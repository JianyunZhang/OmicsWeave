#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
同一实验、同一代谢网络模型下多种示踪剂（如 U-13C 葡萄糖、乳酸钠、甘油、丙氨酸）MFA 结果的整合。

前提
----
在体灌注通常每只动物只输注一种示踪剂：各示踪剂是**互相独立的动物队列**，各自按同一模型拟合。
对同一条反应，各示踪剂测的是**同一个通量**（同一模型、同一单位、同一生理状态），但精度不同：
葡萄糖主要约束糖酵解与 PDH，甘油约束甘油激酶与糖异生入口，丙氨酸约束 GPT / PC，乳酸约束 LDH 与
Cori 循环。示踪剂约束不到的反应，拟合值由约束与初值决定（可辨识性差，或贴在盒约束边界上）。

因此按“示踪剂分层、精度加权”合并，而不是：

- 把各示踪剂的动物混成一组比较：未被约束的估计会稀释有信息的那一个，示踪剂间的系统差异混进组内方差；
- 方向投票（旧流程 ``cross_protocol_consensus``）：不看效应量与精度，且要求多个示踪剂同时显著，
  只有一种示踪剂能解析的反应（如甘油激酶）永远无法通过。

方法（每个比较 × 组织 × 反应）
------------------------------
1. 每个示踪剂内，每只动物一个代表净通量（拟合点 IQR 过滤后的中位数，与单示踪剂分析相同），
   拟合点的方差作为该估计的不确定度；
2. 可辨识性门控：该反应在某示踪剂下为 poorly_constrained / bound_pinned 时，该示踪剂不参与合并；
3. 各示踪剂的组间差 Δ_t 与方差 V_t = (s²_A + f̄_A)/n_A + (s²_B + f̄_B)/n_B，s² 为动物间方差，
   f̄ 为拟合点方差的均值，约束得不好的估计因此自动降权；
4. 逆方差固定效应合并 Δ̂ = Σ w_t Δ_t / Σ w_t（w_t = 1/V_t）。各示踪剂估计的是同一个比较效应，
   固定效应是先验正确的模型；示踪剂间方向相反且异质性显著，说明模型或“同一代谢状态”的假设对某个示踪剂不成立
   （例如示踪剂输注本身扰动了通量），此时不做合并结论，标为 ``multi_tracer_conflict``。异质性在标准化尺度上判断
   （``smd_*``：各示踪剂的 Hedges' g 按 1/Var(g) 合并后的 Cochran Q）：各示踪剂模型的绝对通量可差数十倍，
   绝对尺度的 Q（``heterogeneity_*``，仍保留）主要反映尺度差异。方向一致、只是幅度不一致时不判为冲突，标记
   ``tracer_magnitude_heterogeneity``（每组 3 只时 g 的抽样变异很大）。合并的 ``effect_size`` 是标准化尺度上的合并 g；
5. p 值由示踪剂内的组标签置换校准（分层置换，各示踪剂内组大小不变），不依赖正态近似。示踪剂越多，
   可达的最小 p 越小：单示踪剂 3 vs 3 为 0.1，四个示踪剂各 3 vs 3 为 2/20⁴；
6. BH-FDR 按比较校正；入选规则与单示踪剂分析相同（p、FDR 与 |log2FC| 阈值）。最小可达 p 仍高于阈值时，
   与单示踪剂一样改用效应量（合并的 Hedges' g 与 |log2FC|）。

另外输出各示踪剂对每条反应的约束情况（:func:`tracer_coverage`），以及按条件合并的绝对通量和示踪剂间
一致性（:func:`flux_levels`、:func:`tracer_concordance`），用于发现系统性偏离其余示踪剂的示踪剂。

方法学边界：最理想的做法是在拟合阶段把各示踪剂的 MID 联合拟合到同一组通量（并行标记实验，
Crown & Antoniewicz 2013, PMID 23246523；COMPLETE-MFA, Leighty & Antoniewicz 2013, PMID 24021936）。
在体每只动物只输一种示踪剂时，联合拟合只能在组水平进行，失去以动物为单位的推断；本模块在保留生物学重复的
前提下，合并各示踪剂独立拟合的结果。
"""
from __future__ import annotations

import math
import zlib
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from omicsweave.core.idmap.reaction_pathway_mapping import annotate_flux_dataframe
from omicsweave.core.model import Omics

from .df_analysis import UNCONSTRAINED_LABELS, _log2fc_from_means, _standardized_mean_difference

# tracer_evidence：合并结论的类别
SINGLE = "single_tracer"
CONSISTENT = "multi_tracer_consistent"
CONFLICT = "multi_tracer_conflict"
NOT_IDENTIFIABLE = "not_identifiable"
INSUFFICIENT = "insufficient"
SUPPORTED = (SINGLE, CONSISTENT)
METHOD = "tracer_stratified_permutation"
ALL_TISSUES = "all_tissues"
# 可辨识性从好到差；合并结果的 identifiability 取参与合并的示踪剂中最好的一档
_IDENTIFIABILITY_ORDER = ("well_constrained", "moderate", "undetermined", "poorly_constrained", "bound_pinned")
# 标准化合并权重模式、置换尺度模式与多重检验校正域
SMD_WEIGHTINGS = ("smd", "reliability_smd", "identifiability_smd")
PERMUTATION_SCALES = ("absolute", "relative")
FDR_SCOPES = ("per_tissue", "per_contrast")
_IDENTIFIABILITY_QUALITY = {"well_constrained": 1.0, "moderate": 0.5, "undetermined": 0.25}
# 单个示踪剂内的组划分数不超过该值时全部枚举，否则随机抽取
_MAX_ENUMERATED_SPLITS = 50_000
_TIE_TOLERANCE = 1e-9
_Z_CRITICAL = float(stats.norm.ppf(0.975))


@dataclass(frozen=True)
class PoolingSettings:
    """合并参数：``omics.mfa.analysis.multi_tracer``；判定阈值沿用 ``analysis.differential_flux``。"""

    exclude_unconstrained: bool = True      # poorly_constrained / bound_pinned 的示踪剂不参与合并
    include_fit_variance: bool = True       # 方差中计入每只动物拟合点的方差
    heterogeneity_alpha: float = 0.10       # Cochran Q 的 p 低于该值判为示踪剂间冲突
    min_tracers: int = 1                    # 入选至少需要的参与合并的示踪剂数
    max_exact_permutations: int = 200_000   # 分层置换的组合总数不超过该值时精确枚举
    n_permutations: int = 10_000            # 否则随机置换的次数
    random_state: int = 0
    smd_weighting: str = "smd"              # 标准化权重模式：smd | reliability_smd | identifiability_smd
    permutation_scale: str = "absolute"     # 分层置换统计量尺度：absolute | relative（层内均值幅度归一化）
    n_bootstrap: int = 2000                 # 分层自助法重抽样次数（计算 direction_stability）
    fdr_scope: str = "per_tissue"           # BH-FDR 校正域：per_tissue（按组织分层）| per_contrast（全组织合并）
    # 以下取自 analysis.differential_flux / preprocessing / identifiability
    pvalue_threshold: float = 0.05
    fdr_threshold: float = 0.10
    log2fc_threshold: float = 0.5
    hedges_g_threshold: float = 0.5
    zero_threshold: float = 1.0e-10
    bound_pinned_min_fraction: float = 0.5

    def __post_init__(self) -> None:
        if self.smd_weighting not in SMD_WEIGHTINGS:
            raise ValueError(f"analysis.multi_tracer.smd_weighting 必须为 {' / '.join(SMD_WEIGHTINGS)} 之一，"
                             f"实际为 {self.smd_weighting!r}")
        if self.permutation_scale not in PERMUTATION_SCALES:
            raise ValueError(f"analysis.multi_tracer.permutation_scale 必须为 {' / '.join(PERMUTATION_SCALES)} 之一，"
                             f"实际为 {self.permutation_scale!r}")
        if self.fdr_scope not in FDR_SCOPES:
            raise ValueError(f"analysis.multi_tracer.fdr_scope 必须为 {' / '.join(FDR_SCOPES)} 之一，"
                             f"实际为 {self.fdr_scope!r}")
        if int(self.n_bootstrap) < 0:
            raise ValueError(f"analysis.multi_tracer.n_bootstrap 不能为负数，实际为 {self.n_bootstrap!r}")

    @classmethod
    def from_analysis(cls, analysis: dict[str, Any] | None) -> PoolingSettings:
        analysis = analysis or {}
        own = dict(analysis.get("multi_tracer") or {})
        unknown = sorted(set(own) - set(_OWN_SETTINGS))
        if unknown:
            raise ValueError(f"analysis.multi_tracer 中有无法识别的参数：{'、'.join(unknown)}"
                             f"（可选：{'、'.join(_OWN_SETTINGS)}）")
        diff = analysis.get("differential_flux") or {}
        if "fdr_scope" not in own and "fdr_scope" in diff:
            own["fdr_scope"] = str(diff["fdr_scope"]).lower()
        cohen_d = float(diff.get("cohen_d_threshold", cls.hedges_g_threshold))
        return cls(
            **own,
            pvalue_threshold=float(diff.get("pvalue_threshold", cls.pvalue_threshold)),
            fdr_threshold=float(diff.get("fdr_threshold", cls.fdr_threshold)),
            log2fc_threshold=float(diff.get("log2fc_threshold", cls.log2fc_threshold)),
            hedges_g_threshold=float(diff.get("hedges_g_threshold", cohen_d)),
            zero_threshold=float((analysis.get("preprocessing") or {}).get("zero_threshold", cls.zero_threshold)),
            bound_pinned_min_fraction=float((analysis.get("identifiability") or {}).get(
                "bound_pinned_min_fraction", cls.bound_pinned_min_fraction)),
        )


# analysis.multi_tracer 可设置的参数；其余阈值取自 differential_flux / preprocessing
_OWN_SETTINGS = ("exclude_unconstrained", "include_fit_variance", "heterogeneity_alpha", "min_tracers",
                 "max_exact_permutations", "n_permutations", "random_state", "smd_weighting", "permutation_scale",
                 "n_bootstrap", "fdr_scope")


# ── 各示踪剂对反应的约束 ─────────────────────────────────────────────────────

def tracer_coverage(identifiability: pd.DataFrame, settings: PoolingSettings) -> pd.DataFrame:
    """每个（示踪剂, 组织, 反应）的可辨识性，以及是否参与合并（``informative``）。

    ``identifiability`` 为各示踪剂 ``sample_flux_matrix`` 写出的可辨识性表，另加 ``tracer`` 列。
    ``n_informative_tracers`` 是该反应有多少个示踪剂能约束它，可直接用来看示踪剂之间的互补关系。
    """
    columns = ["tracer", "tissue", "reaction", "identifiability", "median_rel_iqr", "frac_at_bound", "n_samples"]
    table = identifiability[[c for c in columns if c in identifiability.columns]].copy()
    table["identifiability"] = table["identifiability"].fillna("undetermined").astype(str)
    table["informative"] = (~table["identifiability"].isin(UNCONSTRAINED_LABELS)
                            if settings.exclude_unconstrained else True)
    table["n_informative_tracers"] = table.groupby(["tissue", "reaction"])["informative"].transform("sum").astype(int)
    return table.sort_values(["tissue", "reaction", "tracer"]).reset_index(drop=True)


def _with_fit_variance(animals: pd.DataFrame, settings: PoolingSettings) -> pd.DataFrame:
    data = animals.copy()
    if settings.include_fit_variance and "fit_std" in data.columns:
        data["fit_var"] = np.square(pd.to_numeric(data["fit_std"], errors="coerce")).fillna(0.0)
    else:
        data["fit_var"] = 0.0
    return data


# ── 分层置换 ──────────────────────────────────────────────────────────────────

@dataclass
class _Stratum:
    """一个示踪剂在某反应上的数据：每只动物的代表通量、拟合方差、是否属于分子组以及可辨识性分档。"""
    tracer: str
    values: np.ndarray
    fit_var: np.ndarray
    is_num: np.ndarray
    identifiability: str = "moderate"

    @property
    def n(self) -> int:
        return len(self.values)

    @property
    def n_num(self) -> int:
        return int(self.is_num.sum())


@lru_cache(maxsize=64)
def _split_masks(n: int, n_num: int) -> np.ndarray:
    """n 只动物中选 n_num 只为分子组的全部划分，每行一种。"""
    masks = np.zeros((math.comb(n, n_num), n), dtype=bool)
    for row, chosen in enumerate(combinations(range(n), n_num)):
        masks[row, list(chosen)] = True
    masks.setflags(write=False)
    return masks


def _split_terms(stratum: _Stratum, masks: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """每种划分下该示踪剂对合并 z 的贡献 (Δ/V, 1/V)；z = ΣΔ/V / √Σ1/V。

    值先居中，避免平方和相减时的精度损失；V 设一个相对下界，零方差的划分得到极大但有限的 z。
    """
    x = stratum.values - stratum.values.mean()
    f = stratum.fit_var
    n_a = stratum.n_num
    n_b = stratum.n - n_a
    m = masks.astype(float)
    sum_a = m @ x
    sq_a = m @ (x * x)
    mean_a, mean_b = sum_a / n_a, -sum_a / n_b
    var_a = np.maximum((sq_a - n_a * mean_a ** 2) / (n_a - 1), 0.0)
    var_b = np.maximum(((x * x).sum() - sq_a - n_b * mean_b ** 2) / (n_b - 1), 0.0)
    fit_a = (m @ f) / n_a
    fit_b = (f.sum() - m @ f) / n_b
    floor = 1e-12 * (x.var() + f.mean())
    v = np.maximum((var_a + fit_a) / n_a + (var_b + fit_b) / n_b, floor)
    return (mean_a - mean_b) / v, 1.0 / v


def _random_masks(n: int, n_num: int, size: int, rng: np.random.Generator) -> np.ndarray:
    if math.comb(n, n_num) <= _MAX_ENUMERATED_SPLITS:
        masks = _split_masks(n, n_num)
        return masks[rng.integers(0, len(masks), size)]
    chosen = rng.random((size, n)).argsort(axis=1)[:, :n_num]
    masks = np.zeros((size, n), dtype=bool)
    np.put_along_axis(masks, chosen, True, axis=1)
    return masks


def _scale_stratum_relative(stratum: _Stratum, zero_threshold: float) -> _Stratum:
    """层内置换不变尺度归一化：按全部动物 ``mean(|values|)`` 缩放，使不同绝对量纲的示踪剂在置换 z 中可比且保持精确可交换性。"""
    scale = float(np.mean(np.abs(stratum.values))) if len(stratum.values) else 0.0
    if not np.isfinite(scale) or scale <= zero_threshold:
        return stratum
    return _Stratum(stratum.tracer, stratum.values / scale, stratum.fit_var / (scale * scale),
                    stratum.is_num, stratum.identifiability)


def _permutation_test(strata: list[_Stratum], settings: PoolingSettings, seed: int) -> dict[str, Any]:
    """合并 z 的双侧分层置换 p 值：每个示踪剂内独立置换组标签。"""
    active = ([_scale_stratum_relative(s, settings.zero_threshold) for s in strata]
              if settings.permutation_scale == "relative" else strata)
    observed = [_split_terms(s, s.is_num[None, :]) for s in active]
    z_obs = sum(a[0] for a, _ in observed) / math.sqrt(sum(b[0] for _, b in observed))
    total = math.prod(math.comb(s.n, s.n_num) for s in active)
    if total <= settings.max_exact_permutations:
        num, den = np.zeros(1), np.zeros(1)
        for stratum in active:
            a, b = _split_terms(stratum, _split_masks(stratum.n, stratum.n_num))
            num = (num[:, None] + a[None, :]).ravel()
            den = (den[:, None] + b[None, :]).ravel()
        extreme = np.abs(num / np.sqrt(den)) >= abs(z_obs) * (1 - _TIE_TOLERANCE)
        # 各示踪剂两组等大时，整体交换组标签得到 −z，最极端的划分至少有两种
        mirrored = all(2 * s.n_num == s.n for s in active)
        return {"pvalue": float(extreme.mean()), "p_floor": (2.0 if mirrored else 1.0) / total,
                "n_permutations": total, "exact": True, "z": z_obs}

    rng = np.random.default_rng(seed)
    size = settings.n_permutations
    num, den = np.zeros(size), np.zeros(size)
    for stratum in active:
        a, b = _split_terms(stratum, _random_masks(stratum.n, stratum.n_num, size, rng))
        num += a
        den += b
    count = int((np.abs(num / np.sqrt(den)) >= abs(z_obs) * (1 - _TIE_TOLERANCE)).sum())
    return {"pvalue": (1 + count) / (1 + size), "p_floor": 1.0 / (1 + size),
            "n_permutations": size, "exact": False, "z": z_obs}


# ── 合并一个比较 ──────────────────────────────────────────────────────────────

def _best_label(labels) -> str:
    present = set(labels)
    return next((label for label in _IDENTIFIABILITY_ORDER if label in present), "undetermined")


def _describe(values: dict[str, float], digits: int = 4) -> str:
    return ";".join(f"{k}:{v:.{digits}g}" for k, v in values.items())


def _smd_raw_weights(g: np.ndarray, n1: np.ndarray, n2: np.ndarray, *, weighting: str = "smd",
                     reliability_inflation: np.ndarray | None = None,
                     identifiability: list[str] | None = None) -> np.ndarray:
    """标准化合并权重（仅用于有限 g 的子集）：支持 ``smd``、``reliability_smd`` 与 ``identifiability_smd``。"""
    base = 1.0 / ((n1 + n2) / (n1 * n2) + g ** 2 / (2.0 * (n1 + n2)))
    if weighting == "reliability_smd" and reliability_inflation is not None:
        infl = np.asarray(reliability_inflation, dtype=float)
        factor = np.where(np.isfinite(infl) & (infl >= 1.0), infl, 1.0)
        return base / factor
    if weighting == "identifiability_smd" and identifiability is not None:
        quality = np.array([_IDENTIFIABILITY_QUALITY.get(str(lbl), 0.5) for lbl in identifiability], dtype=float)
        return base * quality
    return base


def _pool_standardized(g: np.ndarray, n_num: list[int], n_den: list[int],
                       tracers: list[str] | None = None, *, weighting: str = "smd",
                       reliability_inflation: np.ndarray | None = None,
                       identifiability: list[str] | None = None,
                       g_total: np.ndarray | None = None) -> dict[str, Any]:
    """标准化尺度上的合并（G10：相对尺度）：各示踪剂的 Hedges' g 按所选权重模式加权合并。

    - ``smd``（默认，W2）：经典 Hedges & Olkin 逆方差权重 ``1 / Var(g)``；
    - ``reliability_smd``（W3b）：``1 / (Var(g) · (1 + f̄_fit / s²_p))``，拟合越嘈杂权重越低但不压缩 ``g`` 本身；
    - ``identifiability_smd``（W4）：按可辨识性先验（``well_constrained=1.0, moderate=0.5``）乘 ``1 / Var(g)``。

    同时在提供 ``g_total`` 时计算含拟合方差的全方差合并 Hedges' g（``g_total_var``）作为伴随质控读数。"""
    g = np.asarray(g, dtype=float)
    n1, n2 = np.asarray(n_num, dtype=float), np.asarray(n_den, dtype=float)
    finite = np.isfinite(g)
    out: dict[str, Any] = {"smd": np.nan, "smd_se": np.nan, "smd_heterogeneity_q": np.nan,
                           "smd_heterogeneity_pvalue": np.nan, "smd_i_squared": np.nan,
                           "n_tracers_smd": int(finite.sum()), "smd_weights": "",
                           "smd_weighting": weighting, "g_total_var": np.nan}
    if g_total is not None:
        gt = np.asarray(g_total, dtype=float)
        ok_gt = np.isfinite(gt)
        if ok_gt.any():
            wt = 1.0 / ((n1[ok_gt] + n2[ok_gt]) / (n1[ok_gt] * n2[ok_gt]) + gt[ok_gt] ** 2 / (2.0 * (n1[ok_gt] + n2[ok_gt])))
            out["g_total_var"] = float((wt * gt[ok_gt]).sum() / wt.sum())
    if not finite.any():
        return out
    g_ok, n1_ok, n2_ok = g[finite], n1[finite], n2[finite]
    infl_ok = np.asarray(reliability_inflation, dtype=float)[finite] if reliability_inflation is not None else None
    ident_ok = [lbl for lbl, ok in zip(identifiability, finite) if ok] if identifiability is not None else None
    weight = _smd_raw_weights(g_ok, n1_ok, n2_ok, weighting=weighting,
                              reliability_inflation=infl_ok, identifiability=ident_ok)
    w_total = float(weight.sum())
    pooled = float((weight * g_ok).sum() / w_total)
    out.update(smd=pooled, smd_se=float(math.sqrt(1.0 / w_total)))
    if tracers is not None:
        used_tracers = [t for t, ok in zip(tracers, finite) if ok]
        out["smd_weights"] = _describe(dict(zip(used_tracers, weight / w_total)), 3)
    k = len(g_ok)
    if k > 1:
        q = float((weight * (g_ok - pooled) ** 2).sum())
        out.update(smd_heterogeneity_q=q, smd_heterogeneity_pvalue=float(stats.chi2.sf(q, k - 1)),
                   smd_i_squared=max(0.0, (q - (k - 1)) / q) if q > 0 else 0.0)
    return out


def _stratified_stability(strata: list[_Stratum], scales: np.ndarray, weights: np.ndarray,
                          seed: int, n_bootstrap: int) -> tuple[float, float]:
    """基于激活标准化权重与固定层内尺度因子，计算合并效应的：
    1. ``pooled_loo_stability``：分层留一法中保持与全样本合并效应同号的比例（1.0 = 零反号）；
    2. ``direction_stability``：分层组内自助法（各示踪剂各组独立有放回重抽样）保持同号的比例。"""
    usable = []
    for s, c, w in zip(strata, scales, weights):
        if np.isfinite(c) and c > 0.0 and np.isfinite(w) and w > 0.0:
            a, b = s.values[s.is_num], s.values[~s.is_num]
            usable.append((a, b, float(c), float(w), float(c * (a.mean() - b.mean()))))
    if not usable:
        # 所有层组内零方差（如 [5,5,5] vs [1,1,1]）时，改用等权原始组均值差
        usable = [(s.values[s.is_num], s.values[~s.is_num], 1.0, 1.0,
                   float(s.values[s.is_num].mean() - s.values[~s.is_num].mean())) for s in strata]
    w_sum = sum(item[3] for item in usable)
    if w_sum <= 0.0:
        return np.nan, np.nan
    full = sum(item[3] * item[4] for item in usable) / w_sum
    full_sign = float(np.sign(full)) if np.isfinite(full) else 0.0
    if full_sign == 0.0:
        return 0.0, 0.0

    # 分层留一法
    loo_signs = []
    for idx, (a, b, c, w, _) in enumerate(usable):
        other = sum(item[3] * item[4] for j, item in enumerate(usable) if j != idx)
        for i in range(len(a)):
            val = (other + w * c * float(np.delete(a, i).mean() - b.mean())) / w_sum
            loo_signs.append(float(np.sign(val)))
        for i in range(len(b)):
            val = (other + w * c * float(a.mean() - np.delete(b, i).mean())) / w_sum
            loo_signs.append(float(np.sign(val)))
    loo_stab = float(np.mean(np.asarray(loo_signs) == full_sign)) if loo_signs else np.nan

    # 分层组内自助法（向量化）
    if n_bootstrap <= 0:
        return loo_stab, np.nan
    rng = np.random.default_rng((int(seed) + 1) % 2 ** 32)
    boot_num = np.zeros(int(n_bootstrap), dtype=float)
    for a, b, c, w, _ in usable:
        ba = a[rng.integers(0, len(a), (int(n_bootstrap), len(a)))].mean(axis=1)
        bb = b[rng.integers(0, len(b), (int(n_bootstrap), len(b)))].mean(axis=1)
        boot_num += w * c * (ba - bb)
    boot_stab = float(np.mean(np.sign(boot_num / w_sum) == full_sign))
    return loo_stab, boot_stab


def _pool_reaction(strata: list[_Stratum], settings: PoolingSettings, seed: int) -> dict[str, Any]:
    stats_rows = []
    for s in strata:
        a, b = s.values[s.is_num], s.values[~s.is_num]
        fa, fb = s.fit_var[s.is_num], s.fit_var[~s.is_num]
        n1, n2 = len(a), len(b)
        df = n1 + n2 - 2
        var_a, var_b = a.var(ddof=1), b.var(ddof=1)
        f_a, f_b = float(fa.mean()), float(fb.mean())
        mean_var_a = (var_a + f_a) / n1
        mean_var_b = (var_b + f_b) / n2
        _, g = _standardized_mean_difference(a.mean(), b.mean(), math.sqrt(var_a), math.sqrt(var_b), n1, n2)
        sp2_anim = ((n1 - 1) * var_a + (n2 - 1) * var_b) / df if df > 0 else np.nan
        f_pool = ((n1 - 1) * f_a + (n2 - 1) * f_b) / df if df > 0 else 0.0
        sp2_tot = sp2_anim + f_pool if np.isfinite(sp2_anim) else np.nan
        j_corr = (1.0 - 3.0 / (4.0 * df - 1.0)) if df > 0 and 4.0 * df > 1.0 else np.nan
        scale_c = float(j_corr / math.sqrt(sp2_anim)) if (np.isfinite(sp2_anim) and sp2_anim > 0 and np.isfinite(j_corr)) else np.nan
        g_tot = float(j_corr * (a.mean() - b.mean()) / math.sqrt(sp2_tot)) if (np.isfinite(sp2_tot) and sp2_tot > 0 and np.isfinite(j_corr)) else np.nan
        infl = float(1.0 + f_pool / sp2_anim) if (np.isfinite(sp2_anim) and sp2_anim > 0 and np.isfinite(f_pool)) else np.nan
        stats_rows.append((a.mean(), b.mean(), mean_var_a, mean_var_b, g, scale_c, g_tot, infl))
    mean_a, mean_b, var_mean_a, var_mean_b, g, scale_c, g_tot, infl = (np.array(c, dtype=float) for c in zip(*stats_rows))
    delta = mean_a - mean_b
    variance = var_mean_a + var_mean_b
    # 与置换统计量相同的下界：零方差的示踪剂（如所有动物拟合值相同）仍可参与，但权重有上限
    floor = np.array([1e-12 * (s.values.var() + s.fit_var.mean()) for s in strata])
    weight = 1.0 / np.maximum(variance, floor)
    total = weight.sum()
    diff = float((weight * delta).sum() / total)
    diff_se = math.sqrt(1.0 / total)
    k = len(strata)
    q = float((weight * (delta - diff) ** 2).sum()) if k > 1 else np.nan
    pooled_mean_a = float((weight * mean_a).sum() / total)
    pooled_mean_b = float((weight * mean_b).sum() / total)
    se_a = math.sqrt(float((weight ** 2 * var_mean_a).sum())) / total
    se_b = math.sqrt(float((weight ** 2 * var_mean_b).sum())) / total
    near_zero = min(abs(pooled_mean_a), abs(pooled_mean_b)) <= settings.zero_threshold
    log2fc_se = (np.nan if near_zero else
                 math.sqrt((se_a / pooled_mean_a) ** 2 + (se_b / pooled_mean_b) ** 2) / math.log(2.0))
    tracers = [s.tracer for s in strata]
    idents = [s.identifiability for s in strata]
    n_nums = [s.n_num for s in strata]
    n_dens = [s.n - s.n_num for s in strata]
    perm = _permutation_test(strata, settings, seed)
    standardized = _pool_standardized(
        g, n_nums, n_dens, tracers,
        weighting=settings.smd_weighting,
        reliability_inflation=infl,
        identifiability=idents,
        g_total=g_tot,
    )
    finite_g = np.isfinite(g)
    smd_w_all = np.zeros(len(strata), dtype=float)
    if finite_g.any():
        smd_w_all[finite_g] = _smd_raw_weights(
            g[finite_g], np.asarray(n_nums, float)[finite_g], np.asarray(n_dens, float)[finite_g],
            weighting=settings.smd_weighting,
            reliability_inflation=infl[finite_g],
            identifiability=[lbl for lbl, ok in zip(idents, finite_g) if ok],
        )
    pooled_loo, boot_stab = _stratified_stability(strata, scale_c, smd_w_all, seed, int(settings.n_bootstrap))
    return {
        **standardized,
        "tracer_smd": _describe(dict(zip(tracers, g)), 3),
        "k": k,
        "mean_numerator": pooled_mean_a,
        "mean_denominator": pooled_mean_b,
        "diff": diff,
        "diff_se": diff_se,
        "ci_lower": diff - _Z_CRITICAL * diff_se,
        "ci_upper": diff + _Z_CRITICAL * diff_se,
        "log2fc": _log2fc_from_means(pooled_mean_a, pooled_mean_b, diff, settings.zero_threshold),
        "log2fc_se": log2fc_se,
        # 合并的 Hedges' g 在标准化尺度上按选定权重模式加权（见 _pool_standardized）
        "effect_size": standardized["smd"],
        "statistic": perm["z"],
        "permutation_pvalue": perm["pvalue"],
        "min_achievable_p": perm["p_floor"],
        "n_permutations": perm["n_permutations"],
        "permutation_exact": perm["exact"],
        "permutation_scale": settings.permutation_scale,
        "pooled_loo_stability": pooled_loo,
        "direction_stability": boot_stab,
        "heterogeneity_q": q,
        "heterogeneity_pvalue": float(stats.chi2.sf(q, k - 1)) if k > 1 else np.nan,
        "i_squared": (max(0.0, (q - (k - 1)) / q) if q > 0 else 0.0) if k > 1 else np.nan,
        "n_tracers_up": int((delta > 0).sum()),
        "n_tracers_down": int((delta < 0).sum()),
        "tracer_effects": _describe(dict(zip(tracers, delta))),
        "tracer_weights": _describe(dict(zip(tracers, weight / total)), 3),
    }


def pool_contrast(animals: pd.DataFrame, groups: dict[str, str], numerator: str, denominator: str,
                  coverage: pd.DataFrame, settings: PoolingSettings, seed_key: str = "") -> pd.DataFrame:
    """一个比较：每个（组织, 反应）一行合并结果；多重校正与入选判定见 :func:`to_differential`。

    ``animals`` 每只动物每条反应一行（``sample_id``、``tracer``、``tissue``、``reaction``、``net_flux``、
    ``fit_std``），即各示踪剂 ``sample_flux_matrix`` 的结果加 ``tracer`` 列；``groups`` 为样本 → 所在组
    （分子或分母的水平）。``n_numerator`` / ``n_denominator`` 是均值背后的动物数。
    """
    data = _with_fit_variance(animals[animals["sample_id"].isin(groups)], settings)
    data["is_num"] = data["sample_id"].map(groups).eq(numerator).to_numpy()
    labels = coverage.set_index(["tracer", "tissue", "reaction"])["identifiability"].to_dict()
    informative = set(coverage.loc[coverage["informative"], ["tracer", "tissue", "reaction"]].itertuples(index=False,
                                                                                                         name=None))
    rows = []
    for (tissue, reaction), part in data.groupby(["tissue", "reaction"], sort=True):
        included, excluded, described = [], [], []
        for tracer, sub in part.groupby("tracer", sort=True):
            label = labels.get((str(tracer), tissue, reaction), "undetermined")
            stratum = _Stratum(str(tracer), sub["net_flux"].to_numpy(float), sub["fit_var"].to_numpy(float),
                               sub["is_num"].to_numpy(bool), label)
            n_a, n_b = stratum.n_num, stratum.n - stratum.n_num
            if n_a < 1 or n_b < 1:
                continue
            described.append(stratum)
            contrast_pinned = False
            if settings.exclude_unconstrained and "at_bound" in sub.columns:
                ab = sub["at_bound"].fillna(False).to_numpy(dtype=bool)
                num_ab = float(ab[stratum.is_num].mean()) if n_a > 0 else 0.0
                den_ab = float(ab[~stratum.is_num].mean()) if n_b > 0 else 0.0
                contrast_pinned = max(num_ab, den_ab) >= float(settings.bound_pinned_min_fraction)
            if n_a < 2 or n_b < 2:
                excluded.append(f"{stratum.tracer}:low_n")
            elif (stratum.tracer, tissue, reaction) not in informative:
                excluded.append(f"{stratum.tracer}:{label}")
            elif contrast_pinned:
                excluded.append(f"{stratum.tracer}:contrast_bound_pinned")
            elif stratum.values.var() + stratum.fit_var.mean() <= 0:
                excluded.append(f"{stratum.tracer}:no_variation")
            else:
                included.append(stratum)
        if not described:
            continue
        row: dict[str, Any] = {"tissue": tissue, "reaction": reaction,
                               "n_tracers": len(described), "n_tracers_informative": len(included),
                               "tracers_informative": ";".join(s.tracer for s in included),
                               "tracers_excluded": ";".join(excluded)}
        if included:
            seed = (settings.random_state + zlib.crc32(f"{seed_key}|{tissue}--{reaction}".encode())) % 2 ** 32
            row.update(_pool_reaction(included, settings, seed))
            row["identifiability"] = _best_label(labels.get((s.tracer, tissue, reaction), "undetermined")
                                                 for s in included)
            used = included
        else:
            # 没有可合并的示踪剂：均值只作描述（各示踪剂组均值的简单平均），不做推断
            row.update({"k": 0,
                        "mean_numerator": float(np.mean([s.values[s.is_num].mean() for s in described])),
                        "mean_denominator": float(np.mean([s.values[~s.is_num].mean() for s in described]))})
            row["identifiability"] = _best_label(
                "bound_pinned" if any(e == f"{s.tracer}:contrast_bound_pinned" for e in excluded)
                else labels.get((s.tracer, tissue, reaction), "undetermined")
                for s in described
            )
            used = described
        row["n_numerator"] = sum(s.n_num for s in used)
        row["n_denominator"] = sum(s.n - s.n_num for s in used)
        row["tracer_evidence"] = _classify(row, excluded, settings)
        rows.append(row)
    return pd.DataFrame(rows)


def _classify(row: dict[str, Any], excluded: list[str], settings: PoolingSettings) -> str:
    k = row["k"]
    if k == 0:
        low_n = all(item.endswith(":low_n") for item in excluded)
        return INSUFFICIENT if low_n else NOT_IDENTIFIABLE
    if k == 1:
        return SINGLE
    return CONFLICT if heterogeneous(row, settings) and _opposed(row) else CONSISTENT


def heterogeneous(row: dict[str, Any], settings: PoolingSettings) -> bool:
    """示踪剂间异质性：在标准化尺度上判断（两种以上示踪剂的 g 有定义时），否则退回绝对尺度。"""
    standardized = row.get("smd_heterogeneity_pvalue", np.nan)
    pvalue = standardized if pd.notna(standardized) else row.get("heterogeneity_pvalue", np.nan)
    return bool(pd.notna(pvalue) and pvalue < settings.heterogeneity_alpha)


def _opposed(row: dict[str, Any]) -> bool:
    """参与合并的示踪剂中既有组间差为正的，也有为负的。"""
    return int(row.get("n_tracers_up", 0) or 0) > 0 and int(row.get("n_tracers_down", 0) or 0) > 0


def to_differential(pooled: pd.DataFrame, *, contrast_id: str, contrast_role: str, flux_units: str,
                    settings: PoolingSettings, groups: dict[tuple[str, str], str] | None = None) -> pd.DataFrame:
    """合并结果 → 统一差异结果表（``omicsweave.contract.tables.DIFFERENTIAL_COLUMNS``）。

    - 置换检验在该设计下可达 p 阈值时 ``evidence_level = statistical_test``：p 为置换 p，BH-FDR 在本比较内校正，
      入选需 p、FDR 与 |log2FC| 都过阈值；``groups``（:mod:`omicsweave.mfa.flux_groups`）给出时，化学计量上锁死的
      反应在校正中只计一次；
    - 不可达时与单示踪剂分析一样只给效应量（``effect_size_only``），按合并的 |Hedges' g| 与 |log2FC| 判定；
    - 没有可合并的示踪剂为 ``insufficient``；示踪剂之间冲突（``multi_tracer_conflict``）或参与合并的示踪剂
      少于 ``min_tracers`` 时不入选。
    """
    if pooled.empty:
        return pd.DataFrame()
    table = pooled.copy()
    for column in ("log2fc", "log2fc_se", "effect_size", "statistic", "permutation_pvalue", "min_achievable_p",
                   "diff", "diff_se", "ci_lower", "ci_upper", "heterogeneity_q", "heterogeneity_pvalue", "i_squared"):
        if column not in table.columns:
            table[column] = np.nan
    inferred = table["k"] > 0
    underpowered = inferred & (table["min_achievable_p"] > settings.pvalue_threshold)
    statistical = inferred & ~underpowered
    evidence_level = np.select([statistical, underpowered], ["statistical_test", "effect_size_only"], "insufficient")

    pvalue = table["permutation_pvalue"].where(statistical)
    padj_global = pd.Series(np.nan, index=table.index)
    padj_tissue = pd.Series(np.nan, index=table.index)
    tested = pvalue.notna()
    keys = None
    if groups is not None:
        from .flux_groups import group_keys, grouped_bh
        keys = group_keys(table["tissue"], table["reaction"], groups)
        padj_global = grouped_bh(pvalue, keys)
        for _, part_idx in table.groupby("tissue", sort=False).groups.items():
            padj_tissue.loc[part_idx] = grouped_bh(pvalue.loc[part_idx], keys.loc[part_idx])
    elif tested.any():
        from statsmodels.stats.multitest import multipletests
        padj_global[tested] = multipletests(pvalue[tested], method="fdr_bh")[1]
        for _, part_idx in table.groupby("tissue", sort=False).groups.items():
            sub_tested = tested.loc[part_idx]
            if sub_tested.any():
                idx_valid = part_idx[sub_tested.to_numpy()]
                padj_tissue.loc[idx_valid] = multipletests(pvalue.loc[idx_valid], method="fdr_bh")[1]
    padj = padj_tissue if settings.fdr_scope == "per_tissue" else padj_global

    supported = table["tracer_evidence"].isin(SUPPORTED) & (table["k"] >= settings.min_tracers)
    fc_ok = table["log2fc"].abs().ge(settings.log2fc_threshold)
    p_ok = pvalue.lt(settings.pvalue_threshold) & padj.le(settings.fdr_threshold)
    g_ok = table["effect_size"].abs().ge(settings.hedges_g_threshold)
    selected = supported & ((statistical & p_ok & fc_ok) | (underpowered & g_ok & fc_ok))
    rule_p = (f"permutation p<{settings.pvalue_threshold:g} and fdr_bh<={settings.fdr_threshold:g} "
              f"({settings.fdr_scope}) and |log2FC|>={settings.log2fc_threshold:g}")
    rule_g = f"|Hedges' g|>={settings.hedges_g_threshold:g} and |log2FC|>={settings.log2fc_threshold:g}"
    scope = f"; tracer_evidence in {'/'.join(SUPPORTED)}; informative tracers>={settings.min_tracers}"
    selection_rule = np.select([statistical, underpowered], [rule_p + scope, rule_g + scope], "")

    flags = []
    for i, row in table.iterrows():
        marks = []
        if row["tracer_evidence"] == INSUFFICIENT:
            marks.append("low_n")
        elif row["tracer_evidence"] == NOT_IDENTIFIABLE:
            marks.append(NOT_IDENTIFIABLE)
        elif not np.isfinite(row["effect_size"]):
            marks.append("effect_size_undefined")
        if underpowered[i]:
            marks.append("underpowered")
        if row["tracer_evidence"] == CONFLICT:
            marks.append("tracer_conflict")
        elif row["tracer_evidence"] == SINGLE:
            marks.append(SINGLE)
        elif row["tracer_evidence"] == CONSISTENT and heterogeneous(row, settings):
            marks.append("tracer_magnitude_heterogeneity")     # 方向一致、幅度不一致：不阻止结论，单独标出
        if row["identifiability"] in UNCONSTRAINED_LABELS:
            marks.append(f"identifiability:{row['identifiability']}")
        flags.append(";".join(marks))

    out = pd.DataFrame({
        "contrast_id": contrast_id,
        "contrast_role": contrast_role,
        "tissue": table["tissue"].astype(str),
        "omics": Omics.MFA.value,
        "feature_id": table["reaction"].astype(str),
        "log2fc": table["log2fc"].astype(float),
        "effect_size": table["effect_size"].astype(float),
        "effect_size_type": "hedges_g",
        "statistic": table["statistic"].where(statistical).astype(float),
        "pvalue": pvalue.astype(float),
        "padj": padj.astype(float),
        "padj_global": padj_global.astype(float),
        "fdr_scope": settings.fdr_scope,
        "padj_method": np.where(statistical, "fdr_bh", "none"),
        "n_numerator": table["n_numerator"].astype("int64"),
        "n_denominator": table["n_denominator"].astype("int64"),
        "min_achievable_p": table["min_achievable_p"].astype(float),
        "method": np.where(underpowered, "effect_size_only", METHOD),
        "mean_numerator": table["mean_numerator"].astype(float),
        "mean_denominator": table["mean_denominator"].astype(float),
        "units": flux_units,
        "inference_unit": "biological_replicate",
        "evidence_level": evidence_level,
        "selected": selected.astype(bool),
        "selection_rule": selection_rule,
        "identifiability": table["identifiability"].astype(str),
        "direction_stability": pd.to_numeric(table["direction_stability"], errors="coerce").astype(float)
        if "direction_stability" in table.columns else np.nan,
        "permutation_pvalue": table["permutation_pvalue"].astype(float),
        "flags": flags,
        "flux_id": table["tissue"].astype(str) + "--" + table["reaction"].astype(str),
    })
    for column in ("tracer_evidence", "n_tracers", "n_tracers_informative", "tracers_informative", "tracers_excluded",
                   "diff", "diff_se", "ci_lower", "ci_upper", "log2fc_se", "heterogeneity_q", "heterogeneity_pvalue",
                   "i_squared", "n_tracers_up", "n_tracers_down", "tracer_effects", "tracer_weights",
                   "n_permutations", "permutation_exact", "permutation_scale", "pooled_loo_stability",
                   "smd", "smd_se", "smd_heterogeneity_q", "smd_heterogeneity_pvalue", "smd_i_squared",
                   "n_tracers_smd", "tracer_smd", "smd_weights", "smd_weighting", "g_total_var"):
        if column in table.columns:
            out[column] = table[column].to_numpy()
    if keys is not None:
        out["flux_group"] = [groups.get((str(t), str(r)), "") for t, r in zip(table["tissue"], table["reaction"])]
        out["flux_group_size"] = keys.map(keys.value_counts()).astype(int).to_numpy()
    annotations = annotate_flux_dataframe(pd.DataFrame({"reaction": out["feature_id"]}))
    for column in ("kegg_reaction", "kegg_pathway", "reactome_id", "ec_number"):
        out[column] = annotations[column]
    return out.reset_index(drop=True)


# ── 按条件合并的绝对通量与示踪剂一致性 ────────────────────────────────────────

def flux_levels(animals: pd.DataFrame, conditions: dict[str, str], coverage: pd.DataFrame,
                settings: PoolingSettings) -> tuple[pd.DataFrame, pd.DataFrame]:
    """每个（组织, 反应, 条件）合并各示踪剂的绝对通量，并检查示踪剂之间是否一致。

    同一条件下各示踪剂测的是同一个通量，均值应在误差范围内一致；不一致（Q 检验显著）说明该示踪剂下模型
    失拟，或示踪剂输注本身改变了代谢状态。返回 ``(levels, deviations)``：``deviations`` 每个
    （示踪剂, 组织, 反应, 条件）一行，为该示踪剂相对其余示踪剂合并值的偏离（留一法）。
    只用参与合并（``informative``）且该条件下至少 2 只动物的示踪剂。
    """
    data = _with_fit_variance(animals[animals["sample_id"].isin(conditions)], settings)
    data["condition"] = data["sample_id"].map(conditions)
    keys = ["tissue", "reaction", "condition"]
    per = (data.groupby(["tracer", *keys])
           .agg(n=("net_flux", "size"), mean=("net_flux", "mean"), var=("net_flux", "var"), fit=("fit_var", "mean"))
           .reset_index())
    informative = coverage.loc[coverage["informative"], ["tracer", "tissue", "reaction"]]
    per = per.merge(informative, on=["tracer", "tissue", "reaction"])
    per = per[per["n"] >= 2].copy()
    per["v"] = (per["var"] + per["fit"]) / per["n"]
    per = per[per["v"] > 0].copy()
    if per.empty:
        return pd.DataFrame(), pd.DataFrame()
    per["w"] = 1.0 / per["v"]
    per["wm"] = per["w"] * per["mean"]
    per["sum_w"] = per.groupby(keys)["w"].transform("sum")
    per["sum_wm"] = per.groupby(keys)["wm"].transform("sum")
    per["k"] = per.groupby(keys)["w"].transform("size")
    per["pooled"] = per["sum_wm"] / per["sum_w"]
    per["q_term"] = per["w"] * (per["mean"] - per["pooled"]) ** 2
    per["described"] = per["tracer"] + ":" + per["mean"].map(lambda v: f"{v:.4g}")

    levels = (per.groupby(keys)
              .agg(n_tracers=("tracer", "size"), tracers=("tracer", ";".join), n_animals=("n", "sum"),
                   pooled_flux=("pooled", "first"), sum_w=("sum_w", "first"), heterogeneity_q=("q_term", "sum"),
                   min_tracer_flux=("mean", "min"), max_tracer_flux=("mean", "max"),
                   tracer_means=("described", ";".join))
              .reset_index())
    levels["pooled_se"] = np.sqrt(1.0 / levels.pop("sum_w"))
    multi = levels["n_tracers"] > 1
    dof = levels["n_tracers"] - 1
    q = levels["heterogeneity_q"].where(multi)
    levels["heterogeneity_q"] = q
    levels["heterogeneity_pvalue"] = np.where(multi, stats.chi2.sf(q.fillna(0.0), dof.clip(lower=1)), np.nan)
    levels["i_squared"] = np.where(multi, np.where(q > 0, ((q - dof) / q).clip(lower=0), 0.0), np.nan)
    consistent = levels["heterogeneity_pvalue"] >= settings.heterogeneity_alpha
    levels["tracers_consistent"] = consistent.astype(object).where(multi, None)

    others = per[per["k"] > 1].copy()
    rest_w = others["sum_w"] - others["w"]
    others["others_flux"] = (others["sum_wm"] - others["wm"]) / rest_w
    others["z_deviation"] = (others["mean"] - others["others_flux"]) / np.sqrt(others["v"] + 1.0 / rest_w)
    scale = others["others_flux"].abs()
    others["rel_deviation"] = ((others["mean"] - others["others_flux"]) / scale).where(scale > settings.zero_threshold)
    deviations = others.rename(columns={"mean": "tracer_flux", "n": "n_animals"})[
        ["tracer", *keys, "n_animals", "tracer_flux", "others_flux", "z_deviation", "rel_deviation"]]
    return levels.reset_index(drop=True), deviations.reset_index(drop=True)


def tracer_concordance(deviations: pd.DataFrame) -> pd.DataFrame:
    """按（示踪剂, 组织）汇总留一法偏离：系统性偏离其余示踪剂的示踪剂在这里最明显。

    ``tissue = all_tissues`` 的行汇总全部组织。``frac_beyond_1_96`` 在各示踪剂一致时约为 5%；
    ``median_rel_deviation`` 的符号说明该示踪剂给出的通量系统性偏高还是偏低。只有一个示踪剂偏离时，
    其余示踪剂的留一法参照值也含有它，偏离会部分摊到它们身上，但偏离者的 ``max_abs_z`` 最大。
    """
    if deviations.empty:
        return pd.DataFrame()
    data = pd.concat([deviations, deviations.assign(tissue=ALL_TISSUES)], ignore_index=True)
    data["abs_z"] = data["z_deviation"].abs()
    summary = (data.groupby(["tracer", "tissue"], sort=True)
               .agg(n_compared=("abs_z", "size"), median_abs_z=("abs_z", "median"), max_abs_z=("abs_z", "max"),
                    frac_beyond_1_96=("abs_z", lambda s: float((s > _Z_CRITICAL).mean())),
                    median_rel_deviation=("rel_deviation", "median"))
               .reset_index())
    return summary
