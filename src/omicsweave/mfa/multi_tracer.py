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
   固定效应是先验正确的模型；示踪剂间的异质性（Cochran Q）说明模型或“同一代谢状态”的假设对某个示踪剂
   不成立（例如示踪剂输注本身扰动了通量），此时不做合并结论，标为 ``multi_tracer_conflict``；
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
    # 以下取自 analysis.differential_flux / preprocessing
    pvalue_threshold: float = 0.05
    fdr_threshold: float = 0.10
    log2fc_threshold: float = 0.5
    hedges_g_threshold: float = 0.5
    zero_threshold: float = 1.0e-10

    @classmethod
    def from_analysis(cls, analysis: dict[str, Any] | None) -> PoolingSettings:
        analysis = analysis or {}
        own = dict(analysis.get("multi_tracer") or {})
        unknown = sorted(set(own) - set(_OWN_SETTINGS))
        if unknown:
            raise ValueError(f"analysis.multi_tracer 中有无法识别的参数：{'、'.join(unknown)}"
                             f"（可选：{'、'.join(_OWN_SETTINGS)}）")
        diff = analysis.get("differential_flux") or {}
        cohen_d = float(diff.get("cohen_d_threshold", cls.hedges_g_threshold))
        return cls(
            **own,
            pvalue_threshold=float(diff.get("pvalue_threshold", cls.pvalue_threshold)),
            fdr_threshold=float(diff.get("fdr_threshold", cls.fdr_threshold)),
            log2fc_threshold=float(diff.get("log2fc_threshold", cls.log2fc_threshold)),
            hedges_g_threshold=float(diff.get("hedges_g_threshold", cohen_d)),
            zero_threshold=float((analysis.get("preprocessing") or {}).get("zero_threshold", cls.zero_threshold)),
        )


# analysis.multi_tracer 可设置的参数；其余阈值取自 differential_flux / preprocessing
_OWN_SETTINGS = ("exclude_unconstrained", "include_fit_variance", "heterogeneity_alpha", "min_tracers",
                 "max_exact_permutations", "n_permutations", "random_state")


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
    """一个示踪剂在某反应上的数据：每只动物的代表通量、拟合方差，以及是否属于分子组。"""
    tracer: str
    values: np.ndarray
    fit_var: np.ndarray
    is_num: np.ndarray

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


def _permutation_test(strata: list[_Stratum], settings: PoolingSettings, seed: int) -> dict[str, Any]:
    """合并 z 的双侧分层置换 p 值：每个示踪剂内独立置换组标签。"""
    observed = [_split_terms(s, s.is_num[None, :]) for s in strata]
    z_obs = sum(a[0] for a, _ in observed) / math.sqrt(sum(b[0] for _, b in observed))
    total = math.prod(math.comb(s.n, s.n_num) for s in strata)
    if total <= settings.max_exact_permutations:
        num, den = np.zeros(1), np.zeros(1)
        for stratum in strata:
            a, b = _split_terms(stratum, _split_masks(stratum.n, stratum.n_num))
            num = (num[:, None] + a[None, :]).ravel()
            den = (den[:, None] + b[None, :]).ravel()
        extreme = np.abs(num / np.sqrt(den)) >= abs(z_obs) * (1 - _TIE_TOLERANCE)
        # 各示踪剂两组等大时，整体交换组标签得到 −z，最极端的划分至少有两种
        mirrored = all(2 * s.n_num == s.n for s in strata)
        return {"pvalue": float(extreme.mean()), "p_floor": (2.0 if mirrored else 1.0) / total,
                "n_permutations": total, "exact": True, "z": z_obs}

    rng = np.random.default_rng(seed)
    size = settings.n_permutations
    num, den = np.zeros(size), np.zeros(size)
    for stratum in strata:
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


def _pool_reaction(strata: list[_Stratum], settings: PoolingSettings, seed: int) -> dict[str, Any]:
    stats_rows = []
    for s in strata:
        a, b = s.values[s.is_num], s.values[~s.is_num]
        fa, fb = s.fit_var[s.is_num], s.fit_var[~s.is_num]
        var_a, var_b = a.var(ddof=1), b.var(ddof=1)
        mean_var_a = (var_a + fa.mean()) / len(a)
        mean_var_b = (var_b + fb.mean()) / len(b)
        _, g = _standardized_mean_difference(a.mean(), b.mean(), math.sqrt(var_a), math.sqrt(var_b), len(a), len(b))
        stats_rows.append((a.mean(), b.mean(), mean_var_a, mean_var_b, g))
    mean_a, mean_b, var_mean_a, var_mean_b, g = (np.array(c, dtype=float) for c in zip(*stats_rows))
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
    finite_g = np.isfinite(g)
    pooled_mean_a = float((weight * mean_a).sum() / total)
    pooled_mean_b = float((weight * mean_b).sum() / total)
    se_a = math.sqrt(float((weight ** 2 * var_mean_a).sum())) / total
    se_b = math.sqrt(float((weight ** 2 * var_mean_b).sum())) / total
    near_zero = min(abs(pooled_mean_a), abs(pooled_mean_b)) <= settings.zero_threshold
    log2fc_se = (np.nan if near_zero else
                 math.sqrt((se_a / pooled_mean_a) ** 2 + (se_b / pooled_mean_b) ** 2) / math.log(2.0))
    tracers = [s.tracer for s in strata]
    perm = _permutation_test(strata, settings, seed)
    return {
        "k": k,
        "mean_numerator": pooled_mean_a,
        "mean_denominator": pooled_mean_b,
        "diff": diff,
        "diff_se": diff_se,
        "ci_lower": diff - _Z_CRITICAL * diff_se,
        "ci_upper": diff + _Z_CRITICAL * diff_se,
        "log2fc": _log2fc_from_means(pooled_mean_a, pooled_mean_b, diff, settings.zero_threshold),
        "log2fc_se": log2fc_se,
        "effect_size": (float((weight[finite_g] * g[finite_g]).sum() / weight[finite_g].sum())
                        if finite_g.any() else np.nan),
        "statistic": perm["z"],
        "permutation_pvalue": perm["pvalue"],
        "min_achievable_p": perm["p_floor"],
        "n_permutations": perm["n_permutations"],
        "permutation_exact": perm["exact"],
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
            stratum = _Stratum(str(tracer), sub["net_flux"].to_numpy(float), sub["fit_var"].to_numpy(float),
                               sub["is_num"].to_numpy(bool))
            n_a, n_b = stratum.n_num, stratum.n - stratum.n_num
            if n_a < 1 or n_b < 1:
                continue
            described.append(stratum)
            label = labels.get((stratum.tracer, tissue, reaction), "undetermined")
            if n_a < 2 or n_b < 2:
                excluded.append(f"{stratum.tracer}:low_n")
            elif (stratum.tracer, tissue, reaction) not in informative:
                excluded.append(f"{stratum.tracer}:{label}")
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
            row["identifiability"] = _best_label(labels.get((s.tracer, tissue, reaction), "undetermined")
                                                 for s in described)
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
    return CONFLICT if row["heterogeneity_pvalue"] < settings.heterogeneity_alpha else CONSISTENT


def to_differential(pooled: pd.DataFrame, *, contrast_id: str, contrast_role: str, flux_units: str,
                    settings: PoolingSettings) -> pd.DataFrame:
    """合并结果 → 统一差异结果表（``omicsweave.contract.tables.DIFFERENTIAL_COLUMNS``）。

    - 置换检验在该设计下可达 p 阈值时 ``evidence_level = statistical_test``：p 为置换 p，BH-FDR 在本比较内校正，
      入选需 p、FDR 与 |log2FC| 都过阈值；
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
    padj = pd.Series(np.nan, index=table.index)
    tested = pvalue.notna()
    if tested.any():
        from statsmodels.stats.multitest import multipletests
        padj[tested] = multipletests(pvalue[tested], method="fdr_bh")[1]

    supported = table["tracer_evidence"].isin(SUPPORTED) & (table["k"] >= settings.min_tracers)
    fc_ok = table["log2fc"].abs().ge(settings.log2fc_threshold)
    p_ok = pvalue.lt(settings.pvalue_threshold) & padj.le(settings.fdr_threshold)
    g_ok = table["effect_size"].abs().ge(settings.hedges_g_threshold)
    selected = supported & ((statistical & p_ok & fc_ok) | (underpowered & g_ok & fc_ok))
    rule_p = (f"permutation p<{settings.pvalue_threshold:g} and fdr_bh<={settings.fdr_threshold:g} "
              f"and |log2FC|>={settings.log2fc_threshold:g}")
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
        "permutation_pvalue": table["permutation_pvalue"].astype(float),
        "flags": flags,
        "flux_id": table["tissue"].astype(str) + "--" + table["reaction"].astype(str),
    })
    for column in ("tracer_evidence", "n_tracers", "n_tracers_informative", "tracers_informative", "tracers_excluded",
                   "diff", "diff_se", "ci_lower", "ci_upper", "log2fc_se", "heterogeneity_q", "heterogeneity_pvalue",
                   "i_squared", "n_tracers_up", "n_tracers_down", "tracer_effects", "tracer_weights",
                   "n_permutations", "permutation_exact"):
        if column in table.columns:
            out[column] = table[column]
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
