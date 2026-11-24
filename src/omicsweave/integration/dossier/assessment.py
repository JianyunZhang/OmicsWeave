#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
分轴综合评估（G15，设计文档第 3 节）：不做加和，对每个（反应, 组织, 比较）在四条互不依赖的轴上各给出方向与状态，
再由判定表给出关于通量能说到哪一级（L3–L0、LX）与机制类型，联合结论用合取检验，排序用公开的字典序。

四条轴只读 ``dossier_members`` 中各自的成员行，任何一条轴都不读取其他轴的结果：

- F 通量：各示踪剂的拟合通量与合并估计（本实验）。established 要求至少两种示踪剂可辨识且同向、合并估计通过检验、
  较小组至少 3 只、不是由单只动物驱动；suggestive 要求有可辨识的示踪剂、|Hedges' g| ≥ 1 且逐只去除后方向保持；
  示踪剂之间方向不一致为 conflicting；面板声明的预期示踪剂与数据判定的可辨识性不符时最高为 undetermined；
- C 酶容量：直接酶（转运反应为载体），蛋白优先、RNA 单独作为第二读数。同工酶按丰度占比判断主要成员，主要成员
  方向相反为 mixed；复合体按整体判断（已定量的核心亚基至少一半同向且至少一个通过 FDR），共用亚基不参与；没有蛋白时
  RNA 须在两个来源中同向才到 changed（标 rna_only）。规则 ``rna_single_source_changed`` 打开时，只有一个来源测到
  RNA 的反应在该来源中通过 FDR 即为 changed（另标 rna_single_source，不算复现）；两个来源都测到时仍须同向；
- R 调控：带符号的调控因子与别构效应物，成员变化乘以作用方向得到对活性的指向；由调控因子的量最多到 suggestive，
  测得的效应物变化才到 changed；方向未知或经由转录 / 信号起作用的成员不参与；
- M 底物与产物：产物升指向通量升、底物升指向通量降（一步产物 / 底物比与产物同向）；只表示容量的反应不给方向。

状态共用一套词汇：changed（F 轴为 established）、suggestive、unchanged（等效）、mixed（F 轴为 conflicting）、
undetermined、not_identifiable、insufficient_n、not_detected、not_measured、qc_blocked、no_knowledge。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import t as student_t

from ..mappings.reaction_members import ACTIVITY_MODES

RULES: dict[str, Any] = {"alpha": 0.05, "equivalence_bound": 0.58, "flux_min_abs_g": 1.0,
                         # S4：多示踪剂同向、无异质性、分层置换 p < alpha 且方向稳定时，辅助 Suggestive 路径的 |g| 下限
                         # （设为 None 或 >= flux_min_abs_g 时关闭该路径，仅保留 |g| >= flux_min_abs_g 单路径）
                         "flux_concordant_min_abs_g": 0.8,
                         # S5：通量轴 TOST 等效检验的标准化效应 90% 上界门槛（与 log2FC 等效界联合把关，设为 None 或 <= 0 关闭）
                         "flux_equiv_g_bound": 1.0,
                         "heterogeneity_alpha": 0.10,
                         "direction_stability_min": 0.8, "established_min_group": 3, "minor_share": 0.2,
                         "bootstrap": 2000, "seed": 20261117,
                         # L1 要求 R、M 轴没有相反的指向；suggestive（名义 p）即算指向（设计文档 3.5），可改为 changed
                         "opposition_min_state": "suggestive",
                         # 没有蛋白、只有一个来源测到 RNA 时，该来源通过 FDR 是否判为 changed（默认要求两个来源同向）。
                         # 供没有蛋白组、也没有第二个队列的项目使用；打开后应运行诱饵校准，看 C 轴 changed 是否高于诱饵
                         "rna_single_source_changed": False,
                         # MFA 比较存在整鼠通量系统性尺度漂移（flux_scale_drift）时，F 轴禁止升为 suggestive / changed
                         "flux_scale_drift_gate": True}
MEASURED = ("changed", "nominal", "no_evidence", "equivalent")
MISSING_ORDER = ("qc_blocked", "not_detected", "not_measured", "no_knowledge")
_CAUTION_FLAGS = ("single_subject_driven", "foreign_tissue_marker", "qc_caution")
_LAYER_QC_CAUTION = ("foreign_tissue_marker", "qc_caution")
AXIS_COLUMNS = ["flux_id", "tissue", "contrast_id", "axis", "direction", "state", "effect", "ci_lower", "ci_upper",
                "sign_probability", "pvalue_one_sided", "n_expected", "n_measured", "n_sources", "replicated", "basis",
                "rule", "flags"]
REACTION_COLUMNS = ["flux_id", "tissue", "contrast_id", "reaction_class", "flux_claim", "claim_direction",
                    "n_concordant_axes", "discordant_axes", "mechanism_class", "rho_h_lower", "rho_h_upper",
                    "conjunction_statement", "conjunction_pvalue", "conjunction_padj", "any_pvalue", "completeness",
                    "completeness_by_layer", "testable_by", "phenomena", "rank"]
CLAIM_ORDER = {"L3": 4, "L2": 3, "L1": 2, "L0": 1, "LX": 0}
_RELEVANT = ("flux", "direct_enzyme", "carrier", "regulator", "effector", "substrate", "product", "step_ratio")


@dataclass
class Axis:
    axis: str
    state: str
    direction: int = 0
    effect: float = np.nan
    ci_lower: float = np.nan
    ci_upper: float = np.nan
    sign_probability: float = np.nan
    pvalue_one_sided: float = np.nan
    n_expected: int = 0
    n_measured: int = 0
    n_sources: int = 0
    replicated: bool = False
    basis: str = ""
    rule: str = ""
    flags: list[str] = field(default_factory=list)

    @property
    def pointing(self) -> bool:
        """这条轴给出了方向（状态足以说“升”或“降”）。"""
        return self.direction != 0 and self.state in ("changed", "suggestive")


# ---------------------------------------------------------------------------
# 工具
# ---------------------------------------------------------------------------

def _sign(value: Any) -> int:
    value = pd.to_numeric(value, errors="coerce")
    return 0 if pd.isna(value) or value == 0 else int(np.sign(value))


def one_sided(pvalue: Any, fc: Any, direction: int) -> float:
    """双侧 p → 指定方向的单侧 p。"""
    p = pd.to_numeric(pvalue, errors="coerce")
    if pd.isna(p) or direction == 0:
        return np.nan
    return float(p / 2) if _sign(fc) == direction else float(1 - p / 2)


def simes(pvalues: list[float]) -> float:
    values = sorted(p for p in pvalues if pd.notna(p))
    if not values:
        return np.nan
    n = len(values)
    return float(min(1.0, min(p * n / (i + 1) for i, p in enumerate(values))))


def bh(pvalues: pd.Series) -> pd.Series:
    values = pd.to_numeric(pvalues, errors="coerce")
    valid = values.dropna()
    if valid.empty:
        return pd.Series(np.nan, index=values.index)
    order = valid.sort_values()
    n = len(order)
    adjusted = (order * n / np.arange(1, n + 1)).iloc[::-1].cummin().iloc[::-1].clip(upper=1.0)
    return adjusted.reindex(values.index)


def _activity_sign(row: pd.Series) -> int:
    if row.get("regulator_mode") not in ACTIVITY_MODES:
        return 0
    return {"activates": 1, "inhibits": -1}.get(row.get("regulator_effect"), 0)


def _certainty(row: pd.Series) -> float:
    from .view import certainty
    return certainty(row)


def _missing_state(rows: pd.DataFrame) -> str:
    states = set(rows["cell_state"])
    return next((s for s in MISSING_ORDER if s in states), "not_measured")


def _flags(rows: pd.DataFrame) -> set[str]:
    return {f for value in rows["flags"].fillna("") for f in str(value).split(";") if f}


def _caution_info(rows: pd.DataFrame, target_flags: tuple[str, ...] = _CAUTION_FLAGS) -> tuple[pd.Series, list[str]]:
    """检查给定行集合是否全由带警示标记的行驱动，并返回出现的警示标记列表。"""
    if rows.empty:
        return pd.Series(dtype=bool), []
    raw_flags = rows["flags"].fillna("").astype(str) if "flags" in rows.columns else pd.Series("", index=rows.index)
    qc_conf = (
        rows["qc_status"].astype(str).eq("qc_confounded")
        if "qc_status" in rows.columns else pd.Series(False, index=rows.index)
    )
    masks: dict[str, pd.Series] = {}
    for flag in target_flags:
        hit = raw_flags.str.split(";").map(lambda fs, target=flag: target in fs)
        if flag == "qc_caution":
            hit = hit | qc_conf
        masks[flag] = hit
    any_caution = pd.Series(False, index=rows.index)
    for hit in masks.values():
        any_caution = any_caution | hit
    present = [flag for flag in target_flags if bool(masks[flag].any())]
    return any_caution, present


def _cap_if_cautioned(axis: Axis, driving: pd.DataFrame, target_flags: tuple[str, ...] = _CAUTION_FLAGS) -> None:
    """若支撑 changed 判定的全部成员均带有风险标记，则软封顶降级为 suggestive 并保留警示标记。"""
    if axis.state != "changed" or driving.empty:
        return
    caution_mask, caution_flags = _caution_info(driving, target_flags)
    if bool(caution_mask.all()) and caution_flags:
        axis.state = "suggestive"
        for flag in caution_flags:
            if flag not in axis.flags:
                axis.flags.append(flag)


def _describe(rows: pd.DataFrame, label: str = "member_id") -> str:
    parts = []
    for name, fallback, source, state, fc in zip(rows[label], rows["feature_name"], rows["source"],
                                                 rows["cell_state"], rows["log2fc"]):
        value = f", {fc:+.2f}" if pd.notna(fc) else ""
        parts.append(f"{name or fallback}:{source}({state}{value})")
    return "; ".join(parts)


# ---------------------------------------------------------------------------
# F 轴：通量
# ---------------------------------------------------------------------------

def loo_stability(numerator: np.ndarray, denominator: np.ndarray) -> float:
    """逐只去除一只动物后，组均值差与全部动物时同号的比例；任一组少于 2 只时为缺失。"""
    numerator, denominator = np.asarray(numerator, float), np.asarray(denominator, float)
    numerator, denominator = numerator[np.isfinite(numerator)], denominator[np.isfinite(denominator)]
    if len(numerator) < 2 or len(denominator) < 2:
        return np.nan
    full = np.sign(numerator.mean() - denominator.mean())
    if full == 0:
        return 0.0
    signs = [np.sign(np.delete(numerator, i).mean() - denominator.mean()) for i in range(len(numerator))]
    signs += [np.sign(numerator.mean() - np.delete(denominator, i).mean()) for i in range(len(denominator))]
    return float(np.mean(np.array(signs) == full))


def stratified_pooled_loo(strata: list[tuple[np.ndarray, np.ndarray]]) -> float:
    """多示踪剂分层合并的逐只去除（LOO）方向稳定性（S2）：各可辨识示踪剂按完整层内 Hedges' g 尺度与
    逆方差权重固定合并，逐一去除所有层中的每只动物后检验合并估计是否维持原方向。固定尺度避免 n=2 vs 2
    去一后 df=1 导致 J(1)=0 的退化，也不再被弱贡献示踪剂的单层符号波动否决。"""
    cleaned: list[tuple[np.ndarray, np.ndarray]] = []
    for num, den in strata:
        a = np.asarray(num, float)
        b = np.asarray(den, float)
        a, b = a[np.isfinite(a)], b[np.isfinite(b)]
        if len(a) >= 2 and len(b) >= 2:
            cleaned.append((a, b))
    if not cleaned:
        return np.nan
    if len(cleaned) == 1:
        return loo_stability(*cleaned[0])
    prepared: list[tuple[np.ndarray, np.ndarray, float, float, float]] = []
    for a, b in cleaned:
        na, nb = len(a), len(b)
        df = na + nb - 2
        ss = float(np.sum((a - a.mean()) ** 2) + np.sum((b - b.mean()) ** 2))
        sp = float(np.sqrt(ss / df)) if df > 0 else 0.0
        if not np.isfinite(sp) or sp <= 0.0:
            continue
        j = 1.0 - 3.0 / (4.0 * df - 1.0)
        scale = j / sp
        g = float((a.mean() - b.mean()) * scale)
        var_g = (na + nb) / (na * nb) + (g * g) / (2.0 * (na + nb))
        if not np.isfinite(var_g) or var_g <= 0.0:
            continue
        prepared.append((a, b, scale, g, 1.0 / var_g))
    if not prepared:
        signs = [np.sign(a.mean() - b.mean()) for a, b in cleaned]
        return 1.0 if signs and signs[0] != 0 and len(set(signs)) == 1 else np.nan
    if len(prepared) == 1:
        return loo_stability(prepared[0][0], prepared[0][1])
    total_w = sum(w for *_, w in prepared)
    total_wg = sum(w * g for *_, g, w in prepared)
    full_sign = float(np.sign(total_wg / total_w))
    if full_sign == 0.0:
        return 0.0
    same: list[bool] = []
    for a, b, scale, g, w in prepared:
        base_wg = total_wg - w * g
        mean_a, mean_b = float(a.mean()), float(b.mean())
        for i in range(len(a)):
            g_new = (float(np.delete(a, i).mean()) - mean_b) * scale
            same.append(bool(np.sign((base_wg + w * g_new) / total_w) == full_sign))
        for j in range(len(b)):
            g_new = (mean_a - float(np.delete(b, j).mean())) * scale
            same.append(bool(np.sign((base_wg + w * g_new) / total_w) == full_sign))
    return float(np.mean(same))


def _flux_smd_se(head: pd.Series, identifiable: pd.DataFrame) -> tuple[float, float]:
    """合并标准化效应量 Hedges' g 的标准误与总自由度（用于 S5 通量等效检验）。
    优先读取差异表落盘的 ``smd_se``；历史结果包未存该列时，由可辨识示踪剂的 (n_numerator, n_denominator, effect)
    按逆方差 Hedges' g 合并公式回算。"""
    n_pairs: list[tuple[float, float, float]] = []
    for _, tr in identifiable.iterrows():
        na = pd.to_numeric(tr.get("n_numerator"), errors="coerce")
        nb = pd.to_numeric(tr.get("n_denominator"), errors="coerce")
        gt = pd.to_numeric(tr.get("effect"), errors="coerce")
        if pd.notna(na) and pd.notna(nb) and na >= 2 and nb >= 2:
            n_pairs.append((float(na), float(nb), float(gt) if pd.notna(gt) else 0.0))
    df = sum(na + nb - 2.0 for na, nb, _ in n_pairs)
    if df <= 0:
        ha = pd.to_numeric(head.get("n_numerator"), errors="coerce")
        hb = pd.to_numeric(head.get("n_denominator"), errors="coerce")
        if pd.notna(ha) and pd.notna(hb) and ha >= 2 and hb >= 2:
            df = float(ha + hb - 2.0)
    reported_se = pd.to_numeric(head.get("smd_se"), errors="coerce")
    if pd.notna(reported_se) and reported_se > 0 and df >= 1:
        return float(reported_se), float(df)
    if n_pairs and df >= 1:
        inv_vars = []
        for na, nb, gt in n_pairs:
            vg = (na + nb) / (na * nb) + (gt * gt) / (2.0 * (na + nb))
            if np.isfinite(vg) and vg > 0:
                inv_vars.append(1.0 / vg)
        if inv_vars:
            return float(np.sqrt(1.0 / sum(inv_vars))), float(df)
    ha = pd.to_numeric(head.get("n_numerator"), errors="coerce")
    hb = pd.to_numeric(head.get("n_denominator"), errors="coerce")
    gh = pd.to_numeric(head.get("effect"), errors="coerce")
    if pd.notna(ha) and pd.notna(hb) and ha >= 2 and hb >= 2 and pd.notna(gh):
        vg = (float(ha) + float(hb)) / (float(ha) * float(hb)) + (float(gh) ** 2) / (2.0 * (float(ha) + float(hb)))
        if np.isfinite(vg) and vg > 0:
            return float(np.sqrt(vg)), float(ha + hb - 2.0)
    return np.nan, 0.0


def flux_axis(rows: pd.DataFrame, rules: dict[str, Any], expected: list[str] | None = None,
              stability: dict[str, float] | None = None) -> Axis:
    """``rows``：本实验在该（反应, 组织, 比较）的通量行（各示踪剂与 pooled）；``stability``：各示踪剂由逐只动物的值
    算得的逐只去除方向稳定性（差异结果没有 ``direction_stability`` 时使用）。

    方向不一致（conflicting）要求至少两种可辨识、且各自 |Hedges' g| ≥ ``flux_min_abs_g`` 的示踪剂指向相反；
    合并估计的异质性检验在绝对尺度上进行，不同示踪剂的通量尺度不同（设计文档 3.4），所以异质性 p 单独只记为
    ``tracer_heterogeneity``，不判为冲突，也不允许判为 established。"""
    tracer_rows = rows[rows["tracer"] != "pooled"]
    pooled = rows[rows["tracer"] == "pooled"]
    measured = tracer_rows[~tracer_rows["cell_state"].isin(["not_infused", "not_measured"])]
    has_pooled = not pooled.empty and pooled["cell_state"].iloc[0] != "not_measured"
    axis = Axis("flux", "not_measured", n_expected=len(tracer_rows), n_measured=len(measured), n_sources=1,
                basis=_describe(tracer_rows, "tracer"))
    if measured.empty and not has_pooled:
        axis.rule = "没有示踪数据"
        return axis
    identifiable = measured[~measured["cell_state"].isin(["not_identifiable", "insufficient_n"])]
    if measured.empty:
        identifiable = pooled if has_pooled and pooled["cell_state"].iloc[0] not in ("not_identifiable",
                                                                                     "insufficient_n") else pooled.iloc[0:0]
    if identifiable.empty:
        all_low_n = not measured.empty and measured["cell_state"].eq("insufficient_n").all()
        axis.state = "insufficient_n" if all_low_n or (measured.empty and pooled["cell_state"].iloc[0] ==
                                                        "insufficient_n") else "not_identifiable"
        axis.rule = "所有示踪剂都有一组少于 2 只" if axis.state == "insufficient_n" else "没有示踪剂能约束该反应"
        return axis
    head = pooled.iloc[0] if has_pooled else identifiable.iloc[0]
    signs = {_sign(v) for v in identifiable["log2fc"]} - {0}
    tracer_g = pd.to_numeric(identifiable["effect"], errors="coerce").abs()
    pointing = identifiable[tracer_g >= rules["flux_min_abs_g"]]
    pointing_signs = {_sign(v) for v in pointing["log2fc"]} - {0}
    axis.direction = _sign(head["log2fc"])
    axis.effect = float(head["log2fc"]) if pd.notna(head["log2fc"]) else np.nan
    if pd.notna(head.get("se")) and pd.notna(axis.effect):
        axis.ci_lower, axis.ci_upper = axis.effect - 1.96 * head["se"], axis.effect + 1.96 * head["se"]
    axis.sign_probability = _certainty(head)
    axis.pvalue_one_sided = one_sided(head.get("pvalue"), head["log2fc"], axis.direction)
    heterogeneity = pd.to_numeric(head.get("heterogeneity_pvalue"), errors="coerce")
    heterogeneous = pd.notna(heterogeneity) and heterogeneity < rules["heterogeneity_alpha"]
    if len(pointing) >= 2 and len(pointing_signs) > 1:
        axis.state, axis.direction = "mixed", 0
        axis.rule = "至少两种可辨识、效应明确（|g| ≥ {:g}）的示踪剂方向相反（conflicting）".format(rules["flux_min_abs_g"])
        return axis
    if heterogeneous:
        axis.flags.append("tracer_heterogeneity")
    n_min = pd.to_numeric(pd.Series([head.get("n_numerator"), head.get("n_denominator")]), errors="coerce").min()
    flags = _flags(pooled) | _flags(identifiable)
    single_subject = "single_subject_driven" in flags
    scale_drift = "flux_scale_drift" in flags
    if scale_drift:
        axis.flags.append("flux_scale_drift")
        if rules.get("flux_scale_drift_gate", True):
            axis.state = "undetermined"
            axis.rule = "整鼠通量存在系统性尺度漂移（flux_scale_drift），不作通量变化判定"
            return axis
    reported = pd.to_numeric(head.get("direction_stability"), errors="coerce")
    reported_loo = pd.to_numeric(head.get("pooled_loo_stability"), errors="coerce")
    pooled_computed = pd.to_numeric(stability.get("__pooled__", np.nan), errors="coerce") if stability else np.nan
    computed = [stability[t] for t in identifiable["tracer"] if stability and pd.notna(stability.get(t, np.nan))]
    if pd.notna(reported):
        stable = bool(reported >= rules["direction_stability_min"])
        if pd.notna(reported_loo):
            stable = stable and bool(reported_loo >= 1.0)
    elif pd.notna(pooled_computed):
        # 多示踪剂分层合并的逐只去除（LOO）稳定性（S2）：不再因最弱示踪剂单层反号否决合并结论
        stable = bool(pooled_computed >= 1.0)
    elif computed:
        # 单示踪剂或调用方仅传入分示踪剂 LOO 字典时的兜底口径
        stable = min(computed) >= 1.0
    else:
        stable = False
        axis.flags.append("stability_unknown")
    stable = stable and not single_subject
    fc_sign, g_sign = _sign(head.get("log2fc")), _sign(head.get("effect"))
    if len(identifiable) >= 2 and fc_sign != 0 and g_sign != 0 and fc_sign != g_sign:
        axis.flags.append("scale_direction_discordant")
    passes = pd.notna(head.get("padj")) and head["padj"] < rules["alpha"]
    pval = pd.to_numeric(head.get("pvalue"), errors="coerce")
    g = abs(pd.to_numeric(head.get("effect"), errors="coerce"))
    concordant_g_min = rules.get("flux_concordant_min_abs_g", RULES["flux_concordant_min_abs_g"])
    equiv_g_bound = rules.get("flux_equiv_g_bound", RULES["flux_equiv_g_bound"])
    perm_suggestive = (
        concordant_g_min is not None
        and float(concordant_g_min) < rules["flux_min_abs_g"]
        and len(identifiable) >= 2
        and len(signs) == 1
        and not heterogeneous
        and pd.notna(pval)
        and pval < rules["alpha"]
        and pd.notna(g)
        and g >= float(concordant_g_min)
        and stable
    )
    if len(identifiable) >= 2 and len(signs) == 1 and not heterogeneous and passes \
            and n_min >= rules["established_min_group"] and stable:
        axis.state, axis.rule = "changed", "至少两种可辨识示踪剂同向，合并估计通过检验（established）"
    elif pd.notna(g) and g >= rules["flux_min_abs_g"] and stable:
        axis.state = "suggestive"
        axis.rule = f"有可辨识示踪剂，|g| = {g:.2f} ≥ {rules['flux_min_abs_g']:g}，方向稳定，未达检验阈值"
    elif perm_suggestive:
        axis.state = "suggestive"
        axis.flags.append("multi_tracer_perm_suggestive")
        axis.rule = (f"多示踪剂同向且分层置换 p = {pval:.3g} < {rules['alpha']:g}，"
                     f"|g| = {g:.2f} ≥ {float(concordant_g_min):g}，方向稳定（未过 FDR）")
    else:
        # S5：通量级 TOST 等效性检验（严格把关：log2FC 90% CI 在 ±equivalence_bound 内、无异质性、非单鼠驱动、
        # 任意单示踪剂均无明确效应、且合并 |Hedges' g| 的 90% 上界 < flux_equiv_g_bound）
        smd_se, equiv_df = _flux_smd_se(head, identifiable)
        g_90_hi = float(g + float(student_t.ppf(0.95, equiv_df)) * smd_se) \
            if (pd.notna(g) and pd.notna(smd_se) and equiv_df >= 1) else np.nan
        no_tracer_signal = (
            not identifiable["cell_state"].isin(["changed", "nominal"]).any()
            and not (tracer_g >= rules["flux_min_abs_g"]).any()
        )
        flux_equivalent = (
            equiv_g_bound is not None
            and float(equiv_g_bound) > 0
            and str(head.get("cell_state") or "") == "equivalent"
            and not heterogeneous
            and not single_subject
            and no_tracer_signal
            and pd.notna(g_90_hi)
            and g_90_hi < float(equiv_g_bound)
        )
        if flux_equivalent:
            axis.state = "unchanged"
            axis.direction = 0
            axis.rule = (f"通量等效未变：log2FC 90% 区间在 ±{rules['equivalence_bound']:g} 内，"
                         f"且 |g| 90% 上界 {g_90_hi:.2f} < {float(equiv_g_bound):g}")
        else:
            axis.state = "undetermined"
            reasons: list[str] = []
            if pd.isna(g):
                reasons.append("|g| 缺失")
            elif g < rules["flux_min_abs_g"]:
                reasons.append(f"|g| = {g:.2f} < {rules['flux_min_abs_g']:g}")
            if single_subject:
                reasons.append("单只动物驱动")
            elif "stability_unknown" in axis.flags:
                reasons.append("逐只去除后的方向稳定性未知")
            elif not stable:
                if pd.notna(reported) and reported < rules["direction_stability_min"]:
                    reasons.append(f"方向稳定性 {reported:.2f} < {rules['direction_stability_min']:g}")
                elif pd.notna(reported_loo) and reported_loo < 1.0:
                    reasons.append(f"分层 LOO = {reported_loo:.2f} < 1")
                elif pd.notna(pooled_computed) and pooled_computed < 1.0:
                    reasons.append(f"分层 LOO = {pooled_computed:.2f} < 1")
                else:
                    reasons.append("逐只去除后方向不稳定")
            detail = "，".join(reasons) if reasons else (f"|g| = {g:.2f}" if pd.notna(g) else "")
            axis.rule = f"可辨识，未达 suggestive/changed（{detail}）" if detail else "可辨识，但效应小或不稳定"
    if single_subject:
        axis.flags.append("single_subject_driven")
    if expected:
        # 只比较样本量足够的示踪剂：一组少于 2 只是样本问题，不是可辨识性与预期不符
        evaluable = set(measured.loc[measured["cell_state"] != "insufficient_n", "tracer"])
        constrained, wanted = set(identifiable["tracer"]), set(expected)
        if (constrained - wanted) or ((wanted & evaluable) - constrained):
            axis.flags.append("tracer_expectation_mismatch")
            if axis.state in ("changed", "suggestive", "unchanged"):
                axis.state = "undetermined"
                axis.rule += "；可辨识的示踪剂与面板声明的预期不符，最高为 undetermined"
    return axis


# ---------------------------------------------------------------------------
# C 轴：酶容量
# ---------------------------------------------------------------------------

def _judge_members(rows: pd.DataFrame, rules: dict[str, Any]) -> Axis:
    """一层、一个来源的直接酶 → 容量的状态与方向。"""
    quantified = rows[rows["cell_state"].isin(MEASURED)]
    axis = Axis("capacity", _missing_state(rows), n_expected=len(rows), n_measured=len(quantified),
                basis=_describe(rows))
    if quantified.empty:
        return axis
    fc = pd.to_numeric(quantified["log2fc"], errors="coerce")
    states = quantified["cell_state"]
    if (quantified["member_relation"] == "complex_subunit").any():
        direction = _sign(fc.median())
        same = int((np.sign(fc) == direction).sum())
        axis.effect = float(fc.median())
        changed = states.eq("changed") & (np.sign(fc) == direction)
        if direction and same >= len(quantified) / 2 and changed.any():
            axis.state, axis.direction = "changed", direction
            axis.rule = f"复合体：{same}/{len(quantified)} 个已定量亚基同向，{int(changed.sum())} 个通过 FDR"
            _cap_if_cautioned(axis, quantified[changed])
        elif states.eq("equivalent").all():
            axis.state, axis.rule = "unchanged", "复合体：已定量亚基全部等效"
        elif states.isin(["changed", "nominal"]).any():
            axis.state, axis.direction, axis.rule = "suggestive", direction, "复合体：只有部分亚基有变化证据"
        else:
            axis.state, axis.rule = "undetermined", "复合体：测到，无变化证据"
    elif (quantified["member_relation"] == "pathway_step").any():
        # 集总多步/分枝或功能不可替代的诱导酶反应（如 FA_esterification_c、THERMOGENESIS_UCP、HMGCR）：
        # 不同步骤或不同机制的成员不互相替代，不按可替代同工酶的丰度占比压制
        changed = quantified[states == "changed"]
        dirs = {_sign(v) for v in changed["log2fc"]} - {0}
        if len(dirs) > 1:
            axis.state, axis.effect = "mixed", float(fc.median())
            axis.rule = "多步/集总反应：通过 FDR 的成员方向相反"
        elif dirs:
            direction = next(iter(dirs))
            axis.state, axis.direction = "changed", direction
            axis.effect = float(pd.to_numeric(changed["log2fc"], errors="coerce").median())
            names = "、".join(str(m) for m in changed["member_id"])
            axis.rule = f"多步/诱导反应：{names} 通过 FDR"
            _cap_if_cautioned(axis, changed)
        elif states.eq("equivalent").all():
            axis.state, axis.effect, axis.rule = "unchanged", float(fc.median()), "多步/集总反应：已定量成员全部等效"
        elif states.eq("nominal").any():
            nom = quantified[states == "nominal"]
            nom_dirs = {_sign(v) for v in nom["log2fc"]} - {0}
            if len(nom_dirs) == 1:
                axis.state, axis.direction = "suggestive", next(iter(nom_dirs))
                axis.effect = float(pd.to_numeric(nom["log2fc"], errors="coerce").median())
                axis.rule = "多步/集总反应：部分成员有名义 p"
            else:
                axis.state, axis.effect = "undetermined", float(fc.median())
                axis.rule = "多步/集总反应：名义 p 成员方向不一致"
        else:
            axis.state, axis.effect = "undetermined", float(fc.median())
            axis.rule = "多步/集总反应：测到，既不能判为变化，也不能判为等效"
    else:
        share = pd.to_numeric(quantified["abundance_share"], errors="coerce")
        share = share / share.sum() if share.notna().any() and share.sum() > 0 \
            else pd.Series(1.0 / len(quantified), index=quantified.index)
        share = share.fillna(0.0)
        if "active_abundance_share" in quantified.columns:
            active_share = pd.to_numeric(quantified["active_abundance_share"], errors="coerce")
            eff_share = np.maximum(share, active_share.fillna(share))
        else:
            eff_share = share
        axis.effect = float((share * fc.fillna(0.0)).sum())
        major = share.idxmax()
        # 主要成员：整体或活跃条件占比不低于 minor_share 的成员；成员很多、都低于该值时至少包括占比最高者
        majors = quantified[(eff_share >= rules["minor_share"]) | (quantified.index == major)]
        changed_majors = majors[majors["cell_state"] == "changed"]
        directions = {_sign(v) for v in changed_majors["log2fc"]} - {0}
        if len(directions) > 1:
            axis.state, axis.rule = "mixed", "主要同工酶方向相反（同工酶转换）"
            axis.flags.append("isozyme_switch")
        elif directions:
            axis.state, axis.direction = "changed", next(iter(directions))
            names = "、".join(f"{r.member_id}（占比 {eff_share[i]:.0%}）" for i, r in changed_majors.iterrows())
            axis.rule = f"主要成员 {names} 通过 FDR"
            _cap_if_cautioned(axis, changed_majors)
        elif majors["cell_state"].eq("equivalent").all():
            axis.state, axis.rule = "unchanged", "主要成员全部等效"
        elif states.eq("changed").any():
            changed_all = quantified[states == "changed"]
            c_dirs = {_sign(v) for v in changed_all["log2fc"]} - {0}
            c_dir = next(iter(c_dirs)) if len(c_dirs) == 1 else 0
            major_fcs = pd.to_numeric(majors["log2fc"], errors="coerce").fillna(0.0)
            major_opposing = bool(
                c_dir and (
                    (majors["cell_state"].isin(["changed", "nominal"]) & (np.sign(major_fcs) == -c_dir)).any()
                    or (_sign(quantified.at[major, "log2fc"]) == -c_dir and abs(float(major_fcs.loc[major])) >= 0.3)
                )
            )
            combined_share = float(eff_share[changed_all.index].sum())
            if c_dir and len(changed_all) >= 2 and not major_opposing and combined_share >= rules["minor_share"]:
                axis.state, axis.direction = "changed", c_dir
                axis.effect = float(pd.to_numeric(changed_all["log2fc"], errors="coerce").median())
                names = "、".join(str(m) for m in changed_all["member_id"])
                axis.rule = f"多个同向成员 {names}（合计占比 {combined_share:.0%}）通过 FDR"
                _cap_if_cautioned(axis, changed_all)
            else:
                axis.state, axis.rule = "undetermined", "只有次要同工酶通过 FDR"
                axis.flags.append("minor_isozyme_driven")
        elif majors["cell_state"].eq("nominal").any():
            axis.state, axis.direction = "suggestive", _sign(quantified.at[major, "log2fc"])
            axis.rule = "主要成员只有名义 p"
        else:
            axis.state, axis.rule = "undetermined", "测到，既不能判为变化，也不能判为等效"
    member = quantified.loc[quantified["cell_state"].isin(["changed", "nominal"])]
    axis.pvalue_one_sided = simes([one_sided(r.pvalue, r.log2fc, axis.direction) for r in member.itertuples()]) \
        if axis.direction else np.nan
    return axis


def capacity_axis(rows: pd.DataFrame, rules: dict[str, Any], sources: list[str]) -> Axis:
    """``rows``：该（反应, 组织, 比较）全部来源的成员行。"""
    direct = rows[rows["role"] == "direct_enzyme"]
    if direct[direct["cell_state"] != "no_knowledge"].empty:
        direct = rows[rows["role"].isin(["carrier", "direct_enzyme"])]
    direct = direct[(direct["cell_state"] != "no_knowledge") & (direct["member_relation"].fillna("") != "shared_subunit")]
    if direct.empty:
        return Axis("capacity", "no_knowledge", rule="知识表中没有直接酶或载体")
    judged = {(layer, source): _judge_members(direct[(direct["layer"] == layer) & (direct["source"] == source)], rules)
              for layer in ("proteomics", "transcriptomics") for source in sources
              if not direct[(direct["layer"] == layer) & (direct["source"] == source)].empty}
    protein = [judged[k] for k in judged if k[0] == "proteomics" and judged[k].n_measured]
    rna = [judged[k] for k in judged if k[0] == "transcriptomics" and judged[k].n_measured]
    if protein:
        axis = protein[0]
        axis.n_sources = len(protein)
        axis.replicated = len(protein) >= 2 and all(a.state == "changed" and a.direction == axis.direction
                                                    for a in protein)
        rna_changed = [a for a in rna if a.state == "changed"]
        if rna_changed and not (axis.state == "changed" and all(a.direction == axis.direction for a in rna_changed)):
            axis.flags.append("rna_protein_discordant")
    elif rna:
        changed = [a for a in rna if a.state == "changed"]
        directions = {a.direction for a in changed}
        lead = changed[0] if changed else rna[0]
        axis = Axis("capacity", lead.state, lead.direction, lead.effect, n_expected=sum(a.n_expected for a in rna),
                    n_measured=sum(a.n_measured for a in rna), n_sources=len(rna), basis=" | ".join(a.basis for a in rna),
                    pvalue_one_sided=lead.pvalue_one_sided, flags=["rna_only"] + lead.flags)
        if len(changed) >= 2 and len(directions) == 1:
            axis.state, axis.replicated = "changed", True
            axis.rule = f"没有蛋白数据；RNA 在 {len(changed)} 个来源中同向通过 FDR"
        elif len(directions) > 1:
            axis.state, axis.direction, axis.rule = "mixed", 0, "没有蛋白数据；RNA 在不同来源中方向相反"
        elif changed and len(rna) == 1 and rules.get("rna_single_source_changed", False):
            # 只有一个来源测到 RNA：按规则开关接受单一来源的 FDR；第二个来源测到而没有复现时不适用（走下一支）
            axis.state, axis.direction = "changed", lead.direction
            axis.flags.append("rna_single_source")
            axis.rule = "没有蛋白数据；RNA 只有一个来源，在该来源中通过 FDR（rna_single_source_changed）"
        elif changed or any(a.state == "suggestive" for a in rna):
            axis.state = "suggestive"
            axis.direction = lead.direction or next((a.direction for a in rna if a.direction), 0)
            axis.rule = "没有蛋白数据；RNA 只在一个来源中有变化证据"
        else:
            axis.rule = "没有蛋白数据；" + lead.rule
    else:
        axis = Axis("capacity", _missing_state(direct), n_expected=len(direct), basis=_describe(direct),
                    rule="直接酶在两层都没有测到")
    return axis


# ---------------------------------------------------------------------------
# R 轴：调控；M 轴：底物与产物
# ---------------------------------------------------------------------------

def regulation_axis(rows: pd.DataFrame, rules: dict[str, Any]) -> Axis:
    members = rows[rows["role"].isin(["regulator", "effector"]) & (rows["cell_state"] != "no_knowledge")]
    if members.empty:
        return Axis("regulation", "no_knowledge", rule="知识表中没有调控因子或效应物")
    signs = members.apply(_activity_sign, axis=1)
    signed = members[signs != 0]
    if signed.empty:
        return Axis("regulation", "no_knowledge", n_expected=len(members), basis=_describe(members),
                    rule="调控因子的方向未知或经由转录 / 信号起作用，不判定", flags=["unsigned_regulators"])
    signs = signs[signed.index]
    quantified = signed[signed["cell_state"].isin(MEASURED)]
    axis = Axis("regulation", _missing_state(signed), n_expected=len(signed), n_measured=len(quantified),
                n_sources=signed["source"].nunique(), basis=_describe(signed))
    if quantified.empty:
        axis.rule = "带符号的调控因子都没有测到"
        return axis
    implied = np.sign(pd.to_numeric(quantified["log2fc"], errors="coerce")) * signs[quantified.index]
    evidence = quantified["cell_state"].isin(["changed", "nominal"])
    changed = quantified["cell_state"] == "changed"
    if ((quantified["regulator_effect"] == "inhibits") & evidence
            & (pd.to_numeric(quantified["log2fc"], errors="coerce") > 0)).any():
        axis.flags.append("regulator_inhibitory_up")
    pointing = set(implied[changed].dropna()) - {0} or set(implied[evidence].dropna()) - {0}
    if len(pointing) > 1:
        axis.state, axis.rule = "mixed", "调控因子对活性的指向相反"
    elif pointing:
        axis.direction = int(next(iter(pointing)))
        activity_data = (changed & (quantified["role"] == "effector")).any()
        axis.state = "changed" if activity_data else "suggestive"
        axis.rule = "测得的别构效应物有变化" if activity_data else \
            ("由调控因子的量推断，最多到 suggestive" if changed.any() else "调控因子只有名义 p")
        if axis.state == "changed":
            eff_changed = quantified[changed & (quantified["role"] == "effector") & (implied == axis.direction)]
            _cap_if_cautioned(axis, eff_changed, _LAYER_QC_CAUTION)
        support = quantified[evidence & (implied == axis.direction)]
        axis.effect = float((pd.to_numeric(support["log2fc"], errors="coerce") * signs[support.index]).median())
        axis.pvalue_one_sided = simes([one_sided(r.pvalue, r.log2fc, axis.direction * signs[i])
                                       for i, r in zip(support.index, support.itertuples())])
    elif quantified["cell_state"].eq("equivalent").all():
        axis.state, axis.rule = "unchanged", "带符号的调控因子全部等效"
    else:
        axis.state, axis.rule = "undetermined", "测到，无变化证据"
    return axis


def pool_axis(rows: pd.DataFrame, rules: dict[str, Any], direction_evidence: str) -> Axis:
    members = rows[((rows["layer"] == "metabolomics") & rows["role"].isin(["substrate", "product"]))
                   | ((rows["layer"] == "derived") & (rows["role"] == "step_ratio"))]
    members = members[members["cell_state"] != "no_knowledge"]
    if members.empty:
        return Axis("pool", "no_knowledge", rule="知识表中没有底物或产物")
    quantified = members[members["cell_state"].isin(MEASURED)]
    axis = Axis("pool", _missing_state(members), n_expected=len(members), n_measured=len(quantified),
                n_sources=members["source"].nunique(), basis=_describe(members))
    if quantified.empty:
        axis.rule = "底物与产物都没有测到"
        return axis
    factor = quantified["role"].map({"product": 1, "step_ratio": 1, "substrate": -1})
    implied = np.sign(pd.to_numeric(quantified["log2fc"], errors="coerce")) * factor
    evidence = quantified["cell_state"].isin(["changed", "nominal"])
    changed = quantified["cell_state"] == "changed"
    pointing = set(implied[changed].dropna()) - {0} or set(implied[evidence].dropna()) - {0}
    if direction_evidence != "directional":
        axis.state = "changed" if changed.any() else ("suggestive" if evidence.any() else "undetermined")
        axis.rule, axis.flags = "只表示容量的反应：池的变化不提示净方向", ["capacity_reaction"]
        if axis.state == "changed":
            _cap_if_cautioned(axis, quantified[changed], _LAYER_QC_CAUTION)
        return axis
    if len(pointing) > 1:
        axis.state, axis.rule = "mixed", "底物与产物的指向相反"
    elif pointing:
        axis.direction = int(next(iter(pointing)))
        axis.state = "changed" if changed.any() else "suggestive"
        axis.rule = "产物升 / 底物降指向通量升，反之指向通量降"
        if axis.state == "changed":
            _cap_if_cautioned(axis, quantified[changed & (implied == axis.direction)], _LAYER_QC_CAUTION)
        support = quantified[evidence & (implied == axis.direction)]
        axis.effect = float((pd.to_numeric(support["log2fc"], errors="coerce") * factor[support.index]).median())
        axis.pvalue_one_sided = simes([one_sided(r.pvalue, r.log2fc, axis.direction * factor[i])
                                       for i, r in zip(support.index, support.itertuples())])
    elif quantified["cell_state"].eq("equivalent").all():
        axis.state, axis.rule = "unchanged", "底物与产物全部等效"
    else:
        axis.state, axis.rule = "undetermined", "测到，无变化证据"
    return axis


# ---------------------------------------------------------------------------
# 判定表、调控系数、合取检验
# ---------------------------------------------------------------------------

def decide(f: Axis, c: Axis, r: Axis, m: Axis, rules: dict[str, Any] | None = None,
           semantics: str = "directional") -> dict[str, Any]:
    """判定表（3.5）：结论等级、方向、同向的非通量轴数、反向的轴、现象标签。"""
    if semantics == "not_assessable":
        # 虚拟/血浆/输入输出等不可评估反应（如 BIOMASS_REACTION_*、Salvage_c、*_protein_input）：
        # 保留各轴明细供核对，但绝不允许升为 L1 / L2 / L3 / LX
        return {
            "flux_claim": "L0",
            "claim_direction": 0,
            "n_concordant_axes": 0,
            "concordant": [],
            "discordant_axes": "",
            "phenomena": "not_assessable_reaction",
        }
    rules = {**RULES, **(rules or {})}
    vetoing = ("changed",) if rules["opposition_min_state"] == "changed" else ("changed", "suggestive")
    others = {"capacity": c, "regulation": r, "pool": m}
    if f.state == "mixed":
        claim, direction = "LX", 0
    elif f.state == "changed":
        claim, direction = "L3", f.direction
    elif f.state == "suggestive":
        claim, direction = "L2", f.direction
    elif f.state == "unchanged":
        # 实测通量已通过等效性检验（未变化）：不能再由酶容量外推为 L1 通量变化
        claim, direction = "L0", 0
    else:
        robust = c.state == "changed" and not (set(_CAUTION_FLAGS) & set(c.flags))
        opposite = any(ax.state in vetoing and ax.direction == -c.direction for ax in (r, m))
        claim, direction = ("L1", c.direction) if robust and not opposite else ("L0", 0)
    concordant = [name for name, ax in others.items() if direction and ax.pointing and ax.direction == direction]
    reference = direction or (c.direction if c.pointing else 0)
    discordant = [name for name, ax in others.items() if reference and ax.pointing and ax.direction == -reference]
    phenomena = []
    if c.pointing and c.state == "changed":
        phenomena.append("capacity_up" if c.direction > 0 else "capacity_down")
    if f.state == "unchanged" and c.pointing:
        phenomena.append("flux_unchanged_despite_capacity")
    if c.state == "mixed" or "isozyme_switch" in c.flags:
        phenomena.append("isozyme_switch")
    if "regulator_inhibitory_up" in r.flags:
        phenomena.append("regulator_inhibitory_up")
    if c.pointing and m.pointing and c.direction == -m.direction:
        phenomena.append("capacity_pool_opposed")
    all_axis_flags = c.flags + f.flags + r.flags + m.flags
    for flag in (
        "rna_protein_discordant", "tracer_expectation_mismatch", "rna_single_source",
        "flux_scale_drift", "foreign_tissue_marker", "qc_caution",
    ):
        if flag in all_axis_flags:
            phenomena.append(flag)
    return {"flux_claim": claim, "claim_direction": direction, "n_concordant_axes": len(concordant),
            "concordant": concordant, "discordant_axes": ";".join(discordant), "phenomena": ";".join(phenomena)}


def mechanism(claim: str, f: Axis, c: Axis, rho: tuple[float, float]) -> str:
    """机制类型（只在 L2 及以上给出）：按 C 轴与调控系数 ρ_h 的区间归类。"""
    if claim not in ("L2", "L3"):
        return ""
    if c.state == "unchanged":
        return "metabolic_dominant"
    lower, upper = rho
    if c.state not in ("changed", "suggestive") or pd.isna(lower):
        return "unresolved"
    if upper < 0:
        return "antagonistic"
    if lower <= 1 <= upper:
        return "enzyme_dominant"
    if lower <= 0 <= upper:
        return "metabolic_dominant"
    if 0 < lower and upper < 1:
        return "both"
    return "unresolved"


def regulation_coefficient(flux: pd.DataFrame, enzyme: pd.DataFrame, sides: dict[tuple, str], rules: dict[str, Any],
                           weights: dict[str, float]) -> tuple[float, float]:
    """层级调控系数 ρ_h = Δln E / Δln J 的分层自助区间（3.7）：通量动物与组学动物不是同一批，各层在每一组内独立重抽。

    ``flux``：一种示踪剂的逐只动物通量（线性，要求两组均值为正）；``enzyme``：直接酶成员的逐只动物 log2 值
    （``feature_id`` 区分成员，按 ``weights`` 加权）；``sides``：（层, 样本）→ numerator / denominator。"""
    rng = np.random.default_rng(int(rules["seed"]))

    def groups(frame: pd.DataFrame, layer: str) -> tuple[list, list]:
        side = frame["sample_id"].map(lambda s: sides.get((layer, s), ""))
        return list(frame.loc[side == "numerator", "sample_id"].unique()), \
            list(frame.loc[side == "denominator", "sample_id"].unique())

    f_num, f_den = groups(flux, "mfa")
    e_layer = enzyme["layer"].iloc[0] if not enzyme.empty else ""
    e_num, e_den = groups(enzyme, e_layer)
    if min(len(f_num), len(f_den), len(e_num), len(e_den)) < 2:
        return np.nan, np.nan
    fv = flux.groupby("sample_id")["value"].mean()
    ev = enzyme.pivot_table(index="sample_id", columns="feature_id", values="value", aggfunc="mean")
    w = pd.Series({k: v for k, v in weights.items() if k in ev.columns}, dtype=float)
    if w.empty:
        return np.nan, np.nan
    w = (w / w.sum()).to_numpy()
    columns = [k for k in weights if k in ev.columns]
    b = int(rules["bootstrap"])
    # 向量化的分层自助：每次在各层、各组内有放回地重抽动物
    jn = fv[f_num].to_numpy()[rng.integers(0, len(f_num), (b, len(f_num)))].mean(axis=1)
    jd = fv[f_den].to_numpy()[rng.integers(0, len(f_den), (b, len(f_den)))].mean(axis=1)
    en = ev.loc[e_num, columns].to_numpy()[rng.integers(0, len(e_num), (b, len(e_num)))].mean(axis=1)
    ed = ev.loc[e_den, columns].to_numpy()[rng.integers(0, len(e_den), (b, len(e_den)))].mean(axis=1)
    usable = (jn > 0) & (jd > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        d_ln_j = np.log(np.where(usable, jn, np.nan)) - np.log(np.where(usable, jd, np.nan))
        d_ln_e = np.log(2.0) * np.nansum((en - ed) * w, axis=1)
        rhos = d_ln_e / d_ln_j
    rhos = rhos[np.isfinite(rhos)]
    if len(rhos) < 0.5 * b:
        return np.nan, np.nan
    return float(np.percentile(rhos, 2.5)), float(np.percentile(rhos, 97.5))


# ---------------------------------------------------------------------------
# 汇总
# ---------------------------------------------------------------------------

def _completeness(rows: pd.DataFrame) -> tuple[float, str]:
    relevant = rows[rows["role"].isin(_RELEVANT) & (rows["cell_state"] != "no_knowledge")]
    if relevant.empty:
        return np.nan, ""
    measured = ~relevant["cell_state"].isin(["not_measured", "not_detected", "not_infused"])
    by_layer = "; ".join(f"{layer} {int(measured[part.index].sum())}/{len(part)}"
                         for layer, part in relevant.groupby("layer", sort=True))
    return float(measured.mean()), by_layer


def assess(members: pd.DataFrame, primary: str, rules: dict[str, Any] | None = None,
           expected_tracers: dict[str, list[str]] | None = None, samples: pd.DataFrame | None = None,
           sides: pd.DataFrame | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """``members``：dossier_members（全部来源）；``primary``：本实验的来源标签。返回（axis_assessment, reaction_assessment）。"""
    rules = {**RULES, **(rules or {})}
    expected_tracers = expected_tracers or {}
    # pandas 3 的字符串列默认由 pyarrow 存放，逐格的小筛选每次都要复制（约占九成时间）：先转成普通 object 数组
    members = members.copy()
    for column in members.columns:
        if pd.api.types.is_string_dtype(members[column]) and not pd.api.types.is_object_dtype(members[column]):
            members[column] = members[column].astype(object)
    sources = [primary] + sorted(set(members["source"]) - {primary})
    base = members[members["source"] == primary]
    testable = _testable(base, expected_tracers)
    # 各比较的分子 / 分母样本（本实验）：（层, 样本）→ numerator / denominator
    sides_by_contrast: dict[str, dict[tuple, str]] = {}
    if sides is not None and not sides.empty:
        own_sides = sides[sides["source"] == primary]
        for contrast, part in own_sides.groupby("contrast_id"):
            sides_by_contrast[str(contrast)] = {(s.layer, s.sample_id): s.side for s in part.itertuples()}
    # 先按格子分组一次（逐格布尔筛选在十万行的成员表上是平方级的）
    flux_samples = {key: part for key, part in samples[(samples["source"] == primary) & (samples["layer"] == "mfa")]
                    .groupby(["feature_id", "tissue"])} if samples is not None and not samples.empty else {}
    cells = {key: part for key, part in members.groupby(["flux_id", "tissue", "contrast_id"], sort=False)}
    axis_rows, reaction_rows = [], []
    for (flux_id, tissue, contrast), _ in base.groupby(["flux_id", "tissue", "contrast_id"], sort=True):
        rows = cells[(flux_id, tissue, contrast)]
        own = rows[rows["source"] == primary]
        semantics = str(own["direction_evidence"].dropna().iloc[0]) if own["direction_evidence"].notna().any() else ""
        side_map = sides_by_contrast.get(str(contrast), {})
        mfa_own = own[own["layer"] == "mfa"]
        ident_tracers = set(
            mfa_own.loc[
                (mfa_own["tracer"] != "pooled")
                & ~mfa_own["cell_state"].isin(["not_infused", "not_measured", "not_identifiable", "insufficient_n"]),
                "tracer",
            ].astype(str)
        )
        stability = _flux_stability(flux_samples.get((flux_id, tissue)), side_map, ident_tracers)
        f = flux_axis(mfa_own, rules, expected_tracers.get(flux_id), stability)
        c = capacity_axis(rows, rules, sources)
        r = regulation_axis(rows, rules)
        m = pool_axis(rows, rules, semantics)
        decision = decide(f, c, r, m, rules, semantics=semantics)
        rho = (np.nan, np.nan)
        if decision["flux_claim"] in ("L2", "L3") and samples is not None and c.state in ("changed", "suggestive"):
            rho = _rho(rows, samples, side_map, rules, tissue, f)
        # 联合结论的分量：给出方向的通量轴（L2 / L3），以及与结论同向的非通量轴
        components = ([("flux", f)] if f.pointing else []) + [(n, ax) for n, ax in (("capacity", c), ("regulation", r),
                                                                                     ("pool", m))
                                                              if n in decision["concordant"]]
        word = {1: "up", -1: "down"}
        statement = " & ".join(f"{name} {word[decision['claim_direction'] or ax.direction]}" for name, ax in components) \
            if decision["claim_direction"] and components else ""
        pvals = [ax.pvalue_one_sided for _, ax in components]
        conj = float(max(pvals)) if statement and all(pd.notna(p) for p in pvals) else np.nan
        completeness, by_layer = _completeness(rows)
        reaction_rows.append({
            "flux_id": flux_id, "tissue": tissue, "contrast_id": contrast, "reaction_class": semantics,
            "flux_claim": decision["flux_claim"], "claim_direction": decision["claim_direction"],
            "n_concordant_axes": decision["n_concordant_axes"], "discordant_axes": decision["discordant_axes"],
            "mechanism_class": mechanism(decision["flux_claim"], f, c, rho), "rho_h_lower": rho[0],
            "rho_h_upper": rho[1], "conjunction_statement": statement, "conjunction_pvalue": conj,
            "any_pvalue": simes(pvals) if statement else np.nan, "completeness": completeness,
            "completeness_by_layer": by_layer, "testable_by": ";".join(testable.get((flux_id, tissue), [])),
            "phenomena": decision["phenomena"], "_effect": _lead_effect(decision["flux_claim"], f, c),
            "_replicated": c.replicated})
        for ax in (f, c, r, m):
            axis_rows.append({"flux_id": flux_id, "tissue": tissue, "contrast_id": contrast, "axis": ax.axis,
                              "direction": ax.direction, "state": ax.state, "effect": ax.effect,
                              "ci_lower": ax.ci_lower, "ci_upper": ax.ci_upper, "sign_probability": ax.sign_probability,
                              "pvalue_one_sided": ax.pvalue_one_sided, "n_expected": ax.n_expected,
                              "n_measured": ax.n_measured, "n_sources": ax.n_sources, "replicated": ax.replicated,
                              "basis": ax.basis, "rule": ax.rule, "flags": ";".join(dict.fromkeys(ax.flags))})
    axes = pd.DataFrame(axis_rows, columns=AXIS_COLUMNS)
    reactions = pd.DataFrame(reaction_rows)
    if reactions.empty:
        return axes, pd.DataFrame(columns=REACTION_COLUMNS)
    reactions["conjunction_padj"] = np.nan
    for key, part in reactions.groupby(["tissue", "contrast_id", "conjunction_statement"]):
        if key[2]:
            reactions.loc[part.index, "conjunction_padj"] = bh(part["conjunction_pvalue"])
    reactions["rank"] = 0
    for _, part in reactions.groupby(["tissue", "contrast_id"]):
        keys = pd.DataFrame({"claim": part["flux_claim"].map(CLAIM_ORDER), "axes": part["n_concordant_axes"],
                             "replicated": part["_replicated"].astype(int),
                             "completeness": part["completeness"].fillna(0),
                             "effect": part["_effect"].abs().fillna(0)})
        ordered = keys.sort_values(["claim", "axes", "replicated", "completeness", "effect"], ascending=False,
                                   kind="stable")
        reactions.loc[ordered.index, "rank"] = np.arange(1, len(ordered) + 1)
    return axes, reactions[REACTION_COLUMNS]


def _flux_stability(part: pd.DataFrame | None, sides: dict[tuple, str],
                    identifiable_tracers: set[str] | None = None) -> dict[str, float]:
    """各示踪剂的逐只去除方向稳定性，以及跨可辨识示踪剂的分层合并逐只去除稳定性（``__pooled__``）。"""
    if part is None or part.empty or not sides:
        return {}
    out: dict[str, float] = {}
    pooled_strata: list[tuple[np.ndarray, np.ndarray]] = []
    for tracer, values in part.groupby("tracer"):
        side = values["sample_id"].map(lambda s: sides.get(("mfa", s), ""))
        num = values.loc[side == "numerator", "value"].to_numpy()
        den = values.loc[side == "denominator", "value"].to_numpy()
        out[str(tracer)] = loo_stability(num, den)
        if identifiable_tracers is None or str(tracer) in identifiable_tracers:
            pooled_strata.append((num, den))
    if len(pooled_strata) >= 2:
        pooled_loo = stratified_pooled_loo(pooled_strata)
        if pd.notna(pooled_loo):
            out["__pooled__"] = pooled_loo
    return out


def _lead_effect(claim: str, f: Axis, c: Axis) -> float:
    """排序的最后一个键：最高一级证据的效应量（L2 / L3 取通量，L1 取酶容量）。"""
    if claim in ("L2", "L3", "LX"):
        return f.effect
    return c.effect if pd.notna(c.effect) else np.nan


def _testable(base: pd.DataFrame, expected: dict[str, list[str]]) -> dict[tuple, list[str]]:
    """能检验该反应的示踪剂：该组织中（任一比较）可辨识的示踪剂，以及面板预期但尚未输注的示踪剂。"""
    flux = base[(base["layer"] == "mfa") & (base["tracer"] != "pooled")]
    out: dict[tuple, list[str]] = {}
    for (flux_id, tissue), part in flux.groupby(["flux_id", "tissue"]):
        ok = part[~part["cell_state"].isin(["not_identifiable", "insufficient_n", "not_infused", "not_measured"])]
        found = sorted(set(ok["tracer"]))
        # 尚未输注：该组织在任何比较中都没有这种示踪剂的数据（某个比较缺少某组时也会记 not_infused，不算）
        infused = set(part.loc[part["cell_state"] != "not_infused", "tracer"])
        pending = [t for t in expected.get(flux_id, []) if t not in infused]
        out[(flux_id, tissue)] = found + [f"{t}(not infused)" for t in pending if t not in found]
    return out


def _rho(rows: pd.DataFrame, samples: pd.DataFrame, sides: dict, rules: dict[str, Any], tissue: str,
         f: Axis) -> tuple[float, float]:
    """取最可辨识的一种示踪剂的逐只动物通量与蛋白（没有时 RNA）直接酶，计算 ρ_h 的区间。"""
    flux_rows = rows[(rows["layer"] == "mfa") & (rows["tracer"] != "pooled")
                     & ~rows["cell_state"].isin(["not_identifiable", "insufficient_n", "not_infused", "not_measured"])]
    if flux_rows.empty:
        return np.nan, np.nan
    rank = {"well_constrained": 0, "moderate": 1}
    lead = flux_rows.assign(_r=flux_rows["identifiability"].map(rank).fillna(9)).sort_values("_r").iloc[0]
    flux = samples[(samples["layer"] == "mfa") & (samples["tissue"] == tissue) & (samples["source"] == lead["source"])
                   & (samples["feature_id"] == lead["flux_id"]) & (samples["tracer"] == lead["tracer"])]
    for layer in ("proteomics", "transcriptomics"):
        enzymes = rows[(rows["layer"] == layer) & (rows["role"] == "direct_enzyme") & (rows["source"] == lead["source"])
                       & rows["cell_state"].isin(MEASURED) & (rows["member_relation"].fillna("") != "shared_subunit")]
        if enzymes.empty:
            continue
        shares = pd.to_numeric(enzymes["abundance_share"], errors="coerce").fillna(1.0)
        weights = dict(zip(enzymes["feature_id"].astype(str), shares))
        values = samples[(samples["layer"] == layer) & (samples["tissue"] == tissue)
                         & (samples["source"] == lead["source"]) & samples["feature_id"].isin(list(weights))]
        return regulation_coefficient(flux, values, sides, rules, weights)
    return np.nan, np.nan
