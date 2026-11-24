#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
反应档案（``reaction_dossier``，G13 第一部分）：不以通量为锚，对每个（反应, 组织, 比较）列出各层的每个成员——
酶与转运体（RNA、蛋白）、调控因子、底物与产物、派生比值、别构效应物、各示踪剂的通量——以及它在数据中的
统计量和格子状态。测到的、没测到的、知识表里没有的都各占一行，状态互不混淆（设计见
``docs/arch_design/omicsweave_reaction_dossier_and_assessment.md`` 第 4.1 与第 6 节）。

成员集合、匹配与入选规则（``selected``）都来自 :meth:`FluxEvidenceBuilder.member_matches`，与 ``flux_anchored``
相同；本模式不改变 ``flux_anchored`` 的任何输出。统计量从源结果包的统一结果表按特征 ID 回查。

格子状态（``cell_state``）：

- ``changed``：padj < ``alpha``；``nominal``：p < ``alpha``，padj 未过；
- ``equivalent``：90% 区间（log2FC ± t₀.₉₅·SE，自由度取 ``residual_df``，没有时 n1 + n2 − 2）落在
  ±``equivalence_bound`` 内（双单侧检验）。RNA 的 apeGLM 收缩估计偏向 0，不做等效判定；
- ``no_evidence``：测到，既未判为变化也未判为等效；
- ``not_detected``：知识表中有该成员、该层有数据，但数据中没有这个特征；
- ``not_measured``：该组织、该比较没有这一层的数据；
- ``qc_blocked``：该层在此（组织, 比较）上被 QC 判为不可评估（:mod:`omicsweave.core.reliability`）；
- ``no_knowledge``：知识表中这个角色没有成员（``member_id`` 为空）；
- 通量另有 ``not_infused``（参数 ``tracers`` 中声明、但该组织没有这种示踪剂的数据）、``insufficient_n``
  （较小组少于 2 只）、``not_identifiable``（``poorly_constrained`` 或 ``bound_pinned``）。

参数（``integration.parameters.reaction_dossier``）：``tracers``（应出现的示踪剂）、``reactions`` / ``tissues``
（只列这些；默认全部）、``label``（本实验的来源标签，默认项目 ID）、``rules``（``alpha``、``equivalence_bound``）、
``sources``（其他来源，见 :mod:`omicsweave.integration.dossier.sources`）、``panels``（面板声明，相对实验目录，见
:mod:`omicsweave.integration.dossier.panels`）、``figures``（``formats``、``max_reactions``）；
``pval_cutoff``、``log2fc_cutoff``、``omics_cutoffs``、``fallback_mode``、``isozyme``、``tissue_map``、``qc_imbalance``
没有写时沿用 ``flux_anchored`` 的参数，使入选规则一致。

输出（``integration/reaction_dossier/``）：

- ``tables/dossier_members.parquet``（列见 :data:`COLUMNS`；``source`` 区分来源）；
- ``tables/dossier_samples.parquet``（逐只动物值）与 ``tables/dossier_sample_contrasts.parquet``（各比较的分子 / 分母样本）；
- 有面板时：``tables/panel_elements.parquet``，以及 ``figures/<面板>/<组织>/`` 下的 ``matrix_<比较>``（图 A 与同名的
  TSV / XLSX 导出表）、``dossier_<反应>_<比较>``（图 B）、``tracer_coverage_<比较>``（图 C）；
- ``tables/axis_assessment.parquet`` 与 ``tables/reaction_assessment.parquet``：分轴综合评估（S3，见
  :mod:`omicsweave.integration.dossier.assessment`）；图 A 的评估块即这两张表；
- ``evaluation/assessment_calibration.parquet``：判定规则的校准（S4，E8，见
  :mod:`omicsweave.integration.dossier.calibration`；参数 ``calibration``：``decoys`` 诱饵重复次数，默认 0 不运行；
  ``seed``；``invariance`` 在实际数据上做不变性检查）；面板声明 ``controls`` 时含对照集的符合率；图 A 的图注写出诱饵对照；
- 面板声明 ``branch_nodes`` / ``signed_modules`` 时：``tables/node_split.parquet``、``tables/node_split_samples.parquet``
  与 ``tables/pathway_tests.parquet``（S5，G16：节点分流比、带符号的模块旋转检验、交叉点分析，见
  :mod:`omicsweave.integration.dossier.pathway`；参数 ``pathway``：``rotations``、``seed``），图 D ``node_split_<比较>``
  与同名导出表；
- ``qc/contrast_reliability.parquet``（各来源的可靠性表）。
"""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import t as student_t

from omicsweave.contract.dataset import META_FILE, read_dataset
from omicsweave.contract.export import write_excel_tsv
from omicsweave.core.model import IntegrationMode, Omics
from omicsweave.core.reliability import BLOCKED, QC_DIR, RELIABILITY_FILE, build_reliability
from omicsweave.core.run import RunContext, RunError
from omicsweave.core.tissues import merge_tissue_aliases, resolve_tissue, validate_tissue_aliases
from omicsweave.stats.data import analysis_values

from ..dossier import calibration as cal
from ..dossier import pathway as pw
from ..dossier.assessment import RULES as ASSESSMENT_RULES, assess
from ..evaluation.decoys import _Pool as DecoyPool
from ..mappings.feature_resolver import _canonical_gene_symbol, match_gene_features, match_protein_features
from ..mappings.model_package import DecoyPackage
from ..dossier.panels import Panel, PanelError, load_panels, panel_elements
from ..dossier.samples import _dataset_dir, contrast_samples, dossier_samples
from ..dossier.sources import DossierSource, other_source, primary_source
from ..dossier.view import export_matrix, matrix_view, member_table
from ..integrate import IntegrationRun, mode_dir
from ..mappings import FLUX_GENE_MAPPING, FLUX_METABOLITE_MAPPING, REACTION_ENZYME_MAPPING
from ..mappings.network_layout import LayoutError
from ..mappings.reaction_aliases import assessment_universe
from ..pathway_map import DEFAULTS as PATHWAY_MAP_DEFAULTS, render_maps, settings as pathway_map_settings
from ..mappings.reaction_members import ACTIVITY_MODES
from ..methods.flux_evidence import FluxEvidenceBuilder
from . import flux_anchored as fa

MODE = IntegrationMode.REACTION_DOSSIER
# 与 flux_anchored 共用（没有写时沿用其参数）
SHARED = ("pval_cutoff", "log2fc_cutoff", "omics_cutoffs", "fallback_mode", "isozyme", "tissue_map", "qc_imbalance")
DEFAULTS: dict[str, Any] = {**{key: fa.DEFAULTS[key] for key in SHARED}, "tracers": [], "reactions": [],
                            "tissues": [], "label": "",
                            "rules": {k: v for k, v in ASSESSMENT_RULES.items() if k != "minor_share"},
                            "sources": [], "panels": "", "figures": {"formats": ["pdf", "png"], "max_reactions": 15},
                            # E8：诱饵重复次数（0 为不运行，每次约等于再建一遍成员行）、随机种子、是否在实际数据上做不变性检查
                            "calibration": {"decoys": 0, "seed": 20261117, "invariance": True},
                            # S5：模块旋转检验的旋转次数与随机种子
                            "pathway": {"rotations": 9999, "seed": 20261117},
                            # G18：分组织代谢通路总图（图 E），见 omicsweave.integration.pathway_map
                            "pathway_map": dict(PATHWAY_MAP_DEFAULTS)}
DOSSIER_TABLE = "dossier_members.parquet"
SAMPLES_TABLE = "dossier_samples.parquet"
SIDES_TABLE = "dossier_sample_contrasts.parquet"
ELEMENTS_TABLE = "panel_elements.parquet"
AXIS_TABLE = "axis_assessment.parquet"
ASSESSMENT_TABLE = "reaction_assessment.parquet"
CALIBRATION_TABLE = "assessment_calibration.parquet"
NODE_SPLIT_TABLE = "node_split.parquet"
NODE_SPLIT_SAMPLES = "node_split_samples.parquet"
PATHWAY_TESTS_TABLE = "pathway_tests.parquet"
INVARIANCE_FALLBACK = 30          # 没有面板时，不变性检查取前 30 个反应
Z95 = 1.959963984540054              # ci_lower / ci_upper 的 95% 区间

COLUMNS = [
    "flux_id", "tissue", "contrast_id", "source", "layer", "role", "member_id", "feature_id", "feature_name", "tracer",
    "member_relation", "regulator_effect", "regulator_mode", "direction_member", "mapping_method",
    "log2fc", "log2fc_unshrunk", "effect", "effect_type", "se", "ci_lower", "ci_upper", "pvalue", "padj", "lfsr",
    "n_numerator", "n_denominator", "mean_numerator", "mean_denominator", "abundance_share", "active_abundance_share",
    "cell_state", "state_reason", "selected", "identifiability", "frac_at_bound", "evidence_level", "qc_status",
    "flags", "implied_direction", "direction_evidence", "in_flux_anchored",
    # 通量合并估计的附加读数（S3 的 F 轴与交互诊断用到）：示踪剂间异质性、分层方向稳定性与分层 LOO、
    # 合并权重（标准化 / 绝对尺度）、标准化标准误与方差放大诊断（smd_se / g_total_var）、全局校正 p 值与示踪剂证据类别
    "heterogeneity_pvalue", "direction_stability", "pooled_loo_stability", "tracer_evidence",
    "smd_se", "smd_weights", "tracer_weights", "g_total_var", "smd_weighting", "padj_global", "fdr_scope",
]
_LONG_STATS = {"log2fc": "log2fc", "log2fc_unshrunk": "log2fc_unshrunk", "effect_size": "effect",
               "effect_size_type": "effect_type", "lfc_se": "se", "pvalue": "pvalue", "padj": "padj", "lfsr": "lfsr",
               "n_numerator": "n_numerator", "n_denominator": "n_denominator", "mean_numerator": "mean_numerator",
               "mean_denominator": "mean_denominator", "evidence_level": "evidence_level", "flags": "data_flags",
               "lfc_shrinkage": "lfc_shrinkage", "residual_df": "residual_df", "feature_name": "long_feature_name"}
_FRAME_ID = {"transcriptomics": ("feature_id", "gene", "base_mean"), "proteomics": ("protein", "protein", "abundance"),
             "metabolomics": ("metabolite_id", "compound", "abundance")}
_UNIDENTIFIABLE = ("poorly_constrained", "bound_pinned")
# 没有成员时记 no_knowledge 的角色（每层）
_EXPECTED_ROLES = {"transcriptomics": ("direct_enzyme", "regulator"), "proteomics": ("direct_enzyme", "regulator"),
                   "metabolomics": ("substrate", "product", "effector")}


def _parameters(ctx: RunContext, parameters: dict[str, Any]) -> dict[str, Any]:
    unknown = sorted(set(parameters) - set(DEFAULTS))
    if unknown:
        raise RunError(f"integration.parameters.{MODE.value} 中有无法识别的参数：{'、'.join(unknown)}"
                       f"（可用：{'、'.join(DEFAULTS)}）")
    configured = dict(ctx.experiment.config.integration.parameters.get(IntegrationMode.FLUX_ANCHORED.value, {}) or {})
    inherited = {key: configured[key] for key in SHARED if key in configured and key not in parameters}
    params = {**DEFAULTS, **inherited, **parameters}
    bad_rules = sorted(set(parameters.get("rules") or {}) - set(DEFAULTS["rules"]))
    if bad_rules:
        raise RunError(f"integration.parameters.{MODE.value}.rules 中有无法识别的阈值：{'、'.join(bad_rules)}"
                       f"（可用：{'、'.join(DEFAULTS['rules'])}）")
    params["rules"] = {**DEFAULTS["rules"], **(parameters.get("rules") or {})}
    for bool_rule in ("rna_single_source_changed", "flux_scale_drift_gate"):
        if not isinstance(params["rules"][bool_rule], bool):
            raise RunError(f"integration.parameters.{MODE.value}.rules.{bool_rule} 应为 true 或 false，"
                           f"实际为 {params['rules'][bool_rule]!r}")
    for group in ("calibration", "pathway"):
        bad = sorted(set(parameters.get(group) or {}) - set(DEFAULTS[group]))
        if bad:
            raise RunError(f"integration.parameters.{MODE.value}.{group} 中有无法识别的项：{'、'.join(bad)}"
                           f"（可用：{'、'.join(DEFAULTS[group])}）")
        params[group] = {**DEFAULTS[group], **(parameters.get(group) or {})}
    try:
        params["pathway_map"] = pathway_map_settings(parameters.get("pathway_map"))
    except ValueError as exc:
        raise RunError(f"integration.parameters.{MODE.value}.{exc}") from exc
    params["inherited_from_flux_anchored"] = sorted(inherited)
    return params


def _long_index(table: pd.DataFrame | None, tissue_map: dict[str, str] | None = None) -> pd.DataFrame:
    """统一结果表 → 以（组织, 比较, 特征 ID）为索引的统计量。"""
    if table is None or table.empty:
        return pd.DataFrame()
    out = pd.DataFrame({target: table[source] for source, target in _LONG_STATS.items() if source in table.columns})
    tissues = table["tissue"].astype(str)
    out["tissue"] = tissues.map(lambda t: (tissue_map or {}).get(t, t))
    out["contrast_id"] = table["contrast_id"].astype(str)
    out["feature_id"] = table["feature_id"].astype(str)
    return out.drop_duplicates(["tissue", "contrast_id", "feature_id"]).set_index(["tissue", "contrast_id", "feature_id"])


def _frame_for(layer: str, frames: Any, tissue: str, contrast: str) -> pd.DataFrame | None:
    if layer == "transcriptomics":
        return frames.get(contrast, {}).get(tissue)
    frame = frames.get(tissue)
    return frame if frame is not None and f"log2fc_{contrast}" in frame.columns else None


def _state(row: dict[str, Any], rules: dict[str, float]) -> tuple[str, str]:
    """测到的特征的格子状态与原因。"""
    alpha, bound = float(rules["alpha"]), float(rules["equivalence_bound"])
    if row.get("qc_status") in BLOCKED:
        return "qc_blocked", f"QC：{row['qc_status']}"
    padj, pvalue = row.get("padj"), row.get("pvalue")
    if pd.notna(padj) and padj < alpha:
        return "changed", f"padj {padj:.3g} < {alpha:g}"
    if pd.notna(pvalue) and pvalue < alpha:
        return "nominal", f"p {pvalue:.3g} < {alpha:g}，padj 未过"
    fc, se = row.get("log2fc"), row.get("se")
    shrunk = str(row.get("lfc_shrinkage") or "").strip().lower() not in ("", "none", "nan")
    df = _degrees_of_freedom(row)
    if pd.notna(fc) and pd.notna(se) and se > 0 and not shrunk and df >= 1:
        if abs(fc) + float(student_t.ppf(0.95, df)) * se < bound:
            return "equivalent", f"90% 区间（t，自由度 {df:g}）在 ±{bound:g} 内"
    reason = "收缩估计不做等效判定" if shrunk else ("没有 p 值（只有效应量）" if pd.isna(pvalue) else "")
    return "no_evidence", reason


def _degrees_of_freedom(row: dict[str, Any]) -> float:
    """等效检验的自由度：统一结果表的 ``residual_df``（蛋白组的线性模型；不含 limma 的先验自由度，偏保守），
    否则 n1 + n2 − 2；都没有时为 0（不做等效判定）。"""
    df = pd.to_numeric(row.get("residual_df"), errors="coerce")
    if pd.notna(df):
        return float(df)
    n = pd.to_numeric(pd.Series([row.get("n_numerator"), row.get("n_denominator")]), errors="coerce")
    return float(n.sum() - 2) if n.notna().all() else 0.0


def _implied(role: str, fc: float, activity_sign: int, relation: str = "") -> str:
    """成员的变化对通量方向的指向：酶、转运体、产物升为升，底物升为降，调控因子与效应物按作用方向；
    多个复合体共用的亚基与其余角色不指向。"""
    if pd.isna(fc) or fc == 0 or relation == "shared_subunit":
        return ""
    factor = {"direct_enzyme": 1, "carrier": 1, "product": 1, "flux": 1, "substrate": -1}.get(role, 0)
    if role in ("regulator", "effector"):
        factor = activity_sign
    if not factor:
        return ""
    return "up" if np.sign(fc) * factor > 0 else "down"


def _activity_sign(effect: str, mode: str) -> int:
    if mode not in ACTIVITY_MODES:
        return 0
    return {"activates": 1, "inhibits": -1}.get(effect, 0)


class _Context:
    """一个来源的各层帧、统计量索引、可靠性、稳健性标记与通量表；组织名与比较 ID 都已换成共同名称。"""

    def __init__(self, source: DossierSource, imbalance: float) -> None:
        self.source = source
        self.label = source.label
        tables = source.mapped_tables()
        derived = source.mapped(source.derived, Omics.METABOLOMICS.value)
        empty = pd.DataFrame(columns=["contrast_id", "tissue", "feature_id", "evidence_level"])
        id_types = {source.common_tissue(Omics.TRANSCRIPTOMICS.value, t): v
                    for t, v in fa._transcriptomics_id_types(source.run_dir).items()}
        builders = {
            Omics.TRANSCRIPTOMICS.value: lambda t: fa.transcriptomics_frames(t, id_types),
            Omics.PROTEOMICS.value: fa.proteomics_frames,
            Omics.METABOLOMICS.value: lambda t: fa.metabolomics_frames(t, derived),
        }
        self.frames = {layer: build(tables.get(layer, empty)) for layer, build in builders.items()
                       if layer in source.layers and layer in tables}
        self.long = {layer: _long_index(tables.get(layer)) for layer in self.frames}
        self.derived_long = _long_index(derived)
        # 可靠性与稳健性按来源自己的样本表计算，再换成共同名称
        reliability, self.inputs = build_reliability(source.experiment, source.run_dir, source.reliability_tables(), imbalance)
        self.reliability_table = reliability.assign(
            tissue=[source.common_tissue(str(o), str(t)) for o, t in zip(reliability["omics"], reliability["tissue"])]
        ) if not reliability.empty else reliability
        self.reliability = {
            (str(r.omics), source.common_tissue(str(r.omics), str(r.tissue)), source.contrast_map[str(r.contrast_id)]):
                str(r.status)
            for r in reliability.itertuples() if str(r.contrast_id) in source.contrast_map} \
            if not reliability.empty else {}
        flagged, robustness_inputs = fa.single_subject_features(source.run_dir, source.tables)
        self.flagged = {(o, source.common_tissue(o, t), source.contrast_map[c], f) for o, t, c, f in flagged
                        if c in source.contrast_map}
        self.inputs += robustness_inputs
        self.mfa = _mfa_tables(source) if Omics.MFA.value in source.layers else \
            {"pooled": {}, "by_tracer": {}, "inputs": [], "single_tracer": ""}
        self.inputs += self.mfa["inputs"]


def _omics_rows(builder: FluxEvidenceBuilder, ctx: _Context, layer: str, flux_id: str, tissue: str, contrast: str,
                minor_share: float) -> list[dict[str, Any]]:
    frame = _frame_for(layer, ctx.frames[layer], tissue, contrast)
    members = builder.member_matches(layer, flux_id, frame, contrast)
    id_col, name_col, abundance_col = _FRAME_ID[layer]
    qc_status = ctx.reliability.get((layer, tissue, contrast), "ok" if frame is not None else "")
    rows: list[dict[str, Any]] = []
    for m in members.itertuples(index=False):
        row = {"layer": m.layer, "role": m.role, "member_id": m.member_id, "member_relation": m.member_relation,
               "regulator_effect": m.regulator_effect, "regulator_mode": m.regulator_mode,
               "direction_member": bool(m.direction_member), "mapping_method": m.mapping_method,
               "selected": bool(m.selected), "in_flux_anchored": bool(m.in_flux_anchored), "qc_status": qc_status,
               "tracer": "", "flags": []}
        if frame is None:
            row.update(cell_state="not_measured", state_reason="该组织、该比较没有这一层的数据")
        elif m.row is None:
            row.update(cell_state="not_detected", state_reason="知识表中有该成员，数据中没有对应特征")
        else:
            feature_id = str(frame.at[m.row, id_col]) if id_col in frame.columns else ""
            row.update(feature_id=feature_id, feature_name=str(frame.at[m.row, name_col]),
                       _abundance=pd.to_numeric(frame.at[m.row, abundance_col], errors="coerce")
                       if abundance_col in frame.columns else np.nan)
            long = ctx.derived_long if m.layer == "derived" else ctx.long[layer]
            key = (tissue, contrast, feature_id)
            if not long.empty and key in long.index:
                stats = long.loc[key].to_dict()
                row.update({k: v for k, v in stats.items() if k not in ("data_flags", "long_feature_name")})
                row["flags"] += [f for f in str(stats.get("data_flags") or "").split(";") if f and f != "nan"]
            if (layer, tissue, contrast, feature_id) in ctx.flagged:
                row["flags"].append("single_subject_driven")
        rows.append(row)
    _annotate_members(rows, minor_share)
    present = {r["role"] for r in rows}
    for role in _EXPECTED_ROLES[layer]:
        if role not in present and not (role == "direct_enzyme" and "carrier" in present):
            rows.append({"layer": layer, "role": role, "member_id": "", "cell_state": "no_knowledge",
                         "state_reason": "知识表中这个角色没有成员", "in_flux_anchored": False, "selected": False,
                         "direction_member": False, "tracer": "", "flags": [], "qc_status": qc_status})
    return rows


def _condition_abundances(row: dict[str, Any]) -> tuple[float, float]:
    """成员在分子组与分母组中的条件丰度估计（用于条件诱导型同工酶的活跃占比）。"""
    m_num = pd.to_numeric(pd.Series([row.get("mean_numerator")]), errors="coerce").iloc[0]
    m_den = pd.to_numeric(pd.Series([row.get("mean_denominator")]), errors="coerce").iloc[0]
    if pd.notna(m_num) and pd.notna(m_den) and m_num >= 0 and m_den >= 0 and (m_num + m_den) > 0:
        return float(m_num), float(m_den)
    abund = float(row["_abundance"])
    lfc = pd.to_numeric(pd.Series([row.get("log2fc_unshrunk")]), errors="coerce").iloc[0]
    if pd.isna(lfc):
        lfc = pd.to_numeric(pd.Series([row.get("log2fc")]), errors="coerce").iloc[0]
    if pd.isna(lfc):
        return abund, abund
    ratio = float(np.power(2.0, float(np.clip(lfc, -15.0, 15.0))))
    return 2.0 * abund * ratio / (1.0 + ratio), 2.0 * abund / (1.0 + ratio)


def _annotate_members(rows: list[dict[str, Any]], minor_share: float) -> None:
    """同层同工酶成员的丰度占比（含条件活跃占比，避免强诱导型同工酶被整体均值稀释误杀）与重复特征标记。"""
    isozymes = [r for r in rows if r["role"] in ("direct_enzyme", "carrier") and r.get("member_relation") == "isozyme"
                and pd.notna(r.get("_abundance")) and r["_abundance"] > 0]
    total = sum(r["_abundance"] for r in isozymes)
    cond_pairs = [_condition_abundances(r) for r in isozymes]
    total_num = sum(a for a, _ in cond_pairs)
    total_den = sum(b for _, b in cond_pairs)
    for r, (a_num, a_den) in zip(isozymes, cond_pairs):
        pooled_share = r["_abundance"] / total if total else np.nan
        share_num = a_num / total_num if total_num > 0 else np.nan
        share_den = a_den / total_den if total_den > 0 else np.nan
        candidates = [v for v in (pooled_share, share_num, share_den) if pd.notna(v)]
        active_share = max(candidates) if candidates else np.nan
        r["abundance_share"] = pooled_share
        r["active_abundance_share"] = active_share
        if len(isozymes) >= 2 and pd.notna(active_share) and active_share < minor_share:
            r["flags"].append("minor_isozyme")
    counts = pd.Series([(r["role"], r["member_id"]) for r in rows if r.get("feature_id")]).value_counts()
    for r in rows:
        r.pop("_abundance", None)
        if r.get("feature_id") and counts.get((r["role"], r["member_id"]), 0) > 1:
            r["flags"].append("duplicate_feature")


def _mfa_rows(tables: dict[str, Any], flux_id: str, tissue: str, contrast: str, tracers: list[str],
              reliability: dict, flagged: set | None = None) -> list[dict[str, Any]]:
    pooled, by_tracer = tables["pooled"], tables["by_tracer"]
    qc_status = reliability.get((Omics.MFA.value, tissue, contrast), "")
    base = {"layer": "mfa", "role": "flux", "member_id": flux_id, "feature_id": flux_id, "feature_name": flux_id,
            "direction_member": True, "in_flux_anchored": False, "selected": False, "qc_status": qc_status,
            "flags": []}
    rows = []
    measured = by_tracer.get((tissue, contrast), pd.DataFrame())
    for tracer, part in measured.groupby("tracer", sort=True) if not measured.empty else ():
        hit = part[part["feature_id"] == flux_id]
        if hit.empty:
            rows.append({**base, "tracer": tracer, "cell_state": "not_measured", "state_reason": "该示踪剂的结果中没有该反应"})
        else:
            rows.append({**base, **_mfa_stats(hit.iloc[0]), "tracer": tracer})
    infused = set(measured["tracer"]) if not measured.empty else set()
    for tracer in tracers:
        if tracer not in infused:
            rows.append({**base, "tracer": tracer, "cell_state": "not_infused",
                         "state_reason": "声明的示踪剂在该组织、该比较没有数据"})
    frame = pooled.get((tissue, contrast), pd.DataFrame())
    hit = frame[frame["feature_id"] == flux_id] if not frame.empty else frame
    if hit.empty:
        reason = "该组织、该比较没有通量数据" if frame.empty else "通量结果中没有该反应"
        rows.append({**base, "tracer": "pooled", "cell_state": "not_measured", "state_reason": reason})
    else:
        row = {**base, **_mfa_stats(hit.iloc[0]), "tracer": "pooled", "in_flux_anchored": True}
        if (Omics.MFA.value, tissue, contrast, flux_id) in (flagged or set()):
            row["flags"] = row["flags"] + ["single_subject_driven"]    # E5：入选通量逐只去除后反号或减半
        rows.append(row)
    return rows


def _mfa_stats(r: pd.Series) -> dict[str, Any]:
    """通量行的统计量。净通量在两组间变号时，旧流程的 log2FC 为两组绝对值之比，符号不一定与净通量的变化同向
    （如 −28.5 → +8.0 得 −1.84）；这里按组均值之差定号（幅度不变）并标记 ``net_flux_sign_change``。
    示踪剂间异质性取标准化尺度（``smd_heterogeneity_pvalue``，G10），没有时取绝对尺度。"""
    pvalue = r.get("pvalue")
    if pd.isna(pvalue) and pd.notna(r.get("permutation_pvalue")):
        pvalue = r.get("permutation_pvalue")
    flags = [f for f in str(r.get("flags") or "").split(";") if f and f != "nan"]
    log2fc = r.get("log2fc")
    num, den = pd.to_numeric(r.get("mean_numerator"), errors="coerce"), pd.to_numeric(r.get("mean_denominator"),
                                                                                     errors="coerce")
    if pd.notna(num) and pd.notna(den) and num * den < 0 and pd.notna(log2fc) and np.isfinite(log2fc):
        log2fc = float(np.sign(num - den) * abs(log2fc))
        flags.append("net_flux_sign_change")
    size = pd.to_numeric(r.get("flux_group_size"), errors="coerce")
    if pd.notna(size) and size > 1 and isinstance(r.get("flux_group"), str) and r.get("flux_group"):
        flags.append(f"flux_group:{str(r['flux_group']).split('--', 1)[-1]}")   # 化学计量上锁死，与组内反应同一次检验
    heterogeneity = r.get("smd_heterogeneity_pvalue")
    if pd.isna(heterogeneity):
        heterogeneity = r.get("heterogeneity_pvalue")
    return {"log2fc": log2fc, "effect": r.get("effect_size"), "effect_type": r.get("effect_size_type"),
            "se": r.get("log2fc_se"), "pvalue": pvalue, "padj": r.get("padj"), "lfsr": r.get("lfsr"),
            "n_numerator": r.get("n_numerator"), "n_denominator": r.get("n_denominator"),
            "mean_numerator": r.get("mean_numerator"), "mean_denominator": r.get("mean_denominator"),
            "identifiability": r.get("identifiability"), "frac_at_bound": r.get("frac_at_bound"),
            "evidence_level": r.get("evidence_level"), "selected": bool(r.get("selected", False)),
            "heterogeneity_pvalue": heterogeneity, "direction_stability": r.get("direction_stability"),
            "pooled_loo_stability": r.get("pooled_loo_stability"), "tracer_evidence": r.get("tracer_evidence"),
            "smd_se": r.get("smd_se"), "smd_weights": r.get("smd_weights"), "tracer_weights": r.get("tracer_weights"),
            "g_total_var": r.get("g_total_var"), "smd_weighting": r.get("smd_weighting"),
            "padj_global": r.get("padj_global"), "fdr_scope": r.get("fdr_scope"), "flags": flags}


def _mfa_state(row: dict[str, Any], rules: dict[str, float]) -> tuple[str, str]:
    n = pd.to_numeric(pd.Series([row.get("n_numerator"), row.get("n_denominator")]), errors="coerce")
    if n.notna().all() and n.min() < 2:
        return "insufficient_n", f"较小组 {int(n.min())} 只"
    if str(row.get("identifiability") or "") in _UNIDENTIFIABLE:
        return "not_identifiable", str(row["identifiability"])
    return _state(row, rules)


def _mfa_tables(source: DossierSource) -> dict[str, Any]:
    def split(frame: pd.DataFrame | None) -> dict[tuple[str, str], pd.DataFrame]:
        frame = source.mapped(frame, Omics.MFA.value)
        if frame is None or frame.empty:
            return {}
        frame = frame.assign(feature_id=frame["feature_id"].astype(str))
        return {key: part for key, part in frame.groupby(["tissue", "contrast_id"], sort=True)}

    path = source.run_dir / "omics" / Omics.MFA.value / "tables" / "differential_by_tracer.parquet"
    # 分示踪剂的结果以（比较, 组织, 反应, 示踪剂）为键，不是统一差异结果表的契约，直接读取
    by_tracer = pd.read_parquet(path) if path.is_file() else None
    inputs = [path] if by_tracer is not None else []
    pooled = source.tables.get(Omics.MFA.value)
    single = "" if by_tracer is not None else single_tracer(source)
    if single and pooled is not None and not pooled.empty:
        # 单示踪剂实验：MFA 流程不写分示踪剂表，统一差异表就是这一种示踪剂的结果。按示踪剂另列一行（与 pooled 同值），
        # 否则声明的示踪剂被记为“未输注”，F 轴也取不到按示踪剂、由逐只动物计算的方向稳定性
        by_tracer = pooled.assign(tracer=single)
    return {"pooled": split(pooled), "by_tracer": split(by_tracer), "inputs": inputs, "single_tracer": single}


def single_tracer(source: DossierSource) -> str:
    """来源的 MFA 样本只有一种示踪剂（samples.tsv 的 ``protocol`` 只有一个非空取值）时返回它，否则返回空字符串。
    与 MFA 流程的判断一致：``protocol`` 只有一种取值时按单示踪剂分析，不写分示踪剂表。"""
    table = source.experiment.samples
    if "protocol" not in table.columns:
        return ""
    protocol = table.loc[table["omics"] == Omics.MFA.value, "protocol"].fillna("").astype(str).str.strip()
    named = set(protocol) - {""}
    return next(iter(named)) if len(named) == 1 else ""


def _finish(rows: list[dict[str, Any]], rules: dict[str, float], flux_id: str, tissue: str, contrast: str,
            label: str, direction_evidence: str) -> list[dict[str, Any]]:
    for row in rows:
        if "cell_state" not in row:
            row["cell_state"], row["state_reason"] = (_mfa_state if row["layer"] == "mfa" else _state)(row, rules)
        fc, se = row.get("log2fc"), row.get("se")
        if pd.notna(fc) and pd.notna(se):
            row["ci_lower"], row["ci_upper"] = fc - Z95 * se, fc + Z95 * se
        sign = _activity_sign(row.get("regulator_effect", ""), row.get("regulator_mode", ""))
        row["implied_direction"] = _implied(row["role"], fc, sign, str(row.get("member_relation") or "")) \
            if row["cell_state"] in ("changed", "nominal") else ""
        row["flags"] = ";".join(dict.fromkeys(row.get("flags") or []))
        row.update(flux_id=flux_id, tissue=tissue, contrast_id=contrast, source=label,
                   direction_evidence=direction_evidence)
    return rows


def build_dossier(builder: FluxEvidenceBuilder, ctx: _Context, params: dict[str, Any],
                  reactions: list[str], tissues: set[str] | None = None) -> pd.DataFrame:
    """一个来源的成员行。``tissues``：只列这些（共同名称）组织；None 为全部。"""
    mfa = ctx.mfa
    contrasts = sorted({c for frames in ctx.frames.values() for c in _contrasts(frames)}
                       | {key[1] for key in mfa["pooled"]})
    found = sorted({t for frames in ctx.frames.values() for t in _tissues(frames)} | {key[0] for key in mfa["pooled"]})
    # 声明的示踪剂描述本实验；其他来源只列它自己有数据的示踪剂
    tracers = [str(t) for t in params["tracers"]] if ctx.source.is_primary else []
    minor_share = float((params.get("isozyme") or {}).get("minor_share", 0.2))
    records: list[dict[str, Any]] = []
    for contrast in contrasts:
        for tissue in found:
            if tissues is not None and tissue not in tissues:
                continue
            available = [layer for layer, frames in ctx.frames.items()
                         if _frame_for(layer, frames, tissue, contrast) is not None]
            if not available and (tissue, contrast) not in mfa["pooled"]:
                continue
            for flux_id in reactions:
                rows = []
                for layer in ctx.frames:
                    rows += _omics_rows(builder, ctx, layer, flux_id, tissue, contrast, minor_share)
                if Omics.MFA.value in ctx.source.layers:
                    rows += _mfa_rows(mfa, flux_id, tissue, contrast, tracers, ctx.reliability, ctx.flagged)
                semantics = builder.package.semantics(flux_id).direction_evidence
                records += _finish(rows, params["rules"], flux_id, tissue, contrast, ctx.label, semantics)
    out = pd.DataFrame(records)
    for column in COLUMNS:
        if column not in out.columns:
            out[column] = np.nan
    return out[COLUMNS]


def _decoy_packages(builder: FluxEvidenceBuilder, ctx: _Context, reactions: list[str], tissue: str,
                    rng: np.random.Generator) -> dict[str, DecoyPackage]:
    """一个组织的诱饵映射（E8）：每个反应的直接酶（没有时为载体）与底物 / 产物换成同一组织数据中丰度相近的随机特征
    （E3 的抽样），数目与角色不变；每层一个模型包，各层的诱饵互不影响。"""
    tx = next((per[tissue] for per in ctx.frames.get(Omics.TRANSCRIPTOMICS.value, {}).values() if tissue in per), None)
    frames = {Omics.TRANSCRIPTOMICS.value: tx, Omics.PROTEOMICS.value: ctx.frames.get(Omics.PROTEOMICS.value, {}).get(tissue),
              Omics.METABOLOMICS.value: ctx.frames.get(Omics.METABOLOMICS.value, {}).get(tissue)}
    names = {Omics.TRANSCRIPTOMICS.value: "gene", Omics.PROTEOMICS.value: "protein",
             Omics.METABOLOMICS.value: "compound"}
    picked: dict[str, dict[str, Any]] = {layer: {} for layer in frames}
    for layer, frame in frames.items():
        if frame is None or frame.empty or names[layer] not in frame.columns:
            continue
        pool = DecoyPool(frame, layer)
        for flux_id in reactions:
            rows = builder.member_matches(layer, flux_id, frame)
            matched = rows[rows["row"].notna()]
            exclude = set(matched["row"])
            if layer == Omics.METABOLOMICS.value:
                roles = {}
                for key, role in (("substrates", "substrate"), ("products", "product")):
                    targets = list(matched.loc[(matched["layer"] == "metabolomics") & (matched["role"] == role), "row"])
                    chosen = pool.sample(targets, exclude, rng)
                    exclude |= set(chosen)
                    roles[key] = frame.loc[chosen, names[layer]].astype(str).tolist()
                if roles["substrates"] or roles["products"]:
                    picked[layer][flux_id] = roles
                continue
            direct = matched[(matched["role"] == "direct_enzyme") & (matched["member_relation"] != "shared_subunit")]
            if direct.empty:
                direct = matched[matched["role"] == "carrier"]
            chosen = pool.sample(list(direct["row"]), exclude, rng)
            if chosen:
                picked[layer][flux_id] = frame.loc[chosen, names[layer]].astype(str).tolist()
    base = builder.package
    return {Omics.TRANSCRIPTOMICS.value: DecoyPackage(base, genes=picked[Omics.TRANSCRIPTOMICS.value]),
            Omics.PROTEOMICS.value: DecoyPackage(base, proteins=picked[Omics.PROTEOMICS.value]),
            Omics.METABOLOMICS.value: DecoyPackage(base, metabolites=picked[Omics.METABOLOMICS.value])}


def build_decoy_dossier(builder: FluxEvidenceBuilder, contexts: list[_Context], params: dict[str, Any],
                        reactions: list[str], tissues: set[str] | None, seed: int) -> pd.DataFrame:
    """一次诱饵重复的成员行：组学层换成诱饵映射（同一套诱饵名称用于所有来源），通量行与目标相同，不带派生比值。"""
    rng = np.random.default_rng(seed)
    base = builder.package
    primary = contexts[0]
    found = sorted({t for frames in primary.frames.values() for t in _tissues(frames)})
    minor_share = float((params.get("isozyme") or {}).get("minor_share", 0.2))
    records: list[dict[str, Any]] = []
    try:
        for tissue in found:
            if tissues is not None and tissue not in tissues:
                continue
            packages = _decoy_packages(builder, primary, reactions, tissue, rng)
            for ctx in contexts:
                contrasts = sorted({c for frames in ctx.frames.values() for c in _contrasts(frames)}
                                   | {key[1] for key in ctx.mfa["pooled"]})
                tracers = [str(t) for t in params["tracers"]] if ctx.source.is_primary else []
                for contrast in contrasts:
                    if not [layer for layer, frames in ctx.frames.items()
                            if _frame_for(layer, frames, tissue, contrast) is not None]:
                        continue
                    for flux_id in reactions:
                        rows = []
                        for layer in ctx.frames:
                            builder.package = packages[layer]
                            rows += [r for r in _omics_rows(builder, ctx, layer, flux_id, tissue, contrast, minor_share)
                                     if r["layer"] != "derived"]
                        builder.package = base
                        if Omics.MFA.value in ctx.source.layers:
                            rows += _mfa_rows(ctx.mfa, flux_id, tissue, contrast, tracers, ctx.reliability, ctx.flagged)
                        semantics = base.semantics(flux_id).direction_evidence
                        records += _finish(rows, params["rules"], flux_id, tissue, contrast, ctx.label, semantics)
    finally:
        builder.package = base
    out = pd.DataFrame(records)
    for column in COLUMNS:
        if column not in out.columns:
            out[column] = np.nan
    return out[COLUMNS]


def _contrasts(frames: dict) -> set[str]:
    if frames and isinstance(next(iter(frames.values())), dict):
        return set(frames)
    return {c[len("log2fc_"):] for frame in frames.values() for c in frame.columns if c.startswith("log2fc_")}


def _tissues(frames: dict) -> set[str]:
    if frames and isinstance(next(iter(frames.values())), dict):
        return {t for per in frames.values() for t in per}
    return set(frames)


def _normalize_panels_tissues(panels: dict[str, Panel], aliases: dict[str, str]) -> dict[str, Panel]:
    """将面板声明与对照中的组织名按实验级别名表归一，避免产出同名异写的组织目录。"""
    if not aliases:
        return panels
    for panel in panels.values():
        if panel.tissues:
            panel.tissues = list(dict.fromkeys(resolve_tissue(t, aliases) for t in panel.tissues))
        for control in panel.controls:
            if control.get("tissue"):
                control["tissue"] = resolve_tissue(control["tissue"], aliases)
    return panels


def run(ctx: RunContext, source: Path, parameters: dict[str, Any]) -> IntegrationRun:
    params = _parameters(ctx, parameters)
    try:
        tissue_map = validate_tissue_aliases(params["tissue_map"] or {}, f"integration.parameters.{MODE.value}.tissue_map")
        filter_aliases = merge_tissue_aliases(ctx.experiment.config.design.tissue_aliases, tissue_map,
                                              label=f"integration.parameters.{MODE.value}.tissue_map")
    except ValueError as exc:
        raise RunError(str(exc)) from exc
    label = str(params["label"] or ctx.project.id)
    primary = primary_source(ctx.experiment, source, label, tissue_map)
    sources = [primary]
    for spec in params["sources"] or []:
        other = other_source(ctx.project, dict(spec), set(primary.contrast_map.values()))
        if other.label in {s.label for s in sources}:
            raise RunError(f"reaction_dossier.sources 中的来源标签 {other.label} 重复（本实验的标签为 {label}）")
        sources.append(other)
    panels, panels_path = {}, None
    if params["panels"]:
        panels_path = ctx.experiment.root / str(params["panels"])
        try:
            panels = _normalize_panels_tissues(load_panels(panels_path), filter_aliases)
        except PanelError as exc:
            raise RunError(str(exc)) from exc
    if params["pathway_map"].get("tissues"):
        params["pathway_map"] = {
            **params["pathway_map"],
            "tissues": list(dict.fromkeys(resolve_tissue(t, filter_aliases) for t in params["pathway_map"]["tissues"])),
        }

    scoring = {k: params[k] for k in ("pval_cutoff", "log2fc_cutoff", "omics_cutoffs", "fallback_mode", "isozyme")}
    builder = FluxEvidenceBuilder(SimpleNamespace(experiment_id=ctx.experiment.id, estimand={}, scoring=scoring))
    contexts = [_Context(s, float(params["qc_imbalance"])) for s in sources]
    # 评估集：映射表条目与 MFA 结果反应 ID 的并集，去掉登记为别名的旧命名条目（G18 M0a；MFA 结果用到的 ID 保留）
    universe, excluded_aliases = assessment_universe(
        set(FLUX_GENE_MAPPING) | set(REACTION_ENZYME_MAPPING) | set(FLUX_METABOLITE_MAPPING),
        {f for frame in contexts[0].mfa["pooled"].values() for f in frame["feature_id"]})
    if params["reactions"]:
        reactions, excluded_aliases = [str(r) for r in params["reactions"]], []
    else:
        reactions = universe
    tissues = {resolve_tissue(str(t), filter_aliases) for t in params["tissues"]} or None
    dossier = pd.concat([build_dossier(builder, c, params, reactions, tissues) for c in contexts], ignore_index=True)

    base = mode_dir(ctx, MODE)
    out, qc_dir = base / "tables", base / QC_DIR
    out.mkdir(parents=True)
    qc_dir.mkdir(parents=True)
    dossier.to_parquet(out / DOSSIER_TABLE, index=False)
    pd.concat([c.reliability_table.assign(source=c.label) for c in contexts], ignore_index=True) \
        .to_parquet(qc_dir / RELIABILITY_FILE, index=False)
    samples, sides, dataset_dirs = [], [], []
    for c in contexts:
        table, dirs = dossier_samples(c.source, dossier[dossier["source"] == c.label])
        samples.append(table)
        sides.append(contrast_samples(c.source))
        dataset_dirs += dirs
    samples = pd.concat(samples, ignore_index=True)
    sides = pd.concat(sides, ignore_index=True)
    samples.to_parquet(out / SAMPLES_TABLE, index=False)
    sides.to_parquet(out / SIDES_TABLE, index=False)

    # 分轴综合评估（S3）：预期示踪剂取自面板声明（同一反应出现在多个面板时合并）
    expected = {}
    for panel in panels.values():
        for reaction, tracers in panel.expected_tracers.items():
            expected[reaction] = list(dict.fromkeys(expected.get(reaction, []) + tracers))
    rules = {**params["rules"], "minor_share": float((params.get("isozyme") or {}).get("minor_share", 0.2))}
    axes, reactions_assessed = assess(dossier, label, rules, expected, samples, sides)
    axes.to_parquet(out / AXIS_TABLE, index=False)
    reactions_assessed.to_parquet(out / ASSESSMENT_TABLE, index=False)

    # 判定规则的校准（S4，E8）
    calibration = calibrate(builder, contexts, params, reactions, tissues, dossier, (axes, reactions_assessed),
                            samples, sides, expected, rules, panels)
    (base / "evaluation").mkdir(parents=True)
    calibration.to_parquet(base / "evaluation" / CALIBRATION_TABLE, index=False)

    # 通路层的综合读数（S5，G16）
    pathway = None
    if panels:
        pathway, module_dirs = pathway_readouts(panels, dossier, samples, sides, contexts, params)
        dataset_dirs += module_dirs
        pathway["splits"].to_parquet(out / NODE_SPLIT_TABLE, index=False)
        pathway["split_samples"].to_parquet(out / NODE_SPLIT_SAMPLES, index=False)
        pathway["tests"].to_parquet(out / PATHWAY_TESTS_TABLE, index=False)
    panel_notes = render_panels(base, panels, dossier, samples, sides, contexts, params,
                                (axes, reactions_assessed), calibration, pathway) if panels else []

    # 分组织代谢通路总图（图 E，G18）：不依赖面板
    map_layout = params["pathway_map"]["layout"]
    try:
        map_notes = render_maps(base, base, params["pathway_map"], label,
                                {"整合": ctx.run_dir.name, **{f"来源 {c.label}": c.source.run_id for c in contexts}},
                                members=dossier, axes=axes, assessed=reactions_assessed, calibration=calibration,
                                root=ctx.experiment.root, rules=rules)
    except LayoutError as exc:
        raise RunError(str(exc)) from exc
    # 课题自己的布局文件登记为输入（包内默认布局随代码版本走）
    layout_inputs = [ctx.experiment.root / str(map_layout)] if map_layout else []

    inputs = [p for s in sources for p in s.inputs] + [p for c in contexts for p in c.inputs] \
        + [d / name for d in dict.fromkeys(dataset_dirs) for name in ("matrix.parquet", META_FILE)
           if (d / name).is_file()]
    record = IntegrationRun(MODE, inputs=inputs + ([panels_path] if panels_path else []) + layout_inputs,
                            parameters={**params, "label": label, "excluded_aliases": excluded_aliases,
                                        "sources_resolved": {s.label: {"run": s.run_id, "layers": list(s.layers),
                                                                       "contrast_map": s.contrast_map}
                                                             for s in sources},
                                        "omics": sorted(primary.tables)})
    states = dossier["cell_state"].value_counts()
    record.notes.append(f"反应档案：{len(sources)} 个来源、{dossier['flux_id'].nunique()} 个反应、{len(dossier)} 行；"
                        "格子状态 " + "、".join(f"{state} {n}" for state, n in states.items()))
    if excluded_aliases:
        record.notes.append(f"评估集不含 {len(excluded_aliases)} 个登记为网络反应别名的旧命名条目"
                            f"（reaction_aliases.yaml）：{'、'.join(excluded_aliases)}")
    for ctx in contexts:
        if ctx.mfa.get("single_tracer"):
            record.notes.append(f"来源 {ctx.label} 的 MFA 只有一种示踪剂（{ctx.mfa['single_tracer']}），没有分示踪剂表："
                                "统一差异表即该示踪剂的结果，按示踪剂另列一行")
    if rules["rna_single_source_changed"]:
        single = int(axes["flags"].fillna("").str.contains("rna_single_source").sum())
        record.notes.append(f"规则 rna_single_source_changed 已打开：{single} 个格子的 C 轴由单一来源通过 FDR 的 RNA "
                            "判为 changed（标 rna_single_source，不算复现）")
    declared = set(params["tracers"])
    seen = set(dossier.loc[dossier["layer"].eq("mfa") & dossier["cell_state"].ne("not_infused"), "tracer"]) - {"pooled"}
    if declared - seen:
        record.notes.append(f"声明的示踪剂 {'、'.join(sorted(declared - seen))} 在结果包中没有数据，记为 not_infused")
    if not reactions_assessed.empty:
        claims = reactions_assessed["flux_claim"].value_counts()
        record.notes.append("分轴综合评估：结论等级 " + "、".join(f"{k} {claims.get(k, 0)}"
                                                       for k in ("L3", "L2", "L1", "L0", "LX"))
                            + f"（{len(reactions_assessed)} 个反应 × 组织 × 比较）")
    record.notes += calibration_notes(calibration)
    record.notes += pathway_notes(pathway)
    record.notes += panel_notes
    record.notes += map_notes
    return record


def _module_data(contexts: list[_Context]) -> tuple[Any, list[Path]]:
    """模块检验的数据：按（来源, 层, 组织）读取分析尺度的数据集，基因按与成员表相同的规则匹配。"""
    cache: dict[tuple[str, str, str], pw.ModuleData | None] = {}
    read: list[Path] = []
    by_label = {c.label: c for c in contexts}

    def provide(source: str, layer: str, tissue: str) -> pw.ModuleData | None:
        key = (source, layer, tissue)
        if key in cache:
            return cache[key]
        ctx, out = by_label[source], None
        frames = ctx.frames.get(layer, {})
        frame = next((per[tissue] for per in frames.values() if tissue in per), None) \
            if layer == Omics.TRANSCRIPTOMICS.value else frames.get(tissue)
        directory = _dataset_dir(ctx.source, layer, tissue)
        if frame is not None and directory.is_dir():
            dataset = read_dataset(directory)
            values = analysis_values(dataset)
            matcher = match_gene_features if layer == Omics.TRANSCRIPTOMICS.value else match_protein_features
            id_column = _FRAME_ID[layer][0]

            def match(genes: list[str], frame=frame, matcher=matcher, id_column=id_column) -> dict[str, list[str]]:
                found = matcher(frame, genes)
                by_query: dict[str, list[str]] = {}
                for index in found.indices:
                    by_query.setdefault(found.query_by_index.get(index, ""), []).append(str(frame.at[index, id_column]))
                return {g: by_query.get(_canonical_gene_symbol(g) or g, []) for g in genes}

            # RNA：估计先验时去掉低表达基因（中位 log2(CPM + 1) < 1）
            background = values.index[values.median(axis=1) >= 1.0] if dataset.meta.value_type.value == "counts" \
                else None
            out = pw.ModuleData(values, match, background)
            read.append(directory)
        cache[key] = out
        return out

    return provide, read


def pathway_readouts(panels: dict[str, Panel], dossier: pd.DataFrame, samples: pd.DataFrame, sides: pd.DataFrame,
                     contexts: list[_Context], params: dict[str, Any]) -> tuple[dict[str, pd.DataFrame], list[Path]]:
    """S5：每个面板的节点分流比、模块检验与交叉点分析。"""
    labels = [c.label for c in contexts]
    provide, read = _module_data(contexts)
    tissues = {s: sorted(set(dossier.loc[dossier["source"] == s, "tissue"])) for s in labels}
    contrasts = {s: sorted(set(dossier.loc[dossier["source"] == s, "contrast_id"])) for s in labels}
    settings = params["pathway"]
    splits, split_samples, tests = [], [], []
    for panel in panels.values():
        table, per_animal = pw.node_splits(panel, dossier, samples, sides, labels)
        splits.append(table)
        split_samples.append(per_animal)
        tests.append(pw.signed_module_tests(panel, provide, sides, labels, tissues, contrasts,
                                            int(settings["rotations"]), int(settings["seed"])))
        tests.append(pw.crossover_tests(panel, dossier, labels))
    combine = (lambda parts, columns: pd.concat([p for p in parts if not p.empty], ignore_index=True)
               if any(not p.empty for p in parts) else pd.DataFrame(columns=columns))
    return {"splits": combine(splits, pw.NODE_SPLIT_COLUMNS),
            "split_samples": combine(split_samples, pw.SPLIT_SAMPLE_COLUMNS),
            "tests": combine(tests, pw.TEST_COLUMNS)}, read


def pathway_notes(pathway: dict[str, pd.DataFrame] | None) -> list[str]:
    if not pathway:
        return []
    notes = []
    splits, tests = pathway["splits"], pathway["tests"]
    if not splits.empty:
        computed = splits[splits["state"] == "computed"]
        notes.append(f"节点分流比：{len(splits)} 个（节点 × 出口对 × 层 × 来源 × 组织 × 比较），计算了 {len(computed)} 个，"
                     f"padj < 0.05 的 {int((computed['padj'] < 0.05).sum())} 个；其余状态 "
                     + "、".join(f"{k} {v}" for k, v in splits.loc[splits["state"] != "computed", "state"]
                                .value_counts().items()))
    for test, label in (("signed_module", "模块检验"), ("crossover", "交叉点分析")):
        part = tests[tests["test"] == test] if not tests.empty else tests
        if part.empty:
            continue
        tested = part[part["pvalue"].notna()]
        notes.append(f"{label}：{len(part)} 项，检验了 {len(tested)} 项，padj < 0.05 的 {int((tested['padj'] < 0.05).sum())} 项")
    return notes


def calibrate(builder: FluxEvidenceBuilder, contexts: list[_Context], params: dict[str, Any], reactions: list[str],
              tissues: set[str] | None, dossier: pd.DataFrame, assessed: tuple[pd.DataFrame, pd.DataFrame],
              samples: pd.DataFrame, sides: pd.DataFrame, expected: dict[str, list[str]], rules: dict[str, Any],
              panels: dict[str, Panel]) -> pd.DataFrame:
    """E8：诱饵（参数 calibration.decoys 次）、面板声明的对照集、实际数据上的不变性检查；组别置换记为未运行。"""
    settings = params["calibration"]
    label = contexts[0].label
    parts = []
    if int(settings["decoys"]) > 0:
        runs = []
        for replicate in range(int(settings["decoys"])):
            members = build_decoy_dossier(builder, contexts, params, reactions, tissues,
                                          int(settings["seed"]) + replicate)
            runs.append(assess(members, label, rules, expected, samples, sides))
        parts.append(cal.decoy_summary(assessed, runs))
    controls = [c for p in panels.values() for c in p.expand_controls(sorted(set(dossier["tissue"])),
                                                                       sorted(set(dossier["contrast_id"])))]
    if controls:
        parts.append(cal.control_summary(assessed[1], assessed[0], controls))
    if settings["invariance"]:
        subset = sorted({r for p in panels.values() for r in p.reactions}) or reactions[:INVARIANCE_FALLBACK]
        members = dossier[dossier["flux_id"].isin(subset)]
        parts.append(cal.invariance_checks(members, lambda m: assess(m, label, rules, expected, samples, sides)[0]))
    parts.append(cal.not_run("组别置换（E1）需要按置换后的标签重跑各组学的差异分析，尚未实现"))
    return pd.concat(parts, ignore_index=True)


def calibration_notes(calibration: pd.DataFrame) -> list[str]:
    notes = []
    overall = calibration[(calibration["method"] == "decoy") & (calibration["tissue"] == "")].set_index("metric")
    if not overall.empty:
        notes.append(f"诱饵校准（{int(overall['n_replicates'].iloc[0])} 次）：" + "；".join(
            f"{metric} 目标 {int(overall.at[metric, 'target'])}、诱饵期望 {overall.at[metric, 'decoy_expected']:.1f}"
            for metric in cal.METRICS if metric in overall.index) + "；明细见 evaluation/")
    invariance = calibration[calibration["method"] == "invariance"]
    if not invariance.empty:
        failed = invariance[invariance["passed"] != True]  # noqa: E712
        notes.append("不变性检查：" + ("全部通过" if failed.empty else
                                    "未通过 " + "、".join(f"{r.metric}（{int(r.target)} 处）" for r in failed.itertuples()))
                     + f"（{int(invariance['n_cells'].max())} 个格子）")
    controls = calibration[calibration["metric"] == "concordance"]
    for row in controls.itertuples():
        notes.append(f"对照集（{row.method.replace('control_', '')}）：{row.detail}")
    return notes


def render_panels(base: Path, panels: dict[str, Panel], dossier: pd.DataFrame, samples: pd.DataFrame,
                  sides: pd.DataFrame, contexts: list[_Context], params: dict[str, Any],
                  assessment: tuple[pd.DataFrame, pd.DataFrame] | None = None,
                  calibration: pd.DataFrame | None = None,
                  pathway: dict[str, pd.DataFrame] | None = None) -> list[str]:
    """面板的行表、矩阵视图与导出表、图 A / B / C / D。"""
    from ..plotting import dossier_plots

    figures = {**DEFAULTS["figures"], **(params.get("figures") or {})}
    formats = [str(f) for f in figures["formats"]]
    labels = [c.label for c in contexts]
    runs = {c.label: c.source.run_id for c in contexts}
    elements_all, notes = [], []
    primary = dossier[dossier["source"] == contexts[0].label]
    for panel in panels.values():
        elements = panel_elements(panel)
        elements_all.append(elements)
        tissues = panel.tissues or sorted(set(dossier["tissue"]))
        contrasts = panel.contrasts or sorted(set(dossier["contrast_id"]))
        for chunk_index, chunk in enumerate(_chunks(panel, elements, int(figures["max_reactions"]))):
            suffix = f"_{chunk_index + 1}" if chunk_index or len(panel.reactions) > int(figures["max_reactions"]) else ""
            for tissue in tissues:
                folder = base / "figures" / panel.id / tissue
                for contrast in contrasts:
                    view = matrix_view(panel, chunk, dossier, tissue, contrast, labels, params["tracers"], assessment)
                    if view.empty:
                        continue
                    stem = folder / f"matrix_{contrast}{suffix}"
                    export_matrix(view, member_table(dossier, panel, tissue, contrast), stem)
                    dossier_plots.plot_matrix(view, panel, tissue, contrast, stem, formats, runs,
                                              cal.caption(calibration, tissue, contrast))
        for tissue in panel.tissues or sorted(set(dossier["tissue"])):
            folder = base / "figures" / panel.id / tissue
            for contrast in panel.contrasts or sorted(set(dossier["contrast_id"])):
                dossier_plots.plot_tracer_coverage(primary, panel, tissue, contrast,
                                                   folder / f"tracer_coverage_{contrast}", formats, params["tracers"])
                for reaction in panel.reactions:
                    rows = dossier[(dossier["flux_id"] == reaction) & (dossier["tissue"] == tissue)
                                   & (dossier["contrast_id"] == contrast)]
                    if rows.empty:
                        continue
                    dossier_plots.plot_dossier_page(rows, samples, sides, reaction, tissue, contrast,
                                                    folder / f"dossier_{reaction}_{contrast}", formats, runs)
                if pathway is not None:
                    _render_pathway(pathway, panel, tissue, contrast, folder, sides, formats, runs, dossier_plots)
        notes.append(f"面板 {panel.id}：{len(panel.reactions)} 个反应，图与导出表在 integration/{MODE.value}/figures/{panel.id}/")
    pd.concat(elements_all, ignore_index=True).to_parquet(base / "tables" / ELEMENTS_TABLE, index=False)
    return notes


def _render_pathway(pathway: dict[str, pd.DataFrame], panel: Panel, tissue: str, contrast: str, folder: Path,
                    sides: pd.DataFrame, formats: list[str], runs: dict[str, str], dossier_plots: Any) -> None:
    """图 D（节点分流）与通路层读数的导出表（该面板、组织、比较）。"""
    def here(table: pd.DataFrame) -> pd.DataFrame:
        return table[(table["panel_id"] == panel.id) & (table["tissue"] == tissue)
                     & (table["contrast_id"] == contrast)] if not table.empty else table

    splits, tests = here(pathway["splits"]), here(pathway["tests"])
    folder.mkdir(parents=True, exist_ok=True)
    if not splits.empty:
        per_animal = pathway["split_samples"]
        per_animal = per_animal[(per_animal["panel_id"] == panel.id) & (per_animal["tissue"] == tissue)]
        dossier_plots.plot_node_split(splits, per_animal, sides, panel, tissue, contrast,
                                      folder / f"node_split_{contrast}", formats, runs)
        write_excel_tsv(splits, folder / f"node_split_{contrast}.tsv")
    if not tests.empty:
        write_excel_tsv(tests, folder / f"pathway_tests_{contrast}.tsv")


def _chunks(panel: Panel, elements: pd.DataFrame, limit: int) -> list[pd.DataFrame]:
    """一张矩阵不超过 ``limit`` 个反应：超过时按分段拆成几张（分段本身超过时再按反应拆）。"""
    groups, current, count = [], [], 0
    for segment, reactions in panel.segments:
        if current and count + len(reactions) > limit:
            groups.append(current)
            current, count = [], 0
        current.append(segment)
        count += len(reactions)
    groups.append(current)
    return [elements[elements["segment"].isin(g)].reset_index(drop=True) for g in groups]
