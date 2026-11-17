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
（只列这些；默认全部）、``label``（来源标签，默认项目 ID）、``rules``（``alpha``、``equivalence_bound``）；
``pval_cutoff``、``log2fc_cutoff``、``omics_cutoffs``、``fallback_mode``、``isozyme``、``tissue_map``、``qc_imbalance``
没有写时沿用 ``flux_anchored`` 的参数，使入选规则一致。

输出（``integration/reaction_dossier/``）：``tables/dossier_members.parquet``（列见 :data:`COLUMNS`）、
``qc/contrast_reliability.parquet``。
"""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import t as student_t

from omicsweave.contract.tables import read_differential
from omicsweave.core.model import IntegrationMode, Omics
from omicsweave.core.reliability import BLOCKED, DEFAULT_IMBALANCE, QC_DIR, RELIABILITY_FILE, build_reliability
from omicsweave.core.run import RunContext, RunError

from ..integrate import IntegrationRun, mode_dir
from ..mappings import FLUX_GENE_MAPPING, FLUX_METABOLITE_MAPPING, REACTION_ENZYME_MAPPING
from ..mappings.reaction_members import ACTIVITY_MODES
from ..methods.flux_evidence import FluxEvidenceBuilder
from . import flux_anchored as fa

MODE = IntegrationMode.REACTION_DOSSIER
# 与 flux_anchored 共用（没有写时沿用其参数）
SHARED = ("pval_cutoff", "log2fc_cutoff", "omics_cutoffs", "fallback_mode", "isozyme", "tissue_map", "qc_imbalance")
DEFAULTS: dict[str, Any] = {**{key: fa.DEFAULTS[key] for key in SHARED}, "tracers": [], "reactions": [],
                            "tissues": [], "label": "", "rules": {"alpha": 0.05, "equivalence_bound": 0.58}}
DOSSIER_TABLE = "dossier_members.parquet"
Z95 = 1.959963984540054              # ci_lower / ci_upper 的 95% 区间

COLUMNS = [
    "flux_id", "tissue", "contrast_id", "source", "layer", "role", "member_id", "feature_id", "feature_name", "tracer",
    "member_relation", "regulator_effect", "regulator_mode", "direction_member", "mapping_method",
    "log2fc", "log2fc_unshrunk", "effect", "effect_type", "se", "ci_lower", "ci_upper", "pvalue", "padj", "lfsr",
    "n_numerator", "n_denominator", "mean_numerator", "mean_denominator", "abundance_share",
    "cell_state", "state_reason", "selected", "identifiability", "frac_at_bound", "evidence_level", "qc_status",
    "flags", "implied_direction", "direction_evidence", "in_flux_anchored",
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
    params["rules"] = {**DEFAULTS["rules"], **(parameters.get("rules") or {})}
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
    """一次运行中各层的帧、统计量索引、可靠性与稳健性标记。"""

    def __init__(self, tables: dict[str, pd.DataFrame], source: Path, derived: pd.DataFrame | None,
                 tissue_map: dict[str, str], reliability: pd.DataFrame, flagged: set) -> None:
        empty = pd.DataFrame(columns=["contrast_id", "tissue", "feature_id", "evidence_level"])
        self.frames = {
            Omics.TRANSCRIPTOMICS.value: fa.transcriptomics_frames(tables.get(Omics.TRANSCRIPTOMICS.value, empty),
                                                                  fa._transcriptomics_id_types(source)),
            Omics.PROTEOMICS.value: fa.proteomics_frames(tables.get(Omics.PROTEOMICS.value, empty)),
            Omics.METABOLOMICS.value: fa.metabolomics_frames(tables.get(Omics.METABOLOMICS.value, empty), derived),
        }
        self.long = {omics: _long_index(tables.get(omics)) for omics in self.frames}
        self.derived_long = _long_index(derived)
        self.reliability = {(str(r.omics), str(r.tissue), str(r.contrast_id)): str(r.status)
                            for r in reliability.itertuples()} if not reliability.empty else {}
        self.flagged = flagged
        self.tissue_map = tissue_map


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


def _annotate_members(rows: list[dict[str, Any]], minor_share: float) -> None:
    """同层同工酶成员的丰度占比（次要同工酶标记）与重复特征标记。"""
    isozymes = [r for r in rows if r["role"] in ("direct_enzyme", "carrier") and r.get("member_relation") == "isozyme"
                and pd.notna(r.get("_abundance")) and r["_abundance"] > 0]
    total = sum(r["_abundance"] for r in isozymes)
    for r in isozymes:
        r["abundance_share"] = r["_abundance"] / total if total else np.nan
        if len(isozymes) >= 2 and r["abundance_share"] < minor_share:
            r["flags"].append("minor_isozyme")
    counts = pd.Series([(r["role"], r["member_id"]) for r in rows if r.get("feature_id")]).value_counts()
    for r in rows:
        r.pop("_abundance", None)
        if r.get("feature_id") and counts.get((r["role"], r["member_id"]), 0) > 1:
            r["flags"].append("duplicate_feature")


def _mfa_rows(tables: dict[str, Any], flux_id: str, tissue: str, contrast: str, tracers: list[str],
              reliability: dict) -> list[dict[str, Any]]:
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
        rows.append({**base, **_mfa_stats(hit.iloc[0]), "tracer": "pooled", "in_flux_anchored": True})
    return rows


def _mfa_stats(r: pd.Series) -> dict[str, Any]:
    pvalue = r.get("pvalue")
    if pd.isna(pvalue) and pd.notna(r.get("permutation_pvalue")):
        pvalue = r.get("permutation_pvalue")
    return {"log2fc": r.get("log2fc"), "effect": r.get("effect_size"), "effect_type": r.get("effect_size_type"),
            "se": r.get("log2fc_se"), "pvalue": pvalue, "padj": r.get("padj"), "lfsr": r.get("lfsr"),
            "n_numerator": r.get("n_numerator"), "n_denominator": r.get("n_denominator"),
            "mean_numerator": r.get("mean_numerator"), "mean_denominator": r.get("mean_denominator"),
            "identifiability": r.get("identifiability"), "frac_at_bound": r.get("frac_at_bound"),
            "evidence_level": r.get("evidence_level"), "selected": bool(r.get("selected", False)),
            "flags": [f for f in str(r.get("flags") or "").split(";") if f and f != "nan"]}


def _mfa_state(row: dict[str, Any], rules: dict[str, float]) -> tuple[str, str]:
    n = pd.to_numeric(pd.Series([row.get("n_numerator"), row.get("n_denominator")]), errors="coerce")
    if n.notna().all() and n.min() < 2:
        return "insufficient_n", f"较小组 {int(n.min())} 只"
    if str(row.get("identifiability") or "") in _UNIDENTIFIABLE:
        return "not_identifiable", str(row["identifiability"])
    return _state(row, rules)


def _mfa_tables(tables: dict[str, pd.DataFrame], source: Path, tissue_map: dict[str, str]) -> dict[str, Any]:
    def split(frame: pd.DataFrame | None) -> dict[tuple[str, str], pd.DataFrame]:
        if frame is None or frame.empty:
            return {}
        frame = frame.assign(tissue=frame["tissue"].astype(str).map(lambda t: tissue_map.get(t, t)),
                             contrast_id=frame["contrast_id"].astype(str), feature_id=frame["feature_id"].astype(str))
        return {key: part for key, part in frame.groupby(["tissue", "contrast_id"], sort=True)}

    path = source / "omics" / Omics.MFA.value / "tables" / "differential_by_tracer.parquet"
    # 分示踪剂的结果以（比较, 组织, 反应, 示踪剂）为键，不是统一差异结果表的契约，直接读取
    by_tracer = pd.read_parquet(path) if path.is_file() else None
    return {"pooled": split(tables.get(Omics.MFA.value)), "by_tracer": split(by_tracer),
            "inputs": [path] if by_tracer is not None else []}


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


def build_dossier(builder: FluxEvidenceBuilder, ctx: _Context, mfa: dict[str, Any], params: dict[str, Any],
                  label: str) -> pd.DataFrame:
    tissue_map = ctx.tissue_map
    contrasts = sorted({c for frames in ctx.frames.values() for c in _contrasts(frames)}
                       | {key[1] for key in mfa["pooled"]})
    tissues = sorted({t for frames in ctx.frames.values() for t in _tissues(frames)} | {key[0] for key in mfa["pooled"]})
    if params["tissues"]:
        tissues = [t for t in tissues if t in {tissue_map.get(str(x), str(x)) for x in params["tissues"]}]
    reactions = [str(r) for r in params["reactions"]] or sorted(
        set(FLUX_GENE_MAPPING) | set(REACTION_ENZYME_MAPPING) | set(FLUX_METABOLITE_MAPPING)
        | {f for frame in mfa["pooled"].values() for f in frame["feature_id"]})
    tracers = [str(t) for t in params["tracers"]]
    minor_share = float((params.get("isozyme") or {}).get("minor_share", 0.2))
    records: list[dict[str, Any]] = []
    for contrast in contrasts:
        for tissue in tissues:
            available = [layer for layer, frames in ctx.frames.items()
                         if _frame_for(layer, frames, tissue, contrast) is not None]
            if not available and (tissue, contrast) not in mfa["pooled"]:
                continue
            for flux_id in reactions:
                rows = []
                for layer in ctx.frames:
                    rows += _omics_rows(builder, ctx, layer, flux_id, tissue, contrast, minor_share)
                rows += _mfa_rows(mfa, flux_id, tissue, contrast, tracers, ctx.reliability)
                semantics = builder.package.semantics(flux_id).direction_evidence
                records += _finish(rows, params["rules"], flux_id, tissue, contrast, label, semantics)
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


def run(ctx: RunContext, source: Path, parameters: dict[str, Any]) -> IntegrationRun:
    params = _parameters(ctx, parameters)
    paths = {p.parent.parent.name: p for p in sorted(source.glob("omics/*/tables/differential.parquet"))}
    if not paths:
        raise RunError(f"{MODE.value} 需要至少 1 个组学的差异结果；结果包 {source.name} 中没有。"
                       f"请先运行 omicsweave run（含 differential 步骤）")
    tables = {omics: read_differential(path) for omics, path in paths.items()}
    if Omics.METABOLOMICS.value in tables:
        tables[Omics.METABOLOMICS.value], _ = fa.refresh_metabolite_ids(tables[Omics.METABOLOMICS.value])
    derived_path = source / "omics" / Omics.METABOLOMICS.value / "tables" / fa.DERIVED_TABLE
    derived = read_differential(derived_path) if derived_path.is_file() else None
    tissue_map = {str(k): str(v) for k, v in (params["tissue_map"] or {}).items()}

    reliability, quality_inputs = build_reliability(ctx.experiment, source, tables, float(params["qc_imbalance"]))
    qc_dir = mode_dir(ctx, MODE) / QC_DIR
    qc_dir.mkdir(parents=True)
    reliability.to_parquet(qc_dir / RELIABILITY_FILE, index=False)
    mapped = reliability.copy()
    if not mapped.empty:
        is_mfa = mapped["omics"] == Omics.MFA.value
        mapped.loc[is_mfa, "tissue"] = mapped.loc[is_mfa, "tissue"].map(lambda t: tissue_map.get(t, t))
    flagged, robustness_inputs = fa.single_subject_features(source, tables)

    scoring = {k: params[k] for k in ("pval_cutoff", "log2fc_cutoff", "omics_cutoffs", "fallback_mode", "isozyme")}
    builder = FluxEvidenceBuilder(SimpleNamespace(experiment_id=ctx.experiment.id, estimand={}, scoring=scoring))
    context = _Context(tables, source, derived, tissue_map, mapped, flagged)
    mfa = _mfa_tables(tables, source, tissue_map)
    label = str(params["label"] or ctx.project.id)
    dossier = build_dossier(builder, context, mfa, params, label)

    out = mode_dir(ctx, MODE) / "tables"
    out.mkdir(parents=True)
    dossier.to_parquet(out / DOSSIER_TABLE, index=False)

    record = IntegrationRun(MODE, inputs=list(paths.values()) + ([derived_path] if derived is not None else [])
                            + mfa["inputs"] + quality_inputs + robustness_inputs,
                            parameters={**params, "omics": sorted(paths), "label": label})
    states = dossier["cell_state"].value_counts()
    record.notes.append(f"反应档案：{dossier['flux_id'].nunique()} 个反应、{len(dossier)} 行；格子状态 "
                        + "、".join(f"{state} {n}" for state, n in states.items()))
    declared = set(params["tracers"])
    seen = set(dossier.loc[dossier["layer"].eq("mfa"), "tracer"]) - {"pooled"}
    if declared - seen:
        record.notes.append(f"声明的示踪剂 {'、'.join(sorted(declared - seen))} 在结果包中没有数据，记为 not_infused")
    return record
