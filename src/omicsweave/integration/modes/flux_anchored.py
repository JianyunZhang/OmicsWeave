#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
通量锚定（``flux_anchored``）：以 MFA 差异通量为锚点，对每个（反应, 组织, 比较），按反应 → 基因 / 蛋白 /
代谢物的策展映射，从其他组学的差异结果中提取证据，汇总为多维证据画像（支持、冲突、覆盖、稳定性、设计等级、
证据等级），见 :class:`omicsweave.integration.methods.flux_evidence.FluxEvidenceBuilder`。

读取源结果包中各组学的 ``differential.parquet``，前提：MFA 与至少一个其他组学有差异结果。统一差异结果表按
证据构建器的输入整理，比较按 ``contrast_id`` 对应：

- MFA：每个比较一张表，``feature_id`` 为反应 ID；p 值与 FDR 取主检验路径，Cohen's d 取附加列 ``cohen_d``
  （没有该列时用 ``effect_size``），标准误取附加列 ``log2fc_se``；任一组生物学重复少于 2 个或只给效应量时，
  MFA 显著性判定允许改用效应量（与旧流程的 ``stat_note`` 相同）；
- 转录组：每个（比较, 组织）一张表，基因按 ``feature_id``（类型见数据集的 ``feature_id_type``）匹配；
- 蛋白组、代谢组：每个组织一张宽表，统计量列带比较后缀（``log2fc_<比较>`` 等）。蛋白按 ``feature_id`` 与
  附加列 ``uniprot_id`` 匹配；代谢物按名称（``feature_name``，没有时用 ``feature_id``）与 HMDB / KEGG /
  ChEBI 附加列匹配；
- 各层的有效样本量取 ``n_numerator + n_denominator``。

参数（``integration.parameters.flux_anchored``）：``pval_cutoff`` / ``log2fc_cutoff``（特征显著的阈值）、
``omics_cutoffs``（按组学覆盖阈值）、``fallback_mode``（样本量不足时的效应量判定）、``evidence_profile``
（证据画像的阈值与权重；``mfa_inference_unit`` 默认取 MFA 结果的推断单位），含义同旧流程配置中的
``scoring``；``tissue_map`` 把 MFA 的组织名映射到其他组学的组织名。

证据等级取决于锚点状态（``anchor_state``）：MFA 的证据等级为 insufficient、较小组少于 2 只、不可辨识
（``flags`` / ``identifiability``）、示踪剂方向不一致或 QC 不可评估时，锚点不成立，只给组学之间一致性的等级
（Omics-*，见 :mod:`omicsweave.integration.methods.evidence_profile`）。样本 QC 标记（样本表的 ``qc_flag`` 与各组学的
``qc/sample_quality.parquet``）按当前样本表计算为（组学, 组织, 比较）的可靠性（:mod:`omicsweave.core.reliability`，
参数 ``qc_imbalance``），降权或使该层不可评估。

输出（``integration/flux_anchored/``）：

- ``tables/flux_evidence_profile.parquet``：每个（反应, 组织, 比较）一行的证据画像，另列出各层与通量方向一致的特征；
- ``tables/flux_evidence_table.parquet``：每个（反应, 组织, 比较, 层）一行的证据长表；
- ``tables/flux_pathway_profile.parquet``：按（比较, 组织, 通路）汇总的证据画像；
- ``qc/contrast_reliability.parquet``：本次整合使用的可靠性表。
"""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd

from omicsweave.contract.dataset import META_FILE
from omicsweave.contract.tables import read_differential
from omicsweave.core.model import IntegrationMode, Omics
from omicsweave.core.reliability import DEFAULT_IMBALANCE, QC_DIR, RELIABILITY_FILE, build_reliability, summary_note
from omicsweave.core.run import RunContext, RunError

from ..integrate import IntegrationRun, mode_dir
from ..methods.flux_evidence import FluxEvidenceBuilder

MODE = IntegrationMode.FLUX_ANCHORED
DEFAULTS: dict[str, Any] = {
    "pval_cutoff": 0.05, "log2fc_cutoff": 0.58, "omics_cutoffs": {}, "fallback_mode": {},
    "evidence_profile": {}, "tissue_map": {}, "qc_imbalance": DEFAULT_IMBALANCE,
    "isozyme": {"minor_share": 0.2}, "decoys": {"replicates": 0, "seed": 20261115},
}
# 不传给证据构建器的参数
_MODE_ONLY = ("tissue_map", "qc_imbalance", "decoys")
# 转录组 feature_id_type → 证据构建器按 ID 匹配基因的列；其余类型按基因符号匹配
_GENE_ID_COLUMNS = {"ensembl_gene_id": "ensembl_gene_id", "ensembl_id": "ensembl_gene_id",
                    "ncbi_gene_id": "ncbi_gene_id", "entrez_id": "ncbi_gene_id"}
# 统一结果表中的附加列 → 证据构建器的列
_PROTEIN_IDS = {"uniprot_id": "uniprot_id"}
_METABOLITE_IDS = {"hmdb": "hmdb_id", "kegg_compound": "kegg_compound", "chebi": "chebi_id"}
_MFA_ANNOTATIONS = ("kegg_pathway", "ec_number", "kegg_reaction", "reactome_id")
# 代谢组派生特征（脂质类别总量、质量作用比）的差异结果，与 differential.parquet 同目录
DERIVED_TABLE = "derived_differential.parquet"
# 证据画像中各层与通量方向一致的特征
_CONCORDANT = {"tx_details": ("concordant_genes", "tx_concordant_features"),
               "prot_details": ("concordant_proteins", "prot_concordant_features"),
               "metab_details": ("concordant_metabolites", "metab_concordant_features")}


def _numeric(rows: pd.DataFrame, column: str) -> pd.Series:
    if column not in rows.columns:
        return pd.Series(np.nan, index=rows.index)
    return pd.to_numeric(rows[column], errors="coerce")


def _n_eff(rows: pd.DataFrame) -> pd.Series:
    return _numeric(rows, "n_numerator") + _numeric(rows, "n_denominator")


def mfa_frames(diff: pd.DataFrame, tissue_map: dict[str, str] | None = None) -> dict[str, pd.DataFrame]:
    """MFA 差异结果 → ``{比较: 每行一个（反应, 组织）}``，列名同旧流程的 MFA 差异通量表。"""
    tissue_map = tissue_map or {}
    frames: dict[str, pd.DataFrame] = {}
    for contrast, rows in diff.groupby("contrast_id", sort=True):
        effect_only = rows["method"].astype(str).eq("effect_size_only").to_numpy()
        low_n = (rows[["n_numerator", "n_denominator"]].min(axis=1) < 2).to_numpy()
        frame = pd.DataFrame({
            "flux_id": rows["feature_id"].astype(str),
            "tissue": rows["tissue"].astype(str).map(lambda t: tissue_map.get(t, t)),
            "flux_log2fc": _numeric(rows, "log2fc"),
            "flux_pvalue": _numeric(rows, "pvalue"),
            "flux_fdr": _numeric(rows, "padj"),
            "cohen_d": _numeric(rows, "cohen_d" if "cohen_d" in rows.columns else "effect_size"),
            "stat_note": np.select([low_n & ~effect_only, effect_only],
                                   ["low_n_configured_statistical_test", "fitting_variance_only"],
                                   "biological_replicates"),
            "flux_se": _numeric(rows, "log2fc_se"),
            "n_eff": _n_eff(rows),
            "selection_probability": _numeric(rows, "selection_probability"),
            "direction_stability": _numeric(rows, "direction_stability"),
            "lfsr": _numeric(rows, "lfsr"),
            # 锚点状态的依据（证据画像的 anchor_state）
            "mfa_evidence_level": rows["evidence_level"].astype(str) if "evidence_level" in rows.columns else "",
            "n_min_group": rows[["n_numerator", "n_denominator"]].min(axis=1),
            "mfa_flags": rows["flags"].fillna("").astype(str) if "flags" in rows.columns else "",
        })
        for column in ("identifiability", "tracer_evidence"):
            if column in rows.columns:
                frame[column] = rows[column].fillna("").astype(str).to_numpy()
        for column in _MFA_ANNOTATIONS:
            if column in rows.columns:
                frame[column] = rows[column]
        frames[str(contrast)] = frame.reset_index(drop=True)
    return frames


def transcriptomics_frames(diff: pd.DataFrame, id_types: dict[str, str] | None = None,
                           ) -> dict[str, dict[str, pd.DataFrame]]:
    """转录组差异结果 → ``{比较: {组织: 每行一个基因}}``。``id_types`` 为各组织数据集的 ``feature_id_type``。"""
    id_types = id_types or {}
    frames: dict[str, dict[str, pd.DataFrame]] = {}
    for (contrast, tissue), rows in diff.groupby(["contrast_id", "tissue"], sort=True):
        ids = rows["feature_id"].astype(str)
        names = rows["feature_name"].fillna(ids) if "feature_name" in rows.columns else ids
        frame = pd.DataFrame({
            "gene": names.astype(str),
            "feature_id": ids,                             # 反应档案按它回查统一结果表；匹配不读该列
            "log2fc": _numeric(rows, "log2fc"),
            "pvalue": _numeric(rows, "pvalue"),
            "padj": _numeric(rows, "padj"),
            "evidence_level": rows["evidence_level"] if "evidence_level" in rows.columns else "",
            "lfcSE": _numeric(rows, "lfc_se"),
            "n_eff": _n_eff(rows),
            "lfsr": _numeric(rows, "lfsr"),
            "base_mean": _numeric(rows, "base_mean"),      # 同工酶加权的丰度
        })
        id_column = _GENE_ID_COLUMNS.get(id_types.get(str(tissue), ""))
        if id_column:
            frame[id_column] = ids
        for column in ("n_batches", "batch_ids"):
            if column in rows.columns:
                frame[column] = rows[column]
        frames.setdefault(str(contrast), {})[str(tissue)] = frame.reset_index(drop=True)
    return frames


def _wide_frames(diff: pd.DataFrame, identity: dict[str, str], stats: dict[str, str]) -> dict[str, pd.DataFrame]:
    """每个组织一张宽表：``identity``（源列 → 列名）每个特征一个值，``stats``（源列 → 列名前缀）按比较加后缀。"""
    frames: dict[str, pd.DataFrame] = {}
    for tissue, rows in diff.groupby("tissue", sort=True):
        features = rows.drop_duplicates("feature_id").set_index("feature_id", drop=False)
        frame = pd.DataFrame({target: features[source] for source, target in identity.items()
                              if source in features.columns}, index=features.index)
        for contrast, part in rows.groupby("contrast_id", sort=True):
            part = part.set_index("feature_id")
            for source, prefix in stats.items():
                if source in part.columns:
                    frame[f"{prefix}_{contrast}"] = part[source].reindex(frame.index).to_numpy()
        frames[str(tissue)] = frame.reset_index(drop=True)
    return frames


def proteomics_frames(diff: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """蛋白组差异结果 → ``{组织: 每行一个蛋白}``，列同旧流程的蛋白组宽表（``log2fc_<比较>`` 等）。"""
    diff = diff.assign(n_eff=_n_eff(diff))
    # 同工酶加权的丰度：各比较两组均值（log2 强度）的平均，换回线性尺度
    level = (_numeric(diff, "mean_numerator") + _numeric(diff, "mean_denominator")) / 2
    diff["abundance"] = np.power(2.0, level.groupby([diff["tissue"], diff["feature_id"]]).transform("mean"))
    if "evidence_level" in diff.columns:
        # 旧流程的蛋白组宽表中，经过检验的蛋白记为 inferential
        diff["evidence"] = np.where(diff["evidence_level"].eq("statistical_test"), "inferential",
                                    diff["evidence_level"].astype(str))
    return _wide_frames(diff, {"feature_id": "protein", **_PROTEIN_IDS, "abundance": "abundance"},
                        {"log2fc": "log2fc", "pvalue": "pval", "padj": "padj", "lfc_se": "se", "n_eff": "n_eff",
                         "evidence": "evidence", "lfsr": "lfsr"})


def metabolomics_frames(diff: pd.DataFrame, derived: pd.DataFrame | None = None) -> dict[str, pd.DataFrame]:
    """代谢组差异结果 → ``{组织: 每行一个代谢物}``，列同旧流程的代谢组宽表（``log2fc_<比较>`` 等）。

    ``derived`` 为代谢组派生特征的差异结果（``derived_differential.parquet``：脂质类别总量、质量作用比），
    与测得代谢物并入同一张宽表，以 ``feature_type`` 与 ``related_reaction`` 列区分；``metabolite_id`` 为特征 ID。
    """
    if derived is not None and not derived.empty:
        diff = pd.concat([diff, derived], ignore_index=True)
    ids = diff["feature_id"].astype(str)
    names = diff["feature_name"].fillna(ids) if "feature_name" in diff.columns else ids
    diff = diff.assign(compound=names.astype(str), n_eff=_n_eff(diff))
    # 诱饵映射（E3）按丰度匹配：各比较两组均值（log2 归一化峰面积）的平均，换回线性尺度
    level = (_numeric(diff, "mean_numerator") + _numeric(diff, "mean_denominator")) / 2
    diff["abundance"] = np.power(2.0, level.groupby([diff["tissue"], diff["feature_id"]]).transform("mean"))
    identity = {"compound": "compound", **_METABOLITE_IDS, "feature_id": "metabolite_id", "abundance": "abundance",
                **{c: c for c in ("feature_type", "related_reaction") if c in diff.columns}}
    return _wide_frames(diff, identity, {"log2fc": "log2fc", "pvalue": "pval", "padj": "fdr", "n_eff": "n_eff",
                                         "lfsr": "lfsr"})


_ANNOTATION_IDS = (("hmdb_id", "hmdb"), ("kegg_compound", "kegg_compound"), ("chebi_id", "chebi"))


def refresh_metabolite_ids(diff: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """更正旧结果包中由错配 ID 表写入的代谢物 ID。

    代谢组运行器按名称查代谢物 ID 表写入 HMDB / KEGG / ChEBI 列，而整合按 ID 优先匹配。ID 表中的错配修正后
    （``metabolite_id_mapping.CORRECTED_IDS``），旧结果包里仍是错误的值，例如“Malic acid”写成 D-苹果酸，测得的苹果酸
    因此匹配不到苹果酸脱氢酶。这里只把等于已知错误旧值的 ID 换成修正值，其他注释（含厂商给出的 ID）不动。
    返回 (表, 有更正的行数)。
    """
    from omicsweave.core.idmap.metabolite_id_mapping import corrected_ids

    if diff.empty or "feature_name" not in diff.columns or not any(c in diff.columns for _, c in _ANNOTATION_IDS):
        return diff, 0
    names = diff["feature_name"].astype(str)
    fixes = {name: corrected_ids(name) for name in names.unique()}
    if not any(fixes.values()):
        return diff, 0
    out = diff.copy()
    changed = pd.Series(False, index=diff.index)
    for key, column in _ANNOTATION_IDS:
        if column not in out.columns:
            continue
        current = out[column].astype(object)
        wrong = pd.Series([fixes[n].get(key, ("", ""))[0] for n in names], index=out.index)
        right = pd.Series([fixes[n].get(key, ("", ""))[1] for n in names], index=out.index)
        hit = wrong.ne("") & current.astype(str).eq(wrong)
        out[column] = current.where(~hit, right)
        changed |= hit
    return out, int(changed.sum())


def _transcriptomics_id_types(source: Path) -> dict[str, str]:
    types = {}
    for path in sorted(source.glob(f"omics/{Omics.TRANSCRIPTOMICS.value}/data/*/{META_FILE}")):
        types[path.parent.name] = str(json.loads(path.read_text(encoding="utf-8")).get("feature_id_type", ""))
    return types


ROBUSTNESS_TABLE = "robustness.parquet"
_LAYER_OMICS = {"tx_details": ("transcriptomics", "concordant_genes"), "prot_details": ("proteomics", "concordant_proteins"),
                "metab_details": ("metabolomics", "concordant_metabolites")}


def single_subject_features(source: Path, tables: dict[str, pd.DataFrame]) -> tuple[set, list[Path]]:
    """源结果包中由单只动物驱动的入选特征（E5）：{(组学, 组织, 比较, 特征 ID 或名称)}。"""
    keys, paths = set(), sorted(source.glob(f"omics/*/tables/{ROBUSTNESS_TABLE}"))
    for path in paths:
        omics = path.parent.parent.name
        robust = pd.read_parquet(path)
        robust = robust[robust["single_subject_driven"]]
        names = {}
        diff = tables.get(omics)
        if diff is not None and "feature_name" in diff.columns:
            names = dict(zip(diff["feature_id"].astype(str), diff["feature_name"].astype(str)))
        for contrast, tissue, feature in zip(robust["contrast_id"], robust["tissue"], robust["feature_id"]):
            keys.add((omics, str(tissue), str(contrast), str(feature)))
            if str(feature) in names:
                keys.add((omics, str(tissue), str(contrast), names[str(feature)]))
    return keys, paths


def single_subject_layers(profile: pd.DataFrame, flagged: set) -> pd.Series:
    """每个反应中，与通量方向一致的特征全部由单只动物驱动的层（以分号分隔）。"""
    def layers(row: pd.Series) -> str:
        hit = []
        for details, (omics, key) in _LAYER_OMICS.items():
            features = (row.get(details) or {}).get(key, []) or []
            if features and all((omics, str(row["tissue"]), str(row["comparison"]), str(f)) in flagged
                                for f in features):
                hit.append(omics)
        return ";".join(hit)
    if profile.empty or not flagged:
        return pd.Series("", index=profile.index)
    return profile.apply(layers, axis=1)


def profile_table(profile: pd.DataFrame, roles: dict[str, str]) -> pd.DataFrame:
    """证据画像 → 可写入 parquet 的表：去掉各层明细（字典），补上比较角色与各层方向一致的特征。"""
    if profile.empty:
        return pd.DataFrame(columns=["flux_id", "tissue", "comparison", "contrast_role"])
    out = profile.drop(columns=[c for c in profile.columns if c.endswith("_details")])
    out.insert(int(out.columns.get_loc("comparison")) + 1, "contrast_role", out["comparison"].map(roles))
    for details, (key, column) in _CONCORDANT.items():
        out[column] = profile[details].map(lambda d: ";".join(str(f) for f in (d or {}).get(key, [])))
    # 同工酶加权（G7）：各层容量变化的加权值，以及支持只来自次要同工酶而被撤销的层
    for details, prefix in (("tx_details", "tx"), ("prot_details", "prot")):
        out[f"{prefix}_weighted_capacity_change"] = profile[details].map(
            lambda d: (d or {}).get("weighted_capacity_change", np.nan))
        out[f"{prefix}_major_member"] = profile[details].map(lambda d: (d or {}).get("major_member", ""))
    out["minor_isozyme_driven_layers"] = profile.apply(lambda r: ";".join(
        omics for details, omics in (("tx_details", "transcriptomics"), ("prot_details", "proteomics"))
        if (r[details] or {}).get("minor_isozyme_driven")), axis=1)
    return out


def scoring_config(params: dict[str, Any], mfa: pd.DataFrame) -> dict[str, Any]:
    """证据构建器的 ``scoring``；MFA 推断单位默认取差异结果中的 ``inference_unit``（各行一致时）。"""
    scoring = {k: v for k, v in params.items() if k not in _MODE_ONLY}
    profile = dict(scoring.get("evidence_profile") or {})
    units = mfa["inference_unit"].dropna().unique() if "inference_unit" in mfa.columns else []
    if "mfa_inference_unit" not in profile and len(units) == 1:
        profile["mfa_inference_unit"] = str(units[0])
    scoring["evidence_profile"] = profile
    return scoring


def run(ctx: RunContext, source: Path, parameters: dict[str, Any]) -> IntegrationRun:
    params = {**DEFAULTS, **parameters}
    unknown = sorted(set(parameters) - set(DEFAULTS))
    if unknown:
        raise RunError(f"integration.parameters.{MODE.value} 中有无法识别的参数：{'、'.join(unknown)}"
                       f"（可用：{'、'.join(DEFAULTS)}）")
    paths = {p.parent.parent.name: p for p in sorted(source.glob("omics/*/tables/differential.parquet"))}
    if Omics.MFA.value not in paths or len(paths) < 2:
        found = "、".join(paths) or "没有"
        raise RunError(f"{MODE.value} 需要 MFA 与至少 1 个其他组学的差异结果；结果包 {source.name} 中有：{found}。"
                       f"请先对启用的组学运行 omicsweave run（含 differential 步骤）")
    tables = {omics: read_differential(path) for omics, path in paths.items()}
    reannotated = 0
    if Omics.METABOLOMICS.value in tables:
        tables[Omics.METABOLOMICS.value], reannotated = refresh_metabolite_ids(tables[Omics.METABOLOMICS.value])
    derived_path = source / "omics" / Omics.METABOLOMICS.value / "tables" / DERIVED_TABLE
    derived = read_differential(derived_path) if derived_path.is_file() else None
    empty = pd.DataFrame(columns=["contrast_id", "tissue", "feature_id", "evidence_level"])

    tissue_map = {str(k): str(v) for k, v in (params["tissue_map"] or {}).items()}
    mfa = mfa_frames(tables[Omics.MFA.value], tissue_map)
    tx = transcriptomics_frames(tables.get(Omics.TRANSCRIPTOMICS.value, empty), _transcriptomics_id_types(source))
    prot = proteomics_frames(tables.get(Omics.PROTEOMICS.value, empty))
    metab = metabolomics_frames(tables.get(Omics.METABOLOMICS.value, empty), derived)

    # 样本 QC 标记 → 各层（组织, 比较）的可靠性；按当前样本表计算，写入本次整合的结果包
    reliability, quality_inputs = build_reliability(ctx.experiment, source, tables, float(params["qc_imbalance"]))
    qc_dir = mode_dir(ctx, MODE) / QC_DIR
    qc_dir.mkdir(parents=True)
    reliability.to_parquet(qc_dir / RELIABILITY_FILE, index=False)
    mapped = reliability.copy()
    is_mfa = mapped["omics"] == Omics.MFA.value
    mapped.loc[is_mfa, "tissue"] = mapped.loc[is_mfa, "tissue"].map(lambda t: tissue_map.get(t, t))

    builder = FluxEvidenceBuilder(SimpleNamespace(experiment_id=ctx.experiment.id, estimand={},
                                                  scoring=scoring_config(params, tables[Omics.MFA.value]),
                                                  contrast_reliability=mapped))
    profile = builder.build_all(mfa, tx, prot, metab)
    flagged, robustness_inputs = single_subject_features(source, tables)
    if not profile.empty:
        profile["single_subject_driven_layers"] = single_subject_layers(profile, flagged)
    roles = dict(zip(tables[Omics.MFA.value]["contrast_id"].astype(str),
                     tables[Omics.MFA.value]["contrast_role"].astype(str)))
    out = mode_dir(ctx, MODE) / "tables"
    out.mkdir(parents=True)
    profile_table(profile, roles).to_parquet(out / "flux_evidence_profile.parquet", index=False)
    builder.get_evidence_table().to_parquet(out / "flux_evidence_table.parquet", index=False)
    builder.aggregate_by_pathway(profile).to_parquet(out / "flux_pathway_profile.parquet", index=False)

    decoys = {**DEFAULTS["decoys"], **(params["decoys"] or {})}
    params["decoys"] = decoys
    decoy_summary = None
    if int(decoys["replicates"]) > 0:
        from ..evaluation.decoys import decoy_evaluation

        detail, decoy_summary = decoy_evaluation(builder, mfa, tx, prot, metab, builder.get_evidence_table(),
                                                 replicates=int(decoys["replicates"]), seed=int(decoys["seed"]))
        evaluation = mode_dir(ctx, MODE) / "evaluation"
        evaluation.mkdir(parents=True)
        detail.to_parquet(evaluation / "decoy_evidence.parquet", index=False)
        decoy_summary.to_parquet(evaluation / "decoy_summary.parquet", index=False)

    record = IntegrationRun(MODE, inputs=list(paths.values()) + ([derived_path] if derived is not None else [])
                            + quality_inputs + robustness_inputs,
                            parameters={**params, "omics": sorted(paths),
                                        "derived_metabolomics_features": derived is not None})
    note = summary_note(reliability)
    if note:
        record.notes.append(note)
    if reannotated:
        record.notes.append(f"代谢组差异结果中 {reannotated} 行的代谢物 ID 是代谢物 ID 表修正前写入的错配值，"
                            "已按修正值匹配（重新运行代谢组后不再需要）")
    record.parameters["metabolite_ids_reannotated"] = reannotated
    if decoy_summary is not None and not decoy_summary.empty:
        rows = decoy_summary[decoy_summary["subset"] == "all_directional"]
        record.notes.append("诱饵映射（E3）的假证据比例：" + "；".join(
            f"{r.layer} 目标支持 {r.target_support_rate:.1%}、诱饵支持 {r.decoy_support_rate:.1%}"
            for r in rows.itertuples()) + "；明细见 evaluation/")
    other_tissues = {t for frames in tx.values() for t in frames} | set(prot) | set(metab)
    unmatched = sorted({t for frame in mfa.values() for t in frame["tissue"]} - other_tissues)
    if unmatched:
        record.notes.append(f"MFA 组织 {'、'.join(unmatched)} 在其他组学中没有同名组织，只有 MFA 层证据；"
                            f"名称不同时用 integration.parameters.{MODE.value}.tissue_map 对应")
    contrasts = {c for omics, t in tables.items() if omics != Omics.MFA.value for c in t["contrast_id"].astype(str)}
    missing = sorted(set(mfa) - contrasts)
    if missing:
        record.notes.append(f"比较 {'、'.join(missing)} 只有 MFA 差异结果，其他组学没有对应比较")
    return record
