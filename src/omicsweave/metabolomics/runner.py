#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
代谢组的统一运行器（``omicsweave run --omics metabolomics``）：读取 ``experiment.yaml`` 声明的峰面积表，
按组织预处理，再按其中的比较做差异分析与 KEGG 通路 ORA，写出统一结果表。

计算沿用旧流程的函数（预处理、PLS-DA、Welch t 检验、ORA），只是分组与比较改由样本表和
``design.contrasts`` 给出。样本按数据表中的列顺序参与计算，所以同一份数据与旧流程
（``omicsweave metabolomics run``）的结果逐值一致。富集分析默认仅执行本地 KEGG ORA；代谢物集预排序 GSEA（``gsea_rank_metric`` 默认 Welch t，种子固定为 42）由 ``analysis.run_gsea`` 控制（默认 ``false``，可通过配置或命令行 ``--gsea`` 开启）。旧流程的图件暂未接入。

输入（``omics.metabolomics.input``，``format: peak_table``）::

    tissue_files: {liver: data/liver.xlsx}   # 每个组织一张表（旧版布局）；没有时读取 files，
    files: [data/peaks.tsv]                  #   一张表中可以有多个组织的样本，组织由 samples.tsv 给出
    feature_id_column: Index                 # 可选：特征 ID 列；默认 Index，没有时用名称列
    feature_name_column: Compounds           # 可选：代谢物名称列，用于 KEGG / HMDB 注释

数据表中的样本列名即 samples.tsv 的 ``sample_id``，不同时在 ``data_column`` 列写出数据表中的列名。
``sample_type`` 为 qc / pool 的列参与 QC 过滤（RSD），blank / standard 不进入分析。

输出（``omics/metabolomics/``）::

    data/<组织>/                  预处理后的 log2 归一化峰面积（分析层数据集；只含生物学样本）
    qc/<组织>/                    逐代谢物 QC 指标、样本总信号与归一化因子、预处理摘要
    derived/<组织>/               派生特征（脂质类别总量、质量作用比；log2，定义见 features.parquet），
                                  见 :mod:`omicsweave.metabolomics.derived_features`
    tables/differential.parquet   统一差异结果表
    tables/derived_differential.parquet
                                  派生特征的差异结果（同一 Welch t 检验；按组织 × 比较 × 特征类型单独做 BH 校正，
                                  不并入测得代谢物的 FDR，不要求 VIP）。由 analysis 的 derived_features 开关控制
    tables/enrichment.parquet     KEGG 通路 ORA（显著代谢物，不分上下调）与预排序 GSEA
    figures/                      火山图 volcano_<组织>_<比较>.pdf / .png（``figures`` 步骤，见
                                  :mod:`omicsweave.visualization.volcano`；纵轴为校正后 p）
    details/<组织>/               旧版格式的差异分析宽表与 PLS-DA 验证摘要
"""
from __future__ import annotations

import json
from dataclasses import MISSING, dataclass, fields, replace
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from omicsweave.adapters.enrichment import resolve_run_gsea
from omicsweave.adapters.gseapy_results import prerank_to_contract
from omicsweave.contract.dataset import Dataset, DatasetMeta, ValueType, write_dataset
from omicsweave.contract.tables import empty_enrichment, write_differential, write_enrichment
from omicsweave.core.config import Experiment, Project
from omicsweave.core.idmap.metabolite_id_mapping import annotate_metabolites
from omicsweave.core.model import Omics
from omicsweave.core.run import OmicsRun, RunContext, RunError, contrast_groups, sample_labels
from omicsweave.core.tissues import format_alias_note, select_tissue_file_rows
from omicsweave.visualization.volcano import Thresholds, write_bundle_volcanoes

from .derived_features import UNITS, DerivedFeatures, derive_features
from .data_preprocessing import (
    PreprocessingOptions,
    classify_data_columns,
    normalize_fdr_method,
    preprocess_metabolomics_data,
)
from .differential_analysis import DifferentialAnalysisOptions, analyze_differential_metabolites
from .enrichment._run import calculate_pathway_enrichment
from .experiment_config import MetabolomicsConfig
from .export import derived_differential_to_contract, differential_to_contract, ora_to_contract
from .gsea._run import run_gsea
from .plotting.qc_plots import fit_plsda

STEPS = ("qc", "differential", "enrichment", "figures")
_QC_TYPES = {"qc", "pool"}
_EXCLUDED_TYPES = {"blank", "standard"}
_DISPLAY_ONLY = {"tissue_labels", "plot_colors"}
# 分析参数的默认值与旧版配置加载器相同
_DEFAULTS = {f.name: f.default for f in fields(MetabolomicsConfig)
             if f.default is not MISSING and f.name not in _DISPLAY_ONLY}


@dataclass
class _Block:
    """一个组织的峰面积表：数据列已改名为 sample_id，样本按数据表中的列顺序排列。"""
    tissue: str
    path: Path
    table: pd.DataFrame
    biological: list[str]
    qc: list[str]
    annotations: list[str]
    paths: tuple[Path, ...] = ()

    def __post_init__(self) -> None:
        if not self.paths:
            self.paths = (self.path,)


def _merge_blocks(existing: _Block, incoming: _Block, spec) -> _Block:
    """将同一规范组织来自多张数据表的峰面积表按特征 ID 合并为单一块。"""
    overlap = [s for s in (*incoming.biological, *incoming.qc) if s in set(existing.biological) | set(existing.qc)]
    if overlap:
        raise RunError(f"代谢组 {existing.tissue} 的样本 {'、'.join(overlap[:5])} 在多张数据表中重复出现")
    fid_left, _ = _feature_columns(existing, spec)
    fid_right, _ = _feature_columns(incoming, spec)
    if fid_left != fid_right:
        raise RunError(
            f"代谢组 {existing.tissue} 合并多张数据表时特征 ID 列不一致："
            f"{existing.path.name} 为 {fid_left!r}，{incoming.path.name} 为 {fid_right!r}"
        )
    left_ids = existing.table[fid_left].astype(str)
    right_ids = incoming.table[fid_right].astype(str)
    for label, ids in ((existing.path.name, left_ids), (incoming.path.name, right_ids)):
        if ids.duplicated().any():
            raise RunError(f"代谢组 {existing.tissue}（{label}）的特征 ID 列 {fid_left} 有重复："
                           f"{'、'.join(ids[ids.duplicated()][:5])}")
    left = existing.table.set_axis(left_ids)
    right = incoming.table.set_axis(right_ids)
    ordered = pd.Index(dict.fromkeys([*left.index, *right.index]), name="feature_id")
    annotations = list(dict.fromkeys([*existing.annotations, *incoming.annotations]))
    cols: dict[str, pd.Series] = {}
    for col in annotations:
        if col == fid_left:
            cols[col] = pd.Series(ordered.to_numpy(), index=ordered)
        elif col in left.columns and col in right.columns:
            primary = left[col].reindex(ordered)
            fallback = right[col].reindex(ordered)
            empty = primary.isna() | primary.astype(str).str.strip().eq("")
            cols[col] = primary.where(~empty, fallback)
        elif col in left.columns:
            cols[col] = left[col].reindex(ordered)
        else:
            cols[col] = right[col].reindex(ordered)
    left_samples = [c for c in existing.table.columns if c not in existing.annotations]
    right_samples = [c for c in incoming.table.columns if c not in incoming.annotations]
    for col in left_samples:
        cols[col] = left[col].reindex(ordered)
    for col in right_samples:
        cols[col] = right[col].reindex(ordered)
    merged = pd.DataFrame(cols, index=ordered).reset_index(drop=True)
    return _Block(
        tissue=existing.tissue,
        path=existing.path,
        table=merged,
        biological=[*existing.biological, *incoming.biological],
        qc=[*existing.qc, *incoming.qc],
        annotations=annotations,
        paths=(*existing.paths, *incoming.paths),
    )


def _settings(analysis: dict[str, Any]) -> tuple[dict[str, Any], list[str]]:
    settings = {**_DEFAULTS, **analysis}
    if "fdr_cutoff" not in analysis and "pval_cutoff" in analysis:
        # 与旧版加载器相同：未写 fdr_cutoff 时取 pval_cutoff
        settings["fdr_cutoff"] = settings["pval_cutoff"]
    if str(settings["batch_correction"]).strip().lower() != "none":
        raise RunError("代谢组流程尚不支持批次校正，请设 batch_correction: none，平衡设计或分批分析后做结果级整合")
    unknown = sorted(set(analysis) - set(_DEFAULTS) - _DISPLAY_ONLY)
    return settings, ([f"omics.metabolomics.analysis 中有未使用的参数：{'、'.join(unknown)}"] if unknown else [])


def _read_table(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise RunError(f"找不到代谢组数据表 {path}")
    suffix = path.suffix.lower()
    if suffix in (".xlsx", ".xls"):
        return pd.read_excel(path)
    return pd.read_csv(path, sep="," if suffix == ".csv" else "\t")


def _blocks(project: Project, experiment: Experiment, spec) -> tuple[list[_Block], list[Path], list[str], dict[str, str]]:
    rows = experiment.samples[experiment.samples["omics"] == Omics.METABOLOMICS.value].copy()
    rows["column"] = rows["data_column"].where(rows["data_column"] != "", rows["sample_id"]) \
        if "data_column" in rows.columns else rows["sample_id"]
    rows["type"] = rows["sample_type"].replace("", "sample") if "sample_type" in rows.columns else "sample"
    rows = rows[~rows["type"].isin(_EXCLUDED_TYPES)]
    declared = set(rows["column"])
    applied_aliases = experiment.applied_sample_aliases(Omics.METABOLOMICS)

    tissue_files = getattr(spec, "tissue_files", None) or {}
    if tissue_files:
        resolved_jobs = [
            (str(t), experiment.config.canonical_tissue(str(t), Omics.METABOLOMICS), project.resolve(p))
            for t, p in tissue_files.items()
        ]
        canonical_counts: dict[str, int] = {}
        for raw_key, canon_key, _ in resolved_jobs:
            canonical_counts[canon_key] = canonical_counts.get(canon_key, 0) + 1
            if raw_key != canon_key:
                applied_aliases.setdefault(raw_key, canon_key)
    else:
        resolved_jobs = [(None, None, project.resolve(f)) for f in spec.files]
        canonical_counts = {}
    blocks: dict[str, _Block] = {}
    notes: list[str] = []
    for raw_key, canon_key, path in resolved_jobs:
        table = _read_table(path)
        table.columns = [str(c) for c in table.columns]
        if canon_key is None:
            candidates = rows
        else:
            candidates = select_tissue_file_rows(
                rows, raw_key, canon_key,
                multi_key_canonical=canonical_counts[canon_key] > 1,
                raw_tissues=experiment.raw_tissues,
            )
        present = candidates[candidates["column"].isin(table.columns)]
        metadata = set(classify_data_columns(table.columns).metadata)
        annotations = [c for c in table.columns if c not in declared and c in metadata]
        ignored = [c for c in table.columns if c not in declared and c not in metadata]
        if ignored:
            notes.append(f"{path.name} 中有 {len(ignored)} 列像样本但不在 samples.tsv 中，未参与分析：{'、'.join(ignored[:5])}")
        for block_tissue, part in present.groupby("tissue", sort=False):
            if part["column"].duplicated().any():
                raise RunError(f"代谢组 {block_tissue} 有多个样本对应数据表中的同一列："
                               f"{'、'.join(part.loc[part['column'].duplicated(), 'column'][:5])}")
            by_column = part.set_index("column")
            order = [c for c in table.columns if c in by_column.index]
            biological = [by_column.at[c, "sample_id"] for c in order if by_column.at[c, "type"] not in _QC_TYPES]
            qc = [by_column.at[c, "sample_id"] for c in order if by_column.at[c, "type"] in _QC_TYPES]
            renamed = table[annotations + order].rename(columns=dict(zip(order, by_column.loc[order, "sample_id"])))
            incoming = _Block(block_tissue, path, renamed, biological, qc, annotations)
            blocks[block_tissue] = _merge_blocks(blocks[block_tissue], incoming, spec) if block_tissue in blocks else incoming
    for block in blocks.values():
        if len(block.paths) > 1:
            names = "、".join(p.name for p in block.paths)
            notes.append(f"代谢组 {block.tissue}：已合并来自 {len(block.paths)} 张数据表的样本（{names}，共 {len(block.biological)} 个生物学样本）")
    found_samples = {s for block in blocks.values() for s in (*block.biological, *block.qc)}
    for tissue, part in rows.groupby("tissue", sort=False):
        if tissue not in blocks:
            continue
        absent = [s for s in part["sample_id"] if s not in found_samples]
        if absent:
            notes.append(f"代谢组 {tissue} 的数据表中找不到 {len(absent)} 个样本：{'、'.join(absent[:5])}")
    unread = sorted(set(rows["tissue"]) - set(blocks))
    if unread:
        notes.append(f"samples.tsv 中以下组织的代谢组样本没有对应的数据表：{'、'.join(unread)}")
    return list(blocks.values()), [path for _, _, path in resolved_jobs], notes, applied_aliases


def _feature_columns(block: _Block, spec) -> tuple[str, str | None]:
    name = getattr(spec, "feature_name_column", None) or ("Compounds" if "Compounds" in block.annotations else None)
    feature_id = getattr(spec, "feature_id_column", None) or ("Index" if "Index" in block.annotations else name)
    for label, column in (("feature_id_column", feature_id), ("feature_name_column", name)):
        if column is not None and column not in block.annotations:
            raise RunError(f"{block.path.name} 中没有 {label} 指定的列 {column!r}")
    if feature_id is None:
        raise RunError(f"{block.path.name} 中没有特征 ID 列（Index 或 Compounds），请在 input.feature_id_column 中指定")
    return feature_id, name


def _vip(matrix: pd.DataFrame, samples: list[str], labels: list[str],
         settings: dict[str, Any]) -> tuple[np.ndarray, dict[str, Any]]:
    """PLS-DA 的 VIP；模型未通过交叉验证与置换检验时为空（与旧流程相同）。"""
    result = fit_plsda(matrix[samples].T.astype(float), labels, n_components=2,
                       n_permutations=int(settings["plsda_permutation_num"]),
                       random_state=int(settings["random_seed"]))
    values = result.vip.to_numpy(dtype=float)
    if result.status != "validated":
        values[:] = np.nan
    return values, {"status": result.status, "R2Y": result.r2y, "Q2": result.q2,
                    "CV_Accuracy": result.cv_accuracy, "Permutation_P": result.permutation_pvalue}


def _features(processed: pd.DataFrame, annotations: list[str], ids: pd.Index) -> pd.DataFrame:
    features = processed[annotations].copy()
    for column in features.columns:
        if features[column].dtype == object:
            features[column] = features[column].astype("string")
    features.index = ids
    return features


def _gsea(combined: pd.DataFrame, feature_id: str, contrast_id: str, tissue: str, details: Path,
          settings: dict[str, Any]) -> pd.DataFrame | None:
    """旧流程的代谢物集预排序 GSEA；它以 Index 列为代谢物 ID，特征 ID 列不同时先补上。"""
    table = combined if feature_id == "Index" else combined.assign(Index=combined[feature_id].astype(str))
    result = run_gsea(table, contrast_id, tissue, str(details / "gsea"),
                      min_set_size=int(settings["gsea_min_set_size"]), max_set_size=int(settings["gsea_max_set_size"]),
                      permutation_num=int(settings["gsea_permutation_num"]),
                      rank_metric=str(settings["gsea_rank_metric"]), pathway_scope=str(settings["pathway_scope"]),
                      exclude_overview=bool(settings["exclude_overview_pathways"]))
    return None if result is None or result.empty else result


def _split_kegg_terms(rows: pd.DataFrame) -> pd.DataFrame:
    """旧流程的 GSEA 条目名为 ``ko00010 | 通路名``，拆成 term_id 与 term_name。"""
    parts = rows["term_name"].str.split(" | ", n=1, regex=False)
    return rows.assign(term_id=parts.str[0], term_name=parts.str[1].fillna(rows["term_name"]))


def _derived_differential(derived: DerivedFeatures, tissue: str, groups: dict[str, str], contrast, options, *,
                          fdr_method: str, settings: dict[str, Any]) -> list[pd.DataFrame]:
    """派生特征的一个比较：每种特征类型单独检验与校正（与测得代谢物同一个 Welch t 检验，不要求 VIP）。"""
    levels = (contrast.numerator, contrast.denominator)
    members = {level: [s for s, g in groups.items() if g == level] for level in levels}
    derived_options = replace(options, require_vip_for_significance=False)
    out = []
    for feature_type, definitions in derived.features.groupby("feature_type", sort=True):
        ids = definitions.index
        values = derived.matrix.loc[ids, list(groups)]
        table = pd.concat([pd.DataFrame({"derived_id": ids, "derived_name": definitions["feature_name"].to_numpy()},
                                        index=ids), values], axis=1)
        result = analyze_differential_metabolites(table, groups, [(*levels, contrast.id)], options=derived_options,
                                                  metadata_columns=["derived_id", "derived_name"])
        means = {level: values[cols].mean(axis=1) for level, cols in members.items()}
        rule = (f"padj<{settings['fdr_cutoff']} and |log2FC|>{settings['log2fc_cutoff']}"
                f"（{feature_type} 单独做 BH 校正，不要求 VIP）")
        out.append(derived_differential_to_contract(
            result.contrasts[contrast.id], definitions, units=UNITS[feature_type], selection_rule=rule,
            contrast_id=contrast.id, contrast_role=contrast.role.value, tissue=tissue,
            mean_numerator=means[contrast.numerator], mean_denominator=means[contrast.denominator],
            n_numerator=len(members[contrast.numerator]), n_denominator=len(members[contrast.denominator]),
            fdr_method=fdr_method))
    return out


def run(ctx: RunContext) -> OmicsRun:
    experiment = ctx.experiment
    spec = experiment.config.omics.metabolomics.input
    out = ctx.omics_dir(Omics.METABOLOMICS)
    record = OmicsRun(Omics.METABOLOMICS)
    settings, record.notes = _settings(experiment.config.omics.metabolomics.analysis)
    try:
        do_gsea = resolve_run_gsea(experiment.config.omics.metabolomics.analysis, override=ctx.run_gsea)
    except ValueError as exc:
        raise RunError(str(exc)) from exc
    settings["run_gsea"] = do_gsea
    blocks, record.inputs, notes, applied_aliases = _blocks(ctx.project, experiment, spec)
    record.notes += notes
    if applied_aliases:
        record.notes.append(format_alias_note(Omics.METABOLOMICS.value, applied_aliases))
    if not blocks:
        raise RunError("代谢组数据表中没有找到 samples.tsv 声明的样本")

    fdr_method = str(settings["fdr_method"])
    rule = f"padj<{settings['fdr_cutoff']} and |log2FC|>{settings['log2fc_cutoff']}"
    if settings["require_vip_for_significance"]:
        rule += f" and VIP>={settings['vip_cutoff']}"
    options = DifferentialAnalysisOptions(
        pvalue_cutoff=float(settings["pval_cutoff"]), fdr_cutoff=float(settings["fdr_cutoff"]),
        log2fc_cutoff=float(settings["log2fc_cutoff"]), fdr_method=fdr_method,
        vip_cutoff=float(settings["vip_cutoff"]),
        require_vip_for_significance=bool(settings["require_vip_for_significance"]),
        min_observations_per_group=int(settings["min_observations_per_group"]),
    )
    want_differential = ctx.wants("differential") or ctx.wants("enrichment")
    differential, enrichment, contrasts_run, derived_rows = [], [], [], []
    attempted_enrichment = False
    derived_counts: dict[str, dict[str, int]] = {}

    for block in blocks:
        tissue = block.tissue
        ctx.progress(Omics.METABOLOMICS, "preprocess", tissue=tissue)
        if not block.biological:
            record.notes.append(f"代谢组 {tissue} 没有生物学样本，已跳过")
            continue
        feature_id, name = _feature_columns(block, spec)
        labels = sample_labels(experiment, block.biological)
        pre = preprocess_metabolomics_data(block.table, block.biological, labels, qc_cols=block.qc,
                                           options=PreprocessingOptions(
                                               normalization_method=str(settings["normalization_method"]),
                                               detection_rate=float(settings["detection_rate"]),
                                               imputation_fraction=float(settings["imputation_fraction"]),
                                               qc_filter_enabled=bool(settings["qc_filter_enabled"]),
                                               qc_rsd_cutoff=float(settings["qc_rsd_cutoff"]),
                                               qc_min_detection_rate=float(settings["qc_min_detection_rate"])))
        processed = pre.processed_df
        ids = pd.Index(processed[feature_id].astype(str), name="feature_id")
        if ids.has_duplicates:
            raise RunError(f"代谢组 {tissue} 的特征 ID 列 {feature_id} 有重复：{'、'.join(ids[ids.duplicated()][:5])}")

        matrix = processed[block.biological].set_axis(ids)
        source_file = block.path.name if len(block.paths) == 1 else "; ".join(p.name for p in block.paths)
        meta = DatasetMeta(omics=Omics.METABOLOMICS, tissue=tissue, value_type=ValueType.PEAK_AREA, scale="log2",
                           normalization=pre.summary["normalization_method"], feature_id_type=feature_id,
                           source={"file": source_file, "detection_rule": pre.summary["detection_rule"]})
        write_dataset(out / "data" / tissue, Dataset(matrix, _features(processed, block.annotations, ids), meta))

        derived = None
        if settings["derived_features"]:
            derived = derive_features(processed, block.biological, feature_id_column=feature_id, name_column=name,
                                      min_species=int(settings["lipid_class_min_species"]))
            record.notes += [f"代谢组 {tissue} 的派生特征：{note}" for note in derived.notes]
            if derived.empty:
                derived = None
            else:
                derived_counts[tissue] = {str(k): int(v) for k, v in
                                          derived.features["feature_type"].value_counts().sort_index().items()}
                write_dataset(out / "derived" / tissue, Dataset(derived.matrix, derived.features, DatasetMeta(
                    omics=Omics.METABOLOMICS, tissue=tissue, value_type=ValueType.DERIVED, scale="log2",
                    normalization=pre.summary["normalization_method"], feature_id_type="derived_feature",
                    source={"from": f"data/{tissue}", "module": "omicsweave.metabolomics.derived_features"})))

        if ctx.wants("qc"):
            ctx.progress(Omics.METABOLOMICS, "qc", tissue=tissue)
            qc_dir = out / "qc" / tissue
            qc_dir.mkdir(parents=True, exist_ok=True)
            # QC 指标覆盖检出过滤后、QC 过滤前的全部代谢物，按原表行号取特征 ID
            metrics = pre.qc_metrics.copy()
            metrics.insert(0, "feature_id", block.table.loc[metrics.index, feature_id].astype(str).to_numpy())
            metrics.to_parquet(qc_dir / "feature_metrics.parquet", index=False)
            pre.sample_totals.reset_index().to_parquet(qc_dir / "sample_totals.parquet", index=False)
            (qc_dir / "preprocessing_summary.json").write_text(
                json.dumps(pre.summary, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
        if not want_differential:
            continue

        annotated = annotate_metabolites(processed, name_col=name) if name else processed
        metadata_columns = [c for c in annotated.columns if c not in block.biological and c not in block.qc]
        validation = []
        global_vip = None
        if len(set(labels.values())) >= 2:
            global_vip, info = _vip(processed, block.biological, [labels[s] for s in block.biological], settings)
            validation.append({"comparison": "global_multiclass", **info})
        details = out / "details" / tissue
        details.mkdir(parents=True, exist_ok=True)
        for contrast in experiment.config.design.contrasts:
            ctx.progress(Omics.METABOLOMICS, "differential", tissue=tissue, contrast_id=contrast.id)
            assigned = contrast_groups(experiment, contrast, Omics.METABOLOMICS)
            groups = {s: assigned[s] for s in block.biological if s in assigned}
            levels = (contrast.numerator, contrast.denominator)
            if not all(level in groups.values() for level in levels):
                if groups:
                    record.notes.append(f"比较 {contrast.id} 在代谢组 {tissue} 中有一组没有样本，已跳过")
                continue
            vip, info = _vip(processed, list(groups), list(groups.values()), settings)
            validation.append({"comparison": contrast.id, **info})
            result = analyze_differential_metabolites(annotated, groups, [(*levels, contrast.id)], options=options,
                                                      contrast_vip={contrast.id: vip}, global_vip=global_vip,
                                                      metadata_columns=metadata_columns)
            result.combined.to_csv(details / f"differential_{contrast.id}.csv", index=False, encoding="utf-8-sig")
            members = {level: [s for s, g in groups.items() if g == level] for level in levels}
            means = {level: processed[cols].mean(axis=1) for level, cols in members.items()}
            differential.append(differential_to_contract(
                result.contrasts[contrast.id], contrast_id=contrast.id, contrast_role=contrast.role.value,
                tissue=tissue, feature_id_column=feature_id, feature_name_column=name,
                mean_numerator=means[contrast.numerator], mean_denominator=means[contrast.denominator],
                n_numerator=len(members[contrast.numerator]), n_denominator=len(members[contrast.denominator]),
                fdr_method=fdr_method, selection_rule=rule))
            contrasts_run.append(contrast.id)
            if derived is not None:
                derived_rows += _derived_differential(derived, tissue, groups, contrast, options,
                                                      fdr_method=fdr_method, settings=settings)
            if ctx.wants("enrichment") and "kegg_map" in result.combined.columns:
                ctx.progress(Omics.METABOLOMICS, "enrichment", tissue=tissue, contrast_id=contrast.id)
                attempted_enrichment = True
                ora = calculate_pathway_enrichment(
                    result.combined, significance_col=f"Significant_{contrast.id}", feature_col=feature_id,
                    min_set_size=int(settings["ora_min_set_size"]), max_set_size=int(settings["ora_max_set_size"]),
                    pathway_scope=str(settings["pathway_scope"]),
                    exclude_overview=bool(settings["exclude_overview_pathways"]), fdr_method=fdr_method)
                if not ora.empty:
                    enrichment.append(ora_to_contract(ora, contrast_id=contrast.id, contrast_role=contrast.role.value,
                                                      tissue=tissue, query="all", gene_set_library="KEGG",
                                                      padj_method=normalize_fdr_method(fdr_method)))
                if do_gsea:
                    gsea = _gsea(result.combined, feature_id, contrast.id, tissue, details, settings)
                    if gsea is not None:
                        enrichment.append(prerank_to_contract(
                            gsea, contrast_id=contrast.id, contrast_role=contrast.role.value, tissue=tissue,
                            omics=Omics.METABOLOMICS.value, gene_set_library="KEGG",
                            permutation_num=int(settings["gsea_permutation_num"])).pipe(_split_kegg_terms))
        pd.DataFrame(validation).to_csv(details / "plsda_validation.csv", index=False, encoding="utf-8-sig")

    if differential:
        write_differential(pd.concat(differential, ignore_index=True), out / "tables" / "differential.parquet")
    if derived_rows:
        write_differential(pd.concat(derived_rows, ignore_index=True),
                           out / "tables" / "derived_differential.parquet")
    if enrichment:
        write_enrichment(pd.concat(enrichment, ignore_index=True), out / "tables" / "enrichment.parquet")
    elif ctx.wants("enrichment") and differential:
        write_enrichment(empty_enrichment(), out / "tables" / "enrichment.parquet")
        if attempted_enrichment:
            record.notes.append("代谢组富集分析没有可报告的结果（入选或排序代谢物过少，无显著 KEGG 通路）")
        else:
            record.notes.append("代谢组数据表没有 kegg_map 列或没有可检验的 KEGG 通路，未做 ORA 与 GSEA")
    if ctx.wants("figures"):
        ctx.progress(Omics.METABOLOMICS, "figures")
        # 要求 VIP 时，阈值线外的点也可能未入选；入选规则写在图下方
        record.notes += write_bundle_volcanoes(
            out, experiment.config.design.contrasts, metric="padj",
            thresholds=Thresholds(p=float(settings["fdr_cutoff"]), log2fc=float(settings["log2fc_cutoff"])))
    record.parameters = {
        "tissues": [b.tissue for b in blocks],
        "contrasts": list(dict.fromkeys(contrasts_run)),
        "derived_features": derived_counts,
        "analysis": settings,
    }
    if applied_aliases:
        record.parameters["tissue_aliases_applied"] = applied_aliases
    return record
