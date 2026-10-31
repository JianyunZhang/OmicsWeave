#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
代谢组的统一运行器（``omicsweave run --omics metabolomics``）：读取 ``experiment.yaml`` 声明的峰面积表，
按组织预处理，再按其中的比较做差异分析与 KEGG 通路 ORA，写出统一结果表。

计算沿用旧流程的函数（预处理、PLS-DA、Welch t 检验、ORA），只是分组与比较改由样本表和
``design.contrasts`` 给出。样本按数据表中的列顺序参与计算，所以同一份数据与旧流程
（``omicsweave metabolomics run``）的结果逐值一致。旧流程的图件与 GSEA 暂未接入。

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
    tables/differential.parquet   统一差异结果表
    tables/enrichment.parquet     KEGG 通路 ORA（显著代谢物，不分上下调）
    details/<组织>/               旧版格式的差异分析宽表与 PLS-DA 验证摘要
"""
from __future__ import annotations

import json
from dataclasses import MISSING, dataclass, fields
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from omicsweave.contract.dataset import Dataset, DatasetMeta, ValueType, write_dataset
from omicsweave.contract.tables import write_differential, write_enrichment
from omicsweave.core.config import Experiment, Project
from omicsweave.core.idmap.metabolite_id_mapping import annotate_metabolites
from omicsweave.core.model import Omics
from omicsweave.core.run import OmicsRun, RunContext, RunError, contrast_groups, sample_labels

from .data_preprocessing import (
    PreprocessingOptions,
    classify_data_columns,
    normalize_fdr_method,
    preprocess_metabolomics_data,
)
from .differential_analysis import DifferentialAnalysisOptions, analyze_differential_metabolites
from .enrichment._run import calculate_pathway_enrichment
from .experiment_config import MetabolomicsConfig
from .export import differential_to_contract, ora_to_contract
from .plotting.qc_plots import fit_plsda

STEPS = ("qc", "differential", "enrichment")
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


def _blocks(project: Project, experiment: Experiment, spec) -> tuple[list[_Block], list[Path], list[str]]:
    rows = experiment.samples[experiment.samples["omics"] == Omics.METABOLOMICS.value].copy()
    rows["column"] = rows["data_column"].where(rows["data_column"] != "", rows["sample_id"]) \
        if "data_column" in rows.columns else rows["sample_id"]
    rows["type"] = rows["sample_type"].replace("", "sample") if "sample_type" in rows.columns else "sample"
    rows = rows[~rows["type"].isin(_EXCLUDED_TYPES)]
    declared = set(rows["column"])

    tissue_files = getattr(spec, "tissue_files", None) or {}
    jobs = ([(str(t), project.resolve(p)) for t, p in tissue_files.items()] if tissue_files
            else [(None, project.resolve(f)) for f in spec.files])
    blocks: dict[str, _Block] = {}
    notes: list[str] = []
    for tissue, path in jobs:
        table = _read_table(path)
        table.columns = [str(c) for c in table.columns]
        candidates = rows if tissue is None else rows[rows["tissue"] == tissue]
        present = candidates[candidates["column"].isin(table.columns)]
        if tissue is not None and len(present) < len(candidates):
            absent = candidates.loc[~candidates["column"].isin(table.columns), "sample_id"]
            notes.append(f"代谢组 {tissue} 的数据表中找不到 {len(absent)} 个样本：{'、'.join(absent[:5])}")
        metadata = set(classify_data_columns(table.columns).metadata)
        annotations = [c for c in table.columns if c not in declared and c in metadata]
        ignored = [c for c in table.columns if c not in declared and c not in metadata]
        if ignored:
            notes.append(f"{path.name} 中有 {len(ignored)} 列像样本但不在 samples.tsv 中，未参与分析：{'、'.join(ignored[:5])}")
        for block_tissue, part in present.groupby("tissue", sort=False):
            if part["column"].duplicated().any():
                raise RunError(f"代谢组 {block_tissue} 有多个样本对应数据表中的同一列："
                               f"{'、'.join(part.loc[part['column'].duplicated(), 'column'][:5])}")
            if block_tissue in blocks:
                raise RunError(f"代谢组 {block_tissue} 的样本分布在多张数据表中（{blocks[block_tissue].path.name}、{path.name}）")
            by_column = part.set_index("column")
            order = [c for c in table.columns if c in by_column.index]
            biological = [by_column.at[c, "sample_id"] for c in order if by_column.at[c, "type"] not in _QC_TYPES]
            qc = [by_column.at[c, "sample_id"] for c in order if by_column.at[c, "type"] in _QC_TYPES]
            renamed = table[annotations + order].rename(columns=dict(zip(order, by_column.loc[order, "sample_id"])))
            blocks[block_tissue] = _Block(block_tissue, path, renamed, biological, qc, annotations)
    unread = sorted(set(rows["tissue"]) - set(blocks))
    if unread:
        notes.append(f"samples.tsv 中以下组织的代谢组样本没有对应的数据表：{'、'.join(unread)}")
    return list(blocks.values()), [path for _, path in jobs], notes


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


def run(ctx: RunContext) -> OmicsRun:
    experiment = ctx.experiment
    spec = experiment.config.omics.metabolomics.input
    out = ctx.omics_dir(Omics.METABOLOMICS)
    record = OmicsRun(Omics.METABOLOMICS)
    settings, record.notes = _settings(experiment.config.omics.metabolomics.analysis)
    blocks, record.inputs, notes = _blocks(ctx.project, experiment, spec)
    record.notes += notes
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
    differential, enrichment, contrasts_run = [], [], []

    for block in blocks:
        tissue = block.tissue
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
        meta = DatasetMeta(omics=Omics.METABOLOMICS, tissue=tissue, value_type=ValueType.PEAK_AREA, scale="log2",
                           normalization=pre.summary["normalization_method"], feature_id_type=feature_id,
                           source={"file": block.path.name, "detection_rule": pre.summary["detection_rule"]})
        write_dataset(out / "data" / tissue, Dataset(matrix, _features(processed, block.annotations, ids), meta))

        if ctx.wants("qc"):
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
            if ctx.wants("enrichment") and "kegg_map" in result.combined.columns:
                ora = calculate_pathway_enrichment(
                    result.combined, significance_col=f"Significant_{contrast.id}", feature_col=feature_id,
                    min_set_size=int(settings["ora_min_set_size"]), max_set_size=int(settings["ora_max_set_size"]),
                    pathway_scope=str(settings["pathway_scope"]),
                    exclude_overview=bool(settings["exclude_overview_pathways"]), fdr_method=fdr_method)
                if not ora.empty:
                    enrichment.append(ora_to_contract(ora, contrast_id=contrast.id, contrast_role=contrast.role.value,
                                                      tissue=tissue, query="all", gene_set_library="KEGG",
                                                      padj_method=normalize_fdr_method(fdr_method)))
        pd.DataFrame(validation).to_csv(details / "plsda_validation.csv", index=False, encoding="utf-8-sig")

    if differential:
        write_differential(pd.concat(differential, ignore_index=True), out / "tables" / "differential.parquet")
    if enrichment:
        write_enrichment(pd.concat(enrichment, ignore_index=True), out / "tables" / "enrichment.parquet")
    elif ctx.wants("enrichment") and differential:
        record.notes.append("代谢组数据表没有 kegg_map 列或没有可检验的 KEGG 通路，未做 ORA")
    record.parameters = {
        "tissues": [b.tissue for b in blocks],
        "contrasts": list(dict.fromkeys(contrasts_run)),
        "analysis": settings,
    }
    return record
