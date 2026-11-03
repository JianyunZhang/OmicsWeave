#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
蛋白组的统一运行器（``omicsweave run --omics proteomics``）：读取 ``experiment.yaml`` 声明的 DIA-NN 报告
或蛋白 × 样本矩阵，按组织预处理，再按其中的比较做差异蛋白分析，写出统一结果表。

计算沿用旧流程：报告按 q 值过滤、聚合为蛋白 × 样本矩阵后（矩阵输入直接读取），按组织用全部样本做 log2 与样本中位数居中
（不插补，按总观测数过滤）；每个比较只取该比较的样本，用 limma 方差压缩的线性模型检验
（``analysis.dp_method: wilcoxon`` 时用旧版秩和检验），协变量取 ``design.covariates`` 与
``analysis.design_covariates``。样本按矩阵中的列顺序参与计算，所以同一份报告与旧流程
（``omicsweave proteomics run``）的结果逐值一致。
富集分析（``enrichment`` 步骤）：预排序 GSEA 与旧流程相同（按调节 t 排序，每个基因集库分别检验，默认 KEGG_2019_Mouse）；另做按方向的 ORA（本地计算，背景为有 p 值的蛋白，见 :mod:`omicsweave.adapters.enrichment`）。原始 ``.raw`` 的 DIA-NN 搜库、图件与 ssGSEA 暂未接入。

输入（``omics.proteomics.input``，``format: diann_report``）::

    tissue_files: {liver: diann/liver/report.tsv}   # 每个组织一份报告（旧版布局）；没有时读取 files，
    files: [diann/report.tsv]                       #   一份报告中可以有多个组织的样本，组织由 samples.tsv 给出

样本名由报告 ``Run`` 列按旧流程的规则提取（如 ``..._Liver_86`` → ``Liver86``），即 samples.tsv 的 ``sample_id``；
不同时在 ``data_column`` 中写出提取后的名称。

已在其他软件中汇总好的蛋白 × 样本矩阵（``format: protein_matrix``，``.xlsx`` / ``.csv`` / ``.tsv``）::

    files: [data/protein_matrix.xlsx]   # 同样可用 tissue_files
    sheet: expression_data              # Excel 工作表，默认第一个
    feature_column: Gene.Name           # 特征名列，默认第一列；为空的行改用 uniprot_column
    uniprot_column: Protein.Id          # 可选：UniProt 列，用于富集与整合时的蛋白匹配
    scale: log2                         # 矩阵数值为 log2（默认）或 linear

其余数值列为样本，列名即 ``sample_id``，不同时写在 ``data_column``。不做 q 值过滤；同名特征按
``analysis.protein_aggregation`` 合并，之后的预处理与检验与 DIA-NN 报告相同（已归一化的矩阵会再做一次中位数居中）。

输出（``omics/proteomics/``）::

    data/<组织>/                   log2、中位数居中、未插补的强度（分析层数据集；行 = 蛋白，列 = 样本）
    qc/<组织>/                     样本 QC 指标、蛋白缺失率、重复相关性、PCA 得分与 QC 图
    tables/differential.parquet    统一差异结果表
    tables/enrichment.parquet      ORA 与预排序 GSEA
    details/<组织>/<比较>/          旧版格式的差异蛋白表与设计报告
"""
from __future__ import annotations

from dataclasses import fields
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd

from omicsweave.adapters.enrichment import EnrichmentSettings, ora_rows, prerank_rows
from omicsweave.contract.dataset import Dataset, DatasetMeta, ValueType, write_dataset
from omicsweave.contract.tables import write_differential, write_enrichment
from omicsweave.core.config import Experiment, Project
from omicsweave.core.genesets import GeneSetError, load_libraries
from omicsweave.core.model import Omics
from omicsweave.core.run import OmicsRun, RunContext, RunError, contrast_groups, sample_labels
from omicsweave.core.validate import biological_samples

from .differential_analysis import run_differential_analysis
from .enrichment.gsea import DEFAULT_GENE_SETS, build_ranked_list
from .experiment_config import normalize_fdr_method
from .export import dp_to_contract
from .proteomics_qc import generate_qc_outputs
from .run_proteomics_pipeline import load_and_clean_data, preprocess_data

STEPS = ("qc", "differential", "enrichment")
# 与旧版配置加载器的默认值相同
_DEFAULTS: dict[str, Any] = {
    "pval_cutoff": 0.05, "log2fc_cutoff": 0.58, "fdr_method": "fdr_bh", "dp_method": "limma",
    "min_total_observations": 2, "min_observations_per_group": 2, "imputation_seed": 42,
    "protein_aggregation": "log2_mean", "design_covariates": [],
    "gene_sets": [DEFAULT_GENE_SETS["kegg"]], "run_prerank_gsea": True, "gsea_ranking_metric": "moderated_t",
}
_AGGREGATIONS = ("log2_mean", "linear_mean")
_MATRIX_SCALES = ("log2", "linear")
# 输入格式 → 提示中的名称
_SOURCE_NAMES = {"diann_report": "DIA-NN 报告", "protein_matrix": "蛋白矩阵"}


def _preview(values, limit: int = 5) -> str:
    values = list(map(str, values))
    return "、".join(values[:limit]) + (f" 等 {len(values)} 个" if len(values) > limit else "")


def _samples(experiment: Experiment) -> pd.DataFrame:
    rows = biological_samples(experiment)
    rows = rows[rows["omics"] == Omics.PROTEOMICS.value].drop_duplicates("sample_id").copy()
    rows["column"] = rows["data_column"].where(rows["data_column"] != "", rows["sample_id"]) \
        if "data_column" in rows.columns else rows["sample_id"]
    return rows


def _read_table(path: Path, sheet: str | None) -> pd.DataFrame:
    if path.suffix.lower() in (".xlsx", ".xls"):
        with pd.ExcelFile(path) as book:
            if sheet is not None and sheet not in book.sheet_names:
                raise RunError(f"{path.name} 中没有工作表 {sheet!r}（有：{'、'.join(book.sheet_names)}）")
            table = pd.read_excel(book, sheet_name=sheet if sheet is not None else 0)
    else:
        table = pd.read_csv(path, sep="," if path.suffix.lower() == ".csv" else "\t")
    table.columns = [str(c).strip() for c in table.columns]
    return table


def _aggregate(values: pd.DataFrame, features: pd.Series, aggregation: str) -> pd.DataFrame:
    """同名特征合并，规则与 DIA-NN 报告聚合时相同：log2_mean 为 log2 空间求均值（几何均值）。"""
    values = values.where(values > 0)
    if aggregation == "log2_mean":
        return np.exp2(np.log2(values).groupby(features.to_numpy(), sort=False).mean())
    return values.groupby(features.to_numpy(), sort=False).mean()


def read_protein_matrix(path: Path, spec, aggregation: str,
                        notes: list[str]) -> tuple[pd.DataFrame, dict[str, str]]:
    """``protein_matrix`` 输入 → (蛋白 × 数据列名的线性强度, 蛋白 → UniProt)，与 DIA-NN 报告聚合后的形状相同。"""
    scale = str(getattr(spec, "scale", None) or "log2")
    if scale not in _MATRIX_SCALES:
        raise RunError(f"protein_matrix 的 scale 应为 {' / '.join(_MATRIX_SCALES)}，而不是 {scale!r}")
    table = _read_table(path, getattr(spec, "sheet", None))
    feature_column = getattr(spec, "feature_column", None) or table.columns[0]
    uniprot_column = getattr(spec, "uniprot_column", None)
    for label, column in (("feature_column", feature_column), ("uniprot_column", uniprot_column)):
        if column is not None and column not in table.columns:
            raise RunError(f"{path.name} 中没有 {label} 指定的列 {column!r}")

    def text(column: str | None) -> pd.Series:
        if column is None:
            return pd.Series("", index=table.index)
        return table[column].astype("string").fillna("").str.strip()

    uniprot = text(uniprot_column).str.split(";").str[0]
    features = text(feature_column).where(lambda s: s != "", uniprot)
    unnamed = features == ""
    if unnamed.any():
        notes.append(f"{path.name} 中有 {int(unnamed.sum())} 行没有特征名（{feature_column} 与 UniProt 均为空），已跳过")
    keep = ~unnamed
    samples = [c for c in table.columns
               if c not in (feature_column, uniprot_column) and pd.api.types.is_numeric_dtype(table[c])]
    values = table.loc[keep, samples].astype(float)
    if scale == "log2":
        values = np.exp2(values)
    features = features[keep]
    if features.duplicated().any():
        notes.append(f"{path.name} 中有 {int(features.duplicated().sum())} 行特征名重复，按 {aggregation} 合并")
    matrix = _aggregate(values, features, aggregation)
    matrix.index = matrix.index.astype(str)
    to_uniprot: dict[str, str] = {}
    for feature, accession in zip(features, uniprot[keep]):
        if accession:
            to_uniprot.setdefault(feature, accession)   # 与 DIA-NN 报告相同：同名特征取第一行的 UniProt
    return matrix, to_uniprot


def _read_report(path: Path, spec, aggregation: str, notes: list[str]) -> tuple[pd.DataFrame | None, dict[str, str]]:
    return load_and_clean_data(str(path), protein_aggregation=aggregation)


_LOADERS: dict[str, Callable[..., tuple[pd.DataFrame | None, dict[str, str]]]] = {
    "diann_report": _read_report, "protein_matrix": read_protein_matrix,
}


def _matrices(project: Project, spec, rows: pd.DataFrame, aggregation: str,
              record: OmicsRun) -> list[tuple[str, pd.DataFrame, dict[str, str]]]:
    """按组织返回 (组织, 蛋白 × sample_id 的原始强度, 蛋白 → UniProt)。列顺序与旧流程的矩阵相同。"""
    tissue_files = getattr(spec, "tissue_files", None) or {}
    jobs = ([(str(t), project.resolve(p)) for t, p in tissue_files.items()] if tissue_files
            else [(None, project.resolve(f)) for f in spec.files])
    source = _SOURCE_NAMES[spec.format]
    blocks: dict[str, tuple[pd.DataFrame, dict[str, str]]] = {}
    for tissue, path in jobs:
        if not path.is_file():
            raise RunError(f"找不到{source} {path}")
        record.inputs.append(path)
        matrix, to_uniprot = _LOADERS[spec.format](path, spec, aggregation, record.notes)
        if matrix is None or matrix.empty:
            record.notes.append(f"{path.name} 中没有通过过滤的蛋白定量，已跳过")
            continue
        candidates = rows if tissue is None else rows[rows["tissue"] == tissue]
        declared = dict(zip(candidates["column"], candidates["sample_id"]))
        undeclared = [c for c in matrix.columns if c not in declared]
        if undeclared:
            record.notes.append(f"{path.name} 中有 {len(undeclared)} 个样本不在 samples.tsv 中，未参与分析："
                                f"{_preview(undeclared)}")
        present = [c for c in matrix.columns if c in declared]
        by_tissue = candidates.set_index("column")["tissue"]
        for block_tissue in dict.fromkeys(by_tissue[present]):
            if block_tissue in blocks:
                raise RunError(f"蛋白组 {block_tissue} 的样本分布在多份报告中")
            columns = [c for c in present if by_tissue[c] == block_tissue]
            blocks[block_tissue] = (matrix[columns].rename(columns=declared), to_uniprot)
    found = {s for m, _ in blocks.values() for s in m.columns}
    missing = [s for s in rows["sample_id"] if s not in found]
    if missing:
        record.notes.append(f"{source}中找不到 {len(missing)} 个 samples.tsv 中的样本：{_preview(missing)}")
    return [(t, m, u) for t, (m, u) in blocks.items()]


def _covariate_table(experiment: Experiment, rows: pd.DataFrame, covariates: list[str]) -> pd.DataFrame:
    """协变量：声明为因素时按其类型；未声明的列全部能解析为数值时按连续变量，否则按分类变量。"""
    factors = experiment.config.design.factors
    table = rows.set_index("sample_id")[covariates].copy()
    for column in covariates:
        numeric = pd.to_numeric(table[column], errors="coerce")
        declared = factors.get(column)
        if (declared is not None and declared.type == "numeric") or (declared is None and numeric.notna().all()):
            table[column] = numeric
    return table


def run(ctx: RunContext) -> OmicsRun:
    experiment = ctx.experiment
    settings_block = experiment.config.omics.proteomics
    spec = settings_block.input
    if spec.level == "raw":
        raise RunError("原始 .raw 的 DIA-NN 搜库尚未接入统一运行；已有 DIA-NN 报告时把 omics.proteomics.input 改为 "
                       "{level: table, format: diann_report, files: [报告]}，或暂用 omicsweave proteomics run")
    if spec.format not in _LOADERS:
        raise RunError(f"蛋白组输入格式 {spec.format} 尚未接入统一运行（可用 {'、'.join(_LOADERS)}）")
    settings = {**_DEFAULTS, **settings_block.analysis}
    if settings["protein_aggregation"] not in _AGGREGATIONS:
        raise RunError(f"protein_aggregation 应为 {' / '.join(_AGGREGATIONS)}")
    fdr_method = normalize_fdr_method(str(settings["fdr_method"]))
    pval_cutoff, log2fc_cutoff = float(settings["pval_cutoff"]), float(settings["log2fc_cutoff"])
    rule = f"padj<{pval_cutoff} and |log2FC|>{log2fc_cutoff}"
    out = ctx.omics_dir(Omics.PROTEOMICS)
    record = OmicsRun(Omics.PROTEOMICS)
    rows = _samples(experiment)
    covariates = [c for c in dict.fromkeys([*experiment.config.design.covariates, *settings["design_covariates"]])
                  if c in rows.columns]
    covariate_table = _covariate_table(experiment, rows, covariates)
    blocks = _matrices(ctx.project, spec, rows, str(settings["protein_aggregation"]), record)
    if not blocks:
        raise RunError(f"{_SOURCE_NAMES[spec.format]}中没有找到 samples.tsv 声明的样本")

    tables, contrasts_run, enrichment_inputs = [], [], []
    for tissue, raw, to_uniprot in blocks:
        normalized = preprocess_data(raw, min_total_observations=int(settings["min_total_observations"]))
        if normalized.empty:
            record.notes.append(f"蛋白组 {tissue} 的蛋白都没有达到最少观测数，已跳过")
            continue
        features = pd.DataFrame({"uniprot_id": [to_uniprot.get(p, "") for p in normalized.index]},
                                index=normalized.index)
        meta = DatasetMeta(omics=Omics.PROTEOMICS, tissue=tissue, value_type=ValueType.INTENSITY, scale="log2",
                           normalization="median_centering", feature_id_type="gene_or_protein_group",
                           source={"input_format": spec.format,
                                   "protein_aggregation": settings["protein_aggregation"],
                                   "differential_input": "log2_normalized_unimputed"})
        write_dataset(out / "data" / tissue, Dataset(normalized, features, meta))

        if ctx.wants("qc"):
            labels = sample_labels(experiment, list(normalized.columns))
            qc_meta = covariate_table.loc[normalized.columns].assign(
                condition=[labels[s] for s in normalized.columns])
            generate_qc_outputs(raw.loc[:, normalized.columns], normalized, qc_meta, out / "qc" / tissue,
                                seed=int(settings["imputation_seed"]))
        if not (ctx.wants("differential") or ctx.wants("enrichment")):
            continue
        for contrast in experiment.config.design.contrasts:
            assigned = contrast_groups(experiment, contrast, Omics.PROTEOMICS)
            groups = {s: assigned[s] for s in normalized.columns if s in assigned}
            levels = (contrast.numerator, contrast.denominator)
            if not all(level in groups.values() for level in levels):
                if groups:
                    record.notes.append(f"比较 {contrast.id} 在蛋白组 {tissue} 中有一组没有样本，已跳过")
                continue
            samples = list(groups)
            details = out / "details" / tissue / contrast.id
            metadata = covariate_table.loc[samples].assign(condition=[groups[s] for s in samples])
            combined = run_differential_analysis(
                normalized[samples], groups, [(*levels, contrast.id)], str(details), tissue,
                pval_cutoff, log2fc_cutoff, fdr_method, str(settings["dp_method"]),
                feature_to_uniprot=to_uniprot, sample_metadata=metadata, design_covariates=covariates,
                min_observations_per_group=int(settings["min_observations_per_group"]), flush=lambda _m: None)
            result = pd.read_csv(details / f"dp_{tissue}_{contrast.id}.csv", dtype={"Protein": str, "uniprot_id": str})
            means = {level: normalized[[s for s in samples if groups[s] == level]].mean(axis=1, skipna=True)
                     for level in levels}
            contract = dp_to_contract(result, contrast_id=contrast.id, contrast_role=contrast.role.value,
                                      tissue=tissue, mean_numerator=means[contrast.numerator],
                                      mean_denominator=means[contrast.denominator], selection_rule=rule)
            tables.append(contract)
            enrichment_inputs.append((contrast, tissue, combined, contract))
            contrasts_run.append(contrast.id)

    if tables:
        write_differential(pd.concat(tables, ignore_index=True), out / "tables" / "differential.parquet")
    if ctx.wants("enrichment") and enrichment_inputs:
        _enrichment(ctx, settings, enrichment_inputs, record)
    record.parameters = {**record.parameters, "tissues": [t for t, _, _ in blocks], "covariates": covariates,
                         "contrasts": list(dict.fromkeys(contrasts_run)), "analysis": settings}
    return record


def _enrichment(ctx: RunContext, settings: dict[str, Any], inputs: list, record: OmicsRun) -> None:
    # 与旧流程（enrichment/gsea.py）相同：每个库分别检验；排序蛋白不少于基因集最小大小
    min_size = int(settings.get("gsea_min_set_size", 10))
    defaults = {"gsea_min_ranked": min_size, "gsea_min_set_size": min_size,
                "gsea_seed": int(settings.get("imputation_seed", 42)), "gsea_threads": 1,
                "gsea_joint_libraries": False, "fdr_method": "fdr_bh"}
    enrichment = EnrichmentSettings(**{**defaults, **{f.name: settings[f.name] for f in fields(EnrichmentSettings)
                                                      if f.name in settings}})
    try:
        libraries, files = load_libraries(ctx.project, [str(s) for s in settings["gene_sets"]])
    except GeneSetError as exc:
        record.notes.append(f"蛋白组富集分析未执行：{exc}")
        return
    record.inputs += files
    if settings.get("run_ssgsea"):
        record.notes.append("ssGSEA（样本层面的通路得分）暂未接入统一运行，可用 omicsweave proteomics run")
    parts = []
    for contrast, tissue, combined, contract in inputs:
        keys = dict(contrast_id=contrast.id, contrast_role=contrast.role.value, tissue=tissue,
                    omics=Omics.PROTEOMICS.value)
        parts.append(ora_rows(contract, libraries, enrichment, **keys))
        if settings["run_prerank_gsea"]:
            ranked = build_ranked_list(combined, contrast.id, metric=str(settings["gsea_ranking_metric"]))
            parts.append(prerank_rows(ranked.set_index("Protein")["ranking_score"], libraries, enrichment, **keys))
    parts = [p for p in parts if not p.empty]
    if parts:
        write_enrichment(pd.concat(parts, ignore_index=True),
                         ctx.omics_dir(Omics.PROTEOMICS) / "tables" / "enrichment.parquet")
    else:
        record.notes.append("蛋白组富集分析没有可报告的结果（入选或排序蛋白过少，或基因集与蛋白名不匹配）")
    record.parameters["enrichment"] = {"gene_sets": list(libraries), **enrichment.__dict__}
