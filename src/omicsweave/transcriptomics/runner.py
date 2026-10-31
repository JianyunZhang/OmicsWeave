#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
转录组的统一运行器（``omicsweave run --omics transcriptomics``）：读取 ``experiment.yaml`` 声明的
基因计数矩阵或 kallisto 定量结果，按组织、按比较做差异表达分析，写出统一结果表。

差异分析沿用旧流程：每组至少 2 个生物学重复时用 PyDESeq2（log2FC 经 apeGLM 收缩），设计为
``~[协变量 +] condition``，协变量取 ``design.covariates`` 中在该比较的样本间取值不同的列
（旧流程的 ``~batch + condition`` 即其特例）；样本量不足或拟合失败时退回秩和检验，有协变量时不退回，
不做未校正的检验。样本按计数矩阵的列顺序参与计算，所以同一份计数与旧流程
（``omicsweave transcriptomics run``）的差异结果逐值一致。
原始层 FASTQ 的质控、修剪与 kallisto 定量在阶段 3 接入；富集分析与 GSEA 暂未接入。

输入（``omics.transcriptomics.input``）::

    {level: table, format: counts_matrix, files: [gene_counts.tsv]}   # 基因 × 样本的计数，第一列为基因
    {level: table, format: kallisto, files: [quant/]}                  # quant/<样本>/abundance.tsv
        # 转录本 → 基因：input.tx2gene（transcript_id、gene_id、gene_symbol 三列的表），
        # 或 adapter.transcriptome_fasta（Ensembl cDNA FASTA 的表头）；
        # analysis.gene_aggregation：lengthScaledTPM（默认）或 simple_sum

数据中的样本名即 samples.tsv 的 ``sample_id``，不同时在 ``data_column`` 列写出数据中的名称。

输出（``omics/transcriptomics/``）::

    data/<组织>/                  基因水平计数（分析层数据集；行 = 基因，列 = 样本）
    qc/library_sizes.parquet      每个样本的总计数与检出基因数
    qc/quant_summary.parquet      kallisto 的比对率（kallisto 输入时）
    tables/differential.parquet   统一差异结果表
    details/<组织>/                原工具的差异结果表（含 baseMean、lfcSE、未收缩的 log2FC 等）
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from omicsweave.contract.dataset import Dataset, DatasetMeta, ValueType, empty_features, write_dataset
from omicsweave.contract.tables import write_differential
from omicsweave.core.config import Experiment, Project
from omicsweave.core.model import Omics
from omicsweave.core.run import OmicsRun, RunContext, RunError, contrast_groups
from omicsweave.core.validate import biological_samples

from . import qc_gate
from .export import de_to_contract
from .normalization import AGGREGATION_LENGTH_SCALED_TPM, VALID_AGGREGATION_METHODS
from .rnaseq_pipeline import (
    aggregate_counts_to_gene,
    build_tx2gene_from_fasta,
    merge_kallisto_counts,
    run_pydeseq2_analysis,
)
from .run_wilcoxon_de import run_wilcoxon_analysis

STEPS = ("qc", "differential")


def _preview(values, limit: int = 5) -> str:
    values = list(map(str, values))
    return "、".join(values[:limit]) + (f" 等 {len(values)} 个" if len(values) > limit else "")


def _samples(experiment: Experiment) -> pd.DataFrame:
    """转录组的生物学样本：sample_id、数据中的名称（column）与组织。"""
    rows = biological_samples(experiment)
    rows = rows[rows["omics"] == Omics.TRANSCRIPTOMICS.value].drop_duplicates("sample_id").copy()
    rows["column"] = rows["data_column"].where(rows["data_column"] != "", rows["sample_id"]) \
        if "data_column" in rows.columns else rows["sample_id"]
    return rows


def _read_matrix(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise RunError(f"找不到计数矩阵 {path}")
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    table = pd.read_csv(path, sep=sep, index_col=0)
    table.columns = [str(c) for c in table.columns]
    table.index = table.index.astype(str)
    return table


def _tx2gene(project: Project, spec, adapter: dict) -> tuple[dict, Path]:
    table_path = getattr(spec, "tx2gene", None)
    if table_path:
        path = project.resolve(table_path)
        table = pd.read_csv(path, sep="," if path.suffix.lower() == ".csv" else "\t", dtype=str)
        missing = [c for c in ("transcript_id", "gene_id") if c not in table.columns]
        if missing:
            raise RunError(f"{path.name} 缺少列：{'、'.join(missing)}")
        symbols = table["gene_symbol"] if "gene_symbol" in table.columns else table["gene_id"]
        return ({t: {"gene_id": g, "gene_symbol": s or g}
                 for t, g, s in zip(table["transcript_id"], table["gene_id"], symbols.fillna(""))}, path)
    fasta = adapter.get("transcriptome_fasta")
    if not fasta:
        raise RunError("kallisto 输入需要转录本 → 基因的对应关系：在 input.tx2gene 中给出对应表，"
                       "或在 adapter.transcriptome_fasta 中给出 Ensembl cDNA FASTA")
    path = project.resolve(fasta)
    if not path.is_file():
        raise RunError(f"找不到转录组 FASTA {path}")
    return build_tx2gene_from_fasta(str(path)), path


def _load_counts(ctx: RunContext, rows: pd.DataFrame, record: OmicsRun) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    """返回 (基因 × sample_id 的计数, kallisto 比对率或 None)。"""
    settings = ctx.experiment.config.omics.transcriptomics
    spec = settings.input
    if spec.level == "raw":
        raise RunError("原始层 FASTQ 的质控、修剪与定量在阶段 3 接入统一运行。已有 kallisto 定量结果时，把 "
                       "omics.transcriptomics.input 改为 {level: table, format: kallisto, files: [定量目录]}；"
                       "已有基因计数矩阵时用 format: counts_matrix；也可暂用 omicsweave transcriptomics run")
    if spec.format not in ("counts_matrix", "kallisto"):
        raise RunError(f"转录组输入格式 {spec.format} 尚未接入统一运行（可用 counts_matrix、kallisto）")
    files = [ctx.project.resolve(f) for f in spec.files]
    columns = dict(zip(rows["column"], rows["sample_id"]))

    if spec.format == "counts_matrix":
        if len(files) != 1:
            raise RunError("counts_matrix 目前只读取一个文件，请把各组织的样本放在同一张计数矩阵中")
        table = _read_matrix(files[0])
        present = [c for c in table.columns if c in columns]
        record.inputs.append(files[0])
        counts = table[present].apply(pd.to_numeric, errors="coerce")
        if counts.isna().any().any():
            raise RunError(f"{files[0].name} 的样本列中有非数值或空值")
        quant = None
    else:
        quant_dir = files[0]
        present = [c for c in columns if (quant_dir / c / "abundance.tsv").is_file()]
        record.inputs += [quant_dir / c / "abundance.tsv" for c in present]
        tx_counts, tpm, efflen = merge_kallisto_counts(str(quant_dir), present)
        if tx_counts is None:
            raise RunError(f"{quant_dir} 下没有找到任何样本的 abundance.tsv")
        tx2gene, mapping = _tx2gene(ctx.project, spec, settings.adapter)
        record.inputs.append(mapping)
        method = str(settings.analysis.get("gene_aggregation", AGGREGATION_LENGTH_SCALED_TPM))
        if method not in VALID_AGGREGATION_METHODS:
            raise RunError(f"gene_aggregation 应为 {' / '.join(VALID_AGGREGATION_METHODS)}")
        counts, info = aggregate_counts_to_gene(tx_counts, tx2gene, tpm_df=tpm, efflen_df=efflen, method=method)
        record.parameters["gene_aggregation"] = info.get("gene_aggregation", method)
        quant = qc_gate.collect_quant_summary(quant_dir, present)
        quant.insert(0, "sample_id", quant.pop("sample").map(columns))
    missing = [s for c, s in columns.items() if c not in present]
    if missing:
        record.notes.append(f"转录组数据中找不到 {len(missing)} 个 samples.tsv 中的样本：{_preview(missing)}")
    if not present:
        raise RunError("转录组数据中没有找到 samples.tsv 声明的样本")
    return counts.rename(columns=columns), quant


def _differential(ctx: RunContext, counts: pd.DataFrame, rows: pd.DataFrame, record: OmicsRun) -> list[pd.DataFrame]:
    experiment = ctx.experiment
    design = experiment.config.design
    covariates = [c for c in design.covariates if c in rows.columns]
    by_sample = rows.set_index("sample_id")
    tables, designs = [], {}
    for tissue in dict.fromkeys(rows["tissue"]):
        # 按计数矩阵的列顺序取样本（PyDESeq2 与秩和检验的结果都与此顺序一致）
        tissue_samples = [s for s in counts.columns if by_sample.at[s, "tissue"] == tissue]
        if len(tissue_samples) < 2:
            continue
        details = ctx.omics_dir(Omics.TRANSCRIPTOMICS) / "details" / tissue
        details.mkdir(parents=True, exist_ok=True)
        for contrast in design.contrasts:
            groups = contrast_groups(experiment, contrast, Omics.TRANSCRIPTOMICS)
            samples = [s for s in tissue_samples if s in groups]
            n_num = sum(groups[s] == contrast.numerator for s in samples)
            n_den = len(samples) - n_num
            if n_num < 1 or n_den < 1:
                if samples:
                    record.notes.append(f"比较 {contrast.id} 在转录组 {tissue} 中有一组没有样本，已跳过")
                continue
            meta = pd.DataFrame({"condition": [groups[s] for s in samples]}, index=pd.Index(samples))
            for column in covariates:
                meta[column] = by_sample.loc[samples, column].astype(str)
            active = [c for c in covariates if meta[c].nunique() > 1]
            formula = " + ".join([*active, "condition"])
            label = f"{contrast.id}（{tissue}）"
            result = None
            if n_num >= 2 and n_den >= 2:
                try:
                    result = run_pydeseq2_analysis(counts[samples], meta, design_formula=formula,
                                                   contrast=("condition", contrast.numerator, contrast.denominator))
                except Exception as exc:  # 与旧流程相同：拟合失败时视情况退回秩和检验
                    record.notes.append(f"{label} 的 PyDESeq2 拟合失败：{exc}")
            if result is None:
                if active:
                    record.notes.append(f"{label} 含协变量 {'、'.join(active)}，不能退回未校正的秩和检验，已跳过")
                    continue
                result = run_wilcoxon_analysis(counts[samples], meta, contrast.numerator, contrast.denominator)
            if result is None:
                continue
            designs[f"{contrast.id}/{tissue}"] = f"~{formula}" if "pydeseq2" in set(result["test_method"]) \
                else str(result["test_method"].iloc[0])
            result.rename_axis("gene").to_csv(details / f"de_{contrast.id}.csv")
            tables.append(de_to_contract(result, contrast_id=contrast.id, contrast_role=contrast.role.value,
                                         tissue=tissue))
    record.parameters["designs"] = designs
    return tables


def run(ctx: RunContext) -> OmicsRun:
    out = ctx.omics_dir(Omics.TRANSCRIPTOMICS)
    record = OmicsRun(Omics.TRANSCRIPTOMICS)
    rows = _samples(ctx.experiment)
    counts, quant = _load_counts(ctx, rows, record)
    tissues = rows.set_index("sample_id")["tissue"]
    record.parameters["input_format"] = ctx.experiment.config.omics.transcriptomics.input.format

    id_type = getattr(ctx.experiment.config.omics.transcriptomics.input, "feature_id_type", None) or "gene_symbol"
    for tissue in dict.fromkeys(tissues[counts.columns]):
        samples = [s for s in counts.columns if tissues[s] == tissue]
        matrix = counts[samples]
        meta = DatasetMeta(omics=Omics.TRANSCRIPTOMICS, tissue=str(tissue), value_type=ValueType.COUNTS,
                           scale="linear", feature_id_type=id_type,
                           source={"input_format": record.parameters["input_format"],
                                   **({"gene_aggregation": record.parameters["gene_aggregation"]}
                                      if "gene_aggregation" in record.parameters else {})})
        write_dataset(out / "data" / str(tissue), Dataset(matrix, empty_features(matrix.index), meta))

    if ctx.wants("qc"):
        (out / "qc").mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"sample_id": counts.columns, "tissue": tissues[counts.columns].to_numpy(),
                      "total_counts": counts.sum().to_numpy(),
                      "detected_genes": (counts > 0).sum().to_numpy()}).to_parquet(
            out / "qc" / "library_sizes.parquet", index=False)
        if quant is not None:
            quant.to_parquet(out / "qc" / "quant_summary.parquet", index=False)

    if ctx.wants("differential"):
        tables = _differential(ctx, counts, rows, record)
        if tables:
            write_differential(pd.concat(tables, ignore_index=True), out / "tables" / "differential.parquet")
        record.parameters["contrasts"] = list(dict.fromkeys(t["contrast_id"].iloc[0] for t in tables))
    return record
