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
富集分析（``enrichment`` 步骤）默认对每个比较仅运行本地 ORA（背景为该比较中有 p 值的基因，不上传 Enrichr）；预排序 GSEA（``gseapy.prerank``，log2FC × −log10 p，多个库合并检验）由 ``analysis.run_gsea`` 控制（默认 ``false``，可用配置或命令行 ``--gsea`` 开启，见 :mod:`omicsweave.adapters.enrichment`）。基因集库由 ``analysis.gene_sets`` 指定（本地 GMT 或 Enrichr 库名，见 :mod:`omicsweave.core.genesets`）。
原始层 FASTQ 的质控、修剪与 kallisto 定量在阶段 3 接入。

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
    qc/sample_quality.parquet     通用样本质检结果（表达谱离群、线粒体占比、组织标记保留度、跨组织污染）
    qc/sample_qc_metrics.parquet  逐样本质检诊断指标汇总
    tables/differential.parquet   统一差异结果表
    tables/enrichment.parquet     ORA 与预排序 GSEA
    figures/                      火山图 volcano_<组织>_<比较>.pdf / .png（``figures`` 步骤，见
                                  :mod:`omicsweave.visualization.volcano`；纵轴为校正后 p）
    details/<组织>/                原工具的差异结果表（含 baseMean、lfcSE、未收缩的 log2FC 等）
"""
from __future__ import annotations

from dataclasses import fields
from pathlib import Path

import numpy as np
import pandas as pd

from omicsweave.adapters.enrichment import EnrichmentSettings, ora_rows, prerank_rows, resolve_run_gsea
from omicsweave.contract.dataset import Dataset, DatasetMeta, ValueType, empty_features, write_dataset
from omicsweave.contract.tables import empty_enrichment, write_differential, write_enrichment
from omicsweave.core.config import Experiment, Project
from omicsweave.core.genesets import GeneSetError, load_libraries
from omicsweave.core.model import Omics
from omicsweave.core.run import OmicsRun, RunContext, RunError, contrast_groups
from omicsweave.core.tissues import format_alias_note, select_tissue_file_rows
from omicsweave.core.validate import biological_samples
from omicsweave.visualization.volcano import Thresholds, write_bundle_volcanoes

from . import qc_gate
from .export import LOG2FC_CUTOFF, PADJ_CUTOFF, de_to_contract
from .normalization import AGGREGATION_LENGTH_SCALED_TPM, VALID_AGGREGATION_METHODS
from .rnaseq_pipeline import (
    aggregate_counts_to_gene,
    build_tx2gene_from_fasta,
    merge_kallisto_counts,
    run_pydeseq2_analysis,
)
from .run_wilcoxon_de import run_wilcoxon_analysis
from .sample_qc import SampleQCResult, evaluate_sample_quality, resolve_sample_qc_config

STEPS = ("qc", "differential", "enrichment", "figures")
# 与旧流程（enrichment/_run.py、gsea/_run.py）使用的 Enrichr 库相同
DEFAULT_GENE_SETS = [
    "GO_Biological_Process_2025", "GO_Molecular_Function_2025", "GO_Cellular_Component_2025",
    "KEGG_2019_Mouse", "WikiPathways_2024_Mouse", "Reactome_Pathways_2024",
]


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


def _load_counts(ctx: RunContext, rows: pd.DataFrame,
                 record: OmicsRun) -> tuple[pd.DataFrame, pd.DataFrame | None, dict[str, str]]:
    """返回 (基因 × sample_id 的计数, kallisto 比对率或 None, 生效的组织别名)。"""
    settings = ctx.experiment.config.omics.transcriptomics
    spec = settings.input
    if spec.level == "raw":
        raise RunError("原始层 FASTQ 的质控、修剪与定量在阶段 3 接入统一运行。已有 kallisto 定量结果时，把 "
                       "omics.transcriptomics.input 改为 {level: table, format: kallisto, files: [定量目录]}；"
                       "已有基因计数矩阵时用 format: counts_matrix；也可暂用 omicsweave transcriptomics run")
    if spec.format not in ("counts_matrix", "kallisto"):
        raise RunError(f"转录组输入格式 {spec.format} 尚未接入统一运行（可用 counts_matrix、kallisto）")
    applied_aliases = ctx.experiment.applied_sample_aliases(Omics.TRANSCRIPTOMICS)
    files = [ctx.project.resolve(f) for f in spec.files]

    if spec.format == "counts_matrix":
        tissue_files = getattr(spec, "tissue_files", None) or {}
        if tissue_files:
            jobs = [
                (str(t), ctx.experiment.config.canonical_tissue(str(t), Omics.TRANSCRIPTOMICS), ctx.project.resolve(p))
                for t, p in tissue_files.items()
            ]
            canonical_counts: dict[str, int] = {}
            for raw_key, canon_key, _ in jobs:
                canonical_counts[canon_key] = canonical_counts.get(canon_key, 0) + 1
                if raw_key != canon_key:
                    applied_aliases.setdefault(raw_key, canon_key)
        else:
            if not files:
                raise RunError("counts_matrix 需要至少一张包含样本的计数矩阵")
            jobs = [(None, None, path) for path in files]
            canonical_counts = {}
        parts: list[pd.DataFrame] = []
        seen_samples: set[str] = set()
        for raw_key, canon_key, path in jobs:
            table = _read_matrix(path)
            candidates = rows if canon_key is None else select_tissue_file_rows(
                rows, raw_key, canon_key,
                multi_key_canonical=canonical_counts[canon_key] > 1,
                raw_tissues=ctx.experiment.raw_tissues,
            )
            if candidates["column"].duplicated().any():
                dups = candidates.loc[candidates["column"].duplicated(), "column"].tolist()
                raise RunError(f"{path.name} 中有多个转录组样本对应同一数据列：{_preview(dups)}")
            col_map = dict(zip(candidates["column"], candidates["sample_id"]))
            present_cols = [c for c in table.columns if c in col_map]
            record.inputs.append(path)
            if not present_cols:
                continue
            part = table[present_cols].apply(pd.to_numeric, errors="coerce")
            if part.isna().any().any():
                raise RunError(f"{path.name} 的样本列中有非数值或空值")
            part = part.rename(columns=col_map)
            overlap = [s for s in part.columns if s in seen_samples]
            if overlap:
                raise RunError(f"转录组样本 {_preview(overlap)} 在多张计数矩阵中重复出现")
            seen_samples.update(part.columns)
            parts.append(part)
        if not parts:
            raise RunError("转录组数据中没有找到 samples.tsv 声明的样本")
        if len(parts) == 1:
            counts = parts[0]
        else:
            counts = pd.concat(parts, axis=1)
            if counts.isna().any().any():
                n_missing = int(counts.isna().any(axis=1).sum())
                counts = counts.fillna(0)
                record.notes.append(f"合并多张计数矩阵时有 {n_missing} 个基因未在全部文件中出现，缺失处按 0 填充")
        quant = None
        missing = [s for s in rows["sample_id"] if s not in seen_samples]
        if missing:
            record.notes.append(f"转录组数据中找不到 {len(missing)} 个 samples.tsv 中的样本：{_preview(missing)}")
        return counts, quant, applied_aliases

    columns = dict(zip(rows["column"], rows["sample_id"]))
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
    return counts.rename(columns=columns), quant, applied_aliases


def _flag_foreign_markers(
    contract: pd.DataFrame,
    tissue: str,
    samples: list[str],
    groups: dict[str, str],
    contrast,
    qc_result: SampleQCResult | None,
) -> pd.DataFrame:
    """对保留在 DE 中且含轻度跨组织污染标记样本的比较，给方向吻合的异源组织特异基因标上 ``foreign_tissue_marker``。"""
    if qc_result is None or qc_result.quality.empty or not qc_result.foreign_genes or contract.empty:
        return contract
    q = qc_result.quality
    hit = q[
        (q["metric"] == "cross_tissue_contamination")
        & (q["flag"].fillna("").astype(str).str.strip() != "")
        & (q["tissue"].astype(str) == str(tissue))
        & (q["sample_id"].astype(str).isin(set(samples)))
    ]
    if hit.empty:
        return contract
    pos_genes: set[str] = set()
    neg_genes: set[str] = set()
    for row in hit.itertuples(index=False):
        genes = qc_result.foreign_genes.get((str(tissue), str(row.foreign_tissue)), ())
        if not genes:
            continue
        side = groups.get(str(row.sample_id))
        if side == contrast.numerator:
            pos_genes.update(genes)
        elif side == contrast.denominator:
            neg_genes.update(genes)
    if not pos_genes and not neg_genes:
        return contract
    pos_lower = {g.lower() for g in pos_genes}
    neg_lower = {g.lower() for g in neg_genes}
    fids = contract["feature_id"].astype(str)
    fcs = pd.to_numeric(contract["log2fc"], errors="coerce")
    hit_pos = (fids.isin(pos_genes) | fids.str.lower().isin(pos_lower)) & fcs.gt(0)
    hit_neg = (fids.isin(neg_genes) | fids.str.lower().isin(neg_lower)) & fcs.lt(0)
    mask = (hit_pos | hit_neg).to_numpy()
    if not mask.any():
        return contract
    flags = contract["flags"].fillna("").astype(str) if "flags" in contract.columns else pd.Series("", index=contract.index)
    joined = [
        ";".join(p for p in (f, "foreign_tissue_marker") if p)
        if m and "foreign_tissue_marker" not in f.split(";") else f
        for f, m in zip(flags, mask)
    ]
    return contract.assign(flags=joined)


def _differential(ctx: RunContext, counts: pd.DataFrame, rows: pd.DataFrame,
                  record: OmicsRun, qc_result: SampleQCResult | None = None) -> list[tuple[object, str, pd.DataFrame, pd.DataFrame]]:
    """返回每个比较、每个组织的 (比较, 组织, 原工具结果, 统一差异结果行)。"""
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
            ctx.progress(Omics.TRANSCRIPTOMICS, "differential", tissue=str(tissue), contrast_id=contrast.id)
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
                except ImportError as exc:
                    name = getattr(exc, "name", None) or ""
                    if name.split(".")[0] == "omicsweave":
                        raise
                    pkg = name.split(".")[0] or "pydeseq2"
                    raise RunError(
                        f"{label} 需要 PyDESeq2（缺少或无法导入依赖 {pkg}：{exc}），不能静默退回秩和检验；"
                        f'请安装：pip install "omicsweave[transcriptomics]"（或 pip install {pkg}）'
                    ) from exc
                except Exception as exc:  # 与旧流程相同：仅数值拟合失败时视情况退回秩和检验
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
            result.rename_axis("gene").to_csv(details / f"de_{contrast.id}.csv", encoding="utf-8-sig")
            contract = de_to_contract(result, contrast_id=contrast.id, contrast_role=contrast.role.value,
                                      tissue=tissue)
            contract = _flag_foreign_markers(contract, tissue, samples, groups, contrast, qc_result)
            tables.append((contrast, tissue, result, contract))
    record.parameters["designs"] = designs
    return tables


def run(ctx: RunContext) -> OmicsRun:
    out = ctx.omics_dir(Omics.TRANSCRIPTOMICS)
    record = OmicsRun(Omics.TRANSCRIPTOMICS)
    rows = _samples(ctx.experiment)
    counts, quant, applied_aliases = _load_counts(ctx, rows, record)
    if applied_aliases:
        record.notes.append(format_alias_note(Omics.TRANSCRIPTOMICS.value, applied_aliases))
        record.parameters["tissue_aliases_applied"] = applied_aliases
    tissues = rows.set_index("sample_id")["tissue"]
    record.parameters["input_format"] = ctx.experiment.config.omics.transcriptomics.input.format

    # ── 通用样本质量控制（表达谱离群、线粒体异常、本组织标记保留度、跨组织 RNA 污染） ──
    raw_counts = counts
    qc_cfg = resolve_sample_qc_config(ctx.experiment.config.omics.transcriptomics.analysis.get("sample_qc"))
    manual_flags = (
        dict(zip(rows["sample_id"].astype(str), rows["qc_flag"].fillna("").astype(str)))
        if "qc_flag" in rows.columns else {}
    )
    factor_cols = [c for c in ctx.experiment.config.design.factors if c in rows.columns]
    sample_groups = (
        {
            str(r["sample_id"]): "/".join(str(r.get(c, "") or "") for c in factor_cols)
            for r in rows[["sample_id", *factor_cols]].to_dict("records")
        }
        if factor_cols else None
    )
    qc_result = evaluate_sample_quality(raw_counts, tissues, qc_cfg, manual_flags=manual_flags, groups=sample_groups)

    if ctx.wants("qc") or (qc_cfg.get("enabled", True) and (qc_result.excluded or qc_result.flagged)):
        (out / "qc").mkdir(parents=True, exist_ok=True)
    if ctx.wants("qc"):
        ctx.progress(Omics.TRANSCRIPTOMICS, "qc")
        pd.DataFrame({"sample_id": raw_counts.columns, "tissue": tissues[raw_counts.columns].to_numpy(),
                      "total_counts": raw_counts.sum().to_numpy(),
                      "detected_genes": (raw_counts > 0).sum().to_numpy()}).to_parquet(
            out / "qc" / "library_sizes.parquet", index=False)
        if quant is not None:
            quant.to_parquet(out / "qc" / "quant_summary.parquet", index=False)
        if qc_cfg.get("enabled", True):
            qc_result.metrics.to_parquet(out / "qc" / "sample_qc_metrics.parquet", index=False)
    if qc_cfg.get("enabled", True) and (ctx.wants("qc") or qc_result.excluded or qc_result.flagged):
        qc_result.quality.to_parquet(out / "qc" / "sample_quality.parquet", index=False)

    if qc_cfg.get("filter_excluded", True) and qc_result.excluded:
        kept = [c for c in raw_counts.columns if str(c) not in qc_result.excluded]
        counts = raw_counts[kept]
        record.notes.append(
            f"转录组通用质检过滤 {len(qc_result.excluded)} 个严重异常样本（{_preview(sorted(qc_result.excluded))}），"
            f"不参与差异与下游稳健性分析；明细见 omics/transcriptomics/qc/sample_quality.parquet"
        )
    if qc_result.flagged:
        record.notes.append(
            f"转录组通用质检标记 {len(qc_result.flagged)} 个可疑样本（{_preview(sorted(qc_result.flagged))}），"
            f"已写入 omics/transcriptomics/qc/sample_quality.parquet 供比较可靠性评估"
        )
    if qc_result.excluded or qc_result.flagged:
        record.parameters["sample_qc"] = {
            "enabled": bool(qc_cfg.get("enabled", True)),
            "filter_excluded": bool(qc_cfg.get("filter_excluded", True)),
            "excluded_samples": sorted(qc_result.excluded),
            "flagged_samples": sorted(qc_result.flagged),
        }

    id_type = getattr(ctx.experiment.config.omics.transcriptomics.input, "feature_id_type", None) or "gene_symbol"
    for tissue in dict.fromkeys(tissues[counts.columns]):
        ctx.progress(Omics.TRANSCRIPTOMICS, "dataset", tissue=str(tissue))
        samples = [s for s in counts.columns if tissues[s] == tissue]
        if not samples:
            continue
        matrix = counts[samples]
        meta = DatasetMeta(omics=Omics.TRANSCRIPTOMICS, tissue=str(tissue), value_type=ValueType.COUNTS,
                           scale="linear", feature_id_type=id_type,
                           source={"input_format": record.parameters["input_format"],
                                   **({"gene_aggregation": record.parameters["gene_aggregation"]}
                                      if "gene_aggregation" in record.parameters else {})})
        write_dataset(out / "data" / str(tissue), Dataset(matrix, empty_features(matrix.index), meta))

    if ctx.wants("differential") or ctx.wants("enrichment"):
        results = _differential(ctx, counts, rows, record, qc_result=qc_result) if counts.shape[1] >= 2 else []
        if results:
            write_differential(pd.concat([r[3] for r in results], ignore_index=True),
                               out / "tables" / "differential.parquet")
        record.parameters["contrasts"] = list(dict.fromkeys(r[0].id for r in results))
        if ctx.wants("enrichment") and results:
            _enrichment(ctx, results, record)
    if ctx.wants("figures"):
        ctx.progress(Omics.TRANSCRIPTOMICS, "figures")
        record.notes += write_bundle_volcanoes(out, ctx.experiment.config.design.contrasts, metric="padj",
                                               thresholds=Thresholds(p=PADJ_CUTOFF, log2fc=LOG2FC_CUTOFF))
    return record


def _ranking(result: pd.DataFrame) -> pd.Series:
    """与旧流程（gsea/_run.py）相同：log2FC × −log10(p) 降序，重复基因保留第一个，去掉 0。"""
    table = result[["log2FoldChange", "pvalue"]].dropna()
    metric = (table["log2FoldChange"] * -np.log10(table["pvalue"].clip(lower=1e-300))).sort_values(ascending=False)
    metric = metric[~metric.index.duplicated(keep="first")]
    return metric[metric != 0]


def _enrichment(ctx: RunContext, results: list, record: OmicsRun) -> None:
    analysis = ctx.experiment.config.omics.transcriptomics.analysis
    try:
        do_gsea = resolve_run_gsea(analysis, override=ctx.run_gsea)
    except ValueError as exc:
        raise RunError(str(exc)) from exc
    settings = EnrichmentSettings(
        **{f.name: analysis[f.name] for f in fields(EnrichmentSettings) if f.name in analysis and f.name != "run_gsea"},
        run_gsea=do_gsea,
    )
    try:
        libraries, files = load_libraries(ctx.project, list(analysis.get("gene_sets", DEFAULT_GENE_SETS)))
    except GeneSetError as exc:
        record.notes.append(f"转录组富集分析未执行：{exc}")
        return
    record.inputs += files
    parts = []
    for contrast, tissue, raw, contract in results:
        ctx.progress(Omics.TRANSCRIPTOMICS, "enrichment", tissue=str(tissue), contrast_id=contrast.id)
        keys = dict(contrast_id=contrast.id, contrast_role=contrast.role.value, tissue=tissue,
                    omics=Omics.TRANSCRIPTOMICS.value)
        parts.append(ora_rows(contract, libraries, settings, **keys))
        if settings.run_gsea:
            parts.append(prerank_rows(_ranking(raw), libraries, settings, **keys))
    parts = [p for p in parts if not p.empty]
    enr_path = ctx.omics_dir(Omics.TRANSCRIPTOMICS) / "tables" / "enrichment.parquet"
    if parts:
        write_enrichment(pd.concat(parts, ignore_index=True), enr_path)
    else:
        write_enrichment(empty_enrichment(), enr_path)
        record.notes.append("转录组富集分析没有可报告的结果（入选基因或排序基因过少，或基因集与数据的基因名不匹配）")
    record.parameters["enrichment"] = {"gene_sets": list(libraries), **settings.__dict__}
