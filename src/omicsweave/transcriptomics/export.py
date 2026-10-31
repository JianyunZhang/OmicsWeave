#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
转录组差异表达结果 → 统一差异结果表（``omicsweave.contract.tables.DIFFERENTIAL_COLUMNS``）。

输入为 :func:`omicsweave.transcriptomics.rnaseq_pipeline.run_pydeseq2_analysis`（PyDESeq2，log2FC 经 apeGLM
收缩）或 :func:`omicsweave.transcriptomics.run_wilcoxon_de.run_wilcoxon_analysis`（样本量不足以拟合 DESeq2 时的
秩和检验）的输出，行索引为基因。两者都以 ``padj<0.05`` 且 ``|log2FC|>1`` 作为入选规则。
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from omicsweave.contract.dataset import ContractError
from omicsweave.core.model import Omics

SELECTION_RULE = "padj<0.05 and |log2FC|>1"
# 原工具的检验方法名 → 统一表中的 method
_METHODS = {"pydeseq2": "deseq2", "wilcoxon": "wilcoxon", "wilcoxon_effect_only": "effect_size_only"}
# 原表中保留为附加列的列：原列名 → 统一表中的列名
_EXTRA = {
    "baseMean": "base_mean", "lfcSE": "lfc_se", "log2FoldChange_unshrunk": "log2fc_unshrunk",
    "lfc_shrinkage": "lfc_shrinkage", "design_formula": "design_formula", "batch_adjusted": "batch_adjusted",
    "n_batches": "n_batches", "batch_ids": "batch_ids", "normalization": "normalization",
}


def de_to_contract(results: pd.DataFrame, *, contrast_id: str, contrast_role: str, tissue: str) -> pd.DataFrame:
    """把一个比较、一个组织的差异表达结果转换为统一差异结果表。"""
    missing = [c for c in ("log2FoldChange", "pvalue", "padj", "significant", "test_method", "n_case",
                           "n_control", "evidence_level") if c not in results.columns]
    if missing:
        raise ContractError(f"差异表达结果缺少列：{', '.join(missing)}")
    tested = pd.to_numeric(results["pvalue"], errors="coerce").notna().to_numpy()
    p_floor = pd.to_numeric(results.get("min_achievable_pvalue", pd.Series(np.nan, index=results.index)),
                            errors="coerce")
    flags = np.where(results["evidence_level"] == "effect_size_only", "low_n",
                     np.where(p_floor > 0.05, "underpowered", ""))
    statistic = results["stat"] if "stat" in results.columns else pd.Series(np.nan, index=results.index)
    out = pd.DataFrame({
        "contrast_id": contrast_id,
        "contrast_role": contrast_role,
        "tissue": tissue,
        "omics": Omics.TRANSCRIPTOMICS.value,
        "feature_id": results.index.astype(str),
        "log2fc": pd.to_numeric(results["log2FoldChange"], errors="coerce"),
        "statistic": pd.to_numeric(statistic, errors="coerce"),
        "pvalue": pd.to_numeric(results["pvalue"], errors="coerce"),
        "padj": pd.to_numeric(results["padj"], errors="coerce"),
        # PyDESeq2 与秩和检验都用 Benjamini-Hochberg 校正
        "padj_method": np.where(tested, "fdr_bh", "none"),
        "n_numerator": pd.to_numeric(results["n_case"]).astype("int64"),
        "n_denominator": pd.to_numeric(results["n_control"]).astype("int64"),
        "min_achievable_p": p_floor,
        "method": results["test_method"].map(lambda m: _METHODS.get(str(m), str(m))),
        "inference_unit": "biological_replicate",
        "evidence_level": results["evidence_level"].astype(str),
        "selected": results["significant"].fillna(False).astype(bool),
        "selection_rule": SELECTION_RULE,
        "flags": flags,
    }, index=results.index)
    for source, target in _EXTRA.items():
        if source in results.columns:
            out[target] = results[source]
    return out.reset_index(drop=True)
