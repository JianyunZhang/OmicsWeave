#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
蛋白组差异结果 → 统一差异结果表（``omicsweave.contract.tables.DIFFERENTIAL_COLUMNS``）。

输入为 :func:`omicsweave.proteomics.differential_analysis.run_differential_analysis` 写出的每个比较的表
（``dp_<组织>_<比较>.csv``）：在 log2、按样本中位数居中、未插补的强度上做 limma 方差压缩的线性模型
（或旧版秩和检验）。``n_numerator`` / ``n_denominator`` 是该蛋白实际有观测值的样本数。
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from omicsweave.contract.dataset import ContractError
from omicsweave.core.model import Omics

LOG2_UNITS = "log2 median-centered intensity"
_METHODS = {"limma_moderated_linear_model": "limma", "ols_t_test_unmoderated": "ols_t_unmoderated"}


def _flag(status: str, tested: bool) -> str:
    """未检验的原因：样本或观测不足记为 low_n，设计不可辨识等记为原状态。"""
    if tested:
        return "unmoderated" if status == "ok_unmoderated_t" else ""
    return "low_n" if status.startswith("insufficient") else status


_EXTRA = {"uniprot_id": "uniprot_id", "standard_error": "lfc_se", "residual_df": "residual_df",
          "inference_status": "inference_status", "design_formula": "design_formula"}


def dp_to_contract(table: pd.DataFrame, *, contrast_id: str, contrast_role: str, tissue: str,
                   mean_numerator: pd.Series, mean_denominator: pd.Series, selection_rule: str) -> pd.DataFrame:
    """把一个比较、一个组织的差异蛋白表转换为统一差异结果表；均值按 ``Protein`` 对齐。"""
    missing = [c for c in ("Protein", "log2FoldChange", "pvalue", "padj", "n_numerator", "n_denominator",
                           "evidence_level", "method", "fdr_method", "significant") if c not in table.columns]
    if missing:
        raise ContractError(f"差异蛋白表缺少列：{', '.join(missing)}")
    proteins = table["Protein"].astype(str)
    tested = (table["evidence_level"] == "inferential").to_numpy()
    statistic = table["tstat"] if "tstat" in table.columns else pd.Series(np.nan, index=table.index)
    out = pd.DataFrame({
        "contrast_id": contrast_id,
        "contrast_role": contrast_role,
        "tissue": tissue,
        "omics": Omics.PROTEOMICS.value,
        "feature_id": proteins,
        "log2fc": pd.to_numeric(table["log2FoldChange"], errors="coerce"),
        "statistic": pd.to_numeric(statistic, errors="coerce"),
        "pvalue": pd.to_numeric(table["pvalue"], errors="coerce"),
        "padj": pd.to_numeric(table["padj"], errors="coerce"),
        "padj_method": np.where(tested, table["fdr_method"].astype(str), "none"),
        "n_numerator": pd.to_numeric(table["n_numerator"]).fillna(0).astype("int64"),
        "n_denominator": pd.to_numeric(table["n_denominator"]).fillna(0).astype("int64"),
        "method": table["method"].map(lambda m: _METHODS.get(str(m), str(m))),
        "mean_numerator": mean_numerator.reindex(proteins).to_numpy(dtype=float),
        "mean_denominator": mean_denominator.reindex(proteins).to_numpy(dtype=float),
        "units": LOG2_UNITS,
        "inference_unit": "biological_replicate",
        "evidence_level": np.where(tested, "statistical_test", "insufficient"),
        "selected": table["significant"].fillna(False).astype(bool),
        "selection_rule": selection_rule,
        "flags": [_flag(str(s), t) for s, t in zip(table.get("inference_status", pd.Series("", index=table.index)),
                                                   tested)],
    }, index=table.index)
    for source, target in _EXTRA.items():
        if source in table.columns:
            out[target] = table[source]
    return out.reset_index(drop=True)
