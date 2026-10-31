#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
代谢组结果 → 统一结果表（``omicsweave.contract.tables``）。

- 差异分析：:func:`omicsweave.metabolomics.differential_analysis.analyze_differential_metabolites`
  每个比较一张表（``contrasts[比较名]``），在 log2 归一化峰面积上做 Welch t 检验；
- 通路 ORA：:func:`omicsweave.metabolomics.enrichment._run.calculate_pathway_enrichment` 的输出列
  ``ID / Name / P-value / OddsRatio / SigCount / BgCount / SigTotal / UniverseCount / Members / FDR``。
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from omicsweave.contract.dataset import ContractError
from omicsweave.contract.tables import ORA_QUERIES
from omicsweave.core.model import Omics

from .data_preprocessing import normalize_fdr_method

_REQUIRED = ("ID", "Name", "P-value", "FDR", "SigCount", "BgCount")
_DIFFERENTIAL_REQUIRED = ("log2FoldChange", "statistic", "pvalue", "padj", "significant", "evidence_level")
# 原表中保留为附加列的注释与判定列：原列名 → 统一表中的列名
_DIFFERENTIAL_EXTRA = {
    "VIP_Global": "vip_global", "nominal_significant": "nominal_significant",
    "significant_univariate": "significant_univariate", "significant_with_vip": "significant_with_vip",
    "Formula": "formula", "kegg_map": "kegg_map", "KEGG_Compound": "kegg_compound", "HMDB": "hmdb",
    "ChEBI": "chebi",
}
LOG2_UNITS = "log2 normalized peak area"


def differential_to_contract(table: pd.DataFrame, *, contrast_id: str, contrast_role: str, tissue: str,
                             feature_id_column: str, feature_name_column: str | None,
                             mean_numerator: pd.Series, mean_denominator: pd.Series,
                             n_numerator: int, n_denominator: int, fdr_method: str,
                             selection_rule: str) -> pd.DataFrame:
    """把一个比较的差异表转换为统一差异结果表。均值为两组 log2 归一化峰面积的均值，与 ``table`` 行对齐。"""
    missing = [c for c in (feature_id_column, *_DIFFERENTIAL_REQUIRED) if c not in table.columns]
    if missing:
        raise ContractError(f"代谢组差异表缺少列：{', '.join(missing)}")
    tested = (table["evidence_level"] == "inferential").to_numpy()
    out = pd.DataFrame({
        "contrast_id": contrast_id,
        "contrast_role": contrast_role,
        "tissue": tissue,
        "omics": Omics.METABOLOMICS.value,
        "feature_id": table[feature_id_column].astype(str),
        "log2fc": pd.to_numeric(table["log2FoldChange"], errors="coerce"),
        "statistic": pd.to_numeric(table["statistic"], errors="coerce"),
        "pvalue": pd.to_numeric(table["pvalue"], errors="coerce"),
        "padj": pd.to_numeric(table["padj"], errors="coerce"),
        "padj_method": np.where(tested, normalize_fdr_method(fdr_method), "none"),
        "n_numerator": int(n_numerator),
        "n_denominator": int(n_denominator),
        "method": np.where(tested, "welch_ttest", "insufficient"),
        "mean_numerator": mean_numerator.reindex(table.index).to_numpy(dtype=float),
        "mean_denominator": mean_denominator.reindex(table.index).to_numpy(dtype=float),
        "units": LOG2_UNITS,
        "inference_unit": "biological_replicate",
        "evidence_level": np.where(tested, "statistical_test", "insufficient"),
        "selected": table["significant"].fillna(False).astype(bool),
        "selection_rule": selection_rule,
        "flags": np.where(tested, "", "low_n"),
    }, index=table.index)
    if feature_name_column and feature_name_column != feature_id_column and feature_name_column in table.columns:
        out.insert(out.columns.get_loc("feature_id") + 1, "feature_name", table[feature_name_column].astype(str))
    vip_column = f"VIP_{contrast_id}"
    if vip_column in table.columns:
        out["vip"] = pd.to_numeric(table[vip_column], errors="coerce")
    for source, target in _DIFFERENTIAL_EXTRA.items():
        if source in table.columns and source not in (feature_id_column, feature_name_column):
            out[target] = table[source]
    return out.reset_index(drop=True)


def ora_to_contract(result: pd.DataFrame, *, contrast_id: str, contrast_role: str, tissue: str, query: str,
                    gene_set_library: str = "KEGG", padj_method: str = "fdr_bh") -> pd.DataFrame:
    """``query`` 为参与 ORA 的显著代谢物集合：up / down / all。"""
    if query not in ORA_QUERIES:
        raise ContractError(f"query 应为 {' / '.join(ORA_QUERIES)}")
    missing = [c for c in _REQUIRED if c not in result.columns]
    if missing:
        raise ContractError(f"代谢组 ORA 结果缺少列：{', '.join(missing)}")

    def optional(name: str) -> pd.Series:
        return result[name] if name in result.columns else pd.Series(pd.NA, index=result.index)

    return pd.DataFrame({
        "contrast_id": contrast_id,
        "contrast_role": contrast_role,
        "tissue": tissue,
        "omics": Omics.METABOLOMICS.value,
        "method": "ora",
        "query": query,
        "gene_set_library": gene_set_library,
        "term_id": result["ID"].astype(str),
        "term_name": result["Name"].astype(str),
        "pvalue": pd.to_numeric(result["P-value"], errors="coerce"),
        "padj": pd.to_numeric(result["FDR"], errors="coerce"),
        "padj_method": padj_method,
        "effect_direction": ORA_QUERIES[query],
        "overlap_size": pd.to_numeric(result["SigCount"]).astype("Int64"),
        "set_size": pd.to_numeric(result["BgCount"]).astype("Int64"),
        "query_size": pd.to_numeric(optional("SigTotal")).astype("Int64"),
        "universe_size": pd.to_numeric(optional("UniverseCount")).astype("Int64"),
        "odds_ratio": pd.to_numeric(optional("OddsRatio"), errors="coerce"),
        "hit_features": optional("Members").fillna("").astype(str),
        "tool": "omicsweave.metabolomics.ora",
    }).reset_index(drop=True)
