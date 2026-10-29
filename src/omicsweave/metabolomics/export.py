#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
代谢组通路 ORA 结果 → 统一富集结果表（``omicsweave.contract.tables.ENRICHMENT_COLUMNS``）。

:func:`omicsweave.metabolomics.enrichment._run.calculate_pathway_enrichment` 的输出列：
``ID / Name / P-value / OddsRatio / SigCount / BgCount / SigTotal / UniverseCount / Members / FDR``。
"""
from __future__ import annotations

import pandas as pd

from omicsweave.contract.dataset import ContractError
from omicsweave.contract.tables import ORA_QUERIES
from omicsweave.core.model import Omics

_REQUIRED = ("ID", "Name", "P-value", "FDR", "SigCount", "BgCount")


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
