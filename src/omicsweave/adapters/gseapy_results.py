#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gseapy 结果 → 统一富集结果表（``omicsweave.contract.tables.ENRICHMENT_COLUMNS``）。

- Enrichr（ORA）：``Gene_set / Term / Overlap / P-value / Adjusted P-value / Odds Ratio / Genes``；
  ``Overlap`` 形如 ``5/120``（命中数 / 通路大小）。Enrichr 的校正 p 值为 BH 法。
- prerank（预排序 GSEA）：``Term / ES / NES / NOM p-val / FDR q-val / Tag % / Lead_genes``；
  ``Tag %`` 形如 ``10/40``（leading edge 数 / 通路大小）。同时用多个基因集库时，
  ``Term`` 带 ``<库名>__`` 前缀。``FDR q-val`` 是 GSEA 自己的 FDR，记为 ``gsea_fdr``。
  没有一次置换达到观测富集分数时，gseapy 报告的 ``NOM p-val`` / ``FDR q-val`` 为 0；置换 p 值不应为 0，
  已知置换次数 n 时改为 1/(n+1)，并在 ``flags`` 中记 ``pvalue_floored`` / ``padj_floored``。

既接受 gseapy 返回的 DataFrame，也接受由其保存的 CSV 读回的表。
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from omicsweave.contract.dataset import ContractError
from omicsweave.contract.tables import ORA_QUERIES

_GO_ID = r"\((GO:\d{7})\)\s*$"


def _split_ratio(series: pd.Series) -> tuple[pd.Series, pd.Series]:
    parts = series.astype(str).str.extract(r"^\s*(\d+)\s*/\s*(\d+)\s*$")
    return pd.to_numeric(parts[0]).astype("Int64"), pd.to_numeric(parts[1]).astype("Int64")


def _optional(table: pd.DataFrame, name: str, default) -> pd.Series:
    return table[name] if name in table.columns else pd.Series(default, index=table.index)


def _require(table: pd.DataFrame, columns: list[str], source: str) -> None:
    missing = [c for c in columns if c not in table.columns]
    if missing:
        raise ContractError(f"{source} 结果缺少列：{', '.join(missing)}")


def enrichr_to_contract(results: pd.DataFrame, *, contrast_id: str, contrast_role: str, tissue: str,
                        omics: str, query: str, query_size: int | None = None) -> pd.DataFrame:
    """gseapy Enrichr 结果 → 富集结果表。``query`` 为所用基因列表：up / down / all。"""
    if query not in ORA_QUERIES:
        raise ContractError(f"query 应为 {' / '.join(ORA_QUERIES)}")
    _require(results, ["Gene_set", "Term", "Overlap", "P-value", "Adjusted P-value"], "Enrichr")
    overlap, set_size = _split_ratio(results["Overlap"])
    terms = results["Term"].astype(str)
    return pd.DataFrame({
        "contrast_id": contrast_id,
        "contrast_role": contrast_role,
        "tissue": tissue,
        "omics": omics,
        "method": "ora",
        "query": query,
        "gene_set_library": results["Gene_set"].astype(str),
        "term_id": terms.str.extract(_GO_ID)[0].fillna(""),
        "term_name": terms,
        "pvalue": pd.to_numeric(results["P-value"], errors="coerce"),
        "padj": pd.to_numeric(results["Adjusted P-value"], errors="coerce"),
        "padj_method": "fdr_bh",
        "effect_direction": ORA_QUERIES[query],
        "overlap_size": overlap,
        "set_size": set_size,
        "query_size": pd.array([query_size] * len(results), dtype="Int64"),
        "universe_size": pd.Series(pd.NA, index=results.index, dtype="Int64"),
        "odds_ratio": pd.to_numeric(_optional(results, "Odds Ratio", np.nan), errors="coerce"),
        "es": np.nan,
        "nes": np.nan,
        "hit_features": _optional(results, "Genes", "").fillna("").astype(str),
        "tool": "gseapy.enrichr",
        "flags": "",
    })


def _permutation_floor(values: pd.Series, permutation_num: int | None) -> tuple[pd.Series, pd.Series]:
    """为 0 的置换 p 值 / FDR 改为 1/(n+1)；返回新值与被改动的行。置换次数未知时不改。"""
    floored = (values == 0) & (permutation_num is not None)
    return (values.mask(floored, 1.0 / (permutation_num + 1)) if floored.any() else values), floored


def prerank_to_contract(res2d: pd.DataFrame, *, contrast_id: str, contrast_role: str, tissue: str,
                        omics: str, gene_set_library: str | None = None,
                        permutation_num: int | None = None) -> pd.DataFrame:
    """gseapy prerank 的 ``res2d`` → 富集结果表。

    只用一个基因集库时 ``Term`` 不带库名前缀，需要用 ``gene_set_library`` 指定库名。``permutation_num``
    为本次 GSEA 的置换次数，给出时为 0 的 p 值与 FDR 改为 1/(n+1) 并记入 ``flags``；置换次数未知时保留原值。
    """
    _require(res2d, ["Term", "ES", "NES", "NOM p-val", "FDR q-val"], "prerank")
    terms = res2d["Term"].astype(str)
    prefixed = terms.str.contains("__", regex=False)
    if gene_set_library is None and not prefixed.all():
        raise ContractError("Term 不带基因集库前缀，请用 gene_set_library 指定库名")
    split = terms.str.split("__", n=1, expand=True) if prefixed.any() else None
    library = (split[0].where(prefixed, gene_set_library) if split is not None
               else pd.Series(gene_set_library, index=res2d.index))
    names = split[1].where(prefixed, terms) if split is not None else terms

    nes = pd.to_numeric(res2d["NES"], errors="coerce")
    overlap, set_size = (_split_ratio(res2d["Tag %"]) if "Tag %" in res2d.columns
                         else (pd.Series(pd.NA, index=res2d.index, dtype="Int64"),) * 2)
    pvalue, pvalue_floored = _permutation_floor(pd.to_numeric(res2d["NOM p-val"], errors="coerce"), permutation_num)
    padj, padj_floored = _permutation_floor(pd.to_numeric(res2d["FDR q-val"], errors="coerce"), permutation_num)
    flags = [";".join(mark for mark, hit in (("pvalue_floored", p), ("padj_floored", q)) if hit)
             for p, q in zip(pvalue_floored, padj_floored)]
    return pd.DataFrame({
        "contrast_id": contrast_id,
        "contrast_role": contrast_role,
        "tissue": tissue,
        "omics": omics,
        "method": "gsea_prerank",
        "query": "ranked",
        "gene_set_library": library.astype(str),
        "term_id": names.str.extract(_GO_ID)[0].fillna(""),
        "term_name": names,
        "pvalue": pvalue,
        "padj": padj,
        "padj_method": "gsea_fdr",
        "effect_direction": np.select([nes > 0, nes < 0], ["up", "down"], "none"),
        "overlap_size": overlap,
        "set_size": set_size,
        "query_size": pd.Series(pd.NA, index=res2d.index, dtype="Int64"),
        "universe_size": pd.Series(pd.NA, index=res2d.index, dtype="Int64"),
        "odds_ratio": np.nan,
        "es": pd.to_numeric(res2d["ES"], errors="coerce"),
        "nes": nes,
        "hit_features": _optional(res2d, "Lead_genes", "").fillna("").astype(str),
        "tool": "gseapy.prerank",
        "flags": flags,
    })
