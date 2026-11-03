#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
基于差异结果的基因集富集，输出统一富集结果表（``omicsweave.contract.tables.ENRICHMENT_COLUMNS``）。供转录组、
蛋白组运行器调用；代谢组的 KEGG 代谢物集富集在代谢组模块中实现。

- ORA：入选特征按方向（上调 / 下调）分别检验，背景为该比较中有 p 值的特征，在本地用单侧 Fisher 检验计算
  （:mod:`omicsweave.core.statistics.ora`），不把基因列表上传到 Enrichr 等在线服务；
- 预排序 GSEA：``gseapy.prerank``。多个基因集库可以合在一次检验中（条目名带 ``<库名>__`` 前缀，GSEA 的 FDR
  在全部基因集上计算），也可以每个库分别检验。
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import pandas as pd

from omicsweave.contract.tables import ENRICHMENT_COLUMNS
from omicsweave.core.statistics.ora import over_representation

from .gseapy_results import prerank_to_contract

_GO_ID = r"\((GO:\d{7})\)\s*$"
Libraries = Mapping[str, Mapping[str, list[str]]]


@dataclass(frozen=True)
class EnrichmentSettings:
    # ORA：入选特征总数不足 ora_min_selected 时不做；每个方向至少 ora_min_query 个特征
    ora_min_selected: int = 5
    ora_min_query: int = 3
    ora_min_set_size: int = 5
    ora_max_set_size: int = 2000
    fdr_method: str = "fdr_bh"
    # GSEA：排序特征不足 gsea_min_ranked 个时不做
    gsea_min_ranked: int = 100
    gsea_min_set_size: int = 15
    gsea_max_set_size: int = 500
    gsea_permutation_num: int = 1000
    gsea_seed: int = 42
    gsea_threads: int = 4
    gsea_joint_libraries: bool = True


def _empty() -> pd.DataFrame:
    return pd.DataFrame(columns=[c.name for c in ENRICHMENT_COLUMNS])


def ora_rows(differential: pd.DataFrame, libraries: Libraries, settings: EnrichmentSettings, *,
             contrast_id: str, contrast_role: str, tissue: str, omics: str) -> pd.DataFrame:
    """一个比较、一个组织的统一差异结果表 → ORA 的富集结果行。"""
    selected = differential[differential["selected"]]
    if len(selected) < settings.ora_min_selected:
        return _empty()
    universe = differential.loc[differential["pvalue"].notna(), "feature_id"]
    parts = []
    for query, mask in (("up", selected["log2fc"] > 0), ("down", selected["log2fc"] < 0)):
        features = selected.loc[mask, "feature_id"]
        if len(features) < settings.ora_min_query:
            continue
        for library, gene_sets in libraries.items():
            result = over_representation(features, universe, gene_sets, min_size=settings.ora_min_set_size,
                                         max_size=settings.ora_max_set_size, fdr_method=settings.fdr_method)
            if result.empty:
                continue
            parts.append(pd.DataFrame({
                "contrast_id": contrast_id, "contrast_role": contrast_role, "tissue": tissue, "omics": omics,
                "method": "ora", "query": query, "gene_set_library": library,
                "term_id": result["term_name"].str.extract(_GO_ID)[0].fillna(""),
                "term_name": result["term_name"], "pvalue": result["pvalue"], "padj": result["padj"],
                "padj_method": settings.fdr_method, "effect_direction": query,
                "overlap_size": result["overlap_size"].astype("Int64"), "set_size": result["set_size"].astype("Int64"),
                "query_size": result["query_size"].astype("Int64"),
                "universe_size": result["universe_size"].astype("Int64"),
                "odds_ratio": result["odds_ratio"], "hit_features": result["hit_features"],
                "tool": "omicsweave.ora",
            }))
    return pd.concat(parts, ignore_index=True) if parts else _empty()


def prerank_rows(ranked: pd.Series, libraries: Libraries, settings: EnrichmentSettings, *,
                 contrast_id: str, contrast_role: str, tissue: str, omics: str) -> pd.DataFrame:
    """按统计量降序排列的特征（索引为特征 ID）→ 预排序 GSEA 的富集结果行。"""
    import gseapy

    if len(ranked) < settings.gsea_min_ranked:
        return _empty()
    if settings.gsea_joint_libraries and len(libraries) > 1:
        jobs = [(None, {f"{lib}__{term}": list(members)
                        for lib, sets in libraries.items() for term, members in sets.items()})]
    else:
        jobs = [(lib, {term: list(members) for term, members in sets.items()}) for lib, sets in libraries.items()]
    parts = []
    for library, gene_sets in jobs:
        result = gseapy.prerank(rnk=ranked, gene_sets=gene_sets, outdir=None, min_size=settings.gsea_min_set_size,
                                max_size=settings.gsea_max_set_size, permutation_num=settings.gsea_permutation_num,
                                seed=settings.gsea_seed, threads=settings.gsea_threads, ascending=False,
                                no_plot=True, verbose=False)
        if result.res2d is None or result.res2d.empty:
            continue
        parts.append(prerank_to_contract(result.res2d, contrast_id=contrast_id, contrast_role=contrast_role,
                                         tissue=tissue, omics=omics, gene_set_library=library,
                                         permutation_num=settings.gsea_permutation_num))
    return pd.concat(parts, ignore_index=True) if parts else _empty()
