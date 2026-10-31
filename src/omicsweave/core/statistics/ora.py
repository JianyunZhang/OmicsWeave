#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
过度代表分析（ORA）：单侧 Fisher 精确检验，在本地计算，不调用在线服务。

背景（universe）应取实际参与检验的特征（如差异分析中有 p 值的基因），而不是基因集库里的全部基因；
基因集先与背景取交集，交集大小在 ``[min_size, max_size]`` 内的才参与检验，所有参与检验的基因集
（包括命中数为 0 的）都进入多重检验校正。与代谢组 KEGG ORA 的做法一致。
"""
from __future__ import annotations

from collections.abc import Iterable, Mapping

import numpy as np
import pandas as pd
from scipy import stats

FDR_METHODS = ("fdr_bh", "fdr_by", "bonferroni", "none")
COLUMNS = ["term_name", "pvalue", "padj", "overlap_size", "set_size", "query_size", "universe_size",
           "odds_ratio", "hit_features"]


def adjust_pvalues(pvalues: np.ndarray, method: str = "fdr_bh") -> np.ndarray:
    """多重检验校正（Benjamini-Hochberg / Benjamini-Yekutieli / Bonferroni），只依赖 numpy。"""
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    if method == "none" or n == 0:
        return p.copy()
    if method == "bonferroni":
        return np.minimum(p * n, 1.0)
    if method not in ("fdr_bh", "fdr_by"):
        raise ValueError(f"不支持的校正方法 {method!r}，可选：{'、'.join(FDR_METHODS)}")
    scale = np.sum(1.0 / np.arange(1, n + 1)) if method == "fdr_by" else 1.0
    order = np.argsort(p, kind="stable")
    ranked = p[order] * n * scale / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty(n)
    adjusted[order] = np.minimum(ranked, 1.0)
    return adjusted


def over_representation(query: Iterable[str], universe: Iterable[str], gene_sets: Mapping[str, Iterable[str]], *,
                        min_size: int = 5, max_size: int = 2000, fdr_method: str = "fdr_bh") -> pd.DataFrame:
    """一个基因集库的 ORA；返回按 padj、pvalue 排序的表（列见 :data:`COLUMNS`）。"""
    universe = set(map(str, universe))
    query = set(map(str, query)) & universe
    if min_size < 1 or max_size < min_size:
        raise ValueError("ORA 基因集大小参数无效")
    rows = []
    for term, members in gene_sets.items():
        members = set(map(str, members)) & universe
        if not min_size <= len(members) <= max_size:
            continue
        hits = query & members
        a = len(hits)
        b = len(query) - a
        c = len(members) - a
        d = len(universe) - a - b - c
        odds_ratio, pvalue = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append({"term_name": str(term), "pvalue": float(pvalue), "overlap_size": a, "set_size": len(members),
                     "query_size": len(query), "universe_size": len(universe), "odds_ratio": float(odds_ratio),
                     "hit_features": ";".join(sorted(hits))})
    if not rows:
        return pd.DataFrame(columns=COLUMNS)
    table = pd.DataFrame(rows)
    table["padj"] = adjust_pvalues(table["pvalue"].to_numpy(dtype=float), fdr_method)
    table["odds_ratio"] = table["odds_ratio"].replace([np.inf], np.nan)
    return table[COLUMNS].sort_values(["padj", "pvalue", "term_name"], kind="stable").reset_index(drop=True)
