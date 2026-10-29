#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
基因注释增强 — 为转录组学输出 CSV 添加 Ensembl Gene ID 和 NCBI Gene ID

在不变动现有转录组学流程结构的前提下，对 gene_counts_matrix.csv 和
DE 结果 CSV 增加数据库交叉引用字段，便于与 MFA-gene mapping 字典对接。

使用方式:
    from gene_annotate import annotate_gene_dataframe

    # 增强 gene_counts_matrix (index=gene_symbol)
    enriched = annotate_gene_dataframe(gene_counts_df)
    enriched.to_csv("gene_counts_matrix_annotated.csv")

    # 增强 DE 结果 (index 或 'gene' 列 = gene_symbol)
    enriched = annotate_gene_dataframe(de_results_df)
    enriched.to_csv("de_results_annotated.csv")

依赖:
    omicsweave/core/idmap/mouse_gene_id_mapping.py
"""

from __future__ import annotations

import os
from typing import Optional

import pandas as pd



def _get_annotation_funcs():
    """延迟导入共用映射模块（避免循环依赖）。"""
    from omicsweave.core.idmap.mouse_gene_id_mapping import (
        get_ncbi_gene_id,
        get_ensembl_gene_id,
        get_gene_name,
        resolve_gene_symbol,
    )
    return get_ncbi_gene_id, get_ensembl_gene_id, get_gene_name, resolve_gene_symbol


def annotate_gene_dataframe(
    df: pd.DataFrame,
    gene_col: Optional[str] = None,
    add_gene_name: bool = False,
) -> pd.DataFrame:
    """
    为基因数据的 DataFrame 添加 Ensembl Gene ID, NCBI Gene ID,
    canonical_symbol 和 mapping_status 列。

    Parameters
    ----------
    df : pd.DataFrame
        输入数据框，基因符号位于 index（默认）或 gene_col 指定的列中。
    gene_col : str, optional
        如果基因符号在某个列中（如 DE 结果的 'gene' 列），指定列名。
        如果为 None，则使用 df.index 作为基因符号。
    add_gene_name : bool
        是否也添加 gene_name（基因全名）列。默认 False 以保持输出简洁。

    Returns
    -------
    pd.DataFrame
        增强后的数据框，增加了 gene_symbol, canonical_symbol, ensembl_gene_id,
        ncbi_gene_id, mapping_status 列。
        匹配不到数据库 ID 的行填 None。
    """
    get_ncbi, get_ensembl, get_name, resolve = _get_annotation_funcs()

    if gene_col is not None:
        symbols = df[gene_col].astype(str).tolist()
    else:
        symbols = df.index.astype(str).tolist()

    # 解析基因符号（先别名解析，再查 NCBI/Ensembl）
    canonical_symbols = []
    ncbi_ids = []
    ensembl_ids = []
    mapping_statuses = []

    for s in symbols:
        canonical = resolve(s)
        canonical_symbols.append(canonical)
        ncbi_id = get_ncbi(canonical if canonical else s)
        ensembl_id = get_ensembl(canonical if canonical else s)
        ncbi_ids.append(ncbi_id)
        ensembl_ids.append(ensembl_id)

        if canonical and canonical != s:
            mapping_statuses.append("alias_resolved")
        elif ncbi_id or ensembl_id:
            mapping_statuses.append("direct_match")
        else:
            mapping_statuses.append("unmapped")

    result_df = df.copy()

    # Insert annotation columns after 'significant' for DE results, or at end for matrices
    if gene_col is not None:
        insert_after = 'significant' if 'significant' in result_df.columns else result_df.columns[-1]
        idx = result_df.columns.get_loc(insert_after) + 1
        for i, (col_name, values) in enumerate([
            ('gene_symbol', symbols),
            ('canonical_symbol', canonical_symbols),
            ('ensembl_gene_id', ensembl_ids),
            ('ncbi_gene_id', ncbi_ids),
            ('mapping_status', mapping_statuses),
        ]):
            if col_name not in result_df.columns:
                result_df.insert(idx + i, col_name, values)
    else:
        result_df['ensembl_gene_id'] = ensembl_ids
        result_df['ncbi_gene_id'] = ncbi_ids
        result_df['canonical_symbol'] = canonical_symbols
        result_df['mapping_status'] = mapping_statuses
        # Move annotation columns to end
        cols = list(result_df.columns)
        annot_cols = ['canonical_symbol', 'ensembl_gene_id', 'ncbi_gene_id', 'mapping_status']
        data_cols = [c for c in cols if c not in annot_cols]
        result_df = result_df[data_cols + annot_cols]

    if add_gene_name:
        gene_names = [get_name(s) for s in symbols]
        result_df['gene_name'] = gene_names

    # 统计覆盖率
    matched_ncbi = sum(1 for x in ncbi_ids if x is not None)
    matched_ensembl = sum(1 for x in ensembl_ids if x is not None)
    matched_alias = sum(1 for x in mapping_statuses if x == "alias_resolved")
    total = len(symbols)
    print(f"  基因注释: {matched_ncbi}/{total} 匹配 NCBI ID, "
          f"{matched_ensembl}/{total} 匹配 Ensembl ID, "
          f"{matched_alias}/{total} 别名解析")

    return result_df


def annotate_csv_file(
    csv_path: str,
    gene_col: Optional[str] = None,
    overwrite: bool = True,
) -> str:
    """
    读取 CSV 文件，添加基因注释列，保存回原文件。

    默认覆盖原文件（overwrite=True），将 gene_symbol, canonical_symbol,
    ensembl_gene_id, ncbi_gene_id, mapping_status 列写入主文件。

    同时生成 *_id_mapped.csv（仅含 ID 列，供多组学直接对齐）。

    Parameters
    ----------
    csv_path : str
        输入 CSV 文件路径。
    gene_col : str, optional
        基因符号列名（列式数据）或 None（index 为基因符号）。
    overwrite : bool
        是否覆盖原文件。默认 True。

    Returns
    -------
    str
        输出文件路径。
    """
    df = pd.read_csv(csv_path, index_col=0 if gene_col is None else None)

    if gene_col is not None and gene_col not in df.columns:
        # 如果指定的列不存在，尝试用 index
        if df.index.name == gene_col or gene_col == 'gene':
            gene_col = None

    annotated = annotate_gene_dataframe(df, gene_col=gene_col)

    if overwrite:
        out_path = csv_path
    else:
        base, ext = os.path.splitext(csv_path)
        out_path = f"{base}_annotated{ext}"

    annotated.to_csv(out_path)
    print(f"  注释文件已保存: {out_path}")

    # Generate *_id_mapped.csv for multi-omics alignment
    id_cols = [c for c in ['gene_symbol', 'canonical_symbol', 'ensembl_gene_id', 'ncbi_gene_id', 'mapping_status']
               if c in annotated.columns]
    if id_cols:
        base, ext = os.path.splitext(csv_path)
        mapped_path = f"{base}_id_mapped{ext}"
        annotated[id_cols].to_csv(mapped_path, index=False)
        print(f"  ID 映射文件已保存: {mapped_path}")

    return out_path
