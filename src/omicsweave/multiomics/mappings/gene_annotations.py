#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mouse Gene Annotations — NCBI Gene ID / Ensembl / Alias Resolution

向后兼容重导出模块 — 实际数据已迁移至 omicsweave.common.mouse_gene_id_mapping。
所有现有 import 路径无需修改即可继续使用。

Usage:
    from omicsweave.common.mouse_gene_id_mapping import (
        GENE_ANNOTATIONS, GENE_ALIAS_LOOKUP,
        resolve_gene_symbol, get_ncbi_gene_id, get_ensembl_gene_id,
    )
    # 旧路径仍可用:
    # from omicsweave.multiomics.mappings.gene_annotations import ...
"""
from __future__ import annotations

from omicsweave.common.mouse_gene_id_mapping import (
    GENE_ANNOTATIONS,
    GENE_ALIAS_LOOKUP,
    resolve_gene_symbol,
    get_ncbi_gene_id,
    get_ensembl_gene_id,
    get_gene_name,
    get_gene_info,
)

__all__ = [
    "GENE_ANNOTATIONS",
    "GENE_ALIAS_LOOKUP",
    "resolve_gene_symbol",
    "get_ncbi_gene_id",
    "get_ensembl_gene_id",
    "get_gene_name",
    "get_gene_info",
]
