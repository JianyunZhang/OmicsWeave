#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Metabolite Database ID Annotations (HMDB / KEGG Compound / ChEBI)

向后兼容重导出模块 — 实际数据已迁移至 omicsweave.common.metabolite_id_mapping。
所有现有 import 路径无需修改即可继续使用。

Usage:
    from omicsweave.multiomics.mappings.metabolite_annotations import METABOLITE_DATABASE_IDS
    ids = METABOLITE_DATABASE_IDS.get('Citrate')
    # -> {'hmdb_id': 'HMDB0000094', 'kegg_compound': 'C00158', 'chebi_id': 'CHEBI:16947'}
"""
from __future__ import annotations

from omicsweave.common.metabolite_id_mapping import (
    METABOLITE_DATABASE_IDS,
    lookup_metabolite_ids,
    annotate_metabolites,
    get_hmdb_id,
    get_kegg_compound,
    get_chebi_id,
)

__all__ = [
    "METABOLITE_DATABASE_IDS",
    "lookup_metabolite_ids",
    "annotate_metabolites",
    "get_hmdb_id",
    "get_kegg_compound",
    "get_chebi_id",
]
