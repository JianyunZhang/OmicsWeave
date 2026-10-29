#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mouse Protein UniProt / EC Number Annotation Dictionary

向后兼容重导出模块 — 实际数据已迁移至 omicsweave.common.protein_id_mapping。
所有现有 import 路径无需修改即可继续使用。

Usage:
    from omicsweave.common.protein_id_mapping import (  # 推荐新路径
        MOUSE_PROTEIN_UNIPROT_MAPPING,
        get_uniprot_id, get_ec_for_protein,
    )
    # 旧路径仍可用:
    # from omicsweave.multiomics.mappings.protein_annotations import ...
"""
from __future__ import annotations

from omicsweave.common.mouse_protein_id_mapping import (
    MOUSE_PROTEIN_UNIPROT_MAPPING,
    get_uniprot_id,
    get_ec_for_protein,
)

__all__ = [
    "MOUSE_PROTEIN_UNIPROT_MAPPING",
    "get_uniprot_id",
    "get_ec_for_protein",
]
