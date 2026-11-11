#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""小鼠蛋白注册表：MGI 符号 → UniProtKB/Swiss-Prot 登记号、条目名与 EC 号。

内容由 ``python -m omicsweave.core.idmap.build_registry`` 从 UniProtKB/Swiss-Prot 小鼠条目生成，
按 MGI 登记号对齐到 MGI 符号，存放在 ``data/``；来源与版本见
``omicsweave.core.idmap.mouse_gene_id_mapping.REGISTRY_SOURCES``。不要手工修改。

``MOUSE_PROTEIN_UNIPROT_MAPPING[符号]`` 含 ``uniprot_id``（多个条目时取第一个）、``uniprot_ids``（全部条目）、
``uniprot_entry_name``、``ec_number``、``protein_name``；``UNIPROT_TO_SYMBOL`` 为登记号 → 符号的反查表。

用法::

    from omicsweave.core.idmap.mouse_protein_id_mapping import get_uniprot_id
    get_uniprot_id('Gpd1')   # -> 'P13707'
"""

from __future__ import annotations

from typing import Any, Dict, Optional, Tuple

from ._registry import read_table


def _load() -> tuple[dict, dict]:
    mapping: Dict[str, Dict[str, Any]] = {}
    reverse: Dict[str, str] = {}
    for row in read_table("mouse_protein_registry.tsv.gz"):
        symbol, accession = row["symbol"], row["uniprot_id"]
        reverse[accession] = symbol
        entry = mapping.get(symbol)
        if entry is None:
            mapping[symbol] = {
                "uniprot_id": accession,
                "uniprot_ids": (accession,),
                "uniprot_entry_name": row["entry_name"] or None,
                "ec_number": row["ec_number"] or None,
                "protein_name": row["protein_name"] or None,
            }
        else:
            entry["uniprot_ids"] = (*entry["uniprot_ids"], accession)
    return mapping, reverse


MOUSE_PROTEIN_UNIPROT_MAPPING: Dict[str, Dict[str, Any]]
UNIPROT_TO_SYMBOL: Dict[str, str]
MOUSE_PROTEIN_UNIPROT_MAPPING, UNIPROT_TO_SYMBOL = _load()


def get_uniprot_id(gene_symbol: str) -> Optional[str]:
    """小鼠基因符号的 UniProt 登记号；有多个 Swiss-Prot 条目时返回第一个。"""
    entry = MOUSE_PROTEIN_UNIPROT_MAPPING.get(gene_symbol)
    return entry["uniprot_id"] if entry else None


def get_uniprot_ids(gene_symbol: str) -> Tuple[str, ...]:
    """小鼠基因符号的全部 Swiss-Prot 登记号。"""
    entry = MOUSE_PROTEIN_UNIPROT_MAPPING.get(gene_symbol)
    return tuple(entry["uniprot_ids"]) if entry else ()


def get_ec_for_protein(gene_symbol: str) -> Optional[str]:
    """小鼠基因符号的 EC 号（多个时以 ``"; "`` 连接）。"""
    entry = MOUSE_PROTEIN_UNIPROT_MAPPING.get(gene_symbol)
    return entry["ec_number"] if entry else None
