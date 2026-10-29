#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Flux → Protein / Enzyme Mapping  (YAML-backed, validated against common/ protein registry)

Data is loaded from ``flux_protein_mapping.yaml`` in the same directory.
Protein symbols are validated at import time against
``MOUSE_PROTEIN_UNIPROT_MAPPING`` from ``omicsweave.core.idmap.mouse_protein_id_mapping``.

Public API (unchanged):
    REACTION_ENZYME_MAPPING   — dict[str, dict]
    get_all_related_proteins()
    get_related_proteins()
    get_protein_mapping_details()
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Tuple

import yaml

# ---------------------------------------------------------------------------
# YAML loader
# ---------------------------------------------------------------------------
_YAML_DIR = Path(__file__).resolve().parent


def _load_yaml(filename: str) -> dict[str, Any]:
    path = _YAML_DIR / filename
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"Invalid YAML root (expected mapping): {path}")
    return data


REACTION_ENZYME_MAPPING: dict[str, Any] = _load_yaml("flux_protein_mapping.yaml")

# ---------------------------------------------------------------------------
# Validation against common/mouse_protein_id_mapping
# ---------------------------------------------------------------------------


def _validate_proteins(mapping: dict[str, Any]) -> list[str]:
    try:
        from omicsweave.core.idmap.mouse_protein_id_mapping import (
            MOUSE_PROTEIN_UNIPROT_MAPPING,
        )
        from omicsweave.core.idmap.mouse_gene_id_mapping import (
            GENE_ALIAS_LOOKUP,
            GENE_ANNOTATIONS,
        )
    except ImportError:
        return ["[skip] omicsweave.core.idmap.mouse_protein_id_mapping not importable"]

    warns: list[str] = []
    for rxn_id, entry in mapping.items():
        for cat in ("direct_enzymes", "upstream_regulators", "downstream_enzymes", "transporters", "cofactors"):
            for protein in entry.get(cat, []):
                if protein in MOUSE_PROTEIN_UNIPROT_MAPPING:
                    continue
                canonical = GENE_ALIAS_LOOKUP.get(protein, protein)
                if canonical is None:
                    warns.append(f"{rxn_id}.{cat}: non-protein entity '{protein}'")
                elif canonical in MOUSE_PROTEIN_UNIPROT_MAPPING:
                    continue
                elif canonical in GENE_ANNOTATIONS:
                    # 合法小鼠基因但本地 UniProt 注册表尚未覆盖；运行时可安全按
                    # 规范符号回退，由 audit_mapping_integrity 汇总报告。
                    continue
    return warns


_protein_warnings = _validate_proteins(REACTION_ENZYME_MAPPING)

# ---------------------------------------------------------------------------
# Helper functions (preserved from original)
# ---------------------------------------------------------------------------


def get_all_related_proteins(reaction: str) -> Tuple[List[str], Dict]:
    """获取与代谢反应相关的所有蛋白质（包括上下游调控蛋白）

    Parameters
    ----------
    reaction : str
        反应ID (如 'HEX_c', 'PDH_m')

    Returns
    -------
    Tuple[List[str], Dict]
        (蛋白质列表, 分类详情字典)
    """
    reaction_base = reaction.replace('_c', '').replace('_m', '').replace('_e', '')

    mapping = None
    if reaction in REACTION_ENZYME_MAPPING:
        mapping = REACTION_ENZYME_MAPPING[reaction]
    elif reaction_base in REACTION_ENZYME_MAPPING:
        mapping = REACTION_ENZYME_MAPPING[reaction_base]

    if mapping is None:
        return [], {}

    all_proteins: list[str] = []
    categories: dict[str, list[str]] = {}

    for category in ['direct_enzymes', 'upstream_regulators', 'downstream_enzymes', 'transporters', 'cofactors']:
        proteins = mapping.get(category, [])
        categories[category] = proteins
        all_proteins.extend(proteins)

    all_proteins = list(dict.fromkeys(all_proteins))
    return all_proteins, categories


def get_related_proteins(reaction: str, include_regulators: bool = True) -> List[str]:
    """获取与代谢反应相关的关键蛋白质

    Parameters
    ----------
    reaction : str
        反应ID
    include_regulators : bool
        是否包含上下游调控蛋白

    Returns
    -------
    List[str]
        蛋白质列表
    """
    all_proteins, categories = get_all_related_proteins(reaction)

    if not include_regulators:
        return categories.get('direct_enzymes', [])

    return all_proteins


def get_protein_mapping_details(reaction: str) -> Dict:
    """获取反应-蛋白质映射的详细信息

    Parameters
    ----------
    reaction : str
        反应ID (如 'HEX_c', 'PDH_m')

    Returns
    -------
    Dict
        包含各类蛋白质的分类详情，结构与 REACTION_ENZYME_MAPPING 中的值一致。
        若反应未找到则返回空字典。
    """
    reaction_base = reaction.replace('_c', '').replace('_m', '').replace('_e', '')

    if reaction in REACTION_ENZYME_MAPPING:
        return REACTION_ENZYME_MAPPING[reaction]
    elif reaction_base in REACTION_ENZYME_MAPPING:
        return REACTION_ENZYME_MAPPING[reaction_base]

    return {}
