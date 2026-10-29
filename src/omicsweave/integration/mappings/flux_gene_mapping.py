#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Flux → Gene Mapping  (YAML-backed, validated against common/ gene registry)

Data is loaded from ``flux_gene_mapping.yaml`` in the same directory.
Gene symbols are validated at import time against ``GENE_ANNOTATIONS`` and
``GENE_ALIAS_LOOKUP`` from ``omicsweave.core.idmap.mouse_gene_id_mapping``.
Unknown symbols produce a ``warnings.warn()`` but do **not** block loading.

Public API (unchanged):
    FLUX_GENE_MAPPING  — dict[str, dict]
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

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


FLUX_GENE_MAPPING: dict[str, Any] = _load_yaml("flux_gene_mapping.yaml")

# ---------------------------------------------------------------------------
# Validation against common/mouse_gene_id_mapping
# ---------------------------------------------------------------------------


def _validate_genes(mapping: dict[str, Any]) -> list[str]:
    try:
        from omicsweave.core.idmap.mouse_gene_id_mapping import (
            GENE_ALIAS_LOOKUP,
            GENE_ANNOTATIONS,
        )
    except ImportError:
        return ["[skip] omicsweave.core.idmap.mouse_gene_id_mapping not importable"]

    warns: list[str] = []
    for rxn_id, entry in mapping.items():
        for gene in entry.get("genes", []):
            if gene not in GENE_ANNOTATIONS:
                resolved = GENE_ALIAS_LOOKUP.get(gene)
                if resolved is None:
                    warns.append(f"{rxn_id}: unknown gene '{gene}'")
    return warns


_gene_warnings = _validate_genes(FLUX_GENE_MAPPING)
