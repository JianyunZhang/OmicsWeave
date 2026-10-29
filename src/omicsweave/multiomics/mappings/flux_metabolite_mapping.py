#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Flux → Metabolite Mapping  (YAML-backed, validated against common/ metabolite registry)

Data is loaded from ``flux_metabolite_mapping.yaml`` in the same directory.
Metabolite names are validated at import time via ``lookup_metabolite_ids()``
from ``omicsweave.common.metabolite_id_mapping`` (fuzzy matching).

Public API (unchanged):
    FLUX_METABOLITE_MAPPING  — dict[str, dict]
    get_related_metabolites()
    get_metabolite_roles()
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


FLUX_METABOLITE_MAPPING: dict[str, Any] = _load_yaml("flux_metabolite_mapping.yaml")

# ---------------------------------------------------------------------------
# Validation against common/metabolite_id_mapping
# ---------------------------------------------------------------------------


def _validate_metabolites(mapping: dict[str, Any]) -> list[str]:
    try:
        from omicsweave.common.metabolite_id_mapping import lookup_metabolite_ids
    except ImportError:
        return ["[skip] omicsweave.common.metabolite_id_mapping not importable"]

    warns: list[str] = []

    for rxn_id, entry in mapping.items():
        ids_by_role: dict[str, set[str]] = {}
        for role in ("substrates", "products", "related"):
            role_ids: set[str] = set()
            for met_name in entry.get(role, []):
                result = lookup_metabolite_ids(met_name)
                if result:
                    role_ids.update(str(value) for value in result.values() if value)
            ids_by_role[role] = role_ids
        roles = ("substrates", "products", "related")
        for index, first in enumerate(roles):
            for second in roles[index + 1:]:
                overlap = ids_by_role[first] & ids_by_role[second]
                if overlap:
                    warns.append(
                        f"{rxn_id}: database IDs overlap across {first}/{second}: "
                        + ", ".join(sorted(overlap))
                    )

    return warns


_met_warnings = _validate_metabolites(FLUX_METABOLITE_MAPPING)

# ---------------------------------------------------------------------------
# Helper functions (preserved from original)
# ---------------------------------------------------------------------------


def get_related_metabolites(reaction_name: str) -> list[str] | None:
    """Get list of all related metabolites for a reaction (backward compatible)."""
    for key, mapping in FLUX_METABOLITE_MAPPING.items():
        if key in reaction_name:
            all_metabs = (
                mapping.get("substrates", [])
                + mapping.get("products", [])
                + mapping.get("related", [])
            )
            return all_metabs
    return None


def get_metabolite_roles(reaction_name: str) -> dict[str, list[str]] | None:
    """Get substrate/product/related metabolites for a reaction.

    Returns
    -------
    dict or None
        {"substrates": [...], "products": [...], "related": [...]} or None.
    """
    for key, mapping in FLUX_METABOLITE_MAPPING.items():
        if key in reaction_name:
            return mapping
    return None
