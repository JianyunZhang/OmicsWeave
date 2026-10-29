#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compatibility re-export module.

The canonical reaction registry now lives at:
    omicsweave.common.reaction_pathway_mapping

This module re-exports everything from there so existing import paths
(e.g. ``from omicsweave.multiomics.mappings.reaction_registry import REACTION_REGISTRY``)
continue to work without modification.
"""

from pathlib import Path


from omicsweave.common.reaction_pathway_mapping import (
    REACTION_REGISTRY,
    ReactionInfo,
    get_reaction_info,
    get_kegg_reaction,
    get_kegg_pathway,
    get_reactome_id,
    get_ec_number,
    get_description,
    list_reactions_by_pathway,
    list_reactions_with_ec,
    get_all_pathways,
    annotate_flux_dataframe,
)

__all__ = [
    "REACTION_REGISTRY",
    "ReactionInfo",
    "get_reaction_info",
    "get_kegg_reaction",
    "get_kegg_pathway",
    "get_reactome_id",
    "get_ec_number",
    "get_description",
    "list_reactions_by_pathway",
    "list_reactions_with_ec",
    "get_all_pathways",
    "annotate_flux_dataframe",
]
