"""
Enrichment analysis subpackage for proteomics.

Provides sample-level ssGSEA and comparison-level pre-ranked GSEA.
"""
from .gsea import (
    DEFAULT_GENE_SETS,
    build_ranked_list,
    find_prerank_report,
    find_ssgsea_report,
    run_prerank_gsea,
    run_ssgsea,
    run_ssgsea_from_config,
)

__all__ = [
    "DEFAULT_GENE_SETS",
    "build_ranked_list",
    "find_prerank_report",
    "find_ssgsea_report",
    "run_prerank_gsea",
    "run_ssgsea",
    "run_ssgsea_from_config",
]
