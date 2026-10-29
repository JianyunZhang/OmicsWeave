"""
GO/KEGG Pathway Enrichment Analysis (Enrichr)
===============================================
Self-contained subpackage for ORA enrichment analysis.

Provides:
  - run_enrichment(): public API (callable from other scripts)
"""

from ._run import run_enrichment
