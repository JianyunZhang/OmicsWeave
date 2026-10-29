"""
Pathway visualization module for proteomics.

Regenerates ssGSEA pathway enrichment via the enrichment subpackage,
using the saved imputed expression matrix. Supports multiple gene set libraries.
"""
from __future__ import annotations

import os

import pandas as pd

from . import _config as cfg


def generate_all() -> None:
    """Regenerate ssGSEA pathway analysis for all tissues.

    Requires imputed_matrix.csv saved by the pipeline in 01_dp_results/<tissue>/.
    """
    from ..enrichment.gsea import run_prerank_gsea, run_ssgsea

    cfg.flush_print("Running ssGSEA pathway analysis...")
    for tissue in cfg.TISSUES:
        imp_path = os.path.join(cfg.DP_BASE, tissue, "imputed_matrix.csv")
        if not os.path.exists(imp_path):
            cfg.flush_print(f"  [{tissue}] No imputed matrix found, skipping ssGSEA.")
            continue

        tissue_dir = os.path.join(cfg.PLOT_BASE, tissue)
        os.makedirs(tissue_dir, exist_ok=True)

        imp_df = pd.read_csv(imp_path, index_col=0)
        if cfg.RUN_SSGSEA:
            run_ssgsea(
                imp_df, tissue_dir,
                gene_sets=cfg.GENE_SETS,
                min_size=cfg.GSEA_MIN_SET_SIZE,
                max_size=cfg.GSEA_MAX_SET_SIZE,
                permutation_num=cfg.SSGSEA_PERMUTATION_NUM,
                seed=cfg.GSEA_SEED,
                threads=cfg.GSEA_THREADS,
                flush=cfg.flush_print,
            )

        if cfg.RUN_PRERANK_GSEA:
            results = cfg.load_combined_dp(tissue)
            if results is None:
                cfg.flush_print(
                    f"  [{tissue}] No differential results found, skipping pre-ranked GSEA."
                )
            else:
                run_prerank_gsea(
                    results,
                    tissue_dir,
                    cfg.COMPARISONS,
                    gene_sets=cfg.GENE_SETS,
                    min_size=cfg.GSEA_MIN_SET_SIZE,
                    max_size=cfg.GSEA_MAX_SET_SIZE,
                    permutation_num=cfg.GSEA_PERMUTATION_NUM,
                    seed=cfg.GSEA_SEED,
                    threads=cfg.GSEA_THREADS,
                    flush=cfg.flush_print,
                )

    cfg.flush_print("  Pathway analysis done.")
