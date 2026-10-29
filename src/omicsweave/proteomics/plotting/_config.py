"""
Plotting configuration module.

Populates module-level globals from the experiment config,
providing shared state for all plotting modules.
"""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from ..experiment_config import (
    ProteomicsConfig,
    load_proteomics_config,
    generate_comparisons,
    get_group_colors,
    get_tissue_labels,
    get_comp_labels,
)

# Set plotting style
try:
    plt.style.use("seaborn-whitegrid")
except OSError:
    plt.style.use("seaborn-v0_8-whitegrid")
plt.rcParams["font.sans-serif"] = ["SimHei", "Microsoft YaHei", "Arial"]
plt.rcParams["axes.unicode_minus"] = False

# ============================================================================
# Module-level globals (populated by apply_config)
# ============================================================================
OUTPUT_BASE = ""
DP_BASE = ""         # 01_dp_results
PLOT_BASE = ""       # 02_result_plot
TISSUES: list[str] = []
COMPARISONS: list[tuple[str, str, str]] = []
MOUSE_GROUPS: dict[str, str] = {}
GROUP_COLORS: dict[str, str] = {}
TISSUE_LABELS: dict[str, str] = {}
COMP_LABELS: dict[str, str] = {}

PVAL_CUTOFF = 0.05
LOG2FC_CUTOFF = 0.58
HEATMAP_TOP_N = 50
DP_METHOD = "limma"

GSEA_MIN_SET_SIZE = 10
GSEA_MAX_SET_SIZE = 500
GSEA_PERMUTATION_NUM = 1000
SSGSEA_PERMUTATION_NUM = 0
GSEA_SEED = 42
GSEA_THREADS = 1
RUN_SSGSEA = True
RUN_PRERANK_GSEA = True
GENE_SETS: list[str] = ["KEGG_2019_Mouse"]


def flush_print(msg: str) -> None:
    print(msg, flush=True)


def apply_config(
    experiment: str | None,
    config_path: str | None,
    outdir: str | None,
) -> dict:
    """Load config and populate plotting module globals."""
    global OUTPUT_BASE, DP_BASE, PLOT_BASE, TISSUES, COMPARISONS
    global MOUSE_GROUPS, GROUP_COLORS, TISSUE_LABELS, COMP_LABELS
    global PVAL_CUTOFF, LOG2FC_CUTOFF, HEATMAP_TOP_N, DP_METHOD
    global GSEA_MIN_SET_SIZE, GSEA_MAX_SET_SIZE, GSEA_PERMUTATION_NUM
    global SSGSEA_PERMUTATION_NUM, GSEA_SEED, GSEA_THREADS
    global RUN_SSGSEA, RUN_PRERANK_GSEA, GENE_SETS

    cfg = load_proteomics_config(
        experiment_id=experiment,
        config_path=config_path,
        outdir_override=outdir,
    )

    OUTPUT_BASE = str(cfg.output_dir)
    DP_BASE = os.path.join(OUTPUT_BASE, "01_dp_results")
    PLOT_BASE = os.path.join(OUTPUT_BASE, "02_result_plot")
    TISSUES = list(cfg.tissues)
    COMPARISONS = generate_comparisons(cfg.conditions, cfg.contrast_groups)
    MOUSE_GROUPS = dict(cfg.mouse_groups)
    GROUP_COLORS = get_group_colors(cfg.conditions)
    TISSUE_LABELS = get_tissue_labels(cfg)
    COMP_LABELS = get_comp_labels(cfg.conditions, COMPARISONS)

    PVAL_CUTOFF = cfg.pval_cutoff
    LOG2FC_CUTOFF = cfg.log2fc_cutoff
    HEATMAP_TOP_N = cfg.heatmap_top_n
    DP_METHOD = cfg.dp_method

    GSEA_MIN_SET_SIZE = cfg.gsea_min_set_size
    GSEA_MAX_SET_SIZE = cfg.gsea_max_set_size
    GSEA_PERMUTATION_NUM = cfg.gsea_permutation_num
    SSGSEA_PERMUTATION_NUM = cfg.ssgsea_permutation_num
    GSEA_SEED = cfg.gsea_seed
    GSEA_THREADS = cfg.gsea_threads
    RUN_SSGSEA = cfg.run_ssgsea
    RUN_PRERANK_GSEA = cfg.run_prerank_gsea
    GENE_SETS = list(cfg.gene_sets)

    return {
        "experiment_id": cfg.experiment_id,
        "output_dir": OUTPUT_BASE,
        "_cfg": cfg,  # expose for enrichment module
    }


def load_dp_results(tissue: str, comp_name: str) -> pd.DataFrame | None:
    """Load per-comparison DP CSV from 01_dp_results/<tissue>/."""
    path = os.path.join(DP_BASE, tissue, f"dp_{tissue}_{comp_name}.csv")
    if os.path.exists(path):
        return pd.read_csv(path)
    return None


def load_combined_dp(tissue: str) -> pd.DataFrame | None:
    """Load combined DP CSV from 01_dp_results/<tissue>/."""
    path = os.path.join(DP_BASE, tissue, "differential_proteins.csv")
    if os.path.exists(path):
        return pd.read_csv(path)
    # Backward compat: check old filename
    old_path = os.path.join(DP_BASE, tissue, "differential_expression_wilcoxon.csv")
    if os.path.exists(old_path):
        return pd.read_csv(old_path)
    return None
