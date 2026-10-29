#!/usr/bin/env python3
"""
Shared plotting configuration and utilities for multi-omics visualization.

Populated by ``apply_config()`` with experiment-specific settings (colors,
labels, output directories). All plotting modules read from this module's
module-level variables.

Pattern adapted from: omicsweave/transcriptomics/plotting/_config.py

Author: Claude Code
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")

from omicsweave.visualization.theme import apply_theme

from ..config import (
    MultiOmicsConfig,
    load_multiomics_config,
    generate_comparisons,
    get_group_colors,
    get_tissue_labels,
    get_comp_labels,
    get_plot_colors,
)

# ============================================================================
# Module-level state — populated by apply_config()
# ============================================================================

OUTPUT_DIR: Path = Path(".")
COMMON_DIR: Path = Path(".")   # merged tables + auxiliary CSVs
PLOTS_DIR: Path = Path(".")    # 02_result_plot

TISSUES: list[str] = []
COMPARISONS: list[str] = []
MOUSE_GROUPS: dict[str, str] = {}
GROUP_COLORS: dict[str, str] = {}
TISSUE_LABELS: dict[str, str] = {}
COMP_LABELS: dict[str, str] = {}
PLOT_COLORS: dict[str, Any] = {}

PVAL_CUTOFF: float = 0.05
LOG2FC_CUTOFF: float = 0.58


def apply_config(
    experiment: str | None = None,
    config_path: str | None = None,
    outdir: str | None = None,
) -> dict[str, str]:
    """Load experiment config and populate all module-level variables.

    Parameters
    ----------
    experiment : str, optional
        Experiment ID.
    config_path : str, optional
        Direct path to config.yaml.
    outdir : str, optional
        Override output directory.

    Returns
    -------
    dict
        Summary of loaded configuration.
    """
    global OUTPUT_DIR, COMMON_DIR, PLOTS_DIR
    global TISSUES, COMPARISONS, MOUSE_GROUPS
    global GROUP_COLORS, TISSUE_LABELS, COMP_LABELS, PLOT_COLORS
    global PVAL_CUTOFF, LOG2FC_CUTOFF

    cfg: MultiOmicsConfig = load_multiomics_config(
        experiment_id=experiment,
        config_path=config_path,
        outdir_override=outdir,
    )

    OUTPUT_DIR = cfg.output_dir
    COMMON_DIR = cfg.output_dir / "common"
    PLOTS_DIR = cfg.output_dir / "02_result_plot"

    TISSUES = list(cfg.tissues)
    COMPARISONS = list(cfg.comparisons)
    MOUSE_GROUPS = dict(cfg.mouse_groups)

    comp_tuples = generate_comparisons(cfg.conditions)
    GROUP_COLORS = get_group_colors(cfg.conditions)
    TISSUE_LABELS = get_tissue_labels(cfg)
    COMP_LABELS = get_comp_labels(cfg.conditions, comp_tuples)
    PLOT_COLORS = get_plot_colors(cfg)

    scoring = cfg.scoring
    PVAL_CUTOFF = scoring.get("pval_cutoff", 0.05)
    LOG2FC_CUTOFF = scoring.get("log2fc_cutoff", 0.58)

    return {
        "experiment_id": cfg.experiment_id,
        "config_path": str(cfg.config_path),
        "output_dir": str(cfg.output_dir),
    }


# ============================================================================
# Matplotlib global settings
# ============================================================================

apply_theme(dpi=150, save_dpi=300)


def flush_print(*args: Any, **kwargs: Any) -> None:
    """Print with flush=True for progress updates."""
    print(*args, **kwargs, flush=True)