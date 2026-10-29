"""
Shared configuration and utilities for all plotting modules.
"""

import os
from typing import Any

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')

from omicsweave.visualization.theme import apply_theme
from ..experiment_config import (
    load_transcriptomics_config, config_arg_help, resolve_comparisons,
    get_group_colors, get_tissue_labels, get_comp_labels, get_plot_colors,
)

# ============================================================================
# Module-level state — populated by apply_config()
# ============================================================================
OUTPUT_BASE: str = ""
DE_DIR: str = ""
PLOT_DIR: str = ""
TISSUES: list[str] = []
COMPARISONS: list[str] = []
COMP_TUPLES: list[tuple[str, str, str]] = []
MOUSE_GROUPS: dict[str, str] = {}
GROUP_COLORS: dict[str, str] = {}
TISSUE_LABELS: dict[str, str] = {}
COMP_LABELS: dict[str, str] = {}
PLOT_COLORS: dict[str, Any] = {}

PADJ_CUTOFF: float = 0.05
LOG2FC_CUTOFF: float = 1.0


def apply_config(experiment=None, config_path=None, outdir=None) -> dict:
    """Load experiment config and populate all module-level variables."""
    global OUTPUT_BASE, DE_DIR, PLOT_DIR, TISSUES, COMPARISONS, COMP_TUPLES
    global MOUSE_GROUPS, GROUP_COLORS, TISSUE_LABELS, COMP_LABELS, PLOT_COLORS

    cfg = load_transcriptomics_config(
        experiment_id=experiment,
        config_path=config_path,
        outdir_override=outdir,
    )

    OUTPUT_BASE = os.path.abspath(str(cfg.output_dir))
    os.environ["TRANSCRIPTOMICS_OUTDIR"] = OUTPUT_BASE
    DE_DIR = os.path.join(OUTPUT_BASE, "04_de_results")
    PLOT_DIR = os.path.join(OUTPUT_BASE, "05_result_plot")
    TISSUES = list(cfg.tissues)
    MOUSE_GROUPS = dict(cfg.mouse_groups)

    comp_tuples = resolve_comparisons(cfg)
    COMPARISONS = [name for _, _, name in comp_tuples]
    COMP_TUPLES = list(comp_tuples)

    GROUP_COLORS = get_group_colors(cfg.conditions)
    TISSUE_LABELS = get_tissue_labels(cfg)
    COMP_LABELS = get_comp_labels(cfg.conditions, comp_tuples)
    PLOT_COLORS = get_plot_colors(cfg)

    return {
        "experiment_id": cfg.experiment_id,
        "config_path": str(cfg.config_path),
        "output_dir": OUTPUT_BASE,
    }


def flush_print(*args, **kwargs):
    print(*args, **kwargs, flush=True)


def load_sample_metadata() -> pd.DataFrame | None:
    """
    读取 pipeline 落盘的 sample_metadata.csv（含 condition / tissue / batch）。

    批次列只有在 groups.yaml 配置了 sample_batches 时才存在，因此调用方
    必须容忍其缺失。
    """
    path = os.path.join(OUTPUT_BASE, "sample_metadata.csv")
    if not os.path.exists(path):
        return None
    try:
        return pd.read_csv(path, index_col=0)
    except Exception:
        return None


def get_sample_batches() -> dict[str, str]:
    """返回 {sample_id: batch}；未配置批次时返回空 dict。"""
    meta = load_sample_metadata()
    if meta is None or "batch" not in meta.columns:
        return {}
    batches = meta["batch"].dropna()
    return {str(k): str(v) for k, v in batches.items()}


def load_de_results(tissue: str, comparison: str) -> pd.DataFrame | None:
    """Load a DE result CSV file (PyDESeq2 or Wilcoxon fallback)."""
    for suffix in ("", "_wilcoxon"):
        path = os.path.join(DE_DIR, tissue, f"de_{comparison}{suffix}.csv")
        if os.path.exists(path):
            df = pd.read_csv(path, index_col=0)
            # Annotated files have numeric index; restore gene symbol index
            for col in ('gene_symbol', 'gene'):
                if col in df.columns and not df.index.inferred_type == 'string':
                    df = df.set_index(col)
                    break
            return df
    return None


# 统一主题；转录组图件沿用 200 dpi
apply_theme(dpi=200, save_dpi=200)
