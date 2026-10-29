"""
Heatmap module for proteomics.

Generates top-N protein heatmaps with hierarchical clustering.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from . import _config as cfg


def plot_single(
    df: pd.DataFrame,
    group_map: dict[str, str],
    diff_results: pd.DataFrame,
    output_path: str,
    *,
    group_colors: dict[str, str] | None = None,
    heatmap_top_n: int | None = None,
) -> None:
    """Generate top-N protein heatmap with hierarchical clustering.

    Args:
        df: Expression matrix.
        group_map: Sample name -> group mapping.
        diff_results: Combined DE results DataFrame.
        output_path: Path to save the heatmap.
        group_colors: Optional group -> color mapping.
        heatmap_top_n: Number of top proteins to display.
    """
    top_n = heatmap_top_n or cfg.HEATMAP_TOP_N
    colors = group_colors or cfg.GROUP_COLORS

    significance_cols = [c for c in diff_results.columns if c.startswith("Padj_")]
    if not significance_cols:
        significance_cols = [c for c in diff_results.columns if c.startswith("Pval_")]
    if not significance_cols:
        return

    diff_results = diff_results.copy()
    diff_results["Min_significance"] = diff_results[significance_cols].min(axis=1)
    finite = diff_results["Min_significance"].notna()
    if not finite.any():
        return
    top_proteins = (
        diff_results.loc[finite]
        .sort_values("Min_significance")
        .head(top_n)["Protein"]
        .tolist()
    )

    available = [p for p in top_proteins if p in df.index]
    if not available:
        return
    subset_df = df.loc[available]

    unique_groups = sorted(set(group_map.values()) - {"Unknown"})
    lut = {g: colors.get(g, "#888888") for g in unique_groups}
    col_colors = subset_df.columns.map(group_map).map(lut)

    g = sns.clustermap(
        subset_df, cmap="RdBu_r", col_colors=col_colors,
        z_score=0, center=0, figsize=(12, 12),
    )
    g.savefig(output_path, dpi=150)
    plt.close()


def generate_all() -> None:
    """Regenerate heatmap images and ranked-protein tables for all tissues."""
    from ..experiment_config import get_sample_to_group_map

    cfg.flush_print("Generating heatmap outputs...")
    for tissue in cfg.TISSUES:
        tissue_dir = os.path.join(cfg.PLOT_BASE, tissue)
        os.makedirs(tissue_dir, exist_ok=True)

        res_df = cfg.load_combined_dp(tissue)
        if res_df is None:
            continue

        significance_cols = [c for c in res_df.columns if c.startswith("Padj_")]
        if not significance_cols:
            significance_cols = [c for c in res_df.columns if c.startswith("Pval_")]
        if not significance_cols:
            continue

        res_df = res_df.copy()
        res_df["Min_significance"] = res_df[significance_cols].min(axis=1)
        top_proteins = (
            res_df[res_df["Min_significance"].notna()]
            .sort_values("Min_significance")
            .head(cfg.HEATMAP_TOP_N)
        )

        top_path = os.path.join(tissue_dir, "top_proteins.csv")
        top_proteins[["Protein"] + significance_cols].to_csv(top_path, index=False)

        matrix_path = os.path.join(cfg.DP_BASE, tissue, "imputed_matrix.csv")
        if not os.path.exists(matrix_path):
            cfg.flush_print(f"  [{tissue}] No imputed matrix found, skipping heatmap image.")
            continue
        matrix = pd.read_csv(matrix_path, index_col=0)
        group_map = get_sample_to_group_map(matrix.columns, cfg.MOUSE_GROUPS)
        output_path = os.path.join(tissue_dir, f"heatmap_top{cfg.HEATMAP_TOP_N}.png")
        plot_single(
            matrix,
            group_map,
            res_df,
            output_path,
            group_colors=cfg.GROUP_COLORS,
            heatmap_top_n=cfg.HEATMAP_TOP_N,
        )

    cfg.flush_print("  Heatmap outputs done.")
