#!/usr/bin/env python3
"""
Scoring visualization: score distributions, confidence pie charts.

Adapted from flux_multiomics_integration.py lines 515-576 (concordance analysis).

Author: Claude Code
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from . import _config as cfg


def plot_score_distribution(
    scored_df: pd.DataFrame,
    output_dir: Path | None = None,
) -> list[Path]:
    """Histogram and boxplot of total multi-omics scores.

    Parameters
    ----------
    scored_df : DataFrame
        Output from FluxOmicsScorer.score_all() with columns
        total_score, confidence, comparison, tissue.
    output_dir : Path, optional

    Returns
    -------
    list[Path]
        Saved file paths.
    """
    if scored_df.empty:
        return []

    out = output_dir or cfg.PLOTS_DIR
    out.mkdir(parents=True, exist_ok=True)

    confidence_colors = cfg.PLOT_COLORS.get("confidence", {
        "High": "#1B7837",
        "Medium": "#762A83",
        "Low": "#F4A582",
        "Insufficient": "#D9D9D9",
    })
    saved: list[Path] = []

    # 1. Overall score histogram with confidence bands
    fig, ax = plt.subplots(figsize=(10, 6))
    scores = scored_df["total_score"].values
    ax.hist(scores, bins=30, color="#4C72B0", alpha=0.7, edgecolor="white")

    # Add confidence threshold lines
    thresholds = cfg.PLOT_COLORS  # not ideal but use scoring thresholds
    for label, color in confidence_colors.items():
        pass  # thresholds are in config.scoring, not easily accessible here

    ax.axvline(3.0, color="#1B7837", linestyle="--", alpha=0.7, label="High (3.0)")
    ax.axvline(2.0, color="#762A83", linestyle="--", alpha=0.7, label="Medium (2.0)")
    ax.axvline(1.0, color="#F4A582", linestyle="--", alpha=0.7, label="Low (1.0)")
    ax.set_xlabel("Total Multi-Omics Score")
    ax.set_ylabel("Count")
    ax.set_title("Distribution of Multi-Omics Integration Scores", fontweight="bold")
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    p = out / "score_distribution.png"
    fig.savefig(p, dpi=300)
    plt.close(fig)
    saved.append(p)

    # 2. Confidence level pie chart
    fig2, ax2 = plt.subplots(figsize=(7, 7))
    counts = scored_df["confidence"].value_counts()
    colors = [confidence_colors.get(c, "#999999") for c in counts.index]
    wedges, texts, autotexts = ax2.pie(
        counts.values, labels=counts.index, colors=colors,
        autopct="%1.1f%%", startangle=90, textprops={"fontsize": 11},
    )
    ax2.set_title("Confidence Level Distribution", fontweight="bold")
    fig2.savefig(out / "confidence_pie.png", dpi=300)
    plt.close(fig2)
    saved.append(out / "confidence_pie.png")

    # 3. Score by tissue (boxplot)
    if "tissue" in scored_df.columns:
        tissues = sorted(scored_df["tissue"].unique())
        if len(tissues) > 1:
            fig3, ax3 = plt.subplots(figsize=(max(8, len(tissues) * 1.5), 6))
            data_by_tissue = [
                scored_df.loc[scored_df["tissue"] == t, "total_score"].values
                for t in tissues
            ]
            bp = ax3.boxplot(data_by_tissue, labels=tissues, patch_artist=True)
            tissue_colors = cfg.PLOT_COLORS.get("tissue", {})
            for patch, tissue in zip(bp["boxes"], tissues):
                patch.set_facecolor(tissue_colors.get(tissue, "#4C72B0"))
                patch.set_alpha(0.6)
            ax3.set_ylabel("Total Score")
            ax3.set_title("Score Distribution by Tissue", fontweight="bold")
            ax3.grid(axis="y", alpha=0.3)
            plt.xticks(rotation=30, ha="right")
            plt.tight_layout()
            fig3.savefig(out / "score_by_tissue.png", dpi=300)
            plt.close(fig3)
            saved.append(out / "score_by_tissue.png")

    # 4. Score by comparison (boxplot)
    if "comparison" in scored_df.columns:
        comparisons = sorted(scored_df["comparison"].unique())
        if len(comparisons) > 1:
            fig4, ax4 = plt.subplots(figsize=(max(8, len(comparisons) * 2), 6))
            data_by_comp = [
                scored_df.loc[scored_df["comparison"] == c, "total_score"].values
                for c in comparisons
            ]
            labels = [cfg.COMP_LABELS.get(c, c) for c in comparisons]
            ax4.boxplot(data_by_comp, labels=labels)
            ax4.set_ylabel("Total Score")
            ax4.set_title("Score Distribution by Comparison", fontweight="bold")
            ax4.grid(axis="y", alpha=0.3)
            plt.xticks(rotation=30, ha="right")
            plt.tight_layout()
            fig4.savefig(out / "score_by_comparison.png", dpi=300)
            plt.close(fig4)
            saved.append(out / "score_by_comparison.png")

    print(f"  [scoring] Saved {len(saved)} scoring plots")
    return saved