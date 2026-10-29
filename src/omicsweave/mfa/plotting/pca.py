#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PCA 散点图 + Scree plot

优先从 01_df_results/all_tissues/pca/ 读取 PCA CSV，兼容旧路径
01_df_results/pca/，输出到 02_result_plot/all_tissues/3_pca_clustering/。
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path


def plot_pca(scores_df: pd.DataFrame, tissue: str, out_dir: Path,
             condition_colors: dict, sample_conditions: dict,
             pca_data_dir: Path):
    """PCA scatter plot (PC1 vs PC2) + variance scree bar chart."""
    # 按 tissue 过滤
    if "tissue" not in scores_df.columns:
        return
    t_scores = scores_df[scores_df["tissue"] == tissue].copy()
    if t_scores.empty:
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Scatter plot — 使用 CSV 中的 condition 列
    cond_col = "condition" if "condition" in t_scores.columns else None
    if cond_col:
        for cond in sorted(t_scores[cond_col].unique()):
            mask = t_scores[cond_col] == cond
            color = condition_colors.get(cond, "#999")
            axes[0].scatter(t_scores.loc[mask, "PC1"], t_scores.loc[mask, "PC2"],
                            c=color, label=cond.upper(), s=80, alpha=0.8, edgecolors="white")

    axes[0].set_xlabel("PC1")
    axes[0].set_ylabel("PC2")
    axes[0].set_title(f"PCA - {tissue}", fontweight="bold")
    axes[0].legend()
    axes[0].grid(alpha=0.3)

    # Scree plot — 从 loadings CSV 读取 variance_explained 行
    loadings_path = pca_data_dir / f"flux_pca_loadings_{tissue}.csv"
    if loadings_path.exists():
        loadings_df = pd.read_csv(loadings_path, index_col=0)
        if "variance_explained" in loadings_df.index:
            var_vals = pd.to_numeric(
                loadings_df.loc["variance_explained"], errors="coerce"
            ).dropna().values
            n_pc = min(5, len(var_vals))
            var_vals = var_vals[:n_pc]
            axes[1].bar(range(1, len(var_vals) + 1), var_vals * 100,
                        color="#4C72B0", alpha=0.8)
            axes[1].set_xlabel("Principal Component")
            axes[1].set_ylabel("Variance Explained (%)")
            axes[1].set_title("Scree Plot", fontweight="bold")
            axes[1].set_xticks(range(1, len(var_vals) + 1))
            axes[1].grid(axis="y", alpha=0.3)

    plt.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_dir / f"flux_pca_{tissue}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def generate_pca_plots(output_dir: Path, condition_colors: dict,
                        sample_conditions: dict, tissues: list):
    """从已保存的 PCA scores CSV 生成所有 PCA 图."""
    df_dir = output_dir / "01_df_results" / "all_tissues" / "pca"
    if not df_dir.exists():
        df_dir = output_dir / "01_df_results" / "pca"
    scores_path = df_dir / "flux_pca_scores.csv"
    if not scores_path.exists():
        print("[SKIP pca] 01_df_results/all_tissues/pca/flux_pca_scores.csv 不存在")
        return

    print("\n  Generating PCA plots...")
    scores_df = pd.read_csv(scores_path, index_col=0)
    p_dir = output_dir / "02_result_plot" / "all_tissues" / "3_pca_clustering"

    for tissue in tissues:
        plot_pca(scores_df, tissue, p_dir, condition_colors, sample_conditions, df_dir)
        print(f"    all_tissues/3_pca_clustering/flux_pca_{tissue}.png")
