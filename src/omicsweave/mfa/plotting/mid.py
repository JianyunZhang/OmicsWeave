#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MID 标记分析图（总标记度 + MID 分布）

优先从 01_df_results/all_tissues/mid/ 读取 MID CSV，兼容旧路径
01_df_results/mid/，输出到 02_result_plot/all_tissues/7_mid_analysis/。
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path


def plot_mid_labeling(total_label: pd.DataFrame, metabolites: list,
                       output_dir: Path, condition_colors: dict):
    """Bar chart of total labeling for top metabolites."""
    plot_df = total_label[total_label["metabolite"].isin(metabolites)]
    if plot_df.empty:
        return

    mid_dir = output_dir / "02_result_plot" / "all_tissues" / "7_mid_analysis"
    mid_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(max(10, len(metabolites) * 0.8), 6))
    metabs = sorted(metabolites)
    x = np.arange(len(metabs))
    conds = sorted(plot_df["condition"].unique())
    w = 0.8 / len(conds)

    for ci, cond in enumerate(conds):
        vals = []
        for m in metabs:
            sub = plot_df[(plot_df["metabolite"] == m) & (plot_df["condition"] == cond)]
            vals.append(sub["total_labeling"].values[0] if len(sub) > 0 else 0)
        color = condition_colors.get(cond, "#999")
        ax.bar(x + ci * w, vals, w, label=cond.upper(), color=color, alpha=0.8)

    ax.set_xticks(x + w * (len(conds) - 1) / 2)
    ax.set_xticklabels(metabs, rotation=60, ha="right", fontsize=7)
    ax.set_ylabel("Total Labeling (sum M+1, M+2, ...)")
    ax.set_title("Top Variable Metabolite Labeling by Condition", fontweight="bold")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(mid_dir / "mid_total_labeling.png", dpi=300, bbox_inches="tight")
    plt.close()


def plot_mid_distributions(mid_summary: pd.DataFrame, metabolites: list,
                            output_dir: Path, condition_colors: dict):
    """Grid of grouped bar charts: MID distribution per metabolite."""
    n = len(metabolites)
    cols = min(4, n)
    rows_n = (n + cols - 1) // cols

    mid_dir = output_dir / "02_result_plot" / "all_tissues" / "7_mid_analysis"
    mid_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(rows_n, cols, figsize=(cols * 4, rows_n * 3.5), squeeze=False)

    for idx, metab in enumerate(metabolites):
        ax = axes[idx // cols, idx % cols]
        m_data = mid_summary[mid_summary["metabolite"] == metab]
        if m_data.empty:
            ax.set_visible(False)
            continue

        conds = sorted(m_data["condition"].unique())
        mid_indices = sorted(m_data["mid_index"].unique())
        x = np.arange(len(mid_indices))
        w = 0.8 / len(conds)

        for ci, cond in enumerate(conds):
            vals = [m_data[(m_data["condition"] == cond) & (m_data["mid_index"] == mi)]["mean"].values
                    for mi in mid_indices]
            vals = [v[0] if len(v) > 0 else 0 for v in vals]
            color = condition_colors.get(cond, "#999")
            ax.bar(x + ci * w, vals, w, label=cond.upper(), color=color, alpha=0.8)

        ax.set_xticks(x + w * (len(conds) - 1) / 2)
        ax.set_xticklabels([f"M+{int(i)}" for i in mid_indices], fontsize=7)
        short_name = metab if len(metab) < 25 else metab[:22] + "..."
        ax.set_title(short_name, fontsize=9, fontweight="bold")
        ax.set_ylabel("Fraction")
        if idx == 0:
            ax.legend(fontsize=7)
        ax.grid(axis="y", alpha=0.3)

    for idx in range(n, rows_n * cols):
        axes[idx // cols, idx % cols].set_visible(False)

    plt.suptitle("Mass Isotopomer Distributions", fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(mid_dir / "mid_distributions.png", dpi=300, bbox_inches="tight")
    plt.close()


def generate_mid_plots(output_dir: Path, condition_colors: dict, top_metabolites: int = 20):
    """从已保存的 MID CSV 生成 MID 图."""
    mid_dir = output_dir / "01_df_results" / "all_tissues" / "mid"
    if not mid_dir.exists():
        mid_dir = output_dir / "01_df_results" / "mid"
    total_path = mid_dir / "mid_total_labeling.csv"
    summary_path = mid_dir / "mid_summary.csv"
    if not total_path.exists() or not summary_path.exists():
        print("[SKIP mid] 01_df_results/all_tissues/mid/ CSV 文件不存在")
        return

    print("\n  Generating MID plots...")
    total_label = pd.read_csv(total_path)
    mid_summary = pd.read_csv(summary_path)

    # Find top variable metabolites
    var_per_metab = total_label.groupby("metabolite")["total_labeling"].var()
    top_metabs = var_per_metab.nlargest(top_metabolites).index.tolist()

    if top_metabs:
        plot_mid_labeling(total_label, top_metabs, output_dir, condition_colors)

    unique_metabs = total_label.sort_values("total_labeling", ascending=False)["metabolite"].unique()[:8]
    if len(unique_metabs) > 0:
        plot_mid_distributions(mid_summary, list(unique_metabs), output_dir, condition_colors)

    print("    all_tissues/7_mid_analysis/mid_total_labeling.png")
    print("    all_tissues/7_mid_analysis/mid_distributions.png")
