#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
组织通量概览图（总量 + 均值柱状图）

优先从 01_df_results/all_tissues/overview/ 读取 overview CSV，兼容旧路径
01_df_results/overview/，输出到 02_result_plot/all_tissues/6_tissue_overview/。
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path


def plot_tissue_overview(ov_df: pd.DataFrame, output_dir: Path, condition_colors: dict):
    """Dual panel: total + mean absolute flux per tissue by condition."""
    ov_dir = output_dir / "02_result_plot" / "all_tissues" / "6_tissue_overview"
    ov_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    pivot_total = ov_df.pivot_table(values="total_abs_flux", index="tissue", columns="condition")
    pivot_total.plot(kind="bar", ax=axes[0],
                     color=[condition_colors.get(c, "#999") for c in pivot_total.columns],
                     alpha=0.8, edgecolor="white")
    axes[0].set_title("Total Absolute Flux per Tissue", fontweight="bold")
    axes[0].set_ylabel("Sum |Flux|")
    axes[0].tick_params(axis="x", rotation=30)
    axes[0].legend(title="Condition")
    axes[0].grid(axis="y", alpha=0.3)

    pivot_mean = ov_df.pivot_table(values="mean_abs_flux", index="tissue", columns="condition")
    pivot_mean.plot(kind="bar", ax=axes[1],
                    color=[condition_colors.get(c, "#999") for c in pivot_mean.columns],
                    alpha=0.8, edgecolor="white")
    axes[1].set_title("Mean Absolute Flux per Reaction per Tissue", fontweight="bold")
    axes[1].set_ylabel("Mean |Flux| / reaction")
    axes[1].tick_params(axis="x", rotation=30)
    axes[1].legend(title="Condition")
    axes[1].grid(axis="y", alpha=0.3)

    plt.tight_layout()
    plt.savefig(ov_dir / "tissue_flux_overview.png", dpi=300, bbox_inches="tight")
    plt.close()


def generate_overview_plots(output_dir: Path, condition_colors: dict):
    """从已保存的 overview CSV 生成组织概览图."""
    csv_path = output_dir / "01_df_results" / "all_tissues" / "overview" / "tissue_flux_overview.csv"
    if not csv_path.exists():
        csv_path = output_dir / "01_df_results" / "overview" / "tissue_flux_overview.csv"
    if not csv_path.exists():
        print("[SKIP overview] 01_df_results/all_tissues/overview/tissue_flux_overview.csv 不存在")
        return

    print("\n  Generating tissue overview plot...")
    ov_df = pd.read_csv(csv_path)
    plot_tissue_overview(ov_df, output_dir, condition_colors)
    print("    all_tissues/6_tissue_overview/tissue_flux_overview.png")
