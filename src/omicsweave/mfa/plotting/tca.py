#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TCA 循环柱状图

优先从 01_df_results/all_tissues/tca_cycle/ + 01_df_results/<tissue>/tca_cycle/ 读取 CSV，
兼容旧路径 01_df_results/tca_cycle/。
输出到 02_result_plot/<tissue>/5_tca_cycle/（单组织图）
和 02_result_plot/all_tissues/5_tca_cycle/（跨组织总图）。
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path


def plot_tca_cycle(tca_df: pd.DataFrame, output_dir: Path, condition_colors: dict):
    """Grouped bar chart: mean flux per (tissue x TCA reaction) by condition."""
    pivot = tca_df.pivot_table(values="mean", index=["tissue", "reaction"],
                               columns="condition")
    if pivot.empty:
        return

    tca_dir = output_dir / "02_result_plot" / "all_tissues" / "5_tca_cycle"
    tca_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(max(10, len(pivot) * 0.4), 6))
    pivot.plot(kind="bar", ax=ax,
               color=[condition_colors.get(c, "#999") for c in pivot.columns],
               alpha=0.8, edgecolor="white")
    ax.set_ylabel("Flux")
    ax.set_title("TCA Cycle Fluxes Across Tissues", fontweight="bold")
    labels = [f"{t}\n{r}" for t, r in pivot.index]
    ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=7)
    ax.legend(title="Condition")
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(tca_dir / "tca_cycle_fluxes.png", dpi=300, bbox_inches="tight")
    plt.close()


def generate_tca_plots(output_dir: Path, condition_colors: dict):
    """读取 TCA CSV（兼容新旧目录）并生成按组织分组的 TCA 循环图."""
    df_root = output_dir / "01_df_results"
    if not df_root.exists():
        print("[SKIP tca] 01_df_results/ 目录不存在")
        return

    print("\n  Generating TCA cycle plots...")

    # 收集所有组织的 TCA 数据
    all_frames = []
    generated_tissues = set()

    # 新版按组织目录
    for csv_path in sorted(df_root.glob("*/tca_cycle/tca_cycle_fluxes.csv")):
        tissue = csv_path.parent.parent.name
        df = pd.read_csv(csv_path)
        if "tissue" not in df.columns and not df.empty:
            df = df.copy()
            df["tissue"] = tissue
        if df.empty:
            continue
        all_frames.append(df)

        # 每个组织单独生成 TCA 图
        tca_dir = output_dir / "02_result_plot" / tissue / "5_tca_cycle"
        tca_dir.mkdir(parents=True, exist_ok=True)
        plot_tca_single_tissue(df, tissue, tca_dir, condition_colors)
        generated_tissues.add(tissue)

    # 跨组织目录（新）：单文件包含所有组织
    legacy_csv = df_root / "all_tissues" / "tca_cycle" / "tca_cycle_fluxes.csv"
    # 旧版根层级目录兼容
    if not legacy_csv.exists():
        legacy_csv = df_root / "tca_cycle" / "tca_cycle_fluxes.csv"
    if legacy_csv.exists():
        legacy_df = pd.read_csv(legacy_csv)
        if "tissue" not in legacy_df.columns:
            if not legacy_df.empty and not all_frames:
                all_frames.append(legacy_df)
        else:
            for tissue in sorted(legacy_df["tissue"].dropna().astype(str).unique()):
                if tissue in generated_tissues:
                    continue
                tissue_df = legacy_df[legacy_df["tissue"] == tissue].copy()
                if tissue_df.empty:
                    continue
                all_frames.append(tissue_df)
                tca_dir = output_dir / "02_result_plot" / tissue / "5_tca_cycle"
                tca_dir.mkdir(parents=True, exist_ok=True)
                plot_tca_single_tissue(tissue_df, tissue, tca_dir, condition_colors)
                generated_tissues.add(tissue)

    if not all_frames:
        print("  No TCA cycle data found")
        return

    # 合并所有组织数据生成跨组织总图
    tca_df = pd.concat(all_frames, ignore_index=True)
    plot_tca_cycle(tca_df, output_dir, condition_colors)
    print("  all_tissues/5_tca_cycle/tca_cycle_fluxes.png (combined)")


def plot_tca_single_tissue(tca_df: pd.DataFrame, tissue: str, out_dir: Path,
                            condition_colors: dict):
    """单组织 TCA 循环柱状图。"""
    pivot = tca_df.pivot_table(values="mean", index=["reaction"],
                                columns="condition")
    if pivot.empty:
        return

    fig, ax = plt.subplots(figsize=(max(8, len(pivot) * 0.6), 5))
    pivot.plot(kind="bar", ax=ax,
               color=[condition_colors.get(c, "#999") for c in pivot.columns],
               alpha=0.8, edgecolor="white")
    ax.set_ylabel("Flux")
    ax.set_title(f"TCA Cycle — {tissue}", fontweight="bold")
    ax.set_xticklabels(pivot.index, rotation=45, ha="right", fontsize=9)
    ax.legend(title="Condition")
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_dir / "tca_cycle_fluxes.png", dpi=300, bbox_inches="tight")
    plt.close()
