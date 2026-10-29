#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
通路箱线图：
1) 组织级总览图（全部 pathway 分面 + 全部反应分图）
2) 按 pathway 分类输出（每个 pathway 两类图）

从 01_df_results/pathway_flux/ 读取 pathway_flux_samples CSV，输出到
02_result_plot/<tissue>/4_pathway_flux/。
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from pathlib import Path
import re


def _sanitize_pathway_name(pathway: str) -> str:
    """将 pathway 名称转为稳定的目录/文件名片段。"""
    text = str(pathway).strip().lower()
    text = re.sub(r"[^0-9a-zA-Z_]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text or "unknown_pathway"


def _build_palette(t_df: pd.DataFrame, condition_colors: dict, conditions: list) -> dict:
    """根据数据中实际条件构建颜色映射。"""
    conds_in_data = sorted(t_df["condition"].dropna().astype(str).unique()) if "condition" in t_df.columns else []
    conds = conditions if conditions else conds_in_data
    if not conds:
        conds = conds_in_data
    return {cond: condition_colors.get(cond, "#999") for cond in conds}


def _shorten_label(text: str, max_len: int = 24) -> str:
    """缩短过长反应名，降低 x 轴标签拥挤程度。"""
    txt = str(text)
    if len(txt) <= max_len:
        return txt
    return txt[: max_len - 1] + "…"


def plot_pathway_boxplot(t_df: pd.DataFrame, tissue: str, out_dir: Path,
                          condition_colors: dict, conditions: list,
                          output_filename: str | None = None,
                          title: str | None = None):
    """多面板通路箱线图：按 pathway 分面，x=reaction，hue=condition。"""
    pathways = sorted(t_df["pathway"].unique())
    n_pw = len(pathways)
    if n_pw == 0:
        return

    palette = _build_palette(t_df, condition_colors, conditions)
    # 通路较多时使用双列分面，减少单图横向拉伸
    n_cols = 1 if n_pw <= 2 else 2
    n_rows = int(np.ceil(n_pw / n_cols))
    max_reactions = int(
        max((t_df[t_df["pathway"] == pw]["reaction"].nunique() for pw in pathways), default=1)
    )
    subplot_w = min(14.0, max(8.0, max_reactions * 0.65))
    subplot_h = 4.2
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(subplot_w * n_cols, subplot_h * n_rows),
        squeeze=False,
    )
    axes_flat = axes.flatten()
    legend_handles = None
    legend_labels = None

    for idx, pw in enumerate(pathways):
        ax = axes_flat[idx]
        pw_data = t_df[t_df["pathway"] == pw]
        if pw_data.empty:
            continue

        reactions = sorted(pw_data["reaction"].dropna().astype(str).unique())
        conds_in_pw = sorted(pw_data["condition"].dropna().astype(str).unique()) if "condition" in pw_data.columns else []
        hue_order = [c for c in conditions if c in conds_in_pw] if conditions else conds_in_pw
        if not hue_order:
            hue_order = conds_in_pw

        sns.boxplot(
            data=pw_data,
            x="reaction",
            y="flux",
            hue="condition",
            order=reactions,
            hue_order=hue_order,
            palette=palette,
            width=0.75,
            showfliers=True,
            fliersize=2,
            linewidth=1.0,
            ax=ax,
        )

        tick_fontsize = 8 if len(reactions) <= 12 else 7
        short_labels = [_shorten_label(r, max_len=20) for r in reactions]
        # 显式固定刻度位置后再设置标签，避免 matplotlib 警告
        ticks = ax.get_xticks()
        if len(ticks) == len(short_labels):
            ax.set_xticks(ticks)
            ax.set_xticklabels(short_labels, rotation=35, ha="right", fontsize=tick_fontsize)
        else:
            ax.tick_params(axis="x", labelrotation=35, labelsize=tick_fontsize)
        ax.set_ylabel("Flux (IQR-filtered)")
        ax.set_xlabel("Reaction")
        ax.set_title(f"{pw.replace('_', ' ')}", fontweight="bold")
        ax.grid(axis="y", alpha=0.3)
        if ax.get_legend() is not None:
            if legend_handles is None or legend_labels is None:
                legend_handles, legend_labels = ax.get_legend_handles_labels()
            ax.get_legend().remove()

    # 隐藏空白子图
    for idx in range(n_pw, len(axes_flat)):
        axes_flat[idx].set_visible(False)

    if legend_handles and legend_labels:
        fig.legend(
            legend_handles,
            legend_labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.995),
            ncol=min(len(legend_labels), 4),
            frameon=False,
            title="Conditions",
        )

    fig_title = title or f"Pathway Flux Overview - {tissue}"
    fig.suptitle(fig_title, fontweight="bold", fontsize=14, y=1.02)
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    out_dir.mkdir(parents=True, exist_ok=True)
    save_name = output_filename or f"pathway_flux_{tissue}.png"
    fig.savefig(out_dir / save_name, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_all_reactions_boxplot(t_df: pd.DataFrame, tissue: str, out_dir: Path,
                                condition_colors: dict, conditions: list,
                                output_filename: str | None = None,
                                title: str | None = None):
    """Grid of boxplots: one subplot per reaction, samples individually colored."""
    reactions = sorted(t_df["reaction"].unique())
    if not reactions:
        return

    all_samples_info = t_df[["sample_id", "condition"]].drop_duplicates().sort_values("condition")
    sample_order = all_samples_info["sample_id"].tolist()
    palette = _build_palette(t_df, condition_colors, conditions)

    n_reactions = len(reactions)
    n_cols = min(4, n_reactions)
    n_rows = (n_reactions + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 5, n_rows * 4), squeeze=False)

    for idx, rxn in enumerate(reactions):
        ax = axes[idx // n_cols, idx % n_cols]
        rxn_data = t_df[t_df["reaction"] == rxn]

        box_data = []
        colors = []
        labels = []

        for sid in sample_order:
            sample_vals = rxn_data[rxn_data["sample_id"] == sid]["flux"].values
            if len(sample_vals) > 0:
                box_data.append(sample_vals)
                cond = all_samples_info[all_samples_info["sample_id"] == sid]["condition"].values[0]
                colors.append(palette.get(cond, "#999"))
                labels.append(sid)

        if not box_data:
            ax.set_visible(False)
            continue

        bp = ax.boxplot(box_data, widths=0.6, patch_artist=True, showfliers=True)

        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        for median in bp['medians']:
            median.set_color('black')

        ax.set_xticklabels(labels, rotation=90, fontsize=8)
        ax.set_title(rxn, fontsize=10, fontweight="bold")
        ax.set_ylabel("Flux")
        ax.grid(axis="y", alpha=0.3)

    for idx in range(n_reactions, n_rows * n_cols):
        axes[idx // n_cols, idx % n_cols].set_visible(False)

    fig_title = title or f"All Reactions: Individual Sample Distributions - {tissue}"
    plt.suptitle(fig_title, fontweight="bold", fontsize=16, y=1.02)
    plt.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    save_name = output_filename or f"pathway_flux_{tissue}_individual.png"
    plt.savefig(out_dir / save_name, dpi=300, bbox_inches="tight")
    plt.close()


def generate_by_pathway_plots(t_df: pd.DataFrame, tissue: str, out_dir: Path,
                              condition_colors: dict, conditions: list):
    """按 pathway 分类输出两类图（总反应箱线图 + 单反应拆分箱线图）。"""
    by_pathway_dir = out_dir / "by_pathway"
    pathways = sorted(t_df["pathway"].dropna().astype(str).unique())
    for pathway in pathways:
        pw_df = t_df[t_df["pathway"] == pathway].copy()
        if pw_df.empty:
            continue
        safe_name = _sanitize_pathway_name(pathway)
        pw_out_dir = by_pathway_dir / safe_name

        plot_pathway_boxplot(
            pw_df,
            tissue,
            pw_out_dir,
            condition_colors,
            conditions,
            output_filename=f"pathway_flux_{tissue}_{safe_name}.png",
            title=f"{pathway} - {tissue} (Pathway-level Overview)",
        )
        plot_all_reactions_boxplot(
            pw_df,
            tissue,
            pw_out_dir,
            condition_colors,
            conditions,
            output_filename=f"pathway_flux_{tissue}_{safe_name}_individual.png",
            title=f"{pathway} - {tissue} (Reaction-level Detail)",
        )
        print(f"    {tissue}/4_pathway_flux/by_pathway/{safe_name}/")


def generate_pathway_plots(output_dir: Path, condition_colors: dict,
                            conditions: list, tissues: list):
    """读取 pathway_flux_samples CSV（兼容新旧目录）并生成通路箱线图."""
    df_root = output_dir / "01_df_results"
    generated_tissues = set()

    # 新版按组织目录
    for tissue in tissues:
        csv_path = df_root / tissue / "pathway_flux" / "pathway_flux_samples.csv"
        if not csv_path.exists():
            continue
        sample_df = pd.read_csv(csv_path)
        t_df = sample_df[sample_df["tissue"] == tissue] if "tissue" in sample_df.columns else sample_df
        if t_df.empty:
            continue
        pw_dir = output_dir / "02_result_plot" / tissue / "4_pathway_flux"
        plot_pathway_boxplot(t_df, tissue, pw_dir, condition_colors, conditions)
        plot_all_reactions_boxplot(t_df, tissue, pw_dir, condition_colors, conditions)
        generate_by_pathway_plots(t_df, tissue, pw_dir, condition_colors, conditions)
        print(f"    {tissue}/4_pathway_flux/pathway_flux_{tissue}.png")
        generated_tissues.add(tissue)

    # 旧版根层级目录：单文件包含所有组织
    legacy_csv = df_root / "pathway_flux" / "pathway_flux_samples.csv"
    if not legacy_csv.exists():
        return

    sample_df = pd.read_csv(legacy_csv)
    if "tissue" not in sample_df.columns:
        print("    [SKIP pathway] 根层级 pathway_flux_samples.csv 缺少 tissue 列")
        return

    for tissue in sorted(sample_df["tissue"].dropna().astype(str).unique()):
        if tissue in generated_tissues:
            continue
        t_df = sample_df[sample_df["tissue"] == tissue].copy()
        if t_df.empty:
            continue
        pw_dir = output_dir / "02_result_plot" / tissue / "4_pathway_flux"
        plot_pathway_boxplot(t_df, tissue, pw_dir, condition_colors, conditions)
        plot_all_reactions_boxplot(t_df, tissue, pw_dir, condition_colors, conditions)
        generate_by_pathway_plots(t_df, tissue, pw_dir, condition_colors, conditions)
        print(f"    {tissue}/4_pathway_flux/pathway_flux_{tissue}.png")
