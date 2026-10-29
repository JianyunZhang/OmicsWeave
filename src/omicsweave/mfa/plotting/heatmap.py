#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
差异通量聚类热图

从 01_df_results/heatmap/ 读取 heatmap CSV，输出到 02_result_plot/2_heatmap/。
支持 N<2 低生物学重复模式：自动读取 meta 信息并在图上标注。
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from pathlib import Path


def _extract_tissue_from_heatmap_filename(csv_path: Path, comp_name: str) -> str | None:
    """从旧版根层级 heatmap 文件名提取 tissue。"""
    prefix = "flux_heatmap_raw_"
    suffix = f"_{comp_name}.csv"
    name = csv_path.name
    if name.startswith(prefix) and name.endswith(suffix):
        return name[len(prefix):-len(suffix)]
    return None


def _iter_heatmap_inputs(df_root: Path, tissues: list, comp_name: str):
    """
    迭代热图输入，兼容两种目录布局：
      1) 新版: 01_df_results/<tissue>/heatmap/flux_heatmap_raw_{tissue}_{comp}.csv
      2) 旧版: 01_df_results/heatmap/flux_heatmap_raw_{tissue}_{comp}.csv
    """
    yielded_tissues = set()

    # 新版按组织目录
    for tissue in tissues:
        csv_dir = df_root / tissue / "heatmap"
        csv_path = csv_dir / f"flux_heatmap_raw_{tissue}_{comp_name}.csv"
        if not csv_path.exists():
            continue
        low_n = False
        meta_path = csv_dir / f"heatmap_meta_{comp_name}.csv"
        if meta_path.exists():
            meta = pd.read_csv(meta_path)
            if "effect_size_mode" in meta.columns:
                low_n = bool(meta["effect_size_mode"].iloc[0])
            elif "low_n" in meta.columns:
                low_n = bool(meta["low_n"].iloc[0])
        yielded_tissues.add(tissue)
        yield tissue, csv_path, low_n

    # 旧版根层级目录
    legacy_dir = df_root / "heatmap"
    if not legacy_dir.exists():
        return

    legacy_low_n = False
    legacy_meta_path = legacy_dir / f"heatmap_meta_{comp_name}.csv"
    if legacy_meta_path.exists():
        legacy_meta = pd.read_csv(legacy_meta_path)
        if "effect_size_mode" in legacy_meta.columns:
            legacy_low_n = bool(legacy_meta["effect_size_mode"].iloc[0])
        elif "low_n" in legacy_meta.columns:
            legacy_low_n = bool(legacy_meta["low_n"].iloc[0])

    for csv_path in sorted(legacy_dir.glob(f"flux_heatmap_raw_*_{comp_name}.csv")):
        tissue = _extract_tissue_from_heatmap_filename(csv_path, comp_name)
        if tissue is None or tissue in yielded_tissues:
            continue
        yield tissue, csv_path, legacy_low_n


def _ensure_condition_colors(unique_conds, condition_colors):
    """Ensure every condition has a color; auto-assign from a palette if missing."""
    _PALETTE = [
        "#FFA500", "#008000", "#BA55D3", "#FF4500", "#1E90FF",
        "#FF6347", "#20B2AA", "#DAA520", "#8B008B", "#2E8B57",
        "#CD5C5C", "#4682B4", "#D2691E", "#6A5ACD", "#BDB76B",
    ]
    resolved = dict(condition_colors)
    missing = [c for c in unique_conds if c not in resolved]
    if missing:
        # pick colors not already used
        used = set(resolved.values())
        available = [c for c in _PALETTE if c not in used]
        for i, cond in enumerate(missing):
            resolved[cond] = available[i % len(available)] if available else _PALETTE[i % len(_PALETTE)]
    return resolved


def plot_flux_heatmap(heatmap_z: pd.DataFrame, sample_conditions: dict,
                       condition_colors: dict, comp_name: str, tissue: str,
                       out_dir: Path, low_n: bool = False):
    """Seaborn clustermap for differential flux heatmap."""
    all_sids = list(heatmap_z.columns)

    unique_conds = sorted(set(sample_conditions.get(sid, "unknown") for sid in all_sids))
    resolved_colors = _ensure_condition_colors(unique_conds, condition_colors)

    col_colors = [resolved_colors.get(sample_conditions.get(sid, "unknown"), "#999")
                  for sid in all_sids]

    legend_handles = [plt.Rectangle((0, 0), 1, 1,
                                    facecolor=resolved_colors.get(c, "#999"),
                                    label=c.upper())
                      for c in unique_conds]

    n_rows = len(heatmap_z)
    n_cols = len(all_sids)

    # 宽度：为右侧图例和颜色标尺留足空间
    fig_w = max(14, n_cols * 1.8 + 12)
    fig_h = max(10, n_rows * 0.35)

    g = sns.clustermap(
        heatmap_z,
        col_colors=col_colors,
        cmap="RdBu_r",
        center=0,
        figsize=(fig_w, fig_h),
        xticklabels=True,
        yticklabels=True,
        dendrogram_ratio=(0.12, 0.06),
        cbar_pos=(0.95, 0.35, 0.012, 0.30),
        row_cluster=True,
        col_cluster=False,
    )

    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), fontsize=8, rotation=45, ha="right")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=7)
    g.ax_heatmap.set_xlabel("Samples (Individuals)")
    g.ax_heatmap.set_ylabel("Reaction")

    # 图例放在颜色标尺下方，两者均在热图区域右侧
    g.fig.legend(handles=legend_handles, title="Condition",
                 loc='center left', bbox_to_anchor=(0.955, 0.15),
                 fontsize=9, frameon=True, fancybox=True)

    title = f"Differential Flux Heatmap - {tissue}\n{comp_name}"
    subtitle = "\n(Z-score per reaction, IQR-filtered data)"
    if low_n:
        subtitle += "\n⚠ Effect-size mode: top reactions ranked by |Hedges' g| (biological replicates)"
    plt.suptitle(title + subtitle, fontweight="bold", fontsize=12, y=1.02)

    # 统计方法图注
    if low_n:
        stat_note_text = (
            "Statistical Note: Heatmap shows top reactions ranked by |Hedges' g| effect size "
            "(no significance threshold applied). The effect size is computed from "
            "BIOLOGICAL REPLICATES — each animal contributes one net-flux value (median over "
            "its IQR-filtered model fits); fitting iterations are never treated as replicates. "
            "The p-value path is omitted because the rank test cannot reach the configured "
            "threshold at this sample size. Exploratory — validation with additional "
            "biological replicates is required."
        )
    else:
        stat_note_text = (
            "Statistical Note: Heatmap shows reactions with Mann-Whitney U test P<0.05 on "
            "biological replicates, ranked by |log2FC| (top 50). "
            "Row values are Z-score standardized; each cell is one animal's net flux "
            "(median over its IQR-filtered model fits)."
        )
    g.fig.text(0.5, -0.02, stat_note_text, ha="center", va="top",
               fontsize=6.5, style="italic", wrap=True,
               bbox=dict(boxstyle="round,pad=0.3", facecolor="#f5f5f5", alpha=0.8))
    out_dir.mkdir(parents=True, exist_ok=True)
    g.savefig(out_dir / f"flux_heatmap_{tissue}_{comp_name}.png", dpi=300, bbox_inches="tight")
    plt.close()


def generate_heatmap_plots(output_dir: Path, condition_colors: dict,
                            sample_conditions: dict, tissues: list,
                            comparisons: list):
    """读取 heatmap CSV（兼容新旧目录）并生成按组织分组的热力图."""
    df_root = output_dir / "01_df_results"
    if not df_root.exists():
        print("[SKIP heatmap] 01_df_results/ 目录不存在")
        return

    print("\n  Generating heatmap plots...")
    for cond_a, cond_b in comparisons:
        comp_name = f"{cond_a.upper()}_vs_{cond_b.upper()}"
        for tissue, csv_path, low_n in _iter_heatmap_inputs(df_root, tissues, comp_name):
            raw_df = pd.read_csv(csv_path, index_col=0)
            h_dir = output_dir / "02_result_plot" / tissue / "2_heatmap"
            plot_flux_heatmap(raw_df, sample_conditions, condition_colors,
                               comp_name, tissue, h_dir, low_n=low_n)
            print(f"    {tissue}/2_heatmap/flux_heatmap_{tissue}_{comp_name}.png")
