#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MFA 绘图模块

从已保存的分析表格生成所有可视化图表。支持独立运行（不需要重新计算）。

Usage:
    from omicsweave.mfa.plotting import generate_all_plots
    generate_all_plots(output_dir, condition_colors, tissue_colors, ...)
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .volcano import generate_volcano_plots
from .heatmap import generate_heatmap_plots
from .pca import generate_pca_plots
from .pathway import generate_pathway_plots
from .overview import generate_overview_plots
from .mid import generate_mid_plots
from .qc import generate_qc_plots

# 默认调色板：当配置文件中的颜色不足时自动分配
_PALETTE = [
    "#FFA500", "#008000", "#BA55D3", "#FF4500", "#1E90FF",
    "#FF6347", "#20B2AA", "#DAA520", "#8B008B", "#2E8B57",
    "#CD5C5C", "#4682B4", "#D2691E", "#6A5ACD", "#BDB76B",
]

# 默认组织颜色（覆盖常见组织名称变体）
_DEFAULT_TISSUE_COLORS = {
    "serum": "#FF4500", "heart": "#A83C09", "liver": "#8B008B",
    "kidney": "#663399", "brain": "#FFDAB9", "muscle": "#FA8072",
    "white_adipose": "#FFEFD5",
    "gastroc_muscle": "#FA8072", "sub_q_fat": "#FFEFD5",
    "sq_adipose": "#FFEFD5",
}


def _resolve_colors(
    all_keys_in_data: list[str],
    color_map: dict[str, str],
    defaults: dict[str, str] | None = None,
) -> dict[str, str]:
    """确保每个键都有颜色；先用 defaults 补全，仍缺失时从调色板自动分配。"""
    resolved = dict(color_map)
    # 用默认映射补全
    if defaults:
        for k, v in defaults.items():
            if k not in resolved:
                resolved[k] = v
    # 仍有缺失则从调色板分配
    missing = [c for c in sorted(set(all_keys_in_data)) if c not in resolved]
    if not missing:
        return resolved
    used = set(resolved.values())
    available = [c for c in _PALETTE if c not in used] or _PALETTE
    for i, key in enumerate(missing):
        resolved[key] = available[i % len(available)]
    return resolved


def _infer_conditions(output_dir: Path) -> list[str]:
    """从 PCA scores CSV 推断条件列表。"""
    csv = output_dir / "01_df_results" / "all_tissues" / "pca" / "flux_pca_scores.csv"
    if not csv.exists():
        csv = output_dir / "01_df_results" / "pca" / "flux_pca_scores.csv"
    if not csv.exists():
        return []
    df = pd.read_csv(csv)
    if "condition" not in df.columns:
        return []
    return sorted(df["condition"].unique().tolist())


def generate_all_plots(
    output_dir: Path | str,
    condition_colors: dict | None = None,
    tissue_colors: dict | None = None,
    comparisons: list | None = None,
    tissues: list | None = None,
    conditions: list | None = None,
    sample_conditions: dict | None = None,
    top_metabolites: int = 20,
    analysis_config: dict | None = None,
):
    """
    从已保存的分析表格生成所有 MFA 可视化图表。

    Parameters
    ----------
    output_dir : Path or str
        MFA 分析输出根目录（含 01_df_results/ 子目录）。
    condition_colors : dict, optional
        条件 → 颜色映射。
    tissue_colors : dict, optional
        组织 → 颜色映射。
    comparisons : list of tuple, optional
        比较对列表，如 [("hfd", "lfd")]。为 None 时自动从表格推断。
    tissues : list, optional
        组织列表。为 None 时从表格推断。
    conditions : list, optional
        条件列表。
    sample_conditions : dict, optional
        样本 → 条件映射（热图着色用）。
    top_metabolites : int
        MID 分析中展示的 Top 代谢物数量。
    analysis_config : dict, optional
        config.yaml 中的 analysis 配置块，用于火山图等阈值。
    """
    output_dir = Path(output_dir)
    condition_colors = condition_colors or {}
    tissue_colors = tissue_colors or {}
    analysis_config = analysis_config or {}
    volcano_cfg = analysis_config.get("volcano", {}) or {}

    print("\n" + "=" * 60)
    print("MFA Plot Generation")
    print("=" * 60)

    # 自动推断 comparisons
    if not comparisons:
        comparisons = _infer_comparisons(output_dir)

    # 自动推断 tissues
    if not tissues:
        tissues = _infer_tissues(output_dir)

    # 自动推断 conditions（从 CSV 数据中提取）
    if not conditions:
        conditions = _infer_conditions(output_dir)
    print(f"  [INFO] conditions: {conditions}")
    print(f"  [INFO] condition_colors keys: {list(condition_colors.keys())}")

    # 解析颜色：确保每个条件和组织都有颜色，缺失时自动分配
    if conditions:
        condition_colors = _resolve_colors(conditions, condition_colors)
    if tissues:
        tissue_colors = _resolve_colors(tissues, tissue_colors, _DEFAULT_TISSUE_COLORS)
    print(f"  [INFO] resolved condition_colors: {condition_colors}")
    print(f"  [INFO] resolved tissue_colors: {tissue_colors}")

    # 1. QC plots
    generate_qc_plots(output_dir, condition_colors)

    # 2. Volcano plots
    if comparisons:
        generate_volcano_plots(
            output_dir,
            condition_colors,
            tissue_colors,
            comparisons,
            fc_thresh=float(volcano_cfg.get("log2fc_threshold", 0.5)),
            p_thresh=float(volcano_cfg.get("pvalue_threshold", 0.05)),
            d_thresh=float(volcano_cfg.get("cohen_d_threshold", 0.5)),
            label_top_n=int(volcano_cfg.get("label_top_n", 10)),
            axis_config=volcano_cfg,
        )

    # 3. Heatmap
    if comparisons and sample_conditions:
        generate_heatmap_plots(output_dir, condition_colors, sample_conditions,
                               tissues, comparisons)

    # 4. PCA
    generate_pca_plots(output_dir, condition_colors, sample_conditions, tissues)

    # 5. Pathway flux boxplots
    generate_pathway_plots(output_dir, condition_colors, conditions, tissues)

    # 6. Tissue overview
    generate_overview_plots(output_dir, condition_colors)

    # 7. MID analysis
    generate_mid_plots(output_dir, condition_colors, top_metabolites)

    print("\n" + "=" * 60)
    print("All plots generated!")
    print(f"Output: {output_dir / '02_result_plot'}")
    print("=" * 60)


def _infer_comparisons(output_dir: Path) -> list:
    """从 differential_flux CSV 文件名推断比较对（递归搜索各组织子目录）."""
    diff_root = output_dir / "01_df_results"
    if not diff_root.exists():
        return []
    comps = set()
    for p in sorted(diff_root.glob("**/differential_flux_*.csv")):
        name = p.stem.replace("differential_flux_", "")
        parts = name.split("_vs_")
        if len(parts) == 2:
            comps.add((parts[0].lower(), parts[1].lower()))
    return sorted(comps)


def _infer_tissues(output_dir: Path) -> list:
    """从 PCA scores CSV 推断组织列表."""
    csv = output_dir / "01_df_results" / "all_tissues" / "pca" / "flux_pca_scores.csv"
    if not csv.exists():
        csv = output_dir / "01_df_results" / "pca" / "flux_pca_scores.csv"
    if not csv.exists():
        return []
    df = pd.read_csv(csv)
    return sorted(df["tissue"].unique().tolist()) if "tissue" in df.columns else []
