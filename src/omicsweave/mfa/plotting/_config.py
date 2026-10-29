#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MFA 绘图模块共享配置

存储所有绘图子模块共用的状态（输出目录、配色、条件列表等）。
"""

from __future__ import annotations

from pathlib import Path

from omicsweave.visualization.theme import apply_theme

# 统一主题（MFA 图件在 savefig 处显式指定 dpi）
apply_theme()

# ── 模块级共享状态（由 apply_config 填充） ──────────────────────────

OUTPUT_DIR: Path | None = None
CONDITION_COLORS: dict = {}
TISSUE_COLORS: dict = {}
CONDITIONS: list = []
SAMPLE_CONDITIONS: dict = {}
PATHWAY_GROUPS: dict = {}


def apply_config(
    output_dir: Path | str,
    condition_colors: dict | None = None,
    tissue_colors: dict| None = None,
    conditions: list| None = None,
    sample_conditions: dict| None = None,
    pathway_groups: dict| None = None,
):
    """填充模块级共享变量，在绘图前调用一次."""
    global OUTPUT_DIR, CONDITION_COLORS, TISSUE_COLORS
    global CONDITIONS, SAMPLE_CONDITIONS, PATHWAY_GROUPS

    OUTPUT_DIR = Path(output_dir)
    CONDITION_COLORS = condition_colors or {}
    TISSUE_COLORS = tissue_colors or {}
    CONDITIONS = conditions or []
    SAMPLE_CONDITIONS = sample_conditions or {}
    PATHWAY_GROUPS = pathway_groups or {}


# ── 目录辅助函数 ─────────────────────────────────────────────────────────

_SUBDIR_MAP = {
    "qc": "02_result_plot/all_tissues/0_qc",
    "volcano": "02_result_plot/1_volcano",
    "heatmap": "02_result_plot/2_heatmap",
    "pca": "02_result_plot/all_tissues/3_pca_clustering",
    "pathway": "02_result_plot/4_pathway_flux",
    "tca": "02_result_plot/all_tissues/5_tca_cycle",
    "overview": "02_result_plot/all_tissues/6_tissue_overview",
    "mid": "02_result_plot/all_tissues/7_mid_analysis",
}


def plot_dir(category: str) -> Path:
    """返回指定类别图表的输出目录，并自动创建.

    Args:
        category: "qc" | "volcano" | "heatmap" | "pca" | "pathway" |
                  "tca" | "overview" | "mid"
    """
    if OUTPUT_DIR is None:
        raise RuntimeError("plotting._config 未初始化，请先调用 apply_config()")
    subdir = _SUBDIR_MAP.get(category, category)
    d = OUTPUT_DIR / subdir
    d.mkdir(parents=True, exist_ok=True)
    return d


def df_dir(subdir: str) -> Path:
    """返回 01_df_results/ 下指定子目录的数据路径（只读，不创建）.

    Args:
        subdir: "differential_flux" | "pca" | "heatmap" | "pathway_flux" |
                "tca_cycle" | "qc" | "overview" | "mid"
    """
    if OUTPUT_DIR is None:
        raise RuntimeError("plotting._config 未初始化，请先调用 apply_config()")
    return OUTPUT_DIR / "01_df_results" / subdir
