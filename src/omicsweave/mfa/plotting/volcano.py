#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
火山图 + 组织差异通量汇总图

从 01_df_results/<tissue>/differential_flux/ 读取 CSV，
输出到 02_result_plot/<tissue>/1_volcano/。

支持两种模式：
  - N>=2（biological_replicates）: Y 轴 = -log10(p-value)，P-value 筛选
  - N<2（fitting_variance_only） : Y 轴 = |Cohen's d|，效应量筛选
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

from omicsweave.visualization.palettes import METRO_ATLAS

from . import _config as cfg


def _iter_differential_flux_tables(df_root: Path, comp_name: str):
    """
    迭代指定比较对的差异通量表，兼容两种目录布局：
      1) 新版: 01_df_results/<tissue>/differential_flux/differential_flux_{comp}.csv
      2) 旧版: 01_df_results/differential_flux/differential_flux_{comp}.csv（含 tissue 列）
    """
    pattern = f"**/differential_flux/differential_flux_{comp_name}.csv"
    for csv_path in sorted(df_root.glob(pattern)):
        # 旧版根层级布局：01_df_results/differential_flux/*.csv
        if csv_path.parent.parent == df_root:
            root_df = pd.read_csv(csv_path)
            if root_df.empty:
                continue
            if "tissue" not in root_df.columns:
                print(f"    [SKIP] 根层级 differential_flux 缺少 tissue 列: {csv_path.name}")
                continue
            for tissue, tissue_df in root_df.groupby("tissue", sort=True):
                if tissue_df.empty:
                    continue
                yield str(tissue), tissue_df.copy()
            continue

        # 新版按组织布局
        tissue = csv_path.parent.parent.name
        tissue_df = pd.read_csv(csv_path)
        if tissue_df.empty:
            continue
        yield tissue, tissue_df


def _is_effect_size_mode(df: pd.DataFrame) -> bool:
    """判断本次比较是否走效应量路径（而非 p 值路径）。

    优先读 ``primary_statistic``：它由 df_analysis 直接写出，既覆盖
    N<2 的情形，也覆盖"样本量足够但秩检验最小可达 p 值高于阈值"
    的情形（例如 n=3 vs 3 时 p_min=0.1）。其余为历史布局回退。
    """
    if "primary_statistic" in df.columns:
        return str(df["primary_statistic"].iloc[0]) == "effect_size_only"
    if "stat_note" in df.columns:
        return df["stat_note"].iloc[0] == "fitting_variance_only"
    n_bio_cols = [c for c in df.columns if c.startswith("n_bio_")]
    for col in n_bio_cols:
        if df[col].iloc[0] < 2:
            return True
    return False


def _effect_size_column(df: pd.DataFrame) -> str:
    """效应量列名：优先小样本校正的 Hedges' g，回退 Cohen's d。"""
    if "hedges_g" in df.columns and df["hedges_g"].notna().any():
        return "hedges_g"
    return "cohen_d"


def _axis_note_text() -> str:
    """坐标轴约定说明，避免读者把对称范围误读为数据范围。"""
    return ("Axes: x is forced symmetric about 0 so up- and down-regulation are "
            "visually comparable; its range is the 98th percentile of |log2FC|, "
            "capped at 10 (1024-fold). Net fluxes whose denominator approaches zero "
            "reach |log2FC| above 20, which is a numerical artefact rather than a "
            "stronger biological change. y starts at 0. Open triangles at the "
            "x-limits are off-scale or infinite fold changes, pinned to the edge "
            "rather than dropped.")


def _power_floor_text(df: pd.DataFrame) -> str:
    """秩检验样本量下界的说明片段（无该列时返回空串）。"""
    if "min_achievable_pvalue" not in df.columns:
        return ""
    floor = df["min_achievable_pvalue"].iloc[0]
    if pd.isna(floor):
        return ""
    return (f"Rank-test minimum achievable two-sided p at this sample size "
            f"is {float(floor):.3g} (set by combinatorics, not by the data). ")


def _axis_limits(
    finite_log2fc,
    y_values,
    fc_thresh: float,
    axis_config: dict,
) -> tuple[tuple[float, float], tuple[float, float], float]:
    """计算火山图坐标范围。

    x 轴**默认强制对称于 0**：火山图的全部意义在于比较上调与下调，
    自适应范围会把「一侧数据碰巧更极端」渲染成「一侧变化更多」的视觉
    印象。对称之后，左右两半的面积可直接比较。

    y 轴自 0 起：``-log10(p)`` 的 0 对应 p=1，``|Hedges' g|`` 的 0 对应无
    效应，都是天然基线，浮动基线会夸大微小差异。

    Returns
    -------
    (xlim, ylim, x_edge)
        ``x_edge`` 为无穷大点的绘制位置（x 轴边界）。
    """
    import numpy as _np

    symmetric = bool(axis_config.get("symmetric_x", True))
    pad = float(axis_config.get("axis_pad", 0.08))

    explicit_x = axis_config.get("xlim")
    if explicit_x and len(explicit_x) == 2:
        xlim = (float(explicit_x[0]), float(explicit_x[1]))
    else:
        finite = _np.asarray(finite_log2fc, dtype=float)
        finite = finite[_np.isfinite(finite)]
        # 净通量的分母接近零时 log2FC 会冲到 ±30 以上——那是数值假象，
        # 不是更强的生物学变化。用分位数而不是最大值定范围，否则两三个
        # 这样的点会把全部有意义的变化压成中轴附近的一条窄带。
        # 越界点不会丢失：与 ±inf 一样钉在边界上并用空心三角标出。
        # 经验分界：全库 |log2FC| 的 98 分位为 9.2（约 567 倍），99 分位骤升到
        # 21.7（约 340 万倍）——后者只可能来自接近零的分母。因此在分位数之上
        # 再加一个绝对上限：代谢通量出现 >1000 倍的变化即属数值假象。
        quantile = float(axis_config.get("x_clip_quantile", 0.98))
        max_abs = float(axis_config.get("x_max_abs", 10.0))
        if finite.size:
            span = float(_np.quantile(_np.abs(finite), min(max(quantile, 0.5), 1.0)))
        else:
            span = 0.0
        # 阈值线必须在可视范围内，否则图上看不到判定边界
        span = max(span, abs(fc_thresh) * 1.5)
        span = min(span, max_abs)
        span = span * (1.0 + pad) if span > 0 else 1.0
        xlim = (-span, span) if symmetric else (
            float(_np.min(finite)) - pad * span if finite.size else -span,
            float(_np.max(finite)) + pad * span if finite.size else span,
        )

    explicit_y = axis_config.get("ylim")
    if explicit_y and len(explicit_y) == 2:
        ylim = (float(explicit_y[0]), float(explicit_y[1]))
    else:
        y_arr = _np.asarray(y_values, dtype=float)
        y_arr = y_arr[_np.isfinite(y_arr)]
        y_top = float(_np.max(y_arr)) if y_arr.size else 1.0
        ylim = (0.0, y_top * (1.0 + pad) if y_top > 0 else 1.0)

    x_edge = xlim[1] * (1.0 - pad / 2)
    return xlim, ylim, x_edge


def plot_volcano(df: pd.DataFrame, comp_name: str, tissue: str, out_dir: Path,
                 condition_colors: dict | None = None,
                 fc_thresh: float = 0.5, p_thresh: float = 0.05,
                 d_thresh: float = 0.5, label_top_n: int = 10,
                 axis_config: dict | None = None):
    """
    单组织差异通量火山图。

    坐标轴：x 对称于 0、y 自 0 起（见 ``_axis_limits``）。
    ``log2FC = ±inf`` 的通量（一组通量约为零、另一组不为零）不再被丢弃，
    而是钉在 x 轴边界上并用空心三角标出——这类"有/无"式变化往往是最值得
    看的，静默剔除会让它们从图上消失。
    """
    axis_config = axis_config or {}
    low_n = _is_effect_size_mode(df)
    effect_col = _effect_size_column(df)
    y_col = effect_col if low_n else "pvalue_mannwhitney"

    plot_df = df.dropna(subset=["log2FC", y_col]).copy()
    if plot_df.empty:
        return

    if low_n:
        plot_df["y"] = plot_df[effect_col].abs()
        y_label = "|Hedges' g|" if effect_col == "hedges_g" else "|Cohen's d|"
        y_thresh = d_thresh
    else:
        plot_df["y"] = -np.log10(plot_df["pvalue_mannwhitney"].clip(lower=1e-300))
        y_label = "-log10(p-value)"
        y_thresh = -np.log10(p_thresh)

    plot_df = plot_df[plot_df["y"].notna()]
    if plot_df.empty:
        return

    is_inf = np.isinf(plot_df["log2FC"])
    finite_df = plot_df[~is_inf]

    xlim, ylim, x_edge = _axis_limits(
        finite_df["log2FC"].to_numpy(), plot_df["y"].to_numpy(), fc_thresh, axis_config
    )

    # 超出坐标范围的点与 ±inf 点同等处理：都属于"倍数变化无界"，
    # 钉在边界上而不是丢弃或让它们撑开坐标轴。
    if bool(axis_config.get("show_infinite", True)):
        out_of_range = is_inf | (plot_df["log2FC"].abs() > xlim[1])
    else:
        out_of_range = pd.Series(False, index=plot_df.index)
    edge_df = plot_df[out_of_range]
    finite_df = plot_df[~out_of_range]

    if len(edge_df):
        edge_df = edge_df.copy()
        edge_df["x_plot"] = np.where(edge_df["log2FC"] > 0, x_edge, -x_edge)
        edge_df["y"] = edge_df["y"].clip(upper=ylim[1])

    def _split(frame):
        if frame.empty:
            return frame, frame, frame
        if "is_significant" in frame.columns:
            mask = frame["is_significant"].fillna(False).astype(bool)
        elif low_n:
            mask = frame["y"].ge(y_thresh) & frame["log2FC"].abs().ge(fc_thresh)
        else:
            mask = (frame["pvalue_mannwhitney"] < p_thresh) & frame["log2FC"].abs().ge(fc_thresh)
        return (frame[~mask],
                frame[mask & (frame["log2FC"] > 0)],
                frame[mask & (frame["log2FC"] < 0)])

    non_sig, sig_up, sig_dn = _split(finite_df)
    edge_ns, edge_up, edge_dn = _split(edge_df)

    fig, ax = plt.subplots(figsize=(10, 8))

    ax.scatter(non_sig["log2FC"], non_sig["y"], c=METRO_ATLAS.faint, s=20, alpha=0.55, label="NS")
    ax.scatter(sig_up["log2FC"], sig_up["y"], c=METRO_ATLAS.up, s=30, alpha=0.82,
               label=f"Up ({len(sig_up) + len(edge_up)})")
    ax.scatter(sig_dn["log2FC"], sig_dn["y"], c=METRO_ATLAS.down, s=30, alpha=0.82,
               label=f"Down ({len(sig_dn) + len(edge_dn)})")

    n_edge = len(edge_df)
    if n_edge:
        for frame, color in ((edge_ns, METRO_ATLAS.faint), (edge_up, METRO_ATLAS.up), (edge_dn, METRO_ATLAS.down)):
            if len(frame):
                marker = np.where(frame["log2FC"] > 0, ">", "<")
                for shape in (">", "<"):
                    picked = frame[marker == shape]
                    if len(picked):
                        ax.scatter(picked["x_plot"], picked["y"], s=55,
                                   facecolors="none", edgecolors=color,
                                   linewidths=1.3, marker=shape, zorder=4)
        ax.scatter([], [], s=55, facecolors="none", edgecolors="#555",
                   marker=">", label=f"off-scale |log2FC| ({n_edge})")

    ax.axhline(y_thresh, ls="--", color="grey", alpha=0.5)
    ax.axvline(fc_thresh, ls="--", color="grey", alpha=0.5)
    ax.axvline(-fc_thresh, ls="--", color="grey", alpha=0.5)
    ax.axvline(0, color="#444", linewidth=0.8, alpha=0.6)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)

    top = pd.concat([sig_up.nlargest(label_top_n, "y"),
                     sig_dn.nlargest(label_top_n, "y")])
    for _, r in top.iterrows():
        label = str(r.get("reaction", r.get("flux_id", "")))
        if len(label) > 25:
            label = label[:25]
        ax.annotate(label, (r["log2FC"], r["y"]),
                    fontsize=6, alpha=0.8,
                    textcoords="offset points", xytext=(5, 5))

    ax.set_xlabel("log2(Fold Change)")
    ax.set_ylabel(y_label)

    title = f"Differential Flux Volcano: {comp_name} — {tissue}"
    if low_n:
        n_bio_cols = [c for c in df.columns if c.startswith("n_bio_")]
        n_parts = [f"{c.replace('n_bio_', '')}={int(df[c].iloc[0])}" for c in n_bio_cols]
        title += (f"\n(Effect-size mode — biological replicates: "
                  f"{', '.join(n_parts)})")
    ax.set_title(title, fontweight="bold", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.2)

    if low_n:
        stat_note_text = (
            f"Statistical Note: Y-axis uses |{y_label.strip('|')}| computed from "
            "BIOLOGICAL REPLICATES — each animal contributes one net-flux value "
            "(median over its IQR-filtered model fits); model-fitting iterations are "
            "never treated as replicates. "
            + _power_floor_text(df) +
            "No p-value threshold is applied on this axis. Results are exploratory "
            "and require validation with additional biological replicates. "
            + _axis_note_text()
        )
    else:
        stat_note_text = (
            f"Statistical Note: Y-axis uses -log10(P-value) from Mann-Whitney U test on "
            f"biological replicates. |log2FC|>{fc_thresh} as fold-change cutoff. "
            "FDR correction (Benjamini-Hochberg) applied. "
            + _power_floor_text(df) + _axis_note_text()
        )
    fig.text(0.5, -0.01, stat_note_text, ha="center", va="top",
             fontsize=7, style="italic", wrap=True,
             bbox=dict(boxstyle="round,pad=0.3", facecolor="#f5f5f5", alpha=0.8))

    plt.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_dir / f"volcano_{comp_name}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_tissue_diff_summary(df: pd.DataFrame, comp_name: str, tissue: str, out_dir: Path,
                              tissue_colors: dict | None = None,
                              p_thresh: float = 0.05, d_thresh: float = 0.5):
    """按条件统计差异通量数量的柱状图（单组织）。"""
    low_n = _is_effect_size_mode(df)
    effect_col = _effect_size_column(df)
    if "is_significant" in df.columns:
        sig = df[df["is_significant"].fillna(False).astype(bool)].copy()
    elif low_n:
        sig = df[(df[effect_col].notna()) & (df[effect_col].abs() >= d_thresh)].copy()
    else:
        sig = df[df["pvalue_mannwhitney"] < p_thresh].copy()

    if sig.empty:
        return

    summary = []
    up = (sig["diff"] > 0).sum()
    down = (sig["diff"] < 0).sum()
    summary.append({"tissue": tissue, "up": up, "down": down})

    sdf = pd.DataFrame(summary)
    colors = tissue_colors or cfg.TISSUE_COLORS

    fig, ax = plt.subplots(figsize=(6, 6))
    color = colors.get(tissue, "#999")
    ax.bar([0], [sdf["up"].iloc[0]], 0.6, label="Increased", color=color, alpha=0.8)
    ax.bar([0], [-sdf["down"].iloc[0]], 0.6, label="Decreased", color=color, alpha=0.5, hatch="//")
    ax.axhline(0, color="black", linewidth=0.5)

    ax.set_xticks([0])
    ax.set_xticklabels([tissue], rotation=30, ha="right")
    ax.set_ylabel("Number of Significant Fluxes")
    if "significance_rule" in df.columns and pd.notna(df["significance_rule"].iloc[0]):
        criterion = str(df["significance_rule"].iloc[0])
    else:
        criterion = (f"|{effect_col}|>={d_thresh}" if low_n else f"p<{p_thresh}")
    title = f"Tissue Differential Flux Summary: {comp_name} — {tissue}\n({criterion})"
    ax.set_title(title, fontweight="bold", fontsize=11)
    ax.legend()
    ax.grid(axis="y", alpha=0.3)

    if low_n:
        stat_note_text = (
            f"Statistical Note: Significance defined by |{effect_col}|>={d_thresh} on "
            "biological replicates (one net-flux value per animal). "
            + _power_floor_text(df) +
            "Exploratory — validation with additional replicates required."
        )
    else:
        stat_note_text = (
            f"Statistical Note: Significance defined by Mann-Whitney U test P<{p_thresh} "
            "on biological replicates. FDR correction applied. "
            + _power_floor_text(df)
        )
    fig.text(0.5, -0.01, stat_note_text, ha="center", va="top",
             fontsize=7, style="italic", wrap=True,
             bbox=dict(boxstyle="round,pad=0.3", facecolor="#f5f5f5", alpha=0.8))

    plt.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_dir / f"tissue_diff_summary_{comp_name}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def generate_volcano_plots(output_dir: Path, condition_colors: dict,
                            tissue_colors: dict, comparisons: list,
                            fc_thresh: float = 0.5,
                            p_thresh: float = 0.05,
                            d_thresh: float = 0.5,
                            label_top_n: int = 10,
                            axis_config: dict | None = None):
    """读取 differential_flux CSV（兼容新旧目录）并生成按组织分组的火山图."""
    df_root = output_dir / "01_df_results"
    if not df_root.exists():
        print("[SKIP volcano] 01_df_results/ 目录不存在")
        return

    print("\n  Generating volcano plots...")
    for cond_a, cond_b in comparisons:
        comp_name = f"{cond_a.upper()}_vs_{cond_b.upper()}"
        for tissue, df in _iter_differential_flux_tables(df_root, comp_name):
            v_dir = output_dir / "02_result_plot" / tissue / "1_volcano"
            plot_volcano(
                df, comp_name, tissue, v_dir, condition_colors,
                fc_thresh=fc_thresh, p_thresh=p_thresh,
                d_thresh=d_thresh, label_top_n=label_top_n,
                axis_config=axis_config,
            )
            plot_tissue_diff_summary(
                df, comp_name, tissue, v_dir, tissue_colors,
                p_thresh=p_thresh, d_thresh=d_thresh,
            )
            print(f"    {tissue}/1_volcano/volcano_{comp_name}.png")
