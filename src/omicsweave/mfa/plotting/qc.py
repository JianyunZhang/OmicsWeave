#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
质控图（拟合 loss 柱状图 + 直方图）

优先从 01_df_results/all_tissues/qc/ 读取 QC CSV，兼容旧路径
01_df_results/qc/，输出到 02_result_plot/all_tissues/0_qc/。
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path


def plot_qc(qc_df: pd.DataFrame, output_dir: Path, condition_colors: dict):
    """Bar chart of best-fit loss per sample + histogram by condition."""
    qc_dir = output_dir / "02_result_plot" / "all_tissues" / "0_qc"
    qc_dir.mkdir(parents=True, exist_ok=True)

    n_samples = len(qc_df)
    # 根据样本数自适应画布宽度，保证标签不重叠
    fig_w = max(12, 3 + n_samples * 0.55)
    fig, axes = plt.subplots(1, 2, figsize=(fig_w, 6),
                             gridspec_kw={"width_ratios": [max(2, n_samples / 4), 1]})

    # ── 左图：每样本最优拟合 loss 柱状图 ──
    for cond in sorted(qc_df["condition"].unique()):
        subset = qc_df[qc_df["condition"] == cond]
        color = condition_colors.get(cond, "#999")
        axes[0].bar(
            subset["sample"].str.replace("__", "\n", regex=False),
            subset["min_loss"], alpha=0.8, label=cond, color=color,
        )
    axes[0].set_ylabel("Best-Fit Loss (SSR)")
    axes[0].set_title("Best-Fit Loss per Sample", fontweight="bold")
    axes[0].legend()
    # 自适应字体大小：样本越多字体越小
    tick_fs = max(5, min(9, 140 // n_samples))
    axes[0].tick_params(axis="x", rotation=55, labelsize=tick_fs)
    axes[0].grid(axis="y", alpha=0.3)

    # ── 右图：按条件分组的 loss 分布直方图 + 图注 ──
    for cond in sorted(qc_df["condition"].unique()):
        vals = qc_df[qc_df["condition"] == cond]["min_loss"].values
        color = condition_colors.get(cond, "#999")
        axes[1].hist(vals, bins=max(3, len(vals) // 2), alpha=0.6,
                     label=cond, color=color, edgecolor="white")
    axes[1].set_xlabel("Best-Fit Loss")
    axes[1].set_ylabel("Count")
    axes[1].set_title("Loss Distribution by Condition", fontweight="bold")
    axes[1].legend()
    axes[1].grid(axis="y", alpha=0.3)

    # 右图下方添加图注，说明含义
    n_fits = (int(qc_df["n_fits"].median())
              if "n_fits" in qc_df.columns and qc_df["n_fits"].notna().any()
              else None)
    iterations = f"{n_fits} fitting iterations" if n_fits else "all fitting iterations"
    caption = (
        f"Each bar represents the best-fit loss (min SSR across {iterations}) "
        "per sample.\n"
        "Wider distribution indicates higher heterogeneity among biological "
        "replicates within a condition."
    )
    axes[1].annotate(
        caption, xy=(0.5, 0), xytext=(0.5, -0.22),
        xycoords="axes fraction", textcoords="axes fraction",
        fontsize=7, ha="center", va="top", style="italic",
        wrap=True,
    )

    plt.tight_layout()
    plt.savefig(qc_dir / "qc_fitting_loss.png", dpi=300, bbox_inches="tight")
    plt.close(fig)




def plot_infusion_glycemia(output_dir: Path, condition_colors: dict) -> bool:
    """灌注期血糖时间曲线 + 稳态窗口。

    MFA 假设采样时代谢已达稳态。这张图让"假设成不成立"一眼可见：稳态窗口
    内的曲线应当是平的，持续上升或大幅波动都意味着拟合出的通量描述的是一段
    漂移过程，而非某个确定的代谢状态。
    """
    qc_dir = output_dir / "01_df_results" / "all_tissues" / "qc"
    timecourse_path = qc_dir / "infusion_glycemia_timecourse.csv"
    summary_path = qc_dir / "infusion_steady_state_qc.csv"
    if not timecourse_path.exists():
        return False

    timecourse = pd.read_csv(timecourse_path)
    if timecourse.empty:
        return False
    summary = pd.read_csv(summary_path) if summary_path.exists() else pd.DataFrame()
    passed = (dict(zip(summary["mouse_id"].astype(str), summary["passed"]))
              if not summary.empty else {})

    window = timecourse[timecourse["in_steady_state_window"]]
    window_start = float(window["minute"].min()) if len(window) else None
    window_end = float(window["minute"].max()) if len(window) else None

    conditions = sorted(timecourse["condition"].astype(str).unique())
    fig, ax = plt.subplots(figsize=(11, 6))

    for condition in conditions:
        color = condition_colors.get(condition, "#999")
        subset = timecourse[timecourse["condition"].astype(str) == condition]
        for mouse_id, mouse_rows in subset.groupby("mouse_id"):
            mouse_rows = mouse_rows.sort_values("minute")
            ok = bool(passed.get(str(mouse_id), False))
            ax.plot(
                mouse_rows["minute"], mouse_rows["glucose_mmol_l"],
                marker="o", markersize=3, linewidth=1.4,
                color=color, alpha=0.9 if ok else 0.55,
                linestyle="-" if ok else "--",
            )

    if window_start is not None:
        ax.axvspan(window_start, window_end, color="#4CAF50", alpha=0.08)
        ax.axvline(window_start, color="#4CAF50", linestyle=":", linewidth=1)
        ax.axvline(window_end, color="#4CAF50", linestyle=":", linewidth=1)
        ax.text(
            (window_start + window_end) / 2, ax.get_ylim()[1],
            f"steady-state window [{window_start:.0f}, {window_end:.0f}] min",
            ha="center", va="bottom", fontsize=8, color="#2E7D32",
        )

    handles = [plt.Line2D([], [], color=condition_colors.get(c, "#999"), label=c)
               for c in conditions]
    handles.append(plt.Line2D([], [], color="#555", linestyle="-", label="steady state: pass"))
    handles.append(plt.Line2D([], [], color="#555", linestyle="--", label="steady state: fail"))
    ax.legend(handles=handles, fontsize=8, ncol=2)

    ax.set_xlabel("Time after infusion start (min)")
    ax.set_ylabel("Blood glucose (mmol/L)")
    n_pass = sum(1 for v in passed.values() if v)
    ax.set_title(
        f"Infusion Glycemia Time Course — steady state {n_pass}/{len(passed)} mice"
        if passed else "Infusion Glycemia Time Course",
        fontweight="bold",
    )
    ax.grid(alpha=0.25)

    fig.text(
        0.5, -0.02,
        "Statistical Note: 13C-MFA assumes metabolic and isotopic steady state at "
        "sampling. Fluxes fitted from a drifting glycemia describe the average of a "
        "transient, not a defined physiological state; between-group comparison of "
        "such fluxes is not interpretable.",
        ha="center", va="top", fontsize=7, style="italic", wrap=True,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#f5f5f5", alpha=0.8),
    )

    qc_plot_dir = output_dir / "02_result_plot" / "all_tissues" / "0_qc"
    qc_plot_dir.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(qc_plot_dir / "qc_infusion_glycemia.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return True


def plot_mid_residuals(output_dir: Path) -> bool:
    """实测 vs 模拟 MID 的拟合质量。

    左：全部数据点的实测-模拟散点（理想情况落在 y=x 上）。
    右：失拟最严重的代谢物，定位模型在哪里复现不了实测标记模式。
    """
    mid_dir = output_dir / "01_df_results" / "all_tissues" / "mid"
    long_path = mid_dir / "mid_residuals_long.csv"
    by_metabolite_path = mid_dir / "mid_fit_quality_by_metabolite.csv"
    if not long_path.exists() or not by_metabolite_path.exists():
        return False

    long_df = pd.read_csv(long_path).dropna(subset=["measured", "simulated"])
    by_metabolite = pd.read_csv(by_metabolite_path)
    if long_df.empty or by_metabolite.empty:
        return False

    fig, axes = plt.subplots(1, 2, figsize=(14, 6),
                             gridspec_kw={"width_ratios": [1, 1.3]})

    residual = (long_df["simulated"] - long_df["measured"]).abs()
    axes[0].scatter(long_df["measured"], long_df["simulated"],
                    c=residual, cmap="YlOrRd", s=6, alpha=0.5)
    axes[0].plot([0, 1], [0, 1], color="#333", linestyle="--", linewidth=1)
    axes[0].set_xlabel("Measured MID (fraction)")
    axes[0].set_ylabel("Simulated MID (fraction)")
    axes[0].set_title("Measured vs Simulated MID", fontweight="bold")
    axes[0].grid(alpha=0.25)

    top = by_metabolite.head(15).iloc[::-1]
    labels = [str(m)[:44] for m in top["metabolite"]]
    axes[1].barh(labels, top["median_rmse"], color="#E24A33", alpha=0.85)
    axes[1].set_xlabel("Median RMSE (measured vs simulated)")
    axes[1].set_title("Worst-Fitting Metabolites", fontweight="bold")
    axes[1].tick_params(axis="y", labelsize=7)
    axes[1].axvline(0.05, color="#333", linestyle="--", linewidth=1)
    axes[1].text(0.05, -0.7, "poor-fit threshold 0.05", fontsize=7, color="#333")
    axes[1].grid(axis="x", alpha=0.25)

    fig.text(
        0.5, -0.02,
        "Statistical Note: MID values are fractions summing to 1, so residuals are on "
        "an absolute mole-fraction scale. Typical MS measurement error is 0.003-0.02; "
        "RMSE >= 0.05 means the network cannot reproduce the observed labelling. "
        "Fluxes derived from a model that misfits the data are not interpretable, "
        "regardless of their statistics.",
        ha="center", va="top", fontsize=7, style="italic", wrap=True,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#f5f5f5", alpha=0.8),
    )

    qc_plot_dir = output_dir / "02_result_plot" / "all_tissues" / "0_qc"
    qc_plot_dir.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(qc_plot_dir / "qc_mid_fit_residuals.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return True


def plot_flux_identifiability(output_dir: Path) -> bool:
    """通量可辨识性分布：拟合 IQR 宽度相对于通量自身量级。"""
    path = (output_dir / "01_df_results" / "all_tissues"
            / "flux_matrix" / "flux_identifiability.csv")
    if not path.exists():
        return False

    table = pd.read_csv(path)
    values = pd.to_numeric(table.get("median_rel_iqr"), errors="coerce").dropna()
    if values.empty:
        return False

    fig, axes = plt.subplots(1, 2, figsize=(13, 5),
                             gridspec_kw={"width_ratios": [1.4, 1]})

    clipped = values.clip(upper=5.0)
    axes[0].hist(clipped, bins=50, color="#4682B4", alpha=0.8, edgecolor="white")
    axes[0].axvline(0.2, color="#2E7D32", linestyle="--", linewidth=1.2)
    axes[0].axvline(1.0, color="#C62828", linestyle="--", linewidth=1.2)
    axes[0].text(0.2, axes[0].get_ylim()[1] * 0.95, " well-constrained < 0.2",
                 fontsize=7, color="#2E7D32", va="top")
    axes[0].text(1.0, axes[0].get_ylim()[1] * 0.85, " poorly constrained >= 1.0",
                 fontsize=7, color="#C62828", va="top")
    axes[0].set_xlabel("Median relative IQR of model fits (clipped at 5)")
    axes[0].set_ylabel("Number of fluxes")
    axes[0].set_title("Flux Identifiability Distribution", fontweight="bold")
    axes[0].grid(axis="y", alpha=0.25)

    if "identifiability" in table.columns:
        counts = table["identifiability"].value_counts()
        order = ["well_constrained", "moderate", "poorly_constrained",
                 "bound_pinned", "undetermined"]
        counts = counts.reindex([o for o in order if o in counts.index])
        colors = {"well_constrained": "#2E7D32", "moderate": "#F9A825",
                  "poorly_constrained": "#C62828", "bound_pinned": "#6A1B9A",
                  "undetermined": "#9E9E9E"}
        axes[1].bar(counts.index, counts.values,
                    color=[colors.get(i, "#999") for i in counts.index], alpha=0.85)
        for i, value in enumerate(counts.values):
            axes[1].text(i, value, str(int(value)), ha="center", va="bottom", fontsize=8)
        axes[1].set_ylabel("Number of fluxes")
        axes[1].set_title("Identifiability Classes", fontweight="bold")
        axes[1].tick_params(axis="x", rotation=20, labelsize=8)
        axes[1].grid(axis="y", alpha=0.25)

    fig.text(
        0.5, -0.03,
        "Statistical Note: the fitting iterations are NOT biological replicates, but "
        "they do describe how tightly each flux is pinned down by one animal's data. "
        "A relative IQR near or above 1 means the flux is essentially undetermined; "
        "a between-group difference in such a flux carries no information. "
        "bound_pinned: fits sit on a box bound in most animals - zero spread there "
        "means the optimizer hit the bound, not that the data determine the flux.",
        ha="center", va="top", fontsize=7, style="italic", wrap=True,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#f5f5f5", alpha=0.8),
    )

    qc_plot_dir = output_dir / "02_result_plot" / "all_tissues" / "0_qc"
    qc_plot_dir.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(qc_plot_dir / "qc_flux_identifiability.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return True

def generate_qc_plots(output_dir: Path, condition_colors: dict):
    """从已保存的 QC CSV 生成质控图."""
    print("\n  Generating QC plots...")

    csv_path = output_dir / "01_df_results" / "all_tissues" / "qc" / "qc_fitting_loss.csv"
    if not csv_path.exists():
        csv_path = output_dir / "01_df_results" / "qc" / "qc_fitting_loss.csv"
    if csv_path.exists():
        plot_qc(pd.read_csv(csv_path), output_dir, condition_colors)
        print("    all_tissues/0_qc/qc_fitting_loss.png")
    else:
        print("[SKIP qc] qc_fitting_loss.csv 不存在")

    if plot_infusion_glycemia(output_dir, condition_colors):
        print("    all_tissues/0_qc/qc_infusion_glycemia.png")
    if plot_mid_residuals(output_dir):
        print("    all_tissues/0_qc/qc_mid_fit_residuals.png")
    if plot_flux_identifiability(output_dir):
        print("    all_tissues/0_qc/qc_flux_identifiability.png")
