#!/usr/bin/env python3
"""
Heatmap plots: pathway-condition heatmap, TCA cycle, convergence heatmap,
tissue flux overview.

Adapted from multiomics_integration_analysis.py lines 515-528 (convergence),
711-745 (pathway), 747-799 (TCA), 845-881 (overview).

Author: Claude Code
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from . import _config as cfg


def plot_pathway_condition_heatmap(
    mfa_pathway: pd.DataFrame | None,
    output_dir: Path | None = None,
) -> pd.DataFrame | None:
    """Heatmap of mean flux per pathway, tissue, and condition.

    Parameters
    ----------
    mfa_pathway : DataFrame or None
        MFA pathway data with columns tissue, pathway, condition, mean_flux.
    output_dir : Path, optional

    Returns
    -------
    DataFrame or None
    """
    if mfa_pathway is None:
        print("[SKIP] No MFA pathway data for heatmap")
        return None

    out = output_dir or cfg.PLOTS_DIR
    out.mkdir(parents=True, exist_ok=True)

    df = mfa_pathway.copy()

    # 校验必需的分组列
    required_keys = ["tissue", "pathway", "condition"]
    missing_keys = [c for c in required_keys if c not in df.columns]
    if missing_keys:
        print(
            f"[SKIP] Pathway heatmap: missing required columns {missing_keys}. "
            f"Available columns: {list(df.columns)}"
        )
        return None

    # 自动识别通量大小列（兼容不同命名约定）
    flux_col_candidates = [
        "mean_flux", "flux_mean", "pathway_mean_flux",
        "flux", "mean", "value", "abs_flux", "mean_abs_flux",
    ]
    flux_col = next((c for c in flux_col_candidates if c in df.columns), None)

    # 退而求其次：使用第一个非分组的数值列
    if flux_col is None:
        numeric_cols = [
            c for c in df.columns
            if c not in required_keys and pd.api.types.is_numeric_dtype(df[c])
        ]
        if numeric_cols:
            flux_col = numeric_cols[0]
            print(f"  [info] Using numeric column '{flux_col}' as flux value")

    if flux_col is None:
        print(
            f"[SKIP] Pathway heatmap: no flux value column found. "
            f"Tried {flux_col_candidates}. Available columns: {list(df.columns)}"
        )
        return None

    agg = (
        df.groupby(required_keys)[flux_col]
        .mean()
        .reset_index()
        .rename(columns={flux_col: "mean_flux"})
    )
    agg["label"] = agg["tissue"].astype(str) + " | " + agg["pathway"].astype(str)

    pivot = agg.pivot_table(index="label", columns="condition", values="mean_flux")
    if pivot.empty:
        print("  No data for pathway heatmap")
        return None

    pivot_log = np.log1p(pivot.abs()) * np.sign(pivot)

    fig, ax = plt.subplots(
        figsize=(max(8, pivot.shape[1] * 2), max(8, pivot.shape[0] * 0.35))
    )
    sns.heatmap(
        pivot_log, cmap="RdBu_r", center=0, ax=ax, linewidths=0.3,
        xticklabels=True, yticklabels=True,
        cbar_kws={"label": "log(1+|mean flux|) * sign"},
    )
    ax.set_title("MFA: Mean Pathway Flux by Tissue and Condition", fontweight="bold")
    ax.set_ylabel("")
    ax.set_xlabel("Condition")
    plt.tight_layout()
    fig.savefig(out / "mfa_pathway_condition_heatmap.png", dpi=300)
    plt.close(fig)

    common_dir = cfg.COMMON_DIR
    common_dir.mkdir(parents=True, exist_ok=True)
    pivot.to_csv(common_dir / "mfa_pathway_condition_summary.csv", encoding="utf-8-sig")

    print("  [heatmap] Saved mfa_pathway_condition_heatmap.png + summary CSV")
    return pivot

def plot_tca_cycle(
    mfa_tca: pd.DataFrame | None,
    output_dir: Path | None = None,
) -> pd.DataFrame | None:
    """TCA cycle cross-tissue comparison heatmap and top changes bar chart.

    Parameters
    ----------
    mfa_tca : DataFrame or None
        TCA data with columns tissue, reaction, condition, mean.
    output_dir : Path, optional

    Returns
    -------
    DataFrame or None
    """
    if mfa_tca is None:
        print("[SKIP] No MFA TCA cycle data")
        return None

    out = output_dir or cfg.PLOTS_DIR
    out.mkdir(parents=True, exist_ok=True)

    df = mfa_tca.copy()
    agg = df.groupby(["tissue", "reaction", "condition"])["mean"].mean().reset_index()

    pivot = agg.pivot_table(index=["tissue", "reaction"], columns="condition", values="mean")
    if pivot.empty:
        print("  No TCA data for heatmap")
        return None

    pivot.to_csv(cfg.COMMON_DIR / "mfa_tca_cross_tissue.csv", encoding="utf-8-sig")

    pivot_log = np.log1p(pivot.abs()) * np.sign(pivot)
    fig, ax = plt.subplots(
        figsize=(max(8, pivot.shape[1] * 2.5), max(6, pivot.shape[0] * 0.4))
    )
    sns.heatmap(
        pivot_log, cmap="YlOrRd", ax=ax, linewidths=0.3,
        xticklabels=True, yticklabels=True,
        cbar_kws={"label": "log(1+|flux|)"},
    )
    ax.set_title("TCA Cycle Flux Comparison Across Tissues", fontweight="bold")
    ax.set_ylabel("")
    plt.tight_layout()
    fig.savefig(out / "mfa_tca_cross_tissue_heatmap.png", dpi=300)
    plt.close(fig)

    # Condition difference bar chart
    conditions = sorted(df["condition"].unique())
    if len(conditions) >= 2:
        c1, c2 = conditions[0], conditions[-1]
        p1 = agg[agg["condition"] == c1].set_index(["tissue", "reaction"])["mean"]
        p2 = agg[agg["condition"] == c2].set_index(["tissue", "reaction"])["mean"]
        diff = (p1 - p2).dropna().sort_values()
        if not diff.empty:
            top = pd.concat([diff.head(10), diff.tail(10)]).drop_duplicates()
            fig2, ax2 = plt.subplots(figsize=(10, max(5, len(top) * 0.4)))
            colors = ["#8172B2" if v > 0 else "#C44E52" for v in top.values]
            top.plot.barh(ax=ax2, color=colors, alpha=0.8)
            ax2.set_xlabel(f"Flux Difference ({c1} - {c2})")
            ax2.set_title(f"TCA Cycle: Top Flux Changes ({c1} vs {c2})", fontweight="bold")
            ax2.axvline(0, color="k", linewidth=0.5)
            ax2.grid(axis="x", alpha=0.3)
            plt.tight_layout()
            fig2.savefig(out / "mfa_tca_top_changes.png", dpi=300)
            plt.close(fig2)

    print("  [tca] Saved TCA cycle heatmap + top changes plots")
    return pivot


def plot_convergence_heatmap(
    conv_df: pd.DataFrame,
    output_dir: Path | None = None,
) -> None:
    """Heatmap: tissues vs omics-pair, cell = number of shared pathways.

    Parameters
    ----------
    conv_df : DataFrame
        From PathwayConvergenceAnalyzer with columns tissue, comparison,
        omics_A, omics_B, pathway_term.
    output_dir : Path, optional
    """
    if conv_df.empty:
        return

    out = output_dir or cfg.PLOTS_DIR
    out.mkdir(parents=True, exist_ok=True)

    conv_df = conv_df.copy()
    conv_df["pair"] = conv_df["omics_A"] + " & " + conv_df["omics_B"]
    if "comparison" in conv_df.columns:
        conv_df["context"] = (
            conv_df["tissue"].astype(str)
            + " | "
            + conv_df["comparison"].astype(str)
        )
    else:
        conv_df["context"] = conv_df["tissue"].astype(str)
    pivot = conv_df.groupby(["context", "pair"]).size().unstack(fill_value=0)

    fig, ax = plt.subplots(
        figsize=(max(6, pivot.shape[1] * 2), max(4, pivot.shape[0] * 0.8))
    )
    sns.heatmap(pivot, annot=True, fmt="d", cmap="YlOrRd", ax=ax, linewidths=0.5)
    ax.set_title("Cross-Omics Pathway Convergence (# shared terms)", fontweight="bold")
    ax.set_ylabel("Tissue | Comparison")
    ax.set_xlabel("Omics Pair")
    plt.tight_layout()
    fig.savefig(out / "pathway_convergence_heatmap.png", dpi=300)
    plt.close(fig)


def plot_tissue_flux_overview(
    mfa_overview: pd.DataFrame | None,
    output_dir: Path | None = None,
) -> pd.DataFrame | None:
    """Bar chart of total/mean metabolic flux per tissue and condition.

    Parameters
    ----------
    mfa_overview : DataFrame or None
        With columns tissue, condition, total_abs_flux, mean_abs_flux.
    output_dir : Path, optional

    Returns
    -------
    DataFrame or None
    """
    if mfa_overview is None:
        print("[SKIP] No MFA tissue overview data")
        return None

    out = output_dir or cfg.PLOTS_DIR
    out.mkdir(parents=True, exist_ok=True)

    df = mfa_overview.copy()
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Total absolute flux
    if "total_abs_flux" in df.columns:
        pivot_total = df.pivot_table(
            index="tissue", columns="condition",
            values="total_abs_flux", aggfunc="mean",
        )
        pivot_total.plot(kind="bar", ax=axes[0], alpha=0.8, rot=30)
        axes[0].set_title("Total Absolute Flux per Tissue", fontweight="bold")
        axes[0].set_ylabel("Total |Flux|")
        axes[0].grid(axis="y", alpha=0.3)
        axes[0].legend(title="Condition")

    # Mean absolute flux
    if "mean_abs_flux" in df.columns:
        pivot_mean = df.pivot_table(
            index="tissue", columns="condition",
            values="mean_abs_flux", aggfunc="mean",
        )
        pivot_mean.plot(kind="bar", ax=axes[1], alpha=0.8, rot=30)
        axes[1].set_title("Mean Absolute Flux per Reaction", fontweight="bold")
        axes[1].set_ylabel("Mean |Flux|")
        axes[1].grid(axis="y", alpha=0.3)
        axes[1].legend(title="Condition")

    plt.suptitle("MFA: Tissue Metabolic Activity Overview", fontweight="bold", fontsize=13)
    plt.tight_layout()
    fig.savefig(out / "mfa_tissue_flux_overview.png", dpi=300)
    plt.close(fig)

    common_dir = cfg.COMMON_DIR
    common_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(common_dir / "mfa_tissue_overview_integrated.csv", index=False, encoding="utf-8-sig")

    print("  [overview] Saved tissue flux overview plots + CSV")
    return df
