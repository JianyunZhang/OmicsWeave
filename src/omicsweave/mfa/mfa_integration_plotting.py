#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MFA Integration Plotting

从整合阶段输出的 CSV（integration/tracer_harmonized 与
integration/cross_tracer_scoring）生成整合可视化图表。
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _safe_name(name: str) -> str:
    """将任意字符串转换为文件名安全片段。"""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(name)).strip("_") or "unnamed"


def _clip_text(text: str, max_len: int) -> str:
    s = str(text)
    if len(s) <= max_len:
        return s
    return s[: max(0, max_len - 3)] + "..."


def _load_meta_tables(harmonized_root: Path) -> dict[str, list[tuple[str, pd.DataFrame]]]:
    """加载 tracer_harmonized 下所有 meta 分析 CSV。"""
    result: dict[str, list[tuple[str, pd.DataFrame]]] = {}
    if not harmonized_root.exists():
        return result

    for csv_file in sorted(harmonized_root.glob("*/differential_flux_meta_*.csv")):
        protocol_id = csv_file.parent.name
        comp_name = csv_file.stem.replace("differential_flux_meta_", "")
        df = pd.read_csv(csv_file)
        if df.empty:
            continue
        result.setdefault(protocol_id, []).append((comp_name, df))
    return result


def _plot_meta_forest(
    df: pd.DataFrame,
    protocol_id: str,
    comp_name: str,
    out_file: Path,
    top_n: int = 20,
    max_label_len: int = 46,
) -> bool:
    """绘制 pooled effect 森林图。"""
    required = {"diff_pooled", "ci_lower", "ci_upper"}
    if not required.issubset(df.columns):
        print(f"  [SKIP plot] 缺少森林图所需列: {required - set(df.columns)}")
        return False

    plot_df = df.dropna(subset=["diff_pooled", "ci_lower", "ci_upper"]).copy()
    if plot_df.empty:
        return False

    plot_df = plot_df[np.isfinite(plot_df["diff_pooled"])]
    plot_df = plot_df[np.isfinite(plot_df["ci_lower"]) & np.isfinite(plot_df["ci_upper"])]
    if plot_df.empty:
        return False

    plot_df["abs_diff_pooled"] = plot_df["diff_pooled"].abs()
    plot_df = plot_df.sort_values("abs_diff_pooled", ascending=False).head(top_n).copy()
    if plot_df.empty:
        return False

    reaction_col = "reaction" if "reaction" in plot_df.columns else "flux_id"
    tissue_col = "tissue" if "tissue" in plot_df.columns else None

    labels = []
    for _, row in plot_df.iterrows():
        reaction = _clip_text(row.get(reaction_col, ""), max_label_len)
        if tissue_col:
            labels.append(f"{reaction} | {row.get(tissue_col, '')}")
        else:
            labels.append(reaction)

    y = np.arange(len(plot_df))
    x = plot_df["diff_pooled"].values
    lo = plot_df["ci_lower"].values
    hi = plot_df["ci_upper"].values

    colors = []
    sig_col = "is_significant_pooled" if "is_significant_pooled" in plot_df.columns else None
    for i in range(len(plot_df)):
        xi = x[i]
        is_sig = bool(plot_df.iloc[i][sig_col]) if sig_col else False
        if is_sig and xi >= 0:
            colors.append("#D55E00")
        elif is_sig and xi < 0:
            colors.append("#0072B2")
        else:
            colors.append("#777777")

    fig_h = max(6, 1.8 + len(plot_df) * 0.34)
    fig, ax = plt.subplots(figsize=(12, fig_h))
    for i in range(len(plot_df)):
        ax.plot([lo[i], hi[i]], [y[i], y[i]], color="#999999", linewidth=1.6, zorder=1)
        ax.scatter(x[i], y[i], color=colors[i], s=36, zorder=2)
    ax.axvline(0.0, color="#333333", linestyle="--", linewidth=1.0, alpha=0.8)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Pooled effect diff (95% CI)")
    ax.set_title(f"Meta整合森林图 | {protocol_id} | {comp_name}", fontweight="bold")
    ax.grid(axis="x", alpha=0.25)

    out_file.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return True


def _plot_heterogeneity_distribution(
    protocol_id: str,
    meta_entries: list[tuple[str, pd.DataFrame]],
    out_file: Path,
    heterogeneity_threshold: float,
) -> bool:
    """绘制 I² 分布图。"""
    frames = []
    for comp_name, df in meta_entries:
        if "I_squared" not in df.columns:
            continue
        tmp = df.copy()
        tmp["comparison"] = comp_name
        frames.append(tmp)

    if not frames:
        return False

    merged = pd.concat(frames, ignore_index=True)
    merged = merged.dropna(subset=["I_squared"])
    merged = merged[np.isfinite(merged["I_squared"])]
    if merged.empty:
        return False

    threshold_pct = float(heterogeneity_threshold) * 100.0

    fig, axes = plt.subplots(
        1, 2, figsize=(14, 5), gridspec_kw={"width_ratios": [1.2, 1.4]}
    )

    axes[0].hist(merged["I_squared"], bins=20, color="#4C72B0", alpha=0.75, edgecolor="white")
    axes[0].axvline(threshold_pct, color="#C44E52", linestyle="--", linewidth=1.4)
    axes[0].set_xlabel("I² (%)")
    axes[0].set_ylabel("Count")
    axes[0].set_title(f"I²分布（{protocol_id}）")
    axes[0].grid(axis="y", alpha=0.25)

    comp_order = sorted(merged["comparison"].unique().tolist())
    if len(comp_order) >= 2:
        data = [merged.loc[merged["comparison"] == c, "I_squared"].values for c in comp_order]
        axes[1].boxplot(data, labels=comp_order, showfliers=False)
        axes[1].tick_params(axis="x", rotation=30, labelsize=8)
        axes[1].set_ylabel("I² (%)")
        axes[1].set_title("各比较组异质性分布")
        axes[1].grid(axis="y", alpha=0.25)
    else:
        summary = (
            merged.groupby("comparison", as_index=False)["I_squared"]
            .mean()
            .rename(columns={"I_squared": "mean_i2"})
        )
        axes[1].bar(summary["comparison"], summary["mean_i2"], color="#55A868", alpha=0.8)
        axes[1].tick_params(axis="x", rotation=30, labelsize=8)
        axes[1].set_ylabel("Mean I² (%)")
        axes[1].set_title("平均异质性")
        axes[1].grid(axis="y", alpha=0.25)

    out_file.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return True


def _plot_consensus_score_bar(
    df: pd.DataFrame,
    comp_name: str,
    out_file: Path,
    top_n: int = 25,
    max_label_len: int = 42,
) -> bool:
    """绘制方向一致性得分 Top-N 条形图。"""
    required = {"directional_consensus_score", "consensus_pass"}
    if not required.issubset(df.columns):
        print(f"  [SKIP plot] 缺少一致性条形图所需列: {required - set(df.columns)}")
        return False

    plot_df = df.dropna(subset=["directional_consensus_score"]).copy()
    plot_df = plot_df[np.isfinite(plot_df["directional_consensus_score"])]
    if plot_df.empty:
        return False

    plot_df = plot_df.sort_values("directional_consensus_score", ascending=False).head(top_n).copy()
    if plot_df.empty:
        return False

    reaction_col = "reaction" if "reaction" in plot_df.columns else "flux_id"
    tissue_col = "tissue" if "tissue" in plot_df.columns else None

    labels = []
    for _, row in plot_df.iterrows():
        reaction = _clip_text(row.get(reaction_col, ""), max_label_len)
        if tissue_col:
            labels.append(f"{reaction} | {row.get(tissue_col, '')}")
        else:
            labels.append(reaction)

    y = np.arange(len(plot_df))
    x = plot_df["directional_consensus_score"].values
    colors = ["#2E8B57" if bool(v) else "#C44E52" for v in plot_df["consensus_pass"].values]

    fig_h = max(6, 2 + len(plot_df) * 0.3)
    fig, ax = plt.subplots(figsize=(12, fig_h))
    ax.barh(y, x, color=colors, alpha=0.88)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.05)
    ax.axvline(1.0, color="#333333", linestyle="--", linewidth=1.0, alpha=0.7)
    ax.set_xlabel("Directional consensus score")
    ax.set_title(f"跨Protocol方向一致性 Top-{len(plot_df)} | {comp_name}", fontweight="bold")
    ax.grid(axis="x", alpha=0.25)

    out_file.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return True


def _plot_consensus_summary(summary_file: Path, out_file: Path) -> bool:
    """绘制比较组一致性通过率汇总图。"""
    if not summary_file.exists():
        return False

    df = pd.read_csv(summary_file)
    required = {"comparison", "n_fluxes", "n_consensus_pass"}
    if not required.issubset(df.columns):
        return False
    if df.empty:
        return False

    df = df.copy()
    df["pass_rate"] = np.where(df["n_fluxes"] > 0, df["n_consensus_pass"] / df["n_fluxes"], np.nan)
    df = df.sort_values("pass_rate", ascending=False)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(df["comparison"], df["pass_rate"], color="#4C72B0", alpha=0.85)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Consensus pass rate")
    ax.set_title("跨Protocol一致性通过率汇总", fontweight="bold")
    ax.grid(axis="y", alpha=0.25)
    ax.tick_params(axis="x", rotation=30)

    for i, row in df.reset_index(drop=True).iterrows():
        ax.text(
            i,
            min(1.03, float(row["pass_rate"]) + 0.02),
            f"{int(row['n_consensus_pass'])}/{int(row['n_fluxes'])}",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    out_file.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return True


def generate_integration_plots(output_root: Path, integration_config: dict[str, Any]) -> list[Path]:
    """按 integration.plotting 配置生成整合图表。"""
    plotting_cfg = integration_config.get("plotting", {}) or {}
    if not plotting_cfg.get("enabled", False):
        print("  [integration] integration.plotting.enabled=false，跳过整合绘图")
        return []

    out_subdir = str(plotting_cfg.get("output_subdir", "02_result_plot")).strip() or "02_result_plot"
    plots_cfg = plotting_cfg.get("plots", {}) or {}
    top_n = int(plotting_cfg.get("top_n", 20))
    max_label_len = int(plotting_cfg.get("max_label_len", 46))

    harmonized_root = output_root / "integration" / "tracer_harmonized"
    consensus_root = output_root / "integration" / "cross_tracer_scoring"
    plot_root = output_root / "integration" / out_subdir
    plot_root.mkdir(parents=True, exist_ok=True)

    generated: list[Path] = []

    # 1) 同 protocol Meta 森林图 + I² 分布图
    protocol_meta = _load_meta_tables(harmonized_root)
    if plots_cfg.get("meta_forest", True):
        for protocol_id, entries in sorted(protocol_meta.items()):
            for comp_name, df in entries:
                out_file = (
                    plot_root
                    / "tracer_harmonized"
                    / protocol_id
                    / f"meta_forest_{_safe_name(comp_name)}.png"
                )
                if _plot_meta_forest(
                    df=df,
                    protocol_id=protocol_id,
                    comp_name=comp_name,
                    out_file=out_file,
                    top_n=top_n,
                    max_label_len=max_label_len,
                ):
                    generated.append(out_file)

    if plots_cfg.get("heterogeneity_hist", True):
        het_threshold = (
            integration_config.get("within_protocol", {}) or {}
        ).get("heterogeneity_threshold", 0.50)
        for protocol_id, entries in sorted(protocol_meta.items()):
            out_file = (
                plot_root
                / "tracer_harmonized"
                / protocol_id
                / f"heterogeneity_i2_{_safe_name(protocol_id)}.png"
            )
            if _plot_heterogeneity_distribution(
                protocol_id=protocol_id,
                meta_entries=entries,
                out_file=out_file,
                heterogeneity_threshold=float(het_threshold),
            ):
                generated.append(out_file)

    # 2) 跨 protocol 方向一致性图
    if plots_cfg.get("consensus_bar", True) and consensus_root.exists():
        for csv_file in sorted(consensus_root.glob("directional_consensus_*.csv")):
            comp_name = csv_file.stem.replace("directional_consensus_", "")
            df = pd.read_csv(csv_file)
            if df.empty:
                continue
            out_file = (
                plot_root
                / "cross_tracer_scoring"
                / f"directional_consensus_{_safe_name(comp_name)}.png"
            )
            if _plot_consensus_score_bar(
                df=df,
                comp_name=comp_name,
                out_file=out_file,
                top_n=top_n,
                max_label_len=max_label_len,
            ):
                generated.append(out_file)

        summary_file = consensus_root / "protocol_comparison_summary.csv"
        summary_out = plot_root / "cross_tracer_scoring" / "protocol_comparison_summary.png"
        if _plot_consensus_summary(summary_file, summary_out):
            generated.append(summary_out)

    if generated:
        print(f"  [integration] 整合绘图完成，共生成 {len(generated)} 个图:")
        for p in generated:
            rel = p.relative_to(output_root) if p.is_relative_to(output_root) else p
            print(f"    - {rel}")
    else:
        print("  [integration] 未生成整合图（可能无可用整合结果表）")

    return generated
