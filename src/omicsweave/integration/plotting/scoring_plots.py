#!/usr/bin/env python3
"""
证据画像的图：证据等级分布、支持度与冲突度散点、各组织的支持度分布。

旧版的 0–4 总分直方图与置信度饼图已随加和总分一并删除（阶段 2.8）。
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

from . import _config as cfg

# 锚定等级在前，锚点不成立时的组学一致性等级（Omics-*）在后
GRADE_ORDER = ["Replicated", "Convergent", "Single-layer", "Exploratory", "Conflicting", "Omics-concordant",
               "Omics-single-layer", "Omics-discordant", "No-signal", "Not-assessable"]
_DEFAULT_GRADE_COLORS = {
    "Replicated": "#1B7837", "Convergent": "#5AAE61", "Single-layer": "#A6DBA0", "Exploratory": "#F4A582",
    "Conflicting": "#B2182B", "Omics-concordant": "#4393C3", "Omics-single-layer": "#92C5DE",
    "Omics-discordant": "#D6604D", "No-signal": "#D9D9D9", "Not-assessable": "#F0F0F0",
}


def plot_evidence_profile(
    scored_df: pd.DataFrame,
    output_dir: Path | None = None,
) -> list[Path]:
    """证据画像的汇总图。

    Parameters
    ----------
    scored_df : DataFrame
        :meth:`FluxEvidenceBuilder.build_all` 的输出，含 evidence_grade、support_score、conflict_score、tissue。
    output_dir : Path, optional

    Returns
    -------
    list[Path]
        写出的图件。
    """
    if scored_df.empty:
        return []

    out = output_dir or cfg.PLOTS_DIR
    out.mkdir(parents=True, exist_ok=True)
    colors = {**_DEFAULT_GRADE_COLORS, **cfg.PLOT_COLORS.get("evidence_grade", {})}
    saved: list[Path] = []

    # 1. 证据等级分布
    counts = scored_df["evidence_grade"].value_counts()
    grades = [g for g in GRADE_ORDER if g in counts.index] + [g for g in counts.index if g not in GRADE_ORDER]
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(grades, [counts[g] for g in grades], color=[colors.get(g, "#999999") for g in grades],
           edgecolor="#555555", linewidth=0.5)
    ax.set_ylabel("Reactions")
    ax.set_title("Evidence Grade Distribution", fontweight="bold")
    ax.grid(axis="y", alpha=0.3)
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    path = out / "evidence_grades.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    saved.append(path)

    # 2. 支持度与冲突度
    fig2, ax2 = plt.subplots(figsize=(7, 6))
    for grade in grades:
        rows = scored_df[scored_df["evidence_grade"] == grade]
        ax2.scatter(rows["support_score"], rows["conflict_score"], s=18, alpha=0.75,
                    color=colors.get(grade, "#999999"), label=grade, edgecolor="none")
    ax2.set_xlabel("Support score")
    ax2.set_ylabel("Conflict score")
    ax2.set_xlim(-0.02, 1.02)
    ax2.set_ylim(-0.02, 1.02)
    ax2.set_title("Support vs Conflict per Reaction", fontweight="bold")
    ax2.legend(fontsize=8, frameon=False)
    ax2.grid(alpha=0.3)
    plt.tight_layout()
    path = out / "support_conflict.png"
    fig2.savefig(path, dpi=300)
    plt.close(fig2)
    saved.append(path)

    # 3. 各组织的支持度
    tissues = sorted(scored_df["tissue"].unique()) if "tissue" in scored_df.columns else []
    if len(tissues) > 1:
        fig3, ax3 = plt.subplots(figsize=(max(8, len(tissues) * 1.5), 6))
        bp = ax3.boxplot([scored_df.loc[scored_df["tissue"] == t, "support_score"].values for t in tissues],
                         patch_artist=True)
        ax3.set_xticks(range(1, len(tissues) + 1), [cfg.TISSUE_LABELS.get(t, t) for t in tissues])
        tissue_colors = cfg.PLOT_COLORS.get("tissue", {})
        for patch, tissue in zip(bp["boxes"], tissues):
            patch.set_facecolor(tissue_colors.get(tissue, "#4C72B0"))
            patch.set_alpha(0.6)
        ax3.set_ylabel("Support score")
        ax3.set_title("Support Score by Tissue", fontweight="bold")
        ax3.grid(axis="y", alpha=0.3)
        plt.xticks(rotation=30, ha="right")
        plt.tight_layout()
        path = out / "support_by_tissue.png"
        fig3.savefig(path, dpi=300)
        plt.close(fig3)
        saved.append(path)

    print(f"  [evidence] Saved {len(saved)} evidence profile plots")
    return saved
