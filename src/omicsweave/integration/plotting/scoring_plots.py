#!/usr/bin/env python3
"""
证据画像的图：证据等级分布、支持度与冲突度散点、各组织的支持度分布。

旧版的 0–4 总分直方图与置信度饼图已随加和总分一并删除（阶段 2.8）。
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

from omicsweave.visualization.palettes import EVIDENCE_GRADE_COLORS, get_palette
from omicsweave.visualization.theme import DEFAULT_PRESET, add_economist_tag, themed

from . import _config as cfg

# 锚定等级在前，锚点不成立时的组学一致性等级（Omics-*）在后
GRADE_ORDER = ["Replicated", "Convergent", "Single-layer", "Exploratory", "Conflicting", "Omics-concordant",
               "Omics-single-layer", "Omics-discordant", "No-signal", "Not-assessable"]
# 统一采用色盲安全编辑色板（替换原红绿对比色）
_DEFAULT_GRADE_COLORS = dict(EVIDENCE_GRADE_COLORS)


def plot_evidence_profile(
    scored_df: pd.DataFrame,
    output_dir: Path | None = None,
    *,
    preset: str = DEFAULT_PRESET,
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

    pal = get_palette(preset)
    with themed(preset=preset, save_dpi=300):
        # 1. 证据等级分布
        counts = scored_df["evidence_grade"].value_counts()
        grades = [g for g in GRADE_ORDER if g in counts.index] + [g for g in counts.index if g not in GRADE_ORDER]
        fig, ax = plt.subplots(figsize=(8.5, 4.7))
        if preset == "economist":
            add_economist_tag(fig, x=0.015, y=0.985, width=0.055, height=0.015)
        ax.bar(grades, [counts[g] for g in grades], color=[colors.get(g, pal.neutral) for g in grades],
               edgecolor=pal.ink2, linewidth=0.5)
        ax.set_ylabel("Reactions")
        ax.set_title("Evidence Grade Distribution", loc="left")
        ax.grid(axis="y", alpha=0.35 if preset != "economist" else 0.85)
        plt.xticks(rotation=30, ha="right")
        plt.tight_layout(rect=(0, 0, 1, 0.96 if preset == "economist" else 1.0))
        path = out / "evidence_grades.png"
        fig.savefig(path, dpi=300, facecolor=pal.surface)
        plt.close(fig)
        saved.append(path)

        # 2. 支持度与冲突度
        fig2, ax2 = plt.subplots(figsize=(6.8, 5.7))
        if preset == "economist":
            add_economist_tag(fig2, x=0.015, y=0.985, width=0.065, height=0.014)
        for grade in grades:
            rows = scored_df[scored_df["evidence_grade"] == grade]
            ax2.scatter(rows["support_score"], rows["conflict_score"], s=20, alpha=0.82,
                        color=colors.get(grade, pal.neutral), label=grade, edgecolor="none")
        ax2.set_xlabel("Support score")
        ax2.set_ylabel("Conflict score")
        ax2.set_xlim(-0.02, 1.02)
        ax2.set_ylim(-0.02, 1.02)
        ax2.set_title("Support vs Conflict per Reaction", loc="left")
        ax2.legend(fontsize=8, frameon=False)
        ax2.grid(axis="y", alpha=0.35 if preset != "economist" else 0.85)
        plt.tight_layout(rect=(0, 0, 1, 0.96 if preset == "economist" else 1.0))
        path = out / "support_conflict.png"
        fig2.savefig(path, dpi=300, facecolor=pal.surface)
        plt.close(fig2)
        saved.append(path)

        # 3. 各组织的支持度
        tissues = sorted(scored_df["tissue"].unique()) if "tissue" in scored_df.columns else []
        if len(tissues) > 1:
            fig3, ax3 = plt.subplots(figsize=(max(7.5, len(tissues) * 1.4), 5.3))
            if preset == "economist":
                add_economist_tag(fig3, x=0.015, y=0.985, width=0.055, height=0.014)
            bp = ax3.boxplot([scored_df.loc[scored_df["tissue"] == t, "support_score"].values for t in tissues],
                             patch_artist=True)
            ax3.set_xticks(range(1, len(tissues) + 1), [cfg.TISSUE_LABELS.get(t, t) for t in tissues])
            tissue_colors = cfg.PLOT_COLORS.get("tissue", {})
            for patch, tissue in zip(bp["boxes"], tissues):
                patch.set_facecolor(tissue_colors.get(tissue, pal.down))
                patch.set_alpha(0.75)
            ax3.set_ylabel("Support score")
            ax3.set_title("Support Score by Tissue", loc="left")
            ax3.grid(axis="y", alpha=0.35 if preset != "economist" else 0.85)
            plt.xticks(rotation=30, ha="right")
            plt.tight_layout(rect=(0, 0, 1, 0.96 if preset == "economist" else 1.0))
            path = out / "support_by_tissue.png"
            fig3.savefig(path, dpi=300, facecolor=pal.surface)
            plt.close(fig3)
            saved.append(path)

    print(f"  [evidence] Saved {len(saved)} evidence profile plots")
    return saved
