"""
通路富集结果可视化：柱状图、气泡图。
"""
from __future__ import annotations

import os

import matplotlib.pyplot as plt
import numpy as np

from omicsweave.visualization.palettes import get_palette
from omicsweave.visualization.theme import DEFAULT_PRESET, add_economist_tag, themed

from . import _config  # noqa: F401


def plot_pathway_bar(
    res_df: "pd.DataFrame",
    comp_name: str,
    tissue: str,
    output_dir: str,
    *,
    preset: str = DEFAULT_PRESET,
):
    """通路富集柱状图 (Top 20, -log10 FDR，采用统一出版主题，300 dpi)。"""
    if res_df.empty:
        return
    import pandas as pd

    significance_col = "FDR" if "FDR" in res_df.columns else "P-value"
    plot_df = res_df.nsmallest(20, significance_col).copy()
    plot_df["-log10FDR"] = -np.log10(
        pd.to_numeric(plot_df[significance_col], errors="coerce").clip(lower=1e-300)
    )
    plot_df = plot_df.sort_values("-log10FDR", ascending=True)

    pal = get_palette(preset)
    with themed(preset=preset, save_dpi=300):
        fig, ax = plt.subplots(figsize=(8.5, max(4.6, len(plot_df) * 0.34)))
        if preset == "economist":
            add_economist_tag(fig, x=0.015, y=0.985, width=0.055, height=0.015)
        ax.barh(plot_df["Name"], plot_df["-log10FDR"], color=pal.up if preset == "economist" else pal.down, height=0.68, edgecolor="none")
        display_comp = comp_name.replace("_vs_", " vs ")
        ax.set_title(f"Pathway Enrichment: {display_comp} ({tissue})", loc="left")
        ax.set_xlabel(f"−log$_{{10}}$({significance_col})")
        ax.set_ylabel("")
        ax.grid(axis="x", alpha=0.45)
        ax.grid(axis="y", visible=False)
        plt.tight_layout(rect=(0, 0, 1, 0.96 if preset == "economist" else 1.0))

        out_path = os.path.join(output_dir, tissue)
        os.makedirs(out_path, exist_ok=True)
        plt.savefig(os.path.join(out_path, f"Pathway_Bar_{comp_name}.png"), dpi=300, bbox_inches="tight",
                    facecolor=pal.surface)
        plt.close(fig)


def plot_pathway_dot(
    res_df: "pd.DataFrame",
    comp_name: str,
    tissue: str,
    output_dir: str,
    *,
    preset: str = DEFAULT_PRESET,
):
    """通路富集气泡图 (Rich Factor vs Pathway，采用统一出版主题，300 dpi)。"""
    if res_df.empty:
        return

    import pandas as pd
    significance_col = "FDR" if "FDR" in res_df.columns else "P-value"
    plot_df = res_df.nsmallest(20, significance_col).copy()
    if "RichFactor" not in plot_df.columns:
        plot_df["RichFactor"] = plot_df["SigCount"] / plot_df["BgCount"]
    plot_df["-log10FDR"] = -np.log10(
        pd.to_numeric(plot_df[significance_col], errors="coerce").clip(lower=1e-300)
    )

    pal = get_palette(preset)
    with themed(preset=preset, save_dpi=300):
        fig, ax = plt.subplots(figsize=(8.5, max(4.6, len(plot_df) * 0.36)))
        if preset == "economist":
            add_economist_tag(fig, x=0.015, y=0.985, width=0.055, height=0.015)
        scatter = ax.scatter(
            x=plot_df["RichFactor"],
            y=plot_df["Name"],
            s=plot_df["SigCount"] * 22,
            c=plot_df["-log10FDR"],
            cmap="Reds" if preset == "economist" else "YlOrRd",
            alpha=0.88,
            edgecolors=pal.ink,
            linewidths=0.6,
        )

        cbar = plt.colorbar(scatter, ax=ax, shrink=0.65)
        cbar.set_label(f"−log$_{{10}}$({significance_col})")

        sizes = sorted(set([plot_df["SigCount"].min(), plot_df["SigCount"].max()]))
        handles = [plt.scatter([], [], s=s * 22, c=pal.neutral, alpha=0.75, label=str(s)) for s in sizes]
        ax.legend(handles=handles, title="Count", bbox_to_anchor=(1.05, 1), loc="upper left", frameon=False)

        display_comp = comp_name.replace("_vs_", " vs ")
        ax.set_title(f"Pathway Enrichment: {display_comp} ({tissue})", loc="left")
        ax.set_xlabel("Rich Factor (Sig / Bg)")
        ax.invert_yaxis()
        ax.grid(axis="y", alpha=0.4)
        plt.tight_layout(rect=(0, 0, 1, 0.96 if preset == "economist" else 1.0))

        out_path = os.path.join(output_dir, tissue)
        os.makedirs(out_path, exist_ok=True)
        plt.savefig(os.path.join(out_path, f"Pathway_Dot_{comp_name}.png"), dpi=300, bbox_inches="tight",
                    facecolor=pal.surface)
        plt.close(fig)
