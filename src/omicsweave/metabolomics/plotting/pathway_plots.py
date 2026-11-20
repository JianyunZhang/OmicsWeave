"""
通路富集结果可视化：柱状图、气泡图。
"""
from __future__ import annotations

import os

import matplotlib.pyplot as plt
import numpy as np

from omicsweave.visualization.palettes import NATURE_EDITORIAL
from omicsweave.visualization.theme import themed

from . import _config  # noqa: F401


def plot_pathway_bar(
    res_df: "pd.DataFrame",
    comp_name: str,
    tissue: str,
    output_dir: str,
):
    """通路富集柱状图 (Top 20, -log10 FDR，采用 nature_editorial 单色编辑柱图，300 dpi)。"""
    if res_df.empty:
        return
    import pandas as pd

    significance_col = "FDR" if "FDR" in res_df.columns else "P-value"
    plot_df = res_df.nsmallest(20, significance_col).copy()
    plot_df["-log10FDR"] = -np.log10(
        pd.to_numeric(plot_df[significance_col], errors="coerce").clip(lower=1e-300)
    )
    plot_df = plot_df.sort_values("-log10FDR", ascending=True)

    with themed(preset="nature_editorial", save_dpi=300):
        fig, ax = plt.subplots(figsize=(8.5, max(4.5, len(plot_df) * 0.34)))
        ax.barh(plot_df["Name"], plot_df["-log10FDR"], color=NATURE_EDITORIAL.down, height=0.68, edgecolor="none")
        display_comp = comp_name.replace("_vs_", " vs ")
        ax.set_title(f"Pathway Enrichment: {display_comp} ({tissue})", loc="left")
        ax.set_xlabel(f"−log$_{{10}}$({significance_col})")
        ax.set_ylabel("")
        ax.grid(axis="x", alpha=0.4)
        plt.tight_layout()

        out_path = os.path.join(output_dir, tissue)
        os.makedirs(out_path, exist_ok=True)
        plt.savefig(os.path.join(out_path, f"Pathway_Bar_{comp_name}.png"), dpi=300, bbox_inches="tight")
        plt.close(fig)


def plot_pathway_dot(
    res_df: "pd.DataFrame",
    comp_name: str,
    tissue: str,
    output_dir: str,
):
    """通路富集气泡图 (Rich Factor vs Pathway，采用 nature_editorial，300 dpi)。"""
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

    with themed(preset="nature_editorial", save_dpi=300):
        fig, ax = plt.subplots(figsize=(8.5, max(4.5, len(plot_df) * 0.36)))
        scatter = ax.scatter(
            x=plot_df["RichFactor"],
            y=plot_df["Name"],
            s=plot_df["SigCount"] * 22,
            c=plot_df["-log10FDR"],
            cmap="YlOrRd",
            alpha=0.88,
            edgecolors=NATURE_EDITORIAL.ink,
            linewidths=0.6,
        )

        cbar = plt.colorbar(scatter, ax=ax, shrink=0.65)
        cbar.set_label(f"−log$_{{10}}$({significance_col})")

        sizes = sorted(set([plot_df["SigCount"].min(), plot_df["SigCount"].max()]))
        handles = [plt.scatter([], [], s=s * 22, c=NATURE_EDITORIAL.neutral, alpha=0.75, label=str(s)) for s in sizes]
        ax.legend(handles=handles, title="Count", bbox_to_anchor=(1.05, 1), loc="upper left", frameon=False)

        display_comp = comp_name.replace("_vs_", " vs ")
        ax.set_title(f"Pathway Enrichment: {display_comp} ({tissue})", loc="left")
        ax.set_xlabel("Rich Factor (Sig / Bg)")
        ax.invert_yaxis()
        ax.grid(True, alpha=0.35)

        out_path = os.path.join(output_dir, tissue)
        os.makedirs(out_path, exist_ok=True)
        plt.savefig(os.path.join(out_path, f"Pathway_Dot_{comp_name}.png"), dpi=300, bbox_inches="tight")
        plt.close(fig)
