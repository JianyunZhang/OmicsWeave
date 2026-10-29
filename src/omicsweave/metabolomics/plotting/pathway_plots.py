"""
通路富集结果可视化：柱状图、气泡图。
"""
from __future__ import annotations

import os

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from . import _config  # noqa: F401


def plot_pathway_bar(
    res_df: "pd.DataFrame",
    comp_name: str,
    tissue: str,
    output_dir: str,
):
    """通路富集柱状图 (Top 20, -log10 FDR)。"""
    if res_df.empty:
        return
    import pandas as pd

    significance_col = "FDR" if "FDR" in res_df.columns else "P-value"
    plot_df = res_df.nsmallest(20, significance_col).copy()
    plot_df["-log10FDR"] = -np.log10(
        pd.to_numeric(plot_df[significance_col], errors="coerce").clip(lower=1e-300)
    )

    fig, ax = plt.subplots(figsize=(10, 8))
    sns.barplot(data=plot_df, x="-log10FDR", y="Name", hue="Name", palette="viridis", legend=False, ax=ax)

    display_comp = comp_name.replace("_vs_", " vs ")
    ax.set_title(f"通路富集 (柱状图): {display_comp} ({tissue})")
    ax.set_xlabel(f"-log10({significance_col})")
    ax.set_ylabel("")
    plt.tight_layout()

    out_path = os.path.join(output_dir, tissue)
    os.makedirs(out_path, exist_ok=True)
    plt.savefig(os.path.join(out_path, f"Pathway_Bar_{comp_name}.png"), dpi=150, bbox_inches="tight")
    plt.close()


def plot_pathway_dot(
    res_df: "pd.DataFrame",
    comp_name: str,
    tissue: str,
    output_dir: str,
):
    """通路富集气泡图 (Rich Factor vs Pathway)。"""
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

    fig, ax = plt.subplots(figsize=(10, 8))
    scatter = ax.scatter(
        x=plot_df["RichFactor"],
        y=plot_df["Name"],
        s=plot_df["SigCount"] * 20,
        c=plot_df["-log10FDR"],
        cmap="YlOrRd",
        alpha=0.8,
        edgecolors="grey",
    )

    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label(f"-log10({significance_col})")

    # 大小图例
    sizes = sorted(set([plot_df["SigCount"].min(), plot_df["SigCount"].max()]))
    handles = [plt.scatter([], [], s=s * 20, c="grey", alpha=0.6, label=str(s)) for s in sizes]
    ax.legend(handles=handles, title="Count", bbox_to_anchor=(1.05, 1), loc="upper left")

    display_comp = comp_name.replace("_vs_", " vs ")
    ax.set_title(f"通路富集 (气泡图): {display_comp} ({tissue})")
    ax.set_xlabel("Rich Factor (Sig/Bg)")
    ax.invert_yaxis()
    ax.grid(True, alpha=0.3)

    out_path = os.path.join(output_dir, tissue)
    os.makedirs(out_path, exist_ok=True)
    plt.savefig(os.path.join(out_path, f"Pathway_Dot_{comp_name}.png"), dpi=150, bbox_inches="tight")
    plt.close()
