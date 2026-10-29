"""
PCA plot module for proteomics.

Generates PCA scatter plots from expression matrices.
"""
from __future__ import annotations

import os

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA

from . import _config as cfg


def plot_single(
    df: pd.DataFrame,
    group_map: dict[str, str],
    output_path: str,
    *,
    group_colors: dict[str, str] | None = None,
) -> None:
    """Generate PCA scatter plot.

    Args:
        df: Expression matrix (proteins as index, samples as columns).
        group_map: Sample name -> group mapping.
        output_path: Path to save the PCA plot.
        group_colors: Optional group -> color mapping.
    """
    colors = group_colors or cfg.GROUP_COLORS
    X = df.T
    pca = PCA(n_components=2)
    components = pca.fit_transform(X)
    pca_df = pd.DataFrame(data=components, columns=["PC1", "PC2"])
    pca_df["Sample"] = X.index
    pca_df["Group"] = pca_df["Sample"].map(group_map)
    var = pca.explained_variance_ratio_

    plt.figure(figsize=(10, 8))
    groups_present = sorted(pca_df["Group"].unique())
    palette = {g: colors.get(g, "#888888") for g in groups_present}
    sns.scatterplot(
        data=pca_df, x="PC1", y="PC2",
        hue="Group", s=100, palette=palette,
    )
    for i, txt in enumerate(pca_df["Sample"]):
        plt.annotate(
            txt, (pca_df.PC1[i], pca_df.PC2[i]),
            xytext=(5, 5), textcoords="offset points", fontsize=8,
        )
    plt.title("PCA of Proteomics Data")
    plt.xlabel(f"PC1 ({var[0]:.2%})")
    plt.ylabel(f"PC2 ({var[1]:.2%})")
    plt.legend(title="Group")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


def generate_all() -> None:
    """Regenerate PCA plots for all tissues from saved imputed matrices."""
    cfg.flush_print("Generating PCA plots...")
    from ..experiment_config import get_sample_to_group_map

    for tissue in cfg.TISSUES:
        imp_path = os.path.join(cfg.DP_BASE, tissue, "imputed_matrix.csv")
        if not os.path.exists(imp_path):
            cfg.flush_print(f"  [{tissue}] No imputed matrix found, skipping PCA.")
            continue

        imp_df = pd.read_csv(imp_path, index_col=0)
        group_map = get_sample_to_group_map(imp_df.columns, cfg.MOUSE_GROUPS)

        tissue_dir = os.path.join(cfg.PLOT_BASE, tissue)
        os.makedirs(tissue_dir, exist_ok=True)
        pca_path = os.path.join(tissue_dir, "pca_plot.png")
        plot_single(imp_df, group_map, pca_path)
        cfg.flush_print(f"  [{tissue}] PCA plot regenerated.")

    cfg.flush_print("  PCA plots done.")