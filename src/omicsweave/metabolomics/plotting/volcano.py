"""
火山图绘制。
"""
from __future__ import annotations

import os

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from omicsweave.visualization.palettes import METRO_ATLAS

from . import _config  # noqa: F401


def plot_volcano(
    df,
    comp_name: str,
    tissue: str,
    output_dir: str,
    pval_cutoff: float = 0.05,
    fdr_cutoff: float | None = None,
    log2fc_cutoff: float = 0.58,
    colors: dict | None = None,
):
    """
    绘制火山图。

    Parameters
    ----------
    df : DataFrame — 必须包含 Pval_{comp}, Log2FC_{comp}, Compounds 列
    comp_name : str — 比较名称，如 "HFD_vs_LFD"
    tissue : str — 组织名
    output_dir : str — 输出目录
    """
    p_col = f"Pval_{comp_name}"
    fdr_col = f"FDR_{comp_name}"
    fc_col = f"Log2FC_{comp_name}"
    sig_col = f"Significant_{comp_name}"

    if p_col not in df.columns or fc_col not in df.columns:
        print(f"    跳过火山图 {comp_name}: 缺少列 {p_col} 或 {fc_col}")
        return

    c_up = (colors or {}).get("volcano_up", METRO_ATLAS.up)
    c_down = (colors or {}).get("volcano_down", METRO_ATLAS.down)
    c_ns = (colors or {}).get("volcano_ns", METRO_ATLAS.faint)

    x = df[fc_col]
    significance_metric = fdr_col if fdr_col in df.columns else p_col
    threshold = fdr_cutoff if significance_metric == fdr_col and fdr_cutoff is not None else pval_cutoff
    adjusted = df[significance_metric].clip(lower=1e-300)
    y = -np.log10(adjusted)

    # 显著性分类
    if sig_col in df.columns:
        sig_mask = df[sig_col].fillna(False).astype(bool)
    else:
        sig_mask = (df[significance_metric] < threshold) & (df[fc_col].abs() > log2fc_cutoff)
    point_colors = np.where(sig_mask, np.where(df[fc_col] > 0, c_up, c_down), c_ns)

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.scatter(x, y, c=point_colors.tolist(), alpha=0.6, s=15)

    # 标注 Top 10 显著代谢物
    sig_df = df[sig_mask].copy()
    sig_df["_abs_fc"] = sig_df[fc_col].abs()
    top_labels = sig_df.nsmallest(10, significance_metric)

    for _, row in top_labels.iterrows():
        name = str(row.get("Compounds", row.get("Index", "")))
        if len(name) > 25:
            name = name[:22] + "..."
        ax.text(
            row[fc_col],
            -np.log10(max(row[significance_metric], 1e-300)),
            name,
            fontsize=7,
            ha="right",
            va="bottom",
        )

    # 阈值线
    ax.axhline(-np.log10(threshold), color="grey", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.axvline(log2fc_cutoff, color="grey", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.axvline(-log2fc_cutoff, color="grey", linestyle="--", linewidth=0.8, alpha=0.5)

    display_comp = comp_name.replace("_vs_", " vs ")
    ax.set_title(f"火山图: {display_comp} ({tissue})")
    ax.set_xlabel(f"Log2 Fold Change ({display_comp})")
    ax.set_ylabel("-Log10 FDR" if significance_metric == fdr_col else "-Log10 P-value")

    out_path = os.path.join(output_dir, tissue)
    os.makedirs(out_path, exist_ok=True)
    plt.savefig(os.path.join(out_path, f"Volcano_{comp_name}.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"    火山图保存: {out_path}/Volcano_{comp_name}.png")
