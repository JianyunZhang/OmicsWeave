"""
Heatmap Generation
===================
Generates clustered heatmaps of top DEGs per tissue x comparison.
Samples are filtered and sorted by the two groups in each comparison.

绘图管线（顺序不可调换）:

1. **文库大小归一化** —— median-of-ratios size factor，用该比较的
   *全部* 基因估计后应用到显著基因子集。必须在 ``log2(x+1)`` 之前做：
   测序深度差异是乘性的，先取对数再减就把量纲弄错了。
   size factor 只能用全部基因估计——median-of-ratios 的前提是
   「多数基因不变」，而显著基因按定义正是变了的那批。
2. **按 DE 证据选 Top N** —— padj 升序、|log2FC| 降序。排序依据与任何
   显示变换无关，这是刻意的：曾经的实现在 z-score *之后* 取行方差
   ``nlargest``，而 z-score 使每行方差恒为 1，排序退化成「取索引前 N 行」。
3. **Z-score 仅用于显示** —— 按基因（行）标准化只是为了配色可比，
   它消不掉按样本（列）的深度偏移，所以第 1 步不能省。
"""

import os

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from . import _config as cfg

from ..normalization import median_of_ratios_size_factors


def generate_all():
    """Generate heatmaps of top DEGs per tissue x comparison."""
    out_dir = os.path.join(cfg.PLOT_DIR, '2_heatmap')
    os.makedirs(out_dir, exist_ok=True)
    cfg.flush_print("\n[2/5] Generating Heatmaps...")

    # Load gene counts matrix
    counts_path = os.path.join(cfg.DE_DIR, 'gene_counts_matrix.csv')
    counts_df = pd.read_csv(counts_path, index_col=0, low_memory=False)
    counts_df = counts_df.apply(pd.to_numeric, errors='coerce').fillna(0)
    counts_df = counts_df.loc[:, counts_df.sum() > 0]

    count = 0
    for tissue in cfg.TISSUES:
        tissue_samples = [s for s in counts_df.columns if s.split('_', 1)[1] == tissue]
        if len(tissue_samples) < 4:
            cfg.flush_print(f"  {tissue}: skipped (too few samples)")
            continue

        for num, denom, comp_name in cfg.COMP_TUPLES:
            de_df = cfg.load_de_results(tissue, comp_name)
            if de_df is None:
                continue

            # Filter samples to only the two groups in this comparison
            comp_samples = [
                s for s in tissue_samples
                if cfg.MOUSE_GROUPS.get(s.split('_')[0]) in (num, denom)
            ]
            if len(comp_samples) < 2:
                continue

            # Get significant genes for this comparison
            sig = de_df[
                (de_df['padj'] < cfg.PADJ_CUTOFF)
                & (de_df['log2FoldChange'].abs() > cfg.LOG2FC_CUTOFF)
            ]
            if len(sig) == 0:
                continue

            available_genes = [g for g in sig.index if g in counts_df.index]
            if len(available_genes) == 0:
                continue

            # 1) 文库大小归一化 —— size factor 用该比较的全部基因估计，
            #    再应用到显著基因子集（见模块 docstring）。
            comp_counts = counts_df[comp_samples]
            size_factors = median_of_ratios_size_factors(comp_counts)
            sub_log = np.log2(
                comp_counts.loc[available_genes].div(size_factors, axis=1) + 1
            )

            # 组内零变异的基因 z-score 会得到 0/0，先剔除再选 Top N，
            # 否则它们会占掉名额最后又被 dropna 丢掉。
            sub_log = sub_log[sub_log.std(axis=1) > 0]
            if len(sub_log) < 2:
                continue

            # 2) 按 DE 证据选 Top N（padj 升序，|log2FC| 降序）。
            #    不用行方差：那必须在 z-score 之前算，且 z-score 后恒为 1。
            max_genes = 80
            ranked = (
                sig.loc[sub_log.index]
                .assign(_abs_lfc=lambda d: d['log2FoldChange'].abs())
                .sort_values(['padj', '_abs_lfc'], ascending=[True, False], kind='stable')
            )
            sub_log = sub_log.loc[ranked.index[:max_genes]]

            # 3) Z-score 仅用于显示配色
            sub_z = sub_log.sub(sub_log.mean(axis=1), axis=0).div(sub_log.std(axis=1), axis=0)

            # Sort samples: control group first, then treatment
            sample_groups = []
            for s in comp_samples:
                mouse_id = s.split('_')[0]
                grp = cfg.MOUSE_GROUPS.get(mouse_id, 'Unknown')
                sample_groups.append((grp, s))
            sample_groups.sort(key=lambda x: (x[0] != denom, x[0], x[1]))
            sorted_samples = [s for _, s in sample_groups]

            # Column annotation colors (only the two groups in this comparison)
            comp_colors = {num: cfg.GROUP_COLORS.get(num, '#D32F2F'),
                           denom: cfg.GROUP_COLORS.get(denom, '#1565C0')}
            col_colors = [comp_colors.get(grp, '#999999') for grp, _ in sample_groups]

            # Reorder data columns
            sub_z = sub_z[sorted_samples]

            # Plot
            n_genes = len(sub_z)
            figheight = max(8, n_genes * 0.15 + 2)
            with plt.rc_context({'font.family': 'Arial'}):
                g = sns.clustermap(
                    sub_z,
                    cmap=cfg.PLOT_COLORS.get('heatmap_cmap', 'RdBu_r'),
                    center=0,
                    vmin=cfg.PLOT_COLORS.get('heatmap_vmin', -2),
                    vmax=cfg.PLOT_COLORS.get('heatmap_vmax', 2),
                    col_colors=col_colors,
                    figsize=(max(6, len(sorted_samples) * 0.8 + 2), figheight),
                    xticklabels=True,
                    yticklabels=True,
                    linewidths=0,
                    dendrogram_ratio=(0.12, 0.08),
                    cbar_pos=(0.02, 0.82, 0.03, 0.15),
                    col_cluster=False,
                )
                g.ax_heatmap.set_xlabel('')
                g.ax_heatmap.set_ylabel('')
                g.ax_heatmap.yaxis.tick_right()

                ytick_fontsize = max(6, min(10, 240 // n_genes))
                xtick_fontsize = max(5, min(9, 120 // len(sorted_samples)))
                g.ax_heatmap.tick_params(axis='y', labelsize=ytick_fontsize, length=0)
                g.ax_heatmap.tick_params(axis='x', labelsize=xtick_fontsize,
                                         rotation=270, length=0)

                if hasattr(g, 'ax_col_colors') and g.ax_col_colors is not None:
                    pos_heat = g.ax_heatmap.get_position()
                    pos_col = g.ax_col_colors.get_position()
                    g.ax_col_colors.set_position([
                        pos_col.x0, pos_heat.y1 + 0.012,
                        pos_col.width, pos_col.height,
                    ])

                tissue_label = cfg.TISSUE_LABELS.get(tissue, tissue)
                comp_label = cfg.COMP_LABELS.get(comp_name, comp_name)
                g.fig.suptitle(
                    f'{tissue_label} — {comp_label}\nTop DEGs Heatmap (Z-score, {n_genes} genes)',
                    fontsize=12, y=1.02, fontfamily='Arial',
                )

                legend_elements = [
                    Patch(facecolor=comp_colors[denom], label=denom),
                    Patch(facecolor=comp_colors[num], label=num),
                ]
                g.ax_heatmap.legend(
                    handles=legend_elements, loc='upper left',
                    bbox_to_anchor=(1.02, 1.15), fontsize=8, frameon=False,
                )

                fname = os.path.join(out_dir, f'heatmap_{tissue}_{comp_name}.png')
                g.savefig(fname)
                plt.close(g.fig)
            cfg.flush_print(f"  {tissue}/{comp_name}: {n_genes} genes -> {os.path.basename(fname)}")
            count += 1

    cfg.flush_print(f"  Heatmaps done ({count} plots).")