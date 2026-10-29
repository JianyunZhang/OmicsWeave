"""
Volcano Plot Generation
========================
Generates volcano plots for each tissue x comparison pair, plus a summary panel.
"""

import os

import numpy as np
import matplotlib.pyplot as plt

from . import _config as cfg


def _plot_single(df, tissue, comparison, output_dir):
    """Draw a single volcano plot for one tissue-comparison pair.

    Returns the output filename.
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    # Prepare data
    plot_df = df[['log2FoldChange', 'padj']].dropna().copy()
    plot_df['-log10padj'] = -np.log10(plot_df['padj'].clip(lower=1e-300))

    # Classify genes as up / down / not significant
    up_cond = (plot_df['padj'] < cfg.PADJ_CUTOFF) & (plot_df['log2FoldChange'] > cfg.LOG2FC_CUTOFF)
    dn_cond = (plot_df['padj'] < cfg.PADJ_CUTOFF) & (plot_df['log2FoldChange'] < -cfg.LOG2FC_CUTOFF)
    ns = ~(up_cond | dn_cond)

    up_color = cfg.PLOT_COLORS.get('volcano_up', '#D32F2F')
    dn_color = cfg.PLOT_COLORS.get('volcano_down', '#1565C0')
    ns_color = cfg.PLOT_COLORS.get('volcano_ns', '#BDBDBD')

    # Non-significant first (background)
    ax.scatter(
        plot_df.loc[ns, 'log2FoldChange'],
        plot_df.loc[ns, '-log10padj'],
        c=ns_color, s=4, alpha=0.4, linewidths=0,
        label='Not significant',
    )
    # Up-regulated
    ax.scatter(
        plot_df.loc[up_cond, 'log2FoldChange'],
        plot_df.loc[up_cond, '-log10padj'],
        c=up_color, s=8, alpha=0.7, linewidths=0,
        label=f'Up-regulated ({up_cond.sum()})',
    )
    # Down-regulated
    ax.scatter(
        plot_df.loc[dn_cond, 'log2FoldChange'],
        plot_df.loc[dn_cond, '-log10padj'],
        c=dn_color, s=8, alpha=0.7, linewidths=0,
        label=f'Down-regulated ({dn_cond.sum()})',
    )

    # Significance threshold lines
    ax.axhline(-np.log10(cfg.PADJ_CUTOFF), ls='--', lw=0.8, color='grey', alpha=0.6)
    ax.axvline(cfg.LOG2FC_CUTOFF, ls='--', lw=0.8, color='grey', alpha=0.6)
    ax.axvline(-cfg.LOG2FC_CUTOFF, ls='--', lw=0.8, color='grey', alpha=0.6)

    # Label significant genes (skip those too close to threshold lines for readability)
    sig = plot_df[
        (plot_df['padj'] < cfg.PADJ_CUTOFF)
        & (plot_df['log2FoldChange'].abs() > cfg.LOG2FC_CUTOFF)
    ]
    if len(sig) > 0:
        padj_thresh_y = -np.log10(cfg.PADJ_CUTOFF)
        fc_margin = cfg.LOG2FC_CUTOFF * 0.3
        padj_margin = (sig['-log10padj'].max() - padj_thresh_y) * 0.08
        label_mask = (
            (sig['log2FoldChange'].abs() > cfg.LOG2FC_CUTOFF + fc_margin)
            & (sig['-log10padj'] > padj_thresh_y + padj_margin)
        )
        label_genes = sig[label_mask]
        label_fontsize = 7 if len(label_genes) <= 30 else max(4, 7 - (len(label_genes) - 30) // 30)
        try:
            from adjustText import adjust_text
            texts = [
                ax.text(r['log2FoldChange'], r['-log10padj'], g,
                        fontsize=label_fontsize, ha='center', style='italic')
                for g, r in label_genes.iterrows()
            ]
            adjust_text(texts, arrowprops=dict(arrowstyle='-', color='grey', lw=0.5), ax=ax)
        except ImportError:
            for gene, row in label_genes.iterrows():
                ax.annotate(
                    gene,
                    (row['log2FoldChange'], row['-log10padj']),
                    fontsize=label_fontsize, ha='center', style='italic',
                )

    # Title and axis labels
    tissue_label = cfg.TISSUE_LABELS.get(tissue, tissue)
    comp_label = cfg.COMP_LABELS.get(comparison, comparison)
    ax.set_title(f'{tissue_label} \u2014 {comp_label}')
    ax.set_xlabel('log\u2082(Fold Change)')
    ax.set_ylabel('-log\u2081\u2080(adjusted p-value)')
    ax.legend(loc='upper right', fontsize=8, framealpha=0.9)

    # Symmetric x-axis
    xlim = max(abs(ax.get_xlim()[0]), abs(ax.get_xlim()[1]))
    ax.set_xlim(-xlim, xlim)

    fname = os.path.join(output_dir, f'volcano_{tissue}_{comparison}.png')
    fig.savefig(fname)
    plt.close(fig)
    return fname


def generate_all():
    """Generate volcano plots for all tissue x comparison pairs + summary panel."""
    out_dir = os.path.join(cfg.PLOT_DIR, '1_volcano')
    os.makedirs(out_dir, exist_ok=True)
    cfg.flush_print("\n[1/5] Generating Volcano Plots...")

    count = 0
    for tissue in cfg.TISSUES:
        for comp in cfg.COMPARISONS:
            df = cfg.load_de_results(tissue, comp)
            if df is None:
                continue
            fname = _plot_single(df, tissue, comp, out_dir)
            n_up = ((df['padj'] < cfg.PADJ_CUTOFF) & (df['log2FoldChange'] > cfg.LOG2FC_CUTOFF)).sum()
            n_down = ((df['padj'] < cfg.PADJ_CUTOFF) & (df['log2FoldChange'] < -cfg.LOG2FC_CUTOFF)).sum()
            cfg.flush_print(f"  {tissue}/{comp}: {n_up}\u2191 {n_down}\u2193 -> {os.path.basename(fname)}")
            count += 1

    # Summary multi-panel figure
    cfg.flush_print("  Creating summary panel...")
    n_t = len(cfg.TISSUES)
    n_c = len(cfg.COMPARISONS)
    fig, axes = plt.subplots(n_t, n_c, figsize=(5 * n_c, 4 * n_t), squeeze=False)

    for i, tissue in enumerate(cfg.TISSUES):
        for j, comp in enumerate(cfg.COMPARISONS):
            ax = axes[i, j]
            df = cfg.load_de_results(tissue, comp)
            if df is None:
                ax.text(0.5, 0.5, 'N/A', ha='center', va='center', transform=ax.transAxes)
                ax.set_title(
                    f'{cfg.TISSUE_LABELS.get(tissue, tissue)}\n{cfg.COMP_LABELS.get(comp, comp)}',
                    fontsize=9,
                )
                continue

            plot_df = df[['log2FoldChange', 'padj']].dropna().copy()
            plot_df['-log10p'] = -np.log10(plot_df['padj'].clip(lower=1e-300))

            up = (plot_df['padj'] < cfg.PADJ_CUTOFF) & (plot_df['log2FoldChange'] > cfg.LOG2FC_CUTOFF)
            down = (plot_df['padj'] < cfg.PADJ_CUTOFF) & (plot_df['log2FoldChange'] < -cfg.LOG2FC_CUTOFF)
            ns = ~(up | down)

            ax.scatter(
                plot_df.loc[ns, 'log2FoldChange'], plot_df.loc[ns, '-log10p'],
                c=cfg.PLOT_COLORS.get('volcano_ns', '#BDBDBD'), s=2, alpha=0.3, linewidths=0,
            )
            ax.scatter(
                plot_df.loc[up, 'log2FoldChange'], plot_df.loc[up, '-log10p'],
                c=cfg.PLOT_COLORS.get('volcano_up', '#D32F2F'), s=3, alpha=0.5, linewidths=0,
            )
            ax.scatter(
                plot_df.loc[down, 'log2FoldChange'], plot_df.loc[down, '-log10p'],
                c=cfg.PLOT_COLORS.get('volcano_down', '#1565C0'), s=3, alpha=0.5, linewidths=0,
            )

            ax.axhline(-np.log10(cfg.PADJ_CUTOFF), ls='--', lw=0.5, color='grey')
            ax.axvline(cfg.LOG2FC_CUTOFF, ls='--', lw=0.5, color='grey')
            ax.axvline(-cfg.LOG2FC_CUTOFF, ls='--', lw=0.5, color='grey')

            ax.set_title(
                f'{cfg.TISSUE_LABELS.get(tissue, tissue)}\n{cfg.COMP_LABELS.get(comp, comp)}',
                fontsize=9,
            )
            if j == 0:
                ax.set_ylabel('-log\u2081\u2080(padj)', fontsize=8)
            if i == n_t - 1:
                ax.set_xlabel('log\u2082FC', fontsize=8)
            ax.tick_params(labelsize=7)

    fig.suptitle('Volcano Plots \u2014 All Tissues \u00d7 Comparisons', fontsize=14, y=1.01)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'volcano_summary_panel.png'))
    plt.close(fig)

    cfg.flush_print(f"  Total: {count} volcano plots + 1 summary panel")
