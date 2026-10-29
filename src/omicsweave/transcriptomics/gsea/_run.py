"""
GSEA (Gene Set Enrichment Analysis) — core implementation
==========================================================
Runs pre-ranked GSEA using gseapy for all tissue x comparison pairs.
Ranking metric: log2FoldChange x -log10(pvalue)

Config is loaded from experiment_config.py (shared with other modules).
"""

from __future__ import annotations

import os
from typing import Any

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


from ..experiment_config import (
    load_transcriptomics_config, config_arg_help, resolve_comparisons,
    get_group_colors, get_tissue_labels, get_comp_labels, get_plot_colors,
)

# ============================================================================
# Module-level state — populated by apply_config()
# ============================================================================
OUTPUT_BASE: str = ""
DE_DIR: str = ""
PLOT_DIR: str = ""
TISSUES: list[str] = []
COMPARISONS: list[str] = []
TISSUE_LABELS: dict[str, str] = {}
COMP_LABELS: dict[str, str] = {}
PLOT_COLORS: dict[str, Any] = {}

# Gene set libraries for GSEA
GENE_SET_LIBS = [
    'GO_Biological_Process_2025',
    'KEGG_2019_Mouse',
    'WikiPathways_2024_Mouse',
    'Reactome_Pathways_2024',
]


def apply_config(experiment=None, config_path=None, outdir=None) -> dict:
    """Load experiment config and populate module-level variables."""
    global OUTPUT_BASE, DE_DIR, PLOT_DIR, TISSUES, COMPARISONS
    global TISSUE_LABELS, COMP_LABELS, PLOT_COLORS

    cfg = load_transcriptomics_config(
        experiment_id=experiment,
        config_path=config_path,
        outdir_override=outdir,
    )

    OUTPUT_BASE = os.path.abspath(str(cfg.output_dir))
    os.environ["TRANSCRIPTOMICS_OUTDIR"] = OUTPUT_BASE
    DE_DIR = os.path.join(OUTPUT_BASE, "04_de_results")
    PLOT_DIR = os.path.join(OUTPUT_BASE, "05_result_plot")
    TISSUES = list(cfg.tissues)

    comp_tuples = resolve_comparisons(cfg)
    COMPARISONS = [name for _, _, name in comp_tuples]
    TISSUE_LABELS = get_tissue_labels(cfg)
    COMP_LABELS = get_comp_labels(cfg.conditions, comp_tuples)
    PLOT_COLORS = get_plot_colors(cfg)

    return {
        "experiment_id": cfg.experiment_id,
        "config_path": str(cfg.config_path),
        "output_dir": OUTPUT_BASE,
    }


plt.rcParams.update({
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'figure.dpi': 200,
    'savefig.dpi': 200,
    'savefig.bbox': 'tight',
    'font.family': 'sans-serif',
})


def _flush_print(*args, **kwargs):
    print(*args, **kwargs, flush=True)


def _load_de_results(tissue: str, comparison: str) -> pd.DataFrame | None:
    """Load a DE result CSV file (PyDESeq2 or Wilcoxon fallback)."""
    for suffix in ("", "_wilcoxon"):
        path = os.path.join(DE_DIR, tissue, f"de_{comparison}{suffix}.csv")
        if os.path.exists(path):
            df = pd.read_csv(path, index_col=0)
            for col in ('gene_symbol', 'gene'):
                if col in df.columns and df.index.inferred_type != 'string':
                    df = df.set_index(col)
                    break
            return df
    return None


# ============================================================================
# Core GSEA analysis
# ============================================================================
def _generate_gsea_plots():
    """Run pre-ranked GSEA using gseapy."""
    out_dir = os.path.join(PLOT_DIR, '5_gsea')
    os.makedirs(out_dir, exist_ok=True)
    _flush_print("  Running GSEA prerank analysis...")

    import gseapy as gp

    for tissue in TISSUES:
        for comp in COMPARISONS:
            df = _load_de_results(tissue, comp)
            if df is None:
                continue

            # Ranking metric: log2FC x -log10(pvalue)
            rank_df = df[['log2FoldChange', 'pvalue']].dropna().copy()
            rank_df['rank_metric'] = rank_df['log2FoldChange'] * (
                -np.log10(rank_df['pvalue'].clip(lower=1e-300))
            )
            rank_df = rank_df.sort_values('rank_metric', ascending=False)
            rank_series = rank_df['rank_metric']

            # Remove duplicates and zeros
            rank_series = rank_series[~rank_series.index.duplicated(keep='first')]
            rank_series = rank_series[rank_series != 0]

            if len(rank_series) < 100:
                _flush_print(f"    {tissue}/{comp}: too few ranked genes, skipped")
                continue

            gsea_out = os.path.join(out_dir, f'gsea_{tissue}_{comp}')

            try:
                pre_res = gp.prerank(
                    rnk=rank_series,
                    gene_sets=GENE_SET_LIBS,
                    outdir=gsea_out,
                    min_size=15,
                    max_size=500,
                    permutation_num=1000,
                    seed=42,
                    no_plot=True,
                )

                res = pre_res.res2d
                if res.empty:
                    _flush_print(f"    {tissue}/{comp}: no GSEA results")
                    continue

                # Save full results CSV
                res.to_csv(
                    os.path.join(out_dir, f'gsea_{tissue}_{comp}_results.csv'),
                    index=False,
                )

                # Filter FDR < 0.25, then plot top positive/negative NES
                sig_res = res[res['FDR q-val'].astype(float) < 0.25].copy()
                if sig_res.empty:
                    sig_res = res.head(20)

                sig_res['NES'] = sig_res['NES'].astype(float)
                top_pos = sig_res[sig_res['NES'] > 0].head(10)
                top_neg = sig_res[sig_res['NES'] < 0].head(10)
                plot_terms = pd.concat([top_pos, top_neg])

                if len(plot_terms) > 0:
                    plot_terms = plot_terms.sort_values('NES')

                    fig, ax = plt.subplots(
                        figsize=(10, max(4, len(plot_terms) * 0.35 + 1)),
                    )
                    colors = [
                        PLOT_COLORS.get('gsea_positive', '#D32F2F') if x > 0
                        else PLOT_COLORS.get('gsea_negative', '#1565C0')
                        for x in plot_terms['NES']
                    ]
                    ax.barh(
                        range(len(plot_terms)), plot_terms['NES'],
                        color=colors, alpha=0.8,
                    )
                    ax.set_yticks(range(len(plot_terms)))
                    term_names = [
                        str(t)[:55] + '...' if len(str(t)) > 55 else str(t)
                        for t in plot_terms['Term']
                    ]
                    ax.set_yticklabels(term_names, fontsize=7)
                    ax.set_xlabel('Normalized Enrichment Score (NES)')
                    ax.axvline(0, color='black', lw=0.5)

                    tissue_lbl = TISSUE_LABELS.get(tissue, tissue)
                    comp_lbl = COMP_LABELS.get(comp, comp)
                    ax.set_title(f'GSEA \u2014 {tissue_lbl} | {comp_lbl}')

                    fig.tight_layout()
                    fig.savefig(os.path.join(out_dir, f'gsea_{tissue}_{comp}.png'))
                    plt.close(fig)

                n_sig = (res['FDR q-val'].astype(float) < 0.25).sum()
                _flush_print(
                    f"    {tissue}/{comp}: {len(res)} gene sets, "
                    f"{n_sig} significant (FDR<0.25)",
                )

            except Exception as e:
                _flush_print(f"    {tissue}/{comp}: ERROR - {e}")

    _flush_print(f"  GSEA done. Output: {out_dir}")


# ============================================================================
# Public API
# ============================================================================
def run_gsea(*, experiment: str | None = None,
             config_path: str | None = None,
             outdir: str | None = None) -> dict:
    """Run GSEA analysis with config loading.

    Can be called from other scripts:
        from gsea import run_gsea
        run_gsea(experiment='20250425_cimr_wt_mice_dio')

    Or from CLI:
        python run_gsea.py --experiment <id>
    """
    if any(v is not None for v in (experiment, config_path, outdir)):
        cfg_info = apply_config(experiment, config_path, outdir)
    else:
        cfg_info = {}

    _flush_print("\n" + "=" * 60)
    _flush_print("  GSEA (Gene Set Enrichment Analysis)")
    _flush_print("=" * 60)
    if cfg_info:
        _flush_print(f"Experiment:  {cfg_info.get('experiment_id')}")
        _flush_print(f"DE results:  {DE_DIR}")

    _generate_gsea_plots()
    return cfg_info
