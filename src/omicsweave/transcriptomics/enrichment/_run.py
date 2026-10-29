"""
GO/KEGG Pathway Enrichment Analysis — core implementation
==========================================================
Runs GO/KEGG/WikiPathways/Reactome enrichment via gseapy.enrichr
for all tissue x comparison pairs defined in the experiment config.

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

PADJ_CUTOFF: float = 0.05
LOG2FC_CUTOFF: float = 1.0

# Gene set libraries for enrichment analysis
GENE_SET_LIBS = [
    'GO_Biological_Process_2025',
    'GO_Molecular_Function_2025',
    'GO_Cellular_Component_2025',
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
# Core enrichment analysis
# ============================================================================
def _generate_enrichment_plots():
    """Run GO/KEGG enrichment analysis using gseapy.enrichr."""
    out_dir = os.path.join(PLOT_DIR, '4_enrichment')
    os.makedirs(out_dir, exist_ok=True)
    _flush_print("  Running Enrichr enrichment analysis...")

    import gseapy as gp

    for tissue in TISSUES:
        for comp in COMPARISONS:
            df = _load_de_results(tissue, comp)
            if df is None:
                continue

            sig = df[
                (df['padj'] < PADJ_CUTOFF)
                & (df['log2FoldChange'].abs() > LOG2FC_CUTOFF)
            ]
            if len(sig) < 5:
                _flush_print(f"    {tissue}/{comp}: too few DEGs ({len(sig)}), skipped")
                continue

            up_genes = sig[sig['log2FoldChange'] > 0].index.tolist()
            down_genes = sig[sig['log2FoldChange'] < 0].index.tolist()

            for direction, gene_list in [('up', up_genes), ('down', down_genes)]:
                if len(gene_list) < 3:
                    continue

                try:
                    enr = gp.enrichr(
                        gene_list=gene_list,
                        gene_sets=GENE_SET_LIBS,
                        organism='mouse',
                        outdir=None,
                        no_plot=True,
                    )
                    res = enr.results
                    if res.empty:
                        continue

                    # Save full results CSV
                    res_path = os.path.join(
                        out_dir, f'enrichr_{tissue}_{comp}_{direction}.csv',
                    )
                    res.to_csv(res_path, index=False)

                    # Plot top 15 terms
                    sig_terms = res[res['Adjusted P-value'] < 0.05].head(15)
                    if sig_terms.empty:
                        sig_terms = res.head(15)

                    if len(sig_terms) > 0:
                        fig, ax = plt.subplots(
                            figsize=(10, max(4, len(sig_terms) * 0.35 + 1)),
                        )
                        sig_terms = sig_terms.iloc[::-1]
                        colors = -np.log10(sig_terms['Adjusted P-value'].clip(lower=1e-20))
                        ax.barh(
                            range(len(sig_terms)),
                            -np.log10(sig_terms['Adjusted P-value'].clip(lower=1e-20)),
                            color=plt.get_cmap(
                                PLOT_COLORS.get('enrichment_cmap', 'RdYlBu_r')
                            )(colors / colors.max()) if colors.max() > 0 else '#4CAF50',
                        )
                        ax.set_yticks(range(len(sig_terms)))
                        term_names = [
                            t[:60] + '...' if len(t) > 60 else t
                            for t in sig_terms['Term']
                        ]
                        ax.set_yticklabels(term_names, fontsize=7)
                        ax.set_xlabel('-log\u2081\u2080(Adjusted P-value)')

                        tissue_lbl = TISSUE_LABELS.get(tissue, tissue)
                        comp_lbl = COMP_LABELS.get(comp, comp)
                        dir_lbl = 'Up-regulated' if direction == 'up' else 'Down-regulated'
                        ax.set_title(
                            f'{tissue_lbl} | {comp_lbl} | {dir_lbl} ({len(gene_list)} genes)',
                        )

                        fig.tight_layout()
                        fig.savefig(os.path.join(
                            out_dir, f'enrichr_{tissue}_{comp}_{direction}.png',
                        ))
                        plt.close(fig)

                    _flush_print(
                        f"    {tissue}/{comp}/{direction}: "
                        f"{len(gene_list)} genes, {len(res)} terms",
                    )
                except Exception as e:
                    _flush_print(f"    {tissue}/{comp}/{direction}: ERROR - {e}")

    _flush_print(f"  Enrichment done. Output: {out_dir}")


# ============================================================================
# Public API
# ============================================================================
def run_enrichment(*, experiment: str | None = None,
                   config_path: str | None = None,
                   outdir: str | None = None) -> dict:
    """Run enrichment analysis with config loading.

    Can be called from other scripts:
        from enrichment import run_enrichment
        run_enrichment(experiment='20250425_cimr_wt_mice_dio')

    Or from CLI:
        python run_enrichment.py --experiment <id>
    """
    if any(v is not None for v in (experiment, config_path, outdir)):
        cfg_info = apply_config(experiment, config_path, outdir)
    else:
        cfg_info = {}

    _flush_print("\n" + "=" * 60)
    _flush_print("  GO/KEGG Pathway Enrichment Analysis")
    _flush_print("=" * 60)
    if cfg_info:
        _flush_print(f"Experiment:  {cfg_info.get('experiment_id')}")
        _flush_print(f"DE results:  {DE_DIR}")

    _generate_enrichment_plots()
    return cfg_info
