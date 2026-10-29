"""
Finish GSEA Analysis
====================
Completes GSEA for missing tissue/comparison pairs.
Skips existing results in analysis_output/05_result_plot/5_gsea.

All experiment-specific values (tissues, comparisons, labels) are loaded
from the experiment YAML config. No hardcoded experiment data.

Usage:
  python finish_gsea.py --experiment <id>
  python finish_gsea.py --config <path/to/config.yaml>
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import gseapy as gp

from .experiment_config import (
    load_transcriptomics_config, config_arg_help, resolve_comparisons,
    get_tissue_labels, get_comp_labels,
)

# ============================================================================
# CONFIGURATION — all values loaded from experiment YAML config.
# ============================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Placeholders — will be overwritten by _apply_config()
OUTPUT_BASE = ""
DE_DIR = ""
PLOT_DIR = ""
GSEA_DIR = ""
TISSUES: list[str] = []
COMPARISONS: list[str] = []
TISSUE_LABELS: dict[str, str] = {}
COMP_LABELS: dict[str, str] = {}

# Mouse-specific gene set libraries (non-experiment-specific)
GENE_SET_LIBS = [
    'GO_Biological_Process_2025',
    'KEGG_2019_Mouse',
    'WikiPathways_2024_Mouse',
    'Reactome_Pathways_2024',
]


def _apply_config(experiment: str | None, config_path: str | None, outdir: str | None) -> dict:
    cfg = load_transcriptomics_config(
        experiment_id=experiment,
        config_path=config_path,
        outdir_override=outdir,
    )
    global OUTPUT_BASE, DE_DIR, PLOT_DIR, GSEA_DIR, TISSUES, COMPARISONS
    global TISSUE_LABELS, COMP_LABELS

    OUTPUT_BASE = os.path.abspath(str(cfg.output_dir))
    os.environ["TRANSCRIPTOMICS_OUTDIR"] = OUTPUT_BASE
    DE_DIR = os.path.join(OUTPUT_BASE, "04_de_results")
    PLOT_DIR = os.path.join(OUTPUT_BASE, "05_result_plot")
    GSEA_DIR = os.path.join(PLOT_DIR, "5_gsea")
    TISSUES = list(cfg.tissues)

    # Auto-generate comparisons and labels from config
    comp_tuples = resolve_comparisons(cfg)
    COMPARISONS = [name for _, _, name in comp_tuples]
    TISSUE_LABELS = get_tissue_labels(cfg)
    COMP_LABELS = get_comp_labels(cfg.conditions, comp_tuples)

    return {"experiment_id": cfg.experiment_id, "config_path": str(cfg.config_path), "output_dir": OUTPUT_BASE}


def flush_print(*args, **kwargs):
    print(*args, **kwargs, flush=True)


def load_de_results(tissue, comparison):
    path = os.path.join(DE_DIR, tissue, f"de_{comparison}.csv")
    if not os.path.exists(path):
        path = os.path.join(DE_DIR, tissue, f"de_{comparison}_wilcoxon.csv")
    df = pd.read_csv(path, index_col=0)
    return df


def run_gsea():
    os.makedirs(GSEA_DIR, exist_ok=True)
    flush_print("Starting GSEA completion...")

    for tissue in TISSUES:
        for comp in COMPARISONS:
            # Check if output exists
            res_csv = os.path.join(GSEA_DIR, f'gsea_{tissue}_{comp}_results.csv')
            if os.path.exists(res_csv):
                flush_print(f"  Skipping {tissue}/{comp}: already done")
                continue

            flush_print(f"  Processing {tissue}/{comp}...")

            df = load_de_results(tissue, comp)
            if df is None:
                flush_print(f"    No DE results found for {tissue}/{comp}")
                continue

            # Create ranked gene list: gene -> log2FC * -log10(pvalue)
            rank_df = df[['log2FoldChange', 'pvalue']].dropna().copy()
            rank_df['rank_metric'] = rank_df['log2FoldChange'] * (-np.log10(
                rank_df['pvalue'].clip(lower=1e-300)))
            rank_df = rank_df.sort_values('rank_metric', ascending=False)
            rank_series = rank_df['rank_metric']

            # Remove duplicates and zero
            rank_series = rank_series[~rank_series.index.duplicated(keep='first')]
            rank_series = rank_series[rank_series != 0]

            if len(rank_series) < 100:
                flush_print(f"    Too few ranked genes ({len(rank_series)}), skipping")
                continue

            gsea_out = os.path.join(GSEA_DIR, f'gsea_{tissue}_{comp}')

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
                    verbose=True
                )

                res = pre_res.res2d
                if res.empty:
                    flush_print(f"    No GSEA results found")
                    continue

                # Save results table
                res.to_csv(res_csv, index=False)

                # Plot top enriched & depleted pathways
                sig_res = res[res['FDR q-val'].astype(float) < 0.25].copy()
                if sig_res.empty:
                    sig_res = res.head(20)

                sig_res['NES'] = sig_res['NES'].astype(float)

                top_pos = sig_res[sig_res['NES'] > 0].sort_values('NES', ascending=False).head(10)
                top_neg = sig_res[sig_res['NES'] < 0].sort_values('NES', ascending=True).head(10)
                plot_terms = pd.concat([top_pos, top_neg]).sort_values('NES')

                if len(plot_terms) > 0:
                    fig, ax = plt.subplots(figsize=(10, max(4, len(plot_terms) * 0.35 + 1)))
                    colors = ['#D32F2F' if x > 0 else '#1565C0' for x in plot_terms['NES']]
                    ax.barh(range(len(plot_terms)), plot_terms['NES'], color=colors, alpha=0.8)
                    ax.set_yticks(range(len(plot_terms)))

                    term_names = []
                    for t in plot_terms['Term']:
                        t_str = str(t)
                        if ' (GO:' in t_str:
                            t_str = t_str.split(' (GO:')[0]
                        if len(t_str) > 55:
                            t_str = t_str[:55] + '...'
                        term_names.append(t_str)

                    ax.set_yticklabels(term_names, fontsize=7)
                    ax.set_xlabel('Normalized Enrichment Score (NES)')
                    ax.axvline(0, color='black', lw=0.5)

                    tissue_lbl = TISSUE_LABELS.get(tissue, tissue)
                    comp_lbl = COMP_LABELS.get(comp, comp)
                    ax.set_title(f'GSEA — {tissue_lbl} | {comp_lbl}')

                    fig.tight_layout()
                    plot_png = os.path.join(GSEA_DIR, f'gsea_{tissue}_{comp}.png')
                    fig.savefig(plot_png)
                    plt.close(fig)

                n_sig = (res['FDR q-val'].astype(float) < 0.25).sum()
                flush_print(f"    Done: {len(res)} sets, {n_sig} significant (FDR<0.25)")

            except Exception as e:
                flush_print(f"    ERROR: {e}")

    flush_print("\nGSEA completion finished.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Finish missing GSEA runs driven by experiment config.",
        epilog=config_arg_help(),
    )
    parser.add_argument("--experiment", default=None, help="Experiment id under experiments/<id>/transcriptomics_data/")
    parser.add_argument("--config", default=None, help="Path to transcriptomics run-entry config.yaml or batch.yaml")
    parser.add_argument("--outdir", default=None, help="Override output_dir (transcriptomics_results root)")
    args = parser.parse_args()

    cfg_info = _apply_config(args.experiment, args.config, args.outdir)
    flush_print(f"Experiment: {cfg_info.get('experiment_id')}")
    flush_print(f"Config:     {cfg_info.get('config_path')}")
    flush_print(f"Output:     {cfg_info.get('output_dir')}")

    run_gsea()
