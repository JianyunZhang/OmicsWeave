#!/usr/bin/env python3
"""
Print marker genes from DE result CSVs (config-driven)
=====================================================

For each tissue × comparison pair in the experiment config, prints
a selection of marker genes from the DE results.
No hardcoded tissue/comparison names.

Usage:
  python check_markers.py --experiment <id>
  python check_markers.py --config <path/to/config.yaml>
"""

import os
import argparse
import pandas as pd

from .experiment_config import (
    load_transcriptomics_config, config_arg_help, resolve_comparisons,
)

# Placeholders — set by _apply_config()
BASE: str | None = None
TISSUES: list[str] = []
COMPARISONS: list[str] = []

# Marker gene sets by tissue category (generic — applies to any experiment)
TISSUE_MARKERS: dict[str, list[str]] = {
    "liver": [
        "Srebf1", "Srebf2", "Fasn", "Acaca", "Scd1", "Elovl6", "Hmgcr", "Sqle",
        "Cpt1a", "Cpt2", "Ppara", "Pparg", "Ppargc1a", "Hmgcs2", "Pck1", "G6pc",
        "Tgfb1", "Col1a1", "Col3a1", "Acta2",
        "Il1b", "Tnf", "Ccl2", "Nos2",
    ],
    "heart": [
        "Ppargc1a", "Cpt1b", "Cpt2", "Acadm", "Acadl", "Ucp3",
        "Ndufs1", "Cox4i1", "Atp5f1a",
        "Tnni3", "Tnnt2", "Myh6", "Nppa", "Nppb",
        "Il1b", "Tnf", "Cxcl9", "Stat1", "Isg15",
    ],
    "brain": [
        "Aif1", "Trem2", "Apoe", "Gfap", "Stat1", "Isg15", "Il1b", "Tnf",
    ],
    "kidney": [
        "Havcr1", "Lcn2", "Slc22a6", "Slc34a1", "Ppara", "Ppargc1a",
        "Cpt1a", "Cpt2", "Il1b", "Tnf", "Ccl2", "Isg15",
    ],
    "gast_muscle": [
        "Ppargc1a", "Pdk4", "Cpt1b", "Cpt2", "Ucp3", "Slc2a4",
        "Myog", "Trim63", "Il1b", "Tnf", "Isg15",
    ],
    "sq_adipose": [
        "Pparg", "Adipoq", "Lep", "Plin1", "Cpt1a", "Cpt2",
        "Ppargc1a", "Ucp1", "Tnf", "Il1b", "Ccl2",
    ],
    "go_adipose": [
        "Pparg", "Adipoq", "Lep", "Plin1", "Cpt1a", "Cpt2",
        "Ppargc1a", "Ucp1", "Tnf", "Il1b", "Ccl2",
    ],
}


def show(file_name: str, genes: list[str]) -> None:
    path = os.path.join(BASE, file_name)
    if not os.path.exists(path):
        print(f"\n  {file_name}: MISSING")
        return
    df = pd.read_csv(path).set_index("gene")
    print(f"\n=== {file_name} ===")
    for g in genes:
        if g in df.index:
            r = df.loc[g]
            print(
                f"{g:12s}  log2FC={r['log2FoldChange']:+.3f}  "
                f"padj={r['padj']:.3g}  baseMean={r['baseMean']:.1f}"
            )
        else:
            print(f"{g:12s}  NOT_FOUND")


def main():
    global BASE, TISSUES, COMPARISONS

    parser = argparse.ArgumentParser(
        description="Print marker genes from DE result CSVs (config-driven).",
        epilog=config_arg_help(),
    )
    parser.add_argument("--experiment", default=None, help="Experiment id under experiments/<id>/transcriptomics*/")
    parser.add_argument("--config", default=None, help="Path to transcriptomics run-entry config.yaml or batch.yaml")
    parser.add_argument("--outdir", default=None, help="Override output_dir (transcriptomics_results root)")
    args = parser.parse_args()

    cfg = load_transcriptomics_config(
        experiment_id=args.experiment,
        config_path=args.config,
        outdir_override=args.outdir,
    )
    os.environ["TRANSCRIPTOMICS_OUTDIR"] = os.path.abspath(str(cfg.output_dir))
    BASE = os.path.join(os.environ["TRANSCRIPTOMICS_OUTDIR"], "04_de_results")
    TISSUES = list(cfg.tissues)

    # Auto-generate comparison file suffixes
    comp_tuples = resolve_comparisons(cfg)
    COMPARISONS = [name for _, _, name in comp_tuples]

    # Iterate all tissue × comparison pairs
    for tissue in TISSUES:
        markers = TISSUE_MARKERS.get(tissue, [])
        if not markers:
            continue
        for comp in COMPARISONS:
            fn = os.path.join(tissue, f"de_{comp}.csv")
            show(fn, markers)


if __name__ == "__main__":
    main()
