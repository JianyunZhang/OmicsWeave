#!/usr/bin/env python3
"""
Summarize GSEA results with metabolic keyword filter (config-driven)
===================================================================

Scans all GSEA result CSVs in the output directory matching
{tissue} × {comparison} from the experiment config.
No hardcoded file list.

Usage:
  python summarize_gsea.py --experiment <id>
  python summarize_gsea.py --config <path/to/config.yaml>
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

KEYWORDS = [
    "metabol", "lipid", "cholesterol", "fatty", "beta", "oxid", "mitochond",
    "respir", "electron", "tca", "citrate", "glycol", "glucose", "insulin",
    "ppar", "srebp", "bile", "peroxis", "ketone", "thermogen", "nafld",
]


def keyword_filter(df: pd.DataFrame) -> pd.DataFrame:
    term_l = df["Term"].astype(str).str.lower()
    mask = term_l.str.contains(KEYWORDS[0])
    for k in KEYWORDS[1:]:
        mask = mask | term_l.str.contains(k)
    return df[mask].copy()


def main():
    global BASE, TISSUES, COMPARISONS

    parser = argparse.ArgumentParser(
        description="Summarize GSEA results with metabolic keyword filter (config-driven).",
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
    BASE = os.path.join(os.environ["TRANSCRIPTOMICS_OUTDIR"], "05_result_plot", "5_gsea")
    TISSUES = list(cfg.tissues)

    # Auto-generate comparison file suffixes
    comp_tuples = resolve_comparisons(cfg)
    COMPARISONS = [name for _, _, name in comp_tuples]

    # Dynamically build file list from tissues × comparisons
    files = []
    for tissue in TISSUES:
        for comp in COMPARISONS:
            fn = f"gsea_{tissue}_{comp}_results.csv"
            files.append(fn)

    for fn in files:
        path = os.path.join(BASE, fn)
        print("\n" + "=" * 80)
        print(fn)
        print("=" * 80)
        if not os.path.exists(path):
            print("MISSING")
            continue

        df = pd.read_csv(path)
        df["NES"] = pd.to_numeric(df["NES"], errors="coerce")
        df["FDR q-val"] = pd.to_numeric(df["FDR q-val"], errors="coerce")
        sig = df[df["FDR q-val"] < 0.25].copy()
        print(f"Significant gene sets (FDR<0.25): {len(sig)} / {len(df)}")

        sub = keyword_filter(sig)
        if sub.empty:
            print("No keyword-matched metabolic/mitochondrial sets among significant results.")
            continue

        pos = sub.sort_values("NES", ascending=False).head(12)
        neg = sub.sort_values("NES", ascending=True).head(12)

        print("\nTop +NES (enriched toward UP direction)")
        for _, r in pos.iterrows():
            print(f"  NES={r['NES']:+.2f}  FDR={r['FDR q-val']:.3g}  {r['Term']}")

        print("\nTop -NES (enriched toward DOWN direction)")
        for _, r in neg.iterrows():
            print(f"  NES={r['NES']:+.2f}  FDR={r['FDR q-val']:.3g}  {r['Term']}")


if __name__ == "__main__":
    main()
