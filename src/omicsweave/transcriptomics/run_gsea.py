#!/usr/bin/env python3
"""
GSEA (Gene Set Enrichment Analysis) — thin CLI wrapper.
========================================================
All logic lives in omicsweave/transcriptomics/gsea/.

Usage:
  python run_gsea.py --experiment <id>
  python run_gsea.py --config <path/to/config.yaml>
"""

import os
import argparse


from .experiment_config import config_arg_help
from .gsea import run_gsea


def main():
    parser = argparse.ArgumentParser(
        description="GSEA prerank analysis for RNA-Seq DE results.",
        epilog=config_arg_help(),
    )
    parser.add_argument("--experiment", default=None,
                        help="Experiment id under experiments/<id>/transcriptomics_data/")
    parser.add_argument("--config", default=None,
                        help="Path to transcriptomics run-entry config.yaml or batch.yaml")
    parser.add_argument("--outdir", default=None,
                        help="Override output_dir (transcriptomics_results root)")
    args = parser.parse_args()
    run_gsea(experiment=args.experiment, config_path=args.config, outdir=args.outdir)


if __name__ == "__main__":
    main()
