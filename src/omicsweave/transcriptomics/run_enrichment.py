#!/usr/bin/env python3
"""
GO/KEGG Pathway Enrichment Analysis — thin CLI wrapper.
========================================================
All logic lives in omicsweave/transcriptomics/enrichment/.

Usage:
  python run_enrichment.py --experiment <id>
  python run_enrichment.py --config <path/to/config.yaml>
"""

import os
import argparse


from .experiment_config import config_arg_help
from .enrichment import run_enrichment


def main():
    parser = argparse.ArgumentParser(
        description="GO/KEGG pathway enrichment analysis (Enrichr).",
        epilog=config_arg_help(),
    )
    parser.add_argument("--experiment", default=None,
                        help="Experiment id under experiments/<id>/transcriptomics_data/")
    parser.add_argument("--config", default=None,
                        help="Path to transcriptomics run-entry config.yaml or batch.yaml")
    parser.add_argument("--outdir", default=None,
                        help="Override output_dir (transcriptomics_results root)")
    args = parser.parse_args()
    run_enrichment(experiment=args.experiment, config_path=args.config, outdir=args.outdir)


if __name__ == "__main__":
    main()
