#!/usr/bin/env python3
"""
Re-run enrichment (Step 4) and GSEA (Step 5) (config-driven)
=============================================================

Cleans old results and regenerates enrichment/GSEA plots using
functions from plot_results.py (which reads experiment config).

Usage:
  python rerun_enrichment.py --experiment <id>
  python rerun_enrichment.py --config <path/to/config.yaml>
"""

import os
import argparse
import shutil


from .experiment_config import load_transcriptomics_config, config_arg_help
from .plot_results import (
    generate_enrichment_plots, generate_gsea_plots,
    _apply_config, flush_print,
)


def main():
    parser = argparse.ArgumentParser(
        description="Re-run enrichment and GSEA steps driven by experiment config.",
        epilog=config_arg_help(),
    )
    parser.add_argument("--experiment", default=None, help="Experiment id under experiments/<id>/transcriptomics*/")
    parser.add_argument("--config", default=None, help="Path to transcriptomics run-entry config.yaml or batch.yaml")
    parser.add_argument("--outdir", default=None, help="Override output_dir (transcriptomics_results root)")
    args = parser.parse_args()

    # Apply config to plot_results module globals
    from .plot_results import PLOT_DIR as _plot_dir_before
    cfg_info = _apply_config(args.experiment, args.config, args.outdir)

    # Import PLOT_DIR again after config is applied
    from . import plot_results
    plot_dir = plot_results.PLOT_DIR

    flush_print(f"Experiment: {cfg_info.get('experiment_id')}")
    flush_print(f"Config:     {cfg_info.get('config_path')}")
    flush_print(f"Output:     {plot_dir}")

    # Clean old results
    for sub in ['4_enrichment', '5_gsea']:
        d = os.path.join(plot_dir, sub)
        if os.path.exists(d):
            shutil.rmtree(d)
            flush_print(f"Cleaned: {d}")

    generate_enrichment_plots()
    generate_gsea_plots()

    flush_print("\nDone! Steps 4 & 5 re-generated with mouse-specific databases.")


if __name__ == "__main__":
    main()
