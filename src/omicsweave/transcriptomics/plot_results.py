#!/usr/bin/env python3
"""
RNA-Seq Downstream Visualization (Step 5–7) — thin wrapper.
===========================================================
Step 5 (volcano/heatmap/PCA): logic in plotting/ subpackage
Step 6 (Enrichment):         logic in enrichment/ subpackage
Step 7 (GSEA):               logic in gsea/ subpackage

Usage:
  python plot_results.py --experiment <id>
  python plot_results.py --experiment <id> --skip-enrichment --skip-gsea
"""

import os
import argparse


from .experiment_config import config_arg_help
from .plotting._config import apply_config, flush_print
from .plotting import _config as plot_cfg
from .plotting import generate_volcano_plots, generate_heatmap_plots, generate_pca_plots
from .enrichment import run_enrichment
from .gsea import run_gsea


def main():
    parser = argparse.ArgumentParser(
        description="RNA-Seq downstream plots (volcano/heatmap/PCA). "
                    "Enrichment and GSEA are delegated to separate subpackages.",
        epilog=config_arg_help(),
    )
    parser.add_argument("--experiment", default=None,
                        help="Experiment id under experiments/<id>/transcriptomics_data/")
    parser.add_argument("--config", default=None,
                        help="Path to transcriptomics run-entry config.yaml or batch.yaml")
    parser.add_argument("--outdir", default=None,
                        help="Override output_dir (transcriptomics_results root)")
    parser.add_argument("--skip-enrichment", action="store_true",
                        help="Skip enrichment (run separately via run_enrichment.py)")
    parser.add_argument("--skip-gsea", action="store_true",
                        help="Skip GSEA (run separately via run_gsea.py)")
    args = parser.parse_args()

    apply_config(experiment=args.experiment, config_path=args.config, outdir=args.outdir)

    flush_print("=" * 60)
    flush_print("  RNA-Seq Downstream Visualization (Step 5)")
    flush_print("=" * 60)

    os.makedirs(plot_cfg.PLOT_DIR, exist_ok=True)

    generate_volcano_plots()
    generate_heatmap_plots()
    generate_pca_plots()

    if not args.skip_enrichment:
        run_enrichment(experiment=args.experiment, config_path=args.config, outdir=args.outdir)
    else:
        flush_print("\n  Enrichment: skipped (--skip-enrichment)")

    if not args.skip_gsea:
        run_gsea(experiment=args.experiment, config_path=args.config, outdir=args.outdir)
    else:
        flush_print("  GSEA: skipped (--skip-gsea)")

    flush_print("\n" + "=" * 60)
    flush_print("ALL PLOTS GENERATED")
    flush_print("=" * 60)
    flush_print(f"Output: {plot_cfg.PLOT_DIR}")


if __name__ == "__main__":
    main()
