"""
Proteomics Visualization Pipeline (Config-Driven)
===================================================

Thin orchestrator that delegates to the plotting/ and enrichment/ subpackages.
Generates volcano plots, heatmaps, PCA plots, and ssGSEA pathway results
from pre-computed differential analysis results.

Usage:
    python plot_results.py --experiment <id>
    python plot_results.py --experiment <id> --skip-ssgsea
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path


from .experiment_config import config_arg_help
from .plotting import (
    apply_config,
    generate_volcano_plots,
    generate_heatmap_plots,
    generate_pca_plots,
    generate_pathway_plots,
    flush_print,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Proteomics visualization pipeline (config-driven).",
        epilog=config_arg_help(),
    )
    parser.add_argument("--experiment", default=None, help="Experiment ID")
    parser.add_argument("--config", default=None, help="Path to config.yaml")
    parser.add_argument("--outdir", default=None, help="Override output directory")
    parser.add_argument(
        "--skip-gsea", "--skip-ssgsea", dest="skip_gsea", action="store_true",
        help="Skip both ssGSEA and pre-ranked GSEA",
    )
    args = parser.parse_args()

    cfg_info = apply_config(args.experiment, args.config, args.outdir)

    flush_print(f"Experiment: {cfg_info['experiment_id']}")
    flush_print(f"Output: {cfg_info['output_dir']}")

    # Old outputs may contain raw p-values in Padj columns.  Never redraw them
    # as if they came from the corrected pipeline.
    output_root = Path(cfg_info["output_dir"])
    for tissue in cfg_info["_cfg"].tissues:
        combined = output_root / "01_dp_results" / tissue / "differential_proteins.csv"
        metadata_path = output_root / "01_dp_results" / tissue / "analysis_metadata.json"
        if combined.exists() and not metadata_path.exists():
            raise RuntimeError(
                f"[{tissue}] legacy proteomics results detected. "
                "Re-run run_proteomics_pipeline.py --force before plotting."
            )
        if metadata_path.exists():
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            if str(metadata.get("pipeline_version")) != "6.0":
                raise RuntimeError(
                    f"[{tissue}] results were produced by pipeline "
                    f"{metadata.get('pipeline_version')!r}; re-run with --force."
                )

    generate_volcano_plots()
    generate_heatmap_plots()
    generate_pca_plots()
    if not args.skip_gsea:
        generate_pathway_plots()

    flush_print("\nAll plots generated.")


if __name__ == "__main__":
    main()
