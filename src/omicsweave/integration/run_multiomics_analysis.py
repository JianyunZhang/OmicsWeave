#!/usr/bin/env python3
"""
Multi-omics integration pipeline CLI entry point.

Usage:
    # Via experiment ID (looks for config in experiments/<id>/multiomics_integration_setting/)
    python run_multiomics_analysis.py --experiment 20250425_cimr_wt_mice_dio

    # Via direct config path
    python run_multiomics_analysis.py --config path/to/config.yaml

    # Override output directory
    python run_multiomics_analysis.py --experiment <id> --outdir /custom/output

Pipeline: load config -> load data -> score -> converge -> plot -> report

Output structure:
    multiomics_results/
      01_do_results/<tissue>/do_<tissue>_<comparison>.csv   # per-tissue per-comparison scores
      01_do_results/<tissue>/pathway_<tissue>_<comparison>.csv
      02_result_plot/<tissue>/mfa_volcano_*.png    # per-tissue volcano plots
      02_result_plot/*.png                         # cross-tissue plots
      common/                                      # merged tables + auxiliary CSVs
      multiomics_report.md

Author: Claude Code
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


from omicsweave.integration.config import load_multiomics_config
from omicsweave.integration.loaders import (
    MFALoader,
    TranscriptomicsLoader,
    ProteomicsLoader,
    MetabolomicsLoader,
)
from omicsweave.integration.methods import (
    FluxOmicsScorer,
    PathwayConvergenceAnalyzer,
)
from omicsweave.integration.plotting import (
    apply_config,
    plot_data_availability,
    plot_differential_summary,
    plot_tissue_dashboard,
    plot_mfa_volcanoes,
    plot_pathway_condition_heatmap,
    plot_tca_cycle,
    plot_convergence_heatmap,
    plot_tissue_flux_overview,
    plot_score_distribution,
)
from omicsweave.integration.report import generate_report
from omicsweave.integration.provenance import build_input_audit
from omicsweave.integration.mappings import audit_mapping_integrity


def _save_tissue_csvs(
    df: pd.DataFrame,
    output_dir: Path,
    sub_dir: str,
    group_cols: list[str],
    file_prefix: str = "",
) -> int:
    """Save per-tissue per-comparison CSVs under output_dir/sub_dir/<tissue>/.

    Files are named ``{file_prefix}{tissue}_{comparison}.csv``, e.g.
    ``do_liver_HFD_vs_LFD.csv`` or ``pathway_heart_HSD_vs_LFD.csv``.

    Returns the number of files written.
    """
    if df.empty or len(group_cols) < 2:
        return 0
    n_files = 0
    for (tissue, comp), subset in df.groupby(group_cols[:2]):
        path = output_dir / sub_dir / tissue / f"{file_prefix}{tissue}_{comp}.csv"
        path.parent.mkdir(parents=True, exist_ok=True)
        subset.to_csv(path, index=False, encoding="utf-8-sig")
        n_files += 1
    return n_files


def run_pipeline(
    experiment_id: str | None = None,
    config_path: str | None = None,
    outdir: str | None = None,
) -> None:
    """Execute the full multi-omics integration pipeline.

    Parameters
    ----------
    experiment_id : str, optional
        Experiment ID for config discovery.
    config_path : str, optional
        Direct path to config.yaml.
    outdir : str, optional
        Override output directory.
    """
    # ---- 1. Load configuration ----
    print("=" * 70)
    print("Multi-Omics Integration Pipeline")
    print("=" * 70)

    cfg = load_multiomics_config(
        experiment_id=experiment_id,
        config_path=config_path,
        outdir_override=outdir,
    )
    print(f"Experiment: {cfg.experiment_id}")
    print(f"Output:     {cfg.output_dir}")
    print(f"Tissues:    {', '.join(cfg.tissues)}")
    print(f"Comparisons: {', '.join(cfg.comparisons)}")
    print()

    # Apply plotting config (sets COMMON_DIR, PLOTS_DIR, etc.)
    apply_config(
        experiment=experiment_id,
        config_path=config_path,
        outdir=str(cfg.output_dir) if outdir else None,
    )

    # Shortcut directories
    common_dir = cfg.output_dir / "common"
    plots_dir = cfg.output_dir / "02_result_plot"
    do_dir = cfg.output_dir / "01_do_results"
    for d in [common_dir, plots_dir]:
        d.mkdir(parents=True, exist_ok=True)

    mapping_audit = audit_mapping_integrity()
    mapping_audit.raise_for_errors()
    mapping_audit_rows = [
        {"record_type": "metric", "key": key, "value": value, "message": ""}
        for key, value in sorted(mapping_audit.metrics.items())
    ]
    mapping_audit_rows.extend(
        {"record_type": "warning", "key": "", "value": "", "message": message}
        for message in mapping_audit.warnings
    )
    mapping_audit_df = pd.DataFrame(mapping_audit_rows)
    mapping_audit_df.to_csv(
        common_dir / "mapping_integrity_audit.csv",
        index=False,
        encoding="utf-8-sig",
    )
    print(
        "  Saved: common/mapping_integrity_audit.csv "
        f"({len(mapping_audit.warnings)} warning(s), 0 error(s))"
    )

    # ---- 2. Load data ----
    print("--- Loading Data ---")

    # MFA
    mfa_loader = MFALoader(cfg)
    mfa_data = mfa_loader.load()

    # Transcriptomics
    tx_loader = TranscriptomicsLoader(cfg)
    tx_data = tx_loader.load()

    # Proteomics
    prot_loader = ProteomicsLoader(cfg)
    prot_data = prot_loader.load()

    # Metabolomics
    metab_loader = MetabolomicsLoader(cfg)
    metab_data = metab_loader.load()

    input_audit_df = build_input_audit(
        cfg,
        {
            "mfa": mfa_data.get("diff", {}),
            "transcriptomics": tx_data.get("de", {}),
            "proteomics": prot_data.get("diff", {}),
            "metabolomics": metab_data.get("diff", {}),
        },
    )
    input_audit_df.to_csv(
        common_dir / "input_audit.csv", index=False, encoding="utf-8-sig"
    )
    print("  Saved: common/input_audit.csv")

    print()

    # ---- 3. Scoring ----
    print("--- Scoring Reactions ---")
    scorer = FluxOmicsScorer(cfg)

    scored_df = scorer.score_all(
        mfa_diff=mfa_data["diff"],
        tx_data=tx_data["de"],
        prot_data=prot_data["diff"],
        metab_data=metab_data["diff"],
    )

    if not scored_df.empty:
        # Merged table → common/
        scored_df.to_csv(common_dir / "flux_omics_scores.csv", index=False, encoding="utf-8-sig")
        print(f"  Saved: common/flux_omics_scores.csv ({len(scored_df)} reactions)")

        evidence_df = scorer.get_evidence_table()
        if not evidence_df.empty:
            evidence_df.to_csv(
                common_dir / "evidence_table.csv", index=False, encoding="utf-8-sig"
            )
            print(f"  Saved: common/evidence_table.csv ({len(evidence_df)} layer records)")

        # Per-tissue / per-comparison splits → 01_do_results/
        n_split = _save_tissue_csvs(
            scored_df, cfg.output_dir, "01_do_results",
            ["tissue", "comparison"], file_prefix="do_",
        )
        if n_split:
            print(f"  Saved: 01_do_results/ ({n_split} tissue×comparison score files)")

        # Pathway aggregation
        pathway_df = scorer.aggregate_by_pathway(scored_df)
        if not pathway_df.empty:
            pathway_df.to_csv(common_dir / "pathway_scores.csv", index=False, encoding="utf-8-sig")
            print(f"  Saved: common/pathway_scores.csv ({len(pathway_df)} entries)")

            n_path = _save_tissue_csvs(
                pathway_df, cfg.output_dir, "01_do_results",
                ["tissue", "comparison"], file_prefix="pathway_",
            )
            if n_path:
                print(f"  Saved: 01_do_results/ ({n_path} tissue×comparison pathway files)")
    else:
        pathway_df = pd.DataFrame()
        print("  No reactions scored (no MFA differential data)")

    # ---- 4. Pathway Convergence ----
    print("\n--- Pathway Convergence ---")
    convergence_analyzer = PathwayConvergenceAnalyzer(
        output_dir=cfg.output_dir,
    )
    convergence_df = convergence_analyzer.analyze(
        tissues=cfg.tissues,
        metab_pathway=metab_data["pathway"],
        prot_pathway=prot_data["pathway"],
        tx_gsea=tx_data["gsea"],
    )

    # ---- 5. Plotting ----
    print("\n--- Generating Plots ---")

    # Data availability (significant feature counts)
    plot_data_availability(
        tissues=cfg.tissues,
        mfa_tissues=mfa_data["tissues"],
        metab_tissues=set(metab_data["diff"].keys()),
        prot_tissues=set(prot_data["diff"].keys()),
        tx_tissues=set(tx_data["gsea"].keys()),
        output_dir=plots_dir,
        mfa_diff=mfa_data["diff"],
        metab_diff=metab_data["diff"],
        prot_diff=prot_data["diff"],
        tx_de=tx_data["de"],
        pval_cutoff=scorer.pval_cutoff,
        log2fc_cutoff=scorer.log2fc_cutoff,
    )

    # Differential summary (per-comparison, cross-tissue)
    plot_differential_summary(
        mfa_diff=mfa_data["diff"],
        metab_diff=metab_data["diff"],
        prot_diff=prot_data["diff"],
        tx_de=tx_data["de"],
        tissues=cfg.tissues,
        comparisons=cfg.comparisons,
        output_dir=plots_dir,
    )

    # Tissue dashboards (per-comparison, cross-tissue)
    for comp in cfg.comparisons:
        plot_tissue_dashboard(
            mfa_diff=mfa_data["diff"],
            metab_diff=metab_data["diff"],
            prot_diff=prot_data["diff"],
            tx_gsea=tx_data["gsea"],
            tissues=cfg.tissues,
            comparison=comp,
            comparison_map=cfg.comparison_map,
            output_dir=plots_dir,
        )

    # MFA volcano plots (per-comparison + per-tissue)
    plot_mfa_volcanoes(
        mfa_diff=mfa_data["diff"],
        output_dir=plots_dir,
        tissues=cfg.tissues,
    )

    # MFA heatmaps (cross-tissue)
    plot_pathway_condition_heatmap(
        mfa_pathway=mfa_data["pathway"],
        output_dir=plots_dir,
    )
    plot_tca_cycle(
        mfa_tca=mfa_data["tca"],
        output_dir=plots_dir,
    )
    plot_tissue_flux_overview(
        mfa_overview=mfa_data["overview"],
        output_dir=plots_dir,
    )

    # Convergence heatmap
    if not convergence_df.empty:
        plot_convergence_heatmap(
            conv_df=convergence_df,
            output_dir=plots_dir,
        )

    # Scoring plots (global)
    if not scored_df.empty:
        plot_score_distribution(
            scored_df=scored_df,
            output_dir=plots_dir,
        )

    # ---- 6. Report ----
    print("\n--- Generating Report ---")
    generate_report(
        output_dir=cfg.output_dir,
        experiment_id=cfg.experiment_id,
        scored_df=scored_df if not scored_df.empty else None,
        pathway_df=pathway_df if not pathway_df.empty else None,
        convergence_df=convergence_df if not convergence_df.empty else None,
        pathway_fusion_df=(
            convergence_analyzer.fusion_df
            if not convergence_analyzer.fusion_df.empty else None
        ),
        mfa_diff=mfa_data["diff"],
        tissues=cfg.tissues,
        comparisons=cfg.comparisons,
        estimand=cfg.estimand,
        input_audit_df=input_audit_df,
        mapping_audit_df=mapping_audit_df,
    )

    print("\n" + "=" * 70)
    print("Multi-omics integration complete!")
    print(f"Results: {cfg.output_dir}")
    print("=" * 70)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Multi-Omics Integration Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python run_multiomics_analysis.py --experiment 20250425_cimr_wt_mice_dio\n"
            "  python run_multiomics_analysis.py --config experiments/<id>/multiomics_integration_setting/config.yaml\n"
        ),
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--experiment", type=str,
        help="Experiment ID (looks for config in experiments/<id>/multiomics_integration_setting/)",
    )
    group.add_argument(
        "--config", type=str,
        help="Direct path to config.yaml",
    )
    parser.add_argument(
        "--outdir", type=str, default=None,
        help="Override output directory",
    )

    args = parser.parse_args()
    run_pipeline(
        experiment_id=args.experiment,
        config_path=args.config,
        outdir=args.outdir,
    )


if __name__ == "__main__":
    main()
