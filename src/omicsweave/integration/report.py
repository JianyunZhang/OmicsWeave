#!/usr/bin/env python3
"""
Markdown report generator for multi-omics integration analysis.

Adapted from flux_multiomics_integration.py lines 578-733 (report)
and multiomics_integration_analysis.py lines 912-965 (report).

Generates a self-contained Markdown report with key statistics, confidence
distribution, top hits, pathway summaries, and output file listings.

Author: Claude Code
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd


def generate_report(
    output_dir: Path,
    experiment_id: str,
    scored_df: pd.DataFrame | None = None,
    pathway_df: pd.DataFrame | None = None,
    convergence_df: pd.DataFrame | None = None,
    pathway_fusion_df: pd.DataFrame | None = None,
    mfa_diff: dict[str, pd.DataFrame] | None = None,
    tissues: list[str] | None = None,
    comparisons: list[str] | None = None,
    estimand: dict[str, Any] | None = None,
    input_audit_df: pd.DataFrame | None = None,
    mapping_audit_df: pd.DataFrame | None = None,
) -> Path:
    """Generate a Markdown summary report.

    Parameters
    ----------
    output_dir : Path
        Output directory containing common/, 01_do_results/, 02_result_plot/.
    experiment_id : str
        Experiment identifier.
    scored_df : DataFrame, optional
        Scored reactions from FluxOmicsScorer.
    pathway_df : DataFrame, optional
        Pathway-aggregated scores.
    convergence_df : DataFrame, optional
        Cross-omics pathway convergence results.
    pathway_fusion_df : DataFrame, optional
        Fisher/RRA pathway evidence fusion results.
    mfa_diff : dict, optional
        MFA differential flux tables.
    tissues : list[str], optional
        Tissues analyzed.
    comparisons : list[str], optional
        Comparisons analyzed.

    Returns
    -------
    Path
        Path to the generated report file.
    """
    lines: list[str] = []

    lines.append("# Multi-Omics Integration Report\n\n")
    lines.append(f"**Experiment:** {experiment_id}  \n")
    lines.append(f"**Generated:** {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}  \n")
    if tissues:
        lines.append(f"**Tissues:** {', '.join(sorted(tissues))}  \n")
    if comparisons:
        lines.append(f"**Comparisons:** {', '.join(sorted(comparisons))}  \n")
    if estimand:
        lines.append(f"**Estimand:** {estimand.get('id', 'not specified')}  \n")
    lines.append("\n---\n\n")

    common_dir = output_dir / "common"
    plots_dir = output_dir / "02_result_plot"
    do_dir = output_dir / "01_do_results"

    # ---- Data availability ----
    lines.append("## 1. Data Overview\n\n")
    overview_file = common_dir / "multi_omics_overview.csv"
    if overview_file.exists():
        lines.append("![Data Availability](02_result_plot/multi_omics_overview.png)\n\n")
    if input_audit_df is not None and not input_audit_df.empty:
        lines.append("### Input and Batch Provenance\n\n")
        lines.append(
            "Configured batch counts describe requested inputs; actual batch counts are "
            "derived from loaded result provenance and may be lower when a feature only has "
            "single-batch coverage.\n\n"
        )
        lines.append(
            "| Omics | Config mode | Config batches | Actual source | Actual batches "
            "| Loaded rows | Status |\n"
        )
        lines.append("|---|---|---:|---|---:|---:|---|\n")
        for _, row in input_audit_df.iterrows():
            actual_batches = row.get("actual_n_batches", "")
            actual_batches = "" if pd.isna(actual_batches) else f"{actual_batches:g}"
            lines.append(
                f"| {row.get('omics', '')} | {row.get('configured_mode', '')} "
                f"| {row.get('configured_n_batches', 0)} "
                f"| {row.get('actual_source_modes', '')} | {actual_batches} "
                f"| {row.get('n_loaded_rows', 0)} | {row.get('status', '')} |\n"
            )
        lines.append("\nFull paths and batch IDs are retained in `common/input_audit.csv`.\n\n")
    if mapping_audit_df is not None and not mapping_audit_df.empty:
        warning_rows = mapping_audit_df.loc[
            mapping_audit_df["record_type"].eq("warning")
        ]
        metric_rows = mapping_audit_df.loc[
            mapping_audit_df["record_type"].eq("metric")
        ]
        metric_lookup = dict(zip(metric_rows["key"], metric_rows["value"]))
        lines.append("### Feature Mapping Integrity\n\n")
        lines.append(
            "Feature matching uses standard database IDs first, canonical mouse gene/protein "
            "symbols second, and normalized exact metabolite aliases only when no standard ID "
            "is available. Ambiguous role matches are excluded from scoring.\n\n"
        )
        lines.append(
            f"- Reactions mapped in all three omics layers: "
            f"{metric_lookup.get('all_three_reaction_count', 0)}\n"
        )
        lines.append(
            f"- Metabolite aliases without a standard ID: "
            f"{metric_lookup.get('metabolite_aliases_without_standard_id', 0)}\n"
        )
        lines.append(f"- Non-blocking mapping warnings: {len(warning_rows)}\n\n")
        lines.append(
            "Full metrics and warnings are retained in "
            "`common/mapping_integrity_audit.csv`.\n\n"
        )

    # ---- MFA Summary ----
    if mfa_diff:
        lines.append("## 2. MFA Differential Flux Summary\n\n")
        lines.append("| Comparison | Total fluxes | Significant (p<0.05) | Gate passed (p<0.05 & |FC|>0.5) |\n")
        lines.append("|---|---|---|---|\n")
        for comp, df in sorted(mfa_diff.items()):
            n_total = len(df)
            pcol = "flux_pvalue" if "flux_pvalue" in df.columns else "pvalue_ttest"
            fcol = "flux_log2fc" if "flux_log2fc" in df.columns else "log2FC"
            n_sig_p = int((df[pcol] < 0.05).sum()) if pcol in df.columns else 0
            n_sig_fc = 0
            if pcol in df.columns and fcol in df.columns:
                n_sig_fc = int(((df[pcol] < 0.05) & (df[fcol].abs() > 0.5)).sum())
            lines.append(f"| {comp} | {n_total} | {n_sig_p} | {n_sig_fc} |\n")
        lines.append("\n")

    # ---- Scoring Results ----
    if scored_df is not None and not scored_df.empty:
        lines.append("## 3. Multi-Omics Evidence Profiles\n\n")
        total = len(scored_df)
        lines.append(f"**Total reactions assessed:** {total}\n\n")

        if "evidence_grade" in scored_df.columns:
            grade_counts = scored_df["evidence_grade"].value_counts()
            lines.append("### Evidence Grade Distribution\n\n")
            lines.append("| Evidence grade | Count | Percentage |\n|---|---|---|\n")
            grade_order = [
                "Replicated", "Convergent", "Single-layer", "Exploratory",
                "Conflicting", "No-signal", "Not-assessable",
            ]
            for grade in grade_order:
                n_grade = int(grade_counts.get(grade, 0))
                if n_grade or grade in {"Convergent", "Exploratory", "Conflicting"}:
                    lines.append(f"| {grade} | {n_grade} | {n_grade / total * 100:.1f}% |\n")
            lines.append("\n")

        # 旧评分仅保留用于历史结果复现和现有图表兼容。
        lines.append("### Legacy 0-4 Score Distribution\n\n")
        lines.append("![Score Distribution](02_result_plot/score_distribution.png)\n\n")
        lines.append("![Confidence Pie](02_result_plot/confidence_pie.png)\n\n")

        conf_counts = scored_df["confidence"].value_counts()
        lines.append("| Confidence | Count | Percentage |\n|---|---|---|\n")
        for conf in ["High", "Medium", "Low", "Insufficient"]:
            n = conf_counts.get(conf, 0)
            pct = n / total * 100 if total > 0 else 0
            lines.append(f"| {conf} | {n} | {pct:.1f}% |\n")
        lines.append("\n")

        # Top hits
        lines.append("### Top Multi-Omics Evidence Profiles\n\n")
        lines.append(
            "| Flux ID | Tissue | Comparison | Support | Conflict | Coverage "
            "| Design | Evidence grade | Legacy score |\n"
        )
        lines.append(
            "|---|---|---|---:|---:|---:|---|---|---:|\n"
        )
        for _, row in scored_df.head(20).iterrows():
            lines.append(
                f"| {row.get('flux_id', '')} "
                f"| {row.get('tissue', '')} "
                f"| {row.get('comparison', '')} "
                f"| {row.get('support_score', 0):.3f} "
                f"| {row.get('conflict_score', 0):.3f} "
                f"| {row.get('coverage_score', 0):.3f} "
                f"| {row.get('design_grade', '')} "
                f"| {row.get('evidence_grade', '')} "
                f"| {row.get('legacy_total_score', row.get('total_score', 0)):.3f} |\n"
            )
        lines.append("\n")

    # ---- Pathway Aggregation ----
    if pathway_df is not None and not pathway_df.empty:
        lines.append("## 4. Pathway-Level Analysis\n\n")
        lines.append(
            "| Pathway | Tissue | Comparison | N Reactions | Mean support "
            "| Mean conflict | Mean coverage | Convergent | Conflicting |\n"
        )
        lines.append("|---|---|---|---:|---:|---:|---:|---:|---:|\n")
        for _, row in pathway_df.head(20).iterrows():
            lines.append(
                f"| {row.get('pathway', '')} "
                f"| {row.get('tissue', '')} "
                f"| {row.get('comparison', '')} "
                f"| {row.get('n_reactions', 0)} "
                f"| {row.get('mean_support_score', 0):.3f} "
                f"| {row.get('mean_conflict_score', 0):.3f} "
                f"| {row.get('mean_coverage_score', 0):.3f} "
                f"| {row.get('n_convergent', 0)} "
                f"| {row.get('n_conflicting', 0)} |\n"
            )
        lines.append("\n")

    # ---- Pathway Convergence ----
    if convergence_df is not None and not convergence_df.empty:
        lines.append("## 5. Cross-Omics Pathway Convergence\n\n")
        lines.append("![Convergence Heatmap](02_result_plot/pathway_convergence_heatmap.png)\n\n")
        n_shared = len(convergence_df)
        n_tissues = convergence_df["tissue"].nunique()
        n_comparisons = (
            convergence_df["comparison"].nunique()
            if "comparison" in convergence_df.columns
            else 0
        )
        lines.append(
            f"Found **{n_shared}** shared pathway terms across "
            f"**{n_tissues}** tissues"
            + (f" and **{n_comparisons}** comparisons" if n_comparisons else "")
            + ". Intersections are never pooled across comparisons.\n\n"
        )

    if pathway_fusion_df is not None and not pathway_fusion_df.empty:
        n_fused = len(pathway_fusion_df)
        n_fdr = int((pathway_fusion_df["fused_fdr"] < 0.05).sum())
        lines.append(
            f"Statistical/rank fusion evaluated **{n_fused}** pathways; "
            f"**{n_fdr}** passed fused FDR < 0.05. Omics contributions are "
            "retained in `tables/pathway_fusion.csv`.\n\n"
        )

    # ---- Methodology ----
    lines.append("## 6. Methodology\n\n")
    lines.append("### Multi-dimensional Evidence Strategy\n\n")
    lines.append(
        "The primary result is no longer a fixed 0-4 total. Each omics layer contributes "
        "continuous evidence based on direction, sign certainty, effect strength, mapping "
        "quality and data reliability. Available layers are normalized independently from "
        "coverage, so missing data are not treated as biological absence.\n\n"
    )
    lines.append(
        "The reaction table reports `support_score`, `conflict_score`, `coverage_score`, "
        "bootstrap stability (when biological-replicate values are available), `design_grade` "
        "and `evidence_grade`. EvidenceTable explicitly separates `no_data`, `no_mapping`, "
        "`no_match`, `no_signal` and `opposite`. MFA MCMC fitting points are marked as "
        "pseudo-replicate evidence and cannot receive a confirmatory design grade.\n\n"
    )
    lines.append(
        "The former score matrix remains in `legacy_*` columns solely for backward "
        "compatibility. Its MFA gate and High/Medium/Low labels are not used as the new "
        "cross-experiment confidence claim.\n\n"
    )
    lines.append(
        "Pathway integration additionally writes exact-name overlaps for compatibility and "
        "a Fisher/RRA evidence-fusion table with per-omics source retention and BH FDR.\n\n"
    )

    # ---- Output Files ----
    lines.append("## 7. Output Files\n\n")

    # common/
    if common_dir.exists():
        csv_files = sorted(common_dir.glob("*.csv"))
        if csv_files:
            lines.append("### common/ (merged tables)\n\n")
            for f in csv_files:
                lines.append(f"- `common/{f.name}`\n")
            lines.append("\n")

    # 01_do_results/
    if do_dir.exists():
        tissue_dirs = sorted(d for d in do_dir.iterdir() if d.is_dir())
        if tissue_dirs:
            lines.append("### 01_do_results/ (per-tissue per-comparison)\n\n")
            for td in tissue_dirs:
                csvs = sorted(td.glob("*.csv"))
                if csvs:
                    lines.append(f"**{td.name}/**\n\n")
                    for f in csvs:
                        lines.append(f"- `01_do_results/{td.name}/{f.name}`\n")
                    lines.append("\n")

    # 02_result_plot/
    if plots_dir.exists():
        # Top-level plots (cross-tissue)
        top_pngs = sorted(p for p in plots_dir.glob("*.png"))
        if top_pngs:
            lines.append("### 02_result_plot/ (cross-tissue plots)\n\n")
            for f in top_pngs:
                lines.append(f"- `02_result_plot/{f.name}`\n")
            lines.append("\n")

        # Per-tissue plots
        tissue_plot_dirs = sorted(d for d in plots_dir.iterdir() if d.is_dir())
        if tissue_plot_dirs:
            lines.append("### 02_result_plot/ (per-tissue plots)\n\n")
            for td in tissue_plot_dirs:
                pngs = sorted(td.glob("*.png"))
                if pngs:
                    lines.append(f"**{td.name}/**\n\n")
                    for f in pngs:
                        lines.append(f"- `02_result_plot/{td.name}/{f.name}`\n")
                    lines.append("\n")

    lines.append("---\n\n")
    lines.append("*This report was automatically generated by the Multi-Omics Integration Pipeline.*\n")

    report_path = output_dir / "multiomics_report.md"
    report_path.write_text("".join(lines), encoding="utf-8")
    print(f"[OK] Report saved: {report_path}")
    return report_path
