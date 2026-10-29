"""Config-driven DIA proteomics analysis pipeline (v6.0).

Statistical inference uses log2-normalized, unimputed intensities. Deterministic
down-shifted imputation is retained only for PCA, heatmaps and sample-level
ssGSEA. Comparison-level pathway inference uses pre-ranked GSEA on moderated
t-statistics from all testable proteins.

``--search`` optionally runs DIA-NN 1.9 on Thermo ``.raw`` files first. Without
that flag the pipeline still reads an existing DIA-NN report.
"""

from __future__ import annotations

import argparse
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

from .differential_analysis import (  # noqa: E402
    COMBINED_CSV,
    estimate_limma_priors,
    run_differential_analysis,
)
from .enrichment.gsea import run_prerank_gsea, run_ssgsea  # noqa: E402
from .diann_search import run_tissue_search  # noqa: E402
from .experiment_config import (  # noqa: E402
    ProteomicsConfig,
    build_sample_metadata,
    config_arg_help,
    generate_comparisons,
    get_group_colors,
    get_sample_to_group_map,
    load_proteomics_config,
)
from .plotting.heatmap import plot_single as plot_heatmap  # noqa: E402
from .plotting.pca import plot_single as plot_pca  # noqa: E402
from .plotting.volcano import plot_single_from_combined as plot_volcano  # noqa: E402
from .proteomics_qc import generate_qc_outputs  # noqa: E402


PIPELINE_VERSION = "6.0"
_COMBINED_CSV = COMBINED_CSV
# Backward-compatible public name used by focused tests and downstream imports.
_estimate_limma_priors = estimate_limma_priors


def flush_print(message: str) -> None:
    print(message, flush=True)


def _extract_sample_name(run_name: str) -> str:
    """Extract a stable tissue+animal token from a DIA-NN Run value."""

    name = str(run_name)
    # e.g. ...ileum_20260124_A35 -> ileumA35
    match = re.search(
        r"_([A-Za-z]+(?:_[A-Za-z]+)*)_\d{6,8}_([A-Za-z]*\d+)$", name,
    )
    if match:
        return match.group(1).replace("_", "") + match.group(2)
    # e.g. ...brown_adipose_10 / ...liver_A35
    match = re.search(r"_((?:[A-Za-z]+_)*[A-Za-z]+)_([A-Za-z]*\d+)$", name)
    if match:
        return match.group(1) + match.group(2)
    # e.g. ...Heart52 / ...LiverA35
    match = re.search(r"_([A-Za-z_]+[A-Za-z]\d+)$", name)
    if match:
        return match.group(1).replace("_", "")
    matches = list(re.finditer(r"([A-Za-z]+\d+)", name))
    return matches[-1].group(0) if matches else name


def load_and_clean_data(
    file_path: str,
    *,
    protein_aggregation: str = "log2_mean",
) -> tuple[pd.DataFrame | None, dict[str, str]]:
    """Load a DIA-NN/Spectronaut report and build a protein × sample matrix.

    Multiple protein groups that resolve to one feature are aggregated in log2
    space by default (geometric mean on the original scale).
    """

    flush_print(f"  Loading data from {os.path.basename(file_path)}...")
    usecols = [
        "Run", "Protein.Group", "Genes", "PG.MaxLFQ",
        "PG.Q.Value", "Global.PG.Q.Value",
    ]
    try:
        report = pd.read_csv(file_path, sep="\t", usecols=usecols)
    except ValueError:
        flush_print("  Warning: standard DIA-NN columns not found; trying fallback.")
        report = pd.read_csv(file_path, sep="\t")
        if "PG.Quantity" in report.columns and "PG.MaxLFQ" not in report.columns:
            report["PG.MaxLFQ"] = report["PG.Quantity"]
    except FileNotFoundError:
        flush_print(f"  Error: file not found {file_path}")
        return None, {}

    required = {"Run", "Protein.Group", "PG.MaxLFQ"}
    missing = required - set(report.columns)
    if missing:
        raise ValueError(f"Proteomics report missing required columns: {sorted(missing)}")

    for qvalue_column in ("PG.Q.Value", "Global.PG.Q.Value"):
        if qvalue_column in report.columns:
            qvalues = pd.to_numeric(report[qvalue_column], errors="coerce")
            report = report[qvalues < 0.01].copy()

    report["Sample"] = report["Run"].map(_extract_sample_name)
    if "Genes" in report.columns:
        genes = report["Genes"].replace(r"^\s*$", np.nan, regex=True)
        report["Feature"] = genes.fillna(report["Protein.Group"])
    else:
        report["Feature"] = report["Protein.Group"]
    report["Feature"] = report["Feature"].astype(str)

    feature_to_uniprot: dict[str, str] = {}
    feature_rows = report.dropna(subset=["Protein.Group"]).drop_duplicates("Feature")
    for feature, protein_group in feature_rows[["Feature", "Protein.Group"]].to_numpy():
        feature_to_uniprot[str(feature)] = str(protein_group).split(";")[0]

    quantities = pd.to_numeric(report["PG.MaxLFQ"], errors="coerce")
    report = report[quantities > 0].copy()
    quantities = quantities.loc[report.index]
    if report.empty:
        return pd.DataFrame(), feature_to_uniprot

    if protein_aggregation == "log2_mean":
        report["_aggregation_value"] = np.log2(quantities)
        matrix = report.pivot_table(
            index="Feature",
            columns="Sample",
            values="_aggregation_value",
            aggfunc="mean",
        )
        matrix = np.exp2(matrix)
    elif protein_aggregation == "linear_mean":
        report["_aggregation_value"] = quantities
        matrix = report.pivot_table(
            index="Feature",
            columns="Sample",
            values="_aggregation_value",
            aggfunc="mean",
        )
    else:
        raise ValueError("protein_aggregation must be 'log2_mean' or 'linear_mean'")
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)
    return matrix, feature_to_uniprot


def preprocess_data(
    pivot_df: pd.DataFrame,
    *,
    min_total_observations: int = 2,
) -> pd.DataFrame:
    """Log2-transform and sample-median center without imputing missing values."""

    matrix = pivot_df.apply(pd.to_numeric, errors="coerce")
    matrix = matrix.where(matrix > 0)
    matrix = matrix.dropna(thresh=min_total_observations)
    log_matrix = np.log2(matrix)
    return log_matrix - log_matrix.median(axis=0, skipna=True)


def impute_data(
    df: pd.DataFrame,
    *,
    seed: int = 42,
    shift: float = 1.8,
    width: float = 0.3,
) -> pd.DataFrame:
    """Deterministically impute for visualization/ssGSEA, never for DE tests."""

    rng = np.random.default_rng(seed)
    imputed = df.copy()
    all_values = imputed.to_numpy(dtype=float)
    global_mean = float(np.nanmean(all_values)) if np.isfinite(all_values).any() else 0.0
    global_std = float(np.nanstd(all_values, ddof=1))
    if not np.isfinite(global_std) or global_std <= 0:
        global_std = 1.0

    for column in imputed.columns:
        observed = imputed[column]
        missing = observed.isna()
        if not missing.any():
            continue
        mean = float(observed.mean())
        std = float(observed.std(ddof=1))
        if not np.isfinite(mean):
            mean = global_mean
        if not np.isfinite(std) or std <= 0:
            std = global_std
        imputed.loc[missing, column] = rng.normal(
            loc=mean - shift * std,
            scale=width * std,
            size=int(missing.sum()),
        )
    return imputed


def _analysis_metadata_path(dp_dir: Path) -> Path:
    return dp_dir / "analysis_metadata.json"


def _analysis_signature(
    cfg: ProteomicsConfig,
    tissue: str,
    report_path: Path | None = None,
) -> dict[str, object]:
    if report_path is None:
        report_path = cfg.data_dirs.get(tissue)
    if report_path is not None and report_path.exists():
        source_stat = report_path.stat()
        source_file: dict[str, object] = {
            "path": str(report_path.resolve()),
            "size": int(source_stat.st_size),
            "mtime_ns": int(source_stat.st_mtime_ns),
        }
    else:
        source_file = {"path": str(report_path) if report_path else ""}
    signature = {
        "pipeline_version": PIPELINE_VERSION,
        "dp_method": cfg.dp_method,
        "fdr_method": cfg.fdr_method,
        "pval_cutoff": cfg.pval_cutoff,
        "log2fc_cutoff": cfg.log2fc_cutoff,
        "protein_aggregation": cfg.protein_aggregation,
        "min_total_observations": cfg.min_total_observations,
        "min_observations_per_group": cfg.min_observations_per_group,
        "imputation_seed": cfg.imputation_seed,
        "imputation_shift": cfg.imputation_shift,
        "imputation_width": cfg.imputation_width,
        "heatmap_top_n": cfg.heatmap_top_n,
        "design_covariates": list(cfg.design_covariates),
        "conditions": cfg.conditions,
        "tissue_label": (cfg.tissue_labels or {}).get(tissue, tissue),
        "plot_colors": cfg.plot_colors or {},
        "enrichment": {
            "run_ssgsea": cfg.run_ssgsea,
            "run_prerank_gsea": cfg.run_prerank_gsea,
            "gene_sets": list(cfg.gene_sets),
            "min_set_size": cfg.gsea_min_set_size,
            "max_set_size": cfg.gsea_max_set_size,
            "permutation_num": cfg.gsea_permutation_num,
            "ssgsea_permutation_num": cfg.ssgsea_permutation_num,
            "seed": cfg.gsea_seed,
            "threads": cfg.gsea_threads,
        },
        "comparisons": [
            list(item) for item in generate_comparisons(cfg.conditions, cfg.contrast_groups)
        ],
        "mouse_groups": dict(sorted(cfg.mouse_groups.items())),
        "sample_metadata": {
            key: dict(sorted(values.items()))
            for key, values in sorted(cfg.sample_metadata.items())
        },
        "source_file": source_file,
    }
    # YAML metadata can contain dates or NumPy scalar-like values. Canonicalize
    # to the exact JSON representation that will be read on a later run.
    return json.loads(json.dumps(signature, sort_keys=True, default=str))


def _load_reusable_results(
    dp_dir: Path,
    cfg: ProteomicsConfig,
    report_path: Path | None = None,
) -> pd.DataFrame | None:
    combined_path = dp_dir / COMBINED_CSV
    if not combined_path.exists():
        return None
    metadata_path = _analysis_metadata_path(dp_dir)
    if not metadata_path.exists():
        raise RuntimeError(
            f"Existing results in {dp_dir} predate the corrected v{PIPELINE_VERSION} "
            "statistics. Re-run with --force; they must not be reused."
        )
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if str(metadata.get("pipeline_version")) != PIPELINE_VERSION:
        raise RuntimeError(
            f"Existing results in {dp_dir} were produced by pipeline "
            f"{metadata.get('pipeline_version')!r}; re-run with --force."
        )
    if metadata.get("fdr_method") != cfg.fdr_method:
        raise RuntimeError(
            f"Existing results use FDR method {metadata.get('fdr_method')!r}, "
            f"but config requests {cfg.fdr_method!r}; re-run with --force."
        )
    current_signature = _analysis_signature(cfg, dp_dir.name, report_path)
    previous_signature = metadata.get("analysis_signature")
    if previous_signature != current_signature:
        changed = []
        previous_signature = previous_signature if isinstance(previous_signature, dict) else {}
        for key in sorted(set(previous_signature) | set(current_signature)):
            if previous_signature.get(key) != current_signature.get(key):
                changed.append(key)
        raise RuntimeError(
            f"Existing results in {dp_dir} do not match the current input/config "
            f"({', '.join(changed) or 'missing signature'}); re-run with --force."
        )
    flush_print(f"  [{dp_dir.name}] corrected DP results exist; reusing them.")
    return pd.read_csv(combined_path)


def _write_analysis_metadata(
    dp_dir: Path,
    cfg: ProteomicsConfig,
    tissue: str,
    normalized_matrix: pd.DataFrame,
    report_path: Path | None = None,
) -> None:
    payload = {
        "pipeline_version": PIPELINE_VERSION,
        "experiment_id": cfg.experiment_id,
        "tissue": tissue,
        "dp_method": cfg.dp_method,
        "fdr_method": cfg.fdr_method,
        "differential_input": "log2_normalized_unimputed",
        "protein_aggregation": cfg.protein_aggregation,
        "min_total_observations": cfg.min_total_observations,
        "min_observations_per_group": cfg.min_observations_per_group,
        "design_covariates": cfg.design_covariates,
        "imputation_scope": "visualization_and_ssgsea_only",
        "imputation_seed": cfg.imputation_seed,
        "n_proteins": int(normalized_matrix.shape[0]),
        "n_samples": int(normalized_matrix.shape[1]),
        "analysis_signature": _analysis_signature(cfg, tissue, report_path),
    }
    _analysis_metadata_path(dp_dir).write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8",
    )


def analyze_tissue(
    tissue: str,
    cfg: ProteomicsConfig,
    *,
    force: bool = False,
    skip_gsea: bool = False,
    report_path: Path | None = None,
) -> None:
    """Run load → QC → model → plots → enrichment for one tissue.

    ``report_path`` is set when this invocation ran DIA-NN. That report is then
    the only input. Otherwise the configured ``paths.data_dirs`` report is used.
    """

    flush_print(f"\n{'=' * 60}\nAnalyzing {tissue}\n{'=' * 60}")
    if report_path is None:
        report_path = cfg.data_dirs.get(tissue)
    if report_path is None:
        flush_print(f"  [SKIP] no data path resolved for {tissue}.")
        return

    output_root = Path(cfg.output_dir)
    dp_dir = output_root / "01_dp_results" / tissue
    plot_dir = output_root / "02_result_plot" / tissue
    qc_dir = output_root / "00_qc" / tissue
    dp_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)

    existing_results = None if force else _load_reusable_results(dp_dir, cfg, report_path)
    comparisons = generate_comparisons(cfg.conditions, cfg.contrast_groups)
    group_colors = get_group_colors(cfg.conditions)

    raw_matrix, feature_to_uniprot = load_and_clean_data(
        str(report_path), protein_aggregation=cfg.protein_aggregation,
    )
    if raw_matrix is None or raw_matrix.empty:
        flush_print(f"  [SKIP] no quantified proteins remained for {tissue}.")
        return

    group_map = get_sample_to_group_map(raw_matrix.columns, cfg.mouse_groups)
    unknown_samples = [sample for sample, group in group_map.items() if group == "Unknown"]
    if unknown_samples:
        raise ValueError(
            f"[{tissue}] samples missing from groups.yaml: {', '.join(unknown_samples)}"
        )
    sample_metadata = build_sample_metadata(
        raw_matrix.columns,
        cfg.mouse_groups,
        cfg.sample_metadata,
        tissue=tissue,
    )
    flush_print(f"  Groups: {group_map}")

    normalized = preprocess_data(
        raw_matrix,
        min_total_observations=cfg.min_total_observations,
    )
    if normalized.empty:
        flush_print(f"  [SKIP] all proteins failed the minimum observation filter.")
        return
    normalized.to_csv(dp_dir / "normalized_matrix.csv")

    imputed = impute_data(
        normalized,
        seed=cfg.imputation_seed,
        shift=cfg.imputation_shift,
        width=cfg.imputation_width,
    )
    imputed.to_csv(dp_dir / "imputed_matrix.csv")

    if force or not (qc_dir / "qc_summary.json").exists():
        flush_print("  Generating QC outputs...")
        generate_qc_outputs(
            raw_matrix,
            normalized,
            sample_metadata,
            qc_dir,
            seed=cfg.imputation_seed,
        )
    else:
        flush_print("  QC outputs exist; skipping QC regeneration.")

    pca_path = plot_dir / "pca_plot.png"
    if (force or not pca_path.exists()) and imputed.shape[0] >= 2 and imputed.shape[1] >= 2:
        plot_pca(imputed, group_map, str(pca_path), group_colors=group_colors)

    if existing_results is None:
        results = run_differential_analysis(
            normalized,
            group_map,
            comparisons,
            str(dp_dir),
            tissue,
            cfg.pval_cutoff,
            cfg.log2fc_cutoff,
            cfg.fdr_method,
            cfg.dp_method,
            feature_to_uniprot=feature_to_uniprot,
            sample_metadata=sample_metadata,
            design_covariates=cfg.design_covariates,
            min_observations_per_group=cfg.min_observations_per_group,
            flush=flush_print,
        )
        _write_analysis_metadata(dp_dir, cfg, tissue, normalized, report_path)
    else:
        results = existing_results

    heatmap_path = plot_dir / f"heatmap_top{cfg.heatmap_top_n}.png"
    if force or not heatmap_path.exists():
        plot_heatmap(
            imputed,
            group_map,
            results,
            str(heatmap_path),
            group_colors=group_colors,
            heatmap_top_n=cfg.heatmap_top_n,
        )

    for _numerator, _denominator, comparison in comparisons:
        volcano_path = plot_dir / f"volcano_{comparison}.png"
        if force or not volcano_path.exists():
            plot_volcano(
                results.copy(),
                comparison,
                str(volcano_path),
                pval_cutoff=cfg.pval_cutoff,
                log2fc_cutoff=cfg.log2fc_cutoff,
                method=cfg.dp_method,
            )

    if skip_gsea:
        flush_print("  Skipping ssGSEA and pre-ranked GSEA (--skip-gsea).")
    else:
        if cfg.run_ssgsea:
            run_ssgsea(
                imputed,
                str(plot_dir),
                gene_sets=cfg.gene_sets,
                min_size=cfg.gsea_min_set_size,
                max_size=cfg.gsea_max_set_size,
                permutation_num=cfg.ssgsea_permutation_num,
                seed=cfg.gsea_seed,
                threads=cfg.gsea_threads,
                force=force,
                flush=flush_print,
            )
        if cfg.run_prerank_gsea:
            run_prerank_gsea(
                results,
                str(plot_dir),
                comparisons,
                gene_sets=cfg.gene_sets,
                min_size=cfg.gsea_min_set_size,
                max_size=cfg.gsea_max_set_size,
                permutation_num=cfg.gsea_permutation_num,
                seed=cfg.gsea_seed,
                threads=cfg.gsea_threads,
                force=force,
                flush=flush_print,
            )
    flush_print(f"  Completed {tissue}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=f"Proteomics analysis pipeline v{PIPELINE_VERSION}",
        epilog=config_arg_help(),
    )
    parser.add_argument("--experiment", default=None, help="Experiment ID")
    parser.add_argument("--config", default=None, help="Path to config.yaml")
    parser.add_argument(
        "--experiment-dir", "--experiment_dir", dest="experiment_dir", default=None,
        help="Deprecated compatibility alias: proteomics_data directory containing config.yaml",
    )
    parser.add_argument(
        "--outdir", "--output-dir", "--output_dir", dest="outdir", default=None,
        help="Override output directory",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Recompute statistics. Re-runs DIA-NN only together with --search or --search-only",
    )
    parser.add_argument(
        "--skip-gsea", action="store_true",
        help="Skip both sample-level ssGSEA and comparison-level pre-ranked GSEA",
    )
    parser.add_argument(
        "--tissues", nargs="+", default=None,
        help="Optional subset of configured tissues (useful for validation)",
    )
    parser.add_argument(
        "--search", action="store_true",
        help="Run DIA-NN on Thermo .raw files, then the existing differential analysis",
    )
    parser.add_argument(
        "--search-only", action="store_true",
        help="Run DIA-NN only; do not run QC or differential analysis",
    )
    args = parser.parse_args()
    if args.search and args.search_only:
        parser.error("Use either --search or --search-only, not both.")

    resolved_config = args.config
    if args.experiment_dir:
        if args.experiment or args.config:
            parser.error("--experiment-dir cannot be combined with --experiment/--config")
        experiment_dir = Path(args.experiment_dir).resolve()
        candidate = experiment_dir if experiment_dir.is_file() else experiment_dir / "config.yaml"
        if not candidate.exists():
            parser.error(f"No config.yaml found under --experiment-dir: {experiment_dir}")
        resolved_config = str(candidate)

    cfg = load_proteomics_config(
        experiment_id=args.experiment,
        config_path=resolved_config,
        outdir_override=args.outdir,
    )
    tissues = list(cfg.tissues)
    if args.tissues:
        unknown = sorted(set(args.tissues) - set(tissues))
        if unknown:
            parser.error("Unknown tissues: " + ", ".join(unknown))
        tissues = list(args.tissues)

    do_search = bool(args.search or args.search_only)
    if do_search and not cfg.search_enabled:
        parser.error(
            "--search and --search-only require search.enabled: true, plus "
            "search.diann_exe, search.fasta, and paths.raw_dirs. "
            "--force alone does not run DIA-NN."
        )

    comparisons = generate_comparisons(cfg.conditions, cfg.contrast_groups)
    flush_print(f"Experiment: {cfg.experiment_id}")
    flush_print(f"Output: {cfg.output_dir}")
    flush_print(f"Tissues: {tissues}")
    flush_print(f"Comparisons: {[comparison[2] for comparison in comparisons]}")
    flush_print(
        f"Statistics: {cfg.dp_method}, FDR={cfg.fdr_method}, "
        f"covariates={cfg.design_covariates or 'none'}"
    )
    if do_search:
        flush_print(
            f"Search: DIA-NN {cfg.diann_expected_version}, "
            f"qvalue={cfg.search_qvalue}, force={bool(args.force)}"
        )

    for tissue in tissues:
        # --force recomputes statistics. It re-runs DIA-NN only together with
        # --search or --search-only; a stats-only --force must not repeat a search.
        searched_report = None
        if do_search:
            searched_report = run_tissue_search(tissue, cfg, force=args.force)
        if args.search_only:
            continue
        analyze_tissue(
            tissue,
            cfg,
            force=args.force,
            skip_gsea=args.skip_gsea,
            report_path=searched_report,
        )
    flush_print("\nAll tissues completed.")


if __name__ == "__main__":
    main()
