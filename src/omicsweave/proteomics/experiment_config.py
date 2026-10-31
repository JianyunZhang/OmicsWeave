"""
Proteomics experiment configuration loader.

Mirrors the transcriptomics config pattern (omicsweave/transcriptomics/experiment_config.py).
Provides ProteomicsConfig dataclass, config loader with path resolution, and helper functions
for comparisons, colors, and labels.

Usage:
    from experiment_config import load_proteomics_config, generate_comparisons
    cfg = load_proteomics_config(experiment_id="20250425_cimr_wt_mice_dio")
"""
from __future__ import annotations

import glob
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Sequence

import pandas as pd
import yaml

from omicsweave.core.project import project_dir


# ============================================================================
# Dataclass
# ============================================================================

@dataclass(frozen=True)
class ProteomicsConfig:
    experiment_id: str
    config_path: Path
    output_dir: Path
    # tissue_key -> resolved report TSV path
    data_dirs: dict[str, Path]
    tissues: list[str]
    groups_path: Path
    conditions: dict[str, dict[str, Any]]
    mouse_groups: dict[str, str]
    # Parameters
    pval_cutoff: float
    log2fc_cutoff: float
    heatmap_top_n: int
    fdr_method: str
    dp_method: str  # "limma" or "wilcoxon"
    min_total_observations: int
    min_observations_per_group: int
    imputation_seed: int
    imputation_shift: float
    imputation_width: float
    protein_aggregation: str
    design_covariates: list[str]
    sample_metadata: dict[str, dict[str, Any]]
    # GSEA parameters
    gsea_min_set_size: int
    gsea_max_set_size: int
    gsea_permutation_num: int
    ssgsea_permutation_num: int
    gsea_seed: int
    gsea_threads: int
    run_ssgsea: bool
    run_prerank_gsea: bool
    gene_sets: list[str]  # gene set library names or shorthand keys
    # Optional display
    tissue_labels: dict[str, str] | None = None
    plot_colors: dict[str, Any] | None = None
    contrast_groups: dict[str, list[str]] | None = None
    # Optional Thermo .raw → DIA-NN search. Disabled configs keep reading data_dirs.
    search_enabled: bool = False
    search_engine: str = "diann"
    diann_expected_version: str = "1.9"
    diann_exe: Path | None = None
    search_fasta: Path | None = None
    search_threads: int = 8
    search_qvalue: float = 0.01
    search_reanalyse: bool = True
    search_smart_profiling: bool = True
    raw_dirs: dict[str, Path] = field(default_factory=dict)


# ============================================================================
# Path utilities (same logic as transcriptomics)
# ============================================================================

def repo_root() -> Path:
    """项目目录（实验配置与结果所在仓库），见 :func:`omicsweave.core.project.project_dir`。"""
    return project_dir()


def _is_abs_path_like(p: str) -> bool:
    if not p:
        return False
    return p.startswith("\\\\") or bool(re.match(r"^[A-Za-z]:[\\/]", p))


def _resolve_path(p: str, base_dir: Path) -> Path:
    if _is_abs_path_like(p):
        return Path(p)
    pr = repo_root() / p
    if pr.exists():
        return pr
    return (base_dir / p).resolve()


def _load_yaml(path: Path) -> dict[str, Any]:
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"Invalid YAML root (expected mapping): {path}")
    return data


def _substitute_vars(s: str, vars_map: dict[str, str]) -> str:
    out = s
    for k, v in vars_map.items():
        out = out.replace("${" + k + "}", v)
    return out


_FDR_ALIASES = {
    "bh": "fdr_bh",
    "benjamini-hochberg": "fdr_bh",
    "benjamini_hochberg": "fdr_bh",
    "by": "fdr_by",
    "benjamini-yekutieli": "fdr_by",
    "benjamini_yekutieli": "fdr_by",
}

_VALID_FDR_METHODS = {
    "bonferroni",
    "sidak",
    "holm-sidak",
    "holm",
    "simes-hochberg",
    "hommel",
    "fdr_bh",
    "fdr_by",
    "fdr_tsbh",
    "fdr_tsbky",
}


def normalize_fdr_method(method: str) -> str:
    """Normalize documented aliases to statsmodels names and reject typos.

    Invalid multiple-testing methods are configuration errors.  They must not
    silently fall back to unadjusted p-values.
    """

    normalized = str(method).strip().lower()
    normalized = _FDR_ALIASES.get(normalized, normalized)
    if normalized not in _VALID_FDR_METHODS:
        supported = ", ".join(sorted(_VALID_FDR_METHODS))
        raise ValueError(
            f"Unsupported FDR method '{method}'. Supported methods: {supported}"
        )
    return normalized


def _as_bool(value: Any, *, name: str) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in {"true", "yes", "1", "on"}:
            return True
        if lowered in {"false", "no", "0", "off"}:
            return False
    if isinstance(value, (int, float)) and value in {0, 1}:
        return bool(value)
    raise ValueError(f"{name} must be a boolean, got {value!r}")


# ============================================================================
# Config loader
# ============================================================================

def default_config_path(experiment_id: str) -> Path:
    root = repo_root() / "experiments" / experiment_id
    p1 = root / "proteomics_data" / "config.yaml"
    if p1.exists():
        return p1
    return root / "proteomics" / "config.yaml"


def _find_report_tsv(folder: Path) -> Path:
    """Find the main report TSV inside a DIA-NN output folder.

    Excludes *.pg_matrix.tsv and *.pr_matrix.tsv.
    """
    candidates = sorted(folder.glob("report_*.tsv"))
    candidates = [
        c for c in candidates
        if not c.name.endswith(".pg_matrix.tsv")
        and not c.name.endswith(".pr_matrix.tsv")
    ]
    if not candidates:
        raise FileNotFoundError(
            f"No report_*.tsv found in {folder} "
            f"(found {list(folder.glob('report_*.tsv'))})"
        )
    if len(candidates) > 1:
        raise ValueError(
            f"Multiple main report TSV files found in {folder}: "
            + ", ".join(candidate.name for candidate in candidates)
            + ". Point paths.data_dirs directly to the intended TSV file."
        )
    return candidates[0]


def load_proteomics_config(
    *,
    experiment_id: str | None = None,
    config_path: str | None = None,
    outdir_override: str | None = None,
) -> ProteomicsConfig:
    """Load and resolve proteomics experiment configuration."""
    exp = experiment_id or os.environ.get("EXPERIMENT_ID")
    if config_path:
        cfg_path = Path(config_path).resolve()
    else:
        if not exp:
            raise ValueError(
                "Missing experiment id/config. "
                "Please pass --experiment <id> or --config <path>."
            )
        cfg_path = default_config_path(exp)
    cfg = _load_yaml(cfg_path)
    base_dir = cfg_path.parent

    # Allow config to define/override experiment_id
    if not exp:
        exp = str(cfg.get("experiment_id", "")).strip() or None
    exp = exp or str(cfg.get("experiment_id", "")).strip() or cfg_path.parents[1].name

    paths = cfg.get("paths", {}) or {}
    params = cfg.get("params", {}) or {}

    # ---- output_dir ----
    output_dir_s = str(paths.get("output_dir", "")).strip()
    if outdir_override:
        output_dir_s = outdir_override
    if not output_dir_s:
        raise ValueError(f"paths.output_dir missing in {cfg_path}")
    output_dir_s = _substitute_vars(output_dir_s, {"experiment_id": exp})
    if _is_abs_path_like(output_dir_s):
        output_dir = Path(output_dir_s)
    else:
        # Output paths are a project-level contract and are always relative to
        # the repository root. Input data_dirs remain relative to config.yaml.
        output_dir = (repo_root() / output_dir_s).resolve()

    # ---- data_dirs (tissue -> resolved report TSV path) ----
    data_dirs_raw = paths.get("data_dirs", None)
    data_dirs: dict[str, Path] = {}
    if data_dirs_raw and isinstance(data_dirs_raw, dict):
        # New nested format
        for tissue, folder_rel in data_dirs_raw.items():
            folder_path = _resolve_path(str(folder_rel), base_dir)
            if not folder_path.exists():
                print(f"Warning: data folder not found for {tissue}: {folder_path}")
                continue
            report_tsv = folder_path if folder_path.is_file() else _find_report_tsv(folder_path)
            data_dirs[str(tissue)] = report_tsv
    else:
        # Fallback: legacy flat format (sample_data_<tissue>_folder)
        for key, val in paths.items():
            if key.startswith("sample_data_") and key.endswith("_folder"):
                tissue_key = key.replace("sample_data_", "").replace("_folder", "")
                folder_path = _resolve_path(str(val), base_dir)
                if folder_path.exists():
                    data_dirs[tissue_key] = (
                        folder_path if folder_path.is_file() else _find_report_tsv(folder_path)
                    )

    # ---- tissues list ----
    tissues = list(params.get("tissues", []))
    if not tissues:
        tissues = sorted(data_dirs.keys())
    if not tissues:
        raise ValueError(f"params.tissues missing/empty and no data_dirs in {cfg_path}")

    # ---- parameters ----
    preprocessing = params.get("preprocessing", {}) or {}
    differential_design = params.get("differential_design", {}) or {}
    enrichment = params.get("enrichment", {}) or {}
    if not isinstance(preprocessing, dict):
        raise ValueError("params.preprocessing must be a mapping")
    if not isinstance(differential_design, dict):
        raise ValueError("params.differential_design must be a mapping")
    if not isinstance(enrichment, dict):
        raise ValueError("params.enrichment must be a mapping")

    pval_cutoff = float(params.get("pval_cutoff", 0.05))
    log2fc_cutoff = float(params.get("log2fc_cutoff", 0.58))
    heatmap_top_n = int(params.get("heatmap_top_n", 50))
    fdr_method = normalize_fdr_method(str(params.get("fdr_method", "fdr_bh")))
    dp_method = str(params.get("dp_method", params.get("de_method", "limma"))).lower()
    if dp_method not in ("limma", "wilcoxon"):
        raise ValueError(f"Unsupported dp_method '{dp_method}'. Use 'limma' or 'wilcoxon'.")
    min_total_observations = int(
        preprocessing.get("min_total_observations", params.get("min_total_observations", 2))
    )
    min_observations_per_group = int(
        differential_design.get(
            "min_observations_per_group",
            params.get("min_observations_per_group", 2),
        )
    )
    imputation_seed = int(
        preprocessing.get("imputation_seed", params.get("imputation_seed", 42))
    )
    imputation_shift = float(
        preprocessing.get("imputation_shift", params.get("imputation_shift", 1.8))
    )
    imputation_width = float(
        preprocessing.get("imputation_width", params.get("imputation_width", 0.3))
    )
    protein_aggregation = str(
        preprocessing.get(
            "protein_aggregation",
            params.get("protein_aggregation", "log2_mean"),
        )
    ).strip().lower()
    if protein_aggregation not in {"log2_mean", "linear_mean"}:
        raise ValueError(
            "protein_aggregation must be 'log2_mean' or 'linear_mean'"
        )

    covariates_raw = differential_design.get(
        "covariates", params.get("design_covariates", [])
    )
    if isinstance(covariates_raw, str):
        design_covariates = [covariates_raw]
    elif isinstance(covariates_raw, list):
        design_covariates = [str(value) for value in covariates_raw]
    else:
        raise ValueError("design covariates must be a string or list")
    reserved_covariates = {"condition", "sample", "mouse_id", "tissue"}
    invalid_covariates = sorted(set(design_covariates) & reserved_covariates)
    if invalid_covariates:
        raise ValueError(
            "design covariates contain reserved columns: "
            + ", ".join(invalid_covariates)
        )
    if dp_method == "wilcoxon" and design_covariates:
        raise ValueError(
            "dp_method='wilcoxon' cannot adjust covariates; use dp_method='limma'"
        )

    gsea_min_set_size = int(
        enrichment.get("min_set_size", params.get("gsea_min_set_size", 10))
    )
    gsea_max_set_size = int(
        enrichment.get("max_set_size", params.get("gsea_max_set_size", 500))
    )
    gsea_permutation_num = int(
        enrichment.get("permutation_num", params.get("gsea_permutation_num", 1000))
    )
    ssgsea_permutation_num = int(
        enrichment.get(
            "ssgsea_permutation_num",
            params.get("ssgsea_permutation_num", 0),
        )
    )
    gsea_seed = int(enrichment.get("seed", params.get("gsea_seed", imputation_seed)))
    gsea_threads = int(enrichment.get("threads", params.get("gsea_threads", 1)))
    run_ssgsea = _as_bool(
        enrichment.get("run_ssgsea", params.get("run_ssgsea", True)),
        name="run_ssgsea",
    )
    run_prerank_gsea = _as_bool(
        enrichment.get("run_prerank_gsea", params.get("run_prerank_gsea", True)),
        name="run_prerank_gsea",
    )
    gene_sets_raw = params.get("gene_sets", ["KEGG_2019_Mouse"])
    gene_sets_raw = enrichment.get("gene_sets", gene_sets_raw)
    gene_sets = [str(g) for g in (gene_sets_raw if isinstance(gene_sets_raw, list) else [gene_sets_raw])]

    if not 0 < pval_cutoff <= 1:
        raise ValueError("pval_cutoff must be in (0, 1]")
    if log2fc_cutoff < 0:
        raise ValueError("log2fc_cutoff must be >= 0")
    if heatmap_top_n < 1:
        raise ValueError("heatmap_top_n must be >= 1")
    if min_total_observations < 1:
        raise ValueError("min_total_observations must be >= 1")
    if min_observations_per_group < 2:
        raise ValueError("min_observations_per_group must be >= 2")
    if imputation_shift < 0 or imputation_width <= 0:
        raise ValueError("imputation_shift must be >= 0 and imputation_width > 0")
    if gsea_min_set_size < 1 or gsea_max_set_size < gsea_min_set_size:
        raise ValueError("invalid GSEA gene-set size limits")
    if gsea_permutation_num < 1 or ssgsea_permutation_num < 0:
        raise ValueError("GSEA permutations must be >=1; ssGSEA permutations must be >=0")
    if gsea_threads < 1:
        raise ValueError("gsea_threads must be >= 1")

    # ---- groups ----
    groups_file = str(cfg.get("groups_file", "groups.yaml")).strip() or "groups.yaml"
    groups_path = _resolve_path(groups_file, base_dir)
    groups = _load_yaml(groups_path)

    conditions = groups.get("conditions", {}) or {}
    mouse_groups = groups.get("mouse_groups", {}) or {}
    mouse_groups = {str(k): str(v) for k, v in mouse_groups.items()}

    sample_metadata_raw = groups.get("sample_metadata", {}) or {}
    if not isinstance(sample_metadata_raw, dict):
        raise ValueError(f"sample_metadata must be a mapping in {groups_path}")
    sample_metadata: dict[str, dict[str, Any]] = {}
    for sample_id, values in sample_metadata_raw.items():
        if not isinstance(values, dict):
            raise ValueError(
                f"sample_metadata[{sample_id!r}] must be a mapping in {groups_path}"
            )
        normalized_values = {str(k): v for k, v in values.items()}
        reserved_metadata = sorted(
            set(normalized_values) & {"sample", "mouse_id", "condition", "tissue"}
        )
        if reserved_metadata:
            raise ValueError(
                f"sample_metadata[{sample_id!r}] uses reserved fields: "
                + ", ".join(reserved_metadata)
            )
        sample_metadata[str(sample_id)] = normalized_values

    sample_batches = groups.get("sample_batches", {}) or groups.get("mouse_batches", {}) or {}
    if not isinstance(sample_batches, dict):
        raise ValueError(f"sample_batches must be a mapping in {groups_path}")
    for sample_id, batch in sample_batches.items():
        sample_metadata.setdefault(str(sample_id), {})["batch"] = str(batch)

    contrast_groups_raw = groups.get("contrast_groups", None)
    contrast_groups: dict[str, list[str]] | None = None
    if contrast_groups_raw:
        contrast_groups = {str(k): [str(v) for v in vs] for k, vs in contrast_groups_raw.items()}

    if not conditions:
        raise ValueError(f"conditions missing/empty in {groups_path}")
    if not mouse_groups:
        raise ValueError(f"mouse_groups missing/empty in {groups_path}")
    unknown_condition_codes = sorted(set(mouse_groups.values()) - set(conditions))
    if unknown_condition_codes:
        raise ValueError(
            f"mouse_groups references undefined conditions in {groups_path}: "
            + ", ".join(unknown_condition_codes)
        )

    # ---- optional display ----
    tissue_labels = params.get("tissue_labels", None)
    if tissue_labels:
        tissue_labels = {str(k): str(v) for k, v in tissue_labels.items()}

    plot_colors = params.get("plot_colors", None)

    search_cfg = cfg.get("search", {}) or {}
    if not isinstance(search_cfg, dict):
        raise ValueError(f"search must be a mapping in {cfg_path}")
    search_enabled = _as_bool(search_cfg.get("enabled", False), name="search.enabled")
    search_engine = str(search_cfg.get("engine", "diann")).strip().lower() or "diann"
    diann_expected_version = str(search_cfg.get("expected_version", "1.9")).strip() or "1.9"
    search_threads = int(search_cfg.get("threads", 8))
    search_qvalue = float(search_cfg.get("qvalue", 0.01))
    search_reanalyse = _as_bool(
        search_cfg.get("reanalyse", True), name="search.reanalyse",
    )
    search_smart_profiling = _as_bool(
        search_cfg.get("smart_profiling", True), name="search.smart_profiling",
    )
    if search_threads < 1:
        raise ValueError("search.threads must be >= 1")
    if not 0 < search_qvalue <= 1:
        raise ValueError("search.qvalue must be in (0, 1]")

    raw_dirs_raw = paths.get("raw_dirs", {}) or {}
    if not isinstance(raw_dirs_raw, dict):
        raise ValueError(f"paths.raw_dirs must be a mapping in {cfg_path}")
    raw_dirs: dict[str, Path] = {
        str(tissue): _resolve_path(str(folder), base_dir)
        for tissue, folder in raw_dirs_raw.items()
    }

    diann_exe: Path | None = None
    search_fasta: Path | None = None
    exe_raw = str(search_cfg.get("diann_exe", "")).strip()
    fasta_raw = str(search_cfg.get("fasta", "")).strip()
    if exe_raw:
        diann_exe = _resolve_path(exe_raw, base_dir)
    if fasta_raw:
        search_fasta = _resolve_path(fasta_raw, base_dir)
    if search_enabled:
        if search_engine != "diann":
            raise ValueError(
                f"Unsupported search.engine {search_engine!r}. "
                "This pipeline calls DIA-NN."
            )
        if diann_exe is None or not diann_exe.is_file():
            raise ValueError(
                "search.diann_exe must point to an existing DIA-NN executable "
                f"when search.enabled is true (got {diann_exe})"
            )
        if search_fasta is None or not search_fasta.is_file():
            raise ValueError(
                "search.fasta must point to an existing FASTA file "
                f"when search.enabled is true (got {search_fasta})"
            )
        if not raw_dirs:
            raise ValueError(
                "paths.raw_dirs is required when search.enabled is true"
            )

    return ProteomicsConfig(
        experiment_id=exp,
        config_path=cfg_path,
        output_dir=output_dir,
        data_dirs=data_dirs,
        tissues=tissues,
        groups_path=groups_path,
        conditions=conditions,
        mouse_groups=mouse_groups,
        contrast_groups=contrast_groups,
        pval_cutoff=pval_cutoff,
        log2fc_cutoff=log2fc_cutoff,
        heatmap_top_n=heatmap_top_n,
        fdr_method=fdr_method,
        dp_method=dp_method,
        min_total_observations=min_total_observations,
        min_observations_per_group=min_observations_per_group,
        imputation_seed=imputation_seed,
        imputation_shift=imputation_shift,
        imputation_width=imputation_width,
        protein_aggregation=protein_aggregation,
        design_covariates=design_covariates,
        sample_metadata=sample_metadata,
        gsea_min_set_size=gsea_min_set_size,
        gsea_max_set_size=gsea_max_set_size,
        gsea_permutation_num=gsea_permutation_num,
        ssgsea_permutation_num=ssgsea_permutation_num,
        gsea_seed=gsea_seed,
        gsea_threads=gsea_threads,
        run_ssgsea=run_ssgsea,
        run_prerank_gsea=run_prerank_gsea,
        gene_sets=gene_sets,
        tissue_labels=tissue_labels,
        plot_colors=plot_colors,
        search_enabled=search_enabled,
        search_engine=search_engine,
        diann_expected_version=diann_expected_version,
        diann_exe=diann_exe,
        search_fasta=search_fasta,
        search_threads=search_threads,
        search_qvalue=search_qvalue,
        search_reanalyse=search_reanalyse,
        search_smart_profiling=search_smart_profiling,
        raw_dirs=raw_dirs,
    )


# ============================================================================
# Sample name utilities
# ============================================================================

def extract_sample_name(run_name: str) -> str:
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


def get_sample_to_mouse_id_map(
    columns: list[str] | pd.Index,
    known_mouse_ids: Sequence[str] | dict[str, Any],
) -> dict[str, str | None]:
    """Resolve sample names to configured animal IDs, including IDs like ``A35``.

    The configured IDs are matched as case-insensitive suffixes, longest first.
    Longest-first matching prevents ``A35`` from being accidentally reduced to
    ``35`` when both IDs exist in a project.
    """

    ids = list(known_mouse_ids.keys()) if isinstance(known_mouse_ids, dict) else list(known_mouse_ids)
    ordered_ids = sorted((str(value) for value in ids), key=lambda value: (-len(value), value))
    mapping: dict[str, str | None] = {}
    for column in columns:
        sample = str(column)
        folded = sample.casefold()
        matches = [mouse_id for mouse_id in ordered_ids if folded.endswith(mouse_id.casefold())]
        mapping[sample] = matches[0] if matches else None
    return mapping


def get_sample_to_group_map(
    columns: list[str] | pd.Index,
    mouse_groups: dict[str, str],
) -> dict[str, str]:
    """Map sample names to condition codes using configured animal IDs."""

    mouse_ids = get_sample_to_mouse_id_map(columns, mouse_groups)
    return {
        sample: mouse_groups.get(mouse_id, "Unknown") if mouse_id else "Unknown"
        for sample, mouse_id in mouse_ids.items()
    }


def build_sample_metadata(
    columns: list[str] | pd.Index,
    mouse_groups: dict[str, str],
    configured_metadata: dict[str, dict[str, Any]] | None = None,
    *,
    tissue: str | None = None,
) -> pd.DataFrame:
    """Build an aligned sample metadata table for QC and model covariates."""

    samples = [str(column) for column in columns]
    mouse_ids = get_sample_to_mouse_id_map(samples, mouse_groups)
    configured_metadata = configured_metadata or {}
    rows: list[dict[str, Any]] = []
    for sample in samples:
        mouse_id = mouse_ids[sample]
        condition = mouse_groups.get(mouse_id, "Unknown") if mouse_id else "Unknown"
        row: dict[str, Any] = {
            "sample": sample,
            "mouse_id": mouse_id,
            "condition": condition,
        }
        if tissue is not None:
            row["tissue"] = tissue
        if mouse_id is not None:
            row.update(configured_metadata.get(mouse_id, {}))
        # Exact sample-level metadata overrides animal-level defaults, enabling
        # technical replicates from one animal to carry different batch labels.
        row.update(configured_metadata.get(sample, {}))
        rows.append(row)
    return pd.DataFrame(rows).set_index("sample")


# ============================================================================
# Helper functions (same pattern as transcriptomics)
# ============================================================================

_DEFAULT_PALETTE = [
    "#4C72B0", "#DD8452", "#55A868", "#C44E52",
    "#8172B3", "#937860", "#DA8BC3", "#8C8C8C",
]


def generate_comparisons(
    conditions: dict[str, Any],
    contrast_groups: dict[str, list[str]] | None = None,
) -> list[tuple[str, str, str]]:
    """Generate pairwise comparison list from conditions.

    If ``contrast_groups`` is provided (from groups.yaml), use user-specified
    comparisons only.  Otherwise auto-generate all pairwise comparisons.

    Returns list of (numerator, denominator, comp_name) tuples.
    """
    if contrast_groups:
        comps: list[tuple[str, str, str]] = []
        for label, pair in contrast_groups.items():
            if len(pair) != 2:
                raise ValueError(
                    f"contrast_groups[{label!r}] must contain exactly two conditions"
                )
            a, b = str(pair[0]), str(pair[1])
            undefined = sorted({a, b} - set(conditions))
            if undefined:
                raise ValueError(
                    f"contrast_groups[{label!r}] references undefined conditions: "
                    + ", ".join(undefined)
                )
            if a == b:
                raise ValueError(f"contrast_groups[{label!r}] compares a condition to itself")
            comps.append((a, b, f"{a}_vs_{b}"))
        return comps

    keys = sorted(conditions.keys())
    comps = []
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            comps.append((a, b, f"{a}_vs_{b}"))
    return comps


def get_group_colors(conditions: dict[str, Any]) -> dict[str, str]:
    """Build group -> color mapping from conditions' `color` field."""
    colors = {}
    for i, (key, val) in enumerate(sorted(conditions.items())):
        if isinstance(val, dict) and val.get("color"):
            colors[key] = val["color"]
        else:
            colors[key] = _DEFAULT_PALETTE[i % len(_DEFAULT_PALETTE)]
    return colors


def get_tissue_labels(cfg: ProteomicsConfig) -> dict[str, str]:
    """Return tissue display labels."""
    if cfg.tissue_labels:
        return dict(cfg.tissue_labels)
    return {t: t.replace("_", " ").title() for t in cfg.tissues}


def get_comp_labels(
    conditions: dict[str, Any],
    comparisons: list[tuple[str, str, str]],
) -> dict[str, str]:
    """Build comp_name -> display label mapping."""
    def _disp(code: str) -> str:
        val = conditions.get(code)
        if isinstance(val, dict) and val.get("display"):
            return val["display"]
        return code

    return {name: f"{_disp(a)} vs {_disp(b)}" for a, b, name in comparisons}


def get_cond_map(conditions: dict[str, Any]) -> dict[str, str]:
    """Build condition code -> display name mapping."""
    result = {}
    for key, val in conditions.items():
        if isinstance(val, dict) and val.get("display"):
            result[key] = val["display"]
        else:
            result[key] = key
    return result


_PLOT_COLOR_DEFAULTS: dict[str, Any] = {
    "heatmap_cmap": "RdBu_r",
    "heatmap_vmin": -2,
    "heatmap_vmax": 2,
    "volcano_up": "#D32F2F",
    "volcano_down": "#1565C0",
    "volcano_ns": "#BDBDBD",
    "gsea_positive": "#D32F2F",
    "gsea_negative": "#1565C0",
    "enrichment_cmap": "RdYlBu_r",
    "pca_cmap": "Set2",
}


def get_plot_colors(cfg: ProteomicsConfig) -> dict[str, Any]:
    """Return plot color configuration with defaults."""
    raw = cfg.plot_colors or {}
    return {**_PLOT_COLOR_DEFAULTS, **raw}


def config_arg_help() -> str:
    return (
        "Use experiment config under experiments/<experiment>/proteomics_data/config.yaml. "
        "Pass --experiment <id> or --config <path>. Legacy --experiment_dir and "
        "--output_dir remain compatibility aliases."
    )
