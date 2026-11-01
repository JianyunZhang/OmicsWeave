#!/usr/bin/env python3
"""
Multi-omics integration configuration loader.

Self-contained config system for multi-omics integration analysis.
All parameters (conditions, tissues, comparisons, scoring, etc.) are defined
independent of other omics pipelines.

Config file location: experiments/{experiment_id}/multiomics_integration_setting/config.yaml

Author: Claude Code
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from omicsweave.core.project import project_dir


# ---------------------------------------------------------------------------
# Path resolution utilities (same pattern as transcriptomics experiment_config.py)
# ---------------------------------------------------------------------------

def repo_root() -> Path:
    """项目目录（实验配置与结果所在仓库），见 :func:`omicsweave.core.project.project_dir`。"""
    return project_dir()


def _is_abs_path_like(p: str) -> bool:
    """Check if a path string is absolute (UNC or Windows drive)."""
    if not p:
        return False
    return p.startswith("\\\\") or bool(re.match(r"^[A-Za-z]:[\\/]", p))


def _resolve_path(p: str, base_dir: Path, *, allow_create: bool = False) -> Path:
    """Resolve a path: absolute -> repo-relative -> config-relative.

    Parameters
    ----------
    allow_create : bool
        If True, prefer repo_root() path without checking existence.
        Use for output directories that may not yet exist.
        If False (default), only use repo_root() path when it already exists.
    """
    if _is_abs_path_like(p):
        return Path(p)
    pr = repo_root() / p
    if allow_create or pr.exists():
        return pr
    return (base_dir / p).resolve()


def _resolve_path_or_none(p: str | None, base_dir: Path, *, allow_create: bool = False) -> Path | None:
    """Resolve a path that may be None."""
    if p is None or p == "null" or p == "":
        return None
    return _resolve_path(p, base_dir, allow_create=allow_create)


def _load_yaml(path: Path) -> dict[str, Any]:
    """Load and validate a YAML file."""
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"Invalid YAML root (expected mapping): {path}")
    return data


def _substitute_vars(s: str, vars_map: dict[str, str]) -> str:
    """Replace ${var_name} placeholders in a string."""
    out = s
    for k, v in vars_map.items():
        out = out.replace("${" + k + "}", v)
    return out


def _norm_tissue(name: str) -> str:
    """Normalize tissue name to lowercase_with_underscores."""
    return name.strip().lower().replace(" ", "_").replace("-", "_")


# ---------------------------------------------------------------------------
# Default scoring parameters
# ---------------------------------------------------------------------------

_DEFAULT_SCORING = {
    "pval_cutoff": 0.05,
    "log2fc_cutoff": 0.58,
    "effect_tiers": {
        "log2fc": [0.58, 1.0],
        "cohen_d": [0.8, 1.5],
    },
    # Per-omics significance thresholds (override global pval_cutoff / log2fc_cutoff)
    # If an omics key is present, its pval_cutoff and log2fc_cutoff are used;
    # otherwise the global defaults apply.
    "omics_cutoffs": {
        "mfa": {
            "pval_cutoff": 0.05,
            "log2fc_cutoff": 0.58,
        },
        "transcriptomics": {
            "pval_cutoff": 0.05,
            "log2fc_cutoff": 0.58,
        },
        "proteomics": {
            "pval_cutoff": 0.05,
            "log2fc_cutoff": 0.58,
        },
        "metabolomics": {
            "pval_cutoff": 0.05,
            "log2fc_cutoff": 0.58,
        },
    },
}

_DEFAULT_PLOT_COLORS = {
    "omics": {
        "transcriptomics": "#4C72B0",
        "metabolomics": "#55A868",
        "proteomics": "#C44E52",
        "mfa": "#8172B2",
    },
    "evidence_grade": {
        "Replicated": "#1B7837",
        "Convergent": "#5AAE61",
        "Single-layer": "#A6DBA0",
        "Exploratory": "#F4A582",
        "Conflicting": "#B2182B",
        "No-signal": "#D9D9D9",
        "Not-assessable": "#F0F0F0",
    },
    "tissue": {
        "liver": "#E24A33",
        "heart": "#348ABD",
        "kidney": "#988ED5",
        "brain": "#FBC15E",
        "gast_muscle": "#8EBA42",
        "sq_adipose": "#FFB5B8",
    },
}


# ---------------------------------------------------------------------------
# Config dataclass
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class MultiOmicsConfig:
    """Configuration for multi-omics integration analysis.

    All parameters are self-contained — no dependency on other omics configs.
    """

    experiment_id: str
    config_path: Path
    output_dir: Path

    # Data source directories (all optional — load whatever is available)
    mfa_results: Path | None
    metabolomics_results: Path | None
    proteomics_results: Path | None
    transcriptomics_results: Path | None

    # Experimental conditions (self-defined, not reused from other omics)
    conditions: dict[str, dict[str, Any]]
    mouse_groups: dict[str, str]

    # Tissue and comparison configuration
    tissues: list[str]
    comparisons: list[str]
    mfa_tissue_map: dict[str, str]
    comparison_map: dict[str, str]
    transcriptomics_comparison_map: dict[str, str] | None
    proteomics_comparison_map: dict[str, str] | None
    metabolomics_comparison_map: dict[str, str] | None

    # Estimand and multi-batch input contract (optional, backward compatible)
    estimand: dict[str, Any]
    integration: dict[str, Any]
    harmonization: dict[str, Any]

    # Scoring parameters
    scoring: dict[str, Any]

    # Display
    tissue_labels: dict[str, str] | None
    plot_colors: dict[str, Any] | None


# ---------------------------------------------------------------------------
# Config discovery
# ---------------------------------------------------------------------------

def default_config_path(experiment_id: str) -> Path:
    """Find config.yaml for a given experiment ID."""
    root = repo_root() / "experiments" / experiment_id
    return root / "multiomics_integration_setting" / "config.yaml"


# ---------------------------------------------------------------------------
# Helper functions for derived metadata
# ---------------------------------------------------------------------------

def generate_comparisons(
    conditions: dict[str, dict[str, Any]],
) -> list[tuple[str, str, str]]:
    """Generate all pairwise (cond_a, cond_b, a_vs_b) comparison tuples."""
    keys = sorted(conditions.keys())
    result = []
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            result.append((keys[j], keys[i], f"{keys[j]}_vs_{keys[i]}"))
    return result


def get_group_colors(conditions: dict[str, dict[str, Any]]) -> dict[str, str]:
    """Map condition codes to colors, using config or built-in palette."""
    _DEFAULT_COLORS = [
        "#E24A33", "#348ABD", "#988ED5", "#FBC15E",
        "#8EBA42", "#FFB5B8", "#777777", "#CCCCCC",
    ]
    colors = {}
    for i, (code, info) in enumerate(sorted(conditions.items())):
        colors[code] = info.get("color", _DEFAULT_COLORS[i % len(_DEFAULT_COLORS)])
    return colors


def get_tissue_labels(
    cfg: MultiOmicsConfig,
) -> dict[str, str]:
    """Return display labels for tissues."""
    if cfg.tissue_labels:
        return dict(cfg.tissue_labels)
    return {t: t.replace("_", " ").title() for t in cfg.tissues}


def get_comp_labels(
    conditions: dict[str, dict[str, Any]],
    comparisons: list[tuple[str, str, str]],
) -> dict[str, str]:
    """Map comparison codes to display labels like 'HFD vs LFD'."""
    labels = {}
    for a, b, name in comparisons:
        da = conditions[a].get("display", a)
        db = conditions[b].get("display", b)
        labels[name] = f"{da} vs {db}"
    return labels


def get_cond_map(conditions: dict[str, dict[str, Any]]) -> dict[str, str]:
    """Map condition code -> display name."""
    return {code: info.get("display", code) for code, info in conditions.items()}


def get_plot_colors(cfg: MultiOmicsConfig) -> dict[str, Any]:
    """Merge user plot_colors over defaults."""
    merged = dict(_DEFAULT_PLOT_COLORS)
    if cfg.plot_colors:
        for key, val in cfg.plot_colors.items():
            if isinstance(val, dict) and isinstance(merged.get(key), dict):
                merged[key] = {**merged[key], **val}
            else:
                merged[key] = val
    return merged


# ---------------------------------------------------------------------------
# Main config loader
# ---------------------------------------------------------------------------

def load_multiomics_config(
    *,
    experiment_id: str | None = None,
    config_path: str | None = None,
    outdir_override: str | None = None,
) -> MultiOmicsConfig:
    """Load multi-omics integration config from YAML.

    Args:
        experiment_id: Experiment ID (looks for config in experiments/<id>/multiomics_integration_setting/)
        config_path: Direct path to config.yaml
        outdir_override: Override output directory

    Returns:
        Frozen MultiOmicsConfig instance.
    """
    # Resolve config file
    if config_path:
        cpath = Path(config_path).resolve()
    elif experiment_id:
        cpath = default_config_path(experiment_id)
    elif os.environ.get("EXPERIMENT_ID"):
        cpath = default_config_path(os.environ["EXPERIMENT_ID"])
    else:
        raise ValueError(
            "Must provide --experiment <id> or --config <path>. "
            "Config location: experiments/<id>/multiomics_integration_setting/config.yaml"
        )

    if not cpath.exists():
        raise FileNotFoundError(f"Config file not found: {cpath}")

    config_dir = cpath.parent
    data = _load_yaml(cpath)

    # ---- experiment_id ----
    eid = data.get("experiment_id", experiment_id or config_dir.parent.parent.name)

    # ---- Variable substitution prep ----
    output_dir_str = data.get("paths", {}).get("output_dir", "")
    vars_map = {
        "experiment_id": eid,
        "output_dir": output_dir_str,
    }

    # ---- Resolve paths ----
    paths = data.get("paths", {})

    def _rp(key: str, *, allow_create: bool = False) -> Path | None:
        raw = paths.get(key)
        if raw is None or raw == "null" or raw == "":
            return None
        raw = _substitute_vars(str(raw), vars_map)
        return _resolve_path(raw, config_dir, allow_create=allow_create)

    output_dir = _rp("output_dir", allow_create=True)
    if output_dir is None:
        raise ValueError(f"Missing paths.output_dir in {cpath}")

    if outdir_override:
        output_dir = Path(outdir_override).resolve()

    # Create output dirs
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "tables").mkdir(exist_ok=True)
    (output_dir / "plots").mkdir(exist_ok=True)

    # ---- Conditions and mouse groups (from groups.yaml or inline) ----
    groups_file = config_dir / "groups.yaml"
    if groups_file.exists():
        groups_data = _load_yaml(groups_file)
        conditions_raw = groups_data.get("conditions")
        mouse_groups_raw = groups_data.get("mouse_groups")
    else:
        # Fallback: read from config.yaml (backward compatible)
        conditions_raw = data.get("conditions")
        mouse_groups_raw = data.get("mouse_groups")

    if not conditions_raw or not isinstance(conditions_raw, dict):
        raise ValueError(f"Missing 'conditions' section (check groups.yaml or config.yaml in {config_dir})")

    conditions = {str(k): dict(v) for k, v in conditions_raw.items()}

    mouse_groups = {}
    if mouse_groups_raw and isinstance(mouse_groups_raw, dict):
        mouse_groups = {str(k): str(v) for k, v in mouse_groups_raw.items()}

    # ---- Params ----
    params = data.get("params", {})
    tissues = [_norm_tissue(t) for t in params.get("tissues", [])]
    comparisons = params.get("comparisons", [])
    mfa_tissue_map = params.get("mfa_tissue_map", {})
    comparison_map = params.get("comparison_map", {})
    proteomics_comparison_map = params.get("proteomics_comparison_map", None)
    metabolomics_comparison_map = params.get("metabolomics_comparison_map", None)
    # transcriptomics_comparison_map: explicit name preferred; falls back to comparison_map
    transcriptomics_comparison_map = params.get("transcriptomics_comparison_map", None)

    # Normalize mfa_tissue_map values
    mfa_tissue_map = {str(k): _norm_tissue(str(v)) for k, v in mfa_tissue_map.items()}

    # ---- Scoring ----
    scoring = dict(_DEFAULT_SCORING)
    user_scoring = params.get("scoring", {})
    for k, v in user_scoring.items():
        if isinstance(v, dict) and isinstance(scoring.get(k), dict):
            scoring[k] = {**scoring[k], **v}
        else:
            scoring[k] = v

    # ---- Display ----
    tissue_labels = params.get("tissue_labels")
    if tissue_labels:
        tissue_labels = {str(k): str(v) for k, v in tissue_labels.items()}

    plot_colors = params.get("plot_colors")

    # ---- Estimand / integration ----
    estimand = dict(data.get("estimand", {}) or {})
    integration = dict(data.get("integration", {}) or {})
    harmonization = dict(data.get("harmonization", {}) or {})

    # Resolve optional batch result paths once, so loaders do not need to know
    # where the project config lives. Both ``results`` and ``results_dir`` are accepted.
    inputs = integration.get("inputs", {}) or {}
    resolved_inputs: dict[str, Any] = {}
    for omics_name, input_spec in inputs.items():
        if not isinstance(input_spec, dict):
            resolved_inputs[str(omics_name)] = input_spec
            continue
        spec = dict(input_spec)
        batches = []
        for batch in spec.get("batches", []) or []:
            if not isinstance(batch, dict):
                raise ValueError(f"integration.inputs.{omics_name}.batches entries must be mappings")
            resolved_batch = dict(batch)
            for path_key in ("results", "results_dir"):
                raw_path = resolved_batch.get(path_key)
                if raw_path:
                    raw_path = _substitute_vars(str(raw_path), vars_map)
                    resolved_batch[path_key] = _resolve_path(raw_path, config_dir)
            batches.append(resolved_batch)
        if batches:
            spec["batches"] = batches
        resolved_inputs[str(omics_name)] = spec
    integration["inputs"] = resolved_inputs

    return MultiOmicsConfig(
        experiment_id=eid,
        config_path=cpath,
        output_dir=output_dir,
        mfa_results=_rp("mfa_results"),
        metabolomics_results=_rp("metabolomics_results"),
        proteomics_results=_rp("proteomics_results"),
        transcriptomics_results=_rp("transcriptomics_results"),
        conditions=conditions,
        mouse_groups=mouse_groups,
        tissues=tissues,
        comparisons=comparisons,
        mfa_tissue_map=mfa_tissue_map,
        comparison_map=comparison_map,
        transcriptomics_comparison_map=transcriptomics_comparison_map,
        proteomics_comparison_map=proteomics_comparison_map,
        metabolomics_comparison_map=metabolomics_comparison_map,
        estimand=estimand,
        integration=integration,
        harmonization=harmonization,
        scoring=scoring,
        tissue_labels=tissue_labels,
        plot_colors=plot_colors,
    )
