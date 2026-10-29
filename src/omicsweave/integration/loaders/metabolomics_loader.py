#!/usr/bin/env python3
"""
Metabolomics data loader with normalized column names.

Loads differential analysis and pathway data from metabolomics pipeline output,
normalizes column names to snake_case convention.

Reference:
  - flux_metabolite_cross_analysis.py lines 979-1032
  - multiomics_integration_analysis.py lines 104-138

Author: Claude Code
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pandas as pd

from omicsweave.integration.methods.batch_effect_meta import (
    meta_analyze_feature_batches,
)


# ---------------------------------------------------------------------------
# Column name normalization
# ---------------------------------------------------------------------------

# Patterns: (compiled regex, replacement)
_COL_PATTERNS: list[tuple[re.Pattern, str]] = [
    (re.compile(r"^Compounds$", re.IGNORECASE), "compound"),
    (re.compile(r"^Index$", re.IGNORECASE), "feature_id"),
    (re.compile(r"^Pval_(.+)$", re.IGNORECASE), r"pval_\1"),
    (re.compile(r"^Log2FC_(.+)$", re.IGNORECASE), r"log2fc_\1"),
    (re.compile(r"^FDR_(.+)$", re.IGNORECASE), r"fdr_\1"),
    (re.compile(r"^(?:SE|StdErr)_(.+)$", re.IGNORECASE), r"se_\1"),
    # Standard database ID columns (from metabolomics pipeline output)
    (re.compile(r"^HMDB$", re.IGNORECASE), "hmdb_id"),
    (re.compile(r"^KEGG_Compound$", re.IGNORECASE), "kegg_compound"),
    (re.compile(r"^ChEBI$", re.IGNORECASE), "chebi_id"),
    # These are kept as-is but listed for clarity
    # Formula, VIP_Global, kegg_map
]

# Columns that should pass through unchanged
_KEEP_COLS = {"formula", "vip_global", "kegg_map"}


def _normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize metabolomics DataFrame column names to snake_case convention.

    Rules:
        Compounds        -> compound
        Index            -> feature_id
        HMDB             -> hmdb_id
        KEGG_Compound    -> kegg_compound
        ChEBI            -> chebi_id
        Formula          -> Formula  (kept as-is)
        VIP_Global       -> VIP_Global  (kept as-is)
        kegg_map         -> kegg_map  (kept as-is)
        Pval_{comp}      -> pval_{comp}
        Log2FC_{comp}    -> log2fc_{comp}
        FDR_{comp}       -> fdr_{comp}
    """
    new_cols: list[str] = []
    for col in df.columns:
        renamed = col
        for pat, repl in _COL_PATTERNS:
            m = pat.match(col)
            if m:
                renamed = pat.sub(repl, col)
                break
        new_cols.append(renamed)
    df.columns = new_cols
    return df


# ---------------------------------------------------------------------------
# Tissue directory resolution
# ---------------------------------------------------------------------------

def _norm_tissue(name: str) -> str:
    """Normalize tissue name to lowercase_with_underscores."""
    return name.strip().lower().replace(" ", "_").replace("-", "_")


def _resolve_tissue_dir(base: Path, tissue: str) -> Path | None:
    """Find the tissue sub-directory within metabolomics results.

    Tries the 00_dm_results subdirectory first (current convention),
    then falls back to the root directory (legacy convention).

    Within the resolved scan directory, attempts to match the tissue name
    case-insensitively.
    """
    # Current layout: base/00_dm_results/{tissue}/
    dm_dir = base / "00_dm_results"
    scan_dir = dm_dir if dm_dir.exists() else base

    if not scan_dir.exists():
        return None

    norm = _norm_tissue(tissue)

    # Direct match (lowercase)
    candidate = scan_dir / norm
    if candidate.is_dir():
        return candidate

    # Scan all subdirs case-insensitively
    try:
        for sub in scan_dir.iterdir():
            if sub.is_dir() and _norm_tissue(sub.name) == norm:
                return sub
    except OSError:
        pass

    return None


# ---------------------------------------------------------------------------
# MetabolomicsLoader
# ---------------------------------------------------------------------------

class MetabolomicsLoader:
    """Load metabolomics differential analysis and pathway data.

    Parameters
    ----------
    config : object
        Configuration object with attributes:
        - metabolomics_results : Path | None
            Root directory of metabolomics output.
        - tissues : list[str]
            Normalized (lowercase) tissue names to load.
        - comparisons : list[str]
            Comparison names (e.g. ``["HFD_vs_LFD", "HSD_vs_LFD"]``).
        - metabolomics_comparison_map : dict[str, str] | None
            Maps MFA comparison keys to metabolomics column naming
            (e.g. ``{"HFD_vs_LFD": "HFD-14W_vs_LFD-14W"}``).
            Used to rename metabolomics columns from pipeline format to short names.
    """

    def __init__(self, config: Any) -> None:
        self.results_dir: Path | None = getattr(config, "metabolomics_results", None)
        self.tissues: list[str] = list(getattr(config, "tissues", []))
        self.comparisons: list[str] = list(getattr(config, "comparisons", []))
        self._mfa_comparisons: list[str] = list(self.comparisons)
        # Build reverse map: metabolomics column name → short MFA key
        cmap: dict[str, str] | None = getattr(config, "metabolomics_comparison_map", None)
        if cmap is None:
            cmap = {}
        # Filter out null/None values
        cmap = {k: v for k, v in cmap.items() if v is not None and v != "null"}
        self._reverse_comp_map: dict[str, str] = {v: k for k, v in cmap.items()}
        integration = getattr(config, "integration", {}) or {}
        inputs = integration.get("inputs", {}) if isinstance(integration, dict) else {}
        self._input_spec: dict[str, Any] = dict(inputs.get("metabolomics", {}) or {})
        self._batch_specs: list[dict[str, Any]] = list(
            self._input_spec.get("batches", []) or []
        )
        default_mode = "integrated" if self._batch_specs else "single"
        self._input_mode = str(self._input_spec.get("mode", default_mode)).lower()

    # ---- public API ----

    def load(self) -> dict[str, dict]:
        """Load all metabolomics data.

        Returns
        -------
        dict
            ``{
                "diff": {tissue: DataFrame},
                "pathway": {tissue: {comparison: DataFrame}}
            }``
            Only tissues with successfully loaded data are included.
        """
        if self._input_mode in {"integrated", "auto"} and self._batch_specs:
            return self._load_integrated_batches()
        if self._input_mode == "integrated" and not self._batch_specs:
            raise ValueError(
                "integration.inputs.metabolomics.mode=integrated requires a non-empty batches list"
            )

        if self.results_dir is None:
            print("[WARN] MetabolomicsLoader: metabolomics_results is None, skipping.")
            return {"diff": {}, "pathway": {}}

        rd = Path(self.results_dir)
        if not rd.exists():
            print(f"[WARN] MetabolomicsLoader: results dir not found: {rd}")
            return {"diff": {}, "pathway": {}}

        result = self._load_root(rd)
        print(f"[OK] Metabolomics: loaded {len(result['diff'])} tissue(s)")
        return result

    def _load_root(self, rd: Path) -> dict[str, dict]:
        """读取一个已经独立完成 QC/差异分析的代谢组结果目录。"""
        diff: dict[str, pd.DataFrame] = {}
        pathway: dict[str, dict[str, pd.DataFrame]] = {}

        for tissue in self.tissues:
            tissue_dir = _resolve_tissue_dir(rd, tissue)
            if tissue_dir is None:
                print(f"[WARN] MetabolomicsLoader: tissue dir not found for '{tissue}' in {rd}")
                continue

            # -- differential analysis --
            # Prefer per-comparison files (dm_{tissue}_{comp}.csv), fall back to combined
            dm_files = sorted(tissue_dir.glob(f"dm_{tissue}_*_vs_*.csv"))
            if dm_files:
                # Load per-comparison files and merge into combined DataFrame
                comp_dfs: list[pd.DataFrame] = []
                for dm_file in dm_files:
                    cdf = pd.read_csv(dm_file)
                    # Extract comp_name from filename: dm_{tissue}_{comp_name}.csv
                    prefix = f"dm_{tissue}_"
                    comp_name = dm_file.stem[len(prefix):]
                    cdf = cdf.rename(columns={
                        "log2FoldChange": f"Log2FC_{comp_name}",
                        "statistic": f"Tstat_{comp_name}",
                        "pvalue": f"Pval_{comp_name}",
                        "padj": f"FDR_{comp_name}",
                        "nominal_significant": f"NominalSignificant_{comp_name}",
                        "significant_univariate": f"SignificantUnivariate_{comp_name}",
                        "vip_pass": f"VIPPass_{comp_name}",
                        "significant_with_vip": f"SignificantWithVIP_{comp_name}",
                        "significant": f"Significant_{comp_name}",
                        "evidence_level": f"Evidence_{comp_name}",
                    })
                    comp_dfs.append(cdf)

                # Merge on shared ID columns
                from functools import reduce
                merge_cols = ["Index", "Compounds"]
                optional_cols = ["Formula", "kegg_map", "VIP_Global", "HMDB", "KEGG_Compound", "ChEBI"]
                for opt in optional_cols:
                    if all(opt in cdf.columns for cdf in comp_dfs):
                        merge_cols.append(opt)
                # Keep only merge cols + comparison-specific cols
                kept_dfs = []
                for cdf in comp_dfs:
                    extra = [c for c in cdf.columns if c not in merge_cols]
                    kept_dfs.append(cdf[merge_cols + extra])
                df = reduce(lambda l, r: pd.merge(l, r, on=merge_cols, how="outer"), kept_dfs)

                df = _normalize_columns(df)
                df = self._rename_comparison_columns(df)
                diff[tissue] = df
                print(f"[OK] Metabolomics diff (per-comparison): {tissue} ({len(df)} compounds)")
            else:
                # Fallback: legacy combined file
                diff_file = tissue_dir / "differential_analysis.csv"
                if diff_file.exists():
                    df = pd.read_csv(diff_file)
                    df = _normalize_columns(df)
                    df = self._rename_comparison_columns(df)
                    diff[tissue] = df
                    print(f"[OK] Metabolomics diff (combined): {tissue} ({len(df)} compounds)")
                else:
                    print(f"[WARN] MetabolomicsLoader: no diff file for '{tissue}' at {tissue_dir}")

            # -- pathway files (Pathway_{comparison}.csv) --
            pw: dict[str, pd.DataFrame] = {}
            for pf in sorted(tissue_dir.glob("Pathway_*.csv")):
                comp = pf.stem.replace("Pathway_", "")
                canonical_comp = self._reverse_comp_map.get(comp, comp)
                pw[canonical_comp] = pd.read_csv(pf)
            if pw:
                pathway[tissue] = pw

        return {"diff": diff, "pathway": pathway}

    def _load_integrated_batches(self) -> dict[str, dict]:
        """在各批次上游分析完成后，按比较和代谢物执行结果级随机效应 Meta。"""
        payloads: dict[str, dict[str, dict]] = {}
        for index, spec in enumerate(self._batch_specs, start=1):
            batch_id = str(spec.get("batch_id") or f"batch_{index}")
            if batch_id in payloads:
                raise ValueError(f"Duplicate metabolomics batch_id: {batch_id}")
            raw_path = spec.get("results_dir") or spec.get("results")
            if not raw_path:
                raise ValueError(f"Metabolomics batch '{batch_id}' is missing results_dir/results")
            batch_root = Path(raw_path)
            if not batch_root.exists():
                raise FileNotFoundError(
                    f"Metabolomics batch '{batch_id}' results not found: {batch_root}"
                )
            payloads[batch_id] = self._load_root(batch_root)

        integrated_diff: dict[str, pd.DataFrame] = {}
        integrated_pathway: dict[str, dict[str, pd.DataFrame]] = {}
        all_tissues = sorted({
            tissue
            for payload in payloads.values()
            for tissue in payload["diff"]
        })
        min_batches = int(self._input_spec.get("min_batches", 2))
        protocol_id = str(self._input_spec.get("protocol_id", ""))
        requested_feature = str(self._input_spec.get("feature_column", "compound"))

        for tissue in all_tissues:
            batch_tables = {
                batch_id: payload["diff"][tissue]
                for batch_id, payload in payloads.items()
                if tissue in payload["diff"]
            }
            if not batch_tables:
                continue
            combined = pd.concat(batch_tables.values(), ignore_index=True)
            feature_column = requested_feature
            if feature_column not in combined.columns:
                feature_column = "feature_id" if "feature_id" in combined.columns else "compound"
            if feature_column not in combined.columns:
                raise ValueError(
                    f"Metabolomics batch tables for tissue '{tissue}' lack feature column "
                    f"'{requested_feature}', 'feature_id', and 'compound'"
                )

            statistic_prefixes = (
                "log2fc_", "pval_", "fdr_", "se_", "n_batches_",
                "n_batches_observed_", "n_batches_eligible_", "batch_ids_",
                "eligible_batch_ids_", "heterogeneity_i2_", "heterogeneity_tau2_",
                "q_statistic_", "direction_consistent_", "meta_eligible_",
                "meta_note_", "se_source_", "source_mode_", "protocol_id_",
            )
            metadata_cols = [
                column for column in combined.columns
                if column == feature_column or not column.startswith(statistic_prefixes)
            ]
            integrated = combined[metadata_cols].dropna(subset=[feature_column])
            integrated = integrated.drop_duplicates(feature_column)

            for comparison in self.comparisons:
                effect_column = f"log2fc_{comparison}"
                eligible = {
                    batch_id: table
                    for batch_id, table in batch_tables.items()
                    if effect_column in table.columns and feature_column in table.columns
                }
                meta = meta_analyze_feature_batches(
                    eligible,
                    feature_column=feature_column,
                    effect_column=effect_column,
                    se_candidates=[
                        f"se_{comparison}", "standard_error", "std_error", "se",
                    ],
                    pvalue_candidates=[f"pval_{comparison}", "pvalue", "pval"],
                    min_batches=min_batches,
                )
                if meta.empty:
                    continue
                rename = {
                    "pooled_effect": effect_column,
                    "pooled_se": f"se_{comparison}",
                    "pooled_pvalue": f"pval_{comparison}",
                    "pooled_fdr": f"fdr_{comparison}",
                    "n_batches": f"n_batches_{comparison}",
                    "n_batches_observed": f"n_batches_observed_{comparison}",
                    "n_batches_eligible": f"n_batches_eligible_{comparison}",
                    "batch_ids": f"batch_ids_{comparison}",
                    "eligible_batch_ids": f"eligible_batch_ids_{comparison}",
                    "i_squared": f"heterogeneity_i2_{comparison}",
                    "tau_squared": f"heterogeneity_tau2_{comparison}",
                    "q_statistic": f"q_statistic_{comparison}",
                    "direction_consistent": f"direction_consistent_{comparison}",
                    "meta_eligible": f"meta_eligible_{comparison}",
                    "meta_note": f"meta_note_{comparison}",
                    "se_source": f"se_source_{comparison}",
                }
                keep = [feature_column, *rename]
                current = meta[keep].rename(columns=rename)
                current[f"source_mode_{comparison}"] = "within_omics_batch_meta"
                current[f"protocol_id_{comparison}"] = protocol_id
                integrated = integrated.merge(current, on=feature_column, how="outer")

            if len(integrated) > 0:
                integrated_diff[tissue] = integrated

            tissue_pathways: dict[str, list[pd.DataFrame]] = {}
            for batch_id, payload in payloads.items():
                for comparison, table in payload["pathway"].get(tissue, {}).items():
                    copied = table.copy()
                    copied["batch_id"] = batch_id
                    copied["source_mode"] = "batch_specific_pathway"
                    tissue_pathways.setdefault(comparison, []).append(copied)
            if tissue_pathways:
                integrated_pathway[tissue] = {
                    comparison: pd.concat(frames, ignore_index=True)
                    for comparison, frames in tissue_pathways.items()
                }

        print(
            f"[OK] Metabolomics batch integration: {len(payloads)} batch(es), "
            f"{len(integrated_diff)} tissue(s)"
        )
        return {"diff": integrated_diff, "pathway": integrated_pathway}

    # ---- internal helpers ----

    @staticmethod
    def _extract_conditions(comp_name: str) -> tuple[str, str] | None:
        """Extract the two condition codes from a comparison name like 'A_vs_B'.

        Returns (cond_a, cond_b) or None if not parseable.
        """
        if "_vs_" not in comp_name:
            return None
        parts = comp_name.split("_vs_", 1)
        return parts[0], parts[1]

    @staticmethod
    def _condition_fingerprint(cond: str) -> str:
        """Generate a normalized fingerprint for fuzzy condition matching.

        Strips common age/time suffixes (e.g., '-14W', '-16W', '-8W') to
        enable matching between names like 'HFD-14W-SG-16W' and 'HFD_SG_8W'.
        """
        import re as _re
        # Normalize separators
        s = cond.replace("-", "_")
        # Strip numeric-age suffixes like _14W, _16W, _8W, _4W
        s = _re.sub(r"_\d+W", "", s)
        # Also strip standalone age patterns at end
        s = _re.sub(r"\d+[WwDd]$", "", s)
        return s.lower().strip("_")

    def _auto_map_comparisons(self, file_comp_names: list[str]) -> dict[str, str]:
        """Auto-discover mapping from metabolomics file comp names to MFA keys.

        Uses fuzzy matching on condition codes: strips age/week suffixes and
        compares the core condition identifiers.

        Parameters
        ----------
        file_comp_names : list[str]
            Comparison names extracted from dm_* filenames
            (e.g., ``["HFD-14W-SG-16W_vs_HFD-14W-SHAM-16W"]``).

        Returns
        -------
        dict[str, str]
            ``{file_comp_name: mfa_short_name}`` for matched comparisons.
        """
        result: dict[str, str] = {}

        for fcn in file_comp_names:
            fcond = self._extract_conditions(fcn)
            if fcond is None:
                continue
            fp_a = self._condition_fingerprint(fcond[0])
            fp_b = self._condition_fingerprint(fcond[1])

            for mfa_key in self._mfa_comparisons:
                mcond = self._extract_conditions(mfa_key)
                if mcond is None:
                    continue
                mp_a = self._condition_fingerprint(mcond[0])
                mp_b = self._condition_fingerprint(mcond[1])

                # Check both orderings: (A→A, B→B) or (A→B, B→A)
                if (fp_a == mp_a and fp_b == mp_b) or (fp_a == mp_b and fp_b == mp_a):
                    result[fcn] = mfa_key
                    break

        if result:
            auto_names = [f"{k} -> {v}" for k, v in result.items()]
            print(f"    [auto-map] {auto_names}")

        return result

    def _rename_comparison_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """Rename metabolomics comparison columns to MFA short names.

        Uses explicit ``metabolomics_comparison_map`` first, then falls back
        to auto-discovered mappings from file names.
        """
        # Collect all comparison suffixes from column names
        file_comp_names: set[str] = set()
        prefixes = ("pval_", "log2fc_", "fdr_", "se_")
        for col in df.columns:
            for prefix in prefixes:
                if col.startswith(prefix):
                    file_comp_names.add(col[len(prefix):])

        # Build rename map: explicit map (from config) takes precedence
        rename_map: dict[str, str] = {}

        # Auto-discovered mappings (lower priority)
        auto_map = self._auto_map_comparisons(sorted(file_comp_names))

        for metab_name, short_name in auto_map.items():
            for prefix in prefixes:
                old_col = f"{prefix}{metab_name}"
                new_col = f"{prefix}{short_name}"
                if old_col in df.columns and new_col not in df.columns:
                    rename_map[old_col] = new_col

        # Explicit config map overrides auto-discovered
        for metab_name, short_name in self._reverse_comp_map.items():
            for prefix in prefixes:
                old_col = f"{prefix}{metab_name}"
                new_col = f"{prefix}{short_name}"
                if old_col in df.columns and new_col not in df.columns:
                    rename_map[old_col] = new_col

        if rename_map:
            df = df.rename(columns=rename_map)
            print(f"    [rename] {list(rename_map.values())}")

        return df
