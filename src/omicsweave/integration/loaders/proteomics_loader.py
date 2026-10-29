#!/usr/bin/env python3
"""
Proteomics data loader with normalized column names.

Loads differential expression and pathway data from proteomics pipeline output,
normalizes column names to snake_case convention, and handles comparison name
reversal (e.g., Pval_LFD_vs_HSD vs Pval_HSD_vs_LFD) with automatic log2FC
sign flipping.

Reference: flux_protein_integration_analysis.py lines 679-768

Author: Claude Code
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pandas as pd

from omicsweave.integration.methods.batch_effect_meta import (
    meta_analyze_feature_batches,
)


# ---------------------------------------------------------------------------
# Tissue name mapping: lowercase -> directory casing on disk
# ---------------------------------------------------------------------------

TISSUE_DIR_MAP: dict[str, str] = {
    "brain": "Brain",
    "kidney": "Kidney",
    "liver": "Liver",
    "heart": "Heart",
    "muscle": "Gast_Muscle",
    "gast_muscle": "Gast_Muscle",
    "white_adipose": "SQ_Adipose",
    "sq_adipose": "SQ_Adipose",
}


# ---------------------------------------------------------------------------
# Column name normalization
# ---------------------------------------------------------------------------

# Patterns: (compiled regex, replacement)
_COL_PATTERNS: list[tuple[re.Pattern, str]] = [
    (re.compile(r"^Protein$", re.IGNORECASE), "protein"),
    (re.compile(r"^Pval_(.+)$", re.IGNORECASE), r"pval_\1"),
    (re.compile(r"^Log2FC_(.+)$", re.IGNORECASE), r"log2fc_\1"),
    (re.compile(r"^Padj_(.+)$", re.IGNORECASE), r"padj_\1"),
    (re.compile(r"^(?:SE|StdErr|LfcSE)_(.+)$", re.IGNORECASE), r"se_\1"),
    (re.compile(r"^Evidence_(.+)$", re.IGNORECASE), r"evidence_\1"),
    (re.compile(r"^InferenceStatus_(.+)$", re.IGNORECASE), r"inference_status_\1"),
    (re.compile(r"^NObsNum_(.+)$", re.IGNORECASE), r"n_obs_numerator_\1"),
    (re.compile(r"^NObsDen_(.+)$", re.IGNORECASE), r"n_obs_denominator_\1"),
    (re.compile(r"^ResidualDF_(.+)$", re.IGNORECASE), r"residual_df_\1"),
]


def _normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize proteomics DataFrame column names to snake_case convention.

    Rules:
        Protein          -> protein
        Pval_{comp}      -> pval_{comp}
        Log2FC_{comp}    -> log2fc_{comp}
        Padj_{comp}      -> padj_{comp}
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
# Comparison reversal helper
# ---------------------------------------------------------------------------

def _reverse_comparison(comp: str) -> str:
    """Reverse a comparison string like 'HFD_vs_LFD' -> 'LFD_vs_HFD'."""
    parts = comp.split("_vs_")
    if len(parts) == 2:
        return f"{parts[1]}_vs_{parts[0]}"
    return comp


# ---------------------------------------------------------------------------
# ProteomicsLoader
# ---------------------------------------------------------------------------

class ProteomicsLoader:
    """Load proteomics differential expression and pathway data.

    Parameters
    ----------
    config : object
        Configuration object with attributes:
        - proteomics_results : Path | None
            Root directory of proteomics output.
        - tissues : list[str]
            Normalized (lowercase) tissue names to load.
        - comparisons : list[str]
            Comparison names (e.g. ``["HFD_vs_LFD", "HSD_vs_LFD"]``).
        - comparison_map : dict[str, str]
            Maps short comparison names to condition-code names
            (e.g. ``{"HFD_vs_LFD": "hfd_14w_vs_lfd_14w"}``).
            Used to rename proteomics columns from condition-code to short names.
    """

    def __init__(self, config: Any) -> None:
        self.results_dir: Path | None = getattr(config, "proteomics_results", None)
        self.tissues: list[str] = list(getattr(config, "tissues", []))
        self.comparisons: list[str] = list(getattr(config, "comparisons", []))
        # Build reverse map: condition-code -> short name
        # Prefer proteomics_comparison_map if available; fall back to comparison_map.
        cmap: dict[str, str] = getattr(config, "proteomics_comparison_map", None)
        if cmap is None:
            cmap = getattr(config, "comparison_map", {}) or {}
        # Filter out null/None values (comparisons with no proteomics data)
        cmap = {k: v for k, v in cmap.items() if v is not None and v != "null"}
        self._reverse_comp_map: dict[str, str] = {v: k for k, v in cmap.items()}
        integration = getattr(config, "integration", {}) or {}
        inputs = integration.get("inputs", {}) if isinstance(integration, dict) else {}
        self._input_spec: dict[str, Any] = dict(inputs.get("proteomics", {}) or {})
        self._batch_specs: list[dict[str, Any]] = list(
            self._input_spec.get("batches", []) or []
        )
        default_mode = "integrated" if self._batch_specs else "single"
        self._input_mode = str(self._input_spec.get("mode", default_mode)).lower()
        self._allow_legacy_results = bool(
            self._input_spec.get("allow_legacy_results", False)
        )

    # ---- public API ----

    def load(self) -> dict[str, dict[str, pd.DataFrame | None]]:
        """Load all proteomics data.

        Returns
        -------
        dict
            ``{"diff": {tissue: DataFrame}, "pathway": {tissue: DataFrame | None}}``
            Only tissues with successfully loaded data are included.
        """
        if self._input_mode in {"integrated", "auto"} and self._batch_specs:
            return self._load_integrated_batches()
        if self._input_mode == "integrated" and not self._batch_specs:
            raise ValueError(
                "integration.inputs.proteomics.mode=integrated requires a non-empty batches list"
            )

        if self.results_dir is None:
            print("[WARN] ProteomicsLoader: proteomics_results is None, skipping.")
            return {"diff": {}, "pathway": {}}

        rd = Path(self.results_dir)
        if not rd.exists():
            print(f"[WARN] ProteomicsLoader: results dir not found: {rd}")
            return {"diff": {}, "pathway": {}}

        result = self._load_root(rd)
        print(f"[OK] Proteomics: loaded {len(result['diff'])} tissue(s)")
        return result

    def _load_root(self, rd: Path) -> dict[str, dict[str, pd.DataFrame | None]]:
        """读取一个已经独立完成 QC/差异分析的蛋白质组结果目录。"""
        diff: dict[str, pd.DataFrame] = {}
        pathway: dict[str, pd.DataFrame | None] = {}

        for tissue in self.tissues:
            tissue_dir = self._resolve_tissue_dir(rd, tissue)
            if tissue_dir is None:
                continue

            # -- differential expression --
            diff_file = tissue_dir / "differential_proteins.csv"
            if not diff_file.exists():
                diff_file = tissue_dir / "differential_expression.csv"
            if diff_file.exists():
                metadata_path = tissue_dir / "analysis_metadata.json"
                if not self._allow_legacy_results:
                    if not metadata_path.exists():
                        raise RuntimeError(
                            f"Proteomics results for '{tissue}' have no v6 analysis metadata; "
                            "re-run the corrected proteomics pipeline with --force. "
                            "Set integration.inputs.proteomics.allow_legacy_results=true "
                            "only for explicitly audited historical data."
                        )
                    analysis_metadata = json.loads(
                        metadata_path.read_text(encoding="utf-8")
                    )
                    if str(analysis_metadata.get("pipeline_version")) != "6.0":
                        raise RuntimeError(
                            f"Proteomics results for '{tissue}' were produced by pipeline "
                            f"{analysis_metadata.get('pipeline_version')!r}; expected '6.0'."
                        )
                df = pd.read_csv(diff_file)
                df = _normalize_columns(df)
                df = self._rename_condition_code_columns(df)
                df = self._handle_comparison_reversal(df)
                diff[tissue] = df
                print(f"[OK] Proteomics diff: {tissue} ({len(df)} proteins)")
            else:
                print(f"[WARN] ProteomicsLoader: no diff file for tissue '{tissue}' in {tissue_dir}")

            # -- pathway summary --
            pw_file = tissue_dir / "pathway_summary_all.csv"
            if pw_file.exists():
                pathway[tissue] = pd.read_csv(pw_file)
            else:
                pathway[tissue] = self._load_prerank_pathways(rd, tissue)

        return {"diff": diff, "pathway": pathway}

    def _load_prerank_pathways(
        self,
        root: Path,
        tissue: str,
    ) -> pd.DataFrame | None:
        """Load v6 comparison-level GSEA reports, excluding descriptive ssGSEA."""

        plot_root = root / "02_result_plot"
        if not plot_root.is_dir():
            return None
        tissue_plot_dir = self._resolve_tissue_dir(plot_root, tissue)
        if tissue_plot_dir is None:
            return None
        gsea_root = tissue_plot_dir / "gsea_prerank"
        if not gsea_root.is_dir():
            return None

        frames: list[pd.DataFrame] = []
        for report in sorted(gsea_root.glob("*/*/*prerank*.report.csv")):
            frame = pd.read_csv(report)
            frame["comparison"] = report.parents[1].name
            frame["library"] = report.parent.name
            frame["source_mode"] = "proteomics_prerank_gsea"
            frames.append(frame)
        return pd.concat(frames, ignore_index=True) if frames else None

    def _load_integrated_batches(self) -> dict[str, dict[str, pd.DataFrame | None]]:
        """在各批次上游分析完成后，按比较和蛋白执行结果级随机效应 Meta。"""
        payloads: dict[str, dict[str, dict[str, pd.DataFrame | None]]] = {}
        for index, spec in enumerate(self._batch_specs, start=1):
            batch_id = str(spec.get("batch_id") or f"batch_{index}")
            if batch_id in payloads:
                raise ValueError(f"Duplicate proteomics batch_id: {batch_id}")
            raw_path = spec.get("results_dir") or spec.get("results")
            if not raw_path:
                raise ValueError(f"Proteomics batch '{batch_id}' is missing results_dir/results")
            batch_root = Path(raw_path)
            if not batch_root.exists():
                raise FileNotFoundError(
                    f"Proteomics batch '{batch_id}' results not found: {batch_root}"
                )
            payloads[batch_id] = self._load_root(batch_root)

        integrated_diff: dict[str, pd.DataFrame] = {}
        integrated_pathway: dict[str, pd.DataFrame | None] = {}
        all_tissues = sorted({
            tissue
            for payload in payloads.values()
            for tissue in payload["diff"]
        })
        min_batches = int(self._input_spec.get("min_batches", 2))
        protocol_id = str(self._input_spec.get("protocol_id", ""))

        for tissue in all_tissues:
            batch_tables = {
                batch_id: payload["diff"][tissue]
                for batch_id, payload in payloads.items()
                if tissue in payload["diff"] and payload["diff"][tissue] is not None
            }
            if not batch_tables:
                continue
            combined = pd.concat(batch_tables.values(), ignore_index=True)
            statistic_prefixes = (
                "log2fc_", "pval_", "padj_", "se_", "n_batches_",
                "n_batches_observed_", "n_batches_eligible_", "batch_ids_",
                "eligible_batch_ids_", "heterogeneity_i2_", "heterogeneity_tau2_",
                "q_statistic_", "direction_consistent_", "meta_eligible_",
                "meta_note_", "se_source_", "source_mode_", "protocol_id_",
                "evidence_", "inference_status_", "n_obs_numerator_",
                "n_obs_denominator_", "residual_df_",
            )
            metadata_cols = [
                column for column in combined.columns
                if column == "protein" or not column.startswith(statistic_prefixes)
            ]
            integrated = combined[metadata_cols].dropna(subset=["protein"])
            integrated = integrated.drop_duplicates("protein")

            for comparison in self.comparisons:
                effect_column = f"log2fc_{comparison}"
                eligible = {
                    batch_id: table
                    for batch_id, table in batch_tables.items()
                    if effect_column in table.columns
                }
                meta = meta_analyze_feature_batches(
                    eligible,
                    feature_column="protein",
                    effect_column=effect_column,
                    se_candidates=[
                        f"se_{comparison}", f"lfcSE_{comparison}",
                        "lfcSE", "standard_error", "std_error", "se",
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
                    "pooled_fdr": f"padj_{comparison}",
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
                keep = ["protein", *rename]
                current = meta[keep].rename(columns=rename)
                meta_eligible_col = f"meta_eligible_{comparison}"
                meta_note_col = f"meta_note_{comparison}"
                current[f"evidence_{comparison}"] = (
                    current[meta_eligible_col]
                    .fillna(False)
                    .astype(bool)
                    .map({True: "inferential", False: "insufficient"})
                )
                current[f"inference_status_{comparison}"] = current[meta_note_col]
                current[f"source_mode_{comparison}"] = "within_omics_batch_meta"
                current[f"protocol_id_{comparison}"] = protocol_id
                integrated = integrated.merge(current, on="protein", how="outer")

            if len(integrated) > 0:
                integrated_diff[tissue] = integrated

            pathway_frames = []
            for batch_id, payload in payloads.items():
                table = payload["pathway"].get(tissue)
                if table is not None and not table.empty:
                    table = table.copy()
                    table["batch_id"] = batch_id
                    table["source_mode"] = "batch_specific_pathway"
                    pathway_frames.append(table)
            integrated_pathway[tissue] = (
                pd.concat(pathway_frames, ignore_index=True) if pathway_frames else None
            )

        print(
            f"[OK] Proteomics batch integration: {len(payloads)} batch(es), "
            f"{len(integrated_diff)} tissue(s)"
        )
        return {"diff": integrated_diff, "pathway": integrated_pathway}

    # ---- internal helpers ----

    def _resolve_tissue_dir(self, root: Path, tissue: str) -> Path | None:
        """Find the tissue sub-directory trying multiple casings.

        First resolves the scan directory: checks for ``01_dp_results/``
        subdirectory (current convention), then falls back to root (legacy).

        Then within the scan directory, tries multiple casings:
        1. TISSUE_DIR_MAP lookup (e.g. gast_muscle -> Gast_Muscle)
        2. lowercase (e.g. liver)
        3. Title_Case (e.g. Liver)
        4. original as-is
        5. Scan all subdirectories, match case-insensitively
        """
        # Resolve scan directory: prefer 01_dp_results/ subdirectory
        dp_dir = root / "01_dp_results"
        scan_root = dp_dir if dp_dir.is_dir() else root

        # 1. Known mapping
        mapped = TISSUE_DIR_MAP.get(tissue.lower(), None)
        candidates: list[str] = []
        if mapped:
            candidates.append(mapped)

        # 2. lowercase
        candidates.append(tissue.lower())

        # 3. Title Case (underscores kept, each segment title-cased)
        title = "_".join(s.capitalize() for s in tissue.lower().split("_"))
        candidates.append(title)

        # 4. Original
        candidates.append(tissue)

        for cand in candidates:
            p = scan_root / cand
            if p.is_dir():
                return p

        # 5. Fallback: scan subdirectories, case-insensitive match
        norm = tissue.lower().replace(" ", "_").replace("-", "_")
        try:
            for sub in scan_root.iterdir():
                if sub.is_dir() and sub.name.lower().replace(" ", "_").replace("-", "_") == norm:
                    return sub
        except OSError:
            pass

        print(f"[WARN] ProteomicsLoader: tissue dir not found for '{tissue}' in {scan_root}")
        return None

    def _rename_condition_code_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """Rename condition-code comparison columns to short names.

        Proteomics data may use condition-code comparison names like
        ``pval_hfd_14w_vs_lfd_14w`` instead of the short ``pval_HFD_vs_LFD``.
        This method uses the reverse of ``comparison_map`` to rename them.
        """
        if not self._reverse_comp_map:
            return df

        rename_map: dict[str, str] = {}
        prefixes = (
            "pval_", "log2fc_", "padj_", "se_", "evidence_",
            "inference_status_", "n_obs_numerator_", "n_obs_denominator_",
            "residual_df_",
        )
        for cond_code, short_name in self._reverse_comp_map.items():
            for prefix in prefixes:
                old_col = f"{prefix}{cond_code}"
                new_col = f"{prefix}{short_name}"
                if old_col in df.columns and new_col not in df.columns:
                    rename_map[old_col] = new_col

        if rename_map:
            df = df.rename(columns=rename_map)
            print(f"    [rename] {list(rename_map.values())}")

        return df

    def _handle_comparison_reversal(self, df: pd.DataFrame) -> pd.DataFrame:
        """Ensure comparison columns match the requested direction.

        For each configured comparison (e.g. ``HSD_vs_LFD``), if the column
        ``pval_HSD_vs_LFD`` does not exist but ``pval_LFD_vs_HSD`` does, the
        column is renamed and the ``log2fc_*`` sign is flipped.  This mirrors
        the logic from flux_protein_integration_analysis.py lines 753-768.
        """
        for comp in self.comparisons:
            pval_col = f"pval_{comp}"
            log2fc_col = f"log2fc_{comp}"
            padj_col = f"padj_{comp}"
            evidence_col = f"evidence_{comp}"
            status_col = f"inference_status_{comp}"

            if pval_col in df.columns:
                # Already present in the correct direction
                continue

            rev = _reverse_comparison(comp)
            rev_pval = f"pval_{rev}"
            rev_log2fc = f"log2fc_{rev}"
            rev_padj = f"padj_{rev}"
            rev_se = f"se_{rev}"
            rev_evidence = f"evidence_{rev}"
            rev_status = f"inference_status_{rev}"

            if rev_pval in df.columns:
                # Rename pval
                df = df.rename(columns={rev_pval: pval_col})
                # Flip log2fc sign
                if rev_log2fc in df.columns:
                    df = df.rename(columns={rev_log2fc: log2fc_col})
                    df[log2fc_col] = -df[log2fc_col]
                # Rename padj (no sign flip needed)
                if rev_padj in df.columns:
                    df = df.rename(columns={rev_padj: padj_col})
                if rev_se in df.columns:
                    df = df.rename(columns={rev_se: f"se_{comp}"})
                if rev_evidence in df.columns:
                    df = df.rename(columns={rev_evidence: evidence_col})
                if rev_status in df.columns:
                    df = df.rename(columns={rev_status: status_col})

        return df
