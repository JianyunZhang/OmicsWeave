#!/usr/bin/env python3
"""
MFA (Metabolic Flux Analysis) data loader with normalized column names.

Loads differential flux, pathway flux, TCA cycle, and overview data from
MFA results directories. Applies tissue name normalization and renames
columns to canonical lowercase forms for consistent downstream use.

Usage:
    from omicsweave.integration.loaders import MFALoader
    loader = MFALoader(config)
    data = loader.load()
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd


_META_COL_MAP: dict[str, str] = {
    "log2FC_pooled": "flux_log2fc",
    "se_pooled": "flux_se",
    "cohen_d_pooled": "cohen_d",
    "hedges_g_pooled": "hedges_g",
    "pooled_pvalue": "flux_pvalue",
    "pooled_fdr": "flux_fdr",
    "I_squared": "heterogeneity_i2",
    "tau_squared": "heterogeneity_tau2",
}


class MFALoader:
    """Load MFA results with normalized column names.

    Accepts a config-like object (e.g. ``MultiOmicsConfig``) that provides:
        - ``mfa_results`` (Path | None): path to the MFA results directory
        - ``mfa_tissue_map`` (dict): maps raw MFA tissue names to canonical names

    Example::

        from omicsweave.integration.config import load_multiomics_config
        cfg = load_multiomics_config(experiment_id="20250425_cimr_wt_mice_dio")
        loader = MFALoader(cfg)
        data = loader.load()
        # data["diff"]["HFD_vs_LFD"]  -> DataFrame with normalized columns
        # data["tissues"]              -> {"liver", "heart", ...}
    """

    # Column renaming map: original -> canonical
    _DIFF_COL_MAP: dict[str, str] = {
        "flux_id": "flux_id",
        "tissue": "tissue",
        "log2FC": "flux_log2fc",
        "pvalue_mannwhitney": "flux_pvalue",
        "fdr_mannwhitney": "flux_fdr",
        "cohen_d": "cohen_d",
        "reaction": "reaction",
        # 通路注释列（从 MFA 输出直接透传，无需重命名）
        "ec_number": "ec_number",
        "kegg_reaction": "kegg_reaction",
        "kegg_pathway": "kegg_pathway",
        "reactome_id": "reactome_id",
        "description": "reaction_description",
    }

    def __init__(self, config: Any) -> None:
        self._config = config
        self._results_dir: Path | None = None
        self._tissue_map: dict[str, str] = {}
        self._comparisons: set[str] = set()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def load(self) -> dict[str, Any]:
        """Load all MFA data and return a structured dict.

        Returns:
            {
                "diff":    {comparison_name: DataFrame, ...},
                "pathway": DataFrame | None,
                "tca":     DataFrame | None,
                "overview": DataFrame | None,
                "tissues": set[str],
            }
        """
        # Resolve config attributes (duck-typing)
        mfa_results = getattr(self._config, "mfa_results", None)
        self._tissue_map = getattr(self._config, "mfa_tissue_map", {}) or {}
        self._comparisons = {
            str(value) for value in (getattr(self._config, "comparisons", []) or [])
        }

        if mfa_results is None:
            print("[WARN] MFA results path not configured (mfa_results is None).")
            return self._empty_result()

        self._results_dir = Path(mfa_results)
        if not self._results_dir.exists():
            print(f"[WARN] MFA results directory not found: {self._results_dir}")
            return self._empty_result()
        integrated_dir = self._resolve_integrated_dir(self._results_dir)
        self._results_dir = integrated_dir or self._resolve_results_dir(self._results_dir)

        print(f"[INFO] Loading MFA results from: {self._results_dir}")

        diff_tables = self._load_differential_flux()
        pathway = self._load_pathway()
        tca = self._load_tca()
        overview = self._load_overview()

        # Collect tissues from all loaded differential tables
        tissues: set[str] = set()
        for df in diff_tables.values():
            if "tissue" in df.columns:
                tissues |= set(df["tissue"].dropna().unique())

        n_diff = len(diff_tables)
        n_tissues = len(tissues)
        print(
            f"[OK] MFA results loaded: {n_diff} differential tables, "
            f"{n_tissues} tissues, "
            f"pathway={'Yes' if pathway is not None else 'No'}, "
            f"TCA={'Yes' if tca is not None else 'No'}, "
            f"overview={'Yes' if overview is not None else 'No'}"
        )

        return {
            "diff": diff_tables,
            "pathway": pathway,
            "tca": tca,
            "overview": overview,
            "tissues": tissues,
        }

    # ------------------------------------------------------------------
    # Internal loaders
    # ------------------------------------------------------------------

    def _load_differential_flux(self) -> dict[str, pd.DataFrame]:
        """Load all differential flux CSVs (recursive glob across tissue subdirs) and normalize columns."""
        if self._is_integrated_dir(self._results_dir):
            return self._load_integrated_differential_flux()

        diff_root = self._results_dir / "01_df_results"
        if not diff_root.exists():
            print(f"[WARN] 01_df_results directory not found: {diff_root}")
            return {}

        tables: dict[str, pd.DataFrame] = {}
        for fp in sorted(diff_root.glob("**/differential_flux_*.csv")):
            comp = fp.stem.replace("differential_flux_", "")
            if self._comparisons and comp not in self._comparisons:
                continue
            try:
                df = pd.read_csv(fp)
            except Exception as exc:
                print(f"[WARN] Failed to read {fp.name}: {exc}")
                continue

            df = self._normalize_tissue_col(df)
            df = self._normalize_diff_columns(df)
            # 按比较名合并不同组织的 CSV
            if comp in tables:
                tables[comp] = pd.concat([tables[comp], df], ignore_index=True)
            else:
                tables[comp] = df

        for comp, df in tables.items():
            print(f"  [diff] {comp}: {len(df)} flux entries")

        return tables

    def _load_integrated_differential_flux(self) -> dict[str, pd.DataFrame]:
        """读取同一示踪方案内的批次 Meta 结果，并转换为评分器标准列。"""
        tables: dict[str, pd.DataFrame] = {}
        for fp in sorted(self._results_dir.glob("differential_flux_meta_*.csv")):
            comp = fp.stem.replace("differential_flux_meta_", "")
            if self._comparisons and comp not in self._comparisons:
                continue
            try:
                df = pd.read_csv(fp)
            except Exception as exc:
                print(f"[WARN] Failed to read integrated MFA result {fp.name}: {exc}")
                continue
            df = df.rename(columns={k: v for k, v in _META_COL_MAP.items() if k in df.columns})
            df = self._normalize_tissue_col(df)
            df = self._normalize_diff_columns(df)
            if "flux_pvalue" not in df.columns and "is_significant_pooled" in df.columns:
                # 仅兼容旧产物；新产物必须提供真实 pooled_pvalue/pooled_fdr。
                significant = df["is_significant_pooled"].map(
                    lambda value: str(value).strip().lower() in {"true", "1", "yes"}
                )
                df["flux_pvalue"] = significant.map({True: 0.049, False: 1.0})
                df["flux_fdr"] = df["flux_pvalue"]
                df["pvalue_source"] = "legacy_significance_flag_approximation"
            else:
                df["pvalue_source"] = "pooled_random_effects_meta"
            df["source_mode"] = "within_protocol_meta"
            if "meta_note" in df.columns:
                df["evidence_mode"] = df["meta_note"].map(
                    lambda value: (
                        "biological_replicate"
                        if str(value) in {"random_effects_meta", "high_heterogeneity"}
                        else "descriptive"
                    )
                )
            tables[comp] = df
            print(f"  [meta diff] {comp}: {len(df)} flux entries")
        return tables

    def _load_pathway(self) -> pd.DataFrame | None:
        """Load pathway flux summary (recursive glob across tissue subdirs)."""
        if self._is_integrated_dir(self._results_dir):
            return None
        df_root = self._results_dir / "01_df_results"
        if not df_root.exists():
            return None

        all_frames = []
        for fp in sorted(df_root.glob("**/pathway_flux/pathway_flux_samples.csv")):
            try:
                df = pd.read_csv(fp)
                all_frames.append(df)
            except Exception as exc:
                print(f"[WARN] Failed to read pathway flux: {exc}")
                continue

        if not all_frames:
            print(f"[WARN] No pathway flux files found")
            return None

        df = pd.concat(all_frames, ignore_index=True)
        df = self._normalize_tissue_col(df)
        print(f"  [pathway] {len(df)} entries")
        return df

    def _load_tca(self) -> pd.DataFrame | None:
        """Load TCA cycle fluxes; derive from pathway stats when standalone tca_cycle files are absent."""
        if self._is_integrated_dir(self._results_dir):
            return None
        df_root = self._results_dir / "01_df_results"
        if not df_root.exists():
            return None

        all_frames = []
        # 旧/历史路径：独立 tca_cycle CSV
        for fp in sorted(df_root.glob("**/tca_cycle/tca_cycle_fluxes.csv")):
            try:
                df = pd.read_csv(fp)
                all_frames.append(df)
            except Exception as exc:
                print(f"[WARN] Failed to read TCA cycle: {exc}")
                continue

        # 新路径：不再单独输出 tca_cycle，改从 pathway_flux_reaction_stats 派生
        if not all_frames:
            seen_tissues: set[str] = set()
            for fp in sorted(df_root.glob("*/pathway_flux/pathway_flux_reaction_stats.csv")):
                tissue = fp.parent.parent.name
                if tissue == "all_tissues":
                    continue
                try:
                    df = pd.read_csv(fp)
                except Exception as exc:
                    print(f"[WARN] Failed to read pathway stats for TCA fallback: {exc}")
                    continue
                if "pathway" not in df.columns:
                    continue
                df = df[df["pathway"].astype(str) == "TCA_Cycle"].copy()
                if df.empty:
                    continue
                seen_tissues.add(str(tissue))
                all_frames.append(df)

            cross_fp = df_root / "all_tissues" / "pathway_flux" / "pathway_flux_reaction_stats.csv"
            if cross_fp.exists():
                try:
                    cross_df = pd.read_csv(cross_fp)
                    if "pathway" in cross_df.columns:
                        cross_df = cross_df[cross_df["pathway"].astype(str) == "TCA_Cycle"].copy()
                    if not cross_df.empty:
                        if "tissue" not in cross_df.columns or not seen_tissues:
                            all_frames.append(cross_df)
                        else:
                            for tissue in sorted(cross_df["tissue"].dropna().astype(str).unique()):
                                if tissue in seen_tissues:
                                    continue
                                tissue_df = cross_df[cross_df["tissue"].astype(str) == tissue]
                                if not tissue_df.empty:
                                    all_frames.append(tissue_df.copy())
                except Exception as exc:
                    print(f"[WARN] Failed to read cross-tissue pathway stats for TCA fallback: {exc}")

        if not all_frames:
            print(f"[WARN] No TCA data found (neither tca_cycle files nor pathway TCA_Cycle stats)")
            return None

        df = pd.concat(all_frames, ignore_index=True)
        keep_cols = [c for c in ["tissue", "reaction", "condition", "mean", "std", "n"] if c in df.columns]
        if keep_cols:
            df = df[keep_cols]
        df = self._normalize_tissue_col(df)
        print(f"  [tca] {len(df)} entries")
        return df

    def _load_overview(self) -> pd.DataFrame | None:
        """Load tissue flux overview."""
        if self._is_integrated_dir(self._results_dir):
            return None
        candidates = [
            self._results_dir / "01_df_results" / "all_tissues" / "overview" / "tissue_flux_overview.csv",
            self._results_dir / "01_df_results" / "overview" / "tissue_flux_overview.csv",
        ]
        fp = next((p for p in candidates if p.exists()), None)
        if fp is None:
            print(f"[WARN] Overview file not found in: {candidates[0].parent} or {candidates[1].parent}")
            return None
        try:
            df = pd.read_csv(fp)
        except Exception as exc:
            print(f"[WARN] Failed to read overview: {exc}")
            return None
        df = self._normalize_tissue_col(df)
        print(f"  [overview] {len(df)} entries")
        return df

    # ------------------------------------------------------------------
    # Normalization helpers
    # ------------------------------------------------------------------

    def _normalize_tissue_col(self, df: pd.DataFrame) -> pd.DataFrame:
        """Map tissue values using mfa_tissue_map."""
        if "tissue" in df.columns and self._tissue_map:
            df["tissue"] = df["tissue"].map(
                lambda t: self._tissue_map.get(t, t) if pd.notna(t) else t
            )
        return df

    def _normalize_diff_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """Rename differential flux columns to canonical names.

        Rules:
            - flux_id   -> strip tissue prefix (e.g. "heart--HEX_c" -> "HEX_c")
            - tissue    -> tissue    (values normalized via mfa_tissue_map)
            - log2FC    -> flux_log2fc
            - pvalue_mannwhitney -> flux_pvalue  (prefer mannwhitney)
            - pvalue_ttest       -> flux_pvalue  (fallback)
            - fdr_mannwhitney    -> flux_fdr     (prefer mannwhitney)
            - fdr_ttest          -> flux_fdr     (fallback)
            - cohen_d   -> cohen_d   (kept)
            - reaction  -> reaction  (kept)
        """
        # Strip tissue prefix from flux_id (e.g. "heart--HEX_c" -> "HEX_c")
        if "flux_id" in df.columns:
            df["flux_id"] = (
                df["flux_id"]
                .astype(str)
                .str.replace(r"^[^-]+--", "", regex=True)
            )

        # Build the p-value column: prefer mannwhitney, fallback to ttest
        if "pvalue_mannwhitney" in df.columns:
            df["flux_pvalue"] = df["pvalue_mannwhitney"]
        elif "pvalue_ttest" in df.columns:
            df["flux_pvalue"] = df["pvalue_ttest"]

        # Build the FDR column: prefer mannwhitney, fallback to ttest
        if "fdr_mannwhitney" in df.columns:
            df["flux_fdr"] = df["fdr_mannwhitney"]
        elif "fdr_ttest" in df.columns:
            df["flux_fdr"] = df["fdr_ttest"]

        # Rename log2FC
        if "log2FC" in df.columns:
            df = df.rename(columns={"log2FC": "flux_log2fc"})

        return df

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _empty_result() -> dict[str, Any]:
        """Return an empty result dict when MFA data is not available."""
        return {
            "diff": {},
            "pathway": None,
            "tca": None,
            "overview": None,
            "tissues": set(),
        }

    def _resolve_results_dir(self, root_dir: Path) -> Path:
        """
        解析 MFA 结果目录：
        - 若 root_dir 直接包含 01_df_results，则按单次运行目录处理；
        - 否则尝试从 root_dir/runs/<protocol>/<batch>/ 里选择一个有效运行目录。
        """
        if (root_dir / "01_df_results").exists():
            return root_dir

        runs_root = root_dir / "runs"
        if not runs_root.exists():
            return root_dir

        candidates = [
            p for p in sorted(runs_root.glob("*/*"))
            if p.is_dir() and (p / "01_df_results").exists()
        ]
        if not candidates:
            return root_dir

        pref_protocol = (
            getattr(self._config, "mfa_run_protocol", None)
            or getattr(self._config, "mfa_protocol_id", None)
            or ""
        )
        pref_batch = (
            getattr(self._config, "mfa_run_batch", None)
            or getattr(self._config, "mfa_batch_id", None)
            or ""
        )
        pref_protocol = str(pref_protocol).strip()
        pref_batch = str(pref_batch).strip()
        if pref_protocol and pref_batch:
            preferred = runs_root / pref_protocol / pref_batch
            if preferred in candidates:
                print(f"[INFO] Using configured MFA run: {pref_protocol}/{pref_batch}")
                return preferred

        selected = max(candidates, key=lambda p: p.stat().st_mtime)
        rel = selected.relative_to(root_dir)
        print(f"[INFO] Auto-selected MFA run directory: {rel}")
        return selected

    @staticmethod
    def _is_integrated_dir(path: Path | None) -> bool:
        return bool(path and any(path.glob("differential_flux_meta_*.csv")))

    def _resolve_integrated_dir(self, root_dir: Path) -> Path | None:
        """按配置选择批次 Meta 目录；显式 integrated 模式下缺失即报错。"""
        integration = getattr(self._config, "integration", {}) or {}
        mfa_input = ((integration.get("inputs", {}) or {}).get("mfa", {}) or {})
        mode = str(mfa_input.get("mode", "single_run")).lower()
        if mode not in {"integrated", "auto"}:
            return None

        protocol_id = str(mfa_input.get("protocol_id", "")).strip()
        base = root_dir / "integration" / "tracer_harmonized"
        candidates = []
        if protocol_id:
            candidates.append(base / protocol_id)
        candidates.extend(sorted(p for p in base.glob("*") if p.is_dir()) if base.exists() else [])
        selected = next((p for p in candidates if self._is_integrated_dir(p)), None)
        if selected is not None:
            print(f"[INFO] Using integrated MFA batch results: {selected}")
            return selected
        if mode == "integrated":
            raise FileNotFoundError(
                f"Configured integrated MFA results were not found under: {base}"
            )
        return None
