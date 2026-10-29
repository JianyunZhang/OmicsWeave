#!/usr/bin/env python3
"""
Cross-omics pathway convergence analysis.

For each tissue and comparison, identifies pathway terms shared across multiple omics layers
(metabolomics, proteomics, transcriptomics). Reports pairwise overlaps and
generates convergence heatmaps.

Extracted from multiomics_integration_analysis.py lines 441-528.

Author: Claude Code
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import beta, combine_pvalues


class PathwayConvergenceAnalyzer:
    """Analyze pathway convergence across omics layers.

    Parameters
    ----------
    output_dir : pathlib.Path
        Directory to write output tables and plots.
    pval_cutoff : float
        P-value threshold for pathway significance (default 0.05).
    gsea_fdr_cutoff : float
        FDR threshold for GSEA pathway significance (default 0.25,
        following standard GSEA practice).
    """

    def __init__(
        self,
        output_dir: Any,
        pval_cutoff: float = 0.05,
        gsea_fdr_cutoff: float = 0.25,
    ) -> None:
        self.output_dir = output_dir
        self.pval_cutoff = pval_cutoff
        self.gsea_fdr_cutoff = gsea_fdr_cutoff
        self.evidence_long_df = pd.DataFrame()
        self.fusion_df = pd.DataFrame()

    @staticmethod
    def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
        """Benjamini-Hochberg FDR correction.

        Parameters
        ----------
        pvals : array-like
            Raw p-values.

        Returns
        -------
        np.ndarray
            Adjusted FDR values (same shape).
        """
        pvals = np.asarray(pvals, dtype=float)
        if len(pvals) == 0:
            return np.array([], dtype=float)
        adjusted = np.full(pvals.shape, np.nan, dtype=float)
        finite_mask = np.isfinite(pvals)
        finite = pvals[finite_mask]
        n = len(finite)
        if n == 0:
            return adjusted
        order = np.argsort(finite)
        ranked = finite[order]
        adjusted_ranked = ranked * n / np.arange(1, n + 1)
        adjusted_ranked = np.minimum.accumulate(adjusted_ranked[::-1])[::-1]
        finite_adjusted = np.empty(n, dtype=float)
        finite_adjusted[order] = np.clip(adjusted_ranked, 0.0, 1.0)
        adjusted[finite_mask] = finite_adjusted
        return adjusted

    @staticmethod
    def _normalize_term(value: Any) -> str:
        term = str(value).split("__")[-1].strip().lower()
        return " ".join(term.split())

    @staticmethod
    def _first_probability_column(
        df: pd.DataFrame,
        candidates: tuple[str, ...],
    ) -> str | None:
        for column in candidates:
            if column in df.columns and pd.to_numeric(df[column], errors="coerce").notna().any():
                return column
        return None

    def _evidence_rows(
        self,
        df: pd.DataFrame,
        *,
        tissue: str,
        comparison: str,
        omics: str,
        term_column: str,
        probability_columns: tuple[str, ...],
    ) -> list[dict[str, Any]]:
        if df is None or df.empty or term_column not in df.columns:
            return []
        probability_column = self._first_probability_column(df, probability_columns)
        working = df[[term_column]].copy()
        working["pathway_term"] = working[term_column].map(self._normalize_term)
        working = working[working["pathway_term"] != ""]
        if working.empty:
            return []

        if probability_column:
            probabilities = pd.to_numeric(df.loc[working.index, probability_column], errors="coerce")
            working["probability"] = probabilities.clip(0.0, 1.0)
            working = working.sort_values("probability", na_position="last")
        else:
            working["probability"] = np.nan
        working["rank"] = np.arange(1, len(working) + 1)
        working["rank_fraction"] = working["rank"] / (len(working) + 1.0)

        return [
            {
                "tissue": tissue,
                "comparison": comparison,
                "omics": omics,
                "pathway_term": row["pathway_term"],
                "probability": row["probability"],
                "rank": int(row["rank"]),
                "rank_fraction": float(row["rank_fraction"]),
                "source_probability_column": probability_column or "rank_only",
            }
            for _, row in working.iterrows()
        ]

    def _collect_pathway_evidence(
        self,
        tissues: list[str],
        metab_pathway: dict[str, dict[str, pd.DataFrame]],
        prot_pathway: dict[str, pd.DataFrame | None],
        tx_gsea: dict[str, dict[str, pd.DataFrame]],
    ) -> pd.DataFrame:
        rows: list[dict[str, Any]] = []
        for tissue in sorted(set(tissues)):
            for comparison, df in metab_pathway.get(tissue, {}).items():
                rows.extend(self._evidence_rows(
                    df,
                    tissue=tissue,
                    comparison=str(comparison),
                    omics="metabolomics",
                    term_column="Name",
                    probability_columns=("FDR", "fdr", "P-value", "pvalue", "pval"),
                ))

            prot_df = prot_pathway.get(tissue)
            if prot_df is not None:
                rows.extend(self._evidence_rows(
                    prot_df,
                    tissue=tissue,
                    comparison="pooled",
                    omics="proteomics",
                    term_column="Top Pathway",
                    probability_columns=("FDR", "FDR q-val", "fdr", "pvalue", "pval"),
                ))

            for comparison, df in tx_gsea.get(tissue, {}).items():
                rows.extend(self._evidence_rows(
                    df,
                    tissue=tissue,
                    comparison=str(comparison),
                    omics="transcriptomics",
                    term_column="Term",
                    probability_columns=("FDR q-val", "fdr", "NOM p-val", "pvalue", "pval"),
                ))
        return pd.DataFrame(rows)

    @staticmethod
    def _rra_probability(rank_fractions: np.ndarray) -> float:
        """Robust Rank Aggregation 的 beta 顺序统计量近似。"""
        values = np.sort(np.asarray(rank_fractions, dtype=float))
        values = values[np.isfinite(values)]
        n_values = len(values)
        if n_values < 2:
            return float("nan")
        probabilities = [
            beta.cdf(value, rank, n_values - rank + 1)
            for rank, value in enumerate(values, start=1)
        ]
        return float(min(1.0, min(probabilities) * n_values))

    def _fuse_pathway_evidence(self, evidence_df: pd.DataFrame) -> pd.DataFrame:
        if evidence_df.empty:
            return pd.DataFrame(columns=[
                "tissue", "comparison", "pathway_term", "n_omics", "contributing_omics",
                "fusion_method", "fused_pvalue", "fused_fdr",
            ])

        # comparison 是 estimand 身份的一部分，禁止跨比较池化。
        per_omics_rows: list[dict[str, Any]] = []
        for (tissue, comparison, omics, term), group in evidence_df.groupby(
            ["tissue", "comparison", "omics", "pathway_term"], sort=True
        ):
            probabilities = pd.to_numeric(group["probability"], errors="coerce").dropna()
            if probabilities.empty:
                combined_probability = float("nan")
            else:
                combined_probability = float(min(1.0, probabilities.min() * len(probabilities)))
            per_omics_rows.append({
                "tissue": tissue,
                "comparison": comparison,
                "omics": omics,
                "pathway_term": term,
                "probability": combined_probability,
                "rank_fraction": float(group["rank_fraction"].min()),
            })
        per_omics = pd.DataFrame(per_omics_rows)

        fused_rows: list[dict[str, Any]] = []
        for (tissue, comparison, term), group in per_omics.groupby(
            ["tissue", "comparison", "pathway_term"], sort=True
        ):
            if group["omics"].nunique() < 2:
                continue
            probabilities = pd.to_numeric(group["probability"], errors="coerce").dropna().to_numpy()
            if len(probabilities) >= 2:
                _, fused_probability = combine_pvalues(
                    np.clip(probabilities, np.finfo(float).tiny, 1.0),
                    method="fisher",
                )
                fusion_method = "fisher"
            else:
                fused_probability = self._rra_probability(group["rank_fraction"].to_numpy())
                fusion_method = "rra"

            pvalue_map = {
                row["omics"]: None if pd.isna(row["probability"]) else round(float(row["probability"]), 6)
                for _, row in group.iterrows()
            }
            rank_map = {
                row["omics"]: round(float(row["rank_fraction"]), 6)
                for _, row in group.iterrows()
            }
            fused_rows.append({
                "tissue": tissue,
                "comparison": comparison,
                "pathway_term": term,
                "n_omics": int(group["omics"].nunique()),
                "contributing_omics": ", ".join(sorted(group["omics"].unique())),
                "fusion_method": fusion_method,
                "fused_pvalue": float(fused_probability),
                "omics_probabilities": str(pvalue_map),
                "omics_rank_fractions": str(rank_map),
            })

        fused = pd.DataFrame(fused_rows)
        if fused.empty:
            return pd.DataFrame(columns=[
                "tissue", "comparison", "pathway_term", "n_omics", "contributing_omics",
                "fusion_method", "fused_pvalue", "fused_fdr",
            ])
        fused["fused_fdr"] = np.nan
        for _, indices in fused.groupby(["tissue", "comparison"]).groups.items():
            fused.loc[indices, "fused_fdr"] = self._bh_fdr(
                fused.loc[indices, "fused_pvalue"].to_numpy(dtype=float)
            )
        return fused.sort_values(
            ["fused_fdr", "fused_pvalue", "n_omics"],
            ascending=[True, True, False],
        ).reset_index(drop=True)

    def analyze(
        self,
        tissues: list[str],
        metab_pathway: dict[str, dict[str, pd.DataFrame]],
        prot_pathway: dict[str, pd.DataFrame | None],
        tx_gsea: dict[str, dict[str, pd.DataFrame]],
    ) -> pd.DataFrame:
        """Run pathway convergence analysis across all tissues.

        Parameters
        ----------
        tissues : list[str]
            Tissue names to analyze.
        metab_pathway : dict
            ``{tissue: {comparison: DataFrame}}`` — metabolomics pathway results
            with columns ``Name`` and ``P-value``.
        prot_pathway : dict
            ``{tissue: DataFrame | None}`` — proteomics pathway summary with
            column ``Top Pathway``.
        tx_gsea : dict
            ``{tissue: {comparison: DataFrame}}`` — transcriptomics GSEA results
            with columns ``Term`` and ``FDR q-val``.

        Returns
        -------
        pd.DataFrame
            Columns: tissue, comparison, omics_A, omics_B, pathway_term.
        """
        all_rows: list[dict[str, str]] = []
        tables_dir = self.output_dir / "tables"
        tables_dir.mkdir(parents=True, exist_ok=True)

        self.evidence_long_df = self._collect_pathway_evidence(
            tissues=tissues,
            metab_pathway=metab_pathway,
            prot_pathway=prot_pathway,
            tx_gsea=tx_gsea,
        )
        self.fusion_df = self._fuse_pathway_evidence(self.evidence_long_df)
        self.evidence_long_df.to_csv(
            tables_dir / "pathway_evidence_long.csv", index=False, encoding="utf-8-sig"
        )
        self.fusion_df.to_csv(
            tables_dir / "pathway_fusion.csv", index=False, encoding="utf-8-sig"
        )
        if not self.fusion_df.empty:
            n_fused = int((self.fusion_df["fused_fdr"] < self.pval_cutoff).sum())
            print(
                f"[OK] Pathway evidence fusion: {len(self.fusion_df)} candidates "
                f"({n_fused} FDR<{self.pval_cutoff:g})"
            )

        for tissue in sorted(set(tissues)):
            # comparison 属于 estimand 身份；任何精确名称交集也不得先跨比较合并。
            comparison_sets: dict[str, dict[str, set[str]]] = {}

            # Metabolomics pathways — use BH-corrected FDR when available
            if tissue in metab_pathway:
                for comp, df in metab_pathway[tissue].items():
                    terms: set[str] = set()
                    if "Name" in df.columns and "P-value" in df.columns:
                        try:
                            pvals = df["P-value"].astype(float).values
                            # Apply BH FDR correction
                            fdr_vals = self._bh_fdr(pvals)
                            sig = df[fdr_vals < self.pval_cutoff]
                            terms = {
                                self._normalize_term(term)
                                for term in sig["Name"].dropna()
                            }
                        except (ValueError, TypeError):
                            pass
                    if terms:
                        comparison_sets.setdefault(str(comp), {})[
                            "metabolomics"
                        ] = terms

            # Proteomics pathway summary — filter by FDR when available
            if tissue in prot_pathway and prot_pathway[tissue] is not None:
                df = prot_pathway[tissue]
                if "Top Pathway" in df.columns:
                    comparison_column = next(
                        (column for column in ("comparison", "Comparison") if column in df.columns),
                        None,
                    )
                    groups = (
                        df.groupby(comparison_column, dropna=False)
                        if comparison_column
                        else [("pooled", df)]
                    )
                    for comp, comp_df in groups:
                        probability_column = self._first_probability_column(
                            comp_df, ("FDR", "FDR q-val", "fdr", "pvalue", "pval")
                        )
                        if probability_column is None:
                            # 仅有排序时交给 RRA；不把“存在条目”冒充显著交集。
                            continue
                        try:
                            probabilities = pd.to_numeric(
                                comp_df[probability_column], errors="coerce"
                            )
                            sig = comp_df[probabilities < self.pval_cutoff]
                            terms = {
                                self._normalize_term(term)
                                for term in sig["Top Pathway"].dropna()
                            }
                        except (ValueError, TypeError):
                            terms = set()
                        if terms:
                            comparison_sets.setdefault(str(comp), {})[
                                "proteomics"
                            ] = terms

            # Transcriptomics GSEA
            if tissue in tx_gsea:
                for comp, gdf in tx_gsea[tissue].items():
                    terms: set[str] = set()
                    if "Term" in gdf.columns and "FDR q-val" in gdf.columns:
                        sig = gdf[gdf["FDR q-val"] < self.gsea_fdr_cutoff]
                        raw = sig["Term"].dropna()
                        # Term format: "Library__Pathway Name"
                        terms = {self._normalize_term(term) for term in raw}
                    if terms:
                        comparison_sets.setdefault(str(comp), {})[
                            "transcriptomics"
                        ] = terms

            for comparison, pw_sets in sorted(comparison_sets.items()):
                if len(pw_sets) < 2:
                    continue

                # Pairwise overlaps are computed only inside one estimand context.
                omics_list = sorted(pw_sets.keys())
                for i in range(len(omics_list)):
                    for j in range(i + 1, len(omics_list)):
                        a, b = omics_list[i], omics_list[j]
                        overlap = pw_sets[a] & pw_sets[b]
                        if overlap:
                            for term in sorted(overlap):
                                all_rows.append({
                                    "tissue": tissue,
                                    "comparison": comparison,
                                    "omics_A": a,
                                    "omics_B": b,
                                    "pathway_term": term,
                                })

        if all_rows:
            conv_df = pd.DataFrame(all_rows)

            # Save detailed results
            conv_df.to_csv(tables_dir / "pathway_convergence.csv", index=False, encoding="utf-8-sig")

            # Summary counts
            counts = (
                conv_df
                .groupby(["tissue", "comparison", "omics_A", "omics_B"])
                .size()
                .reset_index(name="n_shared")
            )
            counts.to_csv(tables_dir / "pathway_convergence_counts.csv", index=False, encoding="utf-8-sig")

            print(
                f"[OK] Pathway convergence: {len(all_rows)} shared terms "
                f"across {conv_df['tissue'].nunique()} tissues"
            )
            return conv_df
        else:
            print("[INFO] No cross-omics pathway convergence found")
            return pd.DataFrame(
                columns=["tissue", "comparison", "omics_A", "omics_B", "pathway_term"]
            )
