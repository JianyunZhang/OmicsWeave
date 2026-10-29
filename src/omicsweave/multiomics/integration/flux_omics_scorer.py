#!/usr/bin/env python3
"""
Flux-centric multi-omics evidence profiling with a legacy MFA-gated score.

The historical 0-4 score still uses MFA as both an entry gate and a scoring
layer. It is emitted as ``legacy_*`` fields only. The primary evidence profile
continues evaluating other omics when the MFA gate fails.

The primary output separates support, conflict, coverage, stability, design
grade, and evidence grade. A long-form EvidenceTable retains unavailable,
unmapped, unmatched, no-signal, and opposite-direction states.

Scoring breakdown:
  - MFA (gate + layer): 0–1.0  (significance 0/0.3 + effect tier 0–0.3 + bio annotation 0/0.4)
  - Transcriptomics:    0–1.0  (significance 0/0.4 + direction 0/0.3 + effect tier 0–0.3)
  - Proteomics:         0–1.2  (significance 0/0.5 + direction 0/0.3 + effect tier 0–0.4)
  - Metabolomics:       0–0.8  (significance 0/0.3 + substrate/product-aware direction 0/0.3 + effect tier 0–0.2)
  - Total:              0–4.0

Confidence: High(>=3.0) / Medium(>=2.0) / Low(>=1.0) / Insufficient(<1.0)

The legacy MFA gate uses FDR<0.05 AND |Cohen's d|>0.8 for reproducibility.
This does not repair pseudo-replication; MCMC-fit-point inference is marked as
``pseudo_replicate`` and restricted to exploratory design grade D.

Author: Claude Code
"""

from __future__ import annotations

from collections import Counter
from typing import Any

import numpy as np
import pandas as pd

from ..mappings import (
    FLUX_GENE_MAPPING,
    REACTION_ENZYME_MAPPING,
    FLUX_METABOLITE_MAPPING,
    get_metabolite_roles,
    match_gene_features,
    match_metabolite_roles,
    match_protein_features,
    audit_mapping_integrity,
)
from ..mappings.reaction_registry import get_kegg_pathway
from .evidence_profile import EvidenceProfileScorer, evidence_table


# ---------------------------------------------------------------------------
# Scoring helper functions
# ---------------------------------------------------------------------------

def _effect_tier(log2fc: float, tiers: list[float] | None = None) -> float:
    """Map |log2FC| to a tiered weight for better effect-size discrimination.

    Parameters
    ----------
    log2fc : float
        Log2 fold change value.
    tiers : list of two floats
        [moderate_lower, large_lower] thresholds.
        |fc| < moderate_lower -> 0.0 (trivial)
        moderate_lower <= |fc| < large_lower -> 0.5 (moderate)
        |fc| >= large_lower -> 1.0 (large)
        Default: [0.58, 1.0]
    """
    if tiers is None:
        tiers = [0.58, 1.0]
    if pd.isna(log2fc):
        return 0.0
    abs_fc = abs(log2fc)
    if abs_fc >= tiers[1]:
        return 1.0
    if abs_fc >= tiers[0]:
        return 0.5
    return 0.0


def _confidence_label(score: float, thresholds: dict[str, float]) -> str:
    """Map a numeric score to a confidence label."""
    if score >= thresholds.get("high", 3.0):
        return "High"
    if score >= thresholds.get("medium", 2.0):
        return "Medium"
    if score >= thresholds.get("low", 1.0):
        return "Low"
    return "Insufficient"


def _simes_probability(values: pd.Series) -> float:
    """Simes 聚合，减少“映射对象越多越容易取到极小 p 值”的偏倚。"""
    numeric = pd.to_numeric(values, errors="coerce").dropna().clip(0.0, 1.0)
    if numeric.empty:
        return float("nan")
    ordered = np.sort(numeric.to_numpy(dtype=float))
    ranks = np.arange(1, len(ordered) + 1, dtype=float)
    return float(min(1.0, np.min(ordered * len(ordered) / ranks)))


def _first_numeric_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for column in candidates:
        if column in df.columns and pd.to_numeric(df[column], errors="coerce").notna().any():
            return column
    return None


def _series_numeric_value(row: pd.Series, candidates: list[str]) -> float:
    for key in candidates:
        value = pd.to_numeric(pd.Series([row.get(key, np.nan)]), errors="coerce").iloc[0]
        if pd.notna(value) and np.isfinite(float(value)):
            return float(value)
    return float("nan")


def _mapping_method_counts(methods: dict[Any, str]) -> dict[str, int]:
    """将逐行匹配来源压缩为可序列化审计计数。"""
    return {
        str(method): int(count)
        for method, count in Counter(methods.values()).items()
    }


def _row_evidence_summary(
    matched_df: pd.DataFrame,
    signal_df: pd.DataFrame,
    *,
    effect_col: str | None,
    qvalue_candidates: list[str],
    pvalue_candidates: list[str],
    se_candidates: list[str] | None = None,
    n_batches_candidates: list[str] | None = None,
    batch_ids_candidates: list[str] | None = None,
    i2_candidates: list[str] | None = None,
    tau2_candidates: list[str] | None = None,
    protocol_candidates: list[str] | None = None,
    source_mode_candidates: list[str] | None = None,
) -> dict[str, Any]:
    """从映射后的特征表提取 EvidenceTable 所需的连续统计量。"""
    summary = {
        "representative_effect": float("nan"),
        "se": float("nan"),
        "qvalue": float("nan"),
        "lfsr": float("nan"),
        "n_eff": float("nan"),
        "selection_probability": float("nan"),
        "direction_stability": float("nan"),
    }
    if matched_df.empty:
        return summary

    effect_source = signal_df if not signal_df.empty else matched_df
    if effect_col and effect_col in effect_source.columns:
        effects = pd.to_numeric(effect_source[effect_col], errors="coerce").dropna()
        if not effects.empty:
            summary["representative_effect"] = float(effects.median())

    qvalue_col = _first_numeric_column(matched_df, qvalue_candidates)
    pvalue_col = _first_numeric_column(matched_df, pvalue_candidates)
    probability_col = qvalue_col or pvalue_col
    if probability_col:
        summary["qvalue"] = _simes_probability(matched_df[probability_col])

    lfsr_col = _first_numeric_column(effect_source, ["lfsr", "local_false_sign_rate"])
    if lfsr_col:
        lfsr_values = pd.to_numeric(effect_source[lfsr_col], errors="coerce").dropna()
        if not lfsr_values.empty:
            summary["lfsr"] = float(lfsr_values.median())

    se_col = _first_numeric_column(
        effect_source,
        (se_candidates or [])
        + ["lfcSE", "log2fc_se", "log2FC_se", "standard_error", "std_error", "se"],
    )
    if se_col:
        se_values = pd.to_numeric(effect_source[se_col], errors="coerce").dropna().abs()
        if not se_values.empty:
            summary["se"] = float(se_values.median())

    n_eff_col = _first_numeric_column(
        matched_df,
        ["n_eff", "effective_n", "sample_size", "n_samples", "n"],
    )
    if n_eff_col:
        n_values = pd.to_numeric(matched_df[n_eff_col], errors="coerce").dropna()
        if not n_values.empty:
            summary["n_eff"] = float(n_values.median())
    else:
        case_col = _first_numeric_column(matched_df, ["n_case", "n_treatment"])
        control_col = _first_numeric_column(matched_df, ["n_control", "n_reference"])
        if case_col and control_col:
            totals = (
                pd.to_numeric(matched_df[case_col], errors="coerce")
                + pd.to_numeric(matched_df[control_col], errors="coerce")
            ).dropna()
            if not totals.empty:
                summary["n_eff"] = float(totals.median())

    for output_key, candidates in (
        ("selection_probability", ["selection_probability", "bootstrap_selection_probability"]),
        ("direction_stability", ["direction_stability", "bootstrap_direction_stability"]),
    ):
        column = _first_numeric_column(matched_df, candidates)
        if column:
            values = pd.to_numeric(matched_df[column], errors="coerce").dropna()
            if not values.empty:
                summary[output_key] = float(values.median())
    for output_key, candidates in (
        ("n_batches", (n_batches_candidates or []) + ["n_batches"]),
        ("heterogeneity_i2", (i2_candidates or []) + ["heterogeneity_i2", "I_squared", "i2"]),
        ("heterogeneity_tau2", (tau2_candidates or []) + ["heterogeneity_tau2", "tau_squared"]),
    ):
        column = _first_numeric_column(matched_df, candidates)
        if column:
            values = pd.to_numeric(matched_df[column], errors="coerce").dropna()
            if not values.empty:
                summary[output_key] = float(values.median())
    for output_key, candidates in (
        ("batch_ids", (batch_ids_candidates or []) + ["batch_ids"]),
        ("protocol_id", (protocol_candidates or []) + ["protocol_id"]),
        ("source_mode", (source_mode_candidates or []) + ["source_mode"]),
    ):
        column = next((name for name in candidates if name in matched_df.columns), None)
        if column:
            values = matched_df[column].dropna().astype(str)
            if not values.empty:
                summary[output_key] = values.iloc[0]
    return summary


# ---------------------------------------------------------------------------
# FluxOmicsScorer
# ---------------------------------------------------------------------------

class FluxOmicsScorer:
    """Score each MFA differential flux against multi-omics evidence.

    Parameters
    ----------
    config : object
        Configuration object (``MultiOmicsConfig``) providing:
        - scoring : dict with keys pval_cutoff, log2fc_cutoff, confidence_thresholds, weights
    """

    def __init__(self, config: Any) -> None:
        self.mapping_audit = audit_mapping_integrity()
        self.mapping_audit.raise_for_errors()
        self.experiment_id = str(getattr(config, "experiment_id", ""))
        estimand = getattr(config, "estimand", {}) or {}
        self.estimand_id = str(estimand.get("id", ""))
        scoring = getattr(config, "scoring", {})
        self.pval_cutoff: float = scoring.get("pval_cutoff", 0.05)
        self.log2fc_cutoff: float = scoring.get("log2fc_cutoff", 0.58)

        # Per-omics cutoffs: each omics can override global pval_cutoff / log2fc_cutoff
        global_defaults = {"pval_cutoff": self.pval_cutoff, "log2fc_cutoff": self.log2fc_cutoff}
        omics_cutoffs = scoring.get("omics_cutoffs", {})
        self.omics_cutoffs: dict[str, dict[str, float]] = {}
        for key in ("mfa", "transcriptomics", "proteomics", "metabolomics"):
            self.omics_cutoffs[key] = {**global_defaults, **omics_cutoffs.get(key, {})}

        # Effect tier thresholds from config
        effect_tiers = scoring.get("effect_tiers", {})
        self.log2fc_tiers: list[float] = effect_tiers.get("log2fc", [0.58, 1.0])
        self.cohen_d_tiers: list[float] = effect_tiers.get("cohen_d", [0.8, 1.5])

        self.confidence_thresholds: dict[str, float] = scoring.get(
            "confidence_thresholds", {"high": 3.0, "medium": 2.0, "low": 1.0}
        )
        self.weights: dict[str, float] = scoring.get("weights", {
            "transcriptomics": 1.0,
            "proteomics": 1.2,
            "metabolomics": 0.8,
            "mfa": 1.0,
        })

        # Fallback mode for N=1 / low-n experiments (effect-size-only scoring)
        fallback = scoring.get("fallback_mode", {})
        self.fallback_enabled: bool = fallback.get("enabled", True)
        self.fallback_mfa_gate: dict[str, float] = fallback.get("mfa_gate", {
            "cohen_d_minimum": self.cohen_d_tiers[0],   # 0.8
            "log2fc_minimum": self.log2fc_tiers[0],      # 0.58
        })
        self.fallback_fc_large: float = fallback.get("effect_size_only_thresholds", {}).get(
            "log2fc_large", self.log2fc_tiers[1],        # 1.0
        )
        self.fallback_weight_penalty: float = fallback.get("weight_penalty", 0.6)

        # 新证据画像与旧矩阵解耦。旧分数继续保留，仅用于结果复现。
        self.evidence_profile_scorer = EvidenceProfileScorer(
            scoring.get("evidence_profile", {})
        )
        self._evidence_table = evidence_table([])

    # ---- per-omics scoring methods ----

    def score_transcriptomics(
        self,
        flux_row: pd.Series,
        tx_df: pd.DataFrame | None,
    ) -> tuple[float, dict[str, Any]]:
        """Score a single flux reaction against transcriptomics DE data.

        Only concordant-direction evidence contributes to the score.

        Two evidence modes:
          - statistical_test: padj < cutoff AND |log2fc| > cutoff (standard)
          - effect_size_only: no p-value; |log2fc| >= fallback_fc_large (N=1 fallback)

        Score components (statistical_test, when direction concordance exists):
          - significance: 0.4  |  effect_size_only: 0.2
          - direction:    0.3  |  effect_size_only: 0.2
          - effect:       0–0.3 (same tier)
          - weight penalty: 1.0 | effect_size_only: fallback_weight_penalty
        """
        if tx_df is None or tx_df.empty:
            return 0.0, {"n_genes_checked": 0, "n_sig": 0, "sig_genes": [], "concordant_genes": [], "evidence_mode": "no_data"}

        flux_id = flux_row.get("flux_id", "")
        mapping = FLUX_GENE_MAPPING.get(flux_id, {})
        genes = list(dict.fromkeys(mapping.get("genes", [])))
        gene_roles = mapping.get("gene_roles", {}) or {}
        protein_mapping = REACTION_ENZYME_MAPPING.get(flux_id, {})
        protein_direct = list(protein_mapping.get("direct_enzymes", []))
        if protein_direct:
            direct_genes = set(protein_direct)
        else:
            # 没有蛋白策展条目时才回退到基因映射中的 catalytic。
            direct_genes = set(gene_roles.get("catalytic", []))
        direct_genes.update(gene_roles.get("transporter", []))
        genes = list(dict.fromkeys(genes + list(direct_genes)))
        if not genes:
            return 0.0, {"n_genes_checked": 0, "n_sig": 0, "sig_genes": [], "concordant_genes": [], "evidence_mode": "no_mapping"}

        if not any(
            column in tx_df.columns
            for column in ("gene", "gene_symbol", "canonical_symbol")
        ):
            return 0.0, {"n_genes_checked": 0, "n_sig": 0, "sig_genes": [], "concordant_genes": [], "evidence_mode": "no_data"}

        matched = match_gene_features(tx_df, genes)
        matched_direct = match_gene_features(tx_df, list(direct_genes))
        matched_genes = list(matched.matched_queries)
        matched_direct_genes = list(matched_direct.matched_queries)
        if not matched_genes:
            return 0.0, {
                "n_genes_checked": len(genes),
                "n_genes_found": 0,
                "n_sig": 0,
                "sig_genes": [],
                "concordant_genes": [],
                "evidence_mode": "no_match",
            }

        flux_fc = flux_row.get("flux_log2fc", 0.0)
        if pd.isna(flux_fc):
            flux_fc = 0.0

        sub = tx_df.loc[list(matched.indices)].copy()
        direct_indices = set(matched_direct.indices)

        tx_cut = self.omics_cutoffs["transcriptomics"]
        tx_pval_cut = tx_cut["pval_cutoff"]
        tx_fc_cut = tx_cut["log2fc_cutoff"]

        padj_col = "padj" if "padj" in sub.columns else "pvalue"
        log2fc_col = "log2fc" if "log2fc" in sub.columns else "log2FoldChange"
        evidence_col = "evidence_level" if "evidence_level" in sub.columns else None

        # ── Phase 1: standard statistical significance ──
        sig_mask = pd.Series(True, index=sub.index)
        if padj_col in sub.columns:
            sig_mask &= sub[padj_col].fillna(1.0) < tx_pval_cut
        if log2fc_col in sub.columns:
            sig_mask &= sub[log2fc_col].fillna(0.0).abs() > tx_fc_cut
        # Exclude effect_size_only rows from statistical path
        if evidence_col and evidence_col in sub.columns:
            sig_mask &= sub[evidence_col].fillna("") != "effect_size_only"

        sig_df = sub[sig_mask]
        n_sig = len(sig_df)
        tx_evidence_mode = "statistical_test"

        # ── Phase 2: effect-size-only fallback (N=1) ──
        if n_sig == 0 and self.fallback_enabled and log2fc_col in sub.columns:
            # No statistically significant genes found → try effect-size-only
            # Use the large effect threshold as a stricter substitute for p-value
            fc_fallback_mask = sub[log2fc_col].fillna(0.0).abs() >= self.fallback_fc_large
            # Include effect_size_only rows that were excluded from statistical path
            if evidence_col and evidence_col in sub.columns:
                # These rows have padj=NaN or are marked effect_size_only
                effect_only_mask = sub[evidence_col].fillna("") == "effect_size_only"
                fc_fallback_mask |= (effect_only_mask & sub[log2fc_col].fillna(0.0).abs().ge(self.fallback_fc_large))

            sig_df = sub[fc_fallback_mask]
            n_sig = len(sig_df)
            if n_sig > 0:
                tx_evidence_mode = "effect_size_only"

        gene_column = next(
            column for column in ("gene", "gene_symbol", "canonical_symbol")
            if column in sub.columns
        )
        sig_genes_list = sorted(
            sig_df[gene_column].dropna().astype(str).tolist()
        ) if n_sig > 0 else []
        tx_stat_summary = _row_evidence_summary(
            sub,
            sig_df,
            effect_col=log2fc_col if log2fc_col in sub.columns else None,
            qvalue_candidates=["padj", "qvalue", "fdr"],
            pvalue_candidates=["pvalue", "pval"],
            n_batches_candidates=["n_batches"],
            batch_ids_candidates=["batch_ids"],
            i2_candidates=["heterogeneity_i2", "i2"],
            tau2_candidates=["heterogeneity_tau2", "tau_squared"],
            protocol_candidates=["protocol_id"],
            source_mode_candidates=["source_mode"],
        )

        if n_sig == 0:
            return 0.0, {
                **tx_stat_summary,
                "n_genes_checked": len(genes),
                "n_genes_found": len(matched_genes),
                "n_direct_genes": len(matched_direct_genes),
                "n_sig": 0,
                "sig_genes": [],
                "concordant_genes": [],
                "evidence_mode": "no_significant",
                "direction_source": (
                    "direct_gene" if matched_direct_genes else "indirect_gene"
                ),
            }

        # 只有直接酶/直接转运基因可以提供方向证据。调控和通路邻近基因仅记录为
        # association-only，避免把抑制因子误判为同向支持。
        direct_sig = sig_df.loc[
            [index for index in sig_df.index if index in direct_indices]
        ]
        indirect_sig = sig_df.loc[
            [index for index in sig_df.index if index not in direct_indices]
        ]
        dir_score = 0.0
        max_abs_fc_concordant = 0.0
        n_concordant = 0
        n_opposite = 0
        support_effect = float("nan")
        conflict_effect = float("nan")
        concordant_genes_list: list[str] = []
        direction_source = "direct_gene" if not direct_sig.empty else "indirect_gene"
        has_concordant = False

        if log2fc_col in direct_sig.columns and flux_fc != 0.0:
            direct_fcs = pd.to_numeric(direct_sig[log2fc_col], errors="coerce").dropna()
            same_dir = (
                ((direct_fcs > 0) & (flux_fc > 0))
                | ((direct_fcs < 0) & (flux_fc < 0))
            )
            if same_dir.any():
                has_concordant = True
                concordant_fcs = direct_fcs[same_dir]
                max_abs_fc_concordant = float(concordant_fcs.abs().max())
                n_concordant = int(same_dir.sum())
                support_effect = float(concordant_fcs.abs().median())
                concordant_genes_list = sorted(
                    direct_sig.loc[
                        same_dir[same_dir].index, gene_column
                    ].dropna().astype(str).tolist()
                )
            opposite_fcs = direct_fcs[~same_dir]
            n_opposite = len(opposite_fcs)
            if not opposite_fcs.empty:
                conflict_effect = float(opposite_fcs.abs().median())

        if direct_sig.empty:
            return 0.0, {
                **tx_stat_summary,
                "n_genes_checked": len(genes),
                "n_genes_found": len(matched_genes),
                "n_direct_genes": len(matched_direct_genes),
                "n_sig": n_sig,
                "n_concordant": 0,
                "n_opposite": 0,
                "n_unresolved": n_sig,
                "support_effect": float("nan"),
                "conflict_effect": float("nan"),
                "sig_genes": sig_genes_list,
                "associated_genes": sorted(
                    indirect_sig[gene_column].dropna().astype(str).tolist()
                ),
                "concordant_genes": [],
                "sig_score": 0.0,
                "dir_score": 0.0,
                "effect_score": 0.0,
                "direction_source": "indirect_gene",
                "mapping_method_counts": _mapping_method_counts(
                    matched.match_method_by_index
                ),
                "evidence_mode": tx_evidence_mode,
            }

        # 直接基因仅出现反向变化时，明确记录冲突并不给旧分数。
        if not has_concordant:
            return 0.0, {
                **tx_stat_summary,
                "n_genes_checked": len(genes),
                "n_genes_found": len(matched_genes),
                "n_direct_genes": len(matched_direct_genes),
                "n_sig": n_sig,
                "n_concordant": 0,
                "n_opposite": n_opposite,
                "n_unresolved": len(indirect_sig),
                "support_effect": float("nan"),
                "conflict_effect": conflict_effect,
                "sig_genes": sig_genes_list,
                "associated_genes": sorted(
                    indirect_sig[gene_column].dropna().astype(str).tolist()
                ),
                "concordant_genes": [],
                "sig_score": 0.0,
                "dir_score": 0.0,
                "effect_score": 0.0,
                "direction_override": "opposite",
                "direction_source": "direct_gene",
                "mapping_method_counts": _mapping_method_counts(
                    matched.match_method_by_index
                ),
                "evidence_mode": tx_evidence_mode,
            }

        # Score components — reduced for effect_size_only
        is_fallback = tx_evidence_mode == "effect_size_only"
        sig_score = 0.2 if is_fallback else 0.4
        dir_score = 0.2 if is_fallback else 0.3
        effect_score = _effect_tier(max_abs_fc_concordant, self.log2fc_tiers) * 0.3

        weight = self.weights.get("transcriptomics", 1.0)
        penalty = self.fallback_weight_penalty if is_fallback else 1.0
        total = (sig_score + dir_score + effect_score) * weight * penalty

        details = {
            **tx_stat_summary,
            "n_genes_checked": len(genes),
            "n_genes_found": len(matched_genes),
            "n_direct_genes": len(matched_direct_genes),
            "n_sig": n_sig,
            "n_concordant": n_concordant,
            "n_opposite": n_opposite,
            "n_unresolved": len(indirect_sig),
            "support_effect": support_effect,
            "conflict_effect": conflict_effect,
            "concordant_genes": concordant_genes_list,
            "direction_source": direction_source,
            "sig_genes": sig_genes_list,
            "associated_genes": sorted(
                indirect_sig[gene_column].dropna().astype(str).tolist()
            ),
            "mapping_method_counts": _mapping_method_counts(
                matched.match_method_by_index
            ),
            "sig_score": round(sig_score, 3),
            "dir_score": round(dir_score, 3),
            "effect_score": round(effect_score, 3),
            "evidence_mode": tx_evidence_mode,
            "weight_penalty": penalty,
        }
        return total, details

    def score_proteomics(
        self,
        flux_row: pd.Series,
        prot_df: pd.DataFrame | None,
        comparison: str,
    ) -> tuple[float, dict[str, Any]]:
        """Score a single flux reaction against proteomics DE data.

        Prioritizes direct_enzymes for direction concordance. Falls back to all
        mapped proteins (upstream, transporters, cofactors) with reduced credit.

        Two evidence modes:
          - statistical_test: pval/padj < cutoff AND |log2fc| > cutoff
          - effect_size_only: no p-value; |log2fc| >= fallback_fc_large

        Score components (statistical_test):
          - significance: 0.5 | effect_size_only: 0.25
          - direction:    0.3 (direct) / 0.15 (indirect) | effect_size_only: 0.2 / 0.1
          - effect:       0–0.4 (same tier)
          - weight penalty: 1.0 | effect_size_only: fallback_weight_penalty
        """
        if prot_df is None or prot_df.empty:
            return 0.0, {"n_proteins_checked": 0, "n_sig": 0, "sig_proteins": [], "concordant_proteins": [], "evidence_mode": "no_data"}

        flux_id = flux_row.get("flux_id", "")
        mapping = REACTION_ENZYME_MAPPING.get(flux_id, {})
        direct_enzymes: list[str] = list(mapping.get("direct_enzymes", []))
        directional_proteins = list(direct_enzymes)
        if not directional_proteins:
            # 运输/交换反应通常没有 EC 对应的直接酶，此时转运蛋白就是直接执行者。
            directional_proteins = list(mapping.get("transporters", []))
        directional_proteins = list(dict.fromkeys(directional_proteins))
        all_proteins: list[str] = list(direct_enzymes)
        for key in ("upstream_regulators", "downstream_enzymes",
                     "transporters", "cofactors"):
            all_proteins.extend(mapping.get(key, []))
        all_proteins = list(dict.fromkeys(all_proteins))
        if not all_proteins:
            return 0.0, {"n_proteins_checked": 0, "n_sig": 0, "sig_proteins": [], "concordant_proteins": [], "evidence_mode": "no_mapping"}

        if not any(
            column in prot_df.columns
            for column in ("protein", "Protein", "gene", "gene_symbol", "uniprot_id")
        ):
            return 0.0, {"n_proteins_checked": 0, "n_sig": 0, "sig_proteins": [], "concordant_proteins": [], "evidence_mode": "no_data"}

        matched_result = match_protein_features(prot_df, all_proteins)
        matched_direct_result = match_protein_features(
            prot_df, directional_proteins
        )
        matched = list(matched_result.matched_queries)
        matched_direct = list(matched_direct_result.matched_queries)
        if not matched:
            return 0.0, {
                "n_proteins_checked": len(all_proteins),
                "n_proteins_found": 0,
                "n_sig": 0,
                "sig_proteins": [],
                "concordant_proteins": [],
                "evidence_mode": "no_match",
            }

        sub = prot_df.loc[list(matched_result.indices)].copy()
        direct_indices = set(matched_direct_result.indices)
        protein_column = next(
            (
                column for column in ("protein", "Protein", "gene", "gene_symbol")
                if column in sub.columns
            ),
            "uniprot_id",
        )

        pval_col = f"pval_{comparison}"
        log2fc_col = f"log2fc_{comparison}"
        padj_col = f"padj_{comparison}"
        evidence_col = f"evidence_{comparison}"

        if log2fc_col not in prot_df.columns:
            return 0.0, {
                "n_proteins_checked": len(all_proteins),
                "n_proteins_found": len(matched),
                "n_direct_enzymes": len(matched_direct),
                "n_sig": 0,
                "sig_proteins": [],
                "concordant_proteins": [],
                "mapping_method_counts": _mapping_method_counts(
                    matched_result.match_method_by_index
                ),
                "data_reason": "missing_comparison_effect",
                "evidence_mode": "no_data",
            }

        prot_cut = self.omics_cutoffs["proteomics"]
        prot_pval_cut = prot_cut["pval_cutoff"]
        prot_fc_cut = prot_cut["log2fc_cutoff"]

        # ── Phase 1: standard statistical significance ──
        sig_mask = pd.Series(True, index=sub.index)
        # FDR-adjusted p-values are the primary inferential contract. Raw
        # p-values are used only for audited legacy inputs that have no finite
        # adjusted values at all.
        has_valid_pval = False
        if padj_col in sub.columns and sub[padj_col].notna().any():
            has_valid_pval = True
            sig_mask &= sub[padj_col].fillna(1.0) < prot_pval_cut
        elif pval_col in sub.columns and sub[pval_col].notna().any():
            has_valid_pval = True
            sig_mask &= sub[pval_col].fillna(1.0) < prot_pval_cut
        if evidence_col in sub.columns:
            sig_mask &= sub[evidence_col].astype(str).eq("inferential")
        if log2fc_col in sub.columns:
            sig_mask &= sub[log2fc_col].fillna(0.0).abs() > prot_fc_cut

        sig_df = sub[sig_mask]
        n_sig = len(sig_df)
        prot_evidence_mode = "statistical_test"

        # ── Phase 2: effect-size-only fallback (no valid p-value) ──
        if n_sig == 0 and self.fallback_enabled and not has_valid_pval and log2fc_col in sub.columns:
            fc_fallback_mask = sub[log2fc_col].fillna(0.0).abs() >= self.fallback_fc_large
            sig_df = sub[fc_fallback_mask]
            n_sig = len(sig_df)
            if n_sig > 0:
                prot_evidence_mode = "effect_size_only"

        sig_proteins_list = sorted(
            sig_df[protein_column].dropna().astype(str).tolist()
        ) if n_sig > 0 else []
        prot_stat_summary = _row_evidence_summary(
            sub,
            sig_df,
            effect_col=log2fc_col if log2fc_col in sub.columns else None,
            qvalue_candidates=[padj_col, "padj", "qvalue", "fdr"],
            pvalue_candidates=[pval_col, "pvalue", "pval"],
            se_candidates=[f"se_{comparison}"],
            n_batches_candidates=[f"n_batches_{comparison}"],
            batch_ids_candidates=[f"batch_ids_{comparison}"],
            i2_candidates=[f"heterogeneity_i2_{comparison}"],
            tau2_candidates=[f"heterogeneity_tau2_{comparison}"],
            protocol_candidates=[f"protocol_id_{comparison}"],
            source_mode_candidates=[f"source_mode_{comparison}"],
        )

        if n_sig == 0:
            return 0.0, {
                **prot_stat_summary,
                "n_proteins_checked": len(all_proteins),
                "n_proteins_found": len(matched),
                "n_direct_enzymes": len(matched_direct),
                "n_sig": 0,
                "sig_proteins": [],
                "concordant_proteins": [],
                "evidence_mode": "no_significant",
            }

        flux_fc = flux_row.get("flux_log2fc", 0.0)
        if pd.isna(flux_fc):
            flux_fc = 0.0

        direct_sig = sig_df.loc[
            [index for index in sig_df.index if index in direct_indices]
        ]
        indirect_sig = sig_df.loc[
            [index for index in sig_df.index if index not in direct_indices]
        ]
        dir_score = 0.0
        max_abs_fc_concordant = 0.0
        has_concordant = False
        concordant_proteins_list: list[str] = []
        n_concordant = 0
        n_opposite = 0
        support_effect = float("nan")
        conflict_effect = float("nan")
        direction_source = "direct_enzyme" if not direct_sig.empty else "indirect"
        is_fallback = prot_evidence_mode == "effect_size_only"

        if log2fc_col in direct_sig.columns and flux_fc != 0.0:
            direct_fcs = pd.to_numeric(
                direct_sig[log2fc_col], errors="coerce"
            ).dropna()
            same_dir = (
                ((direct_fcs > 0) & (flux_fc > 0))
                | ((direct_fcs < 0) & (flux_fc < 0))
            )
            if same_dir.any():
                dir_score = 0.2 if is_fallback else 0.3
                has_concordant = True
                concordant_fcs = direct_fcs[same_dir]
                max_abs_fc_concordant = float(concordant_fcs.abs().max())
                support_effect = float(concordant_fcs.abs().median())
                concordant_proteins_list = sorted(
                    direct_sig.loc[
                        same_dir[same_dir].index, protein_column
                    ].dropna().astype(str).tolist()
                )
                n_concordant = int(same_dir.sum())
            opposite_fcs = direct_fcs[~same_dir]
            n_opposite = len(opposite_fcs)
            if not opposite_fcs.empty:
                conflict_effect = float(opposite_fcs.abs().median())

        if direct_sig.empty:
            return 0.0, {
                **prot_stat_summary,
                "n_proteins_checked": len(all_proteins),
                "n_proteins_found": len(matched),
                "n_direct_enzymes": len(matched_direct),
                "n_sig": n_sig,
                "n_concordant": 0,
                "n_opposite": 0,
                "n_unresolved": n_sig,
                "support_effect": float("nan"),
                "conflict_effect": float("nan"),
                "direction_source": "indirect",
                "sig_proteins": sig_proteins_list,
                "associated_proteins": sorted(
                    indirect_sig[protein_column].dropna().astype(str).tolist()
                ),
                "concordant_proteins": [],
                "sig_score": 0.0,
                "dir_score": 0.0,
                "effect_score": 0.0,
                "mapping_method_counts": _mapping_method_counts(
                    matched_result.match_method_by_index
                ),
                "evidence_mode": prot_evidence_mode,
            }

        # 直接执行蛋白仅出现反向变化时，明确记录冲突并不给旧分数。
        if not has_concordant:
            return 0.0, {
                **prot_stat_summary,
                "n_proteins_checked": len(all_proteins),
                "n_proteins_found": len(matched),
                "n_direct_enzymes": len(matched_direct),
                "n_sig": n_sig,
                "n_concordant": 0,
                "n_opposite": n_opposite,
                "n_unresolved": len(indirect_sig),
                "support_effect": float("nan"),
                "conflict_effect": conflict_effect,
                "direction_source": "direct_enzyme",
                "sig_proteins": sig_proteins_list,
                "associated_proteins": sorted(
                    indirect_sig[protein_column].dropna().astype(str).tolist()
                ),
                "concordant_proteins": [],
                "sig_score": 0.0,
                "dir_score": 0.0,
                "effect_score": 0.0,
                "direction_override": "opposite",
                "mapping_method_counts": _mapping_method_counts(
                    matched_result.match_method_by_index
                ),
                "evidence_mode": prot_evidence_mode,
            }

        sig_score = 0.25 if is_fallback else 0.5
        effect_score = _effect_tier(max_abs_fc_concordant, self.log2fc_tiers) * 0.4

        weight = self.weights.get("proteomics", 1.2)
        penalty = self.fallback_weight_penalty if is_fallback else 1.0
        total = (sig_score + dir_score + effect_score) * weight * penalty

        details = {
            **prot_stat_summary,
            "n_proteins_checked": len(all_proteins),
            "n_proteins_found": len(matched),
            "n_direct_enzymes": len(matched_direct),
            "n_sig": n_sig,
            "n_concordant": n_concordant,
            "n_opposite": n_opposite,
            "n_unresolved": len(indirect_sig),
            "support_effect": support_effect,
            "conflict_effect": conflict_effect,
            "direction_source": direction_source,
            "sig_proteins": sig_proteins_list,
            "associated_proteins": sorted(
                indirect_sig[protein_column].dropna().astype(str).tolist()
            ),
            "concordant_proteins": concordant_proteins_list,
            "mapping_method_counts": _mapping_method_counts(
                matched_result.match_method_by_index
            ),
            "sig_score": round(sig_score, 3),
            "dir_score": round(dir_score, 3),
            "effect_score": round(effect_score, 3),
            "evidence_mode": prot_evidence_mode,
            "weight_penalty": penalty,
        }
        return total, details

    def score_metabolomics(
        self,
        flux_row: pd.Series,
        metab_df: pd.DataFrame | None,
        comparison: str,
    ) -> tuple[float, dict[str, Any]]:
        """Score a single flux reaction against metabolomics DE data.

        Uses FLUX_METABOLITE_MAPPING with substrate/product distinction.
        Direction concordance:
          - Product: flux up -> product up (same direction = concordant)
          - Substrate: flux up -> substrate down (opposite direction = concordant)
          - Related: no direction expectation (partial credit)

        Two evidence modes:
          - statistical_test: FDR/pval < cutoff AND |log2fc| > cutoff
          - effect_size_only: no p-value; |log2fc| >= fallback_fc_large

        Score components (statistical_test):
          - significance: 0.3 | effect_size_only: 0.15
          - direction:    0–0.3 | effect_size_only: 0–0.2
          - effect:       0–0.2 (same tier)
          - weight penalty: 1.0 | effect_size_only: fallback_weight_penalty
        """
        if metab_df is None or metab_df.empty:
            return 0.0, {"n_metabolites_checked": 0, "n_sig": 0, "sig_metabolites": [], "concordant_metabolites": [], "evidence_mode": "no_data"}

        flux_id = flux_row.get("flux_id", "")

        roles = get_metabolite_roles(flux_id)
        if roles is None:
            flat = FLUX_METABOLITE_MAPPING.get(flux_id, [])
            if isinstance(flat, list):
                roles = {"substrates": [], "products": [], "related": flat}
            elif isinstance(flat, dict):
                roles = flat
            else:
                return 0.0, {"n_metabolites_checked": 0, "n_sig": 0, "sig_metabolites": [], "concordant_metabolites": [], "evidence_mode": "no_mapping"}

        substrates = roles.get("substrates", [])
        products = roles.get("products", [])
        related = roles.get("related", [])
        all_metabs = substrates + products + related
        n_checked = len(all_metabs)

        if not all_metabs:
            return 0.0, {"n_metabolites_checked": 0, "n_sig": 0, "sig_metabolites": [], "concordant_metabolites": [], "evidence_mode": "no_mapping"}

        if "compound" not in metab_df.columns:
            return 0.0, {"n_metabolites_checked": n_checked, "n_sig": 0, "sig_metabolites": [], "concordant_metabolites": [], "evidence_mode": "no_data"}

        matched_roles = match_metabolite_roles(metab_df, roles)
        n_checked = sum(matched_roles.target_count_by_role.values())
        sub_idx = list(matched_roles.indices_by_role["substrates"])
        prod_idx = list(matched_roles.indices_by_role["products"])
        rel_idx = list(matched_roles.indices_by_role["related"])
        all_idx = list(dict.fromkeys(sub_idx + prod_idx + rel_idx))

        if not all_idx:
            return 0.0, {
                "n_metabolites_checked": n_checked,
                "n_metabolites_found": 0,
                "n_sig": 0,
                "sig_metabolites": [],
                "concordant_metabolites": [],
                "n_ambiguous": len(matched_roles.ambiguous_indices),
                "ambiguous_metabolites": [
                    str(metab_df.at[index, "compound"])
                    for index in matched_roles.ambiguous_indices
                ],
                "evidence_mode": "no_match",
            }

        pval_col = f"pval_{comparison}"
        log2fc_col = f"log2fc_{comparison}"
        fdr_col = f"fdr_{comparison}"

        if log2fc_col not in metab_df.columns:
            return 0.0, {
                "n_metabolites_checked": n_checked,
                "n_metabolites_found": len(all_idx),
                "n_ambiguous": len(matched_roles.ambiguous_indices),
                "n_sig": 0,
                "sig_metabolites": [],
                "concordant_metabolites": [],
                "mapping_method_counts": _mapping_method_counts(
                    matched_roles.match_method_by_index
                ),
                "data_reason": "missing_comparison_effect",
                "evidence_mode": "no_data",
            }

        metab_cut = self.omics_cutoffs["metabolomics"]
        metab_pval_cut = metab_cut["pval_cutoff"]
        metab_fc_cut = metab_cut["log2fc_cutoff"]

        # Track whether valid p-values exist for fallback detection
        _global_has_valid_pval = (
            pval_col in metab_df.columns and metab_df.loc[all_idx, pval_col].notna().any()
        )

        def _is_sig(df_subset: pd.DataFrame) -> tuple[pd.Series, str]:
            """Return (sig_mask, evidence_mode) for a metabolite subset."""
            mask = pd.Series(True, index=df_subset.index)
            has_valid_pval = (
                pval_col in df_subset.columns
                and df_subset[pval_col].notna().any()
            )
            if fdr_col in df_subset.columns:
                if has_valid_pval:
                    mask &= df_subset[fdr_col].fillna(1.0) < metab_pval_cut
                else:
                    mask &= pd.Series(False, index=df_subset.index)
            elif pval_col in df_subset.columns:
                mask &= df_subset[pval_col].fillna(1.0) < metab_pval_cut
            if log2fc_col in df_subset.columns:
                mask &= df_subset[log2fc_col].fillna(0.0).abs() > metab_fc_cut
            return mask, "statistical_test"

        def _is_sig_fallback(df_subset: pd.DataFrame) -> pd.Series:
            """Effect-size-only fallback: use log2FC large threshold only."""
            if log2fc_col not in df_subset.columns:
                return pd.Series(False, index=df_subset.index)
            return df_subset[log2fc_col].fillna(0.0).abs() >= self.fallback_fc_large

        # ── Phase 1: standard significance ──
        all_sub = metab_df.loc[all_idx]
        sig_mask, _ = _is_sig(all_sub)
        sig_df = all_sub[sig_mask]
        n_sig = len(sig_df)
        metab_evidence_mode = "statistical_test"

        # ── Phase 2: effect-size-only fallback ──
        if n_sig == 0 and self.fallback_enabled and not _global_has_valid_pval:
            fb_mask = _is_sig_fallback(all_sub)
            sig_df = all_sub[fb_mask]
            n_sig = len(sig_df)
            if n_sig > 0:
                metab_evidence_mode = "effect_size_only"

        sig_metab_list = sorted(sig_df["compound"].tolist()) if n_sig > 0 else []

        is_fallback = metab_evidence_mode == "effect_size_only"
        sig_score = (0.15 if is_fallback else 0.3) if n_sig > 0 else 0.0

        flux_fc = flux_row.get("flux_log2fc", 0.0)
        if pd.isna(flux_fc):
            flux_fc = 0.0

        dir_score = 0.0
        max_abs_fc = 0.0
        concordant_metabolites_list: list[str] = []

        def _sig_fn(df_subset: pd.DataFrame) -> pd.DataFrame:
            """Get significant rows from subset using current evidence mode."""
            if is_fallback:
                m = _is_sig_fallback(df_subset)
            else:
                m, _ = _is_sig(df_subset)
            return df_subset[m]

        if n_sig > 0 and log2fc_col in sig_df.columns and flux_fc != 0.0:
            # Products: flux up -> product up (same direction = concordant)
            prod_sub = metab_df.loc[prod_idx] if prod_idx else pd.DataFrame()
            if not prod_sub.empty:
                prod_sig = _sig_fn(prod_sub)
                if not prod_sig.empty and log2fc_col in prod_sig.columns:
                    prod_fcs = prod_sig[log2fc_col].dropna()
                    concordant = (
                        ((prod_fcs > 0) & (flux_fc > 0))
                        | ((prod_fcs < 0) & (flux_fc < 0))
                    )
                    if concordant.any():
                        dir_score = 0.2 if is_fallback else 0.3
                        max_abs_fc = max(max_abs_fc, prod_fcs[concordant].abs().max())
                        concordant_metabolites_list.extend(
                            prod_sig.loc[concordant[concordant].index, "compound"].dropna().tolist()
                        )

            # Substrates: flux up -> substrate down (opposite direction = concordant)
            if dir_score < (0.2 if is_fallback else 0.3):
                sub_only = metab_df.loc[sub_idx] if sub_idx else pd.DataFrame()
                if not sub_only.empty:
                    sub_sig = _sig_fn(sub_only)
                    if not sub_sig.empty and log2fc_col in sub_sig.columns:
                        sub_fcs = sub_sig[log2fc_col].dropna()
                        concordant = (
                            ((sub_fcs > 0) & (flux_fc < 0))
                            | ((sub_fcs < 0) & (flux_fc > 0))
                        )
                        if concordant.any():
                            dir_score = 0.2 if is_fallback else 0.3
                            max_abs_fc = max(max_abs_fc, sub_fcs[concordant].abs().max())
                            concordant_metabolites_list.extend(
                                sub_sig.loc[concordant[concordant].index, "compound"].dropna().tolist()
                            )

            # Related: partial credit without direction expectation
            dir_max = 0.2 if is_fallback else 0.3
            if dir_score < dir_max and rel_idx:
                rel_sub = metab_df.loc[rel_idx] if rel_idx else pd.DataFrame()
                if not rel_sub.empty:
                    rel_sig = _sig_fn(rel_sub)
                    if not rel_sig.empty:
                        dir_score = 0.1 if is_fallback else 0.15
                        if log2fc_col in rel_sig.columns:
                            rel_fcs = rel_sig[log2fc_col].dropna()
                            if len(rel_fcs) > 0:
                                max_abs_fc = max(max_abs_fc, rel_fcs.abs().max())

        # General direction fallback
        if dir_score == 0.0 and n_sig > 0 and log2fc_col in sig_df.columns and flux_fc != 0.0:
            metab_fcs = sig_df[log2fc_col].dropna()
            if len(metab_fcs) > 0:
                same_dir = ((metab_fcs > 0) & (flux_fc > 0)) | ((metab_fcs < 0) & (flux_fc < 0))
                if same_dir.any():
                    dir_score = 0.1 if is_fallback else 0.15
                    max_abs_fc = metab_fcs[same_dir].abs().max()
                    concordant_metabolites_list = sorted(
                        sig_df.loc[same_dir[same_dir].index, "compound"].dropna().tolist()
                    )

        # 新证据画像只把明确的底物/产物关系计为方向支持或冲突；related 保持未解析。
        direct_support_idx: list[Any] = []
        direct_conflict_idx: list[Any] = []
        unresolved_idx: list[Any] = []
        if n_sig > 0 and log2fc_col in sig_df.columns and flux_fc != 0.0:
            sub_set = set(sub_idx)
            prod_set = set(prod_idx)
            for row_idx, row_fc in pd.to_numeric(sig_df[log2fc_col], errors="coerce").dropna().items():
                if row_idx in prod_set:
                    expected = (row_fc > 0 and flux_fc > 0) or (row_fc < 0 and flux_fc < 0)
                elif row_idx in sub_set:
                    expected = (row_fc > 0 and flux_fc < 0) or (row_fc < 0 and flux_fc > 0)
                else:
                    unresolved_idx.append(row_idx)
                    continue
                if expected:
                    direct_support_idx.append(row_idx)
                else:
                    direct_conflict_idx.append(row_idx)
        elif n_sig > 0:
            unresolved_idx = list(sig_df.index)

        support_effect = float("nan")
        conflict_effect = float("nan")
        if direct_support_idx:
            support_effect = float(
                pd.to_numeric(sig_df.loc[direct_support_idx, log2fc_col], errors="coerce")
                .dropna().abs().median()
            )
        if direct_conflict_idx:
            conflict_effect = float(
                pd.to_numeric(sig_df.loc[direct_conflict_idx, log2fc_col], errors="coerce")
                .dropna().abs().median()
            )

        metab_stat_summary = _row_evidence_summary(
            all_sub,
            sig_df,
            effect_col=log2fc_col if log2fc_col in all_sub.columns else None,
            qvalue_candidates=[fdr_col, "qvalue", "fdr"],
            pvalue_candidates=[pval_col, "pvalue", "pval"],
            se_candidates=[f"se_{comparison}"],
            n_batches_candidates=[f"n_batches_{comparison}"],
            batch_ids_candidates=[f"batch_ids_{comparison}"],
            i2_candidates=[f"heterogeneity_i2_{comparison}"],
            tau2_candidates=[f"heterogeneity_tau2_{comparison}"],
            protocol_candidates=[f"protocol_id_{comparison}"],
            source_mode_candidates=[f"source_mode_{comparison}"],
        )

        effect_score = _effect_tier(max_abs_fc, self.log2fc_tiers) * 0.2

        weight = self.weights.get("metabolomics", 0.8)
        penalty = self.fallback_weight_penalty if is_fallback else 1.0
        total = (sig_score + dir_score + effect_score) * weight * penalty

        details = {
            **metab_stat_summary,
            "n_metabolites_checked": n_checked,
            "n_metabolites_found": len(all_idx),
            "n_ambiguous": len(matched_roles.ambiguous_indices),
            "ambiguous_metabolites": [
                str(metab_df.at[index, "compound"])
                for index in matched_roles.ambiguous_indices
            ],
            "mapping_method_counts": _mapping_method_counts(
                matched_roles.match_method_by_index
            ),
            "n_sig": n_sig,
            "n_concordant": len(direct_support_idx),
            "n_opposite": len(direct_conflict_idx),
            "n_unresolved": len(unresolved_idx),
            "support_effect": support_effect,
            "conflict_effect": conflict_effect,
            "direction_source": "direct" if direct_support_idx or direct_conflict_idx else "related",
            "sig_metabolites": sig_metab_list,
            "concordant_metabolites": sorted(set(concordant_metabolites_list)),
            "sig_score": round(sig_score, 3),
            "dir_score": round(dir_score, 3),
            "effect_score": round(effect_score, 3),
            "evidence_mode": metab_evidence_mode,
            "weight_penalty": penalty,
        }
        return total, details

    # ---- MFA gate and scoring ----

    def _mfa_gate_passed(self, flux_row: pd.Series) -> tuple[bool, str]:
        """Check if a flux reaction passes the MFA significance gate.

        Returns (passed, evidence_mode) where evidence_mode is one of:
          - "statistical_test"  : standard FDR/pvalue gate passed
          - "effect_size_only"  : no p-value; passed via Cohen's d + log2FC fallback
          - "insufficient"      : gate not passed

        Standard gate: FDR < 0.05 AND |Cohen's d| > 0.8.
        Legacy fallback: pvalue + log2FC.
        Effect-size-only fallback (N=1): |Cohen's d| >= threshold AND |log2FC| >= threshold.
        """
        mfa_cut = self.omics_cutoffs["mfa"]

        # Try FDR-based gate (preferred)
        fdr = flux_row.get("flux_fdr", float("nan"))
        cohen_d = flux_row.get("cohen_d", float("nan"))

        if not pd.isna(fdr) and not pd.isna(cohen_d):
            if float(fdr) < mfa_cut["pval_cutoff"] and abs(float(cohen_d)) > 0.8:
                return True, "statistical_test"
            # FDR available but not significant → check stat_note for effect-size-only context
            if self.fallback_enabled:
                stat_note = flux_row.get("stat_note", "")
                if stat_note in ("fitting_variance_only", "effect_size_only", "low_n_configured_statistical_test"):
                    cd_min = self.fallback_mfa_gate["cohen_d_minimum"]
                    fc_min = self.fallback_mfa_gate["log2fc_minimum"]
                    fcol = "flux_log2fc" if "flux_log2fc" in flux_row.index else "log2FC"
                    fc = flux_row.get(fcol, float("nan"))
                    if (not pd.isna(cohen_d) and abs(float(cohen_d)) >= cd_min
                            and not pd.isna(fc) and abs(float(fc)) >= fc_min):
                        return True, "effect_size_only"
            return False, "insufficient"

        # Legacy: pvalue + log2FC
        pcol = "flux_pvalue" if "flux_pvalue" in flux_row.index else "pvalue_ttest"
        fcol = "flux_log2fc" if "flux_log2fc" in flux_row.index else "log2FC"
        pval = flux_row.get(pcol, float("nan"))
        fc = flux_row.get(fcol, float("nan"))

        if not pd.isna(pval) and not pd.isna(fc):
            if float(pval) < mfa_cut["pval_cutoff"] and abs(float(fc)) > mfa_cut["log2fc_cutoff"]:
                return True, "statistical_test"

        # Effect-size-only fallback: no usable p-value → use Cohen's d + log2FC
        if self.fallback_enabled:
            cd_min = self.fallback_mfa_gate["cohen_d_minimum"]
            fc_min = self.fallback_mfa_gate["log2fc_minimum"]
            if (not pd.isna(cohen_d) and abs(float(cohen_d)) >= cd_min
                    and not pd.isna(fc) and abs(float(fc)) >= fc_min):
                return True, "effect_size_only"

        return False, "insufficient"

    def score_mfa(
        self,
        flux_row: pd.Series,
        mfa_evidence_mode: str = "statistical_test",
    ) -> tuple[float, dict[str, Any]]:
        """Score MFA as an independent scoring layer (0–1.0).

        Only meaningful for reactions that passed the MFA gate.

        Score components (statistical_test mode):
          - significance:  0.5 if FDR < 0.05
          - effect tier:   0–0.5 by |Cohen's d| (0/0.25/0.5)

        Score components (effect_size_only mode):
          - significance:  0.3 (reduced: no statistical test, but effect-size gate passed)
          - effect tier:   0–0.5 by |Cohen's d| (same tiers as statistical mode)
          - weight penalty: 0.6× applied to total
        """
        mfa_cut = self.omics_cutoffs["mfa"]
        is_effect_only = mfa_evidence_mode == "effect_size_only"

        fdr = flux_row.get("flux_fdr", float("nan"))
        cohen_d = flux_row.get("cohen_d", float("nan"))

        # Significance scoring
        sig_score = 0.0
        if is_effect_only:
            # No p-value available, but effect-size gate passed → reduced sig score
            sig_score = 0.3
        else:
            # Standard: FDR-based
            if not pd.isna(fdr) and float(fdr) < mfa_cut["pval_cutoff"]:
                sig_score = 0.5
            else:
                # Fallback to pvalue
                pcol = "flux_pvalue" if "flux_pvalue" in flux_row.index else "pvalue_ttest"
                pval = flux_row.get(pcol, float("nan"))
                if not pd.isna(pval) and float(pval) < mfa_cut["pval_cutoff"]:
                    sig_score = 0.5

        # Effect tier: Cohen's d based (configurable thresholds) — same for both modes
        effect_score = 0.0
        if not pd.isna(cohen_d):
            abs_cd = abs(float(cohen_d))
            if abs_cd >= self.cohen_d_tiers[1]:
                effect_score = 0.5
            elif abs_cd >= self.cohen_d_tiers[0]:
                effect_score = 0.25
        else:
            # Fallback to log2FC tier
            fcol = "flux_log2fc" if "flux_log2fc" in flux_row.index else "log2FC"
            fc = flux_row.get(fcol, float("nan"))
            effect_score = _effect_tier(float(fc) if not pd.isna(fc) else 0.0, self.log2fc_tiers) * 0.5

        weight = self.weights.get("mfa", 1.0)
        penalty = self.fallback_weight_penalty if is_effect_only else 1.0
        total = (sig_score + effect_score) * weight * penalty

        details = {
            "sig_score": round(sig_score, 3),
            "effect_score": round(effect_score, 3),
            "evidence_mode": mfa_evidence_mode,
            "weight_penalty": penalty,
            "se": _series_numeric_value(
                flux_row, ["flux_se", "log2fc_se", "standard_error", "se"]
            ),
            "n_eff": _series_numeric_value(
                flux_row, ["n_eff", "n_biological", "n_bio", "sample_size"]
            ),
            "lfsr": _series_numeric_value(
                flux_row, ["lfsr", "local_false_sign_rate"]
            ),
            "selection_probability": _series_numeric_value(
                flux_row, ["selection_probability", "bootstrap_selection_probability"]
            ),
            "direction_stability": _series_numeric_value(
                flux_row, ["direction_stability", "bootstrap_direction_stability"]
            ),
        }
        return total, details

    # ---- reaction-level aggregation ----

    def score_reaction(
        self,
        flux_row: pd.Series,
        tx_df: pd.DataFrame | None,
        prot_df: pd.DataFrame | None,
        metab_df: pd.DataFrame | None,
        comparison: str,
    ) -> dict[str, Any]:
        """同时计算旧矩阵分数与不受 MFA gate 阻断的新证据画像。"""
        # MFA gate check — returns (passed, evidence_mode)
        gate_passed, mfa_evidence_mode = self._mfa_gate_passed(flux_row)
        base_info = {
            "experiment_id": self.experiment_id,
            "estimand_id": self.estimand_id,
            "flux_id": flux_row.get("flux_id", ""),
            "tissue": flux_row.get("tissue", ""),
            "comparison": comparison,
            "flux_log2fc": flux_row.get("flux_log2fc", float("nan")),
            "flux_pvalue": flux_row.get("flux_pvalue", float("nan")),
            "flux_fdr": flux_row.get("flux_fdr", float("nan")),
            "cohen_d": flux_row.get("cohen_d", float("nan")),
            "flux_se": _series_numeric_value(
                flux_row, ["flux_se", "log2fc_se", "standard_error", "se"]
            ),
            "n_eff": _series_numeric_value(
                flux_row, ["n_eff", "n_biological", "n_bio", "sample_size"]
            ),
            "n_batches": _series_numeric_value(flux_row, ["n_batches"]),
            "batch_ids": flux_row.get("batch_ids", ""),
            "protocol_id": flux_row.get("protocol_id", ""),
            "source_mode": flux_row.get("source_mode", "single_run"),
            "heterogeneity_i2": _series_numeric_value(
                flux_row, ["heterogeneity_i2", "I_squared", "i2"]
            ),
            "heterogeneity_tau2": _series_numeric_value(
                flux_row, ["heterogeneity_tau2", "tau_squared"]
            ),
            # 通路注释透传（从 MFA 输出 → 多组学评分表）
            "kegg_pathway": flux_row.get("kegg_pathway", ""),
            "ec_number": flux_row.get("ec_number", ""),
            "kegg_reaction": flux_row.get("kegg_reaction", ""),
            "reactome_id": flux_row.get("reactome_id", ""),
        }

        # 新画像始终评估其他组学；MFA gate 只约束旧矩阵，避免低功效 MFA 抹除其他证据。
        raw_tx_score, tx_details = self.score_transcriptomics(flux_row, tx_df)
        raw_prot_score, prot_details = self.score_proteomics(flux_row, prot_df, comparison)
        raw_metab_score, metab_details = self.score_metabolomics(flux_row, metab_df, comparison)
        if gate_passed:
            raw_mfa_score, mfa_details = self.score_mfa(flux_row, mfa_evidence_mode)
            mfa_score = raw_mfa_score
            tx_score = raw_tx_score
            prot_score = raw_prot_score
            metab_score = raw_metab_score
        else:
            mfa_score = tx_score = prot_score = metab_score = 0.0
            mfa_details = {
                "sig_score": 0.0,
                "effect_score": 0.0,
                "evidence_mode": mfa_evidence_mode,
                "weight_penalty": 1.0,
                "se": base_info["flux_se"],
                "n_eff": base_info["n_eff"],
                "lfsr": _series_numeric_value(
                    flux_row, ["lfsr", "local_false_sign_rate"]
                ),
                "selection_probability": _series_numeric_value(
                    flux_row, ["selection_probability", "bootstrap_selection_probability"]
                ),
            "direction_stability": _series_numeric_value(
                    flux_row, ["direction_stability", "bootstrap_direction_stability"]
                ),
                "n_batches": base_info["n_batches"],
                "batch_ids": base_info["batch_ids"],
                "protocol_id": base_info["protocol_id"],
                "source_mode": base_info["source_mode"],
                "heterogeneity_i2": base_info["heterogeneity_i2"],
                "heterogeneity_tau2": base_info["heterogeneity_tau2"],
            }

        mfa_details.update({
            "n_batches": base_info["n_batches"],
            "batch_ids": base_info["batch_ids"],
            "protocol_id": base_info["protocol_id"],
            "source_mode": base_info["source_mode"],
            "heterogeneity_i2": base_info["heterogeneity_i2"],
            "heterogeneity_tau2": base_info["heterogeneity_tau2"],
        })

        legacy_total = mfa_score + tx_score + prot_score + metab_score
        legacy_confidence = _confidence_label(legacy_total, self.confidence_thresholds)

        # Collect per-layer evidence modes
        tx_ev = tx_details.get("evidence_mode", "no_data")
        prot_ev = prot_details.get("evidence_mode", "no_data")
        metab_ev = metab_details.get("evidence_mode", "no_data")

        # Overall evidence mode: if any contributing layer is effect_size_only, mark it
        contributing_modes = [mfa_evidence_mode] if gate_passed else []
        if raw_tx_score > 0:
            contributing_modes.append(tx_ev)
        if raw_prot_score > 0:
            contributing_modes.append(prot_ev)
        if raw_metab_score > 0:
            contributing_modes.append(metab_ev)

        if not gate_passed:
            overall_mode = "insufficient"
        elif contributing_modes and all(m == "statistical_test" for m in contributing_modes):
            overall_mode = "statistical_test"
        elif any(m == "statistical_test" for m in contributing_modes):
            overall_mode = "mixed"  # some layers statistical, some effect_size_only
        elif any(m == "effect_size_only" for m in contributing_modes):
            overall_mode = "effect_size_only"
        else:
            overall_mode = "insufficient"

        record = {
            **base_info,
            "mfa_score": round(mfa_score, 3),
            "tx_score": round(tx_score, 3),
            "prot_score": round(prot_score, 3),
            "metab_score": round(metab_score, 3),
            "total_score": round(legacy_total, 3),
            "confidence": legacy_confidence,
            "legacy_mfa_score": round(mfa_score, 3),
            "legacy_tx_score": round(tx_score, 3),
            "legacy_prot_score": round(prot_score, 3),
            "legacy_metab_score": round(metab_score, 3),
            "legacy_total_score": round(legacy_total, 3),
            "legacy_confidence": legacy_confidence,
            "gate_passed": gate_passed,
            "evidence_mode": overall_mode,
            "mfa_evidence_mode": mfa_evidence_mode,
            "tx_evidence_mode": tx_ev,
            "prot_evidence_mode": prot_ev,
            "metab_evidence_mode": metab_ev,
            "mfa_details": mfa_details,
            "tx_details": tx_details,
            "prot_details": prot_details,
            "metab_details": metab_details,
        }
        layer_records = self.evidence_profile_scorer.build_records(record)
        record.update(self.evidence_profile_scorer.summarize(layer_records))
        return record

    # ---- main entry point ----

    def score_all(
        self,
        mfa_diff: dict[str, pd.DataFrame],
        tx_data: dict[str, dict[str, pd.DataFrame]],
        prot_data: dict[str, pd.DataFrame],
        metab_data: dict[str, pd.DataFrame],
    ) -> pd.DataFrame:
        """Score all MFA differential fluxes across all comparisons and tissues.

        Parameters
        ----------
        mfa_diff : dict
            ``{comparison: DataFrame}`` from MFALoader.
        tx_data : dict
            ``{comparison: {tissue: DataFrame}}`` from TranscriptomicsLoader.
        prot_data : dict
            ``{tissue: DataFrame}`` from ProteomicsLoader.
        metab_data : dict
            ``{tissue: DataFrame}`` from MetabolomicsLoader.

        Returns
        -------
        pd.DataFrame
            One row per (flux, tissue, comparison) with scores and details.
        """
        records: list[dict[str, Any]] = []
        layer_records: list[dict[str, Any]] = []

        for comparison, mfa_df in mfa_diff.items():
            if mfa_df.empty:
                continue

            for _, flux_row in mfa_df.iterrows():
                tissue = flux_row.get("tissue", "")

                # Get omics data for this tissue+comparison
                tx_df = None
                if comparison in tx_data and tissue in tx_data[comparison]:
                    tx_df = tx_data[comparison][tissue]

                prot_df = prot_data.get(tissue)
                metab_df = metab_data.get(tissue)

                record = self.score_reaction(
                    flux_row, tx_df, prot_df, metab_df, comparison
                )
                records.append(record)
                layer_records.extend(self.evidence_profile_scorer.build_records(record))

        if not records:
            self._evidence_table = evidence_table([])
            return pd.DataFrame()

        df = pd.DataFrame(records)
        self._evidence_table = evidence_table(layer_records)
        df = df.sort_values(
            ["support_score", "conflict_score", "coverage_score", "legacy_total_score"],
            ascending=[False, True, False, False],
        ).reset_index(drop=True)

        # Summary statistics
        n_total = len(df)
        n_gated = int(df["gate_passed"].sum()) if "gate_passed" in df.columns else n_total
        n_high = (df["confidence"] == "High").sum()
        n_medium = (df["confidence"] == "Medium").sum()
        n_low = (df["confidence"] == "Low").sum()
        n_insuff = (df["confidence"] == "Insufficient").sum()

        print(f"[OK] Scoring complete: {n_total} reactions ({n_gated} passed MFA gate)")
        print(f"     Legacy High: {n_high}, Medium: {n_medium}, Low: {n_low}, Insufficient: {n_insuff}")
        grade_counts = df["evidence_grade"].value_counts().to_dict()
        print(f"     Evidence grades: {grade_counts}")

        return df

    def get_evidence_table(self) -> pd.DataFrame:
        """返回最近一次 ``score_all`` 生成的统一 EvidenceTable 副本。"""
        return self._evidence_table.copy()

    # ---- pathway-level aggregation ----

    def aggregate_by_pathway(
        self,
        scored_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """Aggregate reaction-level scores to pathway level.

        Groups by (comparison, tissue, pathway) and computes:
          - n_reactions: number of scored reactions
          - mean_score, max_score
          - n_high, n_medium, n_low: confidence counts
          - top_reactions: comma-separated list of top reaction IDs
        """
        if scored_df.empty:
            return pd.DataFrame()

        # Add pathway from annotated data, registry, or legacy mapping
        def _get_pathway(row: pd.Series) -> str:
            # 1) 优先使用 MFA 输出中已标注的 kegg_pathway
            kw = row.get("kegg_pathway", "")
            if kw and pd.notna(kw):
                return str(kw)

            # 2) 从 REACTION_REGISTRY 查找
            flux_id = row.get("flux_id", "")
            kegg_pw = get_kegg_pathway(str(flux_id))
            if kegg_pw:
                return kegg_pw

            # 3) 回退到旧的自定义 mapping
            for mapping in (FLUX_GENE_MAPPING, REACTION_ENZYME_MAPPING):
                entry = mapping.get(str(flux_id), {})
                if isinstance(entry, dict) and "pathway" in entry:
                    return entry["pathway"]
                elif isinstance(entry, dict) and "description" in entry:
                    return entry["description"]
            return "Unknown"

        scored_df = scored_df.copy()
        scored_df["pathway"] = scored_df.apply(_get_pathway, axis=1)

        # Aggregate
        def _agg_group(g: pd.DataFrame) -> pd.Series:
            if "support_score" in g.columns:
                top = (
                    g.sort_values(
                        ["support_score", "conflict_score", "coverage_score"],
                        ascending=[False, True, False],
                    )
                    .head(3)["flux_id"]
                    .tolist()
                )
            else:
                top = g.nlargest(3, "total_score")["flux_id"].tolist()
            n_gated = int(g["gate_passed"].sum()) if "gate_passed" in g.columns else len(g)
            result = {
                "n_reactions": len(g),
                "n_gated": n_gated,
                "mean_score": round(g["total_score"].mean(), 3),
                "max_score": round(g["total_score"].max(), 3),
                "n_high": (g["confidence"] == "High").sum(),
                "n_medium": (g["confidence"] == "Medium").sum(),
                "n_low": (g["confidence"] == "Low").sum(),
                "top_reactions": ", ".join(top),
            }
            if "support_score" in g.columns:
                result.update({
                    "mean_support_score": round(g["support_score"].mean(), 4),
                    "median_support_score": round(g["support_score"].median(), 4),
                    "mean_conflict_score": round(g["conflict_score"].mean(), 4),
                    "mean_coverage_score": round(g["coverage_score"].mean(), 4),
                    "n_convergent": int((g["evidence_grade"] == "Convergent").sum()),
                    "n_conflicting": int((g["evidence_grade"] == "Conflicting").sum()),
                    "n_exploratory": int((g["evidence_grade"] == "Exploratory").sum()),
                })
            return pd.Series(result)

        group_columns = ["comparison", "tissue", "pathway"]
        pathway_rows: list[dict[str, Any]] = []
        for group_key, group_df in scored_df.groupby(
            group_columns, sort=True, dropna=False
        ):
            aggregate = _agg_group(group_df).to_dict()
            identity = dict(zip(group_columns, group_key))
            pathway_rows.append({**identity, **aggregate})
        pathway_df = pd.DataFrame(pathway_rows)
        if "mean_support_score" in pathway_df.columns:
            pathway_df = pathway_df.sort_values(
                ["mean_support_score", "mean_conflict_score", "mean_coverage_score"],
                ascending=[False, True, False],
            ).reset_index(drop=True)
        else:
            pathway_df = pathway_df.sort_values("mean_score", ascending=False).reset_index(drop=True)

        print(f"[OK] Pathway aggregation: {len(pathway_df)} pathway entries")
        return pathway_df
