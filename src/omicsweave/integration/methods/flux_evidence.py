#!/usr/bin/env python3
"""
以 MFA 差异通量为锚点的多组学证据构建（flux-anchored evidence builder）。

对每个（反应, 组织, 比较），按反应 → 基因 / 蛋白 / 代谢物的策展映射，从各组学的差异结果中提取证据：
映射与匹配情况、显著特征、与通量方向一致或相反的直接执行者（酶、转运体、底物 / 产物）及其效应量、
统计可信度与数据来源。证据交给 :class:`~.evidence_profile.EvidenceProfileScorer` 汇总为多维证据画像
（支持、冲突、覆盖、稳定性、设计等级、证据等级）与长表 EvidenceTable，这是唯一的整合结果。

旧版的 0–4 加和总分与 High / Medium / Low 置信度已删除（v2 设计文档第 10.3 节，阶段 2.8）。
MFA 显著性判定（FDR < 0.05 且 |Cohen's d| > 0.8；样本量不足时改用效应量判定）只决定 MFA 层能否作为方向锚点，
不影响其他组学的证据评估。拟合点层面的推断标记为 ``pseudo_replicate``，只给探索性设计等级。
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
from omicsweave.core.idmap.reaction_pathway_mapping import get_kegg_pathway
from ..mappings.network_annotation import reaction_semantics
from .evidence_profile import EvidenceProfileScorer, evidence_table


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

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
    n_eff_candidates: list[str] | None = None,
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
        (n_eff_candidates or []) + ["n_eff", "effective_n", "sample_size", "n_samples", "n"],
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
# FluxEvidenceBuilder
# ---------------------------------------------------------------------------

class FluxEvidenceBuilder:
    """为每个 MFA 差异通量构建多组学证据，并汇总为多维证据画像。

    Parameters
    ----------
    config : object
        提供 ``experiment_id``、``estimand``（可选）与 ``scoring``（dict）：``pval_cutoff`` /
        ``log2fc_cutoff`` 与 ``omics_cutoffs`` 决定各组学特征是否显著；``fallback_mode`` 决定样本量
        不足时的效应量判定；``evidence_profile`` 传给证据画像。
    """

    def __init__(self, config: Any) -> None:
        self.mapping_audit = audit_mapping_integrity()
        self.mapping_audit.raise_for_errors()
        self.experiment_id = str(getattr(config, "experiment_id", ""))
        estimand = getattr(config, "estimand", {}) or {}
        self.estimand_id = str(estimand.get("id", ""))
        scoring = getattr(config, "scoring", {}) or {}
        self.pval_cutoff: float = scoring.get("pval_cutoff", 0.05)
        self.log2fc_cutoff: float = scoring.get("log2fc_cutoff", 0.58)

        # 各组学可单独覆盖全局的 pval_cutoff / log2fc_cutoff
        global_defaults = {"pval_cutoff": self.pval_cutoff, "log2fc_cutoff": self.log2fc_cutoff}
        omics_cutoffs = scoring.get("omics_cutoffs", {})
        self.omics_cutoffs: dict[str, dict[str, float]] = {}
        for key in ("mfa", "transcriptomics", "proteomics", "metabolomics"):
            self.omics_cutoffs[key] = {**global_defaults, **omics_cutoffs.get(key, {})}

        # 样本量不足（如每组 1 只）时的效应量判定；阈值默认取 effect_tiers 的下界与上界
        effect_tiers = scoring.get("effect_tiers", {})
        log2fc_tiers: list[float] = effect_tiers.get("log2fc", [0.58, 1.0])
        cohen_d_tiers: list[float] = effect_tiers.get("cohen_d", [0.8, 1.5])
        fallback = scoring.get("fallback_mode", {})
        self.fallback_enabled: bool = fallback.get("enabled", True)
        self.fallback_mfa_gate: dict[str, float] = fallback.get("mfa_gate", {
            "cohen_d_minimum": cohen_d_tiers[0],
            "log2fc_minimum": log2fc_tiers[0],
        })
        self.fallback_fc_large: float = fallback.get("effect_size_only_thresholds", {}).get(
            "log2fc_large", log2fc_tiers[1],
        )

        self.evidence_profile_scorer = EvidenceProfileScorer(scoring.get("evidence_profile", {}))
        self._evidence_table = evidence_table([])

    # ---- 各组学证据 ----

    def transcriptomics_evidence(
        self,
        flux_row: pd.Series,
        tx_df: pd.DataFrame | None,
    ) -> dict[str, Any]:
        """一个反应在转录组中的证据：映射基因的匹配、显著基因与方向。

        显著：padj < cutoff 且 |log2FC| > cutoff；没有统计显著基因时，|log2FC| ≥ fallback 阈值的基因按
        ``effect_size_only`` 记录。只有直接酶 / 转运体基因提供方向证据，调控与通路邻近基因只记为关联。
        """
        if tx_df is None or tx_df.empty:
            return {"n_genes_checked": 0, "n_sig": 0, "sig_genes": [], "concordant_genes": [], "evidence_mode": "no_data"}

        flux_id = flux_row.get("flux_id", "")
        mapping = FLUX_GENE_MAPPING.get(flux_id, {})
        genes = list(dict.fromkeys(mapping.get("genes", [])))
        gene_roles = mapping.get("gene_roles", {}) or {}
        protein_mapping = REACTION_ENZYME_MAPPING.get(flux_id, {})
        protein_direct = list(protein_mapping.get("direct_enzymes", []))
        # 方向成员与蛋白层一致：直接酶（蛋白表优先，没有蛋白条目时用基因表的 catalytic）；
        # 只有没有直接酶的转运 / 交换反应，才由载体充当直接执行者。
        direct_genes = set(protein_direct) if protein_direct else set(gene_roles.get("catalytic", []))
        carriers = set(protein_mapping.get("transporters", [])) | set(gene_roles.get("transporter", []))
        if not direct_genes:
            direct_genes = set(carriers)
        genes = list(dict.fromkeys(genes + sorted(direct_genes) + sorted(carriers)))
        if not reaction_semantics(str(flux_id)).directional:
            direct_genes = set()           # 只表示容量：匹配到的基因全部记为关联，不判支持或冲突
        if not genes:
            return {"n_genes_checked": 0, "n_sig": 0, "sig_genes": [], "concordant_genes": [], "evidence_mode": "no_mapping"}

        if not any(
            column in tx_df.columns
            for column in ("gene", "gene_symbol", "canonical_symbol")
        ):
            return {"n_genes_checked": 0, "n_sig": 0, "sig_genes": [], "concordant_genes": [], "evidence_mode": "no_data"}

        matched = match_gene_features(tx_df, genes)
        matched_direct = match_gene_features(tx_df, list(direct_genes))
        matched_genes = list(matched.matched_queries)
        matched_direct_genes = list(matched_direct.matched_queries)
        if not matched_genes:
            return {
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
            return {
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
            return {
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
                "direction_source": "indirect_gene",
                "mapping_method_counts": _mapping_method_counts(
                    matched.match_method_by_index
                ),
                "evidence_mode": tx_evidence_mode,
            }

        # 直接基因仅出现反向变化时，明确记录冲突并不给旧分数。
        if not has_concordant:
            return {
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
                "direction_override": "opposite",
                "direction_source": "direct_gene",
                "mapping_method_counts": _mapping_method_counts(
                    matched.match_method_by_index
                ),
                "evidence_mode": tx_evidence_mode,
            }

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
            "evidence_mode": tx_evidence_mode,
        }
        return details

    def proteomics_evidence(
        self,
        flux_row: pd.Series,
        prot_df: pd.DataFrame | None,
        comparison: str,
    ) -> dict[str, Any]:
        """一个反应在蛋白组中的证据：方向证据优先取直接酶，没有直接酶时取转运体；
        上游调控、下游酶与辅因子只计入匹配与关联。显著与 ``effect_size_only`` 的判定同转录组。
        """
        if prot_df is None or prot_df.empty:
            return {"n_proteins_checked": 0, "n_sig": 0, "sig_proteins": [], "concordant_proteins": [], "evidence_mode": "no_data"}

        flux_id = flux_row.get("flux_id", "")
        mapping = REACTION_ENZYME_MAPPING.get(flux_id, {})
        direct_enzymes: list[str] = list(mapping.get("direct_enzymes", []))
        directional_proteins = list(direct_enzymes)
        if not directional_proteins:
            # 运输/交换反应通常没有 EC 对应的直接酶，此时转运蛋白就是直接执行者。
            directional_proteins = list(mapping.get("transporters", []))
        directional_proteins = list(dict.fromkeys(directional_proteins))
        if not reaction_semantics(str(flux_id)).directional:
            directional_proteins = []      # 只表示容量：匹配到的蛋白全部记为关联
        all_proteins: list[str] = list(direct_enzymes)
        for key in ("upstream_regulators", "downstream_enzymes",
                     "transporters", "cofactors", "related_isozymes"):
            all_proteins.extend(mapping.get(key, []))
        all_proteins = list(dict.fromkeys(all_proteins))
        if not all_proteins:
            return {"n_proteins_checked": 0, "n_sig": 0, "sig_proteins": [], "concordant_proteins": [], "evidence_mode": "no_mapping"}

        if not any(
            column in prot_df.columns
            for column in ("protein", "Protein", "gene", "gene_symbol", "uniprot_id")
        ):
            return {"n_proteins_checked": 0, "n_sig": 0, "sig_proteins": [], "concordant_proteins": [], "evidence_mode": "no_data"}

        matched_result = match_protein_features(prot_df, all_proteins)
        matched_direct_result = match_protein_features(
            prot_df, directional_proteins
        )
        matched = list(matched_result.matched_queries)
        matched_direct = list(matched_direct_result.matched_queries)
        if not matched:
            return {
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
            return {
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
            n_eff_candidates=[f"n_eff_{comparison}"],
            n_batches_candidates=[f"n_batches_{comparison}"],
            batch_ids_candidates=[f"batch_ids_{comparison}"],
            i2_candidates=[f"heterogeneity_i2_{comparison}"],
            tau2_candidates=[f"heterogeneity_tau2_{comparison}"],
            protocol_candidates=[f"protocol_id_{comparison}"],
            source_mode_candidates=[f"source_mode_{comparison}"],
        )

        if n_sig == 0:
            return {
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
        has_concordant = False
        concordant_proteins_list: list[str] = []
        n_concordant = 0
        n_opposite = 0
        support_effect = float("nan")
        conflict_effect = float("nan")
        direction_source = "direct_enzyme" if not direct_sig.empty else "indirect"

        if log2fc_col in direct_sig.columns and flux_fc != 0.0:
            direct_fcs = pd.to_numeric(
                direct_sig[log2fc_col], errors="coerce"
            ).dropna()
            same_dir = (
                ((direct_fcs > 0) & (flux_fc > 0))
                | ((direct_fcs < 0) & (flux_fc < 0))
            )
            if same_dir.any():
                has_concordant = True
                concordant_fcs = direct_fcs[same_dir]
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
            return {
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
                "mapping_method_counts": _mapping_method_counts(
                    matched_result.match_method_by_index
                ),
                "evidence_mode": prot_evidence_mode,
            }

        # 直接执行蛋白仅出现反向变化时，明确记录冲突并不给旧分数。
        if not has_concordant:
            return {
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
                "direction_override": "opposite",
                "mapping_method_counts": _mapping_method_counts(
                    matched_result.match_method_by_index
                ),
                "evidence_mode": prot_evidence_mode,
            }

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
            "evidence_mode": prot_evidence_mode,
        }
        return details

    def metabolomics_evidence(
        self,
        flux_row: pd.Series,
        metab_df: pd.DataFrame | None,
        comparison: str,
    ) -> dict[str, Any]:
        """一个反应在代谢组中的证据：产物与通量同向、底物与通量反向为支持；关联代谢物不提供方向。
        显著与 ``effect_size_only`` 的判定同转录组。
        """
        if metab_df is None or metab_df.empty:
            return {"n_metabolites_checked": 0, "n_sig": 0, "sig_metabolites": [], "concordant_metabolites": [], "evidence_mode": "no_data"}

        flux_id = flux_row.get("flux_id", "")

        roles = get_metabolite_roles(flux_id)
        if roles is None:
            flat = FLUX_METABOLITE_MAPPING.get(flux_id, [])
            if isinstance(flat, list):
                roles = {"substrates": [], "products": [], "related": flat}
            elif isinstance(flat, dict):
                roles = flat
            else:
                return {"n_metabolites_checked": 0, "n_sig": 0, "sig_metabolites": [], "concordant_metabolites": [], "evidence_mode": "no_mapping"}

        substrates = roles.get("substrates", [])
        products = roles.get("products", [])
        related = roles.get("related", [])
        all_metabs = substrates + products + related
        n_checked = len(all_metabs)

        if not all_metabs:
            return {"n_metabolites_checked": 0, "n_sig": 0, "sig_metabolites": [], "concordant_metabolites": [], "evidence_mode": "no_mapping"}

        if "compound" not in metab_df.columns:
            return {"n_metabolites_checked": n_checked, "n_sig": 0, "sig_metabolites": [], "concordant_metabolites": [], "evidence_mode": "no_data"}

        matched_roles = match_metabolite_roles(metab_df, roles)
        n_checked = sum(matched_roles.target_count_by_role.values())
        sub_idx = list(matched_roles.indices_by_role["substrates"])
        prod_idx = list(matched_roles.indices_by_role["products"])
        rel_idx = list(matched_roles.indices_by_role["related"])
        all_idx = list(dict.fromkeys(sub_idx + prod_idx + rel_idx))

        if not all_idx:
            return {
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
            return {
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

        flux_fc = flux_row.get("flux_log2fc", 0.0)
        if pd.isna(flux_fc):
            flux_fc = 0.0

        # 只有底物 / 产物给出方向：产物与通量同向、底物与通量反向为支持，反之为冲突；关联代谢物保持未解析。
        direct_support_idx: list[Any] = []
        direct_conflict_idx: list[Any] = []
        unresolved_idx: list[Any] = []
        if n_sig > 0 and log2fc_col in sig_df.columns and flux_fc != 0.0:
            # 只表示容量的反应（近平衡、转运、交换）：池大小不提示净方向，底物 / 产物也只记关联
            directional = reaction_semantics(str(flux_id)).directional
            sub_set = set(sub_idx) if directional else set()
            prod_set = set(prod_idx) if directional else set()
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
            n_eff_candidates=[f"n_eff_{comparison}"],
            n_batches_candidates=[f"n_batches_{comparison}"],
            batch_ids_candidates=[f"batch_ids_{comparison}"],
            i2_candidates=[f"heterogeneity_i2_{comparison}"],
            tau2_candidates=[f"heterogeneity_tau2_{comparison}"],
            protocol_candidates=[f"protocol_id_{comparison}"],
            source_mode_candidates=[f"source_mode_{comparison}"],
        )

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
            "concordant_metabolites": sorted(set(
                sig_df.loc[direct_support_idx, "compound"].dropna().astype(str))),
            "evidence_mode": metab_evidence_mode,
        }
        return details

    # ---- MFA ----

    def _mfa_gate_passed(self, flux_row: pd.Series) -> tuple[bool, str]:
        """MFA 显著性判定：通过的反应作为方向锚点（证据画像中 MFA 层为 support）。

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

    def mfa_evidence(
        self,
        flux_row: pd.Series,
        mfa_evidence_mode: str = "statistical_test",
    ) -> dict[str, Any]:
        """MFA 层的连续统计量：标准误、有效样本量、lfsr 与稳定性校准。"""
        return {
            "evidence_mode": mfa_evidence_mode,
            "se": _series_numeric_value(flux_row, ["flux_se", "log2fc_se", "standard_error", "se"]),
            "n_eff": _series_numeric_value(flux_row, ["n_eff", "n_biological", "n_bio", "sample_size"]),
            "lfsr": _series_numeric_value(flux_row, ["lfsr", "local_false_sign_rate"]),
            "selection_probability": _series_numeric_value(
                flux_row, ["selection_probability", "bootstrap_selection_probability"]
            ),
            "direction_stability": _series_numeric_value(
                flux_row, ["direction_stability", "bootstrap_direction_stability"]
            ),
        }

    # ---- 反应层 ----

    def build_record(
        self,
        flux_row: pd.Series,
        tx_df: pd.DataFrame | None,
        prot_df: pd.DataFrame | None,
        metab_df: pd.DataFrame | None,
        comparison: str,
    ) -> dict[str, Any]:
        """一个（反应, 组织, 比较）的各层证据与证据画像；其他组学的证据不受 MFA 显著性判定影响。"""
        gate_passed, mfa_evidence_mode = self._mfa_gate_passed(flux_row)
        semantics = reaction_semantics(str(flux_row.get("flux_id", "")))
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
            "flux_se": _series_numeric_value(flux_row, ["flux_se", "log2fc_se", "standard_error", "se"]),
            "n_eff": _series_numeric_value(flux_row, ["n_eff", "n_biological", "n_bio", "sample_size"]),
            "n_batches": _series_numeric_value(flux_row, ["n_batches"]),
            "batch_ids": flux_row.get("batch_ids", ""),
            "protocol_id": flux_row.get("protocol_id", ""),
            "source_mode": flux_row.get("source_mode", "single_run"),
            "heterogeneity_i2": _series_numeric_value(flux_row, ["heterogeneity_i2", "I_squared", "i2"]),
            "heterogeneity_tau2": _series_numeric_value(flux_row, ["heterogeneity_tau2", "tau_squared"]),
            # 通路注释透传（MFA 输出 → 证据画像）
            "kegg_pathway": flux_row.get("kegg_pathway", ""),
            "ec_number": flux_row.get("ec_number", ""),
            "kegg_reaction": flux_row.get("kegg_reaction", ""),
            "reactome_id": flux_row.get("reactome_id", ""),
            # 网络注释：反应类别与各组学层能否给出方向证据
            "reaction_class": semantics.reaction_class,
            "direction_evidence": semantics.direction_evidence,
        }

        tx_details = self.transcriptomics_evidence(flux_row, tx_df)
        prot_details = self.proteomics_evidence(flux_row, prot_df, comparison)
        metab_details = self.metabolomics_evidence(flux_row, metab_df, comparison)
        mfa_details = self.mfa_evidence(flux_row, mfa_evidence_mode)
        mfa_details.update({key: base_info[key] for key in (
            "n_batches", "batch_ids", "protocol_id", "source_mode", "heterogeneity_i2", "heterogeneity_tau2")})

        tx_ev = tx_details.get("evidence_mode", "no_data")
        prot_ev = prot_details.get("evidence_mode", "no_data")
        metab_ev = metab_details.get("evidence_mode", "no_data")
        # 整体证据路径：MFA 显著时，看给出方向一致证据的各层（代谢组为有显著代谢物）用的是统计检验还是效应量
        contributing_modes = [mfa_evidence_mode] if gate_passed else []
        if tx_details.get("n_concordant", 0) > 0:
            contributing_modes.append(tx_ev)
        if prot_details.get("n_concordant", 0) > 0:
            contributing_modes.append(prot_ev)
        if metab_details.get("n_sig", 0) > 0:
            contributing_modes.append(metab_ev)
        if not gate_passed:
            overall_mode = "insufficient"
        elif contributing_modes and all(m == "statistical_test" for m in contributing_modes):
            overall_mode = "statistical_test"
        elif any(m == "statistical_test" for m in contributing_modes):
            overall_mode = "mixed"
        elif any(m == "effect_size_only" for m in contributing_modes):
            overall_mode = "effect_size_only"
        else:
            overall_mode = "insufficient"

        record = {
            **base_info,
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
        record.update(self.evidence_profile_scorer.summarize(self.evidence_profile_scorer.build_records(record)))
        if not semantics.assessable:
            # 虚拟反应与只在血浆中的反应：组织组学无从评估，不参与证据排序
            record["evidence_grade"] = "Not-assessable"
        return record

    # ---- main entry point ----

    def build_all(
        self,
        mfa_diff: dict[str, pd.DataFrame],
        tx_data: dict[str, dict[str, pd.DataFrame]],
        prot_data: dict[str, pd.DataFrame],
        metab_data: dict[str, pd.DataFrame],
    ) -> pd.DataFrame:
        """为所有比较、组织中的每个 MFA 差异通量构建证据画像。

        Parameters
        ----------
        mfa_diff : ``{比较: DataFrame}``，每行一个（反应, 组织）。
        tx_data : ``{比较: {组织: DataFrame}}``。
        prot_data, metab_data : ``{组织: DataFrame}``，各比较的统计量在 ``log2fc_<比较>`` 等列中。

        Returns
        -------
        pd.DataFrame
            每个（反应, 组织, 比较）一行，按支持度降序、冲突度升序、覆盖度降序排列。
        """
        records: list[dict[str, Any]] = []
        layer_records: list[dict[str, Any]] = []
        for comparison, mfa_df in mfa_diff.items():
            for _, flux_row in mfa_df.iterrows():
                tissue = flux_row.get("tissue", "")
                tx_df = tx_data.get(comparison, {}).get(tissue)
                record = self.build_record(flux_row, tx_df, prot_data.get(tissue), metab_data.get(tissue), comparison)
                records.append(record)
                layer_records.extend(self.evidence_profile_scorer.build_records(record))

        self._evidence_table = evidence_table(layer_records)
        if not records:
            return pd.DataFrame()
        df = pd.DataFrame(records).sort_values(
            ["support_score", "conflict_score", "coverage_score"], ascending=[False, True, False],
        ).reset_index(drop=True)
        n_gated = int(df["gate_passed"].sum())
        print(f"[OK] Flux-anchored evidence: {len(df)} reactions ({n_gated} significant in MFA)")
        print(f"     Evidence grades: {df['evidence_grade'].value_counts().to_dict()}")
        return df

    def get_evidence_table(self) -> pd.DataFrame:
        """返回最近一次 ``build_all`` 生成的统一 EvidenceTable 副本。"""
        return self._evidence_table.copy()

    # ---- pathway-level aggregation ----

    def aggregate_by_pathway(
        self,
        scored_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """把反应层证据画像汇总到通路：按（比较, 组织, 通路）统计反应数、MFA 显著数、平均支持 / 冲突 /
        覆盖度、各证据等级的反应数与支持度最高的 3 个反应。"""
        if scored_df.empty:
            return pd.DataFrame()

        # 通路取 MFA 输出中的 kegg_pathway，其次反应注册表，最后旧版自定义映射
        def _get_pathway(row: pd.Series) -> str:
            kw = row.get("kegg_pathway", "")
            if kw and pd.notna(kw):
                return str(kw)
            flux_id = row.get("flux_id", "")
            kegg_pw = get_kegg_pathway(str(flux_id))
            if kegg_pw:
                return kegg_pw
            for mapping in (FLUX_GENE_MAPPING, REACTION_ENZYME_MAPPING):
                entry = mapping.get(str(flux_id), {})
                if isinstance(entry, dict) and "pathway" in entry:
                    return entry["pathway"]
                elif isinstance(entry, dict) and "description" in entry:
                    return entry["description"]
            return "Unknown"

        scored_df = scored_df.copy()
        scored_df["pathway"] = scored_df.apply(_get_pathway, axis=1)

        def _agg_group(g: pd.DataFrame) -> pd.Series:
            top = (g.sort_values(["support_score", "conflict_score", "coverage_score"],
                                 ascending=[False, True, False]).head(3)["flux_id"].tolist())
            return pd.Series({
                "n_reactions": len(g),
                "n_gated": int(g["gate_passed"].sum()),
                "top_reactions": ", ".join(top),
                "mean_support_score": round(g["support_score"].mean(), 4),
                "median_support_score": round(g["support_score"].median(), 4),
                "mean_conflict_score": round(g["conflict_score"].mean(), 4),
                "mean_coverage_score": round(g["coverage_score"].mean(), 4),
                "n_convergent": int((g["evidence_grade"] == "Convergent").sum()),
                "n_conflicting": int((g["evidence_grade"] == "Conflicting").sum()),
                "n_exploratory": int((g["evidence_grade"] == "Exploratory").sum()),
            })

        group_columns = ["comparison", "tissue", "pathway"]
        pathway_rows: list[dict[str, Any]] = []
        for group_key, group_df in scored_df.groupby(
            group_columns, sort=True, dropna=False
        ):
            aggregate = _agg_group(group_df).to_dict()
            identity = dict(zip(group_columns, group_key))
            pathway_rows.append({**identity, **aggregate})
        pathway_df = pd.DataFrame(pathway_rows).sort_values(
            ["mean_support_score", "mean_conflict_score", "mean_coverage_score"],
            ascending=[False, True, False],
        ).reset_index(drop=True)

        print(f"[OK] Pathway aggregation: {len(pathway_df)} pathway entries")
        return pathway_df
