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

from collections import Counter, defaultdict
from typing import Any

import numpy as np
import pandas as pd

from ..mappings import (
    FLUX_GENE_MAPPING,
    REACTION_ENZYME_MAPPING,
    match_gene_features,
    match_metabolite_roles,
    match_protein_features,
    audit_mapping_integrity,
)
from ..mappings.feature_resolver import (
    IdentifierMatch,
    _canonical_gene_symbol,
    measured_metabolites,
    metabolite_member_of,
    metabolite_members,
)
from ..mappings.flux_protein_mapping import PROTEIN_CATEGORIES
from ..mappings.reaction_members import GENE_ROLE_ROLE, PROTEIN_CATEGORY_ROLE
from omicsweave.core.idmap.reaction_pathway_mapping import get_kegg_pathway
from ..mappings.model_package import DEFAULT_PACKAGE
from ..mappings.network_annotation import lipid_pool_roles
from .evidence_profile import ANCHORED_STATES, EvidenceProfileScorer, evidence_table


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def _derived_indices(flux_id: str, metab_df: pd.DataFrame,
                     pools: dict[str, set[str]] | None = None) -> dict[str, list[Any]]:
    """代谢组派生特征在一个反应中的角色：脂质类别总量按网络中的脂质池作底物 / 产物，质量作用比按
    ``related_reaction`` 记为关联（比值反映氧化还原等状态，不提示净通量方向）。"""
    out: dict[str, list[Any]] = {"substrates": [], "products": [], "related": []}
    if metab_df is None or "feature_type" not in metab_df.columns:
        return out
    types = metab_df["feature_type"].fillna("").astype(str)
    pools = lipid_pool_roles(flux_id) if pools is None else pools
    if pools["substrates"] or pools["products"]:
        ids = metab_df["metabolite_id"].fillna("").astype(str) if "metabolite_id" in metab_df.columns \
            else pd.Series("", index=metab_df.index)
        for index in metab_df.index[types.eq("lipid_class_total")]:
            lipid_class = ids[index].split(":", 1)[-1]
            if lipid_class in pools["substrates"]:
                out["substrates"].append(index)
            elif lipid_class in pools["products"]:
                out["products"].append(index)
    if "related_reaction" in metab_df.columns:
        related = metab_df["related_reaction"].fillna("").astype(str)
        # 质量作用比与由网络生成的一步反应产物 / 底物比（G8）都只记为关联
        ratio = types.isin(["mass_action_ratio", "step_ratio"])
        out["related"] += list(metab_df.index[ratio & related.eq(flux_id)])
    return out

def _implied_counts(signs: pd.Series) -> dict[str, int]:
    """方向成员各自提示的通量方向（+1 升、-1 降）计数；与通量点估计的符号无关，供无锚点时的组学一致性判定。"""
    values = pd.to_numeric(signs, errors="coerce").dropna()
    return {"n_implied_up": int((values > 0).sum()), "n_implied_down": int((values < 0).sum())}


def _isozyme_check(direct: pd.DataFrame, fc_col: str, abundance_col: str, name_col: str, flux_fc: float,
                   concordant: list[str], minor_share: float) -> dict[str, Any]:
    """同工酶加权（G7）：方向成员按该组织中的平均丰度占比加权，给出容量变化的加权值。

    同时计算分子/分母组的条件活跃占比（``max(pooled_share, share_num, share_den)``），避免强诱导型同工酶在整体均值中被稀释误杀。
    与通量同向的成员全部在整体及活跃条件下占比都低于 ``minor_share``，而占比最高的主要同工酶方向相反时，
    记 ``minor_isozyme_driven``（调用方据此不计为支持）。没有丰度列或匹配到的方向成员少于 2 个时不判断。
    """
    if abundance_col not in direct.columns or fc_col not in direct.columns or len(direct) < 2:
        return {}
    abundance = pd.to_numeric(direct[abundance_col], errors="coerce")
    fc = pd.to_numeric(direct[fc_col], errors="coerce")
    usable = abundance.gt(0) & fc.notna()
    if usable.sum() < 2:
        return {}
    a_use = abundance[usable]
    fc_use = fc[usable]
    share = a_use / a_use.sum()
    ratio = np.power(2.0, fc_use.clip(-15.0, 15.0))
    a_num = 2.0 * a_use * ratio / (1.0 + ratio)
    a_den = 2.0 * a_use / (1.0 + ratio)
    share_num = a_num / a_num.sum() if a_num.sum() > 0 else share
    share_den = a_den / a_den.sum() if a_den.sum() > 0 else share
    active_share = np.maximum(share, np.maximum(share_num, share_den))
    names = direct.loc[share.index, name_col].astype(str)
    major = share.idxmax()
    concordant_share = active_share[names.isin(concordant)]
    driven = bool(not concordant_share.empty and concordant_share.max() < minor_share
                  and names[major] not in concordant and fc[major] != 0 and np.sign(fc[major]) == -np.sign(flux_fc))
    return {"weighted_capacity_change": float((share * fc_use).sum()), "major_member": names[major],
            "major_share": round(float(share.max()), 4), "minor_isozyme_driven": driven}


# 反应档案（G13）的成员行：FluxEvidenceBuilder.member_matches 的列
MEMBER_COLUMNS = ["layer", "role", "member_id", "member_relation", "regulator_effect", "regulator_mode",
                  "direction_member", "row", "mapping_method", "selected", "selection_mode", "in_flux_anchored"]
_GENE_COLUMNS = ("gene", "gene_symbol", "canonical_symbol")
_PROTEIN_COLUMNS = ("protein", "Protein", "gene", "gene_symbol", "uniprot_id")


def _metabolite_roles(roles: Any) -> dict[str, list[str]]:
    """模型包给出的代谢物角色 → {substrates, products, related}（旧映射可能只给一个列表，全部记为关联）。"""
    if isinstance(roles, dict):
        return {key: list(roles.get(key, []) or []) for key in ("substrates", "products", "related")}
    return {"substrates": [], "products": [], "related": list(roles or [])}


def _tx_columns(frame: pd.DataFrame) -> tuple[str, str, str | None]:
    """转录组表的 (FDR 列, log2FC 列, 证据级别列)。"""
    return ("padj" if "padj" in frame.columns else "pvalue",
            "log2fc" if "log2fc" in frame.columns else "log2FoldChange",
            "evidence_level" if "evidence_level" in frame.columns else None)


def _non_isozyme_rows(frame: pd.DataFrame, relations: dict[str, str], proteins: bool) -> set[Any]:
    """登记为复合体亚基或共用亚基的成员在数据中的行：同工酶加权只在可互相替代的同工酶之间进行（G17）。"""
    members = [member for member, relation in relations.items() if relation != "isozyme"]
    if not members:
        return set()
    return set((match_protein_features if proteins else match_gene_features)(frame, members).indices)


def _demote_minor_support(details: dict[str, Any], iso: dict[str, Any], concordant_key: str) -> dict[str, Any]:
    """次要同工酶驱动的支持改记为未解析（不判支持）。"""
    details.update(iso)
    if iso.get("minor_isozyme_driven"):
        moved = int(details.get("n_concordant", 0) or 0)
        details.update({"n_concordant": 0, concordant_key: [], "support_effect": float("nan"),
                        "n_unresolved": int(details.get("n_unresolved", 0) or 0) + moved})
    return details


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
    lfsr_candidates: list[str] | None = None,
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

    lfsr_col = _first_numeric_column(effect_source, (lfsr_candidates or []) + ["lfsr", "local_false_sign_rate"])
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
        # 反应知识（成员、语义、脂质池）来自模型包；默认为现有的策展映射（mappings.model_package）
        self.package = getattr(config, "package", None) or DEFAULT_PACKAGE
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

        # 组学 × 组织 × 比较的 QC 可靠性表（omicsweave.core.reliability）；没有时各层不因 QC 降权
        # 同工酶加权：次要同工酶（丰度占比低于该值）单独支持、主要同工酶反向时不计为支持
        self.minor_share: float = float((scoring.get("isozyme") or {}).get("minor_share", 0.2))
        self.evidence_profile_scorer = EvidenceProfileScorer(
            scoring.get("evidence_profile", {}), reliability=getattr(config, "contrast_reliability", None))
        self._evidence_table = evidence_table([])

    # ---- 成员集合与入选规则（证据判定与反应档案共用） ----

    def gene_members(self, flux_id: Any) -> tuple[list[str], set[str]]:
        """转录组层的成员（全部映射基因）与方向成员。"""
        mapping = self.package.gene_mapping(flux_id)
        genes = list(dict.fromkeys(mapping.get("genes", [])))
        gene_roles = mapping.get("gene_roles", {}) or {}
        protein_mapping = self.package.enzyme_mapping(flux_id)
        protein_direct = list(protein_mapping.get("direct_enzymes", []))
        # 方向成员与蛋白层一致：直接酶（蛋白表优先，没有蛋白条目时用基因表的 catalytic）；
        # 只有没有直接酶的转运 / 交换反应，才由载体充当直接执行者。
        direct_genes = set(protein_direct) if protein_direct else set(gene_roles.get("catalytic", []))
        carriers = set(protein_mapping.get("transporters", [])) | set(gene_roles.get("transporter", []))
        if not direct_genes:
            direct_genes = set(carriers)
        genes = list(dict.fromkeys(genes + sorted(direct_genes) + sorted(carriers)))
        relations = self.package.member_relations(str(flux_id))
        # 多个复合体共用的亚基（如 DLD）的变化不只反映这个反应：只记为关联，不作为方向成员（G17）
        direct_genes -= {member for member, relation in relations.items() if relation == "shared_subunit"}
        if not self.package.semantics(str(flux_id)).directional:
            direct_genes = set()           # 只表示容量：匹配到的基因全部记为关联，不判支持或冲突
        return genes, direct_genes

    def protein_members(self, flux_id: Any) -> tuple[list[str], list[str]]:
        """蛋白组层的成员（全部映射蛋白）与方向成员。"""
        mapping = self.package.enzyme_mapping(flux_id)
        direct_enzymes: list[str] = list(mapping.get("direct_enzymes", []))
        directional_proteins = list(direct_enzymes)
        if not directional_proteins:
            # 运输/交换反应通常没有 EC 对应的直接酶，此时转运蛋白就是直接执行者。
            directional_proteins = list(mapping.get("transporters", []))
        relations = self.package.member_relations(str(flux_id))
        # 共用亚基只记为关联，不作为方向成员（G17，同转录组）
        directional_proteins = [p for p in dict.fromkeys(directional_proteins) if relations.get(p) != "shared_subunit"]
        if not self.package.semantics(str(flux_id)).directional:
            directional_proteins = []      # 只表示容量：匹配到的蛋白全部记为关联
        all_proteins: list[str] = list(direct_enzymes)
        for key in PROTEIN_CATEGORIES[1:]:
            all_proteins.extend(mapping.get(key, []))
        return list(dict.fromkeys(all_proteins)), directional_proteins

    def _tx_selection(self, sub: pd.DataFrame, log2fc_col: str, padj_col: str,
                      evidence_col: str | None) -> tuple[pd.Series, str]:
        """转录组的入选：padj < cutoff 且 |log2FC| > cutoff；一个也没有时退回效应量判定。返回 (掩码, 判定方式)。"""
        tx_cut = self.omics_cutoffs["transcriptomics"]
        # ── Phase 1: standard statistical significance ──
        sig_mask = pd.Series(True, index=sub.index)
        if padj_col in sub.columns:
            sig_mask &= sub[padj_col].fillna(1.0) < tx_cut["pval_cutoff"]
        if log2fc_col in sub.columns:
            sig_mask &= sub[log2fc_col].fillna(0.0).abs() > tx_cut["log2fc_cutoff"]
        # Exclude effect_size_only rows from statistical path
        if evidence_col and evidence_col in sub.columns:
            sig_mask &= sub[evidence_col].fillna("") != "effect_size_only"
        # ── Phase 2: effect-size-only fallback (N=1) ──
        if not sig_mask.any() and self.fallback_enabled and log2fc_col in sub.columns:
            # Use the large effect threshold as a stricter substitute for p-value
            fc_fallback_mask = sub[log2fc_col].fillna(0.0).abs() >= self.fallback_fc_large
            if evidence_col and evidence_col in sub.columns:
                # Include effect_size_only rows that were excluded from statistical path
                effect_only_mask = sub[evidence_col].fillna("") == "effect_size_only"
                fc_fallback_mask |= (effect_only_mask & sub[log2fc_col].fillna(0.0).abs().ge(self.fallback_fc_large))
            return fc_fallback_mask, "effect_size_only" if fc_fallback_mask.any() else "statistical_test"
        return sig_mask, "statistical_test"

    def _prot_selection(self, sub: pd.DataFrame, comparison: str) -> tuple[pd.Series, str]:
        """蛋白组的入选：FDR（没有 FDR 时用原始 p）< cutoff、经过检验且 |log2FC| > cutoff；没有有效 p 值时退回效应量判定。"""
        pval_col, log2fc_col = f"pval_{comparison}", f"log2fc_{comparison}"
        padj_col, evidence_col = f"padj_{comparison}", f"evidence_{comparison}"
        prot_cut = self.omics_cutoffs["proteomics"]
        # ── Phase 1: standard statistical significance ──
        sig_mask = pd.Series(True, index=sub.index)
        # FDR-adjusted p-values are the primary inferential contract. Raw
        # p-values are used only for audited legacy inputs that have no finite
        # adjusted values at all.
        has_valid_pval = False
        if padj_col in sub.columns and sub[padj_col].notna().any():
            has_valid_pval = True
            sig_mask &= sub[padj_col].fillna(1.0) < prot_cut["pval_cutoff"]
        elif pval_col in sub.columns and sub[pval_col].notna().any():
            has_valid_pval = True
            sig_mask &= sub[pval_col].fillna(1.0) < prot_cut["pval_cutoff"]
        if evidence_col in sub.columns:
            sig_mask &= sub[evidence_col].astype(str).eq("inferential")
        if log2fc_col in sub.columns:
            sig_mask &= sub[log2fc_col].fillna(0.0).abs() > prot_cut["log2fc_cutoff"]
        # ── Phase 2: effect-size-only fallback (no valid p-value) ──
        if not sig_mask.any() and self.fallback_enabled and not has_valid_pval and log2fc_col in sub.columns:
            fc_fallback_mask = sub[log2fc_col].fillna(0.0).abs() >= self.fallback_fc_large
            return fc_fallback_mask, "effect_size_only" if fc_fallback_mask.any() else "statistical_test"
        return sig_mask, "statistical_test"

    def _metab_selection(self, all_sub: pd.DataFrame, comparison: str) -> tuple[pd.Series, str]:
        """代谢组的入选：有有效 p 值时 FDR < cutoff 且 |log2FC| > cutoff；整组都没有 p 值时退回效应量判定。"""
        pval_col, log2fc_col, fdr_col = f"pval_{comparison}", f"log2fc_{comparison}", f"fdr_{comparison}"
        metab_cut = self.omics_cutoffs["metabolomics"]
        has_valid_pval = pval_col in all_sub.columns and all_sub[pval_col].notna().any()
        # ── Phase 1: standard significance ──
        mask = pd.Series(True, index=all_sub.index)
        if fdr_col in all_sub.columns:
            if has_valid_pval:
                mask &= all_sub[fdr_col].fillna(1.0) < metab_cut["pval_cutoff"]
            else:
                mask &= pd.Series(False, index=all_sub.index)
        elif pval_col in all_sub.columns:
            mask &= all_sub[pval_col].fillna(1.0) < metab_cut["pval_cutoff"]
        if log2fc_col in all_sub.columns:
            mask &= all_sub[log2fc_col].fillna(0.0).abs() > metab_cut["log2fc_cutoff"]
        # ── Phase 2: effect-size-only fallback ──
        if not mask.any() and self.fallback_enabled and not has_valid_pval:
            fallback = all_sub[log2fc_col].fillna(0.0).abs() >= self.fallback_fc_large \
                if log2fc_col in all_sub.columns else pd.Series(False, index=all_sub.index)
            return fallback, "effect_size_only" if fallback.any() else "statistical_test"
        return mask, "statistical_test"

    # ---- 反应档案的成员行（G13） ----

    def member_role(self, flux_id: Any, member: str, layer: str) -> tuple[str, str, str, str]:
        """成员在该层的 (角色, 成员关系, 作用方向, 作用方式)：成员知识表（G17）优先，其次蛋白表的分类，
        转录组再看基因表的角色；都没有时为 associated。"""
        annotated = self.package.member_annotation(str(flux_id))
        if annotated is not None and member in annotated.members:
            known = annotated.members[member]
            return known.role, known.relation, known.effect, known.mode
        protein = self.package.enzyme_mapping(flux_id)
        role = next((PROTEIN_CATEGORY_ROLE[c] for c in PROTEIN_CATEGORIES if member in (protein.get(c) or [])), "")
        if not role and layer == "transcriptomics":
            gene_roles = self.package.gene_mapping(flux_id).get("gene_roles", {}) or {}
            role = next((GENE_ROLE_ROLE.get(r, r) for r, values in gene_roles.items() if member in (values or [])), "")
        role = role or "associated"
        return role, "isozyme" if role in ("direct_enzyme", "carrier") else "", "", ""

    def member_matches(self, layer: str, flux_id: Any, frame: pd.DataFrame | None,
                       comparison: str = "") -> pd.DataFrame:
        """一个反应在一层中的成员行（反应档案，G13），列见 :data:`MEMBER_COLUMNS`。

        成员集合、匹配函数与入选规则（``selected``）都与 ``*_evidence`` 的判定相同，所以 ``in_flux_anchored`` 为真的行
        按层汇总后，与证据长表的映射数、匹配数、入选数一致。每个成员至少一行：没有匹配到时 ``row`` 为 None；
        匹配到多行时每行一条。``frame`` 为 None 时只列成员（该层没有数据）。代谢组另有派生特征（``layer`` 为
        ``derived``）、匹配到多个角色的特征（``ambiguous``），以及成员知识表中的别构效应物；后两者不计入证据判定。
        """
        if layer in ("transcriptomics", "proteomics"):
            return self._identifier_members(layer, flux_id, frame, comparison)
        if layer == "metabolomics":
            return self._metabolite_members(str(flux_id), frame, comparison)
        raise ValueError(f"member_matches 不支持的层：{layer}")

    def _identifier_members(self, layer: str, flux_id: Any, frame: pd.DataFrame | None,
                            comparison: str) -> pd.DataFrame:
        if layer == "transcriptomics":
            members, direct = self.gene_members(flux_id)
            usable = frame is not None and not frame.empty and any(c in frame.columns for c in _GENE_COLUMNS)
            match = match_gene_features(frame, members) if usable and members else IdentifierMatch((), (), {})
        else:
            members, directional = self.protein_members(flux_id)
            direct = set(directional)
            usable = frame is not None and not frame.empty and any(c in frame.columns for c in _PROTEIN_COLUMNS)
            match = match_protein_features(frame, members) if usable and members else IdentifierMatch((), (), {})
        selected: dict[Any, bool] = {}
        mode = ""
        if match.indices:
            sub = frame.loc[list(match.indices)]
            if layer == "transcriptomics":
                padj_col, log2fc_col, evidence_col = _tx_columns(sub)
                mask, mode = self._tx_selection(sub, log2fc_col, padj_col, evidence_col)
            elif f"log2fc_{comparison}" in frame.columns:
                mask, mode = self._prot_selection(sub, comparison)
            else:
                mask, mode = pd.Series(False, index=sub.index), "no_data"
            selected = {index: bool(value) for index, value in mask.items()}
        rows_by_query: dict[str, list[Any]] = defaultdict(list)
        for index in match.indices:
            rows_by_query[match.query_by_index[index]].append(index)
        records, seen = [], set()
        for member in members:
            canonical = _canonical_gene_symbol(member) or str(member)
            if canonical in seen:
                continue
            seen.add(canonical)
            role, relation, effect, regulator_mode = self.member_role(flux_id, str(member), layer)
            for index in rows_by_query.get(canonical) or [None]:
                records.append([layer, role, str(member), relation, effect, regulator_mode, member in direct, index,
                                match.match_method_by_index.get(index, "") if index is not None else "",
                                selected.get(index, False), mode if index is not None else "", True])
        return pd.DataFrame(records, columns=MEMBER_COLUMNS, dtype=object)   # object：行索引与 None 不被转成浮点

    def _metabolite_members(self, flux_id: str, frame: pd.DataFrame | None, comparison: str) -> pd.DataFrame:
        roles = _metabolite_roles(self.package.metabolite_roles(flux_id))
        usable = frame is not None and not frame.empty and "compound" in frame.columns
        derived_idx = _derived_indices(flux_id, frame, self.package.lipid_pools(flux_id)) if usable \
            else {"substrates": [], "products": [], "related": []}
        measured = measured_metabolites(frame) if usable else frame
        match = match_metabolite_roles(measured, roles) if usable else None
        hits = match.indices_by_role if match else {"substrates": (), "products": (), "related": ()}
        all_idx = list(dict.fromkeys([*hits["substrates"], *hits["products"], *hits["related"],
                                      *derived_idx["substrates"], *derived_idx["products"], *derived_idx["related"]]))
        has_effects = usable and f"log2fc_{comparison}" in frame.columns
        selected: dict[Any, bool] = {}
        mode = ""
        if all_idx and has_effects:
            mask, mode = self._metab_selection(frame.loc[all_idx], comparison)
            selected = {index: bool(value) for index, value in mask.items()}
        directional = self.package.semantics(flux_id).directional
        records: list[list[Any]] = []

        def add(layer: str, role: str, member_id: str, index: Any, direction_member: bool, counted: bool,
                method: str = "", effect: str = "", regulator_mode: str = "", chosen: dict | None = None,
                how: str = "") -> None:
            chosen = selected if chosen is None else chosen
            records.append([layer, role, member_id, "", effect, regulator_mode, direction_member, index,
                            method if index is not None else "", bool(chosen.get(index, False)),
                            (how or mode) if index is not None else "", counted])

        def method_of(index: Any) -> str:
            return match.match_method_by_index.get(index, "") if match else ""

        for key, role in (("substrates", "substrate"), ("products", "product"), ("related", "related_metabolite")):
            members = metabolite_members(roles[key])
            owned: dict[str, list[Any]] = defaultdict(list)
            for index in hits[key]:
                owned[metabolite_member_of(measured, index, members) or str(measured.at[index, "compound"])].append(index)
            for member in members:
                rows = owned.pop(member.member_id, [])
                if key == "related" and not rows:
                    continue                 # 关联代谢物只列测到的（知识表的关联条目含通路名等非代谢物）
                for index in rows or [None]:
                    add("metabolomics", role, member.member_id, index, directional and key != "related", True,
                        method_of(index))
            for owner, rows in owned.items():
                for index in rows:
                    add("metabolomics", role, owner, index, directional and key != "related", True, method_of(index))
        for index in (match.ambiguous_indices if match else ()):
            add("metabolomics", "ambiguous", str(measured.at[index, "compound"]), index, False, False,
                "ambiguous_role", chosen={})
        for key, indices in derived_idx.items():
            role = {"substrates": "substrate", "products": "product"}.get(key, "step_ratio")
            for index in indices:
                add("derived", role, str(frame.at[index, "compound"]), index, directional and key != "related", True,
                    "reaction_annotation")
        annotated = self.package.member_annotation(flux_id)
        for effector in (annotated.by_role("effector") if annotated else []):
            found = match_metabolite_roles(measured, {"related": [effector.id]}) if usable else None
            rows = list(found.indices_by_role["related"]) if found else []
            chosen, how = {}, ""
            if rows and has_effects:
                mask, how = self._metab_selection(frame.loc[rows], comparison)
                chosen = {index: bool(value) for index, value in mask.items()}
            for index in rows or [None]:
                add("metabolomics", "effector", effector.id, index, False, False,
                    found.match_method_by_index.get(index, "") if found else "", effector.effect, effector.mode,
                    chosen, how)
        return pd.DataFrame(records, columns=MEMBER_COLUMNS, dtype=object)

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
        genes, direct_genes = self.gene_members(flux_id)
        relations = self.package.member_relations(str(flux_id))
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

        padj_col, log2fc_col, evidence_col = _tx_columns(sub)
        sig_mask, tx_evidence_mode = self._tx_selection(sub, log2fc_col, padj_col, evidence_col)
        sig_df = sub[sig_mask]
        n_sig = len(sig_df)

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
        implied = _implied_counts(direct_sig[log2fc_col] if log2fc_col in direct_sig.columns
                                  else pd.Series(dtype=float))
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
                **implied,
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
            **implied,
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
        excluded = _non_isozyme_rows(tx_df, relations, proteins=False)
        direct_rows = sub.loc[[index for index in sub.index if index in direct_indices and index not in excluded]]
        details = _demote_minor_support(details, _isozyme_check(
            direct_rows, log2fc_col, "base_mean", gene_column, flux_fc, concordant_genes_list, self.minor_share),
            "concordant_genes")
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
        all_proteins, directional_proteins = self.protein_members(flux_id)
        relations = self.package.member_relations(str(flux_id))
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

        sig_mask, prot_evidence_mode = self._prot_selection(sub, comparison)
        sig_df = sub[sig_mask]
        n_sig = len(sig_df)

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
            lfsr_candidates=[f"lfsr_{comparison}"],
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
        implied = _implied_counts(direct_sig[log2fc_col] if log2fc_col in direct_sig.columns
                                  else pd.Series(dtype=float))
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
                **implied,
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
            **implied,
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
        excluded = _non_isozyme_rows(prot_df, relations, proteins=True)
        direct_rows = sub.loc[[index for index in sub.index if index in direct_indices and index not in excluded]]
        details = _demote_minor_support(details, _isozyme_check(
            direct_rows, log2fc_col, "abundance", protein_column, flux_fc, concordant_proteins_list, self.minor_share),
            "concordant_proteins")
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
        derived_idx = _derived_indices(str(flux_id), metab_df, self.package.lipid_pools(str(flux_id)))
        has_derived = any(derived_idx.values())

        roles = self.package.metabolite_roles(flux_id)
        if roles is None or isinstance(roles, list):
            flat = roles if roles is not None else []
            if isinstance(flat, list):
                roles = {"substrates": [], "products": [], "related": flat}
            elif isinstance(flat, dict):
                roles = flat
            elif has_derived:
                roles = {"substrates": [], "products": [], "related": []}
            else:
                return {"n_metabolites_checked": 0, "n_sig": 0, "sig_metabolites": [], "concordant_metabolites": [], "evidence_mode": "no_mapping"}

        substrates = roles.get("substrates", [])
        products = roles.get("products", [])
        related = roles.get("related", [])
        all_metabs = substrates + products + related
        n_checked = len(all_metabs)

        if not all_metabs and not has_derived:
            return {"n_metabolites_checked": 0, "n_sig": 0, "sig_metabolites": [], "concordant_metabolites": [], "evidence_mode": "no_mapping"}

        if "compound" not in metab_df.columns:
            return {"n_metabolites_checked": n_checked, "n_sig": 0, "sig_metabolites": [], "concordant_metabolites": [], "evidence_mode": "no_data"}

        # 派生特征不参与按名称 / ID 的匹配，只按反应注释挂到底物 / 产物 / 关联上
        measured = measured_metabolites(metab_df)
        matched_roles = match_metabolite_roles(measured, roles)
        n_checked = sum(matched_roles.target_count_by_role.values()) + sum(map(len, derived_idx.values()))
        sub_idx = list(matched_roles.indices_by_role["substrates"]) + derived_idx["substrates"]
        prod_idx = list(matched_roles.indices_by_role["products"]) + derived_idx["products"]
        rel_idx = list(matched_roles.indices_by_role["related"]) + derived_idx["related"]
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

        all_sub = metab_df.loc[all_idx]
        sig_mask, metab_evidence_mode = self._metab_selection(all_sub, comparison)
        sig_df = all_sub[sig_mask]
        n_sig = len(sig_df)

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
            directional = self.package.semantics(str(flux_id)).directional
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

        # 不依赖通量符号的方向提示：产物升 → 通量升，底物升 → 通量降（只表示容量的反应不提示方向）
        implied_signs: list[float] = []
        if n_sig > 0 and log2fc_col in sig_df.columns and self.package.semantics(str(flux_id)).directional:
            for row_idx, row_fc in pd.to_numeric(sig_df[log2fc_col], errors="coerce").dropna().items():
                if row_idx in set(prod_idx):
                    implied_signs.append(float(np.sign(row_fc)))
                elif row_idx in set(sub_idx):
                    implied_signs.append(-float(np.sign(row_fc)))

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
            lfsr_candidates=[f"lfsr_{comparison}"],
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
            **_implied_counts(pd.Series(implied_signs, dtype=float)),
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
            "derived_readouts": [f"{metab_df.at[index, 'compound']} ({role})"
                                 for role, indices in derived_idx.items() for index in indices],
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
        semantics = self.package.semantics(str(flux_row.get("flux_id", "")))
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
            # 锚点状态的依据（统一差异结果表中的 MFA 证据等级、较小组的样本量、可辨识性与示踪剂一致性）
            "mfa_evidence_level": flux_row.get("mfa_evidence_level", ""),
            "n_min_group": _series_numeric_value(flux_row, ["n_min_group"]),
            "identifiability": flux_row.get("identifiability", ""),
            "mfa_flags": flux_row.get("mfa_flags", ""),
            "tracer_evidence": flux_row.get("tracer_evidence", ""),
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
                "n_anchor_established": int(g["anchor_state"].isin(ANCHORED_STATES).sum())
                if "anchor_state" in g.columns else 0,
                "n_omics_concordant": int((g["evidence_grade"] == "Omics-concordant").sum()),
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
