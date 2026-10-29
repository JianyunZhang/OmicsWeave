#!/usr/bin/env python3
"""多组学连续证据画像。

该模块只接受标准化后的层级证据，负责把统计可信度、效应强度、映射质量、
数据可靠性和方向关系整合为可解释的多维结果。旧的 0--4 分矩阵不属于本模块。
"""

from __future__ import annotations

from dataclasses import dataclass
from math import erf, sqrt
from typing import Any, Iterable

import numpy as np
import pandas as pd


LAYERS = ("mfa", "transcriptomics", "proteomics", "metabolomics")
ASSESSABLE_STATUSES = {
    "no_signal",
    "support",
    "opposite",
    "mixed",
    "unresolved_signal",
}

_COUNT_KEYS = {
    "mfa": ("n_reactions_checked", "n_reactions_found"),
    "transcriptomics": ("n_genes_checked", "n_genes_found"),
    "proteomics": ("n_proteins_checked", "n_proteins_found"),
    "metabolomics": ("n_metabolites_checked", "n_metabolites_found"),
}

_CONCORDANT_KEYS = {
    "transcriptomics": "concordant_genes",
    "proteomics": "concordant_proteins",
    "metabolomics": "concordant_metabolites",
}


def _as_float(value: Any) -> float:
    """安全转换为有限浮点数；无法转换时返回 NaN。"""
    try:
        result = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return result if np.isfinite(result) else float("nan")


def _format_mapping_methods(value: Any) -> str:
    if not isinstance(value, dict):
        return ""
    return ";".join(
        f"{method}:{int(count)}"
        for method, count in sorted(value.items())
        if int(count) > 0
    )


def _format_ambiguous_features(details: dict[str, Any]) -> str:
    for key in (
        "ambiguous_metabolites",
        "ambiguous_proteins",
        "ambiguous_genes",
    ):
        values = details.get(key)
        if values:
            return ";".join(sorted({str(value) for value in values}))
    return ""


def _clip01(value: float) -> float:
    if not np.isfinite(value):
        return float("nan")
    return float(np.clip(value, 0.0, 1.0))


def _normal_sign_certainty(effect: float, se: float) -> float:
    """由 effect/SE 计算方向可信度代理，范围 0--1。"""
    if not np.isfinite(effect) or not np.isfinite(se) or se <= 0:
        return float("nan")
    z_value = abs(effect / se)
    # 2 * Phi(|z|) - 1；z=0 时为 0，随证据增强趋近 1。
    return _clip01(erf(z_value / sqrt(2.0)))


@dataclass(frozen=True)
class EvidenceProfileConfig:
    """连续证据评分的稳定外部接口。"""

    support_threshold: float
    conflict_threshold: float
    high_conflict_threshold: float
    effect_scales: dict[str, float]
    mapping_quality: dict[str, float]
    mode_reliability: dict[str, float]
    quality_reliability: dict[str, float]
    n_eff_by_layer: dict[str, float]
    mfa_inference_unit: str
    design_grade_override: str | None

    @classmethod
    def from_dict(cls, config: dict[str, Any] | None) -> "EvidenceProfileConfig":
        raw = dict(config or {})
        return cls(
            support_threshold=float(raw.get("support_threshold", 0.15)),
            conflict_threshold=float(raw.get("conflict_threshold", 0.15)),
            high_conflict_threshold=float(raw.get("high_conflict_threshold", 0.25)),
            effect_scales={
                "mfa": 1.0,
                "transcriptomics": 1.0,
                "proteomics": 1.0,
                "metabolomics": 1.0,
                **raw.get("effect_scales", {}),
            },
            mapping_quality={
                "mfa": 1.0,
                "transcriptomics": 0.90,
                "transcriptomics_direct": 1.0,
                "transcriptomics_indirect": 0.65,
                "proteomics_direct": 1.0,
                "proteomics_indirect": 0.65,
                "metabolomics_direct": 0.90,
                "metabolomics_related": 0.50,
                **raw.get("mapping_quality", {}),
            },
            mode_reliability={
                "biological_replicate": 1.0,
                "statistical_test": 0.80,
                "mixed": 0.60,
                "effect_size_only": 0.35,
                "pseudo_replicate": 0.30,
                "descriptive": 0.25,
                **raw.get("mode_reliability", {}),
            },
            quality_reliability={
                "mfa": 1.0,
                "transcriptomics": 1.0,
                "proteomics": 1.0,
                "metabolomics": 1.0,
                **raw.get("quality_reliability", {}),
            },
            n_eff_by_layer={
                str(key): float(value)
                for key, value in raw.get("n_eff_by_layer", {}).items()
            },
            # 当前项目 MFA p 值来自 MCMC 拟合点。上游改为生物学重复推断后，
            # 应在配置中显式设置 biological_replicate。
            mfa_inference_unit=str(raw.get("mfa_inference_unit", "mcmc_fit_points")),
            design_grade_override=raw.get("design_grade"),
        )


class EvidenceProfileScorer:
    """把四层标准化证据整合为 EvidenceTable 和多维证据画像。"""

    def __init__(self, config: dict[str, Any] | None = None) -> None:
        self.config = EvidenceProfileConfig.from_dict(config)

    def build_records(self, scored_record: dict[str, Any]) -> list[dict[str, Any]]:
        """为单个反应构建四行标准化 EvidenceTable 记录。"""
        return [self._build_layer_record(scored_record, layer) for layer in LAYERS]

    def summarize(self, records: Iterable[dict[str, Any]]) -> dict[str, Any]:
        """汇总单个反应的支持、冲突、覆盖、稳定性和证据等级。"""
        table = pd.DataFrame(list(records))
        if table.empty:
            return self._empty_profile()

        assessable = table[table["assessment_status"].isin(ASSESSABLE_STATUSES)]
        n_assessable = len(assessable)
        coverage_score = float(table["coverage_value"].mean())

        if n_assessable:
            support_score = float(assessable["support_evidence"].mean())
            conflict_score = float(assessable["conflict_evidence"].mean())
        else:
            support_score = 0.0
            conflict_score = 0.0

        n_supporting = int(
            (assessable["support_evidence"] >= self.config.support_threshold).sum()
        )
        n_conflicting = int(
            (assessable["conflict_evidence"] >= self.config.conflict_threshold).sum()
        )

        selection_probability = self._nanmean(table["selection_probability"])
        direction_stability = self._nanmean(table["direction_stability"])
        design_grade = self._design_grade(
            assessable,
            selection_probability=selection_probability,
            direction_stability=direction_stability,
        )
        evidence_grade = self._evidence_grade(
            n_assessable=n_assessable,
            n_supporting=n_supporting,
            n_conflicting=n_conflicting,
            support_score=support_score,
            conflict_score=conflict_score,
            design_grade=design_grade,
        )

        return {
            "support_score": round(_clip01(support_score), 4),
            "conflict_score": round(_clip01(conflict_score), 4),
            "coverage_score": round(_clip01(coverage_score), 4),
            "stability_selection_probability": self._round_or_nan(selection_probability),
            "stability_direction": self._round_or_nan(direction_stability),
            "design_grade": design_grade,
            "evidence_grade": evidence_grade,
            "n_assessable_layers": n_assessable,
            "n_supporting_layers": n_supporting,
            "n_conflicting_layers": n_conflicting,
        }

    def _build_layer_record(
        self,
        scored_record: dict[str, Any],
        layer: str,
    ) -> dict[str, Any]:
        details = dict(scored_record.get(f"{self._detail_prefix(layer)}_details", {}) or {})
        mode = str(details.get("evidence_mode", scored_record.get(f"{self._mode_prefix(layer)}_evidence_mode", "no_data")))

        if layer == "mfa":
            return self._build_mfa_record(scored_record, details, mode)

        checked_key, found_key = _COUNT_KEYS[layer]
        n_mapped = int(details.get(checked_key, 0) or 0)
        default_found = n_mapped if mode not in {"no_data", "no_mapping", "no_match"} else 0
        n_matched = int(details.get(found_key, default_found) or 0)
        n_signal = int(details.get("n_sig", 0) or 0)
        concordant_key = _CONCORDANT_KEYS[layer]
        n_support = int(
            details.get("n_concordant", len(details.get(concordant_key, []) or [])) or 0
        )
        if "n_opposite" in details:
            n_conflict = int(details.get("n_opposite", 0) or 0)
        elif details.get("direction_override") == "opposite":
            n_conflict = n_signal
        else:
            n_conflict = max(n_signal - n_support - int(details.get("n_unresolved", 0) or 0), 0)

        status = self._assessment_status(mode, n_signal, n_support, n_conflict)
        support_effect = abs(_as_float(details.get("support_effect")))
        conflict_effect = abs(_as_float(details.get("conflict_effect")))
        representative_effect = _as_float(details.get("representative_effect"))

        if not np.isfinite(support_effect) and n_support:
            support_effect = abs(representative_effect)
        if not np.isfinite(conflict_effect) and n_conflict:
            conflict_effect = abs(representative_effect)

        scale = max(float(self.config.effect_scales.get(layer, 1.0)), np.finfo(float).eps)
        support_strength = self._effect_strength(support_effect, scale)
        conflict_strength = self._effect_strength(conflict_effect, scale)

        qvalue = _as_float(details.get("qvalue"))
        se = _as_float(details.get("se"))
        n_eff = self._n_eff(layer, details.get("n_eff"))
        mapping_quality = self._mapping_quality(layer, details)
        quality_reliability = self._quality_reliability(layer)
        reliability = self._data_reliability(
            mode, n_eff, qvalue, se, quality_reliability
        )
        sign_certainty = self._sign_certainty(
            mode=mode,
            effect=max(support_effect, conflict_effect),
            se=se,
            qvalue=qvalue,
            lfsr=_as_float(details.get("lfsr")),
        )

        signal_denominator = max(n_signal, 1)
        support_fraction = n_support / signal_denominator
        conflict_fraction = n_conflict / signal_denominator
        prevalence = sqrt(n_signal / max(n_matched, 1)) if n_signal else 0.0

        support_evidence = (
            support_fraction
            * prevalence
            * support_strength
            * sign_certainty
            * mapping_quality
            * reliability
        )
        conflict_evidence = (
            conflict_fraction
            * prevalence
            * conflict_strength
            * sign_certainty
            * mapping_quality
            * reliability
        )

        if status == "mixed":
            relative_effect = support_effect - conflict_effect
            direction_relation = "mixed"
        elif status == "opposite":
            relative_effect = -conflict_effect
            direction_relation = "opposite"
        elif status == "support":
            relative_effect = support_effect
            direction_relation = "support"
        else:
            relative_effect = representative_effect
            direction_relation = "unresolved" if n_signal else "none"

        if status == "mixed":
            direction = "mixed"
        elif representative_effect > 0:
            direction = "up"
        elif representative_effect < 0:
            direction = "down"
        else:
            direction = "unknown" if n_signal else "none"

        return {
            **self._identity(scored_record),
            **self._provenance(details),
            "layer": layer,
            "effect": representative_effect,
            "mechanism_aligned_effect": relative_effect,
            "se": se,
            "qvalue": qvalue,
            "n_eff": n_eff,
            "direction": direction,
            "direction_relation": direction_relation,
            "mapping_quality": mapping_quality,
            "mapping_methods": _format_mapping_methods(
                details.get("mapping_method_counts")
            ),
            "n_ambiguous": int(details.get("n_ambiguous", 0) or 0),
            "ambiguous_features": _format_ambiguous_features(details),
            "data_reason": str(details.get("data_reason", "") or ""),
            "quality_reliability": quality_reliability,
            "data_reliability": reliability,
            "evidence_mode": mode,
            "assessment_status": status,
            "sign_certainty": sign_certainty,
            "effect_strength": max(support_strength, conflict_strength),
            "support_evidence": _clip01(support_evidence),
            "conflict_evidence": _clip01(conflict_evidence),
            "coverage_value": 1.0 if status in ASSESSABLE_STATUSES else 0.0,
            "n_mapped": n_mapped,
            "n_matched": n_matched,
            "n_signal": n_signal,
            "n_support": n_support,
            "n_conflict": n_conflict,
            "selection_probability": _as_float(details.get("selection_probability")),
            "direction_stability": _as_float(details.get("direction_stability")),
        }

    def _build_mfa_record(
        self,
        scored_record: dict[str, Any],
        details: dict[str, Any],
        mode: str,
    ) -> dict[str, Any]:
        gate_passed = bool(scored_record.get("gate_passed", False))
        effect = _as_float(scored_record.get("flux_log2fc"))
        se = _as_float(details.get("se", scored_record.get("flux_se")))
        qvalue = _as_float(scored_record.get("flux_fdr"))
        n_eff = self._n_eff("mfa", details.get("n_eff", scored_record.get("n_eff")))

        if self.config.mfa_inference_unit != "biological_replicate":
            mode = "pseudo_replicate" if np.isfinite(qvalue) else "descriptive"
        if np.isfinite(n_eff) and n_eff < 3:
            mode = "descriptive"

        status = "support" if gate_passed else "no_signal"
        n_signal = 1 if gate_passed else 0
        effect_strength = self._effect_strength(
            abs(effect),
            max(float(self.config.effect_scales.get("mfa", 1.0)), np.finfo(float).eps),
        )
        quality_reliability = self._quality_reliability("mfa")
        reliability = self._data_reliability(
            mode, n_eff, qvalue, se, quality_reliability
        )
        sign_certainty = self._sign_certainty(
            mode=mode,
            effect=effect,
            se=se,
            qvalue=qvalue,
            lfsr=_as_float(details.get("lfsr")),
        )
        support_evidence = (
            effect_strength * reliability * sign_certainty if gate_passed else 0.0
        )

        return {
            **self._identity(scored_record),
            **self._provenance(details, fallback=scored_record),
            "layer": "mfa",
            "effect": effect,
            "mechanism_aligned_effect": effect,
            "se": se,
            "qvalue": qvalue,
            "n_eff": n_eff,
            "direction": "up" if effect > 0 else "down" if effect < 0 else "none",
            "direction_relation": "anchor" if gate_passed else "no_signal",
            "mapping_quality": 1.0,
            "quality_reliability": quality_reliability,
            "data_reliability": reliability,
            "evidence_mode": mode,
            "assessment_status": status,
            "sign_certainty": sign_certainty,
            "effect_strength": effect_strength,
            "support_evidence": _clip01(support_evidence),
            "conflict_evidence": 0.0,
            "coverage_value": 1.0,
            "n_mapped": 1,
            "n_matched": 1,
            "n_signal": n_signal,
            "n_support": n_signal,
            "n_conflict": 0,
            "selection_probability": _as_float(details.get("selection_probability")),
            "direction_stability": _as_float(details.get("direction_stability")),
        }

    @staticmethod
    def _detail_prefix(layer: str) -> str:
        return {
            "mfa": "mfa",
            "transcriptomics": "tx",
            "proteomics": "prot",
            "metabolomics": "metab",
        }[layer]

    @staticmethod
    def _mode_prefix(layer: str) -> str:
        return EvidenceProfileScorer._detail_prefix(layer)

    @staticmethod
    def _identity(scored_record: dict[str, Any]) -> dict[str, Any]:
        return {
            "experiment_id": scored_record.get("experiment_id", ""),
            "estimand_id": scored_record.get("estimand_id", ""),
            "flux_id": scored_record.get("flux_id", ""),
            "tissue": scored_record.get("tissue", ""),
            "comparison": scored_record.get("comparison", ""),
        }

    @staticmethod
    def _provenance(
        details: dict[str, Any],
        fallback: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        source = fallback or {}
        return {
            "protocol_id": details.get("protocol_id", source.get("protocol_id", "")),
            "batch_ids": details.get("batch_ids", source.get("batch_ids", "")),
            "n_batches": _as_float(details.get("n_batches", source.get("n_batches"))),
            "source_mode": details.get("source_mode", source.get("source_mode", "single_run")),
            "heterogeneity_i2": _as_float(
                details.get("heterogeneity_i2", source.get("heterogeneity_i2"))
            ),
            "heterogeneity_tau2": _as_float(
                details.get("heterogeneity_tau2", source.get("heterogeneity_tau2"))
            ),
        }

    @staticmethod
    def _assessment_status(
        mode: str,
        n_signal: int,
        n_support: int,
        n_conflict: int,
    ) -> str:
        aliases = {
            "no_significant": "no_signal",
            "insufficient": "no_signal",
        }
        if mode in {"no_data", "no_mapping", "no_match"}:
            return mode
        if n_signal <= 0:
            return aliases.get(mode, "no_signal")
        if n_support > 0 and n_conflict > 0:
            return "mixed"
        if n_support > 0:
            return "support"
        if n_conflict > 0:
            return "opposite"
        return "unresolved_signal"

    def _mapping_quality(self, layer: str, details: dict[str, Any]) -> float:
        explicit = _as_float(details.get("mapping_quality"))
        if np.isfinite(explicit):
            return _clip01(explicit)
        if layer == "proteomics":
            source = details.get("direction_source", "indirect")
            key = "proteomics_direct" if source == "direct_enzyme" else "proteomics_indirect"
            return _clip01(float(self.config.mapping_quality[key]))
        if layer == "transcriptomics":
            source = details.get("direction_source", "indirect_gene")
            key = (
                "transcriptomics_direct"
                if source == "direct_gene" else "transcriptomics_indirect"
            )
            return _clip01(float(self.config.mapping_quality[key]))
        if layer == "metabolomics":
            source = details.get("direction_source", "related")
            key = "metabolomics_direct" if source in {"substrate", "product", "direct"} else "metabolomics_related"
            return _clip01(float(self.config.mapping_quality[key]))
        return _clip01(float(self.config.mapping_quality.get(layer, 1.0)))

    def _data_reliability(
        self,
        mode: str,
        n_eff: float,
        qvalue: float,
        se: float,
        quality_reliability: float = 1.0,
    ) -> float:
        base = float(self.config.mode_reliability.get(mode, 0.50))
        if np.isfinite(n_eff):
            if n_eff < 2:
                n_factor = 0.20
            elif n_eff < 3:
                n_factor = 0.35
            elif n_eff < 6:
                n_factor = 0.75
            else:
                n_factor = 1.0
        else:
            n_factor = 0.70
        uncertainty_factor = 1.0 if np.isfinite(se) or np.isfinite(qvalue) else 0.75
        return _clip01(
            base * n_factor * uncertainty_factor * _clip01(quality_reliability)
        )

    def _quality_reliability(self, layer: str) -> float:
        """返回审核后显式配置的层级 QC 因子；默认不降权。"""
        value = _as_float(self.config.quality_reliability.get(layer, 1.0))
        return _clip01(value) if np.isfinite(value) else 1.0

    @staticmethod
    def _sign_certainty(
        mode: str,
        effect: float,
        se: float,
        qvalue: float,
        lfsr: float,
    ) -> float:
        if np.isfinite(lfsr):
            return _clip01(1.0 - lfsr)

        z_certainty = _normal_sign_certainty(effect, se)
        if np.isfinite(z_certainty):
            return z_certainty

        if mode in {"biological_replicate", "statistical_test", "mixed"} and np.isfinite(qvalue):
            return _clip01(1.0 - qvalue)
        if mode in {"effect_size_only", "pseudo_replicate", "descriptive"}:
            return 0.50
        return 0.50

    @staticmethod
    def _effect_strength(effect: float, scale: float) -> float:
        if not np.isfinite(effect):
            return 0.0
        return _clip01(abs(effect) / scale)

    def _n_eff(self, layer: str, value: Any) -> float:
        parsed = _as_float(value)
        if np.isfinite(parsed):
            return parsed
        return _as_float(self.config.n_eff_by_layer.get(layer))

    def _design_grade(
        self,
        assessable: pd.DataFrame,
        selection_probability: float,
        direction_stability: float,
    ) -> str:
        override = self.config.design_grade_override
        if override in {"A", "B", "C", "D"}:
            return str(override)
        if assessable.empty:
            return "D"

        exploratory_modes = {"effect_size_only", "pseudo_replicate", "descriptive"}
        if assessable["evidence_mode"].isin(exploratory_modes).any():
            return "D"

        known_n = assessable["n_eff"].dropna()
        if not known_n.empty and (known_n < 3).any():
            return "D"
        if len(known_n) < len(assessable) or assessable["se"].isna().any():
            return "C"
        if (
            (known_n >= 6).all()
            and np.isfinite(selection_probability)
            and np.isfinite(direction_stability)
            and selection_probability >= 0.80
            and direction_stability >= 0.80
        ):
            return "A"
        return "B"

    def _evidence_grade(
        self,
        n_assessable: int,
        n_supporting: int,
        n_conflicting: int,
        support_score: float,
        conflict_score: float,
        design_grade: str,
    ) -> str:
        if n_assessable == 0:
            return "Not-assessable"
        if (
            n_supporting > 0
            and n_conflicting > 0
            and support_score >= self.config.high_conflict_threshold
            and conflict_score >= self.config.high_conflict_threshold
        ):
            return "Conflicting"
        if n_supporting == 0:
            return "No-signal"
        if design_grade == "D":
            return "Exploratory"
        if n_supporting >= 2 and n_conflicting == 0:
            return "Convergent"
        if n_supporting == 1 and n_conflicting == 0:
            return "Single-layer"
        if n_conflicting > 0:
            return "Conflicting"
        return "Exploratory"

    @staticmethod
    def _nanmean(series: pd.Series) -> float:
        values = pd.to_numeric(series, errors="coerce").dropna()
        return float(values.mean()) if not values.empty else float("nan")

    @staticmethod
    def _round_or_nan(value: float) -> float:
        return round(float(value), 4) if np.isfinite(value) else float("nan")

    @staticmethod
    def _empty_profile() -> dict[str, Any]:
        return {
            "support_score": 0.0,
            "conflict_score": 0.0,
            "coverage_score": 0.0,
            "stability_selection_probability": float("nan"),
            "stability_direction": float("nan"),
            "design_grade": "D",
            "evidence_grade": "Not-assessable",
            "n_assessable_layers": 0,
            "n_supporting_layers": 0,
            "n_conflicting_layers": 0,
        }


def evidence_table(records: Iterable[dict[str, Any]]) -> pd.DataFrame:
    """把证据记录转换为具有稳定列顺序的 DataFrame。"""
    columns = [
        "experiment_id",
        "estimand_id",
        "flux_id",
        "tissue",
        "comparison",
        "layer",
        "protocol_id",
        "batch_ids",
        "n_batches",
        "source_mode",
        "heterogeneity_i2",
        "heterogeneity_tau2",
        "effect",
        "mechanism_aligned_effect",
        "se",
        "qvalue",
        "n_eff",
        "direction",
        "direction_relation",
        "mapping_quality",
        "mapping_methods",
        "n_ambiguous",
        "ambiguous_features",
        "data_reason",
        "quality_reliability",
        "data_reliability",
        "evidence_mode",
        "assessment_status",
        "sign_certainty",
        "effect_strength",
        "support_evidence",
        "conflict_evidence",
        "coverage_value",
        "n_mapped",
        "n_matched",
        "n_signal",
        "n_support",
        "n_conflict",
        "selection_probability",
        "direction_stability",
    ]
    table = pd.DataFrame(list(records))
    if table.empty:
        return pd.DataFrame(columns=columns)
    for column in columns:
        if column not in table.columns:
            table[column] = np.nan
    return table[columns]
