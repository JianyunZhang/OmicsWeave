#!/usr/bin/env python3
"""
Pathway-centric cross-tissue multi-omics integration scorer.

Unlike the flux-anchored evidence builder (`flux_evidence.py`, which needs
differential flux as the anchor), this scorer is designed for experiments with only
transcriptomics + metabolomics data (no MFA / proteomics). It scores
SYSTEMIC PATHWAYS that recur across multiple tissues by quantifying:

  - tissue breadth          (how many of N tissues show evidence)
  - cross-omics concordance (transcriptomics + metabolomics same direction)
  - effect magnitude        (representative |log2FC|)
  - cross-tissue novelty    (recurrence across >=2 tissues)

The former 0-4.0 composite is retained as a legacy compatibility result only.
Primary output is a curated-evidence profile with support/conflict/coverage,
robust median effect, recurrence, design grade, and evidence grade.

Input: curated pathway evidence (YAML) — see cross_tissue_evidence.yaml schema.
Output: pathway score table + tissue x pathway evidence matrix.

Author: Claude Code
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd


# All tissues the integration considers
DEFAULT_TISSUES = ["heart", "liver", "renal_cortex", "sq_adipose", "gast_muscle"]

# Significance / effect thresholds (same defaults as the flux-anchored evidence builder)
DEFAULT_LOG2FC_CUTOFF = 0.58
DEFAULT_LARGE_LOG2FC = 2.0

# Confidence tiers of this curated cross-tissue scorer (the flux-anchored additive score was removed in phase 2.8)
DEFAULT_CONFIDENCE = {"high": 3.0, "medium": 2.0, "low": 1.0}

# Normalize raw comparison keys (as written in evidence.yaml) to friendly
# "对比组别全称" labels for report rendering. Avoids vague terms like
# mixed/both by listing the actual comparison source(s) per pathway.
_COMPARISON_LABELS = {
    "HFD_Mazdu_vs_HFD_Saline": "Mazdutide_vs_Saline",
    "HFD_Sema_vs_HFD_Saline": "Semaglutide_vs_Saline",
    "HFD_Mazdu_vs_HFD_Sema": "Mazdutide_vs_Semaglutide",
    "HFD_Mazdu_vs_LFD_Mazdu": "Mazdutide_diet(HFD_vs_LFD)",
    "HFD_Sema_vs_LFD_Sema": "Semaglutide_diet(HFD_vs_LFD)",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _confidence_label(score: float, thresholds: dict[str, float]) -> str:
    """Map numeric score to confidence label."""
    if score >= thresholds.get("high", 3.0):
        return "High"
    if score >= thresholds.get("medium", 2.0):
        return "Medium"
    if score >= thresholds.get("low", 1.0):
        return "Low"
    return "Insufficient"


def _sign(x: float | None) -> int:
    """Return -1 / 0 / 1 for a numeric value (None -> 0)."""
    if x is None or pd.isna(x):
        return 0
    if x > 0:
        return 1
    if x < 0:
        return -1
    return 0


# ---------------------------------------------------------------------------
# Scorer
# ---------------------------------------------------------------------------

class CrossTissuePathwayScorer:
    """Score systemic pathways by cross-tissue + cross-omics convergence.

    Parameters
    ----------
    tissues : list[str]
        Tissues considered (default: 5 metabolic tissues).
    log2fc_cutoff : float
        Minimum |log2FC| to count as evidence (default 0.58).
    large_log2fc : float
        |log2FC| threshold for top effect-magnitude tier (default 2.0).
    confidence_thresholds : dict
        High/Medium/Low score boundaries.
    """

    def __init__(
        self,
        tissues: list[str] | None = None,
        *,
        log2fc_cutoff: float = DEFAULT_LOG2FC_CUTOFF,
        large_log2fc: float = DEFAULT_LARGE_LOG2FC,
        confidence_thresholds: dict[str, float] | None = None,
    ) -> None:
        self.tissues = list(tissues) if tissues else list(DEFAULT_TISSUES)
        self.log2fc_cutoff = log2fc_cutoff
        self.large_log2fc = large_log2fc
        self.confidence_thresholds = dict(confidence_thresholds or DEFAULT_CONFIDENCE)

    # ------------------------------------------------------------------
    # Per-dimension scorers
    # ------------------------------------------------------------------

    def _tissue_breadth(self, evidence: dict[str, Any]) -> tuple[float, list[str]]:
        """Score 0-1.5 by number of tissues with any significant evidence."""
        present = []
        for tissue, ev in evidence.items():
            if not isinstance(ev, dict):
                continue
            if self._has_significant_evidence(ev):
                present.append(tissue)
        n = len(present)
        # 0 -> 0, 1 -> 0.3, 2 -> 0.6, 3 -> 0.9, 4 -> 1.2, 5+ -> 1.5
        score = min(n, 5) * 0.3
        return score, present

    @staticmethod
    def _has_significant_evidence(ev: dict[str, Any]) -> bool:
        """A tissue counts as evidence if it has a direction or any omics entry."""
        if ev.get("direction") in ("UP", "DOWN"):
            return True
        for key in ("transcriptomics", "metabolomics"):
            layer = ev.get(key)
            if isinstance(layer, dict) and layer:
                return True
        return False

    def _omics_concordance(self, evidence: dict[str, Any]) -> tuple[float, int]:
        """Score 0-1.5 by tissues where BOTH omics agree in direction.

        Concordance = transcriptomics and metabolomics present AND same sign.
        """
        concordant_tissues = 0
        single_omics_tissues = 0
        for tissue, ev in evidence.items():
            if not isinstance(ev, dict):
                continue
            tx = ev.get("transcriptomics")
            mb = ev.get("metabolomics")
            tx_fc = tx.get("log2fc") if isinstance(tx, dict) else None
            mb_fc = mb.get("log2fc") if isinstance(mb, dict) else None

            has_tx = isinstance(tx, dict) and tx and abs(tx_fc or 0) >= self.log2fc_cutoff
            has_mb = isinstance(mb, dict) and mb and abs(mb_fc or 0) >= self.log2fc_cutoff

            if has_tx and has_mb:
                # Both present: concordant if same sign
                if _sign(tx_fc) == _sign(mb_fc) and _sign(tx_fc) != 0:
                    concordant_tissues += 1
                else:
                    single_omics_tissues += 1  # present but discordant
            elif has_tx or has_mb:
                single_omics_tissues += 1

        if concordant_tissues >= 2:
            score = 1.5
        elif concordant_tissues >= 1:
            score = 1.0
        elif single_omics_tissues >= 1:
            score = 0.5
        else:
            score = 0.0
        return score, concordant_tissues

    def _effect_magnitude(self, evidence: dict[str, Any]) -> tuple[float, float]:
        """Score 0-0.5 by the max representative |log2FC| across all evidence."""
        max_fc = 0.0
        for ev in evidence.values():
            if not isinstance(ev, dict):
                continue
            for layer_key in ("transcriptomics", "metabolomics"):
                layer = ev.get(layer_key)
                if isinstance(layer, dict):
                    fc = layer.get("log2fc")
                    if fc is not None and not pd.isna(fc):
                        max_fc = max(max_fc, abs(fc))
        if max_fc >= self.large_log2fc:
            tier = 0.5
        elif max_fc >= 1.0:
            tier = 0.4
        elif max_fc >= self.log2fc_cutoff:
            tier = 0.3
        elif max_fc > 0:
            tier = 0.1
        else:
            tier = 0.0
        return tier, round(max_fc, 2)

    def _novelty(self, evidence: dict[str, Any]) -> float:
        """旧矩阵的重复组织奖励；仅用于 legacy 结果复现。"""
        n = 0
        for ev in evidence.values():
            if isinstance(ev, dict) and self._has_significant_evidence(ev):
                n += 1
        if n >= 2:
            return 0.5
        if n == 1:
            return 0.2
        return 0.0

    def _evidence_profile(self, evidence: dict[str, Any]) -> dict[str, Any]:
        """为人工策展证据生成不夸大置信度的多维画像。

        YAML 当前没有 SE、qvalue、n_eff 或特征角色，因此只把达到效应阈值的
        条目作为探索性支持；冲突仅接受显式 ``direction_relation: opposite``，
        不根据代谢物正负号猜测底物/产物关系。
        """
        tissue_strengths: list[float] = []
        all_effects: list[float] = []
        explicit_conflicts: list[float] = []
        n_layer_entries = 0
        tissues_with_signal = 0

        for tissue in self.tissues:
            tissue_evidence = evidence.get(tissue, {})
            if not isinstance(tissue_evidence, dict):
                continue
            layer_strengths: list[float] = []
            for layer_key in ("transcriptomics", "metabolomics"):
                layer = tissue_evidence.get(layer_key)
                if not isinstance(layer, dict) or not layer:
                    continue
                n_layer_entries += 1
                effect = layer.get("log2fc")
                if effect is None or pd.isna(effect):
                    continue
                abs_effect = abs(float(effect))
                all_effects.append(abs_effect)
                if abs_effect < self.log2fc_cutoff:
                    continue
                strength = min(abs_effect / self.large_log2fc, 1.0)
                if layer.get("direction_relation") == "opposite":
                    explicit_conflicts.append(strength)
                else:
                    layer_strengths.append(strength)
            if layer_strengths:
                tissues_with_signal += 1
                tissue_strengths.append(float(np.median(layer_strengths)))

        # 人工策展/effect-only 可靠性固定为探索性水平，不能冒充重复实验统计证据。
        curated_reliability = 0.35
        support_score = (
            float(np.median(tissue_strengths)) * curated_reliability
            if tissue_strengths else 0.0
        )
        conflict_score = (
            float(np.median(explicit_conflicts)) * curated_reliability
            if explicit_conflicts else 0.0
        )
        possible_entries = max(len(self.tissues) * 2, 1)
        coverage_score = n_layer_entries / possible_entries
        tissue_recurrence = tissues_with_signal / max(len(self.tissues), 1)
        median_effect = float(np.median(all_effects)) if all_effects else 0.0

        return {
            "support_score": round(float(np.clip(support_score, 0.0, 1.0)), 4),
            "conflict_score": round(float(np.clip(conflict_score, 0.0, 1.0)), 4),
            "coverage_score": round(float(np.clip(coverage_score, 0.0, 1.0)), 4),
            "tissue_recurrence": round(float(np.clip(tissue_recurrence, 0.0, 1.0)), 4),
            "median_abs_log2fc": round(median_effect, 3),
            "stability_selection_probability": float("nan"),
            "stability_direction": float("nan"),
            "design_grade": "D",
            "evidence_grade": "Exploratory" if n_layer_entries else "Not-assessable",
            "evidence_mode": "curated_effect_only" if n_layer_entries else "no_data",
        }

    @staticmethod
    def _comparison_sources(evidence: dict[str, Any]) -> str:
        """Derive the set of comparison sources across tissues, as friendly labels.

        Returns a ' + '-joined string of normalized comparison names (e.g.
        'Mazdutide_vs_Saline + Semaglutide_vs_Saline'). Preserves first-seen
        order across tissues.
        """
        seen: list[str] = []
        for ev in evidence.values():
            if not isinstance(ev, dict):
                continue
            comp = ev.get("comparison")
            if comp and comp not in seen:
                seen.append(comp)
        labels = [_COMPARISON_LABELS.get(c, c) for c in seen]
        return " + ".join(labels)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def score_pathway(self, pathway: dict[str, Any]) -> dict[str, Any]:
        """Score one systemic pathway. Returns dict with score components."""
        evidence = pathway.get("evidence", {})
        breadth, tissues_present = self._tissue_breadth(evidence)
        concordance, n_concordant = self._omics_concordance(evidence)
        magnitude, max_fc = self._effect_magnitude(evidence)
        novelty = self._novelty(evidence)
        profile = self._evidence_profile(evidence)

        # Round components to avoid float artifacts (e.g. 0.8999999999999999 -> 0.9)
        breadth = round(breadth, 2)
        concordance = round(concordance, 2)
        magnitude = round(magnitude, 2)
        novelty = round(novelty, 2)

        total = round(breadth + concordance + magnitude + novelty, 2)
        confidence = _confidence_label(total, self.confidence_thresholds)

        return {
            "id": pathway.get("id", ""),
            "pathway": pathway.get("name", ""),
            "category": pathway.get("category", ""),
            "drug": pathway.get("drug", ""),
            "comparison_sources": self._comparison_sources(evidence),
            "tissue_breadth": breadth,
            "n_tissues": len(tissues_present),
            "tissues_present": ", ".join(tissues_present),
            "omics_concordance": concordance,
            "n_concordant_tissues": n_concordant,
            "effect_magnitude": magnitude,
            "max_log2fc": max_fc,
            "novelty": novelty,
            "total_score": total,
            "confidence": confidence,
            "legacy_tissue_breadth": breadth,
            "legacy_omics_concordance": concordance,
            "legacy_effect_magnitude": magnitude,
            "legacy_novelty": novelty,
            "legacy_total_score": total,
            "legacy_confidence": confidence,
            **profile,
            "description": pathway.get("description", ""),
            "source_reports": ", ".join(pathway.get("source_reports", [])),
        }

    def score_all(self, pathways: list[dict[str, Any]]) -> pd.DataFrame:
        """评估所有通路；优先按新证据画像排序。"""
        rows = [self.score_pathway(p) for p in pathways]
        df = pd.DataFrame(rows)
        df = df.sort_values(
            ["support_score", "conflict_score", "tissue_recurrence", "coverage_score"],
            ascending=[False, True, False, False],
        ).reset_index(drop=True)
        return df

    def build_matrix(self, pathways: list[dict[str, Any]]) -> pd.DataFrame:
        """Build pathway x tissue evidence matrix (direction + top feature)."""
        rows = []
        for p in pathways:
            row: dict[str, Any] = {
                "id": p.get("id", ""),
                "pathway": p.get("name", ""),
            }
            for tissue in self.tissues:
                ev = p.get("evidence", {}).get(tissue, {})
                if isinstance(ev, dict) and self._has_significant_evidence(ev):
                    direction = ev.get("direction", "")
                    # Collect candidate features from both omics layers
                    candidates = []
                    for layer_key in ("transcriptomics", "metabolomics"):
                        layer = ev.get(layer_key)
                        if isinstance(layer, dict) and layer:
                            fname = layer.get("gene") or layer.get("metabolite", "")
                            ffc = layer.get("log2fc")
                            if fname:
                                candidates.append((fname, ffc))
                    # Prefer feature whose sign matches the pathway direction
                    # (avoids confusing "UP gene(-0.6)" for mismatch pathways)
                    feat = ""
                    dir_sign = 1 if direction == "UP" else (-1 if direction == "DOWN" else 0)
                    preferred = None
                    for fname, ffc in candidates:
                        if isinstance(ffc, (int, float)) and _sign(ffc) == dir_sign and dir_sign != 0:
                            preferred = (fname, ffc)
                            break
                    fname, ffc = preferred if preferred else (candidates[0] if candidates else ("", None))
                    if fname:
                        feat = f"{fname}({ffc:+.2f})" if isinstance(ffc, (int, float)) else str(fname)
                    row[tissue] = f"{direction} {feat}".strip()
                else:
                    row[tissue] = "—"
            rows.append(row)
        return pd.DataFrame(rows)
