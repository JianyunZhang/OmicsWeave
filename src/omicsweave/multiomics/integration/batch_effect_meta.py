#!/usr/bin/env python3
"""将独立批次的特征效应量表合并为随机效应 Meta 结果。"""

from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd
from scipy.stats import norm

from omicsweave.common.meta_analysis import random_effects_meta


def _bh_fdr(values: Iterable[float]) -> np.ndarray:
    pvalues = np.asarray(list(values), dtype=float)
    adjusted = np.full(pvalues.shape, np.nan, dtype=float)
    valid = np.isfinite(pvalues)
    if not valid.any():
        return adjusted
    finite = pvalues[valid]
    order = np.argsort(finite)
    ranked = finite[order]
    qvalues = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    qvalues = np.minimum.accumulate(qvalues[::-1])[::-1]
    restored = np.empty_like(qvalues)
    restored[order] = np.clip(qvalues, 0.0, 1.0)
    adjusted[valid] = restored
    return adjusted


def _pick_numeric_column(table: pd.DataFrame, candidates: Iterable[str]) -> str | None:
    for column in candidates:
        if column in table.columns and pd.to_numeric(table[column], errors="coerce").notna().any():
            return column
    return None


def _standard_error(
    table: pd.DataFrame,
    effect_column: str,
    se_candidates: Iterable[str],
    pvalue_candidates: Iterable[str],
) -> tuple[pd.Series, pd.Series]:
    se_column = _pick_numeric_column(table, se_candidates)
    if se_column:
        return (
            pd.to_numeric(table[se_column], errors="coerce").abs(),
            pd.Series("reported", index=table.index, dtype=object),
        )
    pvalue_column = _pick_numeric_column(table, pvalue_candidates)
    effects = pd.to_numeric(table[effect_column], errors="coerce")
    if pvalue_column:
        pvalues = pd.to_numeric(table[pvalue_column], errors="coerce")
        clipped = pvalues.clip(lower=np.finfo(float).tiny, upper=1.0 - np.finfo(float).eps)
        z_values = pd.Series(norm.isf(clipped / 2.0), index=table.index)
        inferred = effects.abs().divide(z_values.where(z_values > 0))
        inferred = inferred.replace([np.inf, -np.inf], np.nan)
        source = pd.Series("pvalue_normal_approximation", index=table.index, dtype=object)
        source[inferred.isna()] = "unavailable"
        return inferred, source
    return (
        pd.Series(np.nan, index=table.index, dtype=float),
        pd.Series("unavailable", index=table.index, dtype=object),
    )


def meta_analyze_feature_batches(
    batch_tables: dict[str, pd.DataFrame],
    *,
    feature_column: str,
    effect_column: str,
    se_candidates: Iterable[str],
    pvalue_candidates: Iterable[str],
    min_batches: int = 2,
) -> pd.DataFrame:
    """按 feature 跨批次合并；单批覆盖仅作为 descriptive passthrough。"""
    rows: list[pd.DataFrame] = []
    for batch_id, table in batch_tables.items():
        if table.empty or feature_column not in table.columns or effect_column not in table.columns:
            continue
        working = table.copy()
        working["_effect"] = pd.to_numeric(working[effect_column], errors="coerce")
        working["_se"], working["_se_source"] = _standard_error(
            working, effect_column, se_candidates, pvalue_candidates
        )
        working["_batch_id"] = str(batch_id)
        working = working.dropna(subset=[feature_column, "_effect"])
        if working.empty:
            continue
        # 同一批次同一 feature 只保留不确定性最小的一行。
        working = working.sort_values("_se", na_position="last").drop_duplicates(feature_column)
        rows.append(working[[feature_column, "_effect", "_se", "_se_source", "_batch_id"]])

    columns = [
        feature_column, "pooled_effect", "pooled_se", "pooled_pvalue", "pooled_fdr",
        "ci_lower", "ci_upper", "n_batches", "n_batches_observed",
        "n_batches_eligible", "batch_ids", "eligible_batch_ids", "i_squared",
        "tau_squared", "q_statistic", "direction_consistent", "meta_eligible",
        "meta_note", "se_source",
    ]
    if not rows:
        return pd.DataFrame(columns=columns)
    joined = pd.concat(rows, ignore_index=True)
    output: list[dict[str, object]] = []
    for feature, group in joined.groupby(feature_column, sort=True):
        observed = group.drop_duplicates("_batch_id")
        observed_ids = sorted(observed["_batch_id"].astype(str).tolist())
        valid = group[np.isfinite(group["_se"]) & (group["_se"] > 0)].copy()
        valid = valid.drop_duplicates("_batch_id")
        eligible_ids = sorted(valid["_batch_id"].astype(str).tolist())
        if valid.empty:
            effects = pd.to_numeric(observed["_effect"], errors="coerce").dropna()
            output.append({
                feature_column: feature,
                "pooled_effect": float(effects.mean()) if not effects.empty else np.nan,
                "pooled_se": np.nan,
                "pooled_pvalue": np.nan,
                "pooled_fdr": np.nan,
                "ci_lower": np.nan,
                "ci_upper": np.nan,
                "n_batches": 0,
                "n_batches_observed": len(observed_ids),
                "n_batches_eligible": 0,
                "batch_ids": ";".join(observed_ids),
                "eligible_batch_ids": "",
                "i_squared": np.nan,
                "tau_squared": np.nan,
                "q_statistic": np.nan,
                "direction_consistent": bool(
                    effects.empty or (effects.ge(0).all() or effects.le(0).all())
                ),
                "meta_eligible": False,
                "meta_note": "not_meta_eligible_missing_uncertainty",
                "se_source": "unavailable",
            })
            continue
        if valid["_batch_id"].nunique() == 1:
            one = valid.iloc[0]
            effect = float(one["_effect"])
            se = float(one["_se"])
            pvalue = float(2.0 * norm.sf(abs(effect / se)))
            output.append({
                feature_column: feature,
                "pooled_effect": effect,
                "pooled_se": se,
                "pooled_pvalue": pvalue,
                "pooled_fdr": np.nan,
                "ci_lower": effect - 1.96 * se,
                "ci_upper": effect + 1.96 * se,
                "n_batches": 1,
                "n_batches_observed": len(observed_ids),
                "n_batches_eligible": 1,
                "batch_ids": ";".join(observed_ids),
                "eligible_batch_ids": str(one["_batch_id"]),
                "i_squared": np.nan,
                "tau_squared": 0.0,
                "q_statistic": 0.0,
                "direction_consistent": True,
                "meta_eligible": False,
                "meta_note": "single_batch_passthrough",
                "se_source": str(one["_se_source"]),
            })
            continue
        if valid["_batch_id"].nunique() < min_batches:
            continue
        meta = random_effects_meta(
            valid["_effect"].to_numpy(dtype=float),
            valid["_se"].to_numpy(dtype=float),
            valid["_batch_id"].astype(str).tolist(),
        )
        output.append({
            feature_column: feature,
            "pooled_effect": meta.pooled_effect,
            "pooled_se": meta.pooled_se,
            "pooled_pvalue": meta.pvalue,
            "pooled_fdr": np.nan,
            "ci_lower": meta.ci_lower,
            "ci_upper": meta.ci_upper,
            "n_batches": meta.k,
            "n_batches_observed": len(observed_ids),
            "n_batches_eligible": meta.k,
            "batch_ids": ";".join(observed_ids),
            "eligible_batch_ids": ";".join(meta.source_ids),
            "i_squared": meta.i_squared,
            "tau_squared": meta.tau_squared,
            "q_statistic": meta.q_statistic,
            "direction_consistent": meta.direction_consistent,
            "meta_eligible": True,
            "meta_note": "random_effects_meta",
            "se_source": ";".join(sorted(valid["_se_source"].astype(str).unique())),
        })
    result = pd.DataFrame(output, columns=columns)
    if not result.empty:
        result["pooled_fdr"] = _bh_fdr(result["pooled_pvalue"])
    return result
