#!/usr/bin/env python3
"""跨组学共用的逆方差 Meta 分析。"""

from __future__ import annotations

from dataclasses import dataclass
from math import erf, sqrt
from typing import Iterable, Sequence

import numpy as np


@dataclass(frozen=True)
class MetaAnalysisResult:
    pooled_effect: float
    pooled_se: float
    ci_lower: float
    ci_upper: float
    pvalue: float
    q_statistic: float
    i_squared: float
    tau_squared: float
    k: int
    source_ids: tuple[str, ...]
    direction_consistent: bool


def _two_sided_normal_pvalue(z_value: float) -> float:
    return float(1.0 - erf(abs(z_value) / sqrt(2.0)))


def random_effects_meta(
    effects: Iterable[float],
    standard_errors: Iterable[float],
    source_ids: Sequence[str] | None = None,
) -> MetaAnalysisResult:
    """DerSimonian–Laird 随机效应 Meta；每个来源必须有有限正 SE。"""
    effect_array = np.asarray(list(effects), dtype=float)
    se_array = np.asarray(list(standard_errors), dtype=float)
    if effect_array.shape != se_array.shape:
        raise ValueError("effects and standard_errors must have the same length")
    labels = list(source_ids or [str(i + 1) for i in range(len(effect_array))])
    if len(labels) != len(effect_array):
        raise ValueError("source_ids must match the number of effects")

    valid = np.isfinite(effect_array) & np.isfinite(se_array) & (se_array > 0)
    effect_array = effect_array[valid]
    se_array = se_array[valid]
    labels = [str(label) for label, keep in zip(labels, valid, strict=True) if keep]
    k = len(effect_array)
    if k == 0:
        raise ValueError("at least one finite effect with a positive standard error is required")

    fixed_weights = 1.0 / np.square(se_array)
    fixed_effect = float(np.sum(fixed_weights * effect_array) / np.sum(fixed_weights))
    q_statistic = float(np.sum(fixed_weights * np.square(effect_array - fixed_effect)))
    if k > 1:
        c_value = float(
            np.sum(fixed_weights)
            - np.sum(np.square(fixed_weights)) / np.sum(fixed_weights)
        )
        tau_squared = max(0.0, (q_statistic - (k - 1)) / c_value) if c_value > 0 else 0.0
        i_squared = max(0.0, (q_statistic - (k - 1)) / q_statistic) * 100.0 if q_statistic > 0 else 0.0
    else:
        tau_squared = 0.0
        i_squared = float("nan")

    random_weights = 1.0 / (np.square(se_array) + tau_squared)
    pooled_effect = float(np.sum(random_weights * effect_array) / np.sum(random_weights))
    pooled_se = float(sqrt(1.0 / np.sum(random_weights)))
    z_value = pooled_effect / pooled_se
    nonzero = effect_array[effect_array != 0]
    return MetaAnalysisResult(
        pooled_effect=pooled_effect,
        pooled_se=pooled_se,
        ci_lower=pooled_effect - 1.96 * pooled_se,
        ci_upper=pooled_effect + 1.96 * pooled_se,
        pvalue=_two_sided_normal_pvalue(z_value),
        q_statistic=q_statistic,
        i_squared=float(i_squared),
        tau_squared=float(tau_squared),
        k=k,
        source_ids=tuple(labels),
        direction_consistent=bool(
            len(nonzero) == 0 or np.all(nonzero > 0) or np.all(nonzero < 0)
        ),
    )
