#!/usr/bin/env python3
"""跨组学共用的生物学重复 bootstrap 与标签置换校准。"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class StabilityResult:
    observed_effect: float
    calibrated_effect_threshold: float
    selection_probability: float
    direction_stability: float
    permutation_pvalue: float
    n_case: int
    n_control: int
    evidence_mode: str

    def to_dict(self) -> dict[str, float | int | str]:
        return asdict(self)


def calibrate_biological_stability(
    case_values: Iterable[float],
    control_values: Iterable[float],
    *,
    n_bootstrap: int = 2000,
    n_permutations: int = 5000,
    alpha: float = 0.05,
    random_state: int | None = 0,
) -> StabilityResult:
    """按独立生物学重复执行组内 bootstrap 和标签置换。

    任一组少于 3 个重复时只返回描述性效应，稳定性和置换字段为 NaN。
    输入不得使用技术重复、MCMC 拟合点或同一样本的多次测量替代生物学重复。
    """
    case = _finite_array(case_values)
    control = _finite_array(control_values)
    observed = float(case.mean() - control.mean()) if len(case) and len(control) else float("nan")

    if len(case) < 3 or len(control) < 3:
        return StabilityResult(
            observed_effect=observed,
            calibrated_effect_threshold=float("nan"),
            selection_probability=float("nan"),
            direction_stability=float("nan"),
            permutation_pvalue=float("nan"),
            n_case=len(case),
            n_control=len(control),
            evidence_mode="descriptive",
        )
    if n_bootstrap < 1 or n_permutations < 1:
        raise ValueError("n_bootstrap 和 n_permutations 必须为正整数")
    if not 0 < alpha < 1:
        raise ValueError("alpha 必须位于 0 和 1 之间")

    rng = np.random.default_rng(random_state)
    combined = np.concatenate([case, control])
    n_case = len(case)
    null_effects = np.empty(n_permutations, dtype=float)
    for index in range(n_permutations):
        permuted = rng.permutation(combined)
        null_effects[index] = permuted[:n_case].mean() - permuted[n_case:].mean()

    calibrated_threshold = float(np.quantile(np.abs(null_effects), 1.0 - alpha))
    permutation_pvalue = float(
        (1 + np.count_nonzero(np.abs(null_effects) >= abs(observed)))
        / (n_permutations + 1)
    )

    bootstrap_effects = np.empty(n_bootstrap, dtype=float)
    for index in range(n_bootstrap):
        case_sample = rng.choice(case, size=len(case), replace=True)
        control_sample = rng.choice(control, size=len(control), replace=True)
        bootstrap_effects[index] = case_sample.mean() - control_sample.mean()

    selection_probability = float(
        np.mean(np.abs(bootstrap_effects) >= calibrated_threshold)
    )
    if observed > 0:
        direction_stability = float(np.mean(bootstrap_effects > 0))
    elif observed < 0:
        direction_stability = float(np.mean(bootstrap_effects < 0))
    else:
        direction_stability = 0.0

    return StabilityResult(
        observed_effect=observed,
        calibrated_effect_threshold=calibrated_threshold,
        selection_probability=selection_probability,
        direction_stability=direction_stability,
        permutation_pvalue=permutation_pvalue,
        n_case=len(case),
        n_control=len(control),
        evidence_mode="biological_replicate",
    )


def _finite_array(values: Iterable[float]) -> np.ndarray:
    array = np.asarray(list(values), dtype=float)
    return array[np.isfinite(array)]
