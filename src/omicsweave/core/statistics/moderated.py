#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
limma 风格的经验贝叶斯方差压缩（Smyth 2004；limma ``fitFDist`` / ``fitFDistUnequalDF1`` 的移植）。

小样本下逐特征的残差方差很不稳定，把各特征的方差向共同先验 ``s0²``（自由度 ``d0``）收缩后，moderated t 的
自由度为 ``d + d0``，检验更稳。最初实现于蛋白组差异分析，现由 :mod:`omicsweave.core.statistics.linear_model`
供各连续型数据层共用。
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import digamma, gammaln, polygamma


def _log_mdigamma(x: np.ndarray | float) -> np.ndarray | float:
    """limma/statmod 的 ``logmdigamma(x) = log(x) - digamma(x)``。"""

    return np.log(x) - digamma(x)


def _trigamma_inverse(value: float) -> float:
    """求解 ``trigamma(y) = value``，迭代式与 limma 相同。"""

    if not np.isfinite(value):
        return 0.0 if value > 0 else np.nan
    if value < 0:
        return np.nan
    if value > 1e7:
        return 1.0 / np.sqrt(value)
    if value < 1e-6:
        return 1.0 / max(value, np.finfo(float).tiny)

    y = 0.5 + 1.0 / value
    for _ in range(51):
        tri = float(polygamma(1, y))
        tetra = float(polygamma(2, y))
        difference = tri * (1.0 - tri / value) / tetra
        y += difference
        if -difference / y < 1e-8:
            break
    return float(y)


def _fit_f_dist_equal_df(variances: np.ndarray, dfs: np.ndarray) -> tuple[float, float]:
    """移植 limma ``fitFDist`` 的等自由度 log-F 矩估计。"""

    x = np.asarray(variances, dtype=float)
    df1 = np.asarray(dfs, dtype=float)
    valid = np.isfinite(x) & np.isfinite(df1) & (x >= 0) & (df1 > 1e-15)
    x = x[valid]
    df1 = df1[valid]
    if x.size < 3:
        positive = x[x > 0]
        return (float(np.median(positive)) if positive.size else 1.0), 0.0

    median_variance = float(np.median(x))
    if median_variance <= 0:
        median_variance = float(np.mean(x[x > 0])) if np.any(x > 0) else 1.0
    x = np.maximum(x, 1e-5 * median_variance)

    half_df = df1 / 2.0
    adjusted_log_variance = np.log(x) + _log_mdigamma(half_df)
    mean_adjusted = float(np.mean(adjusted_log_variance))
    variance_adjusted = float(np.var(adjusted_log_variance, ddof=1))
    variance_adjusted -= float(np.mean(polygamma(1, half_df)))

    if variance_adjusted > 0:
        prior_df = 2.0 * _trigamma_inverse(variance_adjusted)
        prior_variance = float(
            np.exp(mean_adjusted - _log_mdigamma(prior_df / 2.0))
        )
    else:
        prior_df = np.inf
        prior_variance = float(np.mean(x))
    return prior_variance, float(prior_df)


def _fit_f_dist_unequal_df(
    variances: np.ndarray,
    dfs: np.ndarray,
) -> tuple[float, float]:
    """移植 limma >=3.62 ``fitFDistUnequalDF1`` 的 profile MLE 核心。"""

    x = np.asarray(variances, dtype=float)
    df1 = np.asarray(dfs, dtype=float)
    valid = np.isfinite(x) & np.isfinite(df1) & (x >= 0) & (df1 >= 0.01)
    x = x[valid]
    df1 = df1[valid]
    if x.size < 3:
        positive = x[x > 0]
        return (float(np.median(positive)) if positive.size else 1.0), 0.0

    positive = x[x > 0]
    median_variance = float(np.median(positive)) if positive.size else 1.0
    x_positive = np.maximum(x, 1e-12 * median_variance)
    d1 = df1 / 2.0
    adjusted = np.log(x_positive) + _log_mdigamma(d1)
    weights = 1.0 / polygamma(1, d1)
    adjusted_mean = float(np.sum(weights * adjusted) / np.sum(weights))
    d1_times_x = d1 * x_positive

    def minus_twice_log_likelihood(parameter: float) -> float:
        d2 = parameter / (1.0 - parameter)
        d2_s20 = d2 * np.exp(adjusted_mean - _log_mdigamma(d2))
        terms = (
            -(d1 + d2) * np.log1p(d1_times_x / d2_s20)
            - d1 * np.log(d2_s20)
            + gammaln(d1 + d2)
            - gammaln(d2)
        )
        return float(-2.0 * np.sum(terms))

    fit = minimize_scalar(
        minus_twice_log_likelihood,
        bounds=(0.5, 0.9998),
        method="bounded",
        options={"xatol": 1e-10},
    )
    if not fit.success or not np.isfinite(fit.fun):
        return _fit_f_dist_equal_df(x, np.repeat(float(np.median(df1)), len(x)))

    d2 = float(fit.x / (1.0 - fit.x))
    prior_variance = float(np.exp(adjusted_mean - _log_mdigamma(d2)))
    return prior_variance, 2.0 * d2


def estimate_limma_priors(
    variances: np.ndarray,
    dfs: np.ndarray,
) -> tuple[float, float]:
    """估计 limma 的先验方差 ``s0²`` 与先验自由度 ``d0``。

    等自由度时使用官方 ``fitFDist`` 的 log-F 矩估计；缺失值导致自由度
    不等时使用 limma 新版的 profile maximum-likelihood 路径。
    """

    variances = np.asarray(variances, dtype=float)
    dfs = np.asarray(dfs, dtype=float)
    valid = np.isfinite(variances) & np.isfinite(dfs) & (variances >= 0) & (dfs > 0)
    variances = variances[valid]
    dfs = dfs[valid]
    if variances.size < 3:
        positive = variances[variances > 0]
        return (float(np.median(positive)) if positive.size else 1.0), 0.0
    if np.allclose(dfs, dfs[0], rtol=0, atol=1e-12):
        return _fit_f_dist_equal_df(variances, dfs)
    return _fit_f_dist_unequal_df(variances, dfs)
