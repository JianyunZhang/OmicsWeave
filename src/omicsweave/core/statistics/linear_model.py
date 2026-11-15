#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
多特征线性模型与 limma 式方差压缩（G3）。

对特征 × 样本矩阵 ``Y`` 和共同的设计矩阵 ``X``（样本 × 系数）逐特征做最小二乘，再把各特征的残差方差向共同先验
收缩（:func:`omicsweave.core.statistics.moderated.estimate_limma_priors`），给出每个对比（系数的线性组合）的
估计、标准误、moderated t、自由度与 p 值。没有缺失值的特征按矩阵一次求解；有缺失值的特征只用观测到的样本单独
拟合，观测样本使设计不满秩或没有残差自由度时该特征不检验。

用于代谢组、派生特征与通量等连续型数据层（均在 log2 尺度上）；计数数据仍用 PyDESeq2。
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import scipy.linalg
from scipy import stats

from .moderated import estimate_limma_priors


class DesignError(ValueError):
    """设计矩阵不满秩或没有残差自由度。"""


def check_rank(x: np.ndarray, columns: list[str]) -> None:
    """设计矩阵不满秩时说明哪些列可由其他列线性表示（QR 列主元分解中排在秩之后的列）。"""
    n, p = x.shape
    if p == 0:
        raise DesignError("设计矩阵没有列")
    _, r, pivots = scipy.linalg.qr(x, mode="economic", pivoting=True)
    diagonal = np.abs(np.diag(r))
    tolerance = diagonal.max() * max(n, p) * np.finfo(float).eps if diagonal.size else 0.0
    rank = int((diagonal > tolerance).sum())
    if rank < p:
        aliased = [columns[i] for i in pivots[rank:]]
        raise DesignError(f"设计矩阵不满秩（秩 {rank} < {p} 列）：{'、'.join(aliased)} 可由其他列线性表示，"
                          "常见原因是某个因素水平只出现在另一个因素的一个水平中（如处理只在一种饮食中给）")
    if n - p < 1:
        raise DesignError(f"没有残差自由度：{n} 个样本、{p} 个系数")


@dataclass
class ContrastResult:
    """一个对比的逐特征结果，索引为特征。"""
    table: pd.DataFrame          # estimate, se, statistic, df, pvalue, residual_sd, n_obs
    prior_variance: float
    prior_df: float


def fit_contrasts(y: pd.DataFrame, x: np.ndarray, contrasts: dict[str, np.ndarray], *,
                  moderate: bool = True, min_residual_df: int = 1) -> dict[str, ContrastResult]:
    """``y``：特征 × 样本（列顺序与 ``x`` 的行一致）；``contrasts``：名称 → 长度为系数个数的向量。"""
    values = y.to_numpy(dtype=float)
    n_features, n_samples = values.shape
    if x.shape[0] != n_samples:
        raise DesignError(f"设计矩阵有 {x.shape[0]} 行，数据有 {n_samples} 个样本")
    p = x.shape[1]
    vectors = {name: np.asarray(c, dtype=float) for name, c in contrasts.items()}
    estimates = {name: np.full(n_features, np.nan) for name in vectors}
    unscaled = {name: np.full(n_features, np.nan) for name in vectors}
    variance = np.full(n_features, np.nan)
    df = np.full(n_features, np.nan)
    n_obs = np.isfinite(values).sum(axis=1)

    complete = n_obs == n_samples
    if complete.any() and n_samples - p >= min_residual_df:
        xtx_inv = np.linalg.pinv(x.T @ x)
        coef = xtx_inv @ x.T @ values[complete].T                  # p × G
        resid = values[complete] - (x @ coef).T
        variance[complete] = (resid ** 2).sum(axis=1) / (n_samples - p)
        df[complete] = n_samples - p
        for name, c in vectors.items():
            estimates[name][complete] = c @ coef
            unscaled[name][complete] = c @ xtx_inv @ c
    for i in np.flatnonzero(~complete):
        observed = np.isfinite(values[i])
        xo = x[observed]
        if xo.shape[0] - p < min_residual_df or np.linalg.matrix_rank(xo) < p:
            continue
        xtx_inv = np.linalg.pinv(xo.T @ xo)
        coef = xtx_inv @ xo.T @ values[i, observed]
        resid = values[i, observed] - xo @ coef
        df[i] = xo.shape[0] - p
        variance[i] = float(resid @ resid) / df[i]
        for name, c in vectors.items():
            estimates[name][i] = c @ coef
            unscaled[name][i] = c @ xtx_inv @ c

    valid = np.isfinite(variance) & np.isfinite(df) & (df > 0)
    prior_variance, prior_df = float("nan"), 0.0
    posterior_variance, posterior_df = variance.copy(), df.copy()
    if moderate and valid.sum() >= 3:
        prior_variance, prior_df = estimate_limma_priors(variance[valid], df[valid])
        if np.isinf(prior_df):
            posterior_variance[valid] = prior_variance
            posterior_df[valid] = np.inf
        elif prior_df > 0:
            posterior_variance[valid] = (df[valid] * variance[valid] + prior_df * prior_variance) / (df[valid] + prior_df)
            posterior_df[valid] = df[valid] + prior_df
    floor = max((prior_variance if np.isfinite(prior_variance) else np.nanmedian(variance[valid]) if valid.any()
                 else 1.0) * 1e-12, np.finfo(float).tiny)
    posterior_variance = np.where(valid, np.maximum(posterior_variance, floor), np.nan)

    results = {}
    for name in vectors:
        se = np.sqrt(posterior_variance * unscaled[name])
        with np.errstate(divide="ignore", invalid="ignore"):
            t = estimates[name] / se
        pvalue = np.where(np.isinf(posterior_df), 2 * stats.norm.sf(np.abs(t)), 2 * stats.t.sf(np.abs(t), posterior_df))
        table = pd.DataFrame({"estimate": estimates[name], "se": se, "statistic": t, "df": posterior_df,
                              "pvalue": pvalue, "residual_sd": np.sqrt(variance), "n_obs": n_obs}, index=y.index)
        table.loc[~valid, ["se", "statistic", "pvalue"]] = np.nan
        results[name] = ContrastResult(table, prior_variance, prior_df)
    return results


def bh(pvalues: pd.Series | np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg，只校正有限值，其余保留 NaN。"""
    p = np.asarray(pvalues, dtype=float)
    out = np.full(p.shape, np.nan)
    finite = np.isfinite(p)
    if finite.any():
        order = np.argsort(p[finite])
        ranked = p[finite][order] * finite.sum() / np.arange(1, finite.sum() + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        adjusted = np.empty_like(ranked)
        adjusted[order] = np.minimum(ranked, 1.0)
        out[finite] = adjusted
    return out
