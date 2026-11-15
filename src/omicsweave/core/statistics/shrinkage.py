#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
效应量的经验贝叶斯收缩与符号可信度。

- :func:`adaptive_shrinkage`：自适应收缩（Stephens 2017 的 ash，正态混合先验版本）。先验为零点质量加一组以 0 为中心、
  方差成等比网格的正态分布的混合，权重用 EM 估计（零点权重带惩罚，偏保守）。返回后验均值与局部假符号率 lfsr
  （后验中效应为 0 或与估计符号相反的概率的较小者），比 p 值更直接地回答“方向可信吗”；
- :func:`retrodesign`：给定可信的真实效应 D 与标准误 s，按 Gelman 与 Carlin（2014）计算功效、符号错误率（Type S）
  与显著结果的夸大倍数（Type M）；
- :func:`normal_shrinkage`：正态—正态模型下的收缩估计，用作 retrodesign 中“可信的真实效应”。
"""
from __future__ import annotations

import numpy as np
from scipy import stats


def _grid(beta: np.ndarray, se: np.ndarray) -> np.ndarray:
    smallest = max(float(np.min(se)) / 10.0, 1e-8)
    spread = float(np.max(beta ** 2 - se ** 2)) if beta.size else 0.0
    largest = 2.0 * np.sqrt(spread) if spread > 0 else 8.0 * smallest
    largest = max(largest, smallest * 2)
    count = int(np.ceil(np.log(largest / smallest) / np.log(np.sqrt(2.0)))) + 1
    return smallest * np.sqrt(2.0) ** np.arange(count)


def adaptive_shrinkage(beta, se, *, null_penalty: float = 10.0, max_iter: int = 500,
                       tol: float = 1e-8) -> tuple[np.ndarray, np.ndarray]:
    """返回 (后验均值, lfsr)；``beta`` 或 ``se`` 无效（非有限或 se ≤ 0）的位置为 NaN。"""
    beta = np.asarray(beta, dtype=float)
    se = np.asarray(se, dtype=float)
    mean, lfsr = np.full(beta.shape, np.nan), np.full(beta.shape, np.nan)
    ok = np.isfinite(beta) & np.isfinite(se) & (se > 0)
    if ok.sum() < 3:
        return mean, lfsr
    b, s = beta[ok], se[ok]
    sigma = np.concatenate([[0.0], _grid(b, s)])
    variance = sigma[None, :] ** 2 + s[:, None] ** 2
    likelihood = stats.norm.pdf(b[:, None], 0.0, np.sqrt(variance))
    prior = np.concatenate([[null_penalty], np.ones(len(sigma) - 1)])
    weights = np.full(len(sigma), 1.0 / len(sigma))
    for _ in range(max_iter):
        joint = likelihood * weights
        resp = joint / np.maximum(joint.sum(axis=1, keepdims=True), np.finfo(float).tiny)
        updated = resp.sum(axis=0) + prior - 1.0
        updated = np.maximum(updated, 0.0) / np.maximum(updated.sum(), np.finfo(float).tiny)
        if np.max(np.abs(updated - weights)) < tol:
            weights = updated
            break
        weights = updated
    joint = likelihood * weights
    resp = joint / np.maximum(joint.sum(axis=1, keepdims=True), np.finfo(float).tiny)
    with np.errstate(divide="ignore", invalid="ignore"):
        post_var = np.where(sigma[None, :] > 0, 1.0 / (1.0 / sigma[None, :] ** 2 + 1.0 / s[:, None] ** 2), 0.0)
        post_mean = post_var * b[:, None] / s[:, None] ** 2
        negative = np.where(sigma[None, :] > 0, stats.norm.cdf(-post_mean / np.sqrt(post_var)), 0.0)
    p_zero = resp[:, 0]
    p_negative = (resp[:, 1:] * negative[:, 1:]).sum(axis=1)
    p_positive = 1.0 - p_zero - p_negative
    mean[ok] = (resp * post_mean).sum(axis=1)
    lfsr[ok] = np.clip(np.minimum(p_zero + p_positive, p_zero + p_negative), 0.0, 1.0)
    return mean, lfsr


def normal_shrinkage(beta, se) -> np.ndarray:
    """正态—正态收缩：先验方差 τ² 由 E[β̂²] − E[s²] 的矩估计给出（不小于 0）。"""
    beta = np.asarray(beta, dtype=float)
    se = np.asarray(se, dtype=float)
    ok = np.isfinite(beta) & np.isfinite(se) & (se > 0)
    out = np.full(beta.shape, np.nan)
    if not ok.any():
        return out
    tau2 = max(float(np.mean(beta[ok] ** 2) - np.mean(se[ok] ** 2)), 0.0)
    out[ok] = beta[ok] * tau2 / (tau2 + se[ok] ** 2)
    return out


def retrodesign(true_effect, se, alpha: float = 0.05) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """返回 (功效, Type S 符号错误率, Type M 夸大倍数)；真实效应为 0 或无效时为 NaN。"""
    d = np.asarray(true_effect, dtype=float)
    s = np.asarray(se, dtype=float)
    z = stats.norm.isf(alpha / 2.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        lam = d / s
        upper = stats.norm.sf(z - lam)
        lower = stats.norm.cdf(-z - lam)
        power = upper + lower
        type_s = np.where(d > 0, lower, upper) / power
        a, b = z - lam, -z - lam
        positive_part = d * stats.norm.sf(a) + s * stats.norm.pdf(a)
        negative_part = -d * stats.norm.cdf(b) + s * stats.norm.pdf(b)
        exaggeration = (positive_part + negative_part) / power / np.abs(d)
    bad = ~(np.isfinite(d) & np.isfinite(s) & (s > 0) & (d != 0))
    return (np.where(bad, np.nan, power), np.where(bad, np.nan, type_s), np.where(bad, np.nan, exaggeration))
