#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
带符号集合的旋转检验（ROAST，Wu 等 2010；limma ``roast`` 的 ``mean`` 统计量）。

对一个特征集合与一个对比，每个特征在设计矩阵下分解出对比效应 ``e``（按对比的标准误归一，H0 下 ~ N(0, σ²)）与
``d = n − p`` 个残差分量 ``r``；H0 下 ``(e, r)`` 是球对称的。旋转：对集合内所有特征用同一个随机单位向量 ``u``，
令 ``e* = (e, r)·u``、``RSS* = ‖(e, r)‖² − e*²``，从而保留特征间的相关。逐特征 moderated t（先验 ``s0²``、``d0``
来自全体特征，旋转时固定）换算成 z，乘以预先声明的符号后取平均，即集合统计量。

与置换检验相比，旋转在每组 3 只时也给出连续的零分布（置换只有 10 种），适合本课题的设计。符号由知识预先给定，
全部取反时统计量变号、双侧 p 不变（同一随机种子下逐次旋转一一对应）。供反应档案的模块检验（G16）与 G9 共用。
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.linalg
from scipy import stats

ACTIVE_Z = float(np.sqrt(2.0))     # limma 的 PropUp / PropDown：|t| > √2 的成员占比


@dataclass(frozen=True)
class RotationResult:
    statistic: float        # 带符号的 z 的平均（观测）
    pvalue: float           # 双侧：min(1, 2·min(上调 p, 下调 p))
    pvalue_up: float        # 集合按声明的符号上调
    pvalue_down: float      # 集合按声明的符号的反方向变化
    direction: int          # 统计量的符号
    n_members: int
    active_up: float        # 按声明符号同向且 z > √2 的成员占比
    active_down: float
    rotations: int
    df: float               # moderated t 的自由度 d + d0


def _z_from_t(t: np.ndarray, df: float) -> np.ndarray:
    if np.isinf(df):
        return t
    with np.errstate(divide="ignore"):
        tail = stats.t.logsf(np.abs(t), df)
    return np.sign(t) * stats.norm.isf(np.exp(np.maximum(tail, -700.0)))


def signed_rotation_test(y: np.ndarray, x: np.ndarray, contrast: np.ndarray, signs: np.ndarray, *,
                         prior_variance: float = float("nan"), prior_df: float = 0.0, rotations: int = 9999,
                         seed: int | None = None) -> RotationResult:
    """``y``：成员 × 样本（log 尺度，不能有缺失）；``x``：样本 × 系数；``contrast``：系数的线性组合；
    ``signs``：每个成员的 +1 / −1。``prior_variance`` / ``prior_df`` 为全体特征的 limma 先验（没有时 d0 = 0）。"""
    y = np.atleast_2d(np.asarray(y, dtype=float))
    x = np.asarray(x, dtype=float)
    c = np.asarray(contrast, dtype=float)
    signs = np.sign(np.asarray(signs, dtype=float))
    n, p = x.shape
    if y.shape[1] != n or len(signs) != y.shape[0]:
        raise ValueError("y、x、signs 的维度不一致")
    if not np.isfinite(y).all():
        raise ValueError("旋转检验需要完整数据：先去掉有缺失值的成员")
    d = n - np.linalg.matrix_rank(x)
    if d < 1:
        raise ValueError("设计没有残差自由度")
    # 对比方向：c·β̂ = aᵀy，a = X(XᵀX)⁻¹c；e = aᵀy / ‖a‖ 与对比估计同号，H0 下方差为 σ²
    a = x @ np.linalg.pinv(x.T @ x) @ c
    a = a / np.linalg.norm(a)
    q, _ = scipy.linalg.qr(x)                      # 完整 QR：后 n − rank 列张成残差空间
    residual_basis = q[:, n - d:]
    effects = np.column_stack([y @ a, y @ residual_basis])          # G × (d + 1)

    d0 = float(prior_df) if np.isfinite(prior_variance) and prior_df > 0 else 0.0
    s0 = float(prior_variance) if d0 > 0 else 0.0
    df = np.inf if np.isinf(d0) else d + d0

    def set_statistic(e: np.ndarray, rss: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        variance = np.full_like(rss, s0) if np.isinf(d0) else (d0 * s0 + rss) / (d0 + d)
        z = _z_from_t(e / np.sqrt(np.maximum(variance, np.finfo(float).tiny)), df) * signs
        return z.mean(axis=-1), z

    total = (effects ** 2).sum(axis=1)
    observed, z = set_statistic(effects[:, 0], total - effects[:, 0] ** 2)
    rng = np.random.default_rng(seed)
    u = rng.standard_normal((rotations, d + 1))
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    rotated_e = u @ effects.T                                         # rotations × G
    rotated, _ = set_statistic(rotated_e, total[None, :] - rotated_e ** 2)
    tolerance = 1e-12 * max(1.0, abs(float(observed)))
    up = (1 + np.sum(rotated >= observed - tolerance)) / (rotations + 1)
    down = (1 + np.sum(rotated <= observed + tolerance)) / (rotations + 1)
    return RotationResult(statistic=float(observed), pvalue=float(min(1.0, 2 * min(up, down))),
                          pvalue_up=float(up), pvalue_down=float(down), direction=int(np.sign(observed)),
                          n_members=int(y.shape[0]), active_up=float(np.mean(z > ACTIVE_Z)),
                          active_down=float(np.mean(z < -ACTIVE_Z)), rotations=int(rotations), df=float(df))
