#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""跨组学共用：秩检验（Mann-Whitney / Wilcoxon 秩和）的样本量下界。

秩检验的 p 值不是连续的——它由秩的排列组合数决定。在小样本实验里
（MFA、低重复代谢组、n=3/组的转录组）这个离散下界经常**高于**配置的
显著性阈值，此时"0 个显著特征"是方法学上限，不是生物学结论。

本模块把这个界限显式化，供各组学模块在跑检验之前判断：
- ``min_achievable_pvalue``    —— 未校正 p 值的下界
- ``min_achievable_bh_fdr``    —— BH 校正后 FDR 的下界

转录组模块 ``run_wilcoxon_de.py`` 有一份等价的 ``min_achievable_pvalue``
实现（含单元测试 ``test_transcriptomics_fixes.py``）；本模块是 MFA 侧的
共享入口，两者组合数学一致。
"""

from __future__ import annotations

import math

__all__ = ["min_achievable_pvalue", "min_achievable_bh_fdr", "rank_test_power_note"]


def min_achievable_pvalue(n_case: int, n_control: int) -> float:
    """Mann-Whitney 精确检验在给定样本量下可能取到的**最小双侧 p 值**。

    共有 ``C(n1+n2, n1)`` 种秩分配，最极端的一种对应单侧 p = 1/C，
    双侧即 2/C。这个下界与数据无关——纯粹由样本量决定。

    实践含义（MFA 常见样本量）::

        n=1 vs 1  ->  C(2,1)=2    p_min = 1.0    完全无法推断
        n=2 vs 2  ->  C(4,2)=6    p_min = 0.333
        n=3 vs 3  ->  C(6,3)=20   p_min = 0.1    永远达不到 p<0.05
        n=3 vs 4  ->  C(7,3)=35   p_min = 0.057  永远达不到 p<0.05
        n=4 vs 4  ->  C(8,4)=70   p_min = 0.029  p<0.05 刚刚可达
        n=5 vs 5  ->  C(10,5)=252 p_min = 0.0079

    Returns
    -------
    float
        最小双侧 p 值；任一组样本量 < 1 时返回 NaN。
    """
    if n_case < 1 or n_control < 1:
        return float("nan")
    total_arrangements = math.comb(n_case + n_control, n_case)
    return min(1.0, 2.0 / total_arrangements)


def min_achievable_bh_fdr(p_floor: float, n_tests: int, n_at_floor: int) -> float:
    """BH 校正后可能取到的最小 FDR。

    BH 对秩为 ``i``（升序）的 p 值给出 ``p * m / i``，再从大到小做累积
    取最小以保持单调。若有 ``n_at_floor`` 个检验同时压在下界 ``p_floor``
    上，它们中最有利的秩是 ``i = n_at_floor``，因此下界为::

        p_floor * n_tests / n_at_floor

    只有 1 个检验压在下界时惩罚最重（乘以全部检验数），全部检验都压在
    下界时惩罚消失（退化回 ``p_floor``）。

    Parameters
    ----------
    p_floor : float
        ``min_achievable_pvalue`` 的返回值。
    n_tests : int
        参与同一次 BH 校正的检验总数。
    n_at_floor : int
        实际取到（或并列于）下界的检验数量，至少为 1。

    Returns
    -------
    float
        最小可达 FDR；输入非法时返回 NaN。
    """
    if not math.isfinite(p_floor) or n_tests < 1:
        return float("nan")
    n_at_floor = max(1, min(int(n_at_floor), int(n_tests)))
    return min(1.0, float(p_floor) * n_tests / n_at_floor)


def rank_test_power_note(
    n_case: int,
    n_control: int,
    pvalue_threshold: float,
) -> str:
    """生成一行人类可读的说明，用于日志与报告。"""
    p_floor = min_achievable_pvalue(n_case, n_control)
    if not math.isfinite(p_floor):
        return f"n={n_case} vs {n_control}: 样本量不足，无法进行秩检验"
    if p_floor > pvalue_threshold:
        return (
            f"n={n_case} vs {n_control}: 秩检验最小双侧 p = {p_floor:.4g} "
            f"> 阈值 {pvalue_threshold:g}，p 值路径在数学上不可能显著"
        )
    return (
        f"n={n_case} vs {n_control}: 秩检验最小双侧 p = {p_floor:.4g} "
        f"<= 阈值 {pvalue_threshold:g}"
    )
