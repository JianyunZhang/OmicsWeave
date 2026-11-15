#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
可检出效应表（G3，``qc/design_power.parquet``）：每个（组学, 数据集, 组织, 比较）给出样本量、合并组内标准差的中位数
（log2 尺度）、双侧检验在 80% 功效下的最小可检出效应（log2 倍数变化），以及秩检验在该样本量下可达到的最小 p 值。

最小可检出效应按两样本 t 检验近似：``(t_{1−α/2, df} + t_{0.8, df}) · s · √(1/n₁ + 1/n₂)``，``df = n₁ + n₂ − 2``。
计数数据在 log2(CPM + 1) 上计算，只作量级参考（DESeq2 的离散度模型与此不同）。
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from omicsweave.contract.dataset import read_dataset
from omicsweave.core.config import Experiment
from omicsweave.core.statistics.nonparametric import min_achievable_pvalue

from .data import analysis_values, bundle_datasets, pairwise_groups

COLUMNS = ["omics", "kind", "tissue", "contrast_id", "n_numerator", "n_denominator", "n_features",
           "median_residual_sd", "mde_80", "min_achievable_p", "alpha"]
POWER = 0.8


def pooled_sd(first: pd.DataFrame, second: pd.DataFrame) -> pd.Series:
    """逐特征的合并组内标准差（只有一个观测的组不贡献方差）。"""
    n1, n2 = first.notna().sum(axis=1), second.notna().sum(axis=1)
    ss1 = (first.var(axis=1, ddof=1) * (n1 - 1)).fillna(0.0)
    ss2 = (second.var(axis=1, ddof=1) * (n2 - 1)).fillna(0.0)
    df = n1 + n2 - 2
    return np.sqrt((ss1 + ss2) / df.where(df >= 1))


def minimum_detectable_effect(sd: float, n1: int, n2: int, alpha: float, power: float = POWER) -> float:
    df = n1 + n2 - 2
    if df < 1 or not np.isfinite(sd):
        return float("nan")
    return float((stats.t.isf(alpha / 2, df) + stats.t.ppf(power, df)) * sd * np.sqrt(1 / n1 + 1 / n2))


def design_power(experiment: Experiment, run_dir, alpha: float) -> pd.DataFrame:
    rows = []
    for omics, kind, tissue, path in bundle_datasets(run_dir):
        values = analysis_values(read_dataset(path))
        for contrast in experiment.config.design.contrasts:
            groups = pairwise_groups(experiment, contrast, omics, values.columns)
            num = [s for s, g in groups.items() if g == contrast.numerator]
            den = [s for s, g in groups.items() if g == contrast.denominator]
            if not num or not den:
                continue
            sd = pooled_sd(values[num], values[den]).dropna()
            median = float(sd.median()) if not sd.empty else float("nan")
            rows.append((omics, kind, tissue, contrast.id, len(num), len(den), int(len(sd)), median,
                         minimum_detectable_effect(median, len(num), len(den), alpha),
                         float(min_achievable_pvalue(len(num), len(den))), alpha))
    return pd.DataFrame(rows, columns=COLUMNS)
