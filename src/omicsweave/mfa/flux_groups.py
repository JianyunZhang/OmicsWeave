#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
化学计量上锁死的通量组（G10 的一部分）：在检验前把逐次拟合中始终成线性关系的反应合并为一组，组内只计一次检验。

稳态下，无分支线性片段上的反应（如 PFK 与 FBA、TCA 中 AKGD → FUMH）的净通量由质量守恒决定，彼此成比例或只差
一个常数。它们在差异检验中给出完全相同的 t / 秩统计量，分开计数会使 BH 校正的检验数虚增，也会把同一个通量的变化
报告成几个“发现”。

识别方法：把一个模型（一种示踪剂）全部样本的全部拟合（每次拟合一行，净通量 = 正向 − 逆向）按组织堆叠，
两条反应的净通量在所有拟合中的相关系数绝对值不低于 ``min_abs_correlation``（默认 1 − 1e−6）即视为锁死，
按传递闭包合并为组。几百到上千次独立拟合中的精确线性关系只能来自模型结构，而不是数据。所有拟合中都不变的
通量（常贴在约束边界上）不参与配对。多种示踪剂时，一对反应须在每个两者都有变化的模型中都锁死。

``flux_group`` 为组的标识（``组织--首个成员``），单独的反应为空；``slope`` / ``intercept`` 为该反应相对首个成员的
线性关系。:func:`grouped_bh` 在每组只取一个 p 值（组内最小）后做 BH，再把校正值赋给组内所有成员。
"""
from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd

COLUMNS = ["tissue", "reaction", "flux_group", "lead", "n_members", "members", "slope", "intercept", "scope"]
DEFAULT_MIN_ABS_CORRELATION = 1.0 - 1e-6


def net_fits(all_fits: dict[str, pd.DataFrame], samples: Iterable[str] | None = None) -> pd.DataFrame:
    """全部样本的逐次拟合净通量（行 = 拟合，列 = ``组织--反应``）。"""
    wanted = set(samples) if samples is not None else None
    frames = []
    for sample, fits in all_fits.items():
        if wanted is not None and sample not in wanted:
            continue
        columns = [c for c in fits.columns if "--" in str(c) and not str(c).startswith("MIX:")
                   and not str(c).endswith("__R")]
        net = {}
        for column in columns:
            forward = pd.to_numeric(fits[column], errors="coerce")
            reverse = f"{column}__R"
            net[column] = forward - pd.to_numeric(fits[reverse], errors="coerce") if reverse in fits.columns \
                else forward
        frames.append(pd.DataFrame(net))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _locked_pairs(stack: pd.DataFrame, min_abs_correlation: float) -> tuple[set[tuple[str, str]], set[str]]:
    """（锁死的列对, 有变化的列）；列对按字母序。"""
    pairs: set[tuple[str, str]] = set()
    varying: set[str] = set()
    tissues = sorted({str(c).split("--", 1)[0] for c in stack.columns})
    for tissue in tissues:
        block = stack[[c for c in stack.columns if str(c).split("--", 1)[0] == tissue]]
        block = block.dropna(axis=0, how="any")
        if len(block) < 3:
            continue
        spread = block.std(ddof=1)
        scale = block.abs().max().replace(0, 1.0)
        block = block.loc[:, spread > 1e-9 * scale]
        varying |= set(block.columns)
        if block.shape[1] < 2:
            continue
        values = block.to_numpy(dtype=float)
        z = (values - values.mean(axis=0)) / values.std(axis=0, ddof=1)
        corr = (z.T @ z) / (len(z) - 1)
        names = list(block.columns)
        rows, cols = np.nonzero(np.triu(np.abs(corr) >= min_abs_correlation, k=1))
        pairs |= {tuple(sorted((names[i], names[j]))) for i, j in zip(rows, cols)}
    return pairs, varying


def locked_groups(fit_sets: dict[str, dict[str, pd.DataFrame]], *,
                  min_abs_correlation: float = DEFAULT_MIN_ABS_CORRELATION, scope: str = "") -> pd.DataFrame:
    """``fit_sets``：模型（示踪剂）→ 样本 → 拟合表。返回每条组内反应一行（列见 :data:`COLUMNS`）。"""
    stacks = {name: net_fits(fits) for name, fits in fit_sets.items()}
    per_model = {name: _locked_pairs(stack, min_abs_correlation) for name, stack in stacks.items() if not stack.empty}
    candidates = set().union(*(pairs for pairs, _ in per_model.values())) if per_model else set()
    # 每个两者都有变化的模型中都须锁死
    accepted = {pair for pair in candidates
                if all(pair in pairs for pairs, varying in per_model.values() if set(pair) <= varying)}
    parent: dict[str, str] = {}

    def find(x: str) -> str:
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in sorted(accepted):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)
    groups: dict[str, list[str]] = {}
    for column in parent:
        groups.setdefault(find(column), []).append(column)
    reference = next((stack for stack in stacks.values() if not stack.empty), pd.DataFrame())
    rows = []
    for members in groups.values():
        members = sorted(members)
        lead = members[0]
        tissue = lead.split("--", 1)[0]
        for column in members:
            slope, intercept = _line(reference, column, lead)
            rows.append({"tissue": tissue, "reaction": column.split("--", 1)[1],
                         "flux_group": lead, "lead": lead.split("--", 1)[1], "n_members": len(members),
                         "members": ";".join(m.split("--", 1)[1] for m in members), "slope": slope,
                         "intercept": intercept, "scope": scope})
    table = pd.DataFrame(rows, columns=COLUMNS)
    return table.sort_values(["tissue", "flux_group", "reaction"]).reset_index(drop=True)


def _line(stack: pd.DataFrame, column: str, lead: str) -> tuple[float, float]:
    if column == lead:
        return 1.0, 0.0
    if column not in stack.columns or lead not in stack.columns:
        return np.nan, np.nan
    pair = stack[[lead, column]].dropna()
    if len(pair) < 3 or pair[lead].std() == 0:
        return np.nan, np.nan
    slope, intercept = np.polyfit(pair[lead].to_numpy(float), pair[column].to_numpy(float), 1)
    return float(slope), float(intercept)


def group_map(groups: pd.DataFrame) -> dict[tuple[str, str], str]:
    """（组织, 反应）→ 组标识。"""
    if groups is None or groups.empty:
        return {}
    return {(str(r.tissue), str(r.reaction)): str(r.flux_group) for r in groups.itertuples()}


def group_keys(tissues: pd.Series, reactions: pd.Series, groups: dict[tuple[str, str], str] | None) -> pd.Series:
    """每行的检验单位：组内反应为组标识，其余为 ``组织--反应``。"""
    groups = groups or {}
    return pd.Series([groups.get((str(t), str(r)), f"{t}--{r}") for t, r in zip(tissues, reactions)],
                     index=tissues.index)


def grouped_bh(pvalues: pd.Series, keys: pd.Series) -> pd.Series:
    """BH 校正，每个 ``keys`` 只计一次（组内取最小 p）；没有 p 值的行仍为 NaN。"""
    from statsmodels.stats.multitest import multipletests

    p = pd.to_numeric(pvalues, errors="coerce")
    valid = p.notna()
    out = pd.Series(np.nan, index=p.index)
    if not valid.any():
        return out
    representative = p[valid].groupby(keys[valid]).min()
    adjusted = pd.Series(multipletests(representative.to_numpy(float), method="fdr_bh")[1], index=representative.index)
    out[valid] = keys[valid].map(adjusted).to_numpy(float)
    return out


def n_tests(pvalues: pd.Series, keys: pd.Series) -> int:
    p = pd.to_numeric(pvalues, errors="coerce")
    return int(keys[p.notna()].nunique())
