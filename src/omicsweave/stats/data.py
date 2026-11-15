#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""结果包中的分析层数据集、样本信息与设计矩阵。"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from omicsweave.contract.dataset import META_FILE, Dataset, ValueType, read_dataset
from omicsweave.core.config import Experiment
from omicsweave.core.model import Contrast, Factor, Omics
from omicsweave.core.run import contrast_groups
from omicsweave.core.statistics.linear_model import check_rank
from omicsweave.core.validate import biological_samples, matches_level

# 数据集目录：data 为测得特征，derived 为派生特征
KINDS = ("data", "derived")


def bundle_datasets(run_dir: str | Path, kinds: tuple[str, ...] = KINDS) -> list[tuple[str, str, str, Path]]:
    """结果包中的数据集：(组学, 类别, 组织, 目录)。"""
    run_dir = Path(run_dir)
    out = []
    for kind in kinds:
        for meta in sorted(run_dir.glob(f"omics/*/{kind}/*/{META_FILE}")):
            out.append((meta.parent.parent.parent.name, kind, meta.parent.name, meta.parent))
    return out


def analysis_values(dataset: Dataset) -> pd.DataFrame:
    """统计用的 log2 尺度值：计数为 log2(CPM + 1)，线性尺度取 log2（非正值记为缺失），log2 尺度原样。"""
    matrix = dataset.matrix.astype(float)
    if dataset.meta.value_type is ValueType.COUNTS:
        library = matrix.sum(axis=0).replace(0, np.nan)
        return np.log2(matrix / library * 1e6 + 1.0)
    if dataset.meta.scale == "log2":
        return matrix
    return np.log2(matrix.where(matrix > 0))


def sample_table(experiment: Experiment) -> pd.DataFrame:
    """生物学样本（索引为 sample_id），含样本层与动物层的全部列。"""
    return biological_samples(experiment).drop_duplicates("sample_id").set_index("sample_id")


def pairwise_groups(experiment: Experiment, contrast: Contrast, omics: str, columns) -> dict[str, str]:
    """一个比较在数据集列中的样本 → 所在组。"""
    groups = contrast_groups(experiment, contrast, Omics(omics))
    return {s: groups[s] for s in columns if s in groups}


@dataclass
class DesignMatrix:
    x: np.ndarray
    columns: list[str]
    samples: list[str]
    # 报告的系数列 → (项, 各因素的水平)；截距与协变量列不在其中
    coefficients: dict[str, tuple[tuple[str, ...], tuple[str, ...]]] = field(default_factory=dict)
    references: dict[str, str] = field(default_factory=dict)
    formula: str = ""
    # 全部系数列（含协变量与未报告的项）→ (项, 水平)；分类因素 → 数据中出现的水平数
    column_terms: dict[str, tuple[tuple[str, ...], tuple[str, ...]]] = field(default_factory=dict)
    n_levels: dict[str, int] = field(default_factory=dict)

    def contrast(self, label: str) -> np.ndarray:
        """报告一个系数所用的对比向量。

        处理编码下，参与二阶交互的因素，其主效应系数是另一因素取参照水平时的效应；这里改为对另一因素的各水平
        等权平均（主效应系数 + 各交互系数 / 另一因素的水平数），使 pooled_main 的含义与无交互的模型一致。
        与数值因素的交互、三阶及以上的交互不做平均。
        """
        vector = np.zeros(len(self.columns))
        vector[self.columns.index(label)] = 1.0
        term, levels = self.column_terms[label]
        if len(term) == 1 and levels[0]:
            for other_label, (other_term, other_levels) in self.column_terms.items():
                if len(other_term) == 2 and term[0] in other_term:
                    position = other_term.index(term[0])
                    partner = other_term[1 - position]
                    if other_levels[position] == levels[0] and partner in self.n_levels:
                        vector[self.columns.index(other_label)] = 1.0 / self.n_levels[partner]
        return vector


def _levels(values: pd.Series, factor: Factor | None) -> tuple[list[str], str]:
    present = list(dict.fromkeys(values.astype(str)))
    if factor is not None and factor.levels:
        ordered = [lv for lv in factor.levels if lv in present]
    else:
        ordered = sorted(present)
    reference = factor.reference if factor is not None and factor.reference in ordered else ordered[0]
    return [reference] + [lv for lv in ordered if lv != reference], reference


def build_design(table: pd.DataFrame, terms: list[tuple[str, ...]], factors: dict[str, Factor],
                 covariates: list[str]) -> DesignMatrix:
    """``table``：样本 × 变量；``terms``：公式中的项。协变量作为主效应加入（只有一个取值的协变量略去）。

    分类因素按参照水平做处理编码（参照水平取 ``reference``，不在数据中时取第一个出现的水平）；数值因素按原值。
    设计矩阵不满秩时抛出 :class:`~omicsweave.core.statistics.linear_model.DesignError` 并说明原因。
    """
    all_terms = list(terms) + [(c,) for c in covariates if (c,) not in terms and table[c].astype(str).nunique() > 1]
    columns: dict[str, np.ndarray] = {"(Intercept)": np.ones(len(table))}
    dummies: dict[str, dict[str, np.ndarray]] = {}
    references: dict[str, str] = {}
    numeric: dict[str, np.ndarray] = {}
    for name in dict.fromkeys(n for t in all_terms for n in t):
        factor = factors.get(name)
        if factor is not None and factor.type == "numeric":
            numeric[name] = pd.to_numeric(table[name], errors="coerce").to_numpy(dtype=float)
            continue
        levels, reference = _levels(table[name], factor)
        references[name] = reference
        dummies[name] = {lv: (table[name].astype(str) == lv).to_numpy(dtype=float) for lv in levels[1:]}
    coefficients: dict[str, tuple[tuple[str, ...], tuple[str, ...]]] = {}
    column_terms: dict[str, tuple[tuple[str, ...], tuple[str, ...]]] = {}
    for term in all_terms:
        parts = [[("", numeric[n])] if n in numeric else list(dummies[n].items()) for n in term]
        combos = [[]]
        for options in parts:
            combos = [c + [o] for c in combos for o in options]
        for combo in combos:
            label = ":".join(f"{n}[{lv}]" if lv else n for n, (lv, _) in zip(term, combo))
            columns[label] = np.prod([v for _, v in combo], axis=0)
            column_terms[label] = (term, tuple(lv for lv, _ in combo))
            if term in terms:
                coefficients[label] = column_terms[label]
    x = np.column_stack(list(columns.values()))
    check_rank(x, list(columns))
    formula = "~ " + " + ".join(":".join(t) for t in all_terms)
    return DesignMatrix(x, list(columns), list(table.index), coefficients, references, formula, column_terms,
                        {name: len(levels) + 1 for name, levels in dummies.items()})


def select_samples(table: pd.DataFrame, subset: dict[str, str], factors: dict[str, Factor],
                   needed: list[str]) -> pd.DataFrame:
    """按 subset 筛选样本，并去掉模型变量缺失的样本。"""
    for key, value in subset.items():
        table = table[matches_level(table[key], value, factors[key].type == "numeric")]
    present = [c for c in needed if c in table.columns]
    missing = [c for c in needed if c not in table.columns]
    if missing:
        raise KeyError(f"样本表中没有列 {'、'.join(missing)}")
    blank = table[present].astype(str).apply(lambda s: s.str.strip().eq("") | s.eq("nan")).any(axis=1)
    return table[~blank]
