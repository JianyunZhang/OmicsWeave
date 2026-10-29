#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
统一结果表的列规范与读写（v2 设计文档第 8.2 节）。

目前只定义差异结果表 ``differential.parquet``；富集与 MFA 通量结果表在阶段 2 补充。
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import pandas as pd

from omicsweave.core.model import ContrastRole, Omics

from .dataset import ContractError

Kind = Literal["str", "float", "int"]


@dataclass(frozen=True)
class Column:
    name: str
    kind: Kind
    required: bool
    description: str


DIFFERENTIAL_COLUMNS: tuple[Column, ...] = (
    Column("contrast_id", "str", True, "比较 ID，对应 experiment.yaml"),
    Column("contrast_role", "str", True, "disease_effect / treatment_effect / time_effect / other"),
    Column("tissue", "str", True, "组织"),
    Column("omics", "str", True, "组学"),
    Column("feature_id", "str", True, "特征 ID，类型见 dataset.json"),
    Column("feature_name", "str", False, "特征名称"),
    Column("log2fc", "float", True, "分子组相对分母组的 log2 倍数变化"),
    Column("effect_size", "float", False, "效应量，MFA 必填"),
    Column("effect_size_type", "str", False, "如 hedges_g"),
    Column("statistic", "float", False, "检验统计量"),
    Column("pvalue", "float", True, "原始 p 值"),
    Column("padj", "float", True, "校正后 p 值"),
    Column("padj_method", "str", True, "校正方法，如 fdr_bh"),
    Column("n_numerator", "int", True, "分子组实际参与检验的生物学重复数"),
    Column("n_denominator", "int", True, "分母组实际参与检验的生物学重复数"),
    Column("min_achievable_p", "float", False, "秩检验在该样本量下可达到的最小 p 值"),
    Column("method", "str", True, "检验方法"),
    Column("flags", "str", False, "以分号分隔的标记，如 low_n;at_bound"),
)

_KEY = ("contrast_id", "tissue", "omics", "feature_id")


def validate_differential(table: pd.DataFrame) -> list[str]:
    """返回差异结果表不符合规范的地方；空列表表示通过。"""
    problems = []
    missing = [c.name for c in DIFFERENTIAL_COLUMNS if c.required and c.name not in table.columns]
    if missing:
        return [f"缺少必填列：{', '.join(missing)}"]
    for column in DIFFERENTIAL_COLUMNS:
        if column.name not in table.columns:
            continue
        series = table[column.name]
        if column.kind in ("float", "int") and not pd.api.types.is_numeric_dtype(series):
            problems.append(f"列 {column.name} 应为数值")
    if problems:
        return problems
    for name in ("pvalue", "padj", "min_achievable_p"):
        if name in table.columns and ((table[name] < 0) | (table[name] > 1)).any():
            problems.append(f"列 {name} 有超出 [0, 1] 的值")
    for name in ("n_numerator", "n_denominator"):
        if (table[name] < 0).any():
            problems.append(f"列 {name} 有负值")
    bad_roles = set(table["contrast_role"]) - {r.value for r in ContrastRole}
    if bad_roles:
        problems.append(f"contrast_role 取值无法识别：{', '.join(sorted(map(str, bad_roles)))}")
    bad_omics = set(table["omics"]) - {o.value for o in Omics}
    if bad_omics:
        problems.append(f"omics 取值无法识别：{', '.join(sorted(map(str, bad_omics)))}")
    duplicated = table.duplicated(list(_KEY)).sum()
    if duplicated:
        problems.append(f"有 {duplicated} 行在 {' / '.join(_KEY)} 上重复")
    return problems


def write_differential(table: pd.DataFrame, path: str | Path) -> Path:
    problems = validate_differential(table)
    if problems:
        raise ContractError("；".join(problems))
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    ordered = [c.name for c in DIFFERENTIAL_COLUMNS if c.name in table.columns]
    extra = [c for c in table.columns if c not in ordered]
    table[ordered + extra].to_parquet(path, index=False)
    return path


def read_differential(path: str | Path) -> pd.DataFrame:
    table = pd.read_parquet(path)
    problems = validate_differential(table)
    if problems:
        raise ContractError(f"{path}: " + "；".join(problems))
    return table
