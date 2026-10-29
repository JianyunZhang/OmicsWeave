#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
统一结果表的列规范与读写（v2 设计文档第 8.2 节）。

- ``differential.parquet``：所有组学的差异结果。MFA 的行另有必填列（效应量、组均值、单位、推断单位、证据等级），
  使“以生物学重复为推断单位”“样本量不足时只给效应量”等约束在结果表中可见；
- ``enrichment.parquet``：比较层面的通路富集结果，ORA 与预排序 GSEA 共用一张表，以 ``method`` 区分。
  样本层面的 ssGSEA 得分属于数据集（行为通路、列为样本），不放在这里。

表中允许出现规范以外的列（例如保留原工具的输出列以便追溯），写出时排在规范列之后。
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Literal

import numpy as np
import pandas as pd

from omicsweave.core.model import ContrastRole, Omics

from .dataset import ContractError

Kind = Literal["str", "float", "int", "bool"]


@dataclass(frozen=True)
class Column:
    name: str
    kind: Kind
    required: bool
    description: str


EVIDENCE_LEVELS = ("statistical_test", "effect_size_only", "insufficient")
EFFECT_DIRECTIONS = ("up", "down", "mixed", "none")
ENRICHMENT_METHODS = ("ora", "gsea_prerank")
ORA_QUERIES = {"up": "up", "down": "down", "all": "mixed"}  # 查询列表 -> effect_direction

DIFFERENTIAL_COLUMNS: tuple[Column, ...] = (
    Column("contrast_id", "str", True, "比较 ID，对应 experiment.yaml"),
    Column("contrast_role", "str", True, "disease_effect / treatment_effect / time_effect / other"),
    Column("tissue", "str", True, "组织"),
    Column("omics", "str", True, "组学"),
    Column("feature_id", "str", True, "特征 ID，类型见 dataset.json；MFA 为反应 ID"),
    Column("feature_name", "str", False, "特征名称"),
    Column("log2fc", "float", True, "分子组相对分母组的 log2 倍数变化；均值接近 0 或跨 0 时可为空"),
    Column("effect_size", "float", False, "效应量；MFA 必填"),
    Column("effect_size_type", "str", False, "如 hedges_g"),
    Column("statistic", "float", False, "检验统计量"),
    Column("pvalue", "float", True, "主检验的 p 值；只给效应量时为空"),
    Column("padj", "float", True, "主检验的校正后 p 值"),
    Column("padj_method", "str", True, "校正方法，如 fdr_bh；未校正时为 none"),
    Column("n_numerator", "int", True, "分子组实际参与检验的生物学重复数"),
    Column("n_denominator", "int", True, "分母组实际参与检验的生物学重复数"),
    Column("min_achievable_p", "float", False, "秩检验在该样本量下可达到的最小 p 值"),
    Column("method", "str", True, "检验方法或证据路径，如 deseq2、limma、ttest、effect_size_only"),
    Column("mean_numerator", "float", False, "分子组均值；MFA 必填，单位见 units"),
    Column("mean_denominator", "float", False, "分母组均值；MFA 必填"),
    Column("units", "str", False, "均值的单位；MFA 必填，如 nmol/min/g tissue"),
    Column("inference_unit", "str", False, "推断单位，如 biological_replicate；MFA 必填"),
    Column("evidence_level", "str", False, "statistical_test / effect_size_only / insufficient；MFA 必填"),
    Column("selected", "bool", False, "是否按 selection_rule 入选"),
    Column("selection_rule", "str", False, "入选规则的文字描述"),
    Column("identifiability", "str", False, "MFA：通量可辨识性分档，如 well_constrained / bound_pinned"),
    Column("frac_at_bound", "float", False, "MFA：拟合贴在约束边界的动物比例"),
    Column("selection_probability", "float", False, "MFA：生物学稳定性校准的入选概率"),
    Column("direction_stability", "float", False, "MFA：重抽样下方向一致的比例"),
    Column("permutation_pvalue", "float", False, "MFA：置换检验 p 值"),
    Column("flags", "str", False, "以分号分隔的标记，如 low_n;underpowered;effect_size_undefined"),
)
DIFFERENTIAL_KEY = ("contrast_id", "tissue", "omics", "feature_id")
# 只对该组学的行生效的必填列
DIFFERENTIAL_REQUIRED_BY_OMICS: dict[str, tuple[str, ...]] = {
    Omics.MFA.value: ("effect_size", "effect_size_type", "mean_numerator", "mean_denominator",
                      "units", "inference_unit", "evidence_level"),
}

ENRICHMENT_COLUMNS: tuple[Column, ...] = (
    Column("contrast_id", "str", True, "比较 ID"),
    Column("contrast_role", "str", True, "比较角色"),
    Column("tissue", "str", True, "组织"),
    Column("omics", "str", True, "组学"),
    Column("method", "str", True, "ora / gsea_prerank"),
    Column("query", "str", True, "ORA：up / down / all（所用特征列表）；GSEA：ranked"),
    Column("gene_set_library", "str", True, "基因集 / 通路库，如 KEGG_2019_Mouse"),
    Column("term_id", "str", False, "通路或条目 ID，如 ko00010、GO:0006096"),
    Column("term_name", "str", True, "通路或条目名称"),
    Column("pvalue", "float", True, "p 值"),
    Column("padj", "float", True, "校正后 p 值"),
    Column("padj_method", "str", True, "如 fdr_bh；GSEA 自带的 FDR q 值记为 gsea_fdr"),
    Column("effect_direction", "str", True, "up / down / mixed / none：ORA 由查询列表决定，GSEA 由 NES 符号决定"),
    Column("overlap_size", "int", False, "ORA：命中数；GSEA：leading edge 数。ORA 必填"),
    Column("set_size", "int", False, "通路大小（在背景中的成员数）。ORA 必填"),
    Column("query_size", "int", False, "ORA：查询列表大小"),
    Column("universe_size", "int", False, "ORA：背景大小"),
    Column("odds_ratio", "float", False, "ORA：比值比"),
    Column("es", "float", False, "GSEA：富集分数"),
    Column("nes", "float", False, "GSEA：标准化富集分数；GSEA 必填"),
    Column("hit_features", "str", False, "以分号分隔的命中特征（ORA 重叠成员 / GSEA leading edge）"),
    Column("tool", "str", False, "产生结果的工具，如 gseapy.prerank"),
    Column("flags", "str", False, "以分号分隔的标记"),
)
ENRICHMENT_KEY = ("contrast_id", "tissue", "omics", "method", "query", "gene_set_library", "term_name")


def _structure(table: pd.DataFrame, spec: tuple[Column, ...]) -> list[str]:
    """必填列与列类型；有问题时后续的取值检查不再进行。"""
    missing = [c.name for c in spec if c.required and c.name not in table.columns]
    if missing:
        return [f"缺少必填列：{', '.join(missing)}"]
    problems = []
    for column in spec:
        if column.name not in table.columns:
            continue
        series = table[column.name]
        if column.kind in ("float", "int") and not pd.api.types.is_numeric_dtype(series):
            problems.append(f"列 {column.name} 应为数值")
        elif column.kind == "bool" and not pd.api.types.is_bool_dtype(series):
            problems.append(f"列 {column.name} 应为布尔值")
    return problems


def _vocabulary(table: pd.DataFrame, column: str, allowed) -> list[str]:
    if column not in table.columns:
        return []
    unknown = set(table[column].dropna()) - set(allowed)
    return [f"{column} 取值无法识别：{', '.join(sorted(map(str, unknown)))}"] if unknown else []


def _probabilities(table: pd.DataFrame, columns) -> list[str]:
    return [f"列 {name} 有超出 [0, 1] 的值" for name in columns
            if name in table.columns and ((table[name] < 0) | (table[name] > 1)).any()]


def _duplicates(table: pd.DataFrame, key) -> list[str]:
    n = int(table.duplicated(list(key)).sum())
    return [f"有 {n} 行在 {' / '.join(key)} 上重复"] if n else []


def _blank(series: pd.Series) -> pd.Series:
    return series.isna() | (series.astype(str).str.strip() == "")


def _column_or_blank(table: pd.DataFrame, name: str) -> pd.Series:
    return table[name].fillna("").astype(str) if name in table.columns else pd.Series("", index=table.index)


def validate_differential(table: pd.DataFrame) -> list[str]:
    """返回差异结果表不符合规范的地方；空列表表示通过。"""
    problems = _structure(table, DIFFERENTIAL_COLUMNS)
    if problems:
        return problems
    problems += _probabilities(table, ("pvalue", "padj", "min_achievable_p"))
    for name in ("n_numerator", "n_denominator"):
        if (table[name] < 0).any():
            problems.append(f"列 {name} 有负值")
    problems += _vocabulary(table, "contrast_role", [r.value for r in ContrastRole])
    problems += _vocabulary(table, "omics", [o.value for o in Omics])
    problems += _vocabulary(table, "evidence_level", EVIDENCE_LEVELS)

    for omics, columns in DIFFERENTIAL_REQUIRED_BY_OMICS.items():
        rows = table[table["omics"] == omics]
        if rows.empty:
            continue
        missing = [c for c in columns if c not in rows.columns]
        if missing:
            problems.append(f"{omics} 的行缺少必填列：{', '.join(missing)}")
            continue
        for column in columns:
            blank = _blank(rows[column])
            if column == "effect_size":
                # 效应量无定义时允许为空：生物学重复不足，或两组方差均为 0（须标记 effect_size_undefined）
                undefined = _column_or_blank(rows, "flags").str.contains("effect_size_undefined", regex=False)
                blank &= (rows["evidence_level"] != "insufficient") & ~undefined
            if blank.any():
                problems.append(f"{omics} 有 {int(blank.sum())} 行的 {column} 为空")
    return problems + _duplicates(table, DIFFERENTIAL_KEY)


def validate_enrichment(table: pd.DataFrame) -> list[str]:
    """返回富集结果表不符合规范的地方；空列表表示通过。"""
    problems = _structure(table, ENRICHMENT_COLUMNS)
    if problems:
        return problems
    problems += _probabilities(table, ("pvalue", "padj"))
    problems += _vocabulary(table, "contrast_role", [r.value for r in ContrastRole])
    problems += _vocabulary(table, "omics", [o.value for o in Omics])
    problems += _vocabulary(table, "method", ENRICHMENT_METHODS)
    problems += _vocabulary(table, "effect_direction", EFFECT_DIRECTIONS)
    blank_terms = _blank(table["term_name"])
    if blank_terms.any():
        problems.append(f"有 {int(blank_terms.sum())} 行的 term_name 为空")

    ora = table[table["method"] == "ora"]
    if not ora.empty:
        problems += [f"ORA 的 query 取值无法识别：{q}（可选：{' / '.join(ORA_QUERIES)}）"
                     for q in sorted(set(ora["query"]) - set(ORA_QUERIES))]
        for column in ("overlap_size", "set_size"):
            if column not in ora.columns or ora[column].isna().any():
                problems.append(f"ORA 的行需要 {column}")
        if {"overlap_size", "set_size"} <= set(ora.columns) and (ora["overlap_size"] > ora["set_size"]).any():
            problems.append("ORA 有行的 overlap_size 大于 set_size")
        expected = ora["query"].map(ORA_QUERIES)
        if (expected.notna() & (expected != ora["effect_direction"])).any():
            problems.append("ORA 的 effect_direction 应与 query 对应（up→up，down→down，all→mixed）")

    gsea = table[table["method"] == "gsea_prerank"]
    if not gsea.empty:
        if (gsea["query"] != "ranked").any():
            problems.append("GSEA 的 query 应为 ranked")
        if "nes" not in gsea.columns or gsea["nes"].isna().any():
            problems.append("GSEA 的行需要 nes")
        else:
            expected = np.select([gsea["nes"] > 0, gsea["nes"] < 0], ["up", "down"], "none")
            if (expected != gsea["effect_direction"].to_numpy()).any():
                problems.append("GSEA 的 effect_direction 应与 NES 的符号一致")
    return problems + _duplicates(table, ENRICHMENT_KEY)


def _write(table: pd.DataFrame, path: str | Path, spec: tuple[Column, ...],
           validate: Callable[[pd.DataFrame], list[str]]) -> Path:
    problems = validate(table)
    if problems:
        raise ContractError("；".join(problems))
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    ordered = [c.name for c in spec if c.name in table.columns]
    extra = [c for c in table.columns if c not in ordered]
    table[ordered + extra].to_parquet(path, index=False)
    return path


def _read(path: str | Path, validate: Callable[[pd.DataFrame], list[str]]) -> pd.DataFrame:
    table = pd.read_parquet(path)
    problems = validate(table)
    if problems:
        raise ContractError(f"{path}: " + "；".join(problems))
    return table


def write_differential(table: pd.DataFrame, path: str | Path) -> Path:
    return _write(table, path, DIFFERENTIAL_COLUMNS, validate_differential)


def read_differential(path: str | Path) -> pd.DataFrame:
    return _read(path, validate_differential)


def write_enrichment(table: pd.DataFrame, path: str | Path) -> Path:
    return _write(table, path, ENRICHMENT_COLUMNS, validate_enrichment)


def read_enrichment(path: str | Path) -> pd.DataFrame:
    return _read(path, validate_enrichment)
