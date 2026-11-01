#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MFA 差异通量结果 → 统一差异结果表（``omicsweave.contract.tables.DIFFERENTIAL_COLUMNS``）。

:func:`omicsweave.mfa.df_analysis.differential_flux` 每个比较输出一张宽表，组名写在列名里
（``mean_<组>``、``n_bio_<组>``），并同时给出 Welch t 与 Mann-Whitney 两条 p 值路径；实际采用哪条
由 ``primary_pvalue_column`` 决定，样本量不足以达到显著时为空（只给效应量）。转换后：

- ``pvalue`` / ``padj`` 取主路径的 p 值与 FDR；两条路径的原始 p 值另作附加列保留；
- ``method`` 取 ``primary_statistic``（如 ``ttest``、``effect_size_only``）；
- ``selected`` / ``selection_rule`` 取 ``is_significant`` / ``significance_rule``；
- ``flags`` 汇总：生物学重复不足（``low_n``）、效应量无定义（``effect_size_undefined``，如两组方差均为 0）、
  秩检验在该样本量下无法显著（``underpowered``）、通量不可辨识（``identifiability:<分档>``）。
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from omicsweave.contract.dataset import ContractError
from omicsweave.core.model import Omics

from .df_analysis import UNCONSTRAINED_LABELS

_PASSTHROUGH = ("identifiability", "frac_at_bound", "selection_probability", "direction_stability",
                "permutation_pvalue")
_EXTRA = ("flux_id", "diff", "diff_se", "log2fc_se", "cohen_d", "pvalue_ttest", "pvalue_mannwhitney",
          "fdr_ttest", "fdr_mannwhitney", "kegg_reaction", "kegg_pathway", "reactome_id", "ec_number")


def _column(table: pd.DataFrame, name: str, default) -> pd.Series:
    return table[name] if name in table.columns else pd.Series(default, index=table.index)


def _group_column(table: pd.DataFrame, prefix: str, group: str) -> pd.Series:
    """``mean_<组>`` 等列；组名大小写不敏感（比较名为大写，列名沿用配置中的写法）。"""
    matches = [c for c in table.columns if str(c).lower() == f"{prefix}_{group}".lower()]
    if not matches:
        raise ContractError(f"差异通量表缺少列 {prefix}_{group}")
    return table[matches[0]]


def differential_flux_to_contract(table: pd.DataFrame, *, numerator: str, denominator: str, contrast_id: str,
                                  contrast_role: str, flux_units: str) -> pd.DataFrame:
    """把一个比较的差异通量表转换为统一差异结果表。``numerator`` / ``denominator`` 为比较的两组。"""
    for column in ("tissue", "reaction", "log2FC", "hedges_g", "evidence_level", "primary_statistic"):
        if column not in table.columns:
            raise ContractError(f"差异通量表缺少列 {column}")
    primary = _column(table, "primary_pvalue_column", "").fillna("").astype(str)
    fdr_name = primary.str.replace("pvalue_", "fdr_", regex=False)

    def pick(names: pd.Series) -> pd.Series:
        values = [table.at[i, n] if n and n in table.columns else np.nan for i, n in names.items()]
        return pd.to_numeric(pd.Series(values, index=table.index), errors="coerce")

    is_ttest = primary == "pvalue_ttest"
    statistic = pd.to_numeric(_column(table, "t_statistic", np.nan), errors="coerce").where(is_ttest)
    multiple = _column(table, "multiple_testing", "none").fillna("none").astype(str)

    flags = []
    for _, row in table.iterrows():
        marks = []
        if row["evidence_level"] == "insufficient":
            marks.append("low_n")
        elif not np.isfinite(pd.to_numeric(row["hedges_g"], errors="coerce")):
            # 例如两组都贴在约束边界、方差为 0：Hedges' g 为 0/0，无定义
            marks.append("effect_size_undefined")
        if bool(row.get("rank_test_underpowered", False)):
            marks.append("underpowered")
        if row.get("identifiability") in UNCONSTRAINED_LABELS:
            marks.append(f"identifiability:{row['identifiability']}")
        flags.append(";".join(marks))

    out = pd.DataFrame({
        "contrast_id": contrast_id,
        "contrast_role": contrast_role,
        "tissue": table["tissue"].astype(str),
        "omics": Omics.MFA.value,
        "feature_id": table["reaction"].astype(str),
        "log2fc": pd.to_numeric(table["log2FC"], errors="coerce"),
        "effect_size": pd.to_numeric(table["hedges_g"], errors="coerce"),
        "effect_size_type": _column(table, "effect_size_metric", "hedges_g"),
        "statistic": statistic,
        "pvalue": pick(primary),
        "padj": pick(fdr_name),
        "padj_method": multiple.where(primary != "", "none"),
        "n_numerator": pd.to_numeric(_group_column(table, "n_bio", numerator)).astype("int64"),
        "n_denominator": pd.to_numeric(_group_column(table, "n_bio", denominator)).astype("int64"),
        "min_achievable_p": pd.to_numeric(_column(table, "min_achievable_pvalue", np.nan), errors="coerce"),
        "method": table["primary_statistic"].astype(str),
        "mean_numerator": pd.to_numeric(_group_column(table, "mean", numerator), errors="coerce"),
        "mean_denominator": pd.to_numeric(_group_column(table, "mean", denominator), errors="coerce"),
        "units": flux_units,
        "inference_unit": _column(table, "inference_unit", "biological_replicate"),
        "evidence_level": table["evidence_level"].astype(str),
        "selected": _column(table, "is_significant", False).fillna(False).astype(bool),
        "selection_rule": _column(table, "significance_rule", "").fillna("").astype(str),
        "flags": flags,
    }, index=table.index)
    for column in (*_PASSTHROUGH, *_EXTRA):
        if column in table.columns:
            out[column] = table[column]
    return out.reset_index(drop=True)
