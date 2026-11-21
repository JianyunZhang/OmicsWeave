#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
逆转判定（``reversal``）：对 ``design.reversal_pairs`` 中的每个（疾病效应, 处理效应）配对，逐特征判断处理是否逆转了疾病效应。

读取源结果包中各组学的 ``differential.parquet`` 与代谢组派生特征的 ``derived_differential.parquet``，在同一组学、
同一组织中按 ``feature_id`` 对应两个比较。记疾病效应为 d、处理效应为 t（均为 log2 倍数变化），逆转指数
``reversal_index = -t / d``：1 为完全逆转，0 为没有作用，大于 1 为过度逆转，小于 0 为加重。区间按 delta 法
（假定两个比较的估计相互独立；两个比较共用同一组动物时区间偏保守或偏乐观，见文档），标准误取结果表的
``lfc_se`` / ``log2fc_se``，没有时由 p 值反推（``se_source`` 记为 ``from_pvalue``）。

分类（``class``，判定依据写在 ``rule``）：

- ``undetermined``（无法判定）：疾病效应本身未建立（不是统计检验，或 padj ≥ ``alpha``）、处理效应缺失或未检验、
  任一比较在可靠性表中被 QC 判为不可评估，或证据不足以下结论（``inconclusive``）；
- ``reversed`` / ``partially_reversed`` / ``over_reversed``：处理效应显著（padj < ``alpha``）且与疾病方向相反，
  逆转指数的区间包含 1 / 整体小于 1 / 整体大于 1；
- ``aggravated``（加重）：处理效应显著且与疾病方向相同；
- ``not_reversed``（未逆转）：处理效应不显著，并且通过等效性检验（双单侧检验，``1 - 2·alpha`` 区间落在
  ±``equivalence_margin`` × |d| 之内）。低功效的不显著结果归入 undetermined，不会被判为未逆转。

参数（``integration.parameters.reversal``）：``alpha``（默认取 ``project.yaml`` 的 ``defaults.alpha``）、
``equivalence_margin``（等效界限占疾病效应的比例，默认 0.5）、``ci_level``（逆转指数区间，默认 0.95）、
``qc_imbalance``（见 :mod:`omicsweave.core.reliability`）。

输出（``integration/reversal/``）：``tables/reversal.parquet`` 每个（配对, 组学, 组织, 特征）一行；
``tables/reversal_summary.parquet`` 按（配对, 组学, 组织, 特征类型）统计各类数目；``qc/contrast_reliability.parquet``。
通路层面的逆转与析因设计中的交互项需要通路效应的标准误与公式模型，尚未实现。
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import norm

from omicsweave.contract.tables import effect_standard_errors, read_differential
from omicsweave.core.model import IntegrationMode, Omics
from omicsweave.core.reliability import (
    BLOCKED, DEFAULT_IMBALANCE, QC_DIR, RELIABILITY_FILE, build_reliability, summary_note,
)
from omicsweave.core.run import RunContext, RunError
from omicsweave.core.tissues import remap_table_tissues

from ..integrate import IntegrationRun, mode_dir
from .flux_anchored import DERIVED_TABLE

MODE = IntegrationMode.REVERSAL
DEFAULTS: dict[str, Any] = {"alpha": None, "equivalence_margin": 0.5, "ci_level": 0.95,
                            "qc_imbalance": DEFAULT_IMBALANCE}
CLASSES = ("reversed", "partially_reversed", "over_reversed", "aggravated", "not_reversed", "undetermined")
COLUMNS = ["pair_id", "disease_contrast", "treatment_contrast", "omics", "tissue", "feature_id", "feature_name",
           "feature_type", "disease_effect", "disease_se", "disease_padj", "treatment_effect", "treatment_se",
           "treatment_padj", "reversal_index", "ci_lower", "ci_upper", "class", "rule", "se_source"]
standard_errors = effect_standard_errors


def _features(tables: dict[str, pd.DataFrame], derived: pd.DataFrame | None) -> pd.DataFrame:
    parts = [t.assign(feature_type="measured") for t in tables.values()]
    if derived is not None and not derived.empty:
        parts.append(derived.assign(feature_type=derived["feature_type"].fillna("derived").astype(str)))
    table = pd.concat(parts, ignore_index=True)
    if "feature_name" not in table.columns:
        table["feature_name"] = table["feature_id"]
    se, source = standard_errors(table)
    return table.assign(se=se, se_source=source)


def classify(merged: pd.DataFrame, *, alpha: float, equivalence_margin: float, ci_level: float,
             blocked: pd.Series | None = None) -> pd.DataFrame:
    """``merged``：每行一个特征，含 ``*_d``（疾病比较）与 ``*_t``（处理比较）的 log2fc、se、padj、evidence_level。"""
    d, t = merged["log2fc_d"].astype(float), merged["log2fc_t"].astype(float)
    se_d, se_t = merged["se_d"].astype(float), merged["se_t"].astype(float)
    padj_d, padj_t = merged["padj_d"].astype(float), merged["padj_t"].astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        ri = -t / d
        ri_se = np.sqrt(se_t ** 2 / d ** 2 + t ** 2 * se_d ** 2 / d ** 4)
    z = norm.isf((1 - ci_level) / 2)
    lower, upper = ri - z * ri_se, ri + z * ri_se
    ci_ok = np.isfinite(lower) & np.isfinite(upper)

    disease_ok = merged["evidence_level_d"].eq("statistical_test") & padj_d.lt(alpha) & np.isfinite(d) & d.ne(0)
    tested = merged["evidence_level_t"].eq("statistical_test") & np.isfinite(padj_t)
    significant = tested & padj_t.lt(alpha)
    opposite = np.sign(t) == -np.sign(d)
    # 双单侧检验：1 - 2·alpha 区间落在 ±margin 内
    equivalent = tested & np.isfinite(se_t) & ((t.abs() + norm.isf(alpha) * se_t) < equivalence_margin * d.abs())
    blocked = blocked if blocked is not None else pd.Series(False, index=merged.index)

    conditions = [
        blocked.to_numpy(dtype=bool),
        ~disease_ok,
        ~np.isfinite(t),
        significant & opposite & ci_ok & (lower > 1),
        significant & opposite & ci_ok & (upper < 1),
        significant & opposite & ci_ok,
        significant & opposite,
        significant & ~opposite,
        equivalent,
        ~tested,
    ]
    classes = ["undetermined", "undetermined", "undetermined", "over_reversed", "partially_reversed", "reversed",
               "reversed", "aggravated", "not_reversed", "undetermined"]
    rules = ["qc_not_assessable", "disease_effect_not_established", "treatment_effect_missing",
             "treatment_opposite;ri_ci_above_1", "treatment_opposite;ri_ci_below_1",
             "treatment_opposite;ri_ci_includes_1", "treatment_opposite;ri_ci_unavailable",
             "treatment_same_direction", "equivalence_tost", "treatment_not_tested"]
    conditions = [np.asarray(c, dtype=bool) for c in conditions]
    return pd.DataFrame({
        "reversal_index": ri.where(disease_ok & np.isfinite(t)), "ci_lower": lower.where(disease_ok & ci_ok),
        "ci_upper": upper.where(disease_ok & ci_ok),
        "class": np.select(conditions, classes, "undetermined"),
        "rule": np.select(conditions, rules, "inconclusive"),
    }, index=merged.index)


def reversal_table(features: pd.DataFrame, pairs: list[tuple[str, str]], *, alpha: float,
                   equivalence_margin: float, ci_level: float, blocked_keys: set[tuple[str, str, str]]) -> pd.DataFrame:
    keys = ["omics", "tissue", "feature_id"]
    keep = [*keys, "feature_name", "feature_type", "log2fc", "se", "padj", "evidence_level", "se_source"]
    parts = []
    for disease, treatment in pairs:
        d = features.loc[features["contrast_id"] == disease, keep]
        t = features.loc[features["contrast_id"] == treatment, keep]
        merged = d.merge(t.drop(columns=["feature_name", "feature_type"]), on=keys, suffixes=("_d", "_t"))
        if merged.empty:
            continue
        blocked = pd.Series([(o, ti, disease) in blocked_keys or (o, ti, treatment) in blocked_keys
                             for o, ti in zip(merged["omics"], merged["tissue"])], index=merged.index)
        result = classify(merged, alpha=alpha, equivalence_margin=equivalence_margin, ci_level=ci_level,
                          blocked=blocked)
        sources = (merged["se_source_d"].fillna("") + "/" + merged["se_source_t"].fillna("")).str.strip("/")
        parts.append(pd.DataFrame({
            "pair_id": f"{disease}__{treatment}", "disease_contrast": disease, "treatment_contrast": treatment,
            "omics": merged["omics"], "tissue": merged["tissue"], "feature_id": merged["feature_id"],
            "feature_name": merged["feature_name"], "feature_type": merged["feature_type"],
            "disease_effect": merged["log2fc_d"], "disease_se": merged["se_d"], "disease_padj": merged["padj_d"],
            "treatment_effect": merged["log2fc_t"], "treatment_se": merged["se_t"],
            "treatment_padj": merged["padj_t"], **{c: result[c] for c in result.columns}, "se_source": sources,
        }))
    if not parts:
        return pd.DataFrame(columns=COLUMNS)
    return pd.concat(parts, ignore_index=True)[COLUMNS]


def summary_table(table: pd.DataFrame) -> pd.DataFrame:
    keys = ["pair_id", "omics", "tissue", "feature_type"]
    if table.empty:
        return pd.DataFrame(columns=[*keys, "n_features", *CLASSES])
    counts = pd.crosstab([table[k] for k in keys], table["class"]).reindex(columns=list(CLASSES), fill_value=0)
    counts.insert(0, "n_features", counts.sum(axis=1))
    return counts.reset_index()


def run(ctx: RunContext, source: Path, parameters: dict[str, Any]) -> IntegrationRun:
    unknown = sorted(set(parameters) - set(DEFAULTS))
    if unknown:
        raise RunError(f"integration.parameters.{MODE.value} 中有无法识别的参数：{'、'.join(unknown)}"
                       f"（可用：{'、'.join(DEFAULTS)}）")
    params = {**DEFAULTS, **parameters}
    alpha = float(params["alpha"] if params["alpha"] is not None else ctx.project.config.defaults.alpha)
    pairs = [(p.disease, p.treatment) for p in ctx.experiment.config.design.reversal_pairs]
    if not pairs:
        raise RunError(f"{MODE.value} 需要在 experiment.yaml 的 design.reversal_pairs 中声明（疾病效应, 处理效应）配对")
    paths = {p.parent.parent.name: p for p in sorted(source.glob("omics/*/tables/differential.parquet"))}
    if not paths:
        raise RunError(f"结果包 {source.name} 中没有差异结果，请先运行 omicsweave run（含 differential 步骤）")
    tables = {omics: read_differential(path) for omics, path in paths.items()}
    derived_path = source / "omics" / Omics.METABOLOMICS.value / "tables" / DERIVED_TABLE
    derived = read_differential(derived_path) if derived_path.is_file() else None
    try:
        for omics in list(tables):
            aliases = ctx.experiment.config.tissue_aliases_for(omics)
            tables[omics], _ = remap_table_tissues(tables[omics], aliases, context=f"源结果包 {source.name} 的 {omics}")
        if derived is not None:
            derived, _ = remap_table_tissues(
                derived, ctx.experiment.config.tissue_aliases_for(Omics.METABOLOMICS),
                context=f"源结果包 {source.name} 的 metabolomics 派生特征",
            )
    except ValueError as exc:
        raise RunError(str(exc)) from exc

    reliability, quality_inputs = build_reliability(ctx.experiment, source, tables, float(params["qc_imbalance"]))
    blocked_keys = {(o, t, c) for o, t, c in reliability.loc[reliability["status"].isin(BLOCKED),
                                                              ["omics", "tissue", "contrast_id"]].itertuples(index=False)}
    table = reversal_table(_features(tables, derived), pairs, alpha=alpha,
                           equivalence_margin=float(params["equivalence_margin"]),
                           ci_level=float(params["ci_level"]), blocked_keys=blocked_keys)

    base = mode_dir(ctx, MODE)
    (base / "tables").mkdir(parents=True)
    (base / QC_DIR).mkdir(parents=True)
    table.to_parquet(base / "tables" / "reversal.parquet", index=False)
    summary_table(table).to_parquet(base / "tables" / "reversal_summary.parquet", index=False)
    reliability.to_parquet(base / QC_DIR / RELIABILITY_FILE, index=False)

    record = IntegrationRun(MODE, inputs=list(paths.values()) + ([derived_path] if derived is not None else [])
                            + quality_inputs,
                            parameters={**params, "alpha": alpha, "pairs": [list(p) for p in pairs],
                                        "omics": sorted(paths)})
    found = set(table["pair_id"])
    for disease, treatment in pairs:
        if f"{disease}__{treatment}" not in found:
            record.notes.append(f"逆转配对 {disease} / {treatment} 在任何组学中都没有两个比较共有的特征")
    note = summary_note(reliability)
    if note:
        record.notes.append(note)
    return record
