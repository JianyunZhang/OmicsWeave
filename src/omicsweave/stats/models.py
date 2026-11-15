#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
``design.models`` 的公式模型（G3）：在每个组学、每个组织的数据集上拟合，报告合并主效应与交互项。

- 连续型数据层（蛋白、代谢物、代谢组派生特征、通量）用 :func:`~omicsweave.core.statistics.linear_model.fit_contrasts`
  （log2 尺度上的最小二乘加 limma 式方差压缩）。通量只用全部为正值的反应（取 log2）；动物来自不同示踪剂实验时，
  示踪剂（样本表的 ``protocol``）自动作为协变量；合并估计沿用两两比较中该反应的可辨识性与示踪剂一致性标记
  （任一比较判为不可辨识即记入），使证据画像的锚点门控同样适用；
- 转录组计数用 PyDESeq2：被报告的因素作为 ``condition``，模型中的其他主效应与协变量作为调整项；交互项不报告。

每个被报告的系数是一个“比较”：分类因素的主效应为 ``<模型>.<水平>_vs_<参照>``（``estimand_type = pooled_main``），
数值因素为 ``<模型>.<因素>``，交互项为 ``<模型>.<因素[水平]:因素[水平]>``（``interaction``）。合并估计与两两比较并列，
不互相替代；被处理逆转的指标，跨处理臂合并的饮食效应会偏小。设计矩阵不满秩或缺少样本时不拟合，原因写入
``qc/model_fits.parquet``。
"""
from __future__ import annotations

import contextlib
import io

import numpy as np
import pandas as pd

from omicsweave.contract.dataset import ValueType, read_dataset
from omicsweave.core.config import Experiment
from omicsweave.core.model import ModelSpec, Omics
from omicsweave.core.statistics.linear_model import DesignError, bh, fit_contrasts

from .data import DesignMatrix, analysis_values, build_design, bundle_datasets, sample_table, select_samples

ANNOTATION_COLUMNS = ("feature_name", "uniprot_id", "hmdb", "kegg_compound", "chebi")
DERIVED_COLUMNS = ("feature_type", "definition", "related_reaction", "readout", "n_members")
FIT_COLUMNS = ["model_id", "omics", "kind", "tissue", "n_samples", "formula", "status", "reason"]
# 区分示踪剂的样本表列（与 MFA 运行器相同）
TRACER_COLUMN = "protocol"


def _contrast_id(model: ModelSpec, design: DesignMatrix, label: str) -> tuple[str, str, tuple[str, str] | None]:
    """(比较 ID, estimand_type, (分子水平, 参照水平) 或 None)。"""
    term, levels = design.coefficients[label]
    if len(term) == 1 and levels[0]:
        reference = design.references[term[0]]
        return f"{model.id}.{levels[0]}_vs_{reference}", "pooled_main", (levels[0], reference)
    if len(term) == 1:
        return f"{model.id}.{term[0]}", "pooled_main", None
    return f"{model.id}.{label}", "interaction", None


def _cells(table: pd.DataFrame, design: DesignMatrix, label: str) -> tuple[pd.Index, pd.Index]:
    """系数对应的两组样本：主效应为该水平与参照水平，交互项为各因素都取该水平的格子与都取参照水平的格子。"""
    term, levels = design.coefficients[label]
    if not all(levels):
        return table.index, table.index
    hit = np.logical_and.reduce([table[n].astype(str) == lv for n, lv in zip(term, levels)])
    ref = np.logical_and.reduce([table[n].astype(str) == design.references[n] for n in term])
    return table.index[hit], table.index[ref]


def _annotations(differential: pd.DataFrame | None, tissue: str) -> pd.DataFrame:
    if differential is None or differential.empty:
        return pd.DataFrame()
    rows = differential[differential["tissue"].astype(str) == tissue]
    columns = [c for c in ANNOTATION_COLUMNS if c in rows.columns]
    return rows.drop_duplicates("feature_id").set_index("feature_id")[columns] if columns else pd.DataFrame()


def _mfa_quality(differential: pd.DataFrame | None, tissue: str) -> pd.DataFrame:
    """MFA 两两比较中各反应的可辨识性与示踪剂一致性（任一比较判为不可辨识或示踪剂冲突即记入，偏保守）。"""
    if differential is None or differential.empty:
        return pd.DataFrame(columns=["identifiability", "tracer_evidence", "quality_flags"])
    rows = differential[differential["tissue"].astype(str) == tissue]
    flags = rows["flags"].fillna("").astype(str) if "flags" in rows.columns else pd.Series("", index=rows.index)
    tracer = rows["tracer_evidence"].fillna("").astype(str) if "tracer_evidence" in rows.columns \
        else pd.Series("", index=rows.index)
    out = {}
    for feature, part in rows.groupby("feature_id", sort=False):
        tokens = set()
        if flags[part.index].str.split(";").map(lambda f: "not_identifiable" in f).any():
            tokens.add("not_identifiable")
        if tracer[part.index].eq("multi_tracer_conflict").any():
            tokens.add("tracer_conflict")
        identifiability = part["identifiability"].dropna().astype(str) if "identifiability" in part.columns \
            else pd.Series(dtype=str)
        out[str(feature)] = {"identifiability": identifiability.iloc[0] if not identifiability.empty else "",
                             "tracer_evidence": "multi_tracer_conflict" if "tracer_conflict" in tokens else "",
                             "quality_flags": ";".join(sorted(tokens))}
    return pd.DataFrame.from_dict(out, orient="index")


def fit_linear(model: ModelSpec, design: DesignMatrix, table: pd.DataFrame, dataset, omics: str, tissue: str,
               alpha: float, annotations: pd.DataFrame, differential: pd.DataFrame | None = None) -> pd.DataFrame:
    values = analysis_values(dataset)[design.samples]
    raw = dataset.matrix[design.samples].astype(float)
    flags = pd.Series("", index=values.index)
    if dataset.meta.value_type is ValueType.FLUX:
        nonpositive = (raw <= 0).any(axis=1)
        values.loc[nonpositive] = np.nan
        flags[nonpositive] = "nonpositive_flux"
    fits = fit_contrasts(values, design.x, {label: design.contrast(label) for label in design.coefficients})
    parts = []
    for label, fit in fits.items():
        contrast_id, estimand, pair = _contrast_id(model, design, label)
        numerator, denominator = _cells(table, design, label)
        t = fit.table
        tested = t["pvalue"].notna()
        out = pd.DataFrame({
            "contrast_id": contrast_id, "contrast_role": model.role.value, "tissue": tissue, "omics": omics,
            "feature_id": t.index.astype(str), "log2fc": t["estimate"], "lfc_se": t["se"],
            "effect_size": t["estimate"] / t["residual_sd"].replace(0, np.nan),
            "effect_size_type": "standardized_coefficient", "statistic": t["statistic"], "pvalue": t["pvalue"],
            "padj": bh(t["pvalue"]), "padj_method": np.where(tested, "fdr_bh", "none"),
            # 数值因素的斜率没有两组之分：分子记为全部样本，分母记为 0
            "n_numerator": len(numerator),
            "n_denominator": len(denominator) if estimand == "interaction" or pair else 0,
            "n_total": len(table), "method": "moderated_lm",
            "mean_numerator": raw[list(numerator)].mean(axis=1), "mean_denominator": raw[list(denominator)].mean(axis=1),
            "units": dataset.meta.units or ("log2" if dataset.meta.scale == "log2" else dataset.meta.value_type.value),
            "inference_unit": "biological_replicate",
            "evidence_level": np.where(tested, "statistical_test", "insufficient"),
            "estimand_type": estimand, "design_formula": design.formula, "flags": flags.to_numpy(),
        })
        out["selected"] = out["padj"].lt(alpha)
        out["selection_rule"] = f"padj < {alpha:g}"
        if not annotations.empty:
            for column in annotations.columns:
                out[column] = out["feature_id"].map(annotations[column]).to_numpy()
        for column in DERIVED_COLUMNS:
            if column in dataset.features.columns:
                out[column] = out["feature_id"].map(dataset.features[column]).to_numpy()
        if omics == Omics.MFA.value:
            # MFA 行必须给出两组均值；某一水平没有观测的反应无从估计，不写出
            out = out[out["mean_numerator"].notna() & out["mean_denominator"].notna()].copy()
            # 沿用两两比较中该反应的可辨识性与示踪剂一致性，使合并估计同样受锚点门控（G5）
            quality = _mfa_quality(differential, tissue)
            if not quality.empty:
                for column in ("identifiability", "tracer_evidence"):
                    out[column] = out["feature_id"].map(quality[column]).fillna("").to_numpy()
                extra = out["feature_id"].map(quality["quality_flags"]).fillna("")
                out["flags"] = [";".join(p for p in (f, e) if p) for f, e in zip(out["flags"], extra)]
        parts.append(out.reset_index(drop=True))
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def fit_counts(model: ModelSpec, design: DesignMatrix, table: pd.DataFrame, dataset, tissue: str,
               covariates: list[str], notes: list[str]) -> pd.DataFrame:
    from omicsweave.transcriptomics.export import de_to_contract
    from omicsweave.transcriptomics.rnaseq_pipeline import run_pydeseq2_analysis

    counts = dataset.matrix[design.samples].round().astype("int64")
    main = [t[0] for t in model.terms() if len(t) == 1]
    if any(len(t) > 1 for t in model.reported()):
        notes.append(f"模型 {model.id}：转录组不报告交互项（PyDESeq2 的对比只支持主效应）")
    parts = []
    for label, (term, levels) in design.coefficients.items():
        if len(term) > 1 or not levels[0]:
            continue
        factor = term[0]
        others = [n for n in dict.fromkeys(main + covariates) if n != factor and table[n].astype(str).nunique() > 1]
        meta = pd.DataFrame({"condition": table[factor].astype(str)}, index=table.index)
        for name in others:
            meta[name] = table[name].astype(str)
        contrast_id, estimand, (numerator, reference) = _contrast_id(model, design, label)
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                result = run_pydeseq2_analysis(counts, meta, design_formula=" + ".join([*others, "condition"]),
                                               contrast=("condition", numerator, reference))
        except Exception as exc:  # noqa: BLE001 — PyDESeq2 的拟合失败写入提示，不中断其他模型
            notes.append(f"模型 {model.id} 在转录组 {tissue} 的 PyDESeq2 拟合失败：{exc}")
            continue
        if result is None:
            continue
        out = de_to_contract(result, contrast_id=contrast_id, contrast_role=model.role.value, tissue=tissue)
        out["estimand_type"] = estimand
        out["n_total"] = len(table)
        parts.append(out)
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def fit_models(experiment: Experiment, run_dir, differential: dict[str, pd.DataFrame], alpha: float
               ) -> tuple[dict[str, list[pd.DataFrame]], dict[str, list[pd.DataFrame]], pd.DataFrame, list[str]]:
    """返回 (测得特征的模型结果 {组学: [表]}, 派生特征的模型结果, 拟合记录, 提示)。"""
    design_cfg = experiment.config.design
    samples = sample_table(experiment)
    measured: dict[str, list[pd.DataFrame]] = {}
    derived: dict[str, list[pd.DataFrame]] = {}
    fits, notes = [], []
    for model in design_cfg.models:
        wanted = {o.value for o in model.omics}     # 为空时为结果包中的全部组学
        needed = list(dict.fromkeys([n for t in model.terms() for n in t] + design_cfg.covariates))
        for omics, kind, tissue, path in bundle_datasets(run_dir):
            if wanted and omics not in wanted:
                continue
            dataset = read_dataset(path)
            table = samples.loc[[s for s in dataset.matrix.columns if s in samples.index]]
            row = {"model_id": model.id, "omics": omics, "kind": kind, "tissue": tissue, "n_samples": 0,
                   "formula": model.formula, "status": "fit", "reason": ""}
            try:
                table = select_samples(table, model.subset, design_cfg.factors, needed)
                covariates = [c for c in design_cfg.covariates if c not in {n for t in model.terms() for n in t}]
                # MFA：动物来自不同示踪剂实验时，示踪剂（samples.tsv 的 protocol）作为协变量
                if omics == Omics.MFA.value and TRACER_COLUMN in table.columns and TRACER_COLUMN not in covariates \
                        and table[TRACER_COLUMN].astype(str).str.strip().replace("", np.nan).dropna().nunique() > 1:
                    covariates.append(TRACER_COLUMN)
                design = build_design(table, model.terms(), design_cfg.factors, covariates)
            except (DesignError, KeyError, IndexError) as exc:
                row.update(status="not_fit", reason=str(exc), n_samples=len(table))
                fits.append(row)
                continue
            reported = {":".join(t) for t in model.reported()}
            design.coefficients = {k: v for k, v in design.coefficients.items() if ":".join(v[0]) in reported}
            row.update(n_samples=len(table), formula=design.formula)
            if omics == Omics.TRANSCRIPTOMICS.value and dataset.meta.value_type is ValueType.COUNTS:
                result = fit_counts(model, design, table, dataset, tissue, covariates, notes)
            else:
                result = fit_linear(model, design, table, dataset, omics, tissue, alpha,
                                    _annotations(differential.get(omics), tissue), differential.get(omics))
            if result.empty:
                row.update(status="not_fit", reason=row["reason"] or "没有可报告的系数")
            else:
                (derived if kind == "derived" else measured).setdefault(omics, []).append(result)
            fits.append(row)
    return measured, derived, pd.DataFrame(fits, columns=FIT_COLUMNS), list(dict.fromkeys(notes))
