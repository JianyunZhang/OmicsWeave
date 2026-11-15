#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
入选特征的稳健性（E5，``omics/<组学>/tables/robustness.parquet``）。

- 逐只去掉动物：对每个两两比较中入选（``selected``）的特征，在 log2 尺度上（计数为 log2 CPM）用组均值之差估计效应，
  依次去掉一个样本重新估计。去掉某一只后效应反号，或缩小到原来的 ``min_ratio``（默认 0.5）以下，记为
  ``single_subject_driven``，并给出影响最大的样本；
- 符号错误率与夸大倍数（Gelman 与 Carlin 2014）：以正态—正态收缩后的效应作为可信的真实效应，
  结合差异结果的标准误，估计该比较下的功效、Type S 与 Type M。

去掉的是样本（每只动物在一个组织中一个样本）。只剩 1 个样本的组不能再去掉，``n_left_out`` 记录实际重算的次数。
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from omicsweave.contract.dataset import read_dataset
from omicsweave.contract.tables import effect_standard_errors
from omicsweave.core.config import Experiment
from omicsweave.core.statistics.shrinkage import normal_shrinkage, retrodesign

from .data import analysis_values, bundle_datasets, pairwise_groups

FLAG = "single_subject_driven"
DEFAULT_MIN_RATIO = 0.5
COLUMNS = ["contrast_id", "tissue", "omics", "feature_id", "n_numerator", "n_denominator", "n_left_out",
           "full_effect", "loo_min_ratio", "loo_sign_flips", "driver_sample", FLAG, "se", "plausible_effect",
           "power", "type_s", "exaggeration"]


def leave_one_out(first: pd.DataFrame, second: pd.DataFrame) -> tuple[pd.Series, pd.DataFrame]:
    """(全样本的组均值之差, 去掉每个样本后的差：特征 × 被去掉的样本)。"""
    def loo_means(values: pd.DataFrame) -> pd.DataFrame:
        total, count = values.sum(axis=1, min_count=1), values.notna().sum(axis=1)
        out = {}
        for sample in values.columns:
            present = values[sample].notna()
            remaining = count - present.astype(int)
            out[sample] = ((total - values[sample].fillna(0.0)) / remaining.where(remaining > 0))
        return pd.DataFrame(out, index=values.index)

    mean1, mean2 = first.mean(axis=1), second.mean(axis=1)
    loo = pd.concat([loo_means(first).sub(mean2, axis=0), -loo_means(second).sub(mean1, axis=0)], axis=1)
    return mean1 - mean2, loo


def robustness(experiment: Experiment, run_dir, tables: dict[str, pd.DataFrame],
               derived: dict[str, pd.DataFrame], *, alpha: float = 0.05,
               min_ratio: float = DEFAULT_MIN_RATIO) -> dict[str, pd.DataFrame]:
    """返回 {组学: 稳健性表}；``tables`` / ``derived`` 为测得特征与派生特征的差异结果。"""
    contrasts = {c.id: c for c in experiment.config.design.contrasts}
    out: dict[str, list[pd.DataFrame]] = {}
    for omics, kind, tissue, path in bundle_datasets(run_dir):
        source = (derived if kind == "derived" else tables).get(omics)
        if source is None or source.empty:
            continue
        values = analysis_values(read_dataset(path))
        rows_tissue = source[source["tissue"].astype(str) == tissue]
        for contrast_id, rows in rows_tissue.groupby("contrast_id", sort=True):
            contrast = contrasts.get(str(contrast_id))
            if contrast is None:
                continue
            selected = rows[rows["selected"].fillna(False).astype(bool)]
            selected = selected[selected["feature_id"].astype(str).isin(values.index)]
            if selected.empty:
                continue
            groups = pairwise_groups(experiment, contrast, omics, values.columns)
            num = [s for s, g in groups.items() if g == contrast.numerator]
            den = [s for s, g in groups.items() if g == contrast.denominator]
            if not num or not den:
                continue
            ids = selected["feature_id"].astype(str).to_numpy()
            full, loo = leave_one_out(values.loc[ids, num], values.loc[ids, den])
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = loo.div(full, axis=0)
            usable = ratio.notna()
            min_ratio_value = ratio.min(axis=1)
            # 逐行取比值最小的样本；整行都无法重算（如该特征有缺失值）时为空
            filled = ratio.to_numpy(dtype=float, na_value=np.nan)
            any_usable = usable.any(axis=1).to_numpy()
            positions = np.argmin(np.where(np.isnan(filled), np.inf, filled), axis=1)
            driver = pd.Series(np.where(any_usable, ratio.columns.to_numpy()[positions], ""), index=ratio.index)
            se_all, _ = effect_standard_errors(rows)
            plausible = pd.Series(normal_shrinkage(rows["log2fc"], se_all), index=rows.index)
            se = se_all.loc[selected.index].to_numpy()
            d = plausible.loc[selected.index].to_numpy()
            power, type_s, exaggeration = retrodesign(d, se, alpha)
            flagged = (min_ratio_value < min_ratio).fillna(False)
            out.setdefault(omics, []).append(pd.DataFrame({
                "contrast_id": str(contrast_id), "tissue": tissue, "omics": omics, "feature_id": ids,
                "n_numerator": len(num), "n_denominator": len(den), "n_left_out": usable.sum(axis=1).to_numpy(),
                "full_effect": full.to_numpy(), "loo_min_ratio": min_ratio_value.to_numpy(),
                "loo_sign_flips": (ratio <= 0).sum(axis=1).to_numpy(),
                "driver_sample": driver.where(flagged, "").fillna("").astype(str).to_numpy(),
                FLAG: flagged.to_numpy(dtype=bool), "se": se, "plausible_effect": d, "power": power,
                "type_s": type_s, "exaggeration": exaggeration,
            }))
    return {omics: pd.concat(parts, ignore_index=True)[COLUMNS] for omics, parts in out.items()}


def flag_rows(table: pd.DataFrame, robust: pd.DataFrame | None) -> pd.DataFrame:
    """在差异结果中由单只动物驱动的行的 ``flags`` 中加上 ``single_subject_driven``。"""
    if robust is None or robust.empty or table.empty:
        return table
    hits = robust[robust[FLAG]]
    keys = set(zip(hits["contrast_id"], hits["tissue"], hits["feature_id"]))
    mask = [(c, t, f) in keys for c, t, f in zip(table["contrast_id"].astype(str), table["tissue"].astype(str),
                                                   table["feature_id"].astype(str))]
    if not any(mask):
        return table
    flags = table["flags"].fillna("").astype(str) if "flags" in table.columns else pd.Series("", index=table.index)
    return table.assign(flags=[";".join(p for p in (f, FLAG) if p) if m and FLAG not in f.split(";") else f
                               for f, m in zip(flags, mask)])
