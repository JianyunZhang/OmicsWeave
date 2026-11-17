#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
反应档案的逐只动物数据（G13 第二部分）：成员表引用到的特征在各数据集中的逐只动物值，供档案页（图 B）画点。

- 值的尺度与差异分析一致：计数为 log2(CPM + 1)，强度与峰面积为 log2，通量为线性（净通量可为负，不取对数）；
- ``group`` 为该数据集中有变化的动物层因素水平的组合（如 HFD，或 HFD·Sema），``group_order`` 按因素声明的水平排序，
  图的横轴取实验的全部组别，而不只是当前比较的两组；
- ``tracer`` 为 MFA 样本的示踪剂（样本表的 ``protocol``）；``qc_flag`` 取样本表；``at_bound``（逐只动物的通量是否贴在
  约束边界上）目前的 MFA 导出没有提供，为空。

另有比较成员表（:func:`contrast_samples`）：每个比较的分子组与分母组各含哪些样本，图上据此着色。
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from omicsweave.contract.dataset import read_dataset
from omicsweave.core.config import Experiment
from omicsweave.core.model import Omics
from omicsweave.core.run import contrast_groups
from omicsweave.stats.data import analysis_values, sample_table

from .sources import DossierSource

SAMPLE_COLUMNS = ["source", "layer", "tissue", "feature_id", "tracer", "sample_id", "subject_id", "group",
                  "group_order", "value", "value_type", "scale", "units", "at_bound", "qc_flag"]
CONTRAST_COLUMNS = ["source", "layer", "contrast_id", "sample_id", "side", "level"]


def group_labels(experiment: Experiment, samples: list[str]) -> tuple[dict[str, str], dict[str, int]]:
    """样本 → 组标签（在这些样本中有变化的动物层因素的水平，以 · 连接）与组的排序号。"""
    factors = {name: f for name, f in experiment.config.design.factors.items() if f.level == "subject"}
    table = sample_table(experiment)
    values = {s: {f: str(table.at[s, f]) if s in table.index and f in table.columns else "" for f in factors}
              for s in samples}
    varying = [f for f in factors if len({v[f] for v in values.values()}) > 1]
    labels = {s: "·".join(v[f] for f in varying) or "all" for s, v in values.items()}

    def rank(sample: str) -> tuple[int, ...]:
        return tuple(list(factors[f].levels or []).index(values[sample][f])
                     if values[sample][f] in (factors[f].levels or []) else 999 for f in varying)

    ordered = sorted(dict.fromkeys(labels.values()),
                     key=lambda label: min(rank(s) for s, lab in labels.items() if lab == label))
    return labels, {label: i for i, label in enumerate(ordered)}


def _dataset_dir(source: DossierSource, layer: str, tissue: str) -> Path:
    omics, kind = (Omics.METABOLOMICS.value, "derived") if layer == "derived" else (layer, "data")
    return source.run_dir / "omics" / omics / kind / source.raw_tissue(layer if layer != "derived" else omics, tissue)


def dossier_samples(source: DossierSource, members: pd.DataFrame) -> tuple[pd.DataFrame, list[Path]]:
    """``members``：该来源的成员行（dossier_members 的子集）。返回逐只动物表与读取的数据集目录。"""
    measured = members[members["feature_id"].notna() & members["feature_id"].astype(str).ne("")]
    table = sample_table(source.experiment)
    frames, inputs = [], []
    for (layer, tissue), rows in measured.groupby(["layer", "tissue"], sort=True):
        directory = _dataset_dir(source, str(layer), str(tissue))
        if not directory.is_dir():
            continue
        dataset = read_dataset(directory)
        features = [f for f in dict.fromkeys(rows["feature_id"].astype(str)) if f in dataset.matrix.index]
        if not features:
            continue
        flux = dataset.meta.omics is Omics.MFA
        values = dataset.matrix.astype(float) if flux else analysis_values(dataset)
        values = values.loc[features]
        labels, order = group_labels(source.experiment, list(values.columns))
        long = values.rename_axis(index="feature_id", columns=None).reset_index() \
            .melt(id_vars="feature_id", var_name="sample_id", value_name="value")
        meta = table.reindex(long["sample_id"])
        tracer = meta["protocol"].fillna("").astype(str).to_numpy() if flux and "protocol" in meta.columns else ""
        frames.append(pd.DataFrame({
            "source": source.label, "layer": layer, "tissue": tissue, "feature_id": long["feature_id"].astype(str),
            "tracer": tracer, "sample_id": long["sample_id"].astype(str),
            "subject_id": meta["subject_id"].astype(str).to_numpy() if "subject_id" in meta.columns else "",
            "group": long["sample_id"].map(labels), "group_order": long["sample_id"].map(labels).map(order),
            "value": long["value"].to_numpy(dtype=float), "value_type": dataset.meta.value_type.value,
            "scale": "linear" if flux else "log2", "units": dataset.meta.units or "", "at_bound": np.nan,
            "qc_flag": meta["qc_flag"].fillna("").astype(str).to_numpy() if "qc_flag" in meta.columns else "",
        }))
        inputs.append(directory)
    out = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=SAMPLE_COLUMNS)
    return out[SAMPLE_COLUMNS], inputs


def contrast_samples(source: DossierSource) -> pd.DataFrame:
    """每个（层, 比较）的分子组与分母组样本；比较 ID 为本实验的名称。"""
    contrasts = {c.id: c for c in source.experiment.config.design.contrasts}
    rows: list[dict[str, Any]] = []
    for raw, common in source.contrast_map.items():
        contrast = contrasts.get(raw)
        if contrast is None:
            continue
        for layer in source.layers:
            groups = contrast_groups(source.experiment, contrast, Omics(layer))
            for sample, level in groups.items():
                side = "numerator" if level == contrast.numerator else "denominator"
                rows.append({"source": source.label, "layer": layer, "contrast_id": common, "sample_id": sample,
                             "side": side, "level": level})
                if layer == Omics.METABOLOMICS.value:   # 派生特征与测得代谢物是同一批样本
                    rows.append({**rows[-1], "layer": "derived"})
    return pd.DataFrame(rows, columns=CONTRAST_COLUMNS)
