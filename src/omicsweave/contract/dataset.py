#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
分析层数据集：所有统计与整合的唯一输入（v2 设计文档第 6.2 节）。

磁盘结构::

    <dataset_dir>/
    ├── matrix.parquet      # 行 = 特征（feature_id），列 = sample_id
    ├── features.parquet    # 特征注释，索引为 feature_id，与 matrix 行一一对应
    └── dataset.json        # 数值语义与来源（DatasetMeta）
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict, Field, model_validator

from omicsweave.core.model import Omics

MATRIX_FILE = "matrix.parquet"
FEATURES_FILE = "features.parquet"
META_FILE = "dataset.json"


class ValueType(str, Enum):
    COUNTS = "counts"
    TPM = "tpm"
    INTENSITY = "intensity"
    PEAK_AREA = "peak_area"
    CONCENTRATION = "concentration"
    MID_FRACTION = "mid_fraction"
    FLUX = "flux"
    # 由测得特征按样本计算的派生量（如代谢组的脂质类别总量、质量作用比）；定义写在 features 中
    DERIVED = "derived"


class ContractError(ValueError):
    """数据不符合契约。"""


class DatasetMeta(BaseModel):
    model_config = ConfigDict(extra="forbid")

    schema_version: Literal[1] = 1
    omics: Omics
    tissue: str
    value_type: ValueType
    scale: Literal["linear", "log2"]
    normalization: str = "none"
    feature_id_type: str
    units: str | None = None
    # 适配器名称与版本、外部工具版本、输入文件校验值等
    source: dict[str, Any] = Field(default_factory=dict)

    @model_validator(mode="after")
    def _check_semantics(self) -> DatasetMeta:
        if self.value_type in (ValueType.FLUX, ValueType.CONCENTRATION) and not self.units:
            raise ValueError(f"value_type 为 {self.value_type.value} 时必须给出 units")
        if self.value_type is ValueType.COUNTS and self.scale != "linear":
            raise ValueError("counts 只能是 linear 尺度")
        return self


@dataclass
class Dataset:
    matrix: pd.DataFrame
    features: pd.DataFrame
    meta: DatasetMeta

    def check(self) -> None:
        """检查矩阵与特征表的一致性；不一致时抛出 :class:`ContractError`。"""
        problems = []
        if self.matrix.index.has_duplicates:
            problems.append("matrix 的 feature_id 有重复")
        if self.matrix.columns.has_duplicates:
            problems.append("matrix 的 sample_id 有重复")
        if self.features.index.has_duplicates:
            problems.append("features 的 feature_id 有重复")
        if not self.matrix.index.equals(self.features.index):
            problems.append("matrix 与 features 的 feature_id 不一致（集合或顺序不同）")
        non_numeric = [c for c in self.matrix.columns if not pd.api.types.is_numeric_dtype(self.matrix[c])]
        if non_numeric:
            problems.append(f"matrix 中有非数值列：{', '.join(map(str, non_numeric[:5]))}")
        elif self.meta.value_type is ValueType.COUNTS and (self.matrix.to_numpy(dtype=float) < 0).any():
            problems.append("counts 中有负值")
        if problems:
            raise ContractError("；".join(problems))


def write_dataset(directory: str | Path, dataset: Dataset) -> Path:
    dataset.check()
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    matrix = dataset.matrix.copy()
    matrix.index.name = "feature_id"
    matrix.columns = matrix.columns.astype(str)
    features = dataset.features.copy()
    features.index.name = "feature_id"
    matrix.to_parquet(directory / MATRIX_FILE)
    features.to_parquet(directory / FEATURES_FILE)
    (directory / META_FILE).write_text(
        json.dumps(dataset.meta.model_dump(mode="json"), ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return directory


def read_dataset(directory: str | Path) -> Dataset:
    directory = Path(directory)
    for name in (MATRIX_FILE, FEATURES_FILE, META_FILE):
        if not (directory / name).is_file():
            raise ContractError(f"{directory} 中缺少 {name}")
    meta = DatasetMeta.model_validate_json((directory / META_FILE).read_text(encoding="utf-8"))
    dataset = Dataset(
        matrix=pd.read_parquet(directory / MATRIX_FILE),
        features=pd.read_parquet(directory / FEATURES_FILE),
        meta=meta,
    )
    dataset.check()
    return dataset


def empty_features(feature_ids: pd.Index) -> pd.DataFrame:
    """没有注释时使用的空特征表。"""
    return pd.DataFrame(index=pd.Index(feature_ids, name="feature_id"), columns=[], dtype=np.float64)
