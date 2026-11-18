#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
结果包查询接口：按实验与运行读取统一结果表、派生特征、数据集与整合结果，课题脚本不必写死路径或运行 ID。

运行可以指定为 ``"latest"``（默认：该实验最近一次含所需表的完整结果包）或显式的运行 ID。读出的表在
``DataFrame.attrs`` 中带 ``run_id`` 与 ``run_dir``；:class:`Bundle` 本身也有 ``run_id``。

示例::

    from omicsweave.contract.query import open_run, read_table

    diff = read_table("all_tissues", "differential", omics="mfa", project_dir="local/projects/dio_project")
    diff.attrs["run_id"]                                     # 实际读取的运行
    bundle = open_run("all_tissues", "run_20261112_163447_930f", project_dir="local/projects/dio_project")
    profile = bundle.integration("flux_anchored", "flux_evidence_profile")
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from omicsweave.core.config import load_project
from omicsweave.core.project import project_dir as default_project_dir

from .dataset import ContractError, Dataset, read_dataset
from .manifest import MANIFEST_FILE, RunManifest, read_manifest

LATEST = "latest"
# omics/<组学>/tables/ 下的表
OMICS_TABLES = ("differential", "derived_differential", "enrichment")


def _tag(table: pd.DataFrame, bundle: "Bundle") -> pd.DataFrame:
    table.attrs.update(run_id=bundle.run_id, run_dir=str(bundle.path))
    return table


@dataclass(frozen=True)
class Bundle:
    path: Path
    manifest: RunManifest

    @property
    def run_id(self) -> str:
        return self.manifest.run_id

    @property
    def experiment_id(self) -> str:
        return self.manifest.experiment_id

    def tables(self) -> list[str]:
        """结果包中的全部表（Parquet），路径相对于结果包目录。"""
        return sorted(p.relative_to(self.path).as_posix() for p in self.path.rglob("*.parquet"))

    def read(self, relative: str) -> pd.DataFrame:
        path = self.path / relative
        if not path.is_file():
            raise ContractError(f"结果包 {self.run_id} 中没有 {relative}")
        return _tag(pd.read_parquet(path), self)

    def omics(self) -> list[str]:
        return sorted(p.name for p in (self.path / "omics").iterdir() if p.is_dir()) \
            if (self.path / "omics").is_dir() else []

    def omics_table(self, table: str, omics: str | None = None) -> pd.DataFrame:
        """``omics/<组学>/tables/<表>.parquet``；不指定组学时合并全部组学（统一结果表都有 ``omics`` 列）。"""
        names = [omics] if omics else self.omics()
        paths = [self.path / "omics" / name / "tables" / f"{table}.parquet" for name in names]
        found = [p for p in paths if p.is_file()]
        if not found:
            where = f"组学 {omics} " if omics else ""
            raise ContractError(f"结果包 {self.run_id} 中{where}没有 {table}.parquet")
        frames = [pd.read_parquet(p) for p in found]
        return _tag(frames[0] if len(frames) == 1 else pd.concat(frames, ignore_index=True), self)

    def differential(self, omics: str | None = None) -> pd.DataFrame:
        return self.omics_table("differential", omics)

    def derived(self, omics: str = "metabolomics") -> pd.DataFrame:
        return self.omics_table("derived_differential", omics)

    def enrichment(self, omics: str | None = None) -> pd.DataFrame:
        return self.omics_table("enrichment", omics)

    def integration(self, mode: str, table: str) -> pd.DataFrame:
        return self.read(f"integration/{mode}/tables/{table}.parquet")

    def reliability(self, mode: str | None = None) -> pd.DataFrame:
        """组学 × 组织 × 比较的 QC 可靠性：组学结果包的 ``qc/``，或整合模式的 ``integration/<模式>/qc/``。"""
        base = f"integration/{mode}/qc" if mode else "qc"
        return self.read(f"{base}/contrast_reliability.parquet")

    def dataset(self, omics: str, tissue: str, kind: str = "data") -> Dataset:
        """分析层数据集 ``omics/<组学>/<kind>/<组织>/``（``kind`` 为 data，或代谢组派生特征的 derived）。"""
        return read_dataset(self.path / "omics" / omics / kind / tissue)


def list_runs(results_root: str | Path, experiment_id: str) -> list[Path]:
    """实验的完整结果包（有 ``manifest.json``），按运行 ID（即时间）排序。"""
    base = Path(results_root) / experiment_id
    return sorted(p for p in base.glob("*") if (p / MANIFEST_FILE).is_file()) if base.is_dir() else []


def open_run(experiment_id: str, run: str = LATEST, *, project_dir: str | Path | None = None,
             requires: str | None = None) -> Bundle:
    """打开结果包。``run`` 为运行 ID 或 ``"latest"``；``requires`` 是相对结果包的 glob（如
    ``omics/*/tables/differential.parquet``），取最近一次含匹配文件的运行。"""
    project = load_project(project_dir if project_dir is not None else default_project_dir())
    results = project.path("results")
    if run != LATEST:
        path = results / experiment_id / run
        if not (path / MANIFEST_FILE).is_file():
            raise ContractError(f"{path} 不是完整的结果包（没有 {MANIFEST_FILE}）")
        return Bundle(path, read_manifest(path))
    runs = [p for p in list_runs(results, experiment_id) if requires is None or any(p.glob(requires))]
    if not runs:
        need = f"含 {requires} 的" if requires else ""
        raise ContractError(f"{results / experiment_id} 下没有{need}完整结果包")
    return Bundle(runs[-1], read_manifest(runs[-1]))


def read_table(experiment_id: str, table: str, *, omics: str | None = None, mode: str | None = None,
               run: str = LATEST, project_dir: str | Path | None = None) -> pd.DataFrame:
    """读取一张结果表。``mode`` 给出时读整合结果 ``integration/<mode>/tables/<table>.parquet``，否则读组学结果
    （``table`` 为 differential / derived_differential / enrichment，``omics`` 省略时合并全部组学）。
    ``run="latest"`` 时取最近一次含这张表的运行；返回的表在 ``attrs["run_id"]`` 中记录实际读取的运行。"""
    if mode:
        pattern = f"integration/{mode}/tables/{table}.parquet"
        return open_run(experiment_id, run, project_dir=project_dir, requires=pattern).read(pattern)
    if table not in OMICS_TABLES:
        raise ContractError(f"组学结果表应为 {' / '.join(OMICS_TABLES)}；整合结果请指定 mode")
    pattern = f"omics/{omics or '*'}/tables/{table}.parquet"
    return open_run(experiment_id, run, project_dir=project_dir, requires=pattern).omics_table(table, omics)
