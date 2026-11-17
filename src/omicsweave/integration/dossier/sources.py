#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
反应档案的数据来源：本实验的源结果包，以及参数 ``sources`` 中声明的其他项目 / 实验（只读）。

每个来源读入后，组织名与比较 ID 都换成本实验的名称（:meth:`DossierSource.mapped_tables`），档案的其余部分只用
共同名称；读数据集（逐只动物值）时再换回来源自己的组织名（:meth:`DossierSource.raw_tissue`）。

- 本实验：``tissue_map`` 只用于 MFA 的组织名（与 ``flux_anchored`` 相同）；
- 其他来源：``tissue_map`` 用于该来源全部组学的组织名；``contrast_map``（来源的比较 → 本实验的比较）没写时，
  只保留与本实验同名的比较；``layers`` 限定读取的组学。来源的运行按 ``run``（运行 ID 或 latest，取最近一次
  含统一差异结果的运行）解析，项目路径相对本项目目录。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd

from omicsweave.contract.dataset import ContractError
from omicsweave.contract.query import open_run
from omicsweave.contract.tables import read_differential
from omicsweave.core.config import ConfigError, Experiment, Project, load_experiment, load_project
from omicsweave.core.model import Omics
from omicsweave.core.run import RunError

DERIVED_TABLE = "derived_differential.parquet"
SOURCE_KEYS = {"label", "project", "experiment", "run", "layers", "tissue_map", "contrast_map"}
LAYERS = tuple(o.value for o in (Omics.TRANSCRIPTOMICS, Omics.PROTEOMICS, Omics.METABOLOMICS, Omics.MFA))


@dataclass
class DossierSource:
    label: str
    experiment: Experiment
    run_dir: Path
    tissue_map: dict[str, str] = field(default_factory=dict)     # 来源的组织名 → 共同组织名
    contrast_map: dict[str, str] = field(default_factory=dict)   # 来源的比较 → 本实验的比较
    layers: tuple[str, ...] = LAYERS
    is_primary: bool = False                                      # 本实验：tissue_map 只用于 MFA；声明的示踪剂只用于它
    tables: dict[str, pd.DataFrame] = field(default_factory=dict)
    derived: pd.DataFrame | None = None
    inputs: list[Path] = field(default_factory=list)

    @property
    def run_id(self) -> str:
        return self.run_dir.name

    def common_tissue(self, layer: str, tissue: str) -> str:
        if self.is_primary and layer != Omics.MFA.value:
            return tissue
        return self.tissue_map.get(tissue, tissue)

    def raw_tissue(self, layer: str, tissue: str) -> str:
        """共同组织名 → 来源自己的组织名（读数据集用）。"""
        if self.is_primary and layer != Omics.MFA.value:
            return tissue
        return next((raw for raw, common in self.tissue_map.items() if common == tissue), tissue)

    def raw_contrast(self, contrast: str) -> str:
        return next((raw for raw, common in self.contrast_map.items() if common == contrast), contrast)

    def mapped(self, table: pd.DataFrame | None, layer: str) -> pd.DataFrame | None:
        """换成共同的组织名与比较 ID；其他来源只保留映射到本实验的比较。"""
        if table is None or table.empty:
            return table
        out = table[table["contrast_id"].astype(str).isin(self.contrast_map)].copy()
        out["contrast_id"] = out["contrast_id"].astype(str).map(self.contrast_map)
        out["tissue"] = out["tissue"].astype(str).map(lambda t: self.common_tissue(layer, t))
        return out

    def mapped_tables(self) -> dict[str, pd.DataFrame]:
        return {omics: self.mapped(table, omics) for omics, table in self.tables.items()}


def _read_tables(run_dir: Path, layers: tuple[str, ...]) -> tuple[dict[str, pd.DataFrame], pd.DataFrame | None, list[Path]]:
    from omicsweave.integration.modes.flux_anchored import refresh_metabolite_ids

    paths = {p.parent.parent.name: p for p in sorted(run_dir.glob("omics/*/tables/differential.parquet"))
             if p.parent.parent.name in layers}
    tables = {omics: read_differential(path) for omics, path in paths.items()}
    if Omics.METABOLOMICS.value in tables:
        tables[Omics.METABOLOMICS.value], _ = refresh_metabolite_ids(tables[Omics.METABOLOMICS.value])
    derived_path = run_dir / "omics" / Omics.METABOLOMICS.value / "tables" / DERIVED_TABLE
    derived = read_differential(derived_path) if derived_path.is_file() and Omics.METABOLOMICS.value in layers else None
    inputs = list(paths.values()) + ([derived_path] if derived is not None else [])
    return tables, derived, inputs


def primary_source(experiment: Experiment, run_dir: Path, label: str, tissue_map: dict[str, str]) -> DossierSource:
    tables, derived, inputs = _read_tables(run_dir, LAYERS)
    if not tables:
        raise RunError(f"结果包 {run_dir.name} 中没有统一差异结果；请先运行 omicsweave run（含 differential 步骤）")
    contrasts = {str(c) for t in tables.values() for c in t["contrast_id"].astype(str).unique()}
    return DossierSource(label, experiment, run_dir, tissue_map, {c: c for c in contrasts}, LAYERS, True, tables,
                         derived, inputs)


def other_source(project: Project, spec: dict[str, Any], primary_contrasts: set[str]) -> DossierSource:
    unknown = sorted(set(spec) - SOURCE_KEYS)
    if unknown or not spec.get("label") or not spec.get("experiment"):
        raise RunError(f"reaction_dossier.sources 的条目需要 label 与 experiment，可选 {sorted(SOURCE_KEYS)}；"
                       f"无法识别：{unknown}")
    label = str(spec["label"])
    root = project.resolve(spec["project"]) if spec.get("project") else project.root
    try:
        other = load_project(root)
        experiment = load_experiment(other, str(spec["experiment"]))
        bundle = open_run(str(spec["experiment"]), str(spec.get("run", "latest")), project_dir=root,
                          requires="omics/*/tables/differential.parquet")
    except (ConfigError, ContractError, FileNotFoundError) as exc:
        raise RunError(f"来源 {label}：{exc}") from exc
    layers = tuple(str(x) for x in (spec.get("layers") or LAYERS))
    bad = sorted(set(layers) - set(LAYERS))
    if bad:
        raise RunError(f"来源 {label} 的 layers 中有无法识别的组学：{bad}（可选 {list(LAYERS)}）")
    tables, derived, inputs = _read_tables(bundle.path, layers)
    contrasts = {str(c) for t in tables.values() for c in t["contrast_id"].astype(str).unique()}
    contrast_map = {str(k): str(v) for k, v in (spec.get("contrast_map") or {}).items()} or \
        {c: c for c in contrasts & primary_contrasts}
    missing = sorted(set(contrast_map) - contrasts)
    if missing:
        raise RunError(f"来源 {label}（运行 {bundle.path.name}）中没有比较 {missing}")
    if not contrast_map:
        raise RunError(f"来源 {label} 没有与本实验对应的比较：用 contrast_map 声明（来源的比较 → 本实验的比较）")
    tissue_map = {str(k): str(v) for k, v in (spec.get("tissue_map") or {}).items()}
    return DossierSource(label, experiment, bundle.path, tissue_map, contrast_map, layers, False, tables, derived,
                         inputs)
