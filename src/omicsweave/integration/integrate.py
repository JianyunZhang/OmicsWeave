#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
多组学整合入口（``omicsweave integrate``）：读取一次单组学运行（``omicsweave run``）写出的统一结果表，
按整合模式分析，结果写入新的结果包 ``<paths.results>/<实验>/<运行ID>/integration/<模式>/``。

整合只读取数据契约（源结果包中的统一结果表与数据集），不导入各组学模块。运行清单登记源结果包的运行 ID，
以及读取的每个文件的校验值。未指定 ``--mode`` 时运行 ``experiment.yaml`` 中 ``integration.modes`` 列出的模式；
各模式的参数写在 ``integration.parameters.<模式>`` 中。
"""
from __future__ import annotations

import importlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

from omicsweave.contract.manifest import MANIFEST_FILE, RunManifest, file_record, read_manifest, write_manifest
from omicsweave.core.config import EXPERIMENT_FILE, Project
from omicsweave.core.model import IntegrationMode
from omicsweave.core.run import (
    RunContext,
    RunError,
    _expand,
    finish_run,
    is_bundle_output,
    load_validated,
    start_run,
)

# 已实现的整合模式：模式 → "模块:函数"
MODES: dict[IntegrationMode, str] = {
    IntegrationMode.FLUX_ANCHORED: "omicsweave.integration.modes.flux_anchored:run",
    IntegrationMode.PATHWAY_CONVERGENCE: "omicsweave.integration.modes.pathway_convergence:run",
    IntegrationMode.REVERSAL: "omicsweave.integration.modes.reversal:run",
    IntegrationMode.REACTION_DOSSIER: "omicsweave.integration.modes.reaction_dossier:run",
}


@dataclass
class IntegrationRun:
    """一个整合模式的产出；输出文件结束时统一登记。"""
    mode: IntegrationMode
    inputs: list[Path] = field(default_factory=list)
    parameters: dict[str, Any] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)


def mode_dir(ctx: RunContext, mode: IntegrationMode) -> Path:
    return ctx.run_dir / "integration" / mode.value


def source_run(results_root: str | Path, experiment_id: str, run_id: str | None = None) -> Path:
    """源结果包：指定的运行，或最近一次含统一结果表的运行。"""
    base = Path(results_root) / experiment_id
    if run_id:
        candidates = [base / run_id]
    else:
        candidates = sorted(p for p in base.glob("*") if (p / MANIFEST_FILE).is_file()) if base.is_dir() else []
    runs = [p for p in candidates if (p / MANIFEST_FILE).is_file() and any(p.glob("omics/*/tables/*.parquet"))]
    if not runs:
        where = f"{base / run_id}" if run_id else f"{base} 下"
        raise RunError(f"{where}没有含统一结果表的结果包，请先运行 omicsweave run --experiment {experiment_id}")
    return runs[-1]


def _load_mode(mode: IntegrationMode) -> Callable[..., IntegrationRun]:
    module_name, func_name = MODES[mode].split(":")
    return getattr(importlib.import_module(module_name), func_name)


def select_modes(requested: list[str], configured: list[IntegrationMode]) -> list[IntegrationMode]:
    if not requested:
        if not configured:
            raise RunError("没有指定整合模式：用 --mode 指定，或在 experiment.yaml 的 integration.modes 中列出")
        requested = [m.value for m in configured]
    modes = []
    for name in requested:
        try:
            mode = IntegrationMode(name)
        except ValueError:
            raise RunError(f"未知整合模式 {name!r}，可选：{'、'.join(m.value for m in IntegrationMode)}") from None
        if mode not in MODES:
            raise RunError(f"整合模式 {mode.value} 尚未接入，当前可用："
                           f"{'、'.join(m.value for m in MODES)}；也可暂用 omicsweave integration run")
        if mode not in modes:
            modes.append(mode)
    return modes


def integrate_experiment(experiment_id: str, modes: list[str] | None = None, *, source: str | None = None,
                         root: str | Path | Project | None = None, command: list[str] | None = None,
                         run_id: str | None = None, into_source: bool = False,
                         into_existing: bool = False) -> tuple[Path, list[str]]:
    """运行整合模式，返回 (结果包目录, 提示)。``into_source=True`` 时直接写入源结果包并合并清单。"""
    project, experiment = load_validated(experiment_id, root)
    settings = experiment.config.integration
    chosen = select_modes(list(modes or []), settings.modes)
    src = source_run(project.path("results"), experiment.id, source)
    if into_source:
        ctx = RunContext(project, experiment, src.name, src)
    else:
        ctx = start_run(project, experiment, run_id, into_existing=into_existing)
    runs = []
    for mode in chosen:
        try:
            runs.append(_load_mode(mode)(ctx, src, dict(settings.parameters.get(mode, {}))))
        except RunError as exc:
            raise RunError(f"{mode.value} 失败：{exc}（结果包 {ctx.run_dir} 不完整，没有 manifest.json）") from exc
    if into_source:
        existing = read_manifest(src)
        design = experiment.config.design
        candidates = _expand([
            experiment.root / EXPERIMENT_FILE,
            experiment.root / design.subjects,
            experiment.root / design.samples,
            *(p for r in runs for p in r.inputs),
        ])
        input_records = [f for f in existing.files if f.role == "input"]
        seen = {f.path for f in input_records}
        src_resolved = src.resolve()
        for path in candidates:
            if path.resolve().is_relative_to(src_resolved):
                continue
            record = file_record(path, project.root, "input", kind="integration")
            if record.path not in seen:
                seen.add(record.path)
                input_records.append(record)
        outputs = sorted(p for p in src.rglob("*") if p.is_file() and is_bundle_output(p, src))
        merged_params = {
            **existing.parameters,
            "source_run": src.name,
            "modes": [m.value for m in chosen],
            **{r.mode.value: r.parameters for r in runs},
        }
        manifest = RunManifest(
            run_id=existing.run_id,
            project_id=existing.project_id,
            experiment_id=existing.experiment_id,
            created_at=existing.created_at,
            command=existing.command or command or [],
            config_sha256=existing.config_sha256,
            omics=existing.omics,
            parameters=merged_params,
            files=input_records + [file_record(p, src, "output") for p in outputs],
        )
        write_manifest(manifest, src)
    else:
        finish_run(ctx, [p for r in runs for p in r.inputs], omics=[], command=command,
                   parameters={"source_run": src.name, "modes": [m.value for m in chosen],
                               **{r.mode.value: r.parameters for r in runs}},
                   input_kind="integration")
    return ctx.run_dir, [n for r in runs for n in r.notes]
