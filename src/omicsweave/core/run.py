#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
统一运行入口（``omicsweave run``）：按 ``experiment.yaml`` 运行一个或多个组学，生成一个结果包。

结果包目录 ``<paths.results>/<experiment_id>/<run_id>/``（v2 设计文档第 8.1 节）::

    manifest.json     运行 ID、命令、配置哈希、输入输出文件的校验值
    metadata/         subjects / samples（Parquet）与展开后的实验配置
    omics/<组学>/      data/（分析层数据集）、qc/、tables/（统一结果表）、figures/（图件）、details/（旧版格式的明细表）
    qc/               contrast_reliability.parquet：样本 QC 标记对各（组学, 组织, 比较）的影响（见 core.reliability）

每个组学由一个运行器实现（:data:`RUNNERS`），各组学之间互不依赖。运行器只通过 :class:`RunContext`
取得项目、实验和输出目录，不读写进程级全局状态。

``figures`` 步骤由本次运行写出的统一结果表绘图（见 :mod:`omicsweave.visualization`），所以只指定部分步骤时
须同时包含写出差异结果表的 ``differential`` 或 ``enrichment``。
"""
from __future__ import annotations

import importlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

from omicsweave.contract.manifest import MANIFEST_FILE, RunManifest, allocate_run_dir, file_record, write_manifest

from .config import EXPERIMENT_FILE, ConfigError, Experiment, Project, load_experiment, load_project
from .model import Contrast, Omics
from .provenance import sha256_file
from .validate import biological_samples, matches_level, validate_experiment

# 已接入统一运行的组学：组学 → "模块:函数"。模块另以 STEPS 声明可选步骤
RUNNERS: dict[Omics, str] = {
    Omics.MFA: "omicsweave.mfa.runner:run",
    Omics.METABOLOMICS: "omicsweave.metabolomics.runner:run",
    Omics.TRANSCRIPTOMICS: "omicsweave.transcriptomics.runner:run",
    Omics.PROTEOMICS: "omicsweave.proteomics.runner:run",
}
ALL = "all"
METADATA_DIR = "metadata"
FIGURES = "figures"
# 写出 tables/differential.parquet 的步骤；figures 读取它
_TABLE_STEPS = {"differential", "enrichment"}
DELIVERY_META_FILES = frozenset({"index.json", ".index.tmp", ".keep", "summary.xlsx"})
DELIVERY_ONLY_DIRS = frozenset({"candidates", "panels", "pathway_map", "figures", "key_results", "exports"})


def is_bundle_output(path: Path, run_dir: Path) -> bool:
    """是否属于结果包本体输出（排除清单、交付索引、汇总表及交付副本/课题脚本目录）。"""
    rel = path.relative_to(run_dir)
    if len(rel.parts) == 1 and rel.name in {MANIFEST_FILE, *DELIVERY_META_FILES}:
        return False
    if rel.parts[0] in DELIVERY_ONLY_DIRS:
        return False
    return True


class RunError(Exception):
    """运行无法开始，或某个组学运行失败。"""


@dataclass
class RunContext:
    project: Project
    experiment: Experiment
    run_id: str
    run_dir: Path
    # None 表示运行该组学的全部步骤
    steps: tuple[str, ...] | None = None
    # None 表示按各组学配置（默认 False）决定是否运行预排序 GSEA
    run_gsea: bool | None = None

    def omics_dir(self, omics: Omics) -> Path:
        return self.run_dir / "omics" / omics.value

    def wants(self, step: str) -> bool:
        return self.steps is None or step in self.steps


@dataclass
class OmicsRun:
    """一个组学运行器的产出；输出文件不必列出，结束时统一登记结果包内的全部文件。"""
    omics: Omics
    inputs: list[Path] = field(default_factory=list)
    parameters: dict[str, Any] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)


def contrast_groups(experiment: Experiment, contrast: Contrast, omics: Omics) -> dict[str, str]:
    """一个比较涉及的某组学生物学样本 → 所在组（分子或分母的水平）；先按 ``subset`` 筛选。"""
    factors = experiment.config.design.factors
    table = biological_samples(experiment)
    table = table[table["omics"] == omics.value]
    for key, value in contrast.subset.items():
        table = table[matches_level(table[key], value, factors[key].type == "numeric")]
    numeric = factors[contrast.factor].type == "numeric"
    groups: dict[str, str] = {}
    for level in (contrast.numerator, contrast.denominator):
        selected = table[matches_level(table[contrast.factor], level, numeric)]
        groups.update(dict.fromkeys(selected["sample_id"], level))
    return groups


def sample_labels(experiment: Experiment, samples: list[str]) -> dict[str, str]:
    """样本 → 动物层因素水平的组合（如 HFD/GLP1），用作实验分组：

    预处理中按组判断检出率、全局 PLS-DA 等需要“全部样本的一种分组”的地方使用它；
    只有一个因素时就是该因素的水平。没有动物层因素时所有样本为同一组 ``all``。
    """
    factors = [name for name, f in experiment.config.design.factors.items() if f.level == "subject"]
    table = biological_samples(experiment).drop_duplicates("sample_id").set_index("sample_id")
    labels = {}
    for sample in samples:
        values = [str(table.at[sample, f]) for f in factors if sample in table.index and f in table.columns]
        labels[sample] = "/".join(v for v in values if v) or "all"
    return labels


def _load_runner(omics: Omics) -> tuple[Callable[[RunContext], OmicsRun], tuple[str, ...]]:
    module_name, func_name = RUNNERS[omics].split(":")
    try:
        module = importlib.import_module(module_name)
    except ModuleNotFoundError as exc:
        if exc.name and exc.name.split(".")[0] == "omicsweave":
            raise
        # 各组学的分析依赖是可选依赖组，组名与组学名相同
        raise RunError(f"运行 {omics.value} 需要可选依赖 {exc.name}，请安装：pip install \"omicsweave[{omics.value}]\"") from exc
    return getattr(module, func_name), tuple(getattr(module, "STEPS", ()))


def select_omics(experiment: Experiment, requested: list[str]) -> tuple[list[Omics], list[str]]:
    enabled = experiment.config.omics.enabled()
    if ALL in requested:
        chosen = [o for o in enabled if o in RUNNERS]
        skipped = [o.value for o in enabled if o not in RUNNERS]
        notes = [f"以下组学尚未接入统一运行，已跳过（仍可用 omicsweave <组学> run）：{'、'.join(skipped)}"] if skipped else []
        if not chosen:
            raise RunError("实验中没有可统一运行的组学" + (f"（已启用：{'、'.join(o.value for o in enabled)}）" if enabled
                                                     else "（experiment.yaml 中没有启用任何组学）"))
        return chosen, notes
    chosen = []
    for name in requested:
        try:
            omics = Omics(name)
        except ValueError:
            raise RunError(f"未知组学 {name!r}，可选：{'、'.join(o.value for o in Omics)}、{ALL}") from None
        if omics not in enabled:
            raise RunError(f"{name} 未在 experiment.yaml 中启用")
        if omics not in RUNNERS:
            raise RunError(f"{name} 尚未接入统一运行（阶段 2 进行中），请暂用 omicsweave {name} run")
        if omics not in chosen:
            chosen.append(omics)
    return chosen, []


def _write_metadata(ctx: RunContext) -> None:
    directory = ctx.run_dir / METADATA_DIR
    directory.mkdir(parents=True, exist_ok=True)
    ctx.experiment.subjects.to_parquet(directory / "subjects.parquet", index=False)
    ctx.experiment.samples.to_parquet(directory / "samples.parquet", index=False)
    config = ctx.experiment.config.model_dump(mode="json")
    (directory / "experiment.json").write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")


def _expand(paths: list[Path]) -> list[Path]:
    """输入可以是目录（如 Parquet 目录），登记其中的全部文件；按路径去重并保持顺序。"""
    files: list[Path] = []
    for path in paths:
        files += sorted(p for p in path.rglob("*") if p.is_file()) if path.is_dir() else [path]
    return list(dict.fromkeys(p.resolve() for p in files))


def load_validated(experiment_id: str, root: str | Path | Project | None = None) -> tuple[Project, Experiment]:
    """读取项目与实验；配置未通过校验时抛出 :class:`RunError`，不开始运行。"""
    try:
        project = load_project(root)
        experiment = load_experiment(project, experiment_id)
    except ConfigError as exc:
        raise RunError("\n".join(exc.messages)) from exc
    errors = [i.message for i in validate_experiment(experiment) if i.level == "error"]
    if errors:
        raise RunError(f"实验 {experiment.id} 未通过校验（omicsweave validate），请先修正：\n"
                       + "\n".join(f"  - {e}" for e in errors))
    return project, experiment


def start_run(project: Project, experiment: Experiment, run_id: str | None = None,
              steps: list[str] | None = None, *, into_existing: bool = False,
              run_gsea: bool | None = None) -> RunContext:
    """新建结果包目录并写出 metadata/。"""
    base = project.path("results") / experiment.id
    if run_id is None:
        if into_existing:
            raise RunError("into_existing=True 时必须指定 run_id")
        run_dir = allocate_run_dir(base)
        run_id = run_dir.name
    else:
        run_dir = base / run_id
        if into_existing:
            if (run_dir / MANIFEST_FILE).exists() or (run_dir / "omics").exists() or (run_dir / "integration").exists():
                raise RunError(f"{run_dir} 已存在")
        elif run_dir.exists():
            raise RunError(f"{run_dir} 已存在")
    ctx = RunContext(project, experiment, run_id, run_dir, tuple(steps) if steps else None, run_gsea=run_gsea)
    _write_metadata(ctx)
    return ctx


def finish_run(ctx: RunContext, inputs: list[Path], *, omics: list[Omics], parameters: dict[str, Any],
               command: list[str] | None = None, input_kind: str | None = "omics") -> Path:
    """写出运行清单：登记实验配置、样本表与 ``inputs``，以及结果包内的全部输出文件。"""
    experiment = ctx.experiment
    design = experiment.config.design
    files = _expand([experiment.root / EXPERIMENT_FILE, experiment.root / design.subjects,
                     experiment.root / design.samples, *inputs])
    outputs = sorted(p for p in ctx.run_dir.rglob("*") if p.is_file() and is_bundle_output(p, ctx.run_dir))
    manifest = RunManifest(
        run_id=ctx.run_id,
        project_id=ctx.project.id,
        experiment_id=experiment.id,
        command=command or [],
        config_sha256=sha256_file(experiment.root / EXPERIMENT_FILE),
        omics=omics,
        parameters=parameters,
        files=[file_record(p, ctx.project.root, "input", kind=input_kind) for p in files]
        + [file_record(p, ctx.run_dir, "output") for p in outputs],
    )
    return write_manifest(manifest, ctx.run_dir)


def run_experiment(experiment_id: str, omics: list[str] | tuple[str, ...] = (ALL,),
                   steps: list[str] | None = None, root: str | Path | Project | None = None,
                   command: list[str] | None = None, run_id: str | None = None,
                   *, into_existing: bool = False, run_gsea: bool | None = None) -> tuple[Path, list[str]]:
    """运行一个实验，返回 (结果包目录, 提示)。配置未通过校验时不开始运行。"""
    project, experiment = load_validated(experiment_id, root)
    chosen, notes = select_omics(experiment, list(omics))
    runners = {o: _load_runner(o) for o in chosen}
    if steps:
        known = sorted(set().union(*(s for _, s in runners.values())))
        unknown = [s for s in steps if s not in known]
        if unknown:
            raise RunError(f"未知步骤 {'、'.join(unknown)}；{'、'.join(o.value for o in chosen)} 可选：{'、'.join(known)}")
        if FIGURES in steps and not _TABLE_STEPS & set(steps):
            raise RunError(f"{FIGURES} 步骤根据本次运行写出的差异结果表绘图，请与 differential 一起指定，"
                           f"如 --steps differential,{FIGURES}")

    ctx = start_run(project, experiment, run_id, steps, into_existing=into_existing, run_gsea=run_gsea)
    runs: list[OmicsRun] = []
    for o in chosen:
        func, _ = runners[o]
        try:
            runs.append(func(ctx))
        except ModuleNotFoundError as exc:
            if exc.name and exc.name.split(".")[0] == "omicsweave":
                raise
            pkg = exc.name.split(".")[0] if exc.name else str(exc)
            raise RunError(
                f"运行 {o.value} 需要可选依赖 {pkg}，请安装：pip install \"omicsweave[{o.value}]\""
                f"（结果包 {ctx.run_dir} 不完整，没有 manifest.json）"
            ) from exc
        except RunError as exc:
            raise RunError(f"{o.value} 运行失败：{exc}（结果包 {ctx.run_dir} 不完整，没有 manifest.json）") from exc

    # 后处理：派生特征、公式模型、可检出效应、稳健性与样本 QC 可靠性（omicsweave.stats.post）
    if not steps or _TABLE_STEPS & set(steps):
        from omicsweave.stats.post import postprocess
        notes += postprocess(experiment, ctx.run_dir)

    parameters: dict[str, Any] = {"steps": list(steps) if steps else ALL, **{r.omics.value: r.parameters for r in runs}}
    if run_gsea is not None:
        parameters["run_gsea"] = run_gsea
    finish_run(ctx, [p for r in runs for p in r.inputs], omics=chosen, command=command, parameters=parameters)
    return ctx.run_dir, notes + [n for r in runs for n in r.notes]
