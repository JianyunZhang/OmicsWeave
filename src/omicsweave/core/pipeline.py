#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
标准流程入口（``omicsweave pipeline``）：一条命令完成项目的标准分析。

对 ``project.yaml`` 中 ``pipeline.experiments`` 列出的每个实验（未写时为项目中的全部实验）依次运行::

    validate   校验项目与这些实验（有错误时一步也不运行）
    run        单组学分析，生成结果包（同 omicsweave run）
    integrate  读取本次 run 的结果包做多组学整合（同 omicsweave integrate --source <本次 run 的运行 ID>）；
               实验没有启用整合（integration.enabled 与 integration.modes）时跳过
    scripts    pipeline.scripts 中的课题脚本（差异图表、关键结果导出等），占位符换成本次的运行 ID

与逐条运行这些命令相比，各步之间传递的是确切的运行 ID，而不是“最近一次运行”。实验名、脚本文件与占位符在
开始运行前全部检查；某一步失败时停止，并给出续跑命令（``--from integrate|scripts --source <运行ID>``）。

脚本在项目目录下运行，环境变量 ``OMICSWEAVE_PROJECT_DIR`` 指向项目目录；``.py`` 脚本用运行 omicsweave 的
同一个 Python 解释器执行。配置格式见 :class:`omicsweave.core.model.PipelineScript`。
"""
from __future__ import annotations

import os
import shlex
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

from omicsweave.contract.manifest import MANIFEST_FILE, read_manifest

from .config import ConfigError, Experiment, Project, load_experiment, load_project
from .model import PIPELINE_PLACEHOLDER
from .project import PROJECT_DIR_ENV
from .run import ALL, RunError, run_experiment
from .validate import validate_project

# 可以从哪一步开始；validate 总是运行
STAGES = ("run", "integrate", "scripts")
_INTEGRATION_KEYS = {"integration", "integration_dir"}


@dataclass
class Script:
    name: str
    tokens: list[str]     # 拆分后的命令；{project_dir} 与 {experiment} 已替换

    @property
    def needs_integration(self) -> bool:
        return any(set(PIPELINE_PLACEHOLDER.findall(t)) & _INTEGRATION_KEYS for t in self.tokens)


@dataclass
class ExperimentPlan:
    experiment: Experiment
    integrate: bool       # 实验启用了整合
    scripts: list[Script]


@dataclass
class PipelineResult:
    """一个实验的产出：单组学与整合结果包目录（新建的或沿用的；未运行时为 None）与运行过的脚本。"""
    experiment: str
    run_dir: Path | None = None
    integration_dir: Path | None = None
    scripts: list[str] = field(default_factory=list)


def _log(text: str) -> None:
    # 立即输出，保证与脚本子进程的输出按先后顺序出现
    print(text, flush=True)


def _fill(text: str, values: dict[str, str]) -> str:
    """替换 ``values`` 中有的占位符，其余原样保留。"""
    return PIPELINE_PLACEHOLDER.sub(lambda m: values.get(m.group(1), m.group(0)), text)


def _is_python(tokens: list[str]) -> bool:
    return tokens[0].lower().endswith(".py")


def _join(tokens: list[str]) -> str:
    return " ".join(f'"{t}"' if not t or " " in t else t for t in tokens)


def latest_integration(results_root: str | Path, experiment_id: str, source: str) -> Path | None:
    """由单组学结果包 ``source`` 生成的最近一次整合结果包（运行清单的 ``source_run``）。"""
    base = Path(results_root) / experiment_id
    found = [p.parent for p in sorted(base.glob(f"*/{MANIFEST_FILE}"))
             if (p.parent / "integration").is_dir() and read_manifest(p.parent).parameters.get("source_run") == source]
    return found[-1] if found else None


def _targets(project: Project, requested: list[str] | None) -> list[str]:
    known = [d.name for d in project.experiment_dirs()]
    chosen = list(dict.fromkeys(requested or project.config.pipeline.experiments or known))
    if not chosen:
        raise RunError(f"{project.path('experiments')} 下没有找到任何实验")
    unknown = [e for e in chosen if e not in known]
    if unknown:
        raise RunError(f"找不到实验 {'、'.join(unknown)}；本项目的实验：{'、'.join(known) or '无'}")
    for spec in project.config.pipeline.scripts:
        missing = [e for e in spec.experiments if e not in known]
        if missing:
            raise RunError(f"pipeline.scripts 中 {spec.command!r} 的 experiments 含有本项目没有的实验："
                           f"{'、'.join(missing)}")
    return chosen


def _plan(project: Project, targets: list[str]) -> list[ExperimentPlan]:
    """读取实验并展开脚本；脚本文件不存在、需要整合结果包但实验未启用整合时一并报出。"""
    plans, problems = [], []
    for exp_id in targets:
        experiment = load_experiment(project, exp_id)
        settings = experiment.config.integration
        integrate = settings.enabled and bool(settings.modes)
        known = {"project_dir": str(project.root), "experiment": exp_id}
        scripts = []
        for spec in project.config.pipeline.scripts:
            if spec.experiments and exp_id not in spec.experiments:
                continue
            tokens = [_fill(t, known) for t in shlex.split(spec.command)]
            script = Script(spec.name or Path(tokens[0]).stem, tokens)
            path = project.resolve(tokens[0])
            if _is_python(tokens) and not PIPELINE_PLACEHOLDER.search(tokens[0]) and not path.is_file():
                problems.append(f"脚本 {script.name}：找不到 {path}")
            if script.needs_integration and not integrate:
                problems.append(f"脚本 {script.name} 用到整合结果包（{{integration}} / {{integration_dir}}），"
                                f"但实验 {exp_id} 没有启用整合（integration.enabled 与 modes）")
            scripts.append(script)
        plans.append(ExperimentPlan(experiment, integrate, scripts))
    if problems:
        raise RunError("pipeline 配置有误，没有开始运行：\n" + "\n".join(f"  - {p}" for p in dict.fromkeys(problems)))
    return plans


def _run_script(project: Project, exp_id: str, script: Script, tokens: list[str]) -> None:
    argv = [sys.executable, str(project.resolve(tokens[0])), *tokens[1:]] if _is_python(tokens) else tokens
    env = {**os.environ, PROJECT_DIR_ENV: str(project.root)}
    try:
        code = subprocess.run(argv, cwd=project.root, env=env, check=False).returncode
    except OSError as exc:
        raise RunError(f"[{exp_id}] 脚本 {script.name} 无法启动：{exc}") from exc
    if code:
        raise RunError(f"[{exp_id}] 脚本 {script.name} 失败（退出码 {code}）")


def _run_experiment(project: Project, plan: ExperimentPlan, omics: list[str], start: str, source: str | None,
                    dry_run: bool, command: list[str] | None, log: Callable[[str], None]) -> PipelineResult:
    from omicsweave.integration.integrate import integrate_experiment, source_run

    exp_id = plan.experiment.id
    result = PipelineResult(exp_id)
    results_root = project.path("results")

    def step(stage: str, text: str) -> None:
        log(f"[{exp_id}] {stage}：{text}")

    def resume(stage: str) -> str:
        return (f"修正后续跑：omicsweave --project-dir {project.root} pipeline --experiment {exp_id} "
                f"--from {stage} --source {result.run_dir.name}")

    if start == "run":
        step("run", f"omicsweave run --experiment {exp_id} --omics {','.join(omics)}")
        if not dry_run:
            try:
                result.run_dir, notes = run_experiment(exp_id, omics, root=project.root, command=command)
            except RunError as exc:
                raise RunError(f"[{exp_id}] run 失败：{exc}") from exc
            for note in notes:
                log(f"[提示] {note}")
            step("run", f"结果包 {result.run_dir}")
    else:
        result.run_dir = source_run(results_root, exp_id, source)
        step("run", f"沿用结果包 {result.run_dir.name}")
    run_name = result.run_dir.name if result.run_dir else "<run>"

    if not plan.integrate:
        step("integrate", "实验没有启用整合（integration.enabled 与 modes），跳过")
    elif start != "scripts":
        step("integrate", f"omicsweave integrate --experiment {exp_id} --source {run_name}")
        if not dry_run:
            try:
                result.integration_dir, notes = integrate_experiment(exp_id, source=run_name, root=project.root,
                                                                     command=command)
            except RunError as exc:
                raise RunError(f"[{exp_id}] integrate 失败：{exc}\n{resume('integrate')}") from exc
            for note in notes:
                log(f"[提示] {note}")
            step("integrate", f"结果包 {result.integration_dir}")
    else:
        result.integration_dir = latest_integration(results_root, exp_id, run_name)
        if result.integration_dir is None and any(s.needs_integration for s in plan.scripts):
            raise RunError(f"[{exp_id}] 没有由 {run_name} 生成的整合结果包，请从 integrate 开始："
                           f"--from integrate --source {run_name}")
        step("integrate", f"沿用结果包 {result.integration_dir.name}" if result.integration_dir
             else f"没有由 {run_name} 生成的整合结果包")

    if not plan.scripts:
        step("scripts", "没有适用于本实验的课题脚本（project.yaml 的 pipeline.scripts）")
    integration = result.integration_dir
    values = {
        "project_dir": str(project.root), "experiment": exp_id, "run": run_name,
        "run_dir": str(result.run_dir) if result.run_dir else "<run_dir>",
        "integration": integration.name if integration else "<integration>",
        "integration_dir": str(integration) if integration else "<integration_dir>",
    }
    for script in plan.scripts:
        tokens = [_fill(t, values) for t in script.tokens]
        step(script.name, _join(["python", *tokens] if _is_python(tokens) else tokens))
        if dry_run:
            continue
        try:
            _run_script(project, exp_id, script, tokens)
        except RunError as exc:
            raise RunError(f"{exc}\n{resume('scripts')}") from exc
        result.scripts.append(script.name)
    return result


def run_pipeline(root: str | Path | None = None, experiments: list[str] | None = None,
                 omics: list[str] | None = None, *, start: str = "run", source: str | None = None,
                 dry_run: bool = False, command: list[str] | None = None,
                 log: Callable[[str], None] = _log) -> list[PipelineResult]:
    """运行标准流程，返回每个实验的产出。

    ``start`` 为 ``integrate`` 或 ``scripts`` 时沿用已有的单组学结果包（``source``，默认最近一次）；从 ``scripts``
    开始时整合结果包取由该单组学结果包生成的最近一次。``dry_run`` 只做检查并列出各步，不运行、不写文件。
    配置有误或校验未通过时抛出 :class:`RunError`，一步也不运行。
    """
    if start not in STAGES:
        raise RunError(f"未知的起始步骤 {start!r}，可选：{'、'.join(STAGES)}")
    if source and start == "run":
        raise RunError("--source 只用于 --from integrate / scripts（从 run 开始会新建单组学结果包）")
    try:
        project = load_project(root)
    except ConfigError as exc:
        raise RunError("\n".join(exc.messages)) from exc
    targets = _targets(project, experiments)
    if source and len(targets) > 1:
        raise RunError(f"--source 是某一个实验的运行 ID，请用 --experiment 只选一个实验（当前：{'、'.join(targets)}）")

    log(f"项目 {project.root}；实验 {'、'.join(targets)}；从 {start} 开始{'（只列出步骤，不运行）' if dry_run else ''}")
    report = validate_project(project.root, targets)
    for line in report.format().splitlines():
        log(f"[validate] {line}")
    if not report.ok:
        raise RunError("校验未通过（见上），没有开始运行")
    plans = _plan(project, targets)

    results: list[PipelineResult] = []
    for plan in plans:
        try:
            results.append(_run_experiment(project, plan, omics or [ALL], start, source, dry_run, command, log))
        except RunError as exc:
            done = "、".join(r.experiment for r in results)
            raise RunError(f"{exc}" + (f"\n已完成的实验：{done}" if done else "")) from exc
    return results
