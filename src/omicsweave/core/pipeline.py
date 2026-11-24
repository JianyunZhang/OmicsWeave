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
续跑沿用已有结果包时，核对它们登记的输入是否仍与当前一致：单组学结果包对照数据、样本表、参考与 experiment.yaml
中 integration 以外的各节；整合结果包对照 experiment.yaml、样本表与 panels.yaml。有改动时拒绝沿用并说明从哪一步重跑。

脚本在项目目录下运行，环境变量 ``OMICSWEAVE_PROJECT_DIR`` 指向项目目录；``.py`` 脚本用运行 omicsweave 的
同一个 Python 解释器执行。配置格式见 :class:`omicsweave.core.model.PipelineScript`。
"""
from __future__ import annotations

import json
import os
import shlex
import subprocess
import sys
from contextlib import nullcontext
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import yaml

from omicsweave.contract.manifest import MANIFEST_FILE, read_manifest

from .config import EXPERIMENT_FILE, ConfigError, Experiment, Project, load_experiment, load_project
from .delivery import Delivery, validate_delivery_tag
from .model import PIPELINE_PLACEHOLDER, ExperimentConfig, IntegrationMode, Omics
from .provenance import sha256_file
from .project import PROJECT_DIR_ENV
from .retention import CleanupReport, RetentionPolicy, apply_cleanup, format_bytes
from .run import ALL, METADATA_DIR, RunError, run_experiment, select_omics
from .validate import validate_project

# 可以从哪一步开始；validate 总是运行
STAGES = ("run", "integrate", "scripts")
_INTEGRATION_KEYS = {"integration", "integration_dir"}


@dataclass
class Script:
    name: str
    tokens: list[str]     # 拆分后的命令；{project_dir} 与 {experiment} 已替换
    outputs: list[str]

    @property
    def needs_integration(self) -> bool:
        return any(set(PIPELINE_PLACEHOLDER.findall(t)) & _INTEGRATION_KEYS for t in self.tokens)


@dataclass
class ExperimentPlan:
    experiment: Experiment
    integrate: bool       # 实验启用了整合
    scripts: list[Script]
    source: Path | None = None
    integration_source: Path | None = None


@dataclass
class PipelineResult:
    """一个实验的产出：单组学与整合结果包目录（新建的或沿用的；未运行时为 None）与运行过的脚本。"""
    experiment: str
    run_dir: Path | None = None
    integration_dir: Path | None = None
    scripts: list[str] = field(default_factory=list)
    delivery_dir: Path | None = None
    latest_pointer: Path | None = None
    cleanup_report: CleanupReport | None = None


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


def _changed_inputs(project: Project, bundle: Path, skip: frozenset[Path] = frozenset(),
                    *, ignore_kinds: frozenset[str] = frozenset()) -> list[str]:
    """结果包登记的项目内输入文件（实验配置、样本表、数据、参考、panels.yaml）中，内容与当前不同或已不存在的文件。

    结果根目录下的输入（源结果包中的表）按运行 ID 固定，不再核对。
    """
    results = project.path("results").resolve()
    changed = []
    for item in read_manifest(bundle).files:
        if item.role != "input" or item.kind in ignore_kinds:
            continue
        path = (project.root / item.path).resolve()
        if path in skip or path.is_relative_to(results):
            continue
        if not path.is_file() or path.stat().st_size != item.size or sha256_file(path) != item.sha256:
            changed.append(item.path)
    return changed


def _analysis_config_changes(experiment: Experiment, bundle: Path) -> list[str]:
    """``experiment.yaml`` 中影响单组学分析的部分（除 ``integration`` 外各节）与结果包运行时相比改动过的节。

    两边都按当前模型规范化后再比较，模型新增带默认值的字段不算改动；旧配置已无法按当前模型读取时视为全部改动。
    """
    current = experiment.config.model_dump(mode="json")
    try:
        stored = json.loads((bundle / METADATA_DIR / "experiment.json").read_text(encoding="utf-8"))
        stored = ExperimentConfig.model_validate(stored).model_dump(mode="json")
    except (OSError, ValueError):
        return [key for key in current if key != "integration"]
    return [key for key in current if key != "integration" and current[key] != stored.get(key)]


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


def _plan(project: Project, targets: list[str], omics: list[str], start: str, source: str | None,
          integration_run: str | None) -> list[ExperimentPlan]:
    """读取实验并展开脚本；脚本文件不存在、需要整合结果包但实验未启用整合时一并报出。"""
    from omicsweave.integration.integrate import select_modes, source_run

    plans, problems = [], []
    for exp_id in targets:
        experiment = load_experiment(project, exp_id)
        settings = experiment.config.integration
        integrate = settings.enabled and bool(settings.modes)
        modes = select_modes([], settings.modes) if integrate else []
        results_root = project.path("results")
        bundle = source_run(results_root, exp_id, source) if start != "run" else None
        previous = None
        if bundle is None:
            selected, _ = select_omics(experiment, omics)
            for layer in selected:
                omics_config = experiment.config.omics.settings(layer)
                spec = omics_config.input
                if spec.level == "raw" or spec.format == "salmon":
                    problems.append(f"{exp_id} 的 {layer.value} 输入格式 {spec.format} 尚未接入统一运行")
                if spec.format == "counts_matrix" and not spec.files and not getattr(spec, "tissue_files", None):
                    problems.append(f"{exp_id} 的 counts_matrix 需要至少一张包含样本的计数矩阵")
                if spec.format == "kallisto":
                    mapping = getattr(spec, "tx2gene", None) or omics_config.adapter.get("transcriptome_fasta")
                    if not mapping or not project.resolve(mapping).is_file():
                        problems.append(f"{exp_id} 的 kallisto 需要可读的 input.tx2gene 或 adapter.transcriptome_fasta")
            differential = set(selected)
            enrichment = differential - {Omics.MFA}   # MFA 运行器仅写差异通量
        else:
            differential = {Omics(p.parent.parent.name) for p in bundle.glob("omics/*/tables/differential.parquet")}
            enrichment = {Omics(p.parent.parent.name) for p in bundle.glob("omics/*/tables/enrichment.parquet")}
            bundle_manifest = read_manifest(bundle)
            bundle_steps = bundle_manifest.parameters.get("steps")
            if bundle_steps == ALL or (isinstance(bundle_steps, list) and "enrichment" in bundle_steps):
                enrichment |= {o for o in bundle_manifest.omics if o != Omics.MFA}
            # 沿用的单组学结果包须与当前数据、样本表和分析设置一致；只改 integration 一节或整合输入（如 panels）不影响它
            stale = _changed_inputs(
                project, bundle,
                frozenset({(experiment.root / EXPERIMENT_FILE).resolve()}),
                ignore_kinds=frozenset({"integration"}),
            )
            stale += [f"experiment.yaml 的 {key} 一节" for key in _analysis_config_changes(experiment, bundle)]
            if stale:
                problems.append(f"{exp_id} 的单组学结果包 {bundle.name} 已过期：{'、'.join(stale)} 在该次运行之后改过；"
                                f"请从 run 开始（不加 --from / --source）")
        if start != "scripts":
            if IntegrationMode.FLUX_ANCHORED in modes and (Omics.MFA not in differential or len(differential) < 2):
                problems.append(f"{exp_id} 的 flux_anchored 需要 MFA 与至少一个其他组学的差异表")
            if IntegrationMode.PATHWAY_CONVERGENCE in modes and len(enrichment) < 2:
                problems.append(f"{exp_id} 的 pathway_convergence 需要至少两个组学的富集结果")
            if {IntegrationMode.REVERSAL, IntegrationMode.REACTION_DOSSIER}.intersection(modes) and not differential:
                problems.append(f"{exp_id} 的源结果包缺少差异表")
            panel = settings.parameters.get(IntegrationMode.REACTION_DOSSIER, {}).get("panels")
            if IntegrationMode.REACTION_DOSSIER in modes and panel:
                from omicsweave.integration.dossier.panels import PanelError, load_panels
                try:
                    load_panels(experiment.root / panel)
                except (PanelError, OSError, ValueError, yaml.YAMLError) as exc:
                    problems.append(str(exc))
        elif integrate and bundle:
            if integration_run:
                previous = results_root / exp_id / integration_run
                manifest = read_manifest(previous) if (previous / MANIFEST_FILE).is_file() else None
                if (not (previous / "integration").is_dir() or manifest is None or manifest.project_id != project.id
                        or manifest.experiment_id != exp_id or manifest.parameters.get("source_run") != bundle.name):
                    problems.append(f"整合包 {integration_run} 不是由单组学结果包 {bundle.name} 生成的完整整合包")
            else:
                previous = latest_integration(results_root, exp_id, bundle.name)
            # 沿用的整合结果包须与当前 experiment.yaml、样本表与 panels.yaml 一致，否则交付的面板矩阵是旧的
            stale = _changed_inputs(project, previous) if previous and (previous / MANIFEST_FILE).is_file() else []
            if stale:
                problems.append(f"{exp_id} 的整合结果包 {previous.name} 已过期：{'、'.join(stale)} 在整合之后改过；"
                                f"请从 integrate 开始：--from integrate --source {bundle.name}")
        elif integration_run:
            problems.append(f"{exp_id} 未启用整合，无法沿用整合包 {integration_run}")
        known = {"project_dir": str(project.root), "experiment": exp_id}
        scripts = []
        for spec in project.config.pipeline.scripts:
            if spec.experiments and exp_id not in spec.experiments:
                continue
            tokens = [_fill(t, known) for t in shlex.split(spec.command)]
            outputs = [_fill(path, {"experiment": exp_id}) for path in spec.outputs]
            script = Script(spec.name or Path(tokens[0]).stem, tokens, outputs)
            path = project.resolve(tokens[0])
            if _is_python(tokens) and not PIPELINE_PLACEHOLDER.search(tokens[0]) and not path.is_file():
                problems.append(f"脚本 {script.name}：找不到 {path}")
            if script.needs_integration and not integrate:
                problems.append(f"脚本 {script.name} 用到整合结果包（{{integration}} / {{integration_dir}}），"
                                f"但实验 {exp_id} 没有启用整合（integration.enabled 与 modes）")
            scripts.append(script)
        if start == "scripts" and bundle and previous is None and any(s.needs_integration for s in scripts):
            problems.append(f"{exp_id} 没有由 {bundle.name} 生成的整合结果包，请从 integrate 开始")
        plans.append(ExperimentPlan(experiment, integrate, scripts, bundle, previous))
    if problems:
        raise RunError("pipeline 配置有误，没有开始运行：\n" + "\n".join(f"  - {p}" for p in dict.fromkeys(problems)))
    return plans


def _run_script(project: Project, exp_id: str, script: Script, tokens: list[str]) -> None:
    argv = [sys.executable, str(project.resolve(tokens[0])), *tokens[1:]] if _is_python(tokens) else tokens
    env = {**os.environ, PROJECT_DIR_ENV: str(project.root)}
    # 从工作树直接启动时，子解释器也应加载同一份 omicsweave。
    env["PYTHONPATH"] = os.pathsep.join(filter(None, [str(Path(__file__).resolve().parents[2]), env.get("PYTHONPATH")]))
    try:
        code = subprocess.run(argv, cwd=project.root, env=env, check=False).returncode
    except OSError as exc:
        raise RunError(f"[{exp_id}] 脚本 {script.name} 无法启动：{exc}") from exc
    if code:
        raise RunError(f"[{exp_id}] 脚本 {script.name} 失败（退出码 {code}）")


def _run_experiment(project: Project, plan: ExperimentPlan, omics: list[str], start: str,
                    dry_run: bool, command: list[str] | None, log: Callable[[str], None],
                    warnings: list[str], *, tag: str | None = None, pin: bool = False,
                    retention: RetentionPolicy | None = None, run_gsea: bool | None = None) -> PipelineResult:
    from omicsweave.integration.integrate import integrate_experiment

    exp_id = plan.experiment.id
    result = PipelineResult(exp_id, plan.source, plan.integration_source)
    script_keys = [f"script_{i + 1}:{script.name}" for i, script in enumerate(plan.scripts)]
    delivery = None
    if not dry_run:
        try:
            delivery = Delivery(plan.experiment, script_keys, command or [], tag=tag, pinned=pin)
        except (OSError, ValueError) as exc:
            raise RunError(f"无法在 {project.path('results')} 初始化 {exp_id} 的交付目录：{exc}") from exc
        result.delivery_dir = delivery.directory
        delivery.data["validation_warnings"] = warnings
        for key, script in zip(script_keys, plan.scripts):
            delivery.stages[key]["expected_outputs"] = script.outputs
        log(f"[{exp_id}] 本次交付：{delivery.directory}")

    def stage(key: str, success: str = "completed", **details):
        return delivery.stage(key, success, **details) if delivery else nullcontext()

    def step(name: str, message: str) -> None:
        log(f"[{exp_id}] {name}：{message}")

    def resume(at: str) -> str:
        args = ["omicsweave", "--project-dir", str(project.root), "pipeline", "--experiment", exp_id,
                "--results-dir", str(project.path("results")), "--from", at, "--source", result.run_dir.name]
        if at == "scripts" and result.integration_dir:
            args += ["--integration-run", result.integration_dir.name]
        return _join(args)

    current = "run"
    try:
        with stage(current, "reused" if plan.source else "completed"):
            if plan.source:
                step(current, f"沿用结果包 {plan.source.name}")
            else:
                gsea_flag = " --gsea" if run_gsea is True else (" --no-gsea" if run_gsea is False else "")
                step(current, f"omicsweave run --experiment {exp_id} --omics {','.join(omics)}{gsea_flag}")
                if delivery:
                    result.run_dir, notes = run_experiment(
                        exp_id, omics, root=project, command=command,
                        run_id=delivery.directory.name, into_existing=True, run_gsea=run_gsea,
                    )
                    for note in notes:
                        log(f"[提示] {note}")
                    delivery.stages[current]["notes"] = notes
            if delivery:
                delivery.bundle(current, result.run_dir)

        run_name = result.run_dir.name if result.run_dir else "<run>"
        current = "integrate"
        status = "completed" if plan.integrate else "skipped"
        if start == "scripts":
            status = "reused" if plan.integration_source else "skipped"
        with stage(current, status):
            if not plan.integrate:
                step(current, "实验没有启用整合，跳过")
            elif start == "scripts":
                step(current, f"沿用结果包 {result.integration_dir.name}" if result.integration_dir else "没有关联整合包")
            else:
                step(current, f"omicsweave integrate --experiment {exp_id} --source {run_name}")
                if delivery:
                    if plan.source is None:
                        result.integration_dir, notes = integrate_experiment(
                            exp_id, source=run_name, root=project, command=command, into_source=True,
                        )
                    else:
                        result.integration_dir, notes = integrate_experiment(
                            exp_id, source=run_name, root=project, command=command,
                            run_id=delivery.directory.name, into_existing=True,
                        )
                    for note in notes:
                        log(f"[提示] {note}")
                    delivery.stages[current]["notes"] = notes
            if delivery and result.integration_dir:
                delivery.bundle(current, result.integration_dir)

        integration = result.integration_dir
        values = {"experiment": exp_id, "project_dir": str(project.root), "run": run_name,
                  "run_dir": str(result.run_dir) if result.run_dir else "<run_dir>",
                  "integration": integration.name if integration else "<integration>",
                  "integration_dir": str(integration) if integration else "<integration_dir>",
                  "delivery_dir": str(result.delivery_dir) if delivery else "<delivery_dir>"}
        for current, script in zip(script_keys, plan.scripts):
            tokens = [_fill(t, values) for t in script.tokens]
            step(script.name, _join(["python", *tokens] if _is_python(tokens) else tokens))
            if delivery:
                with stage(current, command=tokens, expected_outputs=script.outputs):
                    executable = project.resolve(tokens[0])
                    if executable.is_file():
                        delivery.stages[current]["script_sha256"] = sha256_file(executable)
                    _run_script(project, exp_id, script, tokens)
                    missing = [name for name in script.outputs if not (delivery.directory / name).is_file()]
                    if missing:
                        raise RunError(f"脚本 {script.name} 缺少预期交付文件：{'、'.join(missing)}")
                    result.scripts.append(script.name)
        current = "summary"
        step(current, "复制面板矩阵与图件到 panels/、通路总图到 pathway_map/，写无面板候选表 candidates/，"
                      "生成 summary.xlsx 与 index.json，并更新 latest.json 指针")
        if delivery:
            with stage(current):
                delivery.write_summary()
                delivery.collect_exports()
            result.latest_pointer = delivery.complete()
        if retention is not None and retention.enabled:
            if dry_run:
                step("retention", f"完成后按保留策略自动清理（{retention.describe()}）")
            elif delivery:
                cleanup = apply_cleanup(
                    project.path("results"), retention, [exp_id], dry_run=False,
                    protected_deliveries={exp_id: {delivery.directory.name}},
                )
                result.cleanup_report = cleanup
                if cleanup.removed:
                    step("retention", f"已按保留策略清理 {len(cleanup.removed)} 项，"
                                      f"释放 {format_bytes(cleanup.reclaimed_bytes)}")
                else:
                    step("retention", "已检查保留策略，无需清理")
                for note in cleanup.notes:
                    log(f"[提示] {note}")
                for item, error in cleanup.failed:
                    log(f"[警告] 保留策略未能删除 {item.path}：{error}")
    except Exception as exc:
        message = f"[{exp_id}] {current} 失败：{exc}"
        if result.delivery_dir:
            message += f"；交付记录 {result.delivery_dir / 'index.json'}"
        if result.run_dir:
            message += "\n修正后续跑：" + resume("integrate" if current == "integrate" else "scripts")
        raise RunError(message) from exc
    return result


def run_pipeline(root: str | Path | None = None, experiments: list[str] | None = None,
                 omics: list[str] | None = None, *, start: str = "run", source: str | None = None,
                 results_dir: str | Path | None = None, integration_run: str | None = None,
                 tag: str | None = None, pin: bool = False,
                 retention: RetentionPolicy | None = None, run_gsea: bool | None = None,
                 dry_run: bool = False, command: list[str] | None = None,
                 log: Callable[[str], None] = _log) -> list[PipelineResult]:
    """运行标准流程，返回每个实验的产出。

    ``start`` 为 ``integrate`` 或 ``scripts`` 时沿用已有的单组学结果包（``source``，默认最近一次）；从 ``scripts``
    开始时整合结果包取由该单组学结果包生成的最近一次。``dry_run`` 只做检查并列出各步，不运行、不写文件。
    ``retention`` 未显式传入时读取 ``project.yaml`` 的 ``pipeline.retention``（默认关闭，不自动清理）。
    配置有误或校验未通过时抛出 :class:`RunError`，一步也不运行。
    """
    if start not in STAGES:
        raise RunError(f"未知的起始步骤 {start!r}，可选：{'、'.join(STAGES)}")
    if source and start == "run":
        raise RunError("--source 只用于 --from integrate / scripts（从 run 开始会新建单组学结果包）")
    if tag is not None:
        try:
            tag = validate_delivery_tag(tag)
        except ValueError as exc:
            raise RunError(str(exc)) from exc
    try:
        project = load_project(root)
    except ConfigError as exc:
        raise RunError("\n".join(exc.messages)) from exc
    if results_dir is not None:
        project.config.project.paths.results = str(project.resolve(results_dir).resolve())
    if integration_run:
        if Path(integration_run).name != integration_run or integration_run in (".", ".."):
            raise RunError("--integration-run 应为本实验的运行 ID")
        if start != "scripts" or not source:
            raise RunError("--integration-run 需要 --from scripts 和 --source")
    targets = _targets(project, experiments)
    if source and len(targets) > 1:
        raise RunError(f"--source 是某一个实验的运行 ID，请用 --experiment 只选一个实验（当前：{'、'.join(targets)}）")

    effective_retention = retention if retention is not None else RetentionPolicy.from_settings(project.config.pipeline.retention)
    log(f"项目 {project.root}；结果根目录 {project.path('results')}；实验 {'、'.join(targets)}；"
        f"从 {start} 开始{'（只列出步骤，不运行）' if dry_run else ''}")
    report = validate_project(project, targets)
    for line in report.format().splitlines():
        log(f"[validate] {line}")
    if not report.ok:
        raise RunError("校验未通过（见上），没有开始运行")
    plans = _plan(project, targets, omics or [ALL], start, source, integration_run)

    results: list[PipelineResult] = []
    for plan in plans:
        try:
            warnings = [f"{issue.where}: {issue.message}" for issue in report.warnings
                        if issue.where in ("项目", f"实验 {plan.experiment.id}")]
            results.append(_run_experiment(
                project, plan, omics or [ALL], start, dry_run, command, log, warnings,
                tag=tag, pin=pin, retention=effective_retention, run_gsea=run_gsea,
            ))
        except RunError as exc:
            done = "、".join(r.experiment for r in results)
            raise RunError(f"{exc}" + (f"\n已完成的实验：{done}" if done else "")) from exc
    return results
