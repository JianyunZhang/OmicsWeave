"""统一执行前台分析与后台 Worker 的请求。"""
from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, TypeAdapter, model_validator

from omicsweave.core.progress import NULL_REPORTER, ProgressReporter
from omicsweave.core.path_scope import PathScope, identifier
from omicsweave.core.overrides import TRIALS_DIR, PreparedTrial, TrialOptions, changed_source_outputs, prepare_trial

from .schemas import ErrorCode, ServiceError


class ExecutionRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    kind: Literal["run", "integrate", "pipeline"]
    experiments: list[str]
    omics: list[str] | None = None
    steps: list[str] | None = None
    modes: list[str] | None = None
    source: str | None = None
    start: Literal["run", "integrate", "scripts"] = "run"
    integration_run: str | None = None
    results_dir: Path | None = None
    tag: str | None = None
    pin: bool = False
    dry_run: bool = False
    run_gsea: bool | None = None
    command: list[str] = Field(default_factory=list)
    retention: dict | None = None
    allowed_data_roots: list[Path] | None = None  # None 为原有本机 CLI，列表为受限 MCP 请求
    allow_custom_scripts: bool = True
    trial: TrialOptions | None = None

    @model_validator(mode="after")
    def _one_experiment(self):
        if self.kind != "pipeline" and len(self.experiments) != 1:
            raise ValueError("独立组学运行或整合需要恰好一个实验")
        if self.trial is not None and len(self.experiments) != 1:
            raise ValueError("隔离试验需要显式指定一个实验")
        return self


def _scoped_request(project_dir: Path, request: ExecutionRequest):
    import shlex

    from omicsweave.contract.manifest import read_manifest
    from omicsweave.core.path_scope import ScopeError

    scope = PathScope.create(project_dir, request.allowed_data_roots or ())
    project = scope.project()
    if request.kind == "pipeline" and request.allow_custom_scripts and not request.trial:
        for script in project.config.pipeline.scripts:
            path = scope.path(shlex.split(script.command)[0])
            if path.suffix.lower() != ".py" or not path.is_file():
                raise ScopeError("仅接受显式授权的项目内 Python 课题脚本")
    root = scope.path(request.results_dir or project.path("results"))
    if request.trial:
        scope.path(root / TRIALS_DIR)
    if request.source:
        identifier(request.source)
    if request.integration_run:
        identifier(request.integration_run)
    targets = request.experiments or project.config.pipeline.experiments or [p.name for p in project.experiment_dirs()]
    for exp_id in targets:
        scope.experiment(project, exp_id, execution=request.kind == "run" or (
            request.kind == "pipeline" and request.start == "run"
        ))
        base = scope.tree(root / identifier(exp_id))
        if (request.kind == "pipeline" and request.start != "run") or (request.trial and request.kind == "integrate"):
            if request.trial:
                ids = {value for value in (request.source, request.trial.parent_run_id) if value}
                manifests = [scope.path(base / identifier(value) / "manifest.json") for value in sorted(ids)]
            else:
                manifests = base.glob("*/manifest.json")
            for path in manifests:
                for record in read_manifest(path.parent).files:
                    if record.role == "input":
                        scope.path(record.path, external=True)
    return scope, root


def execute(project: Path, request: ExecutionRequest, reporter: ProgressReporter = NULL_REPORTER) -> dict:
    if request.allowed_data_roots is not None or request.trial:
        from omicsweave.core.genesets import offline_libraries

        if request.allowed_data_roots is not None:
            _scoped_request(project, request)
        with offline_libraries():
            return _execute(project, request, reporter)
    return _execute(project, request, reporter)


def _execute(project: Path, request: ExecutionRequest, reporter: ProgressReporter) -> dict:
    reporter.check_cancelled()
    trial = _prepare_trial(project, request, reporter) if request.trial else None
    trial_data = {"trial": {"results_dir": str(trial.root), **trial.lineage.model_dump(mode="json")}} if trial else {}
    if request.kind == "pipeline":
        from omicsweave.core.pipeline import run_pipeline
        from omicsweave.core.retention import RetentionPolicy

        operations = []
        results = run_pipeline(
            project, request.experiments, request.omics, start=request.start, source=request.source,
            results_dir=request.results_dir, integration_run=request.integration_run, tag=request.tag,
            pin=request.pin, dry_run=request.dry_run, command=request.command, run_gsea=request.run_gsea,
            retention=RetentionPolicy(**request.retention) if request.retention is not None else None,
            log=operations.append if request.dry_run else lambda message: print(message, flush=True),
            reporter=reporter,
            allow_custom_scripts=request.allow_custom_scripts,
            trial=trial,
        )
        # dry_run 按现有 pipeline 预检生成计划，路径可能是待生成的占位对象。
        return {
            "dry_run": request.dry_run,
            "deliveries": TypeAdapter(list[dict]).dump_python([asdict(result) for result in results], mode="json"),
            **trial_data,
            **({"plan": operations} if request.dry_run else {}),
        }

    experiment = request.experiments[0]
    root = project
    if request.results_dir is not None and not trial:
        from omicsweave.core.config import load_project

        root = load_project(project)
        root.config.project.paths.results = str(request.results_dir)
    if request.kind == "run":
        from omicsweave.core.run import ALL, load_validated, run_experiment, select_omics

        if request.dry_run:
            exp = trial.experiment if trial else load_validated(experiment, root)[1]
            chosen, notes = select_omics(exp, request.omics or [ALL])
            return {"dry_run": True, "experiment": experiment, "omics": [layer.value for layer in chosen],
                    "steps": request.steps, "notes": notes, **trial_data}

        if trial:
            trial.reserve()
        directory, notes = run_experiment(
            experiment, request.omics or [ALL], request.steps, root=root, command=request.command,
            run_gsea=request.run_gsea, reporter=reporter,
            prepared=trial.experiment if trial else None, lineage=trial.lineage if trial else None,
            run_id=f"run_{trial.lineage.trial_id}" if trial else None,
        )
    else:
        from omicsweave.integration.integrate import integrate_experiment, select_modes, source_run

        if request.dry_run:
            from omicsweave.core.run import load_validated

            loaded, exp = (trial.experiment.project, trial.experiment) if trial else load_validated(experiment, root)
            chosen = select_modes(request.modes or [], exp.config.integration.modes)
            source = trial.source if trial else source_run(loaded.path("results"), exp.id, request.source)
            return {"dry_run": True, "experiment": experiment, "modes": [mode.value for mode in chosen],
                    "source_run_id": source.name, **trial_data}

        if trial:
            trial.reserve()
        directory, notes = integrate_experiment(
            experiment, request.modes, source=request.source, root=root,
            command=request.command, reporter=reporter,
            prepared=trial.experiment if trial else None, lineage=trial.lineage if trial else None,
            source_bundle=trial.source if trial else None,
            run_id=f"run_{trial.lineage.trial_id}" if trial else None,
        )
    result = {"experiment": experiment, "run_id": directory.name, "run_dir": str(directory), "notes": notes}
    reporter.emit("artifact", stage=request.kind, detail={
        "experiment": experiment, "run_id": directory.name, "run_dir": str(directory),
    })
    return {**result, **trial_data}


def _prepare_trial(project_dir: Path, request: ExecutionRequest, reporter: ProgressReporter) -> PreparedTrial:
    from omicsweave.core.config import load_project
    from omicsweave.core.pipeline import source_changes
    from omicsweave.core.validate import validate_experiment
    from omicsweave.integration.integrate import source_run

    project = load_project(project_dir)
    if request.results_dir is not None:
        project.config.project.paths.results = str(project.resolve(request.results_dir).resolve())
    if request.start == "scripts" or request.integration_run:
        raise ServiceError(ErrorCode.CLI_USAGE_ERROR, "隔离试验从 run 或 integrate 开始")
    reused = request.kind == "integrate" or (request.kind == "pipeline" and request.start == "integrate")
    source = (request.source or request.trial.parent_run_id) if reused else None
    if reused and not source:
        raise ServiceError(ErrorCode.CLI_USAGE_ERROR, "复用结果的试验需要明确指定来源")
    trial = prepare_trial(project, request.experiments[0], request.trial, reused_source=source, job_id=reporter.job_id)
    if request.kind != "pipeline":
        errors = [issue for issue in validate_experiment(trial.experiment) if issue.level == "error"]
        if errors:
            raise ServiceError(ErrorCode.CONFIG_SCHEMA_ERROR, "试验生效配置未通过校验",
                               details={"issues": [asdict(issue) for issue in errors]})
    if trial.source:
        source_run(trial.formal_root, trial.experiment.id, trial.source.name)
        if changed := changed_source_outputs(trial.source):
            raise ServiceError(ErrorCode.STALE_SOURCE_BUNDLE, "复用包的科学输出与正式清单不一致",
                               details={"changed_files": changed, "recommended_start": "run"})
        if request.kind == "integrate":
            stale = source_changes(trial.experiment, trial.source)
            if stale:
                raise ServiceError(ErrorCode.STALE_SOURCE_BUNDLE, "试验要求从 run 阶段开始",
                                   details={"changed_files": stale, "recommended_start": "run"})
    return trial


def submit_execution(project_dir: Path, request: ExecutionRequest) -> dict:
    from omicsweave.core.config import load_project
    from omicsweave.core.jobs import submit

    if request.dry_run:
        raise ServiceError(ErrorCode.CLI_USAGE_ERROR, "请单独运行 --dry-run；预检不会派发后台作业")
    if request.allowed_data_roots is not None:
        from omicsweave.core.file_state import file_lock

        scope = PathScope.create(project_dir, request.allowed_data_roots)
        with file_lock(scope.path(".omicsweave-config.lock")):
            _, root = _scoped_request(project_dir, request)
            request = request.model_copy(update={"results_dir": root, "retention": {"enabled": False}})
            scope.tree(root / "_jobs")
            return submit(root, scope.root, request.model_dump(mode="json"), request.experiments)
    project = load_project(project_dir)
    root = project.resolve(request.results_dir).resolve() if request.results_dir is not None else project.path("results")
    # 新会话内始终使用提交时解析的结果根目录。
    request = request.model_copy(update={"results_dir": root.resolve()})
    experiments = request.experiments or list(project.config.pipeline.experiments)
    return submit(root, project.root, request.model_dump(mode="json"), experiments)


def job_operation(project: Path, root: Path, operation: str, *, job: str | None = None,
                  experiment: str | None = None, active: bool = False, limit: int = 20,
                  tail: int = 50, since_event: int = 0, events_only: bool = False) -> dict:
    from omicsweave.core.jobs import ACTIVE, job_dir, list_status, read_logs, read_status, request_cancel

    if operation == "status" and not job:
        if not 1 <= limit <= 200:
            raise ServiceError(ErrorCode.CLI_USAGE_ERROR, "--limit 必须在 1–200")
        states = [state for state in list_status(root)
                  if state["project_dir"] == str(project)
                  and (not experiment or not state["experiments"] or experiment in state["experiments"])
                  and (not active or state["status"] in ACTIVE)]
        states.sort(key=lambda state: state["created_at"], reverse=True)
        return {"jobs": states[:limit], "total": len(states)}
    try:
        directory = job_dir(root, job)
        state = read_status(directory)
        if state["project_dir"] != str(project):
            raise FileNotFoundError(f"本项目下未找到作业 {job}")
    except FileNotFoundError as exc:
        raise ServiceError(ErrorCode.JOB_NOT_FOUND, str(exc)) from exc
    except (ValueError, KeyError) as exc:
        raise ServiceError(ErrorCode.CLI_USAGE_ERROR, str(exc)) from exc
    if operation == "status":
        return state
    if operation == "cancel":
        state = request_cancel(directory)
        if state["status"] not in ACTIVE and state["status"] != "cancelled":
            raise ServiceError(ErrorCode.JOB_CONFLICT, f"作业已停止，当前状态 {state['status']}")
        return state
    if not 1 <= tail <= 200 or since_event < 0:
        raise ServiceError(ErrorCode.CLI_USAGE_ERROR, "--tail 必须在 1–200，--since-event 必须为非负数")
    return read_logs(directory, tail, since_event, events_only)
