"""将共享应用服务映射到带类型和副作用标注的 MCP 工具。"""
from __future__ import annotations

import base64
import io
from dataclasses import asdict
from pathlib import Path
from typing import Annotated, Literal
from uuid import uuid4

import anyio
from mcp.server.fastmcp import Context, FastMCP
from mcp.types import CallToolResult, ImageContent, TextContent, ToolAnnotations
from pydantic import Field, TypeAdapter

from omicsweave.core.path_scope import identifier
from omicsweave.core.overrides import Intent, TrialOptions
from omicsweave.service import project as project_service
from omicsweave.service.execution import ExecutionRequest, execute, job_operation, submit_execution
from omicsweave.service.project import EditOptions, QCUpdate
from omicsweave.service.query import (
    Filter, ResultView, TableQuery, diagnose_view, diff_payload, query_table, run_digest, trace_entity,
)
from omicsweave.service.schemas import ErrorCode, ServiceError, ServiceResponse, classify_error, json_safe
from omicsweave.service.visualization import render_abstract, render_map

from .security import ServerPolicy

_READ = ToolAnnotations(readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=False)
_COMPUTE = ToolAnnotations(readOnlyHint=False, destructiveHint=False, idempotentHint=False, openWorldHint=False)
_WRITE = ToolAnnotations(readOnlyHint=False, destructiveHint=True, idempotentHint=False, openWorldHint=False)
_CANCEL = ToolAnnotations(readOnlyHint=False, destructiveHint=True, idempotentHint=True, openWorldHint=False)
WaitSeconds = Annotated[float, Field(ge=0, le=30)]
PageSize = Annotated[int, Field(ge=1, le=200)]
Cursor = Annotated[int, Field(ge=0)]
OmicsName = Literal["all", "mfa", "transcriptomics", "proteomics", "metabolomics"]
ModeName = Literal["flux_anchored", "pathway_convergence", "reaction_dossier", "reversal"]


def result_view(policy: ServerPolicy, experiment_id: str, run_id: str, delivery_id: str | None = None,
                trial_id: str | None = None):
    return ResultView(policy.results(experiment_id, trial_id), experiment_id, run_id, delivery_id)


def _preview_image(files: list[str]) -> ImageContent | None:
    from PIL import Image

    pngs = [path for path in files if Path(path).suffix.lower() == ".png"]
    if len(pngs) != 1:
        return None
    with Image.open(pngs[0]) as image:
        image.thumbnail((1280, 1280))
        buffer = io.BytesIO()
        image.save(buffer, format="PNG", optimize=True)
    if buffer.tell() > 512 * 1024:
        return None
    return ImageContent(type="image", data=base64.b64encode(buffer.getvalue()).decode("ascii"), mimeType="image/png")


def register_tools(server: FastMCP, policy: ServerPolicy) -> None:
    scope = policy.scope
    root = scope.root

    def wire(name: str, data=None, *, error: ServiceError | None = None, image: ImageContent | None = None):
        response = ServiceResponse(command=name, project_dir=str(root), ok=error is None, data=json_safe(data),
                                   error=error.detail if error else None)
        content = [TextContent(type="text", text=response.model_dump_json())]
        if image:
            content.append(image)
        return CallToolResult(content=content, structuredContent=response.model_dump(mode="json"), isError=error is not None)

    async def call(name, fn, *, preview_image=False):
        try:
            data = await anyio.to_thread.run_sync(fn)
            image = await anyio.to_thread.run_sync(lambda: _preview_image(data["files"])) if preview_image else None
            return wire(name, data, image=image)
        except Exception as exc:
            return wire(name, error=classify_error(exc))

    async def dispatch(name: str, ctx: Context, experiment_id: str, wait_seconds: float, dry_run: bool,
                       *, overrides=None, intent=None, rationale=None, parent_run_id=None, **kwargs):
        def start():
            trial = None
            if any(value is not None for value in (overrides, intent, rationale, parent_run_id)):
                trial = TrialOptions(overrides=overrides or {}, intent=intent or "sensitivity_check",
                                     rationale=rationale or "", parent_run_id=parent_run_id)
            request = ExecutionRequest(experiments=[identifier(experiment_id)], dry_run=dry_run, command=["mcp", name],
                                       allowed_data_roots=list(scope.data_roots),
                                       allow_custom_scripts=policy.allow_custom_scripts,
                                       retention={"enabled": False}, trial=trial, **kwargs)
            return execute(root, request) if dry_run else submit_execution(root, request)

        try:
            state = await anyio.to_thread.run_sync(start)
            if dry_run or not wait_seconds:
                return wire(name, state)
            from omicsweave.core.jobs import ACTIVE, read_status

            deadline = anyio.current_time() + wait_seconds
            state_path = Path(state["state_path"])
            last_event = 0
            while True:
                state = await anyio.to_thread.run_sync(read_status, state_path.parent)
                seq = state["last_event_seq"]
                if seq > last_event and state["status"] in ACTIVE:
                    last_event = seq
                    await ctx.report_progress(float(seq), message=(state.get("progress") or {}).get("stage"))
                remaining = deadline - anyio.current_time()
                if state["status"] not in ACTIVE or remaining <= 0:
                    break
                await anyio.sleep(min(0.2, remaining))
            return wire(name, state)
        except Exception as exc:
            return wire(name, error=classify_error(exc))

    @server.tool(annotations=_READ)
    async def list_experiments() -> CallToolResult:
        """列出授权项目中的实验。"""
        return await call("list_experiments", lambda: project_service.project_overview(scope))

    @server.tool(annotations=_READ)
    async def validate_experiment(experiment_id: str) -> CallToolResult:
        """校验设计、样本表、路径与强类型参数。"""
        def validate():
            report = project_service.validation(root, [experiment_id], scope=scope)
            return project_service.validation_summary(report)

        return await call("validate_experiment", validate)

    @server.tool(annotations=_READ)
    async def inspect_input_files(experiment_id: str | None = None,
                                  paths: list[str] | None = None) -> CallToolResult:
        """列出指定文件的存在性、类型及字节数。"""
        return await call("inspect_input_files", lambda: project_service.inspect_inputs(scope, experiment_id, paths))

    @server.tool(annotations=_WRITE)
    async def patch_experiment_config(experiment_id: str, patch: dict,
                                      target: Literal["experiment_yaml", "panels_yaml"] = "experiment_yaml",
                                      dry_run: bool = True, reason: str = "",
                                      expected_sha256: str | None = None) -> CallToolResult:
        """预览合并补丁。应用时需要启动授权、预览 SHA-256 和修改依据。"""
        def apply():
            options = EditOptions(dry_run=dry_run, reason=reason, expected_sha256=expected_sha256)
            policy.check_edit(dry_run)
            return project_service.patch_config(scope, experiment_id, patch, options, target)

        return await call("patch_experiment_config", apply)

    @server.tool(annotations=_WRITE)
    async def update_sample_qc(experiment_id: str, updates: list[QCUpdate], dry_run: bool = True,
                               reason: str = "", expected_sha256: str | None = None) -> CallToolResult:
        """预览样本标记及可靠性变化。提交需要启动授权、预览 SHA-256 和修改依据。"""
        def apply():
            options = EditOptions(dry_run=dry_run, reason=reason, expected_sha256=expected_sha256)
            policy.check_edit(dry_run)
            return project_service.update_sample_qc(scope, experiment_id, updates, options)

        return await call("update_sample_qc", apply)

    @server.tool(annotations=_WRITE)
    async def init_or_migrate_experiment(
        experiment_id: str, action: Literal["init", "migrate"] = "init",
        species: Literal["mus_musculus", "homo_sapiens"] = "mus_musculus", dry_run: bool = True,
    ) -> CallToolResult:
        """预览或写入空项目模板。旧配置迁移使用本机 CLI。"""
        def init():
            from omicsweave.core.model import Species
            from omicsweave.core.scaffold import create_project
            from omicsweave.core.file_state import file_lock

            if action == "migrate":
                raise ServiceError(ErrorCode.SECURITY_POLICY_VIOLATION, "旧版路径由研究者通过本机 migrate 预检")
            exp = Path("experiments") / identifier(experiment_id)
            files = [scope.path(name) for name in ("project.yaml", "run_analysis.py", exp / "experiment.yaml",
                                                   exp / "subjects.tsv", exp / "samples.tsv")]
            for folder in (exp, "results", "analysis", "references"):
                scope.path(folder)
            policy.check_edit(dry_run)
            if not dry_run:
                with file_lock(scope.path(".omicsweave-config.lock")):
                    if any(path.exists() for path in files):
                        raise ServiceError(ErrorCode.JOB_CONFLICT, "模板目标已存在")
                    files = create_project(root, "omicsweave_project", Species(species), experiment_id)
            return {"action": action, "dry_run": dry_run, "files_planned_or_written": [str(path) for path in files]}

        return await call("init_or_migrate_experiment", init)

    @server.tool(annotations=_COMPUTE)
    async def submit_pipeline(experiment_id: str, ctx: Context, start: Literal["run", "integrate", "scripts"] = "run",
                              source_run_id: str | None = None, integration_run_id: str | None = None,
                              omics: list[OmicsName] | None = None, tag: str | None = None, pin: bool = False,
                              dry_run: bool = True, wait_seconds: WaitSeconds = 15, overrides: dict | None = None,
                              intent: Intent | None = None, rationale: str | None = None,
                              parent_run_id: str | None = None) -> CallToolResult:
        """预检或启动后台计算。以作业 ID 持续观察阶段进度。"""
        return await dispatch("submit_pipeline", ctx, experiment_id, wait_seconds, dry_run,
                              kind="pipeline", start=start, source=source_run_id, integration_run=integration_run_id,
                              omics=omics, tag=tag, pin=pin, overrides=overrides, intent=intent,
                              rationale=rationale, parent_run_id=parent_run_id)

    @server.tool(annotations=_COMPUTE)
    async def submit_omics_run(experiment_id: str, ctx: Context, omics: list[OmicsName] | None = None,
                               steps: list[Literal["qc", "differential", "enrichment", "figures"]] | None = None,
                               dry_run: bool = True, wait_seconds: WaitSeconds = 15, overrides: dict | None = None,
                               intent: Intent | None = None, rationale: str | None = None,
                               parent_run_id: str | None = None) -> CallToolResult:
        """预检或启动独立的单组学分析。"""
        return await dispatch("submit_omics_run", ctx, experiment_id, wait_seconds, dry_run,
                              kind="run", omics=omics, steps=steps, overrides=overrides, intent=intent,
                              rationale=rationale, parent_run_id=parent_run_id)

    @server.tool(annotations=_COMPUTE)
    async def submit_integration(experiment_id: str, ctx: Context, source_run_id: str | None = None,
                                 modes: list[ModeName] | None = None, dry_run: bool = True,
                                 wait_seconds: WaitSeconds = 15, overrides: dict | None = None,
                                 intent: Intent | None = None, rationale: str | None = None,
                                 parent_run_id: str | None = None) -> CallToolResult:
        """用已有单组学结果包预检或启动整合。"""
        return await dispatch("submit_integration", ctx, experiment_id, wait_seconds, dry_run,
                              kind="integrate", source=source_run_id, modes=modes, overrides=overrides,
                              intent=intent, rationale=rationale, parent_run_id=parent_run_id)

    @server.tool(annotations=_READ)
    async def get_job_status(job_id: str, since_seq: Cursor = 0, tail: PageSize = 20) -> CallToolResult:
        """读取作业状态及增量进度、日志。"""
        def status():
            directory = policy.job_root()
            state = job_operation(root, directory, "status", job=job_id)
            logs = job_operation(root, directory, "logs", job=job_id, since_event=since_seq, tail=tail)
            return {**state, "log": logs}

        return await call("get_job_status", status)

    @server.tool(annotations=_CANCEL)
    async def cancel_job(job_id: str) -> CallToolResult:
        """请求在安全检查点取消作业。"""
        return await call("cancel_job", lambda: job_operation(root, policy.job_root(), "cancel", job=job_id))

    @server.tool(annotations=_READ)
    async def get_run_digest(experiment_id: str, run_id: str = "latest",
                             delivery_id: str | None = None, trial_id: str | None = None) -> CallToolResult:
        """列出结果包的表、组学分层摘要及质控信息。"""
        return await call("get_run_digest", lambda: run_digest(result_view(
            policy, experiment_id, run_id, delivery_id, trial_id,
        )))

    @server.tool(annotations=_READ)
    async def diagnose_run(experiment_id: str, run_id: str = "latest", delivery_id: str | None = None,
                           trial_id: str | None = None, limit: PageSize = 100) -> CallToolResult:
        """只读诊断固定的正式运行或隔离试验。"""
        return await call("diagnose_run", lambda: diagnose_view(
            result_view(policy, experiment_id, run_id, delivery_id, trial_id),
            limit=limit, project_alpha=scope.project().config.defaults.alpha,
        ))

    @server.tool(name="query_table", annotations=_READ)
    async def query(experiment_id: str, table: str, run_id: str = "latest", delivery_id: str | None = None,
                    omics: str | None = None, mode: str | None = None, filters: list[Filter] | None = None,
                    columns: list[str] | None = None, sort_by: str | None = None, ascending: bool = True,
                    limit: PageSize = 50, offset: Cursor = 0, trial_id: str | None = None) -> CallToolResult:
        """从结果 Parquet 读取有字节预算、过滤和投影的分页。"""
        def page():
            request = TableQuery(table=table, omics=omics, mode=mode, filters=filters or [], columns=columns,
                                 sort_by=sort_by, ascending=ascending, limit=limit, offset=offset)
            return query_table(result_view(policy, experiment_id, run_id, delivery_id, trial_id), request)

        return await call("query_table", page)

    @server.tool(name="trace_entity", annotations=_READ)
    async def trace(experiment_id: str, query: str, run_id: str = "latest", tissue: str | None = None,
                    contrast_id: str | None = None, limit: PageSize = 30,
                    trial_id: str | None = None) -> CallToolResult:
        """按确切实体 ID 追溯差异、反应及富集记录。"""
        return await call("trace_entity", lambda: trace_entity(
            result_view(policy, experiment_id, run_id, trial_id=trial_id), query, tissue=tissue,
            contrast=contrast_id, limit=limit,
        ))

    @server.tool(annotations=_COMPUTE)
    async def render_figure(experiment_id: str, kind: Literal["pathway_map", "submaps", "abstract"] = "pathway_map",
                            source_run_id: str | None = None, tissues: list[str] | None = None,
                            contrasts: list[str] | None = None, mode: Literal["atlas", "publication"] = "atlas",
                            style: str = "metro_atlas", color_mode: Literal["semantic", "pathway"] = "semantic",
                            strip_mode: Literal["omics", "axes"] = "omics", regions: list[str] | None = None,
                            submaps: list[str] | None = None,
                            formats: list[Literal["png", "pdf", "svg", "html", "tiff"]] | None = None,
                            return_preview_image: bool = False) -> CallToolResult:
        """由反应档案在本次专用输出目录生成图件。"""
        def render():
            results = policy.results(experiment_id)
            scope.experiment(scope.project(), experiment_id)
            source = identifier(source_run_id) if source_run_id else None
            output = scope.path(results / "_redraws" / identifier(experiment_id) / uuid4().hex)
            selected = formats or ["png"]
            if kind == "abstract":
                return render_abstract(root, experiment=experiment_id, source=source, output=output,
                                       tissues=tissues, contrasts=contrasts, style=style, color_mode=color_mode,
                                       formats=selected)
            overrides = {key: value for key, value in {
                "tissues": tissues, "contrasts": contrasts, "mode": mode, "style": style, "color_mode": color_mode,
                "strip_mode": strip_mode, "regions": regions, "formats": selected,
            }.items() if value is not None}
            overrides["abstract"] = {"enabled": False}
            overrides["submaps"] = {"enabled": kind == "submaps" or bool(submaps),
                                   "panels": submaps or [], "min_signal": "all" if submaps else "assessed"}
            return render_map(root, experiment_id, source, output, overrides)

        return await call("render_figure", render, preview_image=return_preview_image)

    @server.tool(annotations=_READ)
    async def diff_runs(experiment_id: str, run_a: str, run_b: str,
                        trial_a: str | None = None, trial_b: str | None = None) -> CallToolResult:
        """逐表对照两个已存在的结果包。"""
        def diff():
            from omicsweave.contract.diff import diff_bundles
            from omicsweave.contract.export import find_run

            a, b = (find_run(policy.results(experiment_id, trial), experiment_id, identifier(run))
                    for run, trial in ((run_a, trial_a), (run_b, trial_b)))
            return diff_payload(diff_bundles(a, b))

        return await call("diff_runs", diff)

    @server.tool(annotations=_COMPUTE)
    async def export_results(experiment_id: str, run_id: str, format: Literal["xlsx", "tsv"] = "xlsx") -> CallToolResult:
        """将指定结果包导出到本次专用目录。"""
        def export():
            from omicsweave.contract.export import export_bundle, find_run

            directory = policy.results(experiment_id)
            source = find_run(directory, experiment_id, identifier(run_id))
            target = scope.path(directory / "_exports" / experiment_id / uuid4().hex)
            target.mkdir(parents=True)
            written = export_bundle(source, format, target / "results.xlsx" if format == "xlsx" else target)
            return {"run_id": source.name, "output": str(written)}

        return await call("export_results", export)

    @server.tool(annotations=_READ)
    async def cleanup_deliveries(experiment_id: str, keep_last: Annotated[int, Field(ge=1)] = 5) -> CallToolResult:
        """预览保留策略命中项。实际删除使用本机 clean 命令。"""
        def cleanup():
            from omicsweave.core.retention import RetentionPolicy, plan_cleanup

            report = plan_cleanup(policy.results(experiment_id), RetentionPolicy(keep_last=keep_last), [experiment_id])
            return TypeAdapter(dict).dump_python(asdict(report), mode="json")

        return await call("cleanup_deliveries", cleanup)
