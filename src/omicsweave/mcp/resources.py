"""MCP 中可寻址的只读项目与运行元信息。"""
from __future__ import annotations

import json

from mcp.server.fastmcp import FastMCP

from omicsweave.core.param_schemas import schema
from omicsweave.core.path_scope import identifier
from omicsweave.service.execution import job_operation
from omicsweave.service.project import experiment_view, project_overview
from omicsweave.service.query import run_digest

from .security import ServerPolicy
from .tools import result_view


def register_resources(server: FastMCP, policy: ServerPolicy) -> None:
    @server.resource("omicsweave://project/overview", mime_type="application/json")
    def overview() -> dict:
        return project_overview(policy.scope)

    @server.resource("omicsweave://schema/{target}", mime_type="application/json")
    def root_schema(target: str) -> dict:
        return schema(target)

    @server.resource("omicsweave://schema/{family}/{name}", mime_type="application/json")
    def parameter_schema(family: str, name: str) -> dict:
        return schema(f"{family}/{name}")

    @server.resource("omicsweave://experiments/{exp_id}/config", mime_type="application/json")
    def config(exp_id: str) -> dict:
        return experiment_view(policy.scope, exp_id)

    @server.resource("omicsweave://experiments/{exp_id}/latest", mime_type="application/json")
    def latest(exp_id: str) -> dict:
        from omicsweave.core.delivery import build_latest_payload, resolve_latest_delivery

        record, source = resolve_latest_delivery(policy.results(exp_id), exp_id)
        if record is None:
            raise FileNotFoundError(f"实验 {exp_id} 没有已完成的交付")
        index = policy.scope.path(record.index_path)
        saved = json.loads(index.read_text(encoding="utf-8"))
        return {**build_latest_payload(record), "source": source, "stages": saved.get("stages", {})}

    @server.resource("omicsweave://runs/{exp_id}/{run_id}/manifest", mime_type="application/json")
    def manifest(exp_id: str, run_id: str) -> dict:
        from omicsweave.contract.export import find_run
        from omicsweave.contract.manifest import read_manifest

        path = find_run(policy.results(exp_id), exp_id, identifier(run_id))
        return read_manifest(path).model_dump(mode="json")

    @server.resource("omicsweave://runs/{exp_id}/{run_id}/digest", mime_type="application/json")
    def digest(exp_id: str, run_id: str) -> dict:
        return run_digest(result_view(policy, exp_id, run_id))

    @server.resource("omicsweave://jobs/{job_id}/state", mime_type="application/json")
    def job(job_id: str) -> dict:
        return job_operation(policy.scope.root, policy.job_root(), "status", job=job_id)
