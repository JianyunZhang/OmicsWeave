"""MCP 的项目边界、间接路径和写入授权回归。"""
from __future__ import annotations

from datetime import timedelta

import pytest
import yaml

pytest.importorskip("mcp")

import anyio
from mcp.shared.memory import create_connected_server_and_client_session

from omicsweave.mcp.server import build_server
from tests.mfa.test_mfa_runner import EXP, build_two_factor_project


def test_mcp_checks_write_gate_and_indirect_input_paths(tmp_path):
    project = tmp_path / "project"
    project.mkdir()
    build_two_factor_project(project)
    path = project / "experiments" / EXP / "experiment.yaml"
    before = path.read_bytes()
    outside = tmp_path / "outside.tsv"
    outside.write_text("synthetic\n", encoding="utf-8")
    patch = {"experiment": {"description": "test"}}

    async def check():
        server = build_server(project)
        async with create_connected_server_and_client_session(server, read_timeout_seconds=timedelta(seconds=30)) as client:
            preview = await client.call_tool("patch_experiment_config", {"experiment_id": EXP, "patch": patch})
            assert not preview.isError and not preview.structuredContent["data"]["applied"]
            denied = await client.call_tool("patch_experiment_config", {
                "experiment_id": EXP, "patch": patch, "dry_run": False, "reason": "test",
                "expected_sha256": preview.structuredContent["data"]["sha256"],
            })
            assert denied.isError and denied.structuredContent["error"]["code"] == "SECURITY_POLICY_VIOLATION"
            for name, args in [
                ("inspect_input_files", {"paths": [str(outside)]}),
                ("inspect_input_files", {"paths": [r"\\example.invalid\share\synthetic.tsv"]}),
                ("get_run_digest", {"experiment_id": "../other"}),
                ("patch_experiment_config", {"experiment_id": EXP, "patch": {"design": {"subjects": str(outside)}}}),
            ]:
                response = await client.call_tool(name, args)
                assert response.isError and response.structuredContent["error"]["code"] == "SECURITY_POLICY_VIOLATION"
            assert path.read_bytes() == before

            config = yaml.safe_load(before)
            config["design"]["subjects"] = str(outside)
            path.write_text(yaml.safe_dump(config), encoding="utf-8")
            response = await client.call_tool("submit_pipeline", {"experiment_id": EXP, "dry_run": False})
            assert response.isError and response.structuredContent["error"]["code"] == "SECURITY_POLICY_VIOLATION"
            assert not (project / "results" / "_jobs").exists()

    anyio.run(check)
