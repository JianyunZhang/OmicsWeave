"""MCP 的内存协议与真实 stdio 验收。"""
from __future__ import annotations

import json
import os
import shutil
import sys
from datetime import timedelta
from pathlib import Path

import pytest
import yaml

pytest.importorskip("mcp")

import anyio
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
from mcp.shared.memory import create_connected_server_and_client_session
from pydantic import AnyUrl

from omicsweave.mcp.server import build_server

DEMO = Path(__file__).resolve().parents[2] / "examples" / "demo_project"


async def call(session, name, data):
    result = await session.call_tool(name, data, read_timeout_seconds=timedelta(seconds=30))
    assert not result.isError, result.content
    assert result.structuredContent["ok"]
    return result.structuredContent["data"]


def test_mcp_in_memory_demo(tmp_path):
    project = tmp_path / "demo"
    shutil.copytree(DEMO, project, ignore=shutil.ignore_patterns("results", "analysis_results", "__pycache__"))
    path = project / "project.yaml"
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    raw["pipeline"] = {"scripts": [{"command": "analysis/marker.py"}]}
    path.write_text(yaml.safe_dump(raw), encoding="utf-8")
    (project / "analysis").mkdir(exist_ok=True)
    (project / "analysis" / "marker.py").write_text(
        "from pathlib import Path\nPath('analysis/ran').touch()\n", encoding="utf-8",
    )

    async def check():
        server = build_server(project)
        async with create_connected_server_and_client_session(server, read_timeout_seconds=timedelta(seconds=35)) as client:
            names = {tool.name for tool in (await client.list_tools()).tools}
            assert {"validate_experiment", "submit_pipeline", "get_run_digest", "render_figure", "diff_runs"} <= names
            assert len((await client.list_prompts()).prompts) == 4
            response = await client.read_resource(AnyUrl("omicsweave://schema/integration/reversal"))
            assert json.loads(response.contents[0].text)["properties"]["alpha"]
            assert "demo_glp1" in (await call(client, "list_experiments", {}))["experiments"]
            common = {"experiment_id": "demo_glp1"}
            assert (await call(client, "validate_experiment", common))["ok"]
            assert (await call(client, "submit_pipeline", common))["dry_run"]
            assert not (project / "results" / "_jobs").exists()
            submitted = await call(client, "submit_pipeline", {**common, "dry_run": False, "wait_seconds": 0.2})
            seen_progress = 0
            with anyio.fail_after(120):
                while True:
                    state = await call(client, "get_job_status", {"job_id": submitted["job_id"], "since_seq": seen_progress})
                    if state["status"] not in {"queued", "running", "unresponsive"}:
                        break
                    seen_progress = state["log"]["next_event"]
                    await anyio.sleep(0.3)
            assert state["status"] == "completed", state
            assert not (project / "analysis" / "ran").exists()
            digest = await call(client, "get_run_digest", common)
            page = await call(client, "query_table", {**common, "table": "differential", "omics": "transcriptomics", "limit": 10})
            assert len(page["rows"]) == 10 and digest["tables"]
            same = await call(client, "diff_runs", {**common, "run_a": digest["run_id"], "run_b": digest["run_id"]})
            assert same["identical"]
            pointer = project / "results" / "latest.json"
            original = pointer.read_bytes()
            trial_args = {**common, "start": "integrate", "source_run_id": digest["run_id"],
                          "overrides": {"integration.parameters.reversal.equivalence_margin": 0.75},
                          "rationale": "synthetic comparison"}
            preview = await call(client, "submit_pipeline", trial_args)
            assert preview["dry_run"] and not Path(preview["trial"]["results_dir"]).exists()
            trial_job = await call(client, "submit_pipeline", {**trial_args, "dry_run": False, "wait_seconds": 0})
            with anyio.fail_after(80):
                while trial_job["status"] in {"queued", "running", "unresponsive"}:
                    trial_job = await call(client, "get_job_status", {"job_id": trial_job["job_id"]})
                    await anyio.sleep(0.3)
            assert trial_job["status"] == "completed", trial_job
            trial_id = trial_job["result"]["trial"]["trial_id"]
            diagnosed = await call(client, "diagnose_run", {**common, "trial_id": trial_id})
            assert diagnosed["trial_id"] == trial_id and diagnosed["evidence_availability"]
            comparison = await call(client, "diff_runs", {
                **common, "run_a": digest["run_id"], "run_b": diagnosed["run_id"], "trial_b": trial_id,
            })
            assert comparison["config_diff"] and comparison["lineage"]["b"]["parent_run_id"] == digest["run_id"]
            assert pointer.read_bytes() == original

    anyio.run(check)


def test_mcp_stdio_handshake():
    async def check():
        params = StdioServerParameters(
            command=sys.executable, args=["-m", "omicsweave", "--project-dir", str(DEMO), "mcp", "serve"],
            env={k: v for k, v in os.environ.items() if k in {"PYTHONPATH", "PYTHONUTF8", "PYTHONDONTWRITEBYTECODE"}},
        )
        with anyio.fail_after(60):
            async with stdio_client(params) as (read, write):
                async with ClientSession(read, write, read_timeout_seconds=timedelta(seconds=35)) as client:
                    initialized = await client.initialize()
                    assert initialized.serverInfo.name == "OmicsWeave"
                    project = await call(client, "list_experiments", {})
                    assert "demo_glp1" in project["experiments"]

    anyio.run(check)


def test_mcp_renders_and_returns_one_png():
    from tests.integration.test_reaction_dossier import DIET, DossierBundleTest
    from tests.mfa.test_mfa_runner import EXP

    fixture = DossierBundleTest()
    fixture.setUp()
    try:
        _, bundle, _ = fixture.dossier(pathway_map={"enabled": False})

        async def check():
            server = build_server(fixture.root)
            async with create_connected_server_and_client_session(server, read_timeout_seconds=timedelta(seconds=60)) as client:
                payload = await client.call_tool("render_figure", {
                    "experiment_id": EXP, "source_run_id": bundle.name, "kind": "abstract",
                    "tissues": ["liver"], "contrasts": [DIET], "formats": ["png"], "return_preview_image": True,
                })
                assert not payload.isError, payload.content
                assert len(payload.structuredContent["data"]["files"]) == 1
                assert any(part.type == "image" and part.mimeType == "image/png" for part in payload.content)

        anyio.run(check)
    finally:
        fixture.tearDown()
