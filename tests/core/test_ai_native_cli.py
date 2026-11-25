"""A1 从真实命令行入口验证机器协议、状态与结果读取。"""
from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
import yaml

from omicsweave import cli
from omicsweave.contract.manifest import RunManifest, file_record, write_manifest
from omicsweave.core.jobs import read_status
from omicsweave.service.query import ResultView
from tests.mfa.test_mfa_runner import EXP, build_two_factor_project


def call(capsys, root, *args, expected=0):
    with pytest.raises(SystemExit) as exit_info:
        cli.main(["--project-dir", str(root), *args])
    output = capsys.readouterr().out
    assert exit_info.value.code == expected, output
    return json.loads(output)


def wait_status(directory, predicate, timeout=90):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        state = read_status(directory)
        if predicate(state):
            return state
        time.sleep(0.1)
    pytest.fail(f"作业状态超时 {state}")


def _copy_fixture(root):
    build_two_factor_project(root)
    spec = root / "experiments" / EXP / "experiment.yaml"
    doc = yaml.safe_load(spec.read_text(encoding="utf-8"))
    doc["integration"] = {"enabled": True, "modes": ["reversal"]}
    spec.write_text(yaml.safe_dump(doc), encoding="utf-8")
    (root / "analysis").mkdir()
    project = root / "project.yaml"
    doc = yaml.safe_load(project.read_text(encoding="utf-8"))
    doc["pipeline"] = {"scripts": [{"command": "analysis/work.py"}]}
    project.write_text(yaml.safe_dump(doc), encoding="utf-8")
    return root / "analysis" / "work.py"


def test_json_parse_validation_and_backward_compatibility(tmp_path, capsys):
    payload = call(capsys, tmp_path, "run", "--json", expected=2)
    assert payload["error"]["code"] == "CLI_USAGE_ERROR"
    (tmp_path / "project.yaml").write_text("project: [", encoding="utf-8")
    payload = call(capsys, tmp_path, "validate", "--format", "json", expected=10)
    assert payload["error"]["code"] == "CONFIG_SCHEMA_ERROR"
    assert payload["error"]["details"]["issues"]
    payload = call(capsys, tmp_path, "abstract", "--format", "json", "--format", "png", expected=2)
    assert payload["error"]["code"] == "CLI_USAGE_ERROR"


def test_machine_pipeline_preserves_subprocess_streams(tmp_path):
    script = _copy_fixture(tmp_path)
    script.write_text("import os\nos.write(1, b'native child stdout\\n')\nprint('python child stdout')", encoding="utf-8")
    cli_proc = subprocess.run(
        [sys.executable, "-m", "omicsweave.cli", "--project-dir", str(tmp_path), "pipeline", "--json"],
        capture_output=True, timeout=90,
    )
    assert cli_proc.returncode == 0, cli_proc.stderr.decode(errors="replace")
    response = json.loads(cli_proc.stdout)
    assert response["ok"]
    assert b"native child stdout" in cli_proc.stderr and b"python child stdout" in cli_proc.stderr
    assert len(response["data"]["deliveries"]) == 1


def test_background_status_retention_cancel_and_crash(tmp_path, capsys):
    script = _copy_fixture(tmp_path)
    script.write_text("import time\nprint('checkpoint', flush=True)\ntime.sleep(7)", encoding="utf-8")
    submitted = call(capsys, tmp_path, "pipeline", "--background", "--json")["data"]
    directory, job_id = Path(submitted["state_path"]).parent, submitted["job_id"]
    assert submitted["status"] == "queued"
    state = wait_status(directory, lambda s: s["progress"] and s["progress"]["stage"].startswith("script_"))
    previous_heartbeat = state["heartbeat_at"]
    wait_status(directory, lambda s: s["heartbeat_at"] != previous_heartbeat)
    cleanup = call(capsys, tmp_path, "clean", "--drop-failed", "--stale-running-hours", "0", "--json")["data"]
    assert cleanup["to_remove"] == []
    assert any(entry["reason"] == "active_job" for entry in cleanup["kept"])
    called = call(capsys, tmp_path, "cancel", job_id, "--json")["data"]
    assert called["cancel_requested"]
    state = wait_status(directory, lambda s: s["status"] == "cancelled")
    assert state["exit_code"] == 30
    assert not (tmp_path / "results" / "latest.json").exists()
    logs = call(capsys, tmp_path, "logs", job_id, "--tail", "2", "--since-event", "1", "--json")["data"]
    assert [event["seq"] for event in logs["events"]] == [2, 3]
    assert logs["has_more_events"]
    assert "checkpoint" in "\n".join(logs["lines"])

    second = call(capsys, tmp_path, "run", "--experiment", EXP, "--steps", "qc", "--background", "--json")["data"]
    other_dir = Path(second["state_path"]).parent
    running = wait_status(other_dir, lambda s: s["status"] == "running")
    assert running["pid"] == second["pid"] and running["process_identity"]
    os.kill(running["pid"], signal.SIGTERM)
    wait_status(other_dir, lambda s: s["status"] == "stalled", timeout=15)
    assert call(capsys, tmp_path, "status", "--job", second["job_id"], "--json")["data"]["status"] == "stalled"


def test_result_pagination_and_pinned_legacy_latest(tmp_path, capsys):
    root = tmp_path / "results"
    (tmp_path / "project.yaml").write_text(
        "schema_version: 1\nproject: {id: demo, species: mus_musculus}\n", encoding="utf-8",
    )
    run = root / "study" / "run_one"
    layer = run / "omics" / "mfa" / "tables"
    layer.mkdir(parents=True)
    file = layer / "differential.parquet"
    pq.write_table(pa.table({
        "feature_id": ["id", "id-other", "id"], "tissue": ["liver"] * 3,
        "omics": ["mfa"] * 3, "contrast_id": ["a"] * 3, "pvalue": [0.1, float("inf"), 0.3],
        "selected": [True, None, False],
    }), file, row_group_size=1)
    empty_enrichment = layer / "enrichment.parquet"
    pq.write_table(pa.table({"hit_features": pa.nulls(0)}), empty_enrichment)
    manifest = RunManifest(run_id=run.name, project_id="demo", experiment_id="study",
                           files=[file_record(path, run, "output") for path in (file, empty_enrichment)])
    write_manifest(manifest, run)
    (run / "index.json").write_text(json.dumps({
        "status": "completed", "experiment": "study", "bundles": {"run": {"run_id": run.name}},
    }), encoding="utf-8")

    direct = call(capsys, tmp_path, "latest", "--json")
    assert direct["delivery_id"] == run.name and "ok" not in direct
    wrapped = call(capsys, tmp_path, "latest", "--format", "json")
    assert wrapped["data"]["deliveries"][0]["delivery_id"] == run.name
    digest = call(capsys, tmp_path, "inspect", "--experiment", "study", "--json")["data"]
    assert digest["omics_summary"][0]["n_tested"] == 2
    assert digest["tables"][0]["rows"] == 3

    query = ["query", "--experiment", "study", "--table", "differential", "--json"]
    page = call(capsys, tmp_path, *query, "--sort", "pvalue", "--descending", "--limit", "1")["data"]
    assert page["next_offset"] == 1
    assert page["rows"][0]["pvalue"] is None
    page = call(capsys, tmp_path, *query, "--filters", '[{"column":"feature_id","op":"eq","value":"id"}]',
                "--columns", "pvalue", "--offset", "1", "--limit", "1")["data"]
    assert page["rows"] == [{"pvalue": 0.3}] and page["next_offset"] is None
    failure = call(capsys, tmp_path, *query, "--filters", '[{"column":"nope","op":"eq","value":"id"}]', expected=2)
    assert failure["error"]["code"] == "CLI_USAGE_ERROR"
    traced = call(capsys, tmp_path, "trace", "--experiment", "study", "--entity", "id", "--json")["data"]
    assert len(traced["records"]) == 2

    diff = call(capsys, tmp_path, "diff", str(run), str(run), "--json")
    assert diff["data"]["identical"]
    with pytest.raises(Exception, match="无效的结果标识符"):
        ResultView(root, "../study", "latest")


def test_query_byte_budget_paginates_oversized_page(tmp_path):
    from omicsweave.service.query import TableQuery, query_table

    run = tmp_path / "study" / "run_large"
    run.mkdir(parents=True)
    file = run / "large.parquet"
    pq.write_table(pa.table({"long": ["a" * 80_000, "b" * 80_000]}), file)
    write_manifest(RunManifest(run_id=run.name, experiment_id="study", project_id="demo",
                               files=[file_record(file, run, "output")]), run)
    view = ResultView(tmp_path, "study", run.name)
    page = query_table(view, TableQuery(table="large", limit=2))
    assert page["truncated"] and page["next_offset"] == 1
    assert len(page["rows"]) == 1
    next_page = query_table(view, TableQuery(table="large", limit=2, offset=page["next_offset"]))
    assert not next_page["truncated"] and next_page["next_offset"] is None
