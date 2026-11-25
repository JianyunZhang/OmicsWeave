"""真实 CLI 入口的隔离试验、清单兼容和回退选包。"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest
import yaml

from omicsweave import cli
from omicsweave.contract.manifest import read_manifest
from omicsweave.core import delivery
from omicsweave.core.overrides import LINEAGE_PATH, OverrideError, apply_overrides, parse_overrides
from omicsweave.core.retention import RetentionPolicy, plan_cleanup
from tests.mfa.test_mfa_runner import EXP, build_two_factor_project


def _call(capsys, root, *tokens, expected=0):
    with pytest.raises(SystemExit) as outcome:
        cli.main(["--project-dir", str(root), *tokens, "--json"])
    response = json.loads(capsys.readouterr().out)
    assert outcome.value.code == expected, response
    return response


def test_parse_file_and_inline_overrides(tmp_path):
    source = tmp_path / "override.yaml"
    source.write_text("omics:\n  mfa:\n    analysis:\n      flux_groups:\n        enabled: false\n", encoding="utf-8")
    raw = {"integration": {"parameters": {"reaction_dossier": {"pathway_map": {"submaps": False}}}}}
    overrides = parse_overrides(
        ["integration.parameters.reaction_dossier.pathway_map.submaps.min_signal=all"], source,
    )
    result = apply_overrides(raw, overrides)
    submaps = result["integration"]["parameters"]["reaction_dossier"]["pathway_map"]["submaps"]
    assert submaps == {"enabled": False, "min_signal": "all"}
    assert result["omics"]["mfa"]["analysis"]["flux_groups"]["enabled"] is False
    assert raw["integration"]["parameters"]["reaction_dossier"]["pathway_map"]["submaps"] is False
    with pytest.raises(OverrideError, match="重叠"):
        parse_overrides(["omics.mfa.analysis.flux_groups={}"], source)
    with pytest.raises(OverrideError, match="不可"):
        apply_overrides(raw, {"integration.parameters.reaction_dossier.panels": "elsewhere.yaml"})
    source.write_text("design.adaptive_shrinkage: true\ndesign.adaptive_shrinkage: false\n", encoding="utf-8")
    with pytest.raises(OverrideError, match="重复键"):
        parse_overrides([], source)


def test_reaction_grade_changes_align_by_reaction_tissue_and_contrast():
    import pandas as pd

    from omicsweave.contract.diff import compare_frames

    first = pd.DataFrame([{"flux_id": "PC_m", "tissue": "liver", "contrast_id": "diet", "flux_claim": "L1"},
                          {"flux_id": "PDH_m", "tissue": "liver", "contrast_id": "diet", "flux_claim": "L3"}])
    second = first.iloc[::-1].copy()
    second.loc[second.flux_id == "PC_m", "flux_claim"] = "L0"
    path = "integration/reaction_dossier/tables/reaction_assessment.parquet"
    _, _, transitions = compare_frames(path, first, second)
    assert transitions == [{"path": path, "column": "flux_claim", "from": "L1", "to": "L0", "n": 1}]


def test_trial_pipeline_keeps_formal_results_and_links_source(tmp_path, capsys, monkeypatch):
    build_two_factor_project(tmp_path)
    source = tmp_path / "experiments" / EXP / "experiment.yaml"
    document = yaml.safe_load(source.read_text(encoding="utf-8"))
    document["integration"] = {"enabled": True, "modes": ["reversal"]}
    source.write_text(yaml.safe_dump(document), encoding="utf-8")
    formal = _call(capsys, tmp_path, "pipeline", "--experiment", EXP)["data"]["deliveries"][0]
    base = Path(formal["run_dir"])
    pointer = Path(formal["latest_pointer"])
    expected = {p: p.read_bytes() for p in (source, pointer, base / "manifest.json", base / "index.json")}
    trial_args = (
        "pipeline", "--experiment", EXP, "--source", base.name, "--rationale", "合成敏感性对照",
        "--override", "integration.parameters.reversal.equivalence_margin=0.75",
    )
    preview = _call(capsys, tmp_path, *trial_args, "--dry-run")["data"]
    assert preview["dry_run"] and not Path(preview["trial"]["results_dir"]).exists()

    completed = _call(capsys, tmp_path, *trial_args)["data"]
    isolated = Path(completed["deliveries"][0]["integration_dir"])
    trial = completed["trial"]["trial_id"]
    assert isolated.parents[2] == pointer.parent / "_trials"
    lineage = json.loads((isolated / LINEAGE_PATH).read_text(encoding="utf-8"))
    assert lineage["parent_run_id"] == base.name and lineage["source_run_id"] == base.name
    manifest = read_manifest(isolated)
    entry = next(item for item in manifest.files if item.path == LINEAGE_PATH.as_posix())
    assert entry.sha256 == hashlib.sha256((isolated / LINEAGE_PATH).read_bytes()).hexdigest()
    assert manifest.bundle_schema_version == "0.1"
    assert all(path.read_bytes() == data for path, data in expected.items())

    digest = _call(capsys, tmp_path, "inspect", "--experiment", EXP, "--trial", trial)["data"]
    assert digest["run_id"] == isolated.name and digest["lineage"]["parent_run_id"] == base.name
    rows = _call(capsys, tmp_path, "query", "--experiment", EXP, "--trial", trial, "--table", "differential",
                 "--omics", "mfa", "--limit", "2")["data"]["rows"]
    assert rows
    diff = _call(capsys, tmp_path, "diff", base.name, isolated.name, "--experiment", EXP,
                 "--trial-b", trial, expected=1)["data"]
    assert diff["config_diff"] and diff["lineage"]["b"]["trial_id"] == trial

    incompatible = _call(capsys, tmp_path, *trial_args,
                         "--override", "omics.mfa.analysis.flux_groups.enabled=false",
                         "--dry-run", expected=13)
    assert incompatible["error"]["details"]["stale_bundles"][0]["recommended_start"] == "run"
    monkeypatch.setattr(delivery, "read_latest_pointer", lambda *args: None)
    record, mode = delivery.resolve_latest_delivery(pointer.parent, EXP)
    assert mode == "scan" and record.delivery_id == base.name
    plan = plan_cleanup(pointer.parent, RetentionPolicy(keep_last=1, prune_orphan_bundles=True), [EXP])
    assert all(item.path != isolated for item in plan.to_remove)
    assert any(item.path == base and item.reason == "trial_reference" for item in plan.kept)
    table = base / "omics" / "mfa" / "tables" / "differential.parquet"
    table.write_bytes(table.read_bytes() + b"changed test artifact")
    refusal = _call(capsys, tmp_path, *trial_args, "--dry-run", expected=13)
    assert "omics/mfa/tables/differential.parquet" in refusal["error"]["details"]["changed_files"]
