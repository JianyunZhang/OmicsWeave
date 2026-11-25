"""以明确的合成统计与校准记录检验诊断分类。"""
from __future__ import annotations

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
import yaml

from omicsweave.core.diagnostics import inspect_evidence
from omicsweave.core.path_scope import PathScope
from omicsweave.service.project import EditOptions, patch_config
from omicsweave.service.schemas import ServiceError
from tests.mfa.test_mfa_runner import EXP, build_two_factor_project


def test_diagnostics_classifies_measurements_and_missing_calibration(tmp_path):
    from omicsweave.contract.manifest import RunManifest, file_record, write_manifest
    from omicsweave.service.query import ResultView, diagnose_view

    run = tmp_path / "experiment" / "run_diagnostics"
    run.mkdir(parents=True)
    tables = {}

    def write(key, rows):
        path = run / key
        pq.write_table(pa.Table.from_pylist(rows), path)
        tables[key] = path

    subgroup = {"tissue": "liver", "contrast_id": "diet"}
    write("differential.parquet", [{**subgroup, "omics": "mfa", "min_achievable_p": 0.1,
                                   "n_numerator": 3, "n_denominator": 3}])
    write("contrast_reliability.parquet", [{**subgroup, "status": "qc_confounded", "n_usable_numerator": 1}])
    write("sample_quality.parquet", [{"flag": "flag", "tissue": "liver", "metric": "scale"}])
    write("assessment_calibration.parquet", [
        {**subgroup, "method": "decoy", "metric": "L1", "target": 10, "decoy_expected": 9.5,
         "decoy_sd": 1.0, "n_replicates": 20},
        {**subgroup, "tissue": "heart", "method": "decoy", "metric": "L1", "target": 10, "decoy_expected": 10.0,
         "decoy_sd": None, "n_replicates": 1},
    ])
    findings, counts = inspect_evidence(tables, ["输入编码需要复核"], 0.05)
    categories = {finding.code: finding.category for finding in findings}
    assert categories == {
        "L1_CLOSE_TO_DECOYS": "statistical_guardrail", "QC_CONFOUNDED": "statistical_guardrail",
        "RANK_TEST_POWER_LIMIT": "statistical_guardrail", "SAMPLE_QC_FLAG": "sensitivity_explorable",
        "VALIDATION_WARNING": "engineering_fixable",
    }
    assert counts["decoy_strata"] == 1 and counts["decoy_unassessed_strata"] == 1
    write_manifest(RunManifest(run_id=run.name, project_id="demo", experiment_id="experiment",
                               files=[file_record(path, run, "output") for path in tables.values()]), run)
    report = diagnose_view(ResultView(tmp_path, "experiment", run.name))
    assert report["optimization_blocked"] is True and report["blocking_findings"] == 3
    assert report["decoy_calibration_status"] == "partially_assessed"


def test_result_in_trial_requires_post_hoc_panel_with_existing_basis(tmp_path):
    from omicsweave.contract.manifest import RunManifest, write_manifest
    from omicsweave.core.config import load_project

    build_two_factor_project(tmp_path)
    exp = tmp_path / "experiments" / EXP
    exp_yaml = exp / "experiment.yaml"
    document = yaml.safe_load(exp_yaml.read_text(encoding="utf-8"))
    document["integration"] = {"parameters": {"reaction_dossier": {"panels": "panels.yaml"}}}
    exp_yaml.write_text(yaml.safe_dump(document), encoding="utf-8")
    path = exp / "panels.yaml"
    path.write_text("panels: {test: {segments: [{reactions: [PDH_m]}]}}\n", encoding="utf-8")
    before = path.read_bytes()
    seeded = tmp_path / "results" / "_trials" / ("trial_" + "a" * 32) / EXP / "run_observed"
    seeded.mkdir(parents=True)
    write_manifest(RunManifest(run_id=seeded.name, project_id=load_project(tmp_path).id, experiment_id=EXP), seeded)
    scope = PathScope.create(tmp_path)
    patch = {"panels": {"test": {"title": "changed"}}}
    preview = patch_config(scope, EXP, patch, EditOptions(), target="panels_yaml")
    assert preview["guardrail_notices"]
    options = EditOptions(dry_run=False, reason="read result", expected_sha256=preview["sha256"])
    with pytest.raises(ServiceError, match="未通过"):
        patch_config(scope, EXP, patch, options, target="panels_yaml")
    assert path.read_bytes() == before
    patch["selection"] = {"mode": "post_hoc", "basis": "参考 run_observed"}
    applied = patch_config(scope, EXP, patch, options, target="panels_yaml")
    assert applied["applied"]
    assert yaml.safe_load(path.read_bytes())["selection"]["mode"] == "post_hoc"
