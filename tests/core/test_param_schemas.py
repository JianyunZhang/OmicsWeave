"""同步检查参数自省、既有运行默认值和事务性配置编辑。"""
from __future__ import annotations

import csv
import pytest
import yaml

from omicsweave.core.param_schemas import schema
from omicsweave.core.path_scope import PathScope
from omicsweave.core.validate import validate_project
from omicsweave.service.project import EditOptions, QCUpdate, patch_config, update_sample_qc
from omicsweave.service.schemas import ErrorCode, ServiceError
from tests.mfa.test_mfa_runner import EXP, build_two_factor_project


def test_registered_defaults_match_runtime():
    from omicsweave.integration import pathway_map
    from omicsweave.integration.modes import flux_anchored, pathway_convergence, reaction_dossier, reversal
    from omicsweave.metabolomics import runner as metab
    from omicsweave.mfa.runner import FLUX_GROUP_SETTINGS
    from omicsweave.proteomics import runner as proteo
    from omicsweave.transcriptomics.runner import DEFAULT_GENE_SETS

    for target, runtime in [
        ("omics/metabolomics", metab._DEFAULTS), ("omics/proteomics", proteo._DEFAULTS),
        ("omics/mfa", {"flux_groups": FLUX_GROUP_SETTINGS}),
        ("omics/transcriptomics", {"gene_sets": DEFAULT_GENE_SETS}),
        ("integration/flux_anchored", flux_anchored.DEFAULTS), ("integration/reversal", reversal.DEFAULTS),
        ("integration/reaction_dossier", reaction_dossier.DEFAULTS),
        ("integration/pathway_convergence", pathway_convergence.DEFAULTS), ("pathway_map", pathway_map.DEFAULTS),
    ]:
        properties = schema(target)["properties"]
        assert {key: properties[key]["default"] for key in runtime} == runtime
        assert all("x-tunable-class" in prop for prop in properties.values())


def test_typo_is_reported_before_execution(tmp_path):
    build_two_factor_project(tmp_path)
    cfg = tmp_path / "experiments" / EXP / "experiment.yaml"
    raw = yaml.safe_load(cfg.read_text(encoding="utf-8"))
    raw["omics"]["mfa"].setdefault("analysis", {})["flux_groups"] = {"enabledd": False}
    cfg.write_text(yaml.safe_dump(raw), encoding="utf-8")
    report = validate_project(tmp_path)
    assert any(issue.where == "omics.mfa.analysis.flux_groups.enabledd" and issue.code == "CONFIG_SCHEMA_ERROR"
               for issue in report.errors)


def test_patch_atomic_compare_and_swap_and_qc(tmp_path):
    build_two_factor_project(tmp_path)
    scope = PathScope.create(tmp_path)
    experiment = tmp_path / "experiments" / EXP
    path = experiment / "experiment.yaml"
    initial = path.read_bytes()
    patch = {"experiment": {"description": "合成项目验证"}}
    preview = patch_config(scope, EXP, patch, EditOptions())
    assert not preview["applied"] and preview["validation"]["ok"] and path.read_bytes() == initial
    path.write_bytes(initial + b"\n# external edit\n")
    with pytest.raises(ServiceError) as stale:
        patch_config(scope, EXP, patch, EditOptions(dry_run=False, expected_sha256=preview["sha256"], reason="test"))
    assert stale.value.detail.code is ErrorCode.JOB_CONFLICT
    assert path.read_bytes().endswith(b"# external edit\n")
    preview = patch_config(scope, EXP, patch, EditOptions())
    applied = patch_config(scope, EXP, patch,
                           EditOptions(dry_run=False, expected_sha256=preview["sha256"], reason="test"))
    assert applied["applied"] and yaml.safe_load(path.read_text(encoding="utf-8"))["experiment"]["description"] == "合成项目验证"

    sample_path = experiment / "samples.tsv"
    old = sample_path.read_bytes()
    with sample_path.open(encoding="utf-8", newline="") as handle:
        sample = next(csv.DictReader(handle, delimiter="\t"))["sample_id"]
    entries = [QCUpdate(sample_id=sample, qc_flag="flag", qc_reason="合成异常")]
    preview = update_sample_qc(scope, EXP, entries, EditOptions())
    assert preview["recommended_rerun_from"] == "run" and sample_path.read_bytes() == old
    applied = update_sample_qc(scope, EXP, entries,
                               EditOptions(dry_run=False, expected_sha256=preview["sha256"], reason="test"))
    assert applied["applied"] and "合成异常" in sample_path.read_text(encoding="utf-8")
