"""Verify a portable snapshot produced from an actual synthetic omics run."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from omicsweave import cli
from omicsweave.contract.dataset import ContractError, ValueType, read_dataset
from omicsweave.contract.exchange import EXCHANGE_FILE, AnalysisResultExchange, export_exchange, open_exchange


@pytest.fixture(scope="module")
def derived_run():
    from tests.metabolomics.test_derived_features import EXP, RunnerDerivedFeaturesTest

    fixture = RunnerDerivedFeaturesTest
    fixture.setUpClass()
    try:
        yield Path(fixture._tmp.name), EXP, fixture.met.parents[1]
    finally:
        fixture.tearDownClass()


def test_cli_exchanges_derived_data_without_loss(derived_run, capsys):
    project, experiment, run = derived_run
    output = project / "exchange-copy"

    def call(*argv):
        with pytest.raises(SystemExit) as result:
            cli.main(["--project-dir", str(project), "exchange", *argv, "--json"])
        payload = json.loads(capsys.readouterr().out)
        assert result.value.code == 0 and payload["ok"], payload
        return payload["data"]

    built = call("build", "--experiment", experiment, "--run", run.name, "--output", str(output))
    verified = call("verify", str(output))
    assert built["exchange"] == verified["exchange"]
    assert verified["verified_artifacts"] == len(built["exchange"]["artifacts"])
    assert not (output / "manifest.json").exists()
    assert not (output / "metadata" / "samples.parquet").exists()

    snapshot = open_exchange(output)
    derived, = [dataset for dataset in snapshot.manifest.datasets if dataset.metadata.value_type is ValueType.DERIVED]
    original = read_dataset(run / derived.directory)
    imported = snapshot.dataset(derived.directory)
    pd.testing.assert_frame_equal(original.matrix, imported.matrix)
    pd.testing.assert_frame_equal(original.features, imported.features)
    assert imported.meta == original.meta
    assert imported.meta.scale == "log2" and imported.meta.source
    assert not derived.metadata.source
    assert verified["exchange"]["run"]["bundle_schema_version"] == "0.1"


def test_snapshot_detects_changed_files_and_preserves_existing_outputs(derived_run, tmp_path):
    _, _, run = derived_run
    output = tmp_path / "export"
    manifest = export_exchange(run, output)
    original_manifest = (output / EXCHANGE_FILE).read_bytes()
    with pytest.raises(FileExistsError):
        export_exchange(run, output)
    assert (output / EXCHANGE_FILE).read_bytes() == original_manifest
    item = next(item for item in manifest.artifacts if item.kind == "matrix")
    with (output / item.path).open("ab") as data:
        data.write(b"changed synthetic payload")
    with pytest.raises(ContractError, match="校验失败"):
        open_exchange(output)


def test_exchange_schema_is_current():
    root = Path(__file__).resolve().parents[2]
    released = json.loads((root / "schemas/evidence-exchange/v1/exchange.schema.json").read_text(encoding="utf-8"))
    generated = AnalysisResultExchange.model_json_schema()
    assert released == generated
    assert "derived" in released["$defs"]["ValueType"]["enum"]
