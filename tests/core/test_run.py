#!/usr/bin/env python3
"""omicsweave run 的运行框架：校验不通过不开始运行、组学与步骤的选择、结果包与运行清单、命令行入口。"""

from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd
import yaml

from omicsweave import cli
from omicsweave.contract.manifest import read_manifest
from omicsweave.core import run as run_module
from omicsweave.core.run import RunError, run_experiment
from tests.mfa.test_mfa_runner import EXP, build_two_factor_project


def _quiet(func, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()):
        return func(*args, **kwargs)


class RunFrameworkTest(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.root = Path(self._tmp.name).resolve()
        build_two_factor_project(self.root)
        self.exp_dir = self.root / "experiments" / EXP

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def edit_experiment(self, change) -> None:
        path = self.exp_dir / "experiment.yaml"
        doc = yaml.safe_load(path.read_text(encoding="utf-8"))
        change(doc)
        path.write_text(yaml.safe_dump(doc, sort_keys=False), encoding="utf-8")

    def enable_metabolomics(self) -> None:
        """再启用一个尚未接入统一运行的组学（表格层输入，样本齐全，可以通过校验）。"""
        peaks = self.exp_dir / "peaks.tsv"
        peaks.write_text("feature\n", encoding="utf-8")
        self.edit_experiment(lambda doc: doc["omics"].update(metabolomics={
            "enabled": True,
            "input": {"level": "table", "format": "peak_table", "files": [f"experiments/{EXP}/peaks.tsv"]}}))
        samples = pd.read_csv(self.exp_dir / "samples.tsv", sep="\t", dtype=str)
        extra = samples.assign(sample_id="met_" + samples["sample_id"], omics="metabolomics", tissue="liver")
        pd.concat([samples, extra]).to_csv(self.exp_dir / "samples.tsv", sep="\t", index=False)

    def test_validation_errors_stop_before_anything_is_written(self) -> None:
        self.edit_experiment(lambda doc: doc["omics"]["mfa"]["input"].pop("model_id"))
        with self.assertRaisesRegex(RunError, "model_id"):
            run_experiment(EXP, root=self.root)
        self.assertFalse((self.root / "analysis_results").exists())

    def test_omics_selection(self) -> None:
        with self.assertRaisesRegex(RunError, "未在 experiment.yaml 中启用"):
            run_experiment(EXP, ["transcriptomics"], root=self.root)
        with self.assertRaisesRegex(RunError, "未知组学"):
            run_experiment(EXP, ["lipidomics"], root=self.root)
        with mock.patch.dict(run_module.RUNNERS, clear=True):
            with self.assertRaisesRegex(RunError, "尚未接入统一运行"):
                run_experiment(EXP, ["mfa"], root=self.root)
            with self.assertRaisesRegex(RunError, "没有可统一运行的组学"):
                run_experiment(EXP, root=self.root)

    def test_missing_optional_dependency_names_the_extra(self) -> None:
        missing = ModuleNotFoundError("No module named 'sklearn'", name="sklearn")
        with mock.patch.object(run_module.importlib, "import_module", side_effect=missing):
            with self.assertRaisesRegex(RunError, r'omicsweave\[mfa\]'):
                run_experiment(EXP, ["mfa"], root=self.root)

    def test_all_runs_connected_omics_and_reports_the_rest(self) -> None:
        self.enable_metabolomics()
        run_dir, notes = _quiet(run_experiment, EXP, root=self.root)
        self.assertTrue((run_dir / "omics" / "mfa" / "tables" / "differential.parquet").is_file())
        self.assertFalse((run_dir / "omics" / "metabolomics").exists())
        self.assertTrue(any("metabolomics" in n and "尚未接入" in n for n in notes))
        self.assertEqual([o.value for o in read_manifest(run_dir).omics], ["mfa"])

    def test_steps(self) -> None:
        with self.assertRaisesRegex(RunError, "未知步骤 plots"):
            run_experiment(EXP, ["mfa"], steps=["plots"], root=self.root)
        run_dir, _ = _quiet(run_experiment, EXP, ["mfa"], steps=["qc"], root=self.root)
        mfa = run_dir / "omics" / "mfa"
        self.assertTrue((mfa / "qc" / "fitting_loss.parquet").is_file())
        self.assertFalse((mfa / "tables").exists())
        self.assertEqual(read_manifest(run_dir).parameters["steps"], ["qc"])

    def test_each_run_gets_its_own_directory(self) -> None:
        first, _ = _quiet(run_experiment, EXP, ["mfa"], steps=["qc"], root=self.root, run_id="run_a")
        second, _ = _quiet(run_experiment, EXP, ["mfa"], steps=["qc"], root=self.root, run_id="run_b")
        self.assertEqual(first.parent, second.parent)
        self.assertNotEqual(first, second)
        with self.assertRaisesRegex(RunError, "已存在"):
            run_experiment(EXP, ["mfa"], root=self.root, run_id="run_a")

    def test_cli(self) -> None:
        argv = ["--project-dir", str(self.root), "run", "--experiment", EXP, "--omics", "mfa", "--steps", "qc"]
        out = io.StringIO()
        with contextlib.redirect_stdout(out), self.assertRaises(SystemExit) as exit_:
            cli.main(argv)
        self.assertEqual(exit_.exception.code, 0)
        run_dir = Path(out.getvalue().strip().splitlines()[-1].split("：", 1)[1])
        self.assertEqual(read_manifest(run_dir).command, ["omicsweave", *argv])

        err = io.StringIO()
        with contextlib.redirect_stderr(err), self.assertRaises(SystemExit) as exit_:
            cli.main(["--project-dir", str(self.root), "run", "--experiment", EXP, "--omics", "proteomics"])
        self.assertEqual(exit_.exception.code, 1)
        self.assertIn("未在 experiment.yaml 中启用", err.getvalue())


if __name__ == "__main__":
    unittest.main()
