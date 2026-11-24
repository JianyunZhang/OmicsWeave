#!/usr/bin/env python3
"""结果包工具（G11）：查询接口与 omicsweave diff。

验收：对同一输入的两次运行，diff 报告全部相同；修改一个参数后，diff 只报告受影响的表。
"""

from __future__ import annotations

import contextlib
import io
import shutil
import tempfile
import unittest
from pathlib import Path

import pandas as pd
import yaml

from omicsweave.cli import main
from omicsweave.contract.dataset import ContractError
from omicsweave.contract.diff import diff_bundles, format_report, write_report
from omicsweave.contract.manifest import RunManifest, write_manifest
from omicsweave.contract.query import open_run, read_table
from omicsweave.contract.tables import write_differential
from omicsweave.core.model import Omics
from omicsweave.integration.integrate import integrate_experiment
from tests.integration.test_flux_anchored import contract_tables
from tests.mfa.test_mfa_runner import EXP, build_two_factor_project

PROFILE = "integration/flux_anchored/tables/flux_evidence_profile.parquet"
RELIABILITY = "integration/flux_anchored/qc/contrast_reliability.parquet"


def run_cli(argv: list[str]) -> tuple[int, str]:
    out = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(out):
        try:
            main(argv)
        except SystemExit as exc:
            return int(exc.code or 0), out.getvalue()
    return 0, out.getvalue()


class BundleToolsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        cls.root = Path(cls._tmp.name).resolve()
        build_two_factor_project(cls.root)
        cls.source = cls.root / "results" / EXP / "run_20261101_000000_src"
        for name, table in contract_tables().items():
            write_differential(table, cls.source / "omics" / name / "tables" / "differential.parquet")
        write_manifest(RunManifest(run_id=cls.source.name, project_id="toy", experiment_id=EXP,
                                   omics=list(Omics)), cls.source)
        cls.first = cls.integrate(run_id="run_20261102_000000_aaaa")
        cls.second = cls.integrate(run_id="run_20261102_000001_bbbb")
        path = cls.root / "experiments" / EXP / "experiment.yaml"
        doc = yaml.safe_load(path.read_text(encoding="utf-8"))
        doc["integration"] = {"enabled": False, "modes": ["flux_anchored"], "parameters": {
            "flux_anchored": {"evidence_profile": {"support_threshold": 0.9}}}}
        path.write_text(yaml.safe_dump(doc, sort_keys=False), encoding="utf-8")
        cls.changed = cls.integrate(run_id="run_20261102_000002_cccc")

    @classmethod
    def tearDownClass(cls) -> None:
        cls._tmp.cleanup()

    @classmethod
    def integrate(cls, run_id: str) -> Path:
        with contextlib.redirect_stdout(io.StringIO()):
            run_dir, _ = integrate_experiment(EXP, ["flux_anchored"], root=cls.root, source=cls.source.name,
                                              run_id=run_id)
        return run_dir

    def test_same_input_is_identical(self) -> None:
        diff = diff_bundles(self.first, self.second)
        self.assertTrue(diff.identical, format_report(diff))
        self.assertEqual({t.status for t in diff.tables}, {"identical"})
        code, output = run_cli(["--project-dir", str(self.root), "diff", self.first.name, self.second.name])
        self.assertEqual(code, 0, output)
        self.assertIn("逐值相同", output)

    def test_changed_parameter_reports_only_affected_tables(self) -> None:
        diff = diff_bundles(self.first, self.changed)
        status = {t.path: t.status for t in diff.tables}
        self.assertEqual(status[RELIABILITY], "identical")
        self.assertEqual(status["metadata/samples.parquet"], "identical")
        self.assertEqual(status[PROFILE], "values_differ")
        profile = next(t for t in diff.tables if t.path == PROFILE)
        self.assertIn("evidence_grade", profile.columns_differ)
        self.assertNotIn("support_score", profile.columns_differ)     # 阈值只改变等级，不改变连续分数
        moves = diff.transitions[(diff.transitions["path"] == PROFILE)
                                 & (diff.transitions["column"] == "evidence_grade")]
        self.assertFalse(moves.empty)
        self.assertTrue(set(moves["to"]) <= {"No-signal", "Omics-concordant", "Omics-single-layer",
                                             "Omics-discordant", "Exploratory", "Not-assessable"})
        code, output = run_cli(["--project-dir", str(self.root), "diff", self.first.name, self.changed.name,
                                "--experiment", EXP])
        self.assertEqual(code, 1, output)
        self.assertIn("evidence_grade", output)

    def test_rows_columns_and_report_files(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            copy = Path(tmp) / "copy"
            shutil.copytree(self.first, copy)
            table = pd.read_parquet(copy / PROFILE)
            table.drop(index=0).assign(extra=1.0).to_parquet(copy / PROFILE, index=False)
            diff = diff_bundles(self.first, copy)
            profile = next(t for t in diff.tables if t.path == PROFILE)
            self.assertEqual((profile.status, profile.columns_only_b), ("columns_differ", ["extra"]))
            self.assertEqual((profile.rows_only_a, profile.rows_only_b), (1, 0))
            out = write_report(diff, Path(tmp) / "report")
            self.assertEqual(sorted(p.name for p in out.iterdir()),
                             ["columns.tsv", "report.txt", "tables.tsv", "transitions.tsv"])

    def test_query_latest_and_explicit(self) -> None:
        latest = open_run(EXP, project_dir=self.root)
        self.assertEqual(latest.run_id, self.changed.name)
        diff = read_table(EXP, "differential", omics="mfa", project_dir=self.root)
        self.assertEqual(diff.attrs["run_id"], self.source.name)     # 最近一次含组学结果的运行
        self.assertEqual(set(diff["omics"]), {"mfa"})
        everything = read_table(EXP, "differential", project_dir=self.root)
        self.assertEqual(set(everything["omics"]), {o.value for o in Omics})
        profile = read_table(EXP, "flux_evidence_profile", mode="flux_anchored", run=self.first.name,
                             project_dir=self.root)
        self.assertEqual(profile.attrs["run_id"], self.first.name)
        bundle = open_run(EXP, self.second.name, project_dir=self.root)
        self.assertIn(PROFILE, bundle.tables())
        self.assertEqual(len(bundle.reliability("flux_anchored")), len(pd.read_parquet(self.second / RELIABILITY)))
        with self.assertRaises(ContractError):
            read_table(EXP, "enrichment", project_dir=self.root)
        with self.assertRaises(ContractError):
            open_run(EXP, "run_missing", project_dir=self.root)


if __name__ == "__main__":
    unittest.main()
