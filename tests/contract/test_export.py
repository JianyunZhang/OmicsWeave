#!/usr/bin/env python3
"""omicsweave export：结果包中的统一结果表导出为 Excel / TSV，默认取最近一次运行。"""

from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

from omicsweave import cli
from omicsweave.contract import export as export_module
from omicsweave.contract.dataset import ContractError
from omicsweave.contract.export import export_bundle, find_run
from omicsweave.core.run import run_experiment
from tests.mfa.test_mfa_runner import EXP, build_two_factor_project


class ExportTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        cls.root = Path(cls._tmp.name).resolve()
        build_two_factor_project(cls.root)
        with contextlib.redirect_stdout(io.StringIO()):
            cls.first, _ = run_experiment(EXP, ["mfa"], root=cls.root, run_id="run_20261101_000000_aaaa")
            cls.latest, _ = run_experiment(EXP, ["mfa"], root=cls.root, run_id="run_20261101_000001_bbbb")
        cls.results = cls.root / "analysis_results"

    @classmethod
    def tearDownClass(cls) -> None:
        cls._tmp.cleanup()

    def test_find_run(self) -> None:
        self.assertEqual(find_run(self.results, EXP), self.latest)
        self.assertEqual(find_run(self.results, EXP, self.first.name), self.first)
        with self.assertRaisesRegex(ContractError, "不是完整的结果包"):
            find_run(self.results, EXP, "run_missing")
        with self.assertRaisesRegex(ContractError, "omicsweave run"):
            find_run(self.results, "other_exp")

    def test_xlsx(self) -> None:
        path = export_bundle(self.latest)
        self.assertEqual(path, self.latest / "exports" / f"{EXP}_{self.latest.name}.xlsx")
        sheets = pd.read_excel(path, sheet_name=None)
        # MFA 结果包另有逐只动物的拟合分布与通量组（G10）
        self.assertEqual(list(sheets), ["info", "mfa_differential", "mfa_flux_groups", "mfa_sample_flux_fits"])
        info = dict(zip(sheets["info"]["field"], sheets["info"]["value"]))
        self.assertEqual(info["run_id"], self.latest.name)
        self.assertEqual(int(info["rows:mfa_differential"]), len(sheets["mfa_differential"]))
        expected = pd.read_parquet(self.latest / "omics" / "mfa" / "tables" / "differential.parquet")
        self.assertEqual(list(sheets["mfa_differential"].columns), list(expected.columns))

    def test_tsv_and_row_limit(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = export_bundle(self.first, "tsv", Path(tmp) / "out")
            self.assertEqual(sorted(p.name for p in directory.iterdir()),
                             ["info.tsv", "mfa_differential.tsv", "mfa_flux_groups.tsv", "mfa_sample_flux_fits.tsv"])
            table = pd.read_csv(directory / "mfa_differential.tsv", sep="\t")
            self.assertIn("contrast_id", table.columns)
        with mock.patch.object(export_module, "EXCEL_MAX_ROWS", 3):
            with self.assertRaisesRegex(ContractError, "tsv"):
                export_bundle(self.first, "xlsx", self.root / "too_big.xlsx")

    def test_cli(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / "export.xlsx"
            out = io.StringIO()
            with contextlib.redirect_stdout(out), self.assertRaises(SystemExit) as exit_:
                cli.main(["--project-dir", str(self.root), "export", "--experiment", EXP, "--output", str(target)])
            self.assertEqual(exit_.exception.code, 0)
            self.assertTrue(target.is_file())
            self.assertIn("已导出", out.getvalue())
        err = io.StringIO()
        with contextlib.redirect_stderr(err), self.assertRaises(SystemExit) as exit_:
            cli.main(["--project-dir", str(self.root), "export", "--experiment", "other_exp"])
        self.assertEqual(exit_.exception.code, 1)
        self.assertIn("omicsweave run", err.getvalue())


if __name__ == "__main__":
    unittest.main()
