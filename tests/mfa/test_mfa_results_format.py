#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""MFA 求解结果标准格式：模板、校验、Excel / Parquet 读写，以及旧版自研软件格式的转换。

转换的正确性以旧版 MFAAnalyser 的读取结果为准：同一份旧文件，经标准格式还原出的
最优拟合通量矩阵、MID 矩阵与按拟合配对的净通量，必须与旧读取逻辑逐值一致。
"""

from __future__ import annotations

import contextlib
import io
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd
import yaml

from omicsweave import cli
from omicsweave.adapters.mfa import convert_legacy_mfa, legacy_flux_matrix, legacy_mid_matrices
from omicsweave.contract import ContractError
from omicsweave.contract import mfa as mfa_contract
from omicsweave.contract.mfa import MFAInfo, MFAResults, read_mfa, validate_mfa, write_mfa, write_template
from omicsweave.core import project
from omicsweave.mfa.df_analysis import _paired_net_fits
from omicsweave.mfa.mfa_analysis import MFAAnalyser

DEMO = Path(__file__).resolve().parents[2] / "examples" / "demo_project"

INFO = pd.DataFrame({
    "complete_result_name": ["s1", "s2"],
    "sheet_name": ["s1_fit", "s2_fit"],
    "condition": ["control", "treat"],
})


def write_legacy_files(root: Path) -> tuple[Path, Path]:
    flux_path, mid_path = root / "flux_raw_data.xlsx", root / "mid_raw_data.xlsx"
    s1 = pd.DataFrame({
        "Index": [0, 1, 2], "loss": [3.0, 1.0, 2.0],
        "liver--PDH_m": [10.0, 20.0, 30.0], "liver--MDH_m": [5.0, 6.0, 7.0], "liver--MDH_m__R": [1.0, 2.0, 4.0],
        "heart--PDH_m": [8.0, 9.0, 10.0], "FA_usage": [0.5, 0.6, 0.7], "MIX:liver--citrate_0": [0.9, 0.8, 0.7],
    })
    s2 = pd.DataFrame({
        "Index": [0, 1, 2], "loss": [1.0, 4.0, 2.5],
        "liver--PDH_m": [11.0, 21.0, 31.0], "liver--MDH_m": [5.5, 6.5, 7.5], "liver--MDH_m__R": [0.5, 3.0, 2.0],
        "heart--PDH_m": [7.0, 9.5, 12.0], "FA_usage": [0.4, 0.9, 0.8], "MIX:liver--citrate_0": [0.6, 0.5, 0.4],
    })
    with pd.ExcelWriter(flux_path) as writer:
        INFO.to_excel(writer, sheet_name="sheet_information", index=False)
        s1.to_excel(writer, sheet_name="s1_fit", index=False)
        s2.to_excel(writer, sheet_name="s2_fit", index=False)

    def mid_sheet(measured, sims, losses):
        names = ["liver--citrate"] * 3 + ["liver--lactate"] * 2
        sheet = pd.DataFrame({"metabolite_name": names + ["loss"], "mid_index": [0, 1, 2, 0, 1, np.nan],
                              "experiments": measured + [np.nan]})
        for fit, values in enumerate(sims):
            sheet[fit] = values + [losses[fit]]
        return sheet

    m1 = mid_sheet([0.6, 0.3, 0.1, 0.7, 0.3],
                   [[0.5, 0.4, 0.1, 0.6, 0.4], [0.62, 0.28, 0.1, 0.71, 0.29], [0.55, 0.35, 0.1, 0.65, 0.35]],
                   [2.0, 0.5, 1.0])
    m2 = mid_sheet([0.5, 0.35, 0.15, 0.8, 0.2],
                   [[0.52, 0.33, 0.15, 0.79, 0.21], [0.4, 0.4, 0.2, 0.7, 0.3], [0.45, 0.4, 0.15, 0.75, 0.25]],
                   [0.3, 1.2, 0.9])
    with pd.ExcelWriter(mid_path) as writer:
        INFO.assign(sheet_name=["s1_mid", "s2_mid"]).to_excel(writer, sheet_name="sheet_information", index=False)
        m1.to_excel(writer, sheet_name="s1_mid", index=False)
        m2.to_excel(writer, sheet_name="s2_mid", index=False)
    return flux_path, mid_path


def single_estimate_results(**flux_overrides) -> MFAResults:
    flux = pd.DataFrame({
        "sample_id": ["a", "a", "b", "b"],
        "tissue": ["liver"] * 4,
        "reaction_id": ["PDH", "LDH", "PDH", "LDH"],
        "direction": ["net"] * 4,
        "value": [1.0, 2.0, 3.0, 4.0],
        "ci_lower": [0.5, 1.5, 2.5, 3.5],
        "ci_upper": [1.5, 2.5, 3.5, 4.5],
    })
    for column, values in flux_overrides.items():
        flux[column] = values
    return MFAResults(info=MFAInfo(model_id="m", flux_units="nmol/min/g"), flux=flux)


class LegacyConversionTest(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name)
        self.flux_path, self.mid_path = write_legacy_files(self.root)
        self.converted = convert_legacy_mfa(self.flux_path, self.mid_path, model_id="toy", flux_units="nmol/min/g")
        self.analyser = MFAAnalyser(output_dir=str(self.root / "out"))
        with contextlib.redirect_stdout(io.StringIO()):
            self.analyser.load_flux_data(str(self.flux_path))
            self.analyser.load_mid_data(str(self.mid_path))

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def test_converted_results_are_valid(self) -> None:
        errors, warnings = validate_mfa(self.converted.results)
        self.assertEqual(errors, [])
        self.assertEqual(warnings, [])
        self.assertEqual(self.converted.conditions, {"s1": "control", "s2": "treat"})
        self.assertTrue(any("FA_usage" in note for note in self.converted.notes))

    def test_best_fit_flux_matches_legacy_loader(self) -> None:
        expected = self.analyser.flux_df.sort_index().sort_index(axis=1)
        actual = legacy_flux_matrix(self.converted.results).sort_index().sort_index(axis=1)
        # 只比较取值；旧版用 object 索引，pandas 3 默认字符串类型，二者不影响数值
        pd.testing.assert_frame_equal(actual, expected, check_names=False,
                                      check_index_type=False, check_column_type=False)

    def test_mid_matrices_match_legacy_loader(self) -> None:
        simulated, measured = legacy_mid_matrices(self.converted.results)
        for actual, expected in ((simulated, self.analyser.mid_df), (measured, self.analyser.mid_measured_df)):
            pd.testing.assert_frame_equal(actual.sort_index().sort_index(axis=1),
                                          expected.sort_index().sort_index(axis=1),
                                          check_names=False, check_dtype=False, check_index_type=False,
                                          check_column_type=False)

    def test_net_flux_is_paired_by_fit(self) -> None:
        net = self.converted.results.net_flux()
        for sample in ("s1", "s2"):
            rows = net[(net["sample_id"] == sample) & (net["tissue"] == "liver") & (net["reaction_id"] == "MDH_m")]
            expected = _paired_net_fits(self.analyser.flux_all_fits[sample], "liver--MDH_m", "liver--MDH_m__R")
            np.testing.assert_allclose(rows.sort_values("fit_id")["value"].to_numpy(), np.asarray(expected))

    def test_mix_columns_are_not_fluxes(self) -> None:
        self.assertFalse(self.converted.results.flux["reaction_id"].str.contains("citrate").any())


class FormatRoundTripTest(unittest.TestCase):
    def test_excel_and_parquet_round_trip(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            flux_path, mid_path = write_legacy_files(Path(tmp))
            results = convert_legacy_mfa(flux_path, mid_path, model_id="toy", flux_units="nmol/min/g").results
            for target in (Path(tmp) / "standard.xlsx", Path(tmp) / "standard_tables"):
                loaded = read_mfa(write_mfa(results, target))
                self.assertEqual(loaded.info, results.info)
                pd.testing.assert_frame_equal(legacy_flux_matrix(loaded), legacy_flux_matrix(results))
                self.assertEqual(len(loaded.mid), len(results.mid))

    def test_excel_row_limit_points_to_parquet(self) -> None:
        with tempfile.TemporaryDirectory() as tmp, mock.patch.object(mfa_contract, "EXCEL_MAX_ROWS", 3):
            with self.assertRaisesRegex(ContractError, "Parquet"):
                write_mfa(single_estimate_results(), Path(tmp) / "x.xlsx")

    def test_template_has_all_sheets_and_empty_template_fails_check(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = write_template(Path(tmp) / "template.xlsx")
            with pd.ExcelFile(path) as xl:
                sheets = xl.sheet_names
            self.assertEqual(sheets, ["info", "flux", "mid", "fits", "说明"])
            with self.assertRaisesRegex(ContractError, "model_id"):
                read_mfa(path)


class ValidationTest(unittest.TestCase):
    def test_single_estimate_with_ci_is_valid(self) -> None:
        results = single_estimate_results()
        self.assertEqual(validate_mfa(results), ([], []))
        best = results.best_estimates().set_index(["sample_id", "reaction_id"])["value"]
        self.assertEqual(best[("b", "LDH")], 4.0)

    def test_value_outside_ci_is_a_warning(self) -> None:
        _, warnings = validate_mfa(single_estimate_results(value=[1.0, 2.0, 9.0, 4.0]))
        self.assertTrue(any("ci_lower" in w for w in warnings))

    def test_structural_errors(self) -> None:
        cases = {
            "direction 取值无法识别": dict(direction=["net", "net", "both", "net"]),
            "同时给出了 net 与正/逆向": dict(direction=["net", "net", "forward", "net"], reaction_id=["PDH", "LDH", "PDH", "PDH"]),
            "只有 reverse 没有 forward": dict(direction=["net", "net", "reverse", "net"]),
            "重复": dict(reaction_id=["PDH", "PDH", "PDH", "LDH"]),
            "tissue 为空": dict(tissue=["liver", "", "liver", "liver"]),
        }
        for fragment, overrides in cases.items():
            with self.subTest(fragment):
                errors, _ = validate_mfa(single_estimate_results(**overrides))
                self.assertTrue(any(fragment in e for e in errors), errors)

    def test_fits_without_loss_entries(self) -> None:
        results = single_estimate_results(fit_id=["1", "1", "1", "2"])
        results.fits = pd.DataFrame({"sample_id": ["a"], "fit_id": ["1"], "loss": [0.1]})
        errors, _ = validate_mfa(results)
        self.assertTrue(any("找不到 loss" in e for e in errors))


class MfaCliTest(unittest.TestCase):
    # cli.main 会设置全局项目目录；测试后恢复，以免影响 pytest --project-dir 指定的真实数据测试
    def setUp(self) -> None:
        self._saved_project_dir = project._explicit_project_dir

    def tearDown(self) -> None:
        project._explicit_project_dir = self._saved_project_dir

    def run_cli(self, argv: list[str]) -> tuple[int, str]:
        out = io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(out):
            try:
                cli.main(argv)
                return 0, out.getvalue()
            except SystemExit as exc:
                return int(exc.code or 0), out.getvalue()

    def test_convert_then_check(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            flux_path, mid_path = write_legacy_files(Path(tmp))
            out = Path(tmp) / "converted"
            code, output = self.run_cli(["mfa", "convert", "--flux", str(flux_path), "--mid", str(mid_path),
                                         "--model-id", "toy", "--flux-units", "nmol/min/g", "-o", str(out)])
            self.assertEqual(code, 0, output)
            self.assertIn("control", output)
            code, output = self.run_cli(["mfa", "check", str(out)])
            self.assertEqual(code, 0, output)
            code, output = self.run_cli(["mfa", "check", str(flux_path)])
            self.assertEqual(code, 0, output)
            self.assertIn("旧版", output)

    def test_template_command(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / "t.xlsx"
            code, output = self.run_cli(["mfa", "template", str(target)])
            self.assertEqual(code, 0, output)
            self.assertTrue(target.is_file())

    def test_experiment_check_matches_samples_table(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "demo"
            shutil.copytree(DEMO, root)
            exp_dir = root / "experiments" / "demo_glp1"
            (exp_dir / "mfa").mkdir()
            write_legacy_files(exp_dir / "mfa")

            config = yaml.safe_load((exp_dir / "experiment.yaml").read_text(encoding="utf-8"))
            config["omics"]["mfa"] = {"enabled": True, "input": {
                "level": "table", "format": "legacy_mfa_xlsx", "model_id": "toy", "flux_units": "nmol/min/g",
                "files": ["experiments/demo_glp1/mfa/flux_raw_data.xlsx", "experiments/demo_glp1/mfa/mid_raw_data.xlsx"],
            }}
            (exp_dir / "experiment.yaml").write_text(yaml.safe_dump(config, allow_unicode=True), encoding="utf-8")
            samples = exp_dir / "samples.tsv"
            base = samples.read_text(encoding="utf-8")
            samples.write_text(base + "s1\tm01\tmfa\tliver\t4\tb1\tsample\t\ns2\tm05\tmfa\tliver\t4\tb1\tsample\t\n",
                               encoding="utf-8")

            code, output = self.run_cli(["--project-dir", str(root), "validate"])
            self.assertEqual(code, 0, output)
            code, output = self.run_cli(["--project-dir", str(root), "mfa", "check", "--experiment", "demo_glp1"])
            self.assertEqual(code, 0, output)

            samples.write_text(base + "s1\tm01\tmfa\tliver\t4\tb1\tsample\t\n", encoding="utf-8")
            code, output = self.run_cli(["--project-dir", str(root), "mfa", "check", "--experiment", "demo_glp1"])
            self.assertEqual(code, 1)
            self.assertIn("不在 samples.tsv", output)

            del config["omics"]["mfa"]["input"]["model_id"]
            (exp_dir / "experiment.yaml").write_text(yaml.safe_dump(config, allow_unicode=True), encoding="utf-8")
            code, output = self.run_cli(["--project-dir", str(root), "validate"])
            self.assertEqual(code, 1)
            self.assertIn("需要在 mfa.input 中给出 model_id", output)


if __name__ == "__main__":
    unittest.main()
