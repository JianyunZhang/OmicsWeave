#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""回归测试：灌注血糖 QC 应同时识别双语长列名与纯时间表头（``0min``），
并对 ``21..5`` 这类录入笔误规范化而不是静默丢弃。"""

import contextlib
import io
import tempfile
import unittest
from pathlib import Path

import pandas as pd



from omicsweave.mfa.infusion_qc import _parse_timepoint, infusion_steady_state_qc


MINUTES = [0, 30, 60, 90, 120, 150, 180]
# 小鼠 1 平稳；小鼠 2 持续上升，150 min 处带重复小数点笔误；
# 小鼠 3 的 90 min 是无法解析的文本
GLUCOSE = {
    "1": [8.0, 8.1, 8.2, 8.0, 8.1, 8.3, 8.2],
    "2": [8.0, 10.0, 12.0, 14.0, 16.0, "18..5", 20.0],
    "3": [9.0, 9.1, 9.0, "x", 9.2, 9.1, 9.0],
}
QC_CONFIG = {
    "enabled": True,
    "steady_state_window_min": [30, 180],
    "glucose_cv_max": 0.15,
    "glucose_range_max_mmol_l": 2.0,
    "on_fail": "flag",
}
GROUPS = {"1": "ctrl", "2": "treat", "3": "ctrl"}


def _legacy_header(minute):
    return f"{minute} min (blood_glucose_after_infusion_{minute}_min) (mmol/L)"


def _bare_header(minute):
    return f"{minute}min"


def _write_record(path, header_fn):
    data = {"编号 (mice_id)": list(GLUCOSE)}
    data["灌注2h血酮"] = [0.5, 0.6, 0.7]  # 非时间点列，不应被识别
    for i, minute in enumerate(MINUTES):
        data[header_fn(minute)] = [GLUCOSE[m][i] for m in GLUCOSE]
    pd.DataFrame(data).to_excel(path, index=False)


def _run_qc(temp_path, header_fn):
    record = temp_path / f"glycemia_{header_fn.__name__}.xlsx"
    _write_record(record, header_fn)
    out_dir = temp_path / header_fn.__name__
    log = io.StringIO()
    with contextlib.redirect_stdout(log):
        qc = infusion_steady_state_qc(record, out_dir, QC_CONFIG, GROUPS)
    return qc, out_dir / "01_df_results" / "all_tissues" / "qc", log.getvalue()


class ParseTimepointTest(unittest.TestCase):
    def test_both_header_styles(self):
        self.assertEqual(120, _parse_timepoint("blood_glucose_after_infusion_120_min"))
        self.assertEqual(0, _parse_timepoint("0min"))
        self.assertEqual(180, _parse_timepoint("180min"))
        self.assertEqual(15, _parse_timepoint(" 15 MIN "))

    def test_non_timepoint_columns_rejected(self):
        for column in ["灌注2h血酮", "mice_id", "tc/hdl", "120min.1", "fast_120min_note"]:
            self.assertIsNone(_parse_timepoint(column), column)


class InfusionQcHeaderStylesTest(unittest.TestCase):
    def test_legacy_and_bare_headers_give_identical_qc(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as temp_dir:
            temp_path = Path(temp_dir)
            legacy, legacy_dir, _ = _run_qc(temp_path, _legacy_header)
            bare, bare_dir, _ = _run_qc(temp_path, _bare_header)

            self.assertIsNotNone(legacy)
            self.assertIsNotNone(bare)
            for qc_dir in (legacy_dir, bare_dir):
                self.assertTrue((qc_dir / "infusion_steady_state_qc.csv").exists())
                self.assertTrue((qc_dir / "infusion_glycemia_timecourse.csv").exists())
            pd.testing.assert_frame_equal(legacy, bare)

            by_mouse = bare.set_index("mouse_id")
            self.assertEqual(6, by_mouse.loc["1", "n_timepoints"])
            self.assertTrue(by_mouse.loc["1", "passed"])
            self.assertFalse(by_mouse.loc["2", "passed"])

            timecourse = pd.read_csv(bare_dir / "infusion_glycemia_timecourse.csv",
                                     dtype={"mouse_id": str})
            self.assertEqual(sorted(MINUTES), sorted(timecourse["minute"].unique()))

    def test_repeated_decimal_is_normalized_with_warning(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as temp_dir:
            qc, qc_dir, log = _run_qc(Path(temp_dir), _bare_header)

            by_mouse = qc.set_index("mouse_id")
            # 18..5 → 18.5 计入窗口：30–180 min 共 6 点，最大值 20.0、最小值 10.0
            self.assertEqual(6, by_mouse.loc["2", "n_timepoints"])
            self.assertAlmostEqual(10.0, by_mouse.loc["2", "glucose_range_mmol_l"])
            self.assertEqual("150min:normalized_from:18..5", by_mouse.loc["2", "value_notes"])
            self.assertIn("[WARN] 小鼠 2 @ 150 min", log)

            timecourse = pd.read_csv(qc_dir / "infusion_glycemia_timecourse.csv",
                                     dtype={"mouse_id": str}, keep_default_na=False)
            cell = timecourse[(timecourse["mouse_id"] == "2") & (timecourse["minute"] == 150)]
            self.assertAlmostEqual(18.5, float(cell["glucose_mmol_l"].iloc[0]))
            self.assertEqual("normalized_from:18..5", cell["value_note"].iloc[0])

    def test_unparseable_value_is_dropped_with_warning(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as temp_dir:
            qc, _, log = _run_qc(Path(temp_dir), _bare_header)

            by_mouse = qc.set_index("mouse_id")
            self.assertEqual(5, by_mouse.loc["3", "n_timepoints"])
            self.assertEqual("90min:unparseable:x", by_mouse.loc["3", "value_notes"])
            self.assertIn("[WARN] 小鼠 3 @ 90 min", log)
            self.assertEqual("", by_mouse.loc["1", "value_notes"])


if __name__ == "__main__":
    unittest.main()
