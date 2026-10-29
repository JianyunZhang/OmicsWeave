#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""回归测试：灌注血糖 QC 的「动物编号 → 条件」映射。

血糖表按动物编号记录。``20250613_cimr_wt_mice_dio_glp`` 的 ``sample_groups``
以 MFA 样本名为键，动物编号只在 ``mouse_groups`` 里；此前管线只传
``sample_groups``，整张表的动物都被标成 ``unmapped``。这里覆盖：

- 配置加载器读出 ``mouse_groups``；
- 优先 ``mouse_groups``、缺失时退回 ``sample_groups``（20260312 批次的写法）；
- 不在映射中的动物（GLP 批次里改输谷氨酰胺的 85 号）另行报告、不计入判定；
- 映射与血糖表完全没有交集时不拆分，并打印醒目警告；
- 从 ``load_mfa_config`` 经 ``run_from_mfa_config`` 的 run_context 到 QC 的整条链路。
"""

import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd
import yaml



from omicsweave.mfa import run_mfa_analysis
from omicsweave.mfa.infusion_qc import (
    _normalize_mouse_id,
    infusion_steady_state_qc,
    select_animal_groups,
)
from omicsweave.mfa.mfa_experiment_config import load_mfa_config


MINUTES = [0, 30, 60, 90, 120, 150, 180]
# 1、3 平稳；2 持续上升；85 不属于本批次且没有任何血糖值
GLUCOSE = {
    1: [8.0, 8.1, 8.2, 8.0, 8.1, 8.3, 8.2],
    2: [8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0],
    3: [9.0, 9.1, 9.0, 9.1, 9.2, 9.1, 9.0],
    85: [None] * len(MINUTES),
}
QC_CONFIG = {
    "enabled": True,
    "steady_state_window_min": [30, 180],
    "glucose_cv_max": 0.15,
    "glucose_range_max_mmol_l": 2.0,
    "on_fail": "flag",
}
# 动物编号 → 条件（YAML 里写成整数键，加载后应为字符串）
MOUSE_GROUPS = {1: "Ctrl", 2: "Treat", 3: "Ctrl"}
# MFA 样本名 → 条件：与血糖表的动物编号没有任何交集
SAMPLE_NAME_GROUPS = {"ctrl_1": "ctrl", "treat_1": "treat", "ctrl_2": "ctrl"}
EXPECTED_CONDITIONS = {"1": "ctrl", "2": "treat", "3": "ctrl"}


def _write_record(path):
    """写一张血糖表；末尾一行只有备注、没有编号，使编号列被 Excel 读成 float。"""
    ids = list(GLUCOSE) + [None]
    data = {"耳标 (mice_id)": ids, "备注": [""] * len(GLUCOSE) + ["记录人备注"]}
    for i, minute in enumerate(MINUTES):
        data[f"{minute}min"] = [GLUCOSE[m][i] for m in GLUCOSE] + [None]
    pd.DataFrame(data).to_excel(path, index=False)


def _qc_dir(output_dir):
    return Path(output_dir) / "01_df_results" / "all_tissues" / "qc"


def _run_qc(temp_path, groups):
    record = temp_path / "glycemia.xlsx"
    _write_record(record)
    out_dir = temp_path / "out"
    log = io.StringIO()
    with contextlib.redirect_stdout(log):
        qc = infusion_steady_state_qc(record, out_dir, QC_CONFIG, groups)
    return qc, _qc_dir(out_dir), log.getvalue()


def _read_csv(path):
    return pd.read_csv(path, dtype={"mouse_id": str}, keep_default_na=False)


class MouseIdAndSelectionTest(unittest.TestCase):
    def test_integer_valued_float_ids_lose_trailing_zero(self):
        self.assertEqual("86", _normalize_mouse_id(86.0))
        self.assertEqual("86", _normalize_mouse_id(86))
        self.assertEqual("A12", _normalize_mouse_id(" A12 "))
        self.assertEqual("12.5", _normalize_mouse_id(12.5))

    def test_mouse_groups_preferred_over_sample_groups(self):
        chosen = select_animal_groups({"86": "hfd"}, {"hfd_1": "hfd"})
        self.assertEqual({"86": "hfd"}, chosen)

    def test_falls_back_to_sample_groups(self):
        # 20260312 的批次只有 sample_groups，且以动物编号为键
        for empty in (None, {}):
            chosen = select_animal_groups(empty, {"A12": "hfd_sham_8w"})
            self.assertEqual({"A12": "hfd_sham_8w"}, chosen)
        self.assertEqual({}, select_animal_groups(None, None))

    def test_keys_are_stringified(self):
        self.assertEqual({"86": "hfd"}, select_animal_groups({86: "hfd"}, None))


class UnmappedAnimalsTest(unittest.TestCase):
    def test_unmapped_animals_reported_separately(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as temp_dir:
            groups = select_animal_groups(MOUSE_GROUPS, SAMPLE_NAME_GROUPS)
            qc, qc_dir, log = _run_qc(Path(temp_dir), {k: v.lower() for k, v in groups.items()})

            # 编号列被读成 float 也能对上 "1"/"2"/"3"
            self.assertEqual(EXPECTED_CONDITIONS, dict(zip(qc["mouse_id"], qc["condition"])))
            self.assertNotIn("unmapped", set(qc["condition"]))

            unmapped = _read_csv(qc_dir / "infusion_steady_state_unmapped.csv")
            self.assertEqual(["85"], list(unmapped["mouse_id"]))
            self.assertEqual(["unmapped"], list(unmapped["condition"]))
            self.assertIn("[WARN] 1 只动物在血糖表中但不在动物分组中", log)

            main = _read_csv(qc_dir / "infusion_steady_state_qc.csv")
            self.assertEqual({"1", "2", "3"}, set(main["mouse_id"]))
            timecourse = _read_csv(qc_dir / "infusion_glycemia_timecourse.csv")
            self.assertNotIn("85", set(timecourse["mouse_id"]))
            by_condition = _read_csv(qc_dir / "infusion_steady_state_by_condition.csv")
            self.assertEqual({"ctrl", "treat"}, set(by_condition["condition"]))
            self.assertIn("稳态通过: 2/3 只动物", log)

    def test_no_overlap_keeps_everyone_and_warns(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as temp_dir:
            qc, qc_dir, log = _run_qc(Path(temp_dir), SAMPLE_NAME_GROUPS)

            self.assertEqual({"1", "2", "3", "85"}, set(qc["mouse_id"]))
            self.assertEqual({"unmapped"}, set(qc["condition"]))
            self.assertIn("完全没有交集", log)
            self.assertIn("mouse_groups", log)
            # 只写表头，不残留旧结果
            unmapped = _read_csv(qc_dir / "infusion_steady_state_unmapped.csv")
            self.assertTrue(unmapped.empty)
            self.assertIn("mouse_id", unmapped.columns)


class PipelineWiringTest(unittest.TestCase):
    """load_mfa_config → run_from_mfa_config 的 run_context → MFAAnalyser → QC。"""

    def _make_batch(self, root, groups_doc):
        mfa_dir = root / "experiments" / "test_exp" / "mfa_data"
        batch_dir = mfa_dir / "u-13c_glucose_infusion" / "batch_a"
        batch_dir.mkdir(parents=True)
        (mfa_dir / "u-13c_glucose_infusion" / "protocol.yaml").write_text(
            yaml.safe_dump({"infusion_qc_defaults": QC_CONFIG}), encoding="utf-8"
        )
        record = batch_dir / "glycemia.xlsx"
        _write_record(record)
        config = {
            "experiment_id": "test_exp",
            "groups_file": "groups.yaml",
            "paths": {
                "output_dir": str(root / "mfa_results"),
                "glycemia_record_file": str(record),
            },
        }
        cfg_path = batch_dir / "config.yaml"
        cfg_path.write_text(yaml.safe_dump(config), encoding="utf-8")
        (batch_dir / "groups.yaml").write_text(
            yaml.safe_dump(groups_doc, allow_unicode=True), encoding="utf-8"
        )
        return cfg_path

    def _run_qc_only(self, cfg_path):
        cfg = load_mfa_config(config_path=str(cfg_path))

        def qc_only(analyser, comparisons=None, mode="default"):
            analyser.infusion_steady_state_qc()

        with mock.patch.object(run_mfa_analysis, "_run_pipeline", qc_only), \
                contextlib.redirect_stdout(io.StringIO()):
            run_mfa_analysis.run_from_mfa_config(cfg)
        return cfg, _read_csv(_qc_dir(cfg.output_dir) / "infusion_steady_state_qc.csv")

    def test_glp_layout_uses_mouse_groups(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as temp_dir:
            cfg_path = self._make_batch(Path(temp_dir), {
                "conditions": {"ctrl": {"display": "Ctrl"}, "treat": {"display": "Treat"}},
                "sample_groups": SAMPLE_NAME_GROUPS,
                "mouse_groups": MOUSE_GROUPS,
            })
            cfg, qc = self._run_qc_only(cfg_path)

            self.assertEqual(EXPECTED_CONDITIONS, cfg.mouse_groups)
            self.assertEqual(SAMPLE_NAME_GROUPS, cfg.sample_groups)
            self.assertEqual(EXPECTED_CONDITIONS, dict(zip(qc["mouse_id"], qc["condition"])))

    def test_sleeve_gastrectomy_layout_falls_back_to_sample_groups(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as temp_dir:
            cfg_path = self._make_batch(Path(temp_dir), {
                "conditions": {"ctrl": {"display": "Ctrl"}, "treat": {"display": "Treat"}},
                "sample_groups": MOUSE_GROUPS,
            })
            cfg, qc = self._run_qc_only(cfg_path)

            self.assertEqual({}, cfg.mouse_groups)
            self.assertEqual(EXPECTED_CONDITIONS, dict(zip(qc["mouse_id"], qc["condition"])))


if __name__ == "__main__":
    unittest.main()
