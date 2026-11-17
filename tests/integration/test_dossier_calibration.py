#!/usr/bin/env python3
"""判定规则的校准（S4，E8）。

验收（反应档案与分轴综合评估，第 3.11 节与 S4）：
- 诱饵：目标与各次诱饵重复按（组织, 比较）与全部计数，诱饵期望 / 目标为错误比例的估计；
- 对照集：阳性对照看结论或指定轴是否指向预期方向，阴性对照看是否没有给出方向；面板声明错误时报错；
- 不变性：翻转通量的符号，C、R、M 轴不变；删除一整层，不读该层的轴不变；违反时记录格子数；
- 组别置换记为未运行；模式在结果包中写出校准表，图 A 的图注写出诱饵对照。
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from omicsweave.integration.dossier import calibration as cal
from omicsweave.integration.dossier.assessment import assess
from omicsweave.integration.dossier.panels import PanelError, load_panels
from tests.integration import test_dossier_panels as panel_tests     # 只取方法，避免这里再收集一遍它的测试类
from tests.integration.test_dossier_assessment import RULES, TYPICAL, frame
from tests.integration.test_flux_anchored import DIET
from tests.mfa.test_mfa_runner import EXP

PANEL, write_panels = panel_tests.PANEL, panel_tests.write_panels


def axis_rows(flux_id: str, tissue: str = "liver", contrast: str = "C", **states: tuple[str, int]) -> list[dict]:
    """一个（反应, 组织, 比较）的四条轴；没有给出的轴为 undetermined。"""
    return [{"flux_id": flux_id, "tissue": tissue, "contrast_id": contrast, "axis": axis,
             "state": states.get(axis, ("undetermined", 0))[0], "direction": states.get(axis, ("undetermined", 0))[1],
             "effect": np.nan} for axis in ("flux", "capacity", "regulation", "pool")]


def claim(flux_id: str, level: str, direction: int = 0, tissue: str = "liver", contrast: str = "C") -> dict:
    return {"flux_id": flux_id, "tissue": tissue, "contrast_id": contrast, "flux_claim": level,
            "claim_direction": direction}


def assessed(axes: list[list[dict]], claims: list[dict]) -> tuple[pd.DataFrame, pd.DataFrame]:
    return pd.DataFrame([r for rows in axes for r in rows]), pd.DataFrame(claims)


# 目标：肝脏 R1 的 C 轴上升、M 轴同向，结论 L1；R2 不变；心脏 R1 只有 C 轴上升
TARGET = assessed(
    [axis_rows("R1", capacity=("changed", 1), pool=("changed", 1)), axis_rows("R2", capacity=("unchanged", 0)),
     axis_rows("R1", tissue="heart", capacity=("changed", 1))],
    [claim("R1", "L1", 1), claim("R2", "L0"), claim("R1", "L0", tissue="heart")])
# 两次诱饵重复：第一次肝脏 R1 的 C 轴上升但 M 轴反向；第二次什么都没有
DECOYS = [
    assessed([axis_rows("R1", capacity=("changed", 1), pool=("suggestive", -1)), axis_rows("R2"),
              axis_rows("R1", tissue="heart")], [claim("R1", "L0"), claim("R2", "L0"), claim("R1", "L0", tissue="heart")]),
    assessed([axis_rows("R1"), axis_rows("R2"), axis_rows("R1", tissue="heart")],
             [claim("R1", "L0"), claim("R2", "L0"), claim("R1", "L0", tissue="heart")]),
]


class DecoySummaryTest(unittest.TestCase):
    def test_counts_by_group_and_overall(self) -> None:
        table = cal.decoy_summary(TARGET, DECOYS)
        self.assertEqual(list(table.columns), cal.COLUMNS)
        self.assertEqual(len(table), 3 * len(cal.METRICS))                     # 全部 + 两个（组织, 比较）
        rows = table.set_index(["tissue", "contrast_id", "metric"])
        overall = rows.loc[("", "")]
        self.assertEqual(overall.at["capacity_changed", "target"], 2)          # 肝脏 R1 与心脏 R1
        self.assertEqual(overall.at["capacity_changed", "detail"], "1；0")
        self.assertAlmostEqual(overall.at["capacity_changed", "decoy_expected"], 0.5)
        self.assertAlmostEqual(overall.at["capacity_changed", "decoy_sd"], np.std([1, 0], ddof=1))
        self.assertAlmostEqual(overall.at["capacity_changed", "ratio"], 0.25)
        self.assertEqual((overall.at["capacity_and_pool", "target"], overall.at["capacity_and_pool", "decoy_expected"]),
                         (1, 0.0))                                             # 诱饵的 M 轴反向，不算
        self.assertEqual((overall.at["L1", "target"], overall.at["L1", "n_cells"]), (1, 3))
        liver = rows.loc[("liver", "C")]
        self.assertEqual((liver.at["capacity_changed", "target"], liver.at["capacity_changed", "decoy_expected"]),
                         (1, 0.5))
        self.assertAlmostEqual(liver.at["capacity_changed", "ratio"], 0.5)
        heart = rows.loc[("heart", "C")]
        self.assertEqual(heart.at["L1", "target"], 0)
        self.assertTrue(np.isnan(heart.at["L1", "ratio"]))                     # 目标为 0 时不给比例
        self.assertEqual(set(table["n_replicates"]), {2})

    def test_caption(self) -> None:
        table = cal.decoy_summary(TARGET, DECOYS)
        text = cal.caption(table, "liver", "C")
        self.assertTrue(text.startswith("Calibration (decoys ×2): "))
        self.assertIn("C changed target 1, decoy expected 0.5", text)
        self.assertIn("decoys not run", cal.caption(table, "kidney", "C"))
        self.assertEqual(cal.caption(pd.DataFrame(columns=cal.COLUMNS), "liver", "C"), "")


class ControlSummaryTest(unittest.TestCase):
    def test_positive_and_negative_controls(self) -> None:
        axes, reactions = TARGET
        controls = [
            {"reaction": "R1", "kind": "positive", "direction": "up", "axis": "claim", "tissue": "liver", "contrast": "C"},
            {"reaction": "R1", "kind": "positive", "direction": "down", "axis": "capacity", "tissue": "liver",
             "contrast": "C"},
            {"reaction": "R2", "kind": "negative", "axis": "claim", "tissue": "liver", "contrast": "C"},
            {"reaction": "R2", "kind": "negative", "axis": "capacity", "tissue": "liver", "contrast": "C"},
            {"reaction": "R1", "kind": "positive", "direction": "up", "axis": "claim", "tissue": "kidney",
             "contrast": "C"},
        ]
        table = cal.control_summary(reactions, axes, controls)
        single = table[table["metric"] != "concordance"]
        self.assertEqual(list(single["passed"].iloc[:4]), [True, False, True, True])
        self.assertTrue(pd.isna(single["passed"].iloc[4]))
        self.assertEqual(single["detail"].iloc[0], "预期 up，判定 L1 up")
        self.assertEqual(single["detail"].iloc[1], "预期 down，判定 changed up")
        self.assertIn("not_assessed", single["detail"].iloc[4])                # 该组织没有评估：不计入符合率
        concordance = table[table["metric"] == "concordance"].set_index("method")
        self.assertAlmostEqual(concordance.at["control_positive", "target"], 0.5)
        self.assertEqual(concordance.at["control_positive", "detail"], "1/2 个对照符合预期")
        self.assertEqual(concordance.at["control_negative", "target"], 1.0)


class InvarianceTest(unittest.TestCase):
    def run_axes(self, members: pd.DataFrame) -> pd.DataFrame:
        return assess(members, "A", RULES)[0]

    def test_real_assessment_passes(self) -> None:
        table = cal.invariance_checks(frame(*TYPICAL), self.run_axes)
        self.assertEqual(list(table["metric"]), ["flip_flux_sign", "remove_mfa", "remove_proteomics",
                                                 "remove_transcriptomics", "remove_metabolomics"])
        self.assertTrue(table["passed"].all())
        self.assertEqual(set(table["target"]), {0})

    def test_violation_is_counted(self) -> None:
        def leaky(members: pd.DataFrame) -> pd.DataFrame:
            # 故意让 C 轴读取通量：没有通量行时 C 轴改判
            axes = self.run_axes(members)
            if not (members["layer"] == "mfa").any():
                axes.loc[axes["axis"] == "capacity", "state"] = "unchanged"
            return axes

        table = cal.invariance_checks(frame(*TYPICAL), leaky).set_index("metric")
        self.assertFalse(table.at["remove_mfa", "passed"])
        self.assertEqual(table.at["remove_mfa", "target"], 1)
        self.assertEqual(table.at["remove_mfa", "detail"], "R1/liver/C/capacity")
        self.assertTrue(table.at["flip_flux_sign", "passed"])

    def test_permutation_is_recorded_as_not_run(self) -> None:
        row = cal.not_run("E1 未实现").iloc[0]
        self.assertEqual((row["method"], row["detail"]), ("permutation", "E1 未实现"))
        self.assertTrue(pd.isna(row["passed"]))


class ControlDeclarationTest(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.dir = Path(self._tmp.name)

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def load(self, controls) -> object:
        return load_panels(write_panels(self.dir, {"panels": {"p": {
            "tissues": ["liver"], "segments": [{"name": "s", "reactions": ["PC_m", "PDH_m"]}],
            "controls": controls}}}))["p"]

    def test_entries_and_expansion(self) -> None:
        panel = self.load({"positive": [{"reaction": "PC_m", "direction": "up", "axis": "capacity",
                                         "contrast": "HFD_vs_LFD"}],
                           "negative": ["PDH_m"]})
        self.assertEqual(panel.controls, [
            {"reaction": "PC_m", "direction": "up", "axis": "capacity", "contrast": "HFD_vs_LFD", "kind": "positive"},
            {"reaction": "PDH_m", "kind": "negative", "axis": "claim"}])
        expanded = panel.expand_controls(["liver", "heart"], ["HFD_vs_LFD", "HSD_vs_LFD"])
        # 省略的组织取面板声明的（liver），省略的比较取数据中的全部
        self.assertEqual([(c["reaction"], c["tissue"], c["contrast"]) for c in expanded],
                         [("PC_m", "liver", "HFD_vs_LFD"), ("PDH_m", "liver", "HFD_vs_LFD"),
                          ("PDH_m", "liver", "HSD_vs_LFD")])

    def test_declaration_errors(self) -> None:
        cases = {
            "需要 direction": {"positive": ["PC_m"]},
            "axis 应为": {"negative": [{"reaction": "PC_m", "axis": "enzyme"}]},
            "controls 中的反应": {"negative": ["PYK_c"]},
            "controls 应为": {"negative": ["PC_m"], "neutral": ["PDH_m"]},
            "需要 reaction": {"negative": [{"axis": "capacity"}]},
        }
        for message, controls in cases.items():
            with self.subTest(message=message), self.assertRaisesRegex(PanelError, message):
                self.load(controls)


class CalibrationBundleTest(unittest.TestCase):
    """在合成结果包上运行带诱饵与对照的模式。"""
    _bundle = panel_tests.DossierPanelBundleTest
    setUp, tearDown, source_bundle = _bundle.setUp, _bundle.tearDown, _bundle.source_bundle
    set_integration, integrate, run_dossier = _bundle.set_integration, _bundle.integrate, _bundle.run_dossier

    def test_calibration_table_and_notes(self) -> None:
        self.source_bundle()
        doc = {"panels": {name: {**spec, "controls": {
            "positive": [{"reaction": "HEX_c", "direction": "up", "axis": "capacity"}],
            "negative": ["G6Pase_c"]}} for name, spec in PANEL["panels"].items()}}
        write_panels(self.root / "experiments" / EXP, doc)
        base = self.run_dossier(panels="panels.yaml", tracers=["glucose"], label="A", figures={"formats": ["png"]},
                                calibration={"decoys": 2, "seed": 7})
        table = pd.read_parquet(base / "evaluation" / "assessment_calibration.parquet")
        self.assertEqual(list(table.columns), cal.COLUMNS)
        self.assertEqual(set(table["method"]), {"decoy", "control_positive", "control_negative", "invariance",
                                                "permutation"})
        decoys = table[table["method"] == "decoy"]
        self.assertEqual(set(decoys["n_replicates"]), {2})
        self.assertIn(("liver", DIET), set(zip(decoys["tissue"], decoys["contrast_id"])))
        invariance = table[table["method"] == "invariance"]
        self.assertTrue(invariance["passed"].astype(bool).all(), invariance[["metric", "detail"]])
        # 对照的判定与评估表一致
        axes = pd.read_parquet(base / "tables" / "axis_assessment.parquet").set_index(
            ["flux_id", "tissue", "contrast_id", "axis"])
        capacity = axes.loc[("HEX_c", "liver", DIET, "capacity")]
        positive = table[(table["method"] == "control_positive") & (table["metric"] == "HEX_c:capacity")].iloc[0]
        self.assertEqual(bool(positive["passed"]),
                         capacity["state"] in ("changed", "suggestive") and capacity["direction"] == 1)
        self.assertTrue(any(n.startswith("诱饵校准（2 次）") for n in self.notes), self.notes)
        self.assertTrue(any(n.startswith("不变性检查：全部通过") for n in self.notes), self.notes)
        self.assertTrue(cal.caption(table, "liver", DIET).startswith("Calibration (decoys ×2)"))
        # 同一种子重跑，诱饵计数相同
        again = self.run_dossier(panels="panels.yaml", tracers=["glucose"], label="A", figures={"formats": ["png"]},
                                 calibration={"decoys": 2, "seed": 7})
        pd.testing.assert_frame_equal(
            pd.read_parquet(again / "evaluation" / "assessment_calibration.parquet").query("method == 'decoy'"),
            decoys)

    def test_default_runs_no_decoys(self) -> None:
        self.source_bundle()
        base = self.run_dossier(reactions=["HEX_c", "PDH_m"])
        table = pd.read_parquet(base / "evaluation" / "assessment_calibration.parquet")
        self.assertEqual(set(table["method"]), {"invariance", "permutation"})

    def test_unknown_calibration_key(self) -> None:
        self.source_bundle()
        with self.assertRaisesRegex(Exception, "calibration 中有无法识别的项：decoy"):
            self.run_dossier(calibration={"decoy": 1})


if __name__ == "__main__":
    unittest.main()
