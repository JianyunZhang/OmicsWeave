#!/usr/bin/env python3
"""通路层的综合读数（S5，G16）。

验收（反应档案与分轴综合评估，第 5 节 G16）：
- 合成的分支节点中，设定的分流变化在酶层与代谢物层被检出，两个出口同比例变化时不被检出；
- 带符号的模块在符号全部取反后 p 值不变、方向相反；
- 交叉点：上游升、下游降为 restricted，流动方向相反时模式随之互换；两侧为同一代谢物的步骤不检验。
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from omicsweave.core.statistics.rotation import signed_rotation_test
from omicsweave.integration.dossier import pathway as pw
from omicsweave.integration.dossier.panels import PanelError, load_panels
from tests.integration import test_dossier_panels as panel_tests     # 只取方法，避免这里再收集一遍它的测试类
from tests.integration.test_flux_anchored import DIET
from tests.mfa.test_mfa_runner import EXP

write_panels = panel_tests.write_panels
NUM, DEN = [f"n{i}" for i in range(3)], [f"d{i}" for i in range(3)]
NOISE = np.array([0.05, -0.08, 0.03, -0.02, 0.07, -0.05])
G6P_PANEL = {"panels": {"g6p": {
    "tissues": ["liver"], "segments": [{"name": "s", "reactions": ["HEX_c", "G6Pase_c", "PGI_c"]}],
    "branch_nodes": [{"node": "GLC6P_c", "exits": {"release": ["G6Pase_c"], "glycolysis": ["PGI_c"]}}],
    "signed_modules": {"release_over_glycolysis": {"plus": ["G6pc1"], "minus": ["Gpi1"]}}}}}


def member(**fields) -> dict:
    base = {"source": "A", "tissue": "liver", "contrast_id": "C", "flux_id": "G6Pase_c", "layer": "proteomics",
            "role": "direct_enzyme", "member_id": "G6pc1", "feature_id": "G6pc1", "member_relation": "isozyme",
            "cell_state": "no_evidence", "tracer": "", "log2fc": np.nan, "pvalue": np.nan,
            "mean_numerator": np.nan, "mean_denominator": np.nan}
    return {**base, **fields}


def per_animal(layer: str, feature: str, shift: float, base: float = 10.0, tracer: str = "") -> list[dict]:
    """分子组比分母组高 ``shift``（log2，或通量层的倍数）。"""
    rows = []
    for i, sample in enumerate(NUM + DEN):
        value = base + (shift if sample in NUM else 0.0) + NOISE[i]
        if layer == "mfa":
            value = base * (2 ** (shift if sample in NUM else 0.0)) * (1 + NOISE[i])
        rows.append({"source": "A", "layer": layer, "tissue": "liver", "feature_id": feature, "tracer": tracer,
                     "sample_id": f"{tracer}{sample}", "subject_id": sample, "group": "H" if sample in NUM else "L",
                     "group_order": 1 if sample in NUM else 0, "value": value})
    return rows


def sides(*layers: str, tracer: str = "") -> pd.DataFrame:
    return pd.DataFrame([{"source": "A", "layer": layer, "contrast_id": "C", "sample_id": f"{tracer if layer == 'mfa' else ''}{s}",
                          "side": "numerator" if s in NUM else "denominator", "level": "H" if s in NUM else "L"}
                         for layer in layers for s in NUM + DEN])


class RotationTest(unittest.TestCase):
    X = np.column_stack([np.ones(6), [1, 1, 1, 0, 0, 0]])

    def test_sign_flip_keeps_p_and_reverses_direction(self) -> None:
        rng = np.random.default_rng(3)
        y = rng.standard_normal((4, 6))
        signs = np.array([1, 1, -1, -1])
        y[:, :3] += 1.2 * signs[:, None]
        up = signed_rotation_test(y, self.X, np.array([0, 1]), signs, rotations=4999, seed=5)
        down = signed_rotation_test(y, self.X, np.array([0, 1]), -signs, rotations=4999, seed=5)
        self.assertEqual(up.pvalue, down.pvalue)
        self.assertEqual((up.direction, down.direction), (1, -1))
        self.assertAlmostEqual(up.statistic, -down.statistic)
        self.assertEqual((up.pvalue_up, up.pvalue_down), (down.pvalue_down, down.pvalue_up))
        self.assertLess(up.pvalue, 0.05)

    def test_null_is_calibrated(self) -> None:
        rng = np.random.default_rng(0)
        p = np.array([signed_rotation_test(rng.standard_normal((5, 6)), self.X, np.array([0, 1]),
                                           np.array([1, 1, 1, -1, -1]), rotations=499, seed=i).pvalue
                      for i in range(300)])
        self.assertTrue(0.02 <= (p < 0.05).mean() <= 0.09, (p < 0.05).mean())
        self.assertTrue(0.4 <= p.mean() <= 0.6, p.mean())

    def test_prior_and_input_checks(self) -> None:
        y = np.random.default_rng(1).standard_normal((3, 6))
        with_prior = signed_rotation_test(y, self.X, np.array([0, 1]), np.ones(3), prior_variance=1.0, prior_df=8,
                                          rotations=99, seed=1)
        self.assertEqual(with_prior.df, 12.0)                              # d = 4，d0 = 8
        y[0, 0] = np.nan
        with self.assertRaisesRegex(ValueError, "完整数据"):
            signed_rotation_test(y, self.X, np.array([0, 1]), np.ones(3), rotations=99)


class NodeSplitTest(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.panel = load_panels(write_panels(Path(self._tmp.name), G6P_PANEL))["g6p"]

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def members(self) -> pd.DataFrame:
        return pd.DataFrame([
            member(), member(flux_id="PGI_c", member_id="Gpi1", feature_id="Gpi1"),
            # 两个出口共有的特征：从两边去掉
            member(member_id="Shared", feature_id="Shared"), member(flux_id="PGI_c", member_id="Shared",
                                                                   feature_id="Shared"),
            member(layer="metabolomics", role="substrate", member_id="Glucose-6-phosphate", feature_id="G6P"),
            member(layer="metabolomics", role="product", member_id="Glucose", feature_id="GLC"),
            member(flux_id="PGI_c", layer="metabolomics", role="substrate", member_id="Glucose-6-phosphate",
                   feature_id="G6P"),
            member(flux_id="PGI_c", layer="metabolomics", role="product", member_id="Fructose-6-phosphate",
                   feature_id="F6P"),
            member(layer="mfa", role="flux", member_id="G6Pase_c", feature_id="G6Pase_c", tracer="glucose",
                   cell_state="no_evidence"),
            member(flux_id="PGI_c", layer="mfa", role="flux", member_id="PGI_c", feature_id="PGI_c", tracer="glucose",
                   cell_state="not_identifiable"),
        ])

    def run_split(self, release: float, glycolysis: float) -> pd.DataFrame:
        samples = pd.DataFrame(
            per_animal("proteomics", "G6pc1", release) + per_animal("proteomics", "Gpi1", glycolysis)
            + per_animal("proteomics", "Shared", 2.0) + per_animal("metabolomics", "GLC", release)
            + per_animal("metabolomics", "F6P", glycolysis) + per_animal("metabolomics", "G6P", 3.0)
            + per_animal("mfa", "G6Pase_c", release, tracer="glucose")
            + per_animal("mfa", "PGI_c", glycolysis, tracer="glucose"))
        table, per = pw.node_splits(self.panel, self.members(), samples,
                                    pd.concat([sides("proteomics", "metabolomics"), sides("mfa", tracer="glucose")]),
                                    ["A"])
        self.per_animal = per
        return table.set_index("layer")

    def test_split_change_is_detected(self) -> None:
        table = self.run_split(release=1.0, glycolysis=0.0)
        for layer in ("proteomics", "metabolomics"):
            row = table.loc[layer]
            self.assertEqual(row["state"], "computed")
            self.assertAlmostEqual(row["log_ratio_change"], 1.0, delta=0.15)
            self.assertLess(row["pvalue"], 0.01)
        protein = table.loc["proteomics"]
        self.assertEqual((protein["members_a"], protein["members_b"]), ("G6pc1", "Gpi1"))
        self.assertIn("共有、已去掉：Shared", protein["detail"])
        self.assertEqual(table.loc["metabolomics", "members_a"], "Glucose")   # 节点自身（G6P）不计入
        # 通量：一个出口在该示踪剂上不可辨识
        self.assertEqual(table.loc["mfa", "state"], "not_identifiable")
        self.assertIn("PGI_c not_identifiable", table.loc["mfa", "detail"])
        values = self.per_animal[self.per_animal["layer"] == "proteomics"].set_index("sample_id")["value"]
        self.assertEqual(len(values), 6)

    def test_proportional_change_is_not_detected(self) -> None:
        table = self.run_split(release=1.0, glycolysis=1.0)
        for layer in ("proteomics", "metabolomics"):
            self.assertAlmostEqual(table.loc[layer, "log_ratio_change"], 0.0, delta=0.15)
            self.assertGreater(table.loc[layer, "pvalue"], 0.2)

    def test_missing_exit_member(self) -> None:
        members = self.members()
        members = members[~((members["flux_id"] == "PGI_c") & (members["layer"] == "metabolomics")
                            & (members["role"] == "product"))]
        samples = pd.DataFrame(per_animal("metabolomics", "GLC", 1.0))
        table, _ = pw.node_splits(self.panel, members, samples, sides("metabolomics"), ["A"])
        row = table.set_index("layer").loc["metabolomics"]
        self.assertEqual(row["state"], "not_detected")
        self.assertIn("出口 b", row["detail"])


class SignedModuleTest(unittest.TestCase):
    def test_flipped_declaration(self) -> None:
        rng = np.random.default_rng(11)
        genes = [f"g{i}" for i in range(200)] + ["G6pc1", "Gpi1"]
        values = pd.DataFrame(rng.normal(8, 0.3, (len(genes), 6)), index=genes, columns=NUM + DEN)
        values.loc["G6pc1", NUM] += 1.0
        values.loc["Gpi1", NUM] -= 1.0
        data = pw.ModuleData(values, lambda wanted: {g: [g] if g in values.index else [] for g in wanted})
        results = []
        for plus, minus in ((["G6pc1"], ["Gpi1"]), (["Gpi1"], ["G6pc1"])):
            panel = pw.Panel("p", "p", [("s", ["G6Pase_c"])], signed_modules={"m": {"plus": plus, "minus": minus}})
            table = pw.signed_module_tests(panel, lambda s, layer, t: data if layer == "proteomics" else None,
                                           sides("proteomics"), ["A"], {"A": ["liver"]}, {"A": ["C"]}, 1999, 7)
            results.append(table.iloc[0])
        declared, flipped = results
        self.assertEqual(declared["pvalue"], flipped["pvalue"])
        self.assertEqual((declared["direction"], flipped["direction"]), ("as_declared", "opposite"))
        self.assertLess(declared["pvalue"], 0.1)
        self.assertEqual(declared["members"], "G6pc1(+), Gpi1(−)")

    def test_unmatched_members_are_reported(self) -> None:
        values = pd.DataFrame(np.random.default_rng(2).normal(8, 0.3, (20, 6)), index=[f"g{i}" for i in range(20)],
                              columns=NUM + DEN)
        data = pw.ModuleData(values, lambda wanted: {g: [g] if g in values.index else [] for g in wanted})
        panel = pw.Panel("p", "p", [("s", ["G6Pase_c"])], signed_modules={"m": {"plus": ["g1", "Absent"], "minus": []}})
        row = pw.signed_module_tests(panel, lambda s, layer, t: data, sides("proteomics", "transcriptomics"), ["A"],
                                     {"A": ["liver"]}, {"A": ["C"]}, 99, 1).iloc[0]
        self.assertEqual(row["direction"], "not_tested")
        self.assertIn("未匹配：Absent", row["detail"])


class CrossoverTest(unittest.TestCase):
    PANEL = pw.Panel("p", "p", [("s", ["PYR_trans", "G6Pase_c", "PYK_c"])])

    def members(self, flux_means: tuple[float, float] | None = None) -> pd.DataFrame:
        rows = [member(layer="metabolomics", role="substrate", member_id="Glucose-6-phosphate", feature_id="G6P",
                       log2fc=0.8, pvalue=0.002),
                member(layer="metabolomics", role="product", member_id="Glucose", feature_id="GLC", log2fc=-0.9,
                       pvalue=0.004),
                member(flux_id="PYK_c", layer="metabolomics", role="substrate", member_id="Phosphoenolpyruvate",
                       feature_id="PEP", log2fc=0.5, pvalue=0.01),
                member(flux_id="PYR_trans", layer="metabolomics", role="substrate", member_id="Pyruvate",
                       feature_id="PYR", log2fc=0.4, pvalue=0.2),
                member(flux_id="PYR_trans", layer="metabolomics", role="product", member_id="Pyruvate",
                       feature_id="PYR", log2fc=0.4, pvalue=0.2)]
        if flux_means is not None:
            rows.append(member(layer="mfa", role="flux", member_id="G6Pase_c", feature_id="G6Pase_c", tracer="pooled",
                               mean_numerator=flux_means[0], mean_denominator=flux_means[1]))
        return pd.DataFrame(rows)

    def test_patterns(self) -> None:
        table = pw.crossover_tests(self.PANEL, self.members(), ["A"]).set_index("name")
        step = table.loc["G6Pase_c"]
        self.assertEqual(step["direction"], "restricted")                 # 上游升、下游降
        self.assertAlmostEqual(step["pvalue"], 2 * max(0.001, 0.002))     # 两个单侧 p 的最大值 × 2
        self.assertAlmostEqual(step["statistic"], 1.7)
        self.assertIn("assumed_direction", step["flags"])
        self.assertEqual(table.loc["PYR_trans", "direction"], "not_tested")
        self.assertIn("同一个代谢物", table.loc["PYR_trans", "detail"])
        self.assertEqual(table.loc["PYK_c", "direction"], "not_tested")   # 下游没有测得
        # 通量为负（与书写方向相反）：上下游互换，模式随之互换
        reversed_ = pw.crossover_tests(self.PANEL, self.members((-5.0, -3.0)), ["A"]).set_index("name").loc["G6Pase_c"]
        self.assertEqual(reversed_["direction"], "activated")
        self.assertIn("reversed", reversed_["flags"])
        self.assertAlmostEqual(reversed_["pvalue"], step["pvalue"])


class ModuleDeclarationTest(unittest.TestCase):
    def test_errors(self) -> None:
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            for message, module in (("同时在 plus 与 minus 中", {"plus": ["Pcx"], "minus": ["Pcx"]}),
                                    ("没有成员", {"plus": [], "minus": []})):
                doc = {"panels": {"p": {"segments": [{"name": "s", "reactions": ["PC_m"]}],
                                        "signed_modules": {"m": module}}}}
                with self.subTest(message=message), self.assertRaisesRegex(PanelError, message):
                    load_panels(write_panels(Path(tmp), doc))


class PathwayBundleTest(unittest.TestCase):
    """在合成结果包上运行：声明分支节点与模块时写出三张表、图 D 与导出表。"""
    _bundle = panel_tests.DossierPanelBundleTest
    setUp, tearDown, source_bundle = _bundle.setUp, _bundle.tearDown, _bundle.source_bundle
    set_integration, integrate, run_dossier = _bundle.set_integration, _bundle.integrate, _bundle.run_dossier

    def test_tables_figure_and_notes(self) -> None:
        self.source_bundle()
        doc = {"panels": {"g6p": {**G6P_PANEL["panels"]["g6p"], "contrasts": [DIET]}}}
        write_panels(self.root / "experiments" / EXP, doc)
        base = self.run_dossier(panels="panels.yaml", figures={"formats": ["png"]}, pathway={"rotations": 99})
        splits = pd.read_parquet(base / "tables" / "node_split.parquet")
        self.assertEqual(list(splits.columns), pw.NODE_SPLIT_COLUMNS)
        self.assertEqual(set(splits["node"]), {"GLC6P_c"})
        self.assertEqual(set(splits["contrast_id"]), {DIET})
        # 合成结果包没有肝脏的逐只动物数据集：各层如实记为未测
        self.assertTrue(splits["state"].isin(["not_measured", "not_detected", "not_identifiable"]).all(), splits)
        tests = pd.read_parquet(base / "tables" / "pathway_tests.parquet")
        crossover = tests[tests["test"] == "crossover"].set_index("name")
        self.assertEqual(set(crossover.index), {"HEX_c", "G6Pase_c", "PGI_c"})
        folder = base / "figures" / "g6p" / "liver"
        for name in (f"node_split_{DIET}.png", f"node_split_{DIET}.tsv", f"pathway_tests_{DIET}.tsv"):
            self.assertTrue((folder / name).is_file(), name)
        self.assertTrue(any(n.startswith("节点分流比：") for n in self.notes), self.notes)

    def test_unknown_pathway_key(self) -> None:
        self.source_bundle()
        with self.assertRaisesRegex(Exception, "pathway 中有无法识别的项：rotation"):
            self.run_dossier(pathway={"rotation": 10})


if __name__ == "__main__":
    unittest.main()
