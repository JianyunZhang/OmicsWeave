#!/usr/bin/env python3
"""反应档案 S2（G13 第二部分与 G14）：通路面板、逐只动物表、其他来源、矩阵视图与导出表、图。

验收（反应档案与分轴综合评估，第 5 节 G14 与 4.6）：
- 声明中引用的反应不在网络中时报错；
- 面板声明了示踪剂而数据中没有时，矩阵上出现 not_infused 的列；
- 图与导出表共用同一个视图：同一个（反应, 层, 成员）在导出表中的数值与 dossier_members 相同。
"""

from __future__ import annotations

import shutil
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from omicsweave.contract.dataset import Dataset, DatasetMeta, ValueType, write_dataset
from omicsweave.contract.manifest import RunManifest, write_manifest
from omicsweave.contract.tables import write_differential
from omicsweave.core.model import Omics
from omicsweave.integration.dossier.panels import PanelError, load_panels, panel_elements
from omicsweave.integration.dossier.view import matrix_view, summarize, wide_table
from tests.integration.test_flux_anchored import DIET, FluxAnchoredBundleTest
from tests.mfa.test_mfa_runner import EXP, build_two_factor_project

PANEL = {"panels": {"liver_glucose": {
    "title": "toy liver", "tissues": ["liver"], "contrasts": [DIET],
    "segments": [{"name": "uptake", "reactions": ["HEX_c"]}, {"name": "release", "reactions": ["G6Pase_c"]}],
    "expected_tracers": {"G6Pase_c": ["glucose"]}}}}


def write_bundle(root: Path, tables: dict[str, pd.DataFrame], omics: tuple[str, ...]) -> Path:
    """在另一个项目目录中写出源结果包（统一差异结果表与运行清单）。"""
    run_dir = root / "analysis_results" / EXP / "run_20261101_000000_src"
    for name in omics:
        write_differential(tables[name], run_dir / "omics" / name / "tables" / "differential.parquet")
    write_manifest(RunManifest(run_id=run_dir.name, project_id="other", experiment_id=EXP,
                               omics=[Omics(o) for o in omics]), run_dir)
    return run_dir


def write_panels(directory: Path, doc: dict) -> Path:
    path = directory / "panels.yaml"
    path.write_text(yaml.safe_dump(doc, allow_unicode=True, sort_keys=False), encoding="utf-8")
    return path


class PanelTest(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.dir = Path(self._tmp.name)

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def test_elements_share_nodes_between_consecutive_reactions(self) -> None:
        panel = load_panels(write_panels(self.dir, {"panels": {"p": {
            "segments": [{"name": "s", "reactions": ["PYR_trans", "PC_m", "PDH_m"]}],
            "branch_nodes": [{"node": "PYR_m", "exits": {"carboxylation": ["PC_m"], "oxidation": ["PDH_m"]}}]}}}))["p"]
        elements = panel_elements(panel)
        self.assertEqual(list(elements["element_id"]),
                         ["PYR_c", "PYR_trans", "PYR_m", "PC_m", "OAC_m", "PYR_m", "PDH_m", "ACCOA_m"])
        self.assertNotIn("CO2", set(elements["element_id"]))             # 辅因子不作为节点
        self.assertTrue(elements.loc[elements["element_id"] == "PYR_m", "branch_node"].all())
        exits = elements.set_index("element_id")["exit_label"]
        self.assertEqual((exits["PC_m"], exits["PDH_m"]), ("carboxylation", "oxidation"))
        self.assertEqual(elements.set_index("element_id").at["OAC_m", "label"], "oxaloacetate (m)")

    def test_declaration_errors(self) -> None:
        segment = [{"name": "s", "reactions": ["PC_m"]}]
        cases = {
            "不在 MFA 网络中": {"segments": [{"name": "s", "reactions": ["NOT_A_REACTION"]}]},
            "不止一次": {"segments": [{"name": "s", "reactions": ["PC_m", "PC_m"]}]},
            "expected_tracers": {"segments": segment, "expected_tracers": {"PDH_m": ["glucose"]}},
            "无法识别的键": {"segments": segment, "colour": "red"},
            "分支节点": {"segments": segment, "branch_nodes": [{"node": "GLC_c", "exits": {"x": ["PC_m"]}}]},
            "没有 segments": {"title": "empty"},
        }
        for message, spec in cases.items():
            with self.subTest(message=message), self.assertRaisesRegex(PanelError, message):
                load_panels(write_panels(self.dir, {"panels": {"p": spec}}))


def _member(**fields) -> dict:
    base = {"source": "A", "layer": "proteomics", "role": "direct_enzyme", "member_id": "X", "feature_id": "X",
            "feature_name": "X", "member_relation": "isozyme", "regulator_effect": "", "regulator_mode": "",
            "log2fc": np.nan, "se": np.nan, "pvalue": np.nan, "padj": np.nan, "lfsr": np.nan,
            "abundance_share": np.nan, "cell_state": "no_evidence", "selected": False}
    return {**base, **fields}


class SummaryTest(unittest.TestCase):
    def test_states_and_values(self) -> None:
        self.assertEqual(summarize(pd.DataFrame())["display_state"], "no_knowledge")
        self.assertEqual(summarize(pd.DataFrame([_member(cell_state="not_measured", feature_id="")]))["display_state"],
                         "not_measured")
        missing = pd.DataFrame([_member(cell_state="not_detected", feature_id=""),
                                _member(cell_state="not_measured", feature_id="")])
        self.assertEqual(summarize(missing)["display_state"], "not_detected")
        single = summarize(pd.DataFrame([_member(cell_state="changed", selected=True, log2fc=1.0, padj=0.01,
                                                 pvalue=0.001)]))
        self.assertEqual((single["display_state"], single["ring"], single["value"]), ("changed", "selected", 1.0))
        # 同工酶按丰度占比加权：0.9 × 1.0 + 0.1 × (−1.0) = 0.8
        isozymes = pd.DataFrame([_member(member_id="A", log2fc=1.0, abundance_share=0.9, cell_state="changed"),
                                 _member(member_id="B", log2fc=-1.0, abundance_share=0.1)])
        from omicsweave.integration.dossier.view import _enzyme_value
        cell = summarize(isozymes, _enzyme_value)
        self.assertAlmostEqual(cell["value"], 0.8)
        self.assertEqual((cell["display_state"], cell["n_changed"], cell["n_quantified"], cell["ring"]),
                         ("summary", 1, 2, "changed"))


class DossierPanelBundleTest(unittest.TestCase):
    setUp = FluxAnchoredBundleTest.setUp
    tearDown = FluxAnchoredBundleTest.tearDown
    source_bundle = FluxAnchoredBundleTest.source_bundle
    set_integration = FluxAnchoredBundleTest.set_integration
    integrate = FluxAnchoredBundleTest.integrate

    def add_heart_rna(self, root: Path, run_dir: Path) -> None:
        """心脏转录组：样本表加上 12 只动物的 RNA 样本，结果包中写入计数数据集。"""
        samples_path = root / "experiments" / EXP / "samples.tsv"
        samples = pd.read_csv(samples_path, sep="\t")
        subjects = samples["subject_id"].tolist()
        rna = pd.DataFrame({"sample_id": [f"rna_{s}" for s in subjects], "subject_id": subjects,
                            "omics": "transcriptomics", "tissue": "heart"})
        pd.concat([samples, rna], ignore_index=True).to_csv(samples_path, sep="\t", index=False)
        genes = ["Pdha1", "Dlat", "Actb"]
        counts = pd.DataFrame(np.arange(1, 1 + 3 * len(subjects)).reshape(3, -1) * 10.0, index=genes,
                              columns=rna["sample_id"])
        write_dataset(run_dir / "omics" / "transcriptomics" / "data" / "heart", Dataset(
            counts, pd.DataFrame(index=counts.index), DatasetMeta(omics=Omics.TRANSCRIPTOMICS, tissue="heart",
                                                                  value_type=ValueType.COUNTS, scale="linear",
                                                                  feature_id_type="gene_symbol")))

    def run_dossier(self, **parameters) -> Path:
        self.set_integration(enabled=False, modes=["reaction_dossier"], parameters={"reaction_dossier": parameters})
        run_dir, self.notes = self.integrate()
        return run_dir / "integration" / "reaction_dossier"

    def test_samples_and_contrast_sides(self) -> None:
        source = self.source_bundle()
        self.add_heart_rna(self.root, source)
        base = self.run_dossier(reactions=["PDH_m"])
        samples = pd.read_parquet(base / "tables" / "dossier_samples.parquet")
        heart = samples[(samples["layer"] == "transcriptomics") & (samples["tissue"] == "heart")]
        self.assertEqual(set(heart["feature_id"]), {"Pdha1", "Dlat"})           # 只取成员表引用到的特征
        self.assertEqual(len(heart), 2 * 12)
        self.assertEqual(heart["group"].nunique(), 4)                          # diet · treatment 都有变化
        order = heart.drop_duplicates("group").sort_values("group_order")["group"].tolist()
        self.assertEqual(order, ["LFD·Saline", "LFD·GLP1", "HFD·Saline", "HFD·GLP1"])
        self.assertTrue((heart["scale"] == "log2").all())
        sides = pd.read_parquet(base / "tables" / "dossier_sample_contrasts.parquet")
        diet = sides[(sides["layer"] == "transcriptomics") & (sides["contrast_id"] == DIET)]
        self.assertEqual(dict(diet.groupby("side")["level"].first()), {"denominator": "LFD", "numerator": "HFD"})
        self.assertEqual(len(diet), 6)                                         # subset：只含 Saline 的动物

    def test_other_source(self) -> None:
        self.source_bundle()
        other = self.root.parent / f"{self.root.name}_B"
        other.mkdir()
        build_two_factor_project(other)
        try:
            write_bundle(other, self.tables, ("transcriptomics", "mfa"))
            base = self.run_dossier(reactions=["HEX_c", "PDH_m"], tracers=["glucose"], label="A", sources=[{
                "label": "B", "project": str(other), "experiment": EXP, "layers": ["transcriptomics", "mfa"],
                "tissue_map": {"heart": "cardiac"}, "contrast_map": {DIET: DIET}}])
            members = pd.read_parquet(base / "tables" / "dossier_members.parquet")
        finally:
            shutil.rmtree(other, ignore_errors=True)
        b = members[members["source"] == "B"]
        self.assertEqual(set(b["layer"]), {"transcriptomics", "mfa"})
        self.assertEqual(set(b["contrast_id"]), {DIET})
        self.assertIn("cardiac", set(b["tissue"]))                             # 来源的组织名换成共同名称
        self.assertNotIn("heart", set(b["tissue"]))
        # 声明的示踪剂只用于本实验
        self.assertTrue((members[(members["source"] == "A") & (members["tracer"] == "glucose")]["cell_state"]
                         == "not_infused").all())
        self.assertFalse((b["tracer"] == "glucose").any())
        same = members[(members["flux_id"] == "HEX_c") & (members["tissue"] == "liver")
                       & (members["contrast_id"] == DIET) & (members["layer"] == "transcriptomics")
                       & (members["member_id"] == "Hk2")]
        self.assertEqual(sorted(same["source"]), ["A", "B"])
        self.assertEqual(same["log2fc"].nunique(), 1)                          # 同一份合成数据

    def test_panel_figures_and_exports(self) -> None:
        self.source_bundle()
        write_panels(self.root / "experiments" / EXP, PANEL)
        base = self.run_dossier(panels="panels.yaml", tracers=["glucose", "alanine"], label="A",
                                figures={"formats": ["png"]})
        folder = base / "figures" / "liver_glucose" / "liver"
        for name in (f"matrix_{DIET}.png", f"matrix_{DIET}.tsv", f"matrix_{DIET}_members.tsv",
                     "dossier_HEX_c_HFD_vs_LFD.png", "dossier_G6Pase_c_HFD_vs_LFD.png", f"tracer_coverage_{DIET}.png"):
            self.assertTrue((folder / name).is_file(), name)
        elements = pd.read_parquet(base / "tables" / "panel_elements.parquet")
        self.assertEqual(list(elements.loc[elements["element_type"] == "reaction", "element_id"]),
                         ["HEX_c", "G6Pase_c"])
        members = pd.read_parquet(base / "tables" / "dossier_members.parquet")
        panel = load_panels(self.root / "experiments" / EXP / "panels.yaml")["liver_glucose"]
        view = matrix_view(panel, panel_elements(panel), members, "liver", DIET, ["A"], ["glucose", "alanine"])
        # 声明的示踪剂没有数据：矩阵上有 not_infused 的列，预期示踪剂标出
        flux = view[view["block"] == "flux"].set_index(["element_id", "tracer"])
        self.assertEqual(flux.at[("G6Pase_c", "alanine"), "display_state"], "not_infused")
        self.assertTrue(flux.at[("G6Pase_c", "glucose"), "expected_tracer"])
        # 图与导出表共用视图；导出的成员行与 dossier_members 逐值相同
        exported = pd.read_csv(folder / f"matrix_{DIET}.tsv", sep="\t")
        self.assertEqual(len(exported), len(elements))
        self.assertEqual(wide_table(view).shape, exported.shape)
        cell = view[(view["element_id"] == "HEX_c") & (view["column"] == "R:A")].iloc[0]
        self.assertAlmostEqual(exported.loc[exported["element_id"] == "HEX_c", "enzyme R:A log2fc"].iloc[0],
                               round(cell["value"], 4))
        keys = ["source", "flux_id", "layer", "role", "member_id", "feature_id", "tracer"]
        exported_members = pd.read_csv(folder / f"matrix_{DIET}_members.tsv", sep="\t", keep_default_na=False)
        source = members[(members["tissue"] == "liver") & (members["contrast_id"] == DIET)
                         & members["flux_id"].isin(panel.reactions)].copy()
        source[keys] = source[keys].fillna("").astype(str)
        exported_members[keys] = exported_members[keys].astype(str)
        merged = exported_members.merge(source, on=keys, suffixes=("_x", ""))
        self.assertEqual(len(merged), len(exported_members))
        self.assertEqual(len(merged), len(source))
        np.testing.assert_allclose(pd.to_numeric(merged["log2fc_x"], errors="coerce"),
                                   merged["log2fc"].astype(float), equal_nan=True)
        self.assertEqual(list(merged["cell_state_x"]), list(merged["cell_state"]))
        self.assertTrue(any("面板 liver_glucose" in n for n in self.notes))

    def test_panel_errors_stop_the_run(self) -> None:
        self.source_bundle()
        write_panels(self.root / "experiments" / EXP,
                     {"panels": {"p": {"segments": [{"name": "s", "reactions": ["NOT_A_REACTION"]}]}}})
        with self.assertRaisesRegex(Exception, "不在 MFA 网络中"):
            self.run_dossier(panels="panels.yaml")


if __name__ == "__main__":
    unittest.main()
