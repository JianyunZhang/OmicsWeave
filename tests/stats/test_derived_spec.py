#!/usr/bin/env python3
"""声明式派生特征（G8）与脂质简写解析（G7）。

验收（设计文档 G8）：由网络自动生成的 step_ratio 与手写定义逐值一致；成员缺失时该特征不计算并说明原因。
"""

from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from omicsweave.contract.dataset import Dataset, DatasetMeta, ValueType, read_dataset, write_dataset
from omicsweave.contract.tables import read_differential
from omicsweave.core.idmap.lipids import parse_lipid
from omicsweave.core.model import Omics
from omicsweave.metabolomics.derived_features import feature_database_ids, mass_action_ratios
from omicsweave.stats.post import postprocess
from tests.stats.test_models_robustness import StatsCase, subjects

FEATURES = [  # ID、名称、HMDB、KEGG
    ("M01", "Lactic acid", "HMDB0000190", "C00186"), ("M02", "Pyruvic acid", "HMDB0000243", "C00022"),
    ("M03", "TG(16:0_18:1_18:2)", "", ""), ("M04", "TG(18:0_18:1_18:1)", "", ""), ("M05", "TG(52:3)", "", ""),
    ("M06", "PE(16:0/20:4)", "", ""), ("M07", "PE(O-16:0/22:6)", "", ""), ("M08", "PE(P-18:0/20:4)", "", ""),
    ("M09", "PE(18:0/22:6)", "", ""), ("M10", "Glucose", "HMDB0000122", "C00031"),
]
SPECS = [
    {"type": "step_ratio", "omics": "metabolomics", "source": "network"},
    {"id": "lac_pyr", "type": "ratio", "omics": "metabolomics", "numerator": {"database_ids": ["HMDB0000190"]},
     "denominator": {"database_ids": ["HMDB0000243"]}},
    {"id": "pe_ether", "type": "fraction", "omics": "metabolomics", "members": {"lipid_class": ["PE-O", "PE-P"]},
     "denominator": {"lipid_class": ["PE", "PE-O", "PE-P"]}},
    {"id": "tg_chain", "type": "index", "omics": "metabolomics", "members": {"lipid_class": ["TG"]},
     "index": "chain_length"},
    {"id": "tg_db", "type": "index", "omics": "metabolomics", "members": {"lipid_class": ["TG"]},
     "index": "double_bonds"},
    {"id": "cardiolipin", "type": "sum", "omics": "metabolomics", "members": {"lipid_class": ["CL"]}},
]


class DerivedSpecTest(StatsCase):
    def setUp(self) -> None:
        super().setUp()
        exp_dir = self.root / "experiments" / self.experiment().id
        samples = pd.read_csv(exp_dir / "samples.tsv", sep="\t", dtype=str)
        mfa = samples[samples["omics"] == "mfa"]
        met = mfa.assign(sample_id="met_" + mfa["subject_id"], omics="metabolomics", tissue="liver")
        pd.concat([samples, met]).to_csv(exp_dir / "samples.tsv", sep="\t", index=False)
        rng = np.random.default_rng(11)
        meta = subjects()
        values = pd.DataFrame(rng.normal(18, 1, (len(FEATURES), len(meta))), index=[f[0] for f in FEATURES],
                              columns=[f"met_{s}" for s in meta.index])
        values.loc["M01", [f"met_{s}" for s in meta.index[meta["diet"] == "HFD"]]] += 1.5   # HFD 乳酸升高
        self.values = values
        self.annotations = pd.DataFrame({"Compounds": [f[1] for f in FEATURES], "HMDB": [f[2] for f in FEATURES],
                                         "cpd_ID": [f[3] for f in FEATURES]}, index=values.index)
        write_dataset(self.run_dir / "omics" / "metabolomics" / "data" / "liver", Dataset(
            values, self.annotations, DatasetMeta(omics=Omics.METABOLOMICS, tissue="liver",
                                                  value_type=ValueType.PEAK_AREA, scale="log2", feature_id_type="Index")))
        self.set_design()
        path = exp_dir / "experiment.yaml"
        import yaml
        doc = yaml.safe_load(path.read_text(encoding="utf-8"))
        doc["derived"] = SPECS
        path.write_text(yaml.safe_dump(doc, sort_keys=False), encoding="utf-8")
        self.notes = postprocess(self.experiment(), self.run_dir)
        self.derived = read_dataset(self.run_dir / "omics" / "metabolomics" / "derived" / "liver")

    def test_network_step_ratio_matches_handwritten_definition(self) -> None:
        ids = feature_database_ids(self.annotations.assign(_n=self.annotations["Compounds"]), "_n")
        handwritten = mass_action_ratios(self.values, ids, self.annotations["Compounds"])
        expected = handwritten.matrix.loc["ratio:lactate_pyruvate"]
        np.testing.assert_allclose(self.derived.matrix.loc["step_ratio:LDH_c"], expected)
        np.testing.assert_allclose(self.derived.matrix.loc["ratio:lac_pyr"], expected)
        row = self.derived.features.loc["step_ratio:LDH_c"]
        self.assertEqual((row["related_reaction"], row["members"]), ("LDH_c", "M01;M02"))
        self.assertTrue(row["generated_from"].startswith("network:"))

    def test_fraction_and_index(self) -> None:
        linear = np.power(2.0, self.values)
        ether = np.log2(linear.loc[["M07", "M08"]].sum() / linear.loc[["M06", "M07", "M08", "M09"]].sum())
        np.testing.assert_allclose(self.derived.matrix.loc["fraction:pe_ether"], ether)
        tg = ["M03", "M04", "M05"]
        length = pd.Series({m: parse_lipid(n).mean_chain_length for m, n in zip(self.values.index,
                                                                                self.annotations["Compounds"])
                            if m in tg})
        expected = np.log2((linear.loc[tg].mul(length, axis=0)).sum() / linear.loc[tg].sum())
        np.testing.assert_allclose(self.derived.matrix.loc["index:tg_chain"], expected)
        self.assertAlmostEqual(length["M05"], 52 / 3)

    def test_missing_members_are_reported_not_computed(self) -> None:
        self.assertNotIn("sum:cardiolipin", self.derived.matrix.index)
        self.assertTrue(any("cardiolipin" in n and "未计算" in n for n in self.notes), self.notes)

    def test_derived_differential_and_definitions(self) -> None:
        table = read_differential(self.run_dir / "omics" / "metabolomics" / "tables" / "derived_differential.parquet")
        lactate = table[(table["feature_id"] == "step_ratio:LDH_c") & (table["contrast_id"] == "HFD_vs_LFD")].iloc[0]
        self.assertGreater(lactate["log2fc"], 0.8)
        self.assertEqual((lactate["feature_type"], lactate["estimand_type"]), ("step_ratio", "pairwise"))
        self.assertEqual(set(self.derived.features.columns) >= {"definition", "generated_from"}, True)
        # 每个派生类型在“组织 × 比较”内单独校正：只有一个特征的类型，padj 等于 p 值
        one = table[(table["feature_type"] == "fraction") & (table["contrast_id"] == "HFD_vs_LFD")]
        self.assertAlmostEqual(float(one["padj"].iloc[0]), float(one["pvalue"].iloc[0]))


class LipidParserTest(unittest.TestCase):
    def test_shorthand_variants(self) -> None:
        cases = {"TG(16:0_18:1_18:2)": ("TG", 52, 3, "molecular_species"), "TG(52:3)": ("TG", 52, 3, "species"),
                 "PC 34:1": ("PC", 34, 1, "species"), "PE(O-16:0/22:6)": ("PE-O", 38, 6, "molecular_species"),
                 "PE(P-18:0_20:4)": ("PE-P", 38, 4, "molecular_species"), "LPC(16:0)": ("LPC", 16, 0, "molecular_species"),
                 "Cer(d18:1/24:0(2OH))": ("Cer", 42, 1, "molecular_species"), "Cer 18:1;O2/24:0": ("Cer", 42, 1,
                                                                                                   "molecular_species"),
                 "MG(20:4n3)": ("MG", 20, 4, "molecular_species"), "PG(16:0/0:0)[U]": ("PG", 16, 0, "molecular_species")}
        for name, (cls, carbons, double_bonds, level) in cases.items():
            with self.subTest(name=name):
                parsed = parse_lipid(name)
                self.assertEqual((parsed.lipid_class, parsed.carbons, parsed.double_bonds, parsed.level),
                                 (cls, carbons, double_bonds, level))
        self.assertIsNone(parse_lipid("Glucose"))
        self.assertIsNone(parse_lipid("Citric acid"))


if __name__ == "__main__":
    unittest.main()
