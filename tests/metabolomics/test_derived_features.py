#!/usr/bin/env python3
"""代谢组派生特征：脂质类别总量与质量作用比（映射修复阶段 2）。"""

from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from omicsweave.contract.dataset import ValueType, read_dataset
from omicsweave.contract.tables import read_differential
from omicsweave.core.run import run_experiment
from omicsweave.metabolomics.data_preprocessing import apply_fdr_correction
from omicsweave.metabolomics.derived_features import (
    LIPID_CLASS_TOTAL,
    MASS_ACTION_RATIO,
    feature_database_ids,
    lipid_class,
    lipid_class_totals,
    mass_action_ratios,
)

EXP = "derived_exp"


class LipidClassTest(unittest.TestCase):
    def test_class_from_lipid_maps_shorthand(self) -> None:
        cases = {
            "TG(16:0_18:1_18:2)": "TG", "PE(O-16:0_22:6)": "PE-O", "PE(P-18:0/20:4)": "PE-P",
            "Cer(d18:1/24:0)": "Cer", "FFA(16:0)": "FFA", "LPC(16:0/0:0)": "LPC",
            "Cyclo(Pro-Leu)": None, "L-Alanine": None, "Glycerol 3-Phosphate": None, None: None,
        }
        for name, expected in cases.items():
            with self.subTest(name=name):
                self.assertEqual(lipid_class(name), expected)

    def test_totals_are_summed_on_the_linear_scale(self) -> None:
        values = pd.DataFrame({"s1": [1.0, 1.0, 2.0, 5.0, 5.0], "s2": [2.0, 2.0, 3.0, 5.0, 5.0]},
                              index=["t1", "t2", "t3", "p1", "p2"])
        names = pd.Series(["TG(16:0_18:1_18:2)", "TG(16:0_16:0_18:1)", "TG(18:1_18:1_18:1)",
                           "PC(16:0_18:1)", "PC(18:0_20:4)"], index=values.index)

        derived = lipid_class_totals(values, names, min_species=3)

        self.assertEqual(list(derived.matrix.index), ["class_total:TG"])     # PC 只有 2 个分子种
        np.testing.assert_allclose(derived.matrix.loc["class_total:TG"], [3.0, 4.0])   # log2(2+2+4), log2(4+4+8)
        self.assertEqual(derived.features.loc["class_total:TG", "n_members"], 3)
        self.assertEqual(derived.features.loc["class_total:TG", "feature_type"], LIPID_CLASS_TOTAL)


class MassActionRatioTest(unittest.TestCase):
    def test_vendor_ids_are_normalized_and_names_are_a_fallback(self) -> None:
        table = pd.DataFrame({"name": ["Lactic acid", "Pyruvic Acid", "Glycerol 3-phosphate"],
                              "HMDB": ["HMDB00190", "-", ""], "cpd_ID": ["-", "C00022", ""]},
                             index=["a", "b", "c"])

        ids = feature_database_ids(table, "name")

        self.assertEqual(ids["a"], frozenset({"HMDB0000190"}))
        self.assertEqual(ids["b"], frozenset({"C00022"}))
        self.assertIn("HMDB0000126", ids["c"])            # 没有厂商 ID：按名称查

    def test_ratio_uses_the_most_intense_feature_and_reports_missing_components(self) -> None:
        values = pd.DataFrame({"s1": [10.0, 6.0, 8.0], "s2": [11.0, 6.5, 8.0]}, index=["lac", "lac_minor", "pyr"])
        ids = pd.Series({"lac": frozenset({"C00186"}), "lac_minor": frozenset({"HMDB0000190"}),
                         "pyr": frozenset({"C00022"})})
        names = pd.Series({"lac": "Lactic acid", "lac_minor": "L-Lactate", "pyr": "Pyruvic acid"})

        derived = mass_action_ratios(values, ids, names)

        np.testing.assert_allclose(derived.matrix.loc["ratio:lactate_pyruvate"], [2.0, 3.0])
        self.assertEqual(derived.features.loc["ratio:lactate_pyruvate", "related_reaction"], "LDH_c")
        self.assertIn("L-Lactate (lac_minor)", derived.features.loc["ratio:lactate_pyruvate", "definition"])
        self.assertEqual(derived.features["feature_type"].unique().tolist(), [MASS_ACTION_RATIO])
        self.assertTrue(any("acetoacetate" in note for note in derived.notes))


def _peak_table(columns: list[str], hfd: set[str]) -> pd.DataFrame:
    """6 种 TG 在 HFD 中升高 4 倍、乳酸升高 2 倍；3 种 PC、丙酮酸、G3P、DHAP 与 30 个填充代谢物不变。"""
    rng = np.random.default_rng(11)
    rows = ([(f"TG{i}", f"TG(16:0_18:{i}_18:2)", "", "", 4.0) for i in range(6)]
            + [(f"PC{i}", f"PC(16:0_18:{i})", "", "", 1.0) for i in range(3)]
            + [("LAC", "Lactic acid", "HMDB0000190", "C00186", 2.0), ("PYR", "Pyruvic Acid", "HMDB0000243", "C00022", 1.0),
               ("G3P", "Glycerol 3-Phosphate", "HMDB0000126", "C00093", 1.0),
               ("DHAP", "Dihydroxyacetone phosphate", "HMDB0001473", "C00111", 1.0)]
            + [(f"F{i}", f"metabolite_{i}", "", "", 1.0) for i in range(30)])
    table = pd.DataFrame(rows, columns=["peak_id", "name", "HMDB", "cpd_ID", "effect"])
    base = rng.lognormal(12, 0.8, len(table))
    for column in columns:
        effect = table["effect"].to_numpy() if column in hfd else 1.0
        table[column] = base * effect * rng.lognormal(0, 0.05, len(table))
    return table.drop(columns="effect")


class RunnerDerivedFeaturesTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        root = Path(cls._tmp.name).resolve()
        exp_dir = root / "experiments" / EXP
        exp_dir.mkdir(parents=True)
        animals = {f"m{i}": ("LFD" if i <= 4 else "HFD") for i in range(1, 9)}
        hfd = {m for m, g in animals.items() if g == "HFD"}
        _peak_table(list(animals), hfd).to_csv(exp_dir / "peaks.tsv", sep="\t", index=False)
        (root / "project.yaml").write_text("schema_version: 1\nproject: {id: p, species: mus_musculus}\n",
                                           encoding="utf-8")
        (exp_dir / "experiment.yaml").write_text(yaml.safe_dump({
            "schema_version": 1,
            "experiment": {"id": EXP},
            "design": {"factors": {"diet": {"level": "subject", "levels": ["LFD", "HFD"]}},
                       "contrasts": [{"id": "HFD_vs_LFD", "factor": "diet", "numerator": "HFD",
                                      "denominator": "LFD", "role": "disease_effect"}]},
            "omics": {"metabolomics": {"enabled": True, "input": {
                "level": "table", "format": "peak_table", "files": [f"experiments/{EXP}/peaks.tsv"],
                "feature_id_column": "peak_id", "feature_name_column": "name"}}},
        }, sort_keys=False), encoding="utf-8")
        pd.DataFrame({"subject_id": list(animals), "diet": list(animals.values())}).to_csv(
            exp_dir / "subjects.tsv", sep="\t", index=False)
        pd.DataFrame([{"sample_id": m, "subject_id": m, "omics": "metabolomics", "tissue": "liver"}
                      for m in animals]).to_csv(exp_dir / "samples.tsv", sep="\t", index=False)
        with contextlib.redirect_stdout(io.StringIO()):
            run_dir, cls.notes = run_experiment(EXP, ["metabolomics"], steps=["differential"], root=root)
        cls.met = run_dir / "omics" / "metabolomics"
        cls.derived = read_differential(cls.met / "tables" / "derived_differential.parquet").set_index("feature_id")
        cls.measured = read_differential(cls.met / "tables" / "differential.parquet")

    @classmethod
    def tearDownClass(cls) -> None:
        cls._tmp.cleanup()

    def test_derived_dataset(self) -> None:
        dataset = read_dataset(self.met / "derived" / "liver")

        self.assertEqual(dataset.meta.value_type, ValueType.DERIVED)
        self.assertEqual(set(dataset.matrix.index), {"class_total:TG", "class_total:PC", "ratio:lactate_pyruvate",
                                                     "ratio:g3p_dhap"})

    def test_effects_follow_the_simulated_changes(self) -> None:
        self.assertAlmostEqual(self.derived.at["class_total:TG", "log2fc"], 2.0, delta=0.3)
        self.assertAlmostEqual(self.derived.at["class_total:PC", "log2fc"], 0.0, delta=0.3)
        self.assertAlmostEqual(self.derived.at["ratio:lactate_pyruvate", "log2fc"], 1.0, delta=0.3)
        self.assertTrue(self.derived.at["class_total:TG", "selected"])
        self.assertEqual(self.derived.at["class_total:TG", "units"], "log2 summed normalized peak area")
        self.assertEqual(self.derived.at["ratio:g3p_dhap", "related_reaction"], "GPS_c")

    def test_each_feature_type_is_corrected_on_its_own(self) -> None:
        for feature_type, rows in self.derived.groupby("feature_type"):
            with self.subTest(feature_type=feature_type):
                np.testing.assert_allclose(rows["padj"], apply_fdr_correction(rows["pvalue"].to_numpy(), "bh"))

    def test_measured_table_is_unchanged_and_notes_name_missing_ratios(self) -> None:
        self.assertFalse(self.measured["feature_id"].str.contains(":").any())
        self.assertEqual(len(self.measured), 43)
        self.assertTrue(any("acetoacetate" in note for note in self.notes))


if __name__ == "__main__":
    unittest.main()
