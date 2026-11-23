#!/usr/bin/env python3
"""设计感知的统计（G3）与留一稳健性（E5）。

验收（设计文档 G3）：两因素合成数据上，合并主效应的估计与真值一致，功效高于单臂比较；设计矩阵不满秩时拒绝拟合并说明原因。
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from omicsweave.contract.dataset import Dataset, DatasetMeta, ValueType, empty_features, write_dataset
from omicsweave.core.config import load_experiment, load_project
from omicsweave.core.model import ExperimentConfig, Omics
from omicsweave.core.statistics.linear_model import bh, fit_contrasts
from omicsweave.core.statistics.shrinkage import adaptive_shrinkage, retrodesign
from omicsweave.stats.models import fit_models
from omicsweave.stats.post import postprocess
from omicsweave.stats.robustness import robustness
from tests.mfa.test_mfa_runner import EXP, build_two_factor_project

N_FEATURES, N_DIET, N_TREAT = 300, 60, 60
DIET_EFFECT, TREAT_EFFECT, NOISE = 1.0, 0.8, 0.6


def subjects() -> pd.DataFrame:
    """build_two_factor_project 的动物：LFD-Saline m0–m2，LFD-GLP1 m3–m5，HFD-Saline m6–m8，HFD-GLP1 m9–m11。"""
    cells = [(d, t) for d in ("LFD", "HFD") for t in ("Saline", "GLP1")]
    return pd.DataFrame([(f"m{i * 3 + k}", d, t) for i, (d, t) in enumerate(cells) for k in range(3)],
                        columns=["subject_id", "diet", "treatment"]).set_index("subject_id")


def simulate(seed: int = 7) -> pd.DataFrame:
    """log2 强度：前 N_DIET 个特征有饮食效应，随后 N_TREAT 个有处理效应，其余没有效应。"""
    rng = np.random.default_rng(seed)
    meta = subjects()
    diet = (meta["diet"] == "HFD").to_numpy(dtype=float)
    treat = (meta["treatment"] == "GLP1").to_numpy(dtype=float)
    effects = np.zeros((N_FEATURES, len(meta)))
    effects[:N_DIET] += DIET_EFFECT * diet
    effects[N_DIET:N_DIET + N_TREAT] += TREAT_EFFECT * treat
    values = 20 + effects + rng.normal(0, NOISE, effects.shape)
    return pd.DataFrame(values, index=[f"P{i:03d}" for i in range(N_FEATURES)],
                        columns=[f"prot_{s}" for s in meta.index])


class StatsCase(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.root = Path(self._tmp.name).resolve()
        build_two_factor_project(self.root)
        exp_dir = self.root / "experiments" / EXP
        samples = pd.read_csv(exp_dir / "samples.tsv", sep="\t", dtype=str)
        prot = pd.DataFrame({"sample_id": "prot_" + samples["subject_id"], "subject_id": samples["subject_id"],
                             "omics": "proteomics", "tissue": "liver"})
        pd.concat([samples, prot]).to_csv(exp_dir / "samples.tsv", sep="\t", index=False)
        self.run_dir = self.root / "analysis_results" / EXP / "run_test"
        self.write_matrix(simulate())

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def write_matrix(self, matrix: pd.DataFrame) -> None:
        write_dataset(self.run_dir / "omics" / "proteomics" / "data" / "liver", Dataset(
            matrix, empty_features(matrix.index), DatasetMeta(omics=Omics.PROTEOMICS, tissue="liver",
                                                              value_type=ValueType.INTENSITY, scale="log2",
                                                              feature_id_type="protein")))

    def set_design(self, **design) -> None:
        path = self.root / "experiments" / EXP / "experiment.yaml"
        doc = yaml.safe_load(path.read_text(encoding="utf-8"))
        doc["design"].update(design)
        path.write_text(yaml.safe_dump(doc, sort_keys=False), encoding="utf-8")

    def experiment(self):
        return load_experiment(load_project(self.root), EXP)


class FormulaModelTest(StatsCase):
    def models(self, *specs) -> tuple[pd.DataFrame, pd.DataFrame]:
        self.set_design(models=list(specs))
        measured, _, fits, _ = fit_models(self.experiment(), self.run_dir, {}, 0.05)
        rows = pd.concat(measured.get("proteomics", [pd.DataFrame()]), ignore_index=True)
        return rows, fits

    def test_pooled_main_effect_matches_truth_and_has_more_power(self) -> None:
        rows, fits = self.models({"id": "pooled", "formula": "~ diet + treatment", "report": ["diet"]})
        self.assertEqual(set(fits["status"]), {"fit"})
        diet = rows[rows["contrast_id"] == "pooled.HFD_vs_LFD"].set_index("feature_id")
        self.assertEqual(set(diet["estimand_type"]), {"pooled_main"})
        self.assertEqual((diet["n_numerator"].iloc[0], diet["n_denominator"].iloc[0], diet["n_total"].iloc[0]),
                         (6, 6, 12))
        true = diet.index[:N_DIET]
        self.assertAlmostEqual(diet.loc[true, "log2fc"].mean(), DIET_EFFECT, delta=0.1)
        self.assertAlmostEqual(diet.loc[diet.index[N_DIET:], "log2fc"].mean(), 0.0, delta=0.05)

        # 单臂比较：只用 Saline 动物（每组 3 只）
        matrix = simulate()
        saline = [c for c in matrix.columns if c[-2:] in ("m0", "m1", "m2", "m6", "m7", "m8")]
        x = np.column_stack([np.ones(6), [0, 0, 0, 1, 1, 1]])
        single = fit_contrasts(matrix[saline], x, {"diet": np.array([0.0, 1.0])})["diet"].table
        pooled_hits = int((diet.loc[true, "padj"] < 0.05).sum())
        single_hits = int((bh(single["pvalue"])[:N_DIET] < 0.05).sum())
        self.assertGreater(pooled_hits, single_hits)

    def test_interaction_model_reports_averaged_main_effect(self) -> None:
        rows, _ = self.models({"id": "additive", "formula": "~ diet + treatment"},
                              {"id": "factorial", "formula": "~ diet * treatment"})
        estimate = rows.pivot_table(index="feature_id", columns="contrast_id", values="log2fc")
        # 平衡设计中，交互模型对另一因素等权平均的主效应等于加性模型的主效应
        np.testing.assert_allclose(estimate["factorial.HFD_vs_LFD"], estimate["additive.HFD_vs_LFD"], atol=1e-10)
        kinds = rows.drop_duplicates("contrast_id").set_index("contrast_id")["estimand_type"]
        self.assertEqual(kinds["factorial.diet[HFD]:treatment[GLP1]"], "interaction")
        self.assertAlmostEqual(estimate["factorial.diet[HFD]:treatment[GLP1]"].mean(), 0.0, delta=0.1)

    def test_rank_deficient_design_is_refused_with_reason(self) -> None:
        # 去掉 LFD-GLP1 的动物：处理只出现在 HFD 中，交互项不可估计
        matrix = simulate().drop(columns=["prot_m3", "prot_m4", "prot_m5"])
        self.write_matrix(matrix)
        rows, fits = self.models({"id": "factorial", "formula": "~ diet * treatment"})
        self.assertTrue(rows.empty)
        self.assertEqual(fits["status"].tolist(), ["not_fit"])
        self.assertIn("不满秩", fits["reason"].iloc[0])

    def test_mfa_models_adjust_for_tracer_and_keep_identifiability(self) -> None:
        exp_dir = self.root / "experiments" / EXP
        samples = pd.read_csv(exp_dir / "samples.tsv", sep="\t", dtype=str, keep_default_na=False)
        mfa = samples["omics"] == "mfa"
        samples.loc[mfa, "protocol"] = ["glucose" if i % 2 else "lactate" for i in range(int(mfa.sum()))]
        samples.loc[~mfa, "protocol"] = ""
        samples.to_csv(exp_dir / "samples.tsv", sep="\t", index=False)
        rng = np.random.default_rng(2)
        flux = pd.DataFrame(np.exp(rng.normal(1, 0.2, (2, 12))), index=["PDH_m", "CS_m"],
                            columns=[f"m{i}" for i in range(12)])
        write_dataset(self.run_dir / "omics" / "mfa" / "data" / "liver", Dataset(
            flux, empty_features(flux.index), DatasetMeta(omics=Omics.MFA, tissue="liver", value_type=ValueType.FLUX,
                                                          scale="linear", feature_id_type="mfa_reaction", units="nmol")))
        pairwise = pd.DataFrame({"contrast_id": "HFD_vs_LFD", "tissue": "liver", "feature_id": ["PDH_m", "CS_m"],
                                 "flags": ["not_identifiable", ""], "identifiability": ["poorly_constrained", "moderate"]})
        self.set_design(models=[{"id": "pooled", "formula": "~ diet", "omics": ["mfa"]}])
        measured, _, fits, _ = fit_models(self.experiment(), self.run_dir, {"mfa": pairwise}, 0.05)
        self.assertIn("protocol", fits["formula"].iloc[0])
        rows = measured["mfa"][0].set_index("feature_id")
        self.assertIn("not_identifiable", rows.at["PDH_m", "flags"])
        self.assertEqual((rows.at["PDH_m", "identifiability"], rows.at["CS_m", "identifiability"]),
                         ("poorly_constrained", "moderate"))
        self.assertTrue(rows["mean_numerator"].notna().all())

    def test_config_rejects_unknown_terms(self) -> None:
        doc = yaml.safe_load((self.root / "experiments" / EXP / "experiment.yaml").read_text(encoding="utf-8"))
        doc["design"]["models"] = [{"id": "bad", "formula": "~ diet + strain"}]
        with self.assertRaisesRegex(Exception, "strain"):
            ExperimentConfig.model_validate(doc)
        doc["design"]["models"] = [{"id": "bad", "formula": "~ diet", "report": ["treatment"]}]
        with self.assertRaisesRegex(Exception, "report"):
            ExperimentConfig.model_validate(doc)


class RobustnessTest(StatsCase):
    def test_single_animal_driven_feature_is_flagged(self) -> None:
        matrix = simulate()
        # 特征 P250（无真实效应）只在一只 HFD-Saline 动物上极高；P000 是稳定的饮食效应
        matrix.loc["P250", "prot_m6"] += 6.0
        self.write_matrix(matrix)
        diff = pd.DataFrame({
            "contrast_id": "HFD_vs_LFD", "contrast_role": "disease_effect", "tissue": "liver", "omics": "proteomics",
            "feature_id": ["P000", "P250", "P100"], "log2fc": [1.0, 2.0, 0.1], "lfc_se": [0.3, 0.9, 0.3],
            "pvalue": [0.001, 0.04, 0.7], "padj": [0.01, 0.04, 0.9], "padj_method": "fdr_bh", "n_numerator": 3,
            "n_denominator": 3, "method": "limma", "selected": [True, True, False]})
        robust = robustness(self.experiment(), self.run_dir, {"proteomics": diff}, {})["proteomics"]
        rows = robust.set_index("feature_id")
        self.assertEqual(set(rows.index), {"P000", "P250"})
        self.assertFalse(rows.at["P000", "single_subject_driven"])
        self.assertTrue(rows.at["P250", "single_subject_driven"])
        self.assertEqual(rows.at["P250", "driver_sample"], "prot_m6")
        self.assertEqual(rows.at["P250", "n_left_out"], 6)
        self.assertTrue(rows["exaggeration"].ge(1).all())

    def test_feature_without_usable_values_does_not_fail(self) -> None:
        matrix = simulate()
        matrix.loc["P200", ["prot_m0", "prot_m1", "prot_m2"]] = np.nan     # 分母组全部缺失
        self.write_matrix(matrix)
        diff = pd.DataFrame({
            "contrast_id": "HFD_vs_LFD", "contrast_role": "disease_effect", "tissue": "liver", "omics": "proteomics",
            "feature_id": ["P200"], "log2fc": [1.0], "lfc_se": [0.3], "pvalue": [0.001], "padj": [0.01],
            "padj_method": "fdr_bh", "n_numerator": 3, "n_denominator": 0, "method": "limma", "selected": [True]})
        row = robustness(self.experiment(), self.run_dir, {"proteomics": diff}, {})["proteomics"].iloc[0]
        self.assertEqual((row["n_left_out"], row["single_subject_driven"], row["driver_sample"]), (0, False, ""))

    def test_postprocess_writes_power_and_flags(self) -> None:
        matrix = simulate()
        matrix.loc["P250", "prot_m6"] += 6.0
        self.write_matrix(matrix)
        diff = pd.DataFrame({
            "contrast_id": "HFD_vs_LFD", "contrast_role": "disease_effect", "tissue": "liver", "omics": "proteomics",
            "feature_id": ["P000", "P250"], "log2fc": [1.0, 2.0], "lfc_se": [0.3, 0.9], "pvalue": [0.001, 0.04],
            "padj": [0.01, 0.04], "padj_method": "fdr_bh", "n_numerator": 3, "n_denominator": 3, "method": "limma",
            "selected": [True, True]})
        from omicsweave.contract.tables import read_differential, write_differential
        write_differential(diff, self.run_dir / "omics" / "proteomics" / "tables" / "differential.parquet")
        self.set_design(adaptive_shrinkage=True)
        notes = postprocess(self.experiment(), self.run_dir)
        table = read_differential(self.run_dir / "omics" / "proteomics" / "tables" / "differential.parquet")
        self.assertEqual(set(table["estimand_type"]), {"pairwise"})
        self.assertEqual(table.set_index("feature_id").at["P250", "flags"], "single_subject_driven")
        self.assertTrue(any("single_subject_driven" in n for n in notes))
        power = pd.read_parquet(self.run_dir / "qc" / "design_power.parquet").set_index("contrast_id")
        self.assertEqual((power.at["HFD_vs_LFD", "n_numerator"], power.at["HFD_vs_LFD", "n_denominator"]), (3, 3))
        self.assertAlmostEqual(power.at["HFD_vs_LFD", "median_residual_sd"], NOISE, delta=0.15)
        self.assertGreater(power.at["HFD_vs_LFD", "mde_80"], power.at["GLP1_vs_Saline_in_HFD", "mde_80"] * 0.5)
        self.assertTrue((self.run_dir / "qc" / "contrast_reliability.parquet").is_file())

    def test_multi_tracer_mfa_robustness_stratifies_and_ignores_unconstrained_tracers(self) -> None:
        """多示踪剂 MFA 的 E5 稳健性只纳入可辨识且成对的示踪剂并按示踪剂分层，不被不可辨识或未成对示踪剂污染。"""
        exp_dir = self.root / "experiments" / EXP
        subj = pd.DataFrame([
            ("lac_h1", "HFD", "Saline"), ("lac_h2", "HFD", "Saline"),
            ("lac_l1", "LFD", "Saline"), ("lac_l2", "LFD", "Saline"),
            ("gly_h1", "HFD", "Saline"), ("gly_h2", "HFD", "Saline"), ("gly_h3", "HFD", "Saline"),
            ("gly_l1", "LFD", "Saline"), ("gly_l2", "LFD", "Saline"), ("gly_l3", "LFD", "Saline"),
            ("ala_h1", "HFD", "Saline"), ("ala_h2", "HFD", "Saline"), ("ala_h3", "HFD", "Saline"),
            ("ala_l1", "LFD", "Saline"), ("ala_l2", "LFD", "Saline"), ("ala_l3", "LFD", "Saline"),
            ("n15_l1", "LFD", "Saline"), ("n15_l2", "LFD", "Saline"),  # 只有 LFD，未成对
        ], columns=["subject_id", "diet", "treatment"])
        subj.to_csv(exp_dir / "subjects.tsv", sep="\t", index=False)
        proto_of = {"lac": "lactate", "gly": "glycerol", "ala": "alanine", "n15": "alanine_13c15n"}
        samp = pd.DataFrame([
            {"sample_id": s, "subject_id": s, "omics": "mfa", "tissue": "heart", "protocol": proto_of[s[:3]]}
            for s in subj["subject_id"]
        ])
        samp.to_csv(exp_dir / "samples.tsv", sep="\t", index=False)

        cols = list(subj["subject_id"])
        # R_STABLE：仅 lactate 可辨识且极稳定；不可辨识的 alanine 与未成对的 alanine_13c15n 含极端离群值
        r_stable = [154.0, 155.0, 4.0, 7.0,
                    380.0, 170.0, 330.0, 140.0, 85.0, 320.0,
                    0.3, -0.6, 1.0, 12.0, 26.0, 57.0,
                    1.0, 2145.0]
        # R_OUTLIER：仅 glycerol 可辨识且确实由单只动物 gly_h1=270 驱动；其他不可辨识示踪剂呈反向高值（旧逻辑会掩盖）
        r_outlier = [100.0, 100.0, 10.0, 10.0,
                     270.0, 1.0, 1.0, 1.0, 2.0, 1.0,
                     100.0, 100.0, 100.0, 10.0, 10.0, 10.0,
                     10.0, 10.0]
        # R_MULTI：lactate（尺度 ~1）与 glycerol（尺度 ~100）均可辨识且稳定同向
        r_multi = [2.0, 2.2, 0.9, 1.1,
                   210.0, 220.0, 215.0, 100.0, 105.0, 95.0,
                   50.0, -500.0, 50.0, 50.0, 50.0, 50.0,
                   1.0, 5000.0]
        mat = pd.DataFrame([r_stable, r_outlier, r_multi], index=["R_STABLE", "R_OUTLIER", "R_MULTI"], columns=cols)
        write_dataset(self.run_dir / "omics" / "mfa" / "data" / "heart", Dataset(
            mat, empty_features(mat.index),
            DatasetMeta(omics=Omics.MFA, tissue="heart", value_type=ValueType.FLUX, scale="linear",
                        feature_id_type="mfa_reaction", units="nmol"),
        ))
        diff = pd.DataFrame({
            "contrast_id": "HFD_vs_LFD", "contrast_role": "disease_effect", "tissue": "heart", "omics": "mfa",
            "feature_id": ["R_STABLE", "R_OUTLIER", "R_MULTI"],
            "log2fc": [4.6, 3.0, 1.1], "lfc_se": [0.5, 1.2, 0.2], "pvalue": [0.33, 0.1, 0.01],
            "padj": [np.nan, np.nan, 0.03], "padj_method": "fdr_bh",
            "n_numerator": [2, 3, 5], "n_denominator": [2, 3, 5], "method": "tracer_stratified_permutation",
            "selected": [True, True, True],
            "tracers_informative": ["lactate", "glycerol", "glycerol;lactate"],
        })
        rob = robustness(self.experiment(), self.run_dir, {"mfa": diff}, {})["mfa"].set_index("feature_id")
        # R_STABLE 不被不可辨识/未成对示踪剂误伤
        self.assertFalse(bool(rob.at["R_STABLE", "single_subject_driven"]))
        self.assertEqual((int(rob.at["R_STABLE", "n_numerator"]), int(rob.at["R_STABLE", "n_denominator"])), (2, 2))
        self.assertEqual(int(rob.at["R_STABLE", "n_left_out"]), 4)
        self.assertGreater(float(rob.at["R_STABLE", "loo_min_ratio"]), 0.9)
        self.assertEqual(rob.at["R_STABLE", "driver_sample"], "")
        # R_OUTLIER 在有效示踪剂内精准抓出真正驱动离群的动物 gly_h1
        self.assertTrue(bool(rob.at["R_OUTLIER", "single_subject_driven"]))
        self.assertEqual(rob.at["R_OUTLIER", "driver_sample"], "gly_h1")
        self.assertGreaterEqual(int(rob.at["R_OUTLIER", "loo_sign_flips"]), 1)
        # R_MULTI 跨百倍尺度的两个有效示踪剂保持稳定，且总留一数等于 4 + 6 = 10
        self.assertFalse(bool(rob.at["R_MULTI", "single_subject_driven"]))
        self.assertEqual(int(rob.at["R_MULTI", "n_left_out"]), 10)
        self.assertEqual(int(rob.at["R_MULTI", "loo_sign_flips"]), 0)


class ShrinkageTest(unittest.TestCase):
    def test_lfsr_separates_true_effects(self) -> None:
        rng = np.random.default_rng(3)
        truth = np.where(rng.random(1000) < 0.1, 2.0, 0.0)
        se = np.full(1000, 0.5)
        _, lfsr = adaptive_shrinkage(truth + rng.normal(0, se), se)
        self.assertLess(np.median(lfsr[truth > 0]), 0.05)
        self.assertGreater(np.median(lfsr[truth == 0]), 0.5)

    def test_retrodesign_matches_gelman_carlin(self) -> None:
        power, type_s, exaggeration = retrodesign(2.0, 8.1)
        self.assertAlmostEqual(float(power), 0.06, delta=0.005)
        self.assertAlmostEqual(float(type_s), 0.24, delta=0.01)
        self.assertAlmostEqual(float(exaggeration), 9.5, delta=0.2)


if __name__ == "__main__":
    unittest.main()
