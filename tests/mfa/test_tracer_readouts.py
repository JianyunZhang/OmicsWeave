#!/usr/bin/env python3
"""G10 中不依赖新数据的部分：通量组、标准化尺度的示踪剂合并、逐只动物的拟合分布。

验收（通用能力规划 G10）：
- 通量组内的反应只计一次检验；化学计量上锁死（逐次拟合中成线性关系）的反应被识别，独立的反应不被合并；
- 绝对尺度相差数十倍、方向一致的两种示踪剂不判为冲突；方向相反时仍为冲突；
- 结果包写出逐只动物的拟合分布与通量组表；档案按净通量的变化定号，并读入误差线与贴边标记。
"""

from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from omicsweave.contract.tables import read_differential
from omicsweave.core.config import load_experiment, load_project
from omicsweave.core.run import run_experiment
from omicsweave.integration.dossier.samples import dossier_samples
from omicsweave.integration.dossier.sources import primary_source
from omicsweave.integration.modes.reaction_dossier import _mfa_stats
from omicsweave.mfa import flux_groups as fg
from omicsweave.mfa import multi_tracer as mt
from tests.mfa.test_mfa_runner import EXP, build_two_factor_project

UNITS = "nmol/min/g"


def fits(n: int = 40, seed: int = 0) -> dict[str, pd.DataFrame]:
    """两个样本的逐次拟合：B = 2A + 1（锁死），C 与 A 无关，D 在所有拟合中不变（贴边）。"""
    rng = np.random.default_rng(seed)
    out = {}
    for sample in ("s1", "s2"):
        a = rng.normal(10, 2, n)
        out[sample] = pd.DataFrame({"liver--A": a, "liver--B": 2 * a + 1, "liver--B__R": np.full(n, 1.0),
                                    "liver--C": rng.normal(5, 1, n), "liver--D": np.full(n, 1.0),
                                    "heart--A": rng.normal(3, 1, n)})
    return out


class FluxGroupTest(unittest.TestCase):
    def test_locked_reactions_are_grouped(self) -> None:
        table = fg.locked_groups({"glucose": fits()})
        self.assertEqual(set(table["reaction"]), {"A", "B"})          # C 独立、D 不变、心脏的 A 是另一个组织
        row = table.set_index("reaction").loc["B"]
        self.assertEqual((row["flux_group"], row["n_members"], row["members"]), ("liver--A", 2, "A;B"))
        self.assertAlmostEqual(row["slope"], 2.0)                       # 净通量 = 2A + 1 − 1
        self.assertAlmostEqual(row["intercept"], 0.0, places=6)

    def test_pair_must_be_locked_in_every_model(self) -> None:
        other = fits(seed=1)
        for frame in other.values():                                    # 第二个模型中 B 与 A 无关
            frame["liver--B"] = np.random.default_rng(3).normal(20, 2, len(frame))
        self.assertTrue(fg.locked_groups({"glucose": fits(), "lactate": other}).empty)
        for frame in other.values():                                    # 第二个模型中 B 不变：不构成反例
            frame["liver--B"] = 7.0
        self.assertEqual(set(fg.locked_groups({"glucose": fits(), "lactate": other})["reaction"]), {"A", "B"})

    def test_grouped_bh_counts_each_group_once(self) -> None:
        p = pd.Series([0.01, 0.01, 0.02, 0.04, np.nan])
        keys = pd.Series(["g1", "g1", "x", "y", "z"])
        adjusted = fg.grouped_bh(p, keys)
        expected = multipletests([0.01, 0.02, 0.04], method="fdr_bh")[1]
        np.testing.assert_allclose(adjusted[:4], [expected[0], expected[0], expected[1], expected[2]])
        self.assertTrue(np.isnan(adjusted[4]))
        self.assertEqual(fg.n_tests(p, keys), 3)


class StandardizedPoolingTest(unittest.TestCase):
    @staticmethod
    def pooled(spec: dict) -> pd.Series:
        rows, groups = [], {}
        for tracer, (num, den) in spec.items():
            for group, values in (("treat", num), ("control", den)):
                for i, value in enumerate(values):
                    sample = f"{tracer}_{group}{i}"
                    groups[sample] = group
                    rows.append({"sample_id": sample, "tracer": tracer, "tissue": "liver", "reaction": "R1",
                                 "net_flux": value, "fit_std": 0.0})
        coverage = mt.tracer_coverage(pd.DataFrame({"tracer": list(spec), "tissue": "liver", "reaction": "R1",
                                                    "identifiability": "well_constrained"}), mt.PoolingSettings())
        return mt.pool_contrast(pd.DataFrame(rows), groups, "treat", "control", coverage, mt.PoolingSettings()).iloc[0]

    def test_same_direction_on_different_scales_is_consistent(self) -> None:
        # 葡萄糖模型的绝对通量约为甘油模型的 50 倍，两者都上升
        row = self.pooled({"glucose": ([520.0, 480.0, 505.0], [400.0, 420.0, 390.0]),
                           "glycerol": ([10.5, 9.6, 10.1], [8.1, 8.4, 7.7])})
        self.assertLess(row["heterogeneity_pvalue"], 0.05)                # 绝对尺度：Q 主要反映尺度差异
        self.assertGreater(row["smd_heterogeneity_pvalue"], 0.10)
        self.assertEqual(row["tracer_evidence"], mt.CONSISTENT)
        self.assertAlmostEqual(row["effect_size"], row["smd"])
        self.assertGreater(row["smd"], 1.0)

    def test_opposite_directions_still_conflict(self) -> None:
        row = self.pooled({"glucose": ([16.0, 16.2, 15.9], [20.0, 20.1, 19.8]),
                           "lactate": ([28.0, 28.3, 27.9], [20.0, 20.2, 19.9])})
        self.assertLess(row["smd_heterogeneity_pvalue"], 0.10)
        self.assertEqual(row["tracer_evidence"], mt.CONFLICT)

    def test_magnitude_heterogeneity_is_flagged_not_blocking(self) -> None:
        spec = {"glucose": ([12.0, 13.0, 12.5], [10.0, 10.4, 9.9]), "lactate": ([14.0, 11.0, 13.5], [10.5, 9.0, 11.0])}
        row = self.pooled(spec)
        self.assertEqual(row["tracer_evidence"], mt.CONSISTENT)
        table = mt.to_differential(pd.DataFrame([row]), contrast_id="c", contrast_role="other", flux_units=UNITS,
                                   settings=mt.PoolingSettings())
        flags = table.at[0, "flags"]
        if row["smd_heterogeneity_pvalue"] < 0.10:
            self.assertIn("tracer_magnitude_heterogeneity", flags)
        self.assertNotIn("tracer_conflict", flags)


def _add_locked_reaction(root: Path) -> None:
    """在旧版拟合文件中加一条与肝脏 PDH 锁死的反应（CS = 2·PDH + 1）。"""
    path = root / "experiments" / EXP / "mfa" / "flux_raw_data.xlsx"
    sheets = pd.read_excel(path, sheet_name=None)
    for name, frame in sheets.items():
        if name != "sheet_information":
            frame["liver--CS_m"] = 2 * frame["liver--PDH_m"] + 1
    with pd.ExcelWriter(path) as writer:
        for name, frame in sheets.items():
            frame.to_excel(writer, sheet_name=name, index=False)


class BundleTest(unittest.TestCase):
    def test_tables_and_grouped_correction(self) -> None:
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp).resolve()
            build_two_factor_project(root)
            _add_locked_reaction(root)
            with contextlib.redirect_stdout(io.StringIO()):
                run_dir, notes = run_experiment(EXP, ["mfa"], root=root)
            mfa = run_dir / "omics" / "mfa"
            groups = pd.read_parquet(mfa / "tables" / "flux_groups.parquet")
            fits_table = pd.read_parquet(mfa / "tables" / "sample_flux_fits.parquet")
            table = read_differential(mfa / "tables" / "differential.parquet")
            experiment = load_experiment(load_project(root), EXP)
            source = primary_source(experiment, run_dir, "A", {})
            members = pd.DataFrame({"layer": "mfa", "tissue": "liver", "feature_id": ["PDH_m", "CS_m"]})
            samples, _ = dossier_samples(source, members)
        self.assertEqual(set(groups.loc[groups["tissue"] == "liver", "reaction"]), {"CS_m", "PDH_m"})
        self.assertEqual(list(fits_table.columns), ["sample_id", "tissue", "reaction", "tracer", "net_flux", "n_fits",
                                                    "fit_q25", "fit_q75", "fit_p25", "fit_p75", "fit_std", "at_bound"])
        self.assertTrue((fits_table["n_fits"] == 6).all())
        self.assertTrue((fits_table["fit_q25"] <= fits_table["fit_q75"]).all())
        liver = table[(table["tissue"] == "liver") & (table["contrast_id"] == "HFD_vs_LFD")].set_index("feature_id")
        self.assertEqual(liver.at["CS_m", "flux_group"], liver.at["PDH_m", "flux_group"])
        self.assertEqual(int(liver.at["CS_m", "flux_group_size"]), 2)
        tested = liver[liver["padj"].notna()]
        if {"CS_m", "PDH_m"} <= set(tested.index):
            self.assertEqual(tested.at["CS_m", "padj"], tested.at["PDH_m", "padj"])
        self.assertTrue(any(n.startswith("通量组：") for n in notes))
        # 档案读入逐只动物的拟合四分位数与贴边标记
        self.assertTrue(samples["fit_low"].notna().all())
        self.assertTrue((samples["fit_low"] <= samples["fit_high"]).all())
        self.assertTrue(samples["at_bound"].notna().all())


class DossierSignTest(unittest.TestCase):
    def test_sign_change_follows_net_flux(self) -> None:
        row = pd.Series({"log2fc": -1.836, "mean_numerator": 7.976, "mean_denominator": -28.483, "flags": "",
                         "heterogeneity_pvalue": 0.01, "smd_heterogeneity_pvalue": 0.6,
                         "flux_group": "heart--ENO_c", "flux_group_size": 2})
        stats = _mfa_stats(row)
        self.assertAlmostEqual(stats["log2fc"], 1.836)                    # −28.5 → +8.0：净通量上升
        self.assertIn("net_flux_sign_change", stats["flags"])
        self.assertIn("flux_group:ENO_c", stats["flags"])
        self.assertEqual(stats["heterogeneity_pvalue"], 0.6)              # 标准化尺度
        same_sign = _mfa_stats(pd.Series({"log2fc": -0.5, "mean_numerator": 7.0, "mean_denominator": 10.0}))
        self.assertEqual(same_sign["log2fc"], -0.5)
        self.assertNotIn("net_flux_sign_change", same_sign["flags"])


if __name__ == "__main__":
    unittest.main()
