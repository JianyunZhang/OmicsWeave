#!/usr/bin/env python3
"""MFA 通用样本质检、比较级尺度漂移（flux_scale_drift）与多示踪剂单组贴边测试。"""

from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from omicsweave.mfa import multi_tracer as mt
from omicsweave.mfa.sample_qc import (
    SAMPLE_QUALITY_COLUMNS,
    apply_scale_drift_flags,
    evaluate_mfa_quality,
)


def _synthetic_animals(
    group_scales: dict[str, tuple[list[str], float]],
    n_reactions: int = 30,
    seed: int = 11,
) -> pd.DataFrame:
    """构造多组 MFA sample_net_flux 长表（默认 30 条活跃反应，尺度按 group_scales 设定）。"""
    rng = np.random.default_rng(seed)
    base_flux = np.linspace(20.0, 200.0, n_reactions)
    rows = []
    for cond, (sids, scale) in group_scales.items():
        for sid in sids:
            for r_idx, b in enumerate(base_flux):
                val = b * scale * float(rng.lognormal(0.0, 0.04))
                rows.append({
                    "sample_id": sid,
                    "condition": cond,
                    "tissue": "liver",
                    "reaction": f"RXN_{r_idx:02d}",
                    "flux_id": f"liver--RXN_{r_idx:02d}",
                    "net_flux": val,
                    "at_bound": False,
                    "at_upper_bound": False,
                })
    return pd.DataFrame(rows)


class MFASampleQCTest(unittest.TestCase):
    def test_schema_matches_core_reliability(self) -> None:
        self.assertEqual(SAMPLE_QUALITY_COLUMNS, ["sample_id", "metric", "value", "robust_z", "flag", "reason"])

    def test_clean_cohort_has_no_sample_flags_or_scale_drift(self) -> None:
        g1 = [f"s1_{i}" for i in range(4)]
        g2 = [f"s2_{i}" for i in range(4)]
        animals = _synthetic_animals({"G1": (g1, 1.0), "G2": (g2, 1.02)})
        contrasts = {"G2_vs_G1": ({**dict.fromkeys(g1, "G1"), **dict.fromkeys(g2, "G2")}, "G2", "G1")}
        loss = pd.DataFrame({"sample": g1 + g2, "min_loss": [10.0, 10.5, 9.8, 10.2, 10.1, 9.9, 10.4, 10.0]})

        res = evaluate_mfa_quality(animals, contrasts, loss_df=loss)
        self.assertEqual(res.flagged_samples, set())
        self.assertEqual(res.excluded_samples, set())
        self.assertEqual(res.drift_contrasts, set())
        self.assertFalse(bool(res.contrast_scale["scale_drift"].any()))

    def test_detects_sample_scale_ub_and_loss_outliers(self) -> None:
        sids = [f"m{i}" for i in range(8)]
        animals = _synthetic_animals({"ctrl": (sids, 1.0)}, n_reactions=40, seed=5)
        # m0 全局通量放大 2.2 倍 -> flux_scale_factor flag
        animals.loc[animals["sample_id"] == "m0", "net_flux"] *= 2.2
        # m1 有 6 条反应撞上界 3000 -> upper_bound_saturation flag
        m1_mask = (animals["sample_id"] == "m1") & (animals["reaction"].isin([f"RXN_{i:02d}" for i in range(6)]))
        animals.loc[m1_mask, "net_flux"] = 2995.0
        animals.loc[m1_mask, "at_upper_bound"] = True
        animals.loc[m1_mask, "at_bound"] = True
        # m2 拟合残差 loss 暴增 5 倍 -> fitting_loss exclude
        loss = pd.DataFrame({
            "sample": sids,
            "min_loss": [10.0, 10.2, 55.0, 9.8, 10.1, 10.3, 9.9, 10.0],
        })

        res = evaluate_mfa_quality(animals, None, loss_df=loss)
        self.assertIn("m0", res.flagged_samples)
        self.assertIn("m1", res.flagged_samples)
        self.assertIn("m2", res.excluded_samples)
        clean = set(sids) - {"m0", "m1", "m2"}
        self.assertEqual((res.flagged_samples | res.excluded_samples) & clean, set())

    def test_contrast_scale_drift_detects_batch_inflation_and_flags_table(self) -> None:
        w1 = [f"w1_{i}" for i in range(4)]
        w3 = [f"w3_{i}" for i in range(4)]
        w1_ctrl = [f"c1_{i}" for i in range(4)]
        # 3W 整鼠通量相对 1W 系统性放大 1.35 倍（log2 = +0.433）；同周期处理组仅 2 条反应变化
        animals = _synthetic_animals({"W1_Sal": (w1_ctrl, 1.0), "W1_Drug": (w1, 1.0), "W3_Drug": (w3, 1.35)},
                                     n_reactions=30, seed=9)
        animals.loc[animals["sample_id"].isin(w1) & animals["reaction"].isin(["RXN_00", "RXN_01"]), "net_flux"] *= 2.5

        contrasts = {
            "Drug_3W_vs_1W": ({**dict.fromkeys(w1, "W1_Drug"), **dict.fromkeys(w3, "W3_Drug")}, "W3_Drug", "W1_Drug"),
            "Drug_vs_Sal_1W": ({**dict.fromkeys(w1_ctrl, "W1_Sal"), **dict.fromkeys(w1, "W1_Drug")}, "W1_Drug", "W1_Sal"),
        }
        res = evaluate_mfa_quality(animals, contrasts)
        self.assertEqual(res.drift_contrasts, {("Drug_3W_vs_1W", "single")})

        cs = res.contrast_scale.set_index("contrast_id")
        self.assertTrue(bool(cs.at["Drug_3W_vs_1W", "scale_drift"]))
        self.assertEqual(cs.at["Drug_3W_vs_1W", "flag"], "flux_scale_drift")
        self.assertFalse(bool(cs.at["Drug_vs_Sal_1W", "scale_drift"]))

        diff = pd.DataFrame({
            "contrast_id": ["Drug_3W_vs_1W", "Drug_vs_Sal_1W"],
            "tissue": ["liver", "liver"],
            "feature_id": ["RXN_00", "RXN_00"],
            "flags": ["", "existing_flag"],
        })
        flagged = apply_scale_drift_flags(diff, res.contrast_scale)
        self.assertEqual(flagged.loc[0, "flags"], "flux_scale_drift")
        self.assertEqual(flagged.loc[1, "flags"], "existing_flag")


class MultiTracerContrastBoundPinningTest(unittest.TestCase):
    def test_pool_contrast_blocks_tracer_with_group_bound_pinning(self) -> None:
        # 全局覆盖表把 glucose 记为 well_constrained（模拟 28 鼠全局稀释），
        # 但当前比较的 numerator 组有 2/3 只动物 at_bound=True -> 必须在合并时剔除 glucose 并记入 flags
        animals = pd.DataFrame([
            # glucose: num 2/3 pinned, den 0/3 pinned
            {"sample_id": "g_n1", "tracer": "glucose", "tissue": "liver", "reaction": "G6Pase_c", "net_flux": 3000.0, "fit_std": 0.1, "at_bound": True},
            {"sample_id": "g_n2", "tracer": "glucose", "tissue": "liver", "reaction": "G6Pase_c", "net_flux": 2966.0, "fit_std": 5.0, "at_bound": True},
            {"sample_id": "g_n3", "tracer": "glucose", "tissue": "liver", "reaction": "G6Pase_c", "net_flux": 1700.0, "fit_std": 10.0, "at_bound": False},
            {"sample_id": "g_d1", "tracer": "glucose", "tissue": "liver", "reaction": "G6Pase_c", "net_flux": 1100.0, "fit_std": 10.0, "at_bound": False},
            {"sample_id": "g_d2", "tracer": "glucose", "tissue": "liver", "reaction": "G6Pase_c", "net_flux": 1200.0, "fit_std": 10.0, "at_bound": False},
            {"sample_id": "g_d3", "tracer": "glucose", "tissue": "liver", "reaction": "G6Pase_c", "net_flux": 1150.0, "fit_std": 10.0, "at_bound": False},
            # glycerol: 全部未贴边
            {"sample_id": "y_n1", "tracer": "glycerol", "tissue": "liver", "reaction": "G6Pase_c", "net_flux": 900.0, "fit_std": 8.0, "at_bound": False},
            {"sample_id": "y_n2", "tracer": "glycerol", "tissue": "liver", "reaction": "G6Pase_c", "net_flux": 920.0, "fit_std": 8.0, "at_bound": False},
            {"sample_id": "y_n3", "tracer": "glycerol", "tissue": "liver", "reaction": "G6Pase_c", "net_flux": 880.0, "fit_std": 8.0, "at_bound": False},
            {"sample_id": "y_d1", "tracer": "glycerol", "tissue": "liver", "reaction": "G6Pase_c", "net_flux": 700.0, "fit_std": 8.0, "at_bound": False},
            {"sample_id": "y_d2", "tracer": "glycerol", "tissue": "liver", "reaction": "G6Pase_c", "net_flux": 710.0, "fit_std": 8.0, "at_bound": False},
            {"sample_id": "y_d3", "tracer": "glycerol", "tissue": "liver", "reaction": "G6Pase_c", "net_flux": 690.0, "fit_std": 8.0, "at_bound": False},
        ])
        coverage = pd.DataFrame([
            {"tissue": "liver", "reaction": "G6Pase_c", "tracer": "glucose", "identifiability": "well_constrained", "informative": True},
            {"tissue": "liver", "reaction": "G6Pase_c", "tracer": "glycerol", "identifiability": "well_constrained", "informative": True},
        ])
        groups = {
            **dict.fromkeys(["g_n1", "g_n2", "g_n3", "y_n1", "y_n2", "y_n3"], "3W"),
            **dict.fromkeys(["g_d1", "g_d2", "g_d3", "y_d1", "y_d2", "y_d3"], "1W"),
        }
        pooled = mt.pool_contrast(animals, groups, "3W", "1W", coverage, mt.PoolingSettings(), seed_key="3W_vs_1W")
        row = pooled.iloc[0]
        self.assertEqual(row["tracers_informative"], "glycerol")
        self.assertIn("glucose:contrast_bound_pinned", row["tracers_excluded"])

        # 若仅保留 glucose 单示踪剂，该比较直接判为 bound_pinned / not_identifiable
        only_glc = animals[animals["tracer"] == "glucose"].copy()
        cov_glc = coverage[coverage["tracer"] == "glucose"].copy()
        pooled_glc = mt.pool_contrast(only_glc, groups, "3W", "1W", cov_glc, mt.PoolingSettings(), seed_key="3W_vs_1W")
        r0 = pooled_glc.iloc[0]
        self.assertEqual(r0["identifiability"], "bound_pinned")
        self.assertEqual(r0["tracer_evidence"], mt.NOT_IDENTIFIABLE)
        self.assertIn("glucose:contrast_bound_pinned", r0["tracers_excluded"])


if __name__ == "__main__":
    unittest.main()
