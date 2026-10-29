#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""回归测试：盒约束贴边通量不得被标为 well_constrained。

可辨识性此前只看拟合点的相对 IQR。拟合全部贴在盒约束边界（本项目方向通量
下界 1.0；可逆反应两向同贴下界时净通量为 0）时 IQR = 0，于是被判为
``well_constrained``——``20250425_cimr_wt_mice_dio`` 的 data_20251204 中
69 条 well_constrained 里有 39 条如此（肝 PEPCK_c、PDH_m 每只鼠都是 1.0）。
这里覆盖：

- ``_flag_at_bound`` / ``_net_bound_values`` / ``_infer_fit_box_bounds`` 的判定；
- ``sample_net_flux.csv`` 带 ``at_bound``，``flux_identifiability.csv`` 带
  ``frac_at_bound`` 与 ``iqr_identifiability``，``bound_pinned`` 优先于 IQR 分档；
- 上界贴边（3000）只有在从拟合点推断盒约束时才能识别；
- 差异结果带上 ``bound_pinned`` 与 ``frac_at_bound``，
  ``require_identifiable`` 把 ``bound_pinned`` 与 ``poorly_constrained`` 同等剔除；
- 配置校验拒绝越界参数。
"""

import contextlib
import io
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd



from omicsweave.mfa.df_analysis import (
    UNCONSTRAINED_LABELS,
    _flag_at_bound,
    _infer_fit_box_bounds,
    _net_bound_values,
    differential_flux,
    sample_flux_matrix,
)
from omicsweave.mfa.mfa_experiment_config import _validate_analysis_config


LOWER, UPPER = 1.0, 3000.0
N_FITS = 20
CASES = [f"case_{i}" for i in range(3)]
CONTROLS = [f"control_{i}" for i in range(3)]
SAMPLES = CASES + CONTROLS
CONDITIONS = {sid: ("case" if sid.startswith("case") else "control") for sid in SAMPLES}


def _tight(value: float) -> np.ndarray:
    """被数据良好约束的通量：围绕 value 的窄分布（相对 IQR ≈ 0.01）。"""
    return value + np.linspace(-1.0, 1.0, N_FITS)


def _pinned(value: float) -> np.ndarray:
    return np.full(N_FITS, value)


def _build_fits() -> dict[str, pd.DataFrame]:
    """六只动物、单组织 liver 的方向通量拟合点。"""
    me_level = {"case_0": None, "case_1": None, "case_2": 40.0,
                "control_0": 40.0, "control_1": 40.0, "control_2": 40.0}
    diff_level = {"case_0": None, "case_1": None, "case_2": None,
                  "control_0": 40.0, "control_1": 50.0, "control_2": 60.0}
    fits = {}
    for sid in SAMPLES:
        fits[sid] = pd.DataFrame({
            "loss": np.linspace(0.1, 0.2, N_FITS),
            # 每只鼠都贴下界 → bound_pinned（IQR 单独看是 well_constrained）
            "liver--PEPCK_c": _pinned(LOWER),
            # 真正被约束 → well_constrained
            "liver--CS_m": _tight(50.0),
            # 可逆反应两向同贴下界 → 净通量 0 → bound_pinned
            "liver--LDH_c": _pinned(LOWER),
            "liver--LDH_c__R": _pinned(LOWER),
            # 贴上界 → 仅在推断盒约束时可识别
            "liver--PK_c": _pinned(UPPER),
            # 仅 2/6 只贴边 → 不足半数，保留 IQR 分档
            "liver--ME_c": (_pinned(LOWER) if me_level[sid] is None
                            else _tight(me_level[sid])),
            # case 组贴边、control 组被约束 → frac_at_bound = 0.5，且组间差异巨大
            "liver--DIFF_c": (_pinned(LOWER) if diff_level[sid] is None
                              else _tight(diff_level[sid])),
            # MIX 列不是通量，不得参与盒约束推断
            "MIX:liver--glc_0": np.full(N_FITS, 0.0),
        })
    return fits


def _flux_df(fits: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """最优拟合矩阵：取每只鼠的第一个拟合点。"""
    rows = {sid: df.drop(columns=["loss"]).iloc[0] for sid, df in fits.items()}
    return pd.DataFrame.from_dict(rows, orient="index")


def _run_matrix(output_dir: Path, analysis_config=None) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    fits = _build_fits()
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        sample_flux_matrix(_flux_df(fits), fits, CONDITIONS, ["liver"],
                           output_dir, analysis_config=analysis_config)
    matrix_dir = output_dir / "01_df_results" / "all_tissues" / "flux_matrix"
    matrix = pd.read_csv(matrix_dir / "sample_net_flux.csv")
    ident = pd.read_csv(matrix_dir / "flux_identifiability.csv").set_index("reaction")
    return matrix, ident, buffer.getvalue()


class BoundHelpersTest(unittest.TestCase):
    def test_flag_at_bound_uses_absolute_value_and_tolerance(self) -> None:
        values = [1.0, 1.04, 1.06, 0.0, -1.0, 0.03, 50.0, np.nan]
        flags = _flag_at_bound(values, [0.0, 1.0], 0.05)
        self.assertEqual(flags.tolist(),
                         [True, True, False, True, True, True, False, False])

    def test_net_bound_values_cover_both_directions_of_the_box(self) -> None:
        # 不可逆：lb / ub；可逆：两向同侧 → 0，分处两侧 → ub − lb
        self.assertEqual(_net_bound_values([0.0, 1.0], (LOWER, UPPER)),
                         [0.0, 1.0, 2999.0, 3000.0])
        self.assertEqual(_net_bound_values([0.0, 1.0], None), [0.0, 1.0])

    def test_infer_box_ignores_mix_and_loss_columns(self) -> None:
        self.assertEqual(_infer_fit_box_bounds(_build_fits()), (LOWER, UPPER))
        self.assertIsNone(_infer_fit_box_bounds({}))


class SampleFluxMatrixBoundPinnedTest(unittest.TestCase):
    def test_bound_pinned_takes_precedence_over_iqr_bins(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            matrix, ident, stdout = _run_matrix(Path(temp_dir))

        self.assertIn("at_bound", matrix.columns)
        pepck = matrix[matrix["reaction"] == "PEPCK_c"]
        self.assertTrue(pepck["at_bound"].all())

        for reaction in ("PEPCK_c", "LDH_c", "PK_c", "DIFF_c"):
            self.assertEqual(ident.loc[reaction, "identifiability"], "bound_pinned",
                             reaction)
        # 修复前这条会被判为 well_constrained：零离散度 ≠ 被数据决定
        self.assertEqual(ident.loc["PEPCK_c", "iqr_identifiability"], "well_constrained")
        self.assertEqual(ident.loc["LDH_c", "iqr_identifiability"], "undetermined")
        self.assertAlmostEqual(ident.loc["DIFF_c", "frac_at_bound"], 0.5)

        self.assertEqual(ident.loc["CS_m", "identifiability"], "well_constrained")
        self.assertEqual(ident.loc["CS_m", "frac_at_bound"], 0.0)
        # 不足半数贴边：保留 IQR 分档
        self.assertAlmostEqual(ident.loc["ME_c", "frac_at_bound"], 2 / 6)
        self.assertEqual(ident.loc["ME_c", "identifiability"], "well_constrained")

        # 核心不变量：多数动物贴边的通量一条也不能是 well_constrained
        majority_pinned = ident[ident["frac_at_bound"] >= 0.5]
        self.assertFalse(majority_pinned["identifiability"].eq("well_constrained").any())

        self.assertIn("[!] 4/6", stdout)
        self.assertIn("bound_pinned", stdout)
        # PEPCK_c、PK_c、DIFF_c 的拟合 IQR 中位数都≈0；LDH_c 净通量为 0 → undetermined
        self.assertIn("3 条若只看拟合 IQR 会被误判为 well_constrained", stdout)

    def test_upper_bound_requires_inferred_box(self) -> None:
        config = {"identifiability": {"infer_bounds_from_fits": False}}
        with tempfile.TemporaryDirectory() as temp_dir:
            _, ident, _ = _run_matrix(Path(temp_dir), analysis_config=config)
        self.assertEqual(ident.loc["PK_c", "identifiability"], "well_constrained")
        self.assertEqual(ident.loc["PEPCK_c", "identifiability"], "bound_pinned")

    def test_min_fraction_is_configurable(self) -> None:
        config = {"identifiability": {"bound_pinned_min_fraction": 0.3}}
        with tempfile.TemporaryDirectory() as temp_dir:
            _, ident, _ = _run_matrix(Path(temp_dir), analysis_config=config)
        self.assertEqual(ident.loc["ME_c", "identifiability"], "bound_pinned")


class DifferentialFluxIdentifiabilityTest(unittest.TestCase):
    def _run_differential(self, output_dir: Path, require_identifiable: bool) -> pd.DataFrame:
        fits = _build_fits()
        config = {
            "annotation": {"require_standard_reaction_ids": False},
            "differential_flux": {
                "statistic": "effect_size_only",
                "require_identifiable": require_identifiable,
            },
        }
        with contextlib.redirect_stdout(io.StringIO()):
            sample_flux_matrix(_flux_df(fits), fits, CONDITIONS, ["liver"],
                               output_dir, analysis_config=config)
            result = differential_flux(
                _flux_df(fits), fits, CONDITIONS, ["case", "control"], ["liver"],
                output_dir, comparisons=[("case", "control")], analysis_config=config,
            )
        return result["CASE_vs_CONTROL"].set_index("reaction")

    def test_bound_pinned_propagates_and_is_gated_like_poorly_constrained(self) -> None:
        self.assertIn("bound_pinned", UNCONSTRAINED_LABELS)
        self.assertIn("poorly_constrained", UNCONSTRAINED_LABELS)

        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir)
            annotated = self._run_differential(output_dir, require_identifiable=False)
            csv = pd.read_csv(output_dir / "01_df_results" / "liver" / "differential_flux"
                              / "differential_flux_CASE_vs_CONTROL.csv").set_index("reaction")

        self.assertEqual(annotated.loc["DIFF_c", "identifiability"], "bound_pinned")
        self.assertAlmostEqual(annotated.loc["DIFF_c", "frac_at_bound"], 0.5)
        self.assertEqual(csv.loc["DIFF_c", "identifiability"], "bound_pinned")
        self.assertIn("frac_at_bound", csv.columns)
        # 默认只标注不过滤：巨大的组间差异仍判显著
        self.assertTrue(bool(annotated.loc["DIFF_c", "is_significant"]))

        with tempfile.TemporaryDirectory() as temp_dir:
            gated = self._run_differential(Path(temp_dir), require_identifiable=True)
        self.assertFalse(bool(gated.loc["DIFF_c", "is_significant"]))


class CrossProtocolConsensusGateTest(unittest.TestCase):
    def test_bound_pinned_in_any_protocol_fails_consensus(self) -> None:
        from omicsweave.mfa.mfa_integration import cross_protocol_consensus

        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            labels = {"glucose": "well_constrained", "lactate": "bound_pinned"}
            for pid, label in labels.items():
                harmonized = root / "integration" / "tracer_harmonized" / pid
                harmonized.mkdir(parents=True)
                pd.DataFrame({
                    "flux_id": ["liver--PEPCK_c", "liver--CS_m"],
                    "tissue": ["liver", "liver"],
                    "reaction": ["PEPCK_c", "CS_m"],
                    "log2FC_pooled": [1.2, 0.9],
                    "is_significant_pooled": [True, True],
                    "identifiability": [label, "well_constrained"],
                }).to_csv(harmonized / "differential_flux_meta_HFD_vs_LFD.csv", index=False)

            with contextlib.redirect_stdout(io.StringIO()):
                consensus = cross_protocol_consensus(root, list(labels), {}).set_index("reaction")

        self.assertFalse(bool(consensus.loc["PEPCK_c", "consensus_pass"]))
        self.assertTrue(bool(consensus.loc["PEPCK_c", "excluded_poorly_constrained"]))
        self.assertEqual(consensus.loc["PEPCK_c", "n_protocols_well_constrained"], 1)
        self.assertTrue(bool(consensus.loc["CS_m", "consensus_pass"]))


class IdentifiabilityConfigValidationTest(unittest.TestCase):
    def test_rejects_out_of_range_parameters(self) -> None:
        for bad in ({"bound_pinned_min_fraction": 0.0},
                    {"bound_pinned_min_fraction": 1.5},
                    {"bound_tolerance": -0.1},
                    {"bound_values": 1.0},
                    {"bound_values": ["lower"]}):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                _validate_analysis_config({"identifiability": bad})

    def test_accepts_defaults(self) -> None:
        _validate_analysis_config({"identifiability": {
            "bound_values": [0.0, 1.0], "bound_tolerance": 0.05,
            "infer_bounds_from_fits": True, "bound_pinned_min_fraction": 0.5,
        }})


if __name__ == "__main__":
    unittest.main()
