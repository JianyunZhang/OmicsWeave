#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""回归测试：可逆反应净通量必须按拟合迭代逐行配对，再做 IQR 过滤。

旧实现先分别对正向、逆向拟合点做 ``_filter_iqr``（会排序），再逐元素相减，
相当于「第 k 小正向 − 第 k 小逆向」——两者来自不同的拟合迭代。正、逆向
负相关时结果可以反号：``20250425_cimr_wt_mice_dio`` data_20251204 中
hfd_3 kidney--PGI_c 旧值 +703，逐次配对为 −59；2190 个动物×可逆反应中
170 个（8%）反号。

下面的四次拟合是能手算的最小反例（IQR 过滤保留排序后第 2、3 个）::

    迭代  forward  reverse   net
    A        70      100     -30
    B       180      200     -20
    C       290      300     -10
    D      1001        1   +1000   ← 正向最大、逆向最小（负相关）

逐次配对：net 排序 [-30, -20, -10, 1000] → 保留 [-20, -10] → 中位数 −15。
排序配对：forward 保留 [180, 290]，reverse 保留 [100, 200] → [80, 90] → +85。
"""

import contextlib
import io
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd



from omicsweave.mfa.df_analysis import (
    _compute_biological_net_flux,
    _compute_net_flux,
    sample_flux_matrix,
)


FWD = "kidney--PGI_c"
REV = "kidney--PGI_c__R"
FORWARD = [70.0, 180.0, 290.0, 1001.0]
REVERSE = [100.0, 200.0, 300.0, 1.0]


def _fits(forward, reverse) -> pd.DataFrame:
    return pd.DataFrame({
        "loss": np.linspace(0.1, 0.2, len(forward)),
        FWD: forward,
        REV: reverse,
        # 不可逆反应：没有 __R 列，走正向路径
        "kidney--HEX_c": 5.0 + np.arange(len(forward)),
    })


def _flux_df(fits: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = {sid: df.drop(columns=["loss"]).iloc[0] for sid, df in fits.items()}
    return pd.DataFrame.from_dict(rows, orient="index")


class IterationPairedNetFluxTest(unittest.TestCase):
    def setUp(self) -> None:
        self.fits = {"hfd_3": _fits(FORWARD, REVERSE)}
        self.flux_df = _flux_df(self.fits)

    def _net(self, fits, sample_ids, **kwargs) -> np.ndarray:
        return _compute_net_flux(fits, self.flux_df, sample_ids, "kidney", "PGI_c", **kwargs)

    def test_anticorrelated_fits_keep_the_true_sign(self) -> None:
        values = self._net(self.fits, ["hfd_3"])
        np.testing.assert_allclose(values, [-20.0, -10.0])

        sorted_pairing = (np.sort(FORWARD)[1:3] - np.sort(REVERSE)[1:3])
        self.assertTrue((sorted_pairing > 0).all(),
                        "反例本身必须让排序配对给出相反符号")

        animal = _compute_biological_net_flux(
            self.fits, self.flux_df, ["hfd_3"], "kidney", "PGI_c",
        )
        np.testing.assert_allclose(animal, [-15.0])

    def test_iqr_filter_is_applied_to_net_values_not_to_each_direction(self) -> None:
        unfiltered = self._net(self.fits, ["hfd_3"], iqr_filter=False)
        np.testing.assert_allclose(unfiltered, [-30.0, -20.0, -10.0, 1000.0])

    def test_each_animal_is_filtered_separately_when_concatenated(self) -> None:
        nets_2 = np.array([100.0, 110.0, 120.0, 130.0])
        reverse_2 = np.array([50.0, 40.0, 30.0, 20.0])
        fits = dict(self.fits, lfd_1=_fits(reverse_2 + nets_2, reverse_2))

        values = self._net(fits, ["hfd_3", "lfd_1"])
        # 合并后再过滤会得到 [-10, 100, 110, 120]
        np.testing.assert_allclose(values, [-20.0, -10.0, 110.0, 120.0])

        animals = _compute_biological_net_flux(
            fits, self.flux_df, ["hfd_3", "lfd_1"], "kidney", "PGI_c",
        )
        np.testing.assert_allclose(animals, [-15.0, 115.0])

    def test_nan_drops_the_whole_iteration(self) -> None:
        # 旧实现对两列各自 dropna，NaN 位置不同就整体错位一行。
        fits = {"hfd_3": _fits(FORWARD + [5000.0, np.nan], REVERSE + [np.nan, 7.0])}
        values = self._net(fits, ["hfd_3"])
        np.testing.assert_allclose(values, [-20.0, -10.0])

    def test_forward_only_path_is_unchanged(self) -> None:
        values = self._net(self.fits, ["hfd_3"], use_net_flux=False)
        np.testing.assert_allclose(values, [180.0, 290.0])

        irreversible = _compute_net_flux(
            self.fits, self.flux_df, ["hfd_3"], "kidney", "HEX_c",
        )
        np.testing.assert_allclose(irreversible, [6.0, 7.0])

    def test_fallback_without_fits_uses_best_fit_net(self) -> None:
        animal = _compute_biological_net_flux(
            {}, self.flux_df, ["hfd_3"], "kidney", "PGI_c",
        )
        np.testing.assert_allclose(animal, [70.0 - 100.0])

        forward_only = _compute_biological_net_flux(
            {}, self.flux_df, ["hfd_3"], "kidney", "PGI_c", use_net_flux=False,
        )
        np.testing.assert_allclose(forward_only, [70.0])

    def test_sample_flux_matrix_reports_the_paired_net_flux(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            with contextlib.redirect_stdout(io.StringIO()):
                sample_flux_matrix(self.flux_df, self.fits, {"hfd_3": "hfd"},
                                   ["kidney"], Path(temp_dir))
            matrix = pd.read_csv(
                Path(temp_dir) / "01_df_results" / "all_tissues"
                / "flux_matrix" / "sample_net_flux.csv"
            ).set_index("reaction")

        row = matrix.loc["PGI_c"]
        self.assertAlmostEqual(row["net_flux"], -15.0)
        self.assertEqual(row["n_fit_points"], 2)
        self.assertAlmostEqual(row["fit_p25"], -17.5)
        self.assertAlmostEqual(row["fit_p75"], -12.5)


if __name__ == "__main__":
    unittest.main()
