#!/usr/bin/env python3
"""通用火山图：颜色只取决于 selected，纵轴的退回顺序，越界点钉在边界上，文件名与逐字节可重现，线程安全。"""

from __future__ import annotations

import tempfile
import unittest
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from matplotlib.colors import to_hex

from omicsweave.core.model import Contrast
from omicsweave.visualization.figures import figure_stem
from omicsweave.visualization.theme import apply_theme, themed
from omicsweave.visualization.volcano import COLORS, Thresholds, contrast_label, plot_volcano, write_volcanoes


def assert_volcanoes(test: unittest.TestCase, omics_dir: Path) -> None:
    """结果包中一个组学的 figures/ 恰好有差异结果表中每个组织 × 比较的 PDF 与 PNG 火山图。"""
    table = pd.read_parquet(omics_dir / "tables" / "differential.parquet")
    panels = table[["tissue", "contrast_id"]].drop_duplicates().itertuples(index=False)
    expected = {f"{figure_stem('volcano', t, c)}.{fmt}" for t, c in panels for fmt in ("pdf", "png")}
    test.assertTrue(expected)
    test.assertEqual({p.name for p in (omics_dir / "figures").iterdir()}, expected)


def _table(log2fc, *, pvalue=None, padj=None, selected=None, effect_size=None, tissue="liver",
           contrast_id="HFD_vs_LFD") -> pd.DataFrame:
    n = len(log2fc)
    nan = [np.nan] * n
    return pd.DataFrame({
        "contrast_id": contrast_id, "contrast_role": "disease_effect", "tissue": tissue, "omics": "transcriptomics",
        "feature_id": [f"F{i}" for i in range(n)], "feature_name": [f"gene{i}" for i in range(n)],
        "log2fc": log2fc, "pvalue": pvalue if pvalue is not None else nan,
        "padj": padj if padj is not None else nan, "padj_method": "fdr_bh",
        "effect_size": effect_size if effect_size is not None else nan, "effect_size_type": "hedges_g",
        "n_numerator": 4, "n_denominator": 4, "method": "deseq2",
        "selected": selected if selected is not None else [False] * n, "selection_rule": "padj<0.05 and |log2FC|>1",
    })


def _legend(fig) -> list[str]:
    return [t.get_text() for t in fig.axes[0].get_legend().get_texts()]


def _footnote(fig) -> str:
    return fig.axes[0].texts[-1].get_text()


def _points(fig, color: str) -> np.ndarray:
    """某颜色实心圆点的坐标。"""
    found = [c.get_offsets() for c in fig.axes[0].collections
             if len(c.get_facecolors()) and to_hex(c.get_facecolors()[0], keep_alpha=False) == color.lower()]
    return np.vstack(found) if found else np.empty((0, 2))


class PlotVolcanoTest(unittest.TestCase):
    def test_colors_follow_selected_not_thresholds(self) -> None:
        # 第 3 个特征超过两条阈值线但未入选（例如规则还要求 VIP），应画成灰色
        table = _table([2.0, -1.5, 3.0, 0.1], padj=[0.001, 0.01, 0.0001, 0.5], selected=[True, True, False, False])
        with themed():
            fig = plot_volcano(table, thresholds=Thresholds(p=0.05, log2fc=1.0))
        self.assertEqual(_legend(fig), ["Up (1)", "Down (1)", "Not selected (2)"])
        np.testing.assert_allclose(_points(fig, COLORS["up"]), [[2.0, 3.0]])
        np.testing.assert_allclose(_points(fig, COLORS["down"]), [[-1.5, 2.0]])
        self.assertEqual(len(_points(fig, COLORS["ns"])), 2)
        self.assertEqual(fig.axes[0].get_ylabel(), "−log$_{10}$(adjusted p)")
        self.assertIn("Selected: padj<0.05 and |log2FC|>1", _footnote(fig))
        self.assertIn("n = 4 vs 4", _footnote(fig))
        # 标注入选特征，未入选的不标（adjustText 的连线是无文字的注释）
        labels = {t.get_text() for t in fig.axes[0].texts[:-1] if t.get_text()}
        self.assertEqual(labels, {"gene0", "gene1"})

    def test_y_axis_falls_back_to_pvalue_then_effect_size(self) -> None:
        with themed():
            by_p = plot_volcano(_table([1.0, -1.0], pvalue=[0.01, 0.2]))
            by_effect = plot_volcano(_table([1.0, -1.0], effect_size=[-1.2, 0.3]), metric="pvalue",
                                     thresholds=Thresholds(p=0.05, effect_size=0.8))
            nothing = plot_volcano(_table([1.0, -1.0]))
        self.assertEqual(by_p.axes[0].get_ylabel(), "−log$_{10}$(p)")
        self.assertEqual(by_effect.axes[0].get_ylabel(), "|Hedges' g|")
        # 纵轴为效应量时画效应量阈值线，不画 p 值阈值线
        self.assertEqual([line.get_ydata()[0] for line in by_effect.axes[0].lines], [0.8])
        self.assertIsNone(nothing)
        with self.assertRaisesRegex(ValueError, "metric"):
            plot_volcano(_table([1.0]), metric="fdr")

    def test_off_scale_points_are_pinned_and_missing_are_counted(self) -> None:
        table = _table([np.inf, -40.0, 3.0, -2.0, np.nan, 1.0], pvalue=[0.01, 0.02, 0.03, 0.5, 0.01, np.nan],
                       selected=[True, True, True, False, False, False])
        with themed():
            fig = plot_volcano(table, metric="pvalue", x_cap=10)
        ax = fig.axes[0]
        # 横轴半宽取 min(最大有限 |log2FC|, x_cap) = 10，对称于 0
        self.assertAlmostEqual(ax.get_xlim()[0], -ax.get_xlim()[1])
        self.assertAlmostEqual(ax.get_xlim()[1], 10.8)
        pinned = np.vstack([c.get_offsets() for c in ax.collections if c.get_paths() and
                            len(c.get_facecolors()) == 0])
        np.testing.assert_allclose(sorted(pinned[:, 0]), [-10.4, 10.4])
        self.assertIn("Beyond x-axis (2)", _legend(fig))
        self.assertIn("Not shown: 2 features without log2FC or p", _footnote(fig))
        self.assertEqual(ax.get_ylim()[0], 0)

    def test_rank_test_floor_is_stated(self) -> None:
        table = _table([1.0, -1.0], pvalue=[0.1, 0.4]).assign(min_achievable_p=0.1)
        with themed():
            fig = plot_volcano(table, metric="pvalue")
        self.assertIn("Smallest p the rank test can reach at this sample size: 0.1", _footnote(fig))

    def test_contrast_label_and_file_stem(self) -> None:
        contrast = Contrast(id="GLP1_vs_Saline_in_HFD", factor="treatment", numerator="GLP1",
                            denominator="Saline", subset={"diet": "HFD"})
        self.assertEqual(contrast_label(contrast), "GLP1 vs Saline (diet=HFD)")
        self.assertEqual(figure_stem("volcano", "white adipose/epi", "A_vs_B"), "volcano_white_adipose_epi_A_vs_B")


class WriteVolcanoesTest(unittest.TestCase):
    def setUp(self) -> None:
        rng = np.random.default_rng(0)
        n = 300
        log2fc = rng.normal(0, 1.5, n)
        padj = np.clip(rng.uniform(size=n) ** (1 + 2 * np.abs(log2fc)), 1e-12, 1)
        selected = (padj < 0.05) & (np.abs(log2fc) > 1)
        self.table = pd.concat([
            _table(log2fc, padj=padj, pvalue=padj / 2, selected=selected, tissue=tissue, contrast_id=cid)
            for cid in ("B_vs_A", "C_vs_A") for tissue in ("liver", "muscle")
        ] + [_table([0.5], tissue="brain", contrast_id="B_vs_A")], ignore_index=True)
        self.contrasts = [Contrast(id="C_vs_A", factor="g", numerator="C", denominator="A"),
                          Contrast(id="B_vs_A", factor="g", numerator="B", denominator="A")]

    def write(self, directory: Path):
        return write_volcanoes(self.table, directory, contrasts=self.contrasts,
                               thresholds=Thresholds(p=0.05, log2fc=1.0))

    def test_files_order_notes_and_reproducible_bytes(self) -> None:
        rc_before = dict(matplotlib.rcParams)
        with tempfile.TemporaryDirectory() as tmp:
            paths, notes = self.write(Path(tmp, "a"))
            again, _ = self.write(Path(tmp, "b"))
            # 按 contrasts 的顺序，比较内按组织
            self.assertEqual([p.name for p in paths[::2]],
                             ["volcano_liver_C_vs_A.pdf", "volcano_muscle_C_vs_A.pdf",
                              "volcano_liver_B_vs_A.pdf", "volcano_muscle_B_vs_A.pdf"])
            self.assertEqual({p.suffix for p in paths}, {".pdf", ".png"})
            self.assertEqual(notes, ["以下比较 / 组织没有可画的 p 值或效应量，未画火山图：B_vs_A / brain"])
            # 不写入创建时间：重画的图件逐字节一致
            for first, second in zip(paths, again):
                self.assertEqual(first.read_bytes(), second.read_bytes(), first.name)
        self.assertEqual(dict(matplotlib.rcParams), rc_before)

    def test_parallel_threads_match_sequential(self) -> None:
        self.table = self.table[self.table["contrast_id"] == "B_vs_A"]
        with tempfile.TemporaryDirectory() as tmp:
            sequential, _ = self.write(Path(tmp, "seq"))
            with ThreadPoolExecutor(max_workers=2) as pool:
                results = list(pool.map(lambda name: self.write(Path(tmp, name))[0], ["t1", "t2"]))
            for paths in results:
                for expected, actual in zip(sequential, paths):
                    self.assertEqual(expected.read_bytes(), actual.read_bytes(), actual.name)

    def test_global_style_set_by_legacy_modules_does_not_leak(self) -> None:
        # 旧流程的绘图模块在导入时设置全局样式（如代谢组的 seaborn 白色网格）
        self.table = self.table[(self.table["contrast_id"] == "B_vs_A") & (self.table["tissue"] == "liver")]
        with tempfile.TemporaryDirectory() as tmp:
            plain, _ = self.write(Path(tmp, "plain"))
            with matplotlib.rc_context():
                apply_theme(style="seaborn-v0_8-whitegrid")
                styled, _ = self.write(Path(tmp, "styled"))
                with themed():
                    fig = plot_volcano(self.table)
            self.assertFalse(any(line.get_visible() for line in fig.axes[0].xaxis.get_gridlines()))
            for expected, actual in zip(plain, styled):
                self.assertEqual(expected.read_bytes(), actual.read_bytes(), actual.name)


if __name__ == "__main__":
    unittest.main()
