#!/usr/bin/env python3
"""统一绘图主题：字体回退、样式与 rcParams。每个测试都在 rc_context 中运行，结束后恢复全局设置。"""

from __future__ import annotations

import tempfile
import unittest
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from omicsweave.visualization.theme import CJK_FONTS, LATIN_FONTS, apply_theme, font_family, has_cjk_font


class ThemeTest(unittest.TestCase):
    def test_latin_font_first_then_cjk(self) -> None:
        family = font_family()
        self.assertIn(family[0], LATIN_FONTS)
        self.assertTrue(all(name in CJK_FONTS for name in family[1:]))

    def test_apply_theme_settings(self) -> None:
        with plt.rc_context():
            apply_theme(style="seaborn-v0_8-whitegrid", dpi=200, save_dpi=200, rc={"axes.titlesize": 14})
            self.assertEqual(plt.rcParams["font.family"], font_family())
            self.assertEqual(plt.rcParams["savefig.dpi"], 200)
            self.assertEqual(plt.rcParams["axes.titlesize"], 14)
            self.assertEqual(plt.rcParams["pdf.fonttype"], 42)
            self.assertFalse(plt.rcParams["axes.unicode_minus"])
            self.assertTrue(plt.rcParams["axes.grid"])  # 来自白色网格样式

    @unittest.skipUnless(has_cjk_font(), "本机没有中文字体")
    def test_chinese_labels_render_without_missing_glyphs(self) -> None:
        with plt.rc_context(), tempfile.TemporaryDirectory() as tmp, warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            apply_theme()
            fig, ax = plt.subplots()
            ax.set_title("肝脏 Liver：HFD vs LFD")
            ax.set_xlabel("log2FC −1.5")
            fig.savefig(Path(tmp) / "probe.png")
            plt.close(fig)
        missing = [str(w.message) for w in caught if "Glyph" in str(w.message)]
        self.assertEqual(missing, [])


if __name__ == "__main__":
    unittest.main()
