#!/usr/bin/env python3
"""通用主题预设、语义色板与期刊栏宽几何测试（不依赖 pydantic/seaborn）。"""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from omicsweave.visualization.figures import (
    DOUBLE_COL_IN,
    SINGLE_COL_IN,
    journal_figsize,
    save_stem,
)
from omicsweave.visualization.palettes import (
    EVIDENCE_GRADE_COLORS,
    NATURE_EDITORIAL,
    get_palette,
    tier_tokens,
    tint_color,
)
from omicsweave.visualization.theme import (
    DEFAULT_PRESET,
    apply_theme,
    get_preset,
    list_presets,
    themed,
)


class ThemePresetsTest(unittest.TestCase):
    def test_default_preset_is_nature_editorial(self) -> None:
        self.assertEqual(DEFAULT_PRESET, "nature_editorial")
        self.assertEqual(set(list_presets()), {"nature_editorial", "metro_atlas", "slide_dark"})
        preset = get_preset()
        self.assertEqual(preset.palette.surface, "#FFFFFF")
        self.assertFalse(preset.rc["axes.spines.top"])
        self.assertFalse(preset.rc["axes.spines.right"])
        self.assertFalse(preset.rc["axes.grid"])
        with self.assertRaises(ValueError):
            get_preset("unknown_style")

    def test_themed_context_restores_rc_and_isolates_legacy_style(self) -> None:
        rc_before = dict(matplotlib.rcParams)
        with plt.rc_context():
            apply_theme(style="seaborn-v0_8-whitegrid")
            self.assertTrue(plt.rcParams["axes.grid"])
            with themed(preset="nature_editorial") as style:
                self.assertEqual(style.name, "nature_editorial")
                self.assertFalse(plt.rcParams["axes.grid"])
                self.assertEqual(plt.rcParams["savefig.dpi"], 300)
                self.assertEqual(plt.rcParams["pdf.fonttype"], 42)
                self.assertEqual(plt.rcParams["svg.fonttype"], "none")
            with themed(preset="slide_dark") as dark:
                self.assertEqual(dark.palette.surface, "#111827")
                self.assertEqual(plt.rcParams["figure.facecolor"], "#111827")
        self.assertEqual(dict(matplotlib.rcParams), rc_before)

    def test_monochrome_tier_tokens_and_colorblind_safe_grades(self) -> None:
        tok_light = tier_tokens("nature_editorial")
        tok_dark = tier_tokens("slide_dark")
        self.assertEqual([tok_light[k]["dots"] for k in ("L3", "L2", "L1")], [3, 2, 1])
        self.assertEqual([tok_dark[k]["dots"] for k in ("L3", "L2", "L1")], [3, 2, 1])
        # 证据等级不包含高饱和红/蓝方向色，避免色相冲突
        for k in ("L3", "L2", "L1"):
            self.assertNotIn(tok_light[k]["bg"], (NATURE_EDITORIAL.up, NATURE_EDITORIAL.down))
        # scoring_plots 不再使用红绿对比色 (#1B7837)
        self.assertNotEqual(EVIDENCE_GRADE_COLORS["Replicated"], "#1B7837")
        self.assertEqual(get_palette("nature_editorial").up, NATURE_EDITORIAL.up)
        self.assertEqual(len(tint_color("#2166AC", 0.1)), 3)

    def test_journal_figsize_and_deterministic_save_stem(self) -> None:
        w1, h1 = journal_figsize(1, aspect=0.8)
        w2, h2 = journal_figsize(2, aspect=0.65)
        self.assertAlmostEqual(w1, SINGLE_COL_IN, places=4)
        self.assertAlmostEqual(w2, DOUBLE_COL_IN, places=4)
        self.assertGreater(h2, h1)
        with tempfile.TemporaryDirectory() as tmp:
            with themed(preset="nature_editorial"):
                fig, ax = plt.subplots(figsize=(3.5, 2.5))
                ax.plot([0, 1], [0, 1])
                p1 = save_stem(fig, Path(tmp) / "a", formats=("pdf", "png"))
            with themed(preset="nature_editorial"):
                fig, ax = plt.subplots(figsize=(3.5, 2.5))
                ax.plot([0, 1], [0, 1])
                p2 = save_stem(fig, Path(tmp) / "b", formats=("pdf", "png"))
            self.assertEqual(p1[0].read_bytes(), p2[0].read_bytes())
            self.assertEqual(p1[1].read_bytes(), p2[1].read_bytes())


if __name__ == "__main__":
    unittest.main()
