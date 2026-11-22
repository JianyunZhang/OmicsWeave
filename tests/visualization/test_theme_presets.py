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
    ECONOMIST,
    EVIDENCE_GRADE_COLORS,
    METRO_ATLAS,
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
    def test_default_preset_is_metro_atlas_and_economist_available(self) -> None:
        self.assertEqual(DEFAULT_PRESET, "metro_atlas")
        self.assertEqual(
            set(list_presets()),
            {"economist", "nature_editorial", "metro_atlas", "slide_dark", "poster_editorial"},
        )
        preset = get_preset()
        self.assertEqual(preset.name, "metro_atlas")
        self.assertEqual(preset.palette.surface, "#FFFFFF")
        self.assertEqual(preset.palette.up, METRO_ATLAS.up)
        self.assertEqual(preset.palette.down, METRO_ATLAS.down)
        self.assertFalse(preset.rc["axes.spines.top"])
        self.assertFalse(preset.rc["axes.spines.right"])
        self.assertFalse(preset.rc["axes.grid"])

        econ = get_preset("economist")
        self.assertEqual(econ.palette.surface, "#EBEBDF")
        self.assertEqual(econ.palette.up, "#E3120B")
        self.assertEqual(econ.palette.down, "#006BA2")
        self.assertEqual(econ.palette.numerator, "#E3120B")
        self.assertEqual(econ.palette.denominator, "#006BA2")
        self.assertFalse(econ.rc["axes.spines.left"])
        self.assertTrue(econ.rc["axes.spines.bottom"])
        self.assertTrue(econ.rc["axes.grid"])
        self.assertEqual(econ.rc["axes.grid.axis"], "y")

        nat = get_preset("nature_editorial")
        self.assertEqual(nat.palette.surface, "#FFFFFF")
        self.assertFalse(nat.rc["axes.grid"])
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

    def test_poster_preset_uses_horizontal_guides_and_preserves_direction_colors(self) -> None:
        with themed(preset="poster_editorial") as poster:
            self.assertEqual(poster.palette.up, NATURE_EDITORIAL.up)
            self.assertEqual(poster.palette.down, NATURE_EDITORIAL.down)
            self.assertTrue(plt.rcParams["axes.grid"])
            self.assertEqual(plt.rcParams["axes.grid.axis"], "y")
            self.assertEqual(plt.rcParams["axes.titlelocation"], "left")
        with themed(preset="nature_editorial"):
            self.assertFalse(plt.rcParams["axes.grid"])
        with themed(preset="economist"):
            self.assertTrue(plt.rcParams["axes.grid"])
            self.assertEqual(plt.rcParams["figure.facecolor"], ECONOMIST.surface)
        with themed():
            self.assertFalse(plt.rcParams["axes.grid"])
            self.assertEqual(plt.rcParams["figure.facecolor"], METRO_ATLAS.surface)

    def test_monochrome_tier_tokens_and_colorblind_safe_grades(self) -> None:
        tok_econ = tier_tokens("economist")
        tok_light = tier_tokens("nature_editorial")
        tok_dark = tier_tokens("slide_dark")
        self.assertEqual([tok_econ[k]["dots"] for k in ("L3", "L2", "L1")], [3, 2, 1])
        self.assertEqual([tok_light[k]["dots"] for k in ("L3", "L2", "L1")], [3, 2, 1])
        self.assertEqual([tok_dark[k]["dots"] for k in ("L3", "L2", "L1")], [3, 2, 1])
        # 证据等级不包含高饱和红/蓝方向色，避免色相冲突
        for k in ("L3", "L2", "L1"):
            self.assertNotIn(tok_econ[k]["bg"], (ECONOMIST.up, ECONOMIST.down))
            self.assertNotIn(tok_light[k]["bg"], (NATURE_EDITORIAL.up, NATURE_EDITORIAL.down))
        # scoring_plots 不再使用红绿对比色 (#1B7837)
        self.assertNotEqual(EVIDENCE_GRADE_COLORS["Replicated"], "#1B7837")
        self.assertEqual(get_palette().up, METRO_ATLAS.up)
        self.assertEqual(get_palette("economist").up, ECONOMIST.up)
        self.assertEqual(get_palette("nature_editorial").up, NATURE_EDITORIAL.up)
        self.assertEqual(get_palette("poster_editorial"), get_preset("poster_editorial").palette)
        self.assertEqual(len(tint_color("#006BA2", 0.1)), 3)

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
