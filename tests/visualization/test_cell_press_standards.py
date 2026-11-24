#!/usr/bin/env python3
"""Cell Press 出版规范、解耦 JournalProfile、净版子通路图与 1200×1200 方形 GA 合规测试（零 pydantic 依赖）。"""
from __future__ import annotations

import tempfile
import unittest
from io import BytesIO
from pathlib import Path
from unittest.mock import patch

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from PIL import Image

from omicsweave.integration.mappings.network_layout import default_layout
from omicsweave.integration.plotting.pathway_map_html import write_pathway_map_html
from omicsweave.integration.plotting.pathway_map_plots import plot_pathway_map
from omicsweave.visualization.figures import (
    DOUBLE_COL_MM,
    MAX_PAGE_H_MM,
    ONE_HALF_COL_MM,
    SINGLE_COL_MM,
    journal_figsize,
    save_stem,
)
from omicsweave.visualization.graphical_abstract import (
    DEFAULT_GA_LAYOUT,
    GA_LAYOUTS,
    GA_STYLES,
    build_abstract_from_pathway_view,
    render_graphical_abstract,
)
from omicsweave.visualization.journal import (
    CELL_PRESS,
    DEFAULT_JOURNAL_PROFILE,
    NATURE,
    add_panel_label,
    audit_figure_compliance,
    get_journal_profile,
    list_journal_profiles,
)
from omicsweave.visualization.theme import themed
from tests.visualization.test_graphical_abstract import _toy_pathway_view
from tests.visualization.test_pathway_map_modes import _synthetic_map_view


class CellPressStandardsTest(unittest.TestCase):
    def test_journal_profiles_default_cell_press_and_switchable_nature(self) -> None:
        self.assertEqual(DEFAULT_JOURNAL_PROFILE, "cell_press")
        self.assertEqual(set(list_journal_profiles()), {"cell_press", "nature"})
        self.assertEqual(SINGLE_COL_MM, 85.0)
        self.assertEqual(ONE_HALF_COL_MM, 114.0)
        self.assertEqual(DOUBLE_COL_MM, 174.0)
        self.assertEqual(MAX_PAGE_H_MM, 200.0)

        w_cp, _ = journal_figsize(2, aspect=0.6)
        w_nat, _ = journal_figsize(2, aspect=0.6, profile="nature")
        self.assertAlmostEqual(w_cp * 25.4, 174.0, places=3)
        self.assertAlmostEqual(w_nat * 25.4, 183.0, places=3)
        self.assertAlmostEqual(journal_figsize(1.5, profile=CELL_PRESS)[0] * 25.4, 114.0, places=3)
        self.assertAlmostEqual(journal_figsize(1, profile=NATURE)[0] * 25.4, 89.0, places=3)

        with self.assertRaises(ValueError):
            get_journal_profile("unknown_journal")

    def test_themed_applies_cell_press_typography_embedding_and_tiff_eps_export(self) -> None:
        with themed(preset="metro_atlas", journal="cell_press"):
            self.assertEqual(plt.rcParams["pdf.fonttype"], 42)
            self.assertEqual(plt.rcParams["ps.fonttype"], 42)
            self.assertEqual(plt.rcParams["svg.fonttype"], "none")
            self.assertEqual(plt.rcParams["mathtext.fontset"], "stixsans")
            self.assertEqual(plt.rcParams["font.size"], 8.0)
            self.assertGreaterEqual(plt.rcParams["xtick.labelsize"], 6.0)

            fig, ax = plt.subplots(figsize=journal_figsize(1, aspect=0.8, profile="cell_press"))
            ax.plot([0, 1, 2], [0.2, 1.1, 0.8], lw=1.2)
            ax.set_xlabel(r"$\log_2\mathrm{FC}$")
            ax.set_ylabel("Intensity")
            lbl = add_panel_label(ax, "(a).", profile="cell_press")
            self.assertEqual(lbl.get_text(), "A")
            self.assertEqual(lbl.get_fontsize(), 10.0)

            report = audit_figure_compliance(fig, profile="cell_press")
            self.assertTrue(report.passed, f"Unexpected compliance errors: {report.issues}")
            self.assertAlmostEqual(report.width_mm, 85.0, places=1)
            self.assertGreaterEqual(report.min_font_pt or 0.0, 6.0)

            with tempfile.TemporaryDirectory() as tmp:
                out = save_stem(fig, Path(tmp) / "cp_export", formats=("pdf", "png", "tiff", "eps"))
                self.assertEqual([p.suffix for p in out], [".pdf", ".png", ".tiff", ".eps"])
                for p in out:
                    self.assertTrue(p.is_file() and p.stat().st_size > 1000)

        # Nature 面板字母规范为小写加粗
        with themed(preset="metro_atlas", journal="nature"):
            fig_n, ax_n = plt.subplots(figsize=journal_figsize(1, aspect=0.8, profile="nature"))
            lbl_n = add_panel_label(ax_n, "B", profile="nature")
            self.assertEqual(lbl_n.get_text(), "b")
            plt.close(fig_n)

    def test_compliance_auditor_flags_small_fonts_and_oversized_width(self) -> None:
        with themed(preset="metro_atlas", journal="cell_press"):
            fig, ax = plt.subplots(figsize=(9.5, 4.0))
            ax.text(0.5, 0.5, "too small", fontsize=4.2)
            report = audit_figure_compliance(fig, profile="cell_press")
            self.assertFalse(report.passed)
            rules = {i.rule for i in report.issues}
            self.assertIn("column_width", rules)
            self.assertIn("min_font_size", rules)
            plt.close(fig)

    def test_cell_press_square_graphical_abstract_default_is_exact_1200x1200_and_12pt_plus(self) -> None:
        self.assertEqual(DEFAULT_GA_LAYOUT, "cell_press_square")
        self.assertEqual(set(GA_LAYOUTS), {"cell_press_square", "widescreen"})

        spec = build_abstract_from_pathway_view(_toy_pathway_view(), tissue="liver", contrast="HFD_vs_LFD")
        captured_reports = []
        real_save_stem = save_stem

        def spy_save(fig, stem_path, **kwargs):
            captured_reports.append(audit_figure_compliance(fig, profile="cell_press", check_ga=True))
            return real_save_stem(fig, stem_path, **kwargs)

        with tempfile.TemporaryDirectory() as tmp:
            with patch("omicsweave.visualization.graphical_abstract.renderer.save_stem", side_effect=spy_save):
                paths = render_graphical_abstract(
                    spec, Path(tmp) / "ga_cp",
                    style="metro_atlas", formats=("pdf", "png", "tiff"),
                )
            with Image.open(BytesIO((Path(tmp) / "ga_cp.tiff").read_bytes())) as tiff:
                tiff.load()
                self.assertEqual(tiff.size, (1200, 1200))
                self.assertEqual(tiff.tag_v2[259], 5)  # TIFF compression tag 5 = LZW
            self.assertEqual([p.suffix for p in paths], [".pdf", ".png", ".tiff"])
            png_arr = mpimg.imread(paths[1])
            self.assertEqual(png_arr.shape[:2], (1200, 1200), f"Expected exact 1200x1200 px, got {png_arr.shape[:2]}")
            self.assertTrue((Path(tmp) / "ga_cp_caption.txt").is_file())
            self.assertEqual(len(captured_reports), 1)
            rep = captured_reports[0]
            self.assertTrue(rep.passed, f"GA failed Cell Press audit: {rep.issues}")
            self.assertGreaterEqual(rep.min_font_pt or 0.0, 12.0)

            # 所有 5 套皮肤均支持方形与宽屏切换
            for st in GA_STYLES:
                sq = render_graphical_abstract(spec, Path(tmp) / f"sq_{st}", style=st, layout="cell_press_square", formats=("png",))
                arr = mpimg.imread(sq[0])
                self.assertEqual(arr.shape[:2], (1200, 1200), f"Style {st} square GA must be 1200x1200")
            ws = render_graphical_abstract(spec, Path(tmp) / "ws_metro", style="metro_atlas", layout="widescreen", formats=("png",))
            ws_arr = mpimg.imread(ws[0])
            self.assertGreater(ws_arr.shape[1], ws_arr.shape[0])

    def test_publication_submap_enforces_clean_cell_press_width_and_6pt_floor(self) -> None:
        layout = default_layout()
        view, regions = _synthetic_map_view()
        captured = []
        real_save_stem = save_stem

        def spy_save(fig, stem_path, **kwargs):
            captured.append(audit_figure_compliance(fig, profile="cell_press"))
            return real_save_stem(fig, stem_path, **kwargs)

        with tempfile.TemporaryDirectory() as tmp:
            stem = Path(tmp) / "sub_central"
            with patch("omicsweave.integration.plotting.pathway_map_plots.save_stem", side_effect=spy_save):
                out = plot_pathway_map(
                    view, layout, regions,
                    "Liver · HFD_vs_LFD — Central Carbon",
                    ["Run: test_run", "Calibration: ok"],
                    stem,
                    formats=("pdf", "png"),
                    mode="publication",
                    style="metro_atlas",
                    journal="cell_press",
                    focus_regions=["glycolysis_gng", "pyruvate_tca"],
                    panel_label="A",
                )
            self.assertEqual([p.suffix for p in out], [".pdf", ".png"])
            cap_file = Path(tmp) / "sub_central_caption.txt"
            self.assertTrue(cap_file.is_file(), "Clean Cell Press publication panel must emit companion _caption.txt")
            cap_txt = cap_file.read_text(encoding="utf-8")
            self.assertIn("Central Carbon", cap_txt)
            self.assertIn("MFA Pooled IVW", cap_txt)

            self.assertEqual(len(captured), 1)
            rep = captured[0]
            self.assertTrue(rep.passed, f"Publication submap failed Cell Press audit: {rep.issues}")
            self.assertLessEqual(rep.width_mm, 174.5)
            self.assertGreaterEqual(rep.min_font_pt or 0.0, 6.0)

    def test_html_uses_arial_helvetica_stack(self) -> None:
        layout = default_layout()
        view, regions = _synthetic_map_view()
        with tempfile.TemporaryDirectory() as tmp:
            out = write_pathway_map_html(
                view, layout, regions, "Liver · HFD_vs_LFD", ["Run: test"],
                Path(tmp) / "map.html", style="metro_atlas",
            )
            txt = out.read_text(encoding="utf-8")
            self.assertIn('--font-sans: "Arial", "Helvetica Neue", Helvetica', txt)


if __name__ == "__main__":
    unittest.main()
