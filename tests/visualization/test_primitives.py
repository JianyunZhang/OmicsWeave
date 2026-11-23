#!/usr/bin/env python3
"""可复用科学矢量图元测试：三点阶梯徽章、左上角排名角标、耦合轨道底带、降噪四轴微条与贝塞尔箭头。"""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from omicsweave.visualization.palettes import get_palette
from omicsweave.visualization.primitives import (
    OMICS_LAYER_ICONS,
    TIER_DISPLAY_LABELS,
    coupled_badge_label,
    coupled_badge_width,
    draw_axis_strip,
    draw_coupled_badge,
    draw_coupled_track,
    draw_layer_icon,
    draw_organelle_box,
    draw_quad_arrow,
    draw_rank_tab,
    draw_sign_reversal_badge,
    draw_status_cell,
    draw_tier_badge,
    draw_zone_card,
    draw_zone_header,
    tier_badge_text,
    tier_badge_width,
)
from omicsweave.visualization.theme import themed


class PrimitivesTest(unittest.TestCase):
    def test_badge_labels_are_self_explanatory(self) -> None:
        self.assertEqual(TIER_DISPLAY_LABELS, {"L3": "Confirmed", "L2": "Flux", "L1": "Enz-pred"})
        self.assertEqual(tier_badge_text("L3"), "Confirmed")
        self.assertEqual(tier_badge_text("L2"), "Flux")
        self.assertEqual(tier_badge_text("L1"), "Enz-pred")
        self.assertEqual(tier_badge_text("L2", label_style="code"), "L2")
        self.assertEqual(tier_badge_text("L0"), "")
        self.assertEqual(coupled_badge_label(1, total_groups=1), "Linked flux")
        self.assertEqual(coupled_badge_label(2, total_groups=3), "Linked #2")
        self.assertEqual(coupled_badge_label(1, label="Custom"), "Custom")

    def test_badges_and_strips_render_cleanly_across_all_states(self) -> None:
        pal = get_palette("nature_editorial")
        with themed(preset="nature_editorial"):
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.set_xlim(0, 40)
            ax.set_ylim(30, 0)
            ax.set_aspect("equal")

            draw_zone_card(ax, (1, 1, 39, 29), "#35749e", palette=pal)
            draw_zone_header(ax, (2, 3), "Glycolysis", "#35749e", palette=pal)
            draw_organelle_box(ax, (16, 12, 37, 27), "Mitochondrial Matrix", "m", palette=pal, label_pos="bottom_left")

            curves = [{"p0": (4.0, 10.0), "c": (10.0, 8.0), "p1": (16.0, 10.0), "trim0": 1.0, "trim1": 1.0}]
            draw_coupled_track(ax, curves, color=pal.coupled_track, edge_color=pal.coupled_edge)
            draw_quad_arrow(ax, (4.0, 10.0), (10.0, 8.0), (16.0, 10.0), pal.up, 1.8, "-", False, True, 1.0, 1.0, 3.0, halo=True)

            # Tier badges: descriptive labels scale dynamically with text width
            self.assertEqual(draw_tier_badge(ax, 2, 6, "L0", 5.0), 0.0)
            self.assertEqual(tier_badge_width("L0"), 0.0)
            w3 = draw_tier_badge(ax, 2, 6, "L3", 5.0)
            w2 = draw_tier_badge(ax, 2 + w3 + 1, 6, "L2", 5.0)
            w1 = draw_tier_badge(ax, 2 + w3 + w2 + 2, 6, "L1", 5.0)
            self.assertGreater(w3, w1)
            self.assertGreater(w1, w2)
            self.assertGreater(w2, 0)
            # Legacy code style preserves uniform 2-char widths
            self.assertAlmostEqual(
                tier_badge_width("L3", label_style="code"),
                tier_badge_width("L2", label_style="code"),
            )

            # Rank tab, Linked flux badge & net-flux reversal amber badge
            rt_w = draw_rank_tab(ax, 20, 6, 1, 6.0, palette=pal)
            cb_w = draw_coupled_badge(ax, 24, 6, 1, 6.0, palette=pal, total_groups=1)
            cb_w2 = draw_coupled_badge(ax, 24, 9, 2, 6.0, palette=pal, total_groups=3)
            draw_sign_reversal_badge(ax, 18.0, 10.0, radius=0.8, palette=pal, zorder=6.0)
            self.assertGreater(rt_w, 0)
            self.assertAlmostEqual(cb_w, coupled_badge_width(1, total_groups=1))
            self.assertAlmostEqual(cb_w2, coupled_badge_width(2, total_groups=3))

            # Pure-vector layer icons & standalone status cells
            for idx, kind in enumerate(("flux", "protein", "mrna", "metabolite", "regulation")):
                draw_layer_icon(ax, kind, 3.0 + idx * 2.2, 13.5, size=0.9, color=pal.ink, zorder=5.0)
            for idx, (st, dr) in enumerate((
                ("changed", 1), ("changed", -1), ("suggestive", 1), ("suggestive", -1),
                ("mixed", 0), ("unchanged", 0), ("undetermined", 0), ("not_measured", 0),
                ("not_identifiable", 0), ("not_detected", 0), ("qc_blocked", 0),
            )):
                draw_status_cell(ax, 2.0 + idx * 1.4, 15.0, 1.1, st, dr, palette=pal, zorder=5.0)

            # Denoised 4-omics icon strip covering every state
            w_s, h_s = draw_axis_strip(
                ax, 3, 18,
                [("changed", 1), ("suggestive", -1), ("mixed", 0), ("unchanged", 0)],
                header_icons=OMICS_LAYER_ICONS,
                palette=pal,
            )
            self.assertGreater(w_s, 0)
            self.assertGreater(h_s, 0)
            draw_axis_strip(
                ax, 3, 23,
                [("not_measured", 0), ("no_knowledge", 0), ("not_identifiable", 0), ("not_detected", 0)],
                show_header_letters=False, palette=pal,
            )
            draw_axis_strip(
                ax, 18, 18,
                [("insufficient_n", 0), ("undetermined", 0), ("qc_blocked", 0), ("changed", -1)],
                palette=pal,
            )

            # Verify y_up=True + non-square y_scale path used by Figure A/B/C/D matrices
            fig2, ax2 = plt.subplots(figsize=(6, 4))
            ax2.set_xlim(0, 12)
            ax2.set_ylim(0, 10)
            for idx, kind in enumerate(("flux", "protein", "mrna", "metabolite", "regulation")):
                draw_layer_icon(ax2, kind, 1.0 + idx * 1.5, 8.5, size=0.5, color=pal.ink, y_scale=1.25, clip_on=False)
            for idx, (st, dr) in enumerate((
                ("changed", 1), ("suggestive", -1), ("mixed", 0), ("not_identifiable", 0),
                ("not_detected", 0), ("qc_blocked", 0), ("not_measured", 0),
            )):
                draw_status_cell(
                    ax2, 0.8 + idx * 1.2, 5.0, 0.7, st, dr,
                    palette=pal, draw_box=True, y_up=True, y_scale=1.25,
                )
            bw_up = draw_tier_badge(ax2, 1.0, 2.5, "L3", 5.0, scale=0.40, unit_in=0.50, fontsize=5.5, y_scale=1.25)
            self.assertGreater(bw_up, 0)

            with tempfile.TemporaryDirectory() as tmp:
                out = Path(tmp) / "primitives.png"
                out2 = Path(tmp) / "primitives_yup.png"
                fig.savefig(out, dpi=150)
                fig2.savefig(out2, dpi=150)
                self.assertTrue(out.is_file())
                self.assertGreater(out.stat().st_size, 1000)
                self.assertTrue(out2.is_file())
                self.assertGreater(out2.stat().st_size, 1000)
            plt.close(fig)
            plt.close(fig2)

    def test_unchanged_and_undetermined_render_distinct_glyphs(self) -> None:
        pal = get_palette("nature_editorial")
        fig_u, ax_u = plt.subplots()
        draw_status_cell(ax_u, 0.0, 0.0, 1.0, "unchanged", 0, palette=pal)
        fig_d, ax_d = plt.subplots()
        draw_status_cell(ax_d, 0.0, 0.0, 1.0, "undetermined", 0, palette=pal)
        # unchanged draws 1 solid bar; undetermined draws 2 dashed segments in muted color
        self.assertEqual(len(ax_u.lines), 1)
        self.assertEqual(len(ax_d.lines), 2)
        self.assertEqual(ax_u.lines[0].get_color(), pal.ink2)
        self.assertTrue(all(ln.get_color() == pal.muted for ln in ax_d.lines))
        plt.close(fig_u)
        plt.close(fig_d)


if __name__ == "__main__":
    unittest.main()
