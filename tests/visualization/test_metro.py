#!/usr/bin/env python3
"""通用 Beck 地铁拓扑原语（omicsweave.visualization.metro）与莫兰迪分区色单元测试。"""
from __future__ import annotations

import math
import tempfile
import unittest
from pathlib import Path
from xml.etree import ElementTree

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from omicsweave.visualization.metro import (
    draw_metro_coupled_halo,
    draw_metro_loop,
    draw_metro_station,
    draw_metro_track,
    draw_metro_zone,
    fillet_polyline,
    polyline_length,
    route_octilinear,
    sample_polyline_at,
    trim_polyline,
)
from omicsweave.visualization.palettes import (
    METRO_PATHWAY_COLORS,
    VALID_COLOR_MODES,
    pathway_family_color,
)


def _is_octilinear_angle(dx: float, dy: float, tol: float = 1e-2) -> bool:
    ang = abs(math.degrees(math.atan2(dy, dx))) % 45.0
    return min(ang, 45.0 - ang) < tol


class MetroPrimitivesTest(unittest.TestCase):
    def test_route_octilinear_angles_and_parallel_bends(self) -> None:
        # 一般斜向两点：每段必须严格为 0°/45°/90°
        for pref in ("hv_first", "diag_first", "centered"):
            poly = route_octilinear((0.0, 0.0), (14.0, 5.0), prefer=pref)
            self.assertGreaterEqual(len(poly), 3)
            for a, b in zip(poly, poly[1:]):
                self.assertTrue(_is_octilinear_angle(b[0] - a[0], b[1] - a[1]), f"{pref}: {a}->{b}")

        # 正逆平行偏置应分列中轴两侧
        pos = route_octilinear((10.0, 0.0), (10.0, 12.0), bend=1.6)
        neg = route_octilinear((10.0, 0.0), (10.0, 12.0), bend=-1.6)
        self.assertGreater(len(pos), 2)
        self.assertGreater(len(neg), 2)
        pos_mid_x = sum(p[0] for p in pos[1:-1]) / len(pos[1:-1])
        neg_mid_x = sum(p[0] for p in neg[1:-1]) / len(neg[1:-1])
        self.assertNotEqual(math.copysign(1.0, pos_mid_x - 10.0), math.copysign(1.0, neg_mid_x - 10.0))

    def test_fillet_and_trim_preserve_geometry(self) -> None:
        raw = route_octilinear((0.0, 0.0), (12.0, 6.0), prefer="hv_first")
        smooth = fillet_polyline(raw, radius=1.2, arc_steps=6)
        self.assertGreater(len(smooth), len(raw))
        total = polyline_length(smooth)
        trimmed = trim_polyline(smooth, trim0=1.5, trim1=1.5)
        self.assertAlmostEqual(polyline_length(trimmed), total - 3.0, delta=0.15)
        pt, tan = sample_polyline_at(trimmed, polyline_length(trimmed) * 0.5)
        self.assertAlmostEqual(math.hypot(*tan), 1.0, places=5)
        self.assertTrue(math.isfinite(pt[0]) and math.isfinite(pt[1]))

    def test_metro_drawing_primitives_render_pure_vectors(self) -> None:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_xlim(0, 80)
        ax.set_ylim(50, 0)
        try:
            draw_metro_zone(ax, (4, 4, 76, 46), "Mitochondrial Matrix", label_pos="bottom_right")
            p_solid = route_octilinear((10, 14), (34, 24), prefer="centered")
            draw_metro_coupled_halo(ax, [fillet_polyline(p_solid, radius=1.5)])
            draw_metro_track(ax, p_solid, "#C83E37", width=3.8, hollow=False, unit_in=0.075)
            p_hollow = route_octilinear((10, 34), (34, 24), prefer="hv_first")
            draw_metro_track(ax, p_hollow, "#1D5FA7", width=3.4, hollow=True, chevron_bidirectional=True, unit_in=0.075)
            draw_metro_loop(ax, (54, 24), 22, 13, "#C83E37", label="TCA Cycle", subtitle="Flux ↑ · 3 L2+", unit_in=0.075)
            draw_metro_station(ax, (10, 14), 2.2, fill="#FDE8E8", edge="#C83E37", outer_ring="#111827")
            draw_metro_station(ax, (10, 34), 2.2, fill="#F1F5F9", edge="#64748B", slash_color="#CBD5E1")

            with tempfile.TemporaryDirectory() as tmp:
                svg = Path(tmp) / "metro_primitives.svg"
                fig.savefig(svg, format="svg")
                root = ElementTree.parse(svg)
                ns = {"svg": "http://www.w3.org/2000/svg"}
                self.assertIsNotNone(root.find(".//svg:path", ns))
                self.assertIsNone(root.find(".//svg:image", ns))
        finally:
            plt.close(fig)

    def test_pathway_family_color_mapping(self) -> None:
        self.assertEqual(VALID_COLOR_MODES, ("semantic", "pathway"))
        self.assertEqual(pathway_family_color("glycolysis_gng"), METRO_PATHWAY_COLORS["carbohydrate"])
        self.assertEqual(pathway_family_color("tca_cycle"), METRO_PATHWAY_COLORS["respiration"])
        self.assertEqual(pathway_family_color("fatty_acid"), METRO_PATHWAY_COLORS["lipid"])
        self.assertEqual(pathway_family_color("unknown_custom", fallback_accent="#334455"), "#334455")


if __name__ == "__main__":
    unittest.main()
