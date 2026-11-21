#!/usr/bin/env python3
"""独立 Graphical Abstract 子系统测试：MAX_HUBS 护栏、YAML 规格、视图表聚合、三套渲染皮肤与确定性输出。"""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import pandas as pd

from omicsweave.visualization.graphical_abstract import (
    GA_STYLES,
    MAX_HUBS,
    AbstractFlow,
    AbstractHub,
    AbstractSpec,
    build_abstract_from_pathway_view,
    load_abstract_spec,
    render_graphical_abstract,
    validate_spec,
)


def _toy_pathway_view() -> pd.DataFrame:
    return pd.DataFrame([
        # TCA L2+ coupled flux up
        {"element_type": "reaction", "element_id": "CS_m", "region": "pyruvate_tca", "draw_style": "flux",
         "flux_claim": "L2", "flux_group": "ACONT_m", "F_state": "suggestive", "F_direction": 1,
         "C_state": "unchanged", "C_direction": 0, "value": 1.3},
        {"element_type": "reaction", "element_id": "ACONT_m", "region": "pyruvate_tca", "draw_style": "flux",
         "flux_claim": "L2", "flux_group": "ACONT_m", "F_state": "suggestive", "F_direction": 1,
         "C_state": "unchanged", "C_direction": 0, "value": 1.5},
        # Glycolysis L1 enzyme up
        {"element_type": "reaction", "element_id": "HEX_c", "region": "glycolysis_gng", "draw_style": "flux",
         "flux_claim": "L1", "flux_group": "", "F_state": "undetermined", "F_direction": 0,
         "C_state": "changed", "C_direction": 1, "value": float("nan")},
        # Fatty acid omics-only L1
        {"element_type": "aggregate", "element_id": "FAO_m", "region": "fatty_acid", "draw_style": "aggregate",
         "flux_claim": "L1", "flux_group": "", "F_state": "not_measured", "F_direction": 0,
         "C_state": "changed", "C_direction": 1, "value": float("nan")},
        # Metabolite pool
        {"element_type": "metabolite", "element_id": "GLC_c", "region": "", "draw_style": "",
         "flux_claim": "", "flux_group": "", "F_state": "", "F_direction": 0,
         "C_state": "", "C_direction": 0, "value": -1.2},
        {"element_type": "metabolite", "element_id": "BHB_m", "region": "", "draw_style": "",
         "flux_claim": "", "flux_group": "", "F_state": "", "F_direction": 0,
         "C_state": "", "C_direction": 0, "value": 1.1},
    ])


class GraphicalAbstractTest(unittest.TestCase):
    def test_max_hubs_and_integrity_guardrails(self) -> None:
        too_many = tuple(AbstractHub(id=f"h{i}", label=f"H{i}", xy=(10.0 + i, 20.0)) for i in range(MAX_HUBS + 1))
        with self.assertRaisesRegex(ValueError, str(MAX_HUBS)):
            validate_spec(AbstractSpec(title="too big", hubs=too_many))

        with self.assertRaises(ValueError):
            validate_spec(AbstractSpec(title="empty", hubs=()))

        h1 = AbstractHub(id="a", label="A", xy=(10.0, 10.0))
        with self.assertRaises(ValueError):
            validate_spec(AbstractSpec(title="dup", hubs=(h1, h1)))

        with self.assertRaises(ValueError):
            validate_spec(AbstractSpec(title="bad flow", hubs=(h1,), flows=(AbstractFlow("a", "missing"),)))

        with self.assertRaises(ValueError):
            validate_spec(AbstractSpec(title="bad icon", hubs=(AbstractHub("a", "A", (10.0, 10.0), icon="dragon"),)))

    def test_build_from_pathway_view_aggregates_honestly(self) -> None:
        spec = build_abstract_from_pathway_view(_toy_pathway_view(), tissue="liver", contrast="HFD_vs_LFD")
        self.assertLessEqual(len(spec.hubs), MAX_HUBS)
        by_id = {h.id: h for h in spec.hubs}
        self.assertIn("tca_cycle", by_id)
        self.assertEqual(by_id["tca_cycle"].evidence_tier, "L2")
        self.assertEqual(by_id["tca_cycle"].direction, 1)
        self.assertIn("2 Flux-Supported (L2+)", spec.takeaway)
        self.assertIn("2 Enzyme-Predicted (L1)", spec.takeaway)
        self.assertIn("Quantification:", spec.rules_note)
        self.assertIn("0.45", spec.rules_note)

    def test_load_example_yaml_and_render_all_three_styles_deterministically(self) -> None:
        example = Path(__file__).resolve().parents[2] / "examples" / "graphical_abstract.example.yaml"
        spec = load_abstract_spec(example)
        self.assertLessEqual(len(spec.hubs), MAX_HUBS)

        with tempfile.TemporaryDirectory() as tmp:
            for style in GA_STYLES:
                p1 = render_graphical_abstract(spec, Path(tmp) / f"run1_{style}", style=style, formats=("pdf", "png"))
                p2 = render_graphical_abstract(spec, Path(tmp) / f"run2_{style}", style=style, formats=("pdf", "png"))
                self.assertEqual([p.suffix for p in p1], [".pdf", ".png"])
                for a, b in zip(p1, p2):
                    self.assertGreater(a.stat().st_size, 2000)
                    self.assertEqual(a.read_bytes(), b.read_bytes(), f"{style} {a.suffix} not deterministic")

            # metro_atlas 额外支持可选 pathway 通路分区配色模式与纯矢量 SVG 输出
            (svg_path,) = render_graphical_abstract(
                spec, Path(tmp) / "metro_pathway", style="metro_atlas", color_mode="pathway", formats=("svg",),
            )
            self.assertGreater(svg_path.stat().st_size, 2000)
            self.assertNotIn("<image", svg_path.read_text(encoding="utf-8"))

        with self.assertRaises(ValueError):
            render_graphical_abstract(spec, Path("unused"), style="watercolor_neon")
        with self.assertRaises(ValueError):
            render_graphical_abstract(spec, Path("unused"), style="metro_atlas", color_mode="rainbow_neon")


if __name__ == "__main__":
    unittest.main()
