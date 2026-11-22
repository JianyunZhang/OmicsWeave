#!/usr/bin/env python3
"""方案 A 单页交互式 HTML 通路图测试：离线自包含、语义 SVG、经济学人主题与逐只动物散点数据。"""
from __future__ import annotations

import json
import re
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from omicsweave.integration.mappings.network_layout import default_layout
from omicsweave.integration.plotting.pathway_map_html import write_pathway_map_html
from tests.visualization.test_pathway_map_modes import _synthetic_map_view


def _extract_payload_json(html_text: str) -> dict:
    m = re.search(r'<script id="ow-payload" type="application/json">(.*?)</script>', html_text, flags=re.S)
    assert m is not None, "Missing embedded #ow-payload JSON block"
    return json.loads(m.group(1).replace("<\\/", "</"))


class PathwayMapHtmlTest(unittest.TestCase):
    def test_self_contained_economist_html_with_members_and_per_animal_samples(self) -> None:
        layout = default_layout()
        view, regions = _synthetic_map_view()

        members = pd.DataFrame([
            {
                "flux_id": "ACONT_m", "tissue": "liver", "contrast_id": "HFD_vs_LFD", "source": "s1",
                "layer": "mfa", "role": "flux", "member_id": "ACONT_m", "feature_id": "ACONT_m",
                "feature_name": "ACONT_m", "tracer": "glucose", "member_relation": "", "regulator_effect": "",
                "log2fc": 1.4, "padj": 0.03, "pvalue": 0.008, "lfsr": 0.02, "abundance_share": np.nan,
                "mean_numerator": 42.0, "mean_denominator": 16.0, "n_numerator": 3, "n_denominator": 3,
                "cell_state": "suggestive", "state_reason": "", "identifiability": "well_constrained",
            },
            {
                "flux_id": "HEX_c", "tissue": "liver", "contrast_id": "HFD_vs_LFD", "source": "s1",
                "layer": "proteomics", "role": "direct_enzyme", "member_id": "Hk1", "feature_id": "Hk1",
                "feature_name": "Hexokinase 1", "tracer": "", "member_relation": "isozyme", "regulator_effect": "",
                "log2fc": 1.25, "padj": 0.01, "pvalue": 0.002, "lfsr": 0.005, "abundance_share": 0.72,
                "mean_numerator": 12.1, "mean_denominator": 10.85, "n_numerator": 3, "n_denominator": 3,
                "cell_state": "changed", "state_reason": "", "identifiability": "",
            },
        ])
        samples = pd.DataFrame([
            # MFA per-animal with fit IQR & one bound-pinned animal
            {"source": "s1", "layer": "mfa", "tissue": "liver", "feature_id": "ACONT_m", "tracer": "glucose",
             "sample_id": f"mfa_{i}", "subject_id": f"m{i}", "group": "HFD" if i < 3 else ("LFD" if i < 6 else "HSD"),
             "group_order": 1 if i < 3 else (0 if i < 6 else 2),
             "value": (40.0 + i) if i < 3 else (15.0 + i), "scale": "linear", "units": "nmol/min/g",
             "at_bound": (i == 0), "fit_low": 35.0 if i < 3 else 12.0, "fit_high": 45.0 if i < 3 else 20.0,
             "qc_flag": ""}
            for i in range(8)
        ] + [
            {"source": "s1", "layer": "proteomics", "tissue": "liver", "feature_id": "Hk1", "tracer": "",
             "sample_id": f"pro_{i}", "subject_id": f"m{i}", "group": "HFD" if i < 3 else "LFD",
             "group_order": 1 if i < 3 else 0, "value": (12.0 + 0.1 * i) if i < 3 else (10.7 + 0.1 * i),
             "scale": "log2", "units": "", "at_bound": False, "fit_low": np.nan, "fit_high": np.nan, "qc_flag": ""}
            for i in range(6)
        ])
        contrasts = pd.DataFrame([
            {"source": "s1", "layer": "mfa", "contrast_id": "HFD_vs_LFD", "sample_id": f"mfa_{i}",
             "side": "numerator" if i < 3 else "denominator", "level": "HFD" if i < 3 else "LFD"}
            for i in range(6)
        ] + [
            {"source": "s1", "layer": "proteomics", "contrast_id": "HFD_vs_LFD", "sample_id": f"pro_{i}",
             "side": "numerator" if i < 3 else "denominator", "level": "HFD" if i < 3 else "LFD"}
            for i in range(6)
        ])

        with tempfile.TemporaryDirectory() as tmp:
            p1 = write_pathway_map_html(
                view, layout, regions, "Liver · HFD vs LFD", ["Run: demo", "Cal: ok", "Layout: default"],
                Path(tmp) / "map_a.html",
                style="economist", members=members, samples=samples, sample_contrasts=contrasts,
                tissue="liver", contrast_id="HFD_vs_LFD", source="s1",
            )
            p2 = write_pathway_map_html(
                view, layout, regions, "Liver · HFD vs LFD", ["Run: demo", "Cal: ok", "Layout: default"],
                Path(tmp) / "map_b.html",
                style="economist", members=members, samples=samples, sample_contrasts=contrasts,
                tissue="liver", contrast_id="HFD_vs_LFD", source="s1",
            )
            txt = p1.read_text(encoding="utf-8")
            # Reproducible deterministic single-file offline HTML (no external CDN scripts/styles)
            self.assertEqual(p1.read_bytes(), p2.read_bytes())
            self.assertNotIn("<script src=", txt)
            self.assertNotIn('<link rel="stylesheet"', txt)
            self.assertNotIn("https://", txt)
            self.assertIn('data-theme="economist"', txt)
            self.assertIn('data-rxn="ACONT_m"', txt)
            self.assertIn('data-flux-group="ACONT_m"', txt)
            self.assertIn('data-node="GAP_c"', txt)
            self.assertIn('data-agg="FAO_m"', txt)

            payload = _extract_payload_json(txt)
            self.assertEqual(payload["meta"]["initial_style"], "economist")
            self.assertEqual(payload["meta"]["num_label"], "HFD")
            self.assertEqual(payload["meta"]["den_label"], "LFD")
            self.assertTrue(payload["meta"]["has_samples"])
            self.assertEqual(set(payload["palettes"]), {"economist", "metro_atlas", "nature_editorial", "slide_dark"})
            self.assertEqual(payload["palettes"]["economist"]["surface"], "#EBEBDF")
            self.assertEqual(payload["palettes"]["economist"]["numerator"], "#E3120B")
            self.assertEqual(payload["palettes"]["economist"]["denominator"], "#006BA2")

            # Exact parity with map_view + embedded per-animal series
            acont = payload["elements"]["rxn:ACONT_m"]
            self.assertAlmostEqual(acont["value"], 1.4)
            self.assertEqual(acont["flux_group"], "ACONT_m")
            self.assertEqual(len(acont["members"]), 1)
            skey = acont["members"][0]["skey"]
            self.assertIn(skey, payload["samples"])
            ser = payload["samples"][skey]
            self.assertEqual(ser["scale"], "linear")
            self.assertEqual(ser["groups"], ["LFD", "HFD", "HSD"])
            self.assertEqual(len(ser["pts"]), 8)
            # One bound-pinned animal recorded
            self.assertEqual(sum(pt[7] for pt in ser["pts"]), 1)

    def test_html_defaults_to_metro_atlas_with_delicate_pdf_borders(self) -> None:
        layout = default_layout()
        view, regions = _synthetic_map_view()
        with tempfile.TemporaryDirectory() as tmp:
            p = write_pathway_map_html(
                view, layout, regions, "Metro Atlas Default Test", ["Run: demo"],
                Path(tmp) / "metro_default.html",
            )
            txt = p.read_text(encoding="utf-8")
            payload = _extract_payload_json(txt)
            self.assertEqual(payload["meta"]["initial_style"], "metro_atlas")
            self.assertIn('data-theme="metro_atlas"', txt)
            self.assertFalse(payload["meta"]["has_samples"])
            self.assertEqual(payload["samples"], {})
            # Metro Atlas PDF visual parity: region dots, P/R/M squares on aggregates, delicate thin borders
            self.assertIn("ow-region-dot", txt)
            self.assertIn("ow-prm-sq", txt)
            self.assertIn("fitSvgLabels", txt)
            card_strokes = [
                float(m) for m in re.findall(r'class="ow-card-bg"[^>]*stroke-width="([0-9.]+)"', txt)
            ]
            self.assertTrue(card_strokes, "Expected rendered reaction/aggregate card borders")
            self.assertLessEqual(max(card_strokes), 0.18)
            reg_strokes = [
                float(m) for m in re.findall(r'class="ow-region-box[^"]*"[^>]*stroke-width="([0-9.]+)"', txt)
            ]
            self.assertTrue(reg_strokes)
            self.assertLessEqual(max(reg_strokes), 0.10)


if __name__ == "__main__":
    unittest.main()
