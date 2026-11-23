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

    def test_render_maps_omits_figure_e_label_in_outputs(self) -> None:
        """PDF/PNG/HTML 主图、子通路图及 Graphical Abstract 输出文本均不应标注 'Figure E'。"""
        import ast

        root = Path(__file__).resolve().parents[2]
        targets = [
            root / "src/omicsweave/integration/pathway_map.py",
            root / "src/omicsweave/integration/plotting/pathway_map_plots.py",
            root / "src/omicsweave/integration/plotting/pathway_map_html.py",
            root / "src/omicsweave/visualization/graphical_abstract/spec.py",
            root / "src/omicsweave/visualization/graphical_abstract/renderer.py",
        ]
        for path in targets:
            tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
            # Collect docstring AST nodes so only non-docstring runtime strings are checked
            doc_nodes: set[int] = set()
            for node in ast.walk(tree):
                if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                    if (
                        node.body
                        and isinstance(node.body[0], ast.Expr)
                        and isinstance(node.body[0].value, ast.Constant)
                        and isinstance(node.body[0].value.value, str)
                    ):
                        doc_nodes.add(id(node.body[0].value))
            for node in ast.walk(tree):
                if isinstance(node, ast.Constant) and isinstance(node.value, str) and id(node) not in doc_nodes:
                    self.assertNotIn(
                        "figure e",
                        node.value.lower(),
                        f"Unexpected 'Figure E' runtime string in {path.name}:{getattr(node, 'lineno', '?')}: {node.value!r}",
                    )

        # If full runtime dependencies (pydantic) are present, also exercise render_maps end-to-end
        try:
            from unittest.mock import patch
            from omicsweave.integration import pathway_map as pm
        except ModuleNotFoundError:
            return

        view, regions = _synthetic_map_view()
        members = pd.DataFrame([{
            "tissue": "liver", "contrast_id": "HFD_vs_LFD", "source": "U-13C-glucose",
            "layer": "mfa", "flux_id": "CS_m", "cell_state": "suggestive",
            "tracer": "glucose", "n_numerator": 6, "n_denominator": 6,
        }])
        axes = pd.DataFrame([{"tissue": "liver", "contrast_id": "HFD_vs_LFD", "flux_id": "CS_m"}])
        assessed = pd.DataFrame([{"tissue": "liver", "contrast_id": "HFD_vs_LFD", "flux_id": "CS_m"}])
        coverage = pd.DataFrame([{"flux_id": "CS_m", "placement": "edge"}])

        captured: list[str] = []

        def fake_plot(v, lay, reg, title, caption, stem, formats=("pdf", "png", "html"), **kwargs):
            captured.append(str(title))
            if kwargs.get("panel_note"):
                captured.append(str(kwargs["panel_note"]))
            return []

        def fake_ga(spec, stem, **kwargs):
            captured.append(str(spec.provenance))
            return []

        with tempfile.TemporaryDirectory() as tmp:
            dossier_dir = Path(tmp) / "dossier"
            out_dir = Path(tmp) / "out"
            (dossier_dir / "tables").mkdir(parents=True)
            with (
                patch.object(pm, "coverage_table", return_value=coverage),
                patch.object(pm, "map_view", return_value=view),
                patch.object(pm, "region_summary", return_value=regions),
                patch.object(pm, "extra_pathway_fusion", return_value=pd.DataFrame()),
                patch.object(pm, "signed_module_summary", return_value={}),
                patch("omicsweave.integration.plotting.pathway_map_plots.plot_pathway_map", side_effect=fake_plot),
                patch("omicsweave.visualization.graphical_abstract.render_graphical_abstract", side_effect=fake_ga),
            ):
                pm.render_maps(
                    dossier_dir, out_dir, {"formats": ["pdf", "png", "html"]},
                    source="U-13C-glucose", runs={"Integration": "run_test"},
                    members=members, axes=axes, assessed=assessed, calibration=pd.DataFrame(),
                )

        self.assertGreaterEqual(len(captured), 1)
        self.assertTrue(captured[0].startswith("Liver · HFD_vs_LFD"))
        for text in captured:
            self.assertNotIn("figure e", text.lower())

    def test_svg_status_cell_distinguishes_unchanged_and_undetermined(self) -> None:
        from omicsweave.integration.plotting.pathway_map_html import _svg_status_cell
        from omicsweave.visualization.palettes import get_palette

        pal = get_palette("metro_atlas")
        svg_unch = _svg_status_cell(0.0, 0.0, 1.0, "unchanged", 0, pal)
        svg_undet = _svg_status_cell(0.0, 0.0, 1.0, "undetermined", 0, pal)
        self.assertEqual(svg_unch.count("<line "), 1)
        self.assertEqual(svg_undet.count("<line "), 2)
        self.assertIn(pal.ink2, svg_unch)
        self.assertIn(pal.muted, svg_undet)

    def test_flux_gate_checklist_payload_and_drawer_js(self) -> None:
        layout = default_layout()
        view, regions = _synthetic_map_view()
        view["F_flags"] = ""
        view["F_g"] = np.nan
        view["text"] = ""
        hex_mask = (view["element_type"] == "reaction") & (view["element_id"] == "HEX_c")
        view.loc[hex_mask, "F_flags"] = "single_subject_driven"
        view.loc[hex_mask, "F_g"] = 0.66
        view.loc[hex_mask, "text"] = "可辨识，未达 suggestive/changed（|g| = 0.66 < 1，单只动物驱动）"

        members = pd.DataFrame([
            {
                "flux_id": "HEX_c", "tissue": "liver", "contrast_id": "HFD_vs_LFD", "source": "s1",
                "layer": "mfa", "role": "flux", "member_id": "HEX_c", "feature_id": "HEX_c",
                "feature_name": "HEX_c", "tracer": "pooled", "member_relation": "", "regulator_effect": "",
                "log2fc": 0.85, "effect": 0.66, "effect_type": "hedges_g", "smd_se": 0.42, "padj": 0.22,
                "padj_global": 0.42, "fdr_scope": "per_tissue", "pvalue": 0.08,
                "lfsr": np.nan, "abundance_share": np.nan, "mean_numerator": 18.0, "mean_denominator": 10.0,
                "n_numerator": 6, "n_denominator": 6, "cell_state": "no_evidence", "state_reason": "",
                "identifiability": "well_constrained", "heterogeneity_pvalue": 0.61, "direction_stability": 0.72,
                "pooled_loo_stability": 0.83, "smd_weights": "glucose:0.55;lactate:0.45",
                "tracer_weights": "glucose:0.90;lactate:0.10", "g_total_var": 1.42, "smd_weighting": "smd",
                "flags": "single_subject_driven",
            },
        ])
        with tempfile.TemporaryDirectory() as tmp:
            out = write_pathway_map_html(
                view, layout, regions, "Gate Checklist Test", ["Run: demo"],
                Path(tmp) / "gate.html",
                members=members, tissue="liver", contrast_id="HFD_vs_LFD", source="s1",
            )
            txt = out.read_text(encoding="utf-8")
            self.assertIn("renderFluxGateChecklist", txt)
            self.assertIn("Flux Assessment Gate Checklist", txt)
            self.assertIn("Unchanged (Equiv.)", txt)
            self.assertIn("PASS (P+G)", txt)
            payload = _extract_payload_json(txt)
            hx = payload["elements"]["rxn:HEX_c"]
            self.assertEqual(hx["f_flags"], ["single_subject_driven"])
            self.assertAlmostEqual(hx["f_g"], 0.66)
            self.assertIn("|g| = 0.66 < 1", hx["f_rule"])
            m0 = hx["members"][0]
            self.assertAlmostEqual(m0["effect"], 0.66)
            self.assertAlmostEqual(m0["smd_se"], 0.42)
            self.assertAlmostEqual(m0["padj_global"], 0.42)
            self.assertEqual(m0["fdr_scope"], "per_tissue")
            self.assertAlmostEqual(m0["dir_stab"], 0.72)
            self.assertAlmostEqual(m0["loo_stab"], 0.83)
            self.assertAlmostEqual(m0["g_total_var"], 1.42)
            self.assertEqual(m0["smd_weights"], "glucose:0.55;lactate:0.45")
            self.assertEqual(m0["tracer_weights"], "glucose:0.90;lactate:0.10")


if __name__ == "__main__":
    unittest.main()
