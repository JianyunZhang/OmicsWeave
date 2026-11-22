#!/usr/bin/env python3
"""标准化 16 种同位素示踪剂注册表、动态示踪剂证据徽章与多示踪剂 Pooled IVW / Linked flux 溯源测试。"""
from __future__ import annotations

import json
import re
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from omicsweave.integration.mappings.network_layout import default_layout
from omicsweave.integration.plotting.pathway_map_html import write_pathway_map_html
from omicsweave.integration.plotting.pathway_map_plots import (
    _independent_flux_chains,
    _strip_badges,
    plot_pathway_map,
)
from omicsweave.mfa.tracers import (
    TRACER_SPECS,
    canonical_tracer_id,
    format_tracer_badge_label,
    format_tracer_display,
    format_tracer_short,
    list_tracers,
    resolve_tracer,
)
from tests.visualization.test_pathway_map_modes import _synthetic_map_view


class TracersAndBadgesTest(unittest.TestCase):
    def test_all_16_isotope_tracers_registered_and_backward_compatible(self) -> None:
        specs = list_tracers()
        self.assertEqual(len(specs), 16)
        self.assertEqual(len(TRACER_SPECS), 16)

        expected_display = {
            # 4 existing U-13C tracers
            "U-13C-glucose",
            "U-13C-glycerol",
            "U-13C-alanine",
            "U-13C-lactate",
            # 12 future tracers
            "U-13C-palmitic acid",
            "U-15N-alanine",
            "U-13C-15N-alanine",
            "U-13C-glutamine",
            "U-15N-glutamine",
            "U-13C-15N-glutamine",
            "U-13C-glycine",
            "U-15N-glycine",
            "U-13C-15N-glycine",
            "U-13C-leucine",
            "U-15N-leucine",
            "U-13C-15N-leucine",
        }
        self.assertEqual({s.display_name for s in specs}, expected_display)

        # Legacy bare names resolve to U-13C display names seamlessly
        self.assertEqual(format_tracer_display("glucose"), "U-13C-glucose")
        self.assertEqual(format_tracer_display("glycerol"), "U-13C-glycerol")
        self.assertEqual(format_tracer_display("alanine"), "U-13C-alanine")
        self.assertEqual(format_tracer_display("lactate"), "U-13C-lactate")
        self.assertEqual(format_tracer_display("alanine_13c15n"), "U-13C-15N-alanine")
        self.assertEqual(format_tracer_display("palmitate"), "U-13C-palmitic acid")
        self.assertEqual(format_tracer_display("U-13C-palmitic acid"), "U-13C-palmitic acid")
        self.assertEqual(format_tracer_display("15n_glutamine"), "U-15N-glutamine")
        self.assertEqual(format_tracer_display("glutamine_13c15n"), "U-13C-15N-glutamine")
        self.assertEqual(format_tracer_display("glycine_15n"), "U-15N-glycine")
        self.assertEqual(format_tracer_display("leucine_13c15n"), "U-13C-15N-leucine")
        self.assertEqual(format_tracer_display("pooled"), "Pooled (IVW)")

        # Distinct isotope forms NEVER collapse into the same canonical ID
        ala_ids = {
            canonical_tracer_id("U-13C-alanine"),
            canonical_tracer_id("U-15N-alanine"),
            canonical_tracer_id("U-13C-15N-alanine"),
        }
        self.assertEqual(len(ala_ids), 3)
        self.assertEqual(resolve_tracer("U-15N-glutamine").tracked_atoms, ("N",))
        self.assertEqual(resolve_tracer("U-13C-15N-leucine").tracked_atoms, ("C", "N"))
        self.assertEqual(resolve_tracer("U-13C-palmitic acid").entry_node, "PAL_e")
        self.assertEqual(format_tracer_short("U-13C-palmitic acid"), "U-13C-Pal")

    def test_dynamic_tracer_badge_formatting(self) -> None:
        # L2 single tracer → "<short_label> Flux"
        self.assertEqual(format_tracer_badge_label("U-13C-glucose", level="L2"), "U-13C-Glc Flux")
        self.assertEqual(format_tracer_badge_label("palmitate", level="L2"), "U-13C-Pal Flux")
        self.assertEqual(format_tracer_badge_label("U-15N-glutamine", level="L2"), "U-15N-Gln Flux")
        # L2 multi-tracer same isotope vs mixed isotope
        self.assertEqual(format_tracer_badge_label(["glucose", "alanine"], level="L2"), "U-13C-Glc+Ala Flux")
        self.assertEqual(format_tracer_badge_label(["U-13C-alanine", "U-15N-alanine"], level="L2"), "U-13C-Ala+U-15N-Ala Flux")
        # L3 multi-tracer displays exact tracer names
        self.assertEqual(format_tracer_badge_label("U-13C-glucose;U-13C-alanine", level="L3"), "U-13C-Glc+Ala")
        self.assertEqual(format_tracer_badge_label(["U-13C-glutamine", "U-13C-15N-glutamine"], level="L3"), "U-13C-Gln+U-13C15N-Gln")
        self.assertEqual(format_tracer_badge_label(["glucose", "alanine", "lactate"], level="L3"), "U-13C-Glc+Ala+Lac")
        # Fallbacks & L1 unchanged
        self.assertEqual(format_tracer_badge_label("", level="L2"), "Flux")
        self.assertEqual(format_tracer_badge_label("", level="L3"), "Confirmed")
        self.assertEqual(format_tracer_badge_label("glucose", level="L1"), "Enz-pred")

    def test_pathway_map_pdf_png_and_html_render_tracer_badges_and_shared_dof_rules(self) -> None:
        layout = default_layout()
        view, regions = _synthetic_map_view()
        # Attach explicit tracer provenance: ACONT_m & ICDH_m share 1 locked chain (U-13C-glucose, L2),
        # and promote CS_m to L3 multi-tracer (U-13C-glucose + U-13C-alanine) to test both badge styles!
        view["flux_tracers"] = ""
        view["flux_pooling_mode"] = ""
        view.loc[view["element_id"].isin(["ACONT_m", "ICDH_m"]), "flux_tracers"] = "U-13C-glucose"
        view.loc[view["element_id"].isin(["ACONT_m", "ICDH_m"]), "flux_pooling_mode"] = "single_tracer"
        cs_mask = view["element_id"] == "CS_m"
        view.loc[cs_mask, ["draw_style", "F_state", "F_direction", "flux_claim", "claim_direction", "value", "key"]] = [
            "flux", "changed", 1, "L3", 1, 1.15, "key",
        ]
        view.loc[cs_mask, "flux_tracers"] = "U-13C-glucose;U-13C-alanine"
        view.loc[cs_mask, "flux_pooling_mode"] = "multi_tracer_pooled"

        # 3 L2/L3 reactions (CS_m + {ACONT_m, ICDH_m}) = 2 independent flux chains!
        self.assertEqual(_independent_flux_chains(view), 2)

        acont_row = view[view["element_id"] == "ACONT_m"].iloc[0]
        cs_row = view[view["element_id"] == "CS_m"].iloc[0]
        self.assertEqual(_strip_badges(acont_row, {"ACONT_m": 1})["tier_label"], "U-13C-Glc Flux")
        self.assertEqual(_strip_badges(cs_row, {"ACONT_m": 1})["tier_label"], "U-13C-Glc+Ala")

        with tempfile.TemporaryDirectory() as tmp:
            out_files = plot_pathway_map(
                view, layout, regions, "Tracer Badge & Pooled Rules Test", ["Run: demo"],
                Path(tmp) / "map_tracer", formats=("png", "html"), mode="atlas",
            )
            self.assertEqual(len(out_files), 2)
            html_txt = (Path(tmp) / "map_tracer.html").read_text(encoding="utf-8")
            self.assertIn("U-13C-Glc Flux", html_txt)
            self.assertIn("U-13C-Glc+Ala", html_txt)
            self.assertIn("Shared 1 DOF", html_txt)
            self.assertIn("NOT extra confidence", html_txt)
            m = re.search(r'<script id="ow-payload" type="application/json">(.*?)</script>', html_txt, flags=re.S)
            payload = json.loads(m.group(1).replace("<\\/", "</"))
            self.assertEqual(payload["meta"]["counts"]["indep_flux_chains"], 2)
            self.assertEqual(payload["elements"]["rxn:ACONT_m"]["flux_tracers"], ["U-13C-glucose"])
            self.assertEqual(payload["elements"]["rxn:ACONT_m"]["flux_pooling_mode"], "single_tracer")
            self.assertEqual(payload["elements"]["rxn:CS_m"]["flux_tracers"], ["U-13C-glucose", "U-13C-alanine"])
            self.assertEqual(payload["elements"]["rxn:CS_m"]["flux_pooling_mode"], "multi_tracer_pooled")


if __name__ == "__main__":
    unittest.main()
