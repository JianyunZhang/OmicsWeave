#!/usr/bin/env python3
"""反应档案图 A/B/C/D 独立视觉渲染测试（验证统一矢量图标、三点证据徽章与纯矢量图例）。"""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pandas as pd

from omicsweave.integration.dossier.panels import Panel
from omicsweave.integration.plotting.dossier_plots import (
    plot_dossier_page,
    plot_matrix,
    plot_node_split,
    plot_tracer_coverage,
)


def _toy_panel() -> Panel:
    return Panel(
        id="central",
        title="Central Carbon Panel",
        segments=[("Glycolysis", ["HEX_c", "PYK_c"]), ("TCA", ["ACONT_m"])],
        expected_tracers={"HEX_c": ["glucose"], "ACONT_m": ["glucose"]},
    )


class DossierPlotsTest(unittest.TestCase):
    def test_figures_a_b_c_d_render_with_vector_icons_and_badges(self) -> None:
        panel = _toy_panel()

        # Figure A matrix_view covering every state & tier
        mrows = [
            # HEX_c: flux unidentifiable, enzyme changed, L1 Enz-pred
            {"row": 0, "element_type": "reaction", "element_id": "HEX_c", "label": "HEX_c", "segment": "Glycolysis",
             "block": "flux", "column": "glucose", "source": "s1", "layer": "mfa", "tracer": "glucose",
             "expected_tracer": True, "value": np.nan, "certainty": np.nan, "display_state": "not_identifiable",
             "ring": "", "n_changed": 0, "n_quantified": 0, "text": ""},
            {"row": 0, "element_type": "reaction", "element_id": "HEX_c", "label": "HEX_c", "segment": "Glycolysis",
             "block": "enzyme", "column": "P:s1", "source": "s1", "layer": "proteomics", "tracer": "",
             "expected_tracer": False, "value": -1.1, "certainty": 0.95, "display_state": "summary",
             "ring": "selected", "n_changed": 2, "n_quantified": 3, "text": ""},
            {"row": 0, "element_type": "reaction", "element_id": "HEX_c", "label": "HEX_c", "segment": "Glycolysis",
             "block": "enzyme", "column": "R:s1", "source": "s1", "layer": "transcriptomics", "tracer": "",
             "expected_tracer": False, "value": np.nan, "certainty": np.nan, "display_state": "not_detected",
             "ring": "", "n_changed": 0, "n_quantified": 0, "text": ""},
            {"row": 0, "element_type": "reaction", "element_id": "HEX_c", "label": "HEX_c", "segment": "Glycolysis",
             "block": "assessment", "column": "F", "source": "s1", "layer": "", "tracer": "",
             "expected_tracer": False, "value": 0.0, "certainty": np.nan, "display_state": "not_identifiable",
             "ring": "", "n_changed": 0, "n_quantified": 0, "text": ""},
            {"row": 0, "element_type": "reaction", "element_id": "HEX_c", "label": "HEX_c", "segment": "Glycolysis",
             "block": "assessment", "column": "C", "source": "s1", "layer": "", "tracer": "",
             "expected_tracer": False, "value": -1.0, "certainty": np.nan, "display_state": "changed",
             "ring": "", "n_changed": 0, "n_quantified": 0, "text": ""},
            {"row": 0, "element_type": "reaction", "element_id": "HEX_c", "label": "HEX_c", "segment": "Glycolysis",
             "block": "assessment", "column": "level", "source": "s1", "layer": "", "tracer": "",
             "expected_tracer": False, "value": -1.0, "certainty": np.nan, "display_state": "L1",
             "ring": "", "n_changed": 0, "n_quantified": 0, "text": "L1 (1)"},
            # ACONT_m: L2 Flux + mixed + QC
            {"row": 1, "element_type": "reaction", "element_id": "ACONT_m", "label": "ACONT_m", "segment": "TCA",
             "block": "flux", "column": "glucose", "source": "s1", "layer": "mfa", "tracer": "glucose",
             "expected_tracer": True, "value": 1.35, "certainty": 0.92, "display_state": "suggestive",
             "ring": "changed", "n_changed": 0, "n_quantified": 0, "text": ""},
            {"row": 1, "element_type": "reaction", "element_id": "ACONT_m", "label": "ACONT_m", "segment": "TCA",
             "block": "assessment", "column": "F", "source": "s1", "layer": "", "tracer": "",
             "expected_tracer": False, "value": 1.0, "certainty": np.nan, "display_state": "suggestive",
             "ring": "", "n_changed": 0, "n_quantified": 0, "text": ""},
            {"row": 1, "element_type": "reaction", "element_id": "ACONT_m", "label": "ACONT_m", "segment": "TCA",
             "block": "assessment", "column": "C", "source": "s1", "layer": "", "tracer": "",
             "expected_tracer": False, "value": 0.0, "certainty": np.nan, "display_state": "mixed",
             "ring": "", "n_changed": 0, "n_quantified": 0, "text": ""},
            {"row": 1, "element_type": "reaction", "element_id": "ACONT_m", "label": "ACONT_m", "segment": "TCA",
             "block": "assessment", "column": "level", "source": "s1", "layer": "", "tracer": "",
             "expected_tracer": False, "value": 1.0, "certainty": np.nan, "display_state": "L2",
             "ring": "", "n_changed": 0, "n_quantified": 0, "text": "L2 (2)"},
        ]
        mview = pd.DataFrame(mrows)

        # Figure B & C rows
        members = pd.DataFrame([
            {"tissue": "liver", "contrast_id": "HFD_vs_LFD", "source": "s1", "flux_id": "HEX_c",
             "layer": "mfa", "tracer": "glucose", "role": "flux", "member_id": "HEX_c", "feature_id": "HEX_c",
             "feature_name": "HEX_c", "abundance_share": np.nan, "cell_state": "not_identifiable",
             "state_reason": "bound pinned", "identifiability": "bound_pinned", "n_numerator": 3, "n_denominator": 3,
             "log2fc": np.nan, "pvalue": np.nan, "padj": np.nan, "regulator_effect": "",
             "direction_evidence": "forward"},
            {"tissue": "liver", "contrast_id": "HFD_vs_LFD", "source": "s1", "flux_id": "ACONT_m",
             "layer": "mfa", "tracer": "glucose", "role": "flux", "member_id": "ACONT_m", "feature_id": "ACONT_m",
             "feature_name": "ACONT_m", "abundance_share": np.nan, "cell_state": "insufficient_n",
             "state_reason": "n < 2", "identifiability": "moderate", "n_numerator": 1, "n_denominator": 2,
             "log2fc": 0.8, "pvalue": np.nan, "padj": np.nan, "regulator_effect": "",
             "direction_evidence": "forward"},
            {"tissue": "liver", "contrast_id": "HFD_vs_LFD", "source": "s1", "flux_id": "HEX_c",
             "layer": "proteomics", "tracer": "", "role": "direct_enzyme", "member_id": "Hk1", "feature_id": "Hk1",
             "feature_name": "Hk1", "abundance_share": 0.6, "cell_state": "changed",
             "state_reason": "", "identifiability": "", "n_numerator": 3, "n_denominator": 3,
             "log2fc": -1.1, "pvalue": 0.002, "padj": 0.01, "regulator_effect": "",
             "direction_evidence": "forward"},
        ])
        samples = pd.DataFrame([
            {"source": "s1", "layer": "proteomics", "tissue": "liver", "feature_id": "Hk1", "tracer": "",
             "sample_id": f"m{i}", "group": "HFD" if i < 3 else "LFD", "group_order": 1 if i < 3 else 0,
             "value": 10.0 - (1.1 if i < 3 else 0.0) + 0.1 * i, "scale": "log2", "units": ""}
            for i in range(6)
        ])
        sides = pd.DataFrame([
            {"contrast_id": "HFD_vs_LFD", "source": "s1", "layer": "proteomics",
             "sample_id": f"m{i}", "side": "numerator" if i < 3 else "denominator"}
            for i in range(6)
        ])

        # Figure D splits (one computed, one unidentifiable)
        splits = pd.DataFrame([
            {"node": "PYR_m", "exit_a": "PDH_m", "exit_b": "PC_m", "source": "s1", "layer": "mfa", "tracer": "glucose",
             "state": "not_identifiable", "log_ratio_change": np.nan, "pvalue": np.nan, "padj": np.nan,
             "members_a": "PDH_m", "members_b": "PC_m", "detail": "bound pinned"},
            {"node": "PYR_m", "exit_a": "PDH_m", "exit_b": "PC_m", "source": "s1", "layer": "proteomics", "tracer": "",
             "state": "computed", "log_ratio_change": -0.75, "pvalue": 0.02, "padj": 0.04,
             "members_a": "Pdha1", "members_b": "Pcx", "detail": ""},
        ])
        per_anim = pd.DataFrame([
            {"node": "PYR_m", "exit_a": "PDH_m", "exit_b": "PC_m", "source": "s1", "layer": "proteomics", "tracer": "",
             "sample_id": f"m{i}", "group": "HFD" if i < 3 else "LFD", "group_order": 1 if i < 3 else 0,
             "value": (0.2 if i < 3 else 0.95) + 0.05 * i}
            for i in range(6)
        ])

        with tempfile.TemporaryDirectory() as tmp:
            p_a = plot_matrix(mview, panel, "liver", "HFD_vs_LFD", Path(tmp) / "figA", formats=("png",), caption="test")
            p_b = plot_dossier_page(members[members["flux_id"] == "HEX_c"], samples, sides, "HEX_c", "liver", "HFD_vs_LFD", Path(tmp) / "figB", formats=("png",))
            p_c = plot_tracer_coverage(members, panel, "liver", "HFD_vs_LFD", Path(tmp) / "figC", formats=("png",), tracers=("glucose",))
            p_d = plot_node_split(splits, per_anim, sides, panel, "liver", "HFD_vs_LFD", Path(tmp) / "figD", formats=("png",))
            for paths in (p_a, p_b, p_c, p_d):
                self.assertEqual(len(paths), 1)
                self.assertTrue(paths[0].is_file())
                self.assertGreater(paths[0].stat().st_size, 3000)


if __name__ == "__main__":
    unittest.main()
