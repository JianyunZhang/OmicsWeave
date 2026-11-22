#!/usr/bin/env python3
"""图 E 双模式（atlas / publication）与多主题预设独立渲染测试（零 pydantic / seaborn 依赖）。"""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pandas as pd

from omicsweave.integration.mappings.network_layout import default_layout
from omicsweave.integration.plotting.pathway_map_plots import plot_pathway_map


def _synthetic_map_view() -> tuple[pd.DataFrame, pd.DataFrame]:
    layout = default_layout()
    rows = []
    for fid, edge in layout.edges.items():
        if fid in ("CS_m", "ACONT_m", "ICDH_m"):
            rows.append({
                "element_type": "reaction", "element_id": fid, "region": edge.region, "label": fid,
                "draw_style": "flux", "reversible": False, "physio_sign": 1, "sign_change": False,
                "flux_group": "ACONT_m" if fid != "CS_m" else "",
                "F_state": "suggestive", "F_direction": 1, "C_state": "unchanged", "C_direction": 0,
                "R_state": "no_knowledge", "R_direction": 0, "M_state": "undetermined", "M_direction": 0,
                "flux_claim": "L2", "rank": 1 if fid == "ACONT_m" else 2, "key": "key", "value": 1.4,
                "display_state": "suggestive", "ring": "",
            })
        elif fid == "HEX_c":
            rows.append({
                "element_type": "reaction", "element_id": fid, "region": edge.region, "label": fid,
                "draw_style": "flux", "reversible": False, "physio_sign": 1, "sign_change": False,
                "flux_group": "",
                "F_state": "undetermined", "F_direction": 0, "C_state": "changed", "C_direction": 1,
                "R_state": "mixed", "R_direction": 0, "M_state": "not_measured", "M_direction": 0,
                "flux_claim": "L1", "rank": 4, "key": "l1", "value": np.nan,
                "display_state": "undetermined", "ring": "",
            })
        else:
            rows.append({
                "element_type": "reaction", "element_id": fid, "region": edge.region, "label": fid,
                "draw_style": "not_assessed", "reversible": False, "physio_sign": 0, "sign_change": False,
                "flux_group": "",
                "F_state": "not_measured", "F_direction": 0, "C_state": "no_knowledge", "C_direction": 0,
                "R_state": "no_knowledge", "R_direction": 0, "M_state": "no_knowledge", "M_direction": 0,
                "flux_claim": "", "rank": np.nan, "key": "", "value": np.nan,
                "display_state": "not_measured", "ring": "",
            })
    for sp, node in layout.nodes.items():
        if node.kind == "metabolite":
            rows.append({
                "element_type": "metabolite", "element_id": sp, "region": "", "label": node.label or sp,
                "draw_style": "", "value": 1.5 if sp == "GAP_c" else (-1.2 if sp == "GLC_c" else np.nan),
                "display_state": "changed" if sp in ("GAP_c", "GLC_c") else "not_measured",
                "ring": "selected" if sp == "GAP_c" else "",
            })
    for aid, agg in layout.aggregates.items():
        rows.append({
            "element_type": "aggregate", "element_id": aid, "region": agg.region, "label": agg.label or aid,
            "draw_style": "aggregate" if aid == "FAO_m" else "not_assessed",
            "F_state": "not_measured", "F_direction": 0,
            "C_state": "changed" if aid == "FAO_m" else "not_measured", "C_direction": 1 if aid == "FAO_m" else 0,
            "R_state": "no_knowledge", "R_direction": 0, "M_state": "not_measured", "M_direction": 0,
            "flux_claim": "L1" if aid == "FAO_m" else "", "rank": np.nan, "key": "l1" if aid == "FAO_m" else "",
            "P_value": 1.1 if aid == "FAO_m" else np.nan, "P_state": "changed" if aid == "FAO_m" else "not_measured",
            "R_value": np.nan, "R_layer_state": "not_measured",
            "Mx_value": np.nan, "Mx_state": "not_measured",
        })
    view = pd.DataFrame(rows)
    regions = pd.DataFrame([
        {"region": rid, "label": r.label, "kegg": ", ".join(r.kegg), "n_drawn": 2, "n_aggregates": 0,
         "n_flux": 1 if rid in ("pyruvate_tca", "glycolysis_gng") else 0, "n_identifiable": 1,
         "n_omics_only": 0, "n_L2": 3 if rid == "pyruvate_tca" else 0, "n_L1": 1 if rid == "glycolysis_gng" else 0,
         "fusion": ""}
        for rid, r in layout.regions.items()
    ])
    return view, regions


class PathwayMapModesTest(unittest.TestCase):
    def test_publication_and_atlas_modes_render_with_presets(self) -> None:
        layout = default_layout()
        view, regions = _synthetic_map_view()
        with tempfile.TemporaryDirectory() as tmp:
            p_pub = plot_pathway_map(
                view, layout, regions, "Publication Test", ["Line 1"], Path(tmp) / "pub",
                formats=("png",), mode="publication", style="nature_editorial",
            )
            p_dark = plot_pathway_map(
                view, layout, regions, "Dark Slide Test", ["Line 1"], Path(tmp) / "dark",
                formats=("png",), mode="publication", style="slide_dark",
                focus_regions=["glycolysis_gng", "pyruvate_tca"],
            )
            p_atlas = plot_pathway_map(
                view, layout, regions, "Atlas Test", ["Line 1"], Path(tmp) / "atlas",
                formats=("png",), mode="atlas", style="nature_editorial",
            )
            p_metro_sem = plot_pathway_map(
                view, layout, regions, "Metro Semantic", ["Line 1"], Path(tmp) / "metro_sem",
                formats=("png",), mode="publication", style="metro_atlas", color_mode="semantic",
            )
            p_metro_pw = plot_pathway_map(
                view, layout, regions, "Metro Pathway", ["Line 1"], Path(tmp) / "metro_pw",
                formats=("png",), mode="atlas", style="metro_atlas", color_mode="pathway",
            )
            p_axes = plot_pathway_map(
                view, layout, regions, "Metro Axes", ["Line 1"], Path(tmp) / "metro_axes",
                formats=("png",), mode="publication", style="metro_atlas", strip_mode="axes",
                focus_regions=["glycolysis_gng", "pyruvate_tca"],
            )
            p_econ_atlas = plot_pathway_map(
                view, layout, regions, "Economist Atlas", ["Line 1"], Path(tmp) / "econ_atlas",
                formats=("png", "html"), mode="atlas", style="economist",
            )
            self.assertEqual([p.suffix for p in p_econ_atlas], [".png", ".html"])
            self.assertTrue(all(p.is_file() and p.stat().st_size > 5000 for p in p_econ_atlas))

            p_econ_pub = plot_pathway_map(
                view, layout, regions, "Economist Publication", ["Line 1"], Path(tmp) / "econ_pub",
                formats=("png",), mode="publication", style="economist",
                focus_regions=["glycolysis_gng", "pyruvate_tca"],
            )
            p_grey = plot_pathway_map(
                view, layout, regions, "Metro Grey Omics", ["Line 1"], Path(tmp) / "metro_grey",
                formats=("png",), mode="publication", style="metro_atlas", omics_edge_color="grey",
                focus_regions=["glycolysis_gng", "pyruvate_tca"],
            )
            for paths in (p_pub, p_dark, p_atlas, p_metro_sem, p_metro_pw, p_axes, p_econ_pub, p_grey):
                self.assertEqual(len(paths), 1)
                self.assertTrue(paths[0].is_file())
                self.assertGreater(paths[0].stat().st_size, 5000)

        with self.assertRaises(ValueError):
            plot_pathway_map(view, layout, regions, "Bad", [], Path("unused"), mode="cartoon")
        with self.assertRaises(ValueError):
            plot_pathway_map(view, layout, regions, "Bad", [], Path("unused"), style="metro_atlas", color_mode="neon")
        with self.assertRaises(ValueError):
            plot_pathway_map(view, layout, regions, "Bad", [], Path("unused"), strip_mode="invalid")
        with self.assertRaises(ValueError):
            plot_pathway_map(view, layout, regions, "Bad", [], Path("unused"), omics_edge_color="rainbow")


if __name__ == "__main__":
    unittest.main()
