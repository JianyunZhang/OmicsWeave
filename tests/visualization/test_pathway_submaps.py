#!/usr/bin/env python3
"""代谢子通路专题图（pathway_submaps）预设目录、信号筛选与响应式渲染单元测试。"""
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
from omicsweave.integration.plotting.pathway_submaps import (
    ALL_SUBMAP_IDS,
    CORE_SUBMAP_IDS,
    SUBMAP_CATALOG,
    get_submap,
    select_submaps_for_view,
    submap_signal_stats,
)


def _synthetic_submap_view() -> tuple[pd.DataFrame, pd.DataFrame]:
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
                "flux_claim": "L2", "rank": 1, "key": "key", "value": 1.4,
                "display_state": "suggestive", "ring": "",
            })
        elif edge.region in ("glycerolipid", "fatty_acid", "ketone", "amino_acid"):
            rows.append({
                "element_type": "reaction", "element_id": fid, "region": edge.region, "label": fid,
                "draw_style": "omics_only", "reversible": False, "physio_sign": 1, "sign_change": False,
                "flux_group": "",
                "F_state": "not_measured", "F_direction": 0, "C_state": "changed", "C_direction": 1,
                "R_state": "unchanged", "R_direction": 0, "M_state": "not_measured", "M_direction": 0,
                "flux_claim": "L1", "rank": 5, "key": "l1", "value": np.nan,
                "display_state": "changed", "ring": "",
            })
        else:
            # urea_creatine / nucleotide / other 留空（not_assessed），验证自动跳过逻辑
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
                "draw_style": "", "value": 1.2 if sp in ("G3P_c", "CIT_m", "BHB_m", "GLU_c") else np.nan,
                "display_state": "changed" if sp in ("G3P_c", "CIT_m", "BHB_m", "GLU_c") else "not_measured",
                "ring": "selected" if sp in ("G3P_c", "CIT_m") else "",
            })
    for aid, agg in layout.aggregates.items():
        rows.append({
            "element_type": "aggregate", "element_id": aid, "region": agg.region, "label": agg.label or aid,
            "draw_style": "not_assessed",
            "F_state": "not_measured", "F_direction": 0,
            "C_state": "not_measured", "C_direction": 0,
            "R_state": "no_knowledge", "R_direction": 0, "M_state": "not_measured", "M_direction": 0,
            "flux_claim": "", "rank": np.nan, "key": "",
            "P_value": np.nan, "P_state": "not_measured",
            "R_value": np.nan, "R_layer_state": "not_measured",
            "Mx_value": np.nan, "Mx_state": "not_measured",
        })
    view = pd.DataFrame(rows)
    regions = pd.DataFrame([
        {"region": rid, "label": r.label, "kegg": ", ".join(r.kegg), "n_drawn": 2, "n_aggregates": 0,
         "n_flux": 0, "n_identifiable": 0, "n_omics_only": 0, "n_L2": 0, "n_L1": 0, "fusion": ""}
        for rid, r in layout.regions.items()
    ])
    return view, regions


class PathwaySubmapsTest(unittest.TestCase):
    def test_catalog_covers_layout_regions_cleanly(self) -> None:
        layout = default_layout()
        self.assertEqual(
            CORE_SUBMAP_IDS,
            ("central_carbon", "fatty_acid_ketone", "glycerolipid", "amino_acid", "urea_creatine"),
        )
        for spec in SUBMAP_CATALOG:
            for rid in spec.regions:
                self.assertIn(rid, layout.regions, f"{spec.id} references unknown region {rid}")
        with self.assertRaises(ValueError):
            get_submap("nonexistent_submap")

    def test_signal_aware_selection_skips_unassessed_panels(self) -> None:
        layout = default_layout()
        view, _ = _synthetic_submap_view()
        selected, skipped = select_submaps_for_view(view, layout, min_signal="assessed", include_optional_active=True)
        sel_ids = [s.id for s in selected]
        self.assertEqual(sel_ids, ["central_carbon", "fatty_acid_ketone", "glycerolipid", "amino_acid"])
        skip_map = dict(skipped)
        self.assertEqual(skip_map["urea_creatine"], "no_assessed_elements")
        self.assertEqual(skip_map["nucleotide"], "optional_inactive")
        self.assertEqual(skip_map["redox_one_carbon"], "optional_inactive")

        # 显式指定 all + min_signal="all" 时返回全部 7 张子图
        all_sel, all_skip = select_submaps_for_view(view, layout, panels=["all"], min_signal="all")
        self.assertEqual(tuple(s.id for s in all_sel), ALL_SUBMAP_IDS)
        self.assertEqual(all_skip, [])

        cc_stats = submap_signal_stats(view, get_submap("central_carbon"))
        self.assertEqual(cc_stats["n_l2"], 3)
        self.assertTrue(cc_stats["has_active"])

    def test_narrow_and_wide_submaps_render_across_styles(self) -> None:
        layout = default_layout()
        view, regions = _synthetic_submap_view()
        with tempfile.TemporaryDirectory() as tmp:
            # 覆盖极窄单栏（glycerolipid）、中等底部高 y0（urea_creatine）、紧裁线粒体（fatty_acid_ketone）与宽幅（central_carbon）
            cases = [
                ("glycerolipid", "metro_atlas", "semantic"),
                ("urea_creatine", "metro_atlas", "pathway"),
                ("fatty_acid_ketone", "metro_atlas", "semantic"),
                ("central_carbon", "nature_editorial", "semantic"),
            ]
            for sid, style, cm in cases:
                spec = get_submap(sid)
                paths = plot_pathway_map(
                    view, layout, regions,
                    f"Liver · HFD_vs_LFD — {spec.title}",
                    ["Run: test", "Calibration: test"],
                    Path(tmp) / f"{sid}_{style}_{cm}",
                    formats=("png",),
                    mode="publication",
                    style=style,
                    color_mode=cm,
                    focus_regions=list(spec.regions),
                    include_reactions=list(spec.include_reactions) or None,
                    exclude_reactions=list(spec.exclude_reactions) or None,
                    panel_note=f"Sub-pathway panel ({spec.subtitle}).",
                )
                self.assertEqual(len(paths), 1)
                self.assertGreater(paths[0].stat().st_size, 5000)


if __name__ == "__main__":
    unittest.main()
