#!/usr/bin/env python3
"""通路总图（图 E，G18）：别名、布局、视图与出图。

验收（docs/arch_design/omicsweave_pathway_map.md 第 3、4、7 节）：
- 登记为别名的旧命名条目不进入评估集，MFA 结果用到的 ID 除外；
- 包内布局与网络表一致：每个网络反应都画出或写明排除原因，坐标齐全，校验值与网络表相同；
- 图只画视图表：代谢物节点的值与图 A 的池格子相同；箭头方向按两组组均值的符号，变号时标出；
- 成员全部属于画出的网络反应的聚合节点不另画方框；
- reaction_dossier 运行时出图并写出视图表、去向表，交付目录复制图件。
"""
from __future__ import annotations

import contextlib
import io
import re
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import load_workbook

from omicsweave import cli
from omicsweave.core.config import load_experiment, load_project
from omicsweave.core.delivery import Delivery
from omicsweave.integration.dossier.map_view import (
    coverage_table,
    extra_pathway_fusion,
    map_view,
    physiological_sign,
    signed_module_summary,
)
from omicsweave.integration.dossier.panels import load_panels, panel_elements
from omicsweave.integration.dossier.view import matrix_view, node_members
from omicsweave.integration.mappings import FLUX_GENE_MAPPING, FLUX_METABOLITE_MAPPING, REACTION_ENZYME_MAPPING
from omicsweave.integration.mappings.network_annotation import (
    curated_reactions,
    network_reactions,
    pathway_reactions,
    reaction_semantics,
)
from omicsweave.integration.mappings.network_layout import (
    LayoutError, default_layout, load_layout, network_table_sha256, validate_layout)
from omicsweave.integration.mappings.reaction_aliases import assessment_universe, reaction_aliases, validate_aliases
from omicsweave.integration.pathway_map import settings
from tests.integration import test_dossier_panels as _panels
from tests.integration.test_dossier_panels import PANEL, write_panels
from tests.integration.test_flux_anchored import DIET
from tests.mfa.test_mfa_runner import EXP


class AliasTest(unittest.TestCase):
    def test_aliases_are_consistent(self) -> None:
        mapping = set(FLUX_GENE_MAPPING) | set(REACTION_ENZYME_MAPPING) | set(FLUX_METABOLITE_MAPPING)
        self.assertEqual(validate_aliases(network_reactions(), mapping), [])
        self.assertEqual(reaction_aliases()["PK_c"], ("PYK_c",))

    def test_universe_drops_aliases_unless_measured(self) -> None:
        universe, excluded = assessment_universe({"PK_c", "PYK_c", "GAPDH_c", "FAO_m"}, {"GAPDH_c"})
        self.assertEqual(universe, ["FAO_m", "GAPDH_c", "PYK_c"])     # GAPDH_c 是 MFA 结果用到的 ID，保留
        self.assertEqual(excluded, ["PK_c"])

    def test_new_network_members(self) -> None:
        """M0b：原先没有成员的网络反应已登记直接酶或载体。"""
        self.assertIn("Cpt1a", FLUX_GENE_MAPPING["CPT1_c"]["gene_roles"]["catalytic"])
        self.assertEqual(REACTION_ENZYME_MAPPING["CACT_PALCAR"]["transporters"], ["Slc25a20"])
        self.assertEqual(FLUX_METABOLITE_MAPPING["PNP_INO_c"]["products"], ["Hypoxanthine"])
        for flux_id in ("GLND_c", "GLUD_c"):                         # 与碳网络同一种酶，不重复登记
            self.assertNotIn(flux_id, FLUX_GENE_MAPPING)


class LayoutTest(unittest.TestCase):
    def test_package_layout_matches_network(self) -> None:
        layout = default_layout()
        self.assertEqual(validate_layout(layout), [])
        self.assertEqual(layout.network_sha256, network_table_sha256())

    def test_validation_reports_problems(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            tmp = Path(folder) / "layout.yaml"
            lines = default_layout().path.read_text(encoding="utf-8").replace("\r\n", "\n").splitlines()
            lines = [line for line in lines if not line.startswith(("  HEX_c:", "  PYR_m:"))]
            text = "\n".join(lines) + "\n"
            tmp.write_text(text.replace("network_sha256: ", "network_sha256: x"), encoding="utf-8")
            problems = validate_layout(load_layout(tmp))
        joined = "；".join(problems)
        self.assertIn("HEX_c 既不在 edges 中", joined)
        self.assertIn("物种 PYR_m 没有坐标", joined)
        self.assertIn("network_sha256", joined)

    def test_clones_pick_the_nearest_position(self) -> None:
        layout = default_layout()
        geometry = layout.geometry(network_reactions()["GLUD_m"])
        self.assertEqual(dict(geometry["substrates"])["AKG_m"], layout.nodes["AKG_m"].at[1])   # 左侧的克隆，不是 TCA 环上的

    def test_secondary_species_are_notes(self) -> None:
        """转氨的谷氨酸 / α-酮戊二酸不画线：GPT 只连丙酮酸与丙氨酸，次要物种单列。"""
        geometry = default_layout().geometry(network_reactions()["GPT_c"])
        self.assertEqual([s for s, _ in geometry["substrates"]], ["PYR_c"])
        self.assertEqual([s for s, _ in geometry["products"]], ["ALA_c"])
        self.assertEqual(geometry["secondary"], (["GLU_c"], ["AKG_c"]))

    def test_unused_clone_is_reported(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            tmp = Path(folder) / "layout.yaml"
            layout = default_layout()
            gx, gy = layout.nodes["GLN_c"].at[0]
            lines = layout.path.read_text(encoding="utf-8").replace("\r\n", "\n").splitlines()
            lines = [f"  GLN_c: {{at: [[{gx}, {gy}], [5, 99]]}}" if line.startswith("  GLN_c:") else line for line in lines]
            tmp.write_text("\n".join(lines) + "\n", encoding="utf-8")
            problems = validate_layout(load_layout(tmp), check_sha=False)
        self.assertEqual(problems, ["nodes.GLN_c 的位置 [5.0, 99.0] 没有反应连到"])

    def test_settings(self) -> None:
        self.assertEqual(settings({"highlight": {"top_rank": 5}})["highlight"], {"min_claim": "L2", "top_rank": 5})
        with self.assertRaises(ValueError):
            settings({"colour": "red"})
        with self.assertRaises(ValueError):
            settings({"highlight": {"min_claim": "L0"}})


class SignTest(unittest.TestCase):
    def test_physiological_sign(self) -> None:
        rows = lambda *pairs: pd.DataFrame([{"tracer": t, "mean_numerator": a, "mean_denominator": b}
                                            for t, a, b in pairs])
        self.assertEqual(physiological_sign(rows(("pooled", -84.7, -131.9))), (-1, False))   # 反向运行
        self.assertEqual(physiological_sign(rows(("pooled", -28.5, 8.0))), (0, True))        # 变号
        self.assertEqual(physiological_sign(rows(("glucose", 3, 2), ("alanine", 1, 1))), (1, False))
        self.assertEqual(physiological_sign(rows(("glucose", np.nan, np.nan))), (0, False))


class CoverageTest(unittest.TestCase):
    def test_aggregate_covered_by_drawn_steps(self) -> None:
        layout = default_layout()
        members = pd.DataFrame([
            {"flux_id": "TRIOSE_c", "role": "direct_enzyme", "layer": "transcriptomics", "member_id": "Tpi1",
             "cell_state": "no_evidence"},
            {"flux_id": "TPI_c", "role": "direct_enzyme", "layer": "transcriptomics", "member_id": "Tpi1",
             "cell_state": "no_evidence"},
            {"flux_id": "FAO_m", "role": "direct_enzyme", "layer": "transcriptomics", "member_id": "Acox1",
             "cell_state": "no_evidence"},
            {"flux_id": "UREA_cycle", "role": "direct_enzyme", "layer": "transcriptomics", "member_id": "",
             "cell_state": "no_knowledge"},
        ])
        table = coverage_table(layout, members, ["TRIOSE_c", "TPI_c", "FAO_m", "UREA_cycle", "PK_c", "XYZ_c"])
        placement = dict(zip(table["flux_id"], table["placement"]))
        self.assertEqual(placement, {"TRIOSE_c": "covered_by_steps", "TPI_c": "edge", "FAO_m": "aggregate",
                                     "UREA_cycle": "aggregate", "PK_c": "alias", "XYZ_c": "unplaced"})

    def test_layout_pure_english_and_expanded_modules(self) -> None:
        layout = default_layout()
        cjk = re.compile(r"[\u3400-\u4dbf\u4e00-\u9fff\u3000-\u303f\uff00-\uffef]")
        texts = (
            [r.label for r in layout.regions.values()]
            + [str(c.get("label", "")) for c in layout.compartments.values()]
            + [n.label for n in layout.nodes.values() if n.label]
            + [a.label for a in layout.aggregates.values() if a.label]
        )
        self.assertEqual([t for t in texts if cjk.search(t)], [])
        for mod in ("OXPHOS_CI", "OXPHOS_CIII", "OXPHOS_CIV", "OXPHOS_CV", "THERMOGENESIS_UCP",
                    "HIS_metabolism", "PHE_TYR_metabolism", "HYP_metabolism", "BALA_metabolism",
                    "INO_metabolism", "HEME_synthesis"):
            self.assertIn(mod, layout.aggregates)
            self.assertIn(mod, FLUX_GENE_MAPPING)
            self.assertIn(mod, REACTION_ENZYME_MAPPING)
        for edge_id in ("ME_m", "SUC_trans", "FUM_trans", "MAL_trans", "NOS_c", "ORNT_m",
                        "POLYAMINE_synthesis", "FRUCTOSE_METABOLISM", "HEXOSAMINE_HBP",
                        "BCAA_trans", "PROPIONATE_ANAPLEROSIS", "HMGCR", "BILE_acid",
                        "PUFA_EICOSANOID", "PHOSPHOLIPID", "SPHINGOLIPID", "PEROX_FAO",
                        "ETHER_LIPID", "FOLATE_cycle", "MET_cycle", "GSH_synthesis",
                        "PYROGLU", "TAU_metabolism", "ADK_c", "ALLANTOIN",
                        "TRP_metabolism", "NAD_SALVAGE"):
            self.assertIn(edge_id, layout.edges)
            self.assertIn(edge_id, curated_reactions())
            self.assertIn(edge_id, FLUX_GENE_MAPPING)
            self.assertIn(edge_id, REACTION_ENZYME_MAPPING)
            self.assertTrue(curated_reactions()[edge_id].kegg_pathway, edge_id)
            self.assertTrue(curated_reactions()[edge_id].reactome_id, edge_id)

    def test_dual_layer_semantics_and_mfa_precedence(self) -> None:
        mfa = network_reactions()
        cur = curated_reactions()
        merged = pathway_reactions()
        self.assertTrue(set(cur).isdisjoint(set(mfa)))
        self.assertEqual(set(merged), set(mfa) | set(cur))
        # Curated transport / near-equilibrium steps carry capacity semantics; irreversible steps are directional
        for fid in ("SUC_trans", "FUM_trans", "MAL_trans", "ORNT_m", "ADK_c", "FOLATE_cycle"):
            self.assertEqual(reaction_semantics(fid).direction_evidence, "capacity", fid)
        for fid in ("ME_m", "NOS_c", "HMGCR", "PROPIONATE_ANAPLEROSIS", "HEXOSAMINE_HBP"):
            self.assertEqual(reaction_semantics(fid).direction_evidence, "directional", fid)
        # MDH_m must not hijack mitochondrial dicarboxylate carriers (Slc25a10/11)
        self.assertTrue({"Slc25a10", "Slc25a11"}.isdisjoint(REACTION_ENZYME_MAPPING["MDH_m"].get("transporters", [])))

    def test_sidebar_p3_helpers_and_derived_lipid_nodes(self) -> None:
        layout = default_layout()
        fusion = pd.DataFrame([
            {"tissue": "liver", "comparison": DIET, "pathway_term": "citrate cycle (tca cycle)",
             "fused_fdr": 0.01, "n_omics": 3, "omics_directions": "proteomics:up"},
            {"tissue": "liver", "comparison": DIET, "pathway_term": "ppar signaling pathway",
             "fused_fdr": 0.007, "n_omics": 2, "omics_directions": "proteomics:up;transcriptomics:up"},
            {"tissue": "liver", "comparison": DIET, "pathway_term": "retinol metabolism",
             "fused_fdr": 0.40, "n_omics": 2, "omics_directions": "proteomics:up"},
        ])
        extra = extra_pathway_fusion(fusion, layout, "liver", DIET)
        self.assertEqual(list(extra["pathway_term"]), ["ppar signaling pathway"])

        tests = pd.DataFrame([
            {"test": "crossover", "name": "HEX_c", "layer": "metabolomics", "tissue": "liver",
             "contrast_id": DIET, "statistic": 0.1, "pvalue": 0.5, "padj": 0.5, "direction": "concordant",
             "n_members": 2, "members": "GLC -> G6P"},
            {"test": "signed_module", "name": "gluconeogenic_routing", "layer": "proteomics", "tissue": "liver",
             "contrast_id": DIET, "statistic": 2.32, "pvalue": 0.0028, "padj": 0.0028, "direction": "as_declared",
             "n_members": 5, "members": "Pcx(+), Pck1(+)"},
        ])
        sm = signed_module_summary(tests, "liver", DIET)
        self.assertEqual(list(sm["name"]), ["gluconeogenic_routing"])

        rows = pd.DataFrame([
            {"layer": "derived", "role": "product", "feature_id": "class_total:TG", "member_id": "TG class total",
             "log2fc": 2.5, "cell_state": "changed"},
            {"layer": "derived", "role": "substrate", "feature_id": "class_total:FFA", "member_id": "FFA class total",
             "log2fc": -0.3, "cell_state": "no_evidence"},
        ])
        self.assertEqual(list(node_members(rows, "TAG_c", {})["feature_id"]), ["class_total:TG"])
        self.assertEqual(list(node_members(rows, "TGPAL_c", {})["feature_id"]), ["class_total:TG"])
        self.assertEqual(list(node_members(rows, "FA_c", {})["feature_id"]), ["class_total:FFA"])


class PathwayMapBundleTest(unittest.TestCase):
    _bundle = _panels.DossierPanelBundleTest
    setUp, tearDown, source_bundle = _bundle.setUp, _bundle.tearDown, _bundle.source_bundle
    set_integration, integrate, run_dossier = _bundle.set_integration, _bundle.integrate, _bundle.run_dossier

    def test_maps_tables_and_values(self) -> None:
        self.source_bundle()
        write_panels(self.root / "experiments" / EXP, PANEL)
        base = self.run_dossier(panels="panels.yaml", tracers=["glucose"], label="A", figures={"formats": ["png"]},
                                pathway_map={"tissues": ["liver"], "contrasts": [DIET], "formats": ["png"]})
        page = base / "figures" / "pathway_map" / "liver" / f"map_{DIET}"
        self.assertTrue(page.with_suffix(".png").is_file())
        exported = pd.read_csv(page.with_suffix(".tsv"), sep="\t")
        table = pd.read_parquet(base / "tables" / "pathway_map.parquet")
        coverage = pd.read_parquet(base / "tables" / "pathway_map_coverage.parquet")
        self.assertEqual(len(exported), len(table))
        self.assertEqual(set(table["element_type"]), {"reaction", "metabolite", "aggregate"})
        self.assertNotIn("unplaced", set(coverage["placement"]))
        self.assertTrue(any("通路总图（图 E）：1 页" in n for n in self.notes))
        # 代谢物节点与图 A 的池格子取同一个值
        members = pd.read_parquet(base / "tables" / "dossier_members.parquet")
        panel = load_panels(self.root / "experiments" / EXP / "panels.yaml")["liver_glucose"]
        matrix = matrix_view(panel, panel_elements(panel), members, "liver", DIET, ["A"], ["glucose"])
        pools = matrix[(matrix["block"] == "pool") & (matrix["element_type"] == "metabolite")]
        nodes = table[table["element_type"] == "metabolite"].set_index("element_id")
        compared = 0
        for row in pools.itertuples():
            if row.element_id in nodes.index:
                a, b = row.value, nodes.at[row.element_id, "value"]
                self.assertTrue((pd.isna(a) and pd.isna(b)) or abs(a - b) < 1e-12, row.element_id)
                self.assertEqual(row.display_state, nodes.at[row.element_id, "display_state"])
                compared += 1
        self.assertGreater(compared, 0)
        # 反应的等级原样取自评估表
        assessed = pd.read_parquet(base / "tables" / "reaction_assessment.parquet")
        claim = assessed[(assessed["flux_id"] == "HEX_c") & (assessed["tissue"] == "liver")
                         & (assessed["contrast_id"] == DIET)]["flux_claim"].iloc[0]
        hex_row = table[(table["element_id"] == "HEX_c")].iloc[0]
        self.assertEqual(hex_row["flux_claim"], claim)
        # 交付目录 pathway_map/<组织>/ 中有图件（只出 PNG 时为 PNG）与视图表，摘要「通路总图」链接到副本
        delivery = Delivery(load_experiment(load_project(self.root), EXP), [], [])
        delivery.bundle("integrate", base.parents[1])
        delivery.write_summary()
        copies = delivery.directory / "pathway_map" / "liver"
        self.assertTrue((copies / f"map_{DIET}.png").is_file())
        self.assertEqual((copies / f"map_{DIET}.tsv").read_bytes(), page.with_suffix(".tsv").read_bytes())
        workbook = load_workbook(delivery.directory / "summary.xlsx")
        try:
            rows = [[c.value for c in row] for row in workbook["通路总图"].iter_rows(min_row=2)]
        finally:
            workbook.close()
        self.assertIn(["liver", DIET, "图件", f"map_{DIET}.png"], rows)
        # omicsweave map 从整合包重画到新目录：视图表与原页相同，原结果包不变
        before = sorted(p.relative_to(base) for p in base.rglob("*"))
        out = self.root / "redraw"
        with contextlib.redirect_stdout(io.StringIO()), self.assertRaises(SystemExit) as exit_code:
            cli.main(["--project-dir", str(self.root), "map", "--experiment", EXP, "--tissue", "liver",
                      "--format", "png", "--output", str(out)])
        self.assertEqual(exit_code.exception.code, 0)
        redrawn = out / "figures" / "pathway_map" / "liver" / f"map_{DIET}.tsv"
        self.assertEqual(redrawn.read_bytes(), page.with_suffix(".tsv").read_bytes())
        self.assertEqual(sorted(p.relative_to(base) for p in base.rglob("*")), before)

    def test_bad_layout_stops_the_run(self) -> None:
        from omicsweave.core.run import RunError

        self.source_bundle()
        (self.root / "experiments" / EXP / "bad_layout.yaml").write_text("canvas: {width: 10}\n", encoding="utf-8")
        with self.assertRaises(RunError):
            self.run_dossier(reactions=["HEX_c"], pathway_map={"layout": "bad_layout.yaml"})

    def test_each_tracer_map_uses_its_own_statistics_and_animals(self) -> None:
        from tests.visualization.test_pathway_map_html import _extract_payload_json

        source = self.source_bundle()
        pooled = self.tables["mfa"]
        base_rows = pooled[(pooled.contrast_id == DIET) & (pooled.tissue == "liver")]
        up = base_rows.assign(tracer="glucose", log2fc=1., mean_numerator=20., mean_denominator=10.,
                              effect_size=1.5, direction_stability=1., n_numerator=4, n_denominator=3)
        down = up.assign(tracer="lactate", log2fc=-1., mean_numerator=5., effect_size=-1.5, n_numerator=2)
        pd.concat([up, down]).to_parquet(source / "omics/mfa/tables/differential_by_tracer.parquet", index=False)
        base = self.run_dossier(reactions=["HEX_c"], tracers=["glucose", "lactate", "absent"], label="A",
                                pathway_map={"tissues": ["liver"], "contrasts": [DIET], "formats": ["html", "svg"],
                                             "submaps": False, "abstract": False})
        summary = base / "figures/pathway_map/liver"
        combined = _extract_payload_json((summary / f"map_{DIET}.html").read_text(encoding="utf-8"))
        self.assertEqual(combined["elements"]["rxn:HEX_c"]["flux_claim"], "LX")
        self.assertFalse((summary / "by_tracer/absent").exists())
        index = (base / "pathway_maps.html").read_text(encoding="utf-8")
        for name, sign, n in (("glucose", 1, 4), ("lactate", -1, 2)):
            stem = summary / "by_tracer" / name / f"map_{DIET}"
            page = _extract_payload_json(stem.with_suffix(".html").read_text(encoding="utf-8"))
            self.assertEqual(page["meta"]["tracer"], name)
            reaction = page["elements"]["rxn:HEX_c"]
            self.assertEqual({m["tracer"] for m in reaction["members"]}, {name})
            row = pd.read_csv(stem.with_suffix(".tsv"), sep="\t").set_index("element_id").loc["HEX_c"]
            self.assertEqual((row["F_direction"], row["F_state"], row["flux_claim"]), (sign, "suggestive", "L2"))
            self.assertEqual((row["n_numerator"], row["hedges_g"]), (n, sign * 1.5))
            self.assertIn(f"by_tracer/{name}/map_{DIET}.html", index)
        delivery = Delivery(load_experiment(load_project(self.root), EXP), [], [])
        delivery.bundle("integrate", base.parents[1])
        copies = delivery.copy_maps()
        self.assertTrue(all((delivery.directory / r["target"]).is_file() for r in copies))
        self.assertIn(f"by_tracer/lactate/map_{DIET}.html", {r["文件"] for r in copies})
        self.assertIn(f"by_tracer/lactate/map_{DIET}.svg", {r["文件"] for r in copies})
        self.assertIn(f"map_{DIET}.svg", {r["文件"] for r in copies})


if __name__ == "__main__":
    unittest.main()
