#!/usr/bin/env python3
"""无面板候选表：入选只看全局表，面板只作标注；诱饵、QC、单层通路与方向冲突都能追溯。"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd
import yaml

from omicsweave.contract.manifest import RunManifest, file_record, write_manifest
from omicsweave.integration.candidates import (PATHWAY_TOP_N, REACTION_TOP_N, pathway_candidates, panel_info,
                                               reaction_candidates, write_candidates)
from omicsweave.integration.dossier.panels import PanelError, load_panel_selection, load_panels

EXP = "exp"
RUN = "run_20261119_000000_abcd"
# 均为 MFA 网络中的反应；rank 1..12
REACTIONS = ["HEX_c", "PFK_c", "FBA_c", "TPI_c", "GAPD_c", "PGK_c", "PGM_c", "ENO_c", "PYK_c", "PDH_m", "CS_m",
             "LDH_c"]
PANELS = {"selection": {"mode": "post_hoc", "basis": "依据 run_x 的排序"},
          "panels": {"glycolysis": {"title": "糖酵解", "tissues": ["liver"], "contrasts": ["HFD_vs_LFD"],
                                    "segments": [{"name": "下段", "reactions": ["PYK_c", "LDH_c"]}]}}}


def assessment() -> pd.DataFrame:
    rows = []
    for tissue in ("liver", "heart"):
        for rank, flux in enumerate(REACTIONS, start=1):
            # 排名第 12 的 LDH_c 在肝脏为 L1：不在前 10，也应因结论等级入选
            claim = "L1" if (tissue == "liver" and flux in ("HEX_c", "LDH_c")) else ("L2" if rank == 2 else "L0")
            rows.append({"flux_id": flux, "tissue": tissue, "contrast_id": "HFD_vs_LFD", "reaction_class": "capacity",
                         "flux_claim": claim, "claim_direction": {"L1": 1, "L2": -1}.get(claim, 0),
                         "n_concordant_axes": 1, "discordant_axes": "", "mechanism_class": "", "rho_h_lower": None,
                         "rho_h_upper": None, "conjunction_statement": "", "conjunction_pvalue": None,
                         "conjunction_padj": None, "any_pvalue": 0.1, "completeness": 0.9,
                         "completeness_by_layer": "mfa 1/1", "testable_by": "glucose", "phenomena": "", "rank": rank})
    return pd.DataFrame(rows)


def evidence() -> pd.DataFrame:
    """liver：通路 a 两层同向显著（融合）、b 两层方向相反（融合）、c 两层各自显著但融合未过、d..k 只有 RNA 显著。"""
    rows = []

    def add(term, omics, padj, direction, method="gsea_prerank"):
        rows.append({"tissue": "liver", "comparison": "HFD_vs_LFD", "contrast_role": "disease_effect", "omics": omics,
                     "pathway_term": term, "probability": padj, "rank": 1, "rank_fraction": 0.1, "method": method,
                     "query": "", "gene_set_library": "KEGG_2019_Mouse", "term_name": term.title(),
                     "effect_direction": direction, "source_probability_column": "padj"})

    add("a", "proteomics", 0.01, "up"), add("a", "transcriptomics", 0.02, "up")
    add("b", "proteomics", 0.01, "down"), add("b", "transcriptomics", 0.001, "up")
    add("c", "proteomics", 0.04, "up", "ora"), add("c", "transcriptomics", 0.2, "up")
    for i, term in enumerate("defghijk"):
        add(term, "transcriptomics", 0.01 * (i + 1), "up")
        add(term, "proteomics", 0.9, "down")
    add("z", "transcriptomics", 0.6, "up")        # 哪层都不显著，不进候选
    return pd.DataFrame(rows)


def fusion() -> pd.DataFrame:
    fdr = {"a": 0.01, "b": 0.02, "c": 0.3, "z": 0.9, **{t: 0.5 for t in "defghijk"}}
    return pd.DataFrame([{"tissue": "liver", "comparison": "HFD_vs_LFD", "pathway_term": t, "n_omics": 2,
                          "fused_pvalue": v / 10, "fused_fdr": v} for t, v in fdr.items()])


def build(root: Path, panels: dict | None = PANELS, decoys: bool = True) -> Path:
    run = root / "results" / EXP / RUN
    base = run / "integration" / "reaction_dossier"
    for sub in ("tables", "evaluation", "qc"):
        (base / sub).mkdir(parents=True)
    assessment().to_parquet(base / "tables" / "reaction_assessment.parquet", index=False)
    if decoys:
        pd.DataFrame([{"method": "decoy", "tissue": "", "contrast_id": "", "metric": "L1", "target": 2.0,
                       "decoy_expected": 1.5, "decoy_sd": 0.5, "n_replicates": 3},
                      {"method": "decoy", "tissue": "liver", "contrast_id": "HFD_vs_LFD", "metric": "L1", "target": 2.0,
                       "decoy_expected": 1.5, "decoy_sd": 0.5, "n_replicates": 3}]).to_parquet(
            base / "evaluation" / "assessment_calibration.parquet", index=False)
    pd.DataFrame([{"omics": o, "tissue": "liver", "contrast_id": "HFD_vs_LFD", "n_usable_numerator": 3,
                   "n_usable_denominator": n, "status": s, "reason": r}
                  for o, n, s, r in (("proteomics", 3, "ok", ""), ("transcriptomics", 1, "qc_insufficient", "LFD 只剩 1"))]
                 ).to_parquet(base / "qc" / "contrast_reliability.parquet", index=False)
    tables = run / "integration" / "pathway_convergence" / "tables"
    tables.mkdir(parents=True)
    evidence().to_parquet(tables / "pathway_evidence.parquet", index=False)
    fusion().to_parquet(tables / "pathway_fusion.parquet", index=False)

    parameters = {"source_run": "run_src", "pathway_convergence": {"ora_padj_cutoff": 0.05, "gsea_padj_cutoff": 0.25},
                  "reaction_dossier": {"panels": "panels.yaml" if panels else None}}
    files = []
    if panels:
        path = root / "experiments" / EXP / "panels.yaml"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(yaml.safe_dump(panels, allow_unicode=True, sort_keys=False), encoding="utf-8")
        files.append(file_record(path, root, "input"))
    write_manifest(RunManifest(run_id=RUN, project_id="p", experiment_id=EXP, parameters=parameters, files=files), run)
    return run


class CandidateTest(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.root = Path(self._tmp.name).resolve()

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def test_reactions_selected_by_claim_and_rank_with_annotations(self) -> None:
        run = build(self.root)
        info = panel_info(run, self.root)
        full, short = reaction_candidates(run, info.panels)
        self.assertEqual(len(full), 2 * len(REACTIONS))
        liver = short[short["tissue"] == "liver"].set_index("flux_id")
        # 前 10 名全部入选，第 12 名的 L1 因结论等级入选，第 11 名不入选
        self.assertEqual(set(liver.index), set(REACTIONS[:REACTION_TOP_N]) | {"LDH_c"})
        self.assertEqual(liver.at["LDH_c", "selected_by"], "claim")
        self.assertEqual(liver.at["HEX_c", "selected_by"], "claim, top_n")
        self.assertEqual(liver.at["PFK_c", "direction"], "down")
        # 面板只标注声明的组织与比较：心脏同名反应不标
        self.assertEqual(liver.at["PYK_c", "panels"], "glycolysis")
        self.assertEqual(liver.at["LDH_c", "panels"], "glycolysis")
        heart = short[short["tissue"] == "heart"].set_index("flux_id")
        self.assertEqual(heart.at["PYK_c", "panels"], "")
        # 诱饵与 QC 按组织 × 比较合并
        self.assertEqual((liver.at["HEX_c", "cell_L1"], liver.at["HEX_c", "decoy_L1_expected"]), (2.0, 1.5))
        self.assertTrue(pd.isna(heart.at["HEX_c", "decoy_L1_expected"]))
        self.assertEqual(liver.at["HEX_c", "qc_status"], "transcriptomics: qc_insufficient（LFD 只剩 1）")
        self.assertEqual(liver.at["HEX_c", "group_sizes"], "proteomics 3v3; transcriptomics 3v1")

    def test_panels_do_not_change_the_search(self) -> None:
        run = build(self.root)
        other = {"panels": {"tca": {"segments": [{"reactions": ["PDH_m", "CS_m"]}]}}}
        path = self.root / "other.yaml"
        path.write_text(yaml.safe_dump(other), encoding="utf-8")
        results = [reaction_candidates(run, panels) for panels in (None, load_panels(path),
                                                                    panel_info(run, self.root).panels)]
        stripped = [short.drop(columns="panels") for _, short in results]
        for table in stripped[1:]:
            pd.testing.assert_frame_equal(table, stripped[0])
        self.assertEqual(set(results[1][1]["panels"]) - {""}, {"tca"})

    def test_pathway_support_conflicts_and_single_layer(self) -> None:
        run = build(self.root)
        full, short = pathway_candidates(run)
        table = full.set_index("pathway_term")
        self.assertNotIn("z", table.index)
        self.assertEqual(table.at["a", "support"], "fused")
        self.assertEqual(table.at["b", "direction_agreement"], "conflicting")
        self.assertEqual(table.at["b", "layer_directions"], "proteomics:down; transcriptomics:up")
        # c：ORA 阈值 0.05 下蛋白显著、GSEA 阈值 0.25 下 RNA 显著，融合未过
        self.assertEqual(table.at["c", "support"], "multi_layer")
        self.assertEqual(table.at["d", "support"], "single_layer")
        self.assertEqual(table.at["d", "significant_layers"], "transcriptomics")
        self.assertEqual(table.at["d", "layer_padj"], "proteomics:0.9; transcriptomics:0.01")
        self.assertEqual(set(table["matched_by"]), {"normalized_name"})
        # 简表：融合与多层的全部，单层只留 padj 最小的前 N 条
        single = short[short["support"] == "single_layer"]["pathway_term"]
        self.assertEqual(list(single), list("defghijk")[:PATHWAY_TOP_N])
        self.assertEqual(list(short["support"][:3]), ["fused", "fused", "multi_layer"])

    def test_write_candidates_records_selection_and_decoys(self) -> None:
        run = build(self.root)
        result = write_candidates(run, self.root / "out", self.root)
        self.assertEqual([p.name for p in result.files], ["candidates.xlsx", "reactions_full.tsv", "reactions.tsv",
                                                          "pathways_full.tsv", "pathways.tsv"])
        sheets = pd.read_excel(result.files[0], sheet_name=None)
        self.assertEqual(list(sheets), ["说明", "反应候选", "通路候选"])
        readme = dict(zip(sheets["说明"]["项目"], sheets["说明"]["说明"]))
        self.assertEqual(readme["面板选定方式"], "看过结果后选定")
        self.assertEqual(readme["选定依据"], "依据 run_x 的排序")
        self.assertIn("随机水平", readme["诱饵校准"])
        self.assertEqual(result.panels.record()["selection"], "post_hoc")
        liver = result.overview.set_index("tissue").loc["liver"]
        self.assertEqual((liver["reactions_L1"], liver["pathways_fused"], liver["pathways_fused_conflicting"]),
                         (2, 2, 1))
        tsv = (self.root / "out" / "reactions.tsv").read_bytes()
        self.assertTrue(tsv.startswith(b"\xef\xbb\xbf"))

    def test_without_panels_or_decoys_and_with_changed_panels(self) -> None:
        run = build(self.root, panels=None, decoys=False)
        result = write_candidates(run, self.root / "out", self.root)
        self.assertIsNone(result.panels)
        readme = dict(pd.read_excel(result.files[0], sheet_name="说明").to_numpy())
        self.assertIn("没有声明面板", readme["面板"])
        self.assertIn("未运行诱饵校准", readme["诱饵校准"])

        other = self.root / "second"
        run = build(other)
        (other / "experiments" / EXP / "panels.yaml").write_text(yaml.safe_dump(PANELS) + "# 改过\n", encoding="utf-8")
        info = panel_info(run, other)
        self.assertEqual(info.panels, {})
        self.assertIn("改过", info.note)

    def test_panel_selection_is_validated(self) -> None:
        path = self.root / "panels.yaml"
        base = {"panels": PANELS["panels"]}
        path.write_text(yaml.safe_dump(base), encoding="utf-8")
        self.assertEqual(load_panel_selection(path), {"mode": "unspecified", "basis": ""})
        for doc, message in (({**base, "selection": {"mode": "later"}}, "selection"),
                             ({**base, "selection": {"mode": "post_hoc", "why": "x"}}, "selection"),
                             ({**base, "selction": {"mode": "post_hoc"}}, "无法识别的顶层键")):
            path.write_text(yaml.safe_dump(doc, allow_unicode=True), encoding="utf-8")
            with self.subTest(message), self.assertRaisesRegex(PanelError, message):
                load_panels(path)


if __name__ == "__main__":
    unittest.main()
