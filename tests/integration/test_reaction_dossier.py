#!/usr/bin/env python3
"""反应档案（reaction_dossier，G13 第一部分）：成员行、格子状态，以及与 flux_anchored 证据长表的一致性。

验收（反应档案与分轴综合评估，第 5 节 G13）：
- 整层缺失、成员未检出、测到无变化分别得到三种不同的状态；
- 成员表按层聚合后与 flux_evidence_table 的计数一致；
- flux_anchored 的输出不变（由 test_flux_anchored 与队列 A 的结果包对照保证）。
"""

from __future__ import annotations

import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

from omicsweave.contract.manifest import read_manifest
from omicsweave.integration.mappings.feature_resolver import _canonical_gene_symbol
from omicsweave.integration.methods.flux_evidence import FluxEvidenceBuilder
from omicsweave.integration.modes.reaction_dossier import COLUMNS, _implied, _mfa_state, _state
from tests.integration.test_flux_anchored import DIET, DRUG, FluxAnchoredBundleTest
from tests.mfa.test_mfa_runner import EXP

RULES = {"alpha": 0.05, "equivalence_bound": 0.58}


class CellStateTest(unittest.TestCase):
    def test_measured_states(self) -> None:
        n = {"n_numerator": 3, "n_denominator": 3}
        cases = [
            ({"padj": 0.01, "pvalue": 0.001, "log2fc": 1.0, "se": 0.2}, "changed"),
            ({"padj": 0.2, "pvalue": 0.03, "log2fc": 1.0, "se": 0.4}, "nominal"),
            # 90% 区间：0.05 ± t(0.95, 4)·0.1 = ±0.26，在 ±0.58 内
            ({"padj": 0.9, "pvalue": 0.6, "log2fc": 0.05, "se": 0.1, **n}, "equivalent"),
            # 同样的估计，但 SE 大到区间越出界限
            ({"padj": 0.9, "pvalue": 0.6, "log2fc": 0.05, "se": 0.4, **n}, "no_evidence"),
            # 收缩估计偏向 0：不做等效判定
            ({"padj": 0.9, "pvalue": 0.6, "log2fc": 0.05, "se": 0.1, "lfc_shrinkage": "apeglm", **n}, "no_evidence"),
            # 没有自由度（样本数）时不做等效判定
            ({"padj": 0.9, "pvalue": 0.6, "log2fc": 0.05, "se": 0.1}, "no_evidence"),
            ({"padj": 0.001, "pvalue": 0.0001, "log2fc": 2.0, "qc_status": "qc_confounded"}, "qc_blocked"),
        ]
        for row, expected in cases:
            with self.subTest(row=row):
                self.assertEqual(_state(row, RULES)[0], expected)

    def test_flux_states(self) -> None:
        self.assertEqual(_mfa_state({"n_numerator": 1, "n_denominator": 2, "identifiability": "bound_pinned"},
                                    RULES)[0], "insufficient_n")
        self.assertEqual(_mfa_state({"n_numerator": 4, "n_denominator": 3, "identifiability": "bound_pinned"},
                                    RULES)[0], "not_identifiable")
        self.assertEqual(_mfa_state({"n_numerator": 3, "n_denominator": 3, "identifiability": "moderate",
                                     "log2fc": 1.0, "pvalue": np.nan, "padj": np.nan}, RULES)[0], "no_evidence")

    def test_implied_direction(self) -> None:
        self.assertEqual(_implied("substrate", 1.0, 0), "down")
        self.assertEqual(_implied("regulator", 1.0, -1), "down")             # 抑制性激酶升高 → 通量降
        self.assertEqual(_implied("regulator", 1.0, 0), "")                  # 方向未知或经由转录：不指向
        self.assertEqual(_implied("direct_enzyme", 1.0, 0, "shared_subunit"), "")
        self.assertEqual(_implied("transcription_factor", 1.0, 0), "")


class DossierBundleTest(unittest.TestCase):
    """沿用 flux_anchored 测试的合成结果包：饮食比较有四个组学（心脏只有转录组与 MFA），药物比较每组 1 只。
    只借用其辅助方法，不继承（否则会重复运行它的测试）。"""
    setUp = FluxAnchoredBundleTest.setUp
    tearDown = FluxAnchoredBundleTest.tearDown
    source_bundle = FluxAnchoredBundleTest.source_bundle
    set_integration = FluxAnchoredBundleTest.set_integration
    integrate = FluxAnchoredBundleTest.integrate

    def dossier(self, **parameters) -> tuple[pd.DataFrame, Path, list[str]]:
        self.source_bundle()
        if parameters:
            self.set_integration(enabled=False, modes=["flux_anchored", "reaction_dossier"],
                                 parameters={"reaction_dossier": parameters})
        run_dir, notes = self.integrate(["flux_anchored", "reaction_dossier"])
        table = pd.read_parquet(run_dir / "integration" / "reaction_dossier" / "tables" / "dossier_members.parquet")
        return table, run_dir, notes

    @staticmethod
    def cell(table: pd.DataFrame, flux_id: str, tissue: str, contrast: str = DIET) -> pd.DataFrame:
        return table[(table["flux_id"] == flux_id) & (table["tissue"] == tissue) & (table["contrast_id"] == contrast)]

    def test_three_distinct_states(self) -> None:
        table, _, _ = self.dossier()
        self.assertEqual(list(table.columns), COLUMNS)
        heart = self.cell(table, "PDH_m", "heart")
        # 整层缺失：心脏没有蛋白组与代谢组
        self.assertEqual(set(heart.loc[heart["layer"] == "proteomics", "cell_state"]), {"not_measured"})
        self.assertEqual(set(heart.loc[heart["layer"] == "metabolomics", "cell_state"]) - {"no_knowledge"},
                         {"not_measured"})
        tx = heart[heart["layer"] == "transcriptomics"].set_index("member_id")
        self.assertEqual(tx.at["Pdha1", "cell_state"], "changed")          # padj 0.01
        self.assertEqual(tx.at["Dlat", "cell_state"], "no_evidence")       # 测到，无变化证据
        self.assertEqual(tx.at["Dld", "cell_state"], "not_detected")       # 映射了，数据中没有
        self.assertEqual(tx.at["Pdk4", "regulator_effect"], "inhibits")
        self.assertFalse(tx.at["Dld", "direction_member"])                 # 共用亚基
        self.assertTrue(tx.at["Pdha1", "direction_member"])
        self.assertEqual(tx.at["Pdha1", "implied_direction"], "up")

    def test_knowledge_and_effectors(self) -> None:
        table, _, _ = self.dossier()
        pyk = self.cell(table, "PYK_c", "liver")
        effectors = pyk[pyk["role"] == "effector"].set_index("member_id")
        self.assertEqual(effectors.at["Alanine", "cell_state"], "no_evidence")     # 测到（padj 0.8）
        self.assertEqual(effectors.at["Alanine", "regulator_effect"], "inhibits")
        self.assertEqual(effectors.at["Fructose-1,6-bisphosphate", "cell_state"], "not_detected")
        self.assertFalse(effectors["in_flux_anchored"].any())
        pc = self.cell(table, "PC_m", "liver")
        no_knowledge = pc[pc["cell_state"] == "no_knowledge"]
        self.assertEqual(set(zip(no_knowledge["layer"], no_knowledge["role"])),
                         {("transcriptomics", "regulator"), ("proteomics", "regulator")})
        self.assertTrue(no_knowledge["member_id"].eq("").all())

    def test_tracers(self) -> None:
        table, run_dir, notes = self.dossier(tracers=["glucose", "alanine"])
        flux = self.cell(table, "HEX_c", "liver")
        flux = flux[flux["layer"] == "mfa"].set_index("tracer")
        self.assertEqual(flux.at["pooled", "cell_state"], "changed")           # Mann-Whitney FDR 0.01
        self.assertEqual((flux.at["glucose", "cell_state"], flux.at["alanine", "cell_state"]),
                         ("not_infused", "not_infused"))
        self.assertTrue(any("not_infused" in n for n in notes))
        # 没有通量数据的组织 × 比较：只有一行合并的 not_measured
        absent = self.cell(table, "PC_m", "liver", DRUG)
        self.assertIn("not_measured", set(absent.loc[absent["tracer"] == "pooled", "cell_state"]))

    def test_single_named_tracer(self) -> None:
        """MFA 样本只有一种示踪剂（protocol 只有一个取值）时没有分示踪剂表：统一差异表按该示踪剂另列一行，
        声明的这种示踪剂不再记为 not_infused。"""
        path = self.root / "experiments" / EXP / "samples.tsv"
        samples = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        if "protocol" not in samples.columns:
            samples.insert(4, "protocol", "")
        samples.loc[samples["omics"] == "mfa", "protocol"] = "glucose"
        samples.to_csv(path, sep="\t", index=False, lineterminator="\n")
        table, _, notes = self.dossier(tracers=["glucose", "alanine"])
        flux = self.cell(table, "HEX_c", "liver")
        flux = flux[flux["layer"] == "mfa"].set_index("tracer")
        self.assertEqual(flux.at["glucose", "cell_state"], flux.at["pooled", "cell_state"])
        self.assertEqual(flux.at["glucose", "log2fc"], flux.at["pooled", "log2fc"])
        self.assertEqual(flux.at["alanine", "cell_state"], "not_infused")
        self.assertTrue(any("只有一种示踪剂（glucose）" in n for n in notes), notes)
        self.assertFalse(any("glucose 在结果包中没有数据" in n for n in notes), notes)

    def test_counts_match_the_evidence_table(self) -> None:
        """按层汇总后（只计 flux_anchored 也计入的行），映射数、匹配数、入选数与证据长表一致。"""
        table, run_dir, _ = self.dossier()
        evidence = pd.read_parquet(run_dir / "integration" / "flux_anchored" / "tables" / "flux_evidence_table.parquet")
        counted = table[table["in_flux_anchored"].astype(bool)]
        builder = FluxEvidenceBuilder(SimpleNamespace(experiment_id=EXP, estimand={}, scoring={}))
        checked = 0
        for row in evidence[evidence["layer"].isin(["transcriptomics", "proteomics", "metabolomics"])].itertuples():
            layers = ["metabolomics", "derived"] if row.layer == "metabolomics" else [row.layer]
            rows = counted[(counted["flux_id"] == row.flux_id) & (counted["tissue"] == row.tissue)
                           & (counted["contrast_id"] == row.comparison) & counted["layer"].isin(layers)]
            if rows.empty or rows["cell_state"].eq("not_measured").all():
                continue
            found = rows[rows["feature_id"].notna() & rows["feature_id"].ne("")]
            with self.subTest(flux_id=row.flux_id, tissue=row.tissue, contrast=row.comparison, layer=row.layer):
                self.assertEqual(int(rows["selected"].astype(bool).sum()), row.n_signal)
                if row.layer == "metabolomics":
                    self.assertEqual(len(found), row.n_matched)
                else:
                    members = (builder.gene_members if row.layer == "transcriptomics"
                               else builder.protein_members)(row.flux_id)[0]
                    # flux_anchored 按名称计映射数（同一基因的别名如 Gpi 与 Gpi1 计两次），档案按规范符号计一次
                    self.assertEqual(len(members), row.n_mapped)
                    self.assertEqual(rows["member_id"].nunique(),
                                     len({_canonical_gene_symbol(m) or m for m in members}))
                    self.assertEqual(found["member_id"].nunique(), row.n_matched)
                checked += 1
        self.assertGreater(checked, 5)

    def test_parameters_inherit_from_flux_anchored(self) -> None:
        self.source_bundle()
        self.set_integration(enabled=False, modes=["reaction_dossier"], parameters={
            "flux_anchored": {"tissue_map": {"heart": "cardiac"}}, "reaction_dossier": {"reactions": ["PDH_m"]}})
        run_dir, _ = self.integrate()
        table = pd.read_parquet(run_dir / "integration" / "reaction_dossier" / "tables" / "dossier_members.parquet")
        self.assertEqual(set(table["flux_id"]), {"PDH_m"})
        self.assertIn("cardiac", set(table.loc[table["layer"] == "mfa", "tissue"]))
        parameters = read_manifest(run_dir).parameters["reaction_dossier"]
        self.assertEqual(parameters["inherited_from_flux_anchored"], ["tissue_map"])
        with self.assertRaisesRegex(Exception, "无法识别的参数"):
            self.set_integration(enabled=False, modes=["reaction_dossier"],
                                 parameters={"reaction_dossier": {"colours": "red"}})
            self.integrate()


if __name__ == "__main__":
    unittest.main()
