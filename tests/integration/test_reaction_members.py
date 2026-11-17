#!/usr/bin/env python3
"""反应成员知识表（G17）。

验收（反应档案与分轴综合评估，第 5 节 G17）：
- 审计在旧映射上报出 PDH_m 的角色不一致（基因表把 Pdk、Pdp 记为催化成员，蛋白表记为调控因子）；
- 基因表与蛋白表中登记反应的角色由同一张表生成；
- 作用方向为 unknown、或经由转录 / 信号起作用的调控因子不给出活性方向；
- 同工酶加权只在 isozyme 之间进行；共用亚基不作为方向成员。
"""

from __future__ import annotations

import re
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import yaml

from omicsweave.core.idmap.reaction_pathway_mapping import REACTION_REGISTRY
from omicsweave.core.refs import REGISTRIES, build, registry_status
from omicsweave.integration.mappings import FLUX_GENE_MAPPING, REACTION_ENZYME_MAPPING
from omicsweave.integration.mappings.flux_protein_mapping import PROTEIN_CATEGORIES
from omicsweave.integration.mappings.model_package import CuratedNetworkPackage, DecoyPackage
from omicsweave.integration.mappings.reaction_members import (
    REACTION_MEMBERS, audit_reaction_roles, derive_role_tables, load_reaction_members, member_relations,
)
from omicsweave.integration.methods.flux_evidence import FluxEvidenceBuilder

C = "HFD_vs_LFD"
PANEL = ("PYR_trans", "PC_m", "PDH_m", "PEPCK_c", "G6Pase_c", "PYK_c")

# PDH_m 在 2026-11-17 统一角色之前的两张表（基因表把激酶与磷酸酶记在 catalytic 下）
OLD_PDH_GENES = {"PDH_m": {
    "genes": ["Pdha1", "Pdhb", "Dlat", "Dld", "Pdk1", "Pdk2", "Pdk3", "Pdk4", "Pdp1", "Pdp2", "Ppara", "Hif1a"],
    "gene_roles": {"catalytic": ["Pdha1", "Pdhb", "Dlat", "Dld", "Pdk1", "Pdk2", "Pdk3", "Pdk4", "Pdp1", "Pdp2"],
                   "transcription_factor": ["Ppara", "Hif1a"]},
    "pathway": "Pyruvate Metabolism", "description": "Pyruvate dehydrogenase"}}
OLD_PDH_PROTEINS = {"PDH_m": {
    "direct_enzymes": ["Pdha1", "Pdhb", "Dlat", "Dld"],
    "upstream_regulators": ["Pdk1", "Pdk2", "Pdk3", "Pdk4", "Pdp1", "Pdp2"],
    "downstream_enzymes": ["Cs"], "transporters": ["Mpc1", "Mpc2"], "cofactors": [],
    "description": "Pyruvate Dehydrogenase Complex"}}


def _members(entry: dict, key: str) -> set[str]:
    if key == "genes":
        return set(entry.get("genes", []))
    return {m for c in PROTEIN_CATEGORIES for m in entry.get(c, []) or []}


def _parse(text: str):
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
        path = Path(tmp) / "members.yaml"
        path.write_text(text, encoding="utf-8")
        return load_reaction_members(path)


class AnnotationTableTest(unittest.TestCase):
    def test_panel_reactions_are_annotated_and_registered(self) -> None:
        self.assertTrue(set(PANEL) <= set(REACTION_MEMBERS))
        self.assertEqual(sorted(set(REACTION_MEMBERS) - set(REACTION_REGISTRY)), [])

    def test_known_signs_have_verified_sources(self) -> None:
        """方向为激活或抑制的调控因子与效应物都要有文献（PMID）。"""
        for flux_id, annotated in REACTION_MEMBERS.items():
            for member in annotated.members.values():
                with self.subTest(reaction=flux_id, member=member.id):
                    self.assertTrue(all(re.fullmatch(r"PMID:\d+", s) for s in member.sources))
                    if member.effect in ("activates", "inhibits"):
                        self.assertTrue(member.sources)

    def test_activity_sign(self) -> None:
        pdh, pc, g6p, pyk = (REACTION_MEMBERS[r].members for r in ("PDH_m", "PC_m", "G6Pase_c", "PYK_c"))
        self.assertEqual((pdh["Pdk4"].activity_sign, pdh["Pdp1"].activity_sign), (-1, 1))
        self.assertEqual((pdh["Acetyl-CoA"].activity_sign, pc["Acetyl-CoA"].activity_sign), (-1, 1))
        self.assertEqual(pyk["Fructose-1,6-bisphosphate"].activity_sign, 1)
        # unknown，或经由信号 / 转录起作用：只展示，不给方向
        self.assertEqual((pc["Aspartate"].activity_sign, g6p["Insr"].activity_sign), (0, 0))
        self.assertEqual(pdh["Ppara"].activity_sign, 0)

    def test_relations(self) -> None:
        self.assertEqual(member_relations("PDH_m"), {"Pdha1": "complex_subunit", "Pdhb": "complex_subunit",
                                                      "Dlat": "complex_subunit", "Dld": "shared_subunit",
                                                      "Mpc1": "complex_subunit", "Mpc2": "complex_subunit"})
        self.assertEqual(member_relations("PYK_c"), {"Pkm": "isozyme", "Pklr": "isozyme"})
        self.assertEqual(member_relations("HEX_c"), {})                 # 未登记：调用方按同工酶处理
        self.assertEqual(DecoyPackage(CuratedNetworkPackage()).member_relations("PDH_m"), {})

    def test_schema_errors(self) -> None:
        cases = {
            "R:\n  regulators: [Pdk1]\n": "effect",                                   # 调控因子必须有方向
            "R:\n  regulators:\n    Pdk1: {effect: blocks, mode: phosphorylation}\n": "effect",
            "R:\n  regulators:\n    Pdk1: {effect: inhibits, mode: phosphorylation, relation: isozyme}\n": "relation",
            "R:\n  direct_enzymes: [Pkm]\n  downstream_enzymes: [Pkm]\n": "同时属于",
            "R:\n  enzymes: [Pkm]\n": "未知角色",
            "R:\n  direct_enzymes:\n    Pkm: {effect: activates}\n": "effect",
        }
        for text, message in cases.items():
            with self.subTest(text=text), self.assertRaisesRegex(ValueError, message):
                _parse(text)


class RoleTableTest(unittest.TestCase):
    def test_tables_are_generated_from_the_annotation(self) -> None:
        genes, proteins, changed = derive_role_tables(FLUX_GENE_MAPPING, REACTION_ENZYME_MAPPING)
        self.assertEqual(changed, [], "运行 omicsweave refs build reaction_roles 同步两张映射表")
        self.assertEqual(build("reaction_roles", ["--check"]), 0)

    def test_old_pdh_mapping_is_reported_and_regrouped(self) -> None:
        rows = audit_reaction_roles(OLD_PDH_GENES, OLD_PDH_PROTEINS, annotations={})
        mismatched = sorted(r[4] for r in rows if r[7] == "role_mismatch")
        self.assertEqual(mismatched, ["Pdk1", "Pdk2", "Pdk3", "Pdk4", "Pdp1", "Pdp2"])
        self.assertEqual(sum(r[7] == "effect_missing" for r in rows), 6)

        genes, proteins, changed = derive_role_tables(OLD_PDH_GENES, OLD_PDH_PROTEINS)
        self.assertEqual(changed, ["PDH_m"])
        roles = genes["PDH_m"]["gene_roles"]
        self.assertEqual(roles["catalytic"], ["Pdha1", "Pdhb", "Dlat", "Dld"])
        self.assertEqual(roles["regulatory"], ["Pdk1", "Pdk2", "Pdk3", "Pdk4", "Pdp1", "Pdp2"])
        # 只改分组，不增删成员：flux_anchored 匹配到的集合不变
        self.assertEqual(_members(genes["PDH_m"], "genes"), _members(OLD_PDH_GENES["PDH_m"], "genes"))
        self.assertEqual(_members(proteins["PDH_m"], "proteins"), _members(OLD_PDH_PROTEINS["PDH_m"], "proteins"))
        self.assertEqual(audit_reaction_roles(genes, proteins, {"PDH_m": REACTION_MEMBERS["PDH_m"]}), [])

    def test_annotated_reactions_pass_the_audit(self) -> None:
        rows = [r for r in audit_reaction_roles() if r[3] in REACTION_MEMBERS]
        self.assertEqual(rows, [])

    def test_registry_listing(self) -> None:
        status = registry_status().set_index("name")
        self.assertIn("reaction_roles", REGISTRIES)
        self.assertTrue(status.at["reaction_roles", "present"])
        self.assertEqual(status.at["reaction_roles", "n_rows"], len(REACTION_MEMBERS))
        self.assertEqual(status.at["reaction_roles", "generated_on"], "")

    def test_generated_yaml_round_trips(self) -> None:
        """改写后的映射表按 PyYAML 的默认格式写出，再读回内容不变（生成器不会引入格式噪声）。"""
        from omicsweave.integration.mappings.reaction_members import GENE_TABLE, PROTEIN_TABLE
        for path in (GENE_TABLE, PROTEIN_TABLE):
            text = path.read_text(encoding="utf-8").replace("\r\n", "\n")
            dumped = yaml.safe_dump(yaml.safe_load(text), sort_keys=False, allow_unicode=True, default_flow_style=False)
            self.assertEqual(dumped, text, path.name)


def _builder() -> FluxEvidenceBuilder:
    return FluxEvidenceBuilder(SimpleNamespace(experiment_id="toy", estimand={}, scoring={
        "evidence_profile": {"mfa_inference_unit": "biological_replicate"}}))


def _proteins(rows: dict[str, tuple[float, float, float]]) -> pd.DataFrame:
    """{蛋白: (丰度, log2FC, padj)}"""
    return pd.DataFrame({"protein": list(rows), "abundance": [v[0] for v in rows.values()],
                         f"log2fc_{C}": [v[1] for v in rows.values()], f"padj_{C}": [v[2] for v in rows.values()],
                         f"pval_{C}": [v[2] / 2 for v in rows.values()], f"evidence_{C}": "inferential",
                         f"n_eff_{C}": 6})


class ComplexSubunitTest(unittest.TestCase):
    flux = pd.Series({"flux_id": "PDH_m", "tissue": "liver", "flux_log2fc": 0.8})

    def test_shared_subunit_is_not_a_direction_member(self) -> None:
        """只有共用亚基 DLD 与通量同向：不算支持（它也反映其他复合体），也不按“次要同工酶”处理。"""
        details = _builder().proteomics_evidence(self.flux, _proteins({
            "Pdha1": (1000.0, -0.1, 0.8), "Pdhb": (800.0, -0.05, 0.9), "Dlat": (900.0, 0.0, 0.95),
            "Dld": (50.0, 1.2, 0.01)}), C)
        self.assertNotIn("Dld", details["concordant_proteins"])
        self.assertEqual(details["n_concordant"], 0)
        self.assertFalse(details.get("minor_isozyme_driven", False))
        self.assertNotIn("weighted_capacity_change", details)          # 复合体亚基之间不做同工酶加权

    def test_core_subunit_support_is_kept(self) -> None:
        """丰度占比低的核心亚基与通量同向：复合体不按同工酶加权，所以不被当作次要同工酶撤销。"""
        details = _builder().proteomics_evidence(self.flux, _proteins({
            "Pdha1": (1000.0, -0.3, 0.6), "Pdhb": (60.0, 1.0, 0.01), "Dlat": (900.0, -0.2, 0.7)}), C)
        self.assertEqual(details["concordant_proteins"], ["Pdhb"])
        self.assertFalse(details.get("minor_isozyme_driven", False))

    def test_transcriptomics_follows_the_same_relations(self) -> None:
        tx = pd.DataFrame({"gene": ["Pdha1", "Pdhb", "Dlat", "Dld"], "base_mean": [500.0, 400.0, 450.0, 30.0],
                           "log2fc": [-0.1, 0.0, -0.05, 1.5], "padj": [0.8, 0.9, 0.85, 0.001]})
        details = _builder().transcriptomics_evidence(self.flux, tx)
        self.assertNotIn("Dld", details["concordant_genes"])
        self.assertNotIn("weighted_capacity_change", details)


if __name__ == "__main__":
    unittest.main()
