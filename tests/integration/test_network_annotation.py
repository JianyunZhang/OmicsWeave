"""MFA 网络反应注释与方向证据规则（映射修复阶段 1）。"""

from __future__ import annotations

import unittest
from types import SimpleNamespace

import pandas as pd

from omicsweave.core.idmap.reaction_pathway_mapping import REACTION_REGISTRY
from omicsweave.integration.mappings import (
    FLUX_GENE_MAPPING,
    FLUX_METABOLITE_MAPPING,
    REACTION_ENZYME_MAPPING,
)
from omicsweave.integration.mappings.network_annotation import (
    classify,
    network_reactions,
    network_sources,
    parse_module,
    reaction_semantics,
    role_agreement,
    sync_metabolite_roles,
    tissue_side,
)
from omicsweave.integration.methods.flux_evidence import FluxEvidenceBuilder

COMP = "case_vs_control"


class ParseAndClassifyTest(unittest.TestCase):
    def test_parser_reads_literal_reactions_and_skips_generated_ids(self) -> None:
        source = """
name = 'x'
reaction_list = [
    {'id': 'GPS_c', 'sub': [('G3P_c', 'abc')], 'pro': [('DHAP_c', 'abc')], 'reverse': True},
    {'id': 'FA_usage', 'sub': [('FA_c', 'ab')], 'pro': [('ACCOA_m', 'ab')]},
    {'id': f'{name}_exchange', 'sub': [('X_c', 'a')], 'pro': [('X_e', 'a')]},
]
"""
        reactions = {r["flux_id"]: r for r in parse_module(source)}

        self.assertEqual(set(reactions), {"GPS_c", "FA_usage"})
        self.assertEqual(reactions["GPS_c"]["substrates"], ("G3P_c",))
        self.assertTrue(reactions["GPS_c"]["reversible"])
        self.assertFalse(reactions["FA_usage"]["reversible"])

    def test_reaction_classes(self) -> None:
        cases = {
            ("GPS_c", ("G3P_c",), ("DHAP_c",)): "enzymatic",
            ("PYR_trans", ("PYR_c",), ("PYR_m",)): "transport",
            ("ASPGLU_m", ("ASP_m", "GLU_c"), ("ASP_c", "GLU_m")): "transport",
            ("ALA_exchange", ("ALA_e",), ("ALA_c",)): "exchange",
            ("LDH_e", ("PYR_e",), ("LAC_e",)): "extracellular",
            ("GLC_infuse", ("GLC_inf",), ("GLC_e",)): "pseudo",
            ("Salvage_c", ("RIB5P_stock",), ("RIB5P_c",)): "pseudo",
            ("lipid_mobilization_c", ("LIPID",), ("GLYCEROL_c", "FA_c")): "enzymatic",
        }
        for (flux_id, sub, pro), expected in cases.items():
            with self.subTest(flux_id=flux_id):
                self.assertEqual(classify(flux_id, sub, pro), expected)

    def test_tissue_side_roles(self) -> None:
        self.assertEqual(tissue_side(("ALA_e",), ("ALA_c",)), ([], ["ALA"]))
        self.assertEqual(tissue_side(("PYR_c",), ("PYR_m",)), ([], []))
        self.assertEqual(tissue_side(("ASP_c", "AKG_c"), ("OAC_c", "GLU_c")), (["ASP", "AKG"], ["OAC", "GLU"]))


class CommittedAnnotationTest(unittest.TestCase):
    def test_glycerol_and_alanine_reactions_follow_the_network_sign_convention(self) -> None:
        reactions = network_reactions()

        self.assertEqual((reactions["GPS_c"].substrates, reactions["GPS_c"].products), (("G3P_c",), ("DHAP_c",)))
        self.assertEqual(reactions["ALA_exchange"].reaction_class, "exchange")
        self.assertEqual(reactions["ALA_exchange"].tissue_product_names, ("alanine",))
        self.assertEqual(reactions["GLYCEROL_exchange"].tissue_substrate_names, ("glycerol",))
        self.assertEqual(reactions["LDH_e"].reaction_class, "extracellular")

    def test_sources_record_the_network_commit(self) -> None:
        sources = network_sources()

        self.assertEqual(sources["model_id"], "dio_multi_tissue")
        self.assertRegex(sources["networks"]["carbon"]["sha256_of_modules"], r"^[0-9a-f]{64}$")

    def test_semantics(self) -> None:
        expected = {
            "HEX_c": "directional", "PC_m": "directional", "GPS_c": "capacity", "LDH_c": "capacity",
            "GPT_c": "capacity", "PYR_trans": "directional", "ASPGLU_m": "directional", "CIT_trans": "capacity",
            "ALA_exchange": "capacity", "LDH_e": "not_assessable", "GLC_infuse": "not_assessable",
            "BIOMASS_REACTION_3": "not_assessable", "ALA_protein_input": "not_assessable",
        }
        for flux_id, evidence in expected.items():
            with self.subTest(flux_id=flux_id):
                self.assertEqual(reaction_semantics(flux_id).direction_evidence, evidence)

    def test_every_override_names_a_network_reaction(self) -> None:
        import yaml
        from omicsweave.integration.mappings.network_annotation import SEMANTICS_FILE

        overrides = yaml.safe_load(SEMANTICS_FILE.read_text(encoding="utf-8"))["overrides"]
        self.assertFalse(set(overrides) - set(network_reactions()))


class MappingConsistencyTest(unittest.TestCase):
    def test_metabolite_roles_are_in_sync_with_the_network(self) -> None:
        # 改了网络注释或代谢物映射后，请运行 sync_metabolite_roles 重新生成 YAML
        _, changes = sync_metabolite_roles(FLUX_METABOLITE_MAPPING, registered=REACTION_REGISTRY)
        self.assertEqual(changes, {})
        self.assertEqual(set(role_agreement(FLUX_METABOLITE_MAPPING).values()), {"consistent"})

    def test_glycerol_node_entries(self) -> None:
        self.assertEqual(REACTION_ENZYME_MAPPING["G3PP_c"]["direct_enzymes"], ["Pgp"])
        self.assertEqual(REACTION_ENZYME_MAPPING["GK_c"]["direct_enzymes"], ["Gk"])
        self.assertIn("Dgat1", REACTION_ENZYME_MAPPING["FA_esterification_c"]["direct_enzymes"])
        self.assertEqual(REACTION_ENZYME_MAPPING["lipid_mobilization_c"]["direct_enzymes"], ["Pnpla2", "Lipe", "Mgll"])
        self.assertEqual(FLUX_METABOLITE_MAPPING["G3PP_c"]["substrates"], ["glycerol 3-phosphate"])
        self.assertEqual(FLUX_METABOLITE_MAPPING["G3PP_c"]["products"], ["glycerol"])
        self.assertNotIn("Glycerol-3-phosphate", FLUX_METABOLITE_MAPPING["GLYCEROL_exchange"]["products"])

    def test_transport_carriers_and_compartment_specific_isozymes(self) -> None:
        self.assertEqual(REACTION_ENZYME_MAPPING["PYR_trans"]["transporters"], ["Mpc1", "Mpc2"])
        self.assertEqual(REACTION_ENZYME_MAPPING["ME2_c"]["direct_enzymes"], ["Me1"])
        self.assertEqual(REACTION_ENZYME_MAPPING["PEPCK_c"]["related_isozymes"], ["Pck2"])
        self.assertNotIn("Glul", FLUX_GENE_MAPPING["GLND_m"]["gene_roles"].get("catalytic", []))
        self.assertEqual(REACTION_ENZYME_MAPPING["ASPGLU_m"]["transporters"], ["Slc25a12", "Slc25a13"])


class EvidenceBuilderSemanticsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.builder = FluxEvidenceBuilder(SimpleNamespace(scoring={}))

    def proteins(self, rows):
        return pd.DataFrame({"protein": [r[0] for r in rows], "uniprot_id": [r[1] for r in rows],
                             f"log2fc_{COMP}": [r[2] for r in rows], f"padj_{COMP}": [0.001] * len(rows),
                             f"pval_{COMP}": [0.0001] * len(rows)})

    def test_pyruvate_carrier_gives_direction_for_mitochondrial_import(self) -> None:
        flux_row = pd.Series({"flux_id": "PYR_trans", "flux_log2fc": 1.0})

        details = self.builder.proteomics_evidence(flux_row, self.proteins([("Mpc1", "P63030", 0.9)]), COMP)

        self.assertEqual(details["concordant_proteins"], ["Mpc1"])

    def test_near_equilibrium_enzyme_is_association_only(self) -> None:
        flux_row = pd.Series({"flux_id": "GPS_c", "flux_log2fc": 1.0})
        genes = pd.DataFrame({"gene": ["Gpd1"], "log2fc": [-1.5], "padj": [0.001]})

        tx = self.builder.transcriptomics_evidence(flux_row, genes)
        prot = self.builder.proteomics_evidence(flux_row, self.proteins([("Gpd1", "P13707", -1.5)]), COMP)

        for details in (tx, prot):
            self.assertEqual((details["n_concordant"], details["n_opposite"]), (0, 0))
            self.assertEqual(details["n_unresolved"], 1)

    def test_pyruvate_carboxylase_rna_uses_mgi_symbol_and_carriers_are_not_direct(self) -> None:
        flux_row = pd.Series({"flux_id": "PC_m", "flux_log2fc": 1.0})
        genes = pd.DataFrame({"gene": ["Pcx", "Mpc1"], "log2fc": [1.2, -1.3], "padj": [0.001, 0.001]})

        details = self.builder.transcriptomics_evidence(flux_row, genes)

        self.assertEqual(details["concordant_genes"], ["Pcx"])
        self.assertEqual(details["n_opposite"], 0)          # MPC1 下降只记关联，不判为 PC 的冲突

    def test_pseudo_reactions_are_not_assessable(self) -> None:
        flux_row = pd.Series({"flux_id": "GLC_infuse", "tissue": "serum", "flux_log2fc": 2.0,
                              "flux_fdr": 0.001, "cohen_d": 3.0})

        record = self.builder.build_record(flux_row, None, None, None, COMP)

        self.assertEqual(record["reaction_class"], "pseudo")
        self.assertEqual(record["evidence_grade"], "Not-assessable")


if __name__ == "__main__":
    unittest.main()
