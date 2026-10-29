#!/usr/bin/env python3
"""多组学特征精准映射的公共接口回归测试。"""

from __future__ import annotations

import unittest
from types import SimpleNamespace
import tempfile
from pathlib import Path

import pandas as pd

from omicsweave.core.idmap.metabolite_id_mapping import METABOLITE_DATABASE_IDS
from omicsweave.integration.mappings import (
    FLUX_GENE_MAPPING,
    REACTION_ENZYME_MAPPING,
    audit_mapping_integrity,
    match_gene_features,
    match_metabolite_roles,
    match_protein_features,
)
from omicsweave.integration.methods.flux_omics_scorer import FluxOmicsScorer
from omicsweave.integration.loaders.transcriptomics_loader import TranscriptomicsLoader
from omicsweave.integration.loaders.mfa_loader import MFALoader


class GeneFeatureMappingTest(unittest.TestCase):
    def test_gene_alias_case_and_standard_ids_resolve_to_canonical_mouse_gene(self) -> None:
        genes = pd.DataFrame({
            "gene": ["GPI", "unresolved"],
            "ensembl_gene_id": ["", "ENSMUSG00000022143"],
        })

        matched = match_gene_features(genes, ["Gpi1"])

        self.assertEqual(matched.indices, (0, 1))
        self.assertEqual(matched.matched_queries, ("Gpi",))
        self.assertEqual(
            matched.match_method_by_index,
            {0: "canonical_symbol", 1: "ensembl_gene_id"},
        )

    def test_regulatory_mrna_is_association_only_not_direct_flux_support(self) -> None:
        scorer = FluxOmicsScorer(SimpleNamespace(scoring={}))
        flux_row = pd.Series({"flux_id": "HEX_c", "flux_log2fc": 1.0})
        regulatory_only = pd.DataFrame({
            "gene": ["Gckr"],
            "log2fc": [1.5],
            "padj": [0.01],
        })

        score, details = scorer.score_transcriptomics(flux_row, regulatory_only)

        self.assertEqual(score, 0.0)
        self.assertEqual(details["direction_source"], "indirect_gene")
        self.assertEqual(details["n_concordant"], 0)
        self.assertEqual(details["n_opposite"], 0)
        self.assertEqual(details["n_unresolved"], 1)

    def test_tkt2_direct_gene_mapping_excludes_non_ppp_symbols(self) -> None:
        mapping = FLUX_GENE_MAPPING["TKT2_c"]

        self.assertEqual(set(mapping["gene_roles"]["catalytic"]), {"Tkt"})
        self.assertNotIn("Pkp2", mapping["genes"])
        self.assertNotIn("Purat", mapping["genes"])


class ProteinFeatureMappingTest(unittest.TestCase):
    def test_uniprot_id_is_primary_and_gene_alias_is_safe_fallback(self) -> None:
        proteins = pd.DataFrame({
            "protein": ["unresolved", "gpi1"],
            "uniprot_id": ["P70695", ""],
        })

        matched = match_protein_features(proteins, ["Fbp2", "Gpi"])

        self.assertEqual(matched.indices, (0, 1))
        self.assertEqual(matched.matched_queries, ("Fbp2", "Gpi"))
        self.assertEqual(
            matched.match_method_by_index,
            {0: "uniprot_id", 1: "canonical_symbol"},
        )

    def test_regulatory_protein_is_association_only_not_direct_flux_support(self) -> None:
        scorer = FluxOmicsScorer(SimpleNamespace(scoring={}))
        flux_row = pd.Series({"flux_id": "HEX_c", "flux_log2fc": 1.0})
        regulatory_only = pd.DataFrame({
            "protein": ["Gckr"],
            "uniprot_id": ["Q91X44"],
            "pval_case_vs_control": [0.01],
            "log2fc_case_vs_control": [1.5],
        })

        score, details = scorer.score_proteomics(
            flux_row, regulatory_only, "case_vs_control"
        )

        self.assertEqual(score, 0.0)
        self.assertEqual(details["direction_source"], "indirect")
        self.assertEqual(details["n_concordant"], 0)
        self.assertEqual(details["n_opposite"], 0)
        self.assertEqual(details["n_unresolved"], 1)

    def test_symbol_fallback_is_allowed_only_when_query_has_no_uniprot_registry_id(self) -> None:
        proteins = pd.DataFrame({
            "protein": ["Uqcrfs1", "Fbp2"],
            "uniprot_id": ["vendor-accession", "P21550"],
        })

        matched = match_protein_features(proteins, ["Uqcrfs1", "Fbp2"])

        self.assertEqual(matched.indices, (0,))
        self.assertEqual(matched.matched_queries, ("Uqcrfs1",))
        self.assertEqual(
            matched.match_method_by_index,
            {0: "canonical_symbol_no_registry_id"},
        )

    def test_glycerol_phosphate_shuttle_does_not_map_to_glycogen_enzymes(self) -> None:
        protein_mapping = REACTION_ENZYME_MAPPING["GPS_c"]
        gene_mapping = FLUX_GENE_MAPPING["GPS_c"]

        self.assertEqual(set(protein_mapping["direct_enzymes"]), {"Gpd1", "Gpd2"})
        self.assertTrue(
            {"Gys1", "Gys2", "Pygb", "Pygm", "Pygl"}.isdisjoint(
                protein_mapping["direct_enzymes"]
            )
        )
        self.assertEqual(
            set(gene_mapping["gene_roles"]["catalytic"]),
            {"Gpd1", "Gpd2"},
        )

    def test_missing_proteomics_comparison_columns_are_no_data_not_signal(self) -> None:
        scorer = FluxOmicsScorer(SimpleNamespace(scoring={}))
        flux_row = pd.Series({"flux_id": "HEX_c", "flux_log2fc": 1.0})
        proteins = pd.DataFrame({
            "protein": ["Gck"],
            "uniprot_id": ["P52792"],
        })

        score, details = scorer.score_proteomics(
            flux_row, proteins, "missing_vs_comparison"
        )

        self.assertEqual(score, 0.0)
        self.assertEqual(details["evidence_mode"], "no_data")
        self.assertEqual(details["n_sig"], 0)
        self.assertEqual(details["data_reason"], "missing_comparison_effect")


class StandardizedFeatureLoaderTest(unittest.TestCase):
    def test_transcriptomics_loader_reads_tissue_prefixed_multiomics_export(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            tissue_dir = root / "04_de_results" / "liver"
            tissue_dir.mkdir(parents=True)
            pd.DataFrame({
                "gene_symbol": ["Gpi"],
                "ensembl_gene_id": ["ENSMUSG00000022143"],
                "ncbi_gene_id": ["14751"],
                "log2fc": [1.2],
                "pvalue": [0.001],
                "padj": [0.01],
            }).to_csv(
                tissue_dir / "de_liver_for_multiomics_case_vs_control.csv",
                index=False,
            )
            config = SimpleNamespace(
                transcriptomics_results=root,
                tissues=["liver"],
                comparison_map={"case_vs_control": "case_vs_control"},
                transcriptomics_comparison_map={"case_vs_control": "case_vs_control"},
                mfa_tissue_map={},
                integration={},
            )

            loaded = TranscriptomicsLoader(config).load()
            row = loaded["de"]["case_vs_control"]["liver"].iloc[0]

            self.assertEqual(row["gene"], "Gpi")
            self.assertEqual(row["ensembl_gene_id"], "ENSMUSG00000022143")
            self.assertEqual(str(row["ncbi_gene_id"]), "14751")

    def test_mfa_loader_uses_configured_comparisons_as_allowlist(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            diff_dir = root / "01_df_results" / "liver" / "differential_flux"
            diff_dir.mkdir(parents=True)
            for comparison in ("case_vs_control", "historical_vs_control"):
                pd.DataFrame({
                    "flux_id": ["liver--HEX_c"],
                    "tissue": ["liver"],
                    "log2FC": [1.0],
                }).to_csv(
                    diff_dir / f"differential_flux_{comparison}.csv",
                    index=False,
                )
            config = SimpleNamespace(
                mfa_results=root,
                mfa_tissue_map={},
                comparisons=["case_vs_control"],
                integration={},
            )

            loaded = MFALoader(config).load()

            self.assertEqual(set(loaded["diff"]), {"case_vs_control"})


class DatabaseRegistryTest(unittest.TestCase):
    def test_confirmed_hmdb_accessions_are_current_primary_ids(self) -> None:
        self.assertEqual(
            METABOLITE_DATABASE_IDS["Palmitate"]["hmdb_id"],
            "HMDB0000220",
        )
        self.assertEqual(
            METABOLITE_DATABASE_IDS["Inositol 1-phosphate"]["hmdb_id"],
            "HMDB0000213",
        )
        self.assertEqual(
            METABOLITE_DATABASE_IDS["Ribonic acid"]["hmdb_id"],
            "HMDB0000867",
        )

    def test_mapping_integrity_audit_has_no_blocking_errors(self) -> None:
        audit = audit_mapping_integrity()

        self.assertEqual(audit.errors, ())
        self.assertGreater(audit.metrics["metabolite_aliases_without_standard_id"], 0)
        self.assertEqual(audit.metrics["gene_unknown_symbols"], 0)


class EvidenceMappingAuditOutputTest(unittest.TestCase):
    def test_evidence_table_exposes_mapping_provenance_and_ambiguity(self) -> None:
        scorer = FluxOmicsScorer(SimpleNamespace(scoring={}))
        mfa = {
            "case_vs_control": pd.DataFrame({
                "flux_id": ["ACONT_m"],
                "tissue": ["liver"],
                "flux_log2fc": [1.0],
                "flux_fdr": [0.01],
                "cohen_d": [1.2],
            }),
        }
        metabolites = {
            "liver": pd.DataFrame({
                "compound": ["vendor citrate feature"],
                "hmdb_id": ["HMDB0000094"],
                "pval_case_vs_control": [0.001],
                "fdr_case_vs_control": [0.01],
                "log2fc_case_vs_control": [-1.5],
            }),
        }

        scorer.score_all(mfa, {}, {}, metabolites)
        evidence = scorer.get_evidence_table()
        row = evidence.loc[evidence["layer"].eq("metabolomics")].iloc[0]

        self.assertEqual(row["mapping_methods"], "database_id:1")
        self.assertEqual(row["n_ambiguous"], 0)
        self.assertEqual(row["data_reason"], "")


class MetaboliteFeatureMappingTest(unittest.TestCase):
    def test_metabolite_names_require_exact_normalized_aliases(self) -> None:
        compounds = pd.DataFrame({
            "compound": [
                "Citrate",
                "Isocitrate",
                "Piperidine derivative",
                "PI(18:0_20:4)",
                "Phosphoserine",
                "Ser",
            ],
        })
        roles = {
            "substrates": ["Citrate", "Pi", "Ser"],
            "products": ["Isocitrate"],
            "related": [],
        }

        matched = match_metabolite_roles(compounds, roles)

        self.assertEqual(matched.indices_by_role["substrates"], (0, 5))
        self.assertEqual(matched.indices_by_role["products"], (1,))
        self.assertEqual(matched.indices_by_role["related"], ())
        self.assertEqual(matched.ambiguous_indices, ())

    def test_standard_database_id_takes_priority_over_conflicting_name(self) -> None:
        compounds = pd.DataFrame({
            "compound": ["Isocitrate", "unresolved vendor feature"],
            "hmdb_id": ["HMDB0000094", "HMDB0000094"],
        })
        roles = {
            "substrates": ["Citrate"],
            "products": ["Isocitrate"],
            "related": [],
        }

        matched = match_metabolite_roles(compounds, roles)

        self.assertEqual(matched.indices_by_role["substrates"], (0, 1))
        self.assertEqual(matched.indices_by_role["products"], ())
        self.assertEqual(matched.match_method_by_index, {0: "database_id", 1: "database_id"})

    def test_metabolite_synonyms_count_as_one_mapped_biological_target(self) -> None:
        compounds = pd.DataFrame({"compound": ["Citrate"]})
        roles = {
            "substrates": ["Citrate", "Citric acid"],
            "products": [],
            "related": [],
        }

        matched = match_metabolite_roles(compounds, roles)

        self.assertEqual(matched.target_count_by_role["substrates"], 1)
        self.assertEqual(sum(matched.target_count_by_role.values()), 1)

    def test_metabolomics_scorer_rejects_short_alias_substring_false_positives(self) -> None:
        scorer = FluxOmicsScorer(SimpleNamespace(scoring={}))
        flux_row = pd.Series({"flux_id": "GAPD_c", "flux_log2fc": 1.0})
        metabolites = pd.DataFrame({
            "compound": ["Piperidine derivative", "PI(18:0_20:4)"],
            "pval_case_vs_control": [0.001, 0.001],
            "fdr_case_vs_control": [0.01, 0.01],
            "log2fc_case_vs_control": [1.5, 1.5],
        })

        score, details = scorer.score_metabolomics(
            flux_row, metabolites, "case_vs_control"
        )

        self.assertEqual(score, 0.0)
        self.assertEqual(details["evidence_mode"], "no_match")
        self.assertEqual(details["n_metabolites_found"], 0)

    def test_missing_metabolomics_comparison_columns_are_no_data_not_signal(self) -> None:
        scorer = FluxOmicsScorer(SimpleNamespace(scoring={}))
        flux_row = pd.Series({"flux_id": "ACONT_m", "flux_log2fc": 1.0})
        metabolites = pd.DataFrame({
            "compound": ["Citrate", "Isocitrate"],
            "HMDB": ["HMDB0000094", "HMDB0000193"],
        })

        score, details = scorer.score_metabolomics(
            flux_row, metabolites, "missing_vs_comparison"
        )

        self.assertEqual(score, 0.0)
        self.assertEqual(details["evidence_mode"], "no_data")
        self.assertEqual(details["n_sig"], 0)
        self.assertEqual(details["data_reason"], "missing_comparison_effect")


if __name__ == "__main__":
    unittest.main()
