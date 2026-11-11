#!/usr/bin/env python3
"""多组学证据画像的聚焦回归测试。"""

from __future__ import annotations

import tempfile
import unittest
import warnings
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

from omicsweave.core.project import project_dir
from omicsweave.integration.methods.evidence_profile import EvidenceProfileScorer
from omicsweave.integration.methods.cross_tissue_pathway_scorer import CrossTissuePathwayScorer
from omicsweave.integration.methods.flux_evidence import FluxEvidenceBuilder
from omicsweave.integration.methods.pathway_convergence import PathwayConvergenceAnalyzer
from omicsweave.core.statistics.stability import calibrate_biological_stability
from omicsweave.integration.mappings import FLUX_GENE_MAPPING
from omicsweave.integration.run_cross_tissue_integration import (
    run as run_cross_tissue_integration,
)


class EvidenceProfileTest(unittest.TestCase):
    def _base_record(self, protein_mode: str) -> dict:
        protein_details = {
            "evidence_mode": protein_mode,
            "n_proteins_checked": 2,
            "n_proteins_found": 2 if protein_mode == "no_significant" else 0,
            "n_sig": 0,
            "n_concordant": 0,
            "n_opposite": 0,
            "representative_effect": 0.1,
            "qvalue": 0.8,
            "se": 0.2,
            "n_eff": 8,
        }
        return {
            "flux_id": "v_test",
            "tissue": "liver",
            "comparison": "case_vs_control",
            "flux_log2fc": 1.2,
            "flux_fdr": 0.01,
            "flux_se": 0.2,
            "n_eff": 8,
            "gate_passed": True,
            "mfa_evidence_mode": "statistical_test",
            "mfa_details": {"evidence_mode": "statistical_test", "se": 0.2, "n_eff": 8},
            "tx_details": {
                "evidence_mode": "statistical_test",
                "n_genes_checked": 1,
                "n_genes_found": 1,
                "n_sig": 1,
                "n_concordant": 1,
                "n_opposite": 0,
                "support_effect": 1.1,
                "qvalue": 0.01,
                "se": 0.2,
                "n_eff": 8,
            },
            "prot_details": protein_details,
            "metab_details": {
                "evidence_mode": "no_data",
                "n_metabolites_checked": 0,
                "n_sig": 0,
            },
        }

    def test_missing_layer_reduces_coverage_not_available_layer_normalization(self) -> None:
        scorer = EvidenceProfileScorer({"mfa_inference_unit": "biological_replicate"})
        missing_records = scorer.build_records(self._base_record("no_data"))
        measured_records = scorer.build_records(self._base_record("no_significant"))
        missing = scorer.summarize(missing_records)
        measured = scorer.summarize(measured_records)

        self.assertLess(missing["coverage_score"], measured["coverage_score"])
        self.assertGreater(missing["support_score"], measured["support_score"])
        self.assertEqual(missing["evidence_grade"], "Convergent")

    def test_mfa_gate_failure_keeps_other_omics_in_new_profile(self) -> None:
        # 一个可给出方向证据的反应及其直接酶（近平衡、转运反应的酶量只记关联，不能用来测这条路径）
        flux_id, gene = "HEX_c", "Hk2"
        self.assertIn(gene, FLUX_GENE_MAPPING[flux_id]["genes"])
        config = SimpleNamespace(scoring={
            "evidence_profile": {"mfa_inference_unit": "mcmc_fit_points"},
        })
        builder = FluxEvidenceBuilder(config)
        flux_row = pd.Series({
            "flux_id": flux_id,
            "tissue": "liver",
            "flux_log2fc": 1.0,
            "flux_fdr": 0.5,
            "flux_pvalue": 0.5,
            "cohen_d": 0.2,
        })
        tx_df = pd.DataFrame({
            "gene": [gene],
            "log2fc": [1.4],
            "padj": [0.01],
            "lfcSE": [0.2],
            "n_eff": [8],
        })

        result = builder.build_record(
            flux_row, tx_df, None, None, "case_vs_control"
        )
        evidence = builder.evidence_profile_scorer.build_records(result)
        tx_record = next(row for row in evidence if row["layer"] == "transcriptomics")

        self.assertFalse(result["gate_passed"])
        self.assertGreater(result["support_score"], 0.0)
        self.assertEqual(tx_record["assessment_status"], "support")
        self.assertFalse({"legacy_total_score", "total_score", "confidence"} & set(result))

        scored = builder.build_all(
            {"case_vs_control": pd.DataFrame([flux_row])},
            {"case_vs_control": {"liver": tx_df}},
            {},
            {},
        )
        evidence_table = builder.get_evidence_table()
        self.assertEqual(len(scored), 1)
        self.assertEqual(len(evidence_table), 4)
        self.assertTrue({
            "effect", "se", "qvalue", "n_eff", "direction",
            "direction_relation", "assessment_status",
        }.issubset(evidence_table.columns))

    def test_opposite_direction_is_explicit_conflict(self) -> None:
        record = self._base_record("no_data")
        record["tx_details"].update({
            "n_concordant": 0,
            "n_opposite": 1,
            "support_effect": np.nan,
            "conflict_effect": 1.3,
            "direction_override": "opposite",
        })
        scorer = EvidenceProfileScorer({"mfa_inference_unit": "biological_replicate"})
        records = scorer.build_records(record)
        tx_record = next(row for row in records if row["layer"] == "transcriptomics")

        self.assertEqual(tx_record["assessment_status"], "opposite")
        self.assertEqual(tx_record["support_evidence"], 0.0)
        self.assertGreater(tx_record["conflict_evidence"], 0.0)

    def test_reviewed_layer_quality_factor_is_explicit_and_scales_reliability(self) -> None:
        base = self._base_record("no_data")
        default_scorer = EvidenceProfileScorer({
            "mfa_inference_unit": "biological_replicate",
        })
        reviewed_scorer = EvidenceProfileScorer({
            "mfa_inference_unit": "biological_replicate",
            "quality_reliability": {"transcriptomics": 0.5},
        })

        default_tx = next(
            row for row in default_scorer.build_records(base)
            if row["layer"] == "transcriptomics"
        )
        reviewed_tx = next(
            row for row in reviewed_scorer.build_records(base)
            if row["layer"] == "transcriptomics"
        )

        self.assertEqual(reviewed_tx["quality_reliability"], 0.5)
        self.assertAlmostEqual(
            reviewed_tx["data_reliability"],
            default_tx["data_reliability"] * 0.5,
        )

    def test_transcriptomics_mapping_quality_distinguishes_direct_and_regulatory_genes(self) -> None:
        config = SimpleNamespace(scoring={
            "evidence_profile": {"mfa_inference_unit": "biological_replicate"},
        })
        builder = FluxEvidenceBuilder(config)
        flux_row = pd.Series({
            "flux_id": "HEX_c",
            "tissue": "liver",
            "flux_log2fc": 1.0,
            "flux_fdr": 0.01,
            "gate_passed": True,
        })

        direct_details = builder.transcriptomics_evidence(
            flux_row,
            pd.DataFrame({
                "gene": ["Hk1"], "log2fc": [1.2],
                "padj": [0.01], "lfcSE": [0.2], "n_eff": [6],
            }),
        )
        indirect_details = builder.transcriptomics_evidence(
            flux_row,
            pd.DataFrame({
                "gene": ["Hif1a"], "log2fc": [1.2],
                "padj": [0.01], "lfcSE": [0.2], "n_eff": [6],
            }),
        )
        direct_record = builder.evidence_profile_scorer.build_records({
            "flux_id": "HEX_c", "tissue": "liver", "comparison": "case_vs_control",
            "flux_log2fc": 1.0, "gate_passed": True,
            "mfa_details": {}, "tx_details": direct_details,
            "prot_details": {"evidence_mode": "no_data"},
            "metab_details": {"evidence_mode": "no_data"},
        })[1]
        indirect_record = builder.evidence_profile_scorer.build_records({
            "flux_id": "HEX_c", "tissue": "liver", "comparison": "case_vs_control",
            "flux_log2fc": 1.0, "gate_passed": True,
            "mfa_details": {}, "tx_details": indirect_details,
            "prot_details": {"evidence_mode": "no_data"},
            "metab_details": {"evidence_mode": "no_data"},
        })[1]

        self.assertEqual(direct_details["direction_source"], "direct_gene")
        self.assertEqual(indirect_details["direction_source"], "indirect_gene")
        self.assertGreater(
            direct_record["mapping_quality"], indirect_record["mapping_quality"]
        )

    def test_pathway_aggregation_has_no_pandas_apply_deprecation(self) -> None:
        builder = FluxEvidenceBuilder(SimpleNamespace(scoring={}))
        scored = pd.DataFrame({
            "comparison": ["case_vs_control", "case_vs_control"],
            "tissue": ["liver", "liver"],
            "flux_id": ["reaction_a", "reaction_b"],
            "support_score": [0.2, 0.6],
            "conflict_score": [0.0, 0.1],
            "coverage_score": [0.5, 0.75],
            "evidence_grade": ["Exploratory", "Convergent"],
            "gate_passed": [True, False],
        })

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", FutureWarning)
            aggregated = builder.aggregate_by_pathway(scored)

        self.assertEqual(len(aggregated), 1)
        row = aggregated.iloc[0]
        self.assertEqual((row["n_reactions"], row["n_gated"], row["n_convergent"]), (2, 1, 1))
        self.assertEqual(row["top_reactions"], "reaction_b, reaction_a")
        self.assertAlmostEqual(row["mean_support_score"], 0.4)
        self.assertFalse(
            any(issubclass(item.category, FutureWarning) for item in caught)
        )


class StabilityTest(unittest.TestCase):
    def test_bootstrap_and_permutation_use_biological_replicates(self) -> None:
        result = calibrate_biological_stability(
            [2.0, 2.1, 1.9, 2.2],
            [0.0, 0.1, -0.1, 0.2],
            n_bootstrap=300,
            n_permutations=500,
            random_state=7,
        )
        self.assertEqual(result.evidence_mode, "biological_replicate")
        self.assertGreater(result.selection_probability, 0.8)
        self.assertGreater(result.direction_stability, 0.95)

        descriptive = calibrate_biological_stability([1.0], [0.0], random_state=7)
        self.assertEqual(descriptive.evidence_mode, "descriptive")
        self.assertTrue(np.isnan(descriptive.selection_probability))


class CrossTissueProfileTest(unittest.TestCase):
    def test_curated_pathway_uses_median_effect_and_exploratory_grade(self) -> None:
        scorer = CrossTissuePathwayScorer(tissues=["heart", "liver", "kidney"])
        result = scorer.score_pathway({
            "id": "P1",
            "name": "test pathway",
            "evidence": {
                "heart": {"direction": "UP", "transcriptomics": {"log2fc": 1.0}},
                "liver": {"direction": "UP", "transcriptomics": {"log2fc": 1.1}},
                "kidney": {"direction": "UP", "transcriptomics": {"log2fc": 10.0}},
            },
        })

        self.assertEqual(result["median_abs_log2fc"], 1.1)
        self.assertEqual(result["design_grade"], "D")
        self.assertEqual(result["evidence_grade"], "Exploratory")
        self.assertIn("legacy_total_score", result)

    @unittest.skipUnless(
        (project_dir() / "experiments/20250613_cimr_wt_mice_dio_glp/multiomics_integration_setting").is_dir(),
        "需要实验项目数据：设置 OMICSWEAVE_PROJECT_DIR 指向实验项目目录",
    )
    def test_curated_runner_uses_separate_dynamic_report(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir)
            report_path = run_cross_tissue_integration(
                "20250613_cimr_wt_mice_dio_glp",
                output_dir,
            )

            self.assertEqual(report_path.name, "cross_tissue_pathway_report.md")
            self.assertTrue(report_path.exists())
            self.assertFalse((output_dir / "multiomics_report.md").exists())
            report = report_path.read_text(encoding="utf-8")
            self.assertIn("# 跨组织 multiomics 通路整合报告", report)
            self.assertIn("HFD_Mazdu_vs_HFD_Saline", report)
            self.assertIn("来源标识", report)
            self.assertNotIn("底物重新分配五组织图谱", report)


class PathwayFusionTest(unittest.TestCase):
    def test_empty_fusion_schema_retains_comparison_identity(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            analyzer = PathwayConvergenceAnalyzer(Path(temp_dir))
            analyzer.analyze(
                tissues=["liver"],
                metab_pathway={},
                prot_pathway={},
                tx_gsea={},
            )
            self.assertIn("comparison", analyzer.fusion_df.columns)

    def test_bh_fdr_preserves_original_order(self) -> None:
        adjusted = PathwayConvergenceAnalyzer._bh_fdr(np.array([0.04, 0.001, 0.02]))
        np.testing.assert_allclose(adjusted, np.array([0.04, 0.003, 0.03]))

    def test_fisher_fusion_preserves_omics_sources(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            analyzer = PathwayConvergenceAnalyzer(Path(temp_dir))
            convergence = analyzer.analyze(
                tissues=["liver"],
                metab_pathway={
                    "liver": {
                        "case_vs_control": pd.DataFrame({
                            "Name": ["Pathway A", "Pathway B"],
                            "P-value": [0.001, 0.8],
                        })
                    }
                },
                prot_pathway={"liver": None},
                tx_gsea={
                    "liver": {
                        "case_vs_control": pd.DataFrame({
                            "Term": ["Library__Pathway A", "Library__Pathway C"],
                            "FDR q-val": [0.002, 0.9],
                        })
                    }
                },
            )

            self.assertFalse(convergence.empty)
            self.assertEqual(set(convergence["comparison"]), {"case_vs_control"})
            fused = analyzer.fusion_df
            pathway_a = fused[fused["pathway_term"] == "pathway a"].iloc[0]
            self.assertEqual(pathway_a["fusion_method"], "fisher")
            self.assertEqual(pathway_a["n_omics"], 2)
            self.assertIn("metabolomics", pathway_a["contributing_omics"])
            self.assertTrue((Path(temp_dir) / "tables" / "pathway_fusion.csv").exists())

    def test_pathway_fusion_does_not_mix_different_comparisons(self) -> None:
        analyzer = PathwayConvergenceAnalyzer(Path(tempfile.mkdtemp()))
        evidence = pd.DataFrame({
            "tissue": ["liver", "liver"],
            "comparison": ["diet_vs_control", "surgery_vs_sham"],
            "omics": ["metabolomics", "transcriptomics"],
            "pathway_term": ["glycolysis", "glycolysis"],
            "probability": [0.001, 0.002],
            "rank_fraction": [0.1, 0.1],
        })

        fused = analyzer._fuse_pathway_evidence(evidence)

        self.assertTrue(fused.empty)

    def test_exact_name_convergence_does_not_mix_different_comparisons(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            analyzer = PathwayConvergenceAnalyzer(Path(temp_dir))
            convergence = analyzer.analyze(
                tissues=["liver"],
                metab_pathway={
                    "liver": {
                        "diet_vs_control": pd.DataFrame({
                            "Name": ["Glycolysis"],
                            "P-value": [0.001],
                        })
                    }
                },
                prot_pathway={"liver": None},
                tx_gsea={
                    "liver": {
                        "surgery_vs_sham": pd.DataFrame({
                            "Term": ["Library__Glycolysis"],
                            "FDR q-val": [0.001],
                        })
                    }
                },
            )

            self.assertTrue(convergence.empty)
            empty_csv = pd.read_csv(
                Path(temp_dir) / "tables" / "pathway_fusion.csv"
            )
            self.assertIn("comparison", empty_csv.columns)


if __name__ == "__main__":
    unittest.main()
