#!/usr/bin/env python3
"""多批次与跨实验集成的公共接口回归测试。"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import yaml

from omicsweave.core.project import project_dir
from omicsweave.core.statistics.meta_analysis import random_effects_meta
from omicsweave.integration.methods.batch_effect_meta import (
    meta_analyze_feature_batches,
)
from omicsweave.integration.methods.flux_omics_scorer import FluxOmicsScorer
from omicsweave.integration.loaders import (
    MFALoader,
    MetabolomicsLoader,
    ProteomicsLoader,
    TranscriptomicsLoader,
)
from omicsweave.integration.mappings import FLUX_GENE_MAPPING
from omicsweave.integration.project_integration import run_project_integration
from omicsweave.integration.provenance import build_input_audit
from omicsweave.mfa.df_analysis import (
    _compute_biological_net_flux,
    differential_flux,
)
from omicsweave.mfa.mfa_integration import (
    within_protocol_meta_analysis as mfa_within_protocol_meta_analysis,
)
from omicsweave.transcriptomics.transcriptomics_integration import (
    within_protocol_meta_analysis as tx_within_protocol_meta_analysis,
)
from omicsweave.integration.harmonization import (
    align_assay_samples,
    correct_batch_preserving_condition,
    diagnose_residual_batch,
    normalize_assay,
    run_harmonization,
)
from omicsweave.integration.config import load_multiomics_config


class IntegratedResultLoaderTest(unittest.TestCase):
    def test_mfa_loader_prefers_integrated_batch_meta_and_keeps_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            meta_dir = root / "integration" / "tracer_harmonized" / "glucose"
            meta_dir.mkdir(parents=True)
            pd.DataFrame({
                "flux_id": ["liver--HEX_c"],
                "tissue": ["liver"],
                "reaction": ["HEX_c"],
                "log2FC_pooled": [1.2],
                "diff_pooled": [2.4],
                "se_pooled": [0.4],
                "cohen_d_pooled": [1.1],
                "I_squared": [12.5],
                "tau_squared": [0.01],
                "n_batches": [2],
                "batch_ids": ["b1;b2"],
                "meta_note": ["random_effects_meta"],
            }).to_csv(meta_dir / "differential_flux_meta_case_vs_control.csv", index=False)

            cfg = SimpleNamespace(
                mfa_results=root,
                mfa_tissue_map={},
                integration={
                    "inputs": {
                        "mfa": {"mode": "integrated", "protocol_id": "glucose"},
                    }
                },
            )

            loaded = MFALoader(cfg).load()
            row = loaded["diff"]["case_vs_control"].iloc[0]

            self.assertEqual(row["flux_id"], "HEX_c")
            self.assertAlmostEqual(row["flux_log2fc"], 1.2)
            self.assertEqual(row["n_batches"], 2)
            self.assertEqual(row["batch_ids"], "b1;b2")
            self.assertEqual(row["source_mode"], "within_protocol_meta")
            self.assertAlmostEqual(row["heterogeneity_i2"], 12.5)

    def test_transcriptomics_loader_reads_integrated_batch_meta(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            meta_dir = root / "integration" / "within_protocol" / "glucose" / "liver"
            meta_dir.mkdir(parents=True)
            pd.DataFrame({
                "gene": ["G6pc"],
                "pooled_log2FoldChange": [-1.4],
                "pooled_se": [0.3],
                "pooled_pvalue": [0.002],
                "pooled_fdr": [0.01],
                "i2": [0.18],
                "tau_squared": [0.02],
                "n_batches": [3],
                "batch_ids": ["b1,b2,b3"],
                "direction_consistent": [True],
            }).to_csv(meta_dir / "de_case_vs_control_integrated.csv", index=False)

            cfg = SimpleNamespace(
                transcriptomics_results=root,
                tissues=["liver"],
                comparison_map={"case_vs_control": "case_vs_control"},
                transcriptomics_comparison_map={"case_vs_control": "case_vs_control"},
                mfa_tissue_map={},
                integration={
                    "inputs": {
                        "transcriptomics": {
                            "mode": "integrated",
                            "protocol_id": "glucose",
                        },
                    }
                },
            )

            loaded = TranscriptomicsLoader(cfg).load()
            row = loaded["de"]["case_vs_control"]["liver"].iloc[0]

            self.assertEqual(row["gene"], "G6pc")
            self.assertAlmostEqual(row["log2fc"], -1.4)
            self.assertAlmostEqual(row["lfcSE"], 0.3)
            self.assertEqual(row["n_batches"], 3)
            self.assertEqual(row["source_mode"], "within_protocol_meta")

    def test_transcriptomics_pathways_use_canonical_comparison_keys(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            gsea_dir = root / "05_result_plot" / "5_gsea"
            gsea_dir.mkdir(parents=True)
            pd.DataFrame({
                "Term": ["Library__Glycolysis"],
                "FDR q-val": [0.01],
            }).to_csv(
                gsea_dir / "gsea_liver_native_case_vs_control_results.csv",
                index=False,
            )
            pd.DataFrame({
                "Term": ["Library__Unconfigured pathway"],
                "FDR q-val": [0.01],
            }).to_csv(
                gsea_dir / "gsea_liver_unconfigured_vs_control_results.csv",
                index=False,
            )
            cfg = SimpleNamespace(
                transcriptomics_results=root,
                tissues=["liver"],
                comparison_map={"case_vs_control": "native_case_vs_control"},
                transcriptomics_comparison_map={
                    "case_vs_control": "native_case_vs_control"
                },
                mfa_tissue_map={},
                integration={},
            )

            loaded = TranscriptomicsLoader(cfg).load()

            self.assertIn("case_vs_control", loaded["gsea"]["liver"])
            self.assertNotIn("native_case_vs_control", loaded["gsea"]["liver"])
            self.assertNotIn("unconfigured_vs_control", loaded["gsea"]["liver"])

    def test_transcriptomics_loader_backfills_audited_sample_counts(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            de_dir = root / "04_de_results" / "liver"
            de_dir.mkdir(parents=True)
            pd.DataFrame({
                "gene": ["G6pc"],
                "log2FoldChange": [1.1],
                "lfcSE": [0.2],
                "pvalue": [0.01],
                "padj": [0.02],
            }).to_csv(de_dir / "de_native_case_vs_control.csv", index=False)
            pd.DataFrame({
                "sample": ["m1_liver", "m2_liver", "m3_liver", "m4_liver", "m5_liver"],
                "mouse_id": ["m1", "m2", "m3", "m4", "m5"],
                "tissue": ["liver"] * 5,
                "condition": ["native_case", "native_case", "native_case", "control", "control"],
            }).to_csv(root / "sample_metadata.csv", index=False)
            cfg = SimpleNamespace(
                transcriptomics_results=root,
                tissues=["liver"],
                comparison_map={"case_vs_control": "native_case_vs_control"},
                transcriptomics_comparison_map={
                    "case_vs_control": "native_case_vs_control"
                },
                mfa_tissue_map={},
                integration={},
            )

            row = TranscriptomicsLoader(cfg).load()["de"]["case_vs_control"]["liver"].iloc[0]

            self.assertEqual(row["n_case"], 3)
            self.assertEqual(row["n_control"], 2)
            self.assertEqual(row["n_eff"], 5)
            self.assertEqual(row["n_eff_source"], "sample_metadata")


class RandomEffectsMetaTest(unittest.TestCase):
    def test_random_effects_reports_heterogeneity_and_source_labels(self) -> None:
        result = random_effects_meta(
            effects=[0.5, 1.5, 1.0],
            standard_errors=[0.1, 0.1, 0.1],
            source_ids=["batch_a", "batch_b", "batch_c"],
        )

        self.assertEqual(result.k, 3)
        self.assertGreater(result.tau_squared, 0.0)
        self.assertGreater(result.i_squared, 50.0)
        self.assertEqual(result.source_ids, ("batch_a", "batch_b", "batch_c"))
        self.assertAlmostEqual(result.pooled_effect, 1.0, places=6)
        self.assertGreater(result.pooled_se, 0.1)


class WithinOmicsBatchMetaTest(unittest.TestCase):
    def test_missing_uncertainty_is_retained_but_not_meta_eligible(self) -> None:
        result = meta_analyze_feature_batches(
            {
                "batch_a": pd.DataFrame({"feature": ["F1"], "effect": [0.8]}),
                "batch_b": pd.DataFrame({"feature": ["F1"], "effect": [1.2]}),
            },
            feature_column="feature",
            effect_column="effect",
            se_candidates=["se"],
            pvalue_candidates=["pvalue"],
        )

        row = result.iloc[0]
        self.assertAlmostEqual(row["pooled_effect"], 1.0)
        self.assertTrue(pd.isna(row["pooled_se"]))
        self.assertFalse(bool(row["meta_eligible"]))
        self.assertEqual(row["meta_note"], "not_meta_eligible_missing_uncertainty")
        self.assertEqual(row["n_batches_observed"], 2)
        self.assertEqual(row["n_batches_eligible"], 0)
        self.assertEqual(row["batch_ids"], "batch_a;batch_b")

    def test_reported_and_approximated_se_are_auditable(self) -> None:
        result = meta_analyze_feature_batches(
            {
                "batch_a": pd.DataFrame({
                    "feature": ["F1"], "effect": [0.8], "se": [0.2],
                }),
                "batch_b": pd.DataFrame({
                    "feature": ["F1"], "effect": [1.0], "pvalue": [0.01],
                }),
            },
            feature_column="feature",
            effect_column="effect",
            se_candidates=["se"],
            pvalue_candidates=["pvalue"],
        )

        row = result.iloc[0]
        self.assertTrue(bool(row["meta_eligible"]))
        self.assertEqual(row["meta_note"], "random_effects_meta")
        self.assertEqual(row["n_batches_eligible"], 2)
        self.assertIn("reported", row["se_source"])
        self.assertIn("pvalue_normal_approximation", row["se_source"])

    @staticmethod
    def _write_proteomics_batch(root: Path, effect: float, pvalue: float) -> None:
        tissue_dir = root / "01_dp_results" / "Liver"
        tissue_dir.mkdir(parents=True)
        pd.DataFrame({
            "Protein": ["G6pc"],
            "Log2FC_case_vs_control": [effect],
            "Pval_case_vs_control": [pvalue],
        }).to_csv(tissue_dir / "differential_proteins.csv", index=False)
        (tissue_dir / "analysis_metadata.json").write_text(
            '{"pipeline_version": "6.0"}', encoding="utf-8",
        )

    @staticmethod
    def _write_metabolomics_batch(root: Path, effect: float, pvalue: float) -> None:
        tissue_dir = root / "00_dm_results" / "liver"
        tissue_dir.mkdir(parents=True)
        pd.DataFrame({
            "Index": ["glucose"],
            "Compounds": ["Glucose"],
            "Log2FC_case_vs_control": [effect],
            "Pval_case_vs_control": [pvalue],
        }).to_csv(tissue_dir / "differential_analysis.csv", index=False)

    def test_proteomics_loader_meta_analyzes_configured_batches(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            batch_a, batch_b = root / "prot_a", root / "prot_b"
            self._write_proteomics_batch(batch_a, 0.8, 0.01)
            self._write_proteomics_batch(batch_b, 1.0, 0.02)
            cfg = SimpleNamespace(
                proteomics_results=None,
                tissues=["liver"],
                comparisons=["case_vs_control"],
                proteomics_comparison_map={},
                integration={"inputs": {"proteomics": {
                    "mode": "integrated",
                    "batches": [
                        {"batch_id": "p1", "results_dir": batch_a},
                        {"batch_id": "p2", "results_dir": batch_b},
                    ],
                }}},
            )

            row = ProteomicsLoader(cfg).load()["diff"]["liver"].iloc[0]

            self.assertEqual(row["protein"], "G6pc")
            self.assertEqual(row["n_batches_case_vs_control"], 2)
            self.assertEqual(row["batch_ids_case_vs_control"], "p1;p2")
            self.assertEqual(row["source_mode_case_vs_control"], "within_omics_batch_meta")
            self.assertTrue(pd.notna(row["se_case_vs_control"]))

    def test_proteomics_loader_reads_v6_prerank_pathways(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            tissue_dir = root / "01_dp_results" / "liver"
            tissue_dir.mkdir(parents=True)
            pd.DataFrame({
                "Protein": ["G6pc"],
                "Log2FC_case_vs_control": [1.0],
                "Pval_case_vs_control": [0.01],
                "Padj_case_vs_control": [0.02],
                "Evidence_case_vs_control": ["inferential"],
            }).to_csv(tissue_dir / "differential_proteins.csv", index=False)
            (tissue_dir / "analysis_metadata.json").write_text(
                '{"pipeline_version": "6.0"}', encoding="utf-8",
            )
            gsea_dir = (
                root / "02_result_plot" / "liver" / "gsea_prerank"
                / "case_vs_control" / "kegg"
            )
            gsea_dir.mkdir(parents=True)
            pd.DataFrame({
                "Term": ["Glycolysis"], "NES": [1.7], "FDR q-val": [0.01],
            }).to_csv(gsea_dir / "gseapy.gene_set.prerank.report.csv", index=False)
            cfg = SimpleNamespace(
                proteomics_results=root,
                tissues=["liver"],
                comparisons=["case_vs_control"],
                proteomics_comparison_map={},
                integration={},
            )

            payload = ProteomicsLoader(cfg).load()

            self.assertEqual(
                payload["diff"]["liver"].loc[0, "evidence_case_vs_control"],
                "inferential",
            )
            pathway = payload["pathway"]["liver"]
            self.assertIsNotNone(pathway)
            assert pathway is not None
            self.assertEqual(pathway.loc[0, "comparison"], "case_vs_control")
            self.assertEqual(pathway.loc[0, "source_mode"], "proteomics_prerank_gsea")

    def test_metabolomics_loader_meta_analyzes_configured_batches(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            batch_a, batch_b = root / "metab_a", root / "metab_b"
            self._write_metabolomics_batch(batch_a, -0.8, 0.01)
            self._write_metabolomics_batch(batch_b, -1.0, 0.02)
            cfg = SimpleNamespace(
                metabolomics_results=None,
                tissues=["liver"],
                comparisons=["case_vs_control"],
                metabolomics_comparison_map={},
                integration={"inputs": {"metabolomics": {
                    "mode": "integrated",
                    "feature_column": "compound",
                    "batches": [
                        {"batch_id": "m1", "results": batch_a},
                        {"batch_id": "m2", "results": batch_b},
                    ],
                }}},
            )

            row = MetabolomicsLoader(cfg).load()["diff"]["liver"].iloc[0]

            self.assertEqual(row["compound"], "Glucose")
            self.assertEqual(row["n_batches_case_vs_control"], 2)
            self.assertEqual(row["batch_ids_case_vs_control"], "m1;m2")
            self.assertEqual(row["source_mode_case_vs_control"], "within_omics_batch_meta")
            self.assertTrue(pd.notna(row["se_case_vs_control"]))

    def test_batch_provenance_reaches_proteomics_evidence_details(self) -> None:
        cfg = SimpleNamespace(
            experiment_id="experiment_a",
            estimand={"id": "diet_14w"},
            scoring={},
        )
        scorer = FluxOmicsScorer(cfg)
        flux_row = pd.Series({"flux_id": "G6Pase_c", "flux_log2fc": 1.0})
        proteins = pd.DataFrame({
            "protein": ["G6pc"],
            "log2fc_case_vs_control": [0.9],
            "pval_case_vs_control": [0.01],
            "padj_case_vs_control": [0.02],
            "se_case_vs_control": [0.2],
            "n_batches_case_vs_control": [2],
            "batch_ids_case_vs_control": ["p1;p2"],
            "heterogeneity_i2_case_vs_control": [15.0],
            "heterogeneity_tau2_case_vs_control": [0.01],
            "source_mode_case_vs_control": ["within_omics_batch_meta"],
            "protocol_id_case_vs_control": ["dia_v1"],
        })

        _, details = scorer.score_proteomics(
            flux_row, proteins, "case_vs_control"
        )

        self.assertAlmostEqual(details["se"], 0.2)
        self.assertEqual(details["n_batches"], 2)
        self.assertEqual(details["batch_ids"], "p1;p2")
        self.assertAlmostEqual(details["heterogeneity_i2"], 15.0)
        self.assertEqual(details["source_mode"], "within_omics_batch_meta")
        self.assertEqual(details["protocol_id"], "dia_v1")

    def test_proteomics_scoring_prefers_fdr_over_raw_pvalue(self) -> None:
        scorer = FluxOmicsScorer(SimpleNamespace(
            experiment_id="experiment_a", estimand={"id": "diet_14w"}, scoring={},
        ))
        flux_row = pd.Series({"flux_id": "G6Pase_c", "flux_log2fc": 1.0})
        proteins = pd.DataFrame({
            "protein": ["G6pc"],
            "log2fc_case_vs_control": [1.2],
            "pval_case_vs_control": [0.0001],
            "padj_case_vs_control": [0.5],
            "evidence_case_vs_control": ["inferential"],
        })

        score, details = scorer.score_proteomics(
            flux_row, proteins, "case_vs_control",
        )

        self.assertEqual(score, 0.0)
        self.assertEqual(details["n_sig"], 0)
        self.assertEqual(details["evidence_mode"], "no_significant")


class BiologicalInferenceAndBatchMetaTest(unittest.TestCase):
    def test_each_mouse_contributes_one_flux_estimate(self) -> None:
        flux_df = pd.DataFrame(
            {"liver--HEX_c": [1.0, 3.0]},
            index=["mouse_1", "mouse_2"],
        )
        all_fits = {
            "mouse_1": pd.DataFrame({"liver--HEX_c": [0.8, 1.0, 1.2]}),
            "mouse_2": pd.DataFrame({"liver--HEX_c": [2.5, 3.0, 3.5, 4.0]}),
        }

        values = _compute_biological_net_flux(
            all_fits,
            flux_df,
            ["mouse_1", "mouse_2"],
            "liver",
            "HEX_c",
            use_net_flux=True,
            iqr_filter=False,
        )

        self.assertEqual(len(values), 2)
        self.assertAlmostEqual(values[0], 1.0)
        self.assertAlmostEqual(values[1], 3.25)

    def test_differential_flux_exports_log2fc_uncertainty_and_stability(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            sample_ids = [f"case_{i}" for i in range(3)] + [f"control_{i}" for i in range(3)]
            best_values = [4.0, 5.0, 6.0, 1.0, 2.0, 2.5]
            flux_df = pd.DataFrame({"liver--HEX_c": best_values}, index=sample_ids)
            all_fits = {
                sample_id: pd.DataFrame({"liver--HEX_c": [value - 0.1, value, value + 0.1]})
                for sample_id, value in zip(sample_ids, best_values, strict=True)
            }
            sample_conditions = {
                sample_id: ("case" if sample_id.startswith("case") else "control")
                for sample_id in sample_ids
            }

            result = differential_flux(
                flux_df,
                all_fits,
                sample_conditions,
                ["case", "control"],
                ["liver"],
                Path(temp_dir),
                comparisons=[("case", "control")],
                analysis_config={
                    "annotation": {"require_standard_reaction_ids": False},
                    "differential_flux": {
                        "statistic": "ttest",
                        "stability": {
                            "enabled": True,
                            "n_bootstrap": 100,
                            "n_permutations": 100,
                            "random_state": 7,
                            "use_for_selection": False,
                        },
                    },
                },
            )

            row = result["CASE_vs_CONTROL"].iloc[0]
            self.assertEqual(row["inference_unit"], "biological_replicate")
            self.assertTrue(pd.notna(row["log2fc_se"]))
            self.assertGreater(row["log2fc_se"], 0.0)
            self.assertTrue(pd.notna(row["hedges_g"]))
            self.assertTrue(pd.notna(row["selection_probability"]))
            self.assertTrue(pd.notna(row["direction_stability"]))
            self.assertEqual(row["stability_evidence_mode"], "biological_replicate")

    def test_within_protocol_meta_uses_log2fc_and_matching_se_scale(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            for batch_id, effect, se in (("b1", 0.8, 0.2), ("b2", 1.2, 0.3)):
                out = (
                    root / "runs" / "glucose" / batch_id / "01_df_results"
                    / "liver" / "differential_flux"
                )
                out.mkdir(parents=True)
                pd.DataFrame({
                    "flux_id": ["liver--HEX_c"],
                    "tissue": ["liver"],
                    "reaction": ["HEX_c"],
                    "diff": [2.0 + effect],
                    "diff_se": [0.5],
                    "log2FC": [effect],
                    "log2fc_se": [se],
                    "cohen_d": [1.0],
                    "hedges_g": [0.9],
                    "is_significant": [True],
                }).to_csv(out / "differential_flux_CASE_vs_CONTROL.csv", index=False)

            result = mfa_within_protocol_meta_analysis(
                root,
                "glucose",
                ["b1", "b2"],
                {"within_protocol": {"min_batches": 2}},
            )

            row = result["CASE_vs_CONTROL"].iloc[0]
            self.assertEqual(row["effect_scale"], "log2FC")
            # 两批存在异质性，必须使用随机效应而非固定效应权重。
            self.assertAlmostEqual(row["log2FC_pooled"], 0.9375, places=6)
            self.assertAlmostEqual(row["se_pooled"], 0.189983551, places=6)
            self.assertTrue(pd.notna(row["pooled_pvalue"]))
            self.assertEqual(row["n_batches"], 2)
            self.assertEqual(row["meta_note"], "random_effects_meta")

    def test_transcriptomics_meta_uses_percent_i2_and_audits_eligible_batches(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            batch_rows = {
                "b1": pd.DataFrame({
                    "gene": ["GeneMeta", "GenePartial"],
                    "log2FoldChange": [0.8, 1.0],
                    "lfcSE": [0.2, 0.25],
                    "pvalue": [0.001, 0.002],
                    "padj": [0.01, 0.01],
                }),
                "b2": pd.DataFrame({
                    "gene": ["GeneMeta", "GenePartial"],
                    "log2FoldChange": [1.2, 1.3],
                    "lfcSE": [0.3, float("nan")],
                    "pvalue": [0.001, float("nan")],
                    "padj": [0.01, float("nan")],
                }),
            }
            for batch_id, table in batch_rows.items():
                de_dir = root / "runs" / "glucose" / batch_id / "04_de_results" / "liver"
                de_dir.mkdir(parents=True)
                table.to_csv(de_dir / "de_case_vs_control.csv", index=False)

            result = tx_within_protocol_meta_analysis(
                root,
                "glucose",
                ["b1", "b2"],
                {
                    "within_protocol": {
                        "method": "random_effects_meta",
                        "min_batches": 2,
                        # 旧配置比例仍应解释为 50%。
                        "heterogeneity_threshold": 0.5,
                    }
                },
            )

            table = result[("liver", "case_vs_control")].set_index("gene")
            meta = table.loc["GeneMeta"]
            partial = table.loc["GenePartial"]

            self.assertEqual(meta["n_batches_observed"], 2)
            self.assertEqual(meta["n_batches_eligible"], 2)
            self.assertTrue(bool(meta["meta_eligible"]))
            self.assertGreaterEqual(meta["i2"], 0.0)
            self.assertLessEqual(meta["i2"], 100.0)
            self.assertEqual(meta["effect_scale"], "log2FoldChange")
            self.assertIn(meta["meta_note"], {"random_effects_meta", "high_heterogeneity"})

            self.assertEqual(partial["n_batches_observed"], 2)
            self.assertEqual(partial["n_batches_eligible"], 1)
            self.assertFalse(bool(partial["meta_eligible"]))
            self.assertEqual(partial["n_batches"], 1)
            self.assertEqual(partial["meta_note"], "single_batch_passthrough")
            self.assertEqual(partial["eligible_batch_ids"], "b1")


class HarmonizationTest(unittest.TestCase):
    def test_residual_batch_qc_decreases_after_condition_preserving_correction(self) -> None:
        samples = [f"s{i}" for i in range(1, 9)]
        metadata = pd.DataFrame({
            "condition": ["control", "control", "case", "case"] * 2,
            "batch_id": ["b1"] * 4 + ["b2"] * 4,
        }, index=samples)
        condition = np.array([0, 0, 1, 1] * 2, dtype=float)
        batch = np.array([0] * 4 + [1] * 4, dtype=float)
        matrix = pd.DataFrame(
            np.vstack([
                2.0 * condition + 8.0 * batch + np.linspace(0.0, 0.07, 8),
                5.0 * condition + 4.0 * batch + np.linspace(0.0, 0.14, 8),
                6.0 * batch + np.linspace(0.0, 0.21, 8),
            ]),
            index=["f1", "f2", "f3"],
            columns=samples,
        )

        before = diagnose_residual_batch(matrix, metadata, assay="proteomics", stage="normalized")
        corrected = correct_batch_preserving_condition(matrix, metadata)
        after = diagnose_residual_batch(corrected, metadata, assay="proteomics", stage="corrected")

        self.assertEqual(before.summary["status"], "assessable")
        self.assertGreater(before.summary["feature_batch_partial_r2_p90"], 0.8)
        self.assertLess(
            after.summary["feature_batch_partial_r2_p90"],
            before.summary["feature_batch_partial_r2_p90"],
        )
        self.assertGreaterEqual(
            after.summary["suggested_batch_reliability"],
            before.summary["suggested_batch_reliability"],
        )
        self.assertIn("PC1", before.pca_scores.columns)
        self.assertEqual(set(before.feature_qc["feature"]), {"f1", "f2", "f3"})

    def test_union_alignment_marks_missing_assays_and_rejects_metadata_conflict(self) -> None:
        assays = {
            "transcriptomics": pd.DataFrame(
                {"rna_a": [1.0], "rna_b": [2.0]}, index=["g1"]
            ),
            "proteomics": pd.DataFrame({"prot_a": [4.0]}, index=["p1"]),
        }
        sample_map = pd.DataFrame({
            "assay": ["transcriptomics", "transcriptomics", "proteomics"],
            "assay_sample_id": ["rna_a", "rna_b", "prot_a"],
            "primary_sample_id": ["mouse_1", "mouse_2", "mouse_1"],
            "batch_id": ["b1", "b1", "b1"],
            "condition": ["case", "control", "case"],
        })

        result = align_assay_samples(assays, sample_map, min_assays=1)

        self.assertEqual(result.primary_samples, ("mouse_1", "mouse_2"))
        self.assertEqual(list(result.assays["proteomics"].columns), ["mouse_1", "mouse_2"])
        self.assertTrue(pd.isna(result.assays["proteomics"].loc["p1", "mouse_2"]))
        prot_qc = result.missingness[result.missingness["assay"] == "proteomics"].iloc[0]
        self.assertAlmostEqual(prot_qc["missing_fraction"], 0.5)

        conflicted = sample_map.copy()
        conflicted.loc[conflicted["assay"] == "proteomics", "condition"] = "control"
        with self.assertRaisesRegex(ValueError, "inconsistent condition"):
            align_assay_samples(assays, conflicted, min_assays=1)

        missing_batch = sample_map.copy()
        missing_batch.loc[0, "batch_id"] = np.nan
        with self.assertRaisesRegex(ValueError, "missing/empty values"):
            align_assay_samples(assays, missing_batch, min_assays=1)

    def test_sample_map_aligns_assays_and_reports_missingness(self) -> None:
        assays = {
            "transcriptomics": pd.DataFrame(
                {"rna_a": [1.0, 2.0], "rna_b": [2.0, 3.0]}, index=["g1", "g2"]
            ),
            "proteomics": pd.DataFrame(
                {"prot_a": [4.0, float("nan")], "prot_c": [5.0, 6.0]}, index=["p1", "p2"]
            ),
        }
        sample_map = pd.DataFrame({
            "assay": ["transcriptomics", "transcriptomics", "proteomics", "proteomics"],
            "assay_sample_id": ["rna_a", "rna_b", "prot_a", "prot_c"],
            "primary_sample_id": ["mouse_1", "mouse_2", "mouse_1", "mouse_3"],
            "batch_id": ["b1", "b1", "b2", "b2"],
            "condition": ["case", "control", "case", "control"],
        })

        result = align_assay_samples(assays, sample_map, min_assays=2)

        self.assertEqual(result.primary_samples, ("mouse_1",))
        self.assertEqual(list(result.assays["transcriptomics"].columns), ["mouse_1"])
        self.assertEqual(list(result.assays["proteomics"].columns), ["mouse_1"])
        prot_summary = result.missingness[result.missingness["assay"] == "proteomics"].iloc[0]
        self.assertAlmostEqual(prot_summary["missing_fraction"], 0.5)

    def test_batch_correction_preserves_condition_and_rejects_confounding(self) -> None:
        matrix = pd.DataFrame({
            "s1": [0.0, 10.0], "s2": [0.2, 10.2],
            "s3": [5.0, 15.0], "s4": [5.2, 15.2],
            "s5": [0.1, 10.1], "s6": [0.3, 10.3],
            "s7": [5.1, 15.1], "s8": [5.3, 15.3],
        }, index=["f1", "f2"])
        metadata = pd.DataFrame({
            "condition": ["control", "control", "case", "case"] * 2,
            "batch_id": ["b1"] * 4 + ["b2"] * 4,
        }, index=matrix.columns)

        corrected = correct_batch_preserving_condition(matrix, metadata)
        before_effect = matrix.loc["f1", metadata["condition"] == "case"].mean() - matrix.loc["f1", metadata["condition"] == "control"].mean()
        after_effect = corrected.loc["f1", metadata["condition"] == "case"].mean() - corrected.loc["f1", metadata["condition"] == "control"].mean()
        self.assertAlmostEqual(before_effect, after_effect, places=6)

        confounded = metadata.copy()
        confounded["condition"] = ["control"] * 4 + ["case"] * 4
        with self.assertRaisesRegex(ValueError, "confounded"):
            correct_batch_preserving_condition(matrix, confounded)

    def test_assay_specific_normalization_keeps_shape(self) -> None:
        matrix = pd.DataFrame({"s1": [1.0, 4.0], "s2": [2.0, 8.0]}, index=["f1", "f2"])
        normalized = normalize_assay(matrix, "proteomics")
        self.assertEqual(normalized.shape, matrix.shape)
        medians = normalized.median(axis=0)
        self.assertAlmostEqual(float(medians.iloc[0]), float(medians.iloc[1]))

    def test_configured_harmonization_writes_aligned_assays_and_qc(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            rna_path = root / "rna.csv"
            protein_path = root / "protein.csv"
            sample_map_path = root / "sample_map.csv"
            pd.DataFrame({
                "feature": ["GeneA", "GeneB"],
                "rna_1": [100, 50], "rna_2": [140, 45],
                "rna_3": [70, 80], "rna_4": [90, 70],
            }).to_csv(rna_path, index=False)
            pd.DataFrame({
                "feature": ["ProtA", "ProtB"],
                "prot_1": [1000, 500], "prot_2": [1300, 450],
                "prot_3": [700, 800], "prot_4": [850, 700],
            }).to_csv(protein_path, index=False)
            pd.DataFrame({
                "assay": ["transcriptomics"] * 4 + ["proteomics"] * 4,
                "assay_sample_id": [
                    "rna_1", "rna_2", "rna_3", "rna_4",
                    "prot_1", "prot_2", "prot_3", "prot_4",
                ],
                "primary_sample_id": ["m1", "m2", "m3", "m4"] * 2,
                "batch_id": ["b1", "b2", "b1", "b2"] * 2,
                "condition": ["case", "case", "control", "control"] * 2,
            }).to_csv(sample_map_path, index=False)
            config_path = root / "harmonization.yaml"
            config_path.write_text(yaml.safe_dump({
                "output_dir": str(root / "harmonized"),
                "sample_map": str(sample_map_path),
                "min_assays": 2,
                "batch_correction": {"enabled": True},
                "assays": {
                    "transcriptomics": {
                        "path": str(rna_path), "feature_column": "feature",
                    },
                    "proteomics": {
                        "path": str(protein_path), "feature_column": "feature",
                    },
                },
            }), encoding="utf-8")

            result = run_harmonization(config_path)

            self.assertEqual(result.primary_samples, ("m1", "m2", "m3", "m4"))
            out = root / "harmonized"
            self.assertTrue((out / "transcriptomics_harmonized.csv").exists())
            self.assertTrue((out / "proteomics_harmonized.csv").exists())
            self.assertTrue((out / "sample_metadata.csv").exists())
            self.assertTrue((out / "missingness_qc.csv").exists())
            self.assertTrue((out / "batch_qc.csv").exists())
            self.assertTrue((out / "batch_pca_scores.csv").exists())
            self.assertTrue((out / "batch_feature_qc.csv").exists())
            self.assertTrue((out / "harmonization_report.md").exists())
            self.assertEqual(len(list((out / "batch_pca").glob("*.png"))), 4)
            manifest = pd.read_csv(out / "harmonization_manifest.csv")
            self.assertEqual(set(manifest["assay"]), {"transcriptomics", "proteomics"})
            self.assertTrue(manifest["batch_corrected"].all())
            batch_qc = pd.read_csv(out / "batch_qc.csv")
            self.assertEqual(
                set(batch_qc["stage"]), {"normalized", "corrected"}
            )
            report = (out / "harmonization_report.md").read_text(encoding="utf-8")
            self.assertIn("Residual Batch Diagnostics", report)
            self.assertIn("not automatically applied", report)


class ConfigCompatibilityTest(unittest.TestCase):
    @unittest.skipUnless(
        (project_dir() / "experiments/20250425_cimr_wt_mice_dio/multiomics_integration_setting/config.yaml").exists(),
        "需要实验项目数据：设置 OMICSWEAVE_PROJECT_DIR 指向实验项目目录",
    )
    def test_legacy_config_loads_with_empty_integration_contract(self) -> None:
        cfg = load_multiomics_config(
            config_path="experiments/20250425_cimr_wt_mice_dio/multiomics_integration_setting/config.yaml",
            outdir_override=str(Path(tempfile.gettempdir()) / "multiomics_config_compat"),
        )
        self.assertIsInstance(cfg.estimand, dict)
        self.assertIsInstance(cfg.integration, dict)
        self.assertIsInstance(cfg.harmonization, dict)


class InputAuditTest(unittest.TestCase):
    def test_configured_batches_are_not_reported_as_actual_pooled_batches(self) -> None:
        config = SimpleNamespace(
            mfa_results=Path("mfa_results"),
            transcriptomics_results=None,
            proteomics_results=None,
            metabolomics_results=None,
            integration={"inputs": {"mfa": {
                "mode": "integrated",
                "batches": [
                    {"batch_id": "b1", "results": Path("run_1")},
                    {"batch_id": "b2", "results": Path("run_2")},
                ],
            }}},
        )
        loaded_mfa = {
            "case_vs_control": pd.DataFrame({
                "flux_id": ["HEX_c"],
                "n_batches": [1],
                "batch_ids": ["b2"],
                "source_mode": ["within_protocol_meta"],
            })
        }

        audit = build_input_audit(config, {"mfa": loaded_mfa})
        row = audit[audit["omics"] == "mfa"].iloc[0]

        self.assertEqual(row["configured_n_batches"], 2)
        self.assertEqual(row["configured_batch_ids"], "b1;b2")
        self.assertEqual(row["actual_n_batches"], 1)
        self.assertEqual(row["actual_batch_ids"], "b2")


class EvidenceProvenanceTest(unittest.TestCase):
    def test_evidence_table_retains_experiment_batch_and_estimand_identity(self) -> None:
        flux_id, mapping = next(
            (key, value) for key, value in FLUX_GENE_MAPPING.items() if value.get("genes")
        )
        config = SimpleNamespace(
            experiment_id="experiment_a",
            estimand={"id": "diet_effect_14w"},
            integration={},
            scoring={"evidence_profile": {"mfa_inference_unit": "biological_replicate"}},
        )
        scorer = FluxOmicsScorer(config)
        mfa = pd.DataFrame({
            "flux_id": [flux_id],
            "tissue": ["liver"],
            "flux_log2fc": [1.0],
            "flux_fdr": [0.01],
            "cohen_d": [1.2],
            "flux_se": [0.2],
            "n_eff": [8],
            "n_batches": [2],
            "batch_ids": ["mfa_1;mfa_2"],
            "protocol_id": ["glucose"],
            "source_mode": ["within_protocol_meta"],
            "heterogeneity_i2": [8.0],
        })
        tx = pd.DataFrame({
            "gene": [mapping["genes"][0]],
            "log2fc": [1.1],
            "padj": [0.01],
            "lfcSE": [0.25],
            "n_eff": [10],
            "n_batches": [2],
            "batch_ids": ["tx_1,tx_2"],
            "protocol_id": ["glucose"],
            "source_mode": ["within_protocol_meta"],
            "heterogeneity_i2": [4.0],
        })

        scorer.score_all(
            {"case_vs_control": mfa},
            {"case_vs_control": {"liver": tx}},
            {},
            {},
        )
        evidence = scorer.get_evidence_table()
        tx_row = evidence[evidence["layer"] == "transcriptomics"].iloc[0]

        self.assertEqual(tx_row["experiment_id"], "experiment_a")
        self.assertEqual(tx_row["estimand_id"], "diet_effect_14w")
        self.assertEqual(tx_row["batch_ids"], "tx_1,tx_2")
        self.assertEqual(tx_row["n_batches"], 2)
        self.assertEqual(tx_row["source_mode"], "within_protocol_meta")
        self.assertAlmostEqual(tx_row["heterogeneity_i2"], 4.0)


class ProjectIntegrationTest(unittest.TestCase):
    @staticmethod
    def _write_experiment(root: Path, experiment_id: str, effect: float) -> Path:
        result_dir = root / experiment_id / "multiomics_results"
        common_dir = result_dir / "common"
        common_dir.mkdir(parents=True)
        pd.DataFrame({
            "experiment_id": [experiment_id],
            "estimand_id": ["diet_14w"],
            "flux_id": ["HEX_c"],
            "tissue": ["liver"],
            "comparison": ["local_case_vs_control"],
            "layer": ["mfa"],
            "mechanism_aligned_effect": [effect],
            "effect": [effect],
            "se": [0.2],
            "qvalue": [0.01],
            "assessment_status": ["support"],
            "support_evidence": [0.8],
            "conflict_evidence": [0.0],
            "coverage_value": [1.0],
            "design_grade": ["B"],
        }).to_csv(common_dir / "evidence_table.csv", index=False)
        return result_dir

    def test_project_runner_meta_analyzes_comparable_experiments_and_writes_one_report(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            result_a = self._write_experiment(root, "experiment_a", 1.0)
            result_b = self._write_experiment(root, "experiment_b", 1.2)
            output_dir = root / "project_results"
            estimand = {
                "id": "diet_14w",
                "species": "mus_musculus",
                "case": "HFD",
                "control": "LFD",
                "intervention_time": "14w",
                "tracer": "u-13c_glucose",
                "mfa_network": "central_carbon_v1",
                "platforms": {"mfa": "isocor"},
                "primary_endpoint": "flux_log2fc",
            }
            config = {
                "project_id": "diet_replication",
                "output_dir": str(output_dir),
                "experiments": [
                    {
                        "experiment_id": "experiment_a",
                        "results_dir": str(result_a),
                        "estimand": estimand,
                        "comparison_map": {"local_case_vs_control": "HFD_vs_LFD"},
                    },
                    {
                        "experiment_id": "experiment_b",
                        "results_dir": str(result_b),
                        "estimand": estimand,
                        "comparison_map": {"local_case_vs_control": "HFD_vs_LFD"},
                    },
                ],
                "meta_analysis": {"min_experiments": 2, "fdr_threshold": 0.05},
            }
            config_path = root / "project_config.yaml"
            config_path.write_text(yaml.safe_dump(config), encoding="utf-8")

            result = run_project_integration(config_path)

            self.assertTrue(result.report_path.exists())
            self.assertEqual(len(result.meta_results), 1)
            meta_row = result.meta_results.iloc[0]
            self.assertEqual(meta_row["n_experiments"], 2)
            self.assertEqual(meta_row["experiment_ids"], "experiment_a;experiment_b")
            self.assertEqual(meta_row["canonical_comparison"], "HFD_vs_LFD")
            self.assertEqual(meta_row["replication_status"], "replicated")
            self.assertEqual(result.reaction_summary.iloc[0]["evidence_grade"], "Replicated")
            report_text = result.report_path.read_text(encoding="utf-8")
            self.assertIn("experiment_a", report_text)
            self.assertIn("experiment_b", report_text)
            self.assertIn("Replicated", report_text)

    def test_same_estimand_id_with_inconsistent_metadata_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            result_a = self._write_experiment(root, "experiment_a", 1.0)
            result_b = self._write_experiment(root, "experiment_b", 1.2)
            base_estimand = {
                "id": "diet_14w",
                "species": "mus_musculus",
                "case": "HFD",
                "control": "LFD",
                "intervention_time": "14w",
                "tracer": "u-13c_glucose",
                "mfa_network": "central_carbon_v1",
                "platforms": {"mfa": "isocor"},
                "primary_endpoint": "flux_log2fc",
            }
            inconsistent = {**base_estimand, "intervention_time": "8w"}
            config_path = root / "project_config.yaml"
            config_path.write_text(yaml.safe_dump({
                "project_id": "invalid_project",
                "output_dir": str(root / "out"),
                "experiments": [
                    {"experiment_id": "experiment_a", "results_dir": str(result_a), "estimand": base_estimand},
                    {"experiment_id": "experiment_b", "results_dir": str(result_b), "estimand": inconsistent},
                ],
            }), encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "intervention_time"):
                run_project_integration(config_path)

    def test_different_estimands_are_reported_but_not_meta_analyzed_together(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            result_a = self._write_experiment(root, "experiment_a", 1.0)
            result_b = self._write_experiment(root, "experiment_b", 1.2)
            estimand_a = {
                "id": "diet_14w", "species": "mus_musculus", "case": "HFD",
                "control": "LFD", "intervention_time": "14w", "tracer": "u-13c_glucose",
                "mfa_network": "central_carbon_v1", "platforms": {"mfa": "isocor"},
                "primary_endpoint": "flux_log2fc",
            }
            estimand_b = {**estimand_a, "id": "surgery_8w", "case": "SG", "control": "SHAM"}
            config_path = root / "project_config.yaml"
            config_path.write_text(yaml.safe_dump({
                "project_id": "non_comparable_project",
                "output_dir": str(root / "out"),
                "experiments": [
                    {"experiment_id": "experiment_a", "results_dir": str(result_a), "estimand": estimand_a},
                    {"experiment_id": "experiment_b", "results_dir": str(result_b), "estimand": estimand_b},
                ],
            }), encoding="utf-8")

            result = run_project_integration(config_path)

            self.assertTrue(result.meta_results.empty)
            self.assertFalse(result.compatibility_audit["comparable"].any())
            self.assertEqual(set(result.manifest["experiment_id"]), {"experiment_a", "experiment_b"})

    def test_missing_standard_error_is_explicitly_not_meta_eligible(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            result_dir = self._write_experiment(root, "experiment_a", 1.0)
            evidence_path = result_dir / "common" / "evidence_table.csv"
            evidence = pd.read_csv(evidence_path)
            evidence["se"] = float("nan")
            evidence.to_csv(evidence_path, index=False)
            estimand = {
                "id": "diet_14w", "species": "mus_musculus", "case": "HFD",
                "control": "LFD", "intervention_time": "14w", "tracer": "u-13c_glucose",
                "mfa_network": "central_carbon_v1", "platforms": {"mfa": "isocor"},
                "primary_endpoint": "flux_log2fc",
            }
            config_path = root / "project_config.yaml"
            config_path.write_text(yaml.safe_dump({
                "project_id": "missing_se_project",
                "output_dir": str(root / "out"),
                "experiments": [
                    {"experiment_id": "experiment_a", "results_dir": str(result_dir), "estimand": estimand},
                ],
            }), encoding="utf-8")

            result = run_project_integration(config_path)

            self.assertTrue(result.meta_results.empty)
            self.assertEqual(result.manifest.iloc[0]["n_meta_eligible_rows"], 0)
            self.assertFalse(result.evidence_long.iloc[0]["meta_eligible"])

    def test_opposite_but_nonsignificant_layers_are_no_signal_not_conflicting(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            result_dirs = []
            for experiment_id, mfa_effect, tx_effect in (
                ("experiment_a", 0.05, -0.05),
                ("experiment_b", 0.04, -0.04),
            ):
                result_dir = self._write_experiment(root, experiment_id, mfa_effect)
                evidence_path = result_dir / "common" / "evidence_table.csv"
                evidence = pd.read_csv(evidence_path)
                tx = evidence.iloc[0].copy()
                tx["layer"] = "transcriptomics"
                tx["mechanism_aligned_effect"] = tx_effect
                tx["effect"] = tx_effect
                evidence = pd.concat([evidence, tx.to_frame().T], ignore_index=True)
                evidence.to_csv(evidence_path, index=False)
                result_dirs.append(result_dir)
            estimand = {
                "id": "diet_14w", "species": "mus_musculus", "case": "HFD",
                "control": "LFD", "intervention_time": "14w", "tracer": "u-13c_glucose",
                "mfa_network": "central_carbon_v1",
                "platforms": {"mfa": "isocor", "transcriptomics": "kallisto"},
                "primary_endpoint": "flux_log2fc",
            }
            config_path = root / "project_config.yaml"
            config_path.write_text(yaml.safe_dump({
                "project_id": "no_signal_project",
                "output_dir": str(root / "out"),
                "experiments": [
                    {"experiment_id": "experiment_a", "results_dir": str(result_dirs[0]), "estimand": estimand},
                    {"experiment_id": "experiment_b", "results_dir": str(result_dirs[1]), "estimand": estimand},
                ],
                "meta_analysis": {"fdr_threshold": 0.05},
            }), encoding="utf-8")

            result = run_project_integration(config_path)

            self.assertEqual(result.reaction_summary.iloc[0]["evidence_grade"], "No-signal")


if __name__ == "__main__":
    unittest.main()
