"""代谢组学流程 P0-P3 改进项的聚焦回归测试。"""

from __future__ import annotations

import tempfile
import unittest
import warnings
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

from omicsweave.metabolomics.data_preprocessing import (
    PreprocessingOptions,
    apply_fdr_correction,
    classify_data_columns,
    preprocess_metabolomics_data,
)
from omicsweave.metabolomics.differential_analysis import (
    DifferentialAnalysisOptions,
    analyze_differential_metabolites,
)
from omicsweave.metabolomics.enrichment._run import calculate_pathway_enrichment
from omicsweave.metabolomics.experiment_config import (
    MetabolomicsConfig,
    generate_comparisons,
    load_metabolomics_config,
)
from omicsweave.metabolomics.gsea._run import build_ranked_metabolites
from omicsweave.metabolomics.kegg_pathways import get_pathway_name, is_pathway_in_scope
from omicsweave.metabolomics.plotting.qc_plots import fit_plsda
from omicsweave.metabolomics.plotting.heatmap import select_heatmap_features
from omicsweave.metabolomics.metabolomics_pipeline import run_pipeline
from omicsweave.metabolomics.generate_interpretation_report import generate_tissue_summary


class PreprocessingTests(unittest.TestCase):
    def test_column_classification_supports_project_ids_and_qc_samples(self) -> None:
        columns = pd.Index(
            [
                "Index", "Class-I", "A35_sq_adipose", "52-liver",
                "99_liver", "QC01", "QC_02",
            ]
        )

        classified = classify_data_columns(columns, known_mouse_ids={"A35", "52"})

        self.assertEqual(classified.biological, ["A35_sq_adipose", "52-liver"])
        self.assertEqual(classified.qc, ["QC01", "QC_02"])
        self.assertEqual(classified.metadata, ["Index", "Class-I"])
        self.assertEqual(classified.unmapped_biological, ["99_liver"])
        self.assertEqual(classified.mouse_ids["A35_sq_adipose"], "A35")
        self.assertEqual(classified.tissues["A35_sq_adipose"], "sq_adipose")

    def test_qc_rsd_and_detection_gate_filters_unreliable_features(self) -> None:
        frame = pd.DataFrame(
            {
                "Index": ["stable", "unstable", "poor-detection"],
                "A1_liver": [10.0, 10.0, 10.0],
                "A2_liver": [11.0, 11.0, 11.0],
                "B1_liver": [12.0, 12.0, 12.0],
                "B2_liver": [13.0, 13.0, 13.0],
                "QC01": [100.0, 10.0, 5.0],
                "QC02": [101.0, 100.0, 0.0],
                "QC03": [99.0, 1000.0, 0.0],
            }
        )
        samples = ["A1_liver", "A2_liver", "B1_liver", "B2_liver"]
        groups = {samples[0]: "A", samples[1]: "A", samples[2]: "B", samples[3]: "B"}

        result = preprocess_metabolomics_data(
            frame,
            samples,
            groups,
            qc_cols=["QC01", "QC02", "QC03"],
            options=PreprocessingOptions(
                normalization_method="none",
                qc_rsd_cutoff=0.30,
                qc_min_detection_rate=0.80,
            ),
        )

        self.assertEqual(result.processed_df["Index"].tolist(), ["stable"])
        metrics = result.qc_metrics.set_index(frame.loc[result.qc_metrics.index, "Index"])
        self.assertTrue(metrics.loc["stable", "QC_Pass"])
        self.assertFalse(metrics.loc["unstable", "QC_Pass"])
        self.assertFalse(metrics.loc["poor-detection", "QC_Pass"])
        self.assertEqual(result.summary["removed_by_qc"], 2)

    def test_group_detection_and_feature_minimum_imputation_preserve_on_off_signals(self) -> None:
        frame = pd.DataFrame(
            {
                "Index": ["group-specific", "high", "low"],
                "A1": [100.0, 1000.0, 10.0],
                "A2": [100.0, 0.0, 0.0],
                "A3": [100.0, 2000.0, 20.0],
                "B1": [0.0, 1000.0, 10.0],
                "B2": [0.0, 2000.0, 20.0],
                "B3": [0.0, 1000.0, 10.0],
            }
        )
        samples = ["A1", "A2", "A3", "B1", "B2", "B3"]
        groups = {sample: sample[0] for sample in samples}

        result = preprocess_metabolomics_data(
            frame,
            samples,
            groups,
            options=PreprocessingOptions(
                normalization_method="none",
                detection_rate=0.80,
                qc_filter_enabled=False,
            ),
        )
        processed = result.processed_df.set_index("Index")

        self.assertIn("group-specific", processed.index)
        self.assertAlmostEqual(processed.loc["high", "A2"], np.log2(500.0))
        self.assertAlmostEqual(processed.loc["low", "A2"], np.log2(5.0))
        self.assertEqual(result.summary["detection_rule"], "at_least_one_group")

    def test_pqn_removes_uniform_sample_loading_offset(self) -> None:
        frame = pd.DataFrame(
            {
                "Index": ["m1", "m2", "m3", "m4"],
                "Compounds": ["a", "b", "c", "d"],
                "A1_liver": [10.0, 20.0, 40.0, 80.0],
                "A2_liver": [20.0, 40.0, 80.0, 160.0],
                "B1_liver": [40.0, 80.0, 160.0, 320.0],
                "B2_liver": [5.0, 10.0, 20.0, 40.0],
            }
        )
        sample_cols = ["A1_liver", "A2_liver", "B1_liver", "B2_liver"]
        groups = {
            "A1_liver": "A",
            "A2_liver": "A",
            "B1_liver": "B",
            "B2_liver": "B",
        }

        result = preprocess_metabolomics_data(
            frame,
            sample_cols,
            groups,
            options=PreprocessingOptions(normalization_method="pqn"),
        )

        normalized = result.processed_df[sample_cols]
        for sample in sample_cols[1:]:
            np.testing.assert_allclose(
                normalized[sample].to_numpy(),
                normalized[sample_cols[0]].to_numpy(),
                atol=1e-12,
            )
        self.assertEqual(set(result.normalization_factors.index), set(sample_cols))
        self.assertEqual(result.summary["normalization_method"], "pqn")

    def test_fdr_method_is_honoured_and_nan_is_excluded(self) -> None:
        pvalues = np.array([0.01, 0.04, np.nan, 0.20])

        bh = apply_fdr_correction(pvalues, method="bh")
        bonferroni = apply_fdr_correction(pvalues, method="bonferroni")

        np.testing.assert_allclose(
            bh[[0, 1, 3]],
            np.array([0.03, 0.06, 0.20]),
            rtol=1e-12,
        )
        np.testing.assert_allclose(
            bonferroni[[0, 1, 3]],
            np.array([0.03, 0.12, 0.60]),
            rtol=1e-12,
        )
        self.assertTrue(np.isnan(bh[2]))
        with self.assertRaisesRegex(ValueError, "不支持"):
            apply_fdr_correction([0.01, 0.02], method="invented")


class ConfigurationTests(unittest.TestCase):
    def test_new_statistical_defaults_load_without_editing_old_configs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            groups_path = root / "groups.yaml"
            groups_path.write_text(
                "conditions:\n  treatment: {display: Treatment}\n  control: {display: Control, reference: true}\n"
                "mouse_groups:\n  '1': treatment\n  '2': control\n",
                encoding="utf-8",
            )
            config_path = root / "config.yaml"
            config_path.write_text(
                "experiment_id: test\n"
                "paths:\n"
                f"  sample_data_liver: '{(root / 'data.xlsx').as_posix()}'\n"
                f"  output_dir: '{(root / 'out').as_posix()}'\n"
                f"groups_file: '{groups_path.as_posix()}'\n",
                encoding="utf-8",
            )

            cfg = load_metabolomics_config(config_path=str(config_path))

        self.assertEqual(cfg.normalization_method, "pqn")
        self.assertEqual(cfg.gsea_permutation_num, 1000)
        self.assertEqual(cfg.heatmap_selection_method, "variance")
        self.assertEqual(cfg.batch_correction, "none")
        self.assertFalse(cfg.require_vip_for_significance)

    def test_auto_comparison_places_marked_reference_in_denominator(self) -> None:
        conditions = {
            "control": {"display": "Control", "reference": True},
            "treated": {"display": "Treated"},
        }

        self.assertEqual(
            generate_comparisons(conditions),
            [("treated", "control", "Treated_vs_Control")],
        )


class DifferentialAnalysisTests(unittest.TestCase):
    def test_primary_significance_uses_fdr_not_raw_pvalue(self) -> None:
        base = np.array([0.0, 1.0, 2.0])
        rows = []
        for index in range(10):
            delta = 2.6 if index == 0 else 0.0
            rows.append(
                {
                    "Index": f"m{index}",
                    "Compounds": f"metabolite-{index}",
                    "A1": base[0] + delta,
                    "A2": base[1] + delta,
                    "A3": base[2] + delta,
                    "B1": base[0],
                    "B2": base[1],
                    "B3": base[2],
                }
            )
        frame = pd.DataFrame(rows)
        groups = {"A1": "A", "A2": "A", "A3": "A", "B1": "B", "B2": "B", "B3": "B"}

        result = analyze_differential_metabolites(
            frame,
            groups,
            [("A", "B", "A_vs_B")],
            options=DifferentialAnalysisOptions(
                pvalue_cutoff=0.05,
                log2fc_cutoff=0.58,
                fdr_method="bh",
            ),
        )
        contrast = result.contrasts["A_vs_B"].set_index("Index")

        self.assertLess(contrast.loc["m0", "pvalue"], 0.05)
        self.assertGreater(contrast.loc["m0", "padj"], 0.05)
        self.assertTrue(contrast.loc["m0", "nominal_significant"])
        self.assertFalse(contrast.loc["m0", "significant"])
        self.assertEqual(result.design_report.loc[0, "test"], "welch_t_test")

    def test_vip_cutoff_is_reported_and_can_be_required(self) -> None:
        frame = pd.DataFrame(
            {
                "Index": ["low-vip", "high-vip"],
                "A1": [10.0, 20.0], "A2": [10.1, 20.1], "A3": [9.9, 19.9],
                "B1": [0.0, 0.0], "B2": [0.1, 0.1], "B3": [-0.1, -0.1],
            }
        )
        groups = {"A1": "A", "A2": "A", "A3": "A", "B1": "B", "B2": "B", "B3": "B"}
        result = analyze_differential_metabolites(
            frame,
            groups,
            [("A", "B", "A_vs_B")],
            options=DifferentialAnalysisOptions(
                vip_cutoff=1.0,
                require_vip_for_significance=True,
            ),
            contrast_vip={"A_vs_B": np.array([0.5, 1.5])},
        )
        contrast = result.contrasts["A_vs_B"].set_index("Index")

        self.assertTrue(contrast.loc["low-vip", "significant_univariate"])
        self.assertFalse(contrast.loc["low-vip", "significant"])
        self.assertTrue(contrast.loc["high-vip", "vip_pass"])
        self.assertTrue(contrast.loc["high-vip", "significant"])


class PathwayAnalysisTests(unittest.TestCase):
    def test_ora_uses_consistent_annotated_universe_and_adjusts_tests(self) -> None:
        rows = []
        for index in range(20):
            pathways = []
            if index < 5:
                pathways.append("ko00010")
            if 5 <= index < 10:
                pathways.append("ko00020")
            pathways.append("ko01100")
            if index < 8:
                pathways.append("ko05231")
            rows.append(
                {
                    "Index": f"m{index}",
                    "kegg_map": ",".join(pathways),
                    "Significant_A_vs_B": index < 4,
                }
            )
        rows.append(
            {
                "Index": "unannotated-sig",
                "kegg_map": "-",
                "Significant_A_vs_B": True,
            }
        )
        frame = pd.DataFrame(rows)

        result = calculate_pathway_enrichment(
            frame,
            significance_col="Significant_A_vs_B",
            min_set_size=3,
            max_set_size=10,
            pathway_scope="metabolism",
            fdr_method="bh",
        )

        self.assertEqual(set(result["ID"]), {"ko00010", "ko00020"})
        self.assertNotIn("ko01100", set(result["ID"]))
        self.assertNotIn("ko05231", set(result["ID"]))
        self.assertTrue({"P-value", "FDR", "OddsRatio"}.issubset(result.columns))
        self.assertTrue((result["SigTotal"] == 4).all())
        self.assertTrue((result["UniverseCount"] == 10).all())
        glycolysis = result.set_index("ID").loc["ko00010"]
        self.assertEqual(glycolysis["SigCount"], 4)
        self.assertLess(glycolysis["FDR"], 0.05)

    def test_lipid_and_xenobiotic_maps_are_named_without_changing_scope(self) -> None:
        # 广靶 + 脂质组注释中常见、此前显示为 name unavailable 的代谢通路
        expected = {
            "ko00563": "Glycosylphosphatidylinositol (GPI)-anchor biosynthesis",
            "ko00565": "Ether lipid metabolism",
            "ko00980": "Metabolism of xenobiotics by cytochrome P450",
            "ko01040": "Biosynthesis of unsaturated fatty acids",
        }
        for ko_id, name in expected.items():
            self.assertEqual(get_pathway_name(ko_id), name)
        # 名称只用于显示：是否进入统计由编号范围决定，缺名称的通路照常参与
        self.assertEqual(get_pathway_name("ko00999"), "ko00999 (name unavailable)")
        self.assertTrue(is_pathway_in_scope("ko00999"))


class GseaTests(unittest.TestCase):
    def test_ranked_metabolites_drop_nonfinite_values_and_prefer_test_statistic(self) -> None:
        frame = pd.DataFrame(
            {
                "Index": ["m1", "m2", "m3", "m4"],
                "Tstat_A_vs_B": [4.0, np.nan, -2.0, 0.0],
                "Pval_A_vs_B": [0.001, np.nan, 0.02, 1.0],
                "Log2FC_A_vs_B": [2.0, 3.0, -1.0, 0.0],
            }
        )

        ranked = build_ranked_metabolites(frame, "A_vs_B", metric="welch_t")

        self.assertEqual(ranked.index.tolist(), ["m1", "m4", "m3"])
        self.assertTrue(np.isfinite(ranked.to_numpy()).all())
        self.assertEqual(ranked.loc["m1"], 4.0)


class PlsdaTests(unittest.TestCase):
    def test_degenerate_plsda_is_marked_unvalidated(self) -> None:
        features = pd.DataFrame(
            np.ones((6, 4)),
            index=[f"sample-{index}" for index in range(6)],
            columns=[f"m{index}" for index in range(4)],
        )

        result = fit_plsda(
            features,
            ["A", "A", "A", "B", "B", "B"],
            n_permutations=5,
        )

        self.assertEqual(result.status, "degenerate_features")
        self.assertTrue(result.vip.isna().all())
        self.assertTrue(np.isnan(result.q2))

    def test_multiclass_plsda_is_invariant_to_class_name_order_and_validated(self) -> None:
        rng = np.random.default_rng(17)
        labels = np.repeat(["A", "B", "C"], 3)
        matrix = rng.normal(0, 0.15, size=(9, 12))
        matrix[labels == "A", :4] += 2.0
        matrix[labels == "B", 4:8] += 2.0
        matrix[labels == "C", 8:] += 2.0
        features = pd.DataFrame(
            matrix,
            index=[f"sample-{i}" for i in range(9)],
            columns=[f"m{i}" for i in range(12)],
        )
        renamed = pd.Series(labels).map({"A": "Z", "B": "A", "C": "M"}).to_numpy()

        original = fit_plsda(
            features, labels, n_components=2, n_permutations=12, random_state=5,
        )
        relabelled = fit_plsda(
            features, renamed, n_components=2, n_permutations=12, random_state=5,
        )

        np.testing.assert_allclose(original.vip, relabelled.vip, rtol=1e-4, atol=1e-6)
        self.assertEqual(original.status, "validated")
        self.assertTrue(np.isfinite(original.q2))
        self.assertTrue(np.isfinite(original.cv_accuracy))
        self.assertGreaterEqual(original.permutation_pvalue, 0.0)
        self.assertLessEqual(original.permutation_pvalue, 1.0)


class HeatmapTests(unittest.TestCase):
    def test_default_heatmap_selection_is_unsupervised_variance(self) -> None:
        frame = pd.DataFrame(
            {
                "Index": ["stable", "variable", "medium"],
                "Compounds": ["stable", "variable", "medium"],
                "S1": [1.0, -5.0, 0.0],
                "S2": [1.1, 5.0, 1.0],
                "S3": [0.9, -4.0, 0.0],
                "S4": [1.0, 4.0, 1.0],
            }
        )
        first_groups = {"S1": "A", "S2": "A", "S3": "B", "S4": "B"}
        permuted_groups = {"S1": "B", "S2": "A", "S3": "A", "S4": "B"}

        first = select_heatmap_features(
            frame, ["S1", "S2", "S3", "S4"], top_n=1,
            method="variance", groups=first_groups,
        )
        second = select_heatmap_features(
            frame, ["S1", "S2", "S3", "S4"], top_n=1,
            method="variance", groups=permuted_groups,
        )

        self.assertEqual(first["Index"].tolist(), ["variable"])
        self.assertEqual(second["Index"].tolist(), ["variable"])


class PipelineIntegrationTests(unittest.TestCase):
    def test_pipeline_connects_qc_normalization_vip_fdr_and_ora(self) -> None:
        feature_count = 30
        rows = []
        for index in range(feature_count):
            abundance = 1000.0 + index * 20.0
            treatment_multiplier = 4.0 if index < 5 else 1.0
            row = {
                "Index": f"m{index}",
                "Compounds": f"metabolite-{index}",
                "Formula": "C1H2",
                "kegg_map": (
                    "ko00010,ko01100" if index < 5
                    else "ko00020,ko01100" if index < 15
                    else "ko00030,ko01100"
                ),
            }
            for replicate, multiplier in enumerate([0.98, 1.00, 1.02], start=1):
                jitter_a = 1.0 + (((index + replicate) % 5) - 2) * 0.002
                jitter_b = 1.0 + (((index + 2 * replicate) % 7) - 3) * 0.002
                row[f"A{replicate}_liver"] = abundance * treatment_multiplier * multiplier * jitter_a
                row[f"B{replicate}_liver"] = abundance * multiplier * jitter_b
            qc_abundance = abundance * ((treatment_multiplier + 1.0) / 2.0)
            row.update({"QC01": qc_abundance * 0.99, "QC02": qc_abundance, "QC03": qc_abundance * 1.01})
            rows.append(row)
        frame = pd.DataFrame(rows)

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            input_path = root / "metabolites.xlsx"
            frame.to_excel(input_path, index=False)
            groups_path = root / "groups.yaml"
            groups_path.write_text(
                "conditions:\n  treatment: {display: A}\n  control: {display: B, reference: true}\n"
                "contrast_groups:\n  A vs B: [treatment, control]\n"
                "mouse_groups:\n  A1: treatment\n  A2: treatment\n  A3: treatment\n"
                "  B1: control\n  B2: control\n  B3: control\n",
                encoding="utf-8",
            )
            cfg = MetabolomicsConfig(
                experiment_id="integration-test",
                config_path=root / "config.yaml",
                output_dir=root / "results",
                tissue_data_paths={"liver": input_path},
                groups_path=groups_path,
                conditions={
                    "treatment": {"display": "A"},
                    "control": {"display": "B", "reference": True},
                },
                mouse_groups={
                    "A1": "treatment", "A2": "treatment", "A3": "treatment",
                    "B1": "control", "B2": "control", "B3": "control",
                },
                plsda_permutation_num=4,
                gsea_permutation_num=20,
            )

            run_pipeline(cfg, steps=[1, 2, 3, 5, 6])

            output = cfg.output_dir
            differential = pd.read_csv(
                output / "03_differential_analysis" / "liver" / "differential_analysis.csv"
            )
            dm = pd.read_csv(output / "00_dm_results" / "liver" / "dm_liver_A_vs_B.csv")
            pathway = pd.read_csv(output / "05_pathway_enrichment" / "liver" / "Pathway_A_vs_B.csv")

            self.assertTrue((output / "01_data_preprocessing" / "liver" / "normalization_factors.csv").exists())
            self.assertTrue((output / "01_data_preprocessing" / "liver" / "qc_feature_metrics.csv").exists())
            self.assertTrue((output / "02_qc_visualization" / "liver" / "PLS-DA_A_vs_B_validation.csv").exists())
            self.assertIn("VIP_A_vs_B", differential.columns)
            self.assertIn("Significant_A_vs_B", differential.columns)
            self.assertIn("significant_univariate", dm.columns)
            self.assertEqual(dm["significant"].tolist(), differential["Significant_A_vs_B"].tolist())
            self.assertIn("FDR", pathway.columns)
            self.assertNotIn("ko01100", pathway["ID"].tolist())
            self.assertTrue(pathway.set_index("ID").loc["ko00010", "FDR"] < 0.05)
            self.assertTrue(
                (output / "06_gsea" / "liver" / "gsea_A_vs_B_results.csv").exists()
            )


class ReportTests(unittest.TestCase):
    def test_report_separates_fdr_discoveries_from_nominal_findings(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            diff_dir = root / "03_differential_analysis" / "liver"
            diff_dir.mkdir(parents=True)
            pd.DataFrame(
                {
                    "Index": ["nominal-only", "fdr-hit"],
                    "Compounds": ["Nominal", "FDR hit"],
                    "Pval_A_vs_B": [0.01, 0.001],
                    "FDR_A_vs_B": [0.20, 0.01],
                    "Log2FC_A_vs_B": [1.0, 1.2],
                    "Significant_A_vs_B": [False, True],
                }
            ).to_csv(diff_dir / "differential_analysis.csv", index=False)
            cfg = MetabolomicsConfig(
                experiment_id="report-test",
                config_path=root / "config.yaml",
                output_dir=root,
                tissue_data_paths={"liver": root / "input.xlsx"},
                groups_path=root / "groups.yaml",
                conditions={"a": {"display": "A"}, "b": {"display": "B"}},
                mouse_groups={"1": "a", "2": "b"},
            )

            summary = generate_tissue_summary(
                "liver", "Liver", cfg, str(root), [("a", "b", "A_vs_B")],
            )

        self.assertIn("FDR 显著代谢物: **1**", summary)
        self.assertIn("名义发现 (原始 P): **2**", summary)


class DownstreamCompatibilityTests(unittest.TestCase):
    def test_multiomics_loader_names_new_evidence_columns_per_comparison(self) -> None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from omicsweave.integration.loaders.metabolomics_loader import MetabolomicsLoader

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tissue_dir = root / "00_dm_results" / "liver"
            tissue_dir.mkdir(parents=True)
            for comparison in ("A_vs_B", "C_vs_B"):
                pd.DataFrame(
                    {
                        "Index": ["m1"],
                        "Compounds": ["metabolite"],
                        "Formula": ["C1H2"],
                        "kegg_map": ["ko00010"],
                        f"VIP_{comparison}": [1.2],
                        "log2FoldChange": [1.0],
                        "statistic": [3.0],
                        "pvalue": [0.01],
                        "padj": [0.02],
                        "nominal_significant": [True],
                        "significant_univariate": [True],
                        "vip_pass": [True],
                        "significant_with_vip": [True],
                        "significant": [True],
                        "evidence_level": ["inferential"],
                    }
                ).to_csv(
                    tissue_dir / f"dm_liver_{comparison}.csv", index=False,
                )
            config = SimpleNamespace(
                metabolomics_results=root,
                tissues=["liver"],
                comparisons=["A_vs_B", "C_vs_B"],
                metabolomics_comparison_map=None,
                integration={},
            )

            loaded = MetabolomicsLoader(config).load()["diff"]["liver"]

        self.assertIn("SignificantUnivariate_A_vs_B", loaded.columns)
        self.assertIn("SignificantUnivariate_C_vs_B", loaded.columns)
        self.assertFalse(any(column.endswith(("_x", "_y")) for column in loaded.columns))


if __name__ == "__main__":
    unittest.main()
