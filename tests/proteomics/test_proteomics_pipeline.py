"""蛋白质组学流程的聚焦回归测试。"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd


from omicsweave.common.project import project_dir
from omicsweave.proteomics.diann_search import DiannSearchError, run_tissue_search  # noqa: E402
from omicsweave.proteomics.enrichment.gsea import (  # noqa: E402
    build_ranked_list,
    find_ssgsea_report,
    run_prerank_gsea,
    run_ssgsea,
)
from omicsweave.proteomics.experiment_config import (  # noqa: E402
    build_sample_metadata,
    get_sample_to_group_map,
    normalize_fdr_method,
)
from omicsweave.proteomics.plotting.volcano import classify_significance, fdr_cut_on_raw_p  # noqa: E402
from omicsweave.proteomics.proteomics_qc import compute_sample_qc_metrics, generate_qc_outputs  # noqa: E402
from omicsweave.proteomics.run_proteomics_pipeline import (  # noqa: E402
    _estimate_limma_priors,
    impute_data,
    load_and_clean_data,
    main,
    preprocess_data,
    run_differential_analysis,
)


class ConfigAndPreprocessingTests(unittest.TestCase):
    def test_fdr_alias_is_normalized_and_unknown_method_fails(self) -> None:
        self.assertEqual(normalize_fdr_method("bh"), "fdr_bh")
        self.assertEqual(normalize_fdr_method("fdr_bh"), "fdr_bh")
        with self.assertRaisesRegex(ValueError, "Unsupported FDR method"):
            normalize_fdr_method("not-a-method")

    def test_alphanumeric_mouse_ids_are_mapped_by_longest_suffix(self) -> None:
        mouse_groups = {"35": "wrong", "A35": "treated", "52": "control"}
        mapping = get_sample_to_group_map(
            ["LiverA35", "Heart52"], mouse_groups,
        )
        self.assertEqual(mapping, {"LiverA35": "treated", "Heart52": "control"})

        metadata = build_sample_metadata(
            ["LiverA35", "Heart52"],
            mouse_groups,
            {"A35": {"batch": "b1"}, "52": {"batch": "b2"}},
            tissue="liver",
        )
        self.assertEqual(metadata.loc["LiverA35", "mouse_id"], "A35")
        self.assertEqual(metadata.loc["Heart52", "batch"], "b2")

    def test_log_space_protein_aggregation_uses_geometric_mean(self) -> None:
        table = pd.DataFrame({
            "Run": ["run_Liver_52", "run_Liver_52"],
            "Protein.Group": ["P1", "P2"],
            "Genes": ["GeneA", "GeneA"],
            "PG.MaxLFQ": [4.0, 16.0],
            "PG.Q.Value": [0.001, 0.001],
            "Global.PG.Q.Value": [0.001, 0.001],
        })
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "report_test.tsv"
            table.to_csv(path, sep="\t", index=False)
            matrix, mapping = load_and_clean_data(
                str(path), protein_aggregation="log2_mean",
            )

        self.assertIsNotNone(matrix)
        assert matrix is not None
        self.assertAlmostEqual(matrix.loc["GeneA", "Liver52"], 8.0)
        self.assertEqual(mapping["GeneA"], "P1")

    def test_preprocessing_and_imputation_are_reproducible(self) -> None:
        raw = pd.DataFrame(
            {
                "s1": [4.0, 8.0, 2.0],
                "s2": [8.0, np.nan, np.nan],
                "s3": [16.0, 32.0, np.nan],
            },
            index=["p1", "p2", "too_sparse"],
        )
        normalized = preprocess_data(raw, min_total_observations=2)
        self.assertNotIn("too_sparse", normalized.index)

        first = impute_data(normalized, seed=17)
        second = impute_data(normalized, seed=17)
        third = impute_data(normalized, seed=18)
        pd.testing.assert_frame_equal(first, second)
        self.assertFalse(np.allclose(first.to_numpy(), third.to_numpy()))
        self.assertFalse(first.isna().any().any())


class DifferentialAnalysisTests(unittest.TestCase):
    @staticmethod
    def _balanced_expression(n_proteins: int = 60) -> tuple[pd.DataFrame, dict[str, str]]:
        rng = np.random.default_rng(11)
        samples = ["A1", "A2", "A3", "B1", "B2", "B3"]
        matrix = rng.normal(0, 0.35, size=(n_proteins, len(samples)))
        matrix[:8, :3] += 2.0
        df = pd.DataFrame(
            matrix,
            index=[f"P{i:03d}" for i in range(n_proteins)],
            columns=samples,
        )
        groups = {s: ("A" if s.startswith("A") else "B") for s in samples}
        return df, groups

    def test_limma_uses_real_bh_and_writes_inference_provenance(self) -> None:
        expression, groups = self._balanced_expression()
        with tempfile.TemporaryDirectory() as tmp:
            result = run_differential_analysis(
                expression,
                groups,
                [("A", "B", "A_vs_B")],
                tmp,
                "liver",
                pval_cutoff=0.05,
                log2fc_cutoff=0.58,
                fdr_method="fdr_bh",
                dp_method="limma",
                min_observations_per_group=2,
            )
            per_comparison = pd.read_csv(Path(tmp) / "dp_liver_A_vs_B.csv")
            design = pd.read_csv(Path(tmp) / "differential_design_report.csv")

        self.assertIn("Tstat_A_vs_B", result.columns)
        self.assertIn("Evidence_A_vs_B", result.columns)
        self.assertTrue((per_comparison["evidence_level"] == "inferential").all())
        self.assertFalse(
            np.allclose(
                per_comparison["pvalue"], per_comparison["padj"], equal_nan=True,
            )
        )
        self.assertEqual(design.loc[0, "route"], "inferential")
        self.assertGreater(int(per_comparison["significant"].sum()), 0)

    def test_n_less_than_two_is_descriptive_only(self) -> None:
        expression = pd.DataFrame(
            {"A1": [3.0, 2.0], "B1": [1.0, 1.2]},
            index=["P1", "P2"],
        )
        groups = {"A1": "A", "B1": "B"}
        with tempfile.TemporaryDirectory() as tmp:
            run_differential_analysis(
                expression,
                groups,
                [("A", "B", "A_vs_B")],
                tmp,
                "heart",
                pval_cutoff=0.05,
                log2fc_cutoff=0.58,
                fdr_method="fdr_bh",
                dp_method="limma",
                min_observations_per_group=2,
            )
            result = pd.read_csv(Path(tmp) / "dp_heart_A_vs_B.csv")
            design = pd.read_csv(Path(tmp) / "differential_design_report.csv")

        self.assertTrue(result["pvalue"].isna().all())
        self.assertTrue(result["padj"].isna().all())
        self.assertTrue((result["evidence_level"] == "insufficient").all())
        self.assertFalse(result["significant"].any())
        self.assertEqual(design.loc[0, "route"], "insufficient_data")

    def test_tiny_feature_set_falls_back_to_ordinary_ols_t(self) -> None:
        expression = pd.DataFrame(
            {
                "A1": [3.0, 2.2], "A2": [3.2, 2.0], "A3": [2.8, 2.1],
                "B1": [1.0, 2.0], "B2": [1.2, 2.2], "B3": [0.8, 1.9],
            },
            index=["P1", "P2"],
        )
        groups = {sample: sample[0] for sample in expression.columns}
        with tempfile.TemporaryDirectory() as tmp:
            run_differential_analysis(
                expression, groups, [("A", "B", "A_vs_B")], tmp, "liver",
                0.05, 0.58, "fdr_bh", "limma",
            )
            result = pd.read_csv(Path(tmp) / "dp_liver_A_vs_B.csv")

        self.assertTrue(result["pvalue"].notna().all())
        self.assertTrue((result["inference_status"] == "ok_unmoderated_t").all())
        self.assertTrue((result["evidence_level"] == "inferential").all())

    def test_balanced_batch_is_adjusted_but_confounding_is_rejected(self) -> None:
        rng = np.random.default_rng(3)
        samples = [f"S{i}" for i in range(8)]
        expression = pd.DataFrame(
            rng.normal(size=(40, 8)),
            index=[f"P{i}" for i in range(40)],
            columns=samples,
        )
        groups = {s: ("A" if i < 4 else "B") for i, s in enumerate(samples)}
        balanced = pd.DataFrame(
            {
                "condition": [groups[s] for s in samples],
                "batch": ["b1", "b2", "b1", "b2", "b1", "b2", "b1", "b2"],
            },
            index=samples,
        )
        confounded = balanced.copy()
        confounded["batch"] = ["b1"] * 4 + ["b2"] * 4

        with tempfile.TemporaryDirectory() as tmp:
            run_differential_analysis(
                expression, groups, [("A", "B", "A_vs_B")], tmp, "liver",
                0.05, 0.58, "fdr_bh", "limma",
                sample_metadata=balanced,
                design_covariates=["batch"],
                min_observations_per_group=2,
            )
            balanced_report = pd.read_csv(Path(tmp) / "differential_design_report.csv")
        self.assertEqual(balanced_report.loc[0, "route"], "inferential")
        self.assertIn("batch", balanced_report.loc[0, "design_formula"])

        with tempfile.TemporaryDirectory() as tmp:
            run_differential_analysis(
                expression, groups, [("A", "B", "A_vs_B")], tmp, "liver",
                0.05, 0.58, "fdr_bh", "limma",
                sample_metadata=confounded,
                design_covariates=["batch"],
                min_observations_per_group=2,
            )
            rejected = pd.read_csv(Path(tmp) / "dp_liver_A_vs_B.csv")
            rejected_report = pd.read_csv(Path(tmp) / "differential_design_report.csv")
        self.assertEqual(rejected_report.loc[0, "route"], "design_not_identifiable")
        self.assertTrue(rejected["pvalue"].isna().all())

    def test_limma_prior_estimator_is_scale_equivariant(self) -> None:
        variances = np.array([0.2, 0.3, 0.5, 0.7, 1.1, 1.4])
        dfs = np.array([4.0, 4.0, 4.0, 4.0, 4.0, 4.0])
        scale, prior_df = _estimate_limma_priors(variances, dfs)
        scale10, prior_df10 = _estimate_limma_priors(variances * 10.0, dfs)
        self.assertGreater(scale, 0)
        self.assertGreater(prior_df, 0)
        self.assertAlmostEqual(scale10 / scale, 10.0, places=6)
        self.assertAlmostEqual(prior_df10, prior_df, places=6)


class PlotEnrichmentAndQCTests(unittest.TestCase):
    def test_volcano_classification_uses_adjusted_pvalue(self) -> None:
        frame = pd.DataFrame({
            "log2FoldChange": [1.0, -1.0, 1.0],
            "pvalue": [1e-6, 1e-6, 0.2],
            "padj": [0.2, 0.01, 0.3],
        })
        labels = classify_significance(
            frame["log2FoldChange"], frame["padj"], 0.58, 0.05,
        )
        self.assertEqual(labels.tolist(), ["NS", "Down", "NS"])

    def test_fdr_line_sits_at_largest_passing_raw_p(self) -> None:
        # BH is monotone in raw p, so the boundary on a -log10(raw p) axis is
        # the largest raw p that still clears the adjusted cutoff.
        frame = pd.DataFrame({
            "pvalue": [1e-6, 1e-4, 2e-3, 0.04, 0.30],
            "padj": [0.001, 0.008, 0.049, 0.20, 0.80],
        })
        self.assertAlmostEqual(
            fdr_cut_on_raw_p(frame["pvalue"], frame["padj"], 0.05), 2e-3,
        )

    def test_fdr_line_is_none_when_nothing_passes(self) -> None:
        # Panels where the threshold is unreachable must not draw a line that
        # no point can reach; the caller labels them explicitly instead.
        frame = pd.DataFrame({
            "pvalue": [0.01, 0.20, np.nan],
            "padj": [0.90, 0.99, np.nan],
        })
        self.assertIsNone(fdr_cut_on_raw_p(frame["pvalue"], frame["padj"], 0.05))

    def test_fdr_line_ignores_features_without_a_raw_p(self) -> None:
        frame = pd.DataFrame({
            "pvalue": [1e-4, np.nan],
            "padj": [0.01, 0.02],
        })
        self.assertAlmostEqual(
            fdr_cut_on_raw_p(frame["pvalue"], frame["padj"], 0.05), 1e-4,
        )

    def test_ssgsea_detects_real_gseapy_report_name(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "ssgsea_kegg"
            out.mkdir()
            expected = out / "gseapy.gene_set.ssgsea.report.csv"
            expected.write_text("Name,Term,NES\nS1,Path,1\n", encoding="utf-8")
            self.assertEqual(find_ssgsea_report(out), expected)

            with patch("gseapy.ssgsea") as mocked:
                status = run_ssgsea(
                    pd.DataFrame({"S1": [1.0]}, index=["GeneA"]),
                    tmp,
                    gene_sets=["kegg"],
                )
                mocked.assert_not_called()
            self.assertEqual(status.loc[0, "status"], "skipped_existing")

    def test_prerank_uses_all_finite_ranked_proteins(self) -> None:
        results = pd.DataFrame({
            "Protein": [f"G{i}" for i in range(20)],
            "Log2FC_A_vs_B": np.linspace(-2, 2, 20),
            "Tstat_A_vs_B": np.linspace(-4, 4, 20),
            "Pval_A_vs_B": np.linspace(0.001, 0.5, 20),
        })
        ranked = build_ranked_list(results, "A_vs_B", metric="moderated_t")
        self.assertEqual(len(ranked), 20)
        self.assertTrue(ranked.iloc[0, 1] > ranked.iloc[-1, 1])

        def fake_prerank(*, rnk, outdir, **_kwargs):
            self.assertEqual(len(rnk), 20)
            target = Path(outdir)
            target.mkdir(parents=True, exist_ok=True)
            pd.DataFrame({"Term": ["Path"], "NES": [1.5], "FDR q-val": [0.01]}).to_csv(
                target / "gseapy.gene_set.prerank.report.csv", index=False,
            )
            return object()

        with tempfile.TemporaryDirectory() as tmp, patch(
            "gseapy.prerank", side_effect=fake_prerank,
        ):
            status = run_prerank_gsea(
                results,
                tmp,
                [("A", "B", "A_vs_B")],
                gene_sets=["kegg"],
                min_size=3,
                permutation_num=100,
                seed=7,
            )
        self.assertEqual(status.loc[0, "status"], "completed")

    def test_qc_metrics_and_outputs_cover_core_artifacts(self) -> None:
        raw = pd.DataFrame(
            {
                "A1": [1.0, 2.0, np.nan, 4.0],
                "A2": [1.1, 2.1, 3.0, 4.1],
                "B1": [2.0, 3.0, 4.0, np.nan],
                "B2": [2.1, 3.1, 4.1, 5.1],
            },
            index=["P1", "P2", "P3", "P4"],
        )
        metadata = pd.DataFrame(
            {"condition": ["A", "A", "B", "B"], "batch": ["x", "y", "x", "y"]},
            index=raw.columns,
        )
        metrics = compute_sample_qc_metrics(raw)
        self.assertEqual(int(metrics.loc["A1", "n_proteins"]), 3)
        self.assertAlmostEqual(metrics.loc["A1", "missing_pct"], 25.0)

        normalized = np.log2(raw).sub(np.log2(raw).median(axis=0), axis=1)
        with tempfile.TemporaryDirectory() as tmp:
            generate_qc_outputs(raw, normalized, metadata, tmp, seed=9)
            expected = {
                "sample_qc_metrics.csv",
                "protein_missingness.csv",
                "replicate_correlations.csv",
                "intensity_distributions.png",
                "sample_detection_missingness.png",
                "correlation_heatmap.png",
                "pc_covariate_association.csv",
                "qc_summary.json",
            }
            self.assertTrue(expected.issubset({p.name for p in Path(tmp).iterdir()}))


class DiannSearchTests(unittest.TestCase):
    @unittest.skipUnless(
        (project_dir() / "experiments/20250425_cimr_wt_mice_dio/proteomics/config.yaml").exists(),
        "需要实验项目数据：设置 OMICSWEAVE_PROJECT_DIR 指向实验项目目录",
    )
    def test_existing_experiment_leaves_search_disabled(self) -> None:
        from omicsweave.proteomics.experiment_config import load_proteomics_config

        cfg = load_proteomics_config(experiment_id="20250425_cimr_wt_mice_dio")
        self.assertFalse(cfg.search_enabled)
        self.assertEqual(cfg.raw_dirs, {})
        self.assertIsNone(cfg.diann_exe)
        self.assertIsNone(cfg.search_fasta)

    def test_fake_exe_command_matches_library_free_defaults(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _load_search_config(Path(tmp))
            calls: list[list[str]] = []

            report = run_tissue_search(
                "brain", cfg, runner=_fake_diann(calls, version="1.9"),
            )

            self.assertEqual(len(calls), 1)
            command = calls[0]
            raw_name = "MSF250153_zzk_Astral_23min_DIA_100ng_Brain_52.raw"
            self.assertEqual(command[0], str(cfg.diann_exe))
            self.assertEqual(command[command.index("--f") + 1], str(cfg.raw_dirs["brain"] / raw_name))
            self.assertEqual(command[command.index("--fasta") + 1], str(cfg.search_fasta))
            self.assertEqual(command[command.index("--cut") + 1], "K*,R*")
            self.assertEqual(command[command.index("--missed-cleavages") + 1], "1")
            self.assertEqual(command[command.index("--min-pep-len") + 1], "7")
            self.assertEqual(command[command.index("--max-pep-len") + 1], "30")
            self.assertEqual(command[command.index("--min-pr-charge") + 1], "1")
            self.assertEqual(command[command.index("--max-pr-charge") + 1], "4")
            self.assertEqual(command[command.index("--var-mods") + 1], "1")
            self.assertEqual(command[command.index("--qvalue") + 1], "0.01")
            self.assertEqual(command[command.index("--threads") + 1], "4")
            self.assertEqual(command[command.index("--verbose") + 1], "1")
            self.assertEqual(command[command.index("--out") + 1], str(report))
            modifications = [
                command[index + 1]
                for index, token in enumerate(command)
                if token == "--var-mod"
            ]
            self.assertEqual(modifications, [
                "UniMod:35,15.994915,M",
                "UniMod:1,42.010565,*n",
            ])
            for flag in (
                "--fasta-search",
                "--predictor",
                "--gen-spec-lib",
                "--met-excision",
                "--unimod4",
                "--matrices",
                "--reanalyse",
                "--smart-profiling",
            ):
                self.assertIn(flag, command)
            self.assertNotIn("--relaxed-prot-inf", command)
            self.assertNotIn("--lib", command)
            self.assertTrue(report.is_file())
            provenance = json.loads(
                (report.parent / "search_provenance.json").read_text(encoding="utf-8"),
            )
            self.assertEqual(provenance["diann_version"], "1.9")
            self.assertEqual(provenance["samples"], [{
                "raw_file": raw_name,
                "sample": "Brain52",
                "mouse_id": "52",
            }])
            self.assertTrue((report.parent / "diann_cmd.txt").is_file())
            self.assertIn("DIA-NN 1.9", (report.parent / "diann.log").read_text(encoding="utf-8"))

    def test_unmapped_raw_name_does_not_call_the_executable(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _load_search_config(Path(tmp), raw_name="notes.raw")
            calls: list[list[str]] = []
            with self.assertRaisesRegex(ValueError, "notes.raw"):
                run_tissue_search(
                    "brain", cfg, runner=_fake_diann(calls, version="1.9"),
                )
            self.assertEqual(calls, [])
            self.assertFalse((cfg.output_dir / "00_diann_search" / "brain" / "search_provenance.json").exists())

    def test_unchanged_fingerprint_skips_and_force_runs_again(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _load_search_config(Path(tmp))
            calls: list[list[str]] = []
            runner = _fake_diann(calls, version="1.9.2")
            first = run_tissue_search("brain", cfg, runner=runner)
            second = run_tissue_search("brain", cfg, runner=runner)
            self.assertEqual(first, second)
            self.assertEqual(len(calls), 1)

            third = run_tissue_search("brain", cfg, force=True, runner=runner)
            self.assertEqual(third, first)
            self.assertEqual(len(calls), 2)

    def test_unexpected_version_is_not_accepted(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _load_search_config(Path(tmp))
            calls: list[list[str]] = []
            with self.assertRaisesRegex(DiannSearchError, "2.0"):
                run_tissue_search(
                    "brain", cfg, runner=_fake_diann(calls, version="2.0"),
                )
            self.assertEqual(len(calls), 1)
            provenance = cfg.output_dir / "00_diann_search" / "brain" / "search_provenance.json"
            self.assertFalse(provenance.exists())

    def test_force_alone_does_not_launch_search(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cfg_path = _write_search_config(Path(tmp))
            argv = ["pipeline", "--config", str(cfg_path), "--force", "--skip-gsea"]
            with patch("sys.argv", argv), \
                    patch("omicsweave.proteomics.run_proteomics_pipeline.run_tissue_search") as search, \
                    patch("omicsweave.proteomics.run_proteomics_pipeline.analyze_tissue") as analyze:
                main()
            search.assert_not_called()
            analyze.assert_called_once()
            _args, kwargs = analyze.call_args
            self.assertIsNone(kwargs["report_path"])
            self.assertTrue(kwargs["force"])

    def test_search_flag_feeds_the_new_report_and_search_only_stops(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cfg_path = _write_search_config(root)
            report = root / "report_brain.tsv"
            report.write_text("placeholder\n", encoding="utf-8")
            argv = ["pipeline", "--config", str(cfg_path), "--search", "--skip-gsea"]
            with patch("sys.argv", argv), \
                    patch("omicsweave.proteomics.run_proteomics_pipeline.run_tissue_search", return_value=report) as search, \
                    patch("omicsweave.proteomics.run_proteomics_pipeline.analyze_tissue") as analyze:
                main()
            search.assert_called_once()
            _args, search_kwargs = search.call_args
            self.assertFalse(search_kwargs["force"])
            _args, analyze_kwargs = analyze.call_args
            self.assertEqual(analyze_kwargs["report_path"], report)

            analyze.reset_mock()
            search.reset_mock()
            argv[-2] = "--search-only"
            with patch("sys.argv", argv), \
                    patch("omicsweave.proteomics.run_proteomics_pipeline.run_tissue_search", return_value=report) as search, \
                    patch("omicsweave.proteomics.run_proteomics_pipeline.analyze_tissue") as analyze:
                main()
            search.assert_called_once()
            analyze.assert_not_called()


def _write_search_config(root: Path, *, raw_name: str = "MSF250153_zzk_Astral_23min_DIA_100ng_Brain_52.raw") -> Path:
    raw_dir = root / "raw"
    raw_dir.mkdir()
    (raw_dir / raw_name).write_bytes(b"raw")
    exe = root / "fake_diann.exe"
    exe.write_bytes(b"")
    fasta = root / "mouse.fasta"
    fasta.write_text(">sp|P12345|TEST_MOUSE\nACDEFGHIK\n", encoding="utf-8")
    (root / "groups.yaml").write_text(
        "conditions:\n  hfd:\n    display: HFD\nmouse_groups:\n  '52': hfd\n",
        encoding="utf-8",
    )
    output_dir = root / "proteinomics_results"
    config = root / "config.yaml"
    config.write_text(
        "\n".join([
            "experiment_id: diann_search_test",
            "paths:",
            f"  output_dir: {json.dumps(str(output_dir))}",
            "  raw_dirs:",
            "    brain: raw",
            "params:",
            "  tissues: [brain]",
            "search:",
            "  enabled: true",
            "  expected_version: '1.9'",
            "  diann_exe: fake_diann.exe",
            "  fasta: mouse.fasta",
            "  threads: 4",
            "  qvalue: 0.01",
            "  reanalyse: true",
            "  smart_profiling: true",
            "groups_file: groups.yaml",
            "",
        ]),
        encoding="utf-8",
    )
    return config


def _load_search_config(root: Path, **kwargs: str):
    from omicsweave.proteomics.experiment_config import load_proteomics_config

    return load_proteomics_config(config_path=str(_write_search_config(root, **kwargs)))


def _fake_diann(calls: list[list[str]], *, version: str):
    import subprocess

    def runner(command, **_kwargs):
        calls.append(list(command))
        report = Path(command[command.index("--out") + 1])
        report.parent.mkdir(parents=True, exist_ok=True)
        report.write_text("Run\tProtein.Group\tPG.MaxLFQ\n", encoding="utf-8")
        return subprocess.CompletedProcess(
            command, 0, stdout=f"DIA-NN {version}\n", stderr="",
        )

    return runner


if __name__ == "__main__":
    unittest.main()
