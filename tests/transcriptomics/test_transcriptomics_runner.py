#!/usr/bin/env python3
"""转录组统一运行器：与旧流程逐值一致；kallisto 输入、样本量不足时的秩和检验、原始层输入的提示。

核心回归：同一个旧版转录组实验（两个组织、两个批次），旧流程从基因计数矩阵补跑差异表达
（rerun_de_analysis，PyDESeq2 ~batch + condition）与“omicsweave migrate → omicsweave run --omics transcriptomics”
的结果逐值一致；用同一个本地基因集库时，预排序 GSEA 也与旧流程（gsea/_run.py）逐值一致。
"""

from __future__ import annotations

import contextlib
import io
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd
import yaml

from omicsweave.contract.dataset import read_dataset
from omicsweave.contract.tables import read_differential, read_enrichment
from omicsweave.core.genesets import write_gmt
from omicsweave.core.migrate import migrate_experiment, write_migration
from omicsweave.core.project import project_context
from omicsweave.core.run import RunError, run_experiment
from tests.visualization.test_volcano import assert_volcanoes

try:
    import pydeseq2  # noqa: F401
except ImportError:  # pragma: no cover - 依赖可选
    pydeseq2 = None

EXP = "rna_exp"
ANIMALS = {str(i): ("LFD" if i <= 4 else "HFD") for i in range(1, 9)}
BATCHES = {a: ("b1" if int(a) % 2 else "b2") for a in ANIMALS}
TISSUES = ("liver", "heart")
COUNTS = f"analysis_results/{EXP}/transcriptomics_results/04_de_results/gene_counts_matrix.csv"


def gene_counts(columns: list[str], hfd: set[str], b2: set[str], seed: int = 0, n_genes: int = 200) -> pd.DataFrame:
    """负二项计数：前 20 个基因在 HFD 中升高 4 倍，第二批次整体测序深度高 30%。"""
    rng = np.random.default_rng(seed)
    base = rng.lognormal(5, 1.2, n_genes)
    data = {}
    for column in columns:
        mean = base * np.where(np.arange(n_genes) < 20, 4.0 if column in hfd else 1.0, 1.0)
        mean = mean * (1.3 if column in b2 else 1.0)
        data[column] = rng.negative_binomial(5, 5 / (5 + mean))
    return pd.DataFrame(data, index=pd.Index([f"Gene{i:04d}" for i in range(n_genes)], name="gene"))


def build_legacy_experiment(root: Path) -> None:
    tx = root / "experiments" / EXP / "transcriptomics_data"
    tx.mkdir(parents=True)
    samples = [f"{a}_{t}" for t in TISSUES for a in ANIMALS]
    (tx / "config.yaml").write_text(yaml.safe_dump({
        "experiment_id": EXP,
        "paths": {"raw_data_root": "raw/rnaseq",
                  "sample_list_tsv": f"experiments/{EXP}/transcriptomics_data/sample_list.tsv",
                  "output_dir": "analysis_results/${experiment_id}/transcriptomics_results"},
        "reference": {"kallisto_exe": "tools/kallisto.exe", "transcriptome_fasta": "references/tx.fa.gz",
                      "kallisto_index": "references/tx.idx"},
        "params": {"threads": 2, "tissues": list(TISSUES)},
    }), encoding="utf-8")
    (tx / "groups.yaml").write_text(yaml.safe_dump({
        "conditions": {"LFD": {"display": "LFD"}, "HFD": {"display": "HFD"}},
        "mouse_groups": {int(a): g for a, g in ANIMALS.items()},
        "sample_batches": {int(a): b for a, b in BATCHES.items()},
    }), encoding="utf-8")
    (tx / "sample_list.tsv").write_text("sample\tfile1\tfile2\n" + "".join(
        f"{s}\t{s}_R1.fq.gz\t{s}_R2.fq.gz\n" for s in samples), encoding="utf-8")
    raw = root / "raw" / "rnaseq"
    raw.mkdir(parents=True)
    for s in samples:
        for end in ("R1", "R2"):
            (raw / f"{s}_{end}.fq.gz").write_bytes(b"")
    hfd = {s for s in samples if ANIMALS[s.split("_")[0]] == "HFD"}
    b2 = {s for s in samples if BATCHES[s.split("_")[0]] == "b2"}
    path = root / COUNTS
    path.parent.mkdir(parents=True)
    gene_counts(samples, hfd, b2).to_csv(path)


def _quiet(func, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        return func(*args, **kwargs)


def _set_input(root: Path, value: dict, **analysis) -> None:
    path = root / "experiments" / EXP / "experiment.yaml"
    doc = yaml.safe_load(path.read_text(encoding="utf-8"))
    doc["omics"]["transcriptomics"]["input"] = value
    doc["omics"]["transcriptomics"]["analysis"].update(analysis)
    path.write_text(yaml.safe_dump(doc, allow_unicode=True, sort_keys=False), encoding="utf-8")


def write_gene_sets(path: Path) -> Path:
    """DE_SET 为上调的 20 个基因；其余为无差异的基因块。"""
    genes = [f"Gene{i:04d}" for i in range(200)]
    return write_gmt({"DE_SET": genes[:20], "MIXED": genes[10:40], "BLOCK_A": genes[50:80],
                      "BLOCK_B": genes[80:110], "BLOCK_C": genes[110:140]}, path)


@unittest.skipIf(pydeseq2 is None, "pydeseq2 not installed")
class MigratedRunMatchesLegacyTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        from omicsweave.transcriptomics import rerun_de_analysis
        from omicsweave.transcriptomics.gsea import _run as legacy_gsea

        cls._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        cls.root = root = Path(cls._tmp.name).resolve()
        build_legacy_experiment(root)
        cls.gmt = write_gene_sets(root / "gene_sets.gmt")
        # 旧流程：从基因计数矩阵补跑差异表达（与主流程共用 run_pydeseq2_analysis）
        with project_context(root), mock.patch.dict(os.environ), \
                mock.patch.object(sys, "argv", ["rerun_de_analysis", "--experiment", EXP]):
            _quiet(rerun_de_analysis.main)
            # 旧流程的 GSEA：库名换成同一个本地 GMT，其余参数不变
            with mock.patch.object(legacy_gsea, "GENE_SET_LIBS", str(cls.gmt)):
                _quiet(legacy_gsea.apply_config, experiment=EXP)
                _quiet(legacy_gsea._generate_gsea_plots)
            cls.legacy_gsea_dir = Path(legacy_gsea.PLOT_DIR) / "5_gsea"

        write_migration(_quiet(migrate_experiment, EXP, root), root)
        # 迁移结果仍是原始层输入（FASTQ），统一运行尚不支持，应给出改用定量结果的提示
        cls.raw_error = ""
        try:
            run_experiment(EXP, ["transcriptomics"], root=root)
        except RunError as exc:
            cls.raw_error = str(exc)
        _set_input(root, {"level": "table", "format": "counts_matrix", "files": [COUNTS]}, gene_sets=[str(cls.gmt)])
        cls.run_dir, cls.notes = _quiet(run_experiment, EXP, ["transcriptomics"], root=root)
        tables = cls.run_dir / "omics" / "transcriptomics" / "tables"
        cls.table = read_differential(tables / "differential.parquet")
        cls.enrichment = read_enrichment(tables / "enrichment.parquet")

    @classmethod
    def tearDownClass(cls) -> None:
        cls._tmp.cleanup()

    def test_raw_input_points_to_table_formats(self) -> None:
        self.assertIn("阶段 3", self.raw_error)
        self.assertIn("counts_matrix", self.raw_error)

    def test_differential_matches_legacy_pipeline(self) -> None:
        self.assertEqual(set(self.table["contrast_id"]), {"HFD_vs_LFD"})
        for tissue in TISSUES:
            legacy = pd.read_csv(self.root / Path(COUNTS).parent / tissue / "de_HFD_vs_LFD.csv").set_index("gene")
            new = self.table[self.table["tissue"] == tissue].set_index("feature_id")
            self.assertEqual(sorted(new.index), sorted(legacy.index))
            new = new.loc[legacy.index]
            for legacy_column, column in (("log2FoldChange", "log2fc"), ("stat", "statistic"), ("pvalue", "pvalue"),
                                          ("padj", "padj"), ("lfcSE", "lfc_se"),
                                          ("log2FoldChange_unshrunk", "log2fc_unshrunk")):
                np.testing.assert_allclose(new[column].to_numpy(dtype=float),
                                           legacy[legacy_column].to_numpy(dtype=float),
                                           rtol=1e-10, equal_nan=True, err_msg=f"{tissue} {column}")
            self.assertEqual(list(new["selected"]), list(legacy["significant"].astype(bool)))
            self.assertEqual(set(new["design_formula"]), {"~batch + condition"})
            self.assertEqual(set(new["method"]), {"deseq2"})
            self.assertTrue(new.loc["Gene0000", "selected"])
            self.assertTrue((new["n_numerator"] == 4).all() and (new["n_denominator"] == 4).all())

    def test_gsea_matches_legacy_and_ora_finds_the_up_set(self) -> None:
        gsea = self.enrichment[self.enrichment["method"] == "gsea_prerank"]
        for tissue in TISSUES:
            legacy = pd.read_csv(self.legacy_gsea_dir / f"gsea_{tissue}_HFD_vs_LFD_results.csv").set_index("Term")
            new = gsea[gsea["tissue"] == tissue].set_index("term_name").loc[legacy.index]
            # 旧流程保留 gseapy 报告的 0；统一结果表改为 1/(置换次数+1)（默认 1000 次）并在 flags 中注明
            self.assertEqual(list(new["flags"].str.contains("pvalue_floored")), list(legacy["NOM p-val"] == 0))
            legacy[["NOM p-val", "FDR q-val"]] = legacy[["NOM p-val", "FDR q-val"]].replace(0.0, 1 / 1001)
            for legacy_column, column in (("NES", "nes"), ("ES", "es"), ("NOM p-val", "pvalue"),
                                          ("FDR q-val", "padj")):
                np.testing.assert_allclose(new[column].to_numpy(dtype=float),
                                           legacy[legacy_column].to_numpy(dtype=float), rtol=1e-12,
                                           err_msg=f"{tissue} {column}")
            self.assertGreater(new.loc["DE_SET", "nes"], 0)
            self.assertEqual(set(new["gene_set_library"]), {"gene_sets"})
        ora = self.enrichment[(self.enrichment["method"] == "ora") & (self.enrichment["query"] == "up")]
        de_set = ora[ora["term_name"] == "DE_SET"]
        self.assertEqual(sorted(de_set["tissue"]), sorted(TISSUES))
        self.assertTrue((de_set["padj"] < 0.05).all())
        # 背景为有 p 值的基因，而不是基因集库中的全部基因
        self.assertTrue((de_set["universe_size"] <= 200).all())

    def test_bundle_contents(self) -> None:
        tx = self.run_dir / "omics" / "transcriptomics"
        liver = read_dataset(tx / "data" / "liver")
        self.assertEqual(sorted(liver.matrix.columns), sorted(f"{a}_liver" for a in ANIMALS))
        self.assertEqual(liver.meta.value_type.value, "counts")
        sizes = pd.read_parquet(tx / "qc" / "library_sizes.parquet")
        self.assertEqual(len(sizes), 16)
        self.assertEqual(self.notes, [])
        assert_volcanoes(self, tx)


class KallistoInputTest(unittest.TestCase):
    def test_kallisto_quant_with_low_replicates_falls_back_to_rank_sum(self) -> None:
        from omicsweave.transcriptomics.normalization import AGGREGATION_LENGTH_SCALED_TPM
        from omicsweave.transcriptomics.rnaseq_pipeline import aggregate_counts_to_gene, merge_kallisto_counts

        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp).resolve()
            exp_dir = root / "experiments" / EXP
            quant = exp_dir / "quant"
            rng = np.random.default_rng(1)
            transcripts = [f"TX{i:03d}" for i in range(30)]
            tx2gene = pd.DataFrame({"transcript_id": transcripts,
                                    "gene_id": [f"ENSG{i // 2:03d}" for i in range(30)],
                                    "gene_symbol": [f"Gene{i // 2:03d}" for i in range(30)]})
            for name in ("Q1", "Q2", "Q3"):
                (quant / name).mkdir(parents=True)
                length = rng.integers(500, 3000, 30)
                counts = rng.negative_binomial(5, 0.05, 30).astype(float)
                tpm = counts / length
                pd.DataFrame({"target_id": transcripts, "length": length, "eff_length": length - 150.0,
                              "est_counts": counts, "tpm": tpm / tpm.sum() * 1e6}).to_csv(
                    quant / name / "abundance.tsv", sep="\t", index=False)
                (quant / name / "run_info.json").write_text(json.dumps(
                    {"n_processed": 1000, "n_pseudoaligned": 800, "p_pseudoaligned": 80.0}), encoding="utf-8")
            tx2gene.to_csv(exp_dir / "tx2gene.tsv", sep="\t", index=False)
            (root / "project.yaml").write_text("schema_version: 1\nproject: {id: p, species: mus_musculus}\n",
                                               encoding="utf-8")
            (exp_dir / "experiment.yaml").write_text(yaml.safe_dump({
                "schema_version": 1, "experiment": {"id": EXP},
                "design": {"factors": {"diet": {"level": "subject", "levels": ["LFD", "HFD"]}},
                           "contrasts": [{"id": "HFD_vs_LFD", "factor": "diet", "numerator": "HFD",
                                          "denominator": "LFD"}]},
                "omics": {"transcriptomics": {"enabled": True, "input": {
                    "level": "table", "format": "kallisto", "files": [f"experiments/{EXP}/quant"],
                    "tx2gene": f"experiments/{EXP}/tx2gene.tsv"}}},
            }, sort_keys=False), encoding="utf-8")
            pd.DataFrame({"subject_id": ["m1", "m2", "m3"], "diet": ["LFD", "HFD", "HFD"]}).to_csv(
                exp_dir / "subjects.tsv", sep="\t", index=False)
            pd.DataFrame({"sample_id": ["liver_m1", "liver_m2", "liver_m3"], "subject_id": ["m1", "m2", "m3"],
                          "omics": "transcriptomics", "tissue": "liver", "data_column": ["Q1", "Q2", "Q3"]}).to_csv(
                exp_dir / "samples.tsv", sep="\t", index=False)

            run_dir, _ = _quiet(run_experiment, EXP, ["transcriptomics"], steps=["qc", "differential"], root=root)
            tx = run_dir / "omics" / "transcriptomics"
            table = read_differential(tx / "tables" / "differential.parquet")
            dataset = read_dataset(tx / "data" / "liver")
            quant_summary = pd.read_parquet(tx / "qc" / "quant_summary.parquet")

            mapping = {t: {"gene_id": g, "gene_symbol": s} for t, g, s in tx2gene.itertuples(index=False)}
            tx_counts, tpm, efflen = _quiet(merge_kallisto_counts, str(quant), ["Q1", "Q2", "Q3"])
            expected, _ = _quiet(aggregate_counts_to_gene, tx_counts, mapping, tpm_df=tpm, efflen_df=efflen,
                                 method=AGGREGATION_LENGTH_SCALED_TPM)

        pd.testing.assert_frame_equal(
            dataset.matrix, expected.rename(columns={"Q1": "liver_m1", "Q2": "liver_m2", "Q3": "liver_m3"}),
            check_names=False, check_dtype=False, check_index_type=False, check_column_type=False)
        self.assertEqual(list(quant_summary["sample_id"]), ["liver_m1", "liver_m2", "liver_m3"])
        # 1 vs 2：只给效应量，不做检验
        self.assertEqual(set(table["method"]), {"effect_size_only"})
        self.assertEqual(set(table["evidence_level"]), {"effect_size_only"})
        self.assertTrue(table["pvalue"].isna().all() and (table["padj_method"] == "none").all())
        self.assertEqual(set(table["flags"]), {"low_n"})
        self.assertFalse(table["selected"].any())


class PyDESeq2DependencyGuardTest(unittest.TestCase):
    def _build_counts_project(self, root: Path) -> None:
        exp_dir = root / "experiments" / EXP
        data_dir = exp_dir / "data"
        data_dir.mkdir(parents=True)
        samples = ["liver_m1", "liver_m2", "liver_m3", "liver_m4", "liver_m5", "liver_m6"]
        gene_counts(samples, {"liver_m4", "liver_m5", "liver_m6"}, set(), seed=7, n_genes=40).to_csv(
            data_dir / "counts.csv"
        )
        (root / "project.yaml").write_text("schema_version: 1\nproject: {id: p, species: mus_musculus}\n",
                                           encoding="utf-8")
        (exp_dir / "experiment.yaml").write_text(yaml.safe_dump({
            "schema_version": 1, "experiment": {"id": EXP},
            "design": {
                "factors": {"diet": {"level": "subject", "levels": ["LFD", "HFD"]}},
                "contrasts": [{"id": "HFD_vs_LFD", "factor": "diet", "numerator": "HFD", "denominator": "LFD"}],
            },
            "omics": {"transcriptomics": {"enabled": True, "input": {
                "level": "table", "format": "counts_matrix", "files": [f"experiments/{EXP}/data/counts.csv"],
            }}},
        }, sort_keys=False), encoding="utf-8")
        pd.DataFrame({"subject_id": [f"m{i}" for i in range(1, 7)], "diet": ["LFD"] * 3 + ["HFD"] * 3}).to_csv(
            exp_dir / "subjects.tsv", sep="\t", index=False
        )
        pd.DataFrame({"sample_id": samples, "subject_id": [f"m{i}" for i in range(1, 7)],
                      "omics": "transcriptomics", "tissue": "liver"}).to_csv(
            exp_dir / "samples.tsv", sep="\t", index=False
        )

    def test_missing_pydeseq2_fails_fast_without_wilcoxon_fallback(self) -> None:
        from omicsweave.transcriptomics import runner as tx_runner

        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp).resolve()
            self._build_counts_project(root)
            err = ModuleNotFoundError("No module named 'pydeseq2'", name="pydeseq2")
            with mock.patch.object(tx_runner, "run_pydeseq2_analysis", side_effect=err), \
                    mock.patch.object(tx_runner, "run_wilcoxon_analysis") as wilcoxon_mock:
                with self.assertRaises(RunError) as ctx:
                    _quiet(run_experiment, EXP, ["transcriptomics"], steps=["differential"], root=root)
            self.assertIn("pydeseq2", str(ctx.exception))
            self.assertIn("omicsweave[transcriptomics]", str(ctx.exception))
            wilcoxon_mock.assert_not_called()

    def test_numerical_pydeseq2_failure_still_falls_back_to_wilcoxon(self) -> None:
        from omicsweave.transcriptomics import runner as tx_runner

        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp).resolve()
            self._build_counts_project(root)
            with mock.patch.object(tx_runner, "run_pydeseq2_analysis", side_effect=RuntimeError("dispersion failed")):
                run_dir, notes = _quiet(run_experiment, EXP, ["transcriptomics"], steps=["differential"], root=root)
            table = read_differential(run_dir / "omics" / "transcriptomics" / "tables" / "differential.parquet")
            self.assertEqual(set(table["method"]), {"wilcoxon"})
            self.assertTrue(any("dispersion failed" in note for note in notes))

    def test_validate_warns_when_pydeseq2_is_missing(self) -> None:
        from omicsweave.core.validate import validate_project

        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp).resolve()
            self._build_counts_project(root)
            with mock.patch("omicsweave.core.validate.importlib.util.find_spec", return_value=None):
                report = validate_project(root)
            self.assertTrue(report.ok)
            self.assertTrue(any("pydeseq2" in issue.message for issue in report.warnings))


if __name__ == "__main__":
    unittest.main()
