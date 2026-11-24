#!/usr/bin/env python3
"""跨实验复现（G6）。

验收：合成的两个实验中，共享效应被召回，单实验特有效应不被标为复现；基线一致率接近 50%。
"""

from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import norm

from omicsweave.cli import main
from omicsweave.contract.manifest import RunManifest, read_manifest, write_manifest
from omicsweave.contract.tables import write_differential
from omicsweave.core.model import Omics, ProjectConfig
from omicsweave.integration.replication import partial_conjunction, replicate_project

N, SHARED, ONLY_A = 2000, 150, 150
SE = 0.3


def differential(true_effect: np.ndarray, rng: np.random.Generator, omics: str = "transcriptomics",
                 ids: list[str] | None = None) -> pd.DataFrame:
    estimate = true_effect + rng.normal(0, SE, len(true_effect))
    p = 2 * norm.sf(np.abs(estimate) / SE)
    ids = ids or [f"Gene{i}" for i in range(len(true_effect))]
    return pd.DataFrame({"contrast_id": "HFD_vs_LFD", "contrast_role": "disease_effect", "tissue": "liver",
                         "omics": omics, "feature_id": ids, "log2fc": estimate, "lfc_se": SE, "pvalue": p,
                         "padj": np.minimum(p * len(p), 1), "padj_method": "fdr_bh", "n_numerator": 4,
                         "n_denominator": 4, "method": "deseq2"})


class ReplicationTest(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.root = Path(self._tmp.name).resolve()
        rng = np.random.default_rng(5)
        truth_a = np.zeros(N)
        truth_a[:SHARED] = 1.5
        truth_a[SHARED:SHARED + ONLY_A] = 1.5
        truth_b = np.zeros(N)
        truth_b[:SHARED] = 1.5
        self.write_bundle("cohort_a", {"transcriptomics": differential(truth_a, rng),
                                       "proteomics": differential(truth_a[:300], rng, "proteomics",
                                                                  [f"Gene{i}" for i in range(300)])})
        self.write_bundle("cohort_b", {"transcriptomics": differential(truth_b, rng)})
        (self.root / "project.yaml").write_text(yaml.safe_dump({
            "schema_version": 1, "project": {"id": "toy", "species": "mus_musculus"},
            "replication": [{"id": "diet", "cross_layer": [["proteomics", "transcriptomics"]], "sources": [
                {"experiment": "cohort_a", "contrast": "HFD_vs_LFD"},
                {"experiment": "cohort_b", "contrast": "HFD_vs_LFD"}]}]}), encoding="utf-8")

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def write_bundle(self, experiment: str, tables: dict[str, pd.DataFrame]) -> None:
        run_dir = self.root / "results" / experiment / "run_20261101_000000_src"
        for omics, table in tables.items():
            write_differential(table, run_dir / "omics" / omics / "tables" / "differential.parquet")
        write_manifest(RunManifest(run_id=run_dir.name, project_id="toy", experiment_id=experiment,
                                   omics=[Omics(o) for o in tables]), run_dir)

    def run_replication(self) -> tuple[Path, pd.DataFrame, pd.DataFrame]:
        run_dir, _ = replicate_project(self.root)
        tables = run_dir / "integration" / "replication" / "tables"
        return run_dir, pd.read_parquet(tables / "replication.parquet"), pd.read_parquet(
            tables / "replication_baseline.parquet")

    def test_shared_effects_replicate_and_specific_effects_do_not(self) -> None:
        run_dir, result, baseline = self.run_replication()
        same = result[(result["match"] == "same_layer") & (result["layers"] == "transcriptomics")].set_index("key")
        shared = same.loc[[f"Gene{i}" for i in range(SHARED)]]
        only_a = same.loc[[f"Gene{i}" for i in range(SHARED, SHARED + ONLY_A)]]
        null = same.loc[[f"Gene{i}" for i in range(SHARED + ONLY_A, N)]]
        self.assertGreater(shared["replicated"].mean(), 0.9)
        self.assertLess(only_a["replicated"].mean(), 0.05)
        self.assertLess(null["replicated"].mean(), 0.01)
        self.assertTrue((shared["direction"] == "up").mean() > 0.95)
        self.assertAlmostEqual(shared["pooled_effect"].mean(), 1.5, delta=0.1)
        rna = baseline[(baseline["match"] == "same_layer") & (baseline["layers"] == "transcriptomics")].iloc[0]
        self.assertAlmostEqual(rna["sign_concordance_random"], 0.5, delta=0.05)
        self.assertGreater(rna["sign_concordance"], rna["sign_concordance_random"])
        sources = read_manifest(run_dir).parameters["replication"]["diet"]
        self.assertEqual([s["run"] for s in sources], ["run_20261101_000000_src"] * 2)

    def test_cross_layer_correspondence(self) -> None:
        _, result, _ = self.run_replication()
        cross = result[result["match"] == "cross_layer"]
        self.assertEqual(set(cross["layers"]), {"proteomics|transcriptomics"})
        hits = cross.set_index("key")
        self.assertGreater(hits.loc[[f"GENE{i}" for i in range(SHARED)], "replicated"].mean(), 0.8)

    def test_cli_and_config(self) -> None:
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            try:
                main(["--project-dir", str(self.root), "replicate", "--set", "diet"])
            except SystemExit as exc:
                self.assertEqual(exc.code, 0)
        self.assertIn("复现组 diet", out.getvalue())
        with self.assertRaisesRegex(Exception, "min_sources"):
            ProjectConfig.model_validate({"project": {"id": "x", "species": "mus_musculus"}, "replication": [
                {"id": "r", "min_sources": 3, "sources": [{"experiment": "a", "contrast": "c"},
                                                          {"experiment": "b", "contrast": "c"}]}]})

    def test_partial_conjunction(self) -> None:
        p = np.array([[0.01, 0.02], [0.01, 0.9], [0.01, 0.02]])
        e = np.array([[1.0, 2.0], [1.0, 1.0], [1.0, -2.0]])
        pc = partial_conjunction(p, e, 2)
        self.assertAlmostEqual(pc[0], 0.02)          # 两个方向一致：2 × max(单侧 p) = 2 × 0.01
        self.assertGreater(pc[1], 0.5)               # 一个不显著
        self.assertGreater(pc[2], 0.9)               # 方向相反


if __name__ == "__main__":
    unittest.main()
