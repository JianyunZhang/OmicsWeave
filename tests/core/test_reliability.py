#!/usr/bin/env python3
"""样本 QC 标记 → 组学 × 组织 × 比较的可靠性（G2），及其在 omicsweave run 与 flux_anchored 中的作用。"""

from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from omicsweave.contract.manifest import RunManifest, read_manifest, write_manifest
from omicsweave.contract.tables import read_differential, write_differential
from omicsweave.core.config import load_experiment, load_project
from omicsweave.core.model import Omics
from omicsweave.core.reliability import contrast_reliability, sample_flags
from omicsweave.core.run import run_experiment
from omicsweave.core.validate import validate_project
from omicsweave.integration.integrate import integrate_experiment
from tests.integration.test_flux_anchored import DIET, contract_tables
from tests.mfa.test_mfa_runner import EXP, build_two_factor_project

# build_two_factor_project 的动物：LFD-Saline m0–m2，LFD-GLP1 m3–m5，HFD-Saline m6–m8，HFD-GLP1 m9–m11
LFD_SALINE, HFD_SALINE = ["m0", "m1", "m2"], ["m6", "m7", "m8"]


def quiet(func, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()):
        return func(*args, **kwargs)


class ReliabilityCase(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.root = Path(self._tmp.name).resolve()
        build_two_factor_project(self.root)
        self.samples_path = self.root / "experiments" / EXP / "samples.tsv"

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def set_flags(self, flags: dict[str, str], add_rna: bool = False) -> None:
        """给样本加 qc_flag；``add_rna`` 时另为每只动物加一个肝脏转录组样本（ID 为 rna_<动物>）。"""
        samples = pd.read_csv(self.samples_path, sep="\t", dtype=str)
        if add_rna and not samples["omics"].eq("transcriptomics").any():
            rna = samples.assign(sample_id="rna_" + samples["subject_id"], omics="transcriptomics", tissue="liver")
            samples = pd.concat([samples, rna], ignore_index=True)
        samples["qc_flag"] = samples["sample_id"].map(lambda s: flags.get(s, ""))
        samples["qc_reason"] = samples["qc_flag"].map(lambda f: "离群" if f else "")
        samples.to_csv(self.samples_path, sep="\t", index=False)

    def experiment(self):
        return load_experiment(load_project(self.root), EXP)


class ReliabilityRulesTest(ReliabilityCase):
    def reliability(self, flags: dict[str, str]) -> pd.Series:
        self.set_flags(flags)
        exp = self.experiment()
        keys = pd.DataFrame({"omics": ["mfa"], "tissue": ["liver"], "contrast_id": [DIET]})
        return contrast_reliability(exp, keys, sample_flags(exp)).iloc[0]

    def test_no_flags(self) -> None:
        row = self.reliability({})
        self.assertEqual((row["status"], row["reliability"], row["n_numerator"], row["n_denominator"]),
                         ("ok", 1.0, 3, 3))

    def test_balanced_flags_downweight(self) -> None:
        row = self.reliability({"m6": "flag", "m0": "flag"})
        self.assertEqual(row["status"], "downweighted")
        self.assertAlmostEqual(row["reliability"], 4 / 6, places=4)
        self.assertIn("离群", row["reason"])
        row = self.reliability({"m6": "exclude"})
        self.assertEqual((row["status"], row["n_usable_numerator"]), ("downweighted", 2))

    def test_flags_concentrated_in_one_group_are_confounded(self) -> None:
        for flags in ({s: "flag" for s in LFD_SALINE}, {"m6": "flag", "m7": "flag"}):
            with self.subTest(flags=sorted(flags)):
                row = self.reliability(flags)
                self.assertEqual((row["status"], row["reliability"]), ("qc_confounded", 0.0))

    def test_too_few_usable_samples(self) -> None:
        row = self.reliability({"m6": "flag", "m7": "flag", "m0": "flag", "m1": "flag"})
        self.assertEqual((row["status"], row["n_usable_numerator"], row["n_usable_denominator"]),
                         ("qc_insufficient", 1, 1))

    def test_validate_rejects_unknown_flag(self) -> None:
        self.set_flags({"m0": "bad"})
        messages = [i.message for i in validate_project(self.root).errors]
        self.assertTrue(any("qc_flag" in m for m in messages), messages)


class RunReliabilityTest(ReliabilityCase):
    def test_run_writes_table_and_flags_confounded_rows(self) -> None:
        self.set_flags({s: "flag" for s in LFD_SALINE})
        run_dir, notes = quiet(run_experiment, EXP, ["mfa"], steps=["differential"], root=self.root)
        table = pd.read_parquet(run_dir / "qc" / "contrast_reliability.parquet")
        status = table.groupby("contrast_id")["status"].agg(set)
        self.assertEqual(status[DIET], {"qc_confounded"})
        self.assertEqual(status["GLP1_vs_Saline_in_HFD"], {"ok"})
        diff = read_differential(run_dir / "omics" / "mfa" / "tables" / "differential.parquet")
        flagged = diff["flags"].fillna("").str.split(";").map(lambda f: "qc_confounded" in f)
        self.assertTrue(flagged[diff["contrast_id"] == DIET].all())
        self.assertFalse(flagged[diff["contrast_id"] != DIET].any())
        self.assertTrue(any("qc_confounded" in n for n in notes))
        outputs = {f.path for f in read_manifest(run_dir).files if f.role == "output"}
        self.assertIn("qc/contrast_reliability.parquet", outputs)


class IntegrationReliabilityTest(ReliabilityCase):
    def source_bundle(self) -> Path:
        run_dir = self.root / "analysis_results" / EXP / "run_20261101_000000_src"
        for name, table in contract_tables().items():
            write_differential(table, run_dir / "omics" / name / "tables" / "differential.parquet")
        write_manifest(RunManifest(run_id=run_dir.name, project_id="toy", experiment_id=EXP,
                                   omics=list(Omics)), run_dir)
        return run_dir

    def profile(self) -> tuple[pd.DataFrame, pd.DataFrame, Path, list[str]]:
        run_dir, notes = quiet(integrate_experiment, EXP, ["flux_anchored"], root=self.root)
        tables = run_dir / "integration" / "flux_anchored" / "tables"
        profile = pd.read_parquet(tables / "flux_evidence_profile.parquet").set_index(["comparison", "tissue", "flux_id"])
        evidence = pd.read_parquet(tables / "flux_evidence_table.parquet").set_index(
            ["comparison", "tissue", "flux_id", "layer"])
        return profile, evidence, run_dir, notes

    def test_confounded_layer_gives_no_support_or_conflict(self) -> None:
        self.source_bundle()
        self.set_flags({}, add_rna=True)
        before, before_evidence, _, _ = self.profile()
        self.assertEqual(before_evidence.loc[(DIET, "liver", "HEX_c", "transcriptomics"), "assessment_status"], "support")

        # 把 LFD-Saline 的转录组样本全部标记：转录组在该比较中不可评估，其他层不受影响
        self.set_flags({f"rna_{s}": "flag" for s in LFD_SALINE}, add_rna=True)
        after, evidence, run_dir, notes = self.profile()
        tx = evidence.loc[(DIET, "liver", "HEX_c", "transcriptomics")]
        self.assertEqual((tx["assessment_status"], tx["qc_status"]), ("qc_confounded", "qc_confounded"))
        self.assertEqual((tx["support_evidence"], tx["conflict_evidence"]), (0.0, 0.0))
        self.assertIn("离群", tx["data_reason"])
        self.assertEqual(evidence.loc[(DIET, "liver", "HEX_c", "proteomics"), "assessment_status"], "support")
        self.assertEqual(after.loc[(DIET, "liver", "HEX_c"), "n_supporting_layers"],
                         before.loc[(DIET, "liver", "HEX_c"), "n_supporting_layers"] - 1)
        self.assertTrue(any("qc_confounded" in n for n in notes))
        outputs = {f.path for f in read_manifest(run_dir).files if f.role == "output"}
        self.assertIn("integration/flux_anchored/qc/contrast_reliability.parquet", outputs)

    def test_confounded_anchor_is_not_an_anchor(self) -> None:
        self.source_bundle()
        self.set_flags({s: "flag" for s in LFD_SALINE})
        profile, evidence, _, _ = self.profile()
        row = profile.loc[(DIET, "liver", "HEX_c")]
        self.assertEqual(row["anchor_state"], "qc_confounded")
        self.assertEqual(row["evidence_grade"], "Omics-concordant")
        self.assertEqual(evidence.loc[(DIET, "liver", "HEX_c", "mfa"), "assessment_status"], "qc_confounded")

    def test_manual_factor_overrides_value_not_blocking(self) -> None:
        from omicsweave.integration.methods.evidence_profile import EvidenceProfileScorer
        table = pd.DataFrame({"omics": ["transcriptomics"] * 2, "tissue": ["liver"] * 2,
                              "contrast_id": ["a", "b"], "reliability": [0.5, 0.0],
                              "status": ["downweighted", "qc_confounded"], "reason": ["", "x"]})
        scorer = EvidenceProfileScorer({"quality_reliability": {"transcriptomics": 0.9}}, reliability=table)
        self.assertEqual(scorer._layer_quality("transcriptomics", {"tissue": "liver", "comparison": "a"}),
                         (0.9, "downweighted", ""))
        self.assertEqual(scorer._layer_quality("transcriptomics", {"tissue": "liver", "comparison": "b"})[1],
                         "qc_confounded")
        self.assertEqual(EvidenceProfileScorer(reliability=table)._layer_quality(
            "transcriptomics", {"tissue": "liver", "comparison": "a"})[0], 0.5)


if __name__ == "__main__":
    unittest.main()
