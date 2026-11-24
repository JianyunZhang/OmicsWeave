#!/usr/bin/env python3
"""逆转判定（G4）：五类情形各自被正确归类，低功效的不显著结果不会被判为“未逆转”。"""

from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from omicsweave.contract.manifest import RunManifest, read_manifest, write_manifest
from omicsweave.contract.tables import write_differential
from omicsweave.core.model import Omics
from omicsweave.core.run import RunError
from omicsweave.integration.integrate import integrate_experiment
from omicsweave.integration.modes.reversal import classify, standard_errors
from tests.mfa.test_mfa_runner import EXP, build_two_factor_project

DISEASE, TREATMENT = "HFD_vs_LFD", "GLP1_vs_Saline_in_HFD"
# 基因：疾病效应 (log2fc, se, padj)、处理效应 (log2fc, se, padj)、预期分类
CASES = {
    "Full": ((1.0, 0.10, 0.001), (-1.0, 0.15, 0.001), "reversed"),
    "Partial": ((2.0, 0.10, 0.001), (-0.6, 0.10, 0.001), "partially_reversed"),
    "Over": ((0.8, 0.05, 0.001), (-2.4, 0.20, 0.001), "over_reversed"),
    "Worse": ((1.0, 0.10, 0.001), (0.8, 0.10, 0.001), "aggravated"),
    "Same": ((2.0, 0.10, 0.001), (0.05, 0.10, 0.60), "not_reversed"),
    "Noisy": ((1.0, 0.10, 0.001), (-0.2, 0.60, 0.70), "undetermined"),
    "NoDisease": ((0.3, 0.30, 0.50), (-0.3, 0.10, 0.01), "undetermined"),
}


def differential(cases: dict = CASES) -> pd.DataFrame:
    rows = []
    for gene, (disease, treatment, _) in cases.items():
        for contrast, role, (fc, se, padj) in ((DISEASE, "disease_effect", disease),
                                               (TREATMENT, "treatment_effect", treatment)):
            rows.append({"contrast_id": contrast, "contrast_role": role, "tissue": "liver",
                         "omics": "transcriptomics", "feature_id": gene, "feature_name": gene, "log2fc": fc,
                         "lfc_se": se, "pvalue": padj / 2, "padj": padj, "padj_method": "fdr_bh",
                         "n_numerator": 3, "n_denominator": 3, "method": "deseq2",
                         "evidence_level": "statistical_test"})
    return pd.DataFrame(rows)


class ClassifyTest(unittest.TestCase):
    def merged(self, cases: dict = CASES) -> pd.DataFrame:
        table = differential(cases)
        d = table[table["contrast_id"] == DISEASE].set_index("feature_id")
        t = table[table["contrast_id"] == TREATMENT].set_index("feature_id")
        return pd.DataFrame({
            "log2fc_d": d["log2fc"], "se_d": d["lfc_se"], "padj_d": d["padj"], "evidence_level_d": d["evidence_level"],
            "log2fc_t": t["log2fc"], "se_t": t["lfc_se"], "padj_t": t["padj"], "evidence_level_t": t["evidence_level"],
        })

    def test_each_class(self) -> None:
        result = classify(self.merged(), alpha=0.05, equivalence_margin=0.5, ci_level=0.95)
        self.assertEqual(result["class"].to_dict(), {gene: case[2] for gene, case in CASES.items()})
        self.assertEqual(result.loc["NoDisease", "rule"], "disease_effect_not_established")
        self.assertTrue(np.isnan(result.loc["NoDisease", "reversal_index"]))
        self.assertAlmostEqual(result.loc["Full", "reversal_index"], 1.0)

    def test_not_reversed_requires_equivalence(self) -> None:
        """处理效应同样接近 0，但标准误大到区间超出界限时，只能是“无法判定”。"""
        cases = {"Precise": ((2.0, 0.1, 0.001), (0.05, 0.10, 0.6), "not_reversed"),
                 "Imprecise": ((2.0, 0.1, 0.001), (0.05, 1.50, 0.97), "undetermined")}
        result = classify(self.merged(cases), alpha=0.05, equivalence_margin=0.5, ci_level=0.95)
        self.assertEqual(result["class"].to_dict(), {g: c[2] for g, c in cases.items()})
        self.assertEqual(result.loc["Imprecise", "rule"], "inconclusive")

    def test_untested_treatment_is_undetermined(self) -> None:
        merged = self.merged().assign(evidence_level_t="effect_size_only", padj_t=np.nan)
        result = classify(merged, alpha=0.05, equivalence_margin=0.5, ci_level=0.95)
        self.assertEqual(set(result.loc[result.index != "NoDisease", "rule"]), {"treatment_not_tested"})

    def test_standard_error_from_pvalue(self) -> None:
        table = pd.DataFrame({"log2fc": [1.0, 1.0], "pvalue": [0.05, 0.05], "lfc_se": [0.2, np.nan]})
        se, source = standard_errors(table)
        self.assertAlmostEqual(se[0], 0.2)
        self.assertAlmostEqual(se[1], 1.0 / 1.959964, places=4)
        self.assertEqual(list(source), ["reported", "from_pvalue"])


class ReversalModeTest(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.root = Path(self._tmp.name).resolve()
        build_two_factor_project(self.root)
        run_dir = self.root / "results" / EXP / "run_20261101_000000_src"
        write_differential(differential(), run_dir / "omics" / "transcriptomics" / "tables" / "differential.parquet")
        write_manifest(RunManifest(run_id=run_dir.name, project_id="toy", experiment_id=EXP,
                                   omics=[Omics.TRANSCRIPTOMICS]), run_dir)

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def integrate(self):
        with contextlib.redirect_stdout(io.StringIO()):
            return integrate_experiment(EXP, ["reversal"], root=self.root)

    def test_bundle(self) -> None:
        run_dir, notes = self.integrate()
        base = run_dir / "integration" / "reversal"
        table = pd.read_parquet(base / "tables" / "reversal.parquet").set_index("feature_id")
        self.assertEqual(table["class"].to_dict(), {gene: case[2] for gene, case in CASES.items()})
        self.assertEqual(set(table["pair_id"]), {f"{DISEASE}__{TREATMENT}"})
        self.assertEqual(set(table["se_source"]), {"reported/reported"})
        summary = pd.read_parquet(base / "tables" / "reversal_summary.parquet").iloc[0]
        self.assertEqual((summary["n_features"], summary["undetermined"], summary["reversed"]), (7, 2, 1))
        self.assertTrue((base / "qc" / "contrast_reliability.parquet").is_file())
        self.assertEqual(notes, [])
        self.assertEqual(read_manifest(run_dir).parameters["reversal"]["pairs"], [[DISEASE, TREATMENT]])

    def test_requires_reversal_pairs(self) -> None:
        path = self.root / "experiments" / EXP / "experiment.yaml"
        doc = yaml.safe_load(path.read_text(encoding="utf-8"))
        doc["design"]["reversal_pairs"] = []
        path.write_text(yaml.safe_dump(doc, sort_keys=False), encoding="utf-8")
        with self.assertRaisesRegex(RunError, "reversal_pairs"):
            self.integrate()


if __name__ == "__main__":
    unittest.main()
