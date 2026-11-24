#!/usr/bin/env python3
"""同工酶加权（G7）与诱饵映射的假证据率（E3）。

验收（G7）：主要与次要同工酶方向相反的合成反应不被判为支持。
验收（E3）：诱饵集在运行清单中登记随机种子；同一种子的结果可复现。
"""

from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import yaml

from omicsweave.contract.manifest import RunManifest, read_manifest, write_manifest
from omicsweave.contract.tables import write_differential
from omicsweave.core.model import Omics
from omicsweave.integration.evaluation.decoys import _Pool
from omicsweave.integration.integrate import integrate_experiment
from omicsweave.integration.mappings.model_package import CuratedNetworkPackage, DecoyPackage
from omicsweave.integration.methods.flux_evidence import FluxEvidenceBuilder
from tests.integration.test_flux_anchored import contract_tables
from tests.mfa.test_mfa_runner import EXP, build_two_factor_project

C = "HFD_vs_LFD"


def builder(**kwargs) -> FluxEvidenceBuilder:
    return FluxEvidenceBuilder(SimpleNamespace(experiment_id="toy", estimand={}, scoring={
        "evidence_profile": {"mfa_inference_unit": "biological_replicate"}}, **kwargs))


def pyruvate_kinase(pklr_fc: float) -> pd.DataFrame:
    """肝脏丙酮酸激酶：Pklr 为主要同工酶（丰度高），Pkm 为次要同工酶；Pkm 显著升高。"""
    return pd.DataFrame({"protein": ["Pklr", "Pkm"], "abundance": [1000.0, 50.0],
                         f"log2fc_{C}": [pklr_fc, 1.2], f"padj_{C}": [0.4, 0.01], f"pval_{C}": [0.2, 0.001],
                         f"evidence_{C}": "inferential", f"n_eff_{C}": 6})


class IsozymeTest(unittest.TestCase):
    flux = pd.Series({"flux_id": "PYK_c", "tissue": "liver", "flux_log2fc": 0.8})

    def test_minor_isozyme_against_major_is_not_support(self) -> None:
        details = builder().proteomics_evidence(self.flux, pyruvate_kinase(-0.6), C)
        self.assertTrue(details["minor_isozyme_driven"])
        self.assertEqual((details["n_concordant"], details["concordant_proteins"]), (0, []))
        self.assertEqual(details["major_member"], "Pklr")
        self.assertLess(details["weighted_capacity_change"], 0)
        status = builder().evidence_profile_scorer._build_layer_record(
            {"tissue": "liver", "comparison": C, "prot_details": details}, "proteomics")["assessment_status"]
        self.assertNotEqual(status, "support")

    def test_major_isozyme_in_the_same_direction_keeps_support(self) -> None:
        details = builder().proteomics_evidence(self.flux, pyruvate_kinase(0.3), C)
        self.assertFalse(details["minor_isozyme_driven"])
        self.assertEqual(details["concordant_proteins"], ["Pkm"])

    def test_strongly_induced_isozyme_uses_active_condition_share(self) -> None:
        """合并均值占比仅 ~10%、但在活跃条件（上调组）占比超过 20% 的诱导型同工酶不被误判为次要同工酶。"""
        inducible = pd.DataFrame({
            "protein": ["Pklr", "Pkm"],
            "abundance": [880.0, 120.0],
            f"log2fc_{C}": [-0.4, 3.5],
            f"padj_{C}": [0.3, 0.005],
            f"pval_{C}": [0.15, 0.001],
            f"evidence_{C}": "inferential",
            f"n_eff_{C}": 6,
        })
        details = builder().proteomics_evidence(self.flux, inducible, C)
        self.assertFalse(details["minor_isozyme_driven"])
        self.assertEqual(details["concordant_proteins"], ["Pkm"])


class DecoyPackageTest(unittest.TestCase):
    def test_decoy_package_replaces_members_and_keeps_semantics(self) -> None:
        base = CuratedNetworkPackage()
        decoy = DecoyPackage(base, proteins={"PYK_c": ["Actb"]})
        self.assertEqual(decoy.enzyme_mapping("PYK_c"), {"direct_enzymes": ["Actb"]})
        self.assertEqual(decoy.gene_mapping("PYK_c"), {})
        self.assertEqual(decoy.semantics("PYK_c"), base.semantics("PYK_c"))
        details = builder(package=decoy).proteomics_evidence(
            pd.Series({"flux_id": "PYK_c", "tissue": "liver", "flux_log2fc": 0.8}),
            pd.DataFrame({"protein": ["Actb", "Pkm"], f"log2fc_{C}": [1.5, 1.5], f"padj_{C}": [0.01, 0.01],
                          f"evidence_{C}": "inferential"}), C)
        self.assertEqual(details["concordant_proteins"], ["Actb"])

    def test_pool_samples_abundance_matched_non_members(self) -> None:
        frame = pd.DataFrame({"protein": [f"P{i}" for i in range(100)], "abundance": np.arange(1, 101, dtype=float)})
        pool = _Pool(frame, "proteomics")
        rng = np.random.default_rng(0)
        members = [95, 5]
        picked = pool.sample(members, set(members), rng)
        self.assertEqual(len(picked), 2)
        self.assertFalse(set(picked) & set(members))
        self.assertEqual([pool.bins[p] for p in picked], [pool.bins[m] for m in members])


class DecoyBundleTest(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.root = Path(self._tmp.name).resolve()
        build_two_factor_project(self.root)
        source = self.root / "results" / EXP / "run_20261101_000000_src"
        for name, table in contract_tables().items():
            write_differential(table, source / "omics" / name / "tables" / "differential.parquet")
        write_manifest(RunManifest(run_id=source.name, project_id="toy", experiment_id=EXP, omics=list(Omics)), source)
        path = self.root / "experiments" / EXP / "experiment.yaml"
        doc = yaml.safe_load(path.read_text(encoding="utf-8"))
        doc["integration"] = {"enabled": False, "modes": ["flux_anchored"],
                              "parameters": {"flux_anchored": {"decoys": {"replicates": 4, "seed": 7}}}}
        path.write_text(yaml.safe_dump(doc, sort_keys=False), encoding="utf-8")

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def integrate(self) -> tuple[Path, list[str]]:
        with contextlib.redirect_stdout(io.StringIO()):
            return integrate_experiment(EXP, root=self.root)

    def test_decoy_tables_seed_and_reproducibility(self) -> None:
        first, notes = self.integrate()
        second, _ = self.integrate()
        evaluation = first / "integration" / "flux_anchored" / "evaluation"
        detail = pd.read_parquet(evaluation / "decoy_evidence.parquet")
        summary = pd.read_parquet(evaluation / "decoy_summary.parquet")
        self.assertEqual(set(detail["replicate"]), {0, 1, 2, 3})
        self.assertTrue({"transcriptomics", "proteomics", "metabolomics"} <= set(detail["layer"]))
        self.assertTrue(set(summary["subset"]) >= {"all_directional"})
        self.assertTrue(summary["decoy_support_rate"].between(0, 1).all())
        pd.testing.assert_frame_equal(detail, pd.read_parquet(second / "integration" / "flux_anchored" / "evaluation"
                                                              / "decoy_evidence.parquet"))
        self.assertEqual(read_manifest(first).parameters["flux_anchored"]["decoys"], {"replicates": 4, "seed": 7})
        self.assertTrue(any("诱饵" in n for n in notes))


if __name__ == "__main__":
    unittest.main()
