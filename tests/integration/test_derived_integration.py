"""整合读取代谢组派生特征（映射修复阶段 3）：脂质类别总量作底物 / 产物，质量作用比只作关联。"""

from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

from omicsweave.contract.manifest import RunManifest, read_manifest, write_manifest
from omicsweave.contract.tables import write_differential
from omicsweave.core.model import Omics
from omicsweave.integration.integrate import integrate_experiment
from omicsweave.integration.mappings.network_annotation import lipid_pool_roles
from omicsweave.integration.methods.flux_evidence import FluxEvidenceBuilder
from omicsweave.integration.modes.flux_anchored import metabolomics_frames
from tests.integration.test_flux_anchored import DIET, EXP, contract_tables
from tests.mfa.test_mfa_runner import build_two_factor_project

COMP = "case_vs_control"


def derived_table(contrast: str = DIET, tissue: str = "liver") -> pd.DataFrame:
    """派生特征差异表：TG 类别总量升高，乳酸 / 丙酮酸比升高。"""
    rows = [("class_total:TG", "TG class total", "lipid_class_total", "", 2.0, 0.001, 0.01),
            ("class_total:PC", "PC class total", "lipid_class_total", "", 0.1, 0.8, 0.9),
            ("ratio:lactate_pyruvate", "lactate / pyruvate", "mass_action_ratio", "LDH_c", 1.2, 0.002, 0.004)]
    return pd.DataFrame({
        "contrast_id": contrast, "contrast_role": "disease_effect", "tissue": tissue, "omics": "metabolomics",
        "feature_id": [r[0] for r in rows], "feature_name": [r[1] for r in rows], "log2fc": [r[4] for r in rows],
        "statistic": 1.0, "pvalue": [r[5] for r in rows], "padj": [r[6] for r in rows], "padj_method": "fdr_bh",
        "n_numerator": 4, "n_denominator": 4, "method": "welch_ttest", "units": "log2", "inference_unit":
        "biological_replicate", "evidence_level": "statistical_test", "selected": [True, False, True],
        "feature_type": [r[2] for r in rows], "related_reaction": [r[3] for r in rows],
    })


def metab_frame(rows) -> pd.DataFrame:
    """flux_evidence 的代谢组宽表：compound, metabolite_id, feature_type, related_reaction 与按比较的统计列。"""
    return pd.DataFrame({
        "compound": [r[0] for r in rows], "metabolite_id": [r[1] for r in rows],
        "feature_type": [r[2] for r in rows], "related_reaction": [r[3] for r in rows],
        f"log2fc_{COMP}": [r[4] for r in rows], f"pval_{COMP}": [0.001] * len(rows),
        f"fdr_{COMP}": [0.01] * len(rows), f"n_eff_{COMP}": [8] * len(rows),
    })


class LipidPoolTest(unittest.TestCase):
    def test_network_pools_map_to_lipid_class_totals(self) -> None:
        self.assertEqual(lipid_pool_roles("FA_esterification_c"), {"substrates": set(), "products": {"TG"}})
        self.assertEqual(lipid_pool_roles("lipid_mobilization_c"), {"substrates": {"TG"}, "products": set()})
        self.assertEqual(lipid_pool_roles("HEX_c"), {"substrates": set(), "products": set()})


class FramesTest(unittest.TestCase):
    def test_derived_rows_join_the_metabolite_frame_with_their_type(self) -> None:
        frames = metabolomics_frames(contract_tables()["metabolomics"], derived_table())

        liver = frames["liver"].set_index("metabolite_id")
        self.assertEqual(liver.at["class_total:TG", "feature_type"], "lipid_class_total")
        self.assertEqual(liver.at["ratio:lactate_pyruvate", "related_reaction"], "LDH_c")
        self.assertAlmostEqual(liver.at["class_total:TG", f"log2fc_{DIET}"], 2.0)
        self.assertTrue(liver.loc[[i for i in liver.index if ":" not in i], "feature_type"].isna().all())


class EvidenceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.builder = FluxEvidenceBuilder(SimpleNamespace(scoring={}))

    def test_tg_class_total_is_the_product_of_esterification(self) -> None:
        flux_row = pd.Series({"flux_id": "FA_esterification_c", "flux_log2fc": 1.0})
        metab = metab_frame([("TG class total", "class_total:TG", "lipid_class_total", "", 2.0),
                             ("PC class total", "class_total:PC", "lipid_class_total", "", 2.0)])

        details = self.builder.metabolomics_evidence(flux_row, metab, COMP)

        self.assertEqual(details["concordant_metabolites"], ["TG class total"])
        self.assertEqual(details["derived_readouts"], ["TG class total (products)"])

    def test_tg_class_total_is_the_substrate_of_lipolysis(self) -> None:
        flux_row = pd.Series({"flux_id": "lipid_mobilization_c", "flux_log2fc": 1.0})
        metab = metab_frame([("TG class total", "class_total:TG", "lipid_class_total", "", 2.0)])

        details = self.builder.metabolomics_evidence(flux_row, metab, COMP)

        self.assertEqual((details["n_concordant"], details["n_opposite"]), (0, 1))   # 底物与通量同向：冲突

    def test_mass_action_ratio_is_association_only(self) -> None:
        flux_row = pd.Series({"flux_id": "LDH_c", "flux_log2fc": 1.0})
        metab = metab_frame([("lactate / pyruvate", "ratio:lactate_pyruvate", "mass_action_ratio", "LDH_c", 1.2),
                             ("alanine / pyruvate", "ratio:alanine_pyruvate", "mass_action_ratio", "GPT_c", 1.2)])

        details = self.builder.metabolomics_evidence(flux_row, metab, COMP)

        self.assertEqual(details["derived_readouts"], ["lactate / pyruvate (related)"])
        self.assertEqual((details["n_concordant"], details["n_opposite"], details["n_unresolved"]), (0, 0, 1))

    def test_derived_rows_are_not_matched_by_name(self) -> None:
        flux_row = pd.Series({"flux_id": "HEX_c", "flux_log2fc": 1.0})
        metab = metab_frame([("Glucose", np.nan, np.nan, np.nan, -0.8),
                             ("TG class total", "class_total:TG", "lipid_class_total", "", 2.0)])

        details = self.builder.metabolomics_evidence(flux_row, metab, COMP)

        self.assertEqual(details["derived_readouts"], [])
        self.assertEqual(details["concordant_metabolites"], ["Glucose"])


class BundleTest(unittest.TestCase):
    def test_integration_reads_the_derived_table_and_records_it(self) -> None:
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp).resolve()
            build_two_factor_project(root)
            tables = contract_tables()
            run_dir = root / "analysis_results" / EXP / "run_20261101_000000_src"
            for name in ("mfa", "transcriptomics", "proteomics", "metabolomics"):
                write_differential(tables[name], run_dir / "omics" / name / "tables" / "differential.parquet")
            write_differential(derived_table(), run_dir / "omics" / "metabolomics" / "tables"
                               / "derived_differential.parquet")
            write_manifest(RunManifest(run_id=run_dir.name, project_id="toy", experiment_id=EXP,
                                       omics=[Omics(o) for o in ("mfa", "transcriptomics", "proteomics", "metabolomics")]),
                           run_dir)
            with contextlib.redirect_stdout(io.StringIO()):
                out, _ = integrate_experiment(EXP, ["flux_anchored"], root=root)
            manifest = read_manifest(out)

        self.assertTrue(manifest.parameters["flux_anchored"]["derived_metabolomics_features"])
        inputs = {f.path for f in manifest.files if f.role == "input"}
        self.assertIn(f"analysis_results/{EXP}/{run_dir.name}/omics/metabolomics/tables/derived_differential.parquet",
                      inputs)


if __name__ == "__main__":
    unittest.main()
