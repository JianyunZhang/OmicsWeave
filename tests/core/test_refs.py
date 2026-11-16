#!/usr/bin/env python3
"""知识表框架（G7）：RefMet 注册表生成器（离线假响应）、审计与 validate 报告。

验收（G7）：示例项目中登记号与数据不一致的条目被 validate 报出。
"""

from __future__ import annotations

import gzip
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from omicsweave.contract.dataset import Dataset, DatasetMeta, ValueType, write_dataset
from omicsweave.contract.manifest import RunManifest, write_manifest
from omicsweave.core.idmap.build_metabolite_registry import build, resolve, write
from omicsweave.core.idmap.mouse_protein_id_mapping import get_uniprot_ids
from omicsweave.core.model import Omics
from omicsweave.core.refs import REGISTRIES, audit_protein_accessions, registry_status
from omicsweave.core.validate import validate_project
from tests.mfa.test_mfa_runner import EXP, build_two_factor_project

RESPONSES = {
    "refmet/match/Malate": {"refmet_name": "Malic acid", "formula": "C4H6O5", "super_class": "Organic acids"},
    "refmet/name/Malic acid/all": {"refmet_id": "RM0135894", "regno": "37105", "pubchem_cid": "222656"},
    "compound/regno/37105/all": {"hmdb_id": "HMDB0000156", "kegg_id": "C00149", "chebi_id": "6650"},
    "refmet/match/Shuttle": {},
}


def fake_fetch(path: str) -> object:
    return RESPONSES.get(path, {})


class MetaboliteRegistryTest(unittest.TestCase):
    def test_resolve(self) -> None:
        row = resolve("Malate", fake_fetch)
        self.assertEqual((row["refmet_name"], row["hmdb_id"], row["kegg_compound"], row["chebi_id"], row["status"]),
                         ("Malic acid", "HMDB0000156", "C00149", "CHEBI:6650", "matched"))
        self.assertEqual(resolve("Shuttle", fake_fetch)["status"], "unmatched")

    def test_build_cache_and_write(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cache = Path(tmp) / "cache.json"
            rows = build(["Malate", "Shuttle"], fake_fetch, cache=cache, workers=2)
            self.assertIn("refmet/match/Malate", json.loads(cache.read_text(encoding="utf-8")))
            # 第二次运行只用缓存：取数函数一旦被调用就会失败
            again = build(["Malate"], lambda path: (_ for _ in ()).throw(RuntimeError(path)), cache=cache)
            self.assertEqual(again[0]["status"], "matched")
            path = write(rows, Path(tmp))
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                table = pd.read_csv(handle, sep="\t", dtype=str, keep_default_na=False)
            self.assertEqual(table["status"].tolist(), ["matched", "unmatched"])
            sources = json.loads((Path(tmp) / "SOURCES.json").read_text(encoding="utf-8"))
            self.assertEqual(sources["metabolites"]["status"], {"matched": 1, "no_ids": 0, "unmatched": 1})

    def test_registry_listing(self) -> None:
        status = registry_status().set_index("name")
        self.assertEqual(set(status.index), set(REGISTRIES))
        self.assertTrue(status.at["proteins", "present"])


class MetaboliteTableTest(unittest.TestCase):
    """手写代谢物 ID 表中经审计修正的条目（2026-11-16，KEGG / HMDB / ChEBI 逐条核对）。"""

    CORRECTED = {
        "Glyceraldehyde-3-phosphate": {"kegg_compound": "C00118"},
        "HMG-CoA": {"kegg_compound": "C00356", "chebi_id": "CHEBI:11814"},
        "Folate": {"kegg_compound": "C00504"},
        "Malic acid": {"hmdb_id": "HMDB0000156", "kegg_compound": "C00149", "chebi_id": "CHEBI:30797"},
        "2-Phosphoglycerate": {"hmdb_id": "HMDB0003391"},
        "Glycogen": {"hmdb_id": "HMDB0000757"},
        "Heme": {"hmdb_id": "HMDB0003178"},
        "beta-Hydroxybutyrate": {"hmdb_id": "HMDB0000011"},
        "Fructose-2,6-bisphosphate": {"hmdb_id": "HMDB0001047"},
        "Carbamoyl phosphate": {"hmdb_id": "HMDB0001096", "chebi_id": "CHEBI:17672"},
    }

    def test_corrected_entries(self) -> None:
        from omicsweave.core.idmap.metabolite_id_mapping import METABOLITE_DATABASE_IDS, lookup_metabolite_ids
        for name, expected in self.CORRECTED.items():
            with self.subTest(name=name):
                for kind, value in expected.items():
                    self.assertEqual(METABOLITE_DATABASE_IDS[name][kind], value)
        # 错配的编号只留给它们真正对应的化合物
        self.assertEqual(lookup_metabolite_ids("(Rs)-Mevalonic Acid")["kegg_compound"], "C00418")
        self.assertEqual(lookup_metabolite_ids("Homogentisic Acid")["kegg_compound"], "C00544")
        users = {n for n, ids in METABOLITE_DATABASE_IDS.items() if (ids or {}).get("kegg_compound") == "C00354"}
        self.assertFalse(users & {"Glyceraldehyde-3-phosphate"})

    def test_correction_record_matches_the_table(self) -> None:
        from omicsweave.core.idmap.metabolite_id_mapping import CORRECTED_IDS, METABOLITE_DATABASE_IDS
        self.assertEqual(set(CORRECTED_IDS), set(self.CORRECTED))
        for name, fields in CORRECTED_IDS.items():
            for field, (old, new) in fields.items():
                with self.subTest(name=name, field=field):
                    self.assertEqual(METABOLITE_DATABASE_IDS[name][field], new)
                    self.assertNotEqual(old, new)

    def test_no_unreviewed_conflicts_with_refmet(self) -> None:
        from omicsweave.core.refs import ACCEPTED_METABOLITE_VARIANTS, audit_metabolite_table
        self.assertEqual(audit_metabolite_table(), [])
        reviewed = audit_metabolite_table(include_accepted=True)
        self.assertTrue({(r[3], r[7].replace("_conflict", "")) for r in reviewed}
                        <= set(ACCEPTED_METABOLITE_VARIANTS))


class ReannotationTest(unittest.TestCase):
    def test_stale_bundle_ids_follow_the_current_table(self) -> None:
        """旧结果包中由错表写入的苹果酸 ID（D-苹果酸）在整合时按当前表重新注释，与重新运行代谢组一致。"""
        from omicsweave.integration.mappings import get_metabolite_roles, match_metabolite_roles
        from omicsweave.integration.modes.flux_anchored import metabolomics_frames, refresh_metabolite_ids

        stale = pd.DataFrame({
            "contrast_id": "HFD_vs_LFD", "tissue": "liver", "omics": "metabolomics",
            "feature_id": ["MADN0427", "X1"], "feature_name": ["Malic acid", "Unknown peak 7"],
            "hmdb": ["HMDB0031518", "HMDB9999999"], "kegg_compound": ["C00497", ""], "chebi": ["CHEBI:30796", ""],
            "log2fc": [0.5, 0.1], "pvalue": [0.01, 0.5], "padj": [0.02, 0.6], "n_numerator": 3, "n_denominator": 3})
        refreshed, changed = refresh_metabolite_ids(stale)
        self.assertEqual(changed, 1)
        row = refreshed.set_index("feature_id")
        self.assertEqual((row.at["MADN0427", "hmdb"], row.at["MADN0427", "kegg_compound"], row.at["MADN0427", "chebi"]),
                         ("HMDB0000156", "C00149", "CHEBI:30797"))
        self.assertEqual(row.at["X1", "hmdb"], "HMDB9999999")      # 与已知错配无关的注释不动
        roles = get_metabolite_roles("MDH_c")
        for table, expected in ((stale, []), (refreshed, ["Malic acid"])):
            frame = metabolomics_frames(table)["liver"]
            matched = match_metabolite_roles(frame, roles)
            self.assertEqual([frame.at[i, "compound"] for i in matched.indices_by_role["substrates"]], expected)
        self.assertEqual(refresh_metabolite_ids(refreshed)[1], 0)


class ValidateAuditTest(unittest.TestCase):
    def test_validate_reports_accession_mismatch(self) -> None:
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp).resolve()
            build_two_factor_project(root)
            self.assertFalse([i for i in validate_project(root).warnings if "protein_accession" in i.message])
            run_dir = root / "analysis_results" / EXP / "run_20261101_000000_src"
            correct = get_uniprot_ids("Pklr")[0]
            features = pd.DataFrame({"uniprot_id": [correct, "Q00000"]}, index=pd.Index(["Pklr", "Pkm"]))
            write_dataset(run_dir / "omics" / "proteomics" / "data" / "liver", Dataset(
                pd.DataFrame([[1.0, 2.0], [3.0, 4.0]], index=features.index, columns=["s1", "s2"]), features,
                DatasetMeta(omics=Omics.PROTEOMICS, tissue="liver", value_type=ValueType.INTENSITY, scale="log2",
                            feature_id_type="gene")))
            write_manifest(RunManifest(run_id=run_dir.name, project_id="toy", experiment_id=EXP), run_dir)
            rows = audit_protein_accessions(run_dir)
            self.assertEqual([(r[3], r[5], r[7]) for r in rows], [("Pkm", "Q00000", "accession_mismatch")])
            warnings = [i.message for i in validate_project(root).warnings if "protein_accession" in i.message]
            self.assertEqual(len(warnings), 1)
            self.assertIn("Pkm", warnings[0])


if __name__ == "__main__":
    unittest.main()
