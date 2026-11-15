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
