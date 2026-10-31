#!/usr/bin/env python3
"""基因集库（本地 GMT、按库名下载并缓存）与本地 ORA（单侧 Fisher 检验、背景限制、多重检验校正）。"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
from scipy import stats

from omicsweave.core.config import Project
from omicsweave.core.genesets import GeneSetError, cache_path, load_libraries, read_gmt, write_gmt
from omicsweave.core.model import ProjectConfig
from omicsweave.core.statistics.ora import adjust_pvalues, over_representation


def _project(root: Path) -> Project:
    return Project(root, ProjectConfig.model_validate({"project": {"id": "p", "species": "mus_musculus"}}))


class GeneSetLibraryTest(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name)
        self.project = _project(self.root)

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def test_gmt_round_trip_keeps_order_and_drops_duplicates(self) -> None:
        path = write_gmt({"B set": ["g2", "g1"], "A set": ["g3"]}, self.root / "x.gmt")
        path.write_text(path.read_text(encoding="utf-8") + "C set\tdesc\tg4\tg4\t\n", encoding="utf-8")
        self.assertEqual(read_gmt(path), {"B set": ["g2", "g1"], "A set": ["g3"], "C set": ["g4"]})

    def test_local_gmt_and_cached_library(self) -> None:
        write_gmt({"t1": ["a", "b"]}, self.root / "sets" / "custom.gmt")
        downloaded = {"KEGG term": ["Pdk4", "Cpt1a"]}
        with mock.patch("gseapy.get_library", return_value=downloaded) as get_library:
            libraries, files = load_libraries(self.project, ["sets/custom.gmt", "KEGG_2019_Mouse"])
            # 第二次读取缓存，不再下载
            load_libraries(self.project, ["KEGG_2019_Mouse"])
        get_library.assert_called_once_with(name="KEGG_2019_Mouse", organism="Mouse")
        self.assertEqual(libraries, {"custom": {"t1": ["a", "b"]}, "KEGG_2019_Mouse": downloaded})
        self.assertEqual(files[1], cache_path(self.project, "KEGG_2019_Mouse"))
        self.assertTrue(files[1].is_file())

    def test_download_failure_and_missing_file(self) -> None:
        with mock.patch("gseapy.get_library", side_effect=ConnectionError("offline")):
            with self.assertRaisesRegex(GeneSetError, "本地 GMT"):
                load_libraries(self.project, ["Reactome_Pathways_2024"])
        with self.assertRaisesRegex(GeneSetError, "找不到"):
            load_libraries(self.project, ["missing.gmt"])


class OverRepresentationTest(unittest.TestCase):
    def test_matches_fisher_exact_within_universe(self) -> None:
        universe = [f"g{i}" for i in range(100)]
        query = [f"g{i}" for i in range(10)] + ["not_tested"]
        gene_sets = {
            "enriched": [f"g{i}" for i in range(8)] + [f"g{i}" for i in range(50, 55)] + ["outside1", "outside2"],
            "random": [f"g{i}" for i in range(60, 75)],
            "tiny": ["g1", "g2"],
        }
        table = over_representation(query, universe, gene_sets, min_size=5).set_index("term_name")
        self.assertEqual(set(table.index), {"enriched", "random"})   # tiny 小于 min_size
        row = table.loc["enriched"]
        # 不在背景中的基因（not_tested、outside*）不计入
        self.assertEqual((row["overlap_size"], row["set_size"], row["query_size"], row["universe_size"]),
                         (8, 13, 10, 100))
        expected = stats.fisher_exact([[8, 2], [5, 85]], alternative="greater")[1]
        self.assertAlmostEqual(row["pvalue"], expected, places=12)
        self.assertEqual(row["hit_features"], ";".join(sorted(f"g{i}" for i in range(8))))
        np.testing.assert_allclose(table["padj"], adjust_pvalues(table["pvalue"].to_numpy(), "fdr_bh"))
        self.assertEqual(table.loc["random", "overlap_size"], 0)

    def test_adjust_pvalues_matches_statsmodels(self) -> None:
        from statsmodels.stats.multitest import multipletests

        p = np.random.default_rng(0).uniform(size=40) ** 3
        for method in ("fdr_bh", "fdr_by", "bonferroni"):
            np.testing.assert_allclose(adjust_pvalues(p, method), multipletests(p, method=method)[1], rtol=1e-12)
        with self.assertRaises(ValueError):
            adjust_pvalues(p, "holm")


if __name__ == "__main__":
    unittest.main()
