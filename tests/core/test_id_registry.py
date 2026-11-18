"""小鼠基因 / 蛋白 ID 注册表的回归测试。

注册表曾经是手写的：在队列 A 蛋白组中能找到的 185 个符号里有 92 个 UniProt 登记号错误，
337 个基因中有 259 个 Ensembl ID 错误，导致 GPD1、MPC1/2 等蛋白在整合中漏配或误配。
这里的“金标准”登记号同时对照了 UniProt 2026_04 与 DIA 搜库结果自带的登记号。
"""

from __future__ import annotations

import unittest
from pathlib import Path

import pandas as pd

from omicsweave.core.idmap.mouse_gene_id_mapping import (
    ENSEMBL_TO_SYMBOL,
    GENE_ALIAS_LOOKUP,
    GENE_ANNOTATIONS,
    REGISTRY_SOURCES,
    get_ensembl_gene_id,
    resolve_gene_symbol,
)
from omicsweave.core.idmap.mouse_protein_id_mapping import (
    MOUSE_PROTEIN_UNIPROT_MAPPING,
    UNIPROT_TO_SYMBOL,
    get_uniprot_id,
)
from omicsweave.integration.mappings import (
    audit_registry_against_data,
    match_gene_features,
    match_protein_features,
)

REPO = Path(__file__).resolve().parents[2]

UNIPROT_GOLDEN = {
    "Gpd1": "P13707", "Gpd1l": "Q3ULJ0", "Gk": "Q64516", "Pgp": "Q8CHP8", "Mpc1": "P63030", "Mpc2": "Q9D023",
    "Slc25a10": "Q9QZD8", "Slc25a11": "Q9CR62", "Slc25a12": "Q8BH59", "Slc25a13": "Q9QXX4", "Aqp7": "O54794",
    "Aqp9": "Q9JJJ3", "Gpt": "Q8QZR5", "Gpt2": "Q8BGT5", "Gls": "D3Z7P3", "Gls2": "Q571F8", "Slc38a2": "Q8CFE6",
    "Pck2": "Q8BH04", "Pcx": "Q05920", "G6pc1": "P35576", "Ces1d": "Q8VCT4",
}
ENSEMBL_GOLDEN = {
    "Hk1": "ENSMUSG00000037012", "Ldha": "ENSMUSG00000063229", "Gpd1": "ENSMUSG00000023019",
    "Gpi1": "ENSMUSG00000036427", "Pcx": "ENSMUSG00000024892", "Glud1": "ENSMUSG00000021794",
}


class RegistryContentTest(unittest.TestCase):
    def test_uniprot_accessions_match_uniprot_and_measured_data(self) -> None:
        for symbol, accession in UNIPROT_GOLDEN.items():
            with self.subTest(symbol=symbol):
                self.assertEqual(get_uniprot_id(symbol), accession)
                self.assertEqual(UNIPROT_TO_SYMBOL[accession], symbol)

    def test_ensembl_gene_ids_match_ensembl_reference(self) -> None:
        for symbol, gene_id in ENSEMBL_GOLDEN.items():
            with self.subTest(symbol=symbol):
                self.assertEqual(get_ensembl_gene_id(symbol), gene_id)
                self.assertEqual(ENSEMBL_TO_SYMBOL[gene_id], symbol)

    def test_uniprot_and_legacy_names_resolve_to_mgi_symbols(self) -> None:
        expected = {"Pc": "Pcx", "Gpi": "Gpi1", "Gyk": "Gk", "G6pc": "G6pc1", "Atgl": "Pnpla2", "Mct1": "Slc16a1"}
        for alias, symbol in expected.items():
            with self.subTest(alias=alias):
                self.assertEqual(resolve_gene_symbol(alias), symbol)
        self.assertIsNone(resolve_gene_symbol("Oaa"))          # 代谢物，不是基因

    def test_genes_absent_from_mouse_are_not_registered(self) -> None:
        # GLUD2 为类人猿特有；小鼠 Aqp10 为假基因
        for symbol in ("Glud2", "Aqp10"):
            with self.subTest(symbol=symbol):
                self.assertNotIn(symbol, GENE_ANNOTATIONS)
                self.assertNotIn(symbol, GENE_ALIAS_LOOKUP)

    def test_registry_invariants(self) -> None:
        for alias, target in GENE_ALIAS_LOOKUP.items():
            self.assertNotIn(alias, GENE_ANNOTATIONS, f"别名 {alias} 是正式符号")
            if target is not None:
                self.assertIn(target, GENE_ANNOTATIONS, f"别名 {alias} 指向未注册符号 {target}")
        for accession, symbol in UNIPROT_TO_SYMBOL.items():
            self.assertIn(symbol, GENE_ANNOTATIONS, accession)
            self.assertIn(accession, MOUSE_PROTEIN_UNIPROT_MAPPING[symbol]["uniprot_ids"])
        self.assertGreater(len(GENE_ANNOTATIONS), 30000)
        self.assertGreater(len(UNIPROT_TO_SYMBOL), 17000)

    def test_sources_are_recorded(self) -> None:
        self.assertRegex(REGISTRY_SOURCES["uniprot"]["release"], r"^\d{4}_\d{2}$")
        self.assertEqual(len(REGISTRY_SOURCES["uniprot"]["sha256"]), 64)
        self.assertEqual(REGISTRY_SOURCES["ensembl"]["assembly"], "GRCm39")


class FocusProteinMatchingTest(unittest.TestCase):
    """整合中曾漏配或误配的蛋白，按数据自带的登记号必须能正确匹配。"""

    def test_previously_missed_or_mismatched_proteins_match_by_accession(self) -> None:
        symbols = list(UNIPROT_GOLDEN)
        proteins = pd.DataFrame({"protein": symbols, "uniprot_id": [UNIPROT_GOLDEN[s] for s in symbols]})

        matched = match_protein_features(proteins, symbols)

        self.assertEqual(set(matched.matched_queries), set(symbols))
        self.assertEqual(set(matched.match_method_by_index.values()), {"uniprot_id"})

    def test_protein_symbol_and_rna_symbol_of_pyruvate_carboxylase_meet(self) -> None:
        protein = match_protein_features(pd.DataFrame({"protein": ["Pc"], "uniprot_id": ["Q05920"]}), ["Pcx"])
        rna = match_gene_features(pd.DataFrame({"gene": ["Pcx"]}), ["Pc"])

        self.assertEqual(protein.matched_queries, ("Pcx",))
        self.assertEqual(rna.matched_queries, ("Pcx",))

    def test_genes_outside_the_old_337_gene_list_are_matchable(self) -> None:
        genes = pd.DataFrame({"gene": ["Pgp", "Gpd1l", "Slc38a4", "Mogat1"]})

        matched = match_gene_features(genes, ["Pgp", "Gpd1l", "Slc38a4", "Mogat1"])

        self.assertEqual(matched.matched_queries, ("Pgp", "Gpd1l", "Slc38a4", "Mogat1"))


def _example_proteomics() -> Path | None:
    # 示例项目，或本地课题工作区（local/projects，真实课题数据放在这里，不提交）
    pattern = "*/analysis_results/*/run_*/omics/proteomics/tables/differential.parquet"
    hits = sorted([*REPO.glob(f"examples/{pattern}"), *REPO.glob(f"local/projects/{pattern}")])
    return hits[-1] if hits else None


class RegistryAgainstDataTest(unittest.TestCase):
    @unittest.skipUnless(_example_proteomics(),
                         "需要示例项目或本地课题的蛋白组结果包（examples/ 或 local/projects/ 下的 analysis_results）")
    def test_registry_agrees_with_search_engine_accessions(self) -> None:
        audit = audit_registry_against_data(pd.read_parquet(_example_proteomics()))

        self.assertGreater(audit.n_checked, 1000)
        # 剩余不一致来自数据侧：蛋白组以一个成员命名、却给出另一个成员的登记号
        self.assertGreater(audit.agreement, 0.99)

    def test_conflicting_rows_are_reported(self) -> None:
        proteins = pd.DataFrame({"feature_id": ["Gpd1", "Arf1", "Novel1"],
                                 "uniprot_id": ["P13707", "P61205", "A0A000XXX0"]})

        audit = audit_registry_against_data(proteins)

        self.assertEqual((audit.n_checked, audit.n_agree), (2, 1))
        self.assertEqual(audit.conflicts, (("Arf1", "P61205", "Arf3"),))


if __name__ == "__main__":
    unittest.main()
