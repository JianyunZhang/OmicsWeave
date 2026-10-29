#!/usr/bin/env python3
"""统一结果表：富集结果表与 MFA 差异结果的规范，以及从各组学现有输出的转换。

转换器的输入尽量来自真实函数：MFA 的 differential_flux、代谢组的 calculate_pathway_enrichment、
gseapy 的 prerank；Enrichr 需要联网，用与 gseapy 返回值同列名的表代替。
"""

from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from omicsweave.adapters.gseapy_results import enrichr_to_contract, prerank_to_contract
from omicsweave.contract import (
    ContractError,
    read_enrichment,
    validate_differential,
    validate_enrichment,
    write_enrichment,
)
from omicsweave.metabolomics.enrichment._run import calculate_pathway_enrichment
from omicsweave.metabolomics.export import ora_to_contract
from omicsweave.mfa.df_analysis import differential_flux, sample_flux_matrix
from omicsweave.mfa.export import differential_flux_to_contract

CONTRAST = dict(contrast_id="CASE_vs_CONTROL", contrast_role="treatment_effect")
SAMPLES = [f"case_{i}" for i in range(3)] + [f"control_{i}" for i in range(3)]
CONDITIONS = {s: s.split("_")[0] for s in SAMPLES}


def _mfa_fits() -> dict[str, pd.DataFrame]:
    fits = {}
    for i, sid in enumerate(SAMPLES):
        case = sid.startswith("case")
        spread = np.linspace(-1.0, 1.0, 20)
        fits[sid] = pd.DataFrame({
            "loss": np.linspace(0.1, 0.2, 20),
            "liver--CS_m": (80.0 if case else 50.0) + i + spread,
            "liver--LDH_c": 30.0 + i + spread,
            "liver--LDH_c__R": 10.0 + spread,
            # 全部动物贴下界 → bound_pinned
            "liver--PEPCK_c": np.full(20, 1.0),
        })
    return fits


def _run_differential(statistic: str) -> pd.DataFrame:
    fits = _mfa_fits()
    flux_df = pd.DataFrame({sid: df.drop(columns="loss").iloc[0] for sid, df in fits.items()}).T
    config = {"annotation": {"require_standard_reaction_ids": False},
              "differential_flux": {"statistic": statistic}}
    with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()):
        sample_flux_matrix(flux_df, fits, CONDITIONS, ["liver"], Path(tmp), analysis_config=config)
        result = differential_flux(flux_df, fits, CONDITIONS, ["case", "control"], ["liver"], Path(tmp),
                                   comparisons=[("case", "control")], analysis_config=config)
    return result["CASE_vs_CONTROL"]


class MfaDifferentialTest(unittest.TestCase):
    def test_effect_size_only_path(self) -> None:
        legacy = _run_differential("auto")  # 3 vs 3：秩检验无法显著，自动只给效应量
        table = differential_flux_to_contract(legacy, numerator="case", denominator="control",
                                              flux_units="nmol/min/g", **CONTRAST)
        self.assertEqual(validate_differential(table), [])
        self.assertTrue((table["method"] == "effect_size_only").all())
        self.assertTrue(table["pvalue"].isna().all())
        self.assertTrue((table["padj_method"] == "none").all())
        self.assertTrue((table["n_numerator"] == 3).all())
        rows = table.set_index("feature_id")
        self.assertGreater(rows.loc["CS_m", "effect_size"], 0)
        self.assertIn("underpowered", rows.loc["CS_m", "flags"])
        # 所有动物都贴在下界：两组方差为 0，Hedges' g 无定义，必须显式标记
        self.assertTrue(np.isnan(rows.loc["PEPCK_c", "effect_size"]))
        self.assertIn("effect_size_undefined", rows.loc["PEPCK_c", "flags"])
        self.assertIn("identifiability:bound_pinned", rows.loc["PEPCK_c", "flags"])

    def test_ttest_path_takes_primary_pvalue_and_fdr(self) -> None:
        legacy = _run_differential("ttest")
        table = differential_flux_to_contract(legacy, numerator="case", denominator="control",
                                              flux_units="nmol/min/g", **CONTRAST)
        self.assertEqual(validate_differential(table), [])
        np.testing.assert_allclose(table["pvalue"], legacy["pvalue_ttest"].to_numpy())
        np.testing.assert_allclose(table["padj"], legacy["fdr_ttest"].to_numpy())
        self.assertTrue((table["padj_method"] == "fdr_bh").all())

    def test_mfa_rows_need_units_but_insufficient_rows_may_lack_effect_size(self) -> None:
        table = differential_flux_to_contract(_run_differential("auto"), numerator="case",
                                              denominator="control", flux_units="nmol/min/g", **CONTRAST)
        broken = table.assign(units="")
        self.assertTrue(any("units 为空" in p for p in validate_differential(broken)))
        insufficient = table.assign(evidence_level="insufficient", effect_size=np.nan)
        self.assertEqual(validate_differential(insufficient), [])
        no_effect = table.assign(effect_size=np.nan, flags="")
        self.assertTrue(any("effect_size 为空" in p for p in validate_differential(no_effect)))


class MetabolomicsOraTest(unittest.TestCase):
    def test_real_ora_output_converts(self) -> None:
        rows = []
        for index in range(10):
            pathways = ["ko00010"] if index < 5 else ["ko00020"]
            rows.append({"Index": f"m{index}", "kegg_map": ",".join(pathways), "Sig": index < 4})
        result = calculate_pathway_enrichment(pd.DataFrame(rows), significance_col="Sig", min_set_size=3,
                                              max_set_size=10, pathway_scope="metabolism", fdr_method="bh")
        table = ora_to_contract(result, tissue="liver", query="all", **CONTRAST)
        self.assertEqual(validate_enrichment(table), [])
        glycolysis = table.set_index("term_id").loc["ko00010"]
        self.assertEqual((glycolysis["overlap_size"], glycolysis["set_size"]), (4, 5))
        self.assertEqual(glycolysis["effect_direction"], "mixed")


class GseapyTest(unittest.TestCase):
    def test_real_prerank_output_converts(self) -> None:
        import gseapy as gp

        rng = np.random.default_rng(0)
        genes = [f"G{i}" for i in range(200)]
        ranks = pd.Series(np.sort(rng.normal(size=200))[::-1], index=genes)
        sets = {"SET_TOP": genes[:30], "SET_BOTTOM": genes[-30:], "SET_MID": genes[80:110]}
        res = gp.prerank(rnk=ranks, gene_sets=sets, min_size=5, max_size=100, permutation_num=20,
                         seed=1, outdir=None, no_plot=True, threads=1, verbose=False)
        table = prerank_to_contract(res.res2d, tissue="liver", omics="transcriptomics",
                                    gene_set_library="custom", **CONTRAST)
        self.assertEqual(validate_enrichment(table), [])
        rows = table.set_index("term_name")
        self.assertEqual(rows.loc["SET_TOP", "effect_direction"], "up")
        self.assertEqual(rows.loc["SET_BOTTOM", "effect_direction"], "down")
        self.assertEqual(rows.loc["SET_TOP", "set_size"], 30)

    def test_prerank_library_prefix_and_missing_library(self) -> None:
        res2d = pd.DataFrame({"Term": ["KEGG_2019_Mouse__Glycolysis"], "ES": [0.5], "NES": [1.8],
                              "NOM p-val": [0.001], "FDR q-val": [0.01], "Tag %": ["12/60"], "Lead_genes": ["A;B"]})
        table = prerank_to_contract(res2d, tissue="liver", omics="proteomics", **CONTRAST)
        self.assertEqual((table.loc[0, "gene_set_library"], table.loc[0, "term_name"]),
                         ("KEGG_2019_Mouse", "Glycolysis"))
        with self.assertRaises(ContractError):
            prerank_to_contract(res2d.assign(Term="Glycolysis"), tissue="liver", omics="proteomics", **CONTRAST)

    def test_enrichr_columns(self) -> None:
        results = pd.DataFrame({
            "Gene_set": ["GO_Biological_Process_2023", "KEGG_2019_Mouse"],
            "Term": ["glycolytic process (GO:0006096)", "Glycolysis / Gluconeogenesis"],
            "Overlap": ["5/120", "4/67"], "P-value": [1e-4, 2e-3], "Adjusted P-value": [0.01, 0.04],
            "Odds Ratio": [8.1, 6.0], "Combined Score": [70.0, 30.0], "Genes": ["A;B;C;D;E", "A;B;C;D"],
        })
        table = enrichr_to_contract(results, tissue="liver", omics="transcriptomics", query="up", **CONTRAST)
        self.assertEqual(validate_enrichment(table), [])
        self.assertEqual(table.loc[0, "term_id"], "GO:0006096")
        self.assertEqual((table.loc[0, "overlap_size"], table.loc[0, "set_size"]), (5, 120))
        with tempfile.TemporaryDirectory() as tmp:
            loaded = read_enrichment(write_enrichment(table, Path(tmp) / "enrichment.parquet"))
        self.assertEqual(list(loaded["term_id"]), ["GO:0006096", ""])


class EnrichmentRulesTest(unittest.TestCase):
    def _base(self) -> pd.DataFrame:
        results = pd.DataFrame({"Gene_set": ["L"], "Term": ["T"], "Overlap": ["2/10"], "P-value": [0.01],
                                "Adjusted P-value": [0.05]})
        return enrichr_to_contract(results, tissue="liver", omics="transcriptomics", query="down", **CONTRAST)

    def test_rules(self) -> None:
        cases = {
            "overlap_size 大于 set_size": dict(overlap_size=pd.array([20], dtype="Int64")),
            "effect_direction 应与 query 对应": dict(effect_direction="up"),
            "method 取值无法识别": dict(method="fisher"),
            "超出 [0, 1]": dict(padj=1.5),
            "缺少必填列": None,
        }
        for fragment, overrides in cases.items():
            with self.subTest(fragment):
                table = self._base().drop(columns="term_name") if overrides is None else self._base().assign(**overrides)
                self.assertTrue(any(fragment in p for p in validate_enrichment(table)), validate_enrichment(table))
        duplicated = pd.concat([self._base(), self._base()], ignore_index=True)
        self.assertTrue(any("重复" in p for p in validate_enrichment(duplicated)))

    def test_gsea_direction_must_follow_nes(self) -> None:
        res2d = pd.DataFrame({"Term": ["S"], "ES": [0.4], "NES": [1.5], "NOM p-val": [0.01], "FDR q-val": [0.1]})
        table = prerank_to_contract(res2d, tissue="liver", omics="transcriptomics", gene_set_library="L", **CONTRAST)
        self.assertEqual(validate_enrichment(table), [])
        self.assertTrue(any("NES 的符号" in p for p in validate_enrichment(table.assign(effect_direction="down"))))


if __name__ == "__main__":
    unittest.main()
