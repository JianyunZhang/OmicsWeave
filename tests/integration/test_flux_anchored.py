#!/usr/bin/env python3
"""flux_anchored 模式：统一差异结果表整理为证据构建器的输入后，证据与旧流程形状的输入一致；结果包与参数。"""

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

from omicsweave.contract.export import bundle_tables
from omicsweave.contract.manifest import RunManifest, read_manifest, write_manifest
from omicsweave.contract.tables import write_differential
from omicsweave.core.model import Omics
from omicsweave.core.run import RunError
from omicsweave.integration.integrate import integrate_experiment
from omicsweave.integration.methods.flux_evidence import FluxEvidenceBuilder
from omicsweave.integration.modes.flux_anchored import (
    metabolomics_frames, mfa_frames, proteomics_frames, transcriptomics_frames,
)
from omicsweave.metabolomics.export import differential_to_contract
from omicsweave.mfa.export import differential_flux_to_contract
from omicsweave.proteomics.export import dp_to_contract
from omicsweave.transcriptomics.export import de_to_contract
from tests.mfa.test_mfa_runner import EXP, build_two_factor_project

DIET, DRUG = "HFD_vs_LFD", "GLP1_vs_Saline_in_HFD"
ROLES = {DIET: "disease_effect", DRUG: "treatment_effect"}
GROUPS = {DIET: ("HFD", "LFD"), DRUG: ("GLP1", "Saline")}
# 转录组：（比较, 组织）→ 基因、log2FC、padj。药物比较每组 1 只，只给效应量
TX = {
    (DIET, "liver"): [("Hk2", 1.5, 0.001), ("Gck", 0.2, 0.6), ("Aco2", -1.2, 0.01), ("G6pc", 1.4, 0.02),
                      ("Hif1a", 1.0, 0.01), ("Actb", 0.1, 0.9)],
    (DIET, "heart"): [("Pdha1", 1.1, 0.01), ("Dlat", 0.4, 0.3)],
    (DRUG, "liver"): [("Hk2", -1.4, np.nan), ("Aco2", 0.1, np.nan)],
}
# 蛋白组（饮食比较、肝脏）：蛋白、UniProt、log2FC、padj、证据；Aco2 与 ACONT_m 的通量方向相反
PROT = [("Hk2", "", 0.9, 0.01, "inferential"), ("Aco2", "", 0.8, 0.03, "inferential"),
        ("Gck", "P52792", -0.2, 0.7, "inferential"), ("Pdha1", "", 0.6, np.nan, "insufficient")]
# 代谢组（饮食比较、肝脏）：ID、名称、HMDB、KEGG、log2FC、padj
METAB = [("M1", "Glucose", None, "C00031", -0.8, 0.02), ("M2", "Glucose-6-phosphate", None, "C00092", 0.9, 0.01),
         ("M3", "Citrate", "HMDB0000094", "C00158", 1.0, 0.01), ("M4", "Isocitrate", "HMDB0000193", "C00311", 0.7, 0.04),
         ("M5", "Alanine", None, "C00041", 0.1, 0.8)]


def flux_table(contrast: str) -> pd.DataFrame:
    """differential_flux 的输出。饮食比较每组 4 只，主检验为 Mann-Whitney；药物比较每组 1 只，只给效应量。"""
    num, den = GROUPS[contrast]
    if contrast == DIET:
        rows = [("liver", "HEX_c", 1.2, 2.0, 0.002, 0.01), ("liver", "ACONT_m", -0.9, -1.5, 0.004, 0.02),
                ("liver", "G6Pase_c", 0.3, 0.4, 0.40, 0.50), ("heart", "PDH_m", 1.0, 1.3, 0.006, 0.03)]
        n, statistic, primary, level = 4, "mannwhitney", "pvalue_mannwhitney", "statistical_test"
    else:
        rows = [("liver", "HEX_c", -1.1, -1.2, np.nan, np.nan), ("liver", "ACONT_m", 0.2, 0.3, np.nan, np.nan)]
        n, statistic, primary, level = 1, "effect_size_only", "", "effect_size_only"
    table = pd.DataFrame(rows, columns=["tissue", "reaction", "log2FC", "cohen_d", "pvalue_mannwhitney",
                                        "fdr_mannwhitney"])
    return table.assign(
        flux_id=table["tissue"] + "--" + table["reaction"], hedges_g=table["cohen_d"] * 0.9,
        pvalue_ttest=table["pvalue_mannwhitney"] / 2, fdr_ttest=table["fdr_mannwhitney"] / 2, log2fc_se=0.3,
        selection_probability=0.9, direction_stability=0.95, primary_statistic=statistic,
        primary_pvalue_column=primary, evidence_level=level, multiple_testing="fdr_bh",
        is_significant=table["fdr_mannwhitney"] < 0.05,
        **{f"n_bio_{num}": n, f"n_bio_{den}": n, f"mean_{num}": 12.0, f"mean_{den}": 10.0},
    )


def de_results(contrast: str, tissue: str) -> pd.DataFrame:
    """run_pydeseq2_analysis / run_wilcoxon_analysis 的输出，行索引为基因。"""
    table = pd.DataFrame(TX[(contrast, tissue)], columns=["gene", "log2FoldChange", "padj"]).set_index("gene")
    tested = contrast == DIET
    return table.assign(
        pvalue=table["padj"] / 5, lfcSE=0.25, stat=table["log2FoldChange"] / 0.25,
        significant=(table["padj"] < 0.05) & (table["log2FoldChange"].abs() > 1),
        test_method="pydeseq2" if tested else "wilcoxon_effect_only", n_case=4 if tested else 1,
        n_control=4 if tested else 1, evidence_level="statistical_test" if tested else "effect_size_only",
    )


def dp_table() -> pd.DataFrame:
    table = pd.DataFrame(PROT, columns=["Protein", "uniprot_id", "log2FoldChange", "padj", "evidence_level"])
    return table.assign(pvalue=table["padj"] / 4, standard_error=0.2, n_numerator=[4, 4, 4, 1],
                        n_denominator=[4, 4, 4, 1], method="limma_moderated_linear_model", fdr_method="fdr_bh",
                        significant=table["padj"] < 0.05)


def metabolite_table() -> pd.DataFrame:
    table = pd.DataFrame(METAB, columns=["ID", "Compounds", "HMDB", "KEGG_Compound", "log2FoldChange", "padj"])
    return table.assign(pvalue=table["padj"] / 3, statistic=table["log2FoldChange"] * 3,
                        significant=table["padj"] < 0.05, evidence_level="inferential")


def contract_tables() -> dict[str, pd.DataFrame]:
    """各组学的统一差异结果表，由各组学的导出函数生成。"""
    mfa = pd.concat([differential_flux_to_contract(flux_table(c), numerator=GROUPS[c][0], denominator=GROUPS[c][1],
                                                   contrast_id=c, contrast_role=ROLES[c], flux_units="nmol/min/g")
                     for c in (DIET, DRUG)], ignore_index=True)
    tx = pd.concat([de_to_contract(de_results(c, t), contrast_id=c, contrast_role=ROLES[c], tissue=t)
                    for c, t in TX], ignore_index=True)
    dp = dp_table()
    prot = dp_to_contract(dp, contrast_id=DIET, contrast_role=ROLES[DIET], tissue="liver",
                          mean_numerator=pd.Series(20.0, index=dp["Protein"]),
                          mean_denominator=pd.Series(19.0, index=dp["Protein"]), selection_rule="padj<0.05")
    metab_source = metabolite_table()
    metab = differential_to_contract(metab_source, contrast_id=DIET, contrast_role=ROLES[DIET], tissue="liver",
                                     feature_id_column="ID", feature_name_column="Compounds",
                                     mean_numerator=pd.Series(5.0, index=metab_source.index),
                                     mean_denominator=pd.Series(4.0, index=metab_source.index),
                                     n_numerator=4, n_denominator=4, fdr_method="fdr_bh", selection_rule="padj<0.05")
    return {"mfa": mfa, "transcriptomics": tx, "proteomics": prot, "metabolomics": metab}


def legacy_inputs() -> tuple:
    """旧流程加载器给出的形状（MFALoader 等），携带同样的信息。"""
    mfa = {}
    for contrast in (DIET, DRUG):
        t = flux_table(contrast)
        num, den = GROUPS[contrast]
        mfa[contrast] = pd.DataFrame({
            "flux_id": t["reaction"], "tissue": t["tissue"], "flux_log2fc": t["log2FC"],
            "flux_pvalue": t["pvalue_mannwhitney"], "flux_fdr": t["fdr_mannwhitney"], "cohen_d": t["cohen_d"],
            "stat_note": "biological_replicates" if contrast == DIET else "fitting_variance_only",
            "log2fc_se": t["log2fc_se"], "n_eff": t[f"n_bio_{num}"] + t[f"n_bio_{den}"],
            "selection_probability": t["selection_probability"], "direction_stability": t["direction_stability"],
        })
    tx: dict[str, dict[str, pd.DataFrame]] = {}
    for contrast, tissue in TX:
        r = de_results(contrast, tissue)
        tx.setdefault(contrast, {})[tissue] = pd.DataFrame({
            "gene_symbol": r.index, "log2fc": r["log2FoldChange"].to_numpy(), "lfcSE": r["lfcSE"].to_numpy(),
            "pvalue": r["pvalue"].to_numpy(), "padj": r["padj"].to_numpy(),
            "evidence_level": r["evidence_level"].to_numpy(), "n_eff": (r["n_case"] + r["n_control"]).to_numpy(),
        })
    dp = dp_table()
    prot = {"liver": pd.DataFrame({
        "protein": dp["Protein"], "uniprot_id": dp["uniprot_id"], f"log2fc_{DIET}": dp["log2FoldChange"],
        f"pval_{DIET}": dp["pvalue"], f"padj_{DIET}": dp["padj"], f"se_{DIET}": dp["standard_error"],
        f"evidence_{DIET}": dp["evidence_level"], f"n_eff_{DIET}": dp["n_numerator"] + dp["n_denominator"],
    })}
    m = metabolite_table()
    metab = {"liver": pd.DataFrame({
        "compound": m["Compounds"], "hmdb_id": m["HMDB"], "kegg_compound": m["KEGG_Compound"],
        f"log2fc_{DIET}": m["log2FoldChange"], f"pval_{DIET}": m["pvalue"], f"fdr_{DIET}": m["padj"],
        f"n_eff_{DIET}": 8,
    })}
    return mfa, tx, prot, metab


def build(mfa, tx, prot, metab) -> tuple[pd.DataFrame, pd.DataFrame]:
    builder = FluxEvidenceBuilder(SimpleNamespace(
        experiment_id=EXP, estimand={}, scoring={"evidence_profile": {"mfa_inference_unit": "biological_replicate"}}))
    with contextlib.redirect_stdout(io.StringIO()):
        profile = builder.build_all(mfa, tx, prot, metab)
    return profile, builder.get_evidence_table()


def by_key(table: pd.DataFrame, *extra: str) -> pd.DataFrame:
    return table.sort_values(["comparison", "tissue", "flux_id", *extra]).reset_index(drop=True)


class ContractAdapterTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        tables = contract_tables()
        cls.profile, cls.evidence = build(mfa_frames(tables["mfa"]), transcriptomics_frames(tables["transcriptomics"]),
                                          proteomics_frames(tables["proteomics"]),
                                          metabolomics_frames(tables["metabolomics"]))
        cls.legacy_profile, cls.legacy_evidence = build(*legacy_inputs())

    def test_same_evidence_as_legacy_shaped_inputs(self) -> None:
        # 统一结果表另外透传锚点状态的原始依据（旧加载器没有这些列）；由它们得出的 anchor_state 与等级仍须一致
        passthrough = {"mfa_evidence_level", "n_min_group", "mfa_flags", "identifiability", "tracer_evidence"}
        scalar = [c for c in self.legacy_profile.columns if not c.endswith("_details") and c not in passthrough]
        self.assertEqual(sorted(self.profile.columns), sorted(self.legacy_profile.columns))
        pd.testing.assert_frame_equal(by_key(self.profile[scalar]), by_key(self.legacy_profile[scalar]),
                                      check_dtype=False)
        pd.testing.assert_frame_equal(by_key(self.evidence, "layer"), by_key(self.legacy_evidence, "layer"),
                                      check_dtype=False)
        for details, key in (("tx_details", "concordant_genes"), ("prot_details", "concordant_proteins"),
                             ("metab_details", "concordant_metabolites"), ("metab_details", "sig_metabolites")):
            self.assertEqual([d.get(key) for d in by_key(self.profile)[details]],
                             [d.get(key) for d in by_key(self.legacy_profile)[details]])

    def test_evidence_reaches_every_layer(self) -> None:
        rows = self.profile.set_index(["comparison", "tissue", "flux_id"])
        hexokinase = rows.loc[(DIET, "liver", "HEX_c")]
        self.assertTrue(hexokinase["gate_passed"])
        self.assertEqual(hexokinase["tx_details"]["concordant_genes"], ["Hk2"])
        self.assertEqual(hexokinase["prot_details"]["concordant_proteins"], ["Hk2"])
        self.assertEqual(hexokinase["metab_details"]["concordant_metabolites"], ["Glucose", "Glucose-6-phosphate"])
        # 乌头酸酶近平衡：ACO2 与通量反向只记为关联，不再判为冲突
        aconitase = rows.loc[(DIET, "liver", "ACONT_m")]
        self.assertEqual(aconitase["direction_evidence"], "capacity")
        self.assertEqual(aconitase["n_conflicting_layers"], 0)
        self.assertGreaterEqual(aconitase["prot_details"]["n_unresolved"], 1)
        # 药物比较每组 1 只：MFA 与转录组都按效应量判定
        effect_only = rows.loc[(DRUG, "liver", "HEX_c")]
        self.assertEqual((effect_only["mfa_evidence_mode"], effect_only["tx_evidence_mode"]),
                         ("effect_size_only", "effect_size_only"))
        self.assertEqual(effect_only["prot_details"]["data_reason"], "missing_comparison_effect")

    def test_mfa_pvalue_follows_primary_test(self) -> None:
        frame = mfa_frames(contract_tables()["mfa"])[DIET].set_index(["tissue", "flux_id"])
        source = flux_table(DIET).set_index(["tissue", "reaction"])
        np.testing.assert_allclose(frame["flux_fdr"], source["fdr_mannwhitney"])
        np.testing.assert_allclose(frame["flux_se"], 0.3)
        self.assertEqual(set(frame["stat_note"]), {"biological_replicates"})
        self.assertEqual(set(mfa_frames(contract_tables()["mfa"])[DRUG]["stat_note"]), {"fitting_variance_only"})


class FluxAnchoredBundleTest(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.root = Path(self._tmp.name).resolve()
        build_two_factor_project(self.root)
        self.tables = contract_tables()

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def source_bundle(self, omics: tuple[str, ...] = ("mfa", "transcriptomics", "proteomics", "metabolomics"),
                      ) -> Path:
        run_dir = self.root / "analysis_results" / EXP / "run_20261101_000000_src"
        for name in omics:
            write_differential(self.tables[name], run_dir / "omics" / name / "tables" / "differential.parquet")
        write_manifest(RunManifest(run_id=run_dir.name, project_id="toy", experiment_id=EXP,
                                   omics=[Omics(o) for o in omics]), run_dir)
        return run_dir

    def set_integration(self, **integration) -> None:
        path = self.root / "experiments" / EXP / "experiment.yaml"
        doc = yaml.safe_load(path.read_text(encoding="utf-8"))
        # 通路总图（图 E）每页要数秒：合成数据的整合测试默认不出图，test_pathway_map 单独打开
        parameters = dict(integration.get("parameters") or {})
        dossier = parameters.get("reaction_dossier")
        if dossier is None or isinstance(dossier, dict):
            parameters["reaction_dossier"] = {"pathway_map": {"enabled": False}, **(dossier or {})}
            integration = {**integration, "parameters": parameters}
        doc["integration"] = integration
        path.write_text(yaml.safe_dump(doc, sort_keys=False), encoding="utf-8")

    def integrate(self, modes: list[str] | None = None) -> tuple[Path, list[str]]:
        with contextlib.redirect_stdout(io.StringIO()):
            return integrate_experiment(EXP, modes, root=self.root)

    def test_bundle(self) -> None:
        source = self.source_bundle()
        run_dir, notes = self.integrate(["flux_anchored"])
        tables = run_dir / "integration" / "flux_anchored" / "tables"
        profile = pd.read_parquet(tables / "flux_evidence_profile.parquet")
        evidence = pd.read_parquet(tables / "flux_evidence_table.parquet")
        self.assertEqual(len(profile), 6)
        self.assertEqual(len(evidence), 4 * len(profile))
        self.assertFalse([c for c in profile.columns if c.endswith("_details")])
        hexokinase = profile.set_index(["comparison", "tissue", "flux_id"]).loc[(DIET, "liver", "HEX_c")]
        self.assertEqual(hexokinase["contrast_role"], "disease_effect")
        self.assertEqual(hexokinase["metab_concordant_features"], "Glucose;Glucose-6-phosphate")
        # MFA 推断单位取自差异结果（biological_replicate），不降为拟合点层面的 pseudo_replicate
        self.assertNotIn("pseudo_replicate", set(evidence.loc[evidence["layer"] == "mfa", "evidence_mode"]))
        self.assertEqual(notes, [])

        manifest = read_manifest(run_dir)
        self.assertEqual(manifest.parameters["flux_anchored"]["omics"],
                         ["metabolomics", "mfa", "proteomics", "transcriptomics"])
        inputs = {f.path for f in manifest.files if f.role == "input"}
        self.assertIn(f"analysis_results/{EXP}/{source.name}/omics/mfa/tables/differential.parquet", inputs)
        self.assertTrue({"flux_evidence_profile", "flux_evidence_table", "flux_pathway_profile"}
                        <= set(bundle_tables(run_dir)))

    def test_parameters(self) -> None:
        self.source_bundle()
        self.set_integration(enabled=False, modes=["flux_anchored"], parameters={"flux_anchored": {
            "tissue_map": {"heart": "cardiac"}, "evidence_profile": {"mfa_inference_unit": "mcmc_fit_points"}}})
        run_dir, notes = self.integrate()
        tables = run_dir / "integration" / "flux_anchored" / "tables"
        profile = pd.read_parquet(tables / "flux_evidence_profile.parquet")
        evidence = pd.read_parquet(tables / "flux_evidence_table.parquet")
        self.assertIn("cardiac", set(profile["tissue"]))
        self.assertTrue(any("cardiac" in n and "tissue_map" in n for n in notes))
        gated = evidence[(evidence["layer"] == "mfa") & (evidence["comparison"] == DIET)]
        self.assertEqual(set(gated["evidence_mode"]), {"pseudo_replicate"})

        self.set_integration(enabled=False, modes=["flux_anchored"], parameters={"flux_anchored": {"cutoff": 1}})
        with self.assertRaisesRegex(RunError, "无法识别的参数"):
            self.integrate()

    def test_preconditions(self) -> None:
        self.source_bundle(("transcriptomics", "proteomics"))
        with self.assertRaisesRegex(RunError, "需要 MFA"):
            self.integrate(["flux_anchored"])


if __name__ == "__main__":
    unittest.main()
