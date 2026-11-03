#!/usr/bin/env python3
"""omicsweave integrate 与 pathway_convergence 模式：读取统一富集结果表，融合结果与旧版分析器一致。"""

from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import combine_pvalues

from omicsweave import cli
from omicsweave.adapters.gseapy_results import prerank_to_contract
from omicsweave.contract.export import bundle_tables, export_bundle
from omicsweave.contract.manifest import RunManifest, read_manifest, write_manifest
from omicsweave.contract.tables import write_enrichment
from omicsweave.core.model import Omics
from omicsweave.core.run import RunError
from omicsweave.integration.integrate import integrate_experiment
from omicsweave.integration.methods.pathway_convergence import PathwayConvergenceAnalyzer, bh_fdr
from omicsweave.integration.modes.pathway_convergence import evidence_table, fusion_table, overlap_table
from omicsweave.metabolomics.export import ora_to_contract
from tests.mfa.test_mfa_runner import EXP, build_two_factor_project

COMPARISON = "HFD_vs_LFD"
PATHWAYS = ["Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)", "Pyruvate metabolism",
            "Fatty acid degradation", "Oxidative phosphorylation"]


def metabolomics_ora(seed: int) -> pd.DataFrame:
    """代谢组 ORA 结果（calculate_pathway_enrichment 的输出列）；FDR 为 P 值的 BH 校正。"""
    rng = np.random.default_rng(seed)
    pvalues = np.array([0.001, 0.004, 0.2, 0.03, 0.6]) * rng.uniform(0.8, 1.2, 5)
    return pd.DataFrame({
        "ID": [f"ko0{i:04d}" for i in range(5)], "Name": PATHWAYS, "P-value": pvalues,
        "FDR": bh_fdr(pvalues), "OddsRatio": rng.uniform(1, 5, 5), "SigCount": [4, 3, 1, 2, 0],
        "BgCount": [10, 12, 8, 9, 15], "SigTotal": 10, "UniverseCount": 60, "Members": "",
    })


def transcriptomics_gsea(seed: int) -> pd.DataFrame:
    """转录组预排序 GSEA 的 res2d（多个库合并检验，条目带库名前缀）。"""
    rng = np.random.default_rng(seed + 10)
    terms = [f"KEGG_2019_Mouse__{p}" for p in PATHWAYS[:4]] + ["GO_Biological_Process_2025__fatty acid oxidation"]
    nes = np.array([2.1, 1.8, -0.4, -1.9, 1.5])
    return pd.DataFrame({
        "Term": terms, "ES": nes / 3, "NES": nes, "NOM p-val": [0.001, 0.01, 0.7, 0.002, 0.05],
        "FDR q-val": np.array([0.01, 0.06, 0.8, 0.02, 0.3]) * rng.uniform(0.9, 1.1, 5),
        "Tag %": ["10/40"] * 5, "Lead_genes": "",
    })


def contract_rows(tissue: str, seed: int) -> dict[str, pd.DataFrame]:
    return {
        "metabolomics": ora_to_contract(metabolomics_ora(seed), contrast_id=COMPARISON, contrast_role="disease_effect",
                                        tissue=tissue, query="all"),
        "transcriptomics": prerank_to_contract(transcriptomics_gsea(seed), contrast_id=COMPARISON,
                                               contrast_role="disease_effect", tissue=tissue,
                                               omics=Omics.TRANSCRIPTOMICS.value),
    }


class MatchesLegacyAnalyzerTest(unittest.TestCase):
    def test_fusion_and_overlaps(self) -> None:
        tissues = ["liver", "heart"]
        with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()):
            analyzer = PathwayConvergenceAnalyzer(Path(tmp))
            legacy_overlaps = analyzer.analyze(
                tissues,
                metab_pathway={t: {COMPARISON: metabolomics_ora(i)} for i, t in enumerate(tissues)},
                prot_pathway={},
                tx_gsea={t: {COMPARISON: transcriptomics_gsea(i)} for i, t in enumerate(tissues)},
            )
        enrichment = pd.concat([table for i, t in enumerate(tissues) for table in contract_rows(t, i).values()],
                               ignore_index=True)
        evidence = evidence_table(enrichment)
        fusion = fusion_table(evidence)
        keys = ["tissue", "comparison", "pathway_term"]
        legacy = analyzer.fusion_df.set_index(keys).sort_index()
        new = fusion.set_index(keys).sort_index()
        self.assertEqual(list(new.index), list(legacy.index))
        for column in ("fused_pvalue", "fused_fdr"):
            np.testing.assert_allclose(new[column], legacy[column], rtol=1e-12)
        for column in ("n_omics", "contributing_omics", "fusion_method"):
            self.assertEqual(list(new[column]), list(legacy[column]))

        overlaps = overlap_table(enrichment, 0.05, 0.25)
        as_set = lambda t: set(map(tuple, t[["tissue", "comparison", "omics_A", "omics_B", "pathway_term"]].values))
        self.assertEqual(as_set(overlaps), as_set(legacy_overlaps))
        self.assertIn(("liver", COMPARISON, "metabolomics", "transcriptomics", "glycolysis / gluconeogenesis"),
                      as_set(overlaps))

        glycolysis = new.loc[("liver", COMPARISON, "glycolysis / gluconeogenesis")]
        self.assertEqual(glycolysis["omics_directions"], "metabolomics:mixed;transcriptomics:up")
        self.assertEqual(glycolysis["contrast_role"], "disease_effect")


class GseaPermutationFloorTest(unittest.TestCase):
    """没有一次置换达到观测富集分数时 gseapy 报告 p = 0；进入融合前应已改为 1/(置换次数+1)。"""

    def test_zero_gsea_pvalue_is_floored_before_fusion(self) -> None:
        res2d = transcriptomics_gsea(0)
        res2d.loc[0, ["NOM p-val", "FDR q-val"]] = 0.0    # Glycolysis / Gluconeogenesis
        floor = 1 / 1001
        metabolomics = ora_to_contract(metabolomics_ora(0), contrast_id=COMPARISON, contrast_role="disease_effect",
                                       tissue="liver", query="all")

        def fuse(permutation_num: int | None) -> tuple[pd.DataFrame, pd.Series]:
            transcriptomics = prerank_to_contract(res2d, contrast_id=COMPARISON, contrast_role="disease_effect",
                                                  tissue="liver", omics=Omics.TRANSCRIPTOMICS.value,
                                                  permutation_num=permutation_num)
            fusion = fusion_table(evidence_table(pd.concat([metabolomics, transcriptomics], ignore_index=True)))
            return transcriptomics, fusion.set_index("pathway_term").loc["glycolysis / gluconeogenesis"]

        table, glycolysis = fuse(1000)
        self.assertEqual((table.loc[0, "pvalue"], table.loc[0, "padj"]), (floor, floor))
        self.assertEqual(table.loc[0, "flags"], "pvalue_floored;padj_floored")
        self.assertEqual(set(table["flags"].iloc[1:]), {""})
        # 融合 p 值以截断后的输入为下界：Fisher 法对每个输入单调，转录组一侧的概率不低于 1/(n+1)
        metabolomics_p = metabolomics.set_index("term_name").loc["Glycolysis / Gluconeogenesis", "padj"]
        np.testing.assert_allclose(glycolysis["fused_pvalue"],
                                   combine_pvalues([floor, metabolomics_p], method="fisher").pvalue, rtol=1e-12)
        self.assertIn(f"'transcriptomics': {round(floor, 6)}", glycolysis["omics_probabilities"])

        # 置换次数未知时保留 0，融合 p 值只受浮点下限约束（真实数据中出现过约 1e-305）
        table, glycolysis = fuse(None)
        self.assertEqual((table.loc[0, "pvalue"], table.loc[0, "flags"]), (0.0, ""))
        self.assertLess(glycolysis["fused_pvalue"], 1e-300)


class IntegrateTest(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.root = Path(self._tmp.name).resolve()
        build_two_factor_project(self.root)

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def source_bundle(self, omics: list[str], run_id: str = "run_20261101_000000_src") -> Path:
        run_dir = self.root / "analysis_results" / EXP / run_id
        for tissue, seed in (("liver", 0), ("heart", 1)):
            rows = contract_rows(tissue, seed)
            for name in omics:
                path = run_dir / "omics" / name / "tables" / "enrichment.parquet"
                existing = pd.read_parquet(path) if path.is_file() else None
                write_enrichment(pd.concat([t for t in (existing, rows[name]) if t is not None], ignore_index=True),
                                 path)
        write_manifest(RunManifest(run_id=run_id, project_id="toy", experiment_id=EXP,
                                   omics=[Omics(o) for o in omics]), run_dir)
        return run_dir

    def set_integration(self, **integration) -> None:
        path = self.root / "experiments" / EXP / "experiment.yaml"
        doc = yaml.safe_load(path.read_text(encoding="utf-8"))
        doc["integration"] = integration
        path.write_text(yaml.safe_dump(doc, sort_keys=False), encoding="utf-8")

    def test_pathway_convergence_bundle(self) -> None:
        source = self.source_bundle(["metabolomics", "transcriptomics"])
        run_dir, notes = integrate_experiment(EXP, ["pathway_convergence"], root=self.root)
        self.assertNotEqual(run_dir, source)
        tables = run_dir / "integration" / "pathway_convergence" / "tables"
        fusion = pd.read_parquet(tables / "pathway_fusion.parquet")
        self.assertIn("direction_agreement", fusion.columns)
        self.assertEqual(set(fusion["tissue"]), {"liver", "heart"})
        manifest = read_manifest(run_dir)
        self.assertEqual(manifest.parameters["source_run"], source.name)
        self.assertEqual(manifest.parameters["pathway_convergence"]["gsea_padj_cutoff"], 0.25)
        inputs = {f.path for f in manifest.files if f.role == "input"}
        self.assertIn(f"analysis_results/{EXP}/{source.name}/omics/metabolomics/tables/enrichment.parquet", inputs)
        self.assertEqual(notes, [])
        # 导出整合结果包
        self.assertEqual(list(bundle_tables(run_dir)), ["pathway_evidence", "pathway_fusion", "pathway_overlaps"])
        sheets = pd.read_excel(export_bundle(run_dir), sheet_name=None)
        self.assertEqual(list(sheets), ["info", "pathway_evidence", "pathway_fusion", "pathway_overlaps"])

    def test_stale_bundle_with_zero_gsea_padj_is_noted(self) -> None:
        source = self.source_bundle(["metabolomics", "transcriptomics"])
        path = source / "omics" / "transcriptomics" / "tables" / "enrichment.parquet"
        write_enrichment(pd.read_parquet(path).assign(padj=lambda t: t["padj"].mask(t.index == 0, 0.0)), path)
        _, notes = integrate_experiment(EXP, ["pathway_convergence"], root=self.root)
        self.assertTrue(any("1 条 GSEA 结果的 padj 为 0" in note for note in notes), notes)

    def test_configured_modes_and_parameters(self) -> None:
        self.source_bundle(["metabolomics", "transcriptomics"])
        self.set_integration(enabled=False, modes=["pathway_convergence"],
                             parameters={"pathway_convergence": {"gsea_padj_cutoff": 0.01}})
        run_dir, _ = integrate_experiment(EXP, root=self.root)
        overlaps = pd.read_parquet(run_dir / "integration" / "pathway_convergence" / "tables"
                                   / "pathway_overlaps.parquet")
        self.assertNotIn("citrate cycle (tca cycle)", set(overlaps["pathway_term"]))   # GSEA FDR ≈ 0.06
        self.assertEqual(read_manifest(run_dir).parameters["pathway_convergence"]["gsea_padj_cutoff"], 0.01)

        self.set_integration(modes=["pathway_convergence"], parameters={"pathway_convergence": {"cutoff": 1}})
        with self.assertRaisesRegex(RunError, "无法识别的参数"):
            integrate_experiment(EXP, root=self.root)

    def test_preconditions_and_errors(self) -> None:
        with self.assertRaisesRegex(RunError, "omicsweave run"):
            integrate_experiment(EXP, ["pathway_convergence"], root=self.root)
        self.source_bundle(["metabolomics"])
        with self.assertRaisesRegex(RunError, "至少 2 个组学"):
            integrate_experiment(EXP, ["pathway_convergence"], root=self.root)
        with self.assertRaisesRegex(RunError, "尚未接入"):
            integrate_experiment(EXP, ["sample_level"], root=self.root)
        with self.assertRaisesRegex(RunError, "未知整合模式"):
            integrate_experiment(EXP, ["network"], root=self.root)
        with self.assertRaisesRegex(RunError, "没有指定整合模式"):
            integrate_experiment(EXP, root=self.root)

    def test_cli(self) -> None:
        source = self.source_bundle(["metabolomics", "transcriptomics"])
        argv = ["--project-dir", str(self.root), "integrate", "--experiment", EXP, "--mode", "pathway_convergence",
                "--source", source.name]
        out = io.StringIO()
        with contextlib.redirect_stdout(out), self.assertRaises(SystemExit) as exit_:
            cli.main(argv)
        self.assertEqual(exit_.exception.code, 0)
        run_dir = Path(out.getvalue().strip().splitlines()[-1].split("：", 1)[1])
        self.assertEqual(read_manifest(run_dir).command, ["omicsweave", *argv])


if __name__ == "__main__":
    unittest.main()
