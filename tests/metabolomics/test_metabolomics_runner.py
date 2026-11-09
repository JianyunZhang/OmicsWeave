#!/usr/bin/env python3
"""代谢组统一运行器：与旧流程逐值一致，并支持一张表含多个组织、用 data_column 对应列名。

核心回归：同一个旧版代谢组实验（两个组织，含 QC 与空白列），旧流程 run_pipeline 与
“omicsweave migrate → omicsweave run --omics metabolomics”的差异结果、PLS-DA VIP 与 KEGG ORA 逐值一致。
"""

from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from omicsweave.contract.dataset import read_dataset
from omicsweave.contract.tables import read_differential, read_enrichment
from omicsweave.core.migrate import migrate_experiment, write_migration
from omicsweave.core.project import project_context
from omicsweave.core.run import run_experiment
from omicsweave.metabolomics.experiment_config import load_metabolomics_config
from omicsweave.metabolomics.metabolomics_pipeline import run_pipeline
from tests.visualization.test_volcano import assert_volcanoes

EXP = "met_exp"
GROUPS = {str(i): ("LFD" if i <= 5 else "HFD") for i in range(1, 11)}


def peak_table(columns: list[str], hfd_columns: set[str], seed: int, n_features: int = 40) -> pd.DataFrame:
    """合成峰面积表：前 8 个代谢物（KEGG ko00010）在 HFD 中升高；含零值、低检出代谢物、QC 与空白列。"""
    rng = np.random.default_rng(seed)
    base = rng.lognormal(10, 1, n_features)
    table = pd.DataFrame({
        "Index": [f"MW{i:04d}" for i in range(n_features)],
        "Compounds": [f"metabolite_{i}" for i in range(n_features)],
        "Formula": ["C6H12O6"] * n_features,
        "kegg_map": ["ko00010,ko01100"] * 8 + ["ko00020"] * 8 + ["ko00030"] * 8 + [""] * (n_features - 24),
    })
    for column in columns:
        effect = np.where(np.arange(n_features) < 8, 2.5, 1.0) if column in hfd_columns else 1.0
        table[column] = base * effect * rng.lognormal(0, 0.15, n_features)
    biological = [c for c in columns if not c.startswith(("QC", "Blank"))]
    for column in [c for c in columns if c.startswith("QC")]:
        table[column] = table[biological].mean(axis=1) * rng.lognormal(0, 0.05, n_features)
    table.loc[30, biological[:4]] = 0.0          # 少量缺失，需插补
    table.loc[31, biological] = 0.0              # 全部缺失，检出过滤去掉
    table.loc[31, biological[0]] = 1.0
    return table


def build_legacy_experiment(root: Path) -> Path:
    met = root / "experiments" / EXP / "metabolomics_data"
    met.mkdir(parents=True)
    for seed, tissue in enumerate(("liver", "heart")):
        columns = [f"{m}_{tissue}" for m in GROUPS] + ["QC1", "QC2", "QC3", "Blank1"]
        hfd = {f"{m}_{tissue}" for m, g in GROUPS.items() if g == "HFD"}
        peak_table(columns, hfd, seed).to_excel(met / f"{tissue}.xlsx", index=False)
    (met / "config.yaml").write_text(yaml.safe_dump({
        "experiment_id": EXP,
        "paths": {"sample_data_liver": "liver.xlsx", "sample_data_heart": "heart.xlsx",
                  "output_dir": "analysis_results/${experiment_id}/metabolomics_results"},
        "params": {"plsda_permutation_num": 20, "random_seed": 7},
        "groups_file": "groups.yaml",
    }), encoding="utf-8")
    (met / "groups.yaml").write_text(yaml.safe_dump({
        "conditions": {"LFD": {"display": "LFD", "reference": True}, "HFD": {"display": "HFD"}},
        "mouse_groups": GROUPS,
    }), encoding="utf-8")
    return met


def _quiet(func, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()):
        return func(*args, **kwargs)


class MigratedRunMatchesLegacyTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        cls.root = root = Path(cls._tmp.name).resolve()
        build_legacy_experiment(root)
        with project_context(root):
            cfg = load_metabolomics_config(experiment_id=EXP)
            _quiet(run_pipeline, cfg, steps=[1, 2, 3, 5, 6])
        cls.legacy_dir = Path(cfg.output_dir)

        write_migration(_quiet(migrate_experiment, EXP, root), root)
        cls.run_dir, cls.notes = _quiet(run_experiment, EXP, ["metabolomics"], root=root)
        tables = cls.run_dir / "omics" / "metabolomics" / "tables"
        cls.differential = read_differential(tables / "differential.parquet")
        cls.enrichment = read_enrichment(tables / "enrichment.parquet")

    @classmethod
    def tearDownClass(cls) -> None:
        cls._tmp.cleanup()

    def legacy_differential(self, tissue: str) -> pd.DataFrame:
        path = self.legacy_dir / "03_differential_analysis" / tissue / f"dm_{tissue}_HFD_vs_LFD.csv"
        return pd.read_csv(path).set_index("Index")

    def new_differential(self, tissue: str) -> pd.DataFrame:
        return self.differential[self.differential["tissue"] == tissue].set_index("feature_id")

    def test_differential_matches_legacy_pipeline(self) -> None:
        self.assertEqual(set(self.differential["contrast_id"]), {"HFD_vs_LFD"})
        for tissue in ("liver", "heart"):
            legacy, new = self.legacy_differential(tissue), self.new_differential(tissue)
            self.assertEqual(sorted(new.index), sorted(legacy.index))
            new = new.loc[legacy.index]
            for legacy_column, column in (("log2FoldChange", "log2fc"), ("statistic", "statistic"),
                                          ("pvalue", "pvalue"), ("padj", "padj"), ("VIP_HFD_vs_LFD", "vip"),
                                          ("VIP_Global", "vip_global")):
                np.testing.assert_allclose(new[column].to_numpy(dtype=float),
                                           legacy[legacy_column].to_numpy(dtype=float), rtol=1e-12, equal_nan=True,
                                           err_msg=f"{tissue} {column}")
            self.assertEqual(list(new["selected"]), list(legacy["significant"].astype(bool)))
            # 有真实效应的代谢物被检出，且均值之差即 log2FC
            self.assertTrue(new.loc["MW0000", "selected"])
            np.testing.assert_allclose(new["mean_numerator"] - new["mean_denominator"], new["log2fc"], atol=1e-12)
            self.assertTrue((new["n_numerator"] == 5).all() and (new["n_denominator"] == 5).all())
            self.assertNotIn("MW0031", new.index)

    def test_ora_matches_legacy_pipeline(self) -> None:
        ora = self.enrichment[self.enrichment["method"] == "ora"]
        for tissue in ("liver", "heart"):
            legacy = pd.read_csv(self.legacy_dir / "05_pathway_enrichment" / tissue
                                 / "Pathway_HFD_vs_LFD.csv").set_index("ID")
            new = ora[ora["tissue"] == tissue].set_index("term_id").loc[legacy.index]
            np.testing.assert_allclose(new["pvalue"], legacy["P-value"], rtol=1e-12)
            np.testing.assert_allclose(new["padj"], legacy["FDR"], rtol=1e-12)
            self.assertEqual(list(new["overlap_size"]), list(legacy["SigCount"]))
            self.assertLess(new.loc["ko00010", "padj"], 0.05)

    def test_gsea_matches_legacy_pipeline(self) -> None:
        gsea = self.enrichment[self.enrichment["method"] == "gsea_prerank"]
        for tissue in ("liver", "heart"):
            legacy = pd.read_csv(self.legacy_dir / "06_gsea" / tissue / "gsea_HFD_vs_LFD_results.csv")
            legacy["term_id"] = legacy["Term"].str.split(" | ", n=1, regex=False).str[0]
            legacy = legacy.set_index("term_id")
            new = gsea[gsea["tissue"] == tissue].set_index("term_id").loc[legacy.index]
            # 旧流程保留 gseapy 报告的 0；统一结果表改为 1/(置换次数+1)（默认 1000 次）并在 flags 中注明
            self.assertEqual(list(new["flags"].str.contains("pvalue_floored")), list(legacy["NOM p-val"] == 0))
            legacy[["NOM p-val", "FDR q-val"]] = legacy[["NOM p-val", "FDR q-val"]].replace(0.0, 1 / 1001)
            for legacy_column, column in (("NES", "nes"), ("NOM p-val", "pvalue"), ("FDR q-val", "padj")):
                np.testing.assert_allclose(new[column].to_numpy(dtype=float),
                                           legacy[legacy_column].to_numpy(dtype=float), rtol=1e-12,
                                           err_msg=f"{tissue} {column}")
            self.assertGreater(new.loc["ko00010", "nes"], 0)
            self.assertNotIn(" | ", "".join(new["term_name"]))

    def test_dataset_qc_and_column_mapping(self) -> None:
        met = self.run_dir / "omics" / "metabolomics"
        heart = read_dataset(met / "data" / "heart")
        self.assertEqual(sorted(heart.matrix.columns), sorted(f"{m}_heart" for m in GROUPS))
        self.assertEqual((heart.meta.scale, heart.meta.normalization), ("log2", "pqn"))
        self.assertIn("Compounds", heart.features.columns)
        # 两个组织的 QC 列同名，迁移时后出现的改为 metabolomics_QC1，经 data_column 对应回数据表
        qc_samples = {tissue: set(pd.read_parquet(met / "qc" / tissue / "sample_totals.parquet")["sample"])
                      for tissue in ("liver", "heart")}
        self.assertEqual(sorted(name for names in qc_samples.values() for name in names if name.endswith("QC1")),
                         ["QC1", "metabolomics_QC1"])
        metrics = pd.read_parquet(met / "qc" / "heart" / "feature_metrics.parquet")
        self.assertEqual(metrics["feature_id"].iloc[0], "MW0000")

    def test_volcano_for_each_tissue_and_contrast(self) -> None:
        assert_volcanoes(self, self.run_dir / "omics" / "metabolomics")


class PeakTableInputTest(unittest.TestCase):
    def test_one_table_with_two_tissues_and_vendor_column_names(self) -> None:
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp).resolve()
            exp_dir = root / "experiments" / EXP
            exp_dir.mkdir(parents=True)
            animals = {f"m{i}": ("LFD" if i <= 3 else "HFD") for i in range(1, 7)}
            columns = [f"{t[0].upper()}-{m}" for t in ("liver", "heart") for m in animals] + ["7_liver"]
            hfd = {f"{t[0].upper()}-{m}" for t in ("liver", "heart") for m, g in animals.items() if g == "HFD"}
            peak_table(columns, hfd, seed=3).rename(columns={"Index": "peak_id", "Compounds": "name"}).to_csv(
                exp_dir / "peaks.tsv", sep="\t", index=False)
            (root / "project.yaml").write_text("schema_version: 1\nproject: {id: p, species: mus_musculus}\n",
                                               encoding="utf-8")
            (exp_dir / "experiment.yaml").write_text(yaml.safe_dump({
                "schema_version": 1,
                "experiment": {"id": EXP},
                "design": {"factors": {"diet": {"level": "subject", "levels": ["LFD", "HFD"]}},
                           "contrasts": [{"id": "HFD_vs_LFD", "factor": "diet", "numerator": "HFD",
                                          "denominator": "LFD", "role": "disease_effect"}]},
                "omics": {"metabolomics": {"enabled": True, "input": {
                    "level": "table", "format": "peak_table", "files": [f"experiments/{EXP}/peaks.tsv"],
                    "feature_id_column": "peak_id", "feature_name_column": "name"}}},
            }, sort_keys=False), encoding="utf-8")
            pd.DataFrame({"subject_id": list(animals), "diet": list(animals.values())}).to_csv(
                exp_dir / "subjects.tsv", sep="\t", index=False)
            pd.DataFrame([{"sample_id": f"{t}_{m}", "subject_id": m, "omics": "metabolomics", "tissue": t,
                           "data_column": f"{t[0].upper()}-{m}"} for t in ("liver", "heart") for m in animals]).to_csv(
                exp_dir / "samples.tsv", sep="\t", index=False)

            run_dir, notes = _quiet(run_experiment, EXP, ["metabolomics"], steps=["differential"], root=root)
            met = run_dir / "omics" / "metabolomics"
            table = read_differential(met / "tables" / "differential.parquet")
            liver = read_dataset(met / "data" / "liver")
            self.assertFalse((met / "tables" / "enrichment.parquet").exists())
            self.assertFalse((met / "qc").exists())

        self.assertEqual(sorted(set(table["tissue"])), ["heart", "liver"])
        self.assertEqual(table["contrast_role"].unique().tolist(), ["disease_effect"])
        self.assertTrue((table["n_numerator"] == 3).all())
        self.assertEqual(sorted(liver.matrix.columns), sorted(f"liver_{m}" for m in animals))
        self.assertEqual(table.loc[table["feature_id"] == "MW0000", "feature_name"].iloc[0], "metabolite_0")
        self.assertTrue(any("7_liver" in n for n in notes))


if __name__ == "__main__":
    unittest.main()
