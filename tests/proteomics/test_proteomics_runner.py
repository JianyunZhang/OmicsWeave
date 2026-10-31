#!/usr/bin/env python3
"""蛋白组统一运行器：与旧流程逐值一致；一份报告含多个组织、用 data_column 对应样本名；原始层输入的提示。

核心回归：同一个旧版蛋白组实验（两个组织，批次协变量），旧流程 analyze_tissue（limma，~0 + condition + batch）
与“omicsweave migrate → omicsweave run --omics proteomics”的差异结果逐值一致。
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
from omicsweave.core.genesets import write_gmt
from omicsweave.core.migrate import migrate_experiment, write_migration
from omicsweave.core.project import project_context
from omicsweave.core.run import RunError, run_experiment

EXP = "prot_exp"
ANIMALS = {str(i): ("LFD" if i <= 4 else "HFD") for i in range(1, 9)}
BATCHES = {a: ("b1" if int(a) in (1, 2, 5, 6) else "b2") for a in ANIMALS}


def diann_report(runs: dict[str, tuple[str, str]], seed: int, n_proteins: int = 60) -> pd.DataFrame:
    """DIA-NN 长表：runs 为 Run 名 → (分组, 批次)。前 10 个蛋白在 HFD 中升高，第二批次整体偏高；含缺失与 q 值不合格的行。"""
    rng = np.random.default_rng(seed)
    base = rng.normal(20, 2, n_proteins)
    rows = []
    for run, (group, batch) in runs.items():
        for i in range(n_proteins):
            if rng.random() < 0.08:          # 未检出
                continue
            log2 = base[i] + (1.5 if group == "HFD" and i < 10 else 0.0) + (0.4 if batch == "b2" else 0.0)
            rows.append({"Run": run, "Protein.Group": f"P{i:05d}", "Genes": f"Prot{i:03d}",
                         "PG.MaxLFQ": 2 ** (log2 + rng.normal(0, 0.3)),
                         "PG.Q.Value": 0.02 if rng.random() < 0.02 else 0.001, "Global.PG.Q.Value": 0.001})
    return pd.DataFrame(rows)


def write_protein_sets(path: Path) -> Path:
    """UP_SET 含 10 个上调蛋白；其余为无差异的蛋白块。"""
    proteins = [f"Prot{i:03d}" for i in range(60)]
    return write_gmt({"UP_SET": proteins[:15], "BLOCK_A": proteins[20:35], "BLOCK_B": proteins[35:50]}, path)


def _run_name(tissue: str, animal: str) -> str:
    return f"MSF250153_zzk_23min_100ng_{tissue.capitalize()}_{animal}"


def build_legacy_experiment(root: Path) -> None:
    prot = root / "experiments" / EXP / "proteomics_data"
    for seed, tissue in enumerate(("liver", "heart")):
        runs = {_run_name(tissue, a): (g, BATCHES[a]) for a, g in ANIMALS.items()}
        folder = prot / "diann" / tissue
        folder.mkdir(parents=True)
        diann_report(runs, seed).to_csv(folder / f"report_{tissue}.tsv", sep="\t", index=False)
    (prot / "config.yaml").write_text(yaml.safe_dump({
        "experiment_id": EXP,
        "paths": {"output_dir": f"analysis_results/{EXP}/proteomics_results",
                  "data_dirs": {"liver": "diann/liver", "heart": "diann/heart"}},
        "params": {"tissues": ["liver", "heart"], "run_ssgsea": False, "run_prerank_gsea": True,
                   "gene_sets": [str(write_protein_sets(root / "protein_sets.gmt"))],
                   "differential_design": {"covariates": ["batch"]}},
        "groups_file": "groups.yaml",
    }), encoding="utf-8")
    (prot / "groups.yaml").write_text(yaml.safe_dump({
        "conditions": {"LFD": {"display": "LFD"}, "HFD": {"display": "HFD"}},
        "mouse_groups": ANIMALS,
        "sample_batches": BATCHES,
    }), encoding="utf-8")


def _quiet(func, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()):
        return func(*args, **kwargs)


class MigratedRunMatchesLegacyTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        from omicsweave.proteomics.experiment_config import load_proteomics_config
        from omicsweave.proteomics.run_proteomics_pipeline import analyze_tissue

        cls._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        cls.root = root = Path(cls._tmp.name).resolve()
        build_legacy_experiment(root)
        with project_context(root):
            cfg = load_proteomics_config(experiment_id=EXP)
            for tissue in cfg.tissues:
                _quiet(analyze_tissue, tissue, cfg, force=True, skip_gsea=False)
        cls.legacy_dir = Path(cfg.output_dir) / "01_dp_results"
        cls.legacy_plots = Path(cfg.output_dir) / "02_result_plot"

        write_migration(_quiet(migrate_experiment, EXP, root), root)
        cls.run_dir, cls.notes = _quiet(run_experiment, EXP, ["proteomics"], root=root)
        tables = cls.run_dir / "omics" / "proteomics" / "tables"
        cls.table = read_differential(tables / "differential.parquet")
        cls.enrichment = read_enrichment(tables / "enrichment.parquet")

    @classmethod
    def tearDownClass(cls) -> None:
        cls._tmp.cleanup()

    def test_differential_matches_legacy_pipeline(self) -> None:
        self.assertEqual(set(self.table["contrast_id"]), {"HFD_vs_LFD"})
        for tissue in ("liver", "heart"):
            legacy = pd.read_csv(self.legacy_dir / tissue / f"dp_{tissue}_HFD_vs_LFD.csv",
                                 dtype={"Protein": str}).set_index("Protein")
            new = self.table[self.table["tissue"] == tissue].set_index("feature_id")
            self.assertEqual(sorted(new.index), sorted(legacy.index))
            new = new.loc[legacy.index]
            for legacy_column, column in (("log2FoldChange", "log2fc"), ("standard_error", "lfc_se"),
                                          ("tstat", "statistic"), ("pvalue", "pvalue"), ("padj", "padj")):
                np.testing.assert_allclose(new[column].to_numpy(dtype=float),
                                           legacy[legacy_column].to_numpy(dtype=float),
                                           rtol=1e-12, equal_nan=True, err_msg=f"{tissue} {column}")
            self.assertEqual(list(new["n_numerator"]), list(legacy["n_numerator"].fillna(0).astype(int)))
            self.assertEqual(list(new["selected"]), list(legacy["significant"].astype(bool)))
            self.assertEqual(set(new["design_formula"]), {"~0 + condition + batch"})
            self.assertTrue(new.loc["Prot000", "selected"])
            np.testing.assert_allclose((new["mean_numerator"] - new["mean_denominator"]).loc["Prot001"],
                                       new.loc["Prot001", "log2fc"], atol=0.5)

    def test_gsea_matches_legacy_pipeline(self) -> None:
        from omicsweave.proteomics.enrichment.gsea import find_prerank_report

        gsea = self.enrichment[self.enrichment["method"] == "gsea_prerank"]
        self.assertEqual(set(gsea["gene_set_library"]), {"protein_sets"})
        for tissue in ("liver", "heart"):
            report = find_prerank_report(self.legacy_plots / tissue / "gsea_prerank" / "HFD_vs_LFD" / "protein_sets")
            legacy = pd.read_csv(report).set_index("Term")
            new = gsea[gsea["tissue"] == tissue].set_index("term_name").loc[legacy.index]
            # gseapy 写出的报告文件只保留 6 位小数
            for legacy_column, column in (("NES", "nes"), ("NOM p-val", "pvalue"), ("FDR q-val", "padj")):
                np.testing.assert_allclose(new[column].to_numpy(dtype=float),
                                           legacy[legacy_column].to_numpy(dtype=float), rtol=0, atol=1e-6,
                                           err_msg=f"{tissue} {column}")
            self.assertGreater(new.loc["UP_SET", "nes"], 0)

    def test_dataset_and_qc(self) -> None:
        prot = self.run_dir / "omics" / "proteomics"
        liver = read_dataset(prot / "data" / "liver")
        self.assertEqual(sorted(liver.matrix.columns), sorted(f"Liver{a}" for a in ANIMALS))
        self.assertTrue(liver.matrix.isna().any().any())   # 未插补
        self.assertEqual(liver.features.loc["Prot000", "uniprot_id"], "P00000")
        self.assertTrue((prot / "qc" / "liver" / "qc_summary.json").is_file())
        self.assertEqual(self.notes, [])


class ReportInputTest(unittest.TestCase):
    def test_one_report_with_two_tissues_and_data_column(self) -> None:
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp).resolve()
            exp_dir = root / "experiments" / EXP
            exp_dir.mkdir(parents=True)
            animals = {f"m{i}": ("LFD" if i <= 3 else "HFD") for i in range(1, 7)}
            runs = {_run_name(t, str(i)): (g, "b1") for t in ("liver", "heart")
                    for i, g in enumerate(animals.values(), start=1)}
            runs[_run_name("liver", "9")] = ("LFD", "b1")      # 不在 samples.tsv 中
            diann_report(runs, seed=5).to_csv(exp_dir / "report.tsv", sep="\t", index=False)
            (root / "project.yaml").write_text("schema_version: 1\nproject: {id: p, species: mus_musculus}\n",
                                               encoding="utf-8")
            experiment = {
                "schema_version": 1, "experiment": {"id": EXP},
                "design": {"factors": {"diet": {"level": "subject", "levels": ["LFD", "HFD"]}},
                           "contrasts": [{"id": "HFD_vs_LFD", "factor": "diet", "numerator": "HFD",
                                          "denominator": "LFD", "role": "disease_effect"}]},
                "omics": {"proteomics": {"enabled": True, "input": {
                    "level": "raw", "format": "thermo_raw"}}},
            }
            (exp_dir / "experiment.yaml").write_text(yaml.safe_dump(experiment, sort_keys=False), encoding="utf-8")
            pd.DataFrame({"subject_id": list(animals), "diet": list(animals.values())}).to_csv(
                exp_dir / "subjects.tsv", sep="\t", index=False)
            raw_files = []
            rows = []
            for t in ("liver", "heart"):
                for i, m in enumerate(animals, start=1):
                    raw_files.append(exp_dir / f"{t}_{m}.raw")
                    rows.append({"sample_id": f"{t}_{m}", "subject_id": m, "omics": "proteomics", "tissue": t,
                                 "data_column": f"{t.capitalize()}{i}", "files": f"experiments/{EXP}/{t}_{m}.raw"})
            for path in raw_files:
                path.write_bytes(b"")
            pd.DataFrame(rows).to_csv(exp_dir / "samples.tsv", sep="\t", index=False)

            with self.assertRaisesRegex(RunError, "diann_report"):
                run_experiment(EXP, ["proteomics"], root=root)

            experiment["omics"]["proteomics"]["input"] = {"level": "table", "format": "diann_report",
                                                          "files": [f"experiments/{EXP}/report.tsv"]}
            (exp_dir / "experiment.yaml").write_text(yaml.safe_dump(experiment, sort_keys=False), encoding="utf-8")
            run_dir, notes = _quiet(run_experiment, EXP, ["proteomics"], steps=["differential"], root=root)
            table = read_differential(run_dir / "omics" / "proteomics" / "tables" / "differential.parquet")
            heart = read_dataset(run_dir / "omics" / "proteomics" / "data" / "heart")

        self.assertEqual(sorted(set(table["tissue"])), ["heart", "liver"])
        self.assertEqual(set(table["contrast_role"]), {"disease_effect"})
        self.assertEqual(set(table["design_formula"]), {"~0 + condition"})
        self.assertEqual(sorted(heart.matrix.columns), sorted(f"heart_{m}" for m in animals))
        self.assertTrue(any("Liver9" in n for n in notes))


if __name__ == "__main__":
    unittest.main()
