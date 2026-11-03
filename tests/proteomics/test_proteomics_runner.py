#!/usr/bin/env python3
"""蛋白组统一运行器：与旧流程逐值一致；一份报告含多个组织、用 data_column 对应样本名；原始层输入的提示；
蛋白 × 样本矩阵输入（protein_matrix）与同一数据的 DIA-NN 报告结果一致。

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
from omicsweave.core.model import OmicsInput
from omicsweave.core.project import project_context
from omicsweave.core.run import RunError, run_experiment
from omicsweave.proteomics.run_proteomics_pipeline import load_and_clean_data
from omicsweave.proteomics.runner import read_protein_matrix

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
            # 旧流程保留 gseapy 报告的 0；统一结果表改为 1/(置换次数+1)（默认 1000 次）并在 flags 中注明
            self.assertEqual(list(new["flags"].str.contains("pvalue_floored")), list(legacy["NOM p-val"] == 0))
            legacy[["NOM p-val", "FDR q-val"]] = legacy[["NOM p-val", "FDR q-val"]].replace(0.0, 1 / 1001)
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


def _matrix_spec(**extra) -> OmicsInput:
    return OmicsInput(level="table", format="protein_matrix", **extra)


class MatrixInputTest(unittest.TestCase):
    def test_matrix_matches_report(self) -> None:
        """同一份数据写成 DIA-NN 报告与 log2 宽表（Excel，非第一个工作表），差异结果与数据集一致。"""
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp).resolve()
            exp_dir = root / "experiments" / EXP
            exp_dir.mkdir(parents=True)
            animals = {f"m{i}": ("LFD" if i <= 3 else "HFD") for i in range(1, 7)}
            runs = {_run_name("heart", str(i)): (g, "b1") for i, g in enumerate(animals.values(), start=1)}
            diann_report(runs, seed=11).to_csv(exp_dir / "report.tsv", sep="\t", index=False)
            # 报告聚合后的矩阵写成 log2 宽表，样本列名与报告中提取出的样本名相同
            matrix, to_uniprot = _quiet(load_and_clean_data, str(exp_dir / "report.tsv"))
            wide = np.log2(matrix).rename_axis("Gene.Name").reset_index()
            wide.insert(1, "Protein.Id", wide["Gene.Name"].map(to_uniprot))
            with pd.ExcelWriter(exp_dir / "proteins.xlsx") as writer:
                pd.DataFrame({"Sheet": ["expression_data"]}).to_excel(writer, sheet_name="column_descriptions",
                                                                      index=False)
                wide.to_excel(writer, sheet_name="expression_data", index=False)

            (root / "project.yaml").write_text("schema_version: 1\nproject: {id: p, species: mus_musculus}\n",
                                               encoding="utf-8")
            pd.DataFrame({"subject_id": list(animals), "diet": list(animals.values())}).to_csv(
                exp_dir / "subjects.tsv", sep="\t", index=False)
            pd.DataFrame([{"sample_id": f"heart_{m}", "subject_id": m, "omics": "proteomics", "tissue": "heart",
                           "data_column": f"Heart{i}"} for i, m in enumerate(animals, start=1)]).to_csv(
                exp_dir / "samples.tsv", sep="\t", index=False)
            experiment = {
                "schema_version": 1, "experiment": {"id": EXP},
                "design": {"factors": {"diet": {"level": "subject", "levels": ["LFD", "HFD"]}},
                           "contrasts": [{"id": "HFD_vs_LFD", "factor": "diet", "numerator": "HFD",
                                          "denominator": "LFD", "role": "disease_effect"}]},
                "omics": {"proteomics": {"enabled": True}},
            }
            inputs = {
                "diann_report": {"files": [f"experiments/{EXP}/report.tsv"]},
                "protein_matrix": {"files": [f"experiments/{EXP}/proteins.xlsx"], "sheet": "expression_data",
                                   "feature_column": "Gene.Name", "uniprot_column": "Protein.Id"},
            }
            tables, datasets = {}, {}
            for fmt, extra in inputs.items():
                experiment["omics"]["proteomics"]["input"] = {"level": "table", "format": fmt, **extra}
                (exp_dir / "experiment.yaml").write_text(yaml.safe_dump(experiment, sort_keys=False),
                                                         encoding="utf-8")
                run_dir, notes = _quiet(run_experiment, EXP, ["proteomics"], steps=["differential"], root=root)
                self.assertEqual(notes, [], fmt)
                prot = run_dir / "omics" / "proteomics"
                tables[fmt] = read_differential(prot / "tables" / "differential.parquet").set_index("feature_id")
                datasets[fmt] = read_dataset(prot / "data" / "heart")

        report, from_matrix = tables["diann_report"], tables["protein_matrix"].loc[tables["diann_report"].index]
        self.assertEqual(len(from_matrix), len(tables["protein_matrix"]))
        # log2 → 线性 → log2 的舍入误差（约 1e-16）经 limma 先验自由度的迭代估计放大到约 1e-8
        for column in ("log2fc", "statistic", "pvalue", "padj", "mean_numerator", "mean_denominator"):
            np.testing.assert_allclose(from_matrix[column].to_numpy(dtype=float),
                                       report[column].to_numpy(dtype=float), rtol=1e-6, equal_nan=True,
                                       err_msg=column)
        self.assertEqual(list(from_matrix["selected"]), list(report["selected"]))
        self.assertTrue(from_matrix.loc["Prot000", "selected"])
        pd.testing.assert_frame_equal(datasets["protein_matrix"].matrix, datasets["diann_report"].matrix,
                                      rtol=1e-9, check_names=False)
        self.assertEqual(datasets["protein_matrix"].features.loc["Prot000", "uniprot_id"], "P00000")
        self.assertEqual(datasets["protein_matrix"].meta.source["input_format"], "protein_matrix")

    def test_reader_options(self) -> None:
        """线性强度、默认第一列为特征名、空名改用 UniProt、同名合并（几何均值）、非数值列不当作样本。"""
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            path = Path(tmp) / "proteins.csv"
            pd.DataFrame({
                "gene": ["Alb", "Alb", "", ""],
                "accession": ["P07724;Q3", "P07724", "P01942", ""],
                "description": ["serum albumin", "serum albumin", "hemoglobin", "unknown"],
                "s1": [4.0, 16.0, 10.0, 1.0],
                "s2": [2.0, 0.0, 20.0, 1.0],
            }).to_csv(path, index=False)
            notes: list[str] = []
            matrix, to_uniprot = read_protein_matrix(path, _matrix_spec(scale="linear", uniprot_column="accession"),
                                                     "log2_mean", notes)
            with self.assertRaisesRegex(RunError, "scale"):
                read_protein_matrix(path, _matrix_spec(scale="ln"), "log2_mean", [])
            with self.assertRaisesRegex(RunError, "uniprot_column"):
                read_protein_matrix(path, _matrix_spec(uniprot_column="UniProt"), "log2_mean", [])
            with self.assertRaisesRegex(RunError, "工作表"):
                pd.DataFrame({"gene": ["Alb"], "s1": [1.0]}).to_excel(Path(tmp) / "m.xlsx", index=False)
                read_protein_matrix(Path(tmp) / "m.xlsx", _matrix_spec(sheet="expression_data"), "log2_mean", [])

        self.assertEqual(list(matrix.index), ["Alb", "P01942"])
        self.assertEqual(list(matrix.columns), ["s1", "s2"])
        self.assertAlmostEqual(matrix.loc["Alb", "s1"], 8.0)        # sqrt(4 × 16)
        self.assertAlmostEqual(matrix.loc["Alb", "s2"], 2.0)        # 0 视为未检出
        self.assertEqual(to_uniprot, {"Alb": "P07724", "P01942": "P01942"})
        self.assertTrue(any("1 行没有特征名" in n for n in notes))
        self.assertTrue(any("重复" in n for n in notes))


if __name__ == "__main__":
    unittest.main()
