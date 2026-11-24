#!/usr/bin/env python3
"""实验级与组学级组织别名统一（design.tissue_aliases / omics.<组学>.tissue_aliases）。

验证：
1. 配置模型校验：合并、组学级覆写（含显式自映射）、大小写不敏感匹配、链式展开与环路/非法目标拒绝；
2. validate：分属 sq_adipose 与 white_adipose 的样本按同一组织计生物学重复；
3. 四大单组学统一运行（transcriptomics / proteomics / metabolomics / mfa）：同义组织在样本层合并为单一组织、
   仅输出一个规范组织目录，不会同时产出 sq_adipose 与 white_adipose 两个文件夹；
4. 全部整合模式（flux_anchored / pathway_convergence / reversal / reaction_dossier + panels + pathway_map）、
   标准流程交付（pipeline）及跨队列复现（replicate）一致生效；旧结果包在跨组学改名时可直接整合，组学内撞车时明确报错。
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

from omicsweave.contract.dataset import Dataset, DatasetMeta, ValueType, empty_features, write_dataset
from omicsweave.contract.manifest import RunManifest, read_manifest, write_manifest
from omicsweave.contract.mfa import MFAInfo, MFAResults, write_mfa
from omicsweave.contract.tables import read_differential, read_enrichment, write_differential, write_enrichment
from omicsweave.core.config import ConfigError, load_experiment, load_project
from omicsweave.core.genesets import write_gmt
from omicsweave.core.model import ExperimentConfig, Omics
from omicsweave.core.pipeline import run_pipeline
from omicsweave.core.run import RunError, run_experiment
from omicsweave.core.validate import validate_experiment
from omicsweave.integration.integrate import integrate_experiment
from omicsweave.integration.replication import replicate_project
from tests.metabolomics.test_metabolomics_runner import peak_table
from tests.proteomics.test_proteomics_runner import diann_report

EXP = "alias_exp"
DIET = "HFD_vs_LFD"


def _quiet(func, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        return func(*args, **kwargs)


def _write_project(root: Path, *, replication: list | None = None) -> None:
    doc = {
        "schema_version": 1,
        "project": {"id": "alias_proj", "species": "mus_musculus"},
    }
    if replication is not None:
        doc["replication"] = replication
    (root / "project.yaml").write_text(yaml.safe_dump(doc, allow_unicode=True, sort_keys=False), encoding="utf-8")


def _write_tsv(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)


def _make_mfa_results(samples: list[str], tissue_by_sample: dict[str, list[str]], seed: int = 0) -> MFAResults:
    rng = np.random.default_rng(seed)
    flux_rows, fit_rows = [], []
    for sample in samples:
        treated = 1.0 if "hfd" in sample.lower() else 0.0
        for fit_i in range(4):
            fit_id = str(fit_i)
            fit_rows.append({"sample_id": sample, "fit_id": fit_id, "loss": float(rng.uniform(0.5, 2.0))})
            for tissue in tissue_by_sample[sample]:
                flux_rows.append({
                    "sample_id": sample, "tissue": tissue, "reaction_id": "PDH_m", "direction": "net",
                    "fit_id": fit_id, "value": float(10.0 + 4.0 * treated + rng.normal(0, 0.3)),
                })
                flux_rows.append({
                    "sample_id": sample, "tissue": tissue, "reaction_id": "CS_m", "direction": "net",
                    "fit_id": fit_id, "value": float(8.0 - 2.0 * treated + rng.normal(0, 0.2)),
                })
    return MFAResults(
        info=MFAInfo(model_id="central_carbon_v1", flux_units="nmol/min/g", tracer="glucose"),
        flux=pd.DataFrame(flux_rows),
        fits=pd.DataFrame(fit_rows),
    )


class TissueAliasModelTest(unittest.TestCase):
    def test_merge_override_chain_and_case_insensitivity(self) -> None:
        cfg = ExperimentConfig.model_validate({
            "experiment": {"id": "x"},
            "design": {
                "tissue_aliases": {
                    "inguinal_wat": "white_adipose",
                    "white_adipose": "sq_adipose",
                    "cerebral_cortex": "brain",
                    "renal_cortex": "kidney",
                },
            },
            "omics": {
                "mfa": {
                    "tissue_aliases": {
                        "white_adipose": "white_adipose",  # 显式取消 MFA 的实验级改名
                        "muscle": "gast_muscle",
                    },
                },
                "proteomics": {
                    "tissue_aliases": {
                        "SQ_Adipose": "sq_adipose",
                    },
                },
            },
        })
        # 链式展开 + 大小写不敏感解析
        self.assertEqual(cfg.canonical_tissue("inguinal_wat", Omics.TRANSCRIPTOMICS), "sq_adipose")
        self.assertEqual(cfg.canonical_tissue("White_Adipose", Omics.TRANSCRIPTOMICS), "sq_adipose")
        self.assertEqual(cfg.canonical_tissue("cerebral_cortex", Omics.PROTEOMICS), "brain")
        self.assertEqual(cfg.canonical_tissue("renal_cortex", Omics.METABOLOMICS), "kidney")
        # 组学级显式自映射覆盖实验级改名
        self.assertEqual(cfg.canonical_tissue("white_adipose", Omics.MFA), "white_adipose")
        self.assertEqual(cfg.canonical_tissue("muscle", Omics.MFA), "gast_muscle")
        # 实验级链式已展开为 inguinal_wat → sq_adipose；MFA 中保留 white_adipose 不影响 muscle
        self.assertEqual(cfg.canonical_tissue("cerebral_cortex", Omics.MFA), "brain")

    def test_rejects_cycles_conflicts_and_bad_targets(self) -> None:
        # 组内环路
        with self.assertRaisesRegex(Exception, "循环映射"):
            ExperimentConfig.model_validate({
                "experiment": {"id": "x"},
                "design": {"tissue_aliases": {"sq_adipose": "white_adipose", "white_adipose": "sq_adipose"}},
            })
        # 实验级与组学级跨层环路
        with self.assertRaisesRegex(Exception, "循环映射"):
            ExperimentConfig.model_validate({
                "experiment": {"id": "x"},
                "design": {"tissue_aliases": {"white_adipose": "sq_adipose"}},
                "omics": {"mfa": {"tissue_aliases": {"sq_adipose": "white_adipose"}}},
            })
        # 大小写冲突键
        with self.assertRaisesRegex(Exception, "冲突别名"):
            ExperimentConfig.model_validate({
                "experiment": {"id": "x"},
                "design": {"tissue_aliases": {"white_adipose": "sq_adipose", "White_Adipose": "wat"}},
            })
        # 规范名含路径分隔符
        with self.assertRaisesRegex(Exception, "路径分隔符"):
            ExperimentConfig.model_validate({
                "experiment": {"id": "x"},
                "design": {"tissue_aliases": {"white_adipose": "../sq_adipose"}},
            })

    def test_validate_pools_aliased_tissues_for_replicate_counts(self) -> None:
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp)
            _write_project(root)
            exp_dir = root / "experiments" / EXP
            exp_dir.mkdir(parents=True)
            counts_path = exp_dir / "data" / "counts.tsv"
            counts_path.parent.mkdir(parents=True)
            counts_path.write_text("gene\ts1\ts2\ts3\ts4\nG1\t10\t12\t20\t22\n", encoding="utf-8")
            _write_tsv(exp_dir / "subjects.tsv", pd.DataFrame({
                "subject_id": ["m1", "m2", "m3", "m4"],
                "diet": ["LFD", "LFD", "HFD", "HFD"],
            }))
            # LFD 写 sq_adipose、HFD 写 white_adipose；合并后应为 2 vs 2，不报重复不足
            _write_tsv(exp_dir / "samples.tsv", pd.DataFrame({
                "sample_id": ["s1", "s2", "s3", "s4"],
                "subject_id": ["m1", "m2", "m3", "m4"],
                "omics": ["transcriptomics"] * 4,
                "tissue": ["sq_adipose", "sq_adipose", "white_adipose", "white_adipose"],
            }))
            (exp_dir / "experiment.yaml").write_text(yaml.safe_dump({
                "schema_version": 1,
                "experiment": {"id": EXP},
                "design": {
                    "factors": {"diet": {"level": "subject", "levels": ["LFD", "HFD"], "reference": "LFD"}},
                    "contrasts": [{"id": DIET, "factor": "diet", "numerator": "HFD", "denominator": "LFD",
                                   "role": "disease_effect"}],
                    "stratify_by": ["tissue"],
                    "tissue_aliases": {"white_adipose": "sq_adipose"},
                },
                "omics": {
                    "transcriptomics": {
                        "enabled": True,
                        "input": {"level": "table", "format": "counts_matrix",
                                  "files": [f"experiments/{EXP}/data/counts.tsv"]},
                    },
                },
            }, allow_unicode=True, sort_keys=False), encoding="utf-8")

            exp = load_experiment(load_project(root), EXP)
            self.assertEqual(set(exp.samples["tissue"]), {"sq_adipose"})
            issues = validate_experiment(exp)
            self.assertEqual([i for i in issues if i.level == "error"], [])
            self.assertFalse(any("生物学重复不足" in i.message for i in issues), issues)


class SingleOmicsAliasMergeTest(unittest.TestCase):
    """四大单组学运行器：sq_adipose 与 white_adipose 合并为同一组织分析，且只输出一个文件夹。"""

    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.root = Path(self._tmp.name).resolve()
        _write_project(self.root)
        self.exp_dir = self.root / "experiments" / EXP
        self.exp_dir.mkdir(parents=True)
        self.subjects = pd.DataFrame({
            "subject_id": [f"m{i}" for i in range(1, 9)],
            "diet": ["LFD"] * 4 + ["HFD"] * 4,
        })
        _write_tsv(self.exp_dir / "subjects.tsv", self.subjects)

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def _base_exp_yaml(self, omics: dict, *, design_aliases: dict | None = None) -> None:
        doc = {
            "schema_version": 1,
            "experiment": {"id": EXP},
            "design": {
                "factors": {"diet": {"level": "subject", "levels": ["LFD", "HFD"], "reference": "LFD"}},
                "contrasts": [{"id": DIET, "factor": "diet", "numerator": "HFD", "denominator": "LFD",
                               "role": "disease_effect"}],
                "stratify_by": ["tissue"],
                "tissue_aliases": design_aliases or {},
            },
            "omics": omics,
        }
        (self.exp_dir / "experiment.yaml").write_text(
            yaml.safe_dump(doc, allow_unicode=True, sort_keys=False), encoding="utf-8",
        )

    def test_transcriptomics_merges_two_tissue_files_and_mixed_labels(self) -> None:
        rng = np.random.default_rng(1)
        genes = [f"Gene{i:03d}" for i in range(40)]
        # 两份独立的计数矩阵（sq_adipose 文件含前 2 只 LFD + 前 2 只 HFD；white_adipose 文件含后 2 只 LFD + 后 2 只 HFD）
        # 且特意使用同名数据列 col_1..col_4 验证不会串列
        f_sq = self.exp_dir / "data" / "rna_sq.tsv"
        f_wat = self.exp_dir / "data" / "rna_wat.tsv"
        f_sq.parent.mkdir(parents=True, exist_ok=True)
        for path, shift in ((f_sq, 0), (f_wat, 1)):
            data = {}
            for idx, is_hfd in enumerate((False, False, True, True), start=1):
                base = rng.lognormal(4.5, 0.5, len(genes))
                mean = base * np.where(np.arange(len(genes)) < 8, 4.0 if is_hfd else 1.0, 1.0)
                data[f"col_{idx}"] = rng.negative_binomial(5, 5 / (5 + mean))
            pd.DataFrame(data, index=pd.Index(genes, name="gene")).to_csv(path, sep="\t")

        _write_tsv(self.exp_dir / "samples.tsv", pd.DataFrame({
            "sample_id": ["sq_m1", "sq_m2", "sq_m5", "sq_m6", "wat_m3", "wat_m4", "wat_m7", "wat_m8"],
            "subject_id": ["m1", "m2", "m5", "m6", "m3", "m4", "m7", "m8"],
            "omics": ["transcriptomics"] * 8,
            "tissue": ["sq_adipose"] * 4 + ["white_adipose"] * 4,
            "data_column": ["col_1", "col_2", "col_3", "col_4"] * 2,
        }))
        self._base_exp_yaml(
            {
                "transcriptomics": {
                    "enabled": True,
                    "input": {
                        "level": "table",
                        "format": "counts_matrix",
                        "files": [f"experiments/{EXP}/data/rna_sq.tsv", f"experiments/{EXP}/data/rna_wat.tsv"],
                        "tissue_files": {
                            "sq_adipose": f"experiments/{EXP}/data/rna_sq.tsv",
                            "white_adipose": f"experiments/{EXP}/data/rna_wat.tsv",
                        },
                    },
                    "analysis": {"gene_sets": []},
                },
            },
            design_aliases={"white_adipose": "sq_adipose"},
        )

        run_dir, notes = _quiet(run_experiment, EXP, ["transcriptomics"], root=self.root)
        tx_dir = run_dir / "omics" / "transcriptomics"
        self.assertEqual([p.name for p in (tx_dir / "data").iterdir()], ["sq_adipose"])
        self.assertEqual([p.name for p in (tx_dir / "details").iterdir()], ["sq_adipose"])
        self.assertFalse(any("white_adipose" in str(p) for p in tx_dir.rglob("*")))
        diff = read_differential(tx_dir / "tables" / "differential.parquet")
        self.assertEqual(set(diff["tissue"]), {"sq_adipose"})
        self.assertEqual(set(diff["n_numerator"]), {4})
        self.assertEqual(set(diff["n_denominator"]), {4})
        self.assertTrue(any("white_adipose → sq_adipose" in n for n in notes), notes)

    def test_proteomics_merges_two_diann_reports_into_one_tissue_folder(self) -> None:
        sq_runs = {f"MSF_SQ_{m}": ("LFD" if m in ("m1", "m2") else "HFD", "b1") for m in ("m1", "m2", "m5", "m6")}
        wat_runs = {f"MSF_WAT_{m}": ("LFD" if m in ("m3", "m4") else "HFD", "b1") for m in ("m3", "m4", "m7", "m8")}
        p_sq = self.exp_dir / "data" / "report_sq.tsv"
        p_wat = self.exp_dir / "data" / "report_wat.tsv"
        p_sq.parent.mkdir(parents=True, exist_ok=True)
        diann_report(sq_runs, seed=11, n_proteins=30).to_csv(p_sq, sep="\t", index=False)
        diann_report(wat_runs, seed=12, n_proteins=30).to_csv(p_wat, sep="\t", index=False)
        gmt = write_gmt({"UP_SET": [f"Prot{i:03d}" for i in range(15)]}, self.root / "prot.gmt")

        _write_tsv(self.exp_dir / "samples.tsv", pd.DataFrame({
            "sample_id": ["p_sq_1", "p_sq_2", "p_sq_5", "p_sq_6", "p_wat_3", "p_wat_4", "p_wat_7", "p_wat_8"],
            "subject_id": ["m1", "m2", "m5", "m6", "m3", "m4", "m7", "m8"],
            "omics": ["proteomics"] * 8,
            "tissue": ["sq_adipose"] * 4 + ["white_adipose"] * 4,
            "data_column": ["SQm1", "SQm2", "SQm5", "SQm6", "WATm3", "WATm4", "WATm7", "WATm8"],
        }))
        self._base_exp_yaml(
            {
                "proteomics": {
                    "enabled": True,
                    "input": {
                        "level": "table",
                        "format": "diann_report",
                        "files": [f"experiments/{EXP}/data/report_sq.tsv", f"experiments/{EXP}/data/report_wat.tsv"],
                        "tissue_files": {
                            "sq_adipose": f"experiments/{EXP}/data/report_sq.tsv",
                            "white_adipose": f"experiments/{EXP}/data/report_wat.tsv",
                        },
                    },
                    "analysis": {"gene_sets": [str(gmt)]},
                },
            },
            design_aliases={"white_adipose": "sq_adipose"},
        )

        run_dir, notes = _quiet(run_experiment, EXP, ["proteomics"], root=self.root)
        prot_dir = run_dir / "omics" / "proteomics"
        self.assertEqual([p.name for p in (prot_dir / "data").iterdir()], ["sq_adipose"])
        self.assertEqual([p.name for p in (prot_dir / "qc").iterdir()], ["sq_adipose"])
        self.assertEqual([p.name for p in (prot_dir / "details").iterdir()], ["sq_adipose"])
        self.assertFalse(any("white_adipose" in str(p) for p in prot_dir.rglob("*")))
        diff = read_differential(prot_dir / "tables" / "differential.parquet")
        self.assertEqual(set(diff["tissue"]), {"sq_adipose"})
        self.assertEqual(int(diff["n_numerator"].max()), 4)
        self.assertEqual(int(diff["n_denominator"].max()), 4)
        self.assertTrue(any("已合并来自 2 份输入文件" in n for n in notes), notes)

    def test_metabolomics_merges_two_peak_tables_with_shared_qc_column_names(self) -> None:
        f_sq = self.exp_dir / "data" / "met_sq.tsv"
        f_wat = self.exp_dir / "data" / "met_wat.tsv"
        f_sq.parent.mkdir(parents=True, exist_ok=True)
        sq_cols = ["m1_sq", "m2_sq", "m5_sq", "m6_sq", "QC1", "QC2"]
        wat_cols = ["m3_wat", "m4_wat", "m7_wat", "m8_wat", "QC1", "QC2"]
        peak_table(sq_cols, {"m5_sq", "m6_sq"}, seed=21, n_features=32).to_csv(f_sq, sep="\t", index=False)
        peak_table(wat_cols, {"m7_wat", "m8_wat"}, seed=22, n_features=32).to_csv(f_wat, sep="\t", index=False)

        _write_tsv(self.exp_dir / "samples.tsv", pd.DataFrame({
            "sample_id": ["m1_sq", "m2_sq", "m5_sq", "m6_sq", "sq_QC1", "sq_QC2",
                          "m3_wat", "m4_wat", "m7_wat", "m8_wat", "wat_QC1", "wat_QC2"],
            "subject_id": ["m1", "m2", "m5", "m6", "", "", "m3", "m4", "m7", "m8", "", ""],
            "omics": ["metabolomics"] * 12,
            "tissue": ["sq_adipose"] * 6 + ["white_adipose"] * 6,
            "sample_type": ["sample"] * 4 + ["qc", "qc"] + ["sample"] * 4 + ["qc", "qc"],
            "data_column": sq_cols + wat_cols,
        }))
        self._base_exp_yaml(
            {
                "metabolomics": {
                    "enabled": True,
                    "input": {
                        "level": "table",
                        "format": "peak_table",
                        "files": [f"experiments/{EXP}/data/met_sq.tsv", f"experiments/{EXP}/data/met_wat.tsv"],
                        "tissue_files": {
                            "sq_adipose": f"experiments/{EXP}/data/met_sq.tsv",
                            "white_adipose": f"experiments/{EXP}/data/met_wat.tsv",
                        },
                    },
                    "analysis": {"plsda_permutation_num": 10, "random_seed": 7, "derived_features": False},
                },
            },
            design_aliases={"white_adipose": "sq_adipose"},
        )

        run_dir, notes = _quiet(run_experiment, EXP, ["metabolomics"], root=self.root)
        met_dir = run_dir / "omics" / "metabolomics"
        self.assertEqual([p.name for p in (met_dir / "data").iterdir()], ["sq_adipose"])
        self.assertEqual([p.name for p in (met_dir / "qc").iterdir()], ["sq_adipose"])
        self.assertEqual([p.name for p in (met_dir / "details").iterdir()], ["sq_adipose"])
        self.assertFalse(any("white_adipose" in str(p) for p in met_dir.rglob("*")))
        diff = read_differential(met_dir / "tables" / "differential.parquet")
        self.assertEqual(set(diff["tissue"]), {"sq_adipose"})
        self.assertEqual(set(diff["n_numerator"]), {4})
        self.assertEqual(set(diff["n_denominator"]), {4})
        self.assertTrue(any("已合并来自 2 张数据表" in n for n in notes), notes)

    def test_mfa_unifies_tissue_aliases_across_files_and_rejects_same_sample_collision(self) -> None:
        mfa_a = self.exp_dir / "mfa" / "part_a"
        mfa_b = self.exp_dir / "mfa" / "part_b"
        samples_a = ["lfd_1", "lfd_2", "hfd_5", "hfd_6"]
        samples_b = ["lfd_3", "lfd_4", "hfd_7", "hfd_8"]
        write_mfa(_make_mfa_results(samples_a, {s: ["white_adipose", "renal_cortex"] for s in samples_a}, seed=1), mfa_a)
        write_mfa(_make_mfa_results(samples_b, {s: ["sq_adipose", "kidney"] for s in samples_b}, seed=2), mfa_b)

        _write_tsv(self.exp_dir / "samples.tsv", pd.DataFrame({
            "sample_id": samples_a + samples_b,
            "subject_id": ["m1", "m2", "m5", "m6", "m3", "m4", "m7", "m8"],
            "omics": ["mfa"] * 8,
            "tissue": ["multi_tissue"] * 8,
        }))
        self._base_exp_yaml(
            {
                "mfa": {
                    "enabled": True,
                    "input": {
                        "level": "table",
                        "format": "mfa_tables",
                        "files": [f"experiments/{EXP}/mfa/part_a", f"experiments/{EXP}/mfa/part_b"],
                    },
                    "tissue_aliases": {"renal_cortex": "kidney"},
                },
            },
            design_aliases={"white_adipose": "sq_adipose"},
        )

        run_dir, notes = _quiet(run_experiment, EXP, ["mfa"], root=self.root)
        mfa_dir = run_dir / "omics" / "mfa"
        self.assertEqual(sorted(p.name for p in (mfa_dir / "data").iterdir()), ["kidney", "sq_adipose"])
        self.assertFalse(any("white_adipose" in str(p) or "renal_cortex" in str(p) for p in mfa_dir.rglob("*")))
        diff = read_differential(mfa_dir / "tables" / "differential.parquet")
        self.assertEqual(set(diff["tissue"]), {"sq_adipose", "kidney"})
        self.assertEqual(set(diff["n_numerator"]), {4})
        self.assertEqual(set(diff["n_denominator"]), {4})
        self.assertTrue(any("white_adipose → sq_adipose" in n for n in notes), notes)

        # 同一样本同时含 sq_adipose 与 white_adipose 且映射到同一规范名时必须报错
        mfa_dup = self.exp_dir / "mfa" / "dup"
        write_mfa(_make_mfa_results(samples_a, {s: ["sq_adipose", "white_adipose"] for s in samples_a}, seed=3), mfa_dup)
        self._base_exp_yaml(
            {
                "mfa": {
                    "enabled": True,
                    "input": {"level": "table", "format": "mfa_tables", "files": [f"experiments/{EXP}/mfa/dup"]},
                },
            },
            design_aliases={"white_adipose": "sq_adipose"},
        )
        with self.assertRaisesRegex(RunError, "重复通量行"):
            _quiet(run_experiment, EXP, ["mfa"], root=self.root)


class MultiOmicsIntegrationAndReplicationAliasTest(unittest.TestCase):
    """全流程 pipeline + 全整合模式（含旧结果包）+ 跨队列复现的组织别名一致性。"""

    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.root = Path(self._tmp.name).resolve()
        _write_project(self.root)
        self.exp_dir = self.root / "experiments" / EXP
        self.exp_dir.mkdir(parents=True)

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def _build_full_experiment(self) -> None:
        subjects = pd.DataFrame({
            "subject_id": [f"m{i}" for i in range(1, 9)],
            "diet": ["LFD"] * 4 + ["HFD"] * 4,
            "treatment": ["Saline"] * 8,
        })
        _write_tsv(self.exp_dir / "subjects.tsv", subjects)

        # 转录组：sq_adipose + white_adipose（合并为 sq_adipose），另含 cerebral_cortex（归一为 brain）
        rng = np.random.default_rng(101)
        genes = ["Pdha1", "Pdhb", "Cs", "Hk2", "Pkm", "Pfkm", *[f"Gene{i:03d}" for i in range(54)]]
        tx_cols = [f"rna_sq_{i}" for i in (1, 2, 5, 6)] + [f"rna_wat_{i}" for i in (3, 4, 7, 8)] \
            + [f"rna_ctx_{i}" for i in range(1, 9)]
        tx_data = {}
        for col in tx_cols:
            is_hfd = int(col.rsplit("_", 1)[1]) >= 5
            base = rng.lognormal(5.0, 0.4, len(genes))
            mean = base * np.where(np.arange(len(genes)) < 20, 4.0 if is_hfd else 1.0, 1.0)
            tx_data[col] = rng.negative_binomial(6, 6 / (6 + mean))
        tx_path = self.exp_dir / "data" / "rna.tsv"
        tx_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(tx_data, index=pd.Index(genes, name="gene")).to_csv(tx_path, sep="\t")

        # 蛋白组：用 protein_matrix，含 white_adipose 与 brain
        prot_cols = [f"prot_wat_{i}" for i in range(1, 9)] + [f"prot_br_{i}" for i in range(1, 9)]
        prot_table = pd.DataFrame({"Gene.Name": genes, "Protein.Id": [f"P{i:05d}" for i in range(len(genes))]})
        for col in prot_cols:
            is_hfd = int(col.rsplit("_", 1)[1]) >= 5
            prot_table[col] = 20.0 + np.where(np.arange(len(genes)) < 20, 2.0 if is_hfd else 0.0, 0.0) + rng.normal(0, 0.2, len(genes))
        prot_path = self.exp_dir / "data" / "prot.tsv"
        prot_table.to_csv(prot_path, sep="\t", index=False)

        # 代谢组：用 peak_table，含 sq_adipose 与 renal_cortex（归一为 kidney）
        met_sq_cols = [f"met_sq_{i}" for i in range(1, 9)]
        met_rc_cols = [f"met_rc_{i}" for i in range(1, 9)]
        met_sq = peak_table(met_sq_cols, {f"met_sq_{i}" for i in range(5, 9)}, seed=31, n_features=24)
        met_rc = peak_table(met_rc_cols, {f"met_rc_{i}" for i in range(5, 9)}, seed=32, n_features=24)
        met_sq.loc[0, "Compounds"] = "Pyruvate"
        met_sq.loc[1, "Compounds"] = "Citrate"
        met_rc.loc[0, "Compounds"] = "Pyruvate"
        met_rc.loc[1, "Compounds"] = "Citrate"
        met_sq_path = self.exp_dir / "data" / "met_sq.tsv"
        met_rc_path = self.exp_dir / "data" / "met_rc.tsv"
        met_sq.to_csv(met_sq_path, sep="\t", index=False)
        met_rc.to_csv(met_rc_path, sep="\t", index=False)

        # MFA：含 white_adipose、brain、kidney
        mfa_samples = [f"lfd_{i}" for i in range(1, 5)] + [f"hfd_{i}" for i in range(5, 9)]
        mfa_path = self.exp_dir / "mfa" / "results"
        write_mfa(
            _make_mfa_results(mfa_samples, {s: ["white_adipose", "brain", "kidney"] for s in mfa_samples}, seed=41),
            mfa_path,
        )

        gmt = write_gmt({"GLYCOLYSIS_TCA": genes[:20], "BACKGROUND_SET": genes[20:40]}, self.root / "shared.gmt")

        sample_rows = []
        for idx, sid in zip((1, 2, 5, 6), [f"rna_sq_{i}" for i in (1, 2, 5, 6)]):
            sample_rows.append({"sample_id": sid, "subject_id": f"m{idx}", "omics": "transcriptomics", "tissue": "sq_adipose"})
        for idx, sid in zip((3, 4, 7, 8), [f"rna_wat_{i}" for i in (3, 4, 7, 8)]):
            sample_rows.append({"sample_id": sid, "subject_id": f"m{idx}", "omics": "transcriptomics", "tissue": "white_adipose"})
        for idx in range(1, 9):
            sample_rows.append({"sample_id": f"rna_ctx_{idx}", "subject_id": f"m{idx}", "omics": "transcriptomics", "tissue": "cerebral_cortex"})
            sample_rows.append({"sample_id": f"prot_wat_{idx}", "subject_id": f"m{idx}", "omics": "proteomics", "tissue": "white_adipose"})
            sample_rows.append({"sample_id": f"prot_br_{idx}", "subject_id": f"m{idx}", "omics": "proteomics", "tissue": "brain"})
            sample_rows.append({"sample_id": f"met_sq_{idx}", "subject_id": f"m{idx}", "omics": "metabolomics", "tissue": "sq_adipose"})
            sample_rows.append({"sample_id": f"met_rc_{idx}", "subject_id": f"m{idx}", "omics": "metabolomics", "tissue": "renal_cortex"})
        for idx, sid in enumerate(mfa_samples, start=1):
            sample_rows.append({"sample_id": sid, "subject_id": f"m{idx}", "omics": "mfa", "tissue": "multi_tissue"})
        _write_tsv(self.exp_dir / "samples.tsv", pd.DataFrame(sample_rows))

        # panels.yaml 故意写别名 white_adipose / cerebral_cortex，验证交付面板目录也统一到规范组织名
        (self.exp_dir / "panels.yaml").write_text(yaml.safe_dump({
            "panels": {
                "central_carbon": {
                    "title": "Central Carbon",
                    "tissues": ["white_adipose", "cerebral_cortex"],
                    "contrasts": [DIET],
                    "segments": [{"name": "TCA Entry", "reactions": ["PDH_m", "CS_m"]}],
                },
            },
        }, allow_unicode=True, sort_keys=False), encoding="utf-8")

        (self.exp_dir / "experiment.yaml").write_text(yaml.safe_dump({
            "schema_version": 1,
            "experiment": {"id": EXP},
            "design": {
                "factors": {"diet": {"level": "subject", "levels": ["LFD", "HFD"], "reference": "LFD"}},
                "contrasts": [{"id": DIET, "factor": "diet", "numerator": "HFD", "denominator": "LFD",
                               "role": "disease_effect"}],
                "stratify_by": ["tissue"],
                "tissue_aliases": {
                    "white_adipose": "sq_adipose",
                    "cerebral_cortex": "brain",
                    "renal_cortex": "kidney",
                },
            },
            "omics": {
                "transcriptomics": {
                    "enabled": True,
                    "input": {"level": "table", "format": "counts_matrix", "files": [f"experiments/{EXP}/data/rna.tsv"]},
                    "analysis": {"gene_sets": [str(gmt)], "gsea_permutations": 20},
                },
                "proteomics": {
                    "enabled": True,
                    "input": {"level": "table", "format": "protein_matrix", "files": [f"experiments/{EXP}/data/prot.tsv"],
                              "feature_column": "Gene.Name", "uniprot_column": "Protein.Id", "scale": "log2"},
                    "analysis": {"gene_sets": [str(gmt)], "gsea_permutations": 20},
                },
                "metabolomics": {
                    "enabled": True,
                    "input": {
                        "level": "table",
                        "format": "peak_table",
                        "files": [f"experiments/{EXP}/data/met_sq.tsv", f"experiments/{EXP}/data/met_rc.tsv"],
                        "tissue_files": {
                            "sq_adipose": f"experiments/{EXP}/data/met_sq.tsv",
                            "renal_cortex": f"experiments/{EXP}/data/met_rc.tsv",
                        },
                    },
                    "analysis": {"plsda_permutation_num": 10, "gsea_permutation_num": 20, "random_seed": 7,
                                 "derived_features": False},
                },
                "mfa": {
                    "enabled": True,
                    "input": {"level": "table", "format": "mfa_tables", "files": [f"experiments/{EXP}/mfa/results"]},
                },
            },
            "integration": {
                "enabled": True,
                "modes": ["flux_anchored", "pathway_convergence", "reaction_dossier"],
                "parameters": {
                    "reaction_dossier": {
                        "reactions": ["PDH_m", "CS_m"],
                        "panels": "panels.yaml",
                        "figures": {"formats": ["png"], "max_reactions": 2},
                        "pathway_map": {"enabled": True, "formats": ["png"]},
                    },
                },
            },
        }, allow_unicode=True, sort_keys=False), encoding="utf-8")

    def test_end_to_end_pipeline_has_zero_alias_folder_duplication(self) -> None:
        self._build_full_experiment()
        results = _quiet(run_pipeline, self.root, [EXP], log=lambda _m: None)
        self.assertEqual(len(results), 1)
        res = results[0]
        self.assertIsNotNone(res.run_dir)
        self.assertIsNotNone(res.integration_dir)
        self.assertIsNotNone(res.delivery_dir)

        # 整个 results 下不应出现任何别名目录名或文件名
        results_root = self.root / "results"
        bad_paths = [
            p.relative_to(results_root).as_posix()
            for p in results_root.rglob("*")
            if any(alias in p.name for alias in ("white_adipose", "cerebral_cortex", "renal_cortex"))
        ]
        self.assertEqual(bad_paths, [], f"不应输出未归一的别名文件或目录：{bad_paths}")

        # flux_anchored：sq_adipose 对齐全部四层，brain 与 kidney 也成功对齐
        fa_table = pd.read_parquet(res.integration_dir / "integration" / "flux_anchored" / "tables" / "flux_evidence_table.parquet")
        self.assertEqual(set(fa_table["tissue"]), {"sq_adipose", "brain", "kidney"})
        self.assertEqual(set(fa_table.loc[fa_table["tissue"] == "sq_adipose", "layer"]),
                         {"mfa", "transcriptomics", "proteomics", "metabolomics"})

        # pathway_convergence：转录组与蛋白组在 sq_adipose 与 brain 下成功融合
        pc_fusion = pd.read_parquet(res.integration_dir / "integration" / "pathway_convergence" / "tables" / "pathway_fusion.parquet")
        self.assertTrue({"sq_adipose", "brain"}.issubset(set(pc_fusion["tissue"])))

        # reaction_dossier 与交付目录：panels 与 pathway_map 均只输出规范组织文件夹
        self.assertEqual(sorted(p.name for p in (res.delivery_dir / "panels" / "central_carbon").iterdir()),
                         ["brain", "sq_adipose"])
        self.assertEqual(sorted(p.name for p in (res.delivery_dir / "pathway_map").iterdir()),
                         ["brain", "kidney", "sq_adipose"])
        # candidates 面板标注也自动将 panels.yaml 的 white_adipose / cerebral_cortex 归一化为 sq_adipose / brain
        cand_rxns = pd.read_csv(res.delivery_dir / "candidates" / "reactions_full.tsv", sep="\t", keep_default_na=False)
        labeled_tissues = set(cand_rxns.loc[cand_rxns["panels"].str.contains("central_carbon"), "tissue"])
        self.assertEqual(labeled_tissues, {"sq_adipose", "brain"})

    def test_legacy_unaliased_bundle_integrates_across_omics_and_guards_intra_omics_collision(self) -> None:
        """已有旧结果包（单组学各用不同组织名）在加上 tissue_aliases 后：
        - 跨组学异名可在 flux_anchored / pathway_convergence / reaction_dossier（含逐只动物样本读取）正常对齐；
        - 同一组学内出现两个原组织映射到同一规范名时，拒绝静默合并并提示重跑 omicsweave run。
        """
        _write_tsv(self.exp_dir / "subjects.tsv", pd.DataFrame({
            "subject_id": ["m1", "m2", "m3", "m4"], "diet": ["LFD", "LFD", "HFD", "HFD"],
        }))
        _write_tsv(self.exp_dir / "samples.tsv", pd.DataFrame({
            "sample_id": ["r1", "r2", "r3", "r4", "p1", "p2", "p3", "p4", "f1", "f2", "f3", "f4"],
            "subject_id": ["m1", "m2", "m3", "m4"] * 3,
            "omics": ["transcriptomics"] * 4 + ["proteomics"] * 4 + ["mfa"] * 4,
            "tissue": ["cerebral_cortex"] * 4 + ["brain"] * 4 + ["multi_tissue"] * 4,
        }))
        dummy_file = self.exp_dir / "data" / "dummy.tsv"
        dummy_file.parent.mkdir(parents=True, exist_ok=True)
        dummy_file.write_text("x\n1\n", encoding="utf-8")
        dummy_mfa_dir = self.exp_dir / "data" / "dummy_mfa"
        dummy_mfa_dir.mkdir(parents=True, exist_ok=True)

        def write_cfg(aliases: dict) -> None:
            (self.exp_dir / "experiment.yaml").write_text(yaml.safe_dump({
                "schema_version": 1,
                "experiment": {"id": EXP},
                "design": {
                    "factors": {"diet": {"level": "subject", "levels": ["LFD", "HFD"], "reference": "LFD"}},
                    "contrasts": [{"id": DIET, "factor": "diet", "numerator": "HFD", "denominator": "LFD",
                                   "role": "disease_effect"}],
                    "tissue_aliases": aliases,
                },
                "omics": {
                    "transcriptomics": {"enabled": True, "input": {"level": "table", "format": "counts_matrix",
                                                                   "files": [f"experiments/{EXP}/data/dummy.tsv"]}},
                    "proteomics": {"enabled": True, "input": {"level": "table", "format": "protein_matrix",
                                                              "files": [f"experiments/{EXP}/data/dummy.tsv"]}},
                    "mfa": {"enabled": True, "input": {"level": "table", "format": "mfa_tables",
                                                       "files": [f"experiments/{EXP}/data/dummy_mfa"]}},
                },
                "integration": {
                    "enabled": True,
                    "modes": ["flux_anchored", "pathway_convergence", "reaction_dossier"],
                    "parameters": {
                        "reaction_dossier": {"reactions": ["PDH_m"], "pathway_map": {"enabled": False}},
                    },
                },
            }, allow_unicode=True, sort_keys=False), encoding="utf-8")

        # 构造一个旧结果包：RNA 磁盘目录叫 cerebral_cortex、蛋白叫 brain、MFA 叫 brain
        src = self.root / "results" / EXP / "run_20261101_000000_old"
        for omics, tissue, feat, samples in (
            ("transcriptomics", "cerebral_cortex", "Pdha1", ["r1", "r2", "r3", "r4"]),
            ("proteomics", "brain", "Pdha1", ["p1", "p2", "p3", "p4"]),
            ("mfa", "brain", "PDH_m", ["f1", "f2", "f3", "f4"]),
        ):
            diff_row = {
                "contrast_id": [DIET], "contrast_role": ["disease_effect"], "tissue": [tissue], "omics": [omics],
                "feature_id": [feat], "feature_name": [feat], "log2fc": [1.8], "lfc_se": [0.2],
                "pvalue": [0.001], "padj": [0.005], "padj_method": ["fdr_bh"],
                "n_numerator": [2], "n_denominator": [2], "method": ["test"],
                "evidence_level": ["statistical_test"], "selected": [True],
                "mean_numerator": [10.0], "mean_denominator": [4.0],
            }
            if omics == "mfa":
                diff_row.update({
                    "effect_size": [2.0], "effect_size_type": ["hedges_g"],
                    "units": ["nmol/min/g"], "inference_unit": ["biological_replicate"],
                })
            write_differential(pd.DataFrame(diff_row), src / "omics" / omics / "tables" / "differential.parquet")
            mat = pd.DataFrame([[4.0, 4.2, 10.0, 10.2]], index=pd.Index([feat]), columns=samples)
            meta = DatasetMeta(
                omics=Omics(omics), tissue=tissue,
                value_type=ValueType.FLUX if omics == "mfa" else ValueType.INTENSITY,
                scale="linear" if omics == "mfa" else "log2",
                feature_id_type="mfa_reaction" if omics == "mfa" else "gene_symbol",
                units="nmol/min/g" if omics == "mfa" else None,
            )
            write_dataset(src / "omics" / omics / "data" / tissue, Dataset(mat, empty_features(mat.index), meta))
            if omics != "mfa":
                write_enrichment(pd.DataFrame({
                    "contrast_id": [DIET], "contrast_role": ["disease_effect"], "tissue": [tissue], "omics": [omics],
                    "method": ["ora"], "query": ["up"], "gene_set_library": ["KEGG"],
                    "term_id": ["ko00020"], "term_name": ["Citrate cycle (TCA cycle)"],
                    "pvalue": [0.001], "padj": [0.01], "padj_method": ["fdr_bh"], "effect_direction": ["up"],
                    "overlap_size": [5], "set_size": [20],
                }), src / "omics" / omics / "tables" / "enrichment.parquet")
        write_manifest(RunManifest(run_id=src.name, project_id="alias_proj", experiment_id=EXP,
                                   omics=[Omics.TRANSCRIPTOMICS, Omics.PROTEOMICS, Omics.MFA]), src)

        write_cfg({"cerebral_cortex": "brain"})
        int_dir, _ = _quiet(integrate_experiment, EXP, source=src.name, root=self.root)
        # 1) flux_anchored 统一到 brain，覆盖 RNA + 蛋白 + MFA
        fa_ev = pd.read_parquet(int_dir / "integration" / "flux_anchored" / "tables" / "flux_evidence_table.parquet")
        self.assertEqual(set(fa_ev["tissue"]), {"brain"})
        matched_layers = set(fa_ev.loc[fa_ev["n_matched"].fillna(0) > 0, "layer"])
        self.assertTrue({"mfa", "transcriptomics", "proteomics"}.issubset(matched_layers),
                        fa_ev[["layer", "assessment_status", "n_matched", "n_signal"]])
        assessed_layers = set(fa_ev.loc[fa_ev["assessment_status"].isin({"support", "opposite", "mixed", "no_signal", "unresolved_signal"}), "layer"])
        self.assertTrue({"mfa", "transcriptomics", "proteomics"}.issubset(assessed_layers),
                        fa_ev[["layer", "assessment_status", "n_matched"]])
        # 2) pathway_convergence 融合 RNA（原 cerebral_cortex）与蛋白（原 brain）
        pc_fusion = pd.read_parquet(int_dir / "integration" / "pathway_convergence" / "tables" / "pathway_fusion.parquet")
        self.assertEqual(set(pc_fusion["tissue"]), {"brain"})
        self.assertEqual(len(pc_fusion), 1)
        # 3) reaction_dossier 的逐只动物样本表成功回查原磁盘目录 data/cerebral_cortex
        samples_tbl = pd.read_parquet(int_dir / "integration" / "reaction_dossier" / "tables" / "dossier_samples.parquet")
        self.assertEqual(set(samples_tbl["tissue"]), {"brain"})
        self.assertEqual(set(samples_tbl["layer"]), {"mfa", "transcriptomics", "proteomics"})

        # 4) 组学内撞车守卫：若旧结果包中同一组学同时存在 cerebral_cortex 与 brain，拒绝静默合并
        tx_diff_path = src / "omics" / "transcriptomics" / "tables" / "differential.parquet"
        tx_old = read_differential(tx_diff_path)
        write_differential(pd.concat([tx_old, tx_old.assign(tissue="brain")], ignore_index=True), tx_diff_path)
        with self.assertRaisesRegex(RunError, "重新运行 omicsweave run"):
            _quiet(integrate_experiment, EXP, modes=["flux_anchored"], source=src.name, root=self.root)

    def test_replication_respects_experiment_tissue_aliases(self) -> None:
        # 队列 A 的结果包里写 white_adipose，且其 experiment.yaml 声明 white_adipose → sq_adipose；
        # 队列 B 的结果包里写 wat，且 replication source 声明 tissue_map: {wat: sq_adipose}
        for exp_id, tissue, aliases in (("cohort_a", "white_adipose", {"white_adipose": "sq_adipose"}),
                                        ("cohort_b", "wat", {})):
            edir = self.root / "experiments" / exp_id
            edir.mkdir(parents=True, exist_ok=True)
            _write_tsv(edir / "subjects.tsv", pd.DataFrame({"subject_id": ["m1"], "diet": ["LFD"]}))
            _write_tsv(edir / "samples.tsv", pd.DataFrame({
                "sample_id": ["s1"], "subject_id": ["m1"], "omics": ["transcriptomics"], "tissue": [tissue],
            }))
            (edir / "experiment.yaml").write_text(yaml.safe_dump({
                "schema_version": 1, "experiment": {"id": exp_id},
                "design": {"tissue_aliases": aliases},
            }, allow_unicode=True, sort_keys=False), encoding="utf-8")
            run_dir = self.root / "results" / exp_id / "run_20261101_000000_src"
            write_differential(pd.DataFrame({
                "contrast_id": [DIET] * 5, "contrast_role": ["disease_effect"] * 5, "tissue": [tissue] * 5,
                "omics": ["transcriptomics"] * 5, "feature_id": [f"G{i}" for i in range(5)],
                "log2fc": [2.0] * 5, "lfc_se": [0.2] * 5, "pvalue": [1e-4] * 5, "padj": [5e-4] * 5,
                "padj_method": ["fdr_bh"] * 5, "n_numerator": [4] * 5, "n_denominator": [4] * 5, "method": ["deseq2"] * 5,
            }), run_dir / "omics" / "transcriptomics" / "tables" / "differential.parquet")
            write_manifest(RunManifest(run_id=run_dir.name, project_id="alias_proj", experiment_id=exp_id,
                                       omics=[Omics.TRANSCRIPTOMICS]), run_dir)

        _write_project(self.root, replication=[{
            "id": "diet_rep",
            "sources": [
                {"experiment": "cohort_a", "contrast": DIET},
                {"experiment": "cohort_b", "contrast": DIET, "tissue_map": {"wat": "sq_adipose"}},
            ],
        }])
        rep_dir, _ = _quiet(replicate_project, self.root)
        rep = pd.read_parquet(rep_dir / "integration" / "replication" / "tables" / "replication.parquet")
        self.assertEqual(set(rep["tissue"]), {"sq_adipose"})
        self.assertTrue(rep["replicated"].all())

    def test_reversal_respects_experiment_tissue_aliases(self) -> None:
        _write_tsv(self.exp_dir / "subjects.tsv", pd.DataFrame({
            "subject_id": ["m1", "m2", "m3", "m4", "m5", "m6"],
            "group": ["LFD", "LFD", "HFD", "HFD", "HFD_TX", "HFD_TX"],
        }))
        _write_tsv(self.exp_dir / "samples.tsv", pd.DataFrame({
            "sample_id": [f"r{i}" for i in range(1, 7)],
            "subject_id": [f"m{i}" for i in range(1, 7)],
            "omics": ["transcriptomics"] * 6,
            "tissue": ["white_adipose"] * 6,
        }))
        dummy = self.exp_dir / "data" / "dummy.tsv"
        dummy.parent.mkdir(parents=True, exist_ok=True)
        dummy.write_text("x\n1\n", encoding="utf-8")
        (self.exp_dir / "experiment.yaml").write_text(yaml.safe_dump({
            "schema_version": 1,
            "experiment": {"id": EXP},
            "design": {
                "factors": {"group": {"level": "subject", "levels": ["LFD", "HFD", "HFD_TX"], "reference": "LFD"}},
                "contrasts": [
                    {"id": "HFD_vs_LFD", "factor": "group", "numerator": "HFD", "denominator": "LFD", "role": "disease_effect"},
                    {"id": "HFD_TX_vs_HFD", "factor": "group", "numerator": "HFD_TX", "denominator": "HFD", "role": "treatment_effect"},
                ],
                "reversal_pairs": [{"disease": "HFD_vs_LFD", "treatment": "HFD_TX_vs_HFD"}],
                "tissue_aliases": {"white_adipose": "sq_adipose"},
            },
            "omics": {
                "transcriptomics": {"enabled": True, "input": {"level": "table", "format": "counts_matrix",
                                                               "files": [f"experiments/{EXP}/data/dummy.tsv"]}},
            },
            "integration": {"enabled": True, "modes": ["reversal"]},
        }, allow_unicode=True, sort_keys=False), encoding="utf-8")

        src = self.root / "results" / EXP / "run_20261101_000000_rev"
        write_differential(pd.DataFrame({
            "contrast_id": ["HFD_vs_LFD", "HFD_TX_vs_HFD"],
            "contrast_role": ["disease_effect", "treatment_effect"],
            "tissue": ["white_adipose", "white_adipose"],
            "omics": ["transcriptomics", "transcriptomics"],
            "feature_id": ["Pdha1", "Pdha1"],
            "feature_name": ["Pdha1", "Pdha1"],
            "log2fc": [2.0, -2.0],
            "lfc_se": [0.2, 0.2],
            "pvalue": [1e-4, 1e-4],
            "padj": [1e-3, 1e-3],
            "padj_method": ["fdr_bh", "fdr_bh"],
            "n_numerator": [2, 2],
            "n_denominator": [2, 2],
            "method": ["deseq2", "deseq2"],
            "evidence_level": ["statistical_test", "statistical_test"],
        }), src / "omics" / "transcriptomics" / "tables" / "differential.parquet")
        write_manifest(RunManifest(run_id=src.name, project_id="alias_proj", experiment_id=EXP,
                                   omics=[Omics.TRANSCRIPTOMICS]), src)

        int_dir, _ = _quiet(integrate_experiment, EXP, modes=["reversal"], source=src.name, root=self.root)
        rev = pd.read_parquet(int_dir / "integration" / "reversal" / "tables" / "reversal.parquet")
        self.assertEqual(set(rev["tissue"]), {"sq_adipose"})
        self.assertEqual(list(rev["class"]), ["reversed"])


if __name__ == "__main__":
    unittest.main()
