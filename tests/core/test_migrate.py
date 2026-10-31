#!/usr/bin/env python3
"""omicsweave migrate：用合成的旧版四组学配置检查迁移结果，并确认旧配置与数据文件不被改动。"""

from __future__ import annotations

import contextlib
import hashlib
import io
import tempfile
import unittest
from pathlib import Path

import pandas as pd
import yaml

from omicsweave import cli
from omicsweave.core.config import load_experiment, load_project
from omicsweave.core.migrate import (
    MFA_TISSUE,
    REPORT_FILE,
    MigrationError,
    migrate_experiment,
    write_migration,
)
from omicsweave.core.validate import validate_project

EXP = "legacy_exp"
ANIMALS = {"86": "HFD", "87": "HFD", "88": "LFD", "89": "LFD"}


def _yaml(path: Path, doc) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(doc, allow_unicode=True, sort_keys=False), encoding="utf-8")


def _excel(path: Path, sheets: dict[str, pd.DataFrame]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(path) as writer:
        for name, frame in sheets.items():
            frame.to_excel(writer, sheet_name=name, index=False)


def build_legacy_experiment(root: Path) -> None:
    """四个组学的旧版配置：转录组为旧版完整配置，MFA 为扁平布局，另有多组学整合设置。"""
    exp = root / "experiments" / EXP
    # 旧加载器只在文件存在时按项目目录解析相对路径，否则按配置文件所在目录
    for tool in ("tools/kallisto.exe", "references/tx.fa.gz", "references/tx.idx"):
        (root / tool).parent.mkdir(parents=True, exist_ok=True)
        (root / tool).write_bytes(b"")

    tx = exp / "transcriptomics_data"
    _yaml(tx / "config.yaml", {
        "experiment_id": EXP,
        "paths": {"raw_data_root": "raw/rnaseq",
                  "sample_list_tsv": f"experiments/{EXP}/transcriptomics_data/sample_list.tsv",
                  "output_dir": "analysis_results/${experiment_id}/transcriptomics_results"},
        "reference": {"kallisto_exe": "tools/kallisto.exe", "transcriptome_fasta": "references/tx.fa.gz",
                      "kallisto_index": "references/tx.idx"},
        "params": {"threads": 4, "tissues": ["liver"], "strandedness": "reverse"},
    })
    _yaml(tx / "groups.yaml", {
        "conditions": {"LFD": {"display": "LFD"}, "HFD": {"display": "HFD", "color": "#C44E52"}},
        "mouse_groups": {int(k): v for k, v in ANIMALS.items()},
        "sample_batches": {86: "b1", 87: "b2", 88: "b1", 89: "b2"},
    })
    (tx / "sample_list.tsv").write_text(
        "sample\tfile1\tfile2\n" + "".join(f"{a}_liver\t{a}_R1.fq.gz\t{a}_R2.fq.gz\n" for a in ANIMALS),
        encoding="utf-8")

    prot = exp / "proteomics_data"
    _yaml(prot / "config.yaml", {
        "experiment_id": EXP,
        "paths": {"output_dir": f"analysis_results/{EXP}/proteomics_results", "data_dirs": {"liver": "diann/liver"}},
        "params": {"tissues": ["liver"], "pval_cutoff": 0.01},
        "groups_file": "groups.yaml",
    })
    _yaml(prot / "groups.yaml", {"conditions": {"HFD": {}, "LFD": {}}, "mouse_groups": ANIMALS})
    report = prot / "diann" / "liver" / "report_liver.tsv"
    report.parent.mkdir(parents=True)
    report.write_text("Run\tProtein.Group\tPG.MaxLFQ\n" + "".join(
        f"MSF250153_zzk_23min_100ng_Liver_{a}\tP{i}\t1.0\n" for a in ANIMALS for i in range(2)), encoding="utf-8")

    met = exp / "metabolomics_data"
    _excel(met / "liver.xlsx", {"Sheet1": pd.DataFrame(
        columns=["Name", *(f"{a}_liver" for a in ANIMALS), "QC1", "Blank1"])})
    _yaml(met / "config.yaml", {
        "experiment_id": EXP,
        "paths": {"sample_data_liver": "liver.xlsx", "output_dir": f"analysis_results/{EXP}/metabolomics_results"},
        "groups_file": "groups.yaml",
    })
    _yaml(met / "groups.yaml", {
        "conditions": {"LFD": {"display": "LFD", "reference": True}, "HFD": {"display": "HFD"}},
        "mouse_groups": ANIMALS,
        "contrast_groups": {"HFD vs LFD": ["HFD", "LFD"]},
    })

    mfa = exp / "mfa_data"
    names = {f"{g.lower()}_{a}": g for a, g in ANIMALS.items()}
    _excel(mfa / "flux_raw_data.xlsx", {"sheet_information": pd.DataFrame({
        "sheet_name": [f"s{i}" for i in range(len(names))],
        "complete_result_name": list(names),
        "condition": list(names.values()),
    })})
    (mfa / "mid_raw_data.xlsx").write_bytes(b"")
    _yaml(mfa / "config.yaml", {
        "experiment_id": EXP,
        "paths": {"output_dir": "analysis_results/${experiment_id}/mfa_results",
                  "flux_file": f"experiments/{EXP}/mfa_data/flux_raw_data.xlsx",
                  "mid_file": f"experiments/{EXP}/mfa_data/mid_raw_data.xlsx"},
        "comparisons": [["hfd", "lfd"]],
        "analysis": {"differential_flux": {"fdr_threshold": 0.2}},
    })
    _yaml(mfa / "groups.yaml", {
        "conditions": {"lfd": {"color": "#4C72B0"}, "hfd": {"color": "#C44E52"}},
        "sample_groups": {n: g.lower() for n, g in names.items()},
        # 85 号只在血糖记录中出现，没有 MFA 样本
        "mouse_groups": {**{a: g.lower() for a, g in ANIMALS.items()}, "85": "hfd"},
    })

    _yaml(exp / "multiomics_integration_setting" / "config.yaml", {"experiment_id": EXP})


def _fingerprints(root: Path) -> dict[str, str]:
    return {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(root.rglob("*")) if p.is_file()}


def _quiet(func, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()):
        return func(*args, **kwargs)


class MigrateExperimentTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        cls.root = Path(cls._tmp.name).resolve()
        build_legacy_experiment(cls.root)
        cls.before = _fingerprints(cls.root)
        cls.migration = _quiet(migrate_experiment, EXP, cls.root)
        cls.written = write_migration(cls.migration, cls.root)
        cls.project = load_project(cls.root)
        cls.exp = load_experiment(cls.project, EXP)

    @classmethod
    def tearDownClass(cls) -> None:
        cls._tmp.cleanup()

    def samples(self, omics: str) -> pd.DataFrame:
        table = self.exp.samples
        return table[table["omics"] == omics].set_index("sample_id")

    def test_legacy_files_are_untouched(self) -> None:
        after = _fingerprints(self.root)
        new = {p.relative_to(self.root).as_posix() for p in self.written}
        self.assertEqual({k: v for k, v in after.items() if k not in new}, self.before)
        self.assertEqual(new, {"project.yaml", *(f"experiments/{EXP}/{n}" for n in
                                                 ("experiment.yaml", "subjects.tsv", "samples.tsv", REPORT_FILE))})

    def test_design(self) -> None:
        design = self.exp.config.design
        self.assertEqual(design.factors["group"].levels, ["LFD", "HFD"])
        self.assertEqual(design.factors["group"].level, "subject")
        # 四个组学的比较大小写、写法各不相同，合并为同一个
        self.assertEqual([(c.id, c.numerator, c.denominator, c.role.value) for c in design.contrasts],
                         [("HFD_vs_LFD", "HFD", "LFD", "other")])
        self.assertEqual(design.covariates, ["batch"])
        self.assertEqual(dict(zip(self.exp.subjects["subject_id"], self.exp.subjects["group"])), ANIMALS)

    def test_omics_settings(self) -> None:
        omics = self.exp.config.omics
        self.assertEqual({o.value for o in omics.enabled()}, {"transcriptomics", "proteomics", "metabolomics", "mfa"})
        tx = omics.transcriptomics
        self.assertEqual((tx.input.level, tx.input.format, tx.input.strandedness), ("raw", "fastq_gz", "reverse"))
        self.assertEqual(tx.adapter["executable"], "tools/kallisto.exe")
        self.assertEqual(omics.proteomics.input.files, [f"experiments/{EXP}/proteomics_data/diann/liver/report_liver.tsv"])
        self.assertEqual(omics.proteomics.analysis["pval_cutoff"], 0.01)
        self.assertEqual(omics.metabolomics.input.tissue_files, {"liver": f"experiments/{EXP}/metabolomics_data/liver.xlsx"})
        mfa = omics.mfa
        self.assertEqual((mfa.input.format, mfa.input.files), ("legacy_mfa_xlsx", [
            f"experiments/{EXP}/mfa_data/flux_raw_data.xlsx", f"experiments/{EXP}/mfa_data/mid_raw_data.xlsx"]))
        # 生效的分析参数 = 内置默认值 + 批次覆盖
        self.assertEqual(mfa.analysis["differential_flux"]["fdr_threshold"], 0.2)
        self.assertEqual(mfa.analysis["differential_flux"]["pvalue_threshold"], 0.05)
        self.assertEqual(mfa.adapter["batches"][0]["animal_groups"]["85"], "hfd")
        self.assertTrue(self.exp.config.integration.enabled)

    def test_transcriptomics_samples(self) -> None:
        tx = self.samples("transcriptomics")
        self.assertEqual(list(tx.index), [f"{a}_liver" for a in ANIMALS])
        self.assertEqual(tx.loc["86_liver", "files"], "raw/rnaseq/86_R1.fq.gz;raw/rnaseq/86_R2.fq.gz")
        self.assertEqual(dict(zip(tx["subject_id"], tx["batch"])), {"86": "b1", "87": "b2", "88": "b1", "89": "b2"})

    def test_proteomics_and_metabolomics_samples(self) -> None:
        prot = self.samples("proteomics")
        self.assertEqual(dict(zip(prot.index, prot["subject_id"])), {f"Liver{a}": a for a in ANIMALS})
        met = self.samples("metabolomics")
        # 与转录组同名的样本加上组学前缀
        self.assertEqual(set(met.index), {f"metabolomics_{a}_liver" for a in ANIMALS} | {"QC1", "Blank1"})
        self.assertEqual(met.loc["metabolomics_86_liver", "data_column"], "86_liver")
        self.assertEqual(met.loc["QC1", "data_column"], "")
        self.assertEqual(met.loc["QC1", "sample_type"], "qc")
        self.assertEqual(met.loc["QC1", "subject_id"], "")

    def test_mfa_samples_linked_to_animals(self) -> None:
        mfa = self.samples("mfa")
        self.assertEqual(dict(zip(mfa.index, mfa["subject_id"])), {f"{g.lower()}_{a}": a for a, g in ANIMALS.items()})
        self.assertEqual(set(mfa["tissue"]), {MFA_TISSUE})

    def test_report_lists_follow_ups(self) -> None:
        report = (self.root / "experiments" / EXP / REPORT_FILE).read_text(encoding="utf-8")
        for text in ("model_id", "role", "integration.modes", "推断", "组学前缀", "project.yaml"):
            self.assertIn(text, report)

    def test_validation_reports_only_expected_gaps(self) -> None:
        report = validate_project(self.root, [EXP])
        errors = [issue.message for issue in report.errors]
        self.assertTrue(errors)
        for message in errors:
            self.assertRegex(message, "model_id|flux_units|原始文件不存在|modes")


class MigrateLayoutsTest(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.root = Path(self._tmp.name).resolve()

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def test_mfa_protocol_batches(self) -> None:
        mfa = self.root / "experiments" / EXP / "mfa_data"
        _yaml(mfa / "experiment.yaml", {
            "defaults": {"active_protocol": "glc", "active_batch": "b2",
                         "output_root": "analysis_results/${experiment_id}/mfa_results"},
            "tracer_protocols": {"glc": {"enabled": True, "batches": ["b1", "b2"]},
                                 "gln": {"enabled": False, "batches": ["b3"]}},
        })
        for batch, animals in (("b1", ["1", "2"]), ("b2", ["3", "4"]), ("b3", ["5"])):
            _yaml(mfa / ("gln" if batch == "b3" else "glc") / batch / "config.yaml", {"experiment_id": EXP})
            _yaml(mfa / ("gln" if batch == "b3" else "glc") / batch / "groups.yaml", {
                "conditions": {"hfd": {}, "lfd": {}},
                "sample_groups": {a: "hfd" if int(a) % 2 else "lfd" for a in animals},
            })

        migration = _quiet(migrate_experiment, EXP, self.root)
        samples = migration.samples.set_index("sample_id")
        self.assertEqual(list(samples.index), ["1", "2", "3", "4"])
        self.assertEqual(list(samples["batch"]), ["b1", "b1", "b2", "b2"])
        self.assertEqual(set(samples["protocol"]), {"glc"})
        settings = migration.experiment["omics"]["mfa"]
        self.assertEqual([b["batch"] for b in settings["adapter"]["batches"]], ["b1", "b2"])
        self.assertEqual(settings["input"]["files"], [])
        # 未配置 comparisons：与旧流程相同，按条件两两组合
        self.assertEqual([c["id"] for c in migration.experiment["design"]["contrasts"]], ["hfd_vs_lfd"])
        for text in ("停用", "2 个批次", "sample_groups"):
            self.assertIn(text, migration.report)

    def test_existing_project_yaml_is_kept(self) -> None:
        build_legacy_experiment(self.root)
        (self.root / "project.yaml").write_text("schema_version: 1\nproject: {id: lab, species: mus_musculus}\n",
                                                encoding="utf-8")
        migration = _quiet(migrate_experiment, EXP, self.root)
        self.assertIsNone(migration.project_yaml)
        written = write_migration(migration, self.root)
        self.assertNotIn(self.root / "project.yaml", written)

    def test_refuses_to_overwrite_without_force(self) -> None:
        build_legacy_experiment(self.root)
        migration = _quiet(migrate_experiment, EXP, self.root)
        write_migration(migration, self.root)
        with self.assertRaises(FileExistsError):
            write_migration(migration, self.root)
        write_migration(migration, self.root, force=True)

    def test_no_legacy_config(self) -> None:
        (self.root / "experiments" / EXP).mkdir(parents=True)
        with self.assertRaisesRegex(MigrationError, "没有找到旧版组学配置"):
            migrate_experiment(EXP, self.root)
        with self.assertRaisesRegex(MigrationError, "找不到实验目录"):
            migrate_experiment("missing", self.root)


class MigrateCliTest(unittest.TestCase):
    def test_dry_run_then_write(self) -> None:
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp).resolve()
            build_legacy_experiment(root)
            exp_dir = root / "experiments" / EXP
            out = io.StringIO()
            with contextlib.redirect_stdout(out), self.assertRaises(SystemExit) as exit_:
                cli.main(["--project-dir", str(root), "migrate", "--experiment", EXP, "--dry-run"])
            self.assertEqual(exit_.exception.code, 0)
            self.assertIn("预览", out.getvalue())
            self.assertFalse((exp_dir / "experiment.yaml").exists())

            out = io.StringIO()
            with contextlib.redirect_stdout(out), self.assertRaises(SystemExit) as exit_:
                cli.main(["--project-dir", str(root), "migrate", "--experiment", EXP])
            self.assertEqual(exit_.exception.code, 0)
            self.assertTrue((exp_dir / "experiment.yaml").is_file())
            self.assertIn("迁移结果的校验", out.getvalue())

            err = io.StringIO()
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(err), \
                    self.assertRaises(SystemExit) as exit_:
                cli.main(["--project-dir", str(root), "migrate", "--experiment", EXP])
            self.assertEqual(exit_.exception.code, 1)
            self.assertIn("--force", err.getvalue())


if __name__ == "__main__":
    unittest.main()
