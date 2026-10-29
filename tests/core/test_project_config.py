#!/usr/bin/env python3
"""项目模型、配置加载与 validate 规则。以 examples/demo_project 为基线，逐项改坏后检查报错。"""

from __future__ import annotations

import shutil
import tempfile
import unittest
from pathlib import Path

import yaml

from omicsweave.core.config import ConfigError, load_experiment, load_project
from omicsweave.core.model import ContrastRole, ExperimentConfig, Omics
from omicsweave.core.scaffold import create_project
from omicsweave.core.validate import validate_project

DEMO = Path(__file__).resolve().parents[2] / "examples" / "demo_project"
EXP = "demo_glp1"


class DemoProjectCase(unittest.TestCase):
    """每个测试拿到一份 demo 项目的副本。"""

    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name) / "demo"
        shutil.copytree(DEMO, self.root)
        self.exp_dir = self.root / "experiments" / EXP

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def edit_yaml(self, mutate) -> None:
        path = self.exp_dir / "experiment.yaml"
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
        mutate(data)
        path.write_text(yaml.safe_dump(data, allow_unicode=True, sort_keys=False), encoding="utf-8")

    def edit_table(self, name: str, mutate) -> None:
        path = self.exp_dir / name
        lines = path.read_text(encoding="utf-8").splitlines()
        path.write_text("\n".join(mutate(lines)) + "\n", encoding="utf-8")

    def messages(self, level: str) -> list[str]:
        report = validate_project(self.root)
        return [i.message for i in report.issues if i.level == level]

    def assert_message(self, level: str, fragment: str) -> None:
        found = self.messages(level)
        self.assertTrue(any(fragment in m for m in found), f"未找到包含 {fragment!r} 的{level}：{found}")


class DemoBaselineTest(DemoProjectCase):
    def test_demo_project_is_clean(self) -> None:
        report = validate_project(self.root)
        self.assertEqual(report.issues, [])

    def test_loaded_model(self) -> None:
        exp = load_experiment(load_project(self.root), EXP)
        self.assertEqual(exp.config.omics.enabled(), [Omics.TRANSCRIPTOMICS, Omics.METABOLOMICS])
        roles = {c.id: c.role for c in exp.config.design.contrasts}
        self.assertEqual(roles["HFD_vs_LFD"], ContrastRole.DISEASE)
        self.assertEqual(len(exp.samples), 26)


class TableRulesTest(DemoProjectCase):
    def test_duplicate_subject(self) -> None:
        self.edit_table("subjects.tsv", lambda lines: lines + [lines[1]])
        self.assert_message("error", "subject_id 重复")

    def test_unknown_subject_in_samples(self) -> None:
        self.edit_table("samples.tsv", lambda lines: lines + ["rna_liver_x\tm99\ttranscriptomics\tliver\t4\tb1\tsample\t"])
        self.assert_message("error", "在 subjects.tsv 中找不到：m99")

    def test_qc_samples_may_have_no_subject_but_biological_must(self) -> None:
        self.assertNotIn("error", [i.level for i in validate_project(self.root).issues])
        self.edit_table("samples.tsv", lambda lines: [lines[0], lines[1].replace("\tm01\t", "\t\t")] + lines[2:])
        self.assert_message("error", "生物学样本的 subject_id 为空")

    def test_undeclared_factor_level(self) -> None:
        self.edit_table("subjects.tsv", lambda lines: [lines[0], lines[1].replace("LFD", "XFD")] + lines[2:])
        self.assert_message("error", "因素 diet 中有未声明的水平：XFD")

    def test_missing_factor_column(self) -> None:
        self.edit_table("subjects.tsv", lambda lines: [l.rsplit("\t", 1)[0] for l in lines])
        self.edit_yaml(lambda d: d["design"]["factors"].update(sex={"level": "subject"}))
        self.assert_message("error", "因素 sex 声明在 subject 层，但 subjects.tsv 中没有该列")


class ConfigRulesTest(DemoProjectCase):
    def test_unknown_key_is_rejected(self) -> None:
        self.edit_yaml(lambda d: d["design"].update(stratifyby=["tissue"]))
        self.assert_message("error", "design.stratifyby")

    def test_contrast_with_undeclared_factor(self) -> None:
        def mutate(d):
            d["design"]["contrasts"][0]["factor"] = "genotype"
        self.edit_yaml(mutate)
        self.assert_message("error", "引用了未声明的因素 'genotype'")

    def test_reversal_pair_role_mismatch(self) -> None:
        def mutate(d):
            d["design"]["contrasts"][0]["role"] = "other"
        self.edit_yaml(mutate)
        self.assert_message("error", "role 应为 disease_effect")

    def test_unsupported_input_format(self) -> None:
        def mutate(d):
            d["omics"]["mfa"] = {"enabled": True, "input": {"level": "raw", "format": "thermo_raw"}}
        self.edit_yaml(mutate)
        self.assert_message("error", "mfa 的 raw 层不支持格式")

    def test_numeric_levels_from_yaml_become_strings(self) -> None:
        config = ExperimentConfig.model_validate({
            "experiment": {"id": "x"},
            "design": {
                "factors": {"timepoint": {"level": "sample", "type": "numeric", "levels": [4, 14]}},
                "contrasts": [{"id": "t", "factor": "timepoint", "numerator": 14, "denominator": 4}],
            },
        })
        self.assertEqual(config.design.factors["timepoint"].levels, ["4", "14"])
        self.assertEqual(config.design.contrasts[0].numerator, "14")

    def test_reference_variable_substitution(self) -> None:
        project_file = self.root / "project.yaml"
        data = yaml.safe_load(project_file.read_text(encoding="utf-8"))
        data["project"]["references"] = {"counts": "experiments/demo_glp1/data/rna_counts.tsv"}
        project_file.write_text(yaml.safe_dump(data, allow_unicode=True), encoding="utf-8")
        self.edit_yaml(lambda d: d["omics"]["transcriptomics"]["input"].update(files=["${project.references.counts}"]))
        exp = load_experiment(load_project(self.root), EXP)
        self.assertEqual(exp.config.omics.transcriptomics.input.files, ["experiments/demo_glp1/data/rna_counts.tsv"])

        self.edit_yaml(lambda d: d["omics"]["transcriptomics"]["input"].update(files=["${project.references.nope}"]))
        with self.assertRaises(ConfigError):
            load_experiment(load_project(self.root), EXP)


class ReplicateRulesTest(DemoProjectCase):
    def _keep_subjects(self, keep: set[str]) -> None:
        def mutate(lines):
            return [lines[0]] + [l for l in lines[1:] if l.split("\t")[1] in keep or l.split("\t")[1] == ""]
        self.edit_table("samples.tsv", mutate)

    def test_too_few_replicates(self) -> None:
        # HFD-GLP1 只保留 m09
        self._keep_subjects({f"m{i:02d}" for i in range(1, 10)})
        self.assert_message("warning", "GLP1_vs_Saline_in_HFD 在 transcriptomics 中生物学重复不足 2")

    def test_rank_test_floor(self) -> None:
        # 3 vs 3：Mann-Whitney 最小双侧 p = 0.1
        keep = {"m01", "m02", "m03", "m05", "m06", "m07", "m09", "m10", "m11"}
        self._keep_subjects(keep)
        self.assert_message("warning", "3 vs 3，最小 p = 0.1")


class InputRulesTest(DemoProjectCase):
    def test_missing_table_file(self) -> None:
        (self.exp_dir / "data" / "rna_counts.tsv").unlink()
        self.assert_message("error", "transcriptomics 的输入文件不存在")

    def test_enabled_omics_without_samples(self) -> None:
        def mutate(d):
            d["omics"]["mfa"] = {"enabled": True, "input": {"level": "table", "format": "legacy_mfa_xlsx",
                                                            "files": ["experiments/demo_glp1/data/rna_counts.tsv"]}}
        self.edit_yaml(mutate)
        self.assert_message("error", "mfa 已启用，但 samples.tsv 中没有该组学的样本")

    def test_raw_wiff_requires_scan_and_paired_fastq_needs_two_files(self) -> None:
        raw = self.root / "raw"
        raw.mkdir()
        (raw / "a.wiff").write_bytes(b"")
        (raw / "r1.fastq.gz").write_bytes(b"")

        def mutate_yaml(d):
            d["omics"]["metabolomics"]["input"] = {"level": "raw", "format": "sciex_wiff"}
            d["omics"]["transcriptomics"]["input"] = {"level": "raw", "format": "fastq_gz", "layout": "paired"}
        self.edit_yaml(mutate_yaml)

        def mutate_samples(lines):
            out = [lines[0]]
            for line in lines[1:]:
                cells = line.split("\t")
                cells[-1] = "raw/a.wiff" if cells[2] == "metabolomics" else "raw/r1.fastq.gz"
                out.append("\t".join(cells))
            return out
        self.edit_table("samples.tsv", mutate_samples)
        self.assert_message("error", "缺少同名的 .wiff.scan")
        self.assert_message("warning", "为双端测序，但以下样本不是 2 个文件")

    def test_flux_anchored_requires_mfa(self) -> None:
        self.edit_yaml(lambda d: d["integration"].update(modes=["flux_anchored"]))
        self.assert_message("error", "flux_anchored 需要 MFA 数据")


class ScaffoldTest(unittest.TestCase):
    def test_new_project_validates_and_is_not_overwritten(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "proj"
            written = create_project(root, "proj")
            self.assertTrue((root / "project.yaml") in written)
            self.assertEqual(validate_project(root).issues, [])
            with self.assertRaises(FileExistsError):
                create_project(root, "proj")


if __name__ == "__main__":
    unittest.main()
