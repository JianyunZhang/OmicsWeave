#!/usr/bin/env python3
"""omicsweave pipeline：validate → run → integrate → 课题脚本一条命令完成，各步之间传递确切的运行 ID。"""

from __future__ import annotations

import contextlib
import io
import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import yaml

from omicsweave import cli
from omicsweave.contract.manifest import read_manifest
from omicsweave.core.config import ConfigError, load_project
from tests.mfa.test_mfa_runner import EXP, build_two_factor_project

# 课题脚本：把收到的参数与项目目录环境变量写到 analysis/record_<实验>.json
RECORD_SCRIPT = """\
import json, os, sys
from pathlib import Path
experiment, run, integration, run_dir = sys.argv[1:5]
out = Path("analysis") / f"record_{experiment}.json"
history = json.loads(out.read_text(encoding="utf-8")) if out.is_file() else []
history.append({"run": run, "integration": integration, "run_dir": run_dir,
                "project_dir": os.environ.get("OMICSWEAVE_PROJECT_DIR"), "cwd": os.getcwd()})
out.write_text(json.dumps(history), encoding="utf-8")
sys.exit(int(os.environ.get("RECORD_EXIT", "0")))
"""
RECORD = {"name": "record", "command": "analysis/record.py {experiment} {run} {integration} {run_dir}"}


def run_cli(argv: list[str]) -> tuple[int, str]:
    out = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(out):
        try:
            cli.main(argv)
            code = 0
        except SystemExit as exc:
            code = int(exc.code or 0)
    return code, out.getvalue()


class PipelineTest(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.root = Path(self._tmp.name).resolve()
        build_two_factor_project(self.root)
        (self.root / "analysis").mkdir()
        (self.root / "analysis" / "record.py").write_text(RECORD_SCRIPT, encoding="utf-8")
        self.edit("experiments/" + EXP + "/experiment.yaml",
                  lambda doc: doc.update(integration={"enabled": True, "modes": ["reversal"]}))
        self.set_pipeline({"scripts": [RECORD]})

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def edit(self, relative: str, change) -> None:
        path = self.root / relative
        doc = yaml.safe_load(path.read_text(encoding="utf-8"))
        change(doc)
        path.write_text(yaml.safe_dump(doc, sort_keys=False, allow_unicode=True), encoding="utf-8")

    def set_pipeline(self, pipeline: dict) -> None:
        self.edit("project.yaml", lambda doc: doc.update(pipeline=pipeline))

    def pipeline(self, *args: str) -> tuple[int, str]:
        return run_cli(["--project-dir", str(self.root), "pipeline", *args])

    def bundles(self, kind: str | None = None) -> list[Path]:
        """结果包；``kind`` 为 omics / integration 时只取单组学或整合结果包。"""
        base = self.root / "results" / EXP
        found = sorted(p for p in base.iterdir() if (p / "manifest.json").is_file()) if base.is_dir() else []
        return [p for p in found if kind is None or (p / kind).is_dir()]

    def records(self) -> list[dict]:
        return json.loads((self.root / "analysis" / f"record_{EXP}.json").read_text(encoding="utf-8"))

    def test_full_pipeline_passes_exact_run_ids(self) -> None:
        code, output = self.pipeline()
        self.assertEqual(code, 0, output)
        self.assertEqual(len(self.bundles()), 1)
        (run_dir,), (integration_dir,) = self.bundles("omics"), self.bundles("integration")
        self.assertEqual(run_dir, integration_dir)
        self.assertTrue((run_dir / "omics" / "mfa" / "tables" / "differential.parquet").is_file())
        self.assertTrue((run_dir / "integration" / "reversal" / "tables" / "reversal.parquet").is_file())
        self.assertTrue((run_dir / "summary.xlsx").is_file())
        self.assertTrue((run_dir / "index.json").is_file())
        self.assertTrue((self.root / "results" / "latest.json").is_file())
        # 整合读取的是本次 run 的结果包，脚本拿到的也是这个运行 ID
        self.assertEqual(read_manifest(integration_dir).parameters["source_run"], run_dir.name)
        self.assertEqual(self.records(), [{"run": run_dir.name, "integration": integration_dir.name,
                                           "run_dir": str(run_dir), "project_dir": str(self.root),
                                           "cwd": str(self.root)}])
        argv = ["omicsweave", "--project-dir", str(self.root), "pipeline"]
        self.assertEqual(read_manifest(run_dir).command, argv)
        self.assertIn(f"[{EXP}] 整合结果包：{integration_dir}", output)

    def test_resume_from_integrate_and_scripts(self) -> None:
        self.assertEqual(self.pipeline()[0], 0)
        (run_dir,), (first_integration,) = self.bundles("omics"), self.bundles("integration")
        self.assertEqual(run_dir, first_integration)

        code, output = self.pipeline("--from", "scripts")
        self.assertEqual(code, 0, output)
        self.assertEqual(len(self.bundles()), 1)        # 没有新建结果包（仅新增脚本交付目录）
        self.assertEqual(len(list((self.root / "results" / EXP).glob("*/index.json"))), 2)
        self.assertEqual(self.records()[-1]["integration"], first_integration.name)

        # 只给 --source 时从 integrate 开始：新建整合结果包，单组学结果包沿用
        code, output = self.pipeline("--source", run_dir.name)
        self.assertEqual(code, 0, output)
        self.assertEqual(self.bundles("omics"), [run_dir])
        (second_integration,) = set(self.bundles("integration")) - {first_integration}
        self.assertEqual(read_manifest(second_integration).parameters["source_run"], run_dir.name)
        self.assertEqual(self.records()[-1], {**self.records()[0], "integration": second_integration.name})

    def test_reused_bundles_must_match_current_inputs(self) -> None:
        panels_file = self.root / "experiments" / EXP / "panels.yaml"
        initial_panels = yaml.safe_dump({
            "panels": {"p": {"segments": [{"name": "s", "reactions": ["PDH_m"]}]}}
        }, sort_keys=False, allow_unicode=True)
        panels_file.write_text(initial_panels, encoding="utf-8")
        self.edit("experiments/" + EXP + "/experiment.yaml",
                  lambda doc: doc.update(integration={
                      "enabled": True,
                      "modes": ["reversal", "reaction_dossier"],
                      "parameters": {"reaction_dossier": {"panels": "panels.yaml", "pathway_map": {"enabled": False}}},
                  }))
        self.assertEqual(self.pipeline()[0], 0)
        (run_dir,) = self.bundles("omics")
        # 仅改 panels.yaml（整合输入）：--from scripts 判为过期，但仍可沿用同一目录内的单组学结果包从 integrate 续跑
        panels_file.write_text(yaml.safe_dump({
            "panels": {"p": {"segments": [{"name": "s", "reactions": ["CS_m"]}]}}
        }, sort_keys=False, allow_unicode=True), encoding="utf-8")
        code, output = self.pipeline("--from", "scripts", "--dry-run")
        self.assertEqual(code, 1, output)
        self.assertIn("panels.yaml 在整合之后改过", output)
        code, output = self.pipeline("--source", run_dir.name, "--dry-run")
        self.assertEqual(code, 0, output)
        panels_file.write_text(initial_panels, encoding="utf-8")

        experiment = f"experiments/{EXP}/experiment.yaml"
        # 只改 integration 一节（如打开参数）：整合结果包过期，单组学结果包仍可沿用
        self.edit(experiment, lambda doc: doc["integration"].update(parameters={"reversal": {}}))
        code, output = self.pipeline("--from", "scripts", "--dry-run")
        self.assertEqual(code, 1, output)
        self.assertIn(f"已过期：{experiment} 在整合之后改过", output)
        self.assertIn(f"--from integrate --source {run_dir.name}", output)
        code, output = self.pipeline("--source", run_dir.name, "--dry-run")
        self.assertEqual(code, 0, output)
        # 样本表或 integration 以外的设置改过：单组学结果包过期，须从 run 开始
        self.edit(experiment, lambda doc: doc["experiment"].update(description="changed"))
        with open(self.root / "experiments" / EXP / "subjects.tsv", "a", encoding="utf-8") as handle:
            handle.write("\n")
        code, output = self.pipeline("--source", run_dir.name, "--dry-run")
        self.assertEqual(code, 1, output)
        self.assertIn(f"单组学结果包 {run_dir.name} 已过期", output)
        self.assertIn(f"experiments/{EXP}/subjects.tsv", output)
        self.assertIn("experiment.yaml 的 experiment 一节", output)

    def test_dry_run_lists_steps_without_writing(self) -> None:
        code, output = self.pipeline("--dry-run")
        self.assertEqual(code, 0, output)
        self.assertFalse((self.root / "results").exists())
        self.assertIn(f"omicsweave run --experiment {EXP} --omics all", output)
        self.assertIn(f"omicsweave integrate --experiment {EXP} --source <run>", output)
        self.assertIn(f"python analysis/record.py {EXP} <run> <integration> <run_dir>", output)

    def test_configuration_errors_stop_before_running(self) -> None:
        cases = [
            ({"scripts": [{"command": "analysis/missing.py {run}"}]}, [], "找不到"),
            ({"experiments": ["nope"]}, [], "找不到实验 nope"),
            ({"scripts": [{**RECORD, "experiments": ["typo"]}]}, [], "本项目没有的实验：typo"),
            ({"scripts": [RECORD]}, ["--from", "run", "--source", "run_x"], "--source 只用于"),
        ]
        for pipeline, args, message in cases:
            with self.subTest(message):
                self.set_pipeline(pipeline)
                code, output = self.pipeline(*args)
                self.assertEqual(code, 1, output)
                self.assertIn(message, output)
        # 脚本要整合结果包，实验却没有启用整合
        self.set_pipeline({"scripts": [RECORD]})
        self.edit("experiments/" + EXP + "/experiment.yaml",
                  lambda doc: doc.update(integration={"enabled": False}))
        code, output = self.pipeline()
        self.assertEqual(code, 1, output)
        self.assertIn("没有启用整合", output)
        self.assertFalse((self.root / "results").exists())

    def test_unknown_placeholder_is_a_config_error(self) -> None:
        self.set_pipeline({"scripts": [{"command": "analysis/record.py {runid}"}]})
        with self.assertRaisesRegex(ConfigError, r"未知占位符 \{runid\}"):
            load_project(self.root)

    def test_failing_script_gives_resume_command(self) -> None:
        with mock.patch.dict(os.environ, RECORD_EXIT="3"):
            code, output = self.pipeline()
        self.assertEqual(code, 1, output)
        (run_dir,) = self.bundles("omics")
        self.assertIn("脚本 record 失败（退出码 3）", output)
        self.assertIn(f"--from scripts --source {run_dir.name}", output)
        index, = (self.root / "results" / EXP).glob("*/index.json")
        failure = json.loads(index.read_text(encoding="utf-8"))
        self.assertEqual(failure["stages"]["script_1:record"]["status"], "failed")
        self.assertEqual(failure["stages"]["summary"]["status"], "pending")
        self.assertTrue(failure["stages"]["script_1:record"]["script_sha256"])
        self.assertEqual(failure["status"], "failed")
        self.assertIn("--integration-run " + failure["bundles"]["integrate"]["run_id"], output)

    def test_delivery_in_custom_root_is_isolated_and_resumable(self) -> None:
        target = self.root / "中文 结果"
        script = self.root / "analysis" / "交付脚本.py"
        script.write_text("from pathlib import Path\nimport sys\n"
                          "Path(sys.argv[1], '表.txt').write_text(sys.argv[2], encoding='utf-8')\n", encoding="utf-8")
        self.set_pipeline({"scripts": [{"command": "analysis/交付脚本.py {delivery_dir} {run_dir}", "outputs": ["表.txt"]}]})
        code, output = self.pipeline("--results-dir", str(target))
        self.assertEqual(code, 0, output)
        first, = (target / EXP).glob("*/index.json")
        original = first.read_bytes()
        bundle = json.loads(original)["bundles"]
        source = bundle["run"]["run_id"]
        integration = bundle["integrate"]["run_id"]
        self.assertEqual(source, first.parent.name)
        self.assertEqual(integration, first.parent.name)
        args = ("--results-dir", str(target), "--from", "scripts", "--source", source, "--integration-run")
        self.assertEqual(self.pipeline(*args, "wrong-run")[0], 1)
        code, output = self.pipeline(*args, integration)
        self.assertEqual(code, 0, output)
        second, = set(first.parent.parent.glob("*/index.json")) - {first}
        data = json.loads(second.read_bytes())
        self.assertEqual(data["status"], "completed")
        self.assertEqual(data["stages"]["integrate"]["status"], "reused")
        self.assertEqual(first.read_bytes(), original)
        self.assertEqual((second.parent / "表.txt").read_text(encoding="utf-8"), str(target / EXP / source))
        self.assertIn("summary.xlsx", [row["path"] for row in data["artifacts"]])
        self.assertFalse((self.root / "results").exists())
        self.set_pipeline({"scripts": [{"command": "analysis/交付脚本.py {delivery_dir} {run_dir}",
                                        "outputs": ["遗漏.txt"]}]})
        code, output = self.pipeline(*args, integration)
        self.assertEqual(code, 1)
        self.assertIn("遗漏.txt", output)

    def test_dry_run_rejects_unregistered_mode(self) -> None:
        # demo 的两层组学有共同动物，样本设计本身满足 sample_level 的前提。
        demo = Path(__file__).resolve().parents[2] / "examples" / "demo_project"
        root = self.root / "paired"
        shutil.copytree(demo / "experiments", root / "experiments")
        shutil.copyfile(demo / "project.yaml", root / "project.yaml")
        path = root / "experiments/demo_glp1/experiment.yaml"
        config = yaml.safe_load(path.read_text(encoding="utf-8"))
        config["integration"]["modes"].append("sample_level")
        path.write_text(yaml.safe_dump(config, allow_unicode=True), encoding="utf-8")
        code, output = run_cli(["--project-dir", str(root), "pipeline", "--dry-run"])
        self.assertEqual(code, 1)
        self.assertIn("尚未接入", output)
        self.assertFalse((root / "results").exists())


if __name__ == "__main__":
    unittest.main()
