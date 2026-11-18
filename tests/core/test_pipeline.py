#!/usr/bin/env python3
"""omicsweave pipeline：validate → run → integrate → 课题脚本一条命令完成，各步之间传递确切的运行 ID。"""

from __future__ import annotations

import contextlib
import io
import json
import os
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
        base = self.root / "analysis_results" / EXP
        found = sorted(p for p in base.iterdir()) if base.is_dir() else []
        return [p for p in found if kind is None or (p / kind).is_dir()]

    def records(self) -> list[dict]:
        return json.loads((self.root / "analysis" / f"record_{EXP}.json").read_text(encoding="utf-8"))

    def test_full_pipeline_passes_exact_run_ids(self) -> None:
        code, output = self.pipeline()
        self.assertEqual(code, 0, output)
        self.assertEqual(len(self.bundles()), 2)
        (run_dir,), (integration_dir,) = self.bundles("omics"), self.bundles("integration")
        self.assertTrue((run_dir / "omics" / "mfa" / "tables" / "differential.parquet").is_file())
        self.assertTrue((integration_dir / "integration" / "reversal" / "tables" / "reversal.parquet").is_file())
        # 整合读取的是本次 run 的结果包，脚本拿到的也是这两个运行 ID
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

        code, output = self.pipeline("--from", "scripts")
        self.assertEqual(code, 0, output)
        self.assertEqual(len(self.bundles()), 2)        # 没有新建结果包
        self.assertEqual(self.records()[-1]["integration"], first_integration.name)

        # 只给 --source 时从 integrate 开始：新建整合结果包，单组学结果包沿用
        code, output = self.pipeline("--source", run_dir.name)
        self.assertEqual(code, 0, output)
        self.assertEqual(self.bundles("omics"), [run_dir])
        (second_integration,) = set(self.bundles("integration")) - {first_integration}
        self.assertEqual(read_manifest(second_integration).parameters["source_run"], run_dir.name)
        self.assertEqual(self.records()[-1], {**self.records()[0], "integration": second_integration.name})

    def test_dry_run_lists_steps_without_writing(self) -> None:
        code, output = self.pipeline("--dry-run")
        self.assertEqual(code, 0, output)
        self.assertFalse((self.root / "analysis_results").exists())
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
        self.assertFalse((self.root / "analysis_results").exists())

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


if __name__ == "__main__":
    unittest.main()
