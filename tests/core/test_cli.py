#!/usr/bin/env python3
"""命令行：init 与 validate 的退出码和输出。"""

from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
from pathlib import Path

from omicsweave import cli

DEMO = Path(__file__).resolve().parents[2] / "examples" / "demo_project"


def run(argv: list[str]) -> tuple[int, str]:
    out = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(out):
        try:
            cli.main(argv)
            code = 0
        except SystemExit as exc:
            code = int(exc.code or 0)
    return code, out.getvalue()


class CliTest(unittest.TestCase):
    def test_validate_demo_project(self) -> None:
        code, output = run(["--project-dir", str(DEMO), "validate"])
        self.assertEqual(code, 0, output)
        self.assertIn("共 0 个错误，0 个警告", output)

    def test_validate_reports_errors_with_exit_code(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            code, output = run(["--project-dir", tmp, "validate"])
        self.assertEqual(code, 1)
        self.assertIn("找不到", output)

    def test_init_then_validate(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "new_project"
            code, output = run(["--project-dir", str(root), "init", "--id", "p1", "--species", "homo_sapiens"])
            self.assertEqual(code, 0, output)
            self.assertIn("species: homo_sapiens", (root / "project.yaml").read_text(encoding="utf-8"))
            code, output = run(["--project-dir", str(root), "validate"])
            self.assertEqual(code, 0, output)
            code, output = run(["--project-dir", str(root), "init"])
            self.assertEqual(code, 1)
            self.assertIn("已存在", output)

    def test_unknown_argument_for_builtin_command(self) -> None:
        code, _ = run(["--project-dir", str(DEMO), "validate", "--bogus"])
        self.assertEqual(code, 2)


if __name__ == "__main__":
    unittest.main()
