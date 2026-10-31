#!/usr/bin/env python3
"""项目目录的上下文隔离：临时切换后自动恢复，线程之间互不覆盖，命令行运行后不残留设置。"""

from __future__ import annotations

import contextlib
import contextvars
import io
import tempfile
import threading
import unittest
from pathlib import Path

from omicsweave import cli
from omicsweave.core.project import project_context, project_dir, reset_project_dir, set_project_dir

DEMO = Path(__file__).resolve().parents[2] / "examples" / "demo_project"


def default_project_dir() -> Path:
    """未显式指定时的项目目录：在空上下文中求值。"""
    return contextvars.Context().run(project_dir)


class ProjectContextTest(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.a = Path(self._tmp.name, "a").resolve()
        self.b = Path(self._tmp.name, "b").resolve()

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def test_nested_contexts_restore_previous_value(self) -> None:
        before = project_dir()
        with project_context(self.a) as active:
            self.assertEqual(active, self.a)
            with project_context(self.b):
                self.assertEqual(project_dir(), self.b)
            self.assertEqual(project_dir(), self.a)
            with project_context(None):  # None 表示沿用当前设置
                self.assertEqual(project_dir(), self.a)
        self.assertEqual(project_dir(), before)

    def test_set_and_reset_with_token(self) -> None:
        before = project_dir()
        token = set_project_dir(self.a)
        try:
            self.assertEqual(project_dir(), self.a)
        finally:
            reset_project_dir(token)
        self.assertEqual(project_dir(), before)

    def test_threads_do_not_see_each_others_setting(self) -> None:
        barrier = threading.Barrier(2)
        seen: dict[str, Path] = {}

        def worker(name: str, path: Path) -> None:
            with project_context(path):
                barrier.wait()  # 两个线程都已设置好之后再读取
                seen[name] = project_dir()

        threads = [threading.Thread(target=worker, args=("a", self.a)),
                   threading.Thread(target=worker, args=("b", self.b))]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()
        self.assertEqual(seen, {"a": self.a, "b": self.b})

    def test_new_thread_does_not_inherit_callers_setting(self) -> None:
        seen = []
        with project_context(self.a):
            thread = threading.Thread(target=lambda: seen.append(project_dir()))
            thread.start()
            thread.join()
        self.assertEqual(seen, [default_project_dir()])

    def test_cli_does_not_leave_project_dir_behind(self) -> None:
        before = project_dir()
        with contextlib.redirect_stdout(io.StringIO()), self.assertRaises(SystemExit):
            cli.main(["--project-dir", str(DEMO), "validate"])
        self.assertEqual(project_dir(), before)


if __name__ == "__main__":
    unittest.main()
