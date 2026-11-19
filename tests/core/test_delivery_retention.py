#!/usr/bin/env python3
"""时间戳交付命名、latest.json 指针、omicsweave latest / clean、保留策略与入口开关测试。"""

from __future__ import annotations

import contextlib
import io
import json
import os
import re
import tempfile
import unittest
from datetime import datetime, timedelta, timezone
from pathlib import Path
from unittest import mock

import yaml

from omicsweave import cli
from omicsweave.contract.manifest import RunManifest, write_manifest
from omicsweave.core.config import ConfigError, load_experiment, load_project
from omicsweave.core.delivery import (
    DELIVERIES_DIR,
    INDEX_FILE,
    LATEST_POINTER_FILE,
    PIN_FILE,
    SUMMARY_FILE,
    Delivery,
    new_delivery_id,
    open_path,
    read_latest_pointer,
    resolve_latest_delivery,
    validate_delivery_tag,
)
from omicsweave.core.retention import RetentionPolicy, apply_cleanup, plan_cleanup
from omicsweave.core.scaffold import create_project
from tests.mfa.test_mfa_runner import EXP, build_two_factor_project

DELIVERY_ID_RE = re.compile(r"^delivery_\d{8}_\d{6}_[0-9a-f]{4}$")


def run_cli(argv: list[str]) -> tuple[int, str]:
    out = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(out):
        try:
            cli.main(argv)
            code = 0
        except SystemExit as exc:
            code = int(exc.code or 0)
    return code, out.getvalue()


def run_entry(project: Path, argv: list[str], *, retention: RetentionPolicy | None = None,
              allow_clean: bool = False) -> tuple[int, str]:
    out = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(out):
        try:
            cli.project_pipeline_main(project, retention=retention, allow_clean=allow_clean, argv=argv)
            code = 0
        except SystemExit as exc:
            code = int(exc.code or 0)
    return code, out.getvalue()


def _make_fake_delivery(
    results_root: Path,
    exp_id: str,
    delivery_id: str,
    *,
    status: str = "completed",
    created_at: str | None = None,
    completed_at: str | None = None,
    bundles: dict[str, dict] | None = None,
    pinned: bool = False,
    pin_file: bool = False,
    tag: str | None = None,
    mtime: float | None = None,
    with_summary: bool = True,
) -> Path:
    directory = results_root / DELIVERIES_DIR / exp_id / delivery_id
    directory.mkdir(parents=True, exist_ok=True)
    payload = {
        "delivery_id": delivery_id,
        "project": "test_proj",
        "experiment": exp_id,
        "created_at": created_at,
        "completed_at": completed_at,
        "tag": tag,
        "pinned": pinned,
        "status": status,
        "stages": {},
        "bundles": bundles or {},
        "artifacts": [],
    }
    index = directory / INDEX_FILE
    index.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    if with_summary and status == "completed":
        (directory / SUMMARY_FILE).write_bytes(b"fake-xlsx-content")
    if pin_file:
        (directory / PIN_FILE).write_text("pinned\n", encoding="utf-8")
    if mtime is not None:
        os.utime(index, (mtime, mtime))
        os.utime(directory, (mtime, mtime))
    return directory


def _make_fake_bundle(results_root: Path, exp_id: str, run_id: str, *,
                      source_run: str | None = None, pin_file: bool = False) -> Path:
    bundle_dir = results_root / exp_id / run_id
    bundle_dir.mkdir(parents=True, exist_ok=True)
    params = {"source_run": source_run} if source_run else {}
    write_manifest(RunManifest(run_id=run_id, project_id="test_proj", experiment_id=exp_id, parameters=params), bundle_dir)
    (bundle_dir / "payload.bin").write_bytes(b"x" * 2048)
    if pin_file:
        (bundle_dir / PIN_FILE).write_text("keep\n", encoding="utf-8")
    return bundle_dir


class DeliveryNamingAndPointerTest(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.root = Path(self._tmp.name).resolve()
        build_two_factor_project(self.root)
        self.edit(f"experiments/{EXP}/experiment.yaml",
                  lambda doc: doc.update(integration={"enabled": True, "modes": ["reversal"]}))

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def edit(self, relative: str, change) -> None:
        path = self.root / relative
        doc = yaml.safe_load(path.read_text(encoding="utf-8"))
        change(doc)
        path.write_text(yaml.safe_dump(doc, sort_keys=False, allow_unicode=True), encoding="utf-8")

    def test_timestamped_id_and_tag_validation(self) -> None:
        fixed = datetime(2026, 11, 19, 23, 30, 12)
        plain = new_delivery_id(now=fixed)
        self.assertRegex(plain, DELIVERY_ID_RE)
        self.assertTrue(plain.startswith("delivery_20261119_233012_"))

        tagged = new_delivery_id(now=fixed, tag="manuscript_v1.0")
        self.assertTrue(tagged.startswith("delivery_20261119_233012_"))
        self.assertTrue(tagged.endswith("_manuscript_v1.0"))
        self.assertEqual(validate_delivery_tag("  v2-draft  "), "v2-draft")

        for bad in ("", "   ", "../escape", "a/b", "has space", ".hidden", "a" * 65):
            with self.subTest(tag=bad), self.assertRaises(ValueError):
                validate_delivery_tag(bad)

        # 命令行传入非法标签应在运行前拦截，不创建任何结果目录
        code, output = run_cli(["--project-dir", str(self.root), "pipeline", "--tag", "../evil"])
        self.assertEqual(code, 1, output)
        self.assertIn("无效", output)
        self.assertFalse((self.root / "analysis_results").exists())

    def test_delivery_avoids_same_second_collision(self) -> None:
        experiment = load_experiment(load_project(self.root), EXP)
        fixed = datetime(2026, 11, 19, 23, 30, 12)
        d1 = Delivery(experiment, [], [], now=fixed)
        d2 = Delivery(experiment, [], [], now=fixed)
        self.assertNotEqual(d1.directory, d2.directory)
        self.assertRegex(d1.directory.name, DELIVERY_ID_RE)
        self.assertRegex(d2.directory.name, DELIVERY_ID_RE)

    def test_latest_pointer_updated_only_on_success_and_cli_latest(self) -> None:
        # 首次成功运行：生成时间戳目录与 latest.json 指针
        code, output = run_cli(["--project-dir", str(self.root), "pipeline", "--tag", "baseline", "--pin"])
        self.assertEqual(code, 0, output)
        exp_deliveries = self.root / "analysis_results" / DELIVERIES_DIR / EXP
        pointer_path = exp_deliveries / LATEST_POINTER_FILE
        self.assertTrue(pointer_path.is_file())
        self.assertIn(f"[{EXP}] 最新指针：{pointer_path}", output)

        pointer = read_latest_pointer(self.root / "analysis_results", EXP)
        self.assertIsNotNone(pointer)
        first_id = pointer["delivery_id"]
        self.assertTrue(first_id.startswith("delivery_"))
        self.assertTrue(first_id.endswith("_baseline"))
        # 指针内仅存跨平台 POSIX 相对路径
        self.assertEqual(pointer["path"], first_id)
        self.assertEqual(pointer["summary"], f"{first_id}/{SUMMARY_FILE}")
        self.assertEqual(pointer["index"], f"{first_id}/{INDEX_FILE}")
        self.assertTrue(pointer["pinned"])
        self.assertTrue((exp_deliveries / first_id / PIN_FILE).is_file())

        # 配置一个必败脚本重跑：新建失败交付，但不覆盖 latest.json
        (self.root / "analysis").mkdir(exist_ok=True)
        (self.root / "analysis" / "fail.py").write_text("import sys; sys.exit(2)\n", encoding="utf-8")
        self.edit("project.yaml", lambda d: d.update(pipeline={"scripts": [{"command": "analysis/fail.py"}]}))
        code, output = run_cli(["--project-dir", str(self.root), "pipeline", "--from", "scripts"])
        self.assertEqual(code, 1, output)
        after_fail = read_latest_pointer(self.root / "analysis_results", EXP)
        self.assertEqual(after_fail["delivery_id"], first_id)

        # CLI: latest 默认输出 / --path / --target / --json / --list / --open
        code, out_default = run_cli(["--project-dir", str(self.root), "latest"])
        self.assertEqual(code, 0, out_default)
        self.assertIn(first_id, out_default)
        self.assertIn("latest.json 指针", out_default)

        code, out_path = run_cli(["--project-dir", str(self.root), "latest", "--path", "--target", "summary"])
        self.assertEqual(code, 0, out_path)
        self.assertEqual(Path(out_path.strip()), (exp_deliveries / first_id / SUMMARY_FILE))

        code, out_json = run_cli(["--project-dir", str(self.root), "latest", "--json"])
        self.assertEqual(code, 0, out_json)
        parsed = json.loads(out_json)
        self.assertEqual(parsed["delivery_id"], first_id)
        self.assertEqual(parsed["source"], "pointer")
        self.assertTrue(parsed["pointer_exists"])

        code, out_list = run_cli(["--project-dir", str(self.root), "latest", "--list"])
        self.assertEqual(code, 0, out_list)
        self.assertIn("[latest pinned tag=baseline]", out_list)
        self.assertIn("| failed |", out_list)

        with mock.patch("omicsweave.core.delivery.open_path") as mocked_open:
            code, out_open = run_cli(["--project-dir", str(self.root), "latest", "--open"])
            self.assertEqual(code, 0, out_open)
            mocked_open.assert_called_once_with(exp_deliveries / first_id / SUMMARY_FILE)

    def test_latest_fallback_for_legacy_uuid_and_corrupt_pointer(self) -> None:
        results = self.root / "analysis_results"
        t1 = datetime(2026, 11, 18, 10, 0, 0, tzinfo=timezone.utc)
        t2 = datetime(2026, 11, 19, 12, 0, 0, tzinfo=timezone.utc)
        _make_fake_delivery(results, EXP, "767b49fac9fb4257a8848a2a2bd879bf", created_at="", mtime=t1.timestamp())
        _make_fake_delivery(results, EXP, "b6352c7cc6384b769165d6e3168121a6",
                            created_at=t2.isoformat(), mtime=t2.timestamp())

        # 无 latest.json 时自动回退扫描并选出最新的旧版 UUID 交付
        rec, source = resolve_latest_delivery(results, EXP)
        self.assertEqual(source, "scan")
        self.assertIsNotNone(rec)
        self.assertEqual(rec.delivery_id, "b6352c7cc6384b769165d6e3168121a6")

        # 指针含路径穿越或指向不存在目录时拒绝指针并安全回退
        pointer_path = results / DELIVERIES_DIR / EXP / LATEST_POINTER_FILE
        pointer_path.write_text(json.dumps({"path": "../outside"}), encoding="utf-8")
        rec2, source2 = resolve_latest_delivery(results, EXP)
        self.assertEqual(source2, "scan")
        self.assertEqual(rec2.delivery_id, "b6352c7cc6384b769165d6e3168121a6")

    def test_open_path_cross_platform_dispatch(self) -> None:
        target = self.root / "project.yaml"
        with mock.patch("omicsweave.core.delivery.os.name", "nt"), \
                mock.patch("omicsweave.core.delivery.os.startfile", create=True) as startfile:
            open_path(target)
            startfile.assert_called_once_with(str(target.resolve()))
        with mock.patch("omicsweave.core.delivery.os.name", "posix"), \
                mock.patch("omicsweave.core.delivery.sys.platform", "darwin"), \
                mock.patch("omicsweave.core.delivery.subprocess.run") as run_cmd:
            open_path(target)
            run_cmd.assert_called_once_with(["open", str(target.resolve())], check=True)
        with mock.patch("omicsweave.core.delivery.os.name", "posix"), \
                mock.patch("omicsweave.core.delivery.sys.platform", "linux"), \
                mock.patch("omicsweave.core.delivery.subprocess.run") as run_cmd:
            open_path(target)
            run_cmd.assert_called_once_with(["xdg-open", str(target.resolve())], check=True)


class RetentionAndCleanTest(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.root = Path(self._tmp.name).resolve()
        build_two_factor_project(self.root)
        self.results = self.root / "analysis_results"

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def test_policy_and_yaml_validation(self) -> None:
        with self.assertRaises(ValueError):
            RetentionPolicy(keep_last=0)
        with self.assertRaises(ValueError):
            RetentionPolicy(stale_running_hours=-1)

        # project.yaml 中默认关闭 retention
        project = load_project(self.root)
        self.assertFalse(project.config.pipeline.retention.enabled)
        self.assertEqual(project.config.pipeline.retention.keep_last, 5)

        # project.yaml 拼写错误或非法 keep_last 会被拦截
        path = self.root / "project.yaml"
        doc = yaml.safe_load(path.read_text(encoding="utf-8"))
        doc["pipeline"] = {"retention": {"enabled": True, "keep_last": 0}}
        path.write_text(yaml.safe_dump(doc, sort_keys=False, allow_unicode=True), encoding="utf-8")
        with self.assertRaises(ConfigError):
            load_project(self.root)

    def test_cleanup_preserves_latest_pinned_active_and_transitive_bundles(self) -> None:
        now = datetime(2026, 11, 20, 12, 0, 0, tzinfo=timezone.utc)
        # 构造多个 run_* 结果包：
        # run_old_orphan: 孤立，应被清理
        # run_pinned_bundle: 手动 .keep 锁定，应保留
        # run_src_kept: 被 kept integration bundle 传递引用，应保留
        # run_int_kept: 被最新交付引用，应保留
        # run_src_pinned_del: 被 pinned 旧交付引用，应保留
        _make_fake_bundle(self.results, EXP, "run_20261110_000001_aaaa")
        _make_fake_bundle(self.results, EXP, "run_20261110_000002_bbbb", pin_file=True)
        _make_fake_bundle(self.results, EXP, "run_20261111_000001_cccc")
        _make_fake_bundle(self.results, EXP, "run_20261119_000001_dddd")
        _make_fake_bundle(self.results, EXP, "run_20261119_000002_eeee", source_run="run_20261119_000001_dddd")
        # 非结果包目录（如 _redraws）绝不应被误删
        redraws = self.results / EXP / "_redraws" / "custom"
        redraws.mkdir(parents=True)
        (redraws / "map.png").write_bytes(b"img")
        # 独立实验（没有任何交付记录）：孤立清理绝不应清空它
        _make_fake_bundle(self.results, "standalone_exp", "run_20261118_000001_ffff")

        # 构造交付：1 个旧的 pinned、2 个旧的普通成功、1 个最新成功、1 个 failed、1 个超时 running、1 个刚启动 running
        _make_fake_delivery(
            self.results, EXP, "delivery_20261111_000001_1111_milestone",
            created_at=(now - timedelta(days=9)).isoformat(),
            bundles={"run": {"run_id": "run_20261111_000001_cccc"}},
            pin_file=True,
        )
        _make_fake_delivery(
            self.results, EXP, "delivery_20261112_000001_2222",
            created_at=(now - timedelta(days=8)).isoformat(),
            bundles={"run": {"run_id": "run_20261110_000001_aaaa"}},
        )
        _make_fake_delivery(
            self.results, EXP, "delivery_20261115_000001_3333",
            created_at=(now - timedelta(days=5)).isoformat(),
        )
        _make_fake_delivery(
            self.results, EXP, "delivery_20261119_000001_4444",
            created_at=(now - timedelta(days=1)).isoformat(),
            bundles={"integrate": {"run_id": "run_20261119_000002_eeee"}},
        )
        _make_fake_delivery(
            self.results, EXP, "delivery_20261118_000001_5555",
            status="failed", created_at=(now - timedelta(days=2)).isoformat(),
        )
        _make_fake_delivery(
            self.results, EXP, "delivery_20261117_000001_6666",
            status="running", created_at=(now - timedelta(hours=30)).isoformat(),
            mtime=(now - timedelta(hours=30)).timestamp(),
        )
        _make_fake_delivery(
            self.results, EXP, "delivery_20261120_113000_7777",
            status="running", created_at=(now - timedelta(minutes=30)).isoformat(),
            mtime=(now - timedelta(minutes=30)).timestamp(),
        )

        policy = RetentionPolicy(
            enabled=True, keep_last=1, drop_failed=True, prune_orphan_bundles=True, stale_running_hours=24.0,
        )
        plan = plan_cleanup(self.results, policy, [EXP, "standalone_exp"], now=now)
        planned_ids = {(item.kind, item.item_id) for item in plan.to_remove}
        self.assertEqual(planned_ids, {
            ("delivery", "delivery_20261112_000001_2222"),
            ("delivery", "delivery_20261115_000001_3333"),
            ("delivery", "delivery_20261118_000001_5555"),
            ("delivery", "delivery_20261117_000001_6666"),
            ("bundle", "run_20261110_000001_aaaa"),
        })
        self.assertTrue(any("standalone_exp" in note for note in plan.notes))

        # dry_run=True 不删除任何文件
        dry_report = apply_cleanup(self.results, policy, [EXP, "standalone_exp"], dry_run=True, now=now)
        self.assertEqual(dry_report.removed, [])
        self.assertTrue((self.results / DELIVERIES_DIR / EXP / "delivery_20261112_000001_2222").is_dir())

        # 执行清理
        report = apply_cleanup(self.results, policy, [EXP, "standalone_exp"], dry_run=False, now=now)
        self.assertEqual(len(report.removed), 5)
        self.assertEqual(report.failed, [])
        # 应保留的交付与结果包均健在
        for kept_del in (
            "delivery_20261111_000001_1111_milestone",
            "delivery_20261119_000001_4444",
            "delivery_20261120_113000_7777",
        ):
            self.assertTrue((self.results / DELIVERIES_DIR / EXP / kept_del).is_dir())
        for kept_bundle in (
            "run_20261110_000002_bbbb",
            "run_20261111_000001_cccc",
            "run_20261119_000001_dddd",
            "run_20261119_000002_eeee",
        ):
            self.assertTrue((self.results / EXP / kept_bundle).is_dir())
        self.assertTrue((redraws / "map.png").is_file())
        self.assertTrue((self.results / "standalone_exp" / "run_20261118_000001_ffff").is_dir())
        # latest.json 自动同步为最新成功交付
        pointer = read_latest_pointer(self.results, EXP)
        self.assertEqual(pointer["delivery_id"], "delivery_20261119_000001_4444")

    def test_clean_cli_defaults_to_preview_and_applies_only_with_flag(self) -> None:
        now = datetime.now(timezone.utc)
        _make_fake_delivery(self.results, EXP, "delivery_20261110_000001_1111",
                            created_at=(now - timedelta(days=3)).isoformat())
        _make_fake_delivery(self.results, EXP, "delivery_20261111_000001_2222",
                            created_at=(now - timedelta(days=2)).isoformat())
        _make_fake_delivery(self.results, EXP, "delivery_20261112_000001_3333",
                            created_at=(now - timedelta(days=1)).isoformat())

        # 不传 --apply 仅预览，不删文件
        code, output = run_cli(["--project-dir", str(self.root), "clean", "--keep-last", "1"])
        self.assertEqual(code, 0, output)
        self.assertIn("预览模式", output)
        self.assertTrue((self.results / DELIVERIES_DIR / EXP / "delivery_20261110_000001_1111").is_dir())

        # 传 --apply 后执行清理
        code, output = run_cli(["--project-dir", str(self.root), "clean", "--keep-last", "1", "--apply"])
        self.assertEqual(code, 0, output)
        self.assertIn("已执行清理", output)
        self.assertFalse((self.results / DELIVERIES_DIR / EXP / "delivery_20261110_000001_1111").exists())
        self.assertFalse((self.results / DELIVERIES_DIR / EXP / "delivery_20261111_000001_2222").exists())
        self.assertTrue((self.results / DELIVERIES_DIR / EXP / "delivery_20261112_000001_3333").is_dir())

    def test_locked_file_does_not_crash_cleanup(self) -> None:
        now = datetime.now(timezone.utc)
        d_old = _make_fake_delivery(self.results, EXP, "delivery_20261110_000001_1111",
                                    created_at=(now - timedelta(days=2)).isoformat())
        _make_fake_delivery(self.results, EXP, "delivery_20261112_000001_2222",
                            created_at=(now - timedelta(days=1)).isoformat())
        policy = RetentionPolicy(enabled=True, keep_last=1)
        with mock.patch("omicsweave.core.retention.shutil.rmtree", side_effect=PermissionError("Excel locked")):
            report = apply_cleanup(self.results, policy, [EXP], dry_run=False)
        self.assertEqual(len(report.removed), 0)
        self.assertEqual(len(report.failed), 1)
        self.assertEqual(report.failed[0][0].path, d_old)


class EntryFileSwitchTest(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.root = Path(self._tmp.name).resolve()
        build_two_factor_project(self.root)

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def test_entry_switches_disabled_by_default_and_work_when_enabled(self) -> None:
        disabled = RetentionPolicy(enabled=False, keep_last=1)
        # 默认关闭：连续运行 2 次 pipeline，两次交付都保留
        code1, out1 = run_entry(self.root, [], retention=disabled, allow_clean=False)
        self.assertEqual(code1, 0, out1)
        code2, out2 = run_entry(self.root, [], retention=disabled, allow_clean=False)
        self.assertEqual(code2, 0, out2)
        exp_dir = self.root / "analysis_results" / DELIVERIES_DIR / EXP
        deliveries = sorted(p.parent for p in exp_dir.glob("*/index.json"))
        self.assertEqual(len(deliveries), 2)

        # 默认关闭时，clean 命令被拦截且不删除文件；--help 和 latest 仍可用
        code, out = run_entry(self.root, ["clean", "--apply"], retention=disabled, allow_clean=False)
        self.assertEqual(code, 1, out)
        self.assertIn("CLEAN_COMMAND_ENABLED", out)
        self.assertEqual(len(list(exp_dir.glob("*/index.json"))), 2)

        code, out = run_entry(self.root, ["clean", "--help"], retention=disabled, allow_clean=False)
        self.assertEqual(code, 0, out)
        self.assertIn("--apply", out)

        code, out = run_entry(self.root, ["latest", "--path"], retention=disabled, allow_clean=False)
        self.assertEqual(code, 0, out)
        self.assertEqual(Path(out.strip()), deliveries[-1])

        # 在入口开启 CLEAN_COMMAND_ENABLED 后，clean 默认按入口参数预览，加 --apply 才删除
        enabled = RetentionPolicy(enabled=True, keep_last=1, drop_failed=True, prune_orphan_bundles=True)
        code, out = run_entry(self.root, ["clean"], retention=enabled, allow_clean=True)
        self.assertEqual(code, 0, out)
        self.assertIn("预览模式", out)
        self.assertEqual(len(list(exp_dir.glob("*/index.json"))), 2)

        code, out = run_entry(self.root, ["clean", "--apply"], retention=enabled, allow_clean=True)
        self.assertEqual(code, 0, out)
        self.assertEqual(len(list(exp_dir.glob("*/index.json"))), 1)
        # 孤立的旧 run_* 结果包也已联动清理
        self.assertEqual(len(list((self.root / "analysis_results" / EXP).glob("run_*/manifest.json"))), 1)

        # 开启 RETENTION_ENABLED 后，运行 pipeline 结束自动按保留策略清理旧交付
        code3, out3 = run_entry(self.root, [], retention=enabled, allow_clean=True)
        self.assertEqual(code3, 0, out3)
        self.assertIn("已按保留策略清理", out3)
        self.assertEqual(len(list(exp_dir.glob("*/index.json"))), 1)
        self.assertEqual(len(list((self.root / "analysis_results" / EXP).glob("run_*/manifest.json"))), 1)

    def test_scaffold_creates_entry_file_with_switches_off(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            proj = Path(tmp) / "fresh_proj"
            written = create_project(proj, "fresh_proj")
            entry = proj / "run_analysis.py"
            self.assertIn(entry, written)
            text = entry.read_text(encoding="utf-8")
            self.assertIn("RETENTION_ENABLED = False", text)
            self.assertIn("CLEAN_COMMAND_ENABLED = False", text)


if __name__ == "__main__":
    unittest.main()
