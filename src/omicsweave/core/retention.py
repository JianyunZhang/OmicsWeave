"""交付目录与关联结果包的保留策略和清理规划。

默认不在流程中自动启用；只有显式启用（例如在入口文件中开启开关或执行 ``omicsweave clean --apply``）
时才会删除文件。清理始终保护最新一次成功交付、带锁定的交付（``pinned: true`` 或 ``.keep``）以及
正在运行中的交付。
"""
from __future__ import annotations

import shutil
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal

from omicsweave.contract.manifest import MANIFEST_FILE, read_manifest
from .delivery import (
    INDEX_FILE,
    PIN_FILE,
    DeliveryRecord,
    delivery_experiments,
    list_deliveries,
    sync_latest_pointer,
)
from .model import RetentionSettings

REASON_LABELS = {
    "latest": "最新成功交付",
    "pinned": "已锁定（.keep / pinned）",
    "protected": "本次运行保护",
    "within_keep_last": "在保留数量内",
    "keep_all": "不限数量保留",
    "drop_failed_disabled": "保留失败交付",
    "active_running": "仍在运行或未满超时阈值",
    "referenced": "被保留交付引用",
    "exceeds_keep_last": "超出保留数量",
    "failed": "失败交付",
    "stale_running": "中断未完成交付",
    "corrupt": "索引损坏的交付",
    "orphan_bundle": "未被保留交付引用的孤立结果包",
}


@dataclass(frozen=True)
class RetentionPolicy:
    """交付与结果包清理规则；``enabled`` 默认关闭。"""
    enabled: bool = False
    keep_last: int | None = 5
    drop_failed: bool = True
    prune_orphan_bundles: bool = False
    stale_running_hours: float = 24.0

    def __post_init__(self) -> None:
        if self.keep_last is not None and self.keep_last < 1:
            raise ValueError(f"keep_last 至少为 1（当前：{self.keep_last}），不能删除最新一次成功交付")
        if self.stale_running_hours < 0:
            raise ValueError(f"stale_running_hours 不能为负数（当前：{self.stale_running_hours}）")

    @classmethod
    def from_settings(cls, settings: RetentionSettings) -> RetentionPolicy:
        return cls(
            enabled=settings.enabled,
            keep_last=settings.keep_last,
            drop_failed=settings.drop_failed,
            prune_orphan_bundles=settings.prune_orphan_bundles,
            stale_running_hours=settings.stale_running_hours,
        )

    def describe(self) -> str:
        keep_text = f"保留最近 {self.keep_last} 次成功交付" if self.keep_last is not None else "保留全部成功交付"
        failed_text = f"清理失败及超过 {self.stale_running_hours:g}h 的中断交付" if self.drop_failed else "保留失败交付"
        bundle_text = "联动清理孤立结果包" if self.prune_orphan_bundles else "不删除 run_* 结果包"
        return f"{keep_text}；{failed_text}；{bundle_text}"


@dataclass(frozen=True)
class CleanupItem:
    experiment: str
    kind: Literal["delivery", "bundle"]
    item_id: str
    path: Path
    reason: str
    size_bytes: int


@dataclass(frozen=True)
class KeptItem:
    experiment: str
    kind: Literal["delivery", "bundle"]
    item_id: str
    path: Path
    reason: str


@dataclass
class CleanupReport:
    results_root: Path
    policy: RetentionPolicy
    dry_run: bool
    experiments: list[str]
    kept: list[KeptItem] = field(default_factory=list)
    to_remove: list[CleanupItem] = field(default_factory=list)
    removed: list[CleanupItem] = field(default_factory=list)
    failed: list[tuple[CleanupItem, str]] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    @property
    def reclaimable_bytes(self) -> int:
        return sum(item.size_bytes for item in self.to_remove)

    @property
    def reclaimed_bytes(self) -> int:
        return sum(item.size_bytes for item in self.removed)

    @property
    def delivery_removals(self) -> list[CleanupItem]:
        return [item for item in self.to_remove if item.kind == "delivery"]

    @property
    def bundle_removals(self) -> list[CleanupItem]:
        return [item for item in self.to_remove if item.kind == "bundle"]


def directory_size(path: Path) -> int:
    """计算目录占用字节数；遇到瞬态访问错误时忽略该文件。"""
    total = 0
    try:
        for entry in path.rglob("*"):
            try:
                if entry.is_file() and not entry.is_symlink():
                    total += entry.stat().st_size
            except OSError:
                continue
    except OSError:
        return total
    return total


def format_bytes(num_bytes: int) -> str:
    size = float(max(num_bytes, 0))
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024.0 or unit == "TB":
            return f"{int(size)} {unit}" if unit == "B" else f"{size:.2f} {unit}"
        size /= 1024.0
    return f"{size:.2f} TB"


def _collect_protected_bundles(results_root: Path, experiment_id: str,
                               kept_deliveries: list[DeliveryRecord]) -> tuple[dict[str, Path], set[str], dict[str, str]]:
    """找出实验下全部有效结果包及受保护结果包集合（含交付直接引用、source_run 传递引用与 .keep 锁定）。"""
    exp_dir = results_root / experiment_id
    bundles: dict[str, Path] = {}
    if not exp_dir.is_dir():
        return bundles, set(), {}

    exp_resolved = exp_dir.resolve()
    for entry in sorted(exp_dir.iterdir()):
        if entry.is_dir() and not entry.name.startswith((".", "_")) and (entry / MANIFEST_FILE).is_file():
            bundles[entry.name] = entry

    protected: set[str] = set()
    reasons: dict[str, str] = {}

    def mark(run_id: str, reason: str) -> None:
        if run_id in bundles and run_id not in protected:
            protected.add(run_id)
            reasons[run_id] = reason

    for name, path in bundles.items():
        if (path / PIN_FILE).is_file():
            mark(name, "pinned")

    for delivery in kept_deliveries:
        for run_id in delivery.bundles.values():
            if run_id != delivery.delivery_id:
                mark(run_id, "referenced")
        for rel in delivery.bundle_paths.values():
            try:
                resolved = (delivery.directory / rel).resolve()
                if resolved.parent == exp_resolved and resolved.name in bundles and resolved.name != delivery.delivery_id:
                    mark(resolved.name, "referenced")
            except OSError:
                continue

    # 传递保护：若保留了整合/后处理结果包，它依赖的 source_run 与结果包内输入同样必须保留
    queue = list(protected)
    for delivery in kept_deliveries:
        if delivery.delivery_id in bundles and delivery.delivery_id not in queue:
            queue.append(delivery.delivery_id)
    while queue:
        current = queue.pop()
        bundle_dir = bundles[current]
        try:
            manifest = read_manifest(bundle_dir)
        except Exception:
            continue
        source_run = manifest.parameters.get("source_run")
        if isinstance(source_run, str) and source_run != current and source_run in bundles and source_run not in protected:
            mark(source_run, "referenced")
            queue.append(source_run)
        for item in manifest.files:
            if item.role != "input":
                continue
            parts = Path(item.path).parts
            for idx, part in enumerate(parts[:-1]):
                if part == experiment_id and idx + 1 < len(parts):
                    candidate = parts[idx + 1]
                    if candidate != current and candidate in bundles and candidate not in protected:
                        mark(candidate, "referenced")
                        queue.append(candidate)
    return bundles, protected, reasons


def plan_cleanup(
    results_root: str | Path,
    policy: RetentionPolicy,
    experiments: list[str] | None = None,
    *,
    protected_deliveries: dict[str, set[str]] | None = None,
    now: datetime | None = None,
) -> CleanupReport:
    """只读计算清理计划，不修改任何文件。"""
    root = Path(results_root).resolve()
    target_experiments = list(dict.fromkeys(experiments)) if experiments is not None else delivery_experiments(root)
    now_ts = (now or datetime.now(timezone.utc)).timestamp()
    report = CleanupReport(results_root=root, policy=policy, dry_run=True, experiments=target_experiments)
    protected_map = protected_deliveries or {}

    for exp_id in target_experiments:
        records = list_deliveries(root, exp_id)
        exp_protected = protected_map.get(exp_id, set())
        completed = [item for item in records if item.status == "completed"]
        latest_id = completed[-1].delivery_id if completed else None
        recent_ids = (
            {item.delivery_id for item in completed[-policy.keep_last:]}
            if policy.keep_last is not None
            else {item.delivery_id for item in completed}
        )

        kept_deliveries: list[DeliveryRecord] = []
        candidate_removals: list[tuple[DeliveryRecord, str]] = []
        for record in records:
            if record.delivery_id in exp_protected:
                kept_deliveries.append(record)
                report.kept.append(KeptItem(exp_id, "delivery", record.delivery_id, record.directory, "protected"))
                continue
            if record.status == "completed":
                if record.delivery_id == latest_id:
                    kept_deliveries.append(record)
                    report.kept.append(KeptItem(exp_id, "delivery", record.delivery_id, record.directory, "latest"))
                elif record.pinned:
                    kept_deliveries.append(record)
                    report.kept.append(KeptItem(exp_id, "delivery", record.delivery_id, record.directory, "pinned"))
                elif record.delivery_id in recent_ids:
                    reason = "within_keep_last" if policy.keep_last is not None else "keep_all"
                    kept_deliveries.append(record)
                    report.kept.append(KeptItem(exp_id, "delivery", record.delivery_id, record.directory, reason))
                else:
                    candidate_removals.append((record, "exceeds_keep_last"))
                continue

            if record.pinned:
                kept_deliveries.append(record)
                report.kept.append(KeptItem(exp_id, "delivery", record.delivery_id, record.directory, "pinned"))
                continue

            if not policy.drop_failed:
                kept_deliveries.append(record)
                report.kept.append(KeptItem(exp_id, "delivery", record.delivery_id, record.directory, "drop_failed_disabled"))
                continue

            if record.status == "failed":
                candidate_removals.append((record, "failed"))
            else:
                age_hours = max(0.0, (now_ts - record.last_activity_timestamp) / 3600.0)
                if age_hours >= policy.stale_running_hours:
                    reason = "corrupt" if record.status == "corrupt" else "stale_running"
                    candidate_removals.append((record, reason))
                else:
                    kept_deliveries.append(record)
                    report.kept.append(KeptItem(exp_id, "delivery", record.delivery_id, record.directory, "active_running"))

        bundles: dict[str, Path] = {}
        protected_bundles: set[str] = set()
        bundle_reasons: dict[str, str] = {}
        if kept_deliveries and (candidate_removals or policy.prune_orphan_bundles):
            bundles, protected_bundles, bundle_reasons = _collect_protected_bundles(root, exp_id, kept_deliveries)

        for record, reason in candidate_removals:
            if record.delivery_id in protected_bundles:
                kept_deliveries.append(record)
                report.kept.append(KeptItem(
                    exp_id, "delivery", record.delivery_id, record.directory,
                    bundle_reasons[record.delivery_id],
                ))
            else:
                report.to_remove.append(CleanupItem(
                    exp_id, "delivery", record.delivery_id, record.directory,
                    reason, directory_size(record.directory),
                ))

        if policy.prune_orphan_bundles:
            if not kept_deliveries:
                report.notes.append(f"[{exp_id}] 没有可保留的交付记录，跳过孤立结果包清理以防误删独立运行结果包")
                continue
            for run_id, bundle_path in bundles.items():
                if (bundle_path / INDEX_FILE).is_file():
                    continue
                if run_id in protected_bundles:
                    report.kept.append(KeptItem(exp_id, "bundle", run_id, bundle_path, bundle_reasons[run_id]))
                else:
                    report.to_remove.append(CleanupItem(
                        exp_id, "bundle", run_id, bundle_path,
                        "orphan_bundle", directory_size(bundle_path),
                    ))
    return report


def apply_cleanup(
    results_root: str | Path,
    policy: RetentionPolicy,
    experiments: list[str] | None = None,
    *,
    dry_run: bool = True,
    protected_deliveries: dict[str, set[str]] | None = None,
    now: datetime | None = None,
) -> CleanupReport:
    """执行（或预览）清理计划。``dry_run=True`` 时不删除任何文件。"""
    report = plan_cleanup(results_root, policy, experiments, protected_deliveries=protected_deliveries, now=now)
    if dry_run:
        return report

    report.dry_run = False
    root_resolved = report.results_root
    for item in report.to_remove:
        try:
            resolved = item.path.resolve()
            if resolved == root_resolved or not resolved.is_relative_to(root_resolved):
                report.failed.append((item, "目标路径不在结果根目录内，已拒绝删除"))
                continue
            shutil.rmtree(resolved)
            report.removed.append(item)
        except OSError as exc:
            report.failed.append((item, str(exc)))

    for exp_id in report.experiments:
        sync_latest_pointer(root_resolved, exp_id)
    return report


def format_cleanup_report(report: CleanupReport) -> str:
    """将清理计划或执行结果格式化为易读摘要。"""
    mode = "预览（未删除任何文件）" if report.dry_run else "已执行清理"
    lines = [
        f"结果根目录：{report.results_root}",
        f"清理模式：{mode}",
        f"策略：{report.policy.describe()}",
    ]
    if not report.experiments:
        lines.append("未找到任何交付目录，无需清理。")
        return "\n".join(lines)

    for exp_id in report.experiments:
        kept_deliveries = [k for k in report.kept if k.experiment == exp_id and k.kind == "delivery"]
        kept_bundles = [k for k in report.kept if k.experiment == exp_id and k.kind == "bundle"]
        exp_targets = [i for i in report.to_remove if i.experiment == exp_id]
        summary = f"[{exp_id}] 保留交付 {len(kept_deliveries)} 个"
        if report.policy.prune_orphan_bundles:
            summary += f"，保留结果包 {len(kept_bundles)} 个"
        summary += f"；{'待清理' if report.dry_run else '计划清理'} {len(exp_targets)} 项"
        lines.append(summary)
        for item in exp_targets:
            kind_label = "交付" if item.kind == "delivery" else "结果包"
            reason_label = REASON_LABELS.get(item.reason, item.reason)
            action = "将删除" if report.dry_run else ("已删除" if item in report.removed else "删除失败")
            lines.append(
                f"  - [{action}] {kind_label} {item.item_id}（{format_bytes(item.size_bytes)}，原因：{reason_label}）"
            )

    for note in report.notes:
        lines.append(f"[提示] {note}")
    for item, error in report.failed:
        lines.append(f"[警告] 无法删除 {item.path}：{error}")

    if report.dry_run:
        lines.append(
            f"合计：待清理 {len(report.to_remove)} 项（交付 {len(report.delivery_removals)}，"
            f"结果包 {len(report.bundle_removals)}），预计释放 {format_bytes(report.reclaimable_bytes)}"
        )
    else:
        lines.append(
            f"合计：已删除 {len(report.removed)} / {len(report.to_remove)} 项，"
            f"释放空间 {format_bytes(report.reclaimed_bytes)}"
        )
    return "\n".join(lines)
