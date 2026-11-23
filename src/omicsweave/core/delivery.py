"""一次 pipeline 调用的交付清单、latest.json 指针与 Excel 摘要。结果包通过相对路径引用。"""
from __future__ import annotations

import json
import os
import re
import secrets
import shutil
import subprocess
import sys
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Iterator
from urllib.parse import quote

import pandas as pd
import pyarrow.parquet as pq

from omicsweave.contract.manifest import MANIFEST_FILE, read_manifest
from .config import EXPERIMENT_FILE, PROJECT_FILE, Experiment
from .provenance import sha256_file

DELIVERIES_DIR = "_deliveries"
INDEX_FILE = "index.json"
SUMMARY_FILE = "summary.xlsx"
LATEST_POINTER_FILE = "latest.json"
LATEST_SCHEMA_VERSION = 1
PIN_FILE = ".keep"
DELIVERY_TAG_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$")

# 交付目录中反应档案面板产出的子目录，与无面板候选表的子目录
PANELS_DIR = "panels"
CANDIDATES_DIR = "candidates"
# 分组织代谢通路总图（图 E，G18）的副本
MAPS_DIR = "pathway_map"


def validate_delivery_tag(tag: str) -> str:
    """校验交付标签：仅允许字母、数字及 ``_`` / ``.`` / ``-``，长度 1–64，避免跨平台路径问题。"""
    cleaned = tag.strip()
    if not cleaned or not DELIVERY_TAG_PATTERN.match(cleaned) or cleaned in {".", ".."}:
        raise ValueError(
            f"交付标签 {tag!r} 无效：请以字母或数字开头，仅使用字母、数字、_、.、-，且不超过 64 个字符"
        )
    return cleaned


def new_delivery_id(now: datetime | None = None, tag: str | None = None) -> str:
    """形如 ``delivery_20261119_233012_a3f9``（带标签时附 ``_<tag>``）的交付 ID：按时间排序，末尾随机后缀避免同秒冲突。"""
    current = now or datetime.now()
    base = f"delivery_{current:%Y%m%d_%H%M%S}_{secrets.token_hex(2)}"
    if tag is not None:
        return f"{base}_{validate_delivery_tag(tag)}"
    return base


class Delivery:
    def __init__(self, experiment: Experiment, scripts: list[str], command: list[str], *,
                 tag: str | None = None, pinned: bool = False, now: datetime | None = None):
        project = experiment.project
        self.project_root = project.root
        validated_tag = validate_delivery_tag(tag) if tag is not None else None
        parent = project.path("results") / DELIVERIES_DIR / experiment.id
        directory: Path | None = None
        for _ in range(16):
            candidate = parent / new_delivery_id(now=now, tag=validated_tag)
            try:
                candidate.mkdir(parents=True, exist_ok=False)
                directory = candidate
                break
            except FileExistsError:
                continue
        if directory is None:
            raise FileExistsError(f"无法在 {parent} 创建唯一的交付目录")
        self.directory = directory
        if pinned:
            (self.directory / PIN_FILE).write_text("pinned\n", encoding="utf-8")
        self.stages = {"validate": {"status": "completed"}}
        self.stages.update({key: {"status": "pending"} for key in ["run", "integrate", *scripts, "summary"]})
        self.bundles: dict[str, dict] = {}
        self.artifacts: list[dict] = []
        self.data = {
            "delivery_id": self.directory.name, "project": project.id, "experiment": experiment.id,
            "created_at": datetime.now(timezone.utc).isoformat(), "tag": validated_tag, "pinned": bool(pinned),
            "command": command, "status": "running",
            "stages": self.stages, "bundles": self.bundles,
            "configuration": {str(path): sha256_file(path) for path in
                              (project.root / PROJECT_FILE, experiment.root / EXPERIMENT_FILE)},
            "artifacts": self.artifacts,
        }
        self.save()

    def save(self) -> None:
        index = self.directory / INDEX_FILE
        pending = self.directory / ".index.tmp"
        pending.write_text(json.dumps(self.data, ensure_ascii=False, indent=2), encoding="utf-8")
        pending.replace(index)

    def complete(self) -> Path:
        """标记交付完成并原子更新同实验下的 ``latest.json`` 指针。"""
        self.data["status"] = "completed"
        self.data["completed_at"] = datetime.now(timezone.utc).isoformat()
        self.save()
        return self.update_latest_pointer()

    def update_latest_pointer(self) -> Path:
        """写出指向本次交付的跨平台相对路径 JSON 指针（``_deliveries/<实验>/latest.json``）。"""
        delivery_id = self.directory.name
        summary = PurePosixPath(delivery_id, SUMMARY_FILE).as_posix() if (self.directory / SUMMARY_FILE).is_file() else None
        payload = {
            "schema_version": LATEST_SCHEMA_VERSION,
            "project": self.data.get("project"),
            "experiment": self.data.get("experiment"),
            "delivery_id": delivery_id,
            "path": delivery_id,
            "summary": summary,
            "index": PurePosixPath(delivery_id, INDEX_FILE).as_posix(),
            "tag": self.data.get("tag"),
            "pinned": bool(self.data.get("pinned")) or (self.directory / PIN_FILE).is_file(),
            "created_at": self.data.get("created_at"),
            "completed_at": self.data.get("completed_at"),
            "bundles": {stage: info["run_id"] for stage, info in self.bundles.items() if isinstance(info, dict) and info.get("run_id")},
        }
        return write_latest_pointer(self.directory.parent, payload)

    @contextmanager
    def stage(self, name: str, success: str = "completed", **details) -> Iterator[None]:
        record = self.stages[name]
        record.update(status="running", **details)
        self.save()
        try:
            yield
        except BaseException as exc:
            record.update(status="failed", error=str(exc) or type(exc).__name__)
            self.data["status"] = "failed"
            try:
                self.collect_exports()
            except OSError as inventory_error:
                record["inventory_error"] = str(inventory_error)
            self.save()
            raise
        record["status"] = success
        self.save()

    def bundle(self, stage: str, directory: Path) -> None:
        manifest = read_manifest(directory)
        base = Path(os.path.relpath(directory, self.directory))
        self.bundles[stage] = {"run_id": manifest.run_id, "path": base.as_posix(),
                               "config_sha256": manifest.config_sha256, "created_at": manifest.created_at.isoformat(),
                               "manifest_sha256": sha256_file(directory / MANIFEST_FILE)}
        for item in manifest.files:
            if item.role != "output":
                continue
            record = {"stage": stage, "path": (base / item.path).as_posix(),
                      "size": item.size, "sha256": item.sha256}
            if not (directory / item.path).is_file():
                record["missing"] = True
            self.artifacts.append(record)

    def collect_exports(self) -> None:
        self.artifacts[:] = [item for item in self.artifacts if item["stage"] != "delivery"]
        ignored = {self.directory / INDEX_FILE, self.directory / ".index.tmp", self.directory / PIN_FILE}
        for path in sorted(self.directory.rglob("*")):
            if path.is_file() and path not in ignored:
                self.artifacts.append({"stage": "delivery", "path": path.relative_to(self.directory).as_posix(),
                                       "size": path.stat().st_size, "sha256": sha256_file(path)})

    def _panel_rows(self) -> list[dict]:
        """整合结果包中反应档案面板的产出：矩阵（有 XLSX 时只取 XLSX）、图件（有 PDF 时只取 PDF）、节点分流与检验表。"""
        paths = {PurePosixPath(item["path"]) for item in self.artifacts
                 if item["stage"] == "integrate" and not item.get("missing")}
        rows = []
        for path in sorted(paths):
            if path.parts[-6:-3] != ("integration", "reaction_dossier", "figures"):
                continue
            name = path.name
            if name.startswith("matrix_") and path.suffix == ".xlsx":
                kind = "面板矩阵与成员"
            elif name.startswith("matrix_") and path.suffix == ".tsv":
                if path.with_name(path.stem.removesuffix("_members") + ".xlsx") in paths:
                    continue
                kind = "成员变化" if path.stem.endswith("_members") else "面板矩阵"
            elif path.suffix == ".tsv" and name.startswith("node_split_"):
                kind = "节点分流"
            elif path.suffix == ".tsv" and name.startswith("pathway_tests_"):
                kind = "模块与交叉点检验"
            elif name.startswith(("matrix_", "node_split_", "tracer_coverage_")) and (
                    path.suffix == ".pdf" or path.suffix == ".png" and path.with_suffix(".pdf") not in paths):
                kind = "图件"
            else:
                continue
            rows.append({"面板": path.parts[-3], "组织": path.parts[-2], "内容": kind, "文件": name, "source": path})
        return rows

    def copy_panels(self) -> list[dict]:
        """把面板产出复制到交付目录 ``panels/<面板>/<组织>/``：只看交付目录就能拿到感兴趣通路的矩阵与图件。

        逐反应的档案页（图 B）与 PNG 仍只在整合结果包中。
        """
        rows = self._panel_rows()
        for row in rows:
            row["target"] = PurePosixPath(PANELS_DIR, row["面板"], row["组织"], row["文件"])
            destination = self.directory / row["target"]
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2((self.directory / row["source"]).resolve(), destination)
        return rows

    def copy_maps(self) -> list[dict]:
        """把通路总图（图 E）、分专题子通路图与概念摘要图复制到交付目录 ``pathway_map/<组织>/``。"""
        paths = {PurePosixPath(item["path"]) for item in self.artifacts
                 if item["stage"] == "integrate" and not item.get("missing")}
        rows = []
        prefix = ("integration", "reaction_dossier", "figures", "pathway_map")
        for path in sorted(paths):
            # 1) 全景总图 map_<比较> 与概念摘要 abstract_<比较>：.../pathway_map/<组织>/<文件>
            if len(path.parts) >= 6 and path.parts[-6:-2] == prefix:
                if path.name.startswith("map_"):
                    if path.suffix == ".png" and path.with_suffix(".pdf") in paths or path.suffix not in (".pdf", ".png", ".html", ".tsv"):
                        continue
                    tissue = path.parts[-2]
                    kind_label = "视图表" if path.suffix == ".tsv" else ("交互图件" if path.suffix == ".html" else "图件")
                    row = {"组织": tissue, "比较": path.stem.removeprefix("map_"),
                           "内容": kind_label, "文件": path.name,
                           "target": PurePosixPath(MAPS_DIR, tissue, path.name)}
                elif path.name.startswith("abstract_"):
                    if path.suffix == ".png" and path.with_suffix(".pdf") in paths or path.suffix not in (".pdf", ".png", ".tiff", ".tif", ".txt"):
                        continue
                    tissue = path.parts[-2]
                    kind_label = "概念摘要图注" if path.suffix == ".txt" else "概念摘要"
                    row = {"组织": tissue, "比较": path.stem.removeprefix("abstract_").removesuffix("_caption"),
                           "内容": kind_label, "文件": path.name,
                           "target": PurePosixPath(MAPS_DIR, tissue, path.name)}
                else:
                    continue
            # 2) 正文分专题子通路图：.../pathway_map/<组织>/submaps/<比较>_<子图ID>.{pdf,png,tiff} 及配套 _caption.txt
            elif len(path.parts) >= 7 and path.parts[-7:-3] == prefix and path.parts[-2] == "submaps":
                if path.suffix == ".png" and path.with_suffix(".pdf") in paths or path.suffix not in (".pdf", ".png", ".tiff", ".tif", ".txt"):
                    continue
                tissue = path.parts[-3]
                kind_label = "子通路图注" if path.suffix == ".txt" else "子通路图"
                row = {"组织": tissue, "比较": path.stem.removesuffix("_caption"),
                       "内容": kind_label, "文件": f"submaps/{path.name}",
                       "target": PurePosixPath(MAPS_DIR, tissue, "submaps", path.name)}
            else:
                continue
            destination = self.directory / row["target"]
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2((self.directory / path).resolve(), destination)
            rows.append(row)
        return rows

    def _write_map_index(self, writer: pd.ExcelWriter, rows: list[dict]) -> None:
        if not rows:
            return
        pd.DataFrame([{k: row[k] for k in ("组织", "比较", "内容", "文件")} for row in rows]).to_excel(
            writer, sheet_name="通路总图", index=False)
        sheet = writer.sheets["通路总图"]
        sheet.freeze_panes = "A2"
        sheet.auto_filter.ref = sheet.dimensions
        for column, width in zip("ABCD", (18, 24, 10, 40)):
            sheet.column_dimensions[column].width = width
        for index, row in enumerate(rows, start=2):
            cell = sheet.cell(index, 4)
            cell.hyperlink = quote(row["target"].as_posix(), safe="/")
            cell.style = "Hyperlink"

    def _write_panel_index(self, writer: pd.ExcelWriter, rows: list[dict]) -> None:
        if not rows:
            return
        pd.DataFrame([{k: row[k] for k in ("面板", "组织", "内容", "文件")} for row in rows]).to_excel(
            writer, sheet_name="面板结果", index=False)
        sheet = writer.sheets["面板结果"]
        sheet.freeze_panes = "A2"
        sheet.auto_filter.ref = sheet.dimensions
        for column, width in zip("ABCD", (32, 20, 22, 54)):
            sheet.column_dimensions[column].width = width
        for index, row in enumerate(rows, start=2):
            cell = sheet.cell(index, 4)
            cell.hyperlink = quote(row["target"].as_posix(), safe="/")
            cell.style = "Hyperlink"

    def write_candidates(self) -> pd.DataFrame | None:
        """由整合结果包的全局表写出无面板候选表（``candidates/``，见 :mod:`omicsweave.integration.candidates`），
        并把面板声明的选定方式登记到索引；返回按组织 × 比较的概览，没有整合结果包或全局表时为 None。"""
        integrate = self.bundles.get("integrate")
        if not integrate:
            return None
        from omicsweave.integration.candidates import write_candidates

        result = write_candidates((self.directory / integrate["path"]).resolve(), self.directory / CANDIDATES_DIR,
                                  self.project_root)
        if result is None:
            return None
        self.data["panels"] = result.panels.record() if result.panels else None
        return result.overview

    def _write_candidate_overview(self, writer: pd.ExcelWriter, overview: pd.DataFrame) -> None:
        target = f"{CANDIDATES_DIR}/candidates.xlsx"
        overview.assign(明细=target).to_excel(writer, sheet_name="候选概览", index=False)
        sheet = writer.sheets["候选概览"]
        sheet.freeze_panes = "C2"
        sheet.auto_filter.ref = sheet.dimensions
        for row in range(2, len(overview) + 2):
            cell = sheet.cell(row, overview.shape[1] + 1)
            cell.hyperlink = target
            cell.style = "Hyperlink"

    def write_summary(self) -> None:
        """复制面板产出（:meth:`copy_panels`）、写无面板候选表（:meth:`write_candidates`），再写 ``summary.xlsx``。"""
        panels = self.copy_panels()
        maps = self.copy_maps()
        candidates = self.write_candidates()
        self.collect_exports()
        tables, contrasts = [], []
        for item in self.artifacts:
            path = (self.directory / item["path"]).resolve()
            if item.get("missing") or path.suffix != ".parquet" or path.parent.name not in {"tables", "qc", "evaluation"}:
                continue
            tables.append({"阶段": item["stage"], "范围": path.parent.parent.name, "表": path.stem,
                           "记录数": pq.ParquetFile(path).metadata.num_rows, "文件": item["path"]})
            if path.name == "differential.parquet" and item["stage"] == "run":
                columns = ["tissue", "contrast_id", "pvalue", "selected"]
                available = pq.read_schema(path).names
                data = pd.read_parquet(path, columns=[c for c in columns if c in available])
                for (tissue, contrast), group in data.groupby(["tissue", "contrast_id"]):
                    contrasts.append({"组学": path.parent.parent.name, "组织": tissue, "比较": contrast,
                                      "结果行数": len(group), "有 p 值": group.pvalue.notna().sum(),
                                      "按本层规则入选": group.selected.sum() if "selected" in group else None})
        with pd.ExcelWriter(self.directory / "summary.xlsx", engine="openpyxl") as writer:
            pd.DataFrame([{"项目": self.data["project"], "实验": self.data["experiment"],
                           "交付 ID": self.directory.name, "阶段": name, **bundle}
                          for name, bundle in self.bundles.items()]).to_excel(writer, sheet_name="运行来源", index=False)
            pd.DataFrame(contrasts).to_excel(writer, sheet_name="差异结果概览", index=False)
            if candidates is not None:
                self._write_candidate_overview(writer, candidates)
            self._write_panel_index(writer, panels)
            self._write_map_index(writer, maps)
            pd.DataFrame(tables).to_excel(writer, sheet_name="结果表与质控", index=False)
            pd.DataFrame(self.artifacts).to_excel(writer, sheet_name="文件索引", index=False)


# ── 交付查询与 latest.json 指针 ────────────────────────────────────────────

@dataclass(frozen=True)
class DeliveryRecord:
    """一次交付目录的只读摘要，用于 latest 查询、列表展示与保留策略判定。"""
    delivery_id: str
    experiment: str
    directory: Path
    index_path: Path
    summary_path: Path | None
    status: str
    created_at: str | None
    completed_at: str | None
    tag: str | None
    pinned: bool
    bundles: dict[str, str]
    bundle_paths: dict[str, str]
    project: str | None = None
    sort_timestamp: float = 0.0
    last_activity_timestamp: float = 0.0

    def target_path(self, target: str) -> Path:
        if target == "dir":
            return self.directory
        if target == "index":
            return self.index_path
        if target == "summary":
            if self.summary_path is None or not self.summary_path.is_file():
                raise FileNotFoundError(f"交付 {self.delivery_id} 中不存在 {SUMMARY_FILE}")
            return self.summary_path
        raise ValueError(f"未知目标 {target!r}，可选：dir / summary / index")


def _parse_iso_timestamp(value: object) -> float | None:
    if not isinstance(value, str) or not value.strip():
        return None
    text = value.strip()
    if text.endswith("Z"):
        text = text[:-1] + "+00:00"
    try:
        parsed = datetime.fromisoformat(text)
    except ValueError:
        return None
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed.timestamp()


def write_latest_pointer(experiment_deliveries_dir: Path, payload: dict) -> Path:
    """原子写出 ``_deliveries/<实验>/latest.json``，仅使用相对路径以保证跨平台和可搬迁。"""
    experiment_deliveries_dir.mkdir(parents=True, exist_ok=True)
    pointer = experiment_deliveries_dir / LATEST_POINTER_FILE
    pending = experiment_deliveries_dir / f".{LATEST_POINTER_FILE}.tmp.{os.getpid()}.{secrets.token_hex(2)}"
    pending.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    pending.replace(pointer)
    return pointer


def build_latest_payload(record: DeliveryRecord) -> dict:
    """由 :class:`DeliveryRecord` 生成可写入 ``latest.json`` 的跨平台指针字典。"""
    delivery_id = record.directory.name
    summary = PurePosixPath(delivery_id, SUMMARY_FILE).as_posix() if (record.summary_path and record.summary_path.is_file()) else None
    return {
        "schema_version": LATEST_SCHEMA_VERSION,
        "project": record.project,
        "experiment": record.experiment,
        "delivery_id": delivery_id,
        "path": delivery_id,
        "summary": summary,
        "index": PurePosixPath(delivery_id, INDEX_FILE).as_posix(),
        "tag": record.tag,
        "pinned": record.pinned,
        "created_at": record.created_at,
        "completed_at": record.completed_at,
        "bundles": dict(record.bundles),
    }


def read_delivery_record(delivery_dir: Path, experiment_id: str | None = None) -> DeliveryRecord | None:
    """读取单个交付目录的 ``index.json``；非交付目录返回 None，损坏的索引标记为 ``corrupt``。"""
    index_path = delivery_dir / INDEX_FILE
    if not delivery_dir.is_dir() or not index_path.is_file():
        return None
    try:
        index_mtime = index_path.stat().st_mtime
        dir_mtime = delivery_dir.stat().st_mtime
    except OSError:
        return None
    last_activity = max(index_mtime, dir_mtime)
    pinned_file = (delivery_dir / PIN_FILE).is_file()
    summary = delivery_dir / SUMMARY_FILE
    summary_path = summary if summary.is_file() else None
    exp_name = experiment_id or delivery_dir.parent.name
    try:
        data = json.loads(index_path.read_text(encoding="utf-8"))
        if not isinstance(data, dict):
            raise ValueError("index.json 顶层不是对象")
    except (OSError, ValueError):
        return DeliveryRecord(
            delivery_id=delivery_dir.name, experiment=exp_name, directory=delivery_dir,
            index_path=index_path, summary_path=summary_path, status="corrupt",
            created_at=None, completed_at=None, tag=None, pinned=pinned_file,
            bundles={}, bundle_paths={}, project=None,
            sort_timestamp=index_mtime, last_activity_timestamp=last_activity,
        )
    raw_bundles = data.get("bundles") if isinstance(data.get("bundles"), dict) else {}
    bundles: dict[str, str] = {}
    bundle_paths: dict[str, str] = {}
    for stage, info in raw_bundles.items():
        if not isinstance(info, dict):
            continue
        if isinstance(info.get("run_id"), str) and info["run_id"]:
            bundles[stage] = info["run_id"]
        if isinstance(info.get("path"), str) and info["path"]:
            bundle_paths[stage] = info["path"]
    created_at = data.get("created_at") if isinstance(data.get("created_at"), str) and data.get("created_at") else None
    completed_at = data.get("completed_at") if isinstance(data.get("completed_at"), str) and data.get("completed_at") else None
    sort_ts = _parse_iso_timestamp(completed_at) or _parse_iso_timestamp(created_at) or index_mtime
    tag = data.get("tag") if isinstance(data.get("tag"), str) and data.get("tag") else None
    status = data.get("status") if isinstance(data.get("status"), str) and data.get("status") else "unknown"
    project = data.get("project") if isinstance(data.get("project"), str) and data.get("project") else None
    return DeliveryRecord(
        delivery_id=delivery_dir.name,
        experiment=data.get("experiment") if isinstance(data.get("experiment"), str) and data.get("experiment") else exp_name,
        directory=delivery_dir,
        index_path=index_path,
        summary_path=summary_path,
        status=status,
        created_at=created_at,
        completed_at=completed_at,
        tag=tag,
        pinned=bool(data.get("pinned")) or pinned_file,
        bundles=bundles,
        bundle_paths=bundle_paths,
        project=project,
        sort_timestamp=sort_ts,
        last_activity_timestamp=last_activity,
    )


def delivery_experiments(results_root: str | Path) -> list[str]:
    """返回 ``<results_root>/_deliveries/`` 下存在交付目录的实验 ID（按名称排序）。"""
    base = Path(results_root) / DELIVERIES_DIR
    if not base.is_dir():
        return []
    return sorted(entry.name for entry in base.iterdir() if entry.is_dir() and not entry.name.startswith("."))


def list_deliveries(results_root: str | Path, experiment_id: str) -> list[DeliveryRecord]:
    """列出某实验的全部交付记录，按时间从旧到新排序（兼容旧版 UUID 交付目录）。"""
    base = Path(results_root) / DELIVERIES_DIR / experiment_id
    if not base.is_dir():
        return []
    records: list[DeliveryRecord] = []
    for entry in base.iterdir():
        if not entry.is_dir() or entry.name.startswith("."):
            continue
        record = read_delivery_record(entry, experiment_id=experiment_id)
        if record is not None:
            records.append(record)
    records.sort(key=lambda item: (item.sort_timestamp, item.delivery_id))
    return records


def read_latest_pointer(results_root: str | Path, experiment_id: str) -> dict | None:
    """读取 ``_deliveries/<实验>/latest.json``；不存在或格式损坏时返回 None。"""
    pointer = Path(results_root) / DELIVERIES_DIR / experiment_id / LATEST_POINTER_FILE
    if not pointer.is_file():
        return None
    try:
        data = json.loads(pointer.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else None
    except (OSError, ValueError):
        return None


def resolve_latest_delivery(results_root: str | Path, experiment_id: str) -> tuple[DeliveryRecord | None, str]:
    """解析某实验最新一次成功交付，返回 ``(记录, 来源)``，来源为 ``pointer`` / ``scan`` / ``none``。

    优先读取 ``latest.json`` 指针并核对目标目录仍有效且状态为 ``completed``；若指针缺失或已失效，
    自动回退到扫描全部交付目录并取最新一次成功交付（兼容历史 UUID 交付）。
    """
    base = Path(results_root) / DELIVERIES_DIR / experiment_id
    pointer = read_latest_pointer(results_root, experiment_id)
    if pointer and isinstance(pointer.get("path"), str) and pointer["path"].strip():
        try:
            base_resolved = base.resolve()
            target = (base / pointer["path"]).resolve()
            if target.parent == base_resolved and target.is_dir():
                record = read_delivery_record(target, experiment_id=experiment_id)
                if record is not None and record.status == "completed":
                    return record, "pointer"
        except OSError:
            pass
    completed = [item for item in list_deliveries(results_root, experiment_id) if item.status == "completed"]
    if completed:
        return completed[-1], "scan"
    return None, "none"


def sync_latest_pointer(results_root: str | Path, experiment_id: str) -> Path | None:
    """根据当前磁盘上的成功交付重建或清理 ``latest.json``，确保指针始终有效。"""
    base = Path(results_root) / DELIVERIES_DIR / experiment_id
    completed = [item for item in list_deliveries(results_root, experiment_id) if item.status == "completed"]
    pointer = base / LATEST_POINTER_FILE
    if completed:
        return write_latest_pointer(base, build_latest_payload(completed[-1]))
    if pointer.is_file():
        try:
            pointer.unlink()
        except OSError:
            pass
    return None


def open_path(path: Path) -> None:
    """跨平台用系统默认程序打开文件或目录（Windows / macOS / Linux）。"""
    resolved = path.resolve()
    if os.name == "nt":
        os.startfile(str(resolved))  # type: ignore[attr-defined]
    elif sys.platform == "darwin":
        subprocess.run(["open", str(resolved)], check=True)
    else:
        subprocess.run(["xdg-open", str(resolved)], check=True)
