"""一次 pipeline 调用的交付清单、latest.json 指针与 Excel 摘要。结果包通过相对路径引用。"""
from __future__ import annotations

import json
import os
import re
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

from omicsweave.contract.manifest import (
    MANIFEST_FILE,
    allocate_run_dir,
    new_run_id,
    read_manifest,
    validate_run_tag,
)
from .config import EXPERIMENT_FILE, PROJECT_FILE, Experiment
from .file_state import file_lock, write_json
from .provenance import sha256_file

INDEX_FILE = "index.json"
SUMMARY_FILE = "summary.xlsx"
LATEST_POINTER_FILE = "latest.json"
LATEST_SCHEMA_VERSION = 2
PIN_FILE = ".keep"
validate_delivery_tag = validate_run_tag

# 交付目录中反应档案面板产出的子目录，与无面板候选表的子目录
PANELS_DIR = "panels"
CANDIDATES_DIR = "candidates"
# 分组织代谢通路总图（图 E，G18）的副本
MAPS_DIR = "pathway_map"


def new_delivery_id(now: datetime | None = None, tag: str | None = None, *, seq: int = 1) -> str:
    """形如 ``run_20261119_233012``（同秒冲突时附 ``_02``..``_99``，带标签时附 ``_<tag>``）的运行与交付 ID。"""
    return new_run_id(now=now, tag=tag, seq=seq)


def _link_or_copy(src: Path, dst: Path) -> None:
    """优先建立硬链接以零拷贝共享文件，跨文件系统或不支持时回退复制。"""
    src_resolved = src.resolve()
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        try:
            if os.path.samefile(src_resolved, dst):
                return
        except OSError:
            pass
        dst.unlink()
    try:
        os.link(src_resolved, dst)
    except OSError:
        shutil.copy2(src_resolved, dst)


class Delivery:
    def __init__(self, experiment: Experiment, scripts: list[str], command: list[str], *,
                 tag: str | None = None, pinned: bool = False, now: datetime | None = None,
                 job_id: str | None = None, run_tag: str | None = None):
        project = experiment.project
        self.project_root = project.root
        self.results_root = project.path("results")
        self._default_experiment = (
            project.config.pipeline.experiments[0]
            if project.config.pipeline.experiments
            else experiment.id
        )
        validated_tag = validate_delivery_tag(tag) if tag is not None else None
        parent = self.results_root / experiment.id
        self.directory = allocate_run_dir(parent, now=now, tag=run_tag or validated_tag)
        if pinned:
            (self.directory / PIN_FILE).write_text("pinned\n", encoding="utf-8")
        self.stages = {"validate": {"status": "completed"}}
        self.stages.update({key: {"status": "pending"} for key in ["run", "integrate", *scripts, "summary"]})
        self.bundles: dict[str, dict] = {}
        self.artifacts: list[dict] = []
        self.data = {
            "run_id": self.directory.name,
            "delivery_id": self.directory.name, "project": project.id, "experiment": experiment.id,
            "created_at": datetime.now(timezone.utc).isoformat(), "tag": validated_tag, "pinned": bool(pinned),
            "command": command, "status": "running",
            "stages": self.stages, "bundles": self.bundles,
            "configuration": {str(path): sha256_file(path) for path in
                              (project.root / PROJECT_FILE, experiment.root / EXPERIMENT_FILE)},
            "artifacts": self.artifacts,
        }
        if job_id is not None:
            self.data["job_id"] = job_id
        self.save()

    def save(self) -> None:
        index = self.directory / INDEX_FILE
        pending = self.directory / ".index.tmp"
        pending.write_text(json.dumps(self.data, ensure_ascii=False, indent=2), encoding="utf-8")
        pending.replace(index)

    def complete(self) -> Path:
        """标记交付完成并原子更新结果根目录下的 ``latest.json`` 指针。"""
        self.data["status"] = "completed"
        self.data["completed_at"] = datetime.now(timezone.utc).isoformat()
        self.save()
        return self.update_latest_pointer()

    def update_latest_pointer(self) -> Path:
        """写出指向本次交付的跨平台相对路径 JSON 指针（``<结果根>/latest.json``）。"""
        exp_id = self.data["experiment"]
        run_id = self.directory.name
        rel_dir = PurePosixPath(exp_id, run_id)
        summary = (rel_dir / SUMMARY_FILE).as_posix() if (self.directory / SUMMARY_FILE).is_file() else None
        entry = {
            "experiment": exp_id,
            "run_id": run_id,
            "delivery_id": run_id,
            "path": rel_dir.as_posix(),
            "summary": summary,
            "index": (rel_dir / INDEX_FILE).as_posix(),
            "tag": self.data.get("tag"),
            "pinned": bool(self.data.get("pinned")) or (self.directory / PIN_FILE).is_file(),
            "created_at": self.data.get("created_at"),
            "completed_at": self.data.get("completed_at"),
            "bundles": {stage: info["run_id"] for stage, info in self.bundles.items() if isinstance(info, dict) and info.get("run_id")},
        }
        return write_latest_pointer(
            self.results_root,
            entry,
            project_id=self.data.get("project"),
            default_experiment=self._default_experiment,
        )

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
        same_dir = directory.resolve() == self.directory.resolve()
        base = PurePosixPath(".") if same_dir else PurePosixPath(Path(os.path.relpath(directory, self.directory)).as_posix())
        manifest_sha = sha256_file(directory / MANIFEST_FILE)
        self.bundles[stage] = {"run_id": manifest.run_id, "path": base.as_posix(),
                               "config_sha256": manifest.config_sha256, "created_at": manifest.created_at.isoformat(),
                               "manifest_sha256": manifest_sha}
        run_info = self.bundles.get("run")
        same_as_run = isinstance(run_info, dict) and (
            same_dir if run_info.get("path") == "." else
            (self.directory / str(run_info.get("path", ""))).resolve() == directory.resolve()
        )
        if stage == "integrate" and same_as_run and isinstance(run_info, dict):
            run_info["manifest_sha256"] = manifest_sha
        self.artifacts[:] = [item for item in self.artifacts if item["stage"] != stage]
        for item in manifest.files:
            if item.role != "output":
                continue
            top = PurePosixPath(item.path).parts[0] if item.path else ""
            if stage == "run" and top == "integration":
                continue
            if stage == "integrate" and same_as_run and top != "integration":
                continue
            record = {"stage": stage, "path": (base / item.path).as_posix(),
                      "size": item.size, "sha256": item.sha256}
            if not (directory / item.path).is_file():
                record["missing"] = True
            self.artifacts.append(record)

    def collect_exports(self) -> None:
        self.artifacts[:] = [item for item in self.artifacts if item["stage"] != "delivery"]
        ignored_root = {
            self.directory / INDEX_FILE,
            self.directory / ".index.tmp",
            self.directory / PIN_FILE,
            self.directory / MANIFEST_FILE,
        }
        ignored_dirs = {"omics", "integration", "metadata", "qc", "exports"}
        for path in sorted(self.directory.rglob("*")):
            if not path.is_file() or path in ignored_root:
                continue
            rel = path.relative_to(self.directory)
            if rel.parts[0] in ignored_dirs:
                continue
            self.artifacts.append({"stage": "delivery", "path": rel.as_posix(),
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
        """把面板产出链接（或复制）到交付目录 ``panels/<面板>/<组织>/``：只看交付目录就能拿到感兴趣通路的矩阵与图件。

        逐反应的档案页（图 B）与 PNG 仍只在整合结果包中。
        """
        rows = self._panel_rows()
        for row in rows:
            row["target"] = PurePosixPath(PANELS_DIR, row["面板"], row["组织"], row["文件"])
            destination = self.directory / row["target"]
            _link_or_copy(self.directory / row["source"], destination)
        return rows

    def copy_maps(self) -> list[dict]:
        """把通路总图（图 E）、分专题子通路图与概念摘要图链接（或复制）到交付目录 ``pathway_map/<组织>/``。"""
        paths = {PurePosixPath(item["path"]) for item in self.artifacts
                 if item["stage"] == "integrate" and not item.get("missing")}
        rows = []
        prefix = ("integration", "reaction_dossier", "figures", "pathway_map")
        for path in sorted(paths):
            overview = path.parts[-6:-2] == prefix
            by_tracer = path.parts[-8:-4] == prefix and path.parts[-3] == "by_tracer"
            if path.suffix == ".png" and path.with_suffix(".pdf") in paths:
                continue
            if (overview or by_tracer) and path.name.startswith("map_"):
                tissue = path.parts[-4] if by_tracer else path.parts[-2]
                relative = PurePosixPath(*path.parts[-3:]) if by_tracer else PurePosixPath(path.name)
                kind = {".tsv": "视图表", ".html": "交互图件", ".txt": "图注"}.get(path.suffix, "图件")
                contrast = path.stem.removeprefix("map_")
                row = {"组织": tissue, "比较": contrast.removesuffix("_caption") if path.suffix == ".txt" else contrast,
                       "内容": f"{path.parts[-2]} · {kind}" if by_tracer else kind, "文件": relative.as_posix(),
                       "target": PurePosixPath(MAPS_DIR, tissue, relative)}
            elif overview and path.name.startswith("abstract_"):
                tissue = path.parts[-2]
                kind_label = "概念摘要图注" if path.suffix == ".txt" else "概念摘要"
                row = {"组织": tissue, "比较": path.stem.removeprefix("abstract_").removesuffix("_caption"),
                       "内容": kind_label, "文件": path.name, "target": PurePosixPath(MAPS_DIR, tissue, path.name)}
            # 2) 正文分专题子通路图：.../pathway_map/<组织>/submaps/<比较>_<子图ID>.{pdf,png,tiff} 及配套 _caption.txt
            elif len(path.parts) >= 7 and path.parts[-7:-3] == prefix and path.parts[-2] == "submaps":
                tissue = path.parts[-3]
                kind_label = "子通路图注" if path.suffix == ".txt" else "子通路图"
                row = {"组织": tissue, "比较": path.stem.removesuffix("_caption"),
                       "内容": kind_label, "文件": f"submaps/{path.name}",
                       "target": PurePosixPath(MAPS_DIR, tissue, "submaps", path.name)}
            else:
                continue
            destination = self.directory / row["target"]
            _link_or_copy(self.directory / path, destination)
            rows.append(row)
        return rows

    def _write_mfa_index(self, writer: pd.ExcelWriter) -> None:
        """链接到本次运行包中的逐示踪剂和 pooled 统计表及火山图。"""
        run = self.bundles.get("run")
        if not run:
            return
        mfa = PurePosixPath(run["path"], "omics", "mfa")
        paths = {PurePosixPath(item["path"]) for item in self.artifacts
                 if item["stage"] == "run" and not item.get("missing")}
        if mfa / "tables" / "differential_by_tracer.parquet" not in paths:
            return
        rows = []
        for path in sorted(p for p in paths if p.is_relative_to(mfa)):
            parts = path.relative_to(mfa).parts
            if parts == ("tables", "differential_pooled.tsv"):
                scope, kind = "pooled", "差异统计"
            elif len(parts) == 4 and parts[:2] == ("tables", "by_tracer") and parts[-1] == "differential.tsv":
                scope, kind = parts[2], "差异统计"
            elif parts[0] == "figures" and path.name.startswith("volcano_") and path.suffix in (".pdf", ".png"):
                if path.suffix == ".png" and path.with_suffix(".pdf") in paths:
                    continue
                scope = parts[2] if len(parts) == 4 and parts[1] == "by_tracer" else "pooled"
                kind = "火山图"
            else:
                continue
            rows.append({"示踪剂或合并": scope, "内容": kind, "文件": path.name, "target": path})
        if rows:
            self._write_link_index(writer, rows, "MFA统计与图", ["示踪剂或合并", "内容", "文件"], [28, 16, 52])

    def _write_map_index(self, writer: pd.ExcelWriter, rows: list[dict]) -> None:
        if not rows:
            return
        self._write_link_index(writer, rows, "通路总图", ["组织", "比较", "内容", "文件"], [18, 24, 28, 64])

    @staticmethod
    def _write_link_index(writer: pd.ExcelWriter, rows: list[dict], title: str,
                          columns: list[str], widths: list[int]) -> None:
        pd.DataFrame(rows, columns=columns).to_excel(writer, sheet_name=title, index=False)
        sheet = writer.sheets[title]
        sheet.freeze_panes = "A2"
        sheet.auto_filter.ref = sheet.dimensions
        for column, width in zip("ABCD", widths):
            sheet.column_dimensions[column].width = width
        for index, row in enumerate(rows, start=2):
            cell = sheet.cell(index, len(columns))
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
            self._write_mfa_index(writer)
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


def latest_pointer_path(results_root: str | Path) -> Path:
    return Path(results_root) / LATEST_POINTER_FILE


def write_latest_pointer(
    results_root: Path,
    entry: dict,
    *,
    project_id: str | None = None,
    default_experiment: str | None = None,
) -> Path:
    """原子更新 ``<results_root>/latest.json`` 中的实验指针，仅使用相对路径以保证跨平台和可搬迁。"""
    results_root = Path(results_root)
    results_root.mkdir(parents=True, exist_ok=True)
    with file_lock(results_root / ".omicsweave-latest.lock"):
        return _write_latest_pointer(results_root, entry, project_id, default_experiment)


def _write_latest_pointer(results_root: Path, entry: dict, project_id: str | None,
                          default_experiment: str | None) -> Path:
    pointer = latest_pointer_path(results_root)
    existing = read_latest_pointer(results_root) or {}
    experiments = dict(existing.get("experiments")) if isinstance(existing.get("experiments"), dict) else {}
    exp_id = entry["experiment"]
    current = experiments.get(exp_id)
    if isinstance(current, dict):
        current_time = _parse_iso_timestamp(current.get("completed_at"))
        entry_time = _parse_iso_timestamp(entry.get("completed_at"))
        current_path = current.get("path")
        if (current_time is not None and entry_time is not None and current_time > entry_time
                and isinstance(current_path, str) and (results_root / current_path / INDEX_FILE).is_file()):
            return pointer
    experiments[exp_id] = entry
    payload = {
        "schema_version": LATEST_SCHEMA_VERSION,
        "project": project_id or existing.get("project"),
        "default_experiment": default_experiment or existing.get("default_experiment") or exp_id,
        "updated_at": datetime.now(timezone.utc).isoformat(),
        "experiments": dict(sorted(experiments.items())),
    }
    write_json(pointer, payload)
    return pointer


def build_latest_payload(record: DeliveryRecord) -> dict:
    """由 :class:`DeliveryRecord` 生成可写入 ``latest.json`` 的跨平台指针条目字典。"""
    run_id = record.directory.name
    rel_dir = PurePosixPath(record.experiment, run_id)
    summary = (rel_dir / SUMMARY_FILE).as_posix() if (record.summary_path and record.summary_path.is_file()) else None
    return {
        "experiment": record.experiment,
        "run_id": run_id,
        "delivery_id": run_id,
        "path": rel_dir.as_posix(),
        "summary": summary,
        "index": (rel_dir / INDEX_FILE).as_posix(),
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
    """返回 ``<results_root>/`` 下存在交付目录（含 ``index.json``）的实验 ID（按名称排序）。"""
    root = Path(results_root)
    if not root.is_dir():
        return []
    found = []
    for entry in root.iterdir():
        if not entry.is_dir() or entry.name.startswith((".", "_")):
            continue
        if any(
            child.is_dir() and not child.name.startswith((".", "_")) and (child / INDEX_FILE).is_file()
            for child in entry.iterdir()
        ):
            found.append(entry.name)
    return sorted(found)


def list_deliveries(results_root: str | Path, experiment_id: str) -> list[DeliveryRecord]:
    """列出某实验的全部交付记录，按时间从旧到新排序。"""
    base = Path(results_root) / experiment_id
    if not base.is_dir():
        return []
    records: list[DeliveryRecord] = []
    for entry in base.iterdir():
        if not entry.is_dir() or entry.name.startswith((".", "_")):
            continue
        record = read_delivery_record(entry, experiment_id=experiment_id)
        if record is not None:
            records.append(record)
    records.sort(key=lambda item: (item.sort_timestamp, item.delivery_id))
    return records


def read_latest_pointer(results_root: str | Path, experiment_id: str | None = None) -> dict | None:
    """读取 ``<results_root>/latest.json``（或指定实验的条目）；不存在或格式损坏时返回 None。"""
    pointer = latest_pointer_path(results_root)
    if not pointer.is_file():
        return None
    try:
        data = json.loads(pointer.read_text(encoding="utf-8"))
        if not isinstance(data, dict):
            return None
    except (OSError, ValueError):
        return None
    if experiment_id is None:
        return data
    experiments = data.get("experiments")
    if not isinstance(experiments, dict):
        return None
    entry = experiments.get(experiment_id)
    return entry if isinstance(entry, dict) else None


def resolve_latest_delivery(results_root: str | Path, experiment_id: str) -> tuple[DeliveryRecord | None, str]:
    """解析某实验最新一次成功交付，返回 ``(记录, 来源)``，来源为 ``pointer`` / ``scan`` / ``none``。

    优先读取 ``<results_root>/latest.json`` 指针并核对目标目录仍有效且状态为 ``completed``；
    若指针缺失或已失效，自动回退到扫描该实验的全部交付目录并取最新一次成功交付。
    """
    root = Path(results_root)
    exp_dir = root / experiment_id
    entry = read_latest_pointer(root, experiment_id)
    if entry and isinstance(entry.get("path"), str) and entry["path"].strip():
        try:
            exp_resolved = exp_dir.resolve()
            target = (root / entry["path"]).resolve()
            if target.parent == exp_resolved and target.is_dir():
                record = read_delivery_record(target, experiment_id=experiment_id)
                if record is not None and record.status == "completed":
                    return record, "pointer"
        except OSError:
            pass
    completed = [item for item in list_deliveries(root, experiment_id) if item.status == "completed"]
    if completed:
        return completed[-1], "scan"
    return None, "none"


def sync_latest_pointer(results_root: str | Path, experiment_id: str) -> Path | None:
    """根据当前磁盘上的成功交付重建或清理 ``<results_root>/latest.json`` 中的实验条目。"""
    root = Path(results_root)
    if not root.is_dir():
        return None
    with file_lock(root / ".omicsweave-latest.lock"):
        return _sync_latest_pointer(root, experiment_id)


def _sync_latest_pointer(root: Path, experiment_id: str) -> Path | None:
    completed = [item for item in list_deliveries(root, experiment_id) if item.status == "completed"]
    if completed:
        latest_rec = completed[-1]
        return _write_latest_pointer(root, build_latest_payload(latest_rec), latest_rec.project, None)
    pointer = latest_pointer_path(root)
    existing = read_latest_pointer(root)
    if not existing or not isinstance(existing.get("experiments"), dict):
        if pointer.is_file():
            try:
                pointer.unlink()
            except OSError:
                pass
        return None
    experiments = dict(existing["experiments"])
    if experiment_id in experiments:
        del experiments[experiment_id]
        if not experiments:
            try:
                pointer.unlink()
            except OSError:
                pass
            return None
        default_exp = existing.get("default_experiment")
        if default_exp == experiment_id:
            default_exp = next(iter(sorted(experiments)))
        payload = {
            "schema_version": LATEST_SCHEMA_VERSION,
            "project": existing.get("project"),
            "default_experiment": default_exp,
            "updated_at": datetime.now(timezone.utc).isoformat(),
            "experiments": dict(sorted(experiments.items())),
        }
        write_json(pointer, payload)
        return pointer
    return pointer if pointer.is_file() else None


def open_path(path: Path) -> None:
    """跨平台用系统默认程序打开文件或目录（Windows / macOS / Linux）。"""
    resolved = path.resolve()
    if os.name == "nt":
        os.startfile(str(resolved))  # type: ignore[attr-defined]
    elif sys.platform == "darwin":
        subprocess.run(["open", str(resolved)], check=True)
    else:
        subprocess.run(["xdg-open", str(resolved)], check=True)
