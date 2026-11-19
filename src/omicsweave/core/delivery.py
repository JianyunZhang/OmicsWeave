"""一次 pipeline 调用的交付清单与 Excel 摘要。结果包通过相对路径引用。"""
from __future__ import annotations

import json
import os
import shutil
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Iterator
from urllib.parse import quote
from uuid import uuid4

import pandas as pd
import pyarrow.parquet as pq

from omicsweave.contract.manifest import MANIFEST_FILE, read_manifest
from .config import EXPERIMENT_FILE, PROJECT_FILE, Experiment
from .provenance import sha256_file

# 交付目录中反应档案面板产出的子目录，与无面板候选表的子目录
PANELS_DIR = "panels"
CANDIDATES_DIR = "candidates"
# 分组织代谢通路总图（图 E，G18）的副本
MAPS_DIR = "pathway_map"


class Delivery:
    def __init__(self, experiment: Experiment, scripts: list[str], command: list[str]):
        project = experiment.project
        self.project_root = project.root
        self.directory = project.path("results") / "_deliveries" / experiment.id / uuid4().hex
        self.directory.mkdir(parents=True, exist_ok=False)
        self.stages = {"validate": {"status": "completed"}}
        self.stages.update({key: {"status": "pending"} for key in ["run", "integrate", *scripts, "summary"]})
        self.bundles: dict[str, dict] = {}
        self.artifacts: list[dict] = []
        self.data = {
            "delivery_id": self.directory.name, "project": project.id, "experiment": experiment.id,
            "created_at": datetime.now(timezone.utc).isoformat(), "command": command, "status": "running",
            "stages": self.stages, "bundles": self.bundles,
            "configuration": {str(path): sha256_file(path) for path in
                              (project.root / PROJECT_FILE, experiment.root / EXPERIMENT_FILE)},
            "artifacts": self.artifacts,
        }
        self.save()

    def save(self) -> None:
        index = self.directory / "index.json"
        pending = self.directory / ".index.tmp"
        pending.write_text(json.dumps(self.data, ensure_ascii=False, indent=2), encoding="utf-8")
        pending.replace(index)

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
        for path in sorted(self.directory.rglob("*")):
            if path.is_file() and path not in (self.directory / "index.json", self.directory / ".index.tmp"):
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
            shutil.copy2(self.directory / row["source"], destination)
        return rows

    def copy_maps(self) -> list[dict]:
        """把通路总图（图 E）复制到交付目录 ``pathway_map/<组织>/``：每页的 PDF（没有时取 PNG）与同名视图表。"""
        paths = {PurePosixPath(item["path"]) for item in self.artifacts
                 if item["stage"] == "integrate" and not item.get("missing")}
        rows = []
        for path in sorted(paths):
            if path.parts[-6:-2] != ("integration", "reaction_dossier", "figures", "pathway_map") \
                    or not path.name.startswith("map_"):
                continue
            if path.suffix == ".png" and path.with_suffix(".pdf") in paths or path.suffix not in (".pdf", ".png", ".tsv"):
                continue
            row = {"组织": path.parts[-2], "比较": path.stem.removeprefix("map_"),
                   "内容": "视图表" if path.suffix == ".tsv" else "图件", "文件": path.name,
                   "target": PurePosixPath(MAPS_DIR, path.parts[-2], path.name)}
            destination = self.directory / row["target"]
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(self.directory / path, destination)
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

        result = write_candidates(self.directory / integrate["path"], self.directory / CANDIDATES_DIR,
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
            path = self.directory / item["path"]
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
