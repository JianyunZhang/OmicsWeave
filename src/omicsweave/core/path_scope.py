"""应用层路径约束。由可信本机宿主选择项目及可读数据根。"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path, PureWindowsPath
from stat import FILE_ATTRIBUTE_REPARSE_POINT

from .config import EXPERIMENT_FILE, PROJECT_FILE, _parse, _read_yaml, _substitute, load_experiment, load_project
from .model import ID_PATTERN, ExperimentConfig, Omics


class ScopeError(ValueError):
    """请求超出当前进程授权的项目范围。"""


def identifier(value: str) -> str:
    if not re.fullmatch(ID_PATTERN, value) or value.endswith("."):
        raise ScopeError("实验、结果包和作业标识必须是单个目录名")
    return value


def file_component(value: str) -> None:
    if value in (".", "..") or any(char in value for char in "/\\:*?<>|\0"):
        raise ScopeError(f"不能作为目录名的样本或组织标识 {value!r}")


def _indirection(path: Path) -> bool:
    st = path.lstat()
    return path.is_symlink() or bool(getattr(st, "st_file_attributes", 0) & FILE_ATTRIBUTE_REPARSE_POINT)


def _local_path(value: str | Path) -> Path:
    text = str(value)
    windows = PureWindowsPath(text)
    parts = windows.parts[1:] if windows.anchor else windows.parts
    if windows.anchor.startswith("\\\\") or "://" in text or any(":" in part for part in parts):
        raise ScopeError("MCP 使用本机文件路径，拒绝网络路径和替代数据流")
    return Path(value)


@dataclass(frozen=True)
class PathScope:
    root: Path
    data_roots: tuple[Path, ...] = ()

    @classmethod
    def create(cls, project: Path, allowed_data_roots=()) -> PathScope:
        root = _local_path(project).resolve(strict=True)
        data = tuple(_local_path(item).resolve(strict=True) for item in allowed_data_roots)
        if not root.is_dir() or any(not item.is_dir() for item in data):
            raise ScopeError("项目及外部数据根必须是本机目录")
        return cls(root, data)

    def path(self, value: str | Path, *, base: Path | None = None, external: bool = False) -> Path:
        value = _local_path(value)
        candidate = (base or self.root) / value
        resolved = _local_path(candidate.resolve())
        roots = (self.root, *self.data_roots) if external else (self.root,)
        if not any(resolved.is_relative_to(root) for root in roots):
            raise ScopeError(f"路径超出项目授权范围 {candidate}")
        return resolved

    def tree(self, path: Path, *, external: bool = False) -> Path:
        root = self.path(path, external=external)
        if root.is_dir():
            for parent, directories, files in os.walk(root, followlinks=False):
                for name in (*directories, *files):
                    child = Path(parent) / name
                    if _indirection(child):
                        raise ScopeError(f"递归读取或写入目录包含链接 {child}")
        return root

    def project(self):
        self.path(PROJECT_FILE)
        project = load_project(self.root)
        for name in ("experiments", "results", "references", "analysis"):
            self.path(project.path(name))
        for directory in project.config.experiments:
            self.path(directory)
        for reference in project.config.project.references.values():
            self.path(reference, external=True)
        return project

    def experiment(self, project, exp_id: str, *, document: dict | None = None, execution: bool = False):
        from .genesets import cache_path
        from .parameter_defaults import PROTEOMICS, TRANSCRIPTOMICS
        from .validate import _split_files

        exp_dir = self.path(project.path("experiments") / identifier(exp_id))
        path = self.path(exp_dir / EXPERIMENT_FILE)
        raw = _read_yaml(path) if document is None else document
        config = _parse(ExperimentConfig, path, _substitute(raw, project.config.project.references))
        identifier(config.experiment.id)
        self.path(config.design.subjects, base=exp_dir)
        self.path(config.design.samples, base=exp_dir)
        for layer in config.omics.enabled():
            block = config.omics.settings(layer)
            spec = block.input
            if execution and spec.level != "table":
                raise ScopeError(f"MCP 计算需要 {layer.value} 的表格层输入")
            values = [*spec.files, *(getattr(spec, "tissue_files", None) or {}).values(),
                      getattr(spec, "tx2gene", None), block.adapter.get("transcriptome_fasta")]
            for value in filter(None, values):
                self.tree(self.path(value, external=True), external=True)
            libraries = {Omics.TRANSCRIPTOMICS: TRANSCRIPTOMICS["gene_sets"],
                         Omics.PROTEOMICS: PROTEOMICS["gene_sets"]}.get(layer, [])
            for source in block.analysis.get("gene_sets", libraries) or []:
                if not isinstance(source, str):
                    raise ScopeError("基因集须为本地 GMT 路径或缓存库名")
                local = self.path(source, external=True) if source.lower().endswith(".gmt") else self.path(cache_path(project, source))
                if execution and not local.is_file():
                    raise ScopeError(f"MCP 离线计算需要可读 GMT 或已有本地缓存 {source}")

        dossier = config.integration.parameters.get("reaction_dossier", {})
        for reference in (dossier.get("panels"), (dossier.get("pathway_map") or {}).get("layout")):
            if reference:
                self.path(reference, base=exp_dir, external=True)
        if dossier.get("panels"):
            from omicsweave.integration.dossier.panels import _read

            panel_doc = _read(self.path(dossier["panels"], base=exp_dir, external=True))
            for panel_id in (panel_doc.get("panels") or {}):
                file_component(str(panel_id))
        for source in dossier.get("sources") or []:
            root = self.path(source.get("project") or ".")
            if root != self.root:
                raise ScopeError("MCP v1 的二级结果来源须位于同一项目")
            other = self.path(project.path("experiments") / identifier(source["experiment"]))
            other_cfg = self.path(other / EXPERIMENT_FILE)
            other_model = _parse(ExperimentConfig, other_cfg,
                                 _substitute(_read_yaml(other_cfg), project.config.project.references))
            self.path(other_model.design.subjects, base=other)
            self.path(other_model.design.samples, base=other)
            identifier(source.get("run") or "latest")

        experiment = load_experiment(project, exp_dir, document=raw)
        for column in ("sample_id", "tissue"):
            if column in experiment.samples:
                for value in experiment.samples[column].unique():
                    file_component(value)
        if "files" in experiment.samples:
            for value in experiment.samples["files"]:
                for name in _split_files(value):
                    self.path(name, external=True)
        return experiment
