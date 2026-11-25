#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
读取项目配置：``project.yaml``、``experiment.yaml`` 与 ``subjects.tsv`` / ``samples.tsv``
（样本表也可以是 UTF-8 的 ``.csv``，见 :func:`read_table`）。

配置中的相对路径一律相对于项目目录解析。字符串值中的 ``${project.references.<名称>}``
会替换为 ``project.yaml`` 中对应的参考文件路径。
"""
from __future__ import annotations

import codecs
import io
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

import yaml
from pydantic import BaseModel, ValidationError

from .model import ExperimentConfig, Omics, ProjectConfig
from .project import project_dir
from .tissues import resolve_tissue

if TYPE_CHECKING:
    import pandas as pd

PROJECT_FILE = "project.yaml"
EXPERIMENT_FILE = "experiment.yaml"
_VAR = re.compile(r"\$\{project\.references\.([A-Za-z0-9_]+)\}")
# 样本表的分隔符按扩展名决定：.csv 为逗号，其余为制表符
_DELIMITERS = {".csv": ","}
_DELIMITER_NAMES = {"\t": "制表符", ",": "逗号"}


class ConfigError(Exception):
    """配置无法加载。``messages`` 中每一条对应一个具体问题。"""

    def __init__(self, messages: list[str], *, code: str = "CONFIG_SCHEMA_ERROR"):
        super().__init__("\n".join(messages))
        self.messages = messages
        self.code = code


@dataclass
class Project:
    root: Path
    config: ProjectConfig

    @property
    def id(self) -> str:
        return self.config.project.id

    def resolve(self, path: str | Path) -> Path:
        """相对路径按项目目录解析；绝对路径原样返回。"""
        p = Path(path).expanduser()
        return p if p.is_absolute() else self.root / p

    def path(self, name: str) -> Path:
        """``paths`` 中声明的目录：experiments / results / references / analysis。"""
        return self.resolve(getattr(self.config.project.paths, name))

    def experiment_dirs(self) -> list[Path]:
        if self.config.experiments:
            return [self.resolve(p) for p in self.config.experiments]
        base = self.path("experiments")
        if not base.is_dir():
            return []
        return sorted(p for p in base.iterdir() if (p / EXPERIMENT_FILE).is_file())


@dataclass
class Experiment:
    project: Project
    root: Path
    config: ExperimentConfig
    subjects: pd.DataFrame
    samples: pd.DataFrame
    # 组织别名归一前的 sample_id → 原 tissue（仅当配置了 tissue_aliases 时记录）
    raw_tissues: dict[str, str] = field(default_factory=dict)

    @property
    def id(self) -> str:
        return self.config.experiment.id

    def applied_sample_aliases(self, omics: Omics | str) -> dict[str, str]:
        """某组学在 samples.tsv 中实际因 tissue_aliases 发生改名的 ``{原组织名: 规范组织名}``。"""
        if not self.raw_tissues or "sample_id" not in self.samples.columns or "tissue" not in self.samples.columns:
            return {}
        target = omics.value if isinstance(omics, Omics) else str(omics)
        rows = self.samples[self.samples["omics"] == target] if "omics" in self.samples.columns else self.samples
        applied: dict[str, str] = {}
        for sid, canon in zip(rows["sample_id"], rows["tissue"]):
            raw = self.raw_tissues.get(sid, canon)
            if raw and canon and raw != canon:
                applied.setdefault(raw, canon)
        return applied


def _format_validation_error(source: Path, error: ValidationError) -> list[str]:
    messages = []
    for item in error.errors():
        location = ".".join(str(part) for part in item["loc"]) or "(根)"
        messages.append(f"{source.name}: {location}: {item['msg']}")
    return messages


def _read_yaml(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise ConfigError([f"找不到 {path}"], code="INPUT_NOT_FOUND")
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise ConfigError([f"{path.name}: YAML 语法错误：{exc}"]) from exc
    if not isinstance(data, dict):
        raise ConfigError([f"{path.name}: 顶层必须是键值映射"])
    return data


def _parse(model: type[BaseModel], path: Path, data: dict[str, Any]) -> Any:
    try:
        return model.model_validate(data)
    except ValidationError as exc:
        raise ConfigError(_format_validation_error(path, exc)) from exc


def _substitute(value: Any, references: dict[str, str]) -> Any:
    if isinstance(value, str):
        def replace(match: re.Match) -> str:
            key = match.group(1)
            if key not in references:
                raise ConfigError([f"变量 ${{project.references.{key}}} 未在 project.yaml 的 references 中定义"])
            return references[key]
        return _VAR.sub(replace, value)
    if isinstance(value, dict):
        return {k: _substitute(v, references) for k, v in value.items()}
    if isinstance(value, list):
        return [_substitute(v, references) for v in value]
    return value


def _decode_table(path: Path) -> str:
    """UTF-8（可带 BOM）或带 BOM 的 UTF-16（Excel 的“Unicode 文本”）；其他编码报错并说明如何转存。"""
    data = path.read_bytes()
    if data.startswith((codecs.BOM_UTF16_LE, codecs.BOM_UTF16_BE)):
        return data.decode("utf-16")
    try:
        return data.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        line = data[:exc.start].count(b"\n") + 1
        raise ConfigError([
            f"样本表 {path.name} 第 {line} 行不是 UTF-8 编码（Excel 另存为“文本文件(制表符分隔)”时用的是系统编码，"
            "如中文 Windows 的 GBK）。请在 Excel 中另存为“Unicode 文本(*.txt)”后把扩展名改为 .tsv，"
            "或另存为“CSV UTF-8(逗号分隔)(*.csv)”，也可以用文本编辑器转存为 UTF-8"]) from exc


def read_table(path: Path) -> pd.DataFrame:
    """读取样本表：全部按字符串读取，去掉首尾空白，空单元格保留为空字符串。

    ``.csv`` 按逗号分隔，其余（``.tsv``、``.txt`` 等）按制表符分隔。以 ``#`` 开头的行是注释；
    单元格中的 ``#`` 是普通字符。
    """
    import pandas as pd

    if not path.is_file():
        raise ConfigError([f"找不到样本表 {path}"], code="INPUT_NOT_FOUND")
    sep = _DELIMITERS.get(path.suffix.lower(), "\t")
    lines = [line for line in _decode_table(path).splitlines() if not line.lstrip(" ").startswith("#")]
    if not any(line.strip() for line in lines):
        raise ConfigError([f"样本表 {path.name} 是空的，至少需要表头"])
    try:
        table = pd.read_csv(io.StringIO("\n".join(lines)), sep=sep, dtype=str, keep_default_na=False)
    except pd.errors.ParserError as exc:
        raise ConfigError([f"样本表 {path.name} 无法按{_DELIMITER_NAMES[sep]}分隔读取：{exc}"]) from exc
    table.columns = [c.strip() for c in table.columns]
    other = "," if sep == "\t" else "\t"
    if len(table.columns) == 1 and other in table.columns[0]:
        expected = _DELIMITER_NAMES[sep]
        suffix = ".csv" if other == "," else ".tsv"
        raise ConfigError([
            f"样本表 {path.name} 按{expected}分隔只读出一列，表头中却有{_DELIMITER_NAMES[other]}，看起来是"
            f"{_DELIMITER_NAMES[other]}分隔的文件。请改为{expected}分隔，或把扩展名改为 {suffix}"
            "（并相应修改 experiment.yaml 的 design.subjects / design.samples）"])
    return table.apply(lambda col: col.str.strip())


def load_project(root: str | Path | Project | None = None) -> Project:
    """读取项目目录下的 ``project.yaml``；省略 ``root`` 时使用 :func:`project_dir`。"""
    if isinstance(root, Project):
        return root
    root_path = Path(root).expanduser().resolve() if root is not None else project_dir()
    path = root_path / PROJECT_FILE
    return Project(root=root_path, config=_parse(ProjectConfig, path, _read_yaml(path)))


def _has_tissue_aliases(config: ExperimentConfig) -> bool:
    return bool(config.design.tissue_aliases or any(config.omics.settings(o).tissue_aliases for o in Omics))


def apply_sample_tissue_aliases(samples: pd.DataFrame,
                                config: ExperimentConfig) -> tuple[pd.DataFrame, dict[str, str]]:
    """按实验级与组学级 ``tissue_aliases`` 归一 ``samples.tsv`` 的 ``tissue`` 列，并保留原名索引。"""
    if "tissue" not in samples.columns or not _has_tissue_aliases(config):
        return samples, {}
    raw_tissues = dict(zip(samples["sample_id"], samples["tissue"])) if "sample_id" in samples.columns else {}
    aliases_by_omics = {o.value: config.tissue_aliases_for(o) for o in Omics}
    default_aliases = config.tissue_aliases_for(None)
    omics_col = samples["omics"].tolist() if "omics" in samples.columns else [None] * len(samples)
    resolved = [
        resolve_tissue(tissue, aliases_by_omics.get(omics, default_aliases)) if tissue else ""
        for tissue, omics in zip(samples["tissue"].tolist(), omics_col)
    ]
    updated = samples.copy()
    updated["tissue"] = resolved
    return updated, raw_tissues


def load_experiment(project: Project, experiment: str | Path, *, document: dict | None = None) -> Experiment:
    """按实验 ID（``paths.experiments`` 下的目录名）或实验目录路径读取实验。"""
    candidate = Path(experiment)
    # 只有像路径的参数才当作目录，避免实验 ID 恰好与当前目录下的文件夹同名
    is_path = candidate.is_absolute() or len(candidate.parts) > 1
    exp_dir = candidate if is_path else project.path("experiments") / str(experiment)
    path = exp_dir / EXPERIMENT_FILE
    data = _substitute(_read_yaml(path) if document is None else document, project.config.project.references)
    config: ExperimentConfig = _parse(ExperimentConfig, path, data)
    subjects = read_table(exp_dir / config.design.subjects)
    samples, raw_tissues = apply_sample_tissue_aliases(read_table(exp_dir / config.design.samples), config)
    return Experiment(project=project, root=exp_dir, config=config, subjects=subjects, samples=samples,
                      raw_tissues=raw_tissues)
