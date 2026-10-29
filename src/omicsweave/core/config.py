#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
读取项目配置：``project.yaml``、``experiment.yaml`` 与 ``subjects.tsv`` / ``samples.tsv``。

配置中的相对路径一律相对于项目目录解析。字符串值中的 ``${project.references.<名称>}``
会替换为 ``project.yaml`` 中对应的参考文件路径。
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
import yaml
from pydantic import BaseModel, ValidationError

from .model import ExperimentConfig, ProjectConfig
from .project import project_dir

PROJECT_FILE = "project.yaml"
EXPERIMENT_FILE = "experiment.yaml"
_VAR = re.compile(r"\$\{project\.references\.([A-Za-z0-9_]+)\}")


class ConfigError(Exception):
    """配置无法加载。``messages`` 中每一条对应一个具体问题。"""

    def __init__(self, messages: list[str]):
        super().__init__("\n".join(messages))
        self.messages = messages


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

    @property
    def id(self) -> str:
        return self.config.experiment.id


def _format_validation_error(source: Path, error: ValidationError) -> list[str]:
    messages = []
    for item in error.errors():
        location = ".".join(str(part) for part in item["loc"]) or "(根)"
        messages.append(f"{source.name}: {location}: {item['msg']}")
    return messages


def _read_yaml(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise ConfigError([f"找不到 {path}"])
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


def read_table(path: Path) -> pd.DataFrame:
    """读取 TSV 样本表：全部按字符串读取，去掉首尾空白，空单元格保留为空字符串。"""
    if not path.is_file():
        raise ConfigError([f"找不到样本表 {path}"])
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, comment="#")
    table.columns = [c.strip() for c in table.columns]
    return table.apply(lambda col: col.str.strip())


def load_project(root: str | Path | None = None) -> Project:
    """读取项目目录下的 ``project.yaml``；省略 ``root`` 时使用 :func:`project_dir`。"""
    root_path = Path(root).expanduser().resolve() if root is not None else project_dir()
    path = root_path / PROJECT_FILE
    return Project(root=root_path, config=_parse(ProjectConfig, path, _read_yaml(path)))


def load_experiment(project: Project, experiment: str | Path) -> Experiment:
    """按实验 ID（``paths.experiments`` 下的目录名）或实验目录路径读取实验。"""
    candidate = Path(experiment)
    # 只有像路径的参数才当作目录，避免实验 ID 恰好与当前目录下的文件夹同名
    is_path = candidate.is_absolute() or len(candidate.parts) > 1
    exp_dir = candidate if is_path else project.path("experiments") / str(experiment)
    path = exp_dir / EXPERIMENT_FILE
    data = _substitute(_read_yaml(path), project.config.project.references)
    config: ExperimentConfig = _parse(ExperimentConfig, path, data)
    subjects = read_table(exp_dir / config.design.subjects)
    samples = read_table(exp_dir / config.design.samples)
    return Experiment(project=project, root=exp_dir, config=config, subjects=subjects, samples=samples)
