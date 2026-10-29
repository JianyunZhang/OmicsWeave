#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
项目目录定位。

omicsweave 作为独立安装的工具包，不再假设自己位于实验数据仓库内部。配置中的
``experiments/<id>/...``、``analysis_results/<id>/...`` 等相对路径，统一相对于
“项目目录”（存放实验配置、原始数据与分析结果的实验项目仓库）解析。

项目目录按以下顺序确定：

1. :func:`set_project_dir` 显式设置（命令行 ``omicsweave --project-dir <dir> ...``）；
2. 环境变量 ``OMICSWEAVE_PROJECT_DIR``；
3. 从当前工作目录向上查找第一个包含 ``experiments/`` 子目录的目录；
4. 当前工作目录。
"""
from __future__ import annotations

import os
from pathlib import Path

PROJECT_DIR_ENV = "OMICSWEAVE_PROJECT_DIR"
_MARKER_DIR = "experiments"

_explicit_project_dir: Path | None = None


def set_project_dir(path: str | os.PathLike[str] | None) -> None:
    """显式指定项目目录；传入 ``None`` 恢复自动查找。"""
    global _explicit_project_dir
    _explicit_project_dir = Path(path).expanduser().resolve() if path is not None else None


def project_dir() -> Path:
    """返回当前生效的项目目录（见模块说明中的查找顺序）。"""
    if _explicit_project_dir is not None:
        return _explicit_project_dir
    env = os.environ.get(PROJECT_DIR_ENV, "").strip()
    if env:
        return Path(env).expanduser().resolve()
    cwd = Path.cwd().resolve()
    for candidate in (cwd, *cwd.parents):
        if (candidate / _MARKER_DIR).is_dir():
            return candidate
    return cwd
