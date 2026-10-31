#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
项目目录定位。

omicsweave 作为独立安装的工具包，不再假设自己位于实验数据仓库内部。配置中的
``experiments/<id>/...``、``analysis_results/<id>/...`` 等相对路径，统一相对于
“项目目录”（存放实验配置、原始数据与分析结果的实验项目仓库）解析。

项目目录按以下顺序确定：

1. 显式指定：:func:`project_context` / :func:`set_project_dir`（命令行 ``omicsweave --project-dir <dir> ...``）；
2. 环境变量 ``OMICSWEAVE_PROJECT_DIR``；
3. 从当前工作目录向上查找第一个包含 ``experiments/`` 子目录的目录；
4. 当前工作目录。

显式指定的值保存在上下文变量（:mod:`contextvars`）中，而不是模块级全局变量：
每个线程、每个 asyncio 任务各自独立，互不覆盖，为以后在同一进程中并行处理多个项目
（集群作业、服务端）预留。需要临时切换项目时用 :func:`project_context`，退出后自动恢复。
子进程不继承上下文变量，需要时经环境变量传递（见 :data:`PROJECT_DIR_ENV`）。
"""
from __future__ import annotations

import os
from contextlib import contextmanager
from contextvars import ContextVar, Token
from pathlib import Path
from typing import Iterator

PROJECT_DIR_ENV = "OMICSWEAVE_PROJECT_DIR"
_MARKER_DIR = "experiments"

_explicit_project_dir: ContextVar[Path | None] = ContextVar("omicsweave_project_dir", default=None)


def _normalize(path: str | os.PathLike[str] | None) -> Path | None:
    return Path(path).expanduser().resolve() if path is not None else None


def set_project_dir(path: str | os.PathLike[str] | None) -> Token:
    """在当前上下文中显式指定项目目录；传入 ``None`` 恢复自动查找。

    返回的 token 可交给 :func:`reset_project_dir` 恢复原值。只需临时切换时，优先用 :func:`project_context`。
    """
    return _explicit_project_dir.set(_normalize(path))


def reset_project_dir(token: Token) -> None:
    """恢复 :func:`set_project_dir` 之前的值。"""
    _explicit_project_dir.reset(token)


@contextmanager
def project_context(path: str | os.PathLike[str] | None) -> Iterator[Path]:
    """在 ``with`` 块内使用指定的项目目录，退出时恢复原值；``path`` 为 None 时不改变当前设置。"""
    token = _explicit_project_dir.set(_normalize(path)) if path is not None else None
    try:
        yield project_dir()
    finally:
        if token is not None:
            _explicit_project_dir.reset(token)


def project_dir() -> Path:
    """返回当前生效的项目目录（见模块说明中的查找顺序）。"""
    explicit = _explicit_project_dir.get()
    if explicit is not None:
        return explicit
    env = os.environ.get(PROJECT_DIR_ENV, "").strip()
    if env:
        return Path(env).expanduser().resolve()
    cwd = Path.cwd().resolve()
    for candidate in (cwd, *cwd.parents):
        if (candidate / _MARKER_DIR).is_dir():
            return candidate
    return cwd
