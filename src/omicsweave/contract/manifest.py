#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
运行清单 ``manifest.json``：每个结果包（``analysis_results/<experiment_id>/<run_id>/``）一份。

记录运行 ID、命令、配置哈希、软件版本，以及输入输出文件的校验值，用于复现与溯源。
"""
from __future__ import annotations

import hashlib
import platform
import secrets
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field

from omicsweave import __version__
from omicsweave.core.model import Omics

MANIFEST_FILE = "manifest.json"
# 交换契约 v1 在阶段 2 定稿；此前为草案版本
BUNDLE_SCHEMA_VERSION = "0.1"
_CHUNK = 1 << 20


class FileRecord(BaseModel):
    model_config = ConfigDict(extra="forbid")

    path: str          # 相对于项目目录（输入）或结果包目录（输出）
    role: Literal["input", "output"]
    sha256: str
    size: int
    kind: str | None = None


class RunManifest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    bundle_schema_version: str = BUNDLE_SCHEMA_VERSION
    run_id: str
    project_id: str
    experiment_id: str
    created_at: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))
    omicsweave_version: str = __version__
    python_version: str = Field(default_factory=platform.python_version)
    command: list[str] = Field(default_factory=list)
    config_sha256: str | None = None
    omics: list[Omics] = Field(default_factory=list)
    parameters: dict[str, Any] = Field(default_factory=dict)
    files: list[FileRecord] = Field(default_factory=list)


def new_run_id(now: datetime | None = None) -> str:
    """形如 ``run_20261029_153012_a3f9`` 的运行 ID：按时间排序，末尾随机后缀避免同秒冲突。"""
    now = now or datetime.now()
    return f"run_{now:%Y%m%d_%H%M%S}_{secrets.token_hex(2)}"


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(_CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: str | Path, base: str | Path, role: Literal["input", "output"],
                kind: str | None = None) -> FileRecord:
    path, base = Path(path), Path(base)
    try:
        shown = path.resolve().relative_to(base.resolve()).as_posix()
    except ValueError:
        shown = path.resolve().as_posix()
    return FileRecord(path=shown, role=role, sha256=sha256_file(path), size=path.stat().st_size, kind=kind)


def write_manifest(manifest: RunManifest, run_dir: str | Path) -> Path:
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    path = run_dir / MANIFEST_FILE
    path.write_text(manifest.model_dump_json(indent=2), encoding="utf-8")
    return path


def read_manifest(run_dir: str | Path) -> RunManifest:
    return RunManifest.model_validate_json((Path(run_dir) / MANIFEST_FILE).read_text(encoding="utf-8"))
