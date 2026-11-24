#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
运行清单 ``manifest.json``：每个结果包（``results/<experiment_id>/<run_id>/``）一份。

记录运行 ID、命令、配置哈希、软件版本，以及输入输出文件的校验值，用于复现与溯源。
"""
from __future__ import annotations

import platform
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field

from omicsweave import __version__
from omicsweave.core.model import Omics
from omicsweave.core.provenance import sha256_file

MANIFEST_FILE = "manifest.json"
# 交换契约 v1 在阶段 2 定稿；此前为草案版本
BUNDLE_SCHEMA_VERSION = "0.1"
RUN_TAG_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$")


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


def validate_run_tag(tag: str) -> str:
    """校验运行/交付标签：仅允许字母、数字及 ``_`` / ``.`` / ``-``，长度 1–64，避免跨平台路径问题。"""
    cleaned = tag.strip()
    if not cleaned or not RUN_TAG_PATTERN.match(cleaned) or cleaned in {".", ".."}:
        raise ValueError(
            f"交付标签 {tag!r} 无效：请以字母或数字开头，仅使用字母、数字、_、.、-，且不超过 64 个字符"
        )
    return cleaned


def new_run_id(now: datetime | None = None, tag: str | None = None, *, seq: int = 1) -> str:
    """形如 ``run_20261029_153012``（同秒冲突时附 ``_02``..``_99``，带标签时附 ``_<tag>``）的运行 ID：按时间排序。"""
    if seq < 1:
        raise ValueError(f"同秒序号 seq 必须 >= 1，收到 {seq}")
    now = now or datetime.now()
    base = f"run_{now:%Y%m%d_%H%M%S}"
    if seq > 1:
        base = f"{base}_{seq:02d}"
    if tag is not None:
        return f"{base}_{validate_run_tag(tag)}"
    return base


def allocate_run_dir(parent: str | Path, *, now: datetime | None = None, tag: str | None = None) -> Path:
    """在 ``parent`` 下原子创建唯一的运行目录：默认 ``run_YYYYMMDD_HHMMSS``，同秒已存在时自动追加 ``_02``..``_99``。"""
    parent = Path(parent)
    ts = now or datetime.now()
    validated_tag = validate_run_tag(tag) if tag is not None else None
    for seq in range(1, 100):
        candidate = parent / new_run_id(now=ts, tag=validated_tag, seq=seq)
        try:
            candidate.mkdir(parents=True, exist_ok=False)
            return candidate
        except FileExistsError:
            continue
    raise FileExistsError(f"无法在 {parent} 创建唯一的运行目录")


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
