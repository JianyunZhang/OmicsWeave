"""读取 ``data/`` 下由 ``build_registry`` 生成的注册表（gzip 压缩的 TSV 与 SOURCES.json）。"""

from __future__ import annotations

import csv
import gzip
import io
import json
from importlib import resources
from typing import Any

_PACKAGE = "omicsweave.core.idmap"


def read_table(name: str) -> list[dict[str, str]]:
    """读取一张注册表，返回行字典列表；空字段为 ``""``。"""
    raw = resources.files(_PACKAGE).joinpath("data", name).read_bytes()
    with gzip.open(io.BytesIO(raw), "rt", encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def read_sources() -> dict[str, Any]:
    """注册表的来源说明（UniProt 版本、Ensembl 文件、生成日期与汇总计数）。"""
    text = resources.files(_PACKAGE).joinpath("data", "SOURCES.json").read_text(encoding="utf-8")
    return json.loads(text)


def split_ids(value: str) -> tuple[str, ...]:
    return tuple(piece for piece in (value or "").split(";") if piece)
