#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""读取由 RefMet 生成的代谢物注册表（:mod:`.build_metabolite_registry`）；没有生成时为空。"""
from __future__ import annotations

from functools import lru_cache
from importlib import resources

from ._registry import read_table
from .build_metabolite_registry import TABLE_FILE


@lru_cache(maxsize=1)
def registry_by_name() -> dict[str, dict[str, str]]:
    """查询名称（小写、去首尾空白）→ 注册表行。"""
    if not resources.files("omicsweave.core.idmap").joinpath("data", TABLE_FILE).is_file():
        return {}
    return {row["query_name"].strip().lower(): row for row in read_table(TABLE_FILE)}


def lookup(name: str) -> dict[str, str] | None:
    """RefMet 匹配到数据库 ID 的名称返回 ``{hmdb_id, kegg_compound, chebi_id, refmet_name}``，否则 None。"""
    row = registry_by_name().get(str(name).strip().lower())
    if not row or row.get("status") != "matched":
        return None
    return {key: row.get(key, "") for key in ("hmdb_id", "kegg_compound", "chebi_id", "refmet_name")}
