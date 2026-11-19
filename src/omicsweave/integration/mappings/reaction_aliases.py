#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
反应 ID 别名（G18 M0a）：旧 MFA 模型命名的映射表条目与网络反应的对应，见 ``reaction_aliases.yaml`` 的说明。

- :func:`reaction_aliases`：``{别名: (对应的反应, ...)}``；
- :func:`assessment_universe`：反应档案的评估集——映射表条目与 MFA 结果中反应 ID 的并集，去掉别名（MFA 结果用到的除外）。
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Iterable

import yaml

ALIAS_FILE = Path(__file__).with_name("reaction_aliases.yaml")


class AliasError(ValueError):
    """别名表格式错误，或别名与网络反应冲突。"""


@lru_cache(maxsize=1)
def reaction_aliases() -> dict[str, tuple[str, ...]]:
    doc = yaml.safe_load(ALIAS_FILE.read_text(encoding="utf-8")) or {}
    out: dict[str, tuple[str, ...]] = {}
    for alias, spec in (doc.get("aliases") or {}).items():
        targets = spec.get("of") if isinstance(spec, dict) else None
        if not targets or not isinstance(targets, list):
            raise AliasError(f"{ALIAS_FILE.name}：别名 {alias} 缺少 of 列表")
        out[str(alias)] = tuple(str(t) for t in targets)
    return out


def validate_aliases(network_ids: Iterable[str], mapping_ids: Iterable[str]) -> list[str]:
    """别名不能是网络反应；对应的反应必须是网络反应或仍保留的（不是别名的）映射表条目。返回问题列表。"""
    network, mapping = set(network_ids), set(mapping_ids)
    aliases = reaction_aliases()
    problems = []
    for alias, targets in aliases.items():
        if alias in network:
            problems.append(f"{alias} 是网络反应，不能登记为别名")
        for target in targets:
            if target in aliases:
                problems.append(f"{alias} 对应的 {target} 本身是别名")
            elif target not in network and target not in mapping:
                problems.append(f"{alias} 对应的 {target} 既不在网络表中，也不在映射表中")
    return problems


def assessment_universe(mapping_ids: Iterable[str], mfa_ids: Iterable[str]) -> tuple[list[str], list[str]]:
    """（评估集, 被排除的别名）。MFA 结果中出现的 ID 总是保留，旧模型命名的项目因此仍可评估。"""
    aliases = reaction_aliases()
    measured = {str(i) for i in mfa_ids}
    candidates = {str(i) for i in mapping_ids} | measured
    excluded = sorted(i for i in candidates if i in aliases and i not in measured)
    return sorted(candidates - set(excluded)), excluded
