#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
实验级与组学级组织名称归一（``design.tissue_aliases`` / ``omics.<组学>.tissue_aliases``）。

别名表均为 ``{别名: 规范组织名}``，例如::

    design:
      tissue_aliases:
        white_adipose: sq_adipose
        cerebral_cortex: brain
        renal_cortex: kidney

组学级 ``omics.<组学>.tissue_aliases`` 合并并覆写实验级映射；显式自映射（如 ``white_adipose: white_adipose``）
可在某个组学中取消实验级重命名。查找时优先精确匹配，再按大小写不敏感匹配。
"""
from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import pandas as pd

_BAD_PATH_CHARS = ("/", "\\")


def _clean_name(value: Any, kind: str, label: str) -> str:
    text = str(value).strip() if value is not None else ""
    if not text:
        raise ValueError(f"{label} 中的{kind}不能为空")
    return text


def _validate_target(target: str, label: str) -> None:
    if target in (".", "..") or any(ch in target for ch in _BAD_PATH_CHARS):
        raise ValueError(f"{label} 中的规范组织名 {target!r} 不能包含路径分隔符或为 '.' / '..'")


def _resolve_chains(entries: dict[str, str], label: str) -> dict[str, str]:
    """展开链式映射（a → b → c 记为 a → c）并拒绝环路；允许自映射与仅改大小写的映射。"""
    by_fold = {k.casefold(): (k, v) for k, v in entries.items()}
    resolved: dict[str, str] = {}
    for key, target in entries.items():
        chain = [key]
        seen = {key.casefold()}
        current = target
        while True:
            cur_fold = current.casefold()
            # 自映射或仅修改大小写（如 White_Adipose → white_adipose）视为终点
            if cur_fold == chain[-1].casefold():
                break
            nxt = by_fold.get(cur_fold)
            if nxt is None:
                break
            _, next_target = nxt
            if next_target.casefold() == cur_fold:
                current = next_target
                break
            if cur_fold in seen:
                cycle = " → ".join([*chain, current])
                raise ValueError(f"{label} 中存在循环映射：{cycle}")
            chain.append(current)
            seen.add(cur_fold)
            current = next_target
        resolved[key] = current
    return resolved


def validate_tissue_aliases(aliases: Mapping[Any, Any] | None, label: str = "tissue_aliases") -> dict[str, str]:
    """清洗并校验一张组织别名表，返回展开链式映射后的 ``{别名: 规范组织名}``。"""
    if not aliases:
        return {}
    cleaned: dict[str, str] = {}
    seen_fold: dict[str, tuple[str, str]] = {}
    for raw_key, raw_val in aliases.items():
        key = _clean_name(raw_key, "别名键", label)
        target = _clean_name(raw_val, "规范组织名", label)
        _validate_target(target, label)
        fold = key.casefold()
        if fold in seen_fold:
            prev_key, prev_target = seen_fold[fold]
            if prev_target != target:
                raise ValueError(
                    f"{label} 中存在仅大小写不同的冲突别名：{prev_key!r} → {prev_target!r} 与 {key!r} → {target!r}"
                )
            continue
        seen_fold[fold] = (key, target)
        cleaned[key] = target
    return _resolve_chains(cleaned, label)


def merge_tissue_aliases(*maps: Mapping[Any, Any] | None, label: str = "tissue_aliases") -> dict[str, str]:
    """按顺序合并多张组织别名表（后者按大小写不敏感覆写前者同名键），再做环路校验与链式展开。"""
    merged: dict[str, str] = {}
    fold_to_key: dict[str, str] = {}
    for mapping in maps:
        if not mapping:
            continue
        cleaned = validate_tissue_aliases(mapping, label)
        for key, target in cleaned.items():
            fold = key.casefold()
            if fold in fold_to_key and fold_to_key[fold] != key:
                merged.pop(fold_to_key[fold], None)
            merged[key] = target
            fold_to_key[fold] = key
    return _resolve_chains(merged, label) if merged else {}


def _lookup_tables(aliases: Mapping[str, str] | None) -> tuple[dict[str, str], dict[str, str]]:
    if not aliases:
        return {}, {}
    exact = {str(k).strip(): str(v).strip() for k, v in aliases.items() if str(k).strip()}
    lowered = {k.casefold(): v for k, v in exact.items()}
    return exact, lowered


def resolve_tissue(name: Any, aliases: Mapping[str, str] | None) -> str:
    """将单个组织名解析为规范名（精确匹配优先，其次大小写不敏感匹配；空字符串原样返回）。"""
    text = "" if name is None else str(name).strip()
    if not text or not aliases:
        return text
    exact, lowered = _lookup_tables(aliases)
    if text in exact:
        return exact[text]
    return lowered.get(text.casefold(), text)


def remap_tissue_series(series: pd.Series, aliases: Mapping[str, str] | None) -> tuple[pd.Series, dict[str, str]]:
    """批量归一组织列，返回 ``(归一后的 Series, 实际发生改名的 {原名: 规范名})``。"""
    if series.empty or not aliases:
        return series, {}
    exact, lowered = _lookup_tables(aliases)
    applied: dict[str, str] = {}

    def _map(value: Any) -> Any:
        if pd.isna(value):
            return value
        text = str(value).strip()
        if not text:
            return text
        canonical = exact.get(text, lowered.get(text.casefold(), text))
        if canonical != text:
            applied.setdefault(text, canonical)
        return canonical

    mapped = series.map(_map)
    return mapped, applied


def remap_table_tissues(table: pd.DataFrame | None, aliases: Mapping[str, str] | None, *,
                        context: str, allow_merge: bool = False) -> tuple[pd.DataFrame | None, dict[str, str]]:
    """将结果表的 ``tissue`` 列按别名表归一。

    ``allow_merge=False``（默认，用于整合阶段读取已跑完的单组学差异/富集表）时，若同一张表内有多个不同的原组织名
    被映射到同一个规范组织名，则抛出 ``ValueError``，提示先重新运行 ``omicsweave run`` 在样本层面合并分析，
    避免把两份独立统计结果静默拼接或去重丢失。
    """
    if table is None or table.empty or not aliases or "tissue" not in table.columns:
        return table, {}
    raw = table["tissue"].astype(str)
    mapped, applied = remap_tissue_series(raw, aliases)
    if not applied:
        return table, {}
    if not allow_merge:
        by_canonical: dict[str, list[str]] = {}
        for orig, canon in dict.fromkeys(zip(raw.tolist(), mapped.tolist())):
            by_canonical.setdefault(canon, []).append(orig)
        collisions = {canon: raws for canon, raws in by_canonical.items() if len(raws) > 1}
        if collisions:
            detail = "；".join(f"{'、'.join(raws)} → {canon}" for canon, raws in sorted(collisions.items()))
            raise ValueError(
                f"{context} 中存在多个原组织被映射到同一规范组织（{detail}），"
                "不能在整合阶段直接合并已完成的统计表；请先重新运行 omicsweave run 将这些样本作为同一组织合并分析"
            )
    out = table.copy()
    out["tissue"] = mapped
    return out, applied


def format_alias_note(scope: str, applied: Mapping[str, str]) -> str:
    """将实际发生的别名改名格式化为运行提示；未发生改名时返回空字符串。"""
    if not applied:
        return ""
    pairs = "、".join(f"{src} → {dst}" for src, dst in sorted(applied.items()))
    return f"组织别名归一（{scope}）：{pairs}"


def select_tissue_file_rows(rows: pd.DataFrame, raw_key: str, canonical_key: str, *,
                            multi_key_canonical: bool, raw_tissues: Mapping[str, str]) -> pd.DataFrame:
    """为 ``input.tissue_files`` 中的某个条目挑选候选样本行。

    - 若 ``tissue_files`` 中有多个键都映射到同一个规范组织（如同时列了 ``sq_adipose`` 与 ``white_adipose``），
      优先按 ``samples.tsv`` 的原组织名匹配当前文件对应的那批样本（避免不同文件中同名数据列串扰）；
      若样本表中没有该原组织名（例如用户已在 ``samples.tsv`` 中提前写成了规范名），回退到规范组织名匹配。
    - 只有一个键对应到该规范组织时，直接匹配该规范组织下的全部样本。
    """
    canonical_rows = rows[rows["tissue"] == canonical_key]
    if not multi_key_canonical or canonical_rows.empty:
        return canonical_rows
    key_fold = raw_key.strip().casefold()
    matched = canonical_rows[
        canonical_rows["sample_id"].map(lambda sid: str(raw_tissues.get(sid, "")).strip().casefold() == key_fold)
    ]
    return matched if not matched.empty else canonical_rows
