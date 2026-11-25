#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
两个结果包的逐表对照（``omicsweave diff <运行 A> <运行 B>``）。

对两个结果包中的每张 Parquet 表（组学结果、数据集、QC 表、整合结果），按表的主键对齐行后逐格比较，判断是否逐值相同：

- 数值列在 ``atol + rtol·|B|`` 以内视为相同，两边都缺失视为相同；其他列按文本比较；
- 每张表给出状态（``identical`` / ``values_differ`` / ``rows_differ`` / ``columns_differ`` / ``only_in_a`` /
  ``only_in_b``）、行数、只在一边的行数、变化的格数与变化的列；
- 每个变化的数值列给出差值分布（变化格的 |B − A| 的中位数、95% 分位数与最大值）；
- 有等级或分类列（``evidence_grade``、``anchor_state``、``class``、``status``）的表，给出 A → B 的转移计数。

校验值相同的文件直接判为相同，不读取内容。非表格文件（图件、JSON 等）只按校验值统计相同与不同的个数。
``manifest.json`` 与 ``exports/`` 下的导出副本不参与对照。
主键见 :data:`TABLE_KEYS`；未登记的表用其索引（数据集），或用能唯一确定行的全部非数值列，都不行时按行号对齐。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from omicsweave.core.provenance import sha256_file

from .manifest import MANIFEST_FILE
from .tables import DIFFERENTIAL_KEY, ENRICHMENT_KEY

TABLE_KEYS: dict[str, tuple[str, ...]] = {
    "differential": DIFFERENTIAL_KEY,
    "derived_differential": DIFFERENTIAL_KEY,
    "differential_by_tracer": (*DIFFERENTIAL_KEY, "tracer"),
    "enrichment": ENRICHMENT_KEY,
    "flux_evidence_profile": ("comparison", "tissue", "flux_id"),
    "flux_evidence_table": ("comparison", "tissue", "flux_id", "layer"),
    "flux_pathway_profile": ("comparison", "tissue", "pathway"),
    "contrast_reliability": ("omics", "tissue", "contrast_id"),
    "reversal": ("pair_id", "omics", "tissue", "feature_id"),
    "reversal_summary": ("pair_id", "omics", "tissue", "feature_type"),
    "reaction_assessment": ("flux_id", "tissue", "contrast_id"),
    "axis_assessment": ("flux_id", "tissue", "contrast_id", "axis"),
}
CATEGORY_COLUMNS = ("evidence_grade", "anchor_state", "class", "status", "flux_claim", "state")
DEFAULT_ATOL, DEFAULT_RTOL = 1e-12, 1e-9


@dataclass
class TableDiff:
    path: str
    status: str
    key: str = ""
    rows_a: int = 0
    rows_b: int = 0
    rows_only_a: int = 0
    rows_only_b: int = 0
    cells_differ: int = 0
    max_abs_diff: float = float("nan")
    columns_only_a: list[str] = field(default_factory=list)
    columns_only_b: list[str] = field(default_factory=list)
    columns_differ: list[str] = field(default_factory=list)


@dataclass
class BundleDiff:
    a: Path
    b: Path
    tables: list[TableDiff]
    columns: pd.DataFrame        # 每个变化的列：path, column, kind, cells_differ, median/p95/max_abs_diff
    transitions: pd.DataFrame    # path, column, from, to, n
    other_files: dict[str, int]  # 非表格文件：identical / differ / only_in_a / only_in_b

    @property
    def identical(self) -> bool:
        return all(t.status == "identical" for t in self.tables) and not any(
            self.other_files.get(k, 0) for k in ("only_in_a", "only_in_b"))

    def table_frame(self) -> pd.DataFrame:
        rows = [{**t.__dict__, **{k: ";".join(getattr(t, k)) for k in ("columns_only_a", "columns_only_b",
                                                                       "columns_differ")}} for t in self.tables]
        return pd.DataFrame(rows, columns=list(TableDiff.__dataclass_fields__))


def _hashes(run_dir: Path) -> dict[str, str]:
    """结果包内全部文件（运行清单、交付元数据与导出副本除外）的实际校验值；不用清单中记录的值，结果包可能在运行后被改动。"""
    ignored_root = {MANIFEST_FILE, "index.json", ".index.tmp", ".keep", "summary.xlsx"}
    result: dict[str, str] = {}
    for p in sorted(run_dir.rglob("*")):
        if not p.is_file():
            continue
        rel = p.relative_to(run_dir)
        if (len(rel.parts) == 1 and rel.name in ignored_root) or rel.parts[0] == "exports":
            continue
        result[rel.as_posix()] = sha256_file(p)
    return result


def _key(path: str, a: pd.DataFrame, b: pd.DataFrame) -> list[str] | None:
    """行主键：登记的主键 → 非数值列 → None（按索引或行号对齐）。"""
    candidates = TABLE_KEYS.get(Path(path).stem)
    if candidates and all(c in a.columns and c in b.columns for c in candidates):
        return list(candidates)
    if not isinstance(a.index, pd.RangeIndex):
        return None
    text = [c for c in a.columns if c in b.columns and not pd.api.types.is_numeric_dtype(a[c])
            and not pd.api.types.is_bool_dtype(a[c])]
    if text and not a.duplicated(text).any() and not b.duplicated(text).any():
        return text
    return None


def _align(a: pd.DataFrame, b: pd.DataFrame, key: list[str] | None) -> tuple[pd.DataFrame, pd.DataFrame, int, int]:
    if key is None and isinstance(a.index, pd.RangeIndex) and isinstance(b.index, pd.RangeIndex):
        n = min(len(a), len(b))
        return a.iloc[:n], b.iloc[:n], len(a) - n, len(b) - n
    if key is not None:
        a = a.set_index([a[k].astype(str) for k in key]).drop(columns=key)
        b = b.set_index([b[k].astype(str) for k in key]).drop(columns=key)
    if a.index.has_duplicates or b.index.has_duplicates:
        a, b = a.reset_index(drop=True), b.reset_index(drop=True)
        return _align(a, b, None)
    common = a.index.intersection(b.index, sort=False)
    return a.loc[common], b.loc[common], len(a.index.difference(b.index)), len(b.index.difference(a.index))


def _column_diff(x: pd.Series, y: pd.Series, atol: float, rtol: float) -> tuple[np.ndarray, np.ndarray | None, str]:
    """返回 (变化的格, 变化格的 |差值| 或 None, 列类型)。"""
    numeric = (pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y)
               and not pd.api.types.is_bool_dtype(x) and not pd.api.types.is_bool_dtype(y))
    if numeric:
        u, v = x.to_numpy(dtype=float, na_value=np.nan), y.to_numpy(dtype=float, na_value=np.nan)
        both_nan = np.isnan(u) & np.isnan(v)
        with np.errstate(invalid="ignore"):
            close = np.isclose(u, v, atol=atol, rtol=rtol, equal_nan=False) | (u == v)
            changed = ~(close | both_nan)
            delta = np.abs(v - u)[changed]
        return changed, delta[np.isfinite(delta)], "numeric"
    # 缺失值单独比较：pandas 的 str 列转换后仍是 NaN，而 NaN != NaN
    missing_x, missing_y = x.isna().to_numpy(), y.isna().to_numpy()
    u = np.asarray(x.to_numpy(dtype=object), dtype=str)
    v = np.asarray(y.to_numpy(dtype=object), dtype=str)
    return (missing_x != missing_y) | (~missing_x & ~missing_y & (u != v)), None, "text"


def compare_frames(path: str, a: pd.DataFrame, b: pd.DataFrame, atol: float = DEFAULT_ATOL,
                   rtol: float = DEFAULT_RTOL) -> tuple[TableDiff, list[dict], list[dict]]:
    key = _key(path, a, b)
    result = TableDiff(path, "identical", key=";".join(key) if key else ("index" if not isinstance(
        a.index, pd.RangeIndex) else "row"), rows_a=len(a), rows_b=len(b))
    result.columns_only_a = [c for c in a.columns if c not in b.columns]
    result.columns_only_b = [c for c in b.columns if c not in a.columns]
    common_columns = [c for c in a.columns if c in b.columns and c not in (key or [])]
    left, right, result.rows_only_a, result.rows_only_b = _align(a[[*(key or []), *common_columns]],
                                                                 b[[*(key or []), *common_columns]], key)
    columns, transitions = [], []
    maxima = []
    for column in common_columns:
        changed, delta, kind = _column_diff(left[column], right[column], atol, rtol)
        n = int(changed.sum())
        if column in CATEGORY_COLUMNS:
            moved = pd.crosstab(left[column].astype(str).to_numpy(), right[column].astype(str).to_numpy())
            for (src, dst), count in moved.stack().items():
                if count and src != dst:
                    transitions.append({"path": path, "column": column, "from": src, "to": dst, "n": int(count)})
        if not n:
            continue
        result.cells_differ += n
        result.columns_differ.append(column)
        row = {"path": path, "column": column, "kind": kind, "cells_differ": n}
        if delta is not None and delta.size:
            row.update(median_abs_diff=float(np.median(delta)), p95_abs_diff=float(np.quantile(delta, 0.95)),
                       max_abs_diff=float(delta.max()))
            maxima.append(float(delta.max()))
        columns.append(row)
    result.max_abs_diff = max(maxima) if maxima else float("nan")
    if result.columns_only_a or result.columns_only_b:
        result.status = "columns_differ"
    elif result.rows_only_a or result.rows_only_b:
        result.status = "rows_differ"
    elif result.cells_differ:
        result.status = "values_differ"
    return result, columns, transitions


def diff_bundles(a: str | Path, b: str | Path, atol: float = DEFAULT_ATOL, rtol: float = DEFAULT_RTOL) -> BundleDiff:
    a, b = Path(a), Path(b)
    hashes_a, hashes_b = _hashes(a), _hashes(b)
    tables_a = {p for p in hashes_a if p.endswith(".parquet")}
    tables_b = {p for p in hashes_b if p.endswith(".parquet")}
    results, columns, transitions = [], [], []
    for path in sorted(tables_a | tables_b):
        if path not in tables_b:
            results.append(TableDiff(path, "only_in_a"))
            continue
        if path not in tables_a:
            results.append(TableDiff(path, "only_in_b"))
            continue
        if hashes_a.get(path) and hashes_a.get(path) == hashes_b.get(path):
            results.append(TableDiff(path, "identical", key="sha256"))
            continue
        result, cols, moves = compare_frames(path, pd.read_parquet(a / path), pd.read_parquet(b / path), atol, rtol)
        results.append(result)
        columns += cols
        transitions += moves

    other_a = {p: h for p, h in hashes_a.items() if not p.endswith(".parquet")}
    other_b = {p: h for p, h in hashes_b.items() if not p.endswith(".parquet")}
    shared = set(other_a) & set(other_b)
    other = {"identical": sum(other_a[p] == other_b[p] for p in shared),
             "differ": sum(other_a[p] != other_b[p] for p in shared),
             "only_in_a": len(set(other_a) - set(other_b)), "only_in_b": len(set(other_b) - set(other_a))}
    return BundleDiff(a, b, results,
                      pd.DataFrame(columns, columns=["path", "column", "kind", "cells_differ", "median_abs_diff",
                                                     "p95_abs_diff", "max_abs_diff"]),
                      pd.DataFrame(transitions, columns=["path", "column", "from", "to", "n"]), other)


def format_report(diff: BundleDiff, max_columns: int = 6) -> str:
    lines = [f"A：{diff.a}", f"B：{diff.b}", ""]
    counts = pd.Series([t.status for t in diff.tables]).value_counts()
    lines.append("表：" + "，".join(f"{status} {n}" for status, n in counts.items()) if len(counts) else "表：无")
    for t in diff.tables:
        if t.status == "identical":
            continue
        detail = [f"行 {t.rows_a}/{t.rows_b}"] if t.status not in ("only_in_a", "only_in_b") else []
        if t.rows_only_a or t.rows_only_b:
            detail.append(f"只在 A {t.rows_only_a} 行、只在 B {t.rows_only_b} 行")
        if t.cells_differ:
            shown = "、".join(t.columns_differ[:max_columns]) + ("…" if len(t.columns_differ) > max_columns else "")
            detail.append(f"{t.cells_differ} 格变化（{shown}）")
            if np.isfinite(t.max_abs_diff):
                detail.append(f"最大差值 {t.max_abs_diff:.3g}")
        if t.columns_only_a or t.columns_only_b:
            detail.append(f"列只在 A：{'、'.join(t.columns_only_a) or '无'}；只在 B：{'、'.join(t.columns_only_b) or '无'}")
        lines.append(f"  [{t.status}] {t.path}" + (f"：{'；'.join(detail)}" if detail else ""))
    if not diff.transitions.empty:
        lines.append("")
        lines.append("等级 / 分类的转移（A → B）：")
        for (path, column), part in diff.transitions.groupby(["path", "column"], sort=True):
            lines.append(f"  {path} · {column}")
            for row in part.sort_values("n", ascending=False).itertuples(index=False):
                lines.append(f"    {row[2]} → {row[3]}：{row[4]}")
    o = diff.other_files
    lines += ["", f"其他文件（按校验值）：相同 {o['identical']}，不同 {o['differ']}，只在 A {o['only_in_a']}，"
                  f"只在 B {o['only_in_b']}"]
    lines.append("结论：" + ("两个结果包的表逐值相同" if diff.identical else "有差异"))
    return "\n".join(lines)


def write_report(diff: BundleDiff, directory: str | Path) -> Path:
    """写出 tables.tsv、columns.tsv、transitions.tsv 与 report.txt。"""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    from .export import write_excel_tsv
    write_excel_tsv(diff.table_frame(), directory / "tables.tsv")
    write_excel_tsv(diff.columns, directory / "columns.tsv")
    write_excel_tsv(diff.transitions, directory / "transitions.tsv")
    (directory / "report.txt").write_text(format_report(diff), encoding="utf-8")
    return directory
