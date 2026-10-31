#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
旧版自研 MFA 软件导出格式与标准格式之间的转换。

旧格式（``flux_raw_data.xlsx`` / ``mid_raw_data.xlsx``）的约定，与
:meth:`omicsweave.mfa.mfa_analysis.MFAAnalyser.load_flux_data` 的读取逻辑一致：

- ``sheet_information`` 工作表：每个样本一行。``sheet_name`` 指向该样本的工作表，
  ``complete_result_name`` 为样本名，``condition`` 为分组；
- 通量工作表：每行一次拟合。``Index`` 为拟合编号，``loss`` 为残差；通量列名为
  ``组织--反应``，可逆反应的逆向列带 ``__R`` 后缀；``MIX:`` 开头的列是模拟 MID 观测量，不是通量；
- MID 工作表：``metabolite_name`` / ``mid_index`` / ``experiments``（实测）加每次拟合一列模拟值，
  ``metabolite_name == "loss"`` 的行给出每次拟合的残差。MID 只保留最优拟合的模拟值。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from omicsweave.contract.dataset import ContractError
from omicsweave.contract.mfa import MFAInfo, MFAResults, normalize_id

SHEET_INFO = "sheet_information"
REVERSE_SUFFIX = "__R"
UNKNOWN_TISSUE = "unknown"
_NON_FLUX = {"index", "loss"}


@dataclass
class LegacyConversion:
    results: MFAResults
    # 旧文件 sheet_information 中的样本 → 分组，供填写 subjects.tsv / samples.tsv 时参考
    conditions: dict[str, str] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)


def _sample_sheets(path: Path, notes: list[str]) -> list[tuple[str, str, pd.DataFrame]]:
    """按 sheet_information 的顺序返回 (样本名, 分组, 该样本的工作表)。"""
    if not path.is_file():
        raise ContractError(f"找不到 {path}")
    # 一次读入全部工作表：pd.ExcelFile 会一直占用文件句柄，在 Windows 上妨碍随后移动或删除文件
    sheets = pd.read_excel(path, sheet_name=None)
    if SHEET_INFO not in sheets:
        raise ContractError(f"{path.name} 中没有 {SHEET_INFO} 工作表，不是旧版自研软件格式")
    info = sheets[SHEET_INFO]
    if "sheet_name" not in info.columns:
        raise ContractError(f"{path.name} 的 {SHEET_INFO} 缺少 sheet_name 列")
    entries = []
    for _, row in info.iterrows():
        sheet = normalize_id(row["sheet_name"])
        if sheet not in sheets:
            notes.append(f"{path.name}: {SHEET_INFO} 中列出的工作表 {sheet} 不存在，已跳过")
            continue
        sample_id = normalize_id(row.get("complete_result_name", sheet)) or sheet
        entries.append((sample_id, normalize_id(row.get("condition", "")), sheets[sheet]))
    return entries


def detect_legacy_kind(path: str | Path) -> str | None:
    """旧版文件的类型：``flux`` / ``mid``；不是旧版格式时返回 None。"""
    with pd.ExcelFile(path) as xl:
        if SHEET_INFO not in xl.sheet_names:
            return None
        for sheet in xl.sheet_names:
            if sheet != SHEET_INFO:
                header = pd.read_excel(xl, sheet_name=sheet, nrows=0)
                return "mid" if "metabolite_name" in header.columns else "flux"
    return None


def split_flux_column(name: str) -> tuple[str, str, str]:
    """``组织--反应[__R]`` → (组织, 反应, forward/reverse)。没有组织前缀时组织记为 unknown。"""
    tissue, sep, rest = name.partition("--")
    if not sep:
        tissue, rest = UNKNOWN_TISSUE, name
    if rest.endswith(REVERSE_SUFFIX):
        return tissue, rest[: -len(REVERSE_SUFFIX)], "reverse"
    return tissue, rest, "forward"


def read_legacy_flux(path: Path, notes: list[str]) -> tuple[pd.DataFrame, pd.DataFrame | None, dict[str, str]]:
    flux_parts, fit_parts, conditions = [], [], {}
    unknown_tissue: set[str] = set()
    for sample_id, condition, df in _sample_sheets(path, notes):
        if df.empty:
            notes.append(f"{path.name}: 样本 {sample_id} 的工作表为空，已跳过")
            continue
        if "metabolite_name" in df.columns:
            raise ContractError(f"{path.name} 是 MID 文件，不是通量文件")
        conditions[sample_id] = condition
        lower = {str(c).lower(): c for c in df.columns}
        fit_ids = (df[lower["index"]].map(normalize_id) if "index" in lower
                   else pd.Series([str(i) for i in range(len(df))], index=df.index))
        if fit_ids.duplicated().any():
            raise ContractError(f"{path.name}: 样本 {sample_id} 的拟合编号（Index）有重复")
        flux_cols = [c for c in df.columns
                     if str(c).lower() not in _NON_FLUX and not str(c).startswith("MIX:")
                     and pd.api.types.is_numeric_dtype(df[c])]
        long = (df[flux_cols].assign(fit_id=fit_ids.values)
                .melt(id_vars="fit_id", var_name="column", value_name="value")
                .dropna(subset=["value"]))
        parsed = long["column"].astype(str).map(split_flux_column)
        long["tissue"] = parsed.str[0]
        long["reaction_id"] = parsed.str[1]
        long["direction"] = parsed.str[2]
        long["sample_id"] = sample_id
        unknown_tissue |= set(long.loc[long["tissue"] == UNKNOWN_TISSUE, "column"].astype(str))
        flux_parts.append(long[["sample_id", "tissue", "reaction_id", "direction", "fit_id", "value"]])
        if "loss" in lower:
            fit_parts.append(pd.DataFrame({
                "sample_id": sample_id,
                "fit_id": fit_ids.values,
                "loss": pd.to_numeric(df[lower["loss"]], errors="coerce").values,
            }))
    if not flux_parts:
        raise ContractError(f"{path.name} 中没有读到任何通量数据")
    if unknown_tissue:
        notes.append(f"{len(unknown_tissue)} 个通量列没有“组织--”前缀，组织记为 {UNKNOWN_TISSUE}："
                     f"{', '.join(sorted(unknown_tissue)[:5])}")
    fits = pd.concat(fit_parts, ignore_index=True) if fit_parts else None
    return pd.concat(flux_parts, ignore_index=True), fits, conditions


def read_legacy_mid(path: Path, notes: list[str]) -> tuple[pd.DataFrame, dict[str, str]]:
    parts, conditions = [], {}
    for sample_id, condition, df in _sample_sheets(path, notes):
        if df.empty or "metabolite_name" not in df.columns or "mid_index" not in df.columns:
            notes.append(f"{path.name}: 样本 {sample_id} 的工作表为空或缺少 metabolite_name / mid_index，已跳过")
            continue
        conditions[sample_id] = condition
        is_loss = df["metabolite_name"].astype(str) == "loss"
        rows = df[~is_loss]
        iter_cols = [c for c in df.columns
                     if c not in ("metabolite_name", "mid_index", "experiments")
                     and pd.api.types.is_numeric_dtype(df[c])]
        best_col = None
        if is_loss.any() and iter_cols:
            losses = pd.to_numeric(df.loc[is_loss, iter_cols].iloc[0], errors="coerce")
            if losses.notna().any():
                best_col = losses.idxmin()
        names = rows["metabolite_name"].astype(str).str.partition("--")
        has_tissue = names[1] == "--"
        parts.append(pd.DataFrame({
            "sample_id": sample_id,
            "tissue": names[0].where(has_tissue, ""),
            "metabolite": names[2].where(has_tissue, names[0]),
            "isotopologue": pd.to_numeric(rows["mid_index"], errors="coerce"),
            "measured": pd.to_numeric(rows["experiments"], errors="coerce") if "experiments" in rows else float("nan"),
            "simulated": pd.to_numeric(rows[best_col], errors="coerce") if best_col is not None else float("nan"),
            "fit_id": normalize_id(best_col) if best_col is not None else "",
        }))
        if best_col is None:
            notes.append(f"{path.name}: 样本 {sample_id} 没有 loss 行，模拟 MID 留空")
    if not parts:
        raise ContractError(f"{path.name} 中没有读到任何 MID 数据")
    return pd.concat(parts, ignore_index=True), conditions


def convert_legacy_mfa(flux_path: str | Path, mid_path: str | Path | None = None, *, model_id: str,
                       flux_units: str, software: str | None = None, tracer: str | None = None) -> LegacyConversion:
    """读取旧版自研软件导出的 Excel，转换为标准格式。"""
    notes: list[str] = []
    flux, fits, conditions = read_legacy_flux(Path(flux_path), notes)
    mid = None
    if mid_path is not None:
        mid, mid_conditions = read_legacy_mid(Path(mid_path), notes)
        for sample_id, condition in mid_conditions.items():
            conditions.setdefault(sample_id, condition)
    info = MFAInfo(model_id=model_id, flux_units=flux_units, software=software, tracer=tracer)
    return LegacyConversion(results=MFAResults(info=info, flux=flux, mid=mid, fits=fits),
                            conditions=conditions, notes=notes)


def _legacy_columns(table: pd.DataFrame) -> list[str]:
    """旧版通量列名：``组织--反应``，逆向加 ``__R``；没有组织前缀的列保持原样。"""
    return [(r.reaction_id if r.tissue == UNKNOWN_TISSUE else f"{r.tissue}--{r.reaction_id}")
            + (REVERSE_SUFFIX if r.direction == "reverse" else "")
            for r in table.itertuples()]


def legacy_flux_matrix(results: MFAResults) -> pd.DataFrame:
    """还原旧版 ``MFAAnalyser.flux_df``：样本 × ``组织--反应[__R]``，取每个样本的代表值。"""
    best = results.best_estimates()
    best["column"] = _legacy_columns(best)
    return best.pivot(index="sample_id", columns="column", values="value")


def legacy_all_fits(results: MFAResults) -> dict[str, pd.DataFrame]:
    """还原旧版 ``MFAAnalyser.flux_all_fits``：样本 → 每次拟合一行、每个通量一列（有 loss 时另加 loss 列）。

    行按拟合在结果中首次出现的顺序排列；旧版工作表中的 ``Index`` 与 ``MIX:`` 列不参与任何通量计算，不还原。
    """
    flux = results.flux.copy()
    flux["column"] = _legacy_columns(flux)
    losses = (results.fits.set_index(["sample_id", "fit_id"])["loss"]
              if results.fits is not None and not results.fits.empty else None)
    all_fits = {}
    for sample_id, part in flux.groupby("sample_id", sort=False):
        wide = part.pivot(index="fit_id", columns="column", values="value").reindex(pd.unique(part["fit_id"]))
        if losses is not None and sample_id in losses.index.get_level_values(0):
            wide["loss"] = losses.loc[sample_id].reindex(wide.index).to_numpy()
        all_fits[str(sample_id)] = wide.rename_axis(index=None, columns=None).reset_index(drop=True)
    return all_fits


def legacy_mid_matrices(results: MFAResults) -> tuple[pd.DataFrame, pd.DataFrame]:
    """还原旧版 ``MFAAnalyser.mid_df`` 与 ``mid_measured_df``：(metabolite_name, mid_index) × 样本。"""
    if results.mid is None:
        raise ContractError("没有 MID 数据")
    mid = results.mid.copy()
    mid["metabolite_name"] = [f"{t}--{m}" if t else m for t, m in zip(mid["tissue"], mid["metabolite"])]
    mid["mid_index"] = mid["isotopologue"].astype(int)
    index = ["metabolite_name", "mid_index"]
    return (mid.pivot(index=index, columns="sample_id", values="simulated"),
            mid.pivot(index=index, columns="sample_id", values="measured"))
