#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MFA 求解结果的标准格式（v2 设计文档第 6.2 节）。

MFA 求解在外部软件中完成，OmicsWeave 读取其结果。标准格式由四张表组成，有两种载体：

- **Excel**（``omicsweave_mfa_xlsx``）：工作表 ``info`` / ``flux`` / ``mid`` / ``fits``，
  适合每个样本一个估计值（可附置信区间）或少量拟合；
- **Parquet 目录**（``mfa_tables``）：``info.json`` + ``flux.parquet`` / ``mid.parquet`` /
  ``fits.parquet``，列与 Excel 相同。保留全部拟合结果时行数通常超过 Excel 上限，应使用这种载体。

表结构：

- ``info``：``template_version``、``model_id``、``flux_units`` 必填；``software``、``tracer``、``notes`` 可选。
- ``flux``：每行一个（样本, 组织, 反应, 方向, 拟合）的通量值。``direction`` 为
  ``forward`` / ``reverse`` / ``net``：可逆反应分别给出正向与逆向（净通量 = 正向 − 逆向，
  按同一次拟合配对计算），求解软件直接给出净通量时用 ``net``。``fit_id`` 留空表示每个样本只有一个估计值。
- ``mid``：每行一个（样本, 代谢物, 同位素体）的实测与模拟 MID。
- ``fits``：每次拟合的 ``loss``，用于选出最优拟合。

样本与动物、分组的对应关系不写在这里，统一由实验的 ``samples.tsv``（``omics = mfa``）给出。
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict, ValidationError

from omicsweave.core.tissues import remap_tissue_series

from .dataset import ContractError

TEMPLATE_VERSION = 1
# Excel 单个工作表的行数上限（含表头）
EXCEL_MAX_ROWS = 1_048_576
DIRECTIONS = ("forward", "reverse", "net")
_TOL_FRACTION = 0.02
_TOL_SUM = 0.05


class MFAInfo(BaseModel):
    model_config = ConfigDict(extra="forbid")

    template_version: int = TEMPLATE_VERSION
    model_id: str
    flux_units: str
    software: str | None = None
    tracer: str | None = None
    notes: str | None = None


# 列名 -> (是否必填, 类型, 说明)；类型为 str / float / int
FLUX_COLUMNS: dict[str, tuple[bool, str, str]] = {
    "sample_id": (True, "str", "MFA 样本名，需与 samples.tsv 中 omics = mfa 的 sample_id 一致"),
    "tissue": (True, "str", "组织"),
    "reaction_id": (True, "str", "代谢网络模型中的反应 ID（不含方向后缀）"),
    "direction": (True, "str", "forward / reverse / net"),
    "fit_id": (False, "str", "拟合编号；留空表示每个样本只有一个估计值"),
    "value": (True, "float", "通量值，单位见 info.flux_units"),
    "ci_lower": (False, "float", "置信区间下限（单一估计值时可填）"),
    "ci_upper": (False, "float", "置信区间上限"),
}
MID_COLUMNS: dict[str, tuple[bool, str, str]] = {
    "sample_id": (True, "str", "MFA 样本名"),
    "tissue": (False, "str", "组织；多组织模型中区分同名代谢物"),
    "metabolite": (True, "str", "代谢物"),
    "isotopologue": (True, "int", "同位素体序号：0 表示 M+0，1 表示 M+1，依此类推"),
    "measured": (False, "float", "实测 MID 分数"),
    "simulated": (False, "float", "最优拟合的模拟 MID 分数"),
    "sd": (False, "float", "实测标准差"),
    "fit_id": (False, "str", "模拟值所对应的拟合编号"),
}
FITS_COLUMNS: dict[str, tuple[bool, str, str]] = {
    "sample_id": (True, "str", "MFA 样本名"),
    "fit_id": (True, "str", "拟合编号，与 flux.fit_id 对应"),
    "loss": (True, "float", "拟合残差，越小越好"),
}
INFO_FIELDS: dict[str, tuple[bool, str]] = {
    "template_version": (True, f"模板版本，当前为 {TEMPLATE_VERSION}"),
    "model_id": (True, "代谢网络模型标识，如 central_carbon_v1"),
    "flux_units": (True, "通量单位，如 nmol/min/g tissue"),
    "software": (False, "求解软件及版本"),
    "tracer": (False, "规范同位素示踪剂（U-13C / U-15N / U-13C-15N，见 omicsweave.mfa.tracers），如 U-13C-glucose"),
    "notes": (False, "备注"),
}
TABLES = {"flux": FLUX_COLUMNS, "mid": MID_COLUMNS, "fits": FITS_COLUMNS}


def normalize_id(value: Any) -> str:
    """Excel 会把 3 读成 3.0；编号统一成去掉多余小数的字符串，空值为空字符串。"""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value).strip()


def _normalize(table: pd.DataFrame, columns: dict[str, tuple[bool, str, str]]) -> pd.DataFrame:
    table = table.copy()
    table.columns = [str(c).strip() for c in table.columns]
    for name, (_, kind, _) in columns.items():
        if name not in table.columns:
            continue
        if kind == "str":
            table[name] = table[name].map(normalize_id).astype(object)
        else:
            table[name] = pd.to_numeric(table[name], errors="coerce")
    return table


@dataclass
class MFAResults:
    info: MFAInfo
    flux: pd.DataFrame
    mid: pd.DataFrame | None = None
    fits: pd.DataFrame | None = None

    def __post_init__(self) -> None:
        self.flux = _normalize(self.flux, FLUX_COLUMNS)
        if "fit_id" not in self.flux.columns:
            self.flux["fit_id"] = ""
        if self.mid is not None:
            self.mid = _normalize(self.mid, MID_COLUMNS)
        if self.fits is not None:
            self.fits = _normalize(self.fits, FITS_COLUMNS)

    @property
    def samples(self) -> list[str]:
        return sorted(set(self.flux["sample_id"]))

    def best_fit_ids(self) -> dict[str, str]:
        """每个样本 loss 最小的拟合编号；没有 fits 表的样本不出现在结果中。"""
        if self.fits is None or self.fits.empty:
            return {}
        ranked = self.fits.dropna(subset=["loss"]).sort_values(["sample_id", "loss"], kind="stable")
        return dict(ranked.groupby("sample_id")["fit_id"].first())

    def best_estimates(self) -> pd.DataFrame:
        """每个（样本, 组织, 反应, 方向）一个代表值。

        有 loss 时取最优拟合；有多次拟合但没有 loss 时取各次拟合的均值；只有一个估计值时取该值。
        与旧版 ``MFAAnalyser.load_flux_data`` 的取值规则一致。
        """
        keys = ["sample_id", "tissue", "reaction_id", "direction"]
        best = self.best_fit_ids()
        flux = self.flux
        with_best = flux["sample_id"].isin(best)
        picked = flux[with_best & (flux["fit_id"] == flux["sample_id"].map(best))]
        averaged = flux[~with_best].groupby(keys, as_index=False)["value"].mean()
        return pd.concat([picked[keys + ["value"]], averaged], ignore_index=True)

    def net_flux(self) -> pd.DataFrame:
        """按拟合配对计算净通量：同一样本、同一次拟合内正向 − 逆向。

        不可逆反应（没有 reverse）的净通量即正向值；``direction = net`` 的行原样保留。
        不能先分别汇总正、逆向再相减：那样配对的是不同拟合，结果可能偏差很大甚至反号。
        """
        keys = ["sample_id", "tissue", "reaction_id", "fit_id"]
        flux = self.flux
        given = flux[flux["direction"] == "net"][keys + ["value"]]
        forward = flux[flux["direction"] == "forward"].set_index(keys)["value"]
        reverse = flux[flux["direction"] == "reverse"].set_index(keys)["value"]
        net = forward.sub(reverse.reindex(forward.index).fillna(0.0))
        return pd.concat([net.rename("value").reset_index(), given], ignore_index=True)


def validate_mfa(results: MFAResults) -> tuple[list[str], list[str]]:
    """返回（错误, 警告）。错误表示数据无法使用，警告表示数据可能有问题。"""
    errors: list[str] = []
    warnings: list[str] = []

    for name, columns in TABLES.items():
        table = getattr(results, name)
        if table is None:
            continue
        missing = [c for c, (req, _, _) in columns.items() if req and c not in table.columns]
        if missing:
            errors.append(f"{name} 缺少必填列：{', '.join(missing)}")
    if errors:
        return errors, warnings

    flux = results.flux
    if flux.empty:
        errors.append("flux 表为空")
        return errors, warnings
    for column in ("sample_id", "tissue", "reaction_id"):
        blank = flux[column] == ""
        if blank.any():
            errors.append(f"flux 有 {int(blank.sum())} 行的 {column} 为空")
    bad_direction = set(flux["direction"]) - set(DIRECTIONS)
    if bad_direction:
        errors.append(f"flux.direction 取值无法识别：{', '.join(sorted(bad_direction))}（可选：{' / '.join(DIRECTIONS)}）")
    not_finite = ~np.isfinite(flux["value"].to_numpy(dtype=float))
    if not_finite.any():
        errors.append(f"flux.value 有 {int(not_finite.sum())} 个缺失或非数值")
    key = ["sample_id", "tissue", "reaction_id", "direction", "fit_id"]
    duplicated = int(flux.duplicated(key).sum())
    if duplicated:
        errors.append(f"flux 有 {duplicated} 行在 {' / '.join(key)} 上重复")

    reaction_key = ["sample_id", "tissue", "reaction_id"]
    directions = flux.groupby(reaction_key)["direction"].agg(set)
    mixed = directions[directions.map(lambda d: "net" in d and len(d) > 1)]
    if not mixed.empty:
        errors.append(f"{len(mixed)} 个反应同时给出了 net 与正/逆向，例如 {'/'.join(mixed.index[0])}")
    orphan = directions[directions.map(lambda d: "reverse" in d and "forward" not in d)]
    if not orphan.empty:
        errors.append(f"{len(orphan)} 个反应只有 reverse 没有 forward，例如 {'/'.join(orphan.index[0])}")

    if {"ci_lower", "ci_upper"} <= set(flux.columns):
        bounded = flux.dropna(subset=["ci_lower", "ci_upper"])
        outside = bounded[(bounded["value"] < bounded["ci_lower"]) | (bounded["value"] > bounded["ci_upper"])]
        if not outside.empty:
            warnings.append(f"flux 有 {len(outside)} 行的 value 不在 [ci_lower, ci_upper] 内")

    fitted = flux[flux["fit_id"] != ""]
    if not fitted.empty:
        if results.fits is None:
            warnings.append("flux 含多次拟合但没有 fits 表，将以各次拟合的均值作为代表值")
        else:
            known = set(zip(results.fits["sample_id"], results.fits["fit_id"]))
            used = set(zip(fitted["sample_id"], fitted["fit_id"]))
            if used - known:
                errors.append(f"有 {len(used - known)} 个（样本, 拟合）在 fits 表中找不到 loss")
            if results.fits.duplicated(["sample_id", "fit_id"]).any():
                errors.append("fits 表中（sample_id, fit_id）有重复")

    if results.mid is not None and not results.mid.empty:
        mid = results.mid
        if (mid["isotopologue"].isna() | (mid["isotopologue"] < 0)).any():
            errors.append("mid.isotopologue 必须是非负整数")
        value_cols = [c for c in ("measured", "simulated") if c in mid.columns]
        if not value_cols or mid[value_cols].isna().all(axis=1).any():
            errors.append("mid 每行至少需要 measured 或 simulated 之一")
        for column in value_cols:
            values = mid[column].dropna()
            if ((values < -_TOL_FRACTION) | (values > 1 + _TOL_FRACTION)).any():
                warnings.append(f"mid.{column} 有超出 [0, 1] 的值")
        if "measured" in mid.columns:
            group = ["sample_id", "metabolite"] + (["tissue"] if "tissue" in mid.columns else [])
            sums = mid.dropna(subset=["measured"]).groupby(group)["measured"].sum()
            off = sums[(sums - 1).abs() > _TOL_SUM]
            if not off.empty:
                warnings.append(f"{len(off)} 个（样本, 代谢物）的实测 MID 之和偏离 1 超过 {_TOL_SUM}")

    flux_samples = set(flux["sample_id"])
    for name in ("mid", "fits"):
        table = getattr(results, name)
        if table is not None and not table.empty:
            extra = set(table["sample_id"]) - flux_samples
            if extra:
                warnings.append(f"{name} 中有 {len(extra)} 个样本不在 flux 表中")
    return errors, warnings


def merge_mfa(parts: list[MFAResults]) -> MFAResults:
    """合并多组结果（如同一模型下不同示踪剂各自的拟合）。

    只有同一代谢网络模型、同一通量单位的结果才能放在一起分析；样本名不能在两组中重复。
    ``software`` / ``tracer`` / ``notes`` 合并为以分号分隔的不重复取值。
    """
    if not parts:
        raise ContractError("没有可合并的 MFA 结果")
    if len(parts) == 1:
        return parts[0]
    for field, label in (("model_id", "代谢网络模型（model_id）"), ("flux_units", "通量单位（flux_units）")):
        values = sorted({getattr(p.info, field) for p in parts})
        if len(values) > 1:
            raise ContractError(f"多组 MFA 结果的{label}不同：{' / '.join(values)}；只有同一模型、同一单位的结果才能合并分析")
    seen: set[str] = set()
    for part in parts:
        repeated = seen & set(part.samples)
        if repeated:
            raise ContractError(f"样本 {', '.join(sorted(repeated)[:5])} 出现在多组 MFA 结果中，样本名必须唯一")
        seen |= set(part.samples)

    def joined(field: str) -> str | None:
        values = dict.fromkeys(v.strip() for p in parts for v in (getattr(p.info, field) or "").split(";") if v.strip())
        return "; ".join(values) or None

    def stacked(name: str) -> pd.DataFrame | None:
        tables = [getattr(p, name) for p in parts if getattr(p, name) is not None]
        return pd.concat(tables, ignore_index=True) if tables else None

    info = MFAInfo(model_id=parts[0].info.model_id, flux_units=parts[0].info.flux_units,
                   software=joined("software"), tracer=joined("tracer"), notes=joined("notes"))
    return MFAResults(info=info, flux=stacked("flux"), mid=stacked("mid"), fits=stacked("fits"))


def apply_mfa_tissue_aliases(results: MFAResults,
                             aliases: dict[str, str] | None) -> tuple[MFAResults, dict[str, str]]:
    """按组织别名表归一 ``MFAResults`` 中 ``flux`` 与 ``mid`` 的 ``tissue`` 列。

    若同一样本内有多个原组织名映射到同一个规范组织并导致同一反应/MID 冲突，抛出 :class:`ContractError`。
    """
    if not aliases:
        return results, {}
    applied: dict[str, str] = {}
    flux = results.flux
    if not flux.empty and "tissue" in flux.columns:
        orig_tissue = flux["tissue"].astype(str)
        mapped, flux_applied = remap_tissue_series(orig_tissue, aliases)
        if flux_applied:
            applied.update(flux_applied)
            flux = flux.copy()
            flux["tissue"] = mapped
            key = ["sample_id", "tissue", "reaction_id", "direction", "fit_id"]
            if set(key) <= set(flux.columns):
                before_dup = results.flux.duplicated(key)
                after_dup = flux.duplicated(key, keep=False) & ~before_dup
                if after_dup.any():
                    samples = sorted(flux.loc[after_dup, "sample_id"].astype(str).unique())
                    tissues = sorted(flux.loc[after_dup, "tissue"].astype(str).unique())
                    raise ContractError(
                        f"组织别名归一后，同一 MFA 样本（{', '.join(samples[:5])}）在同一规范组织"
                        f"（{', '.join(tissues[:5])}）下出现重复通量行，请检查 tissue_aliases 是否合并了同一样本内的不同组织"
                    )

    mid = results.mid
    if mid is not None and not mid.empty and "tissue" in mid.columns:
        orig_mid = mid["tissue"].astype(str)
        mapped_mid, mid_applied = remap_tissue_series(orig_mid, aliases)
        if mid_applied:
            applied.update(mid_applied)
            mid = mid.copy()
            mid["tissue"] = mapped_mid
            mid_key = ["sample_id", "tissue", "metabolite", "isotopologue"]
            if set(mid_key) <= set(mid.columns):
                before_dup = results.mid.duplicated(mid_key)
                after_dup = mid.duplicated(mid_key, keep=False) & ~before_dup
                if after_dup.any():
                    samples = sorted(mid.loc[after_dup, "sample_id"].astype(str).unique())
                    tissues = sorted(mid.loc[after_dup, "tissue"].astype(str).unique())
                    raise ContractError(
                        f"组织别名归一后，同一 MFA 样本（{', '.join(samples[:5])}）在同一规范组织"
                        f"（{', '.join(tissues[:5])}）下出现重复 MID 行"
                    )

    if not applied:
        return results, {}
    return MFAResults(info=results.info, flux=flux, mid=mid, fits=results.fits), applied


def _parse_info(raw: dict[str, Any], source: str) -> MFAInfo:
    data = {k: v for k, v in raw.items() if v not in (None, "") and not (isinstance(v, float) and np.isnan(v))}
    try:
        return MFAInfo.model_validate(data)
    except ValidationError as exc:
        details = "; ".join(f"{'.'.join(map(str, e['loc']))}: {e['msg']}" for e in exc.errors())
        raise ContractError(f"{source} 的 info 不合法：{details}") from exc


def read_mfa(path: str | Path) -> MFAResults:
    """读取标准格式：``.xlsx`` 文件或 Parquet 目录。"""
    path = Path(path)
    if path.is_dir():
        info_file = path / "info.json"
        if not info_file.is_file():
            raise ContractError(f"{path} 中缺少 info.json")
        info = _parse_info(json.loads(info_file.read_text(encoding="utf-8")), str(path))
        tables = {name: pd.read_parquet(path / f"{name}.parquet") if (path / f"{name}.parquet").is_file() else None
                  for name in TABLES}
        if tables["flux"] is None:
            raise ContractError(f"{path} 中缺少 flux.parquet")
        return MFAResults(info=info, **tables)

    if not path.is_file():
        raise ContractError(f"找不到 {path}")
    sheets = pd.read_excel(path, sheet_name=None)
    if "info" not in sheets or "flux" not in sheets:
        raise ContractError(f"{path.name} 缺少 info 或 flux 工作表；旧版自研软件格式请用 omicsweave mfa convert 转换")
    info_sheet = sheets["info"]
    if not {"key", "value"} <= set(info_sheet.columns):
        raise ContractError(f"{path.name} 的 info 工作表需要 key 与 value 两列")
    info = _parse_info(dict(zip(info_sheet["key"].astype(str).str.strip(), info_sheet["value"])), path.name)
    optional = {name: sheets[name] if name in sheets and not sheets[name].empty else None for name in ("mid", "fits")}
    return MFAResults(info=info, flux=sheets["flux"], **optional)


def write_mfa(results: MFAResults, path: str | Path) -> Path:
    """写出标准格式：路径以 ``.xlsx`` 结尾写 Excel，否则写 Parquet 目录。"""
    errors, _ = validate_mfa(results)
    if errors:
        raise ContractError("；".join(errors))
    path = Path(path)
    tables = {name: getattr(results, name) for name in TABLES if getattr(results, name) is not None}

    if path.suffix.lower() == ".xlsx":
        too_long = [name for name, t in tables.items() if len(t) + 1 > EXCEL_MAX_ROWS]
        if too_long:
            raise ContractError(
                f"{'、'.join(too_long)} 表超过 Excel 的 {EXCEL_MAX_ROWS} 行上限；请改为输出 Parquet 目录（路径不要以 .xlsx 结尾）"
            )
        path.parent.mkdir(parents=True, exist_ok=True)
        info = pd.DataFrame(list(results.info.model_dump(exclude_none=True).items()), columns=["key", "value"])
        with pd.ExcelWriter(path, engine="openpyxl") as writer:
            info.to_excel(writer, sheet_name="info", index=False)
            for name, table in tables.items():
                table.to_excel(writer, sheet_name=name, index=False)
        return path

    path.mkdir(parents=True, exist_ok=True)
    (path / "info.json").write_text(results.info.model_dump_json(indent=2, exclude_none=True), encoding="utf-8")
    for name, table in tables.items():
        table.to_parquet(path / f"{name}.parquet", index=False)
    return path


def write_template(path: str | Path) -> Path:
    """写出空白 Excel 模板：数据工作表只有表头，另有“说明”工作表列出每一列的含义。"""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    guide_rows = [("info", key, "是" if req else "否", desc) for key, (req, desc) in INFO_FIELDS.items()]
    for name, columns in TABLES.items():
        guide_rows += [(name, col, "是" if req else "否", desc) for col, (req, _, desc) in columns.items()]
    guide = pd.DataFrame(guide_rows, columns=["工作表", "列 / 键", "必填", "说明"])
    info = pd.DataFrame({"key": list(INFO_FIELDS), "value": [TEMPLATE_VERSION] + [""] * (len(INFO_FIELDS) - 1)})
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        info.to_excel(writer, sheet_name="info", index=False)
        for name, columns in TABLES.items():
            pd.DataFrame(columns=list(columns)).to_excel(writer, sheet_name=name, index=False)
        guide.to_excel(writer, sheet_name="说明", index=False)
    return path
