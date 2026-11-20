#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
图件的文件名、期刊栏宽几何与确定性保存（v2 设计文档第 8.1、8.3 节）。

- 文件名 ``<图类型>_<组织>_<比较>.<格式>``，如 ``volcano_liver_HFD_vs_LFD.pdf``。
- 默认保存 PDF（矢量，文字可在排版软件中编辑）与 PNG（300 dpi）。
- 不写入创建时间，同一份结果表重画的图件逐字节一致，运行清单中的校验值可以直接比对。
- 自动处理部分 fontTools 版本对中文字体 Type 42 子集化失败的边界情况（回退为 Type 3）。
"""
from __future__ import annotations

import re
from collections.abc import Sequence
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.figure import Figure

FIGURE_FORMATS: tuple[str, ...] = ("pdf", "png")
SAVE_DPI = 300

# 主流期刊印刷栏宽常量（Nature / Cell / Science 体系，单位 mm）
SINGLE_COL_MM = 89.0
ONE_HALF_COL_MM = 120.0
DOUBLE_COL_MM = 183.0
MAX_PAGE_H_MM = 247.0

SINGLE_COL_IN = SINGLE_COL_MM / 25.4
ONE_HALF_COL_IN = ONE_HALF_COL_MM / 25.4
DOUBLE_COL_IN = DOUBLE_COL_MM / 25.4
MAX_PAGE_H_IN = MAX_PAGE_H_MM / 25.4

# 各格式中默认写入时间的元数据项；设为 None 即不写
_NO_DATE = {"pdf": {"CreationDate": None}, "svg": {"Date": None}}


def mm_to_in(mm: float) -> float:
    """毫米转英寸。"""
    return float(mm) / 25.4


def journal_figsize(columns: float = 2.0, aspect: float = 0.68, *, max_height_mm: float = MAX_PAGE_H_MM) -> tuple[float, float]:
    """按期刊栏宽计算标准图幅（英寸）。

    Parameters
    ----------
    columns : float
        ``1`` = 单栏 89 mm、``1.5`` = 半通栏 120 mm、``2`` = 双栏 183 mm；其他正数按双栏比例换算。
    aspect : float
        高宽比 ``height / width``。
    max_height_mm : float
        最大版心高度上限（默认 247 mm）。
    """
    if columns == 1 or columns == 1.0:
        w_mm = SINGLE_COL_MM
    elif columns == 1.5:
        w_mm = ONE_HALF_COL_MM
    elif columns == 2 or columns == 2.0:
        w_mm = DOUBLE_COL_MM
    else:
        w_mm = DOUBLE_COL_MM * (float(columns) / 2.0)
    h_mm = min(w_mm * float(aspect), float(max_height_mm))
    return mm_to_in(w_mm), mm_to_in(h_mm)


def figure_stem(kind: str, *parts: object) -> str:
    """``<图类型>_<部分1>_<部分2>…``；文件名中不宜出现的字符改为下划线。"""
    names = [re.sub(r"[^\w.-]+", "_", str(p)).strip("_") or "NA" for p in (kind, *parts)]
    return "_".join(names)


def save_figure(
    fig: Figure,
    directory: Path,
    stem: str,
    formats: Sequence[str] = FIGURE_FORMATS,
    dpi: int = SAVE_DPI,
    *,
    facecolor: str | None = None,
    close: bool = False,
) -> list[Path]:
    """把 ``fig`` 按 ``formats`` 保存为 ``directory/stem.<格式>``，返回写出的文件。"""
    directory.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    kwargs: dict[str, object] = {"dpi": dpi, "bbox_inches": "tight"}
    if facecolor is not None:
        kwargs["facecolor"] = facecolor
    for fmt in formats:
        fmt_clean = str(fmt).lstrip(".").lower()
        path = directory / f"{stem}.{fmt_clean}"
        try:
            fig.savefig(path, format=fmt_clean, metadata=_NO_DATE.get(fmt_clean), **kwargs)
        except (AttributeError, KeyError, ValueError):
            if fmt_clean != "pdf":
                raise
            # 部分 fontTools 版本对中文字体（微软雅黑、黑体）做 Type 42 子集化时失败：改用 Type 3 嵌入
            with plt.rc_context({"pdf.fonttype": 3}):
                fig.savefig(path, format=fmt_clean, metadata=_NO_DATE.get(fmt_clean), **kwargs)
        paths.append(path)
    if close:
        plt.close(fig)
    return paths


def save_stem(
    fig: Figure,
    stem_path: Path,
    formats: Sequence[str] = FIGURE_FORMATS,
    dpi: int = SAVE_DPI,
    *,
    facecolor: str | None = None,
    close: bool = True,
) -> list[Path]:
    """``save_figure`` 的路径便捷封装：直接传不含扩展名的完整 ``Path``，默认保存后关闭图件。"""
    return save_figure(
        fig,
        stem_path.parent,
        stem_path.name,
        formats=formats,
        dpi=dpi,
        facecolor=facecolor,
        close=close,
    )


__all__ = [
    "DOUBLE_COL_IN",
    "DOUBLE_COL_MM",
    "FIGURE_FORMATS",
    "MAX_PAGE_H_IN",
    "MAX_PAGE_H_MM",
    "ONE_HALF_COL_IN",
    "ONE_HALF_COL_MM",
    "SAVE_DPI",
    "SINGLE_COL_IN",
    "SINGLE_COL_MM",
    "figure_stem",
    "journal_figsize",
    "mm_to_in",
    "save_figure",
    "save_stem",
]
