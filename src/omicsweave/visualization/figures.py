#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
图件的文件名、期刊栏宽几何与确定性保存（v2 设计文档第 8.1、8.3 节）。

- 文件名 ``<图类型>_<组织>_<比较>.<格式>``，如 ``volcano_liver_HFD_vs_LFD.pdf``。
- 默认遵循 ``cell_press`` 出版栏宽（85 / 114 / 174 mm，常规页高上限 200 mm），同时支持 ``profile="nature"``（89 / 120 / 183 mm）切换。
- 默认保存 PDF（嵌入式 Type 42 矢量字体）与 PNG（300 dpi），并支持无损 LZW TIFF、EPS 与 SVG。
- 不写入创建时间，同一份结果表重画的图件逐字节一致，运行清单中的校验值可以直接比对。
- 自动处理部分 fontTools 版本对中文字体 Type 42 子集化失败的边界情况（回退为 Type 3）。
"""
from __future__ import annotations

import os
import re
from collections.abc import Sequence
from io import BytesIO
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from .journal import (
    CELL_PRESS,
    DEFAULT_JOURNAL_PROFILE,
    JournalProfile,
    get_journal_profile,
)

FIGURE_FORMATS: tuple[str, ...] = ("pdf", "png")
SAVE_DPI = CELL_PRESS.default_dpi

# 默认期刊印刷栏宽常量（Cell Press 体系：85 / 114 / 174 mm，常规页高 200 mm；可通过 profile="nature" 切换）
SINGLE_COL_MM = CELL_PRESS.single_col_mm
ONE_HALF_COL_MM = CELL_PRESS.one_half_col_mm
DOUBLE_COL_MM = CELL_PRESS.double_col_mm
MAX_PAGE_H_MM = CELL_PRESS.max_page_h_mm
FULL_PAGE_H_MM = CELL_PRESS.full_page_h_mm

SINGLE_COL_IN = CELL_PRESS.single_col_in
ONE_HALF_COL_IN = CELL_PRESS.one_half_col_in
DOUBLE_COL_IN = CELL_PRESS.double_col_in
MAX_PAGE_H_IN = CELL_PRESS.max_page_h_in
FULL_PAGE_H_IN = CELL_PRESS.full_page_h_in

# 各格式中默认写入时间的元数据项；设为 None 即不写
_NO_DATE = {"pdf": {"CreationDate": None}, "svg": {"Date": None}}


def mm_to_in(mm: float) -> float:
    """毫米转英寸。"""
    return float(mm) / 25.4


def journal_figsize(
    columns: float = 2.0,
    aspect: float = 0.68,
    *,
    profile: str | JournalProfile = DEFAULT_JOURNAL_PROFILE,
    max_height_mm: float | None = None,
) -> tuple[float, float]:
    """按期刊栏宽计算标准图幅（英寸）。

    Parameters
    ----------
    columns : float
        ``1`` = 单栏（Cell Press 85 mm / Nature 89 mm）、``1.5`` = 1.5 栏（114 / 120 mm）、
        ``2`` = 双栏（174 / 183 mm）；其他正数按双栏比例换算。
    aspect : float
        高宽比 ``height / width``。
    profile : str or JournalProfile, default "cell_press"
        期刊出版规格（``"cell_press"`` | ``"nature"``）。
    max_height_mm : float, optional
        最大版心高度上限（省略时取 ``profile.max_page_h_mm``）。
    """
    spec = get_journal_profile(profile)
    w_mm = spec.width_mm(columns)
    cap_h = spec.max_page_h_mm if max_height_mm is None else float(max_height_mm)
    h_mm = min(w_mm * float(aspect), cap_h)
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
    bbox_inches: str | None = "tight",
    pad_inches: float | None = None,
    close: bool = False,
) -> list[Path]:
    """把 ``fig`` 按 ``formats`` 保存为 ``directory/stem.<格式>``，返回写出的文件。

    Parameters
    ----------
    bbox_inches : str or None, default "tight"
        默认裁剪外围白边；对于需要严格像素尺寸的图件（如 Cell Press 1200×1200 px 方形 GA），
        传入 ``bbox_inches=None`` 以精确锁死物理与像素维度。
    """
    directory.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    base_kwargs: dict[str, object] = {"dpi": dpi}
    rc_override: dict[str, object] = {}
    if bbox_inches is not None:
        base_kwargs["bbox_inches"] = bbox_inches
    else:
        # 显式覆写 rcParams["savefig.bbox"] = "tight"，精确锁死物理尺寸与像素（如 1200×1200 px）
        rc_override["savefig.bbox"] = "standard"
    if pad_inches is not None:
        base_kwargs["pad_inches"] = float(pad_inches)
    if facecolor is not None:
        base_kwargs["facecolor"] = facecolor

    for fmt in formats:
        fmt_clean = str(fmt).lstrip(".").lower()
        path = directory / f"{stem}.{fmt_clean}"
        kwargs = dict(base_kwargs)
        if fmt_clean in _NO_DATE:
            kwargs["metadata"] = _NO_DATE[fmt_clean]
        if fmt_clean in ("tiff", "tif"):
            kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        elif fmt_clean in ("jpg", "jpeg"):
            kwargs["pil_kwargs"] = {"quality": 95, "subsampling": 0}
        try:
            with plt.rc_context(rc_override):
                if os.name == "nt" and fmt_clean in ("tif", "tiff"):
                    # Pillow/libtiff 在部分 Windows 运行环境中直接写文件描述符会终止进程。
                    # 同一 LZW 编码器支持内存流，完成编码后再写入目标文件。
                    with BytesIO() as encoded:
                        fig.savefig(encoded, format=fmt_clean, **kwargs)
                        path.write_bytes(encoded.getbuffer())
                else:
                    fig.savefig(path, format=fmt_clean, **kwargs)
        except (AttributeError, KeyError, ValueError):
            if fmt_clean != "pdf":
                raise
            # 部分 fontTools 版本对中文字体（微软雅黑、黑体）做 Type 42 子集化时失败：改用 Type 3 嵌入
            with plt.rc_context({**rc_override, "pdf.fonttype": 3}):
                fig.savefig(path, format=fmt_clean, **kwargs)
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
    bbox_inches: str | None = "tight",
    pad_inches: float | None = None,
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
        bbox_inches=bbox_inches,
        pad_inches=pad_inches,
        close=close,
    )


__all__ = [
    "DOUBLE_COL_IN",
    "DOUBLE_COL_MM",
    "FIGURE_FORMATS",
    "FULL_PAGE_H_IN",
    "FULL_PAGE_H_MM",
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
