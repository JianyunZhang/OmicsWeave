#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
图件的文件名与保存（v2 设计文档第 8.1、8.3 节）。

- 文件名 ``<图类型>_<组织>_<比较>.<格式>``，如 ``volcano_liver_HFD_vs_LFD.pdf``。
- 默认保存 PDF（矢量，文字可在排版软件中编辑）与 PNG（300 dpi）。
- 不写入创建时间，同一份结果表重画的图件逐字节一致，运行清单中的校验值可以直接比对。
"""
from __future__ import annotations

import re
from collections.abc import Sequence
from pathlib import Path

from matplotlib.figure import Figure

FIGURE_FORMATS: tuple[str, ...] = ("pdf", "png")
SAVE_DPI = 300
# 各格式中默认写入时间的元数据项；设为 None 即不写
_NO_DATE = {"pdf": {"CreationDate": None}, "svg": {"Date": None}}


def figure_stem(kind: str, *parts: object) -> str:
    """``<图类型>_<部分1>_<部分2>…``；文件名中不宜出现的字符改为下划线。"""
    names = [re.sub(r"[^\w.-]+", "_", str(p)).strip("_") or "NA" for p in (kind, *parts)]
    return "_".join(names)


def save_figure(fig: Figure, directory: Path, stem: str, formats: Sequence[str] = FIGURE_FORMATS,
                dpi: int = SAVE_DPI) -> list[Path]:
    """把 ``fig`` 按 ``formats`` 保存为 ``directory/stem.<格式>``，返回写出的文件。"""
    directory.mkdir(parents=True, exist_ok=True)
    paths = []
    for fmt in formats:
        path = directory / f"{stem}.{fmt}"
        fig.savefig(path, format=fmt, dpi=dpi, bbox_inches="tight", metadata=_NO_DATE.get(fmt))
        paths.append(path)
    return paths
