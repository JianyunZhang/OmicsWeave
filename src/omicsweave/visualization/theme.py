#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
统一的 matplotlib 外观：字体、字号、分辨率与保存参数。

各组学的绘图模块在导入时调用 :func:`apply_theme`。只有 ``style`` 与分辨率按模块保留原有取值，
其余设置全部统一。

字体：西文用第一个已安装的西文字体（Arial → Liberation Sans → DejaVu Sans），
中文字形回退到已安装的中文字体（Windows 的微软雅黑 / 黑体、macOS 的苹方、Linux 的 Noto / 思源 / 文泉驿）。
matplotlib 只有在 ``font.family`` 直接列出具体字体名时才逐字回退；写成 ``sans-serif`` 只会用其中第一个字体，
中文会显示为方框。未安装的字体不列入，以免产生字体缺失警告。
"""
from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
from matplotlib import font_manager

LATIN_FONTS = ("Arial", "Liberation Sans", "DejaVu Sans")
CJK_FONTS = (
    "Microsoft YaHei", "SimHei",                                       # Windows
    "PingFang SC", "Heiti SC",                                         # macOS
    "Noto Sans CJK SC", "Source Han Sans SC", "WenQuanYi Micro Hei",   # Linux
)

# 统一设置；pdf.fonttype 42 与 svg.fonttype none 使矢量图中的文字可在排版软件中编辑
BASE_RC: dict[str, Any] = {
    "axes.unicode_minus": False,
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
}


def font_family() -> list[str]:
    """已安装字体中的 ``font.family``：一个西文字体在前，随后是全部已安装的中文字体。"""
    installed = {f.name for f in font_manager.fontManager.ttflist}
    latin = next((name for name in LATIN_FONTS if name in installed), "DejaVu Sans")
    return [latin, *(name for name in CJK_FONTS if name in installed)]


def has_cjk_font() -> bool:
    """是否装有可显示中文的字体；没有时中文标签会显示为方框。"""
    return len(font_family()) > 1


def apply_theme(*, style: str | None = None, dpi: int = 150, save_dpi: int = 300,
                rc: dict[str, Any] | None = None) -> None:
    """应用统一主题。

    Parameters
    ----------
    style : str, optional
        先应用的 matplotlib 样式，如 ``seaborn-v0_8-whitegrid``；字体等设置随后覆盖样式中的同名项。
    dpi, save_dpi : int
        屏幕显示与保存时的分辨率；``savefig`` 显式传入 ``dpi`` 时以调用处为准。
    rc : dict, optional
        额外的 rcParams，最后应用。
    """
    if style:
        plt.style.use(style)
    plt.rcParams.update({**BASE_RC, "font.family": font_family(), "figure.dpi": dpi, "savefig.dpi": save_dpi,
                         **(rc or {})})
