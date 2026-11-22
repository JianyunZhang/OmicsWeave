#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Graphical Abstract 统一渲染入口。

专为小规模高层概念摘要（≤ 12 个核心枢纽）设计，与高密度定量绘图完全解耦：
- 支持 ``editorial_ink``（科学手绘插画）、``pixel_abstract``（现代像素风）、``nature_editorial``（极简机制图）、
  ``metro_atlas``（Beck 代谢地铁换乘图，支持 ``color_mode="semantic"`` 默认正文统一红蓝语义色或 ``color_mode="pathway"`` 莫兰迪通路分区色）；
- 默认按双栏宽（约 183 mm / 7.2 in，适配 Cell / Nature Graphical Abstract 比例）渲染并导出 300 dpi + 矢量 PDF/SVG；
- 不写入时间戳，同一输入逐字节可重现。
"""
from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import matplotlib.pyplot as plt

from ..figures import DOUBLE_COL_IN, FIGURE_FORMATS, SAVE_DPI, save_stem
from ..palettes import VALID_COLOR_MODES
from ..theme import themed
from .spec import AbstractSpec, validate_spec
from .styles import GA_STYLES, STYLE_RENDERERS, ga_footer_height

_STYLE_SURFACES = {
    "economist": "#EBEBDF",
    "nature_editorial": "#FFFFFF",
    "editorial_ink": "#FAF8F2",
    "pixel_abstract": "#F4EFE4",
    "metro_atlas": "#FFFFFF",
}


def render_graphical_abstract(
    spec: AbstractSpec,
    stem: Path | str,
    *,
    style: str = "metro_atlas",
    color_mode: str = "semantic",
    formats: Iterable[str] = FIGURE_FORMATS,
    width_in: float = DOUBLE_COL_IN * 1.18,
    dpi: int = SAVE_DPI,
) -> list[Path]:
    """将 :class:`AbstractSpec` 渲染为指定风格的 Graphical Abstract 图件。

    Parameters
    ----------
    spec : AbstractSpec
        摘要图规格（校验 ``1 <= len(hubs) <= 12``）。
    stem : Path or str
        输出路径（不含扩展名）。
    style : ``"editorial_ink"`` | ``"pixel_abstract"`` | ``"nature_editorial"`` | ``"metro_atlas"``
        渲染皮肤，默认 ``editorial_ink``（手绘科学插画风）。
    color_mode : ``"semantic"`` | ``"pathway"``
        色彩语义模式（目前主要用于 ``metro_atlas``）：``"semantic"``（默认）沿用正文统一的红/蓝/灰方向色；
        ``"pathway"`` 采用克制的低饱和莫兰迪通路分区色。
    formats : iterable of str
        输出格式，默认 ``("pdf", "png")``。
    width_in : float
        输出物理宽度（英寸，默认约 8.5 in，可直接用于 Cell/Nature 图文摘要或 PPT 宽屏）。
    dpi : int
        位图输出分辨率（默认 300 dpi）。
    """
    validate_spec(spec)
    if style not in STYLE_RENDERERS:
        raise ValueError(f"未知的 Graphical Abstract 风格 {style!r}（可用：{', '.join(GA_STYLES)}）")
    if color_mode not in VALID_COLOR_MODES:
        raise ValueError(f"未知的 color_mode {color_mode!r}（可用：{', '.join(VALID_COLOR_MODES)}）")

    w_units, h_units = spec.canvas
    header_units = 11.0
    footer_units = ga_footer_height(spec)
    total_h_units = h_units + header_units + footer_units
    height_in = float(width_in) * (total_h_units / float(w_units))
    unit_in = float(width_in) / float(w_units)
    surface = _STYLE_SURFACES[style]

    preset = style if style in ("economist", "metro_atlas") else "nature_editorial"
    with themed(preset=preset, dpi=150, save_dpi=dpi):
        fig = plt.figure(figsize=(float(width_in), height_in), facecolor=surface)
        ax = fig.add_axes([0.0, 0.0, 1.0, 1.0], facecolor=surface)
        ax.set_xlim(0.0, w_units)
        ax.set_ylim(h_units + footer_units, -header_units)
        ax.set_aspect("equal")
        ax.axis("off")

        if style in ("metro_atlas", "economist"):
            STYLE_RENDERERS[style](ax, spec, color_mode=color_mode, unit_in=unit_in)
        else:
            STYLE_RENDERERS[style](ax, spec, unit_in=unit_in)
        return save_stem(fig, Path(stem), formats=tuple(formats), dpi=dpi, facecolor=surface, close=True)


__all__ = ["GA_STYLES", "VALID_COLOR_MODES", "render_graphical_abstract"]
