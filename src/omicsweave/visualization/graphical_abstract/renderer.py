#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Graphical Abstract 统一渲染入口。

专为小规模高层概念摘要（≤ 12 个核心枢纽）设计，与高密度定量绘图完全解耦：
- 支持两种版式策略（``GA_LAYOUTS``）：
  * ``cell_press_square``（**默认**）：严格遵循 Cell Press 官方规范——精确 **1200 × 1200 px @ 300 dpi**
    （4.0 × 4.0 in，``bbox_inches=None``），1:1 正方形单面板，全图 **Arial 12–16 pt**，沿用 ``metro_atlas`` 地铁设计风格，
    同步写出 ``<stem>_caption.txt`` 保存完整溯源与量化细则；
  * ``widescreen``：16:10 宽屏多枢纽概念摘要（适合综述插图、学术报告 PPT 与项目主页）。
- 支持 5 套渲染皮肤：``metro_atlas``（默认）/ ``economist`` / ``nature_editorial`` / ``editorial_ink`` / ``pixel_abstract``；
- 不写入时间戳，同一输入逐字节可重现。
"""
from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import matplotlib.pyplot as plt

from ..figures import DOUBLE_COL_IN, FIGURE_FORMATS, SAVE_DPI, save_stem
from ..journal import DEFAULT_JOURNAL_PROFILE, JournalProfile, get_journal_profile, write_caption_sidecar
from ..palettes import VALID_COLOR_MODES
from ..theme import themed
from .spec import AbstractSpec, validate_spec
from .square import render_cell_press_square
from .styles import GA_STYLES, STYLE_RENDERERS, ga_footer_height

GA_LAYOUTS: tuple[str, ...] = ("cell_press_square", "widescreen")
DEFAULT_GA_LAYOUT = "cell_press_square"

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
    journal: str | JournalProfile = DEFAULT_JOURNAL_PROFILE,
    layout: str | None = None,
    formats: Iterable[str] = FIGURE_FORMATS,
    width_in: float | None = None,
    dpi: int | None = None,
    write_caption: bool | None = None,
) -> list[Path]:
    """将 :class:`AbstractSpec` 渲染为指定风格与版式的 Graphical Abstract 图件。

    Parameters
    ----------
    spec : AbstractSpec
        摘要图规格（校验 ``1 <= len(hubs) <= 12``）。
    stem : Path or str
        输出路径（不含扩展名）。
    style : ``"metro_atlas"`` | ``"economist"`` | ``"nature_editorial"`` | ``"editorial_ink"`` | ``"pixel_abstract"``
        渲染皮肤，默认 ``metro_atlas``。
    color_mode : ``"semantic"`` | ``"pathway"``
        色彩语义模式（``"semantic"`` 默认红/蓝方向色；``"pathway"`` 低饱和通路分区色）。
    journal : str or JournalProfile, default "cell_press"
        期刊出版规格（``"cell_press"`` 默认 1200×1200 方形；``"nature"`` 默认宽屏）。
    layout : ``"cell_press_square"`` | ``"widescreen"`` | None
        显式指定版式；为 ``None`` 时继承 ``get_journal_profile(journal).ga_layout``（默认 ``"cell_press_square"``）。
    formats : iterable of str
        输出格式，默认 ``("pdf", "png")``（另支持 ``"tiff"``, ``"jpg"``, ``"eps"``, ``"svg"``）。
    width_in : float, optional
        宽屏模式下的输出物理宽度（英寸，省略时默认约双栏 1.18×）；方形模式下默认严格遵循 ``journal.ga_size_in``（4.0 in）。
    dpi : int, optional
        位图输出分辨率（默认 300 dpi）。
    write_caption : bool, optional
        是否写出配套 ``<stem>_caption.txt``（默认在 ``cell_press_square`` 且存在 takeaway/provenance 时自动写出）。
    """
    validate_spec(spec)
    if style not in STYLE_RENDERERS:
        raise ValueError(f"未知的 Graphical Abstract 风格 {style!r}（可用：{', '.join(GA_STYLES)}）")
    if color_mode not in VALID_COLOR_MODES:
        raise ValueError(f"未知的 color_mode {color_mode!r}（可用：{', '.join(VALID_COLOR_MODES)}）")

    jprof = get_journal_profile(journal)
    chosen_layout = str(layout or jprof.ga_layout).strip().lower()
    if chosen_layout not in GA_LAYOUTS:
        raise ValueError(f"未知的 Graphical Abstract 版式 {chosen_layout!r}（可用：{', '.join(GA_LAYOUTS)}）")

    out_dpi = int(dpi) if dpi is not None else (jprof.ga_dpi if chosen_layout == "cell_press_square" else SAVE_DPI)
    surface = _STYLE_SURFACES[style]
    preset = style if style in ("economist", "metro_atlas") else "nature_editorial"

    if chosen_layout == "cell_press_square":
        sq_w_in, sq_h_in = jprof.ga_size_in if width_in is None else (float(width_in), float(width_in))
        unit_in = sq_w_in / 100.0
        with themed(preset=preset, journal=jprof, dpi=150, save_dpi=out_dpi):
            fig = plt.figure(figsize=(sq_w_in, sq_h_in), facecolor=surface)
            ax = fig.add_axes([0.0, 0.0, 1.0, 1.0], facecolor=surface)
            render_cell_press_square(ax, spec, style=style, color_mode=color_mode, unit_in=unit_in)
            written = save_stem(
                fig, Path(stem), formats=tuple(formats), dpi=out_dpi,
                facecolor=surface, bbox_inches=None, close=True,
            )
        should_cap = bool(write_caption) if write_caption is not None else bool(spec.takeaway or spec.rules_note or spec.provenance)
        if should_cap:
            cap_lines = [s for s in (spec.takeaway, spec.rules_note, spec.provenance) if s]
            write_caption_sidecar(Path(stem), spec.title, cap_lines, panel_note=spec.subtitle)
        return written

    # Widescreen layout (legacy / reviews / slides)
    w_in = float(width_in) if width_in is not None else DOUBLE_COL_IN * 1.24
    w_units, h_units = spec.canvas
    header_units = 11.0
    footer_units = ga_footer_height(spec)
    total_h_units = h_units + header_units + footer_units
    height_in = w_in * (total_h_units / float(w_units))
    unit_in = w_in / float(w_units)

    with themed(preset=preset, journal=jprof, dpi=150, save_dpi=out_dpi):
        fig = plt.figure(figsize=(w_in, height_in), facecolor=surface)
        ax = fig.add_axes([0.0, 0.0, 1.0, 1.0], facecolor=surface)
        ax.set_xlim(0.0, w_units)
        ax.set_ylim(h_units + footer_units, -header_units)
        ax.set_aspect("equal")
        ax.axis("off")

        if style in ("metro_atlas", "economist"):
            STYLE_RENDERERS[style](ax, spec, color_mode=color_mode, unit_in=unit_in)
        else:
            STYLE_RENDERERS[style](ax, spec, unit_in=unit_in)
        return save_stem(fig, Path(stem), formats=tuple(formats), dpi=out_dpi, facecolor=surface, close=True)


__all__ = [
    "DEFAULT_GA_LAYOUT",
    "GA_LAYOUTS",
    "GA_STYLES",
    "VALID_COLOR_MODES",
    "render_graphical_abstract",
]
