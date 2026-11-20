#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
统一的 matplotlib 出版主题预设引擎（默认 ``nature_editorial``）。

设计原则（见 ``docs/arch_design/omicsweave_visualization_and_graphical_abstract.md`` 第 3 节）：
- 默认主风格为 ``nature_editorial``：纯白背景、开放式坐标轴（无上/右边框、无默认粗网格）、外向细刻度、
  色盲安全色板、300 dpi 保存与可编辑矢量字体（``pdf.fonttype=42``、``svg.fonttype="none"``）。
- 通用图件在 :func:`themed` 上下文中绘图：先复位 matplotlib 默认样式再套用指定预设，退出时完整恢复
  ``rcParams``，并持有可重入锁以保证多线程逐字节一致。
- 字体：西文优先已安装的 Arial / Helvetica / Liberation Sans / DejaVu Sans，中文字形回退到系统已安装字体；
  未安装的字体不列入，避免字形缺失警告。
"""
from __future__ import annotations

import threading
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.style
from matplotlib import font_manager

from .palettes import PALETTES, SemanticPalette, get_palette

LATIN_FONTS = ("Arial", "Helvetica", "Liberation Sans", "DejaVu Sans")
CJK_FONTS = (
    "Microsoft YaHei", "SimHei",                                       # Windows
    "PingFang SC", "Heiti SC",                                         # macOS
    "Noto Sans CJK SC", "Source Han Sans SC", "WenQuanYi Micro Hei",   # Linux
)

DEFAULT_PRESET = "nature_editorial"

# 底层不可变出版基础设置（所有预设与旧式 style 调用均继承）
BASE_RC: dict[str, Any] = {
    "axes.unicode_minus": False,
    "font.size": 9,
    "axes.titlesize": 10.5,
    "axes.labelsize": 9.5,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 8,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
}


@dataclass(frozen=True)
class StylePreset:
    """一套可复用的出版 / 展示主题预设。"""

    name: str
    description: str
    palette: SemanticPalette
    rc: dict[str, Any]


def _preset_rc(palette: SemanticPalette, *, title_weight: str = "bold", grid: bool = False) -> dict[str, Any]:
    return {
        "figure.facecolor": palette.surface,
        "savefig.facecolor": palette.surface,
        "axes.facecolor": palette.surface,
        "axes.edgecolor": palette.ink2,
        "axes.linewidth": 0.7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": grid,
        "grid.color": palette.hairline,
        "grid.linewidth": 0.5,
        "grid.alpha": 0.7,
        "axes.titleweight": title_weight,
        "axes.titlecolor": palette.ink,
        "axes.labelcolor": palette.ink,
        "text.color": palette.ink,
        "xtick.color": palette.ink2,
        "ytick.color": palette.ink2,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.width": 0.65,
        "ytick.major.width": 0.65,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "legend.frameon": False,
        "axes.prop_cycle": plt.cycler(color=list(palette.categorical)),
    }


PRESETS: dict[str, StylePreset] = {
    "nature_editorial": StylePreset(
        name="nature_editorial",
        description="顶刊极简定量风（纯白底、无冗余网格、色盲安全、适合论文正文与主图）",
        palette=PALETTES["nature_editorial"],
        rc=_preset_rc(PALETTES["nature_editorial"]),
    ),
    "metro_atlas": StylePreset(
        name="metro_atlas",
        description="瑞士网格拓扑风（高结构对比度，适合完整代谢通路总图与补充材料）",
        palette=PALETTES["metro_atlas"],
        rc=_preset_rc(PALETTES["metro_atlas"]),
    ),
    "slide_dark": StylePreset(
        name="slide_dark",
        description="暗色高对比演讲风（深石板灰底面，仅高亮显著信号，适合学术大会大屏汇报）",
        palette=PALETTES["slide_dark"],
        rc=_preset_rc(PALETTES["slide_dark"]),
    ),
}


def list_presets() -> tuple[str, ...]:
    """返回全部可用主题预设名。"""
    return tuple(PRESETS)


def get_preset(name: str = DEFAULT_PRESET) -> StylePreset:
    """按预设名获取主题预设对象；未知预设名抛出 ValueError。"""
    if name not in PRESETS:
        raise ValueError(f"未知的绘图主题预设 {name!r}（可用：{', '.join(PRESETS)}）")
    return PRESETS[name]


def font_family() -> list[str]:
    """已安装字体中的 ``font.family``：一个西文字体在前，随后是全部已安装的中文字体。"""
    installed = {f.name for f in font_manager.fontManager.ttflist}
    latin = next((name for name in LATIN_FONTS if name in installed), "DejaVu Sans")
    return [latin, *(name for name in CJK_FONTS if name in installed)]


def has_cjk_font() -> bool:
    """是否装有可显示中文的字体；没有时中文标签会显示为方框。"""
    return len(font_family()) > 1


def theme_rc(
    *,
    preset: str = DEFAULT_PRESET,
    dpi: int = 150,
    save_dpi: int = 300,
    rc: dict[str, Any] | None = None,
    include_preset_rc: bool = True,
) -> dict[str, Any]:
    """构建主题 ``rcParams``；``include_preset_rc=False`` 用于兼容显式指定 matplotlib ``style`` 的旧调用。"""
    spec = get_preset(preset)
    base = {**BASE_RC, **(spec.rc if include_preset_rc else {})}
    return {
        **base,
        "font.family": font_family(),
        "figure.dpi": dpi,
        "savefig.dpi": save_dpi,
        **(rc or {}),
    }


def apply_theme(
    *,
    style: str | None = None,
    preset: str = DEFAULT_PRESET,
    dpi: int = 150,
    save_dpi: int = 300,
    rc: dict[str, Any] | None = None,
) -> None:
    """应用统一主题到进程级 ``plt.rcParams``。

    Parameters
    ----------
    style : str, optional
        显式传入 matplotlib 内置样式（如 ``seaborn-v0_8-whitegrid``）时先套用该样式，再叠加字体/矢量/DPI
        基础设置，不覆写该样式的网格与框线；省略时默认套用 ``preset``（默认 ``nature_editorial``）。
    preset : str, default "nature_editorial"
        OmicsWeave 主题预设名（``nature_editorial`` / ``metro_atlas`` / ``slide_dark``）。
    dpi, save_dpi : int
        屏幕显示与保存分辨率（默认保存 300 dpi）。
    rc : dict, optional
        额外的 rcParams，最后覆写。
    """
    if style:
        plt.style.use(style)
        plt.rcParams.update(theme_rc(preset=preset, dpi=dpi, save_dpi=save_dpi, rc=rc, include_preset_rc=False))
    else:
        plt.rcParams.update(theme_rc(preset=preset, dpi=dpi, save_dpi=save_dpi, rc=rc, include_preset_rc=True))


# 可重入：通用图之间可以嵌套调用
_DRAW_LOCK = threading.RLock()


@contextmanager
def themed(
    *,
    preset: str = DEFAULT_PRESET,
    dpi: int = 150,
    save_dpi: int = 300,
    rc: dict[str, Any] | None = None,
) -> Iterator[StylePreset]:
    """在统一主题预设下绘图并保存，退出时恢复 ``rcParams``；线程安全且可重入。

    进入时先复位为 matplotlib ``default`` 样式再应用 ``preset``，彻底隔离旧模块的全局样式污染。
    上下文返回值 ``as style`` 为当前 :class:`StylePreset`，方便直接读取 ``style.palette``。
    """
    spec = get_preset(preset)
    with (
        _DRAW_LOCK,
        matplotlib.style.context("default"),
        matplotlib.rc_context(theme_rc(preset=preset, dpi=dpi, save_dpi=save_dpi, rc=rc, include_preset_rc=True)),
    ):
        yield spec


__all__ = [
    "BASE_RC",
    "CJK_FONTS",
    "DEFAULT_PRESET",
    "LATIN_FONTS",
    "PRESETS",
    "StylePreset",
    "apply_theme",
    "font_family",
    "get_palette",
    "get_preset",
    "has_cjk_font",
    "list_presets",
    "theme_rc",
    "themed",
]
