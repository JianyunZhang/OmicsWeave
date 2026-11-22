#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Graphical Abstract 专用四套渲染皮肤（严格限定于 ≤ 12 枢纽的高层概念摘要图）。

1. ``nature_editorial``：顶刊极简机制图（纯白底、精密矢量胶囊、色盲安全红蓝流向与无彩色阶梯徽章）；
2. ``editorial_ink``：手绘科学插画风（暖宣纸底、有机双层线粒体嵴/波浪血流、矿物水彩淡晕染 + 精密矢量文字数据层）；
3. ``pixel_abstract``：现代 16-bit 像素概念风（纯矢量硬边像素精灵图标、像素阶梯卡框、正交像素管线与克制矿物调色板）；
4. ``metro_atlas``：Beck 代谢地铁换乘图风（八向 0°/45°/90° 圆角轨道、内嵌方向 chevron、实心 L2+ 主轨 vs 空心白芯 L1 双轨、
   圆角跑道闭环循环与双层圆环换乘站；默认 ``color_mode="semantic"`` 统一采用正文红蓝语义色，可选 ``color_mode="pathway"`` 莫兰迪通路分区色）。
"""
from __future__ import annotations

import math
import textwrap
from typing import Callable

import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, PathPatch, Polygon, Rectangle
from matplotlib.path import Path as MplPath

from ..glyphs import GLYPH_STYLES, draw_glyph
from ..metro import (
    draw_metro_coupled_halo,
    draw_metro_loop,
    draw_metro_station,
    draw_metro_track,
    draw_metro_zone,
    fillet_polyline,
    polyline_length,
    route_octilinear,
    sample_polyline_at,
)
from ..palettes import (
    ECONOMIST,
    ECONOMIST_RED,
    METRO_ATLAS,
    NATURE_EDITORIAL,
    VALID_COLOR_MODES,
    pathway_family_color,
    tint_color,
)
from ..primitives import (
    draw_economist_masthead,
    draw_quad_arrow,
    draw_stealth_head,
    draw_tier_badge,
    estimate_text_width,
    tier_badge_text,
    tier_badge_width,
)
from .spec import AbstractHub, AbstractSpec

GA_STYLES: tuple[str, ...] = (*GLYPH_STYLES, "metro_atlas")

# ---------------------------------------------------------------------------
# 共享布局与配色
# ---------------------------------------------------------------------------

HUB_W = 22.8
HUB_H = 11.2


def _flow_color(direction: int, tier: str, up: str, down: str, muted: str) -> str:
    if direction > 0:
        return up
    if direction < 0:
        return down
    return muted


def _hub_border_and_tint(hub: AbstractHub, up: str, down: str, neutral: str, base: str):
    if hub.direction > 0:
        col = up
        amt = 0.09 if hub.evidence_tier in ("L3", "L2") else 0.055
    elif hub.direction < 0:
        col = down
        amt = 0.09 if hub.evidence_tier in ("L3", "L2") else 0.055
    else:
        col = neutral
        amt = 0.03
    strong = hub.evidence_tier in ("L3", "L2")
    dashed = hub.evidence_tier == "L1"
    return col, tint_color(col, amt, base=base), strong, dashed


def _box_exit_trim(dx: float, dy: float, half_w: float, half_h: float, pad: float = 0.8) -> float:
    """计算从矩形中心沿 (dx, dy) 射出到矩形边框的距离。"""
    norm = math.hypot(dx, dy) or 1.0
    ux, uy = abs(dx) / norm, abs(dy) / norm
    tx = half_w / ux if ux > 1e-6 else math.inf
    ty = half_h / uy if uy > 1e-6 else math.inf
    return min(tx, ty) + pad


def _wrap_rules_lines(spec: AbstractSpec) -> list[str]:
    if not spec.rules_note:
        return []
    wrap_cols = max(84, int((spec.canvas[0] - 8.4) * 0.94))
    return textwrap.wrap(spec.rules_note, wrap_cols)[:3]


def ga_footer_height(spec: AbstractSpec) -> float:
    """计算 Graphical Abstract 底栏精确所需画布高度（量化图例规则卡 + Takeaway + Provenance）。"""
    rules_lines = _wrap_rules_lines(spec)
    legend_box_h = (3.8 + len(rules_lines) * 2.15 + 0.6) if rules_lines else 4.2
    total = 2.0 + legend_box_h + 1.4
    if spec.takeaway:
        total += 5.0 + 1.4
    if spec.provenance:
        total += 3.2
    return total + 1.2


def _draw_ga_bottom_panels(
    ax,
    spec: AbstractSpec,
    *,
    surface: str,
    ink: str,
    ink2: str,
    muted: str,
    card_bg: str,
    card_edge: str,
    up: str,
    preset: str = "nature_editorial",
    unit_in: float = 0.062,
    pixel_mode: bool = False,
    metro_mode: bool = False,
) -> None:
    """绘制所有 Graphical Abstract 共用的量化细则与徽章图例卡 + Takeaway + Provenance。"""
    W, H = spec.canvas
    is_econ = preset == "economist"
    pal = ECONOMIST if is_econ else (METRO_ATLAS if metro_mode else NATURE_EDITORIAL)
    foot_y = H + 2.0

    # 1. 徽章、Linked flux 与聚合量化细则图例卡
    rules_lines = _wrap_rules_lines(spec)
    legend_h = (3.8 + len(rules_lines) * 2.15 + 0.6) if rules_lines else 4.2
    if pixel_mode:
        ax.add_patch(Rectangle((2.2, foot_y), W - 4.4, legend_h, facecolor=card_bg, edgecolor=ink, lw=1.2, zorder=2))
    elif is_econ:
        ax.add_patch(Rectangle((2.2, foot_y), W - 4.4, legend_h, facecolor=card_bg, edgecolor=card_edge, lw=0.65, zorder=2))
        ax.plot([2.2, W - 2.2], [foot_y, foot_y], color=ink, lw=0.95, solid_capstyle="butt", zorder=2.2)
    else:
        ax.add_patch(FancyBboxPatch(
            (2.2, foot_y), W - 4.4, legend_h,
            boxstyle="round,pad=0,rounding_size=1.2",
            facecolor=card_bg, edgecolor=card_edge, lw=0.72, zorder=2,
        ))

    # Row 1: 实景三点证据徽章 + Linked flux 轨道样块
    ry = foot_y + 2.05
    cx = 4.0
    lbl_fs = 5.25
    b_sc = 1.20
    for lvl, desc in (
        ("L3", "Multi-tracer"),
        ("L2", "MFA flux (F sig./sug.)"),
        ("L1", "Enzyme only (C sig./sug.)"),
    ):
        if pixel_mode:
            bw = _draw_pixel_tier_meter(
                ax, cx, ry, lvl, ink=ink, up_col="#D99B26" if lvl != "L1" else ink,
                unit_in=unit_in, zorder=3.2,
            )
        else:
            bw = draw_tier_badge(
                ax, cx, ry, lvl, zorder=3.2, scale=b_sc, preset=preset, unit_in=unit_in,
            )
        ax.text(cx + bw + 0.75, ry, desc, fontsize=lbl_fs, color=ink2, ha="left", va="center", zorder=3.2)
        cx += bw + 0.75 + estimate_text_width(desc, lbl_fs, unit_in) + 2.5

    # Indigo linked-flux conduit sample (explicitly clarifying shared 1 DOF, NOT extra confidence)
    seg_w = 4.2
    linked_txt = "Linked flux (|r|≥1−1e-6; 1 BH DOF, not extra conf.)"
    if cx + seg_w + 1.0 + estimate_text_width(linked_txt, lbl_fs, unit_in) <= W - 3.6:
        if metro_mode:
            draw_metro_coupled_halo(
                ax, [[(cx, ry), (cx + seg_w, ry)]],
                color=pal.coupled_track, edge_color=pal.coupled_edge, width=7.2, zorder=3.0,
            )
            draw_metro_track(
                ax, [(cx, ry), (cx + seg_w, ry)], up,
                width=2.8, hollow=False, surface=card_bg,
                chevrons=True, chevron_spacing=9.0, unit_in=unit_in,
                corner_radius=0.0, halo=False, zorder=3.2,
            )
        else:
            from ..primitives import draw_coupled_track
            draw_coupled_track(
                ax, [{"p0": (cx, ry), "c": (cx + seg_w * 0.5, ry), "p1": (cx + seg_w, ry), "trim0": 0.0, "trim1": 0.0}],
                color=pal.coupled_track, edge_color=pal.coupled_edge, width=6.4, zorder=3.0,
            )
            ax.plot([cx, cx + seg_w], [ry, ry], color=up, lw=1.6, zorder=3.2)
        ax.text(cx + seg_w + 0.85, ry, linked_txt, fontsize=lbl_fs, color=pal.coupled_ink, ha="left", va="center", zorder=3.2)

    # Row 2+: 自动聚合量化统计规则细则
    if rules_lines:
        ax.plot([3.8, W - 3.8], [foot_y + 3.55, foot_y + 3.55], color=card_edge, lw=0.45, alpha=0.75, zorder=2.5)
        for idx, rline in enumerate(rules_lines):
            ax.text(
                4.0, foot_y + 4.95 + idx * 2.15, rline,
                fontsize=5.2, color=ink2, ha="left", va="center", zorder=3.0,
            )

    foot_y += legend_h + 1.4

    # 2. Takeaway 摘要卡
    if spec.takeaway:
        if pixel_mode:
            ax.add_patch(Rectangle((2.2, foot_y), W - 4.4, 5.0, facecolor=card_bg, edgecolor=ink, lw=1.3, zorder=2))
        elif is_econ:
            ax.add_patch(Rectangle((2.2, foot_y), W - 4.4, 5.0, facecolor=card_bg, edgecolor=card_edge, lw=0.65, zorder=2))
            ax.add_patch(Rectangle((2.2, foot_y), 0.75, 5.0, facecolor=ECONOMIST_RED, edgecolor="none", zorder=2.3))
        else:
            ax.add_patch(FancyBboxPatch(
                (2.2, foot_y), W - 4.4, 5.0,
                boxstyle="round,pad=0,rounding_size=1.2",
                facecolor=card_bg, edgecolor=card_edge, lw=0.75, zorder=2,
            ))
        ax.text(4.0, foot_y + 2.5, spec.takeaway, fontsize=6.6, fontweight="medium", color=ink, ha="left", va="center", zorder=3)
        foot_y += 6.4

    # 3. Provenance 来源说明
    if spec.provenance:
        ax.text(2.4, foot_y + 0.4, spec.provenance, fontsize=5.6, color=muted, ha="left", va="top", zorder=3)


def _draw_footer_and_header(
    ax,
    spec: AbstractSpec,
    *,
    surface: str,
    ink: str,
    ink2: str,
    muted: str,
    card_bg: str,
    card_edge: str,
    up: str,
    down: str,
    pixel_mode: bool = False,
    unit_in: float = 0.062,
) -> None:
    W, _H = spec.canvas
    # 顶栏标题
    ax.text(2.2, -7.8, spec.title, fontsize=11.0, fontweight="bold", color=ink, ha="left", va="center")
    if spec.subtitle:
        ax.text(2.2, -3.5, spec.subtitle, fontsize=7.2, color=muted, ha="left", va="center")

    # 右上角极简方向与证据类型图例（与 metro_atlas 语义一致）
    lx = W - 2.4
    for label, col, shape in reversed([
        ("Up (Flux)", up, "up"),
        ("Down", down, "down"),
        ("Enz-pred", ink2, "dashed"),
        ("Stable", muted, "flat"),
    ]):
        tw = len(label) * 1.05 + 3.6
        lx -= tw
        if shape == "up":
            ax.add_patch(Polygon([(lx + 1.0, -7.2), (lx + 1.9, -5.6), (lx + 0.1, -5.6)], closed=True,
                                 facecolor=col, edgecolor="none", zorder=5))
        elif shape == "down":
            ax.add_patch(Polygon([(lx + 1.0, -5.6), (lx + 1.9, -7.2), (lx + 0.1, -7.2)], closed=True,
                                 facecolor=col, edgecolor="none", zorder=5))
        elif shape == "dashed":
            ax.add_line(Line2D([lx + 0.0, lx + 2.2], [-6.4, -6.4], color=col, lw=1.4,
                               linestyle=(0, (2.2, 1.4)), zorder=5))
        else:
            ax.add_line(Line2D([lx + 0.1, lx + 1.9], [-6.4, -6.4], color=col, lw=1.6, zorder=5))
        ax.text(lx + 2.5, -6.4, label, fontsize=6.2, color=ink2, ha="left", va="center")
        lx -= 1.2

    ax.plot([2.2, W - 2.2], [-1.0, -1.0], color=card_edge, lw=0.8 if not pixel_mode else 1.4, zorder=1)

    _draw_ga_bottom_panels(
        ax, spec,
        surface=surface, ink=ink, ink2=ink2, muted=muted,
        card_bg=card_bg, card_edge=card_edge, up=up,
        preset="nature_editorial", unit_in=unit_in, pixel_mode=pixel_mode, metro_mode=False,
    )


# ---------------------------------------------------------------------------
# 皮肤 1：nature_editorial（顶刊极简机制图）
# ---------------------------------------------------------------------------

def render_nature_editorial(ax, spec: AbstractSpec, *, unit_in: float = 0.062) -> None:
    pal = NATURE_EDITORIAL
    ax.set_facecolor(pal.surface)

    # 1. 区室容器
    for comp in spec.compartments:
        x0, y0, x1, y1 = comp.box
        if comp.kind == "blood":
            face, edge, ink = pal.boundary_face, pal.boundary_edge, pal.boundary_ink
        else:
            face, edge, ink = pal.mito_face, pal.mito_edge, pal.mito_ink
        ax.add_patch(FancyBboxPatch(
            (x0, y0), x1 - x0, y1 - y0,
            boxstyle="round,pad=0,rounding_size=2.8",
            facecolor=face, edgecolor=edge, lw=1.15, zorder=0.5,
        ))
        if comp.kind == "mitochondria":
            ax.add_patch(FancyBboxPatch(
                (x0 + 0.9, y0 + 0.9), (x1 - x0) - 1.8, (y1 - y0) - 1.8,
                boxstyle="round,pad=0,rounding_size=2.2",
                facecolor="none", edgecolor=pal.mito_inner, lw=0.65, linestyle=(0, (3.2, 2.0)), zorder=0.6,
            ))
            ax.text((x0 + x1) / 2, y1 - 2.2, comp.label, fontsize=7.0, fontweight="bold", color=ink,
                    ha="center", va="center", zorder=0.9)
        else:
            ax.text((x0 + x1) / 2, y0 + 2.2, comp.label, fontsize=7.0, fontweight="bold", color=ink,
                    ha="center", va="center", zorder=0.9)

    hub_by_id = {h.id: h for h in spec.hubs}

    # 2. 流向箭头
    for fl in spec.flows:
        sh, dh = hub_by_id[fl.source], hub_by_id[fl.target]
        p0, p1 = sh.xy, dh.xy
        dx, dy = p1[0] - p0[0], p1[1] - p0[1]
        length = math.hypot(dx, dy) or 1.0
        nx, ny = -dy / length, dx / length
        ctrl = ((p0[0] + p1[0]) / 2 + nx * fl.curve, (p0[1] + p1[1]) / 2 + ny * fl.curve)
        t0 = _box_exit_trim(ctrl[0] - p0[0], ctrl[1] - p0[1], HUB_W / 2, HUB_H / 2, pad=0.7)
        t1 = _box_exit_trim(p1[0] - ctrl[0], p1[1] - ctrl[1], HUB_W / 2, HUB_H / 2, pad=0.9)

        col = _flow_color(fl.direction, fl.evidence_tier, pal.up, pal.down, pal.edge_grey)
        strong = fl.evidence_tier in ("L3", "L2")
        width = 2.6 if strong else (1.7 if fl.direction else 1.15)
        ls = "-" if (strong or fl.direction == 0) else (0, (3.2, 2.0))

        if fl.coupled:
            from ..primitives import draw_coupled_track
            draw_coupled_track(ax, [{"p0": p0, "c": ctrl, "p1": p1, "trim0": t0, "trim1": t1}],
                               color=pal.coupled_track, edge_color=pal.coupled_edge, width=7.0, zorder=1.8)

        draw_quad_arrow(
            ax, p0, ctrl, p1, col, width, ls, False, True, t0, t1,
            zorder=3.2 if strong else 2.4, halo=strong, surface=pal.surface, arrow_scale=1.35,
        )
        if fl.label:
            mx = 0.25 * p0[0] + 0.5 * ctrl[0] + 0.25 * p1[0]
            my = 0.25 * p0[1] + 0.5 * ctrl[1] + 0.25 * p1[1]
            ax.text(
                mx, my - 0.6, fl.label, fontsize=6.0, style="italic", color=pal.ink2,
                ha="center", va="center", zorder=4.2,
                bbox=dict(boxstyle="round,pad=0.16", facecolor=pal.surface, edgecolor="none", alpha=0.88),
            )

    # 3. 枢纽卡片
    for hub in spec.hubs:
        x, y = hub.xy
        col, fill, strong, dashed = _hub_border_and_tint(hub, pal.up, pal.down, pal.muted, pal.surface)
        ec = col if strong else (pal.ink2 if dashed else pal.card_edge)
        lw = 1.65 if strong else (1.10 if dashed else 0.85)
        ls = (0, (3.0, 1.8)) if dashed else "-"

        ax.add_patch(FancyBboxPatch(
            (x - HUB_W / 2, y - HUB_H / 2), HUB_W, HUB_H,
            boxstyle="round,pad=0,rounding_size=1.8",
            facecolor=fill, edgecolor=ec, lw=lw, linestyle=ls, zorder=5.0,
        ))
        # 左侧方向色条
        if hub.direction:
            ax.add_patch(FancyBboxPatch(
                (x - HUB_W / 2 + 0.35, y - HUB_H / 2 + 0.85), 0.72, HUB_H - 1.70,
                boxstyle="round,pad=0,rounding_size=0.35",
                facecolor=col, edgecolor="none", zorder=5.3,
            ))

        draw_glyph(ax, hub.icon, (x - HUB_W / 2 + 3.35, y), size=4.2, zorder=5.6)

        title_fs = 6.2 if len(hub.label) >= 17 else 6.8
        ax.text(x - HUB_W / 2 + 6.0, y - 1.25, hub.label, fontsize=title_fs, fontweight="bold",
                color=pal.ink, ha="left", va="center", zorder=5.8)
        if hub.subtitle:
            ax.text(x - HUB_W / 2 + 6.0, y + 2.15, hub.subtitle, fontsize=5.9,
                    color=col if hub.direction else pal.muted, fontweight="medium" if hub.direction else "normal",
                    ha="left", va="center", zorder=5.8)
        if hub.evidence_tier in ("L1", "L2", "L3"):
            b_sc = 1.16
            bw = tier_badge_width(hub.evidence_tier, scale=b_sc, unit_in=unit_in)
            draw_tier_badge(ax, x + HUB_W / 2 - bw - 0.85, y - HUB_H / 2 + 1.45, hub.evidence_tier,
                            zorder=6.0, scale=b_sc, preset="nature_editorial", unit_in=unit_in)

    _draw_footer_and_header(
        ax, spec,
        surface=pal.surface, ink=pal.ink, ink2=pal.ink2, muted=pal.muted,
        card_bg=pal.card_bg, card_edge=pal.card_edge, up=pal.up, down=pal.down,
        unit_in=unit_in,
    )


# ---------------------------------------------------------------------------
# 皮肤 2：editorial_ink（手绘科学插画风：有机细胞器墨线淡彩 + 精密数据排版）
# ---------------------------------------------------------------------------

def _organic_rounded_rect_points(
    x0: float, y0: float, x1: float, y1: float, r: float, *, n_per_side: int = 24, jitter: float = 0.24, seed: int = 0
) -> np.ndarray:
    """生成带自然笔触抖动的闭合圆角矩形点列（确定性种子，逐字节可重现）。"""
    rng = np.random.default_rng(seed)
    pts: list[tuple[float, float]] = []
    # 四边 + 四圆角
    corners = [
        (x1 - r, y0 + r, -math.pi / 2, 0.0),
        (x1 - r, y1 - r, 0.0, math.pi / 2),
        (x0 + r, y1 - r, math.pi / 2, math.pi),
        (x0 + r, y0 + r, math.pi, 3 * math.pi / 2),
    ]
    for cx, cy, a0, a1 in corners:
        for t in np.linspace(a0, a1, 10):
            pts.append((cx + r * math.cos(t), cy + r * math.sin(t)))
    arr = np.array(pts, dtype=float)
    # 沿轮廓叠加低频平滑手绘波动（不是锯齿噪点）
    phase = rng.uniform(0, 2 * math.pi, size=3)
    s = np.linspace(0, 2 * math.pi, len(arr), endpoint=False)
    wave = (
        0.55 * np.sin(3 * s + phase[0])
        + 0.30 * np.sin(7 * s + phase[1])
        + 0.15 * np.sin(11 * s + phase[2])
    ) * jitter
    center = np.array([(x0 + x1) / 2.0, (y0 + y1) / 2.0])
    radial = arr - center
    norms = np.hypot(radial[:, 0], radial[:, 1])[:, None]
    norms = np.where(norms < 1e-6, 1.0, norms)
    return arr + (radial / norms) * wave[:, None]


def _draw_mito_cristae(ax, box: tuple[float, float, float, float], stroke: str, wash: str, zorder: float = 0.65) -> None:
    """在线粒体底部与顶部内壁绘制柔和的简笔内嵴（cristae）有机折褶。"""
    x0, y0, x1, y1 = box
    # 顶壁内嵴波浪
    xs_top = np.linspace(x0 + 6.0, x1 - 14.0, 90)
    ys_top = y0 + 1.85 + 1.15 * np.maximum(0.0, np.sin((xs_top - x0) * 0.42)) ** 2
    ax.plot(xs_top, ys_top, color=stroke, lw=0.85, alpha=0.55, solid_capstyle="round", zorder=zorder)
    ax.fill_between(xs_top, y0 + 1.2, ys_top, color=wash, alpha=0.35, zorder=zorder - 0.05)

    # 底壁内嵴波浪
    xs_bot = np.linspace(x0 + 8.0, x1 - 8.0, 90)
    ys_bot = y1 - 1.85 - 1.15 * np.maximum(0.0, np.cos((xs_bot - x0) * 0.38)) ** 2
    ax.plot(xs_bot, ys_bot, color=stroke, lw=0.85, alpha=0.55, solid_capstyle="round", zorder=zorder)
    ax.fill_between(xs_bot, ys_bot, y1 - 1.2, color=wash, alpha=0.35, zorder=zorder - 0.05)


def render_editorial_ink(ax, spec: AbstractSpec, *, unit_in: float = 0.062) -> None:
    # 暖宣纸 + 矿物颜料克制色系
    paper = "#FAF8F2"
    ink = "#1F242D"
    ink2 = "#444B57"
    muted = "#6E7582"
    up = "#BE3E34"       # 朱砂矿红
    down = "#245C94"     # 靛蓝矿青
    card_bg = "#F4F0E6"
    card_edge = "#CFC7B8"

    ax.set_facecolor(paper)
    W, H = spec.canvas

    # 极淡宣纸纹理点彩（确定性随机）
    rng = np.random.default_rng(42)
    fx = rng.uniform(1.0, W - 1.0, 140)
    fy = rng.uniform(0.5, H - 0.5, 140)
    ax.scatter(fx, fy, s=rng.uniform(0.4, 1.5, 140), color="#D8D0C0", alpha=0.35, linewidths=0, zorder=0.1)

    # 1. 有机手绘细胞器容器 + 淡彩晕染
    for idx, comp in enumerate(spec.compartments):
        x0, y0, x1, y1 = comp.box
        if comp.kind == "blood":
            wash, stroke, lbl_col = "#E9EFF4", "#6B8299", "#3B5268"
        else:
            wash, stroke, lbl_col = "#F3EBD9", "#9A8466", "#6B573C"
        outer = _organic_rounded_rect_points(x0, y0, x1, y1, 3.6, jitter=0.32, seed=10 + idx)
        outer2 = _organic_rounded_rect_points(x0 + 0.15, y0 + 0.15, x1 - 0.15, y1 - 0.15, 3.5, jitter=0.22, seed=20 + idx)
        # 多层淡彩水彩底晕
        ax.add_patch(Polygon(outer, closed=True, facecolor=wash, edgecolor="none", alpha=0.72, zorder=0.4))
        ax.add_patch(Polygon(outer, closed=True, facecolor="none", edgecolor=stroke, lw=1.45, alpha=0.88, zorder=0.55))
        ax.add_patch(Polygon(outer2, closed=True, facecolor="none", edgecolor=stroke, lw=0.55, alpha=0.45, zorder=0.58))
        if comp.kind == "mitochondria":
            _draw_mito_cristae(ax, comp.box, stroke=stroke, wash="#E4D5BC", zorder=0.62)
            ax.text((x0 + x1) / 2, y1 - 2.6, comp.label, fontsize=7.2, fontweight="bold", style="italic",
                    color=lbl_col, ha="center", va="center", zorder=0.9,
                    bbox=dict(boxstyle="round,pad=0.22", facecolor=wash, edgecolor="none", alpha=0.85))
        else:
            ax.text((x0 + x1) / 2, y0 + 2.4, comp.label, fontsize=7.2, fontweight="bold", style="italic",
                    color=lbl_col, ha="center", va="center", zorder=0.9)

    hub_by_id = {h.id: h for h in spec.hubs}

    # 2. 精准矢量流向箭头（保留纸墨手绘氛围但绝不抖动失真）
    for fl in spec.flows:
        sh, dh = hub_by_id[fl.source], hub_by_id[fl.target]
        p0, p1 = sh.xy, dh.xy
        dx, dy = p1[0] - p0[0], p1[1] - p0[1]
        length = math.hypot(dx, dy) or 1.0
        nx, ny = -dy / length, dx / length
        ctrl = ((p0[0] + p1[0]) / 2 + nx * fl.curve, (p0[1] + p1[1]) / 2 + ny * fl.curve)
        t0 = _box_exit_trim(ctrl[0] - p0[0], ctrl[1] - p0[1], HUB_W / 2, HUB_H / 2, pad=0.8)
        t1 = _box_exit_trim(p1[0] - ctrl[0], p1[1] - ctrl[1], HUB_W / 2, HUB_H / 2, pad=1.0)

        col = _flow_color(fl.direction, fl.evidence_tier, up, down, "#7A7469")
        strong = fl.evidence_tier in ("L3", "L2")
        width = 2.7 if strong else (1.75 if fl.direction else 1.2)
        ls = "-" if (strong or fl.direction == 0) else (0, (3.4, 2.0))

        if fl.coupled:
            from ..primitives import draw_coupled_track
            draw_coupled_track(ax, [{"p0": p0, "c": ctrl, "p1": p1, "trim0": t0, "trim1": t1}],
                               color="#DCE3F5", edge_color="#5B6AA6", width=7.2, zorder=1.8)

        draw_quad_arrow(
            ax, p0, ctrl, p1, col, width, ls, False, True, t0, t1,
            zorder=3.2 if strong else 2.4, halo=strong, surface=paper, arrow_scale=1.35,
        )
        if fl.label:
            mx = 0.25 * p0[0] + 0.5 * ctrl[0] + 0.25 * p1[0]
            my = 0.25 * p0[1] + 0.5 * ctrl[1] + 0.25 * p1[1]
            ax.text(
                mx, my - 0.6, fl.label, fontsize=6.1, style="italic", color=ink2,
                ha="center", va="center", zorder=4.2,
                bbox=dict(boxstyle="round,pad=0.16", facecolor=paper, edgecolor="none", alpha=0.90),
            )

    # 3. 针管笔双描边 + 淡彩晕染枢纽卡片（内文与徽章精密矢量排版）
    for idx, hub in enumerate(spec.hubs):
        x, y = hub.xy
        col, _, strong, dashed = _hub_border_and_tint(hub, up, down, muted, paper)
        wash = tint_color(col, 0.14 if strong else (0.09 if hub.direction else 0.03), base="#FFFDF8")
        x0, y0, x1, y1 = x - HUB_W / 2, y - HUB_H / 2, x + HUB_W / 2, y + HUB_H / 2
        poly1 = _organic_rounded_rect_points(x0, y0, x1, y1, 2.0, jitter=0.16, seed=100 + idx)
        poly2 = _organic_rounded_rect_points(x0 + 0.18, y0 + 0.22, x1 + 0.18, y1 + 0.22, 2.0, jitter=0.12, seed=200 + idx)
        # 柔和纸面投影
        ax.add_patch(Polygon(poly2, closed=True, facecolor="#E5DFD2", edgecolor="none", alpha=0.55, zorder=4.8))
        # 水彩卡片底
        ax.add_patch(Polygon(poly1, closed=True, facecolor=wash, edgecolor="none", zorder=5.0))
        # 针管笔双轮廓描边
        ec = col if strong else (ink2 if dashed else "#8C8577")
        lw = 1.65 if strong else (1.05 if dashed else 0.90)
        ls = (0, (3.2, 1.8)) if dashed else "-"
        ax.add_patch(Polygon(poly1, closed=True, facecolor="none", edgecolor=ec, lw=lw, linestyle=ls, zorder=5.2))
        if strong:
            ax.add_patch(Polygon(
                _organic_rounded_rect_points(x0 + 0.32, y0 + 0.32, x1 - 0.32, y1 - 0.32, 1.8, jitter=0.10, seed=300 + idx),
                closed=True, facecolor="none", edgecolor=col, lw=0.45, alpha=0.45, zorder=5.25,
            ))

        draw_glyph(ax, hub.icon, (x - HUB_W / 2 + 3.35, y), style="editorial_ink", size=4.2, zorder=5.6)

        title_fs = 6.2 if len(hub.label) >= 17 else 6.8
        ax.text(x - HUB_W / 2 + 6.0, y - 1.25, hub.label, fontsize=title_fs, fontweight="bold",
                color=ink, ha="left", va="center", zorder=5.8)
        if hub.subtitle:
            ax.text(x - HUB_W / 2 + 6.0, y + 2.15, hub.subtitle, fontsize=5.9,
                    color=col if hub.direction else muted, fontweight="medium" if hub.direction else "normal",
                    ha="left", va="center", zorder=5.8)
        if hub.evidence_tier in ("L1", "L2", "L3"):
            b_sc = 1.16
            bw = tier_badge_width(hub.evidence_tier, scale=b_sc, unit_in=unit_in)
            draw_tier_badge(ax, x + HUB_W / 2 - bw - 0.85, y - HUB_H / 2 + 1.45, hub.evidence_tier,
                            zorder=6.0, scale=b_sc, preset="nature_editorial", unit_in=unit_in)

    _draw_footer_and_header(
        ax, spec,
        surface=paper, ink=ink, ink2=ink2, muted=muted,
        card_bg=card_bg, card_edge=card_edge, up=up, down=down,
        unit_in=unit_in,
    )


# ---------------------------------------------------------------------------
# 皮肤 3：pixel_abstract（现代 16-bit 像素概念风 — 严格 ≤ 12 枢纽）
# ---------------------------------------------------------------------------

# 12×12 程序化像素精灵模板（字符映射到语义调色板，零抗锯齿硬边缘）
def _draw_pixel_frame(ax, x0: float, y0: float, w: float, h: float, *, face: str, edge: str, shadow: str, lw: float = 1.4, zorder: float = 5.0) -> None:
    """绘制硬边缘切角 16-bit 像素窗口框（含 1 格像素硬阴影）。"""
    c = 0.75
    # 像素硬阴影
    ax.add_patch(Rectangle((x0 + 0.65, y0 + 0.65), w, h, facecolor=shadow, edgecolor="none", zorder=zorder - 0.15))
    # 八边形像素切角框
    verts = [
        (x0 + c, y0), (x0 + w - c, y0),
        (x0 + w, y0 + c), (x0 + w, y0 + h - c),
        (x0 + w - c, y0 + h), (x0 + c, y0 + h),
        (x0, y0 + h - c), (x0, y0 + c),
    ]
    ax.add_patch(Polygon(verts, closed=True, facecolor=face, edgecolor=edge, lw=lw, joinstyle="miter", zorder=zorder))
    # 顶部 1px 高光边
    ax.plot([x0 + c + 0.3, x0 + w - c - 0.3], [y0 + 0.45, y0 + 0.45], color="#FFFFFF", lw=0.9, alpha=0.85, zorder=zorder + 0.05)


def _pixel_tier_meter_width(tier: str, *, unit_in: float = 0.062) -> float:
    text = tier_badge_text(tier)
    if not text:
        return 0.0
    sq, gap = 0.86, 0.24
    return 0.30 + 3 * sq + 2 * gap + 0.35 + estimate_text_width(text, 4.2, unit_in, bold=True) * 0.96 + 0.38


def _draw_pixel_tier_meter(
    ax,
    x: float,
    y: float,
    tier: str,
    *,
    ink: str,
    up_col: str,
    unit_in: float = 0.062,
    zorder: float = 6.0,
) -> float:
    """3 格像素方块能量条（■■■ Confirmed / ■■□ Flux / ■□□ Enz-pred）。"""
    n_on = {"L3": 3, "L2": 2, "L1": 1}.get(tier, 0)
    if not n_on:
        return 0.0
    sq = 0.86
    gap = 0.24
    bg = "#1E2430" if n_on >= 2 else "#FFFDF7"
    fg = "#FFFDF7" if n_on >= 2 else "#1E2430"
    bw = _pixel_tier_meter_width(tier, unit_in=unit_in)
    bh = 1.50
    ax.add_patch(Rectangle((x, y - bh / 2), bw, bh, facecolor=bg, edgecolor=ink, lw=0.85, zorder=zorder))
    for i in range(3):
        sx = x + 0.30 + i * (sq + gap)
        filled = i < n_on
        ax.add_patch(Rectangle(
            (sx, y - sq / 2), sq, sq,
            facecolor=(up_col if filled else "none"),
            edgecolor=(fg if not filled else "none"),
            lw=0.55, zorder=zorder + 0.1,
        ))
    ax.text(x + 0.30 + 3 * sq + 2 * gap + 0.35, y + 0.02, tier_badge_text(tier), fontsize=4.2, fontweight="bold",
            family="monospace", color=fg, ha="left", va="center", zorder=zorder + 0.12)
    return bw


def _draw_manhattan_pixel_arrow(
    ax, p0: tuple[float, float], p1: tuple[float, float], color: str, width: float, dashed: bool, zorder: float
) -> None:
    """正交像素管线 + 阶梯像素箭头。"""
    dx, dy = p1[0] - p0[0], p1[1] - p0[1]
    t0 = _box_exit_trim(dx, dy, HUB_W / 2, HUB_H / 2, pad=0.6)
    t1 = _box_exit_trim(-dx, -dy, HUB_W / 2, HUB_H / 2, pad=1.3)
    norm = math.hypot(dx, dy) or 1.0
    ux, uy = dx / norm, dy / norm
    s = (p0[0] + ux * t0, p0[1] + uy * t0)
    e = (p1[0] - ux * t1, p1[1] - uy * t1)

    ls = (0, (2.2, 1.6)) if dashed else "-"
    if abs(e[0] - s[0]) > 2.5 and abs(e[1] - s[1]) > 2.5:
        mid_x = round((s[0] + e[0]) / 2.0)
        xs = [s[0], mid_x, mid_x, e[0]]
        ys = [s[1], s[1], e[1], e[1]]
        dir_vec = (1.0 if e[0] >= mid_x else -1.0, 0.0)
    else:
        xs = [s[0], e[0]]
        ys = [s[1], e[1]]
        dir_vec = (ux, uy)

    ax.plot(xs, ys, color="#1E2430", lw=width + 1.1, solid_joinstyle="miter", solid_capstyle="butt", alpha=0.22, zorder=zorder - 0.1)
    ax.plot(xs, ys, color=color, lw=width, linestyle=ls, solid_joinstyle="miter", solid_capstyle="butt", zorder=zorder)
    draw_stealth_head(ax, e, dir_vec, color, size=1.85, zorder=zorder + 0.15)


def render_pixel_abstract(ax, spec: AbstractSpec, *, unit_in: float = 0.062) -> None:
    # 克制矿物系 16-bit 像素调色板（拒绝廉价高饱和霓虹，保持科研演示质感）
    bg = "#F4EFE4"
    ink = "#1E2430"
    ink2 = "#3D4656"
    muted = "#6D7688"
    up = "#D4453B"
    down = "#2B6CB0"
    gold = "#D99B26"
    shadow = "#D4CCB8"
    card_bg = "#FFFDF7"

    ax.set_facecolor(bg)
    W, H = spec.canvas

    # 极淡像素点阵底格（增强像素世界地图秩序感，不干扰阅读）
    gx = np.arange(2, int(W), 4)
    gy = np.arange(2, int(H), 4)
    GXX, GYY = np.meshgrid(gx, gy)
    ax.scatter(GXX.ravel(), GYY.ravel(), s=1.0, color="#DDD5C2", alpha=0.65, marker="s", linewidths=0, zorder=0.1)

    # 1. 像素区室瓦片区
    for comp in spec.compartments:
        x0, y0, x1, y1 = comp.box
        if comp.kind == "blood":
            c_face, c_edge, c_ink = "#E4EDF6", "#5E81AC", "#2E4A6B"
        else:
            c_face, c_edge, c_ink = "#EFE6D2", "#9E845C", "#5E4B2E"
        _draw_pixel_frame(ax, x0, y0, x1 - x0, y1 - y0, face=c_face, edge=c_edge, shadow="#DDD4C0", lw=1.6, zorder=0.5)
        ax.add_patch(Rectangle((x0 + 1.0, y0 + 1.0), (x1 - x0) - 2.0, (y1 - y0) - 2.0,
                               facecolor="none", edgecolor=c_edge, lw=0.75, linestyle=(0, (2.0, 2.0)), zorder=0.6))
        # 居中像素铭牌（线粒体放底部居中，彻底避开顶部 Acyl-CoA 箭头）
        lbl_w = len(comp.label) * 1.18 + 2.8
        lbl_x = (x0 + x1) / 2 - lbl_w / 2
        lbl_y = (y1 - 4.0) if comp.kind == "mitochondria" else (y0 + 1.2)
        ax.add_patch(Rectangle((lbl_x, lbl_y), lbl_w, 2.7, facecolor=c_edge, edgecolor=ink, lw=0.8, zorder=0.85))
        ax.text((x0 + x1) / 2, lbl_y + 1.35, comp.label.upper(), fontsize=5.8, fontweight="bold",
                family="monospace", color="#FFFDF7", ha="center", va="center", zorder=0.9)

    hub_by_id = {h.id: h for h in spec.hubs}

    # 2. 正交像素流向管线
    for fl in spec.flows:
        sh, dh = hub_by_id[fl.source], hub_by_id[fl.target]
        col = _flow_color(fl.direction, fl.evidence_tier, up, down, "#7A8291")
        strong = fl.evidence_tier in ("L3", "L2")
        width = 2.6 if strong else (1.8 if fl.direction else 1.25)
        dashed = bool(fl.direction and not strong)
        _draw_manhattan_pixel_arrow(ax, sh.xy, dh.xy, col, width, dashed, zorder=2.6 if strong else 2.0)
        if fl.label:
            mx, my = (sh.xy[0] + dh.xy[0]) / 2, (sh.xy[1] + dh.xy[1]) / 2
            ax.text(
                mx, my - 0.5, fl.label, fontsize=5.6, family="monospace", color=ink2,
                ha="center", va="center", zorder=4.2,
                bbox=dict(boxstyle="square,pad=0.14", facecolor=bg, edgecolor="none", alpha=0.92),
            )

    # 3. 16-bit 像素枢纽卡片 + 16×16 硬边矢量精灵
    for hub in spec.hubs:
        x, y = hub.xy
        col, _, strong, dashed = _hub_border_and_tint(hub, up, down, muted, card_bg)
        face = tint_color(col, 0.10 if strong else (0.06 if hub.direction else 0.0), base=card_bg)
        edge = col if strong else (ink if dashed else "#6E7688")
        _draw_pixel_frame(
            ax, x - HUB_W / 2, y - HUB_H / 2, HUB_W, HUB_H,
            face=face, edge=edge, shadow=shadow, lw=1.8 if strong else 1.25, zorder=5.0,
        )
        draw_glyph(ax, hub.icon, (x - HUB_W / 2 + 3.6, y), style="pixel_abstract", size=4.8, zorder=5.8)
        title_fs = 6.0 if len(hub.label) >= 17 else 6.5
        ax.text(x - HUB_W / 2 + 6.6, y - 1.25, hub.label, fontsize=title_fs, fontweight="bold",
                color=ink, ha="left", va="center", zorder=5.9)
        if hub.subtitle:
            ax.text(x - HUB_W / 2 + 6.6, y + 2.15, hub.subtitle, fontsize=5.5, family="monospace",
                    color=col if hub.direction else muted, fontweight="bold" if strong else "normal",
                    ha="left", va="center", zorder=5.9)
        if hub.evidence_tier in ("L1", "L2", "L3"):
            bw = _pixel_tier_meter_width(hub.evidence_tier, unit_in=unit_in)
            _draw_pixel_tier_meter(ax, x + HUB_W / 2 - bw - 0.75, y - HUB_H / 2 + 1.55, hub.evidence_tier,
                                   ink=ink, up_col=gold if strong else ink, unit_in=unit_in, zorder=6.1)

    _draw_footer_and_header(
        ax, spec,
        surface=bg, ink=ink, ink2=ink2, muted=muted,
        card_bg=card_bg, card_edge=ink, up=up, down=down, pixel_mode=True,
        unit_in=unit_in,
    )


# ---------------------------------------------------------------------------
# 皮肤 4：metro_atlas（Beck 代谢地铁换乘图 — 八向轨道 + 内嵌 chevron + 跑道闭环 + 圆环换乘站）
# ---------------------------------------------------------------------------

def _is_metro_loop_hub(hub: AbstractHub) -> bool:
    """循环类枢纽（TCA / PPP / Cycle）渲染为标志性跑道闭环。"""
    low_id = hub.id.lower()
    low_lbl = hub.label.lower()
    return hub.icon in ("tca", "ppp") or "cycle" in low_id or "cycle" in low_lbl or "loop" in low_id


def _metro_loop_dims(hub: AbstractHub) -> tuple[float, float]:
    if hub.icon == "tca" or hub.evidence_tier in ("L3", "L2"):
        return 24.4, 13.8
    return 21.6, 10.6


def _metro_hub_color(hub: AbstractHub, pal, color_mode: str) -> str:
    if color_mode == "pathway":
        return pathway_family_color(hub.id, fallback_accent=pathway_family_color(hub.icon))
    if hub.direction > 0:
        return pal.up
    if hub.direction < 0:
        return pal.down
    return pal.edge_grey


def _metro_flow_color(fl, sh: AbstractHub, dh: AbstractHub, pal, color_mode: str) -> str:
    if color_mode == "pathway":
        # 汇入 TCA 的支线优先继承支线来源色（脂肪酸苔绿 / 氨基酸湖蓝 / 糖酵解暖棕），复刻经典代谢地铁分层
        if dh.id == "tca_cycle" and sh.id != "pyruvate_hub":
            return _metro_hub_color(sh, pal, "pathway")
        return _metro_hub_color(dh, pal, "pathway")
    return _flow_color(fl.direction, fl.evidence_tier, pal.up, pal.down, pal.omics_grey)


def _loop_perimeter_anchor(loop_xy: tuple[float, float], dims: tuple[float, float], other_xy: tuple[float, float]) -> tuple[float, float]:
    """计算外部轨道接入圆角跑道闭环的平滑锚点（直壁正交接入或顶/底直段接入）。"""
    cx, cy = loop_xy
    lw, lh = dims
    ox, oy = other_xy
    dx, dy = ox - cx, oy - cy
    if abs(dy) <= lh * 0.42:
        sx = 1.0 if dx >= 0 else -1.0
        ay = max(cy - lh * 0.28, min(cy + lh * 0.28, oy))
        return (cx + sx * lw * 0.5, ay)
    if abs(dx) <= lw * 0.22:
        sy = 1.0 if dy >= 0 else -1.0
        return (cx + dx * 0.4, cy + sy * lh * 0.5)
    sx = 1.0 if dx >= 0 else -1.0
    sy = 1.0 if dy >= 0 else -1.0
    return (cx + sx * lw * 0.23, cy + sy * lh * 0.5)


def _draw_metro_header_and_footer(
    ax,
    spec: AbstractSpec,
    *,
    pal,
    color_mode: str,
    unit_in: float,
) -> None:
    W, H = spec.canvas
    ax.text(2.2, -7.8, spec.title, fontsize=11.0, fontweight="bold", color=pal.ink, ha="left", va="center")
    sub = spec.subtitle
    if color_mode == "pathway":
        mode_tag = "Metro Atlas · Pathway-Family Colors (Solid = Flux · Hollow = Enz-pred)"
    else:
        mode_tag = "Metro Atlas · Semantic Flux Colors (Solid = Flux · Hollow = Enz-pred)"
    ax.text(2.2, -3.5, f"{sub}  ·  {mode_tag}" if sub else mode_tag, fontsize=6.7, color=pal.muted, ha="left", va="center")

    # 右上角真实地铁轨道样本图例（紧凑排列，与左侧主标题保持充足留白）
    lx = W - 2.2
    if color_mode == "pathway":
        items = [
            ("Flux", pal.ink2, False, 3.4),
            ("Enz-pred", pal.ink2, True, 3.2),
            ("Stable", pal.omics_grey, False, 2.2),
        ]
    else:
        items = [
            ("Up (Flux)", pal.up, False, 3.4),
            ("Down", pal.down, False, 3.4),
            ("Enz-pred", pal.ink2, True, 3.1),
            ("Stable", pal.omics_grey, False, 2.2),
        ]
    seg_w = 4.4
    for label, col, hollow, lw in reversed(items):
        tw = len(label) * 0.96 + seg_w + 1.2
        lx -= tw
        draw_metro_track(
            ax, [(lx, -7.6), (lx + seg_w, -7.6)], col,
            width=lw, hollow=hollow, surface=pal.surface,
            chevrons=True, chevron_spacing=10.0, unit_in=unit_in,
            corner_radius=0.0, halo=False, zorder=5.2,
        )
        ax.text(lx + seg_w + 0.75, -7.6, label, fontsize=5.7, color=pal.ink2, ha="left", va="center")
        lx -= 1.2

    ax.plot([2.2, W - 2.2], [-1.0, -1.0], color=pal.card_edge, lw=0.85, zorder=1)

    _draw_ga_bottom_panels(
        ax, spec,
        surface=pal.surface, ink=pal.ink, ink2=pal.ink2, muted=pal.muted,
        card_bg=pal.card_bg, card_edge=pal.card_edge, up=pal.up,
        preset="metro_atlas", unit_in=unit_in, pixel_mode=False, metro_mode=True,
    )


def render_metro_atlas(
    ax,
    spec: AbstractSpec,
    *,
    color_mode: str = "semantic",
    unit_in: float = 0.062,
) -> None:
    """Beck 代谢地铁图渲染器（八向轨道 + 内嵌 chevron + 跑道闭环 + 圆环换乘站）。"""
    if color_mode not in VALID_COLOR_MODES:
        raise ValueError(f"metro_atlas color_mode 应为 {' / '.join(VALID_COLOR_MODES)} 之一，实际为 {color_mode!r}")

    pal = METRO_ATLAS
    ax.set_facecolor(pal.surface)
    W, H = spec.canvas

    # 1. 低存在感浅灰圆角区室底板（100% 视觉焦点让位给地铁轨道与换乘站）
    comp_label_boxes: list[tuple[float, float, float, float]] = []
    for comp in spec.compartments:
        x0, y0, x1, y1 = comp.box
        if comp.kind == "blood":
            draw_metro_zone(
                ax, comp.box, comp.label,
                face="#EEF2F6", edge="#CBD5E1", lw=0.65,
                label_color="#526B84", label_pos="top_center", rounding=3.2, fontsize=6.6, zorder=0.4,
            )
            comp_label_boxes.append(((x0 + x1) / 2 - 10.0, y0, (x0 + x1) / 2 + 10.0, y0 + 4.2))
        elif comp.kind == "mitochondria":
            draw_metro_zone(
                ax, comp.box, comp.label,
                face="#EFEFF3", edge="#D5D9E2", lw=0.75, inner_edge="#CBD2DD",
                label_color="#64748B", label_pos="bottom_center", rounding=3.6, fontsize=6.8, zorder=0.4,
            )
            comp_label_boxes.append(((x0 + x1) / 2 - 15.0, y1 - 4.4, (x0 + x1) / 2 + 15.0, y1))
        else:
            draw_metro_zone(
                ax, comp.box, comp.label,
                face="#EDF5F3", edge="#B8D8D1", lw=0.65,
                label_color="#3D736B", label_pos="top_right", rounding=3.0, fontsize=6.4, zorder=0.4,
            )

    hub_by_id = {h.id: h for h in spec.hubs}
    loop_hubs = {h.id: _metro_loop_dims(h) for h in spec.hubs if _is_metro_loop_hub(h)}
    station_r_map: dict[str, float] = {
        h.id: (2.45 if h.evidence_tier in ("L3", "L2") else 2.15)
        for h in spec.hubs if h.id not in loop_hubs
    }

    # 2. 路由并绘制所有八向地铁轨道（记录每个站台的入射角度，供标签零冲突自动避让）
    incident_angles: dict[str, list[float]] = {h.id: [] for h in spec.hubs}
    loop_junctions: list[tuple[tuple[float, float], str, str]] = []

    for fl in spec.flows:
        sh, dh = hub_by_id[fl.source], hub_by_id[fl.target]
        src_is_loop = sh.id in loop_hubs
        dst_is_loop = dh.id in loop_hubs

        p0 = _loop_perimeter_anchor(sh.xy, loop_hubs[sh.id], dh.xy) if src_is_loop else sh.xy
        p1 = _loop_perimeter_anchor(dh.xy, loop_hubs[dh.id], sh.xy) if dst_is_loop else dh.xy

        prefer = "hv_first" if fl.curve >= 0 else "diag_first"
        # 斜向入环或分叉优先居中或干净单折
        if abs(p1[0] - p0[0]) > 4.0 and abs(p1[1] - p0[1]) > 4.0 and abs(fl.curve) < 1e-4:
            prefer = "centered"
        raw_poly = route_octilinear(p0, p1, bend=0.0, prefer=prefer)
        smooth_poly = fillet_polyline(raw_poly, radius=2.6, arc_steps=8)

        # 记录起点/终点在站台处的真实出射方向角
        if len(smooth_poly) >= 2:
            if not src_is_loop:
                a0 = math.atan2(smooth_poly[1][1] - smooth_poly[0][1], smooth_poly[1][0] - smooth_poly[0][0])
                incident_angles[sh.id].append(a0)
            if not dst_is_loop:
                a1 = math.atan2(smooth_poly[-2][1] - smooth_poly[-1][1], smooth_poly[-2][0] - smooth_poly[-1][0])
                incident_angles[dh.id].append(a1)

        col = _metro_flow_color(fl, sh, dh, pal, color_mode)
        strong = fl.evidence_tier in ("L3", "L2")
        is_l1 = fl.evidence_tier == "L1"
        hollow = is_l1
        width = 4.2 if strong else (3.5 if is_l1 else 2.5)
        trim0 = 0.0 if src_is_loop else (station_r_map[sh.id] + 0.10)
        trim1 = 0.0 if dst_is_loop else (station_r_map[dh.id] + 0.10)

        if fl.coupled:
            draw_metro_coupled_halo(
                ax, [smooth_poly],
                color=pal.coupled_track, edge_color=pal.coupled_edge, width=9.0, zorder=1.8,
            )

        drawn = draw_metro_track(
            ax, raw_poly, col,
            width=width, hollow=hollow, surface=pal.surface,
            chevrons=True, chevron_spacing=7.5, chevron_direction=1,
            unit_in=unit_in, trim0=trim0, trim1=trim1,
            corner_radius=2.6, halo=True, zorder=2.8 if strong else 2.3,
        )

        # 记录入环/出环交点，稍后绘制在跑道环边缘的标志性换乘圆点
        if src_is_loop:
            loop_junctions.append((p0, col, _metro_hub_color(sh, pal, color_mode)))
        if dst_is_loop:
            loop_junctions.append((p1, col, _metro_hub_color(dh, pal, color_mode)))

        if fl.label and len(drawn) >= 2:
            total = polyline_length(drawn)
            mid_pt, tan = sample_polyline_at(drawn, total * 0.50)
            # 沿轨道法向偏移，避开轨道内嵌 chevron（斜向轨道多留一点垂直净空）
            nx, ny = -tan[1], tan[0]
            if ny > 0.35:
                nx, ny = -nx, -ny
            is_diag = abs(tan[0]) > 0.25 and abs(tan[1]) > 0.25
            off = 2.85 if is_diag else 2.05
            lx, ly = mid_pt[0] + nx * off, mid_pt[1] + ny * off
            ax.text(
                lx, ly, fl.label, fontsize=5.9, style="italic", color=pal.ink2,
                ha="center", va="center", zorder=4.4,
                bbox=dict(boxstyle="round,pad=0.14", facecolor=pal.surface, edgecolor="none", alpha=0.90),
            )

    # 3. 绘制跑道闭环枢纽（TCA Cycle / Pentose Phosphate 等）
    for hub in spec.hubs:
        if hub.id not in loop_hubs:
            continue
        lw, lh = loop_hubs[hub.id]
        col = _metro_hub_color(hub, pal, color_mode)
        strong = hub.evidence_tier in ("L3", "L2")
        hollow = hub.evidence_tier == "L1"
        sub_col = (pal.up if hub.direction > 0 else (pal.down if hub.direction < 0 else pal.muted))
        draw_metro_loop(
            ax, hub.xy, lw, lh, col,
            track_width=4.6 if strong else 3.6,
            hollow=hollow,
            surface=pal.surface,
            interior_tint=0.065 if strong else 0.035,
            clockwise=True,
            n_chevrons=6 if lw >= 23.0 else 4,
            label=hub.label,
            subtitle=hub.subtitle,
            label_color=pal.ink if color_mode == "semantic" else col,
            subtitle_color=sub_col,
            label_fontsize=7.6 if strong else 6.8,
            subtitle_fontsize=6.0,
            unit_in=unit_in,
            zorder=3.6 if strong else 3.2,
        )
        if hub.evidence_tier in ("L1", "L2", "L3"):
            b_sc = 1.34
            bw = tier_badge_width(hub.evidence_tier, scale=b_sc, unit_in=unit_in)
            draw_tier_badge(
                ax, hub.xy[0] - bw * 0.5, hub.xy[1] - lh * 0.5 - 2.35,
                hub.evidence_tier, zorder=6.0, scale=b_sc, preset="metro_atlas", unit_in=unit_in,
            )

    # 3b. 跑道闭环边缘换乘圆点（复刻参考图 TCA 环壁圆站）
    for jpt, track_col, loop_col in loop_junctions:
        draw_metro_station(
            ax, jpt, 1.25,
            fill=tint_color(loop_col, 0.22, base=pal.surface),
            edge=loop_col, lw=1.6, surface=pal.surface, halo_pad=0.35, zorder=4.8,
        )

    # 4. 绘制双层圆环换乘站台 + 自动避让轨道的站台铭牌
    for hub in spec.hubs:
        if hub.id in loop_hubs:
            continue
        x, y = hub.xy
        r = station_r_map[hub.id]
        col = _metro_hub_color(hub, pal, color_mode)
        strong = hub.evidence_tier in ("L3", "L2")
        fill = tint_color(col, 0.24 if strong else (0.16 if hub.direction else 0.08), base=pal.surface)
        draw_metro_station(
            ax, (x, y), r,
            fill=fill,
            edge=col,
            outer_ring=pal.ink if strong else None,
            lw=2.2 if strong else 1.85,
            surface=pal.surface,
            halo_pad=0.55,
            zorder=5.2,
        )

        # 选择离入射轨道最远、且不撞区室标题/画布边界的方位放置站台文字块
        inc = incident_angles.get(hub.id, [])
        candidates = [
            # (angle, ox, oy, ha, va, pref_penalty)
            (math.pi / 2, 0.0, r + 1.85, "center", "top", 0.00),      # 下方
            (-math.pi / 2, 0.0, -(r + 2.05), "center", "bottom", 0.03), # 上方
            (0.0, r + 2.3, 0.0, "left", "center", 0.16),               # 右方
            (math.pi, -(r + 2.3), 0.0, "right", "center", 0.18),       # 左方
        ]
        lbl_w = max(len(hub.label), len(hub.subtitle)) * 1.02 + 2.4
        lbl_h = 5.4 if hub.subtitle else 3.0
        best = candidates[0]
        best_score = -1e9
        for ang, ox, oy, ha, va, pref in candidates:
            clearance = min(abs(math.atan2(math.sin(ang - a), math.cos(ang - a))) for a in inc) if inc else math.pi
            cx_l, cy_l = x + ox, y + oy
            bx0 = cx_l - (lbl_w / 2 if ha == "center" else (lbl_w if ha == "right" else 0.0))
            by0 = cy_l - (lbl_h / 2 if va == "center" else (lbl_h if va == "bottom" else 0.0))
            bx1, by1 = bx0 + lbl_w, by0 + lbl_h
            pen = 0.0
            if bx0 < 1.5 or bx1 > W - 1.5 or by0 < 0.8 or by1 > H - 1.0:
                pen += 4.0
            for c0, cy0, c1, cy1 in comp_label_boxes:
                if min(bx1, c1) > max(bx0, c0) and min(by1, cy1) > max(by0, cy0):
                    pen += 3.0
            score = min(clearance, 1.45) - pref - pen
            if score > best_score:
                best_score = score
                best = (ang, ox, oy, ha, va, pref)

        _, ox, oy, ha, va, _ = best
        tx, ty = x + ox, y + oy
        sub_col = pal.up if hub.direction > 0 else (pal.down if hub.direction < 0 else pal.muted)
        if va == "top":
            y_title, y_sub, y_badge = ty + 0.45, ty + 2.55, ty + 4.45
        elif va == "bottom":
            y_title = ty - (2.45 if hub.subtitle else 0.5)
            y_sub = ty - 0.45
            y_badge = y_title - 2.15
        else:
            y_title, y_sub, y_badge = ty - 1.15, ty + 1.25, ty - 3.2

        ax.text(
            tx, y_title, hub.label, fontsize=6.9, fontweight="bold", color=pal.ink,
            ha=ha, va="center", zorder=5.8,
            bbox=dict(boxstyle="round,pad=0.14", facecolor=pal.surface, edgecolor="none", alpha=0.88),
        )
        if hub.subtitle:
            ax.text(
                tx, y_sub, hub.subtitle, fontsize=5.8, color=sub_col,
                fontweight="medium" if hub.direction else "normal",
                ha=ha, va="center", zorder=5.8,
                bbox=dict(boxstyle="round,pad=0.10", facecolor=pal.surface, edgecolor="none", alpha=0.84),
            )
        if hub.evidence_tier in ("L1", "L2", "L3"):
            b_sc = 1.26
            bw = tier_badge_width(hub.evidence_tier, scale=b_sc, unit_in=unit_in)
            bx = tx - bw / 2 if ha == "center" else (tx - bw if ha == "right" else tx)
            draw_tier_badge(
                ax, bx, y_badge, hub.evidence_tier, zorder=6.0, scale=b_sc, preset="metro_atlas", unit_in=unit_in,
            )

    _draw_metro_header_and_footer(ax, spec, pal=pal, color_mode=color_mode, unit_in=unit_in)


# ---------------------------------------------------------------------------
# 皮肤 5：economist（《经济学人》编辑部几何图表风 — 默认主风格）
# ---------------------------------------------------------------------------

def render_economist(ax, spec: AbstractSpec, *, color_mode: str = "semantic", unit_in: float = 0.062) -> None:
    """《经济学人》编辑部视觉风格的 Graphical Abstract（暖灰黄纸底、红标页眉、锐利方卡与红/蓝/灰克制语义）。"""
    pal = ECONOMIST
    W, _H = spec.canvas
    ax.set_facecolor(pal.surface)

    # 1. 标志性左上角红标与编辑部页眉
    draw_economist_masthead(
        ax, 2.2, -10.5, W - 4.4, tag_w=8.8, tag_h=1.45,
        color=ECONOMIST_RED, rule_color="#B8B8AF", y_up=False, zorder=8.0,
    )
    ax.text(2.2, -6.9, spec.title, fontsize=11.4, fontweight="bold", color=pal.ink, ha="left", va="center")
    if spec.subtitle:
        ax.text(2.2, -3.1, spec.subtitle, fontsize=7.3, color=pal.ink2, ha="left", va="center")

    # 右上角编辑部极简方向图例
    lx = W - 2.4
    for label, col, shape in reversed([
        ("Up (Flux)", pal.up, "up"),
        ("Down", pal.down, "down"),
        ("Enz-pred", pal.up, "hollow"),
        ("Stable", pal.edge_grey, "flat"),
    ]):
        tw = len(label) * 1.05 + 3.8
        lx -= tw
        if shape == "up":
            ax.add_patch(Polygon([(lx + 1.0, -7.2), (lx + 1.9, -5.6), (lx + 0.1, -5.6)], closed=True,
                                 facecolor=col, edgecolor="none", zorder=5))
        elif shape == "down":
            ax.add_patch(Polygon([(lx + 1.0, -5.6), (lx + 1.9, -7.2), (lx + 0.1, -7.2)], closed=True,
                                 facecolor=col, edgecolor="none", zorder=5))
        elif shape == "hollow":
            soft_up = tint_color(pal.up, 0.65, base=pal.card_bg)
            draw_metro_track(
                ax, [(lx - 0.15, -6.4), (lx + 2.15, -6.4)], soft_up,
                width=2.1, hollow=True, surface=pal.card_bg, chevrons=False, halo=False, zorder=5.0,
            )
        else:
            ax.add_line(Line2D([lx + 0.1, lx + 2.1], [-6.4, -6.4], color=col, lw=1.8, zorder=5))
        ax.text(lx + 2.85, -6.4, label, fontsize=6.2, color=pal.ink2, ha="left", va="center")
        lx -= 1.2
    ax.plot([2.2, W - 2.2], [-1.0, -1.0], color=pal.ink, lw=0.95, solid_capstyle="butt", zorder=2)

    # 2. 区室容器（方角锐利编辑部面板 + 顶部细黑线）
    for comp in spec.compartments:
        x0, y0, x1, y1 = comp.box
        if comp.kind == "blood":
            face, edge, ink = pal.boundary_face, pal.boundary_edge, pal.boundary_ink
        else:
            face, edge, ink = pal.mito_face, pal.mito_edge, pal.mito_ink
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor=face, edgecolor=pal.card_edge, lw=0.65, zorder=0.5))
        ax.plot([x0, x1], [y0, y0], color=edge, lw=1.25, solid_capstyle="butt", zorder=0.7)
        if comp.kind == "mitochondria":
            ax.add_patch(Rectangle(
                (x0 + 0.8, y0 + 0.8), (x1 - x0) - 1.6, (y1 - y0) - 1.6,
                facecolor="none", edgecolor=pal.mito_inner, lw=0.55, linestyle=(0, (3.0, 2.0)), zorder=0.6,
            ))
            ax.text(x0 + 1.6, y1 - 2.0, comp.label.upper(), fontsize=6.4, fontweight="bold", color=ink,
                    ha="left", va="center", zorder=0.9)
        else:
            ax.text(x0 + 1.6, y0 + 2.0, comp.label.upper(), fontsize=6.4, fontweight="bold", color=ink,
                    ha="left", va="center", zorder=0.9)

    hub_by_id = {h.id: h for h in spec.hubs}

    # 3. 轨道流向（正交 Metro 路线 + Economist 红/蓝/柔色空心/石板灰）
    for fl in spec.flows:
        sh, dh = hub_by_id[fl.source], hub_by_id[fl.target]
        p0, p1 = sh.xy, dh.xy
        raw_pts = route_octilinear(p0, p1, bend=fl.curve * 0.45)
        pts = fillet_polyline(raw_pts, radius=2.0)
        t0 = _box_exit_trim(pts[1][0] - p0[0], pts[1][1] - p0[1], HUB_W / 2, HUB_H / 2, pad=0.65)
        t1 = _box_exit_trim(p1[0] - pts[-2][0], p1[1] - pts[-2][1], HUB_W / 2, HUB_H / 2, pad=0.85)

        strong = fl.evidence_tier in ("L3", "L2")
        hollow = (fl.evidence_tier == "L1") and (fl.direction != 0)
        base_col = _flow_color(fl.direction, fl.evidence_tier, pal.up, pal.down, pal.edge_grey)
        if color_mode == "pathway" and fl.direction == 0:
            base_col = pathway_family_color(sh.icon, fallback_accent=pal.edge_grey)
        col = tint_color(base_col, 0.65, base=pal.card_bg) if hollow else base_col
        width = 3.2 if strong else (2.4 if hollow else (2.0 if fl.direction else 1.35))

        if fl.coupled:
            from ..metro import trim_polyline as _trim_pl
            draw_metro_coupled_halo(
                ax, [_trim_pl(pts, t0, t1)], color=pal.coupled_track, edge_color=pal.coupled_edge,
                width=7.4, zorder=1.8,
            )

        draw_metro_track(
            ax, pts, col,
            width=width, hollow=hollow, surface=pal.card_bg if hollow else pal.surface,
            chevrons=True, chevron_spacing=18.0, chevron_direction=1,
            unit_in=unit_in, trim0=t0, trim1=t1, corner_radius=0.0,
            halo=True, zorder=3.2 if strong else (2.7 if hollow else 2.2),
        )
        if fl.label:
            mid, _ = sample_polyline_at(pts, polyline_length(pts) * 0.5)
            ax.text(
                mid[0], mid[1] - 0.85, fl.label, fontsize=5.8, color=pal.ink2,
                ha="center", va="center", zorder=4.2,
                bbox=dict(boxstyle="square,pad=0.14", facecolor=pal.surface, edgecolor="none", alpha=0.92),
            )

    # 4. 锐利编辑部枢纽卡片（顶部粗色条 + Economist 几何矢量图标）
    for hub in spec.hubs:
        cx, cy = hub.xy
        x0, y0 = cx - HUB_W / 2, cy - HUB_H / 2
        col, _tint_bg, strong, dashed = _hub_border_and_tint(hub, pal.up, pal.down, pal.edge_grey, pal.card_bg)
        card_face = pal.card_bg
        card_ec = pal.ink if strong else (pal.ink2 if dashed else pal.card_edge)
        ax.add_patch(Rectangle(
            (x0, y0), HUB_W, HUB_H,
            facecolor=card_face, edgecolor=card_ec,
            lw=1.15 if strong else (0.90 if dashed else 0.65),
            linestyle=(0, (3.0, 1.8)) if dashed else "-",
            zorder=4.0,
        ))
        # 顶部经济学人粗语义色带（L2+ 实心高饱和；L1 柔色）
        top_bar_col = tint_color(col, 0.65, base=pal.card_bg) if dashed else col
        ax.add_patch(Rectangle(
            (x0, y0), HUB_W, 1.05 if strong else 0.78,
            facecolor=top_bar_col, edgecolor="none", zorder=4.25,
        ))

        draw_glyph(ax, hub.icon, (x0 + 3.4, cy + 0.35), style="economist", size=4.5, zorder=4.6)
        lbl_fs = 6.9 if len(hub.label) <= 15 else (6.3 if len(hub.label) <= 18 else 5.8)
        ax.text(x0 + 6.3, cy - 1.35, hub.label, fontsize=lbl_fs, fontweight="bold", color=pal.ink,
                ha="left", va="center", zorder=4.8)
        if hub.subtitle:
            sub_col = pal.up if hub.direction > 0 else (pal.down if hub.direction < 0 else pal.muted)
            ax.text(x0 + 6.4, cy + 1.55, hub.subtitle, fontsize=5.8, color=sub_col,
                    fontweight="medium" if hub.direction else "normal",
                    ha="left", va="center", zorder=4.8)
        if hub.evidence_tier in ("L1", "L2", "L3"):
            b_sc = 1.24
            bw = tier_badge_width(hub.evidence_tier, scale=b_sc, unit_in=unit_in)
            draw_tier_badge(
                ax, x0 + HUB_W - bw - 0.5, y0 + HUB_H - 0.1, hub.evidence_tier,
                zorder=5.2, scale=b_sc, preset="economist", unit_in=unit_in,
            )

    _draw_ga_bottom_panels(
        ax, spec,
        surface=pal.surface, ink=pal.ink, ink2=pal.ink2, muted=pal.muted,
        card_bg=pal.card_bg, card_edge=pal.card_edge, up=pal.up,
        preset="economist", unit_in=unit_in, pixel_mode=False, metro_mode=True,
    )


STYLE_RENDERERS: dict[str, Callable[..., None]] = {
    "economist": render_economist,
    "nature_editorial": render_nature_editorial,
    "editorial_ink": render_editorial_ink,
    "pixel_abstract": render_pixel_abstract,
    "metro_atlas": render_metro_atlas,
}

__all__ = [
    "GA_STYLES",
    "STYLE_RENDERERS",
    "ga_footer_height",
    "render_economist",
    "render_editorial_ink",
    "render_metro_atlas",
    "render_nature_editorial",
    "render_pixel_abstract",
]
