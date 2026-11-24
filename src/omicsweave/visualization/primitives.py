#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
可复用的科学矢量绘图原语（完全解耦于具体组学与集成业务模块）。

包含：
- 二次贝塞尔曲线、精确边界裁剪（de Casteljau）与隐形矢量箭头（stealth arrowhead）；
- 耦合通量轨道底带（coupled conduit track）：在空间上直观串联化学计量锁死的反应链；
- 无彩色三点证据阶梯徽章（``●●● Confirmed`` / ``●●○ Flux`` / ``●○○ Enz-pred``），把红/蓝高饱和色唯一留给涨跌方向；
- 左上角独立排名角标（``#1``）与双扣链环串联通量徽章（``Linked flux`` / ``Linked #1``）；
- 4 层多组学矢量微图标（Flux / Protein / mRNA / Metab）与纯矢量状态微格（支持 y 向下画布与标准 y 向上非等比坐标系）；
- 学术出版级分区卡片与双层膜细胞器容器。
"""
from __future__ import annotations

import math
from typing import Sequence

from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Ellipse, FancyBboxPatch, PathPatch, Polygon, Rectangle, Wedge
from matplotlib.path import Path as MplPath

from .palettes import ECONOMIST, ECONOMIST_RED, NATURE_EDITORIAL, SemanticPalette, tier_tokens, tint_color

Point = tuple[float, float]


def draw_economist_masthead(
    ax,
    x0: float,
    y_top: float,
    width: float,
    *,
    tag_w: float = 5.4,
    tag_h: float = 1.15,
    color: str = ECONOMIST_RED,
    rule_color: str = "#B8B8AF",
    y_up: bool = False,
    zorder: float = 8.0,
) -> None:
    """在数据坐标系顶部绘制《经济学人》标志性左上角红标与横贯细顶线。"""
    ax.plot(
        [x0, x0 + width],
        [y_top, y_top],
        color=rule_color,
        lw=0.85,
        solid_capstyle="butt",
        zorder=zorder,
    )
    rect_y = (y_top - tag_h) if y_up else y_top
    ax.add_patch(
        Rectangle(
            (x0, rect_y),
            tag_w,
            tag_h,
            facecolor=color,
            edgecolor="none",
            zorder=zorder + 0.1,
        )
    )


# ---------------------------------------------------------------------------
# 字体宽度估算（独立于具体图幅 UNIT）
# ---------------------------------------------------------------------------

def em_units(fontsize: float, unit_inches: float) -> float:
    """将字号（pt）换算为当前画布坐标单位。"""
    return float(fontsize) / 72.0 / float(unit_inches)


def estimate_text_width(text: str, fontsize: float, unit_inches: float, *, bold: bool = False) -> float:
    """估算字符串在画布坐标系下的宽度（用于徽章排版与自动避让）。"""
    scale = 1.10 if bold else 1.0
    total = 0.0
    for ch in str(text):
        if ord(ch) >= 0x2E80:
            total += 1.04
        elif ch.isupper() or ch in "#[]_WMG":
            total += 0.67 * scale
        else:
            total += 0.55 * scale
    return total * em_units(fontsize, unit_inches)


# ---------------------------------------------------------------------------
# 二次贝塞尔曲线与精确边界裁剪
# ---------------------------------------------------------------------------

def quad_point(p0: Point, c: Point, p1: Point, t: float) -> Point:
    u = 1.0 - t
    return (
        u * u * p0[0] + 2.0 * u * t * c[0] + t * t * p1[0],
        u * u * p0[1] + 2.0 * u * t * c[1] + t * t * p1[1],
    )


def quad_tangent(p0: Point, c: Point, p1: Point, t: float) -> Point:
    dx = 2.0 * (1.0 - t) * (c[0] - p0[0]) + 2.0 * t * (p1[0] - c[0])
    dy = 2.0 * (1.0 - t) * (c[1] - p0[1]) + 2.0 * t * (p1[1] - c[1])
    norm = math.hypot(dx, dy) or 1.0
    return dx / norm, dy / norm


def subdivide_quad(p0: Point, c: Point, p1: Point, t0: float, t1: float) -> tuple[Point, Point, Point]:
    """通过 de Casteljau 细分提取二次贝塞尔曲线的 ``[t0, t1]`` 子段。"""
    def lerp(a: Point, b: Point, t: float) -> Point:
        return a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t

    q0 = quad_point(p0, c, p1, t0)
    qc = lerp(c, p1, t0)
    q1 = p1
    if t1 >= 1.0 - 1e-6:
        return q0, qc, q1
    s = (t1 - t0) / max(1.0 - t0, 1e-6)
    r0 = q0
    rc = lerp(q0, qc, s)
    r1 = quad_point(q0, qc, q1, s)
    return r0, rc, r1


def find_t_from_start(p0: Point, c: Point, p1: Point, dist_target: float) -> float:
    if dist_target <= 0:
        return 0.0
    lo, hi = 0.0, 0.48
    for _ in range(14):
        mid = 0.5 * (lo + hi)
        if math.dist(quad_point(p0, c, p1, mid), p0) < dist_target:
            lo = mid
        else:
            hi = mid
    return hi


def find_t_from_end(p0: Point, c: Point, p1: Point, dist_target: float) -> float:
    if dist_target <= 0:
        return 1.0
    lo, hi = 0.52, 1.0
    for _ in range(14):
        mid = 0.5 * (lo + hi)
        if math.dist(quad_point(p0, c, p1, mid), p1) < dist_target:
            hi = mid
        else:
            lo = mid
    return lo


def sample_quad(p0: Point, c: Point, p1: Point, n: int = 12) -> list[Point]:
    return [quad_point(p0, c, p1, i / n) for i in range(n + 1)]


def draw_stealth_head(
    ax,
    tip: Point,
    direction: Point,
    color,
    size: float,
    zorder: float,
    *,
    halo: bool = False,
    surface: str = "#FFFFFF",
) -> None:
    """在数据坐标系中绘制利落的隐形矢量箭头头部（stealth arrowhead）。"""
    ux, uy = direction
    nx, ny = -uy, ux
    length = size
    half_w = size * 0.43
    notch = size * 0.24
    back_l = (tip[0] - ux * length + nx * half_w, tip[1] - uy * length + ny * half_w)
    back_r = (tip[0] - ux * length - nx * half_w, tip[1] - uy * length - ny * half_w)
    inner = (tip[0] - ux * (length - notch), tip[1] - uy * (length - notch))
    verts = [tip, back_l, inner, back_r]
    if halo:
        ax.add_patch(Polygon(
            verts, closed=True, facecolor=surface, edgecolor=surface, lw=1.0,
            joinstyle="round", zorder=zorder - 0.15,
        ))
    ax.add_patch(Polygon(verts, closed=True, facecolor=color, edgecolor="none", zorder=zorder))


def draw_quad_arrow(
    ax,
    p0: Point,
    c: Point,
    p1: Point,
    color,
    width: float,
    style,
    head_start: bool,
    head_end: bool,
    trim_start: float,
    trim_end: float,
    zorder: float,
    *,
    halo: bool = False,
    surface: str = "#FFFFFF",
    arrow_scale: float = 1.0,
) -> None:
    """绘制带精确边界裁剪、可选白色光晕与双向隐形箭头的二次贝塞尔曲线。"""
    chord = math.hypot(p1[0] - p0[0], p1[1] - p0[1])
    if chord < 1e-4:
        return
    arrow_len = max(0.82, min(1.45, (0.72 + width * 0.26) * arrow_scale))
    ts = find_t_from_start(p0, c, p1, trim_start)
    te = find_t_from_end(p0, c, p1, trim_end)
    if te <= ts + 0.04:
        ts, te = 0.2, 0.8
    tip_start = quad_point(p0, c, p1, ts)
    tip_end = quad_point(p0, c, p1, te)
    tan_start = quad_tangent(p0, c, p1, ts)
    tan_end = quad_tangent(p0, c, p1, te)
    ls = find_t_from_start(p0, c, p1, trim_start + (arrow_len * 0.65 if head_start else 0.0))
    le = find_t_from_end(p0, c, p1, trim_end + (arrow_len * 0.65 if head_end else 0.0))
    if le > ls + 0.02:
        q0, qc, q1 = subdivide_quad(p0, c, p1, ls, le)
        path = MplPath([q0, qc, q1], [MplPath.MOVETO, MplPath.CURVE3, MplPath.CURVE3])
        if halo:
            ax.add_patch(PathPatch(
                path, facecolor="none", edgecolor=surface, lw=width + 1.5,
                capstyle="round", zorder=zorder - 0.2,
            ))
        ax.add_patch(PathPatch(
            path, facecolor="none", edgecolor=color, lw=width, linestyle=style,
            capstyle="round", zorder=zorder,
        ))
    if head_end:
        draw_stealth_head(ax, tip_end, tan_end, color, arrow_len, zorder + 0.1, halo=halo, surface=surface)
    if head_start:
        draw_stealth_head(
            ax, tip_start, (-tan_start[0], -tan_start[1]), color, arrow_len, zorder + 0.1, halo=halo, surface=surface
        )


def draw_coupled_track(
    ax,
    curves: Sequence[dict],
    *,
    color: str = NATURE_EDITORIAL.coupled_track,
    edge_color: str = NATURE_EDITORIAL.coupled_edge,
    width: float = 7.2,
    zorder: float = 1.85,
) -> None:
    """在化学计量耦合通量组的反应边下方绘制柔和轨道包络底带，使共检链在视觉上浑然一体。"""
    for cv in curves:
        p0, c, p1 = cv["p0"], cv["c"], cv["p1"]
        ts = find_t_from_start(p0, c, p1, cv.get("trim0", 0.0) * 0.55)
        te = find_t_from_end(p0, c, p1, cv.get("trim1", 0.0) * 0.55)
        if te <= ts + 0.02:
            continue
        q0, qc, q1 = subdivide_quad(p0, c, p1, ts, te)
        path = MplPath([q0, qc, q1], [MplPath.MOVETO, MplPath.CURVE3, MplPath.CURVE3])
        # 外缘微描边 + 内部柔和底带
        ax.add_patch(PathPatch(
            path, facecolor="none", edgecolor=edge_color, lw=width + 0.9,
            alpha=0.22, capstyle="round", zorder=zorder - 0.02,
        ))
        ax.add_patch(PathPatch(
            path, facecolor="none", edgecolor=color, lw=width,
            alpha=0.82, capstyle="round", zorder=zorder,
        ))


# ---------------------------------------------------------------------------
# 语义化矢量徽章：三点证据阶梯、左上角排名角标、耦合链环 G-徽章
# ---------------------------------------------------------------------------
# 自解释无彩色三点证据阶梯徽章（●●● Confirmed / ●●○ Flux / ●○○ Enz-pred）
# ---------------------------------------------------------------------------

TIER_DISPLAY_LABELS: dict[str, str] = {
    "L3": "Confirmed",
    "L2": "Flux",
    "L1": "Enz-pred",
}


def tier_badge_text(level: str, *, label_style: str = "descriptive", label: str | None = None) -> str:
    if level not in ("L1", "L2", "L3"):
        return ""
    if label is not None and str(label).strip():
        return str(label).strip()
    if label_style == "code":
        return level
    return TIER_DISPLAY_LABELS[level]


def tier_badge_width(
    level: str,
    scale: float = 1.0,
    *,
    show_code: bool = True,
    unit_in: float = 0.108,
    label_style: str = "descriptive",
    fontsize: float | None = None,
    label: str | None = None,
) -> float:
    if level not in ("L1", "L2", "L3"):
        return 0.0
    if not show_code:
        return 1.92 * scale
    text = tier_badge_text(level, label_style=label_style, label=label)
    fs = float(fontsize) if fontsize is not None else 3.40 * scale
    dots_w = 1.64 * scale
    txt_w = estimate_text_width(text, fs, unit_in, bold=True) * 0.94
    return dots_w + txt_w + 0.42 * scale


def draw_tier_badge(
    ax,
    x: float,
    y: float,
    level: str,
    zorder: float,
    *,
    scale: float = 1.0,
    preset: str = "nature_editorial",
    show_code: bool = True,
    unit_in: float = 0.108,
    label_style: str = "descriptive",
    fontsize: float | None = None,
    y_scale: float = 1.0,
    label: str | None = None,
) -> float:
    """绘制无彩色三点证据阶梯徽章（``●●● <多示踪剂>`` / ``●●○ <示踪剂> Flux`` / ``●○○ Enz-pred``）。

    采用纯黑白灰材质编码证据强度，绝不与红/蓝涨跌方向冲突。支持通过 ``label`` 覆盖文字（如显示实际约束示踪剂短标）
    及通过 ``y_scale`` 适配非等比坐标轴，返回消耗的 x 宽度。
    """
    if level not in ("L1", "L2", "L3"):
        return 0.0
    tok = tier_tokens(preset)[level]
    text = tier_badge_text(level, label_style=label_style, label=label)
    fs = float(fontsize) if fontsize is not None else 3.40 * scale
    bw = tier_badge_width(
        level, scale=scale, show_code=show_code, unit_in=unit_in, label_style=label_style, fontsize=fs, label=label,
    )
    bh = max(1.14 * scale, em_units(fs, unit_in) * 1.46 if fontsize is not None else 0.0) * y_scale
    lw_scale = min(1.35, (fs / 3.40)) if fontsize is not None else scale
    is_econ = preset == "economist"
    rad = (0.12 if is_econ else 0.38) * scale
    ax.add_patch(FancyBboxPatch(
        (x, y - bh / 2), bw, bh,
        boxstyle=f"round,pad=0,rounding_size={rad}",
        mutation_aspect=y_scale,
        facecolor=tok["bg"], edgecolor=tok["ec"],
        lw=max(0.45, float(tok["lw"]) * lw_scale), linestyle=tok["ls"], zorder=zorder,
    ))
    n_on = int(tok["dots"])
    dot_rx = 0.16 * scale
    dot_ry = dot_rx * y_scale
    dot_step = 0.40 * scale
    dot_x0 = x + 0.42 * scale
    for idx in range(3):
        cx = dot_x0 + idx * dot_step
        filled = idx < n_on
        if filled:
            ax.add_patch(Ellipse(
                (cx, y), 2 * dot_rx, 2 * dot_ry,
                facecolor=tok["dot_on"], edgecolor="none", zorder=zorder + 0.12,
            ))
        else:
            ax.add_patch(Ellipse(
                (cx, y), 2 * dot_rx * 0.92, 2 * dot_ry * 0.92,
                facecolor="none", edgecolor=tok["dot_off"],
                lw=0.45 * lw_scale, zorder=zorder + 0.12,
            ))
    if show_code:
        ax.text(
            dot_x0 + 2 * dot_step + dot_rx + 0.24 * scale,
            y,
            text,
            fontsize=fs,
            fontweight="bold",
            ha="left",
            va="center",
            color=tok["fg"],
            zorder=zorder + 0.15,
        )
    return bw


def rank_tab_width(rank: int, scale: float = 1.0) -> float:
    return (1.68 if rank < 10 else 2.05) * scale


def draw_rank_tab(
    ax,
    card_left: float,
    card_top: float,
    rank: int,
    zorder: float,
    *,
    scale: float = 1.0,
    palette: SemanticPalette = NATURE_EDITORIAL,
) -> float:
    """在反应卡片左上角绘制独立 ``#n`` 排名角标（与右侧证据/耦合徽章彻底物理隔离）。"""
    bw = rank_tab_width(rank, scale)
    bh = 0.92 * scale
    top3 = rank <= 3
    is_econ = palette.name == "economist"
    bg = (ECONOMIST_RED if is_econ else palette.ink) if top3 else palette.card_bg
    ec = (ECONOMIST_RED if is_econ else palette.ink) if top3 else palette.ink2
    fg = "#FFFFFF" if (top3 and is_econ) else (palette.surface if top3 else palette.ink)
    y0 = card_top - bh * 0.62
    x0 = card_left - 0.05 * scale
    ax.add_patch(FancyBboxPatch(
        (x0, y0), bw, bh,
        boxstyle=f"round,pad=0,rounding_size={(0.10 if is_econ else 0.28) * scale}",
        facecolor=bg, edgecolor=ec, lw=(0.55 if top3 else 0.45) * scale, zorder=zorder,
    ))
    ax.text(
        x0 + bw / 2,
        y0 + bh / 2 + 0.02 * scale,
        f"#{rank}",
        fontsize=3.15 * scale,
        fontweight="bold",
        ha="center",
        va="center",
        color=fg,
        zorder=zorder + 0.1,
    )
    return bw


def coupled_badge_label(group_num: int, *, total_groups: int = 1, label: str | None = None) -> str:
    if label is not None:
        return label
    return "Linked flux" if total_groups <= 1 else f"Linked #{group_num}"


def coupled_badge_width(
    group_num: int,
    scale: float = 1.0,
    *,
    total_groups: int = 1,
    label: str | None = None,
    unit_in: float = 0.108,
    fontsize: float | None = None,
) -> float:
    txt = coupled_badge_label(group_num, total_groups=total_groups, label=label)
    fs = float(fontsize) if fontsize is not None else 3.30 * scale
    return 1.62 * scale + estimate_text_width(txt, fs, unit_in, bold=True) * 0.94 + 0.38 * scale


def draw_coupled_badge(
    ax,
    x: float,
    y: float,
    group_num: int,
    zorder: float,
    *,
    scale: float = 1.0,
    palette: SemanticPalette = NATURE_EDITORIAL,
    total_groups: int = 1,
    label: str | None = None,
    unit_in: float = 0.108,
    fontsize: float | None = None,
) -> float:
    """绘制清晰带双环锁链 Icon 与直白文字的串联通量徽章（``🔗 Linked flux`` / ``🔗 Linked #1``）。"""
    txt = coupled_badge_label(group_num, total_groups=total_groups, label=label)
    fs = float(fontsize) if fontsize is not None else 3.30 * scale
    bw = coupled_badge_width(
        group_num, scale=scale, total_groups=total_groups, label=label, unit_in=unit_in, fontsize=fs,
    )
    bh = max(1.14 * scale, em_units(fs, unit_in) * 1.46 if fontsize is not None else 0.0)
    is_econ = palette.name == "economist"
    ax.add_patch(FancyBboxPatch(
        (x, y - bh / 2), bw, bh,
        boxstyle=f"round,pad=0,rounding_size={(0.12 if is_econ else 0.36) * scale}",
        facecolor=palette.coupled_fill, edgecolor=palette.coupled_edge,
        lw=max(0.45, 0.48 * scale), zorder=zorder,
    ))
    # 清晰可读的双扣链环矢量图元
    lh = 0.42 * scale
    lw_link = 0.64 * scale
    lx0 = x + 0.28 * scale
    ax.add_patch(FancyBboxPatch(
        (lx0, y - lh / 2), lw_link, lh,
        boxstyle=f"round,pad=0,rounding_size={lh / 2}",
        facecolor="none", edgecolor=palette.coupled_edge, lw=max(0.45, 0.62 * scale), zorder=zorder + 0.1,
    ))
    ax.add_patch(FancyBboxPatch(
        (lx0 + 0.36 * scale, y - lh / 2), lw_link, lh,
        boxstyle=f"round,pad=0,rounding_size={lh / 2}",
        facecolor="none", edgecolor=palette.coupled_ink, lw=max(0.45, 0.62 * scale), zorder=zorder + 0.12,
    ))
    ax.text(
        x + 1.62 * scale,
        y + 0.02 * scale,
        txt,
        fontsize=fs,
        fontweight="bold",
        ha="left",
        va="center",
        color=palette.coupled_ink,
        zorder=zorder + 0.15,
    )
    return bw


def draw_sign_reversal_badge(
    ax,
    cx: float,
    cy: float,
    radius: float = 0.78,
    *,
    palette: SemanticPalette = NATURE_EDITORIAL,
    lw: float = 0.68,
    zorder: float = 4.5,
    y_scale: float = 1.0,
) -> None:
    """绘制两组间净通量方向反转（Net flux reversal）纯矢量琥珀圆徽章（内部为上下反向双箭头 ``⇄``）。

    与热力学可逆反应（轨道双向箭头 ``< >``）严格区分：仅当实验组与对照组净通量均值异号时出现。
    """
    rx = float(radius)
    ry = rx * y_scale
    ax.add_patch(Ellipse(
        (cx, cy), 2.0 * rx, 2.0 * ry,
        facecolor=palette.sign_change_bg,
        edgecolor=palette.sign_change_edge,
        lw=lw,
        zorder=zorder,
    ))
    arrow_lw = max(0.55, lw * 0.96)
    ink = palette.sign_change_ink

    # 上排向右箭头（→）
    y_top = cy - ry * 0.27
    ax.add_line(Line2D(
        [cx - rx * 0.54, cx + rx * 0.44], [y_top, y_top],
        color=ink, lw=arrow_lw, solid_capstyle="round", zorder=zorder + 0.10,
    ))
    head_r = MplPath(
        [(cx + rx * 0.16, y_top - ry * 0.24), (cx + rx * 0.56, y_top), (cx + rx * 0.16, y_top + ry * 0.24)],
        [MplPath.MOVETO, MplPath.LINETO, MplPath.LINETO],
    )
    ax.add_patch(PathPatch(
        head_r, facecolor="none", edgecolor=ink, lw=arrow_lw,
        capstyle="round", joinstyle="miter", zorder=zorder + 0.12,
    ))

    # 下排向左箭头（←）
    y_bot = cy + ry * 0.27
    ax.add_line(Line2D(
        [cx + rx * 0.54, cx - rx * 0.44], [y_bot, y_bot],
        color=ink, lw=arrow_lw, solid_capstyle="round", zorder=zorder + 0.10,
    ))
    head_l = MplPath(
        [(cx - rx * 0.16, y_bot - ry * 0.24), (cx - rx * 0.56, y_bot), (cx - rx * 0.16, y_bot + ry * 0.24)],
        [MplPath.MOVETO, MplPath.LINETO, MplPath.LINETO],
    )
    ax.add_patch(PathPatch(
        head_l, facecolor="none", edgecolor=ink, lw=arrow_lw,
        capstyle="round", joinstyle="miter", zorder=zorder + 0.12,
    ))


# ---------------------------------------------------------------------------
# 极简纯矢量组学微图标（通量 / 蛋白 / mRNA / 代谢物）与单格状态绘制原语
# ---------------------------------------------------------------------------

AXIS_CODES = ("F", "C", "R", "M")
OMICS_LAYER_ICONS: tuple[str, ...] = ("flux", "protein", "mrna", "metabolite")
AXIS_LAYER_ICONS: tuple[str, ...] = ("flux", "protein", "regulation", "metabolite")


def draw_layer_icon(
    ax,
    kind: str,
    cx: float,
    cy: float,
    size: float,
    color: str,
    *,
    lw: float = 0.62,
    zorder: float = 6.2,
    y_scale: float = 1.0,
    clip_on: bool = True,
) -> None:
    """绘制高可读极简纯矢量微图标：``flux`` / ``protein`` / ``mrna`` / ``metabolite`` / ``regulation``。"""
    rx = size * 0.46
    ry = rx * y_scale
    if kind == "flux":
        # 双层水平流向箭头（通量 / Flux）
        for dy in (-ry * 0.46, ry * 0.46):
            y = cy + dy
            ax.add_line(Line2D(
                [cx - rx * 0.92, cx + rx * 0.85], [y, y],
                color=color, lw=lw, solid_capstyle="round", zorder=zorder, clip_on=clip_on,
            ))
            head = MplPath(
                [(cx + rx * 0.42, y - ry * 0.34), (cx + rx * 0.95, y), (cx + rx * 0.42, y + ry * 0.34)],
                [MplPath.MOVETO, MplPath.LINETO, MplPath.LINETO],
            )
            ax.add_patch(PathPatch(
                head, facecolor="none", edgecolor=color, lw=lw * 0.92,
                capstyle="round", joinstyle="miter", zorder=zorder + 0.02, clip_on=clip_on,
            ))
    elif kind == "protein":
        # 经典生化酶活性口袋球（Protein / Enzyme globule with cleft）
        if abs(y_scale - 1.0) < 1e-4:
            ax.add_patch(Wedge(
                (cx, cy), rx * 0.88, 34, 326,
                facecolor=color, edgecolor="none", zorder=zorder, clip_on=clip_on,
            ))
        else:
            verts = [(cx, cy)] + [
                (
                    cx + rx * 0.88 * math.cos(math.radians(34.0 + 292.0 * k / 28.0)),
                    cy + ry * 0.88 * math.sin(math.radians(34.0 + 292.0 * k / 28.0)),
                )
                for k in range(29)
            ]
            ax.add_patch(Polygon(verts, closed=True, facecolor=color, edgecolor="none", zorder=zorder, clip_on=clip_on))
    elif kind == "mrna":
        # 单链波浪 mRNA + 5' 端帽圆点
        xs = [cx - rx * 0.92 + (1.84 * rx) * (k / 14.0) for k in range(15)]
        ys = [cy + ry * 0.48 * math.sin(k * 2.0 * math.pi / 7.0) for k in range(15)]
        ax.add_line(Line2D(xs, ys, color=color, lw=lw * 1.08, solid_capstyle="round", zorder=zorder, clip_on=clip_on))
        ax.add_patch(Ellipse(
            (xs[0], ys[0]), 2 * rx * 0.22, 2 * ry * 0.22,
            facecolor=color, edgecolor="none", zorder=zorder + 0.03, clip_on=clip_on,
        ))
    elif kind == "metabolite":
        # 规则六边形小分子环（Metabolite ring）
        verts = [
            (cx + rx * 0.88 * math.cos(math.radians(60 * k)), cy + ry * 0.88 * math.sin(math.radians(60 * k)))
            for k in range(6)
        ]
        ax.add_patch(Polygon(verts, closed=True, facecolor="none", edgecolor=color, lw=lw * 1.05, zorder=zorder, clip_on=clip_on))
    elif kind == "regulation":
        # 菱形别构调控开关（Regulation diamond + center dot）
        verts = [(cx, cy - ry * 0.88), (cx + rx * 0.88, cy), (cx, cy + ry * 0.88), (cx - rx * 0.88, cy)]
        ax.add_patch(Polygon(verts, closed=True, facecolor="none", edgecolor=color, lw=lw, zorder=zorder, clip_on=clip_on))
        ax.add_patch(Ellipse(
            (cx, cy), 2 * rx * 0.22, 2 * ry * 0.22,
            facecolor=color, edgecolor="none", zorder=zorder + 0.02, clip_on=clip_on,
        ))


def draw_status_cell(
    ax,
    cell_x: float,
    cells_y: float,
    cell: float,
    state: str,
    direction: int = 0,
    *,
    palette: SemanticPalette = NATURE_EDITORIAL,
    font_scale: float = 1.0,
    zorder: float = 6.1,
    draw_box: bool = False,
    y_up: bool = False,
    y_scale: float = 1.0,
) -> None:
    """绘制单个多组学状态微格（地图卡片、矩阵图与图例共享同一函数，保证 100% 视觉一致）。

    - ``changed`` / ``established``：红/蓝实心三角（显著变化）；
    - ``suggestive`` / ``nominal``：红/蓝空心粗边三角（提示趋势；组学为名义 p<0.05，通量为 |Hedges' g|≥1 且方向稳定）；
    - ``mixed`` / ``conflicting``：红蓝对角斜拼方格（多成员方向冲突）；
    - ``unchanged`` / ``equivalent``：深灰实线短横 ``—``（通过等效性检验，确证无变化）；
    - ``undetermined`` / ``no_evidence``：中灰虚线短横 ``╌``（可辨识/已测到，但效应小、不稳定或证据不足）；
    - ``not_measured`` / ``not_infused`` / ``no_knowledge``：浅灰底 + 淡 45° 斜纹 ``/``；
    - ``not_identifiable``：灰蓝禁用约束圆图标 ``⊘``（替代难以理解的 ``×``）；
    - ``not_detected``：浅灰虚线空心圆 ``◌``（替代 ``?``）；
    - ``qc_blocked`` / ``insufficient_n``：低对比度冷灰双斜线 ``//``（低样本量或无可评估数据，避免满图高显著度警示干扰）。

    ``y_up=True`` 时适配标准 Matplotlib y 轴向上坐标系；``y_scale`` 为 ``x_unit_in / y_unit_in`` 修正纵横比。
    """
    cell_y = cell * y_scale
    if draw_box:
        ax.add_patch(FancyBboxPatch(
            (cell_x, cells_y), cell, cell_y,
            boxstyle=f"round,pad=0,rounding_size={0.14 * (cell / 0.98)}",
            mutation_aspect=y_scale,
            facecolor=palette.surface, edgecolor=palette.card_edge, lw=0.45, zorder=zorder - 0.08,
        ))

    cx = cell_x + cell / 2.0
    cy = cells_y + cell_y / 2.0
    up_s = 1.0 if y_up else -1.0
    st = str(state or "not_measured")
    dr = int(direction or 0)

    if st in ("changed", "established", "suggestive", "nominal") and dr:
        col = palette.up if dr > 0 else palette.down
        solid = st in ("changed", "established")
        rx = cell * 0.31
        ry = cell_y * 0.31
        if dr > 0:
            verts = [
                (cx, cy + up_s * ry),
                (cx + rx * 0.95, cy - up_s * ry * 0.85),
                (cx - rx * 0.95, cy - up_s * ry * 0.85),
            ]
        else:
            verts = [
                (cx, cy - up_s * ry),
                (cx + rx * 0.95, cy + up_s * ry * 0.85),
                (cx - rx * 0.95, cy + up_s * ry * 0.85),
            ]
        ax.add_patch(Polygon(
            verts, closed=True,
            facecolor=col if solid else "none",
            edgecolor=col,
            lw=(0.68 if solid else 0.78) * font_scale,
            zorder=zorder + 0.15,
        ))
    elif st in ("mixed", "conflicting"):
        pad_x = cell * 0.19
        pad_y = cell_y * 0.19
        x_l, x_r = cell_x + pad_x, cell_x + cell - pad_x
        y_top = (cells_y + cell_y - pad_y) if y_up else (cells_y + pad_y)
        y_bot = (cells_y + pad_y) if y_up else (cells_y + cell_y - pad_y)
        ax.add_patch(Polygon(
            [(x_l, y_top), (x_r, y_top), (x_l, y_bot)], closed=True,
            facecolor=palette.up, edgecolor="none", alpha=0.84, zorder=zorder + 0.12,
        ))
        ax.add_patch(Polygon(
            [(x_r, y_top), (x_r, y_bot), (x_l, y_bot)], closed=True,
            facecolor=palette.down, edgecolor="none", alpha=0.84, zorder=zorder + 0.12,
        ))
    elif st in ("unchanged", "equivalent"):
        half = cell * 0.24
        ax.add_line(Line2D(
            [cx - half, cx + half], [cy, cy],
            color=palette.ink2,
            lw=0.86 * font_scale, solid_capstyle="round", zorder=zorder + 0.12,
        ))
    elif st in ("undetermined", "no_evidence"):
        half = cell * 0.27
        gap = cell * 0.06
        for x0_seg, x1_seg in ((cx - half, cx - gap), (cx + gap, cx + half)):
            ax.add_line(Line2D(
                [x0_seg, x1_seg], [cy, cy],
                color=palette.muted,
                lw=0.78 * font_scale, solid_capstyle="round", zorder=zorder + 0.12,
            ))
    elif st in ("not_measured", "not_infused", "no_knowledge"):
        pad_x = cell * 0.07
        pad_y = cell_y * 0.07
        ax.add_patch(Rectangle(
            (cell_x + pad_x, cells_y + pad_y), cell - 2 * pad_x, cell_y - 2 * pad_y,
            facecolor=palette.missing_fill, edgecolor="none", zorder=zorder,
        ))
        y_lo = cells_y + cell_y * (0.24 if y_up else 0.76)
        y_hi = cells_y + cell_y * (0.76 if y_up else 0.24)
        ax.add_line(Line2D(
            [cell_x + cell * 0.24, cell_x + cell * 0.76],
            [y_lo, y_hi],
            color=palette.missing_hatch, lw=0.45 * font_scale, zorder=zorder + 0.06,
        ))
    elif st in ("not_identifiable", "bound_pinned"):
        # 灰蓝不可辨识/受约束斜杠圆 ⊘
        rx = cell * 0.26
        ry = cell_y * 0.26
        ax.add_patch(Ellipse(
            (cx, cy), 2 * rx, 2 * ry, facecolor="none", edgecolor=palette.muted,
            lw=0.62 * font_scale, zorder=zorder + 0.12,
        ))
        dx, dy = rx * 0.70, ry * 0.70
        ax.add_line(Line2D(
            [cx - dx, cx + dx], [cy - up_s * dy, cy + up_s * dy],
            color=palette.muted, lw=0.62 * font_scale, zorder=zorder + 0.14,
        ))
    elif st == "not_detected":
        # 浅灰虚线空心圆 ◌（低于检出限/未检出）
        rx = cell * 0.25
        ry = cell_y * 0.25
        ax.add_patch(Ellipse(
            (cx, cy), 2 * rx, 2 * ry, facecolor="none", edgecolor=palette.muted,
            lw=0.58 * font_scale, linestyle=(0, (1.4, 1.2)), zorder=zorder + 0.12,
        ))
    elif st in ("qc_blocked", "insufficient_n"):
        # 低对比度冷灰双斜线 //（低样本量或无可评估数据，避免满图高显著度琥珀感叹号干扰阅读）
        pad_x = cell * 0.07
        pad_y = cell_y * 0.07
        ax.add_patch(Rectangle(
            (cell_x + pad_x, cells_y + pad_y), cell - 2 * pad_x, cell_y - 2 * pad_y,
            facecolor=palette.missing_fill, edgecolor="none", zorder=zorder,
        ))
        y_lo = cells_y + cell_y * (0.25 if y_up else 0.75)
        y_hi = cells_y + cell_y * (0.75 if y_up else 0.25)
        for dx_off in (-0.09 * cell, 0.09 * cell):
            ax.add_line(Line2D(
                [cell_x + cell * 0.27 + dx_off, cell_x + cell * 0.73 + dx_off],
                [y_lo, y_hi],
                color=palette.muted, lw=0.48 * font_scale, solid_capstyle="round", zorder=zorder + 0.08,
            ))
    else:
        half = cell * 0.22
        ax.add_line(Line2D(
            [cx - half, cx + half], [cy, cy],
            color=palette.muted, lw=0.72 * font_scale, solid_capstyle="round", zorder=zorder + 0.12,
        ))


def draw_axis_strip(
    ax,
    x0: float,
    y_center: float,
    states: Sequence[tuple[str, int]],
    *,
    cell: float = 1.02,
    show_header_letters: bool = True,
    header_icons: Sequence[str] | None = OMICS_LAYER_ICONS,
    header_codes: Sequence[str] = AXIS_CODES,
    header_h: float | None = None,
    zorder: float = 6.0,
    palette: SemanticPalette = NATURE_EDITORIAL,
    font_scale: float = 1.0,
) -> tuple[float, float]:
    """绘制四格多组学/分轴微条（默认表头使用纯矢量微图标：通量 / 蛋白 / mRNA / 代谢物）。

    Returns ``(total_width, total_height)``.
    """
    n_cols = len(states)
    hdr_h = header_h if header_h is not None else ((0.60 if show_header_letters else 0.14) * (cell / 1.02))
    total_w = n_cols * cell
    total_h = cell + hdr_h
    top_y = y_center - total_h / 2.0
    cells_y = top_y + hdr_h

    # 外框容器
    is_econ = palette.name == "economist"
    strip_bg = palette.card_bg if is_econ else palette.surface
    ax.add_patch(FancyBboxPatch(
        (x0, top_y), total_w, total_h,
        boxstyle=f"round,pad=0,rounding_size={(0.08 if is_econ else 0.18) * (cell / 0.98)}",
        facecolor=strip_bg, edgecolor=palette.ink2 if is_econ else palette.card_edge,
        lw=0.45 if is_econ else 0.40, zorder=zorder,
    ))
    # 中性表头底带（纯灰阶，不抢占红蓝语义）
    hdr_bg = tint_color(palette.ink2, 0.11 if is_econ else 0.08, base=strip_bg)
    ax.add_patch(Rectangle(
        (x0 + 0.03, top_y + 0.03), total_w - 0.06, hdr_h - 0.03,
        facecolor=hdr_bg, edgecolor="none", zorder=zorder + 0.05,
    ))

    for i, (state, direction) in enumerate(states):
        cell_x = x0 + i * cell
        cx = cell_x + cell / 2.0

        if show_header_letters:
            if header_icons and i < len(header_icons):
                draw_layer_icon(
                    ax, header_icons[i], cx, top_y + hdr_h * 0.52,
                    size=min(cell * 0.66, hdr_h * 0.86),
                    color=palette.ink2, lw=0.58 * font_scale, zorder=zorder + 0.16,
                )
            elif i < len(header_codes):
                ax.text(
                    cx, top_y + hdr_h * 0.54, header_codes[i],
                    fontsize=2.95 * font_scale, fontweight="bold",
                    ha="center", va="center", color=palette.ink2, zorder=zorder + 0.15,
                )
        if i > 0:
            ax.plot(
                [cell_x, cell_x], [top_y + 0.05, cells_y + cell - 0.04],
                color=palette.hairline, lw=0.30, zorder=zorder + 0.12,
            )

        draw_status_cell(
            ax, cell_x, cells_y, cell, state, direction,
            palette=palette, font_scale=font_scale, zorder=zorder + 0.10, draw_box=False,
        )

    return total_w, total_h


# ---------------------------------------------------------------------------
# 分区底板卡片与细胞器双层膜容器
# ---------------------------------------------------------------------------

def draw_zone_card(
    ax,
    box: tuple[float, float, float, float],
    accent: str,
    *,
    palette: SemanticPalette = NATURE_EDITORIAL,
    tint_fill: float = 0.045,
    tint_edge: float = 0.26,
    rounding: float = 1.6,
    lw: float = 0.65,
    zorder: float = 0.2,
) -> None:
    x0, y0, x1, y1 = box
    is_dark = palette.name == "slide_dark"
    is_econ = palette.name == "economist"
    if is_econ:
        face = tint_color(accent, max(0.025, tint_fill * 0.55), base=palette.card_bg)
        edge = palette.card_edge
        ax.add_patch(FancyBboxPatch(
            (x0, y0), x1 - x0, y1 - y0,
            boxstyle="round,pad=0,rounding_size=0.22",
            facecolor=face, edgecolor=edge, lw=0.55, zorder=zorder,
        ))
        # 《经济学人》多栏图表顶部分隔细黑线
        ax.plot([x0, x1], [y0, y0], color=palette.ink, lw=0.85, solid_capstyle="butt", zorder=zorder + 0.12)
        return
    face = tint_color(accent, tint_fill * (2.2 if is_dark else 1.0), base=palette.surface)
    edge = tint_color(accent, tint_edge * (1.3 if is_dark else 1.0), base=palette.surface)
    ax.add_patch(FancyBboxPatch(
        (x0, y0), x1 - x0, y1 - y0,
        boxstyle=f"round,pad=0,rounding_size={rounding}",
        facecolor=face, edgecolor=edge, lw=lw, zorder=zorder,
    ))


def draw_zone_header(
    ax,
    xy: Point,
    label: str,
    accent: str,
    *,
    palette: SemanticPalette = NATURE_EDITORIAL,
    fontsize: float = 5.4,
    zorder: float = 0.9,
) -> None:
    lx, ly = xy
    if palette.name == "economist":
        # 经济学人锐利方角编辑部短标 + 黑色无衬线粗体分区标题（克制用色，不抢占红/蓝数据焦点）
        bar_col = accent or palette.ink
        ax.add_patch(Rectangle(
            (lx - 0.15, ly - 0.54), 0.36, 1.08,
            facecolor=bar_col, edgecolor="none", zorder=zorder,
        ))
        ax.text(
            lx + 0.58, ly, label,
            fontsize=fontsize, fontweight="bold", color=palette.ink,
            va="center", ha="left", zorder=zorder + 0.05,
        )
        return
    bar_col = accent or palette.ink2
    ax.add_patch(FancyBboxPatch(
        (lx - 0.15, ly - 0.62), 0.42, 1.24,
        boxstyle="round,pad=0,rounding_size=0.18",
        facecolor=bar_col, edgecolor="none", zorder=zorder,
    ))
    ax.text(
        lx + 0.68, ly, label,
        fontsize=fontsize, fontweight="bold", color=palette.ink,
        va="center", ha="left", zorder=zorder + 0.05,
    )


def draw_organelle_box(
    ax,
    box: tuple[float, float, float, float],
    label: str,
    kind: str,
    *,
    palette: SemanticPalette = NATURE_EDITORIAL,
    label_pos: str = "top_right",
    fontsize: float = 5.0,
    zorder: float = 0.5,
) -> None:
    """绘制双层膜细胞器区室（``kind="p"`` 过氧化物酶体，其余为线粒体基质）。"""
    x0, y0, x1, y1 = box
    if kind == "p":
        face, edge, inner, ink = palette.perox_face, palette.perox_edge, palette.perox_inner, palette.perox_ink
        r_out, r_in, inset = 2.0, 1.5, 0.50
        dash = (0, (2.2, 1.8))
    else:
        face, edge, inner, ink = palette.mito_face, palette.mito_edge, palette.mito_inner, palette.mito_ink
        r_out, r_in, inset = 3.2, 2.6, 0.65
        dash = (0, (3.2, 2.0))

    ax.add_patch(FancyBboxPatch(
        (x0, y0), x1 - x0, y1 - y0,
        boxstyle=f"round,pad=0,rounding_size={r_out}",
        facecolor=face, edgecolor=edge, lw=1.20, zorder=zorder,
    ))
    ax.add_patch(FancyBboxPatch(
        (x0 + inset, y0 + inset), (x1 - x0) - 2 * inset, (y1 - y0) - 2 * inset,
        boxstyle=f"round,pad=0,rounding_size={r_in}",
        facecolor="none", edgecolor=inner, lw=0.55, linestyle=dash, zorder=zorder + 0.05,
    ))
    if not label:
        return
    if label_pos == "bottom_left":
        lx, ly, ha = x0 + 1.5, y1 - 1.55, "left"
    elif label_pos == "top_left":
        lx, ly, ha = x0 + 1.5, y0 + 1.55, "left"
    elif label_pos == "bottom_right":
        lx, ly, ha = x1 - 1.5, y1 - 1.55, "right"
    else:
        lx, ly, ha = x1 - 1.5, y0 + 1.45, "right"
    ax.text(
        lx, ly, label,
        fontsize=fontsize, fontweight="bold", ha=ha, va="center", color=ink,
        zorder=zorder + 0.20,
        bbox=dict(boxstyle="round,pad=0.18", facecolor=face, edgecolor="none", alpha=0.88),
    )


__all__ = [
    "AXIS_CODES",
    "AXIS_LAYER_ICONS",
    "OMICS_LAYER_ICONS",
    "Point",
    "TIER_DISPLAY_LABELS",
    "coupled_badge_label",
    "coupled_badge_width",
    "draw_axis_strip",
    "draw_coupled_badge",
    "draw_coupled_track",
    "draw_economist_masthead",
    "draw_layer_icon",
    "draw_organelle_box",
    "draw_quad_arrow",
    "draw_rank_tab",
    "draw_sign_reversal_badge",
    "draw_status_cell",
    "draw_stealth_head",
    "draw_tier_badge",
    "draw_zone_card",
    "draw_zone_header",
    "em_units",
    "estimate_text_width",
    "find_t_from_end",
    "find_t_from_start",
    "quad_point",
    "quad_tangent",
    "rank_tab_width",
    "sample_quad",
    "subdivide_quad",
    "tier_badge_text",
    "tier_badge_width",
]
