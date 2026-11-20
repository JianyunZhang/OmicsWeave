#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
通用 Beck 地铁拓扑几何与代谢轨道图元库（零业务耦合，可供任意项目复用）。

提取自伦敦地铁图 / Metabolic Metro Map 的核心几何语法，并适配学术出版规范：
- 八向圆角折线寻路（0° / 45° / 90° + 光滑圆角倒角）：让高密度网络呈现出极高的秩序感；
- 实心主轨 vs 空心白芯双轨（Solid Track vs Hollow Double-Rail）：天然适配
  「实测通量支撑（L2/L3） vs 仅酶学/转录组预测（L1）」的视觉分层；
- 内嵌方向小箭头（Inline Chevrons）：直接刻在轨道内部，使正逆平行双轨紧贴并排而不打架；
- 圆角跑道闭环（Stadium Loop）：清晰区分循环通路（如 TCA / PPP / β-Oxidation）与线性干线；
- 双层圆环换乘站（Interchange Station）与极淡浅灰区室底板（Metro Zone Plate）。

本模块不依赖 ``omicsweave.integration``、``contract``、``pydantic`` 或 ``seaborn``。
"""
from __future__ import annotations

import math
from typing import Sequence

from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyBboxPatch, PathPatch, Polygon
from matplotlib.path import Path as MplPath

from .palettes import tint_color

Point = tuple[float, float]


# ---------------------------------------------------------------------------
# 八向（0° / 45° / 90°）折线路由、圆角倒角与弧长工具
# ---------------------------------------------------------------------------

def _clean_polyline(points: Sequence[Point], tol: float = 1e-4) -> list[Point]:
    """去除相邻重复点与共线冗余中间点。"""
    dedup: list[Point] = []
    for p in points:
        pt = (float(p[0]), float(p[1]))
        if not dedup or math.dist(dedup[-1], pt) > tol:
            dedup.append(pt)
    if len(dedup) <= 2:
        return dedup
    out: list[Point] = [dedup[0]]
    for i in range(1, len(dedup) - 1):
        a, b, c = out[-1], dedup[i], dedup[i + 1]
        v1 = (b[0] - a[0], b[1] - a[1])
        v2 = (c[0] - b[0], c[1] - b[1])
        n1 = math.hypot(*v1)
        n2 = math.hypot(*v2)
        if n1 < tol or n2 < tol:
            continue
        cross = abs(v1[0] * v2[1] - v1[1] * v2[0]) / (n1 * n2)
        dot = (v1[0] * v2[0] + v1[1] * v2[1]) / (n1 * n2)
        if cross < 1e-3 and dot > 0.0:
            continue
        out.append(b)
    out.append(dedup[-1])
    return out


def _octilinear_segment(
    p0: Point,
    p1: Point,
    *,
    bend: float = 0.0,
    prefer: str = "hv_first",
) -> list[Point]:
    """生成两点之间的单段八向折线（0° / 45° / 90°），支持平行偏置肩线。"""
    dx = p1[0] - p0[0]
    dy = p1[1] - p0[1]
    adx, ady = abs(dx), abs(dy)
    length = math.hypot(dx, dy)
    if length < 1e-6:
        return [p0, p1]

    # 已近似处于 0° / 90° / 45° 方向
    if adx < 0.08 or ady < 0.08 or abs(adx - ady) < 0.08:
        if abs(bend) > 0.12 and length > 1.2:
            ux, uy = dx / length, dy / length
            nx, ny = -uy, ux
            shoulder = min(abs(bend), length * 0.28)
            w1 = (p0[0] + ux * shoulder + nx * bend, p0[1] + uy * shoulder + ny * bend)
            w2 = (p1[0] - ux * shoulder + nx * bend, p1[1] - uy * shoulder + ny * bend)
            return _clean_polyline([p0, w1, w2, p1])
        return [p0, p1]

    sx = 1.0 if dx > 0 else -1.0
    sy = 1.0 if dy > 0 else -1.0

    # 计算两种单弯八向方案的拐点
    if adx > ady:
        w_hv = (p1[0] - sx * ady, p0[1])
        w_diag = (p0[0] + sx * ady, p1[1])
    else:
        w_hv = (p0[0], p1[1] - sy * adx)
        w_diag = (p1[0], p0[1] + sy * adx)

    if abs(bend) > 0.12:
        # 让正负 bend 与贝塞尔法向 (-dy, dx) 方向保持完全一致，使双向平行反应自然分列两侧
        nx, ny = -dy / length, dx / length
        mid_chord = ((p0[0] + p1[0]) * 0.5, (p0[1] + p1[1]) * 0.5)
        dot_hv = (w_hv[0] - mid_chord[0]) * nx + (w_hv[1] - mid_chord[1]) * ny
        use_hv = (dot_hv * bend) >= 0.0
        w = w_hv if use_hv else w_diag
        return _clean_polyline([p0, w, p1])

    if prefer == "centered":
        if adx > ady:
            half = (adx - ady) * 0.5
            w1 = (p0[0] + sx * half, p0[1])
            w2 = (p1[0] - sx * half, p1[1])
        else:
            half = (ady - adx) * 0.5
            w1 = (p0[0], p0[1] + sy * half)
            w2 = (p1[0], p1[1] - sy * half)
        return _clean_polyline([p0, w1, w2, p1])

    if prefer == "diag_first":
        return _clean_polyline([p0, w_diag, p1])
    return _clean_polyline([p0, w_hv, p1])


def route_octilinear(
    p0: Point,
    p1: Point,
    *,
    bend: float = 0.0,
    prefer: str = "hv_first",
    via: Sequence[Point] = (),
) -> list[Point]:
    """生成经过可选中间点 ``via`` 的八向地铁折线。"""
    if not via:
        return _octilinear_segment(p0, p1, bend=bend, prefer=prefer)
    anchors: list[Point] = [p0, *via, p1]
    out: list[Point] = [p0]
    for i in range(len(anchors) - 1):
        seg_pref = prefer if i == 0 else ("diag_first" if prefer == "hv_first" else "hv_first")
        seg = _octilinear_segment(anchors[i], anchors[i + 1], bend=0.0, prefer=seg_pref)
        out.extend(seg[1:])
    return _clean_polyline(out)


def fillet_polyline(
    points: Sequence[Point],
    radius: float = 1.25,
    *,
    arc_steps: int = 6,
    closed: bool = False,
) -> list[Point]:
    """对折线每个内角施加二次贝塞尔光滑倒角，呈现经典地铁圆角转弯。"""
    pts = _clean_polyline(points)
    if len(pts) < 3 or radius <= 1e-4:
        return list(pts)

    n = len(pts)
    out: list[Point] = []
    if not closed:
        out.append(pts[0])

    indices = range(n) if closed else range(1, n - 1)
    for i in indices:
        a = pts[(i - 1) % n]
        b = pts[i]
        c = pts[(i + 1) % n]
        la = math.dist(a, b)
        lc = math.dist(b, c)
        if la < 1e-5 or lc < 1e-5:
            out.append(b)
            continue
        u1 = ((a[0] - b[0]) / la, (a[1] - b[1]) / la)
        u2 = ((c[0] - b[0]) / lc, (c[1] - b[1]) / lc)
        dot = u1[0] * u2[0] + u1[1] * u2[1]
        if dot < -0.996:
            out.append(b)
            continue
        cut = min(radius, 0.42 * la, 0.42 * lc)
        p_in = (b[0] + u1[0] * cut, b[1] + u1[1] * cut)
        p_out = (b[0] + u2[0] * cut, b[1] + u2[1] * cut)
        for step in range(arc_steps + 1):
            t = step / arc_steps
            u = 1.0 - t
            qx = u * u * p_in[0] + 2.0 * u * t * b[0] + t * t * p_out[0]
            qy = u * u * p_in[1] + 2.0 * u * t * b[1] + t * t * p_out[1]
            out.append((qx, qy))

    if not closed:
        out.append(pts[-1])
    elif out:
        out.append(out[0])
    return _clean_polyline(out, tol=1e-5)


def polyline_arclengths(points: Sequence[Point]) -> list[float]:
    """返回折线各顶点的累计弧长列表。"""
    if not points:
        return []
    cum = [0.0]
    for i in range(1, len(points)):
        cum.append(cum[-1] + math.dist(points[i - 1], points[i]))
    return cum


def polyline_length(points: Sequence[Point]) -> float:
    cum = polyline_arclengths(points)
    return cum[-1] if cum else 0.0


def sample_polyline_at(points: Sequence[Point], dist: float) -> tuple[Point, Point]:
    """按弧长采样折线上的坐标与单位切向量 ``(point, unit_tangent)``。"""
    if len(points) < 2:
        p = points[0] if points else (0.0, 0.0)
        return p, (1.0, 0.0)
    cum = polyline_arclengths(points)
    total = cum[-1]
    if total < 1e-9:
        return points[0], (1.0, 0.0)
    target = max(0.0, min(total, float(dist)))
    for i in range(len(points) - 1):
        s0, s1 = cum[i], cum[i + 1]
        if target <= s1 + 1e-9 or i == len(points) - 2:
            seg_len = max(s1 - s0, 1e-9)
            t = max(0.0, min(1.0, (target - s0) / seg_len))
            a, b = points[i], points[i + 1]
            px = a[0] + (b[0] - a[0]) * t
            py = a[1] + (b[1] - a[1]) * t
            tx = (b[0] - a[0]) / seg_len
            ty = (b[1] - a[1]) / seg_len
            return (px, py), (tx, ty)
    a, b = points[-2], points[-1]
    seg_len = math.dist(a, b) or 1.0
    return points[-1], ((b[0] - a[0]) / seg_len, (b[1] - a[1]) / seg_len)


def sample_polyline_points(points: Sequence[Point], n: int = 12) -> list[Point]:
    """沿折线等弧长采样 ``n + 1`` 个点（供碰撞避让计算使用）。"""
    total = polyline_length(points)
    if total < 1e-6 or n <= 0:
        return list(points)
    return [sample_polyline_at(points, total * i / n)[0] for i in range(n + 1)]


def trim_polyline(points: Sequence[Point], trim0: float = 0.0, trim1: float = 0.0) -> list[Point]:
    """从折线首尾分别裁剪指定弧长距离，避免轨道侵入圆环站台内部。"""
    if len(points) < 2:
        return list(points)
    cum = polyline_arclengths(points)
    total = cum[-1]
    if total < 1e-6:
        return list(points)
    t0 = max(0.0, float(trim0))
    t1 = max(0.0, float(trim1))
    if t0 + t1 >= total * 0.90:
        scale = (total * 0.85) / max(t0 + t1, 1e-6)
        t0 *= scale
        t1 *= scale
    s_start = t0
    s_end = total - t1
    if s_end <= s_start + 1e-5:
        mid, _ = sample_polyline_at(points, total * 0.5)
        return [mid, mid]

    p_start, _ = sample_polyline_at(points, s_start)
    p_end, _ = sample_polyline_at(points, s_end)
    out: list[Point] = [p_start]
    for pt, s in zip(points, cum):
        if s_start + 1e-5 < s < s_end - 1e-5:
            out.append(pt)
    out.append(p_end)
    return out


# ---------------------------------------------------------------------------
# 地铁轨道、内嵌箭头、跑道闭环、换乘站与区室底板绘制原语
# ---------------------------------------------------------------------------

def _draw_inline_chevron(
    ax,
    point: Point,
    tangent: Point,
    *,
    half_width: float,
    tail_depth: float,
    color: str,
    lw: float,
    halo_color: str | None = None,
    halo_lw: float = 0.0,
    zorder: float = 3.2,
) -> None:
    """在轨道内部绘制锐利微型方向箭头 ``›``。"""
    px, py = point
    tx, ty = tangent
    norm = math.hypot(tx, ty) or 1.0
    tx, ty = tx / norm, ty / norm
    nx, ny = -ty, tx

    tip = (px + tx * tail_depth * 0.52, py + ty * tail_depth * 0.52)
    w1 = (px - tx * tail_depth * 0.48 + nx * half_width, py - ty * tail_depth * 0.48 + ny * half_width)
    w2 = (px - tx * tail_depth * 0.48 - nx * half_width, py - ty * tail_depth * 0.48 - ny * half_width)
    path = MplPath([w1, tip, w2], [MplPath.MOVETO, MplPath.LINETO, MplPath.LINETO])

    if halo_color and halo_lw > 0:
        ax.add_patch(PathPatch(
            path, facecolor="none", edgecolor=halo_color, lw=halo_lw,
            capstyle="round", joinstyle="miter", zorder=zorder - 0.02,
        ))
    ax.add_patch(PathPatch(
        path, facecolor="none", edgecolor=color, lw=lw,
        capstyle="round", joinstyle="miter", zorder=zorder,
    ))


def draw_metro_coupled_halo(
    ax,
    polylines: Sequence[Sequence[Point]],
    *,
    color: str = "#DBEAFE",
    edge_color: str = "#4F46E5",
    width: float = 8.4,
    zorder: float = 1.85,
) -> None:
    """在地铁轨道下方绘制半透明靛青耦合通量底轨（串联化学计量锁死的反应组）。"""
    for pts in polylines:
        if len(pts) < 2:
            continue
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        ax.add_line(Line2D(
            xs, ys, color=edge_color, lw=width + 1.0, alpha=0.22,
            solid_capstyle="round", solid_joinstyle="round", zorder=zorder,
        ))
        ax.add_line(Line2D(
            xs, ys, color=color, lw=width, alpha=0.88,
            solid_capstyle="round", solid_joinstyle="round", zorder=zorder + 0.05,
        ))


def draw_metro_track(
    ax,
    points: Sequence[Point],
    color,
    *,
    width: float = 3.4,
    hollow: bool = False,
    double_gap_ratio: float = 0.46,
    surface: str = "#FFFFFF",
    chevrons: bool = True,
    chevron_spacing: float = 5.5,
    chevron_direction: int = 1,
    chevron_bidirectional: bool = False,
    unit_in: float = 0.08,
    trim0: float = 0.0,
    trim1: float = 0.0,
    corner_radius: float = 1.25,
    closed: bool = False,
    halo: bool = True,
    zorder: float = 2.6,
    alpha: float = 1.0,
) -> list[Point]:
    """绘制单条 Beck 地铁轨道（支持实心主轨 / 空心白芯双轨 + 等距内嵌方向 chevron）。

    返回实际裁剪平滑后的折线点列。
    """
    if len(points) < 2:
        return list(points)
    smooth = fillet_polyline(points, radius=corner_radius, closed=closed) if corner_radius > 0 else list(points)
    drawn = smooth if closed else trim_polyline(smooth, trim0=trim0, trim1=trim1)
    if len(drawn) < 2:
        return drawn

    xs = [p[0] for p in drawn]
    ys = [p[1] for p in drawn]
    cap = "round" if closed else "butt"

    if halo:
        ax.add_line(Line2D(
            xs, ys, color=surface, lw=width + 1.5, alpha=0.92 * alpha,
            solid_capstyle=cap, solid_joinstyle="round", zorder=zorder - 0.14,
        ))

    ax.add_line(Line2D(
        xs, ys, color=color, lw=width, alpha=alpha,
        solid_capstyle=cap, solid_joinstyle="round", zorder=zorder,
    ))

    if hollow:
        inner_w = max(0.75, width * double_gap_ratio)
        ax.add_line(Line2D(
            xs, ys, color=surface, lw=inner_w, alpha=alpha,
            solid_capstyle=cap, solid_joinstyle="round", zorder=zorder + 0.05,
        ))

    if chevrons:
        total = polyline_length(drawn)
        # 根据物理线宽精确换算数据坐标下的轨道半宽，保证跨幅面（8 in ~ 22 in）内嵌箭头始终严格贴合轨道
        track_w_data = (width / 72.0) / max(float(unit_in), 1e-4)
        half_w = track_w_data * (0.44 if hollow else 0.36)
        tail_d = half_w * 1.35
        chev_lw = max(0.65, width * 0.22)
        chev_col = color if hollow else surface
        chev_halo = surface if hollow else None
        chev_halo_lw = max(1.15, width * 0.42) if hollow else 0.0

        if closed:
            n_chev = max(4, int(round(total / max(chevron_spacing, 1.0))))
            positions = [(i + 0.5) * (total / n_chev) for i in range(n_chev)]
        else:
            margin = max(track_w_data * 1.45, 0.95)
            usable = total - 2.0 * margin
            if usable > 0.35:
                n_chev = max(1, int(round(usable / max(chevron_spacing, 1.0))) + 1)
                if chevron_bidirectional and n_chev == 1 and total > margin * 2.4:
                    n_chev = 2
                if n_chev == 1:
                    positions = [total * 0.5]
                else:
                    positions = [margin + usable * i / (n_chev - 1) for i in range(n_chev)]
            elif total > track_w_data * 1.3:
                positions = [total * 0.5]
            else:
                positions = []

        for idx, s in enumerate(positions):
            pt, tan = sample_polyline_at(drawn, s)
            if chevron_bidirectional:
                sign = -1 if idx < len(positions) / 2 else 1
            else:
                sign = -1 if chevron_direction < 0 else 1
            _draw_inline_chevron(
                ax, pt, (tan[0] * sign, tan[1] * sign),
                half_width=half_w,
                tail_depth=tail_d,
                color=chev_col,
                lw=chev_lw,
                halo_color=chev_halo,
                halo_lw=chev_halo_lw,
                zorder=zorder + 0.12,
            )

    return drawn


def stadium_loop_points(
    center: Point,
    width: float,
    height: float,
    *,
    corner_radius: float | None = None,
    arc_steps: int = 12,
) -> list[Point]:
    """生成顺时针（y 向下画布）圆角跑道闭环点列。"""
    cx, cy = center
    hw, hh = width * 0.5, height * 0.5
    r = min(corner_radius if corner_radius is not None else hh * 0.82, hw * 0.65, hh * 0.95)
    x0, y0, x1, y1 = cx - hw, cy - hh, cx + hw, cy + hh

    # 四个圆角中心与角度区间（y 轴向下，顺时针：顶→右→底→左）
    corners = [
        (x1 - r, y0 + r, -math.pi / 2, 0.0),
        (x1 - r, y1 - r, 0.0, math.pi / 2),
        (x0 + r, y1 - r, math.pi / 2, math.pi),
        (x0 + r, y0 + r, math.pi, 3.0 * math.pi / 2),
    ]
    pts: list[Point] = []
    for ccx, ccy, a0, a1 in corners:
        for k in range(arc_steps + 1):
            ang = a0 + (a1 - a0) * (k / arc_steps)
            pts.append((ccx + r * math.cos(ang), ccy + r * math.sin(ang)))
    pts.append(pts[0])
    return pts


def draw_metro_loop(
    ax,
    center: Point,
    width: float,
    height: float,
    color,
    *,
    track_width: float = 4.2,
    hollow: bool = False,
    surface: str = "#FFFFFF",
    interior_tint: float = 0.045,
    clockwise: bool = True,
    n_chevrons: int = 6,
    label: str = "",
    subtitle: str = "",
    label_color: str | None = None,
    subtitle_color: str | None = None,
    label_fontsize: float = 7.4,
    subtitle_fontsize: float = 5.8,
    unit_in: float = 0.062,
    zorder: float = 3.2,
) -> list[Point]:
    """绘制 Metabolic Metro 标志性圆角跑道闭环（适用于 TCA / PPP / β-Oxidation 等循环通路）。"""
    pts = stadium_loop_points(center, width, height)
    if not clockwise:
        pts = list(reversed(pts))

    # 极浅内部底晕，既保持循环通透感又让标题易读
    inner_fill = tint_color(color, interior_tint, base=surface) if interior_tint > 0 else surface
    ax.add_patch(Polygon(pts[:-1], closed=True, facecolor=inner_fill, edgecolor="none", zorder=zorder - 0.35))

    total = polyline_length(pts)
    spacing = total / max(int(n_chevrons), 2)
    draw_metro_track(
        ax, pts, color,
        width=track_width,
        hollow=hollow,
        surface=surface,
        chevrons=True,
        chevron_spacing=spacing,
        chevron_direction=1,
        unit_in=unit_in,
        corner_radius=0.0,
        closed=True,
        halo=True,
        zorder=zorder,
    )

    cx, cy = center
    if label:
        lcol = label_color or color
        if subtitle:
            ax.text(cx, cy - height * 0.14, label, fontsize=label_fontsize, fontweight="bold",
                    style="italic", color=lcol, ha="center", va="center", zorder=zorder + 0.5)
            ax.text(cx, cy + height * 0.21, subtitle, fontsize=subtitle_fontsize,
                    color=subtitle_color or lcol, fontweight="medium",
                    ha="center", va="center", zorder=zorder + 0.5)
        else:
            ax.text(cx, cy, label, fontsize=label_fontsize, fontweight="bold",
                    style="italic", color=lcol, ha="center", va="center", zorder=zorder + 0.5)
    return pts


def draw_metro_station(
    ax,
    xy: Point,
    radius: float,
    *,
    fill="#FFFFFF",
    edge="#1E293B",
    outer_ring: str | None = None,
    lw: float = 1.45,
    surface: str = "#FFFFFF",
    halo_pad: float = 0.30,
    slash_color: str | None = None,
    glyph: str = "",
    glyph_color: str = "#64748B",
    glyph_fontsize: float = 4.2,
    zorder: float = 4.6,
) -> None:
    """绘制地铁双层圆环换乘站台（白外晕 + 粗色环 + 内部语义填充）。"""
    x, y = xy
    if halo_pad > 0:
        ax.add_patch(Circle((x, y), radius + halo_pad, facecolor=surface, edgecolor="none", zorder=zorder - 0.12))
    if outer_ring:
        ax.add_patch(Circle(
            (x, y), radius + max(0.18, halo_pad * 0.55),
            facecolor="none", edgecolor=outer_ring, lw=max(1.1, lw * 0.85), zorder=zorder + 0.10,
        ))
    ax.add_patch(Circle((x, y), radius, facecolor=fill, edgecolor=edge, lw=lw, zorder=zorder))

    if glyph:
        ax.text(x, y, glyph, fontsize=glyph_fontsize, ha="center", va="center", color=glyph_color, zorder=zorder + 0.15)
    elif slash_color:
        d = radius * 0.44
        ax.add_line(Line2D([x - d, x + d], [y + d, y - d], color=slash_color, lw=0.60, zorder=zorder + 0.15))


def draw_metro_zone(
    ax,
    box: tuple[float, float, float, float],
    label: str,
    *,
    face: str = "#EFF1F5",
    edge: str = "none",
    lw: float = 0.0,
    inner_edge: str | None = None,
    label_color: str = "#64748B",
    label_pos: str = "bottom_left",
    rounding: float = 2.8,
    fontsize: float = 6.2,
    zorder: float = 0.4,
) -> None:
    """绘制低存在感浅灰圆角区室底板（将 100% 视觉权重让给地铁轨道与换乘站）。"""
    x0, y0, x1, y1 = box
    ax.add_patch(FancyBboxPatch(
        (x0, y0), x1 - x0, y1 - y0,
        boxstyle=f"round,pad=0,rounding_size={rounding}",
        facecolor=face, edgecolor=edge, lw=lw, zorder=zorder,
    ))
    if inner_edge:
        pad = min(0.85, (x1 - x0) * 0.04, (y1 - y0) * 0.04)
        ax.add_patch(FancyBboxPatch(
            (x0 + pad, y0 + pad), (x1 - x0) - 2 * pad, (y1 - y0) - 2 * pad,
            boxstyle=f"round,pad=0,rounding_size={max(0.6, rounding - 0.5)}",
            facecolor="none", edgecolor=inner_edge, lw=0.60, linestyle=(0, (3.0, 2.2)), zorder=zorder + 0.06,
        ))
    if not label:
        return
    if label_pos == "top_center":
        lx, ly, ha, va = (x0 + x1) * 0.5, y0 + 2.0, "center", "center"
    elif label_pos == "bottom_center":
        lx, ly, ha, va = (x0 + x1) * 0.5, y1 - 2.0, "center", "center"
    elif label_pos == "top_right":
        lx, ly, ha, va = x1 - 1.4, y0 + 1.35, "right", "top"
    elif label_pos == "top_left":
        lx, ly, ha, va = x0 + 1.4, y0 + 1.35, "left", "top"
    elif label_pos == "bottom_right":
        lx, ly, ha, va = x1 - 1.4, y1 - 0.85, "right", "bottom"
    else:  # bottom_left
        lx, ly, ha, va = x0 + 1.4, y1 - 1.25, "left", "bottom"
    ax.text(
        lx, ly, label, fontsize=fontsize, fontweight="bold", style="italic",
        color=label_color, ha=ha, va=va, zorder=zorder + 0.25,
    )


__all__ = [
    "Point",
    "draw_metro_coupled_halo",
    "draw_metro_loop",
    "draw_metro_station",
    "draw_metro_track",
    "draw_metro_zone",
    "fillet_polyline",
    "polyline_arclengths",
    "polyline_length",
    "route_octilinear",
    "sample_polyline_at",
    "sample_polyline_points",
    "stadium_loop_points",
    "trim_polyline",
]
