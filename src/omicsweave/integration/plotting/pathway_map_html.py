#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
方案 A：单页离线可交互通路探索器（原生语义化 SVG + 内联紧凑 JSON + 零外部依赖原生 JS）。

设计要点：
- 与静态 PDF/PNG 共享 :func:`pathway_map_plots.build_pathway_map_scene` 的全部正交路由、
  空心柔色双轨、避让卡片坐标与锁死通量分组，并按 ``em_units(pt, unit_in)`` 1:1 复刻
  ``metro_atlas`` PDF 的精细边框线宽、字号比例、聚卡片 P/R/M 色块与徽章几何；
- 默认采用学术出版级 ``metro_atlas`` 视觉风格，内置顶部实时主题切换
  （``Metro Atlas`` / ``Economist`` / ``Nature Editorial`` / ``Slide Dark``）；
- 内置浏览器端 ``getBBox()`` 零溢出自校准守护（适配 Windows/macOS/Linux 字体差异），
  确保卡片、徽章、聚合卡与边界胶囊中的任何标签都不会越界；
- 右侧 Reaction Dossier 抽屉内置 **逐只动物/样本散点小图**（按 ``dossier_samples.parquet`` 与
  ``dossier_sample_contrasts.parquet`` 绘制，含 MFA 拟合四分位误差棒、``at_bound`` 空心菱形、
  0 基线与规范横网格/右对齐刻度）及全量成员统计明细表。
"""
from __future__ import annotations

import html
import json
import math
from pathlib import Path
from typing import Any

import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

from omicsweave.visualization.metro import (
    fillet_polyline,
    polyline_length,
    sample_polyline_at,
    trim_polyline,
)
from omicsweave.visualization.palettes import (
    DEFAULT_NORM,
    DEFAULT_VMAX,
    ECONOMIST_RED,
    PALETTES,
    SemanticPalette,
    get_palette,
    pathway_family_color,
    tier_tokens,
    tint_color,
)
from omicsweave.mfa.tracers import format_tracer_badge_label, format_tracer_display
from omicsweave.visualization.primitives import (
    AXIS_LAYER_ICONS,
    OMICS_LAYER_ICONS,
    coupled_badge_label,
    coupled_badge_width,
    em_units,
    tier_badge_text,
    tier_badge_width,
)

from ..mappings.network_annotation import pathway_reactions
from ..mappings.network_layout import Layout
from .pathway_map_plots import (
    AXIS_BAR_H,
    AXIS_HDR_H,
    CELL,
    MIN_EDGE,
    NODE_R,
    NOTE_SIZE,
    UNIT,
    _agg_card_dims,
    _edge_color,
    _enzyme_inferred_direction,
    _extract_strip_states,
    _format_fusion_items,
    _independent_flux_chains,
    _is_compact_label_only,
    _short,
    _species_short,
    _strip_badges,
    _strip_size,
    _text_width,
    build_pathway_map_scene,
)
from .pathway_submaps import SUBMAP_CATALOG


# ---------------------------------------------------------------------------
# Color, Unit & Numeric Utilities
# ---------------------------------------------------------------------------

def _pt(pt_val: float, unit_in: float = UNIT) -> float:
    """将 Matplotlib 点值（pt）精确换算为当前画布数据坐标单位（1:1 对齐 PDF 线宽与字号）。"""
    return em_units(pt_val, unit_in)


def _to_hex(color: Any) -> str:
    if isinstance(color, str) and color.startswith("#") and len(color) in (4, 7):
        return color
    return mcolors.to_hex(color)


def _num(val: Any, digits: int = 4) -> float | None:
    if val is None:
        return None
    try:
        f = float(val)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(f):
        return None
    return round(f, digits)


def _pts_path(pts: list[tuple[float, float]]) -> str:
    if not pts:
        return ""
    parts = [f"M{pts[0][0]:.2f},{pts[0][1]:.2f}"]
    parts.extend(f"L{x:.2f},{y:.2f}" for x, y in pts[1:])
    return "".join(parts)


def _esc(text: Any) -> str:
    return html.escape(str(text or ""), quote=True)


def _palette_dict(name: str) -> dict[str, Any]:
    pal = get_palette(name)
    tok = tier_tokens(name)
    is_metro = name == "metro_atlas"
    is_econ = name == "economist"
    is_dark = name == "slide_dark"
    return {
        "name": pal.name,
        "surface": pal.surface,
        "card_bg": pal.card_bg,
        "strip_card_bg": pal.card_bg if is_econ else pal.surface,
        "card_edge": pal.card_edge,
        "ink": pal.ink,
        "ink2": pal.ink2,
        "muted": pal.muted,
        "hairline": pal.hairline,
        "up": pal.up,
        "down": pal.down,
        "neutral": pal.neutral,
        "diverging_stops": list(pal.diverging_stops),
        "numerator": pal.numerator,
        "denominator": pal.denominator,
        "edge_grey": pal.edge_grey,
        "omics_grey": pal.omics_grey,
        "faint": pal.faint,
        "missing_fill": pal.missing_fill,
        "missing_hatch": pal.missing_hatch,
        "coupled_track": pal.coupled_track,
        "coupled_fill": pal.coupled_fill,
        "coupled_edge": pal.coupled_edge,
        "coupled_ink": pal.coupled_ink,
        "sign_change_bg": pal.sign_change_bg,
        "sign_change_edge": pal.sign_change_edge,
        "sign_change_ink": pal.sign_change_ink,
        "mito_face": "#ECEEF2" if is_metro else pal.mito_face,
        "mito_edge": "#CFD6E0" if is_metro else pal.mito_edge,
        "mito_inner": "#C5CDD8" if is_metro else pal.mito_inner,
        "mito_ink": "#64748B" if is_metro else pal.mito_ink,
        "perox_face": "#EBF4F2" if is_metro else pal.perox_face,
        "perox_edge": "#B4D7CF" if is_metro else pal.perox_edge,
        "perox_inner": "#9CCDC2" if is_metro else pal.perox_inner,
        "perox_ink": "#336D65" if is_metro else pal.perox_ink,
        "boundary_face": pal.boundary_face,
        "boundary_edge": pal.boundary_edge,
        "boundary_ink": pal.boundary_ink,
        "store_face": pal.store_face,
        "store_edge": pal.store_edge,
        "store_ink": pal.store_ink,
        "region_face": "#E5E4D7" if is_econ else ("#F3F5F8" if is_metro else ("#182230" if is_dark else "#F7F8FA")),
        "region_edge": pal.card_edge if is_econ else ("#E2E7EF" if is_metro else pal.card_edge),
        "oxphos_face": "#DDDCD0" if is_econ else ("#E6E9F0" if is_metro else ("#1F2B3D" if is_dark else "#EFF2F7")),
        "oxphos_edge": pal.card_edge if is_econ else ("#CBD2DD" if is_metro else pal.card_edge),
        "red_tag": ECONOMIST_RED if is_econ else pal.up,
        "tier": {
            lvl: {
                "bg": str(tok[lvl]["bg"]),
                "ec": str(tok[lvl]["ec"]),
                "fg": str(tok[lvl]["fg"]),
                "dot_on": str(tok[lvl]["dot_on"]),
                "dot_off": str(tok[lvl]["dot_off"]),
            }
            for lvl in ("L3", "L2", "L1")
        },
    }


# ---------------------------------------------------------------------------
# Pure-SVG Micro-Icons & Status Cells (exact pt-calibrated parity with primitives.py)
# ---------------------------------------------------------------------------

def _svg_layer_icon(kind: str, cx: float, cy: float, size: float, color: str, unit_in: float = UNIT) -> str:
    r = size * 0.44
    sw = _pt(0.56, unit_in)
    if kind == "flux":
        yt, yb = cy - r * 0.46, cy + r * 0.46
        x0, x1 = cx - r * 0.90, cx + r * 0.88
        ah = r * 0.36
        d = (
            f"M{x0:.2f},{yt:.2f}L{x1:.2f},{yt:.2f}"
            f"M{x1 - ah:.2f},{yt - ah * 0.72:.2f}L{x1:.2f},{yt:.2f}L{x1 - ah:.2f},{yt + ah * 0.72:.2f}"
            f"M{x0:.2f},{yb:.2f}L{x1:.2f},{yb:.2f}"
            f"M{x1 - ah:.2f},{yb - ah * 0.72:.2f}L{x1:.2f},{yb:.2f}L{x1 - ah:.2f},{yb + ah * 0.72:.2f}"
        )
        return f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{sw:.3f}" stroke-linecap="round" stroke-linejoin="miter"/>'
    if kind == "protein":
        a1, a2 = math.radians(34), math.radians(326)
        rad = r * 0.92
        x1, y1 = cx + rad * math.cos(a1), cy + rad * math.sin(a1)
        x2, y2 = cx + rad * math.cos(a2), cy + rad * math.sin(a2)
        d = f"M{cx:.2f},{cy:.2f}L{x1:.2f},{y1:.2f}A{rad:.2f},{rad:.2f} 0 1,1 {x2:.2f},{y2:.2f}Z"
        return f'<path d="{d}" fill="{color}"/>'
    if kind == "mrna":
        pts = []
        for i in range(9):
            t = i / 8.0
            px = cx - r * 0.92 + 1.84 * r * t
            py = cy - r * 0.50 * math.sin(2.0 * math.pi * t)
            pts.append(f"{px:.2f},{py:.2f}")
        d = "M" + "L".join(pts)
        p0x, p0y = pts[0].split(",")
        return (
            f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{sw * 1.08:.3f}" stroke-linecap="round"/>'
            f'<circle cx="{p0x}" cy="{p0y}" r="{r * 0.20:.2f}" fill="{color}"/>'
        )
    if kind == "metabolite":
        pts = []
        for i in range(6):
            ang = math.radians(60 * i)
            pts.append(f"{cx + r * 0.88 * math.cos(ang):.2f},{cy + r * 0.88 * math.sin(ang):.2f}")
        return f'<polygon points="{" ".join(pts)}" fill="none" stroke="{color}" stroke-width="{sw * 1.05:.3f}" stroke-linejoin="miter"/>'
    pts = f"{cx:.2f},{cy - r:.2f} {cx + r:.2f},{cy:.2f} {cx:.2f},{cy + r:.2f} {cx - r:.2f},{cy:.2f}"
    return f'<polygon points="{pts}" fill="none" stroke="{color}" stroke-width="{sw:.3f}"/>'


def _svg_status_cell(
    x: float, y: float, cell: float, state: str, direction: int, pal: SemanticPalette,
    *, unit_in: float = UNIT, draw_bg: bool = True, card_base: str | None = None,
) -> str:
    cx, cy = x + cell / 2.0, y + cell / 2.0
    base_col = card_base or pal.surface
    sw_tri = _pt(0.62, unit_in)
    sw_sym = _pt(0.58, unit_in)
    parts: list[str] = []
    if state in ("changed", "suggestive") and direction != 0:
        is_up = direction > 0
        col = pal.up if is_up else pal.down
        fill = col if state == "changed" else "none"
        role = ("cell-up" if is_up else "cell-down") + ("-solid" if state == "changed" else "-hollow")
        w, h = cell * 0.27, cell * 0.23
        if is_up:
            pts = f"{cx:.2f},{cy - h:.2f} {cx - w:.2f},{cy + h:.2f} {cx + w:.2f},{cy + h:.2f}"
        else:
            pts = f"{cx:.2f},{cy + h:.2f} {cx - w:.2f},{cy - h:.2f} {cx + w:.2f},{cy - h:.2f}"
        if draw_bg:
            bg = _to_hex(tint_color(col, 0.14 if state == "changed" else 0.06, base=base_col))
            parts.append(f'<rect class="ow-cell-bg" data-crole="{role}" x="{x + 0.02:.2f}" y="{y + 0.02:.2f}" width="{cell - 0.04:.2f}" height="{cell - 0.04:.2f}" fill="{bg}"/>')
        parts.append(f'<polygon class="ow-cell-sym" data-crole="{role}" points="{pts}" fill="{fill}" stroke="{col}" stroke-width="{sw_tri:.3f}"/>')
    elif state == "mixed":
        d = cell * 0.22
        parts.append(
            f'<path d="M{cx - d:.2f},{cy - d * 0.45:.2f}L{cx + d:.2f},{cy - d * 0.45:.2f}'
            f'M{cx - d:.2f},{cy + d * 0.45:.2f}L{cx + d:.2f},{cy + d * 0.45:.2f}'
            f'M{cx - d * 0.55:.2f},{cy + d * 0.95:.2f}L{cx + d * 0.55:.2f},{cy - d * 0.95:.2f}" '
            f'fill="none" stroke="{pal.ink2}" stroke-width="{sw_sym:.3f}" stroke-linecap="round"/>'
        )
    elif state in ("unchanged", "equivalent"):
        d = cell * 0.22
        parts.append(f'<line x1="{cx - d:.2f}" y1="{cy:.2f}" x2="{cx + d:.2f}" y2="{cy:.2f}" stroke="{pal.ink2}" stroke-width="{sw_sym:.3f}" stroke-linecap="round"/>')
    elif state in ("undetermined", "no_evidence"):
        half = cell * 0.26
        gap = cell * 0.06
        sw_dash = _pt(0.54, unit_in)
        parts.append(f'<line x1="{cx - half:.2f}" y1="{cy:.2f}" x2="{cx - gap:.2f}" y2="{cy:.2f}" stroke="{pal.muted}" stroke-width="{sw_dash:.3f}" stroke-linecap="round"/>')
        parts.append(f'<line x1="{cx + gap:.2f}" y1="{cy:.2f}" x2="{cx + half:.2f}" y2="{cy:.2f}" stroke="{pal.muted}" stroke-width="{sw_dash:.3f}" stroke-linecap="round"/>')
    elif state in ("not_identifiable", "bound_pinned"):
        r = cell * 0.22
        if draw_bg:
            parts.append(f'<rect x="{x + 0.02:.2f}" y="{y + 0.02:.2f}" width="{cell - 0.04:.2f}" height="{cell - 0.04:.2f}" fill="{pal.missing_fill}"/>')
        parts.append(f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{r:.2f}" fill="none" stroke="{pal.muted}" stroke-width="{sw_sym:.3f}"/>')
        dd = r * 0.72
        parts.append(f'<line x1="{cx - dd:.2f}" y1="{cy + dd:.2f}" x2="{cx + dd:.2f}" y2="{cy - dd:.2f}" stroke="{pal.muted}" stroke-width="{sw_sym:.3f}"/>')
    elif state == "not_detected":
        r = cell * 0.22
        if draw_bg:
            parts.append(f'<rect x="{x + 0.02:.2f}" y="{y + 0.02:.2f}" width="{cell - 0.04:.2f}" height="{cell - 0.04:.2f}" fill="{pal.missing_fill}"/>')
        parts.append(
            f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{r:.2f}" fill="none" stroke="{pal.muted}" '
            f'stroke-width="{sw_sym:.3f}" stroke-dasharray="{_pt(1.2, unit_in):.3f} {_pt(1.0, unit_in):.3f}"/>'
        )
    elif state in ("insufficient_n", "qc_blocked"):
        d = cell * 0.21
        off = cell * 0.09
        sw_dbl = _pt(0.50, unit_in)
        if draw_bg:
            parts.append(f'<rect x="{x + 0.02:.2f}" y="{y + 0.02:.2f}" width="{cell - 0.04:.2f}" height="{cell - 0.04:.2f}" fill="{pal.missing_fill}"/>')
        for dx_off in (-off, off):
            parts.append(
                f'<line x1="{cx - d + dx_off:.2f}" y1="{cy + d:.2f}" x2="{cx + d + dx_off:.2f}" y2="{cy - d:.2f}" '
                f'stroke="{pal.muted}" stroke-width="{sw_dbl:.3f}" stroke-linecap="round"/>'
            )
    else:  # not_measured / no_knowledge / not_infused
        d = cell * 0.19
        if draw_bg:
            parts.append(f'<rect x="{x + 0.02:.2f}" y="{y + 0.02:.2f}" width="{cell - 0.04:.2f}" height="{cell - 0.04:.2f}" fill="{pal.missing_fill}"/>')
        parts.append(f'<line x1="{cx - d:.2f}" y1="{cy + d:.2f}" x2="{cx + d:.2f}" y2="{cy - d:.2f}" stroke="{pal.missing_hatch}" stroke-width="{_pt(0.52, unit_in):.3f}"/>')
    return "".join(parts)


def _svg_tier_badge(
    x: float, y: float, level: str, scale: float, pal: SemanticPalette, preset: str, unit_in: float,
    *, label: str | None = None,
) -> tuple[str, float]:
    if level not in ("L1", "L2", "L3"):
        return "", 0.0
    tok = tier_tokens(preset)[level]
    text = tier_badge_text(level, label=label)
    bw = tier_badge_width(level, scale=scale, unit_in=unit_in, label=label)
    bh = 1.14 * scale
    # Exact PDF pt -> data-unit conversion (0.90 factor absorbs web sans-serif bold metrics)
    fs = _pt(3.40 * scale * 0.90, unit_in)
    rx = (0.12 if preset == "economist" else 0.38) * scale
    sw = _pt(float(tok["lw"]) * scale, unit_in)
    dash = f' stroke-dasharray="{_pt(2.2, unit_in):.3f} {_pt(1.2, unit_in):.3f}"' if level == "L1" else ""
    parts = [
        f'<g class="ow-tier-badge" data-tier="{level}" data-scale="{scale:.3f}">',
        f'<rect class="ow-tb-bg" x="{x:.2f}" y="{y - bh / 2:.2f}" width="{bw:.2f}" height="{bh:.2f}" rx="{rx:.2f}" '
        f'fill="{tok["bg"]}" stroke="{tok["ec"]}" stroke-width="{sw:.3f}"{dash}/>',
    ]
    n_on = int(tok["dots"])
    dot_r = 0.155 * scale
    dot_step = 0.395 * scale
    dot_x0 = x + 0.41 * scale
    for idx in range(3):
        cx = dot_x0 + idx * dot_step
        if idx < n_on:
            parts.append(f'<circle class="ow-tb-dot-on" cx="{cx:.2f}" cy="{y:.2f}" r="{dot_r:.2f}" fill="{tok["dot_on"]}"/>')
        else:
            parts.append(
                f'<circle class="ow-tb-dot-off" cx="{cx:.2f}" cy="{y:.2f}" r="{dot_r * 0.92:.2f}" '
                f'fill="none" stroke="{tok["dot_off"]}" stroke-width="{_pt(0.45 * scale, unit_in):.3f}"/>'
            )
    tx = dot_x0 + 2 * dot_step + dot_r + 0.22 * scale
    parts.append(
        f'<text class="ow-tb-txt" x="{tx:.2f}" y="{y + 0.02 * scale:.2f}" font-size="{fs:.3f}" '
        f'font-weight="700" fill="{tok["fg"]}" dominant-baseline="central">{_esc(text)}</text></g>'
    )
    return "".join(parts), bw


def _svg_coupled_badge(
    x: float, y: float, group_num: int, total_groups: int, scale: float, pal: SemanticPalette, preset: str, unit_in: float,
) -> tuple[str, float]:
    txt = coupled_badge_label(group_num, total_groups=total_groups)
    bw = coupled_badge_width(group_num, scale=scale, total_groups=total_groups, unit_in=unit_in)
    bh = 1.14 * scale
    fs = _pt(3.30 * scale * 0.88, unit_in)
    rx = (0.12 if preset == "economist" else 0.36) * scale
    sw_box = _pt(0.48 * scale, unit_in)
    sw_link = _pt(0.62 * scale, unit_in)
    lh = 0.40 * scale
    lw_link = 0.62 * scale
    lx0 = x + 0.27 * scale
    svg = (
        f'<g class="ow-coupled-badge" data-scale="{scale:.3f}">'
        f'<rect class="ow-cb-bg" x="{x:.2f}" y="{y - bh / 2:.2f}" width="{bw:.2f}" height="{bh:.2f}" rx="{rx:.2f}" '
        f'fill="{pal.coupled_fill}" stroke="{pal.coupled_edge}" stroke-width="{sw_box:.3f}"/>'
        f'<rect x="{lx0:.2f}" y="{y - lh / 2:.2f}" width="{lw_link:.2f}" height="{lh:.2f}" rx="{lh / 2:.2f}" '
        f'fill="none" stroke="{pal.coupled_edge}" stroke-width="{sw_link:.3f}"/>'
        f'<rect x="{lx0 + 0.35 * scale:.2f}" y="{y - lh / 2:.2f}" width="{lw_link:.2f}" height="{lh:.2f}" rx="{lh / 2:.2f}" '
        f'fill="none" stroke="{pal.coupled_ink}" stroke-width="{sw_link:.3f}"/>'
        f'<text class="ow-cb-txt" x="{x + 1.42 * scale:.2f}" y="{y + 0.02 * scale:.2f}" font-size="{fs:.3f}" '
        f'font-weight="700" fill="{pal.coupled_ink}" dominant-baseline="central">{_esc(txt)}</text></g>'
    )
    return svg, bw


def _svg_sign_reversal_badge(cx: float, cy: float, radius: float, pal: SemanticPalette, unit_in: float) -> str:
    """纯矢量琥珀圆反转徽章（与 primitives.draw_sign_reversal_badge 完全一致，无字体偏移风险）。"""
    r = float(radius)
    sw_ring = _pt(0.68, unit_in)
    sw_arr = _pt(0.65, unit_in)
    yt = cy - r * 0.27
    yb = cy + r * 0.27
    d = (
        f"M{cx - r * 0.54:.2f},{yt:.2f}L{cx + r * 0.44:.2f},{yt:.2f}"
        f"M{cx + r * 0.16:.2f},{yt - r * 0.24:.2f}L{cx + r * 0.56:.2f},{yt:.2f}L{cx + r * 0.16:.2f},{yt + r * 0.24:.2f}"
        f"M{cx + r * 0.54:.2f},{yb:.2f}L{cx - r * 0.44:.2f},{yb:.2f}"
        f"M{cx - r * 0.16:.2f},{yb - r * 0.24:.2f}L{cx - r * 0.56:.2f},{yb:.2f}L{cx - r * 0.16:.2f},{yb + r * 0.24:.2f}"
    )
    return (
        f'<g class="ow-sign-rev">'
        f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{r:.2f}" fill="{pal.sign_change_bg}" '
        f'stroke="{pal.sign_change_edge}" stroke-width="{sw_ring:.3f}"/>'
        f'<path d="{d}" fill="none" stroke="{pal.sign_change_ink}" stroke-width="{sw_arr:.3f}" '
        f'stroke-linecap="round" stroke-linejoin="miter"/>'
        f"</g>"
    )


def _svg_chevrons(
    pts: list[tuple[float, float]], direction: int, bidirectional: bool, stroke: str, sw: float, unit_in: float,
) -> str:
    total = polyline_length(pts)
    if total < 1.2:
        return ""
    if bidirectional:
        dists = [(total * 0.28, -1), (total * 0.72, 1)] if total >= 3.2 else [(total * 0.5, 1)]
    else:
        spacing = 4.8
        n = max(1, min(4, int(round(total / spacing))))
        dists = [(total * (i + 1) / (n + 1), direction) for i in range(n)]
    chev_len = max(0.32, min(0.50, sw * 1.15))
    chev_half = chev_len * 0.66
    chev_sw = max(_pt(0.70, unit_in), sw * 0.30)
    segs: list[str] = []
    for dist, dsign in dists:
        (px, py), (tx, ty) = sample_polyline_at(pts, dist)
        ux, uy = tx * dsign, ty * dsign
        nx, ny = -uy, ux
        tip = (px + ux * chev_len * 0.45, py + uy * chev_len * 0.45)
        w1 = (px - ux * chev_len * 0.55 + nx * chev_half, py - uy * chev_len * 0.55 + ny * chev_half)
        w2 = (px - ux * chev_len * 0.55 - nx * chev_half, py - uy * chev_len * 0.55 - ny * chev_half)
        segs.append(f"M{w1[0]:.2f},{w1[1]:.2f}L{tip[0]:.2f},{tip[1]:.2f}L{w2[0]:.2f},{w2[1]:.2f}")
    return (
        f'<path class="ow-edge-chev" d="{"".join(segs)}" fill="none" stroke="{stroke}" '
        f'stroke-width="{chev_sw:.3f}" stroke-linecap="round" stroke-linejoin="miter"/>'
    )


# ---------------------------------------------------------------------------
# Semantic SVG Scene Serializer (1:1 visual finesse parity with metro_atlas PDF)
# ---------------------------------------------------------------------------

def _build_svg(
    scene: dict[str, Any],
    layout: Layout,
    *,
    palette: SemanticPalette,
    preset: str,
    color_mode: str,
    strip_mode: str,
    omics_edge_color: str,
    top_rank: int,
    mode: str,
) -> tuple[str, tuple[float, float, float, float]]:
    """把共享场景几何直接序列化为与 ``metro_atlas`` PDF 1:1 精细对齐的语义化 SVG 字符串。"""
    cmap = palette.diverging_cmap()
    pdict = _palette_dict(preset)
    reactions: pd.DataFrame = scene["reactions"]
    aggs: pd.DataFrame = scene["aggs"]
    nodes_df: pd.DataFrame = scene["nodes_df"]
    active_regions: set[str] = scene["active_regions"]
    active_species: set[str] = scene["active_species"]
    active_node_positions = scene["active_node_positions"]
    comp_boxes: dict[str, tuple[float, float, float, float]] = scene["comp_boxes"]
    geometries: dict[str, dict] = scene["geometries"]
    group_indices: dict[str, int] = scene["group_indices"]
    multi_groups: set[int] = scene["multi_groups"]
    node_labels: dict[tuple[str, int], tuple[float, float, str, str]] = scene["node_labels"]
    positions: dict[str, tuple[float, float]] = scene["positions"]
    bx0, by0, bx1, by1 = scene["bounds"]
    font_scale: float = scene["font_scale"]
    unit_in: float = scene["unit_in"]
    is_econ = preset == "economist"
    card_bg = pdict["strip_card_bg"]

    vb = (bx0 - 1.2, by0 - 1.2, (bx1 - bx0) + 2.4, (by1 - by0) + 2.4)
    svg: list[str] = [
        f'<svg id="ow-svg" xmlns="http://www.w3.org/2000/svg" viewBox="{vb[0]:.2f} {vb[1]:.2f} {vb[2]:.2f} {vb[3]:.2f}" '
        f'preserveAspectRatio="xMidYMid meet">',
        f'<rect id="ow-canvas-bg" x="{vb[0] - 80:.2f}" y="{vb[1] - 80:.2f}" width="{vb[2] + 160:.2f}" height="{vb[3] + 160:.2f}" fill="{palette.surface}"/>',
    ]

    # 1. Regions (delicate PDF-calibrated border; Economist top black rule hidden in metro_atlas)
    reg_rx = 0.25 if is_econ else 2.0
    reg_lw = _pt(0.45 if is_econ else 0.50, unit_in)
    rule_lw = _pt(0.80, unit_in)
    rule_disp = "inline" if is_econ else "none"
    svg.append('<g id="ow-layer-regions">')
    for rid, region in layout.regions.items():
        if rid not in active_regions or region.box is None or rid == "oxphos":
            continue
        x0, y0, x1, y1 = region.box
        svg.append(
            f'<g class="ow-region" data-region="{_esc(rid)}">'
            f'<rect class="ow-region-box" x="{x0:.2f}" y="{y0:.2f}" width="{x1 - x0:.2f}" height="{y1 - y0:.2f}" '
            f'rx="{reg_rx:.2f}" fill="{pdict["region_face"]}" stroke="{pdict["region_edge"]}" stroke-width="{reg_lw:.3f}"/>'
            f'<line class="ow-region-rule" x1="{x0:.2f}" y1="{y0:.2f}" x2="{x1:.2f}" y2="{y0:.2f}" '
            f'stroke="{palette.ink}" stroke-width="{rule_lw:.3f}" style="display:{rule_disp}"/>'
            f"</g>"
        )
    svg.append("</g>")

    # 2. Compartments (Mitochondria / Peroxisome)
    svg.append('<g id="ow-layer-compartments">')
    comp_lbl_fs = _pt(5.0 * font_scale, unit_in)
    for cid, box in comp_boxes.items():
        comp = layout.compartments[cid]
        x0, y0, x1, y1 = box
        if cid == "p":
            face, edge, inner, ink = pdict["perox_face"], pdict["perox_edge"], pdict["perox_inner"], pdict["perox_ink"]
            lx, ly, anchor = x1 - 1.3, y0 + 1.35, "end"
            crx = 0.35 if is_econ else 2.4
            clw = _pt(0.75, unit_in)
        else:
            face, edge, inner, ink = pdict["mito_face"], pdict["mito_edge"], pdict["mito_inner"], pdict["mito_ink"]
            if "oxphos" in active_regions:
                lx, ly, anchor = (x0 + x1) * 0.5, y0 + 1.85, "middle"
            else:
                lx, ly, anchor = x1 - 1.3, y1 - 1.15, "end"
            crx = 0.35 if is_econ else 2.6
            clw = _pt(0.80, unit_in)
        svg.append(
            f'<g class="ow-comp" data-comp="{_esc(cid)}">'
            f'<rect class="ow-comp-box" x="{x0:.2f}" y="{y0:.2f}" width="{x1 - x0:.2f}" height="{y1 - y0:.2f}" '
            f'rx="{crx:.2f}" fill="{face}" stroke="{edge}" stroke-width="{clw:.3f}"/>'
            f'<rect class="ow-comp-inner" x="{x0 + 0.65:.2f}" y="{y0 + 0.65:.2f}" width="{(x1 - x0) - 1.30:.2f}" height="{(y1 - y0) - 1.30:.2f}" '
            f'rx="{max(0.2, crx - 0.35):.2f}" fill="none" stroke="{inner}" stroke-width="{_pt(0.55, unit_in):.3f}" '
            f'stroke-dasharray="{_pt(3.5, unit_in):.2f} {_pt(2.4, unit_in):.2f}"/>'
            f'<text class="ow-comp-lbl" x="{lx:.2f}" y="{ly:.2f}" font-size="{comp_lbl_fs:.3f}" font-weight="700" '
            f'font-style="italic" fill="{ink}" text-anchor="{anchor}" dominant-baseline="central">{_esc(comp.get("label", ""))}</text>'
            f"</g>"
        )
    if "oxphos" in active_regions and "oxphos" in layout.regions and layout.regions["oxphos"].box is not None:
        x0, y0, x1, y1 = layout.regions["oxphos"].box
        ox_rx = 0.25 if is_econ else 1.6
        svg.append(
            f'<g class="ow-region" data-region="oxphos">'
            f'<rect class="ow-region-box ow-oxphos-box" x="{x0:.2f}" y="{y0:.2f}" width="{x1 - x0:.2f}" height="{y1 - y0:.2f}" '
            f'rx="{ox_rx:.2f}" fill="{pdict["oxphos_face"]}" stroke="{pdict["oxphos_edge"]}" stroke-width="{_pt(0.55, unit_in):.3f}"/>'
            f'<line class="ow-region-rule" x1="{x0:.2f}" y1="{y0:.2f}" x2="{x1:.2f}" y2="{y0:.2f}" '
            f'stroke="{palette.ink2}" stroke-width="{rule_lw:.3f}" style="display:{rule_disp}"/>'
            f"</g>"
        )
    svg.append("</g>")

    # 3. Region Headings (Metro colored dot in metro_atlas; square bar in economist)
    svg.append('<g id="ow-layer-region-titles">')
    reg_title_fs = _pt(5.4 * font_scale, unit_in)
    reg_halo_sw = _pt(2.0, unit_in)
    for rid, region in layout.regions.items():
        if rid not in active_regions:
            continue
        lx, ly = region.label_xy
        metro_accent = pathway_family_color(rid, fallback_accent=region.accent) if color_mode == "pathway" else (region.accent or palette.ink2)
        econ_accent = pathway_family_color(rid, fallback_accent=region.accent) if color_mode == "pathway" else palette.ink
        in_mito = rid in ("pyruvate_tca", "ketone", "oxphos")
        plate_bg = pdict["mito_face"] if in_mito else pdict["region_face"]
        tx_val = (lx + 0.90) if is_econ else (lx + 1.45)
        svg.append(
            f'<g class="ow-region-title" data-region="{_esc(rid)}" data-in-mito="{1 if in_mito else 0}" '
            f'data-accent-metro="{_esc(metro_accent)}" data-accent-econ="{_esc(econ_accent)}" style="cursor:pointer">'
            f'<circle class="ow-region-dot" cx="{lx + 0.55:.2f}" cy="{ly:.2f}" r="0.50" '
            f'fill="{metro_accent}" stroke="{palette.surface}" stroke-width="{_pt(0.45, unit_in):.3f}" '
            f'style="display:{"none" if is_econ else "inline"}"/>'
            f'<rect class="ow-region-accent" x="{lx + 0.15:.2f}" y="{ly - 0.46:.2f}" width="0.35" height="0.92" '
            f'fill="{econ_accent}" style="display:{"inline" if is_econ else "none"}"/>'
            f'<text class="ow-region-text" x="{tx_val:.2f}" data-x-metro="{lx + 1.45:.2f}" data-x-econ="{lx + 0.90:.2f}" '
            f'y="{ly:.2f}" font-size="{reg_title_fs:.3f}" font-weight="700" fill="{palette.ink}" '
            f'dominant-baseline="central" paint-order="stroke" stroke="{plate_bg}" stroke-width="{reg_halo_sw:.3f}" '
            f'stroke-linejoin="round">{_esc(region.label)}</text>'
            f"</g>"
        )
    svg.append("</g>")

    # 4. Linked Flux Conduits + Reaction Tracks
    draw_order = (
        reactions["draw_style"].map({"not_assessed": 0, "omics_only": 1, "flux": 2}).fillna(0)
        + reactions["key"].eq("l1").astype(int)
        + reactions["value"].notna().astype(int) * 2
        + reactions["key"].eq("key").astype(int) * 2
    )
    ordered_rxns = reactions.assign(_z=draw_order).sort_values("_z", kind="stable")

    # Conduits (exact 8.2pt fill / 9.2pt edge)
    cond_edge_sw = _pt(9.2, unit_in)
    cond_fill_sw = _pt(8.2, unit_in)
    svg.append('<g id="ow-layer-conduits">')
    for _, row in ordered_rxns.iterrows():
        fid = str(row["element_id"])
        fg = str(row["flux_group"] or "").strip()
        gi = group_indices.get(fg) if (fg and row["F_state"] in ("changed", "suggestive")) else None
        if gi not in multi_groups:
            continue
        mpaths = geometries[fid].get("metro_paths") or []
        for mp in mpaths:
            pts = fillet_polyline(trim_polyline(mp["points"], mp["trim0"], mp["trim1"]), radius=1.15)
            d = _pts_path(pts)
            if d:
                svg.append(
                    f'<g class="ow-conduit" data-flux-group="{_esc(fg)}" data-rxn="{_esc(fid)}">'
                    f'<path class="ow-conduit-edge" d="{d}" fill="none" stroke="{palette.coupled_edge}" '
                    f'stroke-width="{cond_edge_sw:.3f}" stroke-linecap="round" stroke-linejoin="round" opacity="0.24"/>'
                    f'<path class="ow-conduit-fill" d="{d}" fill="none" stroke="{palette.coupled_track}" '
                    f'stroke-width="{cond_fill_sw:.3f}" stroke-linecap="round" stroke-linejoin="round" opacity="0.86"/>'
                    f"</g>"
                )
    svg.append("</g>")

    # Edges (exact pt widths matching _draw_edge)
    halo_extra = _pt(1.9, unit_in)
    note_fs = _pt(NOTE_SIZE * font_scale, unit_in)
    svg.append('<g id="ow-layer-edges">')
    for _, row in ordered_rxns.iterrows():
        fid = str(row["element_id"])
        geom = geometries[fid]
        draw_style = str(row["draw_style"] or "")
        claim = str(row["flux_claim"] or "")
        fg = str(row["flux_group"] or "").strip()
        strong = row["key"] == "key"
        colored = draw_style == "flux" and pd.notna(row["value"])
        enz_dir, enz_mag = _enzyme_inferred_direction(row)
        is_enz = (not colored) and (enz_dir != 0)
        hollow = (draw_style == "omics_only") or is_enz
        is_na = draw_style == "not_assessed"
        physio = int(row["physio_sign"]) if pd.notna(row["physio_sign"]) else 0
        sign_chg = bool(row["sign_change"])
        both = sign_chg or (physio == 0 and bool(row["reversible"]))
        chev_dir = -1 if (not both and physio < 0) else 1

        col_hex = _to_hex(
            _edge_color(
                row, palette, cmap, metro=True, color_mode=color_mode,
                omics_edge_color=omics_edge_color, layout=layout,
            )
        )
        if colored:
            val_num = float(row["value"])
            edge_role = "mfa_dir"
        elif is_enz:
            val_num = enz_dir * max(0.4, enz_mag)
            edge_role = "enz_up" if enz_dir > 0 else "enz_down"
        elif draw_style == "flux":
            val_num = None
            edge_role = "mfa_grey"
        elif draw_style == "omics_only":
            val_num = None
            edge_role = "omics_grey"
        else:
            val_num = None
            edge_role = "not_assessed"

        if is_na:
            sw = _pt(0.95, unit_in)
        elif hollow:
            sw = _pt(2.85 if (row["key"] == "l1" or is_enz) else 2.20, unit_in)
        else:
            sw = _pt(3.70 if strong else (3.00 if colored else 2.15), unit_in)
        core_sw = max(_pt(0.75, unit_in), sw * 0.44)
        val_attr = f"{val_num:.3f}" if val_num is not None else ""

        svg.append(
            f'<g class="ow-edge" data-rxn="{_esc(fid)}" data-region="{_esc(row["region"])}" '
            f'data-claim="{_esc(claim)}" data-draw="{_esc(draw_style)}" data-erole="{edge_role}" '
            f'data-val="{val_attr}" data-flux-group="{_esc(fg)}" data-sign-change="{1 if sign_chg else 0}" '
            f'data-hollow="{1 if hollow else 0}">'
        )
        mpaths = geom.get("metro_paths") or []
        for mp in mpaths:
            pts = fillet_polyline(trim_polyline(mp["points"], mp["trim0"], mp["trim1"]), radius=1.15)
            d = _pts_path(pts)
            if not d:
                continue
            svg.append(f'<path class="ow-hit" d="{d}" fill="none" stroke="transparent" stroke-width="2.2" stroke-linecap="round"/>')
            if not is_na:
                svg.append(
                    f'<path class="ow-edge-halo" d="{d}" fill="none" stroke="{palette.surface}" '
                    f'stroke-width="{sw + halo_extra:.3f}" stroke-linecap="round" stroke-linejoin="round" opacity="0.90"/>'
                )
            svg.append(
                f'<path class="ow-edge-main" d="{d}" fill="none" stroke="{col_hex}" '
                f'stroke-width="{sw:.3f}" stroke-linecap="round" stroke-linejoin="round"/>'
            )
            if hollow:
                svg.append(
                    f'<path class="ow-edge-core" d="{d}" fill="none" stroke="{palette.surface}" '
                    f'stroke-width="{core_sw:.3f}" stroke-linecap="round" stroke-linejoin="round"/>'
                )
            if not is_na:
                chev_col = col_hex if hollow else palette.surface
                svg.append(_svg_chevrons(pts, chev_dir, both, chev_col, sw, unit_in))
        if len(mpaths) > 1 and not is_na:
            hx, hy = geom["hub"]
            svg.append(
                f'<circle class="ow-edge-hub" cx="{hx:.2f}" cy="{hy:.2f}" r="0.34" '
                f'fill="{col_hex}" stroke="{palette.surface}" stroke-width="{_pt(0.45, unit_in):.3f}"/>'
            )
        if sign_chg:
            hx, hy = geom["hub"]
            svg.append(_svg_sign_reversal_badge(hx, hy, 0.78, palette, unit_in))
        if geom.get("note"):
            hx, hy = geom["hub"]
            svg.append(
                f'<text class="ow-edge-note" x="{hx:.2f}" y="{hy + 0.55:.2f}" font-size="{note_fs:.3f}" '
                f'font-style="italic" fill="{palette.muted}" text-anchor="middle" dominant-baseline="hanging" '
                f'paint-order="stroke" stroke="{palette.surface}" stroke-width="{_pt(1.6, unit_in):.3f}" '
                f'stroke-linejoin="round">{_esc(geom["note"])}</text>'
            )
        svg.append("</g>")
    svg.append("</g>")

    # 5. Aggregate Cards (exact PDF parity: delicate borders, adaptive title, 4-omics strip + P/R/M squares)
    header_icons = OMICS_LAYER_ICONS if strip_mode == "omics" else AXIS_LAYER_ICONS
    svg.append('<g id="ow-layer-aggregates">')
    for _, row in aggs.iterrows():
        aid = str(row["element_id"])
        agg = layout.aggregates[aid]
        x, y = agg.xy
        w, h = _agg_card_dims(agg)
        faint = row["draw_style"] == "not_assessed"
        is_key = row["key"] == "key"
        is_l1 = row["key"] == "l1"
        claim = str(row["flux_claim"]) if row["flux_claim"] in ("L1", "L2", "L3") else ""
        ec = palette.ink if is_key else (palette.ink2 if is_l1 else (palette.faint if faint else palette.card_edge))
        lw = _pt(1.25 if is_key else (0.88 if is_l1 else 0.58), unit_in)
        dash = f' stroke-dasharray="{_pt(2.4, unit_in):.3f} {_pt(1.3, unit_in):.3f}"' if is_l1 else ""
        arx = 0.16 if is_econ else 0.55
        svg.append(
            f'<g class="ow-agg" data-agg="{_esc(aid)}" data-region="{_esc(row["region"])}" data-claim="{_esc(claim)}">'
            f'<rect class="ow-card-bg" x="{x - w / 2:.2f}" y="{y - h / 2:.2f}" width="{w:.2f}" height="{h:.2f}" '
            f'rx="{arx:.2f}" fill="{card_bg}" stroke="{ec}" stroke-width="{lw:.3f}"{dash}/>'
        )
        title_str = str(row["label"])
        title_y = y - 1.02
        if claim:
            b_sc = 0.78 * font_scale
            badge_w = tier_badge_width(claim, scale=b_sc, unit_in=unit_in)
            badge_x = x + w / 2 - badge_w - 0.36
            avail = w - badge_w - 1.10
            t_pt = 3.85 * font_scale
            bold_t = is_key or is_l1
            while t_pt > 2.65 and _text_width(title_str, t_pt, bold=bold_t, unit=unit_in) > avail:
                t_pt -= 0.16
            svg.append(
                f'<text class="ow-card-title" data-max-x="{badge_x - 0.22:.2f}" x="{x - w / 2 + 0.44:.2f}" '
                f'y="{title_y:.2f}" font-size="{_pt(t_pt * 0.93, unit_in):.3f}" '
                f'font-weight="{700 if bold_t else 600}" fill="{palette.muted if faint else palette.ink}" '
                f'dominant-baseline="central">{_esc(title_str)}</text>'
            )
            b_svg, _ = _svg_tier_badge(badge_x, title_y, claim, b_sc, palette, preset, unit_in)
            svg.append(b_svg)
        else:
            t_pt = 4.0 * font_scale
            while t_pt > 2.80 and _text_width(title_str, t_pt, bold=False, unit=unit_in) > (w - 0.90):
                t_pt -= 0.16
            svg.append(
                f'<text class="ow-card-title" data-max-w="{w - 0.72:.2f}" x="{x:.2f}" y="{title_y:.2f}" '
                f'font-size="{_pt(t_pt * 0.93, unit_in):.3f}" font-weight="600" '
                f'fill="{palette.muted if faint else palette.ink}" text-anchor="middle" dominant-baseline="central">{_esc(title_str)}</text>'
            )

        # Bottom row: exact PDF layout (4-omics strip + P/R/M icon & log2FC color squares)
        compact = w < 10.5
        mc = 0.72 if compact else 0.82
        m_hdr = 0.44 if (is_key or is_l1) else 0.14
        strip_w = 4 * mc
        sq_w = 0.84 if compact else 0.96
        prm_step = 1.55 if compact else 2.15
        prm_w = 2 * prm_step + sq_w + 0.65
        gap = max(0.55, min(1.15, w - strip_w - prm_w - 1.0))
        content_w = strip_w + gap + prm_w
        left_x = x - content_w / 2.0
        row_y = y + 0.68

        states = _extract_strip_states(row, strip_mode=strip_mode)
        total_h = mc + m_hdr
        st_top = row_y - total_h / 2.0
        cells_y = st_top + m_hdr
        hdr_bg = _to_hex(tint_color(palette.ink2, 0.09, base=card_bg))
        svg.append(
            f'<rect class="ow-strip-frame" x="{left_x:.2f}" y="{st_top:.2f}" width="{strip_w:.2f}" height="{total_h:.2f}" '
            f'rx="0.08" fill="{card_bg}" stroke="{palette.card_edge}" stroke-width="{_pt(0.40, unit_in):.3f}"/>'
            f'<rect class="ow-strip-hdr" x="{left_x + 0.02:.2f}" y="{st_top + 0.02:.2f}" width="{strip_w - 0.04:.2f}" height="{max(0.04, m_hdr - 0.02):.2f}" fill="{hdr_bg}"/>'
        )
        for idx, (st, dr) in enumerate(states):
            cx_i = left_x + idx * mc
            if (is_key or is_l1) and idx < len(header_icons):
                svg.append(_svg_layer_icon(header_icons[idx], cx_i + mc / 2, st_top + m_hdr * 0.52, min(mc * 0.62, m_hdr * 0.82), palette.ink2, unit_in))
            if idx > 0:
                svg.append(f'<line x1="{cx_i:.2f}" y1="{st_top + 0.03:.2f}" x2="{cx_i:.2f}" y2="{cells_y + mc - 0.03:.2f}" stroke="{palette.hairline}" stroke-width="{_pt(0.30, unit_in):.3f}"/>')
            svg.append(_svg_status_cell(cx_i, cells_y, mc, st, dr, palette, unit_in=unit_in, card_base=card_bg))

        prm_x0 = left_x + strip_w + gap + 0.62
        for i, (icon_kind, value_col, state_col) in enumerate((
            ("protein", "P_value", "P_state"),
            ("mrna", "R_value", "R_layer_state"),
            ("metabolite", "Mx_value", "Mx_state"),
        )):
            sx = prm_x0 + i * prm_step
            val_raw = row[value_col] if value_col in row.index else np.nan
            val_n = _num(val_raw, 3) if pd.notna(val_raw) else None
            st_str = str(row[state_col]) if (state_col in row.index and pd.notna(row[state_col])) else "not_measured"
            face = _to_hex(cmap(DEFAULT_NORM(val_n))) if val_n is not None else palette.missing_fill
            svg.append(_svg_layer_icon(icon_kind, sx - 0.34, row_y, 0.58 if compact else 0.68, palette.muted, unit_in))
            svg.append(
                f'<rect class="ow-prm-sq" data-val="{val_n if val_n is not None else ""}" '
                f'x="{sx:.2f}" y="{row_y - sq_w / 2:.2f}" width="{sq_w:.2f}" height="{sq_w:.2f}" rx="0.14" '
                f'fill="{face}" stroke="{palette.card_edge}" stroke-width="{_pt(0.35, unit_in):.3f}"/>'
            )
            if val_n is None:
                svg.append(_svg_status_cell(sx, row_y - sq_w / 2, sq_w, st_str, 0, palette, unit_in=unit_in, draw_bg=False, card_base=card_bg))
        svg.append("</g>")
    svg.append("</g>")

    # 6. Metabolite Stations / Boundary Pills / Store Pills (exact PDF radius, thin strokes & crisp labels)
    svg.append('<g id="ow-layer-nodes">')
    node_lbl_fs = _pt(4.6 * font_scale * 0.95, unit_in)
    node_halo_sw = _pt(2.0, unit_in)
    station_r = NODE_R * 1.05
    for species, node in layout.nodes.items():
        if species not in active_species:
            continue
        nrow = nodes_df.loc[species] if species in nodes_df.index else None
        for idx, (nx, ny) in enumerate(node.at):
            if active_node_positions is not None and (species, idx) not in active_node_positions:
                continue
            if node.kind == "metabolite":
                val = _num(nrow["value"]) if (nrow is not None and pd.notna(nrow["value"])) else None
                measured = val is not None
                state = str(nrow["display_state"] or "") if nrow is not None else "not_measured"
                ring = str(nrow["ring"] or "") if nrow is not None else ""
                sig = ring in ("selected", "changed")
                fill = _to_hex(cmap(DEFAULT_NORM(val))) if measured else palette.missing_fill
                edge = palette.ink if sig else (palette.ink2 if measured else palette.omics_grey)
                lw = _pt(1.35 if sig else (1.10 if measured else 0.85), unit_in)
                lx, ly, ha, va = node_labels.get((species, idx), (nx, ny + NODE_R + 0.55, "center", "top", None))[:4]
                anchor = {"left": "start", "right": "end"}.get(ha, "middle")
                vbase = {"top": "hanging", "bottom": "ideographic"}.get(va, "central")
                short_lbl = node.label or _species_short(species, str(nrow["label"]) if nrow is not None else species)
                svg.append(
                    f'<g class="ow-node" data-node="{_esc(species)}" data-val="{val if measured else ""}" '
                    f'data-state="{_esc(state)}" data-ring="{_esc(ring)}">'
                    f'<circle class="ow-node-halo" cx="{nx:.2f}" cy="{ny:.2f}" r="{station_r + 0.34:.2f}" fill="{palette.surface}"/>'
                    f'<circle class="ow-node-dot" cx="{nx:.2f}" cy="{ny:.2f}" r="{station_r:.2f}" fill="{fill}" stroke="{edge}" stroke-width="{lw:.3f}"/>'
                )
                if sig:
                    svg.append(
                        f'<circle class="ow-node-ring" cx="{nx:.2f}" cy="{ny:.2f}" r="{station_r + 0.20:.2f}" '
                        f'fill="none" stroke="{palette.ink}" stroke-width="{_pt(1.10, unit_in):.3f}"/>'
                    )
                if not measured:
                    if state in ("not_measured", "no_knowledge"):
                        d = station_r * 0.44
                        svg.append(f'<line x1="{nx - d:.2f}" y1="{ny + d:.2f}" x2="{nx + d:.2f}" y2="{ny - d:.2f}" stroke="{palette.missing_hatch}" stroke-width="{_pt(0.58, unit_in):.3f}"/>')
                    elif state == "not_detected":
                        svg.append(
                            f'<circle cx="{nx:.2f}" cy="{ny:.2f}" r="{station_r * 0.44:.2f}" fill="none" stroke="{palette.muted}" '
                            f'stroke-width="{_pt(0.62 * font_scale, unit_in):.3f}" stroke-dasharray="{_pt(1.5, unit_in):.3f} {_pt(1.2, unit_in):.3f}"/>'
                        )
                    elif state in ("qc_blocked", "insufficient_n"):
                        d = station_r * 0.38
                        off = station_r * 0.16
                        for dx_off in (-off, off):
                            svg.append(
                                f'<line x1="{nx - d + dx_off:.2f}" y1="{ny + d:.2f}" x2="{nx + d + dx_off:.2f}" y2="{ny - d:.2f}" '
                                f'stroke="{palette.muted}" stroke-width="{_pt(0.54 * font_scale, unit_in):.3f}" stroke-linecap="round"/>'
                            )
                svg.append(
                    f'<text class="ow-node-lbl" x="{lx:.2f}" y="{ly:.2f}" font-size="{node_lbl_fs:.3f}" '
                    f'font-weight="{700 if sig else 500}" fill="{palette.ink if (sig or measured) else palette.muted}" '
                    f'text-anchor="{anchor}" dominant-baseline="{vbase}" paint-order="stroke" stroke="{palette.surface}" '
                    f'stroke-width="{node_halo_sw:.3f}" stroke-linejoin="round">{_esc(short_lbl)}</text>'
                    f"</g>"
                )
            elif node.kind == "boundary":
                name = node.label or _species_short(species)
                text = f"{name} (Blood)" if species.endswith("_e") else name
                fs_pt = 3.85 * font_scale
                bw = max(3.2, _text_width(text, fs_pt, bold=True, unit=unit_in) + 1.10)
                bh = 1.44
                brx = 0.22 if is_econ else 0.72
                svg.append(
                    f'<g class="ow-boundary" data-node="{_esc(species)}">'
                    f'<rect class="ow-pill-bg" x="{nx - bw / 2:.2f}" y="{ny - bh / 2:.2f}" width="{bw:.2f}" height="{bh:.2f}" rx="{brx:.2f}" '
                    f'fill="{palette.boundary_face}" stroke="{palette.boundary_edge}" stroke-width="{_pt(0.85, unit_in):.3f}"/>'
                    f'<text class="ow-pill-txt" x="{nx:.2f}" y="{ny:.2f}" font-size="{_pt(fs_pt * 0.92, unit_in):.3f}" font-weight="700" '
                    f'fill="{palette.boundary_ink}" text-anchor="middle" dominant-baseline="central">{_esc(text)}</text>'
                    f"</g>"
                )
            else:  # store / sink
                text = node.label or _species_short(species)
                fs_pt = 3.95 * font_scale
                bw = max(3.8, _text_width(text, fs_pt, bold=True, unit=unit_in) + 1.20)
                bh = 1.48
                srx = 0.20 if is_econ else 0.72
                svg.append(
                    f'<g class="ow-store" data-node="{_esc(species)}">'
                    f'<rect class="ow-pill-bg" x="{nx - bw / 2:.2f}" y="{ny - bh / 2:.2f}" width="{bw:.2f}" height="{bh:.2f}" rx="{srx:.2f}" '
                    f'fill="{palette.store_face}" stroke="{palette.store_edge}" stroke-width="{_pt(0.85, unit_in):.3f}"/>'
                    f'<text class="ow-pill-txt" x="{nx:.2f}" y="{ny:.2f}" font-size="{_pt(fs_pt * 0.92, unit_in):.3f}" font-weight="700" '
                    f'fill="{palette.store_ink}" text-anchor="middle" dominant-baseline="central">{_esc(text)}</text>'
                    f"</g>"
                )
    svg.append("</g>")

    # 7. Reaction Cards / 4-Omics Strips (exact pt-calibrated borders, font sizes & zero-overflow widths)
    svg.append('<g id="ow-layer-strips">')
    total_groups = max(len(group_indices), 1)
    compact_fs = _pt(3.85 * font_scale, unit_in)
    for _, row in reactions.iterrows():
        fid = str(row["element_id"])
        x0, y0 = positions[fid]
        hub = geometries[fid]["hub"]
        claim = str(row["flux_claim"] or "")
        fg = str(row["flux_group"] or "").strip()
        draw_style = str(row["draw_style"] or "")
        info = _strip_badges(row, group_indices, top_rank=top_rank, show_rank=(mode != "publication"))
        label = info["label"]

        if _is_compact_label_only(row, mode):
            is_na = draw_style == "not_assessed"
            svg.append(
                f'<g class="ow-strip ow-strip-compact" data-rxn="{_esc(fid)}" data-region="{_esc(row["region"])}" '
                f'data-claim="{_esc(claim)}" data-draw="{_esc(draw_style)}" data-flux-group="{_esc(fg)}">'
                f'<text class="ow-rxn-lbl" x="{x0 + 0.15:.2f}" y="{y0:.2f}" font-size="{compact_fs:.3f}" '
                f'fill="{palette.muted if is_na else palette.ink2}" dominant-baseline="central" '
                f'paint-order="stroke" stroke="{palette.surface}" stroke-width="{_pt(1.8, unit_in):.3f}" '
                f'stroke-linejoin="round">{_esc(label)}</text>'
                f"</g>"
            )
            continue

        w, h = _strip_size(row, group_indices, top_rank=top_rank, mode=mode, font_scale=font_scale, unit_in=unit_in)
        is_key = row["key"] == "key"
        is_l1 = row["key"] == "l1"
        is_hi = row["key"] in ("key", "l1", "top")
        bold = is_key or is_l1
        ec = palette.ink if is_key else (palette.ink2 if is_l1 else palette.card_edge)
        lw = _pt(1.25 if is_key else (0.88 if is_l1 else 0.48), unit_in)
        dash = f' stroke-dasharray="{_pt(2.4, unit_in):.3f} {_pt(1.3, unit_in):.3f}"' if is_l1 else ""
        crx = 0.12 if is_econ else 0.38
        card_left = x0 - 0.12
        card_top = y0 - h / 2.0

        # Precompute exact inner layout so card_w mathematically encloses every label & badge
        st_x = x0 + 0.08
        rx_x = st_x + 4 * CELL + 0.34
        lbl_pt = 4.35 * font_scale
        lbl_w = _text_width(label, lbl_pt, bold=bold, unit=unit_in)
        bx = rx_x + lbl_w + 0.48
        b_scale = 1.02 * font_scale
        badge_parts: list[str] = []
        content_right = rx_x + lbl_w
        if info["tier"] in ("L1", "L2", "L3"):
            bsvg, bw = _svg_tier_badge(
                bx, y0, info["tier"], b_scale, palette, preset, unit_in,
                label=info.get("tier_label") or None,
            )
            badge_parts.append(bsvg)
            content_right = bx + bw
            bx += bw + 0.34
        if info["group"] is not None:
            csvg, cw = _svg_coupled_badge(bx, y0, int(info["group"]), total_groups, b_scale, palette, preset, unit_in)
            badge_parts.append(csvg)
            content_right = bx + cw
        card_w = max(w + 0.12, (content_right + 0.38) - card_left)

        # Leader line if strip shifted away from hub
        end = min(((x0, y0), (card_left + card_w, y0)), key=lambda p: math.dist(p, hub))
        leader_svg = ""
        if math.dist(end, hub) > 2.0:
            leader_svg = (
                f'<line class="ow-leader" x1="{hub[0]:.2f}" y1="{hub[1]:.2f}" x2="{end[0]:.2f}" y2="{end[1]:.2f}" '
                f'stroke="{palette.omics_grey}" stroke-width="{_pt(0.40, unit_in):.3f}" '
                f'stroke-dasharray="{_pt(1.4, unit_in):.3f} {_pt(1.4, unit_in):.3f}"/>'
                f'<circle cx="{hub[0]:.2f}" cy="{hub[1]:.2f}" r="0.20" fill="{palette.edge_grey}" '
                f'stroke="{palette.surface}" stroke-width="{_pt(0.25, unit_in):.3f}"/>'
            )

        parts = [
            f'<g class="ow-strip" data-rxn="{_esc(fid)}" data-region="{_esc(row["region"])}" '
            f'data-claim="{_esc(claim)}" data-draw="{_esc(draw_style)}" data-flux-group="{_esc(fg)}" '
            f'data-tier-style="{"key" if is_key else ("l1" if is_l1 else "std")}">',
            leader_svg,
            f'<rect class="ow-card-bg" x="{card_left:.2f}" y="{card_top:.2f}" width="{card_w:.2f}" height="{h:.2f}" '
            f'rx="{crx:.2f}" fill="{card_bg}" stroke="{ec}" stroke-width="{lw:.3f}"{dash} opacity="0.98"/>',
        ]

        # Top-left #rank tab (exact draw_rank_tab geometry)
        if info["rank"] is not None:
            rk = int(info["rank"])
            r_sc = 0.96 * font_scale
            rw = (1.68 if rk < 10 else 2.05) * r_sc
            rh = 0.92 * r_sc
            top3 = rk <= 3
            rbg = (ECONOMIST_RED if is_econ else palette.ink) if top3 else palette.card_bg
            rec = rbg if top3 else palette.ink2
            rfg = "#FFFFFF" if (top3 and is_econ) else (palette.surface if top3 else palette.ink)
            ry0 = card_top - rh * 0.62
            rx0 = card_left - 0.05 * r_sc
            parts.append(
                f'<rect class="ow-rank-tab" data-top3="{1 if top3 else 0}" x="{rx0:.2f}" y="{ry0:.2f}" width="{rw:.2f}" height="{rh:.2f}" '
                f'rx="{(0.10 if is_econ else 0.28) * r_sc:.2f}" fill="{rbg}" stroke="{rec}" '
                f'stroke-width="{_pt((0.55 if top3 else 0.45) * r_sc, unit_in):.3f}"/>'
                f'<text class="ow-rank-txt" data-top3="{1 if top3 else 0}" x="{rx0 + rw / 2:.2f}" y="{ry0 + rh / 2 + 0.02 * r_sc:.2f}" '
                f'font-size="{_pt(3.15 * r_sc, unit_in):.3f}" font-weight="700" fill="{rfg}" '
                f'text-anchor="middle" dominant-baseline="central">#{rk}</text>'
            )

        # 4-omics icon strip (exact draw_axis_strip proportions)
        states = _extract_strip_states(row, strip_mode=strip_mode)
        hdr_h = AXIS_HDR_H if is_hi else AXIS_BAR_H
        total_h = CELL + hdr_h
        st_top = y0 - total_h / 2.0
        cells_y = st_top + hdr_h
        hdr_bg = _to_hex(tint_color(palette.ink2, 0.09, base=card_bg))
        parts.append(
            f'<rect class="ow-strip-frame" x="{st_x:.2f}" y="{st_top:.2f}" width="{4 * CELL:.2f}" height="{total_h:.2f}" '
            f'rx="0.08" fill="{card_bg}" stroke="{palette.card_edge}" stroke-width="{_pt(0.42, unit_in):.3f}"/>'
            f'<rect class="ow-strip-hdr" x="{st_x + 0.02:.2f}" y="{st_top + 0.02:.2f}" width="{4 * CELL - 0.04:.2f}" height="{max(0.04, hdr_h - 0.02):.2f}" fill="{hdr_bg}"/>'
        )
        for i, (st, dr) in enumerate(states):
            cx_i = st_x + i * CELL
            if is_hi and i < len(header_icons):
                parts.append(_svg_layer_icon(header_icons[i], cx_i + CELL / 2, st_top + hdr_h * 0.52, min(CELL * 0.64, hdr_h * 0.84), palette.ink2, unit_in))
            if i > 0:
                parts.append(f'<line x1="{cx_i:.2f}" y1="{st_top + 0.03:.2f}" x2="{cx_i:.2f}" y2="{cells_y + CELL - 0.03:.2f}" stroke="{palette.hairline}" stroke-width="{_pt(0.32, unit_in):.3f}"/>')
            parts.append(_svg_status_cell(cx_i, cells_y, CELL, st, dr, palette, unit_in=unit_in, card_base=card_bg))

        # Reaction label + right-side badges
        parts.append(
            f'<text class="ow-rxn-lbl" x="{rx_x:.2f}" y="{y0:.2f}" '
            f'font-size="{_pt(lbl_pt * 0.94, unit_in):.3f}" font-weight="{700 if bold else 500}" '
            f'fill="{palette.ink if bold else palette.ink2}" dominant-baseline="central">{_esc(label)}</text>'
        )
        parts.extend(badge_parts)
        parts.append("</g>")
        svg.append("".join(parts))
    svg.append("</g>")

    svg.append("</svg>")
    return "\n".join(svg), vb


# ---------------------------------------------------------------------------
# Compact Data Payload Extraction (Elements + Members + Per-Animal Samples)
# ---------------------------------------------------------------------------

def _extract_payload(
    view: pd.DataFrame,
    layout: Layout,
    regions: pd.DataFrame,
    scene: dict[str, Any],
    *,
    title: str,
    caption: list[str],
    style: str,
    color_mode: str,
    strip_mode: str,
    omics_edge_color: str,
    vb: tuple[float, float, float, float],
    extra_fusion: pd.DataFrame | None = None,
    signed_modules: pd.DataFrame | None = None,
    members: pd.DataFrame | None = None,
    samples: pd.DataFrame | None = None,
    sample_contrasts: pd.DataFrame | None = None,
    tissue: str = "",
    contrast_id: str = "",
    source: str = "",
) -> dict[str, Any]:
    rxn_lookup = pathway_reactions()
    reactions_df: pd.DataFrame = scene["reactions"]
    aggs_df: pd.DataFrame = scene["aggs"]
    nodes_df: pd.DataFrame = scene["nodes_df"]
    active_regions: set[str] = scene["active_regions"]

    # 1. Submap quick-focus bounding boxes
    _SHORT_SUBMAP_LABELS = {
        "central_carbon": "Central Carbon",
        "fatty_acid_ketone": "Fatty Acid & Ketone",
        "glycerolipid": "Glycerolipid",
        "amino_acid": "Amino Acids",
        "urea_creatine": "Urea & Creatine",
        "nucleotide": "Nucleotide",
        "redox_one_carbon": "Redox & 1C",
    }
    submap_views: list[dict[str, Any]] = [
        {"id": "all", "label": "Full Atlas", "box": [round(v, 2) for v in vb]},
    ]
    for sub in SUBMAP_CATALOG:
        pts_x: list[float] = []
        pts_y: list[float] = []
        for rid in sub.regions:
            if rid in active_regions and rid in layout.regions and layout.regions[rid].box is not None:
                rx0, ry0, rx1, ry1 = layout.regions[rid].box
                pts_x.extend([rx0, rx1])
                pts_y.extend([ry0, ry1])
        if pts_x and pts_y:
            pad = 2.8
            submap_views.append({
                "id": sub.id,
                "label": _SHORT_SUBMAP_LABELS.get(sub.id, sub.title.split(" / ")[0]),
                "full_title": sub.title,
                "regions": list(sub.regions),
                "box": [
                    round(min(pts_x) - pad, 2),
                    round(min(pts_y) - pad, 2),
                    round((max(pts_x) - min(pts_x)) + 2 * pad, 2),
                    round((max(pts_y) - min(pts_y)) + 2 * pad, 2),
                ],
            })

    # 2. Regions dashboard
    region_list: list[dict[str, Any]] = []
    for row in regions.itertuples():
        rid = str(row.region)
        reg_meta = layout.regions.get(rid)
        rbox = list(reg_meta.box) if (reg_meta and reg_meta.box is not None) else None
        sig_items, ns_note = _format_fusion_items(getattr(row, "fusion", ""))
        region_list.append({
            "id": rid,
            "label": "Extra-Network Modules" if rid == "other" else str(row.label),
            "kegg": str(getattr(row, "kegg", "") or ""),
            "n_drawn": int(getattr(row, "n_drawn", 0) or 0),
            "n_flux": int(getattr(row, "n_flux", 0) or 0),
            "n_identifiable": int(getattr(row, "n_identifiable", 0) or 0),
            "n_omics_only": int(getattr(row, "n_omics_only", 0) or 0),
            "n_aggregates": int(getattr(row, "n_aggregates", 0) or 0),
            "n_L2": int(getattr(row, "n_L2", 0) or 0),
            "n_L1": int(getattr(row, "n_L1", 0) or 0),
            "fusion_sig": [{"text": t, "strong": bool(s)} for t, s in sig_items[:3]],
            "fusion_note": ns_note or "",
            "box": [round(v, 2) for v in rbox] if rbox else None,
        })

    # 3. Elements metadata (reactions, aggregates, metabolites)
    elements: dict[str, dict[str, Any]] = {}
    search_index: list[dict[str, Any]] = []

    for _, row in reactions_df.iterrows():
        fid = str(row["element_id"])
        nr = rxn_lookup.get(fid)
        eq = ""
        if nr is not None:
            arrow = " ⇌ " if nr.reversible else " → "
            subs = " + ".join(_species_short(s) for s in nr.substrates)
            prods = " + ".join(_species_short(p) for p in nr.products)
            eq = f"{subs}{arrow}{prods}"
        enz_dir, enz_mag = _enzyme_inferred_direction(row)
        hub = scene["geometries"][fid]["hub"]
        raw_ft = str(row["flux_tracers"]).strip() if ("flux_tracers" in row.index and pd.notna(row.get("flux_tracers"))) else ""
        flux_tracers_list = [format_tracer_display(t) for t in raw_ft.split(";") if t.strip()]
        flux_pool_mode = str(row["flux_pooling_mode"]).strip() if ("flux_pooling_mode" in row.index and pd.notna(row.get("flux_pooling_mode"))) else ""
        claim_str = str(row.get("flux_claim") or "")
        tier_label = format_tracer_badge_label(flux_tracers_list, level=claim_str) if claim_str in ("L1", "L2", "L3") else ""
        el = {
            "kind": "reaction",
            "id": fid,
            "short": _short(fid),
            "label": str(row.get("label") or _short(fid)),
            "region": str(row.get("region") or ""),
            "region_label": layout.regions[str(row["region"])].label if str(row.get("region")) in layout.regions else "",
            "equation": eq,
            "draw_style": str(row.get("draw_style") or ""),
            "reversible": bool(row.get("reversible")),
            "physio_sign": int(row["physio_sign"]) if pd.notna(row.get("physio_sign")) else 0,
            "sign_change": bool(row.get("sign_change")),
            "flux_group": str(row.get("flux_group") or "").strip(),
            "flux_tracers": flux_tracers_list,
            "flux_pooling_mode": flux_pool_mode,
            "tier_label": tier_label,
            "kegg_pathway": str(row.get("kegg_pathway") or ""),
            "kegg_reaction": str(row.get("kegg_reaction") or ""),
            "ec_number": str(row.get("ec_number") or ""),
            "flux_claim": claim_str,
            "rank": int(row["rank"]) if pd.notna(row.get("rank")) else None,
            "value": _num(row.get("value"), 3),
            "enz_dir": int(enz_dir),
            "enz_mag": _num(enz_mag, 3),
            "axes": {
                "F": {"state": str(row.get("F_state") or ""), "dir": int(row["F_direction"]) if pd.notna(row.get("F_direction")) else 0},
                "C": {"state": str(row.get("C_state") or ""), "dir": int(row["C_direction"]) if pd.notna(row.get("C_direction")) else 0},
                "R": {"state": str(row.get("R_state") or ""), "dir": int(row["R_direction"]) if pd.notna(row.get("R_direction")) else 0},
                "M": {"state": str(row.get("M_state") or ""), "dir": int(row["M_direction"]) if pd.notna(row.get("M_direction")) else 0},
            },
            "layers": {
                "P": {"val": _num(row.get("P_value"), 3), "state": str(row.get("P_state") or "")},
                "R": {"val": _num(row.get("R_value"), 3), "state": str(row.get("R_layer_state") or "")},
                "M": {"val": _num(row.get("Mx_value"), 3), "state": str(row.get("Mx_state") or "")},
            },
            "f_rule": str(row.get("text") or "").strip() if pd.notna(row.get("text")) else "",
            "f_flags": [
                f.strip() for f in str(row.get("F_flags") or "").split(";")
                if f.strip() and f.strip() != "nan"
            ] if pd.notna(row.get("F_flags")) else [],
            "f_g": _num(row.get("F_g"), 2),
            "xy": [round(hub[0], 2), round(hub[1], 2)],
            "members": [],
        }
        elements[f"rxn:{fid}"] = el
        search_index.append({
            "key": f"rxn:{fid}",
            "kind": "Reaction",
            "title": f"{_short(fid)} ({fid})",
            "sub": f"{el['region_label']} · {el['flux_claim'] or el['draw_style']}",
            "tokens": f"{fid} {_short(fid)} {el['label']} {el['ec_number']} {el['kegg_reaction']} {eq}".lower(),
        })

    for _, row in aggs_df.iterrows():
        aid = str(row["element_id"])
        agg = layout.aggregates[aid]
        el = {
            "kind": "aggregate",
            "id": aid,
            "short": str(row.get("label") or aid),
            "label": str(row.get("label") or aid),
            "region": str(row.get("region") or ""),
            "region_label": layout.regions[str(row["region"])].label if str(row.get("region")) in layout.regions else "",
            "equation": "Multi-step / pathway module (omics-aggregated)",
            "draw_style": str(row.get("draw_style") or ""),
            "reversible": False,
            "physio_sign": 0,
            "sign_change": False,
            "flux_group": "",
            "kegg_pathway": str(row.get("kegg_pathway") or ""),
            "kegg_reaction": str(row.get("kegg_reaction") or ""),
            "ec_number": str(row.get("ec_number") or ""),
            "flux_claim": str(row.get("flux_claim") or ""),
            "rank": int(row["rank"]) if pd.notna(row.get("rank")) else None,
            "value": _num(row.get("P_value"), 3),
            "enz_dir": int(row["C_direction"]) if pd.notna(row.get("C_direction")) else 0,
            "enz_mag": _num(abs(float(row["P_value"])) if pd.notna(row.get("P_value")) else None, 3),
            "axes": {
                "F": {"state": str(row.get("F_state") or "not_measured"), "dir": 0},
                "C": {"state": str(row.get("C_state") or ""), "dir": int(row["C_direction"]) if pd.notna(row.get("C_direction")) else 0},
                "R": {"state": str(row.get("R_state") or ""), "dir": int(row["R_direction"]) if pd.notna(row.get("R_direction")) else 0},
                "M": {"state": str(row.get("M_state") or ""), "dir": int(row["M_direction"]) if pd.notna(row.get("M_direction")) else 0},
            },
            "layers": {
                "P": {"val": _num(row.get("P_value"), 3), "state": str(row.get("P_state") or "")},
                "R": {"val": _num(row.get("R_value"), 3), "state": str(row.get("R_layer_state") or "")},
                "M": {"val": _num(row.get("Mx_value"), 3), "state": str(row.get("Mx_state") or "")},
            },
            "xy": [round(agg.xy[0], 2), round(agg.xy[1], 2)],
            "members": [],
        }
        elements[f"agg:{aid}"] = el
        search_index.append({
            "key": f"agg:{aid}",
            "kind": "Module",
            "title": f"{el['label']} ({aid})",
            "sub": f"{el['region_label']} · {el['flux_claim'] or 'aggregate'}",
            "tokens": f"{aid} {el['label']} {el['region_label']}".lower(),
        })

    for sp, node in layout.nodes.items():
        if sp not in scene["active_species"] or node.kind != "metabolite":
            continue
        nrow = nodes_df.loc[sp] if sp in nodes_df.index else None
        xy = node.at[0]
        short_lbl = _species_short(sp, node.label)
        el = {
            "kind": "metabolite",
            "id": sp,
            "short": short_lbl,
            "label": str(nrow["label"]) if (nrow is not None and pd.notna(nrow.get("label"))) else (node.label or sp),
            "region": "",
            "region_label": "Metabolite Pool (whole-tissue measurement)",
            "equation": "",
            "draw_style": "metabolite",
            "value": _num(nrow.get("value"), 3) if nrow is not None else None,
            "display_state": str(nrow.get("display_state") or "not_measured") if nrow is not None else "not_measured",
            "ring": str(nrow.get("ring") or "") if nrow is not None else "",
            "n_changed": int(nrow["n_changed"]) if (nrow is not None and pd.notna(nrow.get("n_changed"))) else 0,
            "n_quantified": int(nrow["n_quantified"]) if (nrow is not None and pd.notna(nrow.get("n_quantified"))) else 0,
            "members_text": str(nrow.get("members") or "") if nrow is not None else "",
            "xy": [round(xy[0], 2), round(xy[1], 2)],
            "members": [],
        }
        elements[f"node:{sp}"] = el
        search_index.append({
            "key": f"node:{sp}",
            "kind": "Metabolite",
            "title": f"{short_lbl} ({sp})",
            "sub": el["label"],
            "tokens": f"{sp} {short_lbl} {el['label']} {el['members_text']}".lower(),
        })

    # 4. Member details & Per-animal samples (strictly filtered to on-map elements!)
    needed_sample_keys: set[tuple[str, str, str]] = set()
    num_label, den_label = "Numerator (Exp)", "Denominator (Ctrl)"

    if members is not None and not members.empty:
        m = members.copy()
        if tissue and "tissue" in m.columns:
            m = m[m["tissue"].astype(str) == str(tissue)]
        if contrast_id and "contrast_id" in m.columns:
            m = m[m["contrast_id"].astype(str) == str(contrast_id)]
        if source and "source" in m.columns:
            m = m[m["source"].astype(str) == str(source)]
        m_valid = m[m["cell_state"].astype(str) != "no_knowledge"].copy()

        # Attach reaction & aggregate members
        on_map_fids = set(reactions_df["element_id"].astype(str)) | set(aggs_df["element_id"].astype(str))
        m_rxn = m_valid[m_valid["flux_id"].astype(str).isin(on_map_fids)]
        for fid, grp in m_rxn.groupby("flux_id", sort=False):
            fid_s = str(fid)
            el_key = f"rxn:{fid_s}" if f"rxn:{fid_s}" in elements else (f"agg:{fid_s}" if f"agg:{fid_s}" in elements else "")
            if not el_key:
                continue
            extra_tokens: list[str] = []
            # Sort: MFA first (well-constrained/directional before insufficient_n/not_identifiable), then P/R/M/derived
            layer_ord = {"mfa": 0, "proteomics": 1, "transcriptomics": 2, "metabolomics": 3, "derived": 4}
            state_ord = {"changed": 0, "suggestive": 1, "unchanged": 2, "equivalent": 2, "undetermined": 3, "no_evidence": 3, "insufficient_n": 4, "not_identifiable": 5, "bound_pinned": 5, "not_detected": 6, "not_measured": 7}
            ident_ord = {"well_constrained": 0, "moderate": 1, "poorly_constrained": 2, "": 3}
            grp = grp.assign(
                _lo=grp["layer"].map(layer_ord).fillna(9),
                _so=grp["cell_state"].astype(str).map(state_ord).fillna(4),
                _io=grp["identifiability"].fillna("").astype(str).map(ident_ord).fillna(3),
                _sh=pd.to_numeric(grp.get("abundance_share"), errors="coerce").fillna(-1.0),
                _fc=pd.to_numeric(grp.get("log2fc"), errors="coerce").abs().fillna(-1.0),
            ).sort_values(["_lo", "_so", "_io", "_sh", "_fc"], ascending=[True, True, True, False, False], kind="stable")

            mem_records: list[dict[str, Any]] = []
            for r in grp.itertuples():
                layer = str(getattr(r, "layer", "") or "")
                tracer = str(getattr(r, "tracer", "") or "")
                feat_id = str(getattr(r, "feature_id", "") or "")
                mem_id = str(getattr(r, "member_id", "") or "")
                feat_name = str(getattr(r, "feature_name", "") or "")
                disp_name = format_tracer_display(tracer) if layer == "mfa" else (mem_id or feat_name or feat_id)
                skey = f"{layer}|{feat_id}|{tracer if layer == 'mfa' else ''}" if feat_id else ""
                if feat_id:
                    needed_sample_keys.add((layer, feat_id, tracer if layer == "mfa" else ""))
                if disp_name:
                    extra_tokens.append(disp_name)
                if feat_name and feat_name != disp_name:
                    extra_tokens.append(feat_name)
                mem_records.append({
                    "layer": layer,
                    "role": str(getattr(r, "role", "") or ""),
                    "name": disp_name,
                    "feature_id": feat_id,
                    "feature_name": feat_name,
                    "tracer": tracer,
                    "relation": str(getattr(r, "member_relation", "") or ""),
                    "reg_effect": str(getattr(r, "regulator_effect", "") or ""),
                    "log2fc": _num(getattr(r, "log2fc", None), 3),
                    "effect": _num(getattr(r, "effect", None), 2),
                    "effect_type": str(getattr(r, "effect_type", "") or "") if pd.notna(getattr(r, "effect_type", None)) else "",
                    "smd_se": _num(getattr(r, "smd_se", None), 3),
                    "padj": _num(getattr(r, "padj", None), 5),
                    "padj_global": _num(getattr(r, "padj_global", None), 5),
                    "fdr_scope": str(getattr(r, "fdr_scope", "") or "") if pd.notna(getattr(r, "fdr_scope", None)) else "",
                    "pvalue": _num(getattr(r, "pvalue", None), 5),
                    "lfsr": _num(getattr(r, "lfsr", None), 4),
                    "share": _num(getattr(r, "abundance_share", None), 3),
                    "mean_num": _num(getattr(r, "mean_numerator", None), 3),
                    "mean_den": _num(getattr(r, "mean_denominator", None), 3),
                    "n_num": int(r.n_numerator) if pd.notna(getattr(r, "n_numerator", np.nan)) else None,
                    "n_den": int(r.n_denominator) if pd.notna(getattr(r, "n_denominator", np.nan)) else None,
                    "state": str(getattr(r, "cell_state", "") or ""),
                    "reason": str(getattr(r, "state_reason", "") or ""),
                    "ident": str(getattr(r, "identifiability", "") or ""),
                    "het_p": _num(getattr(r, "heterogeneity_pvalue", None), 4),
                    "dir_stab": _num(getattr(r, "direction_stability", None), 3),
                    "loo_stab": _num(getattr(r, "pooled_loo_stability", None), 3),
                    "smd_weights": str(getattr(r, "smd_weights", "") or "") if pd.notna(getattr(r, "smd_weights", None)) else "",
                    "tracer_weights": str(getattr(r, "tracer_weights", "") or "") if pd.notna(getattr(r, "tracer_weights", None)) else "",
                    "g_total_var": _num(getattr(r, "g_total_var", None), 3),
                    "smd_weighting": str(getattr(r, "smd_weighting", "") or "") if pd.notna(getattr(r, "smd_weighting", None)) else "",
                    "flags": str(getattr(r, "flags", "") or "") if pd.notna(getattr(r, "flags", None)) else "",
                    "skey": skey,
                })
            elements[el_key]["members"] = mem_records
            if extra_tokens:
                for item in search_index:
                    if item["key"] == el_key:
                        item["tokens"] += " " + " ".join(dict.fromkeys(extra_tokens)).lower()
                        break

        # Attach metabolite node members (so clicking a metabolite circle also shows per-animal dots!)
        try:
            from ..dossier.view import node_members as _node_members
        except Exception:
            _node_members = None
        for sp in list(scene["active_species"]):
            el_key = f"node:{sp}"
            if el_key not in elements:
                continue
            try:
                if _node_members is not None:
                    nm = _node_members(m, sp)
                else:
                    sp_base = sp.rsplit("_", 1)[0].lower()
                    nm = m[
                        m["layer"].astype(str).isin(("metabolomics", "derived"))
                        & (
                            m["member_id"].fillna("").astype(str).str.lower().eq(sp_base)
                            | m["feature_name"].fillna("").astype(str).str.lower().eq(sp_base)
                        )
                    ]
            except Exception:
                continue
            if nm.empty:
                continue
            nm = nm[nm["cell_state"].astype(str) != "no_knowledge"].drop_duplicates(["layer", "feature_id"])
            m_recs: list[dict[str, Any]] = []
            for r in nm.itertuples():
                layer = str(getattr(r, "layer", "") or "")
                feat_id = str(getattr(r, "feature_id", "") or "")
                feat_name = str(getattr(r, "feature_name", "") or "") or str(getattr(r, "member_id", "") or "") or feat_id
                skey = f"{layer}|{feat_id}|" if feat_id else ""
                if feat_id:
                    needed_sample_keys.add((layer, feat_id, ""))
                m_recs.append({
                    "layer": layer,
                    "role": "pool",
                    "name": feat_name,
                    "feature_id": feat_id,
                    "feature_name": feat_name,
                    "tracer": "",
                    "relation": "",
                    "reg_effect": "",
                    "log2fc": _num(getattr(r, "log2fc", None), 3),
                    "padj": _num(getattr(r, "padj", None), 5),
                    "pvalue": _num(getattr(r, "pvalue", None), 5),
                    "lfsr": _num(getattr(r, "lfsr", None), 4),
                    "share": None,
                    "mean_num": _num(getattr(r, "mean_numerator", None), 3),
                    "mean_den": _num(getattr(r, "mean_denominator", None), 3),
                    "n_num": int(r.n_numerator) if pd.notna(getattr(r, "n_numerator", np.nan)) else None,
                    "n_den": int(r.n_denominator) if pd.notna(getattr(r, "n_denominator", np.nan)) else None,
                    "state": str(getattr(r, "cell_state", "") or ""),
                    "reason": str(getattr(r, "state_reason", "") or ""),
                    "ident": "",
                    "skey": skey,
                })
            elements[el_key]["members"] = m_recs

    # 5. Per-animal samples extraction (compact tuple layout)
    sample_series: dict[str, dict[str, Any]] = {}
    side_map: dict[tuple[str, str], str] = {}
    if sample_contrasts is not None and not sample_contrasts.empty:
        sc = sample_contrasts.copy()
        if contrast_id and "contrast_id" in sc.columns:
            sc = sc[sc["contrast_id"].astype(str) == str(contrast_id)]
        if source and "source" in sc.columns:
            sc = sc[sc["source"].astype(str) == str(source)]
        for r in sc.itertuples():
            s_val = str(getattr(r, "side", "") or "")
            side_map[(str(getattr(r, "layer", "")), str(getattr(r, "sample_id", "")))] = s_val
            lvl = str(getattr(r, "level", "") or "")
            if lvl:
                if s_val == "numerator":
                    num_label = lvl
                elif s_val == "denominator":
                    den_label = lvl

    if samples is not None and not samples.empty and needed_sample_keys:
        sp_df = samples.copy()
        if tissue and "tissue" in sp_df.columns:
            sp_df = sp_df[sp_df["tissue"].astype(str) == str(tissue)]
        if source and "source" in sp_df.columns:
            sp_df = sp_df[sp_df["source"].astype(str) == str(source)]
        needed_fids = {fid for _, fid, _ in needed_sample_keys}
        sp_df = sp_df[sp_df["feature_id"].astype(str).isin(needed_fids)]
        if not sp_df.empty:
            sp_df = sp_df.assign(
                _tracer_clean=sp_df["tracer"].fillna("").astype(str) if "tracer" in sp_df.columns else "",
                _go=pd.to_numeric(sp_df["group_order"], errors="coerce").fillna(999).astype(int) if "group_order" in sp_df.columns else 0,
            )
            for (layer, feat_id, tr_req) in sorted(needed_sample_keys):
                part = sp_df[(sp_df["layer"].astype(str) == layer) & (sp_df["feature_id"].astype(str) == feat_id)]
                if layer == "mfa" and tr_req and tr_req != "pooled":
                    part = part[part["_tracer_clean"] == tr_req]
                part = part[pd.to_numeric(part["value"], errors="coerce").notna()]
                if part.empty:
                    continue
                groups_df = part[["group", "_go"]].drop_duplicates().sort_values(["_go", "group"])
                groups_list = [str(g) for g in groups_df["group"]]
                g_idx = {g: i for i, g in enumerate(groups_list)}
                pts: list[list[Any]] = []
                for r in part.sort_values(["_go", "sample_id"]).itertuples():
                    sid = str(getattr(r, "sample_id", "") or "")
                    subj = str(getattr(r, "subject_id", "") or sid)
                    grp_name = str(getattr(r, "group", "") or "")
                    side_full = side_map.get((layer, sid), "other")
                    side_code = "n" if side_full == "numerator" else ("d" if side_full == "denominator" else "o")
                    at_b = 1 if (pd.notna(getattr(r, "at_bound", None)) and bool(getattr(r, "at_bound", False))) else 0
                    pts.append([
                        subj,
                        sid,
                        g_idx.get(grp_name, 0),
                        side_code,
                        _num(getattr(r, "value", None), 4),
                        _num(getattr(r, "fit_low", None), 4),
                        _num(getattr(r, "fit_high", None), 4),
                        at_b,
                        str(getattr(r, "qc_flag", "") or ""),
                    ])
                scale = str(part["scale"].iloc[0]) if "scale" in part.columns else "log2"
                units = str(part["units"].iloc[0]) if ("units" in part.columns and str(part["units"].iloc[0]) not in ("", "TBD", "nan")) else ""
                sample_series[f"{layer}|{feat_id}|{tr_req if layer == 'mfa' else ''}"] = {
                    "scale": scale,
                    "units": units,
                    "groups": groups_list,
                    "pts": pts,
                }

    # 6. Top ranked reactions for overview card
    top_rxns: list[dict[str, Any]] = []
    ranked = reactions_df[reactions_df["rank"].notna()].sort_values("rank").head(12)
    for _, r in ranked.iterrows():
        fid = str(r["element_id"])
        top_rxns.append({
            "key": f"rxn:{fid}",
            "rank": int(r["rank"]),
            "id": fid,
            "short": _short(fid),
            "claim": str(r.get("flux_claim") or ""),
            "value": _num(r.get("value"), 2),
            "region": layout.regions[str(r["region"])].label if str(r.get("region")) in layout.regions else str(r.get("region") or ""),
        })

    # Extra fusion pathways (top 6)
    extra_list: list[dict[str, Any]] = []
    if extra_fusion is not None and not extra_fusion.empty:
        for r in extra_fusion.head(6).itertuples():
            extra_list.append({
                "term": str(getattr(r, "pathway_term", "")),
                "fdr": _num(getattr(r, "fused_fdr", None), 4),
                "dirs": str(getattr(r, "omics_directions", "") or ""),
            })

    claims = view.loc[view["flux_claim"].isin(("L1", "L2", "L3")), "flux_claim"].value_counts()
    indep_chains = _independent_flux_chains(view)
    return {
        "meta": {
            "title": title,
            "tracer": str(view["tracer"].iloc[0]) if "tracer" in view and not view.empty else "",
            "tissue": tissue,
            "contrast": contrast_id,
            "caption": list(caption),
            "initial_style": style if style in ("metro_atlas", "economist", "nature_editorial", "slide_dark") else "metro_atlas",
            "color_mode": color_mode,
            "strip_mode": strip_mode,
            "omics_edge_color": omics_edge_color,
            "num_label": num_label,
            "den_label": den_label,
            "has_samples": bool(sample_series),
            "counts": {
                "L3": int(claims.get("L3", 0)),
                "L2": int(claims.get("L2", 0)),
                "L1": int(claims.get("L1", 0)),
                "indep_flux_chains": int(indep_chains),
                "linked_groups": len(scene["multi_groups"]),
            },
            "vb": [round(v, 2) for v in vb],
        },
        "palettes": {k: _palette_dict(k) for k in ("metro_atlas", "economist", "nature_editorial", "slide_dark")},
        "submaps": submap_views,
        "regions": region_list,
        "top_reactions": top_rxns,
        "extra_fusion": extra_list,
        "elements": elements,
        "samples": sample_series,
        "search": search_index,
    }


# ---------------------------------------------------------------------------
# Offline Self-Contained Single-File HTML Template (Metro Atlas Default + Theme Switcher)
# ---------------------------------------------------------------------------

_CSS = r"""
:root {
  --surface: #FFFFFF;
  --card: #F8FAFC;
  --card-edge: #DCE3EC;
  --ink: #0F172A;
  --ink2: #334155;
  --muted: #64748B;
  --hairline: #E2E8F0;
  --up: #C81E1E;
  --down: #1D4ED8;
  --neutral: #94A3B8;
  --coupled-bg: #EEF4FF;
  --coupled-edge: #3B82F6;
  --coupled-ink: #1D4ED8;
  --red-tag: #C81E1E;
  --radius: 5px;
  --font-sans: "Arial", "Helvetica Neue", Helvetica, "Liberation Sans", "Microsoft YaHei", "PingFang SC", sans-serif;
  --font-mono: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
}
* { box-sizing: border-box; }
html, body {
  margin: 0; padding: 0; height: 100%; width: 100%;
  background: var(--surface); color: var(--ink);
  font-family: var(--font-sans); font-size: 13px; line-height: 1.42;
  overflow: hidden;
}
#app { display: flex; flex-direction: column; height: 100vh; width: 100vw; }

/* Masthead: delicate academic hairline by default; Economist red tag & heavy rule when active */
.masthead {
  flex: 0 0 auto; background: var(--surface);
  border-bottom: 1px solid var(--card-edge);
  padding: 0 18px 8px 18px; position: relative;
}
.red-tag-bar { height: 5px; border-top: none; position: relative; margin-bottom: 3px; }
.red-tag {
  position: absolute; top: -1px; left: 0; width: 44px; height: 8px;
  background: var(--red-tag); display: none;
}
body[data-theme="economist"] .masthead { border-bottom: 1.5px solid var(--ink); }
body[data-theme="economist"] .red-tag-bar { height: 10px; border-top: 1px solid var(--card-edge); margin-bottom: 5px; }
body[data-theme="economist"] .red-tag { display: block; }

.mast-row { display: flex; align-items: center; justify-content: space-between; gap: 12px; flex-wrap: wrap; }
.mast-title-wrap { min-width: 240px; }
.mast-title { margin: 0; font-size: 16.5px; font-weight: 800; letter-spacing: -0.015em; color: var(--ink); }
.mast-sub { font-size: 11.5px; color: var(--muted); margin-top: 2px; display: flex; gap: 12px; flex-wrap: wrap; }
.mast-pills { display: flex; align-items: center; gap: 6px; flex-wrap: wrap; }
.stat-pill {
  display: inline-flex; align-items: center; gap: 5px;
  padding: 2px 8px; font-size: 11px; font-weight: 700;
  border-radius: var(--radius); border: 1px solid var(--card-edge); background: var(--card); color: var(--ink);
}
.stat-pill.l3 { background: var(--ink); color: var(--surface); border-color: var(--ink); }
.stat-pill.l2 { background: #475569; color: #FFFFFF; border-color: #475569; }
.stat-pill.l1 { border-style: dashed; border-color: var(--ink2); color: var(--ink); background: var(--surface); }
.stat-pill.link { background: var(--coupled-bg); border-color: var(--coupled-edge); color: var(--coupled-ink); }

/* Control Toolbar */
.toolbar {
  flex: 0 0 auto; display: flex; align-items: center; gap: 6px; flex-wrap: wrap;
  padding: 5px 16px; background: var(--card); border-bottom: 1px solid var(--card-edge);
}
.tb-group { display: inline-flex; align-items: center; gap: 3px; }
.tb-sep { width: 1px; height: 16px; background: var(--card-edge); margin: 0 2px; }
.tb-label { font-size: 10px; font-weight: 700; text-transform: uppercase; letter-spacing: 0.04em; color: var(--muted); margin-right: 2px; }
.btn {
  appearance: none; background: var(--surface); color: var(--ink);
  border: 1px solid var(--card-edge); border-radius: var(--radius);
  padding: 2.5px 7px; font-size: 11px; font-weight: 600; cursor: pointer;
  transition: background 0.12s, border-color 0.12s, color 0.12s;
}
.btn:hover { border-color: var(--ink2); }
.btn.active { background: var(--ink); color: var(--surface); border-color: var(--ink); }
.btn.red-accent.active { background: var(--up); color: #fff; border-color: var(--up); }
select.ctrl-select {
  background: var(--surface); color: var(--ink); border: 1px solid var(--card-edge);
  border-radius: var(--radius); padding: 2.5px 6px; font-size: 11px; font-weight: 600; cursor: pointer;
}
.search-wrap { position: relative; min-width: 175px; }
.search-input {
  width: 100%; background: var(--surface); color: var(--ink);
  border: 1px solid var(--card-edge); border-radius: var(--radius);
  padding: 2.5px 22px 2.5px 7px; font-size: 11px; outline: none;
}
.search-input:focus { border-color: var(--ink); }
.search-results {
  display: none; position: absolute; top: 100%; left: 0; right: 0; z-index: 40;
  max-height: 290px; overflow-y: auto; background: var(--card);
  border: 1px solid var(--card-edge); box-shadow: 0 6px 18px rgba(15, 23, 42, 0.10);
  border-radius: var(--radius); margin-top: 2px;
}
.search-results.open { display: block; }
.sr-item {
  padding: 6px 9px; cursor: pointer; border-bottom: 1px solid var(--hairline);
  display: flex; justify-content: space-between; align-items: baseline; gap: 6px;
}
.sr-item:hover { background: var(--surface); }
.sr-title { font-weight: 700; font-size: 11.5px; color: var(--ink); }
.sr-sub { font-size: 10.5px; color: var(--muted); white-space: nowrap; overflow: hidden; text-overflow: ellipsis; max-width: 125px; }

/* Workspace Split */
.workspace { flex: 1 1 auto; display: flex; min-height: 0; position: relative; overflow: hidden; }
.stage { flex: 1 1 auto; position: relative; overflow: hidden; background: var(--surface); cursor: grab; user-select: none; }
.stage.panning { cursor: grabbing; }
#ow-svg { width: 100%; height: 100%; display: block; }

/* SVG Interactivity (delicate pt-calibrated highlight strokes) */
.ow-edge, .ow-node, .ow-agg, .ow-strip { cursor: pointer; transition: opacity 0.14s ease; }
.ow-dimmed { opacity: 0.14 !important; }
.ow-selected .ow-card-bg { stroke: var(--up) !important; stroke-width: 0.24 !important; }
.ow-selected .ow-node-dot { stroke: var(--up) !important; stroke-width: 0.28 !important; }
.ow-selected .ow-edge-main { filter: drop-shadow(0 0 0.30px var(--up)); }
.ow-group-pulse .ow-conduit-fill { opacity: 1 !important; stroke-width: 1.50 !important; }
.ow-group-pulse .ow-conduit-edge { opacity: 0.72 !important; stroke-width: 1.85 !important; }

/* Floating Zoom & Legend Overlay */
.zoom-ctrls {
  position: absolute; left: 14px; bottom: 14px; z-index: 15;
  display: flex; flex-direction: column; gap: 3px;
}
.zoom-ctrls .btn { width: 30px; height: 28px; padding: 0; font-size: 14px; display: flex; align-items: center; justify-content: center; background: var(--card); }
.legend-float {
  position: absolute; left: 54px; bottom: 14px; z-index: 16; max-width: 660px;
  background: var(--card); border: 1px solid var(--card-edge); border-top: 2px solid var(--ink2);
  border-radius: var(--radius); padding: 9px 12px; font-size: 11px; color: var(--ink2);
  box-shadow: 0 6px 18px rgba(15, 23, 42, 0.08); display: none;
}
.legend-float.open { display: block; }
.legend-grid { display: grid; grid-template-columns: repeat(2, minmax(220px, 1fr)); gap: 6px 14px; }
.lg-row { display: flex; align-items: center; gap: 7px; }
.lg-swatch { flex: 0 0 32px; height: 10px; display: inline-flex; align-items: center; justify-content: center; }

/* Hover Tooltip */
.tooltip {
  position: fixed; z-index: 60; pointer-events: none; display: none; max-width: 310px;
  background: var(--surface); color: var(--ink);
  border: 1px solid var(--card-edge); border-top: 2.5px solid var(--up);
  box-shadow: 0 8px 20px rgba(15, 23, 42, 0.12);
  border-radius: var(--radius); padding: 8px 10px; font-size: 11.5px;
}
body[data-theme="economist"] .tooltip { background: var(--card); border-color: var(--ink); border-top: 3px solid var(--red-tag); box-shadow: none; }
.tt-head { display: flex; justify-content: space-between; align-items: center; gap: 8px; margin-bottom: 3px; }
.tt-title { font-weight: 800; font-size: 12.5px; }
.tt-eq { font-size: 11px; color: var(--ink2); margin-bottom: 4px; }
.tt-meta { font-size: 10.5px; color: var(--muted); }

/* Right Dossier Drawer */
.drawer {
  flex: 0 0 445px; width: 445px; max-width: 46vw;
  background: var(--card); border-left: 1px solid var(--card-edge);
  display: flex; flex-direction: column; min-height: 0; overflow: hidden;
  transition: width 0.16s ease, flex-basis 0.16s ease;
}
body[data-theme="economist"] .drawer { border-left: 1.5px solid var(--ink); }
.drawer.collapsed { flex-basis: 0; width: 0; border-left: none; }
.drawer-head {
  flex: 0 0 auto; padding: 11px 14px 9px 14px;
  border-bottom: 1px solid var(--card-edge); background: var(--surface);
}
.dh-top { display: flex; justify-content: space-between; align-items: flex-start; gap: 8px; }
.dh-title { margin: 0; font-size: 15px; font-weight: 800; color: var(--ink); line-height: 1.25; }
.dh-sub { font-size: 11.5px; color: var(--ink2); margin-top: 3px; }
.dh-tags { display: flex; gap: 5px; flex-wrap: wrap; margin-top: 6px; }
.tag {
  font-size: 10.5px; font-weight: 700; padding: 1px 6px; border-radius: var(--radius);
  border: 1px solid var(--card-edge); background: var(--surface); color: var(--ink2);
}
.tag.up { background: var(--up); color: #fff; border-color: var(--up); }
.tag.down { background: var(--down); color: #fff; border-color: var(--down); }
.tag.linked { background: var(--coupled-bg); color: var(--coupled-ink); border-color: var(--coupled-edge); cursor: pointer; }
.drawer-body { flex: 1 1 auto; overflow-y: auto; padding: 12px 14px 18px 14px; }
.sec-title {
  font-size: 10.5px; font-weight: 800; text-transform: uppercase; letter-spacing: 0.045em;
  color: var(--ink2); border-top: 1px solid var(--card-edge); padding-top: 7px; margin: 12px 0 7px 0;
  display: flex; justify-content: space-between; align-items: center;
}
body[data-theme="economist"] .sec-title { color: var(--ink); border-top: 1.5px solid var(--ink); }
.sec-title:first-child { margin-top: 0; border-top: none; padding-top: 0; }

/* 4-Omics Axis Grid */
.axis-grid { display: grid; grid-template-columns: repeat(4, 1fr); gap: 6px; margin-bottom: 8px; }
.axis-card {
  background: var(--surface); border: 1px solid var(--card-edge); border-radius: var(--radius);
  padding: 6px 7px; text-align: center;
}
.ac-lbl { font-size: 10px; font-weight: 700; color: var(--muted); text-transform: uppercase; }
.ac-val { font-size: 12.5px; font-weight: 800; margin-top: 2px; }
.ac-st { font-size: 10px; color: var(--ink2); margin-top: 1px; }

/* Flux Assessment Gate Checklist */
.gate-box {
  background: var(--surface); border: 1px solid var(--card-edge);
  border-radius: var(--radius); padding: 7px 8px; margin-bottom: 8px;
}
body[data-theme="economist"] .gate-box { border-top: 1.5px solid var(--ink); }
.gate-hdr {
  display: flex; justify-content: space-between; align-items: center;
  font-size: 10.5px; font-weight: 800; color: var(--ink); margin-bottom: 5px;
}
.gate-row {
  display: flex; align-items: baseline; justify-content: space-between; gap: 6px;
  font-size: 10.5px; padding: 2.5px 0; border-top: 1px solid var(--hairline);
}
.gate-pill {
  display: inline-block; font-size: 9.5px; font-weight: 700; padding: 0.5px 5px;
  border-radius: var(--radius); white-space: nowrap; flex-shrink: 0;
}
.gate-pass { background: #E6F4EA; color: #137333; }
.gate-fail { background: #FCE8E6; color: #A50E0E; }
.gate-warn { background: #FEF7E0; color: #9A6700; }
.gate-na { background: var(--hairline); color: var(--muted); }
body[data-theme="slide_dark"] .gate-pass { background: #163826; color: #81C995; }
body[data-theme="slide_dark"] .gate-fail { background: #431B1B; color: #F28B82; }
body[data-theme="slide_dark"] .gate-warn { background: #3E2F10; color: #FDD663; }
.gate-meta { font-size: 10px; color: var(--muted); margin-top: 4px; line-height: 1.4; }

/* Per-Animal Scatter Small Multiples */
.layer-tabs { display: flex; gap: 4px; flex-wrap: wrap; margin-bottom: 8px; }
.scatter-grid { display: grid; grid-template-columns: repeat(2, 1fr); gap: 8px; }
.scatter-grid.single-col { grid-template-columns: 1fr; }
.sc-card {
  background: var(--surface); border: 1px solid var(--card-edge);
  border-radius: var(--radius); padding: 6px 7px 5px 7px;
}
body[data-theme="economist"] .sc-card { border-top: 1.5px solid var(--ink); }
.sc-head { font-size: 11px; font-weight: 800; color: var(--ink); white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
.sc-sub { font-size: 10px; color: var(--muted); margin-bottom: 3px; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
.sc-svg { width: 100%; height: 112px; display: block; overflow: visible; }
.sc-legend {
  display: flex; align-items: center; gap: 10px; flex-wrap: wrap;
  font-size: 10.5px; color: var(--ink2); margin-bottom: 7px;
}
.sc-dot { width: 8px; height: 8px; border-radius: 50%; display: inline-block; margin-right: 3px; vertical-align: middle; }

/* Member Table */
table.mem-tbl { width: 100%; border-collapse: collapse; font-size: 11px; }
table.mem-tbl th {
  text-align: left; font-size: 10px; font-weight: 700; color: var(--muted);
  border-bottom: 1px solid var(--card-edge); padding: 4px 4px; cursor: pointer;
}
body[data-theme="economist"] table.mem-tbl th { border-bottom: 1.5px solid var(--ink); }
table.mem-tbl td { border-bottom: 1px solid var(--hairline); padding: 4px 4px; vertical-align: middle; }
table.mem-tbl tr:hover td { background: var(--surface); }
.fc-bar-wrap { display: inline-flex; align-items: center; gap: 4px; }
.fc-bar { height: 7px; display: inline-block; min-width: 2px; }

/* Overview Coverage Bars */
.reg-row { padding: 6px 0; border-bottom: 1px solid var(--hairline); cursor: pointer; }
.reg-row:hover { background: var(--surface); margin: 0 -6px; padding: 6px 6px; }
.reg-top { display: flex; justify-content: space-between; font-size: 11.5px; font-weight: 700; }
.reg-bar { height: 6px; display: flex; background: var(--hairline); margin: 4px 0 3px 0; border-radius: 2px; overflow: hidden; }
.reg-stat { font-size: 10.5px; color: var(--muted); display: flex; justify-content: space-between; }
.rank-item {
  display: flex; align-items: center; justify-content: space-between; gap: 6px;
  padding: 5px 4px; border-bottom: 1px solid var(--hairline); cursor: pointer; font-size: 11.5px;
}
.rank-item:hover { background: var(--surface); }
.rk-badge {
  display: inline-block; min-width: 22px; text-align: center; font-size: 10px; font-weight: 800;
  padding: 1px 4px; background: var(--ink); color: var(--surface); border-radius: var(--radius);
}
.rk-badge.top3 { background: var(--up); color: #fff; }

/* Provenance Footer */
.footer {
  flex: 0 0 auto; background: var(--surface); border-top: 1px solid var(--card-edge);
  padding: 4px 18px; font-size: 10.5px; color: var(--muted);
  display: flex; justify-content: space-between; gap: 12px; flex-wrap: wrap;
}
"""

_JS = r"""
(function() {
  "use strict";
  const DATA = JSON.parse(document.getElementById("ow-payload").textContent);
  const svg = document.getElementById("ow-svg");
  const stage = document.getElementById("ow-stage");
  const tooltip = document.getElementById("ow-tooltip");
  const drawer = document.getElementById("ow-drawer");
  const drawerHead = document.getElementById("ow-drawer-head");
  const drawerBody = document.getElementById("ow-drawer-body");

  let state = {
    theme: DATA.meta.initial_style || "metro_atlas",
    omicsColor: DATA.meta.omics_edge_color || "tinted",
    filter: "all",
    selectedKey: null,
    activeScatterLayer: "key",
    vb: [...DATA.meta.vb],
    baseVb: [...DATA.meta.vb],
    dragging: false,
    dragStart: null,
    dragVb: null,
    memSort: { col: "abs_fc", desc: true },
  };

  /* ---------- Color Helpers (Exact Semantic Diverging Ramp) ---------- */
  function hexToRgb(h) {
    const c = (h || "#888888").replace("#", "");
    const n = parseInt(c.length === 3 ? c.split("").map(x => x + x).join("") : c, 16);
    return [(n >> 16) & 255, (n >> 8) & 255, n & 255];
  }
  function rgbToHex(r, g, b) {
    return "#" + [r, g, b].map(v => Math.max(0, Math.min(255, Math.round(v))).toString(16).padStart(2, "0")).join("");
  }
  function mixHex(c1, c2, t) {
    const a = hexToRgb(c1), b = hexToRgb(c2);
    return rgbToHex(a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t);
  }
  function divergingColor(val, pal) {
    if (val === null || val === undefined || isNaN(val)) return pal.card_bg;
    const stops = pal.diverging_stops;
    const t = Math.max(0, Math.min(1, (val + 2.0) / 4.0));
    const scaled = t * (stops.length - 1);
    const idx = Math.min(stops.length - 2, Math.floor(scaled));
    return mixHex(stops[idx], stops[idx + 1], scaled - idx);
  }
  function enzTintColor(dir, mag, pal) {
    const base = dir > 0 ? pal.up : pal.down;
    const frac = 0.56 + 0.16 * Math.min(Math.abs(mag || 1.0) / 2.0, 1.0);
    return mixHex(pal.surface, base, frac);
  }
  function fmtNum(v, d) {
    if (v === null || v === undefined || isNaN(v)) return "—";
    return Number(v).toFixed(d === undefined ? 2 : d);
  }
  function fmtSigned(v, d) {
    if (v === null || v === undefined || isNaN(v)) return "—";
    const n = Number(v);
    return (n >= 0 ? "+" : "") + n.toFixed(d === undefined ? 2 : d);
  }
  function fmtP(v) {
    if (v === null || v === undefined || isNaN(v)) return "—";
    const n = Number(v);
    return n < 0.001 ? n.toExponential(1) : n.toFixed(3);
  }
  function esc(s) {
    return String(s == null ? "" : s).replace(/[&<>"']/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
  }

  /* ---------- Cross-Platform Zero-Overflow SVG Label Auto-Fit ---------- */
  function fitSvgLabels() {
    try {
      // 1. Ensure tier & coupled badge text sits comfortably inside its pill
      svg.querySelectorAll(".ow-tier-badge, .ow-coupled-badge").forEach(b => {
        const rect = b.querySelector(".ow-tb-bg, .ow-cb-bg");
        const txt = b.querySelector(".ow-tb-txt, .ow-cb-txt");
        if (!rect || !txt) return;
        const rb = rect.getBBox();
        const tb = txt.getBBox();
        const padR = 0.18;
        const maxRight = rb.x + rb.width - padR;
        if (tb.width > 0 && tb.x + tb.width > maxRight) {
          const avail = maxRight - tb.x;
          const curFs = parseFloat(txt.getAttribute("font-size")) || 0.40;
          const ratio = Math.max(0.80, avail / tb.width);
          txt.setAttribute("font-size", (curFs * ratio).toFixed(3));
          const tb2 = txt.getBBox();
          if (tb2.x + tb2.width > maxRight) {
            const needW = (tb2.x + tb2.width + padR) - rb.x;
            rect.setAttribute("width", needW.toFixed(2));
          }
        }
      });

      // 2. Ensure reaction strip labels never collide with badges or overflow card right edge
      svg.querySelectorAll(".ow-strip:not(.ow-strip-compact)").forEach(g => {
        const card = g.querySelector(".ow-card-bg");
        const lbl = g.querySelector(".ow-rxn-lbl");
        if (!card || !lbl) return;
        const lb = lbl.getBBox();
        let cursorX = lb.x + lb.width + 0.36;
        const badges = g.querySelectorAll(".ow-tier-badge, .ow-coupled-badge");
        badges.forEach(b => {
          b.removeAttribute("transform");
          const rect = b.querySelector(".ow-tb-bg, .ow-cb-bg");
          if (!rect) return;
          const rx = parseFloat(rect.getAttribute("x")) || 0;
          const rw = parseFloat(rect.getAttribute("width")) || 0;
          if (rx < cursorX) {
            const dx = cursorX - rx;
            b.setAttribute("transform", `translate(${dx.toFixed(2)},0)`);
            cursorX = cursorX + rw + 0.30;
          } else {
            cursorX = rx + rw + 0.30;
          }
        });
        const cardX = parseFloat(card.getAttribute("x")) || 0;
        const cardW = parseFloat(card.getAttribute("width")) || 0;
        const rightContent = badges.length ? (cursorX + 0.08) : (lb.x + lb.width + 0.38);
        if (rightContent > cardX + cardW) {
          card.setAttribute("width", (rightContent - cardX).toFixed(2));
        }
      });

      // 3. Ensure aggregate titles never collide with their tier badge or box edges
      svg.querySelectorAll(".ow-agg").forEach(g => {
        const title = g.querySelector(".ow-card-title");
        if (!title) return;
        const tb = title.getBBox();
        if (tb.width <= 0) return;
        const curFs = parseFloat(title.getAttribute("font-size")) || 0.44;
        const maxXAttr = title.getAttribute("data-max-x");
        const maxWAttr = title.getAttribute("data-max-w");
        if (maxXAttr) {
          const maxX = parseFloat(maxXAttr);
          if (tb.x + tb.width > maxX) {
            const ratio = Math.max(0.68, (maxX - tb.x) / tb.width);
            title.setAttribute("font-size", (curFs * ratio).toFixed(3));
          }
        } else if (maxWAttr) {
          const maxW = parseFloat(maxWAttr);
          if (tb.width > maxW) {
            const ratio = Math.max(0.68, maxW / tb.width);
            title.setAttribute("font-size", (curFs * ratio).toFixed(3));
          }
        }
      });

      // 4. Ensure boundary & store pills always enclose their labels with symmetric padding
      svg.querySelectorAll(".ow-boundary, .ow-store").forEach(g => {
        const rect = g.querySelector(".ow-pill-bg");
        const txt = g.querySelector(".ow-pill-txt");
        if (!rect || !txt) return;
        const tb = txt.getBBox();
        const curW = parseFloat(rect.getAttribute("width")) || 0;
        const needW = tb.width + 0.96;
        if (needW > curW) {
          const cx = (parseFloat(rect.getAttribute("x")) || 0) + curW / 2;
          rect.setAttribute("x", (cx - needW / 2).toFixed(2));
          rect.setAttribute("width", needW.toFixed(2));
        }
      });
    } catch (_) {
      /* getBBox unavailable in non-visual environments */
    }
  }

  /* ---------- Live Theme & Edge-Style Re-Rendering ---------- */
  function applyTheme(themeName) {
    state.theme = themeName;
    const pal = DATA.palettes[themeName] || DATA.palettes.metro_atlas;
    const isEcon = themeName === "economist";
    const cardSurface = pal.strip_card_bg || pal.surface;
    document.body.setAttribute("data-theme", themeName);
    const r = document.documentElement.style;
    r.setProperty("--surface", pal.surface);
    r.setProperty("--card", pal.card_bg);
    r.setProperty("--card-edge", pal.card_edge);
    r.setProperty("--ink", pal.ink);
    r.setProperty("--ink2", pal.ink2);
    r.setProperty("--muted", pal.muted);
    r.setProperty("--hairline", pal.hairline);
    r.setProperty("--up", pal.up);
    r.setProperty("--down", pal.down);
    r.setProperty("--neutral", pal.neutral);
    r.setProperty("--coupled-bg", pal.coupled_fill);
    r.setProperty("--coupled-edge", pal.coupled_edge);
    r.setProperty("--coupled-ink", pal.coupled_ink);
    r.setProperty("--red-tag", pal.red_tag);
    r.setProperty("--radius", isEcon ? "2px" : "5px");

    const bgRect = document.getElementById("ow-canvas-bg");
    if (bgRect) bgRect.setAttribute("fill", pal.surface);

    svg.querySelectorAll(".ow-region-box").forEach(el => {
      const isOx = el.classList.contains("ow-oxphos-box");
      el.setAttribute("fill", isOx ? pal.oxphos_face : pal.region_face);
      el.setAttribute("stroke", isOx ? pal.oxphos_edge : pal.region_edge);
      el.setAttribute("rx", isEcon ? "0.25" : (isOx ? "1.6" : "2.0"));
    });
    svg.querySelectorAll(".ow-region-rule").forEach(el => {
      el.style.display = isEcon ? "inline" : "none";
      el.setAttribute("stroke", pal.ink);
    });
    svg.querySelectorAll(".ow-region-title").forEach(g => {
      const dot = g.querySelector(".ow-region-dot");
      const bar = g.querySelector(".ow-region-accent");
      const txt = g.querySelector(".ow-region-text");
      const inMito = g.getAttribute("data-in-mito") === "1";
      const plateBg = inMito ? pal.mito_face : pal.region_face;
      if (dot) {
        dot.style.display = isEcon ? "none" : "inline";
        dot.setAttribute("stroke", pal.surface);
      }
      if (bar) {
        bar.style.display = isEcon ? "inline" : "none";
        bar.setAttribute("fill", isEcon ? pal.ink : (g.getAttribute("data-accent-metro") || pal.ink2));
      }
      if (txt) {
        txt.setAttribute("x", isEcon ? txt.getAttribute("data-x-econ") : txt.getAttribute("data-x-metro"));
        txt.setAttribute("fill", pal.ink);
        txt.setAttribute("stroke", plateBg);
      }
    });

    svg.querySelectorAll(".ow-comp").forEach(g => {
      const isP = g.getAttribute("data-comp") === "p";
      const box = g.querySelector(".ow-comp-box");
      const inn = g.querySelector(".ow-comp-inner");
      const lbl = g.querySelector(".ow-comp-lbl");
      const crx = isEcon ? 0.35 : (isP ? 2.4 : 2.6);
      if (box) {
        box.setAttribute("fill", isP ? pal.perox_face : pal.mito_face);
        box.setAttribute("stroke", isP ? pal.perox_edge : pal.mito_edge);
        box.setAttribute("rx", crx.toFixed(2));
      }
      if (inn) {
        inn.setAttribute("stroke", isP ? pal.perox_inner : pal.mito_inner);
        inn.setAttribute("rx", Math.max(0.2, crx - 0.35).toFixed(2));
      }
      if (lbl) lbl.setAttribute("fill", isP ? pal.perox_ink : pal.mito_ink);
    });

    svg.querySelectorAll(".ow-conduit-fill").forEach(el => el.setAttribute("stroke", pal.coupled_track));
    svg.querySelectorAll(".ow-conduit-edge").forEach(el => el.setAttribute("stroke", pal.coupled_edge));

    recolorEdges(pal);

    svg.querySelectorAll(".ow-node").forEach(g => {
      const valStr = g.getAttribute("data-val");
      const st = g.getAttribute("data-state");
      const ring = g.getAttribute("data-ring");
      const sig = ring === "selected" || ring === "changed";
      const halo = g.querySelector(".ow-node-halo");
      const dot = g.querySelector(".ow-node-dot");
      const oring = g.querySelector(".ow-node-ring");
      const lbl = g.querySelector(".ow-node-lbl");
      const measured = valStr !== "" && !isNaN(Number(valStr));
      if (halo) halo.setAttribute("fill", pal.surface);
      if (dot) {
        dot.setAttribute("fill", measured ? divergingColor(Number(valStr), pal) : pal.missing_fill);
        dot.setAttribute("stroke", sig ? pal.ink : (measured ? pal.ink2 : pal.omics_grey));
      }
      if (oring) oring.setAttribute("stroke", pal.ink);
      if (lbl) {
        lbl.setAttribute("fill", (sig || measured) ? pal.ink : pal.muted);
        lbl.setAttribute("stroke", pal.surface);
      }
    });
    svg.querySelectorAll(".ow-boundary .ow-pill-bg").forEach(el => {
      el.setAttribute("fill", pal.boundary_face);
      el.setAttribute("stroke", pal.boundary_edge);
      el.setAttribute("rx", isEcon ? "0.22" : "0.72");
    });
    svg.querySelectorAll(".ow-boundary .ow-pill-txt").forEach(el => el.setAttribute("fill", pal.boundary_ink));
    svg.querySelectorAll(".ow-store .ow-pill-bg").forEach(el => {
      el.setAttribute("fill", pal.store_face);
      el.setAttribute("stroke", pal.store_edge);
      el.setAttribute("rx", isEcon ? "0.20" : "0.72");
    });
    svg.querySelectorAll(".ow-store .ow-pill-txt").forEach(el => el.setAttribute("fill", pal.store_ink));

    svg.querySelectorAll(".ow-strip .ow-card-bg").forEach(el => {
      const st = el.closest(".ow-strip").getAttribute("data-tier-style");
      el.setAttribute("fill", cardSurface);
      el.setAttribute("stroke", st === "key" ? pal.ink : (st === "l1" ? pal.ink2 : pal.card_edge));
      el.setAttribute("rx", isEcon ? "0.12" : "0.38");
    });
    svg.querySelectorAll(".ow-agg .ow-card-bg").forEach(el => {
      el.setAttribute("fill", cardSurface);
      el.setAttribute("rx", isEcon ? "0.16" : "0.55");
    });
    svg.querySelectorAll(".ow-strip-frame").forEach(el => {
      el.setAttribute("fill", cardSurface);
      el.setAttribute("stroke", pal.card_edge);
    });
    svg.querySelectorAll(".ow-strip-hdr").forEach(el => el.setAttribute("fill", mixHex(cardSurface, pal.ink2, 0.09)));
    svg.querySelectorAll(".ow-prm-sq").forEach(el => {
      const v = el.getAttribute("data-val");
      el.setAttribute("fill", (v !== "" && !isNaN(Number(v))) ? divergingColor(Number(v), pal) : pal.missing_fill);
      el.setAttribute("stroke", pal.card_edge);
    });
    svg.querySelectorAll(".ow-rank-tab").forEach(el => {
      const top3 = el.getAttribute("data-top3") === "1";
      const bg = top3 ? (isEcon ? pal.red_tag : pal.ink) : pal.card_bg;
      el.setAttribute("fill", bg);
      el.setAttribute("stroke", top3 ? bg : pal.ink2);
    });
    svg.querySelectorAll(".ow-rank-txt").forEach(el => {
      const top3 = el.getAttribute("data-top3") === "1";
      el.setAttribute("fill", top3 ? (isEcon ? "#FFFFFF" : pal.surface) : pal.ink);
    });
    svg.querySelectorAll(".ow-strip-compact .ow-rxn-lbl").forEach(el => {
      el.setAttribute("stroke", pal.surface);
      el.setAttribute("fill", el.closest(".ow-strip").getAttribute("data-draw") === "not_assessed" ? pal.muted : pal.ink2);
    });
    svg.querySelectorAll(".ow-strip:not(.ow-strip-compact) .ow-rxn-lbl").forEach(el => {
      const st = el.closest(".ow-strip").getAttribute("data-tier-style");
      el.setAttribute("fill", (st === "key" || st === "l1") ? pal.ink : pal.ink2);
    });
    svg.querySelectorAll(".ow-card-title").forEach(el => el.setAttribute("fill", pal.ink));
    svg.querySelectorAll(".ow-tier-badge").forEach(g => {
      const t = pal.tier[g.getAttribute("data-tier")];
      const sc = parseFloat(g.getAttribute("data-scale")) || 1.0;
      if (!t) return;
      const bg = g.querySelector(".ow-tb-bg");
      if (bg) { bg.setAttribute("fill", t.bg); bg.setAttribute("stroke", t.ec); bg.setAttribute("rx", ((isEcon ? 0.12 : 0.38) * sc).toFixed(2)); }
      g.querySelectorAll(".ow-tb-dot-on").forEach(c => c.setAttribute("fill", t.dot_on));
      g.querySelectorAll(".ow-tb-dot-off").forEach(c => c.setAttribute("stroke", t.dot_off));
      const tx = g.querySelector(".ow-tb-txt");
      if (tx) tx.setAttribute("fill", t.fg);
    });
    svg.querySelectorAll(".ow-coupled-badge").forEach(g => {
      const sc = parseFloat(g.getAttribute("data-scale")) || 1.0;
      const bg = g.querySelector(".ow-cb-bg");
      if (bg) { bg.setAttribute("fill", pal.coupled_fill); bg.setAttribute("stroke", pal.coupled_edge); bg.setAttribute("rx", ((isEcon ? 0.12 : 0.36) * sc).toFixed(2)); }
      const tx = g.querySelector(".ow-cb-txt");
      if (tx) tx.setAttribute("fill", pal.coupled_ink);
    });
    svg.querySelectorAll(".ow-cell-sym").forEach(el => {
      const cr = el.getAttribute("data-crole") || "";
      const col = cr.includes("up") ? pal.up : pal.down;
      el.setAttribute("stroke", col);
      el.setAttribute("fill", cr.includes("solid") ? col : "none");
    });
    svg.querySelectorAll(".ow-cell-bg").forEach(el => {
      const cr = el.getAttribute("data-crole") || "";
      const col = cr.includes("up") ? pal.up : pal.down;
      el.setAttribute("fill", mixHex(cardSurface, col, cr.includes("solid") ? 0.14 : 0.06));
    });

    renderLegendFloat(pal);
    if (state.selectedKey) renderDrawer(state.selectedKey); else renderOverview();
  }

  function recolorEdges(pal) {
    svg.querySelectorAll(".ow-edge").forEach(g => {
      const role = g.getAttribute("data-erole");
      const valStr = g.getAttribute("data-val");
      const hollow = g.getAttribute("data-hollow") === "1";
      let col = pal.edge_grey;
      if (role === "mfa_dir" && valStr !== "") {
        const v = Number(valStr);
        const eff = Math.sign(v) * Math.max(Math.abs(v), 0.8);
        col = divergingColor(eff, pal);
      } else if (role === "enz_up" || role === "enz_down") {
        if (state.omicsColor === "tinted") {
          col = enzTintColor(role === "enz_up" ? 1 : -1, Math.abs(Number(valStr || 1)), pal);
        } else {
          col = pal.omics_grey;
        }
      } else if (role === "omics_grey") {
        col = pal.omics_grey;
      } else if (role === "not_assessed") {
        col = pal.faint;
      }
      g.querySelectorAll(".ow-edge-main").forEach(p => p.setAttribute("stroke", col));
      g.querySelectorAll(".ow-edge-halo").forEach(p => p.setAttribute("stroke", pal.surface));
      g.querySelectorAll(".ow-edge-core").forEach(p => p.setAttribute("stroke", pal.surface));
      g.querySelectorAll(".ow-edge-chev").forEach(p => p.setAttribute("stroke", hollow ? col : pal.surface));
      g.querySelectorAll(".ow-edge-hub").forEach(c => { c.setAttribute("fill", col); c.setAttribute("stroke", pal.surface); });
      g.querySelectorAll(".ow-edge-note").forEach(t => { t.setAttribute("fill", pal.muted); t.setAttribute("stroke", pal.surface); });
    });
  }

  /* ---------- Viewport Pan & Zoom ---------- */
  function setViewBox(v) {
    state.vb = v;
    svg.setAttribute("viewBox", v.map(n => n.toFixed(2)).join(" "));
  }
  function animateViewBox(target, ms) {
    const start = [...state.vb];
    const t0 = performance.now();
    const dur = ms || 280;
    function step(now) {
      const p = Math.min(1, (now - t0) / dur);
      const e = 1 - Math.pow(1 - p, 3);
      setViewBox(start.map((s, i) => s + (target[i] - s) * e));
      if (p < 1) requestAnimationFrame(step);
    }
    requestAnimationFrame(step);
  }
  function zoomAt(factor, clientX, clientY) {
    const rect = svg.getBoundingClientRect();
    const rx = (clientX - rect.left) / rect.width;
    const ry = (clientY - rect.top) / rect.height;
    const [x, y, w, h] = state.vb;
    const nw = Math.max(14, Math.min(state.baseVb[2] * 1.4, w * factor));
    const nh = nw * (h / w);
    setViewBox([x + (w - nw) * rx, y + (h - nh) * ry, nw, nh]);
  }
  function focusPoint(xy, span) {
    const s = span || 36;
    const ratio = state.baseVb[3] / state.baseVb[2];
    const w = s, h = s * ratio;
    animateViewBox([xy[0] - w / 2, xy[1] - h / 2, w, h], 300);
  }

  stage.addEventListener("wheel", e => {
    e.preventDefault();
    zoomAt(e.deltaY > 0 ? 1.15 : 0.87, e.clientX, e.clientY);
  }, { passive: false });
  stage.addEventListener("mousedown", e => {
    if (e.button !== 0) return;
    state.dragging = true;
    state.dragStart = [e.clientX, e.clientY];
    state.dragVb = [...state.vb];
    stage.classList.add("panning");
  });
  window.addEventListener("mousemove", e => {
    if (!state.dragging) return;
    const rect = svg.getBoundingClientRect();
    const dx = (e.clientX - state.dragStart[0]) / rect.width * state.dragVb[2];
    const dy = (e.clientY - state.dragStart[1]) / rect.height * state.dragVb[3];
    setViewBox([state.dragVb[0] - dx, state.dragVb[1] - dy, state.dragVb[2], state.dragVb[3]]);
  });
  window.addEventListener("mouseup", () => {
    state.dragging = false;
    stage.classList.remove("panning");
  });

  /* ---------- Filtering & Group Highlighting ---------- */
  function applyFilter(mode) {
    state.filter = mode;
    document.querySelectorAll("[data-filter]").forEach(b => b.classList.toggle("active", b.getAttribute("data-filter") === mode));
    svg.querySelectorAll(".ow-edge, .ow-strip").forEach(el => {
      const claim = el.getAttribute("data-claim") || "";
      const fg = el.getAttribute("data-flux-group") || "";
      const sc = el.getAttribute("data-sign-change") === "1";
      let keep = true;
      if (mode === "l2") keep = (claim === "L3" || claim === "L2");
      else if (mode === "l1") keep = (claim === "L1");
      else if (mode === "linked") keep = Boolean(fg);
      else if (mode === "reversal") keep = sc;
      el.classList.toggle("ow-dimmed", !keep);
    });
    svg.querySelectorAll(".ow-agg").forEach(el => {
      const claim = el.getAttribute("data-claim") || "";
      let keep = true;
      if (mode === "l2") keep = (claim === "L3" || claim === "L2");
      else if (mode === "l1") keep = (claim === "L1");
      else if (mode === "linked" || mode === "reversal") keep = false;
      el.classList.toggle("ow-dimmed", !keep);
    });
  }

  function highlightFluxGroup(fg, on) {
    svg.querySelectorAll(".ow-group-pulse").forEach(el => el.classList.remove("ow-group-pulse"));
    if (!fg || !on) return;
    svg.querySelectorAll(`[data-flux-group="${CSS.escape(fg)}"]`).forEach(el => el.classList.add("ow-group-pulse"));
  }

  function selectElement(elKey, shouldZoom) {
    svg.querySelectorAll(".ow-selected").forEach(e => e.classList.remove("ow-selected"));
    state.selectedKey = elKey;
    state.activeScatterLayer = "key";
    if (!elKey || !DATA.elements[elKey]) {
      renderOverview();
      return;
    }
    drawer.classList.remove("collapsed");
    const [kind, id] = elKey.split(":");
    const sel = kind === "rxn" ? `[data-rxn="${CSS.escape(id)}"]` : (kind === "agg" ? `[data-agg="${CSS.escape(id)}"]` : `[data-node="${CSS.escape(id)}"]`);
    svg.querySelectorAll(sel).forEach(e => e.classList.add("ow-selected"));
    const el = DATA.elements[elKey];
    if (shouldZoom && el.xy) focusPoint(el.xy, 38);
    if (el.flux_group) highlightFluxGroup(el.flux_group, true); else highlightFluxGroup("", false);
    renderDrawer(elKey);
  }

  /* ---------- Per-Animal Scatter Mini-Plot SVG Renderer (Economist Style!) ---------- */
  function seededJitter(seedStr) {
    let h = 2166136261;
    for (let i = 0; i < seedStr.length; i++) {
      h ^= seedStr.charCodeAt(i);
      h = Math.imul(h, 16777619);
    }
    return ((h >>> 0) % 1000) / 1000 - 0.5;
  }

  function renderScatterSvg(mem, series, pal) {
    const W = 192, H = 112;
    const ml = 8, mr = 34, mt = 8, mb = 22;
    const pw = W - ml - mr, ph = H - mt - mb;
    const groups = series.groups || [];
    const pts = series.pts || [];
    if (!pts.length) return `<div style="font-size:10.5px;color:var(--muted);padding:12px 2px;">No individual animal values</div>`;

    let vals = pts.map(p => p[4]).filter(v => v !== null && !isNaN(v));
    pts.forEach(p => {
      if (p[5] !== null && !isNaN(p[5])) vals.push(p[5]);
      if (p[6] !== null && !isNaN(p[6])) vals.push(p[6]);
    });
    let vmin = Math.min(...vals), vmax = Math.max(...vals);
    if (series.scale === "linear") {
      if (vmin > 0 && vmax > 0 && vmin / vmax < 0.25) vmin = 0;
      if (vmax < 0 && vmin < 0 && vmax / vmin < 0.25) vmax = 0;
    }
    let span = vmax - vmin;
    if (span < 1e-6) { const pad = Math.max(1.0, Math.abs(vmin) * 0.15); vmin -= pad; vmax += pad; span = vmax - vmin; }
    else { vmin -= span * 0.10; vmax += span * 0.12; span = vmax - vmin; }

    const yOf = v => mt + ph - ((v - vmin) / span) * ph;
    const xOf = gi => ml + (gi + 0.5) * (pw / Math.max(groups.length, 1));

    // Clean 4-tick horizontal Economist gridlines (right-aligned numeric ticks!)
    const ticks = [0, 1, 2, 3].map(i => vmin + (span * i) / 3);
    let out = `<svg class="sc-svg" viewBox="0 0 ${W} ${H}">`;
    ticks.forEach(tv => {
      const yy = yOf(tv);
      const tlbl = Math.abs(tv) >= 1000 ? tv.toExponential(1) : (Math.abs(tv) >= 10 ? tv.toFixed(0) : tv.toFixed(1));
      out += `<line x1="${ml}" y1="${yy.toFixed(1)}" x2="${W - mr + 2}" y2="${yy.toFixed(1)}" stroke="${pal.hairline}" stroke-width="0.8"/>`;
      out += `<text x="${W - mr + 5}" y="${yy.toFixed(1)}" font-size="8.5" fill="${pal.muted}" dominant-baseline="central">${tlbl}</text>`;
    });
    // Zero reference line for linear MFA net flux (Sarah Leo honesty guideline!)
    if (series.scale === "linear" && vmin <= 0 && vmax >= 0) {
      const y0 = yOf(0);
      out += `<line x1="${ml}" y1="${y0.toFixed(1)}" x2="${W - mr + 2}" y2="${y0.toFixed(1)}" stroke="${pal.ink}" stroke-width="1.0" stroke-dasharray="2.5 2"/>`;
    }
    // Economist charcoal bottom baseline
    out += `<line x1="${ml}" y1="${mt + ph}" x2="${W - mr + 2}" y2="${mt + ph}" stroke="${pal.ink}" stroke-width="1.15"/>`;

    // Group medians & animal points
    const byG = groups.map(() => []);
    pts.forEach(p => { if (byG[p[2]]) byG[p[2]].push(p); });
    const colW = pw / Math.max(groups.length, 1);

    byG.forEach((gpts, gi) => {
      const gx = xOf(gi);
      const gname = groups[gi] || "";
      const shortG = gname.length > 9 ? gname.slice(0, 8) + "…" : gname;
      out += `<line x1="${gx.toFixed(1)}" y1="${mt + ph}" x2="${gx.toFixed(1)}" y2="${mt + ph + 3}" stroke="${pal.ink}" stroke-width="0.9"/>`;
      out += `<text x="${gx.toFixed(1)}" y="${mt + ph + 12}" font-size="8.5" font-weight="600" fill="${pal.ink2}" text-anchor="middle">${esc(shortG)}</text>`;
      if (!gpts.length) return;
      const sortedV = gpts.map(p => p[4]).sort((a, b) => a - b);
      const med = sortedV[Math.floor(sortedV.length / 2)];
      const ym = yOf(med);
      const bw = Math.min(14, colW * 0.36);
      out += `<line x1="${(gx - bw).toFixed(1)}" y1="${ym.toFixed(1)}" x2="${(gx + bw).toFixed(1)}" y2="${ym.toFixed(1)}" stroke="${pal.ink}" stroke-width="1.6"/>`;

      gpts.forEach(p => {
        const [subj, sid, , sideCode, val, fl, fh, atB, qcf] = p;
        const jx = gx + seededJitter(sid + mem.name) * Math.min(14, colW * 0.42);
        const py = yOf(val);
        const col = sideCode === "n" ? pal.numerator : (sideCode === "d" ? pal.denominator : pal.neutral);
        const tip = `${subj} (${sid}) | ${gname} | ${fmtNum(val, 3)}${series.units ? " " + series.units : ""}` +
                    (fl !== null && fh !== null ? ` | fit IQR [${fmtNum(fl, 2)}, ${fmtNum(fh, 2)}]` : "") +
                    (atB ? " | at box bound" : "") + (qcf ? ` | QC: ${qcf}` : "");
        if (fl !== null && fh !== null) {
          out += `<line x1="${jx.toFixed(1)}" y1="${yOf(fl).toFixed(1)}" x2="${jx.toFixed(1)}" y2="${yOf(fh).toFixed(1)}" stroke="${col}" stroke-width="1.0" opacity="0.65"/>`;
        }
        if (atB) {
          // Hollow diamond for bound-pinned MFA animal
          const d = 3.6;
          out += `<polygon class="ow-pt" data-tip="${esc(tip)}" points="${jx.toFixed(1)},${(py - d).toFixed(1)} ${(jx + d).toFixed(1)},${py.toFixed(1)} ${jx.toFixed(1)},${(py + d).toFixed(1)} ${(jx - d).toFixed(1)},${py.toFixed(1)}" fill="${pal.surface}" stroke="${col}" stroke-width="1.4"/>`;
        } else {
          out += `<circle class="ow-pt" data-tip="${esc(tip)}" cx="${jx.toFixed(1)}" cy="${py.toFixed(1)}" r="3.2" fill="${col}" stroke="${qcf ? pal.ink : pal.surface}" stroke-width="${qcf ? 1.1 : 0.7}"${qcf ? ' stroke-dasharray="1.5 1"' : ''}/>`;
        }
      });
    });
    out += "</svg>";
    return out;
  }

  /* ---------- Right Drawer: Overview & Reaction Dossier ---------- */
  function renderOverview() {
    const pal = DATA.palettes[state.theme];
    const c = DATA.meta.counts;
    drawerHead.innerHTML = `
      <div class="dh-top">
        <div>
          <h2 class="dh-title">Pathway Atlas Dashboard</h2>
          <div class="dh-sub">Click any reaction track, metabolite station, or module card to inspect per-animal scatter plots &amp; member evidence.</div>
        </div>
        <button class="btn" id="btn-collapse-drawer" title="Collapse panel">Hide</button>
      </div>
      <div class="dh-tags">
        <span class="tag">${c.L3} Confirmed</span>
        <span class="tag">${c.L2} Flux${(c.L2 + c.L3) > c.indep_flux_chains ? ` (${c.indep_flux_chains} indep.)` : ""}</span>
        <span class="tag">${c.L1} Enz-pred</span>
        <span class="tag linked">${c.linked_groups} Linked chain${c.linked_groups === 1 ? "" : "s"} (shared 1 DOF)</span>
      </div>`;
    document.getElementById("btn-collapse-drawer").onclick = () => drawer.classList.add("collapsed");

    let html = `<div class="sec-title"><span>Top Ranked Reactions</span><span>Click to zoom</span></div>`;
    if (DATA.top_reactions.length) {
      DATA.top_reactions.forEach(tr => {
        const vcol = tr.value > 0 ? "var(--up)" : (tr.value < 0 ? "var(--down)" : "var(--muted)");
        html += `
          <div class="rank-item" data-jump="${esc(tr.key)}">
            <div><span class="rk-badge ${tr.rank <= 3 ? "top3" : ""}">#${tr.rank}</span> <strong>${esc(tr.short)}</strong> <span style="color:var(--muted);font-size:10.5px;">${esc(tr.region)}</span></div>
            <div style="display:flex;gap:6px;align-items:center;">
              <span style="font-weight:700;color:${vcol};">${tr.value !== null ? fmtSigned(tr.value, 2) : "—"}</span>
              <span class="tag">${esc(tr.claim || "—")}</span>
            </div>
          </div>`;
      });
    } else {
      html += `<div style="color:var(--muted);font-size:11.5px;">No ranked reactions in this contrast.</div>`;
    }

    html += `<div class="sec-title"><span>Pathway Coverage &amp; Fusion</span><span>Click region</span></div>`;
    DATA.regions.forEach(rg => {
      const tot = Math.max(rg.n_drawn, 1);
      const wId = (rg.n_identifiable / tot) * 100;
      const wUn = (Math.max(0, rg.n_flux - rg.n_identifiable) / tot) * 100;
      const wOm = (rg.n_omics_only / tot) * 100;
      html += `
        <div class="reg-row" data-zoom-region="${esc(rg.id)}">
          <div class="reg-top"><span>${esc(rg.label)}</span><span style="font-weight:600;color:var(--ink2);">L2+ ${rg.n_L2} · L1 ${rg.n_L1}</span></div>
          <div class="reg-bar">
            <div style="width:${wId.toFixed(1)}%;background:${pal.ink2};"></div>
            <div style="width:${wUn.toFixed(1)}%;background:${pal.edge_grey};"></div>
            <div style="width:${wOm.toFixed(1)}%;background:${pal.omics_grey};"></div>
          </div>
          <div class="reg-stat"><span>Flux ${rg.n_flux}/${rg.n_drawn} (${rg.n_identifiable} ident.) · Omics ${rg.n_omics_only}</span><span>${esc(rg.kegg.split(",")[0] || "")}</span></div>
          ${rg.fusion_sig.map(f => `<div style="font-size:10.5px;margin-top:2px;color:${f.strong ? "var(--up)" : "var(--ink2)"};">${f.strong ? "●" : "○"} ${esc(f.text)}</div>`).join("")}
        </div>`;
    });

    if (DATA.extra_fusion.length) {
      html += `<div class="sec-title"><span>Extra-Network Fusion (FDR &lt; 0.25)</span><span>Exploratory</span></div>`;
      DATA.extra_fusion.forEach(ef => {
        html += `<div style="padding:4px 0;border-bottom:1px solid var(--hairline);font-size:11px;display:flex;justify-content:space-between;gap:8px;">
          <span>${esc(ef.term)}</span><span style="color:${ef.fdr < 0.05 ? "var(--up)" : "var(--muted)"};font-weight:700;white-space:nowrap;">FDR ${fmtP(ef.fdr)}</span>
        </div>`;
      });
    }

    drawerBody.innerHTML = html;
    drawerBody.querySelectorAll("[data-jump]").forEach(el => {
      el.onclick = () => selectElement(el.getAttribute("data-jump"), true);
    });
    drawerBody.querySelectorAll("[data-zoom-region]").forEach(el => {
      el.onclick = () => zoomToRegion(el.getAttribute("data-zoom-region"));
    });
  }

  function renderFluxGateChecklist(el) {
    if (!el || el.kind !== "reaction" || el.draw_style !== "flux" || !el.axes || !el.axes.F) return "";
    const fState = String(el.axes.F.state || "");
    const fFlags = Array.isArray(el.f_flags) ? el.f_flags : [];
    const fRule = String(el.f_rule || "");
    const mfaMems = (el.members || []).filter(m => m.layer === "mfa");
    const pooled = mfaMems.find(m => m.tracer === "pooled") || null;
    const unident = new Set(["not_infused", "not_measured", "not_identifiable", "insufficient_n"]);
    const ident = mfaMems.filter(m => m.tracer !== "pooled" && !unident.has(m.state));
    const nIdent = ident.length || (el.flux_tracers ? el.flux_tracers.length : 0);
    const head = pooled || ident[0] || null;
    const pill = (cls, txt) => `<span class="gate-pill gate-${cls}">${esc(txt)}</span>`;

    const rows = [];
    // Gate 1: Identifiable tracers
    if (nIdent >= 2) {
      rows.push(["Identifiable tracers", `${nIdent} tracers (Sig. eligible)`, "pass", "PASS"]);
    } else if (nIdent === 1) {
      rows.push(["Identifiable tracers", "1 tracer (Suggestive max)", "warn", "1 ONLY"]);
    } else {
      rows.push(["Identifiable tracers", "0 identifiable tracers", "fail", "BLOCKED"]);
    }

    // Gate 2: Direction concordance
    if (fState === "mixed") {
      rows.push(["Tracer concordance", "Opposing clear effects (|g| >= 1)", "fail", "CONFLICT"]);
    } else if (nIdent >= 2 && ident.length >= 2) {
      const signs = new Set(ident.map(m => (m.log2fc > 0 ? 1 : (m.log2fc < 0 ? -1 : 0))).filter(s => s !== 0));
      if (signs.size <= 1 && signs.size > 0) {
        rows.push(["Tracer concordance", "All identifiable tracers agree", "pass", "PASS"]);
      } else {
        rows.push(["Tracer concordance", "Minor tracer differs in sign", "warn", "PARTIAL"]);
      }
    } else if (nIdent === 1) {
      rows.push(["Tracer concordance", "Single tracer (no cross-check)", "na", "N/A"]);
    } else {
      rows.push(["Tracer concordance", "No identifiable tracer", "na", "N/A"]);
    }

    // Gate 3: Standardized effect size |Hedges' g|
    const gVal = (el.f_g !== null && el.f_g !== undefined) ? el.f_g : (head ? head.effect : null);
    const gVar = pooled && pooled.g_total_var !== null && pooled.g_total_var !== undefined ? pooled.g_total_var : null;
    const varTxt = gVar !== null ? ` · Var=${fmtNum(gVar, 2)}` : "";
    const isPermSugg = fFlags.includes("multi_tracer_perm_suggestive");
    if (gVal !== null && gVal !== undefined) {
      if (fState === "unchanged") {
        rows.push(["Effect size |Hedges' g|", `g = ${fmtSigned(gVal, 2)} (TOST 90% CI < 1.0)${varTxt}`, "pass", "EQUIV"]);
      } else if (Math.abs(gVal) >= 1.0) {
        rows.push(["Effect size |Hedges' g|", `g = ${fmtSigned(gVal, 2)} (|g| >= 1.0)${varTxt}`, "pass", "PASS"]);
      } else if (isPermSugg) {
        rows.push(["Effect size |Hedges' g|", `g = ${fmtSigned(gVal, 2)} (k>=2 & p<0.05 & |g|>=0.8)${varTxt}`, "pass", "PASS (P+G)"]);
      } else {
        rows.push(["Effect size |Hedges' g|", `g = ${fmtSigned(gVal, 2)} (|g| < 1.0)${varTxt}`, "fail", "|g| < 1"]);
      }
    } else {
      rows.push(["Effect size |Hedges' g|", "Not available", "na", "N/A"]);
    }

    // Gate 4: Pooled FDR / p-value & heterogeneity
    const padj = head ? head.padj : null;
    const padjGlobal = head ? head.padj_global : null;
    const fdrScope = head && head.fdr_scope ? String(head.fdr_scope) : "";
    const pval = head ? head.pvalue : null;
    const hetP = head ? head.het_p : null;
    const nMin = head && head.n_num !== null && head.n_den !== null ? Math.min(head.n_num, head.n_den) : null;
    const hetBad = fFlags.includes("tracer_heterogeneity") || (hetP !== null && hetP < 0.10);
    const statParts = [];
    if (padj !== null) {
      const scopeTag = fdrScope === "per_tissue" ? "(tissue)" : (fdrScope === "per_contrast" ? "(global)" : "");
      statParts.push(`FDR${scopeTag}=${fmtP(padj)}`);
    }
    if (padjGlobal !== null && fdrScope === "per_tissue" && padjGlobal !== padj) {
      statParts.push(`FDR(global)=${fmtP(padjGlobal)}`);
    }
    if (pval !== null) statParts.push(`p=${fmtP(pval)}`);
    if (hetP !== null) statParts.push(`Q p=${fmtP(hetP)}`);
    if (nMin !== null) statParts.push(`n_min=${nMin}`);
    if (padj !== null && padj < 0.05 && !hetBad && (nMin === null || nMin >= 3)) {
      rows.push(["Significance & Q test", statParts.join(" · "), "pass", "SIG"]);
    } else if (hetBad) {
      rows.push(["Significance & Q test", statParts.join(" · "), "warn", "HETERO"]);
    } else if (pval !== null && pval < 0.05) {
      rows.push(["Significance & Q test", statParts.join(" · "), "warn", "NOMINAL"]);
    } else if (statParts.length) {
      rows.push(["Significance & Q test", statParts.join(" · "), "fail", "NS"]);
    } else {
      rows.push(["Significance & Q test", "No pooled test statistics", "na", "N/A"]);
    }

    // Gate 5: Directional stability & robustness
    const memFlags = head && head.flags ? String(head.flags) : "";
    const singleSubj = fFlags.includes("single_subject_driven") || memFlags.includes("single_subject_driven");
    const mismatch = fFlags.includes("tracer_expectation_mismatch");
    const stabUnk = fFlags.includes("stability_unknown");
    const dirStab = pooled ? pooled.dir_stab : (head ? head.dir_stab : null);
    const looStab = pooled ? pooled.loo_stab : (head ? head.loo_stab : null);
    const stabInfo = [];
    if (dirStab !== null && dirStab !== undefined) stabInfo.push(`boot=${fmtNum(dirStab, 2)}`);
    if (looStab !== null && looStab !== undefined) stabInfo.push(`LOO=${fmtNum(looStab, 2)}`);
    if (mismatch) {
      rows.push(["Stability & expectation", "Identifiability mismatches panel expectation", "fail", "MISMATCH"]);
    } else if (singleSubj) {
      rows.push(["Stability & expectation", `Single-subject driven${stabInfo.length ? " (" + stabInfo.join(", ") + ")" : ""}`, "fail", "1-SUBJ"]);
    } else if (stabUnk) {
      rows.push(["Stability & expectation", "Per-animal stability unavailable", "warn", "UNKNOWN"]);
    } else if ((dirStab !== null && dirStab < 0.8) || (looStab !== null && looStab < 1.0) || fRule.includes("LOO") || fRule.includes("不稳定")) {
      rows.push(["Stability & expectation", stabInfo.length ? stabInfo.join(" · ") : "Direction flips on LOO", "fail", "UNSTABLE"]);
    } else if (fState === "changed" || fState === "suggestive" || (dirStab !== null && dirStab >= 0.8 && (looStab === null || looStab >= 1.0))) {
      rows.push(["Stability & expectation", stabInfo.length ? stabInfo.join(" · ") : "Passed stratified LOO", "pass", "STABLE"]);
    } else {
      rows.push(["Stability & expectation", stabInfo.length ? stabInfo.join(" · ") : "Not blocking", "na", "N/A"]);
    }

    const stateTag = fState === "changed" ? pill("pass", "Sig. (Established)")
      : (fState === "suggestive" ? pill("warn", "Suggestive")
      : (fState === "unchanged" ? pill("pass", "Unchanged (Equiv.)")
      : (fState === "mixed" ? pill("fail", "Conflicting")
      : pill("na", fState.replace(/_/g, " ") || "undetermined"))));
    let out = `<div class="gate-box"><div class="gate-hdr"><span>Flux Assessment Gate Checklist</span>${stateTag}</div>`;
    rows.forEach(([lbl, detail, cls, badge]) => {
      out += `<div class="gate-row"><span><strong>${esc(lbl)}:</strong> <span style="color:var(--ink2);">${esc(detail)}</span></span>${pill(cls, badge)}</div>`;
    });
    const meta = [];
    if (pooled && pooled.smd_weights) {
      const modeLbl = pooled.smd_weighting ? ` (${pooled.smd_weighting})` : "";
      meta.push(`<strong>SMD weights${esc(modeLbl)}:</strong> ${esc(pooled.smd_weights)}`);
    }
    if (pooled && pooled.tracer_weights) {
      meta.push(`<strong>Abs IVW weights:</strong> ${esc(pooled.tracer_weights)}`);
    }
    if (fFlags.length) {
      meta.push(`<strong>Flags:</strong> ${esc(fFlags.join(", "))}`);
    }
    if (fRule) {
      meta.push(`<strong>Rule:</strong> ${esc(fRule)}`);
    }
    if (meta.length) {
      out += `<div class="gate-meta">${meta.join("<br>")}</div>`;
    }
    out += `</div>`;
    return out;
  }

  function renderDrawer(elKey) {
    const el = DATA.elements[elKey];
    if (!el) return renderOverview();
    const pal = DATA.palettes[state.theme];
    const dirClass = el.value > 0 || el.enz_dir > 0 ? "up" : ((el.value < 0 || el.enz_dir < 0) ? "down" : "");
    const claimDots = { L3: "●●●", L2: "●●○", L1: "●○○" };
    const claimFallback = { L3: "Confirmed", L2: "Flux", L1: "Enz-pred" };
    const claimBadgeTxt = el.flux_claim && claimDots[el.flux_claim]
      ? `${claimDots[el.flux_claim]} ${el.tier_label || claimFallback[el.flux_claim]}`
      : "";

    let tags = "";
    if (claimBadgeTxt) tags += `<span class="tag">${esc(claimBadgeTxt)}</span>`;
    if (el.value !== null && el.value !== undefined) tags += `<span class="tag ${dirClass}">log2FC ${fmtSigned(el.value, 2)}</span>`;
    if (el.rank) tags += `<span class="tag">Rank #${el.rank}</span>`;
    if (el.flux_group) tags += `<span class="tag linked" id="btn-hl-linked" title="Stoichiometrically locked chain (shared 1 DOF, NOT extra confidence) — click to highlight">&#128279; Linked flux: ${esc(el.flux_group)} (1 DOF)</span>`;
    if (el.sign_change) tags += `<span class="tag" style="background:${pal.sign_change_bg};color:${pal.sign_change_ink};border-color:${pal.sign_change_edge};">&#8644; Net-flux sign reversal</span>`;
    if (el.ec_number) tags += `<span class="tag">EC ${esc(el.ec_number)}</span>`;
    if (el.kegg_reaction) tags += `<span class="tag">KEGG ${esc(el.kegg_reaction)}</span>`;

    drawerHead.innerHTML = `
      <div class="dh-top">
        <div>
          <h2 class="dh-title">${esc(el.short)} <span style="font-weight:500;color:var(--muted);font-size:12px;">${el.id !== el.short ? esc(el.id) : ""}</span></h2>
          <div class="dh-sub">${esc(el.region_label)}${el.equation ? " · " + esc(el.equation) : ""}</div>
        </div>
        <div style="display:flex;gap:4px;">
          <button class="btn" id="btn-back-overview">Overview</button>
          <button class="btn" id="btn-collapse-drawer">Hide</button>
        </div>
      </div>
      <div class="dh-tags">${tags}</div>`;

    document.getElementById("btn-back-overview").onclick = () => selectElement(null, false);
    document.getElementById("btn-collapse-drawer").onclick = () => drawer.classList.add("collapsed");
    const hlBtn = document.getElementById("btn-hl-linked");
    if (hlBtn) hlBtn.onclick = () => highlightFluxGroup(el.flux_group, true);

    let body = "";
    // 1. Four-omics / four-axis summary strip
    if (el.kind !== "metabolite" && el.axes) {
      const cells = [
        { lbl: "Flux (MFA)", val: el.value, st: el.axes.F.state, dr: el.axes.F.dir },
        { lbl: "Protein (P)", val: el.layers.P.val, st: el.layers.P.state || el.axes.C.state, dr: el.axes.C.dir },
        { lbl: "mRNA (R)", val: el.layers.R.val, st: el.layers.R.state, dr: el.layers.R.val > 0 ? 1 : (el.layers.R.val < 0 ? -1 : 0) },
        { lbl: "Metab (M)", val: el.layers.M.val, st: el.layers.M.state || el.axes.M.state, dr: el.axes.M.dir },
      ];
      body += `<div class="axis-grid">`;
      cells.forEach(c => {
        const col = (c.st === "changed" || c.st === "suggestive") ? (c.dr > 0 ? "var(--up)" : (c.dr < 0 ? "var(--down)" : "var(--ink)")) : "var(--ink2)";
        body += `<div class="axis-card">
          <div class="ac-lbl">${c.lbl}</div>
          <div class="ac-val" style="color:${col};">${c.val !== null && c.val !== undefined ? fmtSigned(c.val, 2) : "—"}</div>
          <div class="ac-st">${esc(String(c.st || "—").replace(/_/g, " "))}</div>
        </div>`;
      });
      body += `</div>`;
      body += renderFluxGateChecklist(el);
      if (el.flux_tracers && el.flux_tracers.length) {
        const tList = el.flux_tracers.join(", ");
        let modeDesc = "";
        if (el.flux_pooling_mode === "multi_tracer_pooled") {
          modeDesc = `<strong>Pooled IVW (${el.flux_tracers.length} identifiable tracers):</strong> ${esc(tList)} — combined via tracer-stratified inverse-variance fixed-effects pooling (unconstrained / low-n tracers gated out).`;
        } else if (el.flux_pooling_mode === "single_tracer") {
          modeDesc = `<strong>Active MFA Tracer:</strong> ${esc(tList)} (one identifiable tracer with sufficient animals).`;
        } else if (el.flux_pooling_mode === "multi_tracer_conflict") {
          modeDesc = `<strong>Multi-Tracer Conflict:</strong> ${esc(tList)} have opposing clear effects (|g| ≥ 1); no pooled flux direction assigned.`;
        }
        if (modeDesc) {
          body += `<div style="font-size:10.5px;color:var(--ink2);background:var(--surface);border-left:3px solid var(--ink2);padding:5px 8px;margin-bottom:8px;">${modeDesc}</div>`;
        }
      }
      if (el.flux_claim === "L1") {
        body += `<div style="font-size:10.5px;color:var(--ink2);background:var(--surface);border-left:3px solid var(--up);padding:5px 8px;margin-bottom:8px;">
          <strong>Enz-pred (L1):</strong> Direction inferred from enzyme capacity (Protein/mRNA), not measured MFA flux. Verify against decoy baseline before citing.
        </div>`;
      }
      if (el.flux_group) {
        body += `<div style="font-size:10.5px;color:var(--coupled-ink);background:var(--coupled-bg);border-left:3px solid var(--coupled-edge);padding:5px 8px;margin-bottom:8px;">
          <strong>Linked flux (${esc(el.flux_group)} — Shared 1 DOF, NOT extra confidence):</strong> Per-fit net fluxes maintain |r| ≥ 1−1e-6 across all independent fits (steady-state mass balance on an unbranched segment). Counted as 1 single hypothesis in grouped BH-FDR to prevent multiple-testing inflation and duplicate discovery counting.
        </div>`;
      }
    } else if (el.kind === "metabolite") {
      body += `<div style="font-size:11px;color:var(--ink2);background:var(--surface);padding:6px 8px;margin-bottom:8px;border:1px solid var(--card-edge);">
        Whole-tissue metabolomics pool: <strong>${el.value !== null ? "log2FC " + fmtSigned(el.value, 2) : esc(el.display_state)}</strong>
        ${el.n_quantified ? ` (${el.n_changed}/${el.n_quantified} changed)` : ""}. Compartment nodes share one tissue pool value.
      </div>`;
    }

    // 2. Per-Animal / Per-Sample Scatter Mini-Plots (Economist Style!)
    const mems = el.members || [];
    const memsWithSamples = mems.filter(m => m.skey && DATA.samples[m.skey]);
    const layerLabels = { key: "Key Multi-Omics", mfa: "MFA Flux", proteomics: "Protein", transcriptomics: "mRNA", metabolomics: "Metabolite", derived: "Ratio" };
    const availLayers = ["key"];
    ["mfa", "proteomics", "transcriptomics", "metabolomics", "derived"].forEach(l => {
      if (memsWithSamples.some(m => m.layer === l)) availLayers.push(l);
    });
    if (!availLayers.includes(state.activeScatterLayer)) state.activeScatterLayer = "key";

    body += `<div class="sec-title"><span>Per-Animal / Sample Distribution</span><span>bar = median</span></div>`;
    body += `<div class="sc-legend">
      <span><i class="sc-dot" style="background:${pal.numerator};"></i>${esc(DATA.meta.num_label)} (Exp)</span>
      <span><i class="sc-dot" style="background:${pal.denominator};"></i>${esc(DATA.meta.den_label)} (Ctrl)</span>
      <span><i class="sc-dot" style="background:${pal.neutral};"></i>Other groups</span>
      <span style="color:var(--muted);">◇ bound-pinned</span>
    </div>`;

    if (memsWithSamples.length) {
      body += `<div class="layer-tabs">`;
      availLayers.forEach(l => {
        const cnt = l === "key" ? Math.min(6, memsWithSamples.length) : memsWithSamples.filter(m => m.layer === l).length;
        body += `<button class="btn ${state.activeScatterLayer === l ? "active" : ""}" data-sc-layer="${l}">${layerLabels[l]} (${cnt})</button>`;
      });
      body += `</div>`;

      let chosen = [];
      if (state.activeScatterLayer === "key") {
        // Pick up to 2 per layer so user sees Flux + Protein + mRNA + Metabolite together!
        ["mfa", "proteomics", "transcriptomics", "metabolomics", "derived"].forEach(l => {
          const sub = memsWithSamples.filter(m => m.layer === l && (l !== "mfa" || m.tracer !== "pooled"));
          const fallback = sub.length ? sub : memsWithSamples.filter(m => m.layer === l);
          chosen.push(...fallback.slice(0, l === "mfa" ? 2 : 2));
        });
        chosen = chosen.slice(0, 6);
      } else {
        chosen = memsWithSamples.filter(m => m.layer === state.activeScatterLayer).slice(0, 12);
      }

      const gridCls = chosen.length === 1 ? "scatter-grid single-col" : "scatter-grid";
      body += `<div class="${gridCls}">`;
      chosen.forEach(m => {
        const ser = DATA.samples[m.skey];
        const tag = { mfa: "MFA", proteomics: "Prot", transcriptomics: "mRNA", metabolomics: "Metab", derived: "Ratio" }[m.layer] || m.layer;
        const statParts = [];
        if (m.log2fc !== null) statParts.push(`FC ${fmtSigned(m.log2fc, 2)}`);
        if (m.padj !== null) statParts.push(`q=${fmtP(m.padj)}`);
        else if (m.pvalue !== null) statParts.push(`p=${fmtP(m.pvalue)}`);
        if (m.share !== null) statParts.push(`${Math.round(m.share * 100)}%`);
        body += `<div class="sc-card">
          <div class="sc-head" title="${esc(m.name)} (${tag})">${esc(m.name)} <span style="font-weight:500;color:var(--muted);">· ${tag}</span></div>
          <div class="sc-sub">${esc(statParts.join(" · ") || m.state)}${ser.units ? " · " + esc(ser.units) : ""}</div>
          ${renderScatterSvg(m, ser, pal)}
        </div>`;
      });
      body += `</div>`;
    } else {
      body += `<div style="font-size:11px;color:var(--muted);padding:8px;background:var(--surface);border:1px solid var(--card-edge);">
        ${DATA.meta.has_samples ? "No quantified individual animal series linked to this element." : "Per-animal sample table (dossier_samples.parquet) was not bundled with this export."}
      </div>`;
    }

    // 3. Full Member Evidence Table
    body += `<div class="sec-title"><span>Mapped Multi-Omics Members (${mems.length})</span><span>Click header to sort</span></div>`;
    if (mems.length) {
      const sorted = [...mems].sort((a, b) => {
        const sc = state.memSort.col;
        const av = sc === "abs_fc" ? Math.abs(a.log2fc ?? -1) : (sc === "padj" ? (a.padj ?? a.pvalue ?? 99) : String(a[sc] || ""));
        const bv = sc === "abs_fc" ? Math.abs(b.log2fc ?? -1) : (sc === "padj" ? (b.padj ?? b.pvalue ?? 99) : String(b[sc] || ""));
        if (typeof av === "string") return state.memSort.desc ? bv.localeCompare(av) : av.localeCompare(bv);
        return state.memSort.desc ? bv - av : av - bv;
      });
      body += `<table class="mem-tbl"><thead><tr>
        <th data-sort="layer">Layer</th><th data-sort="name">Member</th>
        <th data-sort="abs_fc">log2FC</th><th data-sort="padj">padj/p</th><th>State</th>
      </tr></thead><tbody>`;
      const lShort = { mfa: "Flux", proteomics: "Prot", transcriptomics: "mRNA", metabolomics: "Metab", derived: "Ratio" };
      sorted.forEach(m => {
        const fc = m.log2fc;
        const barW = fc !== null ? Math.min(34, Math.abs(fc) * 14) : 0;
        const barCol = fc > 0 ? pal.up : pal.down;
        const pv = m.padj !== null ? fmtP(m.padj) : (m.pvalue !== null ? fmtP(m.pvalue) + "*" : "—");
        const gSub = (m.layer === "mfa" && m.effect !== null && m.effect !== undefined)
          ? `<div style="font-size:9.5px;color:var(--muted);font-family:var(--font-mono);">g=${fmtSigned(m.effect, 2)}</div>`
          : "";
        body += `<tr>
          <td><span class="tag">${lShort[m.layer] || m.layer}</span></td>
          <td title="${esc(m.role)} ${esc(m.relation)} ${esc(m.reason)}"><strong>${esc(m.name)}</strong>${m.share !== null ? ` <span style="color:var(--muted);font-size:10px;">(${Math.round(m.share * 100)}%)</span>` : ""}</td>
          <td><span class="fc-bar-wrap">${fc !== null ? `<span class="fc-bar" style="width:${barW.toFixed(0)}px;background:${barCol};"></span>` : ""}<span style="color:${fc !== null ? barCol : "var(--muted)"};font-weight:700;">${fmtSigned(fc, 2)}</span></span>${gSub}</td>
          <td style="font-family:var(--font-mono);font-size:10px;">${pv}</td>
          <td style="color:var(--muted);font-size:10px;">${esc(m.state.replace(/_/g, " "))}</td>
        </tr>`;
      });
      body += `</tbody></table>`;
    } else {
      body += `<div style="font-size:11px;color:var(--muted);">No mapped omics members in knowledge base.</div>`;
    }

    drawerBody.innerHTML = body;
    drawerBody.querySelectorAll("[data-sc-layer]").forEach(btn => {
      btn.onclick = () => { state.activeScatterLayer = btn.getAttribute("data-sc-layer"); renderDrawer(elKey); };
    });
    drawerBody.querySelectorAll("th[data-sort]").forEach(th => {
      th.onclick = () => {
        const col = th.getAttribute("data-sort");
        if (state.memSort.col === col) state.memSort.desc = !state.memSort.desc;
        else state.memSort = { col, desc: col !== "padj" && col !== "name" };
        renderDrawer(elKey);
      };
    });
    drawerBody.querySelectorAll(".ow-pt").forEach(pt => {
      pt.addEventListener("mouseenter", e => showRawTooltip(e, pt.getAttribute("data-tip")));
      pt.addEventListener("mousemove", moveTooltip);
      pt.addEventListener("mouseleave", hideTooltip);
    });
  }

  /* ---------- Tooltip Helpers ---------- */
  function showElementTooltip(e, elKey) {
    const el = DATA.elements[elKey];
    if (!el) return;
    const claimDots = { L3: "●●●", L2: "●●○", L1: "●○○" };
    const claimFallback = { L3: "Confirmed", L2: "Flux", L1: "Enz-pred" };
    const claimText = el.flux_claim && claimDots[el.flux_claim]
      ? `${claimDots[el.flux_claim]} ${el.tier_label || claimFallback[el.flux_claim]}`
      : (el.flux_claim || "");
    const fluxSrc = (el.flux_tracers && el.flux_tracers.length)
      ? (el.flux_pooling_mode === "multi_tracer_pooled"
          ? `Pooled IVW: ${el.flux_tracers.join(" + ")}`
          : (el.flux_pooling_mode === "single_tracer"
              ? `Tracer: ${el.flux_tracers[0]} (single-tracer)`
              : `Conflict: ${el.flux_tracers.join(" vs ")}`))
      : "";
    let valHtml = "";
    if (el.value !== null && el.value !== undefined) {
      const vc = el.value > 0 ? "var(--up)" : (el.value < 0 ? "var(--down)" : "var(--ink2)");
      valHtml = `<span style="font-weight:800;color:${vc};">log2FC ${fmtSigned(el.value, 2)}</span>`;
    } else if (el.enz_dir) {
      const vc = el.enz_dir > 0 ? "var(--up)" : "var(--down)";
      valHtml = `<span style="font-weight:700;color:${vc};">Enz-pred ${el.enz_dir > 0 ? "▲ Up" : "▼ Down"}</span>`;
    }
    const topMem = (el.members || []).filter(m => m.log2fc !== null).slice(0, 3).map(m => `${m.name} (${fmtSigned(m.log2fc, 2)})`).join(", ");
    tooltip.innerHTML = `
      <div class="tt-head"><span class="tt-title">${esc(el.short)}</span>${valHtml}</div>
      ${el.equation ? `<div class="tt-eq">${esc(el.equation)}</div>` : ""}
      <div class="tt-meta">${esc(el.region_label)}${claimText ? " · " + esc(claimText) : ""}${el.flux_group ? " · Linked (1 DOF): " + esc(el.flux_group) : ""}</div>
      ${fluxSrc ? `<div class="tt-meta" style="margin-top:2px;color:var(--ink);font-weight:600;">${esc(fluxSrc)}</div>` : ""}
      ${topMem ? `<div class="tt-meta" style="margin-top:3px;color:var(--ink2);">Top: ${esc(topMem)}</div>` : ""}
      <div class="tt-meta" style="margin-top:3px;">Click for per-animal plots &amp; full dossier</div>`;
    tooltip.style.display = "block";
    moveTooltip(e);
  }
  function showRawTooltip(e, text) {
    tooltip.innerHTML = `<div style="font-size:11px;font-weight:600;">${esc(text)}</div>`;
    tooltip.style.display = "block";
    moveTooltip(e);
  }
  function moveTooltip(e) {
    const pad = 14;
    let lx = e.clientX + pad, ly = e.clientY + pad;
    if (lx + 320 > window.innerWidth) lx = e.clientX - 320;
    if (ly + 140 > window.innerHeight) ly = e.clientY - 120;
    tooltip.style.left = Math.max(6, lx) + "px";
    tooltip.style.top = Math.max(6, ly) + "px";
  }
  function hideTooltip() { tooltip.style.display = "none"; }

  /* ---------- Interactive Legend Overlay ---------- */
  function renderLegendFloat(pal) {
    const softUp = enzTintColor(1, 1.2, pal);
    const softDn = enzTintColor(-1, 1.2, pal);
    const el = document.getElementById("ow-legend-float");
    el.innerHTML = `
      <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:6px;">
        <strong>Visual Encoding &amp; Scientific Rules</strong>
        <button class="btn" id="btn-close-lg" style="padding:1px 6px;">×</button>
      </div>
      <div class="legend-grid">
        <div class="lg-row"><span class="lg-swatch"><svg width="32" height="10"><line x1="1" y1="5" x2="31" y2="5" stroke="${pal.up}" stroke-width="4"/></svg></span><span><strong>Solid Red / Blue rail:</strong> ${DATA.meta.tracer ? "Net flux assessed within the displayed tracer cohort." : "MFA net flux assessed across informative tracers (k=1: single-tracer evidence)."}</span></div>
        <div class="lg-row"><span class="lg-swatch"><svg width="32" height="10"><line x1="1" y1="5" x2="15" y2="5" stroke="${softUp}" stroke-width="4"/><line x1="2" y1="5" x2="14" y2="5" stroke="${pal.card_bg}" stroke-width="1.8"/><line x1="17" y1="5" x2="31" y2="5" stroke="${softDn}" stroke-width="4"/><line x1="18" y1="5" x2="30" y2="5" stroke="${pal.card_bg}" stroke-width="1.8"/></svg></span><span><strong>Hollow soft red/blue rail:</strong> Enz-pred (L1) inferred dir. (≤72% sat., lower conf.; never metabolite-only)</span></div>
        <div class="lg-row"><span class="lg-swatch"><svg width="32" height="10"><line x1="1" y1="5" x2="31" y2="5" stroke="${pal.edge_grey}" stroke-width="3"/><text x="7" y="8" font-size="9" fill="${pal.card_bg}">&lt;</text><text x="22" y="8" font-size="9" fill="${pal.card_bg}">&gt;</text></svg></span><span><strong>Opposing &lt; &gt; grey rail:</strong> Reversible reaction without single net direction (or multi-tracer conflict)</span></div>
        <div class="lg-row"><span class="lg-swatch"><svg width="32" height="10"><line x1="2" y1="5" x2="30" y2="5" stroke="${pal.coupled_track}" stroke-width="8" stroke-linecap="round"/><line x1="2" y1="5" x2="30" y2="5" stroke="${pal.up}" stroke-width="3"/></svg></span><span><strong>Linked flux (Shared 1 DOF — NOT extra confidence):</strong> Per-fit |r| ≥ 1−1e-6 across fits (mass balance); 1 hypothesis in grouped BH-FDR</span></div>
        <div class="lg-row"><span class="lg-swatch"><span style="display:inline-block;width:16px;height:16px;border-radius:50%;background:${pal.sign_change_bg};border:1px solid ${pal.sign_change_edge};color:${pal.sign_change_ink};text-align:center;line-height:14px;font-size:10px;font-weight:700;">&#8644;</span></span><span><strong>Amber &#8644; badge:</strong> Net flux reverses direction between experimental &amp; control groups</span></div>
        <div class="lg-row"><span class="lg-swatch"><svg width="32" height="12"><line x1="2" y1="6" x2="12" y2="6" stroke="${pal.ink2}" stroke-width="2" stroke-linecap="round"/><line x1="19" y1="6" x2="24" y2="6" stroke="${pal.muted}" stroke-width="1.8" stroke-linecap="round"/><line x1="26" y1="6" x2="31" y2="6" stroke="${pal.muted}" stroke-width="1.8" stroke-linecap="round"/></svg></span><span><strong>Cell states (solid vs dashed):</strong> Solid <code>—</code> = equivalence-tested <em>Unchanged</em>; dashed <code>--</code> = <em>Undetermined</em> (weak / unstable / insufficient evidence). Flux: <em>Sig.</em> = ≥2 tracers &amp; FDR&lt;0.05; <em>Suggestive</em> = |Hedges' g|≥1 &amp; LOO-stable</span></div>
        <div class="lg-row"><span class="lg-swatch"><i class="sc-dot" style="background:${pal.numerator};"></i><i class="sc-dot" style="background:${pal.denominator};"></i></span><span><strong>Per-animal scatter:</strong> ${esc(DATA.meta.num_label)} (Exp) vs ${esc(DATA.meta.den_label)} (Ctrl) · Grey = Other · ◇ = bound-pinned</span></div>
      </div>
      <div style="margin-top:7px;padding-top:6px;border-top:1px solid var(--hairline);font-size:10.5px;color:var(--muted);line-height:1.4;">
        <strong>Pooling &amp; Isotope Rule:</strong> Unconstrained (poorly_constrained / bound_pinned) and low-n tracers are gated out before inverse-variance fixed-effects pooling (w = 1/[(s²<sub>animal</sub> + f̄<sub>fit</sub>)/n], stratified permutation p). Supported isotopes: U-13C, U-15N, and U-13C-15N (kept in separate strata).
      </div>`;
    document.getElementById("btn-close-lg").onclick = () => el.classList.remove("open");
  }

  /* ---------- Region Zoom & SVG/PNG Export ---------- */
  function zoomToRegion(rid) {
    const rg = DATA.regions.find(r => r.id === rid);
    if (!rg || !rg.box) return;
    const [x0, y0, x1, y1] = rg.box;
    const pad = 3.0;
    animateViewBox([x0 - pad, y0 - pad, (x1 - x0) + 2 * pad, (y1 - y0) + 2 * pad], 300);
  }

  function exportSvgFile() {
    const clone = svg.cloneNode(true);
    clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
    const blob = new Blob([new XMLSerializer().serializeToString(clone)], { type: "image/svg+xml;charset=utf-8" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = `pathway_map_${DATA.meta.tissue || "atlas"}_${DATA.meta.contrast || "view"}_${state.theme}.svg`;
    a.click();
    URL.revokeObjectURL(a.href);
  }

  function exportPngFile() {
    const clone = svg.cloneNode(true);
    clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
    const [,, vw, vh] = state.vb;
    const scale = Math.min(2400 / vw, 1800 / vh);
    const cw = Math.round(vw * scale), ch = Math.round(vh * scale);
    clone.setAttribute("width", cw);
    clone.setAttribute("height", ch);
    const url = URL.createObjectURL(new Blob([new XMLSerializer().serializeToString(clone)], { type: "image/svg+xml;charset=utf-8" }));
    const img = new Image();
    img.onload = () => {
      const canvas = document.createElement("canvas");
      canvas.width = cw; canvas.height = ch;
      const ctx = canvas.getContext("2d");
      ctx.fillStyle = DATA.palettes[state.theme].surface;
      ctx.fillRect(0, 0, cw, ch);
      ctx.drawImage(img, 0, 0, cw, ch);
      URL.revokeObjectURL(url);
      const a = document.createElement("a");
      a.href = canvas.toDataURL("image/png");
      a.download = `pathway_map_${DATA.meta.tissue || "atlas"}_${DATA.meta.contrast || "view"}_${state.theme}.png`;
      a.click();
    };
    img.src = url;
  }

  /* ---------- Wire All DOM & SVG Event Handlers ---------- */
  svg.querySelectorAll(".ow-edge, .ow-strip").forEach(el => {
    const fid = el.getAttribute("data-rxn");
    const fg = el.getAttribute("data-flux-group");
    el.addEventListener("mouseenter", e => {
      showElementTooltip(e, "rxn:" + fid);
      if (fg) highlightFluxGroup(fg, true);
    });
    el.addEventListener("mousemove", moveTooltip);
    el.addEventListener("mouseleave", () => {
      hideTooltip();
      if (!state.selectedKey) highlightFluxGroup("", false);
    });
    el.addEventListener("click", e => { e.stopPropagation(); selectElement("rxn:" + fid, false); });
  });
  svg.querySelectorAll(".ow-agg").forEach(el => {
    const aid = el.getAttribute("data-agg");
    el.addEventListener("mouseenter", e => showElementTooltip(e, "agg:" + aid));
    el.addEventListener("mousemove", moveTooltip);
    el.addEventListener("mouseleave", hideTooltip);
    el.addEventListener("click", e => { e.stopPropagation(); selectElement("agg:" + aid, false); });
  });
  svg.querySelectorAll(".ow-node").forEach(el => {
    const sp = el.getAttribute("data-node");
    el.addEventListener("mouseenter", e => showElementTooltip(e, "node:" + sp));
    el.addEventListener("mousemove", moveTooltip);
    el.addEventListener("mouseleave", hideTooltip);
    el.addEventListener("click", e => { e.stopPropagation(); selectElement("node:" + sp, false); });
  });
  svg.querySelectorAll(".ow-region-title").forEach(el => {
    el.addEventListener("click", e => { e.stopPropagation(); zoomToRegion(el.getAttribute("data-region")); });
  });

  // Toolbar bindings
  const themeSel = document.getElementById("sel-theme");
  themeSel.value = state.theme;
  themeSel.onchange = () => applyTheme(themeSel.value);

  const subWrap = document.getElementById("tb-submaps");
  DATA.submaps.forEach((sm, idx) => {
    const b = document.createElement("button");
    b.className = "btn" + (idx === 0 ? " active" : "");
    b.setAttribute("data-submap", sm.id);
    b.textContent = sm.label;
    b.title = sm.full_title || sm.label;
    b.onclick = () => {
      subWrap.querySelectorAll(".btn").forEach(x => x.classList.remove("active"));
      b.classList.add("active");
      animateViewBox(sm.box, 300);
    };
    subWrap.appendChild(b);
  });

  document.querySelectorAll("[data-filter]").forEach(btn => {
    btn.onclick = () => applyFilter(btn.getAttribute("data-filter"));
  });

  const omicsBtn = document.getElementById("btn-omics-color");
  function syncOmicsBtn() {
    omicsBtn.textContent = state.omicsColor === "tinted" ? "Enz-pred: Tinted" : "Enz-pred: Grey";
    omicsBtn.classList.toggle("active", state.omicsColor === "tinted");
  }
  syncOmicsBtn();
  omicsBtn.onclick = () => {
    state.omicsColor = state.omicsColor === "tinted" ? "grey" : "tinted";
    syncOmicsBtn();
    recolorEdges(DATA.palettes[state.theme]);
  };

  document.getElementById("btn-legend").onclick = () => document.getElementById("ow-legend-float").classList.toggle("open");
  document.getElementById("btn-drawer-toggle").onclick = () => drawer.classList.toggle("collapsed");
  document.getElementById("btn-export-svg").onclick = exportSvgFile;
  document.getElementById("btn-export-png").onclick = exportPngFile;
  document.getElementById("btn-z-in").onclick = () => { const r = svg.getBoundingClientRect(); zoomAt(0.78, r.left + r.width / 2, r.top + r.height / 2); };
  document.getElementById("btn-z-out").onclick = () => { const r = svg.getBoundingClientRect(); zoomAt(1.28, r.left + r.width / 2, r.top + r.height / 2); };
  document.getElementById("btn-z-reset").onclick = () => animateViewBox(state.baseVb, 260);

  // Instant search across reactions, genes, proteins, metabolites, EC
  const sInput = document.getElementById("ow-search");
  const sMenu = document.getElementById("ow-search-results");
  sInput.addEventListener("input", () => {
    const q = sInput.value.trim().toLowerCase();
    if (!q) { sMenu.classList.remove("open"); return; }
    const hits = DATA.search.filter(item => item.tokens.includes(q)).slice(0, 14);
    if (!hits.length) {
      sMenu.innerHTML = `<div class="sr-item"><span class="sr-sub">No match for "${esc(q)}"</span></div>`;
    } else {
      sMenu.innerHTML = hits.map(h => `
        <div class="sr-item" data-hit="${esc(h.key)}">
          <span class="sr-title">${esc(h.title)}</span>
          <span class="sr-sub">${esc(h.kind)} · ${esc(h.sub)}</span>
        </div>`).join("");
      sMenu.querySelectorAll("[data-hit]").forEach(item => {
        item.onclick = () => {
          sMenu.classList.remove("open");
          sInput.value = "";
          selectElement(item.getAttribute("data-hit"), true);
        };
      });
    }
    sMenu.classList.add("open");
  });
  document.addEventListener("click", e => { if (!e.target.closest(".search-wrap")) sMenu.classList.remove("open"); });
  document.addEventListener("keydown", e => {
    if (e.key === "/" && document.activeElement !== sInput) { e.preventDefault(); sInput.focus(); }
    else if (e.key === "Escape") { sMenu.classList.remove("open"); document.getElementById("ow-legend-float").classList.remove("open"); }
    else if (e.key === "0" && document.activeElement !== sInput) animateViewBox(state.baseVb, 260);
  });

  // Boot with initial theme + cross-platform zero-overflow label fit + URL deep-link params
  const params = new URLSearchParams(window.location.search);
  const qTheme = params.get("theme");
  if (qTheme && DATA.palettes[qTheme]) {
    state.theme = qTheme;
    themeSel.value = qTheme;
  }
  applyTheme(state.theme);
  fitSvgLabels();
  if (document.fonts && document.fonts.ready) {
    document.fonts.ready.then(fitSvgLabels);
  }
  if (params.get("drawer") === "0") {
    drawer.classList.add("collapsed");
  }
  const qFocus = params.get("focus");
  if (qFocus) {
    const sm = DATA.submaps.find(s => s.id === qFocus);
    if (sm) {
      setViewBox(sm.box);
      subWrap.querySelectorAll(".btn").forEach(x => x.classList.toggle("active", x.getAttribute("data-submap") === sm.id));
    } else {
      const rg = DATA.regions.find(r => r.id === qFocus);
      if (rg && rg.box) setViewBox([rg.box[0] - 2.5, rg.box[1] - 2.5, (rg.box[2] - rg.box[0]) + 5.0, (rg.box[3] - rg.box[1]) + 5.0]);
    }
  }
  const qZoom = params.get("zoom");
  if (qZoom) {
    const nums = qZoom.split(",").map(Number);
    if (nums.length === 4 && nums.every(n => isFinite(n))) setViewBox(nums);
  }
  const qSelect = params.get("select");
  if (qSelect && DATA.elements[qSelect]) {
    selectElement(qSelect, false);
  }
})();
"""


def write_pathway_map_html(
    view: pd.DataFrame,
    layout: Layout,
    regions: pd.DataFrame,
    title: str,
    caption: list[str],
    output_path: Path | str,
    *,
    scene: dict[str, Any] | None = None,
    extra_fusion: pd.DataFrame | None = None,
    signed_modules: pd.DataFrame | None = None,
    top_rank: int = 10,
    mode: str = "atlas",
    style: str = "metro_atlas",
    color_mode: str = "semantic",
    strip_mode: str = "omics",
    omics_edge_color: str = "tinted",
    focus_regions: list[str] | None = None,
    include_reactions: list[str] | None = None,
    exclude_reactions: list[str] | None = None,
    members: pd.DataFrame | None = None,
    samples: pd.DataFrame | None = None,
    sample_contrasts: pd.DataFrame | None = None,
    tissue: str = "",
    contrast_id: str = "",
    source: str = "",
) -> Path:
    """生成单文件离线可交互通路探索器 HTML（方案 A：原生语义化 SVG + Metro Atlas 默认学术皮肤 + 逐只动物散点小图）。"""
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    pal_name = style if style in PALETTES else "metro_atlas"
    palette = get_palette(pal_name)

    if scene is None:
        scene = build_pathway_map_scene(
            view, layout,
            mode=mode, style=pal_name, top_rank=top_rank,
            focus_regions=focus_regions, include_reactions=include_reactions, exclude_reactions=exclude_reactions,
        )

    svg_markup, vb = _build_svg(
        scene, layout,
        palette=palette, preset=pal_name, color_mode=color_mode,
        strip_mode=strip_mode, omics_edge_color=omics_edge_color,
        top_rank=top_rank, mode=mode,
    )
    payload = _extract_payload(
        view, layout, regions, scene,
        title=title, caption=caption, style=pal_name,
        color_mode=color_mode, strip_mode=strip_mode, omics_edge_color=omics_edge_color,
        vb=vb, extra_fusion=extra_fusion, signed_modules=signed_modules,
        members=members, samples=samples, sample_contrasts=sample_contrasts,
        tissue=tissue, contrast_id=contrast_id, source=source,
    )
    # Safe JSON embedding (escape </script> to prevent premature tag close)
    json_str = json.dumps(payload, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")

    c = payload["meta"]["counts"]
    sub_line = caption[0] if caption else ""
    cal_line = caption[1] if len(caption) > 1 else ""
    lay_line = caption[2] if len(caption) > 2 else ""

    doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>{_esc(title)}</title>
<style>{_CSS}</style>
</head>
<body data-theme="{_esc(payload['meta']['initial_style'])}">
<div id="app">
  <header class="masthead">
    <div class="red-tag-bar"><span class="red-tag" aria-hidden="true"></span></div>
    <div class="mast-row">
      <div class="mast-title-wrap">
        <h1 class="mast-title">{_esc(title)}</h1>
        <div class="mast-sub">
          <span>{_esc(sub_line)}</span>
        </div>
      </div>
      <div class="mast-pills">
        <span class="stat-pill l3">&#9679;&#9679;&#9679; Multi-tracer: {c['L3']}</span>
        <span class="stat-pill l2">&#9679;&#9679;&#9675; Flux: {c['L2']}{f" ({c['indep_flux_chains']} indep.)" if (c['L2'] + c['L3']) > c['indep_flux_chains'] else ""}</span>
        <span class="stat-pill l1">&#9679;&#9675;&#9675; Enz-pred: {c['L1']}</span>
        <span class="stat-pill link" title="Stoichiometrically locked chains sharing 1 DOF in grouped BH-FDR (not extra confidence)">&#128279; Linked: {c['linked_groups']} (1 DOF)</span>
      </div>
    </div>
  </header>

  <nav class="toolbar" aria-label="Pathway map controls">
    <div class="tb-group">
      <span class="tb-label">Style</span>
      <select id="sel-theme" class="ctrl-select" title="Switch visual theme live">
        <option value="metro_atlas">Metro Atlas (Default)</option>
        <option value="economist">Economist (Editorial)</option>
        <option value="nature_editorial">Nature Editorial</option>
        <option value="slide_dark">Slide Dark</option>
      </select>
    </div>
    <span class="tb-sep"></span>
    <div class="tb-group" id="tb-submaps"><span class="tb-label">Focus</span></div>
    <span class="tb-sep"></span>
    <div class="tb-group">
      <span class="tb-label">Filter</span>
      <button class="btn active" data-filter="all">All</button>
      <button class="btn" data-filter="l2">Flux (L2+)</button>
      <button class="btn" data-filter="l1">Enz-pred (L1)</button>
      <button class="btn" data-filter="linked">Linked</button>
      <button class="btn" data-filter="reversal">&#8644; Flip</button>
    </div>
    <span class="tb-sep"></span>
    <div class="tb-group">
      <button class="btn red-accent" id="btn-omics-color" title="Toggle soft red/blue vs neutral grey for Enz-pred (L1) edges">Enz-pred: Tinted</button>
    </div>
    <span class="tb-sep"></span>
    <div class="search-wrap">
      <input id="ow-search" class="search-input" type="search" placeholder="Search reaction, gene, metabolite, EC (/)"/>
      <div id="ow-search-results" class="search-results"></div>
    </div>
    <div class="tb-group" style="margin-left:auto;">
      <button class="btn" id="btn-legend">Legend</button>
      <button class="btn" id="btn-export-svg">SVG</button>
      <button class="btn" id="btn-export-png">PNG</button>
      <button class="btn" id="btn-drawer-toggle">Dossier</button>
    </div>
  </nav>

  <main class="workspace">
    <section class="stage" id="ow-stage">
      {svg_markup}
      <div class="zoom-ctrls">
        <button class="btn" id="btn-z-in" title="Zoom in">+</button>
        <button class="btn" id="btn-z-out" title="Zoom out">&minus;</button>
        <button class="btn" id="btn-z-reset" title="Reset full view (0)">&#8634;</button>
      </div>
      <div class="legend-float" id="ow-legend-float"></div>
    </section>
    <aside class="drawer" id="ow-drawer">
      <div class="drawer-head" id="ow-drawer-head"></div>
      <div class="drawer-body" id="ow-drawer-body"></div>
    </aside>
  </main>

  <footer class="footer">
    <span>{_esc(cal_line)}</span>
    <span>{_esc(lay_line)}</span>
  </footer>
</div>
<div class="tooltip" id="ow-tooltip"></div>
<script id="ow-payload" type="application/json">{json_str}</script>
<script>{_JS}</script>
</body>
</html>
"""
    out.write_text(doc, encoding="utf-8")
    return out


__all__ = ["write_pathway_map_html"]
