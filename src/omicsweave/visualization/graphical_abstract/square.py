#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cell Press 官方 1:1 方形 Graphical Abstract 渲染策略模块（零业务耦合，可跨项目复用）。

严格遵循 Cell Press Graphical Abstract 规范：
- 物理与像素尺寸：精确 **1200 × 1200 px @ 300 dpi**（4.0 × 4.0 in，``bbox_inches=None`` 无裁切漂移）；
- 单面板叙事流：逻辑画布 ``100 × 100``，采用上→下、左→右清晰生理流向（循环/胞质底物 → 线粒体核心代谢结局）；
- 字体硬底线：全图统一无衬线 **Arial / Helvetica**，所有可见文字严格处于 **12–16 pt** 区间
  （标题 14–15 pt Bold，枢纽与图例 12–13 pt，彻底杜绝缩略图糊字）；
- 沿用 ``metro_atlas`` Beck 地铁视觉语法（八向圆角轨道、内嵌 chevron、实心 L2+ vs 空心 L1 双轨、
  TCA 圆角跑道闭环与双层换乘圆站），同时兼容 ``economist`` / ``nature_editorial`` / ``editorial_ink`` / ``pixel_abstract``。
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from matplotlib.patches import FancyBboxPatch

from ..glyphs import GLYPH_STYLES, draw_glyph
from ..metro import (
    draw_metro_coupled_halo,
    draw_metro_loop,
    draw_metro_station,
    draw_metro_track,
    draw_metro_zone,
    fillet_polyline,
    route_octilinear,
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
from ..primitives import draw_economist_masthead, estimate_text_width
from .spec import AbstractFlow, AbstractHub, AbstractSpec

# ---------------------------------------------------------------------------
# 方形画布核心槽位（100 × 100 坐标系，保证 12 pt Arial 下所有铭牌零碰撞）
# ---------------------------------------------------------------------------

# 3 列 × 2 行宽间距地铁网：列中心 x = 17, 50, 84（配合简明铭牌保证 12 pt Bold 零接触）
_SQUARE_SLOTS: dict[str, tuple[float, float]] = {
    "top_left": (17.0, 29.0),
    "top_mid": (50.0, 29.0),
    "top_right": (84.0, 29.0),
    "bot_left": (17.0, 60.0),
    "bot_mid": (51.0, 60.0),
    "bot_right": (85.0, 60.0),
}

# 标准代谢枢纽别名 → (优先槽位, 简明 Cell Press 12 pt 铭牌，控制在 ≤10 字符确保列间净空 ≥8 units)
_CANONICAL_SLOT_RULES: dict[str, tuple[str, str]] = {
    "blood_glc": ("top_left", "Blood Glc"),
    "glycerolipid": ("top_left", "TAG / G3P"),
    "glycolysis_gng": ("top_mid", "Glycolysis"),
    "glycolysis": ("top_mid", "Glycolysis"),
    "ppp": ("top_mid", "PPP"),
    "fatty_acid": ("top_right", "Fatty Acid"),
    "pyruvate_hub": ("bot_left", "Pyruvate"),
    "pyruvate": ("bot_left", "Pyruvate"),
    "amino_acid": ("bot_left", "Amino Acid"),
    "tca_cycle": ("bot_mid", "TCA Cycle"),
    "tca": ("bot_mid", "TCA Cycle"),
    "ketone": ("bot_right", "Ketones"),
    "oxphos": ("bot_right", "OXPHOS"),
}

_TIER_PRIORITY: dict[str, int] = {"L3": 5, "L2": 4, "L1": 3, "omics": 2, "": 1}
_CORE_BACKBONE: set[str] = {"blood_glc", "glycolysis_gng", "glycolysis", "fatty_acid", "pyruvate_hub", "pyruvate", "tca_cycle", "tca"}

_SQUARE_SURFACES: dict[str, str] = {
    "metro_atlas": "#FFFFFF",
    "nature_editorial": "#FFFFFF",
    "economist": ECONOMIST.surface,
    "editorial_ink": "#FAF8F2",
    "pixel_abstract": "#F4EFE4",
}


@dataclass(frozen=True)
class _SquareHub:
    hub: AbstractHub
    xy: tuple[float, float]
    short_label: str
    is_loop: bool
    slot: str


def _short_square_label(hub: AbstractHub) -> str:
    if hub.id in _CANONICAL_SLOT_RULES:
        return _CANONICAL_SLOT_RULES[hub.id][1]
    lbl = str(hub.label).strip()
    if len(lbl) <= 13:
        return lbl
    # 优先在分隔符处截短，避免生硬截断单词
    for sep in (" & ", " / ", " - ", " — ", ": "):
        if sep in lbl:
            first = lbl.split(sep, 1)[0].strip()
            if 3 <= len(first) <= 13:
                return first
    words = lbl.split()
    out = ""
    for w in words:
        cand = f"{out} {w}".strip()
        if len(cand) <= 13:
            out = cand
        else:
            break
    return out or lbl[:12]


def _hub_score(hub: AbstractHub) -> tuple[int, int, float, int]:
    return (
        _TIER_PRIORITY.get(hub.evidence_tier, 1),
        1 if hub.direction != 0 else 0,
        float(abs(hub.magnitude)),
        1 if hub.id in _CORE_BACKBONE else 0,
    )


def _is_loop_hub(hub: AbstractHub) -> bool:
    return hub.icon in ("tca", "ppp") or hub.id in ("tca_cycle", "tca")


def _select_and_place_square_hubs(spec: AbstractSpec) -> list[_SquareHub]:
    """将任意 AbstractSpec（自动聚合或跨课题 YAML）映射到无重叠 100×100 Cell Press 方形网格（≤ 6 核心枢纽）。"""
    canonical_hits = sum(1 for h in spec.hubs if h.id in _CANONICAL_SLOT_RULES)
    placed: list[_SquareHub] = []

    if canonical_hits >= max(2, len(spec.hubs) // 2):
        # 按槽位竞争选出每个生理槽位证据最强/最核心的枢纽，空槽再从剩余高证据枢纽补位
        by_slot: dict[str, list[AbstractHub]] = {k: [] for k in _SQUARE_SLOTS}
        leftovers: list[AbstractHub] = []
        for hub in spec.hubs:
            if hub.id in _CANONICAL_SLOT_RULES:
                slot = _CANONICAL_SLOT_RULES[hub.id][0]
                by_slot[slot].append(hub)
            else:
                leftovers.append(hub)

        chosen: dict[str, AbstractHub] = {}
        overflow: list[AbstractHub] = list(leftovers)
        for slot, cands in by_slot.items():
            if not cands:
                continue
            # 同一槽位中：主干核心（如 Glycolysis / Pyruvate / TCA）在同级时优先，证据更强（L2+）胜出
            cands_sorted = sorted(
                cands,
                key=lambda h: (
                    1 if (h.id in ("tca_cycle", "tca", "glycolysis_gng", "glycolysis") and h.evidence_tier != "") else 0,
                    *_hub_score(h),
                ),
                reverse=True,
            )
            chosen[slot] = cands_sorted[0]
            overflow.extend(cands_sorted[1:])

        # 若有空槽位且存在活跃剩余枢纽，则择优补入空槽
        empty_slots = [s for s in _SQUARE_SLOTS if s not in chosen]
        overflow_sorted = sorted(overflow, key=_hub_score, reverse=True)
        for slot, hub in zip(empty_slots, overflow_sorted):
            if _TIER_PRIORITY.get(hub.evidence_tier, 1) >= 2 or hub.direction != 0:
                chosen[slot] = hub

        for slot in _SQUARE_SLOTS:
            hub = chosen.get(slot)
            if hub is None:
                continue
            placed.append(_SquareHub(
                hub=hub,
                xy=_SQUARE_SLOTS[slot],
                short_label=_short_square_label(hub),
                is_loop=(slot == "bot_mid" and _is_loop_hub(hub)),
                slot=slot,
            ))
        return placed

    # 通用跨课题 YAML 回退：按原始 (y, x) 空间秩序排入 3×2 方形网格（最多 6 个核心枢纽）
    top_hubs = sorted(spec.hubs, key=_hub_score, reverse=True)[:6]
    # 按原始 y 分为上下两排，排内按 x 从左到右
    ys = sorted(h.xy[1] for h in top_hubs)
    mid_y = ys[len(ys) // 2] if ys else 40.0
    row1 = sorted([h for h in top_hubs if h.xy[1] < mid_y], key=lambda h: h.xy[0])
    row2 = sorted([h for h in top_hubs if h.xy[1] >= mid_y], key=lambda h: h.xy[0])
    while len(row1) > 3:
        row2.insert(0, row1.pop())
    while len(row2) > 3:
        row1.append(row2.pop(0))

    def _row_xs(n: int) -> list[float]:
        if n == 1:
            return [52.0]
        if n == 2:
            return [28.0, 76.0]
        return [18.0, 52.0, 84.0]

    for row_idx, (row, y_val) in enumerate(((row1, 28.0), (row2, 60.0))):
        for col_idx, (hub, x_val) in enumerate(zip(row, _row_xs(len(row)))):
            slot_name = f"{'top' if row_idx == 0 else 'bot'}_{col_idx}"
            placed.append(_SquareHub(
                hub=hub,
                xy=(x_val, y_val),
                short_label=_short_square_label(hub),
                is_loop=(_is_loop_hub(hub) and len(row) <= 3 and col_idx == (1 if len(row) == 3 else 0)),
                slot=slot_name,
            ))
    return placed


def _square_title(spec: AbstractSpec, unit_in: float, *, max_w: float | None = None) -> tuple[str, float]:
    """裁剪并适配 Cell Press 方形顶栏标题，确保字号严格处于 12.0–15.0 pt 且零越界、不断词。"""
    raw = str(spec.title or "Metabolic Rewiring").strip()
    for suffix in (" — Graphical Abstract", " - Graphical Abstract", ": Graphical Abstract"):
        if raw.endswith(suffix):
            raw = raw[:-len(suffix)].strip()
    budget = float(max_w) if max_w is not None else (64.0 if spec.contrast else 90.0)
    # 优先完整保留标题，其次优雅凝练 "<Tissue> Metabolic Rewiring" → "<Tissue> Metabolism"
    candidates = [raw]
    if raw.endswith(" Metabolic Rewiring"):
        candidates.append(raw[: -len(" Metabolic Rewiring")] + " Metabolism")
    for cand in candidates:
        for fs in (14.5, 13.5, 12.5, 12.0):
            if estimate_text_width(cand, fs, unit_in, bold=True) <= budget:
                return cand, fs
    words = raw.split()
    trimmed = ""
    for w in words:
        cand = f"{trimmed} {w}".strip()
        if estimate_text_width(cand, 12.0, unit_in, bold=True) <= budget:
            trimmed = cand
        else:
            break
    return (trimmed or raw[:18]), 12.0


def _palette_for_style(style: str):
    if style == "economist":
        return ECONOMIST
    if style == "metro_atlas":
        return METRO_ATLAS
    return NATURE_EDITORIAL


def _hub_color(hub: AbstractHub, pal, color_mode: str) -> str:
    if color_mode == "pathway":
        return pathway_family_color(hub.id)
    if hub.direction > 0:
        return pal.up
    if hub.direction < 0:
        return pal.down
    return pal.omics_grey


def _flow_color(flow: AbstractFlow, src: AbstractHub, dst: AbstractHub, pal, color_mode: str) -> str:
    if color_mode == "pathway":
        return pathway_family_color(dst.id if dst.direction else src.id)
    if flow.direction > 0:
        return pal.up
    if flow.direction < 0:
        return pal.down
    return pal.omics_grey


def _loop_anchor(cx: float, cy: float, lw: float, lh: float, other: tuple[float, float]) -> tuple[float, float]:
    ox, oy = other
    dx, dy = ox - cx, oy - cy
    if abs(dy) <= lh * 0.45:
        sx = 1.0 if dx >= 0 else -1.0
        return (cx + sx * lw * 0.5, max(cy - lh * 0.25, min(cy + lh * 0.25, oy)))
    if abs(dx) <= lw * 0.28:
        sy = 1.0 if dy >= 0 else -1.0
        return (cx + dx * 0.35, cy + sy * lh * 0.5)
    sx = 1.0 if dx >= 0 else -1.0
    sy = 1.0 if dy >= 0 else -1.0
    return (cx + sx * lw * 0.26, cy + sy * lh * 0.5)


def render_cell_press_square(
    ax,
    spec: AbstractSpec,
    *,
    style: str = "metro_atlas",
    color_mode: str = "semantic",
    unit_in: float = 0.04,
) -> None:
    """在 ``100 × 100`` 画布（4.0 × 4.0 in @ 300 dpi = 1200 × 1200 px）渲染 Cell Press 方形 GA。

    全图所有文字严格处于 **12.0–15.0 pt** 无衬线 Arial/Helvetica，保留 ``metro_atlas`` 地铁换乘美学。
    """
    if color_mode not in VALID_COLOR_MODES:
        raise ValueError(f"color_mode 应为 {' / '.join(VALID_COLOR_MODES)} 之一，实际为 {color_mode!r}")

    pal = _palette_for_style(style)
    surface = _SQUARE_SURFACES.get(style, pal.surface)
    is_econ = style == "economist"
    use_glyph_cards = style in ("editorial_ink", "pixel_abstract")

    ax.set_facecolor(surface)
    ax.set_xlim(0.0, 100.0)
    ax.set_ylim(100.0, 0.0)
    ax.set_aspect("equal")
    ax.axis("off")

    # ── 1. 顶部极简标题栏（12–14.5 pt Arial，无冗长副标题）────────────────
    if is_econ:
        draw_economist_masthead(
            ax, 4.0, 1.8, 92.0,
            tag_w=10.5, tag_h=2.0, color=ECONOMIST_RED, rule_color="#B8B8AF", y_up=False, zorder=8.0,
        )
        title_y = 8.4
        rule_y = 14.2
    else:
        title_y = 7.4
        rule_y = 13.4

    c_width = 0.0
    if spec.contrast:
        c_txt = str(spec.contrast).strip()
        c_fs = 12.0
        c_width = estimate_text_width(c_txt, c_fs, unit_in, bold=True) + 4.4
        cx0 = 96.0 - c_width
        ax.add_patch(FancyBboxPatch(
            (cx0, title_y - 3.3), c_width, 6.6,
            boxstyle=f"round,pad=0,rounding_size={0.5 if is_econ else 2.2}",
            facecolor=pal.card_bg, edgecolor=pal.card_edge, lw=0.75, zorder=5.5,
        ))
        ax.text(
            cx0 + c_width * 0.5, title_y, c_txt,
            fontsize=c_fs, fontweight="bold", color=pal.ink2,
            ha="center", va="center", zorder=5.8,
        )

    title_budget = (92.0 - c_width - 2.6) if c_width > 0 else 92.0
    title_txt, title_fs = _square_title(spec, unit_in, max_w=title_budget)
    ax.text(
        4.0, title_y, title_txt,
        fontsize=title_fs, fontweight="bold", color=pal.ink,
        ha="left", va="center", zorder=6.0,
    )

    ax.plot([4.0, 96.0], [rule_y, rule_y], color=pal.card_edge, lw=0.80, zorder=1.5)

    # ── 2. 选择并排布核心枢纽（≤ 6 枢纽，确保 12 pt 零重叠）────────────────
    sq_hubs = _select_and_place_square_hubs(spec)
    sq_by_id = {sh.hub.id: sh for sh in sq_hubs}
    has_mito = any(sh.hub.compartment == "mitochondria" or sh.slot.startswith("bot") for sh in sq_hubs)
    has_blood = any(sh.hub.compartment == "blood" and sh.slot == "top_left" for sh in sq_hubs)

    # ── 3. 低存在感细胞器区室底板（12 pt 角落标签，绝不与站台打架）─────────
    if has_blood:
        draw_metro_zone(
            ax, (4.0, 16.5, 33.0, 41.0), "",
            face="#EEF2F6" if not is_econ else "#E3E2D5",
            edge="#CBD5E1" if not is_econ else pal.card_edge,
            lw=0.70, rounding=3.2, zorder=0.4,
        )
    if has_mito:
        mito_face = "#F1F3F7" if not is_econ else "#E5E4D7"
        mito_edge = "#CFD6E0" if not is_econ else pal.card_edge
        mito_inner = "#C3CCD8" if not is_econ else "#B8B8AF"
        draw_metro_zone(
            ax, (4.0, 44.5, 96.0, 83.0), "",
            face=mito_face, edge=mito_edge, lw=0.85,
            inner_edge=mito_inner, rounding=4.2, zorder=0.4,
        )
        ax.text(
            93.5, 80.0, "Mitochondria",
            fontsize=12.0, fontweight="bold", color=pal.muted,
            ha="right", va="bottom", zorder=0.8,
        )

    # ── 4. 地铁八向轨道（实心 = Flux L2+；空心白芯双轨 = Enz-pred L1）──────
    loop_w, loop_h = 26.5, 19.0
    station_r = 3.2
    loop_junctions: list[tuple[tuple[float, float], str]] = []

    for fl in spec.flows:
        if fl.source not in sq_by_id or fl.target not in sq_by_id:
            continue
        sh, dh = sq_by_id[fl.source], sq_by_id[fl.target]
        p0 = _loop_anchor(*sh.xy, loop_w, loop_h, dh.xy) if sh.is_loop else sh.xy
        p1 = _loop_anchor(*dh.xy, loop_w, loop_h, sh.xy) if dh.is_loop else dh.xy

        prefer = "hv_first" if fl.curve >= 0 else "diag_first"
        if abs(p1[0] - p0[0]) > 8.0 and abs(p1[1] - p0[1]) > 8.0:
            prefer = "hv_first"
        raw_poly = route_octilinear(p0, p1, bend=0.0, prefer=prefer)
        smooth_poly = fillet_polyline(raw_poly, radius=3.5, arc_steps=8)

        col = _flow_color(fl, sh.hub, dh.hub, pal, color_mode)
        strong = fl.evidence_tier in ("L3", "L2")
        is_l1 = fl.evidence_tier == "L1"
        track_w = 5.2 if strong else (4.4 if is_l1 else 2.8)
        trim0 = 0.0 if sh.is_loop else (station_r + 0.35)
        trim1 = 0.0 if dh.is_loop else (station_r + 0.35)

        if fl.coupled:
            draw_metro_coupled_halo(
                ax, [smooth_poly],
                color=pal.coupled_track, edge_color=pal.coupled_edge, width=10.5, zorder=1.6,
            )

        draw_metro_track(
            ax, raw_poly, col,
            width=track_w, hollow=is_l1, surface=surface,
            chevrons=True, chevron_spacing=12.0, chevron_direction=1,
            unit_in=unit_in, trim0=trim0, trim1=trim1,
            corner_radius=3.5, halo=True, zorder=2.8 if strong else 2.2,
        )
        if sh.is_loop:
            loop_junctions.append((p0, _hub_color(sh.hub, pal, color_mode)))
        if dh.is_loop:
            loop_junctions.append((p1, _hub_color(dh.hub, pal, color_mode)))

    # ── 5. 跑道闭环枢纽（如 TCA Cycle）与边缘换乘圆点──────────────────────
    for sh in sq_hubs:
        if not sh.is_loop:
            continue
        col = _hub_color(sh.hub, pal, color_mode)
        strong = sh.hub.evidence_tier in ("L3", "L2")
        hollow = sh.hub.evidence_tier == "L1"
        arrow = " ↑" if sh.hub.direction > 0 else (" ↓" if sh.hub.direction < 0 else "")
        sub_txt = f"Flux{arrow}" if strong else (f"Enz{arrow}" if hollow else ("Stable" if not arrow else f"Pool{arrow}"))
        sub_col = pal.up if sh.hub.direction > 0 else (pal.down if sh.hub.direction < 0 else pal.muted)

        draw_metro_loop(
            ax, sh.xy, loop_w, loop_h, col,
            track_width=5.4 if strong else 4.4,
            hollow=hollow,
            surface=surface,
            interior_tint=0.08 if strong else 0.04,
            clockwise=True,
            n_chevrons=4,
            label="",
            subtitle="",
            unit_in=unit_in,
            zorder=3.6,
        )
        ax.text(
            sh.xy[0], sh.xy[1] - 2.4, sh.short_label,
            fontsize=12.5, fontweight="bold", color=pal.ink if color_mode == "semantic" else col,
            ha="center", va="center", zorder=4.4,
        )
        ax.text(
            sh.xy[0], sh.xy[1] + 3.2, sub_txt,
            fontsize=12.0, fontweight="bold", color=sub_col,
            ha="center", va="center", zorder=4.4,
        )

    for jpt, lcol in loop_junctions:
        draw_metro_station(
            ax, jpt, 1.55,
            fill=tint_color(lcol, 0.25, base=surface),
            edge=lcol, lw=1.6, surface=surface, halo_pad=0.45, zorder=4.6,
        )

    # ── 6. 换乘圆站与 12 pt 无重叠铭牌──────────────────────────────────────
    for sh in sq_hubs:
        if sh.is_loop:
            continue
        x, y = sh.xy
        col = _hub_color(sh.hub, pal, color_mode)
        strong = sh.hub.evidence_tier in ("L3", "L2")
        fill = tint_color(col, 0.26 if strong else (0.16 if sh.hub.direction else 0.08), base=surface)

        if use_glyph_cards and style in GLYPH_STYLES:
            ax.add_patch(FancyBboxPatch(
                (x - 4.4, y - 4.4), 8.8, 8.8,
                boxstyle=f"round,pad=0,rounding_size={0.5 if style == 'pixel_abstract' else 2.2}",
                facecolor=fill, edgecolor=col, lw=2.0 if strong else 1.5, zorder=5.0,
            ))
            draw_glyph(
                ax, sh.hub.icon, (x, y),
                style=style, size=6.2, zorder=5.4,
            )
        else:
            draw_metro_station(
                ax, (x, y), station_r,
                fill=fill,
                edge=col,
                outer_ring=pal.ink if strong else None,
                lw=2.4 if strong else 1.9,
                surface=surface,
                halo_pad=0.65,
                zorder=5.0,
            )

        arrow = " ↑" if sh.hub.direction > 0 else (" ↓" if sh.hub.direction < 0 else "")
        lbl_txt = f"{sh.short_label}{arrow}"
        lbl_col = pal.ink if sh.hub.direction == 0 else (pal.up if sh.hub.direction > 0 else pal.down)
        # 顶排标签置于站台上方（y = 22.2），底排置于站台下方（y = 67.2），底色与所在区室无缝融合
        if sh.slot.startswith("top"):
            ly, va = y - 6.8, "center"
            lbl_bg = ("#EEF2F6" if not is_econ else "#E3E2D5") if (has_blood and sh.slot == "top_left") else surface
        else:
            ly, va = y + 7.2, "center"
            lbl_bg = ("#F1F3F7" if not is_econ else "#E5E4D7") if has_mito else surface
        # 保证左右边界留白 ≥ 2.5 逻辑单位
        est_w = estimate_text_width(lbl_txt, 12.0, unit_in, bold=True)
        lx = max(3.0 + est_w * 0.5, min(97.0 - est_w * 0.5, x))
        ax.text(
            lx, ly, lbl_txt,
            fontsize=12.0, fontweight="bold", color=lbl_col,
            ha="center", va=va, zorder=5.8,
            bbox=dict(boxstyle="round,pad=0.14", facecolor=lbl_bg, edgecolor="none", alpha=0.92),
        )

    # ── 7. 底部极简 1 行地铁图例条（严格 12.0 pt Arial，零微小脚注）─────────
    ax.plot([4.0, 96.0], [86.0, 86.0], color=pal.card_edge, lw=0.80, zorder=1.5)
    leg_y = 92.6
    if color_mode == "pathway":
        leg_items = [
            ("Flux", pal.ink2, False, 4.4),
            ("Enz-pred", pal.ink2, True, 4.0),
            ("Stable", pal.omics_grey, False, 2.6),
        ]
    else:
        leg_items = [
            ("Flux ↑", pal.up, False, 4.4),
            ("Flux ↓", pal.down, False, 4.4),
            ("Enz-pred", pal.ink2, True, 4.0),
        ]
    seg_w = 8.4
    gap = 4.6
    widths = [seg_w + 1.8 + estimate_text_width(lbl, 12.0, unit_in, bold=True) for lbl, _, _, _ in leg_items]
    total_w = sum(widths) + gap * (len(widths) - 1)
    cx = max(4.0, (100.0 - total_w) * 0.5)
    for (lbl, col, hollow, lw), item_w in zip(leg_items, widths):
        draw_metro_track(
            ax, [(cx, leg_y), (cx + seg_w, leg_y)], col,
            width=lw, hollow=hollow, surface=surface,
            chevrons=True, chevron_spacing=12.0, unit_in=unit_in,
            corner_radius=0.0, halo=False, zorder=5.2,
        )
        ax.text(
            cx + seg_w + 1.8, leg_y, lbl,
            fontsize=12.0, fontweight="bold", color=pal.ink2,
            ha="left", va="center", zorder=5.4,
        )
        cx += item_w + gap


__all__ = ["render_cell_press_square"]
