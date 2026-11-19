#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
图 E（G18）：分组织代谢通路总图，只画 :func:`omicsweave.integration.dossier.map_view.map_view` 的表。

编码（设计文档 ``omicsweave_pathway_map.md`` 第 4 节）：

- 反应：实线 = 该组织有拟合通量，虚线 = 只有组学，细点线 = 未评估。箭头在生理净方向上（两组变号时双向并标 ⇄）。
  边的颜色只表示 F 轴，换算到生理方向（红 = 通量增加、蓝 = 减少，色阶与图 A 相同）；F 轴不是 suggestive / changed 时为灰；
- 反应中心旁的四格小条为 F、C、R、M 四轴（实心三角 = changed / established，空心 = suggestive，其余为状态符号；
  F、M 两轴换算到生理方向，C、R 为酶容量与活性）。L2 及以上加实线粗框，L1 加虚线框，排名前若干的反应标名次；
- 代谢物节点：填色为代谢组 log2FC（通过 FDR 加粗描边）；? 未检出、/ 没有代谢组、– 没有可映射的特征；
- 聚合节点（方框）：只有组学证据，P、R、M 三格，不画箭头；
- 侧栏为各区域的覆盖计数与通路融合结果（只作探索排序），图注写来源、样本数、诱饵比与布局校验值。
"""
from __future__ import annotations

import math
import textwrap
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Rectangle

from omicsweave.visualization.theme import themed

from ..mappings.network_annotation import network_reactions
from ..mappings.network_layout import Layout
from .dossier_plots import AXIS_GLYPHS, DIVERGING, HAIRLINE, INK, INK2, MUTED, NORM, SURFACE, VMAX, _save

UNIT = 0.12                      # 英寸 / 画布单位
SIDEBAR = 36.0                   # 侧栏宽度（画布单位）
FOOTER = 17.0                    # 图例与图注高度
NODE_R = 0.85
CELL = 0.95
NOTE_SIZE = 3.2                  # 次要底物 / 产物注记的字号
MIN_EDGE = 0.8                   # F 轴有指向的边着色时的最小 |log2FC|
UP, DOWN = "#e34948", "#2a78d6"
EDGE_GREY, OMICS_GREY, FAINT = "#8f8d86", "#b9b7b0", "#d9d8d2"
COMPARTMENT_FACE, COMPARTMENT_EDGE = "#f6f5f1", "#c9c7bf"
NODE_GLYPHS = {"not_detected": "?", "not_measured": "/", "no_knowledge": "–", "qc_blocked": "Q"}
STRIP_AXES = (("F", "F_state", "F_direction"), ("C", "C_state", "C_direction"), ("R", "R_state", "R_direction"),
              ("M", "M_state", "M_direction"))


def _short(flux_id: str) -> str:
    for suffix in ("_c", "_m"):
        if flux_id.endswith(suffix):
            return flux_id[: -len(suffix)]
    return flux_id


def _shift(point: tuple[float, float], normal: tuple[float, float], offset: float) -> tuple[float, float]:
    return point[0] + normal[0] * offset, point[1] + normal[1] * offset


def _normal(geometry: dict) -> tuple[float, float]:
    ends = [p for _, p in geometry["substrates"]] + [p for _, p in geometry["products"]]
    if len(ends) < 2:
        return 0.0, 0.0
    a = geometry["substrates"][0][1] if geometry["substrates"] else geometry["hub"]
    b = geometry["products"][0][1] if geometry["products"] else geometry["hub"]
    dx, dy = b[0] - a[0], b[1] - a[1]
    norm = math.hypot(dx, dy) or 1.0
    return -dy / norm, dx / norm


def _edge_color(row: pd.Series) -> str:
    if row["draw_style"] == "not_assessed":
        return FAINT
    if row["draw_style"] == "omics_only":
        return OMICS_GREY
    value = row["value"]
    if pd.isna(value):
        return EDGE_GREY
    # F 轴有指向的边至少取色阶上 ±MIN_EDGE 的颜色：效应很小的边否则接近中性色，看起来像未评估；数值以视图表为准
    return DIVERGING(NORM(math.copysign(max(abs(value), MIN_EDGE), value)))


def _segment(ax, start, end, color, width, style, head: bool, shrink_end: float) -> None:
    ax.annotate("", xy=end, xytext=start, zorder=2,
                arrowprops=dict(arrowstyle="-|>" if head else "-", color=color, lw=width, linestyle=style,
                                shrinkA=0, shrinkB=shrink_end, mutation_scale=5.5, joinstyle="miter"))


def _edge_geometry(layout: Layout, flux_id: str) -> dict:
    """一条边画出来的几何：中心（含平行偏移）、折点、各端点（含偏移），以及用于避让的线段。"""
    edge = layout.edges[flux_id]
    geometry = layout.geometry(network_reactions()[flux_id])
    normal = _normal(geometry)
    hub = _shift(geometry["hub"], normal, edge.offset)
    via = [_shift(p, normal, edge.offset) for p in edge.via]
    substrates = [_shift(p, normal, edge.offset) for _, p in geometry["substrates"]]
    products = [_shift(p, normal, edge.offset) for _, p in geometry["products"]]
    segments = []
    for point in substrates:
        path = [hub] + via[::-1] + [point]
        segments += list(zip(path[:-1], path[1:]))
    segments += [(hub, p) for p in products]
    ins, outs = geometry["secondary"]
    name = lambda species: species.rsplit("_", 1)[0]
    note = ""
    if ins or outs:
        note = (", ".join(name(s) for s in dict.fromkeys(ins)) or "") + " → " + (", ".join(name(s) for s in dict.fromkeys(outs)) or "")
    return {"hub": hub, "via": via, "substrates": substrates, "products": products, "segments": segments,
            "note": note.strip()}


def _note_rect(geometry: dict) -> tuple | None:
    if not geometry["note"]:
        return None
    hx, hy = geometry["hub"]
    w = _text_width(geometry["note"], NOTE_SIZE)
    return (hx - w / 2, hy + 0.45, hx + w / 2, hy + 1.15)


# ---------------------------------------------------------------------------
# 四格小条与反应名的避让放置
# ---------------------------------------------------------------------------

def _em(fontsize: float) -> float:
    """字号（磅）对应的画布单位。"""
    return fontsize / 72 / UNIT


def _text_width(text: str, fontsize: float) -> float:
    return sum(0.56 if ord(ch) < 0x2E80 else 1.0 for ch in str(text)) * _em(fontsize)


def _strip_text(row: pd.Series) -> str:
    label = _short(row["element_id"])
    if row["draw_style"] == "not_assessed":
        return label
    extra = ""
    if row["flux_claim"] in ("L1", "L2", "L3"):
        extra = f" {row['flux_claim']}"
    if pd.notna(row["rank"]) and row["key"] in ("key", "l1", "top"):
        extra += f" #{int(row['rank'])}"
    # 锁死的通量组只在 F 轴有指向时标出（提醒同组的红边是同一次检验）；其余反应的组见视图表的 flux_group 列
    if row["flux_group"] and row["F_state"] in ("changed", "suggestive"):
        extra += f" [{_short(str(row['flux_group']))}]"
    return label + extra


def _strip_size(row: pd.Series) -> tuple[float, float]:
    if row["draw_style"] == "not_assessed":
        return _text_width(_strip_text(row), 3.6) + 0.2, 0.8
    return 4 * CELL + 0.3 + _text_width(_strip_text(row), 3.8) + 0.2, CELL + 0.3


def _overlap(a: tuple, b: tuple) -> float:
    dx = min(a[2], b[2]) - max(a[0], b[0])
    dy = min(a[3], b[3]) - max(a[1], b[1])
    return dx * dy if dx > 0 and dy > 0 else 0.0


def _segment_inside(segment: tuple, rect: tuple, samples: int = 12) -> float:
    (x0, y0), (x1, y1) = segment
    length = math.hypot(x1 - x0, y1 - y0)
    inside = sum(1 for i in range(samples + 1)
                 if rect[0] <= x0 + (x1 - x0) * i / samples <= rect[2] and rect[1] <= y0 + (y1 - y0) * i / samples <= rect[3])
    return length * inside / (samples + 1)


def _obstacles(view: pd.DataFrame, layout: Layout) -> list[tuple[tuple, float]]:
    """（矩形, 权重）：代谢物节点与名称、边界与库存、聚合方框、区域标题。"""
    out = []
    labels = view[view["element_type"] == "metabolite"].set_index("element_id")["label"].to_dict()
    for species, node in layout.nodes.items():
        for x, y in node.at:
            if node.kind == "metabolite":
                out.append(((x - NODE_R - 0.1, y - NODE_R - 0.1, x + NODE_R + 0.1, y + NODE_R + 0.1), 3.0))
                w = _text_width(labels.get(species, species), 3.8)
                out.append(((x - w / 2, y + NODE_R + 0.2, x + w / 2, y + NODE_R + 0.95), 2.0))
            elif node.kind == "boundary":
                w = _text_width(node.label or species, 3.4) + 1.2
                out.append(((x - 0.7, y - 0.7, x + 0.7, y + 0.7), 2.0))
                out.append(((x - w / 2, y + 1.2, x + w / 2, y + 1.85), 1.0))
            else:
                out.append(((x - 1.8, y - 0.8, x + 1.8, y + 0.8), 2.0))
    for agg in layout.aggregates.values():
        x, y = agg.xy
        out.append(((x - 4.8, y - 1.7, x + 4.8, y + 1.7), 2.0))
    for region in layout.regions.values():
        x, y = region.label_xy
        out.append(((x, y - 0.5, x + _text_width(region.label, 5.0), y + 0.5), 1.5))
    for comp in layout.compartments.values():
        x1, y0 = comp["box"][2], comp["box"][1]
        out.append(((x1 - 1 - _text_width(comp.get("label", ""), 4.6), y0 + 0.7, x1 - 1, y0 + 1.7), 1.0))
    return out


def _place_strips(view: pd.DataFrame, layout: Layout, geometries: dict[str, dict]) -> dict[str, tuple[float, float]]:
    """每条边的小条左端中点。布局写了 ``strip`` 时照用；否则在中心四周的候选位置中取与节点、名称、方框、已放好的小条
    及边的重叠最少的一个。重要的反应（L2 及以上、L1、排名靠前）先放。"""
    obstacles = _obstacles(view, layout) + [(r, 2.0) for g in geometries.values() if (r := _note_rect(g))]
    segments = [s for g in geometries.values() for s in g["segments"]]
    placed: list[tuple] = []
    out: dict[str, tuple[float, float]] = {}
    priority = {"key": 0, "l1": 1, "top": 2}
    rows = view[view["element_type"] == "reaction"].copy()
    rows["order"] = rows["key"].map(priority).fillna(3) + rows["draw_style"].eq("not_assessed") * 1
    for _, row in rows.sort_values("order", kind="stable").iterrows():
        flux_id = row["element_id"]
        hx, hy = geometries[flux_id]["hub"]
        w, h = _strip_size(row)
        explicit = layout.edges[flux_id].strip
        if explicit is not None:
            candidates = [explicit]
        else:
            gap = 0.9
            candidates = [(gap, 0), (gap, -1.5), (gap, 1.5), (-gap - w, 0), (-gap - w, -1.5), (-gap - w, 1.5),
                          (-w / 2, -1.7), (-w / 2, 1.7), (gap, -2.8), (gap, 2.8), (-gap - w, -2.8), (-gap - w, 2.8)]
        best, best_score = None, math.inf
        for index, (dx, dy) in enumerate(candidates):
            rect = (hx + dx, hy + dy - h / 2, hx + dx + w, hy + dy + h / 2)
            score = index * 0.02
            score += sum(_overlap(rect, r) * weight for r, weight in obstacles)
            score += sum(_overlap(rect, r) * 4.0 for r in placed)
            score += sum(_segment_inside(s, rect) for s in segments) * 0.35
            if rect[0] < -1 or rect[2] > layout.width + 1 or rect[1] < 0 or rect[3] > layout.height + 1:
                score += 50
            if score < best_score:
                best, best_score = (dx, dy, rect), score
        dx, dy, rect = best
        placed.append(rect)
        out[flux_id] = (hx + dx, hy + dy)
    return out


def _draw_edge(ax, layout: Layout, row: pd.Series, geometry: dict) -> tuple[float, float]:
    hub = geometry["hub"]
    style = {"flux": "-", "omics_only": (0, (3, 1.6)), "not_assessed": (0, (1, 1.6))}[row["draw_style"]]
    strong = row["key"] == "key"
    width = 0.5 if row["draw_style"] == "not_assessed" else (1.6 if strong else 1.0)
    color = _edge_color(row)
    physio = int(row["physio_sign"]) if pd.notna(row["physio_sign"]) else 0
    both = bool(row["sign_change"]) or (physio == 0 and bool(row["reversible"]))
    node_pt = NODE_R * UNIT * 72 + 1.0
    head_substrate, head_product = both or physio < 0, both or physio >= 0
    for point in geometry["substrates"]:
        # 底物一侧：中心 → 折点 → 底物；箭头在底物端时只给最后一段加箭头
        path = [hub] + geometry["via"][::-1] + [point]
        for i, (a, b) in enumerate(zip(path[:-1], path[1:])):
            last = i == len(path) - 2
            _segment(ax, a, b, color, width, style, head_substrate and last, node_pt if last and head_substrate else 0)
    for point in geometry["products"]:
        _segment(ax, hub, point, color, width, style, head_product, node_pt if head_product else 0)
    if bool(row["sign_change"]):
        ax.text(hub[0], hub[1] - 0.9, "↔", fontsize=4.5, ha="center", va="center", color=INK2, zorder=4)
    if geometry["note"]:
        ax.text(hub[0], hub[1] + 0.45, geometry["note"], fontsize=NOTE_SIZE, ha="center", va="top", color=MUTED,
                style="italic", zorder=4)
    return hub


def _draw_leader(ax, hub: tuple[float, float], position: tuple[float, float], row: pd.Series) -> None:
    """小条不紧挨反应中心时，画一条细引线连到小条较近的一端。"""
    w, _ = _strip_size(row)
    x0, y0 = position
    end = min(((x0, y0), (x0 + w, y0)), key=lambda p: math.dist(p, hub))
    if math.dist(end, hub) > 1.2:
        ax.add_line(Line2D([hub[0], end[0]], [hub[1], end[1]], color=MUTED, lw=0.3, linestyle=(0, (1, 1)), zorder=4))
    ax.add_patch(plt.Circle(hub, 0.18, facecolor=MUTED, edgecolor="none", zorder=4))


def _draw_strip(ax, position: tuple[float, float], row: pd.Series) -> None:
    x0, y0 = position
    if row["draw_style"] == "not_assessed":
        ax.text(x0, y0, _strip_text(row), fontsize=3.6, ha="left", va="center", color=MUTED, zorder=5)
        return
    for i, (_, state_col, dir_col) in enumerate(STRIP_AXES):
        x = x0 + i * CELL
        ax.add_patch(Rectangle((x, y0 - CELL / 2), CELL, CELL, facecolor="white", edgecolor=HAIRLINE, lw=0.3, zorder=5))
        state = str(row[state_col]) if pd.notna(row[state_col]) else "not_measured"
        direction = int(row[dir_col]) if pd.notna(row[dir_col]) else 0
        cx, cy = x + CELL / 2, y0
        if state in ("changed", "suggestive") and direction:
            color = UP if direction > 0 else DOWN
            ax.scatter([cx], [cy], s=5.5, marker="^" if direction > 0 else "v", linewidths=0.5, zorder=6,
                       facecolors=color if state == "changed" else "none", edgecolors=color)
        else:
            ax.text(cx, cy, AXIS_GLYPHS.get(state, "·"), fontsize=3.4, ha="center", va="center", zorder=6,
                    color=INK2 if state in ("unchanged", "mixed") else MUTED)
    frame = {"key": ("-", 1.1), "l1": ((0, (2, 1)), 0.7)}.get(row["key"])
    if frame:
        ax.add_patch(Rectangle((x0 - 0.15, y0 - CELL / 2 - 0.15), 4 * CELL + 0.3, CELL + 0.3, facecolor="none",
                               edgecolor=INK, lw=frame[1], linestyle=frame[0], zorder=6))
    ax.text(x0 + 4 * CELL + 0.3, y0, _strip_text(row), fontsize=3.8, ha="left", va="center", zorder=6,
            color=INK if row["key"] in ("key", "l1") else INK2, fontweight="bold" if row["key"] == "key" else "normal",
            bbox=dict(boxstyle="square,pad=0.05", facecolor=SURFACE, edgecolor="none", alpha=0.75))


def _draw_node(ax, layout: Layout, species: str, row: pd.Series | None) -> None:
    node = layout.nodes[species]
    for x, y in node.at:
        if node.kind == "boundary":
            ax.add_patch(Rectangle((x - 0.55, y - 0.55), 1.1, 1.1, facecolor="white", edgecolor=MUTED, lw=0.5,
                                   zorder=3))
            name = node.label or species.rsplit("_", 1)[0]
            ax.text(x, y + 1.25, f"{name}（血）" if species.endswith("_e") else name, fontsize=3.4, ha="center",
                    va="top", color=MUTED, style="italic", zorder=4)
            continue
        if node.kind == "store":
            ax.add_patch(FancyBboxPatch((x - 1.6, y - 0.6), 3.2, 1.2, boxstyle="round,pad=0.1,rounding_size=0.4",
                                        facecolor="#eceae4", edgecolor=MUTED, lw=0.4, zorder=3))
            ax.text(x, y, node.label or species, fontsize=3.4, ha="center", va="center", color=INK2, zorder=4)
            continue
        state = str(row["display_state"]) if row is not None and pd.notna(row["display_state"]) else "not_measured"
        value = row["value"] if row is not None else np.nan
        face = DIVERGING(NORM(value)) if pd.notna(value) else SURFACE
        ring = row["ring"] if row is not None and isinstance(row["ring"], str) else ""
        edge_color, lw = (INK, 1.1) if ring in ("selected", "changed") else (MUTED, 0.45)
        ax.add_patch(plt.Circle((x, y), NODE_R, facecolor=face, edgecolor=edge_color, lw=lw, zorder=4))
        if pd.isna(value) and state in NODE_GLYPHS:
            ax.text(x, y, NODE_GLYPHS[state], fontsize=3.8, ha="center", va="center", color=MUTED, zorder=5)
        label = row["label"] if row is not None else species
        ax.text(x, y + NODE_R + 0.25, str(label), fontsize=3.8, ha="center", va="top", color=INK, zorder=5)


def _draw_aggregate(ax, layout: Layout, row: pd.Series) -> None:
    agg = layout.aggregates[row["element_id"]]
    x, y = agg.xy
    w, h = 9.2, 3.0
    faint = row["draw_style"] == "not_assessed"
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.05,rounding_size=0.5",
                                facecolor="white", edgecolor=FAINT if faint else OMICS_GREY, lw=0.6,
                                linestyle=(0, (3, 1.6)), zorder=3))
    frame = {"key": ("-", 1.1), "l1": ((0, (2, 1)), 0.8)}.get(row["key"])
    if frame:
        ax.add_patch(FancyBboxPatch((x - w / 2 - 0.2, y - h / 2 - 0.2), w + 0.4, h + 0.4,
                                    boxstyle="round,pad=0.05,rounding_size=0.6", facecolor="none", edgecolor=INK,
                                    lw=frame[1], linestyle=frame[0], zorder=3))
    title = str(row["label"]) + (f"  {row['flux_claim']}" if row["flux_claim"] in ("L1", "L2", "L3") else "")
    ax.text(x, y - 0.55, title, fontsize=3.6, ha="center", va="center", color=MUTED if faint else INK, zorder=4)
    for i, (code, value_col, state_col) in enumerate((("P", "P_value", "P_state"), ("R", "R_value", "R_layer_state"),
                                                     ("M", "Mx_value", "Mx_state"))):
        cx = x - 2.4 + i * 2.4
        value = row[value_col]
        state = str(row[state_col]) if pd.notna(row[state_col]) else "not_measured"
        face = DIVERGING(NORM(value)) if pd.notna(value) else "white"
        ax.add_patch(Rectangle((cx - 0.45, y + 0.15), 0.9, 0.9, facecolor=face, edgecolor=HAIRLINE, lw=0.3, zorder=4))
        glyph = "" if pd.notna(value) else NODE_GLYPHS.get(state, AXIS_GLYPHS.get(state, ""))
        ax.text(cx - 0.6, y + 0.6, code, fontsize=3.0, ha="right", va="center", color=MUTED, zorder=5)
        if glyph:
            ax.text(cx, y + 0.6, glyph, fontsize=3.2, ha="center", va="center", color=MUTED, zorder=5)


def _sidebar(ax, regions: pd.DataFrame, x0: float, y0: float) -> None:
    ax.text(x0, y0, "各区域读数（计数描述完整度，不是评分）", fontsize=5.2, fontweight="bold", color=INK, va="top")
    y = y0 + 2.2
    for row in regions.itertuples():
        if row.region == "other":
            continue
        ax.text(x0, y, f"{row.label}  {row.kegg}", fontsize=4.4, fontweight="bold", color=INK, va="top")
        y += 1.5
        ax.text(x0, y, f"画出 {row.n_drawn} · 有通量 {row.n_flux} · 可辨识 {row.n_identifiable} · 只有组学 "
                       f"{row.n_omics_only} · L2+ {row.n_L2} · L1 {row.n_L1}", fontsize=3.8, color=INK2, va="top")
        y += 1.4
        fusion = row.fusion or "通路融合：无匹配的多层通路行"
        for line in textwrap.wrap(("通路融合（探索）：" if row.fusion else "") + fusion, 64)[:3]:
            ax.text(x0, y, line, fontsize=3.5, color=MUTED, va="top")
            y += 1.2
        y += 0.9


def _legend(ax, y0: float, width: float, caption: list[str]) -> None:
    # 色阶
    xs = np.linspace(0, 14, 57)
    for i, x in enumerate(xs[:-1]):
        value = -VMAX + 2 * VMAX * i / (len(xs) - 2)
        ax.add_patch(Rectangle((1 + x, y0), xs[1] - xs[0], 0.9, facecolor=DIVERGING(NORM(value)), edgecolor="none"))
    for v in (-2, 0, 2):
        ax.text(1 + (v + VMAX) / (2 * VMAX) * 14, y0 + 1.5, f"{v:+d}" if v else "0", fontsize=3.8, ha="center",
                va="top", color=INK2)
    ax.text(1, y0 - 0.5, "log2FC（节点：代谢组；边：F 轴，按生理方向，|log2FC| < 0.8 的边按 0.8 着色）", fontsize=3.8,
            color=INK2, va="bottom")
    items = [((0, ()), "-", EDGE_GREY, "有拟合通量"), (None, (0, (3, 1.6)), OMICS_GREY, "只有组学（该组织模型没有这一步）"),
             (None, (0, (1, 1.6)), FAINT, "未评估（没有通量，也没有成员知识）")]
    for i, (_, style, color, text) in enumerate(items):
        x = 18 + i * 17
        ax.add_line(Line2D([x, x + 3], [y0 + 0.4, y0 + 0.4], color=color, lw=1.0, linestyle=style))
        ax.text(x + 3.5, y0 + 0.4, text, fontsize=3.8, va="center", color=INK2)
    glyphs = "四格小条 F C R M：▲▼ 实心 changed（F：established）· 空心 suggestive · = 不变 · ≠ 混合 · · 未定 · " \
             "x 不可辨识 · ! 样本不足 · / 未测 · ? 未检出 · – 无知识 ｜ 实线框 L2 及以上，虚线框 L1，#n 排名，[…] 锁死的通量组"
    ax.text(18, y0 + 2.2, glyphs, fontsize=3.6, va="center", color=INK2)
    ax.text(18, y0 + 3.6, "箭头为生理净方向；↔ 两组间净通量变号。F、M 两轴的方向已换算到生理方向；C、R 为酶容量与活性。"
                          "红边多不等于通路上调：锁死的反应是同一次检验。", fontsize=3.6, va="center", color=INK2)
    y = y0 + 5.4
    for line in caption:
        for wrapped in textwrap.wrap(line, 230):
            ax.text(1, y, wrapped, fontsize=3.7, va="top", color=INK2)
            y += 1.25


def plot_pathway_map(view: pd.DataFrame, layout: Layout, regions: pd.DataFrame, title: str, caption: list[str],
                     stem: Path, formats: Iterable[str] = ("pdf", "png")) -> list[Path]:
    width, height = layout.width + SIDEBAR, layout.height + FOOTER
    with themed(dpi=150, save_dpi=300):
        fig = plt.figure(figsize=(width * UNIT, (height + 5) * UNIT))
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_xlim(-1, width)
        ax.set_ylim(height, -5)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.text(0, -3.2, title, fontsize=7.5, fontweight="bold", color=INK, va="center")
        for comp in layout.compartments.values():
            x0, y0, x1, y1 = comp["box"]
            ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0, boxstyle="round,pad=0,rounding_size=2.5",
                                        facecolor=COMPARTMENT_FACE, edgecolor=COMPARTMENT_EDGE, lw=1.0, zorder=0))
            ax.text(x1 - 1, y0 + 1.2, comp.get("label", ""), fontsize=4.6, ha="right", va="center", color=MUTED)
        for region in layout.regions.values():
            ax.text(*region.label_xy, region.label, fontsize=5.0, fontweight="bold", color=INK2, va="center", zorder=1)
        reactions = view[view["element_type"] == "reaction"]
        geometries = {flux_id: _edge_geometry(layout, flux_id) for flux_id in reactions["element_id"]}
        positions = _place_strips(view, layout, geometries)
        for _, row in reactions.iterrows():
            _draw_edge(ax, layout, row, geometries[row["element_id"]])
        for _, row in reactions.iterrows():
            _draw_leader(ax, geometries[row["element_id"]]["hub"], positions[row["element_id"]], row)
            _draw_strip(ax, positions[row["element_id"]], row)
        nodes = view[view["element_type"] == "metabolite"].set_index("element_id")
        for species in layout.nodes:
            _draw_node(ax, layout, species, nodes.loc[species] if species in nodes.index else None)
        for _, row in view[view["element_type"] == "aggregate"].iterrows():
            _draw_aggregate(ax, layout, row)
        _sidebar(ax, regions, layout.width + 2, 0)
        _legend(ax, layout.height + 2.5, width, caption)
        return _save(fig, stem, formats)
