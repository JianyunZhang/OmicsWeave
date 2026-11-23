#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
反应档案的图（G14）：

- 图 A :func:`plot_matrix`：通路证据矩阵，只画 :func:`omicsweave.integration.dossier.view.matrix_view` 的表；
- 图 B :func:`plot_dossier_page`：一个（反应, 组织, 比较）的档案页，每个成员一个面板，画出每只动物的点；
- 图 C :func:`plot_tracer_coverage`：面板中各反应 × 示踪剂的可辨识性；
- 图 D :func:`plot_node_split`：分支节点的分流比（各层的逐只动物点图并排）。

编码（设计文档 4.5）：颜色为带符号的效应量（log2FC，蓝降红升，零点为中性灰，色阶固定并截断于 ±2）；点的大小为
方向的确定程度；通过 FDR 的格子加描边（同时达到入选规则的描边更粗）；不可辨识、未测、未检出等状态各有固定的
符号，不画成空白或零效应的颜色；图例出现在每张图上。档案页横轴为数据中的全部组别，当前比较的分子组与分母组
着色，其余组为灰色。
"""
from __future__ import annotations

import math
import re
import textwrap
import unicodedata
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

from omicsweave.mfa.tracers import format_tracer_display, format_tracer_short
from omicsweave.visualization.figures import save_stem
from omicsweave.visualization.palettes import (
    DEFAULT_NORM,
    DEFAULT_VMAX,
    IDENTIFIABILITY_COLORS,
    METRO_ATLAS,
    SemanticPalette,
    get_palette,
)
from omicsweave.visualization.primitives import (
    draw_layer_icon,
    draw_status_cell,
    draw_tier_badge,
    estimate_text_width,
    tier_badge_width,
)
from omicsweave.visualization.theme import DEFAULT_PRESET, add_economist_tag, style_economist_axes, themed

from ..dossier.panels import Panel

INK = METRO_ATLAS.ink
INK2 = METRO_ATLAS.ink2
MUTED = METRO_ATLAS.muted
HAIRLINE = METRO_ATLAS.hairline
SURFACE = METRO_ATLAS.surface
NUMERATOR = METRO_ATLAS.numerator
DENOMINATOR = METRO_ATLAS.denominator
DIVERGING = METRO_ATLAS.diverging_cmap()
VMAX = DEFAULT_VMAX
NORM = DEFAULT_NORM
_SPECIAL_DATA_STATES = frozenset({
    "not_identifiable", "insufficient_n", "not_measured", "not_infused",
    "not_detected", "no_knowledge", "qc_blocked",
})
# 兼容保留字典导出名，视觉渲染统一使用 primitives 纯矢量状态图标
GLYPHS = {"equivalent": "—", "not_identifiable": "⊘", "insufficient_n": "⚠", "not_measured": "/",
          "not_infused": "/", "not_detected": "◌", "no_knowledge": "/", "qc_blocked": "⚠"}
BLOCK_TITLES = {"flux": "FLUX by tracer", "enzyme": "ENZYME", "regulator": "REGULATOR", "pool": "POOL",
                "assessment": "ASSESSMENT"}
AXIS_GLYPHS = {"unchanged": "—", "mixed": "≠", "undetermined": "╌", "not_identifiable": "⊘", "insufficient_n": "⚠",
               "not_detected": "◌", "not_measured": "/", "qc_blocked": "⚠", "no_knowledge": "/"}
IDENTIFIABILITY = dict(IDENTIFIABILITY_COLORS)
_HEADER_LAYER_ICONS = {"P": "protein", "R": "mrna", "M": "metabolite"}
_ASSESSMENT_COL_ICONS = {"F": "flux", "C": "protein", "R": "regulation", "M": "metabolite"}
_CONCORDANT_RE = re.compile(r"\((\d+)\)")
_ABS_G_RE = re.compile(r"\|g\|\s*=\s*([0-9]+(?:\.[0-9]+)?)")
_ABS_G_LT_RE = re.compile(r"\|g\|\s*=\s*([0-9]+(?:\.[0-9]+)?)\s*<")


def _undetermined_tag(cell: pd.Series) -> str:
    """F 轴 undetermined 时的精简门控原因标签（纯 ASCII，避免字库缺字），排在虚线下方。"""
    if str(cell.get("column") or "") != "F" or str(cell.get("display_state") or "") != "undetermined":
        return ""
    flags = {f.strip() for f in str(cell.get("flags") or "").split(";") if f.strip() and f.strip() != "nan"}
    text = str(cell.get("text") or "")
    if "tracer_expectation_mismatch" in flags:
        return "mismatch"
    if "single_subject_driven" in flags or "单只动物驱动" in text:
        return "1-subj"
    if "stability_unknown" in flags or "稳定性未知" in text:
        return "LOO?"
    m_lt = _ABS_G_LT_RE.search(text)
    if m_lt:
        return f"g={float(m_lt.group(1)):.2f}"
    m_g = _ABS_G_RE.search(text)
    if m_g:
        g_val = float(m_g.group(1))
        return f"g={g_val:.2f}" if g_val < 1.0 else "LOO"
    if "LOO" in text or "方向稳定性" in text or "不稳定" in text:
        return "LOO"
    if "|g| 缺失" in text:
        return "g=NA"
    return ""


def _size(certainty: float) -> float:
    if certainty is None or not np.isfinite(certainty):
        return 70.0
    return 25.0 + 175.0 * float(np.clip((certainty - 0.5) / 0.5, 0, 1))


def _save(fig: plt.Figure, stem: Path, formats: Iterable[str], *, surface: str = SURFACE) -> list[Path]:
    return save_stem(fig, stem, formats=tuple(formats), dpi=300, facecolor=surface, close=True)


def _draw_cell(
    ax: plt.Axes,
    x: float,
    y: float,
    cell: pd.Series,
    *,
    y_scale: float = 1.25,
    pal: SemanticPalette = METRO_ATLAS,
    diverging=None,
) -> None:
    """一个格子：点在上，log2FC 写在点的下方；特殊/缺失状态使用统一纯矢量图标。"""
    cmap = diverging if diverging is not None else pal.diverging_cmap()
    state = str(cell["display_state"] or "")
    if state in _SPECIAL_DATA_STATES:
        s = 0.66
        draw_status_cell(
            ax, x - s / 2, y - s * y_scale / 2, s, state, 0,
            palette=pal, font_scale=1.18, zorder=3.0, draw_box=False, y_up=True, y_scale=y_scale,
        )
        return
    value, cy = cell["value"], y + 0.13
    if state == "equivalent":
        ax.scatter([x], [cy], s=_size(cell["certainty"]), color=pal.card_bg, edgecolors=pal.hairline, linewidths=0.6,
                   zorder=2)
        ax.plot([x - 0.14, x + 0.14], [cy, cy], color=pal.ink2, lw=1.0, solid_capstyle="round", zorder=3)
    else:
        face = cmap(NORM(value)) if pd.notna(value) else "none"
        ring = cell["ring"]
        edge, width = (pal.ink, 1.6) if ring == "selected" else (pal.ink, 0.8) if ring == "changed" else (pal.muted, 0.4)
        ax.scatter([x], [cy], s=_size(cell["certainty"]), color=face, edgecolors=edge, linewidths=width, zorder=2)
    if pd.notna(value):
        ax.text(x, y - 0.3, f"{value:+.2f}", ha="center", va="center", fontsize=6.0, color=pal.ink2, zorder=3)
    if state == "summary" and cell["n_quantified"]:
        ax.text(x + 0.3, y + 0.38, f"{int(cell['n_changed'])}/{int(cell['n_quantified'])}", ha="left", va="bottom",
                fontsize=6.0, color=pal.ink2, zorder=3)


def _draw_assessment(
    ax: plt.Axes,
    x: float,
    y: float,
    cell: pd.Series,
    *,
    y_scale: float = 1.25,
    pal: SemanticPalette = METRO_ATLAS,
    preset: str = DEFAULT_PRESET,
) -> None:
    """评估块的一格：四条轴统一纯矢量图标，结论等级使用自解释三点徽章（Confirmed / Flux / Enz-pred）。"""
    state = str(cell["display_state"] or "")
    direction = int(cell["value"]) if pd.notna(cell["value"]) else 0
    if cell["column"] == "level":
        if state in ("L1", "L2", "L3"):
            b_sc, b_fs = 0.42, 6.0
            bw = tier_badge_width(state, scale=b_sc, unit_in=0.50, fontsize=b_fs)
            m = _CONCORDANT_RE.search(str(cell.get("text") or ""))
            suffix = f" ({m.group(1)})" if m else ""
            suf_w = estimate_text_width(suffix, 6.0, 0.50) if suffix else 0.0
            bx = x - (bw + suf_w) / 2.0
            draw_tier_badge(
                ax, bx, y, state, zorder=3.2, scale=b_sc, preset=preset,
                unit_in=0.50, fontsize=b_fs, y_scale=y_scale,
            )
            if suffix:
                ax.text(bx + bw + 0.06, y, suffix, ha="left", va="center", fontsize=6.0, color=pal.muted, zorder=3.2)
        elif state == "LX":
            ax.text(x, y, "Conflict", ha="center", va="center", fontsize=6.2, fontweight="bold", color=pal.ink2, zorder=3.1)
        else:
            ax.plot([x - 0.16, x + 0.16], [y, y], color=pal.muted, lw=0.95, solid_capstyle="round", zorder=3.1)
        return
    tag = _undetermined_tag(cell)
    s = 0.62 if tag else 0.66
    cy = (y + 0.10) if tag else y
    draw_status_cell(
        ax, x - s / 2, cy - s * y_scale / 2, s, state, direction,
        palette=pal, font_scale=1.22, zorder=3.0, draw_box=False, y_up=True, y_scale=y_scale,
    )
    if tag:
        ax.text(x, y - 0.30, tag, ha="center", va="center", fontsize=6.0, color=pal.muted, zorder=3.1)


def _short_label(block: str, column: str, tracer: str) -> str:
    """列的短名：通量列为规范同位素示踪剂短标（如 U-13C-Glc / Pooled），证据等级列为 Evidence，其余为层/轴代码。"""
    if block == "flux":
        return format_tracer_short(str(tracer))
    if block == "assessment" and column == "level":
        return "Evidence"
    return str(column).split(":")[0]


def plot_matrix(view: pd.DataFrame, panel: Panel, tissue: str, contrast: str, stem: Path,
                formats: Iterable[str] = ("pdf", "png"), runs: dict[str, str] | None = None,
                caption: str = "", *, style: str = DEFAULT_PRESET) -> list[Path]:
    """图 A：只画 ``view``（矩阵视图）中的格子。"""
    pal = get_palette(style)
    diverging = pal.diverging_cmap()
    rows = view[["row", "element_type", "element_id", "label", "segment"]].drop_duplicates().sort_values("row")
    row_y = {r: -i for i, r in enumerate(rows["row"])}
    sources = list(dict.fromkeys(view["source"]))
    # 每块内的列按来源排在一起（来源写在上方的来源行），来源内保持视图中的顺序
    columns = view[["block", "column", "source", "tracer"]].drop_duplicates(["block", "column"])
    columns = columns.assign(_s=columns["source"].map(sources.index), _o=range(len(columns)))
    x, positions, labels, block_spans, source_spans = 0.0, {}, {}, {}, []
    for block in ("flux", "enzyme", "regulator", "pool", "assessment"):
        part = columns[columns["block"] == block].sort_values(["_s", "_o"])
        if part.empty:
            continue
        start = x
        for source, group in part.groupby("_s", sort=True):
            first = x
            for r in group.itertuples():
                step = 2.8 if (block, r.column) == ("assessment", "level") else 1.0
                positions[(block, r.column)] = x + (step - 1.0) / 2
                labels[(block, r.column)] = _short_label(block, r.column, r.tracer)
                x += step
            if block != "assessment":                # 评估块综合全部来源，不写来源行
                source_spans.append((first, x - 1.0, sources[source]))
            x += 0.25                                # 来源之间留一点空
        x -= 0.25
        block_spans[block] = (start, x - 1.0)
        x += 0.8
    last = x - 1.8                                   # 最后一列的位置；分段名写在网格右侧之外
    right = last + 3.2
    x_unit_in = 0.50
    row_h = 0.42
    top_in = 2.08 if style == "economist" else 1.95
    bot_in = 2.25
    width = max(8.2, 2.6 + x_unit_in * (right + 0.7))
    x_span = right + 0.7
    eff_x_in = (width - 2.6) / x_span
    y_scale = eff_x_in / row_h
    height = top_in + row_h * len(rows) + bot_in
    with themed(preset=style):
        fig = plt.figure(figsize=(width, height), facecolor=pal.surface)
        if style == "economist":
            add_economist_tag(fig, x=0.015, y=1 - 0.08 / height, width=min(0.065, 0.52 / width),
                              height=0.09 / height, rule=True, rule_right=0.985)
        ax = fig.add_axes((2.4 / width, bot_in / height, (width - 2.6) / width, (height - top_in - bot_in) / height))
        ax.set_facecolor(pal.surface)
        ax.set_xlim(-0.7, right)
        ax.set_ylim(-len(rows) + 0.4, 0.8)
        ax.axis("off")
        previous_segment = None
        stripe_bg = "#E2E1D5" if style == "economist" else "#f4f3ef"
        for _, row in rows.iterrows():
            y = row_y[row["row"]]
            reaction = row["element_type"] == "reaction"
            ax.text(-0.9, y, f"{row['label']}", ha="right", va="center", fontsize=8 if reaction else 7,
                    color=pal.ink if reaction else pal.ink2, fontweight="bold" if reaction else "normal",
                    transform=ax.transData, clip_on=False)
            if reaction:
                ax.add_patch(plt.Rectangle((-0.7, y - 0.5), last + 1.2, 1.0, color=stripe_bg, zorder=0, lw=0))
            if row["segment"] != previous_segment:
                ax.plot([-0.7, last + 0.5], [y + 0.5, y + 0.5],
                        color=pal.ink if style == "economist" else pal.hairline,
                        lw=0.75 if style == "economist" else 0.6, zorder=1)
                ax.text(last + 0.8, y + 0.45, row["segment"], ha="left", va="top", fontsize=6.5, color=pal.muted)
                previous_segment = row["segment"]
        for key, cx in positions.items():
            blk, _col = key
            name = labels[key]
            icon_kind = (
                _ASSESSMENT_COL_ICONS.get(name) if blk == "assessment"
                else (_HEADER_LAYER_ICONS.get(name) if blk != "flux" else None)
            )
            if icon_kind:
                draw_layer_icon(
                    ax, icon_kind, cx, 0.96, size=0.38, color=pal.ink2, lw=0.72,
                    zorder=3.0, y_scale=y_scale, clip_on=False,
                )
                ax.text(cx, 1.38, name, ha="center", va="bottom", fontsize=6.2, fontweight="bold", color=pal.ink2, clip_on=False)
            elif name == "Evidence":
                ax.text(cx, 0.78, name, ha="center", va="bottom", fontsize=6.5, fontweight="bold", color=pal.ink, clip_on=False)
            else:
                ax.text(cx, 0.68, name, ha="center", va="bottom", rotation=90 if len(name) > 2 else 0, fontsize=6.5,
                        color=pal.ink2, clip_on=False)
        for start, end, source in source_spans:      # 来源行：来源名居中，下划细线
            ax.plot([start - 0.35, end + 0.35], [2.15, 2.15], color=pal.muted, lw=0.6, clip_on=False)
            ax.text((start + end) / 2, 2.22, source, ha="center", va="bottom", fontsize=7, color=pal.ink,
                    clip_on=False)
        for block, (start, end) in block_spans.items():
            ax.text((start + end) / 2, 2.85, BLOCK_TITLES[block], ha="center", va="bottom", fontsize=7.5, color=pal.ink,
                    fontweight="bold", clip_on=False)
        for _, cell in view.iterrows():
            cx, cy = positions[(cell["block"], cell["column"])], row_y[cell["row"]]
            if cell["block"] == "flux" and bool(cell["expected_tracer"]):
                ax.add_patch(plt.Rectangle((cx - 0.45, cy - 0.45), 0.9, 0.9, fill=False, lw=0.6, ls=(0, (2, 1.5)),
                                           ec=pal.muted, zorder=1))
            if cell["block"] == "assessment":
                _draw_assessment(ax, cx, cy, cell, y_scale=y_scale, pal=pal, preset=style)
            else:
                _draw_cell(ax, cx, cy, cell, y_scale=y_scale, pal=pal, diverging=diverging)
        title_y = (1 - 0.32 / height) if style == "economist" else (1 - 0.28 / height)
        fig.text(0.015, title_y, f"{panel.title}  ·  {tissue}  ·  {contrast}", fontsize=10.5, color=pal.ink,
                 fontweight="bold", va="top")
        _matrix_legend(fig, width, height, runs or {}, caption=caption, pal=pal, preset=style, diverging=diverging)
        return _save(fig, stem, formats, surface=pal.surface)


def _matrix_legend(
    fig: plt.Figure,
    width: float,
    height: float,
    runs: dict[str, str],
    caption: str = "",
    *,
    pal: SemanticPalette = METRO_ATLAS,
    preset: str = DEFAULT_PRESET,
    diverging=None,
) -> None:
    """图 A 统一纯矢量视觉图例（与图 E / 子图 100% 同一套色块、状态图标与三点徽章）。"""
    cmap = diverging if diverging is not None else pal.diverging_cmap()
    # 顶部横排：Colorbar + 点大小/描边图例
    cax = fig.add_axes((0.35 / width, 1.62 / height, 1.55 / width, 0.13 / height))
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=NORM, cmap=cmap), cax=cax, orientation="horizontal")
    bar.set_ticks([-2, -1, 0, 1, 2])
    bar.ax.tick_params(labelsize=6.0, colors=pal.ink2, length=2)
    bar.outline.set_edgecolor(pal.hairline)
    cax.set_title("log2FC (regulators: × activity sign)", fontsize=6.4, color=pal.ink2, pad=3)

    handles = [plt.scatter([], [], s=_size(c), color="#d9d8d3", edgecolors=pal.muted, linewidths=0.4)
               for c in (0.6, 0.8, 1.0)]
    handles += [plt.scatter([], [], s=85, color="#d9d8d3", edgecolors=pal.ink, linewidths=w) for w in (0.8, 1.6)]
    handles += [plt.scatter([], [], s=85, color="none", edgecolors=pal.muted, linewidths=0.6)]
    labels = ["P(sign) 0.6", "0.8", "1.0", "FDR < 0.05", "FDR & effect", "unsigned"]
    fig.legend(handles, labels, loc="lower left", bbox_to_anchor=(2.10 / width, 1.48 / height), ncol=6, fontsize=6.0,
               frameon=False, handletextpad=0.25, columnspacing=0.9)

    # 下方等比英寸坐标图例画布（1 unit = 1 inch），直接渲染真实矢量图标与三点徽章
    lw_in = width - 0.60
    lh_in = 1.38
    lax = fig.add_axes((0.30 / width, 0.06 / height, lw_in / width, lh_in / height))
    lax.set_facecolor(pal.surface)
    lax.set_xlim(0.0, lw_in)
    lax.set_ylim(0.0, lh_in)
    lax.axis("off")

    sw = 0.165
    # Row 1 (y = 1.18): 方向与状态矢量图标（显式区分 Unchanged 等效未变 与 Undetermined 证据不足/不稳定）
    y1 = 1.18
    x = 0.02
    for st_pair, desc in (
        ((("changed", 1), ("changed", -1)), "Sig. (FDR<0.05)"),
        ((("suggestive", 1), ("suggestive", -1)), "Suggestive (p<0.05 / |g|≥1)"),
        ((("mixed", 0),), "Mixed"),
        ((("unchanged", 0),), "Unchanged (equiv.)"),
        ((("undetermined", 0),), "Undetermined"),
        ((("not_measured", 0),), "Not measured"),
        ((("not_identifiable", 0),), "Unident. flux"),
        ((("not_detected", 0),), "Not detected"),
        ((("qc_blocked", 0),), "QC / Low-n"),
    ):
        item_w = len(st_pair) * (sw + 0.035) + 0.04 + estimate_text_width(desc, 6.0, 1.0) + 0.18
        if x > 0.05 and x + item_w > lw_in:
            y1 -= 0.23
            x = 0.02
        for st, dr in st_pair:
            draw_status_cell(
                lax, x, y1 - sw / 2, sw, st, dr,
                palette=pal, font_scale=0.95, zorder=2.5, draw_box=True, y_up=True, y_scale=1.0,
            )
            x += sw + 0.035
        lax.text(x + 0.03, y1, desc, fontsize=6.0, color=pal.ink2, va="center", ha="left")
        x += 0.04 + estimate_text_width(desc, 6.0, 1.0) + 0.18

    # Row 2 (y = y1 - 0.26): 自解释三点证据徽章 + 4 层图标说明
    y2 = y1 - 0.26
    x = 0.02
    for lvl, desc in (
        ("L3", "Multi-tracer"),
        ("L2", "MFA supported"),
        ("L1", "Enzyme only"),
    ):
        bw = draw_tier_badge(
            lax, x, y2, lvl, zorder=2.5, scale=0.18, preset=preset,
            unit_in=1.0, fontsize=6.0, y_scale=1.0,
        )
        lax.text(x + bw + 0.07, y2, desc, fontsize=6.0, color=pal.ink2, va="center", ha="left")
        x += bw + 0.07 + estimate_text_width(desc, 6.0, 1.0) + 0.24
    # 附带 4 层微图标含义
    for ikind, ilbl in (("flux", "Flux"), ("protein", "Prot"), ("mrna", "mRNA"), ("metabolite", "Metab")):
        iw = 0.16 + estimate_text_width(ilbl, 6.0, 1.0) + 0.16
        if x + iw > lw_in:
            break
        draw_layer_icon(lax, ikind, x + 0.07, y2, size=0.15, color=pal.ink, lw=0.72, zorder=2.5)
        lax.text(x + 0.17, y2, ilbl, fontsize=6.0, color=pal.ink2, va="center", ha="left")
        x += iw

    # Row 3 & 4: 来源、说明（含 Flux 轴判定口径）与校准注记
    y3 = max(0.32, y2 - 0.25)
    note_parts = [
        "Flux: Sig.=≥2 tracers & FDR<0.05; Sugg.=|g|≥1 or (≥2 & p<0.05 & |g|≥0.8); —=equiv.; dashed=undet (tag: gate)",
        "k/n: changed/quantified",
        "dashed box: expected tracer",
        "number under dot: log2FC",
    ]
    if runs:
        note_parts.insert(0, "Sources: " + "; ".join(f"{k} = {v}" for k, v in runs.items()))
    lax.text(0.02, y3, "  ·  ".join(note_parts), fontsize=6.0, color=pal.muted, va="center", ha="left")
    if caption:
        lax.text(0.02, max(0.08, y3 - 0.20), caption, fontsize=6.0, color=pal.ink2, va="center", ha="left")


# ---------------------------------------------------------------------------
# 图 B
# ---------------------------------------------------------------------------

def _page_items(rows: pd.DataFrame) -> list[tuple[str, list[pd.Series]]]:
    """档案页的分区与每个分区的成员行。"""
    def pick(mask: pd.Series) -> list[pd.Series]:
        part = rows[mask & rows["cell_state"].ne("no_knowledge")].copy()
        part["_share"] = pd.to_numeric(part["abundance_share"], errors="coerce").fillna(-1)
        part = part.sort_values(["source", "_share", "member_id"], ascending=[True, False, True], kind="stable")
        return [r for _, r in part.iterrows()]

    mfa = rows["layer"].eq("mfa")
    # 有分示踪剂结果的来源画各示踪剂；只有合并估计的来源（如没有分示踪剂表）画合并的一格
    per_tracer = set(rows.loc[mfa & rows["tracer"].ne("pooled") & rows["cell_state"].ne("not_infused"), "source"])
    flux = mfa & (rows["tracer"].ne("pooled") | ~rows["source"].isin(per_tracer))
    enzyme = rows["role"].isin(["direct_enzyme", "carrier"])
    sections = [
        ("Flux by tracer (per animal)", pick(flux)),
        ("Enzyme · protein", pick(enzyme & rows["layer"].eq("proteomics"))),
        ("Enzyme · RNA", pick(enzyme & rows["layer"].eq("transcriptomics"))),
        ("Regulators and effectors", pick(rows["role"].eq("regulator") | rows["role"].eq("effector"))),
        ("Substrates and products", pick(rows["role"].isin(["substrate", "product"]) & rows["layer"].eq("metabolomics"))),
        ("Product / substrate ratio", pick(rows["layer"].eq("derived"))),
    ]
    return [(title, items) for title, items in sections if items]


def _member_samples(row: pd.Series, samples: pd.DataFrame) -> pd.DataFrame:
    if not isinstance(row.get("feature_id"), str) or not row["feature_id"]:
        return samples.iloc[0:0]
    part = samples[(samples["source"] == row["source"]) & (samples["layer"] == row["layer"])
                   & (samples["tissue"] == row["tissue"]) & (samples["feature_id"] == row["feature_id"])]
    if row["layer"] == "mfa" and row["tracer"] != "pooled":
        part = part[part["tracer"] == row["tracer"]]
    return part


# 图 B 面板标题每行最多的显示宽度（ASCII 记 1、全角字符记 2）：6 pt 字号下约 0.046 英寸 / 字符，面板列距 1.75 英寸，
# 留出间隙后约 34 个字符；超过时在分隔符处折行，否则会压到右侧面板的标题上
TITLE_WIDTH = 34
TITLE_MAX_LINES = 4


def _text_width(text: str) -> int:
    return sum(2 if unicodedata.east_asian_width(ch) in "WF" else 1 for ch in text)


def _hard_wrap(text: str, width: int) -> list[str]:
    """单个片段超过宽度时折行：优先在空格处断开，没有空格时按显示宽度逐字断开。"""
    if _text_width(text) <= width:
        return [text]
    lines, current = [], ""
    for word in text.split(" "):
        candidate = f"{current} {word}" if current else word
        if _text_width(candidate) <= width:
            current = candidate
            continue
        if current:
            lines.append(current)
        current = ""
        for ch in word:
            if _text_width(current + ch) > width:
                lines.append(current)
                current = ""
            current += ch
    if current:
        lines.append(current)
    return lines


def _pack(parts: list[str], sep: str, width: int = TITLE_WIDTH) -> list[str]:
    """把片段用 ``sep`` 依次拼接，超过 ``width`` 时另起一行；片段本身过长时折行。"""
    lines, current = [], ""
    for part in (p for p in parts if p):
        candidate = f"{current}{sep}{part}" if current else part
        if current and _text_width(candidate) > width:
            lines.append(current)
            current = part
        else:
            current = candidate
    if current:
        lines.append(current)
    return [piece for line in lines for piece in _hard_wrap(line, width)]


def _panel_title(row: pd.Series) -> str:
    """成员名与来源一行、状态与统计量一行；任一行超过面板宽度时在分隔符处折行，最多 ``TITLE_MAX_LINES`` 行。"""
    name = row["member_id"] or row.get("feature_name") or ""
    if row["layer"] == "mfa":
        name = format_tracer_display(str(row["tracer"]))
    tag = {"proteomics": "P", "transcriptomics": "R", "metabolomics": "M", "derived": "ratio", "mfa": "flux"}
    origin = f"{tag.get(row['layer'], row['layer'])}:{row['source']}"
    if row["role"] in ("regulator", "effector") and row.get("regulator_effect"):
        origin += f" ({row['regulator_effect']})"
    stats = [str(row["cell_state"])]
    if pd.notna(row.get("log2fc")):
        stats.append(f"log2FC {row['log2fc']:+.2f}")
    if pd.notna(row.get("padj")):
        stats.append(f"q {row['padj']:.2g}")
    elif pd.notna(row.get("pvalue")):
        stats.append(f"p {row['pvalue']:.2g}")
    if row["layer"] == "mfa" and isinstance(row.get("identifiability"), str) and row["identifiability"]:
        stats.append(row["identifiability"])
    lines = _pack([str(name), origin], " · ") + _pack(stats, "  ")
    if len(lines) > TITLE_MAX_LINES:
        lines = lines[:TITLE_MAX_LINES - 1] + [textwrap.shorten(" ".join(lines[TITLE_MAX_LINES - 1:]), TITLE_WIDTH,
                                                                placeholder="…")]
    return "\n".join(lines)


def _draw_member(
    ax: plt.Axes,
    row: pd.Series,
    data: pd.DataFrame,
    sides: dict[str, str],
    *,
    panel_aspect: float = 1.5,
    pal: SemanticPalette = METRO_ATLAS,
    preset: str = DEFAULT_PRESET,
) -> None:
    is_econ = preset == "economist"
    ax.set_title(_panel_title(row), fontsize=6, color=pal.ink, loc="left", pad=3)
    if is_econ:
        style_economist_axes(ax, right_y=False, zero_line=False, horizontal_grid=True)
        ax.tick_params(labelsize=6.0, colors=pal.ink2)
    else:
        ax.tick_params(labelsize=6.0, colors=pal.ink2, length=2)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        for spine in ("left", "bottom"):
            ax.spines[spine].set_color(pal.hairline)
    if data.empty or data["value"].notna().sum() == 0:
        ax.grid(False)
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.0)
        s = 0.22
        draw_status_cell(
            ax, 0.5 - s / 2, 0.68 - s * panel_aspect / 2, s, str(row["cell_state"] or "not_measured"), 0,
            palette=pal, font_scale=1.15, zorder=3.0, draw_box=True, y_up=True, y_scale=panel_aspect,
        )
        reason = str(row.get("state_reason") or row["cell_state"]).replace("_", " ")
        text = textwrap.fill(reason.strip(), 16, break_long_words=True)
        ax.text(0.5, 0.28, text, ha="center", va="center", fontsize=6.0, color=pal.muted)
        ax.set_xticks([])
        ax.set_yticks([])
        return
    ax.ticklabel_format(axis="y", useOffset=False, style="plain")
    groups = data[["group", "group_order"]].drop_duplicates().sort_values("group_order")
    positions = {g: i for i, g in enumerate(groups["group"])}
    rng = np.random.default_rng(0)
    for group, part in data.groupby("group", sort=False):
        xs = positions[group] + rng.uniform(-0.12, 0.12, len(part))
        colors = [pal.numerator if sides.get(s) == "numerator" else pal.denominator if sides.get(s) == "denominator"
                  else pal.neutral for s in part["sample_id"]]
        if row["layer"] == "mfa" and "fit_low" in part.columns:
            # 通量：误差线为该动物全部拟合点的四分位距；贴在盒约束边界上的动物画空心点
            low, high = pd.to_numeric(part["fit_low"], errors="coerce"), pd.to_numeric(part["fit_high"], errors="coerce")
            spread = low.notna() & high.notna()
            ax.vlines(xs[spread.to_numpy()], low[spread], high[spread], colors=[c for c, s in zip(colors, spread) if s],
                      linewidths=0.65, alpha=0.75, zorder=2)
            pinned = part["at_bound"].eq(True).to_numpy(dtype=bool) if "at_bound" in part.columns \
                else np.zeros(len(part), dtype=bool)
            ax.scatter(xs[~pinned], part["value"].to_numpy()[~pinned], s=14, color=np.array(colors, dtype=object)[~pinned],
                       edgecolors=pal.surface, linewidths=0.5, zorder=3)
            ax.scatter(xs[pinned], part["value"].to_numpy()[pinned], s=16, facecolors=pal.surface,
                       edgecolors=np.array(colors, dtype=object)[pinned], linewidths=0.95, zorder=3)
        else:
            ax.scatter(xs, part["value"], s=14, color=colors, edgecolors=pal.surface, linewidths=0.5, zorder=3)
        median = part["value"].median()
        if pd.notna(median):
            ax.plot([positions[group] - 0.26, positions[group] + 0.26], [median, median],
                    color=pal.ink if is_econ else pal.ink2, lw=1.0 if is_econ else 0.8, zorder=2.5)
    values = data["value"].dropna()
    center, span = float(values.median()), float(values.max() - values.min())
    if span <= 1e-6 * max(1.0, abs(center)):        # 全部相同（如都贴在边界上）：不放大浮点误差
        pad = max(1.0, 0.1 * abs(center))
        ax.set_ylim(center - pad, center + pad)
    scale = data["scale"].iloc[0]
    if is_econ and scale == "linear":
        ymin, ymax = ax.get_ylim()
        if ymin <= 0.0 <= ymax:
            ax.axhline(0.0, color=pal.ink, lw=0.75, ls="-", zorder=1.6)
    ax.set_xticks(range(len(positions)))
    ax.set_xticklabels(list(positions), rotation=45 if max(map(len, positions)) > 4 else 0, ha="right"
                       if max(map(len, positions)) > 4 else "center")
    ax.set_xlim(-0.6, len(positions) - 0.4)
    units = data["units"].iloc[0] if "units" in data and data["units"].iloc[0] not in ("", "TBD") else ""
    ax.set_ylabel("net flux" + (f" ({units})" if units else " (units TBD)") if scale == "linear" else "log2 value",
                  fontsize=6.0, color=pal.ink2)


def plot_dossier_page(rows: pd.DataFrame, samples: pd.DataFrame, sides: pd.DataFrame, reaction: str, tissue: str,
                      contrast: str, stem: Path, formats: Iterable[str] = ("pdf", "png"),
                      runs: dict[str, str] | None = None, ncols: int = 6,
                      *, style: str = DEFAULT_PRESET) -> list[Path]:
    """图 B：``rows`` 为该（反应, 组织, 比较）的 dossier_members（全部来源）。"""
    pal = get_palette(style)
    sections = _page_items(rows)
    side_map = {(s.source, s.layer, s.sample_id): s.side
                for s in sides[sides["contrast_id"] == contrast].itertuples()} if not sides.empty else {}
    order = list(dict.fromkeys(rows["source"]))
    # 每个分区内，每个来源另起一行（行首写来源名），两个来源的成员上下对齐
    grouped = [(title, [(source, [r for r in items if r["source"] == source]) for source in order
                        if any(r["source"] == source for r in items)]) for title, items in sections]

    def extra(items: list[pd.Series]) -> float:
        """组名较长（旋转后的刻度标签较高）的行多留一些底部空间。"""
        longest = max((len(str(g)) for r in items for g in _member_samples(r, samples)["group"].unique()), default=0)
        return 0.0 if longest <= 6 else min(0.55, 0.045 * longest)

    width, panel_h, left = 0.6 + 1.75 * ncols, 1.75, 1.05
    row_heights = [math.ceil(len(items) / ncols) * panel_h + extra(items) for _, groups in grouped for _, items in groups]
    top_pad = 1.02 if style == "economist" else 0.9
    height = top_pad + sum(row_heights) + 0.32 * len(sections) + 0.6
    with themed(preset=style):
        fig = plt.figure(figsize=(width, height), facecolor=pal.surface)
        if style == "economist":
            add_economist_tag(fig, x=0.01, y=1 - 0.06 / height, width=min(0.055, 0.55 / width),
                              height=0.085 / height, rule=True, rule_right=0.99)
        y = height - top_pad
        for title, groups in grouped:
            fig.text(0.01, y / height, title, fontsize=8, color=pal.ink, fontweight="bold", va="top")
            y -= 0.32
            for source, items in groups:
                fig.text(0.18 / width, (y - panel_h / 2 + 0.2) / height, source, fontsize=9, color=pal.ink2,
                         fontweight="bold", ha="center", va="center")
                for i, row in enumerate(items):
                    r, c = divmod(i, ncols)
                    top = y - r * panel_h
                    ax = fig.add_axes(((left + c * 1.75) / width, (top - panel_h + 0.55) / height, 1.2 / width,
                                       (panel_h - 0.95) / height))
                    ax.set_facecolor(pal.surface)
                    sides_here = {sample: side for (src, layer, sample), side in side_map.items()
                                  if src == row["source"] and layer == row["layer"]}
                    _draw_member(
                        ax, row, _member_samples(row, samples), sides_here,
                        panel_aspect=1.2 / (panel_h - 0.95), pal=pal, preset=style,
                    )
                y -= math.ceil(len(items) / ncols) * panel_h + extra(items)
        direction = rows["direction_evidence"].dropna().iloc[0] if rows["direction_evidence"].notna().any() else ""
        ty = 1 - (0.28 if style == "economist" else 0.25) / height
        sy = 1 - (0.58 if style == "economist" else 0.55) / height
        fig.text(0.01, ty, f"{reaction}  ·  {tissue}  ·  {contrast}", fontsize=10.5, color=pal.ink,
                 fontweight="bold", va="top")
        fig.text(0.01, sy, f"reaction semantics: {direction}   ·   per-animal points; bar = group median",
                 fontsize=6.5, color=pal.ink2, va="top")
        handles = [Line2D([], [], marker="o", ls="", color=pal.numerator, markersize=4.5),
                   Line2D([], [], marker="o", ls="", color=pal.denominator, markersize=4.5),
                   Line2D([], [], marker="o", ls="", color=pal.neutral, markersize=4.5)]
        fig.legend(handles, [f"{contrast}: numerator (exp)", "denominator (ctrl)", "other groups"], loc="upper right",
                   bbox_to_anchor=(0.99, ty + 0.04 / height), ncol=3, fontsize=6, frameon=False)
        if runs:
            fig.text(0.01, 0.12 / height, "Sources: " + "; ".join(f"{k} = {v}" for k, v in runs.items())
                     + "   values: log2 (CPM + 1 for RNA), net flux linear; flux whiskers = IQR of that animal's fits, "
                       "hollow = at a box-constraint bound (from sample_flux_fits; absent in older bundles)",
                     fontsize=6.0, color=pal.muted)
        return _save(fig, stem, formats, surface=pal.surface)


# ---------------------------------------------------------------------------
# 图 D
# ---------------------------------------------------------------------------

SPLIT_LAYER_NAMES = {"mfa": "flux", "proteomics": "enzyme P", "transcriptomics": "enzyme R", "metabolomics": "metabolite"}
SPLIT_STATES = {"not_identifiable": "⊘", "insufficient_n": "⚠", "not_measured": "/", "not_detected": "◌",
                "same_members": "—"}


def _split_title(row: pd.Series) -> str:
    head = SPLIT_LAYER_NAMES.get(row["layer"], row["layer"]) + (f" · {format_tracer_display(str(row['tracer']))}" if row["tracer"] else "")
    if row["state"] == "computed":
        stats = f"Δ log2(a/b) {row['log_ratio_change']:+.2f}  p {row['pvalue']:.2g}"
        if pd.notna(row.get("padj")):
            stats += f"  q {row['padj']:.2g}"
    else:
        stats = str(row["state"]).replace("_", " ")
    members = f"a: {row['members_a'] or '—'} | b: {row['members_b'] or '—'}"
    return f"{head}\n{stats}\n{textwrap.shorten(members, 46, placeholder='…')}"


def _draw_split(
    ax: plt.Axes,
    row: pd.Series,
    data: pd.DataFrame,
    sides: dict[str, str],
    *,
    panel_aspect: float = 1.37,
    pal: SemanticPalette = METRO_ATLAS,
    preset: str = DEFAULT_PRESET,
) -> None:
    is_econ = preset == "economist"
    ax.set_title(_split_title(row), fontsize=6.0, color=pal.ink, loc="left", pad=3)
    if is_econ:
        style_economist_axes(ax, right_y=False, zero_line=False, horizontal_grid=True)
        ax.tick_params(labelsize=6.0, colors=pal.ink2)
    else:
        ax.tick_params(labelsize=6.0, colors=pal.ink2, length=2)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        for spine in ("left", "bottom"):
            ax.spines[spine].set_color(pal.hairline)
    if data.empty or data["value"].notna().sum() == 0:
        ax.grid(False)
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.0)
        st = str(row["state"] or "not_measured")
        icon_st = "unchanged" if st == "same_members" else st
        s = 0.21
        draw_status_cell(
            ax, 0.5 - s / 2, 0.68 - s * panel_aspect / 2, s, icon_st, 0,
            palette=pal, font_scale=1.15, zorder=3.0, draw_box=True, y_up=True, y_scale=panel_aspect,
        )
        detail = str(row.get("detail") or st).replace("_", " ")
        text = textwrap.fill(detail.strip(), 22, break_long_words=False)
        ax.text(0.5, 0.27, text, ha="center", va="center", fontsize=6.0, color=pal.muted)
        ax.set_xticks([])
        ax.set_yticks([])
        return
    groups = data[["group", "group_order"]].drop_duplicates().sort_values("group_order")
    positions = {g: i for i, g in enumerate(groups["group"])}
    rng = np.random.default_rng(0)
    for group, part in data.groupby("group", sort=False):
        xs = positions[group] + rng.uniform(-0.12, 0.12, len(part))
        colors = [pal.numerator if sides.get(s) == "numerator" else pal.denominator if sides.get(s) == "denominator"
                  else pal.neutral for s in part["sample_id"]]
        ax.scatter(xs, part["value"], s=14, color=colors, edgecolors=pal.surface, linewidths=0.5, zorder=3)
        median = part["value"].median()
        if pd.notna(median):
            ax.plot([positions[group] - 0.26, positions[group] + 0.26], [median, median],
                    color=pal.ink if is_econ else pal.ink2, lw=1.0 if is_econ else 0.8, zorder=2.5)
    ax.set_xticks(range(len(positions)))
    longest = max(map(len, positions))
    ax.set_xticklabels(list(positions), rotation=45 if longest > 4 else 0, ha="right" if longest > 4 else "center")
    ax.set_xlim(-0.6, len(positions) - 0.4)
    ax.set_ylabel("log2(a / b)", fontsize=6.0, color=pal.ink2)


def plot_node_split(splits: pd.DataFrame, per_animal: pd.DataFrame, sides: pd.DataFrame, panel: Panel, tissue: str,
                    contrast: str, stem: Path, formats: Iterable[str] = ("pdf", "png"),
                    runs: dict[str, str] | None = None, *, style: str = DEFAULT_PRESET) -> list[Path]:
    """图 D：每个（节点, 出口对, 来源）一行，各层（通量按示踪剂）的逐只动物分流比并排。``splits`` 为该面板、组织、
    比较的 node_split 行；不能计算的格子保留面板并写明原因。"""
    pal = get_palette(style)
    layer_rank = {layer: i for i, layer in enumerate(SPLIT_LAYER_NAMES)}
    columns = sorted({(r.layer, r.tracer) for r in splits.itertuples()},
                     key=lambda c: (layer_rank.get(c[0], 9), c[1]))
    keys = list(dict.fromkeys((r.node, r.exit_a, r.exit_b, r.source) for r in splits.itertuples()))
    side_map = {(s.source, s.layer, s.sample_id): s.side
                for s in sides[sides["contrast_id"] == contrast].itertuples()} if not sides.empty else {}
    panel_w, panel_h, left = 1.85, 2.0, 1.45
    panel_aspect = (panel_w - 0.55) / (panel_h - 1.05)
    width = max(6.8, left + panel_w * len(columns) + 0.3)
    top_pad = 1.10 if style == "economist" else 1.0
    height = top_pad + panel_h * len(keys) + 0.72
    with themed(preset=style):
        fig = plt.figure(figsize=(width, height), facecolor=pal.surface)
        if style == "economist":
            add_economist_tag(fig, x=0.01, y=1 - 0.06 / height, width=min(0.065, 0.50 / width),
                              height=0.085 / height, rule=True, rule_right=0.99)
        for i, (node, exit_a, exit_b, source) in enumerate(keys):
            top = height - top_pad - i * panel_h
            fig.text(0.1 / width, (top - panel_h / 2 + 0.25) / height, f"{node}\na = {exit_a}\nb = {exit_b}\n{source}",
                     fontsize=7, color=pal.ink, va="center", ha="left", linespacing=1.4)
            for j, (layer, tracer) in enumerate(columns):
                hit = splits[(splits["node"] == node) & (splits["exit_a"] == exit_a) & (splits["exit_b"] == exit_b)
                             & (splits["source"] == source) & (splits["layer"] == layer) & (splits["tracer"] == tracer)]
                ax = fig.add_axes(((left + j * panel_w) / width, (top - panel_h + 0.5) / height, (panel_w - 0.55) / width,
                                   (panel_h - 1.05) / height))
                ax.set_facecolor(pal.surface)
                if hit.empty:              # 该来源没有这一层（或这种示踪剂）：保留面板并写明，不留空白
                    placeholder = pd.Series({"layer": layer, "tracer": tracer, "state": "not_measured",
                                             "members_a": "", "members_b": "",
                                             "detail": f"not in source {source}"})
                    _draw_split(ax, placeholder, per_animal.iloc[0:0], {}, panel_aspect=panel_aspect, pal=pal, preset=style)
                    continue
                data = per_animal[(per_animal["node"] == node) & (per_animal["exit_a"] == exit_a)
                                  & (per_animal["exit_b"] == exit_b) & (per_animal["source"] == source)
                                  & (per_animal["layer"] == layer) & (per_animal["tracer"] == tracer)]
                here = {sample: side for (src, lay, sample), side in side_map.items() if src == source and lay == layer}
                _draw_split(
                    ax, hit.iloc[0],
                    data if hit.iloc[0]["state"] == "computed" else data.iloc[0:0],
                    here, panel_aspect=panel_aspect, pal=pal, preset=style,
                )
        ty = 1 - (0.28 if style == "economist" else 0.25) / height
        sy = 1 - (0.58 if style == "economist" else 0.55) / height
        fig.text(0.01, ty, f"{panel.title}  ·  {tissue}  ·  {contrast}  ·  node split",
                 fontsize=10.5, color=pal.ink, fontweight="bold", va="top")
        fig.text(0.01, sy, "per animal: log2 of exit a over exit b (enzymes and metabolites summed on the "
                 "linear scale; flux = net flux of exit reactions consuming the node); bar = group median",
                 fontsize=6.5, color=pal.ink2, va="top")
        handles = [Line2D([], [], marker="o", ls="", color=c, markersize=4.5)
                   for c in (pal.numerator, pal.denominator, pal.neutral)]
        fig.legend(handles, [f"{contrast}: numerator (exp)", "denominator (ctrl)", "other groups"], loc="upper right",
                   bbox_to_anchor=(0.99, ty + 0.04 / height), ncol=3, fontsize=6, frameon=False)
        # 底部统一纯矢量状态图例
        lax = fig.add_axes((0.12 / width, 0.28 / height, (width - 0.24) / width, 0.34 / height))
        lax.set_facecolor(pal.surface)
        lax.set_xlim(0.0, width - 0.24)
        lax.set_ylim(0.0, 0.34)
        lax.axis("off")
        lx, sw = 0.02, 0.145
        for st, desc in (
            ("not_identifiable", "Unidentifiable"),
            ("insufficient_n", "Low-n"),
            ("not_measured", "Not measured"),
            ("not_detected", "Not detected"),
            ("unchanged", "Same members"),
        ):
            draw_status_cell(lax, lx, 0.17 - sw / 2, sw, st, 0, palette=pal, font_scale=0.9, zorder=2.5, draw_box=True, y_up=True)
            lax.text(lx + sw + 0.04, 0.17, desc, fontsize=6.0, color=pal.muted, va="center", ha="left")
            lx += sw + 0.06 + estimate_text_width(desc, 6.0, 1.0) + 0.18
        if runs:
            fig.text(0.01, 0.10 / height, "Sources: " + "; ".join(f"{k} = {v}" for k, v in runs.items()),
                     fontsize=6.0, color=pal.muted)
        return _save(fig, stem, formats, surface=pal.surface)


# ---------------------------------------------------------------------------
# 图 C
# ---------------------------------------------------------------------------

def plot_tracer_coverage(rows: pd.DataFrame, panel: Panel, tissue: str, contrast: str, stem: Path,
                         formats: Iterable[str] = ("pdf", "png"), tracers: Iterable[str] = (),
                         *, style: str = DEFAULT_PRESET) -> list[Path]:
    """图 C：面板反应 × 示踪剂的可辨识性（颜色），格内为该比较两组的动物数（分子 vs 分母）。

    ``rows``：本实验在该（组织, 比较）的分示踪剂通量行（dossier_members 中 ``layer`` 为 mfa、``tracer`` 不是 pooled），
    与图 A 的通量列是同一批行。"""
    pal = get_palette(style)
    part = rows[(rows["tissue"] == tissue) & (rows["contrast_id"] == contrast) & rows["layer"].eq("mfa")
                & rows["tracer"].ne("pooled") & rows["flux_id"].isin(panel.reactions)]
    order = list(dict.fromkeys([*tracers, *sorted(set(part["tracer"]))]))
    reactions = panel.reactions
    n_cols = max(len(order), 1)
    width = max(5.6, 2.2 + 0.78 * n_cols)
    height = (2.02 if style == "economist" else 1.85) + 0.34 * len(reactions)
    ax_w_in = width - 1.9
    ax_h_in = height - (1.82 if style == "economist" else 1.65)
    x_u = ax_w_in / n_cols
    y_u = ax_h_in / max(len(reactions), 1)
    y_scale = x_u / y_u
    with themed(preset=style):
        fig = plt.figure(figsize=(width, height), facecolor=pal.surface)
        if style == "economist":
            add_economist_tag(fig, x=0.015, y=1 - 0.06 / height, width=min(0.08, 0.45 / width),
                              height=0.08 / height, rule=True, rule_right=0.985)
        ax = fig.add_axes((1.6 / width, 1.15 / height, ax_w_in / width, ax_h_in / height))
        ax.set_facecolor(pal.surface)
        ax.grid(False)
        ax.set_xlim(-0.5, n_cols - 0.5)          # 没有示踪剂时保留一列宽的空图
        ax.set_ylim(len(reactions) - 0.5, -0.5)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([format_tracer_display(str(t)) for t in order], fontsize=6.4, color=pal.ink2)
        ax.set_yticks(range(len(reactions)))
        ax.set_yticklabels(reactions, fontsize=7, color=pal.ink)
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        cells = part.drop_duplicates(["flux_id", "tracer"]).set_index(["flux_id", "tracer"]) if not part.empty \
            else pd.DataFrame()
        for i, reaction in enumerate(reactions):
            for j, tracer in enumerate(order):
                key = (reaction, tracer)
                cell = cells.loc[key] if not cells.empty and key in cells.index else None
                if cell is None or cell["cell_state"] in ("not_infused", "not_measured"):
                    s = min(0.32, 0.58 / max(y_scale, 1.0))
                    draw_status_cell(
                        ax, j - s / 2, i - s * y_scale / 2, s, "not_measured", 0,
                        palette=pal, font_scale=1.0, zorder=2.2, draw_box=False, y_up=False, y_scale=y_scale,
                    )
                    continue
                level = cell["identifiability"] if isinstance(cell["identifiability"], str) else ""
                ax.add_patch(plt.Rectangle((j - 0.45, i - 0.4), 0.9, 0.8, color=IDENTIFIABILITY.get(level, pal.hairline),
                                           lw=0))
                n = (cell.get("n_numerator"), cell.get("n_denominator"))
                text = f"{int(n[0])} vs {int(n[1])}" if all(pd.notna(v) for v in n) else ""
                ink = pal.card_bg if level in ("well_constrained", "moderate") else pal.ink2
                flag_state = (
                    "insufficient_n" if cell["cell_state"] == "insufficient_n"
                    else ("bound_pinned" if level == "bound_pinned" else "")
                )
                if flag_state:
                    s = min(0.24, 0.54 / max(y_scale, 1.0))
                    draw_status_cell(
                        ax, j - 0.37, i - s * y_scale / 2, s, flag_state, 0,
                        palette=pal, font_scale=0.95, zorder=3.0, draw_box=False, y_up=False, y_scale=y_scale,
                    )
                    ax.text(j + 0.08, i, text, ha="center", va="center", fontsize=6.0, color=ink, zorder=3.0)
                else:
                    ax.text(j, i, text, ha="center", va="center", fontsize=6.0, color=ink, zorder=3.0)
        for i, reaction in enumerate(reactions):       # 面板声明中预期能约束该反应的示踪剂
            for tracer in panel.expected_tracers.get(reaction, []):
                if tracer in order:
                    j = order.index(tracer)
                    ax.add_patch(plt.Rectangle((j - 0.48, i - 0.44), 0.96, 0.88, fill=False, lw=0.9,
                                               ls=(0, (2, 1.5)), ec=pal.ink2, zorder=3))
        fig.text(0.015, 1 - (0.24 if style == "economist" else 0.18) / height,
                 f"Tracer coverage  ·  {panel.title}  ·  {tissue}  ·  {contrast}", fontsize=9.5,
                 color=pal.ink, fontweight="bold", va="top")
        handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in IDENTIFIABILITY.values()]
        fig.legend(handles, [k.replace("_", " ") for k in IDENTIFIABILITY], loc="lower left",
                   bbox_to_anchor=(0.015, 0.06 / height), ncol=4, fontsize=6.0, frameon=False)
        # 纯矢量图标图例条（替代 x / ! / / 文本）
        lax = fig.add_axes((0.12 / width, 0.52 / height, (width - 0.24) / width, 0.36 / height))
        lax.set_facecolor(pal.surface)
        lax.set_xlim(0.0, width - 0.24)
        lax.set_ylim(0.0, 0.36)
        lax.axis("off")
        lx, sw = 0.02, 0.15
        for st, desc in (
            ("not_measured", "Not infused / measured"),
            ("bound_pinned", "Bound-pinned"),
            ("insufficient_n", "Low-n (<2)"),
        ):
            draw_status_cell(lax, lx, 0.18 - sw / 2, sw, st, 0, palette=pal, font_scale=0.92, zorder=2.5, draw_box=True, y_up=True)
            lax.text(lx + sw + 0.04, 0.18, desc, fontsize=6.0, color=pal.muted, va="center", ha="left")
            lx += sw + 0.06 + estimate_text_width(desc, 6.0, 1.0) + 0.16
        lax.add_patch(Rectangle((lx, 0.18 - sw * 0.42), sw, sw * 0.84, fill=False, lw=0.85, ls=(0, (2, 1.5)), ec=pal.ink2))
        lax.text(lx + sw + 0.05, 0.18, "Expected tracer  ·  a vs b = num vs den", fontsize=6.0, color=pal.muted, va="center", ha="left")
        return _save(fig, stem, formats, surface=pal.surface)
