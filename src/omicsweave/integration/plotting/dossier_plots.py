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
import textwrap
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D

from omicsweave.visualization.theme import themed

from ..dossier.panels import Panel

INK, INK2, MUTED, HAIRLINE, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#fcfcfb"
NUMERATOR, DENOMINATOR = "#eb6834", "#2a78d6"            # 分类色第 2、1 位
DIVERGING = LinearSegmentedColormap.from_list(
    "dossier_diverging", ["#104281", "#2a78d6", "#9ec5f4", "#f0efec", "#f5b5b4", "#e34948", "#8c1f1f"])
VMAX = 2.0
NORM = Normalize(-VMAX, VMAX, clip=True)
GLYPHS = {"equivalent": "=", "not_identifiable": "x", "insufficient_n": "!", "not_measured": "/",
          "not_infused": "/", "not_detected": "?", "no_knowledge": "–", "qc_blocked": "QC"}
GLYPH_LEGEND = [("=", "equivalent (90% CI within ±0.58)"), ("x", "not identifiable"), ("!", "insufficient n"),
                ("/", "not measured / not infused"), ("?", "mapped, not detected"), ("–", "no member in knowledge"),
                ("QC", "QC blocked")]
BLOCK_TITLES = {"flux": "FLUX by tracer", "enzyme": "ENZYME", "regulator": "REGULATOR", "pool": "POOL",
                "assessment": "ASSESSMENT"}
# 评估块：形状表示方向与强度（实心 = changed / established，空心 = suggestive），颜色与色阶的升降一致
AXIS_GLYPHS = {"unchanged": "=", "mixed": "≠", "undetermined": "·", "not_identifiable": "x", "insufficient_n": "!",
               "not_detected": "?", "not_measured": "/", "qc_blocked": "QC", "no_knowledge": "–"}
AXIS_LEGEND = [("filled ▲▼", "changed (flux: established)"), ("hollow", "suggestive"), ("=", "unchanged"),
               ("≠", "mixed / tracers conflict"), ("·", "undetermined"), ("L3–L0, LX", "claim level (n concordant axes)")]
IDENTIFIABILITY = {"well_constrained": "#1c5cab", "moderate": "#5598e7", "poorly_constrained": "#b7d3f6",
                   "bound_pinned": "#e1e0d9"}


def _size(certainty: float) -> float:
    if certainty is None or not np.isfinite(certainty):
        return 70.0
    return 25.0 + 175.0 * float(np.clip((certainty - 0.5) / 0.5, 0, 1))


def _save(fig: plt.Figure, stem: Path, formats: Iterable[str]) -> list[Path]:
    stem.parent.mkdir(parents=True, exist_ok=True)
    paths = []
    for fmt in formats:
        path = stem.with_suffix(f".{fmt}")
        try:
            fig.savefig(path, bbox_inches="tight", facecolor=SURFACE)
        except (AttributeError, KeyError, ValueError):
            if fmt != "pdf":
                raise
            # 部分 fontTools 版本对中文字体（微软雅黑、黑体）做 Type 42 子集化时失败：改用 Type 3 嵌入（文字不可编辑）
            with plt.rc_context({"pdf.fonttype": 3}):
                fig.savefig(path, bbox_inches="tight", facecolor=SURFACE)
        paths.append(path)
    plt.close(fig)
    return paths


def _draw_cell(ax: plt.Axes, x: float, y: float, cell: pd.Series) -> None:
    """一个格子：点在上，log2FC 写在点的下方；没有测到的状态只画符号。"""
    state = cell["display_state"]
    if state in GLYPHS and state != "equivalent":
        ax.text(x, y, GLYPHS[state], ha="center", va="center", fontsize=7.5, color=MUTED)
        return
    value, cy = cell["value"], y + 0.13
    if state == "equivalent":
        ax.scatter([x], [cy], s=_size(cell["certainty"]), color="#f0efec", edgecolors=HAIRLINE, linewidths=0.6,
                   zorder=2)
        ax.text(x, cy, "=", ha="center", va="center", fontsize=7, color=INK2, zorder=3)
    else:
        face = DIVERGING(NORM(value)) if pd.notna(value) else "none"
        ring = cell["ring"]
        edge, width = (INK, 1.6) if ring == "selected" else (INK, 0.8) if ring == "changed" else (MUTED, 0.4)
        ax.scatter([x], [cy], s=_size(cell["certainty"]), color=face, edgecolors=edge, linewidths=width, zorder=2)
    if pd.notna(value):
        ax.text(x, y - 0.3, f"{value:+.2f}", ha="center", va="center", fontsize=5.2, color=INK2, zorder=3)
    if state == "summary" and cell["n_quantified"]:
        ax.text(x + 0.3, y + 0.38, f"{int(cell['n_changed'])}/{int(cell['n_quantified'])}", ha="left", va="bottom",
                fontsize=5, color=INK2, zorder=3)


def _draw_assessment(ax: plt.Axes, x: float, y: float, cell: pd.Series) -> None:
    """评估块的一格：四条轴为方向符号，结论等级为文字。"""
    state, direction = cell["display_state"], int(cell["value"]) if pd.notna(cell["value"]) else 0
    if cell["column"] == "level":
        strong = state in ("L1", "L2", "L3")
        ax.text(x, y, str(cell["text"]), ha="center", va="center", fontsize=7.5 if strong else 7,
                color=INK if strong or state == "LX" else MUTED, fontweight="bold" if strong else "normal")
        return
    if state in ("changed", "suggestive") and direction:
        # 用标记而不是字形画三角：实心 = changed，空心 = suggestive，大小一致（字体回退会让空心字形变小）
        color = "#e34948" if direction > 0 else "#2a78d6"
        ax.scatter([x], [y], s=55, marker="^" if direction > 0 else "v", linewidths=1.2, zorder=3,
                   facecolors=color if state == "changed" else "none", edgecolors=color)
        return
    ax.text(x, y, AXIS_GLYPHS.get(state, "·"), ha="center", va="center", fontsize=8,
            color=INK2 if state in ("unchanged", "mixed") else MUTED)


def _short_label(block: str, column: str, tracer: str) -> str:
    """列的短名：通量列为示踪剂，其余为层的代码（P、R、M）；来源写在上方的来源行。"""
    return str(tracer) if block == "flux" else str(column).split(":")[0]


def plot_matrix(view: pd.DataFrame, panel: Panel, tissue: str, contrast: str, stem: Path,
                formats: Iterable[str] = ("pdf", "png"), runs: dict[str, str] | None = None,
                caption: str = "") -> list[Path]:
    """图 A：只画 ``view``（矩阵视图）中的格子。"""
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
                step = 1.8 if (block, r.column) == ("assessment", "level") else 1.0
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
    width = 2.6 + 0.5 * (right + 0.7)
    row_h = 0.40
    height = 1.9 + row_h * len(rows) + 1.9
    with themed():
        fig = plt.figure(figsize=(width, height), facecolor=SURFACE)
        ax = fig.add_axes((2.4 / width, 1.9 / height, (width - 2.6) / width, (height - 3.8) / height))
        ax.set_facecolor(SURFACE)
        ax.set_xlim(-0.7, right)
        ax.set_ylim(-len(rows) + 0.4, 0.8)
        ax.axis("off")
        previous_segment = None
        for _, row in rows.iterrows():
            y = row_y[row["row"]]
            reaction = row["element_type"] == "reaction"
            ax.text(-0.9, y, f"{row['label']}", ha="right", va="center", fontsize=8 if reaction else 7,
                    color=INK if reaction else INK2, fontweight="bold" if reaction else "normal",
                    transform=ax.transData, clip_on=False)
            if reaction:
                ax.add_patch(plt.Rectangle((-0.7, y - 0.5), last + 1.2, 1.0, color="#f4f3ef", zorder=0, lw=0))
            if row["segment"] != previous_segment:
                ax.plot([-0.7, last + 0.5], [y + 0.5, y + 0.5], color=HAIRLINE, lw=0.6, zorder=1)
                ax.text(last + 0.8, y + 0.45, row["segment"], ha="left", va="top", fontsize=6.5, color=MUTED)
                previous_segment = row["segment"]
        for key, cx in positions.items():
            name = labels[key]
            ax.text(cx, 0.65, name, ha="center", va="bottom", rotation=90 if len(name) > 2 else 0, fontsize=6.5,
                    color=INK2)
        for start, end, source in source_spans:      # 来源行：来源名居中，下划细线
            ax.plot([start - 0.35, end + 0.35], [2.05, 2.05], color=MUTED, lw=0.6, clip_on=False)
            ax.text((start + end) / 2, 2.12, source, ha="center", va="bottom", fontsize=7, color=INK,
                    clip_on=False)
        for block, (start, end) in block_spans.items():
            ax.text((start + end) / 2, 2.75, BLOCK_TITLES[block], ha="center", va="bottom", fontsize=7.5, color=INK,
                    fontweight="bold", clip_on=False)
        for _, cell in view.iterrows():
            cx, cy = positions[(cell["block"], cell["column"])], row_y[cell["row"]]
            if cell["block"] == "flux" and bool(cell["expected_tracer"]):
                ax.add_patch(plt.Rectangle((cx - 0.45, cy - 0.45), 0.9, 0.9, fill=False, lw=0.6, ls=(0, (2, 1.5)),
                                           ec=MUTED, zorder=1))
            if cell["block"] == "assessment":
                _draw_assessment(ax, cx, cy, cell)
            else:
                _draw_cell(ax, cx, cy, cell)
        fig.text(0.01, 1 - 0.35 / height, f"{panel.title}  ·  {tissue}  ·  {contrast}", fontsize=10, color=INK,
                 fontweight="bold", va="top")
        _matrix_legend(fig, width, height, runs or {})
        if caption:                                  # 校准结果（E8）写入图注
            fig.text(0.04, 0.04 / height, caption, fontsize=6, color=INK2)
        return _save(fig, stem, formats)


def _matrix_legend(fig: plt.Figure, width: float, height: float, runs: dict[str, str]) -> None:
    cax = fig.add_axes((0.04, 0.9 / height, 1.6 / width, 0.12 / height))
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=NORM, cmap=DIVERGING), cax=cax, orientation="horizontal")
    bar.set_ticks([-2, -1, 0, 1, 2])
    bar.ax.tick_params(labelsize=6, colors=INK2, length=2)
    bar.outline.set_edgecolor(HAIRLINE)
    cax.set_title("log2FC (regulators: × activity sign)", fontsize=6.5, color=INK2, pad=3)
    handles = [plt.scatter([], [], s=_size(c), color="#d9d8d3", edgecolors=MUTED, linewidths=0.4)
               for c in (0.6, 0.8, 1.0)]
    handles += [plt.scatter([], [], s=90, color="#d9d8d3", edgecolors=INK, linewidths=w) for w in (0.8, 1.6)]
    handles += [plt.scatter([], [], s=90, color="none", edgecolors=MUTED, linewidths=0.6)]
    labels = ["P(sign) 0.6", "0.8", "1.0", "FDR < 0.05", "FDR and effect threshold",
              "measured, no signed regulator"]
    fig.legend(handles, labels, loc="lower left", bbox_to_anchor=(2.1 / width, 0.55 / height), ncol=6, fontsize=6,
               frameon=False, handletextpad=0.3, columnspacing=1.0)
    glyphs = "   ".join(f"{g} {label}" for g, label in GLYPH_LEGEND)
    fig.text(0.04, 0.38 / height, glyphs + "   k/n changed / quantified members   dashed box: expected tracer",
             fontsize=6, color=INK2)
    fig.text(0.04, 0.27 / height, "ASSESSMENT  F flux · C enzyme capacity · R regulation · M substrates/products:   "
             + "   ".join(f"{g} {label}" for g, label in AXIS_LEGEND), fontsize=6, color=INK2)
    if runs:
        fig.text(0.04, 0.15 / height, "Sources: " + "; ".join(f"{k} = {v}" for k, v in runs.items())
                 + "   P protein · R RNA · M metabolite (pool: metabolite or product/substrate ratio)"
                 + "   number under each point: log2FC (cell summary for multi-member cells)",
                 fontsize=6, color=MUTED)


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


def _panel_title(row: pd.Series) -> str:
    name = row["member_id"] or row.get("feature_name") or ""
    if row["layer"] == "mfa":
        name = str(row["tracer"])
    tag = {"proteomics": "P", "transcriptomics": "R", "metabolomics": "M", "derived": "ratio", "mfa": "flux"}
    head = f"{name} · {tag.get(row['layer'], row['layer'])}:{row['source']}"
    if row["role"] in ("regulator", "effector") and row.get("regulator_effect"):
        head += f" ({row['regulator_effect']})"
    stats = []
    if pd.notna(row.get("log2fc")):
        stats.append(f"log2FC {row['log2fc']:+.2f}")
    if pd.notna(row.get("padj")):
        stats.append(f"q {row['padj']:.2g}")
    elif pd.notna(row.get("pvalue")):
        stats.append(f"p {row['pvalue']:.2g}")
    if row["layer"] == "mfa" and isinstance(row.get("identifiability"), str) and row["identifiability"]:
        stats.append(row["identifiability"])
    return head + "\n" + row["cell_state"] + ("  " + "  ".join(stats) if stats else "")


def _draw_member(ax: plt.Axes, row: pd.Series, data: pd.DataFrame, sides: dict[str, str]) -> None:
    ax.set_title(_panel_title(row), fontsize=6, color=INK, loc="left", pad=3)
    ax.tick_params(labelsize=5.5, colors=INK2, length=2)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(HAIRLINE)
    if data.empty or data["value"].notna().sum() == 0:
        reason = str(row.get("state_reason") or row["cell_state"])
        text = textwrap.fill(f"{GLYPHS.get(row['cell_state'], '')} {reason}".strip(), 14, break_long_words=True)
        ax.text(0.5, 0.5, text, ha="center", va="center", fontsize=5.5, color=MUTED, transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])
        return
    ax.ticklabel_format(axis="y", useOffset=False, style="plain")
    groups = data[["group", "group_order"]].drop_duplicates().sort_values("group_order")
    positions = {g: i for i, g in enumerate(groups["group"])}
    rng = np.random.default_rng(0)
    for group, part in data.groupby("group", sort=False):
        xs = positions[group] + rng.uniform(-0.12, 0.12, len(part))
        colors = [NUMERATOR if sides.get(s) == "numerator" else DENOMINATOR if sides.get(s) == "denominator"
                  else MUTED for s in part["sample_id"]]
        ax.scatter(xs, part["value"], s=12, color=colors, edgecolors=SURFACE, linewidths=0.5, zorder=3)
        median = part["value"].median()
        if pd.notna(median):
            ax.plot([positions[group] - 0.25, positions[group] + 0.25], [median, median], color=INK2, lw=0.8, zorder=2)
    values = data["value"].dropna()
    center, span = float(values.median()), float(values.max() - values.min())
    if span <= 1e-6 * max(1.0, abs(center)):        # 全部相同（如都贴在边界上）：不放大浮点误差
        pad = max(1.0, 0.1 * abs(center))
        ax.set_ylim(center - pad, center + pad)
    ax.set_xticks(range(len(positions)))
    ax.set_xticklabels(list(positions), rotation=45 if max(map(len, positions)) > 4 else 0, ha="right"
                       if max(map(len, positions)) > 4 else "center")
    ax.set_xlim(-0.6, len(positions) - 0.4)
    scale = data["scale"].iloc[0]
    units = data["units"].iloc[0] if "units" in data and data["units"].iloc[0] not in ("", "TBD") else ""
    ax.set_ylabel("net flux" + (f" ({units})" if units else " (units TBD)") if scale == "linear" else "log2 value",
                  fontsize=5.5, color=INK2)


def plot_dossier_page(rows: pd.DataFrame, samples: pd.DataFrame, sides: pd.DataFrame, reaction: str, tissue: str,
                      contrast: str, stem: Path, formats: Iterable[str] = ("pdf", "png"),
                      runs: dict[str, str] | None = None, ncols: int = 6) -> list[Path]:
    """图 B：``rows`` 为该（反应, 组织, 比较）的 dossier_members（全部来源）。"""
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
    height = 0.9 + sum(row_heights) + 0.32 * len(sections) + 0.6
    with themed():
        fig = plt.figure(figsize=(width, height), facecolor=SURFACE)
        y = height - 0.9
        for title, groups in grouped:
            fig.text(0.01, y / height, title, fontsize=8, color=INK, fontweight="bold", va="top")
            y -= 0.32
            for source, items in groups:
                fig.text(0.18 / width, (y - panel_h / 2 + 0.2) / height, source, fontsize=9, color=INK2,
                         fontweight="bold", ha="center", va="center")
                for i, row in enumerate(items):
                    r, c = divmod(i, ncols)
                    top = y - r * panel_h
                    ax = fig.add_axes(((left + c * 1.75) / width, (top - panel_h + 0.55) / height, 1.2 / width,
                                       (panel_h - 0.95) / height))
                    ax.set_facecolor(SURFACE)
                    sides_here = {sample: side for (src, layer, sample), side in side_map.items()
                                  if src == row["source"] and layer == row["layer"]}
                    _draw_member(ax, row, _member_samples(row, samples), sides_here)
                y -= math.ceil(len(items) / ncols) * panel_h + extra(items)
        direction = rows["direction_evidence"].dropna().iloc[0] if rows["direction_evidence"].notna().any() else ""
        fig.text(0.01, 1 - 0.25 / height, f"{reaction}  ·  {tissue}  ·  {contrast}", fontsize=10, color=INK,
                 fontweight="bold", va="top")
        fig.text(0.01, 1 - 0.55 / height, f"reaction semantics: {direction}   ·   per-animal points; bar = group median",
                 fontsize=6.5, color=INK2, va="top")
        handles = [Line2D([], [], marker="o", ls="", color=NUMERATOR, markersize=4),
                   Line2D([], [], marker="o", ls="", color=DENOMINATOR, markersize=4),
                   Line2D([], [], marker="o", ls="", color=MUTED, markersize=4)]
        fig.legend(handles, [f"{contrast}: numerator", "denominator", "other groups"], loc="upper right",
                   bbox_to_anchor=(0.99, 1 - 0.2 / height), ncol=3, fontsize=6, frameon=False)
        if runs:
            fig.text(0.01, 0.12 / height, "Sources: " + "; ".join(f"{k} = {v}" for k, v in runs.items())
                     + "   values: log2 (CPM + 1 for RNA), net flux linear; at-bound per animal not exported by MFA",
                     fontsize=5.5, color=MUTED)
        return _save(fig, stem, formats)


# ---------------------------------------------------------------------------
# 图 D
# ---------------------------------------------------------------------------

SPLIT_LAYER_NAMES = {"mfa": "flux", "proteomics": "enzyme P", "transcriptomics": "enzyme R", "metabolomics": "metabolite"}
SPLIT_STATES = {"not_identifiable": "x", "insufficient_n": "!", "not_measured": "/", "not_detected": "?",
                "same_members": "="}


def _split_title(row: pd.Series) -> str:
    head = SPLIT_LAYER_NAMES.get(row["layer"], row["layer"]) + (f" · {row['tracer']}" if row["tracer"] else "")
    if row["state"] == "computed":
        stats = f"Δ log2(a/b) {row['log_ratio_change']:+.2f}  p {row['pvalue']:.2g}"
        if pd.notna(row.get("padj")):
            stats += f"  q {row['padj']:.2g}"
    else:
        stats = row["state"]
    members = f"a: {row['members_a'] or '—'} | b: {row['members_b'] or '—'}"
    return f"{head}\n{stats}\n{textwrap.shorten(members, 46, placeholder='…')}"


def _draw_split(ax: plt.Axes, row: pd.Series, data: pd.DataFrame, sides: dict[str, str]) -> None:
    ax.set_title(_split_title(row), fontsize=5.8, color=INK, loc="left", pad=3)
    ax.tick_params(labelsize=5.5, colors=INK2, length=2)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(HAIRLINE)
    if data.empty or data["value"].notna().sum() == 0:
        text = textwrap.fill(f"{SPLIT_STATES.get(row['state'], '')} {row.get('detail') or row['state']}".strip(), 24,
                             break_long_words=False)
        ax.text(0.5, 0.5, text, ha="center", va="center", fontsize=5.5, color=MUTED, transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])
        return
    groups = data[["group", "group_order"]].drop_duplicates().sort_values("group_order")
    positions = {g: i for i, g in enumerate(groups["group"])}
    rng = np.random.default_rng(0)
    for group, part in data.groupby("group", sort=False):
        xs = positions[group] + rng.uniform(-0.12, 0.12, len(part))
        colors = [NUMERATOR if sides.get(s) == "numerator" else DENOMINATOR if sides.get(s) == "denominator"
                  else MUTED for s in part["sample_id"]]
        ax.scatter(xs, part["value"], s=12, color=colors, edgecolors=SURFACE, linewidths=0.5, zorder=3)
        median = part["value"].median()
        if pd.notna(median):
            ax.plot([positions[group] - 0.25, positions[group] + 0.25], [median, median], color=INK2, lw=0.8, zorder=2)
    ax.set_xticks(range(len(positions)))
    longest = max(map(len, positions))
    ax.set_xticklabels(list(positions), rotation=45 if longest > 4 else 0, ha="right" if longest > 4 else "center")
    ax.set_xlim(-0.6, len(positions) - 0.4)
    ax.set_ylabel("log2(a / b)", fontsize=5.5, color=INK2)


def plot_node_split(splits: pd.DataFrame, per_animal: pd.DataFrame, sides: pd.DataFrame, panel: Panel, tissue: str,
                    contrast: str, stem: Path, formats: Iterable[str] = ("pdf", "png"),
                    runs: dict[str, str] | None = None) -> list[Path]:
    """图 D：每个（节点, 出口对, 来源）一行，各层（通量按示踪剂）的逐只动物分流比并排。``splits`` 为该面板、组织、
    比较的 node_split 行；不能计算的格子保留面板并写明原因。"""
    layer_rank = {layer: i for i, layer in enumerate(SPLIT_LAYER_NAMES)}
    columns = sorted({(r.layer, r.tracer) for r in splits.itertuples()},
                     key=lambda c: (layer_rank.get(c[0], 9), c[1]))
    keys = list(dict.fromkeys((r.node, r.exit_a, r.exit_b, r.source) for r in splits.itertuples()))
    side_map = {(s.source, s.layer, s.sample_id): s.side
                for s in sides[sides["contrast_id"] == contrast].itertuples()} if not sides.empty else {}
    panel_w, panel_h, left = 1.85, 2.0, 1.45
    width, height = left + panel_w * len(columns) + 0.3, 1.0 + panel_h * len(keys) + 0.5
    with themed():
        fig = plt.figure(figsize=(width, height), facecolor=SURFACE)
        for i, (node, exit_a, exit_b, source) in enumerate(keys):
            top = height - 1.0 - i * panel_h
            fig.text(0.1 / width, (top - panel_h / 2 + 0.25) / height, f"{node}\na = {exit_a}\nb = {exit_b}\n{source}",
                     fontsize=7, color=INK, va="center", ha="left", linespacing=1.4)
            for j, (layer, tracer) in enumerate(columns):
                hit = splits[(splits["node"] == node) & (splits["exit_a"] == exit_a) & (splits["exit_b"] == exit_b)
                             & (splits["source"] == source) & (splits["layer"] == layer) & (splits["tracer"] == tracer)]
                ax = fig.add_axes(((left + j * panel_w) / width, (top - panel_h + 0.5) / height, (panel_w - 0.55) / width,
                                   (panel_h - 1.05) / height))
                ax.set_facecolor(SURFACE)
                if hit.empty:              # 该来源没有这一层（或这种示踪剂）：保留面板并写明，不留空白
                    placeholder = pd.Series({"layer": layer, "tracer": tracer, "state": "not_measured",
                                             "members_a": "", "members_b": "",
                                             "detail": f"not in source {source}"})
                    _draw_split(ax, placeholder, per_animal.iloc[0:0], {})
                    continue
                data = per_animal[(per_animal["node"] == node) & (per_animal["exit_a"] == exit_a)
                                  & (per_animal["exit_b"] == exit_b) & (per_animal["source"] == source)
                                  & (per_animal["layer"] == layer) & (per_animal["tracer"] == tracer)]
                here = {sample: side for (src, lay, sample), side in side_map.items() if src == source and lay == layer}
                _draw_split(ax, hit.iloc[0], data if hit.iloc[0]["state"] == "computed" else data.iloc[0:0], here)
        fig.text(0.01, 1 - 0.25 / height, f"{panel.title}  ·  {tissue}  ·  {contrast}  ·  node split",
                 fontsize=10, color=INK, fontweight="bold", va="top")
        fig.text(0.01, 1 - 0.55 / height, "per animal: log2 of exit a over exit b (enzymes and metabolites summed on the "
                 "linear scale; flux = net flux of exit reactions consuming the node); bar = group median",
                 fontsize=6.5, color=INK2, va="top")
        handles = [Line2D([], [], marker="o", ls="", color=c, markersize=4) for c in (NUMERATOR, DENOMINATOR, MUTED)]
        fig.legend(handles, [f"{contrast}: numerator", "denominator", "other groups"], loc="upper right",
                   bbox_to_anchor=(0.99, 1 - 0.2 / height), ncol=3, fontsize=6, frameon=False)
        legend = "   ".join(f"{g} {s.replace('_', ' ')}" for s, g in SPLIT_STATES.items())
        fig.text(0.01, 0.3 / height, legend + "   labelled-isotope split awaits G10", fontsize=5.5, color=MUTED)
        if runs:
            fig.text(0.01, 0.12 / height, "Sources: " + "; ".join(f"{k} = {v}" for k, v in runs.items()),
                     fontsize=5.5, color=MUTED)
        return _save(fig, stem, formats)


# ---------------------------------------------------------------------------
# 图 C
# ---------------------------------------------------------------------------

def plot_tracer_coverage(rows: pd.DataFrame, panel: Panel, tissue: str, contrast: str, stem: Path,
                         formats: Iterable[str] = ("pdf", "png"), tracers: Iterable[str] = ()) -> list[Path]:
    """图 C：面板反应 × 示踪剂的可辨识性（颜色），格内为该比较两组的动物数（分子 vs 分母）。

    ``rows``：本实验在该（组织, 比较）的分示踪剂通量行（dossier_members 中 ``layer`` 为 mfa、``tracer`` 不是 pooled），
    与图 A 的通量列是同一批行。"""
    part = rows[(rows["tissue"] == tissue) & (rows["contrast_id"] == contrast) & rows["layer"].eq("mfa")
                & rows["tracer"].ne("pooled") & rows["flux_id"].isin(panel.reactions)]
    order = list(dict.fromkeys([*tracers, *sorted(set(part["tracer"]))]))
    reactions = panel.reactions
    width, height = 2.2 + 0.75 * max(len(order), 1), 1.6 + 0.32 * len(reactions)
    with themed():
        fig = plt.figure(figsize=(width, height), facecolor=SURFACE)
        ax = fig.add_axes((1.6 / width, 1.0 / height, (width - 1.9) / width, (height - 1.5) / height))
        ax.set_xlim(-0.5, max(len(order), 1) - 0.5)          # 没有示踪剂时保留一列宽的空图
        ax.set_ylim(len(reactions) - 0.5, -0.5)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order, fontsize=7, color=INK2)
        ax.set_yticks(range(len(reactions)))
        ax.set_yticklabels(reactions, fontsize=7, color=INK)
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
                    ax.text(j, i, "/", ha="center", va="center", fontsize=8, color=MUTED)
                    continue
                level = cell["identifiability"] if isinstance(cell["identifiability"], str) else ""
                ax.add_patch(plt.Rectangle((j - 0.45, i - 0.4), 0.9, 0.8, color=IDENTIFIABILITY.get(level, HAIRLINE),
                                           lw=0))
                n = (cell.get("n_numerator"), cell.get("n_denominator"))
                text = f"{int(n[0])} vs {int(n[1])}" if all(pd.notna(v) for v in n) else ""
                if cell["cell_state"] == "insufficient_n":
                    text = "! " + text
                elif level == "bound_pinned":
                    text = "x " + text
                ink = SURFACE if level in ("well_constrained", "moderate") else INK2
                ax.text(j, i, text, ha="center", va="center", fontsize=6, color=ink)
        for i, reaction in enumerate(reactions):       # 面板声明中预期能约束该反应的示踪剂
            for tracer in panel.expected_tracers.get(reaction, []):
                if tracer in order:
                    j = order.index(tracer)
                    ax.add_patch(plt.Rectangle((j - 0.48, i - 0.44), 0.96, 0.88, fill=False, lw=0.9,
                                               ls=(0, (2, 1.5)), ec=INK2, zorder=3))
        fig.text(0.01, 1 - 0.2 / height, f"Tracer coverage  ·  {panel.title}  ·  {tissue}  ·  {contrast}", fontsize=9,
                 color=INK, fontweight="bold", va="top")
        handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in IDENTIFIABILITY.values()]
        fig.legend(handles, [k.replace("_", " ") for k in IDENTIFIABILITY], loc="lower left",
                   bbox_to_anchor=(0.01, 0.0), ncol=4, fontsize=6, frameon=False)
        fig.text(0.01, 0.62 / height, "a vs b = animals in numerator vs denominator   / not infused or not measured"
                 "   x bound-pinned   ! fewer than 2 animals in a group   dashed box: expected tracer",
                 fontsize=6, color=MUTED)
        return _save(fig, stem, formats)
