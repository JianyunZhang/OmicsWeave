#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
火山图：由统一差异结果表（``differential.parquet``，v2 设计文档第 8.2 节）画出，各组学共用，每个组织 × 比较一张。

- 横轴 log2FC，对称于 0，上调与下调的面积可以直接比较；纵轴为 −log10(校正后 p) 或 −log10(p)，自 0 起。
  一个面板中没有任何 p 值时（如 MFA 只给效应量的比较），纵轴改为 |效应量|。
- 点的颜色只由表中的 ``selected`` 与 log2FC 的符号决定，所以图与结果表的入选结果始终一致。阈值线只是参考：
  入选规则（``selection_rule``）写在图下方；规则中有阈值线表达不了的条件（如 VIP、同时要求 p 与 FDR、
  稳定性校准）时，线外的点也可能未入选。
- log2FC 为 ±inf 或超出 ``x_cap`` 的点钉在横轴边界上，用空心三角标出，不丢弃；log2FC 或纵轴指标为空的特征
  不画，个数写在图下方。
- 标注纵轴最高的若干个入选特征；装有 adjustText 时自动避让。

图中文字为英文，便于直接用于论文。
"""
from __future__ import annotations

import textwrap
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

from omicsweave.contract.tables import read_differential
from omicsweave.core.model import Contrast

from .figures import FIGURE_FORMATS, figure_stem, save_figure
from .theme import themed

KIND = "volcano"
COLORS = {"up": "#D32F2F", "down": "#1565C0", "ns": "#BDBDBD"}
_LEGEND = {"up": "Up", "down": "Down", "ns": "Not selected"}
_Y_LABELS = {"padj": "−log$_{10}$(adjusted p)", "pvalue": "−log$_{10}$(p)"}
_Y_NAMES = {"padj": "adjusted p", "pvalue": "p", "effect_size": "effect size"}
_EFFECT_NAMES = {"hedges_g": "Hedges' g", "cohen_d": "Cohen's d"}
_P_FLOOR = 1e-300
_LABEL_CHARS = 24
_RC = {"axes.spines.top": False, "axes.spines.right": False, "legend.fontsize": 8}


@dataclass(frozen=True)
class Thresholds:
    """参考线的位置，为 None 的不画；入选以表中的 ``selected`` 为准。"""
    p: float | None = None            # 纵轴为 ``metric`` 时的水平线
    log2fc: float | None = None       # 竖线 ±log2fc
    effect_size: float | None = None  # 纵轴为 |效应量| 时的水平线


def contrast_label(contrast: Contrast) -> str:
    """如 ``GLP1 vs Saline (diet=HFD)``。"""
    label = f"{contrast.numerator} vs {contrast.denominator}"
    if contrast.subset:
        label += " (" + ", ".join(f"{k}={v}" for k, v in contrast.subset.items()) + ")"
    return label


def _numeric(panel: pd.DataFrame, column: str) -> pd.Series:
    if column not in panel.columns:
        return pd.Series(np.nan, index=panel.index)
    return pd.to_numeric(panel[column], errors="coerce").astype(float)


def _y_axis(panel: pd.DataFrame, metric: str) -> tuple[pd.Series, str, str] | None:
    """纵轴 (取值, 所用列, 轴标题)：先用 ``metric``，没有任何值时依次退回 p 值与 |效应量|。"""
    for column in dict.fromkeys((metric, "pvalue")):
        values = _numeric(panel, column)
        if values.notna().any():
            return -np.log10(values.clip(lower=_P_FLOOR)), column, _Y_LABELS[column]
    values = _numeric(panel, "effect_size")
    if values.notna().any():
        kinds = panel["effect_size_type"].dropna().astype(str) if "effect_size_type" in panel else pd.Series()
        kind = kinds.iloc[0] if len(kinds) else ""
        return values.abs(), "effect_size", f"|{_EFFECT_NAMES.get(kind, kind or 'effect size')}|"
    return None


def _span(values: pd.Series) -> str:
    """生物学重复数：各特征相同时写一个数，不同时写范围（缺失值较多的蛋白组常见）。"""
    counts = pd.to_numeric(values, errors="coerce").dropna().astype(int)
    if counts.empty:
        return "?"
    low, high = counts.min(), counts.max()
    return str(low) if low == high else f"{low}–{high}"


def _short(name: str) -> str:
    return name if len(name) <= _LABEL_CHARS else name[:_LABEL_CHARS - 1] + "…"


def _label_points(ax: Axes, x: np.ndarray, y: np.ndarray, names: list[str]) -> None:
    try:
        from adjustText import adjust_text
    except ImportError:
        adjust_text = None
    if adjust_text is None:
        for xi, yi, name in zip(x, y, names):
            ax.annotate(name, (xi, yi), xytext=(3, 2), textcoords="offset points", fontsize=7, zorder=4)
        return
    texts = [ax.text(xi, yi, name, fontsize=7, zorder=4) for xi, yi, name in zip(x, y, names)]
    adjust_text(texts, ax=ax, arrowprops={"arrowstyle": "-", "color": "0.5", "lw": 0.5})


def plot_volcano(panel: pd.DataFrame, *, title: str = "", metric: str = "padj",
                 thresholds: Thresholds = Thresholds(), x_cap: float | None = None, n_labels: int = 10,
                 colors: dict[str, str] | None = None) -> Figure | None:
    """一个组织 × 比较的火山图；没有可画的点时返回 None。

    Parameters
    ----------
    panel : DataFrame
        统一差异结果表中该组织、该比较的行。
    metric : ``"padj"`` 或 ``"pvalue"``
        纵轴首选的指标，应与该组学入选规则所用的指标相同。
    thresholds : Thresholds
        参考线；``thresholds.p`` 只在纵轴确为 ``metric`` 时画出。
    x_cap : float, optional
        横轴半宽的上限；超出的点钉在边界上。MFA 的净通量在分母接近 0 时 |log2FC| 可达 20 以上，属于数值假象。
    n_labels : int
        标注的入选特征个数。

    应在 :func:`omicsweave.visualization.theme.themed` 中调用并保存（:func:`write_volcanoes` 已这样做）。
    """
    if metric not in _Y_LABELS:
        raise ValueError(f"metric 应为 {' / '.join(_Y_LABELS)}，而不是 {metric!r}")
    axis = _y_axis(panel, metric)
    if axis is None:
        return None
    y_all, y_column, y_label = axis
    x_all = _numeric(panel, "log2fc")
    shown = (x_all.notna() & np.isfinite(y_all)).to_numpy()
    if not shown.any():
        return None
    rows = panel.loc[shown]
    x, y = x_all.to_numpy()[shown], y_all.to_numpy()[shown]
    selected = (rows["selected"].eq(True).fillna(False).astype(bool).to_numpy() if "selected" in rows
                else np.zeros(len(rows), dtype=bool))
    category = np.where(selected & (x > 0), "up", np.where(selected & (x < 0), "down", "ns"))

    # 横轴对称于 0；超出半宽（或为 ±inf）的点钉在边界上
    finite = np.isfinite(x)
    half = float(np.abs(x[finite]).max()) if finite.any() else 1.0
    if x_cap is not None:
        half = min(half, x_cap)
    if thresholds.log2fc:
        half = max(half, thresholds.log2fc * 1.5)
    half = half or 1.0
    off = ~finite | (np.abs(x) > half)
    x_plot = np.where(off, np.sign(x) * half * 1.04, x)

    fig = Figure(figsize=(6, 5))
    FigureCanvasAgg(fig)
    ax = fig.add_subplot()
    colors = {**COLORS, **(colors or {})}
    size = 5 if len(x) > 2000 else 12 if len(x) > 300 else 22
    rasterized = len(x) > 500  # 点多时矢量图中只把散点栅格化，文字与坐标轴仍为矢量
    for cat in ("ns", "down", "up"):  # 未入选的点在最下层
        inside, outside = (category == cat) & ~off, (category == cat) & off
        if inside.any():
            ax.scatter(x_plot[inside], y[inside], s=size, c=colors[cat], alpha=0.5 if cat == "ns" else 0.85,
                       linewidths=0, rasterized=rasterized, zorder=1 if cat == "ns" else 2)
        for side, marker in ((1, ">"), (-1, "<")):
            pinned = outside & (np.sign(x) == side)
            if pinned.any():
                ax.scatter(x_plot[pinned], y[pinned], s=size * 2.5, marker=marker, facecolors="none",
                           edgecolors=colors[cat], linewidths=0.8, zorder=3)

    line = {"color": "0.45", "lw": 0.8, "ls": "--", "zorder": 0}
    hline = None
    if y_column == metric and thresholds.p:
        hline = -np.log10(thresholds.p)
    elif y_column == "effect_size" and thresholds.effect_size:
        hline = thresholds.effect_size
    if hline is not None:
        ax.axhline(hline, **line)
    if thresholds.log2fc:
        for v in (-thresholds.log2fc, thresholds.log2fc):
            ax.axvline(v, **line)

    if n_labels > 0 and selected.any():
        names = rows["feature_name"].fillna("").astype(str) if "feature_name" in rows else pd.Series("", index=rows.index)
        names = names.where(names.str.strip() != "", rows["feature_id"].astype(str))
        frame = pd.DataFrame({"x": x_plot, "y": y, "abs_x": np.abs(x_plot), "name": names.to_numpy()})[selected]
        top = frame.sort_values(["y", "abs_x"], ascending=False, kind="mergesort").head(n_labels)
        _label_points(ax, top["x"].to_numpy(), top["y"].to_numpy(), [_short(n) for n in top["name"]])

    y_top = max(float(y.max()), hline or 0.0)
    ax.set_xlim(-half * 1.08, half * 1.08)
    ax.set_ylim(0, y_top * 1.06 if y_top > 0 else 1.0)
    ax.set_xlabel("log$_2$ fold change")
    ax.set_ylabel(y_label)
    if title:
        ax.set_title(title, loc="left")

    handles = [Line2D([], [], marker="o", ls="", color=colors[c], markersize=5,
                      label=f"{_LEGEND[c]} ({int((category == c).sum()):,})") for c in ("up", "down", "ns")]
    if off.any():
        handles.append(Line2D([], [], marker=">", ls="", markerfacecolor="none", markeredgecolor="0.4",
                              markersize=5, label=f"Beyond x-axis ({int(off.sum()):,})"))
    # 图例放在坐标区右侧，不遮挡点和标签
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1), frameon=False, borderaxespad=0)

    notes = [f"n = {_span(panel.get('n_numerator', pd.Series()))} vs "
             f"{_span(panel.get('n_denominator', pd.Series()))}"]
    rules = sorted(set(panel["selection_rule"].dropna().astype(str)) - {""}) if "selection_rule" in panel else []
    if rules:
        notes.append("Selected: " + "; ".join(rules))
    hidden = int((~shown).sum())
    if hidden:
        notes.append(f"Not shown: {hidden:,} features without log2FC or {_Y_NAMES[y_column]}")
    floor = _numeric(panel, "min_achievable_p").min()
    if y_column != "effect_size" and np.isfinite(floor):
        notes.append(f"Smallest p the rank test can reach at this sample size: {floor:.3g}")
    ax.annotate("\n".join(textwrap.fill(n, 100) for n in notes), xy=(0, 0), xycoords="axes fraction",
                xytext=(0, -36), textcoords="offset points", va="top", ha="left", fontsize=7, color="0.35")
    return fig


def write_volcanoes(table: pd.DataFrame, directory: Path, *, contrasts: Iterable[Contrast] = (),
                    metric: str = "padj", thresholds: Thresholds = Thresholds(), x_cap: float | None = None,
                    n_labels: int = 10, formats: Sequence[str] = FIGURE_FORMATS) -> tuple[list[Path], list[str]]:
    """为统一差异结果表中的每个组织 × 比较画一张火山图，写到 ``directory``；返回 (写出的文件, 提示)。

    ``contrasts`` 给出图题中的比较说明（如 ``HFD vs LFD``）与图的顺序；表中其余比较按比较 ID 排在后面。
    """
    labels = {c.id: contrast_label(c) for c in contrasts}
    order = {cid: i for i, cid in enumerate(labels)}
    panels = dict(tuple(table.groupby(["contrast_id", "tissue"], sort=False)))
    keys = sorted(panels, key=lambda k: (order.get(k[0], len(order)), str(k[0]), str(k[1])))
    paths: list[Path] = []
    skipped: list[str] = []
    for contrast_id, tissue in keys:
        with themed(rc=_RC):
            fig = plot_volcano(panels[(contrast_id, tissue)], title=f"{tissue} · {labels.get(contrast_id, contrast_id)}",
                               metric=metric, thresholds=thresholds, x_cap=x_cap, n_labels=n_labels)
            if fig is None:
                skipped.append(f"{contrast_id} / {tissue}")
                continue
            paths += save_figure(fig, directory, figure_stem(KIND, tissue, contrast_id), formats)
    notes = [f"以下比较 / 组织没有可画的 p 值或效应量，未画火山图：{'、'.join(skipped)}"] if skipped else []
    return paths, notes


def write_bundle_volcanoes(omics_dir: Path, contrasts: Iterable[Contrast] = (), **options) -> list[str]:
    """结果包中一个组学的火山图：读取 ``tables/differential.parquet``，写到 ``figures/``；返回提示。

    ``options`` 同 :func:`write_volcanoes`。各组学运行器的 ``figures`` 步骤调用它。
    """
    path = omics_dir / "tables" / "differential.parquet"
    if not path.is_file():
        return [f"{omics_dir.name}：没有差异结果表，未画火山图"]
    _, notes = write_volcanoes(read_differential(path), omics_dir / "figures", contrasts=contrasts, **options)
    return [f"{omics_dir.name}：{note}" for note in notes]
