#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
期刊出版规格（``JournalProfile``）与自动化合规审计引擎（零业务耦合，可跨项目复用）。

将「视觉风格（``metro_atlas`` / ``economist`` / ``nature_editorial``）」与
「期刊印刷排版规范（``cell_press`` / ``nature``）」正交解耦：
- 默认出版规格为 ``cell_press``（Cell Press Figure Guidelines）：
  * 标准栏宽：单栏 **85 mm**、1.5 栏 **114 mm**、双栏 **174 mm**，常规页高上限 **200 mm**（整页 **225 mm**）；
  * 字体体系：统一无衬线 **Arial / Helvetica**，数学符号强制无衬线 ``stixsans``，PDF/EPS 嵌入 Type 42 字体；
  * 字号底线：最终印刷尺寸下正文图任意可读文本 **≥ 6.0 pt**（主体 7–8 pt，面板字母 **A, B, C** 10 pt 粗体大写）；
  * 线宽底线：主要边框/轴线 **≥ 0.50 pt**，辅助细线不低于 **0.25 pt**；
  * 正文净版排版：出版子图默认不在图内嵌入标题栏与运行哈希图注，配套输出 ``<stem>_caption.txt``；
  * Graphical Abstract：默认严格 **1200 × 1200 px（300 dpi，4.0 × 4.0 in）** 1:1 正方形单面板，**Arial 12–16 pt**。
- 保留 ``nature``（89 / 120 / 183 mm，页高 247 mm，小写加粗面板标号 ``a, b, c``）随时切换。
"""
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.text import Text

DEFAULT_JOURNAL_PROFILE = "cell_press"

# 无衬线字体识别关键词（用于合规审计）
_SANS_TOKENS: tuple[str, ...] = (
    "sans",
    "arial",
    "helvetica",
    "liberation",
    "dejavu sans",
    "noto sans",
    "yahei",
    "simhei",
    "pingfang",
    "heiti",
    "source han sans",
    "wenquanyi",
    "inter",
    "roboto",
    "calibri",
    "segoe",
    "avenir",
    "myriad",
)


@dataclass(frozen=True)
class JournalProfile:
    """一套完整、不可变的期刊印刷排版与格式契约。"""

    name: str
    display_name: str
    # 物理栏宽与页高（mm）
    single_col_mm: float
    one_half_col_mm: float
    double_col_mm: float
    max_page_h_mm: float
    full_page_h_mm: float
    # 字体与排版（pt）
    font_families: tuple[str, ...]
    mathtext_fontset: str
    min_font_pt: float
    body_font_pt: float
    tick_font_pt: float
    legend_font_pt: float
    title_font_pt: float
    panel_label_pt: float
    panel_label_uppercase: bool
    # 线宽（pt）
    min_linewidth_pt: float
    hairline_lw_pt: float
    axis_lw_pt: float
    # 色彩与输出
    color_space: str
    default_dpi: int
    line_art_dpi: int
    embed_fonts: bool
    clean_publication_framing: bool
    # Graphical Abstract 契约
    ga_layout: str
    ga_width_px: int
    ga_height_px: int
    ga_dpi: int
    ga_min_font_pt: float
    ga_max_font_pt: float

    @property
    def single_col_in(self) -> float:
        return self.single_col_mm / 25.4

    @property
    def one_half_col_in(self) -> float:
        return self.one_half_col_mm / 25.4

    @property
    def double_col_in(self) -> float:
        return self.double_col_mm / 25.4

    @property
    def max_page_h_in(self) -> float:
        return self.max_page_h_mm / 25.4

    @property
    def full_page_h_in(self) -> float:
        return self.full_page_h_mm / 25.4

    @property
    def ga_size_in(self) -> tuple[float, float]:
        return (self.ga_width_px / float(self.ga_dpi), self.ga_height_px / float(self.ga_dpi))

    def width_mm(self, columns: float = 2.0) -> float:
        """按栏数（1 / 1.5 / 2 或自定义比例）返回物理宽度（mm）。"""
        if columns == 1 or columns == 1.0:
            return self.single_col_mm
        if columns == 1.5:
            return self.one_half_col_mm
        if columns == 2 or columns == 2.0:
            return self.double_col_mm
        return self.double_col_mm * (float(columns) / 2.0)

    def width_in(self, columns: float = 2.0) -> float:
        return self.width_mm(columns) / 25.4

    def clamp_font(self, pt: float, *, floor: float | None = None) -> float:
        """确保最终印刷字号不低于当前期刊规范底线。"""
        min_pt = self.min_font_pt if floor is None else float(floor)
        return max(float(pt), min_pt)

    def clamp_lw(self, lw: float, *, floor: float | None = None) -> float:
        """确保线条宽度不低于当前期刊最小印刷线宽底线。"""
        min_lw = self.min_linewidth_pt if floor is None else float(floor)
        return max(float(lw), min_lw)

    def format_panel_label(self, label: str) -> str:
        """按期刊规范格式化面板角标（Cell Press：大写无句点 ``A``；Nature：小写 ``a``）。"""
        cleaned = str(label).strip().strip("().: ")
        return cleaned.upper() if self.panel_label_uppercase else cleaned.lower()

    def rc_overrides(self) -> dict[str, Any]:
        """返回与该期刊规范对应的 matplotlib ``rcParams`` 增量。"""
        rc: dict[str, Any] = {
            "font.size": self.body_font_pt,
            "axes.titlesize": self.title_font_pt,
            "axes.labelsize": self.body_font_pt,
            "xtick.labelsize": self.tick_font_pt,
            "ytick.labelsize": self.tick_font_pt,
            "legend.fontsize": self.legend_font_pt,
            "axes.linewidth": self.axis_lw_pt,
            "xtick.major.width": max(self.hairline_lw_pt, self.axis_lw_pt * 0.92),
            "ytick.major.width": max(self.hairline_lw_pt, self.axis_lw_pt * 0.92),
            "grid.linewidth": self.hairline_lw_pt,
            "mathtext.fontset": self.mathtext_fontset,
            "mathtext.default": "regular",
        }
        if self.embed_fonts:
            rc.update({
                "pdf.fonttype": 42,
                "ps.fonttype": 42,
                "svg.fonttype": "none",
            })
        return rc


CELL_PRESS = JournalProfile(
    name="cell_press",
    display_name="Cell Press",
    single_col_mm=85.0,
    one_half_col_mm=114.0,
    double_col_mm=174.0,
    max_page_h_mm=200.0,
    full_page_h_mm=225.0,
    font_families=("Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"),
    mathtext_fontset="stixsans",
    min_font_pt=6.0,
    body_font_pt=8.0,
    tick_font_pt=7.0,
    legend_font_pt=7.0,
    title_font_pt=9.0,
    panel_label_pt=10.0,
    panel_label_uppercase=True,
    min_linewidth_pt=0.25,
    hairline_lw_pt=0.50,
    axis_lw_pt=0.65,
    color_space="RGB",
    default_dpi=300,
    line_art_dpi=600,
    embed_fonts=True,
    clean_publication_framing=True,
    ga_layout="cell_press_square",
    ga_width_px=1200,
    ga_height_px=1200,
    ga_dpi=300,
    ga_min_font_pt=12.0,
    ga_max_font_pt=16.0,
)

NATURE = JournalProfile(
    name="nature",
    display_name="Nature Portfolio",
    single_col_mm=89.0,
    one_half_col_mm=120.0,
    double_col_mm=183.0,
    max_page_h_mm=247.0,
    full_page_h_mm=247.0,
    font_families=("Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"),
    mathtext_fontset="stixsans",
    min_font_pt=5.5,
    body_font_pt=8.0,
    tick_font_pt=7.0,
    legend_font_pt=7.0,
    title_font_pt=9.0,
    panel_label_pt=8.0,
    panel_label_uppercase=False,
    min_linewidth_pt=0.25,
    hairline_lw_pt=0.45,
    axis_lw_pt=0.60,
    color_space="RGB",
    default_dpi=300,
    line_art_dpi=600,
    embed_fonts=True,
    clean_publication_framing=False,
    ga_layout="widescreen",
    ga_width_px=2550,
    ga_height_px=1700,
    ga_dpi=300,
    ga_min_font_pt=6.0,
    ga_max_font_pt=12.0,
)

JOURNAL_PROFILES: dict[str, JournalProfile] = {
    "cell_press": CELL_PRESS,
    "nature": NATURE,
}


def list_journal_profiles() -> tuple[str, ...]:
    """返回全部可用期刊出版规格名称。"""
    return tuple(JOURNAL_PROFILES)


def get_journal_profile(profile: str | JournalProfile = DEFAULT_JOURNAL_PROFILE) -> JournalProfile:
    """按名称（或直接传入 :class:`JournalProfile` 实例）获取期刊出版规格。"""
    if isinstance(profile, JournalProfile):
        return profile
    key = str(profile).strip().lower()
    if key not in JOURNAL_PROFILES:
        raise ValueError(f"未知的期刊出版规格 {profile!r}（可用：{', '.join(JOURNAL_PROFILES)}）")
    return JOURNAL_PROFILES[key]


def add_panel_label(
    target: Axes | Figure,
    label: str,
    *,
    profile: str | JournalProfile = DEFAULT_JOURNAL_PROFILE,
    x: float | None = None,
    y: float | None = None,
    fontsize: float | None = None,
    color: str = "#111827",
) -> Text:
    """按期刊规范在面板左上角添加标准面板字母（Cell Press：粗体大写 ``A``；Nature：粗体小写 ``a``）。"""
    spec = get_journal_profile(profile)
    text = spec.format_panel_label(label)
    fs = float(fontsize) if fontsize is not None else spec.panel_label_pt
    if isinstance(target, Axes):
        px = -0.10 if x is None else float(x)
        py = 1.04 if y is None else float(y)
        return target.text(
            px, py, text,
            transform=target.transAxes,
            fontsize=fs, fontweight="bold", color=color,
            ha="left", va="bottom", clip_on=False, zorder=10.0,
        )
    px = 0.008 if x is None else float(x)
    py = 0.992 if y is None else float(y)
    return target.text(
        px, py, text,
        fontsize=fs, fontweight="bold", color=color,
        ha="left", va="top", zorder=10.0,
    )


def write_caption_sidecar(
    stem: Path | str,
    title: str,
    caption: Sequence[str] = (),
    *,
    panel_note: str = "",
    method_note: str = "",
) -> Path:
    """在净版出版模式下将完整标题、统计方法与溯源元数据写入同名 ``<stem>_caption.txt``。"""
    base = Path(stem)
    out_path = base.with_name(f"{base.name}_caption.txt")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    if title:
        lines.append(str(title).strip())
    if panel_note:
        lines.append(str(panel_note).strip())
    if method_note:
        lines.append(str(method_note).strip())
    for item in caption:
        text = str(item).strip()
        if text:
            lines.append(text)
    out_path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
    return out_path


# ---------------------------------------------------------------------------
# 自动化出版合规审计（字号、字体族、物理栏宽、线宽底线与方形 GA 核验）
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ComplianceIssue:
    severity: str   # "error" | "warning"
    rule: str
    message: str


@dataclass(frozen=True)
class ComplianceReport:
    profile: str
    passed: bool
    width_mm: float
    height_mm: float
    min_font_pt: float | None
    max_font_pt: float | None
    min_linewidth_pt: float | None
    issues: tuple[ComplianceIssue, ...]


def _is_sans_family(families: Sequence[str] | str) -> bool:
    if isinstance(families, str):
        names = [families]
    else:
        names = [str(f) for f in families]
    for name in names:
        low = name.lower()
        if any(tok in low for tok in _SANS_TOKENS):
            return True
    return False


def audit_figure_compliance(
    fig: Figure,
    *,
    profile: str | JournalProfile = DEFAULT_JOURNAL_PROFILE,
    check_ga: bool = False,
    check_width: bool = True,
    font_tol_pt: float = 0.05,
) -> ComplianceReport:
    """对 Matplotlib :class:`Figure` 执行期刊出版合规预检。

    检查项：
    1. 物理幅面：标准图件宽度不超过双栏上限（Cell Press 174 mm）；GA 核验精确 1:1 正方形（1200×1200 px @ 300 dpi）；
    2. 最小字号：遍历图内全部可见 ``Text``，核验 ``fontsize >= min_font_pt``（GA 为 ``>= ga_min_font_pt``）；
    3. 无衬线字体：核验全部可见文本均属于无衬线字体族（Arial / Helvetica 等）；
    4. 最小线宽：遍历 ``Line2D`` 与带描边 ``Patch``，核验 ``linewidth >= min_linewidth_pt``。
    """
    spec = get_journal_profile(profile)
    w_in, h_in = (float(fig.get_figwidth()), float(fig.get_figheight()))
    w_mm, h_mm = w_in * 25.4, h_in * 25.4
    issues: list[ComplianceIssue] = []

    if check_ga and spec.ga_layout == "cell_press_square":
        exp_w_in, exp_h_in = spec.ga_size_in
        if abs(w_in - exp_w_in) > 0.02 or abs(h_in - exp_h_in) > 0.02:
            issues.append(ComplianceIssue(
                "error", "ga_dimensions",
                f"GA 尺寸 {w_in:.2f}×{h_in:.2f} in 与 {spec.display_name} 要求 {exp_w_in:.2f}×{exp_h_in:.2f} in "
                f"（{spec.ga_width_px}×{spec.ga_height_px} px @ {spec.ga_dpi} dpi）不一致",
            ))
    elif check_width:
        if w_mm > spec.double_col_mm + 1.0:
            issues.append(ComplianceIssue(
                "error", "column_width",
                f"图宽 {w_mm:.1f} mm 超出 {spec.display_name} 双栏上限 {spec.double_col_mm:.1f} mm",
            ))
        if h_mm > spec.full_page_h_mm + 1.0:
            issues.append(ComplianceIssue(
                "warning", "page_height",
                f"图高 {h_mm:.1f} mm 超出 {spec.display_name} 整页高度上限 {spec.full_page_h_mm:.1f} mm",
            ))

    min_allowed_font = spec.ga_min_font_pt if check_ga else spec.min_font_pt
    font_sizes: list[float] = []
    seen_non_sans = False
    small_samples: list[str] = []

    texts: list[Text] = list(fig.texts)
    for ax in fig.axes:
        texts.extend(ax.texts)
        texts.extend([ax.title, ax.xaxis.label, ax.yaxis.label])
        if ax.axison:
            texts.extend(ax.get_xticklabels())
            texts.extend(ax.get_yticklabels())
        leg = ax.get_legend()
        if leg is not None:
            texts.extend(leg.get_texts())
            if leg.get_title() is not None:
                texts.append(leg.get_title())
    for leg in fig.legends:
        texts.extend(leg.get_texts())

    for t in texts:
        if not t.get_visible() or not str(t.get_text()).strip():
            continue
        fs = float(t.get_fontsize())
        font_sizes.append(fs)
        if fs + font_tol_pt < min_allowed_font:
            if len(small_samples) < 5:
                small_samples.append(f"{t.get_text()!r} ({fs:.2f} pt)")
        if not _is_sans_family(t.get_fontfamily()):
            seen_non_sans = True

    if small_samples:
        issues.append(ComplianceIssue(
            "error", "min_font_size",
            f"检测到低于 {spec.display_name} 下限 {min_allowed_font:.1f} pt 的文字：{', '.join(small_samples)}",
        ))
    if seen_non_sans:
        issues.append(ComplianceIssue(
            "error", "font_family",
            f"检测到非无衬线字体，{spec.display_name} 要求统一使用 Arial / Helvetica 等无衬线体",
        ))

    lws: list[float] = []
    thin_count = 0
    for ax in fig.axes:
        for line in ax.lines:
            if not line.get_visible():
                continue
            lw = float(line.get_linewidth() or 0.0)
            if lw > 0:
                lws.append(lw)
                if lw + 1e-4 < spec.min_linewidth_pt:
                    thin_count += 1
        for patch in ax.patches:
            if not patch.get_visible():
                continue
            ec = patch.get_edgecolor()
            lw = float(patch.get_linewidth() or 0.0)
            has_edge = lw > 0 and not (isinstance(ec, tuple) and len(ec) == 4 and ec[3] == 0)
            if has_edge:
                lws.append(lw)
                if lw + 1e-4 < spec.min_linewidth_pt:
                    thin_count += 1
    if thin_count:
        issues.append(ComplianceIssue(
            "warning", "min_linewidth",
            f"存在 {thin_count} 处描边线宽小于 {spec.min_linewidth_pt:.2f} pt（最小 {min(lws):.2f} pt）",
        ))

    passed = not any(item.severity == "error" for item in issues)
    return ComplianceReport(
        profile=spec.name,
        passed=passed,
        width_mm=round(w_mm, 2),
        height_mm=round(h_mm, 2),
        min_font_pt=round(min(font_sizes), 2) if font_sizes else None,
        max_font_pt=round(max(font_sizes), 2) if font_sizes else None,
        min_linewidth_pt=round(min(lws), 3) if lws else None,
        issues=tuple(issues),
    )


__all__ = [
    "CELL_PRESS",
    "ComplianceIssue",
    "ComplianceReport",
    "DEFAULT_JOURNAL_PROFILE",
    "JOURNAL_PROFILES",
    "JournalProfile",
    "NATURE",
    "add_panel_label",
    "audit_figure_compliance",
    "get_journal_profile",
    "list_journal_profiles",
    "write_caption_sidecar",
]
