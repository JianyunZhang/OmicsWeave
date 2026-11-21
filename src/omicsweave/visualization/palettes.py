#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
统一语义色板与发散色阶（色盲安全，按出版预设组织）。

设计原则（见 ``docs/arch_design/omicsweave_visualization_and_graphical_abstract.md`` 第 3.2 节）：
- 高饱和红 / 蓝唯一保留给生物学方向（上调 / 下调）；
- 证据等级（L3 / L2 / L1）采用无彩色黑白灰阶材质，彻底消除与涨跌方向的色相冲突；
- 分类色板取自色盲友好的 Okabe-Ito / ColorBrewer Dark2 调和体系；
- 不依赖任何业务模块或 seaborn，供所有组学模块与外部项目直接复用。
"""
from __future__ import annotations

from dataclasses import dataclass

import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap, Normalize


@dataclass(frozen=True)
class SemanticPalette:
    """一套完整的出版级语义配色令牌。"""

    name: str
    # 画布与排版阶层
    surface: str                 # 整图背景
    card_bg: str                 # 面板 / 侧栏 / 图例卡片底色
    card_edge: str               # 卡片描边
    ink: str                     # 主文字与显著描边
    ink2: str                    # 次文字
    muted: str                   # 注释、未检出、弱化标签
    hairline: str                # 细分隔线
    # 方向与连续效应量（唯一高饱和红 / 蓝通道）
    up: str
    down: str
    neutral: str
    diverging_stops: tuple[str, ...]
    # 对比分组（图 B / D 等分子 vs 分母）
    numerator: str
    denominator: str
    # 骨架与缺测层级
    edge_grey: str               # 已拟合 MFA 但无指向的边
    omics_grey: str              # 仅组学虚线边
    faint: str                   # 未评估骨架
    missing_fill: str            # 缺测微格底色
    missing_hatch: str           # 缺测微格细斜纹色
    # 耦合通量组（与红 / 蓝涨跌区分的靛青通道）
    coupled_track: str
    coupled_fill: str
    coupled_edge: str
    coupled_ink: str
    # 变号标记
    sign_change_bg: str
    sign_change_edge: str
    sign_change_ink: str
    # 细胞器区室
    mito_face: str
    mito_edge: str
    mito_inner: str
    mito_ink: str
    perox_face: str
    perox_edge: str
    perox_inner: str
    perox_ink: str
    # 边界与库存节点
    boundary_face: str
    boundary_edge: str
    boundary_ink: str
    store_face: str
    store_edge: str
    store_ink: str
    # 分类色板（色盲友好）
    categorical: tuple[str, ...]

    def diverging_cmap(self) -> LinearSegmentedColormap:
        """固定连续发散色阶（蓝降 → 极浅中灰白 → 红升）。"""
        return LinearSegmentedColormap.from_list(f"ow_{self.name}_diverging", list(self.diverging_stops))


# ── 三套核心出版 / 展示配色 ──────────────────────────────────────────────

NATURE_EDITORIAL = SemanticPalette(
    name="nature_editorial",
    surface="#FFFFFF",
    card_bg="#F8F9FA",
    card_edge="#D8DEE6",
    ink="#111827",
    ink2="#374151",
    muted="#6B7280",
    hairline="#E5E7EB",
    up="#C83E37",
    down="#2166AC",
    neutral="#9CA3AF",
    diverging_stops=("#1B4F8A", "#3B82C4", "#A8C9EA", "#F7F7F7", "#F4B6A8", "#D9483B", "#8F1D24"),
    numerator="#D95F02",
    denominator="#2166AC",
    edge_grey="#6B7280",
    omics_grey="#9CA3AF",
    faint="#D1D5DB",
    missing_fill="#F3F4F6",
    missing_hatch="#D1D5DB",
    coupled_track="#E0E7FF",
    coupled_fill="#EEF2FF",
    coupled_edge="#6366F1",
    coupled_ink="#3730A3",
    sign_change_bg="#FFFBEB",
    sign_change_edge="#B45309",
    sign_change_ink="#92400E",
    mito_face="#FAF8F3",
    mito_edge="#A8A29E",
    mito_inner="#D6D3D1",
    mito_ink="#57534E",
    perox_face="#F0FDFA",
    perox_edge="#5EEAD4",
    perox_inner="#99F6E4",
    perox_ink="#0F766E",
    boundary_face="#F1F5F9",
    boundary_edge="#94A3B8",
    boundary_ink="#334155",
    store_face="#F5F3EF",
    store_edge="#A8A29E",
    store_ink="#44403C",
    categorical=("#2166AC", "#D95F02", "#1B9E77", "#7570B3", "#E7298A", "#66A61E", "#A6761D", "#666666"),
)

METRO_ATLAS = SemanticPalette(
    name="metro_atlas",
    surface="#FFFFFF",
    card_bg="#F8FAFC",
    card_edge="#CBD5E1",
    ink="#0F172A",
    ink2="#334155",
    muted="#64748B",
    hairline="#E2E8F0",
    up="#C83E37",
    down="#1D5FA7",
    neutral="#94A3B8",
    diverging_stops=("#17467B", "#3174B8", "#9EC3E8", "#F8FAFC", "#F4B4A6", "#D44236", "#8B1A22"),
    numerator="#D95F02",
    denominator="#1D5FA7",
    edge_grey="#64748B",
    omics_grey="#94A3B8",
    faint="#CBD5E1",
    missing_fill="#F1F5F9",
    missing_hatch="#CBD5E1",
    coupled_track="#DBEAFE",
    coupled_fill="#EFF6FF",
    coupled_edge="#4F46E5",
    coupled_ink="#312E81",
    sign_change_bg="#FEF3C7",
    sign_change_edge="#B45309",
    sign_change_ink="#78350F",
    mito_face="#F8F6F0",
    mito_edge="#94A3B8",
    mito_inner="#CBD5E1",
    mito_ink="#475569",
    perox_face="#ECFDF5",
    perox_edge="#34D399",
    perox_inner="#A7F3D0",
    perox_ink="#065F46",
    boundary_face="#EFF6FF",
    boundary_edge="#60A5FA",
    boundary_ink="#1E3A8A",
    store_face="#F1F5F9",
    store_edge="#94A3B8",
    store_ink="#334155",
    categorical=("#1D5FA7", "#D95F02", "#0D9488", "#6D28D9", "#DB2777", "#65A30D", "#B45309", "#475569"),
)

SLIDE_DARK = SemanticPalette(
    name="slide_dark",
    surface="#111827",
    card_bg="#1F2937",
    card_edge="#374151",
    ink="#F9FAFB",
    ink2="#D1D5DB",
    muted="#9CA3AF",
    hairline="#374151",
    up="#FB7185",
    down="#38BDF8",
    neutral="#6B7280",
    diverging_stops=("#0284C7", "#38BDF8", "#7DD3FC", "#1F2937", "#FDA4AF", "#FB7185", "#E11D48"),
    numerator="#FB923C",
    denominator="#38BDF8",
    edge_grey="#6B7280",
    omics_grey="#4B5563",
    faint="#374151",
    missing_fill="#1F2937",
    missing_hatch="#374151",
    coupled_track="#312E81",
    coupled_fill="#1E1B4B",
    coupled_edge="#818CF8",
    coupled_ink="#C7D2FE",
    sign_change_bg="#451A03",
    sign_change_edge="#F59E0B",
    sign_change_ink="#FDE68A",
    mito_face="#1C1917",
    mito_edge="#57534E",
    mito_inner="#44403C",
    mito_ink="#D6D3D1",
    perox_face="#042F2E",
    perox_edge="#14B8A6",
    perox_inner="#115E59",
    perox_ink="#5EEAD4",
    boundary_face="#1E293B",
    boundary_edge="#475569",
    boundary_ink="#CBD5E1",
    store_face="#262626",
    store_edge="#525252",
    store_ink="#D4D4D4",
    categorical=("#38BDF8", "#FB923C", "#34D399", "#A78BFA", "#F472B6", "#A3E635", "#FBBF24", "#9CA3AF"),
)

PALETTES: dict[str, SemanticPalette] = {
    "nature_editorial": NATURE_EDITORIAL,
    "metro_atlas": METRO_ATLAS,
    "slide_dark": SLIDE_DARK,
    "poster_editorial": NATURE_EDITORIAL,
}

# 默认全局常量（对应 nature_editorial，供所有模块统一导入）
DIRECTIONAL: dict[str, str] = {
    "up": NATURE_EDITORIAL.up,
    "down": NATURE_EDITORIAL.down,
    "ns": "#BDBDBD",
    "neutral": NATURE_EDITORIAL.neutral,
}

# 证据等级无彩色材质令牌（彻底剥离色相冲突）
TIER_TOKENS: dict[str, dict[str, object]] = {
    "L3": {
        "dots": 3,
        "bg": "#111827",
        "ec": "#111827",
        "fg": "#FFFFFF",
        "dot_on": "#FFFFFF",
        "dot_off": "#4B5563",
        "lw": 0.55,
        "ls": "-",
        "label": "Confirmed Flux",
    },
    "L2": {
        "dots": 2,
        "bg": "#374151",
        "ec": "#374151",
        "fg": "#FFFFFF",
        "dot_on": "#FFFFFF",
        "dot_off": "#6B7280",
        "lw": 0.55,
        "ls": "-",
        "label": "Flux-Supported",
    },
    "L1": {
        "dots": 1,
        "bg": "#FFFFFF",
        "ec": "#64748B",
        "fg": "#334155",
        "dot_on": "#334155",
        "dot_off": "#CBD5E1",
        "lw": 0.55,
        "ls": (0, (1.8, 1.0)),
        "label": "Enzyme-Predicted",
    },
}

# 暗色背景下的证据等级反转令牌
TIER_TOKENS_DARK: dict[str, dict[str, object]] = {
    "L3": {
        "dots": 3,
        "bg": "#F9FAFB",
        "ec": "#F9FAFB",
        "fg": "#111827",
        "dot_on": "#111827",
        "dot_off": "#9CA3AF",
        "lw": 0.55,
        "ls": "-",
        "label": "Confirmed Flux",
    },
    "L2": {
        "dots": 2,
        "bg": "#D1D5DB",
        "ec": "#D1D5DB",
        "fg": "#111827",
        "dot_on": "#111827",
        "dot_off": "#6B7280",
        "lw": 0.55,
        "ls": "-",
        "label": "Flux-Supported",
    },
    "L1": {
        "dots": 1,
        "bg": "#1F2937",
        "ec": "#9CA3AF",
        "fg": "#E5E7EB",
        "dot_on": "#E5E7EB",
        "dot_off": "#4B5563",
        "lw": 0.55,
        "ls": (0, (1.8, 1.0)),
        "label": "Enzyme-Predicted",
    },
}

# 色盲安全的多组学证据等级色板（替换旧版 scoring_plots 的红绿配色）
EVIDENCE_GRADE_COLORS: dict[str, str] = {
    "Replicated": "#0F5C6E",
    "Convergent": "#2A8C96",
    "Single-layer": "#7EC4C8",
    "Exploratory": "#E6A266",
    "Conflicting": "#C83E37",
    "Omics-concordant": "#2166AC",
    "Omics-single-layer": "#6BAED6",
    "Omics-discordant": "#D95F02",
    "No-signal": "#D1D5DB",
    "Not-assessable": "#ECEFF3",
}

# 示踪剂可辨识性单色序列（清晰区分 4 档）
IDENTIFIABILITY_COLORS: dict[str, str] = {
    "well_constrained": "#1B4F8A",
    "moderate": "#4A89C7",
    "poorly_constrained": "#B9D4EE",
    "bound_pinned": "#E5E7EB",
}

# 地铁模式（metro_atlas）颜色语义开关：
# - "semantic"（默认）：与正文 nature_editorial 统一，红/蓝/灰唯一编码生物学涨跌方向；
# - "pathway"（可选）：低饱和莫兰迪学术通路分区色，适用于无实验对比数据的概念综述或导览图。
VALID_COLOR_MODES: tuple[str, ...] = ("semantic", "pathway")

METRO_PATHWAY_COLORS: dict[str, str] = {
    "carbohydrate": "#C87A28",   # 琥珀暖棕（Carbohydrate / Glycolysis / GNG / Glycogen）
    "respiration": "#B84A42",    # 砖红（Pyruvate / TCA / Cellular Respiration）
    "pentose": "#945294",        # 绛紫（Pentose Phosphate）
    "lipid": "#4A8A5C",          # 苔绿（Fatty Acid / Glycerolipid / Peroxisome）
    "ketone": "#B56275",         # 烟粉褐（Ketone Bodies / β-Oxidation branch）
    "amino": "#387CA0",          # 湖蓝（Amino Acids / BCAA / Urea）
    "nucleotide": "#7B8548",     # 橄榄绿（Nucleotide / Cofactor / Other）
    "redox": "#5E6E82",          # 石板灰蓝（OXPHOS / Redox）
    "circulation": "#4679A3",    # 静脉蓝（Blood / Circulation / Boundary）
    "default": "#5E6E82",
}

_PATHWAY_FAMILY_MAP: dict[str, str] = {
    # Figure E region IDs & GA hub / icon IDs
    "glycolysis_gng": "carbohydrate",
    "glycogen": "carbohydrate",
    "glucose": "carbohydrate",
    "carbohydrate": "carbohydrate",
    "ppp": "pentose",
    "pentose": "pentose",
    "pyruvate_tca": "respiration",
    "pyruvate_hub": "respiration",
    "pyruvate": "respiration",
    "tca_cycle": "respiration",
    "tca": "respiration",
    "respiration": "respiration",
    "ketone": "ketone",
    "fatty_acid": "lipid",
    "glycerolipid": "lipid",
    "peroxisome": "lipid",
    "lipid": "lipid",
    "amino_acid": "amino",
    "urea_creatine": "amino",
    "amino": "amino",
    "nucleotide": "nucleotide",
    "other": "nucleotide",
    "oxphos": "redox",
    "atp": "redox",
    "redox": "redox",
    "blood_glc": "circulation",
    "blood": "circulation",
    "circulation": "circulation",
}


def pathway_family_color(key: str, *, fallback_accent: str = "") -> str:
    """将通路区域 ID、概念枢纽 ID 或图标类型映射为克制的学术莫兰迪地铁分区色。"""
    norm = str(key or "").strip().lower()
    if norm in _PATHWAY_FAMILY_MAP:
        return METRO_PATHWAY_COLORS[_PATHWAY_FAMILY_MAP[norm]]
    for token, fam in _PATHWAY_FAMILY_MAP.items():
        if token in norm:
            return METRO_PATHWAY_COLORS[fam]
    if fallback_accent:
        return fallback_accent
    return METRO_PATHWAY_COLORS["default"]

DEFAULT_VMAX = 2.0
DEFAULT_NORM = Normalize(-DEFAULT_VMAX, DEFAULT_VMAX, clip=True)


def get_palette(name: str = "nature_editorial") -> SemanticPalette:
    """按预设名返回对应的语义配色对象；未知名称抛出 ValueError。"""
    if name not in PALETTES:
        raise ValueError(f"未知的配色预设 {name!r}（可用：{', '.join(PALETTES)}）")
    return PALETTES[name]


def tier_tokens(preset: str = "nature_editorial") -> dict[str, dict[str, object]]:
    """按背景明暗返回无彩色证据等级令牌。"""
    return TIER_TOKENS_DARK if preset == "slide_dark" else TIER_TOKENS


def tint_color(hex_color: str, amount: float = 0.08, *, base: str = "#FFFFFF") -> tuple[float, float, float]:
    """将主题色与底色按比例混合，生成柔和低饱和的分区底色或边框色。"""
    if not hex_color:
        return mcolors.to_rgb(base)
    r, g, b = mcolors.to_rgb(hex_color)
    br, bg, bb = mcolors.to_rgb(base)
    return (
        br * (1.0 - amount) + r * amount,
        bg * (1.0 - amount) + g * amount,
        bb * (1.0 - amount) + b * amount,
    )
