#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
正文分专题代谢子通路预设目录与信号感知筛选器（零 contract / pydantic / seaborn 依赖）。

设计定位：
- 与全景图 E（``map_<比较>``）共享同一份 ``pathway_map`` 视图表与同一套 ``network_layout.yaml`` 坐标，
  零数据冗余、数值 100% 一致；
- 面向正文单栏 / 双栏机制图：将 13 个区域按经典代谢专题组织为 5 张默认核心子图 + 2 张可选扩展子图，
  让每张正文子图聚焦 13–52 条反应，字号自然达到 6.5–8 pt 印刷可读标准；
- 内置信号感知跳过机制（``min_signal``）：避免在缺少对应通路数据的组织/比较中生成大量全灰空壳图。
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import pandas as pd

from ..mappings.network_layout import Layout

VALID_MIN_SIGNALS: tuple[str, ...] = ("assessed", "active", "all")


@dataclass(frozen=True)
class SubmapSpec:
    """单个正文代谢专题子图的声明式规格。"""

    id: str
    title: str
    subtitle: str
    regions: tuple[str, ...]
    core: bool = True
    include_reactions: tuple[str, ...] = ()
    exclude_reactions: tuple[str, ...] = ()


# 5 张默认核心子图 + 2 张可选扩展子图（按生物学叙事顺序列出）
SUBMAP_CATALOG: tuple[SubmapSpec, ...] = (
    SubmapSpec(
        id="central_carbon",
        title="Carbohydrate & Central Carbon Metabolism",
        subtitle="Glycogen · Glycolysis / GNG · Pentose Phosphate · Pyruvate Hub & TCA",
        regions=("glycogen", "glycolysis_gng", "ppp", "pyruvate_tca"),
        core=True,
    ),
    SubmapSpec(
        id="fatty_acid_ketone",
        title="Fatty Acid, Ketone & Peroxisomal Metabolism",
        subtitle="Long-Chain FA · Lipid Mobilization · Ketogenesis · Peroxisomal β-Ox",
        regions=("fatty_acid", "ketone", "peroxisome"),
        core=True,
        # lipid_mobilization_c 在全局画布中物理坐标位于右侧脂滴/FFA 区（x=151），自然归入本子图
        include_reactions=("lipid_mobilization_c",),
    ),
    SubmapSpec(
        id="glycerolipid",
        title="Glycerol & Glycerolipid Metabolism",
        subtitle="G3P · LPA · PA · DAG · TAG · Phospholipid Synthesis & Lipolysis",
        regions=("glycerolipid",),
        core=True,
        # 剔除位于远端右侧脂滴区的 lipid_mobilization_c（已归入 fatty_acid_ketone），保持紧凑单栏幅面
        exclude_reactions=("lipid_mobilization_c",),
    ),
    SubmapSpec(
        id="amino_acid",
        title="Amino Acid & BCAA Metabolism",
        subtitle="Glutamate / Aspartate / Serine · Transamination & BCAA Anaplerosis",
        regions=("amino_acid",),
        core=True,
    ),
    SubmapSpec(
        id="urea_creatine",
        title="Urea Cycle, Polyamines & Creatine",
        subtitle="Arginine / Ornithine / Citrulline · Polyamine & Creatine Synthesis",
        regions=("urea_creatine",),
        core=True,
    ),
    SubmapSpec(
        id="nucleotide",
        title="Purine & Pyrimidine Nucleotide Metabolism",
        subtitle="De Novo & Salvage Purine / Pyrimidine Pathways",
        regions=("nucleotide",),
        core=False,
    ),
    SubmapSpec(
        id="redox_one_carbon",
        title="One-Carbon, Glutathione & Redox Sub-Networks",
        subtitle="Glutathione · Folate 1C · Methionine / Cysteine · NAD(P) & Tryptophan Modules",
        regions=("other",),
        core=False,
    ),
)

SUBMAP_BY_ID: dict[str, SubmapSpec] = {spec.id: spec for spec in SUBMAP_CATALOG}
CORE_SUBMAP_IDS: tuple[str, ...] = tuple(spec.id for spec in SUBMAP_CATALOG if spec.core)
ALL_SUBMAP_IDS: tuple[str, ...] = tuple(spec.id for spec in SUBMAP_CATALOG)


def get_submap(panel_id: str) -> SubmapSpec:
    """按子图 ID 获取预设规格；未知 ID 抛出 ValueError。"""
    if panel_id not in SUBMAP_BY_ID:
        raise ValueError(f"未知的代谢子图 ID {panel_id!r}（可用：{', '.join(ALL_SUBMAP_IDS)}）")
    return SUBMAP_BY_ID[panel_id]


def list_submaps(*, core_only: bool = False) -> tuple[SubmapSpec, ...]:
    """列出全部（或仅 5 张核心）代谢子图预设。"""
    return tuple(s for s in SUBMAP_CATALOG if (s.core or not core_only))


def submap_signal_stats(view: pd.DataFrame, spec: SubmapSpec) -> dict[str, int | bool]:
    """统计指定子图在当前 ``view`` 中的覆盖与活跃信号数。"""
    is_el = view["element_type"].isin(("reaction", "aggregate"))
    in_reg = view["region"].isin(spec.regions)
    if spec.include_reactions:
        in_reg = in_reg | (view["element_type"].eq("reaction") & view["element_id"].isin(spec.include_reactions))
    if spec.exclude_reactions:
        in_reg = in_reg & ~(view["element_type"].eq("reaction") & view["element_id"].isin(spec.exclude_reactions))
    sub = view[is_el & in_reg]
    if sub.empty:
        return {
            "n_total": 0, "n_assessed": 0, "n_flux": 0, "n_l2": 0, "n_l1": 0,
            "n_directional": 0, "has_assessed": False, "has_active": False,
        }
    draw_style = sub["draw_style"].fillna("").astype(str)
    claims = sub["flux_claim"].fillna("").astype(str)
    f_state = sub["F_state"].fillna("").astype(str) if "F_state" in sub.columns else pd.Series("", index=sub.index)
    c_state = sub["C_state"].fillna("").astype(str) if "C_state" in sub.columns else pd.Series("", index=sub.index)

    n_assessed = int((draw_style != "not_assessed").sum())
    n_flux = int((draw_style == "flux").sum())
    n_l2 = int(claims.isin(("L2", "L3")).sum())
    n_l1 = int((claims == "L1").sum())
    n_directional = int((f_state.isin(("changed", "suggestive")) | c_state.isin(("changed", "suggestive"))).sum())
    return {
        "n_total": int(len(sub)),
        "n_assessed": n_assessed,
        "n_flux": n_flux,
        "n_l2": n_l2,
        "n_l1": n_l1,
        "n_directional": n_directional,
        "has_assessed": n_assessed > 0,
        "has_active": bool(n_l2 > 0 or n_l1 > 0 or n_directional > 0),
    }


def select_submaps_for_view(
    view: pd.DataFrame,
    layout: Layout,
    *,
    panels: Sequence[str] | None = None,
    min_signal: str = "assessed",
    include_optional_active: bool = True,
) -> tuple[list[SubmapSpec], list[tuple[str, str]]]:
    """根据配置与当前视图的真实信号筛选应渲染的正文子图。

    Returns
    -------
    (selected, skipped)
        ``selected`` 为按顺序待渲染的 :class:`SubmapSpec` 列表；
        ``skipped`` 为 ``(submap_id, reason)`` 跳过说明，便于写入运行日志。
    """
    if min_signal not in VALID_MIN_SIGNALS:
        raise ValueError(f"submaps.min_signal 应为 {' / '.join(VALID_MIN_SIGNALS)} 之一，实际为 {min_signal!r}")

    explicit = list(panels or ())
    if explicit:
        if any(str(p).strip().lower() == "all" for p in explicit):
            candidates = list(SUBMAP_CATALOG)
        else:
            candidates = [get_submap(str(p).strip()) for p in dict.fromkeys(explicit)]
        check_optional = False
    else:
        candidates = [s for s in SUBMAP_CATALOG if s.core]
        check_optional = bool(include_optional_active)

    selected: list[SubmapSpec] = []
    skipped: list[tuple[str, str]] = []

    def _present_regions(spec: SubmapSpec) -> tuple[str, ...]:
        return tuple(r for r in spec.regions if r in layout.regions)

    for spec in candidates:
        valid_regs = _present_regions(spec)
        if not valid_regs:
            skipped.append((spec.id, "layout_missing_regions"))
            continue
        effective = SubmapSpec(
            spec.id, spec.title, spec.subtitle, valid_regs, spec.core,
            include_reactions=spec.include_reactions, exclude_reactions=spec.exclude_reactions,
        )
        stats = submap_signal_stats(view, effective)
        if min_signal == "all":
            selected.append(effective)
        elif min_signal == "active":
            if stats["has_active"]:
                selected.append(effective)
            else:
                skipped.append((spec.id, "no_active_signal"))
        else:  # "assessed"
            if stats["has_assessed"]:
                selected.append(effective)
            else:
                skipped.append((spec.id, "no_assessed_elements"))

    if check_optional:
        for spec in SUBMAP_CATALOG:
            if spec.core:
                continue
            valid_regs = _present_regions(spec)
            if not valid_regs:
                continue
            effective = SubmapSpec(
            spec.id, spec.title, spec.subtitle, valid_regs, spec.core,
            include_reactions=spec.include_reactions, exclude_reactions=spec.exclude_reactions,
        )
            stats = submap_signal_stats(view, effective)
            if stats["has_active"]:
                selected.append(effective)
            else:
                skipped.append((spec.id, "optional_inactive"))

    return selected, skipped


__all__ = [
    "ALL_SUBMAP_IDS",
    "CORE_SUBMAP_IDS",
    "SUBMAP_BY_ID",
    "SUBMAP_CATALOG",
    "SubmapSpec",
    "VALID_MIN_SIGNALS",
    "get_submap",
    "list_submaps",
    "select_submaps_for_view",
    "submap_signal_stats",
]
