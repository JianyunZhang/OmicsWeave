#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Graphical Abstract（概念摘要图）的数据契约、护栏校验与数据驱动聚合器。

核心约束：
- 严格限定为**小规模高层概念摘要（``MAX_HUBS = 12``）**，禁止将 100+ 反应的明细网络塞进手绘或像素风；
- 支持两种输入路径：
  1. :func:`build_abstract_from_pathway_view`：直接从图 E 视图表聚合出 6–10 个核心代谢枢纽（所有方向、
     证据等级与注记均严格溯源至真实评估结果）；
  2. :func:`load_abstract_spec` / :func:`spec_from_dict`：任意课题可编写轻量 YAML 跨项目复用。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from ..glyphs import glyph_ids

MAX_HUBS = 12
VALID_ICONS = glyph_ids(category="concept")
VALID_COMPARTMENTS = ("cytosol", "mitochondria", "peroxisome", "blood", "extracellular")
TIER_RANK = {"L3": 3, "L2": 2, "L1": 1, "omics": 0, "": -1}


@dataclass(frozen=True)
class AbstractCompartment:
    """高层区室容器（如线粒体基质、细胞外血浆带）。"""

    id: str
    label: str
    box: tuple[float, float, float, float]   # (x0, y0, x1, y1)，y 轴向下
    kind: str = "mitochondria"               # mitochondria | peroxisome | blood | cytosol


@dataclass(frozen=True)
class AbstractHub:
    """概念摘要图中的单个核心通路/功能枢纽。"""

    id: str
    label: str
    xy: tuple[float, float]
    direction: int = 0                       # +1 上调 / -1 下调 / 0 无明显方向或混合
    magnitude: float = 0.0                   # 0.0–2.0 归一化强度
    evidence_tier: str = ""                  # "L3" | "L2" | "L1" | "omics" | ""
    subtitle: str = ""                       # 极简数据溯源副标（如 "Flux ↑ · 3 L2"）
    compartment: str = "cytosol"
    icon: str = "generic"
    n_l2: int = 0
    n_l1: int = 0


@dataclass(frozen=True)
class AbstractFlow:
    """概念摘要图中两个枢纽之间的宏观代谢/调控流向。"""

    source: str
    target: str
    label: str = ""
    direction: int = 0
    magnitude: float = 0.0
    evidence_tier: str = ""
    curve: float = 0.0
    coupled: bool = False


@dataclass(frozen=True)
class AbstractSpec:
    """完整的 Graphical Abstract 声明式规格。"""

    title: str
    hubs: tuple[AbstractHub, ...]
    flows: tuple[AbstractFlow, ...] = ()
    compartments: tuple[AbstractCompartment, ...] = ()
    subtitle: str = ""
    contrast: str = ""
    takeaway: str = ""
    rules_note: str = ""
    provenance: str = ""
    canvas: tuple[float, float] = (120.0, 74.0)
    metadata: dict[str, Any] = field(default_factory=dict)


DEFAULT_AGGREGATION_RULES: str = (
    "Quantification: Hub tier = max claim (Confirmed > Flux > Enz-pred); "
    "Dir. priority: Flux (F sig./sug. majority) > Enzyme (C majority) > Pool (|median log2FC| ≥ 0.45); "
    "hits = L2+ count (else L1); Track = higher-tier end."
)


def validate_spec(spec: AbstractSpec) -> None:
    """校验摘要图规格：强制 ``1 <= len(hubs) <= MAX_HUBS`` 且边端点均已定义。"""
    if not spec.hubs:
        raise ValueError("Graphical Abstract 至少需要 1 个枢纽节点（hubs 不能为空）")
    if len(spec.hubs) > MAX_HUBS:
        raise ValueError(
            f"Graphical Abstract 仅用于小规模概念摘要（最多 {MAX_HUBS} 个枢纽，当前 {len(spec.hubs)} 个）；"
            "大尺度全量网络请使用标准通路总图（mode='atlas' 或 'publication'）"
        )
    ids = [h.id for h in spec.hubs]
    if len(set(ids)) != len(ids):
        raise ValueError("Graphical Abstract 中存在重复的 hub id")
    id_set = set(ids)
    for hub in spec.hubs:
        if hub.icon not in VALID_ICONS:
            raise ValueError(f"hub {hub.id!r} 的 icon {hub.icon!r} 不支持（可用：{', '.join(VALID_ICONS)}）")
    for flow in spec.flows:
        if flow.source not in id_set or flow.target not in id_set:
            raise ValueError(f"flow {flow.source!r} -> {flow.target!r} 引用了未定义的 hub")


def spec_from_dict(doc: dict[str, Any]) -> AbstractSpec:
    """从字典构建并校验 :class:`AbstractSpec`。"""
    raw_canvas = doc.get("canvas") or [120.0, 74.0]
    if isinstance(raw_canvas, dict):
        canvas = (float(raw_canvas.get("width", 120.0)), float(raw_canvas.get("height", 74.0)))
    else:
        canvas = (float(raw_canvas[0]), float(raw_canvas[1]))

    compartments = tuple(
        AbstractCompartment(
            id=str(c["id"]),
            label=str(c.get("label", c["id"])),
            box=tuple(float(v) for v in c["box"]),
            kind=str(c.get("kind", "mitochondria")),
        )
        for c in (doc.get("compartments") or ())
    )
    hubs = tuple(
        AbstractHub(
            id=str(h["id"]),
            label=str(h.get("label", h["id"])),
            xy=(float(h["xy"][0]), float(h["xy"][1])),
            direction=int(h.get("direction", 0)),
            magnitude=float(h.get("magnitude", 0.0)),
            evidence_tier=str(h.get("evidence_tier", "")),
            subtitle=str(h.get("subtitle", "")),
            compartment=str(h.get("compartment", "cytosol")),
            icon=str(h.get("icon", "generic")),
            n_l2=int(h.get("n_l2", 0)),
            n_l1=int(h.get("n_l1", 0)),
        )
        for h in (doc.get("hubs") or ())
    )
    flows = tuple(
        AbstractFlow(
            source=str(f["source"]),
            target=str(f["target"]),
            label=str(f.get("label", "")),
            direction=int(f.get("direction", 0)),
            magnitude=float(f.get("magnitude", 0.0)),
            evidence_tier=str(f.get("evidence_tier", "")),
            curve=float(f.get("curve", 0.0)),
            coupled=bool(f.get("coupled", False)),
        )
        for f in (doc.get("flows") or ())
    )
    spec = AbstractSpec(
        title=str(doc.get("title", "Metabolic Rewiring Summary")),
        subtitle=str(doc.get("subtitle", "")),
        contrast=str(doc.get("contrast", "")),
        takeaway=str(doc.get("takeaway", "")),
        rules_note=str(doc.get("rules_note", "")),
        provenance=str(doc.get("provenance", "")),
        hubs=hubs,
        flows=flows,
        compartments=compartments,
        canvas=canvas,
        metadata=dict(doc.get("metadata") or {}),
    )
    validate_spec(spec)
    return spec


def load_abstract_spec(path: str | Path) -> AbstractSpec:
    """从 YAML 文件读取 Graphical Abstract 规格。"""
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"找不到 Graphical Abstract 规格文件：{p}")
    doc = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    return spec_from_dict(doc)


# ---------------------------------------------------------------------------
# 从图 E 视图表自动聚合为 8–10 枢纽的高层代谢摘要图
# ---------------------------------------------------------------------------

# 经典中心代谢概念布局（画布 138 x 80，充足间距保证边注记与枢纽标题零截断）
_CANONICAL_HUBS: tuple[dict[str, Any], ...] = (
    {
        "id": "blood_glc",
        "label": "Circulating Fuel",
        "regions": (),
        "metabolites": ("GLC_c", "LAC_c", "BHB_c"),
        "xy": (15.0, 17.5),
        "compartment": "blood",
        "icon": "blood",
    },
    {
        "id": "glycolysis_gng",
        "label": "Glycolysis / GNG",
        "regions": ("glycolysis_gng", "glycogen"),
        "metabolites": ("GLC6P_c", "FRU6P_c", "GAP_c", "PEP_c"),
        "xy": (47.0, 17.5),
        "compartment": "cytosol",
        "icon": "glucose",
    },
    {
        "id": "ppp",
        "label": "Pentose Phosphate",
        "regions": ("ppp",),
        "metabolites": ("RIB5P_c", "SED7P_c", "E4P_c"),
        "xy": (81.5, 13.5),
        "compartment": "cytosol",
        "icon": "ppp",
    },
    {
        "id": "glycerolipid",
        "label": "G3P & Glycerolipids",
        "regions": ("glycerolipid",),
        "metabolites": ("G3P_c", "TAG_c", "DAG_c"),
        "xy": (16.0, 45.0),
        "compartment": "cytosol",
        "icon": "glycerolipid",
    },
    {
        "id": "pyruvate_hub",
        "label": "Pyruvate Hub",
        "regions": ("pyruvate_tca",),
        "rxn_filter": ("PYR_trans", "PDH_m", "PC_m", "LDH_c", "PEPCK_c", "ME_m", "ME2_c"),
        "metabolites": ("PYR_c", "PYR_m"),
        "xy": (47.0, 45.0),
        "compartment": "mitochondria",
        "icon": "pyruvate",
    },
    {
        "id": "tca_cycle",
        "label": "TCA Cycle",
        "regions": ("pyruvate_tca",),
        "rxn_filter": ("CS_m", "ACONT_m", "ICDH_m", "AKGD_m", "SUCOAS_m", "SUCD_m", "FUMH_m", "MDH_m"),
        "metabolites": ("CIT_m", "ICIT_m", "AKG_m", "SUC_m", "FUM_m", "MAL_m"),
        "xy": (83.0, 47.0),
        "compartment": "mitochondria",
        "icon": "tca",
    },
    {
        "id": "ketone",
        "label": "Ketone Bodies",
        "regions": ("ketone",),
        "metabolites": ("BHB_m", "ACAC_m"),
        "xy": (118.5, 45.0),
        "compartment": "mitochondria",
        "icon": "ketone",
    },
    {
        "id": "fatty_acid",
        "label": "Fatty Acid & β-Ox",
        "regions": ("fatty_acid", "peroxisome"),
        "metabolites": ("PAL_c", "TGPAL_c", "CHOL_c", "FA_c"),
        "xy": (118.5, 17.5),
        "compartment": "cytosol",
        "icon": "lipid",
    },
    {
        "id": "amino_acid",
        "label": "Amino Acids & BCAA",
        "regions": ("amino_acid",),
        "metabolites": ("GLU_c", "GLN_c", "ASP_c", "ALA_c", "SER_c"),
        "xy": (47.0, 67.5),
        "compartment": "cytosol",
        "icon": "amino",
    },
    {
        "id": "oxphos",
        "label": "OXPHOS & Redox",
        "regions": ("oxphos", "other"),
        "metabolites": ("GSH_c", "NAD_POOL_c"),
        "xy": (118.5, 67.5),
        "compartment": "mitochondria",
        "icon": "atp",
    },
)

_CANONICAL_FLOWS: tuple[tuple[str, str, str, float], ...] = (
    ("blood_glc", "glycolysis_gng", "Glc uptake", 0.0),
    ("glycolysis_gng", "ppp", "G6P shunt", -2.2),
    ("glycolysis_gng", "glycerolipid", "DHAP → G3P", 2.5),
    ("glycolysis_gng", "pyruvate_hub", "Pyruvate", 0.0),
    ("pyruvate_hub", "tca_cycle", "Ac-CoA / OAA", 0.0),
    ("tca_cycle", "ketone", "Ac-CoA", 0.0),
    ("fatty_acid", "tca_cycle", "β-Ox / Citrate", 3.5),
    ("fatty_acid", "ketone", "Acyl-CoA", 0.0),
    ("amino_acid", "tca_cycle", "Anaplerosis", -2.5),
    ("tca_cycle", "oxphos", "NADH / FADH2", 2.0),
)


def _summarize_hub_rows(
    view: pd.DataFrame,
    spec_def: dict[str, Any],
) -> tuple[int, float, str, str, int, int, bool]:
    """从图 E 视图行汇总单个概念枢纽的 (direction, magnitude, tier, subtitle, n_l2, n_l1, has_data)。"""
    rxns = view[view["element_type"].isin(("reaction", "aggregate"))]
    mets = view[view["element_type"] == "metabolite"]

    regions = spec_def["regions"]
    rxn_filter = spec_def.get("rxn_filter")
    if rxn_filter:
        sub_r = rxns[rxns["element_id"].isin(rxn_filter)]
    elif regions:
        sub_r = rxns[rxns["region"].isin(regions)]
    else:
        sub_r = rxns.iloc[0:0]

    sub_m = mets[mets["element_id"].isin(spec_def.get("metabolites", ()))]

    # 统计 L2+/L1 等级
    claims = sub_r["flux_claim"].fillna("").astype(str) if not sub_r.empty else pd.Series(dtype=str)
    n_l3 = int((claims == "L3").sum())
    n_l2 = int(claims.isin(("L2", "L3")).sum())
    n_l1 = int((claims == "L1").sum())
    n_flux = int((sub_r["draw_style"] == "flux").sum()) if not sub_r.empty else 0
    n_omics = int((sub_r["draw_style"].isin(("omics_only", "aggregate"))).sum()) if not sub_r.empty else 0
    m_vals = pd.to_numeric(sub_m["value"], errors="coerce").dropna() if not sub_m.empty else pd.Series(dtype=float)

    has_data = bool(n_flux > 0 or n_omics > 0 or not m_vals.empty)
    if not has_data:
        return 0, 0.0, "", "No data", 0, 0, False

    tier = "L3" if n_l3 > 0 else ("L2" if n_l2 > 0 else ("L1" if n_l1 > 0 else ("omics" if (n_omics or not m_vals.empty) else "")))

    # 1) 优先取 F 轴显著/提示方向（最硬通量证据）
    f_sig = sub_r[sub_r["F_state"].isin(("changed", "suggestive"))] if not sub_r.empty else sub_r
    c_sig = sub_r[sub_r["C_state"].isin(("changed", "suggestive"))] if not sub_r.empty else sub_r

    direction = 0
    magnitude = 0.0
    parts: list[str] = []

    if not f_sig.empty:
        dirs = pd.to_numeric(f_sig["F_direction"], errors="coerce").fillna(0)
        vals = pd.to_numeric(f_sig["value"], errors="coerce").abs().dropna()
        pos, neg = int((dirs > 0).sum()), int((dirs < 0).sum())
        direction = 1 if pos > neg else (-1 if neg > pos else 0)
        magnitude = float(np.clip(vals.median() if not vals.empty else 1.2, 0.8, 2.0))
        arrow = "↑" if direction > 0 else ("↓" if direction < 0 else "↕")
        parts.append(f"Flux {arrow}")
    elif not c_sig.empty:
        dirs = pd.to_numeric(c_sig["C_direction"], errors="coerce").fillna(0)
        pos, neg = int((dirs > 0).sum()), int((dirs < 0).sum())
        direction = 1 if pos > neg else (-1 if neg > pos else 0)
        magnitude = 1.05 if direction else 0.4
        arrow = "↑" if direction > 0 else ("↓" if direction < 0 else "↕")
        parts.append(f"Enzyme {arrow}")
    elif not m_vals.empty:
        med = float(m_vals.median())
        if abs(med) >= 0.45:
            direction = 1 if med > 0 else -1
            magnitude = float(np.clip(abs(med), 0.5, 2.0))
            parts.append(f"Pool {'↑' if direction > 0 else '↓'}")
        else:
            parts.append("Pool stable")
    else:
        parts.append("Measured · stable")

    if n_l2 > 0:
        parts.append(f"{n_l2} {'hit' if n_l2 == 1 else 'hits'}")
    elif n_l1 > 0:
        parts.append(f"{n_l1} {'hit' if n_l1 == 1 else 'hits'}")

    return direction, magnitude, tier, " · ".join(parts), n_l2, n_l1, True


def build_abstract_from_pathway_view(
    view: pd.DataFrame,
    regions: pd.DataFrame | None = None,
    *,
    title: str = "",
    tissue: str = "",
    contrast: str = "",
    provenance: str = "",
    max_hubs: int = 10,
) -> AbstractSpec:
    """将图 E 的全量视图表自动聚合成 ≤ 10 枢纽的 Graphical Abstract 规格。

    所有枢纽的方向、强度与证据层级均严格来自真实评估表，绝不虚构未测信号。
    """
    cap = min(max(int(max_hubs), 3), MAX_HUBS)
    hubs: list[AbstractHub] = []
    hub_map: dict[str, AbstractHub] = {}

    for hdef in _CANONICAL_HUBS:
        direction, magnitude, tier, subtitle, n_l2, n_l1, has_data = _summarize_hub_rows(view, hdef)
        if not has_data and hdef["id"] != "blood_glc":
            continue
        hub = AbstractHub(
            id=hdef["id"],
            label=hdef["label"],
            xy=hdef["xy"],
            direction=direction,
            magnitude=magnitude,
            evidence_tier=tier,
            subtitle=subtitle,
            compartment=hdef["compartment"],
            icon=hdef["icon"],
            n_l2=n_l2,
            n_l1=n_l1,
        )
        hubs.append(hub)
        hub_map[hub.id] = hub
        if len(hubs) >= cap:
            break

    # 若出现极端空白表，保底返回已汇总的枢纽
    if not hubs:
        hubs = [
            AbstractHub(id="central_metabolism", label="Central Metabolism", xy=(60.0, 38.0), subtitle="No assessed rows")
        ]
        hub_map = {h.id: h for h in hubs}

    # 检查 TCA 耦合通量组标记
    tca_rows = view[(view["element_type"] == "reaction") & (view["region"] == "pyruvate_tca")
                    & view["F_state"].isin(("changed", "suggestive"))]
    tca_coupled = bool(not tca_rows.empty and tca_rows["flux_group"].dropna().astype(str).str.strip().ne("").any())

    flows: list[AbstractFlow] = []
    for src, dst, label, curve in _CANONICAL_FLOWS:
        if src not in hub_map or dst not in hub_map:
            continue
        sh, dh = hub_map[src], hub_map[dst]
        # 流向优先取证据更强的一端（但绝不无中生有）
        primary = dh if TIER_RANK.get(dh.evidence_tier, -1) >= TIER_RANK.get(sh.evidence_tier, -1) else sh
        flows.append(AbstractFlow(
            source=src,
            target=dst,
            label=label,
            direction=primary.direction if primary.evidence_tier in ("L3", "L2", "L1") else 0,
            magnitude=primary.magnitude if primary.evidence_tier in ("L3", "L2", "L1") else 0.0,
            evidence_tier=primary.evidence_tier,
            curve=curve,
            coupled=(src == "pyruvate_hub" and dst == "tca_cycle" and tca_coupled),
        ))

    compartments: list[AbstractCompartment] = [
        AbstractCompartment("blood", "Circulation", (2.5, 7.0, 27.5, 28.0), kind="blood"),
        AbstractCompartment("mito", "Mitochondrial Matrix", (33.5, 32.5, 133.5, 77.5), kind="mitochondria"),
    ]

    # 生成诚实的 1 行 Takeaway 摘要
    rxn_view = view[view["element_type"].isin(("reaction", "aggregate"))]
    total_l2 = int(rxn_view["flux_claim"].isin(("L2", "L3")).sum()) if not rxn_view.empty else 0
    total_l1 = int((rxn_view["flux_claim"] == "L1").sum()) if not rxn_view.empty else 0
    l2_hubs = [h.label for h in hubs if h.n_l2 > 0]
    if l2_hubs:
        takeaway = (
            f"Summary: {total_l2} Flux-Supported (L2+) hits in {', '.join(l2_hubs[:3])}; "
            f"{total_l1} Enzyme-Predicted (L1) candidates (verify L1 against decoy baseline)."
        )
    else:
        takeaway = (
            f"Summary: 0 Flux-Supported (L2+) hits; {total_l1} Enzyme-Predicted (L1) candidates "
            "(omics-only inference — verify against decoy baseline)."
        )

    tissue_nice = str(tissue).replace("_", " ").title() if tissue else "Tissue"
    contrast_nice = str(contrast).replace("_vs_", " vs ") if contrast else ""
    header_title = title or f"{tissue_nice} Metabolic Rewiring — Graphical Abstract"
    sub = f"{contrast_nice}  ·  High-level Multi-Omics & MFA Conceptual Summary (≤{MAX_HUBS} Hubs)" if contrast_nice else (
        f"High-level Multi-Omics & MFA Conceptual Summary (≤{MAX_HUBS} Hubs)"
    )

    spec = AbstractSpec(
        title=header_title,
        subtitle=sub,
        contrast=contrast_nice,
        takeaway=takeaway,
        rules_note=DEFAULT_AGGREGATION_RULES,
        provenance=provenance,
        hubs=tuple(hubs),
        flows=tuple(flows),
        compartments=tuple(compartments),
        canvas=(138.0, 80.0),
        metadata={"tissue": tissue, "contrast": contrast, "n_l2": total_l2, "n_l1": total_l1},
    )
    validate_spec(spec)
    return spec


__all__ = [
    "AbstractCompartment",
    "AbstractFlow",
    "AbstractHub",
    "AbstractSpec",
    "MAX_HUBS",
    "VALID_COMPARTMENTS",
    "VALID_ICONS",
    "build_abstract_from_pathway_view",
    "load_abstract_spec",
    "spec_from_dict",
    "validate_spec",
]
