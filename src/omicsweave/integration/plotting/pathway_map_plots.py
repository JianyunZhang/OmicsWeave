#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure E (G18): Tissue Multi-Omics Pathway Atlas — renders strictly from
:func:`omicsweave.integration.dossier.map_view.map_view`.

Refactored onto the decoupled :mod:`omicsweave.visualization` foundation:
- Default style: ``metro_atlas`` (and ``nature_editorial`` / ``slide_dark`` / ``editorial_ink`` presets; pure white
  surface, colorblind-safe red/blue reserved strictly for biological direction, monochrome 3-dot tier ladders
  ``●●● Confirmed`` / ``●●○ Flux`` / ``●○○ Enz-pred``, isolated top-left ``#n`` rank tabs, translucent indigo
  coupled-flux conduit tracks + ``Linked flux`` badges, and 4-omics pure-vector icon strips ``Flux | Protein | mRNA | Metab``).
- Dual-mode output:
  * ``mode="atlas"``: full 13-region diagnostic / supplemental atlas with right dashboard sidebar and
    4-column visual legend.
  * ``mode="publication"``: focused journal-ready sub-map cropping automatically to active regions
    (or user-specified ``focus_regions``), folding inactive grey regions into a compact summary strip,
    hiding exploratory rank tabs and non-key micro-strips, and boosting typography for double-column print.
"""
from __future__ import annotations

import math
import re
import textwrap
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyBboxPatch, Polygon, Rectangle

from omicsweave.visualization.figures import DOUBLE_COL_IN, ONE_HALF_COL_IN, save_stem
from omicsweave.visualization.metro import (
    draw_metro_coupled_halo,
    draw_metro_station,
    draw_metro_track,
    draw_metro_zone,
    fillet_polyline,
    polyline_length,
    route_octilinear,
    sample_polyline_at,
    sample_polyline_points,
    trim_polyline,
)
from omicsweave.visualization.palettes import (
    DEFAULT_NORM,
    DEFAULT_VMAX,
    VALID_COLOR_MODES,
    SemanticPalette,
    get_palette,
    pathway_family_color,
    tier_tokens,
    tint_color,
)
from omicsweave.visualization.primitives import (
    AXIS_LAYER_ICONS,
    OMICS_LAYER_ICONS,
    coupled_badge_width,
    draw_axis_strip,
    draw_coupled_badge,
    draw_coupled_track,
    draw_layer_icon,
    draw_organelle_box,
    draw_quad_arrow,
    draw_rank_tab,
    draw_sign_reversal_badge,
    draw_status_cell,
    draw_tier_badge,
    draw_zone_card,
    draw_zone_header,
    em_units,
    estimate_text_width,
    quad_point,
    sample_quad,
    tier_badge_width,
)
from omicsweave.visualization.theme import DEFAULT_PRESET, themed

from ..mappings.network_annotation import pathway_reactions
from ..mappings.network_layout import Layout

UNIT = 0.108                     # inches per canvas unit in atlas mode (~22.5 x 15.6 in)
SIDEBAR = 44.0                   # right dashboard sidebar width (atlas mode)
HEADER = 6.2                     # top header bar height
FOOTER = 19.6                    # bottom 4-column legend + provenance caption height
NODE_R = 1.05                    # metabolite node radius
CELL = 1.04                      # single cell width/height in the 4-layer omics micro-strip
AXIS_HDR_H = 0.62                # icon header height on highlighted cards
AXIS_BAR_H = 0.14                # slim neutral header bar on standard cards
NOTE_SIZE = 3.6                  # secondary substrate/product note font size
MIN_EDGE = 0.8                   # floor |log2FC| for coloring directional F-axis edges

VALID_MODES = ("atlas", "publication")
VALID_STRIP_MODES = ("omics", "axes")
STRIP_AXIS_COLS = (
    ("F_state", "F_direction"),
    ("C_state", "C_direction"),
    ("R_state", "R_direction"),
    ("M_state", "M_direction"),
)

# Concise biochemical abbreviations for canvas clarity (full standard names remain in pathway_map.tsv)
SHORT_LABELS: dict[str, str] = {
    "GLC": "Glc", "GLC6P": "G6P", "FRU6P": "F6P", "FRU16BP": "F1,6BP",
    "DHAP": "DHAP", "GAP": "GAP", "13BPG": "1,3BPG", "3PG": "3PG", "2PG": "2PG",
    "PEP": "PEP", "PYR": "Pyr", "LAC": "Lac", "GLC1P": "G1P",
    "RUL5P": "Ru5P", "RIB5P": "R5P", "XYL5P": "Xu5P", "SED7P": "S7P", "E4P": "E4P",
    "SER": "Ser", "GLY": "Gly",
    "G3P": "G3P", "GLYCEROL": "Glycerol", "LPA": "LPA", "PA": "PA", "DAG": "DAG", "MAG": "MAG", "TAG": "TAG",
    "ACE": "Acetate", "ACCOA": "Ac-CoA", "OAC": "OAA", "CIT": "Cit", "ICIT": "iCit",
    "AKG": "α-KG", "SUCCOA": "Suc-CoA", "SUC": "Suc", "FUM": "Fum", "MAL": "Mal",
    "ASP": "Asp", "GLU": "Glu", "GLN": "Gln", "ALA": "Ala", "ASN": "Asn",
    "ACACCOA": "AcAc-CoA", "ACAC": "AcAc", "BHB": "β-HB",
    "MALCOA": "Mal-CoA", "STE": "Stearate", "OLE": "Oleate",
    "PAL": "Palmitate", "PALO": "Palmitoleate", "TGPAL": "TG-Pal",
    "PALCOA": "Pal-CoA", "PALCAR": "Pal-Car", "ACCAR": "Ac-Car", "FA": "FFA",
    "CP": "Carb-P", "ORN": "Orn", "CITR": "Citr", "ARG": "Arg",
    "P5C": "P5C", "PRO": "Pro", "GAA": "GAA", "CRE": "Cr", "PCR": "PCr", "CRN": "Crn",
    "IMP": "IMP", "AMP": "AMP", "ADP": "ADP", "ADO": "Ado", "INO": "Ino", "GMP": "GMP", "GUO": "Guo",
    "UMP": "UMP", "CMP": "CMP", "CYTD": "Cyd", "URI": "Urd", "DCYTD": "dCyd",
    "FRU": "Fru", "UDPGLCNAC": "UDP-GlcNAc",
    "VAL": "Val", "LEU": "Leu", "ILE": "Ile", "PROPCAR": "Prop-Car",
    "CHOL": "Chol", "CE": "CE", "BILE": "Bile Acids",
    "LIN": "Linoleate", "AA": "Arachidonate",
    "PL": "Phospholipids", "SPH": "Sphingolipids", "ETHER": "Ether Lipids",
    "MET": "Met", "HCY": "Hcy", "CYS": "Cys", "GSH": "GSH", "OXOPRO": "PyroGlu", "TAU": "Taurine",
    "POLYAMINE": "Polyamines", "ALLANTOIN": "Allantoin",
    "TRP": "Trp", "KYN": "Kyn", "NAM": "NAM", "NAD_POOL": "NAD(P)",
}


def _short(flux_id: str) -> str:
    for suffix in ("_c", "_m"):
        if flux_id.endswith(suffix):
            return flux_id[: -len(suffix)]
    return flux_id


def _species_short(species: str, fallback: str = "") -> str:
    base = species.rsplit("_", 1)[0]
    if base in SHORT_LABELS:
        return SHORT_LABELS[base]
    return fallback or base


def _text_width(text: str, fontsize: float, bold: bool = False, unit: float = UNIT) -> float:
    return estimate_text_width(text, fontsize, unit, bold=bold)


def _node_clearance(layout: Layout, species: str) -> float:
    kind = layout.nodes[species].kind
    if kind == "metabolite":
        return NODE_R + 0.30
    if kind == "boundary":
        return 1.38
    return 1.98


def _agg_card_dims(agg) -> tuple[float, float]:
    """Adaptive aggregate card (width, height) with generous badge room."""
    if agg.region == "oxphos":
        if agg.flux_id in ("OXPHOS_CV", "THERMOGENESIS_UCP"):
            return 11.8, 4.0
        return 9.1, 4.0
    if agg.region == "other":
        return 13.4, 4.05
    w = max(12.2, min(14.4, _text_width(agg.label or agg.flux_id, 4.2, bold=True) + 4.8))
    return w, 4.0


# ---------------------------------------------------------------------------
# Edge Geometry Computation
# ---------------------------------------------------------------------------

def _edge_color(
    row: pd.Series,
    palette: SemanticPalette,
    cmap,
    *,
    metro: bool = False,
    color_mode: str = "semantic",
    layout: Layout | None = None,
) -> str:
    if row["draw_style"] == "not_assessed":
        return palette.faint
    if metro and color_mode == "pathway":
        reg_id = str(row.get("region", "") or "")
        fallback = layout.regions[reg_id].accent if (layout is not None and reg_id in layout.regions) else ""
        return pathway_family_color(reg_id, fallback_accent=fallback)
    if row["draw_style"] == "omics_only":
        return "#5B6C82" if metro else palette.omics_grey
    value = row["value"]
    if pd.isna(value):
        return palette.edge_grey
    return cmap(DEFAULT_NORM(math.copysign(max(abs(value), MIN_EDGE), value)))


def _edge_geometry(layout: Layout, flux_id: str, *, metro: bool = False) -> dict:
    edge = layout.edges[flux_id]
    geometry = layout.geometry(pathway_reactions()[flux_id])
    subs = geometry["substrates"]
    pros = geometry["products"]
    curves: list[dict] = []
    metro_paths: list[dict] = []
    incident: list[tuple[str, tuple[float, float], float]] = []

    def add_incident(species: str, pt: tuple[float, float], towards: tuple[float, float]) -> None:
        angle = math.atan2(towards[1] - pt[1], towards[0] - pt[0])
        incident.append((species, pt, angle))

    if metro:
        if len(subs) == 1 and len(pros) == 1 and not edge.via:
            (s_sp, p0), (p_sp, p1) = subs[0], pros[0]
            bow = edge.offset * 1.65 + edge.curve
            if edge.hub is not None and abs(bow) < 1e-4:
                raw = route_octilinear(p0, p1, via=(edge.hub,), prefer="hv_first")
            else:
                raw = route_octilinear(p0, p1, bend=bow, prefer="hv_first")
            smooth = fillet_polyline(raw, radius=1.15, arc_steps=6)
            t0, t1 = _node_clearance(layout, s_sp), _node_clearance(layout, p_sp)
            drawn = trim_polyline(smooth, t0, t1)
            hub, _ = sample_polyline_at(drawn, polyline_length(drawn) * 0.5)
            metro_paths.append({"points": raw, "smooth": smooth, "drawn": drawn, "trim0": t0, "trim1": t1, "role": "single"})
            add_incident(s_sp, p0, smooth[1] if len(smooth) > 1 else p1)
            add_incident(p_sp, p1, smooth[-2] if len(smooth) > 1 else p0)
        elif len(subs) == 1 and len(pros) == 1 and edge.via:
            (s_sp, p0), (p_sp, p1) = subs[0], pros[0]
            via_pts = [*edge.via, geometry["hub"]]
            raw = route_octilinear(p0, p1, via=via_pts, prefer="hv_first")
            smooth = fillet_polyline(raw, radius=1.15, arc_steps=6)
            t0, t1 = _node_clearance(layout, s_sp), _node_clearance(layout, p_sp)
            drawn = trim_polyline(smooth, t0, t1)
            hub = geometry["hub"]
            metro_paths.append({"points": raw, "smooth": smooth, "drawn": drawn, "trim0": t0, "trim1": t1, "role": "single"})
            add_incident(s_sp, p0, smooth[1] if len(smooth) > 1 else p1)
            add_incident(p_sp, p1, smooth[-2] if len(smooth) > 1 else p0)
        else:
            hub = geometry["hub"]
            for s_sp, p in subs:
                raw = route_octilinear(p, hub, prefer="hv_first")
                smooth = fillet_polyline(raw, radius=1.0, arc_steps=5)
                t0 = _node_clearance(layout, s_sp)
                drawn = trim_polyline(smooth, t0, 0.0)
                metro_paths.append({"points": raw, "smooth": smooth, "drawn": drawn, "trim0": t0, "trim1": 0.0, "role": "sub"})
                add_incident(s_sp, p, smooth[1] if len(smooth) > 1 else hub)
            for p_sp, p in pros:
                raw = route_octilinear(hub, p, prefer="diag_first")
                smooth = fillet_polyline(raw, radius=1.0, arc_steps=5)
                t1 = _node_clearance(layout, p_sp)
                drawn = trim_polyline(smooth, 0.0, t1)
                metro_paths.append({"points": raw, "smooth": smooth, "drawn": drawn, "trim0": 0.0, "trim1": t1, "role": "pro"})
                add_incident(p_sp, p, smooth[-2] if len(smooth) > 1 else hub)
        samples = [pt for mp in metro_paths for pt in sample_polyline_points(mp["smooth"], 10)]
        ins, outs = geometry["secondary"]
        note = ""
        if ins or outs:
            left = ", ".join(_species_short(s) for s in dict.fromkeys(ins))
            right = ", ".join(_species_short(s) for s in dict.fromkeys(outs))
            note = f"{left} → {right}".strip()
        return {"hub": hub, "curves": curves, "metro_paths": metro_paths, "samples": samples, "incident": incident, "note": note}

    if len(subs) == 1 and len(pros) == 1 and not edge.via:
        (s_sp, p0), (p_sp, p1) = subs[0], pros[0]
        dx, dy = p1[0] - p0[0], p1[1] - p0[1]
        length = math.hypot(dx, dy) or 1.0
        nx, ny = -dy / length, dx / length
        bow = edge.offset * 1.65 + edge.curve
        mid_chord = ((p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2)
        if abs(bow) > 1e-4:
            ctrl = (mid_chord[0] + nx * bow, mid_chord[1] + ny * bow)
        elif edge.hub is not None:
            ctrl = (2 * edge.hub[0] - mid_chord[0], 2 * edge.hub[1] - mid_chord[1])
        else:
            ctrl = mid_chord
        hub = quad_point(p0, ctrl, p1, 0.5)
        curves.append({
            "p0": p0, "c": ctrl, "p1": p1,
            "trim0": _node_clearance(layout, s_sp), "trim1": _node_clearance(layout, p_sp),
            "role": "single",
        })
        add_incident(s_sp, p0, ctrl if math.dist(p0, ctrl) > 1e-3 else p1)
        add_incident(p_sp, p1, ctrl if math.dist(p1, ctrl) > 1e-3 else p0)
    elif len(subs) == 1 and len(pros) == 1 and edge.via:
        (s_sp, p0), (p_sp, p1) = subs[0], pros[0]
        pts = [p0, *edge.via, geometry["hub"], p1]
        anchors = [pts[0]] + [
            ((pts[i][0] + pts[i + 1][0]) / 2, (pts[i][1] + pts[i + 1][1]) / 2)
            for i in range(1, len(pts) - 2)
        ] + [pts[-1]]
        ctrls = pts[1:-1]
        for idx, ctrl in enumerate(ctrls):
            a, b = anchors[idx], anchors[idx + 1]
            role = "single" if len(ctrls) == 1 else ("sub" if idx == 0 else ("pro" if idx == len(ctrls) - 1 else "mid"))
            curves.append({
                "p0": a, "c": ctrl, "p1": b,
                "trim0": _node_clearance(layout, s_sp) if idx == 0 else 0.0,
                "trim1": _node_clearance(layout, p_sp) if idx == len(ctrls) - 1 else 0.0,
                "role": role,
            })
        hub = geometry["hub"]
        add_incident(s_sp, p0, pts[1])
        add_incident(p_sp, p1, pts[-2])
    else:
        hub = geometry["hub"]
        if subs and pros:
            sc = (sum(p[0] for _, p in subs) / len(subs), sum(p[1] for _, p in subs) / len(subs))
            pc = (sum(p[0] for _, p in pros) / len(pros), sum(p[1] for _, p in pros) / len(pros))
            ux, uy = pc[0] - sc[0], pc[1] - sc[1]
            un = math.hypot(ux, uy) or 1.0
            ux, uy = ux / un, uy / un
        else:
            ux, uy = 1.0, 0.0
        for s_sp, p in subs:
            d = math.dist(p, hub)
            bend = 0.35 * d if len(subs) > 1 else 0.0
            ctrl = (hub[0] - ux * bend, hub[1] - uy * bend) if bend > 0 else ((p[0] + hub[0]) / 2, (p[1] + hub[1]) / 2)
            curves.append({
                "p0": p, "c": ctrl, "p1": hub,
                "trim0": _node_clearance(layout, s_sp), "trim1": 0.0, "role": "sub",
            })
            add_incident(s_sp, p, ctrl if math.dist(p, ctrl) > 1e-3 else hub)
        for p_sp, p in pros:
            d = math.dist(hub, p)
            bend = 0.35 * d if len(pros) > 1 else 0.0
            ctrl = (hub[0] + ux * bend, hub[1] + uy * bend) if bend > 0 else ((hub[0] + p[0]) / 2, (hub[1] + p[1]) / 2)
            curves.append({
                "p0": hub, "c": ctrl, "p1": p,
                "trim0": 0.0, "trim1": _node_clearance(layout, p_sp), "role": "pro",
            })
            add_incident(p_sp, p, ctrl if math.dist(p, ctrl) > 1e-3 else hub)

    samples = [pt for cv in curves for pt in sample_quad(cv["p0"], cv["c"], cv["p1"], 10)]
    ins, outs = geometry["secondary"]
    note = ""
    if ins or outs:
        left = ", ".join(_species_short(s) for s in dict.fromkeys(ins))
        right = ", ".join(_species_short(s) for s in dict.fromkeys(outs))
        note = f"{left} → {right}".strip()
    return {"hub": hub, "curves": curves, "samples": samples, "incident": incident, "note": note}


# ---------------------------------------------------------------------------
# Smart Collision-Free Label & Micro-Strip Placement
# ---------------------------------------------------------------------------

def _node_label_placements(
    layout: Layout,
    geometries: dict[str, dict],
    bounds: tuple[float, float, float, float],
    font_scale: float = 1.0,
    unit_in: float = UNIT,
) -> dict[tuple[str, int], tuple[float, float, str, str, tuple]]:
    angles_by_pos: dict[tuple[str, int], list[float]] = {}
    for geom in geometries.values():
        for sp, pt, ang in geom["incident"]:
            node = layout.nodes[sp]
            idx = min(range(len(node.at)), key=lambda i: math.dist(node.at[i], pt))
            angles_by_pos.setdefault((sp, idx), []).append(ang)

    dirs = [
        (math.pi / 2, 0.0, NODE_R + 0.48, "center", "top", 0.00),
        (0.0, NODE_R + 0.50, 0.0, "left", "center", 0.04),
        (math.pi, -(NODE_R + 0.50), 0.0, "right", "center", 0.05),
        (-math.pi / 2, 0.0, -(NODE_R + 0.48), "center", "bottom", 0.08),
        (math.pi / 4, NODE_R * 0.90 + 0.3, NODE_R * 0.90 + 0.25, "left", "top", 0.14),
        (3 * math.pi / 4, -(NODE_R * 0.90 + 0.3), NODE_R * 0.90 + 0.25, "right", "top", 0.15),
        (-math.pi / 4, NODE_R * 0.90 + 0.3, -(NODE_R * 0.90 + 0.25), "left", "bottom", 0.17),
        (-3 * math.pi / 4, -(NODE_R * 0.90 + 0.3), -(NODE_R * 0.90 + 0.25), "right", "bottom", 0.18),
    ]
    bx_min, by_min, bx_max, by_max = bounds
    fs = 4.6 * font_scale
    out: dict[tuple[str, int], tuple[float, float, str, str, tuple]] = {}
    for sp, node in layout.nodes.items():
        if node.kind != "metabolite":
            continue
        text = node.label or _species_short(sp)
        w = _text_width(text, fs, unit=unit_in) + 0.38
        h = em_units(fs, unit_in) * 1.28
        for idx, (x, y) in enumerate(node.at):
            inc = angles_by_pos.get((sp, idx), [])
            best = None
            best_score = -1e9
            for ang, ox, oy, ha, va, pref in dirs:
                clearance = min(abs(math.atan2(math.sin(ang - a), math.cos(ang - a))) for a in inc) if inc else math.pi
                lx, ly = x + ox, y + oy
                bx0 = lx - (w / 2 if ha == "center" else (w if ha == "right" else 0.0))
                by0 = ly - (h / 2 if va == "center" else (h if va == "bottom" else 0.0))
                bbox = (bx0, by0, bx0 + w, by0 + h)
                oob = -6.0 if (bbox[0] < bx_min + 0.4 or bbox[2] > bx_max - 0.4 or bbox[1] < by_min + 0.4 or bbox[3] > by_max - 0.4) else 0.0
                score = min(clearance, 1.35) - pref + oob
                if score > best_score:
                    best_score = score
                    best = (lx, ly, ha, va, bbox)
            out[(sp, idx)] = best
    return out


def _flux_group_indices(view: pd.DataFrame) -> dict[str, int]:
    """Assign integer indices (1, 2, ...) to coupled flux groups with directional F-axis evidence."""
    rxns = view[(view["element_type"] == "reaction") & view["F_state"].isin(("changed", "suggestive"))]
    groups: list[str] = []
    for g in rxns["flux_group"].dropna():
        gs = str(g).strip()
        if gs and gs not in groups:
            groups.append(gs)
    return {g: i + 1 for i, g in enumerate(groups)}


def _strip_badges(
    row: pd.Series,
    group_indices: dict[str, int],
    *,
    top_rank: int = 10,
    show_rank: bool = True,
) -> dict:
    """Return structured badge descriptors for a reaction card."""
    label = _short(row["element_id"])
    if row["draw_style"] == "not_assessed":
        return {"label": label, "tier": "", "rank": None, "group": None}
    tier = str(row["flux_claim"]) if row["flux_claim"] in ("L1", "L2", "L3") else ""
    rank = None
    if show_rank and pd.notna(row["rank"]):
        r_int = int(row["rank"])
        if r_int <= top_rank and row["key"] in ("key", "l1", "top"):
            rank = r_int
    group = None
    if row["flux_group"] and row["F_state"] in ("changed", "suggestive"):
        group = group_indices.get(str(row["flux_group"]).strip())
    return {"label": label, "tier": tier, "rank": rank, "group": group}


def _is_compact_label_only(row: pd.Series, mode: str) -> bool:
    """In publication mode, non-key reactions collapse to clean enzyme labels for print legibility."""
    if row["draw_style"] == "not_assessed":
        return True
    if mode == "publication" and row["key"] not in ("key", "l1") and row["F_state"] not in ("changed", "suggestive"):
        return True
    return False


def _val_sign(val) -> int:
    if pd.isna(val):
        return 0
    f = float(val)
    return 1 if f > 0 else (-1 if f < 0 else 0)


def _extract_strip_states(row: pd.Series, strip_mode: str = "omics") -> list[tuple[str, int]]:
    """提取四格微条的 ``[(state, direction), ...]``：

    - ``strip_mode="omics"``（默认）：4 层组学 ``[Flux | Protein | mRNA | Metab]``；
    - ``strip_mode="axes"``：S3 四机制轴 ``[Flux | Enzyme Capacity | Regulation | Pool]``。
    """
    if strip_mode == "axes":
        return [
            (
                str(row[sc]) if (sc in row.index and pd.notna(row[sc])) else "not_measured",
                int(row[dc]) if (dc in row.index and pd.notna(row[dc])) else 0,
            )
            for sc, dc in STRIP_AXIS_COLS
        ]

    # 1. Flux
    f_st = str(row["F_state"]) if ("F_state" in row.index and pd.notna(row["F_state"])) else "not_measured"
    f_dir = int(row["F_direction"]) if ("F_direction" in row.index and pd.notna(row["F_direction"])) else 0

    # 2. Protein & 3. mRNA（若视图由 map_view 生成则直接取真实分层读数；合成最小表回退 C_state）
    has_layer_cols = (
        "P_state" in row.index and pd.notna(row["P_state"])
    ) or (
        "R_layer_state" in row.index and pd.notna(row["R_layer_state"])
    )
    if has_layer_cols:
        p_st = str(row["P_state"]) if ("P_state" in row.index and pd.notna(row["P_state"])) else "not_measured"
        p_dir = _val_sign(row["P_value"]) if "P_value" in row.index else 0
        if p_st == "summary":
            p_st = "changed" if abs(float(row["P_value"] or 0.0)) >= 0.58 else "unchanged"
        r_st = str(row["R_layer_state"]) if ("R_layer_state" in row.index and pd.notna(row["R_layer_state"])) else "not_measured"
        r_dir = _val_sign(row["R_value"]) if "R_value" in row.index else 0
        if r_st == "summary":
            r_st = "changed" if abs(float(row["R_value"] or 0.0)) >= 0.58 else "unchanged"
    else:
        p_st = str(row["C_state"]) if ("C_state" in row.index and pd.notna(row["C_state"])) else "not_measured"
        p_dir = int(row["C_direction"]) if ("C_direction" in row.index and pd.notna(row["C_direction"])) else 0
        r_st, r_dir = "not_measured", 0

    # 4. Metabolite（优先取沿生理方向的底物/产物评估 M_state，未测时回退 Mx_state）
    m_st = str(row["M_state"]) if ("M_state" in row.index and pd.notna(row["M_state"])) else "not_measured"
    m_dir = int(row["M_direction"]) if ("M_direction" in row.index and pd.notna(row["M_direction"])) else 0
    if m_st in ("not_measured", "no_knowledge") and "Mx_state" in row.index and pd.notna(row["Mx_state"]):
        mx_st = str(row["Mx_state"])
        if mx_st not in ("not_measured", "no_knowledge"):
            m_st = mx_st if mx_st != "summary" else "unchanged"
            m_dir = _val_sign(row["Mx_value"]) if "Mx_value" in row.index else 0

    return [(f_st, f_dir), (p_st, p_dir), (r_st, r_dir), (m_st, m_dir)]


def _strip_size(
    row: pd.Series,
    group_indices: dict[str, int],
    *,
    top_rank: int = 10,
    mode: str = "atlas",
    font_scale: float = 1.0,
    unit_in: float = UNIT,
) -> tuple[float, float]:
    show_rank = mode != "publication"
    info = _strip_badges(row, group_indices, top_rank=top_rank, show_rank=show_rank)
    if _is_compact_label_only(row, mode):
        fs = 3.9 * font_scale
        return _text_width(info["label"], fs, unit=unit_in) + 0.55, 0.98 * font_scale
    is_hi = row["key"] in ("key", "l1", "top")
    bold = row["key"] in ("key", "l1")
    lbl_fs = 4.35 * font_scale
    badge_scale = 1.02 * font_scale
    total_groups = max(len(group_indices), 1)
    w = 4 * CELL + 0.46 + _text_width(info["label"], lbl_fs, bold=bold, unit=unit_in) + 0.48
    gap = 0.34
    if info["tier"]:
        w += gap + tier_badge_width(info["tier"], scale=badge_scale, unit_in=unit_in)
    if info["group"] is not None:
        w += gap + coupled_badge_width(info["group"], scale=badge_scale, total_groups=total_groups, unit_in=unit_in)
    w += 0.38
    h = CELL + (AXIS_HDR_H + 0.44 if is_hi else AXIS_BAR_H + 0.42)
    return w, h


def _overlap(a: tuple, b: tuple) -> float:
    dx = min(a[2], b[2]) - max(a[0], b[0])
    dy = min(a[3], b[3]) - max(a[1], b[1])
    return dx * dy if dx > 0 and dy > 0 else 0.0


def _points_inside(points: list[tuple[float, float]], rect: tuple) -> int:
    x0, y0, x1, y1 = rect
    return sum(1 for px, py in points if x0 <= px <= x1 and y0 <= py <= y1)


def _note_rect(geometry: dict) -> tuple | None:
    if not geometry["note"]:
        return None
    hx, hy = geometry["hub"]
    w = _text_width(geometry["note"], NOTE_SIZE) + 0.35
    return (hx - w / 2, hy + 0.4, hx + w / 2, hy + 1.25)


def _obstacles(
    layout: Layout,
    node_labels: dict,
    active_species: set[str],
    active_aggs: set[str],
    active_regions: set[str],
    *,
    metro: bool = False,
    comp_boxes: dict[str, tuple[float, float, float, float]] | None = None,
) -> list[tuple[tuple, float]]:
    out = []
    for species, node in layout.nodes.items():
        if species not in active_species:
            continue
        for idx, (x, y) in enumerate(node.at):
            if node.kind == "metabolite":
                out.append(((x - NODE_R - 0.30, y - NODE_R - 0.30, x + NODE_R + 0.30, y + NODE_R + 0.30), 4.6))
                if (species, idx) in node_labels:
                    out.append((node_labels[(species, idx)][4], 3.8))
            elif node.kind == "boundary":
                out.append(((x - 2.7, y - 0.95, x + 2.7, y + 0.95), 2.8))
            else:
                out.append(((x - 2.8, y - 0.95, x + 2.8, y + 0.95), 2.8))
    for aid, agg in layout.aggregates.items():
        if aid not in active_aggs:
            continue
        x, y = agg.xy
        w, h = _agg_card_dims(agg)
        out.append(((x - w / 2 - 0.30, y - h / 2 - 0.30, x + w / 2 + 0.30, y + h / 2 + 0.30), 3.6))
    for rid, region in layout.regions.items():
        if rid not in active_regions:
            continue
        x, y = region.label_xy
        out.append(((x - 0.4, y - 0.95, x + _text_width(region.label, 5.5, bold=True) + 1.5, y + 0.95), 2.8))
    boxes = comp_boxes if comp_boxes is not None else {cid: tuple(c["box"]) for cid, c in layout.compartments.items()}
    for cid, box in boxes.items():
        comp = layout.compartments[cid]
        x0, y0, x1, y1 = box
        lw = _text_width(comp.get("label", ""), 5.0, bold=True) + 1.6
        if cid == "p":
            out.append(((x1 - lw - 0.6, y0 + 0.5, x1 - 0.6, y0 + 2.4), 2.2))
        elif metro and "oxphos" in active_regions:
            cx_m = (x0 + x1) * 0.5
            out.append(((cx_m - lw * 0.5 - 0.8, y0 + 0.6, cx_m + lw * 0.5 + 0.8, y0 + 3.4), 4.2))
        elif metro and not ({"pyruvate_tca", "ketone"} & active_regions):
            out.append(((x0 + 0.8, y0 + 0.6, x0 + lw + 1.0, y0 + 2.8), 3.0))
        elif metro:
            out.append(((x1 - lw - 1.0, y1 - 2.6, x1 - 0.6, y1 - 0.5), 3.0))
        else:
            # 线粒体标签移至左下角避开 BHB_trans
            out.append(((x0 + 0.6, y1 - 2.5, x0 + lw + 0.8, y1 - 0.6), 2.2))
    return out


def _place_strips(
    reactions: pd.DataFrame,
    layout: Layout,
    geometries: dict[str, dict],
    node_labels: dict,
    group_indices: dict[str, int],
    *,
    active_species: set[str],
    active_aggs: set[str],
    active_regions: set[str],
    bounds: tuple[float, float, float, float],
    top_rank: int = 10,
    mode: str = "atlas",
    font_scale: float = 1.0,
    unit_in: float = UNIT,
    metro: bool = False,
    comp_boxes: dict[str, tuple[float, float, float, float]] | None = None,
) -> dict[str, tuple[float, float]]:
    obstacles = _obstacles(
        layout, node_labels, active_species, active_aggs, active_regions,
        metro=metro, comp_boxes=comp_boxes,
    ) + [(r, 2.4) for g in geometries.values() if (r := _note_rect(g))]
    all_samples = [pt for g in geometries.values() for pt in g["samples"]]
    placed: list[tuple] = []
    out: dict[str, tuple[float, float]] = {}
    bx0, by0, bx1, by1 = bounds
    priority = {"key": 0, "l1": 1, "top": 2}
    rows = reactions.copy()
    rows["order"] = (
        rows["key"].map(priority).fillna(3)
        + rows["draw_style"].eq("omics_only") * 0.5
        + rows["draw_style"].eq("not_assessed") * 2.0
    )
    for _, row in rows.sort_values("order", kind="stable").iterrows():
        flux_id = row["element_id"]
        hx, hy = geometries[flux_id]["hub"]
        w, h = _strip_size(row, group_indices, top_rank=top_rank, mode=mode, font_scale=font_scale, unit_in=unit_in)
        explicit = layout.edges[flux_id].strip
        if explicit is not None:
            candidates = [explicit]
        elif _is_compact_label_only(row, mode):
            g = 0.65
            candidates = [
                (g, 0.0), (-g - w, 0.0), (-w / 2, -1.10), (-w / 2, 1.10),
                (g, -1.25), (g, 1.25), (-g - w, -1.25), (-g - w, 1.25),
            ]
        else:
            g = 0.90
            candidates = [
                (g, 0.0), (-g - w, 0.0), (-w / 2, -1.75), (-w / 2, 1.75),
                (g, -1.65), (g, 1.65), (-g - w, -1.65), (-g - w, 1.65),
                (g + 1.35, 0.0), (-g - 1.35 - w, 0.0), (-w / 2, -3.00), (-w / 2, 3.00),
                (g, -3.10), (g, 3.10), (-g - w, -3.10), (-g - w, 3.10),
            ]
        best, best_score = None, math.inf
        for idx, (dx, dy) in enumerate(candidates):
            # 若左上角带 #n 角标，顶部多预留 0.55 避免撞上相邻元素
            top_pad = 0.55 if (mode != "publication" and row["key"] in ("key", "l1", "top")) else 0.0
            rect = (hx + dx, hy + dy - h / 2 - top_pad, hx + dx + w, hy + dy + h / 2)
            dist_pen = math.hypot(dx + w / 2, dy) * 0.03 + idx * 0.015
            score = dist_pen
            score += sum(_overlap(rect, r) * weight for r, weight in obstacles)
            score += sum(_overlap(rect, r) * 5.8 for r in placed)
            score += _points_inside(all_samples, rect) * 0.16
            if rect[0] < bx0 + 0.4 or rect[2] > bx1 - 0.3 or rect[1] < by0 + 0.4 or rect[3] > by1 - 0.3:
                score += 80.0
            if score < best_score:
                best, best_score = (dx, dy, rect), score
        dx, dy, rect = best
        placed.append(rect)
        out[flux_id] = (hx + dx, hy + dy)
    return out


# ---------------------------------------------------------------------------
# Active Region Selection & Viewport Computation for Publication Mode
# ---------------------------------------------------------------------------

def _select_active_regions(
    view: pd.DataFrame,
    layout: Layout,
    focus_regions: Sequence[str] | None,
) -> tuple[set[str], list[str]]:
    """Return (active_region_ids, folded_inactive_labels)."""
    if focus_regions:
        unknown = [r for r in focus_regions if r not in layout.regions]
        if unknown:
            raise ValueError(f"pathway_map.regions 包含未声明的区域：{', '.join(unknown)}")
        active = set(focus_regions)
    else:
        rxns = view[view["element_type"] == "reaction"]
        aggs = view[view["element_type"] == "aggregate"]
        active = set()
        for rid in layout.regions:
            sub_r = rxns[rxns["region"] == rid]
            sub_a = aggs[aggs["region"] == rid]
            has_flux = bool((sub_r["draw_style"] == "flux").any())
            has_l2 = bool(sub_r["flux_claim"].isin(("L2", "L3")).any() or sub_a["flux_claim"].isin(("L2", "L3")).any())
            has_l1 = bool(sub_r["flux_claim"].eq("L1").any() or sub_a["flux_claim"].eq("L1").any())
            # 出版模式自动聚焦：保留有 MFA 通量或 L2+ 证据的核心区；若仅有零星 L1 的大块外围区也可排除以防稀疏
            if has_flux or has_l2 or (has_l1 and rid not in ("other", "nucleotide", "urea_creatine")):
                active.add(rid)
        if not active:
            active = set(layout.regions)

    folded = [layout.regions[r].label for r in layout.regions if r not in active]
    return active, folded


def _compute_clipped_compartments(
    layout: Layout,
    active_regions: set[str],
    active_node_positions: set[tuple[str, int]] | None,
    active_aggs: set[str],
    *,
    mode: str = "atlas",
) -> dict[str, tuple[float, float, float, float]]:
    """计算可见细胞器容器的精确范围；子图模式下自动紧裁线粒体空白侧，避免拖入大块空白。"""
    out: dict[str, tuple[float, float, float, float]] = {}
    for cid, comp in layout.compartments.items():
        mx0, my0, mx1, my1 = (float(v) for v in comp["box"])
        if cid == "p":
            if "peroxisome" in active_regions:
                out[cid] = (mx0, my0, mx1, my1)
            continue

        # 线粒体容器：仅当含线粒体相关区域或有活跃节点落在其内部时绘制
        mito_regs = {"pyruvate_tca", "ketone", "oxphos"}
        if mode == "atlas" or mito_regs <= active_regions:
            if mito_regs & active_regions:
                out[cid] = (mx0, my0, mx1, my1)
            continue

        in_xs: list[float] = []
        in_ys: list[float] = []
        if active_node_positions is not None:
            for sp, idx in active_node_positions:
                x, y = layout.nodes[sp].at[idx]
                if mx0 - 1.2 <= x <= mx1 + 1.2 and my0 - 1.2 <= y <= my1 + 1.2:
                    in_xs.append(x)
                    in_ys.append(y)
        for aid in active_aggs:
            x, y = layout.aggregates[aid].xy
            w, h = _agg_card_dims(layout.aggregates[aid])
            if mx0 - 1.0 <= x <= mx1 + 1.0 and my0 - 1.0 <= y <= my1 + 1.0:
                in_xs.extend([x - w / 2, x + w / 2])
                in_ys.extend([y - h / 2, y + h / 2])
        for rid in mito_regs & active_regions:
            if rid in layout.regions:
                lx, ly = layout.regions[rid].label_xy
                if mx0 - 1.0 <= lx <= mx1 + 1.0 and my0 - 1.0 <= ly <= my1 + 1.0:
                    in_xs.extend([lx, lx + _text_width(layout.regions[rid].label, 5.4, bold=True) + 2.0])
                    in_ys.append(ly)

        if not in_xs or not in_ys:
            if mito_regs & active_regions:
                out[cid] = (mx0, my0, mx1, my1)
            continue

        cx0 = max(mx0, min(in_xs) - 4.8)
        cx1 = min(mx1, max(in_xs) + 6.2)
        cy0 = max(my0, min(in_ys) - 4.2)
        cy1 = min(my1, max(in_ys) + 4.8)
        if cx1 - cx0 < 16.0:
            mid = (cx0 + cx1) * 0.5
            cx0, cx1 = max(mx0, mid - 8.0), min(mx1, mid + 8.0)
        if cy1 - cy0 < 14.0:
            mid = (cy0 + cy1) * 0.5
            cy0, cy1 = max(my0, mid - 7.0), min(my1, mid + 7.0)
        out[cid] = (cx0, cy0, cx1, cy1)
    return out


def _compute_publication_viewport(
    layout: Layout,
    active_regions: set[str],
    active_node_positions: set[tuple[str, int]],
    active_aggs: set[str],
    comp_boxes: dict[str, tuple[float, float, float, float]],
) -> tuple[float, float, float, float]:
    xs: list[float] = []
    ys: list[float] = []
    for rid in active_regions:
        reg = layout.regions[rid]
        if reg.box is not None:
            xs.extend([reg.box[0], reg.box[2]])
            ys.extend([reg.box[1], reg.box[3]])
        xs.append(reg.label_xy[0])
        ys.append(reg.label_xy[1])
    for sp, idx in active_node_positions:
        x, y = layout.nodes[sp].at[idx]
        xs.append(x)
        ys.append(y)
    for aid in active_aggs:
        x, y = layout.aggregates[aid].xy
        w, h = _agg_card_dims(layout.aggregates[aid])
        xs.extend([x - w / 2, x + w / 2])
        ys.extend([y - h / 2, y + h / 2])
    for bx0, by0, bx1, by1 in comp_boxes.values():
        xs.extend([bx0, bx1])
        ys.extend([by0, by1])
    if not xs or not ys:
        return 0.0, 0.0, layout.width, layout.height
    pad = 3.0
    vx0 = max(0.0, min(xs) - pad)
    vy0 = max(0.0, min(ys) - pad)
    vx1 = min(layout.width, max(xs) + pad + 2.5)
    vy1 = min(layout.height, max(ys) + pad)
    # 极窄单栏子图（如 glycerolipid）适度向右拓宽留白，便于四轴微条与紧凑图例优雅舒展
    min_vw = 44.0
    if vx1 - vx0 < min_vw:
        need = min_vw - (vx1 - vx0)
        add_right = min(layout.width - vx1, need * 0.80)
        vx1 += add_right
        rem = need - add_right
        if rem > 0:
            vx0 = max(0.0, vx0 - rem)
    return vx0, vy0, vx1, vy1


# ---------------------------------------------------------------------------
# Rendering: Zones, Organelles, Coupled Tracks, Edges, Strips, Nodes, Cards
# ---------------------------------------------------------------------------

def _draw_regions_and_compartments(
    ax,
    layout: Layout,
    active_regions: set[str],
    *,
    palette: SemanticPalette,
    font_scale: float = 1.0,
    metro: bool = False,
    color_mode: str = "semantic",
    comp_boxes: dict[str, tuple[float, float, float, float]] | None = None,
) -> None:
    boxes = comp_boxes if comp_boxes is not None else {
        cid: tuple(comp["box"]) for cid, comp in layout.compartments.items()
        if (cid == "p" and "peroxisome" in active_regions)
        or (cid != "p" and ({"pyruvate_tca", "ketone", "oxphos"} & active_regions))
    }
    if metro:
        # 地铁风格：区室退为极淡无抢戏的浅灰圆角底板（100% 视觉聚焦于八向轨道与换乘站）
        for rid, region in layout.regions.items():
            if rid not in active_regions:
                continue
            if region.box is not None and rid != "oxphos":
                draw_metro_zone(
                    ax, region.box, "",
                    face="#F3F5F8", edge="#E2E7EF", lw=0.50, rounding=2.0, zorder=0.2,
                )
        for cid, box in boxes.items():
            comp = layout.compartments[cid]
            if cid == "p":
                draw_metro_zone(
                    ax, box, comp.get("label", ""),
                    face="#EBF4F2", edge="#B4D7CF", lw=0.75, inner_edge="#9CCDC2",
                    label_color="#336D65", label_pos="top_right", rounding=2.4,
                    fontsize=5.0 * font_scale, zorder=0.5,
                )
            else:
                if "oxphos" in active_regions:
                    mito_pos = "top_center"
                elif not ({"pyruvate_tca", "ketone"} & active_regions):
                    mito_pos = "top_left"
                else:
                    mito_pos = "bottom_right"
                draw_metro_zone(
                    ax, box, comp.get("label", ""),
                    face="#ECEEF2", edge="#CFD6E0", lw=0.80, inner_edge="#C5CDD8",
                    label_color="#64748B", label_pos=mito_pos, rounding=2.6,
                    fontsize=5.0 * font_scale, zorder=0.5,
                )
        if "oxphos" in active_regions and "oxphos" in layout.regions and layout.regions["oxphos"].box is not None:
            reg = layout.regions["oxphos"]
            draw_metro_zone(ax, reg.box, "", face="#E6E9F0", edge="#CBD2DD", lw=0.55, rounding=1.6, zorder=0.62)
        for rid, region in layout.regions.items():
            if rid not in active_regions:
                continue
            accent = pathway_family_color(rid, fallback_accent=region.accent) if color_mode == "pathway" else (region.accent or palette.ink2)
            lx, ly = region.label_xy
            # 地铁线路圆点 + 清晰无衬线线路名（带极淡柔光底衬，避免跨区室轨道割裂标题）
            plate_bg = "#ECEEF2" if rid in ("pyruvate_tca", "ketone", "oxphos") else "#F3F5F8"
            ax.add_patch(Circle((lx + 0.55, ly), 0.52, facecolor=accent, edgecolor=palette.surface, lw=0.45, zorder=3.85))
            ax.text(
                lx + 1.45, ly, region.label, fontsize=5.4 * font_scale, fontweight="bold",
                color=accent if color_mode == "pathway" else palette.ink, ha="left", va="center", zorder=3.85,
                bbox=dict(boxstyle="round,pad=0.14", facecolor=plate_bg, edgecolor="none", alpha=0.88),
            )
        return

    # 1. Soft zone cards (standard nature_editorial / slide_dark)
    for rid, region in layout.regions.items():
        if rid not in active_regions:
            continue
        if region.box is not None and rid != "oxphos":
            draw_zone_card(ax, region.box, region.accent, palette=palette, tint_fill=0.045, tint_edge=0.25, zorder=0.2)

    # 2. Organelle compartments (Mitochondria label placed bottom-left so BHB_trans never crosses it!)
    for cid, box in boxes.items():
        comp = layout.compartments[cid]
        draw_organelle_box(
            ax,
            box,
            comp.get("label", ""),
            kind=cid,
            palette=palette,
            label_pos="top_right" if cid == "p" else "bottom_left",
            fontsize=5.0 * font_scale,
            zorder=0.5,
        )

    # 2b. Inner OXPHOS sub-zone inside Mitochondria
    if "oxphos" in active_regions and "oxphos" in layout.regions and layout.regions["oxphos"].box is not None:
        reg = layout.regions["oxphos"]
        draw_zone_card(ax, reg.box, reg.accent, palette=palette, tint_fill=0.06, tint_edge=0.28, rounding=1.4, lw=0.6, zorder=0.62)

    # 3. Zone headers with colored accent bar
    for rid, region in layout.regions.items():
        if rid not in active_regions:
            continue
        draw_zone_header(ax, region.label_xy, region.label, region.accent, palette=palette, fontsize=5.4 * font_scale, zorder=0.9)


def _draw_edge(
    ax,
    row: pd.Series,
    geometry: dict,
    *,
    palette: SemanticPalette,
    cmap,
    layout: Layout | None = None,
    is_coupled_member: bool = False,
    font_scale: float = 1.0,
    metro: bool = False,
    color_mode: str = "semantic",
    unit_in: float = UNIT,
) -> None:
    draw_style = row["draw_style"]
    strong = row["key"] == "key"
    colored = draw_style == "flux" and pd.notna(row["value"])
    zorder = 1.6 if draw_style == "not_assessed" else (3.6 if strong else (3.1 if colored else 2.2))
    color = _edge_color(row, palette, cmap, metro=metro, color_mode=color_mode, layout=layout)
    physio = int(row["physio_sign"]) if pd.notna(row["physio_sign"]) else 0
    both = bool(row["sign_change"]) or (physio == 0 and bool(row["reversible"]))

    if metro and geometry.get("metro_paths"):
        mpaths = geometry["metro_paths"]
        hollow = draw_style == "omics_only"
        is_na = draw_style == "not_assessed"
        if is_na:
            width = 0.95
        elif hollow:
            width = 2.75 if row["key"] == "l1" else 2.20
        else:
            width = 3.70 if strong else (3.00 if colored else 2.15)
        chev_dir = -1 if (not both and physio < 0) else 1

        if is_coupled_member:
            draw_metro_coupled_halo(
                ax, [mp["drawn"] for mp in mpaths],
                color=palette.coupled_track, edge_color=palette.coupled_edge,
                width=8.2, zorder=1.9,
            )

        for mp in mpaths:
            draw_metro_track(
                ax, mp["points"], color,
                width=width,
                hollow=hollow,
                surface=palette.surface,
                chevrons=(not is_na),
                chevron_spacing=4.8,
                chevron_direction=chev_dir,
                chevron_bidirectional=both,
                unit_in=unit_in,
                trim0=mp["trim0"],
                trim1=mp["trim1"],
                corner_radius=1.15,
                halo=(not is_na),
                zorder=zorder,
            )
        if len(mpaths) > 1 and not is_na:
            ax.add_patch(Circle(
                geometry["hub"], 0.36, facecolor=color, edgecolor=palette.surface, lw=0.45, zorder=zorder + 0.18,
            ))
    else:
        style = {"flux": "-", "omics_only": (0, (3.2, 2.0)), "not_assessed": (0, (1.2, 1.8))}[draw_style]
        width = 0.62 if draw_style == "not_assessed" else (2.25 if strong else (1.60 if colored else 1.05))
        head_sub = both or physio < 0
        head_pro = both or physio >= 0
        halo = strong or colored

        if is_coupled_member:
            draw_coupled_track(
                ax, geometry["curves"],
                color=palette.coupled_track, edge_color=palette.coupled_edge,
                width=7.4, zorder=1.9,
            )

        for cv in geometry["curves"]:
            role = cv["role"]
            hs = head_sub if role in ("single", "sub") else False
            he = head_pro if role in ("single", "pro") else False
            draw_quad_arrow(
                ax, cv["p0"], cv["c"], cv["p1"], color, width, style, hs, he,
                cv["trim0"], cv["trim1"], zorder=zorder, halo=halo, surface=palette.surface,
            )

    hub = geometry["hub"]
    if bool(row["sign_change"]):
        draw_sign_reversal_badge(ax, hub[0], hub[1], radius=0.78, palette=palette, lw=0.68, zorder=4.5)
    if geometry["note"]:
        ax.text(
            hub[0], hub[1] + 0.45, geometry["note"], fontsize=NOTE_SIZE * font_scale, ha="center", va="top",
            color=palette.muted, style="italic", zorder=4.2,
            bbox=dict(boxstyle="round,pad=0.06", facecolor=palette.surface, edgecolor="none", alpha=0.84),
        )


def _draw_leader(
    ax,
    hub: tuple[float, float],
    position: tuple[float, float],
    row: pd.Series,
    group_indices: dict[str, int],
    *,
    palette: SemanticPalette,
    top_rank: int = 10,
    mode: str = "atlas",
    font_scale: float = 1.0,
    unit_in: float = UNIT,
) -> None:
    if _is_compact_label_only(row, mode):
        return
    w, _ = _strip_size(row, group_indices, top_rank=top_rank, mode=mode, font_scale=font_scale, unit_in=unit_in)
    x0, y0 = position
    end = min(((x0, y0), (x0 + w, y0)), key=lambda p: math.dist(p, hub))
    if math.dist(end, hub) > 2.0:
        ax.add_line(Line2D([hub[0], end[0]], [hub[1], end[1]], color=palette.omics_grey, lw=0.40,
                           linestyle=(0, (1.4, 1.4)), zorder=4.1))
        ax.add_patch(Circle(hub, 0.20, facecolor=palette.edge_grey, edgecolor=palette.surface, lw=0.25, zorder=4.2))


def _draw_strip(
    ax,
    position: tuple[float, float],
    row: pd.Series,
    group_indices: dict[str, int],
    *,
    palette: SemanticPalette,
    preset: str = DEFAULT_PRESET,
    top_rank: int = 10,
    mode: str = "atlas",
    strip_mode: str = "omics",
    font_scale: float = 1.0,
    unit_in: float = UNIT,
) -> None:
    x0, y0 = position
    show_rank = mode != "publication"
    info = _strip_badges(row, group_indices, top_rank=top_rank, show_rank=show_rank)
    label = info["label"]

    if _is_compact_label_only(row, mode):
        is_na = row["draw_style"] == "not_assessed"
        ax.text(
            x0 + 0.15, y0, label, fontsize=3.85 * font_scale, ha="left", va="center",
            color=palette.muted if is_na else palette.ink2, zorder=5.2,
            bbox=dict(boxstyle="round,pad=0.07", facecolor=palette.surface, edgecolor="none", alpha=0.78),
        )
        return

    w, h = _strip_size(row, group_indices, top_rank=top_rank, mode=mode, font_scale=font_scale, unit_in=unit_in)
    is_key = row["key"] == "key"
    is_l1 = row["key"] == "l1"
    is_hi = row["key"] in ("key", "l1", "top")
    z = 8.0 if is_key else (7.0 if is_l1 else 5.8)

    if is_key:
        bg, ec, lw, ls = palette.surface, palette.ink, 1.25, "-"
    elif is_l1:
        bg, ec, lw, ls = palette.surface, palette.ink2, 0.88, (0, (2.4, 1.3))
    else:
        bg, ec, lw, ls = palette.surface, palette.card_edge, 0.48, "-"

    card_left = x0 - 0.12
    card_top = y0 - h / 2
    ax.add_patch(FancyBboxPatch(
        (card_left, card_top), w + 0.12, h,
        boxstyle="round,pad=0,rounding_size=0.38",
        facecolor=bg, edgecolor=ec, lw=lw, linestyle=ls, alpha=0.98, zorder=z,
    ))

    # Top-left isolated #n rank tab (never collides with right-side tier / coupled badges)
    if info["rank"] is not None:
        draw_rank_tab(ax, card_left, card_top, info["rank"], zorder=z + 0.55, scale=0.96 * font_scale, palette=palette)

    # 4-layer omics (Flux | Protein | mRNA | Metab) or 4-axis micro-strip with crisp pure-vector header icons
    states = _extract_strip_states(row, strip_mode=strip_mode)
    header_icons = OMICS_LAYER_ICONS if strip_mode == "omics" else AXIS_LAYER_ICONS
    strip_x = x0 + 0.08
    draw_axis_strip(
        ax, strip_x, y0, states,
        cell=CELL,
        show_header_letters=is_hi,
        header_icons=header_icons,
        header_h=AXIS_HDR_H if is_hi else AXIS_BAR_H,
        zorder=z + 0.12,
        palette=palette,
        font_scale=font_scale,
    )

    # Reaction label + self-explanatory Tier badge (●●● Confirmed / ●●○ Flux / ●○○ Enz-pred) + Linked flux badge
    tx = strip_x + 4 * CELL + 0.34
    bold = is_key or is_l1
    lbl_fs = 4.35 * font_scale
    ax.text(
        tx, y0, label, fontsize=lbl_fs, ha="left", va="center", zorder=z + 0.30,
        color=palette.ink if bold else palette.ink2, fontweight="bold" if bold else "normal",
    )
    bx = tx + _text_width(label, lbl_fs, bold=bold, unit=unit_in) + 0.48
    b_scale = 1.02 * font_scale
    total_groups = max(len(group_indices), 1)
    if info["tier"]:
        bx += draw_tier_badge(
            ax, bx, y0, info["tier"], zorder=z + 0.42, scale=b_scale, preset=preset, unit_in=unit_in,
        ) + 0.34
    if info["group"] is not None:
        draw_coupled_badge(
            ax, bx, y0, info["group"], zorder=z + 0.42, scale=b_scale, palette=palette,
            total_groups=total_groups, unit_in=unit_in,
        )


def _draw_unmeasured_node_marker(
    ax,
    x: float,
    y: float,
    radius: float,
    state: str,
    *,
    palette: SemanticPalette,
    font_scale: float = 1.0,
    zorder: float = 4.75,
) -> None:
    """在未量化代谢物圆站内部绘制统一的纯矢量状态图标（替代裸 ASCII ``?`` / ``Q``）。"""
    if state == "not_detected":
        ax.add_patch(Circle(
            (x, y), radius * 0.44, facecolor="none", edgecolor=palette.muted,
            lw=0.62 * font_scale, linestyle=(0, (1.5, 1.2)), zorder=zorder,
        ))
    elif state in ("qc_blocked", "insufficient_n"):
        r = radius * 0.48
        verts = [(x, y - r), (x + r * 0.96, y + r * 0.78), (x - r * 0.96, y + r * 0.78)]
        ax.add_patch(Polygon(
            verts, closed=True, facecolor=palette.sign_change_bg,
            edgecolor=palette.sign_change_edge, lw=0.60 * font_scale, zorder=zorder,
        ))
        ax.add_line(Line2D([x, x], [y - r * 0.34, y + r * 0.16], color=palette.sign_change_ink, lw=0.55 * font_scale, zorder=zorder + 0.03))
        ax.add_patch(Circle((x, y + r * 0.44), r * 0.10, facecolor=palette.sign_change_ink, edgecolor="none", zorder=zorder + 0.03))
    elif state in ("not_measured", "no_knowledge"):
        d = radius * 0.44
        ax.add_line(Line2D([x - d, x + d], [y + d, y - d], color=palette.missing_hatch, lw=0.58, zorder=zorder))


def _draw_node(
    ax,
    layout: Layout,
    species: str,
    row: pd.Series | None,
    node_labels: dict,
    *,
    palette: SemanticPalette,
    cmap,
    font_scale: float = 1.0,
    unit_in: float = UNIT,
    active_positions: set[tuple[str, int]] | None = None,
    metro: bool = False,
) -> None:
    node = layout.nodes[species]
    for idx, (x, y) in enumerate(node.at):
        if active_positions is not None and (species, idx) not in active_positions:
            continue
        if node.kind == "boundary":
            name = node.label or _species_short(species)
            text = f"{name} (Blood)" if species.endswith("_e") else name
            fs = 3.85 * font_scale
            w = max(3.2, _text_width(text, fs, unit=unit_in) + 1.0)
            ax.add_patch(FancyBboxPatch(
                (x - w / 2, y - 0.72), w, 1.44,
                boxstyle=f"round,pad=0,rounding_size={0.72 if metro else 0.55}",
                facecolor=palette.boundary_face, edgecolor=palette.boundary_edge,
                lw=0.85 if metro else 0.60, zorder=4.2,
            ))
            ax.text(x, y, text, fontsize=fs, ha="center", va="center", color=palette.boundary_ink,
                    fontweight="bold" if metro else "medium", zorder=4.3)
            continue
        if node.kind == "store":
            text = node.label or _species_short(species)
            fs = 3.95 * font_scale
            w = max(3.8, _text_width(text, fs, unit=unit_in) + 1.1)
            ax.add_patch(FancyBboxPatch(
                (x - w / 2, y - 0.74), w, 1.48,
                boxstyle=f"round,pad=0,rounding_size={0.72 if metro else 0.50}",
                facecolor=palette.store_face, edgecolor=palette.store_edge,
                lw=0.85 if metro else 0.60, zorder=4.2,
            ))
            ax.text(x, y, text, fontsize=fs, ha="center", va="center", color=palette.store_ink,
                    fontweight="bold" if metro else "medium", zorder=4.3)
            continue

        state = str(row["display_state"]) if row is not None and pd.notna(row["display_state"]) else "not_measured"
        value = row["value"] if row is not None else np.nan
        measured = pd.notna(value)
        face = cmap(DEFAULT_NORM(value)) if measured else palette.missing_fill
        ring = row["ring"] if row is not None and isinstance(row["ring"], str) else ""
        sig = ring in ("selected", "changed")

        if metro:
            draw_metro_station(
                ax, (x, y), NODE_R * 1.05,
                fill=face,
                edge=palette.ink if sig else (palette.ink2 if measured else palette.omics_grey),
                outer_ring=palette.ink if sig else None,
                lw=1.35 if sig else (1.10 if measured else 0.85),
                surface=palette.surface,
                halo_pad=0.34,
                zorder=4.6,
            )
            if not measured:
                _draw_unmeasured_node_marker(ax, x, y, NODE_R * 1.05, state, palette=palette, font_scale=font_scale, zorder=4.76)
        else:
            ax.add_patch(Circle((x, y), NODE_R + 0.24, facecolor=palette.surface, edgecolor="none", zorder=4.4))
            if sig:
                ax.add_patch(Circle((x, y), NODE_R + 0.16, facecolor="none", edgecolor=palette.ink, lw=1.20, zorder=4.65))
            ax.add_patch(Circle(
                (x, y), NODE_R,
                facecolor=face,
                edgecolor=palette.ink if sig else (palette.ink2 if measured else palette.faint),
                lw=0.80 if sig else 0.55,
                zorder=4.55,
            ))
            if not measured:
                _draw_unmeasured_node_marker(ax, x, y, NODE_R, state, palette=palette, font_scale=font_scale, zorder=4.72)

        lx, ly, ha, va, _ = node_labels[(species, idx)]
        short_label = node.label or _species_short(species, str(row["label"]) if row is not None else species)
        ax.text(
            lx, ly, short_label, fontsize=4.6 * font_scale, ha=ha, va=va,
            color=palette.ink if (sig or measured) else palette.muted,
            fontweight="bold" if sig else "medium", zorder=5.0,
            bbox=dict(boxstyle="round,pad=0.10", facecolor=palette.surface, edgecolor="none", alpha=0.88),
        )


def _draw_aggregate(
    ax,
    layout: Layout,
    row: pd.Series,
    *,
    palette: SemanticPalette,
    cmap,
    preset: str = DEFAULT_PRESET,
    strip_mode: str = "omics",
    font_scale: float = 1.0,
    unit_in: float = UNIT,
) -> None:
    agg = layout.aggregates[row["element_id"]]
    x, y = agg.xy
    w, h = _agg_card_dims(agg)
    faint = row["draw_style"] == "not_assessed"
    is_key = row["key"] == "key"
    is_l1 = row["key"] == "l1"

    ec = palette.ink if is_key else (palette.ink2 if is_l1 else (palette.faint if faint else palette.card_edge))
    lw = 1.25 if is_key else (0.88 if is_l1 else 0.58)
    ls = "-" if is_key else ((0, (2.4, 1.3)) if is_l1 else "-")
    ax.add_patch(FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle="round,pad=0,rounding_size=0.55",
        facecolor=palette.surface, edgecolor=ec, lw=lw, linestyle=ls, zorder=3.5,
    ))

    title = str(row["label"])
    claim = str(row["flux_claim"]) if row["flux_claim"] in ("L1", "L2", "L3") else ""
    title_y = y - 1.02
    if claim:
        b_scale = 0.78 * font_scale
        badge_w = tier_badge_width(claim, scale=b_scale, unit_in=unit_in)
        avail = w - badge_w - 1.05
        t_fs = 3.85 * font_scale
        while t_fs > 2.80 and _text_width(title, t_fs, bold=(is_key or is_l1), unit=unit_in) > avail:
            t_fs -= 0.16
        ax.text(
            x - w / 2 + 0.45, title_y, title, fontsize=t_fs, ha="left", va="center",
            color=palette.muted if faint else palette.ink, fontweight="bold" if (is_key or is_l1) else "medium",
            zorder=4.0,
        )
        draw_tier_badge(
            ax, x + w / 2 - badge_w - 0.38, title_y, claim,
            zorder=4.15, scale=b_scale, preset=preset, unit_in=unit_in,
        )
    else:
        t_fs = 4.0 * font_scale if _text_width(title, 4.0 * font_scale, unit=unit_in) <= w - 0.8 else 3.45 * font_scale
        ax.text(
            x, title_y, title, fontsize=t_fs, ha="center", va="center",
            color=palette.muted if faint else palette.ink, fontweight="medium", zorder=4.0,
        )

    compact = w < 10.5
    mc = 0.72 if compact else 0.82
    m_hdr = 0.44 if (is_key or is_l1) else 0.14
    strip_w = 4 * mc
    sq_w = 0.84 if compact else 0.96
    prm_step = 1.55 if compact else 2.15
    prm_w = 2 * prm_step + sq_w + 0.65
    gap = max(0.55, min(1.15, w - strip_w - prm_w - 1.0))
    content_w = strip_w + gap + prm_w
    left_x = x - content_w / 2
    row_y = y + 0.68

    states = _extract_strip_states(row, strip_mode=strip_mode)
    header_icons = OMICS_LAYER_ICONS if strip_mode == "omics" else AXIS_LAYER_ICONS
    draw_axis_strip(
        ax, left_x, row_y, states,
        cell=mc, show_header_letters=(is_key or is_l1), header_icons=header_icons, header_h=m_hdr,
        zorder=4.0, palette=palette, font_scale=0.82 * font_scale,
    )

    prm_x0 = left_x + strip_w + gap + 0.62
    for i, (icon_kind, value_col, state_col) in enumerate((
        ("protein", "P_value", "P_state"),
        ("mrna", "R_value", "R_layer_state"),
        ("metabolite", "Mx_value", "Mx_state"),
    )):
        sx = prm_x0 + i * prm_step
        value = row[value_col]
        state = str(row[state_col]) if pd.notna(row[state_col]) else "not_measured"
        measured = pd.notna(value)
        face = cmap(DEFAULT_NORM(value)) if measured else palette.missing_fill
        draw_layer_icon(
            ax, icon_kind, sx - 0.34, row_y,
            size=0.58 if compact else 0.68, color=palette.muted, lw=0.55 * font_scale, zorder=4.15,
        )
        ax.add_patch(FancyBboxPatch(
            (sx, row_y - sq_w / 2), sq_w, sq_w,
            boxstyle="round,pad=0,rounding_size=0.14",
            facecolor=face, edgecolor=palette.card_edge, lw=0.35, zorder=4.0,
        ))
        if not measured:
            draw_status_cell(
                ax, sx, row_y - sq_w / 2, sq_w, state, 0,
                palette=palette, font_scale=0.78 * font_scale, zorder=4.12, draw_box=False,
            )


# ---------------------------------------------------------------------------
# Right Dashboard Sidebar (Atlas Mode)
# ---------------------------------------------------------------------------

_FUSION_RE = re.compile(r"^(.*?):\s*FDR\s+([0-9.eE+-]+)\s*\((.*?)\)$")


def _compact_dirs(text: str) -> str:
    rep = {
        "proteomics": "P", "transcriptomics": "R", "metabolomics": "M",
        ":up": "↑", ":down": "↓", ":mixed": "↕", "方向未定": "mixed",
    }
    out = str(text or "")
    for k, v in rep.items():
        out = out.replace(k, v)
    return out.replace(";", " ").strip() or "mixed"


def _format_fusion_items(raw: str) -> tuple[list[tuple[str, bool]], str]:
    if not raw:
        return [], ""
    sig_items: list[tuple[str, bool]] = []
    min_ns = math.inf
    n_ns = 0
    for part in [p.strip() for p in re.split(r";\s*(?![^()]*\))", raw) if p.strip()]:
        m = _FUSION_RE.match(part)
        if not m:
            sig_items.append((part, False))
            continue
        term, fdr_s, dirs = m.group(1).strip(), m.group(2), m.group(3)
        try:
            fdr = float(fdr_s)
        except ValueError:
            fdr = 1.0
        if fdr < 0.25:
            sig_items.append((f"{term}  FDR {fdr:.2g} ({_compact_dirs(dirs)})", fdr < 0.05))
        else:
            n_ns += 1
            min_ns = min(min_ns, fdr)
    ns_note = f"{n_ns} other pathway(s) FDR ≥ 0.25 (min {min_ns:.2g})" if n_ns else ""
    return sig_items, ns_note


def _sidebar(
    ax,
    regions: pd.DataFrame,
    layout: Layout,
    x0: float,
    y0: float,
    *,
    palette: SemanticPalette,
    extra_fusion: pd.DataFrame | None = None,
    signed_modules: pd.DataFrame | None = None,
) -> None:
    panel_w = SIDEBAR - 2.2
    ax.add_patch(FancyBboxPatch(
        (x0 - 0.8, y0), panel_w + 0.8, layout.height,
        boxstyle="round,pad=0,rounding_size=1.8",
        facecolor=palette.card_bg, edgecolor=palette.card_edge, lw=0.7, zorder=0.5,
    ))
    ax.text(x0 + 0.6, y0 + 1.15, "Pathway Coverage & Multi-Omics Fusion", fontsize=5.6, fontweight="bold",
            color=palette.ink, va="top")
    ax.text(x0 + 0.6, y0 + 2.85, "Coverage tracks evidence completeness; Fusion FDR is exploratory",
            fontsize=3.65, color=palette.muted, va="top")

    bar_cols = (palette.ink2, palette.edge_grey, palette.omics_grey, palette.hairline)
    bar_legend = list(zip(bar_cols, ("Ident. Flux", "Unident.", "Omics-Only", "Unassessed")))
    bx = x0 + 0.6
    for color, label in bar_legend:
        ax.add_patch(Rectangle((bx, y0 + 4.35), 1.3, 0.62, facecolor=color, edgecolor="none", zorder=1))
        ax.text(bx + 1.65, y0 + 4.66, label, fontsize=3.35, color=palette.ink2, va="center")
        bx += 1.8 + _text_width(label, 3.35) + 0.85
    ax.plot([x0 + 0.4, x0 + panel_w - 0.8], [y0 + 5.65, y0 + 5.65], color=palette.card_edge, lw=0.55, zorder=1)

    has_sm = signed_modules is not None and not signed_modules.empty
    has_ex = extra_fusion is not None and not extra_fusion.empty
    sm_h = (3.4 + min(len(signed_modules), 4) * 1.48) if has_sm else 0.0
    ex_h = (3.5 + min(len(extra_fusion), 5) * 1.75) if has_ex else 0.0
    reserved_bottom = (sm_h + 1.0 if has_sm else 0.0) + (ex_h + 1.0 if has_ex else 0.0)
    max_reg_y = y0 + layout.height - reserved_bottom - 1.2

    y = y0 + 6.45
    bar_w = panel_w - 2.0
    for row in regions.itertuples():
        if y > max_reg_y - 4.2:
            break
        reg_meta = layout.regions.get(row.region)
        accent = reg_meta.accent if reg_meta and reg_meta.accent else palette.ink2
        n_agg = int(getattr(row, "n_aggregates", 0) or 0)

        ax.add_patch(Circle((x0 + 1.0, y + 0.50), 0.42, facecolor=accent, edgecolor="none", zorder=1.5))
        reg_title = "Extra-Network Modules" if row.region == "other" else str(row.label)
        ax.text(x0 + 1.9, y, reg_title, fontsize=4.55, fontweight="bold", color=palette.ink, va="top")
        if row.kegg:
            kegg_str = str(row.kegg)
            if len(kegg_str) > 26:
                kegg_str = ", ".join(kegg_str.split(", ")[:3]) + "…"
            ax.text(x0 + panel_w - 0.9, y + 0.12, kegg_str, fontsize=3.25, color=palette.muted, va="top", ha="right")
        y += 1.52

        if int(row.n_drawn) > 0:
            total = max(int(row.n_drawn), 1)
            n_id = int(row.n_identifiable)
            n_unid = max(0, int(row.n_flux) - n_id)
            n_om = int(row.n_omics_only)
            n_na = max(0, total - n_id - n_unid - n_om)
            cur_x = x0 + 0.6
            for count, color in zip((n_id, n_unid, n_om, n_na), bar_cols):
                if count > 0:
                    seg_w = bar_w * count / total
                    ax.add_patch(Rectangle((cur_x, y), seg_w, 0.58, facecolor=color, edgecolor="none", zorder=1.2))
                    cur_x += seg_w
            y += 0.92
            stat = f"Flux {row.n_flux}/{row.n_drawn} ({row.n_identifiable} ident.) · Omics {row.n_omics_only}"
            if n_agg > 0:
                stat += f" · +{n_agg} mod."
        else:
            ax.add_patch(Rectangle((x0 + 0.6, y), bar_w, 0.58, facecolor=palette.omics_grey, edgecolor="none", zorder=1.2))
            if int(row.n_L1) > 0 and n_agg > 0:
                ax.add_patch(Rectangle((x0 + 0.6, y), bar_w * min(1.0, int(row.n_L1) / n_agg), 0.58,
                                       facecolor=accent, edgecolor="none", alpha=0.65, zorder=1.3))
            y += 0.92
            stat = f"Omics-only modules: {n_agg}"

        ax.text(x0 + 0.6, y, stat, fontsize=3.6, color=palette.ink2, va="top")
        claim_text = f"L2+ {row.n_L2} · L1 {row.n_L1}"
        ax.text(x0 + panel_w - 0.9, y, claim_text, fontsize=3.6, va="top", ha="right",
                color=palette.ink if (row.n_L2 or row.n_L1) else palette.muted,
                fontweight="bold" if row.n_L2 else "normal")
        y += 1.30

        sig_items, ns_note = _format_fusion_items(row.fusion)
        for text, strong_sig in sig_items[:2]:
            for li, wrapped in enumerate(textwrap.wrap(text, 70)[:1]):
                if len(text) > 70:
                    wrapped = wrapped.rstrip() + "…"
                bullet = ("● " if strong_sig else "○ ") if li == 0 else "  "
                ax.text(x0 + 0.9, y, bullet + wrapped, fontsize=3.45,
                        color=palette.up if strong_sig else palette.ink,
                        fontweight="bold" if strong_sig else "normal", va="top")
                y += 1.12
        if not sig_items and ns_note:
            ax.text(x0 + 0.6, y, ns_note, fontsize=3.3, color=palette.muted, va="top")
            y += 1.05
        elif not sig_items and not ns_note:
            ax.text(x0 + 0.6, y, "Fusion: no matched multi-omics term", fontsize=3.3, color=palette.muted, va="top")
            y += 1.05

        y += 0.55
        ax.plot([x0 + 0.6, x0 + panel_w - 0.9], [y - 0.22, y - 0.22], color=palette.hairline, lw=0.45, zorder=1)
        y += 0.28

    cur_bot = y0 + layout.height - 0.9
    if has_ex:
        ex_top = cur_bot - ex_h
        ax.add_patch(FancyBboxPatch(
            (x0 + 0.2, ex_top), panel_w - 1.0, ex_h,
            boxstyle="round,pad=0,rounding_size=1.1",
            facecolor=tint_color(palette.up, 0.04, base=palette.surface),
            edgecolor=tint_color(palette.up, 0.22, base=palette.surface), lw=0.6, zorder=1.0,
        ))
        ax.text(x0 + 0.9, ex_top + 0.9, "Top Extra-Network Pathways (Fusion, FDR < 0.25)",
                fontsize=4.3, fontweight="bold", color=palette.ink, va="top", zorder=1.2)
        ey = ex_top + 2.55
        for erow in extra_fusion.head(5).itertuples():
            fdr = float(erow.fused_fdr)
            strong = fdr < 0.05
            term = str(erow.pathway_term)
            if len(term) > 48:
                term = term[:47].rstrip() + "…"
            dirs = _compact_dirs(erow.omics_directions)
            bullet = "● " if strong else "○ "
            ax.text(x0 + 0.9, ey, f"{bullet}{term}", fontsize=3.55, va="top",
                    color=palette.up if strong else palette.ink, fontweight="bold" if strong else "normal", zorder=1.2)
            ax.text(x0 + panel_w - 1.3, ey, f"FDR {fdr:.2g} ({dirs})", fontsize=3.4, va="top", ha="right",
                    color=palette.up if strong else palette.ink2, zorder=1.2)
            ey += 1.72
        cur_bot = ex_top - 0.9

    if has_sm:
        sm_top = cur_bot - sm_h
        ax.add_patch(FancyBboxPatch(
            (x0 + 0.2, sm_top), panel_w - 1.0, sm_h,
            boxstyle="round,pad=0,rounding_size=1.1",
            facecolor=tint_color(palette.down, 0.045, base=palette.surface),
            edgecolor=tint_color(palette.down, 0.24, base=palette.surface), lw=0.6, zorder=1.0,
        ))
        ax.text(x0 + 0.9, sm_top + 0.9, "Signed Module Rotation Tests (G16)",
                fontsize=4.3, fontweight="bold", color=palette.ink, va="top", zorder=1.2)
        sy = sm_top + 2.55
        layer_tag = {"proteomics": "Protein", "transcriptomics": "RNA", "metabolomics": "Metab"}
        for srow in signed_modules.head(4).itertuples():
            padj = float(srow.padj) if pd.notna(srow.padj) else 1.0
            stat = float(srow.statistic) if pd.notna(srow.statistic) else 0.0
            sig = padj < 0.05
            mod_name = str(srow.name).replace("_", " ")
            lay = layer_tag.get(str(srow.layer), str(srow.layer))
            dir_s = str(srow.direction).replace("_", " ")
            bullet = "● " if sig else "○ "
            ax.text(x0 + 0.9, sy, f"{bullet}{mod_name} [{lay}]", fontsize=3.5, va="top",
                    color=palette.up if sig else palette.ink, fontweight="bold" if sig else "normal", zorder=1.2)
            ax.text(x0 + panel_w - 1.3, sy, f"z={stat:+.2f}, padj={padj:.2g} ({dir_s})", fontsize=3.35,
                    va="top", ha="right", color=palette.up if sig else palette.ink2, zorder=1.2)
            sy += 1.48


# ---------------------------------------------------------------------------
# Bottom Legends (4-Column Atlas Legend & Compact Publication Legend)
# ---------------------------------------------------------------------------

def _legend_atlas(
    ax,
    y0: float,
    total_width: float,
    caption: list[str],
    *,
    palette: SemanticPalette,
    cmap,
    preset: str = DEFAULT_PRESET,
    metro: bool = False,
    color_mode: str = "semantic",
    strip_mode: str = "omics",
    unit_in: float = UNIT,
) -> None:
    panel_w = total_width - 2.0
    legend_h = 13.8
    ax.add_patch(FancyBboxPatch(
        (1.0, y0), panel_w, legend_h,
        boxstyle="round,pad=0,rounding_size=1.4",
        facecolor=palette.card_bg, edgecolor=palette.card_edge, lw=0.7, zorder=1,
    ))

    c1 = 2.6
    c2 = 1.0 + panel_w * 0.175
    c3 = 1.0 + panel_w * 0.410
    c4 = 1.0 + panel_w * 0.775
    for cx in (c2 - 1.5, c3 - 1.5, c4 - 1.5):
        ax.plot([cx, cx], [y0 + 1.0, y0 + legend_h - 1.0], color=palette.card_edge, lw=0.55, zorder=1.5)

    # Col 1: Color Scale
    ax.text(c1, y0 + 1.0, "1. Color Scale (log2FC)", fontsize=4.7, fontweight="bold", color=palette.ink, va="top")
    bar_len = (c2 - c1) - 4.6
    xs = np.linspace(0, bar_len, 61)
    for i, x in enumerate(xs[:-1]):
        value = -DEFAULT_VMAX + 2 * DEFAULT_VMAX * i / (len(xs) - 2)
        ax.add_patch(Rectangle((c1 + x, y0 + 3.2), xs[1] - xs[0] + 0.05, 1.35,
                               facecolor=cmap(DEFAULT_NORM(value)), edgecolor="none", zorder=2))
    ax.add_patch(Rectangle((c1, y0 + 3.2), bar_len, 1.35, facecolor="none", edgecolor=palette.card_edge, lw=0.45, zorder=2.2))
    for v, lab in ((-2, "-2 (Down)"), (0, "0"), (2, "+2 (Up)")):
        ax.text(c1 + (v + DEFAULT_VMAX) / (2 * DEFAULT_VMAX) * bar_len, y0 + 4.95, lab,
                fontsize=3.75, ha="center", va="top", color=palette.ink2)
    if metro and color_mode == "pathway":
        col1_note = (
            "Node fill: tissue metabolite or lipid-class total\n"
            "Track color: Morandi pathway family\n"
            "(check 4-icon strip for Up/Down direction)"
        )
    else:
        col1_note = (
            "Node fill: tissue metabolite or lipid-class total\n"
            "Edge color: MFA flux effect along net flow\n"
            "(directional |log2FC| < 0.8 floored at ±0.8)"
        )
    ax.text(c1, y0 + 7.2, col1_note, fontsize=3.7, color=palette.ink2, va="top", linespacing=1.45)

    # Col 2: Edges & Linked Flux Conduit
    ax.text(c2, y0 + 1.0, "2. Reaction Tracks & Net Flow" if metro else "2. Reaction Edges & Net Flow",
            fontsize=4.7, fontweight="bold", color=palette.ink, va="top")
    if metro:
        metro_samples = [
            (y0 + 3.05, palette.edge_grey, 2.2, False, True, False, "MFA flux fitted (solid rail; grey = undecided)"),
            (y0 + 4.88, palette.up, 3.2, False, True, False,
             "MFA flux directional (solid rail; red up / blue down)" if color_mode == "semantic" else "MFA flux active (solid rail; color = pathway family)"),
            (y0 + 6.71, "#5B6C82", 2.6, True, True, False, "Omics-only / Enzyme-predicted (hollow double-rail)"),
            (y0 + 8.54, palette.edge_grey, 2.2, False, True, True, "Reversible reaction (no single net dir.; opposing < >)"),
            (y0 + 10.37, palette.faint, 1.0, False, False, False, "Not assessed skeleton (no MFA & no mapped members)"),
        ]
        for ey, col, lw, hollow, chev, bidir, desc in metro_samples:
            draw_metro_track(
                ax, [(c2, ey), (c2 + 6.4, ey)], col,
                width=lw, hollow=hollow, surface=palette.card_bg,
                chevrons=chev, chevron_spacing=10.0, chevron_bidirectional=bidir,
                unit_in=unit_in, corner_radius=0.0, halo=False, zorder=2.5,
            )
            ax.text(c2 + 7.4, ey, desc, fontsize=3.6, color=palette.ink2, va="center")
        draw_metro_coupled_halo(
            ax, [[(c2, y0 + 12.22), (c2 + 6.4, y0 + 12.22)]],
            color=palette.coupled_track, edge_color=palette.coupled_edge, width=7.0, zorder=2.2,
        )
        draw_metro_track(
            ax, [(c2, y0 + 12.22), (c2 + 6.4, y0 + 12.22)], palette.up,
            width=2.8, hollow=False, surface=palette.card_bg,
            chevrons=True, chevron_spacing=10.0, unit_in=unit_in,
            corner_radius=0.0, halo=False, zorder=2.6,
        )
        ax.text(c2 + 7.4, y0 + 12.22, "Indigo conduit: stoichiometrically linked steps (shared MFA test)",
                fontsize=3.55, color=palette.coupled_ink, va="center")
    else:
        edge_samples = [
            (y0 + 3.05, "-", palette.edge_grey, 1.2, False, True, "MFA flux fitted (grey = undecided / unident.)"),
            (y0 + 4.88, "-", palette.up, 1.9, False, True, "MFA flux directional (red up / blue down)"),
            (y0 + 6.71, (0, (3.2, 2.0)), palette.omics_grey, 1.05, False, True, "Omics-only / Enzyme-predicted edge (no MFA yet)"),
            (y0 + 8.54, "-", palette.edge_grey, 1.2, True, True, "Reversible reaction (no single net direction)"),
            (y0 + 10.37, (0, (1.2, 1.8)), palette.faint, 0.72, False, True, "Not assessed (no MFA & no mapped members)"),
        ]
        for ey, st, col, lw, hs, he, desc in edge_samples:
            draw_quad_arrow(ax, (c2, ey), (c2 + 3.2, ey), (c2 + 6.4, ey), col, lw, st, hs, he, 0.0, 0.0,
                            zorder=2.5, surface=palette.card_bg)
            ax.text(c2 + 7.4, ey, desc, fontsize=3.6, color=palette.ink2, va="center")
        draw_coupled_track(ax, [{"p0": (c2, y0 + 12.22), "c": (c2 + 3.2, y0 + 12.22), "p1": (c2 + 6.4, y0 + 12.22), "trim0": 0, "trim1": 0}],
                           color=palette.coupled_track, edge_color=palette.coupled_edge, width=6.4, zorder=2.2)
        draw_quad_arrow(ax, (c2, y0 + 12.22), (c2 + 3.2, y0 + 12.22), (c2 + 6.4, y0 + 12.22), palette.up, 1.6, "-", False, True, 0.0, 0.0,
                        zorder=2.6, surface=palette.card_bg)
        ax.text(c2 + 7.4, y0 + 12.22, "Indigo conduit: stoichiometrically linked steps (shared MFA test)",
                fontsize=3.55, color=palette.coupled_ink, va="center")

    # Col 3: 4-Layer Icon Strip, Real Visual State Swatches & Self-Explanatory Badges
    layer_defs = (
        (("flux", "Flux"), ("protein", "Protein"), ("mrna", "mRNA"), ("metabolite", "Metab"))
        if strip_mode == "omics" else
        (("flux", "Flux"), ("protein", "Enzyme"), ("regulation", "Regul."), ("metabolite", "Pool"))
    )
    ax.text(
        c3, y0 + 1.0,
        "3. Multi-Omics Strip Icons, Cell States & Evidence Badges" if strip_mode == "omics" else "3. Four-Axis Strip Icons, Cell States & Evidence Badges",
        fontsize=4.7, fontweight="bold", color=palette.ink, va="top",
    )
    # Row 3a: 4 enlarged layer header icons + directional colored cell swatches
    ry1 = y0 + 3.75
    ix = c3 + 0.3
    for ikind, ilbl in layer_defs:
        draw_layer_icon(ax, ikind, ix + 0.55, ry1, size=1.12, color=palette.ink, lw=0.75, zorder=2.5)
        ax.text(ix + 1.45, ry1, ilbl, fontsize=3.65, fontweight="bold", color=palette.ink2, va="center")
        ix += 1.55 + _text_width(ilbl, 3.65, bold=True, unit=unit_in) + 1.45
    # Directional swatches right beside layer icons
    ix += 0.8
    sw = 1.25
    for st_pair, desc in (
        ((("changed", 1), ("changed", -1)), "Sig. (FDR<0.05)"),
        ((("suggestive", 1), ("suggestive", -1)), "Suggestive (p<0.05)"),
        ((("mixed", 0),), "Mixed"),
    ):
        for st, dr in st_pair:
            draw_status_cell(ax, ix, ry1 - sw / 2, sw, st, dr, palette=palette, font_scale=1.1, zorder=2.5, draw_box=True)
            ix += sw + 0.28
        ax.text(ix + 0.25, ry1, desc, fontsize=3.55, color=palette.ink2, va="center")
        ix += 0.35 + _text_width(desc, 3.55, unit=unit_in) + 1.35

    # Row 3b: Non-directional / Missing / Alert icon swatches (every icon drawn directly!)
    ry2 = y0 + 6.40
    ix = c3 + 0.3
    for st, desc in (
        ("unchanged", "Unchanged"),
        ("not_measured", "Not measured"),
        ("not_identifiable", "Unidentifiable flux"),
        ("not_detected", "Not detected"),
        ("qc_blocked", "QC / Low-n"),
    ):
        draw_status_cell(ax, ix, ry2 - sw / 2, sw, st, 0, palette=palette, font_scale=1.1, zorder=2.5, draw_box=True)
        ax.text(ix + sw + 0.45, ry2, desc, fontsize=3.55, color=palette.ink2, va="center")
        ix += sw + 0.55 + _text_width(desc, 3.55, unit=unit_in) + 1.55
    draw_sign_reversal_badge(ax, ix + 0.72, ry2, radius=0.72, palette=palette, lw=0.65, zorder=2.5)
    ax.text(ix + 1.85, ry2, "Net flux reversal (group sign flip)", fontsize=3.55, color=palette.sign_change_ink,
            fontweight="medium", va="center")

    # Row 3c: Self-explanatory 3-dot tier badges (●●● Confirmed / ●●○ Flux / ●○○ Enz-pred)
    by1 = y0 + 9.15
    bx = c3 + 0.3
    for lvl, txt in (
        ("L3", "Multi-tracer confirmed"),
        ("L2", "MFA flux supported"),
        ("L1", "Enzyme-predicted (dashed; verify decoy)"),
    ):
        bx += draw_tier_badge(ax, bx, by1, lvl, zorder=2.5, scale=1.08, preset=preset, unit_in=unit_in) + 0.65
        ax.text(bx, by1, txt, fontsize=3.6, color=palette.ink2, va="center")
        bx += _text_width(txt, 3.6, unit=unit_in) + 2.0

    # Row 3d: Corner rank tab + Linked flux badge with clear statistical meaning
    by2 = y0 + 11.70
    bx = c3 + 0.3
    draw_rank_tab(ax, bx + 0.1, by2 + 0.25, 1, zorder=2.5, scale=1.12, palette=palette)
    ax.text(bx + 2.55, by2, "Corner #n: Top-N rank", fontsize=3.6, color=palette.ink2, va="center")
    bx2 = bx + 2.55 + _text_width("Corner #n: Top-N rank", 3.6, unit=unit_in) + 2.8
    bx2 += draw_coupled_badge(ax, bx2, by2, 1, zorder=2.5, scale=1.06, palette=palette, total_groups=1, unit_in=unit_in) + 0.7
    ax.text(bx2, by2, "Stoichiometrically linked reactions sharing one MFA test (non-independent)",
            fontsize=3.55, color=palette.ink2, va="center")

    # Col 4: Nodes, Organelles & Modules
    ax.text(c4, y0 + 1.0, "4. Nodes, Organelles & Modules", fontsize=4.7, fontweight="bold", color=palette.ink, va="top")
    ax.add_patch(Circle((c4 + 1.4, y0 + 3.55), 0.92, facecolor=cmap(DEFAULT_NORM(0.8)), edgecolor=palette.ink2, lw=0.55, zorder=2.2))
    ax.add_patch(Circle((c4 + 4.3, y0 + 3.55), 1.08, facecolor="none", edgecolor=palette.ink, lw=1.15, zorder=2.3))
    ax.add_patch(Circle((c4 + 4.3, y0 + 3.55), 0.92, facecolor=cmap(DEFAULT_NORM(1.7)), edgecolor=palette.ink, lw=0.75, zorder=2.2))
    ax.add_patch(Circle((c4 + 7.2, y0 + 3.55), 0.92, facecolor=palette.missing_fill, edgecolor=palette.faint, lw=0.55, zorder=2.2))
    _draw_unmeasured_node_marker(ax, c4 + 7.2, y0 + 3.55, 0.92, "not_detected", palette=palette, zorder=2.4)
    ax.text(c4 + 9.0, y0 + 3.55, "Metabolite (double ring = FDR sig.; ◌ = not detected)", fontsize=3.55, color=palette.ink2, va="center")

    ax.add_patch(FancyBboxPatch((c4 + 0.4, y0 + 5.55), 5.6, 1.35, boxstyle="round,pad=0,rounding_size=0.55",
                                facecolor=palette.boundary_face, edgecolor=palette.boundary_edge, lw=0.55, zorder=2.2))
    ax.text(c4 + 3.2, y0 + 6.22, "Glc (Blood)", fontsize=3.35, ha="center", va="center", color=palette.boundary_ink, zorder=2.3)
    ax.add_patch(FancyBboxPatch((c4 + 6.6, y0 + 5.55), 4.0, 1.35, boxstyle="round,pad=0,rounding_size=0.50",
                                facecolor=palette.store_face, edgecolor=palette.store_edge, lw=0.55, zorder=2.2))
    ax.text(c4 + 8.6, y0 + 6.22, "Store", fontsize=3.35, ha="center", va="center", color=palette.store_ink, zorder=2.3)
    ax.text(c4 + 11.4, y0 + 6.22, "Plasma boundary / storage pool", fontsize=3.55, color=palette.ink2, va="center")

    ax.add_patch(FancyBboxPatch((c4 + 0.4, y0 + 7.65), 4.8, 1.35, boxstyle="round,pad=0,rounding_size=0.45",
                                facecolor=palette.mito_face, edgecolor=palette.mito_edge, lw=0.9, zorder=2.2))
    ax.text(c4 + 2.8, y0 + 8.32, "Mito", fontsize=3.3, ha="center", va="center", color=palette.mito_ink, fontweight="bold", zorder=2.3)
    ax.add_patch(FancyBboxPatch((c4 + 5.8, y0 + 7.65), 4.8, 1.35, boxstyle="round,pad=0,rounding_size=0.45",
                                facecolor=palette.perox_face, edgecolor=palette.perox_edge, lw=0.9, zorder=2.2))
    ax.text(c4 + 8.2, y0 + 8.32, "Perox", fontsize=3.3, ha="center", va="center", color=palette.perox_ink, fontweight="bold", zorder=2.3)
    ax.text(c4 + 11.4, y0 + 8.32, "Mitochondrial matrix / Peroxisome", fontsize=3.55, color=palette.ink2, va="center")

    ax.add_patch(FancyBboxPatch((c4 + 0.4, y0 + 9.85), 10.6, 2.55, boxstyle="round,pad=0,rounding_size=0.45",
                                facecolor=palette.surface, edgecolor=palette.card_edge, lw=0.55, zorder=2.2))
    ax.text(c4 + 5.7, y0 + 10.42, "Omics Module", fontsize=3.3, ha="center", va="center", color=palette.ink, zorder=2.3)
    draw_axis_strip(ax, c4 + 0.85, y0 + 11.50, [("not_measured", 0), ("changed", 1), ("unchanged", 0), ("not_measured", 0)],
                    cell=0.66, show_header_letters=False, header_h=0.14, zorder=2.3, palette=palette, font_scale=0.8)
    for i, (ikind, col) in enumerate((("protein", cmap(DEFAULT_NORM(1.1))), ("mrna", cmap(DEFAULT_NORM(-0.9))), ("metabolite", palette.missing_fill))):
        cx = c4 + 5.05 + i * 2.0
        draw_layer_icon(ax, ikind, cx - 0.34, y0 + 11.50, size=0.58, color=palette.muted, lw=0.52, zorder=2.4)
        ax.add_patch(Rectangle((cx, y0 + 11.12), 0.76, 0.76, facecolor=col, edgecolor=palette.card_edge, lw=0.3, zorder=2.3))
    ax.text(c4 + 11.8, y0 + 11.10, "Isolated module (4-icon strip + Prot / mRNA / Metab)", fontsize=3.55, color=palette.ink2, va="center")

    # Compact Provenance & Decoy Calibration Caption
    y = y0 + legend_h + 1.05
    for line in caption:
        for wrapped in textwrap.wrap(line, 275):
            ax.text(1.4, y, wrapped, fontsize=3.95, va="top", color=palette.ink2)
            y += 1.48


def _flow_rows(widths: Sequence[float], max_w: float, gap: float) -> tuple[list[tuple[int, float]], int]:
    """将一组定宽图例项按 ``max_w`` 自动换行，返回 ``([(row_idx, x_offset), ...], n_rows)``。"""
    placed: list[tuple[int, float]] = []
    row = 0
    x = 0.0
    for w in widths:
        if x > 0.0 and x + w > max_w:
            row += 1
            x = 0.0
        placed.append((row, x))
        x += w + gap
    return placed, (row + 1 if widths else 0)


def _header_pills(view: pd.DataFrame, palette: SemanticPalette, preset: str) -> list[tuple[str, str, str, str]]:
    rxns = view[view["element_type"] == "reaction"]
    n_flux = int((rxns["draw_style"] == "flux").sum())
    n_l2 = int(view["flux_claim"].isin(["L2", "L3"]).sum())
    n_l1 = int((view["flux_claim"] == "L1").sum())
    n_agg = int((view["element_type"] == "aggregate").sum())
    tok = tier_tokens(preset)
    pills = [
        (f"Flux: {n_l2}", tok["L2"]["bg"], tok["L2"]["fg"], "L2"),
        (f"Enz-pred: {n_l1}", tok["L1"]["bg"], tok["L1"]["fg"], "L1"),
    ]
    edge_txt = f"MFA Flux: {n_flux}/{len(rxns)}" + (f" | Mod: {n_agg}" if n_agg else "")
    pills.append((edge_txt, palette.card_bg, palette.ink2, ""))
    return pills


def _header_geometry(
    title: str,
    view: pd.DataFrame,
    avail_w: float,
    *,
    palette: SemanticPalette,
    preset: str,
    font_scale: float,
    unit_in: float,
) -> tuple[float, bool, float, float]:
    """计算顶栏高度、是否换行堆叠及适配字号 ``(header_h, stacked, title_fs, pill_fs)``。"""
    pills = _header_pills(view, palette, preset)
    title_fs = 8.4 * font_scale
    pill_fs = 4.15 * font_scale

    def _pills_total_w(pfs: float) -> float:
        total = 0.0
        for text, _, _, icon in pills:
            icon_pad = 2.05 * (pfs / 4.15) if icon else 0.0
            total += _text_width(text, pfs, bold=True, unit=unit_in) + 1.8 + icon_pad + 1.0
        return max(0.0, total - 1.0)

    tw = _text_width(title, title_fs, bold=True, unit=unit_in)
    pw = _pills_total_w(pill_fs)
    if tw + pw + 4.0 <= avail_w:
        return HEADER, False, title_fs, pill_fs

    # 窄/中幅子图自动切换为清爽两行顶栏，避免主标题与统计药丸碰撞
    while title_fs > 5.2 and _text_width(title, title_fs, bold=True, unit=unit_in) > avail_w - 2.0:
        title_fs -= 0.35
    while pill_fs > 3.2 and _pills_total_w(pill_fs) > avail_w - 2.0:
        pill_fs -= 0.20
    return 9.2, True, title_fs, pill_fs


def _pub_legend_plan(
    vw: float,
    caption: list[str],
    folded_regions: list[str],
    panel_note: str,
    *,
    font_scale: float,
    unit_in: float,
    metro: bool,
    color_mode: str,
    strip_mode: str = "omics",
) -> dict:
    """预先计算出版/子图紧凑图例的精确换行与总高度（全部状态与图标均以真实图元呈现）。"""
    fs = 3.75 * font_scale
    inner_w = max(24.0, vw - 4.0)
    bar_len = min(20.0, max(11.0, vw * 0.22))
    sw = 1.28

    # Row-group A: Colorbar + Track/Arrow types + Reversible + Linked conduit + Sign reversal + FDR-sig metabolite
    cb_w = max(bar_len, _text_width("log2FC (−2 → +2)", fs * 0.90, bold=True, unit=unit_in)) + 1.0
    if metro:
        flux_lbl = "MFA Flux (solid)" if color_mode == "semantic" else "MFA Flux (pathway)"
        a_items = [
            ("colorbar", "", cb_w),
            ("metro_flux", flux_lbl, 5.8 + _text_width(flux_lbl, fs, unit=unit_in)),
            ("metro_omics", "Omics-only (hollow)", 5.8 + _text_width("Omics-only (hollow)", fs, unit=unit_in)),
            ("metro_rev", "Reversible (no net dir.)", 5.8 + _text_width("Reversible (no net dir.)", fs, unit=unit_in)),
            ("metro_coupled", "Linked flux conduit", 5.8 + _text_width("Linked flux conduit", fs, unit=unit_in)),
            ("sign_reversal", "Net flux reversal", 2.3 + _text_width("Net flux reversal", fs, unit=unit_in)),
            ("metab", "FDR-sig. metab.", 2.6 + _text_width("FDR-sig. metab.", fs, unit=unit_in)),
        ]
    else:
        a_items = [
            ("colorbar", "", cb_w),
            ("arrow_flux", "MFA Flux Δ", 5.6 + _text_width("MFA Flux Δ", fs, unit=unit_in)),
            ("arrow_rev", "Reversible (no net dir.)", 5.6 + _text_width("Reversible (no net dir.)", fs, unit=unit_in)),
            ("arrow_coupled", "Linked flux conduit", 5.6 + _text_width("Linked flux conduit", fs, unit=unit_in)),
            ("sign_reversal", "Net flux reversal", 2.3 + _text_width("Net flux reversal", fs, unit=unit_in)),
            ("metab", "FDR-sig. metab.", 2.6 + _text_width("FDR-sig. metab.", fs, unit=unit_in)),
        ]
    a_pos, a_rows = _flow_rows([w for _, _, w in a_items], inner_w, gap=2.2)

    # Row-group B: 4 Omics Layer Icons + Real Colored Cell State Swatches (Solid ▲▼, Hollow △▽, Mixed, —, /, ⊘, ◌, ⚠)
    layer_defs = (
        (("flux", "Flux"), ("protein", "Protein"), ("mrna", "mRNA"), ("metabolite", "Metab"))
        if strip_mode == "omics" else
        (("flux", "Flux"), ("protein", "Enzyme"), ("regulation", "Regul."), ("metabolite", "Pool"))
    )
    b_items: list[tuple[str, str, float]] = []
    for ikind, ilbl in layer_defs:
        b_items.append((f"icon:{ikind}", ilbl, 1.55 + _text_width(ilbl, fs, bold=True, unit=unit_in)))
    b_items.extend([
        ("cell_sig", "Sig. (FDR<0.05)", 2 * sw + 0.65 + _text_width("Sig. (FDR<0.05)", fs, unit=unit_in)),
        ("cell_sug", "Suggestive (p<0.05)", 2 * sw + 0.65 + _text_width("Suggestive (p<0.05)", fs, unit=unit_in)),
        ("cell_mixed", "Mixed", sw + 0.50 + _text_width("Mixed", fs, unit=unit_in)),
        ("cell_unchanged", "Unchanged", sw + 0.50 + _text_width("Unchanged", fs, unit=unit_in)),
        ("cell_na", "Not measured", sw + 0.50 + _text_width("Not measured", fs, unit=unit_in)),
        ("cell_unident", "Unident. flux", sw + 0.50 + _text_width("Unident. flux", fs, unit=unit_in)),
        ("cell_nodet", "Not detected", sw + 0.50 + _text_width("Not detected", fs, unit=unit_in)),
        ("cell_qc", "QC / Low-n", sw + 0.50 + _text_width("QC / Low-n", fs, unit=unit_in)),
    ])
    b_pos, b_rows = _flow_rows([w for _, _, w in b_items], inner_w, gap=1.8)

    # Row-group C: Self-explanatory Tier Badges + Linked Flux Badge
    c_scale = 1.02 * font_scale
    linked_desc = "Shared MFA test (non-independent)"
    c_items = [
        ("L3", "Multi-tracer", tier_badge_width("L3", scale=c_scale, unit_in=unit_in) + 0.55 + _text_width("Multi-tracer", fs, unit=unit_in)),
        ("L2", "MFA supported", tier_badge_width("L2", scale=c_scale, unit_in=unit_in) + 0.55 + _text_width("MFA supported", fs, unit=unit_in)),
        ("L1", "Enzyme only", tier_badge_width("L1", scale=c_scale, unit_in=unit_in) + 0.55 + _text_width("Enzyme only", fs, unit=unit_in)),
        ("LINK", linked_desc, coupled_badge_width(1, scale=c_scale, total_groups=1, unit_in=unit_in) + 0.55 + _text_width(linked_desc, fs, unit=unit_in)),
    ]
    c_pos, c_rows = _flow_rows([w for _, _, w in c_items], inner_w, gap=2.0)

    # Disclosure note (submap scope or folded regions)
    char_w = max(em_units(fs * 0.88, unit_in) * 0.56, 0.18)
    wrap_cols = max(36, int(inner_w / char_w))
    note_text = panel_note or (
        ("Inactive / unassessed regions folded from print view: " + " · ".join(folded_regions))
        if folded_regions else ""
    )
    note_lines = textwrap.wrap(note_text, wrap_cols)[:2] if note_text else []

    cap_wrap_cols = max(40, int((vw - 2.0) / max(em_units(fs * 0.86, unit_in) * 0.56, 0.18)))
    cap_lines: list[str] = []
    for line in caption[:2]:
        cap_lines.extend(textwrap.wrap(line, cap_wrap_cols))

    row_a = 3.20
    row_b = 2.65
    row_c = 2.85
    legend_h = (
        1.5
        + a_rows * row_a + 0.45
        + b_rows * row_b + 0.50
        + c_rows * row_c
        + (0.45 + len(note_lines) * 1.80 if note_lines else 0.0)
        + 0.9
    )
    footer_h = 1.6 + legend_h + 1.0 + len(cap_lines) * 1.65 + 0.8
    return {
        "fs": fs, "bar_len": bar_len, "sw": sw,
        "row_a": row_a, "row_b": row_b, "row_c": row_c,
        "a_items": a_items, "a_pos": a_pos, "a_rows": a_rows,
        "b_items": b_items, "b_pos": b_pos, "b_rows": b_rows,
        "c_items": c_items, "c_pos": c_pos, "c_rows": c_rows,
        "note_lines": note_lines, "cap_lines": cap_lines,
        "legend_h": legend_h, "footer_h": footer_h,
    }


def _legend_publication(
    ax,
    vx0: float,
    y0: float,
    vw: float,
    plan: dict,
    *,
    palette: SemanticPalette,
    cmap,
    preset: str = DEFAULT_PRESET,
    font_scale: float = 1.25,
    unit_in: float = UNIT,
) -> None:
    """响应式出版/子图紧凑图例（带真实组学图标、实心/空心色格、状态图标与自解释徽章）。"""
    legend_h = plan["legend_h"]
    fs = plan["fs"]
    bar_len = plan["bar_len"]
    sw = plan["sw"]
    row_a, row_b, row_c = plan["row_a"], plan["row_b"], plan["row_c"]

    ax.add_patch(FancyBboxPatch(
        (vx0 + 0.5, y0), vw - 1.0, legend_h,
        boxstyle="round,pad=0,rounding_size=1.2",
        facecolor=palette.card_bg, edgecolor=palette.card_edge, lw=0.7, zorder=1,
    ))

    # Row-group A: Colorbar + Tracks + Metabolite station
    y_a0 = y0 + 2.45
    for (kind, lbl, _), (r_idx, x_off) in zip(plan["a_items"], plan["a_pos"]):
        x = vx0 + 2.0 + x_off
        y = y_a0 + r_idx * row_a
        if kind == "colorbar":
            xs = np.linspace(0, bar_len, 41)
            for i, xi in enumerate(xs[:-1]):
                val = -DEFAULT_VMAX + 2 * DEFAULT_VMAX * i / (len(xs) - 2)
                ax.add_patch(Rectangle((x + xi, y - 0.55), xs[1] - xs[0] + 0.05, 1.15,
                                       facecolor=cmap(DEFAULT_NORM(val)), edgecolor="none", zorder=2))
            ax.add_patch(Rectangle((x, y - 0.55), bar_len, 1.15, facecolor="none", edgecolor=palette.card_edge, lw=0.45, zorder=2.2))
            ax.text(x, y - 0.85, "log2FC (−2 → +2)", fontsize=fs * 0.88, fontweight="bold", color=palette.ink, va="bottom")
        elif kind == "metro_flux":
            draw_metro_track(
                ax, [(x, y), (x + 5.0, y)], palette.up,
                width=3.2, hollow=False, surface=palette.card_bg,
                chevrons=True, chevron_spacing=9.0, unit_in=unit_in,
                corner_radius=0.0, halo=False, zorder=2.5,
            )
            ax.text(x + 5.7, y, lbl, fontsize=fs, color=palette.ink2, va="center")
        elif kind == "metro_omics":
            draw_metro_track(
                ax, [(x, y), (x + 5.0, y)], "#5B6C82",
                width=2.6, hollow=True, surface=palette.card_bg,
                chevrons=True, chevron_spacing=9.0, unit_in=unit_in,
                corner_radius=0.0, halo=False, zorder=2.5,
            )
            ax.text(x + 5.7, y, lbl, fontsize=fs, color=palette.ink2, va="center")
        elif kind == "metro_rev":
            draw_metro_track(
                ax, [(x, y), (x + 5.0, y)], palette.edge_grey,
                width=2.4, hollow=False, surface=palette.card_bg,
                chevrons=True, chevron_spacing=9.0, chevron_bidirectional=True,
                unit_in=unit_in, corner_radius=0.0, halo=False, zorder=2.5,
            )
            ax.text(x + 5.7, y, lbl, fontsize=fs, color=palette.ink2, va="center")
        elif kind == "metro_coupled":
            draw_metro_coupled_halo(
                ax, [[(x, y), (x + 5.0, y)]],
                color=palette.coupled_track, edge_color=palette.coupled_edge, width=7.0, zorder=2.2,
            )
            draw_metro_track(
                ax, [(x, y), (x + 5.0, y)], palette.up,
                width=2.8, hollow=False, surface=palette.card_bg,
                chevrons=True, chevron_spacing=9.0, unit_in=unit_in,
                corner_radius=0.0, halo=False, zorder=2.5,
            )
            ax.text(x + 5.7, y, lbl, fontsize=fs, color=palette.coupled_ink, va="center")
        elif kind == "arrow_flux":
            draw_quad_arrow(ax, (x, y), (x + 2.4, y), (x + 4.8, y),
                            palette.up, 1.9, "-", False, True, 0.0, 0.0, zorder=2.5, surface=palette.card_bg)
            ax.text(x + 5.5, y, lbl, fontsize=fs, color=palette.ink2, va="center")
        elif kind == "arrow_rev":
            draw_quad_arrow(ax, (x, y), (x + 2.4, y), (x + 4.8, y),
                            palette.edge_grey, 1.4, "-", True, True, 0.0, 0.0, zorder=2.5, surface=palette.card_bg)
            ax.text(x + 5.5, y, lbl, fontsize=fs, color=palette.ink2, va="center")
        elif kind == "arrow_coupled":
            draw_coupled_track(ax, [{"p0": (x, y), "c": (x + 2.4, y), "p1": (x + 4.8, y), "trim0": 0, "trim1": 0}],
                               color=palette.coupled_track, edge_color=palette.coupled_edge, width=6.4, zorder=2.2)
            draw_quad_arrow(ax, (x, y), (x + 2.4, y), (x + 4.8, y),
                            palette.up, 1.5, "-", False, True, 0.0, 0.0, zorder=2.5, surface=palette.card_bg)
            ax.text(x + 5.5, y, lbl, fontsize=fs, color=palette.coupled_ink, va="center")
        elif kind == "sign_reversal":
            draw_sign_reversal_badge(ax, x + 0.85, y, radius=0.82, palette=palette, lw=0.70, zorder=2.5)
            ax.text(x + 2.15, y, lbl, fontsize=fs, color=palette.sign_change_ink, fontweight="medium", va="center")
        elif kind == "metab":
            ax.add_patch(Circle((x + 0.9, y), 0.90, facecolor=cmap(DEFAULT_NORM(1.5)), edgecolor=palette.ink, lw=0.75, zorder=2.3))
            ax.add_patch(Circle((x + 0.9, y), 1.08, facecolor="none", edgecolor=palette.ink, lw=1.1, zorder=2.4))
            ax.text(x + 2.5, y, lbl, fontsize=fs, color=palette.ink2, va="center")

    # Row-group B: 4 Layer Icons + Real Visual Cell Status Swatches
    y_b0 = y0 + 1.5 + plan["a_rows"] * row_a + 1.45
    for (kind, lbl, _), (r_idx, x_off) in zip(plan["b_items"], plan["b_pos"]):
        x = vx0 + 2.0 + x_off
        y = y_b0 + r_idx * row_b
        if kind.startswith("icon:"):
            ikind = kind.split(":", 1)[1]
            draw_layer_icon(ax, ikind, x + 0.58, y, size=1.15, color=palette.ink, lw=0.78, zorder=2.5)
            ax.text(x + 1.45, y, lbl, fontsize=fs, fontweight="bold", color=palette.ink2, va="center")
        elif kind == "cell_sig":
            draw_status_cell(ax, x, y - sw / 2, sw, "changed", 1, palette=palette, font_scale=font_scale, zorder=2.5, draw_box=True)
            draw_status_cell(ax, x + sw + 0.22, y - sw / 2, sw, "changed", -1, palette=palette, font_scale=font_scale, zorder=2.5, draw_box=True)
            ax.text(x + 2 * sw + 0.62, y, lbl, fontsize=fs, color=palette.ink2, va="center")
        elif kind == "cell_sug":
            draw_status_cell(ax, x, y - sw / 2, sw, "suggestive", 1, palette=palette, font_scale=font_scale, zorder=2.5, draw_box=True)
            draw_status_cell(ax, x + sw + 0.22, y - sw / 2, sw, "suggestive", -1, palette=palette, font_scale=font_scale, zorder=2.5, draw_box=True)
            ax.text(x + 2 * sw + 0.62, y, lbl, fontsize=fs, color=palette.ink2, va="center")
        else:
            st_map = {
                "cell_mixed": "mixed",
                "cell_unchanged": "unchanged",
                "cell_na": "not_measured",
                "cell_unident": "not_identifiable",
                "cell_nodet": "not_detected",
                "cell_qc": "qc_blocked",
            }
            draw_status_cell(ax, x, y - sw / 2, sw, st_map[kind], 0, palette=palette, font_scale=font_scale, zorder=2.5, draw_box=True)
            ax.text(x + sw + 0.45, y, lbl, fontsize=fs, color=palette.ink2, va="center")

    # Row-group C: Self-explanatory Tier Badges + Linked Flux Badge
    y_c0 = y0 + 1.5 + plan["a_rows"] * row_a + 0.45 + plan["b_rows"] * row_b + 1.65
    c_scale = 1.02 * font_scale
    for (kind, lbl, _), (r_idx, x_off) in zip(plan["c_items"], plan["c_pos"]):
        x = vx0 + 2.0 + x_off
        y = y_c0 + r_idx * row_c
        if kind in ("L3", "L2", "L1"):
            bw = draw_tier_badge(ax, x, y, kind, zorder=2.5, scale=c_scale, preset=preset, unit_in=unit_in)
            ax.text(x + bw + 0.55, y, lbl, fontsize=fs, color=palette.ink2, va="center")
        elif kind == "LINK":
            bw = draw_coupled_badge(ax, x, y, 1, zorder=2.5, scale=c_scale, palette=palette, total_groups=1, unit_in=unit_in)
            ax.text(x + bw + 0.55, y, lbl, fontsize=fs, color=palette.ink2, va="center")

    # Disclosure note
    if plan["note_lines"]:
        y_n = y_c0 + (plan["c_rows"] - 1) * row_c + 2.25
        for line in plan["note_lines"]:
            ax.text(vx0 + 2.0, y_n, line, fontsize=fs * 0.88, color=palette.muted, va="center")
            y_n += 1.80

    # Short calibration / provenance footer
    y = y0 + legend_h + 1.0
    for wrapped in plan["cap_lines"]:
        ax.text(vx0 + 1.0, y, wrapped, fontsize=fs * 0.86, va="top", color=palette.ink2)
        y += 1.65


def _draw_header(
    ax,
    title: str,
    view: pd.DataFrame,
    x_left: float,
    x_right: float,
    *,
    y_base: float = 0.0,
    palette: SemanticPalette,
    preset: str = DEFAULT_PRESET,
    font_scale: float = 1.0,
    unit_in: float = UNIT,
    stacked: bool = False,
    title_fs: float | None = None,
    pill_fs: float | None = None,
) -> None:
    t_fs = title_fs if title_fs is not None else 8.4 * font_scale
    p_fs = pill_fs if pill_fs is not None else 4.15 * font_scale
    pills = _header_pills(view, palette, preset)
    tok = tier_tokens(preset)

    if stacked:
        ax.text(x_left + 1.0, y_base - 6.7, title, fontsize=t_fs, fontweight="bold",
                color=palette.ink, va="center", ha="left")
        py = y_base - 3.1
        px = x_left + 1.0
        dot_scale = p_fs / 4.15
        for text, bg, fg, icon in pills:
            icon_pad = 2.05 * dot_scale if icon else 0.0
            tw = _text_width(text, p_fs, bold=True, unit=unit_in)
            pw = tw + 1.8 + icon_pad
            ax.add_patch(FancyBboxPatch(
                (px, py - 1.12), pw, 2.24, boxstyle="round,pad=0,rounding_size=0.85",
                facecolor=bg, edgecolor=palette.card_edge if icon != "L2" else "none",
                lw=0.6 if icon != "L2" else 0.0,
                linestyle=(0, (2.0, 1.2)) if icon == "L1" else "-",
                zorder=2,
            ))
            if icon in ("L2", "L1"):
                n_on = 2 if icon == "L2" else 1
                dot_col, off_col = tok[icon]["dot_on"], tok[icon]["dot_off"]
                for di in range(3):
                    cx = px + 0.78 + di * 0.48 * dot_scale
                    if di < n_on:
                        ax.add_patch(Circle((cx, py), 0.17 * dot_scale, facecolor=dot_col, edgecolor="none", zorder=2.2))
                    else:
                        ax.add_patch(Circle((cx, py), 0.15 * dot_scale, facecolor="none", edgecolor=off_col, lw=0.45, zorder=2.2))
            ax.text(px + 0.78 + icon_pad, py, text, fontsize=p_fs, fontweight="bold", color=fg,
                    ha="left", va="center", zorder=2.2)
            px += pw + 1.0
    else:
        ax.text(x_left + 1.0, y_base - 3.5, title, fontsize=t_fs, fontweight="bold",
                color=palette.ink, va="center", ha="left")
        px = x_right - 1.0
        dot_scale = p_fs / 4.15
        for text, bg, fg, icon in reversed(pills):
            icon_pad = 2.05 * dot_scale if icon else 0.0
            tw = _text_width(text, p_fs, bold=True, unit=unit_in)
            pw = tw + 1.8 + icon_pad
            ax.add_patch(FancyBboxPatch(
                (px - pw, y_base - 4.70), pw, 2.30, boxstyle="round,pad=0,rounding_size=0.9",
                facecolor=bg, edgecolor=palette.card_edge if icon != "L2" else "none",
                lw=0.6 if icon != "L2" else 0.0,
                linestyle=(0, (2.0, 1.2)) if icon == "L1" else "-",
                zorder=2,
            ))
            if icon in ("L2", "L1"):
                n_on = 2 if icon == "L2" else 1
                dot_col, off_col = tok[icon]["dot_on"], tok[icon]["dot_off"]
                for di in range(3):
                    cx = px - pw + 0.78 + di * 0.48 * dot_scale
                    if di < n_on:
                        ax.add_patch(Circle((cx, y_base - 3.55), 0.17 * dot_scale, facecolor=dot_col, edgecolor="none", zorder=2.2))
                    else:
                        ax.add_patch(Circle((cx, y_base - 3.55), 0.15 * dot_scale, facecolor="none", edgecolor=off_col, lw=0.45, zorder=2.2))
            ax.text(px - pw + 0.78 + icon_pad, y_base - 3.55, text, fontsize=p_fs, fontweight="bold", color=fg,
                    ha="left", va="center", zorder=2.2)
            px -= pw + 1.0
    ax.plot([x_left + 1.0, x_right - 1.0], [y_base - 1.05, y_base - 1.05], color=palette.card_edge, lw=0.85, zorder=1)


# ---------------------------------------------------------------------------
# Main Entry Point
# ---------------------------------------------------------------------------

def plot_pathway_map(
    view: pd.DataFrame,
    layout: Layout,
    regions: pd.DataFrame,
    title: str,
    caption: list[str],
    stem: Path,
    formats: Iterable[str] = ("pdf", "png"),
    extra_fusion: pd.DataFrame | None = None,
    signed_modules: pd.DataFrame | None = None,
    top_rank: int = 10,
    *,
    mode: str = "atlas",
    style: str = "metro_atlas",
    color_mode: str = "semantic",
    strip_mode: str = "omics",
    focus_regions: Sequence[str] | None = None,
    include_reactions: Sequence[str] | None = None,
    exclude_reactions: Sequence[str] | None = None,
    panel_note: str = "",
) -> list[Path]:
    """Render Figure E (Per-Tissue Multi-Omics Pathway Atlas or Focused Sub-Pathway Panel).

    Parameters
    ----------
    mode : ``"atlas"`` | ``"publication"``
        ``"atlas"`` renders the full 13-region network + right dashboard sidebar;
        ``"publication"`` crops to active regions (or ``focus_regions``) with responsive single-/double-column
        typography, smart organelle clipping, and hides non-key micro-strips and exploratory rank tabs.
    style : ``"metro_atlas"`` | ``"nature_editorial"`` | ``"slide_dark"``
        Theme preset name from :mod:`omicsweave.visualization.theme` (default ``"metro_atlas"``).
    color_mode : ``"semantic"`` | ``"pathway"``
        Color encoding for ``style="metro_atlas"``: ``"semantic"`` (default) keeps the unified main-text
        red/blue directional scale; ``"pathway"`` uses muted Morandi pathway-family colors.
    focus_regions : sequence of str, optional
        Explicit region IDs to include when ``mode="publication"`` (empty = auto-select active regions).
    panel_note : str, optional
        Custom concise disclosure note shown inside the publication/sub-map legend bar.
    """
    if mode not in VALID_MODES:
        raise ValueError(f"pathway_map mode 应为 {' / '.join(VALID_MODES)} 之一，实际为 {mode!r}")
    if color_mode not in VALID_COLOR_MODES:
        raise ValueError(f"pathway_map color_mode 应为 {' / '.join(VALID_COLOR_MODES)} 之一，实际为 {color_mode!r}")
    if strip_mode not in VALID_STRIP_MODES:
        raise ValueError(f"pathway_map strip_mode 应为 {' / '.join(VALID_STRIP_MODES)} 之一，实际为 {strip_mode!r}")
    palette = get_palette(style)
    cmap = palette.diverging_cmap()
    metro = style == "metro_atlas"

    all_reactions = view[view["element_type"] == "reaction"]
    all_aggs = view[view["element_type"] == "aggregate"]
    nodes_df = view[view["element_type"] == "metabolite"].set_index("element_id")

    if mode == "publication":
        active_regions, folded_regions = _select_active_regions(view, layout, focus_regions)
        rxn_mask = all_reactions["region"].isin(active_regions)
        if include_reactions:
            rxn_mask = rxn_mask | all_reactions["element_id"].isin(include_reactions)
        if exclude_reactions:
            rxn_mask = rxn_mask & ~all_reactions["element_id"].isin(exclude_reactions)
        reactions = all_reactions[rxn_mask]
        aggs = all_aggs[all_aggs["region"].isin(active_regions)]
        active_aggs = set(aggs["element_id"])
        geometries = {fid: _edge_geometry(layout, fid, metro=metro) for fid in reactions["element_id"]}
        active_species: set[str] = set()
        active_node_positions: set[tuple[str, int]] = set()
        for fid in reactions["element_id"]:
            geom = layout.geometry(pathway_reactions()[fid])
            for sp, pt in geom["substrates"] + geom["products"]:
                active_species.add(sp)
                node = layout.nodes[sp]
                idx = min(range(len(node.at)), key=lambda i: math.dist(node.at[i], pt))
                active_node_positions.add((sp, idx))
        comp_boxes = _compute_clipped_compartments(
            layout, active_regions, active_node_positions, active_aggs, mode="publication",
        )
        vx0, vy0, vx1, vy1 = _compute_publication_viewport(
            layout, active_regions, active_node_positions, active_aggs, comp_boxes,
        )
        vw, vh = vx1 - vx0, vy1 - vy0
        # 按子图跨度自适应单栏 / 1.5 栏 / 双栏物理幅面，保持文字、节点与四轴微条比例和谐
        if vw >= 95.0:
            fig_w = max(DOUBLE_COL_IN * 1.48, vw * 0.092)
        elif vw >= 58.0:
            fig_w = max(DOUBLE_COL_IN * 1.15, vw * 0.116)
        else:
            fig_w = max(ONE_HALF_COL_IN * 1.24, vw * 0.145)
        unit_in = fig_w / (vw + 1.0)
        font_scale = max(1.04, min(1.26, 1.16 * (unit_in / 0.105)))

        header_view = pd.concat([reactions, aggs], ignore_index=True)
        header_h, hdr_stacked, hdr_title_fs, hdr_pill_fs = _header_geometry(
            title, header_view, vw,
            palette=palette, preset=style, font_scale=font_scale, unit_in=unit_in,
        )
        legend_plan = _pub_legend_plan(
            vw, caption, folded_regions, panel_note,
            font_scale=font_scale, unit_in=unit_in, metro=metro, color_mode=color_mode, strip_mode=strip_mode,
        )
        pub_footer = legend_plan["footer_h"]
        fig_h = fig_w * ((vh + header_h + pub_footer) / (vw + 1.0))
        x_lim = (vx0 - 0.5, vx1 + 0.5)
        y_lim = (vy1 + pub_footer + 0.5, vy0 - header_h)
        y_base = vy0
        bounds = (vx0, vy0, vx1, vy1)
    else:
        active_regions = set(layout.regions)
        folded_regions = []
        reactions = all_reactions
        aggs = all_aggs
        active_aggs = set(aggs["element_id"])
        active_species = set(layout.nodes)
        active_node_positions = None
        comp_boxes = _compute_clipped_compartments(layout, active_regions, None, active_aggs, mode="atlas")
        geometries = {fid: _edge_geometry(layout, fid, metro=metro) for fid in reactions["element_id"]}
        font_scale = 1.0
        unit_in = UNIT
        total_w = layout.width + SIDEBAR
        total_h = layout.height + FOOTER
        header_view = view
        header_h, hdr_stacked, hdr_title_fs, hdr_pill_fs = HEADER, False, 8.4, 4.15
        legend_plan = None
        fig_w = total_w * UNIT
        fig_h = (total_h + header_h) * UNIT
        x_lim = (-0.5, total_w + 0.5)
        y_lim = (total_h + 0.5, -header_h)
        y_base = 0.0
        bounds = (0.0, 0.0, layout.width, layout.height)

    with themed(preset=style, dpi=150, save_dpi=300):
        fig = plt.figure(figsize=(fig_w, fig_h), facecolor=palette.surface)
        ax = fig.add_axes([0, 0, 1, 1], facecolor=palette.surface)
        ax.set_xlim(*x_lim)
        ax.set_ylim(*y_lim)
        ax.set_aspect("equal")
        ax.axis("off")

        _draw_header(
            ax, title, header_view,
            x_lim[0] + 0.5, x_lim[1] - 0.5,
            y_base=y_base, palette=palette, preset=style, font_scale=font_scale, unit_in=unit_in,
            stacked=hdr_stacked, title_fs=hdr_title_fs, pill_fs=hdr_pill_fs,
        )
        _draw_regions_and_compartments(
            ax, layout, active_regions,
            palette=palette, font_scale=font_scale, metro=metro, color_mode=color_mode, comp_boxes=comp_boxes,
        )

        group_indices = _flux_group_indices(reactions)
        # 仅当同一 flux_group 至少有 2 条方向性边时才画轨道包络（单边无需耦合轨道）
        grp_counts: dict[int, int] = {}
        for _, r in reactions.iterrows():
            if r["flux_group"] and r["F_state"] in ("changed", "suggestive"):
                gi = group_indices.get(str(r["flux_group"]).strip())
                if gi is not None:
                    grp_counts[gi] = grp_counts.get(gi, 0) + 1
        multi_groups = {gi for gi, cnt in grp_counts.items() if cnt >= 2}

        node_labels = _node_label_placements(layout, geometries, bounds, font_scale=font_scale, unit_in=unit_in)
        positions = _place_strips(
            reactions, layout, geometries, node_labels, group_indices,
            active_species=active_species, active_aggs=active_aggs, active_regions=active_regions,
            bounds=bounds, top_rank=top_rank, mode=mode, font_scale=font_scale, unit_in=unit_in,
            metro=metro, comp_boxes=comp_boxes,
        )

        draw_order = (
            reactions["draw_style"].map({"not_assessed": 0, "omics_only": 1, "flux": 2}).fillna(0)
            + reactions["value"].notna().astype(int) * 2
            + reactions["key"].eq("key").astype(int) * 2
        )
        for _, row in reactions.assign(_z=draw_order).sort_values("_z", kind="stable").iterrows():
            gi = group_indices.get(str(row["flux_group"]).strip()) if (row["flux_group"] and row["F_state"] in ("changed", "suggestive")) else None
            _draw_edge(
                ax, row, geometries[row["element_id"]],
                palette=palette, cmap=cmap, layout=layout,
                is_coupled_member=(gi in multi_groups),
                font_scale=font_scale, metro=metro, color_mode=color_mode, unit_in=unit_in,
            )

        for _, row in aggs.iterrows():
            _draw_aggregate(
                ax, layout, row, palette=palette, cmap=cmap, preset=style,
                strip_mode=strip_mode, font_scale=font_scale, unit_in=unit_in,
            )

        for species in layout.nodes:
            if species not in active_species:
                continue
            _draw_node(
                ax, layout, species,
                nodes_df.loc[species] if species in nodes_df.index else None,
                node_labels, palette=palette, cmap=cmap, font_scale=font_scale,
                unit_in=unit_in, active_positions=active_node_positions, metro=metro,
            )

        for _, row in reactions.iterrows():
            _draw_leader(
                ax, geometries[row["element_id"]]["hub"], positions[row["element_id"]], row, group_indices,
                palette=palette, top_rank=top_rank, mode=mode, font_scale=font_scale, unit_in=unit_in,
            )
            _draw_strip(
                ax, positions[row["element_id"]], row, group_indices,
                palette=palette, preset=style, top_rank=top_rank, mode=mode,
                strip_mode=strip_mode, font_scale=font_scale, unit_in=unit_in,
            )

        if mode == "atlas":
            _sidebar(
                ax, regions, layout, layout.width + 2.4, 1.0,
                palette=palette, extra_fusion=extra_fusion, signed_modules=signed_modules,
            )
            _legend_atlas(
                ax, layout.height + 1.8, layout.width + SIDEBAR, caption,
                palette=palette, cmap=cmap, preset=style, metro=metro,
                color_mode=color_mode, strip_mode=strip_mode, unit_in=unit_in,
            )
        else:
            _legend_publication(
                ax, vx0, vy1 + 1.6, vw, legend_plan,
                palette=palette, cmap=cmap, preset=style, font_scale=font_scale, unit_in=unit_in,
            )

        return save_stem(fig, stem, formats=tuple(formats), dpi=300, facecolor=palette.surface, close=True)
