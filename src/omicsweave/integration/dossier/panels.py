#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
通路面板（G14）：课题目录中的 ``panels.yaml`` 声明有序的反应、分支节点、预期示踪剂、预先定号的模块与对照集；
本模块解析声明并生成面板的行（反应与代谢物节点交替，``panel_elements``）。

代谢物节点由网络化学计量补全：每个反应上方是它的底物，下方是产物；相邻两个反应中，前一个的产物恰是后一个的
底物时，节点只出现一次。辅因子（CO2、ATP、NAD(H)、CoA 等）不作为节点。

格式（顶层键 ``panels``，每个面板）::

    title, tissues, contrasts                      # 标题；只为这些组织 / 比较出图（省略时全部）
    segments: [{name, reactions: [...]}]           # 自上而下的阅读顺序
    branch_nodes: [{node: PYR_m, exits: {出口名: [反应...]}}]
    expected_tracers: {反应: [示踪剂...]}
    signed_modules: {模块: {plus: [...], minus: [...]}}
    controls: {positive: [...], negative: [...]}   # 校准（S4）的对照：反应 ID，或 {reaction, direction, tissue, contrast, axis}
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
import yaml

from ..mappings.network_annotation import network_reactions, split_metabolite

PANEL_KEYS = {"title", "tissues", "contrasts", "segments", "branch_nodes", "expected_tracers", "signed_modules",
              "controls"}
ELEMENT_COLUMNS = ["panel_id", "order", "element_type", "element_id", "label", "segment", "branch_node",
                   "exit_label"]
COFACTORS = {"CO2", "HCO3", "ATP", "ADP", "AMP", "GTP", "GDP", "UTP", "UDP", "NAD", "NADH", "NADP", "NADPH", "FAD",
             "FADH2", "COA", "PI", "PPI", "H2O", "H", "O2", "NH3", "NH4", "Q", "QH2"}


class PanelError(ValueError):
    pass


@dataclass
class Panel:
    id: str
    title: str
    segments: list[tuple[str, list[str]]]
    tissues: list[str] = field(default_factory=list)
    contrasts: list[str] = field(default_factory=list)
    branch_nodes: list[dict[str, Any]] = field(default_factory=list)
    expected_tracers: dict[str, list[str]] = field(default_factory=dict)
    signed_modules: dict[str, dict[str, list[str]]] = field(default_factory=dict)
    controls: list[dict[str, Any]] = field(default_factory=list)   # 展开后的对照（每个组织、比较一条）

    def expand_controls(self, tissues: Iterable[str], contrasts: Iterable[str]) -> list[dict[str, Any]]:
        """对照中省略的组织 / 比较取面板声明的（面板也没写时取数据中的全部）。"""
        out = []
        for control in self.controls:
            for tissue in [control["tissue"]] if control.get("tissue") else (self.tissues or list(tissues)):
                for contrast in [control["contrast"]] if control.get("contrast") else (self.contrasts or list(contrasts)):
                    out.append({**control, "tissue": tissue, "contrast": contrast})
        return out

    @property
    def reactions(self) -> list[str]:
        return [r for _, reactions in self.segments for r in reactions]


def _main_species(species: tuple[str, ...]) -> list[str]:
    return [s for s in species if split_metabolite(s)[0].upper() not in COFACTORS]


@lru_cache(maxsize=1)
def species_labels() -> dict[str, str]:
    """网络物种 → 可读名称：取各反应的组织层名称（与非辅因子物种一一对应时）；没有时用物种名的基名。"""
    names: dict[str, str] = {}
    for reaction in network_reactions().values():
        for species, tissue_names in ((reaction.substrates, reaction.tissue_substrate_names),
                                      (reaction.products, reaction.tissue_product_names)):
            main = _main_species(species)
            if tissue_names and len(main) == len(tissue_names):
                for s, name in zip(main, tissue_names):
                    names.setdefault(split_metabolite(s)[0], name)
    return names


def species_label(species: str) -> str:
    base, compartment = split_metabolite(species)
    name = species_labels().get(base, base)
    return f"{name} ({compartment})" if compartment else name


def load_panels(path: str | Path) -> dict[str, Panel]:
    path = Path(path)
    if not path.is_file():
        raise PanelError(f"找不到面板声明 {path}")
    doc = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    raw = doc.get("panels") if isinstance(doc, dict) else None
    if not isinstance(raw, dict) or not raw:
        raise PanelError(f"{path.name}：顶层应为 panels: {{面板 ID: 声明}}")
    network = network_reactions()
    panels: dict[str, Panel] = {}
    for panel_id, spec in raw.items():
        where = f"{path.name} 的面板 {panel_id}"
        if not isinstance(spec, dict):
            raise PanelError(f"{where}：应为映射")
        unknown = sorted(set(spec) - PANEL_KEYS)
        if unknown:
            raise PanelError(f"{where}：无法识别的键 {unknown}，可选 {sorted(PANEL_KEYS)}")
        segments = []
        for i, segment in enumerate(spec.get("segments") or []):
            if not isinstance(segment, dict) or not segment.get("reactions"):
                raise PanelError(f"{where}：第 {i + 1} 段应为 {{name, reactions: [...]}}")
            segments.append((str(segment.get("name") or f"段 {i + 1}"), [str(r) for r in segment["reactions"]]))
        if not segments:
            raise PanelError(f"{where}：没有 segments")
        panel = Panel(str(panel_id), str(spec.get("title") or panel_id), segments,
                      [str(t) for t in spec.get("tissues") or []], [str(c) for c in spec.get("contrasts") or []],
                      list(spec.get("branch_nodes") or []),
                      {str(k): [str(t) for t in v] for k, v in (spec.get("expected_tracers") or {}).items()},
                      {str(k): {s: [str(m) for m in (v or {}).get(s, [])] for s in ("plus", "minus")}
                       for k, v in (spec.get("signed_modules") or {}).items()},
                      _controls(spec.get("controls") or {}, where))
        missing = [r for r in panel.reactions if r not in network]
        if missing:
            raise PanelError(f"{where}：反应 {missing} 不在 MFA 网络中")
        duplicated = sorted({r for r in panel.reactions if panel.reactions.count(r) > 1})
        if duplicated:
            raise PanelError(f"{where}：反应 {duplicated} 出现了不止一次")
        stray = sorted(set(panel.expected_tracers) - set(panel.reactions))
        if stray:
            raise PanelError(f"{where}：expected_tracers 中的反应 {stray} 不在面板中")
        stray = sorted({c["reaction"] for c in panel.controls} - set(panel.reactions))
        if stray:
            raise PanelError(f"{where}：controls 中的反应 {stray} 不在面板中")
        nodes = {s for r in panel.reactions for s in (*network[r].substrates, *network[r].products)}
        for node in panel.branch_nodes:
            exits = {str(r) for rs in (node.get("exits") or {}).values() for r in rs}
            if str(node.get("node")) not in nodes or not exits <= set(panel.reactions):
                raise PanelError(f"{where}：分支节点 {node.get('node')} 不是面板反应的代谢物，或出口反应不在面板中")
        panels[panel.id] = panel
    return panels


CONTROL_AXES = ("claim", "flux", "capacity", "regulation", "pool")


def _controls(spec: Any, where: str) -> list[dict[str, Any]]:
    """``controls: {positive: [...], negative: [...]}``。每条为反应 ID，或 {reaction, direction（阳性必填：up / down），
    tissue, contrast, axis（claim 或某条轴，默认 claim）}。"""
    if not isinstance(spec, dict) or set(spec) - {"positive", "negative"}:
        raise PanelError(f"{where}：controls 应为 {{positive: [...], negative: [...]}}")
    out = []
    for kind in ("positive", "negative"):
        for item in spec.get(kind) or []:
            entry = {"reaction": item} if isinstance(item, str) else dict(item or {})
            entry = {**entry, "kind": kind, "axis": str(entry.get("axis") or "claim")}
            if not entry.get("reaction"):
                raise PanelError(f"{where}：controls.{kind} 的条目需要 reaction")
            if kind == "positive" and entry.get("direction") not in ("up", "down"):
                raise PanelError(f"{where}：阳性对照 {entry['reaction']} 需要 direction: up 或 down")
            if entry["axis"] not in CONTROL_AXES:
                raise PanelError(f"{where}：对照 {entry['reaction']} 的 axis 应为 {CONTROL_AXES} 之一")
            out.append(entry)
    return out


def panel_elements(panel: Panel) -> pd.DataFrame:
    """面板的行：代谢物节点（``metabolite``，ID 为网络物种）与反应（``reaction``）按阅读顺序交替。"""
    network = network_reactions()
    branch = {str(n.get("node")): n for n in panel.branch_nodes}
    exits = {str(r): str(name) for n in panel.branch_nodes for name, rs in (n.get("exits") or {}).items() for r in rs}
    rows: list[list[Any]] = []

    def add(kind: str, element: str, label: str, segment: str) -> None:
        rows.append([panel.id, len(rows), kind, element, label, segment,
                     kind == "metabolite" and element in branch, exits.get(element, "") if kind == "reaction" else ""])

    previous: set[str] = set()
    for segment, reactions in panel.segments:
        for reaction in reactions:
            entry = network[reaction]
            for species in _main_species(entry.substrates):
                if species not in previous:
                    add("metabolite", species, species_label(species), segment)
            add("reaction", reaction, reaction, segment)
            products = _main_species(entry.products)
            for species in products:
                add("metabolite", species, species_label(species), segment)
            previous = set(products)
    return pd.DataFrame(rows, columns=ELEMENT_COLUMNS)
