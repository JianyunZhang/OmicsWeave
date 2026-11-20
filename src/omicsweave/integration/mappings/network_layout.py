#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
通路总图的布局（G18 M1，设计文档 ``docs/arch_design/omicsweave_pathway_map.md`` 第 3 节）：读取与校验
``data/network_layout.yaml``。

- 节点是网络物种（同一代谢物在不同区室各占一个节点）。一个物种可以有多个位置（克隆）：连接到它的反应取离反应中心
  最近的一个，用来避免 GLU_c、ASP_c 这类枢纽代谢物拉出很长的线；
- 边是网络反应：底物 → 反应中心 → 产物；中心默认取各端点的重心，可用 ``hub`` 指定，``via`` 为中心前的折点；
- 聚合节点是不在网络中、只有组学证据的评估 ID，画成方框；
- ``exclude`` 中的 ID 不上图，必须写原因。

校验（:func:`validate_layout`）：网络表的每个反应都在 ``edges`` 或 ``exclude`` 中；每条边用到的物种都有坐标
（``hide_species`` 除外）；区域与区室都有声明；``network_sha256`` 与包内网络表一致。
"""
from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable

import yaml

from .network_annotation import DATA_DIR, TABLE_FILE, NetworkReaction, network_reactions, pathway_reactions

LAYOUT_FILE = DATA_DIR / "network_layout.yaml"
NODE_KINDS = ("metabolite", "boundary", "store")
EXCLUDE_REASONS = ("pseudo", "serum_only", "alias", "duplicate_of", "covered_by_steps", "not_drawn")


class LayoutError(ValueError):
    """布局文件格式错误或与网络表不一致。"""


@dataclass(frozen=True)
class Node:
    species: str
    at: tuple[tuple[float, float], ...]
    kind: str = "metabolite"
    label: str = ""


@dataclass(frozen=True)
class Edge:
    flux_id: str
    region: str
    hub: tuple[float, float] | None = None
    via: tuple[tuple[float, float], ...] = ()
    offset: float = 0.0
    strip: tuple[float, float] | None = None    # 四轴小条左端中点相对反应中心的位置；不写时由绘图按避让自动放置
    secondary: tuple[str, ...] = ()             # 次要底物 / 产物（转氨的谷氨酸与 α-酮戊二酸、合成中的谷氨酰胺等）：不画线，在中心旁写名称
    curve: float = 0.0                          # 弧线法向弯曲度（画布单位；0 表示直线或仅由 offset 拱弯）


@dataclass(frozen=True)
class Aggregate:
    flux_id: str
    region: str
    xy: tuple[float, float]
    label: str = ""


@dataclass(frozen=True)
class Region:
    id: str
    label: str
    label_xy: tuple[float, float]
    kegg: tuple[str, ...] = ()
    terms: tuple[str, ...] = ()
    box: tuple[float, float, float, float] | None = None   # 可选分区柔和底板 [x0, y0, x1, y1]
    accent: str = ""                                       # 可选分区主题色（十六进制）


@dataclass(frozen=True)
class Layout:
    path: Path
    model_id: str
    network_sha256: str
    width: float
    height: float
    compartments: dict[str, dict[str, Any]]
    regions: dict[str, Region]
    hide_species: frozenset[str]
    nodes: dict[str, Node]
    edges: dict[str, Edge]
    aggregates: dict[str, Aggregate]
    exclude: dict[str, dict[str, Any]]
    sha256: str = ""
    notes: list[str] = field(default_factory=list)

    def endpoint(self, species: str, near: tuple[float, float]) -> tuple[float, float]:
        """物种离 ``near`` 最近的位置。"""
        return min(self.nodes[species].at, key=lambda p: math.dist(p, near))

    def geometry(self, reaction: NetworkReaction) -> dict[str, Any]:
        """一条边的几何：{"hub", "substrates": [(物种, 点)], "products": [...], "secondary": ([底物], [产物])}；
        隐藏的物种不画，次要物种只写名称。"""
        edge = self.edges[reaction.flux_id]
        skip = self.hide_species | set(edge.secondary)
        subs = [s for s in reaction.substrates if s not in skip]
        pros = [s for s in reaction.products if s not in skip]
        sides = list(dict.fromkeys(subs)) + list(dict.fromkeys(pros))
        if edge.hub is not None:
            hub = edge.hub
        else:
            hub = _centroid([self.nodes[s].at[0] for s in sides])
            for _ in range(3):        # 克隆：先取第一个位置的重心，再改取离重心最近的位置
                hub = _centroid([self.endpoint(s, hub) for s in sides])
        return {"hub": hub,
                "substrates": [(s, self.endpoint(s, hub)) for s in dict.fromkeys(subs)],
                "products": [(s, self.endpoint(s, hub)) for s in dict.fromkeys(pros)],
                "secondary": ([s for s in reaction.substrates if s in edge.secondary],
                              [s for s in reaction.products if s in edge.secondary])}


def _centroid(points: list[tuple[float, float]]) -> tuple[float, float]:
    if not points:
        return (0.0, 0.0)
    return (sum(p[0] for p in points) / len(points), sum(p[1] for p in points) / len(points))


def _point(value: Any, where: str) -> tuple[float, float]:
    if not (isinstance(value, (list, tuple)) and len(value) == 2 and all(isinstance(v, (int, float)) for v in value)):
        raise LayoutError(f"{where}：坐标应为 [x, y]，实际为 {value!r}")
    return float(value[0]), float(value[1])


def network_table_sha256() -> str:
    path = DATA_DIR / TABLE_FILE
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest() if path.exists() else ""


def load_layout(path: str | Path | None = None) -> Layout:
    path = Path(path) if path else LAYOUT_FILE
    if not path.is_file():
        raise LayoutError(f"找不到布局文件 {path}")
    raw = path.read_bytes()
    doc = yaml.safe_load(raw.decode("utf-8")) or {}
    canvas = doc.get("canvas") or {}
    regions = {}
    for rid, spec in (doc.get("regions") or {}).items():
        spec = spec or {}
        raw_box = spec.get("box")
        box = None
        if raw_box is not None:
            if not (isinstance(raw_box, (list, tuple)) and len(raw_box) == 4
                    and all(isinstance(v, (int, float)) for v in raw_box)):
                raise LayoutError(f"regions.{rid}.box：应为 [x0, y0, x1, y1]，实际为 {raw_box!r}")
            box = tuple(float(v) for v in raw_box)
        regions[str(rid)] = Region(str(rid), str(spec.get("label", rid)), _point(spec.get("label_xy"), f"regions.{rid}"),
                                   tuple(spec.get("kegg") or ()), tuple(spec.get("terms") or ()),
                                   box, str(spec.get("accent", "")))
    nodes = {}
    for species, spec in (doc.get("nodes") or {}).items():
        if isinstance(spec, list):
            spec = {"at": [spec]}
        spec = dict(spec or {})
        if "xy" in spec:
            spec["at"] = [spec.pop("xy")]
        at = tuple(_point(p, f"nodes.{species}") for p in spec.get("at") or [])
        if not at:
            raise LayoutError(f"nodes.{species}：没有坐标")
        kind = str(spec.get("kind", "metabolite"))
        if kind not in NODE_KINDS:
            raise LayoutError(f"nodes.{species}：kind 应为 {'、'.join(NODE_KINDS)}，实际为 {kind}")
        nodes[str(species)] = Node(str(species), at, kind, str(spec.get("label", "")))
    edges = {}
    for flux_id, spec in (doc.get("edges") or {}).items():
        spec = spec or {}
        hub = _point(spec["hub"], f"edges.{flux_id}.hub") if spec.get("hub") is not None else None
        via = tuple(_point(p, f"edges.{flux_id}.via") for p in spec.get("via") or [])
        strip = _point(spec["strip"], f"edges.{flux_id}.strip") if spec.get("strip") is not None else None
        edges[str(flux_id)] = Edge(str(flux_id), str(spec.get("region", "")), hub, via, float(spec.get("offset", 0.0)),
                                   strip, tuple(str(s) for s in spec.get("secondary") or ()),
                                   float(spec.get("curve", 0.0)))
    aggregates = {str(k): Aggregate(str(k), str((v or {}).get("region", "")), _point((v or {}).get("xy"), f"aggregates.{k}"),
                                    str((v or {}).get("label", "")))
                  for k, v in (doc.get("aggregates") or {}).items()}
    exclude = {str(k): dict(v or {}) for k, v in (doc.get("exclude") or {}).items()}
    return Layout(path, str(doc.get("model_id", "")), str(doc.get("network_sha256", "")),
                  float(canvas.get("width", 100)), float(canvas.get("height", 100)),
                  {str(k): dict(v or {}) for k, v in (doc.get("compartments") or {}).items()}, regions,
                  frozenset(str(s) for s in doc.get("hide_species") or ()), nodes, edges, aggregates, exclude,
                  hashlib.sha256(raw.replace(b"\r\n", b"\n")).hexdigest())


@lru_cache(maxsize=4)
def default_layout() -> Layout:
    return load_layout(LAYOUT_FILE)


def validate_layout(layout: Layout, reactions: dict[str, NetworkReaction] | None = None,
                    check_sha: bool = True) -> list[str]:
    """返回问题列表（空表示通过）。MFA 网络表中的每个反应都须在 edges 或 exclude 中；
    edges 同时允许使用 ``curated_reactions.tsv`` 中的策展扩展反应（便于后续 MFA 扩展无缝承接）。"""
    mfa_reactions = network_reactions() if reactions is None else reactions
    all_reactions = pathway_reactions() if reactions is None else reactions
    problems: list[str] = []
    if check_sha and layout.network_sha256 != network_table_sha256():
        problems.append("network_sha256 与包内 network_reactions.tsv 不一致：网络表改过，需要核对坐标后更新校验值")
    for flux_id in mfa_reactions:
        if flux_id not in layout.edges and flux_id not in layout.exclude:
            problems.append(f"网络反应 {flux_id} 既不在 edges 中，也不在 exclude 中")
    for flux_id in layout.edges:
        if flux_id not in all_reactions:
            problems.append(f"edges.{flux_id} 不是网络反应或策展扩展反应")
            continue
        if flux_id in layout.exclude:
            problems.append(f"{flux_id} 同时在 edges 与 exclude 中")
        edge = layout.edges[flux_id]
        if edge.region not in layout.regions:
            problems.append(f"edges.{flux_id}：区域 {edge.region!r} 没有声明")
        reaction = all_reactions[flux_id]
        for species in edge.secondary:
            if species not in reaction.substrates + reaction.products:
                problems.append(f"edges.{flux_id}：次要物种 {species} 不是该反应的底物或产物")
        for species in reaction.substrates + reaction.products:
            if species not in layout.hide_species and species not in edge.secondary and species not in layout.nodes:
                problems.append(f"edges.{flux_id}：物种 {species} 没有坐标")
    for flux_id, agg in layout.aggregates.items():
        if flux_id in mfa_reactions or flux_id in layout.edges:
            problems.append(f"aggregates.{flux_id} 是网络反应，应画成边")
        if agg.region not in layout.regions:
            problems.append(f"aggregates.{flux_id}：区域 {agg.region!r} 没有声明")
    for flux_id, spec in layout.exclude.items():
        reason = spec.get("reason")
        if reason not in EXCLUDE_REASONS:
            problems.append(f"exclude.{flux_id}：原因应为 {'、'.join(EXCLUDE_REASONS)} 之一，实际为 {reason!r}")
    used: dict[str, set] = {}
    for flux_id, edge in layout.edges.items():
        reaction = all_reactions.get(flux_id)
        skip = layout.hide_species | set(edge.secondary)
        if reaction is not None and all(s in layout.nodes for s in reaction.substrates + reaction.products
                                        if s not in skip):
            geometry = layout.geometry(reaction)
            for species, point in geometry["substrates"] + geometry["products"]:
                used.setdefault(species, set()).add(point)
    for species, node in layout.nodes.items():
        if species not in used:
            problems.append(f"nodes.{species} 没有被任何画出的反应用到")
        else:
            for point in node.at:
                if point not in used[species]:
                    problems.append(f"nodes.{species} 的位置 {list(point)} 没有反应连到")
    return problems


def unplaced_ids(layout: Layout, assessed: Iterable[str]) -> list[str]:
    """评估集中既不是边、也不是聚合节点、也没有排除原因的 ID。"""
    return sorted(i for i in set(map(str, assessed))
                  if i not in layout.edges and i not in layout.aggregates and i not in layout.exclude)
