#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
反应成员知识表（G17）：成员的角色、直接酶之间的关系（同工酶 / 复合体亚基 / 共用亚基），以及调控因子与别构效应物对
酶活性的作用方向。数据在同目录的 ``reaction_members.yaml``，只约束登记在其中的反应。

- 角色统一：登记反应在基因表（``gene_roles``）与蛋白表（各分类）中的角色由本表生成（:func:`sync_role_tables`，
  即 ``omicsweave refs build reaction_roles``）。生成只改角色分组，不增删两张表的成员，所以 ``flux_anchored``
  匹配到的成员集合不变；两张表都没有的成员只记在本表（如别构效应物），供反应档案展示。
- 成员关系：同工酶加权只在 ``isozyme`` 之间进行；``shared_subunit`` 不作为方向成员（见 ``flux_evidence``）。
- 审计（:func:`audit_reaction_roles`，``omicsweave refs audit``）：两张表对同一成员的角色不一致；登记反应的表与本表
  不同步、或表中有本表没有的成员；调控因子没有作用方向；效应物名称查不到代谢物 ID。
"""
from __future__ import annotations

import argparse
import copy
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping

import yaml

from .flux_protein_mapping import PROTEIN_CATEGORIES

_DIR = Path(__file__).resolve().parent
ANNOTATION_FILE = _DIR / "reaction_members.yaml"
GENE_TABLE = _DIR / "flux_gene_mapping.yaml"
PROTEIN_TABLE = _DIR / "flux_protein_mapping.yaml"


@dataclass(frozen=True)
class RoleSpec:
    role: str               # 反应档案中的角色名
    protein_category: str   # flux_protein_mapping.yaml 的分类；效应物为空（代谢物不在蛋白表）
    gene_role: str          # flux_gene_mapping.yaml 的 gene_roles 键；效应物为空


# 键为 reaction_members.yaml 中的角色键
ROLES: dict[str, RoleSpec] = {
    "direct_enzymes": RoleSpec("direct_enzyme", "direct_enzymes", "catalytic"),
    "carriers": RoleSpec("carrier", "transporters", "transporter"),
    "regulators": RoleSpec("regulator", "upstream_regulators", "regulatory"),
    "effectors": RoleSpec("effector", "", ""),
    "transcription_factors": RoleSpec("transcription_factor", "transcription_factors", "transcription_factor"),
    "upstream_enzymes": RoleSpec("upstream_enzyme", "upstream_enzymes", "upstream"),
    "downstream_enzymes": RoleSpec("downstream_enzyme", "downstream_enzymes", "downstream"),
    "related_isozymes": RoleSpec("related_isozyme", "related_isozymes", "related_isozyme"),
    "cofactors": RoleSpec("cofactor", "cofactors", "cofactor"),
    "pathway_neighbors": RoleSpec("pathway_neighbor", "pathway_neighbors", "pathway_neighbor"),
}
PROTEIN_CATEGORY_ROLE = {s.protein_category: s.role for s in ROLES.values() if s.protein_category}
GENE_ROLE_ROLE = {s.gene_role: s.role for s in ROLES.values() if s.gene_role}
_ROLE_SPEC = {s.role: s for s in ROLES.values()}
GENE_ROLE_ORDER = tuple(s.gene_role for s in ROLES.values() if s.gene_role)

RELATIONS = ("isozyme", "complex_subunit", "shared_subunit", "pathway_step")
EFFECTS = ("activates", "inhibits", "unknown")
ACTIVITY_MODES = ("phosphorylation", "dephosphorylation", "allosteric", "protein_interaction")
MODES = (*ACTIVITY_MODES, "transcription", "signaling", "unknown")
_SIGNED_ROLES = ("regulator", "effector")
_RELATION_ROLES = ("direct_enzyme", "carrier")
_ATTRIBUTES = {"relation", "effect", "mode", "sources", "note"}


@dataclass(frozen=True)
class Member:
    id: str
    role: str
    relation: str = ""
    effect: str = ""
    mode: str = ""
    sources: tuple[str, ...] = ()
    note: str = ""

    @property
    def activity_sign(self) -> int:
        """对酶活性的作用方向：激活 +1、抑制 −1；方向未知、或经由酶量 / 其他通路起作用（转录、信号）时为 0。"""
        if self.role not in _SIGNED_ROLES or self.mode not in ACTIVITY_MODES:
            return 0
        return {"activates": 1, "inhibits": -1}.get(self.effect, 0)


@dataclass(frozen=True)
class ReactionMembers:
    flux_id: str
    members: dict[str, Member] = field(default_factory=dict)

    def by_role(self, role: str) -> list[Member]:
        return [m for m in self.members.values() if m.role == role]

    def relations(self) -> dict[str, str]:
        return {m.id: m.relation for m in self.members.values() if m.role in _RELATION_ROLES}


def _parse_reaction(flux_id: str, entry: Any) -> ReactionMembers:
    if not isinstance(entry, Mapping):
        raise ValueError(f"{flux_id}：应为 {{角色: 成员}}")
    members: dict[str, Member] = {}
    for key, values in entry.items():
        if key not in ROLES:
            raise ValueError(f"{flux_id}：未知角色 {key!r}，可选：{'、'.join(ROLES)}")
        role = ROLES[key].role
        items = values.items() if isinstance(values, Mapping) else ((v, None) for v in (values or []))
        for member_id, attrs in items:
            member_id = str(member_id)
            attrs = dict(attrs or {})
            unknown = set(attrs) - _ATTRIBUTES
            if unknown:
                raise ValueError(f"{flux_id}.{member_id}：未知属性 {sorted(unknown)}")
            if member_id in members:
                raise ValueError(f"{flux_id}.{member_id}：同时属于 {members[member_id].role} 与 {role}")
            relation = str(attrs.get("relation") or ("isozyme" if role in _RELATION_ROLES else ""))
            if relation and (role not in _RELATION_ROLES or relation not in RELATIONS):
                raise ValueError(f"{flux_id}.{member_id}：relation {relation!r} 只用于直接酶与转运体，可选：{RELATIONS}")
            effect, mode = str(attrs.get("effect") or ""), str(attrs.get("mode") or "")
            if role in _SIGNED_ROLES:
                if effect not in EFFECTS or mode not in MODES:
                    raise ValueError(f"{flux_id}.{member_id}：{role} 必须写 effect（{EFFECTS}）与 mode（{MODES}）")
            elif effect or mode:
                raise ValueError(f"{flux_id}.{member_id}：effect / mode 只用于调控因子与效应物")
            sources = tuple(str(s) for s in (attrs.get("sources") or ()))
            members[member_id] = Member(member_id, role, relation, effect, mode, sources, str(attrs.get("note") or ""))
    return ReactionMembers(flux_id, members)


def load_reaction_members(path: str | Path = ANNOTATION_FILE) -> dict[str, ReactionMembers]:
    data = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
    if not isinstance(data, Mapping):
        raise ValueError(f"{path}：根节点应为 {{反应: 角色}}")
    return {str(k): _parse_reaction(str(k), v) for k, v in data.items()}


REACTION_MEMBERS: dict[str, ReactionMembers] = load_reaction_members()


def reaction_members(flux_id: str) -> ReactionMembers | None:
    return REACTION_MEMBERS.get(flux_id)


def member_relations(flux_id: str) -> dict[str, str]:
    """登记反应中直接酶与转运体的关系；未登记的反应返回空（调用方按同工酶处理）。"""
    annotated = REACTION_MEMBERS.get(flux_id)
    return annotated.relations() if annotated else {}


# ---------------------------------------------------------------------------
# 由本表生成两张映射表的角色分组
# ---------------------------------------------------------------------------

def _dedupe(values: Iterable[Any]) -> list[Any]:
    return list(dict.fromkeys(values))


def _protein_entry(entry: dict[str, Any], annotated: ReactionMembers) -> dict[str, Any]:
    regrouped: dict[str, list[str]] = defaultdict(list)
    for category in PROTEIN_CATEGORIES:
        for member in entry.get(category, []) or []:
            known = annotated.members.get(str(member))
            target = _ROLE_SPEC[known.role].protein_category if known and _ROLE_SPEC[known.role].protein_category \
                else category       # 本表没有的成员留在原分类，由审计报告
            if member not in regrouped[target]:
                regrouped[target].append(member)
    out: dict[str, Any] = {}
    added = [c for c in PROTEIN_CATEGORIES if regrouped.get(c) and c not in entry]
    for key, value in entry.items():
        if key == "description":
            out.update({c: regrouped[c] for c in added})
            added = []
        out[key] = list(regrouped.get(key, [])) if key in PROTEIN_CATEGORIES else value
    out.update({c: regrouped[c] for c in added})
    return out


def _gene_entry(entry: dict[str, Any], annotated: ReactionMembers) -> dict[str, Any]:
    roles = entry.get("gene_roles", {}) or {}
    ordered = _dedupe([m for values in roles.values() for m in (values or [])] + list(entry.get("genes", []) or []))
    original_roles = defaultdict(list)
    for gene_role, values in roles.items():
        for member in values or []:
            original_roles[member].append(gene_role)
    regrouped: dict[str, list[str]] = defaultdict(list)
    for member in ordered:
        known = annotated.members.get(str(member))
        targets = [_ROLE_SPEC[known.role].gene_role] if known and _ROLE_SPEC[known.role].gene_role \
            else original_roles.get(member, [])
        for target in targets:
            if member not in regrouped[target]:
                regrouped[target].append(member)
    keys = [k for k in roles if regrouped.get(k)] + [k for k in GENE_ROLE_ORDER if regrouped.get(k) and k not in roles]
    keys += [k for k in regrouped if regrouped[k] and k not in keys]
    out = dict(entry)
    out["gene_roles"] = {k: regrouped[k] for k in keys}
    return out


def derive_role_tables(gene_map: Mapping[str, Any], protein_map: Mapping[str, Any],
                       annotations: Mapping[str, ReactionMembers] | None = None
                       ) -> tuple[dict[str, Any], dict[str, Any], list[str]]:
    """返回按本表改写角色分组后的 (基因表, 蛋白表, 有改动的反应)。未登记的反应原样保留。"""
    annotations = REACTION_MEMBERS if annotations is None else annotations
    genes, proteins = copy.deepcopy(dict(gene_map)), copy.deepcopy(dict(protein_map))
    changed: list[str] = []
    for flux_id, annotated in annotations.items():
        before = (genes.get(flux_id), proteins.get(flux_id))
        if flux_id in proteins:
            proteins[flux_id] = _protein_entry(proteins[flux_id], annotated)
        if flux_id in genes:
            genes[flux_id] = _gene_entry(genes[flux_id], annotated)
        if (genes.get(flux_id), proteins.get(flux_id)) != before:
            changed.append(flux_id)
    return genes, proteins, changed


def _dump_yaml(path: Path, data: Mapping[str, Any]) -> None:
    raw = path.read_bytes()
    text = yaml.safe_dump(dict(data), sort_keys=False, allow_unicode=True, default_flow_style=False)
    newline = "\r\n" if b"\r\n" in raw else "\n"
    path.write_bytes(text.replace("\n", newline).encode("utf-8"))


def sync_role_tables(write: bool = True) -> list[str]:
    """按本表改写两张映射表中登记反应的角色分组；返回有改动的反应。``write=False`` 时只比较。"""
    def read(path: Path) -> dict[str, Any]:
        return yaml.safe_load(path.read_text(encoding="utf-8"))

    genes, proteins, changed = derive_role_tables(read(GENE_TABLE), read(PROTEIN_TABLE))
    if write and changed:
        _dump_yaml(GENE_TABLE, genes)
        _dump_yaml(PROTEIN_TABLE, proteins)
    return changed


# ---------------------------------------------------------------------------
# 审计
# ---------------------------------------------------------------------------

def _table_roles(gene_entry: Mapping[str, Any] | None, protein_entry: Mapping[str, Any] | None
                 ) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    """成员 → 角色集合（统一词汇），分别取自基因表与蛋白表。"""
    gene_roles: dict[str, set[str]] = defaultdict(set)
    for gene_role, values in ((gene_entry or {}).get("gene_roles", {}) or {}).items():
        for member in values or []:
            gene_roles[str(member)].add(GENE_ROLE_ROLE.get(gene_role, gene_role))
    protein_roles: dict[str, set[str]] = defaultdict(set)
    for category in PROTEIN_CATEGORIES:
        for member in (protein_entry or {}).get(category, []) or []:
            protein_roles[str(member)].add(PROTEIN_CATEGORY_ROLE[category])
    return gene_roles, protein_roles


def _row(check: str, flux_id: str, member: str, first: str, second: str, issue: str, omics: str = "") -> tuple:
    # 列与 omicsweave.core.refs.AUDIT_COLUMNS 一致：check, omics, tissue, feature_id, feature_name, data_value, registry_value, issue
    return (check, omics, "", flux_id, member, first, second, issue)


def audit_reaction_roles(gene_map: Mapping[str, Any] | None = None, protein_map: Mapping[str, Any] | None = None,
                         annotations: Mapping[str, ReactionMembers] | None = None) -> list[tuple]:
    """知识表自身的检查，不需要数据。

    - ``reaction_roles / role_mismatch``：未登记的反应中，同一成员在基因表与蛋白表的角色不同（data_value 为基因表，
      registry_value 为蛋白表）；
    - ``reaction_roles / out_of_sync``：登记反应中，表的角色与本表不同（registry_value 为本表）；
    - ``reaction_roles / unannotated_member``：登记反应中，表里有、本表没有的成员；
    - ``reaction_roles / effector_unresolved``：效应物名称查不到代谢物 ID；
    - ``reaction_roles / alias_duplicate``：同一张表中指向同一基因的不同名称（registry_value 为规范符号）；
    - ``regulator_effect / effect_missing``：未登记反应中的调控因子（没有作用方向）。
    """
    from .flux_gene_mapping import FLUX_GENE_MAPPING
    from .flux_protein_mapping import REACTION_ENZYME_MAPPING

    gene_map = FLUX_GENE_MAPPING if gene_map is None else gene_map
    protein_map = REACTION_ENZYME_MAPPING if protein_map is None else protein_map
    annotations = REACTION_MEMBERS if annotations is None else annotations
    rows: list[tuple] = []
    for flux_id in sorted(set(gene_map) | set(protein_map)):
        gene_roles, protein_roles = _table_roles(gene_map.get(flux_id), protein_map.get(flux_id))
        rows += _alias_duplicates(flux_id, gene_map.get(flux_id), protein_map.get(flux_id))
        annotated = annotations.get(flux_id)
        if annotated is None:
            for member in sorted(set(gene_roles) & set(protein_roles)):
                if gene_roles[member] != protein_roles[member]:
                    rows.append(_row("reaction_roles", flux_id, member, "基因表 " + "、".join(sorted(gene_roles[member])),
                                     "蛋白表 " + "、".join(sorted(protein_roles[member])), "role_mismatch"))
            for member in sorted(m for m in set(gene_roles) | set(protein_roles)
                                 if "regulator" in gene_roles.get(m, set()) | protein_roles.get(m, set())):
                rows.append(_row("regulator_effect", flux_id, member, "regulator", "", "effect_missing"))
            continue
        for label, table in (("基因表", gene_roles), ("蛋白表", protein_roles)):
            for member, roles in sorted(table.items()):
                known = annotated.members.get(member)
                if known is None:
                    rows.append(_row("reaction_roles", flux_id, member, f"{label} " + "、".join(sorted(roles)), "",
                                     "unannotated_member"))
                elif roles != {known.role}:
                    rows.append(_row("reaction_roles", flux_id, member, f"{label} " + "、".join(sorted(roles)),
                                     f"本表 {known.role}", "out_of_sync"))
    rows += _unresolved_effectors(annotations)
    return rows


def _alias_duplicates(flux_id: str, gene_entry: Mapping[str, Any] | None,
                      protein_entry: Mapping[str, Any] | None) -> list[tuple]:
    """同一张表中指向同一基因的不同名称（如 Gpi 与 Gpi1）：匹配时合并，但 flux_anchored 的映射数按名称计两次。"""
    from .feature_resolver import _canonical_gene_symbol

    rows = []
    for label, members in (("基因表", (gene_entry or {}).get("genes", []) or []),
                           ("蛋白表", [m for c in PROTEIN_CATEGORIES for m in (protein_entry or {}).get(c, []) or []])):
        names: dict[str, list[str]] = defaultdict(list)
        for member in dict.fromkeys(str(m) for m in members):
            names[_canonical_gene_symbol(member) or member].append(member)
        rows += [_row("reaction_roles", flux_id, "、".join(aliases), label, canonical, "alias_duplicate")
                 for canonical, aliases in sorted(names.items()) if len(aliases) > 1]
    return rows


def _unresolved_effectors(annotations: Mapping[str, ReactionMembers]) -> list[tuple]:
    from omicsweave.core.idmap.metabolite_id_mapping import lookup_metabolite_ids

    return [_row("reaction_roles", flux_id, m.id, "effector", "", "effector_unresolved", omics="metabolomics")
            for flux_id, annotated in sorted(annotations.items()) for m in annotated.by_role("effector")
            if not lookup_metabolite_ids(m.id)]


# ---------------------------------------------------------------------------
# omicsweave refs build reaction_roles [--check]
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="omicsweave refs build reaction_roles",
                                     description="由 reaction_members.yaml 生成基因表与蛋白表中登记反应的角色分组")
    parser.add_argument("--check", action="store_true", help="只比较，不写文件；不同步时退出码为 1")
    args = parser.parse_args(argv or [])
    changed = sync_role_tables(write=not args.check)
    if not changed:
        print(f"两张映射表与 {ANNOTATION_FILE.name} 已同步（登记 {len(REACTION_MEMBERS)} 个反应）")
        return 0
    verb = "需要改写" if args.check else "已改写"
    print(f"{verb}角色分组的反应（{len(changed)} 个）：{'、'.join(changed)}")
    return 1 if args.check else 0


__all__ = ["ROLES", "RELATIONS", "EFFECTS", "MODES", "ACTIVITY_MODES", "Member", "ReactionMembers", "REACTION_MEMBERS",
           "load_reaction_members", "reaction_members", "member_relations", "derive_role_tables", "sync_role_tables",
           "audit_reaction_roles", "main"]
