#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
“模型包”接口（G7）：证据构建器需要的反应知识——反应 → 基因 / 蛋白 / 代谢物的成员、反应语义（能否给出方向证据）
与脂质池——由一个模型包提供，而不是直接读模块级的映射表。

- :class:`CuratedNetworkPackage`：现有的策展映射与由 MFA 网络生成的反应语义（默认）；
- :class:`DecoyPackage`：把各反应的方向成员换成随机特征的诱饵映射（E3），反应语义沿用原包。

其他网络或基因组尺度模型的基因—反应规则可以实现同一接口后接入，证据构建器不必修改。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Protocol

from . import FLUX_GENE_MAPPING, FLUX_METABOLITE_MAPPING, REACTION_ENZYME_MAPPING, get_metabolite_roles
from .network_annotation import ReactionSemantics, lipid_pool_roles, network_reactions, reaction_semantics
from .reaction_members import member_relations


class ModelPackage(Protocol):
    id: str

    def gene_mapping(self, flux_id: str) -> dict[str, Any]: ...

    def enzyme_mapping(self, flux_id: str) -> dict[str, Any]: ...

    def member_relations(self, flux_id: str) -> dict[str, str]:
        """直接酶 / 转运体之间的关系（isozyme、complex_subunit、shared_subunit）；没有登记的成员按同工酶处理。"""
        ...

    def metabolite_roles(self, flux_id: str) -> dict[str, list[str]] | list[str] | None: ...

    def semantics(self, flux_id: str) -> ReactionSemantics: ...

    def lipid_pools(self, flux_id: str) -> dict[str, set[str]]: ...


@dataclass(frozen=True)
class CuratedNetworkPackage:
    id: str = "curated_mfa_network"

    def gene_mapping(self, flux_id: str) -> dict[str, Any]:
        return FLUX_GENE_MAPPING.get(flux_id, {})

    def enzyme_mapping(self, flux_id: str) -> dict[str, Any]:
        return REACTION_ENZYME_MAPPING.get(flux_id, {})

    def member_relations(self, flux_id: str) -> dict[str, str]:
        return member_relations(flux_id)

    def metabolite_roles(self, flux_id: str):
        roles = get_metabolite_roles(flux_id)
        return roles if roles is not None else FLUX_METABOLITE_MAPPING.get(flux_id)

    def semantics(self, flux_id: str) -> ReactionSemantics:
        return reaction_semantics(flux_id)

    def lipid_pools(self, flux_id: str) -> dict[str, set[str]]:
        return lipid_pool_roles(flux_id)

    def reactions(self):
        return network_reactions()


@dataclass
class DecoyPackage:
    """诱饵映射：``genes`` / ``proteins`` / ``metabolites`` 为 {反应: 诱饵成员}。代谢物诱饵保留底物 / 产物角色
    （{反应: {"substrates": [...], "products": [...]}}）。没有给出诱饵的层视为没有映射。"""
    base: ModelPackage
    genes: dict[str, list[str]] = field(default_factory=dict)
    proteins: dict[str, list[str]] = field(default_factory=dict)
    metabolites: dict[str, dict[str, list[str]]] = field(default_factory=dict)
    id: str = "decoy"

    def gene_mapping(self, flux_id: str) -> dict[str, Any]:
        members = self.genes.get(flux_id, [])
        return {"genes": members, "gene_roles": {"catalytic": members}} if members else {}

    def enzyme_mapping(self, flux_id: str) -> dict[str, Any]:
        members = self.proteins.get(flux_id, [])
        return {"direct_enzymes": members} if members else {}

    def member_relations(self, flux_id: str) -> dict[str, str]:
        return {}                          # 诱饵成员互不相关，按同工酶处理

    def metabolite_roles(self, flux_id: str):
        roles = self.metabolites.get(flux_id)
        return {"substrates": roles.get("substrates", []), "products": roles.get("products", []), "related": []} \
            if roles else None

    def semantics(self, flux_id: str) -> ReactionSemantics:
        return self.base.semantics(flux_id)

    def lipid_pools(self, flux_id: str) -> dict[str, set[str]]:
        return {"substrates": set(), "products": set()}


DEFAULT_PACKAGE = CuratedNetworkPackage()
