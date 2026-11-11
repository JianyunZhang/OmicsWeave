#!/usr/bin/env python3
"""多组学反应映射的一致性审计。"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Any

import pandas as pd

from omicsweave.core.idmap.metabolite_id_mapping import lookup_metabolite_ids
from omicsweave.core.idmap.mouse_gene_id_mapping import (
    GENE_ALIAS_LOOKUP,
    GENE_ANNOTATIONS,
)
from omicsweave.core.idmap.mouse_protein_id_mapping import (
    MOUSE_PROTEIN_UNIPROT_MAPPING,
    UNIPROT_TO_SYMBOL,
)
from omicsweave.core.idmap.reaction_pathway_mapping import REACTION_REGISTRY

from .flux_gene_mapping import FLUX_GENE_MAPPING
from .flux_metabolite_mapping import FLUX_METABOLITE_MAPPING
from .flux_protein_mapping import REACTION_ENZYME_MAPPING


PROTEIN_ROLES = (
    "direct_enzymes",
    "upstream_regulators",
    "downstream_enzymes",
    "transporters",
    "cofactors",
    "related_isozymes",
)
METABOLITE_ROLES = ("substrates", "products", "related")


@dataclass(frozen=True)
class MappingIntegrityAudit:
    """映射完整性检查结果。errors 可用于阻断 pipeline。"""

    errors: tuple[str, ...]
    warnings: tuple[str, ...]
    metrics: dict[str, int]

    def raise_for_errors(self) -> None:
        if self.errors:
            raise ValueError("多组学映射完整性检查失败:\n- " + "\n- ".join(self.errors))


def _canonical_gene(symbol: str) -> str | None:
    from .feature_resolver import _canonical_gene_symbol

    return _canonical_gene_symbol(symbol)


def _database_ids(name: str) -> set[str]:
    identifiers = lookup_metabolite_ids(name)
    if not identifiers:
        return set()
    return {str(value) for value in identifiers.values() if value}


def audit_mapping_integrity() -> MappingIntegrityAudit:
    """检查反应键、标识符、角色歧义和注册表内部冲突。"""
    errors: list[str] = []
    warnings: list[str] = []
    metrics: Counter[str] = Counter()

    mappings: dict[str, dict[str, Any]] = {
        "gene": FLUX_GENE_MAPPING,
        "protein": REACTION_ENZYME_MAPPING,
        "metabolite": FLUX_METABOLITE_MAPPING,
    }
    for layer, mapping in mappings.items():
        orphaned = sorted(set(mapping) - set(REACTION_REGISTRY))
        if orphaned:
            errors.append(f"{layer} 映射含未注册反应: {', '.join(orphaned)}")
        metrics[f"{layer}_reaction_count"] = len(mapping)

    gene_unknown: set[str] = set()
    gene_duplicate_entries = 0
    for reaction, entry in FLUX_GENE_MAPPING.items():
        genes = list(entry.get("genes", []))
        gene_duplicate_entries += len(genes) - len(set(genes))
        for gene in genes:
            if _canonical_gene(str(gene)) is None:
                gene_unknown.add(str(gene))
                warnings.append(f"{reaction}: 未注册基因符号 {gene}")
    metrics["gene_unknown_symbols"] = len(gene_unknown)
    metrics["gene_duplicate_entries"] = gene_duplicate_entries

    protein_unknown: set[str] = set()
    protein_missing_uniprot: set[str] = set()
    protein_non_gene: set[str] = set()
    for reaction, entry in REACTION_ENZYME_MAPPING.items():
        for role in PROTEIN_ROLES:
            values = list(entry.get(role, []))
            metrics["protein_duplicate_entries"] += len(values) - len(set(values))
            for protein in values:
                canonical = _canonical_gene(str(protein))
                if canonical is None:
                    if protein in GENE_ALIAS_LOOKUP and GENE_ALIAS_LOOKUP[protein] is None:
                        protein_non_gene.add(str(protein))
                        errors.append(f"{reaction}.{role}: 非蛋白实体 {protein}")
                    else:
                        protein_unknown.add(str(protein))
                    continue
                if canonical not in MOUSE_PROTEIN_UNIPROT_MAPPING:
                    protein_missing_uniprot.add(canonical)
    metrics["protein_unknown_symbols"] = len(protein_unknown)
    metrics["protein_symbols_without_uniprot"] = len(protein_missing_uniprot)
    metrics["protein_non_gene_entities"] = len(protein_non_gene)
    if protein_unknown:
        warnings.append(
            "蛋白映射中有未纳入本地基因注册表的符号: "
            + ", ".join(sorted(protein_unknown))
        )
    if protein_missing_uniprot:
        warnings.append(
            f"{len(protein_missing_uniprot)} 个规范小鼠蛋白尚无本地 UniProt 注释"
        )

    accessions: dict[str, set[str]] = defaultdict(set)
    for symbol, info in MOUSE_PROTEIN_UNIPROT_MAPPING.items():
        accession = str(info.get("uniprot_id") or "").strip().upper()
        canonical = _canonical_gene(symbol) or symbol
        if accession:
            accessions[accession].add(canonical)
    accession_conflicts = {
        accession: symbols
        for accession, symbols in accessions.items()
        if len(symbols) > 1
    }
    metrics["uniprot_accession_conflicts"] = len(accession_conflicts)
    for accession, symbols in sorted(accession_conflicts.items()):
        errors.append(
            f"UniProt {accession} 同时映射到不同基因: {', '.join(sorted(symbols))}"
        )

    metabolite_aliases: set[str] = set()
    metabolite_without_ids: set[str] = set()
    metabolite_duplicate_entries = 0
    for reaction, entry in FLUX_METABOLITE_MAPPING.items():
        ids_by_role: dict[str, set[str]] = {}
        for role in METABOLITE_ROLES:
            aliases = [str(value).strip() for value in entry.get(role, [])]
            if any(not alias for alias in aliases):
                errors.append(f"{reaction}.{role}: 含空代谢物名称")
            metabolite_duplicate_entries += len(aliases) - len(set(aliases))
            metabolite_aliases.update(aliases)
            role_ids: set[str] = set()
            for alias in aliases:
                identifiers = _database_ids(alias)
                if identifiers:
                    role_ids.update(identifiers)
                else:
                    metabolite_without_ids.add(alias)
            ids_by_role[role] = role_ids
        for index, first in enumerate(METABOLITE_ROLES):
            for second in METABOLITE_ROLES[index + 1:]:
                overlap = ids_by_role[first] & ids_by_role[second]
                if overlap:
                    errors.append(
                        f"{reaction}: 标准代谢物 ID 同时属于 {first}/{second}: "
                        + ", ".join(sorted(overlap))
                    )
    metrics["metabolite_unique_aliases"] = len(metabolite_aliases)
    metrics["metabolite_aliases_without_standard_id"] = len(metabolite_without_ids)
    metrics["metabolite_duplicate_entries"] = metabolite_duplicate_entries
    if metabolite_without_ids:
        warnings.append(
            f"{len(metabolite_without_ids)} 个代谢物别名仅可进行规范化精确名称匹配"
        )

    metrics["all_three_reaction_count"] = len(
        set(FLUX_GENE_MAPPING)
        & set(REACTION_ENZYME_MAPPING)
        & set(FLUX_METABOLITE_MAPPING)
    )
    return MappingIntegrityAudit(
        errors=tuple(dict.fromkeys(errors)),
        warnings=tuple(dict.fromkeys(warnings)),
        metrics=dict(metrics),
    )


@dataclass(frozen=True)
class RegistryDataAudit:
    """注册表与一份蛋白组数据自带登记号的对照。

    只统计“符号能解析为规范基因、且登记号在注册表中”的行；``conflicts`` 为数据行的符号与
    注册表中该登记号所属基因不一致的行（数据行符号, 登记号, 注册表中的基因）。
    """

    n_rows: int
    n_checked: int
    n_agree: int
    conflicts: tuple[tuple[str, str, str], ...]
    n_unregistered_accession: int

    @property
    def agreement(self) -> float:
        return self.n_agree / self.n_checked if self.n_checked else float("nan")


def audit_registry_against_data(
    proteins: pd.DataFrame,
    symbol_column: str = "feature_id",
    accession_column: str = "uniprot_id",
) -> RegistryDataAudit:
    """检查注册表的 UniProt 登记号是否与蛋白组数据自带的“基因符号—登记号”一致。

    蛋白组搜库结果中的登记号与基因名来自同一份 FASTA，是最直接的对照；不一致说明注册表或数据标注有误。
    蛋白组（``;`` 分隔的多个蛋白）只看第一个成员。
    """
    from .feature_resolver import _canonical_gene_symbol

    pairs = proteins[[symbol_column, accession_column]].dropna().drop_duplicates()
    n_checked = n_agree = n_unregistered = 0
    conflicts: list[tuple[str, str, str]] = []
    for symbol_value, accession_value in pairs.itertuples(index=False):
        symbol = str(symbol_value).split(";")[0].strip()
        accession = str(accession_value).split(";")[0].strip().upper()
        canonical = _canonical_gene_symbol(symbol)
        if not canonical or not accession:
            continue
        owner = UNIPROT_TO_SYMBOL.get(accession)
        if owner is None:
            n_unregistered += 1
            continue
        n_checked += 1
        if owner == canonical:
            n_agree += 1
        else:
            conflicts.append((symbol, accession, owner))
    return RegistryDataAudit(len(pairs), n_checked, n_agree, tuple(sorted(conflicts)), n_unregistered)


__all__ = ["MappingIntegrityAudit", "audit_mapping_integrity", "RegistryDataAudit", "audit_registry_against_data"]
