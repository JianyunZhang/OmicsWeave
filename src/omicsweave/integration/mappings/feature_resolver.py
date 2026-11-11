#!/usr/bin/env python3
"""跨组学特征标识符解析。

该模块是 scorer 与各类手工映射表之间的统一 seam。调用方只提交数据表和
目标映射；名称规范化、标准 ID 解析和歧义隔离都隐藏在模块内部。
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Any, Mapping, Sequence
import re
import unicodedata
import weakref

import pandas as pd

from omicsweave.core.idmap.mouse_gene_id_mapping import (
    ENSEMBL_TO_SYMBOL,
    GENE_ALIAS_LOOKUP,
    GENE_ANNOTATIONS,
    NCBI_TO_SYMBOL,
)
from omicsweave.core.idmap.mouse_protein_id_mapping import (
    MOUSE_PROTEIN_UNIPROT_MAPPING,
    UNIPROT_TO_SYMBOL,
)
from omicsweave.core.idmap.metabolite_id_mapping import lookup_metabolite_ids


METABOLITE_ROLES = ("substrates", "products", "related")
METABOLITE_ID_COLUMNS = (
    "hmdb_id", "HMDB",
    "kegg_compound", "KEGG_Compound",
    "chebi_id", "ChEBI",
)


@dataclass(frozen=True)
class MetaboliteRoleMatch:
    """一次代谢物角色匹配的不可变结果。"""

    indices_by_role: dict[str, tuple[Any, ...]]
    ambiguous_indices: tuple[Any, ...]
    match_method_by_index: dict[Any, str]
    target_count_by_role: dict[str, int]


@dataclass(frozen=True)
class IdentifierMatch:
    """基因或蛋白标识符匹配结果。"""

    indices: tuple[Any, ...]
    matched_queries: tuple[str, ...]
    match_method_by_index: dict[Any, str]


def _unique_casefold(pairs: Sequence[tuple[str, str | None]]) -> dict[str, str | None]:
    """忽略大小写的查找表；大小写折叠后指向不同目标的键不收录（只能精确匹配）。"""
    targets: dict[str, set[str | None]] = defaultdict(set)
    for key, target in pairs:
        targets[key.casefold()].add(target)
    return {key: next(iter(values)) for key, values in targets.items() if len(values) == 1}


_GENE_SYMBOL_LOWER = _unique_casefold([(symbol, symbol) for symbol in GENE_ANNOTATIONS])
_GENE_ALIAS_LOWER = _unique_casefold(list(GENE_ALIAS_LOOKUP.items()))
_ENSEMBL_TO_GENE = {gene_id.upper(): symbol for gene_id, symbol in ENSEMBL_TO_SYMBOL.items()}
_NCBI_TO_GENE = dict(NCBI_TO_SYMBOL)
_GENE_FEATURE_CACHE: dict[int, tuple[weakref.ReferenceType, Any]] = {}
_PROTEIN_FEATURE_CACHE: dict[int, tuple[weakref.ReferenceType, Any]] = {}
_METABOLITE_FEATURE_CACHE: dict[int, tuple[weakref.ReferenceType, Any]] = {}


def _get_cached(
    cache: dict[int, tuple[weakref.ReferenceType, Any]],
    dataframe: pd.DataFrame,
) -> Any | None:
    record = cache.get(id(dataframe))
    if record is not None and record[0]() is dataframe:
        return record[1]
    return None


def _set_cached(
    cache: dict[int, tuple[weakref.ReferenceType, Any]],
    dataframe: pd.DataFrame,
    value: Any,
) -> None:
    key = id(dataframe)

    def _cleanup(reference: weakref.ReferenceType) -> None:
        current = cache.get(key)
        if current is not None and current[0] is reference:
            cache.pop(key, None)

    cache[key] = (weakref.ref(dataframe, _cleanup), value)


def _normalize_name(value: Any) -> str:
    """规范化名称，但不删除具有生物学含义的立体化学前缀。"""
    if value is None or pd.isna(value):
        return ""
    text = unicodedata.normalize("NFKC", str(value)).strip().casefold()
    text = text.replace("–", "-").replace("—", "-").replace("−", "-")
    return re.sub(r"\s+", " ", text)


def _normalize_database_id(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    text = unicodedata.normalize("NFKC", str(value)).strip().upper()
    if not text or text in {"-", "NA", "N/A", "NONE", "NAN", "NULL"}:
        return ""
    hmdb = re.fullmatch(r"HMDB0*(\d+)", text)
    if hmdb:
        return f"HMDB{int(hmdb.group(1)):07d}"
    kegg = re.fullmatch(r"C(\d{5})", text)
    if kegg:
        return f"C{kegg.group(1)}"
    chebi = re.fullmatch(r"(?:CHEBI:?)?(\d+)", text)
    if chebi and (text.startswith("CHEBI") or text.isdigit()):
        return f"CHEBI:{int(chebi.group(1))}"
    return ""


def _database_ids_from_value(value: Any) -> set[str]:
    if value is None or pd.isna(value):
        return set()
    return {
        normalized
        for token in re.split(r"[;|,]", str(value))
        if (normalized := _normalize_database_id(token))
    }


def _metabolite_target_key(alias: Any) -> str:
    identifiers = lookup_metabolite_ids(str(alias))
    if identifiers:
        for field in ("hmdb_id", "kegg_compound", "chebi_id"):
            normalized = _normalize_database_id(identifiers.get(field))
            if normalized:
                return normalized
    return "name:" + _normalize_name(alias)


def _canonical_gene_symbol(value: Any) -> str | None:
    text = "" if value is None or (not isinstance(value, str) and pd.isna(value)) else str(value).strip()
    if text in GENE_ANNOTATIONS:
        return text
    if text in GENE_ALIAS_LOOKUP:
        return GENE_ALIAS_LOOKUP[text]
    normalized = _normalize_name(value)
    if not normalized:
        return None
    canonical = _GENE_SYMBOL_LOWER.get(normalized)
    if canonical:
        return canonical
    alias_target = _GENE_ALIAS_LOWER.get(normalized)
    return alias_target or None


def _normalize_ncbi_gene_id(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    if re.fullmatch(r"\d+\.0", text):
        text = text[:-2]
    return text


def _gene_feature_index(dataframe: pd.DataFrame) -> dict[str, tuple[tuple[Any, str], ...]]:
    cached = _get_cached(_GENE_FEATURE_CACHE, dataframe)
    if cached is not None:
        return cached

    symbol_columns = [
        column for column in ("gene", "gene_symbol", "canonical_symbol")
        if column in dataframe.columns
    ]
    ensembl_columns = [
        column for column in ("ensembl_gene_id", "ensembl_id")
        if column in dataframe.columns
    ]
    ncbi_columns = [
        column for column in ("ncbi_gene_id", "gene_id")
        if column in dataframe.columns
    ]
    index: dict[str, list[tuple[Any, str]]] = defaultdict(list)
    for row_index, row in dataframe.iterrows():
        row_gene: str | None = None
        method = ""
        for column in ensembl_columns:
            value = row.get(column)
            if value is not None and not pd.isna(value):
                row_gene = _ENSEMBL_TO_GENE.get(str(value).strip().upper())
                if row_gene:
                    method = "ensembl_gene_id"
                    break
        if row_gene is None:
            for column in ncbi_columns:
                row_gene = _NCBI_TO_GENE.get(_normalize_ncbi_gene_id(row.get(column)))
                if row_gene:
                    method = "ncbi_gene_id"
                    break
        if row_gene is None:
            for column in symbol_columns:
                row_gene = _canonical_gene_symbol(row.get(column))
                if row_gene:
                    method = "canonical_symbol"
                    break
        if row_gene:
            index[row_gene].append((row_index, method))
    frozen = {gene: tuple(rows) for gene, rows in index.items()}
    _set_cached(_GENE_FEATURE_CACHE, dataframe, frozen)
    return frozen


def match_gene_features(
    dataframe: pd.DataFrame,
    genes: Sequence[str],
) -> IdentifierMatch:
    """用规范小鼠基因符号、Ensembl 或 NCBI Gene ID 精确匹配 mRNA。"""
    canonical_queries = tuple(dict.fromkeys(
        canonical
        for gene in genes
        if (canonical := _canonical_gene_symbol(gene))
    ))
    query_set = set(canonical_queries)
    if dataframe.empty or not query_set:
        return IdentifierMatch((), (), {})

    feature_index = _gene_feature_index(dataframe)
    selected: dict[Any, tuple[str, str]] = {}
    for query in canonical_queries:
        for row_index, method in feature_index.get(query, ()):
            selected[row_index] = (query, method)
    matched_indices = [index for index in dataframe.index if index in selected]
    methods: dict[Any, str] = {}
    matched_canonical: set[str] = set()
    for row_index in matched_indices:
        row_gene, method = selected[row_index]
        methods[row_index] = method
        matched_canonical.add(row_gene)

    matched_queries = tuple(
        query for query in canonical_queries if query in matched_canonical
    )
    return IdentifierMatch(tuple(matched_indices), matched_queries, methods)


def _uniprot_ids_from_value(value: Any) -> set[str]:
    if value is None or pd.isna(value):
        return set()
    pieces = re.split(r"[;,|\s]+", str(value).strip().upper())
    return {piece for piece in pieces if piece and piece not in {"SP", "TR"}}


def _protein_feature_index(dataframe: pd.DataFrame) -> dict[str, Any]:
    cached = _get_cached(_PROTEIN_FEATURE_CACHE, dataframe)
    if cached is not None:
        return cached

    protein_columns = [
        column for column in ("protein", "Protein", "gene", "gene_symbol")
        if column in dataframe.columns
    ]
    uniprot_columns = [
        column for column in ("uniprot_id", "UniProt", "accession")
        if column in dataframe.columns
    ]
    accession_index: dict[str, list[Any]] = defaultdict(list)
    symbol_index: dict[str, list[tuple[Any, bool]]] = defaultdict(list)
    for row_index, row in dataframe.iterrows():
        row_accessions: set[str] = set()
        for column in uniprot_columns:
            row_accessions.update(_uniprot_ids_from_value(row.get(column)))
        for accession in row_accessions:
            accession_index[accession].append(row_index)
        canonical = next(
            (
                resolved for column in protein_columns
                if (resolved := _canonical_gene_symbol(row.get(column)))
            ),
            None,
        )
        if canonical:
            # None：数据行没有登记号；空集：登记号都不在注册表中；否则为注册表中这些登记号所属的基因
            owners = (
                {UNIPROT_TO_SYMBOL[a] for a in row_accessions if a in UNIPROT_TO_SYMBOL}
                if row_accessions else None
            )
            symbol_index[canonical].append((row_index, owners))
    frozen = {
        "accessions": {
            accession: tuple(rows) for accession, rows in accession_index.items()
        },
        "symbols": {
            symbol: tuple(rows) for symbol, rows in symbol_index.items()
        },
    }
    _set_cached(_PROTEIN_FEATURE_CACHE, dataframe, frozen)
    return frozen


def match_protein_features(
    dataframe: pd.DataFrame,
    proteins: Sequence[str],
) -> IdentifierMatch:
    """用 UniProt 登记号优先、规范小鼠基因符号回退匹配蛋白。

    数据行的登记号与注册表一致时按登记号匹配。按符号回退只在数据行没有可信冲突时进行：
    数据行若带有注册表登记在另一个基因名下的 Swiss-Prot 登记号，说明符号与登记号矛盾，不按符号匹配；
    数据行的登记号不在注册表中（TrEMBL、亚型、厂商编号）时允许按符号回退。
    """
    canonical_queries = tuple(dict.fromkeys(
        canonical
        for protein in proteins
        if (canonical := _canonical_gene_symbol(protein))
    ))
    query_set = set(canonical_queries)
    if dataframe.empty or not query_set:
        return IdentifierMatch((), (), {})

    query_by_uniprot: dict[str, set[str]] = {}
    for query in canonical_queries:
        info = MOUSE_PROTEIN_UNIPROT_MAPPING.get(query, {})
        for accession in info.get("uniprot_ids") or (info.get("uniprot_id"),):
            accession = str(accession or "").strip().upper()
            if accession:
                query_by_uniprot.setdefault(accession, set()).add(query)

    feature_index = _protein_feature_index(dataframe)
    selected_queries: dict[Any, set[str]] = defaultdict(set)
    methods: dict[Any, str] = {}
    for accession, queries in query_by_uniprot.items():
        for row_index in feature_index["accessions"].get(accession, ()):
            selected_queries[row_index].update(queries)
            methods[row_index] = "uniprot_id"
    for query in canonical_queries:
        for row_index, row_owners in feature_index["symbols"].get(query, ()):
            if row_owners is None:
                method = "canonical_symbol"
            elif not row_owners:
                method = "canonical_symbol_unregistered_accession"
            elif query in row_owners:
                method = "uniprot_id"
            else:
                continue  # 登记号属于另一个基因：符号与登记号矛盾
            selected_queries[row_index].add(query)
            methods.setdefault(row_index, method)

    matched_indices = [
        row_index for row_index in dataframe.index
        if len(selected_queries.get(row_index, set())) == 1
    ]
    matched_canonical = {
        next(iter(selected_queries[row_index])) for row_index in matched_indices
    }
    methods = {row_index: methods[row_index] for row_index in matched_indices}

    matched_queries = tuple(
        query for query in canonical_queries if query in matched_canonical
    )
    return IdentifierMatch(tuple(matched_indices), matched_queries, methods)


def _metabolite_feature_index(
    dataframe: pd.DataFrame,
    feature_column: str,
) -> dict[str, Any]:
    cache = _get_cached(_METABOLITE_FEATURE_CACHE, dataframe)
    if cache is None:
        cache = {}
        _set_cached(_METABOLITE_FEATURE_CACHE, dataframe, cache)
    if feature_column in cache:
        return cache[feature_column]

    available_id_columns = [
        column for column in METABOLITE_ID_COLUMNS if column in dataframe.columns
    ]
    id_index: dict[str, list[Any]] = defaultdict(list)
    name_index: dict[str, list[Any]] = defaultdict(list)
    for row_index, value in dataframe[feature_column].items():
        row_ids: set[str] = set()
        for column in available_id_columns:
            row_ids.update(_database_ids_from_value(dataframe.at[row_index, column]))
        if row_ids:
            for identifier in row_ids:
                id_index[identifier].append(row_index)
        else:
            normalized = _normalize_name(value)
            if normalized:
                name_index[normalized].append(row_index)
    prepared = {
        "ids": {identifier: tuple(rows) for identifier, rows in id_index.items()},
        "names": {name: tuple(rows) for name, rows in name_index.items()},
    }
    cache[feature_column] = prepared
    return prepared


def match_metabolite_roles(
    dataframe: pd.DataFrame,
    roles: Mapping[str, Sequence[str]],
    *,
    feature_column: str = "compound",
) -> MetaboliteRoleMatch:
    """按规范化后的完整名称匹配底物、产物和关联代谢物。

    不执行任意子串搜索。同一数据行若命中多个角色，会从所有角色结果中移除并
    放入 ``ambiguous_indices``，避免同时被解释为底物和产物。
    """
    empty = {role: () for role in METABOLITE_ROLES}
    if feature_column not in dataframe.columns or dataframe.empty:
        counts = {
            role: len({
                key for alias in roles.get(role, ())
                if (key := _metabolite_target_key(alias)) != "name:"
            })
            for role in METABOLITE_ROLES
        }
        return MetaboliteRoleMatch(empty, (), {}, counts)

    aliases = {
        role: {
            normalized
            for value in roles.get(role, ())
            if (normalized := _normalize_name(value))
        }
        for role in METABOLITE_ROLES
    }
    target_count_by_role = {
        role: len({
            key for alias in roles.get(role, ())
            if (key := _metabolite_target_key(alias)) != "name:"
        })
        for role in METABOLITE_ROLES
    }
    role_database_ids: dict[str, set[str]] = {
        role: set() for role in METABOLITE_ROLES
    }
    for role in METABOLITE_ROLES:
        for alias in roles.get(role, ()):
            identifiers = lookup_metabolite_ids(str(alias))
            if identifiers:
                role_database_ids[role].update(
                    normalized
                    for value in identifiers.values()
                    if (normalized := _normalize_database_id(value))
                )

    feature_index = _metabolite_feature_index(dataframe, feature_column)
    role_hits: dict[str, set[Any]] = {
        role: set() for role in METABOLITE_ROLES
    }
    method_by_index: dict[Any, str] = {}
    for role in METABOLITE_ROLES:
        for identifier in role_database_ids[role]:
            for row_index in feature_index["ids"].get(identifier, ()):
                role_hits[role].add(row_index)
                method_by_index[row_index] = "database_id"
        for alias in aliases[role]:
            for row_index in feature_index["names"].get(alias, ()):
                role_hits[role].add(row_index)
                method_by_index.setdefault(row_index, "exact_name")

    hits_by_index = {
        row_index: [
            role for role in METABOLITE_ROLES if row_index in role_hits[role]
        ]
        for row_index in dataframe.index
        if any(row_index in role_hits[role] for role in METABOLITE_ROLES)
    }

    ambiguous = tuple(
        row_index for row_index, matched_roles in hits_by_index.items()
        if len(matched_roles) > 1
    )
    ambiguous_set = set(ambiguous)
    indices_by_role = {
        role: tuple(
            row_index for row_index, matched_roles in hits_by_index.items()
            if row_index not in ambiguous_set and role in matched_roles
        )
        for role in METABOLITE_ROLES
    }
    methods = {
        row_index: method_by_index[row_index]
        for row_index in hits_by_index
        if row_index not in ambiguous_set
    }
    return MetaboliteRoleMatch(
        indices_by_role,
        ambiguous,
        methods,
        target_count_by_role,
    )


__all__ = [
    "IdentifierMatch",
    "MetaboliteRoleMatch",
    "match_gene_features",
    "match_metabolite_roles",
    "match_protein_features",
]
