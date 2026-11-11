#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""小鼠基因注册表：MGI 符号 → Ensembl / NCBI Gene / MGI ID，以及别名解析。

内容由 ``python -m omicsweave.core.idmap.build_registry`` 从 Ensembl GRCm39 cDNA 序列头与
UniProtKB/Swiss-Prot 小鼠条目生成，存放在 ``data/``，来源与版本见 ``REGISTRY_SOURCES``
（``data/SOURCES.json``）。不要手工修改这些表；需要更正时改上游或改 ``build_registry.CURATED_ALIASES`` 后重新生成。

每个 ``GENE_ANNOTATIONS[符号]`` 含：
    - ncbi_gene_id / ncbi_gene_ids : NCBI Gene ID（来自 UniProt 交叉引用；多个时 ncbi_gene_id 取第一个）
    - ensembl_id / ensembl_ids     : Ensembl 基因 ID（不含版本号）
    - mgi_id                       : MGI 登记号
    - gene_name, gene_type, chromosome
    - aliases                      : 解析到该符号的别名

``GENE_ALIAS_LOOKUP`` 把非规范写法映射到 MGI 符号；值为 ``None`` 的条目表示“不是基因”（代谢物、辅因子等）。
只收录无歧义的别名：别名本身若是另一个基因的正式符号，或同时指向多个基因，构建时即被丢弃。

用法::

    from omicsweave.core.idmap.mouse_gene_id_mapping import resolve_gene_symbol, get_ensembl_gene_id
    resolve_gene_symbol('Pc')        # -> 'Pcx'
    get_ensembl_gene_id('Gpd1')      # -> 'ENSMUSG00000023019'
"""

from __future__ import annotations

from typing import Any, Dict, Optional

from ._registry import read_sources, read_table, split_ids

REGISTRY_SOURCES: Dict[str, Any] = read_sources()


def _load() -> tuple[dict, dict, dict, dict]:
    annotations: Dict[str, Dict[str, Any]] = {}
    for row in read_table("mouse_gene_registry.tsv.gz"):
        ensembl = split_ids(row["ensembl_gene_ids"])
        ncbi = split_ids(row["ncbi_gene_ids"])
        annotations[row["symbol"]] = {
            "ncbi_gene_id": ncbi[0] if ncbi else None,
            "ncbi_gene_ids": ncbi,
            "ensembl_id": ensembl[0] if ensembl else None,
            "ensembl_ids": ensembl,
            "mgi_id": row["mgi_id"] or None,
            "gene_name": row["gene_name"] or None,
            "gene_type": row["gene_type"] or None,
            "chromosome": row["chromosome"] or None,
            "aliases": [],
        }
    aliases: Dict[str, Optional[str]] = {}
    for row in read_table("mouse_gene_aliases.tsv.gz"):
        target = row["symbol"] or None
        aliases[row["alias"]] = target
        if target in annotations:
            annotations[target]["aliases"].append(row["alias"])
    ensembl_to: Dict[str, str] = {}
    ncbi_owners: Dict[str, set[str]] = {}
    for symbol, info in annotations.items():
        for gene_id in info["ensembl_ids"]:
            ensembl_to[gene_id] = symbol
        for gene_id in info["ncbi_gene_ids"]:
            ncbi_owners.setdefault(gene_id, set()).add(symbol)
    ncbi_to = {gene_id: next(iter(owners)) for gene_id, owners in ncbi_owners.items() if len(owners) == 1}
    return annotations, aliases, ensembl_to, ncbi_to


GENE_ANNOTATIONS: Dict[str, Dict[str, Any]]
GENE_ALIAS_LOOKUP: Dict[str, Optional[str]]
ENSEMBL_TO_SYMBOL: Dict[str, str]
NCBI_TO_SYMBOL: Dict[str, str]
GENE_ANNOTATIONS, GENE_ALIAS_LOOKUP, ENSEMBL_TO_SYMBOL, NCBI_TO_SYMBOL = _load()


def resolve_gene_symbol(symbol: str) -> Optional[str]:
    """把基因符号解析为 MGI 规范符号。

    正式符号原样返回；别名返回其目标（“不是基因”的条目返回 ``None``）；
    注册表中没有的符号原样返回（视为调用方已给出规范写法）。
    """
    if symbol in GENE_ANNOTATIONS:
        return symbol
    if symbol in GENE_ALIAS_LOOKUP:
        return GENE_ALIAS_LOOKUP[symbol]
    return symbol


def _info(symbol: str) -> Dict[str, Any]:
    canonical = resolve_gene_symbol(symbol)
    if canonical is None:
        return {}
    return GENE_ANNOTATIONS.get(canonical, {})


def get_ncbi_gene_id(symbol: str) -> Optional[str]:
    """NCBI Gene ID（先做别名解析）。"""
    return _info(symbol).get("ncbi_gene_id")


def get_ensembl_gene_id(symbol: str) -> Optional[str]:
    """Ensembl 基因 ID（先做别名解析）。"""
    return _info(symbol).get("ensembl_id")


def get_gene_name(symbol: str) -> Optional[str]:
    """基因全名（先做别名解析）。"""
    return _info(symbol).get("gene_name")


def get_gene_info(symbol: str) -> Dict[str, Any]:
    """完整注释字典；不是基因的条目返回 ``{'ncbi_gene_id': None, 'note': 'non-gene entry'}``，未知符号返回空字典。"""
    if resolve_gene_symbol(symbol) is None:
        return {"ncbi_gene_id": None, "note": "non-gene entry"}
    return _info(symbol)
