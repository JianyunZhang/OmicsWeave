#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
omicsweave.core — 跨组学共享的基础模块。

提供所有组学分析工具（代谢组学、蛋白质组学、转录组学、代谢通量分析、多组学整合）
共用的标准化数据库 ID 映射字典和工具函数，避免代码重复，便于统一维护。

子模块:
    project                          — 项目目录定位
    safety                           — pandas 数值安全检查（numexpr 静默算错的防护）
    statistics.nonparametric         — 秩检验在给定样本量下可达到的最小 p 值 / FDR
    statistics.meta_analysis         — 随机效应 Meta 分析
    statistics.stability             — 生物学稳定性校准
    idmap.metabolite_id_mapping      — 代谢物 → HMDB / KEGG Compound / ChEBI
    idmap.mouse_protein_id_mapping   — 基因符号 → UniProt / EC Number
    idmap.mouse_gene_id_mapping      — 基因符号 → NCBI Gene ID / Ensembl ID（含别名解析）
    idmap.reaction_pathway_mapping   — MFA 反应 ID → KEGG Reaction / Pathway / Reactome / EC
"""

from importlib import import_module


# 轻量 CLI 与后台运行时按需导入统计和 ID 注释模块。
_MODULE_EXPORTS = {
    'idmap.metabolite_id_mapping': (
        'METABOLITE_DATABASE_IDS',
        'lookup_metabolite_ids',
        'annotate_metabolites',
        'get_hmdb_id',
        'get_kegg_compound',
        'get_chebi_id',
    ),
    'idmap.mouse_protein_id_mapping': (
        'MOUSE_PROTEIN_UNIPROT_MAPPING',
        'get_uniprot_id',
        'get_ec_for_protein',
    ),
    'idmap.mouse_gene_id_mapping': (
        'GENE_ANNOTATIONS',
        'GENE_ALIAS_LOOKUP',
        'resolve_gene_symbol',
        'get_ncbi_gene_id',
        'get_ensembl_gene_id',
        'get_gene_name',
        'get_gene_info',
    ),
    'idmap.reaction_pathway_mapping': (
        'REACTION_REGISTRY',
        'ReactionInfo',
        'get_reaction_info',
        'get_kegg_reaction',
        'get_kegg_pathway',
        'get_reactome_id',
        'get_ec_number',
        'get_description',
        'list_reactions_by_pathway',
        'list_reactions_with_ec',
        'get_all_pathways',
        'annotate_flux_dataframe',
    ),
    'statistics.meta_analysis': (
        'MetaAnalysisResult',
        'random_effects_meta',
    ),
    'statistics.stability': (
        'StabilityResult',
        'calibrate_biological_stability',
    ),
}
_EXPORTS = {name: module for module, names in _MODULE_EXPORTS.items() for name in names}
__all__ = list(_EXPORTS)


def __getattr__(name: str):
    if name not in _EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(import_module(f".{_EXPORTS[name]}", __name__), name)
    globals()[name] = value
    return value


def __dir__():
    return sorted(set(globals()).union(__all__))
