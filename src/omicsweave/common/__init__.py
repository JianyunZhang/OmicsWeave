#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
omicsweave.common — 跨组学工具共享模块。

提供所有组学分析工具（代谢组学、蛋白质组学、转录组学、代谢通量分析、多组学整合）
共用的标准化数据库 ID 映射字典和工具函数，避免代码重复，便于统一维护。

模块:
    metabolite_id_mapping      — 代谢物 → HMDB / KEGG Compound / ChEBI
    protein_id_mapping         — 基因符号 → UniProt / EC Number
    gene_id_mapping            — 基因符号 → NCBI Gene ID / Ensembl ID（含别名解析）
    reaction_pathway_mapping   — MFA 反应 ID → KEGG Reaction / Pathway / Reactome / EC
"""

from .metabolite_id_mapping import (
    METABOLITE_DATABASE_IDS,
    lookup_metabolite_ids,
    annotate_metabolites,
    get_hmdb_id,
    get_kegg_compound,
    get_chebi_id,
)
from .mouse_protein_id_mapping import (
    MOUSE_PROTEIN_UNIPROT_MAPPING,
    get_uniprot_id,
    get_ec_for_protein,
)
from .mouse_gene_id_mapping import (
    GENE_ANNOTATIONS,
    GENE_ALIAS_LOOKUP,
    resolve_gene_symbol,
    get_ncbi_gene_id,
    get_ensembl_gene_id,
    get_gene_name,
    get_gene_info,
)
from .reaction_pathway_mapping import (
    REACTION_REGISTRY,
    ReactionInfo,
    get_reaction_info,
    get_kegg_reaction,
    get_kegg_pathway,
    get_reactome_id,
    get_ec_number,
    get_description,
    list_reactions_by_pathway,
    list_reactions_with_ec,
    get_all_pathways,
    annotate_flux_dataframe,
)
from .meta_analysis import MetaAnalysisResult, random_effects_meta
from .stability import StabilityResult, calibrate_biological_stability

__all__ = [
    # 代谢物
    "METABOLITE_DATABASE_IDS",
    "lookup_metabolite_ids",
    "annotate_metabolites",
    "get_hmdb_id",
    "get_kegg_compound",
    "get_chebi_id",
    # 蛋白质
    "MOUSE_PROTEIN_UNIPROT_MAPPING",
    "get_uniprot_id",
    "get_ec_for_protein",
    # 基因
    "GENE_ANNOTATIONS",
    "GENE_ALIAS_LOOKUP",
    "resolve_gene_symbol",
    "get_ncbi_gene_id",
    "get_ensembl_gene_id",
    "get_gene_name",
    "get_gene_info",
    # 代谢反应通路
    "REACTION_REGISTRY",
    "ReactionInfo",
    "get_reaction_info",
    "get_kegg_reaction",
    "get_kegg_pathway",
    "get_reactome_id",
    "get_ec_number",
    "get_description",
    "list_reactions_by_pathway",
    "list_reactions_with_ec",
    "get_all_pathways",
    "annotate_flux_dataframe",
    # 统计整合
    "MetaAnalysisResult",
    "random_effects_meta",
    "StabilityResult",
    "calibrate_biological_stability",
]
