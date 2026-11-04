"""
KEGG 通路名称映射字典及查询函数。

从 run_metabolomics_analysis.py 中提取，独立为可复用模块。
"""
from __future__ import annotations

import re

KEGG_PATHWAYS: dict[str, str] = {
    'ko00010': 'Glycolysis / Gluconeogenesis',
    'ko00020': 'Citrate cycle (TCA cycle)',
    'ko00030': 'Pentose phosphate pathway',
    'ko00040': 'Pentose and glucuronate interconversions',
    'ko00051': 'Fructose and mannose metabolism',
    'ko00052': 'Galactose metabolism',
    'ko00053': 'Ascorbate and aldarate metabolism',
    'ko00061': 'Fatty acid biosynthesis',
    'ko00062': 'Fatty acid elongation',
    'ko00071': 'Fatty acid degradation',
    'ko00072': 'Synthesis and degradation of ketone bodies',
    'ko00100': 'Steroid biosynthesis',
    'ko00120': 'Primary bile acid biosynthesis',
    'ko00130': 'Ubiquinone and other terpenoid-quinone biosynthesis',
    'ko00140': 'Steroid hormone biosynthesis',
    'ko00190': 'Oxidative phosphorylation',
    'ko00220': 'Arginine biosynthesis',
    'ko00230': 'Purine metabolism',
    'ko00232': 'Caffeine metabolism',
    'ko00240': 'Pyrimidine metabolism',
    'ko00250': 'Alanine, aspartate and glutamate metabolism',
    'ko00260': 'Glycine, serine and threonine metabolism',
    'ko00270': 'Cysteine and methionine metabolism',
    'ko00280': 'Valine, leucine and isoleucine degradation',
    'ko00290': 'Valine, leucine and isoleucine biosynthesis',
    'ko00300': 'Lysine biosynthesis',
    'ko00310': 'Lysine degradation',
    'ko00330': 'Arginine and proline metabolism',
    'ko00340': 'Histidine metabolism',
    'ko00350': 'Tyrosine metabolism',
    'ko00360': 'Phenylalanine metabolism',
    'ko00380': 'Tryptophan metabolism',
    'ko00400': 'Phenylalanine, tyrosine and tryptophan biosynthesis',
    'ko00410': 'beta-Alanine metabolism',
    'ko00430': 'Taurine and hypotaurine metabolism',
    'ko00440': 'Phosphonate and phosphinate metabolism',
    'ko00450': 'Selenocompound metabolism',
    'ko00470': 'D-Amino acid metabolism',
    'ko00480': 'Glutathione metabolism',
    'ko00500': 'Starch and sucrose metabolism',
    'ko00510': 'N-Glycan biosynthesis',
    'ko00515': 'Mannose type O-glycan biosynthesis',
    'ko00520': 'Amino sugar and nucleotide sugar metabolism',
    'ko00524': 'Neomycin, kanamycin and gentamicin biosynthesis',
    'ko00561': 'Glycerolipid metabolism',
    'ko00562': 'Inositol phosphate metabolism',
    'ko00563': 'Glycosylphosphatidylinositol (GPI)-anchor biosynthesis',
    'ko00564': 'Glycerophospholipid metabolism',
    'ko00565': 'Ether lipid metabolism',
    'ko00590': 'Arachidonic acid metabolism',
    'ko00591': 'Linoleic acid metabolism',
    'ko00592': 'alpha-Linolenic acid metabolism',
    'ko00600': 'Sphingolipid metabolism',
    'ko00620': 'Pyruvate metabolism',
    'ko00630': 'Glyoxylate and dicarboxylate metabolism',
    'ko00640': 'Propanoate metabolism',
    'ko00650': 'Butanoate metabolism',
    'ko00670': 'One carbon pool by folate',
    'ko00710': 'Carbon fixation in photosynthetic organisms',
    'ko00720': 'Carbon fixation pathways in prokaryotes',
    'ko00730': 'Thiamine metabolism',
    'ko00740': 'Riboflavin metabolism',
    'ko00750': 'Vitamin B6 metabolism',
    'ko00760': 'Nicotinate and nicotinamide metabolism',
    'ko00770': 'Pantothenate and CoA biosynthesis',
    'ko00780': 'Biotin metabolism',
    'ko00785': 'Lipoic acid metabolism',
    'ko00790': 'Folate biosynthesis',
    'ko00830': 'Retinol metabolism',
    'ko00860': 'Porphyrin and chlorophyll metabolism',
    'ko00900': 'Terpenoid backbone biosynthesis',
    'ko00910': 'Nitrogen metabolism',
    'ko00920': 'Sulfur metabolism',
    'ko00970': 'Aminoacyl-tRNA biosynthesis',
    'ko00980': 'Metabolism of xenobiotics by cytochrome P450',
    'ko00982': 'Drug metabolism - cytochrome P450',
    'ko00983': 'Drug metabolism - other enzymes',
    'ko01040': 'Biosynthesis of unsaturated fatty acids',
    'ko01100': 'Metabolic pathways',
    'ko01200': 'Carbon metabolism',
    'ko01210': '2-Oxocarboxylic acid metabolism',
    'ko01212': 'Fatty acid metabolism',
    'ko01230': 'Biosynthesis of amino acids',
    'ko01232': 'Nucleotide metabolism',
    'ko01240': 'Biosynthesis of cofactors',
    'ko01250': 'Biosynthesis of nucleotide sugars',
    'ko01521': 'EGFR tyrosine kinase inhibitor resistance',
    'ko01522': 'Endocrine resistance',
    'ko01523': 'Antifolate resistance',
    'ko04141': 'Protein processing in endoplasmic reticulum',
    'ko04142': 'Lysosome',
    'ko04146': 'Peroxisome',
    'ko04920': 'Adipocytokine signaling pathway',
    'ko04910': 'Insulin signaling pathway',
    'ko04260': 'Cardiac muscle contraction',
    'ko04923': 'Regulation of lipolysis in adipocytes',
    'ko04975': 'Fat digestion and absorption',
    'ko04714': 'Thermogenesis',
    'ko04723': 'Retrograde endocannabinoid signaling',
    'ko04931': 'Insulin resistance',
    'ko04977': 'Vitamin digestion and absorption',
    'ko04979': 'Cholesterol metabolism',
    'ko05231': 'Choline metabolism in cancer',
    'ko05417': 'Lipid and atherosclerosis',
}


# KEGG BRITE br08901 中 Metabolism > Global and overview maps 的全部条目
KEGG_OVERVIEW_PATHWAYS = {
    "ko01100",  # Metabolic pathways
    "ko01110",  # Biosynthesis of secondary metabolites
    "ko01120",  # Microbial metabolism in diverse environments
    "ko01200",  # Carbon metabolism
    "ko01210",  # 2-Oxocarboxylic acid metabolism
    "ko01212",  # Fatty acid metabolism
    "ko01220",  # Degradation of aromatic compounds
    "ko01230",  # Biosynthesis of amino acids
    "ko01232",  # Nucleotide metabolism
    "ko01240",  # Biosynthesis of cofactors
    "ko01250",  # Biosynthesis of nucleotide sugars
    "ko01310",  # Nitrogen cycle
    "ko01320",  # Sulfur cycle
}

# ko015xx 编号落在 ko00001–ko01999 内，但属于 Human Diseases > Drug resistance
# （如 ko01521 EGFR TKI resistance、ko01522 Endocrine resistance、ko01523 Antifolate resistance）
KEGG_DRUG_RESISTANCE_RANGE = range(1500, 1600)


def is_pathway_in_scope(
    ko_id: str,
    scope: str = "metabolism",
    *,
    exclude_overview: bool = True,
) -> bool:
    """判断通路是否进入统计集合；默认仅分析 KEGG Metabolism 类。

    ``metabolism``：ko00001–ko01999，去掉 ko015xx 耐药图。``all``：任意 ``koNNNNN``。

    ``exclude_overview`` 与 scope 无关，去掉 KEGG 归为 Global and overview maps 的全部通路，
    包括 ko012xx（Carbon metabolism、Biosynthesis of amino acids 等）。这些图是具体代谢图的并集，
    同一批代谢物会被重复检验；集合大（真实小鼠数据中 20–70 个已注释代谢物），只增加 BH 分母而
    不提供独立证据，显著时也只能回到组成它的具体通路解释。需要时设为 False 保留。
    """
    pathway_id = str(ko_id).strip()
    if exclude_overview and pathway_id in KEGG_OVERVIEW_PATHWAYS:
        return False
    normalized_scope = str(scope).strip().lower()
    if normalized_scope == "all":
        return bool(re.fullmatch(r"ko\d{5}", pathway_id))
    if normalized_scope != "metabolism":
        raise ValueError("pathway_scope 仅支持 'metabolism' 或 'all'")
    match = re.fullmatch(r"ko(\d{5})", pathway_id)
    if not match:
        return False
    number = int(match.group(1))
    return 0 < number < 2000 and number not in KEGG_DRUG_RESISTANCE_RANGE


def get_pathway_name(ko_id: str) -> str:
    """查询 KEGG 通路名称；未知名称显式标记，避免误作已注释结果。"""
    return KEGG_PATHWAYS.get(ko_id, f"{ko_id} (name unavailable)")


def build_metabolite_sets(kegg_map_series: "pd.Series", index_series: "pd.Series") -> dict[str, set[str]]:
    """
    从 kegg_map 列构建 pathway_id -> set(metabolite_index) 映射。

    Parameters
    ----------
    kegg_map_series : pd.Series
        每行为逗号分隔的 KEGG pathway ID（如 "ko00010,ko01100"）
    index_series : pd.Series
        对应的代谢物 Index 列

    Returns
    -------
    dict[str, set[str]]
        {pathway_id: {metabolite_index_1, metabolite_index_2, ...}}
    """
    import pandas as pd
    pathway_map: dict[str, set[str]] = {}
    for idx_val, kegg_str in zip(index_series, kegg_map_series):
        if pd.isna(kegg_str) or str(kegg_str).strip() == "":
            continue
        for ko in str(kegg_str).split(","):
            ko = ko.strip()
            if not ko or ko == "-":
                continue
            if ko not in pathway_map:
                pathway_map[ko] = set()
            pathway_map[ko].add(str(idx_val))
    return pathway_map
