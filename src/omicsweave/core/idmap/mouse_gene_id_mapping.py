#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mouse Gene ID Mapping — NCBI Gene ID / Ensembl / Alias Resolution

Central reference for all omics analysis tools (transcriptomics, proteomics,
metabolomics, multi-omics integration). Maps canonical MGI gene symbols to
standardized database identifiers.

Each entry provides:
    - ncbi_gene_id   : NCBI Gene (Entrez) ID (e.g. '15275')
    - ensembl_id     : Ensembl Gene ID (e.g. 'ENSMUSG00000030278')
    - gene_name      : Full gene name / description
    - gene_type      : protein_coding, ncRNA, etc.
    - chromosome     : Chromosome location
    - aliases        : Known alternative / historical symbols

Usage:
    from omicsweave.core.idmap.mouse_gene_id_mapping import (
        GENE_ANNOTATIONS, GENE_ALIAS_LOOKUP,
        resolve_gene_symbol, get_ncbi_gene_id, get_ensembl_gene_id,
    )
    ncbi_id = get_ncbi_gene_id('Hk1')  # -> '15275'

Gene symbols follow MGI (Mouse Genome Informatics) nomenclature.
All IDs are for Mus musculus (taxid: 10090).

References:
    - MGI:    https://www.informatics.jax.org/
    - NCBI:   https://www.ncbi.nlm.nih.gov/gene/
    - Ensembl: https://www.ensembl.org/Mus_musculus/
"""

from typing import Dict, List, Optional


# ══════════════════════════════════════════════════════════════════════════════
# ALIAS / NON-STANDARD SYMBOL RESOLUTION
# Maps non-standard gene symbols to their canonical MGI symbols.
# None = not a gene (metabolite, cofactor, or other non-gene entry)
# ══════════════════════════════════════════════════════════════════════════════
GENE_ALIAS_LOOKUP: Dict[str, Optional[str]] = {
    # ── Malformed gene symbols (wrong capitalization or numbering) ──
    'Pdh1': 'Pdha1',
    'Pdh1a': 'Pdha1',
    'Ogdh2': 'Ogdh',
    'Suc1': 'Suclg1',
    'Suc2': 'Sucla2',
    'Pgls': 'Pgl',
    'Gpk': 'Gyk',
    'Acats1': 'Acat1',
    'Tkt1': 'Tkt',
    'Rpi1': 'Rpia',
    'Ugpdh': 'Ugp2',
    'Gpdt1': 'Gpam',
    'Gpi1': 'Gpi',
    'Bpgm1': 'Bpgm',

    # ── Protein name aliases → canonical gene symbols ──
    'Atgl': 'Pnpla2',
    'Hsl': 'Lipe',
    'Mct1': 'Slc16a1',
    'Fatp1': 'Slc27a1',
    'Fatp2': 'Slc27a2',
    'Fatp3': 'Slc27a3',
    'Fatp4': 'Slc27a4',
    'Cd147': 'Bsg',
    'Nrf2': 'Nfe2l2',  # common alias
    'Chrebp': 'Mlxipl',  # ChREBP = MLXIPL protein

    # ── Nuclear receptor / TF naming errors ──
    'Ppar1a': 'Ppara',
    'Ppar1b': 'Ppard',
    'Ppar1g': 'Pparg',
    'Ppar2': 'Ppard',
    'Ppar2g': 'Pparg',
    'Pparg1a': 'Ppargc1a',
    'Pparg1b': 'Ppargc1b',
    'Pparg2': 'Pparg',
    'Srebp1': 'Srebf1',
    'Srebp1c': 'Srebf1',
    'Lxr1a': 'Nr1h3',
    'Fx1r': 'Nr1h4',

    # ── AMPK naming variants ──
    'Ampk1a': 'Prkaa1',
    'Ampk1b': 'Prkaa2',

    # ── Non-gene entries (metabolites/cofactors mistakenly in gene lists) ──
    'Oaa': None,        # Oxaloacetate = metabolite
    'Nadh': None,       # NADH = cofactor
    'Nadu': None,       # NAD+ abbreviation = cofactor
    'Biotin': None,     # Vitamin/cofactor
    'Carnitine': None,  # Metabolite
    'Pp1': None,        # Protein phosphatase 1 = protein complex, not gene
    'Ppk1': None,       # Not a valid mouse gene symbol

    # ── Other aliases ──
    'Gcn2': 'Eif2ak4', # GCN2 = EIF2AK4
    'Phd1': 'Egln2',   # PHD1 = EGLN2
    'Phd2': 'Egln1',   # PHD2 = EGLN1
    'Nampt1': 'Nampt', # Nampt1 is an alias for Nampt
    'Purat': None,     # purine salvage enzyme (underspecified), not a valid gene
}


# ══════════════════════════════════════════════════════════════════════════════
# GENE ANNOTATIONS
# Key = canonical MGI gene symbol
# NCBI Gene IDs will be populated from batch NCBI query
# ══════════════════════════════════════════════════════════════════════════════
GENE_ANNOTATIONS: Dict[str, Dict[str, Optional[str]]] = {

    # ══════════════════════════════════════════════════════════════════════════
    # GLYCOLYSIS / GLUCONEOGENESIS
    # ══════════════════════════════════════════════════════════════════════════
    'Hk1': {'ncbi_gene_id': '15275', 'ensembl_id': 'ENSMUSG00000030278', 'gene_name': 'hexokinase 1', 'gene_type': 'protein_coding', 'chromosome': '10', 'aliases': ['Hk-1']},
    'Hk2': {'ncbi_gene_id': '15277', 'ensembl_id': 'ENSMUSG00000020805', 'gene_name': 'hexokinase 2', 'gene_type': 'protein_coding', 'chromosome': '6', 'aliases': []},
    'Hk3': {'ncbi_gene_id': '15204', 'ensembl_id': 'ENSMUSG00000029416', 'gene_name': 'hexokinase 3', 'gene_type': 'protein_coding', 'chromosome': '3', 'aliases': []},
    'Gck': {'ncbi_gene_id': '103988', 'ensembl_id': 'ENSMUSG00000041798', 'gene_name': 'glucokinase', 'gene_type': 'protein_coding', 'chromosome': '11', 'aliases': ['Hk4', 'Hk-4', 'Gk']},
    'Gckr': {'ncbi_gene_id': '145953', 'ensembl_id': 'ENSMUSG00000027048', 'gene_name': 'glucokinase regulator', 'gene_type': 'protein_coding', 'chromosome': '5', 'aliases': ['Gkr']},
    'Hkdc1': {'ncbi_gene_id': '216379', 'ensembl_id': 'ENSMUSG00000096143', 'gene_name': 'hexokinase domain containing 1', 'gene_type': 'protein_coding', 'chromosome': '10', 'aliases': []},
    'Gpi': {'ncbi_gene_id': '14751', 'ensembl_id': 'ENSMUSG00000022143', 'gene_name': 'glucose-6-phosphate isomerase', 'gene_type': 'protein_coding', 'chromosome': '7', 'aliases': ['Gpi1']},
    'Pfkp': {'ncbi_gene_id': '56418', 'ensembl_id': 'ENSMUSG00000040655', 'gene_name': 'phosphofructokinase, platelet', 'gene_type': 'protein_coding', 'chromosome': '13', 'aliases': []},
    'Pfkm': {'ncbi_gene_id': '18642', 'ensembl_id': 'ENSMUSG00000033071', 'gene_name': 'phosphofructokinase, muscle', 'gene_type': 'protein_coding', 'chromosome': '9', 'aliases': []},
    'Pfkl': {'ncbi_gene_id': '56399', 'ensembl_id': 'ENSMUSG00000030027', 'gene_name': 'phosphofructokinase, liver', 'gene_type': 'protein_coding', 'chromosome': '10', 'aliases': []},
    'Pfkfb1': {'ncbi_gene_id': '56320', 'ensembl_id': 'ENSMUSG00000026858', 'gene_name': '6-phosphofructo-2-kinase/fructose-2,6-bisphosphatase 1', 'gene_type': 'protein_coding', 'chromosome': '4', 'aliases': []},
    'Pfkfb2': {'ncbi_gene_id': '56321', 'ensembl_id': 'ENSMUSG00000022472', 'gene_name': '6-phosphofructo-2-kinase/fructose-2,6-bisphosphatase 2', 'gene_type': 'protein_coding', 'chromosome': '1', 'aliases': []},
    'Pfkfb3': {'ncbi_gene_id': '56322', 'ensembl_id': 'ENSMUSG00000026639', 'gene_name': '6-phosphofructo-2-kinase/fructose-2,6-bisphosphatase 3', 'gene_type': 'protein_coding', 'chromosome': '2', 'aliases': []},
    'Pfkfb4': {'ncbi_gene_id': '270198', 'ensembl_id': 'ENSMUSG00000045750', 'gene_name': '6-phosphofructo-2-kinase/fructose-2,6-bisphosphatase 4', 'gene_type': 'protein_coding', 'chromosome': '9', 'aliases': []},
    'Aldoa': {'ncbi_gene_id': '11668', 'ensembl_id': 'ENSMUSG00000030195', 'gene_name': 'aldolase, fructose-bisphosphate A', 'gene_type': 'protein_coding', 'chromosome': '4', 'aliases': []},
    'Aldob': {'ncbi_gene_id': '230767', 'ensembl_id': 'ENSMUSG00000056470', 'gene_name': 'aldolase, fructose-bisphosphate B', 'gene_type': 'protein_coding', 'chromosome': '9', 'aliases': []},
    'Aldoc': {'ncbi_gene_id': '11670', 'ensembl_id': 'ENSMUSG00000030947', 'gene_name': 'aldolase, fructose-bisphosphate C', 'gene_type': 'protein_coding', 'chromosome': '11', 'aliases': []},
    'Aldoart1': {'ncbi_gene_id': '97823', 'ensembl_id': 'ENSMUSG00000076500', 'gene_name': 'aldolase A retrogene 1', 'gene_type': 'protein_coding', 'chromosome': '4', 'aliases': []},
    'Tpi1': {'ncbi_gene_id': '21991', 'ensembl_id': 'ENSMUSG00000023305', 'gene_name': 'triosephosphate isomerase 1', 'gene_type': 'protein_coding', 'chromosome': '1', 'aliases': []},
    'Gapdh': {'ncbi_gene_id': '14433', 'ensembl_id': 'ENSMUSG00000016614', 'gene_name': 'glyceraldehyde-3-phosphate dehydrogenase', 'gene_type': 'protein_coding', 'chromosome': '7', 'aliases': []},
    'Gapdhs': {'ncbi_gene_id': '238627', 'ensembl_id': 'ENSMUSG00000041474', 'gene_name': 'glyceraldehyde-3-phosphate dehydrogenase, spermatogenic', 'gene_type': 'protein_coding', 'chromosome': '7', 'aliases': []},
    'Pgk1': {'ncbi_gene_id': '18654', 'ensembl_id': 'ENSMUSG00000053690', 'gene_name': 'phosphoglycerate kinase 1', 'gene_type': 'protein_coding', 'chromosome': 'X', 'aliases': []},
    'Pgk2': {'ncbi_gene_id': '18655', 'ensembl_id': 'ENSMUSG00000029483', 'gene_name': 'phosphoglycerate kinase 2', 'gene_type': 'protein_coding', 'chromosome': '12', 'aliases': []},
    'Pgam1': {'ncbi_gene_id': '18740', 'ensembl_id': 'ENSMUSG00000031636', 'gene_name': 'phosphoglycerate mutase 1', 'gene_type': 'protein_coding', 'chromosome': '5', 'aliases': []},
    'Pgam2': {'ncbi_gene_id': '18741', 'ensembl_id': 'ENSMUSG00000031635', 'gene_name': 'phosphoglycerate mutase 2', 'gene_type': 'protein_coding', 'chromosome': '5', 'aliases': []},
    'Bpgm': {'ncbi_gene_id': '19883', 'ensembl_id': 'ENSMUSG00000020236', 'gene_name': 'bisphosphoglycerate mutase', 'gene_type': 'protein_coding', 'chromosome': '6', 'aliases': []},
    'Eno1': {'ncbi_gene_id': '13806', 'ensembl_id': 'ENSMUSG00000031707', 'gene_name': 'enolase 1, (alpha)', 'gene_type': 'protein_coding', 'chromosome': '4', 'aliases': []},
    'Eno2': {'ncbi_gene_id': '13807', 'ensembl_id': 'ENSMUSG00000015226', 'gene_name': 'enolase 2, (gamma)', 'gene_type': 'protein_coding', 'chromosome': '6', 'aliases': []},
    'Eno3': {'ncbi_gene_id': '13808', 'ensembl_id': 'ENSMUSG00000026277', 'gene_name': 'enolase 3, (beta)', 'gene_type': 'protein_coding', 'chromosome': '11', 'aliases': []},
    'Pkm': {'ncbi_gene_id': '18746', 'ensembl_id': 'ENSMUSG00000031414', 'gene_name': 'pyruvate kinase M', 'gene_type': 'protein_coding', 'chromosome': '9', 'aliases': ['Pkm2']},
    'Pklr': {'ncbi_gene_id': '18745', 'ensembl_id': 'ENSMUSG00000060316', 'gene_name': 'pyruvate kinase L/R', 'gene_type': 'protein_coding', 'chromosome': '3', 'aliases': []},
    'Ldha': {'ncbi_gene_id': '16828', 'ensembl_id': 'ENSMUSG00000029518', 'gene_name': 'lactate dehydrogenase A', 'gene_type': 'protein_coding', 'chromosome': '7', 'aliases': []},
    'Ldhb': {'ncbi_gene_id': '16832', 'ensembl_id': 'ENSMUSG00000030458', 'gene_name': 'lactate dehydrogenase B', 'gene_type': 'protein_coding', 'chromosome': '6', 'aliases': []},
    'Ldhc': {'ncbi_gene_id': '16833', 'ensembl_id': 'ENSMUSG00000036218', 'gene_name': 'lactate dehydrogenase C', 'gene_type': 'protein_coding', 'chromosome': '12', 'aliases': []},
    'Fbp1': {'ncbi_gene_id': '14121', 'ensembl_id': 'ENSMUSG00000029153', 'gene_name': 'fructose-bisphosphatase 1', 'gene_type': 'protein_coding', 'chromosome': '13', 'aliases': []},
    'Fbp2': {'ncbi_gene_id': '14122', 'ensembl_id': 'ENSMUSG00000024159', 'gene_name': 'fructose-bisphosphatase 2', 'gene_type': 'protein_coding', 'chromosome': '9', 'aliases': []},
    'Pck1': {'ncbi_gene_id': '18539', 'ensembl_id': 'ENSMUSG00000027512', 'gene_name': 'phosphoenolpyruvate carboxykinase 1, cytosolic', 'gene_type': 'protein_coding', 'chromosome': '2', 'aliases': ['Pepck']},
    'Pck2': {'ncbi_gene_id': '72300', 'ensembl_id': 'ENSMUSG00000023698', 'gene_name': 'phosphoenolpyruvate carboxykinase 2, mitochondrial', 'gene_type': 'protein_coding', 'chromosome': '14', 'aliases': []},
    'Pc': {'ncbi_gene_id': '18563', 'ensembl_id': 'ENSMUSG00000024406', 'gene_name': 'pyruvate carboxylase', 'gene_type': 'protein_coding', 'chromosome': '19', 'aliases': []},
    'G6pc': {'ncbi_gene_id': '14377', 'ensembl_id': 'ENSMUSG00000024471', 'gene_name': 'glucose-6-phosphatase, catalytic subunit', 'gene_type': 'protein_coding', 'chromosome': '11', 'aliases': []},
    'G6pc2': {'ncbi_gene_id': '26375', 'ensembl_id': 'ENSMUSG00000022511', 'gene_name': 'glucose-6-phosphatase 2', 'gene_type': 'protein_coding', 'chromosome': '2', 'aliases': []},
    'G6pc3': {'ncbi_gene_id': '26374', 'ensembl_id': 'ENSMUSG00000020116', 'gene_name': 'glucose-6-phosphatase 3', 'gene_type': 'protein_coding', 'chromosome': '10', 'aliases': []},
    'Pgm1': {'ncbi_gene_id': '18805', 'ensembl_id': 'ENSMUSG00000001440', 'gene_name': 'phosphoglucomutase 1', 'gene_type': 'protein_coding', 'chromosome': '5', 'aliases': []},
    'Pgm2': {'ncbi_gene_id': '218693', 'ensembl_id': 'ENSMUSG00000026363', 'gene_name': 'phosphoglucomutase 2', 'gene_type': 'protein_coding', 'chromosome': '13', 'aliases': []},
    'Pgm3': {'ncbi_gene_id': '103968', 'ensembl_id': 'ENSMUSG00000030948', 'gene_name': 'phosphoglucomutase 3', 'gene_type': 'protein_coding', 'chromosome': '11', 'aliases': []},
    'Pgm5': {'ncbi_gene_id': '228349', 'ensembl_id': 'ENSMUSG00000032395', 'gene_name': 'phosphoglucomutase 5', 'gene_type': 'protein_coding', 'chromosome': '5', 'aliases': []},
    'Gys1': {'ncbi_gene_id': '14936', 'ensembl_id': 'ENSMUSG00000030156', 'gene_name': 'glycogen synthase 1, muscle', 'gene_type': 'protein_coding', 'chromosome': '15', 'aliases': []},
    'Gys2': {'ncbi_gene_id': '14937', 'ensembl_id': 'ENSMUSG00000061046', 'gene_name': 'glycogen synthase 2, liver', 'gene_type': 'protein_coding', 'chromosome': '6', 'aliases': []},
    'Pygb': {'ncbi_gene_id': '19077', 'ensembl_id': 'ENSMUSG00000026597', 'gene_name': 'glycogen phosphorylase, brain', 'gene_type': 'protein_coding', 'chromosome': '15', 'aliases': []},
    'Pygm': {'ncbi_gene_id': '19079', 'ensembl_id': 'ENSMUSG00000030842', 'gene_name': 'glycogen phosphorylase, muscle associated', 'gene_type': 'protein_coding', 'chromosome': '19', 'aliases': []},
    'Pygl': {'ncbi_gene_id': '19078', 'ensembl_id': 'ENSMUSG00000029647', 'gene_name': 'glycogen phosphorylase, liver', 'gene_type': 'protein_coding', 'chromosome': '12', 'aliases': []},
    'Ugp2': {'ncbi_gene_id': '22337', 'ensembl_id': 'ENSMUSG00000024443', 'gene_name': 'UDP-glucose pyrophosphorylase 2', 'gene_type': 'protein_coding', 'chromosome': '4', 'aliases': []},
    'Gsk3a': {'ncbi_gene_id': '60449', 'ensembl_id': 'ENSMUSG00000030834', 'gene_name': 'glycogen synthase kinase 3 alpha', 'gene_type': 'protein_coding', 'chromosome': '7', 'aliases': []},
    'Gsk3b': {'ncbi_gene_id': '56637', 'ensembl_id': 'ENSMUSG00000023234', 'gene_name': 'glycogen synthase kinase 3 beta', 'gene_type': 'protein_coding', 'chromosome': '16', 'aliases': []},
    'Ppp1r3a': {'ncbi_gene_id': '19046', 'ensembl_id': 'ENSMUSG00000026002', 'gene_name': 'protein phosphatase 1, regulatory subunit 3A', 'gene_type': 'protein_coding', 'chromosome': '6', 'aliases': []},
    'Ppp1r3b': {'ncbi_gene_id': '69505', 'ensembl_id': 'ENSMUSG00000029998', 'gene_name': 'protein phosphatase 1, regulatory subunit 3B', 'gene_type': 'protein_coding', 'chromosome': '9', 'aliases': []},
    'Ppp1r3c': {'ncbi_gene_id': '71748', 'ensembl_id': 'ENSMUSG00000026508', 'gene_name': 'protein phosphatase 1, regulatory subunit 3C', 'gene_type': 'protein_coding', 'chromosome': '5', 'aliases': []},

    # ══════════════════════════════════════════════════════════════════════════
    # TCA CYCLE / PDH
    # ══════════════════════════════════════════════════════════════════════════
    'Cs': {'ncbi_gene_id': '12995', 'ensembl_id': 'ENSMUSG00000005366', 'gene_name': 'citrate synthase', 'gene_type': 'protein_coding', 'chromosome': '10', 'aliases': []},
    'Aco1': {'ncbi_gene_id': '11428', 'ensembl_id': 'ENSMUSG00000020266', 'gene_name': 'aconitase 1', 'gene_type': 'protein_coding', 'chromosome': '4', 'aliases': ['Irp1']},
    'Aco2': {'ncbi_gene_id': '11429', 'ensembl_id': 'ENSMUSG00000022677', 'gene_name': 'aconitase 2, mitochondrial', 'gene_type': 'protein_coding', 'chromosome': '15', 'aliases': []},
    'Idh3a': {'ncbi_gene_id': '15906', 'ensembl_id': 'ENSMUSG00000005330', 'gene_name': 'isocitrate dehydrogenase 3 (NAD+) alpha', 'gene_type': 'protein_coding', 'chromosome': '1', 'aliases': []},
    'Idh3b': {'ncbi_gene_id': '52519', 'ensembl_id': 'ENSMUSG00000032140', 'gene_name': 'isocitrate dehydrogenase 3 (NAD+) beta', 'gene_type': 'protein_coding', 'chromosome': '3', 'aliases': []},
    'Idh3g': {'ncbi_gene_id': '15907', 'ensembl_id': 'ENSMUSG00000029584', 'gene_name': 'isocitrate dehydrogenase 3 (NAD+) gamma', 'gene_type': 'protein_coding', 'chromosome': '1', 'aliases': []},
    'Idh1': {'ncbi_gene_id': '15926', 'ensembl_id': 'ENSMUSG00000025899', 'gene_name': 'isocitrate dehydrogenase 1 (NADP+), soluble', 'gene_type': 'protein_coding', 'chromosome': '1', 'aliases': []},
    'Idh2': {'ncbi_gene_id': '269951', 'ensembl_id': 'ENSMUSG00000055970', 'gene_name': 'isocitrate dehydrogenase 2 (NADP+), mitochondrial', 'gene_type': 'protein_coding', 'chromosome': '7', 'aliases': []},
    'Ogdh': {'ncbi_gene_id': '18293', 'ensembl_id': 'ENSMUSG00000020545', 'gene_name': 'oxoglutarate dehydrogenase', 'gene_type': 'protein_coding', 'chromosome': '14', 'aliases': []},
    'Dlst': {'ncbi_gene_id': '71744', 'ensembl_id': 'ENSMUSG00000005318', 'gene_name': 'dihydrolipoamide S-succinyltransferase', 'gene_type': 'protein_coding', 'chromosome': '12', 'aliases': []},
    'Dld': {'ncbi_gene_id': '13382', 'ensembl_id': 'ENSMUSG00000028381', 'gene_name': 'dihydrolipoamide dehydrogenase', 'gene_type': 'protein_coding', 'chromosome': '10', 'aliases': []},
    'Sdha': {'ncbi_gene_id': '66945', 'ensembl_id': 'ENSMUSG00000023860', 'gene_name': 'succinate dehydrogenase complex, subunit A', 'gene_type': 'protein_coding', 'chromosome': '12', 'aliases': []},
    'Sdhb': {'ncbi_gene_id': '67680', 'ensembl_id': 'ENSMUSG00000009548', 'gene_name': 'succinate dehydrogenase complex, subunit B', 'gene_type': 'protein_coding', 'chromosome': '4', 'aliases': []},
    'Sdhc': {'ncbi_gene_id': '66052', 'ensembl_id': 'ENSMUSG00000005320', 'gene_name': 'succinate dehydrogenase complex, subunit C', 'gene_type': 'protein_coding', 'chromosome': '1', 'aliases': []},
    'Sdhd': {'ncbi_gene_id': '66946', 'ensembl_id': 'ENSMUSG00000028331', 'gene_name': 'succinate dehydrogenase complex, subunit D', 'gene_type': 'protein_coding', 'chromosome': '11', 'aliases': []},
    'Fh1': {'ncbi_gene_id': '14194', 'ensembl_id': 'ENSMUSG00000030087', 'gene_name': 'fumarate hydratase 1', 'gene_type': 'protein_coding', 'chromosome': '15', 'aliases': []},
    'Mdh1': {'ncbi_gene_id': '17449', 'ensembl_id': 'ENSMUSG00000020343', 'gene_name': 'malate dehydrogenase 1, NAD (soluble)', 'gene_type': 'protein_coding', 'chromosome': '11', 'aliases': []},
    'Mdh2': {'ncbi_gene_id': '17448', 'ensembl_id': 'ENSMUSG00000014808', 'gene_name': 'malate dehydrogenase 2, NAD (mitochondrial)', 'gene_type': 'protein_coding', 'chromosome': '5', 'aliases': []},
    'Suclg1': {'ncbi_gene_id': '20909', 'ensembl_id': 'ENSMUSG00000020501', 'gene_name': 'succinyl-CoA ligase, GDP-forming, alpha subunit', 'gene_type': 'protein_coding', 'chromosome': '6', 'aliases': []},
    'Sucla2': {'ncbi_gene_id': '56480', 'ensembl_id': 'ENSMUSG00000037201', 'gene_name': 'succinyl-CoA ligase, ADP-forming, beta subunit', 'gene_type': 'protein_coding', 'chromosome': '13', 'aliases': []},
    'Pdha1': {'ncbi_gene_id': '18597', 'ensembl_id': 'ENSMUSG00000025867', 'gene_name': 'pyruvate dehydrogenase E1 alpha 1', 'gene_type': 'protein_coding', 'chromosome': '7', 'aliases': []},
    'Pdhb': {'ncbi_gene_id': '18598', 'ensembl_id': 'ENSMUSG00000031802', 'gene_name': 'pyruvate dehydrogenase E1 beta', 'gene_type': 'protein_coding', 'chromosome': '15', 'aliases': []},
    'Dlat': {'ncbi_gene_id': '108667', 'ensembl_id': 'ENSMUSG00000022561', 'gene_name': 'dihydrolipoamide S-acetyltransferase', 'gene_type': 'protein_coding', 'chromosome': '9', 'aliases': []},
    'Pdk1': {'ncbi_gene_id': '228026', 'ensembl_id': 'ENSMUSG00000027122', 'gene_name': 'pyruvate dehydrogenase kinase, isozyme 1', 'gene_type': 'protein_coding', 'chromosome': '2', 'aliases': []},
    'Pdk2': {'ncbi_gene_id': '241398', 'ensembl_id': 'ENSMUSG00000031150', 'gene_name': 'pyruvate dehydrogenase kinase, isozyme 2', 'gene_type': 'protein_coding', 'chromosome': '11', 'aliases': []},
    'Pdk3': {'ncbi_gene_id': '241399', 'ensembl_id': 'ENSMUSG00000032884', 'gene_name': 'pyruvate dehydrogenase kinase, isozyme 3', 'gene_type': 'protein_coding', 'chromosome': '12', 'aliases': []},
    'Pdk4': {'ncbi_gene_id': '27273', 'ensembl_id': 'ENSMUSG00000029517', 'gene_name': 'pyruvate dehydrogenase kinase, isozyme 4', 'gene_type': 'protein_coding', 'chromosome': '1', 'aliases': []},
    'Pdp1': {'ncbi_gene_id': '76831', 'ensembl_id': 'ENSMUSG00000031874', 'gene_name': 'pyruvate dehydrogenase phosphatase catalytic subunit 1', 'gene_type': 'protein_coding', 'chromosome': '1', 'aliases': []},
    'Pdp2': {'ncbi_gene_id': '320629', 'ensembl_id': 'ENSMUSG00000040266', 'gene_name': 'pyruvate dehydrogenase phosphatase catalytic subunit 2', 'gene_type': 'protein_coding', 'chromosome': '16', 'aliases': []},
    'Mpc1': {'ncbi_gene_id': '67774', 'ensembl_id': 'ENSMUSG00000028862', 'gene_name': 'mitochondrial pyruvate carrier 1', 'gene_type': 'protein_coding', 'chromosome': '5', 'aliases': []},
    'Mpc2': {'ncbi_gene_id': '66589', 'ensembl_id': 'ENSMUSG00000028436', 'gene_name': 'mitochondrial pyruvate carrier 2', 'gene_type': 'protein_coding', 'chromosome': '10', 'aliases': []},

    # ── Placeholder: additional genes to be completed from NCBI batch query ──
    # The following sections will be populated after the NCBI batch query completes.
    # Format is maintained for easy merge.

    # ══════════════════════════════════════════════════════════════════════════
    # PPP / OXIDATIVE STRESS
    # ══════════════════════════════════════════════════════════════════════════
    'G6pd': {'ncbi_gene_id': '14381', 'ensembl_id': 'ENSMUSG00000031400', 'gene_name': 'glucose-6-phosphate dehydrogenase', 'gene_type': 'protein_coding', 'chromosome': 'X', 'aliases': []},
    'Pgd': {'ncbi_gene_id': '18261', 'ensembl_id': 'ENSMUSG00000027601', 'gene_name': 'phosphogluconate dehydrogenase', 'gene_type': 'protein_coding', 'chromosome': '4', 'aliases': []},
    'Pgl': {'ncbi_gene_id': '66097', 'ensembl_id': 'ENSMUSG00000025836', 'gene_name': '6-phosphogluconolactonase', 'gene_type': 'protein_coding', 'chromosome': '18', 'aliases': ['Pgls']},
    'Tkt': {'ncbi_gene_id': '21843', 'ensembl_id': 'ENSMUSG00000028335', 'gene_name': 'transketolase', 'gene_type': 'protein_coding', 'chromosome': '14', 'aliases': []},
    'Taldo1': {'ncbi_gene_id': '21346', 'ensembl_id': 'ENSMUSG00000026149', 'gene_name': 'transaldolase 1', 'gene_type': 'protein_coding', 'chromosome': '3', 'aliases': []},
    'Taldo2': {'ncbi_gene_id': '241250', 'ensembl_id': 'ENSMUSG00000042848', 'gene_name': 'transaldolase 2', 'gene_type': 'protein_coding', 'chromosome': '13', 'aliases': []},
    'Rpia': {'ncbi_gene_id': '103833', 'ensembl_id': 'ENSMUSG00000024507', 'gene_name': 'ribose 5-phosphate isomerase A', 'gene_type': 'protein_coding', 'chromosome': '15', 'aliases': ['Rpi1']},
    'Rpe': {'ncbi_gene_id': '20020', 'ensembl_id': 'ENSMUSG00000032217', 'gene_name': 'ribulose-5-phosphate-3-epimerase', 'gene_type': 'protein_coding', 'chromosome': '7', 'aliases': []},
    'Prps1': {'ncbi_gene_id': '19149', 'ensembl_id': 'ENSMUSG00000035570', 'gene_name': 'phosphoribosyl pyrophosphate synthetase 1', 'gene_type': 'protein_coding', 'chromosome': '4', 'aliases': []},
    'Prps2': {'ncbi_gene_id': '19150', 'ensembl_id': 'ENSMUSG00000037934', 'gene_name': 'phosphoribosyl pyrophosphate synthetase 2', 'gene_type': 'protein_coding', 'chromosome': 'X', 'aliases': []},
    'Prpsap1': {'ncbi_gene_id': '108014', 'ensembl_id': 'ENSMUSG00000031425', 'gene_name': 'phosphoribosyl pyrophosphate synthetase-associated protein 1', 'gene_type': 'protein_coding', 'chromosome': '5', 'aliases': []},
    'Sod1': {'ncbi_gene_id': '20655', 'ensembl_id': 'ENSMUSG00000022960', 'gene_name': 'superoxide dismutase 1, soluble', 'gene_type': 'protein_coding', 'chromosome': '16', 'aliases': []},
    'Sod2': {'ncbi_gene_id': '20656', 'ensembl_id': 'ENSMUSG00000023364', 'gene_name': 'superoxide dismutase 2, mitochondrial', 'gene_type': 'protein_coding', 'chromosome': '17', 'aliases': []},
    'Cat': {'ncbi_gene_id': '12359', 'ensembl_id': 'ENSMUSG00000027136', 'gene_name': 'catalase', 'gene_type': 'protein_coding', 'chromosome': '2', 'aliases': []},
    'Gpx1': {'ncbi_gene_id': '14775', 'ensembl_id': 'ENSMUSG00000043970', 'gene_name': 'glutathione peroxidase 1', 'gene_type': 'protein_coding', 'chromosome': '9', 'aliases': []},
    'Gclc': {'ncbi_gene_id': '14629', 'ensembl_id': 'ENSMUSG00000033298', 'gene_name': 'glutamate-cysteine ligase, catalytic subunit', 'gene_type': 'protein_coding', 'chromosome': '9', 'aliases': []},
    'Gclm': {'ncbi_gene_id': '14630', 'ensembl_id': 'ENSMUSG00000023855', 'gene_name': 'glutamate-cysteine ligase, modifier subunit', 'gene_type': 'protein_coding', 'chromosome': '4', 'aliases': []},
    'Gsr': {'ncbi_gene_id': '14782', 'ensembl_id': 'ENSMUSG00000031547', 'gene_name': 'glutathione reductase', 'gene_type': 'protein_coding', 'chromosome': '8', 'aliases': []},
    'Prdx1': {'ncbi_gene_id': '18477', 'ensembl_id': 'ENSMUSG00000028679', 'gene_name': 'peroxiredoxin 1', 'gene_type': 'protein_coding', 'chromosome': '4', 'aliases': []},
    'Prdx3': {'ncbi_gene_id': '18479', 'ensembl_id': 'ENSMUSG00000005370', 'gene_name': 'peroxiredoxin 3', 'gene_type': 'protein_coding', 'chromosome': '1', 'aliases': []},
    'Nfe2l2': {'ncbi_gene_id': '18024', 'ensembl_id': 'ENSMUSG00000020844', 'gene_name': 'NFE2 like bZIP transcription factor 2', 'gene_type': 'protein_coding', 'chromosome': '7', 'aliases': ['Nrf2']},
    'Keap1': {'ncbi_gene_id': '27066', 'ensembl_id': 'ENSMUSG00000032214', 'gene_name': 'Kelch-like ECH-associated protein 1', 'gene_type': 'protein_coding', 'chromosome': '9', 'aliases': []},
    'Nfs1': {'ncbi_gene_id': '18041', 'ensembl_id': 'ENSMUSG00000027618', 'gene_name': 'nitrogen fixation gene 1 (S. cerevisiae)', 'gene_type': 'protein_coding', 'chromosome': '2', 'aliases': []},
    'Iscu': {'ncbi_gene_id': '70694', 'ensembl_id': 'ENSMUSG00000024431', 'gene_name': 'iron-sulfur cluster assembly enzyme', 'gene_type': 'protein_coding', 'chromosome': '5', 'aliases': []},
    'Fxn': {'ncbi_gene_id': '14297', 'ensembl_id': 'ENSMUSG00000031705', 'gene_name': 'frataxin', 'gene_type': 'protein_coding', 'chromosome': '19', 'aliases': []},

    'Acaca': {'ncbi_gene_id': '107476', 'ensembl_id': 'ENSMUSG00000020532', 'gene_name': 'acetyl-CoA carboxylase alpha', 'gene_type': 'protein_coding', 'chromosome': '11', 'aliases': []},
    'Acacb': {'ncbi_gene_id': '100705', 'ensembl_id': 'ENSMUSG00000042010', 'gene_name': 'acetyl-CoA carboxylase beta', 'gene_type': 'protein_coding', 'chromosome': '5', 'aliases': []},
    'Fasn': {'ncbi_gene_id': '14104', 'ensembl_id': 'ENSMUSG00000025153', 'gene_name': 'fatty acid synthase', 'gene_type': 'protein_coding', 'chromosome': '11', 'aliases': []},
    'Scd1': {'ncbi_gene_id': '20249', 'ensembl_id': 'ENSMUSG00000037071', 'gene_name': 'stearoyl-CoA desaturase 1', 'gene_type': 'protein_coding', 'chromosome': '19', 'aliases': []},
    'Scd2': {'ncbi_gene_id': '20250', 'ensembl_id': 'ENSMUSG00000025202', 'gene_name': 'stearoyl-CoA desaturase 2', 'gene_type': 'protein_coding', 'chromosome': '19', 'aliases': []},
    'Fads1': {'ncbi_gene_id': '24390', 'ensembl_id': 'ENSMUSG00000046333', 'gene_name': 'fatty acid desaturase 1', 'gene_type': 'protein_coding', 'chromosome': '19', 'aliases': []},
    'Fads2': {'ncbi_gene_id': '56473', 'ensembl_id': 'ENSMUSG00000046334', 'gene_name': 'fatty acid desaturase 2', 'gene_type': 'protein_coding', 'chromosome': '19', 'aliases': []},
    'Fads3': {'ncbi_gene_id': '74563', 'ensembl_id': 'ENSMUSG00000046335', 'gene_name': 'fatty acid desaturase 3', 'gene_type': 'protein_coding', 'chromosome': '19', 'aliases': []},
    'Elovl1': {'ncbi_gene_id': '54325', 'ensembl_id': 'ENSMUSG00000006390', 'gene_name': 'ELOVL fatty acid elongase 1', 'gene_type': 'protein_coding', 'chromosome': '4', 'aliases': []},
    'Elovl2': {'ncbi_gene_id': '72090', 'ensembl_id': 'ENSMUSG00000031787', 'gene_name': 'ELOVL fatty acid elongase 2', 'gene_type': 'protein_coding', 'chromosome': '13', 'aliases': []},
    'Elovl3': {'ncbi_gene_id': '12686', 'ensembl_id': 'ENSMUSG00000031754', 'gene_name': 'ELOVL fatty acid elongase 3', 'gene_type': 'protein_coding', 'chromosome': '19', 'aliases': []},
    'Elovl4': {'ncbi_gene_id': '109565', 'ensembl_id': 'ENSMUSG00000045568', 'gene_name': 'ELOVL fatty acid elongase 4', 'gene_type': 'protein_coding', 'chromosome': '9', 'aliases': []},
    'Elovl5': {'ncbi_gene_id': '76030', 'ensembl_id': 'ENSMUSG00000032115', 'gene_name': 'ELOVL fatty acid elongase 5', 'gene_type': 'protein_coding', 'chromosome': '9', 'aliases': []},
    'Elovl6': {'ncbi_gene_id': '230235', 'ensembl_id': 'ENSMUSG00000041220', 'gene_name': 'ELOVL fatty acid elongase 6', 'gene_type': 'protein_coding', 'chromosome': '3', 'aliases': []},
    'Elovl7': {'ncbi_gene_id': '53622', 'ensembl_id': 'ENSMUSG00000027867', 'gene_name': 'ELOVL fatty acid elongase 7', 'gene_type': 'protein_coding', 'chromosome': '13', 'aliases': []},
    'Acly': {'ncbi_gene_id': '20943', 'ensembl_id': 'ENSMUSG00000026056', 'gene_name': 'ATP citrate lyase', 'gene_type': 'protein_coding', 'chromosome': '11', 'aliases': []},
    'Acss1': {'ncbi_gene_id': '68738', 'ensembl_id': 'ENSMUSG00000031683', 'gene_name': 'acyl-CoA synthetase short-chain family member 1', 'gene_type': 'protein_coding', 'chromosome': '2', 'aliases': []},
    'Acss2': {'ncbi_gene_id': '20920', 'ensembl_id': 'ENSMUSG00000027605', 'gene_name': 'acyl-CoA synthetase short-chain family member 2', 'gene_type': 'protein_coding', 'chromosome': '10', 'aliases': []},
    'Cpt1a': {'ncbi_gene_id': '12894', 'ensembl_id': 'ENSMUSG00000024900', 'gene_name': 'carnitine palmitoyltransferase 1a', 'gene_type': 'protein_coding', 'chromosome': '19', 'aliases': []},
    'Cpt1b': {'ncbi_gene_id': '12895', 'ensembl_id': 'ENSMUSG00000060093', 'gene_name': 'carnitine palmitoyltransferase 1b', 'gene_type': 'protein_coding', 'chromosome': '15', 'aliases': []},
    'Cpt1c': {'ncbi_gene_id': '78070', 'ensembl_id': 'ENSMUSG00000043261', 'gene_name': 'carnitine palmitoyltransferase 1c', 'gene_type': 'protein_coding', 'chromosome': '7', 'aliases': []},
    'Cpt2': {'ncbi_gene_id': '12922', 'ensembl_id': 'ENSMUSG00000022447', 'gene_name': 'carnitine palmitoyltransferase 2', 'gene_type': 'protein_coding', 'chromosome': '5', 'aliases': []},
    'Acadm': {'ncbi_gene_id': '11398', 'ensembl_id': 'ENSMUSG00000022014', 'gene_name': 'acyl-CoA dehydrogenase medium chain', 'gene_type': 'protein_coding', 'chromosome': '3', 'aliases': []},
    'Acadl': {'ncbi_gene_id': '11399', 'ensembl_id': 'ENSMUSG00000026030', 'gene_name': 'acyl-CoA dehydrogenase long chain', 'gene_type': 'protein_coding', 'chromosome': '13', 'aliases': []},
    'Acads': {'ncbi_gene_id': '11400', 'ensembl_id': 'ENSMUSG00000029379', 'gene_name': 'acyl-CoA dehydrogenase short chain', 'gene_type': 'protein_coding', 'chromosome': '12', 'aliases': []},
    'Acadvl': {'ncbi_gene_id': '11360', 'ensembl_id': 'ENSMUSG00000038568', 'gene_name': 'acyl-CoA dehydrogenase very long chain', 'gene_type': 'protein_coding', 'chromosome': '12', 'aliases': []},
    'Acad9': {'ncbi_gene_id': '26039', 'ensembl_id': 'ENSMUSG00000027162', 'gene_name': 'acyl-CoA dehydrogenase family member 9', 'gene_type': 'protein_coding', 'chromosome': '19', 'aliases': []},
    'Acad10': {'ncbi_gene_id': '68147', 'ensembl_id': 'ENSMUSG00000042644', 'gene_name': 'acyl-CoA dehydrogenase family member 10', 'gene_type': 'protein_coding', 'chromosome': '5', 'aliases': []},
    'Hadh': {'ncbi_gene_id': '15107', 'ensembl_id': 'ENSMUSG00000027984', 'gene_name': 'hydroxyacyl-CoA dehydrogenase', 'gene_type': 'protein_coding', 'chromosome': '12', 'aliases': []},
    'Hadha': {'ncbi_gene_id': '17058', 'ensembl_id': 'ENSMUSG00000035238', 'gene_name': 'hydroxyacyl-CoA dehydrogenase trifunctional multienzyme complex subunit alpha', 'gene_type': 'protein_coding', 'chromosome': '8', 'aliases': []},
    'Hadhb': {'ncbi_gene_id': '16961', 'ensembl_id': 'ENSMUSG00000035237', 'gene_name': 'hydroxyacyl-CoA dehydrogenase trifunctional multienzyme complex subunit beta', 'gene_type': 'protein_coding', 'chromosome': '8', 'aliases': []},
    'Ehhadh': {'ncbi_gene_id': '74147', 'ensembl_id': 'ENSMUSG00000032255', 'gene_name': 'enoyl-CoA hydratase and 3-hydroxyacyl CoA dehydrogenase', 'gene_type': 'protein_coding', 'chromosome': '16', 'aliases': []},
    'Mlycd': {'ncbi_gene_id': '56690', 'ensembl_id': 'ENSMUSG00000074064', 'gene_name': 'malonyl-CoA decarboxylase', 'gene_type': 'protein_coding', 'chromosome': '8', 'aliases': []},
    'Bdh1': {'ncbi_gene_id': '71911', 'ensembl_id': 'ENSMUSG00000030201', 'gene_name': '3-hydroxybutyrate dehydrogenase 1', 'gene_type': 'protein_coding', 'chromosome': '16', 'aliases': []},
    'Hmgcs2': {'ncbi_gene_id': '50762', 'ensembl_id': 'ENSMUSG00000029287', 'gene_name': '3-hydroxy-3-methylglutaryl-CoA synthase 2', 'gene_type': 'protein_coding', 'chromosome': '3', 'aliases': []},
    'Hmgcl': {'ncbi_gene_id': '17897', 'ensembl_id': 'ENSMUSG00000022183', 'gene_name': '3-hydroxy-3-methylglutaryl-CoA lyase', 'gene_type': 'protein_coding', 'chromosome': '4', 'aliases': []},
    'Oxct1': {'ncbi_gene_id': '67041', 'ensembl_id': 'ENSMUSG00000022186', 'gene_name': '3-oxoacid CoA transferase 1', 'gene_type': 'protein_coding', 'chromosome': '15', 'aliases': []},
    'Oxct2': {'ncbi_gene_id': '17071', 'ensembl_id': 'ENSMUSG00000031711', 'gene_name': '3-oxoacid CoA transferase 2', 'gene_type': 'protein_coding', 'chromosome': '9', 'aliases': []},
    'Acat1': {'ncbi_gene_id': '110446', 'ensembl_id': 'ENSMUSG00000032047', 'gene_name': 'acetyl-CoA acetyltransferase 1', 'gene_type': 'protein_coding', 'chromosome': '11', 'aliases': []},
    'Glud1': {'ncbi_gene_id': '14410', 'ensembl_id': 'ENSMUSG00000021794', 'gene_name': 'glutamate dehydrogenase 1', 'gene_type': 'protein_coding', 'chromosome': '14', 'aliases': []},
    'Glud2': {'ncbi_gene_id': '9778', 'ensembl_id': 'ENSMUSG00000025824', 'gene_name': 'glutamate dehydrogenase 2', 'gene_type': 'protein_coding', 'chromosome': '8', 'aliases': []},
    'Gls': {'ncbi_gene_id': '14627', 'ensembl_id': 'ENSMUSG00000026103', 'gene_name': 'glutaminase', 'gene_type': 'protein_coding', 'chromosome': '1', 'aliases': []},
    'Gls2': {'ncbi_gene_id': '65080', 'ensembl_id': 'ENSMUSG00000028513', 'gene_name': 'glutaminase 2 (liver)', 'gene_type': 'protein_coding', 'chromosome': '10', 'aliases': []},
    'Glul': {'ncbi_gene_id': '14645', 'ensembl_id': 'ENSMUSG00000026473', 'gene_name': 'glutamate-ammonia ligase (glutamine synthetase)', 'gene_type': 'protein_coding', 'chromosome': '1', 'aliases': []},
    'Got1': {'ncbi_gene_id': '14718', 'ensembl_id': 'ENSMUSG00000031668', 'gene_name': 'glutamic-oxaloacetic transaminase 1', 'gene_type': 'protein_coding', 'chromosome': '15', 'aliases': []},
    'Got2': {'ncbi_gene_id': '14719', 'ensembl_id': 'ENSMUSG00000023459', 'gene_name': 'glutamic-oxaloacetic transaminase 2', 'gene_type': 'protein_coding', 'chromosome': '8', 'aliases': []},
    'Gpt': {'ncbi_gene_id': '14851', 'ensembl_id': 'ENSMUSG00000022546', 'gene_name': 'glutamic pyruvic transaminase', 'gene_type': 'protein_coding', 'chromosome': '15', 'aliases': []},
    'Gpt2': {'ncbi_gene_id': '20593', 'ensembl_id': 'ENSMUSG00000038238', 'gene_name': 'glutamic pyruvate transaminase 2', 'gene_type': 'protein_coding', 'chromosome': '8', 'aliases': []},
    'Me1': {'ncbi_gene_id': '14422', 'ensembl_id': 'ENSMUSG00000050240', 'gene_name': 'malic enzyme 1 (NADP+ dependent, cytosolic)', 'gene_type': 'protein_coding', 'chromosome': '9', 'aliases': []},
    'Me2': {'ncbi_gene_id': '14184', 'ensembl_id': 'ENSMUSG00000026993', 'gene_name': 'malic enzyme 2 (NAD+ dependent, mitochondrial)', 'gene_type': 'protein_coding', 'chromosome': '18', 'aliases': []},
    'Me3': {'ncbi_gene_id': '228960', 'ensembl_id': 'ENSMUSG00000039068', 'gene_name': 'malic enzyme 3 (NADP+ dependent, mitochondrial)', 'gene_type': 'protein_coding', 'chromosome': '17', 'aliases': []},
    'Cps1': {'ncbi_gene_id': '19185', 'ensembl_id': 'ENSMUSG00000025271', 'gene_name': 'carbamoyl-phosphate synthetase 1', 'gene_type': 'protein_coding', 'chromosome': '4', 'aliases': []},
    'Otc': {'ncbi_gene_id': '18418', 'ensembl_id': 'ENSMUSG00000031144', 'gene_name': 'ornithine transcarbamylase', 'gene_type': 'protein_coding', 'chromosome': '1', 'aliases': []},
    'Ass1': {'ncbi_gene_id': '11841', 'ensembl_id': 'ENSMUSG00000025289', 'gene_name': 'argininosuccinate synthetase 1', 'gene_type': 'protein_coding', 'chromosome': '11', 'aliases': []},
    'Phgdh': {'ncbi_gene_id': '236539', 'ensembl_id': 'ENSMUSG00000032767', 'gene_name': 'phosphoglycerate dehydrogenase', 'gene_type': 'protein_coding', 'chromosome': '3', 'aliases': []},
    'Psat1': {'ncbi_gene_id': '107280', 'ensembl_id': 'ENSMUSG00000029238', 'gene_name': 'phosphoserine aminotransferase 1', 'gene_type': 'protein_coding', 'chromosome': '19', 'aliases': []},
    'Psph': {'ncbi_gene_id': '100080', 'ensembl_id': 'ENSMUSG00000026242', 'gene_name': 'phosphoserine phosphatase', 'gene_type': 'protein_coding', 'chromosome': '11', 'aliases': []},
    'Shmt1': {'ncbi_gene_id': '20485', 'ensembl_id': 'ENSMUSG00000020543', 'gene_name': 'serine hydroxymethyltransferase 1', 'gene_type': 'protein_coding', 'chromosome': '12', 'aliases': []},
    'Shmt2': {'ncbi_gene_id': '108091', 'ensembl_id': 'ENSMUSG00000003345', 'gene_name': 'serine hydroxymethyltransferase 2', 'gene_type': 'protein_coding', 'chromosome': '5', 'aliases': []},
    'Gldc': {'ncbi_gene_id': '105497', 'ensembl_id': 'ENSMUSG00000024827', 'gene_name': 'glycine decarboxylase', 'gene_type': 'protein_coding', 'chromosome': '4', 'aliases': []},
    'Amt': {'ncbi_gene_id': '10680', 'ensembl_id': 'ENSMUSG00000032607', 'gene_name': 'aminomethyltransferase', 'gene_type': 'protein_coding', 'chromosome': '9', 'aliases': []},
    'Gpd1': {'ncbi_gene_id': '14555', 'ensembl_id': 'ENSMUSG00000024162', 'gene_name': 'glycerol-3-phosphate dehydrogenase 1 (cytosolic)', 'gene_type': 'protein_coding', 'chromosome': '15', 'aliases': []},
    'Gpd2': {'ncbi_gene_id': '14556', 'ensembl_id': 'ENSMUSG00000029991', 'gene_name': 'glycerol-3-phosphate dehydrogenase 2 (mitochondrial)', 'gene_type': 'protein_coding', 'chromosome': '6', 'aliases': []},
    'Gk': {'ncbi_gene_id': '14933', 'ensembl_id': 'ENSMUSG00000024077', 'gene_name': 'glycerol kinase', 'gene_type': 'protein_coding', 'chromosome': '7', 'aliases': []},
    'Gk2': {'ncbi_gene_id': '243775', 'ensembl_id': 'ENSMUSG00000041141', 'gene_name': 'glycerol kinase 2', 'gene_type': 'protein_coding', 'chromosome': '5', 'aliases': []},
    'Gyk': {'ncbi_gene_id': '14935', 'ensembl_id': 'ENSMUSG00000031483', 'gene_name': 'glycerol kinase', 'gene_type': 'protein_coding', 'chromosome': '5', 'aliases': []},
    'Gpam': {'ncbi_gene_id': '14710', 'ensembl_id': 'ENSMUSG00000027470', 'gene_name': 'glycerol-3-phosphate acyltransferase, mitochondrial', 'gene_type': 'protein_coding', 'chromosome': '19', 'aliases': []},
    'Agpat1': {'ncbi_gene_id': '56161', 'ensembl_id': 'ENSMUSG00000034254', 'gene_name': '1-acylglycerol-3-phosphate O-acyltransferase 1', 'gene_type': 'protein_coding', 'chromosome': '9', 'aliases': []},
    'Agpat2': {'ncbi_gene_id': '66101', 'ensembl_id': 'ENSMUSG00000031776', 'gene_name': '1-acylglycerol-3-phosphate O-acyltransferase 2', 'gene_type': 'protein_coding', 'chromosome': '11', 'aliases': []},
    'Agpat3': {'ncbi_gene_id': '11874', 'ensembl_id': 'ENSMUSG00000029373', 'gene_name': '1-acylglycerol-3-phosphate O-acyltransferase 3', 'gene_type': 'protein_coding', 'chromosome': '1', 'aliases': []},
    'Dgat1': {'ncbi_gene_id': '14574', 'ensembl_id': 'ENSMUSG00000022555', 'gene_name': 'diacylglycerol O-acyltransferase 1', 'gene_type': 'protein_coding', 'chromosome': '17', 'aliases': []},
    'Dgat2': {'ncbi_gene_id': '9451', 'ensembl_id': 'ENSMUSG00000022636', 'gene_name': 'diacylglycerol O-acyltransferase 2', 'gene_type': 'protein_coding', 'chromosome': '7', 'aliases': []},
    'Pnpla2': {'ncbi_gene_id': '66853', 'ensembl_id': 'ENSMUSG00000025509', 'gene_name': 'patatin-like phospholipase domain containing 2 (ATGL)', 'gene_type': 'protein_coding', 'chromosome': '7', 'aliases': []},
    'Lipe': {'ncbi_gene_id': '16890', 'ensembl_id': 'ENSMUSG00000055980', 'gene_name': 'lipase E, hormone sensitive', 'gene_type': 'protein_coding', 'chromosome': '5', 'aliases': []},
    'Mgll': {'ncbi_gene_id': '23945', 'ensembl_id': 'ENSMUSG00000029683', 'gene_name': 'monoglyceride lipase', 'gene_type': 'protein_coding', 'chromosome': '6', 'aliases': []},
    'Lpl': {'ncbi_gene_id': '16956', 'ensembl_id': 'ENSMUSG00000029285', 'gene_name': 'lipoprotein lipase', 'gene_type': 'protein_coding', 'chromosome': '8', 'aliases': []},
    'Hlcs': {'ncbi_gene_id': '110948', 'ensembl_id': 'ENSMUSG00000024466', 'gene_name': 'holocarboxylase synthetase', 'gene_type': 'protein_coding', 'chromosome': '16', 'aliases': []},
    'Mcc1': {'ncbi_gene_id': '21325', 'ensembl_id': 'ENSMUSG00000029681', 'gene_name': 'methylcrotonoyl-CoA carboxylase 1 (alpha)', 'gene_type': 'protein_coding', 'chromosome': '3', 'aliases': []},
    'Chka': {'ncbi_gene_id': '12660', 'ensembl_id': 'ENSMUSG00000024843', 'gene_name': 'choline kinase alpha', 'gene_type': 'protein_coding', 'chromosome': '19', 'aliases': []},
    'Chkb': {'ncbi_gene_id': '21139', 'ensembl_id': 'ENSMUSG00000028716', 'gene_name': 'choline kinase beta', 'gene_type': 'protein_coding', 'chromosome': '7', 'aliases': []},
    'Pcyt1a': {'ncbi_gene_id': '16432', 'ensembl_id': 'ENSMUSG00000033565', 'gene_name': 'phosphate cytidylyltransferase 1, choline, alpha isoform', 'gene_type': 'protein_coding', 'chromosome': '16', 'aliases': []},
    'Pcyt2': {'ncbi_gene_id': '16972', 'ensembl_id': 'ENSMUSG00000038546', 'gene_name': 'phosphate cytidylyltransferase 2, ethanolamine', 'gene_type': 'protein_coding', 'chromosome': '9', 'aliases': []},
    'Fabp1': {'ncbi_gene_id': '14080', 'ensembl_id': 'ENSMUSG00000024750', 'gene_name': 'fatty acid binding protein 1, liver', 'gene_type': 'protein_coding', 'chromosome': '6', 'aliases': []},
    'Fabp2': {'ncbi_gene_id': '14079', 'ensembl_id': 'ENSMUSG00000023057', 'gene_name': 'fatty acid binding protein 2, intestinal', 'gene_type': 'protein_coding', 'chromosome': '3', 'aliases': []},
    'Fabp3': {'ncbi_gene_id': '14077', 'ensembl_id': 'ENSMUSG00000028795', 'gene_name': 'fatty acid binding protein 3, muscle and heart', 'gene_type': 'protein_coding', 'chromosome': '4', 'aliases': []},
    'Fabp4': {'ncbi_gene_id': '11770', 'ensembl_id': 'ENSMUSG00000062515', 'gene_name': 'fatty acid binding protein 4, adipocyte', 'gene_type': 'protein_coding', 'chromosome': '3', 'aliases': []},
    'Fabp5': {'ncbi_gene_id': '16592', 'ensembl_id': 'ENSMUSG00000027222', 'gene_name': 'fatty acid binding protein 5, epidermal', 'gene_type': 'protein_coding', 'chromosome': '3', 'aliases': []},
    'Fabp6': {'ncbi_gene_id': '11769', 'ensembl_id': 'ENSMUSG00000015023', 'gene_name': 'fatty acid binding protein 6, ileal', 'gene_type': 'protein_coding', 'chromosome': '16', 'aliases': []},
    'Fabp7': {'ncbi_gene_id': '12140', 'ensembl_id': 'ENSMUSG00000046132', 'gene_name': 'fatty acid binding protein 7, brain', 'gene_type': 'protein_coding', 'chromosome': '1', 'aliases': []},
    'Cd36': {'ncbi_gene_id': '12491', 'ensembl_id': 'ENSMUSG00000022944', 'gene_name': 'CD36 molecule', 'gene_type': 'protein_coding', 'chromosome': '5', 'aliases': []},
    'Slc2a1': {'ncbi_gene_id': '20525', 'ensembl_id': 'ENSMUSG00000028645', 'gene_name': 'solute carrier family 2 member 1 (GLUT1)', 'gene_type': 'protein_coding', 'chromosome': '4', 'aliases': []},
    'Slc2a2': {'ncbi_gene_id': '20526', 'ensembl_id': 'ENSMUSG00000027690', 'gene_name': 'solute carrier family 2 member 2 (GLUT2)', 'gene_type': 'protein_coding', 'chromosome': '3', 'aliases': []},
    'Slc2a3': {'ncbi_gene_id': '20527', 'ensembl_id': 'ENSMUSG00000028782', 'gene_name': 'solute carrier family 2 member 3 (GLUT3)', 'gene_type': 'protein_coding', 'chromosome': '4', 'aliases': []},
    'Slc2a4': {'ncbi_gene_id': '20528', 'ensembl_id': 'ENSMUSG00000028639', 'gene_name': 'solute carrier family 2 member 4 (GLUT4)', 'gene_type': 'protein_coding', 'chromosome': '3', 'aliases': []},
    'Slc2a5': {'ncbi_gene_id': '20529', 'ensembl_id': 'ENSMUSG00000031750', 'gene_name': 'solute carrier family 2 member 5 (GLUT5)', 'gene_type': 'protein_coding', 'chromosome': '10', 'aliases': []},
    'Slc2a8': {'ncbi_gene_id': '20530', 'ensembl_id': 'ENSMUSG00000035933', 'gene_name': 'solute carrier family 2 member 8 (GLUT8)', 'gene_type': 'protein_coding', 'chromosome': '4', 'aliases': []},
    'Slc5a1': {'ncbi_gene_id': '20460', 'ensembl_id': 'ENSMUSG00000033855', 'gene_name': 'solute carrier family 5 member 1 (SGLT1)', 'gene_type': 'protein_coding', 'chromosome': '7', 'aliases': []},
    'Slc5a2': {'ncbi_gene_id': '20461', 'ensembl_id': 'ENSMUSG00000033856', 'gene_name': 'solute carrier family 5 member 2 (SGLT2)', 'gene_type': 'protein_coding', 'chromosome': '7', 'aliases': []},
    'Slc5a3': {'ncbi_gene_id': '20536', 'ensembl_id': 'ENSMUSG00000035037', 'gene_name': 'solute carrier family 5 member 3 (SMIT1)', 'gene_type': 'protein_coding', 'chromosome': '16', 'aliases': []},
    'Slc16a1': {'ncbi_gene_id': '20503', 'ensembl_id': 'ENSMUSG00000026836', 'gene_name': 'solute carrier family 16 member 1 (MCT1)', 'gene_type': 'protein_coding', 'chromosome': '3', 'aliases': []},
    'Slc16a3': {'ncbi_gene_id': '80879', 'ensembl_id': 'ENSMUSG00000024456', 'gene_name': 'solute carrier family 16 member 3 (MCT4)', 'gene_type': 'protein_coding', 'chromosome': '8', 'aliases': []},
    'Slc16a6': {'ncbi_gene_id': '10464', 'ensembl_id': 'ENSMUSG00000036658', 'gene_name': 'solute carrier family 16 member 6', 'gene_type': 'protein_coding', 'chromosome': '11', 'aliases': []},
    'Slc16a7': {'ncbi_gene_id': '20504', 'ensembl_id': 'ENSMUSG00000037015', 'gene_name': 'solute carrier family 16 member 7 (MCT2)', 'gene_type': 'protein_coding', 'chromosome': '10', 'aliases': []},
    'Slc1a1': {'ncbi_gene_id': '20510', 'ensembl_id': 'ENSMUSG00000024935', 'gene_name': 'solute carrier family 1 member 1 (EAAT3)', 'gene_type': 'protein_coding', 'chromosome': '19', 'aliases': []},
    'Slc1a2': {'ncbi_gene_id': '20511', 'ensembl_id': 'ENSMUSG00000005089', 'gene_name': 'solute carrier family 1 member 2 (EAAT2 / GLT-1)', 'gene_type': 'protein_coding', 'chromosome': '2', 'aliases': []},
    'Slc1a3': {'ncbi_gene_id': '20512', 'ensembl_id': 'ENSMUSG00000005268', 'gene_name': 'solute carrier family 1 member 3 (EAAT1 / GLAST)', 'gene_type': 'protein_coding', 'chromosome': '15', 'aliases': []},
    'Slc1a4': {'ncbi_gene_id': '55963', 'ensembl_id': 'ENSMUSG00000020142', 'gene_name': 'solute carrier family 1 member 4 (ASCT1)', 'gene_type': 'protein_coding', 'chromosome': '11', 'aliases': []},
    'Slc1a5': {'ncbi_gene_id': '20514', 'ensembl_id': 'ENSMUSG00000001920', 'gene_name': 'solute carrier family 1 member 5 (ASCT2)', 'gene_type': 'protein_coding', 'chromosome': '7', 'aliases': []},
    'Slc1a6': {'ncbi_gene_id': '20515', 'ensembl_id': 'ENSMUSG00000005357', 'gene_name': 'solute carrier family 1 member 6 (EAAT4)', 'gene_type': 'protein_coding', 'chromosome': '10', 'aliases': []},
    'Slc1a7': {'ncbi_gene_id': '230982', 'ensembl_id': 'ENSMUSG00000035803', 'gene_name': 'solute carrier family 1 member 7 (EAAT5)', 'gene_type': 'protein_coding', 'chromosome': '4', 'aliases': []},
    'Slc7a1': {'ncbi_gene_id': '17874', 'ensembl_id': 'ENSMUSG00000041313', 'gene_name': 'solute carrier family 7 member 1 (CAT1)', 'gene_type': 'protein_coding', 'chromosome': '17', 'aliases': []},
    'Slc7a2': {'ncbi_gene_id': '17875', 'ensembl_id': 'ENSMUSG00000028655', 'gene_name': 'solute carrier family 7 member 2 (CAT2)', 'gene_type': 'protein_coding', 'chromosome': '15', 'aliases': []},
    'Slc7a5': {'ncbi_gene_id': '20539', 'ensembl_id': 'ENSMUSG00000007576', 'gene_name': 'solute carrier family 7 member 5 (LAT1)', 'gene_type': 'protein_coding', 'chromosome': '5', 'aliases': []},
    'Slc7a11': {'ncbi_gene_id': '26570', 'ensembl_id': 'ENSMUSG00000029204', 'gene_name': 'solute carrier family 7 member 11 (xCT)', 'gene_type': 'protein_coding', 'chromosome': '7', 'aliases': []},
    'Slc38a1': {'ncbi_gene_id': '56362', 'ensembl_id': 'ENSMUSG00000029001', 'gene_name': 'solute carrier family 38 member 1 (SNAT1)', 'gene_type': 'protein_coding', 'chromosome': '6', 'aliases': []},
    'Slc38a2': {'ncbi_gene_id': '56363', 'ensembl_id': 'ENSMUSG00000029183', 'gene_name': 'solute carrier family 38 member 2 (SNAT2)', 'gene_type': 'protein_coding', 'chromosome': '11', 'aliases': []},
    'Slc38a3': {'ncbi_gene_id': '22325', 'ensembl_id': 'ENSMUSG00000041856', 'gene_name': 'solute carrier family 38 member 3 (SNAT3)', 'gene_type': 'protein_coding', 'chromosome': '16', 'aliases': []},
    'Slc38a5': {'ncbi_gene_id': '105727', 'ensembl_id': 'ENSMUSG00000031163', 'gene_name': 'solute carrier family 38 member 5 (SNAT5)', 'gene_type': 'protein_coding', 'chromosome': '5', 'aliases': []},
    'Slc3a2': {'ncbi_gene_id': '17254', 'ensembl_id': 'ENSMUSG00000010097', 'gene_name': 'solute carrier family 3 member 2 (CD98hc / 4F2hc)', 'gene_type': 'protein_coding', 'chromosome': '19', 'aliases': []},
    'Slc6a9': {'ncbi_gene_id': '14664', 'ensembl_id': 'ENSMUSG00000044336', 'gene_name': 'solute carrier family 6 member 9 (GlyT1)', 'gene_type': 'protein_coding', 'chromosome': '4', 'aliases': []},
    'Slc25a1': {'ncbi_gene_id': '55997', 'ensembl_id': 'ENSMUSG00000006201', 'gene_name': 'solute carrier family 25 member 1 (CIC)', 'gene_type': 'protein_coding', 'chromosome': '16', 'aliases': []},
    'Slc25a10': {'ncbi_gene_id': '56663', 'ensembl_id': 'ENSMUSG00000025408', 'gene_name': 'solute carrier family 25 member 10 (DIC)', 'gene_type': 'protein_coding', 'chromosome': '11', 'aliases': []},
    'Slc25a11': {'ncbi_gene_id': '67841', 'ensembl_id': 'ENSMUSG00000021417', 'gene_name': 'solute carrier family 25 member 11 (OGC)', 'gene_type': 'protein_coding', 'chromosome': '10', 'aliases': []},
    'Slc25a12': {'ncbi_gene_id': '14636', 'ensembl_id': 'ENSMUSG00000024690', 'gene_name': 'solute carrier family 25 member 12 (Aralar1)', 'gene_type': 'protein_coding', 'chromosome': '6', 'aliases': []},
    'Slc25a13': {'ncbi_gene_id': '13932', 'ensembl_id': 'ENSMUSG00000024073', 'gene_name': 'solute carrier family 25 member 13 (Citrin)', 'gene_type': 'protein_coding', 'chromosome': '2', 'aliases': []},
    'Slc27a1': {'ncbi_gene_id': '26457', 'ensembl_id': 'ENSMUSG00000031808', 'gene_name': 'solute carrier family 27 member 1 (FATP1)', 'gene_type': 'protein_coding', 'chromosome': '7', 'aliases': []},
    'Slc27a2': {'ncbi_gene_id': '26458', 'ensembl_id': 'ENSMUSG00000027812', 'gene_name': 'solute carrier family 27 member 2 (FATP2)', 'gene_type': 'protein_coding', 'chromosome': '2', 'aliases': []},
    'Slc27a3': {'ncbi_gene_id': '74152', 'ensembl_id': 'ENSMUSG00000027935', 'gene_name': 'solute carrier family 27 member 3 (FATP3)', 'gene_type': 'protein_coding', 'chromosome': '7', 'aliases': []},
    'Slc27a4': {'ncbi_gene_id': '56729', 'ensembl_id': 'ENSMUSG00000038853', 'gene_name': 'solute carrier family 27 member 4 (FATP4)', 'gene_type': 'protein_coding', 'chromosome': '3', 'aliases': []},
    'Slc27a5': {'ncbi_gene_id': '26459', 'ensembl_id': 'ENSMUSG00000041005', 'gene_name': 'solute carrier family 27 member 5 (FATP5)', 'gene_type': 'protein_coding', 'chromosome': '8', 'aliases': []},
    'Slc27a6': {'ncbi_gene_id': '225579', 'ensembl_id': 'ENSMUSG00000024600', 'gene_name': 'solute carrier family 27 member 6 (FATP6)', 'gene_type': 'protein_coding', 'chromosome': '18', 'aliases': []},
    'Slc22a5': {'ncbi_gene_id': '20519', 'ensembl_id': 'ENSMUSG00000021638', 'gene_name': 'solute carrier family 22 member 5 (OCTN2)', 'gene_type': 'protein_coding', 'chromosome': '10', 'aliases': []},
    'Slc37a1': {'ncbi_gene_id': '22467', 'ensembl_id': 'ENSMUSG00000023787', 'gene_name': 'solute carrier family 37 member 1', 'gene_type': 'protein_coding', 'chromosome': '14', 'aliases': []},
    'Slc37a2': {'ncbi_gene_id': '106068', 'ensembl_id': 'ENSMUSG00000032122', 'gene_name': 'solute carrier family 37 member 2', 'gene_type': 'protein_coding', 'chromosome': '9', 'aliases': []},
    'Slc37a4': {'ncbi_gene_id': '21564', 'ensembl_id': 'ENSMUSG00000040136', 'gene_name': 'solute carrier family 37 member 4 (G6PT1)', 'gene_type': 'protein_coding', 'chromosome': '21', 'aliases': []},
    'Slc5a8': {'ncbi_gene_id': '69005', 'ensembl_id': 'ENSMUSG00000031642', 'gene_name': 'solute carrier family 5 member 8 (SMCT1)', 'gene_type': 'protein_coding', 'chromosome': '10', 'aliases': []},
    'Slc4a1': {'ncbi_gene_id': '20393', 'ensembl_id': 'ENSMUSG00000029006', 'gene_name': 'solute carrier family 4 member 1 (AE1 / Band 3)', 'gene_type': 'protein_coding', 'chromosome': '11', 'aliases': []},
    'Aqp3': {'ncbi_gene_id': '11828', 'ensembl_id': 'ENSMUSG00000028435', 'gene_name': 'aquaporin 3', 'gene_type': 'protein_coding', 'chromosome': '5', 'aliases': []},
    'Aqp7': {'ncbi_gene_id': '11832', 'ensembl_id': 'ENSMUSG00000030957', 'gene_name': 'aquaporin 7', 'gene_type': 'protein_coding', 'chromosome': '4', 'aliases': []},
    'Aqp9': {'ncbi_gene_id': '32901', 'ensembl_id': 'ENSMUSG00000032204', 'gene_name': 'aquaporin 9', 'gene_type': 'protein_coding', 'chromosome': '9', 'aliases': []},
    'Aqp10': {'ncbi_gene_id': '57140', 'ensembl_id': 'ENSMUSG00000039997', 'gene_name': 'aquaporin 10', 'gene_type': 'protein_coding', 'chromosome': '3', 'aliases': []},
    'Slc17a6': {'ncbi_gene_id': '140919', 'ensembl_id': 'ENSMUSG00000031618', 'gene_name': 'solute carrier family 17 member 6 (VGLUT2)', 'gene_type': 'protein_coding', 'chromosome': '10', 'aliases': []},
    'Slc17a7': {'ncbi_gene_id': '109838', 'ensembl_id': 'ENSMUSG00000045376', 'gene_name': 'solute carrier family 17 member 7 (VGLUT1)', 'gene_type': 'protein_coding', 'chromosome': '10', 'aliases': []},
    'Slc17a8': {'ncbi_gene_id': '140920', 'ensembl_id': 'ENSMUSG00000041622', 'gene_name': 'solute carrier family 17 member 8 (VGLUT3)', 'gene_type': 'protein_coding', 'chromosome': '5', 'aliases': []},
    'Insr': {'ncbi_gene_id': '16337', 'ensembl_id': 'ENSMUSG00000005534', 'gene_name': 'insulin receptor', 'gene_type': 'protein_coding', 'chromosome': '17', 'aliases': []},
    'Irs1': {'ncbi_gene_id': '16367', 'ensembl_id': 'ENSMUSG00000032711', 'gene_name': 'insulin receptor substrate 1', 'gene_type': 'protein_coding', 'chromosome': '14', 'aliases': []},
    'Irs2': {'ncbi_gene_id': '384783', 'ensembl_id': 'ENSMUSG00000038894', 'gene_name': 'insulin receptor substrate 2', 'gene_type': 'protein_coding', 'chromosome': '9', 'aliases': []},
    'Pik3r1': {'ncbi_gene_id': '18708', 'ensembl_id': 'ENSMUSG00000024368', 'gene_name': 'phosphoinositide-3-kinase regulatory subunit 1', 'gene_type': 'protein_coding', 'chromosome': '4', 'aliases': []},
    'Akt1': {'ncbi_gene_id': '11651', 'ensembl_id': 'ENSMUSG00000001729', 'gene_name': 'thymoma viral proto-oncogene (AKT serine/threonine kinase 1)', 'gene_type': 'protein_coding', 'chromosome': '12', 'aliases': []},
    'Akt2': {'ncbi_gene_id': '11652', 'ensembl_id': 'ENSMUSG00000000666', 'gene_name': 'thymoma viral proto-oncogene (AKT serine/threonine kinase 2)', 'gene_type': 'protein_coding', 'chromosome': '7', 'aliases': []},
    'Foxo1': {'ncbi_gene_id': '56458', 'ensembl_id': 'ENSMUSG00000044167', 'gene_name': 'forkhead box O1', 'gene_type': 'protein_coding', 'chromosome': '4', 'aliases': []},
    'Mtor': {'ncbi_gene_id': '56717', 'ensembl_id': 'ENSMUSG00000028991', 'gene_name': 'mechanistic target of rapamycin kinase', 'gene_type': 'protein_coding', 'chromosome': '12', 'aliases': []},
    'Rptor': {'ncbi_gene_id': '74370', 'ensembl_id': 'ENSMUSG00000027875', 'gene_name': 'regulatory associated protein of MTOR, complex 1', 'gene_type': 'protein_coding', 'chromosome': '12', 'aliases': []},
    'Rheb': {'ncbi_gene_id': '19739', 'ensembl_id': 'ENSMUSG00000028945', 'gene_name': 'Ras homolog enriched in brain', 'gene_type': 'protein_coding', 'chromosome': '1', 'aliases': []},
    'Mlst8': {'ncbi_gene_id': '105849', 'ensembl_id': 'ENSMUSG00000054941', 'gene_name': 'MTOR associated protein, LST8 homolog', 'gene_type': 'protein_coding', 'chromosome': '17', 'aliases': []},
    'Deptor': {'ncbi_gene_id': '58963', 'ensembl_id': 'ENSMUSG00000025348', 'gene_name': 'DEP domain containing MTOR-interacting protein', 'gene_type': 'protein_coding', 'chromosome': '15', 'aliases': []},
    'Tsc1': {'ncbi_gene_id': '64930', 'ensembl_id': 'ENSMUSG00000031714', 'gene_name': 'TSC complex subunit 1', 'gene_type': 'protein_coding', 'chromosome': '11', 'aliases': []},
    'Tsc2': {'ncbi_gene_id': '22084', 'ensembl_id': 'ENSMUSG00000038953', 'gene_name': 'TSC complex subunit 2', 'gene_type': 'protein_coding', 'chromosome': '7', 'aliases': []},
    'Rraga': {'ncbi_gene_id': '23370', 'ensembl_id': 'ENSMUSG00000031867', 'gene_name': 'Ras-related GTP binding A', 'gene_type': 'protein_coding', 'chromosome': '9', 'aliases': []},
    'Rragb': {'ncbi_gene_id': '73858', 'ensembl_id': 'ENSMUSG00000022090', 'gene_name': 'Ras-related GTP binding B', 'gene_type': 'protein_coding', 'chromosome': '14', 'aliases': []},
    'Rragc': {'ncbi_gene_id': '26986', 'ensembl_id': 'ENSMUSG00000035846', 'gene_name': 'Ras-related GTP binding C', 'gene_type': 'protein_coding', 'chromosome': '10', 'aliases': []},
    'Rragd': {'ncbi_gene_id': '29197', 'ensembl_id': 'ENSMUSG00000025304', 'gene_name': 'Ras-related GTP binding D', 'gene_type': 'protein_coding', 'chromosome': '9', 'aliases': []},
    'Igf1': {'ncbi_gene_id': '16000', 'ensembl_id': 'ENSMUSG00000020053', 'gene_name': 'insulin-like growth factor 1', 'gene_type': 'protein_coding', 'chromosome': '16', 'aliases': []},
    'Igf1r': {'ncbi_gene_id': '16001', 'ensembl_id': 'ENSMUSG00000005993', 'gene_name': 'insulin-like growth factor I receptor', 'gene_type': 'protein_coding', 'chromosome': '11', 'aliases': []},
    'Egfr': {'ncbi_gene_id': '13759', 'ensembl_id': 'ENSMUSG00000020122', 'gene_name': 'epidermal growth factor receptor', 'gene_type': 'protein_coding', 'chromosome': '9', 'aliases': []},
    'Hif1a': {'ncbi_gene_id': '15251', 'ensembl_id': 'ENSMUSG00000021109', 'gene_name': 'hypoxia inducible factor 1 alpha subunit', 'gene_type': 'protein_coding', 'chromosome': '19', 'aliases': []},
    'Tp53': {'ncbi_gene_id': '22059', 'ensembl_id': 'ENSMUSG00000059552', 'gene_name': 'transformation related protein 53 (p53)', 'gene_type': 'protein_coding', 'chromosome': '11', 'aliases': []},
    'Myc': {'ncbi_gene_id': '17869', 'ensembl_id': 'ENSMUSG00000022346', 'gene_name': 'myelocytomatosis oncogene (c-Myc)', 'gene_type': 'protein_coding', 'chromosome': '15', 'aliases': []},
    'Mycn': {'ncbi_gene_id': '18109', 'ensembl_id': 'ENSMUSG00000037169', 'gene_name': 'v-myc myelocytomatosis viral related oncogene, neuroblastoma derived', 'gene_type': 'protein_coding', 'chromosome': '17', 'aliases': []},
    'Mlxipl': {'ncbi_gene_id': '58805', 'ensembl_id': 'ENSMUSG00000004022', 'gene_name': 'MLX interacting protein like (ChREBP)', 'gene_type': 'protein_coding', 'chromosome': '5', 'aliases': []},
    'Mlx': {'ncbi_gene_id': '21428', 'ensembl_id': 'ENSMUSG00000017801', 'gene_name': 'MAX-like protein X', 'gene_type': 'protein_coding', 'chromosome': '9', 'aliases': []},
    'Srebf1': {'ncbi_gene_id': '20787', 'ensembl_id': 'ENSMUSG00000020538', 'gene_name': 'sterol regulatory element binding transcription factor 1', 'gene_type': 'protein_coding', 'chromosome': '11', 'aliases': []},
    'Srebf2': {'ncbi_gene_id': '20788', 'ensembl_id': 'ENSMUSG00000021880', 'gene_name': 'sterol regulatory element binding transcription factor 2', 'gene_type': 'protein_coding', 'chromosome': '11', 'aliases': []},
    'Ppara': {'ncbi_gene_id': '19013', 'ensembl_id': 'ENSMUSG00000022383', 'gene_name': 'peroxisome proliferator activated receptor alpha', 'gene_type': 'protein_coding', 'chromosome': '18', 'aliases': []},
    'Ppard': {'ncbi_gene_id': '19015', 'ensembl_id': 'ENSMUSG00000025106', 'gene_name': 'peroxisome proliferator activated receptor delta', 'gene_type': 'protein_coding', 'chromosome': '17', 'aliases': []},
    'Pparg': {'ncbi_gene_id': '19016', 'ensembl_id': 'ENSMUSG00000000349', 'gene_name': 'peroxisome proliferator activated receptor gamma', 'gene_type': 'protein_coding', 'chromosome': '9', 'aliases': []},
    'Ppargc1a': {'ncbi_gene_id': '19017', 'ensembl_id': 'ENSMUSG00000029167', 'gene_name': 'peroxisome proliferative activated receptor gamma coactivator 1 alpha (PGC-1a)', 'gene_type': 'protein_coding', 'chromosome': '4', 'aliases': []},
    'Ppargc1b': {'ncbi_gene_id': '170826', 'ensembl_id': 'ENSMUSG00000031916', 'gene_name': 'peroxisome proliferative activated receptor gamma coactivator 1 beta (PGC-1b)', 'gene_type': 'protein_coding', 'chromosome': '18', 'aliases': []},
    'Nr1h3': {'ncbi_gene_id': '22259', 'ensembl_id': 'ENSMUSG00000025868', 'gene_name': 'nuclear receptor subfamily 1 group H member 3 (LXR-alpha)', 'gene_type': 'protein_coding', 'chromosome': '15', 'aliases': []},
    'Nr1h4': {'ncbi_gene_id': '20186', 'ensembl_id': 'ENSMUSG00000025937', 'gene_name': 'nuclear receptor subfamily 1 group H member 4 (FXR)', 'gene_type': 'protein_coding', 'chromosome': '1', 'aliases': []},
    'Hnf4a': {'ncbi_gene_id': '15378', 'ensembl_id': 'ENSMUSG00000017950', 'gene_name': 'hepatic nuclear factor 4 alpha', 'gene_type': 'protein_coding', 'chromosome': '11', 'aliases': []},
    'Creb1': {'ncbi_gene_id': '12912', 'ensembl_id': 'ENSMUSG00000025959', 'gene_name': 'cAMP responsive element binding protein 1', 'gene_type': 'protein_coding', 'chromosome': '9', 'aliases': []},
    'Crtc2': {'ncbi_gene_id': '67955', 'ensembl_id': 'ENSMUSG00000027936', 'gene_name': 'CREB regulated transcription coactivator 2', 'gene_type': 'protein_coding', 'chromosome': '3', 'aliases': []},
    'Nrf1': {'ncbi_gene_id': '18188', 'ensembl_id': 'ENSMUSG00000001586', 'gene_name': 'nuclear respiratory factor 1', 'gene_type': 'protein_coding', 'chromosome': '8', 'aliases': []},
    'Tfam': {'ncbi_gene_id': '21780', 'ensembl_id': 'ENSMUSG00000003923', 'gene_name': 'transcription factor A, mitochondrial', 'gene_type': 'protein_coding', 'chromosome': '9', 'aliases': []},
    'Zeb1': {'ncbi_gene_id': '21417', 'ensembl_id': 'ENSMUSG00000036540', 'gene_name': 'zinc finger E-box binding homeobox 1', 'gene_type': 'protein_coding', 'chromosome': '5', 'aliases': []},
    'Atf4': {'ncbi_gene_id': '11911', 'ensembl_id': 'ENSMUSG00000042406', 'gene_name': 'activating transcription factor 4', 'gene_type': 'protein_coding', 'chromosome': '15', 'aliases': []},
    'Prkaa1': {'ncbi_gene_id': '105787', 'ensembl_id': 'ENSMUSG00000050697', 'gene_name': 'protein kinase AMP-activated catalytic subunit alpha 1', 'gene_type': 'protein_coding', 'chromosome': '7', 'aliases': []},
    'Prkaa2': {'ncbi_gene_id': '108079', 'ensembl_id': 'ENSMUSG00000028518', 'gene_name': 'protein kinase AMP-activated catalytic subunit alpha 2', 'gene_type': 'protein_coding', 'chromosome': '7', 'aliases': []},
    'Prkag1': {'ncbi_gene_id': '19082', 'ensembl_id': 'ENSMUSG00000005956', 'gene_name': 'protein kinase AMP-activated non-catalytic subunit gamma 1', 'gene_type': 'protein_coding', 'chromosome': '11', 'aliases': []},
    'Egln1': {'ncbi_gene_id': '17068', 'ensembl_id': 'ENSMUSG00000030846', 'gene_name': 'egl-9 family hypoxia-inducible factor 1 (PHD2)', 'gene_type': 'protein_coding', 'chromosome': '11', 'aliases': []},
    'Egln2': {'ncbi_gene_id': '17070', 'ensembl_id': 'ENSMUSG00000031796', 'gene_name': 'egl-9 family hypoxia-inducible factor 1a (PHD1)', 'gene_type': 'protein_coding', 'chromosome': '10', 'aliases': []},
    'Vhl': {'ncbi_gene_id': '22346', 'ensembl_id': 'ENSMUSG00000031659', 'gene_name': 'von Hippel-Lindau tumor suppressor (pVHL)', 'gene_type': 'protein_coding', 'chromosome': '17', 'aliases': []},
    'Hnrnpa1': {'ncbi_gene_id': '15382', 'ensembl_id': 'ENSMUSG00000020379', 'gene_name': 'heterogeneous nuclear ribonucleoprotein A1', 'gene_type': 'protein_coding', 'chromosome': '10', 'aliases': []},
    'Hnrnpa2b1': {'ncbi_gene_id': '53379', 'ensembl_id': 'ENSMUSG00000090633', 'gene_name': 'heterogeneous nuclear ribonucleoprotein A2/B1', 'gene_type': 'protein_coding', 'chromosome': '10', 'aliases': []},
    'Ptbp1': {'ncbi_gene_id': '19205', 'ensembl_id': 'ENSMUSG00000025158', 'gene_name': 'polypyrimidine tract binding protein 1', 'gene_type': 'protein_coding', 'chromosome': '7', 'aliases': []},
    'Srsf1': {'ncbi_gene_id': '110809', 'ensembl_id': 'ENSMUSG00000042367', 'gene_name': 'serine and arginine-rich splicing factor 1', 'gene_type': 'protein_coding', 'chromosome': '4', 'aliases': []},
    'Srsf3': {'ncbi_gene_id': '20383', 'ensembl_id': 'ENSMUSG00000028660', 'gene_name': 'serine and arginine-rich splicing factor 3', 'gene_type': 'protein_coding', 'chromosome': '10', 'aliases': []},
    'Eif4e': {'ncbi_gene_id': '13684', 'ensembl_id': 'ENSMUSG00000028156', 'gene_name': 'eukaryotic translation initiation factor 4E', 'gene_type': 'protein_coding', 'chromosome': '11', 'aliases': []},
    'Eif4ebp1': {'ncbi_gene_id': '13685', 'ensembl_id': 'ENSMUSG00000031428', 'gene_name': 'eukaryotic translation initiation factor 4E binding protein 1', 'gene_type': 'protein_coding', 'chromosome': '12', 'aliases': []},
    'Eif2ak4': {'ncbi_gene_id': '27103', 'ensembl_id': 'ENSMUSG00000029110', 'gene_name': 'eukaryotic translation initiation factor 2 alpha kinase 4 (GCN2)', 'gene_type': 'protein_coding', 'chromosome': '10', 'aliases': []},
    'Gde1': {'ncbi_gene_id': '56209', 'ensembl_id': 'ENSMUSG00000033917', 'gene_name': 'glycerophosphodiester phosphodiesterase 1', 'gene_type': 'protein_coding', 'chromosome': '7', 'aliases': []},
    'Agl': {'ncbi_gene_id': '14634', 'ensembl_id': 'ENSMUSG00000020209', 'gene_name': 'amylo-1,6-glucosidase, 4-alpha-glucanotransferase', 'gene_type': 'protein_coding', 'chromosome': '1', 'aliases': []},
    'Cad': {'ncbi_gene_id': '22258', 'ensembl_id': 'ENSMUSG00000032237', 'gene_name': 'carbamoyl-phosphate synthetase 2, aspartate transcarbamylase, and dihydroorotase', 'gene_type': 'protein_coding', 'chromosome': '16', 'aliases': []},
    'Ppat': {'ncbi_gene_id': '231327', 'ensembl_id': 'ENSMUSG00000042844', 'gene_name': 'phosphoribosyl pyrophosphate amidotransferase', 'gene_type': 'protein_coding', 'chromosome': '5', 'aliases': []},
    'Gart': {'ncbi_gene_id': '14738', 'ensembl_id': 'ENSMUSG00000029486', 'gene_name': 'phosphoribosylglycinamide formyltransferase', 'gene_type': 'protein_coding', 'chromosome': '5', 'aliases': []},
    'Bsg': {'ncbi_gene_id': '12215', 'ensembl_id': 'ENSMUSG00000023175', 'gene_name': 'basigin (CD147)', 'gene_type': 'protein_coding', 'chromosome': '5', 'aliases': []},
    'Pgam5': {'ncbi_gene_id': '72599', 'ensembl_id': 'ENSMUSG00000043073', 'gene_name': 'PGAM family member 5, mitochondrial serine/threonine protein phosphatase', 'gene_type': 'protein_coding', 'chromosome': '7', 'aliases': []},
    'Fgf21': {'ncbi_gene_id': '56636', 'ensembl_id': 'ENSMUSG00000030829', 'gene_name': 'fibroblast growth factor 21', 'gene_type': 'protein_coding', 'chromosome': '17', 'aliases': []},
    'Fgfr1': {'ncbi_gene_id': '14182', 'ensembl_id': 'ENSMUSG00000031565', 'gene_name': 'fibroblast growth factor receptor 1', 'gene_type': 'protein_coding', 'chromosome': '12', 'aliases': []},
    'Klb': {'ncbi_gene_id': '83379', 'ensembl_id': 'ENSMUSG00000029126', 'gene_name': 'klotho beta', 'gene_type': 'protein_coding', 'chromosome': '6', 'aliases': []},
    'Kat2a': {'ncbi_gene_id': '14534', 'ensembl_id': 'ENSMUSG00000031468', 'gene_name': 'K(lysine) acetyltransferase 2A (GCN5)', 'gene_type': 'protein_coding', 'chromosome': '9', 'aliases': []},
    'Ep300': {'ncbi_gene_id': '328572', 'ensembl_id': 'ENSMUSG00000024494', 'gene_name': 'E1A binding protein p300', 'gene_type': 'protein_coding', 'chromosome': '1', 'aliases': []},
    'Crebbp': {'ncbi_gene_id': '12914', 'ensembl_id': 'ENSMUSG00000022521', 'gene_name': 'CREB binding protein (CBP)', 'gene_type': 'protein_coding', 'chromosome': '14', 'aliases': []},
    'Nampt': {'ncbi_gene_id': '59027', 'ensembl_id': 'ENSMUSG00000025366', 'gene_name': 'nicotinamide phosphoribosyltransferase', 'gene_type': 'protein_coding', 'chromosome': '14', 'aliases': []},
    'Nmnat1': {'ncbi_gene_id': '17716', 'ensembl_id': 'ENSMUSG00000028920', 'gene_name': 'nicotinamide nucleotide adenylyltransferase 1', 'gene_type': 'protein_coding', 'chromosome': '13', 'aliases': []},
    'Nmnat2': {'ncbi_gene_id': '103858', 'ensembl_id': 'ENSMUSG00000031307', 'gene_name': 'nicotinamide nucleotide adenylyltransferase 2', 'gene_type': 'protein_coding', 'chromosome': '9', 'aliases': []},
    'Uqcc2': {'ncbi_gene_id': '67991', 'ensembl_id': 'ENSMUSG00000021946', 'gene_name': 'ubiquinol-cytochrome c reductase complex assembly factor 2', 'gene_type': 'protein_coding', 'chromosome': '12', 'aliases': []},
    'Cyc1': {'ncbi_gene_id': '66269', 'ensembl_id': 'ENSMUSG00000024129', 'gene_name': 'cytochrome c-1', 'gene_type': 'protein_coding', 'chromosome': '4', 'aliases': []},
    'Uqcrfs1': {'ncbi_gene_id': '66694', 'ensembl_id': 'ENSMUSG00000044140', 'gene_name': 'ubiquinol-cytochrome c reductase, Rieske iron-sulfur polypeptide 1', 'gene_type': 'protein_coding', 'chromosome': '7', 'aliases': []},
    'Car2': {'ncbi_gene_id': '12349', 'ensembl_id': 'ENSMUSG00000027562', 'gene_name': 'carbonic anhydrase 2', 'gene_type': 'protein_coding', 'chromosome': '9', 'aliases': []},
    'Car9': {'ncbi_gene_id': '12357', 'ensembl_id': 'ENSMUSG00000031924', 'gene_name': 'carbonic anhydrase 9', 'gene_type': 'protein_coding', 'chromosome': '8', 'aliases': []},
    'Pkp2': {'ncbi_gene_id': '21656', 'ensembl_id': 'ENSMUSG00000025400', 'gene_name': 'plakophilin 2', 'gene_type': 'protein_coding', 'chromosome': '10', 'aliases': []},
    'Ogdhl': {'ncbi_gene_id': '239017', 'ensembl_id': 'ENSMUSG00000021913', 'gene_name': 'oxoglutarate dehydrogenase-like', 'gene_type': 'protein_coding', 'chromosome': '14', 'aliases': []},


}


# ══════════════════════════════════════════════════════════════════════════════
# Lookup helpers
# ══════════════════════════════════════════════════════════════════════════════

def resolve_gene_symbol(symbol: str) -> Optional[str]:
    """Resolve a gene symbol to its canonical MGI form.

    Returns the canonical symbol, or None if the symbol represents
    a non-gene entry (metabolite, cofactor, etc.).

    If the symbol is not in GENE_ALIAS_LOOKUP and is not a known
    non-gene, returns the symbol unchanged (assumed canonical).
    """
    if symbol in GENE_ALIAS_LOOKUP:
        return GENE_ALIAS_LOOKUP[symbol]
    return symbol  # assume canonical if not in lookup


def get_ncbi_gene_id(symbol: str) -> Optional[str]:
    """Return NCBI Gene ID for a gene symbol (after alias resolution)."""
    canonical = resolve_gene_symbol(symbol)
    if canonical is None:
        return None
    info = GENE_ANNOTATIONS.get(canonical, {})
    return info.get('ncbi_gene_id')


def get_ensembl_gene_id(symbol: str) -> Optional[str]:
    """Return Ensembl Gene ID for a gene symbol (after alias resolution)."""
    canonical = resolve_gene_symbol(symbol)
    if canonical is None:
        return None
    info = GENE_ANNOTATIONS.get(canonical, {})
    return info.get('ensembl_id')


def get_gene_name(symbol: str) -> Optional[str]:
    """Return full gene name for a gene symbol (after alias resolution)."""
    canonical = resolve_gene_symbol(symbol)
    if canonical is None:
        return None
    info = GENE_ANNOTATIONS.get(canonical, {})
    return info.get('gene_name')


def get_gene_info(symbol: str) -> Dict[str, Optional[str]]:
    """Return full annotation dict for a gene, or empty dict."""
    canonical = resolve_gene_symbol(symbol)
    if canonical is None:
        return {'ncbi_gene_id': None, 'note': 'non-gene entry'}
    return GENE_ANNOTATIONS.get(canonical, {})
