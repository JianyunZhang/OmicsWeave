#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Metabolic Reaction Pathway Mapping — KEGG / Reactome / EC Number Annotations

Central lookup for all MFA flux reaction IDs used in this project,
mapping each custom reaction ID to standard public database identifiers.

This module is the single source of truth for reaction pathway annotations
and is shared across the MFA pipeline, multi-omics integration, and any
other analysis tool that needs to map flux reactions to pathways.

Each entry provides:
    - kegg_reaction : KEGG Reaction ID (e.g. 'R00771')
    - kegg_pathway  : KEGG Pathway map ID (e.g. 'map00010')
    - reactome_id   : Reactome stable ID (e.g. 'R-MMU-70349')
    - ec_number     : EC number(s) for the catalysed step
    - description   : Short English description

Mouse (Mus musculus) Reactome prefix = 'R-MMU-'.
'aggregate' reactions (e.g. GLYCOLYSIS, BIOMASS) lack single EC / KEGG reaction IDs.

References:
    - KEGG Reaction:     https://www.genome.jp/dbget-bin/www_bfind?reaction
    - KEGG Pathway:      https://www.genome.jp/kegg/pathway.html
    - Reactome:          https://reactome.org/content/schema/Reaction
    - BRENDA EC:         https://www.brenda-enzymes.org/

Usage:
    from omicsweave.core.idmap.reaction_pathway_mapping import (
        REACTION_REGISTRY,
        get_reaction_info,
        get_kegg_pathway,
        get_ec_number,
        annotate_flux_dataframe,
    )

    # Look up a single reaction
    info = get_reaction_info('HEX_c')
    # -> {'kegg_reaction': 'R00771', 'kegg_pathway': 'map00010', ...}

    # Annotate a MFA differential flux DataFrame
    df = annotate_flux_dataframe(df, reaction_col='reaction')
    # DataFrame now includes ec_number, kegg_pathway, reactome_id columns
"""

from typing import Any, Dict, List, Optional

import pandas as pd

# ─── type alias ──────────────────────────────────────────────────────────────
ReactionInfo = Dict[str, Optional[str]]

# ══════════════════════════════════════════════════════════════════════════════
# GLYCOLYSIS / GLUCONEOGENESIS
# KEGG map00010   Reactome R-MMU-70171
# ══════════════════════════════════════════════════════════════════════════════
REACTION_REGISTRY: Dict[str, ReactionInfo] = {

    # ── Glycolysis ────────────────────────────────────────────────────────────
    'GLC_exchange': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00010',
        'reactome_id': 'R-MMU-428969',
        'ec_number': None,
        'description': 'Glucose transport (GLUT/SGLT family)',
    },
    'GLC_infuse': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00010',
        'reactome_id': 'R-MMU-428969',
        'ec_number': None,
        'description': 'Glucose infusion input (mapped to glucose transport/exchange context)',
    },
    'GLYCOLYSIS': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00010',
        'reactome_id': 'R-MMU-70171',
        'ec_number': None,
        'description': 'Glycolysis (aggregate pathway)',
    },
    'GLYCOLYSIS_c': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00010',
        'reactome_id': 'R-MMU-70171',
        'ec_number': None,
        'description': 'Glycolysis (aggregate pathway, cytosolic naming alias)',
    },
    'HEX_c': {
        'kegg_reaction': 'R00771',
        'kegg_pathway': 'map00010',
        'reactome_id': 'R-MMU-70493',
        'ec_number': '2.7.1.1',
        'description': 'Hexokinase / Glucokinase: Glc + ATP -> G6P + ADP',
    },
    'PGI_c': {
        'kegg_reaction': 'R02739',
        'kegg_pathway': 'map00010',
        'reactome_id': 'R-MMU-70507',
        'ec_number': '5.3.1.9',
        'description': 'Phosphoglucose isomerase: G6P <=> F6P',
    },
    'PFK_c': {
        'kegg_reaction': 'R04779',
        'kegg_pathway': 'map00010',
        'reactome_id': 'R-MMU-70524',
        'ec_number': '2.7.1.11',
        'description': 'Phosphofructokinase-1: F6P + ATP -> F1,6BP + ADP',
    },
    'FBA_c': {
        'kegg_reaction': 'R01070',
        'kegg_pathway': 'map00010',
        'reactome_id': 'R-MMU-70550',
        'ec_number': '4.1.2.13',
        'description': 'Fructose-bisphosphate aldolase: F1,6BP <=> DHAP + GAP',
    },
    'TPI_c': {
        'kegg_reaction': 'R01015',
        'kegg_pathway': 'map00010',
        'reactome_id': 'R-MMU-70561',
        'ec_number': '5.3.1.1',
        'description': 'Triosephosphate isomerase: DHAP <=> GAP',
    },
    'TRIOSE_c': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00010',
        'reactome_id': 'R-MMU-70561',
        'ec_number': '5.3.1.1',
        'description': 'Triose phosphate metabolism (TPI + downstream)',
    },
    'GAPDH_c': {
        'kegg_reaction': 'R01061',
        'kegg_pathway': 'map00010',
        'reactome_id': 'R-MMU-70586',
        'ec_number': '1.2.1.12',
        'description': 'GAPDH: GAP + Pi + NAD+ <=> 1,3BPG + NADH',
    },
    'GAPD_c': {
        'kegg_reaction': 'R01061',
        'kegg_pathway': 'map00010',
        'reactome_id': 'R-MMU-70586',
        'ec_number': '1.2.1.12',
        'description': 'GAPDH (alias): GAP + Pi + NAD+ <=> 1,3BPG + NADH',
    },
    'PGK_c': {
        'kegg_reaction': 'R01830',
        'kegg_pathway': 'map00010',
        'reactome_id': 'R-MMU-70594',
        'ec_number': '2.7.2.3',
        'description': 'Phosphoglycerate kinase: 1,3BPG + ADP <=> 3PG + ATP',
    },
    'PGM_c': {
        'kegg_reaction': 'R01518',
        'kegg_pathway': 'map00010',
        'reactome_id': 'R-MMU-70599',
        'ec_number': '5.4.2.12',
        'description': 'Phosphoglycerate mutase: 3PG <=> 2PG',
    },
    'ENO_c': {
        'kegg_reaction': 'R00558',
        'kegg_pathway': 'map00010',
        'reactome_id': 'R-MMU-70608',
        'ec_number': '4.2.1.11',
        'description': 'Enolase: 2PG <=> PEP + H2O',
    },
    'PK_c': {
        'kegg_reaction': 'R00200',
        'kegg_pathway': 'map00010',
        'reactome_id': 'R-MMU-70613',
        'ec_number': '2.7.1.40',
        'description': 'Pyruvate kinase: PEP + ADP -> Pyr + ATP',
    },
    'PYK_c': {
        'kegg_reaction': 'R00200',
        'kegg_pathway': 'map00010',
        'reactome_id': 'R-MMU-70613',
        'ec_number': '2.7.1.40',
        'description': 'Pyruvate kinase (alias of PK_c)',
    },
    'LDH_c': {
        'kegg_reaction': 'R00703',
        'kegg_pathway': 'map00010',
        'reactome_id': 'R-MMU-71373',
        'ec_number': '1.1.1.27',
        'description': 'Lactate dehydrogenase: Pyr + NADH <=> Lac + NAD+',
    },
    'LDH_e': {
        'kegg_reaction': 'R00703',
        'kegg_pathway': 'map00010',
        'reactome_id': 'R-MMU-71373',
        'ec_number': '1.1.1.27',
        'description': 'Lactate dehydrogenase (extracellular naming alias)',
    },
    'LAC_exchange': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00010',
        'reactome_id': 'R-MMU-428969',
        'ec_number': None,
        'description': 'Lactate transport (MCT family)',
    },
    'LAC_trans': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00010',
        'reactome_id': 'R-MMU-428969',
        'ec_number': None,
        'description': 'Lactate transmembrane transport',
    },

    # ── Gluconeogenesis ───────────────────────────────────────────────────────
    'FBPase_c': {
        'kegg_reaction': 'R00752',
        'kegg_pathway': 'map00010',
        'reactome_id': 'R-MMU-70449',
        'ec_number': '3.1.3.11',
        'description': 'Fructose-1,6-bisphosphatase: F1,6BP + H2O -> F6P + Pi',
    },
    'PEPCK_c': {
        'kegg_reaction': 'R00341',
        'kegg_pathway': 'map00010',
        'reactome_id': 'R-MMU-72342',
        'ec_number': '4.1.1.32',
        'description': 'PEPCK: OAA + GTP <=> PEP + CO2 + GDP',
    },
    'PC_m': {
        'kegg_reaction': 'R00345',
        'kegg_pathway': 'map00010',
        'reactome_id': 'R-MMU-72338',
        'ec_number': '6.4.1.1',
        'description': 'Pyruvate carboxylase: Pyr + CO2 + ATP -> OAA + ADP + Pi',
    },
    'G6Pase_c': {
        'kegg_reaction': 'R00299',
        'kegg_pathway': 'map00010',
        'reactome_id': 'R-MMU-70353',
        'ec_number': '3.1.3.9',
        'description': 'Glucose-6-phosphatase: G6P + H2O -> Glc + Pi',
    },
    'G6Pase': {
        'kegg_reaction': 'R00299',
        'kegg_pathway': 'map00010',
        'reactome_id': 'R-MMU-70353',
        'ec_number': '3.1.3.9',
        'description': 'Glucose-6-phosphatase (alias)',
    },

    # ── Glycogen Metabolism ───────────────────────────────────────────────────
    # KEGG map00500   Reactome R-MMU-3322077
    'PGMT_c': {
        'kegg_reaction': 'R02367',
        'kegg_pathway': 'map00500',
        'reactome_id': 'R-MMU-3322077',
        'ec_number': '5.4.2.2',
        'description': 'Phosphoglucomutase: G1P <=> G6P',
    },
    'GLYCOGENOLYSIS': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00500',
        'reactome_id': 'R-MMU-3322077',
        'ec_number': '2.4.1.1',
        'description': 'Glycogenolysis (aggregate)',
    },
    'GLYCOGENOLYSIS_c': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00500',
        'reactome_id': 'R-MMU-3322077',
        'ec_number': '2.4.1.1',
        'description': 'Glycogenolysis (cytosolic)',
    },

    # ══════════════════════════════════════════════════════════════════════════
    # PENTOSE PHOSPHATE PATHWAY
    # KEGG map00030   Reactome R-MMU-71336
    # ══════════════════════════════════════════════════════════════════════════
    'G6PDH_c': {
        'kegg_reaction': 'R00835',
        'kegg_pathway': 'map00030',
        'reactome_id': 'R-MMU-71336',
        'ec_number': '1.1.1.49',
        'description': 'G6PDH: G6P + NADP+ -> 6PGL + NADPH',
    },
    'G6PDH2R_PGL_GND_c': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00030',
        'reactome_id': 'R-MMU-71336',
        'ec_number': '1.1.1.49',
        'description': 'Oxidative PPP aggregate (G6PDH + PGL + 6PGD)',
    },
    'PGL_c': {
        'kegg_reaction': 'R02035',
        'kegg_pathway': 'map00030',
        'reactome_id': 'R-MMU-71339',
        'ec_number': '3.1.1.31',
        'description': '6-Phosphogluconolactonase: 6PGL + H2O -> 6PG',
    },
    'GND_c': {
        'kegg_reaction': 'R01528',
        'kegg_pathway': 'map00030',
        'reactome_id': 'R-MMU-71342',
        'ec_number': '1.1.1.44',
        'description': '6-Phosphogluconate dehydrogenase: 6PG + NADP+ -> Ru5P + CO2 + NADPH',
    },
    'TKT1_c': {
        'kegg_reaction': 'R01641',
        'kegg_pathway': 'map00030',
        'reactome_id': 'R-MMU-71350',
        'ec_number': '2.2.1.1',
        'description': 'Transketolase 1: R5P + Xu5P <=> S7P + GAP',
    },
    'TKT2_c': {
        'kegg_reaction': 'R01827',
        'kegg_pathway': 'map00030',
        'reactome_id': 'R-MMU-71355',
        'ec_number': '2.2.1.1',
        'description': 'Transketolase 2: E4P + Xu5P <=> F6P + GAP',
    },
    'TALA_c': {
        'kegg_reaction': 'R01826',
        'kegg_pathway': 'map00030',
        'reactome_id': 'R-MMU-71352',
        'ec_number': '2.2.1.2',
        'description': 'Transaldolase: S7P + GAP <=> E4P + F6P',
    },
    'RPI_c': {
        'kegg_reaction': 'R01056',
        'kegg_pathway': 'map00030',
        'reactome_id': 'R-MMU-71346',
        'ec_number': '5.3.1.6',
        'description': 'Ribose-5-phosphate isomerase: R5P <=> Ru5P',
    },
    'RPE_c': {
        'kegg_reaction': 'R01529',
        'kegg_pathway': 'map00030',
        'reactome_id': 'R-MMU-71348',
        'ec_number': '5.1.3.1',
        'description': 'Ribulose-5-phosphate epimerase: Ru5P <=> Xu5P',
    },

    # ══════════════════════════════════════════════════════════════════════════
    # TCA CYCLE
    # KEGG map00020   Reactome R-MMU-71406
    # ══════════════════════════════════════════════════════════════════════════
    'CS_m': {
        'kegg_reaction': 'R00351',
        'kegg_pathway': 'map00020',
        'reactome_id': 'R-MMU-71403',
        'ec_number': '2.3.3.1',
        'description': 'Citrate synthase: OAA + AcCoA + H2O -> Cit + CoA',
    },
    'ACONT_m': {
        'kegg_reaction': 'R01325',
        'kegg_pathway': 'map00020',
        'reactome_id': 'R-MMU-71406',
        'ec_number': '4.2.1.3',
        'description': 'Aconitase: Cit <=> Isocit',
    },
    'ICDH_m': {
        'kegg_reaction': 'R00709',  # NAD+ dependent
        'kegg_pathway': 'map00020',
        'reactome_id': 'R-MMU-71412',
        'ec_number': '1.1.1.41',
        'description': 'Isocitrate dehydrogenase (NAD+): Isocit + NAD+ => aKG + CO2 + NADH',
    },
    'AKGDH_m': {
        'kegg_reaction': 'R00267',  # same as AKGD_m
        'kegg_pathway': 'map00020',
        'reactome_id': 'R-MMU-71419',
        'ec_number': '1.2.4.2',
        'description': 'aKG dehydrogenase: aKG + CoA + NAD+ -> SuccCoA + CO2 + NADH',
    },
    'AKGD_m': {
        'kegg_reaction': 'R00267',
        'kegg_pathway': 'map00020',
        'reactome_id': 'R-MMU-71419',
        'ec_number': '1.2.4.2',
        'description': 'aKG dehydrogenase (alias)',
    },
    'SUCOAS_m': {
        'kegg_reaction': 'R00405',  # reverse direction
        'kegg_pathway': 'map00020',
        'reactome_id': 'R-MMU-71424',
        'ec_number': '6.2.1.5',
        'description': 'Succinyl-CoA synthetase: SuccCoA + ADP + Pi <=> Succ + ATP + CoA',
    },
    'SUCD_m': {
        'kegg_reaction': 'R02164',
        'kegg_pathway': 'map00020',
        'reactome_id': 'R-MMU-71428',
        'ec_number': '1.3.5.1',
        'description': 'Succinate dehydrogenase (Complex II): Succ + Q <=> Fum + QH2',
    },
    'SUC_trans': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00020',
        'reactome_id': 'R-MMU-372843',
        'ec_number': None,
        'description': 'Succinate transport (SLC25A10 dicarboxylate carrier)',
    },
    'FRD_m': {
        'kegg_reaction': 'R02164',  # reverse of SUCD
        'kegg_pathway': 'map00020',
        'reactome_id': 'R-MMU-71428',
        'ec_number': '1.3.5.1',
        'description': 'Fumarate reductase (reverse SDH)',
    },
    'FUM_m': {
        'kegg_reaction': 'R01082',
        'kegg_pathway': 'map00020',
        'reactome_id': 'R-MMU-71433',
        'ec_number': '4.2.1.2',
        'description': 'Fumarase: Fum + H2O <=> Mal',
    },
    'FUM1_m': {
        'kegg_reaction': 'R01082',
        'kegg_pathway': 'map00020',
        'reactome_id': 'R-MMU-71433',
        'ec_number': '4.2.1.2',
        'description': 'Fumarase (alias)',
    },
    'FUMH_m': {
        'kegg_reaction': 'R01082',
        'kegg_pathway': 'map00020',
        'reactome_id': 'R-MMU-71433',
        'ec_number': '4.2.1.2',
        'description': 'Fumarase (alias)',
    },
    'FUM_trans': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00020',
        'reactome_id': None,
        'ec_number': None,
        'description': 'Fumarate transport',
    },
    'MDH_m': {
        'kegg_reaction': 'R00342',
        'kegg_pathway': 'map00020',
        'reactome_id': 'R-MMU-71437',
        'ec_number': '1.1.1.37',
        'description': 'Malate dehydrogenase (mitochondrial): Mal + NAD+ <=> OAA + NADH',
    },
    'MDH2_m': {
        'kegg_reaction': 'R00342',
        'kegg_pathway': 'map00020',
        'reactome_id': 'R-MMU-71437',
        'ec_number': '1.1.1.37',
        'description': 'Malate dehydrogenase (mitochondrial, explicit)',
    },
    'MDH_c': {
        'kegg_reaction': 'R00342',
        'kegg_pathway': 'map00020',
        'reactome_id': 'R-MMU-71437',
        'ec_number': '1.1.1.37',
        'description': 'Malate dehydrogenase (cytosolic)',
    },
    'MDH1_c': {
        'kegg_reaction': 'R00342',
        'kegg_pathway': 'map00020',
        'reactome_id': 'R-MMU-71437',
        'ec_number': '1.1.1.37',
        'description': 'Malate dehydrogenase (cytosolic, explicit)',
    },
    'CIT_m': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00020',
        'reactome_id': 'R-MMU-75849',
        'ec_number': None,
        'description': 'Citrate mitochondrial pool / transport (SLC25A1 context)',
    },
    'CIT_trans': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00020',
        'reactome_id': 'R-MMU-75849',
        'ec_number': None,
        'description': 'Citrate transport (SLC25A1 exchanges CIT and MAL across mitochondrial membrane)',
    },
    'MAL_trans': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00020',
        'reactome_id': 'R-MMU-198440',
        'ec_number': None,
        'description': 'Malate transport (SLC25A11 oxoglutarate/malate carrier)',
    },

    # ══════════════════════════════════════════════════════════════════════════
    # PYRUVATE DEHYDROGENASE & PYRUVATE METABOLISM
    # KEGG map00020 / map00620   Reactome R-MMU-71406
    # ══════════════════════════════════════════════════════════════════════════
    'PDH_m': {
        'kegg_reaction': 'R00209',
        'kegg_pathway': 'map00020',
        'reactome_id': 'R-MMU-71399',
        'ec_number': '1.2.4.1',
        'description': 'Pyruvate dehydrogenase: Pyr + CoA + NAD+ -> AcCoA + CO2 + NADH',
    },
    'PYR_exchange': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00620',
        'reactome_id': 'R-MMU-428969',
        'ec_number': None,
        'description': 'Pyruvate transport (MPC complex)',
    },
    'PYR_trans': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00620',
        'reactome_id': 'R-MMU-428969',
        'ec_number': None,
        'description': 'Pyruvate transmembrane transport',
    },

    # ══════════════════════════════════════════════════════════════════════════
    # FATTY ACID METABOLISM
    # KEGG map00071 (FA degradation) / map00061 (FA biosynthesis)
    # Reactome R-MMU-77288 / R-MMU-75105
    # ══════════════════════════════════════════════════════════════════════════
    'FA_exchange': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00071',
        'reactome_id': 'R-MMU-5659128',
        'ec_number': None,
        'description': 'Fatty acid uptake / transport',
    },
    'FA_usage': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00071',
        'reactome_id': 'R-MMU-77288',
        'ec_number': None,
        'description': 'Fatty acid utilization (aggregate)',
    },
    'FAO_m': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00071',
        'reactome_id': 'R-MMU-77288',
        'ec_number': None,
        'description': 'Fatty acid beta-oxidation (mitochondrial, aggregate)',
    },
    'CPT1': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00071',
        'reactome_id': 'R-MMU-77289',
        'ec_number': '2.3.1.21',
        'description': 'Carnitine palmitoyltransferase 1 (rate-limiting FAO step)',
    },
    'CARNITINE_shuttle': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00071',
        'reactome_id': 'R-MMU-77289',
        'ec_number': None,
        'description': 'Carnitine shuttle (CPT1 + CACT + CPT2)',
    },
    'ACC_c': {
        'kegg_reaction': 'R00742',
        'kegg_pathway': 'map00061',
        'reactome_id': 'R-MMU-75105',
        'ec_number': '6.4.1.2',
        'description': 'Acetyl-CoA carboxylase: AcCoA + CO2 + ATP -> MalonylCoA + ADP + Pi',
    },
    'FAS_c': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00061',
        'reactome_id': 'R-MMU-75109',
        'ec_number': '2.3.1.85',
        'description': 'Fatty acid synthase (aggregate: AcCoA + MalCoA -> Palmitate)',
    },
    'ACOT_c': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00061',
        'reactome_id': 'R-MMU-5696103',
        'ec_number': '3.1.2.20',
        'description': 'Acyl-CoA thioesterase: Acyl-CoA + H2O -> FFA + CoA',
    },

    # ══════════════════════════════════════════════════════════════════════════
    # KETONE BODY METABOLISM
    # KEGG map00072   Reactome R-MMU-77288
    # ══════════════════════════════════════════════════════════════════════════
    'BHBDH_m': {
        'kegg_reaction': 'R01361',
        'kegg_pathway': 'map00072',
        'reactome_id': 'R-MMU-77296',
        'ec_number': '1.1.1.30',
        'description': 'beta-Hydroxybutyrate dehydrogenase: AcAc + NADH <=> BHB + NAD+',
    },
    'BHB_trans': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00072',
        'reactome_id': 'R-MMU-428969',
        'ec_number': None,
        'description': 'BHB transmembrane transport (MCT)',
    },
    'BHB_exchange': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00072',
        'reactome_id': 'R-MMU-428969',
        'ec_number': None,
        'description': 'BHB exchange',
    },
    'BHB_input': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00072',
        'reactome_id': 'R-MMU-428969',
        'ec_number': None,
        'description': 'BHB uptake',
    },
    'BHB_output': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00072',
        'reactome_id': 'R-MMU-428969',
        'ec_number': None,
        'description': 'BHB release',
    },
    'ACAC_trans': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00072',
        'reactome_id': 'R-MMU-428969',
        'ec_number': None,
        'description': 'Acetoacetate transport',
    },
    'ACAC_production': {
        'kegg_reaction': 'R01360',
        'kegg_pathway': 'map00072',
        'reactome_id': 'R-MMU-77295',
        'ec_number': '4.1.3.4',
        'description': 'HMG-CoA lyase: HMG-CoA <=> AcCoA + Acetoacetate',
    },
    'ACAC_production_m': {
        'kegg_reaction': 'R01360',
        'kegg_pathway': 'map00072',
        'reactome_id': 'R-MMU-77295',
        'ec_number': '4.1.3.4',
        'description': 'HMG-CoA lyase mitochondrial naming alias',
    },
    'ACACCOA_lyase_m': {
        'kegg_reaction': 'R01360',
        'kegg_pathway': 'map00072',
        'reactome_id': 'R-MMU-77295',
        'ec_number': '4.1.3.4',
        'description': 'Acetoacetyl-CoA lyase / HMG-CoA lyase naming alias in ketogenesis branch',
    },
    'SCOT_m': {
        'kegg_reaction': 'R01780',
        'kegg_pathway': 'map00072',
        'reactome_id': 'R-MMU-77310',
        'ec_number': '2.8.3.5',
        'description': 'SCOT: SuccCoA + 3-oxo acid <=> Succ + 3-oxoacyl-CoA',
    },

    # ══════════════════════════════════════════════════════════════════════════
    # AMINO ACID METABOLISM
    # ══════════════════════════════════════════════════════════════════════════

    # ── Glutamate / Glutamine ────────────────────────────────────────────────
    'GLU_trans': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00251',
        'reactome_id': 'R-MMU-428562',
        'ec_number': None,
        'description': 'Glutamate transport (EAAT family)',
    },
    'GLU_exchange': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00251',
        'reactome_id': 'R-MMU-428562',
        'ec_number': None,
        'description': 'Glutamate exchange',
    },
    'GLUD_m': {
        'kegg_reaction': 'R00243',
        'kegg_pathway': 'map00251',
        'reactome_id': 'R-MMU-70593',
        'ec_number': '1.4.1.3',
        'description': 'Glutamate dehydrogenase: Glu + H2O + NAD+ <=> aKG + NH3 + NADH',
    },
    'GLNS_c': {
        'kegg_reaction': 'R00253',
        'kegg_pathway': 'map00251',
        'reactome_id': 'R-MMU-70602',
        'ec_number': '6.3.1.2',
        'description': 'Glutamine synthetase: Glu + NH3 + ATP -> Gln + ADP + Pi',
    },
    'GLND_m': {
        'kegg_reaction': 'R00256',
        'kegg_pathway': 'map00251',
        'reactome_id': 'R-MMU-70603',
        'ec_number': '3.5.1.2',
        'description': 'Glutaminase: Gln + H2O -> Glu + NH3',
    },
    'GLN_exchange': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00251',
        'reactome_id': 'R-MMU-428562',
        'ec_number': None,
        'description': 'Glutamine exchange',
    },
    'GLN_trans': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00251',
        'reactome_id': 'R-MMU-428562',
        'ec_number': None,
        'description': 'Glutamine transport',
    },

    # ── Aspartate / Alanine Transaminases ────────────────────────────────────
    'ASPTA_m': {
        'kegg_reaction': 'R00397',
        'kegg_pathway': 'map00250',
        'reactome_id': 'R-MMU-70508',
        'ec_number': '2.6.1.1',
        'description': 'Aspartate transaminase (mito): Glu + OAA <=> aKG + Asp',
    },
    'ASPTA_c': {
        'kegg_reaction': 'R00397',
        'kegg_pathway': 'map00250',
        'reactome_id': 'R-MMU-70508',
        'ec_number': '2.6.1.1',
        'description': 'Aspartate transaminase (cyto): Glu + OAA <=> aKG + Asp',
    },
    'ASPGLU_m': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00250',
        'reactome_id': 'R-MMU-70508',
        'ec_number': '2.6.1.1',
        'description': 'Aspartate-glutamate metabolism (malate-aspartate shuttle)',
    },
    'ASPGLu_m': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00250',
        'reactome_id': 'R-MMU-70508',
        'ec_number': '2.6.1.1',
        'description': 'Aspartate-glutamate metabolism (case variant)',
    },
    'ASP_exchange': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00250',
        'reactome_id': 'R-MMU-428562',
        'ec_number': None,
        'description': 'Aspartate exchange',
    },
    'GPT_c': {
        'kegg_reaction': 'R00258',
        'kegg_pathway': 'map00250',
        'reactome_id': 'R-MMU-70512',
        'ec_number': '2.6.1.2',
        'description': 'Alanine transaminase: Glu + Pyr <=> aKG + Ala',
    },
    'ALA_exchange': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00250',
        'reactome_id': 'R-MMU-428562',
        'ec_number': None,
        'description': 'Alanine exchange',
    },
    'AKGMAL_m': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00250',
        'reactome_id': 'R-MMU-198440',
        'ec_number': None,
        'description': 'aKG-malate interconversion (malate-aspartate shuttle via SLC25A11)',
    },
    'AKG_trans': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00250',
        'reactome_id': 'R-MMU-198440',
        'ec_number': None,
        'description': 'Alpha-ketoglutarate transport (SLC25A11 oxoglutarate/malate carrier)',
    },

    # ── Other Amino Acids ────────────────────────────────────────────────────
    'SER_exchange': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00260',
        'reactome_id': 'R-MMU-428562',
        'ec_number': None,
        'description': 'Serine exchange',
    },
    'GLY_exchange': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00260',
        'reactome_id': 'R-MMU-428562',
        'ec_number': None,
        'description': 'Glycine exchange',
    },
    'BCAA_trans': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00280',
        'reactome_id': 'R-MMU-428562',
        'ec_number': None,
        'description': 'Branched-chain amino acid transport',
    },
    'HIS_metabolism': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00340',
        'reactome_id': 'R-MMU-71064',
        'ec_number': None,
        'description': 'Histidine metabolism',
    },
    'PRO_metabolism': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00330',
        'reactome_id': 'R-MMU-70688',
        'ec_number': None,
        'description': 'Proline metabolism (catabolism via PRODH/ALDH4A1)',
    },
    'TRP_metabolism': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00380',
        'reactome_id': 'R-MMU-71240',
        'ec_number': None,
        'description': 'Tryptophan metabolism (kynurenine pathway)',
    },
    'MET_cycle': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00270',
        'reactome_id': 'R-MMU-1614635',
        'ec_number': None,
        'description': 'Methionine cycle / sulfur amino acid metabolism (MAT/CBS/CTH)',
    },
    'PHE_TYR_metabolism': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00360',
        'reactome_id': 'R-MMU-8963691',
        'ec_number': None,
        'description': 'Phenylalanine and tyrosine metabolism',
    },
    'BALA_metabolism': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00410',
        'reactome_id': 'R-MMU-71291',
        'ec_number': None,
        'description': 'Beta-alanine metabolism (broad amino acid metabolism context)',
    },
    'ARG_metabolism': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00330',
        'reactome_id': 'R-MMU-70635',
        'ec_number': None,
        'description': 'Arginine metabolism (urea cycle context: ASS/ARG1/ARG2)',
    },
    'ORN_metabolism': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00330',
        'reactome_id': 'R-MMU-70635',
        'ec_number': None,
        'description': 'Ornithine metabolism (urea cycle context: ODC/OTC)',
    },
    'TAU_metabolism': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00430',
        'reactome_id': 'R-MMU-1614635',
        'ec_number': None,
        'description': 'Taurine metabolism (cysteine -> CSA -> hypotaurine -> taurine)',
    },
    'HYP_metabolism': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00330',
        'reactome_id': 'R-MMU-70688',
        'ec_number': None,
        'description': 'Hydroxyproline metabolism (proline catabolism context, ALDH4A1)',
    },
    'INO_metabolism': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00220',
        'reactome_id': 'R-MMU-15869',
        'ec_number': None,
        'description': 'Inosine / purine metabolism (nucleotide metabolism)',
    },

    # ── Urea Cycle ───────────────────────────────────────────────────────────
    'CPS1_m': {
        'kegg_reaction': 'R00149',
        'kegg_pathway': 'map00330',
        'reactome_id': 'R-MMU-70589',
        'ec_number': '6.3.4.16',
        'description': 'Carbamoyl phosphate synthetase I: NH3 + HCO3- + 2ATP -> CP + 2ADP + Pi',
    },
    'UREA_cycle': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00330',
        'reactome_id': 'R-MMU-71064',
        'ec_number': None,
        'description': 'Urea cycle (aggregate)',
    },

    # ══════════════════════════════════════════════════════════════════════════
    # GLYCEROL METABOLISM
    # KEGG map00561   Reactome R-MMU-70363
    # ══════════════════════════════════════════════════════════════════════════
    'GPS_c': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00561',
        'reactome_id': 'R-MMU-70363',
        'ec_number': None,
        'description': 'Glycerol-3-phosphate shuttle (aggregate)',
    },
    'GLYCEROL_exchange': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00561',
        'reactome_id': 'R-MMU-428969',
        'ec_number': None,
        'description': 'Glycerol transport (AQP family)',
    },
    'GLYCK_c': {
        'kegg_reaction': 'R00847',
        'kegg_pathway': 'map00561',
        'reactome_id': 'R-MMU-70365',
        'ec_number': '2.7.1.30',
        'description': 'Glycerol kinase: Glycerol + ATP -> G3P + ADP',
    },
    'GK_c': {
        'kegg_reaction': 'R00847',
        'kegg_pathway': 'map00561',
        'reactome_id': 'R-MMU-70365',
        'ec_number': '2.7.1.30',
        'description': 'Glycerol kinase (alias)',
    },
    'G3PP_c': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00561',
        'reactome_id': None,
        'ec_number': '3.1.3.-',
        'description': 'Glycerol-3-phosphate phosphatase branch (aggregate / model-specific)',
    },
    'G3PDH_c': {
        'kegg_reaction': 'R00842',
        'kegg_pathway': 'map00561',
        'reactome_id': 'R-MMU-70363',
        'ec_number': '1.1.1.8',
        'description': 'Glycerol-3-phosphate dehydrogenase (cytosolic): DHAP + NADH <=> G3P + NAD+',
    },
    'G3PD_m': {
        'kegg_reaction': 'R00849',
        'kegg_pathway': 'map00561',
        'reactome_id': 'R-MMU-70363',
        'ec_number': '1.1.5.3',
        'description': 'Glycerol-3-phosphate dehydrogenase (mitochondrial): G3P + Q <=> DHAP + QH2',
    },

    # ══════════════════════════════════════════════════════════════════════════
    # ACETATE / ACETYL-CoA METABOLISM
    # ══════════════════════════════════════════════════════════════════════════
    'ACE_exchange': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00620',
        'reactome_id': 'R-MMU-428969',
        'ec_number': None,
        'description': 'Acetate exchange',
    },
    'ACE_trans': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00620',
        'reactome_id': 'R-MMU-428969',
        'ec_number': None,
        'description': 'Acetate transport',
    },
    'ACITL_c': {
        'kegg_reaction': 'R00352',
        'kegg_pathway': 'map00020',
        'reactome_id': 'R-MMU-70403',
        'ec_number': '2.3.3.8',
        'description': 'ATP citrate lyase: Cit + ATP + CoA -> AcCoA + OAA + ADP + Pi',
    },
    'ACSS2_m': {
        'kegg_reaction': 'R00235',
        'kegg_pathway': 'map00620',
        'reactome_id': None,
        'ec_number': '6.2.1.1',
        'description': 'Acetyl-CoA synthetase (ACSS2-like mapping): Acetate + CoA + ATP -> Acetyl-CoA + AMP + PPi',
    },

    # ══════════════════════════════════════════════════════════════════════════
    # LIPID SYNTHESIS / MOBILIZATION
    # ══════════════════════════════════════════════════════════════════════════
    'LIPID_c': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00564',
        'reactome_id': 'R-MMU-75109',
        'ec_number': None,
        'description': 'Triglyceride synthesis (aggregate)',
    },
    'lipid_mobilization': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00561',
        'reactome_id': 'R-MMU-73931',
        'ec_number': None,
        'description': 'Lipid mobilization / lipolysis (aggregate)',
    },
    'FA_esterification_c': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00561',
        'reactome_id': None,
        'ec_number': None,
        'description': 'Fatty acid esterification to triacylglycerol (lumped FA_c -> TG_store; '
                       'Kennedy pathway GPAT / AGPAT / lipin / DGAT, plus MGAT)',
    },
    'lipid_mobilization_c': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00561',
        'reactome_id': 'R-MMU-73931',
        'ec_number': None,
        'description': 'Lipid mobilization / lipolysis (cytosolic naming alias)',
    },
    'PHOSPHOLIPID': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00564',
        'reactome_id': 'R-MMU-1483257',
        'ec_number': None,
        'description': 'Phospholipid metabolism (PC/PA/PS/PI/CL synthesis)',
    },
    'SPHINGOLIPID': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00600',
        'reactome_id': 'R-MMU-428157',
        'ec_number': None,
        'description': 'Sphingolipid metabolism (de novo biosynthesis + catabolism)',
    },
    'HMGCR': {
        'kegg_reaction': 'R02082',
        'kegg_pathway': 'map00100',
        'reactome_id': 'R-MMU-194692',
        'ec_number': '1.1.1.34',
        'description': 'HMG-CoA reductase (cholesterol synthesis rate-limiting)',
    },
    'BILE_acid': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00120',
        'reactome_id': 'R-MMU-193366',
        'ec_number': None,
        'description': 'Bile acid metabolism (aggregate)',
    },

    # ══════════════════════════════════════════════════════════════════════════
    # MALIC ENZYME
    # ══════════════════════════════════════════════════════════════════════════
    'ME_c': {
        'kegg_reaction': 'R00216',
        'kegg_pathway': 'map00620',
        'reactome_id': 'R-MMU-71437',
        'ec_number': '1.1.1.40',
        'description': 'Malic enzyme (NADP+): Mal + NADP+ <=> Pyr + CO2 + NADPH',
    },
    'ME_m': {
        'kegg_reaction': 'R00214',
        'kegg_pathway': 'map00620',
        'reactome_id': 'R-MMU-71437',
        'ec_number': '1.1.1.38',
        'description': 'Malic enzyme (mitochondrial, NAD+): Mal + NAD+ <=> Pyr + CO2 + NADH',
    },
    'ME1_c': {
        'kegg_reaction': 'R00216',
        'kegg_pathway': 'map00620',
        'reactome_id': 'R-MMU-71437',
        'ec_number': '1.1.1.40',
        'description': 'Malic enzyme 1 (cytosolic, NADP+)',
    },
    'ME2_c': {
        'kegg_reaction': 'R00214',
        'kegg_pathway': 'map00620',
        'reactome_id': 'R-MMU-71437',
        'ec_number': '1.1.1.38',
        'description': 'Malic enzyme 2 (mitochondrial, NAD+)',
    },
    'ME2_m': {
        'kegg_reaction': 'R00214',
        'kegg_pathway': 'map00620',
        'reactome_id': 'R-MMU-71437',
        'ec_number': '1.1.1.38',
        'description': 'Malic enzyme 2 (mitochondrial, explicit)',
    },

    # ══════════════════════════════════════════════════════════════════════════
    # BIOMASS / GROWTH
    # ══════════════════════════════════════════════════════════════════════════
    'BIOMASS_REACTION_1': {
        'kegg_reaction': None,
        'kegg_pathway': None,
        'reactome_id': None,
        'ec_number': None,
        'description': 'Biomass synthesis (pseudo-reaction)',
    },
    'BIOMASS_REACTION_2': {
        'kegg_reaction': None,
        'kegg_pathway': None,
        'reactome_id': None,
        'ec_number': None,
        'description': 'Biomass synthesis (pseudo-reaction variant 2)',
    },
    'BIOMASS_REACTION_3': {
        'kegg_reaction': None,
        'kegg_pathway': None,
        'reactome_id': None,
        'ec_number': None,
        'description': 'Biomass synthesis (pseudo-reaction variant 3)',
    },
    'BIOMASS_REACTION_4': {
        'kegg_reaction': None,
        'kegg_pathway': None,
        'reactome_id': None,
        'ec_number': None,
        'description': 'Biomass synthesis (pseudo-reaction variant 4)',
    },
    'BIOMASS_REACTION_5': {
        'kegg_reaction': None,
        'kegg_pathway': None,
        'reactome_id': None,
        'ec_number': None,
        'description': 'Biomass synthesis (pseudo-reaction variant 5)',
    },
    'BIOMASS_REACTION_6': {
        'kegg_reaction': None,
        'kegg_pathway': None,
        'reactome_id': None,
        'ec_number': None,
        'description': 'Biomass synthesis (pseudo-reaction, 6-component)',
    },
    'BIOMASS_REACTION_7': {
        'kegg_reaction': None,
        'kegg_pathway': None,
        'reactome_id': None,
        'ec_number': None,
        'description': 'Biomass synthesis (pseudo-reaction variant 7)',
    },

    # ══════════════════════════════════════════════════════════════════════════
    # MISCELLANEOUS / NUCLEOTIDE / COFACTOR METABOLISM
    # ══════════════════════════════════════════════════════════════════════════
    'NP_c': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00240',
        'reactome_id': None,
        'ec_number': '2.4.2.1',
        'description': 'Nucleoside phosphorylase (purine nucleoside + Pi <=> base + R1P)',
    },
    'Salvage_c': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00230',
        'reactome_id': 'R-MMU-74217',
        'ec_number': None,
        'description': 'Nucleotide salvage pathway (purine salvage)',
    },
    'PHGDH_PSAT_PSP_c': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00260',
        'reactome_id': None,
        'ec_number': '1.1.1.95;2.6.1.52;3.1.3.3',
        'description': 'De novo serine biosynthesis aggregate (PHGDH + PSAT + PSP)',
    },
    'PROS_m': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00330',
        'reactome_id': 'R-MMU-70688',
        'ec_number': '1.5.1.2',
        'description': 'Proline synthesis / turnover aggregate (PYCR, mitochondrial)',
    },
    'ADK_c': {
        'kegg_reaction': 'R00127',
        'kegg_pathway': 'map00230',
        'reactome_id': None,
        'ec_number': '2.7.4.3',
        'description': 'Adenylate kinase: 2ADP <=> ATP + AMP',
    },
    'CK_c': {
        'kegg_reaction': 'R01878',
        'kegg_pathway': 'map00230',
        'reactome_id': None,
        'ec_number': '2.7.3.2',
        'description': 'Creatine kinase: Cr + ATP <=> PCr + ADP',
    },
    'ALLANTOIN': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00260',
        'reactome_id': 'R-MMU-74259',
        'ec_number': None,
        'description': 'Allantoin metabolism (purine catabolism)',
    },
    'GSH_synthesis': {
        'kegg_reaction': 'R00497',
        'kegg_pathway': 'map00480',
        'reactome_id': None,
        'ec_number': '6.3.2.3',
        'description': 'Glutathione synthase (step 2): gamma-Glu-Cys + Gly + ATP -> GSH + ADP + Pi',
    },
    'HEME_synthesis': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00860',
        'reactome_id': 'R-MMU-189439',
        'ec_number': None,
        'description': 'Heme biosynthesis (aggregate)',
    },
    'POLYAMINE_synthesis': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00330',
        'reactome_id': 'R-MMU-351200',
        'ec_number': None,
        'description': 'Polyamine synthesis (putrescine -> spermidine/spermine via SRM/SMS)',
    },
    'FOLATE_cycle': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00670',
        'reactome_id': 'R-MMU-196757',
        'ec_number': None,
        'description': 'Folate / one-carbon cycle (folate and pterine metabolism)',
    },
    'PYROGLU': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00251',
        'reactome_id': 'R-MMU-174403',
        'ec_number': None,
        'description': 'Pyroglutamate (5-oxoproline) metabolism (glutathione recycling context)',
    },

    # ── G18 M0b (2026-11-19): network reactions without earlier annotation ──────
    # KEGG reaction / map / EC checked with the KEGG REST API; lumped reactions have no single reaction ID
    'ACSL_PAL_c': {
        'kegg_reaction': 'R01280',
        'kegg_pathway': 'map00071',
        'reactome_id': None,
        'ec_number': '6.2.1.3',
        'description': 'Long-chain acyl-CoA synthetase: palmitate + CoA -> palmitoyl-CoA',
    },
    'CPT1_c': {
        'kegg_reaction': 'R01923',
        'kegg_pathway': 'map00071',
        'reactome_id': None,
        'ec_number': '2.3.1.21',
        'description': 'Carnitine palmitoyltransferase 1: palmitoyl-CoA -> palmitoylcarnitine (outer membrane)',
    },
    'CACT_PALCAR': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00071',
        'reactome_id': None,
        'ec_number': None,
        'description': 'Carnitine-acylcarnitine translocase: palmitoylcarnitine into mitochondria',
    },
    'CPT2_m': {
        'kegg_reaction': 'R01923',
        'kegg_pathway': 'map00071',
        'reactome_id': None,
        'ec_number': '2.3.1.21',
        'description': 'Carnitine palmitoyltransferase 2: palmitoylcarnitine -> palmitoyl-CoA (inner membrane)',
    },
    'CRAT_m': {
        'kegg_reaction': 'R02396',
        'kegg_pathway': None,
        'reactome_id': None,
        'ec_number': '2.3.1.7',
        'description': 'Carnitine acetyltransferase: acetyl-CoA <-> acetylcarnitine',
    },
    'PAL_beta_oxidation': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00071',
        'reactome_id': None,
        'ec_number': None,
        'description': 'Mitochondrial beta-oxidation of palmitate to acetyl-CoA (lumped)',
    },
    'ELOVL6_c': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00062',
        'reactome_id': None,
        'ec_number': '2.3.1.199',
        'description': 'Fatty acid elongase 6: palmitate (C16) -> stearate (C18)',
    },
    'SCD_PAL_c': {
        'kegg_reaction': None,
        'kegg_pathway': 'map01040',
        'reactome_id': None,
        'ec_number': '1.14.19.1',
        'description': 'Stearoyl-CoA desaturase: palmitate -> palmitoleate',
    },
    'SCD_STE_c': {
        'kegg_reaction': 'R02222',
        'kegg_pathway': 'map01040',
        'reactome_id': None,
        'ec_number': '1.14.19.1',
        'description': 'Stearoyl-CoA desaturase: stearate -> oleate',
    },
    'GPAT_c': {
        'kegg_reaction': 'R00851',
        'kegg_pathway': 'map00561',
        'reactome_id': None,
        'ec_number': '2.3.1.15',
        'description': 'Glycerol-3-phosphate acyltransferase: G3P -> lysophosphatidate',
    },
    'AGPAT_c': {
        'kegg_reaction': 'R02241',
        'kegg_pathway': 'map00561',
        'reactome_id': None,
        'ec_number': '2.3.1.51',
        'description': '1-Acylglycerol-3-phosphate acyltransferase: LPA -> phosphatidate',
    },
    'LIPIN_c': {
        'kegg_reaction': 'R02239',
        'kegg_pathway': 'map00561',
        'reactome_id': None,
        'ec_number': '3.1.3.4',
        'description': 'Phosphatidate phosphatase (lipin): PA -> DAG',
    },
    'DGAT_c': {
        'kegg_reaction': 'R02251',
        'kegg_pathway': 'map00561',
        'reactome_id': None,
        'ec_number': '2.3.1.20',
        'description': 'Diacylglycerol acyltransferase: DAG -> TAG',
    },
    'ATGL_c': {
        'kegg_reaction': 'R02250',
        'kegg_pathway': 'map00561',
        'reactome_id': None,
        'ec_number': '3.1.1.3',
        'description': 'Adipose triglyceride lipase: TAG -> DAG',
    },
    'HSL_c': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00561',
        'reactome_id': None,
        'ec_number': '3.1.1.79',
        'description': 'Hormone-sensitive lipase: DAG -> MAG',
    },
    'MGL_c': {
        'kegg_reaction': 'R01351',
        'kegg_pathway': 'map00561',
        'reactome_id': None,
        'ec_number': '3.1.1.23',
        'description': 'Monoacylglycerol lipase: MAG -> glycerol',
    },
    'PAL_esterification_c': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00561',
        'reactome_id': None,
        'ec_number': None,
        'description': 'Palmitate esterification to triacylglycerol (lumped; members as FA_esterification_c)',
    },
    'PAL_TG_hydrolysis_c': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00561',
        'reactome_id': None,
        'ec_number': None,
        'description': 'Hydrolysis of palmitate-containing triacylglycerol (lumped lipolysis)',
    },
    'GLYCOGEN_SYNTHESIS_c': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00500',
        'reactome_id': None,
        'ec_number': None,
        'description': 'Glycogen synthesis from glucose-1-phosphate (UGP2, glycogen synthase, branching enzyme; lumped)',
    },
    'CPS1_c': {
        'kegg_reaction': 'R00149',
        'kegg_pathway': 'map00220',
        'reactome_id': None,
        'ec_number': '6.3.4.16',
        'description': 'Carbamoyl-phosphate synthetase 1: NH3 -> carbamoyl phosphate',
    },
    'OTC_c': {
        'kegg_reaction': 'R01398',
        'kegg_pathway': 'map00220',
        'reactome_id': None,
        'ec_number': '2.1.3.3',
        'description': 'Ornithine transcarbamylase: ornithine + carbamoyl phosphate -> citrulline',
    },
    'ASS_ASL_c': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00220',
        'reactome_id': None,
        'ec_number': '6.3.4.5, 4.3.2.1',
        'description': 'Argininosuccinate synthetase + lyase: citrulline + aspartate -> arginine (lumped)',
    },
    'ARGN_c': {
        'kegg_reaction': 'R00551',
        'kegg_pathway': 'map00220',
        'reactome_id': None,
        'ec_number': '3.5.3.1',
        'description': 'Arginase: arginine -> ornithine + urea',
    },
    'GATM_c': {
        'kegg_reaction': 'R00565',
        'kegg_pathway': 'map00330',
        'reactome_id': None,
        'ec_number': '2.1.4.1',
        'description': 'Glycine amidinotransferase: arginine + glycine -> ornithine + guanidinoacetate',
    },
    'GAMT_c': {
        'kegg_reaction': 'R01883',
        'kegg_pathway': 'map00330',
        'reactome_id': None,
        'ec_number': '2.1.1.2',
        'description': 'Guanidinoacetate N-methyltransferase: guanidinoacetate -> creatine',
    },
    'OAT_c': {
        'kegg_reaction': 'R00667',
        'kegg_pathway': 'map00330',
        'reactome_id': None,
        'ec_number': '2.6.1.13',
        'description': 'Ornithine aminotransferase: ornithine -> P5C + glutamate',
    },
    'P5CS_c': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00330',
        'reactome_id': None,
        'ec_number': '2.7.2.11, 1.2.1.41',
        'description': 'P5C synthase: glutamate -> 1-pyrroline-5-carboxylate (lumped)',
    },
    'PYCR_c': {
        'kegg_reaction': 'R01248',
        'kegg_pathway': 'map00330',
        'reactome_id': None,
        'ec_number': '1.5.1.2',
        'description': 'P5C reductase: P5C -> proline',
    },
    'PRODH_c': {
        'kegg_reaction': 'R01253',
        'kegg_pathway': 'map00330',
        'reactome_id': None,
        'ec_number': '1.5.5.2',
        'description': 'Proline dehydrogenase: proline -> P5C',
    },
    'P5CDH_c': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00330',
        'reactome_id': None,
        'ec_number': '1.2.1.88',
        'description': 'P5C dehydrogenase: P5C -> glutamate',
    },
    'ASNS_c': {
        'kegg_reaction': 'R00578',
        'kegg_pathway': 'map00250',
        'reactome_id': None,
        'ec_number': '6.3.5.4',
        'description': 'Asparagine synthetase: aspartate + glutamine -> asparagine + glutamate',
    },
    'ASPG_c': {
        'kegg_reaction': 'R00485',
        'kegg_pathway': 'map00250',
        'reactome_id': None,
        'ec_number': '3.5.1.1',
        'description': 'Asparaginase: asparagine -> aspartate + NH3',
    },
    'SHMT_c': {
        'kegg_reaction': 'R00945',
        'kegg_pathway': 'map00260',
        'reactome_id': None,
        'ec_number': '2.1.2.1',
        'description': 'Serine hydroxymethyltransferase: serine <-> glycine',
    },
    'PURINE_DE_NOVO_c': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00230',
        'reactome_id': None,
        'ec_number': None,
        'description': 'De novo purine synthesis to IMP (lumped)',
    },
    'ADSS_ADSL_c': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00230',
        'reactome_id': None,
        'ec_number': '6.3.4.4, 4.3.2.2',
        'description': 'Adenylosuccinate synthase + lyase: IMP + aspartate -> AMP (lumped)',
    },
    'IMPDH_GMPS_c': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00230',
        'reactome_id': None,
        'ec_number': '1.1.1.205, 6.3.5.2',
        'description': 'IMP dehydrogenase + GMP synthase: IMP -> GMP (lumped)',
    },
    'AMPD_c': {
        'kegg_reaction': 'R00181',
        'kegg_pathway': 'map00230',
        'reactome_id': None,
        'ec_number': '3.5.4.6',
        'description': 'AMP deaminase: AMP -> IMP + NH3',
    },
    'ADA_c': {
        'kegg_reaction': 'R01560',
        'kegg_pathway': 'map00230',
        'reactome_id': None,
        'ec_number': '3.5.4.4',
        'description': 'Adenosine deaminase: adenosine -> inosine + NH3',
    },
    'NT5_AMP_c': {
        'kegg_reaction': 'R00183',
        'kegg_pathway': 'map00230',
        'reactome_id': None,
        'ec_number': '3.1.3.5',
        'description': "5'-Nucleotidase: AMP -> adenosine",
    },
    'NT5_GMP_c': {
        'kegg_reaction': 'R01227',
        'kegg_pathway': 'map00230',
        'reactome_id': None,
        'ec_number': '3.1.3.5',
        'description': "5'-Nucleotidase: GMP -> guanosine",
    },
    'PNP_INO_c': {
        'kegg_reaction': 'R01863',
        'kegg_pathway': 'map00230',
        'reactome_id': None,
        'ec_number': '2.4.2.1',
        'description': 'Purine nucleoside phosphorylase: inosine -> hypoxanthine',
    },
    'PNP_GUO_c': {
        'kegg_reaction': 'R02147',
        'kegg_pathway': 'map00230',
        'reactome_id': None,
        'ec_number': '2.4.2.1',
        'description': 'Purine nucleoside phosphorylase: guanosine -> guanine',
    },
    'PYRIMIDINE_DE_NOVO_c': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00240',
        'reactome_id': None,
        'ec_number': None,
        'description': 'De novo pyrimidine synthesis to UMP (lumped)',
    },
    'CTPS_c': {
        'kegg_reaction': 'R00573',
        'kegg_pathway': 'map00240',
        'reactome_id': None,
        'ec_number': '6.3.4.2',
        'description': 'CTP synthase: UTP + glutamine -> CTP + glutamate (network: UMP -> CMP)',
    },
    'NT5_UMP_c': {
        'kegg_reaction': 'R00963',
        'kegg_pathway': 'map00240',
        'reactome_id': None,
        'ec_number': '3.1.3.5',
        'description': "5'-Nucleotidase: UMP -> uridine",
    },
    'NT5_CMP_c': {
        'kegg_reaction': 'R01548',
        'kegg_pathway': 'map00240',
        'reactome_id': None,
        'ec_number': '3.1.3.5',
        'description': "5'-Nucleotidase: CMP -> cytidine",
    },
    'CDA_c': {
        'kegg_reaction': 'R01878',
        'kegg_pathway': 'map00240',
        'reactome_id': None,
        'ec_number': '3.5.4.5',
        'description': 'Cytidine deaminase: cytidine -> uridine + NH3',
    },
    'RNR_DCYTD_c': {
        'kegg_reaction': 'R02024',
        'kegg_pathway': 'map00240',
        'reactome_id': None,
        'ec_number': '1.17.4.1',
        'description': 'Ribonucleotide reductase: CDP -> dCDP (network: CMP -> deoxycytidine, lumped)',
    },
    'URI_degradation_c': {
        'kegg_reaction': 'R01876',
        'kegg_pathway': 'map00240',
        'reactome_id': None,
        'ec_number': '2.4.2.3',
        'description': 'Uridine phosphorylase: uridine -> uracil (start of pyrimidine degradation)',
    },
    'DCYTD_degradation_c': {
        'kegg_reaction': None,
        'kegg_pathway': 'map00240',
        'reactome_id': None,
        'ec_number': None,
        'description': 'Deoxycytidine catabolism (cytidine deaminase, thymidine phosphorylase; lumped)',
    },
}


# ══════════════════════════════════════════════════════════════════════════════
# Lookup helpers
# ══════════════════════════════════════════════════════════════════════════════

REACTION_ALIASES: Dict[str, str] = {
    # Naming aliases observed in project outputs / model variants
    'ACAC_production_m': 'ACAC_production',
    'ACACCOA_lyase_m': 'ACAC_production',
    'GLYCOLYSIS_c': 'GLYCOLYSIS',
    'GLC_infuse': 'GLC_exchange',
    'LDH_e': 'LDH_c',
    'lipid_mobilization_c': 'lipid_mobilization',
    'BIOMASS_REACTION_2': 'BIOMASS_REACTION_1',
    'BIOMASS_REACTION_3': 'BIOMASS_REACTION_1',
    'BIOMASS_REACTION_4': 'BIOMASS_REACTION_1',
    'BIOMASS_REACTION_5': 'BIOMASS_REACTION_1',
    'BIOMASS_REACTION_7': 'BIOMASS_REACTION_1',
}

_REGISTRY_LOWER = {k.lower(): k for k in REACTION_REGISTRY}
_ALIAS_LOWER = {k.lower(): v for k, v in REACTION_ALIASES.items()}


def _normalize_reaction_id(reaction_id: str) -> str:
    """Normalize raw reaction identifier to the core model reaction key."""
    rid = str(reaction_id).strip()
    if not rid:
        return rid
    if rid.startswith('MIX:'):
        rid = rid[4:]
    if '--' in rid:
        rid = rid.split('--', 1)[1]
    rid = rid.strip()
    if rid.endswith('__R'):
        rid = rid[:-3]
    return rid


def resolve_reaction_id(reaction_id: str) -> Optional[str]:
    """Resolve raw/custom reaction ID to canonical REACTION_REGISTRY key."""
    if reaction_id is None:
        return None
    rid = _normalize_reaction_id(reaction_id)
    if not rid:
        return None

    if rid in REACTION_REGISTRY:
        return rid

    alias_target = REACTION_ALIASES.get(rid)
    if alias_target in REACTION_REGISTRY:
        return alias_target

    rid_lower = rid.lower()
    canonical = _REGISTRY_LOWER.get(rid_lower)
    if canonical:
        return canonical

    alias_target = _ALIAS_LOWER.get(rid_lower)
    if alias_target in REACTION_REGISTRY:
        return alias_target

    return None


def get_reaction_info(reaction_id: str) -> ReactionInfo:
    """Return full annotation dict for a reaction, or empty dict if unknown."""
    resolved = resolve_reaction_id(reaction_id)
    if not resolved:
        return {}
    return REACTION_REGISTRY.get(resolved, {})


def get_kegg_reaction(reaction_id: str) -> Optional[str]:
    """Return KEGG Reaction ID (e.g. 'R00771') or None."""
    return get_reaction_info(reaction_id).get('kegg_reaction')


def get_kegg_pathway(reaction_id: str) -> Optional[str]:
    """Return KEGG Pathway map ID (e.g. 'map00010') or None."""
    return get_reaction_info(reaction_id).get('kegg_pathway')


def get_reactome_id(reaction_id: str) -> Optional[str]:
    """Return Reactome stable ID (e.g. 'R-MMU-70493') or None."""
    return get_reaction_info(reaction_id).get('reactome_id')


def get_ec_number(reaction_id: str) -> Optional[str]:
    """Return EC number string or None."""
    return get_reaction_info(reaction_id).get('ec_number')


def get_description(reaction_id: str) -> str:
    """Return English description or empty string."""
    return get_reaction_info(reaction_id).get('description', '')


def list_unmapped_reactions(reaction_ids: List[str]) -> pd.DataFrame:
    """Return a diagnostics table for reaction IDs that fail registry mapping."""
    rows = []
    for raw in sorted({str(r).strip() for r in reaction_ids if r is not None and str(r).strip()}):
        normalized = _normalize_reaction_id(raw)
        resolved = resolve_reaction_id(raw)
        rows.append({
            'raw_reaction': raw,
            'normalized_reaction': normalized,
            'mapped_reaction': resolved,
            'is_mapped': bool(resolved),
        })
    return pd.DataFrame(rows)


def list_reactions_by_pathway(kegg_pathway: str) -> List[str]:
    """Return all reaction IDs belonging to a KEGG pathway map."""
    return [
        rid for rid, info in REACTION_REGISTRY.items()
        if info.get('kegg_pathway') == kegg_pathway
    ]


def list_reactions_with_ec() -> List[str]:
    """Return all reaction IDs that have an EC number annotation."""
    return [
        rid for rid, info in REACTION_REGISTRY.items()
        if info.get('ec_number') is not None
    ]


def get_all_pathways() -> Dict[str, List[str]]:
    """Return {kegg_pathway_id: [reaction_ids]} for pathway-level grouping."""
    pathways: Dict[str, List[str]] = {}
    for rid, info in REACTION_REGISTRY.items():
        pw = info.get('kegg_pathway')
        if pw:
            pathways.setdefault(pw, []).append(rid)
    return {k: sorted(v) for k, v in pathways.items()}


# ══════════════════════════════════════════════════════════════════════════════
# DataFrame annotation helper (for MFA / multi-omics pipeline output)
# ══════════════════════════════════════════════════════════════════════════════

_ANNOTATION_COLUMNS = ['kegg_reaction', 'kegg_pathway', 'reactome_id', 'ec_number', 'description']


def annotate_flux_dataframe(df: pd.DataFrame, reaction_col: str = 'reaction') -> pd.DataFrame:
    """Annotate a MFA flux DataFrame with public database pathway IDs.

    Given a DataFrame with a ``reaction_col`` containing MFA reaction IDs
    (e.g. 'HEX_c', 'CS_m'), adds columns for kegg_reaction, kegg_pathway,
    reactome_id, ec_number, and description.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame (e.g. differential flux results).
    reaction_col : str
        Column name containing the reaction ID (default: 'reaction').

    Returns
    -------
    pd.DataFrame
        Annotated DataFrame with 5 extra columns appended.
    """
    if reaction_col not in df.columns:
        return df

    annotations = df[reaction_col].map(lambda rid: get_reaction_info(str(rid)) if pd.notna(rid) else {})
    for col in _ANNOTATION_COLUMNS:
        df[col] = annotations.map(lambda info: info.get(col) if isinstance(info, dict) else None)
    return df
