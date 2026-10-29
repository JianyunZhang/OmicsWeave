#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mouse Protein UniProt / EC Number Annotation Dictionary

Central reference for all omics analysis tools. Maps canonical MGI gene symbols
to their corresponding mouse (Mus musculus) UniProt Swiss-Prot accessions,
entry names, and EC numbers.

Usage:
    from omicsweave.common.protein_id_mapping import (
        MOUSE_PROTEIN_UNIPROT_MAPPING,
        get_uniprot_id, get_ec_for_protein,
    )
    info = MOUSE_PROTEIN_UNIPROT_MAPPING.get('Hk1')
    # -> {'uniprot_id': 'P17710', 'uniprot_entry_name': 'HK1_MOUSE', 'ec_number': '2.7.1.1'}

Sources:
    - UniProtKB/Swiss-Prot (reviewed) mouse entries
    - BRENDA enzyme database for EC numbers
    - NCBI Gene for mouse-specific gene symbols
"""

from typing import Dict, Optional

MOUSE_PROTEIN_UNIPROT_MAPPING: Dict[str, Dict[str, Optional[str]]] = {
    # =========================================================================
    # Glycolysis / Gluconeogenesis
    # =========================================================================
    'Gck': {'uniprot_id': 'P52792', 'uniprot_entry_name': 'HK4_MOUSE', 'ec_number': '2.7.1.1'},
    'Hk1': {'uniprot_id': 'P17710', 'uniprot_entry_name': 'HK1_MOUSE', 'ec_number': '2.7.1.1'},
    'Hk2': {'uniprot_id': 'O08528', 'uniprot_entry_name': 'HK2_MOUSE', 'ec_number': '2.7.1.1'},
    'Hk3': {'uniprot_id': 'Q3TRM8', 'uniprot_entry_name': 'HK3_MOUSE', 'ec_number': '2.7.1.1'},
    'Hkdc1': {'uniprot_id': 'Q6PUI6', 'uniprot_entry_name': 'HKDC1_MOUSE', 'ec_number': '2.7.1.1'},
    'Gckr': {'uniprot_id': 'Q91X44', 'uniprot_entry_name': 'GCKR_MOUSE', 'ec_number': None},
    'Gpi': {'uniprot_id': 'P06745', 'uniprot_entry_name': 'G6PI_MOUSE', 'ec_number': '5.3.1.9'},
    'Gpi1': {'uniprot_id': 'P06745', 'uniprot_entry_name': 'G6PI_MOUSE', 'ec_number': '5.3.1.9'},
    'Pfkp': {'uniprot_id': 'Q9WUA3', 'uniprot_entry_name': 'PFKP_MOUSE', 'ec_number': '2.7.1.11'},
    'Pfkl': {'uniprot_id': 'P12382', 'uniprot_entry_name': 'PFKL_MOUSE', 'ec_number': '2.7.1.11'},
    'Pfkm': {'uniprot_id': 'P47857', 'uniprot_entry_name': 'PFKM_MOUSE', 'ec_number': '2.7.1.11'},
    'Pfkfb1': {'uniprot_id': 'A2AFM9', 'uniprot_entry_name': 'PFKFB1_MOUSE', 'ec_number': '2.7.1.105'},
    'Pfkfb2': {'uniprot_id': 'P70265', 'uniprot_entry_name': 'F262_MOUSE', 'ec_number': '2.7.1.105; 3.1.3.46'},
    'Pfkfb3': {'uniprot_id': 'A7UAK8', 'uniprot_entry_name': 'PFKFB3_MOUSE', 'ec_number': '2.7.1.105'},
    'Pfkfb4': {'uniprot_id': 'Q3UMR5', 'uniprot_entry_name': 'PFKFB4_MOUSE', 'ec_number': '2.7.1.105'},
    'Aldoa': {'uniprot_id': 'P05064', 'uniprot_entry_name': 'ALDOA_MOUSE', 'ec_number': '4.1.2.13'},
    'Aldob': {'uniprot_id': 'Q91Y97', 'uniprot_entry_name': 'ALDOB_MOUSE', 'ec_number': '4.1.2.13'},
    'Aldoc': {'uniprot_id': 'P05063', 'uniprot_entry_name': 'ALDOC_MOUSE', 'ec_number': '4.1.2.13'},
    'Tpi1': {'uniprot_id': 'P17751', 'uniprot_entry_name': 'TPIS_MOUSE', 'ec_number': '5.3.1.1'},
    'Gapdh': {'uniprot_id': 'P16858', 'uniprot_entry_name': 'G3P_MOUSE', 'ec_number': '1.2.1.12'},
    'Gapdhs': {'uniprot_id': 'Q9JLR4', 'uniprot_entry_name': 'GPDS_MOUSE', 'ec_number': '1.2.1.12'},
    'Pgk1': {'uniprot_id': 'P09411', 'uniprot_entry_name': 'PGK1_MOUSE', 'ec_number': '2.7.2.3'},
    'Pgk2': {'uniprot_id': 'P09042', 'uniprot_entry_name': 'PGK2_MOUSE', 'ec_number': '2.7.2.3'},
    'Pgam1': {'uniprot_id': 'Q9DBJ1', 'uniprot_entry_name': 'PGM1_MOUSE', 'ec_number': '5.4.2.12'},
    'Pgam2': {'uniprot_id': 'P15259', 'uniprot_entry_name': 'PGM2_MOUSE', 'ec_number': '5.4.2.12'},
    'Eno1': {'uniprot_id': 'P17182', 'uniprot_entry_name': 'ENOA_MOUSE', 'ec_number': '4.2.1.11'},
    'Eno2': {'uniprot_id': 'P17183', 'uniprot_entry_name': 'ENOG_MOUSE', 'ec_number': '4.2.1.11'},
    'Eno3': {'uniprot_id': 'P21550', 'uniprot_entry_name': 'ENOB_MOUSE', 'ec_number': '4.2.1.11'},
    'Pkm': {'uniprot_id': 'P52480', 'uniprot_entry_name': 'KPYM_MOUSE', 'ec_number': '2.7.1.40'},
    'Pklr': {'uniprot_id': 'P53657', 'uniprot_entry_name': 'KPYR_MOUSE', 'ec_number': '2.7.1.40'},
    'Ldha': {'uniprot_id': 'P06151', 'uniprot_entry_name': 'LDHA_MOUSE', 'ec_number': '1.1.1.27'},
    'Ldhb': {'uniprot_id': 'P16125', 'uniprot_entry_name': 'LDHB_MOUSE', 'ec_number': '1.1.1.27'},
    'Ldhc': {'uniprot_id': 'P00336', 'uniprot_entry_name': 'LDHC_MOUSE', 'ec_number': '1.1.1.27'},
    'Pgm1': {'uniprot_id': 'Q9D0F9', 'uniprot_entry_name': 'PGM1_MOUSE', 'ec_number': '5.4.2.2'},
    'Pgm2': {'uniprot_id': 'Q8VCT4', 'uniprot_entry_name': 'PGM2_MOUSE', 'ec_number': '5.4.2.2'},
    'Pgm3': {'uniprot_id': 'Q9EQ09', 'uniprot_entry_name': 'PGM3_MOUSE', 'ec_number': '5.4.2.2'},
    'Ugp2': {'uniprot_id': 'Q91ZV5', 'uniprot_entry_name': 'UGP2_MOUSE', 'ec_number': '2.7.7.9'},
    'Gys1': {'uniprot_id': 'Q9Z1E4', 'uniprot_entry_name': 'GYS1_MOUSE', 'ec_number': '2.4.1.11'},
    'Gys2': {'uniprot_id': 'Q9WUB5', 'uniprot_entry_name': 'GYS2_MOUSE', 'ec_number': '2.4.1.11'},
    'G6pc': {'uniprot_id': 'P35576', 'uniprot_entry_name': 'G6PC1_MOUSE', 'ec_number': '3.1.3.9'},
    'G6pc2': {'uniprot_id': 'Q9JL03', 'uniprot_entry_name': 'G6PC2_MOUSE', 'ec_number': '3.1.3.9'},
    'G6pc3': {'uniprot_id': 'Q9CYD6', 'uniprot_entry_name': 'G6PC3_MOUSE', 'ec_number': '3.1.3.9'},
    'Pck1': {'uniprot_id': 'Q9Z2V4', 'uniprot_entry_name': 'PCK1_MOUSE', 'ec_number': '4.1.1.32'},
    'Pck2': {'uniprot_id': 'Q8BYR3', 'uniprot_entry_name': 'PCK2_MOUSE', 'ec_number': '4.1.1.32'},
    'Fbp1': {'uniprot_id': 'Q9QXD6', 'uniprot_entry_name': 'FBP1_MOUSE', 'ec_number': '3.1.3.11'},
    'Fbp2': {'uniprot_id': 'P70695', 'uniprot_entry_name': 'F16P2_MOUSE', 'ec_number': '3.1.3.11'},
    'Pygb': {'uniprot_id': 'Q9WUB3', 'uniprot_entry_name': 'PYGB_MOUSE', 'ec_number': '2.4.1.1'},
    'Pygm': {'uniprot_id': 'Q9WUB2', 'uniprot_entry_name': 'PYGM_MOUSE', 'ec_number': '2.4.1.1'},
    'Pygl': {'uniprot_id': 'O70172', 'uniprot_entry_name': 'PYGL_MOUSE', 'ec_number': '2.4.1.1'},
    'Gbe1': {'uniprot_id': 'Q9R1I0', 'uniprot_entry_name': 'GBE1_MOUSE', 'ec_number': '2.4.1.181'},
    'Gsk3a': {'uniprot_id': 'Q2NL51', 'uniprot_entry_name': 'GSK3A_MOUSE', 'ec_number': '2.7.11.26'},
    'Gsk3b': {'uniprot_id': 'Q9WV60', 'uniprot_entry_name': 'GSK3B_MOUSE', 'ec_number': '2.7.11.26'},
    'Phka1': {'uniprot_id': 'P42816', 'uniprot_entry_name': 'PHKA1_MOUSE', 'ec_number': '2.7.11.19'},
    'Phka2': {'uniprot_id': 'O08795', 'uniprot_entry_name': 'PHKA2_MOUSE', 'ec_number': '2.7.11.19'},
    'Phkb': {'uniprot_id': 'Q7SIH1', 'uniprot_entry_name': 'PHKB_MOUSE', 'ec_number': None},
    'Phkg1': {'uniprot_id': 'O88664', 'uniprot_entry_name': 'PHKG1_MOUSE', 'ec_number': '2.7.11.19'},
    'Phkg2': {'uniprot_id': 'Q9WUA4', 'uniprot_entry_name': 'PHKG2_MOUSE', 'ec_number': '2.7.11.19'},

    # =========================================================================
    # TCA Cycle
    # =========================================================================
    'Cs': {'uniprot_id': 'Q9CZU6', 'uniprot_entry_name': 'CISY_MOUSE', 'ec_number': '2.3.3.1'},
    'Aco1': {'uniprot_id': 'P28271', 'uniprot_entry_name': 'ACON1_MOUSE', 'ec_number': '4.2.1.3'},
    'Aco2': {'uniprot_id': 'Q99KI0', 'uniprot_entry_name': 'ACON2_MOUSE', 'ec_number': '4.2.1.3'},
    'Idh1': {'uniprot_id': 'O88844', 'uniprot_entry_name': 'IDHC_MOUSE', 'ec_number': '1.1.1.42'},
    'Idh2': {'uniprot_id': 'P54071', 'uniprot_entry_name': 'IDHP_MOUSE', 'ec_number': '1.1.1.42'},
    'Idh3a': {'uniprot_id': 'Q9D6R2', 'uniprot_entry_name': 'IDH3A_MOUSE', 'ec_number': '1.1.1.41'},
    'Idh3b': {'uniprot_id': 'Q9QYX4', 'uniprot_entry_name': 'IDH3B_MOUSE', 'ec_number': '1.1.1.41'},
    'Idh3g': {'uniprot_id': 'Q9JHI5', 'uniprot_entry_name': 'IDH3G_MOUSE', 'ec_number': '1.1.1.41'},
    'Ogdh': {'uniprot_id': 'Q60597', 'uniprot_entry_name': 'ODO1_MOUSE', 'ec_number': '1.2.4.2'},
    'Ogdhl': {'uniprot_id': 'Q99KI7', 'uniprot_entry_name': 'ODOHL_MOUSE', 'ec_number': None},
    'Dlst': {'uniprot_id': 'Q9D2G2', 'uniprot_entry_name': 'ODO2_MOUSE', 'ec_number': '2.3.1.61'},
    'Dld': {'uniprot_id': 'O08749', 'uniprot_entry_name': 'DLDH_MOUSE', 'ec_number': '1.8.1.4'},
    'Suclg1': {'uniprot_id': 'Q9WUM5', 'uniprot_entry_name': 'SUCLG1_MOUSE', 'ec_number': '6.2.1.5'},
    'Suclg2': {'uniprot_id': 'Q9Z2I2', 'uniprot_entry_name': 'SUCLG2_MOUSE', 'ec_number': '6.2.1.5'},
    'Sucla2': {'uniprot_id': 'Q9JHR9', 'uniprot_entry_name': 'SUCLA2_MOUSE', 'ec_number': '6.2.1.5'},
    'Sdha': {'uniprot_id': 'Q8K2B3', 'uniprot_entry_name': 'SDHA_MOUSE', 'ec_number': '1.3.5.1'},
    'Sdhb': {'uniprot_id': 'Q9CQA3', 'uniprot_entry_name': 'SDHB_MOUSE', 'ec_number': '1.3.5.1'},
    'Sdhc': {'uniprot_id': 'Q9CXV1', 'uniprot_entry_name': 'SDHC_MOUSE', 'ec_number': '1.3.5.1'},
    'Sdhd': {'uniprot_id': 'Q9ERS2', 'uniprot_entry_name': 'SDHD_MOUSE', 'ec_number': '1.3.5.1'},
    'Sdhaf1': {'uniprot_id': 'Q9DAQ0', 'uniprot_entry_name': 'SDAF1_MOUSE', 'ec_number': None},
    'Sdhaf2': {'uniprot_id': 'Q8C6W8', 'uniprot_entry_name': 'SDAF2_MOUSE', 'ec_number': None},
    'Fh1': {'uniprot_id': 'P97807', 'uniprot_entry_name': 'FUMH_MOUSE', 'ec_number': '4.2.1.2'},
    'Mdh1': {'uniprot_id': 'P14152', 'uniprot_entry_name': 'MDHC_MOUSE', 'ec_number': '1.1.1.37'},
    'Mdh2': {'uniprot_id': 'P08249', 'uniprot_entry_name': 'MDHM_MOUSE', 'ec_number': '1.1.1.37'},
    'Ndufs1': {'uniprot_id': 'Q91VD9', 'uniprot_entry_name': 'NDUS1_MOUSE', 'ec_number': '7.1.1.2'},
    'Ndufv1': {'uniprot_id': 'Q8C0Q2', 'uniprot_entry_name': 'NDUV1_MOUSE', 'ec_number': '7.1.1.2'},

    # =========================================================================
    # Amino Acid Metabolism
    # =========================================================================
    'Got1': {'uniprot_id': 'P05201', 'uniprot_entry_name': 'AATC_MOUSE', 'ec_number': '2.6.1.1'},
    'Got2': {'uniprot_id': 'P05202', 'uniprot_entry_name': 'AATM_MOUSE', 'ec_number': '2.6.1.1'},
    'Gpt': {'uniprot_id': 'P24288', 'uniprot_entry_name': 'ALAT1_MOUSE', 'ec_number': '2.6.1.2'},
    'Gpt2': {'uniprot_id': 'Q8R4Y4', 'uniprot_entry_name': 'ALAT2_MOUSE', 'ec_number': '2.6.1.2'},
    'Glud1': {'uniprot_id': 'P26443', 'uniprot_entry_name': 'DHE3_MOUSE', 'ec_number': '1.4.1.3'},
    'Glul': {'uniprot_id': 'P15105', 'uniprot_entry_name': 'GLNA_MOUSE', 'ec_number': '6.3.1.2'},
    'Gls': {'uniprot_id': 'Q9Z2Z6', 'uniprot_entry_name': 'GLSK_MOUSE', 'ec_number': '3.5.1.2'},
    'Gls2': {'uniprot_id': 'Q9QXN8', 'uniprot_entry_name': 'GLS2_MOUSE', 'ec_number': '3.5.1.2'},
    'Phgdh': {'uniprot_id': 'Q61763', 'uniprot_entry_name': 'SERAD_MOUSE', 'ec_number': '1.1.1.95'},
    'Psat1': {'uniprot_id': 'Q8BGQ7', 'uniprot_entry_name': 'PSAT1_MOUSE', 'ec_number': '2.6.1.52'},
    'Gldc': {'uniprot_id': 'Q9QXF2', 'uniprot_entry_name': 'GLDC_MOUSE', 'ec_number': '1.4.4.2'},
    'Amt': {'uniprot_id': 'Q9CT18', 'uniprot_entry_name': 'AMT_MOUSE', 'ec_number': '2.1.2.10'},
    'Cps1': {'uniprot_id': 'Q8C196', 'uniprot_entry_name': 'CPSM_MOUSE', 'ec_number': '6.3.4.16'},
    'Nags': {'uniprot_id': 'Q8BMA4', 'uniprot_entry_name': 'NAGS_MOUSE', 'ec_number': '2.7.2.8'},
    'Otc': {'uniprot_id': 'P11725', 'uniprot_entry_name': 'OTC_MOUSE', 'ec_number': '2.1.3.3'},

    # =========================================================================
    # Pyruvate Dehydrogenase Complex
    # =========================================================================
    'Pdha1': {'uniprot_id': 'P35486', 'uniprot_entry_name': 'ODPA_MOUSE', 'ec_number': '1.2.4.1'},
    'Pdhb': {'uniprot_id': 'Q9D051', 'uniprot_entry_name': 'ODPB_MOUSE', 'ec_number': '1.2.4.1'},
    'Dlat': {'uniprot_id': 'Q8BMF4', 'uniprot_entry_name': 'ODP2_MOUSE', 'ec_number': '2.3.1.12'},
    'Pdk1': {'uniprot_id': 'Q3U1D8', 'uniprot_entry_name': 'PDK1_MOUSE', 'ec_number': '2.7.11.2'},
    'Pdk2': {'uniprot_id': 'Q8BKJ9', 'uniprot_entry_name': 'PDK2_MOUSE', 'ec_number': '2.7.11.2'},
    'Pdk4': {'uniprot_id': 'Q166J5', 'uniprot_entry_name': 'PDK4_MOUSE', 'ec_number': '2.7.11.2'},
    'Pc': {'uniprot_id': 'Q05920', 'uniprot_entry_name': 'PYC_MOUSE', 'ec_number': '6.4.1.1'},
    'Mpc1': {'uniprot_id': 'Q8VC05', 'uniprot_entry_name': 'MPC1_MOUSE', 'ec_number': None},
    'Mpc2': {'uniprot_id': 'Q3UZH3', 'uniprot_entry_name': 'MPC2_MOUSE', 'ec_number': None},

    # =========================================================================
    # Fatty Acid / Lipid Metabolism
    # =========================================================================
    'Cpt1a': {'uniprot_id': 'P97742', 'uniprot_entry_name': 'CPT1A_MOUSE', 'ec_number': '2.3.1.21'},
    'Cpt1b': {'uniprot_id': 'Q924X2', 'uniprot_entry_name': 'CPT1B_MOUSE', 'ec_number': '2.3.1.21'},
    'Cpt1c': {'uniprot_id': 'Q8NDH6', 'uniprot_entry_name': 'CPT1C_MOUSE', 'ec_number': None},
    'Cpt2': {'uniprot_id': 'P52825', 'uniprot_entry_name': 'CPT2_MOUSE', 'ec_number': '2.3.1.21'},
    'Acadm': {'uniprot_id': 'P45952', 'uniprot_entry_name': 'ACADM_MOUSE', 'ec_number': '1.3.8.1'},
    'Acadvl': {'uniprot_id': 'P50544', 'uniprot_entry_name': 'ACADV_MOUSE', 'ec_number': '1.3.8.9'},
    'Hadha': {'uniprot_id': 'Q8BMS1', 'uniprot_entry_name': 'ECHA_MOUSE', 'ec_number': '2.3.1.16'},
    'Hadhb': {'uniprot_id': 'Q99JY0', 'uniprot_entry_name': 'ECHB_MOUSE', 'ec_number': '2.3.1.16'},
    'Acsl1': {'uniprot_id': 'P41233', 'uniprot_entry_name': 'ACSL1_MOUSE', 'ec_number': '6.2.1.3'},
    'Acsl3': {'uniprot_id': 'Q99NB1', 'uniprot_entry_name': 'ACSL3_MOUSE', 'ec_number': '6.2.1.3'},
    'Acsl4': {'uniprot_id': 'O35551', 'uniprot_entry_name': 'ACSL4_MOUSE', 'ec_number': '6.2.1.3'},
    'Acss1': {'uniprot_id': 'Q9CRB4', 'uniprot_entry_name': 'ACS1_MOUSE', 'ec_number': '6.2.1.1'},
    'Acss2': {'uniprot_id': 'Q99NB7', 'uniprot_entry_name': 'ACS2_MOUSE', 'ec_number': '6.2.1.1'},
    'Acly': {'uniprot_id': 'Q91V92', 'uniprot_entry_name': 'ACLY_MOUSE', 'ec_number': '2.3.3.8'},
    'Acaca': {'uniprot_id': 'Q5SWU9', 'uniprot_entry_name': 'ACACA_MOUSE', 'ec_number': '6.4.1.2'},
    'Acacb': {'uniprot_id': 'E9Q4Z2', 'uniprot_entry_name': 'ACACB_MOUSE', 'ec_number': '6.4.1.2'},
    'Fasn': {'uniprot_id': 'P19096', 'uniprot_entry_name': 'FAS_MOUSE', 'ec_number': '2.3.1.85'},
    'Scd1': {'uniprot_id': 'P13516', 'uniprot_entry_name': 'SCD1_MOUSE', 'ec_number': '1.14.19.1'},
    'Scd2': {'uniprot_id': 'P13011', 'uniprot_entry_name': 'SCD2_MOUSE', 'ec_number': '1.14.19.1'},
    'Fads1': {'uniprot_id': 'Q920L5', 'uniprot_entry_name': 'FADS1_MOUSE', 'ec_number': '1.14.19.44'},
    'Fads2': {'uniprot_id': 'Q9Z0R1', 'uniprot_entry_name': 'FADS2_MOUSE', 'ec_number': '1.14.19.44'},
    'Srebf1': {'uniprot_id': 'Q9WTN3', 'uniprot_entry_name': 'SRBP1_MOUSE', 'ec_number': None},
    'Srebf2': {'uniprot_id': 'Q3U1N2', 'uniprot_entry_name': 'SRBP2_MOUSE', 'ec_number': None},
    'Chrebp': {'uniprot_id': 'Q8R1D0', 'uniprot_entry_name': 'MLXI_MOUSE', 'ec_number': None},
    'Cd36': {'uniprot_id': 'Q08852', 'uniprot_entry_name': 'CD36_MOUSE', 'ec_number': None},
    'Dgat1': {'uniprot_id': 'Q9Z2A7', 'uniprot_entry_name': 'DGAT1_MOUSE', 'ec_number': '2.3.1.20'},
    'Dgat2': {'uniprot_id': 'Q9DCV3', 'uniprot_entry_name': 'DGAT2_MOUSE', 'ec_number': '2.3.1.20'},
    'Pnpla2': {'uniprot_id': 'Q8BJ56', 'uniprot_entry_name': 'PLPL2_MOUSE', 'ec_number': '3.1.1.3'},
    'Lipe': {'uniprot_id': 'P54310', 'uniprot_entry_name': 'LIPE_MOUSE', 'ec_number': '3.1.1.79'},
    'Mgll': {'uniprot_id': 'O35678', 'uniprot_entry_name': 'MGLL_MOUSE', 'ec_number': '3.1.1.23'},

    # =========================================================================
    # Ketone Body Metabolism
    # =========================================================================
    'Hmgcs2': {'uniprot_id': 'P54869', 'uniprot_entry_name': 'HMCS2_MOUSE', 'ec_number': '2.3.3.10'},
    'Hmgcl': {'uniprot_id': 'P38060', 'uniprot_entry_name': 'HMCL_MOUSE', 'ec_number': '4.1.3.4'},
    'Bdh1': {'uniprot_id': 'Q9JHI4', 'uniprot_entry_name': 'BDH_MOUSE', 'ec_number': '1.1.1.30'},
    'Bdh2': {'uniprot_id': 'Q9D894', 'uniprot_entry_name': 'BDH2_MOUSE', 'ec_number': '1.1.1.30'},
    'Oxct1': {'uniprot_id': 'Q9D0K2', 'uniprot_entry_name': 'SCOT1_MOUSE', 'ec_number': '2.8.3.5'},

    # =========================================================================
    # Pentose Phosphate Pathway
    # =========================================================================
    'G6pd': {'uniprot_id': 'Q790Y8', 'uniprot_entry_name': 'G6PDX_MOUSE', 'ec_number': '1.1.1.49'},
    'Pgl': {'uniprot_id': 'Q9JLK3', 'uniprot_entry_name': '6PGL_MOUSE', 'ec_number': '3.1.1.31'},
    'Rgn': {'uniprot_id': 'Q9JLJ8', 'uniprot_entry_name': 'REGN_MOUSE', 'ec_number': None},
    'Pgd': {'uniprot_id': 'Q9DCD0', 'uniprot_entry_name': '6PGD_MOUSE', 'ec_number': '1.1.1.44'},
    'Tkt': {'uniprot_id': 'P40142', 'uniprot_entry_name': 'TKT_MOUSE', 'ec_number': '2.2.1.1'},
    'Taldo1': {'uniprot_id': 'Q93092', 'uniprot_entry_name': 'TALDO_MOUSE', 'ec_number': '2.2.1.2'},
    'Rpi': {'uniprot_id': 'Q8VBT5', 'uniprot_entry_name': 'RPIA_MOUSE', 'ec_number': '5.3.1.6'},
    'Rpe': {'uniprot_id': 'Q8CAX6', 'uniprot_entry_name': 'RPE_MOUSE', 'ec_number': '5.1.3.1'},

    # =========================================================================
    # Glycerol Metabolism
    # =========================================================================
    'Gpd1': {'uniprot_id': 'P28847', 'uniprot_entry_name': 'GPD1_MOUSE', 'ec_number': '1.1.1.8'},
    'Gpd2': {'uniprot_id': 'Q64521', 'uniprot_entry_name': 'GPD2_MOUSE', 'ec_number': '1.1.5.3'},
    'Gk': {'uniprot_id': 'Q64516', 'uniprot_entry_name': 'GLPK_MOUSE', 'ec_number': '2.7.1.30'},
    'Gyk': {'uniprot_id': 'Q9QK41', 'uniprot_entry_name': 'GYK2_MOUSE', 'ec_number': '2.7.1.30'},

    # =========================================================================
    # Malic Enzyme
    # =========================================================================
    'Me1': {'uniprot_id': 'P06801', 'uniprot_entry_name': 'MDHC_MOUSE', 'ec_number': '1.1.1.40'},
    'Me2': {'uniprot_id': 'Q9ERT9', 'uniprot_entry_name': 'MDHM_MOUSE', 'ec_number': '1.1.1.38'},
    'Me3': {'uniprot_id': 'Q8BWF0', 'uniprot_entry_name': 'MDHP_MOUSE', 'ec_number': '1.1.1.40'},

    # =========================================================================
    # Transporters (GLUT, MCT, Mitochondrial carriers, etc.)
    # =========================================================================
    'Slc2a1': {'uniprot_id': 'P17809', 'uniprot_entry_name': 'GTR1_MOUSE', 'ec_number': None},
    'Slc2a2': {'uniprot_id': 'P14246', 'uniprot_entry_name': 'GTR2_MOUSE', 'ec_number': None},
    'Slc2a3': {'uniprot_id': 'P20472', 'uniprot_entry_name': 'GTR3_MOUSE', 'ec_number': None},
    'Slc2a4': {'uniprot_id': 'P14142', 'uniprot_entry_name': 'GTR4_MOUSE', 'ec_number': None},
    'Slc2a5': {'uniprot_id': 'Q9WV38', 'uniprot_entry_name': 'GTR5_MOUSE', 'ec_number': None},
    'Slc2a8': {'uniprot_id': 'Q9JID3', 'uniprot_entry_name': 'GTR8_MOUSE', 'ec_number': None},
    'Slc16a1': {'uniprot_id': 'P53986', 'uniprot_entry_name': 'MOT1_MOUSE', 'ec_number': None},
    'Slc16a3': {'uniprot_id': 'O70417', 'uniprot_entry_name': 'MOT3_MOUSE', 'ec_number': None},
    'Slc16a7': {'uniprot_id': 'Q8VBT6', 'uniprot_entry_name': 'MOT2_MOUSE', 'ec_number': None},
    'Slc25a1': {'uniprot_id': 'P32089', 'uniprot_entry_name': 'SCMC1_MOUSE', 'ec_number': None},
    'Slc25a10': {'uniprot_id': 'Q9DCS9', 'uniprot_entry_name': 'DIC_MOUSE', 'ec_number': None},
    'Slc25a11': {'uniprot_id': 'P56514', 'uniprot_entry_name': 'OGC_MOUSE', 'ec_number': None},
    'Slc25a12': {'uniprot_id': 'Q8BH29', 'uniprot_entry_name': 'AGC1_MOUSE', 'ec_number': None},
    'Slc25a13': {'uniprot_id': 'Q9QYZ4', 'uniprot_entry_name': 'AGC2_MOUSE', 'ec_number': None},
    'Slc25a2': {'uniprot_id': 'Q8C7K3', 'uniprot_entry_name': 'ORC2_MOUSE', 'ec_number': None},
    'Slc25a15': {'uniprot_id': 'O35722', 'uniprot_entry_name': 'ORC1_MOUSE', 'ec_number': None},
    'Slc27a1': {'uniprot_id': 'Q60628', 'uniprot_entry_name': 'SL27A1_MOUSE', 'ec_number': None},
    'Slc27a2': {'uniprot_id': 'O08912', 'uniprot_entry_name': 'SL27A2_MOUSE', 'ec_number': None},
    'Slc27a4': {'uniprot_id': 'Q60644', 'uniprot_entry_name': 'SL27A4_MOUSE', 'ec_number': None},
    'Slc27a5': {'uniprot_id': 'Q8VCM7', 'uniprot_entry_name': 'SL27A5_MOUSE', 'ec_number': None},
    'Slc1a1': {'uniprot_id': 'P51906', 'uniprot_entry_name': 'EAA3_MOUSE', 'ec_number': None},
    'Slc1a2': {'uniprot_id': 'P43006', 'uniprot_entry_name': 'EAA2_MOUSE', 'ec_number': None},
    'Slc1a3': {'uniprot_id': 'P56564', 'uniprot_entry_name': 'EAA1_MOUSE', 'ec_number': None},
    'Slc1a5': {'uniprot_id': 'Q8C3J4', 'uniprot_entry_name': 'AAAT_MOUSE', 'ec_number': None},
    'Slc5a8': {'uniprot_id': 'Q8R4I0', 'uniprot_entry_name': 'S5A8_MOUSE', 'ec_number': None},
    'Slc7a1': {'uniprot_id': 'P28568', 'uniprot_entry_name': 'CAT1_MOUSE', 'ec_number': None},
    'Slc7a2': {'uniprot_id': 'P30684', 'uniprot_entry_name': 'CAT2_MOUSE', 'ec_number': None},
    'Slc7a5': {'uniprot_id': 'Q60629', 'uniprot_entry_name': 'LAT1_MOUSE', 'ec_number': None},
    'Slc7a8': {'uniprot_id': 'Q9QXW9', 'uniprot_entry_name': 'LAT2_MOUSE', 'ec_number': None},
    'Slc38a1': {'uniprot_id': 'Q8K0T5', 'uniprot_entry_name': 'S38A1_MOUSE', 'ec_number': None},
    'Slc38a2': {'uniprot_id': 'Q8BH59', 'uniprot_entry_name': 'S38A2_MOUSE', 'ec_number': None},
    'Slc38a3': {'uniprot_id': 'Q8VCT8', 'uniprot_entry_name': 'S38A3_MOUSE', 'ec_number': None},
    'Slc6a9': {'uniprot_id': 'P46460', 'uniprot_entry_name': 'GLYT1_MOUSE', 'ec_number': None},
    'Slc37a1': {'uniprot_id': 'Q8C0N9', 'uniprot_entry_name': 'S37A1_MOUSE', 'ec_number': None},
    'Slc37a4': {'uniprot_id': 'Q64471', 'uniprot_entry_name': 'S37A4_MOUSE', 'ec_number': None},
    'Aqp3': {'uniprot_id': 'Q3U7U3', 'uniprot_entry_name': 'AQP3_MOUSE', 'ec_number': None},
    'Aqp7': {'uniprot_id': 'Q6P582', 'uniprot_entry_name': 'AQP7_MOUSE', 'ec_number': None},
    'Aqp9': {'uniprot_id': 'Q9QUS3', 'uniprot_entry_name': 'AQP9_MOUSE', 'ec_number': None},
}


def get_uniprot_id(gene_symbol: str) -> Optional[str]:
    """Return UniProt accession for a mouse gene symbol, or None."""
    entry = MOUSE_PROTEIN_UNIPROT_MAPPING.get(gene_symbol)
    return entry['uniprot_id'] if entry else None


def get_ec_for_protein(gene_symbol: str) -> Optional[str]:
    """Return EC number for a mouse gene symbol, or None."""
    entry = MOUSE_PROTEIN_UNIPROT_MAPPING.get(gene_symbol)
    return entry['ec_number'] if entry else None
