#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LIPID MAPS 简写的解析（Liebisch 等 2020；语法参照 Goslin，Kopczynski 等 2020）。

支持厂商表中常见的写法::

    TG(16:0_18:1_18:2)   TG(52:3)          TG 16:0_18:1_18:2   PC(16:0/18:1)   PC 34:1
    PE(O-16:0/22:6)      PE(P-18:0_20:4)   LPC(16:0)           Cer(d18:1/24:0) SM(d18:1/16:0)
    Cer 18:1;O2/24:0     CE(18:2)          FFA(18:2)           CAR(16:0)       DG(16:0_18:1)

返回类别（醚脂按 O- / P- 记为 ``PE-O`` / ``PE-P``，与 :func:`omicsweave.metabolomics.derived_features.lipid_class` 一致）、
各链的碳数与双键数、总碳数与总双键数，以及注释层级（只有总组成时为 ``species``，有各链组成时为 ``molecular_species``）。
"""
from __future__ import annotations

import re
from dataclasses import dataclass

# 每个分子中的酰基 / 烷基链数（鞘脂计鞘氨醇碱基与 N-酰基两条链）
CHAINS_PER_CLASS = {
    "TG": 3, "DG": 2, "MG": 1, "PC": 2, "PE": 2, "PS": 2, "PG": 2, "PI": 2, "PA": 2, "BMP": 2, "CL": 4,
    "LPC": 1, "LPE": 1, "LPS": 1, "LPG": 1, "LPI": 1, "LPA": 1, "CE": 1, "FFA": 1, "FA": 1, "CAR": 1, "ACar": 1,
    "Cer": 2, "SM": 2, "HexCer": 2, "Hex1Cer": 2, "Hex2Cer": 2, "LacCer": 2, "SPH": 1, "MGDG": 2, "DGDG": 2,
    "WE": 2, "NAE": 1, "PMe": 2, "PEt": 2,
}
_HEAD = re.compile(r"^\s*(?P<cls>[A-Za-z][A-Za-z0-9]*?)\s*(?:\((?P<inner>.*)\)(?:\[[^\]]*\])?|\s+(?P<bare>\S.*))\s*$")
# 链上的修饰（羟基 (OH) / (2OH)、新命名法的 ;O2）不影响碳数与双键数
# 链后缀 n3 / e（ω 位置、醚键）与尾随的 [U] 等标注忽略；N- 为 N-酰基链（如 LNAPE 的 N-20:3）；厂商表偶有把 : 写成 . 的（PI(17:0_20.4)）
_CHAIN = re.compile(r"^(?P<ether>[OP]-)?(?:N-)?(?P<base>[dtm])?(?P<c>\d{1,2})[:.](?P<db>\d{1,2})(?:n-?\d+|e|p)?"
                    r"(?P<rest>(?:;[^/_()]*|\([^()]*\))*)$")


@dataclass(frozen=True)
class LipidName:
    lipid_class: str
    chains: tuple[tuple[int, int], ...]
    ether: str | None
    level: str

    @property
    def carbons(self) -> int:
        return sum(c for c, _ in self.chains)

    @property
    def double_bonds(self) -> int:
        return sum(d for _, d in self.chains)

    @property
    def n_chains(self) -> int:
        return CHAINS_PER_CLASS.get(self.lipid_class.split("-")[0], len(self.chains) or 1)

    @property
    def mean_chain_length(self) -> float:
        return self.carbons / max(self.n_chains, 1)


def parse_lipid(name: object) -> LipidName | None:
    """解析一个脂质名称；不是脂质简写时返回 None。"""
    if name is None:
        return None
    match = _HEAD.match(str(name))
    if not match:
        return None
    body = match["inner"] if match["inner"] is not None else match["bare"]
    pieces = [p.strip() for p in re.split(r"[/_]", body) if p.strip()]
    if not pieces:
        return None
    chains, ether = [], None
    for index, piece in enumerate(pieces):
        chain = _CHAIN.match(piece)
        if not chain:
            return None
        if chain["ether"] and ether is None:
            ether = chain["ether"][0]
        chains.append((int(chain["c"]), int(chain["db"])))
    cls = match["cls"]
    expected = CHAINS_PER_CLASS.get(cls)
    level = "molecular_species" if len(chains) > 1 or expected == 1 else "species"
    return LipidName(cls + (f"-{ether}" if ether else ""), tuple(chains), ether, level)


def lipid_class(name: object) -> str | None:
    parsed = parse_lipid(name)
    return parsed.lipid_class if parsed else None
