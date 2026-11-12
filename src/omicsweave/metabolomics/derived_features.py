#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
代谢组的派生特征：脂质类别总量与质量作用比。两者都在预处理后的 log2 归一化峰面积上按样本计算，
再与测得的代谢物用同一个 Welch t 检验，但各自单独做多重检验校正（不并入测得代谢物的 FDR）。

脂质类别总量（``feature_type = lipid_class_total``）
    按 LIPID MAPS 简写（``TG(16:0_18:1_18:2)``、``PE(O-16:0_22:6)``、``Cer(d18:1/24:0)``）从名称解析类别，
    醚脂按 O- / P- 单列（PE-O、PE-P）；同一类别的分子种在线性尺度上求和后取 log2。
    分子种“同向个数”不是检验，类别总量才是一个可检验的特征。分子种少于 ``min_species`` 的类别不计算。

质量作用比（``feature_type = mass_action_ratio``）
    近平衡反应的产物 / 底物比，如乳酸 / 丙酮酸（LDH，胞质 NADH / NAD⁺）、G3P / DHAP（GPD1）、
    β-羟丁酸 / 乙酰乙酸（BHBDH，线粒体 NADH / NAD⁺；Williamson 1967）。按样本取 log2(分子) − log2(分母)。
    广靶代谢组是相对定量，不同化合物的峰面积不能直接比较，所以比值本身没有绝对意义；但两组之间比值的
    log2 倍数变化等于两者 log2 倍数变化之差，离子化效率在其中抵消，组间比较有效。比值也混合了胞质与线粒体，
    不能换算成游离 NAD⁺ / NADH。

组分按数据表自带的 HMDB / KEGG ID 匹配（``HMDB``、``cpd_ID``、``KEGG_Compound`` 列）；没有这些列时按名称查 ID。
一个组分匹配到多个特征（同一代谢物的不同离子或方法）时，取中位强度最高的一个，其余记在 ``definition`` 中。
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from omicsweave.core.idmap.metabolite_id_mapping import lookup_metabolite_ids

LIPID_CLASS_TOTAL = "lipid_class_total"
MASS_ACTION_RATIO = "mass_action_ratio"
UNITS = {
    LIPID_CLASS_TOTAL: "log2 summed normalized peak area",
    MASS_ACTION_RATIO: "log2 ratio of normalized peak areas",
}
_LIPID_NAME = re.compile(r"^(?P<cls>[A-Za-z][A-Za-z0-9]*)\((?P<ether>[OP]-)?[dtm]?\d{1,2}:\d{1,2}")
_VENDOR_ID_COLUMNS = ("HMDB", "hmdb", "hmdb_id", "cpd_ID", "KEGG_Compound", "kegg_compound", "KEGG")


@dataclass(frozen=True)
class Component:
    name: str
    ids: frozenset[str]


@dataclass(frozen=True)
class RatioDefinition:
    id: str
    name: str
    numerator: Component
    denominator: Component
    reaction: str
    readout: str


def _component(name: str, *ids: str) -> Component:
    return Component(name, frozenset(ids))


# LC-MS 不分 D / L 与外消旋体；厂商库在不同批次中对同一个峰会给出外消旋或未指定构型的 ID
# （例如峰 MADN1568 “Glycerol 3-Phosphate”在一个批次为 HMDB0000126 / C00093，另一个批次为 HMDB0252849 / C03189）
LACTATE = _component("lactate", "HMDB0000190", "C00186", "HMDB0001311", "C00256", "HMDB0144295", "C01432")
PYRUVATE = _component("pyruvate", "HMDB0000243", "C00022")
G3P = _component("glycerol 3-phosphate", "HMDB0000126", "C00093", "HMDB0252849", "C03189")
DHAP = _component("dihydroxyacetone phosphate", "HMDB0001473", "C00111")
ALANINE = _component("alanine", "HMDB0000161", "C00041")
GLUTAMATE = _component("glutamate", "HMDB0000148", "C00025", "C00302")
AKG = _component("alpha-ketoglutarate", "HMDB0000208", "C00026")
BHB = _component("3-hydroxybutyrate", "HMDB0000011", "HMDB0000357", "C01089")
ACETOACETATE = _component("acetoacetate", "HMDB0000060", "C00164")
# 同一个苹果酸峰（MADN0427）在不同批次的厂商注释中写作 L- 或 D-
MALATE = _component("malate", "HMDB0000156", "HMDB0031518", "C00149", "C00497")
OXALOACETATE = _component("oxaloacetate", "HMDB0000223", "C00036")

RATIOS: tuple[RatioDefinition, ...] = (
    RatioDefinition("lactate_pyruvate", "lactate / pyruvate", LACTATE, PYRUVATE, "LDH_c",
                    "cytosolic NADH / NAD+ (LDH near equilibrium)"),
    RatioDefinition("g3p_dhap", "glycerol 3-phosphate / DHAP", G3P, DHAP, "GPS_c",
                    "cytosolic NADH / NAD+ (GPD1 near equilibrium)"),
    RatioDefinition("bhb_acetoacetate", "3-hydroxybutyrate / acetoacetate", BHB, ACETOACETATE, "BHBDH_m",
                    "mitochondrial NADH / NAD+ (BHBDH near equilibrium)"),
    RatioDefinition("malate_oxaloacetate", "malate / oxaloacetate", MALATE, OXALOACETATE, "MDH_c",
                    "NADH / NAD+ via MDH (cytosol and mitochondria mixed)"),
    RatioDefinition("alanine_pyruvate", "alanine / pyruvate", ALANINE, PYRUVATE, "GPT_c",
                    "ALT mass action, with glutamate / alpha-ketoglutarate"),
    RatioDefinition("glutamate_akg", "glutamate / alpha-ketoglutarate", GLUTAMATE, AKG, "GPT_c",
                    "transaminase / GLUD mass action, with alanine / pyruvate"),
)


@dataclass
class DerivedFeatures:
    """派生特征：``matrix`` 行为派生特征、列为样本（log2）；``features`` 为定义（与 matrix 行对齐）。"""
    matrix: pd.DataFrame
    features: pd.DataFrame
    notes: list[str] = field(default_factory=list)

    @property
    def empty(self) -> bool:
        return self.matrix.empty


FEATURE_COLUMNS = ["feature_type", "feature_name", "definition", "n_members", "members", "related_reaction",
                   "readout", "units"]


def lipid_class(name: object) -> str | None:
    """从 LIPID MAPS 简写解析脂质类别；醚脂记为 ``PE-O`` / ``PE-P``。不是脂质分子种时返回 None。"""
    match = _LIPID_NAME.match(str(name).strip()) if name is not None else None
    if not match:
        return None
    return match["cls"] + (f"-{match['ether'][0]}" if match["ether"] else "")


def lipid_class_totals(values: pd.DataFrame, names: pd.Series, *, min_species: int = 3) -> DerivedFeatures:
    """``values``：行为特征、列为样本的 log2 值；``names``：与行对齐的名称。"""
    classes = names.reindex(values.index).map(lipid_class)
    rows, definitions = [], []
    for cls, members in classes.dropna().groupby(classes.dropna(), sort=True):
        ids = list(members.index)
        if len(ids) < min_species:
            continue
        linear = np.power(2.0, values.loc[ids].astype(float))
        rows.append(pd.Series(np.log2(linear.sum(axis=0, min_count=1)), name=f"class_total:{cls}"))
        definitions.append({
            "feature_type": LIPID_CLASS_TOTAL, "feature_name": f"{cls} class total",
            "definition": f"log2 of the summed normalized peak areas of {len(ids)} {cls} species",
            "n_members": len(ids), "members": ";".join(map(str, ids)), "related_reaction": "",
            "readout": f"{cls} pool size", "units": UNITS[LIPID_CLASS_TOTAL],
        })
    return _assemble(rows, definitions, values.columns)


def _normalize_id(value: object) -> str | None:
    text = str(value).strip().upper()
    hmdb = re.fullmatch(r"HMDB0*(\d+)", text)
    if hmdb:
        return f"HMDB{int(hmdb.group(1)):07d}"
    kegg = re.fullmatch(r"C\d{5}", text)
    return text if kegg else None


def feature_database_ids(table: pd.DataFrame, name_column: str | None) -> pd.Series:
    """每个特征的 HMDB / KEGG ID 集合：优先用数据表自带的 ID 列，没有时按名称查 ID。"""
    columns = [c for c in _VENDOR_ID_COLUMNS if c in table.columns]
    out = {}
    for index, row in table.iterrows():
        ids = set()
        for column in columns:
            for piece in re.split(r"[;,|\s]+", str(row[column]) if pd.notna(row[column]) else ""):
                if (normalized := _normalize_id(piece)):
                    ids.add(normalized)
        if not ids and name_column and name_column in table.columns:
            found = lookup_metabolite_ids(str(row[name_column])) or {}
            ids = {n for key in ("hmdb_id", "kegg_compound") if (n := _normalize_id(found.get(key) or ""))}
        out[index] = frozenset(ids)
    return pd.Series(out, dtype=object)


def _pick(component: Component, ids: pd.Series, values: pd.DataFrame, names: pd.Series) -> tuple[object, list]:
    hits = [index for index, feature_ids in ids.items() if feature_ids & component.ids]
    if not hits:
        return None, []
    medians = values.loc[hits].median(axis=1)
    chosen = medians.idxmax()
    others = [h for h in hits if h != chosen]
    return chosen, [f"{names.get(h, h)} ({h})" for h in others]


def mass_action_ratios(values: pd.DataFrame, ids: pd.Series, names: pd.Series,
                       ratios: tuple[RatioDefinition, ...] = RATIOS) -> DerivedFeatures:
    """``values``：行为特征、列为样本的 log2 值；``ids``：每个特征的 HMDB / KEGG ID 集合（与行对齐）。"""
    ids = ids.reindex(values.index).apply(lambda v: v if isinstance(v, frozenset) else frozenset())
    names = names.reindex(values.index)
    rows, definitions, notes = [], [], []
    for ratio in ratios:
        numerator, num_others = _pick(ratio.numerator, ids, values, names)
        denominator, den_others = _pick(ratio.denominator, ids, values, names)
        missing = [c.name for c, hit in ((ratio.numerator, numerator), (ratio.denominator, denominator))
                   if hit is None]
        if missing:
            notes.append(f"{ratio.name}：未测到 {'、'.join(missing)}，未计算")
            continue
        series = values.loc[numerator].astype(float) - values.loc[denominator].astype(float)
        rows.append(series.rename(f"ratio:{ratio.id}"))
        extra = "; ".join(filter(None, [
            f"other {ratio.numerator.name} features not used: {', '.join(num_others)}" if num_others else "",
            f"other {ratio.denominator.name} features not used: {', '.join(den_others)}" if den_others else ""]))
        definitions.append({
            "feature_type": MASS_ACTION_RATIO, "feature_name": ratio.name,
            "definition": (f"log2({names.get(numerator, numerator)}) - log2({names.get(denominator, denominator)})"
                           + (f"; {extra}" if extra else "")),
            "n_members": 2, "members": f"{numerator};{denominator}", "related_reaction": ratio.reaction,
            "readout": ratio.readout, "units": UNITS[MASS_ACTION_RATIO],
        })
    derived = _assemble(rows, definitions, values.columns)
    derived.notes = notes
    return derived


def _assemble(rows: list[pd.Series], definitions: list[dict], samples: pd.Index) -> DerivedFeatures:
    if not rows:
        return DerivedFeatures(pd.DataFrame(columns=samples, dtype=float),
                               pd.DataFrame(columns=FEATURE_COLUMNS))
    matrix = pd.DataFrame(rows).reindex(columns=samples)
    matrix.index.name = "feature_id"
    features = pd.DataFrame(definitions, index=matrix.index)[FEATURE_COLUMNS]
    return DerivedFeatures(matrix, features)


def derive_features(processed: pd.DataFrame, samples: list[str], *, feature_id_column: str,
                    name_column: str | None, min_species: int = 3) -> DerivedFeatures:
    """从一个组织的预处理结果（行为特征，含注释列与样本列）计算全部派生特征。"""
    index = pd.Index(processed[feature_id_column].astype(str), name="feature_id")
    values = processed[samples].set_axis(index)
    names = (processed[name_column].astype(str) if name_column else processed[feature_id_column].astype(str))
    names = names.set_axis(index)
    lipids = lipid_class_totals(values, names, min_species=min_species)
    ratios = mass_action_ratios(values, feature_database_ids(processed.set_axis(index), name_column), names)
    parts = [d for d in (lipids, ratios) if not d.empty]
    if not parts:
        return DerivedFeatures(pd.DataFrame(columns=samples, dtype=float), pd.DataFrame(columns=FEATURE_COLUMNS),
                               ratios.notes)
    return DerivedFeatures(pd.concat([d.matrix for d in parts]), pd.concat([d.features for d in parts]),
                           ratios.notes)
