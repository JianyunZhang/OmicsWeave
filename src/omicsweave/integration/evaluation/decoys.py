#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
假证据率（E3）：用诱饵映射估计“映射层与判定规则制造证据”的倾向（目标—诱饵策略，Elias 与 Gygi 2007）。

对每个可给出方向证据的（反应, 组织, 比较）和每个组学层，把该层匹配到的方向成员（直接酶 / 转运体，代谢物为底物与
产物）换成同一组织数据中丰度相近（同一丰度十分位）的随机特征，数目与角色不变，再用同一个证据构建器判定该层是否
“支持”。诱饵的支持率是在映射与真实生物学无关时也会出现的支持率；它与目标映射支持率之比估计支持证据中的假证据比例。

随机种子与重复次数登记在运行清单中（``integration.parameters.flux_anchored.decoys``）。只表示容量的反应（近平衡、
转运、交换）与通量估计为 0 的反应不参与。
"""
from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from ..mappings import match_metabolite_roles
from ..mappings.feature_resolver import measured_metabolites
from ..mappings.model_package import DecoyPackage
from ..methods.evidence_profile import ANCHORED_STATES

LAYERS = {"transcriptomics": ("tx_details", "n_direct_genes", "base_mean"),
          "proteomics": ("prot_details", "n_direct_enzymes", "abundance"),
          "metabolomics": ("metab_details", None, "abundance")}
NAME_COLUMNS = {"transcriptomics": ("gene", "gene_symbol", "canonical_symbol"), "proteomics": ("protein", "Protein"),
                "metabolomics": ("compound",)}
DETAIL_COLUMNS = ["replicate", "comparison", "tissue", "flux_id", "layer", "anchor_state", "n_members",
                  "target_status", "decoy_status"]


def _name_column(frame: pd.DataFrame, layer: str) -> str | None:
    return next((c for c in NAME_COLUMNS[layer] if c in frame.columns), None)


def _deciles(frame: pd.DataFrame, column: str) -> pd.Series:
    values = pd.to_numeric(frame[column], errors="coerce") if column in frame.columns else pd.Series(np.nan,
                                                                                                    index=frame.index)
    if values.notna().sum() < 10:
        return pd.Series(0, index=frame.index)
    return pd.qcut(values.rank(method="first"), 10, labels=False).fillna(-1).astype(int)


class _Pool:
    """一张数据表的诱饵候选：按丰度十分位分组（派生特征不作候选）。"""

    def __init__(self, frame: pd.DataFrame, layer: str):
        candidates = frame.index
        if "feature_type" in frame.columns:
            candidates = frame.index[frame["feature_type"].isna()]
        bins = _deciles(frame.loc[candidates], LAYERS[layer][2])
        self.bins = bins
        self.all = np.asarray(candidates, dtype=object)
        self.groups = {b: np.asarray(idx, dtype=object) for b, idx in bins.groupby(bins).groups.items()}

    def sample(self, targets: list[Any], exclude: set, rng: np.random.Generator) -> list[Any]:
        chosen: list[Any] = []
        for target in targets:
            group = self.groups.get(self.bins.get(target), self.all)
            for source in (group, self.all):
                pick = None
                for _ in range(50):
                    candidate = source[int(rng.integers(len(source)))] if len(source) else None
                    if candidate is not None and candidate not in exclude and candidate not in chosen:
                        pick = candidate
                        break
                if pick is not None:
                    chosen.append(pick)
                    break
        return chosen


def _metabolite_members(builder, flux_id: str, frame: pd.DataFrame) -> dict[str, list[Any]]:
    roles = builder.package.metabolite_roles(flux_id)
    if not isinstance(roles, dict):
        return {"substrates": [], "products": []}
    measured = measured_metabolites(frame)
    matched = match_metabolite_roles(measured, roles)
    return {role: list(matched.indices_by_role[role]) for role in ("substrates", "products")}


def _direct_members(layer: str, builder, flux_row: pd.Series, frame: pd.DataFrame, details: dict) -> list[Any]:
    """真实映射在该层匹配到的方向成员（数据行索引）。"""
    from ..mappings import match_gene_features, match_protein_features

    flux_id = str(flux_row.get("flux_id", ""))
    enzymes = builder.package.enzyme_mapping(flux_id)
    if layer == "proteomics":
        members = list(enzymes.get("direct_enzymes", [])) or list(enzymes.get("transporters", []))
        return list(match_protein_features(frame, members).indices) if members else []
    genes = builder.package.gene_mapping(flux_id)
    roles = genes.get("gene_roles", {}) or {}
    direct = list(enzymes.get("direct_enzymes", [])) or list(roles.get("catalytic", []))
    if not direct:
        direct = list(enzymes.get("transporters", [])) + list(roles.get("transporter", []))
    return list(match_gene_features(frame, direct).indices) if direct else []


def _evaluate(builder, layer: str, flux_row: pd.Series, frame: pd.DataFrame, comparison: str) -> str:
    if layer == "transcriptomics":
        details = builder.transcriptomics_evidence(flux_row, frame)
    elif layer == "proteomics":
        details = builder.proteomics_evidence(flux_row, frame, comparison)
    else:
        details = builder.metabolomics_evidence(flux_row, frame, comparison)
    record = {"tissue": flux_row.get("tissue", ""), "comparison": comparison, LAYERS[layer][0]: details}
    return str(builder.evidence_profile_scorer._build_layer_record(record, layer)["assessment_status"])


def decoy_evaluation(builder, mfa: dict[str, pd.DataFrame], tx: dict[str, dict[str, pd.DataFrame]],
                     prot: dict[str, pd.DataFrame], metab: dict[str, pd.DataFrame], evidence: pd.DataFrame,
                     *, replicates: int, seed: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    """返回 (逐条诱饵结果, 按层与锚点状态汇总)。``evidence`` 为真实映射的证据长表。"""
    base = builder.package
    target = evidence.set_index(["comparison", "tissue", "flux_id", "layer"])
    rows = []
    pools: dict[tuple[str, str, int], _Pool] = {}
    try:
        for replicate in range(replicates):
            rng = np.random.default_rng(seed + replicate)
            for comparison, mfa_df in mfa.items():
                for _, flux_row in mfa_df.iterrows():
                    flux_id, tissue = str(flux_row.get("flux_id", "")), str(flux_row.get("tissue", ""))
                    fc = flux_row.get("flux_log2fc", 0.0)
                    if not base.semantics(flux_id).directional or pd.isna(fc) or fc == 0:
                        continue
                    frames = {"transcriptomics": tx.get(comparison, {}).get(tissue),
                              "proteomics": prot.get(tissue), "metabolomics": metab.get(tissue)}
                    for layer, frame in frames.items():
                        if frame is None or frame.empty or (comparison, tissue, flux_id, layer) not in target.index:
                            continue
                        name_column = _name_column(frame, layer)
                        if name_column is None:
                            continue
                        builder.package = base
                        pool = pools.setdefault((layer, tissue, id(frame)), _Pool(frame, layer))
                        if layer == "metabolomics":
                            roles = _metabolite_members(builder, flux_id, frame)
                            members = roles["substrates"] + roles["products"]
                            if not members:
                                continue
                            subs = pool.sample(roles["substrates"], set(members), rng)
                            prods = pool.sample(roles["products"], set(members) | set(subs), rng)
                            decoy = DecoyPackage(base, metabolites={flux_id: {
                                "substrates": frame.loc[subs, name_column].astype(str).tolist(),
                                "products": frame.loc[prods, name_column].astype(str).tolist()}})
                        else:
                            members = _direct_members(layer, builder, flux_row, frame, {})
                            if not members:
                                continue
                            picked = pool.sample(members, set(members), rng)
                            names = frame.loc[picked, name_column].astype(str).tolist()
                            decoy = DecoyPackage(base, **({"genes": {flux_id: names}} if layer == "transcriptomics"
                                                          else {"proteins": {flux_id: names}}))
                        builder.package = decoy
                        status = _evaluate(builder, layer, flux_row, frame, comparison)
                        key = (comparison, tissue, flux_id, layer)
                        rows.append((replicate, comparison, tissue, flux_id, layer,
                                     str(target.loc[key, "anchor_state"]), len(members),
                                     str(target.loc[key, "assessment_status"]), status))
    finally:
        builder.package = base
    detail = pd.DataFrame(rows, columns=DETAIL_COLUMNS)
    return detail, summarize(detail)


def summarize(detail: pd.DataFrame) -> pd.DataFrame:
    """按层（及锚点是否成立）汇总：目标支持率、诱饵支持率与两者之比（估计的假证据比例，上限 1）。"""
    columns = ["layer", "subset", "n_cells", "n_replicates", "target_support_rate", "decoy_support_rate",
               "decoy_conflict_rate", "false_evidence_ratio"]
    if detail.empty:
        return pd.DataFrame(columns=columns)
    detail = detail.assign(anchored=detail["anchor_state"].isin(ANCHORED_STATES))
    out = []
    for layer, part in detail.groupby("layer", sort=True):
        for subset, rows in (("all_directional", part), ("anchored", part[part["anchored"]])):
            if rows.empty:
                continue
            first = rows[rows["replicate"] == rows["replicate"].min()]
            target_rate = float((first["target_status"] == "support").mean())
            decoy_rate = float((rows["decoy_status"] == "support").mean())
            conflict_rate = float(rows["decoy_status"].isin(["opposite", "mixed"]).mean())
            ratio = min(decoy_rate / target_rate, 1.0) if target_rate > 0 else float("nan")
            out.append((layer, subset, len(first), int(rows["replicate"].nunique()), target_rate, decoy_rate,
                        conflict_rate, ratio))
    return pd.DataFrame(out, columns=columns)
