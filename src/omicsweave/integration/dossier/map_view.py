#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
通路总图（图 E，G18）的视图：一个（组织, 比较）在总图布局上的每个元素一行。图只画这张表，导出表就是这张表。

元素有三种：

- ``reaction``：网络反应（布局中的边）。``draw_style`` 为 ``flux``（该组织有拟合通量）、``omics_only``（没有通量、
  有组学成员）或 ``not_assessed``（两者都没有）。``physio_sign`` 为生理净方向相对网络书写方向的符号（两组组均值同号时取
  该符号；变号时为 0 并记 ``sign_change``；没有通量时为 0）。F、M 两轴的方向换算到生理方向（乘以 ``physio_sign``）；
  C、R 两轴为酶容量与活性，不换号。边的颜色值 ``value`` 只在 F 轴为 suggestive 或 changed 时给出，取合并估计的
  log2FC 的绝对值乘以生理方向上的变化方向；
- ``metabolite``：网络物种节点。取值规则与图 A 的池格子相同（:func:`omicsweave.integration.dossier.view.node_members`
  选出的特征，:func:`omicsweave.integration.dossier.view.summarize` 汇总），所以两张图的数值一致；代谢组测全组织，
  同一代谢物的各区室节点共用一个值；
- ``aggregate``：不在网络中、只有组学的评估 ID（方框）。P、R 格的规则同图 A 的酶格子，M 格为底物与产物；
  成员（直接酶与载体）全部被画出的网络反应覆盖时不上图，记为 ``covered_by_steps``（见 :func:`coverage_table`）。

方向与状态取自 S3 的 ``axis_assessment`` 与 ``reaction_assessment``，不重新判定。
"""
from __future__ import annotations

from typing import Any, Iterable

import numpy as np
import pandas as pd

from omicsweave.core.idmap.reaction_pathway_mapping import REACTION_REGISTRY

from ..mappings.network_annotation import NetworkReaction, network_reactions
from ..mappings.network_layout import Layout, unplaced_ids
from ..mappings.reaction_aliases import reaction_aliases
from .panels import species_labels
from .view import _enzyme_value, _plain_value, node_members, summarize

AXES = ("flux", "capacity", "regulation", "pool")
POINTING = ("changed", "suggestive")
MAP_COLUMNS = ["tissue", "contrast_id", "element_type", "element_id", "region", "label", "draw_style", "reversible",
               "physio_sign", "sign_change", "flux_group", "kegg_pathway", "kegg_reaction", "ec_number",
               "F_state", "F_direction", "C_state", "C_direction", "R_state", "R_direction", "M_state", "M_direction",
               "flux_claim", "claim_direction", "n_concordant_axes", "rank", "key",
               "value", "certainty", "display_state", "ring", "n_changed", "n_quantified",
               "P_value", "P_state", "R_value", "R_layer_state", "Mx_value", "Mx_state", "members", "text"]


def _sign(value: Any) -> int:
    value = pd.to_numeric(value, errors="coerce")
    return 0 if pd.isna(value) or value == 0 else (1 if value > 0 else -1)


def physiological_sign(flux_rows: pd.DataFrame) -> tuple[int, bool]:
    """（生理净方向相对书写方向的符号, 两组间是否变号）。优先取合并估计的组均值，没有时取各示踪剂多数。"""
    usable = flux_rows.dropna(subset=["mean_numerator", "mean_denominator"])
    if usable.empty:
        return 0, False
    pooled = usable[usable["tracer"] == "pooled"]
    rows = pooled if not pooled.empty else usable
    signs = [(_sign(r.mean_numerator), _sign(r.mean_denominator)) for r in rows.itertuples()]
    if any(a and b and a != b for a, b in signs):
        return 0, True
    votes = [a or b for a, b in signs if a or b]
    if not votes:
        return 0, False
    total = sum(votes)
    return (0, False) if total == 0 else (1 if total > 0 else -1, False)


def _flux_group(flux_rows: pd.DataFrame) -> str:
    for flags in flux_rows["flags"].dropna():
        for flag in (flags if not isinstance(flags, str) else flags.split(";")):
            if str(flag).startswith("flux_group:"):
                return str(flag).split(":", 1)[1]
    return ""


def _registry(flux_id: str) -> dict[str, Any]:
    entry = REACTION_REGISTRY.get(flux_id) or {}
    return {"kegg_pathway": entry.get("kegg_pathway") or "", "kegg_reaction": entry.get("kegg_reaction") or "",
            "ec_number": entry.get("ec_number") or ""}


def _axis(axes: pd.DataFrame, flux_id: str) -> dict[str, tuple[str, int, float, str]]:
    here = axes[axes["flux_id"] == flux_id]
    out = {}
    for axis in AXES:
        row = here[here["axis"] == axis]
        out[axis] = ("not_measured", 0, np.nan, "") if row.empty else \
            (str(row["state"].iloc[0]), int(row["direction"].iloc[0] or 0), float(pd.to_numeric(row["effect"].iloc[0],
                                                                                                    errors="coerce")),
             str(row["rule"].iloc[0] or ""))
    return out


def _executors(rows: pd.DataFrame) -> pd.DataFrame:
    chosen = rows[rows["role"] == "direct_enzyme"]
    if chosen[chosen["cell_state"] != "no_knowledge"].empty:
        chosen = rows[rows["role"].isin(("carrier", "direct_enzyme"))]
    return chosen[chosen["member_relation"].fillna("") != "shared_subunit"]


def _layer_cells(rows: pd.DataFrame) -> dict[str, Any]:
    """聚合节点的 P、R、M 三格。"""
    prot = summarize(_executors(rows[rows["layer"] == "proteomics"]), _enzyme_value)
    rna = summarize(_executors(rows[rows["layer"] == "transcriptomics"]), _enzyme_value)
    met_rows = rows[(rows["layer"] == "metabolomics") & rows["role"].isin(("substrate", "product"))]
    measured = met_rows[met_rows["feature_id"].notna() & met_rows["feature_id"].astype(str).ne("")]
    met = summarize(measured.drop_duplicates("feature_id") if not measured.empty
                    else met_rows.drop_duplicates("member_id"), _plain_value)
    return {"P_value": prot["value"], "P_state": prot["display_state"], "R_value": rna["value"],
            "R_layer_state": rna["display_state"], "Mx_value": met["value"], "Mx_state": met["display_state"]}


def coverage_table(layout: Layout, members: pd.DataFrame, assessed: Iterable[str]) -> pd.DataFrame:
    """评估集中每个 ID 在总图上的去向：edge、aggregate、covered_by_steps、alias、exclude（原因）或 unplaced。"""
    executors = members[members["role"].isin(("direct_enzyme", "carrier"))
                        & members["layer"].isin(("proteomics", "transcriptomics"))
                        & (members["cell_state"] != "no_knowledge")]
    by_reaction = executors.groupby("flux_id")["member_id"].apply(
        lambda s: {m for m in s.dropna().astype(str) if m.strip()})
    drawn = set().union(*(by_reaction.get(r, set()) for r in layout.edges)) if layout.edges else set()
    aliases = reaction_aliases()
    rows = []
    assessed = sorted(set(map(str, assessed)))
    for flux_id in assessed:
        if flux_id in layout.edges:
            rows.append((flux_id, "edge", layout.edges[flux_id].region, ""))
        elif flux_id in layout.aggregates:
            own = by_reaction.get(flux_id, set())
            if own and own <= drawn:
                rows.append((flux_id, "covered_by_steps", layout.aggregates[flux_id].region,
                             f"{len(own)} 个直接酶 / 载体全部属于画出的网络反应"))
            else:
                rows.append((flux_id, "aggregate", layout.aggregates[flux_id].region,
                             f"{len(own & drawn)} / {len(own)} 个直接酶 / 载体与画出的网络反应共有" if own else ""))
        elif flux_id in aliases:
            rows.append((flux_id, "alias", "", "对应 " + "、".join(aliases[flux_id])))
        elif flux_id in layout.exclude:
            spec = layout.exclude[flux_id]
            rows.append((flux_id, f"exclude:{spec.get('reason')}", "", str(spec.get("of", ""))))
    for flux_id in unplaced_ids(layout, [i for i in assessed if i not in aliases]):
        rows.append((flux_id, "unplaced", "", "布局中没有这个评估 ID"))
    return pd.DataFrame(rows, columns=["flux_id", "placement", "region", "detail"])


def map_view(layout: Layout, members: pd.DataFrame, axes: pd.DataFrame, assessed: pd.DataFrame, tissue: str,
             contrast: str, source: str, coverage: pd.DataFrame | None = None, top_rank: int = 10,
             min_claim: str = "L2", reactions: dict[str, NetworkReaction] | None = None) -> pd.DataFrame:
    """``members``：dossier_members；``axes`` / ``assessed``：S3 的两张评估表；``source``：本实验的来源标签。"""
    reactions = network_reactions() if reactions is None else reactions
    order = {"L3": 4, "L2": 3, "L1": 2, "L0": 1, "LX": 0}
    rows = members[(members["tissue"] == tissue) & (members["contrast_id"] == contrast) & (members["source"] == source)]
    here_axes = axes[(axes["tissue"] == tissue) & (axes["contrast_id"] == contrast)]
    here_claims = assessed[(assessed["tissue"] == tissue) & (assessed["contrast_id"] == contrast)].set_index("flux_id")
    by_id = dict(tuple(rows.groupby("flux_id")))
    empty_rows = rows.iloc[0:0]
    labels = species_labels()
    records: list[dict[str, Any]] = []

    def claim_fields(flux_id: str, physio: int) -> dict[str, Any]:
        if flux_id not in here_claims.index:
            return {"flux_claim": "", "claim_direction": 0, "n_concordant_axes": 0, "rank": np.nan, "key": ""}
        claim = here_claims.loc[flux_id]
        level = str(claim["flux_claim"])
        direction = int(claim["claim_direction"] or 0)
        if level in ("L2", "L3") and physio:          # 通量等级：方向换算到生理方向；L1 为酶容量的方向，不换号
            direction *= physio
        rank = pd.to_numeric(claim.get("rank"), errors="coerce")
        key = "key" if order.get(level, 0) >= order[min_claim] else ("l1" if level == "L1" else "")
        if not key and pd.notna(rank) and rank <= top_rank:
            key = "top"
        return {"flux_claim": level, "claim_direction": direction,
                "n_concordant_axes": int(claim.get("n_concordant_axes") or 0), "rank": rank, "key": key}

    for flux_id, edge in layout.edges.items():
        reaction = reactions[flux_id]
        part = by_id.get(flux_id, empty_rows)
        flux_rows = part[part["layer"] == "mfa"]
        has_flux = bool((~flux_rows["cell_state"].isin(["not_measured", "not_infused"])).any())
        knowledge = part[(part["layer"] != "mfa") & (part["cell_state"] != "no_knowledge")]
        style = "flux" if has_flux else ("omics_only" if not knowledge.empty else "not_assessed")
        physio, change = physiological_sign(flux_rows) if has_flux else (0, False)
        ax = _axis(here_axes, flux_id)
        f_state, f_dir, f_effect, f_rule = ax["flux"]
        f_phys = f_dir * physio if physio else f_dir
        value = abs(f_effect) * f_phys if f_state in POINTING and pd.notna(f_effect) else np.nan
        m_dir = ax["pool"][1] * physio if physio else ax["pool"][1]
        records.append({"tissue": tissue, "contrast_id": contrast, "element_type": "reaction", "element_id": flux_id,
                        "region": edge.region, "label": flux_id, "draw_style": style, "reversible": reaction.reversible,
                        "physio_sign": physio, "sign_change": change, "flux_group": _flux_group(flux_rows),
                        **_registry(flux_id),
                        "F_state": f_state, "F_direction": f_phys, "C_state": ax["capacity"][0],
                        "C_direction": ax["capacity"][1], "R_state": ax["regulation"][0],
                        "R_direction": ax["regulation"][1], "M_state": ax["pool"][0], "M_direction": m_dir,
                        **claim_fields(flux_id, physio), "value": value, "certainty": np.nan,
                        "display_state": f_state, "ring": "", "n_changed": 0, "n_quantified": 0, "text": f_rule})

    keys: dict[str, str] = {}
    met_rows = rows[rows["layer"] == "metabolomics"]
    for species, node in layout.nodes.items():
        if node.kind != "metabolite":
            continue
        cell = summarize(node_members(met_rows, species, keys)) if not met_rows.empty else \
            {"value": np.nan, "certainty": np.nan, "n_changed": 0, "n_quantified": 0, "ring": "", "members": "",
             "display_state": "not_measured"}
        base = species.rsplit("_", 1)[0]
        records.append({"tissue": tissue, "contrast_id": contrast, "element_type": "metabolite", "element_id": species,
                        "region": "", "label": node.label or labels.get(base, base), "draw_style": "",
                        **{k: cell.get(k) for k in ("value", "certainty", "display_state", "ring", "n_changed",
                                                    "n_quantified", "members")}})

    coverage = coverage_table(layout, members, here_claims.index) if coverage is None else coverage
    hidden = set(coverage.loc[coverage["placement"] == "covered_by_steps", "flux_id"])
    for flux_id, agg in layout.aggregates.items():
        if flux_id in hidden:
            continue
        part = by_id.get(flux_id, empty_rows)
        cells = _layer_cells(part) if not part.empty else \
            {"P_value": np.nan, "P_state": "not_measured", "R_value": np.nan, "R_layer_state": "not_measured",
             "Mx_value": np.nan, "Mx_state": "not_measured"}
        ax = _axis(here_axes, flux_id)
        records.append({"tissue": tissue, "contrast_id": contrast, "element_type": "aggregate", "element_id": flux_id,
                        "region": agg.region, "label": agg.label or flux_id,
                        "draw_style": "aggregate" if not part.empty else "not_assessed", **_registry(flux_id),
                        "C_state": ax["capacity"][0], "C_direction": ax["capacity"][1], **claim_fields(flux_id, 0),
                        **cells, "display_state": ax["capacity"][0]})
    view = pd.DataFrame(records)
    for column in MAP_COLUMNS:
        if column not in view.columns:
            view[column] = np.nan
    return view[MAP_COLUMNS]


def region_summary(view: pd.DataFrame, layout: Layout, fusion: pd.DataFrame | None, tissue: str,
                   contrast: str) -> pd.DataFrame:
    """侧栏：每个区域的覆盖计数（不是评分）与通路融合结果（按区域声明的通路名匹配，只作探索排序）。"""
    out = []
    reactions = view[view["element_type"] == "reaction"]
    for rid, region in layout.regions.items():
        part = reactions[reactions["region"] == rid]
        if part.empty and rid != "other":
            continue
        row = {"region": rid, "label": region.label, "kegg": ", ".join(region.kegg), "n_drawn": len(part),
               "n_flux": int((part["draw_style"] == "flux").sum()),
               "n_identifiable": int((part["draw_style"].eq("flux")
                                      & ~part["F_state"].isin(["not_identifiable", "insufficient_n", "not_measured"])).sum()),
               "n_omics_only": int((part["draw_style"] == "omics_only").sum()),
               "n_L2": int(part["flux_claim"].isin(["L2", "L3"]).sum()), "n_L1": int((part["flux_claim"] == "L1").sum()),
               "fusion": ""}
        if fusion is not None and not fusion.empty and region.terms:
            hit = fusion[(fusion["tissue"] == tissue) & (fusion["comparison"] == contrast)
                         & fusion["pathway_term"].str.lower().isin([t.lower() for t in region.terms])]
            row["fusion"] = "; ".join(f"{r.pathway_term}: FDR {r.fused_fdr:.2g} ({r.omics_directions or '方向未定'})"
                                      for r in hit.itertuples())
        out.append(row)
    return pd.DataFrame(out)
