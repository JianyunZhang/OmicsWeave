#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
通路证据矩阵的视图（图 A 与导出表共用，G14）：一个（面板, 组织, 比较）的每一行（反应或代谢物节点）× 每一列
（各块、各来源、各示踪剂）的格子，全部由 ``dossier_members`` 汇总得到。图只画这张表，导出表就是这张表与
成员行，所以同一个（反应, 层, 成员）在图、表与档案页中的数值相同。

列分为五块：通量（各示踪剂与合并估计，每个有通量数据的来源一组）、酶（P 蛋白、R RNA，每个来源一列）、
调控（P、R 与 M 别构效应物）、池（M：代谢物节点的测得值；反应行为一步产物 / 底物比）、评估（四条轴 F、C、R、M
的状态与方向，以及结论等级；原样取自 S3 的 axis_assessment 与 reaction_assessment）。

格子的汇总只用于展示，不是判定（判定见 S3 的四轴）：

- 单个成员：显示它自己的状态与 log2FC；
- 多个成员（``display_state`` 为 ``summary``）：酶取同工酶按丰度占比加权的 log2FC，没有占比（如复合体）时取中位数；
  共用亚基（如 DLD）不计入酶的格子（值、个数与描边），只在成员行与档案页中列出；调控取“对酶活性的指向”——log2FC 乘以作用方向（激活 +1、抑制 −1）的中位数，方向未知的调控因子
  只计入个数；``n_changed / n_quantified`` 为通过 FDR 的成员数 / 测到的成员数，描边表示至少一个成员通过 FDR；
- 没有测到的格子取最能说明原因的状态：全部 not_measured、QC 阻断、未检出，或知识表没有成员。

``certainty``（点的大小）：1 − lfsr；没有时为 Φ(|log2FC| / SE)；再没有时为 1 − p / 2，即方向判断正确的近似概率。
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy.stats import norm

from omicsweave.contract.export import write_excel_tsv
from ..mappings.feature_resolver import metabolite_members
from ..mappings.reaction_members import ACTIVITY_MODES
from .panels import Panel, species_labels
from ..mappings.network_annotation import split_metabolite

MEASURED = ("changed", "nominal", "no_evidence", "equivalent")
FLUX_STATES = ("not_identifiable", "insufficient_n", "not_infused")
VIEW_COLUMNS = ["panel_id", "tissue", "contrast_id", "row", "element_type", "element_id", "label", "segment",
                "block", "column", "source", "layer", "tracer", "expected_tracer", "value", "certainty",
                "display_state", "n_changed", "n_quantified", "ring", "members", "text"]
AXIS_LABELS = {"flux": "F", "capacity": "C", "regulation": "R", "pool": "M"}
LAYER_CODE = {"proteomics": "P", "transcriptomics": "R", "metabolomics": "M"}
MEMBER_EXPORT = ["source", "flux_id", "layer", "role", "member_id", "feature_id", "feature_name", "tracer",
                 "member_relation", "regulator_effect", "regulator_mode", "log2fc", "log2fc_unshrunk", "se", "pvalue",
                 "padj", "lfsr", "n_numerator", "n_denominator", "abundance_share", "cell_state", "state_reason",
                 "selected", "identifiability", "flags", "implied_direction"]


def certainty(row: pd.Series | dict) -> float:
    lfsr, fc, se, p = (row.get(k) for k in ("lfsr", "log2fc", "se", "pvalue"))
    if pd.notna(lfsr):
        return float(1 - lfsr)
    if pd.notna(fc) and pd.notna(se) and se > 0:
        return float(norm.cdf(abs(fc) / se))
    if pd.notna(p):
        return float(1 - p / 2)
    return np.nan


def _activity_sign(row: pd.Series) -> int:
    if row.get("regulator_mode") not in ACTIVITY_MODES:
        return 0
    return {"activates": 1, "inhibits": -1}.get(row.get("regulator_effect"), 0)


def _enzyme_value(rows: pd.DataFrame) -> float:
    rows = rows[rows["member_relation"].fillna("") != "shared_subunit"]
    fc = pd.to_numeric(rows["log2fc"], errors="coerce")
    share = pd.to_numeric(rows.get("abundance_share"), errors="coerce") if "abundance_share" in rows else None
    if share is not None and share.notna().any():
        usable = share.notna() & fc.notna()
        return float((fc[usable] * share[usable]).sum() / share[usable].sum()) if usable.any() else np.nan
    return float(fc.median()) if fc.notna().any() else np.nan


def _regulator_value(rows: pd.DataFrame) -> float:
    signs = rows.apply(_activity_sign, axis=1)
    fc = pd.to_numeric(rows["log2fc"], errors="coerce")
    signed = (fc * signs)[signs != 0].dropna()
    return float(signed.median()) if not signed.empty else np.nan


def _plain_value(rows: pd.DataFrame) -> float:
    fc = pd.to_numeric(rows["log2fc"], errors="coerce").dropna()
    return float(fc.median()) if not fc.empty else np.nan


def summarize(rows: pd.DataFrame, value_fn=_plain_value, single_value_fn=None) -> dict[str, Any]:
    """一组成员行 → 一个格子。"""
    empty = {"value": np.nan, "certainty": np.nan, "n_changed": 0, "n_quantified": 0, "ring": "", "members": ""}
    if rows.empty or rows["cell_state"].eq("no_knowledge").all():
        return {**empty, "display_state": "no_knowledge"}
    rows = rows[rows["cell_state"] != "no_knowledge"]
    quantified = rows[rows["cell_state"].isin(MEASURED)]
    members = ";".join(f"{r.member_id or r.feature_name}({r.cell_state}"
                       + (f", {r.log2fc:+.2f}" if pd.notna(r.log2fc) else "") + ")" for r in rows.itertuples())
    if quantified.empty:
        states = list(dict.fromkeys(rows["cell_state"]))
        for state in ("qc_blocked", *FLUX_STATES, "not_detected", "not_measured"):
            if state in states:
                return {**empty, "display_state": state, "members": members}
        return {**empty, "display_state": states[0], "members": members}
    n_changed = int(quantified["cell_state"].eq("changed").sum())
    out = {"n_changed": n_changed, "n_quantified": len(quantified), "members": members}
    if len(quantified) == 1 and len(rows) == 1:
        row = quantified.iloc[0]
        value = (single_value_fn or (lambda r: r["log2fc"]))(row)
        ring = "selected" if row["cell_state"] == "changed" and bool(row["selected"]) else \
            ("changed" if row["cell_state"] == "changed" else "")
        return {**out, "display_state": row["cell_state"], "value": value, "certainty": certainty(row), "ring": ring}
    values = quantified.apply(certainty, axis=1)
    return {**out, "display_state": "summary", "value": value_fn(quantified),
            "certainty": float(values.median()) if values.notna().any() else np.nan,
            "ring": "changed" if n_changed else ""}


def _signed_single(row: pd.Series) -> float:
    sign = _activity_sign(row)
    return float(row["log2fc"]) * sign if sign else np.nan


def _metabolite_key(name: str) -> str:
    members = metabolite_members([name])
    return members[0].key if members else f"name:{name.lower()}"


def species_key(species: str, keys: dict[str, str]) -> str:
    """网络物种（如 ``PYR_m``）→ 代谢物的标准键；``keys`` 为调用方持有的缓存。区室不区分（数据中同一个特征）。"""
    base = split_metabolite(species)[0]
    return keys.setdefault(base, _metabolite_key(species_labels().get(base, base)))


def node_members(rows: pd.DataFrame, species: str, keys: dict[str, str]) -> pd.DataFrame:
    """成员行中属于某个代谢物节点的底物 / 产物行：有测得特征时按特征去重，否则按成员去重（保留未检出等状态）。"""
    key = species_key(species, keys)
    part = rows[(rows["layer"] == "metabolomics") & rows["role"].isin(("substrate", "product"))]
    part = part[part["member_id"].map(lambda m: keys.setdefault(f"member:{m}", _metabolite_key(str(m)))) == key]
    measured = part[part["feature_id"].notna() & part["feature_id"].astype(str).ne("")]
    return measured.drop_duplicates("feature_id") if not measured.empty else part.drop_duplicates("member_id")


def matrix_view(panel: Panel, elements: pd.DataFrame, members: pd.DataFrame, tissue: str, contrast: str,
                sources: list[str], tracers: Iterable[str] = (),
                assessment: tuple[pd.DataFrame, pd.DataFrame] | None = None) -> pd.DataFrame:
    """``members``：dossier_members（全部来源）；``sources``：来源标签的顺序（本实验在前）；``assessment``：
    （axis_assessment, reaction_assessment），有时加上评估块（四条轴与结论等级，原样取自评估表）。"""
    cell_rows = members[(members["tissue"] == tissue) & (members["contrast_id"] == contrast)
                        & members["flux_id"].isin(panel.reactions)]
    records: list[dict[str, Any]] = []
    layers_by_source = {s: set(cell_rows.loc[cell_rows["source"] == s, "layer"]) for s in sources}
    mfa_sources = [s for s in sources if (cell_rows.loc[(cell_rows["source"] == s) & (cell_rows["layer"] == "mfa"),
                                                        "cell_state"] != "not_measured").any()]
    declared = list(dict.fromkeys(tracers))
    flux_columns = []
    for source in mfa_sources:
        present = cell_rows.loc[(cell_rows["source"] == source) & (cell_rows["layer"] == "mfa"), "tracer"]
        # 声明的示踪剂描述本实验（第一个来源）；其他来源只列它自己有数据的示踪剂
        lead = declared if source == sources[0] else []
        order = lead + sorted(set(present) - set(lead) - {"pooled"}) + ["pooled"]
        flux_columns += [(source, tracer) for tracer in order]
    label = (lambda tracer, source: tracer if len(mfa_sources) == 1 else f"{tracer}·{source}")
    keys = {}

    def add(element: pd.Series, block: str, column: str, source: str, layer: str, tracer: str, cell: dict,
            expected: bool = False) -> None:
        records.append({"panel_id": panel.id, "tissue": tissue, "contrast_id": contrast, "row": int(element["order"]),
                        "element_type": element["element_type"], "element_id": element["element_id"],
                        "label": element["label"], "segment": element["segment"], "block": block, "column": column,
                        "source": source, "layer": layer, "tracer": tracer, "expected_tracer": expected, **cell})

    for _, element in elements.iterrows():
        if element["element_type"] == "reaction":
            reaction = element["element_id"]
            rows = cell_rows[cell_rows["flux_id"] == reaction]
            for source, tracer in flux_columns:
                part = rows[(rows["source"] == source) & (rows["layer"] == "mfa") & (rows["tracer"] == tracer)]
                expected = tracer in panel.expected_tracers.get(reaction, [])
                add(element, "flux", label(tracer, source), source, "mfa", tracer, summarize(part), expected)
            for block, roles, fn, single in (("enzyme", ("direct_enzyme",), _enzyme_value, None),
                                             ("regulator", ("regulator",), _regulator_value, _signed_single)):
                for layer in ("proteomics", "transcriptomics"):
                    for source in sources:
                        if layer not in layers_by_source[source]:
                            continue
                        part = rows[(rows["source"] == source) & (rows["layer"] == layer)]
                        chosen = part[part["role"].isin(roles)]
                        if block == "enzyme":
                            if chosen[chosen["cell_state"] != "no_knowledge"].empty:
                                chosen = part[part["role"].isin(("carrier", "direct_enzyme"))]   # 转运反应：载体即执行者
                            # 共用亚基不参与方向判断，不计入格子（成员明细与档案页仍列出）
                            chosen = chosen[chosen["member_relation"].fillna("") != "shared_subunit"]
                        add(element, block, f"{LAYER_CODE[layer]}:{source}", source, layer, "",
                            summarize(chosen, fn, single))
            for source in sources:
                if "metabolomics" not in layers_by_source[source]:
                    continue
                part = rows[(rows["source"] == source) & (rows["layer"] == "metabolomics")]
                add(element, "regulator", f"M:{source}", source, "metabolomics", "",
                    summarize(part[part["role"] == "effector"], _regulator_value, _signed_single))
                ratio = rows[(rows["source"] == source) & (rows["layer"] == "derived") & (rows["role"] == "step_ratio")]
                add(element, "pool", f"M:{source}", source, "derived", "", summarize(ratio))
            if assessment is not None:
                _assessment_cells(add, element, assessment, reaction, tissue, contrast, sources[0])
        else:
            for source in sources:
                if "metabolomics" not in layers_by_source[source]:
                    continue
                part = node_members(cell_rows[cell_rows["source"] == source], element["element_id"], keys)
                add(element, "pool", f"M:{source}", source, "metabolomics", "", summarize(part))
    view = pd.DataFrame(records)
    for column in VIEW_COLUMNS:
        if column not in view.columns:
            view[column] = np.nan
    return view[VIEW_COLUMNS]


def _assessment_cells(add, element: pd.Series, assessment: tuple[pd.DataFrame, pd.DataFrame], reaction: str,
                      tissue: str, contrast: str, source: str) -> None:
    axes, reactions = assessment
    empty = {"value": np.nan, "certainty": np.nan, "n_changed": 0, "n_quantified": 0, "ring": "", "members": ""}
    here = axes[(axes["flux_id"] == reaction) & (axes["tissue"] == tissue) & (axes["contrast_id"] == contrast)]
    for axis, label in AXIS_LABELS.items():
        row = here[here["axis"] == axis]
        if row.empty:
            add(element, "assessment", label, source, "", "", {**empty, "display_state": "not_measured", "text": ""})
            continue
        row = row.iloc[0]
        add(element, "assessment", label, source, "", "",
            {**empty, "display_state": row["state"], "value": float(row["direction"]), "text": row["rule"],
             "members": row["basis"]})
    claim = reactions[(reactions["flux_id"] == reaction) & (reactions["tissue"] == tissue)
                      & (reactions["contrast_id"] == contrast)]
    if claim.empty:
        return
    claim = claim.iloc[0]
    level = claim["flux_claim"] + (f" ({int(claim['n_concordant_axes'])})" if claim["flux_claim"] in ("L1", "L2", "L3")
                                   else "")
    add(element, "assessment", "level", source, "", "",
        {**empty, "display_state": claim["flux_claim"], "value": float(claim["claim_direction"]), "text": level,
         "members": claim["phenomena"]})


def wide_table(view: pd.DataFrame) -> pd.DataFrame:
    """导出用的宽表：行与图相同，每个格子三列（值、状态、通过数 / 测到数）。"""
    index = ["row", "element_type", "element_id", "label", "segment"]
    out = view[index].drop_duplicates().set_index("row")
    for (block, column), cell in view.groupby(["block", "column"], sort=False):
        cell = cell.set_index("row")
        name = f"{block} {column}"
        if block == "assessment":                 # 评估块：方向（+1 / −1 / 0）、状态；结论等级另写文字
            out[f"{name} direction"] = cell["value"]
            out[f"{name} state"] = cell["display_state"]
            if column == "level":
                out[f"{name} text"] = cell["text"]
            continue
        out[f"{name} log2fc"] = cell["value"].round(4)
        out[f"{name} state"] = cell["display_state"]
        out[f"{name} k/n"] = cell.apply(lambda r: f"{int(r.n_changed)}/{int(r.n_quantified)}"
                                        if r.n_quantified else "", axis=1)
    return out.reset_index()


def member_table(members: pd.DataFrame, panel: Panel, tissue: str, contrast: str) -> pd.DataFrame:
    order = {r: i for i, r in enumerate(panel.reactions)}
    rows = members[(members["tissue"] == tissue) & (members["contrast_id"] == contrast)
                   & members["flux_id"].isin(panel.reactions)].copy()
    rows["_order"] = rows["flux_id"].map(order)
    return rows.sort_values(["_order", "source", "layer", "role", "member_id"], kind="stable")[MEMBER_EXPORT]


def export_matrix(view: pd.DataFrame, members: pd.DataFrame, stem: Path) -> list[Path]:
    """写出 ``<stem>.tsv``（宽表）、``<stem>_members.tsv``（成员行）与 ``<stem>.xlsx``（两张表各一页）。"""
    stem.parent.mkdir(parents=True, exist_ok=True)
    wide = wide_table(view)
    paths = [stem.with_suffix(".tsv"), stem.parent / f"{stem.name}_members.tsv"]
    write_excel_tsv(wide, paths[0])
    write_excel_tsv(members, paths[1])
    try:
        with pd.ExcelWriter(stem.with_suffix(".xlsx")) as writer:
            wide.to_excel(writer, sheet_name="matrix", index=False)
            members.to_excel(writer, sheet_name="members", index=False)
        paths.append(stem.with_suffix(".xlsx"))
    except ImportError:                     # 没有安装 openpyxl 时只写 TSV
        pass
    return paths
