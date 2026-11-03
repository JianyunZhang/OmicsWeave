#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
通路收敛（``pathway_convergence``）：同一组织、同一比较下，多个组学的通路富集结果在通路层面是否一致。

读取源结果包中各组学的 ``enrichment.parquet``，通路按名称匹配（去掉库名前缀、小写、合并空白），比较按
``contrast_id`` 区分，不跨比较合并。前提：至少两个组学有富集结果。

- 证据：每条富集结果一行，概率取其 padj，在（组织, 比较, 组学）内按概率排序得到秩。GSEA 的置换 p 值 / FDR
  为 0 时已在生成富集结果表时改为 1/(置换次数+1)（``flags`` 记 ``pvalue_floored`` / ``padj_floored``）；
  旧结果包里仍为 0 的会使融合 p 值趋近 0，此时在运行记录中提示重新运行；
- 融合：同一组学内的多条结果（不同方法、方向或库）按 Bonferroni 合并，至少两个组学有该通路时组学间用
  Fisher 法合并，各（组织, 比较）内 BH 校正（与旧流程相同，见
  :func:`omicsweave.integration.methods.pathway_convergence.fuse_pathway_evidence`）。另列出各组学中该通路
  最显著那条结果的效应方向，以及方向是否一致；
- 共享：各组学显著的通路（ORA 的 padj < ``ora_padj_cutoff``，GSEA 的 FDR < ``gsea_padj_cutoff``）两两取交集。

输出（``integration/pathway_convergence/tables/``）：``pathway_evidence.parquet``、``pathway_fusion.parquet``、
``pathway_overlaps.parquet``。
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from omicsweave.contract.tables import read_enrichment
from omicsweave.core.model import IntegrationMode
from omicsweave.core.run import RunContext, RunError

from ..integrate import IntegrationRun, mode_dir
from ..methods.pathway_convergence import fuse_pathway_evidence, normalize_term

MODE = IntegrationMode.PATHWAY_CONVERGENCE
# 与旧流程相同：ORA 用 padj < 0.05；GSEA 按惯例用 FDR < 0.25
DEFAULTS: dict[str, Any] = {"ora_padj_cutoff": 0.05, "gsea_padj_cutoff": 0.25}


def evidence_table(enrichment: pd.DataFrame) -> pd.DataFrame:
    """统一富集结果表 → 证据长表（列与 ``fuse_pathway_evidence`` 的输入一致，另保留来源信息）。"""
    working = enrichment.assign(
        comparison=enrichment["contrast_id"].astype(str),
        pathway_term=enrichment["term_name"].map(normalize_term),
        probability=pd.to_numeric(enrichment["padj"], errors="coerce").clip(0.0, 1.0),
    )
    working = working[working["pathway_term"] != ""]
    parts = []
    for _, group in working.groupby(["tissue", "comparison", "omics"], sort=True):
        group = group.sort_values("probability", na_position="last")
        parts.append(group.assign(rank=np.arange(1, len(group) + 1),
                                  rank_fraction=np.arange(1, len(group) + 1) / (len(group) + 1.0)))
    columns = ["tissue", "comparison", "contrast_role", "omics", "pathway_term", "probability", "rank",
               "rank_fraction", "method", "query", "gene_set_library", "term_name", "effect_direction"]
    if not parts:
        return pd.DataFrame(columns=columns)
    return pd.concat(parts, ignore_index=True).assign(source_probability_column="padj")[
        [*columns, "source_probability_column"]]


def _agreement(directions: list[str]) -> str:
    signed = {d for d in directions if d in ("up", "down")}
    if len(signed) == 2:
        return "conflicting"
    if len(signed) == 1 and sum(d in ("up", "down") for d in directions) >= 2:
        return "consistent"
    return "undetermined"


def fusion_table(evidence: pd.DataFrame) -> pd.DataFrame:
    """融合结果，并补上比较角色、各组学的效应方向与方向一致性。"""
    fused = fuse_pathway_evidence(evidence)
    if fused.empty:
        return fused.assign(contrast_role=pd.Series(dtype=str), omics_directions=pd.Series(dtype=str),
                            direction_agreement=pd.Series(dtype=str))
    best = (evidence.sort_values("probability", na_position="last")
            .drop_duplicates(["tissue", "comparison", "pathway_term", "omics"]))
    directions: dict[tuple, dict[str, str]] = {}
    for row in best.itertuples(index=False):
        directions.setdefault((row.tissue, row.comparison, row.pathway_term), {})[row.omics] = row.effect_direction
    roles = evidence.drop_duplicates("comparison").set_index("comparison")["contrast_role"]
    keys = list(zip(fused["tissue"], fused["comparison"], fused["pathway_term"]))
    fused["contrast_role"] = fused["comparison"].map(roles)
    fused["omics_directions"] = [";".join(f"{o}:{d}" for o, d in sorted(directions[k].items())) for k in keys]
    fused["direction_agreement"] = [_agreement(list(directions[k].values())) for k in keys]
    return fused


def overlap_table(enrichment: pd.DataFrame, ora_cutoff: float, gsea_cutoff: float) -> pd.DataFrame:
    """各组学显著通路的两两交集。"""
    cutoff = np.where(enrichment["method"] == "gsea_prerank", gsea_cutoff, ora_cutoff)
    significant = enrichment[pd.to_numeric(enrichment["padj"], errors="coerce") < cutoff]
    rows = []
    for (tissue, contrast, role), group in significant.groupby(["tissue", "contrast_id", "contrast_role"], sort=True):
        sets = {omics: set(g["term_name"].map(normalize_term)) - {""} for omics, g in group.groupby("omics")}
        names = sorted(sets)
        for i, a in enumerate(names):
            for b in names[i + 1:]:
                rows += [{"tissue": tissue, "comparison": contrast, "contrast_role": role, "omics_A": a,
                          "omics_B": b, "pathway_term": term} for term in sorted(sets[a] & sets[b])]
    return pd.DataFrame(rows, columns=["tissue", "comparison", "contrast_role", "omics_A", "omics_B",
                                       "pathway_term"])


def run(ctx: RunContext, source: Path, parameters: dict[str, Any]) -> IntegrationRun:
    params = {**DEFAULTS, **parameters}
    unknown = sorted(set(parameters) - set(DEFAULTS))
    if unknown:
        raise RunError(f"integration.parameters.{MODE.value} 中有无法识别的参数：{'、'.join(unknown)}"
                       f"（可用：{'、'.join(DEFAULTS)}）")
    paths = sorted(source.glob("omics/*/tables/enrichment.parquet"))
    if len(paths) < 2:
        found = "、".join(p.parent.parent.name for p in paths) or "没有"
        raise RunError(f"{MODE.value} 需要至少 2 个组学的富集结果；结果包 {source.name} 中只有：{found}。"
                       f"请先对启用的组学运行 omicsweave run（含 enrichment 步骤）")
    enrichment = pd.concat([read_enrichment(p) for p in paths], ignore_index=True)

    evidence = evidence_table(enrichment)
    fusion = fusion_table(evidence)
    overlaps = overlap_table(enrichment, float(params["ora_padj_cutoff"]), float(params["gsea_padj_cutoff"]))
    tables = mode_dir(ctx, MODE) / "tables"
    tables.mkdir(parents=True)
    evidence.to_parquet(tables / "pathway_evidence.parquet", index=False)
    fusion.to_parquet(tables / "pathway_fusion.parquet", index=False)
    overlaps.to_parquet(tables / "pathway_overlaps.parquet", index=False)

    record = IntegrationRun(MODE, inputs=paths, parameters={**params, "omics": sorted(enrichment["omics"].unique())})
    zero = (enrichment["method"] == "gsea_prerank") & (pd.to_numeric(enrichment["padj"], errors="coerce") == 0)
    if zero.any():
        record.notes.append(f"有 {int(zero.sum())} 条 GSEA 结果的 padj 为 0（结果包早于置换下限 1/(n+1) 的修正），"
                            "相关通路的融合 p 值会趋近 0；请重新运行 omicsweave run 的 enrichment 步骤")
    if fusion.empty:
        record.notes.append("没有在两个及以上组学中同时出现的通路（通路名不同、或富集结果过少）")
    return record
