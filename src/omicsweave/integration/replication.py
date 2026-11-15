#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
跨实验、跨队列的复现（G6，``omicsweave replicate``）：``project.yaml`` 的 ``replication`` 中每个复现组列出估计对象
相同的若干比较（可来自其他项目），对每个特征、派生特征与通路判断效应是否复现。

- 对应：同层复现按特征对应（转录组、蛋白组与派生特征按 ``feature_id``，测得代谢物按规范化名称，因为各队列的厂商
  特征编号不同；通路按规范化的条目名称，只用预排序 GSEA）；``cross_layer`` 中声明的组学对（如一个队列的蛋白对另一个
  队列的 RNA）按基因符号对应，只用于两个来源的复现组；
- 检验：方向性的部分合取检验（Benjamini 与 Heller 2008）。每个来源的双侧 p 值换成两个方向的单侧 p 值，在每个方向上
  对“至少 u 个来源有该方向效应”做 Simes 型合并 ``min_{i≥u} (n−u+1)/(i−u+1)·p_(i)``，取两个方向中较小者并乘 2；
  在（复现组, 对应方式, 组学, 组织, 类别）内 BH 校正得到 ``replication_fdr``，小于 ``alpha`` 记为 ``replicated``；
- 合并：同层特征有标准误时给出随机效应合并估计与 I²（:func:`~omicsweave.core.statistics.meta_analysis.random_effects_meta`）；
- 基线：每组对应的全部特征中效应同号的比例，以及把一个来源的效应在特征间随机重排后的同号比例（应接近 50%）。

结果写入项目级结果包 ``<paths.results>/_project/<运行ID>/integration/replication/tables/``：``replication.parquet`` 与
``replication_baseline.parquet``；运行清单登记每个来源实际读取的运行。
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from omicsweave.contract.dataset import ContractError
from omicsweave.contract.manifest import RunManifest, file_record, new_run_id, write_manifest
from omicsweave.contract.query import open_run
from omicsweave.contract.tables import effect_standard_errors
from omicsweave.core.config import Project, load_project
from omicsweave.core.model import ReplicationSet, ReplicationSource
from omicsweave.core.statistics.linear_model import bh
from omicsweave.core.statistics.meta_analysis import random_effects_meta

PROJECT_RUNS = "_project"
SEED = 20261115


class ReplicationError(ValueError):
    pass


def _normalize(text: object) -> str:
    return re.sub(r"[^0-9a-z]+", "", str(text).lower())


@dataclass
class SourceTable:
    index: int
    source: ReplicationSource
    run_id: str
    inputs: list[Path]
    rows: pd.DataFrame        # kind, layer, tissue, key, symbol, feature_id, effect, se, pvalue


def load_source(project: Project, index: int, source: ReplicationSource) -> SourceTable:
    root = project.resolve(source.project) if source.project else project.root
    try:
        bundle = open_run(source.experiment, source.run, project_dir=root,
                          requires="omics/*/tables/differential.parquet")
    except ContractError as exc:
        raise ReplicationError(f"来源 {index}（{source.experiment}）：{exc}") from exc
    inputs = sorted(bundle.path.glob("omics/*/tables/*differential.parquet"))
    parts = []
    for path in inputs:
        table = pd.read_parquet(path)
        table = table[table["contrast_id"].astype(str) == source.contrast]
        if "estimand_type" in table.columns:
            table = table[table["estimand_type"].fillna("pairwise") == "pairwise"]
        if table.empty:
            continue
        derived = "feature_type" in table.columns and path.name.startswith("derived")
        se, _ = effect_standard_errors(table)
        names = table["feature_name"].fillna(table["feature_id"]) if "feature_name" in table.columns \
            else table["feature_id"]
        omics = table["omics"].astype(str)
        key = np.where(derived, table["feature_id"].astype(str),
                       np.where(omics == "metabolomics", names.map(_normalize), table["feature_id"].astype(str)))
        parts.append(pd.DataFrame({
            "kind": "derived" if derived else "feature", "layer": omics, "tissue": table["tissue"].astype(str),
            "key": key, "symbol": names.astype(str).str.upper(), "feature_id": table["feature_id"].astype(str),
            "effect": pd.to_numeric(table["log2fc"], errors="coerce"), "se": se,
            "pvalue": pd.to_numeric(table["pvalue"], errors="coerce")}))
    enrichment = sorted(bundle.path.glob("omics/*/tables/enrichment.parquet"))
    for path in enrichment:
        table = pd.read_parquet(path)
        table = table[(table["contrast_id"].astype(str) == source.contrast) & (table["method"] == "gsea_prerank")]
        if table.empty:
            continue
        table = table.sort_values("pvalue").drop_duplicates(["tissue", "omics", "term_name"])
        parts.append(pd.DataFrame({
            "kind": "pathway", "layer": table["omics"].astype(str), "tissue": table["tissue"].astype(str),
            "key": table["term_name"].map(_normalize), "symbol": "", "feature_id": table["term_name"].astype(str),
            "effect": pd.to_numeric(table["nes"], errors="coerce"), "se": np.nan,
            "pvalue": pd.to_numeric(table["pvalue"], errors="coerce")}))
    inputs += enrichment
    if not parts:
        raise ReplicationError(f"来源 {index}（{source.experiment} 的运行 {bundle.run_id}）中没有比较 {source.contrast} 的结果")
    rows = pd.concat(parts, ignore_index=True)
    rows["tissue"] = rows["tissue"].map(lambda t: source.tissue_map.get(t, t))
    rows = rows.dropna(subset=["effect"]).drop_duplicates(["kind", "layer", "tissue", "key"])
    return SourceTable(index, source, bundle.run_id, inputs, rows)


def partial_conjunction(pvalues: np.ndarray, effects: np.ndarray, u: int) -> np.ndarray:
    """方向性的部分合取 p 值。``pvalues`` / ``effects``：特征 × 来源（缺失为 NaN，按 p = 1 处理）。"""
    p = np.where(np.isfinite(pvalues), np.clip(pvalues, 0, 1), 1.0)
    e = np.where(np.isfinite(effects), effects, 0.0)
    n = p.shape[1]
    up = np.where(e > 0, p / 2, 1 - p / 2)
    down = np.where(e < 0, p / 2, 1 - p / 2)

    def simes(one_sided: np.ndarray) -> np.ndarray:
        ordered = np.sort(one_sided, axis=1)
        weights = np.array([(n - u + 1) / (i - u + 1) for i in range(u, n + 1)])
        return np.min(ordered[:, u - 1:] * weights, axis=1)

    return np.minimum(2 * np.minimum(simes(up), simes(down)), 1.0)


def _match(tables: list[SourceTable], cross_layer: list[list[str]]) -> list[tuple[str, str, pd.DataFrame]]:
    """(对应方式, 组学, 宽表)：宽表每行一个对应上的特征，列 effect_<i>、se_<i>、pvalue_<i>、feature_id_<i>。"""
    keys = ["kind", "layer", "tissue", "key"]
    out = []
    merged = None
    for table in tables:
        part = table.rows.set_index(keys)[["feature_id", "effect", "se", "pvalue"]].add_suffix(f"_{table.index}")
        merged = part if merged is None else merged.join(part, how="outer")
    for (kind, layer), rows in merged.groupby(level=["kind", "layer"], sort=True):
        out.append(("same_layer", layer, rows.reset_index()))
    if len(tables) == 2:
        first, second = tables
        for pair in cross_layer:
            for a, b in (tuple(pair), tuple(reversed(pair))):
                left = first.rows[(first.rows["layer"] == a) & (first.rows["kind"] == "feature")]
                right = second.rows[(second.rows["layer"] == b) & (second.rows["kind"] == "feature")]
                joined = left.merge(right, on=["tissue", "symbol"], suffixes=(f"_{first.index}", f"_{second.index}"))
                if joined.empty:
                    continue
                joined = joined.drop_duplicates(["tissue", "symbol"]).assign(kind="feature", key=joined["symbol"])
                out.append(("cross_layer", f"{a}|{b}", joined))
    return out


def replicate_set(project: Project, spec: ReplicationSet) -> tuple[pd.DataFrame, pd.DataFrame, list[SourceTable]]:
    tables = [load_source(project, i + 1, s) for i, s in enumerate(spec.sources)]
    n = len(tables)
    u = spec.min_sources or n
    results, baselines = [], []
    rng = np.random.default_rng(SEED)
    for match, layer, wide in _match(tables, [[o.value for o in pair] for pair in spec.cross_layer]):
        effects = np.column_stack([pd.to_numeric(wide.get(f"effect_{t.index}"), errors="coerce") for t in tables])
        pvalues = np.column_stack([pd.to_numeric(wide.get(f"pvalue_{t.index}"), errors="coerce") for t in tables])
        present = np.isfinite(effects).sum(axis=1)
        keep = present >= u
        wide, effects, pvalues = wide[keep].reset_index(drop=True), effects[keep], pvalues[keep]
        if wide.empty:
            continue
        signs = np.sign(np.where(np.isfinite(effects), effects, 0.0))
        positive, negative = (signs > 0).sum(axis=1), (signs < 0).sum(axis=1)
        direction = np.where(positive >= u, "up", np.where(negative >= u, "down", "mixed"))
        out = wide[["tissue", "kind", "key"]].copy()
        out.insert(0, "set_id", spec.id)
        out.insert(1, "match", match)
        out.insert(2, "layers", layer)
        for t in tables:
            for column in ("feature_id", "effect", "pvalue"):
                out[f"{column}_{t.index}"] = wide.get(f"{column}_{t.index}")
        out["n_sources"] = present[keep]
        out["direction"] = direction
        out["sign_consistent"] = (positive == present[keep]) | (negative == present[keep])
        out["pc_pvalue"] = partial_conjunction(pvalues, effects, u)
        out["replication_fdr"] = np.nan
        for _, idx in out.groupby(["tissue", "kind"]).groups.items():
            out.loc[idx, "replication_fdr"] = bh(out.loc[idx, "pc_pvalue"])
        out["replicated"] = out["replication_fdr"].lt(spec.alpha) & out["direction"].ne("mixed")
        if match == "same_layer":
            ses = np.column_stack([pd.to_numeric(wide.get(f"se_{t.index}"), errors="coerce") for t in tables])
            pooled = [random_effects_meta(e[np.isfinite(e) & np.isfinite(s)], s[np.isfinite(e) & np.isfinite(s)])
                      if (np.isfinite(e) & np.isfinite(s) & (s > 0)).sum() >= 2 else None
                      for e, s in zip(effects, ses)]
            out["pooled_effect"] = [r.pooled_effect if r else np.nan for r in pooled]
            out["pooled_se"] = [r.pooled_se if r else np.nan for r in pooled]
            out["pooled_pvalue"] = [r.pvalue if r else np.nan for r in pooled]
            out["i_squared"] = [r.i_squared if r else np.nan for r in pooled]
        results.append(out)
        for (tissue, kind), rows in out.groupby(["tissue", "kind"], sort=True):
            idx = rows.index
            both = np.isfinite(effects[idx]).all(axis=1)
            e = effects[idx][both]
            observed = float(np.mean(np.sign(e[:, 0])[:, None] == np.sign(e[:, 1:]))) if len(e) else np.nan
            shuffled = e.copy()
            if len(e):
                shuffled[:, 1:] = shuffled[rng.permutation(len(e)), 1:]
            random = float(np.mean(np.sign(shuffled[:, 0])[:, None] == np.sign(shuffled[:, 1:]))) if len(e) else np.nan
            baselines.append((spec.id, match, layer, tissue, kind, int(len(rows)), observed, random,
                              int(rows["replicated"].sum())))
    columns = ["set_id", "match", "layers", "tissue", "kind", "n_features", "sign_concordance",
               "sign_concordance_random", "n_replicated"]
    return (pd.concat(results, ignore_index=True) if results else pd.DataFrame(),
            pd.DataFrame(baselines, columns=columns), tables)


def replicate_project(root: str | Path | None = None, set_ids: list[str] | None = None,
                      command: list[str] | None = None, run_id: str | None = None) -> tuple[Path, list[str]]:
    project = load_project(root)
    sets = [s for s in project.config.replication if not set_ids or s.id in set_ids]
    if not sets:
        known = "、".join(s.id for s in project.config.replication) or "无"
        raise ReplicationError(f"project.yaml 的 replication 中没有要运行的复现组（已声明：{known}）")
    run_id = run_id or new_run_id()
    run_dir = project.path("results") / PROJECT_RUNS / run_id
    tables_dir = run_dir / "integration" / "replication" / "tables"
    tables_dir.mkdir(parents=True)
    results, baselines, inputs, sources, notes = [], [], [], {}, []
    for spec in sets:
        result, baseline, tables = replicate_set(project, spec)
        results.append(result)
        baselines.append(baseline)
        inputs += [p for t in tables for p in t.inputs]
        sources[spec.id] = [{"project": t.source.project or ".", "experiment": t.source.experiment,
                             "contrast": t.source.contrast, "run": t.run_id} for t in tables]
        if not result.empty:
            hits = result[result["replicated"]]
            notes.append(f"复现组 {spec.id}：{len(result)} 个对应上的特征 / 通路，{len(hits)} 个复现"
                         f"（{', '.join(f'{k} {v}' for k, v in hits['kind'].value_counts().items())}）")
    pd.concat(results, ignore_index=True).to_parquet(tables_dir / "replication.parquet", index=False)
    pd.concat(baselines, ignore_index=True).to_parquet(tables_dir / "replication_baseline.parquet", index=False)
    outputs = sorted(p for p in run_dir.rglob("*") if p.is_file())
    write_manifest(RunManifest(
        run_id=run_id, project_id=project.id, experiment_id=PROJECT_RUNS, command=command or [],
        parameters={"replication": sources, "seed": SEED},
        files=[file_record(p, project.root, "input") for p in dict.fromkeys(inputs)]
        + [file_record(p, run_dir, "output") for p in outputs]), run_dir)
    return run_dir, notes
