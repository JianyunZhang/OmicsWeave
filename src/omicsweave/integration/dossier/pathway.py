#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
通路层的综合读数（S5，G16；设计文档 3.8、3.9）：节点分流比、预先定号的模块检验、交叉点分析。三者都只在面板声明的
范围内进行（分支节点、带符号的模块、面板的反应顺序），不由数据挑选。

节点分流比（``tables/node_split.parquet``，逐只动物值在 ``tables/node_split_samples.parquet``）
    对面板 ``branch_nodes`` 中每个节点的每一对出口（a, b），在同一只动物内取两个出口之比，再比较两组：

    - 酶（蛋白、RNA）：两个出口反应的直接酶（转运反应为载体；共用亚基不计入），各自在线性尺度上求和后取 log2，
      相减（与 G8 的 ``ratio`` 定义相同）。两个出口共有的特征从两边去掉；
    - 代谢物：两个出口反应在节点另一侧的主要代谢物（去掉辅因子；一个代谢物有多个特征时取中位强度最高的一个），
      同样求和后相减。同一个代谢物出现在两边时去掉；
    - 通量：每种示踪剂，两个出口中以节点为底物的反应的逐只动物净通量求和，取 log2 之比。两组中有任何一只动物的
      和不为正，或任一出口反应在该示踪剂、该比较上不可辨识（``not_identifiable``、``insufficient_n``）时不计算；
    - 标记层（两个出口产物的标记之比）待 G10。

    组间差异用 :func:`~omicsweave.core.statistics.linear_model.fit_contrasts`（两组线性模型；分流比特征很少，不做
    方差压缩），``padj`` 在（面板, 层, 来源, 组织, 比较）内做 BH。``log_ratio_change`` 为 log2 尺度。

预先定号的模块检验（``tables/pathway_tests.parquet``，``test = signed_module``）
    面板 ``signed_modules`` 中每个模块的 plus / minus 基因，在蛋白层与 RNA 层分别做带符号的旋转检验
    （:func:`~omicsweave.core.statistics.rotation.signed_rotation_test`，limma 先验取该数据集全部完整特征）。一个基因
    匹配到多个特征时取中位强度最高的一个；有缺失值的成员不参与。``direction`` 为 ``as_declared``（按声明的符号重排）或
    ``opposite``；全部符号取反时 p 不变、方向相反。

交叉点分析（``test = crossover``）
    面板中每个只有一个主要底物、一个主要产物的反应，比较上游与下游代谢物的变化：上游升、下游降为 ``restricted``
    （该步骤受限），上游降、下游升为 ``activated``。p 值为两种模式中较小的合取 p（两个单侧 p 的最大值）乘以 2；
    一个代谢物有多个特征时单侧 p 取 Simes 合并。流动方向取通量：该反应可辨识且两组均值同号时按其符号（负号表示与
    书写方向相反），否则按书写方向假定（``assumed_direction``）。两侧是同一个代谢物（区室在数据中不可区分）、
    多底物或多产物的步骤不检验，原因写入 ``detail``。分支节点与辅因子偶联的步骤标注在 ``flags`` 中：交叉点定理只
    适用于无分支、无守恒基团的线性片段，这里的结果只作为假说来源。
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Iterable

import numpy as np
import pandas as pd

from omicsweave.core.statistics.linear_model import fit_contrasts
from omicsweave.core.statistics.rotation import signed_rotation_test

from ..mappings.network_annotation import network_reactions, split_metabolite
from .assessment import bh, one_sided, simes
from .panels import COFACTORS, Panel, _main_species
from .view import node_members, species_key

NODE_SPLIT_COLUMNS = ["panel_id", "node", "exit_a", "exit_b", "layer", "source", "tissue", "contrast_id", "tracer",
                      "state", "log_ratio_change", "se", "pvalue", "padj", "n_numerator", "n_denominator",
                      "members_a", "members_b", "detail"]
SPLIT_SAMPLE_COLUMNS = ["panel_id", "node", "exit_a", "exit_b", "layer", "source", "tissue", "tracer", "sample_id",
                        "subject_id", "group", "group_order", "value"]
TEST_COLUMNS = ["panel_id", "test", "name", "layer", "source", "tissue", "contrast_id", "statistic", "pvalue", "padj",
                "direction", "n_members", "members", "flags", "detail"]
ENZYME_LAYERS = ("proteomics", "transcriptomics")
SPLIT_LAYERS = ("mfa", "proteomics", "transcriptomics", "metabolomics")
BLOCKED_FLUX = ("not_identifiable", "insufficient_n", "not_measured", "not_infused")
MIN_GROUP = 2


def _measured(rows: pd.DataFrame) -> pd.DataFrame:
    return rows[rows["feature_id"].notna() & rows["feature_id"].astype(str).ne("")]


def _log2_sum(values: pd.DataFrame) -> pd.Series:
    """各特征在线性尺度上求和后取 log2（与 G8 ``ratio`` 的定义相同；缺失的特征不计入）。"""
    return np.log2(np.power(2.0, values.astype(float)).sum(axis=0, min_count=1))


def _scope(panel: Panel, values: Iterable[str], declared: list[str]) -> list[str]:
    found = sorted(set(values))
    return [v for v in found if v in declared] if declared else found


def exit_pairs(panel: Panel) -> list[tuple[str, tuple[str, list[str]], tuple[str, list[str]]]]:
    """面板声明的每个分支节点的每一对出口（按声明顺序）。"""
    pairs = []
    for node in panel.branch_nodes:
        exits = [(str(name), [str(r) for r in reactions]) for name, reactions in (node.get("exits") or {}).items()]
        pairs += [(str(node["node"]), exits[i], exits[j]) for i in range(len(exits)) for j in range(i + 1, len(exits))]
    return pairs


def _other_side(node: str, reaction: str) -> list[str]:
    """出口反应在节点另一侧的主要代谢物（节点是底物时取产物，节点是产物时取底物）。"""
    entry = network_reactions()[reaction]
    if node in entry.substrates:
        return _main_species(entry.products)
    if node in entry.products:
        return _main_species(entry.substrates)
    return []


def _consumers(node: str, reactions: list[str]) -> list[str]:
    return [r for r in reactions if node in network_reactions()[r].substrates]


def _wide(values: pd.DataFrame) -> pd.DataFrame:
    """逐只动物长表 → 特征 × 样本。"""
    if values.empty:
        return pd.DataFrame()
    return values.drop_duplicates(["feature_id", "sample_id"]).set_index(["feature_id", "sample_id"])["value"] \
        .astype(float).unstack("sample_id")


@dataclass
class _Side:
    features: dict[str, str]           # 特征 → 成员名
    reason: str = ""


def _enzyme_side(rows: pd.DataFrame, reactions: list[str], layer: str) -> _Side:
    """出口反应的执行者，规则与图 A 的酶格子相同：直接酶；反应没有直接酶时取载体；共用亚基不计入。"""
    features: dict[str, str] = {}
    for reaction in reactions:
        part = rows[(rows["flux_id"] == reaction) & (rows["layer"] == layer)]
        chosen = part[part["role"] == "direct_enzyme"]
        if chosen[chosen["cell_state"] != "no_knowledge"].empty:
            chosen = part[part["role"].isin(("carrier", "direct_enzyme"))]
        chosen = _measured(chosen[chosen["member_relation"].fillna("") != "shared_subunit"])
        features.update(zip(chosen["feature_id"].astype(str), chosen["member_id"].astype(str)))
    return _Side(features)


def _metabolite_side(rows: pd.DataFrame, node: str, reactions: list[str], wide: pd.DataFrame,
                     keys: dict[str, str]) -> _Side:
    features: dict[str, str] = {}
    node_key = species_key(node, keys)
    for reaction in reactions:
        reaction_rows = rows[rows["flux_id"] == reaction]
        for species in _other_side(node, reaction):
            if species_key(species, keys) == node_key:
                continue
            part = _measured(node_members(reaction_rows, species, keys))
            ids = [f for f in part["feature_id"].astype(str) if f in wide.index]
            if ids:
                best = str(wide.loc[ids].median(axis=1).idxmax())
                features[best] = str(part.loc[part["feature_id"].astype(str) == best, "member_id"].iloc[0])
    return _Side(features)


def _test(split: pd.Series, numerator: list[str], denominator: list[str]) -> dict[str, Any]:
    num = [s for s in numerator if s in split.index and np.isfinite(split[s])]
    den = [s for s in denominator if s in split.index and np.isfinite(split[s])]
    out = {"n_numerator": len(num), "n_denominator": len(den)}
    if len(num) < MIN_GROUP or len(den) < MIN_GROUP:
        return {**out, "state": "insufficient_n", "detail": f"有值的动物 {len(num)} vs {len(den)}（每组至少 {MIN_GROUP} 只）"}
    samples = num + den
    x = np.column_stack([np.ones(len(samples)), [1.0] * len(num) + [0.0] * len(den)])
    table = fit_contrasts(split[samples].to_frame().T, x, {"effect": np.array([0.0, 1.0])}, moderate=False)["effect"].table
    row = table.iloc[0]
    return {**out, "state": "computed", "log_ratio_change": float(row["estimate"]), "se": float(row["se"]),
            "pvalue": float(row["pvalue"]), "detail": ""}


def node_splits(panel: Panel, members: pd.DataFrame, samples: pd.DataFrame, sides: pd.DataFrame,
                sources: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """返回（node_split 表, 逐只动物的分流比）。"""
    rows_out: list[dict[str, Any]] = []
    sample_out: list[pd.DataFrame] = []
    keys: dict[str, str] = {}
    pairs = exit_pairs(panel)
    if not pairs:
        return pd.DataFrame(columns=NODE_SPLIT_COLUMNS), pd.DataFrame(columns=SPLIT_SAMPLE_COLUMNS)
    for source in sources:
        src_rows = members[members["source"] == source]
        for tissue in _scope(panel, src_rows["tissue"], panel.tissues):
            rows = src_rows[src_rows["tissue"] == tissue]
            contrasts = _scope(panel, rows["contrast_id"], panel.contrasts)
            values = samples[(samples["source"] == source) & (samples["tissue"] == tissue)]
            for node, (name_a, reactions_a), (name_b, reactions_b) in pairs:
                base = {"panel_id": panel.id, "node": node, "exit_a": name_a, "exit_b": name_b, "source": source,
                        "tissue": tissue}
                for layer in SPLIT_LAYERS:
                    if not (rows["layer"] == layer).any():
                        continue
                    layer_values = values[values["layer"] == layer]
                    tracers = sorted(set(layer_values["tracer"].dropna()) - {""}) if layer == "mfa" else [""]
                    for tracer in tracers or [""]:
                        here = layer_values[layer_values["tracer"] == tracer] if layer == "mfa" else layer_values
                        wide = _wide(here)
                        split, side_a, side_b, problem = _split(layer, node, reactions_a, reactions_b, rows, wide, keys)
                        meta = {**base, "layer": layer, "tracer": tracer,
                                "members_a": ", ".join(sorted(set(side_a.features.values()))),
                                "members_b": ", ".join(sorted(set(side_b.features.values())))}
                        if split is not None:
                            info = here.drop_duplicates("sample_id").set_index("sample_id")
                            sample_out.append(pd.DataFrame({
                                **{k: meta[k] for k in ("panel_id", "node", "exit_a", "exit_b", "layer", "source",
                                                        "tissue", "tracer")},
                                "sample_id": split.index.astype(str),
                                "subject_id": info["subject_id"].reindex(split.index).to_numpy(),
                                "group": info["group"].reindex(split.index).to_numpy(),
                                "group_order": info["group_order"].reindex(split.index).to_numpy(),
                                "value": split.to_numpy(dtype=float)}))
                        for contrast in contrasts:
                            record = {**meta, "contrast_id": contrast}
                            if problem:
                                rows_out.append({**record, "state": problem[0],
                                                 "detail": "；".join(filter(None, [problem[1], side_a.reason]))})
                                continue
                            blocked = _flux_blocked(rows, contrast, tracer, reactions_a, reactions_b, node) \
                                if layer == "mfa" else None
                            if blocked:
                                rows_out.append({**record, "state": "not_identifiable", "detail": blocked})
                                continue
                            side = sides[(sides["source"] == source) & (sides["layer"] == layer)
                                         & (sides["contrast_id"] == contrast)]
                            numerator = list(side.loc[side["side"] == "numerator", "sample_id"].astype(str))
                            denominator = list(side.loc[side["side"] == "denominator", "sample_id"].astype(str))
                            result = _test(split, numerator, denominator)
                            result["detail"] = "；".join(filter(None, [result.get("detail"), side_a.reason]))
                            rows_out.append({**record, **result})
    table = pd.DataFrame(rows_out)
    for column in NODE_SPLIT_COLUMNS:
        if column not in table.columns:
            table[column] = np.nan
    table = table[NODE_SPLIT_COLUMNS]
    if not table.empty:
        table["padj"] = table.groupby(["panel_id", "layer", "source", "tissue", "contrast_id"], group_keys=False)[
            "pvalue"].apply(bh)
    split_samples = pd.concat(sample_out, ignore_index=True) if sample_out else pd.DataFrame(columns=SPLIT_SAMPLE_COLUMNS)
    return table, split_samples[SPLIT_SAMPLE_COLUMNS]


def _split(layer: str, node: str, reactions_a: list[str], reactions_b: list[str], rows: pd.DataFrame,
           wide: pd.DataFrame, keys: dict[str, str]) -> tuple[pd.Series | None, _Side, _Side, tuple[str, str] | None]:
    """一个层（与示踪剂）的逐只动物分流比；不能计算时返回（None, …, (状态, 原因)）。"""
    if wide.empty:
        return None, _Side({}), _Side({}), ("not_measured", "该组织没有这一层的逐只动物数据")
    if layer == "mfa":
        side_a = _Side({r: r for r in _consumers(node, reactions_a) if r in wide.index})
        side_b = _Side({r: r for r in _consumers(node, reactions_b) if r in wide.index})
        missing = [name for name, side in (("a", side_a), ("b", side_b)) if not side.features]
        if missing:
            return None, side_a, side_b, ("not_measured", "出口中没有以该节点为底物、且有逐只动物通量的反应")
        sum_a, sum_b = wide.loc[list(side_a.features)].sum(axis=0, min_count=1), \
            wide.loc[list(side_b.features)].sum(axis=0, min_count=1)
        valid = sum_a.notna() & sum_b.notna()
        if (sum_a[valid] <= 0).any() or (sum_b[valid] <= 0).any():
            return None, side_a, side_b, ("not_identifiable", "有动物的出口净通量不为正，比值没有定义")
        return np.log2(sum_a / sum_b), side_a, side_b, None
    if layer == "metabolomics":
        side_a = _metabolite_side(rows, node, reactions_a, wide, keys)
        side_b = _metabolite_side(rows, node, reactions_b, wide, keys)
    else:
        side_a, side_b = _enzyme_side(rows, reactions_a, layer), _enzyme_side(rows, reactions_b, layer)
    shared = set(side_a.features) & set(side_b.features)
    note = f"两个出口共有、已去掉：{', '.join(sorted({side_a.features[f] for f in shared}))}" if shared else ""
    side_a = _Side({f: m for f, m in side_a.features.items() if f not in shared}, note)
    side_b = _Side({f: m for f, m in side_b.features.items() if f not in shared}, note)
    if shared and (not side_a.features or not side_b.features):
        return None, side_a, side_b, ("same_members", "两个出口只有共同的特征，比值恒为 0")
    a = [f for f in side_a.features if f in wide.index]
    b = [f for f in side_b.features if f in wide.index]
    if not a or not b:
        empty = "a" if not a else "b"
        return None, side_a, side_b, ("not_detected", f"出口 {empty} 在这一层没有测得的成员")
    return _log2_sum(wide.loc[a]) - _log2_sum(wide.loc[b]), side_a, side_b, None


def _flux_blocked(rows: pd.DataFrame, contrast: str, tracer: str, reactions_a: list[str], reactions_b: list[str],
                  node: str) -> str:
    part = rows[(rows["layer"] == "mfa") & (rows["contrast_id"] == contrast)]
    reasons = []
    for reaction in _consumers(node, reactions_a) + _consumers(node, reactions_b):
        hit = part[(part["flux_id"] == reaction) & (part["tracer"] == tracer)]
        if hit.empty:          # 来源没有分示踪剂的结果（如队列 B）：取合并估计的状态
            hit = part[(part["flux_id"] == reaction) & (part["tracer"] == "pooled")]
        state = str(hit["cell_state"].iloc[0]) if not hit.empty else "not_measured"
        if state in BLOCKED_FLUX:
            reasons.append(f"{reaction} {state}")
    return "；".join(reasons)


# ---------------------------------------------------------------------------
# 预先定号的模块检验
# ---------------------------------------------------------------------------

@dataclass
class ModuleData:
    """一个（来源, 层, 组织）的分析尺度数据（特征 × 样本）、基因 → 特征的匹配函数，以及估计 limma 先验时用的
    特征（``background``；None 为全部。RNA 去掉低表达基因，否则低计数基因的方差会抬高先验）。"""
    values: pd.DataFrame
    match: Callable[[list[str]], dict[str, list[str]]]
    background: pd.Index | None = None


def module_signs(panel: Panel) -> dict[str, dict[str, int]]:
    return {name: {**{g: 1 for g in spec.get("plus", [])}, **{g: -1 for g in spec.get("minus", [])}}
            for name, spec in panel.signed_modules.items()}


def signed_module_tests(panel: Panel, data: Callable[[str, str, str], ModuleData | None], sides: pd.DataFrame,
                        sources: list[str], tissues: dict[str, list[str]], contrasts: dict[str, list[str]],
                        rotations: int, seed: int) -> pd.DataFrame:
    """``data(source, layer, tissue)``：模块检验用的数据；``tissues`` / ``contrasts``：每个来源的组织与比较。"""
    rows: list[dict[str, Any]] = []
    modules = module_signs(panel)
    if not modules:
        return pd.DataFrame(columns=TEST_COLUMNS)
    for source in sources:
        for layer in ENZYME_LAYERS:
            for tissue in _scope(panel, tissues.get(source, []), panel.tissues):
                found = data(source, layer, tissue)
                if found is None:
                    continue
                for contrast in _scope(panel, contrasts.get(source, []), panel.contrasts):
                    side = sides[(sides["source"] == source) & (sides["layer"] == layer)
                                 & (sides["contrast_id"] == contrast)]
                    numerator = [s for s in side.loc[side["side"] == "numerator", "sample_id"].astype(str)
                                 if s in found.values.columns]
                    denominator = [s for s in side.loc[side["side"] == "denominator", "sample_id"].astype(str)
                                   if s in found.values.columns]
                    base = {"panel_id": panel.id, "test": "signed_module", "layer": layer, "source": source,
                            "tissue": tissue, "contrast_id": contrast}
                    if len(numerator) < MIN_GROUP or len(denominator) < MIN_GROUP:
                        rows += [{**base, "name": name, "direction": "not_tested",
                                  "detail": f"该比较的样本 {len(numerator)} vs {len(denominator)}"} for name in modules]
                        continue
                    samples = numerator + denominator
                    x = np.column_stack([np.ones(len(samples)), [1.0] * len(numerator) + [0.0] * len(denominator)])
                    values = found.values[samples]
                    pool = values if found.background is None else values.loc[values.index.intersection(found.background)]
                    complete = pool[np.isfinite(pool).all(axis=1) & (pool.var(axis=1) > 0)]
                    prior = fit_contrasts(complete, x, {"effect": np.array([0.0, 1.0])})["effect"]
                    for name, signs in modules.items():
                        rows.append({**base, "name": name, **_module_row(found, values, signs, x, prior, rotations,
                                                                          seed)})
    table = pd.DataFrame(rows)
    for column in TEST_COLUMNS:
        if column not in table.columns:
            table[column] = np.nan
    table = table[TEST_COLUMNS]
    if not table.empty:
        table["padj"] = table.groupby(["panel_id", "test", "layer", "source", "tissue", "contrast_id"],
                                      group_keys=False)["pvalue"].apply(bh)
    return table


def _module_row(found: ModuleData, values: pd.DataFrame, signs: dict[str, int], x: np.ndarray, prior: Any,
                rotations: int, seed: int) -> dict[str, Any]:
    matched = found.match(list(signs))
    chosen, missing, incomplete = {}, [], []
    for gene, sign in signs.items():
        ids = [f for f in matched.get(gene, []) if f in values.index]
        if not ids:
            missing.append(gene)
            continue
        best = str(values.loc[ids].median(axis=1).idxmax())
        if not np.isfinite(values.loc[best]).all():
            incomplete.append(gene)
            continue
        chosen[gene] = (best, sign)
    members = ", ".join(f"{g}({'+' if s > 0 else '−'})" for g, (_, s) in chosen.items())
    notes = ([f"未匹配：{'、'.join(missing)}"] if missing else []) + \
        ([f"有缺失值未参与：{'、'.join(incomplete)}"] if incomplete else [])
    if len(chosen) < 2:
        return {"direction": "not_tested", "n_members": len(chosen), "members": members,
                "detail": "；".join(["参与检验的成员少于 2 个"] + notes)}
    y = values.loc[[f for f, _ in chosen.values()]].to_numpy(dtype=float)
    result = signed_rotation_test(y, x, np.array([0.0, 1.0]), np.array([s for _, s in chosen.values()]),
                                  prior_variance=prior.prior_variance, prior_df=prior.prior_df, rotations=rotations,
                                  seed=seed)
    return {"statistic": result.statistic, "pvalue": result.pvalue,
            "direction": "as_declared" if result.direction > 0 else "opposite" if result.direction < 0 else "none",
            "n_members": result.n_members, "members": members,
            "detail": "；".join([f"按声明方向 z > √2 的成员占 {result.active_up:.0%}、反向 {result.active_down:.0%}；"
                                f"旋转 {result.rotations} 次，自由度 {result.df:.1f}"] + notes)}


# ---------------------------------------------------------------------------
# 交叉点分析
# ---------------------------------------------------------------------------

def _node_pvalues(part: pd.DataFrame) -> tuple[float, float, float]:
    """（中位 log2FC, 上升的单侧 p, 下降的单侧 p）；多个特征时单侧 p 取 Simes 合并。"""
    fc = pd.to_numeric(part["log2fc"], errors="coerce")
    up = simes([one_sided(p, f, 1) for p, f in zip(part["pvalue"], fc)])
    down = simes([one_sided(p, f, -1) for p, f in zip(part["pvalue"], fc)])
    return float(fc.median()), up, down


def _flow(rows: pd.DataFrame, reaction: str) -> tuple[int, str]:
    """流动方向：+1 按书写方向，−1 相反；第二项为依据（data / assumed）。"""
    hit = rows[(rows["flux_id"] == reaction) & (rows["layer"] == "mfa") & (rows["tracer"] == "pooled")]
    if not hit.empty and str(hit["cell_state"].iloc[0]) not in BLOCKED_FLUX:
        means = pd.to_numeric(hit[["mean_numerator", "mean_denominator"]].iloc[0], errors="coerce")
        if means.notna().all() and (np.sign(means) == np.sign(means.iloc[0])).all() and means.iloc[0] != 0:
            return int(np.sign(means.iloc[0])), "data"
    return 1, "assumed"


def crossover_tests(panel: Panel, members: pd.DataFrame, sources: list[str]) -> pd.DataFrame:
    keys: dict[str, str] = {}
    branch = {str(n.get("node")) for n in panel.branch_nodes}
    usage: dict[str, int] = {}
    network = network_reactions()
    for reaction in panel.reactions:
        for species in {*_main_species(network[reaction].substrates), *_main_species(network[reaction].products)}:
            usage[species] = usage.get(species, 0) + 1
    rows_out: list[dict[str, Any]] = []
    for source in sources:
        src = members[(members["source"] == source) & members["flux_id"].isin(panel.reactions)]
        if not (src["layer"] == "metabolomics").any():
            continue
        for tissue in _scope(panel, src["tissue"], panel.tissues):
            for contrast in _scope(panel, src.loc[src["tissue"] == tissue, "contrast_id"], panel.contrasts):
                rows = src[(src["tissue"] == tissue) & (src["contrast_id"] == contrast)]
                for reaction in panel.reactions:
                    rows_out.append({"panel_id": panel.id, "test": "crossover", "name": reaction,
                                     "layer": "metabolomics", "source": source, "tissue": tissue,
                                     "contrast_id": contrast,
                                     **_crossover(reaction, rows, keys, branch, usage)})
    table = pd.DataFrame(rows_out)
    for column in TEST_COLUMNS:
        if column not in table.columns:
            table[column] = np.nan
    table = table[TEST_COLUMNS]
    if not table.empty:
        table["padj"] = table.groupby(["panel_id", "test", "source", "tissue", "contrast_id"],
                                      group_keys=False)["pvalue"].apply(bh)
    return table


def _crossover(reaction: str, rows: pd.DataFrame, keys: dict[str, str], branch: set[str],
               usage: dict[str, int]) -> dict[str, Any]:
    entry = network_reactions()[reaction]
    substrates, products = _main_species(entry.substrates), _main_species(entry.products)
    flags = []
    if any(split_metabolite(s)[0].upper() in COFACTORS for s in (*entry.substrates, *entry.products)):
        flags.append("cofactor_coupled")
    if len(substrates) != 1 or len(products) != 1:
        return {"direction": "not_tested", "flags": ";".join(flags),
                "detail": f"不是单底物、单产物的步骤（主要底物 {len(substrates)} 个、产物 {len(products)} 个）"}
    substrate, product = substrates[0], products[0]
    if species_key(substrate, keys) == species_key(product, keys):
        return {"direction": "not_tested", "flags": ";".join(flags),
                "detail": f"{substrate} 与 {product} 是同一个代谢物（数据不区分区室）"}
    if {substrate, product} & branch or max(usage.get(substrate, 0), usage.get(product, 0)) > 2:
        flags.append("branch")
    sign, basis = _flow(rows, reaction)
    upstream, downstream = (substrate, product) if sign > 0 else (product, substrate)
    if basis == "assumed":
        flags.append("assumed_direction")
    elif sign < 0:
        flags.append("reversed")
    reaction_rows = rows[rows["flux_id"] == reaction]
    parts = {name: _measured(node_members(reaction_rows, species, keys)) for name, species in
             (("upstream", upstream), ("downstream", downstream))}
    absent = [name for name, part in parts.items() if part.empty or pd.to_numeric(part["pvalue"], errors="coerce")
              .notna().sum() == 0]
    members = f"{upstream} → {downstream}"
    if absent:
        return {"direction": "not_tested", "members": members, "flags": ";".join(flags),
                "detail": f"{'、'.join(absent)} 的代谢物没有测得或没有检验结果"}
    (fc_up, up_u, down_u), (fc_down, up_w, down_w) = _node_pvalues(parts["upstream"]), \
        _node_pvalues(parts["downstream"])
    restricted, activated = max(up_u, down_w), max(down_u, up_w)
    pattern = "restricted" if fc_up > 0 > fc_down else "activated" if fc_up < 0 < fc_down else "concordant"
    return {"statistic": fc_up - fc_down, "pvalue": min(1.0, 2 * min(restricted, activated)), "direction": pattern,
            "n_members": len(parts["upstream"]) + len(parts["downstream"]), "members": members,
            "flags": ";".join(flags),
            "detail": f"上游 {fc_up:+.2f}、下游 {fc_down:+.2f}（log2FC）；流动方向{'按通量' if basis == 'data' else '按书写方向假定'}"}
