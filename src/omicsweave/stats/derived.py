#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
声明式派生特征（G8，experiment.yaml 的 ``derived``）。

在每个组织的分析层数据集（log2 尺度；计数为 log2 CPM）上按样本计算，写入 ``omics/<组学>/derived/<组织>/``
（与运行器已写出的派生特征合并，同 ID 的替换），再对每个两两比较做两组 moderated t 检验，各派生类型在
“组织 × 比较”内单独做 BH 校正，追加到 ``tables/derived_differential.parquet``。

成员不足 ``min_members`` 或未测到时该特征不计算，原因写入运行提示。``features.parquet`` 记录 ``feature_type``、
``definition``、``generated_from``，成员带有鉴定可信等级列（``identification_level`` / ``msi_level`` / ``Level``）时
取成员中的最差等级。

``step_ratio`` 的 ``source: network`` 由 MFA 网络注释表（``network_reactions.tsv``）生成：只取去掉辅因子后恰有一个底物和
一个产物的反应，代谢物按 :func:`~omicsweave.integration.mappings.match_metabolite_roles` 匹配（与 flux_anchored 相同），
一个代谢物命中多个特征时取中位强度最高的一个。
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from omicsweave.contract.dataset import Dataset, DatasetMeta, ValueType, read_dataset, write_dataset
from omicsweave.contract.tables import read_differential, write_differential
from omicsweave.core.config import Experiment
from omicsweave.core.idmap.lipids import parse_lipid
from omicsweave.core.model import DerivedSpec, FeatureSelector, Omics
from omicsweave.core.statistics.linear_model import bh, fit_contrasts

from .data import analysis_values, bundle_datasets, pairwise_groups

FEATURE_COLUMNS = ["feature_type", "feature_name", "definition", "n_members", "members", "related_reaction",
                   "readout", "units", "generated_from", "identification_level"]
_NAME_COLUMNS = ("feature_name", "Compounds", "compound", "name", "gene_symbol")
_LEVEL_COLUMNS = ("identification_level", "msi_level", "MSI level", "Level")
# 不作为一步反应底物 / 产物的辅因子与小分子
COFACTORS = {"nad", "nadh", "nadp", "nadph", "atp", "adp", "amp", "gtp", "gdp", "coa", "h2o", "co2", "pi", "ppi",
             "h", "o2", "fad", "fadh2", "nh3", "nh4"}
UNITS = {"sum": "log2 summed abundance", "ratio": "log2 ratio", "step_ratio": "log2 product / substrate ratio",
         "fraction": "log2 fraction", "index": "log2 abundance-weighted index", "module_score": "mean z-score"}


def _names(dataset: Dataset) -> pd.Series:
    for column in _NAME_COLUMNS:
        if column in dataset.features.columns:
            return dataset.features[column].fillna(pd.Series(dataset.features.index, index=dataset.features.index)
                                                   ).astype(str)
    return pd.Series(dataset.matrix.index.astype(str), index=dataset.matrix.index)


def _database_ids(dataset: Dataset, names: pd.Series) -> pd.Series:
    from omicsweave.metabolomics.derived_features import feature_database_ids
    table = dataset.features.assign(_name=names)
    return feature_database_ids(table, "_name")


def select(selector: FeatureSelector, dataset: Dataset, names: pd.Series) -> list[str]:
    """按 ID / 名称、脂质类别或数据库 ID 选取特征（并集，保持数据集中的顺序）。"""
    wanted_ids = set(selector.features)
    classes = set(selector.lipid_class)
    chosen = set()
    if wanted_ids:
        chosen |= {f for f in dataset.matrix.index if f in wanted_ids or names.get(f) in wanted_ids}
    if classes:
        chosen |= {f for f in dataset.matrix.index if (p := parse_lipid(names.get(f))) and p.lipid_class in classes}
    if selector.database_ids:
        from omicsweave.metabolomics.derived_features import _normalize_id
        targets = {n for i in selector.database_ids if (n := _normalize_id(i))}
        ids = _database_ids(dataset, names)
        chosen |= {f for f, found in ids.items() if found & targets}
    return [f for f in dataset.matrix.index if f in chosen]


def _log2_sum(values: pd.DataFrame) -> pd.Series:
    return np.log2(np.power(2.0, values.astype(float)).sum(axis=0, min_count=1))


def _best(ids: list[str], values: pd.DataFrame) -> str:
    return str(values.loc[ids].median(axis=1).idxmax())


def _level(dataset: Dataset, members: list[str]) -> float:
    for column in _LEVEL_COLUMNS:
        if column in dataset.features.columns:
            levels = pd.to_numeric(dataset.features.loc[members, column], errors="coerce").dropna()
            if not levels.empty:
                return float(levels.max())
    return float("nan")


def _definition(spec_type: str, name: str, text: str, members: list[str], dataset: Dataset, *,
                reaction: str = "", readout: str = "", origin: str = "") -> dict:
    return {"feature_type": spec_type, "feature_name": name, "definition": text, "n_members": len(members),
            "members": ";".join(members), "related_reaction": reaction, "readout": readout, "units": UNITS[spec_type],
            "generated_from": origin, "identification_level": _level(dataset, members)}


def network_step_ratios(dataset: Dataset, values: pd.DataFrame, names: pd.Series
                        ) -> tuple[list[tuple[str, pd.Series, dict]], list[str]]:
    from omicsweave.integration.mappings import get_metabolite_roles, match_metabolite_roles
    from omicsweave.integration.mappings.network_annotation import network_reactions

    frame = pd.DataFrame({"compound": names}, index=dataset.matrix.index)
    for source, target in (("HMDB", "hmdb_id"), ("cpd_ID", "kegg_compound"), ("KEGG_Compound", "kegg_compound"),
                           ("ChEBI", "chebi_id")):
        if source in dataset.features.columns and target not in frame.columns:
            frame[target] = dataset.features[source]
    out, notes = [], []
    for flux_id, reaction in sorted(network_reactions().items()):
        substrates = {n for n in reaction.tissue_substrate_names if n.lower() not in COFACTORS}
        products = {n for n in reaction.tissue_product_names if n.lower() not in COFACTORS}
        roles = get_metabolite_roles(flux_id)
        if len(substrates) != 1 or len(products) != 1 or not roles:
            continue
        sub_aliases = [a for a in roles.get("substrates", []) if a.lower() in {s.lower() for s in substrates}
                       or _same(a, next(iter(substrates)))]
        prod_aliases = [a for a in roles.get("products", []) if a.lower() in {p.lower() for p in products}
                        or _same(a, next(iter(products)))]
        if not sub_aliases or not prod_aliases:
            continue
        matched = match_metabolite_roles(frame, {"substrates": sub_aliases, "products": prod_aliases, "related": []})
        sub_hits = [str(i) for i in matched.indices_by_role["substrates"]]
        prod_hits = [str(i) for i in matched.indices_by_role["products"]]
        if not sub_hits or not prod_hits:
            continue
        substrate, product = _best(sub_hits, values), _best(prod_hits, values)
        name = f"{next(iter(products))} / {next(iter(substrates))}"
        out.append((f"step_ratio:{flux_id}", values.loc[product] - values.loc[substrate], _definition(
            "step_ratio", name, f"log2({names.get(product, product)}) - log2({names.get(substrate, substrate)})",
            [product, substrate], dataset, reaction=flux_id, readout=f"{flux_id} product / substrate",
            origin=f"network:{reaction.network}:{flux_id}")))
    return out, notes


def _same(alias: str, canonical: str) -> bool:
    normalize = lambda s: "".join(ch for ch in s.lower() if ch.isalnum())  # noqa: E731
    return normalize(alias) == normalize(canonical)


def compute(spec: DerivedSpec, dataset: Dataset, values: pd.DataFrame, names: pd.Series
            ) -> tuple[list[tuple[str, pd.Series, dict]], list[str]]:
    """返回 ([(ID, 按样本的值, 定义)], 未计算的原因)。"""
    if spec.type == "step_ratio" and spec.source == "network":
        return network_step_ratios(dataset, values, names)
    feature_id = f"{spec.type}:{spec.id}"
    origin = f"config:{spec.id}"

    def members(selector: FeatureSelector | None, label: str) -> list[str] | str:
        found = select(selector, dataset, names) if selector is not None else []
        if len(found) < spec.min_members:
            return f"{spec.id}：{label}只匹配到 {len(found)} 个特征（至少 {spec.min_members} 个），未计算"
        return found

    if spec.type in ("ratio", "step_ratio", "fraction"):
        top = members(spec.numerator if spec.type != "fraction" else spec.members, "分子")
        bottom = members(spec.denominator, "分母")
        problems = [m for m in (top, bottom) if isinstance(m, str)]
        if problems:
            return [], problems
        if spec.type == "step_ratio":
            top, bottom = [_best(top, values)], [_best(bottom, values)]
        series = _log2_sum(values.loc[top]) - _log2_sum(values.loc[bottom])
        text = f"log2(sum of {len(top)}) - log2(sum of {len(bottom)})"
        return [(feature_id, series, _definition(spec.type, spec.id, text, list(top) + list(bottom), dataset,
                                                 origin=origin))], []
    found = members(spec.members, "成员")
    if isinstance(found, str):
        return [], [found]
    if spec.type == "sum":
        return [(feature_id, _log2_sum(values.loc[found]),
                 _definition("sum", spec.id, f"log2 of the summed abundance of {len(found)} features", found,
                             dataset, origin=origin))], []
    if spec.type == "index":
        parsed = {f: parse_lipid(names.get(f)) for f in found}
        parsed = {f: p for f, p in parsed.items() if p is not None}
        if len(parsed) < spec.min_members:
            return [], [f"{spec.id}：能解析出链组成的脂质只有 {len(parsed)} 个，未计算"]
        ids = list(parsed)
        x = pd.Series({f: (p.mean_chain_length if spec.index == "chain_length" else p.double_bonds)
                       for f, p in parsed.items()})
        weights = np.power(2.0, values.loc[ids].astype(float))
        index = weights.mul(x, axis=0).sum(axis=0, min_count=1) / weights.sum(axis=0, min_count=1)
        # 双键数可以为 0，取 log2(1 + 指数)；平均链长总为正，取 log2(指数)
        series = np.log2(index) if spec.index == "chain_length" else np.log2(1.0 + index)
        what = "mean acyl chain length" if spec.index == "chain_length" else "1 + mean double bonds"
        return [(feature_id, series, _definition("index", spec.id, f"log2 abundance-weighted {what} over {len(ids)} "
                                                 "lipids", ids, dataset, origin=origin))], []
    # module_score：各成员在样本间标准化后取平均
    block = values.loc[found].astype(float)
    z = block.sub(block.mean(axis=1), axis=0).div(block.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    return [(feature_id, z.mean(axis=0, skipna=True),
             _definition("module_score", spec.id, f"mean z-score of {len(found)} features", found, dataset,
                         origin=origin))], []


def _differential(experiment: Experiment, omics: str, tissue: str, matrix: pd.DataFrame, features: pd.DataFrame,
                  alpha: float) -> pd.DataFrame:
    parts = []
    for contrast in experiment.config.design.contrasts:
        groups = pairwise_groups(experiment, contrast, omics, matrix.columns)
        num = [s for s, g in groups.items() if g == contrast.numerator]
        den = [s for s, g in groups.items() if g == contrast.denominator]
        if not num or not den:
            continue
        samples = num + den
        x = np.column_stack([np.ones(len(samples)), [1.0] * len(num) + [0.0] * len(den)])
        fit = fit_contrasts(matrix[samples], x, {"effect": np.array([0.0, 1.0])})["effect"].table
        padj = pd.Series(np.nan, index=fit.index)
        for _, ids in features.groupby("feature_type").groups.items():
            padj.loc[ids] = bh(fit.loc[ids, "pvalue"])
        tested = fit["pvalue"].notna()
        out = pd.DataFrame({
            "contrast_id": contrast.id, "contrast_role": contrast.role.value, "tissue": tissue, "omics": omics,
            "feature_id": fit.index.astype(str), "feature_name": features["feature_name"].reindex(fit.index),
            "log2fc": fit["estimate"], "lfc_se": fit["se"], "statistic": fit["statistic"], "pvalue": fit["pvalue"],
            "padj": padj, "padj_method": np.where(tested, "fdr_bh", "none"), "n_numerator": len(num),
            "n_denominator": len(den), "method": "moderated_lm", "mean_numerator": matrix[num].mean(axis=1),
            "mean_denominator": matrix[den].mean(axis=1), "units": features["units"].reindex(fit.index),
            "inference_unit": "biological_replicate",
            "evidence_level": np.where(tested, "statistical_test", "insufficient"),
        })
        out["selected"] = out["padj"].lt(alpha)
        out["selection_rule"] = f"padj < {alpha:g}（各派生类型单独校正）"
        for column in ("feature_type", "definition", "related_reaction", "readout", "n_members", "generated_from"):
            out[column] = features[column].reindex(fit.index).to_numpy()
        parts.append(out.reset_index(drop=True))
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def write_derived(experiment: Experiment, run_dir: str | Path) -> list[str]:
    run_dir = Path(run_dir)
    alpha = float(experiment.project.config.defaults.alpha)
    notes: list[str] = []
    by_omics: dict[str, list[DerivedSpec]] = {}
    for spec in experiment.config.derived:
        by_omics.setdefault(spec.omics.value, []).append(spec)
    for omics, specs in by_omics.items():
        rows_all = []
        for _, kind, tissue, path in [d for d in bundle_datasets(run_dir, ("data",)) if d[0] == omics]:
            dataset = read_dataset(path)
            values, names = analysis_values(dataset), _names(dataset)
            computed = []
            for spec in specs:
                found, problems = compute(spec, dataset, values, names)
                computed += found
                notes += [f"{omics} {tissue}：{p}" for p in problems]
            if not computed:
                continue
            matrix = pd.DataFrame({fid: series for fid, series, _ in computed}).T.reindex(columns=dataset.matrix.columns)
            matrix.index.name = "feature_id"
            features = pd.DataFrame([definition for _, _, definition in computed], index=matrix.index)[FEATURE_COLUMNS]
            target = run_dir / "omics" / omics / "derived" / tissue
            if (target / "dataset.json").is_file():
                existing = read_dataset(target)
                keep = existing.matrix.index.difference(matrix.index, sort=False)
                matrix = pd.concat([existing.matrix.loc[keep], matrix])
                features = pd.concat([existing.features.loc[keep], features]).reindex(columns=list(dict.fromkeys(
                    [*existing.features.columns, *FEATURE_COLUMNS])))
                meta = existing.meta
            else:
                meta = DatasetMeta(omics=Omics(omics), tissue=tissue, value_type=ValueType.DERIVED, scale="log2",
                                   normalization=dataset.meta.normalization, feature_id_type="derived_feature",
                                   source={"from": f"data/{tissue}", "module": "omicsweave.stats.derived"})
            write_dataset(target, Dataset(matrix, features, meta))
            new_ids = [fid for fid, _, _ in computed]
            rows_all.append(_differential(experiment, omics, tissue, matrix.loc[new_ids], features.loc[new_ids], alpha))
        rows_all = [r for r in rows_all if not r.empty]
        if not rows_all:
            continue
        new = pd.concat(rows_all, ignore_index=True)
        counts = new.drop_duplicates(["tissue", "feature_id"])["feature_type"].value_counts()
        path = run_dir / "omics" / omics / "tables" / "derived_differential.parquet"
        if path.is_file():
            old = read_differential(path)
            keys = set(zip(new["contrast_id"], new["tissue"], new["feature_id"]))
            old = old[[k not in keys for k in zip(old["contrast_id"], old["tissue"], old["feature_id"])]]
            new = pd.concat([old, new], ignore_index=True)
        write_differential(new, path)
        notes.append(f"{omics} 声明式派生特征（组织 × 特征）：" + "、".join(f"{t} {int(n)}" for t, n in counts.items()))
    return notes
