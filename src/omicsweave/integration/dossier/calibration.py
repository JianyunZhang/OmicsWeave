#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
判定规则的校准（E8，设计文档 3.11）：给分轴综合评估的结论一个可估计的错误率。

- 诱饵映射：把每个反应在每个组织中的直接酶（转运反应为载体）与底物 / 产物换成同一组织数据中丰度相近的随机特征
  （E3 的抽样），重跑成员行与四轴判定；统计诱饵达到“C 轴 changed”“C 轴 changed 且 M 轴同向”与 L1 的个数，
  与目标对照。诱饵没有调控因子知识（R 轴为 no_knowledge）、不带派生比值；通量行与目标相同。
- 对照集：面板声明中的阳性对照（已知方向）与阴性对照（预期不变），统计判定与预期的符合率。
- 不变性：翻转通量的符号，C、R、M 轴不变；删除一整层，不读该层的轴不变。在实际数据上运行（面板反应，没有面板时
  取前 30 个反应），违反时记录格子数。
- 组别置换（E1）需要按置换后的标签重跑各组学的差异分析，尚未实现，在结果表中记为未运行。

输出 ``evaluation/assessment_calibration.parquet``（列见 :data:`COLUMNS`），图 A 的图注写出该组织、该比较的诱饵对照。
"""
from __future__ import annotations

from typing import Any, Callable, Iterable

import numpy as np
import pandas as pd

COLUMNS = ["method", "tissue", "contrast_id", "metric", "n_cells", "target", "decoy_expected", "decoy_sd",
           "n_replicates", "ratio", "passed", "detail"]
# 每条轴读取的层：删除一整层时，不读该层的轴应当不变
AXIS_LAYERS = {"flux": {"mfa"}, "capacity": {"proteomics", "transcriptomics"},
               "regulation": {"proteomics", "transcriptomics", "metabolomics"}, "pool": {"metabolomics", "derived"}}
LAYERS = ("mfa", "proteomics", "transcriptomics", "metabolomics", "derived")


def _flags(axes: pd.DataFrame, reactions: pd.DataFrame) -> pd.DataFrame:
    """每个（反应, 组织, 比较）是否达到各项：L1、C 轴 changed、C 轴 changed 且 M 轴同向。参数顺序与
    :func:`~omicsweave.integration.dossier.assessment.assess` 的返回值相同。"""
    keys = ["flux_id", "tissue", "contrast_id"]
    wide = axes.drop_duplicates(keys + ["axis"]).set_index(keys + ["axis"])[["state", "direction"]].unstack("axis")
    capacity_changed = wide[("state", "capacity")].eq("changed")
    pool_same = wide[("state", "pool")].isin(["changed", "suggestive"]) & \
        (wide[("direction", "pool")] == wide[("direction", "capacity")]) & (wide[("direction", "capacity")] != 0)
    out = pd.DataFrame({"capacity_changed": capacity_changed, "capacity_and_pool": capacity_changed & pool_same})
    claims = reactions.set_index(keys)["flux_claim"]
    out["L1"] = claims.reindex(out.index).eq("L1")
    return out.reset_index()


METRICS = ("L1", "capacity_changed", "capacity_and_pool")


def decoy_summary(target: tuple[pd.DataFrame, pd.DataFrame],
                  decoys: list[tuple[pd.DataFrame, pd.DataFrame]]) -> pd.DataFrame:
    """目标与各次诱饵重复的计数：按（组织, 比较）与全部。"""
    target_flags = _flags(*target)
    decoy_flags = [_flags(*run) for run in decoys]
    rows = []
    groups = [("", "", target_flags.index)] + [(t, c, part.index) for (t, c), part in
                                                 target_flags.groupby(["tissue", "contrast_id"], sort=True)]
    for tissue, contrast, index in groups:
        keys = target_flags.loc[index, ["flux_id", "tissue", "contrast_id"]]
        for metric in METRICS:
            observed = int(target_flags.loc[index, metric].sum())
            counts = []
            for flags in decoy_flags:
                merged = keys.merge(flags, on=["flux_id", "tissue", "contrast_id"], how="left")
                counts.append(int(merged[metric].fillna(False).astype(bool).sum()))
            expected = float(np.mean(counts)) if counts else np.nan
            rows.append({"method": "decoy", "tissue": tissue, "contrast_id": contrast, "metric": metric,
                         "n_cells": len(index), "target": observed, "decoy_expected": expected,
                         "decoy_sd": float(np.std(counts, ddof=1)) if len(counts) > 1 else np.nan,
                         "n_replicates": len(counts),
                         "ratio": min(expected / observed, 1.0) if observed and counts else np.nan,
                         "passed": np.nan, "detail": "；".join(str(c) for c in counts)})
    return pd.DataFrame(rows, columns=COLUMNS)


def control_summary(reactions: pd.DataFrame, axes: pd.DataFrame, controls: Iterable[dict[str, Any]]) -> pd.DataFrame:
    """对照集：阳性对照看指定的结论或轴是否指向预期方向；阴性对照看是否没有给出方向。"""
    rows = []
    claims = reactions.set_index(["flux_id", "tissue", "contrast_id"])
    axis_table = axes.set_index(["flux_id", "tissue", "contrast_id", "axis"])
    for control in controls:
        key = (control["reaction"], control["tissue"], control["contrast"])
        target = control.get("axis", "claim")
        if key not in claims.index:
            observed, correct = "not_assessed", np.nan
        elif target == "claim":
            claim = claims.loc[key]
            direction = int(claim["claim_direction"])
            observed = f"{claim['flux_claim']} {'up' if direction > 0 else 'down' if direction < 0 else ''}".strip()
            correct = (direction == {"up": 1, "down": -1}[control["direction"]]) if control["kind"] == "positive" \
                else claim["flux_claim"] in ("L0",)
        else:
            axis = axis_table.loc[(*key, target)]
            direction = int(axis["direction"])
            pointing = axis["state"] in ("changed", "suggestive") and direction != 0
            observed = f"{axis['state']} {'up' if direction > 0 else 'down' if direction < 0 else ''}".strip()
            correct = (pointing and direction == {"up": 1, "down": -1}[control["direction"]]) \
                if control["kind"] == "positive" else not pointing
        rows.append({"method": f"control_{control['kind']}", "tissue": control["tissue"],
                     "contrast_id": control["contrast"], "metric": f"{control['reaction']}:{target}", "n_cells": 1,
                     "target": np.nan, "decoy_expected": np.nan, "decoy_sd": np.nan, "n_replicates": 0, "ratio": np.nan,
                     "passed": correct, "detail": f"预期 {control.get('direction', 'no change')}，判定 {observed}"})
    for kind in ("positive", "negative"):
        part = [r for r in rows if r["method"] == f"control_{kind}"]
        assessed = [bool(r["passed"]) for r in part if pd.notna(r["passed"])]    # 没有评估的对照不计入符合率
        if part:
            rows.append({"method": f"control_{kind}", "tissue": "", "contrast_id": "", "metric": "concordance",
                         "n_cells": len(part), "target": float(np.mean(assessed)) if assessed else np.nan,
                         "decoy_expected": np.nan, "decoy_sd": np.nan, "n_replicates": 0, "ratio": np.nan,
                         "passed": np.nan, "detail": f"{sum(assessed)}/{len(assessed)} 个对照符合预期"})
    return pd.DataFrame(rows, columns=COLUMNS)


def _axis_values(axes: pd.DataFrame) -> pd.DataFrame:
    return axes.set_index(["flux_id", "tissue", "contrast_id", "axis"])[["state", "direction", "effect"]]


def invariance_checks(members: pd.DataFrame, run: Callable[[pd.DataFrame], pd.DataFrame]) -> pd.DataFrame:
    """在实际数据上检查各轴互不依赖。``run``：成员行 → axis_assessment。"""
    base = _axis_values(run(members))
    rows = []

    def compare(name: str, changed: pd.DataFrame, axes_kept: Iterable[str]) -> None:
        after = _axis_values(changed)
        kept = base[base.index.get_level_values("axis").isin(list(axes_kept))]
        aligned = after.reindex(kept.index)
        same = (kept["state"] == aligned["state"]) & (kept["direction"] == aligned["direction"]) & \
            ((kept["effect"] == aligned["effect"]) | (kept["effect"].isna() & aligned["effect"].isna()))
        violations = int((~same).sum())
        cells = kept.index.droplevel("axis").unique()
        rows.append({"method": "invariance", "tissue": "", "contrast_id": "", "metric": name,
                     "n_cells": len(cells), "target": violations, "decoy_expected": np.nan, "decoy_sd": np.nan,
                     "n_replicates": 0, "ratio": np.nan, "passed": violations == 0,
                     "detail": "、".join(f"{f}/{t}/{c}/{a}" for f, t, c, a in kept.index[~same][:5])})

    flipped = members.copy()
    is_flux = flipped["layer"] == "mfa"
    for column in ("log2fc", "effect", "ci_lower", "ci_upper"):
        if column in flipped.columns:
            flipped.loc[is_flux, column] = -pd.to_numeric(flipped.loc[is_flux, column], errors="coerce")
    compare("flip_flux_sign", run(flipped), ("capacity", "regulation", "pool"))
    for layer in LAYERS:
        if (members["layer"] == layer).any():
            kept = [axis for axis, layers in AXIS_LAYERS.items() if layer not in layers]
            compare(f"remove_{layer}", run(members[members["layer"] != layer]), kept)
    return pd.DataFrame(rows, columns=COLUMNS)


def not_run(reason: str) -> pd.DataFrame:
    return pd.DataFrame([{"method": "permutation", "tissue": "", "contrast_id": "", "metric": "all", "n_cells": 0,
                          "target": np.nan, "decoy_expected": np.nan, "decoy_sd": np.nan, "n_replicates": 0,
                          "ratio": np.nan, "passed": np.nan, "detail": reason}], columns=COLUMNS)


def caption(calibration: pd.DataFrame, tissue: str, contrast: str) -> str:
    """图 A 的图注：该组织、该比较的诱饵对照。"""
    if calibration is None or calibration.empty:
        return ""
    part = calibration[(calibration["method"] == "decoy") & (calibration["tissue"] == tissue)
                       & (calibration["contrast_id"] == contrast)].set_index("metric")
    if part.empty:
        return "Calibration: decoys not run (calibration.decoys = 0)"
    labels = {"L1": "L1", "capacity_changed": "C changed", "capacity_and_pool": "C changed & M concordant"}
    return f"Calibration (decoys ×{int(part['n_replicates'].iloc[0])}): " + "; ".join(
        f"{labels[m]} target {int(part.at[m, 'target'])}, decoy expected {part.at[m, 'decoy_expected']:.1f}"
        for m in METRICS if m in part.index)
