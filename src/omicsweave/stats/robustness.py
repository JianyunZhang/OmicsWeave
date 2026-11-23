#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
入选特征的稳健性（E5，``omics/<组学>/tables/robustness.parquet``）。

- 逐只去掉动物：对每个两两比较中入选（``selected``）的特征，在 log2 尺度上（计数为 log2 CPM）用组均值之差估计效应，
  依次去掉一个样本重新估计。去掉某一只后效应反号，或缩小到原来的 ``min_ratio``（默认 0.5）以下，记为
  ``single_subject_driven``，并给出影响最大的样本；
- 符号错误率与夸大倍数（Gelman 与 Carlin 2014）：以正态—正态收缩后的效应作为可信的真实效应，
  结合差异结果的标准误，估计该比较下的功效、Type S 与 Type M。

去掉的是样本（每只动物在一个组织中一个样本）。只剩 1 个样本的组不能再去掉，``n_left_out`` 记录实际重算的次数。
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

from omicsweave.contract.dataset import read_dataset
from omicsweave.contract.tables import effect_standard_errors
from omicsweave.core.config import Experiment
from omicsweave.core.model import Omics
from omicsweave.core.statistics.shrinkage import normal_shrinkage, retrodesign

from .data import analysis_values, bundle_datasets, pairwise_groups

FLAG = "single_subject_driven"
DEFAULT_MIN_RATIO = 0.5
TRACER_COVERAGE = Path("omics") / Omics.MFA.value / "tables" / "tracer_coverage.parquet"
COLUMNS = ["contrast_id", "tissue", "omics", "feature_id", "n_numerator", "n_denominator", "n_left_out",
           "full_effect", "loo_min_ratio", "loo_sign_flips", "driver_sample", FLAG, "se", "plausible_effect",
           "power", "type_s", "exaggeration"]


def leave_one_out(first: pd.DataFrame, second: pd.DataFrame) -> tuple[pd.Series, pd.DataFrame]:
    """(全样本的组均值之差, 去掉每个样本后的差：特征 × 被去掉的样本)。"""
    def loo_means(values: pd.DataFrame) -> pd.DataFrame:
        total, count = values.sum(axis=1, min_count=1), values.notna().sum(axis=1)
        out = {}
        for sample in values.columns:
            present = values[sample].notna()
            remaining = count - present.astype(int)
            out[sample] = ((total - values[sample].fillna(0.0)) / remaining.where(remaining > 0))
        return pd.DataFrame(out, index=values.index)

    mean1, mean2 = first.mean(axis=1), second.mean(axis=1)
    loo = pd.concat([loo_means(first).sub(mean2, axis=0), -loo_means(second).sub(mean1, axis=0)], axis=1)
    return mean1 - mean2, loo


def _mfa_protocols(experiment: Experiment, columns) -> dict[str, str] | None:
    """多示踪剂 MFA 的样本 → 示踪剂（protocol）映射；单示踪剂或未标注 protocol 时返回 None。"""
    table = experiment.samples
    if "protocol" not in table.columns:
        return None
    mfa = table[(table["omics"] == Omics.MFA.value) & table["sample_id"].isin(list(columns))]
    if mfa.empty:
        return None
    mapping = {
        str(sid): str(proto).strip()
        for sid, proto in zip(mfa["sample_id"], mfa["protocol"].fillna(""))
        if str(proto).strip()
    }
    return mapping if len(set(mapping.values())) >= 2 else None


def _mfa_informative_tracers(run_dir: Path, tissue: str, selected: pd.DataFrame) -> dict[str, set[str]] | None:
    """各入选 MFA 反应参与合并的可辨识示踪剂集合；优先读差异表的 ``tracers_informative``，再读 ``tracer_coverage.parquet``。"""
    if "tracers_informative" in selected.columns:
        out: dict[str, set[str]] = {}
        for fid, text in zip(selected["feature_id"].astype(str), selected["tracers_informative"]):
            if pd.isna(text):
                out[fid] = set()
            else:
                out[fid] = {t.strip() for t in str(text).split(";") if t.strip()}
        return out
    cov_path = Path(run_dir) / TRACER_COVERAGE
    if not cov_path.is_file():
        return None
    cov = pd.read_parquet(cov_path)
    part = cov[cov["tissue"].astype(str) == str(tissue)]
    if part.empty:
        return None
    out = {str(rxn): set() for rxn in part["reaction"].astype(str).unique()}
    for r in part[part["informative"].fillna(False).astype(bool)].itertuples(index=False):
        out.setdefault(str(r.reaction), set()).add(str(r.tracer).strip())
    return out


def _stratified_mfa_loo_one(
    row_values: pd.Series,
    num: list[str],
    den: list[str],
    protocols: dict[str, str],
    allowed: set[str] | None,
    min_ratio: float,
) -> tuple[int, int, int, float, float, int, str, bool]:
    """单条反应在多示踪剂 MFA 下的分层留一稳健性（仅纳入可辨识且每组 ≥ 2 只的示踪剂）。

    - 单一有效示踪剂（k = 1）：直接在该示踪剂内做组均值差留一法，不混入其他不可辨识示踪剂样本；
    - 多有效示踪剂（k ≥ 2）：以全样本各示踪剂的 Hedges' g 尺度因子 c_t = J(df_t)/s_{p,t} 与逆方差权重
      w_t = 1/Var(g_t) 固定标准化尺度，逐只去除任一示踪剂中的 1 只动物重算合并效应 E_{-i}，求比值 E_{-i}/E_full。
      固定全样本 s_{p,t} 可避免 n=2 示踪剂去 1 只后自由度 df=1 时 Hedges 校正归零的奇异性，且 k=1 时退化为精确均值差比值。
    """
    by_tracer: dict[str, tuple[list[str], list[float], list[str], list[float]]] = {}
    for side_idx, samples in ((0, num), (1, den)):
        for sid in samples:
            tracer = protocols.get(str(sid), "")
            if not tracer or (allowed is not None and tracer not in allowed):
                continue
            val = float(row_values[sid]) if pd.notna(row_values[sid]) else np.nan
            if not np.isfinite(val):
                continue
            entry = by_tracer.setdefault(tracer, ([], [], [], []))
            if side_idx == 0:
                entry[0].append(str(sid))
                entry[1].append(val)
            else:
                entry[2].append(str(sid))
                entry[3].append(val)

    strata = [
        (t, a_ids, np.asarray(a_vals, float), b_ids, np.asarray(b_vals, float))
        for t, (a_ids, a_vals, b_ids, b_vals) in sorted(by_tracer.items())
        if len(a_vals) >= 2 and len(b_vals) >= 2
    ]
    n_num = sum(len(s[2]) for s in strata)
    n_den = sum(len(s[4]) for s in strata)
    if not strata:
        return n_num, n_den, 0, np.nan, np.nan, 0, "", False

    if len(strata) == 1:
        _, a_ids, a, b_ids, b = strata[0]
        full = float(a.mean() - b.mean())
        if not np.isfinite(full) or full == 0.0:
            return n_num, n_den, 0, full, np.nan, 0, "", False
        loo_sids = list(a_ids) + list(b_ids)
        loo_vals = np.array(
            [float(np.delete(a, i).mean() - b.mean()) for i in range(len(a))]
            + [float(a.mean() - np.delete(b, i).mean()) for i in range(len(b))],
            dtype=float,
        )
        ratios = loo_vals / full
    else:
        weighted = []
        for t, a_ids, a, b_ids, b in strata:
            n1, n2 = len(a), len(b)
            df = n1 + n2 - 2
            sp2 = ((n1 - 1) * a.var(ddof=1) + (n2 - 1) * b.var(ddof=1)) / df
            if not np.isfinite(sp2) or sp2 <= 0.0:
                continue
            c = (1.0 - 3.0 / (4.0 * df - 1.0)) / math.sqrt(sp2)
            delta = float(a.mean() - b.mean())
            g = c * delta
            w = 1.0 / ((n1 + n2) / (n1 * n2) + g * g / (2.0 * (n1 + n2)))
            if np.isfinite(g) and np.isfinite(w) and w > 0.0:
                weighted.append((t, a_ids, a, b_ids, b, c, w, g))
        if not weighted:
            # 全部示踪剂组内零方差时的回退：各示踪剂等权均值差
            deltas = np.array([float(s[2].mean() - s[4].mean()) for s in strata], dtype=float)
            full = float(deltas.mean())
            if not np.isfinite(full) or full == 0.0:
                return n_num, n_den, 0, full, np.nan, 0, "", False
            loo_sids, loo_list = [], []
            for idx, (_, a_ids, a, b_ids, b) in enumerate(strata):
                other = float(deltas.sum() - deltas[idx])
                for i, sid in enumerate(a_ids):
                    loo_sids.append(sid)
                    loo_list.append((other + float(np.delete(a, i).mean() - b.mean())) / len(strata))
                for i, sid in enumerate(b_ids):
                    loo_sids.append(sid)
                    loo_list.append((other + float(a.mean() - np.delete(b, i).mean())) / len(strata))
            ratios = np.asarray(loo_list, dtype=float) / full
        else:
            w_sum = sum(item[6] for item in weighted)
            full = float(sum(item[6] * item[7] for item in weighted) / w_sum)
            if not np.isfinite(full) or full == 0.0:
                return n_num, n_den, 0, full, np.nan, 0, "", False
            loo_sids, loo_list = [], []
            for idx, (_, a_ids, a, b_ids, b, c, w, _) in enumerate(weighted):
                other = sum(item[6] * item[7] for j, item in enumerate(weighted) if j != idx)
                for i, sid in enumerate(a_ids):
                    d_loo = float(np.delete(a, i).mean() - b.mean())
                    loo_sids.append(sid)
                    loo_list.append((other + w * c * d_loo) / w_sum)
                for i, sid in enumerate(b_ids):
                    d_loo = float(a.mean() - np.delete(b, i).mean())
                    loo_sids.append(sid)
                    loo_list.append((other + w * c * d_loo) / w_sum)
            ratios = np.asarray(loo_list, dtype=float) / full

    finite = np.isfinite(ratios)
    if not finite.any():
        return n_num, n_den, 0, full, np.nan, 0, "", False
    min_idx = int(np.argmin(np.where(finite, ratios, np.inf)))
    min_r = float(ratios[min_idx])
    flips = int((ratios[finite] <= 0.0).sum())
    flagged = bool(min_r < min_ratio)
    driver = loo_sids[min_idx] if flagged else ""
    return n_num, n_den, int(finite.sum()), full, min_r, flips, driver, flagged


def robustness(experiment: Experiment, run_dir, tables: dict[str, pd.DataFrame],
               derived: dict[str, pd.DataFrame], *, alpha: float = 0.05,
               min_ratio: float = DEFAULT_MIN_RATIO) -> dict[str, pd.DataFrame]:
    """返回 {组学: 稳健性表}；``tables`` / ``derived`` 为测得特征与派生特征的差异结果。"""
    contrasts = {c.id: c for c in experiment.config.design.contrasts}
    out: dict[str, list[pd.DataFrame]] = {}
    for omics, kind, tissue, path in bundle_datasets(run_dir):
        source = (derived if kind == "derived" else tables).get(omics)
        if source is None or source.empty:
            continue
        values = analysis_values(read_dataset(path))
        protocols = _mfa_protocols(experiment, values.columns) if omics == Omics.MFA.value and kind != "derived" else None
        rows_tissue = source[source["tissue"].astype(str) == tissue]
        for contrast_id, rows in rows_tissue.groupby("contrast_id", sort=True):
            contrast = contrasts.get(str(contrast_id))
            if contrast is None:
                continue
            selected = rows[rows["selected"].fillna(False).astype(bool)]
            selected = selected[selected["feature_id"].astype(str).isin(values.index)]
            if selected.empty:
                continue
            groups = pairwise_groups(experiment, contrast, omics, values.columns)
            num = [s for s, g in groups.items() if g == contrast.numerator]
            den = [s for s, g in groups.items() if g == contrast.denominator]
            if not num or not den:
                continue
            ids = selected["feature_id"].astype(str).to_numpy()
            se_all, _ = effect_standard_errors(rows)
            plausible = pd.Series(normal_shrinkage(rows["log2fc"], se_all), index=rows.index)
            se = se_all.loc[selected.index].to_numpy()
            d = plausible.loc[selected.index].to_numpy()
            power, type_s, exaggeration = retrodesign(d, se, alpha)

            if protocols is not None:
                informative_map = _mfa_informative_tracers(Path(run_dir), tissue, selected)
                strat_rows = [
                    _stratified_mfa_loo_one(
                        values.loc[fid], num, den, protocols,
                        informative_map.get(fid) if informative_map is not None else None,
                        min_ratio,
                    )
                    for fid in ids
                ]
                n_num_arr, n_den_arr, n_lo_arr, full_arr, min_r_arr, flips_arr, driver_arr, flagged_arr = zip(*strat_rows)
                out.setdefault(omics, []).append(pd.DataFrame({
                    "contrast_id": str(contrast_id), "tissue": tissue, "omics": omics, "feature_id": ids,
                    "n_numerator": np.asarray(n_num_arr, dtype=int),
                    "n_denominator": np.asarray(n_den_arr, dtype=int),
                    "n_left_out": np.asarray(n_lo_arr, dtype=int),
                    "full_effect": np.asarray(full_arr, dtype=float),
                    "loo_min_ratio": np.asarray(min_r_arr, dtype=float),
                    "loo_sign_flips": np.asarray(flips_arr, dtype=int),
                    "driver_sample": np.asarray(driver_arr, dtype=str),
                    FLAG: np.asarray(flagged_arr, dtype=bool), "se": se, "plausible_effect": d,
                    "power": power, "type_s": type_s, "exaggeration": exaggeration,
                }))
                continue

            full, loo = leave_one_out(values.loc[ids, num], values.loc[ids, den])
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = loo.div(full, axis=0)
            usable = ratio.notna()
            min_ratio_value = ratio.min(axis=1)
            # 逐行取比值最小的样本；整行都无法重算（如该特征有缺失值）时为空
            filled = ratio.to_numpy(dtype=float, na_value=np.nan)
            any_usable = usable.any(axis=1).to_numpy()
            positions = np.argmin(np.where(np.isnan(filled), np.inf, filled), axis=1)
            driver = pd.Series(np.where(any_usable, ratio.columns.to_numpy()[positions], ""), index=ratio.index)
            flagged = (min_ratio_value < min_ratio).fillna(False)
            out.setdefault(omics, []).append(pd.DataFrame({
                "contrast_id": str(contrast_id), "tissue": tissue, "omics": omics, "feature_id": ids,
                "n_numerator": len(num), "n_denominator": len(den), "n_left_out": usable.sum(axis=1).to_numpy(),
                "full_effect": full.to_numpy(), "loo_min_ratio": min_ratio_value.to_numpy(),
                "loo_sign_flips": (ratio <= 0).sum(axis=1).to_numpy(),
                "driver_sample": driver.where(flagged, "").fillna("").astype(str).to_numpy(),
                FLAG: flagged.to_numpy(dtype=bool), "se": se, "plausible_effect": d, "power": power,
                "type_s": type_s, "exaggeration": exaggeration,
            }))
    return {omics: pd.concat(parts, ignore_index=True)[COLUMNS] for omics, parts in out.items()}


def flag_rows(table: pd.DataFrame, robust: pd.DataFrame | None) -> pd.DataFrame:
    """在差异结果中由单只动物驱动的行的 ``flags`` 中加上 ``single_subject_driven``。"""
    if robust is None or robust.empty or table.empty:
        return table
    hits = robust[robust[FLAG]]
    keys = set(zip(hits["contrast_id"], hits["tissue"], hits["feature_id"]))
    mask = [(c, t, f) in keys for c, t, f in zip(table["contrast_id"].astype(str), table["tissue"].astype(str),
                                                   table["feature_id"].astype(str))]
    if not any(mask):
        return table
    flags = table["flags"].fillna("").astype(str) if "flags" in table.columns else pd.Series("", index=table.index)
    return table.assign(flags=[";".join(p for p in (f, FLAG) if p) if m and FLAG not in f.split(";") else f
                               for f, m in zip(flags, mask)])
