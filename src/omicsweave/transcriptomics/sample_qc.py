#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
转录组通用样本质量控制（``omics/transcriptomics/qc/sample_quality.parquet``）。

面向任意单组织或多组织 RNA-seq 实验，不硬编码任何样本编号、组织名称或基因符号，全部由计数矩阵的数据分布自动推导：

1. ``cross_tissue_contamination``（多组织自动启用）：
   由各组织表达中位数自动提取每个（目标组织 T，异源组织 F）的成对高特异标记基因，计算每个样本中异源组织标记相对同组织中位数的
   对数抬升（``dlog``）、稳健 Z 分数（MAD）与估计的异源 RNA 质量占比（``frac_sum``），识别解剖取材混入或文库交叉污染。
2. ``own_marker_retention``（多组织自动启用）：
   由各组织表达中位数自动提取本组织特异标记，计算每个样本相对本组织中位表达的保留比（``retention``），识别组织身份严重缺失
   （如皮下脂肪混入淋巴结、心肌细胞核转录本严重耗竭等；即使整组同时失效也能检出）。
3. ``transcriptome_correlation``（单组织与多组织均适用）：
   在表达基因上计算样本与同组织中位表达谱的相关系数，以及高变基因上的组内成对中位相关系数，结合 MAD 稳健 Z 分数识别离群样本。
4. ``mitochondrial_fraction``（检测到 ``mt-`` / ``MT-`` 基因时自动启用）：
   统计线粒体基因计数占比，识别极端线粒体富集（如 >70% 读数被线粒体占据）或相对同组织中位数极端偏离的样本。

输出的 ``sample_quality.parquet`` 严格遵守 :mod:`omicsweave.core.reliability` 契约（``sample_id``、``metric``、``value``、
``robust_z``、``flag``、``reason``），并在 ``filter_excluded=True`` 时于差异分析前过滤 ``flag == "exclude"`` 的严重异常样本。
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd

# 与 omicsweave.core.reliability.SAMPLE_QUALITY_COLUMNS 保持一致的契约列
SAMPLE_QUALITY_COLUMNS = ["sample_id", "metric", "value", "robust_z", "flag", "reason"]
MT_GENE_PATTERN = re.compile(r"^(?:mt-|MT-)")

OUTPUT_COLUMNS = [*SAMPLE_QUALITY_COLUMNS, "tissue", "foreign_tissue", "detail"]

DEFAULT_SAMPLE_QC_CONFIG: dict[str, Any] = {
    "enabled": True,
    "filter_excluded": True,
    "min_samples_per_tissue": 4,
    "min_expressed_genes": 50,
    "marker_min_log_cpm": 6.0,
    "marker_min_global_gap": 3.0,
    "marker_min_pair_gap": 5.0,
    "marker_max_target_log_cpm": 3.5,
    "marker_own_min_gap": 4.0,
    "marker_top_n": 25,
    "marker_min_count": 5,
    "contam_flag_log2": 1.5,
    "contam_exclude_log2": 2.5,
    "contam_min_robust_z": 6.0,
    "contam_flag_fraction": 0.0015,
    "contam_exclude_fraction": 0.005,
    "marker_retention_flag": 0.70,
    "marker_retention_flag_z": -2.5,
    "marker_retention_exclude": 0.35,
    "joint_degradation_retention_max": 0.72,
    "joint_degradation_retention_z": -2.0,
    "joint_degradation_var_corr_max": 0.90,
    "joint_degradation_var_corr_z": 2.5,
    "corr_exclude_min": 0.86,
    "corr_exclude_robust_z": 5.0,
    "corr_flag_min": 0.92,
    "corr_flag_robust_z": 4.5,
    "var_corr_top_n": 1000,
    "var_corr_flag_min": 0.90,
    "var_corr_flag_robust_z": 4.5,
    "mt_exclude_max": 0.88,
    "mt_exclude_relative_min": 0.60,
    "mt_exclude_robust_z": 5.0,
    "mt_flag_robust_z": 4.0,
    "mt_flag_min_delta": 0.25,
    "mt_flag_low": False,
}


@dataclass(frozen=True)
class SampleQCResult:
    """通用转录组样本质检产出。"""

    quality: pd.DataFrame          # sample_quality.parquet（逐样本 × 指标，含 flag 与 reason）
    metrics: pd.DataFrame          # sample_qc_metrics.parquet（逐样本汇总诊断表）
    excluded: set[str]             # 需过滤的严重异常样本（flag == "exclude"）
    flagged: set[str]              # 标记但不过滤的样本（flag == "flag"，且不在 excluded 中）
    parameters: dict[str, Any]
    foreign_genes: dict[tuple[str, str], tuple[str, ...]] = field(default_factory=dict)


def resolve_sample_qc_config(raw: dict[str, Any] | None) -> dict[str, Any]:
    """合并用户配置与默认值。"""
    cfg = dict(DEFAULT_SAMPLE_QC_CONFIG)
    if isinstance(raw, dict):
        cfg.update(raw)
    return cfg


def _robust_z(series: pd.Series, floor: float) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    valid = values.dropna()
    if len(valid) < 2:
        return pd.Series(0.0, index=series.index, dtype=float)
    median = float(valid.median())
    mad = float(np.median(np.abs(valid.to_numpy(dtype=float) - median)))
    scale = max(1.4826 * mad, float(floor))
    return (values - median) / scale


def _cpm_and_log_cpm(counts: pd.DataFrame) -> tuple[pd.Series, pd.DataFrame, pd.DataFrame]:
    numeric = counts.apply(pd.to_numeric, errors="coerce").fillna(0.0).clip(lower=0.0)
    lib_size = numeric.sum(axis=0)
    cpm = numeric.div(lib_size.replace(0.0, np.nan), axis=1).fillna(0.0) * 1.0e6
    log_cpm = np.log2(cpm + 1.0)
    return lib_size, cpm, log_cpm


def _derive_own_markers(med_log: pd.DataFrame, cfg: dict[str, Any]) -> dict[str, list[str]]:
    """按组织间表达中位数差自动推导每个组织的特异标记基因。"""
    tissues = list(med_log.columns)
    if len(tissues) < 2:
        return {}
    min_log = float(cfg["marker_min_log_cpm"])
    min_gap = float(cfg["marker_own_min_gap"])
    top_n = int(cfg["marker_top_n"])
    min_count = int(cfg["marker_min_count"])
    markers: dict[str, list[str]] = {}
    for tissue in tissues:
        others = [t for t in tissues if t != tissue]
        gap = med_log[tissue] - med_log[others].max(axis=1)
        selected = gap[(med_log[tissue] >= min_log) & (gap >= min_gap)].sort_values(ascending=False)
        chosen = [str(g) for g in selected.head(top_n).index]
        if len(chosen) >= min_count:
            markers[tissue] = chosen
    return markers


def _derive_pair_marker_sets(
    med_log: pd.DataFrame,
    cfg: dict[str, Any],
) -> tuple[dict[tuple[str, str], list[str]], dict[tuple[str, str], tuple[str, ...]]]:
    """自动推导 (目标组织 T, 异源组织 F) 的成对特异标记（用于样本污染打分与特征级标记）。"""
    tissues = list(med_log.columns)
    if len(tissues) < 2:
        return {}, {}
    min_log = float(cfg["marker_min_log_cpm"])
    min_global_gap = float(cfg["marker_min_global_gap"])
    min_pair_gap = float(cfg["marker_min_pair_gap"])
    max_target = float(cfg["marker_max_target_log_cpm"])
    top_n = int(cfg["marker_top_n"])
    min_count = int(cfg["marker_min_count"])

    other_max = {f: med_log[[o for o in tissues if o != f]].max(axis=1) for f in tissues}
    pairs: dict[tuple[str, str], list[str]] = {}
    foreign_genes: dict[tuple[str, str], tuple[str, ...]] = {}
    for target in tissues:
        for foreign in tissues:
            if foreign == target:
                continue
            gap_all = med_log[foreign] - other_max[foreign]
            gap_pair = med_log[foreign] - med_log[target]
            mask = (
                (med_log[foreign] >= min_log)
                & (gap_all >= min_global_gap)
                & (gap_pair >= min_pair_gap)
                & (med_log[target] <= max_target)
            )
            ordered = [str(g) for g in gap_pair[mask].sort_values(ascending=False).index]
            if len(ordered) >= min_count:
                pairs[(target, foreign)] = ordered[:top_n]
                foreign_genes[(target, foreign)] = tuple(ordered)
    return pairs, foreign_genes


def _derive_pair_markers(med_log: pd.DataFrame, cfg: dict[str, Any]) -> dict[tuple[str, str], list[str]]:
    """自动推导 (目标组织 T, 异源组织 F) 的成对特异标记：在 F 中高表达、全组织特异、在 T 中近零。"""
    return _derive_pair_marker_sets(med_log, cfg)[0]


def evaluate_sample_quality(
    counts: pd.DataFrame,
    tissues: pd.Series | dict[str, str],
    config: dict[str, Any] | None = None,
    *,
    manual_flags: dict[str, str] | None = None,
    groups: pd.Series | dict[str, str] | None = None,
) -> SampleQCResult:
    """对多组织或单组织计数矩阵执行通用样本质检。

    Parameters
    ----------
    counts : pd.DataFrame
        基因 × 样本的计数矩阵（列为 ``sample_id``）。
    tissues : pd.Series or dict[str, str]
        ``sample_id -> tissue`` 映射。
    config : dict, optional
        覆盖 :data:`DEFAULT_SAMPLE_QC_CONFIG` 的参数。
    manual_flags : dict[str, str], optional
        来自 ``samples.tsv`` 的人工 ``qc_flag``（仅用于合并 ``excluded`` 集合，不重复写入质量表）。
    groups : pd.Series or dict[str, str], optional
        ``sample_id -> 实验分组标签`` 映射；提供时，高变基因成对相关离群判定仅对同实验组内部
        同样偏离的样本生效，避免强药物或饮食表型重塑被误判为技术离群。
    """
    cfg = resolve_sample_qc_config(config)
    sample_ids = [str(c) for c in counts.columns]
    tissue_map = {s: str(tissues[s]) for s in sample_ids if s in tissues}
    groups_dict = groups.to_dict() if isinstance(groups, pd.Series) else (groups or {})
    group_map = {str(s): str(groups_dict[s]).strip() for s in sample_ids if s in groups_dict and str(groups_dict[s]).strip()}
    manual = {str(k): str(v).strip().lower() for k, v in (manual_flags or {}).items() if str(v).strip()}

    if not cfg.get("enabled", True) or counts.empty or not sample_ids:
        empty_q = pd.DataFrame(columns=OUTPUT_COLUMNS)
        empty_m = pd.DataFrame(columns=["sample_id", "tissue", "qc_status"])
        excluded = {s for s in sample_ids if manual.get(s) == "exclude"}
        flagged = {s for s in sample_ids if manual.get(s) and manual.get(s) != "exclude"} - excluded
        return SampleQCResult(empty_q, empty_m, excluded, flagged, cfg)

    lib_size, cpm, log_cpm = _cpm_and_log_cpm(counts)
    by_tissue: dict[str, list[str]] = {}
    for sid in sample_ids:
        t = tissue_map.get(sid, "")
        if t:
            by_tissue.setdefault(t, []).append(sid)

    min_n = int(cfg["min_samples_per_tissue"])
    min_genes = int(cfg["min_expressed_genes"])

    # 仅使用样本数达标的组织构建跨组织中位参考，避免单样本组织污染标记基因
    eligible_tissues = [t for t, sids in by_tissue.items() if len(sids) >= min_n]
    med_log = pd.DataFrame({t: log_cpm[by_tissue[t]].median(axis=1) for t in eligible_tissues}) if eligible_tissues else pd.DataFrame(index=log_cpm.index)
    med_cpm = pd.DataFrame({t: cpm[by_tissue[t]].median(axis=1) for t in eligible_tissues}) if eligible_tissues else pd.DataFrame(index=cpm.index)

    own_markers = _derive_own_markers(med_log, cfg) if len(eligible_tissues) >= 2 else {}
    pair_markers, foreign_genes = _derive_pair_marker_sets(med_log, cfg) if len(eligible_tissues) >= 2 else ({}, {})

    mt_mask = counts.index.astype(str).map(lambda g: bool(MT_GENE_PATTERN.match(g))).to_numpy(dtype=bool)
    has_mt = bool(mt_mask.any())
    mt_frac_all = (
        counts.iloc[mt_mask].apply(pd.to_numeric, errors="coerce").fillna(0.0).sum(axis=0)
        / lib_size.replace(0.0, np.nan)
    ).fillna(0.0) if has_mt else pd.Series(np.nan, index=counts.columns, dtype=float)

    quality_rows: list[dict[str, Any]] = []
    metric_summary: dict[str, dict[str, Any]] = {
        s: {
            "sample_id": s,
            "tissue": tissue_map.get(s, ""),
            "library_size": float(lib_size[s]),
            "detected_genes": int((counts[s] > 0).sum()),
            "mt_fraction": float(mt_frac_all[s]) if has_mt else np.nan,
            "corr_median": np.nan,
            "corr_median_rz": np.nan,
            "var_pairwise_corr": np.nan,
            "var_pairwise_corr_rz": np.nan,
            "own_marker_retention": np.nan,
            "own_marker_retention_rz": np.nan,
            "top_foreign_tissue": "",
            "foreign_log2_enrichment": np.nan,
            "foreign_enrichment_rz": np.nan,
            "foreign_rna_fraction": np.nan,
            "qc_status": "ok",
            "qc_reasons": "",
        }
        for s in sample_ids
    }

    for tissue, sids in by_tissue.items():
        if len(sids) < min_n:
            continue

        # ── 1. 转录组整体与高变基因相关系数 ─────────────────────────────────
        within_ok: dict[str, bool] = {}
        ref_log = log_cpm[sids].median(axis=1)
        expr_genes = ref_log.index[ref_log >= 1.0]
        if len(expr_genes) >= min_genes:
            sub_expr = log_cpm.loc[expr_genes, sids]
            ref_expr = ref_log.loc[expr_genes]
            r_all = pd.Series({s: float(sub_expr[s].corr(ref_expr)) for s in sids}, dtype=float)
            dist_all = 1.0 - r_all
            rz_all = _robust_z(dist_all, floor=0.005)

            # 高变基因成对中位相关系数（能捕捉 2 只同组离群但彼此相似的样本偏离大群体）
            top_n_var = min(int(cfg["var_corr_top_n"]), len(expr_genes))
            top_var = sub_expr.var(axis=1).nlargest(top_n_var).index
            corr_mat = sub_expr.loc[top_var].corr()
            r_var_pw = pd.Series(
                {s: float(corr_mat.loc[s].drop(s).median()) for s in sids},
                dtype=float,
            )
            dist_var = 1.0 - r_var_pw
            rz_var = _robust_z(dist_var, floor=0.005)
            for s in sids:
                grp = group_map.get(s)
                peers = [p for p in sids if p != s and group_map.get(p) == grp] if grp else []
                if peers:
                    r_within = float(corr_mat.loc[s, peers].median())
                    within_ok[s] = bool(np.isfinite(r_within) and r_within >= float(cfg["var_corr_flag_min"]))
                else:
                    within_ok[s] = False

            for s in sids:
                metric_summary[s]["corr_median"] = float(r_all[s])
                metric_summary[s]["corr_median_rz"] = float(rz_all[s])
                metric_summary[s]["var_pairwise_corr"] = float(r_var_pw[s])
                metric_summary[s]["var_pairwise_corr_rz"] = float(rz_var[s])

                flag, reason, value, rz = "", "", float(r_all[s]), float(rz_all[s])
                detail = f"r_all={r_all[s]:.4f} (z={rz_all[s]:+.1f}); r_var_pw={r_var_pw[s]:.4f} (z={rz_var[s]:+.1f})"
                if np.isfinite(r_all[s]) and r_all[s] < float(cfg["corr_exclude_min"]) and rz_all[s] >= float(cfg["corr_exclude_robust_z"]):
                    flag = "exclude"
                    reason = f"表达谱严重离群（r={r_all[s]:.3f}, z={rz_all[s]:.1f}）"
                elif np.isfinite(r_all[s]) and r_all[s] < float(cfg["corr_flag_min"]) and rz_all[s] >= float(cfg["corr_flag_robust_z"]):
                    flag = "flag"
                    reason = f"表达谱偏离同组织中位谱（r={r_all[s]:.3f}, z={rz_all[s]:.1f}）"
                elif (
                    not within_ok.get(s, False)
                    and np.isfinite(r_var_pw[s])
                    and r_var_pw[s] < float(cfg["var_corr_flag_min"])
                    and rz_var[s] >= float(cfg["var_corr_flag_robust_z"])
                ):
                    flag = "flag"
                    value, rz = float(r_var_pw[s]), float(rz_var[s])
                    reason = f"高变基因表达谱偏离同组织群体（r_var={r_var_pw[s]:.3f}, z={rz_var[s]:.1f}）"
                quality_rows.append({
                    "sample_id": s, "metric": "transcriptome_correlation",
                    "value": value, "robust_z": rz, "flag": flag, "reason": reason,
                    "tissue": tissue, "foreign_tissue": "", "detail": detail,
                })

        # ── 2. 线粒体读数占比 ─────────────────────────────────────────────
        if has_mt:
            mt_t = mt_frac_all[sids]
            med_mt = float(mt_t.median())
            rz_mt = _robust_z(mt_t, floor=0.02)
            for s in sids:
                val = float(mt_t[s])
                rz = float(rz_mt[s])
                delta = val - med_mt
                flag, reason = "", ""
                if val >= float(cfg["mt_exclude_max"]):
                    flag = "exclude"
                    reason = f"线粒体读数占比极端偏高（{val:.1%} >= {float(cfg['mt_exclude_max']):.0%}）"
                elif (
                    val >= float(cfg["mt_exclude_relative_min"])
                    and rz >= float(cfg["mt_exclude_robust_z"])
                    and delta >= float(cfg["mt_flag_min_delta"])
                ):
                    flag = "exclude"
                    reason = f"线粒体读数占比相对同组织中位数（{med_mt:.1%}）极端偏高（{val:.1%}, z={rz:+.1f}）"
                elif (
                    rz >= float(cfg["mt_flag_robust_z"]) and delta >= float(cfg["mt_flag_min_delta"])
                ) or (
                    bool(cfg.get("mt_flag_low", False))
                    and rz <= -float(cfg["mt_flag_robust_z"])
                    and delta <= -float(cfg["mt_flag_min_delta"])
                ):
                    flag = "flag"
                    direction = "偏高" if delta > 0 else "偏低"
                    reason = f"线粒体读数占比相对同组织中位数（{med_mt:.1%}）异常{direction}（{val:.1%}, z={rz:+.1f}）"
                quality_rows.append({
                    "sample_id": s, "metric": "mitochondrial_fraction",
                    "value": val, "robust_z": rz, "flag": flag, "reason": reason,
                    "tissue": tissue, "foreign_tissue": "",
                    "detail": f"mt={val:.1%} (tissue_median={med_mt:.1%}, z={rz:+.1f})",
                })

        # ── 3. 本组织特异标记保留度（多组织） ──────────────────────────────
        own_m = own_markers.get(tissue, [])
        if own_m:
            denom = med_cpm.loc[own_m, tissue].replace(0.0, np.nan)
            ret = pd.Series({s: float((cpm.loc[own_m, s] / denom).median()) for s in sids}, dtype=float)
            ret_log = np.log2(ret.clip(lower=1.0e-4))
            rz_ret = _robust_z(ret_log, floor=0.15)
            for s in sids:
                val, rz = float(ret[s]), float(rz_ret[s])
                metric_summary[s]["own_marker_retention"] = val
                metric_summary[s]["own_marker_retention_rz"] = rz
                flag, reason = "", ""
                if np.isfinite(val) and val < float(cfg["marker_retention_exclude"]):
                    flag = "exclude"
                    reason = f"本组织特异标记严重缺失（保留比 {val:.2f} < {float(cfg['marker_retention_exclude']):.2f}）"
                elif np.isfinite(val) and val < float(cfg["marker_retention_flag"]) and rz <= float(cfg["marker_retention_flag_z"]):
                    flag = "flag"
                    reason = f"本组织特异标记明显偏低（保留比 {val:.2f}, z={rz:+.1f}）"
                else:
                    r_var_val = metric_summary[s]["var_pairwise_corr"]
                    r_var_z = metric_summary[s]["var_pairwise_corr_rz"]
                    if (
                        not within_ok.get(s, False)
                        and np.isfinite(val) and np.isfinite(r_var_val) and np.isfinite(r_var_z)
                        and val < float(cfg["joint_degradation_retention_max"])
                        and rz <= float(cfg["joint_degradation_retention_z"])
                        and r_var_val < float(cfg["joint_degradation_var_corr_max"])
                        and r_var_z >= float(cfg["joint_degradation_var_corr_z"])
                    ):
                        flag = "flag"
                        reason = (
                            f"本组织特异标记偏低且高变表达谱偏离群体"
                            f"（保留比 {val:.2f}, z={rz:+.1f}；r_var={r_var_val:.3f}, z={r_var_z:+.1f}）"
                        )
                quality_rows.append({
                    "sample_id": s, "metric": "own_marker_retention",
                    "value": val, "robust_z": rz, "flag": flag, "reason": reason,
                    "tissue": tissue, "foreign_tissue": "",
                    "detail": f"retention={val:.3f} (z={rz:+.1f}, n_markers={len(own_m)})",
                })

        # ── 4. 跨组织 RNA 污染检测（多组织） ───────────────────────────────
        foreign_candidates = [f for f in eligible_tissues if f != tissue and (tissue, f) in pair_markers]
        if foreign_candidates:
            best_by_sample: dict[str, tuple[int, float, float, float, str, str, str]] = {}
            for foreign in foreign_candidates:
                fm = pair_markers[(tissue, foreign)]
                ref_f_cpm = med_cpm.loc[fm, foreign]
                ref_f_sum = float(ref_f_cpm.sum())
                if ref_f_sum <= 0.0:
                    continue
                ref_t_log = med_log.loc[fm, tissue]
                dlog = pd.Series({s: float((log_cpm.loc[fm, s] - ref_t_log).median()) for s in sids}, dtype=float)
                frac_sum = pd.Series({s: float(cpm.loc[fm, s].sum() / ref_f_sum) for s in sids}, dtype=float)
                rz_dl = _robust_z(dlog, floor=0.15)
                for s in sids:
                    dl, fr, rz = float(dlog[s]), float(frac_sum[s]), float(rz_dl[s])
                    flag, reason, sev = "", "", 0
                    if (
                        dl >= float(cfg["contam_exclude_log2"])
                        and rz >= float(cfg["contam_min_robust_z"])
                        and fr >= float(cfg["contam_exclude_fraction"])
                    ):
                        flag, sev = "exclude", 2
                        reason = f"疑似严重混入 {foreign} 组织 RNA（异源标记抬升 +{dl:.2f} log2，占比 {fr:.2%}，z={rz:.1f}）"
                    elif (
                        dl >= float(cfg["contam_flag_log2"])
                        and rz >= float(cfg["contam_min_robust_z"])
                        and fr >= float(cfg["contam_flag_fraction"])
                    ):
                        flag, sev = "flag", 1
                        reason = f"疑似混入 {foreign} 组织 RNA（异源标记抬升 +{dl:.2f} log2，占比 {fr:.2%}，z={rz:.1f}）"
                    prev = best_by_sample.get(s)
                    score = (sev, fr, dl)
                    if prev is None or score > (prev[0], prev[1], prev[2]):
                        best_by_sample[s] = (sev, fr, dl, rz, foreign, flag, reason)
                    # 若触发 flag/exclude，单独记录每一对被触发的异源组织以便溯源
                    if flag:
                        quality_rows.append({
                            "sample_id": s, "metric": "cross_tissue_contamination",
                            "value": fr, "robust_z": rz, "flag": flag, "reason": reason,
                            "tissue": tissue, "foreign_tissue": foreign,
                            "detail": f"foreign={foreign}; dlog={dl:+.2f}; frac={fr:.4f}; z={rz:+.1f}; n_markers={len(fm)}",
                        })
            for s, (sev, fr, dl, rz, foreign, flag, reason) in best_by_sample.items():
                metric_summary[s]["top_foreign_tissue"] = foreign
                metric_summary[s]["foreign_log2_enrichment"] = dl
                metric_summary[s]["foreign_enrichment_rz"] = rz
                metric_summary[s]["foreign_rna_fraction"] = fr
                if sev == 0:
                    quality_rows.append({
                        "sample_id": s, "metric": "cross_tissue_contamination",
                        "value": fr, "robust_z": rz, "flag": "", "reason": "",
                        "tissue": tissue, "foreign_tissue": foreign,
                        "detail": f"top_foreign={foreign}; dlog={dl:+.2f}; frac={fr:.4f}; z={rz:+.1f}",
                    })

    quality_df = pd.DataFrame(quality_rows, columns=OUTPUT_COLUMNS)
    auto_excluded = set(quality_df.loc[quality_df["flag"] == "exclude", "sample_id"].astype(str))
    auto_flagged = set(quality_df.loc[quality_df["flag"] == "flag", "sample_id"].astype(str)) - auto_excluded
    excluded = auto_excluded | {s for s in sample_ids if manual.get(s) == "exclude"}
    flagged = (auto_flagged | {s for s in sample_ids if manual.get(s) and manual.get(s) != "exclude"}) - excluded

    reasons_by_sample: dict[str, list[str]] = {}
    for row in quality_df[quality_df["flag"].astype(str).ne("")].itertuples(index=False):
        reasons_by_sample.setdefault(str(row.sample_id), []).append(str(row.reason))
    for s in sample_ids:
        if s in excluded:
            metric_summary[s]["qc_status"] = "exclude"
        elif s in flagged:
            metric_summary[s]["qc_status"] = "flag"
        metric_summary[s]["qc_reasons"] = "；".join(dict.fromkeys(reasons_by_sample.get(s, [])))

    metrics_df = pd.DataFrame([metric_summary[s] for s in sample_ids])
    return SampleQCResult(quality_df, metrics_df, excluded, flagged, cfg, foreign_genes)
