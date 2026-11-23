#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MFA 通用样本与比较级质检（``omics/mfa/qc/sample_quality.parquet``、``contrast_scale_qc.parquet``）。

面向任意单示踪剂或多示踪剂 13C-MFA 实验，不硬编码任何样本 ID、分组名或反应名：

1. 样本级质量控制（``sample_quality.parquet``，遵守 :mod:`omicsweave.core.reliability` 契约）：
   - ``fitting_loss``：识别拟合残差（loss）显著高于同示踪剂队列中位数的离群样本；
   - ``flux_scale_factor``：在全组织活跃反应上计算每只动物相对同示踪剂队列中位通量的 ``med_log2_scale`` 与整鼠总绝对通量比
     ``log2_total_flux_ratio``，识别单只动物整体通量尺度异常放大或塌缩；
   - ``upper_bound_saturation``：统计每只动物贴近推断盒约束上界的反应数与占比，识别优化器大面积撞顶的样本。
2. 比较级尺度漂移控制（``contrast_scale_qc.parquet``）：
   - 对每个（比较, 示踪剂），比较分子组与分母组的整鼠总通量对数均值差（``delta_log2_total_flux``）与全反应中位尺度差
     （``delta_median_log2_scale``），当两者同向且同时超过阈值时标为 ``flux_scale_drift``（例如跨时长批次 3W vs 1W 导致的
     系统性通量放大），并将该标记写入统一差异表 ``flags``，防止尺度漂移伪造成通量 L2/L3 结论。
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

# 与 omicsweave.core.reliability.SAMPLE_QUALITY_COLUMNS 保持一致
SAMPLE_QUALITY_COLUMNS = ["sample_id", "metric", "value", "robust_z", "flag", "reason"]
OUTPUT_QUALITY_COLUMNS = [*SAMPLE_QUALITY_COLUMNS, "tracer", "detail"]
CONTRAST_SCALE_COLUMNS = [
    "contrast_id", "tracer", "n_numerator", "n_denominator", "n_active_reactions",
    "delta_log2_total_flux", "delta_median_log2_scale",
    "mean_ub_numerator", "mean_ub_denominator", "ub_count_diff",
    "scale_drift", "ub_asymmetry", "flag", "reason",
]

DEFAULT_MFA_QC_CONFIG: dict[str, Any] = {
    "enabled": True,
    "min_samples_per_tracer": 4,
    "active_reaction_min_flux": 5.0,
    "upper_bound_rel_tolerance": 0.02,
    "loss_flag_robust_z": 4.0,
    "loss_flag_min_ratio": 1.5,
    "loss_exclude_robust_z": 8.0,
    "loss_exclude_min_ratio": 3.0,
    "scale_flag_log2": 0.60,
    "scale_flag_robust_z": 3.0,
    "total_scale_flag_log2": 0.65,
    "total_scale_flag_robust_z": 2.5,
    "ub_flag_min_count": 4,
    "ub_flag_min_fraction": 0.008,
    "ub_flag_robust_z": 3.0,
    "contrast_scale_min_reactions": 10,
    "contrast_scale_total_log2": 0.30,
    "contrast_scale_median_log2": 0.15,
    "contrast_ub_asym_min_diff": 1.0,
}


@dataclass(frozen=True)
class MFAQualityResult:
    """MFA 通用质检结果。"""

    sample_quality: pd.DataFrame
    sample_metrics: pd.DataFrame
    contrast_scale: pd.DataFrame
    flagged_samples: set[str]
    excluded_samples: set[str]
    drift_contrasts: set[tuple[str, str]]   # {(contrast_id, tracer)}
    parameters: dict[str, Any]


def resolve_mfa_qc_config(raw: dict[str, Any] | None) -> dict[str, Any]:
    cfg = dict(DEFAULT_MFA_QC_CONFIG)
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


def _prepare_animals(animals: pd.DataFrame) -> pd.DataFrame:
    df = animals.copy()
    if "tracer" not in df.columns:
        df["tracer"] = "single"
    df["tracer"] = df["tracer"].fillna("single").astype(str)
    if "flux_id" not in df.columns:
        df["flux_id"] = df["tissue"].astype(str) + "--" + df["reaction"].astype(str)
    df["abs_flux"] = pd.to_numeric(df["net_flux"], errors="coerce").abs()
    return df


def _infer_upper_bound_flags(sub: pd.DataFrame, rel_tol: float) -> pd.Series:
    if "at_upper_bound" in sub.columns:
        return sub["at_upper_bound"].fillna(False).astype(bool)
    vals = sub["abs_flux"].dropna()
    if vals.empty:
        return pd.Series(False, index=sub.index)
    ub = float(vals.max())
    # 避免把极小通量（如 [0, 1] 玩具夹具）误认为盒约束上界
    if not np.isfinite(ub) or ub < 100.0:
        return pd.Series(False, index=sub.index)
    return sub["abs_flux"].ge(ub * (1.0 - float(rel_tol))).fillna(False)


def _sample_loss_map(loss_df: pd.DataFrame | None) -> dict[str, float]:
    if loss_df is None or loss_df.empty:
        return {}
    sid_col = "sample_id" if "sample_id" in loss_df.columns else ("sample" if "sample" in loss_df.columns else None)
    val_col = next((c for c in ("min_loss", "mean_loss", "loss") if c in loss_df.columns), None)
    if sid_col is None or val_col is None:
        return {}
    grouped = loss_df.groupby(loss_df[sid_col].astype(str))[val_col].apply(
        lambda s: float(pd.to_numeric(s, errors="coerce").min())
    )
    return {k: v for k, v in grouped.items() if np.isfinite(v)}


def evaluate_mfa_quality(
    animals: pd.DataFrame | None,
    contrast_groups_by_id: dict[str, tuple[dict[str, str], str, str]] | None = None,
    loss_df: pd.DataFrame | None = None,
    config: dict[str, Any] | None = None,
) -> MFAQualityResult:
    """评估 MFA 样本级质量与比较级全局通量尺度漂移。

    Parameters
    ----------
    animals : pd.DataFrame
        ``sample_flux_matrix`` 产出的长表（每样本 × 反应一行），多示踪剂时含 ``tracer`` 列。
    contrast_groups_by_id : dict[str, tuple[dict[str, str], str, str]], optional
        ``{contrast_id: ({sample_id: level}, numerator, denominator)}``。
    loss_df : pd.DataFrame, optional
        ``quality_control`` 产出的拟合 loss 汇总表。
    config : dict, optional
        覆盖 :data:`DEFAULT_MFA_QC_CONFIG` 的配置。
    """
    cfg = resolve_mfa_qc_config(config)
    if not cfg.get("enabled", True) or animals is None or animals.empty:
        return MFAQualityResult(
            pd.DataFrame(columns=OUTPUT_QUALITY_COLUMNS),
            pd.DataFrame(),
            pd.DataFrame(columns=CONTRAST_SCALE_COLUMNS),
            set(), set(), set(), cfg,
        )

    df = _prepare_animals(animals)
    loss_map = _sample_loss_map(loss_df)
    min_n = int(cfg["min_samples_per_tracer"])
    min_flux = float(cfg["active_reaction_min_flux"])
    ub_rel = float(cfg["upper_bound_rel_tolerance"])

    quality_rows: list[dict[str, Any]] = []
    metric_rows: list[dict[str, Any]] = []
    contrast_rows: list[dict[str, Any]] = []

    for tracer, sub in df.groupby("tracer", sort=True):
        sub = sub.copy()
        sub["_at_ub"] = _infer_upper_bound_flags(sub, ub_rel)
        sample_ids = sorted(sub["sample_id"].astype(str).unique())

        rxn_med = sub.groupby("flux_id")["abs_flux"].median()
        active_rxns = rxn_med.index[rxn_med >= min_flux]
        n_active = int(len(active_rxns))

        per_s = sub.groupby(sub["sample_id"].astype(str)).agg(
            total_abs_flux=("abs_flux", "sum"),
            n_reactions=("flux_id", "nunique"),
            n_at_upper_bound=("_at_ub", "sum"),
            frac_at_upper_bound=("_at_ub", "mean"),
        ).reindex(sample_ids)

        if n_active > 0:
            sub_act = sub[sub["flux_id"].isin(active_rxns)].copy()
            sub_act["_log2_rel"] = np.log2(
                (sub_act["abs_flux"].to_numpy(dtype=float) + 1.0e-6)
                / (sub_act["flux_id"].map(rxn_med).to_numpy(dtype=float) + 1.0e-6)
            )
            per_s["med_log2_scale"] = sub_act.groupby(sub_act["sample_id"].astype(str))["_log2_rel"].median().reindex(sample_ids)
        else:
            per_s["med_log2_scale"] = 0.0

        med_tot = float(per_s["total_abs_flux"].median()) if per_s["total_abs_flux"].notna().any() else np.nan
        if np.isfinite(med_tot) and med_tot > 0:
            per_s["log2_total_flux_ratio"] = np.log2(per_s["total_abs_flux"].clip(lower=1.0e-9) / med_tot)
        else:
            per_s["log2_total_flux_ratio"] = 0.0

        per_s["min_loss"] = [loss_map.get(s, np.nan) for s in sample_ids]

        rz_med_scale = _robust_z(per_s["med_log2_scale"], floor=0.10)
        rz_tot_scale = _robust_z(per_s["log2_total_flux_ratio"], floor=0.10)
        rz_ub = _robust_z(per_s["n_at_upper_bound"], floor=1.0)
        rz_loss = _robust_z(per_s["min_loss"], floor=0.05) if per_s["min_loss"].notna().sum() >= 2 else pd.Series(np.nan, index=sample_ids)
        med_loss = float(per_s["min_loss"].median()) if per_s["min_loss"].notna().any() else np.nan

        can_stat_flag = len(sample_ids) >= min_n

        for s in sample_ids:
            ms = float(per_s.at[s, "med_log2_scale"])
            lt = float(per_s.at[s, "log2_total_flux_ratio"])
            n_ub = int(per_s.at[s, "n_at_upper_bound"])
            f_ub = float(per_s.at[s, "frac_at_upper_bound"])
            lv = float(per_s.at[s, "min_loss"]) if pd.notna(per_s.at[s, "min_loss"]) else np.nan
            z_ms = float(rz_med_scale.at[s])
            z_lt = float(rz_tot_scale.at[s])
            z_ub = float(rz_ub.at[s])
            z_lv = float(rz_loss.at[s]) if pd.notna(rz_loss.at[s]) else np.nan

            # 1. 全局通量尺度因子
            s_flag, s_reason, s_val, s_rz = "", "", ms, z_ms
            if can_stat_flag and n_active >= int(cfg["contrast_scale_min_reactions"]):
                if abs(ms) >= float(cfg["scale_flag_log2"]) and abs(z_ms) >= float(cfg["scale_flag_robust_z"]):
                    s_flag = "flag"
                    s_reason = f"单鼠全反应中位通量尺度异常（log2={ms:+.2f}, z={z_ms:+.1f}）"
                elif abs(lt) >= float(cfg["total_scale_flag_log2"]) and abs(z_lt) >= float(cfg["total_scale_flag_robust_z"]):
                    s_flag = "flag"
                    s_val, s_rz = lt, z_lt
                    s_reason = f"单鼠整鼠总通量尺度异常（log2={lt:+.2f}, z={z_lt:+.1f}）"
            quality_rows.append({
                "sample_id": s, "metric": "flux_scale_factor",
                "value": s_val, "robust_z": s_rz, "flag": s_flag, "reason": s_reason,
                "tracer": str(tracer),
                "detail": f"med_log2_scale={ms:+.3f} (z={z_ms:+.1f}); log2_total_ratio={lt:+.3f} (z={z_lt:+.1f})",
            })

            # 2. 上界贴边饱和度
            ub_flag, ub_reason = "", ""
            if (
                can_stat_flag
                and n_ub >= int(cfg["ub_flag_min_count"])
                and f_ub >= float(cfg["ub_flag_min_fraction"])
                and z_ub >= float(cfg["ub_flag_robust_z"])
            ):
                ub_flag = "flag"
                ub_reason = f"过多反应贴近盒约束上界（{n_ub} 条，占比 {f_ub:.2%}，z={z_ub:+.1f}）"
            quality_rows.append({
                "sample_id": s, "metric": "upper_bound_saturation",
                "value": float(n_ub), "robust_z": z_ub, "flag": ub_flag, "reason": ub_reason,
                "tracer": str(tracer),
                "detail": f"n_at_ub={n_ub}/{int(per_s.at[s, 'n_reactions'])} ({f_ub:.3%}, z={z_ub:+.1f})",
            })

            # 3. 拟合残差 loss
            if np.isfinite(lv) and np.isfinite(med_loss) and med_loss > 0:
                l_flag, l_reason = "", ""
                ratio = lv / med_loss
                if can_stat_flag and z_lv >= float(cfg["loss_exclude_robust_z"]) and ratio >= float(cfg["loss_exclude_min_ratio"]):
                    l_flag = "exclude"
                    l_reason = f"MFA 拟合残差极端偏高（loss={lv:.3g}, {ratio:.1f}x 中位数, z={z_lv:+.1f}）"
                elif can_stat_flag and z_lv >= float(cfg["loss_flag_robust_z"]) and ratio >= float(cfg["loss_flag_min_ratio"]):
                    l_flag = "flag"
                    l_reason = f"MFA 拟合残差偏高（loss={lv:.3g}, {ratio:.1f}x 中位数, z={z_lv:+.1f}）"
                quality_rows.append({
                    "sample_id": s, "metric": "fitting_loss",
                    "value": lv, "robust_z": z_lv, "flag": l_flag, "reason": l_reason,
                    "tracer": str(tracer),
                    "detail": f"min_loss={lv:.4g} (median={med_loss:.4g}, ratio={ratio:.2f}x, z={z_lv:+.1f})",
                })

            metric_rows.append({
                "sample_id": s, "tracer": str(tracer),
                "total_abs_flux": float(per_s.at[s, "total_abs_flux"]),
                "log2_total_flux_ratio": lt, "log2_total_flux_rz": z_lt,
                "med_log2_scale": ms, "med_log2_scale_rz": z_ms,
                "n_at_upper_bound": n_ub, "frac_at_upper_bound": f_ub, "ub_rz": z_ub,
                "min_loss": lv, "loss_rz": z_lv,
            })

        # ── 比较级整鼠尺度漂移与上界不对称评估 ────────────────────────────────
        if contrast_groups_by_id:
            min_rxns = int(cfg["contrast_scale_min_reactions"])
            tot_cut = float(cfg["contrast_scale_total_log2"])
            med_cut = float(cfg["contrast_scale_median_log2"])
            ub_cut = float(cfg["contrast_ub_asym_min_diff"])
            for cid, (groups, num_level, den_level) in contrast_groups_by_id.items():
                ns = [s for s, g in groups.items() if g == num_level and s in per_s.index]
                ds = [s for s, g in groups.items() if g == den_level and s in per_s.index]
                if len(ns) < 2 or len(ds) < 2:
                    continue
                d_tot = float(per_s.loc[ns, "log2_total_flux_ratio"].mean() - per_s.loc[ds, "log2_total_flux_ratio"].mean())
                d_med = float(per_s.loc[ns, "med_log2_scale"].mean() - per_s.loc[ds, "med_log2_scale"].mean())
                ub_n = float(per_s.loc[ns, "n_at_upper_bound"].mean())
                ub_d = float(per_s.loc[ds, "n_at_upper_bound"].mean())
                ub_diff = ub_n - ub_d

                same_sign = (d_tot > 0 and d_med > 0) or (d_tot < 0 and d_med < 0)
                scale_drift = bool(
                    n_active >= min_rxns
                    and same_sign
                    and abs(d_tot) >= tot_cut
                    and abs(d_med) >= med_cut
                )
                ub_asym = bool(abs(ub_diff) >= ub_cut and max(ub_n, ub_d) >= 1.0)
                marks, reasons = [], []
                if scale_drift:
                    marks.append("flux_scale_drift")
                    reasons.append(f"整鼠通量存在系统性尺度漂移（总通量差 {d_tot:+.2f} log2，中位尺度差 {d_med:+.2f} log2）")
                if ub_asym:
                    marks.append("upper_bound_asymmetry")
                    reasons.append(f"两组上界贴边不对称（平均贴顶反应数 {ub_n:.2f} vs {ub_d:.2f}）")
                contrast_rows.append({
                    "contrast_id": str(cid), "tracer": str(tracer),
                    "n_numerator": len(ns), "n_denominator": len(ds), "n_active_reactions": n_active,
                    "delta_log2_total_flux": d_tot, "delta_median_log2_scale": d_med,
                    "mean_ub_numerator": ub_n, "mean_ub_denominator": ub_d, "ub_count_diff": ub_diff,
                    "scale_drift": scale_drift, "ub_asymmetry": ub_asym,
                    "flag": ";".join(marks), "reason": "；".join(reasons),
                })

    quality_df = pd.DataFrame(quality_rows, columns=OUTPUT_QUALITY_COLUMNS)
    metrics_df = pd.DataFrame(metric_rows)
    contrast_df = pd.DataFrame(contrast_rows, columns=CONTRAST_SCALE_COLUMNS)
    excluded = set(quality_df.loc[quality_df["flag"] == "exclude", "sample_id"].astype(str))
    flagged = set(quality_df.loc[quality_df["flag"] == "flag", "sample_id"].astype(str)) - excluded
    drift = {
        (str(r.contrast_id), str(r.tracer))
        for r in contrast_df[contrast_df["scale_drift"]].itertuples(index=False)
    } if not contrast_df.empty else set()
    return MFAQualityResult(quality_df, metrics_df, contrast_df, flagged, excluded, drift, cfg)


def apply_scale_drift_flags(
    table: pd.DataFrame,
    contrast_scale: pd.DataFrame,
    *,
    tracer_col: str | None = None,
) -> pd.DataFrame:
    """把 ``flux_scale_drift`` 标记追加到差异通量表的 ``flags`` 列。"""
    if table is None or table.empty or contrast_scale is None or contrast_scale.empty:
        return table
    drifted = contrast_scale[contrast_scale["scale_drift"]]
    if drifted.empty:
        return table

    out = table.copy()
    existing = out["flags"].fillna("").astype(str) if "flags" in out.columns else pd.Series("", index=out.index)
    if tracer_col and tracer_col in out.columns:
        keys = set(zip(drifted["contrast_id"].astype(str), drifted["tracer"].astype(str)))
        hit = [
            (str(c), str(t)) in keys
            for c, t in zip(out["contrast_id"].astype(str), out[tracer_col].astype(str))
        ]
    else:
        # 单示踪剂或所有示踪剂均漂移的合并比较
        by_contrast = contrast_scale.groupby("contrast_id")["scale_drift"].all()
        drifted_cids = set(by_contrast[by_contrast].index.astype(str))
        hit = [str(c) in drifted_cids for c in out["contrast_id"].astype(str)]

    if not any(hit):
        return table
    out["flags"] = [
        ";".join(p for p in (f, "flux_scale_drift") if p) if m and "flux_scale_drift" not in f.split(";") else f
        for f, m in zip(existing, hit)
    ]
    return out
