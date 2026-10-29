#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MFA Multi-Batch / Multi-Protocol Integration Engine

对同一 protocol 下的多个 batch 执行随机效应 Meta 分析，
跨 protocol 执行方向一致性评分，输出整合后的差异通量表和报告。

Statistical method: DerSimonian-Laird random-effects estimator.
No external dependencies beyond numpy/scipy.

Usage:
    from mfa_integration import within_protocol_meta_analysis, cross_protocol_consensus

    results = within_protocol_meta_analysis(output_root, protocol_id, batch_ids, config)
"""

from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd
from scipy.stats import norm

from omicsweave.core.statistics.meta_analysis import random_effects_meta

# ─── 常量 ──────────────────────────────────────────────────────────────────

_INTEGRATION_DIR = "integration"
_HARMONIZED_SUBDIR = "tracer_harmonized"
_CROSS_TRACER_SUBDIR = "cross_tracer_scoring"

_DIFF_ANNOTATION_COLS = ["kegg_reaction", "kegg_pathway", "reactome_id", "ec_number", "description"]

# 拟合质量列：跨 protocol 表决时需要知道某条通量在各 protocol 下是否真的
# 被数据约束——未被约束的通量投票没有意义。
_QUALITY_COLS = ["identifiability", "median_rel_iqr", "frac_at_bound"]
# 视为"未被数据决定"的标签，与 df_analysis.UNCONSTRAINED_LABELS 保持同步：
# bound_pinned 指多数动物的拟合贴在盒约束边界，离散度为 0 但并非被约束。
_UNCONSTRAINED_LABELS = ("poorly_constrained", "bound_pinned")

# 由 discover_batch_differential_results() 填充：{batch_id: {tissue: 生成日期}}
# 供设计报告如实记录「保留但未纳入」的组织。
_LEGACY_TISSUES: dict[str, dict[str, str]] = {}


# ─── 差异通量结果发现 ─────────────────────────────────────────────────────

def normalize_tissue_aliases(
    df: pd.DataFrame,
    aliases: dict[str, str],
) -> pd.DataFrame:
    """把同一生理隔室的不同叫法归一到规范名。

    不同批次可能用 ``serum`` / ``plasma`` 指同一份循环血样本。名字不统一
    时，``flux_id`` 也随之不同（``serum--LDH_e`` vs ``plasma--LDH_e``），
    两个批次就永远无法在 meta 分析里对齐。别名在**载入后、分组前**统一
    重写 ``tissue`` 与 ``flux_id``。
    """
    if not aliases or df.empty or "tissue" not in df.columns:
        return df

    lowered = {str(k).lower(): str(v) for k, v in aliases.items()}
    original = df["tissue"].astype(str)
    canonical = original.str.lower().map(lowered).fillna(original)
    if canonical.equals(original):
        return df

    df = df.copy()
    df["tissue_original"] = original
    df["tissue"] = canonical
    if "reaction" in df.columns:
        df["flux_id"] = canonical + "--" + df["reaction"].astype(str)
    return df


def _current_run_tissues(base: Path) -> set[str] | None:
    """本次运行实际产出的组织集合（来自样本级通量矩阵）。

    流程是增量/幂等的，不会删除上一版数据留下的组织目录。若不加区分地
    glob，整合会把**过期组织**的旧 CSV 一并读进来。样本级通量矩阵每次
    运行整体重写，因此可作为「当前有效组织」的权威来源。
    """
    matrix_path = base / "all_tissues" / "flux_matrix" / "sample_net_flux.csv"
    if not matrix_path.exists():
        return None
    try:
        tissues = pd.read_csv(matrix_path, usecols=["tissue"])["tissue"]
    except (ValueError, KeyError):
        return None
    return {str(t) for t in tissues.unique()}


def discover_batch_differential_results(
    output_root: Path,
    protocol_id: str,
    batch_id: str,
) -> dict[tuple[str, str], Path]:
    """发现某个 batch 的差异通量 CSV 文件。

    key 必须是 **(comparison, tissue)**：同一个 comparison 在每个组织下
    都有一份 CSV，早期版本只用 comparison 作 key，导致同一批次里除最后
    一个组织外**全部被覆盖丢弃**，整合实际只看到了一个组织。

    Returns
    -------
    dict mapping (comp_name, tissue) to Path, e.g. ("HFD_vs_LFD", "liver").
    """
    base = output_root / "runs" / protocol_id / batch_id / "01_df_results"
    if not base.exists():
        return {}

    valid_tissues = _current_run_tissues(base)
    results: dict[tuple[str, str], Path] = {}
    stale: dict[str, str] = {}

    # 扫描 <tissue>/differential_flux/differential_flux_*.csv
    for csv_file in sorted(base.glob("*/differential_flux/differential_flux_*.csv")):
        tissue = csv_file.parent.parent.name
        if valid_tissues is not None and tissue not in valid_tissues:
            # 本次运行未产出该组织：文件是更早一版输入数据 + 更早一版代码
            # 的产物。文件保留在磁盘上，但不能与本次结果混进同一个统计
            # 池——两者的推断单位和阈值都不同。
            stale.setdefault(
                tissue,
                pd.Timestamp(csv_file.stat().st_mtime, unit="s").strftime("%Y-%m-%d"),
            )
            continue
        comp_name = csv_file.stem.replace("differential_flux_", "")
        results[(comp_name, tissue)] = csv_file

    # 兼容旧版扁平结构：differential_flux/*.csv（组织从表内 tissue 列取）
    if not results:
        for csv_file in sorted(base.glob("differential_flux/differential_flux_*.csv")):
            comp_name = csv_file.stem.replace("differential_flux_", "")
            results[(comp_name, "")] = csv_file

    if stale:
        detail = ", ".join(f"{t}(生成于 {d})" for t, d in sorted(stale.items()))
        print(f"    [discover] {batch_id}: {len(stale)} 个组织来自更早一版运行，"
              f"文件保留但不参与统计整合: {detail}")
        print(f"               如需纳入，请把这些组织放回 flux_raw_data.xlsx "
              f"并用当前版本重跑该批次。")

    _LEGACY_TISSUES[batch_id] = dict(stale)
    return results


def _load_differential_flux(csv_path: Path) -> pd.DataFrame:
    """加载差异通量 CSV，统一 key 列。"""
    df = pd.read_csv(csv_path)
    # 确保关键列存在
    for col in ("flux_id", "tissue", "reaction", "diff"):
        if col not in df.columns:
            return pd.DataFrame()
    return df


# ══════════════════════════════════════════════════════════════════════════
# DerSimonian-Laird 随机效应 Meta 分析
# ══════════════════════════════════════════════════════════════════════════

def _dersimonian_laird(
    estimates: np.ndarray,
    ses: np.ndarray,
) -> dict[str, float]:
    """DerSimonian-Laird 随机效应 Meta 分析。

    Parameters
    ----------
    estimates : array of effect sizes (diff)
    ses : array of standard errors

    Returns
    -------
    dict with: pooled_effect, pooled_se, ci_lower, ci_upper,
               Q_statistic, I_squared, tau_squared, k
    """
    k = len(estimates)
    if k == 0:
        return {}

    if k == 1:
        se = float(ses[0])
        est = float(estimates[0])
        return {
            "pooled_effect": est,
            "pooled_se": se,
            "ci_lower": est - 1.96 * se,
            "ci_upper": est + 1.96 * se,
            "Q_statistic": 0.0,
            "I_squared": np.nan,
            "tau_squared": 0.0,
            "k": 1,
        }

    # 固定效应权重
    w = 1.0 / (ses ** 2)
    w_sum = np.sum(w)

    # 固定效应合并效应
    theta_fe = np.sum(w * estimates) / w_sum

    # Cochran's Q
    Q = np.sum(w * (estimates - theta_fe) ** 2)

    # τ² (between-study variance)
    w2_sum = np.sum(w ** 2)
    denominator = w_sum - w2_sum / w_sum
    tau2 = max(0.0, (Q - k + 1) / denominator) if denominator > 0 else 0.0

    # 随机效应权重
    w_star = 1.0 / (ses ** 2 + tau2)
    w_star_sum = np.sum(w_star)

    # 随机效应合并效应
    theta_re = np.sum(w_star * estimates) / w_star_sum
    se_re = np.sqrt(1.0 / w_star_sum)

    # I²
    I2 = max(0.0, (Q - k + 1) / Q) * 100 if Q > 0 else 0.0

    return {
        "pooled_effect": float(theta_re),
        "pooled_se": float(se_re),
        "ci_lower": float(theta_re - 1.96 * se_re),
        "ci_upper": float(theta_re + 1.96 * se_re),
        "Q_statistic": float(Q),
        "I_squared": float(I2),
        "tau_squared": float(tau2),
        "k": k,
    }


def _compute_se_from_stats(row: pd.Series) -> float:
    """从差异通量统计列计算标准误。

    优先使用 std_A/n_A + std_B/n_B 公式。
    回退：如果无法计算则返回 NaN。
    """
    # 尝试从条件统计列计算
    std_cols_a = [c for c in row.index if c.startswith("std_") and c != "std"]
    mean_cols_a = [c for c in row.index if c.startswith("mean_") and c != "mean"]
    n_cols = [c for c in row.index if c.startswith("n_points_")]

    if len(std_cols_a) >= 2 and len(n_cols) >= 2:
        try:
            se_a = float(row[std_cols_a[0]]) / np.sqrt(max(1, float(row[n_cols[0]])))
            se_b = float(row[std_cols_a[1]]) / np.sqrt(max(1, float(row[n_cols[1]])))
            return np.sqrt(se_a ** 2 + se_b ** 2)
        except (ValueError, TypeError, ZeroDivisionError):
            pass

    # 回退：从 cohen_d 和 diff 估算
    diff = row.get("diff", np.nan)
    cohen_d = row.get("cohen_d", np.nan)
    if pd.notna(diff) and pd.notna(cohen_d) and abs(cohen_d) > 1e-10:
        return abs(float(diff)) / abs(float(cohen_d))

    return np.nan


def _compute_log2fc_se(row: pd.Series) -> tuple[float, str]:
    """返回与 ``log2FC`` 同尺度的 SE 及其来源。

    来源分三档：

    - ``reported`` —— 由组内方差经 delta method 直接算出，可用；
    - ``pvalue_normal_approximation`` —— 从 p 值反解（Cochrane Handbook
      §6.3 认可的插补），前提是该 p 值是**真实的检验结果**；
    - ``pvalue_at_sample_size_floor`` —— 反解所用的 p 值正好压在秩检验的
      样本量下界上（如 n=1 vs 3 时 p 恒为 0.5）。这种 p 值不含任何精度
      信息，由它反解出的 SE 只是 |效应量| 的常数倍，**不构成方差估计**。
    """
    explicit = pd.to_numeric(
        pd.Series([row.get("log2fc_se", np.nan)]), errors="coerce"
    ).iloc[0]
    if pd.notna(explicit) and np.isfinite(float(explicit)) and float(explicit) > 0:
        return abs(float(explicit)), "reported"

    effect = pd.to_numeric(
        pd.Series([row.get("log2FC", np.nan)]), errors="coerce"
    ).iloc[0]
    primary = str(row.get("primary_pvalue_column", "") or "")
    candidates = [primary, "pvalue_ttest", "pvalue_mannwhitney", "pvalue"]
    for column in dict.fromkeys(name for name in candidates if name):
        pvalue = pd.to_numeric(
            pd.Series([row.get(column, np.nan)]), errors="coerce"
        ).iloc[0]
        if (
            pd.notna(effect) and np.isfinite(float(effect)) and float(effect) != 0
            and pd.notna(pvalue) and np.isfinite(float(pvalue))
            and 0 < float(pvalue) < 1
        ):
            z_value = float(norm.isf(float(pvalue) / 2.0))
            if np.isfinite(z_value) and z_value > 0:
                floor = pd.to_numeric(
                    pd.Series([row.get("min_achievable_pvalue", np.nan)]),
                    errors="coerce",
                ).iloc[0]
                at_floor = (
                    pd.notna(floor)
                    and float(pvalue) <= float(floor) + 1e-12
                )
                source = ("pvalue_at_sample_size_floor" if at_floor
                          else "pvalue_normal_approximation")
                return abs(float(effect)) / z_value, source
    return np.nan, "unavailable"


def _bh_fdr(pvalues: pd.Series) -> pd.Series:
    values = pd.to_numeric(pvalues, errors="coerce")
    output = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.notna() & np.isfinite(values)
    if not valid.any():
        return output
    finite = values.loc[valid].clip(0.0, 1.0)
    order = finite.sort_values().index
    ranked = finite.loc[order].to_numpy(dtype=float)
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output.loc[order] = np.clip(adjusted, 0.0, 1.0)
    return output


# ══════════════════════════════════════════════════════════════════════════
# 批次覆盖度与整合路线判定
# ══════════════════════════════════════════════════════════════════════════

ROUTE_META = "random_effects_meta"
ROUTE_UNION = "union_single_batch"
ROUTE_DISJOINT = "union_disjoint_fluxes"
ROUTE_LEGACY = "excluded_legacy_run"


def read_noise_floor(output_root: Path, protocol_id: str) -> tuple[float, float]:
    """读取桥接 QC 得到的技术噪声地板 ``(p95, median)``（log2 尺度）。

    只有在 sample_id 确实标识同一只动物时（``interpretation ==
    technical_noise_floor``）才返回数值；否则返回 NaN——那种情况下的
    数字含生物学差异，拿它当阈值会把真实效应一并滤掉。
    """
    path = (output_root / _INTEGRATION_DIR / "bridging_qc" / protocol_id
            / "bridge_noise_floor.csv")
    if not path.exists():
        return float("nan"), float("nan")
    try:
        floor_df = pd.read_csv(path)
    except Exception:
        return float("nan"), float("nan")
    if "interpretation" in floor_df.columns:
        floor_df = floor_df[floor_df["interpretation"] == "technical_noise_floor"]
    if floor_df.empty:
        return float("nan"), float("nan")
    p95 = pd.to_numeric(floor_df.get("noise_floor_log2fc"), errors="coerce").max()
    med = pd.to_numeric(floor_df.get("median_abs_log2_ratio"), errors="coerce").max()
    return float(p95) if pd.notna(p95) else float("nan"), \
           float(med) if pd.notna(med) else float("nan")


def _se_quality(sources) -> str:
    """按 Meta 分析惯例给合并结果的方差来源分档。

    从 p 值反解 SE 是 Cochrane Handbook §6.3 认可的插补，**前提是那个
    p 值是真实的检验结果**。当 p 值压在秩检验的样本量下界上（n=1 vs 3
    时恒为 0.5），它不含精度信息，反解出的 SE 只是 |效应量| 的常数倍。
    这种"方差"不能支撑显著性主张，因此单独分出一档。
    """
    unique = set(str(s) for s in sources)
    if not unique or unique <= {"unavailable"}:
        return "unavailable"
    if "pvalue_at_sample_size_floor" in unique:
        return "imputed_from_floor_pvalue"
    if unique <= {"reported"}:
        return "all_reported"
    return "partially_imputed"


def build_integration_design_report(
    harmonized_dir: Path,
    protocol_id: str,
    batch_results: dict[str, dict[str, pd.DataFrame]],
    min_batches: int,
) -> pd.DataFrame:
    """写出跨批次整合的**安全性产物**：每个 (comparison x tissue) 的路线判定。

    这张表回答一个在合并之前必须先问清楚的问题：这些批次到底是彼此的
    重复，还是互补的组织板块？

    - 两个批次都覆盖同一 (comparison, tissue) 且反应名可对齐
      -> ``random_effects_meta``，随机效应合并有意义；
    - 只有一个批次覆盖
      -> ``union_single_batch``，只能并集拼接。它不是 pooled effect，
         不应与真正合并过的结果混在同一个 FDR 池里；
    - 多个批次覆盖同一组织但反应集不相交
      -> ``union_disjoint_fluxes``，通常是模型网络或命名不一致，
         需要人工核对而不是自动合并。

    与转录组的 ``joint_design_report.csv`` 作用相同：先读它，再决定
    整合结果能信到什么程度。
    """
    coverage: dict[tuple[str, str], dict[str, set]] = {}
    for batch_id, comp_dict in batch_results.items():
        for comp_name, df in comp_dict.items():
            if "tissue" not in df.columns:
                continue
            for tissue, group in df.groupby("tissue"):
                coverage.setdefault((comp_name, str(tissue)), {})[batch_id] = set(
                    group["flux_id"].astype(str)
                )

    rows = []
    for (comp_name, tissue), per_batch in sorted(coverage.items()):
        batch_ids = sorted(per_batch)
        flux_sets = [per_batch[b] for b in batch_ids]
        shared = set.intersection(*flux_sets) if len(flux_sets) > 1 else set(flux_sets[0])
        union = set.union(*flux_sets)

        if len(batch_ids) < min_batches:
            route = ROUTE_UNION
            note = "only_one_batch_covers_this_tissue_x_comparison"
        elif len(shared) == 0:
            route = ROUTE_DISJOINT
            note = "batches_share_the_tissue_but_no_common_flux_id"
        else:
            route = ROUTE_META
            note = "ok"

        rows.append({
            "protocol_id": protocol_id,
            "comparison": comp_name,
            "tissue": tissue,
            "n_batches": len(batch_ids),
            "batch_ids": ";".join(batch_ids),
            "n_flux_union": len(union),
            "n_flux_shared": len(shared),
            "flux_overlap_fraction": round(len(shared) / len(union), 4) if union else 0.0,
            "route": route,
            "note": note,
        })

    # 保留在磁盘上、但因来自更早一版运行而未纳入统计的组织，如实记录，
    # 不做静默丢弃。
    for batch_id, legacy in sorted(_LEGACY_TISSUES.items()):
        for tissue, generated in sorted(legacy.items()):
            rows.append({
                "protocol_id": protocol_id,
                "comparison": "(all)",
                "tissue": tissue,
                "n_batches": 0,
                "batch_ids": batch_id,
                "n_flux_union": 0,
                "n_flux_shared": 0,
                "flux_overlap_fraction": 0.0,
                "route": ROUTE_LEGACY,
                "note": f"retained_on_disk_from_earlier_run_{generated}",
            })

    design_df = pd.DataFrame(rows)
    if design_df.empty:
        return design_df

    harmonized_dir.mkdir(parents=True, exist_ok=True)
    design_df.to_csv(
        harmonized_dir / f"integration_design_report_{protocol_id}.csv", index=False
    )

    # ── 批次层面的组织重叠诊断 ────────────────────────────────────────
    batch_tissues = {
        batch_id: {
            str(t)
            for df in comp_dict.values()
            if "tissue" in df.columns
            for t in df["tissue"].unique()
        }
        for batch_id, comp_dict in batch_results.items()
    }
    if len(batch_tissues) > 1:
        shared_tissues = set.intersection(*batch_tissues.values())
        all_tissues = set.union(*batch_tissues.values())
        print(f"    [design] 组织覆盖: "
              + ", ".join(f"{b}={len(t)}" for b, t in sorted(batch_tissues.items()))
              + f" | 共有 {len(shared_tissues)}/{len(all_tissues)}")
        if not shared_tissues:
            print("    [design] !! 各批次组织集完全不重叠 —— 这些批次是**互补的组织板块**,")
            print("             不是彼此的技术重复。随机效应 Meta 分析在此不适用,")
            print("             整合退化为并集拼接(union),结果不构成跨批次验证。")
        elif len(shared_tissues) < len(all_tissues) / 2:
            print(f"    [design] !  组织重叠有限（共有 {sorted(shared_tissues)}）,"
                  f" 大部分 (comparison x tissue) 只能走 union 路线。")

    route_counts = design_df["route"].value_counts().to_dict()
    print(f"    [design] 路线判定: {route_counts}")
    print(f"    [design] → integration_design_report_{protocol_id}.csv")
    return design_df


# ─── 同 Protocol 内 Meta 分析 ──────────────────────────────────────────────

def within_protocol_meta_analysis(
    output_root: Path,
    protocol_id: str,
    batch_ids: list[str],
    integration_config: dict[str, Any],
) -> dict[str, pd.DataFrame]:
    """对同一 protocol 的多个 batch 执行随机效应 Meta 分析。

    Parameters
    ----------
    output_root : Path
        mfa_results 根目录
    protocol_id : str
    batch_ids : list[str]
        有数据的 batch ID 列表
    integration_config : dict
        experiment.yaml 中 integration 段

    Returns
    -------
    dict mapping comp_name to DataFrame with meta-analysis results.
    """
    noise_floor_p95, noise_floor_median = read_noise_floor(output_root, protocol_id)
    if np.isfinite(noise_floor_p95):
        print(f"  [integration] 技术噪声地板: |log2FC| P95={noise_floor_p95:.3f}"
              f" (≈{2 ** noise_floor_p95:.1f}x), "
              f"中位数={noise_floor_median:.3f}"
              f" (≈{2 ** noise_floor_median:.1f}x)")

    tissue_aliases = integration_config.get("tissue_aliases", {}) or {}
    if tissue_aliases:
        print(f"  [integration] 组织别名: {tissue_aliases}")
    wp_config = integration_config.get("within_protocol", {}) or {}
    min_batches = wp_config.get("min_batches", 2)
    heterogeneity_threshold = wp_config.get("heterogeneity_threshold", 0.50)
    pvalue_threshold = float(wp_config.get("pvalue_threshold", 0.05))
    fdr_threshold = float(wp_config.get("fdr_threshold", 0.10))
    log2fc_threshold = float(wp_config.get("log2fc_threshold", 0.0))

    # 发现每个 batch 的差异通量文件
    # 旧版组织登记表按 protocol 重置，避免上一个 protocol 的记录串进来。
    _LEGACY_TISSUES.clear()
    batch_results: dict[str, dict[str, pd.DataFrame]] = {}
    all_comps: set[str] = set()

    for bid in batch_ids:
        comp_paths = discover_batch_differential_results(output_root, protocol_id, bid)
        batch_dfs: dict[str, pd.DataFrame] = {}
        frames_by_comp: dict[str, list[pd.DataFrame]] = {}
        for (comp_name, _tissue), csv_path in comp_paths.items():
            df = _load_differential_flux(csv_path)
            if df.empty:
                continue
            df = normalize_tissue_aliases(df, tissue_aliases)
            df["_batch_id"] = bid
            uncertainty = df.apply(_compute_log2fc_se, axis=1)
            df["_se"] = [value[0] for value in uncertainty]
            df["_se_source"] = [value[1] for value in uncertainty]
            df["_effect"] = pd.to_numeric(df.get("log2FC"), errors="coerce")
            # 同一 comparison 的多个组织必须纵向拼接，而不是互相覆盖。
            frames_by_comp.setdefault(comp_name, []).append(df)
            all_comps.add(comp_name)
        for comp_name, frames in frames_by_comp.items():
            batch_dfs[comp_name] = pd.concat(frames, ignore_index=True)
        if batch_dfs:
            batch_results[bid] = batch_dfs

    if not all_comps:
        print(f"  [integration] 无差异通量数据可整合 (protocol={protocol_id})")
        return {}

    n_data_batches = len(batch_results)
    print(f"  [integration] protocol={protocol_id}: "
          f"{n_data_batches} batch(es) with data, {len(all_comps)} comparison(s)")

    # 输出目录
    harmonized_dir = output_root / _INTEGRATION_DIR / _HARMONIZED_SUBDIR / protocol_id
    harmonized_dir.mkdir(parents=True, exist_ok=True)

    # 合并之前先判定路线：这些批次是重复，还是互补的组织板块？
    build_integration_design_report(
        harmonized_dir, protocol_id, batch_results, min_batches
    )

    meta_results = {}

    for comp_name in sorted(all_comps):
        # 收集该 comparison 在各 batch 中的数据
        comp_frames = []
        for bid, comp_dict in batch_results.items():
            if comp_name in comp_dict:
                comp_frames.append(comp_dict[comp_name])

        if not comp_frames:
            continue

        # 合并所有 batch 的数据
        all_data = pd.concat(comp_frames, ignore_index=True)

        # 按 (flux_id, tissue) 分组进行 Meta 分析
        group_key = ["flux_id", "tissue"]
        meta_rows = []

        for (fid, tis), grp in all_data.groupby(group_key):
            # Meta 的效应量和 SE 必须都在 log2FC 尺度。
            valid = grp.dropna(subset=["_effect", "_se"])
            valid = valid[valid["_se"] > 0]

            if len(valid) == 0:
                continue

            # 取第一行的注释列和 reaction 名
            first = valid.iloc[0]
            row_base = {
                "flux_id": fid,
                "tissue": tis,
                "reaction": first.get("reaction", fid.split("--")[-1] if "--" in fid else fid),
                "n_batches": len(valid),
                "n_batches_observed": int(grp["_batch_id"].nunique()),
                "effect_scale": "log2FC",
            }

            # 添加注释列与拟合质量列
            for col in _DIFF_ANNOTATION_COLS + _QUALITY_COLS:
                row_base[col] = first.get(col, np.nan)

            if len(valid) == 1:
                # 单 batch 直通
                row = row_base.copy()
                row["diff_pooled"] = float(valid.iloc[0]["diff"])
                row["se_pooled"] = float(valid.iloc[0]["_se"])
                ci_se = float(valid.iloc[0]["_se"])
                pooled_effect = float(valid.iloc[0]["_effect"])
                row["ci_lower"] = pooled_effect - 1.96 * ci_se
                row["ci_upper"] = pooled_effect + 1.96 * ci_se
                row["log2FC_pooled"] = pooled_effect
                row["cohen_d_pooled"] = float(valid.iloc[0].get("cohen_d", np.nan))
                row["hedges_g_pooled"] = float(valid.iloc[0].get("hedges_g", np.nan))
                row["Q_statistic"] = 0.0
                row["I_squared"] = np.nan
                row["tau_squared"] = 0.0
                row["pooled_pvalue"] = float(
                    2.0 * norm.sf(abs(pooled_effect / ci_se))
                )
                row["is_significant_pooled"] = False
                row["meta_note"] = "single_batch_passthrough"
                row["integration_route"] = ROUTE_UNION
                # 该 p 值来自单个批次的原始检验，不是 pooled effect。
                # 之所以看起来"像"合并结果，只是因为 SE 是由原 p 值反解的，
                # 再正向算回去必然复现原值——这不构成任何跨批次证据。
                row["pvalue_origin"] = "single_batch_original"
                row["batch_ids"] = valid.iloc[0]["_batch_id"]
                row["se_source"] = str(valid.iloc[0]["_se_source"])
                row["se_quality"] = _se_quality([valid.iloc[0]["_se_source"]])
                row["direction_consistent_across_batches"] = False
                meta_rows.append(row)
            elif len(valid) >= min_batches:
                # 随机效应 Meta 分析
                estimates = valid["_effect"].values.astype(float)
                ses = valid["_se"].values.astype(float)

                # 过滤 NaN
                mask = np.isfinite(estimates) & np.isfinite(ses) & (ses > 0)
                if mask.sum() < min_batches:
                    continue

                labels = valid.loc[mask, "_batch_id"].astype(str).tolist()
                ma = random_effects_meta(estimates[mask], ses[mask], labels)

                row = row_base.copy()
                row["n_batches"] = ma.k
                row["diff_pooled"] = float(
                    pd.to_numeric(valid.get("diff"), errors="coerce").mean()
                )
                row["se_pooled"] = ma.pooled_se
                row["ci_lower"] = ma.ci_lower
                row["ci_upper"] = ma.ci_upper
                row["log2FC_pooled"] = ma.pooled_effect
                row["pooled_pvalue"] = ma.pvalue

                cohen_d_vals = valid["cohen_d"].dropna().values
                row["cohen_d_pooled"] = float(np.mean(cohen_d_vals)) if len(cohen_d_vals) > 0 else np.nan
                hedges_g_vals = pd.to_numeric(
                    valid.get("hedges_g", pd.Series(dtype=float)), errors="coerce"
                ).dropna().values
                row["hedges_g_pooled"] = (
                    float(np.mean(hedges_g_vals)) if len(hedges_g_vals) > 0 else np.nan
                )

                row["Q_statistic"] = ma.q_statistic
                row["I_squared"] = ma.i_squared
                row["tau_squared"] = ma.tau_squared
                row["is_significant_pooled"] = False

                # 异质性标记
                if pd.notna(ma.i_squared) and ma.i_squared / 100 > heterogeneity_threshold:
                    row["meta_note"] = "high_heterogeneity"
                else:
                    row["meta_note"] = "random_effects_meta"

                row["integration_route"] = ROUTE_META
                row["pvalue_origin"] = "random_effects_pooled"
                row["batch_ids"] = ";".join(ma.source_ids)
                sources = valid.loc[mask, "_se_source"].astype(str).tolist()
                row["se_source"] = ";".join(sorted(set(sources)))
                row["se_quality"] = _se_quality(sources)
                row["direction_consistent_across_batches"] = bool(ma.direction_consistent)
                meta_rows.append(row)

        if not meta_rows:
            continue

        meta_df = pd.DataFrame(meta_rows)
        if "integration_route" not in meta_df.columns:
            meta_df["integration_route"] = ROUTE_UNION

        # FDR **按路线分池**：真正合并过的 pooled effect 与仅单批次直通的
        # 结果是两类不同强度的证据，混在同一个 BH 池里会让后者搭上前者的
        # 便车（或反之被稀释）。
        meta_df["pooled_fdr"] = np.nan
        for route, idx in meta_df.groupby("integration_route").groups.items():
            meta_df.loc[idx, "pooled_fdr"] = _bh_fdr(
                meta_df.loc[idx, "pooled_pvalue"]
            )
        meta_df["fdr_pool"] = meta_df["integration_route"]
        meta_df["fdr_pool_size"] = meta_df.groupby("integration_route")[
            "pooled_pvalue"
        ].transform("size")

        meta_df["is_significant_pooled"] = (
            meta_df["pooled_pvalue"].le(pvalue_threshold)
            & meta_df["pooled_fdr"].le(fdr_threshold)
            & meta_df["log2FC_pooled"].abs().ge(log2fc_threshold)
        ).fillna(False)

        if "se_quality" not in meta_df.columns:
            meta_df["se_quality"] = "unavailable"
        if "direction_consistent_across_batches" not in meta_df.columns:
            meta_df["direction_consistent_across_batches"] = False

        # 跨批次证据分两档：
        #   is_cross_batch_supported —— 真正合并过、pooled 显著，且方差估计
        #     不是从下界 p 值反解来的（见 _se_quality）。这是最强的一档。
        #   cross_batch_direction_agree —— 只要求两批次效应同向。当样本量
        #     不足以给出可信方差时，方向一致仍是有意义的弱证据。
        # 技术噪声地板：同一只动物在两个批次里各拟合一次的差异，不含
        # 生物学信号。合并效应量若没有超过它，就无法与技术重现性区分。
        meta_df["noise_floor_log2fc_p95"] = noise_floor_p95
        meta_df["noise_floor_log2fc_median"] = noise_floor_median
        abs_pooled = meta_df["log2FC_pooled"].abs()
        meta_df["above_noise_floor_median"] = (
            abs_pooled.ge(noise_floor_median) if np.isfinite(noise_floor_median)
            else pd.Series(True, index=meta_df.index)
        ).fillna(False)
        meta_df["above_noise_floor_p95"] = (
            abs_pooled.ge(noise_floor_p95) if np.isfinite(noise_floor_p95)
            else pd.Series(True, index=meta_df.index)
        ).fillna(False)

        meta_df["is_cross_batch_supported"] = (
            meta_df["is_significant_pooled"]
            & meta_df["integration_route"].eq(ROUTE_META)
            & ~meta_df["se_quality"].eq("imputed_from_floor_pvalue")
            & ~meta_df["se_quality"].eq("unavailable")
            & meta_df["above_noise_floor_p95"]
        )
        meta_df["cross_batch_direction_agree"] = (
            meta_df["integration_route"].eq(ROUTE_META)
            & meta_df["direction_consistent_across_batches"]
        )

        # 排序
        sort_col = "diff_pooled" if "diff_pooled" in meta_df.columns else None
        if sort_col:
            meta_df = meta_df.reindex(
                meta_df[sort_col].abs().sort_values(ascending=False).index
            ).reset_index(drop=True)

        # 保存
        out_path = harmonized_dir / f"differential_flux_meta_{comp_name}.csv"
        meta_df.to_csv(out_path, index=False)
        n_meta_route = int(meta_df["integration_route"].eq(ROUTE_META).sum())
        n_union_route = int(meta_df["integration_route"].eq(ROUTE_UNION).sum())
        print(f"    → {out_path.name}: {len(meta_df)} fluxes "
              f"({n_meta_route} meta-analysed, {n_union_route} union-only), "
              f"{int(meta_df['is_significant_pooled'].sum())} significant, "
              f"{int(meta_df['is_cross_batch_supported'].sum())} cross-batch supported, "
              f"{int(meta_df['cross_batch_direction_agree'].sum())} direction-agree")
        if np.isfinite(noise_floor_p95) and n_meta_route:
            meta_only = meta_df[meta_df["integration_route"].eq(ROUTE_META)]
            print(f"      [NOISE] {int(meta_only['above_noise_floor_p95'].sum())}"
                  f"/{n_meta_route} 条合并效应超过技术噪声地板 P95")
        if n_meta_route and not int(meta_df["is_cross_batch_supported"].sum()):
            floor_rows = int(
                meta_df["integration_route"].eq(ROUTE_META).__and__(
                    meta_df["se_quality"].eq("imputed_from_floor_pvalue")
                ).sum()
            )
            if floor_rows:
                print(f"      [SE] {floor_rows}/{n_meta_route} 条合并结果的方差来自"
                      f"下界 p 值反解，按惯例不支撑显著性主张")

        meta_results[comp_name] = meta_df

    # 生成 meta 分析汇总表
    _write_meta_summary(harmonized_dir, protocol_id, meta_results, n_data_batches)

    return meta_results


def _write_meta_summary(
    harmonized_dir: Path,
    protocol_id: str,
    meta_results: dict[str, pd.DataFrame],
    n_batches: int,
) -> None:
    """生成同 protocol 内 meta 分析的汇总统计表。"""
    rows = []
    for comp_name, df in meta_results.items():
        n_total = len(df)
        n_sig = int(df["is_significant_pooled"].sum())
        n_high_het = int((df["meta_note"] == "high_heterogeneity").sum())
        n_passthrough = int((df["meta_note"] == "single_batch_passthrough").sum())
        n_meta = n_total - n_passthrough

        avg_I2 = df["I_squared"].dropna().mean() if n_meta > 0 else np.nan
        avg_tau2 = df["tau_squared"].mean() if "tau_squared" in df.columns else np.nan

        n_cross_batch = (
            int(df["is_cross_batch_supported"].sum())
            if "is_cross_batch_supported" in df.columns else 0
        )
        rows.append({
            "protocol_id": protocol_id,
            "comparison": comp_name,
            "n_data_batches": n_batches,
            "n_fluxes": n_total,
            "n_significant": n_sig,
            "n_cross_batch_supported": n_cross_batch,
            "n_meta_analyzed": n_meta,
            "n_passthrough": n_passthrough,
            "n_high_heterogeneity": n_high_het,
            "avg_I_squared": round(avg_I2, 2) if pd.notna(avg_I2) else np.nan,
            "avg_tau_squared": round(avg_tau2, 4) if pd.notna(avg_tau2) else np.nan,
        })

    if rows:
        summary_df = pd.DataFrame(rows)
        summary_df.to_csv(harmonized_dir / f"meta_analysis_summary_{protocol_id}.csv", index=False)


# ══════════════════════════════════════════════════════════════════════════
# 跨 Protocol 方向一致性分析
# ══════════════════════════════════════════════════════════════════════════

def cross_protocol_consensus(
    output_root: Path,
    protocol_ids: list[str],
    integration_config: dict[str, Any],
) -> Optional[pd.DataFrame]:
    """跨协议方向一致性分析。

    Parameters
    ----------
    output_root : Path
    protocol_ids : list[str]
        有数据且 enabled 的协议列表。
    integration_config : dict
        experiment.yaml 中 integration 段。

    Returns
    -------
    DataFrame with consensus results, or None if insufficient data.
    """
    if len(protocol_ids) < 2:
        print("  [integration] 跨 protocol 整合需要 ≥2 个有数据的 protocol，跳过")
        return None

    cp_config = integration_config.get("cross_protocol", {}) or {}
    min_protocols = cp_config.get("min_protocols", 2)
    require_same_direction = cp_config.get("require_same_direction", True)
    # 仅凭方向一致就"通过"，会让两个 protocol 里都只是噪声的通量也过关——
    # 符号本身只有两种取值，随机也有一半概率一致。因此默认还要求：
    # 至少 min_protocols_significant 个 protocol 判定该通量显著，且该通量
    # 不属于 poorly_constrained / bound_pinned（exclude_poorly_constrained
    # 沿用旧键名，覆盖 _UNCONSTRAINED_LABELS 全部标签）。
    min_protocols_significant = int(
        cp_config.get("min_protocols_significant", min_protocols)
    )
    exclude_poorly_constrained = bool(
        cp_config.get("exclude_poorly_constrained", True)
    )

    # 加载各 protocol 的 harmonized 结果
    protocol_data: dict[str, dict[str, pd.DataFrame]] = {}
    for pid in protocol_ids:
        harmonized_dir = output_root / _INTEGRATION_DIR / _HARMONIZED_SUBDIR / pid
        if not harmonized_dir.exists():
            continue
        comp_dfs = {}
        for csv_file in harmonized_dir.glob("differential_flux_meta_*.csv"):
            comp_name = csv_file.stem.replace("differential_flux_meta_", "")
            comp_dfs[comp_name] = pd.read_csv(csv_file)
        if comp_dfs:
            protocol_data[pid] = comp_dfs

    if len(protocol_data) < min_protocols:
        print(f"  [integration] 仅有 {len(protocol_data)} 个 protocol 有 harmonized 数据，"
              f"需要 ≥{min_protocols}，跳过跨 protocol 分析")
        return None

    # 收集所有 comparison
    all_comps: set[str] = set()
    for comp_dict in protocol_data.values():
        all_comps.update(comp_dict.keys())

    consensus_rows = []

    for comp_name in sorted(all_comps):
        # 收集各 protocol 的数据
        frames = []
        for pid, comp_dict in protocol_data.items():
            if comp_name in comp_dict:
                df = comp_dict[comp_name].copy()
                df["_protocol_id"] = pid
                frames.append(df)

        if len(frames) < min_protocols:
            continue

        merged = pd.concat(frames, ignore_index=True)

        # 按 (flux_id, tissue) 分组
        for (fid, tis), grp in merged.groupby(["flux_id", "tissue"]):
            if len(grp) < min_protocols:
                continue

            # 各 protocol 的方向
            log2fc_vals = grp["log2FC_pooled"].dropna().values
            if len(log2fc_vals) < min_protocols:
                continue

            signs = np.sign(log2fc_vals)
            n_up = int(np.sum(signs > 0))
            n_down = int(np.sum(signs < 0))
            n_agree = max(n_up, n_down)
            n_protocols = len(signs)

            consensus_score = n_agree / n_protocols

            # 各 protocol 是否给出了实质证据
            if "is_significant_pooled" in grp.columns:
                n_significant = int(
                    grp["is_significant_pooled"].fillna(False).astype(bool).sum()
                )
            else:
                n_significant = 0
            if "identifiability" in grp.columns:
                identifiability_values = grp["identifiability"].astype(str)
                unconstrained = identifiability_values.isin(_UNCONSTRAINED_LABELS)
                n_well_constrained = int((~unconstrained).sum())
                any_poorly_constrained = bool(unconstrained.any())
            else:
                n_well_constrained = n_protocols
                any_poorly_constrained = False

            # 多数方向
            if n_up > n_down:
                direction = "up"
            elif n_down > n_up:
                direction = "down"
            else:
                direction = "mixed"

            # 通过条件：方向一致 + 足够多的 protocol 给出显著证据 +
            # （默认）该通量在任一 protocol 中都不是未被约束的
            consensus_pass = True
            if require_same_direction and consensus_score < 1.0:
                consensus_pass = False
            if n_significant < min_protocols_significant:
                consensus_pass = False
            if exclude_poorly_constrained and any_poorly_constrained:
                consensus_pass = False

            first = grp.iloc[0]
            row = {
                "flux_id": fid,
                "tissue": tis,
                "reaction": first.get("reaction", ""),
                "comparison": comp_name,
                "n_protocols": n_protocols,
                "n_agree": n_agree,
                "n_up": n_up,
                "n_down": n_down,
                "directional_consensus_score": round(consensus_score, 3),
                "direction": direction,
                "n_protocols_significant": n_significant,
                "n_protocols_well_constrained": n_well_constrained,
                "min_protocols_significant": min_protocols_significant,
                "excluded_poorly_constrained": bool(
                    exclude_poorly_constrained and any_poorly_constrained
                ),
                "consensus_pass": consensus_pass,
                "protocol_details": ";".join(
                    f"{r['_protocol_id']}:{'up' if r['log2FC_pooled'] > 0 else 'down'}"
                    for _, r in grp.iterrows() if pd.notna(r.get("log2FC_pooled"))
                ),
            }

            # 注释列与拟合质量列
            for col in _DIFF_ANNOTATION_COLS + _QUALITY_COLS:
                row[col] = first.get(col, np.nan)

            consensus_rows.append(row)

    if not consensus_rows:
        return None

    consensus_df = pd.DataFrame(consensus_rows)

    # 输出
    scoring_dir = output_root / _INTEGRATION_DIR / _CROSS_TRACER_SUBDIR
    scoring_dir.mkdir(parents=True, exist_ok=True)

    for comp_name in sorted(consensus_df["comparison"].unique()):
        comp_df = consensus_df[consensus_df["comparison"] == comp_name]
        out_path = scoring_dir / f"directional_consensus_{comp_name}.csv"
        comp_df.to_csv(out_path, index=False)
        n_pass = int(comp_df["consensus_pass"].sum())
        n_direction = int(comp_df["directional_consensus_score"].ge(1.0).sum())
        print(f"    → {out_path.name}: {len(comp_df)} fluxes, "
              f"{n_direction} direction-consistent, {n_pass} pass consensus "
              f"(需 >={min_protocols_significant} 个 protocol 显著"
              f"{'，且排除未约束通量' if exclude_poorly_constrained else ''})")

    # Protocol 对比汇总
    summary_rows = []
    for comp_name in sorted(consensus_df["comparison"].unique()):
        comp_df = consensus_df[consensus_df["comparison"] == comp_name]
        summary_rows.append({
            "comparison": comp_name,
            "n_fluxes": len(comp_df),
            "n_consensus_pass": int(comp_df["consensus_pass"].sum()),
            "n_direction_up": int((comp_df["direction"] == "up").sum()),
            "n_direction_down": int((comp_df["direction"] == "down").sum()),
            "n_direction_mixed": int((comp_df["direction"] == "mixed").sum()),
            "avg_consensus_score": round(comp_df["directional_consensus_score"].mean(), 3),
        })
    if summary_rows:
        pd.DataFrame(summary_rows).to_csv(
            scoring_dir / "protocol_comparison_summary.csv", index=False)

    return consensus_df


# ══════════════════════════════════════════════════════════════════════════
# 整合分析 Markdown 报告
# ══════════════════════════════════════════════════════════════════════════


def _design_report_section(harmonized_dir: Path, protocol_id: str) -> list[str]:
    """整合路线判定小节：先说清批次之间到底是什么关系。"""
    path = harmonized_dir / f"integration_design_report_{protocol_id}.csv"
    if not path.exists():
        return []

    design = pd.read_csv(path)
    if design.empty:
        return []

    lines = ["### 整合路线判定（先读这一节）", ""]
    counts = design["route"].value_counts()
    n_meta = int(counts.get(ROUTE_META, 0))
    n_union = int(counts.get(ROUTE_UNION, 0))
    n_disjoint = int(counts.get(ROUTE_DISJOINT, 0))
    total = len(design)

    lines.append(f"共 {total} 个 (comparison x tissue) 组合：")
    lines.append("")
    lines.append(f"- `{ROUTE_META}`（可做随机效应合并）: **{n_meta}**")
    lines.append(f"- `{ROUTE_UNION}`（仅单批次覆盖，只能并集拼接）: **{n_union}**")
    lines.append(f"- `{ROUTE_DISJOINT}`（同组织但反应集不相交，需人工核对）: **{n_disjoint}**")
    n_legacy = int(counts.get(ROUTE_LEGACY, 0))
    if n_legacy:
        lines.append(f"- `{ROUTE_LEGACY}`（更早一版运行的产物，文件保留但不参与统计）: **{n_legacy}**")
    lines.append("")
    if n_legacy:
        lines.append("> 标为 `excluded_legacy_run` 的组织目录仍完整保留在磁盘上，")
        lines.append("> 只是它们由更早一版输入数据 + 更早一版代码生成，推断单位与")
        lines.append("> 显著性阈值都与本次不同，混入同一统计池会破坏可比性。")
        lines.append("> 如需纳入，请把这些组织放回 `flux_raw_data.xlsx` 后重跑该批次。")
        lines.append("")

    if n_meta == 0:
        lines.append("> **本 protocol 下没有任何组合能真正做跨批次合并。**")
        lines.append("> 各批次覆盖的是互补的组织板块，不是彼此的技术重复；")
        lines.append("> 下方所有 pooled 结果实为单批次直通（union），")
        lines.append("> 不构成跨批次验证，`is_cross_batch_supported` 恒为 False。")
        lines.append("")

    tissue_view = (
        design.groupby(["tissue", "route"], as_index=False)
        .agg(n_comparisons=("comparison", "nunique"),
             n_batches=("n_batches", "max"),
             n_flux_shared=("n_flux_shared", "max"),
             n_flux_union=("n_flux_union", "max"))
        .sort_values(["route", "tissue"])
    )
    lines.append(tissue_view.to_markdown(index=False))
    lines.append("")
    return lines


def _bridging_section(output_root: Path, protocol_id: str) -> list[str]:
    """桥接样本噪声地板小节。"""
    path = (output_root / _INTEGRATION_DIR / "bridging_qc" / protocol_id
            / "bridge_noise_floor.csv")
    if not path.exists():
        return []

    floor_df = pd.read_csv(path)
    if floor_df.empty:
        return []

    lines = ["### 桥接样本一致性与技术噪声地板", ""]
    usable = floor_df[floor_df["n_pairs"] > 0]
    if usable.empty:
        lines.append("各批次之间没有共同的 (样本, 通量) 组合，无法估计技术噪声地板。")
        lines.append("")
        lines.append("> 建议：后续批次固定纳入 2-3 个桥接样本（同一动物分装或混合 QC 池），")
        lines.append("> 否则跨批次差异中的技术成分永远无法与生物学成分分离。")
        lines.append("")
        return lines

    display_cols = [c for c in ["batch_a", "batch_b", "n_bridge_samples",
                                "n_bridge_fluxes", "n_pairs", "n_usable_pairs",
                                "sign_agreement_rate", "median_abs_log2_ratio",
                                "noise_floor_log2fc"]
                    if c in usable.columns]
    lines.append(usable[display_cols].to_markdown(index=False))
    lines.append("")

    floors = usable["noise_floor_log2fc"].dropna()
    identity = str(usable.get("identity_basis", pd.Series([""])).iloc[0])
    if len(floors):
        worst = float(floors.max())
        if identity == "group_index_label":
            lines.append("> **这不是技术噪声地板。** 上游 `complete_result_name` 是")
            lines.append("> 「组名+序号」标签（如 `glucose__male__hfd_sham_8w__1`），")
            lines.append("> 不含动物 ID，两个批次里的同名样本未必是同一只动物。")
            lines.append("")
            lines.append(f"按同名配对得到的 |log2 比值| 在 95% 分位处为 **{worst:.3f}**"
                         f"（约 {2 ** worst:.2f} 倍），它是")
            lines.append("「技术重现性 + 生物学个体差异 + 可能的配对错误」三者之和的**上界**。")
            lines.append("")
            lines.append("要把它变成真正可用的技术噪声地板，需要二者之一：")
            lines.append("")
            lines.append("1. 在 `sheet_information.complete_result_name` 中写入动物 ID")
            lines.append("   （与 `groups.yaml` 的 `sample_groups` 键一致），使跨批次配对可验证；")
            lines.append("2. 后续批次固定纳入 2-3 个真正的桥接样本（同一动物分装或混合 QC 池）。")
            lines.append("")
        else:
            medians = pd.to_numeric(
                usable.get("median_abs_log2_ratio"), errors="coerce"
            ).dropna()
            typical = float(medians.max()) if len(medians) else float("nan")
            lines.append("**解读**：同一只动物的同一组织在两个批次中被独立拟合，")
            lines.append("两次估计之差不含任何生物学信号——它就是技术重现性本身。")
            lines.append("")
            if pd.notna(typical):
                lines.append(f"- 典型差异（中位数）：|log2FC| **{typical:.3f}**"
                             f"，约 **{2 ** typical:.2f} 倍**")
            lines.append(f"- 重复性限（P95，ISO 5725 意义上的临界差）："
                         f"|log2FC| **{worst:.3f}**，约 **{2 ** worst:.2f} 倍**")
            lines.append("")
            lines.append(f"因此 **|log2FC| < {worst:.3f} 的跨批次合并效应无法与技术")
            lines.append("重现性区分**，不应作为结论依据。整合结果中的")
            lines.append("`above_noise_floor_p95` / `above_noise_floor_median` 两列记录了")
            lines.append("每条通量是否越过该门槛，`is_cross_batch_supported` 已把 P95 门槛")
            lines.append("纳入判定。")
            lines.append("")
            if pd.notna(typical) and typical >= 0.58:
                lines.append(f"> **注意**：中位技术差异（{2 ** typical:.1f} 倍）已经高于")
                lines.append("> 本实验配置的生物学效应阈值 |log2FC| ≥ 0.58（1.5 倍）。")
                lines.append("> 这意味着在当前拟合重现性下，跨批次比较基本无法分辨")
                lines.append("> 该量级的生物学效应——优先排查拟合收敛性与测量重复性，")
                lines.append("> 而不是继续在结果层面做统计。")
                lines.append("")
    return lines

def generate_integration_report(
    output_root: Path,
    protocol_ids: list[str],
    integration_config: dict[str, Any],
    experiment_id: str,
) -> Path:
    """生成整合分析 Markdown 报告。

    输出: <output_root>/integration/integration_report.md
    """
    lines = [
        f"# MFA 多批次/多协议整合分析报告",
        f"",
        f"**实验**: `{experiment_id}`",
        f"",
    ]

    # 协议概览
    lines.append("## 协议与批次概览")
    lines.append("")
    lines.append("| Protocol | 数据批次 | 整合方式 |")
    lines.append("|---|---|---|")

    for pid in protocol_ids:
        harmonized_dir = output_root / _INTEGRATION_DIR / _HARMONIZED_SUBDIR / pid
        summary_file = harmonized_dir / f"meta_analysis_summary_{pid}.csv" if harmonized_dir.exists() else None
        if summary_file and summary_file.exists():
            summ = pd.read_csv(summary_file)
            n_batches = int(summ["n_data_batches"].iloc[0]) if len(summ) > 0 else 0
            n_meta = int(summ["n_meta_analyzed"].sum()) if "n_meta_analyzed" in summ else 0
            n_passthrough = int(summ["n_passthrough"].sum()) if "n_passthrough" in summ else 0
            if n_meta > 0:
                method = "随机效应 Meta 分析"
            elif n_passthrough > 0:
                method = "仅单批次覆盖（未形成 pooled effect）"
            else:
                method = "无可合并结果"
        else:
            n_batches = 0
            method = "无数据"
        lines.append(f"| {pid} | {n_batches} | {method} |")

    lines.append("")

    # 各协议 Meta 分析结果
    for pid in protocol_ids:
        harmonized_dir = output_root / _INTEGRATION_DIR / _HARMONIZED_SUBDIR / pid
        if not harmonized_dir.exists():
            continue

        lines.append(f"## Protocol: {pid}")
        lines.append("")

        lines.extend(_design_report_section(harmonized_dir, pid))
        lines.extend(_bridging_section(output_root, pid))

        # 汇总表
        summary_file = harmonized_dir / f"meta_analysis_summary_{pid}.csv"
        if summary_file.exists():
            summ = pd.read_csv(summary_file)
            lines.append("### 汇总统计")
            lines.append("")
            lines.append(summ.to_markdown(index=False))
            lines.append("")

        # 各 comparison 的 top fluxes
        for csv_file in sorted(harmonized_dir.glob("differential_flux_meta_*.csv")):
            comp_name = csv_file.stem.replace("differential_flux_meta_", "")
            df = pd.read_csv(csv_file)

            lines.append(f"### 比较组: {comp_name}")
            lines.append("")
            lines.append(f"- 总通量数: {len(df)}")
            sig = df[df["is_significant_pooled"] == True]
            lines.append(f"- 显著通量数: {len(sig)}")

            n_het = int((df["meta_note"] == "high_heterogeneity").sum())
            if n_het > 0:
                lines.append(f"- 高异质性通量数 (I²>{integration_config.get('within_protocol', {}).get('heterogeneity_threshold', 0.50)*100:.0f}%): {n_het}")

            lines.append("")

            # Top 10 显著上调/下调
            if len(sig) > 0:
                sig_sorted = sig.reindex(sig["diff_pooled"].abs().sort_values(ascending=False).index)
                top = sig_sorted.head(10)
                display_cols = [c for c in ["reaction", "tissue", "diff_pooled", "se_pooled",
                                             "log2FC_pooled", "I_squared", "meta_note"]
                                if c in top.columns]
                lines.append("**Top 显著通量:**")
                lines.append("")
                lines.append(top[display_cols].to_markdown(index=False))
                lines.append("")

    # 跨 protocol 一致性
    scoring_dir = output_root / _INTEGRATION_DIR / _CROSS_TRACER_SUBDIR
    consensus_files = list(scoring_dir.glob("directional_consensus_*.csv")) if scoring_dir.exists() else []
    if consensus_files:
        lines.append("## 跨 Protocol 方向一致性")
        lines.append("")
        for csv_file in sorted(consensus_files):
            comp_name = csv_file.stem.replace("directional_consensus_", "")
            df = pd.read_csv(csv_file)
            n_pass = int(df["consensus_pass"].sum())
            lines.append(f"### {comp_name}")
            lines.append(f"- 总通量数: {len(df)}")
            lines.append(f"- 方向一致通量数: {n_pass}")
            lines.append(f"- 平均一致性分数: {df['directional_consensus_score'].mean():.3f}")
            lines.append("")

            pass_df = df[df["consensus_pass"] == True]
            if len(pass_df) > 0:
                top_pass = pass_df.head(10)
                display_cols = [c for c in ["reaction", "tissue", "direction",
                                            "directional_consensus_score", "n_protocols"]
                                if c in top_pass.columns]
                lines.append("**一致性最高的通量:**")
                lines.append("")
                lines.append(top_pass[display_cols].to_markdown(index=False))
                lines.append("")

    lines.append("---")
    lines.append(f"*报告生成时间: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}*")

    # 写入
    report_path = output_root / _INTEGRATION_DIR / "integration_report.md"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"  [integration] 报告已写入: {report_path}")

    return report_path
