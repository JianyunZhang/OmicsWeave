#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MFA 跨批次桥接样本一致性检查（技术噪声地板估计）

当同一只动物的同一组织在两个批次里各被拟合一次（不同上机日期、不同
仪器、或仅仅是两次独立的模型运行），两次结果之差**不含任何生物学信号**
——它纯粹是技术 + 拟合的重现性。这个差值的分布就是跨批次比较的
**噪声地板（noise floor）**：任何小于它的跨批次效应量都不可信。

本模块从各批次的 ``01_df_results/all_tissues/flux_matrix/sample_net_flux.csv``
中找出重复出现的 (sample_id, flux_id)，逐对量化差异，并输出：

- ``bridge_sample_overview.csv``  哪些样本 / 组织 / 通量在多个批次中重复
- ``bridge_flux_pairs.csv``       每个桥接对的两次估计与其对数比
- ``bridge_noise_floor.csv``      每个批次对的噪声地板汇总

设计要点
--------
桥接样本不是"多出来的重复"，而是**唯一可测量的技术重现性来源**。若实验
里本来没有刻意留桥接样本，任何被两个批次共同覆盖的组织都会被自动识别
出来充当这一角色。

**前提：sample_id 必须能跨批次唯一标识同一只动物。** 上游
``sheet_information.complete_result_name`` 若写成
``glucose__male__hfd_sham_8w__1`` 这种「组名+序号」标签，两个批次里的
同名样本未必是同一只鼠——此时按同名配对得到的差异里混着生物学个体差异，
**不能**当作技术噪声地板。``classify_sample_identity()`` 会自动识别这种
情况，把结果标注为 ``upper_bound_confounded_with_biology`` 并在日志与
报告中明确说明，而不是给出一个看起来像噪声地板的数字。

对数比只在两次估计**同号且量级足够**时计算——净通量在零附近时
``log2(a/b)`` 会爆炸，那属于数值噪声而非技术噪声。被排除的对会单独计数
并报告，不会被悄悄丢弃。
"""

from __future__ import annotations

from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

_INTEGRATION_DIR = "integration"
_BRIDGING_SUBDIR = "bridging_qc"

# 相对量级门槛：两次估计的较大绝对值必须达到该批次对中位绝对通量的这个比例，
# 否则对数比不可解释（分母接近零）。
_DEFAULT_MAGNITUDE_FRACTION = 0.05

# 噪声地板取对数比绝对值的这个分位数。
_DEFAULT_FLOOR_QUANTILE = 0.95


IDENTITY_ANIMAL = "explicit_sample_id"
IDENTITY_GROUP_INDEX = "group_index_label"


def classify_sample_identity(matrices: dict[str, pd.DataFrame]) -> str:
    """判断 ``sample_id`` 到底是动物标识，还是「组名+序号」标签。

    上游 ``sheet_information.complete_result_name`` 常写成
    ``glucose__male__hfd_sham_8w__1`` 这种形式——它编码的是**组别和组内
    序号**，不是动物 ID。两个批次里的 ``hfd_sham_8w__1`` 未必是同一只鼠，
    此时按 sample_id 配对得到的差异里混着生物学个体差异，**不能**当作
    技术噪声地板。

    判据：sample_id 中是否包含该样本自己的 condition 字符串。
    """
    total = 0
    encoded = 0
    for df in matrices.values():
        if "condition" not in df.columns:
            continue
        pairs = df[["sample_id", "condition"]].drop_duplicates()
        for sample_id, condition in pairs.itertuples(index=False):
            condition = str(condition).strip().lower()
            if not condition or condition == "unknown":
                continue
            total += 1
            if condition in str(sample_id).strip().lower():
                encoded += 1
    if total == 0:
        return IDENTITY_ANIMAL
    return IDENTITY_GROUP_INDEX if encoded / total >= 0.8 else IDENTITY_ANIMAL


def _sample_flux_path(output_root: Path, protocol_id: str, batch_id: str) -> Path:
    return (
        output_root / "runs" / protocol_id / batch_id / "01_df_results"
        / "all_tissues" / "flux_matrix" / "sample_net_flux.csv"
    )


def load_sample_flux_matrices(
    output_root: Path,
    protocol_id: str,
    batch_ids: list[str],
) -> dict[str, pd.DataFrame]:
    """读取各批次的样本级净通量表。缺失的批次会被跳过并提示。"""
    matrices: dict[str, pd.DataFrame] = {}
    for batch_id in batch_ids:
        path = _sample_flux_path(output_root, protocol_id, batch_id)
        if not path.exists():
            print(f"    [bridge] {batch_id}: 缺少 sample_net_flux.csv"
                  f"（该批次需以当前版本重跑）")
            continue
        df = pd.read_csv(path)
        required = {"sample_id", "flux_id", "net_flux"}
        if not required.issubset(df.columns):
            print(f"    [bridge] {batch_id}: sample_net_flux.csv 缺少必需列 "
                  f"{sorted(required - set(df.columns))}")
            continue
        df["sample_id"] = df["sample_id"].astype(str)
        matrices[batch_id] = df
    return matrices


def apply_tissue_aliases(
    matrices: dict[str, pd.DataFrame],
    aliases: dict[str, str],
) -> dict[str, pd.DataFrame]:
    """把同一隔室的不同叫法归一，否则 ``serum--X`` 与 ``plasma--X``
    永远配不上对，该隔室就白白失去了桥接价值。"""
    if not aliases:
        return matrices
    lowered = {str(k).lower(): str(v) for k, v in aliases.items()}
    out = {}
    for batch_id, df in matrices.items():
        if "tissue" not in df.columns:
            out[batch_id] = df
            continue
        df = df.copy()
        original = df["tissue"].astype(str)
        canonical = original.str.lower().map(lowered).fillna(original)
        df["tissue"] = canonical
        if "reaction" in df.columns:
            df["flux_id"] = canonical + "--" + df["reaction"].astype(str)
        out[batch_id] = df
    return out


def _pair_frame(
    df_a: pd.DataFrame,
    df_b: pd.DataFrame,
    batch_a: str,
    batch_b: str,
) -> pd.DataFrame:
    """构造两个批次之间的 (sample_id, flux_id) 桥接对。"""
    keys = ["sample_id", "flux_id"]
    cols = keys + ["condition", "tissue", "reaction", "net_flux", "fit_std"]
    left = df_a[[c for c in cols if c in df_a.columns]].copy()
    right = df_b[[c for c in cols if c in df_b.columns]].copy()

    merged = left.merge(right, on=keys, how="inner", suffixes=("_a", "_b"))
    if merged.empty:
        return merged

    merged["batch_a"] = batch_a
    merged["batch_b"] = batch_b
    for base in ("condition", "tissue", "reaction"):
        col_a, col_b = f"{base}_a", f"{base}_b"
        if col_a in merged.columns:
            merged[base] = merged[col_a]
            merged = merged.drop(columns=[c for c in (col_a, col_b) if c in merged.columns])

    merged["abs_diff"] = merged["net_flux_a"] - merged["net_flux_b"]
    denominator = merged[["net_flux_a", "net_flux_b"]].abs().max(axis=1)
    merged["rel_diff"] = np.where(
        denominator > 0, merged["abs_diff"] / denominator, np.nan
    )
    merged["sign_agree"] = (
        np.sign(merged["net_flux_a"]) == np.sign(merged["net_flux_b"])
    )
    return merged


def _annotate_log2_ratio(
    pairs: pd.DataFrame,
    magnitude_fraction: float,
) -> tuple[pd.DataFrame, float]:
    """在同号且量级足够的对上计算 log2 比，返回 (标注后的表, 量级门槛)。"""
    magnitudes = pairs[["net_flux_a", "net_flux_b"]].abs().to_numpy().ravel()
    magnitudes = magnitudes[np.isfinite(magnitudes)]
    median_magnitude = float(np.median(magnitudes)) if len(magnitudes) else 0.0
    magnitude_floor = magnitude_fraction * median_magnitude

    large_enough = (
        pairs[["net_flux_a", "net_flux_b"]].abs().min(axis=1) >= magnitude_floor
    )
    usable = pairs["sign_agree"] & large_enough & (pairs["net_flux_b"] != 0)

    ratio = np.where(
        usable,
        np.abs(pairs["net_flux_a"]) / np.abs(pairs["net_flux_b"]).replace(0, np.nan),
        np.nan,
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        pairs = pairs.copy()
        pairs["log2_ratio"] = np.log2(ratio)
    pairs["usable_for_noise_floor"] = usable
    pairs.loc[~usable, "log2_ratio"] = np.nan
    pairs["exclusion_reason"] = np.where(
        usable,
        "",
        np.where(~pairs["sign_agree"], "sign_reversal", "below_magnitude_floor"),
    )
    return pairs, magnitude_floor


def bridging_sample_qc(
    output_root: Path,
    protocol_id: str,
    batch_ids: list[str],
    integration_config: dict[str, Any] | None = None,
) -> pd.DataFrame | None:
    """跨批次桥接样本一致性检查。

    Parameters
    ----------
    output_root : Path
        ``mfa_results`` 根目录。
    protocol_id : str
    batch_ids : list[str]
        该 protocol 下有数据的批次。
    integration_config : dict, optional
        读取 ``bridging`` 段：``magnitude_fraction``、``floor_quantile``。

    Returns
    -------
    pd.DataFrame or None
        噪声地板汇总表；无桥接样本时返回 None。
    """
    config = (integration_config or {}).get("bridging", {}) or {}
    magnitude_fraction = float(
        config.get("magnitude_fraction", _DEFAULT_MAGNITUDE_FRACTION)
    )
    floor_quantile = float(config.get("floor_quantile", _DEFAULT_FLOOR_QUANTILE))
    # 显式声明 sample_id 是否可跨批次唯一标识同一只动物；
    # 不声明时自动判别。
    declared_identity = config.get("sample_id_is_animal_id", None)

    if len(batch_ids) < 2:
        print(f"    [bridge] protocol={protocol_id}: 仅 1 个批次，无跨批次桥接可查")
        return None

    matrices = load_sample_flux_matrices(output_root, protocol_id, batch_ids)
    if len(matrices) < 2:
        print(f"    [bridge] protocol={protocol_id}: 可用的样本级通量表 < 2，跳过")
        return None

    matrices = apply_tissue_aliases(
        matrices, (integration_config or {}).get("tissue_aliases", {}) or {}
    )

    if declared_identity is None:
        identity_basis = classify_sample_identity(matrices)
    else:
        identity_basis = IDENTITY_ANIMAL if bool(declared_identity) else IDENTITY_GROUP_INDEX

    if identity_basis == IDENTITY_GROUP_INDEX:
        interpretation = "upper_bound_confounded_with_biology"
        print("    [bridge] !! sample_id 是「组名+序号」标签，不是动物 ID。")
        print("             两个批次里的同名样本未必是同一只动物，因此下面的数字")
        print("             **不是技术噪声地板**，而是「技术 + 生物学个体差异 +")
        print("             可能的配对错误」的上界。要得到真正的噪声地板，需要在")
        print("             complete_result_name 中写入动物 ID，或纳入真正的桥接样本。")
    else:
        interpretation = "technical_noise_floor"

    out_dir = output_root / _INTEGRATION_DIR / _BRIDGING_SUBDIR / protocol_id
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── 样本覆盖总览 ──────────────────────────────────────────────────────
    overview_rows = []
    for batch_id, df in matrices.items():
        for sample_id, group in df.groupby("sample_id"):
            overview_rows.append({
                "protocol_id": protocol_id,
                "batch_id": batch_id,
                "sample_id": sample_id,
                "condition": group["condition"].iloc[0] if "condition" in group else "",
                "n_tissues": group["tissue"].nunique() if "tissue" in group else np.nan,
                "n_fluxes": group["flux_id"].nunique(),
            })
    overview = pd.DataFrame(overview_rows)
    if not overview.empty:
        batches_per_sample = overview.groupby("sample_id")["batch_id"].nunique()
        overview["n_batches_for_sample"] = overview["sample_id"].map(batches_per_sample)
        overview["is_bridging_sample"] = overview["n_batches_for_sample"] >= 2
        overview.to_csv(out_dir / "bridge_sample_overview.csv", index=False)

    # ── 逐批次对构造桥接 ──────────────────────────────────────────────────
    all_pairs = []
    floor_rows = []

    for batch_a, batch_b in combinations(sorted(matrices), 2):
        pairs = _pair_frame(matrices[batch_a], matrices[batch_b], batch_a, batch_b)
        if pairs.empty:
            print(f"    [bridge] {batch_a} x {batch_b}: 无共同 (sample_id, flux_id)")
            floor_rows.append({
                "protocol_id": protocol_id,
                "batch_a": batch_a,
                "batch_b": batch_b,
                "n_bridge_samples": 0,
                "n_bridge_fluxes": 0,
                "n_pairs": 0,
                "n_usable_pairs": 0,
                "sign_agreement_rate": np.nan,
                "median_abs_log2_ratio": np.nan,
                "p75_abs_log2_ratio": np.nan,
                "noise_floor_log2fc": np.nan,
                "floor_quantile": floor_quantile,
                "magnitude_floor": np.nan,
                "identity_basis": identity_basis,
                "interpretation": interpretation,
                "note": "no_shared_sample_flux_pairs",
            })
            continue

        pairs, magnitude_floor = _annotate_log2_ratio(pairs, magnitude_fraction)
        all_pairs.append(pairs)

        usable = pairs[pairs["usable_for_noise_floor"]]
        abs_log2 = usable["log2_ratio"].abs().dropna()

        floor_rows.append({
            "protocol_id": protocol_id,
            "batch_a": batch_a,
            "batch_b": batch_b,
            "n_bridge_samples": int(pairs["sample_id"].nunique()),
            "n_bridge_fluxes": int(pairs["flux_id"].nunique()),
            "n_pairs": int(len(pairs)),
            "n_usable_pairs": int(len(abs_log2)),
            "sign_agreement_rate": float(pairs["sign_agree"].mean()),
            "median_abs_log2_ratio": float(abs_log2.median()) if len(abs_log2) else np.nan,
            "p75_abs_log2_ratio": float(abs_log2.quantile(0.75)) if len(abs_log2) else np.nan,
            "noise_floor_log2fc": (
                float(abs_log2.quantile(floor_quantile)) if len(abs_log2) else np.nan
            ),
            "floor_quantile": floor_quantile,
            "magnitude_floor": magnitude_floor,
            "identity_basis": identity_basis,
            "interpretation": interpretation,
            "note": "ok" if len(abs_log2) else "no_usable_pairs",
        })

    if not floor_rows:
        return None

    if all_pairs:
        pairs_df = pd.concat(all_pairs, ignore_index=True)
        pairs_df = pairs_df.reindex(
            pairs_df["log2_ratio"].abs().sort_values(ascending=False, na_position="last").index
        ).reset_index(drop=True)
        pairs_df.to_csv(out_dir / "bridge_flux_pairs.csv", index=False)

    floor_df = pd.DataFrame(floor_rows)
    floor_df.to_csv(out_dir / "bridge_noise_floor.csv", index=False)

    # ── 控制台摘要 ────────────────────────────────────────────────────────
    for _, row in floor_df.iterrows():
        if row["n_pairs"] == 0:
            continue
        print(f"    [bridge] {row['batch_a']} x {row['batch_b']}: "
              f"{int(row['n_bridge_samples'])} 个桥接样本 x "
              f"{int(row['n_bridge_fluxes'])} 条通量 = {int(row['n_pairs'])} 对")
        print(f"             方向一致率 {row['sign_agreement_rate']:.1%}, "
              f"|log2 ratio| 中位数 {row['median_abs_log2_ratio']:.3f}")
        if pd.notna(row["noise_floor_log2fc"]):
            floor = float(row["noise_floor_log2fc"])
            label = ("噪声地板" if interpretation == "technical_noise_floor"
                     else "重现性上界（含生物学差异，非技术噪声地板）")
            print(f"             {label} (P{floor_quantile * 100:.0f}) = "
                  f"|log2FC| {floor:.3f}  (≈ {2 ** floor:.2f}x)")

    print(f"    [bridge] → {out_dir}")
    return floor_df
