"""
QC 指标采集与阈值门控
======================

此前 ``experiment.yaml`` 里的 ``analysis_defaults.qc`` /
``analysis_defaults.alignment`` 与 ``batch.yaml`` 里的 ``batch_qc`` 都是
「声明但不执行」——没有任何代码读取，QC 纯粹是描述性的，不合格样本
不会被拦下。本模块负责把这些阈值真正落地。

同时把 kallisto 的比对率从 stderr 日志里救出来：``run_info.json``
含结构化的 ``p_pseudoaligned``，之前只被 print 掉，从未落盘，导致
``integration.weighting.quality_signals`` 里声明的 ``mapping_rate``
根本无数据可用。
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

import pandas as pd

# on_qc_fail 支持的策略
POLICY_FLAG = "flag"        # 仅告警，样本仍参与分析
POLICY_EXCLUDE = "exclude"  # 不合格样本移出下游 DE
POLICY_FAIL = "fail"        # 直接中止流程
VALID_POLICIES = (POLICY_FLAG, POLICY_EXCLUDE, POLICY_FAIL)


DEFAULT_QC_THRESHOLDS: dict[str, float] = {
    "min_q20_pct": 0.0,
    "min_q30_pct": 0.0,
    "max_n_content_pct": 100.0,
    "min_pseudoaligned_pct": 0.0,
    "min_total_reads": 0,
}


def read_kallisto_run_info(quant_dir: str | Path, sample_id: str) -> dict[str, Any]:
    """
    读取单样本的 kallisto ``run_info.json``。

    Returns
    -------
    dict
        含 ``n_processed`` / ``n_pseudoaligned`` / ``p_pseudoaligned``；
        文件缺失或损坏时相应值为 None。
    """
    info_path = Path(quant_dir) / sample_id / "run_info.json"
    record: dict[str, Any] = {
        "sample": sample_id,
        "n_processed": None,
        "n_pseudoaligned": None,
        "p_pseudoaligned": None,
        "n_unique": None,
        "kallisto_version": None,
        "index_version": None,
        "run_info_found": False,
    }
    if not info_path.exists():
        return record

    try:
        with open(info_path, encoding="utf-8") as handle:
            doc = json.load(handle)
    except (json.JSONDecodeError, OSError):
        return record

    record.update({
        "n_processed": doc.get("n_processed"),
        "n_pseudoaligned": doc.get("n_pseudoaligned"),
        "p_pseudoaligned": doc.get("p_pseudoaligned"),
        "n_unique": doc.get("n_unique"),
        "kallisto_version": doc.get("kallisto_version"),
        "index_version": doc.get("index_version"),
        "run_info_found": True,
    })
    return record


def collect_quant_summary(quant_dir: str | Path, sample_ids: list[str]) -> pd.DataFrame:
    """汇总所有样本的比对率，返回可直接落盘的表。"""
    rows = [read_kallisto_run_info(quant_dir, sid) for sid in sample_ids]
    return pd.DataFrame(rows)


def summarize_qc_reports(qc_reports: dict[str, dict[str, Any]]) -> pd.DataFrame:
    """
    把每样本 R1/R2 的 QC JSON 压成「每样本一行」的表。

    取 R1/R2 中较差的一侧作为样本级指标——门控应当由最差的一端决定，
    否则 R2 质量崩溃会被 R1 掩盖。
    """
    rows = []
    for sample_id, report in qc_reports.items():
        ends = [report[end] for end in ("R1", "R2") if end in report]
        if not ends:
            continue
        rows.append({
            "sample": sample_id,
            "total_reads": min(e.get("total_reads", 0) for e in ends),
            "mean_quality": min(e.get("mean_quality", 0) for e in ends),
            "q20_pct": min(e.get("q20_pct", 0) for e in ends),
            "q30_pct": min(e.get("q30_pct", 0) for e in ends),
            "gc_content_pct": sum(e.get("gc_content_pct", 0) for e in ends) / len(ends),
            "n_content_pct": max(e.get("n_content_pct", 0) for e in ends),
        })
    return pd.DataFrame(rows)


def evaluate_qc_gate(
    qc_summary: pd.DataFrame,
    quant_summary: pd.DataFrame | None,
    thresholds: dict[str, Any],
) -> pd.DataFrame:
    """
    按阈值逐样本判定，返回门控结果表。

    只检查阈值中实际配置了的项；未配置的项不参与判定（用
    ``DEFAULT_QC_THRESHOLDS`` 的中性值），这样配置文件里省略某项时
    不会意外把所有样本判为不合格。

    Returns
    -------
    pd.DataFrame
        列：sample, qc_pass, qc_failures(以 ";" 连接), 以及参与判定的指标值。
    """
    merged = qc_summary.copy()
    if quant_summary is not None and not quant_summary.empty:
        cols = [c for c in ("sample", "p_pseudoaligned", "n_processed") if c in quant_summary.columns]
        merged = merged.merge(quant_summary[cols], on="sample", how="left")

    effective = dict(DEFAULT_QC_THRESHOLDS)
    for key, value in (thresholds or {}).items():
        if value is not None and key in effective:
            effective[key] = float(value)

    # (阈值键, 数据列, 比较方向) —— "min" 表示低于阈值即不合格
    checks = [
        ("min_q20_pct", "q20_pct", "min"),
        ("min_q30_pct", "q30_pct", "min"),
        ("min_total_reads", "total_reads", "min"),
        ("max_n_content_pct", "n_content_pct", "max"),
        ("min_pseudoaligned_pct", "p_pseudoaligned", "min"),
    ]

    failures: list[str] = []
    for _, row in merged.iterrows():
        reasons = []
        for thr_key, col, direction in checks:
            limit = effective[thr_key]
            if col not in merged.columns:
                continue
            value = row.get(col)
            if value is None or pd.isna(value):
                continue
            if direction == "min" and limit > 0 and float(value) < limit:
                reasons.append(f"{col}={float(value):.2f}<{limit:g}")
            elif direction == "max" and limit < 100.0 and float(value) > limit:
                reasons.append(f"{col}={float(value):.2f}>{limit:g}")
        failures.append(";".join(reasons))

    merged["qc_failures"] = failures
    merged["qc_pass"] = merged["qc_failures"] == ""
    return merged


def apply_qc_policy(
    gate_df: pd.DataFrame,
    policy: str,
) -> tuple[list[str], str]:
    """
    根据 ``on_qc_fail`` 策略决定如何处理不合格样本。

    Returns
    -------
    (excluded_samples, resolved_policy)
        ``excluded_samples`` 仅在 ``exclude`` 策略下非空。
        ``fail`` 策略由调用方在拿到返回值后中止。
    """
    resolved = str(policy or POLICY_FLAG).strip().lower()
    if resolved not in VALID_POLICIES:
        resolved = POLICY_FLAG

    failed = gate_df.loc[~gate_df["qc_pass"], "sample"].astype(str).tolist()
    if resolved == POLICY_EXCLUDE:
        return failed, resolved
    return [], resolved
