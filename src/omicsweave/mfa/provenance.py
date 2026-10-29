#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MFA 运行溯源与跨批次可合并性断言

两个批次的差异通量能不能放进同一个 Meta 分析，取决于它们背后的技术链是否
一致：同一种示踪剂、同一套代谢网络、同一类提取方法。这些信息此前只以
注释形式躺在 ``batch.yaml`` 里（``detection_device``、
``metabolites_extraction_protocal``），没有任何代码读取，也就无从阻止把
两条不可比的技术链悄悄合并。

本模块为每次批次运行写出 ``run_provenance.json``，并在整合阶段比对各批次的
指纹，对不一致项按 ``integration.provenance_mismatch_policy``
（``fail`` | ``warn`` | ``ignore``，默认 ``warn``）作出反应。

与转录组的 ``run_provenance.json`` 作用相同。
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime
from pathlib import Path
from typing import Any

_PROVENANCE_FILE = "run_provenance.json"

# 这些字段不一致意味着两个批次测的不是同一件事，合并前必须人工确认。
_CRITICAL_FIELDS = [
    "tracer_type",
    "tracer_isotope",
    "reaction_set_hash",
    "extraction_protocol",
]
# 这些字段不一致只是提示信息（仪器换代、操作人不同等）。
_ADVISORY_FIELDS = [
    "detection_device",
    "detection_division",
    "n_fit_iterations",
]


def _reaction_set_hash(reactions) -> str:
    """代谢网络指纹：反应名集合的稳定哈希。

    网络拓扑不同的两次拟合，同名通量的数值含义也不同——合并它们是无声的
    错误，因此把反应集合本身纳入指纹。
    """
    unique = sorted({str(r) for r in reactions})
    digest = hashlib.sha256("|".join(unique).encode("utf-8")).hexdigest()
    return digest[:16]


def write_run_provenance(
    output_dir: Path | str,
    *,
    experiment_id: str,
    protocol_id: str | None,
    batch_id: str | None,
    batch_metadata: dict[str, Any] | None,
    tracer: dict[str, Any] | None,
    reactions,
    tissues,
    conditions,
    n_samples: int,
    n_fit_iterations: int | None,
    analysis_config: dict[str, Any] | None,
) -> Path:
    """写出本次批次运行的溯源指纹。"""
    metadata = dict(batch_metadata or {})
    tracer = dict(tracer or {})
    diff_cfg = (analysis_config or {}).get("differential_flux", {}) or {}

    record = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "experiment_id": experiment_id,
        "protocol_id": protocol_id,
        "batch_id": batch_id,
        "tracer_type": tracer.get("type"),
        "tracer_isotope": tracer.get("isotope"),
        "tracer_enrichment_atom_percent": tracer.get("isotopic_enrichment_atom_percent"),
        "detection_device": metadata.get("detection_device"),
        "detection_division": metadata.get("detection_division"),
        "detection_date": str(metadata.get("detection_date", "")) or None,
        "extraction_protocol": metadata.get("metabolites_extraction_protocal"),
        "reaction_set_hash": _reaction_set_hash(reactions),
        "n_reactions": len({str(r) for r in reactions}),
        "tissues": sorted({str(t) for t in tissues}),
        "conditions": sorted({str(c) for c in conditions}),
        "n_samples": int(n_samples),
        "n_fit_iterations": n_fit_iterations,
        "analysis_thresholds": {
            "statistic": diff_cfg.get("statistic"),
            "pvalue_threshold": diff_cfg.get("pvalue_threshold"),
            "fdr_threshold": diff_cfg.get("fdr_threshold"),
            "log2fc_threshold": diff_cfg.get("log2fc_threshold"),
            "cohen_d_threshold": diff_cfg.get("cohen_d_threshold"),
        },
    }

    path = Path(output_dir) / _PROVENANCE_FILE
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(record, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"  Saved {_PROVENANCE_FILE} "
          f"(network {record['reaction_set_hash']}, "
          f"{record['n_reactions']} reactions, "
          f"extraction={record['extraction_protocol']})")
    return path


def read_run_provenance(
    output_root: Path,
    protocol_id: str,
    batch_id: str,
) -> dict[str, Any] | None:
    path = output_root / "runs" / protocol_id / batch_id / _PROVENANCE_FILE
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def assert_batch_compatibility(
    output_root: Path,
    protocol_id: str,
    batch_ids: list[str],
    policy: str = "warn",
) -> list[dict[str, Any]]:
    """比对同一 protocol 下各批次的技术链指纹。

    Returns
    -------
    list of dict
        每个不一致字段一条记录（字段名、各批次取值、严重程度）。

    Raises
    ------
    RuntimeError
        当 ``policy == "fail"`` 且存在关键字段不一致时。
    """
    policy = str(policy or "warn").lower()
    if policy == "ignore":
        return []

    records = {}
    missing = []
    for batch_id in batch_ids:
        record = read_run_provenance(output_root, protocol_id, batch_id)
        if record is None:
            missing.append(batch_id)
        else:
            records[batch_id] = record

    if missing:
        print(f"    [provenance] {len(missing)} 个批次缺少 {_PROVENANCE_FILE}"
              f"（需以当前版本重跑）: {missing}")
    if len(records) < 2:
        return []

    mismatches = []
    for field in _CRITICAL_FIELDS + _ADVISORY_FIELDS:
        values = {b: r.get(field) for b, r in records.items()}
        distinct = {json.dumps(v, ensure_ascii=False, sort_keys=True)
                    for v in values.values()}
        if len(distinct) <= 1:
            continue
        mismatches.append({
            "field": field,
            "severity": "critical" if field in _CRITICAL_FIELDS else "advisory",
            "values": values,
        })

    if not mismatches:
        print("    [provenance] 各批次技术链指纹一致")
        return []

    critical = [m for m in mismatches if m["severity"] == "critical"]
    for mismatch in mismatches:
        tag = "!!" if mismatch["severity"] == "critical" else " -"
        detail = ", ".join(f"{b}={v}" for b, v in mismatch["values"].items())
        print(f"    [provenance] {tag} {mismatch['field']} 不一致: {detail}")

    if critical:
        print("    [provenance] 关键字段不一致意味着两个批次测的不是同一件事，")
        print("                 合并结果需人工确认后才可使用。")
        if policy == "fail":
            fields = ", ".join(m["field"] for m in critical)
            raise RuntimeError(
                f"批次技术链不兼容（{fields}）；"
                f"设置 integration.provenance_mismatch_policy: warn 可降级为提示。"
            )

    return mismatches
