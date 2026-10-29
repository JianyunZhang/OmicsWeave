#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
运行溯源：文件指纹、``run_provenance.json`` 的读写，以及跨批次的一致性比对。

跨批次合并的前提是技术链一致：同一份参考、同一版本的工具、同一套代谢网络与示踪剂。
各组学在每次批次运行时写出 ``run_provenance.json``（字段由各组学决定），整合阶段用
:func:`compare_provenance` 比对：

- **关键字段**（critical）不一致或缺少溯源文件，说明批次之间不可比，结论为不一致；
- **提示字段**（advisory）不一致只作提示，例如仪器换代、拟合次数不同。

各组学的字段与采集逻辑见 ``omicsweave.transcriptomics.provenance``、``omicsweave.mfa.provenance``。
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

PROVENANCE_FILE = "run_provenance.json"

# 超过此大小的文件只做「大小 + 首尾分块」指纹，避免整文件哈希拖慢流程（如 kallisto 索引）
FULL_HASH_MAX_BYTES = 512 * 1024 * 1024
_CHUNK_BYTES = 4 * 1024 * 1024


def sha256_file(path: str | Path) -> str:
    """整个文件的 sha256。"""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(_CHUNK_BYTES), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_fingerprint(path: str | Path) -> dict[str, Any]:
    """文件指纹：路径、是否存在、大小与 sha256。

    不超过 :data:`FULL_HASH_MAX_BYTES` 的文件做完整哈希（``hash_mode = sha256_full``）；
    更大的文件哈希「大小 + 首尾各一块」（``sha256_head_tail``），在保持区分度的同时避免读完整个文件。
    """
    p = Path(path)
    record: dict[str, Any] = {"path": str(p), "exists": p.exists()}
    if not p.exists() or not p.is_file():
        return record

    size = p.stat().st_size
    record["size_bytes"] = size
    if size <= FULL_HASH_MAX_BYTES:
        record["hash_mode"] = "sha256_full"
        record["sha256"] = sha256_file(p)
        return record

    digest = hashlib.sha256()
    digest.update(str(size).encode("ascii"))
    with open(p, "rb") as handle:
        digest.update(handle.read(_CHUNK_BYTES))
        handle.seek(max(0, size - _CHUNK_BYTES))
        digest.update(handle.read(_CHUNK_BYTES))
    record["hash_mode"] = "sha256_head_tail"
    record["sha256"] = digest.hexdigest()
    return record


def write_run_provenance(output_dir: str | Path, record: dict[str, Any]) -> Path:
    """把溯源记录写入 ``<output_dir>/run_provenance.json``，返回文件路径。"""
    path = Path(output_dir) / PROVENANCE_FILE
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding="utf-8")
    return path


def load_run_provenance(run_dir: str | Path) -> dict[str, Any] | None:
    """读取某次运行的溯源记录；文件不存在或无法解析时返回 None。"""
    path = Path(run_dir) / PROVENANCE_FILE
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return None


def compare_provenance(
    provenances: dict[str, dict[str, Any] | None],
    critical: Iterable[str],
    advisory: Iterable[str] = (),
) -> dict[str, Any]:
    """比对多个批次的溯源记录。

    Parameters
    ----------
    provenances : dict
        ``{batch_id: 溯源记录或 None}``；None 表示缺少溯源文件。
    critical, advisory : iterable of str
        关键字段与提示字段。

    Returns
    -------
    dict
        - ``consistent``：没有关键字段不一致，且没有批次缺少溯源文件；
        - ``mismatches``：``{字段: {batch_id: 取值}}``，含关键与提示字段；
        - ``severity``：``{字段: "critical" | "advisory"}``；
        - ``missing``：缺少溯源文件的批次（无法判定，视为不一致）；
        - ``n_compared``：参与比对的批次数。
    """
    critical, advisory = tuple(critical), tuple(advisory)
    missing = sorted(bid for bid, record in provenances.items() if not record)
    present = {bid: record for bid, record in provenances.items() if record}

    mismatches: dict[str, dict[str, Any]] = {}
    severity: dict[str, str] = {}
    for field in (*critical, *advisory):
        values = {bid: record.get(field) for bid, record in present.items()}
        distinct = {json.dumps(v, sort_keys=True, ensure_ascii=False) for v in values.values()}
        if len(distinct) > 1:
            mismatches[field] = values
            severity[field] = "critical" if field in critical else "advisory"

    return {
        "consistent": "critical" not in severity.values() and not missing,
        "mismatches": mismatches,
        "severity": severity,
        "missing": missing,
        "n_compared": len(present),
    }


def format_provenance_report(result: dict[str, Any]) -> str:
    """把 :func:`compare_provenance` 的结果渲染成可读文本。"""
    if result.get("consistent") and not result.get("mismatches"):
        return f"参考与技术参数在 {result.get('n_compared', 0)} 个批次间一致。"

    lines = []
    if result.get("missing"):
        lines.append(f"缺少 {PROVENANCE_FILE} 的批次（无法校验）: " + ", ".join(result["missing"]))
    severity = result.get("severity") or {}
    for field, values in (result.get("mismatches") or {}).items():
        rendered = "; ".join(f"{bid}={value!r}" for bid, value in sorted(values.items()))
        note = "（提示）" if severity.get(field) == "advisory" else ""
        lines.append(f"字段 {field} 不一致{note}: {rendered}")
    return "\n".join(lines)
