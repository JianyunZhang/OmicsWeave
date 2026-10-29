"""
运行溯源信息采集与跨批次一致性校验
====================================

跨批次合并的前提是**技术链路一致**。如果两个批次用了不同版本的转录组
参考、不同的 kallisto 版本、或不同的链特异性设置，它们的基因层计数在
定义上就不可比——此时做 Meta 分析只会把技术差异当成生物学差异。

``capture_run_provenance()`` 在每次运行时把实际使用的参考文件指纹、
工具版本和技术参数写入 ``run_provenance.json``；
``compare_provenance()`` 在整合阶段比对各批次，发现不一致时由调用方
按策略拒绝整合。
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
from pathlib import Path
from typing import Any

# 超过此大小的文件只做「首尾分块 + 大小」指纹，避免整文件哈希拖慢流程
_FULL_HASH_MAX_BYTES = 512 * 1024 * 1024
_CHUNK_BYTES = 4 * 1024 * 1024

# 参与跨批次一致性判定的字段（技术链路的核心）
CRITICAL_PROVENANCE_FIELDS = (
    "transcriptome_fingerprint",
    "kallisto_version",
    "strandedness",
    "gene_aggregation",
)


def file_fingerprint(path: str | Path) -> dict[str, Any]:
    """
    计算文件指纹。

    小文件用完整 sha256；超大文件（如 kallisto index）用「大小 + 首尾
    分块哈希」，在保持区分度的同时避免读完整个文件。
    """
    p = Path(path)
    record: dict[str, Any] = {"path": str(p), "exists": p.exists()}
    if not p.exists() or not p.is_file():
        return record

    size = p.stat().st_size
    record["size_bytes"] = size

    digest = hashlib.sha256()
    if size <= _FULL_HASH_MAX_BYTES:
        record["hash_mode"] = "sha256_full"
        with open(p, "rb") as handle:
            for chunk in iter(lambda: handle.read(_CHUNK_BYTES), b""):
                digest.update(chunk)
    else:
        record["hash_mode"] = "sha256_head_tail"
        digest.update(str(size).encode("ascii"))
        with open(p, "rb") as handle:
            digest.update(handle.read(_CHUNK_BYTES))
            handle.seek(max(0, size - _CHUNK_BYTES))
            digest.update(handle.read(_CHUNK_BYTES))

    record["sha256"] = digest.hexdigest()
    return record


def detect_kallisto_version(kallisto_exe: str | Path) -> str | None:
    """调用 ``kallisto version`` 获取版本串；不可用时返回 None。"""
    try:
        result = subprocess.run(
            [str(kallisto_exe), "version"],
            capture_output=True, text=True, timeout=30,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    text = (result.stdout or "").strip() or (result.stderr or "").strip()
    return text or None


def capture_run_provenance(
    output_dir: str | Path,
    *,
    experiment_id: str,
    protocol_id: str | None,
    batch_id: str | None,
    transcriptome_fasta: str | Path,
    kallisto_exe: str | Path,
    kallisto_index: str | Path,
    strandedness: str,
    gene_aggregation: str,
    technical_params: dict[str, Any] | None = None,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """采集本次运行的溯源信息并写入 ``run_provenance.json``。"""
    provenance: dict[str, Any] = {
        "experiment_id": experiment_id,
        "protocol_id": protocol_id,
        "batch_id": batch_id,
        "transcriptome_fingerprint": file_fingerprint(transcriptome_fasta).get("sha256"),
        "transcriptome_path": str(transcriptome_fasta),
        "kallisto_index_path": str(kallisto_index),
        "kallisto_version": detect_kallisto_version(kallisto_exe),
        "strandedness": strandedness,
        "gene_aggregation": gene_aggregation,
        "technical_params": dict(technical_params or {}),
    }
    if extra:
        provenance.update(extra)

    out_path = Path(output_dir) / "run_provenance.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as handle:
        json.dump(provenance, handle, indent=2, ensure_ascii=False)

    provenance["_path"] = str(out_path)
    return provenance


def load_run_provenance(run_dir: str | Path) -> dict[str, Any] | None:
    """读取某次运行的溯源信息；不存在返回 None。"""
    path = Path(run_dir) / "run_provenance.json"
    if not path.exists():
        return None
    try:
        with open(path, encoding="utf-8") as handle:
            return json.load(handle)
    except (json.JSONDecodeError, OSError):
        return None


def compare_provenance(
    provenances: dict[str, dict[str, Any] | None],
    fields: tuple[str, ...] = CRITICAL_PROVENANCE_FIELDS,
) -> dict[str, Any]:
    """
    比对多个批次的溯源信息。

    Parameters
    ----------
    provenances : dict
        ``{batch_id: provenance_dict_or_None}``。

    Returns
    -------
    dict
        ``consistent`` 为整体结论；``mismatches`` 列出不一致字段及各批次取值；
        ``missing`` 列出缺少溯源文件的批次（无法判定，视为不一致）。
    """
    missing = sorted([bid for bid, prov in provenances.items() if not prov])
    present = {bid: prov for bid, prov in provenances.items() if prov}

    mismatches: dict[str, dict[str, Any]] = {}
    for field in fields:
        values = {bid: prov.get(field) for bid, prov in present.items()}
        distinct = {json.dumps(v, sort_keys=True, ensure_ascii=False) for v in values.values()}
        if len(distinct) > 1:
            mismatches[field] = values

    return {
        "consistent": not mismatches and not missing,
        "mismatches": mismatches,
        "missing": missing,
        "n_compared": len(present),
    }


def format_provenance_report(result: dict[str, Any]) -> str:
    """把 ``compare_provenance()`` 的结果渲染成可读文本。"""
    if result.get("consistent"):
        return f"参考与技术参数在 {result.get('n_compared', 0)} 个批次间一致。"

    lines = []
    if result.get("missing"):
        lines.append(
            "缺少 run_provenance.json 的批次（无法校验）: "
            + ", ".join(result["missing"])
        )
    for field, values in (result.get("mismatches") or {}).items():
        rendered = "; ".join(f"{bid}={value!r}" for bid, value in sorted(values.items()))
        lines.append(f"字段 {field} 不一致: {rendered}")
    return "\n".join(lines)
