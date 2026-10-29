"""
转录组运行溯源
==============

跨批次合并的前提是**技术链路一致**。如果两个批次用了不同版本的转录组
参考、不同的 kallisto 版本、或不同的链特异性设置，它们的基因层计数在
定义上就不可比——此时做 Meta 分析只会把技术差异当成生物学差异。

``capture_run_provenance()`` 在每次运行时把实际使用的参考文件指纹、
工具版本和技术参数写入 ``run_provenance.json``；整合阶段用
:func:`omicsweave.core.provenance.compare_provenance` 按
:data:`CRITICAL_PROVENANCE_FIELDS` 比对各批次，发现不一致时由调用方按策略拒绝整合。
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any

from omicsweave.core.provenance import file_fingerprint, write_run_provenance

# 参与跨批次一致性判定的字段（技术链路的核心）
CRITICAL_PROVENANCE_FIELDS = (
    "transcriptome_fingerprint",
    "kallisto_version",
    "strandedness",
    "gene_aggregation",
)


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

    out_path = write_run_provenance(output_dir, provenance)
    provenance["_path"] = str(out_path)
    return provenance
