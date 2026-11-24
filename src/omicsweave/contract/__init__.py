"""
数据契约：分析层数据集、统一结果表与运行清单的模型、读写和校验。

各组学模块只通过这一层读写分析数据；同一套模型导出的 JSON Schema 供下游（如 TheraWeave）使用。
交换契约 v1 在阶段 2 定稿，当前为草案（``BUNDLE_SCHEMA_VERSION``）。
"""
from typing import Any

from .dataset import ContractError, Dataset, DatasetMeta, ValueType, empty_features, read_dataset, write_dataset
from .manifest import (
    BUNDLE_SCHEMA_VERSION,
    FileRecord,
    RunManifest,
    allocate_run_dir,
    file_record,
    new_run_id,
    read_manifest,
    write_manifest,
)
from .tables import (
    DIFFERENTIAL_COLUMNS,
    ENRICHMENT_COLUMNS,
    read_differential,
    read_enrichment,
    validate_differential,
    validate_enrichment,
    write_differential,
    write_enrichment,
)


def export_json_schema() -> dict[str, dict[str, Any]]:
    """返回数据集元数据与运行清单的 JSON Schema。"""
    return {
        "dataset": DatasetMeta.model_json_schema(),
        "manifest": RunManifest.model_json_schema(),
    }


__all__ = [
    "BUNDLE_SCHEMA_VERSION",
    "ContractError",
    "DIFFERENTIAL_COLUMNS",
    "ENRICHMENT_COLUMNS",
    "Dataset",
    "DatasetMeta",
    "FileRecord",
    "RunManifest",
    "ValueType",
    "allocate_run_dir",
    "empty_features",
    "export_json_schema",
    "file_record",
    "new_run_id",
    "read_dataset",
    "read_differential",
    "read_enrichment",
    "read_manifest",
    "validate_differential",
    "validate_enrichment",
    "write_dataset",
    "write_differential",
    "write_enrichment",
    "write_manifest",
]
