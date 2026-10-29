#!/usr/bin/env python3
"""多组学输入模式、批次来源和结果覆盖审计。"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd


OMICS_INPUTS = (
    ("mfa", "mfa_results"),
    ("transcriptomics", "transcriptomics_results"),
    ("proteomics", "proteomics_results"),
    ("metabolomics", "metabolomics_results"),
)


def _iter_tables(value: Any):
    if isinstance(value, pd.DataFrame):
        yield value
    elif isinstance(value, Mapping):
        for nested in value.values():
            yield from _iter_tables(nested)


def _split_ids(value: Any) -> set[str]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return set()
    return {
        item.strip()
        for item in re.split(r"[;,]", str(value))
        if item.strip()
    }


def _values_from_columns(tables: list[pd.DataFrame], prefixes: tuple[str, ...]) -> set[str]:
    values: set[str] = set()
    for table in tables:
        for column in table.columns:
            if column in prefixes or column.startswith(tuple(f"{prefix}_" for prefix in prefixes)):
                for value in table[column].dropna():
                    values.update(_split_ids(value))
    return values


def _numeric_max(tables: list[pd.DataFrame], prefix: str) -> float:
    values: list[float] = []
    for table in tables:
        for column in table.columns:
            if column == prefix or column.startswith(f"{prefix}_"):
                numeric = pd.to_numeric(table[column], errors="coerce")
                values.extend(numeric[np.isfinite(numeric)].tolist())
    return float(max(values)) if values else np.nan


def build_input_audit(
    config: Any,
    loaded_inputs: Mapping[str, Any],
) -> pd.DataFrame:
    """汇总配置输入与实际载入结果；不会把配置批次数当成实际 pooled 批次数。"""
    integration = getattr(config, "integration", {}) or {}
    input_specs = integration.get("inputs", {}) if isinstance(integration, dict) else {}
    rows: list[dict[str, Any]] = []

    for omics, path_attribute in OMICS_INPUTS:
        spec = dict(input_specs.get(omics, {}) or {})
        batch_specs = list(spec.get("batches", []) or [])
        configured_mode = str(
            spec.get("mode", "integrated" if batch_specs else "single")
        ).lower()
        configured_ids = [
            str(batch.get("batch_id") or f"batch_{index}")
            for index, batch in enumerate(batch_specs, start=1)
        ]
        configured_paths = []
        base_path = getattr(config, path_attribute, None)
        if base_path:
            configured_paths.append(str(Path(base_path)))
        for batch in batch_specs:
            path = batch.get("results_dir") or batch.get("results")
            if path:
                configured_paths.append(str(Path(path)))

        tables = list(_iter_tables(loaded_inputs.get(omics, {})))
        actual_modes = _values_from_columns(tables, ("source_mode",))
        actual_batch_ids = _values_from_columns(tables, ("batch_ids",))
        actual_protocols = _values_from_columns(tables, ("protocol_id",))
        actual_n_batches = _numeric_max(tables, "n_batches")
        if not actual_modes and tables:
            actual_modes = {"single_results_directory"}

        rows.append({
            "omics": omics,
            "configured_mode": configured_mode,
            "configured_n_batches": len(batch_specs),
            "configured_batch_ids": ";".join(configured_ids),
            "configured_paths": ";".join(dict.fromkeys(configured_paths)),
            "n_loaded_tables": len(tables),
            "n_loaded_rows": int(sum(len(table) for table in tables)),
            "actual_source_modes": ";".join(sorted(actual_modes)),
            "actual_n_batches": actual_n_batches,
            "actual_batch_ids": ";".join(sorted(actual_batch_ids)),
            "actual_protocol_ids": ";".join(sorted(actual_protocols)),
            "status": "loaded" if tables else "not_available",
        })
    return pd.DataFrame(rows)


__all__ = ["build_input_audit"]
