#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
把结果包中的统一结果表导出为 Excel 或 TSV，便于查看与共享（``omicsweave export``）。

Parquet 是结果包的正式格式；导出文件是派生副本，默认写在结果包的 ``exports/`` 下，不登记在运行清单中。
Excel 每个工作表最多 1,048,576 行，超出时请改用 TSV。
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from .dataset import ContractError
from .manifest import MANIFEST_FILE, RunManifest, read_manifest

EXPORT_DIR = "exports"
EXCEL_MAX_ROWS = 1_048_576
FORMATS = ("xlsx", "tsv")


def write_excel_tsv(table: pd.DataFrame, path: Path, **options) -> None:
    """写出便于 Excel 识别编码的制表符文件，保留列值与调用方的精度设置。"""
    table.to_csv(path, sep="\t", encoding="utf-8-sig", index=False, **options)


def find_run(results_root: str | Path, experiment_id: str, run_id: str | None = None) -> Path:
    """实验的结果包目录；不指定 ``run_id`` 时取最近一次运行（运行 ID 按时间排序）。"""
    base = Path(results_root) / experiment_id
    if run_id:
        run_dir = base / run_id
        if not (run_dir / MANIFEST_FILE).is_file():
            raise ContractError(f"{run_dir} 不是完整的结果包（没有 {MANIFEST_FILE}）")
        return run_dir
    runs = sorted(p for p in base.glob("*") if (p / MANIFEST_FILE).is_file()) if base.is_dir() else []
    if not runs:
        raise ContractError(f"{base} 下没有完整的结果包，请先运行 omicsweave run --experiment {experiment_id}")
    return runs[-1]


def bundle_tables(run_dir: str | Path) -> dict[str, pd.DataFrame]:
    """结果包中的结果表：组学结果名为 ``<组学>_<表名>``；整合结果的表名已带模式特有的前缀（如 pathway_fusion），
    直接用表名，与其他表重名时再加 ``<模式>_``。"""
    run_dir = Path(run_dir)
    tables = {f"{path.parent.parent.name}_{path.stem}": pd.read_parquet(path)
              for path in sorted(run_dir.glob("omics/*/tables/*.parquet"))}
    for path in sorted(run_dir.glob("integration/*/tables/*.parquet")):
        name = path.stem if path.stem not in tables else f"{path.parent.parent.name}_{path.stem}"
        tables[name] = pd.read_parquet(path)
    return tables


def _info(manifest: RunManifest, tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = [
        ("run_id", manifest.run_id),
        ("project_id", manifest.project_id),
        ("experiment_id", manifest.experiment_id),
        ("created_at", manifest.created_at.isoformat()),
        ("omicsweave_version", manifest.omicsweave_version),
        ("python_version", manifest.python_version),
        ("command", " ".join(manifest.command)),
        ("config_sha256", manifest.config_sha256 or ""),
        ("omics", ", ".join(o.value for o in manifest.omics)),
        *((f"rows:{name}", str(len(table))) for name, table in tables.items()),
    ]
    return pd.DataFrame(rows, columns=["field", "value"])


def export_bundle(run_dir: str | Path, fmt: str = "xlsx", output: str | Path | None = None) -> Path:
    """导出结果表，返回写出的 Excel 文件或 TSV 目录。"""
    if fmt not in FORMATS:
        raise ContractError(f"导出格式应为 {' / '.join(FORMATS)}")
    run_dir = Path(run_dir)
    manifest = read_manifest(run_dir)
    tables = bundle_tables(run_dir)
    if not tables:
        raise ContractError(f"{run_dir} 中没有结果表（omics/<组学>/tables/ 或 integration/<模式>/tables/）")
    info = _info(manifest, tables)
    stem = f"{manifest.experiment_id}_{manifest.run_id}"

    if fmt == "xlsx":
        too_large = [name for name, table in tables.items() if len(table) + 1 > EXCEL_MAX_ROWS]
        if too_large:
            raise ContractError(f"{'、'.join(too_large)} 超过 Excel 的行数上限，请改用 --format tsv")
        path = Path(output) if output else run_dir / EXPORT_DIR / f"{stem}.xlsx"
        path.parent.mkdir(parents=True, exist_ok=True)
        with pd.ExcelWriter(path) as writer:
            info.to_excel(writer, sheet_name="info", index=False)
            for name, table in tables.items():
                table.to_excel(writer, sheet_name=name[:31], index=False)
        return path

    directory = Path(output) if output else run_dir / EXPORT_DIR / stem
    directory.mkdir(parents=True, exist_ok=True)
    write_excel_tsv(info, directory / "info.tsv")
    for name, table in tables.items():
        write_excel_tsv(table, directory / f"{name}.tsv")
    return directory
