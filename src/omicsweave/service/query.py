"""有界、只读的 Parquet 摘要、切片和实体溯源。"""
from __future__ import annotations

import json
import re
from collections import defaultdict
from dataclasses import asdict
from datetime import date, datetime, time, timedelta
from decimal import Decimal
from itertools import islice
from pathlib import Path
from typing import Any, Literal

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.parquet as pq
from pydantic import BaseModel, Field

from omicsweave.contract.manifest import MANIFEST_FILE, read_manifest
from omicsweave.core.delivery import read_delivery_record, resolve_latest_delivery
from omicsweave.core.overrides import TRIALS_DIR, changed_source_outputs, formal_bundle, lineage_for, trial_root
from omicsweave.core.provenance import sha256_file

from .schemas import ErrorCode, ServiceError, json_safe

MAX_ROWS = 200
MAX_BYTES = 120 * 1024  # 留出响应信封余量
MAX_SORT_WINDOW = 20_000
BATCH_SIZE = 1024


class Filter(BaseModel):
    column: str
    op: Literal["eq", "ne", "gt", "gte", "lt", "lte", "in", "contains", "is_null", "not_null"]
    value: Any = None


class TableQuery(BaseModel):
    table: str
    omics: str | None = None
    mode: str | None = None
    filters: list[Filter] = Field(default_factory=list, max_length=32)
    columns: list[str] | None = None
    sort_by: str | None = None
    ascending: bool = True
    limit: int = Field(default=50, ge=1, le=MAX_ROWS)
    offset: int = Field(default=0, ge=0)


def diff_payload(result, output: Path | None = None) -> dict:
    from omicsweave.contract.diff import write_report

    return {"a": str(result.a), "b": str(result.b), "tables": [asdict(table) for table in result.tables],
            "columns": result.columns.to_dict("records"), "transitions": result.transitions.to_dict("records"),
            "other_files": result.other_files, "identical": result.identical,
            **configuration_diff(result.a, result.b),
            "report_dir": str(write_report(result, output)) if output else None}


def configuration_diff(a: Path, b: Path) -> dict:
    def metadata(bundle: Path):
        config_path = _within(bundle / "metadata" / "experiment.json", bundle)
        config = json.loads(config_path.read_text(encoding="utf-8")) if config_path.is_file() else None
        lineage = lineage_for(bundle)
        return config, lineage.model_dump(mode="json") if lineage else None

    def compare(left, right, prefix):
        for key in sorted(left.keys() | right.keys()):
            path = f"{prefix}.{key}" if prefix else key
            old, new = left.get(key), right.get(key)
            if key not in left and isinstance(new, dict) and new:
                yield from compare({}, new, path)
            elif key not in right and isinstance(old, dict) and old:
                yield from compare(old, {}, path)
            elif key not in left:
                yield {"path": path, "op": "add", "after": new}
            elif key not in right:
                yield {"path": path, "op": "remove", "before": old}
            elif isinstance(old, dict) and isinstance(new, dict):
                yield from compare(old, new, path)
            elif type(old) is not type(new) or old != new:
                yield {"path": path, "op": "replace", "before": old, "after": new}

    left, left_lineage = metadata(a)
    right, right_lineage = metadata(b)
    changes = list(compare(left, right, "")) if isinstance(left, dict) and isinstance(right, dict) else None
    return {"config_diff": changes, "config_identical": not changes if changes is not None else None,
            "lineage": {"a": left_lineage, "b": right_lineage}}


def _identifier(value: str) -> str:
    if not value or value in {".", ".."} or any(c in value for c in "/\\:*?[]\0"):
        raise ServiceError(ErrorCode.SECURITY_POLICY_VIOLATION, f"无效的结果标识符 {value!r}")
    return value


def _within(path: Path, parent: Path) -> Path:
    if not path.resolve().is_relative_to(parent.resolve()):
        raise ServiceError(ErrorCode.SECURITY_POLICY_VIOLATION, f"结果文件超出允许的目录 {path}")
    return path


def _open_manifest(directory: Path):
    _within(directory / MANIFEST_FILE, directory)
    try:
        return read_manifest(directory)
    except FileNotFoundError:
        raise
    except (ValueError, OSError) as exc:
        raise ServiceError(ErrorCode.CONTRACT_ERROR, f"无法读取结果清单 {directory}: {exc}") from exc


class ResultView:
    """一次读取固定的交付及其声明的结果包。"""

    def __init__(self, root: Path, experiment: str, run: str, delivery_id: str | None = None,
                 trial_id: str | None = None):
        root = root.resolve()
        formal = None
        if trial_id:
            formal, root = root, trial_root(root, trial_id)
        elif root.parent.name == TRIALS_DIR:
            formal, trial_id = root.parent.parent, root.name
            trial_root(formal, trial_id)
        base = _within(root / _identifier(experiment), root)
        if delivery_id:
            directory = _within(base / _identifier(delivery_id), base)
            _within(directory / "index.json", directory)
            delivery = read_delivery_record(directory, experiment)
            if delivery is None or delivery.status != "completed":
                raise FileNotFoundError(f"未找到已完成的交付 {delivery_id}")
        elif run == "latest":
            delivery, _ = resolve_latest_delivery(root, experiment)
        else:
            delivery = None

        self.experiment = experiment
        self.delivery = delivery
        self.tables: dict[str, Path] = {}
        self.figures: dict[str, Path] = {}
        self.manifests = {}
        self.directories: dict[str, Path] = {}
        self.trial_id = trial_id
        self.lineage = None
        self.reused_source: Path | None = None
        if delivery:
            _within(delivery.directory, base)
            _within(delivery.index_path, delivery.directory)
            candidates = [delivery.directory.name] if trial_id else list(dict.fromkeys(delivery.bundles.values()))
        elif run != "latest":
            candidates = [_identifier(run)]
        else:
            candidates = sorted((p.name for p in base.glob("*")
                                 if (p / MANIFEST_FILE).is_file() and not (p / "index.json").is_file()), reverse=True)[:1]
        if not candidates:
            raise FileNotFoundError(f"{base} 下没有可读取的结果包")

        for name in candidates:
            manifest = self._add_bundle(base, name)
            source = manifest.parameters.get("source_run")
            if trial_id:
                self.lineage = lineage_for(self.directories[name])
                if (self.lineage is None or self.lineage.trial_id != trial_id
                        or self.lineage.experiment_id != experiment or self.lineage.project_id != manifest.project_id):
                    raise ServiceError(ErrorCode.CONTRACT_ERROR, "结果包缺少匹配的试验谱系")
                if self.lineage.source_run_id:
                    if source != self.lineage.source_run_id:
                        raise ServiceError(ErrorCode.CONTRACT_ERROR, "清单与试验复用来源不一致")
                    source_dir = formal_bundle(formal, experiment, source)
                    if sha256_file(source_dir / MANIFEST_FILE) != self.lineage.source_manifest_sha256:
                        raise ServiceError(ErrorCode.STALE_SOURCE_BUNDLE, "试验复用的来源清单已变化")
                    self._add_bundle(source_dir.parent, source)
                    self.reused_source = source_dir
            elif isinstance(source, str) and source != name:
                source_dir = _within(base / _identifier(source), base)
                if (source_dir / MANIFEST_FILE).is_file():
                    self._add_bundle(base, source)
        self.run_id = candidates[0]
        if delivery:
            for path in (delivery.directory / "candidates").glob("*.parquet"):
                key = f"candidates/{path.name}"
                self.tables[key] = _within(path, delivery.directory)

    def _add_bundle(self, base: Path, name: str):
        directory = _within(base / _identifier(name), base)
        manifest = _open_manifest(directory)
        if manifest.experiment_id != self.experiment or manifest.run_id != name:
            raise ServiceError(ErrorCode.CONTRACT_ERROR, f"结果清单与所请求的实验或运行 ID 不匹配 {directory}")
        self.manifests[name] = manifest
        self.directories[name] = directory
        for record in manifest.files:
            if record.role != "output":
                continue
            relative = Path(record.path)
            if relative.suffix.lower() not in {".parquet", ".png", ".pdf", ".svg"}:
                continue
            path = _within(directory / relative, directory)
            key = f"{name}/{relative.as_posix()}"
            if relative.suffix.lower() == ".parquet":
                self.tables[key] = path
            else:
                self.figures[key] = path
        return manifest

    def table(self, request: TableQuery) -> tuple[str, Path]:
        name = request.table if request.table.endswith(".parquet") else f"{request.table}.parquet"
        matches = [(key, path) for key, path in self.tables.items() if
                   (key == name or key.endswith(f"/{name}")) and
                   (not request.omics or f"/omics/{request.omics}/" in key) and
                   (not request.mode or f"/integration/{request.mode}/" in key)]
        if not matches:
            raise FileNotFoundError(f"结果包中没有表 {name}")
        if len(matches) > 1:
            raise ServiceError(ErrorCode.CLI_USAGE_ERROR, f"表名 {name} 有多个匹配，请指定组学、整合模式或完整表名",
                               details={"matches": [key for key, _ in matches][:MAX_ROWS]})
        key, path = matches[0]
        self.check_reused_source(path)
        return key, path

    def check_reused_source(self, only: Path | None = None) -> None:
        if self.reused_source is not None and (only is None or only.is_relative_to(self.reused_source)):
            if changed := changed_source_outputs(self.reused_source, only=only):
                raise ServiceError(ErrorCode.STALE_SOURCE_BUNDLE, "试验复用的科学输出与清单不一致",
                                   details={"changed_files": changed})


def _encoded(value) -> bytes:
    def standard(item):
        if isinstance(item, (datetime, date, time)):
            return item.isoformat()
        if isinstance(item, (Decimal, timedelta)):
            return str(item)
        if isinstance(item, bytes):
            return item.hex()
        raise TypeError(f"不支持的 Parquet 值类型 {type(item).__name__}")

    return json.dumps(json_safe(value), ensure_ascii=True, allow_nan=False, default=standard,
                      separators=(",", ":")).encode("utf-8")


def _page(rows, metadata: dict) -> tuple[list[dict], bool]:
    budget = MAX_BYTES - len(_encoded(metadata)) - 1024
    page = []
    for row in rows:
        encoded = _encoded(row)
        if len(encoded) + 1 > budget:
            if not page:
                raise ServiceError(ErrorCode.CLI_USAGE_ERROR, "单行或结果元数据超过响应限制，请使用 --columns 裁剪")
            return page, True
        budget -= len(encoded) + 1
        page.append(json.loads(encoded))
    return page, False


def _column(schema: pa.Schema, name: str):
    if schema.get_field_index(name) == -1:
        raise ServiceError(ErrorCode.CLI_USAGE_ERROR, f"未知或不唯一的表字段 {name!r}")
    return schema.field(name).type


def _filter(schema: pa.Schema, condition: Filter):
    dtype = _column(schema, condition.column)
    field = ds.field(condition.column)
    value = condition.value
    op = condition.op
    if op in {"is_null", "not_null"}:
        return field.is_null() if op == "is_null" else field.is_valid()
    if op == "contains":
        if not isinstance(value, str) or not (pa.types.is_string(dtype) or pa.types.is_large_string(dtype)):
            raise ServiceError(ErrorCode.CLI_USAGE_ERROR, "contains 仅接受字符串列及字符串值")
        return pc.match_substring(field, value)
    if op == "in":
        if not isinstance(value, list) or not 1 <= len(value) <= 200 or any(v is None for v in value):
            raise ServiceError(ErrorCode.CLI_USAGE_ERROR, "in 需要 1–200 个非空值")
        return field.isin(pa.array(value).cast(dtype, safe=True))
    if value is None or isinstance(value, (list, dict)):
        raise ServiceError(ErrorCode.CLI_USAGE_ERROR, "比较条件需要非空标量；空值请用 is_null")
    scalar = pa.scalar(value)
    if pa.types.is_integer(dtype) and not pa.types.is_integer(scalar.type):
        raise ServiceError(ErrorCode.CLI_USAGE_ERROR, f"列 {condition.column!r} 需要整数条件")
    scalar = scalar.cast(dtype, safe=True)
    return {"eq": field.__eq__, "ne": field.__ne__, "gt": field.__gt__, "gte": field.__ge__,
            "lt": field.__lt__, "lte": field.__le__}[op](scalar)


def query_table(view: ResultView, request: TableQuery) -> dict:
    key, path = view.table(request)
    return read_table_page(view.run_id, key, path, request)


def read_table_page(run_id: str, key: str, path: Path, request: TableQuery) -> dict:
    """分页读取已由调用方授权的结果表。"""
    try:
        dataset = ds.dataset(str(path), format="parquet")
    except pa.ArrowException as exc:
        raise ServiceError(ErrorCode.CONTRACT_ERROR, f"无法读取 Parquet 表 {key}: {exc}") from exc
    columns = request.columns or dataset.schema.names
    for name in columns:
        _column(dataset.schema, name)
    if len(columns) != len(set(columns)):
        raise ServiceError(ErrorCode.CLI_USAGE_ERROR, "列名不能重复")
    condition = None
    try:
        for clause in request.filters:
            expr = _filter(dataset.schema, clause)
            condition = expr if condition is None else condition & expr
        selected = list(columns)
        if request.sort_by:
            _column(dataset.schema, request.sort_by)
            if request.sort_by not in selected:
                selected.append(request.sort_by)
            if request.offset + request.limit > MAX_SORT_WINDOW:
                raise ServiceError(ErrorCode.CLI_USAGE_ERROR, f"排序窗口不能超过 {MAX_SORT_WINDOW} 行，请缩小过滤范围")
        matched = dataset.count_rows(filter=condition)
    except (pa.ArrowException, TypeError, ValueError, OverflowError) as exc:
        raise ServiceError(ErrorCode.CLI_USAGE_ERROR, f"无法应用查询条件 {exc}") from exc

    metadata = {"run_id": run_id, "table_path": key, "total_rows_before_filter": dataset.count_rows(),
                "total_rows_matched": matched, "offset": request.offset, "limit": request.limit, "columns": columns}
    windows = dataset.to_batches(columns=selected, filter=condition, batch_size=BATCH_SIZE)
    needed = request.offset + request.limit
    if request.sort_by:
        # 扫描所有匹配行以保证全局顺序，最多保留 offset + limit 行。
        order_column = "__row_order"
        while order_column in selected:
            order_column += "_"
        keys = [(request.sort_by, "ascending" if request.ascending else "descending"), (order_column, "ascending")]
        top = None
        position = 0
        for batch in windows:
            window = pa.Table.from_batches([batch]).append_column(
                order_column, pa.array(range(position, position + batch.num_rows), type=pa.int64()),
            )
            position += batch.num_rows
            if top is not None:
                window = pa.concat_tables([top, window])
            try:
                top = window.sort_by(keys).slice(0, needed)
            except pa.ArrowNotImplementedError as exc:
                raise ServiceError(ErrorCode.CLI_USAGE_ERROR, f"不支持按 {request.sort_by} 列排序") from exc
        rows = top.slice(request.offset, request.limit).select(columns).to_pylist() if top is not None else []
    else:
        rows = []
        position = 0
        for batch in windows:
            end = position + batch.num_rows
            if end > request.offset:
                start = max(0, request.offset - position)
                rows.extend(batch.slice(start, min(end, needed) - max(position, request.offset)).to_pylist())
            if end >= needed:
                break
            position = end
    page, truncated = _page(rows, metadata)
    next_offset = request.offset + len(page)
    return {**metadata, "rows": page, "truncated": truncated,
            "next_offset": next_offset if page and next_offset < matched else None}


def diagnose_view(view: ResultView, *, limit: int = 100, project_alpha: float | None = None) -> dict:
    from omicsweave.core.diagnostics import inspect_evidence, inspect_registry

    if not 1 <= limit <= MAX_ROWS:
        raise ServiceError(ErrorCode.CLI_USAGE_ERROR, "--limit 必须在 1–200")
    view.check_reused_source()
    warnings = []
    if view.delivery:
        index = _within(view.delivery.index_path, view.delivery.directory)
        saved = json.loads(index.read_text(encoding="utf-8"))
        warnings = [str(warning) for warning in (saved.get("validation_warnings") or [])]
    alpha = view.lineage.project_alpha if view.lineage else project_alpha if project_alpha is not None else 0.05
    findings, availability = inspect_evidence(view.tables, warnings, alpha)
    findings += inspect_registry(view.directories)
    findings.sort(key=lambda item: item.category != "statistical_guardrail")
    blocked = sum(item.category == "statistical_guardrail" for item in findings)
    calibration_status = (
        "partially_assessed" if availability["decoy_strata"] and availability["decoy_unassessed_strata"]
        else "assessed" if availability["decoy_strata"]
        else "unassessed"
    )
    report = {"experiment_id": view.experiment, "run_id": view.run_id, "trial_id": view.trial_id,
              "delivery_id": view.delivery.delivery_id if view.delivery else None,
              "alpha": alpha, "alpha_source": "trial" if view.lineage else
              "current_project_defaults" if project_alpha is not None else "default_screen",
              "evidence_availability": availability, "optimization_blocked": blocked > 0,
              "decoy_calibration_status": calibration_status,
              "blocking_findings": blocked, "total_findings": len(findings)}
    page, _ = _page((item.payload() for item in findings[:limit]), report)
    return {**report, "findings": page, "truncated": len(page) < len(findings)}


def run_digest(view: ResultView) -> dict:
    view.check_reused_source()
    tables, strata, candidates, health = [], [], {}, []
    for key, path in sorted(view.tables.items()):
        parquet = pq.ParquetFile(path)
        schema = parquet.schema_arrow
        tables.append({"name": key, "rows": parquet.metadata.num_rows, "columns": schema.names[:64],
                       "column_count": len(schema.names)})
        stem = path.stem
        if stem == "differential":
            columns = [name for name in ("tissue", "contrast_id", "omics", "pvalue", "selected") if name in schema.names]
            grouped = defaultdict(lambda: {"rows": 0, "n_tested": 0, "n_selected": 0 if "selected" in columns else None})
            if not {"tissue", "contrast_id", "omics"} <= set(columns):
                continue
            for batch in parquet.iter_batches(batch_size=BATCH_SIZE, columns=columns):
                for row in batch.to_pylist():
                    item = grouped[(row["omics"], row["tissue"], row["contrast_id"])]
                    item["rows"] += 1
                    pv = row.get("pvalue")
                    item["n_tested"] += int(pv is not None and 0 <= pv <= 1)
                    if row.get("selected") is True:
                        item["n_selected"] += 1
            strata.extend({"table": key, "omics": layer, "tissue": tissue, "contrast_id": contrast, **counts}
                          for (layer, tissue, contrast), counts in grouped.items())
        elif stem in {"reaction_assessment", "pathway_fusion"}:
            field = next((name for name in ("flux_claim", "direction_agreement") if name in schema.names), None)
            if field:
                counts = defaultdict(int)
                for batch in parquet.iter_batches(batch_size=BATCH_SIZE, columns=[field]):
                    for scalar in batch.column(0):
                        counts[str(scalar.as_py())] += 1
                candidates[key] = {"field": field, "counts": dict(counts)}
        elif stem in {"assessment_calibration", "contrast_reliability"}:
            rows = next(parquet.iter_batches(batch_size=20), None)
            health.append({"table": key, "rows": rows.to_pylist() if rows is not None else [],
                           "total_rows": parquet.metadata.num_rows})
    warnings = []
    if view.delivery:
        index = _within(view.delivery.index_path, view.delivery.directory)
        warnings = json.loads(index.read_text(encoding="utf-8")).get("validation_warnings") or []
    data = {"experiment_id": view.experiment, "run_id": view.run_id, "bundles": list(view.manifests),
            "delivery_id": view.delivery.delivery_id if view.delivery else None,
            **({"trial_id": view.trial_id, "lineage": view.lineage.model_dump(mode="json")} if view.lineage else {})}
    sections = {"tables": tables, "omics_summary": strata, "figures": sorted(view.figures),
                "candidate_overview": [{"table": key, **row} for key, row in candidates.items()],
                "qc_and_calibration_preview": health, "validation_warnings": warnings}
    # 单一预算涵盖全部分区，并报告每组未返回的条目数量。
    flat = ({"section": section, "value": row} for section, rows in sections.items() for row in rows[:MAX_ROWS])
    page, _ = _page(flat, data)
    data.update({section: [] for section in sections})
    for row in page:
        data[row["section"]].append(row["value"])
    data["omitted"] = {section: len(rows) - len(data[section]) for section, rows in sections.items()
                       if len(rows) > len(data[section])}
    return data


_TRACE_FIELDS = {
    "differential": ("feature_id", "feature_name"),
    "derived_differential": ("feature_id", "feature_name", "related_reaction"),
    "differential_by_tracer": ("feature_id", "feature_name"),
    "enrichment": ("hit_features",),
    "dossier_members": ("flux_id", "feature_id", "feature_name", "gene_symbol", "member_id"),
    "reaction_assessment": ("flux_id", "reaction_id"),
    "axis_assessment": ("flux_id", "reaction_id"),
}


def trace_entity(view: ResultView, entity: str, *, tissue: str | None = None, contrast: str | None = None,
                 limit: int = 30) -> dict:
    if not entity.strip() or len(entity) > 256 or not 1 <= limit <= MAX_ROWS:
        raise ServiceError(ErrorCode.CLI_USAGE_ERROR, "实体不能为空或超过 256 字符，limit 必须在 1–200")
    datasets = []
    for key, path in sorted(view.tables.items()):
        fields = _TRACE_FIELDS.get(path.stem)
        if fields is None:
            continue
        view.check_reused_source(path)
        dataset = ds.dataset(str(path), format="parquet")
        expr = None
        for field in fields:
            if field not in dataset.schema.names:
                continue
            dtype = dataset.schema.field(field).type
            # 空的标准表可将整列写成 Arrow null 类型。
            if pa.types.is_null(dtype):
                continue
            value = ds.field(field).cast(pa.string())
            match = pc.match_substring_regex(value, rf"(^|;){re.escape(entity)}(;|$)") \
                if field == "hit_features" else value == entity
            expr = match if expr is None else expr | match
        if expr is None:
            continue
        for field, value in (("tissue", tissue), ("contrast_id", contrast)):
            if value is not None and field in dataset.schema.names:
                expr = expr & (ds.field(field) == value)
        rows = dataset.count_rows(filter=expr)
        if rows:
            datasets.append((key, dataset, expr, rows))
    metadata = {"entity": entity, "run_id": view.run_id, "experiment_id": view.experiment,
                "matches": [{"table": key, "rows": rows} for key, _, _, rows in datasets][:MAX_ROWS]}
    records = ({"table": key, "row": row}
               for key, dataset, expr, _ in datasets
               for row in dataset.head(limit, filter=expr).to_pylist())
    page, truncated = _page(islice(records, limit), metadata)
    return {**metadata, "records": page, "truncated": truncated or sum(item[3] for item in datasets) > len(page)}
