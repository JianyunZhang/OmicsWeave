"""项目发现、校验和事务性的配置与样本 QC 接口。"""
from __future__ import annotations

import csv
import difflib
import hashlib
import io
import re
from dataclasses import asdict, replace
from pathlib import Path
from typing import Callable, Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field

from omicsweave.core.config import EXPERIMENT_FILE
from omicsweave.core.file_state import file_lock, replace_bytes
from omicsweave.core.path_scope import PathScope, file_component

from .schemas import ErrorCode, ServiceError


class EditOptions(BaseModel):
    model_config = ConfigDict(extra="forbid")

    dry_run: bool = True
    reason: str = Field(default="", max_length=2000)
    expected_sha256: str | None = Field(default=None, pattern=r"^[0-9a-f]{64}$")


class QCUpdate(BaseModel):
    model_config = ConfigDict(extra="forbid")

    sample_id: str = Field(min_length=1, max_length=256)
    qc_flag: Literal["flag", "exclude", ""]
    qc_reason: str = Field(default="", max_length=2000)


def project_overview(scope: PathScope) -> dict:
    project = scope.project()
    return {"project_id": project.id, "title": project.config.project.title,
            "species": project.config.project.species.value,
            "experiments": [scope.path(p).name for p in project.experiment_dirs()],
            "pipeline_scripts": [item.name or f"script_{i + 1}"
                                 for i, item in enumerate(project.config.pipeline.scripts)]}


def validation(root: Path, experiments: list[str] | None = None, *, scope: PathScope | None = None):
    from omicsweave.core.validate import validate_project

    if scope:
        project = scope.project()
        targets = experiments or [p.name for p in project.experiment_dirs()]
        for exp_id in targets:
            scope.experiment(project, exp_id)
            scope.tree(project.path("results") / exp_id)
    return validate_project(root, experiments)


def validation_summary(report) -> dict:
    return {"ok": report.ok, "issues": [asdict(issue) for issue in report.issues],
            "error_count": len(report.errors), "warning_count": len(report.warnings)}


def experiment_view(scope: PathScope, experiment: str) -> dict:
    from omicsweave.core.validate import biological_samples

    exp = scope.experiment(scope.project(), experiment)
    samples = biological_samples(exp)
    matrix = samples.groupby(["omics", "tissue"], dropna=False).size().rename("samples").reset_index()
    return {"experiment_id": exp.id, "experiment_config": exp.config.model_dump(mode="json"),
            "subject_count": len(exp.subjects), "sample_count": len(exp.samples),
            "sample_matrix": matrix.to_dict("records")}


def inspect_inputs(scope: PathScope, experiment: str | None = None, paths: list[str] | None = None) -> dict:
    if not experiment and not paths:
        raise ServiceError(ErrorCode.CLI_USAGE_ERROR, "请指定实验 ID 或文件路径")
    inputs = [(None, scope.path(value, external=True)) for value in paths or []]
    if experiment:
        project = scope.project()
        exp = scope.experiment(project, experiment)
        for layer in exp.config.omics.enabled():
            spec = exp.config.omics.settings(layer).input
            values = list(dict.fromkeys([*spec.files, *(getattr(spec, "tissue_files", None) or {}).values()]))
            inputs.extend((layer.value, scope.path(value, external=True)) for value in values)
    rows = []
    for layer, path in inputs[:200]:
        # 只报告文件信息，避免在发现时读取大矩阵或返回实验记录。
        rows.append({"path": str(path), "omics": layer, "exists": path.exists(),
                     "kind": "directory" if path.is_dir() else "file",
                     "size_bytes": path.stat().st_size if path.is_file() else None})
    return {"experiment_id": experiment, "files": rows, "total_files": len(inputs), "truncated": len(inputs) > len(rows)}


def _merge(document: dict, patch: dict) -> dict:
    result = dict(document)
    for key, value in patch.items():
        if value is None:
            result.pop(key, None)
        elif isinstance(value, dict):
            before = result.get(key)
            result[key] = _merge(before if isinstance(before, dict) else {}, value)
        else:
            result[key] = value
    return result


def _apply(scope: PathScope, exp, path: Path, old: bytes, new: bytes, options: EditOptions,
           panel_guard: Callable[[], list[str]] | None = None) -> dict:
    from omicsweave.core.jobs import ACTIVE, list_status

    digest = hashlib.sha256(old).hexdigest()
    if options.dry_run or new == old:
        return {"applied": False, "dry_run": options.dry_run, "sha256": digest}
    if not options.reason.strip() or options.expected_sha256 != digest:
        raise ServiceError(ErrorCode.JOB_CONFLICT, "提交需要非空理由及当前预览返回的 expected_sha256")
    with file_lock(scope.root / ".omicsweave-config.lock"):
        if scope.path(path) != path or path.read_bytes() != old:
            raise ServiceError(ErrorCode.JOB_CONFLICT, "目标文件已变更，请重新预览")
        root = scope.path(exp.project.path("results"))
        if any(job["status"] in ACTIVE and (not job["experiments"] or exp.id in job["experiments"])
               for job in list_status(root)):
            raise ServiceError(ErrorCode.JOB_CONFLICT, "该实验存在活动作业")
        if panel_guard is not None and (notices := panel_guard()):
            raise ServiceError(ErrorCode.VALIDATION_FAILED, "面板来源需要重新确认", details={"notices": notices})
        replace_bytes(path, new)
    return {"applied": True, "dry_run": False, "sha256": digest,
            "written_sha256": hashlib.sha256(new).hexdigest(), "reason": options.reason}


def _impact(stage: str | None) -> dict:
    return {"invalidates_omics_bundles": stage == "run", "invalidates_integration_bundles": stage is not None,
            "recommended_rerun_from": stage}


def _panel_notices(exp, selection: dict[str, str], declared: bool) -> list[str]:
    from omicsweave.contract.manifest import read_manifest
    from omicsweave.core.overrides import TRIALS_DIR

    results = exp.project.path("results").resolve()
    paths = [*(results / exp.id).glob("*/manifest.json"),
             *(results / TRIALS_DIR).glob(f"trial_*/{exp.id}/*/manifest.json")]
    published = [path for path in paths if path.resolve() == path and path.is_file()]
    if not published:
        return []
    if not declared or selection["mode"] != "post_hoc":
        return ["已有正式或试验结果。面板变更需要显式登记 selection.mode=post_hoc"]
    requested = {word.rstrip(".") for word in re.findall(r"[A-Za-z0-9_.-]+", selection["basis"])}
    for path in published:
        if path.parent.name not in requested:
            continue
        try:
            manifest = read_manifest(path.parent)
        except (OSError, ValueError):
            continue
        if (manifest.project_id, manifest.experiment_id, manifest.run_id) == (exp.project.id, exp.id, path.parent.name):
            return []
    return ["selection.basis 需要引用本实验已有且清单有效的运行 ID，可同时附上解释"]


def patch_config(scope: PathScope, experiment: str, patch: dict, options: EditOptions,
                 target: Literal["experiment_yaml", "panels_yaml"] = "experiment_yaml") -> dict:
    from omicsweave.core.validate import ValidationReport, validate_experiment
    from omicsweave.integration.dossier.panels import PanelError, load_panel_selection, load_panels

    project = scope.project()
    exp = scope.experiment(project, experiment)
    scope.tree(project.path("results") / exp.id)
    if target == "panels_yaml":
        configured = exp.config.integration.parameters.get("reaction_dossier", {}).get("panels")
        if not configured:
            raise ServiceError(ErrorCode.INPUT_NOT_FOUND, "实验未配置面板文件")
        path = scope.path(configured, base=exp.root)
    else:
        path = scope.path(exp.root / EXPERIMENT_FILE)
    old = path.read_bytes()
    text = old.decode("utf-8-sig")
    source = yaml.safe_load(text)
    if not isinstance(source, dict) or not patch:
        raise ServiceError(ErrorCode.CLI_USAGE_ERROR, "补丁和原始 YAML 须为非空映射")
    updated = _merge(source, patch)
    changes = [key for key in sorted(source.keys() | updated.keys()) if source.get(key) != updated.get(key)]
    panel_guard = None
    if target == "panels_yaml":
        try:
            panels = load_panels(path, document=updated)
            for panel_id in panels:
                file_component(panel_id)
            selection = load_panel_selection(path, document=updated)
        except PanelError as exc:
            raise ServiceError(ErrorCode.CONFIG_SCHEMA_ERROR, str(exc)) from exc
        if changes:
            declared = (patch.get("selection") or {}).get("mode") == "post_hoc"
            panel_guard = lambda: _panel_notices(exp, selection, declared)
        candidate = exp
    else:
        candidate = scope.experiment(project, experiment, document=updated)
        if candidate.id != exp.id:
            raise ServiceError(ErrorCode.CLI_USAGE_ERROR, "补丁不能更改实验 ID")
        old_panels = exp.config.integration.parameters.get("reaction_dossier", {}).get("panels")
        new_panels = candidate.config.integration.parameters.get("reaction_dossier", {}).get("panels")
        if new_panels and new_panels != old_panels:
            panel_path = scope.path(new_panels, base=candidate.root, external=True)

            def panel_guard() -> list[str]:
                return _panel_notices(exp, load_panel_selection(panel_path), True)
    notices = panel_guard() if panel_guard else []
    issues = ValidationReport(validate_experiment(candidate))
    if not options.dry_run and (not issues.ok or notices):
        raise ServiceError(ErrorCode.VALIDATION_FAILED, "配置预检未通过",
                           details={"issues": [asdict(issue) for issue in issues.issues], "notices": notices})
    rendered = yaml.safe_dump(updated, sort_keys=False, allow_unicode=True) if changes else text
    yaml_diff = "".join(difflib.unified_diff(text.splitlines(keepends=True), rendered.splitlines(keepends=True),
                                            fromfile=path.name, tofile=path.name))
    new = rendered.encode("utf-8-sig" if old.startswith(b"\xef\xbb\xbf") else "utf-8")
    stage = ("integrate" if target == "panels_yaml" or changes == ["integration"] else "run") if changes else None
    return {**_apply(scope, exp, path, old, new, options, panel_guard), "path": str(path), "yaml_diff": yaml_diff,
            "validation": validation_summary(issues), "guardrail_notices": notices,
            "invalidation_impact": _impact(stage)}


def update_sample_qc(scope: PathScope, experiment: str, updates: list[QCUpdate], options: EditOptions) -> dict:
    import pandas as pd

    from omicsweave.core.config import apply_sample_tissue_aliases
    from omicsweave.core.reliability import contrast_reliability, sample_flags
    from omicsweave.core.validate import ValidationReport, validate_experiment

    project = scope.project()
    exp = scope.experiment(project, experiment)
    path = scope.path(exp.config.design.samples, base=exp.root)
    scope.tree(project.path("results") / exp.id)
    old = path.read_bytes()
    try:
        text = old.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise ServiceError(ErrorCode.CONFIG_SCHEMA_ERROR, "QC 写入需要 UTF-8 样本表") from exc
    comments = [line for line in text.splitlines(keepends=True) if line.lstrip().startswith("#")]
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    reader = csv.DictReader(io.StringIO("".join(line for line in text.splitlines(keepends=True)
                                             if line.strip() and not line.lstrip().startswith("#"))), delimiter=sep)
    columns = list(reader.fieldnames or [])
    rows = list(reader)
    if len(set(columns)) != len(columns) or any(None in row or None in row.values() for row in rows):
        raise ServiceError(ErrorCode.CONFIG_SCHEMA_ERROR, "样本表列名重复或行宽与表头不一致")
    by_id = {row.get("sample_id"): row for row in rows}
    requested = {item.sample_id: item for item in updates}
    if not 1 <= len(updates) <= 200 or len(requested) != len(updates):
        raise ServiceError(ErrorCode.CLI_USAGE_ERROR, "提供 1–200 个互不重复的样本更新")
    missing = sorted(set(requested) - set(by_id))
    if missing or len(by_id) != len(rows) or "sample_id" not in columns:
        raise ServiceError(ErrorCode.VALIDATION_FAILED, "样本表 ID 缺失、重复或无法匹配", details={"missing": missing})
    changed = []
    for sample, item in requested.items():
        if item.qc_flag and not item.qc_reason.strip():
            raise ServiceError(ErrorCode.CLI_USAGE_ERROR, f"标记样本 {sample} 需要给出依据")
        before = {key: by_id[sample].get(key, "") for key in ("qc_flag", "qc_reason")}
        after = {"qc_flag": item.qc_flag, "qc_reason": item.qc_reason.strip() if item.qc_flag else ""}
        if before != after:
            by_id[sample].update(after)
            changed.append({"sample_id": sample, "before": before, "after": after})
    columns.extend(key for key in ("qc_flag", "qc_reason") if key not in columns)
    frame = pd.DataFrame(rows, columns=columns).fillna("")
    samples, aliases = apply_sample_tissue_aliases(frame, exp.config)
    candidate = replace(exp, samples=samples, raw_tissues=aliases)
    issues = ValidationReport(validate_experiment(candidate))
    if not options.dry_run and not issues.ok:
        raise ServiceError(ErrorCode.VALIDATION_FAILED, "更新后的样本表未通过校验",
                           details={"issues": [asdict(issue) for issue in issues.issues]})

    keys = samples[["omics", "tissue"]].drop_duplicates()
    keys = keys[keys["omics"].isin([layer.value for layer in exp.config.omics.enabled()])]
    keys = keys.merge(pd.DataFrame({"contrast_id": [c.id for c in exp.config.design.contrasts]}), how="cross")
    projected = []
    if changed and not keys.empty:
        columns_key = ["omics", "tissue", "contrast_id"]
        before = contrast_reliability(exp, keys, sample_flags(exp))
        after = contrast_reliability(candidate, keys, sample_flags(candidate))
        delta = before[columns_key + ["status"]].merge(after[columns_key + ["status"]],
                                                        on=columns_key, suffixes=("_before", "_after"))
        projected = delta[delta.status_before != delta.status_after].to_dict("records")
    output = io.StringIO()
    output.writelines(comments)
    writer = csv.DictWriter(output, columns, delimiter=sep, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    new = output.getvalue().encode("utf-8-sig" if old.startswith(b"\xef\xbb\xbf") else "utf-8") if changed else old
    return {**_apply(scope, exp, path, old, new, options), "changed_rows": changed,
            "projected_reliability_changes": projected, "validation": validation_summary(issues),
            "recommended_rerun_from": "run" if changed else None}
