"""隔离试验的参数、目录和运行谱系。正式 v0.1 清单仍使用原有字段。"""
from __future__ import annotations

import json
import re
from copy import deepcopy
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Literal
from uuid import uuid4

import yaml
from pydantic import BaseModel, ConfigDict, Field

from omicsweave.contract.manifest import MANIFEST_FILE, read_manifest

from .config import EXPERIMENT_FILE, Experiment, Project, load_experiment
from .file_state import write_json
from .path_scope import identifier
from .provenance import sha256_file

TRIALS_DIR = "_trials"
TRIAL_RECORD = "trial.json"
LINEAGE_PATH = Path("metadata", "lineage.json")
TRIAL_ID = re.compile(r"^trial_[a-f0-9]{32}\Z")
Intent = Literal["sensitivity_check", "post_hoc_exploration", "qc_refinement", "visualization_tune", "baseline", "manual"]


class OverrideError(ValueError):
    """试验参数或基线引用不符合约定。"""


class TrialOptions(BaseModel):
    model_config = ConfigDict(extra="forbid")
    overrides: dict[str, Any] = Field(default_factory=dict)
    intent: Intent = "sensitivity_check"
    rationale: str = Field(min_length=1, max_length=2000)
    parent_run_id: str | None = None


class RunLineage(BaseModel):
    model_config = ConfigDict(extra="forbid")
    schema_version: Literal[1] = 1
    trial_id: str = Field(pattern=TRIAL_ID.pattern.replace(r"\Z", "$"))
    project_id: str
    experiment_id: str
    intent: Intent
    rationale: str
    job_id: str | None = None
    parent_run_id: str
    parent_manifest_sha256: str
    source_run_id: str | None = None
    source_manifest_sha256: str | None = None
    base_config_sha256: str
    project_alpha: float
    overrides: dict[str, Any]


def _add_override(results: dict, parts: tuple[str, ...], value: Any) -> None:
    if not 1 <= len(parts) <= 32 or any(not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_-]*", key) for key in parts):
        raise OverrideError("覆盖路径使用点分隔的配置键")
    if isinstance(value, dict) and value:
        for key, child in value.items():
            if not isinstance(key, str):
                raise OverrideError("覆盖的映射键须为字符串")
            _add_override(results, (*parts, *key.split(".")), child)
        return
    name = ".".join(parts)
    if any(name == other or name.startswith(f"{other}.") or other.startswith(f"{name}.") for other in results):
        raise OverrideError(f"覆盖路径重复或重叠 {name}")
    try:
        # 拒绝 YAML 的日期、集合和非有限浮点数等不可交换值。
        results[name] = json.loads(json.dumps(value, allow_nan=False))
    except (ValueError, TypeError) as exc:
        raise OverrideError(f"覆盖值须可序列化为 JSON {name}") from exc


def flatten_overrides(document: dict) -> dict:
    if not isinstance(document, dict):
        raise OverrideError("覆盖文件与覆盖参数须为映射")
    result: dict = {}
    for key, value in document.items():
        if not isinstance(key, str):
            raise OverrideError("覆盖路径须为字符串")
        _add_override(result, tuple(key.split(".")), value)
    return result


def parse_overrides(items: list[str] | None, filename: Path | None) -> dict:
    class UniqueKeysLoader(yaml.SafeLoader):
        def construct_mapping(self, node, deep=False):
            mapping = super().construct_mapping(node, deep=deep)
            if len(mapping) != len(node.value):
                raise OverrideError("覆盖 YAML 包含重复键")
            return mapping

    try:
        result = flatten_overrides(yaml.load(filename.read_text(encoding="utf-8-sig"), Loader=UniqueKeysLoader)) if filename else {}
        for item in items or []:
            key, separator, value = item.partition("=")
            if not separator or not value.strip():
                raise OverrideError("--override 使用 key.path=yaml_value")
            _add_override(result, tuple(key.strip().split(".")), yaml.load(value, Loader=UniqueKeysLoader))
    except yaml.YAMLError as exc:
        raise OverrideError(f"无法解析覆盖参数 {exc}") from exc
    return result


def _check_path(parts: list[str]) -> None:
    from .param_schemas import INTEGRATION_SCHEMAS, OMICS_SCHEMAS

    changed = parts[3:]
    if parts[:2] == ["design", "adaptive_shrinkage"]:
        return
    if len(parts) < 4:
        raise OverrideError(f"只能覆盖已登记的分析方法、整合和展示参数 {'.'.join(parts)}")
    omics = parts[0] == "omics" and parts[1] in OMICS_SCHEMAS and parts[2] == "analysis"
    integration = parts[:2] == ["integration", "parameters"] and parts[2] in INTEGRATION_SCHEMAS
    model = OMICS_SCHEMAS[parts[1]] if omics else INTEGRATION_SCHEMAS[parts[2]] if integration else None
    if (model is None or changed[0] not in model.model_fields
            or any(name in {"gene_sets", "sources", "panels", "layout"} for name in changed)):
        raise OverrideError(f"该路径不可用于试验覆盖 {'.'.join(parts)}")


def apply_overrides(document: dict, overrides: dict) -> dict:
    effective = deepcopy(document)
    for name, value in overrides.items():
        parts = name.split(".")
        _check_path(parts)
        current = effective
        for key in parts[:-1]:
            child = current.get(key)
            if not isinstance(child, dict):
                # 接受已有的 submaps: false / abstract: false 简写。
                child = {"enabled": child} if isinstance(child, bool) and key in {"submaps", "abstract"} else {}
                current[key] = child
            current = child
        current[parts[-1]] = value
    return effective


def trial_root(root: Path, trial_id: str) -> Path:
    if not TRIAL_ID.fullmatch(trial_id):
        raise OverrideError("无效的试验 ID")
    directory = root.resolve() / TRIALS_DIR / trial_id
    if directory.resolve() != directory:
        raise OverrideError("试验目录必须位于结果根内，且不能重定向")
    return directory


def formal_bundle(root: Path, experiment: str, run_id: str) -> Path:
    base = root.resolve() / identifier(experiment)
    bundle = base / identifier(run_id)
    if bundle.resolve() != bundle or (bundle / MANIFEST_FILE).resolve() != bundle / MANIFEST_FILE:
        raise OverrideError("试验的来源必须是同一实验下的正式结果包")
    return bundle


def changed_source_outputs(bundle: Path, only: Path | None = None) -> list[str]:
    """核对复用的科学数据和配置快照，图件及交付副本由原索引管理。"""
    changed = []
    for record in read_manifest(bundle).files:
        path = bundle / record.path
        if record.role != "output" or path.suffix.lower() not in {".parquet", ".json"} or (only and path != only):
            continue
        if (not path.resolve().is_relative_to(bundle.resolve()) or not path.is_file()
                or path.stat().st_size != record.size or sha256_file(path) != record.sha256):
            changed.append(record.path)
    return changed


@dataclass(frozen=True)
class PreparedTrial:
    lineage: RunLineage
    formal_root: Path
    root: Path
    experiment: Experiment
    source: Path | None

    def reserve(self) -> None:
        from .jobs import root_lock

        with root_lock(self.formal_root):
            if trial_root(self.formal_root, self.lineage.trial_id) != self.root:
                raise OverrideError("试验输出目录已变化")
            for run_id, digest in ((self.lineage.parent_run_id, self.lineage.parent_manifest_sha256),
                                   (self.lineage.source_run_id, self.lineage.source_manifest_sha256)):
                if run_id and sha256_file(formal_bundle(
                    self.formal_root, self.experiment.id, run_id,
                ) / MANIFEST_FILE) != digest:
                    raise OverrideError("试验的基线清单已变化，请重新预检")
            current = self.experiment.root / EXPERIMENT_FILE
            if sha256_file(current) != self.lineage.base_config_sha256:
                raise OverrideError("实验配置已变化，请重新预检")
            self.root.mkdir(parents=True, exist_ok=False)
            write_json(self.root / TRIAL_RECORD, self.lineage.model_dump(mode="json"))


def prepare_trial(project: Project, experiment: str, options: TrialOptions, *,
                  reused_source: str | None, job_id: str | None = None) -> PreparedTrial:
    from hashlib import sha256

    if not options.rationale.strip():
        raise OverrideError("试验需要填写理由")
    root = project.path("results").resolve()
    parent_id = options.parent_run_id or reused_source
    if not parent_id:
        raise OverrideError("试验需要明确的 --parent-run 或 --source")

    def baseline(run_id: str) -> Path:
        path = formal_bundle(root, experiment, run_id)
        manifest = read_manifest(path)
        if (manifest.project_id, manifest.experiment_id, manifest.run_id) != (project.id, experiment, run_id):
            raise OverrideError("试验的来源清单应与项目、实验和运行 ID 一致")
        return path

    parent = baseline(parent_id)
    source = baseline(reused_source) if reused_source else None
    overrides = flatten_overrides(options.overrides)
    path = project.path("experiments") / identifier(experiment) / EXPERIMENT_FILE
    raw = path.read_bytes()
    document = yaml.safe_load(raw.decode("utf-8-sig"))
    if not isinstance(document, dict):
        raise OverrideError("实验配置须为映射")
    trial_id = f"trial_{uuid4().hex}"
    directory = trial_root(root, trial_id)
    trial_project = replace(project, config=project.config.model_copy(deep=True))
    trial_project.config.project.paths.results = str(directory)
    effective = load_experiment(trial_project, experiment, document=apply_overrides(document, overrides))
    if effective.id != experiment:
        raise OverrideError("试验要求配置实验 ID 与目录一致")
    lineage = RunLineage(
        trial_id=trial_id, project_id=project.id, experiment_id=experiment, intent=options.intent,
        rationale=options.rationale.strip(), job_id=job_id, parent_run_id=parent_id,
        parent_manifest_sha256=sha256_file(parent / MANIFEST_FILE), source_run_id=reused_source,
        source_manifest_sha256=sha256_file(source / MANIFEST_FILE) if source else None,
        base_config_sha256=sha256(raw).hexdigest(), project_alpha=project.config.defaults.alpha,
        overrides=overrides,
    )
    return PreparedTrial(lineage, root, directory, effective, source)


def lineage_for(bundle: Path) -> RunLineage | None:
    path = bundle / LINEAGE_PATH
    if not path.is_file():
        return None
    if not path.resolve().is_relative_to(bundle.resolve()):
        raise OverrideError("谱系文件超出结果目录")
    return RunLineage.model_validate_json(path.read_bytes())


def trial_references(root: Path, experiment: str) -> set[str]:
    """登记过的试验保护其基线及复用包，便于审计和重复分析。"""
    references: set[str] = set()
    for path in (root / TRIALS_DIR).glob(f"trial_*/{TRIAL_RECORD}"):
        try:
            if path.resolve() != path or trial_root(root, path.parent.name) != path.parent:
                continue
            lineage = RunLineage.model_validate_json(path.read_bytes())
        except (ValueError, OSError):
            continue
        if lineage.experiment_id == experiment:
            references.add(lineage.parent_run_id)
            if lineage.source_run_id:
                references.add(lineage.source_run_id)
    return references
