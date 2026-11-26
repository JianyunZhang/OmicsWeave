"""Versioned, read-only result snapshots for file and MCP exchange."""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from hashlib import sha256
from pathlib import Path, PurePosixPath, PureWindowsPath
from typing import Annotated, Literal

from pydantic import AfterValidator, BaseModel, ConfigDict, Field, ValidationError, model_validator

from omicsweave.core.file_state import write_json
from omicsweave.core.model import Omics
from omicsweave.core.overrides import LINEAGE_PATH, RunLineage
from omicsweave.core.provenance import sha256_file

from .dataset import FEATURES_FILE, MATRIX_FILE, META_FILE, ContractError, Dataset, DatasetMeta, read_dataset
from .manifest import MANIFEST_FILE, FileRecord, RunManifest

EXCHANGE_FILE = "exchange.json"
EXCHANGE_FORMAT = "omicsweave.analysis-result.exchange"
EXCHANGE_VERSION = 1
_OMICS_TABLES = {"differential", "derived_differential", "differential_by_tracer", "enrichment"}
_INTEGRATION_TABLES = {
    "flux_anchored": {"flux_evidence_profile", "flux_pathway_profile"},
    "pathway_convergence": {"pathway_fusion"},
    "reaction_dossier": {"reaction_assessment", "axis_assessment"},
    "reversal": {"reversal", "reversal_summary"},
}
_LAYERS = {layer.value for layer in Omics}
_DATASET_FILES = {META_FILE: "dataset_metadata", MATRIX_FILE: "matrix", FEATURES_FILE: "features"}
Kind = Literal["table", "dataset_metadata", "matrix", "features"]
Digest = Annotated[str, Field(pattern=r"^[a-f0-9]{64}$")]


def relative_path(value: str) -> str:
    path = PurePosixPath(value)
    if (not path.parts or path.is_absolute() or path.as_posix() != value
            or any(part == ".." or part.endswith((" ", ".")) or any(c in part for c in '\\:*?"<>|')
                   or any(ord(c) < 32 for c in part) or PureWindowsPath(part).is_reserved()
                   for part in path.parts)):
        raise ValueError("交换路径必须是跨平台的规范相对路径")
    return value


RelativePath = Annotated[str, AfterValidator(relative_path)]


def is_result_table(value: str) -> bool:
    path = PurePosixPath(value)
    parts = path.parts
    if len(parts) != 4 or parts[2] != "tables" or path.suffix != ".parquet":
        return False
    namespace, layer = parts[:2]
    return ((namespace == "omics" and layer in _LAYERS and path.stem in _OMICS_TABLES)
            or (namespace == "integration" and path.stem in _INTEGRATION_TABLES.get(layer, set())))


def _dataset_path(value: str) -> bool:
    parts = PurePosixPath(value).parts
    return len(parts) >= 4 and parts[0] == "omics" and parts[1] in _LAYERS and parts[2] in {"data", "derived"}


def _file(root: Path, relative: str) -> Path:
    path = root / relative_path(relative)
    if path.resolve() != path or not path.is_file():
        raise ContractError(f"交换文件不存在或经过重定向 {relative}")
    return path


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid")


class ExchangeRun(_Strict):
    project_id: str
    experiment_id: str
    run_id: str
    created_at: datetime
    omicsweave_version: str
    bundle_schema_version: str
    manifest_sha256: Digest
    config_sha256: Digest | None = None


class ExchangeLineage(_Strict):
    trial_id: str
    intent: str
    parent_run_id: str
    parent_manifest_sha256: Digest
    source_run_id: str | None = None
    source_manifest_sha256: Digest | None = None
    lineage_sha256: Digest


class ExchangeArtifact(_Strict):
    path: RelativePath
    kind: Kind
    size: int = Field(ge=0)
    sha256: Digest


class ExchangeDataset(_Strict):
    directory: RelativePath
    metadata: DatasetMeta


class AnalysisResultExchange(_Strict):
    # 交换配置文件和 RunManifest 分开版本化。
    format: Literal["omicsweave.analysis-result.exchange"]
    version: Literal[1]
    run: ExchangeRun
    lineage: ExchangeLineage | None = None
    datasets: list[ExchangeDataset]
    artifacts: list[ExchangeArtifact] = Field(min_length=1)

    @model_validator(mode="after")
    def _inventory(self) -> AnalysisResultExchange:
        actual = {item.path: item.kind for item in self.artifacts}
        if len({path.casefold() for path in actual}) != len(self.artifacts):
            raise ValueError("交换文件路径重复")
        expected = {path: kind for path, kind in actual.items() if kind == "table" and is_result_table(path)}
        directories = {item.directory.casefold() for item in self.datasets}
        if len(directories) != len(self.datasets):
            raise ValueError("数据集目录重复")
        for item in self.datasets:
            if (not _dataset_path(item.directory) or item.metadata.source
                    or PurePosixPath(item.directory).parts[1] != item.metadata.omics.value):
                raise ValueError("数据集目录或公开元信息不符合交换契约")
            expected.update({f"{item.directory}/{filename}": kind for filename, kind in _DATASET_FILES.items()})
        if actual != expected:
            raise ValueError("交换清单与实际声明的表和数据集不匹配")
        return self


def _source_bytes(root: Path, record: FileRecord) -> bytes:
    contents = _file(root, record.path).read_bytes()
    if len(contents) != record.size or sha256(contents).hexdigest() != record.sha256:
        raise ContractError(f"源结果校验失败 {record.path}")
    return contents


def _copy_artifact(source: Path, target: Path, record: ExchangeArtifact) -> None:
    digest = sha256()
    size = 0
    destination = target / record.path
    destination.parent.mkdir(parents=True, exist_ok=True)
    with _file(source, record.path).open("rb") as reader, destination.open("xb") as writer:
        while block := reader.read(1024 * 1024):
            writer.write(block)
            digest.update(block)
            size += len(block)
    if size != record.size or digest.hexdigest() != record.sha256:
        raise ContractError(f"源结果校验失败 {record.path}")


def export_exchange(source: Path, output: Path) -> AnalysisResultExchange:
    """发布一个固定结果包的受限快照，清单在所有文件复制并检验后写入。"""
    source, output = source.resolve(strict=True), output.resolve()
    if output.is_relative_to(source) or output.exists():
        raise FileExistsError("交换输出须为源结果包之外的全新目录")
    raw = _file(source, MANIFEST_FILE).read_bytes()
    manifest = RunManifest.model_validate_json(raw)
    records = {}
    for item in manifest.files:
        if item.role != "output":
            continue
        name = PurePosixPath(item.path).name
        if not (is_result_table(item.path) or (_dataset_path(item.path) and name in _DATASET_FILES)
                or item.path == LINEAGE_PATH.as_posix()):
            continue
        path = relative_path(item.path)
        if path in records:
            raise ContractError(f"源结果清单包含重复文件 {path}")
        records[path] = item

    datasets = []
    selected = {path: "table" for path in records if is_result_table(path)}
    for path in sorted(records):
        rel = PurePosixPath(path)
        if rel.name != META_FILE or not _dataset_path(rel.parent.as_posix()):
            continue
        meta = DatasetMeta.model_validate_json(_source_bytes(source, records[path]))
        for filename, kind in _DATASET_FILES.items():
            artifact = (rel.parent / filename).as_posix()
            if artifact not in records:
                raise ContractError(f"源结果数据集缺少已登记的产物 {artifact}")
            selected[artifact] = kind
        datasets.append(ExchangeDataset(directory=rel.parent.as_posix(), metadata=meta.model_copy(update={"source": {}})))
    if not selected:
        raise ContractError("结果包没有可交换的标准结果表或数据集")

    lineage = None
    if (source / LINEAGE_PATH).exists():
        if LINEAGE_PATH.as_posix() not in records:
            raise ContractError("源谱系需要包含在结果清单中")
        record = records[LINEAGE_PATH.as_posix()]
        parsed = RunLineage.model_validate_json(_source_bytes(source, record))
        if (parsed.project_id, parsed.experiment_id) != (manifest.project_id, manifest.experiment_id):
            raise ContractError("试验谱系与源结果包身份不一致")
        lineage = ExchangeLineage(trial_id=parsed.trial_id, intent=parsed.intent, parent_run_id=parsed.parent_run_id,
                                  parent_manifest_sha256=parsed.parent_manifest_sha256,
                                  source_run_id=parsed.source_run_id, source_manifest_sha256=parsed.source_manifest_sha256,
                                  lineage_sha256=record.sha256)
    report = AnalysisResultExchange(
        format=EXCHANGE_FORMAT, version=EXCHANGE_VERSION,
        run=ExchangeRun(project_id=manifest.project_id, experiment_id=manifest.experiment_id,
                        run_id=manifest.run_id, created_at=manifest.created_at,
                        omicsweave_version=manifest.omicsweave_version, bundle_schema_version=manifest.bundle_schema_version,
                        manifest_sha256=sha256(raw).hexdigest(), config_sha256=manifest.config_sha256),
        lineage=lineage, datasets=datasets,
        artifacts=[ExchangeArtifact(path=path, kind=kind, size=records[path].size, sha256=records[path].sha256)
                   for path, kind in sorted(selected.items())],
    )
    output.mkdir(parents=True, exist_ok=False)
    for item in report.artifacts:
        _copy_artifact(source, output, item)
    if sha256_file(_file(source, MANIFEST_FILE)) != report.run.manifest_sha256:
        raise ContractError("导出期间源清单发生变更")
    write_json(output / EXCHANGE_FILE, report.model_dump(mode="json"))
    return report


@dataclass(frozen=True)
class ExchangeSnapshot:
    root: Path
    manifest: AnalysisResultExchange

    def artifact_path(self, artifact: ExchangeArtifact) -> Path:
        path = _file(self.root, artifact.path)
        if path.stat().st_size != artifact.size or sha256_file(path) != artifact.sha256:
            raise ContractError(f"交换文件校验失败 {artifact.path}")
        return path

    def dataset(self, directory: str) -> Dataset:
        description = next((item for item in self.manifest.datasets if item.directory == directory), None)
        if description is None:
            raise ContractError("交换包未登记该数据集")
        for item in self.manifest.artifacts:
            if item.path in {f"{directory}/{name}" for name in _DATASET_FILES}:
                self.artifact_path(item)
        result = read_dataset(self.root / directory)
        if result.meta.model_copy(update={"source": {}}) != description.metadata:
            raise ContractError("数据集语义与交换清单不一致")
        return result


def open_exchange(directory: Path) -> ExchangeSnapshot:
    """校验交换清单、登记文件和数据集元信息，可用于独立消费者。"""
    root = directory.resolve(strict=True)
    try:
        manifest = AnalysisResultExchange.model_validate_json(_file(root, EXCHANGE_FILE).read_bytes())
    except ValidationError as exc:
        raise ContractError("交换清单不符合当前版本的 Schema") from exc
    snapshot = ExchangeSnapshot(root, manifest)
    for item in manifest.artifacts:
        snapshot.artifact_path(item)
    for dataset in manifest.datasets:
        meta = DatasetMeta.model_validate_json(_file(root, f"{dataset.directory}/{META_FILE}").read_bytes())
        if meta.model_copy(update={"source": {}}) != dataset.metadata:
            raise ContractError(f"数据集语义与交换清单不一致 {dataset.directory}")
    return snapshot
