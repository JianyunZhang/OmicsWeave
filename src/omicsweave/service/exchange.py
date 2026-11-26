"""选择确切运行并构建可携带的只读结果快照。"""
from __future__ import annotations

from pathlib import Path

from omicsweave.contract.exchange import export_exchange, open_exchange
from omicsweave.contract.manifest import read_manifest
from omicsweave.core.overrides import formal_bundle, trial_root
from omicsweave.core.path_scope import PathScope, identifier

from .schemas import ErrorCode, ServiceError


def build_exchange(project_dir: Path, experiment: str, run_id: str, directory: Path,
                   trial_id: str | None = None) -> dict:
    if run_id == "latest":
        raise ServiceError(ErrorCode.CLI_USAGE_ERROR, "交换需要指定确切的运行 ID")
    scope = PathScope.create(project_dir)
    project = scope.project()
    root = project.path("results")
    if trial_id:
        root = trial_root(root, trial_id)
    source = formal_bundle(root, identifier(experiment), identifier(run_id))
    manifest = read_manifest(source)
    if (manifest.project_id, manifest.experiment_id, manifest.run_id) != (project.id, experiment, run_id):
        raise ServiceError(ErrorCode.CONTRACT_ERROR, "运行清单与所请求的结果包不一致")
    output = scope.path(directory)
    # 保留正式与试验结果根的原有扫描语义。
    if output.is_relative_to(project.path("results").resolve() / experiment):
        raise ServiceError(ErrorCode.SECURITY_POLICY_VIOLATION, "请将交换包写在正式实验结果树之外")
    try:
        exported = export_exchange(source, output)
    except FileExistsError as exc:
        raise ServiceError(ErrorCode.JOB_CONFLICT, str(exc)) from exc
    return {"directory": str(output), "exchange": exported.model_dump(mode="json")}


def verify_exchange(directory: Path) -> dict:
    snapshot = open_exchange(directory)
    return {"directory": str(snapshot.root), "verified_artifacts": len(snapshot.manifest.artifacts),
            "exchange": snapshot.manifest.model_dump(mode="json")}
