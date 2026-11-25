"""固定 MCP 服务的授权范围。"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from omicsweave.core.path_scope import PathScope, identifier
from omicsweave.core.overrides import trial_root
from omicsweave.service.schemas import ErrorCode, ServiceError


@dataclass(frozen=True)
class ServerPolicy:
    scope: PathScope
    allow_writes: bool = False
    allow_custom_scripts: bool = False

    def check_edit(self, dry_run: bool) -> None:
        if not dry_run and not self.allow_writes:
            raise ServiceError(ErrorCode.SECURITY_POLICY_VIOLATION, "配置写入需要服务器启动时启用 --allow-writes")

    def results(self, experiment: str | None = None, trial: str | None = None) -> Path:
        directory = self.scope.path(self.scope.project().path("results"))
        if experiment is not None:
            self.scope.tree(directory / identifier(experiment))
        if trial is not None:
            directory = self.scope.path(trial_root(directory, trial))
            if experiment is not None:
                self.scope.tree(directory / identifier(experiment))
        return directory

    def job_root(self) -> Path:
        root = self.results()
        self.scope.tree(root / "_jobs")
        return root
