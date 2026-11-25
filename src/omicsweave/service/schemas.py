"""版本化的机器响应与应用错误。"""
from __future__ import annotations

import math
from datetime import datetime, timezone
from enum import Enum
from typing import Any, Literal

from pydantic import BaseModel, Field, ValidationError


class ErrorCode(str, Enum):
    CLI_USAGE_ERROR = "CLI_USAGE_ERROR"
    CONFIG_SCHEMA_ERROR = "CONFIG_SCHEMA_ERROR"
    VALIDATION_FAILED = "VALIDATION_FAILED"
    INPUT_NOT_FOUND = "INPUT_NOT_FOUND"
    STALE_SOURCE_BUNDLE = "STALE_SOURCE_BUNDLE"
    OPTIONAL_DEP_MISSING = "OPTIONAL_DEP_MISSING"
    SECURITY_POLICY_VIOLATION = "SECURITY_POLICY_VIOLATION"
    EXECUTION_FAILED = "EXECUTION_FAILED"
    CONTRACT_ERROR = "CONTRACT_ERROR"
    JOB_NOT_FOUND = "JOB_NOT_FOUND"
    JOB_CONFLICT = "JOB_CONFLICT"
    JOB_CANCELLED = "JOB_CANCELLED"


EXIT_CODES = {
    ErrorCode.CLI_USAGE_ERROR: 2,
    ErrorCode.CONFIG_SCHEMA_ERROR: 10,
    ErrorCode.VALIDATION_FAILED: 11,
    ErrorCode.INPUT_NOT_FOUND: 12,
    ErrorCode.STALE_SOURCE_BUNDLE: 13,
    ErrorCode.OPTIONAL_DEP_MISSING: 14,
    ErrorCode.SECURITY_POLICY_VIOLATION: 15,
    ErrorCode.EXECUTION_FAILED: 20,
    ErrorCode.CONTRACT_ERROR: 21,
    ErrorCode.JOB_NOT_FOUND: 22,
    ErrorCode.JOB_CONFLICT: 22,
    ErrorCode.JOB_CANCELLED: 30,
}


class ErrorDetail(BaseModel):
    code: ErrorCode
    message: str
    hint: str | None = None
    retryable: bool = False
    details: dict[str, Any] = Field(default_factory=dict)


class ServiceError(Exception):
    def __init__(self, code: ErrorCode, message: str, *, hint: str | None = None,
                 retryable: bool = False, details: dict | None = None):
        super().__init__(message)
        self.detail = ErrorDetail(code=code, message=message, hint=hint,
                                  retryable=retryable, details=details or {})

    @property
    def exit_code(self) -> int:
        return EXIT_CODES[self.detail.code]


class ServiceWarning(BaseModel):
    code: str
    where: str
    message: str


class ServiceResponse(BaseModel):
    ok: bool
    schema_version: Literal[1] = 1
    command: str | None
    project_dir: str | None
    timestamp: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))
    data: Any = None
    warnings: list[ServiceWarning] = Field(default_factory=list)
    error: ErrorDetail | None = None


def json_safe(value):
    """保留数值缺失，递归处理 Parquet 和统计结果中的非有限浮点数。"""
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {k: json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    return value


def classify_error(exc: Exception) -> ServiceError:
    if isinstance(exc, ServiceError):
        return exc
    if isinstance(exc, ValidationError):
        return ServiceError(ErrorCode.CLI_USAGE_ERROR, "请求参数校验失败",
                            details={"field_errors": exc.errors(include_url=False, include_context=False, include_input=False)})

    from omicsweave.contract.dataset import ContractError
    from omicsweave.core.config import ConfigError
    from omicsweave.core.overrides import OverrideError
    from omicsweave.core.path_scope import ScopeError
    from omicsweave.core.progress import JobCancelledError
    from omicsweave.core.run import RunError

    cause: BaseException | None = exc
    while cause is not None:
        if isinstance(cause, ScopeError):
            return ServiceError(ErrorCode.SECURITY_POLICY_VIOLATION, str(cause))
        if isinstance(cause, OverrideError):
            return ServiceError(ErrorCode.CLI_USAGE_ERROR, str(cause))
        if isinstance(cause, ConfigError):
            return ServiceError(ErrorCode(cause.code), str(exc), details={"field_errors": cause.messages})
        if isinstance(cause, RunError) and cause.code:
            return ServiceError(ErrorCode(cause.code), str(exc), details=cause.details,
                                retryable=cause.code == "STALE_SOURCE_BUNDLE")
        if isinstance(cause, JobCancelledError):
            return ServiceError(ErrorCode.JOB_CANCELLED, str(exc))
        if isinstance(cause, FileNotFoundError):
            return ServiceError(ErrorCode.INPUT_NOT_FOUND, str(exc))
        if isinstance(cause, ModuleNotFoundError) and cause.name and not cause.name.startswith("omicsweave"):
            return ServiceError(ErrorCode.OPTIONAL_DEP_MISSING, str(exc), details={"module": cause.name})
        if isinstance(cause, ContractError):
            return ServiceError(ErrorCode.CONTRACT_ERROR, str(exc))
        cause = cause.__cause__
    return ServiceError(ErrorCode.EXECUTION_FAILED, str(exc) or type(exc).__name__)
