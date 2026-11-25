"""运行器共享的进度和协作取消接口。"""
from __future__ import annotations

from contextlib import contextmanager
from typing import Protocol


class JobCancelledError(Exception):
    """Worker 已到达可以停止的位置。"""


class ProgressReporter(Protocol):
    job_id: str | None

    def emit(self, event: str, *, stage: str, detail: dict | None = None) -> None: ...

    def check_cancelled(self) -> None: ...


class NullReporter:
    job_id = None

    def emit(self, event: str, *, stage: str, detail: dict | None = None) -> None:
        pass

    def check_cancelled(self) -> None:
        pass


NULL_REPORTER = NullReporter()


@contextmanager
def progress_stage(reporter: ProgressReporter, stage: str, **detail):
    reporter.check_cancelled()
    reporter.emit("stage_start", stage=stage, detail=detail)
    yield
    reporter.emit("stage_end", stage=stage, detail=detail)
