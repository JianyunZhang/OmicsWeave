"""内部作业入口，argv 仅接收已创建的作业目录。"""
from __future__ import annotations

import sys
import traceback
from pathlib import Path

from omicsweave.core.jobs import JobReporter, _read
from omicsweave.core.project import project_context
from omicsweave.service.execution import ExecutionRequest, execute
from omicsweave.service.schemas import ErrorCode, classify_error


def main(directory: Path) -> int:
    with JobReporter(directory) as reporter:
        try:
            spec = _read(directory / "spec.json")
            if spec["schema_version"] != 1:
                raise ValueError("无法读取此版本的作业规格")
            root = Path(spec["project_dir"])
            with project_context(root):
                result = execute(root, ExecutionRequest.model_validate(spec["request"]), reporter)
        except Exception as exc:
            failure = classify_error(exc)
            cancelled = failure.detail.code == ErrorCode.JOB_CANCELLED
            if not cancelled:
                traceback.print_exc()
            reporter.finish("cancelled" if cancelled else "failed", failure.exit_code,
                            error=failure.detail.model_dump(mode="json"))
            return failure.exit_code
        reporter.finish("completed", 0, result=result)
        return 0


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1])))
