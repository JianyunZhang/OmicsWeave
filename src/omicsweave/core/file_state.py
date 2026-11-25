"""本地跨进程文件锁及 JSON 状态写入。"""
from __future__ import annotations

import errno
import json
import os
import time
from contextlib import contextmanager
from pathlib import Path
from uuid import uuid4


@contextmanager
def file_lock(path: Path):
    """使用 OS 锁，进程退出时自动释放。锁文件本身始终保留。"""
    with path.open("a+b") as handle:
        if os.name == "nt":
            import msvcrt

            deadline = time.monotonic() + 30
            while True:
                try:
                    handle.seek(0)
                    msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
                    break
                except OSError as exc:
                    if exc.errno not in (errno.EACCES, errno.EAGAIN) or time.monotonic() >= deadline:
                        raise
                    time.sleep(0.02)
        else:
            import fcntl

            fcntl.flock(handle, fcntl.LOCK_EX)
        yield


def replace_bytes(path: Path, contents: bytes) -> None:
    """以同目录临时文件替换目标，读者只会看到完整文件。"""
    pending = path.with_name(f".{path.name}.{uuid4().hex}.tmp")
    try:
        pending.write_bytes(contents)
        if path.exists():
            pending.chmod(path.stat().st_mode & 0o777)
        for attempt in range(10):
            try:
                pending.replace(path)
                return
            except PermissionError:
                if os.name != "nt" or attempt == 9:
                    raise
                time.sleep(0.01)
    finally:
        pending.unlink(missing_ok=True)


def write_json(path: Path, data: dict) -> None:
    replace_bytes(path, json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False).encode("utf-8"))
