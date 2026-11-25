"""纯文件系统后台作业，状态由 Worker 单写，查询只生成实时视图。"""
from __future__ import annotations

import json
import os
import re
import socket
import subprocess
import sys
import threading
from collections import deque
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

from .file_state import file_lock, write_json
from .progress import JobCancelledError

ACTIVE = frozenset({"queued", "running", "unresponsive"})
JOB_ID = re.compile(r"job_[0-9a-f]{32}\Z")
HEARTBEAT_SECONDS = 5
STARTUP_GRACE_SECONDS = 60


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _read(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"{path} 应为 JSON 对象")
    return value


def root_lock(root: Path):
    return file_lock(root / ".omicsweave-jobs.lock")


def job_dir(root: Path, job_id: str) -> Path:
    if not JOB_ID.fullmatch(job_id):
        raise ValueError("无效的 job_id")
    base = (root.resolve() / "_jobs")
    path = base / job_id
    if path.resolve() != path:
        raise ValueError("作业目录不能通过符号链接重定向")
    if not (path / "state.json").is_file():
        raise FileNotFoundError(f"未找到作业 {job_id}")
    return path


def _windows_process(pid: int) -> tuple[bool | None, str | None]:
    import ctypes
    from ctypes import wintypes as win

    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel.OpenProcess.argtypes = (win.DWORD, win.BOOL, win.DWORD)
    kernel.OpenProcess.restype = win.HANDLE
    kernel.CloseHandle.argtypes = (win.HANDLE,)
    kernel.WaitForSingleObject.argtypes = (win.HANDLE, win.DWORD)
    kernel.WaitForSingleObject.restype = win.DWORD
    kernel.GetProcessTimes.argtypes = (win.HANDLE, *([ctypes.POINTER(win.FILETIME)] * 4))
    handle = kernel.OpenProcess(0x1000 | 0x00100000, False, pid)  # 查询身份与等待句柄
    if not handle:
        return (False if ctypes.get_last_error() == 87 else None), None
    try:
        status = kernel.WaitForSingleObject(handle, 0)
        if status == 0:
            return False, None
        times = [win.FILETIME() for _ in range(4)]
        if status != 258 or not kernel.GetProcessTimes(handle, *(ctypes.byref(t) for t in times)):
            return None, None
        born = times[0]
        return True, f"{born.dwHighDateTime}:{born.dwLowDateTime}"
    finally:
        kernel.CloseHandle(handle)


def _process(pid: int) -> tuple[bool | None, str | None]:
    if pid <= 0:
        return False, None
    if os.name == "nt":
        return _windows_process(pid)
    if sys.platform.startswith("linux"):
        try:
            fields = Path(f"/proc/{pid}/stat").read_text().rsplit(")", 1)[1].split()
            boot = Path("/proc/sys/kernel/random/boot_id").read_text().strip()
            return fields[0] not in {"Z", "X"}, f"{boot}:{fields[19]}"
        except (OSError, IndexError):
            pass
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False, None
    except PermissionError:
        return None, None
    return True, None


def submit(root: Path, project: Path, request: dict, experiments: list[str]) -> dict:
    root = root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    with root_lock(root):
        directory = root / "_jobs" / f"job_{uuid4().hex}"
        if directory.resolve() != directory:
            raise ValueError("作业目录不能通过符号链接重定向")
        directory.mkdir(parents=True)
        state = {
            "schema_version": 1, "job_id": directory.name, "kind": request["kind"],
            "project_dir": str(project), "experiments": experiments, "host": socket.gethostname(),
            "status": "queued", "pid": None, "created_at": now_iso(), "started_at": None,
            "heartbeat_at": None, "finished_at": None, "exit_code": None, "progress": None,
            "last_event_seq": 0, "artifacts": {}, "error": None, "result": None,
        }
        write_json(directory / "spec.json", {"schema_version": 1, "project_dir": str(project), "request": request})
        write_json(directory / "state.json", state)
        env = {**os.environ, "PYTHONIOENCODING": "utf-8", "PYTHONUNBUFFERED": "1"}
        env["PYTHONPATH"] = os.pathsep.join(filter(None, [str(Path(__file__).resolve().parents[2]), env.get("PYTHONPATH")]))
        options = (
            {"creationflags": subprocess.DETACHED_PROCESS | subprocess.CREATE_NEW_PROCESS_GROUP}
            if os.name == "nt" else {"start_new_session": True}
        )
        try:
            with (directory / "output.log").open("ab") as output:
                proc = subprocess.Popen([sys.executable, "-m", "omicsweave._worker", str(directory)],
                                        env=env, stdin=subprocess.DEVNULL, stdout=output,
                                        stderr=subprocess.STDOUT, **options)
        except OSError as exc:
            state.update(status="failed", finished_at=now_iso(), exit_code=20,
                         error={"code": "EXECUTION_FAILED", "message": str(exc)})
            write_json(directory / "state.json", state)
            raise
        threading.Thread(target=proc.wait, daemon=True).start()
        _, identity = _process(proc.pid)
        write_json(directory / "launch.json", {"pid": proc.pid, "process_identity": identity})
        return {**state, "pid": proc.pid, "state_path": str(directory / "state.json")}


def read_status(directory: Path) -> dict:
    state = _read(directory / "state.json")
    state["state_path"] = str(directory / "state.json")
    state["cancel_requested"] = (directory / "cancel.flag").exists()
    if state["status"] not in ACTIVE or state["host"] != socket.gethostname():
        return state

    if state["pid"] is None:
        try:
            state.update(_read(directory / "launch.json"))
        except FileNotFoundError:
            pass
    age = (datetime.now(timezone.utc) - datetime.fromisoformat(state["heartbeat_at"] or state["created_at"])).total_seconds()
    state["heartbeat_age_seconds"] = max(0.0, age)
    if state["pid"] is None:
        if age > STARTUP_GRACE_SECONDS:
            state["status"] = "stalled"
        return state

    live, current = _process(state["pid"])
    identity = state.get("process_identity")
    if live is False or (identity and current and current != identity):
        # 进程可能在读到运行态后刚落盘终态并退出，重新读取一次。
        final = _read(directory / "state.json")
        if final["status"] not in ACTIVE:
            return {**state, **final}
        state["status"] = "stalled"
    elif age > (STARTUP_GRACE_SECONDS if state["status"] == "queued" else 3 * HEARTBEAT_SECONDS):
        state["status"] = "unresponsive"
    return state


def list_status(root: Path):
    for directory in (root / "_jobs").glob("job_*"):
        try:
            yield read_status(job_dir(root, directory.name))
        except (OSError, ValueError, KeyError):
            continue


def request_cancel(directory: Path) -> dict:
    state = read_status(directory)
    if state["status"] in ACTIVE:
        (directory / "cancel.flag").touch(exist_ok=True)
        state["cancel_requested"] = True
    return state


def read_logs(directory: Path, tail: int, since_event: int, events_only: bool) -> dict:
    payload = {"job_id": directory.name, "lines": [], "events": [], "has_more_events": False}
    path = directory / "events.jsonl"
    if path.is_file():
        with path.open("rb") as handle:
            for raw in handle:
                if not raw.endswith(b"\n"):
                    break
                event = json.loads(raw)
                if event["seq"] > since_event:
                    if len(payload["events"]) == tail:
                        payload["has_more_events"] = True
                        break
                    payload["events"].append(event)
    log = directory / "output.log"
    if not events_only and log.is_file():
        with log.open("rb") as handle:
            handle.seek(max(0, log.stat().st_size - 32768))
            if handle.tell():
                handle.readline()
            payload["lines"] = list(deque(
                handle.read(32768).decode("utf-8", errors="replace").splitlines(),
                maxlen=tail,
            ))
    payload["next_event"] = payload["events"][-1]["seq"] if payload["events"] else since_event
    return payload


class JobReporter:
    def __init__(self, directory: Path):
        self.directory = directory
        self.job_id = directory.name
        self.state = _read(directory / "state.json")
        self.lock = threading.Lock()
        self.stop = threading.Event()
        self.thread = threading.Thread(target=self._heartbeat, daemon=True)

    def __enter__(self):
        pid = os.getpid()
        _, identity = _process(pid)
        self.state.update(status="running", pid=pid, process_identity=identity,
                          started_at=now_iso(), heartbeat_at=now_iso())
        self.emit("job_start", stage="startup")
        self.thread.start()
        return self

    def __exit__(self, *_):
        self.stop.set()
        self.thread.join()

    def _heartbeat(self) -> None:
        while not self.stop.wait(HEARTBEAT_SECONDS):
            with self.lock:
                if self.state["status"] == "running":
                    self.state["heartbeat_at"] = now_iso()
                    self._save()

    def _save(self) -> None:
        write_json(self.directory / "state.json", self.state)

    def emit(self, event: str, *, stage: str, detail: dict | None = None) -> None:
        with self.lock:
            self._emit(event, stage, detail or {})

    def _emit(self, event: str, stage: str, detail: dict) -> None:
        seq = self.state["last_event_seq"] + 1
        record = {"seq": seq, "timestamp": now_iso(), "event": event, "stage": stage, "detail": detail}
        with (self.directory / "events.jsonl").open("a", encoding="utf-8", newline="\n") as handle:
            handle.write(json.dumps(record, ensure_ascii=False, allow_nan=False) + "\n")
        self.state.update(last_event_seq=seq, progress={"stage": stage, **detail})
        if event == "artifact":
            self.state["artifacts"].setdefault(detail["experiment"], {}).update(detail)
        self._save()

    def check_cancelled(self) -> None:
        if (self.directory / "cancel.flag").exists():
            raise JobCancelledError("已响应取消请求")

    def finish(self, status: str, exit_code: int, *, result: dict | None = None, error: dict | None = None) -> None:
        with self.lock:
            self.state.update(status=status, exit_code=exit_code, finished_at=now_iso(), result=result, error=error)
            self._emit("job_end", status, {})
