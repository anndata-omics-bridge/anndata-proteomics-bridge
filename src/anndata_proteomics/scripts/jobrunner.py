"""Project-agnostic background-job runner for marimo apps (no marimo import).

Distilled from MetaboStatsHub `src/cd_convert/prolfqua_runner.py`. Drop into a project as
e.g. `yourpkg/jobrunner.py`, then drive it from a marimo file (see marimo_wiring.py).

The contract: start a command as a background subprocess streaming stdout+stderr to a log
file; inspect it (poll + log tail) without side effects; tail/zip outputs; terminate the whole
process tree; and serve outputs over a real static HTTP server so links open in the browser.
"""

from __future__ import annotations

import hashlib
import io
import os
import signal
import socket
import subprocess
import sys
import zipfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from time import time
from urllib.parse import quote


@dataclass
class Job:
    """A running (or finished) background subprocess job."""

    command: tuple[str, ...]
    process: subprocess.Popen
    log_file: Path
    output_dir: Path
    run_key: str | None = None
    started_at: float = 0.0
    preexisting: frozenset[Path] = field(default_factory=frozenset)


@dataclass(frozen=True)
class JobStatus:
    """Immutable snapshot of a job, built by inspect_job (safe to render)."""

    command: tuple[str, ...]
    returncode: int | None
    running: bool
    log_file: Path
    log_text: str
    output_dir: Path
    reports: tuple[Path, ...]
    run_key: str | None = None

    @property
    def success(self) -> bool:
        return self.returncode == 0


@dataclass
class StaticServer:
    """A `python -m http.server` process rooted at an output directory."""

    root_dir: Path
    host: str
    port: int
    process: subprocess.Popen

    @property
    def base_url(self) -> str:
        return f"http://{self.host}:{self.port}"

    @property
    def running(self) -> bool:
        return self.process.poll() is None


def make_run_key(*parts: object) -> str:
    """sha256 fingerprint of the inputs that define a run (detects selection changes)."""
    digest = hashlib.sha256()
    for part in parts:
        digest.update(str(part).encode("utf-8"))
        digest.update(b"\0")
    return digest.hexdigest()


def start_job(
    command: Sequence[str],
    output_dir: Path | str,
    *,
    log_file: Path | str | None = None,
    run_key: str | None = None,
    env: Mapping[str, str] | None = None,
    report_glob: str = "*.html",
    popen=subprocess.Popen,
) -> Job:
    """Launch `command` in the background, streaming stdout+stderr to a log file."""
    outdir = Path(output_dir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    log_path = (Path(log_file).expanduser().resolve()
                if log_file is not None else outdir / "console.log")
    log_path.parent.mkdir(parents=True, exist_ok=True)

    # Own process group/session so terminate_job can kill the whole tree (wrapper + children).
    group_kwargs: dict[str, object] = (
        {"creationflags": subprocess.CREATE_NEW_PROCESS_GROUP} if os.name == "nt"
        else {"start_new_session": True}
    )
    command = tuple(str(part) for part in command)
    preexisting = (frozenset(p.resolve() for p in outdir.rglob(report_glob) if p.is_file())
                   if outdir.exists() else frozenset())
    started = time()
    with log_path.open("w", encoding="utf-8") as handle:
        handle.write("$ " + " ".join(command) + "\n\n")
        handle.flush()
        process = popen(
            list(command),
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
            env=dict(env) if env is not None else None,
            **group_kwargs,
        )
    return Job(command=command, process=process, log_file=log_path, output_dir=outdir,
               run_key=run_key, started_at=started, preexisting=preexisting)


def read_text_tail(path: Path | str, max_log_chars: int = 40000) -> str:
    """Return the tail of a UTF-8 text file (with a truncation marker if clipped)."""
    file_path = Path(path).expanduser()
    if not file_path.exists():
        return ""
    text = file_path.read_text(encoding="utf-8", errors="replace")
    if len(text) <= max_log_chars:
        return text
    return "... log truncated ...\n" + text[-max_log_chars:]


def inspect_job(job: Job, *, report_glob: str = "*.html", max_log_chars: int = 40000) -> JobStatus:
    """Poll the job and read its log tail. Pure: no state mutation, safe to call every refresh."""
    returncode = job.process.poll()
    reports: tuple[Path, ...] = ()
    if job.output_dir.exists():
        found = sorted(job.output_dir.rglob(report_glob))
        # "fresh" = not present at start, or rewritten after start (mtime-granularity-immune).
        reports = tuple(
            p for p in found
            if p.resolve() not in job.preexisting or p.stat().st_mtime >= job.started_at
        )
    return JobStatus(
        command=job.command,
        returncode=returncode,
        running=returncode is None,
        log_file=job.log_file,
        log_text=read_text_tail(job.log_file, max_log_chars=max_log_chars),
        output_dir=job.output_dir,
        reports=reports,
        run_key=job.run_key,
    )


def _signal_group(process: subprocess.Popen, *, force: bool) -> bool:
    pid = getattr(process, "pid", None)
    if not isinstance(pid, int):
        return False
    if os.name == "nt":
        args = ["taskkill", "/T", "/PID", str(pid)] if not force else ["taskkill", "/F", "/T", "/PID", str(pid)]
        try:
            subprocess.run(args, capture_output=True, check=False)
            return True
        except OSError:
            return False
    try:
        pgid = os.getpgid(pid)
    except OSError:
        return False
    if pgid != pid:          # child is not its own group leader → caller falls back
        return False
    try:
        os.killpg(pgid, signal.SIGKILL if force else signal.SIGTERM)
        return True
    except OSError:
        return False


def terminate_job(job: Job | None, timeout: float = 5.0) -> bool:
    """Terminate a still-running job and its child tree. No-op (False) if already done/None."""
    if job is None or job.process.poll() is not None:
        return False
    proc = job.process
    if not _signal_group(proc, force=False):
        proc.terminate()
    try:
        proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        if not _signal_group(proc, force=True):
            proc.kill()
        try:
            proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            pass
    return True


def find_available_port(host: str = "127.0.0.1") -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.bind((host, 0))
        return int(sock.getsockname()[1])


def start_static_server(root_dir: Path | str, host: str = "127.0.0.1",
                        port: int | None = None, popen=subprocess.Popen) -> StaticServer:
    """Serve `root_dir` over HTTP so report/log links open in the browser (not file://)."""
    root = Path(root_dir).expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"output directory not found: {root}")
    port = port or find_available_port(host)
    process = popen(
        [sys.executable, "-m", "http.server", str(port), "--bind", host, "--directory", str(root)],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, text=True,
    )
    return StaticServer(root_dir=root, host=host, port=port, process=process)


def path_to_url(path: Path | str, server: StaticServer) -> str:
    """Map a file under the server root to its HTTP URL."""
    rel = Path(path).expanduser().resolve().relative_to(server.root_dir.expanduser().resolve())
    return f"{server.base_url}/" + "/".join(quote(part) for part in rel.parts)


def zip_dir_to_bytes(directory: Path | str, exclude_suffixes: Sequence[str] = ()) -> bytes:
    """Zip a directory in memory for mo.download (so logs/outputs are grabbable)."""
    root = Path(directory).expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"output directory not found: {root}")
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, mode="w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(root.rglob("*")):
            if path.is_file() and not (exclude_suffixes and path.name.endswith(tuple(exclude_suffixes))):
                zf.write(path, arcname=str(path.relative_to(root)))
    return buffer.getvalue()
