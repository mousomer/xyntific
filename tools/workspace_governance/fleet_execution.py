"""Bounded execution evidence without retaining or printing consumer output."""

from __future__ import annotations

import hashlib
import math
import os
import selectors
import signal
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Any

from tools.workspace_governance.fleet_support import (
    MAKE_ENVIRONMENT_VARIABLES,
    execution_policy,
)

COLLECTOR_ID = "workspace-governance-fleet-verify-v2"


def safe_environment(verification: dict[str, Any] | None = None) -> dict[str, str]:
    """Allow declared ambient names; values are never included in diagnostics."""
    policy = execution_policy(verification)
    environment = {
        name: os.environ[name] for name in policy["inherit"] if name in os.environ
    }
    environment["PATH"] = ":".join(policy["path"])
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    return environment


def make_environment(environment: dict[str, str]) -> dict[str, str]:
    return {
        name: environment[name]
        for name in MAKE_ENVIRONMENT_VARIABLES
        if name in environment
    }


def _execution_limits(timeout: float, diagnostic_limit: int) -> None:
    if type(diagnostic_limit) is not int or not 0 <= diagnostic_limit <= 128 * 1024:
        raise ValueError(
            "Diagnostic capture limit must be an integer from zero to 128 KiB"
        )
    if (
        isinstance(timeout, bool)
        or not isinstance(timeout, (int, float))
        or not math.isfinite(timeout)
        or timeout <= 0
    ):
        raise ValueError("Execution timeout must be finite and positive")


def _start_process(argv: list[str], root: Path, environment: dict[str, str]):
    startup = (
        tempfile.TemporaryDirectory(prefix="governance-zsh-startup-")
        if Path(argv[0]).name == "zsh"
        else None
    )
    if startup is not None:
        environment = {**environment, "ZDOTDIR": startup.name}
    try:
        process = subprocess.Popen(
            argv,
            cwd=root,
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            start_new_session=True,
        )
    except BaseException:
        if startup is not None:
            startup.cleanup()
        raise
    return process, startup


def run_digest(
    argv: list[str],
    root: Path,
    environment: dict[str, str],
    timeout: float,
    *,
    diagnostic_limit: int = 0,
) -> dict[str, Any]:
    """Drain both pipes in bounded chunks; kill only this process group on timeout."""
    _execution_limits(timeout, diagnostic_limit)
    started = time.monotonic()
    deadline = started + timeout
    streams = {
        name: {"hash": hashlib.sha256(), "bytes": 0, "preview": bytearray()}
        for name in ("stdout", "stderr")
    }
    timed_out = False
    process, startup = _start_process(argv, root, environment)
    with selectors.DefaultSelector() as selector:
        for name, stream in (("stdout", process.stdout), ("stderr", process.stderr)):
            os.set_blocking(stream.fileno(), False)
            selector.register(stream, selectors.EVENT_READ, name)
        while selector.get_map() or process.poll() is None:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                timed_out = True
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                break
            for key, _ in selector.select(min(remaining, 0.25)):
                data = os.read(key.fileobj.fileno(), 64 * 1024)
                if data:
                    streams[key.data]["hash"].update(data)
                    streams[key.data]["bytes"] += len(data)
                    streams[key.data]["preview"].extend(
                        data[
                            : max(
                                0, diagnostic_limit - len(streams[key.data]["preview"])
                            )
                        ]
                    )
                else:
                    selector.unregister(key.fileobj)
    process.wait()
    process.stdout.close()
    process.stderr.close()
    if startup is not None:
        startup.cleanup()
    return {
        "exit_code": process.returncode,
        "timed_out": timed_out,
        "duration_seconds": round(time.monotonic() - started, 6),
        **{f"{name}_bytes": info["bytes"] for name, info in streams.items()},
        **{
            f"{name}_sha256": info["hash"].hexdigest() for name, info in streams.items()
        },
        **(
            {
                f"{name}_preview": bytes(info["preview"])
                for name, info in streams.items()
            }
            if diagnostic_limit
            else {}
        ),
    }
