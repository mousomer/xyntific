"""Bounded secret-pattern checks over Git-visible text; semantic review remains required."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
from pathlib import Path

MAX_BYTES = 2 * 1024 * 1024
MAX_FINDINGS = 50
PATTERNS = (
    ("private-key", re.compile(r"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----")),
    ("aws-access-key", re.compile(r"\b(?:AKIA|ASIA)[A-Z0-9]{16}\b")),
    ("github-token", re.compile(r"\bgh[pousr]_[A-Za-z0-9]{36,}\b")),
    ("openai-token", re.compile(r"\bsk-(?:proj-)?[A-Za-z0-9_-]{40,}\b")),
)


def read_text(path: Path) -> tuple[str, str | None]:
    try:
        with path.open("rb") as stream:
            data = stream.read(MAX_BYTES + 1)
        if len(data) > MAX_BYTES:
            return "", "size-bound"
        if b"\0" in data:
            return "", "binary"
        return data.decode("utf-8"), None
    except UnicodeError:
        return "", "binary-or-encoding"
    except OSError:
        return "", "unreadable-text"


def bounded_append(collection: list, item: dict) -> None:
    if len(collection) < MAX_FINDINGS:
        collection.append(item)


def scan(root: Path) -> dict:
    root = root.resolve()
    listed = subprocess.run(
        [
            "git",
            "-C",
            str(root),
            "ls-files",
            "-z",
            "--cached",
            "--others",
            "--exclude-standard",
        ],
        capture_output=True,
        check=False,
        timeout=30,
    )
    if listed.returncode:
        raise ValueError("secret scan requires an observable Git checkout")
    findings, skipped, errors = [], [], []
    inspected = 0
    finding_count = 0
    skipped_count = 0
    error_count = 0
    for raw in sorted(set(listed.stdout.split(b"\0")) - {b""}):
        relative = raw.decode("utf-8", errors="surrogateescape")
        path = root / relative
        if not path.exists() and not path.is_symlink():
            continue  # A tracked deletion contains no working-tree bytes.
        if path.is_symlink() or not path.resolve().is_relative_to(root):
            skipped_count += 1
            bounded_append(skipped, {"path": relative, "reason": "symlink"})
            continue
        if not path.is_file():
            skipped_count += 1
            bounded_append(skipped, {"path": relative, "reason": "non-regular"})
            continue
        content, reason = read_text(path)
        if reason is not None:
            collection = errors if reason == "unreadable-text" else skipped
            error_count += int(reason == "unreadable-text")
            skipped_count += int(reason != "unreadable-text")
            bounded_append(collection, {"path": relative, "reason": reason})
            continue
        inspected += 1
        for pattern_id, pattern in PATTERNS:
            for match in pattern.finditer(content):
                finding_count += 1
                bounded_append(
                    findings,
                    {
                        "pattern_id": pattern_id,
                        "path": relative,
                        "line": content.count("\n", 0, match.start()) + 1,
                    },
                )
    return {
        "profile": "git-visible-text-patterns-v1",
        "inspected_files": inspected,
        "findings": findings[:MAX_FINDINGS],
        "finding_count": finding_count,
        "skipped": skipped[:MAX_FINDINGS],
        "skipped_count": skipped_count,
        "errors": errors[:MAX_FINDINGS],
        "error_count": error_count,
        "limits": "Known patterns only; binary, symlink and oversized files require other evidence. Ignored files and Git history are outside this check.",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    args = parser.parse_args()
    try:
        result = scan(args.root)
    except (OSError, ValueError, subprocess.TimeoutExpired):
        print(
            json.dumps(
                {
                    "status": "unavailable",
                    "reason": "Git-visible inputs could not be inspected",
                }
            )
        )
        return 2
    print(json.dumps(result, indent=2, ensure_ascii=True))
    return 1 if result["finding_count"] else 2 if result["error_count"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
