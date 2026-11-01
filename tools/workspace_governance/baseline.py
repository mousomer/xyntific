"""Standalone governance baseline with bounded, content-safe diagnostic output."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from tools.workspace_governance.diagnostics import execute, render

CHECKS = (
    ("environment_doctor", ["./gov", "doctor"]),
    ("workspace_governance", ["./gov", "check", "--generated", "--summary", "--json"]),
    ("policy_compliance", ["/bin/bash", "scripts/check_policy_compliance.sh", "."]),
    ("git_sanitation", ["/bin/bash", "scripts/check_git_sanitation.sh", "."]),
    ("template_drift", ["/bin/bash", "scripts/check_policy_template_drift.sh", "."]),
)


def run_baseline(root: Path, *, timeout: float = 120) -> int:
    checks = [
        *CHECKS,
        (
            "pack_sanitation",
            [
                sys.executable,
                "-B",
                "tools/workspace_governance/sanitation.py",
                "--root",
                ".",
            ],
        ),
    ]
    # Output is private process data, parsed through a closed diagnostic grammar.
    for name, argv in checks:
        report, status = execute(
            name, argv, root, environment=dict(os.environ), timeout=timeout
        )
        print(render(report), flush=True)
        if status:
            return status if status > 0 else 1
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--timeout", type=float, default=120)
    args = parser.parse_args(argv)
    if not 0 < args.timeout <= 7200:
        parser.error("timeout must be positive and bounded")
    return run_baseline(args.root.resolve(), timeout=args.timeout)


if __name__ == "__main__":
    if __package__ in {None, ""}:
        raise SystemExit(
            "Invoke the baseline as a module with the selected runtime root"
        )
    raise SystemExit(main())
