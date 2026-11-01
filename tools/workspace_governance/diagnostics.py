"""Bounded checker classifications; raw output and matched contents never become diagnostics.

The caller owns execution and its exit semantics. This module translates only closed
legacy headers or known structured metadata. Locations require a declared allowlist;
unknown output remains an observation count and digest, never displayed text.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import sys
from collections import Counter
from pathlib import Path, PurePosixPath
from typing import Any

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

MAX_INPUT_BYTES = 1024 * 1024
MAX_DIAGNOSTICS = 12
BASELINE_LOCATIONS = frozenset(
    {
        ".workspace_policy_version.json",
        "scripts/verify.sh",
        "scripts/check_git_sanitation.sh",
        "scripts/check_policy_compliance.sh",
        "scripts/check_policy_template_drift.sh",
        ".policy/policy_template_hashes.json",
        "docs",
        "scripts",
        ".policy",
    }
)
CHECKERS = {
    "policy_compliance": "Restore the declared required policy files/directories and track required policy sources.",
    "git_sanitation": "Review the reported path class against the owning sanitation policy; preserve its warning semantics.",
    "template_drift": "Restore approved template bytes or review and regenerate the template hash manifest through its owner.",
    "pack_sanitation": "Remove the reported secret pattern from the working tree without copying matched content into logs; resolve unreadable inputs.",
    "workspace_check": "Repair the declared manifest, authority, generated-output or pack relationship identified by the rule.",
    "environment": "Select and certify the declared interpreter and tools through the machine overlay.",
    "generated_outputs": "Regenerate through the declared renderer, preserve protected inputs and rerun its read-only check.",
}
ALIASES = {
    "scripts/check_policy_compliance.sh": "policy_compliance",
    "scripts/check_git_sanitation.sh": "git_sanitation",
    "scripts/check_policy_template_drift.sh": "template_drift",
}
GENERIC_ACTION = "Review and rerun the declared owning checker; repair its reported rule using the consumer's authority."

GOVERNANCE_RULES = {
    "BROKEN_REFERENCE": "Restore the missing or invalid declared source, pointer, skill adapter or renderer target.",
    "CONFLICTING_VALUE": "Reconcile the conflicting declaration or generated output with its owning authority; preserve protected inputs.",
    "PACK_DRIFT": "Restore the reviewed pack bytes and lock, or publish and adopt a reviewed successor release.",
    "DUPLICATE_AUTHORITY": "Choose one exclusive owner for the duplicated scope and repair its references.",
    "AMBIGUOUS_AUTHORITY": "Resolve overlapping routes or authority selection through the owning manifest.",
    "ENVIRONMENT_MISMATCH": "Repair the declared interpreter, runtime floor, package binding or tool selection.",
}
WARNING_HEADERS = {
    "WARN: absolute Unix paths found in git-visible text files": (
        "sanitation.absolute-unix",
        "Review absolute Unix path references through the sanitation owner.",
    ),
    "WARN: absolute Windows paths found in git-visible text files": (
        "sanitation.absolute-windows",
        "Review absolute Windows path references through the sanitation owner.",
    ),
    "WARN: '../' patterns found in git-visible text files": (
        "sanitation.parent-traversal",
        "Review parent traversal references against approved confined-path rules.",
    ),
    "WARN: filenames with spaces detected (tracked or untracked, excluding ignored)": (
        "sanitation.spaced-filenames",
        "Review filenames with spaces against the existing repository sanitation policy.",
    ),
}
PATTERN_IDS = frozenset(
    {"private-key", "aws-access-key", "github-token", "openai-token"}
)
SKIP_REASONS = frozenset(
    {
        "binary",
        "binary-or-encoding",
        "size-bound",
        "symlink",
        "non-regular",
        "unreadable-text",
    }
)


def _relative(value: Any) -> str | None:
    if not isinstance(value, str) or not value or len(value) > 240 or "\\" in value:
        return None
    path = PurePosixPath(value)
    if (
        path.is_absolute()
        or value != path.as_posix()
        or any(part in {"..", ".git"} for part in path.parts)
    ):
        return None
    return value


def _count(value: Any) -> int:
    return value if type(value) is int and 0 <= value <= 2**63 - 1 else 0


def summarize(  # noqa: C901 - closed legacy-header and structured-metadata classification transaction
    checker: str,
    returncode: int,
    stdout: bytes = b"",
    stderr: bytes = b"",
    *,
    allowed_locations=(),
    max_diagnostics: int = MAX_DIAGNOSTICS,
) -> dict[str, Any]:
    """Classify bounded known metadata without changing the checker's return code."""
    checker = ALIASES.get(checker, checker.replace("-", "_"))
    if (
        not re.fullmatch(r"[a-z][a-z0-9_]{0,63}", checker)
        or type(returncode) is not int
        or not 1 <= max_diagnostics <= MAX_DIAGNOSTICS
    ):
        raise ValueError("Unknown checker, invalid exit status or diagnostic bound")
    if not isinstance(stdout, bytes) or not isinstance(stderr, bytes):
        raise TypeError("Checker output must be captured bytes")
    action = CHECKERS.get(checker, GENERIC_ACTION)
    locations = BASELINE_LOCATIONS | frozenset(
        value for item in allowed_locations if (value := _relative(item))
    )
    data = stdout + b"\n" + stderr
    rows, counts = [], Counter()
    sampled_reasons = Counter()
    unknown = 0

    def add(rule, severity, remediation, location=None, line=None):
        counts[severity] += 1
        row = {"rule_id": rule, "severity": severity, "remediation": remediation}
        if location in locations and _relative(location):
            row["location"] = location
        if type(line) is int and 0 < line <= 2**31 - 1:
            row["line"] = line
        if len(rows) < max_diagnostics:
            rows.append(row)

    truncated = len(data) > MAX_INPUT_BYTES
    parsed = None
    if not truncated:
        try:
            parsed = json.loads(stdout.decode("utf-8"))
        except (ValueError, UnicodeError):
            pass
    metrics = {}
    if isinstance(parsed, dict) and parsed.get("kind") == "bounded-checker-diagnostics":
        items = parsed.get("diagnostics", [])
        for item in items if isinstance(items, list) else []:
            if isinstance(item, dict) and item.get("rule_id") in GOVERNANCE_RULES:
                add(
                    item["rule_id"],
                    "error",
                    GOVERNANCE_RULES[item["rule_id"]],
                    item.get("location"),
                    item.get("line"),
                )
    elif (
        checker == "pack_sanitation"
        and isinstance(parsed, dict)
        and parsed.get("profile") == "git-visible-text-patterns-v1"
    ):
        metrics = {
            name: _count(parsed.get(name))
            for name in (
                "inspected_files",
                "finding_count",
                "skipped_count",
                "error_count",
            )
        }
        findings = parsed.get("findings", [])
        errors = parsed.get("errors", [])
        skipped = parsed.get("skipped", [])
        for item in findings[:MAX_DIAGNOSTICS] if isinstance(findings, list) else []:
            if isinstance(item, dict) and item.get("pattern_id") in PATTERN_IDS:
                add(
                    "secret." + item["pattern_id"],
                    "error",
                    action,
                    item.get("path"),
                    item.get("line"),
                )
        for item in errors[:MAX_DIAGNOSTICS] if isinstance(errors, list) else []:
            if isinstance(item, dict) and item.get("reason") == "unreadable-text":
                add(
                    "scan.unreadable",
                    "error",
                    "Repair read access or record the protected-input coverage gap.",
                    item.get("path"),
                )
        for item in skipped if isinstance(skipped, list) else []:
            if isinstance(item, dict) and item.get("reason") in SKIP_REASONS:
                sampled_reasons[item["reason"]] += 1
        if metrics.get("finding_count", 0) and counts["error"] == 0:
            add("scan.findings", "error", action)
        if metrics.get("error_count", 0) and not errors:
            add(
                "scan.unreadable",
                "error",
                "Repair read access or record the protected-input coverage gap.",
            )
    elif isinstance(parsed, (dict, list)):
        items = parsed if isinstance(parsed, list) else parsed.get("diagnostics", [])
        if isinstance(items, list):
            for item in items:
                if isinstance(item, dict) and item.get("code") in GOVERNANCE_RULES:
                    sources = item.get("sources", [])
                    source = (
                        next(
                            (
                                value.partition("#")[0]
                                for value in sources
                                if isinstance(value, str)
                                and value.partition("#")[0] in locations
                            ),
                            None,
                        )
                        if isinstance(sources, list)
                        else None
                    )
                    add(item["code"], "error", GOVERNANCE_RULES[item["code"]], source)
    if parsed is None or stderr:
        plain = data if parsed is None else stderr
        for line in (
            plain[:MAX_INPUT_BYTES].decode("utf-8", errors="replace").splitlines()
        ):
            if line in WARNING_HEADERS:
                rule, rule_action = WARNING_HEADERS[line]
                add(rule, "warning", rule_action)
                continue
            match = re.fullmatch(
                r"Missing required (?:policy )?(file|directory): (.+)", line
            )
            if match:
                add("policy.missing-" + match[1], "error", action, match[2])
                continue
            match = re.fullmatch(
                r"Required policy file must be tracked in git: (.+)", line
            )
            if match:
                add("policy.untracked-source", "error", action, match[1])
                continue
            if checker == "template_drift":
                if (
                    line
                    == "FAIL: missing template hash manifest: .policy/policy_template_hashes.json"
                ):
                    add(
                        "template.missing-manifest",
                        "error",
                        action,
                        ".policy/policy_template_hashes.json",
                    )
                    continue
                if line.startswith(
                    "FAIL: unsupported hash algorithm in .policy/policy_template_hashes.json: "
                ):
                    add(
                        "template.hash-algorithm",
                        "error",
                        "Use the declared SHA-256 template manifest format.",
                        ".policy/policy_template_hashes.json",
                    )
                    continue
                match = re.fullmatch(r"FAIL: missing template-tracked file: (.+)", line)
                if match:
                    add("template.missing-source", "error", action, match[1])
                    continue
                match = re.fullmatch(
                    r"FAIL: template drift detected for (.+) \([^\n]*\)", line
                )
                if match:
                    add("template.content-drift", "error", action, match[1])
                    continue
            if re.fullmatch(
                r"WARN: sanitation checks reported \d+ issue class\(es\)", line
            ):
                continue
            if line.startswith(("WARN:", "WARNING:")):
                add("checker.warning", "warning", action)
                continue
            if line:
                unknown += 1
    if returncode and counts["error"] == 0:
        add("checker.failed", "error", action)
    return {
        "schema_version": 1,
        "kind": "bounded-checker-diagnostics",
        "checker": checker,
        "exit_status": returncode,
        "status": "failed"
        if returncode
        else "passed_with_warnings"
        if counts["warning"]
        else "passed",
        "diagnostics": rows,
        "error_count": counts["error"],
        "warning_count": counts["warning"],
        "suppressed_diagnostics": max(0, sum(counts.values()) - len(rows)),
        "unclassified_output_lines": unknown,
        "metrics": metrics,
        "sampled_skip_reasons": dict(sorted(sampled_reasons.items())),
        "output_bytes": len(stdout) + len(stderr),
        "output_sha256": hashlib.sha256(stdout + stderr).hexdigest(),
        "classification_input_truncated": truncated,
    }


def render(summary: dict[str, Any]) -> str:
    """At most fourteen lines; escaped locations and fixed vocabulary only."""
    header = f"{summary['checker']}: {summary['status']} (exit={summary['exit_status']}; errors={summary['error_count']}; warnings={summary['warning_count']})"
    if summary["metrics"]:
        header += "; " + "; ".join(
            f"{key}={value}" for key, value in summary["metrics"].items()
        )
    lines = [header]
    for item in summary["diagnostics"]:
        where = (
            " location=" + json.dumps(item["location"], ensure_ascii=True)
            if "location" in item
            else ""
        )
        if "line" in item:
            where += ":" + str(item["line"])
        lines.append(
            f"{item['severity']}: {item['rule_id']}{where}; {item['remediation']}"
        )
    if summary["suppressed_diagnostics"]:
        lines.append(
            f"additional_diagnostics={summary['suppressed_diagnostics']} (bounded output)"
        )
    return "\n".join(lines)


def execute(
    checker: str,
    argv: list[str],
    root: Path,
    *,
    environment=None,
    timeout: float = 120,
    allowed_locations=(),
) -> tuple[dict[str, Any], int]:
    """Run a governed checker with bounded in-memory diagnostics and no raw log files."""
    from tools.workspace_governance.fleet_execution import run_digest

    summarize(checker, 0, allowed_locations=allowed_locations)
    if (
        not isinstance(argv, list)
        or not argv
        or any(not isinstance(arg, str) or "\0" in arg for arg in argv)
    ):
        raise ValueError("Checker argv must be a nonempty string list")
    if (
        isinstance(timeout, bool)
        or not isinstance(timeout, (int, float))
        or not math.isfinite(timeout)
        or not 0 < timeout <= 86400
    ):
        raise ValueError("Checker timeout is outside its finite bound")
    # Preserve the already-certified source binding; refuse startup/Make injection.
    prohibited = {
        "BASH_ENV",
        "ENV",
        "SHELLOPTS",
        "BASHOPTS",
        "CDPATH",
        "IFS",
        "ZDOTDIR",
        "MAKEFLAGS",
        "MFLAGS",
        "MAKEFILES",
        "GNUMAKEFLAGS",
    }
    source = dict(os.environ if environment is None else environment)
    effective = {
        name: value for name, value in source.items() if name not in prohibited
    }
    try:
        result = run_digest(argv, root, effective, timeout, diagnostic_limit=128 * 1024)
    except OSError:
        report = summarize(checker, 2, allowed_locations=allowed_locations)
        report["diagnostics"][0] = {
            "rule_id": "checker.unavailable",
            "severity": "error",
            "remediation": "Restore the declared checker executable and working directory, then rerun the owning check.",
        }
        report["execution_unavailable"] = True
        return report, 2
    code = 124 if result["timed_out"] else result["exit_code"]
    if code < 0:
        code = 128 - code
    stdout, stderr = (
        result.pop("stdout_preview", b""),
        result.pop("stderr_preview", b""),
    )
    report = summarize(
        checker, code, stdout, stderr, allowed_locations=allowed_locations
    )
    if result["timed_out"]:
        fallback = next(
            (
                row
                for row in report["diagnostics"]
                if row["rule_id"] == "checker.failed"
            ),
            None,
        )
        if fallback is not None:
            report["diagnostics"].remove(fallback)
            report["error_count"] -= 1
        report["error_count"] += 1
        if len(report["diagnostics"]) < MAX_DIAGNOSTICS:
            report["diagnostics"].append(
                {
                    "rule_id": "checker.timeout",
                    "severity": "error",
                    "remediation": "Review the declared checker timeout and incomplete execution, then rerun the owning check.",
                }
            )
        else:
            report["suppressed_diagnostics"] += 1
    report.update(
        execution=result,
        classification_input_truncated=(
            report["classification_input_truncated"]
            or len(stdout) < result["stdout_bytes"]
            or len(stderr) < result["stderr_bytes"]
        ),
        counts_scope="bounded captured prefix"
        if len(stdout) < result["stdout_bytes"] or len(stderr) < result["stderr_bytes"]
        else "complete checker output",
    )
    return report, code


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checker", required=True)
    parser.add_argument("--exit-status", type=int)
    parser.add_argument("--allowed-location", action="append", default=[])
    parser.add_argument("--timeout", type=float, default=120)
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--run", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    try:
        if args.run is not None:
            if args.exit_status is not None:
                parser.error("Execution and captured-input modes are exclusive")
            result, code = execute(
                args.checker,
                args.run,
                Path.cwd(),
                timeout=args.timeout,
                allowed_locations=args.allowed_location,
            )
        else:
            if args.exit_status is None:
                parser.error("Captured-input mode requires --exit-status")
            data = sys.stdin.buffer.read(MAX_INPUT_BYTES + 1)
            result = summarize(
                args.checker,
                args.exit_status,
                data,
                allowed_locations=args.allowed_location,
            )
            code = 0
    except (TypeError, ValueError):
        parser.error("Invalid bounded diagnostic contract")
    print(
        json.dumps(result, ensure_ascii=True, separators=(",", ":"))
        if args.json
        else render(result)
    )
    return code


if __name__ == "__main__":
    raise SystemExit(main())
