"""Guarded fleet plan/apply/audit/undo. Machine registries remain private."""

from __future__ import annotations

import argparse
import base64
import json
import math
import re
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path, PurePosixPath
from typing import Any
from uuid import uuid4

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tools.workspace_governance.fleet_execution import (
    COLLECTOR_ID,
    make_environment,
    run_digest,
    safe_environment,
)
from tools.workspace_governance.fleet_support import (
    BEGIN,
    END,
    EXECUTION_POLICY,
    PROTECTED,
    SCHEMA,
    FleetError,
    atomic,
    canonical,
    confined,
    digest,
    execution_policy,
    git_identity,
    instruction_aliases,
    load,
    registry_rows,
    relative,
    release_files,
    remove,
    rendered,
    resolved_verification,
    state,
    validate_verification_argv,
)

__all__ = [
    "BEGIN",
    "END",
    "FleetError",
    "apply",
    "atomic",
    "audit",
    "canonical",
    "collect_acceptance",
    "collect_baselines",
    "digest",
    "git_identity",
    "main",
    "merge_evidence",
    "plan",
    "relative",
    "release_files",
    "undo",
]


def plan(registry_path: Path, release_path: Path) -> dict[str, Any]:
    registry = load(registry_path)
    release, identity, payload = release_files(release_path)
    output: dict[str, Any] = {
        "schema_version": SCHEMA,
        "kind": "fleet-adoption-plan",
        "registry_sha256": digest(canonical(registry)),
        "registry_path": str(registry_path.absolute()),
        "release_path": str(release_path.absolute()),
        "release_manifest_sha256": identity,
        "scope": registry.get("scope", {}),
        "checkouts": [],
    }
    modes = {entry["path"]: entry.get("mode", 0o644) for entry in release["files"]}
    for row in registry_rows(registry):
        item = {
            "id": row["id"],
            "path": row["path"],
            "repository_id": row["repository_id"],
            "integration_mode": row["integration_mode"],
            "operations": [],
            "errors": [],
        }
        output["checkouts"].append(item)
        if row["integration_mode"] in {"external", "blocked"}:
            item["status"] = (
                "external" if row["integration_mode"] == "external" else "blocked"
            )
            continue
        try:
            if row.get("classification") in PROTECTED:
                raise FleetError(
                    "Protected repository classification cannot be mutated"
                )
            root = Path(row["path"]).absolute()
            item["git"] = git_identity(root)
            item["instruction_aliases"] = instruction_aliases(row, root)
            desired = {name: (data, modes[name]) for name, data in payload.items()}
            for name, entry in rendered(row, root, release, identity).items():
                if name in desired:
                    raise FleetError(f"Release and binding both own {name}")
                desired[name] = entry
            for name, (data, mode) in desired.items():
                current = state(confined(root, name))
                after = {"sha256": digest(data), "mode": mode}
                if current == after:
                    continue
                section_owned = name in {
                    "AGENTS.md",
                    "CLAUDE.md",
                    row["verification"]["runner"]["path"],
                }
                replace = row.get("replace_sha256", {}).get(name)
                if (
                    current["sha256"] is not None
                    and not section_owned
                    and replace != current["sha256"]
                ):
                    raise FleetError(f"Unapproved preexisting file: {name}")
                item["operations"].append(
                    {
                        "path": name,
                        "before": current,
                        "after": after,
                        "content_base64": base64.b64encode(data).decode(),
                    }
                )
            item["status"] = "ready" if item["operations"] else "unchanged"
        except (
            FleetError,
            OSError,
            KeyError,
            UnicodeError,
            subprocess.TimeoutExpired,
        ) as exc:
            item["status"] = "blocked"
            item["errors"].append(str(exc))
            item["operations"] = []
    return output


def restore(root: Path, operation: dict[str, Any]) -> None:
    path = confined(root, operation["path"])
    if state(path) != operation["after"]:
        raise FleetError(f"Rollback refuses changed file: {operation['path']}")
    before = operation["before"]
    if before["sha256"] is None:
        remove(path)
    else:
        data = base64.b64decode(operation["backup_base64"], validate=True)
        if digest(data) != before["sha256"]:
            raise FleetError("Rollback backup hash mismatch")
        atomic(path, data, before["mode"])


def pending_operations(
    item: dict[str, Any],
    approved: dict[str, tuple[bytes, int]],
    root: Path,
    row: dict[str, Any],
) -> list[tuple[dict[str, Any], bytes, str | None]]:
    pending = []
    paths = [operation["path"] for operation in item["operations"]]
    if len(paths) != len(set(paths)):
        raise FleetError("Plan contains duplicate managed paths")
    for operation in item["operations"]:
        path = confined(root, operation["path"])
        current = state(path)
        if current == operation["after"]:
            continue
        section_owned = operation["path"] in {
            "AGENTS.md",
            "CLAUDE.md",
            row["verification"]["runner"]["path"],
        }
        if (
            current["sha256"] is not None
            and not section_owned
            and row.get("replace_sha256", {}).get(operation["path"])
            != current["sha256"]
        ):
            raise FleetError(
                f"Registry does not authorize replacing: {operation['path']}"
            )
        if current != operation["before"]:
            raise FleetError(f"Preimage changed: {operation['path']}")
        data = base64.b64decode(operation["content_base64"], validate=True)
        if digest(data) != operation["after"]["sha256"] or approved.get(
            operation["path"]
        ) != (data, operation["after"]["mode"]):
            raise FleetError("Plan content differs from approved release or binding")
        backup = base64.b64encode(path.read_bytes()).decode() if path.exists() else None
        pending.append((operation, data, backup))
    return pending


def rollback_writes(record: dict[str, Any], root: Path) -> None:
    for saved in reversed(record["writes"]):
        try:
            if state(confined(root, saved["path"])) == saved["after"]:
                restore(root, saved)
                saved["status"] = "rolled_back"
        except (FleetError, OSError) as error:
            record["errors"].append(str(error))


def apply_checkout(
    item: dict[str, Any],
    row: dict[str, Any],
    record: dict[str, Any],
    receipt: dict[str, Any],
    receipt_path: Path,
    release: dict[str, Any],
    identity: str,
    payload: dict[str, bytes],
) -> None:
    root = Path(item["path"])
    try:
        if (
            row["path"] != item["path"]
            or row.get("classification") in PROTECTED
            or row["integration_mode"] not in {"direct", "inherited_checkout"}
        ):
            raise FleetError("Plan location differs from approved registry")
        if git_identity(root)["common_dir_sha256"] != item["git"]["common_dir_sha256"]:
            raise FleetError("Repository identity changed since planning")
        aliases = instruction_aliases(row, root)
        if aliases != item.get("instruction_aliases", {}):
            raise FleetError("Instruction alias changed since planning")
        record["instruction_aliases"] = aliases
        approved = {
            entry["path"]: (payload[entry["path"]], entry.get("mode", 0o644))
            for entry in release["files"]
        }
        approved.update(rendered(row, root, release, identity))
        pending = pending_operations(item, approved, root, row)
        for operation, data, backup in pending:
            if instruction_aliases(row, root) != aliases:
                raise FleetError("Instruction alias changed during application")
            saved = {k: operation[k] for k in ("path", "before", "after")}
            saved["backup_base64"] = backup
            saved["status"] = "pending"
            record["writes"].append(saved)
            atomic(receipt_path, canonical(receipt))
            atomic(confined(root, operation["path"]), data, operation["after"]["mode"])
            saved["status"] = "written"
            atomic(receipt_path, canonical(receipt))
        if instruction_aliases(row, root) != aliases:
            raise FleetError("Instruction alias changed during application")
        record["status"] = "installed" if pending else "unchanged"
    except (FleetError, OSError, KeyError, subprocess.TimeoutExpired) as error:
        record["status"] = "failed"
        record["errors"].append(str(error))
        rollback_writes(record, root)


def apply(plan_path: Path, receipt_path: Path) -> dict[str, Any]:
    declaration = load(plan_path)
    if declaration.get("kind") != "fleet-adoption-plan":
        raise FleetError("Not an adoption plan")
    release, identity, payload = release_files(
        Path(declaration["release_path"]), declaration["release_manifest_sha256"]
    )
    registry = load(Path(declaration["registry_path"]))
    if digest(canonical(registry)) != declaration["registry_sha256"]:
        raise FleetError("Registry changed since planning")
    rows = {row["id"]: row for row in registry_rows(registry)}
    planned_ids = [row["id"] for row in declaration["checkouts"]]
    if len(planned_ids) != len(set(planned_ids)) or set(planned_ids) != set(rows):
        raise FleetError(
            "Plan does not account for every registry checkout exactly once"
        )
    if receipt_path.exists():
        raise FleetError("Receipt already exists; preserve previous undo evidence")
    receipt: dict[str, Any] = {
        "schema_version": SCHEMA,
        "kind": "fleet-adoption-receipt",
        "plan_sha256": digest(canonical(declaration)),
        "checkouts": [],
    }
    atomic(receipt_path, canonical(receipt))
    for item in declaration["checkouts"]:
        record = {
            "id": item["id"],
            "path": item["path"],
            "status": item["status"],
            "writes": [],
            "errors": [],
        }
        receipt["checkouts"].append(record)
        if item["status"] in {"ready", "unchanged"}:
            apply_checkout(
                item,
                rows[item["id"]],
                record,
                receipt,
                receipt_path,
                release,
                identity,
                payload,
            )
        atomic(receipt_path, canonical(receipt))
    return receipt


def undo(receipt_path: Path) -> dict[str, Any]:
    receipt = load(receipt_path)
    if receipt.get("kind") != "fleet-adoption-receipt":
        raise FleetError("Not an adoption receipt")
    for record in reversed(receipt["checkouts"]):
        root = Path(record["path"])
        aliases = record.get("instruction_aliases", {})
        try:
            expected = {name: item["target"] for name, item in aliases.items()}
            if instruction_aliases({"instruction_aliases": expected}, root) != aliases:
                raise FleetError("Undo refuses changed instruction alias")
        except (FleetError, OSError, KeyError) as exc:
            record["errors"].append(str(exc))
            record["status"] = "undo_refused"
            atomic(receipt_path, canonical(receipt))
            continue
        for operation in reversed(record["writes"]):
            if operation["status"] in {"undone", "rolled_back"}:
                continue
            try:
                current = state(confined(root, operation["path"]))
                if current == operation["before"]:
                    operation["status"] = "undone"
                else:
                    restore(root, operation)
                    operation["status"] = "undone"
            except (FleetError, OSError) as exc:
                record["errors"].append(str(exc))
                record["status"] = "undo_refused"
        if record["writes"] and all(
            op["status"] in {"undone", "rolled_back"} for op in record["writes"]
        ):
            record["status"] = "undone"
        atomic(receipt_path, canonical(receipt))
    return receipt


def installed_metadata(
    row: dict[str, Any],
    root: Path,
    release: dict[str, Any],
    identity: str,
    payload: dict[str, bytes],
) -> tuple[dict[str, Any], str]:
    git = git_identity(root)
    for entry in release["files"]:
        name = entry["path"]
        expected = {"sha256": digest(payload[name]), "mode": entry.get("mode", 0o644)}
        if state(confined(root, name)) != expected:
            raise FleetError(f"Release file absent or drifted: {name}")
    binding = confined(root, ".governance/fleet.json")
    if load(binding).get("release_manifest_sha256") != identity:
        raise FleetError("Binding release identity differs")
    for name, (data, mode) in rendered(row, root, release, identity).items():
        if state(confined(root, name)) != {"sha256": digest(data), "mode": mode}:
            raise FleetError(f"Managed binding drifted: {name}")
    return git, digest(binding.read_bytes())


EVIDENCE_SCHEMA = 2
EVIDENCE_KIND = "fleet-verification-evidence"
EVIDENCE_TRUST = "collector-reported"


def verification_timeout(
    verification: dict[str, Any], key: str, default: float
) -> float:
    value = verification.get(key, default)
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value <= 0
    ):
        raise FleetError("Verification timeout must be a finite positive number")
    return float(value)


def verification_snapshot(git: dict[str, Any], binding: str) -> dict[str, Any]:
    return {
        name: git[name]
        for name in ("common_dir_sha256", "head", "head_state", "state_sha256")
    } | {"binding_sha256": binding}


def finish_verification_identity(record: dict[str, Any], root: Path) -> None:
    git = git_identity(root)
    binding = state(confined(root, ".governance/fleet.json"))["sha256"]
    current = verification_snapshot(git, binding)
    record.update(
        post_verification=current,
        input_state_stable=record["pre_verification"] == current,
        head=git["head"],
        head_state=git["head_state"],
        state_sha256=git["state_sha256"],
        binding_sha256=binding,
        completed_at=datetime.now(UTC).isoformat(),
    )


def baseline_command(row: dict[str, Any], root: Path) -> list[str]:
    name = relative(row["verification"]["runner"].get("governance_runner", "gov"))
    return [str(confined(root, name))] + (["check"] if name == "gov" else [])


def execution_contract(
    verification: dict[str, Any], root: Path, argv: list[str]
) -> tuple[list[str], dict[str, str], dict[str, Any]]:
    executed, assignments, tool = resolved_verification(verification, root, argv=argv)
    environment = safe_environment(verification)
    environment.update(assignments)
    contract = {
        "policy": EXECUTION_POLICY,
        "declaration_sha256": digest(canonical(execution_policy(verification))),
        "tool": tool,
        "inherited_names": sorted(
            name
            for name in execution_policy(verification)["inherit"]
            if name in environment
        ),
    }
    return executed, environment, contract


def collect_checkout(
    row: dict[str, Any],
    record: dict[str, Any],
    release: dict[str, Any],
    identity: str,
    payload: dict[str, bytes],
    *,
    baseline_only: bool = False,
) -> None:
    if row["integration_mode"] == "external":
        record["status"] = "not_applicable"
        return
    if (
        row["integration_mode"] not in {"direct", "inherited_checkout"}
        or row.get("classification") in PROTECTED
    ):
        raise FleetError("Checkout is not eligible for verification execution")
    root = Path(row["path"])
    verification = row["verification"]
    validate_verification_argv(verification, root)
    before_git, before_binding = installed_metadata(
        row, root, release, identity, payload
    )
    record["pre_verification"] = verification_snapshot(before_git, before_binding)
    declared_baseline = baseline_command(row, root)
    actual_baseline, environment, contract = execution_contract(
        verification, root, declared_baseline
    )
    record.update(
        governance_argv=declared_baseline,
        governance_executed_argv=actual_baseline,
        governance_execution=contract,
        make_environment=make_environment(environment),
    )
    baseline = run_digest(
        actual_baseline,
        root,
        environment,
        verification_timeout(verification, "governance_timeout_seconds", 120),
    )
    record.update({f"governance_{key}": value for key, value in baseline.items()})
    record["baseline_status"] = (
        "passed"
        if baseline["exit_code"] == 0 and not baseline["timed_out"]
        else "failed"
    )
    if record["baseline_status"] != "passed":
        finish_verification_identity(record, root)
        raise FleetError("Required governance baseline failed or timed out")
    if baseline_only:
        finish_verification_identity(record, root)
        installed_metadata(row, root, release, identity, payload)
        if not record["input_state_stable"]:
            raise FleetError(
                "Governance baseline changed Git-visible inputs or the binding"
            )
        record.update(
            status="baseline_passed", product_execution="not_requested", exit_code=None
        )
        return
    if verification.get("coverage") == "governance-only":
        finish_verification_identity(record, root)
        record.update(
            verification_coverage="governance-only",
            product_execution="unavailable",
            exit_code=None,
        )
        raise FleetError(
            "Product acceptance unavailable: verification coverage is governance-only"
        )
    actual, environment, contract = execution_contract(
        verification, root, verification["argv"]
    )
    record.update(
        argv=verification["argv"],
        executed_argv=actual,
        execution=contract,
        product_execution="requested",
    )
    record.update(
        run_digest(
            actual,
            root,
            environment,
            verification_timeout(verification, "timeout_seconds", 600),
        )
    )
    finish_verification_identity(record, root)
    if not record["input_state_stable"]:
        raise FleetError(
            "Verification changed Git-visible inputs or the governance binding"
        )
    installed_metadata(row, root, release, identity, payload)
    if record["exit_code"] != 0 or record["timed_out"]:
        raise FleetError("Declared product verification failed or timed out")
    record["status"] = "passed"


def collect_acceptance(
    registry_path: Path,
    release_path: Path,
    *,
    checkout_ids: list[str] | None = None,
    output_path: Path | None = None,
    baseline_only: bool = False,
) -> dict[str, Any]:
    """Collector-reported attestations; this UID can rewrite collector and proof."""
    registry = load(registry_path)
    release, identity, payload = release_files(release_path)
    rows = registry_rows(registry)
    selected = (
        set(checkout_ids) if checkout_ids is not None else {row["id"] for row in rows}
    )
    if checkout_ids is not None and len(checkout_ids) != len(selected):
        raise FleetError("Verification selected duplicate checkouts")
    if not selected or not selected.issubset({row["id"] for row in rows}):
        raise FleetError("Verification selected no checkout or an unknown checkout")
    if output_path is not None and any(
        output_path.absolute()
        .resolve()
        .is_relative_to(Path(row["path"]).absolute().resolve())
        for row in rows
    ):
        raise FleetError(
            "Verification evidence must remain outside consumer repositories"
        )
    mode = "baseline" if baseline_only else "product"
    result = {
        "schema_version": EVIDENCE_SCHEMA,
        "kind": EVIDENCE_KIND,
        "trust": EVIDENCE_TRUST,
        "collector_id": COLLECTOR_ID,
        "execution_policy": EXECUTION_POLICY,
        "release_manifest_sha256": identity,
        "registry_sha256": digest(canonical(registry)),
        "run_id": uuid4().hex,
        "mode": mode,
        "selection": [row["id"] for row in rows if row["id"] in selected],
        "started_at": datetime.now(UTC).isoformat(),
        "completed_at": None,
        "completed": False,
        "checkouts": [],
        "failures": [],
        "ok": False,
    }
    if output_path is not None:
        atomic(output_path, canonical(result))
    for row in rows:
        if row["id"] not in selected:
            continue
        record = {
            "id": row["id"],
            "release_manifest_sha256": identity,
            "collector_id": COLLECTOR_ID,
            "execution_policy": EXECUTION_POLICY,
            "run_id": result["run_id"],
            "mode": mode,
            "status": "failed",
            "errors": [],
            "started_at": datetime.now(UTC).isoformat(),
        }
        result["checkouts"].append(record)
        try:
            collect_checkout(
                row, record, release, identity, payload, baseline_only=baseline_only
            )
        except (
            FleetError,
            OSError,
            KeyError,
            UnicodeError,
            ValueError,
            subprocess.TimeoutExpired,
        ) as error:
            record["errors"].append(str(error))
            result["failures"].append(row["id"])
        finally:
            record["completed_at"] = datetime.now(UTC).isoformat()
        if output_path is not None:
            atomic(output_path, canonical(result))
    result.update(
        completed=True,
        completed_at=datetime.now(UTC).isoformat(),
        ok=not result["failures"],
    )
    if output_path is not None:
        atomic(output_path, canonical(result))
    return result


def collect_baselines(
    registry_path: Path, release_path: Path, **kwargs: Any
) -> dict[str, Any]:
    return collect_acceptance(registry_path, release_path, baseline_only=True, **kwargs)


def _timestamp(value: Any) -> datetime:
    if not isinstance(value, str):
        raise FleetError("Evidence timestamp must be an ISO UTC string")
    try:
        parsed = datetime.fromisoformat(value)
    except ValueError as error:
        raise FleetError("Malformed evidence timestamp") from error
    if parsed.tzinfo is None or parsed.utcoffset() != UTC.utcoffset(parsed):
        raise FleetError("Evidence timestamp must have UTC timezone")
    return parsed


def _hash(value: Any) -> bool:
    return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None


def _execution_record(record: dict[str, Any], prefix: str = "") -> None:
    for name in ("stdout_sha256", "stderr_sha256"):
        if not _hash(record.get(prefix + name)):
            raise FleetError("Evidence lacks valid execution output identities")
    for name in ("stdout_bytes", "stderr_bytes"):
        value = record.get(prefix + name)
        if type(value) is not int or value < 0:
            raise FleetError("Evidence output sizes must be nonnegative integers")
    if (
        type(record.get(prefix + "exit_code")) is not int
        or type(record.get(prefix + "timed_out")) is not bool
    ):
        raise FleetError("Evidence execution outcome has invalid types")
    duration = record.get(prefix + "duration_seconds")
    if (
        isinstance(duration, bool)
        or not isinstance(duration, (int, float))
        or not math.isfinite(duration)
        or duration < 0
    ):
        raise FleetError("Evidence execution duration must be finite and nonnegative")


def _snapshot_record(value: Any) -> None:
    keys = {"common_dir_sha256", "head", "head_state", "state_sha256", "binding_sha256"}
    if (
        not isinstance(value, dict)
        or set(value) != keys
        or any(not _hash(value.get(name)) for name in keys - {"head", "head_state"})
    ):
        raise FleetError("Evidence snapshot identity is malformed")
    if value["head_state"] == "unborn":
        if value["head"] is not None:
            raise FleetError("Unborn evidence snapshot has a commit")
    elif (
        value["head_state"] != "commit"
        or not isinstance(value["head"], str)
        or re.fullmatch(r"[0-9a-f]{40}|[0-9a-f]{64}", value["head"]) is None
    ):
        raise FleetError("Evidence snapshot commit identity is malformed")


def _contract_record(value: Any) -> None:
    if (
        not isinstance(value, dict)
        or set(value) != {"policy", "declaration_sha256", "tool", "inherited_names"}
        or value["policy"] != EXECUTION_POLICY
        or not _hash(value["declaration_sha256"])
    ):
        raise FleetError("Evidence execution contract is malformed")
    names = value["inherited_names"]
    if (
        not isinstance(names, list)
        or any(not isinstance(name, str) for name in names)
        or names != sorted(set(names))
    ):
        raise FleetError("Evidence inherited environment declaration is malformed")
    tool = value["tool"]
    if (
        not isinstance(tool, dict)
        or set(tool) != {"path", "resolved_path", "sha256", "mode"}
        or not _hash(tool["sha256"])
        or type(tool["mode"]) is not int
        or not 0 <= tool["mode"] <= 0o7777
        or any(
            not isinstance(tool[key], str) or not Path(tool[key]).is_absolute()
            for key in ("path", "resolved_path")
        )
    ):
        raise FleetError("Evidence executable identity is malformed")


def _evidence_header(
    value: dict[str, Any], registry: dict[str, Any], identity: str
) -> tuple[dict[str, Any], list[dict[str, Any]], datetime, datetime]:
    expected = {
        "schema_version": EVIDENCE_SCHEMA,
        "kind": EVIDENCE_KIND,
        "trust": EVIDENCE_TRUST,
        "collector_id": COLLECTOR_ID,
        "execution_policy": EXECUTION_POLICY,
        "release_manifest_sha256": identity,
        "registry_sha256": digest(canonical(registry)),
        "completed": True,
    }
    if (
        not isinstance(value, dict)
        or type(value.get("schema_version")) is not int
        or any(value.get(name) != item for name, item in expected.items())
        or type(value.get("completed")) is not bool
    ):
        raise FleetError(
            "Evidence envelope is historic, incomplete, inconsistent or uses different pins"
        )
    if (
        not isinstance(value.get("run_id"), str)
        or re.fullmatch(r"[0-9a-f]{32}", value["run_id"]) is None
    ):
        raise FleetError("Evidence run identity is invalid")
    if value.get("mode") not in {"product", "baseline", "mixed"}:
        raise FleetError("Evidence execution mode is invalid")
    started, completed = (
        _timestamp(value.get("started_at")),
        _timestamp(value.get("completed_at")),
    )
    if completed < started:
        raise FleetError("Evidence run completed before it started")
    selection = value.get("selection")
    known = {row["id"] for row in registry_rows(registry)}
    if (
        not isinstance(selection, list)
        or not selection
        or any(not isinstance(item, str) for item in selection)
        or len(set(selection)) != len(selection)
        or not set(selection).issubset(known)
    ):
        raise FleetError(
            "Evidence selection contains missing, duplicate or unknown IDs"
        )
    records = value.get("checkouts")
    if (
        not isinstance(records, list)
        or any(not isinstance(record, dict) for record in records)
        or [record.get("id") for record in records] != selection
    ):
        raise FleetError(
            "Evidence does not finalize every selected checkout exactly once"
        )
    return expected, records, started, completed


def _evidence_composition(
    value: dict[str, Any],
    registry: dict[str, Any],
    identity: str,
    records: list[dict[str, Any]],
    *,
    composed: bool,
) -> set[str]:
    sources = value.get("source_runs")
    if sources is None:
        if value["mode"] == "mixed":
            raise FleetError("Mixed evidence requires source run envelopes")
        return {value["run_id"]}
    if composed or not isinstance(sources, list) or not sources:
        raise FleetError("Invalid or nested evidence composition")
    combined, modes, source_ids = [], set(), set()
    for source in sources:
        combined.extend(validate_evidence(source, registry, identity, composed=True))
        if source["run_id"] in source_ids:
            raise FleetError("Evidence composition contains duplicate source runs")
        source_ids.add(source["run_id"])
        modes.add(source["mode"])
    if combined != records or value["mode"] != (
        next(iter(modes)) if len(modes) == 1 else "mixed"
    ):
        raise FleetError("Evidence composition differs from its source records")
    return source_ids


def _evidence_record_header(
    record: dict[str, Any],
    value: dict[str, Any],
    expected: dict[str, Any],
    source_ids: set[str],
    started: datetime,
    completed: datetime,
) -> str:
    if (
        record.get("run_id") not in source_ids
        or record.get("mode") not in {"product", "baseline"}
        or (value.get("source_runs") is None and record["mode"] != value["mode"])
    ):
        raise FleetError("Record execution run or mode differs from its envelope")
    for name in ("release_manifest_sha256", "collector_id", "execution_policy"):
        if record.get(name) != expected[name]:
            raise FleetError("Evidence record provenance differs from its envelope")
    first, last = (
        _timestamp(record.get("started_at")),
        _timestamp(record.get("completed_at")),
    )
    if not started <= first <= last <= completed:
        raise FleetError("Evidence record timestamps are inconsistent")
    status, errors = record.get("status"), record.get("errors")
    if (
        status not in {"passed", "baseline_passed", "failed", "not_applicable"}
        or not isinstance(errors, list)
        or any(not isinstance(error, str) or not error for error in errors)
    ):
        raise FleetError("Evidence record outcome is invalid")
    if status == "failed" and not errors:
        raise FleetError("Failed evidence lacks its failure reason")
    if status != "failed" and errors:
        raise FleetError("Passing evidence contains failure reasons")
    return status


def _passing_execution_proof(record: dict[str, Any], status: str) -> None:
    if status not in {"passed", "baseline_passed"}:
        return
    _execution_record(record, "governance_")
    if (
        record.get("baseline_status") != "passed"
        or record["governance_exit_code"] != 0
        or record["governance_timed_out"]
        or type(record.get("input_state_stable")) is not bool
        or not record["input_state_stable"]
        or record.get("pre_verification") != record.get("post_verification")
    ):
        raise FleetError("Passing evidence lacks a stable passing baseline")
    if status == "passed":
        _execution_record(record)
        if (
            record["mode"] != "product"
            or record["exit_code"] != 0
            or record["timed_out"]
            or record.get("product_execution") != "requested"
        ):
            raise FleetError("Passing product evidence has inconsistent execution")
    elif (
        record["mode"] != "baseline"
        or record.get("exit_code") is not None
        or record.get("product_execution") != "not_requested"
    ):
        raise FleetError("Baseline-only evidence claims product execution")


def _evidence_execution_proof(record: dict[str, Any], status: str) -> None:
    if "baseline_status" in record:
        _execution_record(record, "governance_")
        _contract_record(record.get("governance_execution"))
        if record.get("make_environment") != {}:
            raise FleetError(
                "Evidence retains unsafe or missing Make environment provenance"
            )
        observed = (
            "passed"
            if record["governance_exit_code"] == 0
            and not record["governance_timed_out"]
            else "failed"
        )
        if record["baseline_status"] != observed:
            raise FleetError("Evidence baseline outcome is inconsistent")
    if record.get("exit_code") is not None:
        _execution_record(record)
        _contract_record(record.get("execution"))
    for name in ("pre_verification", "post_verification"):
        if name in record:
            _snapshot_record(record[name])
    _passing_execution_proof(record, status)


def validate_evidence(
    value: dict[str, Any],
    registry: dict[str, Any],
    identity: str,
    *,
    composed: bool = False,
) -> list[dict[str, Any]]:
    """Validate consistency, never authenticate this same-user execution claim."""
    expected, records, started, completed = _evidence_header(value, registry, identity)
    source_ids = _evidence_composition(
        value, registry, identity, records, composed=composed
    )
    failures = []
    for record in records:
        status = _evidence_record_header(
            record, value, expected, source_ids, started, completed
        )
        _evidence_execution_proof(record, status)
        if status == "failed":
            failures.append(record["id"])
    if (
        value.get("failures") != failures
        or type(value.get("ok")) is not bool
        or value["ok"] != (not failures)
    ):
        raise FleetError("Evidence aggregate outcome differs from its records")
    return records


def merge_evidence(
    registry_path: Path, release_path: Path, sources: list[dict[str, Any]]
) -> dict[str, Any]:
    """Compose finalized source runs unchanged, including honest failed records."""
    registry = load(registry_path)
    _release, identity, _payload = release_files(release_path)
    if not sources:
        raise FleetError("Evidence composition requires finalized source runs")
    records = [
        record
        for source in sources
        for record in validate_evidence(source, registry, identity, composed=True)
    ]
    modes = {source["mode"] for source in sources}
    result = {
        "schema_version": EVIDENCE_SCHEMA,
        "kind": EVIDENCE_KIND,
        "trust": EVIDENCE_TRUST,
        "collector_id": COLLECTOR_ID,
        "execution_policy": EXECUTION_POLICY,
        "release_manifest_sha256": identity,
        "registry_sha256": digest(canonical(registry)),
        "run_id": uuid4().hex,
        "mode": next(iter(modes)) if len(modes) == 1 else "mixed",
        "selection": [record["id"] for record in records],
        "started_at": min((source["started_at"] for source in sources), key=_timestamp),
        "completed_at": datetime.now(UTC).isoformat(),
        "completed": True,
        "checkouts": records,
        "failures": [
            record["id"] for record in records if record["status"] == "failed"
        ],
        "ok": all(source["ok"] for source in sources),
        "source_runs": sources,
    }
    validate_evidence(result, registry, identity)
    return result


def baseline_matches(
    item: dict[str, Any], row: dict[str, Any], record: dict[str, Any], identity: str
) -> bool:
    if item.get("id") != row["id"]:
        return False
    root = Path(row["path"])
    expected = verification_snapshot(record["git"], record["binding_sha256"])
    declared = baseline_command(row, root)
    executed, _environment, contract = execution_contract(
        row["verification"], root, declared
    )
    return (
        item.get("id") == row["id"]
        and item.get("release_manifest_sha256") == identity
        and item.get("baseline_status") == "passed"
        and item.get("governance_exit_code") == 0
        and item.get("governance_timed_out") is False
        and item.get("input_state_stable") is True
        and item.get("pre_verification") == expected
        and item.get("post_verification") == expected
        and item.get("governance_argv") == declared
        and item.get("governance_executed_argv") == executed
        and item.get("governance_execution") == contract
    )


def acceptance_matches(
    item: dict[str, Any], row: dict[str, Any], record: dict[str, Any], identity: str
) -> bool:
    if (
        item.get("id") != row["id"]
        or row["verification"].get("coverage") == "governance-only"
        or item.get("mode") != "product"
    ):
        return False
    executed, _environment, contract = execution_contract(
        row["verification"], Path(row["path"]), row["verification"]["argv"]
    )
    return (
        baseline_matches(item, row, record, identity)
        and item.get("status") == "passed"
        and item.get("errors") == []
        and item.get("exit_code") == 0
        and item.get("timed_out") is False
        and item.get("argv") == row["verification"]["argv"]
        and item.get("executed_argv") == executed
        and item.get("execution") == contract
    )


def observe_governance(record: dict[str, Any], root: Path, row: dict[str, Any]) -> None:
    declared = [str(confined(root, "gov")), "check"]
    actual, environment, contract = execution_contract(
        row["verification"], root, declared
    )
    before = git_identity(root)
    checked = run_digest(
        actual,
        root,
        environment,
        verification_timeout(row["verification"], "structural_timeout_seconds", 120),
    )
    after = git_identity(root)
    record.update(
        structural_argv=declared,
        structural_executed_argv=actual,
        structural_execution=contract,
        structural_check="passed"
        if checked["exit_code"] == 0 and not checked["timed_out"] and before == after
        else "failed",
        **{f"structural_{name}": value for name, value in checked.items()},
    )
    if record["structural_check"] != "passed":
        raise FleetError(
            "Structural governance check failed, timed out or changed inputs"
        )


def _retired_external_observation(root: Path, policy: dict[str, Any]) -> dict[str, Any]:
    lifecycle = policy["lifecycle"]
    if not isinstance(lifecycle, dict) or set(lifecycle) != {
        "state",
        "recorded_at",
        "evidence_ref",
    }:
        raise FleetError("Invalid external absence lifecycle")
    reference = lifecycle["evidence_ref"]
    head = policy.get("observed_git_head")
    if (
        lifecycle["state"] != "retired_absent"
        or not isinstance(reference, str)
        or not re.fullmatch(r"[A-Za-z0-9_./-]+", reference)
        or PurePosixPath(reference).is_absolute()
        or reference != PurePosixPath(reference).as_posix()
        or any(part in {".", ".."} for part in reference.split("/"))
        or policy.get("observed_root_validity") != "resolved"
        or not isinstance(head, str)
        or not re.fullmatch(r"[0-9a-f]{40}|[0-9a-f]{64}", head)
    ):
        raise FleetError("Invalid external absence identity or evidence reference")
    try:
        recorded = datetime.fromisoformat(lifecycle["recorded_at"])
        if recorded.utcoffset() is None or recorded > datetime.now(UTC):
            raise ValueError("invalid observation time")
    except (TypeError, ValueError) as error:
        raise FleetError("Invalid external absence observation time") from error
    try:
        root.lstat()
    except FileNotFoundError:
        return {
            "location_exists": False,
            "observed_at": datetime.now(UTC).isoformat(),
            "policy_release_manifest_sha256": policy["release_manifest_sha256"],
            "lifecycle": lifecycle,
        }
    raise FleetError("Retired external location reappeared; owner review required")


def _observe_external(
    root: Path, record: dict[str, Any], policy: dict[str, Any], identity: str
) -> None:
    if (
        policy.get("mode") != "read-only"
        or not policy.get("owner")
        or not policy.get("reason")
        or policy.get("release_manifest_sha256") != identity
    ):
        raise FleetError(
            "External integration lacks location or declared read-only policy"
        )
    if "lifecycle" in policy:
        record.update(
            external_policy=policy,
            status="governed_absent",
            observation=_retired_external_observation(root, policy),
        )
        return
    if not root.exists():
        raise FleetError(
            "External integration lacks location or declared read-only policy"
        )
    record.update(
        external_policy=policy,
        status="governed_read_only",
        observation={
            "location_exists": True,
            "policy_release_manifest_sha256": identity,
        },
    )


def audit_checkout(
    row: dict[str, Any],
    record: dict[str, Any],
    release: dict[str, Any],
    identity: str,
    payload: dict[str, bytes],
    evidence: list[dict[str, Any]],
    run_checks: bool,
) -> None:
    root = Path(row["path"])
    if row["integration_mode"] == "blocked":
        raise FleetError("Location remains blocked")
    if row["integration_mode"] != "external" and row.get("classification") in PROTECTED:
        raise FleetError(
            "Protected repository classification cannot be directly adopted"
        )
    if row["integration_mode"] == "external":
        _observe_external(root, record, row.get("external_policy", {}), identity)
        return
    git, binding = installed_metadata(row, root, release, identity, payload)
    record.update(installed=True, git=git, binding_sha256=binding)
    if run_checks:
        observe_governance(record, root, row)
    baseline = any(baseline_matches(item, row, record, identity) for item in evidence)
    accepted = any(acceptance_matches(item, row, record, identity) for item in evidence)
    record["governance_baseline"] = "passed" if baseline else "not_observed"
    if accepted:
        record["product_acceptance"] = "collector_reported_passed"
    if row["verification"].get("coverage") == "governance-only":
        record["verification_coverage"] = "governance-only"
        raise FleetError(
            "Product acceptance unavailable: verification coverage is governance-only"
        )
    if row.get("acceptance_required", True) and not accepted:
        raise FleetError("Product acceptance evidence missing or stale")
    if not baseline:
        raise FleetError("Declared governance baseline has not been observed")
    record["status"] = "accepted" if accepted else "installed_only"


def audit(
    registry_path: Path,
    release_path: Path,
    *,
    run_checks: bool = False,
    evidence_path: Path | None = None,
    baseline_evidence_path: Path | None = None,
) -> dict[str, Any]:
    registry = load(registry_path)
    release, identity, payload = release_files(release_path)
    evidence = []
    evidence_errors = []
    for selected_path, baseline_only in (
        (evidence_path, False),
        (baseline_evidence_path, True),
    ):
        if selected_path is None:
            continue
        try:
            value = json.loads(selected_path.read_text())
            selected_records = validate_evidence(value, registry, identity)
            if baseline_only and any(
                record["mode"] != "baseline" for record in selected_records
            ):
                raise FleetError("Baseline evidence source claims product execution")
            evidence.extend(selected_records)
        except (FleetError, OSError, ValueError, TypeError, KeyError) as error:
            evidence_errors.append(str(error))
    result = {
        "schema_version": SCHEMA,
        "kind": "fleet-adoption-audit",
        "release_manifest_sha256": identity,
        "registry_sha256": digest(canonical(registry)),
        "evidence_trust": EVIDENCE_TRUST,
        "completed_at": datetime.now(UTC).isoformat(),
        "product_execution": "not_requested",
        "coverage_complete": bool(
            registry.get("scope", {}).get("whole_computer_coverage")
        ),
        "checkouts": [],
        "failures": [],
    }
    if evidence_errors:
        result["evidence_errors"] = evidence_errors
        result["failures"].append("invalid_verification_evidence")
    for row in registry_rows(registry):
        record = {
            "id": row["id"],
            "path": row["path"],
            "integration_mode": row["integration_mode"],
            "installed": False,
            "structural_check": "not_run",
            "governance_baseline": "not_observed",
            "product_acceptance": "not_run",
            "errors": [],
        }
        result["checkouts"].append(record)
        try:
            audit_checkout(
                row, record, release, identity, payload, evidence, run_checks
            )
        except (
            FleetError,
            OSError,
            KeyError,
            UnicodeError,
            ValueError,
            subprocess.TimeoutExpired,
        ) as error:
            record["errors"].append(str(error))
            record["status"] = "incomplete"
            result["failures"].append(row["id"])
    if not result["coverage_complete"]:
        result["failures"].append("scan_scope_incomplete")
    if registry.get("scope", {}).get("scan_errors") or registry.get("scope", {}).get(
        "unresolved"
    ):
        result["coverage_complete"] = False
        result["failures"].append("unresolved_scan_evidence")
    result["ok"] = not result["failures"]
    return result


def arguments(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("plan", "audit", "verify", "baseline"):
        command = sub.add_parser(name)
        command.add_argument("--registry", required=True, type=Path)
        command.add_argument("--release", required=True, type=Path)
        command.add_argument(
            "--output", type=Path, required=name in {"plan", "verify", "baseline"}
        )
        if name == "audit":
            command.add_argument("--run-checks", action="store_true")
            command.add_argument("--evidence", type=Path)
            command.add_argument("--baseline-evidence", type=Path)
        if name in {"verify", "baseline"}:
            command.add_argument("--checkout", action="append", dest="checkout_ids")
        if name == "verify":
            command.add_argument("--baseline-only", action="store_true")
    command = sub.add_parser("apply")
    command.add_argument("--plan", required=True, type=Path)
    command.add_argument("--receipt", required=True, type=Path)
    command = sub.add_parser("undo")
    command.add_argument("--receipt", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = arguments(argv)
    try:
        if args.command == "plan":
            result = plan(args.registry, args.release)
            failed = any(item["status"] == "blocked" for item in result["checkouts"])
        elif args.command == "apply":
            result = apply(args.plan, args.receipt)
            failed = any(
                item["status"] in {"failed", "blocked"} for item in result["checkouts"]
            )
        elif args.command in {"verify", "baseline"}:
            result = collect_acceptance(
                args.registry,
                args.release,
                checkout_ids=args.checkout_ids,
                output_path=args.output,
                baseline_only=args.command == "baseline"
                or getattr(args, "baseline_only", False),
            )
            failed = not result["ok"]
        elif args.command == "undo":
            result = undo(args.receipt)
            failed = any(
                item["status"] == "undo_refused" for item in result["checkouts"]
            )
        else:
            result = audit(
                args.registry,
                args.release,
                run_checks=args.run_checks,
                evidence_path=args.evidence,
                baseline_evidence_path=args.baseline_evidence,
            )
            failed = not result["ok"]
        if getattr(args, "output", None):
            atomic(args.output, canonical(result))
        if args.command in {"plan", "apply", "undo", "verify", "baseline"}:
            # Plans/receipts contain preserved consumer bytes; keep them out of logs.
            summary = {
                "kind": result.get("kind"),
                "checkouts": [
                    {key: row[key] for key in ("id", "status", "errors") if key in row}
                    for row in result["checkouts"]
                ],
            }
            print(canonical(summary).decode(), end="")
        elif not getattr(args, "output", None):
            print(canonical(result).decode(), end="")
        return int(failed)
    except (FleetError, OSError, KeyError, json.JSONDecodeError) as exc:
        print(json.dumps({"error": str(exc)}), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
