from __future__ import annotations

import argparse
import json
import shlex
import sys
import time
from contextlib import redirect_stdout
from pathlib import Path
from typing import Any

if __package__ in {None, ""}:
    root_hint = Path(__file__).resolve().parents[3]
    sys.path.insert(0, str(root_hint))

from tools.workspace_governance import diagnostics as safe_diagnostics
from tools.workspace_governance import telemetry
from tools.workspace_governance.resolver.core import GovernanceError, GovernanceResolver
from tools.workspace_governance.validators.core import (
    binding_environment,
    doctor,
    interpreter_candidate,
    load_manifest_json,
    pack_hash,
    resolve_execution_mode,
    resolve_interpreter,
)


def _emit(value: object, *, as_json: bool) -> None:
    if as_json:
        print(json.dumps(value, indent=2, sort_keys=True))
        return
    if isinstance(value, list):
        for item in value:
            if isinstance(item, dict):
                print(f"{item['code']}: {item['fact']}")
                print(f"  owner: {item['owner']}")
                print(f"  sources: {', '.join(item['sources'])}")
                print(f"  reason: {item['reason']}")
                print(f"  repair: {item['repair']}")
            else:
                print(item)
        return
    print(json.dumps(value, indent=2, sort_keys=True))


def _sync(root: Path, resolver: GovernanceResolver) -> dict[str, object]:
    workspace, _, _ = resolver.load()
    lock_rel = workspace["governance_pack"]["lock"]
    lock_path = root / lock_rel
    if not lock_path.is_file():
        raise ValueError(f"pack lock does not exist: {lock_rel}")
    existing = load_manifest_json(lock_path)
    pack_rel = existing["pack_path"]
    pack_root = root / pack_rel
    if not pack_root.is_dir():
        raise ValueError(f"locked pack path does not exist: {pack_rel}")
    manifest = load_manifest_json(pack_root / "MANIFEST.json")
    digest, files = pack_hash(pack_root)
    version = (pack_root / manifest["version_file"]).read_text(encoding="utf-8").strip()
    lock = {
        "schema_version": 1,
        "pack_path": pack_rel,
        "pack_name": manifest["name"],
        "version": version,
        "revision": manifest["revision"],
        "content_sha256": digest,
        "files": files,
        "lock_algorithm": manifest["lock_algorithm"],
    }
    lock_path.write_text(
        json.dumps(lock, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return lock


def _env(
    root: Path,
    resolver: GovernanceResolver,
    *,
    as_json: bool,
    allow_missing_interpreter: bool = False,
    observed: dict[str, Any] | None = None,
) -> int:
    """Emit the execution environment for callers that are not Python.

    Shell cannot import the resolver, and `resolve_python_env.sh` prints one
    interpreter and rejects arguments, so a second output would break its
    contract. This is the one place a non-Python caller learns both the
    interpreter and the source binding, keeping that derivation out of every
    script that starts a subprocess.
    """
    _, project, _ = resolver.load()
    local, local_tiers = resolver.local_overlay()
    interpreter, _, issues = resolve_interpreter(
        root, project, local, local_tiers=local_tiers
    )
    execution_mode, _, mode_issues = resolve_execution_mode(project, None, local)
    issues.extend(mode_issues)
    if (
        interpreter is None
        and allow_missing_interpreter
        and execution_mode == "installed"
    ):
        candidate, selection_reason = interpreter_candidate(
            root, project, local, local_tiers=local_tiers
        )
        if selection_reason == "repository-local environment":
            interpreter, issues = candidate, mode_issues
    if issues or interpreter is None:
        if observed is not None:
            observed["diagnostics"] = [item.to_dict() for item in issues]
        with redirect_stdout(sys.stderr):
            _emit([item.to_dict() for item in issues], as_json=as_json)
        return 1
    # The mode travels under the variable the project declares for it, so a
    # child process reads the mode back through the same override it set.
    assignments = {
        "PYTHON_BIN": str(interpreter),
    }
    if project["environment"].get("profile", "python") == "python":
        assignments[project["environment"]["execution_mode"]["override"]] = (
            execution_mode
        )
    assignments.update(binding_environment(root, project, None, local))
    if as_json:
        _emit(assignments, as_json=True)
        return 0
    for name, value in assignments.items():
        print(f"export {name}={shlex.quote(value)}")
    return 0


def _doctor_command(
    root: Path,
    resolver: GovernanceResolver,
    args: argparse.Namespace,
    observed: dict[str, Any],
) -> int:
    _, project, _ = resolver.load()
    local, local_tiers = resolver.local_overlay()
    result, issues = doctor(
        root,
        project,
        local,
        route=args.route,
        local_tiers=local_tiers,
        project_source=resolver.project_source(),
    )
    observed["diagnostics"] = [item.to_dict() for item in issues]
    if args.print_interpreter and not issues:
        print(result["interpreter"])
    elif args.print_interpreter:
        with redirect_stdout(sys.stderr):
            _emit([item.to_dict() for item in issues], as_json=False)
    else:
        payload = {**result, "diagnostics": [item.to_dict() for item in issues]}
        _emit(payload, as_json=args.json)
    return 1 if issues else 0


def _record(
    args: argparse.Namespace,
    argv: list[str],
    root: Path,
    status: int,
    started: float,
    observed: dict[str, Any],
) -> None:
    """Record this execution as a first-hand observation, if enabled."""
    resolver = GovernanceResolver.for_root(root)
    settings = resolver.telemetry_settings()

    def event() -> dict[str, object]:
        if settings.key_directory is None:
            raise telemetry.TelemetryError("telemetry store has no directory")
        key = telemetry.telemetry_key(settings.key_directory)
        workspace = load_manifest_json(resolver.workspace_path)
        try:
            project_id = resolver.load()[1]["project"]["id"]
        except (OSError, ValueError, TypeError, KeyError):
            project_id = None
        return telemetry.build_command_execution(
            root=root,
            workspace_id=workspace["workspace_id"],
            project_id=project_id,
            producer_version=telemetry.producer_version(),
            program="gov",
            arguments=telemetry.recorded_gov_arguments(argv, key),
            key=key,
            exit_status=status,
            duration_ms=int((time.monotonic() - started) * 1000),
            governance=telemetry.governance_enrichment(
                args.command, observed["diagnostics"], key, observed["resolution"]
            ),
        )

    telemetry.record(settings, event)


def main(argv: list[str] | None = None) -> int:
    started = time.monotonic()
    original_argv = list(sys.argv[1:] if argv is None else argv)
    raw_argv = list(original_argv)
    as_json = "--json" in raw_argv
    raw_argv = [item for item in raw_argv if item != "--json"]
    parser = argparse.ArgumentParser(prog="gov")
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--json",
        action="store_true",
        help="emit JSON (accepted before or after the subcommand)",
    )
    sub = parser.add_subparsers(dest="command", required=True)
    check_parser = sub.add_parser("check")
    check_parser.add_argument(
        "--summary",
        action="store_true",
        help="emit bounded rule IDs and repair classes without raw checker output",
    )
    check_parser.add_argument(
        "--generated",
        action="store_true",
        help="also execute generated-output checks in the certified consumer environment",
    )
    sub.add_parser("resolve").add_argument(
        "--mode", choices=("LOCAL_FIX", "FEATURE", "STRUCTURAL_CHANGE")
    )
    explain = sub.add_parser("explain")
    explain.add_argument("query", nargs="?")
    explain.add_argument("--task")
    doctor_parser = sub.add_parser("doctor")
    doctor_parser.add_argument("--print-interpreter", action="store_true")
    doctor_parser.add_argument("--route")
    sub.add_parser("sync")
    env_parser = sub.add_parser("env")
    env_parser.add_argument("--allow-missing-interpreter", action="store_true")
    telemetry_parser = sub.add_parser("telemetry")
    telemetry_parser.add_argument("telemetry_action", choices=("status", "selftest"))
    args = parser.parse_args(raw_argv)
    args.json = as_json or args.json
    root = args.root.resolve()
    observed: dict[str, Any] = {"diagnostics": [], "resolution": None}
    status = _dispatch(args, root, explain, observed)
    if args.command == "telemetry":
        return status
    try:
        _record(args, original_argv, root, status, started, observed)
    except Exception as exc:  # noqa: BLE001 - final telemetry-isolation boundary
        # Recording is observational only. The command status above is final;
        # do not catch BaseException so Ctrl-C and process termination remain
        # visible to their callers.
        print(f"gov: telemetry not recorded ({type(exc).__name__})", file=sys.stderr)
    return status


def _telemetry_command(resolver, root, args):
    settings = resolver.telemetry_settings()
    if args.telemetry_action == "status":
        _emit(telemetry.status(settings), as_json=args.json)
        return 0
    payload, exit_status = telemetry.selftest(settings, root)
    _emit(payload, as_json=args.json)
    return exit_status


def _check_command(resolver, args, observed):
    issues = resolver.check(generated=args.generated)
    observed["diagnostics"] = [item.to_dict() for item in issues]
    if getattr(args, "summary", False):
        allowed = [
            source.partition("#")[0] for item in issues for source in item.sources
        ]
        payload = safe_diagnostics.summarize(
            "workspace_check",
            1 if issues else 0,
            json.dumps([item.to_dict() for item in issues]).encode(),
            allowed_locations=allowed,
        )
        if args.json:
            _emit(payload, as_json=True)
        else:
            print(safe_diagnostics.render(payload))
    else:
        _emit(
            [item.to_dict() for item in issues] if issues else {"status": "ok"},
            as_json=args.json,
        )
    return 1 if issues else 0


def _dispatch(
    args: argparse.Namespace,
    root: Path,
    explain: argparse.ArgumentParser,
    observed: dict[str, Any],
) -> int:
    try:
        resolver = GovernanceResolver.for_root(root)
        if args.command == "telemetry":
            return _telemetry_command(resolver, root, args)
        if args.command == "check":
            return _check_command(resolver, args, observed)
        if args.command == "resolve":
            observed["resolution"] = resolver.resolve(mode=args.mode)
            _emit(observed["resolution"], as_json=args.json)
            return 0
        if args.command == "explain":
            if args.query is None and args.task is None:
                explain.print_usage(sys.stderr)
                return 2
            observed["resolution"] = resolver.explain(args.query, task=args.task)
            _emit(observed["resolution"], as_json=args.json)
            return 0
        if args.command == "sync":
            _emit(_sync(root, resolver), as_json=args.json)
            return 0
        if args.command == "env":
            return _env(
                root,
                resolver,
                as_json=args.json,
                allow_missing_interpreter=args.allow_missing_interpreter,
                observed=observed,
            )
        return _doctor_command(root, resolver, args, observed)
    except (
        GovernanceError,
        OSError,
        ValueError,
        TypeError,
        KeyError,
        StopIteration,
        ImportError,
    ) as exc:
        diagnostics = exc.diagnostics if isinstance(exc, GovernanceError) else []
        observed["diagnostics"] = [item.to_dict() for item in diagnostics]
        if diagnostics:
            _emit([item.to_dict() for item in diagnostics], as_json=args.json)
        else:
            payload = {
                "status": "ENVIRONMENT_INVALID"
                if args.command == "doctor"
                else "invalid",
                "diagnostics": [
                    {
                        "code": "ENVIRONMENT_MISMATCH"
                        if isinstance(exc, ImportError) or args.command == "doctor"
                        else "BROKEN_REFERENCE",
                        "fact": args.command,
                        "owner": "environment"
                        if args.command == "doctor"
                        else "project",
                        "sources": [],
                        "reason": str(exc),
                        "repair": "repair the referenced metadata or dependency",
                    }
                ],
            }
            observed["diagnostics"] = payload["diagnostics"]
            with redirect_stdout(
                sys.stderr if getattr(args, "print_interpreter", False) else sys.stdout
            ):
                _emit(payload, as_json=args.json)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
