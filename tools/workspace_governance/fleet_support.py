"""Guarded local-fleet adoption; registry and receipts are private machine data.

A pinned release authenticates bytes; it does not infer publisher trust or policy
approval. Existing verification and instructions retain their original bodies.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import secrets
import shlex
import stat
import subprocess
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path, PurePosixPath
from typing import Any

SCHEMA = 1
BEGIN = "<!-- BEGIN workspace-governance managed dispatch -->"
END = "<!-- END workspace-governance managed dispatch -->"
VERIFY_BEGIN = "# BEGIN workspace-governance managed verification"
VERIFY_END = "# END workspace-governance managed verification"
PROTECTED = {
    "governance_distribution_source",
    "application_private_git_store",
    "repository_metadata_store",
    "dependency_or_tool_cache",
    "dependency_or_extension_checkout",
    "registered_submodule",
    "experiment_generated_checkout",
    "test_fixture_repository",
    "backup_checkout",
    "package_manager_tool_source",
}
MODES = {"direct", "inherited_checkout", "external", "blocked"}
MAKE_ENVIRONMENT_VARIABLES = ("MAKEFLAGS", "MFLAGS", "MAKEFILES", "GNUMAKEFLAGS")
EXECUTION_POLICY = "workspace-governance-execution-v2"
DEFAULT_TOOL_PATH = ("/usr/bin", "/bin", "/usr/sbin", "/sbin", "/opt/homebrew/bin")
DEFAULT_TOOLS = {
    "sh": "/bin/sh",
    "bash": "/bin/bash",
    "zsh": "/bin/zsh",
    "make": "/usr/bin/make",
    "gmake": "/opt/homebrew/bin/gmake",
}
DEFAULT_INHERITED_NAMES = frozenset(
    {
        "HOME",
        "USER",
        "LOGNAME",
        "LANG",
        "LC_ALL",
        "LC_CTYPE",
        "TERM",
        "TMPDIR",
        "TMP",
        "TEMP",
        "XDG_CONFIG_HOME",
        "XDG_STATE_HOME",
        "XDG_CACHE_HOME",
        "CODEX_MODE",
        "CI",
    }
)
PROHIBITED_ENVIRONMENT_NAMES = frozenset(
    {
        *MAKE_ENVIRONMENT_VARIABLES,
        "PATH",
        "BASH_ENV",
        "ENV",
        "ZDOTDIR",
        "SHELLOPTS",
        "BASHOPTS",
        "CDPATH",
        "IFS",
        "PYTHONPATH",
        "PYTHONHOME",
        "PYTHONINSPECT",
        "LD_PRELOAD",
        "LD_LIBRARY_PATH",
        "DYLD_INSERT_LIBRARIES",
        "DYLD_LIBRARY_PATH",
        "DYLD_FRAMEWORK_PATH",
    }
)


class FleetError(ValueError):
    """An invalid declaration or unsafe write."""


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2) + "\n").encode()


def load(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict) or value.get("schema_version") != SCHEMA:
        raise FleetError(f"{path.name}: expected schema_version {SCHEMA}")
    return value


def relative(value: str) -> str:
    path = PurePosixPath(value)
    if not value or path.is_absolute() or any(p in {"..", ".git"} for p in path.parts):
        raise FleetError(f"Unsafe relative path: {value!r}")
    if "\\" in value or value != path.as_posix() or value == ".":
        raise FleetError(f"Noncanonical relative path: {value!r}")
    return value


def confined(root: Path, name: str) -> Path:
    path = root / relative(name)
    for component in (path, *path.parents):
        if component == root:
            break
        if component.is_symlink():
            raise FleetError(f"Symlink in managed path: {name}")
    if not path.resolve().is_relative_to(root.resolve()):
        raise FleetError(f"Path escapes checkout: {name}")
    return path


def state(path: Path) -> dict[str, Any]:
    if path.is_symlink():
        raise FleetError(f"Managed file is a symlink: {path.name}")
    if not path.exists():
        return {"sha256": None, "mode": None}
    if not path.is_file():
        raise FleetError(f"Managed path is not a regular file: {path.name}")
    return {
        "sha256": digest(path.read_bytes()),
        "mode": stat.S_IMODE(path.stat().st_mode),
    }


@contextmanager
def parent_descriptor(path: Path, *, create: bool = True) -> Iterator[int]:
    """Open every directory without following symlinks, including during writes."""
    descriptor = os.open(path.anchor, os.O_RDONLY | os.O_DIRECTORY)
    try:
        for component in path.absolute().parent.parts[1:]:
            try:
                next_descriptor = os.open(
                    component,
                    os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                    dir_fd=descriptor,
                )
            except FileNotFoundError:
                if not create:
                    raise
                os.mkdir(component, 0o755, dir_fd=descriptor)
                next_descriptor = os.open(
                    component,
                    os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                    dir_fd=descriptor,
                )
            os.close(descriptor)
            descriptor = next_descriptor
        yield descriptor
    finally:
        os.close(descriptor)


def atomic(path: Path, data: bytes, mode: int = 0o600) -> None:
    path = path.absolute()
    temporary = ".fleet-" + secrets.token_hex(12)
    with parent_descriptor(path) as directory:
        descriptor = os.open(
            temporary,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
            0o600,
            dir_fd=directory,
        )
        try:
            with os.fdopen(descriptor, "wb") as stream:
                stream.write(data)
                stream.flush()
                os.fchmod(stream.fileno(), mode)
                os.fsync(stream.fileno())
            os.replace(temporary, path.name, src_dir_fd=directory, dst_dir_fd=directory)
            os.fsync(directory)
        finally:
            try:
                os.unlink(temporary, dir_fd=directory)
            except FileNotFoundError:
                pass


def remove(path: Path) -> None:
    with parent_descriptor(path.absolute()) as directory:
        os.unlink(path.name, dir_fd=directory)


def git_query(root: Path, *args: str, optional: bool = False) -> bytes | None:
    environment = {
        key: value
        for key, value in os.environ.items()
        if key in DEFAULT_INHERITED_NAMES
    }
    environment.update(PATH=":".join(DEFAULT_TOOL_PATH), GIT_PAGER="cat", PAGER="cat")
    result = subprocess.run(
        [
            "/usr/bin/git",
            "-c",
            "core.fsmonitor=false",
            "-c",
            "core.untrackedCache=false",
            "-C",
            str(root),
            *args,
        ],
        capture_output=True,
        check=False,
        timeout=15,
        env=environment,
    )
    if result.returncode:
        if optional:
            return None
        raise FleetError("Unresolved Git checkout")
    return result.stdout


def untracked_entry(root: Path, name: str) -> dict[str, Any]:
    """Hash visible bytes with bounded memory; never follow an untracked symlink."""
    path = root / name.rstrip("/")
    with parent_descriptor(path.absolute(), create=False) as directory:
        info = os.stat(path.name, dir_fd=directory, follow_symlinks=False)
        record = {
            "path": name,
            "type": stat.S_IFMT(info.st_mode),
            "mode": stat.S_IMODE(info.st_mode),
        }
        if stat.S_ISLNK(info.st_mode):
            record["sha256"] = digest(
                os.fsencode(os.readlink(path.name, dir_fd=directory))
            )
        elif stat.S_ISREG(info.st_mode):
            descriptor = os.open(
                path.name, os.O_RDONLY | os.O_NOFOLLOW, dir_fd=directory
            )
            with os.fdopen(descriptor, "rb") as stream:
                content = hashlib.sha256()
                for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                    content.update(chunk)
            record["sha256"] = content.hexdigest()
        elif stat.S_ISDIR(info.st_mode):
            raise FleetError(
                "Git-visible directory identity unsupported; declare a child checkout scope"
            )
        else:
            record["device"] = info.st_rdev
    return record


def git_head(root: Path) -> tuple[str | None, str]:
    head = git_query(root, "rev-parse", "--verify", "HEAD", optional=True)
    if head is not None:
        return head.decode().strip(), "commit"
    symbolic = git_query(root, "symbolic-ref", "--quiet", "HEAD", optional=True)
    if symbolic is None:
        raise FleetError("Unresolved Git HEAD")
    branch = symbolic.decode().strip()
    if (
        git_query(root, "show-ref", "--verify", "--quiet", branch, optional=True)
        is not None
    ):
        raise FleetError("Unresolved existing Git HEAD")
    return None, "unborn"


def git_identity(root: Path) -> dict[str, Any]:
    if root.is_symlink() or not root.is_dir():
        raise FleetError("Checkout root is missing or symlinked")
    top = os.fsdecode(git_query(root, "rev-parse", "--show-toplevel")).strip()
    if Path(top).resolve() != root.resolve():
        raise FleetError("Declared location is not a Git checkout root")
    common = Path(os.fsdecode(git_query(root, "rev-parse", "--git-common-dir")).strip())
    if not common.is_absolute():
        common = root / common
    head, head_state = git_head(root)
    status = git_query(root, "status", "--porcelain", "-z", "--untracked-files=all")
    delta = (
        git_query(
            root,
            "diff",
            "--binary",
            "--no-ext-diff",
            "--no-textconv",
            "--no-color",
            "--full-index",
            "--no-renames",
            "--src-prefix=a/",
            "--dst-prefix=b/",
            "HEAD",
        )
        if head
        else (
            git_query(
                root,
                "diff",
                "--binary",
                "--no-ext-diff",
                "--no-textconv",
                "--no-color",
                "--full-index",
                "--no-renames",
                "--src-prefix=a/",
                "--dst-prefix=b/",
                "--cached",
            )
            + git_query(
                root,
                "diff",
                "--binary",
                "--no-ext-diff",
                "--no-textconv",
                "--no-color",
                "--full-index",
                "--no-renames",
                "--src-prefix=a/",
                "--dst-prefix=b/",
            )
        )
    )
    names = git_query(root, "ls-files", "--others", "--exclude-standard", "-z")
    visible = [
        untracked_entry(root, os.fsdecode(name)) for name in names.split(b"\0") if name
    ]
    snapshot = hashlib.sha256()
    for data in (status, delta, canonical(visible)):
        snapshot.update(len(data).to_bytes(8, "big"))
        snapshot.update(data)
    return {
        "common_dir_sha256": digest(os.fsencode(common.resolve())),
        "head": head,
        "head_state": head_state,
        "state_sha256": snapshot.hexdigest(),
        "dirty_entries": len([entry for entry in status.split(b"\0") if entry]),
        "untracked_entries": len(visible),
    }


def _runtime_release_files(
    source_root: Path, entries: list[dict[str, Any]], expected_sources: set[Path]
) -> dict[str, bytes]:
    payload = {}
    for entry in entries:
        name = relative(entry["path"])
        if name in payload:
            raise FleetError("Duplicate release destination")
        source = confined(source_root, entry["source"])
        data = source.read_bytes()
        if digest(data) != entry["sha256"]:
            raise FleetError("Release byte mismatch")
        if entry.get("mode", 0o644) not in {0o644, 0o755}:
            raise FleetError("Unsupported release mode")
        if state(source)["mode"] != entry.get("mode", 0o644):
            raise FleetError("Release source mode mismatch")
        expected_sources.add(source)
        payload[name] = data
    if not payload:
        raise FleetError("Release has no files")
    return payload


def _component_release_files(
    source_root: Path, name: str, component: dict[str, Any], expected_sources: set[Path]
) -> None:
    if (
        not isinstance(component, dict)
        or type(component.get("schema_version")) is not int
        or component["schema_version"] != 1
    ):
        raise FleetError("Unsupported release component")
    files = component.get("files")
    if (
        not isinstance(files, list)
        or not files
        or digest(canonical(files)) != component.get("identity_sha256")
    ):
        raise FleetError("Release component file identity mismatch")
    destinations = set()
    for entry in files:
        destination = relative(entry["path"])
        if destination in destinations or not entry["source"].startswith(name + "/"):
            raise FleetError("Duplicate or misplaced release component file")
        destinations.add(destination)
        source = confined(source_root, entry["source"])
        if (
            source in expected_sources
            or entry.get("mode", 0o644) not in {0o644, 0o755}
            or state(source)
            != {"sha256": entry["sha256"], "mode": entry.get("mode", 0o644)}
        ):
            raise FleetError("Release component bytes or modes differ")
        expected_sources.add(source)


def _release_source_set(source_root: Path, expected_sources: set[Path]) -> None:
    for base, directories, names in os.walk(source_root, followlinks=False):
        if any((Path(base) / name).is_symlink() for name in directories):
            raise FleetError("Symlink directory in immutable release")
        for name in names:
            if Path(base) / name not in expected_sources:
                raise FleetError("Undeclared file in immutable release")


def release_files(
    path: Path, expected: str | None = None
) -> tuple[dict[str, Any], str, dict[str, bytes]]:
    if path.is_symlink() or not path.is_file():
        raise FleetError("Release manifest must be a regular file")
    release = json.loads(path.read_text(encoding="utf-8"))
    if (
        not isinstance(release, dict)
        or type(release.get("schema_version")) is not int
        or release["schema_version"] not in {1, 2}
    ):
        raise FleetError("Unsupported release schema")
    identity = digest(canonical(release))
    if expected is not None and identity != expected:
        raise FleetError("Release manifest changed since planning")
    for key in ("name", "version", "revision"):
        if not isinstance(release.get(key), str) or not release[key]:
            raise FleetError("Release identity is incomplete")
    root_name = release.get("release_root", ".")
    source_root = (
        path.parent.resolve()
        if root_name == "."
        else confined(path.parent.resolve(), root_name)
    )
    expected_sources = {path.resolve()}
    payload = _runtime_release_files(
        source_root, release.get("files", []), expected_sources
    )
    for name in ("verification", "distribution"):
        if name in release:
            if release["schema_version"] != 2:
                raise FleetError("Release components require schema two")
            _component_release_files(source_root, name, release[name], expected_sources)
    _release_source_set(source_root, expected_sources)
    return release, identity, payload


def registry_rows(registry: dict[str, Any]) -> list[dict[str, Any]]:
    rows = registry.get("checkouts")
    if not isinstance(rows, list):
        raise FleetError("Registry requires explicit checkouts, not a discovery census")
    ids: set[str] = set()
    paths: set[str] = set()
    for row in rows:
        if not isinstance(row.get("id"), str) or not row["id"] or row["id"] in ids:
            raise FleetError("Missing or duplicate checkout id")
        path = str(Path(row.get("path", "")).absolute())
        if not row.get("path") or not Path(row["path"]).is_absolute() or path in paths:
            raise FleetError("Missing or duplicate checkout path")
        ids.add(row["id"])
        paths.add(path)
        if row.get("integration_mode") not in MODES:
            raise FleetError(f"{row['id']}: explicit integration_mode required")
        if not row.get("repository_id"):
            raise FleetError(f"{row['id']}: repository identity required")
    return rows


def managed_block(
    original: str, start: str, end: str, body: str, *, after_shebang: bool = False
) -> str:
    if original.count(start) != original.count(end) or original.count(start) > 1:
        raise FleetError("Malformed or duplicate managed section")
    block = f"{start}\n{body.rstrip()}\n{end}\n"
    if start in original:
        first, end_start = original.index(start), original.index(end)
        if end_start <= first:
            raise FleetError("Managed section ends before it begins")
        for position, marker in ((first, start), (end_start, end)):
            if (position and original[position - 1] != "\n") or original[
                position + len(marker) : position + len(marker) + 1
            ] not in {"", "\n"}:
                raise FleetError("Managed markers must occupy complete lines")
        last = end_start + len(end)
        if original[last : last + 1] == "\n":
            last += 1
        return original[:first] + block + original[last:]
    if after_shebang and original.startswith("#!"):
        first, separator, rest = original.partition("\n")
        return first + separator + block + rest
    return original.rstrip("\n") + ("\n\n" if original else "") + block


def shell_runner(text: str, name: str, governance_runner: str = "gov") -> str:
    if not text.startswith("#!") or not re.search(
        r"(?:sh|bash|zsh)(?:\s|$)", text.splitlines()[0]
    ):
        raise FleetError("Shell runner binding requires a shell shebang")
    upward = "/".join([".."] * len(PurePosixPath(name).parent.parts)) or "."
    hook = (
        'GOVERNANCE_REPO_ROOT="$(CDPATH= cd -- "$(dirname -- "$0")/'
        + upward
        + '" && pwd)"\n'
        '"$GOVERNANCE_REPO_ROOT/'
        + relative(governance_runner)
        + '"'
        + (" check" if governance_runner == "gov" else "")
        + ' || exit "$?"'
    )
    return managed_block(text, VERIFY_BEGIN, VERIFY_END, hook, after_shebang=True)


def _make_declaration_tokens(runner: dict[str, Any], name: str) -> set[str]:
    """Public-safe tokens are reviewed boundaries, never evaluated by planning."""
    values = runner.get(name, [])
    if not isinstance(values, list) or any(
        not isinstance(value, str) for value in values
    ):
        raise FleetError(f"Make {name} requires a list of relative tokens")
    if len(values) != len(set(values)):
        raise FleetError(f"Make {name} has duplicate tokens")
    for value in values:
        symbolic = re.sub(r"\$\([A-Za-z_][A-Za-z0-9_]*\)", "variable", value)
        if (
            not symbolic
            or not re.fullmatch(r"[A-Za-z0-9_.*?/-]+", symbolic)
            or symbolic.startswith("/")
            or ".." in PurePosixPath(symbolic).parts
        ):
            raise FleetError(
                f"Make {name} requires public-safe relative tokens: {value}"
            )
    return set(values)


def _make_build_rules(runner: dict[str, Any]) -> set[str]:
    values = runner.get("build_rules", [])
    if not isinstance(values, list) or any(
        not isinstance(value, str) for value in values
    ):
        raise FleetError(
            "Make build_rules requires a list of relative rule expressions"
        )
    if len(values) != len(set(values)):
        raise FleetError("Make build_rules has duplicate expressions")
    for value in values:
        expression = value
        filtered = re.fullmatch(
            r"\$\(filter ([A-Za-z0-9_.%/-]+(?: [A-Za-z0-9_.%/-]+)*),\$\([A-Za-z_][A-Za-z0-9_]*\)\)",
            value,
        )
        if filtered:
            expression = filtered.group(1)
        else:
            expression = re.sub(r"\$\([A-Za-z_][A-Za-z0-9_]*\)", "variable", value)
        if not expression or any(
            not re.fullmatch(r"[A-Za-z0-9_.%/-]+", token)
            or token.startswith("/")
            or ".." in PurePosixPath(token).parts
            for token in expression.split()
        ):
            raise FleetError(
                f"Make build_rules requires public-safe simple build expressions: {value}"
            )
    return set(values)


def _make_unmanaged_text(text: str) -> str:
    if VERIFY_BEGIN not in text and VERIFY_END not in text:
        return text
    # Reuse marker validation, then discard only the old owned section for
    # discovery. Project-owned additions to the gate remain visible.
    validated = managed_block(text, VERIFY_BEGIN, VERIFY_END, "")
    first = validated.index(VERIFY_BEGIN)
    last = validated.index(VERIFY_END) + len(VERIFY_END)
    if validated[last : last + 1] == "\n":
        last += 1
    return validated[:first] + validated[last:]


def _make_logical_lines(text: str) -> Iterator[str]:
    lines = iter(text.splitlines())
    for raw in lines:
        while raw.endswith("\\"):
            try:
                raw = raw[:-1] + " " + next(lines).lstrip()
            except StopIteration as exc:
                raise FleetError("Make has an unterminated continued line") from exc
        yield raw


def _make_validate_statement(line: str) -> None:
    if re.search(r"\$[({]\s*(?:eval|call|guile|shell|file)\b", line):
        raise FleetError(
            "Make effectful eval/call/guile/shell/file cannot establish a literal verification graph"
        )
    if line.startswith(".IGNORE"):
        raise FleetError("Make .IGNORE can bypass governance rejection")
    if line.startswith(".SECONDEXPANSION"):
        raise FleetError("Make secondary expansion is unsupported for verification")
    if re.match(r"(?:(?:override|export|private)\s+)*\.RECIPEPREFIX\b", line):
        raise FleetError("Make recipe prefix changes are unsupported")
    if re.match(r"(?:ifeq|ifneq|ifdef|ifndef|else|endif|define|endef)\b", line):
        raise FleetError(
            "Make conditional or defined rules are unsupported for verification"
        )


def _make_metadata_line(line: str, includes: set[str], used_includes: set[str]) -> bool:
    included = re.match(r"(-include|sinclude|include)\s+(.+)$", line)
    if included:
        directive, expression = included.groups()
        if directive == "include":
            raise FleetError(
                "Make include rules are unsupported; isolate build metadata"
            )
        if expression not in includes:
            raise FleetError(
                f"Make optional include requires reviewed build_includes: {expression}"
            )
        used_includes.add(expression)
        return True
    if re.match(
        r"(?:(?:override|export|private)\s+)*(?:MAKEFLAGS|MFLAGS|GNUMAKEFLAGS|MAKEFILES|MAKEOVERRIDES|SHELL)\s*(?:[?:+!]?=)",
        line,
    ):
        raise FleetError(
            "Make control variable assignments may bypass governance rejection"
        )
    if re.match(
        r"(?:(?:override|export|private)\s+)*\.[A-Za-z_][A-Za-z0-9_]*\s*(?:[?:+!]?=)",
        line,
    ):
        raise FleetError(
            "Make special variable assignments may change the verification graph"
        )
    # Ordinary assignments are opaque consumer build facts. Target-specific
    # assignments are recorded as rule errors and rejected if reached.
    return bool(
        re.match(
            r"(?:override\s+|export\s+|private\s+)?[A-Za-z0-9_.-]+\s*(?:[?:+!]?=)", line
        )
    )


def _make_rule_targets(
    left: str, build_rules: set[str], used_rules: set[str]
) -> list[str]:
    if "\\" in left:
        raise FleetError("Make escaped target names are unsupported for verification")
    if "$" in left or "%" in left:
        if left not in build_rules:
            raise FleetError(
                f"Make dynamic or pattern rule requires reviewed build_rules: {left}"
            )
        used_rules.add(left)
    return left.split() if "$" not in left else [left]


def _make_rule_error(left: str, right: str, prerequisite: str) -> str | None:
    error = "grouped-target" if left.endswith("&") else None
    if right.startswith(":"):
        error = "double-colon"
    elif ":" in right.split(";", 1)[0]:
        error = "static-pattern or target-specific"
    elif re.search(r"(?:^|\s)[A-Za-z0-9_]+\s*[?+!]?=", right):
        error = "target-specific assignment"
    if prerequisite.count("|") > 1:
        error = "multiple order-only separators"
    return error


def _make_record_rule(
    graph: dict[str, dict[str, Any]],
    line: str,
    build_rules: set[str],
    used_rules: set[str],
) -> list[str]:
    if ":" not in line:
        raise FleetError(
            f"Make unsupported top-level statement may change verification rules: {line}"
        )
    left, right = line.split(":", 1)
    left = left.strip()
    names = _make_rule_targets(left, build_rules, used_rules)
    prerequisite, separator, _recipe = right.partition(";")
    error = _make_rule_error(left, right, prerequisite)
    dependencies = prerequisite.replace("|", " ").split()
    current = []
    for name in names:
        if not name or name == ".PHONY":
            continue
        node = graph.setdefault(
            name, {"dependencies": [], "recipe": False, "error": None}
        )
        node["dependencies"].extend(dependencies)
        node["recipe"] = node["recipe"] or bool(separator)
        node["error"] = node["error"] or error
        current.append(name)
    return current


def _make_require_used(
    declaration: str, declared: set[str], used: set[str], *, kind: str
) -> None:
    unused = declared - used
    if unused:
        raise FleetError(
            f"Make {declaration} contains unused {kind}: " + ", ".join(sorted(unused))
        )


def _make_graph(
    text: str, build: set[str], includes: set[str], build_rules: set[str]
) -> tuple[dict[str, dict[str, Any]], set[str]]:
    graph: dict[str, dict[str, Any]] = {}
    used_includes: set[str] = set()
    used_rules: set[str] = set()
    current: list[str] = []
    for raw in _make_logical_lines(text):
        if raw.startswith("\t"):
            for name in current:
                graph[name]["recipe"] = True
            continue
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        current = []
        _make_validate_statement(line)
        if _make_metadata_line(line, includes, used_includes):
            continue
        current = _make_record_rule(graph, line, build_rules, used_rules)
    _make_require_used("build_rules", build_rules, used_rules, kind="expressions")
    _make_require_used("build_includes", includes, used_includes, kind="expressions")
    return graph, used_includes


def _make_verification_targets(runner: dict[str, Any]) -> list[str]:
    targets = runner.get("targets", [])
    if (
        not isinstance(targets, list)
        or not targets
        or any(
            not isinstance(target, str) or not re.fullmatch(r"[A-Za-z0-9_.-]+", target)
            for target in targets
        )
        or len(targets) != len(set(targets))
    ):
        raise FleetError("Make binding requires explicit unique safe targets")
    return targets


def _make_validate_gate(
    graph: dict[str, dict[str, Any]], targets: list[str], gate: str, text: str
) -> None:
    for target in targets:
        if target not in graph:
            raise FleetError(f"Declared Make target is absent: {target}")
    if gate in graph and VERIFY_BEGIN not in text:
        raise FleetError("Existing workspace-governance-check target is unmanaged")
    if gate in graph and graph[gate]["recipe"]:
        raise FleetError(
            "Project-owned governance gate recipe would replace the managed recipe"
        )
    graph.setdefault(gate, {"dependencies": [], "recipe": False, "error": None})


def _make_literal_node(graph: dict[str, dict[str, Any]], name: str) -> dict[str, Any]:
    if not re.fullmatch(r"[A-Za-z0-9_./-]+", name):
        raise FleetError(
            f"Make dynamic prerequisite requires reviewed build_prerequisites: {name}"
        )
    if name not in graph:
        raise FleetError(
            f"Make verification prerequisite requires an explicit literal rule or build_prerequisites: {name}"
        )
    node = graph[name]
    if node["error"]:
        raise FleetError(
            f"Make {node['error']} rule is unsupported in verification: {name}"
        )
    return node


def _make_closure(
    graph: dict[str, dict[str, Any]],
    roots: list[str],
    build: set[str],
    gate: str,
    used_build: set[str],
    *,
    governance: bool = False,
) -> set[str]:
    visited: set[str] = set()
    active: set[str] = set()

    def visit(name: str) -> None:
        if name in active:
            raise FleetError(f"Make verification dependency cycle at {name}")
        if name in visited:
            return
        if name == gate and not governance:
            return
        if name in build:
            if governance:
                raise FleetError(f"Make build exclusion overlaps gate ancestry: {name}")
            used_build.add(name)
            return
        node = _make_literal_node(graph, name)
        active.add(name)
        for dependency in node["dependencies"]:
            visit(dependency)
        active.remove(name)
        visited.add(name)

    for root in roots:
        visit(root)
    return visited


def make_runner(text: str, runner: dict[str, Any]) -> str:
    targets = _make_verification_targets(runner)
    gate = "workspace-governance-check"
    build = _make_declaration_tokens(runner, "build_prerequisites")
    includes = _make_declaration_tokens(runner, "build_includes")
    build_rules = _make_build_rules(runner)
    if set(targets) & (build | {gate}):
        raise FleetError(
            "Make build boundary or governance gate cannot be a verification root"
        )
    unmanaged = _make_unmanaged_text(text)
    graph, _ = _make_graph(unmanaged, build, includes, build_rules)
    _make_validate_gate(graph, targets, gate, text)
    used_build: set[str] = set()
    governance = _make_closure(graph, [gate], build, gate, used_build, governance=True)
    if set(targets) & governance:
        raise FleetError(
            "Make governance gate ancestry contains a product verification root"
        )
    product = _make_closure(graph, targets, build, gate, used_build)
    overlap = {name for name in product & governance if graph[name]["recipe"]}
    if overlap:
        raise FleetError(
            "Make product recipe overlaps governance gate ancestry: "
            + ", ".join(sorted(overlap))
        )
    guarded = set(targets) | {name for name in product if graph[name]["recipe"]}
    _make_require_used("build_prerequisites", build, used_build, kind="tokens")
    hook = (
        ".PHONY: workspace-governance-check\nworkspace-governance-check:\n\t@./"
        + relative(runner.get("governance_runner", "gov"))
        + (" check" if runner.get("governance_runner", "gov") == "gov" else "")
        + "\n\n"
        + " ".join(sorted(guarded))
        + ": | workspace-governance-check"
    )
    return managed_block(text, VERIFY_BEGIN, VERIFY_END, hook)


def prohibited_environment_name(name: str) -> bool:
    return name in PROHIBITED_ENVIRONMENT_NAMES or name.startswith(
        ("BASH_FUNC_", "GIT_", "LD_", "DYLD_")
    )


def execution_policy(verification: dict[str, Any] | None = None) -> dict[str, Any]:
    """Explicit approved inputs, not a sandbox against an actor with this UID."""
    declaration = (verification or {}).get("execution_environment", {})
    if not isinstance(declaration, dict) or set(declaration) - {
        "version",
        "inherit",
        "path",
        "approved_tools",
    }:
        raise FleetError("Invalid execution environment declaration")
    if declaration.get("version", EXECUTION_POLICY) != EXECUTION_POLICY:
        raise FleetError("Unsupported execution environment policy")
    inherited = declaration.get("inherit", sorted(DEFAULT_INHERITED_NAMES))
    if (
        not isinstance(inherited, list)
        or any(
            not isinstance(name, str)
            or not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", name)
            or prohibited_environment_name(name)
            for name in inherited
        )
        or len(set(inherited)) != len(inherited)
    ):
        raise FleetError("Unsafe execution environment inheritance")
    paths = declaration.get("path", list(DEFAULT_TOOL_PATH))
    if (
        not isinstance(paths, list)
        or not paths
        or any(
            not isinstance(path, str)
            or not Path(path).is_absolute()
            or ":" in path
            or str(Path(path)) != path
            or ".." in Path(path).parts
            for path in paths
        )
    ):
        raise FleetError("Execution PATH requires approved absolute directories")
    tools = dict(DEFAULT_TOOLS)
    declared_tools = declaration.get("approved_tools", {})
    if not isinstance(declared_tools, dict) or set(declared_tools) - set(DEFAULT_TOOLS):
        raise FleetError("Unknown approved shell or Make executable")
    for name, path in declared_tools.items():
        if (
            not isinstance(path, str)
            or not Path(path).is_absolute()
            or str(Path(path)) != path
            or ".." in Path(path).parts
        ):
            raise FleetError("Approved executable requires an absolute path")
        tools[name] = path
    return {
        "version": EXECUTION_POLICY,
        "inherit": sorted(inherited),
        "path": paths,
        "approved_tools": tools,
    }


def approved_executable(
    value: str, verification: dict[str, Any] | None = None
) -> tuple[str, dict[str, Any]]:
    policy = execution_policy(verification)
    name = Path(value).name
    expected = policy["approved_tools"].get(name)
    if expected is None or value not in {name, expected}:
        raise FleetError("Executable differs from the approved execution policy")
    path = Path(expected)
    if not path.is_file() or not os.access(path, os.X_OK):
        raise FleetError("Approved executable is unavailable")
    return expected, {
        "path": expected,
        "resolved_path": str(path.resolve()),
        "sha256": digest(path.read_bytes()),
        "mode": stat.S_IMODE(path.stat().st_mode),
    }


def command_without_env(argv: list[str]) -> list[str]:
    if Path(argv[0]).name != "env":
        return argv
    if argv[0] not in {"env", "/usr/bin/env"}:
        raise FleetError("Foreign env executable is forbidden")
    index = 1
    while index < len(argv) and re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*=.*", argv[index]):
        if prohibited_environment_name(argv[index].split("=", 1)[0]):
            raise FleetError(
                "Environment assignment may replace or bypass the verifier"
            )
        index += 1
    if index < len(argv) and argv[index] == "--":
        index += 1
    if index >= len(argv) or argv[index].startswith("-"):
        raise FleetError("Verification env argv requires an explicit command")
    return argv[index:]


def runner_path_matches(value: str, name: str, root: Path) -> bool:
    path = Path(value)
    if path.is_absolute():
        return path == root / name
    return path.as_posix().removeprefix("./") == name


def shell_command_matches(command: list[str], name: str, root: Path) -> bool:
    if runner_path_matches(command[0], name, root) and (
        Path(command[0]).is_absolute() or "/" in command[0]
    ):
        return True
    if Path(command[0]).name not in {"sh", "bash", "zsh"}:
        return False
    index = 1
    while index < len(command) and command[index] in {"-e", "-u", "-eu", "-ue", "--"}:
        index += 1
    return index < len(command) and runner_path_matches(command[index], name, root)


def make_command_targets(
    command: list[str], runner: dict[str, Any], root: Path
) -> list[str]:
    targets = []
    filenames = {path.name for path in root.iterdir()}
    makefile = next(
        (
            name
            for name in ("GNUmakefile", "makefile", "Makefile")
            if name in filenames and (root / name).is_file()
        ),
        "Makefile",
    )
    index = 1
    explicit_file = False
    while index < len(command):
        argument = command[index]
        if argument in {"-f", "--file"} and index + 1 < len(command):
            if explicit_file:
                raise FleetError("Repeated Makefile arguments are forbidden")
            explicit_file = True
            index += 1
            makefile = command[index]
            if not runner_path_matches(makefile, runner["path"], root):
                raise FleetError("Verification makefile differs from the hooked runner")
        elif (
            argument in {"-j", "--jobs"}
            and index + 1 < len(command)
            and command[index + 1].isdigit()
        ):
            index += 1
        elif argument in {"-s", "--silent", "--no-print-directory"} or re.fullmatch(
            r"(?:-j|--jobs=)[0-9]+", argument
        ):
            pass
        elif argument.startswith("-") or "=" in argument:
            raise FleetError(
                "Verification make flags may not bypass or replace the runner"
            )
        else:
            targets.append(argument)
        index += 1
    if not runner_path_matches(makefile, runner["path"], root):
        raise FleetError("Verification makefile differs from the hooked runner")
    return targets


def validate_verification_argv(verification: dict[str, Any], root: Path) -> None:
    argv = verification.get("argv")
    if (
        not isinstance(argv, list)
        or not argv
        or any(not isinstance(arg, str) or not arg for arg in argv)
    ):
        raise FleetError("An explicit verification argv is required")
    runner = verification.get("runner")
    if not runner:
        raise FleetError("A real verification runner binding is required")
    command = command_without_env(argv)
    if runner["kind"] == "shell":
        matched = shell_command_matches(command, relative(runner["path"]), root)
        if matched:
            resolve_shell_command(command, root, verification)
    elif runner["kind"] == "make" and Path(command[0]).name in {"make", "gmake"}:
        approved_executable(command[0], verification)
        targets = make_command_targets(command, runner, root)
        matched = bool(targets) and set(targets).issubset(runner.get("targets", []))
    else:
        matched = False
    if not matched:
        raise FleetError("Verification argv does not invoke the hooked runner")


def resolve_shell_command(
    command: list[str], root: Path, verification: dict[str, Any]
) -> tuple[list[str], dict[str, Any]]:
    """Execute scripts with an approved shell, bypassing PATH-based shebang lookup."""
    if Path(command[0]).name in {"sh", "bash", "zsh"}:
        executable, identity = approved_executable(command[0], verification)
        index = 1
        while index < len(command) and command[index] in {
            "-e",
            "-u",
            "-eu",
            "-ue",
            "--",
        }:
            index += 1
        if index >= len(command):
            raise FleetError("Shell argv requires the declared script")
        name = (
            str(Path(command[index]).relative_to(root))
            if Path(command[index]).is_absolute()
            else command[index].removeprefix("./")
        )
        script = confined(root, name)
        return [
            executable,
            *command[1:index],
            str(script),
            *command[index + 1 :],
        ], identity
    name = (
        str(Path(command[0]).relative_to(root))
        if Path(command[0]).is_absolute()
        else command[0].removeprefix("./")
    )
    script = confined(root, name)
    text = (
        script.read_text()
        if script.is_file()
        else verification.get("runner", {}).get("seed_content", "")
    )
    if not text.startswith("#!"):
        raise FleetError("Declared runner requires a supported shell shebang")
    parts = shlex.split(text.splitlines()[0][2:])
    if parts and parts[0] == "/usr/bin/env":
        parts = parts[1:]
    if not parts or len(parts) != 1 or Path(parts[0]).name not in {"sh", "bash", "zsh"}:
        raise FleetError("Declared runner requires an approved shell shebang")
    executable, identity = approved_executable(parts[0], verification)
    return [executable, str(script), *command[1:]], identity


def resolved_verification(
    verification: dict[str, Any], root: Path, *, argv: list[str] | None = None
) -> tuple[list[str], dict[str, str], dict[str, Any]]:
    """Return executed argv, explicit safe assignments, and executable identity."""
    declared = verification["argv"] if argv is None else argv
    command = command_without_env(declared)
    assignments = {}
    if Path(declared[0]).name == "env":
        for argument in declared[1:]:
            if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*=.*", argument):
                break
            key, value = argument.split("=", 1)
            assignments[key] = value
    if Path(command[0]).name in {"make", "gmake"}:
        executable, identity = approved_executable(command[0], verification)
        return [executable, *command[1:]], assignments, identity
    executed, identity = resolve_shell_command(command, root, verification)
    return executed, assignments, identity


def instruction_aliases(row: dict[str, Any], root: Path) -> dict[str, Any]:
    declared = row.get("instruction_aliases", {})
    if not isinstance(declared, dict) or set(declared) - {"CLAUDE.md"}:
        raise FleetError("Only the bounded CLAUDE.md instruction alias is supported")
    result = {}
    path = root / "CLAUDE.md"
    if declared:
        target = declared["CLAUDE.md"]
        if (
            target not in {"AGENTS.md", "./AGENTS.md"}
            or not path.is_symlink()
            or os.readlink(path) != target
        ):
            raise FleetError("Declared instruction alias is missing or changed")
        authority = confined(root, "AGENTS.md")
        if not authority.is_file():
            raise FleetError("Instruction alias requires a regular AGENTS.md authority")
        info = path.lstat()
        result["CLAUDE.md"] = {
            "target": target,
            "type": stat.S_IFMT(info.st_mode),
            "mode": stat.S_IMODE(info.st_mode),
        }
    elif path.is_symlink():
        raise FleetError(
            "Declare the bounded CLAUDE.md instruction alias before adoption"
        )
    return result


def verification_text(verification: dict[str, Any], root: Path) -> tuple[str, str, int]:
    runner = verification["runner"]
    name = relative(runner["path"])
    path = confined(root, name)
    seed = runner.get("seed_content")
    if seed is not None and (
        runner["kind"] != "shell" or verification.get("coverage") != "governance-only"
    ):
        raise FleetError("New runners require shell seed and governance-only coverage")
    if path.is_file():
        return name, path.read_text(), state(path)["mode"] or 0o644
    if seed is None:
        raise FleetError(f"Verification runner is absent: {name}")
    if not isinstance(seed, str) or runner.get("seed_mode", 0o755) not in {
        0o644,
        0o755,
    }:
        raise FleetError("Invalid shell runner seed")
    return name, seed, runner.get("seed_mode", 0o755)


def verification_binding(row: dict[str, Any], root: Path) -> tuple[str, bytes, int]:
    verification = row.get("verification", {})
    validate_verification_argv(verification, root)
    name, text, mode = verification_text(verification, root)
    runner = verification["runner"]
    if runner["kind"] == "shell":
        integration = runner.get("integration", "managed-hook")
        if integration == "native":
            if runner.get("native_sha256") != digest(text.encode()):
                raise FleetError("Reviewed native verifier body has drifted")
            if text.count(VERIFY_BEGIN) != 1 or text.count(VERIFY_END) != 1:
                raise FleetError(
                    "Native verifier requires one complete managed declaration"
                )
            # Validate grammar, without replacing the reviewed native graph.
            managed_block(text, VERIFY_BEGIN, VERIFY_END, "native graph declaration")
        elif integration == "managed-hook":
            text = shell_runner(text, name, runner.get("governance_runner", "gov"))
        else:
            raise FleetError("Unsupported governance verification integration")
    elif runner["kind"] == "make":
        text = make_runner(text, runner)
    else:
        raise FleetError("Unsupported verification runner kind")
    return name, text.encode(), mode


def rendered(
    row: dict[str, Any], root: Path, release: dict[str, Any], identity: str
) -> dict[str, tuple[bytes, int]]:
    body = (
        "Use `./gov check` before work and `./gov resolve` for governed context.\n"
        "The adopted pack is pinned in `.governance/fleet.json`. Existing local\n"
        "instructions and stronger verification duties continue to apply.\n"
        "Run the declared project verifier; its managed hook checks the pack first."
    )
    generated: dict[str, tuple[bytes, int]] = {}
    aliases = instruction_aliases(row, root)
    for name in ("AGENTS.md", "CLAUDE.md"):
        if name in aliases:
            continue
        path = confined(root, name)
        text = path.read_text() if path.exists() else ""
        generated[name] = (
            managed_block(text, BEGIN, END, body).encode(),
            state(path)["mode"] or 0o644,
        )
    verification = row.get("verification", {})
    public_verification = {
        key: value
        for key, value in verification.items()
        if key != "execution_environment"
    }
    metadata = {
        "schema_version": SCHEMA,
        "release": {k: release[k] for k in ("name", "version", "revision")},
        "release_manifest_sha256": identity,
        "profile": row.get("profile", "tools"),
        "integration_mode": row["integration_mode"],
        "verification": public_verification,
        "authorities": list(dict.fromkeys(row.get("authorities", []))),
        "instruction_aliases": aliases,
        "skill_adapters": row.get("skill_adapters", []),
    }
    if "execution_environment" in verification:
        metadata["execution_environment_binding"] = {
            "version": execution_policy(verification)["version"],
            "source": "private-registry",
            "sha256": digest(canonical(verification["execution_environment"])),
        }
    generated[".governance/fleet.json"] = (canonical(metadata), 0o644)
    for entry in row.get("managed_files", []):
        name = relative(entry["path"])
        if name in generated:
            raise FleetError(f"Duplicate generated destination: {name}")
        generated[name] = (entry["content"].encode(), entry.get("mode", 0o644))
    name, data, mode = verification_binding(row, root)
    if name in generated:
        raise FleetError(f"Binding and verification both own {name}")
    generated[name] = (data, mode)
    return generated
