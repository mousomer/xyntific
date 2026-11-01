from __future__ import annotations

import fnmatch
import hashlib
import json
import os
import re
import shutil
import stat
import subprocess
import tomllib
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from tools.workspace_governance.resolver.roles import (
    ARTIFACT_ROLES,
    REQUIRED_PACK_ROOTS,
    WORK_TYPES,
    WRITE_COMPATIBILITY,
    RoleDeclarationError,
    declaration_path_defect,
    default_project_manifest,
    fixed_root_roles,
    role_bootstrap_lists,
)

DIAGNOSTIC_CLASSES = {
    "DUPLICATE_AUTHORITY",
    "CONFLICTING_VALUE",
    "AMBIGUOUS_AUTHORITY",
    "BROKEN_REFERENCE",
    "PACK_DRIFT",
    "ENVIRONMENT_MISMATCH",
}
SCHEMA_FILES = {
    "workspace": "schemas/workspace.schema.json",
    "project": "schemas/project.schema.json",
    "local": "schemas/workspace-local.schema.json",
}

OBLIGATION_COMPARATORS = {
    "required_all_of": "set_superset",
    "allowed_any_of": "set_subset",
    "trigger_coverage": "set_superset",
    "forbidden_set": "set_superset",
    "permitted_set": "set_subset",
    "exact_set": "equality",
    "minimum": "numeric_gte",
    "maximum": "numeric_lte",
    "required_boolean": "false_to_true_only",
    "identity": "equality",
    "unknown": "human_review",
}
# The modes the resolver implements. Not a per-project choice, so it is code
# rather than manifest: a project may say which variable carries the mode and
# what it defaults to, never invent a third meaning for it.
EXECUTION_MODES = ("source", "installed")
OVERLAY_REASONS = {
    "local": "approved local overlay",
    "workspace": "inherited workspace overlay",
}
SUPPORTED_COMMANDS = {
    "check",
    "resolve",
    "explain",
    "doctor",
    "sync",
    "env",
    "telemetry",
}
PYTHON_APPLICATION_FIELDS = {
    "python_requires",
    "dependency_authority",
    "execution_mode",
    "package_import",
    "editable_source",
    "critical_packages",
}
# Diagnostics name the manifest the workspace declares; this label stands in
# only when the workspace names no default project at all.
UNRESOLVED_PROJECT_SOURCE = "project manifest"


@dataclass(frozen=True)
class Diagnostic:
    code: str
    fact: str
    owner: str
    sources: tuple[str, ...]
    reason: str
    repair: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _project_source(workspace: dict[str, Any]) -> str:
    return default_project_manifest(workspace) or UNRESOLVED_PROJECT_SOURCE


def _diag(
    code: str,
    fact: str,
    owner: str,
    sources: list[str],
    reason: str,
    repair: str,
) -> Diagnostic:
    return Diagnostic(code, fact, owner, tuple(sources), reason, repair)


def load_manifest_json(path: Path) -> dict[str, Any]:
    def reject(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"duplicate JSON key: {key}")
            result[key] = value
        return result

    value = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=reject)
    if not isinstance(value, dict):
        raise TypeError("top-level value must be an object")
    return value


def _type_matches(value: object, expected: str) -> bool:
    return {
        "object": isinstance(value, dict),
        "array": isinstance(value, list),
        "string": isinstance(value, str),
        "integer": isinstance(value, int) and not isinstance(value, bool),
        "boolean": isinstance(value, bool),
    }.get(expected, True)


def _schema_issues(  # noqa: C901 - bounded recursive schema evaluator
    value: object,
    schema: dict[str, Any],
    *,
    path: str,
    source: str,
) -> list[Diagnostic]:
    issues: list[Diagnostic] = []

    def fail(reason: str) -> None:
        issues.append(
            _diag(
                "CONFLICTING_VALUE",
                path or "manifest",
                "schema",
                [source],
                reason,
                "change the manifest to satisfy its versioned schema",
            )
        )

    expected = schema.get("type")
    if isinstance(expected, str) and not _type_matches(value, expected):
        fail(f"expected {expected}, got {type(value).__name__}")
        return issues
    if "const" in schema and value != schema["const"]:
        fail(f"value must equal {schema['const']!r}")
    if "enum" in schema and value not in schema["enum"]:
        fail(f"value is not one of {schema['enum']!r}")
    if isinstance(value, str) and len(value) < schema.get("minLength", 0):
        fail("string is shorter than minLength")
    if (
        isinstance(value, int)
        and not isinstance(value, bool)
        and "minimum" in schema
        and value < schema["minimum"]
    ):
        fail(f"integer is below minimum {schema['minimum']}")
    if isinstance(value, list):
        if len(value) < schema.get("minItems", 0):
            fail("array has too few items")
        if "maxItems" in schema and len(value) > schema["maxItems"]:
            fail("array has too many items")
        if schema.get("uniqueItems") and len(
            {json.dumps(v, sort_keys=True) for v in value}
        ) != len(value):
            fail("array items must be unique")
        item_schema = schema.get("items")
        if isinstance(item_schema, dict):
            for index, item in enumerate(value):
                issues.extend(
                    _schema_issues(
                        item, item_schema, path=f"{path}[{index}]", source=source
                    )
                )
    if isinstance(value, dict):
        required = schema.get("required", [])
        for key in required:
            if key not in value:
                issues.append(
                    _diag(
                        "BROKEN_REFERENCE",
                        f"{path}.{key}".strip("."),
                        "schema",
                        [source],
                        "required field is missing",
                        "add the field required by the versioned schema",
                    )
                )
        properties = schema.get("properties", {})
        additional = schema.get("additionalProperties", True)
        for key, nested in value.items():
            nested_schema = (
                properties.get(key) if isinstance(properties, dict) else None
            )
            if nested_schema is None and additional is False:
                fail(f"undeclared field is not allowed: {key}")
            elif nested_schema is None and isinstance(additional, dict):
                nested_schema = additional
            if isinstance(nested_schema, dict):
                issues.extend(
                    _schema_issues(
                        nested,
                        nested_schema,
                        path=f"{path}.{key}".strip("."),
                        source=source,
                    )
                )
        if len(value) < schema.get("minProperties", 0):
            fail("object has too few properties")
        if "maxProperties" in schema and len(value) > schema["maxProperties"]:
            fail("reserved object must remain empty")
    return issues


def _schema_field_usage_issues(schema: dict[str, Any], source: str) -> list[Diagnostic]:
    issues: list[Diagnostic] = []
    allowed = {"consumed", "validation-only", "reserved-for-future-use"}
    registry = load_manifest_json(Path(source).parent.parent / "field-usage.json")

    def visit(node: object, path: str) -> None:
        if not isinstance(node, dict):
            return
        properties = node.get("properties")
        if isinstance(properties, dict):
            for key, child in properties.items():
                child_path = f"{path}.{key}".strip(".")
                usage = (
                    child.get("x-governance-use") if isinstance(child, dict) else None
                )
                behavior = registry.get(child.get("x-governance-behavior"), {})
                reserved_empty = usage != "reserved-for-future-use" or (
                    child.get("type") == "object"
                    and child.get("maxProperties") == 0
                    or child.get("type") == "array"
                    and child.get("maxItems") == 0
                )
                if (
                    usage not in allowed
                    or behavior.get("classification") != usage
                    or not behavior.get("handler")
                    or not behavior.get("test")
                    or not reserved_empty
                ):
                    issues.append(
                        _diag(
                            "BROKEN_REFERENCE",
                            child_path,
                            "schema",
                            [source],
                            "declared field lacks a valid use classification",
                            "classify it as consumed, validation-only, or reserved-for-future-use",
                        )
                    )
                visit(child, child_path)
        visit(node.get("items"), f"{path}[]")
        additional = node.get("additionalProperties")
        if isinstance(additional, dict):
            visit(additional, f"{path}.*")

    visit(schema, "")
    return issues


def _load_schemas(
    pack_root: Path,
) -> tuple[dict[str, dict[str, Any]], list[Diagnostic]]:
    schemas: dict[str, dict[str, Any]] = {}
    issues: list[Diagnostic] = []
    for layer, rel in SCHEMA_FILES.items():
        try:
            schema = load_manifest_json(pack_root / rel)
            usage_issues = _schema_field_usage_issues(schema, str(pack_root / rel))
        except (OSError, ValueError, TypeError) as exc:
            issues.append(
                _diag(
                    "BROKEN_REFERENCE",
                    f"schema:{layer}",
                    "shared_pack",
                    [rel],
                    str(exc),
                    "restore the versioned schema",
                )
            )
            continue
        schemas[layer] = schema
        issues.extend(usage_issues)
    return schemas, issues


def _json_pointer(document: object, pointer: str) -> object:
    if not pointer.startswith("/"):
        raise ValueError("JSON pointer must begin with /")
    value = document
    for token in pointer[1:].split("/"):
        if not isinstance(value, dict):
            raise KeyError(token)
        value = value[token.replace("~1", "/").replace("~0", "~")]
    return value


def _split_reference(root: Path, reference: str) -> tuple[Path, str]:
    raw_path, marker, pointer = reference.partition("#")
    path = root / raw_path
    return path, pointer if marker else ""


def _route_facade(project: dict[str, Any]) -> dict[str, Any]:
    authorities = {item["authority_id"]: item for item in project["authorities"]}
    facade: dict[str, Any] = {}
    for route_id, route in project["routes"].items():
        entry: dict[str, Any] = {}
        if route["dispatch_paths"]:
            entry["dispatch_paths"] = route["dispatch_paths"]
        entry["authority_keys"] = sorted(
            authorities[item]["legacy_key"] for item in route["authority_refs"]
        )
        entry["typical_verification_requirements"] = route[
            "typical_verification_requirements"
        ]
        facade[route_id] = entry
    return facade


def _valid_authority_source(
    root: Path, entry: dict[str, Any], by_id: dict[str, Any]
) -> bool:
    source, pointer = _split_reference(root, entry["source"])
    kind = entry["source_type"]
    valid = source.is_dir() if kind == "directory" else source.is_file()
    if kind == "json_pointer" and valid:
        try:
            _json_pointer(load_manifest_json(source), pointer)
            valid = bool(pointer)
        except (OSError, ValueError, KeyError, TypeError):
            valid = False
    elif pointer:
        valid = False
    if kind == "alias":
        target = by_id.get(entry.get("alias_of"))
        valid = bool(
            target
            and target["source_type"] != "alias"
            and target["source"] == entry["source"]
            and not entry["exclusive"]
            and not entry["canonical_governance"]
        )
    elif "alias_of" in entry or entry["authority_type"] == "alias":
        valid = False
    return valid


def _authority_graph_issues(
    root: Path,
    project: dict[str, Any],
    project_source: str = UNRESOLVED_PROJECT_SOURCE,
) -> list[Diagnostic]:
    issues: list[Diagnostic] = []
    authorities = project["authorities"]
    by_id = {item["authority_id"]: item for item in authorities}
    sources: dict[str, list[str]] = {}
    for entry in authorities:
        valid = _valid_authority_source(root, entry, by_id)
        if not valid:
            issues.append(
                _diag(
                    "BROKEN_REFERENCE",
                    "authority:" + entry["authority_id"],
                    "project",
                    [entry["source"]],
                    "invalid typed authority source or alias",
                    "correct the source type; aliases must be nonexclusive references",
                )
            )
        if entry["exclusive"]:
            path, pointer = _split_reference(root, entry["source"])
            identity = path.resolve().as_posix() + ("#" + pointer if pointer else "")
            sources.setdefault(identity, []).append(entry["authority_id"])
    for source, identities in sources.items():
        if len(identities) > 1:
            issues.append(
                _diag(
                    "AMBIGUOUS_AUTHORITY",
                    "authority.source",
                    "project",
                    [source],
                    "exclusive identities share a source: " + ", ".join(identities),
                    "use a nonexclusive alias or distinct typed source",
                )
            )
    return issues + _canonical_owner_issues(root, project, by_id, project_source)


def _canonical_owner_issues(
    root: Path, project: dict[str, Any], by_id: dict[str, Any], project_source: str
) -> list[Diagnostic]:
    issues: list[Diagnostic] = []
    authorities = project["authorities"]
    try:
        ref = project["canonical_owner_set"]
        if "owners" in ref:
            owner_set = ref["owners"]
        else:
            owner_source = by_id[ref["authority_ref"]]["source"]
            owner_set = _json_pointer(
                load_manifest_json(root / owner_source), ref["json_pointer"]
            )
        if (
            not isinstance(owner_set, dict)
            or not owner_set
            or not all(isinstance(v, str) for v in owner_set.values())
        ):
            raise ValueError(
                "canonical owner set must be a nonempty domain-to-source mapping"
            )
        canonical = [a for a in authorities if a["canonical_governance"]]
        for domain, source in owner_set.items():
            matches = [
                a
                for a in canonical
                if a["source"] == source
                and a["authority_type"] == "human"
                and a["exclusive"]
                and a["source_type"] == "file"
            ]
            if len(matches) != 1:
                issues.append(
                    _diag(
                        "BROKEN_REFERENCE",
                        "canonical_owner:" + domain,
                        "project",
                        [source],
                        "canonical owner requires exactly one applicable human authority",
                        "restore the canonical authority identity and flag",
                    )
                )
        if any(a["source"] not in owner_set.values() for a in canonical):
            raise ValueError("canonical authority is outside the defined owner set")
    except (OSError, ValueError, KeyError, TypeError) as exc:
        issues.append(
            _diag(
                "BROKEN_REFERENCE",
                "canonical_owner_set",
                "project",
                [project_source],
                str(exc),
                "restore the existing owner-set reference",
            )
        )
    return issues


def _semantic_manifest_issues(  # noqa: C901 - cross-field integrity transaction
    root: Path, workspace: dict[str, Any], project: dict[str, Any]
) -> list[Diagnostic]:
    source = _project_source(workspace)
    issues = _binding_contract_issues(project, source)
    if issues:
        return issues
    authorities = project["authorities"]
    by_id: dict[str, list[dict[str, Any]]] = {}
    for entry in authorities:
        by_id.setdefault(entry["authority_id"], []).append(entry)
    for authority_id, entries in sorted(by_id.items()):
        if len(entries) > 1:
            issues.append(
                _diag(
                    "DUPLICATE_AUTHORITY",
                    authority_id,
                    "project",
                    [item["source"] for item in entries],
                    "stable authority ID is declared more than once",
                    "deduplicate declarations",
                )
            )
        for entry in entries:
            if not _split_reference(root, entry["source"])[0].exists():
                issues.append(
                    _diag(
                        "BROKEN_REFERENCE",
                        f"authority:{authority_id}",
                        "project",
                        [entry["source"]],
                        "authority source does not exist",
                        "restore or correct the authority source",
                    )
                )
    for required in sorted(set(project.get("required_authorities", [])) - set(by_id)):
        issues.append(
            _diag(
                "BROKEN_REFERENCE",
                f"authority:{required}",
                "project",
                [source],
                "required stable authority ID is missing",
                "restore the stable authority ID",
            )
        )
    scopes: dict[str, list[dict[str, Any]]] = {}
    for entry in authorities:
        if entry["exclusive"]:
            scopes.setdefault(entry["scope"], []).append(entry)
        alias = entry.get("alias_of")
        if alias is not None and alias not in by_id:
            issues.append(
                _diag(
                    "BROKEN_REFERENCE",
                    f"authority:{entry['authority_id']}",
                    "project",
                    [source],
                    f"alias target is missing: {alias}",
                    "restore or correct the alias",
                )
            )
    for scope, entries in scopes.items():
        if len(entries) > 1:
            issues.append(
                _diag(
                    "AMBIGUOUS_AUTHORITY",
                    scope,
                    "project",
                    [item["source"] for item in entries],
                    "exclusive scope has multiple owners",
                    "resolve the authority conflict explicitly",
                )
            )
    reached: set[str] = set()
    requires_facade = any(
        surface["kind"] == "compatibility_facade"
        for surface in project["generated_surfaces"]
    )
    for route_id, route in project["routes"].items():
        for authority_id in route["authority_refs"]:
            if authority_id not in by_id:
                issues.append(
                    _diag(
                        "BROKEN_REFERENCE",
                        f"route:{route_id}",
                        "project",
                        [source],
                        f"route references missing authority {authority_id}",
                        "correct the typed authority reference",
                    )
                )
            reached.add(authority_id)
            if (
                requires_facade
                and authority_id in by_id
                and "legacy_key" not in by_id[authority_id][0]
            ):
                issues.append(
                    _diag(
                        "BROKEN_REFERENCE",
                        f"route:{route_id}",
                        "project",
                        [source],
                        f"routed authority lacks compatibility identity: {authority_id}",
                        "add its legacy_key before exposing it through the route facade",
                    )
                )
        for path in route["dispatch_paths"]:
            if not (root / path).exists():
                issues.append(
                    _diag(
                        "BROKEN_REFERENCE",
                        f"route:{route_id}",
                        "project",
                        [path],
                        "dispatch path does not exist",
                        "restore or correct the dispatch path",
                    )
                )
        for tool_id in route.get("environment_tools", []):
            if tool_id not in project["environment"]["route_tools"]:
                issues.append(
                    _diag(
                        "BROKEN_REFERENCE",
                        f"route:{route_id}",
                        "project",
                        [source],
                        f"route references missing environment tool {tool_id}",
                        "declare the route tool or correct the reference",
                    )
                )
    aliases = {
        item["authority_id"]: item["alias_of"]
        for item in authorities
        if "alias_of" in item
    }
    reached.update(alias for alias, target in aliases.items() if target in reached)
    for entry in authorities:
        if (
            entry["authority_type"] not in {"project_metadata"}
            and entry["authority_id"] not in reached
        ):
            issues.append(
                _diag(
                    "BROKEN_REFERENCE",
                    f"authority:{entry['authority_id']}",
                    "project",
                    [entry["source"]],
                    "authority is unreachable from every project route",
                    "add it to an applicable route or declare a reachable alias",
                )
            )
    for profile_id, profile in project["execution"]["profiles"].items():
        for route in profile["routes"]:
            if route not in project["routes"]:
                issues.append(
                    _diag(
                        "BROKEN_REFERENCE",
                        f"profile:{profile_id}",
                        "project",
                        [source],
                        f"unknown route {route}",
                        "correct the route reference",
                    )
                )
    for scenario in project["execution"]["representative_scenarios"]:
        for route in scenario["routes"]:
            if route not in project["routes"]:
                issues.append(
                    _diag(
                        "BROKEN_REFERENCE",
                        f"scenario:{scenario['id']}",
                        "project",
                        [source],
                        f"unknown route {route}",
                        "correct the route reference",
                    )
                )
    selectable = {
        route
        for profile in project["execution"]["profiles"].values()
        for route in profile["routes"]
    } | {
        route
        for scenario in project["execution"]["representative_scenarios"]
        for route in scenario["routes"]
    }
    for route_id in sorted(set(project["routes"]) - selectable):
        issues.append(
            _diag(
                "BROKEN_REFERENCE",
                f"route:{route_id}",
                "project",
                [source],
                "route is unreachable from every execution profile and scenario",
                "add the route to a profile or declare a representative scenario",
            )
        )
    if project["execution"]["default_mode"] not in project["execution"]["profiles"]:
        issues.append(
            _diag(
                "BROKEN_REFERENCE",
                "execution.default_mode",
                "project",
                [source],
                "default mode has no profile",
                "declare the profile or correct the default",
            )
        )
    for fact in ("canonical", "full"):
        if "command" in project["verification"][fact]:
            continue
        authority_id = project["verification"][fact]["authority_ref"]
        if authority_id not in by_id:
            issues.append(
                _diag(
                    "BROKEN_REFERENCE",
                    f"verification.{fact}",
                    "project",
                    [source],
                    f"verification fact references missing authority {authority_id}",
                    "correct the authority reference",
                )
            )
    for fact in ("python_requires", "dependency_authority"):
        if project["environment"].get("profile", "python") != "python":
            continue
        authority_id = project["environment"][fact]["authority_ref"]
        if authority_id not in by_id:
            issues.append(
                _diag(
                    "BROKEN_REFERENCE",
                    f"environment.{fact}",
                    "project",
                    [source],
                    f"environment fact references missing authority {authority_id}",
                    "correct the authority reference",
                )
            )
    dependency = project["environment"].get("dependency_authority")
    target = by_id.get(dependency["authority_ref"], []) if dependency else []
    if dependency and target and dependency["source"] != target[0]["source"]:
        issues.append(
            _diag(
                "CONFLICTING_VALUE",
                "environment.dependency_authority",
                "project",
                [dependency["source"], target[0]["source"]],
                "dependency source disagrees with authority",
                "use the referenced dependency authority source",
            )
        )
    default_id = workspace["defaults"]["project"]
    members = [item for item in workspace["projects"] if item["id"] == default_id]
    if len(members) != 1:
        issues.append(
            _diag(
                "AMBIGUOUS_AUTHORITY",
                "workspace.default_project",
                "workspace",
                [".governance/workspace.json"],
                "default project must resolve to exactly one member",
                "declare one unique default project",
            )
        )
    return issues


def _declared_path_problem(
    root: Path, path: str, pack_rel: str | None
) -> tuple[str, str] | None:
    defect = declaration_path_defect(path)
    if defect is not None:
        return "CONFLICTING_VALUE", defect
    if pack_rel is not None and (path == pack_rel or path.startswith(pack_rel + "/")):
        return "CONFLICTING_VALUE", "the pack declares the roles of its own files"
    if not (root / path).is_file():
        return "BROKEN_REFERENCE", "declared artifact does not exist"
    return None


def _role_issue(
    source: str, code: str, fact: str, reason: str, repair: str
) -> Diagnostic:
    return _diag(code, fact, "project", [source], reason, repair)


def _role_bootstrap_issues(
    root: Path, lists: dict[str, list[str]], pack_rel: str | None, source: str
) -> list[Diagnostic]:
    issues: list[Diagnostic] = []
    seen: set[str] = set()
    for key, paths in lists.items():
        for index, path in enumerate(paths):
            fact = f"role_bootstrap.{key}[{index}]"
            problem = _declared_path_problem(root, path, pack_rel)
            if problem is not None:
                issues.append(
                    _role_issue(
                        source,
                        problem[0],
                        fact,
                        problem[1],
                        "declare one existing file",
                    )
                )
            elif path in seen:
                issues.append(
                    _role_issue(
                        source,
                        "CONFLICTING_VALUE",
                        fact,
                        "path is declared more than once in role_bootstrap",
                        "give each bootstrap path exactly one role",
                    )
                )
            seen.add(path)
    return issues


def _artifact_role_issues(
    root: Path, workspace: dict[str, Any], project: dict[str, Any], pack_root: Path
) -> list[Diagnostic]:
    """Hold each role declaration to one existing project file.

    The pack owns its own files' roles and a bootstrap root keeps its fixed
    role, so a project declaration can relabel neither to move a write out of
    the Manifest surface.
    """
    try:
        pack_rel: str | None = (
            pack_root.resolve().relative_to(root.resolve()).as_posix()
        )
    except ValueError:
        pack_rel = None
    source = _project_source(workspace)
    lists = role_bootstrap_lists(project)
    issues = _role_bootstrap_issues(root, lists, pack_rel, source)
    if issues:
        return issues
    try:
        fixed = fixed_root_roles(workspace, project, source)
    except RoleDeclarationError as exc:
        return [
            _role_issue(
                source,
                "CONFLICTING_VALUE",
                "role_bootstrap",
                str(exc),
                "declare valid roots",
            )
        ]
    for index, path in enumerate(lists["bookkeeping"]):
        if path in fixed:
            issues.append(
                _role_issue(
                    source,
                    "CONFLICTING_VALUE",
                    f"role_bootstrap.bookkeeping[{index}]",
                    f"a bootstrap root keeps its fixed role {fixed[path][0]}",
                    "remove the root from the bookkeeping list",
                )
            )
    for path, role in project.get("artifact_roles", {}).items():
        fact = f"artifact_roles.{path}"
        problem = _declared_path_problem(root, path, pack_rel)
        if problem is not None:
            issues.append(
                _role_issue(
                    source,
                    problem[0],
                    fact,
                    problem[1],
                    "declare one existing project file",
                )
            )
        elif path in fixed and fixed[path][0] != role:
            issues.append(
                _role_issue(
                    source,
                    "CONFLICTING_VALUE",
                    fact,
                    f"a bootstrap root keeps its fixed role {fixed[path][0]}",
                    "declare the root's fixed role or remove the declaration",
                )
            )
    return issues


def _normalized_facade(value: object) -> dict[str, Any] | None:
    if not isinstance(value, dict):
        return None
    result = {}
    for key, route in value.items():
        if not isinstance(route, dict):
            return None
        keys = route.get("authority_keys")
        if not isinstance(keys, list) or not all(
            isinstance(item, str) for item in keys
        ):
            return None
        result[key] = {**route, "authority_keys": sorted(keys)}
    return result


def _binding_contract_issues(project: dict[str, Any], source: str) -> list[Diagnostic]:
    """Keep optional representations exclusive without weakening old profiles."""
    defects: list[tuple[str, str]] = []
    for fact in ("canonical", "full"):
        spec = project["verification"][fact]
        if set(spec) not in ({"command"}, {"authority_ref", "json_pointer"}):
            defects.append(
                (
                    f"verification.{fact}",
                    "declare either a direct command or one authority/json-pointer reference",
                )
            )
    if set(project["canonical_owner_set"]) not in (
        {"owners"},
        {"authority_ref", "json_pointer"},
    ):
        defects.append(
            (
                "canonical_owner_set",
                "declare either direct owners or one authority/json-pointer reference",
            )
        )
    environment = project["environment"]
    present = PYTHON_APPLICATION_FIELDS & environment.keys()
    if environment.get("profile", "python") == "python":
        if present != PYTHON_APPLICATION_FIELDS:
            defects.append(
                (
                    "environment",
                    "Python application profile requires all application certification fields",
                )
            )
    elif present:
        defects.append(
            (
                "environment",
                "tools profile cannot declare Python application certification fields",
            )
        )
    defects.extend(_surface_contract_defects(project))
    return [
        _diag(
            "CONFLICTING_VALUE",
            fact,
            "project",
            [source],
            reason,
            "declare one complete supported binding",
        )
        for fact, reason in defects
    ]


def _surface_contract_defects(project: dict[str, Any]) -> list[tuple[str, str]]:
    defects: list[tuple[str, str]] = []
    for surface in project["generated_surfaces"]:
        if surface["kind"] == "compatibility_facade":
            if (
                "transform" not in surface
                or "check_command" in surface
                or "timeout_seconds" in surface
            ):
                defects.append(
                    (
                        surface["surface_id"],
                        "compatibility facade requires its transform only",
                    )
                )
        elif "check_command" not in surface or "transform" in surface:
            defects.append(
                (
                    surface["surface_id"],
                    "generated output requires a check command instead of a transform",
                )
            )
    return defects


def _protected_check_paths(root: Path, project: dict[str, Any]) -> set[str]:
    """Snapshot versioned/nonignored files and explicit governance inputs.

    This detects mutation of protected files; it is not an OS sandbox and does
    not certify ignored runtime/cache files or Git administrative metadata.
    """
    try:
        tracked = subprocess.run(
            ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"],
            cwd=root,
            capture_output=True,
            check=False,
        )
    except OSError:
        tracked = None
    if tracked is not None and tracked.returncode == 0:
        paths = {os.fsdecode(path) for path in tracked.stdout.split(b"\0") if path}
    else:
        paths = set()
        excluded = {
            ".git",
            "__pycache__",
            ".pytest_cache",
            ".ruff_cache",
            ".venv",
            "venv",
            "env",
        }
        for directory, directories, files in os.walk(root, followlinks=False):
            paths.update(
                (Path(directory) / name).relative_to(root).as_posix()
                for name in directories
                if (Path(directory) / name).is_symlink()
            )
            directories[:] = [name for name in directories if name not in excluded]
            paths.update(
                (Path(directory) / name).relative_to(root).as_posix() for name in files
            )
    paths = {
        path
        for path in paths
        if "__pycache__" not in Path(path).parts and not path.endswith(".pyc")
    }
    for authority in project["authorities"]:
        relative = authority["source"].partition("#")[0]
        path = root / relative
        paths.add(relative)
        if path.is_dir():
            paths.update(
                item.relative_to(root).as_posix()
                for item in path.rglob("*")
                if item.is_file() or item.is_symlink()
            )
    for surface in project["generated_surfaces"]:
        paths.update(surface[key].partition("#")[0] for key in ("source", "target"))
    workspace = load_manifest_json(root / ".governance/workspace.json")
    lock = load_manifest_json(root / workspace["governance_pack"]["lock"])
    paths.update(
        {
            ".governance/workspace.json",
            workspace["governance_pack"]["lock"],
            _project_source(workspace),
        }
    )
    paths.update(lock["pack_path"] + "/" + relative for relative in lock["files"])
    return paths


def _protected_check_snapshot(
    root: Path, project: dict[str, Any]
) -> dict[str, tuple[Any, ...]]:
    snapshot: dict[str, tuple[Any, ...]] = {}
    for relative in _protected_check_paths(root, project):
        path = root / relative
        try:
            info = path.lstat()
        except FileNotFoundError:
            snapshot[relative] = ("missing",)
            continue
        metadata = (
            info.st_mode,
            info.st_size,
            info.st_mtime_ns,
            info.st_ctime_ns,
            info.st_ino,
        )
        if stat.S_ISLNK(info.st_mode):
            if path.is_file():
                target_info = path.stat()
                with path.open("rb") as stream:
                    target_digest = hashlib.file_digest(stream, "sha256").hexdigest()
                snapshot[relative] = (
                    *metadata,
                    os.readlink(path),
                    target_info.st_mtime_ns,
                    target_info.st_ctime_ns,
                    target_info.st_ino,
                    target_digest,
                )
            else:
                snapshot[relative] = (*metadata, os.readlink(path))
        elif stat.S_ISREG(info.st_mode):
            with path.open("rb") as stream:
                snapshot[relative] = (
                    *metadata,
                    hashlib.file_digest(stream, "sha256").hexdigest(),
                )
        else:
            snapshot[relative] = metadata
    return snapshot


def _regular_generated_paths(root: Path, paths: tuple[str, ...]) -> bool:
    return all(
        (root / relative).is_file()
        and not any(
            (root / part).is_symlink()
            for part in (Path(relative), *Path(relative).parents)
        )
        for relative in paths
    )


def _generated_output_defect(root: Path, surface: dict[str, Any]) -> str | None:
    source, target = surface["source"], surface["target"]
    command = surface["check_command"]
    if any(declaration_path_defect(path) for path in (source, target)):
        return (
            "generated output source and target must be repository-relative file paths"
        )
    elif not all((root / path).is_file() for path in (source, target)):
        return "generated output references a missing renderer or target"
    elif not _regular_generated_paths(root, (source, target)):
        return "generated output renderer and target must be regular files without symlink path components"
    elif not all(
        (root / path).resolve().is_relative_to(root.resolve())
        for path in (source, target)
    ):
        return "generated output paths must remain inside the checkout"
    elif not (
        command[0] == source
        or len(command) > 1
        and command[0] == "{python}"
        and command[1] == source
    ):
        return (
            "check command must invoke its declared renderer directly or with {python}"
        )
    return None


def _generated_mutation_reason(
    root: Path,
    project: dict[str, Any],
    target: str,
    before: dict[str, tuple[Any, ...]],
) -> str | None:
    try:
        after = _protected_check_snapshot(root, project)
        changed = sorted(
            path
            for path in before.keys() | after.keys()
            if before.get(path) != after.get(path)
        )
        confined = _regular_generated_paths(root, (target,)) and (
            root / target
        ).resolve().is_relative_to(root.resolve())
        if changed or not confined:
            digest = hashlib.sha256(
                b"\0".join(os.fsencode(path) for path in changed)
            ).hexdigest()
            return f"generated-output check modified protected inputs or outputs (changed_paths={len(changed)}; paths_sha256={digest}); target remains confined={confined}"
    except (OSError, ValueError, KeyError, TypeError):
        return "generated-output check invalidated protected input or output metadata"
    return None


def _generated_contract(root, project, surface, local, local_tiers):
    command = surface["check_command"]
    if command[0] == "{python}":
        interpreter, _, issues = resolve_interpreter(
            root, project, local, local_tiers=local_tiers
        )
        if issues or interpreter is None:
            return None, issues
        argv = [str(interpreter), *command[1:]]
    else:
        argv = [str(root / command[0]), *command[1:]]
    _, _, issues = resolve_execution_mode(project, local=local)
    if issues:
        return None, issues
    environment = {
        **os.environ,
        **binding_environment(root, project, local=local),
        "PYTHONDONTWRITEBYTECODE": "1",
    }
    contract = (
        tuple(argv),
        str(root.resolve()),
        tuple(sorted(environment.items())),
        surface.get("timeout_seconds", 120),
    )
    return contract, []


def _generated_output_group(root, project, surfaces, contract):
    """Run one exact contract, preserving protection and attribution for all outputs."""
    argv, cwd, environment, timeout = contract
    before = _protected_check_snapshot(root, project)
    reason = None
    try:
        result = subprocess.run(
            list(argv),
            cwd=Path(cwd),
            capture_output=True,
            timeout=timeout,
            check=False,
            env=dict(environment),
        )
        if result.returncode:
            outputs = "; ".join(
                f"{name}_bytes={len(data)}; {name}_sha256={hashlib.sha256(data).hexdigest()}"
                for name, data in (("stdout", result.stdout), ("stderr", result.stderr))
            )
            reason = f"generated output check failed with exit {result.returncode}; {outputs}"
    except subprocess.TimeoutExpired:
        reason = "generated output check exceeded its declared timeout"
    except OSError as exc:
        reason = f"generated output check could not execute (errno={exc.errno})"
    mutation = next(
        (
            issue
            for surface in surfaces
            if (
                issue := _generated_mutation_reason(
                    root, project, surface["target"], before
                )
            )
        ),
        None,
    )
    reason = mutation or reason
    if reason is None:
        return []
    return [
        _diag(
            "CONFLICTING_VALUE",
            surface["surface_id"],
            "project",
            [surface["source"], surface["target"]],
            reason,
            "regenerate through the declared renderer and rerun its read-only check",
        )
        for surface in surfaces
    ]


def _generated_output_issue(
    root, project, surface, local=None, *, local_tiers=None, generated=True
):
    reason = _generated_output_defect(root, surface)
    if reason:
        return [
            _diag(
                "BROKEN_REFERENCE",
                surface["surface_id"],
                "project",
                [surface["source"], surface["target"]],
                reason,
                "declare existing checkout-local renderer and output paths",
            )
        ]
    if not generated:
        return []
    contract, issues = _generated_contract(root, project, surface, local, local_tiers)
    return issues or _generated_output_group(root, project, [surface], contract)


def _generated_surface_issues(  # noqa: C901 - declaration preflight, exact grouping and facade transaction
    root: Path,
    project: dict[str, Any],
    local: dict[str, Any] | None = None,
    *,
    local_tiers: dict[str, str] | None = None,
    generated: bool = False,
) -> list[Diagnostic]:
    issues: list[Diagnostic] = []
    outputs = [
        surface
        for surface in project["generated_surfaces"]
        if surface["kind"] == "generated_output"
    ]
    # Validate every declaration before any renderer can run.
    for surface in outputs:
        issues.extend(
            _generated_output_issue(root, project, surface, local, generated=False)
        )
    if generated and not issues:
        groups = {}
        for surface in outputs:
            contract, defects = _generated_contract(
                root, project, surface, local, local_tiers
            )
            issues.extend(defects)
            if contract is not None:
                groups.setdefault(contract, []).append(surface)
        if not issues:
            for contract, surfaces in groups.items():
                issues.extend(
                    _generated_output_group(root, project, surfaces, contract)
                )
    for surface in project["generated_surfaces"]:
        if surface["kind"] == "generated_output":
            continue
        source_path, source_pointer = _split_reference(root, surface["source"])
        target_path, target_pointer = _split_reference(root, surface["target"])
        if not source_path.is_file() or not target_path.is_file():
            issues.append(
                _diag(
                    "BROKEN_REFERENCE",
                    surface["surface_id"],
                    "project",
                    [surface["source"], surface["target"]],
                    "generated relationship references a missing file",
                    "restore the declared source and target",
                )
            )
            continue
        try:
            source_value = _json_pointer(
                load_manifest_json(source_path), source_pointer
            )
            target_value = _json_pointer(
                load_manifest_json(target_path), target_pointer
            )
        except (ValueError, KeyError, TypeError) as exc:
            issues.append(
                _diag(
                    "BROKEN_REFERENCE",
                    surface["surface_id"],
                    "project",
                    [surface["source"], surface["target"]],
                    f"generated relationship pointer failed: {exc}",
                    "correct the source or target pointer",
                )
            )
            continue
        if surface["transform"] == "authority-legacy-key-route-v1":
            expected = _route_facade(project)
            normalized = _normalized_facade(target_value)
            if source_value != project["routes"] or normalized != expected:
                issues.append(
                    _diag(
                        "CONFLICTING_VALUE",
                        surface["surface_id"],
                        "project",
                        [surface["source"], surface["target"]],
                        "legacy route compatibility facade diverges from canonical project routes",
                        "regenerate or update the compatibility facade from project routes",
                    )
                )
    return issues


def _skill_adapter_issues(root: Path, project: dict[str, Any]) -> list[Diagnostic]:  # noqa: C901 - bounded adapter identity transaction
    """Tie discoverable adapters to canonical bytes and the adopted release."""
    issues = []
    seen = set()
    for entry in project.get("skill_adapters", []):
        path, target = entry["path"], entry["target"]
        reason = None
        if path in seen:
            reason = "duplicate skill adapter path"
        seen.add(path)
        if any(declaration_path_defect(item) for item in (path, target)):
            reason = "skill adapter paths must be confined repository-relative files"
        elif not _regular_generated_paths(root, (path, target)) or not all(
            (root / item).resolve().is_relative_to(root.resolve())
            for item in (path, target)
        ):
            reason = "skill adapter and canonical target must be confined regular files"
        else:
            expected_target = (
                "tools/workspace_governance/skills/" + entry["name"] + "/SKILL.md"
            )
            text = (root / path).read_text(encoding="utf-8")
            canonical = (root / target).read_text(encoding="utf-8")
            canonical_names = re.findall(r"^name: ([^\n]+)$", canonical, re.MULTILINE)
            descriptions = re.findall(
                r"^description: ([^\n]+)$", canonical, re.MULTILINE
            )
            names = re.findall(r"^name: ([^\n]+)$", text, re.MULTILINE)
            links = re.findall(r"\[canonical shared skill\]\(([^)]+)\)", text)
            pins = re.findall(r"Canonical source SHA-256: `([0-9a-f]{64})`", text)
            binding = root / ".governance/fleet.json"
            if not re.fullmatch(r"[0-9a-f]{64}", entry["sha256"]) or not re.fullmatch(
                r"[0-9a-f]{64}", entry["release_manifest_sha256"]
            ):
                reason = (
                    "skill adapter digests must be full lowercase SHA-256 identities"
                )
            elif (
                target != expected_target
                or names != [entry["name"]]
                or len(links) != 1
                or canonical_names != [entry["name"]]
                or len(descriptions) != 1
            ):
                reason = (
                    "skill adapter name or canonical target disagrees with declaration"
                )
            elif (
                Path(links[0]).is_absolute()
                or (root / path).parent.joinpath(links[0]).resolve()
                != (root / target).resolve()
            ):
                reason = "skill adapter discovery link disagrees with canonical target"
            elif text != (
                "---\nname: "
                + entry["name"]
                + "\ndescription: "
                + descriptions[0]
                + "\n---\n\nRead and apply the [canonical shared skill]("
                + links[0]
                + ") before performing this workflow.\n\n"
                + "Canonical source SHA-256: `"
                + entry["sha256"]
                + "`. This discovery adapter does not own policy.\n"
            ):
                reason = (
                    "skill adapter body diverges from its bounded discovery contract"
                )
            elif (
                pins != [entry["sha256"]]
                or hashlib.sha256((root / target).read_bytes()).hexdigest()
                != entry["sha256"]
            ):
                reason = "skill adapter canonical source hash is stale"
            elif (
                binding.is_symlink()
                or not _regular_generated_paths(root, (".governance/fleet.json",))
                or load_manifest_json(binding).get("release_manifest_sha256")
                != entry["release_manifest_sha256"]
            ):
                reason = (
                    "skill adapter release differs from the adopted release binding"
                )
        if reason:
            issues.append(
                _diag(
                    "BROKEN_REFERENCE",
                    "skill_adapter:" + path,
                    "project",
                    [path, target],
                    reason,
                    "regenerate the adapter contract from the reviewed adopted release",
                )
            )
    return issues


def validate_manifests(
    root: Path,
    workspace: dict[str, Any],
    project: dict[str, Any],
    local: dict[str, Any] | None,
    *,
    pack_root: Path,
    local_documents: list[tuple[dict[str, Any], str]] | None = None,
    local_tiers: dict[str, str] | None = None,
    generated: bool = False,
) -> list[Diagnostic]:
    schemas, issues = _load_schemas(pack_root)
    source = _project_source(workspace)
    if local_documents is None:
        local_documents = (
            [(local, ".governance/workspace.local.json")] if local is not None else []
        )
    layers: list[tuple[str, dict[str, Any] | None, str]] = [
        ("workspace", workspace, ".governance/workspace.json"),
        ("project", project, source),
        *(("local", payload, source) for payload, source in local_documents),
    ]
    for layer, value, source in layers:
        if value is not None and layer in schemas:
            issues.extend(_schema_issues(value, schemas[layer], path="", source=source))
    if issues:
        return issues
    issues.extend(_semantic_manifest_issues(root, workspace, project))
    issues.extend(_skill_adapter_issues(root, project))
    issues.extend(_artifact_role_issues(root, workspace, project, pack_root))
    issues.extend(_authority_graph_issues(root, project, source))
    if not issues:
        issues.extend(
            _generated_surface_issues(
                root, project, local, local_tiers=local_tiers, generated=generated
            )
        )
    for rel in (
        project["sanitation"]["secret_scanner"],
        project["sanitation"]["repository_check"],
    ):
        if not (root / rel).is_file():
            issues.append(
                _diag(
                    "BROKEN_REFERENCE",
                    "sanitation",
                    "project",
                    [rel],
                    "declared sanitation entrypoint does not exist",
                    "restore or correct the entrypoint",
                )
            )
    machine_path = re.compile(r'(?:^|[\s"\'])(?:/Users/|/home/|[A-Za-z]:[\\/])')
    for rel, payload in (
        (".governance/workspace.json", workspace),
        (source, project),
    ):
        if machine_path.search(json.dumps(payload)):
            issues.append(
                _diag(
                    "CONFLICTING_VALUE",
                    "machine_local_path",
                    "local",
                    [rel],
                    "tracked manifest contains a machine-local absolute path",
                    "move the path to the ignored local overlay",
                )
            )
    return issues


def _pack_files(pack_root: Path, exclusions: list[str]) -> list[str]:
    files: list[str] = []
    for path in pack_root.rglob("*"):
        if not path.is_file():
            continue
        rel = path.relative_to(pack_root).as_posix()
        if any(
            fnmatch.fnmatch(rel, pattern) or fnmatch.fnmatch(f"/{rel}", pattern)
            for pattern in exclusions
        ):
            continue
        files.append(rel)
    return sorted(files)


def content_hash(
    root: Path, files: list[str], algorithm: str = "sha256-path-and-content-v1"
) -> str:
    if algorithm != "sha256-path-and-content-v1":
        raise ValueError(f"unsupported pack hash algorithm: {algorithm}")
    digest = hashlib.sha256()
    for rel in sorted(files):
        digest.update(rel.encode() + b"\0")
        digest.update((root / rel).read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


def pack_hash(pack_root: Path) -> tuple[str, list[str]]:
    manifest = load_manifest_json(pack_root / "MANIFEST.json")
    files = _pack_files(pack_root, manifest["pack_hash_excludes"])
    return content_hash(pack_root, files, manifest["lock_algorithm"]), files


def _pack_metadata_issues(
    manifest: dict[str, Any], files: list[str]
) -> list[Diagnostic]:
    """Validate the pack's own role, root, and comparison metadata.

    These declarations make the shared vocabulary, the pack's own artifact
    identities, and the tables enforcement reads explicit. They do not classify
    project artifacts; that is `resolver.roles` over project declarations.
    """
    issues: list[Diagnostic] = []
    issues.extend(_runtime_metadata_issues(manifest))
    work_types = manifest.get("work_types")
    if not isinstance(work_types, list) or sorted(work_types) != sorted(WORK_TYPES):
        issues.append(
            _diag(
                "PACK_DRIFT",
                "MANIFEST.work_types",
                "shared_pack",
                ["MANIFEST.json"],
                "pack must declare exactly Planning, Coding, and Manifest",
                "restore the versioned work-type vocabulary",
            )
        )
    declarations = manifest.get("artifact_roles")
    if not isinstance(declarations, list):
        declarations = []
    declared_paths: list[str] = []
    invalid_roles = False
    for declaration in declarations:
        if not isinstance(declaration, dict):
            invalid_roles = True
            continue
        path, role = declaration.get("path"), declaration.get("role")
        if not isinstance(path, str) or not isinstance(role, str):
            invalid_roles = True
            continue
        declared_paths.append(path)
        invalid_roles |= role not in ARTIFACT_ROLES
    if (
        invalid_roles
        or set(declared_paths) != set(files)
        or len(declared_paths) != len(set(declared_paths))
    ):
        issues.append(
            _diag(
                "PACK_DRIFT",
                "MANIFEST.artifact_roles",
                "shared_pack",
                ["MANIFEST.json"],
                "pack artifact-role declarations must assign one known role to every packed file",
                "declare each packed artifact exactly once with a supported role",
            )
        )
    if manifest.get("obligation_comparators") != OBLIGATION_COMPARATORS:
        issues.append(
            _diag(
                "PACK_DRIFT",
                "MANIFEST.obligation_comparators",
                "shared_pack",
                ["MANIFEST.json"],
                "pack comparator metadata must retain every directional descriptor mapping",
                "restore the versioned obligation-comparator metadata",
            )
        )
    if manifest.get("write_compatibility") != WRITE_COMPATIBILITY:
        issues.append(
            _diag(
                "PACK_DRIFT",
                "MANIFEST.write_compatibility",
                "shared_pack",
                ["MANIFEST.json"],
                "pack write-compatibility metadata must match the table enforcement reads",
                "restore the versioned write-compatibility table",
            )
        )
    roots = manifest.get("bootstrap_roots")
    if (
        not isinstance(roots, list)
        or not all(isinstance(root, str) for root in roots)
        or len(roots) != len(set(roots))
        or not set(roots) <= set(files)
        or not REQUIRED_PACK_ROOTS <= set(roots)
    ):
        issues.append(
            _diag(
                "PACK_DRIFT",
                "MANIFEST.bootstrap_roots",
                "shared_pack",
                ["MANIFEST.json"],
                "pack bootstrap roots must be packed files that keep the enforcement roots",
                "declare each root once and retain "
                + ", ".join(sorted(REQUIRED_PACK_ROOTS)),
            )
        )
    return issues


def _runtime_metadata_issues(manifest: dict[str, Any]) -> list[Diagnostic]:
    runtime = manifest.get("runtime")
    requirement = runtime.get("python_requires") if isinstance(runtime, dict) else None
    floor = (
        re.fullmatch(r">=\s*(\d+(?:\.\d+)*)", requirement)
        if isinstance(requirement, str)
        else None
    )
    if (
        not isinstance(runtime, dict)
        or set(runtime) != {"python_requires"}
        or floor is None
        or tuple(map(int, floor[1].split("."))) < (3, 11)
    ):
        return [
            _diag(
                "PACK_DRIFT",
                "MANIFEST.runtime",
                "shared_pack",
                ["MANIFEST.json"],
                "pack runtime must declare a Python >= version floor of at least 3.11",
                "restore the supported pack runtime contract",
            )
        ]
    return []


def validate_pack(root: Path, workspace: dict[str, Any]) -> list[Diagnostic]:
    lock_rel = workspace["governance_pack"]["lock"]
    try:
        lock = load_manifest_json(root / lock_rel)
        pack_root = root / lock["pack_path"]
        manifest = load_manifest_json(pack_root / "MANIFEST.json")
        version = (
            (pack_root / manifest["version_file"]).read_text(encoding="utf-8").strip()
        )
        actual, files = pack_hash(pack_root)
    except (OSError, ValueError, TypeError, KeyError) as exc:
        return [
            _diag(
                "PACK_DRIFT",
                "governance_pack",
                "shared_pack",
                [lock_rel],
                str(exc),
                "restore a coherent pack and lock",
            )
        ]
    usage = manifest.get("field_usage")
    allowed_usage = {"consumed", "validation-only", "reserved-for-future-use"}
    if (
        not isinstance(usage, dict)
        or set(usage) != set(manifest)
        or set(usage.values()) - allowed_usage
    ):
        return [
            _diag(
                "PACK_DRIFT",
                "MANIFEST.field_usage",
                "shared_pack",
                [f"{lock['pack_path']}/MANIFEST.json"],
                "pack manifest fields are not completely classified",
                "classify every manifest field",
            )
        ]
    # Report broken metadata and a stale lock together: repairing one must not
    # be what first reveals the other.
    issues = _pack_metadata_issues(manifest, files)
    identity = {
        "pack_name": manifest["name"],
        "version": manifest["version"],
        "revision": manifest["revision"],
        "content_sha256": actual,
        "files": files,
        "lock_algorithm": manifest["lock_algorithm"],
    }
    for key, expected in identity.items():
        if lock.get(key) != expected:
            issues.append(
                _diag(
                    "PACK_DRIFT",
                    f"governance_pack.{key}",
                    "shared_pack",
                    [lock_rel, f"{lock['pack_path']}/MANIFEST.json"],
                    "pack lock identity does not match the pack manifest/content",
                    "run gov sync after reviewing the pack update",
                )
            )
            break
    if issues:
        return issues
    if (
        version != manifest["version"]
        or workspace["governance_pack"]["required_schema"] != manifest["schema_version"]
    ):
        return [
            _diag(
                "PACK_DRIFT",
                "governance_pack.version",
                "shared_pack",
                [lock_rel],
                "pack version or required schema disagrees",
                "install a compatible pack or update the workspace contract",
            )
        ]
    if (
        set(manifest["commands"]) != SUPPORTED_COMMANDS
        or set(manifest["diagnostics"]) != DIAGNOSTIC_CLASSES
    ):
        return [
            _diag(
                "PACK_DRIFT",
                "MANIFEST.capabilities",
                "shared_pack",
                [f"{lock['pack_path']}/MANIFEST.json"],
                "declared commands or diagnostics do not match implementation",
                "align the manifest with implemented capabilities",
            )
        ]
    return []


def _python_specifier(root: Path, project: dict[str, Any]) -> tuple[str, str]:
    environment = project["environment"]
    if environment.get("profile", "python") == "tools":
        workspace = load_manifest_json(root / ".governance/workspace.json")
        lock = load_manifest_json(root / workspace["governance_pack"]["lock"])
        source = lock["pack_path"] + "/MANIFEST.json"
        manifest = load_manifest_json(root / source)
        return manifest.get("runtime", {}).get("python_requires", ">=3.11"), source
    authority_id = environment["python_requires"]["authority_ref"]
    authority = next(
        item for item in project["authorities"] if item["authority_id"] == authority_id
    )
    payload = tomllib.loads((root / authority["source"]).read_text(encoding="utf-8"))
    value: object = payload
    for token in environment["python_requires"]["source_field"].split("."):
        value = value[token]
    if not isinstance(value, str):
        raise TypeError("Python requirement must be a string")
    return value, authority["source"]


def interpreter_candidate(
    root: Path,
    project: dict[str, Any],
    local: dict[str, Any] | None,
    environ: dict[str, str] | None = None,
    *,
    local_tiers: dict[str, str] | None = None,
) -> tuple[Path, str]:
    """Select the governed interpreter path before environment validation."""
    env = os.environ if environ is None else environ
    override = project["environment"]["interpreter_override"]
    if env.get(override):
        raw, reason = env[override], f"explicit override {override}"
    elif local and local.get("interpreter"):
        tier = (local_tiers or {}).get("interpreter", "local")
        raw = local["interpreter"]
        reason = OVERLAY_REASONS.get(tier, OVERLAY_REASONS["local"])
    else:
        raw, reason = (
            project["environment"]["preferred"],
            "repository-local environment",
        )
    path = Path(raw).expanduser()
    if not path.is_absolute():
        path = root / path
    return path, reason


def resolve_interpreter(
    root: Path,
    project: dict[str, Any],
    local: dict[str, Any] | None,
    environ: dict[str, str] | None = None,
    *,
    local_tiers: dict[str, str] | None = None,
) -> tuple[Path | None, str, list[Diagnostic]]:
    path, reason = interpreter_candidate(
        root, project, local, environ, local_tiers=local_tiers
    )
    if not path.is_file() or not os.access(path, os.X_OK):
        return (
            None,
            reason,
            [
                _diag(
                    "ENVIRONMENT_MISMATCH",
                    "python.interpreter",
                    "environment",
                    [str(path)],
                    "highest-priority approved interpreter is not executable",
                    "repair or remove that override; no fallback is permitted",
                )
            ],
        )
    try:
        specifier_text, source = _python_specifier(root, project)
        if project["environment"].get("profile", "python") == "tools":
            floor = re.fullmatch(r">=\s*(\d+(?:\.\d+)*)", specifier_text)
            if floor is None:
                raise ValueError(
                    "pack runtime requires an explicit >= Python version floor"
                )
            probe_code = (
                "import platform,runpy; v=platform.python_version(); print(v); "
                "floor=runpy.run_path("
                + repr(str(Path(__file__).resolve().parents[1] / "runtime_floor.py"))
                + ")['numeric_python_floor']; "
                "raise SystemExit(0 if floor(v," + repr(specifier_text) + ") else 1)"
            )
        else:
            probe_code = (
                "import platform; from packaging.specifiers import SpecifierSet; "
                "from packaging.version import Version; v=platform.python_version(); "
                "print(v); raise SystemExit(0 if Version(v) in SpecifierSet("
                + repr(specifier_text)
                + ") else 1)"
            )
        probe = subprocess.run(
            [str(path), "-c", probe_code],
            cwd=root,
            text=True,
            capture_output=True,
            timeout=20,
            check=False,
        )
        version = probe.stdout.strip()
        if probe.returncode != 0:
            return (
                None,
                reason,
                [
                    _diag(
                        "ENVIRONMENT_MISMATCH",
                        "python.version",
                        "project_metadata",
                        [str(path), source],
                        f"interpreter version {version or '<unknown>'} cannot certify {specifier_text}; packaging must be installed in the selected environment",
                        "select an approved interpreter satisfying project metadata",
                    )
                ],
            )
    except (
        OSError,
        ValueError,
        TypeError,
        KeyError,
        StopIteration,
        ImportError,
        subprocess.TimeoutExpired,
    ) as exc:
        return (
            None,
            reason,
            [
                _diag(
                    "ENVIRONMENT_MISMATCH",
                    "python.version",
                    "environment",
                    [str(path)],
                    str(exc),
                    "repair project metadata or the selected interpreter",
                )
            ],
        )
    return path.absolute(), reason, []


def source_binding_path(root: Path, project: dict[str, Any]) -> Path:
    """Where this checkout's importable source lives.

    Derived from the declared editable source rather than spelled out in each
    caller, so the binding is one governed fact instead of a shell convention
    repeated wherever a Python subprocess is started.
    """
    return (root / project["environment"]["editable_source"]).resolve().parent


def binding_environment(
    root: Path,
    project: dict[str, Any],
    environ: dict[str, str] | None = None,
    local: dict[str, Any] | None = None,
) -> dict[str, str]:
    """Variables that bind execution to this checkout, for the declared mode.

    Source mode prepends this checkout's source so imports resolve here rather
    than through whichever checkout an installed distribution was built from.
    Installed mode adds nothing: the distribution is the binding, and injecting
    a path would make `editable_install` unable to observe a wrong one.
    """
    env = os.environ if environ is None else environ
    if project["environment"].get("profile", "python") == "tools":
        return {}
    mode, _, issues = resolve_execution_mode(project, env, local)
    if issues or mode != "source":
        return {}
    binding = str(source_binding_path(root, project))
    inherited = env.get("PYTHONPATH", "")
    return {
        "PYTHONPATH": os.pathsep.join((binding, inherited)) if inherited else binding
    }


def resolve_execution_mode(
    project: dict[str, Any],
    environ: dict[str, str] | None = None,
    local: dict[str, Any] | None = None,
) -> tuple[str, str, list[Diagnostic]]:
    """Declare whether the project runs from a checkout or an installed dist.

    This is declared and never inferred. Whether the project package is
    imported from a checkout or from an installed distribution decides which dependency
    metadata is truthful, and reading that back from whichever source happens
    to answer would make the check agree with any environment it runs in --
    including one holding another worktree's stale build artifacts.
    """
    env = os.environ if environ is None else environ
    if project["environment"].get("profile", "python") == "tools":
        return "tools", "declared tools profile", []
    spec = project["environment"]["execution_mode"]
    override = spec["override"]
    if env.get(override):
        mode, reason = env[override], f"explicit override {override}"
    elif local and local.get("execution_mode"):
        # A machine whose shared environment owns no source binding is in source
        # mode for every checkout on it, so the declaration belongs beside the
        # interpreter it describes rather than in each shell that starts a gate.
        mode, reason = local["execution_mode"], "approved local overlay"
    else:
        mode, reason = spec["default"], "declared project default"
    if mode not in EXECUTION_MODES:
        return (
            mode,
            reason,
            [
                _diag(
                    "CONFLICTING_VALUE",
                    "environment.execution_mode",
                    "environment",
                    [mode, override],
                    "execution mode is not one this resolver implements",
                    f"select one of {', '.join(EXECUTION_MODES)}",
                )
            ],
        )
    return mode, reason, []


def _route_environment(
    root: Path,
    project: dict[str, Any],
    local: dict[str, Any] | None,
    env: dict[str, str],
    route: str | None,
    project_source: str,
) -> tuple[dict[str, str], list[Diagnostic]]:
    observed: dict[str, str] = {}
    issues: list[Diagnostic] = []
    if route is None:
        return observed, issues
    route_spec = project["routes"].get(route)
    if route_spec is None:
        return observed, [
            _diag(
                "BROKEN_REFERENCE",
                f"route:{route}",
                "project",
                [project_source],
                "requested route is unknown",
                "select a declared route",
            )
        ]
    local_tools = local.get("tool_paths", {}) if local else {}
    for tool_id in route_spec.get("environment_tools", []):
        spec = project["environment"]["route_tools"][tool_id]
        executable = (
            local_tools.get(tool_id)
            or env.get(spec["override"])
            or next(
                (
                    shutil.which(item, path=env.get("PATH"))
                    for item in spec["candidates"]
                    if shutil.which(item, path=env.get("PATH"))
                ),
                None,
            )
        )
        if not executable:
            issues.append(
                _diag(
                    "ENVIRONMENT_MISMATCH",
                    f"tool:{tool_id}",
                    "environment",
                    [spec["override"]],
                    "route-required tool is unavailable",
                    "configure local tool_paths, the approved override, or install it",
                )
            )
            continue
        try:
            version = subprocess.run(
                [executable, *spec["version_args"]],
                cwd=root,
                text=True,
                capture_output=True,
                timeout=20,
                check=False,
                env=env,
            )
            lines = (version.stdout or version.stderr).strip().splitlines()
            if version.returncode or not lines:
                raise ValueError("route tool version probe failed")
            observed[tool_id] = lines[0]
        except (OSError, ValueError, subprocess.TimeoutExpired) as exc:
            issues.append(
                _diag(
                    "ENVIRONMENT_MISMATCH",
                    f"tool:{tool_id}",
                    "environment",
                    [executable],
                    str(exc),
                    "repair the tool",
                )
            )
    return observed, issues


def doctor(  # noqa: C901 - environment probes remain one deterministic transaction
    root: Path,
    project: dict[str, Any],
    local: dict[str, Any] | None,
    environ: dict[str, str] | None = None,
    route: str | None = None,
    *,
    local_tiers: dict[str, str] | None = None,
    project_source: str = UNRESOLVED_PROJECT_SOURCE,
) -> tuple[dict[str, Any], list[Diagnostic]]:
    env = os.environ if environ is None else environ
    interpreter, reason, issues = resolve_interpreter(
        root, project, local, env, local_tiers=local_tiers
    )
    if project["environment"].get("profile", "python") == "tools":
        requirement, requirement_source = _python_specifier(root, project)
        observed, tool_issues = _route_environment(
            root, project, local, env, route, project_source
        )
        issues.extend(tool_issues)
        return {
            "status": "ENVIRONMENT_INVALID" if issues else "ok",
            "interpreter": str(interpreter) if interpreter else None,
            "selection_reason": reason,
            "profile": "tools",
            "python_requires": requirement,
            "python_requires_source": requirement_source,
            "route_tools": observed,
        }, issues
    execution_mode, execution_mode_reason, mode_issues = resolve_execution_mode(
        project, env, local
    )
    issues.extend(mode_issues)
    try:
        requirement, requirement_source = _python_specifier(root, project)
    except (
        OSError,
        ValueError,
        TypeError,
        KeyError,
        StopIteration,
        ImportError,
    ) as exc:
        requirement, requirement_source = None, None
        issues.append(
            _diag(
                "ENVIRONMENT_MISMATCH",
                "python.metadata",
                "project_metadata",
                [],
                str(exc),
                "repair the Python requirement authority",
            )
        )
    result: dict[str, Any] = {
        "status": "ENVIRONMENT_INVALID" if issues else "ok",
        "interpreter": str(interpreter) if interpreter else None,
        "selection_reason": reason,
        "execution_mode": execution_mode,
        "execution_mode_reason": execution_mode_reason,
        "python_requires": requirement,
        "python_requires_source": requirement_source,
        "dependency_authority": project["environment"]["dependency_authority"][
            "source"
        ],
    }
    if interpreter is None:
        return result, issues
    package_import = project["environment"]["package_import"]
    modules = list(
        dict.fromkeys([package_import, *project["environment"]["critical_packages"]])
    )
    # Importability and installed-distribution ownership are different claims,
    # and source mode turns on the second, so the probe reports both.
    probe_code = (
        "import importlib,importlib.metadata,json,pathlib,platform; "
        "mods=[importlib.import_module(n) for n in " + repr(modules) + "]; "
        "root=importlib.import_module(" + repr(package_import) + "); "
        "names={(d.metadata['Name'] or '').lower().replace('_','-') "
        "for d in importlib.metadata.distributions()}; "
        "print(json.dumps({'version':platform.python_version(),"
        "'package':str(pathlib.Path(root.__file__).resolve()),"
        "'critical_packages':" + repr(modules) + ","
        "'distribution':" + repr(package_import.lower().replace("_", "-")) + " in names"
        "}))"
    )
    # The probe must run under the same binding the gate will use, or doctor
    # answers a question about a different environment than the one verified.
    binding = binding_environment(root, project, env, local)
    result["source_binding"] = binding.get("PYTHONPATH")
    try:
        probe = subprocess.run(
            [str(interpreter), "-c", probe_code],
            cwd=root,
            text=True,
            capture_output=True,
            timeout=20,
            check=False,
            # A supplied environment is a complete probe environment.  Merging
            # ambient variables back in here made `doctor(environ=...)` answer
            # a different question from its caller's explicit environment.
            env={**env, **binding},
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        probe = None
        issues.append(
            _diag(
                "ENVIRONMENT_MISMATCH",
                "python.probe",
                "environment",
                [str(interpreter)],
                str(exc),
                "repair the selected environment",
            )
        )
    if probe is not None and probe.returncode != 0:
        issues.append(
            _diag(
                "ENVIRONMENT_MISMATCH",
                "project.import",
                "environment",
                [str(interpreter)],
                "project or critical package import failed",
                "install current-checkout development dependencies",
            )
        )
    elif probe is not None:
        try:
            data = json.loads(probe.stdout)
            actual = Path(data["package"]).parent
        except (ValueError, KeyError, TypeError) as exc:
            issues.append(
                _diag(
                    "ENVIRONMENT_MISMATCH",
                    "python.probe",
                    "environment",
                    [str(interpreter)],
                    str(exc),
                    "repair the environment probe output",
                )
            )
            result["status"] = "ENVIRONMENT_INVALID"
            return result, issues
        result.update(data)
        expected = (root / project["environment"]["editable_source"]).resolve()
        if actual != expected:
            # Under source mode the injected binding precedes site-packages, so
            # a mismatch here means the binding itself is wrong rather than that
            # some other checkout owns an install. Say which, or the repair is
            # the opposite of the actual fault.
            source_bound = execution_mode == "source"
            issues.append(
                _diag(
                    "ENVIRONMENT_MISMATCH",
                    "source_binding" if source_bound else "editable_install",
                    "project",
                    [str(actual), str(expected)],
                    "bound source is not this checkout"
                    if source_bound
                    else "editable install targets another checkout",
                    "correct environment.editable_source or the injected binding"
                    if source_bound
                    else "reinstall editable from the current checkout",
                )
            )
        elif execution_mode == "source" and data.get("distribution"):
            # An installed distribution keeps one checkout authoritative for
            # every other, which is the ownership source mode exists to remove.
            # Importability alone cannot see it: the binding wins the import and
            # the distribution stays, silently, as a second answer.
            issues.append(
                _diag(
                    "ENVIRONMENT_MISMATCH",
                    "source_neutrality",
                    "environment",
                    [str(interpreter)],
                    "source mode requires no installed project distribution",
                    "uninstall the project distribution from the shared environment",
                )
            )
        elif execution_mode == "source" and not binding:
            # Importing this checkout without a declared binding means an
            # installed distribution answered, which is the cross-checkout
            # ownership source mode exists to remove.
            issues.append(
                _diag(
                    "ENVIRONMENT_MISMATCH",
                    "source_binding",
                    "environment",
                    [str(expected)],
                    "source mode resolved no binding, so a distribution answered",
                    "declare execution_mode source only with a resolvable binding",
                )
            )
    observed, tool_issues = _route_environment(
        root, project, local, env, route, project_source
    )
    issues.extend(tool_issues)
    if observed:
        result["route_tools"] = observed
    result["status"] = "ENVIRONMENT_INVALID" if issues else "ok"
    return result, issues
