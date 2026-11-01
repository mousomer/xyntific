"""Artifact-role resolution for the three-type work ontology.

Everything here is a pure function of governance documents that the caller has
already loaded. Keeping file access out lets one implementation resolve both a
treatment captured earlier (the base of a work segment) and the treatment on
disk now, and lets the validators reuse it without an import cycle.

Roles come from declarations and from a fixed bootstrap projection of declared
governance facts. Nothing is inferred from routes, dispatch paths, extensions,
or directory names; a path no rule names is ``unclassified``.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import PurePosixPath
from typing import Any

WORKSPACE_MANIFEST = ".governance/workspace.json"

WORK_TYPES = ("Planning", "Coding", "Manifest")
ARTIFACT_ROLES = (
    "governance_treatment",
    "executable_machinery",
    "product_authority",
    "planning_document",
    "generated",
    "bookkeeping",
)
ROLE_PROVENANCE = ("explicit_declared", "bootstrap_projected", "unclassified")
COMPATIBILITY_RESULTS = ("compatible", "contradiction", "conditional", "unclassified")

# The write-compatibility table of the work-type decision. A cell is either
# unconditional or names the condition the decision attaches to it. Conditional
# cells are reported as such, never folded into a pass or a contradiction,
# because they are where an observed write needs a reviewer's judgement.
WRITE_COMPATIBILITY: dict[str, dict[str, str]] = {
    "Planning": {
        "governance_treatment": "contradiction",
        "executable_machinery": "contradiction",
        "product_authority": "allowed",
        "planning_document": "allowed",
        "generated": "contradiction_unless_justified",
        "bookkeeping": "allowed_if_subordinate",
    },
    "Coding": {
        "governance_treatment": "contradiction",
        "executable_machinery": "allowed",
        "product_authority": "allowed_with_authority_obligations",
        "planning_document": "allowed_if_subordinate",
        "generated": "allowed_through_owning_generator",
        "bookkeeping": "allowed",
    },
    "Manifest": {
        "governance_treatment": "allowed",
        "executable_machinery": "contradiction",
        "product_authority": "contradiction_unless_justified",
        "planning_document": "allowed_if_documenting_manifest_decision",
        "generated": "allowed_through_owning_generator",
        "bookkeeping": "allowed",
    },
}
WRITE_CONDITIONS = frozenset(
    cell
    for row in WRITE_COMPATIBILITY.values()
    for cell in row.values()
    if cell not in {"allowed", "contradiction"}
)

# Pack files that every release must keep among its bootstrap roots: the
# declaration of the pack's own surface and the enforcement that reads it. A
# release may add roots; dropping one of these would let a change to the
# enforcement escape judgement by its base version.
REQUIRED_PACK_ROOTS = frozenset(
    {"MANIFEST.json", "resolver/roles.py", "validators/core.py"}
)

_DRIVE = re.compile(r"^[A-Za-z]:")
_PATTERN_CHARACTERS = frozenset("*?[]#\\")


class RoleDeclarationError(ValueError):
    """A governance document cannot be read as role declarations."""


def declaration_path_defect(path: object) -> str | None:
    """Why ``path`` cannot name one repository file, or None when it can."""
    if not isinstance(path, str) or not path or path != path.strip():
        return "path must be a non-empty string without padding"
    if _PATTERN_CHARACTERS & set(path):
        return "path must name one literal file, not a pattern or pointer"
    if path.startswith("/") or _DRIVE.match(path):
        return "path must be repository-relative"
    if any(part in {"", ".", ".."} for part in path.split("/")):
        return "path must be normalized, without empty, '.' or '..' segments"
    return None


def normalize_repo_path(path: str) -> str:
    """Normalize an observed path; refuse one that leaves the repository."""
    normalized = PurePosixPath(path).as_posix().removeprefix("./").rstrip("/")
    if (
        not normalized
        or normalized == "."
        or normalized.startswith("/")
        or ".." in normalized.split("/")
    ):
        raise RoleDeclarationError(f"path is not repository-relative: {path!r}")
    return normalized


@dataclass(frozen=True)
class RoleResolution:
    role: str
    provenance: str
    source: str | None

    def to_dict(self) -> dict[str, object]:
        return {
            "role": self.role,
            "role_provenance": self.provenance,
            "role_source": self.source,
        }

    @classmethod
    def from_dict(cls, value: object) -> RoleResolution:
        if not isinstance(value, dict):
            raise RoleDeclarationError("role resolution must be an object")
        role, provenance = value.get("role"), value.get("role_provenance")
        source = value.get("role_source")
        if (
            role not in ARTIFACT_ROLES
            or provenance not in ROLE_PROVENANCE[:2]
            or not isinstance(source, str)
        ):
            raise RoleDeclarationError("role resolution is invalid")
        return cls(role, provenance, source)


UNCLASSIFIED = RoleResolution("unclassified", "unclassified", None)


@dataclass(frozen=True)
class TreatmentDocuments:
    """The loaded governance documents that decide artifact roles."""

    project_path: str
    project: dict[str, Any]
    workspace: dict[str, Any]
    pack_manifest: dict[str, Any]
    pack_path: str


@dataclass(frozen=True)
class RoleIndex:
    explicit: dict[str, RoleResolution]
    bootstrap: dict[str, RoleResolution]
    bootstrap_directories: tuple[tuple[str, RoleResolution], ...]
    roots: frozenset[str]

    def resolve(self, path: str) -> RoleResolution:
        normalized = normalize_repo_path(path)
        if normalized in self.explicit:
            return self.explicit[normalized]
        if normalized in self.bootstrap:
            return self.bootstrap[normalized]
        for directory, resolution in self.bootstrap_directories:
            if normalized.startswith(f"{directory}/"):
                return resolution
        return UNCLASSIFIED

    def is_root(self, path: str) -> bool:
        return normalize_repo_path(path) in self.roots

    def to_dict(self) -> dict[str, object]:
        return {
            "explicit": {k: v.to_dict() for k, v in self.explicit.items()},
            "bootstrap": {k: v.to_dict() for k, v in self.bootstrap.items()},
            "bootstrap_directories": [
                [directory, resolution.to_dict()]
                for directory, resolution in self.bootstrap_directories
            ],
            "roots": sorted(self.roots),
        }

    @classmethod
    def from_dict(cls, value: object) -> RoleIndex:
        if not isinstance(value, dict):
            raise RoleDeclarationError("role index must be an object")
        explicit, bootstrap = value.get("explicit"), value.get("bootstrap")
        directories, roots = value.get("bootstrap_directories"), value.get("roots")
        if (
            not isinstance(explicit, dict)
            or not isinstance(bootstrap, dict)
            or not isinstance(directories, list)
            or not isinstance(roots, list)
            or not all(isinstance(root, str) for root in roots)
        ):
            raise RoleDeclarationError("role index is incomplete")
        pairs = []
        for pair in directories:
            if not isinstance(pair, list) or len(pair) != 2:
                raise RoleDeclarationError("role index directory is invalid")
            pairs.append(
                (normalize_repo_path(pair[0]), RoleResolution.from_dict(pair[1]))
            )
        return cls(
            explicit={
                normalize_repo_path(k): RoleResolution.from_dict(v)
                for k, v in explicit.items()
            },
            bootstrap={
                normalize_repo_path(k): RoleResolution.from_dict(v)
                for k, v in bootstrap.items()
            },
            bootstrap_directories=tuple(pairs),
            roots=frozenset(normalize_repo_path(root) for root in roots),
        )


def compatibility(
    work_type: str,
    role: str,
    table: dict[str, dict[str, str]] | None = None,
) -> tuple[str, str | None]:
    """Judge one write: (result, condition) from the compatibility table."""
    if work_type not in WORK_TYPES:
        raise RoleDeclarationError(
            f"work_type must be exactly one of: {', '.join(WORK_TYPES)}"
        )
    if role == "unclassified":
        return "unclassified", None
    if role not in ARTIFACT_ROLES:
        raise RoleDeclarationError(f"unknown artifact role: {role!r}")
    cell = (WRITE_COMPATIBILITY if table is None else table)[work_type][role]
    if cell == "allowed":
        return "compatible", None
    if cell == "contradiction":
        return "contradiction", None
    return "conditional", cell


def validate_compatibility_table(table: object) -> dict[str, dict[str, str]]:
    """Accept a captured table only if it covers every type and role exactly."""
    known = {"allowed", "contradiction", *WRITE_CONDITIONS}
    if (
        not isinstance(table, dict)
        or set(table) != set(WORK_TYPES)
        or any(
            not isinstance(row, dict)
            or set(row) != set(ARTIFACT_ROLES)
            or not set(row.values()) <= known
            for row in table.values()
        )
    ):
        raise RoleDeclarationError("compatibility table is incomplete or unknown")
    return {work_type: dict(row) for work_type, row in table.items()}


def project_manifest_paths(workspace: dict[str, Any]) -> list[tuple[int, str]]:
    """(member index, repository path) of each project manifest in this checkout."""
    paths = []
    for index, member in enumerate(workspace.get("projects", [])):
        if not isinstance(member, dict):
            continue
        joined = PurePosixPath(str(member.get("repository", "."))) / str(
            member.get("manifest", "")
        )
        try:
            paths.append((index, normalize_repo_path(joined.as_posix())))
        except RoleDeclarationError:
            continue  # another repository's manifest is not an artifact here
    return paths


def default_project_manifest(workspace: dict[str, Any]) -> str | None:
    """Repository path of the default project's manifest, when one is declared.

    Several members sharing the default identity still yield the first, so
    validation can load it and report the ambiguity instead of losing it.
    """
    defaults = workspace.get("defaults")
    default_id = defaults.get("project") if isinstance(defaults, dict) else None
    members = [
        member
        for member in workspace.get("projects", [])
        if isinstance(member, dict) and member.get("id") == default_id
    ]
    if not members:
        return None
    joined = PurePosixPath(str(members[0].get("repository", "."))) / str(
        members[0].get("manifest", "")
    )
    try:
        return normalize_repo_path(joined.as_posix())
    except RoleDeclarationError:
        return None


def role_bootstrap_lists(project: dict[str, Any]) -> dict[str, list[str]]:
    """The project's declared instruction roots and bookkeeping records."""
    declared = project.get("role_bootstrap", {})
    if not isinstance(declared, dict):
        raise RoleDeclarationError("role_bootstrap must be an object")
    lists = {}
    for key in ("instruction_roots", "bookkeeping"):
        values = declared.get(key, [])
        if not isinstance(values, list) or not all(isinstance(v, str) for v in values):
            raise RoleDeclarationError(f"role_bootstrap.{key} must list paths")
        lists[key] = values
    return lists


def fixed_root_roles(
    workspace: dict[str, Any], project: dict[str, Any], project_path: str
) -> dict[str, tuple[str, str]]:
    """Roots whose role no declaration may change: path -> (role, source)."""
    fixed: dict[str, tuple[str, str]] = {
        WORKSPACE_MANIFEST: ("governance_treatment", "workspace manifest")
    }
    for index, path in project_manifest_paths(workspace):
        fixed.setdefault(
            path,
            ("governance_treatment", f"{WORKSPACE_MANIFEST}#/projects/{index}"),
        )
    lock = workspace.get("governance_pack", {}).get("lock")
    if isinstance(lock, str):
        fixed.setdefault(
            normalize_repo_path(lock),
            ("generated", f"{WORKSPACE_MANIFEST}#/governance_pack/lock"),
        )
    for index, path in enumerate(role_bootstrap_lists(project)["instruction_roots"]):
        fixed.setdefault(
            normalize_repo_path(path),
            (
                "governance_treatment",
                f"{project_path}#/role_bootstrap/instruction_roots/{index}",
            ),
        )
    return fixed


def _json_pointer_token(value: str) -> str:
    return value.replace("~", "~0").replace("/", "~1")


def _explicit_roles(documents: TreatmentDocuments) -> dict[str, RoleResolution]:
    pack_prefix = f"{documents.pack_path}/"
    explicit: dict[str, RoleResolution] = {}
    project_roles = documents.project.get("artifact_roles", {})
    if not isinstance(project_roles, dict):
        raise RoleDeclarationError("project artifact_roles must be an object")
    for raw_path, role in project_roles.items():
        if role not in ARTIFACT_ROLES:
            raise RoleDeclarationError(f"invalid project role declaration: {raw_path}")
        path = normalize_repo_path(raw_path)
        if path.startswith(pack_prefix):
            continue  # the pack owns the roles of its own files
        explicit[path] = RoleResolution(
            role,
            "explicit_declared",
            f"{documents.project_path}#/artifact_roles/{_json_pointer_token(raw_path)}",
        )
    declarations = documents.pack_manifest.get("artifact_roles")
    if not isinstance(declarations, list):
        raise RoleDeclarationError("pack manifest artifact_roles must be a list")
    for index, declaration in enumerate(declarations):
        if not isinstance(declaration, dict):
            raise RoleDeclarationError(f"pack artifact_roles[{index}] is invalid")
        relative, role = declaration.get("path"), declaration.get("role")
        if not isinstance(relative, str) or role not in ARTIFACT_ROLES:
            raise RoleDeclarationError(f"pack artifact_roles[{index}] is invalid")
        explicit[normalize_repo_path(pack_prefix + relative)] = RoleResolution(
            role,
            "explicit_declared",
            f"{pack_prefix}MANIFEST.json#/artifact_roles/{index}",
        )
    return explicit


def _authority_roles(
    documents: TreatmentDocuments, bootstrap: dict[str, RoleResolution]
) -> list[tuple[str, RoleResolution]]:
    authorities = documents.project.get("authorities")
    if not isinstance(authorities, list):
        raise RoleDeclarationError("project authorities must be a list")
    directories = []
    for index, authority in enumerate(authorities):
        if not isinstance(authority, dict) or not isinstance(
            authority.get("source"), str
        ):
            raise RoleDeclarationError(f"authorities[{index}] is invalid")
        source_type = authority.get("source_type")
        if source_type not in {"file", "directory"}:
            continue
        if authority.get("canonical_governance") is True or (
            authority.get("authority_type") == "compatibility"
        ):
            role = "governance_treatment"
        elif authority.get("authority_type") == "human":
            role = "product_authority"
        else:
            continue
        resolution = RoleResolution(
            role,
            "bootstrap_projected",
            f"{documents.project_path}#/authorities/{index} "
            f"({authority.get('authority_id')})",
        )
        path = normalize_repo_path(authority["source"])
        if source_type == "directory":
            directories.append((path, resolution))
        else:
            bootstrap.setdefault(path, resolution)
    return directories


def build_role_index(documents: TreatmentDocuments) -> RoleIndex:
    """Resolve every declared and bootstrap-projected role for one treatment.

    Precedence: the pack's declarations for its own files; then the fixed roles
    of bootstrap roots, which no project declaration can change; then project
    declarations; then the bootstrap projection of declared governance facts;
    then declared authority directories, longest first.
    """
    project, workspace = documents.project, documents.workspace
    fixed = fixed_root_roles(workspace, project, documents.project_path)
    explicit = _explicit_roles(documents)
    for path, (role, _) in fixed.items():
        if path in explicit and explicit[path].role != role:
            del explicit[path]

    bootstrap = {
        path: RoleResolution(role, "bootstrap_projected", source)
        for path, (role, source) in fixed.items()
    }
    for index, path in enumerate(role_bootstrap_lists(project)["bookkeeping"]):
        bootstrap.setdefault(
            normalize_repo_path(path),
            RoleResolution(
                "bookkeeping",
                "bootstrap_projected",
                f"{documents.project_path}#/role_bootstrap/bookkeeping/{index}",
            ),
        )
    directories = _authority_roles(documents, bootstrap)
    for index, surface in enumerate(project.get("generated_surfaces", [])):
        target = surface.get("target") if isinstance(surface, dict) else None
        # A pointer target is one generated value inside a larger file, which
        # does not make the whole file generated.
        if isinstance(target, str) and "#" not in target:
            bootstrap.setdefault(
                normalize_repo_path(target),
                RoleResolution(
                    "generated",
                    "bootstrap_projected",
                    f"{documents.project_path}#/generated_surfaces/{index}",
                ),
            )

    pack_roots = documents.pack_manifest.get("bootstrap_roots")
    if not isinstance(pack_roots, list) or not all(
        isinstance(root, str) for root in pack_roots
    ):
        raise RoleDeclarationError("pack manifest bootstrap_roots must list paths")
    roots = set(fixed) | {
        normalize_repo_path(f"{documents.pack_path}/{root}") for root in pack_roots
    }
    return RoleIndex(
        explicit=dict(sorted(explicit.items())),
        bootstrap=dict(sorted(bootstrap.items())),
        bootstrap_directories=tuple(
            sorted(directories, key=lambda item: (-len(item[0]), item[0]))
        ),
        roots=frozenset(roots),
    )


__all__ = [
    "ARTIFACT_ROLES",
    "COMPATIBILITY_RESULTS",
    "REQUIRED_PACK_ROOTS",
    "ROLE_PROVENANCE",
    "WORKSPACE_MANIFEST",
    "WORK_TYPES",
    "WRITE_COMPATIBILITY",
    "WRITE_CONDITIONS",
    "RoleDeclarationError",
    "RoleIndex",
    "RoleResolution",
    "TreatmentDocuments",
    "build_role_index",
    "compatibility",
    "declaration_path_defect",
    "default_project_manifest",
    "fixed_root_roles",
    "normalize_repo_path",
    "project_manifest_paths",
    "role_bootstrap_lists",
    "validate_compatibility_table",
]
