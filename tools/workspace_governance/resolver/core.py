from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from tools.workspace_governance.resolver.roles import (
    WORKSPACE_MANIFEST,
    RoleIndex,
    TreatmentDocuments,
    build_role_index,
    normalize_repo_path,
)
from tools.workspace_governance.validators.core import (
    UNRESOLVED_PROJECT_SOURCE,
    Diagnostic,
    load_manifest_json,
    validate_manifests,
    validate_pack,
)

USER_OVERLAY_DIR = "workspace-governance"
# A workspace id becomes a filename, so it may not contain anything that could
# address a different directory. The shell bootstrap resolver applies the same
# rule; both derive one path and must agree on what is addressable.
SAFE_WORKSPACE_ID = re.compile(r"^[A-Za-z0-9._-]+$")


def user_overlay_path(
    workspace_id: str, environ: dict[str, str] | None = None
) -> Path | None:
    """Locate the machine-wide overlay shared by every checkout of a workspace.

    Keyed by workspace identity rather than checkout path so a worktree
    inherits it wherever it lives, and honouring XDG_CONFIG_HOME so a test or
    an operator can substitute an isolated configuration root. Returns None for
    an identity that cannot safely become a filename.
    """
    if not SAFE_WORKSPACE_ID.match(workspace_id):
        return None
    env = os.environ if environ is None else environ
    base = env.get("XDG_CONFIG_HOME") or str(Path.home() / ".config")
    return Path(base).expanduser() / USER_OVERLAY_DIR / f"{workspace_id}.local.json"


class GovernanceError(ValueError):
    def __init__(self, diagnostics: list[Diagnostic]):
        super().__init__(diagnostics[0].reason if diagnostics else "governance error")
        self.diagnostics = diagnostics


@dataclass(frozen=True)
class GovernanceResolver:
    root: Path
    workspace_path: Path
    # None when the workspace names no default project; load() then fails
    # with that reason instead of guessing a conventional location.
    project_path: Path | None
    local_path: Path | None
    user_local_path: Path | None = None

    @classmethod
    def for_root(cls, root: Path) -> GovernanceResolver:
        root = root.resolve()
        workspace_path = root / WORKSPACE_MANIFEST
        project_path: Path | None = None
        user_local: Path | None = None
        try:
            workspace = load_manifest_json(workspace_path)
            default_id = workspace["defaults"]["project"]
            members = [
                item for item in workspace["projects"] if item["id"] == default_id
            ]
            # Several members sharing the default identity still load the
            # first, so validation reports the ambiguity rather than a missing
            # manifest.
            if members:
                project_path = root / members[0]["repository"] / members[0]["manifest"]
            inherited = user_overlay_path(str(workspace["workspace_id"]))
            if inherited is not None and inherited.is_file():
                user_local = inherited
        except (OSError, ValueError, TypeError, KeyError):
            # check() owns the bounded diagnostic; construction remains total.
            pass
        local = root / ".governance/workspace.local.json"
        return cls(
            root,
            workspace_path,
            project_path,
            local if local.exists() else None,
            user_local,
        )

    def local_documents(self) -> list[tuple[dict[str, Any], str, str]]:
        """Return overlays lowest-precedence first as (payload, tier, source)."""
        documents: list[tuple[dict[str, Any], str, str]] = []
        for path, tier in (
            (self.user_local_path, "workspace"),
            (self.local_path, "local"),
        ):
            if path is not None:
                documents.append((load_manifest_json(path), tier, self._source(path)))
        return documents

    def local_overlay(self) -> tuple[dict[str, Any] | None, dict[str, str]]:
        """Layer overlays per key so a repository setting overrides an inherited one.

        Merging by key rather than by document keeps a repository overlay that
        names only `tool_paths` from silently discarding the inherited
        interpreter.
        """
        merged: dict[str, Any] = {}
        tiers: dict[str, str] = {}
        for payload, tier, _ in self.local_documents():
            for key, value in payload.items():
                merged[key], tiers[key] = value, tier
        return (merged or None), tiers

    def _project_manifest(self) -> Path:
        if self.project_path is None:
            raise ValueError("the workspace declares no default project manifest")
        return self.project_path

    def project_source(self) -> str:
        """The project manifest as diagnostics name it."""
        if self.project_path is None:
            return UNRESOLVED_PROJECT_SOURCE
        return self._source(self.project_path)

    def load(self) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any] | None]:
        return (
            load_manifest_json(self.workspace_path),
            load_manifest_json(self._project_manifest()),
            self.local_overlay()[0],
        )

    def treatment_documents(self) -> TreatmentDocuments:
        """Load the documents that decide artifact roles in this checkout.

        Local overlays are deliberately not read: they select machine
        locations and can never change what an artifact is.
        """
        workspace = load_manifest_json(self.workspace_path)
        lock = load_manifest_json(self.root / workspace["governance_pack"]["lock"])
        pack_path = normalize_repo_path(lock["pack_path"])
        return TreatmentDocuments(
            project_path=self.project_source(),
            project=load_manifest_json(self._project_manifest()),
            workspace=workspace,
            pack_manifest=load_manifest_json(self.root / pack_path / "MANIFEST.json"),
            pack_path=pack_path,
        )

    def telemetry_settings(self, environ: dict[str, str] | None = None):
        """Whether and where this checkout's workspace records telemetry.

        Only the machine's local overlays can enable recording; an unreadable
        workspace or overlay leaves it off rather than failing the command.
        """
        # Imported here: the telemetry module builds on this one.
        from tools.workspace_governance.telemetry import resolve_settings

        try:
            workspace = load_manifest_json(self.workspace_path)
            local, _ = self.local_overlay()
        except (OSError, ValueError, TypeError, KeyError):
            return resolve_settings(None, None, environ)
        workspace_id = workspace.get("workspace_id")
        return resolve_settings(
            local, workspace_id if isinstance(workspace_id, str) else None, environ
        )

    def role_index(self) -> RoleIndex:
        """Artifact roles under the treatment currently in this checkout."""
        return build_role_index(self.treatment_documents())

    def _source(self, path: Path) -> str:
        try:
            return path.relative_to(self.root).as_posix()
        except ValueError:
            return str(path)

    def check(self, *, generated: bool = False) -> list[Diagnostic]:
        try:
            documents = self.local_documents()
        except (OSError, ValueError, TypeError, KeyError) as exc:
            # Name the overlays themselves: a malformed inherited overlay must
            # not be reported against the tracked workspace manifests.
            return [
                Diagnostic(
                    "BROKEN_REFERENCE",
                    "manifest",
                    "local",
                    tuple(
                        self._source(path)
                        for path in (self.user_local_path, self.local_path)
                        if path is not None
                    ),
                    str(exc),
                    "restore valid strict JSON in the local interpreter overlays",
                )
            ]
        try:
            workspace, project, local = self.load()
            lock = load_manifest_json(self.root / workspace["governance_pack"]["lock"])
            pack_root = self.root / lock["pack_path"]
        except (OSError, ValueError, TypeError, KeyError) as exc:
            return [
                Diagnostic(
                    "BROKEN_REFERENCE",
                    "manifest",
                    "workspace/project",
                    (
                        self._source(self.workspace_path),
                        self.project_source(),
                    ),
                    str(exc),
                    "restore valid strict JSON manifests and a resolvable pack lock",
                )
            ]
        return validate_manifests(
            self.root,
            workspace,
            project,
            local,
            pack_root=pack_root,
            local_documents=[(payload, source) for payload, _, source in documents],
            local_tiers={
                key: tier for payload, tier, _ in documents for key in payload
            },
            generated=generated,
        ) + validate_pack(self.root, workspace)

    def _verification_fact(
        self, spec: dict[str, Any], authority_map: dict[str, Any], fact: str
    ) -> dict[str, Any]:
        if "command" in spec:
            return {
                "value": spec["command"],
                "owner": "project",
                "source": f"{self.project_source()}#/verification/{fact}/command",
            }
        authority_id, pointer = spec["authority_ref"], spec["json_pointer"]
        authority = authority_map[authority_id]
        value: Any = load_manifest_json(self.root / authority["source"])
        try:
            for token in pointer.removeprefix("/").split("/"):
                value = value[token.replace("~1", "/").replace("~0", "~")]
        except (KeyError, TypeError) as exc:
            raise GovernanceError(
                [
                    Diagnostic(
                        "BROKEN_REFERENCE",
                        pointer,
                        authority_id,
                        (authority["source"],),
                        f"authority JSON pointer does not resolve: {exc}",
                        "repair the reference or restore the authority field",
                    )
                ]
            ) from exc
        return {
            "value": value,
            "owner": authority_id,
            "source": f"{authority['source']}#{pointer}",
        }

    def resolve(
        self, *, task: str | None = None, mode: str | None = None
    ) -> dict[str, Any]:
        issues = self.check()
        if issues:
            raise GovernanceError(issues)
        workspace, project, local = self.load()
        project_source = self.project_source()
        workspace_source = self._source(self.workspace_path)
        profiles = project["execution"]["profiles"]
        selected = mode or project["execution"]["default_mode"]
        matched_scenario: str | None = None
        matched_scenario_priority: int | None = None
        scenario_routes: list[str] | None = None
        if task:
            task_lower = task.lower()
            applicable = [
                scenario
                for scenario in project["execution"]["representative_scenarios"]
                if all(token.lower() in task_lower for token in scenario["match_all"])
            ]
            if applicable:
                highest_priority = max(scenario["priority"] for scenario in applicable)
                highest = [
                    scenario
                    for scenario in applicable
                    if scenario["priority"] == highest_priority
                ]
                if len(highest) != 1:
                    scenario_ids = sorted(scenario["id"] for scenario in highest)
                    raise GovernanceError(
                        [
                            Diagnostic(
                                "AMBIGUOUS_AUTHORITY",
                                "execution.representative_scenarios",
                                "project",
                                (project_source,),
                                (
                                    "multiple applicable scenarios share highest "
                                    f"priority {highest_priority}: "
                                    f"{', '.join(scenario_ids)}"
                                ),
                                "assign distinct priorities to overlapping scenarios",
                            )
                        ]
                    )
                scenario = highest[0]
                (
                    selected,
                    matched_scenario,
                    matched_scenario_priority,
                    scenario_routes,
                ) = (
                    scenario["mode"],
                    scenario["id"],
                    scenario["priority"],
                    scenario["routes"],
                )
        if selected not in profiles:
            raise GovernanceError(
                [
                    Diagnostic(
                        "BROKEN_REFERENCE",
                        "execution.mode",
                        "project",
                        (project_source,),
                        f"unknown execution mode {selected}",
                        "select a declared execution profile",
                    )
                ]
            )

        authority_map = {item["authority_id"]: item for item in project["authorities"]}
        route_ids = scenario_routes or profiles[selected]["routes"]
        routes = {route_id: project["routes"][route_id] for route_id in route_ids}
        authority_ids = sorted(
            {item for route in routes.values() for item in route["authority_refs"]}
        )

        def entry(value: Any, owner: str, source: str) -> dict[str, Any]:
            return {"value": value, "owner": owner, "source": source}

        entries = {
            "workspace.id": entry(
                workspace["workspace_id"], "workspace", workspace_source
            ),
            "project.id": entry(project["project"]["id"], "project", project_source),
            "project.name": entry(
                project["project"]["name"], "project", project_source
            ),
            "execution.mode": entry(
                selected, "project", f"{project_source}#/execution"
            ),
            "execution.exploration": entry(
                profiles[selected]["exploration"],
                "project",
                f"{project_source}#/execution/profiles/{selected}",
            ),
            "execution.final_evidence_rule": entry(
                profiles[selected]["final_evidence_rule"],
                "project",
                f"{project_source}#/execution/profiles/{selected}",
            ),
            "verification.canonical": self._verification_fact(
                project["verification"]["canonical"], authority_map, "canonical"
            ),
            "verification.full": self._verification_fact(
                project["verification"]["full"], authority_map, "full"
            ),
            "verification.targeted": entry(
                project["verification"]["targeted"],
                "project",
                f"{project_source}#/verification/targeted",
            ),
            "verification.rule": entry(
                project["verification"]["rule"],
                "project",
                f"{project_source}#/verification/rule",
            ),
            "routes": entry(routes, "project", f"{project_source}#/routes"),
            "authorities": entry(
                [authority_map[item] for item in authority_ids],
                "project",
                f"{project_source}#/authorities",
            ),
            "matched_scenario": entry(
                matched_scenario,
                "project",
                f"{project_source}#/execution/representative_scenarios",
            ),
            "matched_scenario_priority": entry(
                matched_scenario_priority,
                "project",
                f"{project_source}#/execution/representative_scenarios",
            ),
            "local_overlay": entry(
                self._source(self.local_path) if local and self.local_path else None,
                "local",
                self._source(self.local_path)
                if self.local_path
                else ".governance/workspace.local.json",
            ),
        }
        return {"schema_version": 1, "entries": entries}

    def explain(
        self, query: str | None = None, *, task: str | None = None
    ) -> dict[str, Any]:
        resolved = self.resolve(task=task)
        if task:
            return resolved
        if query is None:
            raise ValueError("explain requires a key or --task")
        if query.startswith("authority:"):
            authority_id = query.split(":", 1)[1]
            _, project, _ = self.load()
            for value in project["authorities"]:
                if value["authority_id"] == authority_id:
                    source = self.project_source()
                    return {
                        "schema_version": 1,
                        "entries": {
                            query: {
                                "value": value,
                                "owner": "project",
                                "source": f"{source}#/authorities",
                            }
                        },
                    }
            raise GovernanceError(
                [
                    Diagnostic(
                        "BROKEN_REFERENCE",
                        query,
                        "project",
                        (self.project_source(),),
                        "unknown authority ID",
                        "use a declared stable authority ID",
                    )
                ]
            )
        if query not in resolved["entries"]:
            raise GovernanceError(
                [
                    Diagnostic(
                        "BROKEN_REFERENCE",
                        query,
                        "resolved",
                        ("ephemeral resolution",),
                        "unknown resolved key",
                        "request a key emitted by gov resolve",
                    )
                ]
            )
        return {"schema_version": 1, "entries": {query: resolved["entries"][query]}}


__all__ = ["Diagnostic", "GovernanceError", "GovernanceResolver"]
