"""Export immutable local policy-kit distributions without consumer-specific facts."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import tempfile
from fnmatch import fnmatch
from pathlib import Path, PurePosixPath

BASELINE_FILES = (
    (
        "scripts/check_git_sanitation_template.sh",
        ".governance/baseline/check_git_sanitation.sh",
    ),
    (
        "scripts/check_policy_compliance_template.sh",
        ".governance/baseline/check_policy_compliance.sh",
    ),
    (
        "scripts/check_policy_template_drift_template.sh",
        ".governance/baseline/check_policy_template_drift.sh",
    ),
    (
        "scripts/update_policy_template_hashes_template.sh",
        ".governance/baseline/update_policy_template_hashes.sh",
    ),
    ("templates/verify_template.sh", ".governance/baseline/verify_template.sh"),
    ("sanitation_baseline.gitignore", ".governance/baseline/sanitation.gitignore"),
    ("POLICY_VERSION.json", ".governance/baseline/legacy-policy-version.json"),
    ("required_files_manifest.json", ".governance/baseline/legacy-required-files.json"),
)


def canonical_bytes(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2) + "\n").encode()


def relative_path(value: str) -> str:
    path = PurePosixPath(value)
    if (
        not value
        or path.is_absolute()
        or any(p in {".", ".."} for p in value.split("/"))
        or "\\" in value
    ):
        raise ValueError("invalid relative export path")
    return path.as_posix()


def source_bytes(root: Path, relative: str) -> bytes:
    path = root / relative_path(relative)
    if path.is_symlink() or not path.resolve().is_relative_to(root.resolve()):
        raise ValueError("export source escapes declared root")
    if not path.is_file():
        raise ValueError("export source is not a regular file: " + relative)
    return path.read_bytes()


def git_revision(root: Path) -> str | None:
    result = subprocess.run(
        ["git", "-C", str(root), "rev-parse", "HEAD"],
        capture_output=True,
        timeout=15,
        check=False,
    )
    return result.stdout.decode().strip() if result.returncode == 0 else None


def legacy_baseline_root(policy_kit: Path) -> Path:
    retained = policy_kit / "compatibility" / "1.1.3"
    return retained if retained.is_dir() else policy_kit


def collect_payload(pack: Path, policy_kit: Path) -> tuple[dict[str, bytes], dict]:
    manifest = json.loads(source_bytes(pack, "MANIFEST.json"))
    declared = [relative_path(item["path"]) for item in manifest["artifact_roles"]]
    if len(declared) != len(set(declared)):
        raise ValueError("duplicate declared pack artifact")
    exclusions = manifest.get(
        "pack_hash_excludes", ["**/__pycache__/**", "**/*.pyc", "**/.DS_Store"]
    )
    actual = set()
    for path in pack.rglob("*"):
        if not path.is_file():
            continue
        relative = path.relative_to(pack).as_posix()
        if not any(
            fnmatch(relative, pattern) or fnmatch("/" + relative, pattern)
            for pattern in exclusions
        ):
            actual.add(relative)
    if set(declared) != actual:
        raise ValueError("pack export does not match declared artifact set")
    payload = {
        "tools/workspace_governance/" + rel: source_bytes(pack, rel)
        for rel in sorted(declared)
    }
    for source, destination in BASELINE_FILES:
        payload[destination] = source_bytes(legacy_baseline_root(policy_kit), source)
    return payload, manifest


def collect_component(root: Path, prefix: str) -> tuple[dict[str, bytes], dict]:
    """Separate kit-only verification/control files from consumer payloads."""
    if root.is_symlink() or not root.is_dir():
        raise ValueError("component root must be a regular directory")
    root = root.resolve()
    payload = {}
    rows = []
    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root).as_posix()
        if (
            any(
                part in {"__pycache__", ".pytest_cache", ".ruff_cache", ".git"}
                for part in path.relative_to(root).parts
            )
            or path.suffix == ".pyc"
            or path.name == ".DS_Store"
        ):
            continue
        if path.is_symlink():
            raise ValueError("symlink in release component")
        if not path.is_file():
            continue
        data = source_bytes(root, relative)
        mode = path.stat().st_mode & 0o777
        if mode not in {0o644, 0o755}:
            raise ValueError("unsupported component source mode")
        source = prefix + "/" + relative
        payload[source] = data
        rows.append(
            {
                "path": relative,
                "source": source,
                "sha256": hashlib.sha256(data).hexdigest(),
                "mode": mode,
            }
        )
    if not rows:
        raise ValueError("release component is empty")
    return payload, {
        "schema_version": 1,
        "files": rows,
        "identity_sha256": hashlib.sha256(canonical_bytes(rows)).hexdigest(),
    }


def declared_rows(manifest: dict) -> list[dict]:
    rows = list(manifest["files"])
    for name in ("verification", "distribution"):
        rows.extend(manifest.get(name, {}).get("files", []))
    return rows


def verify_existing_release(output: Path, manifest: dict) -> None:
    existing_path = output / "release.json"
    if json.loads(source_bytes(output, "release.json")) != manifest:
        raise ValueError(
            "immutable release output already exists with a different identity"
        )
    expected = {existing_path}
    for item in declared_rows(manifest):
        data = source_bytes(output, item["source"])
        path = output / item["source"]
        if (
            hashlib.sha256(data).hexdigest() != item["sha256"]
            or path.stat().st_mode & 0o777 != item["mode"]
        ):
            raise ValueError("existing immutable release bytes or modes have drifted")
        expected.add(path)
    for path in output.rglob("*"):
        if path.is_symlink():
            raise ValueError("symlink in immutable release")
        if path.is_file() and path not in expected:
            raise ValueError("undeclared file in immutable release")


def stage_release(
    output: Path, manifest: dict, source_payload: dict[str, bytes]
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".governance-export-", dir=output.parent
    ) as staging:
        staged = Path(staging) / "release"
        staged.mkdir()
        for item in declared_rows(manifest):
            target = staged / item["source"]
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(source_payload[item["source"]])
            target.chmod(item["mode"])
        (staged / "release.json").write_bytes(canonical_bytes(manifest))
        if output.exists() or output.is_symlink():
            raise ValueError("immutable release output appeared during export")
        os.rename(staged, output)


def build_manifest(
    pack_manifest: dict, policy_kit: Path, version: str, files: list
) -> dict:
    return {
        "schema_version": 1,
        "name": "workspace/policy-kit",
        "version": version,
        "revision": pack_manifest["revision"],
        "release_root": ".",
        "files": files,
        "provenance": {
            "publisher": "workspace-policy-owner",
            "source_identity": "workspace/policy-kit",
            "runtime_name": pack_manifest["name"],
            "runtime_version": pack_manifest["version"],
            "source_state": "retained-export-bytes; incidental source history is separate publication provenance",
            "payload_identity": hashlib.sha256(canonical_bytes(files)).hexdigest(),
            "legacy_policy_version": json.loads(
                source_bytes(legacy_baseline_root(policy_kit), "POLICY_VERSION.json")
            )["policy_version"],
        },
    }


def payload_mode(destination: str, pack: Path, policy_kit: Path) -> int:
    prefix = "tools/workspace_governance/"
    if destination.startswith(prefix):
        source = pack / destination[len(prefix) :]
    else:
        baseline_sources = {target: name for name, target in BASELINE_FILES}
        source = legacy_baseline_root(policy_kit) / baseline_sources[destination]
    mode = source.stat().st_mode & 0o777
    if mode not in {0o644, 0o755}:
        raise ValueError("unsupported runtime source mode")
    return mode


def export_release(
    pack: Path,
    policy_kit: Path,
    output: Path,
    version: str = "1.2.0",
    *,
    verification_root: Path | None = None,
    distribution_root: Path | None = None,
) -> dict:
    payload, pack_manifest = collect_payload(pack.resolve(), policy_kit.resolve())
    files = [
        {
            "source": "payload/" + destination,
            "path": destination,
            "sha256": hashlib.sha256(data).hexdigest(),
            "mode": payload_mode(destination, pack, policy_kit),
        }
        for destination, data in sorted(payload.items())
    ]
    manifest = build_manifest(pack_manifest, policy_kit, version, files)
    source_payload = {"payload/" + name: data for name, data in payload.items()}
    for name, component_root in (
        ("verification", verification_root),
        ("distribution", distribution_root),
    ):
        if component_root is not None:
            component_payload, declaration = collect_component(component_root, name)
            source_payload.update(component_payload)
            manifest[name] = declaration
            manifest["schema_version"] = 2
    identity = hashlib.sha256(canonical_bytes(manifest)).hexdigest()
    if output.is_symlink():
        raise ValueError("immutable release output must not be a symlink")
    existed = output.exists()
    if existed:
        verify_existing_release(output, manifest)
    else:
        stage_release(output, manifest, source_payload)
    return {
        "release": str(output / "release.json"),
        "identity": identity,
        "files": len(files),
        "changed": not existed,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pack", type=Path, required=True)
    parser.add_argument("--policy-kit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--version", default="1.2.0")
    parser.add_argument("--verification-root", type=Path)
    parser.add_argument("--distribution-root", type=Path)
    args = parser.parse_args()
    try:
        result = export_release(
            args.pack,
            args.policy_kit,
            args.output,
            args.version,
            verification_root=args.verification_root,
            distribution_root=args.distribution_root,
        )
    except (OSError, ValueError, KeyError, subprocess.TimeoutExpired) as error:
        parser.exit(1, "release export failed: " + str(error) + "\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
