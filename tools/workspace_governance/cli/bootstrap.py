"""Dependency-light startup; deliberately compatible with pre-tomllib Python."""

import json
import re
import runpy
import sys
from pathlib import Path


def main(check_only=False):
    root = Path(__file__).resolve().parents[3]
    sources = ["tools/workspace_governance/MANIFEST.json#/runtime/python_requires"]
    try:
        # Governance has its own runtime. Native and documentation consumers do
        # not need Python application metadata merely to load the pack.
        manifest = json.loads((Path(__file__).parents[1] / "MANIFEST.json").read_text())
        requirement = manifest.get("runtime", {}).get("python_requires", ">=3.11")
        if check_only:
            # Existing environment-creation callers explicitly ask to certify
            # the application's bootstrap before creating its environment.
            text = (root / "pyproject.toml").read_text(encoding="utf-8")
            project = re.search(r"(?ms)^\[project\]\s*$(.*?)(?=^\[|\Z)", text)
            declared = (
                re.search(r'requires-python\s*=\s*["\']([^"\']+)', project.group(1))
                if project
                else None
            )
            if declared is None:
                raise ValueError(
                    "Cannot read project.requires-python lower bound for bootstrap"
                )
            requirement = declared.group(1)
            sources = ["pyproject.toml"]
        # Application specifiers may contain upper bounds and exclusions;
        # this preliminary check only extracts the established lower bound.
        floor = (
            re.search(r">=\s*(\d+(?:\.\d+)*)", requirement)
            if check_only
            else re.fullmatch(r">=\s*(\d+(?:\.\d+)*)", requirement)
        )
        if floor is None:
            raise ValueError("Bootstrap requires an explicit >= Python version floor")
        numeric_floor = runpy.run_path(
            str(Path(__file__).parents[1] / "runtime_floor.py")
        )["numeric_python_floor"]
        if not numeric_floor(
            ".".join(str(n) for n in sys.version_info[:3]), ">=" + floor.group(1)
        ):
            raise ValueError("Bootstrap Python does not satisfy " + requirement)
        if check_only:
            return 0
        sys.path.insert(0, str(root))
        runpy.run_module("tools.workspace_governance.cli.gov", run_name="__main__")
    except (ImportError, OSError, ValueError) as exc:
        print(
            json.dumps(
                {
                    "status": "ENVIRONMENT_INVALID",
                    "diagnostics": [
                        {
                            "code": "ENVIRONMENT_MISMATCH",
                            "fact": "bootstrap.interpreter",
                            "owner": "environment",
                            "sources": sources,
                            "reason": str(exc),
                            "repair": "Select a supported interpreter using the launcher bootstrap override",
                        }
                    ],
                }
            ),
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
