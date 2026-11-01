#!/usr/bin/env bash
set -euo pipefail
GOVERNANCE_REPO_ROOT="$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)"
cd "$GOVERNANCE_REPO_ROOT"
GOVERNED_ENV="$(./gov env)"
eval "$GOVERNED_ENV"
exec "$PYTHON_BIN" -I -B -c 'import pathlib,runpy,sys,types; root=pathlib.Path.cwd(); namespace=types.ModuleType("tools"); namespace.__package__="tools"; namespace.__path__=[str(root/"tools")]; sys.modules["tools"]=namespace; sys.path.insert(0,str(root)); runpy.run_module("tools.workspace_governance.baseline",run_name="__main__")' --root .
