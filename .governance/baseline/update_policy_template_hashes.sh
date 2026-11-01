#!/usr/bin/env bash
set -euo pipefail

REPO_PATH="${1:-.}"
cd "$REPO_PATH"

if git rev-parse --show-toplevel >/dev/null 2>&1; then
  REPO_ROOT="$(git rev-parse --show-toplevel)"
else
  REPO_ROOT="$(pwd)"
fi
cd "$REPO_ROOT"

HASH_MANIFEST="$REPO_ROOT/.policy/policy_template_hashes.json"
if [[ ! -f "$HASH_MANIFEST" ]]; then
  echo "FAIL: missing template hash manifest: .policy/policy_template_hashes.json" >&2
  exit 1
fi

PYTHON_BIN="${PYTHON_BIN:-python3}"
if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
  if command -v python >/dev/null 2>&1; then
    PYTHON_BIN="python"
  else
    echo "FAIL: python3/python not found" >&2
    exit 1
  fi
fi

"$PYTHON_BIN" - <<'PY'
import hashlib
import json
from pathlib import Path

manifest_path = Path('.policy/policy_template_hashes.json')
data = json.loads(manifest_path.read_text(encoding='utf-8'))
templates = data.get('templates', {})
updated = {}
for rel_path in templates:
    p = Path(rel_path)
    if not p.is_file():
        raise SystemExit(f'FAIL: missing template-tracked file: {rel_path}')
    updated[rel_path] = hashlib.sha256(p.read_bytes()).hexdigest()

data['templates'] = updated
manifest_path.write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8')
print(f'Updated {manifest_path} ({len(updated)} hashes)')
PY
