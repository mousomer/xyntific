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
  PYTHON_BIN="python"
fi

"$PYTHON_BIN" - <<'PY'
import hashlib
import json
from pathlib import Path
import sys

repo_root = Path(".")
manifest_path = repo_root / ".policy" / "policy_template_hashes.json"
data = json.loads(manifest_path.read_text(encoding="utf-8"))
algo = str(data.get("hash_algorithm", "sha256")).strip().lower()
templates = data.get("templates", {})

if algo != "sha256":
    print(f"FAIL: unsupported hash algorithm in {manifest_path}: {algo}", file=sys.stderr)
    sys.exit(1)

failed = False
for rel_path, expected in templates.items():
    p = repo_root / rel_path
    if not p.is_file():
        print(f"FAIL: missing template-tracked file: {rel_path}", file=sys.stderr)
        failed = True
        continue
    actual = hashlib.sha256(p.read_bytes()).hexdigest()
    if actual != str(expected).strip().lower():
        print(
            f"FAIL: template drift detected for {rel_path} "
            f"({actual} != {expected})",
            file=sys.stderr,
        )
        failed = True

if failed:
    sys.exit(1)

print(f"PASS: policy template hashes match ({len(templates)} files)")
PY
