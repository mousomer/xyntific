#!/usr/bin/env bash
set -euo pipefail
GOVERNANCE_REPO_ROOT="$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)"
exec make -C "$GOVERNANCE_REPO_ROOT" check-governance
