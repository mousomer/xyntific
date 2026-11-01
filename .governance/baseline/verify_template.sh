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

echo "Running policy compliance (template baseline)..."
"$REPO_ROOT/scripts/check_policy_compliance.sh" "$REPO_ROOT"

if [[ -x "$REPO_ROOT/scripts/check_policy_template_drift.sh" ]]; then
  echo "Running policy template drift check..."
  "$REPO_ROOT/scripts/check_policy_template_drift.sh" "$REPO_ROOT"
fi

if [[ -x "$REPO_ROOT/scripts/check_policy_compliance_repo.sh" ]]; then
  echo "Running policy compliance (repo extensions)..."
  "$REPO_ROOT/scripts/check_policy_compliance_repo.sh" "$REPO_ROOT"
fi

echo "Running lint..."
# <lint command>

echo "Running tests..."
# <test command>

echo "Verification complete."
