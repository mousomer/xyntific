#!/usr/bin/env bash
set -euo pipefail

echo "Checking required files and directories..."

REQUIRED_FILES=(
  ".workspace_policy_version.json"
  "scripts/verify.sh"
  "scripts/check_git_sanitation.sh"
  "scripts/check_policy_compliance.sh"
  "scripts/check_policy_template_drift.sh"
  ".policy/policy_template_hashes.json"
)

REQUIRED_DIRS=(
  "docs"
  "scripts"
  ".policy"
)

FAILED=0

for file in "${REQUIRED_FILES[@]}"; do
  if [ ! -f "$file" ]; then
    echo "Missing required file: $file"
    FAILED=1
  fi
done

for dir in "${REQUIRED_DIRS[@]}"; do
  if [ ! -d "$dir" ]; then
    echo "Missing required directory: $dir"
    FAILED=1
  fi
done

if [ "$FAILED" -ne 0 ]; then
  echo "Policy compliance check failed."
  exit 1
fi

echo "Policy compliance check complete."
