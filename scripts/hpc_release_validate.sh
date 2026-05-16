#!/usr/bin/env bash
# Release validation entrypoint intended for HPC or another machine with enough
# local resources. Do not run expensive validation on a constrained laptop.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

PYTHON_BIN="${PYTHON_BIN:-python3}"
if [[ "$PYTHON_BIN" == "python3" && -x "$REPO_ROOT/.hpc-venv/bin/python" ]]; then
  PYTHON_BIN="$REPO_ROOT/.hpc-venv/bin/python"
fi

for arg in "$@"; do
  case "$arg" in
    -h|--help)
      echo "Usage: bash scripts/hpc_release_validate.sh"
      exit 0
      ;;
    *)
      echo "Unknown argument: $arg" >&2
      exit 2
      ;;
  esac
done

echo "[1/3] Python regression tests"
"$PYTHON_BIN" -m unittest discover -s tests -p 'test_review_fixes.py'

echo "[2/3] Whitespace/conflict-marker diff check"
git diff --check
git diff --cached --check

echo "[3/3] Release hygiene check"
"$PYTHON_BIN" - <<'PY'
from pathlib import Path
import subprocess
import sys

repo = Path.cwd()
forbidden = (".claude/", "__pycache__/")
max_bytes = 50 * 1024 * 1024

tracked = subprocess.check_output(["git", "ls-files"], text=True).splitlines()
offenders = []
for rel in tracked:
    path = repo / rel
    if not path.exists():
        continue
    if any(part in rel for part in forbidden):
        offenders.append(rel)
    elif path.stat().st_size > max_bytes:
        offenders.append(f"{rel} ({path.stat().st_size} bytes)")

if offenders:
    print("Release hygiene failures:")
    for item in offenders:
        print(f"  - {item}")
    sys.exit(1)
print("Release hygiene OK")
PY
