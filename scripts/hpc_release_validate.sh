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

RUN_V8_SUMMARY=0

for arg in "$@"; do
  case "$arg" in
    --v8-summary)
      RUN_V8_SUMMARY=1
      ;;
    -h|--help)
      echo "Usage: bash scripts/hpc_release_validate.sh [--v8-summary]"
      exit 0
      ;;
    *)
      echo "Unknown argument: $arg" >&2
      exit 2
      ;;
  esac
done

echo "[1/5] Python regression tests"
"$PYTHON_BIN" -m unittest discover -s tests -p 'test_review_fixes.py'

echo "[2/5] Whitespace/conflict-marker diff check"
git diff --check
git diff --cached --check

echo "[3/5] Release hygiene check"
"$PYTHON_BIN" - <<'PY'
import json
from pathlib import Path
import subprocess
import sys

repo = Path.cwd()
forbidden = (".claude/", "v8/bridge/geo_cache/", "__pycache__/")
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

prov_dir = repo / "v8" / "provenance" / "runs"
if prov_dir.exists():
    for manifest in sorted(prov_dir.glob("*.json")):
        text = manifest.read_text()
        rel = manifest.relative_to(repo)
        if "TBD_HPC_RUN" in text:
            offenders.append(f"{rel} contains TBD_HPC_RUN")
        try:
            data = json.loads(text)
        except json.JSONDecodeError as exc:
            offenders.append(f"{rel} is invalid JSON: {exc}")
            continue
        if data.get("status") == "draft":
            offenders.append(f"{rel} remains draft")
        if data.get("status") in {"hpc_validated", "release_candidate", "frozen"}:
            if not data.get("generated_at"):
                offenders.append(f"{rel} missing generated_at")
            validation = data.get("validation", {})
            if validation.get("local_tests_run") is not False:
                offenders.append(f"{rel} should record local_tests_run=false")
            if validation.get("hpc_tests_run") is not True:
                offenders.append(f"{rel} missing hpc_tests_run=true")

if offenders:
    print("Release hygiene failures:")
    for item in offenders:
        print(f"  - {item}")
    sys.exit(1)
print("Release hygiene OK")
PY

echo "[4/5] v8 provenance and beta metadata validation"
"$PYTHON_BIN" scripts/validate_v8_provenance.py

if [[ "$RUN_V8_SUMMARY" -eq 1 ]]; then
  echo "[5/5] Regenerate v8 summary"
  "$PYTHON_BIN" v8/RESULTS_SUMMARY.py
else
  echo "[5/5] Skipping v8 summary regeneration; pass --v8-summary to run it"
fi
