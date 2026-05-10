#!/usr/bin/env bash
# HPC entrypoint for v8 DECOMPOSE. Audits raw analog cache readiness, then runs
# factorial analog models, Mars extrapolation, bounded-dose sensitivity, and
# optional bootstrap CIs.
set -euo pipefail

RUN_BOOTSTRAP=1
AUDIT_ONLY=0
ALLOW_INCOMPLETE_CACHE=0

for arg in "$@"; do
  case "$arg" in
    --audit-only)
      AUDIT_ONLY=1
      ;;
    --allow-incomplete-cache)
      ALLOW_INCOMPLETE_CACHE=1
      ;;
    --skip-bootstrap)
      RUN_BOOTSTRAP=0
      ;;
    -h|--help)
      echo "Usage: bash scripts/hpc_v8_decompose.sh [--audit-only] [--allow-incomplete-cache] [--skip-bootstrap]"
      exit 0
      ;;
    *)
      echo "Unknown argument: $arg" >&2
      exit 2
      ;;
  esac
done

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

PYTHON_BIN="${PYTHON_BIN:-python3}"
if [[ "$PYTHON_BIN" == "python3" && -x "$REPO_ROOT/.hpc-venv/bin/python" ]]; then
  PYTHON_BIN="$REPO_ROOT/.hpc-venv/bin/python"
fi

echo "[v8 DECOMPOSE 0/5] Raw analog cache audit"
"$PYTHON_BIN" v8/decompose/raw_cache_audit.py

if [[ "$AUDIT_ONLY" -eq 1 ]]; then
  echo "v8 DECOMPOSE audit complete."
  exit 0
fi

if [[ "$ALLOW_INCOMPLETE_CACHE" -ne 1 ]]; then
  "$PYTHON_BIN" - <<'PY'
import json
import sys
from pathlib import Path

audit_path = Path("v8/decompose/evaluation/raw_cache_audit.json")
audit = json.loads(audit_path.read_text())
if not audit.get("full_rerun_ready"):
    print(
        "DECOMPOSE raw cache is incomplete; refusing full rerun. "
        "Restore the missing OSDR caches listed in raw_cache_audit.json, "
        "or pass --allow-incomplete-cache for an explicitly partial/debug run.",
        file=sys.stderr,
    )
    sys.exit(3)
PY
fi

echo "[v8 DECOMPOSE 1/5] Factorial analog decomposition"
"$PYTHON_BIN" v8/decompose/factorial_analog.py

echo "[v8 DECOMPOSE 2/5] Mars extrapolation"
"$PYTHON_BIN" v8/decompose/mars_extrapolate.py

echo "[v8 DECOMPOSE 3/5] Mars bounded-dose sensitivity"
"$PYTHON_BIN" v8/decompose/mars_saturation_sensitivity.py

if [[ "$RUN_BOOTSTRAP" -eq 1 ]]; then
  echo "[v8 DECOMPOSE 4/5] Mars bootstrap confidence intervals"
  "$PYTHON_BIN" v8/decompose/mars_bootstrap_ci.py
else
  echo "[v8 DECOMPOSE 4/5] Skipping Mars bootstrap confidence intervals"
fi

echo "v8 DECOMPOSE complete. Update v8/provenance/runs/decompose_*.json with HPC metadata and checksums."
