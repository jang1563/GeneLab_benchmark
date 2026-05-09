#!/usr/bin/env bash
# HPC entrypoint for v8 DECOMPOSE. Runs factorial analog models, Mars
# extrapolation, and optional bootstrap CIs.
set -euo pipefail

RUN_BOOTSTRAP=1

for arg in "$@"; do
  case "$arg" in
    --skip-bootstrap)
      RUN_BOOTSTRAP=0
      ;;
    -h|--help)
      echo "Usage: bash scripts/hpc_v8_decompose.sh [--skip-bootstrap]"
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

echo "[v8 DECOMPOSE 1/3] Factorial analog decomposition"
"$PYTHON_BIN" v8/decompose/factorial_analog.py

echo "[v8 DECOMPOSE 2/3] Mars extrapolation"
"$PYTHON_BIN" v8/decompose/mars_extrapolate.py

if [[ "$RUN_BOOTSTRAP" -eq 1 ]]; then
  echo "[v8 DECOMPOSE 3/3] Mars bootstrap confidence intervals"
  "$PYTHON_BIN" v8/decompose/mars_bootstrap_ci.py
else
  echo "[v8 DECOMPOSE 3/3] Skipping Mars bootstrap confidence intervals"
fi

echo "v8 DECOMPOSE complete. Update v8/provenance/runs/decompose_*.json with HPC metadata and checksums."
