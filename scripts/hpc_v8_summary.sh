#!/usr/bin/env bash
# HPC entrypoint for regenerating v8 consolidated summaries after pillar runs.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

PYTHON_BIN="${PYTHON_BIN:-python3}"
if [[ "$PYTHON_BIN" == "python3" && -x "$REPO_ROOT/.hpc-venv/bin/python" ]]; then
  PYTHON_BIN="$REPO_ROOT/.hpc-venv/bin/python"
fi

echo "[v8 SUMMARY] Consolidated results table"
"$PYTHON_BIN" v8/RESULTS_SUMMARY.py

echo "v8 SUMMARY complete. Update v8/provenance/runs/summary_consolidated_results.json with HPC metadata and checksums."
