#!/usr/bin/env bash
# HPC entrypoint for v8 causal integration. Requires regenerated DECOMPOSE
# betas and INTERVENE tissue signatures.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

PYTHON_BIN="${PYTHON_BIN:-python3}"
if [[ "$PYTHON_BIN" == "python3" && -x "$REPO_ROOT/.hpc-venv/bin/python" ]]; then
  PYTHON_BIN="$REPO_ROOT/.hpc-venv/bin/python"
fi

echo "[v8 CAUSAL] Build ICP DAG summaries"
"$PYTHON_BIN" v8/causal/build_dag.py

echo "v8 CAUSAL complete. Update v8/provenance/runs/causal_icp_dag.json with HPC metadata and checksums."
