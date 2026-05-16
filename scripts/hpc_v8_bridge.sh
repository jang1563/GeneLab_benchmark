#!/usr/bin/env bash
# HPC entrypoint for v8 BRIDGE. Intended for machines with SpaceOmicsBench,
# processed fGSEA outputs, and enough resources for CV/bootstrap recomputation.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

PYTHON_BIN="${PYTHON_BIN:-python3}"
if [[ "$PYTHON_BIN" == "python3" && -x "$REPO_ROOT/.hpc-venv/bin/python" ]]; then
  PYTHON_BIN="$REPO_ROOT/.hpc-venv/bin/python"
fi

SOB_ROOT="${SPACEOMICS_ROOT:-"$REPO_ROOT/../SpaceOmicsBench"}"
if [[ ! -d "$SOB_ROOT" ]]; then
  echo "SPACEOMICS_ROOT does not resolve to a SpaceOmicsBench checkout: $SOB_ROOT" >&2
  echo "Set SPACEOMICS_ROOT before running this script." >&2
  exit 2
fi

export SPACEOMICS_ROOT="$SOB_ROOT"

echo "[v8 BRIDGE 1/4] Aggregate mouse -> human NES conservation"
"$PYTHON_BIN" v8/bridge/link_spaceomicsbench.py

echo "[v8 BRIDGE 2/4] Per-tissue mouse NES Spearman matrix"
"$PYTHON_BIN" v8/bridge/tissue_nes_bridge.py

echo "[v8 BRIDGE 3/4] Supervised conservation with mouse NES features"
"$PYTHON_BIN" v8/bridge/supervised_conservation.py

echo "[v8 BRIDGE 4/4] Leakage audit for supervised conservation"
"$PYTHON_BIN" v8/bridge/leakage_audit.py

echo "v8 BRIDGE complete. Update v8/provenance/runs/bridge_*.json with HPC metadata and checksums."
