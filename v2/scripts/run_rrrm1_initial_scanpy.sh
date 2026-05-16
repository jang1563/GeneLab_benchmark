#!/bin/bash
set -euo pipefail

CONDA_ENV_BASE="${CONDA_PREFIX:-$HOME/miniconda3}/envs/seurat.v5.R.4.3.3.python.3.11.nfcore"
PYTHON_BIN="${CONDA_ENV_BASE}/bin/python"
SCRIPT="${SCRATCH_DIR:?Set SCRATCH_DIR}/rrrm1_scrna/rrrm1_initial_scanpy.py"

if [[ ! -x "${PYTHON_BIN}" ]]; then
    echo "ERROR: Python env not found at ${PYTHON_BIN}"
    exit 1
fi

if [[ ! -f "${SCRIPT}" ]]; then
    echo "ERROR: Initial scanpy script not found at ${SCRIPT}"
    exit 1
fi

export LD_LIBRARY_PATH="${CONDA_ENV_BASE}/lib:${LD_LIBRARY_PATH:-}"
export MPLBACKEND=Agg

exec "${PYTHON_BIN}" "${SCRIPT}" "$@"
