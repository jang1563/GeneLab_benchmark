#!/bin/bash
# HPC Setup for v7: install PyG + scPRINT-2 in perturb_seq_new conda env.
# Run this ONCE on the login node before submitting SLURM jobs.
# Usage: bash v7/scripts/hpc_setup_v7.sh

set -eo pipefail

BASE_DIR="${GENELAB_ROOT:?Set GENELAB_ROOT}"
MODEL_DIR="$BASE_DIR/v7/models/scprint2"
mkdir -p "$MODEL_DIR"

echo "=========================================="
echo "v7 HPC Setup: PyG + scPRINT-2"
echo "=========================================="

# Activate conda env
source ${CONDA_PREFIX:-$HOME/miniconda3}/etc/profile.d/conda.sh
conda activate perturb_seq_new

PYTHON=$(which python)
echo "Python: $PYTHON"
echo "PyTorch: $($PYTHON -c 'import torch; print(torch.__version__)')"
TORCH_VER=$($PYTHON -c 'import torch; print(torch.__version__.split("+")[0])')
CUDA_VER=$($PYTHON -c 'import torch; print(torch.version.cuda or "cpu")' 2>/dev/null || echo "cpu")
echo "CUDA: $CUDA_VER"

if [ "$CUDA_VER" = "cpu" ]; then
    PYG_FLAVOR="cpu"
else
    PYG_FLAVOR="cu$(echo "$CUDA_VER" | tr -d '.')"
fi

# ── 1. Install PyTorch Geometric ───────────────────────────────────────────────
echo ""
echo ">>> Installing PyTorch Geometric..."
if $PYTHON -c "import torch_geometric" 2>/dev/null; then
    echo "  PyG already installed: $($PYTHON -c 'import torch_geometric; print(torch_geometric.__version__)')"
else
    # Match the published wheel index to the active torch/cuda build.
    PYG_URL="https://data.pyg.org/whl/torch-${TORCH_VER}+${PYG_FLAVOR}.html"
    echo "  Installing from: $PYG_URL"
    "$PYTHON" -m pip install torch-scatter torch-sparse torch-cluster torch-spline-conv \
        -f "$PYG_URL" --quiet
    "$PYTHON" -m pip install torch-geometric --quiet
    echo "  Installed: $($PYTHON -c 'import torch_geometric; print(torch_geometric.__version__)')"
fi

# ── 2. Install AnnData (needed for scPRINT) ──────────────────────────────────
echo ""
echo ">>> Checking AnnData..."
if $PYTHON -c "import anndata" 2>/dev/null; then
    echo "  AnnData OK: $($PYTHON -c 'import anndata; print(anndata.__version__)')"
else
    "$PYTHON" -m pip install anndata --quiet
fi

# ── 3. Install scPRINT-2 ──────────────────────────────────────────────────────
echo ""
echo ">>> Installing scPRINT-2..."
if $PYTHON -c "import scprint2" 2>/dev/null; then
    echo "  scPRINT-2 already installed"
else
    "$PYTHON" -m pip install scprint2 huggingface_hub --quiet
    echo "  Installed scPRINT-2"
fi

# ── 4. Verify scPRINT-2 package + checkpoint helper ──────────────────────────
echo ""
echo ">>> Verifying scPRINT-2 package..."
$PYTHON - << 'PYEOF'
try:
    from scprint2 import scPRINT2
    print("  scPRINT-2 importable: OK")
    print(f"  Model class: {scPRINT2.__name__}")
except Exception as e:
    print(f"  WARNING: scPRINT-2 import failed: {e}")
    raise
PYEOF

# ── 5. Download default checkpoint on the login node if missing ──────────────
echo ""
echo ">>> Checking default scPRINT-2 checkpoint..."
if [ -f "$MODEL_DIR/medium-v1.5.ckpt" ]; then
    echo "  Found: $MODEL_DIR/medium-v1.5.ckpt"
else
    echo "  Downloading medium-v1.5.ckpt to $MODEL_DIR"
    $PYTHON - << PYEOF
from pathlib import Path

model_dir = Path("$MODEL_DIR")
filename = "medium-v1.5.ckpt"
repo_ids = ["jkobject/scPRINT", "jkobject/scPRINT-2"]

for repo_id in repo_ids:
    try:
        from huggingface_hub import hf_hub_download
        out = hf_hub_download(repo_id=repo_id, filename=filename, local_dir=model_dir)
        print(f"  Downloaded from {repo_id}: {out}")
        break
    except Exception as exc:
        print(f"  Download attempt failed for {repo_id}: {exc}")
else:
    print("  WARNING: could not pre-download checkpoint.")
    print("  The benchmark script will try again at runtime if networking is available.")
PYEOF
fi

# ── 6. Verify all imports for GNN script ─────────────────────────────────────
echo ""
echo ">>> Verifying GNN imports..."
$PYTHON - << 'PYEOF'
import torch
import torch_geometric
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from torch_geometric.nn import GATConv, GCNConv, global_mean_pool
import numpy as np
import pandas as pd
import scipy
import sklearn
print(f"  torch:         {torch.__version__}")
print(f"  torch_geometric: {torch_geometric.__version__}")
print(f"  numpy:         {np.__version__}")
print(f"  sklearn:       {sklearn.__version__}")
print(f"  GPU available: {torch.cuda.is_available()}")
if torch.cuda.is_available():
    print(f"  GPU:           {torch.cuda.get_device_name(0)}")
    print(f"  GPU Memory:    {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
PYEOF

echo ""
echo "=========================================="
echo "Setup complete! You can now submit SLURM jobs."
echo "  bash v7/scripts/hpc_submit_v7.sh"
echo "  bash v7/scripts/hpc_submit_v7.sh --gnn-only"
echo "  bash v7/scripts/hpc_submit_v7.sh --scprint-only"
echo "=========================================="
