#!/bin/bash
#SBATCH --job-name=rrrm2_convert
#SBATCH --partition=scu-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=4:00:00
#SBATCH --output=/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark/v3/logs/rrrm2_convert_%a.log
#SBATCH --error=/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark/v3/logs/rrrm2_convert_%a.err
#SBATCH --array=0-3

# GLDS array: 402, 403, 404, 405
GLDS_ARRAY=(402 403 404 405)
GLDS=${GLDS_ARRAY[$SLURM_ARRAY_TASK_ID]}

echo "=== Task $SLURM_ARRAY_TASK_ID: GLDS-$GLDS ==="
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "Memory limit: ${SLURM_MEM_PER_NODE}MB"

SCRIPT_DIR="/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark/v3/scripts"
LOG_DIR="/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark/v3/logs"
mkdir -p "$LOG_DIR"

# --- Step 1: R extraction (Seurat → MTX/CSV) ---
echo ""
echo "=== Step 1: R extraction ==="
export LD_LIBRARY_PATH=""
source /home/fs01/jak4013/miniconda3/miniconda3/etc/profile.d/conda.sh
set -eo pipefail
conda activate seurat.v5.R.4.3.3

Rscript "${SCRIPT_DIR}/convert_rrrm2_seurat_to_h5ad.R" "$GLDS"
R_EXIT=$?

if [ $R_EXIT -ne 0 ]; then
    echo "ERROR: R extraction failed for GLDS-$GLDS (exit code $R_EXIT)"
    exit 1
fi

echo "R extraction succeeded."

# --- Step 2: Python conversion (MTX → h5ad) ---
echo ""
echo "=== Step 2: Python conversion ==="

# Use the seurat.v5 env's Python with reticulate
# First check if anndata is available, install if not
python3 -c "import anndata" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "Installing anndata and scipy..."
    pip install --quiet anndata scipy
fi

python3 "${SCRIPT_DIR}/mtx_to_h5ad.py" "$GLDS"
PY_EXIT=$?

if [ $PY_EXIT -ne 0 ]; then
    echo "WARNING: Python h5ad conversion failed (exit code $PY_EXIT)"
    echo "MTX/metadata files still available for manual conversion."
    # Don't exit with error — MTX extraction succeeded
fi

echo ""
echo "=== GLDS-$GLDS conversion complete ==="
echo "Date: $(date)"
