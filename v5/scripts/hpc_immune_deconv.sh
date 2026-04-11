#!/bin/bash
#SBATCH --job-name=v5_immune
#SBATCH --array=0-7
#SBATCH --partition=scu-cpu
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -eo pipefail

BASE_DIR="${GENELAB_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"
LOG_DIR="$BASE_DIR/v5/logs"
mkdir -p "$LOG_DIR"
cd "$BASE_DIR"

TISSUES=(liver gastrocnemius kidney thymus eye skin lung colon)
TISSUE=${TISSUES[$SLURM_ARRAY_TASK_ID]}

echo "=== v5 Phase 1: Immune Deconvolution ==="
echo "Tissue: $TISSUE"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "Repo: $BASE_DIR"
echo "Log dir: $LOG_DIR"

# Activate conda
source ~/.bashrc
conda activate perturb_seq_new

python -u "$BASE_DIR/v5/scripts/immune_deconv.py" --tissue "$TISSUE"

echo "=== Done: $TISSUE ==="
