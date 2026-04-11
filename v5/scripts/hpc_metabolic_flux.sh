#!/bin/bash
#SBATCH --job-name=v5_flux
#SBATCH --array=0-5
#SBATCH --partition=scu-cpu
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -eo pipefail

BASE_DIR="${GENELAB_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"
LOG_DIR="$BASE_DIR/v5/logs"
mkdir -p "$LOG_DIR"
cd "$BASE_DIR"

# 6 LOMO tissues only (lung/colon excluded — single-mission)
TISSUES=(liver gastrocnemius kidney thymus eye skin)
TISSUE=${TISSUES[$SLURM_ARRAY_TASK_ID]}

echo "=== v5 Phase 3: Metabolic Flux Modeling ==="
echo "Tissue: $TISSUE"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "Repo: $BASE_DIR"
echo "Log dir: $LOG_DIR"

# Activate conda
source ~/.bashrc
conda activate perturb_seq_new

python -u "$BASE_DIR/v5/scripts/metabolic_flux.py" --tissue "$TISSUE"

echo "=== Done: $TISSUE ==="
