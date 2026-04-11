#!/bin/bash
#SBATCH --job-name=scprint2
#SBATCH --array=0-5               # 6 tissues
#SBATCH --partition=scu-gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

# scPRINT-2 Zero-Shot Benchmark — v7 GeneLabBench
# Uses Python API (not CLI) with a local checkpoint + LaminDB+bionty
#
# Results in: GeneLab_benchmark/v7/evaluation/SCPRINT2_{tissue}.json

set -eo pipefail

BASE_DIR="${GENELAB_ROOT:?Set GENELAB_ROOT}"
LOG_DIR="$BASE_DIR/v7/logs"
mkdir -p "$LOG_DIR"

TISSUES=(liver gastrocnemius kidney thymus eye skin)
TISSUE=${TISSUES[$SLURM_ARRAY_TASK_ID]}

echo "Job ID: $SLURM_JOB_ID"
echo "Array task: $SLURM_ARRAY_TASK_ID → tissue=$TISSUE"
echo "Node: $SLURMD_NODENAME"
echo "GPU: $CUDA_VISIBLE_DEVICES"
echo "Repo: $BASE_DIR"
echo "Log dir: $LOG_DIR"
date

source ${CONDA_PREFIX:-$HOME/miniconda3}/etc/profile.d/conda.sh
conda activate perturb_seq_new

# Ensure LaminDB instance is connected (populated on login node)
python -c "import lamindb as ln; ln.connect('scprint_v2')" 2>/dev/null || true

SCRIPT_DIR="$BASE_DIR/v7/unified"
CKPT_PATH="${SCPRINT_CKPT_PATH:-$BASE_DIR/v7/models/scprint2/medium-v1.5.ckpt}"
if [ ! -f "$CKPT_PATH" ] && [ -f "$BASE_DIR/v7/models/scprint2/medium-v1.5_fixed.ckpt" ]; then
    CKPT_PATH="$BASE_DIR/v7/models/scprint2/medium-v1.5_fixed.ckpt"
fi
cd "$BASE_DIR"

python -u "$SCRIPT_DIR/scprint2_benchmark.py" \
    --tissue "$TISSUE" \
    --ckpt-path "$CKPT_PATH" \
    --n-bootstrap 1000 \
    --n-perm 100 \
    --seed 42 \
    --embedding-cache-dir "v7/evaluation/embeddings_cache"

echo "Done: $TISSUE"
date
