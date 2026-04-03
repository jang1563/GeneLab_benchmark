#!/bin/bash
#SBATCH --job-name=v4_tabnet
#SBATCH --partition=scu-gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --gres=gpu:1
#SBATCH --time=4:00:00
#SBATCH --output=v4/logs/M1_tabnet_%A_%a.out
#SBATCH --error=v4/logs/M1_tabnet_%A_%a.err
#SBATCH --array=0-7  # 8 tissues × 1 GPU method = 8 jobs

# GeneLabBench v4: TabNet evaluation (GPU)
# Usage: sbatch v4/scripts/hpc_multi_method_gpu.sh [--features gene]

set -eo pipefail
source ~/miniconda3/miniconda3/etc/profile.d/conda.sh
conda activate perturb_seq_new

cd ${GENELAB_ROOT:?Set GENELAB_ROOT to your project directory}

FEATURES=${1:-gene}

TISSUES=(liver gastrocnemius kidney thymus eye skin lung colon)
TISSUE=${TISSUES[$SLURM_ARRAY_TASK_ID]}

echo "=========================================="
echo "Job: ${SLURM_ARRAY_TASK_ID} → ${TISSUE} / tabnet / ${FEATURES}"
echo "Time: $(date)"
echo "Node: $(hostname)"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'N/A')"
echo "=========================================="

python -u v4/scripts/multi_method_eval.py \
    --tissue "$TISSUE" \
    --method tabnet \
    --features "$FEATURES" \
    --n-bootstrap 2000 \
    --n-perm 1000

echo "Done: $(date)"
