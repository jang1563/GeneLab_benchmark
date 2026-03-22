#!/bin/bash
#SBATCH --job-name=v4_multi_method
#SBATCH --partition=scu-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --output=v4/logs/M1_%A_%a.out
#SBATCH --error=v4/logs/M1_%A_%a.err
#SBATCH --array=0-55  # 8 tissues × 7 CPU methods = 56 jobs (indices 0-55)

# GeneLabBench v4: Multi-method evaluation (CPU methods)
# Usage: sbatch v4/scripts/hpc_multi_method_cpu.sh [--features gene]

set -eo pipefail
source ~/miniconda3/miniconda3/etc/profile.d/conda.sh
conda activate perturb_seq_new

cd /athena/cayuga_0003/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark

# Feature type from argument (default: gene)
FEATURES=${1:-gene}

# Define arrays
TISSUES=(liver gastrocnemius kidney thymus eye skin lung colon)
METHODS=(elasticnet_lr pca_lr rf xgb svm_rbf knn mlp)

N_TISSUES=${#TISSUES[@]}
N_METHODS=${#METHODS[@]}

# Map SLURM_ARRAY_TASK_ID to tissue × method
TISSUE_IDX=$((SLURM_ARRAY_TASK_ID / N_METHODS))
METHOD_IDX=$((SLURM_ARRAY_TASK_ID % N_METHODS))

TISSUE=${TISSUES[$TISSUE_IDX]}
METHOD=${METHODS[$METHOD_IDX]}

echo "=========================================="
echo "Job: ${SLURM_ARRAY_TASK_ID} → ${TISSUE} / ${METHOD} / ${FEATURES}"
echo "Time: $(date)"
echo "Node: $(hostname)"
echo "=========================================="

# Run evaluation
python -u v4/scripts/multi_method_eval.py \
    --tissue "$TISSUE" \
    --method "$METHOD" \
    --features "$FEATURES" \
    --n-bootstrap 2000 \
    --n-perm 1000

echo "Done: $(date)"
