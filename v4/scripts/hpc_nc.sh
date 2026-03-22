#!/bin/bash
#SBATCH --job-name=v4_nc
#SBATCH --partition=scu-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=v4/logs/NC_%A_%a.out
#SBATCH --error=v4/logs/NC_%A_%a.err
#SBATCH --array=0-127  # NC1(64) + NC2(64) = 128 jobs

# GeneLabBench v4: Negative Controls (NC1 shuffled labels + NC2 random features)
# NC1: 8 tissues × 8 methods = 64 (indices 0-63)
# NC2: 8 tissues × 8 methods = 64 (indices 64-127)

set -eo pipefail
source ~/miniconda3/miniconda3/etc/profile.d/conda.sh
conda activate perturb_seq_new

cd /athena/cayuga_0003/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark

TISSUES=(liver gastrocnemius kidney thymus eye skin lung colon)
METHODS=(elasticnet_lr pca_lr rf xgb svm_rbf knn mlp tabnet)

N_TISSUES=${#TISSUES[@]}
N_METHODS=${#METHODS[@]}
NC1_SIZE=$((N_TISSUES * N_METHODS))  # 64

TASK_ID=${SLURM_ARRAY_TASK_ID}

if [ $TASK_ID -lt $NC1_SIZE ]; then
    # NC1: Shuffled labels
    TISSUE_IDX=$((TASK_ID / N_METHODS))
    METHOD_IDX=$((TASK_ID % N_METHODS))
    TISSUE=${TISSUES[$TISSUE_IDX]}
    METHOD=${METHODS[$METHOD_IDX]}

    echo "=========================================="
    echo "NC1 Shuffled Labels: ${TISSUE} / ${METHOD}"
    echo "Job: ${TASK_ID}, Time: $(date), Node: $(hostname)"
    echo "=========================================="

    python -u v4/scripts/nc_shuffled_labels.py \
        --tissue "$TISSUE" \
        --method "$METHOD" \
        --n-bootstrap 2000 \
        --n-perm 1000
else
    # NC2: Random features
    IDX=$((TASK_ID - NC1_SIZE))
    TISSUE_IDX=$((IDX / N_METHODS))
    METHOD_IDX=$((IDX % N_METHODS))
    TISSUE=${TISSUES[$TISSUE_IDX]}
    METHOD=${METHODS[$METHOD_IDX]}

    echo "=========================================="
    echo "NC2 Random Features: ${TISSUE} / ${METHOD}"
    echo "Job: ${TASK_ID}, Time: $(date), Node: $(hostname)"
    echo "=========================================="

    python -u v4/scripts/nc_random_features.py \
        --tissue "$TISSUE" \
        --method "$METHOD" \
        --n-bootstrap 2000 \
        --n-perm 1000
fi

echo "Done: $(date)"
