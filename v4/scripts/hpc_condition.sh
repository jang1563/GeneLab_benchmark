#!/bin/bash
#SBATCH --job-name=v4_cond
#SBATCH --partition=scu-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=v4/logs/COND_%A_%a.out
#SBATCH --error=v4/logs/COND_%A_%a.err
#SBATCH --array=0-127  # D3(96) + D5(32) = 128 jobs

# GeneLabBench v4: Condition Prediction (D3 mission ID + D5 hardware)
# D3: 6 tissues × 8 methods × 2 features = 96 (indices 0-95)
# D5: 2 tissues × 8 methods × 2 features = 32 (indices 96-127)

set -eo pipefail
source ~/miniconda3/miniconda3/etc/profile.d/conda.sh
conda activate perturb_seq_new

cd ${GENELAB_ROOT:?Set GENELAB_ROOT to your project directory}

D3_TISSUES=(liver gastrocnemius kidney thymus eye skin)
D5_TISSUES=(liver thymus)
METHODS=(elasticnet_lr pca_lr rf xgb svm_rbf knn mlp tabnet)
FEATURES=(gene pathway_hallmark)

N_D3_TISSUES=${#D3_TISSUES[@]}
N_D5_TISSUES=${#D5_TISSUES[@]}
N_METHODS=${#METHODS[@]}
N_FEATURES=${#FEATURES[@]}

D3_SIZE=$((N_D3_TISSUES * N_METHODS * N_FEATURES))  # 96

TASK_ID=${SLURM_ARRAY_TASK_ID}

if [ $TASK_ID -lt $D3_SIZE ]; then
    # D3: Mission ID prediction
    TISSUE_IDX=$((TASK_ID / (N_METHODS * N_FEATURES)))
    REMAINDER=$((TASK_ID % (N_METHODS * N_FEATURES)))
    METHOD_IDX=$((REMAINDER / N_FEATURES))
    FEAT_IDX=$((REMAINDER % N_FEATURES))

    TISSUE=${D3_TISSUES[$TISSUE_IDX]}
    METHOD=${METHODS[$METHOD_IDX]}
    FEAT=${FEATURES[$FEAT_IDX]}

    echo "=========================================="
    echo "D3 Mission ID: ${TISSUE} / ${METHOD} / ${FEAT}"
    echo "Job: ${TASK_ID}, Time: $(date), Node: $(hostname)"
    echo "=========================================="

    python -u v4/scripts/condition_prediction_v4.py \
        --task D3 \
        --tissue "$TISSUE" \
        --method "$METHOD" \
        --features "$FEAT" \
        --n-bootstrap 2000 \
        --n-perm 1000
else
    # D5: Hardware prediction
    IDX=$((TASK_ID - D3_SIZE))
    TISSUE_IDX=$((IDX / (N_METHODS * N_FEATURES)))
    REMAINDER=$((IDX % (N_METHODS * N_FEATURES)))
    METHOD_IDX=$((REMAINDER / N_FEATURES))
    FEAT_IDX=$((REMAINDER % N_FEATURES))

    TISSUE=${D5_TISSUES[$TISSUE_IDX]}
    METHOD=${METHODS[$METHOD_IDX]}
    FEAT=${FEATURES[$FEAT_IDX]}

    echo "=========================================="
    echo "D5 Hardware: ${TISSUE} / ${METHOD} / ${FEAT}"
    echo "Job: ${TASK_ID}, Time: $(date), Node: $(hostname)"
    echo "=========================================="

    python -u v4/scripts/condition_prediction_v4.py \
        --task D5 \
        --tissue "$TISSUE" \
        --method "$METHOD" \
        --features "$FEAT" \
        --n-bootstrap 2000 \
        --n-perm 1000
fi

echo "Done: $(date)"
