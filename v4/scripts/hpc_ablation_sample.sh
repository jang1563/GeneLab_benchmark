#!/bin/bash
#SBATCH --job-name=v4_abl_s
#SBATCH --partition=scu-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=v4/logs/ABL_S_%A_%a.out
#SBATCH --error=v4/logs/ABL_S_%A_%a.err
#SBATCH --array=0-399  # 5 fractions × 2 methods × 4 tissues × 10 repeats = 400

# GeneLabBench v4: Sample Size Ablation
# 5 fractions × 2 methods × 4 tissues × 10 repeats = 400 jobs

set -eo pipefail
source ~/miniconda3/miniconda3/etc/profile.d/conda.sh
conda activate perturb_seq_new

cd ${GENELAB_ROOT:?Set GENELAB_ROOT to your project directory}

TISSUES=(liver thymus skin gastrocnemius)
METHODS=(pca_lr elasticnet_lr)
FRACTIONS=(0.2 0.4 0.6 0.8 1.0)

N_TISSUES=${#TISSUES[@]}
N_METHODS=${#METHODS[@]}
N_FRACTIONS=${#FRACTIONS[@]}
N_REPEATS=10

TASK_ID=${SLURM_ARRAY_TASK_ID}

# Map: tissue × method × fraction × repeat
# Total per tissue: N_METHODS × N_FRACTIONS × N_REPEATS = 100
# Total: 4 × 100 = 400

TISSUE_IDX=$((TASK_ID / (N_METHODS * N_FRACTIONS * N_REPEATS)))
REMAINDER=$((TASK_ID % (N_METHODS * N_FRACTIONS * N_REPEATS)))
METHOD_IDX=$((REMAINDER / (N_FRACTIONS * N_REPEATS)))
REMAINDER2=$((REMAINDER % (N_FRACTIONS * N_REPEATS)))
FRAC_IDX=$((REMAINDER2 / N_REPEATS))
REPEAT=$((REMAINDER2 % N_REPEATS))

TISSUE=${TISSUES[$TISSUE_IDX]}
METHOD=${METHODS[$METHOD_IDX]}
FRAC=${FRACTIONS[$FRAC_IDX]}

echo "=========================================="
echo "Sample Ablation: ${TISSUE} / ${METHOD} / frac=${FRAC} / repeat=${REPEAT}"
echo "Job: ${TASK_ID}, Time: $(date), Node: $(hostname)"
echo "=========================================="

python -u v4/scripts/ablation_sample_size.py \
    --tissue "$TISSUE" \
    --method "$METHOD" \
    --fraction "$FRAC" \
    --repeat "$REPEAT" \
    --n-bootstrap 2000 \
    --n-perm 1000

echo "Done: $(date)"
