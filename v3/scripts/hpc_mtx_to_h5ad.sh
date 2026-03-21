#!/bin/bash
#SBATCH --job-name=mtx2h5ad
#SBATCH --partition=scu-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=1:00:00
#SBATCH --output=/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark/v3/logs/mtx2h5ad_%a.log
#SBATCH --error=/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark/v3/logs/mtx2h5ad_%a.err
#SBATCH --array=0-3

GLDS_ARRAY=(402 403 404 405)
GLDS=${GLDS_ARRAY[$SLURM_ARRAY_TASK_ID]}

echo "=== Task $SLURM_ARRAY_TASK_ID: GLDS-$GLDS ==="
echo "Node: $(hostname)"
echo "Date: $(date)"

SCRIPT_DIR="/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark/v3/scripts"

export LD_LIBRARY_PATH=""
source /home/fs01/jak4013/miniconda3/miniconda3/etc/profile.d/conda.sh
set -eo pipefail
conda activate scgpt_env_new

python3 "${SCRIPT_DIR}/mtx_to_h5ad.py" "$GLDS"

echo ""
echo "=== GLDS-$GLDS h5ad conversion complete ==="
echo "Date: $(date)"
