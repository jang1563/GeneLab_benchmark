#!/bin/bash
#SBATCH --job-name=f5_rrrm2
#SBATCH --partition=scu-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=8:00:00
#SBATCH --output=/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark/v3/logs/f5_rrrm2.log
#SBATCH --error=/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark/v3/logs/f5_rrrm2.err

echo "=== F5 RRRM-2 Benchmark ==="
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "Memory: $(free -h | head -2)"

WORK_DIR="/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark"
SCRIPT="${WORK_DIR}/v3/scripts/f5_rrrm2_benchmark.py"
DATA_DIR="${WORK_DIR}/v3/data/rrrm2"

# Activate conda env with anndata, scipy, sklearn, gseapy, statsmodels
export LD_LIBRARY_PATH=""
source /home/fs01/jak4013/miniconda3/miniconda3/etc/profile.d/conda.sh
set -eo pipefail
conda activate scgpt_env_new

# Check dependencies
python3 -u -c "
import anndata, scipy, numpy, pandas, gseapy
from sklearn.linear_model import LogisticRegression
from statsmodels.stats.multitest import multipletests
print('All dependencies OK')
print(f'  anndata={anndata.__version__}')
print(f'  gseapy={gseapy.__version__}')
"

echo ""
echo "=== Running F5 tasks ==="
echo "Data dir: ${DATA_DIR}"
ls -lh ${DATA_DIR}/*.h5ad

# Run all F5 tasks (F5a, F5b, F5c run on all data; F5d needs F5b first; F5e is BM-specific)
python3 -u "${SCRIPT}" \
    --tasks F5a,F5b,F5c,F5d,F5e \
    --data-dir "${DATA_DIR}"

echo ""
echo "=== F5 Complete ==="
echo "Date: $(date)"
ls -lh ${WORK_DIR}/v3/evaluation/F5*.json
