#!/bin/bash
#SBATCH --job-name=e4_multisp
#SBATCH --partition=scu-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark/v3/logs/e4_multispecies.log
#SBATCH --error=/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark/v3/logs/e4_multispecies.err

echo "=== E4 Multi-Species NES Concordance ==="
echo "Node: $(hostname)"
echo "Date: $(date)"

WORK_DIR="/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark"

export LD_LIBRARY_PATH=""
source /home/fs01/jak4013/miniconda3/miniconda3/etc/profile.d/conda.sh
set -eo pipefail
conda activate scgpt_env_new

# Check gseapy
python3 -u -c "import gseapy; print(f'gseapy {gseapy.__version__}')"

# Step 1: Download mouse data from OSDR (if not already present)
echo ""
echo "=== Step 1: Download mouse data ==="
python3 -u "${WORK_DIR}/v3/scripts/download_mouse_e4.py"

# Step 2: Also need ensembl symbol map (copy from processed/ if exists)
if [ ! -f "${WORK_DIR}/processed/ensembl_symbol_map.csv" ]; then
    echo "WARNING: ensembl_symbol_map.csv not found — ENSMUSG→symbol mapping will fail"
fi

# Step 3: Run E4
echo ""
echo "=== Step 2: Run E4 ==="
python3 -u "${WORK_DIR}/v3/scripts/e4_multispecies_nes.py"

echo ""
echo "=== E4/E5 Complete ==="
echo "Date: $(date)"
ls -lh ${WORK_DIR}/v3/evaluation/E4*.json ${WORK_DIR}/v3/evaluation/E5*.json 2>/dev/null
