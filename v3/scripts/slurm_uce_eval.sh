#!/bin/bash
#SBATCH --job-name=uce_eval
#SBATCH --partition=scu-gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=4:00:00
#SBATCH --output=v3/logs/uce_eval_%j.out
#SBATCH --error=v3/logs/uce_eval_%j.err

set -euo pipefail

PROJECT_DIR="/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark"
CONDA_BASE="/home/fs01/jak4013/miniconda3/miniconda3"
UCE_ENV="$CONDA_BASE/envs/uce_env"
DEFAULT_UCE_ARGS="--all --nlayers 4 --batch-size 25"

echo "=== UCE Foundation Model Evaluation ==="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'N/A')"
echo "Project dir: $PROJECT_DIR"
echo "Conda base: $CONDA_BASE"
echo "UCE env: $UCE_ENV"

source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$UCE_ENV"

which python
python -V
python - <<'PY'
import shutil
import sys
required = ["anndata", "scanpy", "uce"]
missing = []
for name in required:
    try:
        __import__(name)
    except Exception:
        missing.append(name)
print("python_executable", sys.executable)
print("uce_cli", shutil.which("uce-eval-single-anndata"))
if missing:
    raise SystemExit(f"Missing Python packages: {missing}")
if shutil.which("uce-eval-single-anndata") is None:
    raise SystemExit("Missing required CLI: uce-eval-single-anndata")
PY

cd "$PROJECT_DIR"

# Ensure output dirs exist
mkdir -p v3/evaluation v3/logs

read -r -a UCE_ARGV <<< "${UCE_ARGS:-$DEFAULT_UCE_ARGS}"
echo "UCE args: ${UCE_ARGV[*]}"

python -u v3/scripts/uce_eval.py "${UCE_ARGV[@]}" 2>&1

echo "=== UCE evaluation complete: $(date) ==="
