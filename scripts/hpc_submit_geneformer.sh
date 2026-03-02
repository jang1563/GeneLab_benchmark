#!/bin/bash
# ── GeneLab Benchmark: Geneformer LOMO Fine-Tuning (Cornell Cayuga HPC) ───────
#
# Cluster: Cayuga (Cornell CAC) — Slurm v25.05.0
# VPN required; contact scu@med.cornell.edu for access.
#
# GPU nodes (scu-gpu partition covers g0001-g0003):
#   g0001:       128 CPUs, 1024GB RAM, 4× NVIDIA A100 80GB PCIe
#   g0002-g0003: 128 CPUs, 1024GB RAM, 4× NVIDIA A40  48GB PCIe each
#   g0004:       128 CPUs, 1032GB RAM, 4× NVIDIA H100 (restricted — NOT in scu-gpu)
#
# Max time: 7 days on scu-gpu partition
# GRES syntax: --gres=gpu:a40:<n>  or  --gres=gpu:a100:<n>
#
# Storage:
#   Home dir:  /home/fs01/<cwid>     (NFS, NOT backed up — small files only)
#   Athena:    /athena/cayuga_XXXX/scratch/$USER/  ← use this for all job I/O
#   Lab link:  /athena/<labname>     → /athena/cayuga_XXXX
#
# ONE-TIME SETUP (run on login node, NOT in a job):
#   1. Install Miniconda in Athena scratch (NOT in $HOME):
#      mkdir -p /athena/cayuga_XXXX/scratch/$USER/miniconda3
#      wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
#           -O /tmp/miniconda.sh
#      bash /tmp/miniconda.sh -b -u -p /athena/cayuga_XXXX/scratch/$USER/miniconda3
#      source /athena/cayuga_XXXX/scratch/$USER/miniconda3/etc/profile.d/conda.sh
#
#   2. Create geneformer environment (Python 3.11, CUDA-compatible torch):
#      conda create -n geneformer python=3.11 -y
#      conda activate geneformer
#      pip install torch --index-url https://download.pytorch.org/whl/cu124
#      pip install transformers datasets accelerate scikit-learn safetensors tqdm huggingface_hub
#
#   3. Copy project to Athena scratch:
#      cp -r /path/to/GeneLab_benchmark /athena/cayuga_XXXX/scratch/$USER/
#
#   4. Pre-tokenize on login node (CPU only, fast ~5 min):
#      cd /athena/cayuga_XXXX/scratch/$USER/GeneLab_benchmark
#      conda activate geneformer
#      python scripts/geneformer_tokenize.py --task A4 --model-version v1
#      # Tokenized data saved to tasks/A4_thymus_lomo/fold_*/geneformer_tokens/v1/
#
# SUBMIT ARRAY JOB (4 LOMO folds simultaneously):
#   sbatch scripts/hpc_submit_geneformer.sh
#
# CHECK JOB STATUS:
#   squeue -u $USER
#   sacct -j <JOBID> --format=JobID,State,ExitCode,Elapsed,MaxRSS
#
# INTERACTIVE DEBUG SESSION (A40, 1h):
#   srun -p scu-gpu --gres=gpu:a40:1 --mem=32G --time=01:00:00 --pty bash
#
# NOTE: fold_RR-23_holdout is the HELD-OUT evaluation fold — excluded from LOMO array.
# ─────────────────────────────────────────────────────────────────────────────

#SBATCH --job-name=geneformer_lomo
#SBATCH --partition=scu-gpu
#SBATCH --output=logs/geneformer_%A_%a.log
#SBATCH --error=logs/geneformer_%A_%a.err
#SBATCH --time=04:00:00                   # 4h per fold; increase to 08:00:00 if needed
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G                         # A40/A100 nodes: 1024GB total — 32G is safe
#SBATCH --gres=gpu:a40:1                 # 1x A40 48GB. Swap to gpu:a100:1 for A100.
#SBATCH --array=0-3                      # 4 LOMO folds (RR-23 is held-out, excluded)

# ── Configuration ─────────────────────────────────────────────────────────────
# 4 LOMO folds only — RR-23 is the held-out evaluation fold and must NOT be
# included in cross-validation training.
FOLDS=("MHU-1" "MHU-2" "RR-6" "RR-9")
FOLD=${FOLDS[$SLURM_ARRAY_TASK_ID]}
TASK="A4"
MODEL_VERSION="v1"
EPOCHS=10
BATCH_SIZE=16        # A40 48GB: 16 comfortable; try 32 if memory allows
LR="2e-5"
FREEZE_LAYERS=4      # Freeze bottom 4/6 BERT layers (top-2 + head trainable)
                     # Recommended for small-n (n≈30-60): limits overfitting
                     # Use 0 for full fine-tuning if n>200; use 6 for head-only
SEED=42

# ── EDIT THESE: Cayuga account details ────────────────────────────────────────
# Replace cayuga_XXXX with your actual project ID (e.g. cayuga_0123)
ATHENA_SCRATCH="/athena/cayuga_XXXX/scratch/${USER}"
PROJECT_DIR="${ATHENA_SCRATCH}/GeneLab_benchmark"
CONDA_BASE="${ATHENA_SCRATCH}/miniconda3"   # Miniconda in Athena, NOT $HOME
CONDA_ENV="geneformer"
# ─────────────────────────────────────────────────────────────────────────────

echo "======================================"
echo "Geneformer LOMO Fine-Tuning — Cayuga"
echo "Job ID: ${SLURM_JOB_ID} | Array: ${SLURM_ARRAY_TASK_ID}"
echo "Fold: ${FOLD} | Task: ${TASK} | Model: Geneformer-${MODEL_VERSION}"
echo "Epochs: ${EPOCHS} | Batch: ${BATCH_SIZE} | LR: ${LR} | freeze_layers: ${FREEZE_LAYERS}"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader 2>/dev/null \
    || echo "nvidia-smi unavailable"
echo "======================================"

# Activate conda (installed in Athena, NOT $HOME)
source "${CONDA_BASE}/etc/profile.d/conda.sh" \
    || { echo "ERROR: conda not found at ${CONDA_BASE}"; exit 1; }
conda activate "${CONDA_ENV}" \
    || { echo "ERROR: conda env '${CONDA_ENV}' not found"; exit 1; }

cd "${PROJECT_DIR}" || { echo "ERROR: project dir not found: ${PROJECT_DIR}"; exit 1; }
mkdir -p logs

echo "Working dir: $(pwd)"
echo "Python: $(which python) — $(python --version)"
echo "PyTorch: $(python -c 'import torch; print(torch.__version__, "CUDA:", torch.cuda.is_available())')"

# Step 1: Tokenize fold if not already done
# (Safe to re-run; will skip if tokenized data already exists)
echo ""
echo "[Step 1] Tokenizing fold ${FOLD} (will skip if already done)..."
python scripts/geneformer_tokenize.py \
    --task "${TASK}" \
    --fold "${FOLD}" \
    --model-version "${MODEL_VERSION}" \
    2>&1 | tee "logs/tokenize_${FOLD}.log"

# Step 2: Fine-tune on GPU
echo ""
echo "[Step 2] Fine-tuning fold ${FOLD} on GPU..."
python scripts/geneformer_finetune.py \
    --task "${TASK}" \
    --fold "${FOLD}" \
    --model-version "${MODEL_VERSION}" \
    --device cuda \
    --epochs "${EPOCHS}" \
    --batch-size "${BATCH_SIZE}" \
    --lr "${LR}" \
    --freeze-layers "${FREEZE_LAYERS}" \
    --seed "${SEED}" \
    2>&1 | tee "logs/finetune_${FOLD}.log"

EXIT_CODE=$?
echo ""
echo "======================================"
echo "Done: Fold ${FOLD} — exit code: ${EXIT_CODE}"
echo "======================================"
exit ${EXIT_CODE}
