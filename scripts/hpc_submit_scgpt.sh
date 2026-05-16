#!/bin/bash
# ── GeneLab Benchmark: scGPT LOMO Fine-Tuning (Cornell Cayuga HPC) ────────────
#
# Cluster: Cayuga (Cornell CAC) — Slurm v25.05.0
#
# GPU nodes (scu-gpu partition):
#   g0001:       4× NVIDIA A100 80GB PCIe  → gpu:a100:1
#   g0002-g0003: 4× NVIDIA A40 48GB PCIe   → gpu:a40:1 (default)
#
# ENVIRONMENT CRITICAL NOTE:
#   scGPT env has CXXABI_1.3.15 incompatibility with login node libstdc++.
#   FIX: Add conda lib to LD_LIBRARY_PATH (done in this script).
#   This is needed for IPython → sqlite3 → _sqlite3 linkage.
#
# Conda env: scgpt_env_new (Python 3.10, scgpt 0.2.4, flash_attn 1.0.4)
# Conda path: ${CONDA_PREFIX:-$HOME/miniconda3}
#
# ONE-TIME SETUP (already done as of 2026-03-09):
#   [x] conda env scgpt_env_new: Python 3.10, scgpt 0.2.4, flash_attn 1.0.4
#   [x] Pretrained model: models/scgpt_whole_human/ (vocab.json + best_model.pt)
#   [x] Ortholog table: data/mouse/ensembl_mouse_human_orthologs_with_symbols.tsv
#
# USAGE:
#   # Pilot: thymus only (4 folds, ~4h each)
#   sbatch --array=0-3 scripts/hpc_submit_scgpt.sh
#
#   # Custom tissue:
#   sbatch --export=TASK=A1,TASK_DIR=A1_liver_lomo,FOLDS_STR="RR-1 RR-3 RR-6 RR-8 RR-9 MHU-2" \
#          --array=0-5 scripts/hpc_submit_scgpt.sh
#
#   # Check status:
#   squeue -u $USER
#   sacct -j <JOBID> --format=JobID,State,ExitCode,Elapsed,MaxRSS
#
#   # Interactive debug (A40, 1h):
#   srun -p scu-gpu --gres=gpu:a40:1 --mem=32G --time=01:00:00 --pty bash
# ─────────────────────────────────────────────────────────────────────────────

#SBATCH --job-name=scgpt_lomo
#SBATCH --partition=scu-gpu
#SBATCH --output=logs/scgpt_%A_%a.log
#SBATCH --error=logs/scgpt_%A_%a.err
#SBATCH --time=06:00:00        # 6h per fold (generous; actual ~1-2h)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G              # 32G CPU RAM (model is ~200MB, data is small)
#SBATCH --gres=gpu:a40:1       # 1x A40 48GB. Swap to gpu:a100:1 if needed.
# NOTE: --array set by sbatch command (e.g. --array=0-3 for thymus)

# ── Configuration (overridable via sbatch --export) ───────────────────────────
TASK="${TASK:-A4}"
TASK_DIR="${TASK_DIR:-A4_thymus_lomo}"
MODEL_VERSION="${MODEL_VERSION:-whole_human}"
EPOCHS="${EPOCHS:-10}"
BATCH_SIZE="${BATCH_SIZE:-8}"      # A40 48GB: 8 comfortable for 12L/512d model
LR="${LR:-1e-4}"
FREEZE_LAYERS="${FREEZE_LAYERS:-10}"  # Freeze bottom 10/12 layers (top2 + head)
SEED="${SEED:-42}"

# Default LOMO folds (thymus A4)
if [ -z "${FOLDS_STR}" ]; then
    FOLDS=("MHU-1" "MHU-2" "RR-6" "RR-9")
else
    read -ra FOLDS <<< "${FOLDS_STR}"
fi
FOLD="${FOLDS[$SLURM_ARRAY_TASK_ID]}"

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT_DIR="${GENELAB_ROOT:?Set GENELAB_ROOT to your project directory}"
CONDA_BASE="${CONDA_PREFIX:-$HOME/miniconda3}"
CONDA_ENV="scgpt_env_new"
CONDA_ENV_PATH="${CONDA_BASE}/envs/${CONDA_ENV}"
PYTHON="${CONDA_ENV_PATH}/bin/python"

# ── CRITICAL: Fix libstdc++ compatibility on Cayuga ───────────────────────────
# scgpt_env_new uses newer ICU library (libicui18n.so.78) that requires
# CXXABI_1.3.15, which is NOT in the system /lib64/libstdc++.so.6 (too old).
# Solution: prepend conda env lib so conda's libstdc++ is found first.
export LD_LIBRARY_PATH="${CONDA_ENV_PATH}/lib:${LD_LIBRARY_PATH}"

echo "======================================"
echo "scGPT LOMO Fine-Tuning — Cayuga"
echo "Job ID: ${SLURM_JOB_ID} | Array: ${SLURM_ARRAY_TASK_ID}"
echo "Task: ${TASK} | Task Dir: ${TASK_DIR}"
echo "Fold: ${FOLD} | Model: scGPT-${MODEL_VERSION}"
echo "Epochs: ${EPOCHS} | Batch: ${BATCH_SIZE} | LR: ${LR} | freeze: ${FREEZE_LAYERS}/12"
echo "LD_LIBRARY_PATH: ${LD_LIBRARY_PATH:0:80}..."
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader 2>/dev/null \
    || echo "nvidia-smi unavailable"
echo "======================================"

# Validate fold
if [ -z "${FOLD}" ]; then
    echo "ERROR: FOLD is empty. SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}, FOLDS=(${FOLDS[*]})"
    exit 1
fi

# Activate conda (via source .bashrc which has conda init)
source "${CONDA_BASE}/etc/profile.d/conda.sh" || {
    echo "ERROR: Cannot source conda profile"
    exit 1
}
conda activate "${CONDA_ENV}" || {
    echo "ERROR: conda env '${CONDA_ENV}' not found"
    exit 1
}

cd "${PROJECT_DIR}" || { echo "ERROR: project dir not found: ${PROJECT_DIR}"; exit 1; }
mkdir -p "${PROJECT_DIR}/logs"

echo "Working dir: $(pwd)"
echo "Python: $(which python) — $(python --version)"
echo "PyTorch: $(python -c 'import torch; print(torch.__version__, "CUDA:", torch.cuda.is_available(), "GPUs:", torch.cuda.device_count())')"
echo "scGPT: $(python -c 'import scgpt; print(scgpt.__version__)')"

# ── Step 1: Tokenize fold if not already done ──────────────────────────────────
echo ""
echo "[Step 1] Tokenizing fold ${FOLD} (will skip if already done)..."
python scripts/scgpt_tokenize.py \
    --task "${TASK}" \
    --task-dir "${TASK_DIR}" \
    --fold "${FOLD}" \
    --model-version "${MODEL_VERSION}" \
    --model-base "${PROJECT_DIR}/models" \
    --data-base "${PROJECT_DIR}/data" \
    --tasks-base "${PROJECT_DIR}/tasks" \
    2>&1 | tee "logs/scgpt_tokenize_${TASK}_${FOLD}.log"

TOKENIZE_EXIT=$?
if [ ${TOKENIZE_EXIT} -ne 0 ]; then
    echo "ERROR: Tokenization failed (exit ${TOKENIZE_EXIT})"
    exit ${TOKENIZE_EXIT}
fi

# ── Step 2: Fine-tune on GPU ───────────────────────────────────────────────────
echo ""
echo "[Step 2] Fine-tuning fold ${FOLD} on GPU..."
python scripts/scgpt_finetune.py \
    --task "${TASK}" \
    --fold "${FOLD}" \
    --model-version "${MODEL_VERSION}" \
    --model-base "${PROJECT_DIR}/models" \
    --tasks-base "${PROJECT_DIR}/tasks" \
    --eval-base "${PROJECT_DIR}/evaluation" \
    --device cuda \
    --epochs "${EPOCHS}" \
    --batch-size "${BATCH_SIZE}" \
    --lr "${LR}" \
    --freeze-layers "${FREEZE_LAYERS}" \
    --seed "${SEED}" \
    2>&1 | tee "logs/scgpt_finetune_${TASK}_${FOLD}.log"

EXIT_CODE=$?
echo ""
echo "======================================"
echo "Done: ${TASK} fold ${FOLD} — exit code: ${EXIT_CODE}"
echo "======================================"
exit ${EXIT_CODE}
