#!/bin/bash
# ── GeneLab Benchmark: Submit Geneformer LOMO for All 6 Tissues ───────────────
#
# Master orchestration script for Cayuga HPC.
# Submits one Slurm array job per tissue, each array element = one LOMO fold.
#
# Usage:
#   bash scripts/hpc_submit_all_tissues.sh                  # submit all 6 tissues
#   bash scripts/hpc_submit_all_tissues.sh --tissue A4      # single tissue only
#   bash scripts/hpc_submit_all_tissues.sh --tokenize-only  # tokenize on login node
#   bash scripts/hpc_submit_all_tissues.sh --dry-run        # show sbatch commands only
#   bash scripts/hpc_submit_all_tissues.sh --aggregate      # collect results after HPC
#
# Prerequisites:
#   - conda activate mouse_gf
#   - PYTHONPATH includes Mouse-Geneformer source
#   - Run from project root directory
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ── Project paths ────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"
HPC_SCRIPT="${SCRIPT_DIR}/hpc_submit_geneformer.sh"

# ── Tissue Configuration ────────────────────────────────────────────────────
# Format: TASK|TASK_DIR|FOLD1 FOLD2 ...|N_FOLDS
declare -a TISSUES=(
    "A1|A1_liver_lomo|RR-1 RR-3 RR-6 RR-8 RR-9 MHU-2|6"
    "A2|A2_gastrocnemius_lomo|RR-1 RR-5 RR-9|3"
    "A3|A3_kidney_lomo|RR-1 RR-3 RR-7|3"
    "A4|A4_thymus_lomo|RR-6 RR-9 MHU-1 MHU-2|4"
    "A5|A5_skin_lomo|RR-6 RR-7 MHU-2|3"
    "A6|A6_eye_lomo|RR-1 RR-3 TBD|3"
)

# ── Shared hyperparameters ───────────────────────────────────────────────────
MODEL_VERSION="mouse_gf"
EPOCHS=10
BATCH_SIZE=16
LR="2e-5"
FREEZE_LAYERS=4
SEED=42

# ── Parse arguments ──────────────────────────────────────────────────────────
FILTER_TISSUE=""
TOKENIZE_ONLY=false
DRY_RUN=false
AGGREGATE=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --tissue)
            FILTER_TISSUE="$2"
            shift 2
            ;;
        --tokenize-only)
            TOKENIZE_ONLY=true
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --aggregate)
            AGGREGATE=true
            shift
            ;;
        -h|--help)
            echo "Usage: bash $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --tissue TASK    Submit only this tissue (e.g., A4)"
            echo "  --tokenize-only  Tokenize all tissues on login node (no GPU)"
            echo "  --dry-run        Print sbatch commands without submitting"
            echo "  --aggregate      Collect results + compare with baselines"
            echo "  -h, --help       Show this help"
            echo ""
            echo "Tissues: A1 (liver), A2 (gastrocnemius), A3 (kidney),"
            echo "         A4 (thymus), A5 (skin), A6 (eye)"
            echo ""
            echo "Total: 22 LOMO folds across 6 tissues."
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# ── Aggregate mode ───────────────────────────────────────────────────────────
if [ "${AGGREGATE}" = true ]; then
    echo "=== Collecting Geneformer results across all tissues ==="
    cd "${PROJECT_DIR}"
    python scripts/aggregate_geneformer_results.py
    exit $?
fi

# ── Main loop ────────────────────────────────────────────────────────────────
echo "============================================================"
echo "GeneLab Benchmark — Geneformer Multi-Tissue Deployment"
echo "Model: Mouse-Geneformer (${MODEL_VERSION})"
echo "Config: epochs=${EPOCHS}, batch=${BATCH_SIZE}, lr=${LR}, freeze=${FREEZE_LAYERS}"
if [ -n "${FILTER_TISSUE}" ]; then
    echo "Filter: ${FILTER_TISSUE} only"
fi
if [ "${TOKENIZE_ONLY}" = true ]; then
    echo "Mode: TOKENIZE ONLY (login node, no GPU)"
fi
if [ "${DRY_RUN}" = true ]; then
    echo "Mode: DRY RUN (no submission)"
fi
echo "============================================================"
echo ""

mkdir -p "${PROJECT_DIR}/logs"

SUBMITTED=0
TOTAL_FOLDS=0

for entry in "${TISSUES[@]}"; do
    IFS='|' read -r TASK TASK_DIR FOLDS_STR N_FOLDS <<< "${entry}"

    # Filter by tissue if specified
    if [ -n "${FILTER_TISSUE}" ] && [ "${TASK}" != "${FILTER_TISSUE}" ]; then
        continue
    fi

    ARRAY_MAX=$((N_FOLDS - 1))
    echo "── ${TASK}: ${TASK_DIR} (${N_FOLDS} folds: ${FOLDS_STR}) ──"

    if [ "${TOKENIZE_ONLY}" = true ]; then
        # Tokenize on login node (CPU, no sbatch)
        echo "  Tokenizing ${TASK} (all folds)..."
        if [ "${DRY_RUN}" = true ]; then
            echo "  [DRY RUN] python scripts/geneformer_tokenize.py --task ${TASK} --task-dir ${TASK_DIR} --model-version ${MODEL_VERSION}"
        else
            python scripts/geneformer_tokenize.py \
                --task "${TASK}" \
                --task-dir "${TASK_DIR}" \
                --model-version "${MODEL_VERSION}" \
                2>&1 | tee "logs/tokenize_${TASK}_all.log" || true
        fi
    else
        # Submit array job via sbatch
        # Use export + --export=ALL to preserve shell environment (PATH, conda, etc.)
        # Slurm's --export=VAR=val strips all other env vars; --export=ALL keeps them.
        if [ "${DRY_RUN}" = true ]; then
            echo "  [DRY RUN] export TASK=${TASK} TASK_DIR=${TASK_DIR} FOLDS_STR='${FOLDS_STR}'"
            echo "  sbatch --export=ALL --array=0-${ARRAY_MAX} --job-name=gf_${TASK} ${HPC_SCRIPT}"
        else
            echo "  Submitting ${N_FOLDS} array jobs..."
            export TASK TASK_DIR FOLDS_STR MODEL_VERSION EPOCHS BATCH_SIZE LR FREEZE_LAYERS SEED
            sbatch --export=ALL \
                --array=0-"${ARRAY_MAX}" \
                --job-name="gf_${TASK}" \
                "${HPC_SCRIPT}"
        fi
    fi

    SUBMITTED=$((SUBMITTED + 1))
    TOTAL_FOLDS=$((TOTAL_FOLDS + N_FOLDS))
    echo ""
done

echo "============================================================"
if [ "${TOKENIZE_ONLY}" = true ]; then
    echo "Tokenization complete: ${SUBMITTED} tissues, ${TOTAL_FOLDS} folds."
    echo ""
    echo "Next: submit GPU jobs:"
    echo "  bash scripts/hpc_submit_all_tissues.sh"
elif [ "${DRY_RUN}" = true ]; then
    echo "Dry run: ${SUBMITTED} tissues, ${TOTAL_FOLDS} total folds."
    echo "Remove --dry-run to actually submit."
else
    echo "Submitted: ${SUBMITTED} tissues, ${TOTAL_FOLDS} total GPU jobs."
    echo ""
    echo "Monitor: squeue -u \$USER"
    echo "After completion: bash scripts/hpc_submit_all_tissues.sh --aggregate"
fi
echo "============================================================"
