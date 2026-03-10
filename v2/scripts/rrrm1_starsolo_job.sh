#!/bin/bash
#SBATCH --job-name=rrrm1_starsolo
#SBATCH --partition=scu-gpu
#SBATCH --gres=gpu:0
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/athena/masonlab/scratch/users/jak4013/rrrm1_scrna/logs/starsolo_%A_%a.out
#SBATCH --error=/athena/masonlab/scratch/users/jak4013/rrrm1_scrna/logs/starsolo_%A_%a.err
#SBATCH --array=0-3
#
# rrrm1_starsolo_job.sh
# STARsolo alignment for RRRM-1 scRNA-seq (OSD-918/920/924/934)
# STAR 2.7.11b, GRCm39-2024-A index (pre-built)
# 10x Chromium 3' v3 chemistry (16bp CB + 12bp UMI)
#
# Submit: /opt/ohpc/pub/software/slurm/24.05.2/bin/sbatch rrrm1_starsolo_job.sh

set -euo pipefail

# ── Array index → OSD mapping ──────────────────────────────────────────────
OSD_NUMS=(918 920 924 934)
TISSUES=(blood eye muscle skin)
OSD_N="${OSD_NUMS[$SLURM_ARRAY_TASK_ID]}"
TISSUE="${TISSUES[$SLURM_ARRAY_TASK_ID]}"

echo "=== RRRM-1 STARsolo: OSD-${OSD_N} (${TISSUE}) ==="
echo "Job: ${SLURM_JOB_ID}, Array: ${SLURM_ARRAY_TASK_ID}"
date

# ── Paths ──────────────────────────────────────────────────────────────────
SCRATCH="/athena/masonlab/scratch/users/jak4013/rrrm1_scrna"
FASTQ_DIR="${SCRATCH}/OSD-${OSD_N}/fastq"
OUT_DIR="${SCRATCH}/OSD-${OSD_N}/starsolo"
STAR_BIN="/opt/ohpc/pub/software/STAR/2.7.11b/bin/Linux_x86_64/STAR"
GENOME_DIR="/athena/masonlab/scratch/users/jak4013/reference/cellranger/refdata-gex-GRCm39-2024-A/star"
WHITELIST="/athena/masonlab/scratch/users/jak4013/software/cellranger-6.1.1/lib/python/cellranger/barcodes/3M-february-2018.txt.gz"

mkdir -p "${OUT_DIR}" "${SCRATCH}/logs"

# ── Check FASTQs exist ─────────────────────────────────────────────────────
if [[ ! -d "${FASTQ_DIR}" ]] || [[ -z "$(ls "${FASTQ_DIR}"/*.fastq.gz 2>/dev/null)" ]]; then
    echo "ERROR: No FASTQs found at ${FASTQ_DIR}"
    echo "Run: bash rrrm1_osdr_download.sh ${OSD_N}"
    exit 1
fi

# ── Pair R1 (cDNA) and R2 (CB+UMI) files ─────────────────────────────────
# For 10x Chromium: R1 = CB+UMI (28bp), R2 = cDNA read
# STAR convention: --readFilesIn cDNA_reads CB_UMI_reads
R2_FILES=$(ls "${FASTQ_DIR}"/*_R1_*.fastq.gz 2>/dev/null || ls "${FASTQ_DIR}"/*_R1.fastq.gz 2>/dev/null || true)
R1_FILES=$(ls "${FASTQ_DIR}"/*_R2_*.fastq.gz 2>/dev/null || ls "${FASTQ_DIR}"/*_R2.fastq.gz 2>/dev/null || true)

# Fallback: all R1 for CB+UMI, all R2 for cDNA
if [[ -z "$R1_FILES" ]] || [[ -z "$R2_FILES" ]]; then
    echo "Using glob pattern for FASTQ pairing..."
    R1_FILES=$(ls "${FASTQ_DIR}"/*R1*.fastq.gz | tr '\n' ',')
    R2_FILES=$(ls "${FASTQ_DIR}"/*R2*.fastq.gz | tr '\n' ',')
fi

# Remove trailing comma
R1_FILES=$(echo "$R1_FILES" | sed 's/,$//')
R2_FILES=$(echo "$R2_FILES" | sed 's/,$//')

echo "cDNA reads (R2): ${R2_FILES}"
echo "CB+UMI reads (R1): ${R1_FILES}"

# ── STARsolo ──────────────────────────────────────────────────────────────
echo ""
echo "Running STARsolo..."

${STAR_BIN} \
    --runThreadN 16 \
    --genomeDir "${GENOME_DIR}" \
    --readFilesIn ${R2_FILES} ${R1_FILES} \
    --readFilesCommand zcat \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist "${WHITELIST}" \
    --soloCBstart 1 --soloCBlen 16 \
    --soloUMIstart 17 --soloUMIlen 12 \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --soloUMIfiltering MultiGeneUMI_CR \
    --soloCellFilter EmptyDrops_CR \
    --soloFeatures Gene GeneFull \
    --soloOutDir "${OUT_DIR}" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
    --outFileNamePrefix "${OUT_DIR}/" \
    --twopassMode Basic \
    2>&1 | tee "${SCRATCH}/logs/starsolo_${OSD_N}.log"

echo ""
echo "=== STARsolo done for OSD-${OSD_N} (${TISSUE}) ==="
echo "Output: ${OUT_DIR}"
ls "${OUT_DIR}/Solo.out/GeneFull/filtered/" 2>/dev/null || ls "${OUT_DIR}/Solo.out/Gene/filtered/" 2>/dev/null || echo "No filtered output dir found"
date
