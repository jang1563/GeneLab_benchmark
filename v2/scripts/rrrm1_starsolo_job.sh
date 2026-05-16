#!/bin/bash
#SBATCH --job-name=rrrm1_starsolo
#SBATCH --partition=scu-gpu
#SBATCH --gres=gpu:0
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --array=0-3%1
#
# rrrm1_starsolo_job.sh
# STARsolo alignment for RRRM-1 scRNA-seq (OSD-918/920/924/934)
# STAR 2.7.3a, GRCm39-2024-A index (pre-built)
# 10x Chromium 3' v3 chemistry (R1=16bp CB + 12bp UMI, R2=cDNA)
#
# Submit: /opt/ohpc/pub/software/slurm/24.05.2/bin/sbatch rrrm1_starsolo_job.sh

set -euo pipefail

OSD_NUMS=(918 920 924 934)
TISSUES=(blood eye muscle skin)
OSD_N="${OSD_NUMS[$SLURM_ARRAY_TASK_ID]}"
TISSUE="${TISSUES[$SLURM_ARRAY_TASK_ID]}"

echo "=== RRRM-1 STARsolo: OSD-${OSD_N} (${TISSUE}) ==="
echo "Job: ${SLURM_JOB_ID}, Array: ${SLURM_ARRAY_TASK_ID}"
date

SCRATCH="${SCRATCH_DIR:?Set SCRATCH_DIR}/rrrm1_scrna"
FASTQ_DIR="${SCRATCH}/OSD-${OSD_N}/fastq"
OUT_DIR="${SCRATCH}/OSD-${OSD_N}/starsolo"
STAR_BIN="/opt/ohpc/pub/software/STAR/2.7.3a/bin/Linux_x86_64_static/STAR"
GENOME_DIR="${SCRATCH_DIR:?Set SCRATCH_DIR}/reference/cellranger/refdata-gex-GRCm39-2024-A/star"
WHITELIST="${SCRATCH_DIR:?Set SCRATCH_DIR}/reference/3M-february-2018.txt"

mkdir -p "${OUT_DIR}" "${SCRATCH}/logs"

if [[ ! -d "${FASTQ_DIR}" ]] || [[ -z "$(ls "${FASTQ_DIR}"/*.fastq.gz 2>/dev/null)" ]]; then
    echo "ERROR: No FASTQs found at ${FASTQ_DIR}"
    echo "Run: bash rrrm1_osdr_download.sh ${OSD_N}"
    exit 1
fi

# 10x convention on these files:
#   R1 = barcode+UMI (28 bp)
#   R2 = cDNA read
# STAR convention:
#   --readFilesIn cDNA_reads CB_UMI_reads
CDNA_FILES=$(find "${FASTQ_DIR}" -maxdepth 1 -type f \( -name "*_R2_*.fastq.gz" -o -name "*_R2.fastq.gz" \) | sort | paste -sd, -)
CBUMI_FILES=$(find "${FASTQ_DIR}" -maxdepth 1 -type f \( -name "*_R1_*.fastq.gz" -o -name "*_R1.fastq.gz" \) | sort | paste -sd, -)

if [[ -z "${CDNA_FILES}" ]] || [[ -z "${CBUMI_FILES}" ]]; then
    echo "ERROR: Could not build paired R2/R1 FASTQ lists in ${FASTQ_DIR}"
    exit 1
fi

N_CDNA=$(echo "${CDNA_FILES}" | tr ',' '\n' | wc -l)
N_CBUMI=$(echo "${CBUMI_FILES}" | tr ',' '\n' | wc -l)
if [[ "${N_CDNA}" -ne "${N_CBUMI}" ]]; then
    echo "ERROR: cDNA file count (${N_CDNA}) != CB+UMI file count (${N_CBUMI})"
    exit 1
fi

echo "cDNA reads (R2): ${CDNA_FILES}"
echo "CB+UMI reads (R1): ${CBUMI_FILES}"
echo "Paired libraries: ${N_CDNA}"

TWOPASS_MODE="Basic"
if [[ "${OSD_N}" == "934" ]]; then
    # OSD-934 overflows STAR 2.7.3a 2-pass junction insertion limits.
    TWOPASS_MODE="None"
fi

# ── STARsolo ──────────────────────────────────────────────────────────────
echo ""
echo "Running STARsolo..."

"${STAR_BIN}" \
    --runThreadN 16 \
    --genomeDir "${GENOME_DIR}" \
    --readFilesIn "${CDNA_FILES}" "${CBUMI_FILES}" \
    --readFilesCommand zcat \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist "${WHITELIST}" \
    --soloCBstart 1 --soloCBlen 16 \
    --soloUMIstart 17 --soloUMIlen 12 \
    --soloBarcodeReadLength 0 \
    --soloCBmatchWLtype 1MM_multi_pseudocounts \
    --soloUMIfiltering MultiGeneUMI \
    --soloCellFilter CellRanger2.2 3000 0.99 10 \
    --soloFeatures Gene GeneFull \
    --outSAMtype None \
    --outFileNamePrefix "${OUT_DIR}/" \
    --twopassMode "${TWOPASS_MODE}" \
    2>&1 | tee "${SCRATCH}/logs/starsolo_${OSD_N}.log"

echo ""
echo "=== STARsolo done for OSD-${OSD_N} (${TISSUE}) ==="
echo "Output: ${OUT_DIR}"
ls "${OUT_DIR}/Solo.out/GeneFull/filtered/" 2>/dev/null || ls "${OUT_DIR}/Solo.out/Gene/filtered/" 2>/dev/null || echo "No filtered output dir found"
date
