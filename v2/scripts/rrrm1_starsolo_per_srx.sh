#!/bin/bash
#SBATCH --job-name=rrrm1_per_srx
#SBATCH --partition=scu-gpu
#SBATCH --gres=gpu:0
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=6:00:00
#SBATCH --output=/athena/masonlab/scratch/users/jak4013/rrrm1_scrna/logs/per_srx_%A_%a.out
#SBATCH --error=/athena/masonlab/scratch/users/jak4013/rrrm1_scrna/logs/per_srx_%A_%a.err
#SBATCH --array=0-31%4
#
# rrrm1_starsolo_per_srx.sh
# Per-animal STARsolo alignment for RRRM-1 scRNA-seq
# 32 tasks: 4 tissues × 8 animals (4 GC + 4 FLT each)
# Output: per-SRX Solo.out matrices with condition labels
#
# Submit: /opt/ohpc/pub/software/slurm/24.05.2/bin/sbatch rrrm1_starsolo_per_srx.sh

set -euo pipefail

# ── Lookup tables (ordered: OSD-918 ×8, OSD-920 ×8, OSD-924 ×8, OSD-934 ×8) ──
OSD_NUMS=(
  918 918 918 918 918 918 918 918
  920 920 920 920 920 920 920 920
  924 924 924 924 924 924 924 924
  934 934 934 934 934 934 934 934
)
TISSUES=(
  blood blood blood blood blood blood blood blood
  eye   eye   eye   eye   eye   eye   eye   eye
  muscle muscle muscle muscle muscle muscle muscle muscle
  skin  skin  skin  skin  skin  skin  skin  skin
)
SRX_IDS=(
  SRX28491856 SRX28491916 SRX28491975 SRX28491923
  SRX28491934 SRX28491993 SRX28491866 SRX28491877
  SRX28492002 SRX28492003 SRX28492004 SRX28492005
  SRX28492006 SRX28492007 SRX28492008 SRX28492009
  SRX28491958 SRX28491959 SRX28491961 SRX28491962
  SRX28491963 SRX28491964 SRX28491965 SRX28491966
  SRX28491990 SRX28491859 SRX28491991 SRX28491992
  SRX28491994 SRX28491995 SRX28491996 SRX28491997
)
CONDITIONS=(
  GC GC GC GC FLT FLT FLT FLT
  GC GC GC GC FLT FLT FLT FLT
  GC GC GC GC FLT FLT FLT FLT
  GC GC GC GC FLT FLT FLT FLT
)
AGE_MONTHS=(
  3 8 3 8 3 8 3 8
  3 8 3 8 3 8 3 8
  3 8 3 8 3 8 3 8
  3 8 3 8 3 8 3 8
)
SOURCE_NAMES=(
  HC_LAR_09 HC_LAR_19 HC_LAR_10 HC_LAR_20
  FL_LAR_09 FL_LAR_19 FL_LAR_10 FL_LAR_20
  HC_LAR_09 HC_LAR_19 HC_LAR_10 HC_LAR_20
  FL_LAR_09 FL_LAR_19 FL_LAR_10 FL_LAR_20
  HC_LAR_09 HC_LAR_19 HC_LAR_10 HC_LAR_20
  FL_LAR_09 FL_LAR_19 FL_LAR_10 FL_LAR_20
  HC_LAR_09 HC_LAR_19 HC_LAR_10 HC_LAR_20
  FL_LAR_09 FL_LAR_19 FL_LAR_10 FL_LAR_20
)

TASK="${SLURM_ARRAY_TASK_ID}"
OSD_N="${OSD_NUMS[$TASK]}"
TISSUE="${TISSUES[$TASK]}"
SRX="${SRX_IDS[$TASK]}"
CONDITION="${CONDITIONS[$TASK]}"
AGE="${AGE_MONTHS[$TASK]}"
SOURCE="${SOURCE_NAMES[$TASK]}"

echo "=== RRRM-1 per-SRX STARsolo ==="
echo "Task:      $TASK"
echo "OSD:       OSD-${OSD_N}  tissue=${TISSUE}"
echo "SRX:       ${SRX}  condition=${CONDITION}  age=${AGE}mo  source=${SOURCE}"
date

SCRATCH="/athena/masonlab/scratch/users/jak4013/rrrm1_scrna"
FASTQ_DIR="${SCRATCH}/OSD-${OSD_N}/fastq"
OUT_DIR="${SCRATCH}/OSD-${OSD_N}/starsolo_per_srx/${SRX}"
STAR_BIN="/opt/ohpc/pub/software/STAR/2.7.3a/bin/Linux_x86_64_static/STAR"
GENOME_DIR="/athena/masonlab/scratch/users/jak4013/reference/cellranger/refdata-gex-GRCm39-2024-A/star"
WHITELIST="/athena/masonlab/scratch/users/jak4013/reference/3M-february-2018.txt"

mkdir -p "${OUT_DIR}" "${SCRATCH}/logs"

# Find this SRX's R1 and R2 FASTQs
# File naming: GLDS-746_scRNA-Seq_SRX{ID}_R1_raw.fastq.gz
GLDS_ID="GLDS-746"
# For OSD-920, might use different GLDS ID — auto-detect
R1=$(find "${FASTQ_DIR}" -maxdepth 1 -name "*${SRX}_R1*.fastq.gz" | sort | head -1)
R2=$(find "${FASTQ_DIR}" -maxdepth 1 -name "*${SRX}_R2*.fastq.gz" | sort | head -1)

if [[ -z "${R1}" ]] || [[ -z "${R2}" ]]; then
    echo "ERROR: Could not find FASTQs for ${SRX} in ${FASTQ_DIR}"
    echo "Available files:"
    ls "${FASTQ_DIR}" | grep -i "${SRX}" || echo "  (none matching ${SRX})"
    exit 1
fi

echo "R1 (CB+UMI): ${R1}"
echo "R2 (cDNA):   ${R2}"

# STAR soloType CB_UMI_Simple: readFilesIn <cDNA=R2> <CB_UMI=R1>
TWOPASS_MODE="Basic"
if [[ "${OSD_N}" == "934" ]]; then
    TWOPASS_MODE="None"
fi

echo ""
echo "Running STARsolo for ${SRX} ..."

"${STAR_BIN}" \
    --runThreadN 16 \
    --genomeDir "${GENOME_DIR}" \
    --readFilesIn "${R2}" "${R1}" \
    --readFilesCommand zcat \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist "${WHITELIST}" \
    --soloCBstart 1 --soloCBlen 16 \
    --soloUMIstart 17 --soloUMIlen 12 \
    --soloBarcodeReadLength 0 \
    --soloCBmatchWLtype 1MM_multi_pseudocounts \
    --soloUMIfiltering MultiGeneUMI \
    --soloCellFilter CellRanger2.2 3000 0.99 10 \
    --soloFeatures GeneFull \
    --outSAMtype None \
    --outFileNamePrefix "${OUT_DIR}/" \
    --twopassMode "${TWOPASS_MODE}" \
    2>&1 | tee "${SCRATCH}/logs/starsolo_${SRX}.log"

# Write condition metadata file alongside STARsolo output
META_FILE="${OUT_DIR}/sample_metadata.txt"
echo "srx=${SRX}" > "${META_FILE}"
echo "osd=OSD-${OSD_N}" >> "${META_FILE}"
echo "tissue=${TISSUE}" >> "${META_FILE}"
echo "condition=${CONDITION}" >> "${META_FILE}"
echo "age_months=${AGE}" >> "${META_FILE}"
echo "source_name=${SOURCE}" >> "${META_FILE}"

echo ""
echo "=== Done: ${SRX} (OSD-${OSD_N} ${TISSUE} ${CONDITION}) ==="
echo "Output: ${OUT_DIR}"
ls "${OUT_DIR}/Solo.out/GeneFull/filtered/" 2>/dev/null || echo "WARNING: No filtered output dir"
date
