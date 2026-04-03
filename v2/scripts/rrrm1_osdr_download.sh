#!/bin/bash
# rrrm1_osdr_download.sh
# Download RRRM-1 scRNA-seq FASTQs from OSDR for a given OSD number
#
# Usage: bash rrrm1_osdr_download.sh <OSD_NUMBER>
# Example: bash rrrm1_osdr_download.sh 934
#
# Target OSD numbers (scRNA-seq only):
#   OSD-918  Blood (peripheral blood)
#   OSD-920  Eye (extraocular muscle, retina, cornea, lens, optic nerve)
#   OSD-924  Limb muscle (hindlimb + forelimb)
#   OSD-934  Skin (dorsal, FACS-sorted keratinocytes)

set -euo pipefail

OSD_N="${1:-}"
if [[ -z "$OSD_N" ]]; then
    echo "Usage: $0 <OSD_NUMBER>"
    exit 1
fi

BASE_URL="https://osdr.nasa.gov"
OSDR_API="${BASE_URL}/osdr/data/osd/files/${OSD_N}"
DOWNLOAD_BASE="${BASE_URL}/geode-py/ws/studies/OSD-${OSD_N}/download?source=datamanager&file="
SCRATCH="${SCRATCH_DIR:?Set SCRATCH_DIR}/rrrm1_scrna"
OUT_DIR="${SCRATCH}/OSD-${OSD_N}/fastq"

mkdir -p "${OUT_DIR}"

echo "=== Downloading OSD-${OSD_N} FASTQs to ${OUT_DIR} ==="

# Get file list from OSDR API
FILELIST=$(curl -s "${OSDR_API}" | python3 -c "
import json, sys
d = json.load(sys.stdin)
studies = d.get('studies', {})
for k, v in studies.items():
    for f in v.get('study_files', []):
        fname = f.get('file_name', '')
        if fname.endswith('.fastq.gz') or fname.endswith('.fastq.gz '):
            print(fname)
")

if [[ -z "$FILELIST" ]]; then
    echo "ERROR: No FASTQ files found for OSD-${OSD_N}"
    exit 1
fi

echo "Files to download:"
echo "$FILELIST"
echo ""

N_FILES=$(echo "$FILELIST" | wc -l)
echo "Total: ${N_FILES} FASTQ files"
echo ""

# Download each file
DOWNLOADED=0
SKIPPED=0
while IFS= read -r fname; do
    fname=$(echo "$fname" | tr -d ' ')
    OUT_FILE="${OUT_DIR}/${fname}"

    if [[ -f "$OUT_FILE" ]]; then
        echo "  SKIP (exists): ${fname}"
        ((SKIPPED++)) || true
        continue
    fi

    URL="${DOWNLOAD_BASE}${fname}"
    echo "  Downloading: ${fname}"
    curl -L -o "${OUT_FILE}" "${URL}" --silent --show-error \
        --retry 3 --retry-delay 10
    ((DOWNLOADED++)) || true
done <<< "$FILELIST"

echo ""
echo "=== Done: ${DOWNLOADED} downloaded, ${SKIPPED} skipped ==="
echo "Output: ${OUT_DIR}"
ls -lh "${OUT_DIR}"
