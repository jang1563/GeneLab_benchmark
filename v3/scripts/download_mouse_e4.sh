#!/bin/bash
# Download mouse tissue data needed for E4 from OSDR
# Run on Cayuga: bash download_mouse_e4.sh

set -eo pipefail
BASE="${SCRATCH_DIR:?Set SCRATCH_DIR}/huggingface/benchmark/GeneLab_benchmark/data/mouse"
OSDR="https://osdr.nasa.gov/osdr/data/osd/files"

download_pair() {
    local tissue=$1 mission=$2 glds_num=$3 suffix=$4
    local dir="${BASE}/${tissue}/${mission}"
    mkdir -p "$dir"

    local glds="GLDS-${glds_num}"
    local nc="${glds}_rna_seq_Normalized_Counts${suffix}.csv"
    local st="${glds}_rna_seq_SampleTable${suffix}.csv"

    if [ -f "${dir}/${nc}" ]; then
        echo "EXISTS: ${tissue}/${mission}/${nc}"
        return
    fi

    echo "Downloading: ${tissue}/${mission} (${glds})..."

    # OSDR file URL pattern
    local osd_url="${OSDR}/${glds_num}"
    wget -q -O "${dir}/${nc}" "${osd_url}/${nc}" 2>/dev/null || \
        curl -sfL -o "${dir}/${nc}" "${osd_url}/${nc}" 2>/dev/null || \
        echo "FAILED: ${nc}"

    wget -q -O "${dir}/${st}" "${osd_url}/${st}" 2>/dev/null || \
        curl -sfL -o "${dir}/${st}" "${osd_url}/${st}" 2>/dev/null || \
        echo "FAILED: ${st}"

    ls -lh "${dir}/${nc}" "${dir}/${st}" 2>/dev/null
}

echo "=== Downloading mouse tissue data for E4 ==="

# liver
download_pair "liver" "RR-1" "48" "_GLbulkRNAseq"
download_pair "liver" "RR-3" "137" ""

# thymus
download_pair "thymus" "MHU-2" "289" "_GLbulkRNAseq"
download_pair "thymus" "RR-6" "244" "_GLbulkRNAseq"

# kidney
download_pair "kidney" "RR-1" "102" "_GLbulkRNAseq"
download_pair "kidney" "RR-3" "163" "_GLbulkRNAseq"

# eye
download_pair "eye" "RR-1" "100" ""

# skin
download_pair "skin" "RR-7" "254" ""
download_pair "skin" "RR-6" "243" "_GLbulkRNAseq"

# gastrocnemius
download_pair "gastrocnemius" "RR-1" "101" "_GLbulkRNAseq"
download_pair "gastrocnemius" "RR-9" "326" "_GLbulkRNAseq"

echo ""
echo "=== Download complete ==="
find "${BASE}" -name "*.csv" -newer "${BASE}" | wc -l
echo "new CSV files"
