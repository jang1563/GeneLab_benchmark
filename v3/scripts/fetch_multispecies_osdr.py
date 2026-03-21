#!/usr/bin/env python3
"""Download multi-species spaceflight RNA-seq data from OSDR for Phase 1 (E4/E5).

Datasets:
  - OSD-207: Drosophila melanogaster, whole body (female), FLT vs GC, 4 genotypes
  - OSD-37: Arabidopsis thaliana, seedling pool, FLT vs GC, 4 ecotypes
  - OSD-120: Arabidopsis thaliana, root, FLT vs GC × Light, 3 genotypes
"""
import os
import json
import urllib.request
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent.parent
DATA_DIR = BASE_DIR / "data" / "multispecies"

# Files to download per dataset
DATASETS = {
    "OSD-207": {
        "species": "drosophila",
        "files": [
            "GLDS-207_rna_seq_Normalized_Counts_GLbulkRNAseq.csv",
            "GLDS-207_rna_seq_SampleTable_GLbulkRNAseq.csv",
        ],
    },
    "OSD-37": {
        "species": "arabidopsis",
        "files": [
            "GLDS-37_rna_seq_Normalized_Counts_GLbulkRNAseq.csv",
            "GLDS-37_rna_seq_SampleTable_GLbulkRNAseq.csv",
        ],
    },
    "OSD-120": {
        "species": "arabidopsis",
        "files": [
            "GLDS-120_rna_seq_Normalized_Counts_GLbulkRNAseq.csv",
            "GLDS-120_rna_seq_SampleTable_GLbulkRNAseq.csv",
        ],
    },
}


def get_download_urls(osd_id: str) -> dict:
    """Get file download URLs from OSDR API."""
    num = osd_id.split("-")[1]
    api_url = f"https://osdr.nasa.gov/osdr/data/osd/files/{num}"
    print(f"  Fetching file list from {api_url}...")

    req = urllib.request.Request(api_url)
    with urllib.request.urlopen(req, timeout=30) as resp:
        data = json.loads(resp.read().decode())

    study_files = data["studies"][osd_id]["study_files"]
    url_map = {}
    base_url = "https://osdr.nasa.gov"
    for sf in study_files:
        url = sf["remote_url"]
        # OSDR API returns relative URLs — prepend base
        if url.startswith("/"):
            url = base_url + url
        url_map[sf["file_name"]] = url
    return url_map


def download_file(url: str, dest: Path):
    """Download a file from URL to destination."""
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  Already exists: {dest.name} ({dest.stat().st_size / 1e6:.1f} MB)")
        return

    print(f"  Downloading {dest.name}...")
    urllib.request.urlretrieve(url, str(dest))
    print(f"  Saved: {dest.name} ({dest.stat().st_size / 1e6:.1f} MB)")


def main():
    for osd_id, info in DATASETS.items():
        species = info["species"]
        out_dir = DATA_DIR / species
        out_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n=== {osd_id} ({species}) ===")

        try:
            url_map = get_download_urls(osd_id)
        except Exception as e:
            print(f"  ERROR fetching URLs: {e}")
            continue

        for fname in info["files"]:
            if fname not in url_map:
                print(f"  WARNING: {fname} not found in OSDR API response")
                # Try with different naming patterns
                continue

            dest = out_dir / fname
            try:
                download_file(url_map[fname], dest)
            except Exception as e:
                print(f"  ERROR downloading {fname}: {e}")

    print("\n=== Download complete ===")
    # Summary
    for osd_id, info in DATASETS.items():
        species = info["species"]
        out_dir = DATA_DIR / species
        files = list(out_dir.glob(f"GLDS-{osd_id.split('-')[1]}*"))
        print(f"  {osd_id}: {len(files)} files in {out_dir}")


if __name__ == "__main__":
    main()
