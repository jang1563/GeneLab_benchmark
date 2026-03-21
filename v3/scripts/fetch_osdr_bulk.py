#!/usr/bin/env python3
"""
fetch_osdr_bulk.py — OSDR에서 GeneLab processed bulk RNA-seq 파일 다운로드.

Downloads Normalized_Counts + SampleTable CSVs for Phase 4 & 5 datasets.

Usage:
    python3 v3/scripts/fetch_osdr_bulk.py                  # Download all
    python3 v3/scripts/fetch_osdr_bulk.py --phase 4        # Phase 4 only
    python3 v3/scripts/fetch_osdr_bulk.py --osd 248        # Single dataset
    python3 v3/scripts/fetch_osdr_bulk.py --dry-run        # Preview only

Output:
    data/mouse/{tissue}/{mission}/GLDS-{id}_rna_seq_*.csv
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path
from urllib.request import urlopen, Request, urlretrieve
from urllib.error import URLError, HTTPError

BASE_DIR = Path(__file__).resolve().parent.parent.parent  # GeneLab_benchmark/
DATA_DIR = BASE_DIR / "data" / "mouse"
OSDR_FILES_API = "https://osdr.nasa.gov/osdr/data/osd/files/{osd_id}"
OSDR_DOWNLOAD_BASE = "https://osdr.nasa.gov"

# Files to download per dataset
TARGET_FILE_PATTERNS = [
    "Normalized_Counts_GLbulkRNAseq.csv",
    "Normalized_Counts_rRNArm_GLbulkRNAseq.csv",
    "SampleTable_GLbulkRNAseq.csv",
    "RSEM_Unnormalized_Counts_GLbulkRNAseq.csv",
    "bulkRNASeq_v2_runsheet.csv",
    "bulkRNASeq_v1_runsheet.csv",
    # Fallbacks for older GeneLab pipeline versions
    "Normalized_Counts.csv",
    "SampleTable.csv",
]

# Dataset → tissue/mission mapping
DATASETS = {
    # Phase 4: New bulk RNA-seq
    248: {"tissue": "lung",  "mission": "RR-6",  "phase": 4},
    247: {"tissue": "colon", "mission": "RR-6",  "phase": 4},
    253: {"tissue": "kidney", "mission": "RR-7", "phase": 4},  # Already downloaded
    240: {"tissue": "skin",  "mission": "RR-5_dorsal",  "phase": 4},
    241: {"tissue": "skin",  "mission": "RR-5_femoral", "phase": 4},
    243: {"tissue": "skin",  "mission": "RR-6_dorsal",  "phase": 4},
    # Phase 5: Radiation
    202: {"tissue": "retina_brain", "mission": "LDR_HLU", "phase": 5},
    211: {"tissue": "spleen",      "mission": "LDR_HLU", "phase": 5},
    237: {"tissue": "skin_rad",    "mission": "LDR_HLU", "phase": 5},
}


def fetch_file_list(osd_id):
    """Fetch file list from OSDR files API."""
    url = OSDR_FILES_API.format(osd_id=osd_id)
    req = Request(url, headers={"User-Agent": "GeneLabBench/3.0"})
    try:
        with urlopen(req, timeout=30) as resp:
            raw = json.loads(resp.read().decode("utf-8"))
        key = f"OSD-{osd_id}"
        if key not in raw.get("studies", {}):
            return []
        return raw["studies"][key].get("study_files", [])
    except (URLError, HTTPError) as e:
        print(f"  ERROR fetching file list: {e}")
        return []


def find_target_files(file_list):
    """Find files matching our target patterns."""
    matches = []
    for f in file_list:
        fname = f.get("file_name", "")
        for pattern in TARGET_FILE_PATTERNS:
            if fname.endswith(pattern) or pattern in fname:
                matches.append({
                    "file_name": fname,
                    "file_size": f.get("file_size", 0),
                    "remote_url": f.get("remote_url", ""),
                })
                break
    return matches


def download_file(remote_url, output_path):
    """Download a file from OSDR."""
    full_url = OSDR_DOWNLOAD_BASE + remote_url
    output_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        req = Request(full_url, headers={"User-Agent": "GeneLabBench/3.0"})
        with urlopen(req, timeout=300) as resp:
            with open(output_path, "wb") as out:
                while True:
                    chunk = resp.read(65536)
                    if not chunk:
                        break
                    out.write(chunk)
        return True
    except (URLError, HTTPError) as e:
        print(f"    DOWNLOAD ERROR: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(description="Download OSDR bulk RNA-seq files")
    parser.add_argument("--phase", type=int, help="Download only this phase (4 or 5)")
    parser.add_argument("--osd", type=int, help="Download only this OSD ID")
    parser.add_argument("--dry-run", action="store_true", help="Preview only")
    parser.add_argument("--force", action="store_true", help="Re-download existing files")
    args = parser.parse_args()

    print("=" * 70)
    print("GeneLabBench v3 — OSDR Bulk RNA-seq Download")
    print("=" * 70)

    datasets = DATASETS
    if args.osd:
        datasets = {args.osd: DATASETS[args.osd]}
    elif args.phase:
        datasets = {k: v for k, v in DATASETS.items() if v["phase"] == args.phase}

    total_downloaded = 0
    total_skipped = 0
    total_failed = 0

    for osd_id, info in datasets.items():
        tissue = info["tissue"]
        mission = info["mission"]
        phase = info["phase"]
        output_dir = DATA_DIR / tissue / mission

        print(f"\n--- OSD-{osd_id} → {tissue}/{mission} (Phase {phase}) ---")

        # Check if already downloaded
        if output_dir.exists() and not args.force:
            existing = list(output_dir.glob("*.csv"))
            if existing:
                print(f"  Already exists ({len(existing)} CSV files). Use --force to re-download.")
                total_skipped += 1
                continue

        # Fetch file list
        file_list = fetch_file_list(osd_id)
        if not file_list:
            print(f"  No files found on OSDR")
            total_failed += 1
            continue

        # Find target files
        targets = find_target_files(file_list)
        if not targets:
            print(f"  No matching processed files found in {len(file_list)} total files")
            # Show available categories for debugging
            categories = set(f.get("category", "") for f in file_list)
            print(f"  Available categories: {categories}")
            total_failed += 1
            continue

        print(f"  Found {len(targets)} target files:")
        for t in targets:
            size_mb = t["file_size"] / 1_000_000
            print(f"    {t['file_name']} ({size_mb:.1f} MB)")

        if args.dry_run:
            continue

        # Download each file
        for t in targets:
            output_path = output_dir / t["file_name"]
            if output_path.exists() and not args.force:
                print(f"    SKIP (exists): {t['file_name']}")
                continue

            print(f"    Downloading: {t['file_name']}...", end=" ", flush=True)
            if download_file(t["remote_url"], output_path):
                actual_size = output_path.stat().st_size
                print(f"OK ({actual_size / 1_000_000:.1f} MB)")
                total_downloaded += 1
            else:
                total_failed += 1

        time.sleep(0.5)  # Rate limit

    print(f"\n{'=' * 70}")
    print(f"SUMMARY: Downloaded={total_downloaded}, Skipped={total_skipped}, Failed={total_failed}")


if __name__ == "__main__":
    main()
