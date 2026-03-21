#!/usr/bin/env python3
"""Download mouse tissue normalized counts + sample tables from OSDR for E4.

Run on Cayuga:
  python3 download_mouse_e4.py
"""
import json
import os
import urllib.request
from pathlib import Path

BASE = Path("/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark/data/mouse")
OSDR_API = "https://osdr.nasa.gov/osdr/data/osd/files/{osd_num}"
OSDR_DL = "https://osdr.nasa.gov"

# tissue → mission → (glds_num, osd_num, count_suffix, sample_suffix)
# osd_num may differ from glds_num
DATASETS = {
    "liver": {
        "RR-1": ("GLDS-48", 48, "_GLbulkRNAseq"),
        "RR-3": ("GLDS-137", 137, ""),
    },
    "thymus": {
        "MHU-2": ("GLDS-289", 289, "_GLbulkRNAseq"),
        "RR-6": ("GLDS-244", 244, "_GLbulkRNAseq"),
    },
    "kidney": {
        "RR-1": ("GLDS-102", 102, "_GLbulkRNAseq"),
        "RR-3": ("GLDS-163", 163, "_GLbulkRNAseq"),
    },
    "eye": {
        "RR-1": ("GLDS-100", 100, ""),
    },
    "skin": {
        "RR-7": ("GLDS-254", 254, ""),
        "RR-6": ("GLDS-243", 243, "_GLbulkRNAseq"),
    },
    "gastrocnemius": {
        "RR-1": ("GLDS-101", 101, "_GLbulkRNAseq"),
        "RR-9": ("GLDS-326", 326, "_GLbulkRNAseq"),
    },
}


def get_download_url(osd_num: int, filename: str) -> str:
    """Query OSDR API to get the actual download URL for a file."""
    api_url = OSDR_API.format(osd_num=osd_num)
    try:
        with urllib.request.urlopen(api_url, timeout=30) as resp:
            data = json.loads(resp.read())
    except Exception as e:
        print(f"  API error for OSD-{osd_num}: {e}")
        return ""

    studies = data.get("studies", {})
    for sid, sdata in studies.items():
        for f in sdata.get("study_files", []):
            if f.get("file_name") == filename:
                rel_url = f.get("remote_url", "")
                if rel_url:
                    return OSDR_DL + rel_url
    return ""


def download_file(url: str, dest: Path) -> bool:
    """Download a file from URL to dest."""
    try:
        urllib.request.urlretrieve(url, str(dest))
        size = dest.stat().st_size
        print(f"  OK: {dest.name} ({size / 1e6:.1f} MB)")
        return True
    except Exception as e:
        print(f"  FAIL: {dest.name} — {e}")
        if dest.exists():
            dest.unlink()
        return False


def main():
    print("=== Downloading mouse tissue data for E4 ===\n")
    total_ok = 0
    total_fail = 0

    for tissue, missions in DATASETS.items():
        for mission, (glds, osd_num, suffix) in missions.items():
            out_dir = BASE / tissue / mission
            out_dir.mkdir(parents=True, exist_ok=True)

            # Normalized counts
            nc_name = f"{glds}_rna_seq_Normalized_Counts{suffix}.csv"
            nc_path = out_dir / nc_name
            if nc_path.exists() and nc_path.stat().st_size > 1000:
                print(f"EXISTS: {tissue}/{mission}/{nc_name}")
                total_ok += 1
            else:
                print(f"Downloading {tissue}/{mission}/{nc_name}...")
                url = get_download_url(osd_num, nc_name)
                if url:
                    if download_file(url, nc_path):
                        total_ok += 1
                    else:
                        total_fail += 1
                else:
                    print(f"  URL not found for {nc_name}")
                    total_fail += 1

            # Sample table
            st_name = f"{glds}_rna_seq_SampleTable{suffix}.csv"
            st_path = out_dir / st_name
            if st_path.exists() and st_path.stat().st_size > 100:
                print(f"EXISTS: {tissue}/{mission}/{st_name}")
                total_ok += 1
            else:
                print(f"Downloading {tissue}/{mission}/{st_name}...")
                url = get_download_url(osd_num, st_name)
                if url:
                    if download_file(url, st_path):
                        total_ok += 1
                    else:
                        total_fail += 1
                else:
                    print(f"  URL not found for {st_name}")
                    total_fail += 1

    print(f"\n=== Done: {total_ok} OK, {total_fail} failed ===")


if __name__ == "__main__":
    main()
