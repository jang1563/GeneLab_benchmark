#!/usr/bin/env python3
"""
fetch_osdr.py — GeneLab_benchmark Phase 1: Data Download

Downloads GeneLab standard processed files (normalized counts, sample tables,
differential expression) for verified OSDR studies.

Priority order (Phase 1 first):
  1. Liver (A1, B1): OSD-48, 137, 245, 379, 242, 686
  2. Gastrocnemius (A2): OSD-101, 401, 326
  3. Kidney (A3): OSD-102, 163, 253
  4. Thymus (A4): OSD-244, 289, 421
  5. Eye (A6): OSD-100, 194, 397
  6. Skin (A5, supplementary): OSD-238, 239, 243

File targets per study (GeneLab standard pipeline outputs):
  - *_Normalized_Counts_GLbulkRNAseq.csv  (or *_Normalized_Counts.csv)
  - *_SampleTable_GLbulkRNAseq.csv        (sample metadata)
  - *_differential_expression_GLbulkRNAseq.csv (DESeq2 results)
  - *_contrasts_GLbulkRNAseq.csv          (contrast definitions)
  - *_bulkRNASeq_v*_runsheet.csv          (full metadata)
  - *_RSEM_Unnormalized_Counts*.csv        (raw counts)

Download method: GEODE endpoint (from SpaceOmicsBench fetch_genelab.py patterns)
  URL: https://osdr.nasa.gov/geode-py/ws/studies/{OSD_ID}/download
  Params: source=datamanager, file={filename}

Usage:
  python fetch_osdr.py                          # download all Phase 1 studies
  python fetch_osdr.py --tissue liver           # liver only
  python fetch_osdr.py --osd OSD-48 OSD-137    # specific studies
  python fetch_osdr.py --list                   # list what would be downloaded
  python fetch_osdr.py --check                  # check which files already exist
"""

import json
import time
import argparse
import requests
from pathlib import Path
from datetime import datetime

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data" / "mouse"
GLDS_VERIFIED_JSON = BASE_DIR / "GLDS_verified.json"
GEODE_DOWNLOAD = "https://osdr.nasa.gov/geode-py/ws/studies"
BIODATA_API = "https://visualization.osdr.nasa.gov/biodata/api/v2"

# ── Download targets per study ─────────────────────────────────────────────────
# Ordered by priority. The script tries each pattern and downloads what exists.
FILE_PRIORITY_PATTERNS = [
    # GeneLab v2+ standard pipeline outputs (preferred)
    "_rna_seq_Normalized_Counts_GLbulkRNAseq.csv",
    "_rna_seq_Normalized_Counts_rRNArm_GLbulkRNAseq.csv",
    "_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv",
    "_rna_seq_SampleTable_GLbulkRNAseq.csv",
    "_rna_seq_differential_expression_GLbulkRNAseq.csv",
    "_rna_seq_contrasts_GLbulkRNAseq.csv",
    "_rna_seq_bulkRNASeq_v2_runsheet.csv",
    "_rna_seq_bulkRNASeq_v1_runsheet.csv",
    # Older GeneLab pipeline outputs (v1)
    "_rna_seq_Normalized_Counts.csv",
    "_rna_seq_Unnormalized_Counts.csv",
    "_rna_seq_SampleTable.csv",
    "_rna_seq_contrasts.csv",
    "_rna_seq_differential_expression.csv",
]

# ── Study download manifest ────────────────────────────────────────────────────
# Maps OSD ID → tissue, mission, GLDS prefix (for filename construction)
DOWNLOAD_MANIFEST = {
    # ── Liver: Phase 1 priority ────────────────────────────────────────────────
    "OSD-48": {
        "tissue": "liver", "mission": "RR-1", "glds_prefix": "GLDS-48",
        "phase": 1, "task": "A1",
    },
    "OSD-137": {
        "tissue": "liver", "mission": "RR-3", "glds_prefix": "GLDS-137",
        "phase": 1, "task": "A1",
    },
    "OSD-245": {
        "tissue": "liver", "mission": "RR-6", "glds_prefix": "GLDS-245",
        "phase": 1, "task": "A1",
    },
    "OSD-379": {
        "tissue": "liver", "mission": "RR-8", "glds_prefix": "GLDS-379",
        "phase": 1, "task": "A1",
    },
    "OSD-242": {
        "tissue": "liver", "mission": "RR-9", "glds_prefix": "GLDS-242",
        "phase": 1, "task": "A1",
    },
    "OSD-686": {
        "tissue": "liver", "mission": "MHU-2", "glds_prefix": "GLDS-617",
        "phase": 1, "task": "A1",
        "notes": "GLDS-617 RNA-seq files. Has 3 groups: uG, 1G-centrifuge, GC",
    },
    # ── Gastrocnemius ─────────────────────────────────────────────────────────
    "OSD-101": {
        "tissue": "gastrocnemius", "mission": "RR-1", "glds_prefix": "GLDS-101",
        "phase": 2, "task": "A2",
    },
    "OSD-401": {
        "tissue": "gastrocnemius", "mission": "RR-5", "glds_prefix": "GLDS-401",
        "phase": 2, "task": "A2",
    },
    "OSD-326": {
        "tissue": "gastrocnemius", "mission": "RR-9", "glds_prefix": "GLDS-326",
        "phase": 2, "task": "A2",
    },
    # ── Kidney ────────────────────────────────────────────────────────────────
    "OSD-102": {
        "tissue": "kidney", "mission": "RR-1", "glds_prefix": "GLDS-102",
        "phase": 2, "task": "A3",
    },
    "OSD-163": {
        "tissue": "kidney", "mission": "RR-3", "glds_prefix": "GLDS-163",
        "phase": 2, "task": "A3",
    },
    "OSD-253": {
        "tissue": "kidney", "mission": "RR-7", "glds_prefix": "GLDS-253",
        "phase": 2, "task": "A3",
    },
    # ── Thymus ────────────────────────────────────────────────────────────────
    "OSD-244": {
        "tissue": "thymus", "mission": "RR-6", "glds_prefix": "GLDS-244",
        "phase": 2, "task": "A4",
    },
    "OSD-289": {
        "tissue": "thymus", "mission": "MHU-2", "glds_prefix": "GLDS-289",
        "phase": 2, "task": "A4",
    },
    "OSD-421": {
        "tissue": "thymus", "mission": "RR-9", "glds_prefix": "GLDS-421",
        "phase": 2, "task": "A4",
    },
    # ── Eye / Retina ──────────────────────────────────────────────────────────
    "OSD-100": {
        "tissue": "eye", "mission": "RR-1", "glds_prefix": "GLDS-100",
        "phase": 2, "task": "A6",
    },
    "OSD-194": {
        "tissue": "eye", "mission": "RR-3", "glds_prefix": "GLDS-194",
        "phase": 2, "task": "A6",
    },
    "OSD-397": {
        "tissue": "eye", "mission": "TBD", "glds_prefix": "GLDS-397",
        "phase": 2, "task": "A6",
        "notes": "Mission TBD. Normalized counts available.",
    },
    # ── Skin (supplementary) ─────────────────────────────────────────────────
    "OSD-238": {
        "tissue": "skin", "mission": "MHU-2 (dorsal)", "glds_prefix": "GLDS-238",
        "phase": 3, "task": "A5",
    },
    "OSD-239": {
        "tissue": "skin", "mission": "MHU-2 (femoral)", "glds_prefix": "GLDS-239",
        "phase": 3, "task": "A5",
    },
    "OSD-243": {
        "tissue": "skin", "mission": "RR-6", "glds_prefix": "GLDS-243",
        "phase": 3, "task": "A5",
    },
    "OSD-254": {
        "tissue": "skin", "mission": "RR-7", "glds_prefix": "GLDS-254",
        "phase": 3, "task": "A5",
        "notes": "Mixed strain (C57BL/6J + C3H/HeJ). Use C57BL/6J subset only for Track 2a.",
    },
    # ── HU analog ────────────────────────────────────────────────────────────
    "OSD-295": {
        "tissue": "soleus", "mission": "HU", "glds_prefix": "GLDS-295",
        "phase": 3, "task": "B6",
        "notes": "Hindlimb Unloading RNA-seq. B6 candidate.",
    },
    # ── Category J (pipeline comparison) — needs GLDS-168 ────────────────────
    "OSD-168": {
        "tissue": "liver", "mission": "RR-1+3 (merged)", "glds_prefix": "GLDS-168",
        "phase": 3, "task": "J1",
        "notes": "Category J only. Duplicate of OSD-48+137. Pipeline version comparison.",
    },
}


# ── API helpers ────────────────────────────────────────────────────────────────

def list_files(osd_id: str) -> dict:
    """List all files for an OSD study via Biodata API."""
    url = f"{BIODATA_API}/dataset/{osd_id}/files/"
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        return data.get(osd_id, data).get("files", {})
    except Exception as e:
        print(f"    [WARN] Could not list files for {osd_id}: {e}")
        return {}


def download_file(osd_id: str, filename: str, outpath: Path,
                  timeout: int = 300, chunk_size: int = 8192) -> bool:
    """
    Download a single file via GEODE endpoint.
    Follows S3 presigned URL redirect automatically.
    """
    url = f"{GEODE_DOWNLOAD}/{osd_id}/download"
    params = {"source": "datamanager", "file": filename}

    try:
        resp = requests.get(url, params=params, stream=True,
                            timeout=timeout, allow_redirects=True)
        resp.raise_for_status()

        outpath.parent.mkdir(parents=True, exist_ok=True)
        total_bytes = 0
        with open(outpath, "wb") as f:
            for chunk in resp.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
                    total_bytes += len(chunk)

        size_kb = total_bytes / 1024
        print(f"      ✓ {filename} ({size_kb:.0f} KB)")
        return True

    except requests.exceptions.HTTPError as e:
        if e.response.status_code == 404:
            return False  # File not found — try next pattern
        print(f"      ✗ HTTP {e.response.status_code}: {filename}")
        return False
    except Exception as e:
        print(f"      ✗ Error downloading {filename}: {e}")
        return False


# ── Core download logic ────────────────────────────────────────────────────────

def get_output_dir(info: dict) -> Path:
    """Determine output directory: data/mouse/{tissue}/{mission}/"""
    tissue = info["tissue"].replace(" ", "_").replace("/", "_")
    mission = info["mission"].replace(" ", "_").replace("/", "_").replace("+", "_")
    return DATA_DIR / tissue / mission


def download_study(osd_id: str, info: dict, dry_run: bool = False,
                   force: bool = False) -> dict:
    """
    Download all processable files for a single study.
    Returns summary dict with counts.
    """
    glds_prefix = info["glds_prefix"]
    outdir = get_output_dir(info)
    osd_num = osd_id.replace("OSD-", "")
    notes = info.get("notes", "")

    print(f"\n  [{osd_id}] {info['tissue']} / {info['mission']} (task: {info['task']})")
    if notes:
        print(f"    Note: {notes}")

    if dry_run:
        print(f"    → would download to: {outdir}")
        return {"osd_id": osd_id, "downloaded": 0, "skipped": 0, "failed": 0}

    # Get file listing from API to discover actual filenames
    all_files = list_files(osd_id)
    if not all_files:
        print(f"    [WARN] No files found for {osd_id}")
        return {"osd_id": osd_id, "downloaded": 0, "skipped": 0, "failed": 1}

    downloaded = 0
    skipped = 0
    failed = 0

    # Try each priority pattern
    for pattern in FILE_PRIORITY_PATTERNS:
        # Build candidate filename
        target_filename = f"{glds_prefix}{pattern}"

        # Check if file exists in study
        if target_filename not in all_files:
            continue

        outpath = outdir / target_filename
        if outpath.exists() and not force:
            print(f"      ↷ {target_filename} (already exists)")
            skipped += 1
            continue

        success = download_file(osd_id, target_filename, outpath)
        if success:
            downloaded += 1
        else:
            failed += 1

        time.sleep(0.2)  # Brief pause between files

    if downloaded == 0 and skipped == 0:
        print(f"    [WARN] No matching processed files found for {glds_prefix}")
        # Fallback: try to find normalized counts from full file listing
        norm_candidates = [
            f for f in all_files
            if "normalized" in f.lower() and f.endswith(".csv")
        ]
        if norm_candidates:
            print(f"    Fallback: trying {len(norm_candidates)} normalized count file(s)...")
            for fname in norm_candidates[:3]:
                outpath = outdir / fname
                if not outpath.exists() or force:
                    success = download_file(osd_id, fname, outpath)
                    if success:
                        downloaded += 1
                time.sleep(0.2)

    return {
        "osd_id": osd_id,
        "tissue": info["tissue"],
        "mission": info["mission"],
        "task": info["task"],
        "outdir": str(outdir),
        "downloaded": downloaded,
        "skipped": skipped,
        "failed": failed,
    }


def write_download_log(results: list[dict]) -> None:
    """Write download summary to data/DOWNLOAD_LOG.json."""
    log_path = DATA_DIR.parent / "DOWNLOAD_LOG.json"
    log = {
        "generated_at": datetime.now().isoformat(),
        "results": results,
        "summary": {
            "total_downloaded": sum(r["downloaded"] for r in results),
            "total_skipped": sum(r["skipped"] for r in results),
            "total_failed": sum(r["failed"] for r in results),
        },
    }
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(json.dumps(log, indent=2))
    print(f"\n✓ Download log: {log_path}")


# ── CLI ────────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="Download OSDR processed RNA-seq files for GeneLab_benchmark"
    )
    parser.add_argument(
        "--osd", nargs="+", default=None,
        help="Specific OSD IDs to download (e.g. OSD-48 OSD-137)"
    )
    parser.add_argument(
        "--tissue", type=str, default=None,
        help="Filter by tissue (e.g. liver, kidney, thymus)"
    )
    parser.add_argument(
        "--phase", type=int, default=None,
        help="Download only studies from a specific phase (1, 2, or 3)"
    )
    parser.add_argument(
        "--task", type=str, default=None,
        help="Download only studies for a specific task (e.g. A1, A2)"
    )
    parser.add_argument(
        "--list", action="store_true",
        help="List what would be downloaded without downloading"
    )
    parser.add_argument(
        "--check", action="store_true",
        help="Check which files already exist locally"
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Re-download files even if they already exist"
    )
    parser.add_argument(
        "--delay", type=float, default=1.0,
        help="Delay in seconds between study downloads. Default: 1.0"
    )
    return parser.parse_args()


def filter_manifest(args) -> dict:
    """Apply CLI filters to DOWNLOAD_MANIFEST."""
    manifest = DOWNLOAD_MANIFEST

    if args.osd:
        manifest = {k: v for k, v in manifest.items() if k in args.osd}
    if args.tissue:
        manifest = {k: v for k, v in manifest.items()
                    if v["tissue"] == args.tissue}
    if args.phase:
        manifest = {k: v for k, v in manifest.items()
                    if v.get("phase") == args.phase}
    if args.task:
        manifest = {k: v for k, v in manifest.items()
                    if v.get("task") == args.task}

    return manifest


def check_existing(manifest: dict) -> None:
    """Report which files already exist locally."""
    print("\n=== Existing Files Check ===")
    for osd_id, info in manifest.items():
        outdir = get_output_dir(info)
        glds_prefix = info["glds_prefix"]
        existing = list(outdir.glob(f"{glds_prefix}*.csv")) if outdir.exists() else []
        status = f"✓ {len(existing)} files" if existing else "✗ none"
        print(f"  {osd_id} ({info['tissue']}/{info['mission']}): {status}")
        for f in existing:
            print(f"    {f.name}")


def main():
    args = parse_args()
    manifest = filter_manifest(args)

    if not manifest:
        print("No studies match the specified filters.")
        return

    print("=" * 60)
    print("GeneLab_benchmark — Phase 1 Data Download")
    print(f"Studies to process: {len(manifest)}")
    print(f"Output directory:   {DATA_DIR}")
    print("=" * 60)

    if args.check:
        check_existing(manifest)
        return

    if args.list:
        print("\nWould download:")
        for osd_id, info in manifest.items():
            outdir = get_output_dir(info)
            print(f"  {osd_id} ({info['tissue']} / {info['mission']}) → {outdir}")
        return

    # Sort by phase (Phase 1 first)
    sorted_manifest = sorted(manifest.items(), key=lambda x: x[1].get("phase", 99))

    results = []
    for osd_id, info in sorted_manifest:
        result = download_study(osd_id, info, dry_run=False, force=args.force)
        results.append(result)
        time.sleep(args.delay)

    write_download_log(results)

    # Final summary
    total_dl = sum(r["downloaded"] for r in results)
    total_sk = sum(r["skipped"] for r in results)
    total_fa = sum(r["failed"] for r in results)

    print("\n" + "=" * 60)
    print(f"✓ Downloaded: {total_dl} files")
    print(f"↷ Skipped (exist): {total_sk} files")
    print(f"✗ Failed: {total_fa} files")
    print("=" * 60)

    if total_fa > 0:
        print("\nFailed studies — check manually:")
        for r in results:
            if r["failed"] > 0:
                print(f"  {r['osd_id']} ({r['tissue']}/{r['mission']})")


if __name__ == "__main__":
    main()
