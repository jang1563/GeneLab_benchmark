#!/usr/bin/env python3
"""
catalog_datasets.py — GeneLab_benchmark Phase 1: GLDS ID Verification

Queries the NASA OSDR Biodata API to:
  1. Verify each candidate GLDS ID exists
  2. Confirm bulk RNA-seq files are present (not microarray/spatial/scRNA)
  3. Extract control type (GC/VC/BC) from metadata if available
  4. Flag assay type (bulk RNA-seq vs spatial vs single-cell vs microarray)
  5. Search for Hindlimb Unloading (HU) studies

Output:
  - DATA_CATALOG.md  (human-readable table)
  - GLDS_verified.json (machine-readable, used by fetch_osdr.py)

Usage:
  python catalog_datasets.py               # verify all candidate GLDS IDs
  python catalog_datasets.py --search-hu   # also search for HU analog studies
  python catalog_datasets.py --osd OSD-48  # check single study (OSD or GLDS prefix)

API reference:
  Biodata API: https://visualization.osdr.nasa.gov/biodata/api/v2/
  Based on patterns from SpaceOmicsBench/v2_public/scripts/fetch_genelab.py

Design decisions applied:
  - GLDS-168 flagged as DUPLICATE (see DESIGN_DECISIONS.md DD-05)
  - Assay type filtering: bulk RNA-seq only for Category A
  - BALB/c strain flag for Track 2a/2b separation (DD-06)
"""

import json
import time
import argparse
import requests
from pathlib import Path
from datetime import datetime

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent.parent
DATA_CATALOG_MD = BASE_DIR / "DATA_CATALOG.md"
GLDS_VERIFIED_JSON = BASE_DIR / "GLDS_verified.json"

# ── API endpoints (from SpaceOmicsBench verified patterns) ─────────────────────
BIODATA_API = "https://visualization.osdr.nasa.gov/biodata/api/v2"
OSDR_SEARCH_API = "https://osdr.nasa.gov/osdr/data/osd/meta"
OSDR_FILES_API = "https://osdr.nasa.gov/osdr/data/osd/files"

# ── Candidate GLDS IDs from PLAN.md v0.5 ──────────────────────────────────────
# Status legend: "primary" = included in main analysis, "duplicate" = excluded
# Track: "2a" = C57BL/6J only, "2b" = all strains
CANDIDATE_DATASETS = [
    # ── Liver (A1, B1) ────────────────────────────────────────────────────────
    {"osd_id": "OSD-48",  "tissue": "liver",         "mission": "RR-1",  "track": "2a", "notes": "primary"},
    {"osd_id": "OSD-137", "tissue": "liver",         "mission": "RR-3",  "track": "2a", "notes": "primary"},
    {"osd_id": "OSD-168", "tissue": "liver",         "mission": "RR-1+3","track": "2a", "notes": "DUPLICATE of OSD-48+137 — Category J only (DESIGN_DECISIONS DD-05)"},
    {"osd_id": "OSD-245", "tissue": "liver",         "mission": "RR-6",  "track": "2a", "notes": "primary"},
    {"osd_id": "OSD-379", "tissue": "liver",         "mission": "RR-8",  "track": "2a", "notes": "primary"},
    {"osd_id": "OSD-242", "tissue": "liver",         "mission": "RR-9",  "track": "2a", "notes": "primary"},
    {"osd_id": "OSD-686", "tissue": "liver",         "mission": "MHU-2", "track": "2a", "notes": "primary — GLDS-617 RNA-seq files (uG vs GC vs 1G-centrifuge, n=3/group)"},
    {"osd_id": "OSD-617", "tissue": "liver",         "mission": "MHU-2", "track": "exc","notes": "EXCLUDE — cytology/estrous cycle only, NOT RNA-seq. RNA-seq is in OSD-686"},
    {"osd_id": "OSD-173", "tissue": "liver",         "mission": "STS-135","track": "qc","notes": "C57BL/6CR (C57BL/6 subline), n=2/group — FAILS QC minimum (n<3/group)"},
    {"osd_id": "OSD-25",  "tissue": "liver",         "mission": "STS-135","track": "exc","notes": "MICROARRAY confirmed — not RNA-seq. Completely excluded"},
    # ── Gastrocnemius (A2, B4) ────────────────────────────────────────────────
    {"osd_id": "OSD-101", "tissue": "gastrocnemius", "mission": "RR-1",  "track": "2a", "notes": "primary"},
    {"osd_id": "OSD-326", "tissue": "gastrocnemius", "mission": "RR-9",  "track": "2a", "notes": "primary"},
    {"osd_id": "OSD-401", "tissue": "gastrocnemius", "mission": "RR-5",  "track": "2a", "notes": "primary"},
    # ── Soleus (supplementary) ────────────────────────────────────────────────
    {"osd_id": "OSD-104", "tissue": "soleus",        "mission": "RR-1",  "track": "2a", "notes": "supplementary — check mission count"},
    {"osd_id": "OSD-638", "tissue": "soleus",        "mission": "MHU-8", "track": "2a", "notes": "supplementary — unverified"},
    # ── Kidney (A3, B5) ───────────────────────────────────────────────────────
    {"osd_id": "OSD-102", "tissue": "kidney",        "mission": "RR-1",  "track": "2a", "notes": "primary"},
    {"osd_id": "OSD-163", "tissue": "kidney",        "mission": "RR-3",  "track": "2a", "notes": "primary"},
    {"osd_id": "OSD-253", "tissue": "kidney",        "mission": "RR-7",  "track": "2a", "notes": "primary"},
    {"osd_id": "OSD-674", "tissue": "kidney",        "mission": "?",     "track": "2a", "notes": "unverified — mission TBD"},
    # ── Thymus (A4) ───────────────────────────────────────────────────────────
    {"osd_id": "OSD-4",   "tissue": "thymus",        "mission": "STS-118","track": "2a","notes": "shuttle mission — verify assay type"},
    {"osd_id": "OSD-244", "tissue": "thymus",        "mission": "RR-6",  "track": "2a", "notes": "primary"},
    {"osd_id": "OSD-289", "tissue": "thymus",        "mission": "MHU-2", "track": "2a", "notes": "primary"},
    {"osd_id": "OSD-421", "tissue": "thymus",        "mission": "RR-9",  "track": "2a", "notes": "primary"},
    # ── Skin (A5) ─────────────────────────────────────────────────────────────
    {"osd_id": "OSD-238", "tissue": "skin",          "mission": "MHU-2", "track": "2a", "notes": "primary"},
    {"osd_id": "OSD-239", "tissue": "skin",          "mission": "MHU-2", "track": "2a", "notes": "primary — check if same study as OSD-238"},
    {"osd_id": "OSD-243", "tissue": "skin",          "mission": "RR-6",  "track": "2a", "notes": "primary"},
    {"osd_id": "OSD-689", "tissue": "skin",          "mission": "RR-8",  "track": "exc","notes": "EXCLUDE — scRNA-seq (GLDS-689), not bulk RNA-seq"},
    {"osd_id": "OSD-793", "tissue": "skin",          "mission": "RRRM-1","track": "exc","notes": "EXCLUDE — targetRNAseq (NanoString panel), not standard bulk mRNA"},
    # ── Eye / Retina (A6) ─────────────────────────────────────────────────────
    {"osd_id": "OSD-100", "tissue": "eye",           "mission": "RR-1",  "track": "2a", "notes": "primary"},
    {"osd_id": "OSD-194", "tissue": "eye",           "mission": "RR-3",  "track": "2a", "notes": "primary"},
    {"osd_id": "OSD-397", "tissue": "eye/retina",    "mission": "TBD",   "track": "2a", "notes": "bulk RNA-seq confirmed (n=9 FLT, n=7 GC). Mission TBD — verify ISA metadata"},
    {"osd_id": "OSD-664", "tissue": "eye",           "mission": "MHU-8", "track": "exc","notes": "EXCLUDE — western blot + immunofluorescence only, NOT RNA-seq"},
    {"osd_id": "OSD-87",  "tissue": "eye",           "mission": "?",     "track": "exc","notes": "EXCLUDE — microarray"},
    # ── HU Analog studies (for B6) ───────────────────────────────────────────
    {"osd_id": "OSD-295", "tissue": "soleus",        "mission": "HU",    "track": "HU", "notes": "Hindlimb Unloading soleus RNA-seq — B6 candidate"},
    {"osd_id": "OSD-335", "tissue": "liver",         "mission": "HU+HZE","track": "exc","notes": "EXCLUDE: miRNA-Seq (not mRNA bulk RNA-seq)"},
    # ── Known non-bulk (Category A exclusions — confirm) ─────────────────────
    {"osd_id": "OSD-270", "tissue": "heart",         "mission": "RR-3",  "track": "exc","notes": "EXCLUDE from A: Visium spatial transcriptomics"},
    {"osd_id": "OSD-596", "tissue": "heart",         "mission": "NG-11", "track": "exc","notes": "heart — verify assay type, ≤2 missions"},
]

# ── HU analog search terms ─────────────────────────────────────────────────────
HU_QUERY_TERMS = [
    "hindlimb unloading",
    "tail suspension",
    "antiorthostatic suspension",
]

# ── Assay type keywords for bulk RNA-seq detection ────────────────────────────
BULK_RNASEQ_KEYWORDS = [
    "normalized_counts", "normalized-counts",
    "rna-seq_Normalized", "rnaseq_normalized",
    "DESeq2", "edgeR", "limma",
    "_counts_", "_count_matrix",
    "Gene_Expression_Processed",
    "rna-seq_", "bulk_rna",
]
SPATIAL_KEYWORDS = ["visium", "spatial", "GeoMx", "geomx", "nanostring"]
SINGLE_CELL_KEYWORDS = ["scrna", "snrna", "10x", "10X", "single_cell", "singlecell"]
MICROARRAY_KEYWORDS = ["affymetrix", "agilent", "microarray", "array_"]


# ── API functions ──────────────────────────────────────────────────────────────

def get_osd_files(osd_id: str, retries: int = 3, delay: float = 1.0) -> dict:
    """
    List files for an OSD study using the Biodata API.
    Returns dict of {filename: file_info} or {} on failure.
    """
    url = f"{BIODATA_API}/dataset/{osd_id}/files/"
    for attempt in range(retries):
        try:
            resp = requests.get(url, timeout=30)
            resp.raise_for_status()
            data = resp.json()
            osd_data = data.get(osd_id, data)
            files = osd_data.get("files", {})
            if isinstance(files, dict):
                return files
            return {}
        except requests.exceptions.HTTPError as e:
            if e.response.status_code == 404:
                return {}  # Study does not exist
            if attempt < retries - 1:
                time.sleep(delay)
        except Exception:
            if attempt < retries - 1:
                time.sleep(delay)
    return {}


def get_osd_metadata(osd_id: str) -> dict:
    """
    Get study metadata (title, factors, organism, assay type) via OSDR API.
    Returns dict or {} on failure.
    """
    numeric_id = osd_id.replace("OSD-", "")
    url = f"{OSDR_SEARCH_API}/{numeric_id}"
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        studies = data.get("hits", {}).get("hits", [])
        if studies:
            return studies[0].get("_source", {})
        return {}
    except Exception:
        return {}


def detect_assay_type(filenames: list[str]) -> str:
    """
    Classify the assay type from file names.
    Returns: 'bulk_rnaseq' | 'spatial' | 'single_cell' | 'microarray' | 'unknown'
    """
    names_str = " ".join(filenames).lower()

    for kw in SPATIAL_KEYWORDS:
        if kw.lower() in names_str:
            return "spatial"
    for kw in SINGLE_CELL_KEYWORDS:
        if kw.lower() in names_str:
            return "single_cell"
    for kw in MICROARRAY_KEYWORDS:
        if kw.lower() in names_str:
            return "microarray"
    for kw in BULK_RNASEQ_KEYWORDS:
        if kw.lower() in names_str:
            return "bulk_rnaseq"
    return "unknown"


def detect_control_types(filenames: list[str], metadata: dict) -> list[str]:
    """
    Infer available control types (GC/VC/BC) from file names and metadata.
    Ground Control (GC) = flight hardware on Earth
    Vivarium Control (VC) = standard cage
    Basal Control (BC) = pre-launch
    """
    names_str = " ".join(filenames).lower()
    meta_str = json.dumps(metadata).lower()
    combined = names_str + " " + meta_str

    controls = []
    if any(kw in combined for kw in ["ground_control", "groundcontrol", "grd_ctrl", "gc_"]):
        controls.append("GC")
    if any(kw in combined for kw in ["vivarium", "viv_ctrl", "vc_", "cage_control"]):
        controls.append("VC")
    if any(kw in combined for kw in ["basal", "baseline", "pre_flight", "preflight", "bc_"]):
        controls.append("BC")

    # Fall back: if no keywords found, check if we see both flight and control conditions
    # (common in GeneLab standard processed files)
    return controls if controls else ["unknown"]


def search_hu_studies(max_results: int = 20) -> list[dict]:
    """
    Search OSDR for Hindlimb Unloading mouse RNA-seq studies.
    Uses OSDR search API with text queries.
    """
    hu_studies = []
    seen_ids = set()
    search_url = "https://osdr.nasa.gov/osdr/data/osd/meta/search"

    for term in HU_QUERY_TERMS:
        params = {
            "term": term,
            "from": 0,
            "size": max_results,
            "type": "cgene",
        }
        try:
            resp = requests.get(search_url, params=params, timeout=30)
            if resp.status_code != 200:
                continue
            data = resp.json()
            hits = data.get("hits", {}).get("hits", [])
            for hit in hits:
                src = hit.get("_source", {})
                osd_id = src.get("Accession", "")
                if osd_id and osd_id not in seen_ids:
                    seen_ids.add(osd_id)
                    hu_studies.append({
                        "osd_id": osd_id,
                        "title": src.get("Study Title", ""),
                        "organism": src.get("Organism", ""),
                        "assay": src.get("Assay", ""),
                        "search_term": term,
                    })
            time.sleep(0.5)
        except Exception:
            continue

    return hu_studies


# ── Verification logic ─────────────────────────────────────────────────────────

def verify_dataset(candidate: dict, verbose: bool = True) -> dict:
    """
    Verify a single candidate GLDS/OSD study.
    Returns enriched result dict.
    """
    osd_id = candidate["osd_id"]
    if verbose:
        print(f"  Checking {osd_id} ({candidate['tissue']}, {candidate['mission']})...", end=" ", flush=True)

    files = get_osd_files(osd_id)
    filenames = list(files.keys())

    if not filenames:
        status = "❌ NOT FOUND"
        assay_type = "unknown"
        control_types = []
        n_files = 0
        bulk_rnaseq_files = []
    else:
        assay_type = detect_assay_type(filenames)
        metadata = get_osd_metadata(osd_id)
        control_types = detect_control_types(filenames, metadata)
        n_files = len(filenames)

        # Identify bulk RNA-seq processed files
        bulk_rnaseq_files = [
            f for f in filenames
            if any(kw.lower() in f.lower() for kw in BULK_RNASEQ_KEYWORDS)
        ]

        if candidate["track"] == "exc":
            status = "⛔ EXCLUDED (design)"
        elif candidate["track"] == "qc":
            status = "⚠️ FAILS QC (n<3/group)"
        elif candidate["track"] == "HU":
            status = f"🔬 HU ANALOG ({assay_type})"
        elif assay_type == "bulk_rnaseq":
            status = "✅ Verified (bulk RNA-seq)"
        elif assay_type == "spatial":
            status = "⚠️ SPATIAL (not Category A)"
        elif assay_type == "single_cell":
            status = "⚠️ SINGLE-CELL (not Category A)"
        elif assay_type == "microarray":
            status = "⚠️ MICROARRAY (not Category A)"
        else:
            status = "⚠️ UNKNOWN assay type — manual check"

        if "DUPLICATE" in candidate.get("notes", ""):
            status = "🔁 DUPLICATE — Category J only"

    if verbose:
        print(status)

    result = {
        **candidate,
        "status": status,
        "assay_type": assay_type,
        "control_types_detected": control_types,
        "n_files": n_files,
        "bulk_rnaseq_files_found": len(bulk_rnaseq_files),
        "bulk_rnaseq_file_examples": bulk_rnaseq_files[:3],  # first 3 examples
        "verified_at": datetime.utcnow().isoformat(),
    }
    return result


# ── Output generation ──────────────────────────────────────────────────────────

def write_data_catalog_md(results: list[dict], hu_results: list[dict]) -> None:
    """Generate DATA_CATALOG.md from verification results."""
    lines = [
        "# DATA_CATALOG.md — GeneLab_benchmark",
        f"**자동 생성**: `catalog_datasets.py` — {datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC')}",
        "**PLAN.md 버전**: v0.5",
        "",
        "> ⚠️ Status legend:",
        "> ✅ Verified (bulk RNA-seq) | ❌ NOT FOUND | ⚠️ Wrong assay type",
        "> 🔁 DUPLICATE (excluded from Category A) | ⛔ EXCLUDED by design",
        "",
        "---",
        "",
        "## Candidate Datasets (PLAN.md v0.5)",
        "",
        "| OSD ID | Tissue | Mission | Track | Status | Assay | Control Types | Bulk RNA-seq Files | Notes |",
        "|---|---|---|---|---|---|---|---|---|",
    ]

    for r in results:
        ctrl = ", ".join(r.get("control_types_detected", [])) or "?"
        n_bulk = r.get("bulk_rnaseq_files_found", 0)
        lines.append(
            f"| {r['osd_id']} | {r['tissue']} | {r['mission']} | {r['track']} "
            f"| {r['status']} | {r['assay_type']} | {ctrl} | {n_bulk} | {r.get('notes', '')} |"
        )

    # Summary stats
    verified = [r for r in results if "✅" in r["status"]]
    not_found = [r for r in results if "❌" in r["status"]]
    wrong_assay = [r for r in results if "⚠️" in r["status"]]
    excluded = [r for r in results if "⛔" in r["status"] or "🔁" in r["status"]]

    lines += [
        "",
        "---",
        "",
        "## Summary",
        "",
        f"- ✅ Verified (bulk RNA-seq): **{len(verified)}**",
        f"- ❌ Not Found: **{len(not_found)}**",
        f"- ⚠️ Wrong assay type: **{len(wrong_assay)}**",
        f"- ⛔/🔁 Excluded by design: **{len(excluded)}**",
        "",
        "---",
        "",
    ]

    if hu_results:
        lines += [
            "## Hindlimb Unloading (HU) Analog Studies Found",
            "",
            "| OSD ID | Title | Organism | Assay | Search Term |",
            "|---|---|---|---|---|",
        ]
        for h in hu_results:
            title_short = h['title'][:60] + "..." if len(h['title']) > 60 else h['title']
            lines.append(
                f"| {h['osd_id']} | {title_short} | {h['organism']} | {h['assay']} | {h['search_term']} |"
            )
        lines += ["", "---", ""]

    lines += [
        "## Notes on Excluded Studies",
        "",
        "### GLDS-168 / OSD-168 — DUPLICATE",
        "OSD-168 is a reprocessed combination of OSD-48 (RR-1) + OSD-137 (RR-3) liver data.",
        "Including in Category A would cause sample duplication across LOMO train/test folds.",
        "**Use only in Category J (J1): pipeline version comparison (v1 vs v2 vs v3).**",
        "Reference: DESIGN_DECISIONS.md DD-05",
        "",
        "### OSD-25 (STS-135) — Track 2b only",
        "Uses BALB/c mouse strain, not C57BL/6J. Excluded from Track 2a (primary analysis).",
        "Track 2b (secondary, supplementary): included as cross-strain experiment.",
        "Reference: DESIGN_DECISIONS.md DD-06",
        "",
        "### OSD-270 (Heart, Visium) — NOT Category A",
        "Visium spatial transcriptomics, not bulk RNA-seq. Category F (v2.0) only.",
        "Reference: PLAN.md v0.5 Section 3.2",
    ]

    DATA_CATALOG_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\n✓ DATA_CATALOG.md written → {DATA_CATALOG_MD}")


def write_verified_json(results: list[dict], hu_results: list[dict]) -> None:
    """Write machine-readable GLDS_verified.json for use by fetch_osdr.py."""
    output = {
        "generated_at": datetime.utcnow().isoformat(),
        "plan_version": "0.5",
        "datasets": {
            r["osd_id"]: {
                "tissue": r["tissue"],
                "mission": r["mission"],
                "track": r["track"],
                "status": r["status"],
                "assay_type": r["assay_type"],
                "control_types_detected": r.get("control_types_detected", []),
                "bulk_rnaseq_files_found": r.get("bulk_rnaseq_files_found", 0),
                "bulk_rnaseq_file_examples": r.get("bulk_rnaseq_file_examples", []),
                "notes": r.get("notes", ""),
                "include_category_a": "✅" in r.get("status", ""),
                "include_category_j": "🔁" in r.get("status", "") or "✅" in r.get("status", ""),
            }
            for r in results
        },
        "hu_analog_studies": hu_results,
    }
    GLDS_VERIFIED_JSON.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(f"✓ GLDS_verified.json written → {GLDS_VERIFIED_JSON}")


# ── CLI ────────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="Verify GLDS/OSD IDs via OSDR API and generate DATA_CATALOG.md"
    )
    parser.add_argument(
        "--osd", type=str, default=None,
        help="Check a single OSD ID (e.g. OSD-48). Skip full catalog run."
    )
    parser.add_argument(
        "--search-hu", action="store_true",
        help="Also search for Hindlimb Unloading analog studies."
    )
    parser.add_argument(
        "--delay", type=float, default=0.5,
        help="Delay (seconds) between API calls to avoid rate limiting. Default: 0.5"
    )
    parser.add_argument(
        "--no-output", action="store_true",
        help="Skip writing output files (for testing)."
    )
    return parser.parse_args()


def main():
    args = parse_args()

    if args.osd:
        # Single-study quick check mode
        osd_id = args.osd if args.osd.startswith("OSD-") else f"OSD-{args.osd}"
        candidate = {"osd_id": osd_id, "tissue": "?", "mission": "?", "track": "?", "notes": "manual check"}
        result = verify_dataset(candidate, verbose=True)
        print(json.dumps(result, indent=2))
        return

    print("=" * 60)
    print("GeneLab_benchmark — OSDR Catalog Verification")
    print(f"Checking {len(CANDIDATE_DATASETS)} candidate studies...")
    print("=" * 60)

    results = []
    for candidate in CANDIDATE_DATASETS:
        result = verify_dataset(candidate, verbose=True)
        results.append(result)
        time.sleep(args.delay)

    hu_results = []
    if args.search_hu:
        print("\nSearching for Hindlimb Unloading (HU) analog studies...")
        hu_results = search_hu_studies()
        print(f"  Found {len(hu_results)} HU candidate studies.")

    if not args.no_output:
        write_data_catalog_md(results, hu_results)
        write_verified_json(results, hu_results)

    # Print summary
    print("\n" + "=" * 60)
    verified = [r for r in results if "✅" in r["status"]]
    not_found = [r for r in results if "❌" in r["status"]]
    wrong = [r for r in results if "⚠️" in r["status"]]
    print(f"✅ Verified (bulk RNA-seq): {len(verified)}")
    print(f"❌ Not Found:              {len(not_found)}")
    print(f"⚠️  Wrong assay/unknown:   {len(wrong)}")
    print("=" * 60)

    if not_found:
        print("\nNot found (check OSD IDs):")
        for r in not_found:
            print(f"  {r['osd_id']} — {r['tissue']} {r['mission']}")

    if wrong:
        print("\nWrong assay type (manual check needed):")
        for r in wrong:
            print(f"  {r['osd_id']} — {r['tissue']} {r['mission']}: {r['assay_type']}")


if __name__ == "__main__":
    main()
