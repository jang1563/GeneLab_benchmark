#!/usr/bin/env python3
"""
catalog_v3_datasets.py — Automated verification of v3 candidate dataset metadata using the OSDR API.

Usage:
    python3 v3/scripts/catalog_v3_datasets.py

Output:
    v3/docs/DATA_CATALOG_V3.md — Verification results in Markdown
    v3/processed/dataset_catalog_v3.json — Raw JSON output
"""

import json
import sys
import time
from pathlib import Path
from urllib.request import urlopen, Request
from urllib.error import URLError, HTTPError

BASE_DIR = Path(__file__).resolve().parent.parent.parent  # GeneLab_benchmark/
OUTPUT_MD = BASE_DIR / "v3" / "docs" / "DATA_CATALOG_V3.md"
OUTPUT_JSON = BASE_DIR / "v3" / "processed" / "dataset_catalog_v3.json"

# OSDR API
OSDR_META_API = "https://osdr.nasa.gov/osdr/data/osd/meta/{osd_id}"
OSDR_FILES_API = "https://osdr.nasa.gov/osdr/data/osd/files/{osd_id}"

# Datasets to verify
DATASETS = [
    # Phase 1: Multi-Species
    {"osd_id": 207, "expected": "Drosophila head RNA-seq, spaceflight", "phase": 1, "group": "Multi-Species"},
    {"osd_id": 7,   "expected": "Arabidopsis RNA-seq, spaceflight (ABRS)", "phase": 1, "group": "Multi-Species"},
    {"osd_id": 37,  "expected": "Arabidopsis RNA-seq, spaceflight (BRIC)", "phase": 1, "group": "Multi-Species"},
    {"osd_id": 38,  "expected": "Arabidopsis RNA-seq + proteomics", "phase": 1, "group": "Multi-Species"},
    {"osd_id": 120, "expected": "Arabidopsis RNA-seq, spaceflight (ISS)", "phase": 1, "group": "Multi-Species"},
    {"osd_id": 208, "expected": "Arabidopsis root RNA-seq (NOT spaceflight)", "phase": 1, "group": "Multi-Species"},
    {"osd_id": 113, "expected": "C. elegans spaceflight (microarray?)", "phase": 1, "group": "C. elegans"},
    {"osd_id": 112, "expected": "C. elegans spaceflight", "phase": 1, "group": "C. elegans"},
    {"osd_id": 41,  "expected": "C. elegans spaceflight", "phase": 1, "group": "C. elegans"},
    # Phase 2: Spatial
    {"osd_id": 270, "expected": "RR-3 heart Visium spatial", "phase": 2, "group": "Spatial"},
    # Phase 3: RRRM-2 scRNA-seq
    {"osd_id": 402, "expected": "RRRM-2 femur bone marrow scRNA-seq", "phase": 3, "group": "RRRM-2 scRNA"},
    {"osd_id": 403, "expected": "RRRM-2 humerus bone marrow scRNA-seq", "phase": 3, "group": "RRRM-2 scRNA"},
    {"osd_id": 404, "expected": "RRRM-2 PBMCs scRNA-seq", "phase": 3, "group": "RRRM-2 scRNA"},
    {"osd_id": 405, "expected": "RRRM-2 spleen scRNA-seq", "phase": 3, "group": "RRRM-2 scRNA"},
    # Phase 4: Existing Extension
    {"osd_id": 253, "expected": "Kidney RR-7 RNA-seq (extends kidney)", "phase": 4, "group": "New Bulk"},
    {"osd_id": 240, "expected": "Dorsal skin RR-5 RNA-seq", "phase": 4, "group": "New Bulk"},
    {"osd_id": 241, "expected": "Femoral skin RR-5 RNA-seq", "phase": 4, "group": "New Bulk"},
    {"osd_id": 243, "expected": "Dorsal skin RR-6 RNA-seq", "phase": 4, "group": "New Bulk"},
    {"osd_id": 248, "expected": "Lung RR-6 RNA-seq", "phase": 4, "group": "New Bulk"},
    {"osd_id": 247, "expected": "Colon RR-6 RNA-seq (new tissue)", "phase": 4, "group": "New Bulk"},
    # Phase 5: Radiation
    {"osd_id": 202, "expected": "Retina+brain, LDR+HLU, RNA-seq", "phase": 5, "group": "Radiation"},
    {"osd_id": 211, "expected": "Spleen, LDR+HLU, RNA-seq", "phase": 5, "group": "Radiation"},
    {"osd_id": 237, "expected": "Skin, LDR+HLU, RNA-seq", "phase": 5, "group": "Radiation"},
    {"osd_id": 73,  "expected": "Radiation, microarray (verify)", "phase": 5, "group": "Radiation (microarray)"},
    {"osd_id": 109, "expected": "Radiation, microarray (verify)", "phase": 5, "group": "Radiation (microarray)"},
    # Phase 6: Multimodal
    {"osd_id": 804, "expected": "RR-1 MicroCT bone imaging", "phase": 6, "group": "Multimodal"},
]


def fetch_json(url, timeout=30):
    """Fetch JSON from URL with error handling."""
    try:
        req = Request(url, headers={"User-Agent": "GeneLabBench/3.0"})
        with urlopen(req, timeout=timeout) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except (URLError, HTTPError, json.JSONDecodeError) as e:
        return {"_error": str(e)}


def extract_metadata(osd_id):
    """Extract study metadata from OSDR meta API."""
    info = {
        "osd_id": f"OSD-{osd_id}",
        "title": "",
        "description": "",
        "organism": "",
        "tissue": "",
        "assay_types": [],
        "sample_count": 0,
        "factors": [],
        "factor_values": {},
        "first_data_file": "",
        "release_date": "",
        "status": "error",
    }

    # Fetch metadata
    raw = fetch_json(OSDR_META_API.format(osd_id=osd_id))
    if "_error" in raw:
        info["status"] = f"API error: {raw['_error']}"
        return info

    key = f"OSD-{osd_id}"
    if "study" not in raw or key not in raw["study"]:
        info["status"] = "No study data in response"
        return info

    osd_data = raw["study"][key]
    studies = osd_data.get("studies", [])
    if not studies:
        info["status"] = "No studies array"
        return info

    s = studies[0]

    # Title + description
    info["title"] = s.get("title", "")
    info["description"] = s.get("description", "")[:300]
    info["release_date"] = s.get("publicReleaseDate", "")

    # Build @id → name lookups
    char_cat_lookup = {}  # @id → characteristicType annotationValue
    for cc in s.get("characteristicCategories", []):
        cid = cc.get("@id", "")
        ctype = cc.get("characteristicType", {}).get("annotationValue", "")
        if cid and ctype:
            char_cat_lookup[cid] = ctype

    factor_lookup = {}  # @id → factorName
    for f in s.get("factors", []):
        fid = f.get("@id", "")
        fname = f.get("factorName", "")
        if fid and fname:
            factor_lookup[fid] = fname

    # Organism + tissue from samples (sources often have empty characteristics)
    samples = s.get("materials", {}).get("samples", [])
    info["sample_count"] = len(samples)

    if samples:
        sample0 = samples[0]
        for ch in sample0.get("characteristics", []):
            cat_id = ch.get("category", {}).get("@id", "")
            cat_name = char_cat_lookup.get(cat_id, "")
            val = ch.get("value", {})
            val_name = val.get("annotationValue", "") if isinstance(val, dict) else str(val)
            if cat_name.lower() == "organism":
                info["organism"] = val_name
            elif "tissue" in cat_name.lower() or cat_name.lower() in ("organism part", "organ"):
                info["tissue"] = val_name

        # Factor values from ALL samples to get unique values
        for sample in samples:
            for fv in sample.get("factorValues", []):
                fcat_id = fv.get("category", {}).get("@id", "")
                fname = factor_lookup.get(fcat_id, "")
                fval = fv.get("value", {})
                fval_name = fval.get("annotationValue", "") if isinstance(fval, dict) else str(fval)
                if fname:
                    info["factor_values"].setdefault(fname, set()).add(fval_name)

    # Convert sets to sorted lists for JSON serialization
    for k, v in info["factor_values"].items():
        info["factor_values"][k] = sorted(v)

    # Factors
    for f in s.get("factors", []):
        fname = f.get("factorName", "")
        ftype = f.get("factorType", {}).get("annotationValue", "")
        info["factors"].append(f"{fname} ({ftype})" if ftype else fname)

    # Assays
    for a in s.get("assays", []):
        filename = a.get("filename", "")
        # Detect assay type from filename
        fl = filename.lower()
        if "rna-seq" in fl or "rna-sequencing" in fl:
            info["assay_types"].append("RNA-seq")
        elif "spatial-transcriptom" in fl:
            info["assay_types"].append("Spatial (Visium)")
        elif "microarray" in fl or "genechip" in fl or "affymetrix" in fl:
            info["assay_types"].append("Microarray")
        elif "mass-spectrometry" in fl or "proteom" in fl:
            info["assay_types"].append("Proteomics")
        elif "methylation" in fl or "methyl" in fl:
            info["assay_types"].append("Methylation")
        elif "metabol" in fl:
            info["assay_types"].append("Metabolomics")
        elif "atac" in fl:
            info["assay_types"].append("ATAC-seq")
        elif "metagenomic" in fl or "amplicon" in fl:
            info["assay_types"].append("Metagenomics")
        elif "micro-computed-tomography" in fl or "microct" in fl:
            info["assay_types"].append("MicroCT")
        elif "image" in fl or "photography" in fl:
            info["assay_types"].append("Imaging")
        else:
            info["assay_types"].append(f"Unknown ({filename[:60]})")

        # Get first data file
        data_files = a.get("dataFiles", [])
        if data_files and not info["first_data_file"]:
            info["first_data_file"] = data_files[0].get("name", "")

    info["status"] = "ok"
    return info


def determine_go_nogo(info, expected):
    """Determine if dataset meets Go criteria."""
    assays = [a.lower() for a in info["assay_types"]]
    has_rnaseq = any("rna-seq" in a for a in assays)
    has_microarray = any("microarray" in a for a in assays)
    has_samples = info["sample_count"] >= 4  # minimum

    if info["status"] != "ok":
        return "❌ API error"

    if "microarray" in expected.lower():
        # We expected microarray
        if has_microarray:
            return "✅ Microarray (as expected)"
        elif has_rnaseq:
            return "✅ RNA-seq (better than expected!)"
        else:
            return "⚠️ Unknown assay"

    # We expected RNA-seq or similar
    if has_rnaseq:
        if has_samples:
            return "✅ RNA-seq"
        else:
            return "⚠️ RNA-seq but low n"
    elif has_microarray:
        return "⚠️ Microarray (not RNA-seq)"
    elif "scrna" in expected.lower() or "scRNA" in expected:
        # Check if any assay mentions single-cell
        return "⚠️ Check scRNA-seq files"
    elif "visium" in expected.lower() or "spatial" in expected.lower():
        return "⚠️ Check spatial files"
    elif "microct" in expected.lower():
        return "⚠️ Non-omics (MicroCT)"
    else:
        return "⚠️ Unknown assay type"


def main():
    print("=" * 80)
    print("GeneLabBench v3 — OSDR Dataset Catalog Verification")
    print("=" * 80)

    results = []

    for ds in DATASETS:
        osd_id = ds["osd_id"]
        expected = ds["expected"]
        phase = ds["phase"]
        group = ds["group"]

        print(f"\n  OSD-{osd_id} (expected: {expected})")
        info = extract_metadata(osd_id)
        info["expected"] = expected
        info["phase"] = phase
        info["group"] = group

        # Determine Go/No-Go
        info["go_nogo"] = determine_go_nogo(info, expected)

        print(f"    Title: {info['title'][:80]}")
        print(f"    Organism: {info['organism']}")
        print(f"    Assays: {info['assay_types']}")
        print(f"    Samples: {info['sample_count']}")
        print(f"    Factors: {info['factors']}")
        print(f"    Factor values: {dict(list(info['factor_values'].items())[:3])}")
        print(f"    Go/No-Go: {info['go_nogo']}")

        results.append(info)
        time.sleep(0.3)  # rate limit

    # Save JSON
    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_JSON, "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False, default=str)
    print(f"\n\nJSON saved to: {OUTPUT_JSON}")

    # Generate markdown
    generate_markdown(results)
    print(f"Markdown saved to: {OUTPUT_MD}")

    # Print summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    go_count = sum(1 for r in results if r["go_nogo"].startswith("✅"))
    warn_count = sum(1 for r in results if r["go_nogo"].startswith("⚠️"))
    fail_count = sum(1 for r in results if r["go_nogo"].startswith("❌"))
    print(f"  ✅ Go: {go_count}  |  ⚠️ Warning: {warn_count}  |  ❌ Fail: {fail_count}")
    print(f"  Total datasets checked: {len(results)}")


def generate_markdown(results):
    """Generate DATA_CATALOG_V3.md."""
    lines = [
        "# GeneLabBench v3 — Data Catalog",
        "",
        "**Generated**: 2026-03-18",
        "**Method**: OSDR API automated verification (`v3/scripts/catalog_v3_datasets.py`)",
        "",
        "## Summary",
        "",
        "| OSD ID | Organism | Assay | n | Title | Phase | Go? |",
        "|--------|----------|-------|---|-------|-------|-----|",
    ]

    for r in sorted(results, key=lambda x: (x["phase"], x["osd_id"])):
        assay = ", ".join(r["assay_types"]) if r["assay_types"] else "?"
        title_short = r["title"][:45] + ("..." if len(r["title"]) > 45 else "")
        org = r["organism"][:20] if r["organism"] else "?"
        lines.append(
            f"| {r['osd_id']} | {org} | {assay} | {r['sample_count']} | {title_short} | {r['phase']} | {r['go_nogo'][:3]} |"
        )

    # Group by phase
    by_phase = {}
    for r in results:
        by_phase.setdefault(r["phase"], []).append(r)

    for phase in sorted(by_phase.keys()):
        lines.extend(["", f"## Phase {phase} Details", ""])
        for r in by_phase[phase]:
            lines.append(f"### {r['osd_id']} — {r['go_nogo']}")
            lines.append(f"- **Title**: {r['title']}")
            lines.append(f"- **Organism**: {r['organism']}")
            lines.append(f"- **Assay types**: {', '.join(r['assay_types']) if r['assay_types'] else 'Unknown'}")
            lines.append(f"- **Sample count**: {r['sample_count']}")
            lines.append(f"- **Factors**: {', '.join(r['factors']) if r['factors'] else 'N/A'}")
            if r["factor_values"]:
                for fname, fvals in r["factor_values"].items():
                    lines.append(f"  - {fname}: {', '.join(fvals)}")
            lines.append(f"- **Expected**: {r['expected']}")
            lines.append(f"- **Description**: {r['description']}")
            if r["first_data_file"]:
                lines.append(f"- **First data file**: `{r['first_data_file']}`")
            lines.append("")

    OUTPUT_MD.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_MD, "w") as f:
        f.write("\n".join(lines))


if __name__ == "__main__":
    main()
