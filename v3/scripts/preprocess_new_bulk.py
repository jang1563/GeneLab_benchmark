#!/usr/bin/env python3
"""
preprocess_new_bulk.py — Phase 4 신규 bulk RNA-seq 데이터 전처리.

v1 quality_filter.py 패턴을 따름:
  - Normalized counts → log2(x+1) → samples × genes CSV
  - SampleTable → metadata CSV with standardized columns

Handles:
  - Complex condition strings (RR-6 format: "Space.Flight....30.day...On.Earth...")
  - Simple condition strings (RR-5 format: "Ground.Control", "Space.Flight")
  - Multiple strains (C57BL/6J Track 2a, BAL-TAL Track 2b)

Usage:
    python3 v3/scripts/preprocess_new_bulk.py                    # All Phase 4
    python3 v3/scripts/preprocess_new_bulk.py --tissue lung      # Single tissue
    python3 v3/scripts/preprocess_new_bulk.py --dry-run          # Preview

Output:
    processed/A_detection/{tissue}/{tissue}_{mission}_log2_norm.csv
    processed/A_detection/{tissue}/{tissue}_{mission}_metadata.csv
"""

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent.parent  # GeneLab_benchmark/
DATA_DIR = BASE_DIR / "data" / "mouse"
PROCESSED_DIR = BASE_DIR / "processed" / "A_detection"

# Dataset definitions for Phase 4
PHASE4_DATASETS = [
    # New tissues
    {
        "osd_id": 248, "tissue": "lung", "mission": "RR-6",
        "data_dir": "lung/RR-6",
        "counts_file": "GLDS-248_rna_seq_Normalized_Counts_GLbulkRNAseq.csv",
        "sample_table": "GLDS-248_rna_seq_SampleTable_GLbulkRNAseq.csv",
    },
    {
        "osd_id": 247, "tissue": "colon", "mission": "RR-6",
        "data_dir": "colon/RR-6",
        "counts_file": "GLDS-247_rna_seq_Normalized_Counts_GLbulkRNAseq.csv",
        "sample_table": "GLDS-247_rna_seq_SampleTable_GLbulkRNAseq.csv",
    },
    # Extended skin (RR-5 dorsal — BAL-TAL strain → Track 2b)
    {
        "osd_id": 240, "tissue": "skin", "mission": "RR-5_dorsal",
        "data_dir": "skin/RR-5_dorsal",
        "counts_file": "GLDS-240_rna_seq_Normalized_Counts.csv",
        "sample_table": "GLDS-240_rna_seq_SampleTable.csv",
        "strain": "BAL-TAL",
    },
    # Extended skin (RR-5 femoral — BAL-TAL strain → Track 2b)
    {
        "osd_id": 241, "tissue": "skin", "mission": "RR-5_femoral",
        "data_dir": "skin/RR-5_femoral",
        "counts_file": "GLDS-241_rna_seq_Normalized_Counts.csv",
        "sample_table": "GLDS-241_rna_seq_SampleTable.csv",
        "strain": "BAL-TAL",
    },
    # Extended skin (RR-6 dorsal — C57BL/6J)
    {
        "osd_id": 243, "tissue": "skin", "mission": "RR-6_dorsal",
        "data_dir": "skin/RR-6_dorsal",
        "counts_file": "GLDS-243_rna_seq_Normalized_Counts_GLbulkRNAseq.csv",
        "sample_table": "GLDS-243_rna_seq_SampleTable_GLbulkRNAseq.csv",
    },
]


def parse_condition(condition_str):
    """Parse OSDR condition string into standardized metadata fields.

    Handles both formats:
      - Simple: "Ground.Control", "Space.Flight"
      - Complex: "Space.Flight....30.day...On.Earth...Upon.euthanasia"
    """
    meta = {
        "label": "",
        "duration_days": "",
        "euthanasia_location": "",
        "dissection_condition": "",
    }

    # Normalize separator
    cond = condition_str.replace("...", ".").replace("..", ".")

    # Extract spaceflight label
    cond_lower = cond.lower()
    if "space.flight" in cond_lower or "spaceflight" in cond_lower:
        meta["label"] = "Flight"
    elif "ground.control" in cond_lower:
        meta["label"] = "Ground"
    elif "basal.control" in cond_lower:
        meta["label"] = "Basal"
    elif "vivarium" in cond_lower:
        meta["label"] = "Vivarium"
    else:
        meta["label"] = condition_str

    # Extract duration
    dur_match = re.search(r"(\d+)[\._\s]*day", cond, re.IGNORECASE)
    if dur_match:
        meta["duration_days"] = dur_match.group(1)

    # Extract euthanasia location
    if "on.iss" in cond_lower:
        meta["euthanasia_location"] = "ISS"
    elif "on.earth" in cond_lower:
        meta["euthanasia_location"] = "Earth"

    # Extract dissection condition
    if "carcass" in cond_lower:
        meta["dissection_condition"] = "Carcass"
    elif "upon.euthanasia" in cond_lower or "lar" in cond_lower:
        meta["dissection_condition"] = "LAR"

    return meta


def process_dataset(ds, dry_run=False):
    """Process a single dataset: counts + sample table → log2_norm + metadata."""
    tissue = ds["tissue"]
    mission = ds["mission"]
    osd_id = ds["osd_id"]

    counts_path = DATA_DIR / ds["data_dir"] / ds["counts_file"]
    st_path = DATA_DIR / ds["data_dir"] / ds["sample_table"]

    out_dir = PROCESSED_DIR / tissue
    out_log2 = out_dir / f"{tissue}_{mission}_log2_norm.csv"
    out_meta = out_dir / f"{tissue}_{mission}_metadata.csv"

    print(f"\n--- OSD-{osd_id}: {tissue}/{mission} ---")

    if not counts_path.exists():
        print(f"  SKIP: counts file not found: {counts_path}")
        return None
    if not st_path.exists():
        print(f"  SKIP: sample table not found: {st_path}")
        return None

    # Load counts (genes × samples)
    counts = pd.read_csv(counts_path, index_col=0)
    print(f"  Counts: {counts.shape} (raw)")

    # Ensure genes are rows (ENSMUSG in index)
    if not str(counts.index[0]).startswith("ENSMUSG"):
        if str(counts.columns[0]).startswith("ENSMUSG"):
            counts = counts.T
            print(f"  Transposed → {counts.shape}")
        else:
            print(f"  WARNING: No ENSMUSG IDs found in rows or columns!")

    # Filter to protein-coding genes (ENSMUSG only)
    ensmusg_mask = [str(g).startswith("ENSMUSG") for g in counts.index]
    counts = counts.loc[ensmusg_mask]
    print(f"  After ENSMUSG filter: {counts.shape}")

    # log2(x+1) transform
    counts_numeric = counts.apply(pd.to_numeric, errors="coerce").fillna(0)
    log2_counts = np.log2(counts_numeric + 1)

    # Transpose to samples × genes (v1 convention for log2_norm)
    log2_T = log2_counts.T
    print(f"  log2_norm shape: {log2_T.shape} (samples × genes)")

    # Load sample table
    st = pd.read_csv(st_path)
    # Handle unnamed first column (sample names)
    if "Unnamed: 0" in st.columns:
        st = st.rename(columns={"Unnamed: 0": "sample_name"})
    elif st.columns[0] != "sample_name":
        st = st.rename(columns={st.columns[0]: "sample_name"})

    # Parse conditions into metadata
    meta_rows = []
    for _, row in st.iterrows():
        sample = row["sample_name"]
        condition = str(row.get("condition", ""))
        parsed = parse_condition(condition)
        parsed["sample_name"] = sample
        parsed["mission"] = mission
        parsed["osd_id"] = f"OSD-{osd_id}"
        parsed["tissue"] = tissue
        parsed["condition_raw"] = condition
        parsed["strain"] = ds.get("strain", "C57BL/6J")
        meta_rows.append(parsed)

    meta = pd.DataFrame(meta_rows)
    meta = meta.set_index("sample_name")

    # Align counts and metadata
    common = sorted(set(log2_T.index) & set(meta.index))
    if len(common) < len(log2_T):
        print(f"  WARNING: {len(log2_T) - len(common)} samples not in sample table")
    log2_T = log2_T.loc[common]
    meta = meta.loc[common]

    # Summary
    label_counts = meta["label"].value_counts().to_dict()
    print(f"  Samples: {len(meta)} | Labels: {label_counts}")

    if dry_run:
        return {
            "tissue": tissue,
            "mission": mission,
            "n_samples": len(meta),
            "n_genes": log2_T.shape[1],
            "labels": label_counts,
        }

    # Save
    out_dir.mkdir(parents=True, exist_ok=True)
    log2_T.to_csv(out_log2)
    meta.to_csv(out_meta)
    print(f"  Saved: {out_log2.name}, {out_meta.name}")

    return {
        "tissue": tissue,
        "mission": mission,
        "n_samples": len(meta),
        "n_genes": log2_T.shape[1],
        "labels": label_counts,
        "output_log2": str(out_log2),
        "output_meta": str(out_meta),
    }


def main():
    parser = argparse.ArgumentParser(description="Preprocess Phase 4 bulk RNA-seq")
    parser.add_argument("--tissue", help="Process only this tissue")
    parser.add_argument("--dry-run", action="store_true", help="Preview only")
    args = parser.parse_args()

    print("=" * 70)
    print("GeneLabBench v3 — Phase 4 Bulk RNA-seq Preprocessing")
    print("=" * 70)

    datasets = PHASE4_DATASETS
    if args.tissue:
        datasets = [d for d in datasets if d["tissue"] == args.tissue]

    results = []
    for ds in datasets:
        r = process_dataset(ds, dry_run=args.dry_run)
        if r:
            results.append(r)

    # Save summary
    print("\n" + "=" * 70)
    print("PREPROCESSING SUMMARY")
    print("=" * 70)
    for r in results:
        strain_note = ""
        if any(d.get("strain") == "BAL-TAL" for d in PHASE4_DATASETS
               if d["tissue"] == r["tissue"] and d["mission"] == r["mission"]):
            strain_note = " [BAL-TAL]"
        print(f"  {r['tissue']}/{r['mission']}{strain_note}: "
              f"{r['n_samples']} samples × {r['n_genes']} genes | {r['labels']}")

    if not args.dry_run:
        summary_path = PROCESSED_DIR / "phase4_preprocess_summary.json"
        with open(summary_path, "w") as f:
            json.dump(results, f, indent=2, default=str)
        print(f"\nSummary saved to: {summary_path}")


if __name__ == "__main__":
    main()
