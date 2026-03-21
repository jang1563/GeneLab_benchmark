#!/usr/bin/env python3
"""
preprocess_radiation.py — Phase 5 radiation × HLU 2×2 데이터 전처리.

Datasets:
  OSD-202 (GLDS-202): Brain, 2×2×2 (radiation × HLU × time), n=40
  OSD-211 (GLDS-211): Spleen, 2×2 (radiation × HLU), n=21
  OSD-237 (GLDS-237): Skin, 2×2 (radiation × HLU), n=21

Output:
  processed/R_radiation/{tissue}/{tissue}_log2_norm.csv
  processed/R_radiation/{tissue}/{tissue}_metadata.csv

Metadata columns:
  sample_name, radiation (0/1), hlu (0/1), timepoint, condition_raw, tissue, osd_id
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent.parent
DATA_DIR = BASE_DIR / "data" / "mouse"
PROCESSED_DIR = BASE_DIR / "processed" / "R_radiation"

DATASETS = [
    {
        "osd_id": 202, "glds_id": 202, "tissue": "brain",
        "data_dir": "retina_brain/LDR_HLU",
        "counts_file": "GLDS-202_rna_seq_Normalized_Counts_GLbulkRNAseq.csv",
        "sample_table": "GLDS-202_rna_seq_SampleTable_GLbulkRNAseq.csv",
        "runsheet": "GLDS-202_rna_seq_bulkRNASeq_v2_runsheet.csv",
        "has_runsheet": True,
    },
    {
        "osd_id": 211, "glds_id": 211, "tissue": "spleen",
        "data_dir": "spleen/LDR_HLU",
        "counts_file": "GLDS-211_rna_seq_Normalized_Counts.csv",
        "sample_table": "GLDS-211_rna_seq_SampleTable.csv",
        "has_runsheet": False,
    },
    {
        "osd_id": 237, "glds_id": 237, "tissue": "skin_rad",
        "data_dir": "skin_rad/LDR_HLU",
        "counts_file": "GLDS-237_rna_seq_Normalized_Counts.csv",
        "sample_table": "GLDS-237_rna_seq_SampleTable.csv",
        "has_runsheet": False,
    },
]


def parse_condition_202(row):
    """Parse OSD-202 condition using runsheet factor values."""
    rad_val = row.get("Factor Value[Ionizing Radiation]", "")
    hlu_val = row.get("Factor Value[Hindlimb Unloading]", "")
    time_val = row.get("Factor Value[Time of Sample Collection After Treatment]", "")

    rad_str = str(rad_val).lower()
    radiation = 1 if "cobalt" in rad_str or "gamma" in rad_str else 0
    hlu = 1 if "unloaded" in str(hlu_val).lower() else 0
    timepoint = str(time_val).strip()

    return {"radiation": radiation, "hlu": hlu, "timepoint": timepoint}


def parse_condition_211(condition_str):
    """Parse OSD-211 abbreviated condition: HLLC_IRC, HLLC_IR, HLU_IRC, HLU_IR."""
    cond = str(condition_str).strip()
    hlu = 1 if cond.startswith("HLU") else 0
    radiation = 1 if cond.endswith("_IR") and not cond.endswith("_IRC") else 0
    return {"radiation": radiation, "hlu": hlu, "timepoint": ""}


def parse_condition_237(condition_str):
    """Parse OSD-237 long condition strings."""
    cond = str(condition_str).lower().replace("...", ".").replace("..", ".")
    hlu = 1 if "unloaded" in cond else 0
    radiation = 1 if "irradiated.with" in cond or "0.04.gy" in cond else 0
    return {"radiation": radiation, "hlu": hlu, "timepoint": ""}


def process_dataset(ds):
    """Process a single radiation dataset."""
    tissue = ds["tissue"]
    osd_id = ds["osd_id"]
    data_path = DATA_DIR / ds["data_dir"]

    counts_path = data_path / ds["counts_file"]
    st_path = data_path / ds["sample_table"]

    print(f"\n--- OSD-{osd_id}: {tissue} ---")

    if not counts_path.exists():
        print(f"  SKIP: counts not found: {counts_path}")
        return None

    # Load counts
    counts = pd.read_csv(counts_path, index_col=0)
    if not str(counts.index[0]).startswith("ENSMUSG"):
        if str(counts.columns[0]).startswith("ENSMUSG"):
            counts = counts.T
    ensmusg_mask = [str(g).startswith("ENSMUSG") for g in counts.index]
    counts = counts.loc[ensmusg_mask]
    counts_numeric = counts.apply(pd.to_numeric, errors="coerce").fillna(0)
    log2_counts = np.log2(counts_numeric + 1)
    log2_T = log2_counts.T  # samples × genes
    print(f"  log2_norm: {log2_T.shape} (samples × genes)")

    # Load sample table + parse conditions
    st = pd.read_csv(st_path)
    if "Unnamed: 0" in st.columns:
        st = st.rename(columns={"Unnamed: 0": "sample_name"})
    elif st.columns[0] != "sample_name":
        st = st.rename(columns={st.columns[0]: "sample_name"})

    # For OSD-202, merge with runsheet for clean factor values
    if ds.get("has_runsheet"):
        rs_path = data_path / ds["runsheet"]
        rs = pd.read_csv(rs_path)
        # Merge on sample name
        rs = rs.rename(columns={"Sample Name": "sample_name"})
        st = st.merge(rs[["sample_name",
                           "Factor Value[Ionizing Radiation]",
                           "Factor Value[Hindlimb Unloading]",
                           "Factor Value[Time of Sample Collection After Treatment]"]],
                       on="sample_name", how="left")

    # Parse conditions into metadata
    meta_rows = []
    for _, row in st.iterrows():
        sample = row["sample_name"]
        condition = str(row.get("condition", ""))

        if ds.get("has_runsheet"):
            parsed = parse_condition_202(row)
        elif osd_id == 211:
            parsed = parse_condition_211(condition)
        elif osd_id == 237:
            parsed = parse_condition_237(condition)
        else:
            parsed = {"radiation": 0, "hlu": 0, "timepoint": ""}

        parsed["sample_name"] = sample
        parsed["condition_raw"] = condition
        parsed["tissue"] = tissue
        parsed["osd_id"] = f"OSD-{osd_id}"

        # Derived group label
        rad_str = "Rad" if parsed["radiation"] else "Ctrl"
        hlu_str = "HLU" if parsed["hlu"] else "NL"
        parsed["group"] = f"{hlu_str}_{rad_str}"

        meta_rows.append(parsed)

    meta = pd.DataFrame(meta_rows).set_index("sample_name")

    # Align
    common = sorted(set(log2_T.index) & set(meta.index))
    log2_T = log2_T.loc[common]
    meta = meta.loc[common]

    # Summary
    group_counts = meta["group"].value_counts().to_dict()
    print(f"  Samples: {len(meta)} | Groups: {group_counts}")

    # Save
    out_dir = PROCESSED_DIR / tissue
    out_dir.mkdir(parents=True, exist_ok=True)
    log2_T.to_csv(out_dir / f"{tissue}_log2_norm.csv")
    meta.to_csv(out_dir / f"{tissue}_metadata.csv")
    print(f"  Saved to: {out_dir}")

    return {
        "tissue": tissue,
        "osd_id": osd_id,
        "n_samples": len(meta),
        "n_genes": log2_T.shape[1],
        "groups": group_counts,
    }


def main():
    print("=" * 70)
    print("GeneLabBench v3 — Phase 5 Radiation Preprocessing")
    print("=" * 70)

    results = []
    for ds in DATASETS:
        r = process_dataset(ds)
        if r:
            results.append(r)

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    for r in results:
        print(f"  {r['tissue']} (OSD-{r['osd_id']}): {r['n_samples']} samples × {r['n_genes']} genes")
        print(f"    Groups: {r['groups']}")

    summary_path = PROCESSED_DIR / "phase5_preprocess_summary.json"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with open(summary_path, "w") as f:
        json.dump(results, f, indent=2, default=str)


if __name__ == "__main__":
    main()
