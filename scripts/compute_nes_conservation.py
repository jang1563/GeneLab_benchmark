#!/usr/bin/env python3
"""
compute_nes_conservation.py — Compute cross-mission NES Spearman r and update JSON.

For each tissue, computes pairwise Spearman correlation of NES scores (Hallmark)
across fGSEA mission outputs. Updates evaluation/NES_conservation_vs_transfer.json.

Usage:
  python scripts/compute_nes_conservation.py --tissue skin
  python scripts/compute_nes_conservation.py --all
"""
import json
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr
from itertools import combinations
from datetime import datetime

BASE_DIR  = Path(__file__).resolve().parent.parent
FGSEA_DIR = BASE_DIR / "processed" / "fgsea"
NES_JSON  = BASE_DIR / "evaluation" / "NES_conservation_vs_transfer.json"

# Transfer AUROC values sourced from B_cross_mission_summary.json (pca_lr)
TRANSFER_AUROC = {
    "liver":         {"auroc": 0.577,  "ci": [0.492, 0.666],  "n": 6},
    "gastrocnemius": {"auroc": 0.801,  "ci": [0.653, 0.944],  "n": 3},
    "kidney":        {"auroc": 0.555,  "ci": [0.397, 0.681],  "n": 3},
    "thymus":        {"auroc": 0.860,  "ci": [0.763, 0.953],  "n": 4},
    "eye":           {"auroc": 0.754,  "ci": [0.688, 0.838],  "n": 3},
    "skin":          {"auroc": 0.7718, "ci": [0.6914, 0.8338], "n": 3},
}

NOTES = {
    "gastrocnemius": "Only 2 missions have fGSEA data (RR-5 has no DGE file).",
    "skin": (
        "RR-7 (GLDS-254) excluded: only normalized counts available, no DGE. "
        "fGSEA run on 3 sources: RR-6, MHU-2_dorsal (GLDS-238), "
        "MHU-2_femoral (GLDS-239). Transfer AUROC (B5) uses merged MHU-2 + RR-7."
    ),
}


def compute_tissue(tissue: str) -> dict | None:
    """Compute NES pairwise Spearman r for one tissue."""
    tissue_dir = FGSEA_DIR / tissue
    if not tissue_dir.exists():
        print(f"  [SKIP] {tissue}: fGSEA directory not found ({tissue_dir})")
        return None

    missions = sorted({
        f.stem.replace("_fgsea_hallmark", "")
        for f in tissue_dir.glob("*_fgsea_hallmark.csv")
    })
    if len(missions) < 2:
        print(f"  [SKIP] {tissue}: need >= 2 missions, found {missions}")
        return None

    nes = {}
    for m in missions:
        csv = tissue_dir / f"{m}_fgsea_hallmark.csv"
        df = pd.read_csv(csv)
        nes[m] = df.set_index("pathway")["NES"]

    pairs = {}
    for m1, m2 in combinations(sorted(nes), 2):
        common = nes[m1].index.intersection(nes[m2].index)
        if len(common) < 10:
            print(f"  [WARN] {tissue}: {m1}/{m2} only {len(common)} common pathways, skipping")
            continue
        r, _ = spearmanr(nes[m1][common], nes[m2][common])
        pairs[f"{m1}_{m2}"] = round(float(r), 3)

    if not pairs:
        return None

    mean_r = round(float(np.mean(list(pairs.values()))), 3)
    transfer = TRANSFER_AUROC.get(tissue)
    if transfer is None:
        print(f"  [WARN] {tissue}: no transfer AUROC data in TRANSFER_AUROC table")
        return None

    entry = {
        "nes_mean_r":          mean_r,
        "n_missions_nes":      len(nes),
        "nes_pairs":           pairs,
        "transfer_auroc":      transfer["auroc"],
        "transfer_ci":         transfer["ci"],
        "n_missions_transfer": transfer["n"],
    }
    if tissue in NOTES:
        entry["note"] = NOTES[tissue]
    return entry


def main():
    parser = argparse.ArgumentParser(
        description="Compute NES Spearman r and update NES_conservation_vs_transfer.json"
    )
    parser.add_argument("--tissue", nargs="+", help="Tissue(s) to compute")
    parser.add_argument("--all", action="store_true", help="Recompute all tissues")
    args = parser.parse_args()

    if not NES_JSON.exists():
        print(f"[ERROR] JSON not found: {NES_JSON}")
        return

    existing = json.loads(NES_JSON.read_text())
    data = existing.get("data", {})

    if args.all:
        tissues = list(TRANSFER_AUROC.keys())
    elif args.tissue:
        tissues = args.tissue
    else:
        print("[ERROR] Specify --tissue <name> or --all")
        return

    for tissue in tissues:
        print(f"\n  Processing: {tissue}")
        result = compute_tissue(tissue)
        if result:
            data[tissue] = result
            print(f"    nes_mean_r = {result['nes_mean_r']:.3f}  "
                  f"(n_pairs={len(result['nes_pairs'])}, "
                  f"n_missions={result['n_missions_nes']})")
            for pair_key, r_val in sorted(result["nes_pairs"].items()):
                print(f"      {pair_key}: r={r_val:.3f}")

    # Update correlation section (rank order: NES mean_r vs transfer AUROC)
    tissues_with_data = {
        t: d for t, d in data.items()
        if isinstance(d, dict) and "nes_mean_r" in d and "transfer_auroc" in d
    }
    sorted_by_nes = sorted(
        tissues_with_data.items(),
        key=lambda x: x[1]["nes_mean_r"],
        reverse=True
    )

    existing["data"] = data
    existing["correlation"] = {
        "tissues_ordered_by_nes": [t for t, _ in sorted_by_nes],
        "nes_values":   [d["nes_mean_r"]     for _, d in sorted_by_nes],
        "auroc_values": [d["transfer_auroc"] for _, d in sorted_by_nes],
        "note": existing.get("correlation", {}).get("note", (
            "Gastrocnemius outlier: only 2 fGSEA missions (RR-5 no DGE). "
            "Excluding gastrocnemius and skin outliers where applicable: "
            "rank-order correlation between NES conservation and transfer AUROC "
            "remains strong (see key_finding)."
        )),
    }
    existing["timestamp"] = datetime.now().isoformat()

    NES_JSON.write_text(json.dumps(existing, indent=2))
    print(f"\n  Updated: {NES_JSON.relative_to(BASE_DIR)}")
    print(f"  Tissues in JSON: {list(data.keys())}")


if __name__ == "__main__":
    main()
