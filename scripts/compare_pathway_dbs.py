#!/usr/bin/env python3
"""
compare_pathway_dbs.py — Multi-DB pathway comparison analysis for V1 paper.

Generates:
  1. AUROC comparison table: tissue × DB (6×4)
  2. NES conservation comparison: tissue × DB (6×4)
  3. MitoCarta deep dive: OXPHOS Complex I-V per tissue
  4. Cross-DB concordance analysis

Usage:
    python scripts/compare_pathway_dbs.py
"""

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent.parent
EVAL_DIR = ROOT / "evaluation"
FGSEA_DIR = ROOT / "processed" / "fgsea"

TISSUES = ["liver", "gastrocnemius", "kidney", "thymus", "eye", "skin"]
TASK_IDS = {"liver": "A1", "gastrocnemius": "A2", "kidney": "A3",
            "thymus": "A4", "skin": "A5", "eye": "A6"}
DBS = ["hallmark", "kegg", "reactome", "mitocarta"]


# ── 1. AUROC Comparison ──────────────────────────────────────────────────────

def collect_auroc():
    """Collect PCA-LR and LR mean AUROC across tissues × DBs."""
    results = {}
    for model in ["pca_lr", "lr"]:
        table = {}
        for tissue in TISSUES:
            task_id = TASK_IDS[tissue]
            row = {}
            for db in DBS:
                fpath = EVAL_DIR / f"{task_id}_pathway_{db}_results.json"
                if fpath.exists():
                    data = json.load(open(fpath))
                    row[db] = {
                        "mean_auroc": data[model]["mean_auroc"],
                        "std_auroc": data[model]["std_auroc"],
                        "mean_ci_lower": data[model]["mean_ci_lower"],
                        "mean_perm_p": data[model]["mean_perm_p"],
                        "go": data[model]["go"],
                        "n_folds": data[model]["n_folds"],
                    }
                else:
                    row[db] = None
            table[tissue] = row
        results[model] = table
    return results


# ── 2. NES Conservation ──────────────────────────────────────────────────────

def collect_nes_conservation():
    """Collect NES conservation (nes_mean_r) across tissues × DBs."""
    table = {}
    for tissue in TISSUES:
        row = {}
        for db in DBS:
            fpath = EVAL_DIR / f"NES_conservation_{db}.json"
            if fpath.exists():
                data = json.load(open(fpath))
                if tissue in data["data"]:
                    entry = data["data"][tissue]
                    row[db] = {
                        "nes_mean_r": entry.get("nes_mean_r"),
                        "n_missions": entry.get("n_missions_nes"),
                        "n_pathways": entry.get("n_pathways"),
                    }
                else:
                    row[db] = None
            else:
                row[db] = None
        table[tissue] = row
    return table


# ── 3. MitoCarta Deep Dive ────────────────────────────────────────────────────

def mitocarta_deep_dive():
    """Analyze MitoCarta OXPHOS sub-complex NES across tissues."""
    oxphos_keywords = {
        "Complex_I": ["COMPLEX_I", "NADH"],
        "Complex_II": ["COMPLEX_II", "SUCCINATE"],
        "Complex_III": ["COMPLEX_III", "CYTOCHROME_B"],
        "Complex_IV": ["COMPLEX_IV", "CYTOCHROME_C_OXIDASE"],
        "Complex_V": ["COMPLEX_V", "ATP_SYNTHASE"],
        "OXPHOS_general": ["OXPHOS", "OXIDATIVE_PHOSPHORYLATION"],
    }

    results = {}
    for tissue in TISSUES:
        tissue_dir = FGSEA_DIR / tissue
        mitocarta_files = sorted(tissue_dir.glob("*_fgsea_mitocarta.csv"))
        if not mitocarta_files:
            continue

        tissue_data = {}
        for fpath in mitocarta_files:
            mission = fpath.stem.replace("_fgsea_mitocarta", "")
            df = pd.read_csv(fpath)
            if "pathway" not in df.columns:
                continue

            mission_data = {}
            for complex_name, keywords in oxphos_keywords.items():
                mask = df["pathway"].str.upper().apply(
                    lambda p: any(k in p.upper() for k in keywords)
                )
                matching = df[mask]
                if not matching.empty:
                    mission_data[complex_name] = {
                        "pathways": matching["pathway"].tolist(),
                        "NES": matching["NES"].tolist(),
                        "padj": matching["padj"].tolist(),
                        "mean_NES": float(matching["NES"].mean()),
                    }
            tissue_data[mission] = mission_data
        results[tissue] = tissue_data
    return results


# ── 4. Cross-DB Concordance ───────────────────────────────────────────────────

def cross_db_concordance():
    """Compare NES direction concordance between Hallmark OXPHOS and MitoCarta."""
    concordance = {}
    for tissue in TISSUES:
        tissue_dir = FGSEA_DIR / tissue
        hallmark_files = sorted(tissue_dir.glob("*_fgsea_hallmark.csv"))
        mitocarta_files = sorted(tissue_dir.glob("*_fgsea_mitocarta.csv"))

        hallmark_missions = {f.stem.replace("_fgsea_hallmark", ""): f for f in hallmark_files}
        mitocarta_missions = {f.stem.replace("_fgsea_mitocarta", ""): f for f in mitocarta_files}

        common_missions = sorted(set(hallmark_missions) & set(mitocarta_missions))
        if not common_missions:
            continue

        tissue_concordance = {}
        for mission in common_missions:
            hm_df = pd.read_csv(hallmark_missions[mission])
            mc_df = pd.read_csv(mitocarta_missions[mission])

            # Find HALLMARK_OXIDATIVE_PHOSPHORYLATION
            hm_oxphos = hm_df[hm_df["pathway"].str.contains("OXIDATIVE_PHOSPHORYLATION", case=False)]
            if hm_oxphos.empty:
                continue
            hm_nes = float(hm_oxphos["NES"].iloc[0])

            # Find all MitoCarta OXPHOS-related
            mc_oxphos = mc_df[mc_df["pathway"].str.upper().str.contains(
                "OXPHOS|COMPLEX_I|COMPLEX_II|COMPLEX_III|COMPLEX_IV|COMPLEX_V|NADH|ATP_SYNTHASE|SUCCINATE|CYTOCHROME"
            )]
            if mc_oxphos.empty:
                continue

            # Direction concordance
            mc_directions = (mc_oxphos["NES"] > 0).astype(int)
            hm_direction = 1 if hm_nes > 0 else 0
            n_concordant = int((mc_directions == hm_direction).sum())
            n_total = len(mc_directions)

            tissue_concordance[mission] = {
                "hallmark_OXPHOS_NES": hm_nes,
                "mitocarta_OXPHOS_pathways": mc_oxphos["pathway"].tolist(),
                "mitocarta_NES": mc_oxphos["NES"].tolist(),
                "n_concordant": n_concordant,
                "n_total": n_total,
                "concordance_rate": round(n_concordant / n_total, 3) if n_total > 0 else None,
            }
        concordance[tissue] = tissue_concordance
    return concordance


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("=== Multi-DB Pathway Comparison Analysis ===\n")

    # 1. AUROC
    print("1. Collecting AUROC results...")
    auroc = collect_auroc()

    # 2. NES Conservation
    print("2. Collecting NES conservation...")
    nes_conservation = collect_nes_conservation()

    # 3. MitoCarta Deep Dive
    print("3. Running MitoCarta deep dive...")
    mitocarta = mitocarta_deep_dive()

    # 4. Cross-DB Concordance
    print("4. Computing cross-DB concordance...")
    concordance = cross_db_concordance()

    # Save combined output
    output = {
        "auroc": auroc,
        "nes_conservation": nes_conservation,
        "mitocarta_deep_dive": mitocarta,
        "cross_db_concordance": concordance,
    }

    out_path = EVAL_DIR / "multi_db_comparison.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nSaved → {out_path}")

    # Print summary tables
    print("\n" + "=" * 70)
    print("AUROC Summary (PCA-LR)")
    print("=" * 70)
    print(f"{'Tissue':<15}", end="")
    for db in DBS:
        print(f"{db:>12}", end="")
    print(f"{'best_db':>12}")
    print("-" * 75)
    for tissue in TISSUES:
        print(f"{tissue:<15}", end="")
        best_db, best_auroc = None, -1
        for db in DBS:
            entry = auroc["pca_lr"][tissue].get(db)
            if entry:
                val = entry["mean_auroc"]
                print(f"{val:>12.3f}", end="")
                if val > best_auroc:
                    best_auroc = val
                    best_db = db
            else:
                print(f"{'---':>12}", end="")
        print(f"{best_db:>12}")

    print("\n" + "=" * 70)
    print("NES Conservation (nes_mean_r)")
    print("=" * 70)
    print(f"{'Tissue':<15}", end="")
    for db in DBS:
        print(f"{db:>12}", end="")
    print()
    print("-" * 63)
    for tissue in TISSUES:
        print(f"{tissue:<15}", end="")
        for db in DBS:
            entry = nes_conservation[tissue].get(db)
            if entry and entry.get("nes_mean_r") is not None:
                print(f"{entry['nes_mean_r']:>12.3f}", end="")
            else:
                print(f"{'---':>12}", end="")
        print()

    # MitoCarta summary
    print("\n" + "=" * 70)
    print("MitoCarta OXPHOS Deep Dive")
    print("=" * 70)
    for tissue in TISSUES:
        if tissue not in mitocarta:
            continue
        print(f"\n  {tissue.upper()}:")
        for mission, complexes in mitocarta[tissue].items():
            oxphos_vals = []
            for cname in ["Complex_I", "Complex_II", "Complex_III", "Complex_IV", "Complex_V"]:
                if cname in complexes:
                    nes = complexes[cname]["mean_NES"]
                    oxphos_vals.append(f"{cname}={nes:+.2f}")
            if oxphos_vals:
                print(f"    {mission}: {', '.join(oxphos_vals)}")

    # Concordance summary
    print("\n" + "=" * 70)
    print("Hallmark OXPHOS vs MitoCarta Concordance")
    print("=" * 70)
    for tissue in TISSUES:
        if tissue not in concordance or not concordance[tissue]:
            continue
        rates = [v["concordance_rate"] for v in concordance[tissue].values() if v.get("concordance_rate") is not None]
        if rates:
            mean_rate = np.mean(rates)
            print(f"  {tissue:<15} mean concordance: {mean_rate:.1%} ({len(rates)} missions)")


if __name__ == "__main__":
    main()
