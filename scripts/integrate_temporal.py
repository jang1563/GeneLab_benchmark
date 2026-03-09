#!/usr/bin/env python3
"""Integrate v2 temporal analysis results (T1-T3) into main evaluation framework.

Reads from v2/evaluation/T_temporal_summary.json and individual result files,
outputs summary JSONs to evaluation/ following the established schema pattern.

Usage:
    python scripts/integrate_temporal.py
    python scripts/integrate_temporal.py --dry-run
"""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
V2_EVAL = PROJECT_ROOT / "v2" / "evaluation"
EVAL_DIR = PROJECT_ROOT / "evaluation"


def load_v2_summary():
    """Load the comprehensive T_temporal_summary.json."""
    path = V2_EVAL / "T_temporal_summary.json"
    if not path.exists():
        print(f"ERROR: {path} not found", file=sys.stderr)
        sys.exit(1)
    with open(path) as f:
        return json.load(f)


def build_t1_summary(data):
    """Build T1 ISS-T vs LAR summary JSON."""
    t1 = data["T1"]

    # Extract key metrics for each sub-task
    subtasks = {}
    for key in ["T1a_RR-6_liver", "T1b_RR-8_liver", "T1c_RR-6_thymus"]:
        entry = t1[key]
        conditions = entry["conditions"]
        subtask = {
            "task": entry["task"],
            "tissue": entry["tissue"],
            "mission": entry["mission"],
            "confound_warning": entry.get("confound_warning", ""),
            "conditions": {}
        }
        for cond_key, cond_val in conditions.items():
            subtask["conditions"][cond_key] = {
                "auroc": cond_val["auroc"],
                "ci_low": cond_val["ci_low"],
                "ci_high": cond_val["ci_high"],
                "perm_p": cond_val["perm_p"],
                "n_total": cond_val["n_total"],
                "cv_method": cond_val["cv_method"],
            }
        if "interpretation" in entry:
            subtask["interpretation"] = entry["interpretation"]
        subtasks[key] = subtask

    # Cross-mission transfer
    t1d = t1["T1d_cross_mission"]
    transfers = {}
    for tkey, tval in t1d["transfers"].items():
        transfers[tkey] = {
            "auroc": tval["auroc"],
            "ci_low": tval["ci_low"],
            "ci_high": tval["ci_high"],
            "n_train": tval["n_train"],
            "n_test": tval["n_test"],
        }
    subtasks["T1d_cross_mission"] = {
        "task": "T1d",
        "tissue": "liver",
        "description": "Cross-mission temporal transfer",
        "transfers": transfers,
    }

    return {
        "timestamp": datetime.now().isoformat(),
        "description": "T1: ISS-T vs LAR sacrifice timing classification",
        "design_decisions": ["DD-18"],
        "method": "PCA-LR, RepeatedStratifiedKFold(5x10), gene and pathway features",
        "source": "v2/evaluation/T_temporal_summary.json",
        "subtasks": subtasks,
        "verdict": "Preservation artifact dominates biological timing effect",
        "evidence": (
            "GC AUROC >= FLT AUROC in most cases "
            "(RR-8: GC=0.973 > FLT=0.930, excess=-0.043). "
            "BSL near-perfect separation confirms preservation is primary driver."
        ),
    }


def build_t2_summary(data):
    """Build T2 LAR Recovery summary JSON."""
    t2 = data["T2"]
    missions = {}

    for key in ["T2_RR-6", "T2_RR-8"]:
        entry = t2[key]
        mission = entry["mission"]

        # Summarize pathway recovery
        recovery = entry["pathway_recovery_summary"]
        n_recovering = recovery["n_recovering"]
        n_total = recovery["n_total"]
        mean_recovery = recovery["mean"]

        missions[mission] = {
            "task": "T2",
            "tissue": entry["tissue"],
            "mission": mission,
            "pca_recovery_ratio": entry["pca_recovery_ratio"],
            "pca_variance_explained_top3": entry["pca_variance_explained"],
            "pathway_recovery": {
                "mean": mean_recovery,
                "median": recovery["median"],
                "n_recovering": n_recovering,
                "n_total": n_total,
            },
            "classification_recovery": entry["classification_recovery"],
            "groups": entry["groups"],
        }

    return {
        "timestamp": datetime.now().isoformat(),
        "description": "T2: LAR recovery signature analysis",
        "design_decisions": ["DD-20"],
        "method": "PCA distance ratios + direction-aware pathway recovery (DD-20)",
        "source": "v2/evaluation/T_temporal_summary.json",
        "missions": missions,
        "verdict": "LAR samples show partial-to-strong recovery toward baseline",
        "evidence": (
            "RR-6: PCA R=0.842 (partial), 12/26 recovering. "
            "RR-8: PCA R=0.652 (stronger recovery), 25/27 recovering with overshoot. "
            "Classification: FLT_LAR flight_prob 0.185 (RR-6) / 0.404 (RR-8)."
        ),
    }


def build_t3_summary(data):
    """Build T3 Age x Spaceflight summary JSON."""
    t3 = data["T3"]
    subtasks = t3["subtasks"]

    # Extract key subtask results
    results = {}
    for key in ["T3a_gene", "T3a_pathway",
                 "T3b_FLT_gene", "T3b_FLT_pathway",
                 "T3b_GC_gene", "T3b_GC_pathway",
                 "T3b_VIV_gene", "T3b_VIV_pathway",
                 "T3d_OLD_gene", "T3d_OLD_pathway",
                 "T3d_YNG_gene", "T3d_YNG_pathway"]:
        entry = subtasks[key]
        result = {
            "description": entry["description"],
            "auroc": entry["auroc"],
            "ci_low": entry["ci_low"],
            "ci_high": entry["ci_high"],
            "perm_p": entry["perm_p"],
        }
        # Add sample counts (different keys for different subtasks)
        if "n_old" in entry:
            result["n_old"] = entry["n_old"]
            result["n_yng"] = entry["n_yng"]
        elif "n_flight" in entry:
            result["n_flight"] = entry["n_flight"]
            result["n_gc"] = entry["n_gc"]
            result["age_group"] = entry["age_group"]
        if "condition" in entry:
            result["condition"] = entry["condition"]
        results[key] = result

    # ANOVA summary
    anova = subtasks["T3c_anova"]
    results["T3c_anova"] = {
        "n_pathways": anova["n_pathways"],
        "n_significant_interaction": anova["n_significant_interaction"],
        "n_samples": anova["n_samples"],
        "timing_filter": anova["timing_filter"],
        "min_fdr_interaction": min(
            p["fdr_interaction"] for p in anova["pathways"]
        ),
    }

    return {
        "timestamp": datetime.now().isoformat(),
        "description": "T3: Age x Spaceflight interaction analysis (RR-8 liver)",
        "design_decisions": ["DD-19"],
        "method": "PCA-LR RSKF for classification, two-way ANOVA for pathway interaction",
        "source": "v2/evaluation/T_temporal_summary.json",
        "tissue": t3["tissue"],
        "mission": t3["mission"],
        "results": results,
        "aging_hypothesis": t3["aging_hypothesis"],
        "verdict": "Spaceflight amplifies aging: SUPPORTED",
        "evidence": (
            f"OLD AUROC={t3['aging_hypothesis']['OLD_auroc']:.3f} vs "
            f"YNG AUROC={t3['aging_hypothesis']['YNG_auroc']:.3f} "
            f"(delta=+{t3['aging_hypothesis']['delta']:.3f}). "
            "ANOVA: 0/50 significant interactions at FDR<0.05 (underpowered)."
        ),
    }


def verify_values(t1_json, t2_json, t3_json, v2_data):
    """Verify key values match v2 originals."""
    errors = []

    # T1 checks — both missions
    for task_key, cond_key in [
        ("T1a_RR-6_liver", "FLT_gene"),
        ("T1b_RR-8_liver", "GC_gene"),
        ("T1c_RR-6_thymus", "FLT_gene"),
    ]:
        val = t1_json["subtasks"][task_key]["conditions"][cond_key]["auroc"]
        v2_val = v2_data["T1"][task_key]["conditions"][cond_key]["auroc"]
        if val != v2_val:
            errors.append(f"{task_key} {cond_key} AUROC: {val} != {v2_val}")

    # T1d cross-mission transfer
    t1d_val = t1_json["subtasks"]["T1d_cross_mission"]["transfers"]
    for tkey in t1d_val:
        if t1d_val[tkey]["auroc"] != v2_data["T1"]["T1d_cross_mission"]["transfers"][tkey]["auroc"]:
            errors.append(f"T1d {tkey} AUROC mismatch")

    # T2 checks — both missions
    for mission in ["RR-6", "RR-8"]:
        t2_r = t2_json["missions"][mission]["pca_recovery_ratio"]
        v2_r = v2_data["T2"][f"T2_{mission}"]["pca_recovery_ratio"]
        if t2_r != v2_r:
            errors.append(f"T2 {mission} PCA ratio: {t2_r} != {v2_r}")

    # T3 checks — hypothesis + key subtasks
    t3_old = t3_json["aging_hypothesis"]["OLD_auroc"]
    v2_old = v2_data["T3"]["aging_hypothesis"]["OLD_auroc"]
    if t3_old != v2_old:
        errors.append(f"T3 OLD AUROC: {t3_old} != {v2_old}")

    t3_delta = t3_json["aging_hypothesis"]["delta"]
    v2_delta = v2_data["T3"]["aging_hypothesis"]["delta"]
    if t3_delta != v2_delta:
        errors.append(f"T3 delta: {t3_delta} != {v2_delta}")

    for key in ["T3a_gene", "T3b_FLT_gene", "T3d_YNG_gene"]:
        val = t3_json["results"][key]["auroc"]
        v2_val = v2_data["T3"]["subtasks"][key]["auroc"]
        if val != v2_val:
            errors.append(f"T3 {key} AUROC: {val} != {v2_val}")

    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true",
                        help="Print summary without writing files")
    args = parser.parse_args()

    v2_data = load_v2_summary()
    print(f"Loaded v2 temporal summary (version={v2_data['version']})")

    t1_json = build_t1_summary(v2_data)
    t2_json = build_t2_summary(v2_data)
    t3_json = build_t3_summary(v2_data)

    # Verify
    errors = verify_values(t1_json, t2_json, t3_json, v2_data)
    if errors:
        print("VERIFICATION ERRORS:", file=sys.stderr)
        for e in errors:
            print(f"  - {e}", file=sys.stderr)
        sys.exit(1)
    print("Verification passed: all values match v2 originals")

    # Summary
    print(f"\nT1: {len(t1_json['subtasks'])} subtasks, verdict: {t1_json['verdict']}")
    print(f"T2: {len(t2_json['missions'])} missions, verdict: {t2_json['verdict']}")
    print(f"T3: {len(t3_json['results'])} results, verdict: {t3_json['verdict']}")

    if args.dry_run:
        print("\n[DRY RUN] Would write to:")
        for name in ["T1_temporal_timing_summary.json",
                      "T2_recovery_summary.json",
                      "T3_age_spaceflight_summary.json"]:
            print(f"  {EVAL_DIR / name}")
        return

    # Write output files
    outputs = {
        "T1_temporal_timing_summary.json": t1_json,
        "T2_recovery_summary.json": t2_json,
        "T3_age_spaceflight_summary.json": t3_json,
    }
    for name, data in outputs.items():
        path = EVAL_DIR / name
        with open(path, "w") as f:
            json.dump(data, f, indent=2)
        print(f"Wrote {path} ({path.stat().st_size} bytes)")

    print("\nDone. Run with --dry-run to preview without writing.")


if __name__ == "__main__":
    main()
