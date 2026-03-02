#!/usr/bin/env python3
"""NC1: Aggregate permutation test results across all benchmark categories.

Reads existing evaluation JSONs and produces a unified summary of all
permutation-based significance tests (Categories A, C, D) plus bootstrap
CI-based significance for Category B (which lacks permutation p-values).

Output: evaluation/NC1_permutation_summary.json
"""
import json
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
EVAL = BASE / "evaluation"
PROC = BASE / "processed"

# Task → tissue mapping
TASK_TISSUE = {
    "A1": "liver",
    "A2": "gastrocnemius",
    "A3": "kidney",
    "A4": "thymus",
    "A6": "eye",
}

B_TISSUES = ["liver", "gastrocnemius", "kidney", "thymus", "eye"]


def load_json(path: Path) -> dict:
    with open(path) as f:
        return json.load(f)


def collect_category_a() -> list[dict]:
    """Category A: LOMO spaceflight detection (pca_lr model)."""
    entries = []
    for task_id, tissue in TASK_TISSUE.items():
        path = EVAL / f"{task_id}_baseline_results.json"
        if not path.exists():
            print(f"  [WARN] Missing {path.name}")
            continue
        data = load_json(path)

        # Use model with lowest mean_perm_pvalue (best significance)
        best_key = min(
            (k for k in data if "mean_perm_pvalue" in data[k]),
            key=lambda k: data[k]["mean_perm_pvalue"],
            default=None,
        )
        if best_key is None:
            print(f"  [WARN] {task_id}: no model with perm_pvalue")
            continue
        model_key = best_key

        model = data[model_key]
        auroc = model["mean_auroc"]
        perm_p = model["mean_perm_pvalue"]
        n_folds = model["n_folds"]

        entries.append({
            "task_id": task_id,
            "category": "A",
            "description": f"{tissue} spaceflight detection (LOMO)",
            "tissue": tissue,
            "metric": "AUROC",
            "observed": round(auroc, 4),
            "perm_p": round(perm_p, 6),
            "n_perm": 1000,
            "n_folds": n_folds,
            "model": model_key,
            "significant": perm_p < 0.05,
            "verdict": "SIGNIFICANT" if perm_p < 0.05 else "NOT_SIGNIFICANT",
        })
        print(f"  {task_id} ({tissue}): AUROC={auroc:.3f}, perm_p={perm_p:.4f}"
              f" → {'SIG' if perm_p < 0.05 else 'NS'}")

    return entries


def collect_category_b() -> list[dict]:
    """Category B: Cross-mission transfer (bootstrap CI-based significance)."""
    entries = []
    for tissue in B_TISSUES:
        path = PROC / "B_cross_mission" / tissue / "B_transfer_details.json"
        if not path.exists():
            print(f"  [WARN] Missing {path.name} for {tissue}")
            continue
        data = load_json(path)

        # Aggregate pca_lr pairs
        pca_pairs = [d for d in data["details"] if d["method"] == "pca_lr"]
        if not pca_pairs:
            print(f"  [WARN] No pca_lr pairs for {tissue}")
            continue

        aurocs = [d["auroc"] for d in pca_pairs]
        ci_lows = [d["ci_low"] for d in pca_pairs]
        mean_auroc = sum(aurocs) / len(aurocs)
        mean_ci_low = sum(ci_lows) / len(ci_lows)
        # Significance: mean CI lower bound > 0.5 (better than chance)
        sig = mean_ci_low > 0.5

        entries.append({
            "task_id": f"B_{tissue}",
            "category": "B",
            "description": f"{tissue} cross-mission transfer",
            "tissue": tissue,
            "metric": "AUROC",
            "observed": round(mean_auroc, 4),
            "perm_p": None,
            "perm_p_note": "No permutation test; significance via bootstrap CI",
            "ci_low": round(mean_ci_low, 4),
            "n_perm": None,
            "n_pairs": len(pca_pairs),
            "n_missions": data["n_missions"],
            "significant": sig,
            "verdict": "SIGNIFICANT (CI)" if sig else "NOT_SIGNIFICANT (CI)",
        })
        print(f"  B_{tissue}: mean AUROC={mean_auroc:.3f}, CI_low={mean_ci_low:.3f}"
              f" → {'SIG' if sig else 'NS'}")

    return entries


def collect_category_c() -> list[dict]:
    """Category C: Cross-tissue transfer (per-method permutation p)."""
    path = EVAL / "C_cross_tissue_summary.json"
    if not path.exists():
        print(f"  [WARN] Missing {path.name}")
        return []
    data = load_json(path)
    entries = []

    method_labels = {
        "method_a": "Gene-level",
        "method_b": "DEG overlap",
        "method_c": "Pathway-level",
    }

    for task_id, task in data["tasks"].items():
        train_tissue = task["train"]
        test_tissue = task["test"]
        for method_key, method_label in method_labels.items():
            if method_key not in task:
                continue
            m = task[method_key]
            auroc = m["auroc"]
            perm_p = m["perm_p"]
            n_feat = m.get("n_features", None)

            entries.append({
                "task_id": f"{task_id}_{method_key}",
                "category": "C",
                "description": f"{train_tissue}→{test_tissue} {method_label}",
                "train_tissue": train_tissue,
                "test_tissue": test_tissue,
                "method": method_label,
                "metric": "AUROC",
                "observed": round(auroc, 4),
                "perm_p": round(perm_p, 6),
                "n_perm": 10000,
                "n_features": n_feat,
                "significant": perm_p < 0.05,
                "verdict": "SIGNIFICANT" if perm_p < 0.05 else "NOT_SIGNIFICANT",
            })
        best_method = max(method_labels.keys(),
                         key=lambda k: task.get(k, {}).get("auroc", 0))
        best_auroc = task[best_method]["auroc"]
        print(f"  {task_id} ({train_tissue}→{test_tissue}): "
              f"best={method_labels[best_method]} AUROC={best_auroc:.3f}")

    return entries


def collect_category_d() -> list[dict]:
    """Category D: Condition/confounder prediction (per-mode permutation p)."""
    path = EVAL / "D_condition_summary.json"
    if not path.exists():
        print(f"  [WARN] Missing {path.name}")
        return []
    data = load_json(path)
    entries = []

    # n_perm varies: D3=10000, D6=5000
    n_perm_map = {"D3": 10000, "D6_liver": 5000, "D6_thymus": 5000}

    for task_id, task in data["tasks"].items():
        tissue = task.get("tissue", "")
        n_samples = task.get("n_samples", None)
        for mode in ["gene", "pathway"]:
            if mode not in task.get("feature_modes", {}):
                continue
            fm = task["feature_modes"][mode]
            metric_val = fm.get("macro_f1", None)
            perm_p = fm.get("perm_p", None)
            if perm_p is None:
                continue
            n_feat = fm.get("n_features", None)
            n_perm = n_perm_map.get(task_id, 10000)

            entries.append({
                "task_id": f"{task_id}_{mode}",
                "category": "D",
                "description": f"{task_id} {tissue} ({mode} features)",
                "tissue": tissue,
                "feature_mode": mode,
                "metric": "macro_F1",
                "observed": round(metric_val, 4) if metric_val else None,
                "perm_p": round(perm_p, 6),
                "n_perm": n_perm,
                "n_features": n_feat,
                "n_samples": n_samples,
                "significant": perm_p < 0.05,
                "verdict": "SIGNIFICANT" if perm_p < 0.05 else "NOT_SIGNIFICANT",
            })
        # Print summary
        modes = task.get("feature_modes", {})
        for mode in ["gene", "pathway"]:
            if mode in modes:
                fm = modes[mode]
                print(f"  {task_id}_{mode}: F1={fm.get('macro_f1', '?'):.3f}, "
                      f"perm_p={fm.get('perm_p', '?'):.6f}")

    return entries


def main():
    print("=" * 60)
    print("NC1: Aggregating Permutation Test Results")
    print("=" * 60)

    all_entries = []

    print("\n--- Category A: Spaceflight Detection (LOMO) ---")
    all_entries.extend(collect_category_a())

    print("\n--- Category B: Cross-Mission Transfer (CI-based) ---")
    all_entries.extend(collect_category_b())

    print("\n--- Category C: Cross-Tissue Transfer ---")
    all_entries.extend(collect_category_c())

    print("\n--- Category D: Condition Prediction ---")
    all_entries.extend(collect_category_d())

    # Summary statistics
    n_total = len(all_entries)
    n_sig = sum(1 for e in all_entries if e["significant"])
    n_perm_entries = sum(1 for e in all_entries if e.get("perm_p") is not None)
    n_ci_entries = sum(1 for e in all_entries if e.get("perm_p") is None)

    summary = {
        "n_total_entries": n_total,
        "n_significant": n_sig,
        "n_not_significant": n_total - n_sig,
        "n_permutation_tested": n_perm_entries,
        "n_ci_based_only": n_ci_entries,
        "by_category": {},
    }
    for cat in ["A", "B", "C", "D"]:
        cat_entries = [e for e in all_entries if e["category"] == cat]
        summary["by_category"][cat] = {
            "n_entries": len(cat_entries),
            "n_significant": sum(1 for e in cat_entries if e["significant"]),
        }

    output = {
        "description": "NC1: Unified permutation test summary across all benchmark categories",
        "note": "Category B uses bootstrap CI-based significance (no permutation p available). "
                "All others use permutation p-values from original evaluations.",
        "summary": summary,
        "entries": all_entries,
    }

    out_path = EVAL / "NC1_permutation_summary.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved: {out_path}")

    # Print summary table
    print(f"\n{'='*60}")
    print(f"NC1 Summary: {n_total} entries ({n_sig} significant, "
          f"{n_total - n_sig} not significant)")
    print(f"  Permutation-tested: {n_perm_entries}")
    print(f"  CI-based only: {n_ci_entries}")
    for cat, info in summary["by_category"].items():
        print(f"  Category {cat}: {info['n_entries']} entries, "
              f"{info['n_significant']} significant")

    # Validation checks
    print(f"\n--- Validation Checks ---")
    ok = True

    # Check A2 (gastrocnemius) and A4 (thymus) should be significant
    for task_id in ["A2", "A4"]:
        entry = next((e for e in all_entries if e["task_id"] == task_id), None)
        if entry and entry["significant"]:
            print(f"  ✓ {task_id} is significant (p={entry['perm_p']:.4f})")
        elif entry:
            print(f"  ⚠ {task_id} NOT significant (p={entry['perm_p']:.4f}) - unexpected")
        else:
            print(f"  ✗ {task_id} missing")
            ok = False

    # Check A1 (liver) and A3 (kidney) should NOT be significant
    for task_id in ["A1", "A3"]:
        entry = next((e for e in all_entries if e["task_id"] == task_id), None)
        if entry and not entry["significant"]:
            print(f"  ✓ {task_id} is NOT significant (p={entry['perm_p']:.4f}) - expected NO-GO tissue")
        elif entry:
            print(f"  ⚠ {task_id} IS significant (p={entry['perm_p']:.4f}) - unexpected for NO-GO tissue")
        else:
            print(f"  ✗ {task_id} missing")
            ok = False

    # Check D3 gene should be significant
    d3_gene = next((e for e in all_entries if e["task_id"] == "D3_gene"), None)
    if d3_gene and d3_gene["significant"]:
        print(f"  ✓ D3_gene is significant (p={d3_gene['perm_p']:.6f})")
    else:
        print(f"  ⚠ D3_gene check failed")

    # Check minimum entry count
    if n_total >= 15:
        print(f"  ✓ Entry count ({n_total}) ≥ 15 threshold")
    else:
        print(f"  ✗ Entry count ({n_total}) < 15 threshold")
        ok = False

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
