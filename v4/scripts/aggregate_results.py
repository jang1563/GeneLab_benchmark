#!/usr/bin/env python3
"""Aggregate and summarize all 256 v4 multi-method evaluation results."""

import json
import numpy as np
from pathlib import Path

eval_dir = Path(__file__).resolve().parent.parent / "evaluation"
tissues = ["liver", "gastrocnemius", "kidney", "thymus", "eye", "skin", "lung", "colon"]
methods = ["elasticnet_lr", "pca_lr", "rf", "xgb", "svm_rbf", "knn", "mlp", "tabnet"]
features = ["gene", "pathway_hallmark", "pathway_kegg", "combined"]
method_labels = {"elasticnet_lr": "EN-LR", "pca_lr": "PCA-LR", "rf": "RF", "xgb": "XGB",
                 "svm_rbf": "SVM", "knn": "kNN", "mlp": "MLP", "tabnet": "TabNet"}

# Load all results
all_results = {}
for t in tissues:
    for f in features:
        for m in methods:
            fname = f"M1_{t}_{f}_{m}.json"
            fp = eval_dir / fname
            if fp.exists():
                data = json.loads(fp.read_text())
                all_results[(t, f, m)] = data

print(f"Loaded: {len(all_results)}/256 evaluations")
print()

# ==== TABLE 1: Gene features ====
print("TABLE 1: Gene Features — AUROC")
mlabs = [method_labels[m] for m in methods]
header = f"{'Tissue':15s}" + "".join(f"{l:>9s}" for l in mlabs)
print(header)
print("-" * len(header))
for t in tissues:
    row = f"{t:15s}"
    for m in methods:
        key = (t, "gene", m)
        if key in all_results:
            a = all_results[key]["mean_auroc"]
            p = all_results[key]["mean_perm_pvalue"]
            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
            row += f"  {a:.3f}{sig:<2s}"
        else:
            row += "      N/A"
    print(row)

print("-" * len(header))
for label, tset in [("Mean(LOMO6)", tissues[:6]), ("Mean(all8)", tissues)]:
    row = f"{label:15s}"
    for m in methods:
        vals = [all_results[(t,"gene",m)]["mean_auroc"] for t in tset if (t,"gene",m) in all_results]
        row += f"  {np.mean(vals):.3f}  "
    print(row)

# ==== TABLE 2: Pathway comparison ====
print(f"\nTABLE 2: Gene vs Pathway (PCA-LR)")
print(f"{'Tissue':15s}{'Gene':>8s}{'Hallmark':>10s}{'KEGG':>8s}{'Combined':>10s}  {'Best':>8s}")
print("-" * 62)
for t in tissues:
    row = f"{t:15s}"
    best_f, best_a = "gene", 0
    for f in features:
        key = (t, f, "pca_lr")
        if key in all_results:
            a = all_results[key]["mean_auroc"]
            row += f"{a:10.3f}"
            if a > best_a:
                best_a = a
                best_f = f
        else:
            row += "       N/A"
    row += f"  {best_f:>8s}"
    print(row)

# ==== TABLE 3: Best overall per tissue ====
print(f"\nTABLE 3: Best Method × Feature per Tissue")
print(f"{'Tissue':15s}{'Method':>12s}{'Feature':>18s}{'AUROC':>8s}{'p':>8s}")
print("-" * 65)
for t in tissues:
    best_key, best_auroc = None, -1
    for f in features:
        for m in methods:
            key = (t, f, m)
            if key in all_results:
                a = all_results[key]["mean_auroc"]
                if a > best_auroc:
                    best_auroc = a
                    best_key = key
    if best_key:
        p = all_results[best_key]["mean_perm_pvalue"]
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
        t2, f2, m2 = best_key
        print(f"{t2:15s}{method_labels[m2]:>12s}{f2:>18s}  {best_auroc:.4f}{sig:<3s}  {p:.4f}")

# ==== Verification ====
print(f"\n==== VERIFICATION ====")
leaks = [k for k in all_results if all_results[k]["mean_auroc"] > 0.99]
print(f"Data leakage (AUROC>0.99): {len(leaks)} {'PASS' if not leaks else 'FAIL!'}")

thymus_sig = sum(1 for m in methods if ("thymus","gene",m) in all_results
                 and all_results[("thymus","gene",m)]["mean_perm_pvalue"] < 0.05)
print(f"Thymus sig methods (gene, p<0.05): {thymus_sig} {'PASS' if thymus_sig >= 2 else 'CHECK'}")

n_sig = sum(1 for k in all_results if all_results[k]["mean_perm_pvalue"] < 0.05)
print(f"Total significant (p<0.05): {n_sig}/{len(all_results)}")

# ==== Save aggregated summary ====
summary = {}
for t in tissues:
    summary[t] = {}
    for f in features:
        summary[t][f] = {}
        for m in methods:
            key = (t, f, m)
            if key in all_results:
                d = all_results[key]
                summary[t][f][m] = {
                    "auroc": d["mean_auroc"],
                    "std": d["std_auroc"],
                    "ci_lower": d["mean_ci_lower"],
                    "perm_p": d["mean_perm_pvalue"],
                    "n_folds": d["n_folds"],
                    "cv": d["cv_strategy"],
                    "wall_time": d["wall_time_sec"],
                }

out_file = eval_dir / "M1_summary.json"
out_file.write_text(json.dumps(summary, indent=2))
print(f"\nSaved: {out_file}")
