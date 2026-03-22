#!/usr/bin/env python3
"""
aggregate_nc_ablation.py — GeneLabBench v4: Aggregate NC + Ablation Results

Collects all NC1, NC2, D3, D5, and ablation results into summary JSONs.

Usage:
  python aggregate_nc_ablation.py [--eval-dir v4/evaluation/]
"""

import json
import argparse
import numpy as np
from pathlib import Path
from datetime import datetime

import sys
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from v4_utils import TISSUE_MISSIONS, V4_EVAL_DIR
from classifier_registry import CLASSIFIERS


def load_json(path):
    """Load JSON file, return None if missing."""
    if path.exists():
        return json.loads(path.read_text())
    return None


def aggregate_nc(eval_dir):
    """Aggregate NC1 + NC2 results."""
    tissues = list(TISSUE_MISSIONS.keys())
    methods = list(CLASSIFIERS.keys())

    nc1_results = {}
    nc2_results = {}

    for tissue in tissues:
        nc1_results[tissue] = {}
        nc2_results[tissue] = {}
        for method in methods:
            # NC1
            r = load_json(eval_dir / f"NC1_shuffled_{tissue}_{method}.json")
            if r:
                nc1_results[tissue][method] = {
                    "mean_auroc": r["mean_auroc"],
                    "std_auroc": r["std_auroc"],
                    "n_folds": r["n_folds"],
                    "mean_perm_pvalue": r.get("mean_perm_pvalue"),
                }

            # NC2
            r = load_json(eval_dir / f"NC2_random_{tissue}_{method}.json")
            if r:
                nc2_results[tissue][method] = {
                    "mean_auroc": r["mean_auroc"],
                    "std_auroc": r["std_auroc"],
                    "n_folds": r["n_folds"],
                    "mean_perm_pvalue": r.get("mean_perm_pvalue"),
                }

    # Validation
    nc1_pass = 0
    nc1_warn = 0
    nc1_total = 0
    nc2_pass = 0
    nc2_warn = 0
    nc2_total = 0

    for tissue in tissues:
        for method in methods:
            if method in nc1_results[tissue]:
                nc1_total += 1
                auroc = nc1_results[tissue][method]["mean_auroc"]
                if 0.35 <= auroc <= 0.65:
                    nc1_pass += 1
                else:
                    nc1_warn += 1
                    print(f"  NC1 WARN: {tissue}/{method} AUROC={auroc:.3f}")

            if method in nc2_results[tissue]:
                nc2_total += 1
                auroc = nc2_results[tissue][method]["mean_auroc"]
                if 0.35 <= auroc <= 0.65:
                    nc2_pass += 1
                else:
                    nc2_warn += 1
                    print(f"  NC2 WARN: {tissue}/{method} AUROC={auroc:.3f}")

    summary = {
        "experiment": "negative_controls_summary",
        "NC1_shuffled_labels": {
            "total": nc1_total,
            "pass": nc1_pass,
            "warn": nc1_warn,
            "results": nc1_results,
        },
        "NC2_random_features": {
            "total": nc2_total,
            "pass": nc2_pass,
            "warn": nc2_warn,
            "results": nc2_results,
        },
        "timestamp": datetime.now().isoformat(),
    }

    print(f"\nNC1: {nc1_pass}/{nc1_total} PASS, {nc1_warn} WARN")
    print(f"NC2: {nc2_pass}/{nc2_total} PASS, {nc2_warn} WARN")

    return summary


def aggregate_condition(eval_dir):
    """Aggregate D3 + D5 results."""
    d3_results = {}
    d5_results = {}

    # D3: 6 LOMO tissues × 8 methods × 2 features
    d3_tissues = ["liver", "gastrocnemius", "kidney", "thymus", "eye", "skin"]
    methods = list(CLASSIFIERS.keys())
    features = ["gene", "pathway_hallmark"]

    for tissue in d3_tissues:
        d3_results[tissue] = {}
        for method in methods:
            d3_results[tissue][method] = {}
            for feat in features:
                r = load_json(eval_dir / f"D3_v4_{tissue}_{method}_{feat}.json")
                if r:
                    d3_results[tissue][method][feat] = {
                        "macro_f1": r["macro_f1"],
                        "accuracy": r["accuracy"],
                        "perm_pvalue": r["perm_pvalue"],
                        "n_classes": r["n_classes"],
                    }

    # D5: liver, thymus × 8 methods × 2 features
    d5_tissues = ["liver", "thymus"]
    for tissue in d5_tissues:
        d5_results[tissue] = {}
        for method in methods:
            d5_results[tissue][method] = {}
            for feat in features:
                r = load_json(eval_dir / f"D5_v4_{tissue}_{method}_{feat}.json")
                if r:
                    d5_results[tissue][method][feat] = {
                        "macro_f1": r["macro_f1"],
                        "accuracy": r["accuracy"],
                        "perm_pvalue": r["perm_pvalue"],
                    }

    # J5: gene vs pathway comparison
    j5 = {}
    for tissue in d3_tissues:
        j5[tissue] = {}
        for method in methods:
            gene_r = d3_results.get(tissue, {}).get(method, {}).get("gene")
            path_r = d3_results.get(tissue, {}).get(method, {}).get("pathway_hallmark")
            if gene_r and path_r:
                j5[tissue][method] = {
                    "gene_f1": gene_r["macro_f1"],
                    "pathway_f1": path_r["macro_f1"],
                    "diff": round(path_r["macro_f1"] - gene_r["macro_f1"], 4),
                    "winner": "pathway" if path_r["macro_f1"] > gene_r["macro_f1"] else "gene",
                }

    summary = {
        "experiment": "condition_prediction_summary",
        "D3_mission_id": d3_results,
        "D5_hardware": d5_results,
        "J5_gene_vs_pathway": j5,
        "timestamp": datetime.now().isoformat(),
    }

    # Print summary table
    print("\nD3 Mission ID Prediction (macro-F1):")
    print(f"  {'tissue':<15s}", end="")
    for m in methods:
        print(f" {m:>12s}", end="")
    print()
    for tissue in d3_tissues:
        print(f"  {tissue:<15s}", end="")
        for m in methods:
            r = d3_results.get(tissue, {}).get(m, {}).get("gene")
            val = f"{r['macro_f1']:.3f}" if r else "  ---"
            print(f" {val:>12s}", end="")
        print()

    return summary


def aggregate_ablation(eval_dir):
    """Aggregate all ablation results."""
    tissues = list(TISSUE_MISSIONS.keys())

    # Feature count ablation
    feat_results = {}
    top_k_values = [100, 500, 1000, 2000, 5000, 10000, "all"]
    for tissue in tissues:
        feat_results[tissue] = {}
        for method in ["pca_lr", "elasticnet_lr"]:
            feat_results[tissue][method] = {}
            for k in top_k_values:
                k_str = str(k)
                r = load_json(eval_dir / f"ABL_feat_{tissue}_{method}_k{k_str}.json")
                if r:
                    feat_results[tissue][method][k_str] = {
                        "mean_auroc": r["mean_auroc"],
                        "std_auroc": r["std_auroc"],
                        "top_k_actual": r.get("top_k_actual"),
                    }

    # PCA component ablation
    pca_results = {}
    nc_values = [5, 10, 20, 50, 100, 200]
    for tissue in tissues:
        pca_results[tissue] = {}
        for nc in nc_values:
            r = load_json(eval_dir / f"ABL_pca_{tissue}_nc{nc}.json")
            if r:
                pca_results[tissue][str(nc)] = {
                    "mean_auroc": r["mean_auroc"],
                    "std_auroc": r["std_auroc"],
                    "mean_explained_var": r.get("mean_explained_var"),
                }

    # Sample size ablation
    sample_results = {}
    sample_tissues = ["liver", "thymus", "skin", "gastrocnemius"]
    fractions = [0.2, 0.4, 0.6, 0.8, 1.0]
    for tissue in sample_tissues:
        sample_results[tissue] = {}
        for method in ["pca_lr", "elasticnet_lr"]:
            sample_results[tissue][method] = {}
            for frac in fractions:
                frac_str = f"{frac:.1f}".replace(".", "")
                repeat_aurocs = []
                for rep in range(10):
                    r = load_json(eval_dir /
                                  f"ABL_sample_{tissue}_{method}_f{frac_str}_r{rep}.json")
                    if r:
                        repeat_aurocs.append(r["mean_auroc"])

                if repeat_aurocs:
                    sample_results[tissue][method][str(frac)] = {
                        "mean_auroc": round(float(np.mean(repeat_aurocs)), 4),
                        "std_auroc": round(float(np.std(repeat_aurocs)), 4),
                        "n_repeats": len(repeat_aurocs),
                    }

    # Bootstrap ablation
    boot_results = {}
    for tissue in tissues:
        r = load_json(eval_dir / f"ABL_boot_{tissue}.json")
        if r:
            boot_results[tissue] = {}
            for level_key, level_data in r.get("levels", {}).items():
                boot_results[tissue][level_key] = {
                    "mean_auroc": level_data["mean_auroc"],
                    "mean_ci_width": level_data.get("mean_ci_width"),
                }

    summary = {
        "experiment": "ablation_summary",
        "S1_feature_count": feat_results,
        "S2_pca_components": pca_results,
        "S3_sample_size": sample_results,
        "S4_bootstrap_stability": boot_results,
        "timestamp": datetime.now().isoformat(),
    }

    # Print summary
    print("\nFeature Count Ablation (PCA-LR mean AUROC):")
    print(f"  {'tissue':<15s}", end="")
    for k in top_k_values:
        print(f" {'k='+str(k):>8s}", end="")
    print()
    for tissue in tissues:
        print(f"  {tissue:<15s}", end="")
        for k in top_k_values:
            r = feat_results.get(tissue, {}).get("pca_lr", {}).get(str(k))
            val = f"{r['mean_auroc']:.3f}" if r else "  ---"
            print(f" {val:>8s}", end="")
        print()

    print("\nPCA Components Ablation (mean AUROC):")
    print(f"  {'tissue':<15s}", end="")
    for nc in nc_values:
        print(f" {'nc='+str(nc):>8s}", end="")
    print()
    for tissue in tissues:
        print(f"  {tissue:<15s}", end="")
        for nc in nc_values:
            r = pca_results.get(tissue, {}).get(str(nc))
            val = f"{r['mean_auroc']:.3f}" if r else "  ---"
            print(f" {val:>8s}", end="")
        print()

    print("\nBootstrap Stability (mean CI width):")
    for tissue in tissues:
        if tissue in boot_results:
            widths = []
            for level_key in ["100", "500", "1000", "2000", "5000"]:
                w = boot_results[tissue].get(level_key, {}).get("mean_ci_width")
                widths.append(f"{w:.3f}" if w else "---")
            print(f"  {tissue:<15s}: {' → '.join(widths)}")

    return summary


def main():
    parser = argparse.ArgumentParser(description="Aggregate NC + Ablation Results")
    parser.add_argument("--eval-dir", type=str, default=None)
    args = parser.parse_args()

    eval_dir = Path(args.eval_dir) if args.eval_dir else V4_EVAL_DIR

    print("=" * 70)
    print("GeneLabBench v4: NC + Ablation Aggregation")
    print("=" * 70)

    # NC
    print("\n── Negative Controls ──")
    nc_summary = aggregate_nc(eval_dir)
    (eval_dir / "NC_summary.json").write_text(json.dumps(nc_summary, indent=2))
    print(f"  Saved: {eval_dir / 'NC_summary.json'}")

    # Condition prediction
    print("\n── Condition Prediction (D3/D5) ──")
    cond_summary = aggregate_condition(eval_dir)
    (eval_dir / "D_condition_v4_summary.json").write_text(
        json.dumps(cond_summary, indent=2))
    print(f"  Saved: {eval_dir / 'D_condition_v4_summary.json'}")

    # Ablation
    print("\n── Ablation Studies ──")
    abl_summary = aggregate_ablation(eval_dir)
    (eval_dir / "ABL_summary.json").write_text(json.dumps(abl_summary, indent=2))
    print(f"  Saved: {eval_dir / 'ABL_summary.json'}")

    # Regression tests
    print("\n── Regression Tests ──")
    m1_ref = {}
    for tissue in TISSUE_MISSIONS:
        r = load_json(eval_dir / f"M1_{tissue}_gene_pca_lr.json")
        if r:
            m1_ref[tissue] = r["mean_auroc"]

    # Check feature K=all matches M1
    feat_check = abl_summary.get("S1_feature_count", {})
    for tissue in TISSUE_MISSIONS:
        if tissue not in m1_ref:
            continue
        r = feat_check.get(tissue, {}).get("pca_lr", {}).get("all")
        if r:
            diff = abs(r["mean_auroc"] - m1_ref[tissue])
            status = "PASS" if diff < 0.002 else "FAIL"
            print(f"  K=all {tissue:<15s}: {r['mean_auroc']:.4f} vs M1={m1_ref[tissue]:.4f} "
                  f"(diff={diff:.4f}) {status}")

    # Check PCA nc=50 matches M1
    pca_check = abl_summary.get("S2_pca_components", {})
    for tissue in TISSUE_MISSIONS:
        if tissue not in m1_ref:
            continue
        r = pca_check.get(tissue, {}).get("50")
        if r:
            diff = abs(r["mean_auroc"] - m1_ref[tissue])
            status = "PASS" if diff < 0.002 else "FAIL"
            print(f"  nc=50 {tissue:<15s}: {r['mean_auroc']:.4f} vs M1={m1_ref[tissue]:.4f} "
                  f"(diff={diff:.4f}) {status}")

    # Check bootstrap AUROC matches M1
    boot_check = abl_summary.get("S4_bootstrap_stability", {})
    for tissue in TISSUE_MISSIONS:
        if tissue not in m1_ref:
            continue
        r = boot_check.get(tissue, {}).get("2000")
        if r:
            diff = abs(r["mean_auroc"] - m1_ref[tissue])
            status = "PASS" if diff < 0.002 else "FAIL"
            print(f"  boot  {tissue:<15s}: {r['mean_auroc']:.4f} vs M1={m1_ref[tissue]:.4f} "
                  f"(diff={diff:.4f}) {status}")

    print(f"\nDone. Results in {eval_dir}")


if __name__ == "__main__":
    main()
