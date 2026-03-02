#!/usr/bin/env python3
"""
gene_vs_pathway_comparison.py — GeneLab_benchmark: Category J5

Systematic comparison of gene-level vs pathway-level feature representations
across all benchmark categories.

Category A (Spaceflight Detection, LOMO):
  - Gene: log2 normalized counts → PCA → LR → AUROC
  - Pathway: GSVA Hallmark scores (50) → LR → AUROC
  (FRESH evaluation with consistent methodology)

Category C (Cross-Tissue Transfer):
  - Gene: Method A (common genes → PCA → LR)
  - Pathway: Method C (fGSEA top-20 → GSVA → LR)
  (EXISTING results from cross_tissue_transfer.py)

Category D (Condition Prediction):
  - Gene: log2 norm → PCA → LR
  - Pathway: GSVA Hallmark → LR
  (EXISTING results from condition_prediction.py)

Output:
  evaluation/J5_gene_vs_pathway.json   # Full comparison table
  (Also prints formatted table for paper)

Usage:
  python scripts/gene_vs_pathway_comparison.py
  python scripts/gene_vs_pathway_comparison.py --skip-category-a  # use existing only
"""

import json
import argparse
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")

from utils import (
    load_metadata, load_gene_features, load_pathway_features,
    align_features_with_meta, TISSUE_MISSIONS,
)

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent.parent
PROCESSED_DIR = BASE_DIR / "processed" / "A_detection"
PATHWAY_DIR = BASE_DIR / "processed" / "pathway_scores"
RESULTS_DIR = BASE_DIR / "evaluation"

# ── Config ─────────────────────────────────────────────────────────────────────
FLIGHT_LABEL = "Flight"
GROUND_LABELS = {"GC", "VC"}
VARIANCE_PERCENTILE = 0.25


# ── Category A: Spaceflight Detection LOMO ──────────────────────────────────

def lomo_auroc(tissue, feature_mode="gene"):
    """Run LOMO evaluation for a tissue. Returns mean AUROC and per-fold AUROCs."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    from sklearn.metrics import roc_auc_score

    meta = load_metadata(tissue)
    if feature_mode == "gene":
        feat = load_gene_features(tissue)
    else:
        feat = load_pathway_features(tissue)

    if feat is None:
        return None, {}

    feat_aligned, meta_aligned = align_features_with_meta(feat, meta)

    # Binary labels: Flight=1, Ground=0
    labels_raw = meta_aligned["label"]
    binary = pd.Series(np.nan, index=meta_aligned.index)
    binary[labels_raw == FLIGHT_LABEL] = 1
    binary[labels_raw.isin(GROUND_LABELS)] = 0
    valid = binary.notna()
    feat_aligned = feat_aligned[valid]
    meta_aligned = meta_aligned[valid]
    binary = binary[valid].astype(int)

    missions = meta_aligned["mission"].values
    unique_missions = sorted(set(missions))

    fold_aurocs = {}
    for test_mission in unique_missions:
        test_mask = missions == test_mission
        train_mask = ~test_mask

        y_train = binary.values[train_mask]
        y_test = binary.values[test_mask]

        if len(np.unique(y_test)) < 2 or len(np.unique(y_train)) < 2:
            continue
        if sum(test_mask) < 3:
            continue

        X_train = feat_aligned.values[train_mask].astype(float)
        X_test = feat_aligned.values[test_mask].astype(float)

        # Handle NaN
        X_train = np.nan_to_num(X_train, nan=0.0)
        X_test = np.nan_to_num(X_test, nan=0.0)

        if feature_mode == "gene":
            # Variance filter (train only)
            var = np.var(X_train, axis=0)
            threshold = np.percentile(var, VARIANCE_PERCENTILE * 100)
            keep = var > threshold
            if keep.sum() < 50:
                keep = np.ones(X_train.shape[1], dtype=bool)
            X_train = X_train[:, keep]
            X_test = X_test[:, keep]

            # Scale + PCA
            scaler = StandardScaler()
            X_train = scaler.fit_transform(X_train)
            X_test = scaler.transform(X_test)
            n_comp = min(50, X_train.shape[0] - 1, X_train.shape[1])
            pca = PCA(n_components=n_comp, random_state=42)
            X_train = pca.fit_transform(X_train)
            X_test = pca.transform(X_test)
        else:
            # Pathway: just scale
            scaler = StandardScaler()
            X_train = scaler.fit_transform(X_train)
            X_test = scaler.transform(X_test)

        clf = LogisticRegression(
            solver="lbfgs", class_weight="balanced",
            max_iter=2000, C=1.0, random_state=42,
        )
        clf.fit(X_train, y_train)

        try:
            y_score = clf.predict_proba(X_test)[:, 1]
            auroc = roc_auc_score(y_test, y_score)
        except Exception:
            auroc = np.nan

        fold_aurocs[test_mission] = auroc

    if fold_aurocs:
        mean_auroc = np.nanmean(list(fold_aurocs.values()))
    else:
        mean_auroc = np.nan

    return mean_auroc, fold_aurocs


def run_category_a_comparison():
    """Run gene vs pathway LOMO for all tissues."""
    print("=" * 70)
    print("J5 Category A: Spaceflight Detection (LOMO)")
    print("=" * 70)

    results = {}

    for tissue, missions in TISSUE_MISSIONS.items():
        print(f"\n  {tissue} ({len(missions)} missions):")

        gene_mean, gene_folds = lomo_auroc(tissue, "gene")
        pathway_mean, pathway_folds = lomo_auroc(tissue, "pathway")

        print(f"    Gene    AUROC: {gene_mean:.4f}" if gene_mean else "    Gene    AUROC: N/A")
        print(f"    Pathway AUROC: {pathway_mean:.4f}" if pathway_mean else "    Pathway AUROC: N/A")

        if gene_mean and pathway_mean:
            diff = pathway_mean - gene_mean
            print(f"    Diff (P-G):    {diff:+.4f}")

        results[tissue] = {
            "gene_mean_auroc": gene_mean,
            "gene_per_fold": gene_folds,
            "pathway_mean_auroc": pathway_mean,
            "pathway_per_fold": pathway_folds,
        }

    return results


# ── Collect Existing Results ─────────────────────────────────────────────────

def load_category_c_results():
    """Load Method A (gene) vs Method C (pathway) from C summary."""
    f = RESULTS_DIR / "C_cross_tissue_summary.json"
    if not f.exists():
        return {}
    with open(f) as fh:
        data = json.load(fh)

    results = {}
    for task_id, task_data in data.get("tasks", {}).items():
        gene_auroc = task_data.get("method_a", {}).get("auroc")
        pathway_auroc = task_data.get("method_c", {}).get("auroc")
        results[task_id] = {
            "train": task_data.get("train"),
            "test": task_data.get("test"),
            "gene_auroc": gene_auroc,
            "pathway_auroc": pathway_auroc,
        }
    return results


def load_category_d_results():
    """Load gene vs pathway F1 from D summary."""
    f = RESULTS_DIR / "D_condition_summary.json"
    if not f.exists():
        return {}
    with open(f) as fh:
        data = json.load(fh)

    return data.get("j5_comparison", {})


# ── Build J5 Summary Table ───────────────────────────────────────────────────

def build_j5_table(cat_a, cat_c, cat_d):
    """Build unified comparison table."""
    rows = []

    # Category A
    for tissue, data in cat_a.items():
        gene = data.get("gene_mean_auroc")
        pathway = data.get("pathway_mean_auroc")
        if gene is not None and pathway is not None:
            rows.append({
                "category": "A",
                "task": f"A_{tissue}",
                "description": f"Spaceflight detection ({tissue}, LOMO)",
                "metric": "AUROC",
                "gene": round(gene, 4),
                "pathway": round(pathway, 4),
                "diff": round(pathway - gene, 4),
                "winner": "pathway" if pathway > gene else "gene",
            })

    # Category C
    for task_id, data in cat_c.items():
        gene = data.get("gene_auroc")
        pathway = data.get("pathway_auroc")
        if gene is not None and pathway is not None:
            rows.append({
                "category": "C",
                "task": task_id,
                "description": f"Cross-tissue {data.get('train','?')}→{data.get('test','?')}",
                "metric": "AUROC",
                "gene": round(gene, 4),
                "pathway": round(pathway, 4),
                "diff": round(pathway - gene, 4),
                "winner": "pathway" if pathway > gene else "gene",
            })

    # Category D
    for task_id, data in cat_d.items():
        gene = data.get("gene_macro_f1")
        pathway = data.get("pathway_macro_f1")
        if gene is not None and pathway is not None:
            rows.append({
                "category": "D",
                "task": task_id,
                "description": task_id.replace("_", " "),
                "metric": "macro-F1",
                "gene": round(gene, 4),
                "pathway": round(pathway, 4),
                "diff": round(pathway - gene, 4),
                "winner": "pathway" if pathway > gene else "gene",
            })

    return rows


def print_j5_table(rows):
    """Print formatted table for paper."""
    print("\n" + "=" * 90)
    print("J5: Gene-level vs Pathway-level Feature Representation — Full Comparison")
    print("=" * 90)
    print(f"{'Category':<4} {'Task':<20} {'Metric':<8} {'Gene':>8} {'Pathway':>8} {'Diff':>8} {'Winner':>10}")
    print("-" * 90)

    for row in rows:
        print(f"{row['category']:<4} {row['task']:<20} {row['metric']:<8} "
              f"{row['gene']:>8.4f} {row['pathway']:>8.4f} {row['diff']:>+8.4f} "
              f"{row['winner']:>10}")

    # Summary statistics
    pathway_wins = sum(1 for r in rows if r["winner"] == "pathway")
    gene_wins = sum(1 for r in rows if r["winner"] == "gene")
    mean_diff = np.mean([r["diff"] for r in rows])

    print("-" * 90)
    print(f"Summary: Pathway wins {pathway_wins}/{len(rows)}, "
          f"Gene wins {gene_wins}/{len(rows)}")
    print(f"Mean diff (pathway - gene): {mean_diff:+.4f}")

    # Per-category summary
    for cat in ["A", "C", "D"]:
        cat_rows = [r for r in rows if r["category"] == cat]
        if cat_rows:
            pw = sum(1 for r in cat_rows if r["winner"] == "pathway")
            md = np.mean([r["diff"] for r in cat_rows])
            print(f"  Category {cat}: pathway wins {pw}/{len(cat_rows)}, "
                  f"mean diff = {md:+.4f}")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="J5: Gene vs Pathway comparison")
    parser.add_argument("--skip-category-a", action="store_true",
                        help="Skip Category A LOMO (use only existing C/D results)")
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Category A: fresh LOMO evaluation
    cat_a = {}
    if not args.skip_category_a:
        cat_a = run_category_a_comparison()

    # Category C: existing results
    cat_c = load_category_c_results()
    print(f"\nCategory C: Loaded {len(cat_c)} cross-tissue pairs")

    # Category D: existing results
    cat_d = load_category_d_results()
    print(f"Category D: Loaded {len(cat_d)} tasks")

    # Build table
    rows = build_j5_table(cat_a, cat_c, cat_d)
    print_j5_table(rows)

    # Save
    output = {
        "j5_table": rows,
        "category_a_details": cat_a,
        "category_c_source": "C_cross_tissue_summary.json",
        "category_d_source": "D_condition_summary.json",
        "summary": {
            "total_comparisons": len(rows),
            "pathway_wins": sum(1 for r in rows if r["winner"] == "pathway"),
            "gene_wins": sum(1 for r in rows if r["winner"] == "gene"),
            "mean_diff": float(np.mean([r["diff"] for r in rows])) if rows else 0,
        },
        "timestamp": datetime.now().isoformat(),
    }

    out_file = RESULTS_DIR / "J5_gene_vs_pathway.json"
    with open(out_file, "w") as f:
        json.dump(output, f, indent=2, default=str)

    print(f"\nSaved to {out_file}")


if __name__ == "__main__":
    main()
