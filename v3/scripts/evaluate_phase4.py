#!/usr/bin/env python3
"""
evaluate_phase4.py — Phase 4 spaceflight detection benchmark.

Tasks:
  A7: Lung (single-mission RR-6, stratified 5-fold CV)
  A7b: Colon (single-mission RR-6, stratified 5-fold CV)
  A8: Extended Skin (existing LOMO + new missions: RR-5 dorsal/femoral)

Methodology: v1 PCA-LR baseline (DD-01, DD-03 compliant)
  - Variance filter (top 75%, train-only)
  - StandardScaler (train-fit, test-transform)
  - PCA (n_components=20)
  - LogisticRegression

Metrics:
  - AUROC per fold
  - Bootstrap 95% CI (n=2000)
  - Permutation p-value (n=1000)

Usage:
    python3 v3/scripts/evaluate_phase4.py              # All tasks
    python3 v3/scripts/evaluate_phase4.py --task A7    # Single task

Output:
    v3/evaluation/A7_lung.json
    v3/evaluation/A7b_colon.json
    v3/evaluation/A8_extended_skin.json
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

BASE_DIR = Path(__file__).resolve().parent.parent.parent
PROCESSED_DIR = BASE_DIR / "processed" / "A_detection"
OUTPUT_DIR = BASE_DIR / "v3" / "evaluation"

VARIANCE_PERCENTILE = 0.25  # Keep top 75%
N_PCA = 20
N_BOOTSTRAP = 2000
N_PERM = 1000
SEED = 42


def load_single_mission(tissue, mission):
    """Load log2_norm + metadata for a single mission."""
    log2_path = PROCESSED_DIR / tissue / f"{tissue}_{mission}_log2_norm.csv"
    meta_path = PROCESSED_DIR / tissue / f"{tissue}_{mission}_metadata.csv"

    log2 = pd.read_csv(log2_path, index_col=0)
    meta = pd.read_csv(meta_path, index_col=0)

    # Ensure samples × genes orientation
    if log2.shape[0] < log2.shape[1]:
        pass  # Already samples × genes
    else:
        log2 = log2.T

    # Filter to ENSMUSG columns only
    ensmusg_cols = [c for c in log2.columns if str(c).startswith("ENSMUSG")]
    log2 = log2[ensmusg_cols]
    log2 = log2.apply(pd.to_numeric, errors="coerce").fillna(0)

    # Align
    common = sorted(set(log2.index) & set(meta.index))
    log2 = log2.loc[common]
    meta = meta.loc[common]

    return log2, meta


def make_binary_labels(meta):
    """Convert metadata labels to binary (Flight=1, Ground=0).
    Excludes Basal, Vivarium, AG, etc.
    Handles both v1 labels (GC) and Phase 4 labels (Ground).
    """
    FLIGHT_LABELS = {"Flight"}
    GROUND_LABELS = {"Ground", "GC", "Ground Control"}

    is_flight = meta["label"].isin(FLIGHT_LABELS)
    is_ground = meta["label"].isin(GROUND_LABELS)
    mask = is_flight | is_ground
    meta_bin = meta.loc[mask].copy()
    y = is_flight.loc[mask].astype(int)
    return meta_bin, y


def variance_filter(train_X, percentile=VARIANCE_PERCENTILE):
    """Train-only variance filter."""
    gene_var = train_X.var(axis=0)
    threshold = gene_var.quantile(percentile)
    selected = gene_var[gene_var >= threshold].index.tolist()
    return selected


def bootstrap_auroc(y_true, y_score, n_boot=N_BOOTSTRAP, seed=SEED):
    """Bootstrap 95% CI for AUROC."""
    rng = np.random.default_rng(seed)
    n = len(y_true)
    scores = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if len(np.unique(y_true[idx])) < 2:
            continue
        scores.append(roc_auc_score(y_true[idx], y_score[idx]))
    if not scores:
        return (float("nan"), float("nan"))
    return float(np.percentile(scores, 2.5)), float(np.percentile(scores, 97.5))


def permutation_p(y_true, y_score, n_perm=N_PERM, seed=SEED):
    """Permutation test p-value."""
    rng = np.random.default_rng(seed)
    observed = roc_auc_score(y_true, y_score)
    null_scores = []
    for _ in range(n_perm):
        y_perm = rng.permutation(y_true)
        if len(np.unique(y_perm)) < 2:
            continue
        null_scores.append(roc_auc_score(y_perm, y_score))
    null_arr = np.array(null_scores)
    p = float((np.sum(null_arr >= observed) + 1) / (len(null_arr) + 1))
    return p


def evaluate_single_mission_cv(tissue, mission, n_splits=5):
    """Evaluate spaceflight detection with stratified k-fold CV."""
    log2, meta = load_single_mission(tissue, mission)
    meta_bin, y = make_binary_labels(meta)
    X = log2.loc[meta_bin.index]

    print(f"  Data: {X.shape[0]} samples × {X.shape[1]} genes")
    print(f"  Labels: Flight={int(y.sum())}, Ground={int(len(y)-y.sum())}")

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=SEED)
    fold_results = []

    for fold_idx, (train_idx, test_idx) in enumerate(skf.split(X, y)):
        train_X = X.iloc[train_idx]
        test_X = X.iloc[test_idx]
        train_y = y.iloc[train_idx].values
        test_y = y.iloc[test_idx].values

        # Variance filter (train-only)
        selected = variance_filter(train_X)
        train_X_f = train_X[selected]
        test_X_f = test_X[selected]

        # Standardize
        scaler = StandardScaler()
        train_X_s = scaler.fit_transform(train_X_f)
        test_X_s = scaler.transform(test_X_f)

        # PCA
        n_comp = min(N_PCA, train_X_s.shape[0] - 1, train_X_s.shape[1])
        pca = PCA(n_components=n_comp, random_state=SEED)
        train_X_p = pca.fit_transform(train_X_s)
        test_X_p = pca.transform(test_X_s)

        # Logistic Regression
        lr = LogisticRegression(max_iter=1000, random_state=SEED)
        lr.fit(train_X_p, train_y)
        y_score = lr.predict_proba(test_X_p)[:, 1]

        # Metrics
        if len(np.unique(test_y)) < 2:
            print(f"    Fold {fold_idx}: SKIP (single class in test)")
            continue

        auroc = roc_auc_score(test_y, y_score)
        ci_lo, ci_hi = bootstrap_auroc(test_y, y_score)
        perm_p = permutation_p(test_y, y_score)

        fold_results.append({
            "fold": fold_idx,
            "n_train": len(train_y),
            "n_test": len(test_y),
            "n_genes_selected": len(selected),
            "n_pca": n_comp,
            "auroc": round(auroc, 4),
            "ci_lower": round(ci_lo, 4),
            "ci_upper": round(ci_hi, 4),
            "perm_p": round(perm_p, 4),
        })
        print(f"    Fold {fold_idx}: AUROC={auroc:.3f} [{ci_lo:.3f}, {ci_hi:.3f}] p={perm_p:.4f}")

    if not fold_results:
        return None

    mean_auroc = np.mean([r["auroc"] for r in fold_results])
    mean_ci_lo = np.mean([r["ci_lower"] for r in fold_results])
    mean_perm_p = np.mean([r["perm_p"] for r in fold_results])

    return {
        "tissue": tissue,
        "mission": mission,
        "method": "PCA-LR (5-fold stratified CV)",
        "n_folds": len(fold_results),
        "mean_auroc": round(mean_auroc, 4),
        "mean_ci_lower": round(mean_ci_lo, 4),
        "mean_perm_p": round(mean_perm_p, 4),
        "fold_results": fold_results,
    }


def evaluate_extended_skin():
    """Evaluate skin with new missions added to LOMO.

    Existing skin missions (v1): MHU-2_(dorsal), MHU-2_(femoral), RR-6, RR-7
    New missions (v3): RR-5_dorsal (BAL-TAL), RR-5_femoral (BAL-TAL)

    Note: RR-6_dorsal (Phase 4, OSD-243) is excluded because v1 RR-6 is
    already from the same OSD-243 dataset. skin_all_missions is also excluded.
    """
    print("\n=== A8: Extended Skin ===")

    # Missions to exclude from glob results
    EXCLUDE_MISSIONS = {
        "all_missions",   # v1 combined file — would duplicate all data
        "RR-6_dorsal",    # Phase 4 OSD-243 — duplicate of v1 RR-6 (also OSD-243)
    }

    # Load all available skin missions
    skin_dir = PROCESSED_DIR / "skin"
    missions = []
    for f in sorted(skin_dir.glob("skin_*_log2_norm.csv")):
        mission = f.stem.replace("skin_", "").replace("_log2_norm", "")
        if mission in EXCLUDE_MISSIONS:
            print(f"  Excluding: {mission}")
            continue
        missions.append(mission)

    print(f"  Missions: {missions}")

    # Load and concatenate all missions
    all_X = []
    all_meta = []
    loaded_missions = []
    for mission in missions:
        log2, meta = load_single_mission("skin", mission)
        meta_bin, y = make_binary_labels(meta)
        if len(meta_bin) == 0:
            print(f"  {mission}: SKIP (no Flight/Ground samples)")
            continue
        X = log2.loc[meta_bin.index]
        meta_bin = meta_bin.copy()
        meta_bin["_y"] = y
        # Override mission column with filename-extracted name
        # (v1 metadata uses "MHU-2 (dorsal)" but filename gives "MHU-2_(dorsal)")
        meta_bin["mission"] = mission
        all_X.append(X)
        all_meta.append(meta_bin)
        loaded_missions.append(mission)
    missions = loaded_missions

    if not all_X:
        print("  No skin data found!")
        return None

    # Find common genes
    common_genes = set(all_X[0].columns)
    for X in all_X[1:]:
        common_genes &= set(X.columns)
    common_genes = sorted(common_genes)
    print(f"  Common genes across {len(missions)} missions: {len(common_genes)}")

    # Concatenate
    combined_X = pd.concat([X[common_genes] for X in all_X])
    combined_meta = pd.concat(all_meta)
    combined_y = combined_meta["_y"]

    print(f"  Combined: {combined_X.shape[0]} samples × {combined_X.shape[1]} genes")
    print(f"  Labels: Flight={int(combined_y.sum())}, Ground={int(len(combined_y)-combined_y.sum())}")

    # LOMO: leave one mission out
    fold_results = []
    for test_mission in missions:
        test_mask = combined_meta["mission"] == test_mission
        if test_mask.sum() < 3:
            print(f"    {test_mission}: SKIP (n={test_mask.sum()} < 3)")
            continue

        train_X = combined_X.loc[~test_mask]
        test_X = combined_X.loc[test_mask]
        train_y = combined_y.loc[~test_mask].values
        test_y = combined_y.loc[test_mask].values

        if len(np.unique(test_y)) < 2:
            print(f"    {test_mission}: SKIP (single class)")
            continue

        # Pipeline
        selected = variance_filter(train_X)
        train_X_f = train_X[selected]
        test_X_f = test_X[selected]

        scaler = StandardScaler()
        train_X_s = scaler.fit_transform(train_X_f)
        test_X_s = scaler.transform(test_X_f)

        n_comp = min(N_PCA, train_X_s.shape[0] - 1, train_X_s.shape[1])
        pca = PCA(n_components=n_comp, random_state=SEED)
        train_X_p = pca.fit_transform(train_X_s)
        test_X_p = pca.transform(test_X_s)

        lr = LogisticRegression(max_iter=1000, random_state=SEED)
        lr.fit(train_X_p, train_y)
        y_score = lr.predict_proba(test_X_p)[:, 1]

        auroc = roc_auc_score(test_y, y_score)
        ci_lo, ci_hi = bootstrap_auroc(test_y, y_score)
        perm_p = permutation_p(test_y, y_score)

        strain = combined_meta.loc[test_mask, "strain"].iloc[0] if "strain" in combined_meta.columns else "?"
        fold_results.append({
            "test_mission": test_mission,
            "strain": strain,
            "n_train": len(train_y),
            "n_test": len(test_y),
            "n_genes_selected": len(selected),
            "auroc": round(auroc, 4),
            "ci_lower": round(ci_lo, 4),
            "ci_upper": round(ci_hi, 4),
            "perm_p": round(perm_p, 4),
        })
        print(f"    {test_mission} ({strain}): AUROC={auroc:.3f} [{ci_lo:.3f}, {ci_hi:.3f}] p={perm_p:.4f}")

    if not fold_results:
        return None

    mean_auroc = np.mean([r["auroc"] for r in fold_results])
    return {
        "tissue": "skin",
        "method": "PCA-LR (LOMO, extended)",
        "n_missions": len(missions),
        "n_folds": len(fold_results),
        "mean_auroc": round(mean_auroc, 4),
        "fold_results": fold_results,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--task", choices=["A7", "A7b", "A8", "all"], default="all")
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    results = {}

    if args.task in ("A7", "all"):
        print("\n=== A7: Lung Spaceflight Detection ===")
        r = evaluate_single_mission_cv("lung", "RR-6")
        if r:
            results["A7"] = r
            with open(OUTPUT_DIR / "A7_lung.json", "w") as f:
                json.dump(r, f, indent=2)
            print(f"  → Mean AUROC: {r['mean_auroc']}")

    if args.task in ("A7b", "all"):
        print("\n=== A7b: Colon Spaceflight Detection ===")
        r = evaluate_single_mission_cv("colon", "RR-6")
        if r:
            results["A7b"] = r
            with open(OUTPUT_DIR / "A7b_colon.json", "w") as f:
                json.dump(r, f, indent=2)
            print(f"  → Mean AUROC: {r['mean_auroc']}")

    if args.task in ("A8", "all"):
        r = evaluate_extended_skin()
        if r:
            results["A8"] = r
            with open(OUTPUT_DIR / "A8_extended_skin.json", "w") as f:
                json.dump(r, f, indent=2)
            print(f"  → Mean AUROC: {r['mean_auroc']}")

    # Summary
    print("\n" + "=" * 70)
    print("PHASE 4 EVALUATION SUMMARY")
    print("=" * 70)
    for task_id, r in results.items():
        print(f"  {task_id}: mean AUROC = {r['mean_auroc']:.3f} ({r['method']})")


if __name__ == "__main__":
    main()
