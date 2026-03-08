#!/usr/bin/env python3
"""
condition_prediction.py — GeneLab_benchmark: Category D Condition Prediction

D3: Mission ID prediction (liver, 6-class)
  - Tests how distinguishable missions are in transcriptomic space
  - High accuracy ≈ strong mission batch effects (confounder quantification)
  - Evaluation: RepeatedStratifiedKFold (5-fold × 10 repeats)

D4: Strain effect (thymus GC only, C57BL/6J vs C57BL/6CR)
  - Exploratory: minority class n=3
  - Evaluation: Leave-One-Out CV

D5: Hardware effect (RR vs MHU, derived from mission)
  - Collinear with D3 — interpret as upper bound
  - Implemented for liver and thymus

D6: Artificial gravity effect (MHU-2, 3-class: uG vs AG vs GC)
  - Separates microgravity from other spaceflight stressors
  - Evaluation: Leave-One-Out CV (tiny sample size, n=9)
  - Implemented for both liver and thymus MHU-2

Each task evaluated with both feature representations:
  - Gene-level: log2 normalized counts (variance-filtered within train fold)
  - Pathway-level: GSVA Hallmark scores (50 pathways, no filtering needed)

Metrics: macro-F1 (primary), accuracy, per-class F1, confusion matrix

Output:
  processed/D_condition/
    D3_liver_results.json
    D4_strain_results.json
    D5_liver_results.json
    D5_thymus_results.json
    D6_liver_results.json
    D6_thymus_results.json
  evaluation/
    D_condition_summary.json  (merged — single-task runs preserve prior results)

Usage:
  python scripts/condition_prediction.py --task D3
  python scripts/condition_prediction.py --task D6
  python scripts/condition_prediction.py --all
  python scripts/condition_prediction.py --all --no-bootstrap
"""

import json
import argparse
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from collections import Counter

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")

from utils import (
    load_metadata, load_gene_features, load_pathway_features,
    align_features_with_meta, TISSUE_MISSIONS as ALL_TISSUE_MISSIONS,
)

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent.parent
PROCESSED_DIR = BASE_DIR / "processed" / "A_detection"
PATHWAY_DIR = BASE_DIR / "processed" / "pathway_scores"
D_OUTPUT_DIR = BASE_DIR / "processed" / "D_condition"
RESULTS_DIR = BASE_DIR / "evaluation"

# ── Config ─────────────────────────────────────────────────────────────────────
VARIANCE_PERCENTILE = 0.25   # keep top 75% by variance (train-only)
N_BOOTSTRAP = 2000
CI_ALPHA = 0.05
N_SPLITS = 5
N_REPEATS = 10

# D tasks only use liver and thymus
TISSUE_MISSIONS = {
    "liver": ALL_TISSUE_MISSIONS["liver"],
    "thymus": ALL_TISSUE_MISSIONS["thymus"],
}

# ── Task Definitions ──────────────────────────────────────────────────────────
TASKS = {
    "D3": {
        "tissue": "liver",
        "target": "mission",
        "description": "Mission ID prediction (6-class, liver)",
        "cv_method": "repeated_stratified_kfold",
    },
    "D4": {
        "tissue": "thymus",
        "target": "strain",
        "description": "Strain effect (thymus GC, C57BL/6J vs C57BL/6CR, exploratory n=3)",
        "cv_method": "loo",
    },
    "D5_liver": {
        "tissue": "liver",
        "target": "hardware",
        "description": "Hardware effect (liver, RR vs MHU, collinear with D3)",
        "cv_method": "loo",
    },
    "D5_thymus": {
        "tissue": "thymus",
        "target": "hardware",
        "description": "Hardware effect (thymus, RR vs MHU, collinear with D3)",
        "cv_method": "repeated_stratified_kfold",
    },
    "D6_liver": {
        "tissue": "liver",
        "target": "gravity",
        "description": "Artificial gravity (liver MHU-2, 3-class)",
        "cv_method": "loo",
        "missions_filter": ["MHU-2"],
    },
    "D6_thymus": {
        "tissue": "thymus",
        "target": "gravity",
        "description": "Artificial gravity (thymus MHU-2, 3-class)",
        "cv_method": "loo",
        "missions_filter": ["MHU-2"],
    },
}

# Hardware mapping (derived from mission name)
HARDWARE_MAP = {
    "RR-1": "RR", "RR-3": "RR", "RR-5": "RR",
    "RR-6": "RR", "RR-7": "RR", "RR-8": "RR", "RR-9": "RR",
    "MHU-1": "MHU", "MHU-2": "MHU",
    "TBD": "unknown",
}


# ── Label Extraction ──────────────────────────────────────────────────────────

def get_mission_labels(meta):
    """Extract mission ID labels for D3."""
    return meta["mission"].values, sorted(meta["mission"].unique())


def get_strain_labels(meta):
    """Extract strain labels for D4. Binary: C57BL/6J vs C57BL/6CR."""
    strains = meta["strain"].values
    classes = sorted(set(strains))
    return strains, classes


def get_hardware_labels(meta):
    """Extract hardware labels for D5 (derived from mission name)."""
    hw = meta["mission"].map(HARDWARE_MAP).values
    # Remove unknowns
    valid = hw != "unknown"
    if not valid.all():
        print(f"  [INFO] Removed {(~valid).sum()} samples with unknown hardware")
    classes = sorted(set(hw[valid]))
    return hw, classes, valid


def get_gravity_labels(meta):
    """Extract gravity condition labels for D6.
    Maps: Flight/uG → 'uG', AG → 'AG', GC → 'GC'
    """
    label_map = {"Flight": "uG", "AG": "AG", "GC": "GC"}
    raw = meta["label"].values
    mapped = np.array([label_map.get(v, v) for v in raw])
    # Verify all mapped
    valid = np.isin(mapped, ["uG", "AG", "GC"])
    if not valid.all():
        bad = set(mapped[~valid])
        raise ValueError(f"Unexpected labels for D6: {bad}")
    classes = ["uG", "AG", "GC"]
    return mapped, classes


# ── Statistical Utilities ────────────────────────────────────────────────────

def bootstrap_macro_f1_ci(y_true, y_pred, n_boot=None, alpha=CI_ALPHA,
                          random_state=42):
    """Bootstrap 95% CI for macro-F1."""
    from sklearn.metrics import f1_score

    if n_boot is None:
        n_boot = N_BOOTSTRAP
    if n_boot <= 0:
        return np.nan, np.nan

    rng = np.random.RandomState(random_state)
    n = len(y_true)
    boot_f1s = []
    for _ in range(n_boot):
        idx = rng.randint(0, n, size=n)
        y_t = y_true[idx]
        y_p = y_pred[idx]
        if len(np.unique(y_t)) < 2:
            continue
        try:
            boot_f1s.append(f1_score(y_t, y_p, average="macro", zero_division=0))
        except Exception:
            continue

    if len(boot_f1s) < 100:
        return np.nan, np.nan

    lo = np.percentile(boot_f1s, 100 * alpha / 2)
    hi = np.percentile(boot_f1s, 100 * (1 - alpha / 2))
    return lo, hi


def permutation_test_macro_f1(y_true, y_pred, n_perm=10000, random_state=42):
    """Permutation test: is macro-F1 significantly better than chance?"""
    from sklearn.metrics import f1_score

    observed = f1_score(y_true, y_pred, average="macro", zero_division=0)
    rng = np.random.RandomState(random_state)
    count_ge = 0
    for _ in range(n_perm):
        perm = rng.permutation(y_true)
        pf1 = f1_score(perm, y_pred, average="macro", zero_division=0)
        if pf1 >= observed:
            count_ge += 1
    return (count_ge + 1) / (n_perm + 1)


# ── Core Classification ─────────────────────────────────────────────────────

def classify_repeated_stratified_kfold(X, y, feature_mode="gene",
                                       n_splits=N_SPLITS, n_repeats=N_REPEATS):
    """Repeated Stratified K-Fold classification for D3.
    Returns per-sample predictions (aggregated across repeats).
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    from sklearn.model_selection import RepeatedStratifiedKFold

    rskf = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats,
                                   random_state=42)
    n_samples = len(y)
    classes = sorted(np.unique(y))
    n_classes = len(classes)

    # Accumulate prediction probabilities across repeats
    prob_accum = np.zeros((n_samples, n_classes))
    count_accum = np.zeros(n_samples)

    for train_idx, test_idx in rskf.split(X, y):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train = y[train_idx]

        # Feature processing
        if feature_mode == "gene":
            # Variance filter on train
            var = np.var(X_train, axis=0)
            threshold = np.percentile(var, VARIANCE_PERCENTILE * 100)
            keep = var > threshold
            if keep.sum() < 50:
                keep = np.ones(X_train.shape[1], dtype=bool)
            X_train = X_train[:, keep]
            X_test = X_test[:, keep]

            # PCA for dimensionality reduction
            n_components = min(50, X_train.shape[0] - 1, X_train.shape[1])
            scaler = StandardScaler()
            X_train = scaler.fit_transform(X_train)
            X_test = scaler.transform(X_test)
            pca = PCA(n_components=n_components, random_state=42)
            X_train = pca.fit_transform(X_train)
            X_test = pca.transform(X_test)
        else:
            # Pathway features: just scale, no PCA needed (only ~50 features)
            scaler = StandardScaler()
            X_train = scaler.fit_transform(X_train)
            X_test = scaler.transform(X_test)

        # Logistic Regression (multinomial)
        clf = LogisticRegression(
            solver="lbfgs", class_weight="balanced",
            max_iter=2000, C=1.0, random_state=42,
        )
        clf.fit(X_train, y_train)
        proba = clf.predict_proba(X_test)

        # Map to correct class indices
        for i, cls in enumerate(clf.classes_):
            cls_idx = classes.index(cls)
            prob_accum[test_idx, cls_idx] += proba[:, i]
        count_accum[test_idx] += 1

    # Average probabilities and predict
    prob_avg = prob_accum / count_accum[:, None]
    y_pred = np.array([classes[i] for i in np.argmax(prob_avg, axis=1)])

    return y_pred, prob_avg, classes


def classify_loo(X, y, feature_mode="gene"):
    """Leave-One-Out CV for small datasets (D6).
    Returns per-sample predictions.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import LeaveOneOut

    loo = LeaveOneOut()
    classes = sorted(np.unique(y))
    n_samples = len(y)
    y_pred = np.empty(n_samples, dtype=y.dtype)
    y_proba = np.zeros((n_samples, len(classes)))

    for train_idx, test_idx in loo.split(X, y):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train = y[train_idx]

        # Scale
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)

        if feature_mode == "gene":
            # Variance filter on train
            var = np.var(X_train, axis=0)
            threshold = np.percentile(var, VARIANCE_PERCENTILE * 100)
            keep = var > threshold
            if keep.sum() < 20:
                keep = np.ones(X_train.shape[1], dtype=bool)
            X_train_f = X_train[:, keep]
            X_test_f = X_test[:, keep]

            # PCA
            from sklearn.decomposition import PCA
            n_comp = min(5, X_train_f.shape[0] - 1, X_train_f.shape[1])
            pca = PCA(n_components=n_comp, random_state=42)
            X_train_f = pca.fit_transform(X_train_f)
            X_test_f = pca.transform(X_test_f)
        else:
            X_train_f = X_train
            X_test_f = X_test

        clf = LogisticRegression(
            solver="lbfgs", class_weight="balanced",
            max_iter=2000, C=1.0, random_state=42,
        )
        clf.fit(X_train_f, y_train)
        y_pred[test_idx] = clf.predict(X_test_f)

        proba = clf.predict_proba(X_test_f)
        for i, cls in enumerate(clf.classes_):
            cls_idx = classes.index(cls)
            y_proba[test_idx, cls_idx] = proba[:, i]

    return y_pred, y_proba, classes


# ── Task Runners ─────────────────────────────────────────────────────────────

def run_d3(do_bootstrap=True):
    """Run D3: Mission ID prediction on liver."""
    from sklearn.metrics import (f1_score, accuracy_score,
                                 classification_report, confusion_matrix)

    print("=" * 70)
    print("D3: Mission ID Prediction (liver, 6-class)")
    print("=" * 70)

    tissue = "liver"
    meta = load_metadata(tissue)

    # Load features
    gene_feat = load_gene_features(tissue)
    pathway_feat = load_pathway_features(tissue)

    # Class distribution
    mission_counts = Counter(meta["mission"])
    print(f"\nSamples: {len(meta)}")
    print("Mission distribution:")
    for m in sorted(mission_counts):
        print(f"  {m}: {mission_counts[m]} samples ({100*mission_counts[m]/len(meta):.1f}%)")

    results = {"task": "D3", "tissue": tissue, "n_samples": len(meta),
               "class_distribution": dict(mission_counts), "feature_modes": {}}

    for mode, feat_df in [("gene", gene_feat), ("pathway", pathway_feat)]:
        if feat_df is None:
            print(f"\n  [{mode}] No features available, skipping")
            continue

        print(f"\n  [{mode}] Features: {feat_df.shape[1]}")

        # Align
        feat_aligned, meta_aligned = align_features_with_meta(feat_df, meta)
        print(f"  [{mode}] Aligned samples: {len(feat_aligned)}")

        X = feat_aligned.values.astype(float)
        y, classes = get_mission_labels(meta_aligned)

        # Check NaN
        nan_mask = np.isnan(X)
        if nan_mask.any():
            nan_cols = nan_mask.any(axis=0).sum()
            print(f"  [{mode}] Warning: {nan_cols} features with NaN, imputing with 0")
            X = np.nan_to_num(X, nan=0.0)

        # Classify
        print(f"  [{mode}] Running {N_SPLITS}-fold × {N_REPEATS}-repeat stratified CV...")
        y_pred, y_proba, classes = classify_repeated_stratified_kfold(
            X, y, feature_mode=mode, n_splits=N_SPLITS, n_repeats=N_REPEATS
        )

        macro_f1 = f1_score(y, y_pred, average="macro", zero_division=0)
        accuracy = accuracy_score(y, y_pred)
        per_class = classification_report(y, y_pred, output_dict=True,
                                          zero_division=0)
        cm = confusion_matrix(y, y_pred, labels=classes).tolist()

        print(f"  [{mode}] macro-F1: {macro_f1:.4f}")
        print(f"  [{mode}] accuracy: {accuracy:.4f}")

        # Per-class F1
        for cls in classes:
            if cls in per_class:
                f1 = per_class[cls]["f1-score"]
                print(f"    {cls}: F1={f1:.3f}")

        # Bootstrap CI
        ci_low, ci_high = np.nan, np.nan
        if do_bootstrap:
            ci_low, ci_high = bootstrap_macro_f1_ci(y, y_pred)
            print(f"  [{mode}] 95% CI: [{ci_low:.4f}, {ci_high:.4f}]")

        # Permutation test
        perm_p = permutation_test_macro_f1(y, y_pred)
        print(f"  [{mode}] permutation p-value: {perm_p:.6f}")

        mode_result = {
            "n_features": feat_aligned.shape[1],
            "macro_f1": macro_f1,
            "accuracy": accuracy,
            "ci_low": ci_low,
            "ci_high": ci_high,
            "perm_p": perm_p,
            "per_class_f1": {cls: per_class.get(cls, {}).get("f1-score", 0)
                            for cls in classes},
            "confusion_matrix": cm,
            "classes": classes,
        }
        results["feature_modes"][mode] = mode_result

    return results


def run_d6(tissue, do_bootstrap=True):
    """Run D6: Artificial gravity prediction for a specific tissue's MHU-2."""
    from sklearn.metrics import (f1_score, accuracy_score,
                                 classification_report, confusion_matrix)

    task_name = f"D6_{tissue}"
    print("=" * 70)
    print(f"D6: Artificial Gravity ({tissue} MHU-2, 3-class)")
    print("=" * 70)

    meta = load_metadata(tissue)

    # Filter to MHU-2 only
    meta_mhu2 = meta[meta["mission"] == "MHU-2"].copy()
    if len(meta_mhu2) == 0:
        print(f"  No MHU-2 samples found for {tissue}")
        return None

    # Filter to uG/AG/GC labels only
    valid_labels = {"Flight", "AG", "GC"}
    meta_mhu2 = meta_mhu2[meta_mhu2["label"].isin(valid_labels)]

    print(f"\nMHU-2 samples: {len(meta_mhu2)}")
    label_counts = Counter(meta_mhu2["label"])
    for lab in sorted(label_counts):
        print(f"  {lab}: {label_counts[lab]}")

    results = {"task": task_name, "tissue": tissue, "n_samples": len(meta_mhu2),
               "label_distribution": dict(label_counts), "feature_modes": {}}

    # Gene features
    gene_feat = load_gene_features(tissue)
    # Pathway features (MHU-2 GSVA)
    pathway_feat = load_pathway_features(tissue)

    for mode, feat_df in [("gene", gene_feat), ("pathway", pathway_feat)]:
        if feat_df is None:
            print(f"\n  [{mode}] No features available, skipping")
            continue

        print(f"\n  [{mode}] Total features available: {feat_df.shape[1]}")

        # Align with MHU-2 metadata
        try:
            feat_aligned, meta_aligned = align_features_with_meta(feat_df, meta_mhu2)
        except ValueError as e:
            print(f"  [{mode}] Alignment failed: {e}")
            continue

        print(f"  [{mode}] Aligned MHU-2 samples: {len(feat_aligned)}")

        X = feat_aligned.values.astype(float)
        y, classes = get_gravity_labels(meta_aligned)

        # Check NaN
        nan_mask = np.isnan(X)
        if nan_mask.any():
            X = np.nan_to_num(X, nan=0.0)

        if len(np.unique(y)) < 2:
            print(f"  [{mode}] Too few classes after filtering, skipping")
            continue

        # LOO-CV
        print(f"  [{mode}] Running LOO-CV (n={len(y)})...")
        y_pred, y_proba, classes = classify_loo(X, y, feature_mode=mode)

        macro_f1 = f1_score(y, y_pred, average="macro", zero_division=0)
        accuracy = accuracy_score(y, y_pred)
        per_class = classification_report(y, y_pred, output_dict=True,
                                          zero_division=0)
        cm = confusion_matrix(y, y_pred, labels=classes).tolist()

        print(f"  [{mode}] macro-F1: {macro_f1:.4f}")
        print(f"  [{mode}] accuracy: {accuracy:.4f}")
        for cls in classes:
            if cls in per_class:
                f1 = per_class[cls]["f1-score"]
                print(f"    {cls}: F1={f1:.3f}")

        # Bootstrap CI (may be unstable with n=9)
        ci_low, ci_high = np.nan, np.nan
        if do_bootstrap and len(y) >= 6:
            ci_low, ci_high = bootstrap_macro_f1_ci(y, y_pred)
            if not np.isnan(ci_low):
                print(f"  [{mode}] 95% CI: [{ci_low:.4f}, {ci_high:.4f}]")

        # Permutation test
        perm_p = permutation_test_macro_f1(y, y_pred, n_perm=5000)
        print(f"  [{mode}] permutation p-value: {perm_p:.6f}")

        mode_result = {
            "n_features": feat_aligned.shape[1],
            "macro_f1": macro_f1,
            "accuracy": accuracy,
            "ci_low": ci_low,
            "ci_high": ci_high,
            "perm_p": perm_p,
            "per_class_f1": {cls: per_class.get(cls, {}).get("f1-score", 0)
                            for cls in classes},
            "confusion_matrix": cm,
            "classes": classes,
            "y_true": y.tolist(),
            "y_pred": y_pred.tolist(),
        }
        results["feature_modes"][mode] = mode_result

    return results


# ── D4: Strain Effect ────────────────────────────────────────────────────────

def run_d4(do_bootstrap=True):
    """Run D4: Strain effect prediction on thymus GC samples only.

    EXPLORATORY ANALYSIS: C57BL/6CR has only n=3 samples (MHU-1 GC).
    Statistical power is extremely limited. Results should NOT be used
    for significance claims.
    """
    from sklearn.metrics import (f1_score, accuracy_score,
                                 classification_report, confusion_matrix)

    print("=" * 70)
    print("D4: Strain Effect (thymus GC, EXPLORATORY n_minority=3)")
    print("=" * 70)

    tissue = "thymus"
    meta = load_metadata(tissue)

    # Filter to GC samples only (to control for experimental condition)
    gc_labels = {"GC", "Ground Control", "Ground"}
    meta_gc = meta[meta["label"].isin(gc_labels)].copy()

    if "strain" not in meta_gc.columns:
        print("  [ERROR] No 'strain' column in metadata")
        return None

    strain_counts = Counter(meta_gc["strain"])
    print(f"\nGC samples by strain: {len(meta_gc)}")
    for s in sorted(strain_counts):
        print(f"  {s}: {strain_counts[s]}")

    # Warn about minority class
    min_count = min(strain_counts.values())
    if min_count <= 5:
        print(f"\n  ⚠ WARNING: minority class has only {min_count} samples.")
        print(f"    Results are EXPLORATORY ONLY. No significance claims possible.")

    results = {
        "task": "D4", "tissue": tissue, "n_samples": len(meta_gc),
        "sample_filter": "GC only (Ground Control)",
        "class_distribution": dict(strain_counts),
        "warning": f"Exploratory analysis: minority class n={min_count}. "
                   f"Statistical power extremely limited.",
        "feature_modes": {},
    }

    gene_feat = load_gene_features(tissue)
    pathway_feat = load_pathway_features(tissue)

    for mode, feat_df in [("gene", gene_feat), ("pathway", pathway_feat)]:
        if feat_df is None:
            print(f"\n  [{mode}] No features available, skipping")
            continue

        print(f"\n  [{mode}] Features: {feat_df.shape[1]}")

        try:
            feat_aligned, meta_aligned = align_features_with_meta(feat_df, meta_gc)
        except ValueError as e:
            print(f"  [{mode}] Alignment failed: {e}")
            continue

        print(f"  [{mode}] Aligned samples: {len(feat_aligned)}")

        X = feat_aligned.values.astype(float)
        y, classes = get_strain_labels(meta_aligned)
        X = np.nan_to_num(X, nan=0.0)

        if len(np.unique(y)) < 2:
            print(f"  [{mode}] Only one strain class present, skipping")
            continue

        # LOO-CV
        print(f"  [{mode}] Running LOO-CV (n={len(y)})...")
        y_pred, y_proba, classes = classify_loo(X, y, feature_mode=mode)

        macro_f1 = f1_score(y, y_pred, average="macro", zero_division=0)
        accuracy = accuracy_score(y, y_pred)
        per_class = classification_report(y, y_pred, output_dict=True,
                                          zero_division=0)
        cm = confusion_matrix(y, y_pred, labels=classes).tolist()

        print(f"  [{mode}] macro-F1: {macro_f1:.4f}")
        print(f"  [{mode}] accuracy: {accuracy:.4f}")
        for cls in classes:
            if cls in per_class:
                print(f"    {cls}: F1={per_class[cls]['f1-score']:.3f}")

        ci_low, ci_high = np.nan, np.nan
        if do_bootstrap and len(y) >= 6:
            ci_low, ci_high = bootstrap_macro_f1_ci(y, y_pred)
            if not np.isnan(ci_low):
                print(f"  [{mode}] 95% CI: [{ci_low:.4f}, {ci_high:.4f}]")

        perm_p = permutation_test_macro_f1(y, y_pred, n_perm=5000)
        print(f"  [{mode}] permutation p-value: {perm_p:.6f}")

        results["feature_modes"][mode] = {
            "n_features": feat_aligned.shape[1],
            "macro_f1": macro_f1,
            "accuracy": accuracy,
            "ci_low": ci_low,
            "ci_high": ci_high,
            "perm_p": perm_p,
            "per_class_f1": {cls: per_class.get(cls, {}).get("f1-score", 0)
                            for cls in classes},
            "confusion_matrix": cm,
            "classes": classes,
            "y_true": y.tolist(),
            "y_pred": y_pred.tolist(),
        }

    return results


# ── D5: Hardware Effect ──────────────────────────────────────────────────────

def run_d5(tissue, do_bootstrap=True):
    """Run D5: Hardware effect prediction (RR vs MHU, derived from mission).

    WARNING: Hardware is collinear with mission ID (D3). D5 AUROC should
    be interpreted as an upper bound of D3, not independently.
    """
    from sklearn.metrics import (f1_score, accuracy_score,
                                 classification_report, confusion_matrix)

    task_name = f"D5_{tissue}"
    print("=" * 70)
    print(f"D5: Hardware Effect ({tissue}, RR vs MHU)")
    print("=" * 70)

    meta = load_metadata(tissue)

    # Add hardware column
    meta = meta.copy()
    meta["hardware"] = meta["mission"].map(HARDWARE_MAP)
    meta = meta[meta["hardware"] != "unknown"]

    hw_counts = Counter(meta["hardware"])
    print(f"\nSamples by hardware: {len(meta)}")
    for hw in sorted(hw_counts):
        print(f"  {hw}: {hw_counts[hw]}")

    min_count = min(hw_counts.values())
    collinear_note = ("Hardware is collinear with mission ID (D3). "
                      "D5 macro-F1 should be interpreted as an upper bound of D3.")

    results = {
        "task": task_name, "tissue": tissue, "n_samples": len(meta),
        "class_distribution": dict(hw_counts),
        "collinearity_warning": collinear_note,
        "feature_modes": {},
    }

    gene_feat = load_gene_features(tissue)
    pathway_feat = load_pathway_features(tissue)

    for mode, feat_df in [("gene", gene_feat), ("pathway", pathway_feat)]:
        if feat_df is None:
            continue

        print(f"\n  [{mode}] Features: {feat_df.shape[1]}")

        try:
            feat_aligned, meta_aligned = align_features_with_meta(feat_df, meta)
        except ValueError as e:
            print(f"  [{mode}] Alignment failed: {e}")
            continue

        print(f"  [{mode}] Aligned samples: {len(feat_aligned)}")

        X = feat_aligned.values.astype(float)
        y = meta_aligned["hardware"].values
        X = np.nan_to_num(X, nan=0.0)

        classes = sorted(np.unique(y))
        if len(classes) < 2:
            print(f"  [{mode}] Only one hardware class present, skipping")
            continue

        # Choose CV method based on minority class size
        if min_count <= 15:
            print(f"  [{mode}] Running LOO-CV (n={len(y)}, minority={min_count})...")
            y_pred, y_proba, classes = classify_loo(X, y, feature_mode=mode)
            n_perm = 5000
        else:
            print(f"  [{mode}] Running {N_SPLITS}-fold × {N_REPEATS}-repeat CV...")
            y_pred, y_proba, classes = classify_repeated_stratified_kfold(
                X, y, feature_mode=mode)
            n_perm = 10000

        macro_f1 = f1_score(y, y_pred, average="macro", zero_division=0)
        accuracy = accuracy_score(y, y_pred)
        per_class = classification_report(y, y_pred, output_dict=True,
                                          zero_division=0)
        cm = confusion_matrix(y, y_pred, labels=classes).tolist()

        print(f"  [{mode}] macro-F1: {macro_f1:.4f}")
        print(f"  [{mode}] accuracy: {accuracy:.4f}")
        for cls in classes:
            if cls in per_class:
                print(f"    {cls}: F1={per_class[cls]['f1-score']:.3f}")

        ci_low, ci_high = np.nan, np.nan
        if do_bootstrap:
            ci_low, ci_high = bootstrap_macro_f1_ci(y, y_pred)
            if not np.isnan(ci_low):
                print(f"  [{mode}] 95% CI: [{ci_low:.4f}, {ci_high:.4f}]")

        perm_p = permutation_test_macro_f1(y, y_pred, n_perm=n_perm)
        print(f"  [{mode}] permutation p-value: {perm_p:.6f}")

        results["feature_modes"][mode] = {
            "n_features": feat_aligned.shape[1],
            "macro_f1": macro_f1,
            "accuracy": accuracy,
            "ci_low": ci_low,
            "ci_high": ci_high,
            "perm_p": perm_p,
            "per_class_f1": {cls: per_class.get(cls, {}).get("f1-score", 0)
                            for cls in classes},
            "confusion_matrix": cm,
            "classes": classes,
            "y_true": y.tolist(),
            "y_pred": y_pred.tolist(),
        }

    return results


# ── Summary ──────────────────────────────────────────────────────────────────

def build_summary(all_results):
    """Build D category summary with J5 gene-vs-pathway comparison.

    Merges new results into existing summary (if present) so that
    running a single task (e.g. --task D6) does not overwrite results
    from other tasks.
    """
    summary_path = RESULTS_DIR / "D_condition_summary.json"

    # Load existing summary to preserve prior results
    if summary_path.exists():
        with open(summary_path) as f:
            summary = json.load(f)
        # Ensure required keys exist
        summary.setdefault("tasks", {})
        summary.setdefault("j5_comparison", {})
    else:
        summary = {"tasks": {}, "j5_comparison": {}}

    # Update config with current run parameters
    summary["config"] = {
        "variance_percentile": VARIANCE_PERCENTILE,
        "n_splits": N_SPLITS,
        "n_repeats": N_REPEATS,
        "n_bootstrap": N_BOOTSTRAP,
        "timestamp": datetime.now().isoformat(),
    }

    for task_name, res in all_results.items():
        if res is None:
            continue
        summary["tasks"][task_name] = res

        # J5: Gene vs Pathway comparison
        modes = res.get("feature_modes", {})
        if "gene" in modes and "pathway" in modes:
            gene_f1 = modes["gene"]["macro_f1"]
            path_f1 = modes["pathway"]["macro_f1"]
            summary["j5_comparison"][task_name] = {
                "gene_macro_f1": gene_f1,
                "pathway_macro_f1": path_f1,
                "diff_pathway_minus_gene": path_f1 - gene_f1,
                "winner": "pathway" if path_f1 > gene_f1 else "gene",
            }

    return summary


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    global N_BOOTSTRAP, N_SPLITS, N_REPEATS

    parser = argparse.ArgumentParser(
        description="Category D: Condition Prediction")
    parser.add_argument("--task", choices=["D3", "D4", "D5", "D6"],
                        help="Run specific task")
    parser.add_argument("--all", action="store_true",
                        help="Run all tasks")
    parser.add_argument("--no-bootstrap", action="store_true",
                        help="Skip bootstrap CI (faster)")
    parser.add_argument("--n-repeats", type=int, default=10,
                        help="Number of CV repeats for D3")
    args = parser.parse_args()

    if not args.task and not args.all:
        parser.print_help()
        return

    do_bootstrap = not args.no_bootstrap
    if args.no_bootstrap:
        N_BOOTSTRAP = 0
    N_REPEATS = args.n_repeats

    # Create output dirs
    D_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    all_results = {}

    # D3: Mission ID
    if args.task == "D3" or args.all:
        res = run_d3(do_bootstrap=do_bootstrap)
        all_results["D3"] = res
        # Save per-task result
        with open(D_OUTPUT_DIR / "D3_liver_results.json", "w") as f:
            json.dump(res, f, indent=2, default=str)

    # D4: Strain Effect
    if args.task == "D4" or args.all:
        res = run_d4(do_bootstrap=do_bootstrap)
        all_results["D4"] = res
        if res is not None:
            with open(D_OUTPUT_DIR / "D4_strain_results.json", "w") as f:
                json.dump(res, f, indent=2, default=str)

    # D5: Hardware Effect
    if args.task == "D5" or args.all:
        for tissue in ["liver", "thymus"]:
            res = run_d5(tissue, do_bootstrap=do_bootstrap)
            task_key = f"D5_{tissue}"
            all_results[task_key] = res
            if res is not None:
                with open(D_OUTPUT_DIR / f"D5_{tissue}_results.json", "w") as f:
                    json.dump(res, f, indent=2, default=str)

    # D6: Artificial Gravity
    if args.task == "D6" or args.all:
        for tissue in ["liver", "thymus"]:
            res = run_d6(tissue, do_bootstrap=do_bootstrap)
            task_key = f"D6_{tissue}"
            all_results[task_key] = res
            if res is not None:
                with open(D_OUTPUT_DIR / f"D6_{tissue}_results.json", "w") as f:
                    json.dump(res, f, indent=2, default=str)

    # Build and save summary
    summary = build_summary(all_results)
    with open(RESULTS_DIR / "D_condition_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    # Print J5 comparison
    print("\n" + "=" * 70)
    print("J5: Gene-level vs Pathway-level Comparison (Category D)")
    print("=" * 70)
    for task, comp in summary.get("j5_comparison", {}).items():
        print(f"\n  {task}:")
        print(f"    Gene    macro-F1: {comp['gene_macro_f1']:.4f}")
        print(f"    Pathway macro-F1: {comp['pathway_macro_f1']:.4f}")
        print(f"    Diff (P-G):       {comp['diff_pathway_minus_gene']:+.4f}")
        print(f"    Winner:           {comp['winner']}")

    print(f"\nResults saved to {RESULTS_DIR / 'D_condition_summary.json'}")


if __name__ == "__main__":
    main()
