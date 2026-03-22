#!/usr/bin/env python3
"""
condition_prediction_v4.py — GeneLabBench v4: D3/D5 Multi-Method Extension

Extends v1 condition_prediction.py to all 8 classifiers.

D3: Mission ID prediction (multi-class)
  - 6 LOMO tissues with mission labels
  - RepeatedStratifiedKFold (5-fold × 10 repeats)
  - Metric: macro-F1 (primary), accuracy, per-class F1

D5: Hardware type prediction (binary: RR vs MHU)
  - Only tissues with BOTH RR and MHU missions: liver, thymus
  - RepeatedStratifiedKFold (5-fold × 10 repeats)

Each: 8 methods × {gene, pathway_hallmark} features

Usage:
  python condition_prediction_v4.py --task D3 --tissue liver --method pca_lr --features gene
  python condition_prediction_v4.py --task D5 --tissue liver --method all --features gene
"""

import json
import time
import argparse
import warnings
import traceback
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from sklearn.metrics import f1_score, accuracy_score, confusion_matrix
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.preprocessing import LabelEncoder

warnings.filterwarnings("ignore")

import sys
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from v4_utils import (
    TISSUE_MISSIONS, LOMO_TISSUES,
    BASE_DIR, V4_EVAL_DIR,
    load_metadata, load_gene_features, align_features_with_meta,
)
from classifier_registry import (
    CLASSIFIERS, CPU_METHODS, GPU_METHODS,
    get_classifier, adapt_pca_components, fit_tabnet_with_eval,
)

# ── Config ────────────────────────────────────────────────────────────────────

PATHWAY_DIR = BASE_DIR / "processed" / "pathway_scores"

HARDWARE_MAP = {
    "RR-1": "RR", "RR-3": "RR", "RR-5": "RR",
    "RR-6": "RR", "RR-7": "RR", "RR-8": "RR", "RR-9": "RR",
    "MHU-1": "MHU", "MHU-2": "MHU",
    "TBD": "unknown",
}

# D3 applicable tissues: multi-mission LOMO tissues
D3_TISSUES = LOMO_TISSUES  # liver, gastro, kidney, thymus, eye, skin

# D5 applicable tissues: must have BOTH RR and MHU missions
D5_TISSUES = ["liver", "thymus"]

N_SPLITS = 5
N_REPEATS = 10
N_BOOTSTRAP = 2000
N_PERM = 1000


# ── Feature Loading ──────────────────────────────────────────────────────────

def load_pathway_features(tissue, db="hallmark"):
    """Load GSVA pathway scores."""
    all_scores = []
    for mission in TISSUE_MISSIONS.get(tissue, []):
        f = PATHWAY_DIR / tissue / f"{mission}_gsva_{db}.csv"
        if not f.exists():
            continue
        scores = pd.read_csv(f, index_col=0)
        all_scores.append(scores)

    if not all_scores:
        return None

    combined = pd.concat(all_scores)
    if "mission" in combined.columns:
        combined = combined.drop(columns=["mission"])
    combined = combined[~combined.index.duplicated(keep="first")]
    return combined


def load_features(tissue, feature_type):
    """Load features by type."""
    if feature_type == "gene":
        return load_gene_features(tissue)
    elif feature_type == "pathway_hallmark":
        return load_pathway_features(tissue, "hallmark")
    else:
        raise ValueError(f"Unknown feature_type: {feature_type}")


# ── Statistical Utilities ────────────────────────────────────────────────────

def bootstrap_macro_f1_ci(y_true, y_pred, n_boot=N_BOOTSTRAP, alpha=0.05, seed=42):
    """Bootstrap 95% CI for macro-F1."""
    rng = np.random.default_rng(seed)
    n = len(y_true)
    boot_f1s = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        yt, yp = y_true[idx], y_pred[idx]
        if len(np.unique(yt)) < 2:
            continue
        boot_f1s.append(f1_score(yt, yp, average="macro", zero_division=0))

    if len(boot_f1s) < 100:
        return float("nan"), float("nan")

    lo = float(np.percentile(boot_f1s, 100 * alpha / 2))
    hi = float(np.percentile(boot_f1s, 100 * (1 - alpha / 2)))
    return lo, hi


def permutation_pvalue_f1(y_true, y_pred, n_perm=N_PERM, seed=42):
    """Permutation p-value for macro-F1."""
    observed = f1_score(y_true, y_pred, average="macro", zero_division=0)
    rng = np.random.default_rng(seed)
    count_ge = 0
    for _ in range(n_perm):
        y_perm = rng.permutation(y_true)
        pf1 = f1_score(y_perm, y_pred, average="macro", zero_division=0)
        if pf1 >= observed:
            count_ge += 1
    return (count_ge + 1) / (n_perm + 1)


# ── Core Evaluation ──────────────────────────────────────────────────────────

def evaluate_condition(tissue, method_key, feature_type, task="D3",
                       n_splits=N_SPLITS, n_repeats=N_REPEATS,
                       n_bootstrap=N_BOOTSTRAP, n_perm=N_PERM,
                       seed=42, verbose=True):
    """Run D3 or D5 evaluation with one method on one tissue."""
    t_start = time.time()
    label, _ = get_classifier(method_key, seed)

    if verbose:
        print(f"\n{'='*60}")
        print(f"{task} | Tissue: {tissue} | Method: {label} | Features: {feature_type}")
        print(f"CV: {n_splits}-fold × {n_repeats}-repeat stratified")
        print(f"{'='*60}")

    # Load features
    feat_df = load_features(tissue, feature_type)
    if feat_df is None:
        print(f"  [SKIP] No {feature_type} features for {tissue}")
        return None

    meta = load_metadata(tissue)
    feat_df, meta = align_features_with_meta(feat_df, meta)

    # Extract labels based on task
    if task == "D3":
        y_raw = meta["mission"].values
    elif task == "D5":
        meta = meta.copy()
        meta["hardware"] = meta["mission"].map(HARDWARE_MAP)
        valid = meta["hardware"] != "unknown"
        meta = meta[valid]
        feat_df = feat_df[valid]
        y_raw = meta["hardware"].values
    else:
        raise ValueError(f"Unknown task: {task}")

    classes = sorted(np.unique(y_raw))
    n_classes = len(classes)

    if n_classes < 2:
        print(f"  [SKIP] Only {n_classes} class(es)")
        return None

    # Encode labels to integers for classifiers
    le = LabelEncoder()
    y = le.fit_transform(y_raw)
    X = feat_df.values.astype(np.float32)
    X = np.nan_to_num(X, nan=0.0)

    if verbose:
        from collections import Counter
        dist = Counter(y_raw)
        print(f"  Samples: {len(y)}, Classes: {n_classes}")
        for cls in classes:
            print(f"    {cls}: {dist[cls]}")

    # Repeated Stratified K-Fold
    rskf = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats,
                                   random_state=seed)

    # Accumulate predictions across repeats
    prob_accum = np.zeros((len(y), n_classes))
    count_accum = np.zeros(len(y))
    train_times = []

    # For gene features, apply PCA preprocessing to reduce dimensionality
    # before classifiers that don't include PCA internally.
    # This matches v1 condition_prediction.py design and avoids infeasible
    # SAGA convergence on 22K raw features (50-fold CV).
    # Methods with internal PCA (pca_lr) or tree-based (rf, xgb) don't need this.
    NEEDS_GENE_PCA = {"elasticnet_lr", "svm_rbf", "knn", "mlp", "tabnet"}
    apply_gene_pca = (feature_type == "gene" and method_key in NEEDS_GENE_PCA)

    for train_idx, test_idx in rskf.split(X, y):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train = y[train_idx]

        # Dimensionality reduction for gene features (matches v1 approach)
        if apply_gene_pca:
            from sklearn.preprocessing import StandardScaler as SS
            from sklearn.decomposition import PCA
            n_comp = min(50, X_train.shape[0] - 1, X_train.shape[1])
            scaler_pre = SS()
            X_train = scaler_pre.fit_transform(X_train)
            X_test = scaler_pre.transform(X_test)
            pca_pre = PCA(n_components=n_comp, random_state=seed)
            X_train = pca_pre.fit_transform(X_train)
            X_test = pca_pre.transform(X_test)

        _, model = get_classifier(method_key, seed)
        adapt_pca_components(model, X_train.shape[0], X_train.shape[1])

        t0 = time.time()
        try:
            if method_key == "tabnet":
                n_val = max(2, int(0.1 * len(y_train)))
                fit_tabnet_with_eval(
                    model, X_train[:-n_val], y_train[:-n_val].astype(np.int64),
                    X_train[-n_val:], y_train[-n_val:].astype(np.int64),
                    max_epochs=200, patience=20,
                )
            else:
                model.fit(X_train, y_train)
        except Exception as e:
            if verbose:
                print(f"  [ERROR] fold: {e}")
            continue
        train_times.append(time.time() - t0)

        # Predict
        if hasattr(model, "predict_proba"):
            proba = model.predict_proba(X_test)
        else:
            # For models without predict_proba, use one-hot from predict
            pred = model.predict(X_test)
            proba = np.zeros((len(pred), n_classes))
            for i, p in enumerate(pred):
                proba[i, p] = 1.0

        # Map model classes to our class order
        if hasattr(model, "classes_"):
            model_classes = model.classes_
        elif hasattr(model, "named_steps"):
            # Pipeline: get classes from last step
            last = model.named_steps.get("clf", model)
            model_classes = getattr(last, "classes_", np.arange(n_classes))
        else:
            model_classes = np.arange(n_classes)

        for i, cls in enumerate(model_classes):
            if i < proba.shape[1]:
                prob_accum[test_idx, cls] += proba[:, i]
        count_accum[test_idx] += 1

    # Average predictions and compute final metrics
    valid_mask = count_accum > 0
    if not valid_mask.all():
        print(f"  [WARN] {(~valid_mask).sum()} samples never in test set")

    prob_avg = prob_accum[valid_mask] / count_accum[valid_mask, None]
    y_pred = np.argmax(prob_avg, axis=1)
    y_true = y[valid_mask]

    # Convert back to original labels for reporting
    y_true_labels = le.inverse_transform(y_true)
    y_pred_labels = le.inverse_transform(y_pred)

    macro_f1 = float(f1_score(y_true, y_pred, average="macro", zero_division=0))
    accuracy = float(accuracy_score(y_true, y_pred))
    cm = confusion_matrix(y_true, y_pred).tolist()

    # Per-class F1
    per_class_f1 = {}
    for cls_idx, cls_name in enumerate(classes):
        binary_true = (y_true_labels == cls_name).astype(int)
        binary_pred = (y_pred_labels == cls_name).astype(int)
        per_class_f1[cls_name] = round(float(f1_score(binary_true, binary_pred, zero_division=0)), 4)

    # Bootstrap CI
    ci_lower, ci_upper = bootstrap_macro_f1_ci(y_true, y_pred, n_bootstrap, seed=seed)

    # Permutation test
    perm_p = permutation_pvalue_f1(y_true, y_pred, n_perm, seed=seed)

    total_time = time.time() - t_start

    result = {
        "experiment": task,
        "model": label,
        "model_key": method_key,
        "tissue": tissue,
        "feature_type": feature_type,
        "cv_strategy": f"repeated_stratified_{n_splits}x{n_repeats}",
        "n_classes": n_classes,
        "classes": classes,
        "class_distribution": {cls: int((y_raw == cls).sum()) for cls in classes},
        "n_samples": len(y),
        "macro_f1": round(macro_f1, 4),
        "accuracy": round(accuracy, 4),
        "ci_lower": round(ci_lower, 4) if not np.isnan(ci_lower) else None,
        "ci_upper": round(ci_upper, 4) if not np.isnan(ci_upper) else None,
        "perm_pvalue": round(perm_p, 4),
        "per_class_f1": per_class_f1,
        "confusion_matrix": cm,
        "mean_train_time_s": round(float(np.mean(train_times)), 2) if train_times else 0,
        "n_bootstrap": n_bootstrap,
        "n_perm": n_perm,
        "seed": seed,
        "wall_time_sec": round(total_time, 1),
        "timestamp": datetime.now().isoformat(),
    }

    if verbose:
        print(f"\n  macro-F1 = {macro_f1:.4f}")
        print(f"  accuracy = {accuracy:.4f}")
        if ci_lower is not None:
            print(f"  CI: [{ci_lower:.4f}, {ci_upper:.4f}]")
        print(f"  perm p = {perm_p:.4f}")
        for cls, f1val in per_class_f1.items():
            print(f"    {cls}: F1={f1val:.3f}")

    return result


def main():
    parser = argparse.ArgumentParser(description="D3/D5 Multi-Method Condition Prediction")
    parser.add_argument("--task", required=True, choices=["D3", "D5"])
    parser.add_argument("--tissue", required=True)
    parser.add_argument("--method", required=True,
                        help="Classifier key (or 'all' / 'cpu' / 'gpu')")
    parser.add_argument("--features", required=True,
                        choices=["gene", "pathway_hallmark"])
    parser.add_argument("--n-splits", type=int, default=N_SPLITS)
    parser.add_argument("--n-repeats", type=int, default=N_REPEATS)
    parser.add_argument("--n-bootstrap", type=int, default=N_BOOTSTRAP)
    parser.add_argument("--n-perm", type=int, default=N_PERM)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output-dir", type=str, default=None)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    # Validate tissue for task
    if args.task == "D3" and args.tissue not in D3_TISSUES:
        print(f"[ERROR] D3 requires multi-mission LOMO tissue. Valid: {D3_TISSUES}")
        return
    if args.task == "D5" and args.tissue not in D5_TISSUES:
        print(f"[ERROR] D5 requires tissue with both RR+MHU. Valid: {D5_TISSUES}")
        return

    verbose = not args.quiet
    output_dir = Path(args.output_dir) if args.output_dir else V4_EVAL_DIR

    if args.method == "all":
        methods = list(CLASSIFIERS.keys())
    elif args.method == "cpu":
        methods = CPU_METHODS
    elif args.method == "gpu":
        methods = GPU_METHODS
    else:
        methods = [args.method]

    for method_key in methods:
        if method_key not in CLASSIFIERS:
            print(f"[WARN] Unknown method: {method_key}")
            continue

        try:
            result = evaluate_condition(
                tissue=args.tissue,
                method_key=method_key,
                feature_type=args.features,
                task=args.task,
                n_splits=args.n_splits,
                n_repeats=args.n_repeats,
                n_bootstrap=args.n_bootstrap,
                n_perm=args.n_perm,
                seed=args.seed,
                verbose=verbose,
            )
            if result:
                fname = f"{args.task}_v4_{args.tissue}_{method_key}_{args.features}.json"
                output_dir.mkdir(parents=True, exist_ok=True)
                (output_dir / fname).write_text(json.dumps(result, indent=2))
                if verbose:
                    print(f"  Saved: {output_dir / fname}")
        except Exception as e:
            print(f"[ERROR] {args.tissue}/{method_key}: {e}")
            traceback.print_exc()


if __name__ == "__main__":
    main()
