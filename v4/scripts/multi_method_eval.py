#!/usr/bin/env python3
"""
multi_method_eval.py — GeneLabBench v4: Core multi-method evaluation engine

Runs any classifier from classifier_registry.py on any tissue with LOMO-CV
or stratified k-fold, reporting AUROC + bootstrap CI + permutation p-value.

Feature types: gene, pathway_hallmark, pathway_kegg, combined (gene+hallmark)

Usage:
  # Single evaluation
  python multi_method_eval.py --tissue liver --method pca_lr --features gene

  # All methods for one tissue
  python multi_method_eval.py --tissue liver --method all --features gene

  # Specific method + feature via SLURM array
  python multi_method_eval.py --tissue liver --method rf --features pathway_hallmark

  # Use v1 task folds (for reproduction)
  python multi_method_eval.py --tissue liver --method pca_lr --features gene --use-task-folds

  # High-precision permutation for best method
  python multi_method_eval.py --tissue liver --method pca_lr --features gene --n-perm 5000
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

warnings.filterwarnings("ignore")

# ── Imports from v4 modules ──────────────────────────────────────────────────
import sys
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
sys.path.insert(0, str(SCRIPT_DIR.parent.parent / "scripts"))  # for scripts/utils.py

from v4_utils import (
    TISSUE_MISSIONS, LOMO_TISSUES, KFOLD_TISSUES, TASK_MAP,
    BASE_DIR, V4_EVAL_DIR,
    load_metadata, load_gene_features, align_features_with_meta,
    get_folds, get_folds_from_task_dir,
    bootstrap_auroc, permutation_pvalue, save_evaluation_result,
)
from classifier_registry import (
    CLASSIFIERS, CPU_METHODS, GPU_METHODS,
    get_classifier, adapt_pca_components, fit_tabnet_with_eval,
)


# ── Pathway Loading ──────────────────────────────────────────────────────────

PATHWAY_DIR = BASE_DIR / "processed" / "pathway_scores"

def load_pathway_features(tissue, db="hallmark"):
    """Load GSVA pathway scores. Returns samples × pathways DataFrame."""
    all_scores = []
    for mission in TISSUE_MISSIONS.get(tissue, []):
        f = PATHWAY_DIR / tissue / f"{mission}_gsva_{db}.csv"
        if not f.exists():
            continue
        scores = pd.read_csv(f, index_col=0)
        all_scores.append(scores)

    if not all_scores:
        raise FileNotFoundError(
            f"No pathway scores found for tissue={tissue}, db={db}. "
            f"Checked: {PATHWAY_DIR / tissue}/"
        )

    combined = pd.concat(all_scores)
    if "mission" in combined.columns:
        combined = combined.drop(columns=["mission"])

    # Drop duplicates (MHU-2 fix)
    combined = combined[~combined.index.duplicated(keep="first")]

    return combined


def load_features(tissue, feature_type):
    """Load features by type. Returns (features_df, feature_type_label)."""
    if feature_type == "gene":
        return load_gene_features(tissue), "gene"
    elif feature_type == "pathway_hallmark":
        return load_pathway_features(tissue, "hallmark"), "pathway_hallmark"
    elif feature_type == "pathway_kegg":
        return load_pathway_features(tissue, "kegg"), "pathway_kegg"
    elif feature_type == "combined":
        genes = load_gene_features(tissue)
        pathways = load_pathway_features(tissue, "hallmark")
        # Align samples
        common_samples = sorted(set(genes.index) & set(pathways.index))
        if len(common_samples) < 5:
            raise ValueError(f"Too few common samples for combined features: {len(common_samples)}")
        combined = pd.concat([genes.loc[common_samples], pathways.loc[common_samples]], axis=1)
        return combined, "combined"
    else:
        raise ValueError(f"Unknown feature_type: {feature_type}")


# ── Core Evaluation ──────────────────────────────────────────────────────────

def evaluate_single(tissue, method_key, feature_type,
                    n_bootstrap=2000, n_perm=1000,
                    use_task_folds=False, seed=42, verbose=True):
    """Run one method on one tissue with one feature type.

    Returns result dict matching v1 A1_baseline_results.json schema.
    """
    t_start = time.time()
    label, _ = get_classifier(method_key, seed)

    if verbose:
        print(f"\n{'='*60}")
        print(f"Tissue: {tissue} | Method: {label} ({method_key})")
        print(f"Features: {feature_type} | Bootstrap: {n_bootstrap} | Perm: {n_perm}")
        cv_type = "LOMO" if tissue in LOMO_TISSUES else "5-fold stratified"
        print(f"CV: {cv_type}")
        if use_task_folds:
            print(f"Folds: from tasks/ directory (v1 reproduction)")
        print(f"{'='*60}")

    # Load features
    if feature_type == "gene" and use_task_folds and tissue in TASK_MAP:
        # Use pre-existing task folds (exact v1 reproduction)
        folds = get_folds_from_task_dir(tissue)
        if verbose:
            print(f"Loaded {len(folds)} folds from tasks/{TASK_MAP[tissue]}/")
    else:
        # Generate folds from raw data
        features = load_features(tissue, feature_type)[0]
        meta = load_metadata(tissue)
        features, meta = align_features_with_meta(features, meta)
        folds = get_folds(tissue, gene_features=features, meta=meta)
        if verbose:
            print(f"Generated {len(folds)} folds from raw data "
                  f"({features.shape[0]} samples × {features.shape[1]} features)")

    fold_results = []
    peak_memory_mb = 0

    for fold in folds:
        fold_name = fold["fold_name"]
        test_mission = fold["test_mission"]

        X_train = fold["train_X"]
        y_train = fold["train_y"]
        X_test = fold["test_X"]
        y_test = fold["test_y"]

        if len(np.unique(y_test)) < 2:
            if verbose:
                print(f"  [SKIP] {test_mission}: only one class in test set")
            continue

        # Build fresh model for each fold
        _, model = get_classifier(method_key, seed)

        # Adapt PCA components for small training sets
        adapt_pca_components(model, X_train.shape[0], X_train.shape[1])

        # Train
        t0 = time.time()
        try:
            if method_key == "tabnet":
                # TabNet needs special handling
                # Use 10% of training for validation
                n_val = max(2, int(0.1 * len(y_train)))
                fit_tabnet_with_eval(
                    model, X_train[:-n_val], y_train[:-n_val].astype(np.int64),
                    X_train[-n_val:], y_train[-n_val:].astype(np.int64),
                    max_epochs=200, patience=20
                )
            else:
                model.fit(X_train, y_train)
        except Exception as e:
            if verbose:
                print(f"  [ERROR] {test_mission}: training failed — {e}")
                traceback.print_exc()
            continue
        train_time = time.time() - t0

        # Predict
        if hasattr(model, "predict_proba"):
            y_score = model.predict_proba(X_test)[:, 1]
        else:
            y_score = model.decision_function(X_test)

        # Metrics
        from sklearn.metrics import roc_auc_score
        auroc = float(roc_auc_score(y_test, y_score))
        mean_bs, lower_bs, upper_bs = bootstrap_auroc(y_test, y_score, n_bootstrap)
        pval = permutation_pvalue(y_test, y_score, n_perm)

        result = {
            "test_mission": test_mission,
            "fold_name": fold_name,
            "auroc": round(auroc, 4),
            "bootstrap_mean": round(mean_bs, 4),
            "ci_lower": round(lower_bs, 4),
            "ci_upper": round(upper_bs, 4),
            "perm_pvalue": round(pval, 4),
            "n_train": fold["n_train"],
            "n_test": fold["n_test"],
            "n_features": X_train.shape[1],
            "n_flight_test": int(y_test.sum()),
            "n_ground_test": int((y_test == 0).sum()),
            "train_time_s": round(train_time, 2),
        }
        fold_results.append(result)

        if verbose:
            print(f"  {test_mission}: AUROC={auroc:.3f} "
                  f"[{lower_bs:.3f}–{upper_bs:.3f}] "
                  f"p={pval:.3f} (n={fold['n_test']}, t={train_time:.1f}s)")

        # Track memory
        try:
            import resource
            mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024  # KB → MB on Linux
            peak_memory_mb = max(peak_memory_mb, mem)
        except Exception:
            pass

    if not fold_results:
        if verbose:
            print(f"  [ERROR] No valid folds for {tissue}/{method_key}/{feature_type}")
        return None

    # Aggregate across folds
    aurocs = [r["auroc"] for r in fold_results]
    pvals = [r["perm_pvalue"] for r in fold_results]
    ci_lowers = [r["ci_lower"] for r in fold_results]

    total_time = time.time() - t_start

    summary = {
        "model": label,
        "model_key": method_key,
        "tissue": tissue,
        "feature_type": feature_type,
        "cv_strategy": "LOMO" if tissue in LOMO_TISSUES else "5-fold_stratified",
        "n_folds": len(fold_results),
        "mean_auroc": round(float(np.mean(aurocs)), 4),
        "std_auroc": round(float(np.std(aurocs)), 4),
        "min_auroc": round(float(np.min(aurocs)), 4),
        "max_auroc": round(float(np.max(aurocs)), 4),
        "mean_ci_lower": round(float(np.mean(ci_lowers)), 4),
        "mean_perm_pvalue": round(float(np.mean(pvals)), 4),
        "n_bootstrap": n_bootstrap,
        "n_perm": n_perm,
        "wall_time_sec": round(total_time, 1),
        "peak_memory_mb": round(peak_memory_mb, 1),
        "timestamp": datetime.now().isoformat(),
        "folds": fold_results,
    }

    if verbose:
        print(f"\n  ── {label} Summary ──")
        print(f"  Mean AUROC = {summary['mean_auroc']:.4f} ± {summary['std_auroc']:.4f}")
        print(f"  CI lower = {summary['mean_ci_lower']:.4f}")
        print(f"  Perm p = {summary['mean_perm_pvalue']:.4f}")
        print(f"  Wall time: {total_time:.1f}s")

    return summary


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="GeneLabBench v4 Multi-Method Evaluation")
    parser.add_argument("--tissue", required=True,
                        choices=list(TISSUE_MISSIONS.keys()),
                        help="Tissue to evaluate")
    parser.add_argument("--method", required=True,
                        help="Classifier key (or 'all' / 'cpu' / 'gpu')")
    parser.add_argument("--features", required=True,
                        choices=["gene", "pathway_hallmark", "pathway_kegg", "combined", "all"],
                        help="Feature type")
    parser.add_argument("--n-bootstrap", type=int, default=2000,
                        help="Bootstrap iterations for CI (default: 2000)")
    parser.add_argument("--n-perm", type=int, default=1000,
                        help="Permutation iterations for p-value (default: 1000)")
    parser.add_argument("--use-task-folds", action="store_true",
                        help="Use v1 pre-split task/ folds (gene features only)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output-dir", type=str, default=None,
                        help="Output directory (default: v4/evaluation/)")
    parser.add_argument("--quiet", action="store_true")

    args = parser.parse_args()
    verbose = not args.quiet
    output_dir = Path(args.output_dir) if args.output_dir else V4_EVAL_DIR

    # Resolve method list
    if args.method == "all":
        methods = list(CLASSIFIERS.keys())
    elif args.method == "cpu":
        methods = CPU_METHODS
    elif args.method == "gpu":
        methods = GPU_METHODS
    else:
        methods = [args.method]

    # Resolve feature list
    if args.features == "all":
        feature_types = ["gene", "pathway_hallmark", "pathway_kegg", "combined"]
    else:
        feature_types = [args.features]

    # Run evaluations
    all_results = {}
    for feature_type in feature_types:
        for method_key in methods:
            if method_key not in CLASSIFIERS:
                print(f"[WARN] Unknown method: {method_key}")
                continue

            result_key = f"{method_key}_{feature_type}"
            try:
                result = evaluate_single(
                    tissue=args.tissue,
                    method_key=method_key,
                    feature_type=feature_type,
                    n_bootstrap=args.n_bootstrap,
                    n_perm=args.n_perm,
                    use_task_folds=args.use_task_folds,
                    seed=args.seed,
                    verbose=verbose,
                )
                if result:
                    all_results[result_key] = result

                    # Save individual result
                    fname = f"M1_{args.tissue}_{feature_type}_{method_key}.json"
                    output_dir.mkdir(parents=True, exist_ok=True)
                    (output_dir / fname).write_text(json.dumps(result, indent=2))
                    if verbose:
                        print(f"  Saved: {output_dir / fname}")

            except FileNotFoundError as e:
                print(f"[SKIP] {args.tissue}/{feature_type}/{method_key}: {e}")
            except Exception as e:
                print(f"[ERROR] {args.tissue}/{feature_type}/{method_key}: {e}")
                traceback.print_exc()

    # Save combined results for this tissue
    if all_results:
        combined_fname = f"M1_{args.tissue}_combined_results.json"
        (output_dir / combined_fname).write_text(json.dumps(all_results, indent=2))
        if verbose:
            print(f"\nCombined results: {output_dir / combined_fname}")
            print(f"\nTotal evaluations: {len(all_results)}")
            for k, v in all_results.items():
                print(f"  {k}: AUROC={v['mean_auroc']:.4f}")


if __name__ == "__main__":
    main()
