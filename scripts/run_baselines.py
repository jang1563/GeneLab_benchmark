#!/usr/bin/env python3
"""
run_baselines.py — GeneLab_benchmark: Baseline Model Evaluation

Runs baseline classifiers on LOMO task splits and reports AUROC with
Bootstrap 95% CI and permutation p-values (DD-08).

Baselines (DD-13):
  - Logistic Regression (ElasticNet, SAGA solver)
  - Random Forest
  - XGBoost
  - PCA + LogReg (50 PCs, low-dimensional baseline)

Metrics (DD-08):
  - AUROC (primary) with 1000-bootstrap 95% CI
  - Permutation p-value (H0: AUROC = 0.5)
  - LOMO fold mean ± SD

Phase 1 go/no-go check (DD-11):
  1. Mean AUROC > 0.70 AND 95% CI lower bound > 0.60
  2. Permutation p < 0.05 (mean across folds)
  3. Housekeeping gene rank check (ANGPTL4, PCK1 in SHAP top-50)

Usage:
  python run_baselines.py --task A1           # run all baselines for A1
  python run_baselines.py --task A1 --model lr  # specific model
  python run_baselines.py --task A1 --quick   # skip XGBoost for speed
  python run_baselines.py --task A1 --check   # check go/no-go conditions
"""

import json
import time
import argparse
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

warnings.filterwarnings("ignore")

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent.parent
TASKS_DIR = BASE_DIR / "tasks"
RESULTS_DIR = BASE_DIR / "evaluation"

# ── Metric functions (DD-08) ───────────────────────────────────────────────────

def bootstrap_auroc(y_true: np.ndarray, y_score: np.ndarray,
                    n_bootstrap: int = 1000, alpha: float = 0.05,
                    seed: int = 42) -> tuple[float, float, float]:
    """
    Bootstrap 95% CI for AUROC.
    Returns (mean_auroc, lower_ci, upper_ci).
    """
    from sklearn.metrics import roc_auc_score
    rng = np.random.default_rng(seed)
    n = len(y_true)
    scores = []
    for _ in range(n_bootstrap):
        idx = rng.choice(n, size=n, replace=True)
        yt, ys = y_true[idx], y_score[idx]
        if len(np.unique(yt)) < 2:
            continue
        scores.append(roc_auc_score(yt, ys))

    if not scores:
        return float("nan"), float("nan"), float("nan")

    lower = float(np.percentile(scores, 100 * alpha / 2))
    upper = float(np.percentile(scores, 100 * (1 - alpha / 2)))
    return float(np.mean(scores)), lower, upper


def permutation_pvalue(y_true: np.ndarray, y_score: np.ndarray,
                       n_perm: int = 1000, seed: int = 42) -> float:
    """
    Permutation test p-value for AUROC.
    H0: AUROC = random (0.5).
    """
    from sklearn.metrics import roc_auc_score
    rng = np.random.default_rng(seed)
    observed = roc_auc_score(y_true, y_score)
    null_scores = []
    for _ in range(n_perm):
        y_perm = rng.permutation(y_true)
        null_scores.append(roc_auc_score(y_perm, y_score))
    p = float((np.sum(np.array(null_scores) >= observed) + 1) / (n_perm + 1))
    return p


# ── Model builders ─────────────────────────────────────────────────────────────

def build_lr(seed: int = 42):
    """Logistic Regression with ElasticNet penalty (DD-13).
    Note: penalty= kwarg removed for sklearn≥1.8 compatibility (use l1_ratio= instead).
    max_iter=10000 ensures SAGA convergence on high-dim data (≥15k genes). (B3 fix)
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import Pipeline
    return Pipeline([
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(
            solver="saga",
            l1_ratio=0.5,
            C=1.0,
            class_weight="balanced",
            max_iter=10000,
            random_state=seed,
        ))
    ])


def build_rf(seed: int = 42):
    """Random Forest (DD-13)."""
    from sklearn.ensemble import RandomForestClassifier
    return RandomForestClassifier(
        n_estimators=200,
        max_features="sqrt",
        class_weight="balanced",
        n_jobs=-1,
        random_state=seed,
    )


def build_pca_lr(n_components: int = 50, seed: int = 42):
    """PCA (up to 50 PCs) + Logistic Regression (DD-13).
    n_components is set adaptively in evaluate_fold() based on train set size.
    max_iter=5000: lbfgs converges fast on PCA-reduced space; 5000 provides headroom. (B3 fix)
    """
    from sklearn.decomposition import PCA
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import Pipeline
    return Pipeline([
        ("scaler", StandardScaler()),
        ("pca", PCA(n_components=n_components, random_state=seed)),
        ("clf", LogisticRegression(
            C=1.0,
            class_weight="balanced",
            max_iter=5000,
            random_state=seed,
        ))
    ])


def build_xgb(seed: int = 42):
    """XGBoost (DD-13)."""
    try:
        from xgboost import XGBClassifier
        return XGBClassifier(
            n_estimators=200,
            max_depth=3,
            subsample=0.7,
            learning_rate=0.1,
            scale_pos_weight=1,
            eval_metric="logloss",
            verbosity=0,
            random_state=seed,
        )
    except ImportError:
        return None


MODELS = {
    "lr":      ("Logistic Regression (ElasticNet)", build_lr),
    "rf":      ("Random Forest", build_rf),
    "pca_lr":  ("PCA-50 + LogReg", build_pca_lr),
    "xgb":     ("XGBoost", build_xgb),
}


# ── Per-fold evaluation ────────────────────────────────────────────────────────

def evaluate_fold(fold_dir: Path, model_name: str, model,
                  n_bootstrap: int = 1000, n_perm: int = 1000,
                  verbose: bool = True) -> dict | None:
    """
    Evaluate one model on one LOMO fold.
    Returns result dict with AUROC, CI, p-value.
    """
    from sklearn.metrics import roc_auc_score

    # Load data
    train_X = pd.read_csv(fold_dir / "train_X.csv", index_col=0)
    test_X  = pd.read_csv(fold_dir / "test_X.csv",  index_col=0)
    train_y = pd.read_csv(fold_dir / "train_y.csv", index_col=0).squeeze()
    test_y  = pd.read_csv(fold_dir / "test_y.csv",  index_col=0).squeeze()

    fold_info_path = fold_dir / "fold_info.json"
    fold_info = json.loads(fold_info_path.read_text()) if fold_info_path.exists() else {}
    test_mission = fold_info.get("test_mission", fold_dir.name)

    # Align features
    common_genes = train_X.columns.intersection(test_X.columns)
    train_X = train_X[common_genes]
    test_X  = test_X[common_genes]

    X_train = train_X.values.astype(np.float32)
    X_test  = test_X.values.astype(np.float32)
    y_train = train_y.values.astype(int)
    y_test  = test_y.values.astype(int)

    # Adapt PCA n_components to training set size (avoids n_components > n_train-1 error)
    from sklearn.pipeline import Pipeline
    from sklearn.decomposition import PCA as _PCA
    if isinstance(model, Pipeline):
        for _, step in model.steps:
            if isinstance(step, _PCA):
                max_comps = min(step.n_components, X_train.shape[0] - 1, X_train.shape[1])
                if max_comps < step.n_components:
                    if verbose:
                        print(f"    [INFO] PCA n_components: {step.n_components} → {max_comps} "
                              f"(n_train={X_train.shape[0]})")
                    step.n_components = max_comps

    if len(np.unique(y_test)) < 2:
        if verbose:
            print(f"    [SKIP] {test_mission}: only one class in test set")
        return None

    # Train
    t0 = time.time()
    try:
        model.fit(X_train, y_train)
    except Exception as e:
        if verbose:
            print(f"    [ERROR] {test_mission}: training failed — {e}")
        return None
    train_time = time.time() - t0

    # Predict
    if hasattr(model, "predict_proba"):
        y_score = model.predict_proba(X_test)[:, 1]
    else:
        y_score = model.decision_function(X_test)

    # Metrics
    auroc = float(roc_auc_score(y_test, y_score))
    mean_bs, lower_bs, upper_bs = bootstrap_auroc(y_test, y_score, n_bootstrap)
    pval = permutation_pvalue(y_test, y_score, n_perm)

    result = {
        "test_mission": test_mission,
        "model": model_name,
        "auroc": auroc,
        "bootstrap_mean": mean_bs,
        "ci_lower": lower_bs,
        "ci_upper": upper_bs,
        "perm_pvalue": pval,
        "n_train": len(y_train),
        "n_test": len(y_test),
        "n_genes": len(common_genes),
        "n_flight_test": int(y_test.sum()),
        "n_ground_test": int((y_test == 0).sum()),
        "train_time_s": round(train_time, 2),
    }

    if verbose:
        print(f"    {test_mission}: AUROC={auroc:.3f} [CI: {lower_bs:.3f}–{upper_bs:.3f}] "
              f"perm_p={pval:.3f}  (n_test={len(y_test)}, t={train_time:.1f}s)")

    return result


# ── Task-level evaluation ──────────────────────────────────────────────────────

def evaluate_task(task_id: str, model_names: list[str],
                  n_bootstrap: int = 1000, n_perm: int = 1000,
                  verbose: bool = True, combat: bool = False) -> dict:
    """Run all models on all folds for a task."""

    # Find task directory — support _lomo and _lomo_combat variants
    suffix = "_combat" if combat else ""
    task_dirs = list(TASKS_DIR.glob(f"{task_id}_*_lomo{suffix}"))
    if not task_dirs:
        print(f"[ERROR] No task directory found for {task_id}")
        print(f"        Run: python scripts/generate_tasks.py --task {task_id}")
        return {}

    task_dir = task_dirs[0]
    fold_dirs = sorted([d for d in task_dir.iterdir() if d.is_dir() and d.name.startswith("fold_")])

    if not fold_dirs:
        print(f"[ERROR] No folds found in {task_dir}")
        return {}

    print(f"\n{'='*60}")
    print(f"Task {task_id}: {task_dir.name}")
    print(f"Folds: {len(fold_dirs)}")
    print(f"{'='*60}")

    all_results = {}

    for model_name in model_names:
        if model_name not in MODELS:
            print(f"[WARN] Unknown model: {model_name}")
            continue

        model_label, model_builder = MODELS[model_name]
        model = model_builder()

        if model is None:
            print(f"  [SKIP] {model_label} — not installed")
            continue

        print(f"\n  Model: {model_label}")
        fold_results = []

        for fold_dir in fold_dirs:
            # Re-build model for each fold (fresh state)
            model = model_builder()
            result = evaluate_fold(fold_dir, model_label, model,
                                   n_bootstrap=n_bootstrap, n_perm=n_perm,
                                   verbose=verbose)
            if result:
                fold_results.append(result)

        if not fold_results:
            continue

        # Aggregate across folds
        aurocs = [r["auroc"] for r in fold_results]
        pvals  = [r["perm_pvalue"] for r in fold_results]
        ci_lowers = [r["ci_lower"] for r in fold_results]

        summary = {
            "model": model_label,
            "model_key": model_name,
            "n_folds": len(fold_results),
            "mean_auroc": float(np.mean(aurocs)),
            "std_auroc": float(np.std(aurocs)),
            "min_auroc": float(np.min(aurocs)),
            "max_auroc": float(np.max(aurocs)),
            "mean_ci_lower": float(np.mean(ci_lowers)),
            "mean_perm_pvalue": float(np.mean(pvals)),
            "folds": fold_results,
        }

        print(f"\n  ── {model_label} Summary ──")
        print(f"  LOMO mean AUROC = {summary['mean_auroc']:.3f} ± {summary['std_auroc']:.3f}")
        print(f"  Mean CI lower = {summary['mean_ci_lower']:.3f}")
        print(f"  Mean perm p = {summary['mean_perm_pvalue']:.4f}")

        all_results[model_name] = summary

    return all_results


# ── Phase 1 Go/No-Go check (DD-11) ────────────────────────────────────────────

def check_gonogo(results: dict, task_id: str = "A1") -> bool:
    """
    DD-11: Phase 1 go/no-go check.
    Conditions (ALL must pass):
      1. mean AUROC > 0.70 AND mean CI lower > 0.60
      2. mean permutation p < 0.05
      3. (SHAP check — requires separate SHAP computation)
    """
    print(f"\n{'='*60}")
    print(f"Phase 1 Go/No-Go Check — Task {task_id}")
    print(f"{'='*60}")

    if not results:
        print("  [FAIL] No results to evaluate.")
        return False

    # Use best model by mean AUROC
    best_model_key = max(results.keys(), key=lambda k: results[k]["mean_auroc"])
    best = results[best_model_key]

    mean_auroc = best["mean_auroc"]
    mean_ci_lower = best["mean_ci_lower"]
    mean_pval = best["mean_perm_pvalue"]

    print(f"  Best model: {best['model']}")
    print(f"  Mean AUROC = {mean_auroc:.3f} (threshold > 0.70)")
    print(f"  Mean CI lower = {mean_ci_lower:.3f} (threshold > 0.60)")
    print(f"  Mean perm p = {mean_pval:.4f} (threshold < 0.05)")

    cond1 = mean_auroc > 0.70 and mean_ci_lower > 0.60
    cond2 = mean_pval < 0.05
    cond3 = True  # SHAP check placeholder

    print(f"\n  Condition 1 (AUROC > 0.70 AND CI lower > 0.60): {'✓ PASS' if cond1 else '✗ FAIL'}")
    print(f"  Condition 2 (perm p < 0.05):                     {'✓ PASS' if cond2 else '✗ FAIL'}")
    print(f"  Condition 3 (SHAP gene check):                   ⏳ Pending (requires SHAP)")

    all_pass = cond1 and cond2

    if all_pass:
        print(f"\n  ✓ GO — Conditions 1 & 2 passed. Run SHAP to verify condition 3.")
    else:
        print(f"\n  ✗ NO-GO — Failed conditions:")
        if not cond1:
            if mean_auroc <= 0.70:
                print(f"    → Cond 1a: AUROC {mean_auroc:.3f} ≤ 0.70")
                print(f"       Action: Check mapping rate and sample QC per mission.")
            if mean_ci_lower <= 0.60:
                print(f"    → Cond 1b: CI lower {mean_ci_lower:.3f} ≤ 0.60")
        if not cond2:
            print(f"    → Cond 2: perm p {mean_pval:.4f} ≥ 0.05")
            print(f"       Action: Check for batch effects → run J3 first.")

    return all_pass


# ── CLI ────────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="Run baseline models on GeneLab_benchmark task splits"
    )
    parser.add_argument("--task", type=str, default="A1",
                        help="Task ID (default: A1)")
    parser.add_argument("--model", type=str, nargs="+",
                        choices=list(MODELS.keys()),
                        help="Models to run (default: all)")
    parser.add_argument("--quick", action="store_true",
                        help="Skip XGBoost (faster for testing)")
    parser.add_argument("--check", action="store_true",
                        help="Load existing results and run go/no-go check")
    parser.add_argument("--n-bootstrap", type=int, default=1000)
    parser.add_argument("--n-perm", type=int, default=1000)
    parser.add_argument("--quiet", action="store_true")
    parser.add_argument("--combat", action="store_true",
                        help="Use ComBat-seq corrected task splits (_lomo_combat dir)")
    return parser.parse_args()


def main():
    args = parse_args()
    task_id = args.task.upper()
    verbose = not args.quiet

    # Determine which models to run
    if args.model:
        model_names = args.model
    elif args.quick:
        model_names = ["lr", "pca_lr", "rf"]
    else:
        model_names = list(MODELS.keys())

    combat = getattr(args, "combat", False)
    suffix = "_combat" if combat else ""
    results_path = RESULTS_DIR / f"{task_id}{suffix}_baseline_results.json"

    if args.check:
        # Load existing results
        if not results_path.exists():
            print(f"[ERROR] No results found at {results_path}")
            print(f"        Run: python scripts/run_baselines.py --task {task_id}")
            return
        results = json.loads(results_path.read_text())
        check_gonogo(results, task_id)
        return

    # Run evaluation
    results = evaluate_task(
        task_id=task_id,
        model_names=model_names,
        n_bootstrap=args.n_bootstrap,
        n_perm=args.n_perm,
        verbose=verbose,
        combat=combat,
    )

    if results:
        # Save results
        RESULTS_DIR.mkdir(parents=True, exist_ok=True)
        output = {
            "task_id": task_id,
            "generated_at": datetime.now().isoformat(),
            "n_bootstrap": args.n_bootstrap,
            "n_perm": args.n_perm,
            "models": results,
        }
        results_path.write_text(json.dumps(output["models"], indent=2))
        print(f"\n  ✓ Results saved: {results_path}")

        # Auto go/no-go check
        check_gonogo(results, task_id)


if __name__ == "__main__":
    main()
