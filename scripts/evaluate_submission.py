#!/usr/bin/env python3
"""
evaluate_submission.py — GeneLab_benchmark: Submission Evaluator

Evaluates an external model submission against ground truth labels.
Computes AUROC, bootstrap CI, and permutation p-value per fold and overall.

Submission format: JSON (see docs/submission_format.md)

Usage:
  python evaluate_submission.py --submission path/to/submission.json --task A4
  python evaluate_submission.py --submission path/to/submission.json --task A4 --validate-only
  python evaluate_submission.py --submission path/to/submission.json --task A4 --output results.json
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

# ── Paths ─────────────────────────────────────────────────────────────────────
REPO_ROOT = Path(__file__).resolve().parent.parent
TASKS_DIR = REPO_ROOT / "tasks"

# ── Go/No-Go Thresholds (DD-11) ───────────────────────────────────────────────
AUROC_THRESHOLD = 0.700
CI_LOWER_THRESHOLD = 0.500
PERM_P_THRESHOLD = 0.050

# ── Evaluation Parameters ─────────────────────────────────────────────────────
N_BOOTSTRAP = 2000
N_PERM = 1000
RANDOM_SEED = 42


# ── Helper Functions ──────────────────────────────────────────────────────────

def bootstrap_auroc_ci(y_true: np.ndarray, y_score: np.ndarray,
                        n_bootstrap: int = N_BOOTSTRAP,
                        seed: int = RANDOM_SEED):
    """Bootstrap 95% CI for AUROC."""
    if n_bootstrap <= 0:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    n = len(y_true)
    boot_scores = []
    for _ in range(n_bootstrap):
        idx = rng.integers(0, n, n)
        if len(np.unique(y_true[idx])) < 2:
            continue
        boot_scores.append(roc_auc_score(y_true[idx], y_score[idx]))
    if not boot_scores:
        return float("nan"), float("nan")
    ci_lo = float(np.percentile(boot_scores, 2.5))
    ci_hi = float(np.percentile(boot_scores, 97.5))
    return ci_lo, ci_hi


def permutation_pvalue(y_true: np.ndarray, y_score: np.ndarray,
                        n_perm: int = N_PERM,
                        seed: int = RANDOM_SEED):
    """Permutation p-value for AUROC (H0: AUROC = 0.5). Pseudocount applied."""
    rng = np.random.default_rng(seed)
    observed = roc_auc_score(y_true, y_score)
    null_scores = []
    for _ in range(n_perm):
        y_perm = rng.permutation(y_true)
        null_scores.append(roc_auc_score(y_perm, y_score))
    null_arr = np.array(null_scores)
    p = float((np.sum(null_arr >= observed) + 1) / (n_perm + 1))
    return p


# ── Validation ────────────────────────────────────────────────────────────────

def get_task_dir(task_id: str, task_dir_name=None) -> Path:
    """Return the task directory for a given task ID with deterministic resolution."""
    if task_dir_name:
        task_dir = TASKS_DIR / task_dir_name
        if not task_dir.exists() or not task_dir.is_dir():
            raise ValueError(f"Task directory not found: {task_dir}")
        if not task_dir.name.startswith(f"{task_id}_"):
            raise ValueError(
                f"--task-dir '{task_dir.name}' does not match task_id '{task_id}'"
            )
        return task_dir

    candidates = sorted(
        t for t in TASKS_DIR.iterdir()
        if t.is_dir() and t.name.startswith(f"{task_id}_")
    )
    if not candidates:
        raise ValueError(f"Task '{task_id}' not found in {TASKS_DIR}")
    if len(candidates) > 1:
        names = ", ".join(t.name for t in candidates)
        raise ValueError(
            f"Ambiguous task_id '{task_id}': {names}. Use --task-dir to select one."
        )
    return candidates[0]


def validate_submission(submission: dict, task_id: str,
                        task_dir_name=None) -> list[str]:
    """
    Validate submission structure. Returns list of error messages (empty = valid).
    """
    errors = []

    # Required top-level fields
    for field in ["task_id", "model_name", "predictions"]:
        if field not in submission:
            errors.append(f"Missing required field: '{field}'")

    if errors:
        return errors

    # task_id match
    if submission["task_id"] != task_id:
        errors.append(
            f"task_id mismatch: submission has '{submission['task_id']}' but evaluating '{task_id}'"
        )

    # model_name format
    model_name = submission.get("model_name", "")
    if not model_name or len(model_name) > 50:
        errors.append(f"model_name must be non-empty and ≤50 chars (got '{model_name}')")
    import re
    if not re.match(r'^[\w\-\.]+$', model_name):
        errors.append(f"model_name must be alphanumeric/underscore/hyphen/dot only: '{model_name}'")

    # Check predictions structure
    if not isinstance(submission.get("predictions"), dict):
        errors.append("'predictions' must be a dict")
        return errors

    # Load expected folds from task
    try:
        task_dir = get_task_dir(task_id, task_dir_name=task_dir_name)
    except ValueError as e:
        errors.append(str(e))
        return errors

    # Support both Category A (fold_ prefix) and Category B (pair_ prefix)
    fold_prefix = "pair_" if task_id.startswith("B") else "fold_"
    all_folds = sorted([d.name for d in task_dir.iterdir()
                        if d.is_dir() and d.name.startswith(fold_prefix)])
    # Holdout folds are optional (labels are private; Tier 1 = LOMO only)
    expected_folds = [f for f in all_folds if "holdout" not in f]
    optional_folds = [f for f in all_folds if "holdout" in f]

    submitted_folds = sorted(submission["predictions"].keys())

    # Check all LOMO folds present
    for fold in expected_folds:
        if fold not in submission["predictions"]:
            errors.append(f"Missing fold in predictions: '{fold}'")

    for fold in submitted_folds:
        if fold not in expected_folds and fold not in optional_folds:
            errors.append(f"Unknown fold in submission: '{fold}' (expected: {expected_folds})")

    # Check sample IDs and probability values per fold (LOMO + optional holdout if submitted)
    folds_to_check = expected_folds + [f for f in optional_folds if f in submission["predictions"]]
    for fold in folds_to_check:
        if fold not in submission["predictions"]:
            continue
        fold_preds = submission["predictions"][fold]
        if not isinstance(fold_preds, dict):
            errors.append(f"[{fold}] predictions must be a dict of sample_id -> probability")
            continue
        fold_dir = task_dir / fold
        test_y_path = fold_dir / "test_y.csv"
        if not test_y_path.exists():
            errors.append(f"Cannot find {test_y_path}")
            continue

        y_df = pd.read_csv(test_y_path, index_col=0)
        expected_samples = set(y_df.index.astype(str))
        submitted_samples = set(fold_preds.keys())

        missing = expected_samples - submitted_samples
        extra = submitted_samples - expected_samples

        if missing:
            errors.append(f"[{fold}] Missing predictions for {len(missing)} samples: {sorted(missing)[:5]}...")
        if extra:
            errors.append(f"[{fold}] Unknown sample IDs: {sorted(extra)[:5]}...")

        # Check probability values
        for sample_id, prob in fold_preds.items():
            if not isinstance(prob, (int, float)):
                errors.append(f"[{fold}] Non-numeric probability for '{sample_id}': {prob}")
                break
            if not (0.0 <= float(prob) <= 1.0):
                errors.append(f"[{fold}] Probability out of [0,1] for '{sample_id}': {prob}")
                break

    return errors


# ── Evaluation ────────────────────────────────────────────────────────────────

def evaluate_fold(fold_dir: Path, fold_predictions: dict,
                  n_bootstrap: int = N_BOOTSTRAP,
                  n_perm: int = N_PERM) -> dict:
    """Evaluate one fold. Returns per-fold metrics dict."""
    y_df = pd.read_csv(fold_dir / "test_y.csv", index_col=0)
    y_true = y_df.iloc[:, 0].values.astype(float)
    sample_ids = y_df.index.astype(str).tolist()

    y_score = np.array([float(fold_predictions[sid]) for sid in sample_ids])

    if len(np.unique(y_true)) < 2:
        return {
            "fold": fold_dir.name,
            "n_test": len(y_true),
            "n_flight": int((y_true == 1).sum()),
            "n_ground": int((y_true == 0).sum()),
            "auroc": None,
            "ci_lower": None,
            "ci_upper": None,
            "perm_p": None,
            "error": "Only one class in test set",
        }

    auroc = float(roc_auc_score(y_true, y_score))
    ci_lo, ci_hi = bootstrap_auroc_ci(y_true, y_score, n_bootstrap=n_bootstrap)
    perm_p = permutation_pvalue(y_true, y_score, n_perm=n_perm)

    return {
        "fold": fold_dir.name,
        "n_test": len(y_true),
        "n_flight": int((y_true == 1).sum()),
        "n_ground": int((y_true == 0).sum()),
        "auroc": auroc,
        "ci_lower": ci_lo,
        "ci_upper": ci_hi,
        "perm_p": perm_p,
    }


def evaluate_submission_full(submission: dict, task_id: str,
                              verbose: bool = True,
                              n_bootstrap: int = N_BOOTSTRAP,
                              n_perm: int = N_PERM,
                              task_dir_name=None) -> dict:
    """Run full evaluation. Returns results dict."""
    task_dir = get_task_dir(task_id, task_dir_name=task_dir_name)
    # Support Category A (fold_) and Category B (pair_); holdout folds excluded
    fold_prefix = "pair_" if task_id.startswith("B") else "fold_"
    fold_dirs = sorted([d for d in task_dir.iterdir()
                        if d.is_dir() and d.name.startswith(fold_prefix)
                        and "holdout" not in d.name])

    fold_results = []
    for fold_dir in fold_dirs:
        fold_name = fold_dir.name
        fold_preds = submission["predictions"].get(fold_name, {})
        result = evaluate_fold(fold_dir, fold_preds,
                               n_bootstrap=n_bootstrap, n_perm=n_perm)
        fold_results.append(result)
        if verbose:
            auroc_str = f"{result['auroc']:.3f}" if result["auroc"] is not None else "N/A"
            ci_str = (f"[{result['ci_lower']:.3f}, {result['ci_upper']:.3f}]"
                      if result["ci_lower"] is not None else "[N/A]")
            perm_str = f"{result['perm_p']:.3f}" if result["perm_p"] is not None else "N/A"
            print(f"  {fold_name}: AUROC={auroc_str} {ci_str} perm_p={perm_str} "
                  f"(n={result['n_test']}: {result['n_flight']}F + {result['n_ground']}G)")

    # Aggregate across folds
    valid_aurocs = [r["auroc"] for r in fold_results if r["auroc"] is not None]
    valid_perm_ps = [r["perm_p"] for r in fold_results if r["perm_p"] is not None]

    if valid_aurocs:
        mean_auroc = float(np.mean(valid_aurocs))
        std_auroc = float(np.std(valid_aurocs))
        mean_perm_p = float(np.mean(valid_perm_ps))

        # Pool all test samples for overall CI
        all_y_true, all_y_score = [], []
        for fold_dir in fold_dirs:
            y_df = pd.read_csv(fold_dir / "test_y.csv", index_col=0)
            y_true = y_df.iloc[:, 0].values.astype(float)
            sample_ids = y_df.index.astype(str).tolist()
            fold_name = fold_dir.name
            fold_preds = submission["predictions"].get(fold_name, {})
            y_score = np.array([float(fold_preds.get(sid, 0.5)) for sid in sample_ids])
            all_y_true.extend(y_true.tolist())
            all_y_score.extend(y_score.tolist())

        all_y_true = np.array(all_y_true)
        all_y_score = np.array(all_y_score)
        overall_auroc = float(roc_auc_score(all_y_true, all_y_score))
        overall_ci_lo, overall_ci_hi = bootstrap_auroc_ci(
            all_y_true, all_y_score, n_bootstrap=n_bootstrap
        )

        # Go/No-Go (based on fold mean CI lower from per-fold CI means)
        mean_ci_lower = float(np.mean([r["ci_lower"] for r in fold_results
                                        if r["ci_lower"] is not None]))
        go_nogo = (mean_auroc > AUROC_THRESHOLD
                   and mean_ci_lower > CI_LOWER_THRESHOLD
                   and mean_perm_p < PERM_P_THRESHOLD)
    else:
        mean_auroc = std_auroc = mean_perm_p = None
        overall_auroc = overall_ci_lo = overall_ci_hi = None
        mean_ci_lower = None
        go_nogo = False

    results = {
        "task_id": task_id,
        "model_name": submission["model_name"],
        "model_description": submission.get("model_description", ""),
        "submission_date": submission.get("submission_date", ""),
        "n_folds": len(fold_results),
        "fold_results": fold_results,
        "summary": {
            "mean_auroc": mean_auroc,
            "std_auroc": std_auroc,
            "mean_ci_lower": mean_ci_lower,
            "overall_auroc_pooled": overall_auroc,
            "overall_ci_lower": overall_ci_lo,
            "overall_ci_upper": overall_ci_hi,
            "mean_perm_p": mean_perm_p,
            "go_nogo": go_nogo,
        },
        "thresholds": {
            "auroc": AUROC_THRESHOLD,
            "ci_lower": CI_LOWER_THRESHOLD,
            "perm_p": PERM_P_THRESHOLD,
        },
    }
    return results


# ── CLI ───────────────────────────────────────────────────────────────────────

def print_summary(results: dict) -> None:
    s = results["summary"]
    task_id = results["task_id"]
    is_cross_mission = task_id.startswith("B")

    print()
    print("=" * 60)
    print(f"  Task:   {task_id}")
    print(f"  Model:  {results['model_name']}")
    label = "Pairs" if is_cross_mission else "Folds"
    print(f"  {label}: {results['n_folds']}")
    print("-" * 60)

    if s["mean_auroc"] is not None:
        print(f"  Mean AUROC:       {s['mean_auroc']:.4f} ± {s['std_auroc']:.4f}")
        print(f"  Mean CI lower:    {s['mean_ci_lower']:.4f}")
        print(f"  Pooled AUROC:     {s['overall_auroc_pooled']:.4f} "
              f"[{s['overall_ci_lower']:.4f}, {s['overall_ci_upper']:.4f}]")
        print(f"  Mean perm_p:      {s['mean_perm_p']:.4f}")
        print("-" * 60)

        if is_cross_mission:
            # Category B: no single GO/NO-GO; show transfer pattern summary (DD-17)
            fold_results = results.get("fold_results", [])
            large_pairs = [r for r in fold_results if r["n_test"] >= 10]
            sig_pairs   = [r for r in fold_results if r.get("perm_p") is not None
                           and r["perm_p"] < 0.05]
            hi_auroc    = [r for r in fold_results if r.get("auroc") is not None
                           and r["auroc"] >= 0.70]
            print(f"  Transfer Pattern Summary (DD-17):")
            print(f"    AUROC ≥ 0.70 pairs: {len(hi_auroc)}/{results['n_folds']}")
            print(f"    perm_p < 0.05 pairs: {len(sig_pairs)}/{results['n_folds']}")
            if large_pairs:
                large_aurocs = [r["auroc"] for r in large_pairs if r["auroc"] is not None]
                large_sig    = [r for r in large_pairs if r.get("perm_p") is not None
                                and r["perm_p"] < 0.05]
                print(f"    Large pairs (n≥10): {len(large_pairs)}  "
                      f"mean AUROC={np.mean(large_aurocs):.3f}  "
                      f"sig={len(large_sig)}/{len(large_pairs)}")
            print(f"  (Category B: no single GO/NO-GO — use per-pair table above)")
        else:
            status = "✓ GO" if s["go_nogo"] else "✗ NO-GO"
            print(f"  Go/No-Go:         {status}")
            if not s["go_nogo"]:
                if s["mean_auroc"] <= AUROC_THRESHOLD:
                    print(f"    ✗ Mean AUROC {s['mean_auroc']:.3f} ≤ threshold {AUROC_THRESHOLD}")
                if s["mean_ci_lower"] <= CI_LOWER_THRESHOLD:
                    print(f"    ✗ CI lower {s['mean_ci_lower']:.3f} ≤ threshold {CI_LOWER_THRESHOLD}")
                if s["mean_perm_p"] >= PERM_P_THRESHOLD:
                    print(f"    ✗ perm_p {s['mean_perm_p']:.3f} ≥ threshold {PERM_P_THRESHOLD}")
    else:
        print(f"  No valid {'pairs' if is_cross_mission else 'folds'} to aggregate.")

    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate a GeneLab benchmark submission."
    )
    parser.add_argument("--submission", required=True,
                        help="Path to submission JSON file")
    parser.add_argument("--task", required=True,
                        help="Task ID (e.g., A2, A4)")
    parser.add_argument("--task-dir", default=None,
                        help="Explicit task directory name under tasks/ (recommended for A1 variants)")
    parser.add_argument("--validate-only", action="store_true",
                        help="Only validate submission format, do not evaluate")
    parser.add_argument("--output", default=None,
                        help="Save evaluation results to JSON file")
    parser.add_argument("--n-bootstrap", type=int, default=N_BOOTSTRAP,
                        help=f"Bootstrap iterations (default: {N_BOOTSTRAP})")
    parser.add_argument("--n-perm", type=int, default=N_PERM,
                        help=f"Permutation iterations (default: {N_PERM})")
    args = parser.parse_args()
    if args.n_bootstrap < 0:
        print("Error: --n-bootstrap must be >= 0", file=sys.stderr)
        sys.exit(1)
    if args.n_perm < 0:
        print("Error: --n-perm must be >= 0", file=sys.stderr)
        sys.exit(1)

    # Load submission
    sub_path = Path(args.submission)
    if not sub_path.exists():
        print(f"Error: submission file not found: {sub_path}", file=sys.stderr)
        sys.exit(1)

    with open(sub_path) as f:
        submission = json.load(f)

    print(f"\nGeneLab Benchmark Evaluator")
    print(f"Task: {args.task} | Submission: {sub_path.name}")
    print()

    # Validate
    print("Validating submission format...")
    errors = validate_submission(submission, args.task, task_dir_name=args.task_dir)
    if errors:
        print(f"  ✗ Validation FAILED ({len(errors)} errors):")
        for e in errors:
            print(f"    - {e}")
        sys.exit(1)
    else:
        print("  ✓ Validation passed")

    if args.validate_only:
        return

    # Evaluate
    n_bootstrap = args.n_bootstrap
    n_perm = args.n_perm

    print(f"\nEvaluating {len(submission['predictions'])} fold(s)...")
    results = evaluate_submission_full(submission, args.task, verbose=True,
                                       n_bootstrap=n_bootstrap, n_perm=n_perm,
                                       task_dir_name=args.task_dir)

    print_summary(results)

    # Save
    if args.output:
        out_path = Path(args.output)
        with open(out_path, "w") as f:
            json.dump(results, f, indent=2)
        print(f"\nResults saved to: {out_path}")
    else:
        # Default output path
        out_name = (f"evaluation/submission_{submission['model_name']}"
                    f"_{args.task}_eval.json")
        out_path = REPO_ROOT / out_name
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w") as f:
            json.dump(results, f, indent=2)
        print(f"\nResults saved to: {out_path}")


if __name__ == "__main__":
    main()
