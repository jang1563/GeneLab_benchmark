#!/usr/bin/env python3
"""
generate_submission.py — Generate PCA-LR baseline submission JSONs.

Produces per-sample prediction files in the standard submission format used by
evaluate_submission.py.  Replicates the build_pca_lr() configuration from
run_baselines.py so results are directly comparable to existing baseline files.

A-tasks (LOMO):
  Loads each fold_*_test directory, trains PCA-LR on train split, saves
  predicted Flight probabilities for test samples.

B-tasks (cross-mission):
  Loads the tissue-wide log2-normalised matrix from processed/A_detection/,
  applies the same variance filter as cross_mission_transfer.py (DD-03),
  trains PCA-LR for each directed pair, saves predicted probabilities.

Output:
  evaluation/submission_PCALR_baseline_{task_id}.json

Usage:
  python scripts/generate_submission.py --task A5
  python scripts/generate_submission.py --task A1 --task-dir A1_liver_lomo
  python scripts/generate_submission.py --task A4 --fold RR-9
  python scripts/generate_submission.py --task A5 A6 B5 B6
  python scripts/generate_submission.py --all
"""
import json
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

BASE_DIR    = Path(__file__).resolve().parent.parent
TASKS_DIR   = BASE_DIR / "tasks"
PROC_DIR    = BASE_DIR / "processed" / "A_detection"
RESULTS_DIR = BASE_DIR / "evaluation"

# DD-03: top 75th-percentile variance filter (train-only)
VARIANCE_PERCENTILE = 0.25  # discard bottom 25%
FLIGHT_LABEL  = "Flight"
GROUND_LABELS = {"GC", "VC"}

A_TASKS = {
    "A1": "liver_lomo",
    "A2": "gastrocnemius_lomo",
    "A4": "thymus_lomo",
    "A5": "skin_lomo",
    "A6": "eye_lomo",
}
B_TASKS = {
    "B2": "gastrocnemius",
    "B4": "thymus",
    "B5": "skin",
    "B6": "eye",
}
ALL_TASKS = list(A_TASKS) + list(B_TASKS)


def resolve_a_task_dir(task_id: str, tissue_suffix: str,
                       task_dir_name: str | None = None) -> Path:
    """
    Resolve A-task directory deterministically.
    If task_dir_name is given, use it with task_id consistency checks.
    """
    if task_dir_name:
        task_dir = TASKS_DIR / task_dir_name
        if not task_dir.exists() or not task_dir.is_dir():
            raise FileNotFoundError(f"Task directory not found: {task_dir}")
        if not task_dir.name.startswith(f"{task_id}_"):
            raise ValueError(
                f"--task-dir '{task_dir.name}' does not match task_id '{task_id}'"
            )
        return task_dir

    expected = TASKS_DIR / f"{task_id}_{tissue_suffix}"
    if expected.exists() and expected.is_dir():
        return expected

    candidates = sorted(
        d for d in TASKS_DIR.iterdir()
        if d.is_dir() and d.name.startswith(f"{task_id}_")
    )
    if not candidates:
        raise FileNotFoundError(
            f"No task directory found for {task_id}. Use --task-dir explicitly."
        )
    if len(candidates) > 1:
        names = ", ".join(d.name for d in candidates)
        raise ValueError(
            f"Ambiguous task_id '{task_id}': {names}. Use --task-dir to select one."
        )
    return candidates[0]


def resolve_a_fold_dirs(task_dir: Path, fold_name: str | None = None) -> list[Path]:
    """
    Resolve A-task folds with exact matching only.
    Accepts:
      --fold RR-9           -> fold_RR-9_test
      --fold fold_RR-9_test -> fold_RR-9_test
    """
    if fold_name:
        if fold_name.startswith("fold_"):
            candidate_names = [fold_name]
            if not fold_name.endswith("_test"):
                candidate_names.append(f"{fold_name}_test")
        else:
            candidate_names = [f"fold_{fold_name}_test"]

        matches = [task_dir / name for name in candidate_names
                   if (task_dir / name).exists() and (task_dir / name).is_dir()]
        if not matches:
            available = sorted(
                d.name for d in task_dir.iterdir()
                if d.is_dir() and d.name.startswith("fold_") and d.name.endswith("_test")
            )
            raise FileNotFoundError(
                f"No fold found matching '{fold_name}' in {task_dir}. "
                f"Available test folds: {available}"
            )
        fold_dir = matches[0]
        if not fold_dir.name.endswith("_test"):
            raise ValueError(
                f"Fold '{fold_dir.name}' is not a *_test fold. "
                "Submission generation only supports fold_*_test."
            )
        return [fold_dir]

    fold_dirs = sorted([
        d for d in task_dir.iterdir()
        if d.is_dir() and d.name.startswith("fold_") and d.name.endswith("_test")
    ])
    if not fold_dirs:
        raise FileNotFoundError(f"No fold_*_test directories in {task_dir}")
    return fold_dirs


def build_pca_lr(n_train: int) -> Pipeline:
    """PCA-LR — matches run_baselines.py build_pca_lr() exactly."""
    n_comps = min(50, n_train - 1)
    return Pipeline([
        ("scaler", StandardScaler()),
        ("pca",    PCA(n_components=n_comps, random_state=42)),
        ("clf",    LogisticRegression(
            C=1.0,
            class_weight="balanced",
            max_iter=5000,
            random_state=42,
        )),
    ])


def generate_a_submission(task_id: str, tissue_suffix: str,
                          task_dir_name: str | None = None,
                          fold_name: str | None = None) -> dict:
    """
    LOMO: train PCA-LR on each fold_*_test, return {fold_key: {sample: prob}}.
    Only _test folds are included (holdout folds without test_y are skipped).
    """
    task_dir = resolve_a_task_dir(task_id, tissue_suffix, task_dir_name=task_dir_name)
    fold_dirs = resolve_a_fold_dirs(task_dir, fold_name=fold_name)

    print(f"  Task directory: {task_dir.name}")
    if fold_name:
        print(f"  Fold filter: {fold_name} (exact match)")

    predictions = {}
    for fold_dir in fold_dirs:
        train_X = pd.read_csv(fold_dir / "train_X.csv", index_col=0)
        test_X  = pd.read_csv(fold_dir / "test_X.csv",  index_col=0)
        train_y = pd.read_csv(fold_dir / "train_y.csv", index_col=0).squeeze()

        # Align columns (should already match, but be safe)
        common = train_X.columns.intersection(test_X.columns)
        train_X = train_X[common]
        test_X  = test_X[common]

        model = build_pca_lr(len(train_X))
        model.fit(train_X.values.astype(np.float32),
                  train_y.values.astype(int))
        probs = model.predict_proba(test_X.values.astype(np.float32))[:, 1]

        fold_key = fold_dir.name  # e.g. "fold_RR-6_test"
        predictions[fold_key] = {
            sid: round(float(p), 8)
            for sid, p in zip(test_X.index, probs)
        }
        n_flt = int(train_y.sum())
        n_gnd = int((train_y == 0).sum())
        print(f"    {fold_key}: n_train={len(train_X)} ({n_flt}F+{n_gnd}G), "
              f"n_test={len(test_X)}, n_features={len(common)}")

    return predictions


def generate_b_submission(tissue: str) -> dict:
    """
    Cross-mission: train PCA-LR for each directed pair, return predictions.
    Applies DD-03 variance filter on train set only.
    """
    expr_path = PROC_DIR / tissue / f"{tissue}_all_missions_log2_norm.csv"
    meta_path = PROC_DIR / tissue / f"{tissue}_all_missions_metadata.csv"

    if not expr_path.exists():
        raise FileNotFoundError(f"Expression matrix not found: {expr_path}")

    expr = pd.read_csv(expr_path, index_col=0)
    meta = pd.read_csv(meta_path, index_col=0)

    # Drop non-numeric expression columns if any leaked in
    non_gene = [c for c in expr.columns if c in {"mission", "osd_id", "label"}]
    if non_gene:
        expr = expr.drop(columns=non_gene)
    expr = expr.select_dtypes(include=[np.number])

    # Binary labels: Flight=1, GC/VC=0, AG/BC excluded
    binary = pd.Series(np.nan, index=meta.index)
    binary[meta["label"] == FLIGHT_LABEL] = 1
    binary[meta["label"].isin(GROUND_LABELS)] = 0

    valid = ~binary.isna()
    expr   = expr[valid]
    meta   = meta[valid]
    binary = binary[valid]

    missions = sorted(meta["mission"].unique())
    predictions = {}

    for train_m in missions:
        for test_m in missions:
            if train_m == test_m:
                continue
            train_mask = meta["mission"] == train_m
            test_mask  = meta["mission"] == test_m

            X_tr = expr[train_mask]
            X_te = expr[test_mask]
            y_tr = binary[train_mask]

            n_flt = int((y_tr == 1).sum())
            n_gnd = int((y_tr == 0).sum())
            if n_flt < 3 or n_gnd < 3:
                print(f"    [SKIP] pair_{train_m}_{test_m}: "
                      f"n_flt={n_flt}, n_gnd={n_gnd} (< 3)")
                continue

            # DD-03 variance filter (train-only)
            gene_var = X_tr.var(axis=0)
            threshold = gene_var.quantile(VARIANCE_PERCENTILE)
            genes = gene_var[gene_var >= threshold].index
            X_tr = X_tr[genes]
            X_te = X_te[genes]

            model = build_pca_lr(len(X_tr))
            model.fit(X_tr.values, y_tr.values.astype(int))
            probs = model.predict_proba(X_te.values)[:, 1]

            pair_key = f"pair_{train_m}_{test_m}"
            predictions[pair_key] = {
                sid: float(p)
                for sid, p in zip(X_te.index, probs)
            }
            print(f"    {pair_key}: n_train={len(X_tr)} ({n_flt}F+{n_gnd}G), "
                  f"n_test={len(X_te)}, n_genes={len(genes)}")

    return predictions


def save_submission(task_id: str, predictions: dict,
                    tissue_label: str) -> Path:
    """Save evaluation/submission_PCALR_baseline_{task_id}.json."""
    is_b = task_id.startswith("B")
    description = (
        f"StandardScaler + PCA(min(50,n_train-1)) + LogReg(L2, C=1.0, balanced), "
        + (
            f"PCA-LR cross-mission transfer, {tissue_label}"
            if is_b else
            f"LOMO Tier 1 baseline — {tissue_label}"
        )
    )
    out = {
        "task_id": task_id,
        "model_name": "PCA-LR_baseline_v1",
        "model_description": description,
        "tier": "1",
        "submission_date": datetime.now().strftime("%Y-%m-%d"),
        "predictions": predictions,
    }
    path = RESULTS_DIR / f"submission_PCALR_baseline_{task_id}.json"
    path.write_text(json.dumps(out, indent=2))
    return path


def main():
    parser = argparse.ArgumentParser(
        description="Generate PCA-LR baseline submission JSONs"
    )
    parser.add_argument(
        "--all", action="store_true",
        help="Generate submissions for all supported tasks"
    )
    parser.add_argument(
        "--task", nargs="+", default=["all"],
        help=f"Task IDs. Use 'all' or list. Available: {ALL_TASKS}",
    )
    parser.add_argument(
        "--task-dir", default=None,
        help="Explicit task directory under tasks/ (Category A only; single task mode)"
    )
    parser.add_argument(
        "--fold", default=None,
        help="Exact A-task fold selector: RR-9 or fold_RR-9_test (single task mode)"
    )
    args = parser.parse_args()

    if args.all or "all" in args.task:
        tasks = ALL_TASKS
    else:
        tasks = []
        for t in args.task:
            t_up = t.upper()
            if t_up in ALL_TASKS:
                tasks.append(t_up)
            else:
                print(f"[WARN] Unknown task: {t}. Available: {ALL_TASKS}")
        if not tasks:
            print("No valid tasks. Exiting.")
            return

    if args.task_dir is not None:
        if len(tasks) != 1 or tasks[0] not in A_TASKS:
            parser.error("--task-dir can only be used with a single Category A task")
    if args.fold is not None:
        if len(tasks) != 1 or tasks[0] not in A_TASKS:
            parser.error("--fold can only be used with a single Category A task")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    for task_id in tasks:
        print(f"\n{'='*60}")
        print(f"Generating: {task_id}")
        print(f"{'='*60}")

        try:
            if task_id in A_TASKS:
                tissue_suffix = A_TASKS[task_id]
                preds = generate_a_submission(
                    task_id,
                    tissue_suffix,
                    task_dir_name=args.task_dir if len(tasks) == 1 else None,
                    fold_name=args.fold if len(tasks) == 1 else None,
                )
                tissue_label = tissue_suffix.replace("_lomo", "")
            else:
                tissue = B_TASKS[task_id]
                preds = generate_b_submission(tissue)
                tissue_label = tissue

            if not preds:
                print(f"  [WARN] No predictions generated for {task_id}")
                continue

            out_path = save_submission(task_id, preds, tissue_label)
            print(f"\n  Saved: {out_path.name}  "
                  f"({len(preds)} {'folds' if task_id in A_TASKS else 'pairs'})")

        except Exception as e:
            print(f"  [ERROR] {task_id}: {e}")
            import traceback
            traceback.print_exc()

    print("\nDone.")


if __name__ == "__main__":
    main()
