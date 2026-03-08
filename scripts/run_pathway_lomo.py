#!/usr/bin/env python3
"""
run_pathway_lomo.py — LOMO evaluation using GSVA pathway scores as features.

Replicates run_baselines.py logic but uses per-mission GSVA enrichment scores
(samples × 50 pathways) instead of gene-level features.

Usage:
    python scripts/run_pathway_lomo.py --tissue kidney --db hallmark
    python scripts/run_pathway_lomo.py --tissue eye --db hallmark
    python scripts/run_pathway_lomo.py --tissue kidney --db kegg
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT = Path(__file__).parent.parent
PATHWAY_DIR = ROOT / "processed" / "pathway_scores"
TASKS_DIR = ROOT / "tasks"
EVAL_DIR = ROOT / "evaluation"

TASK_MAP = {
    "kidney": ("A3", "A3_kidney_lomo"),
    "eye":    ("A6", "A6_eye_lomo"),
    "gastrocnemius": ("A2", "A2_gastrocnemius_lomo"),
    "thymus": ("A4", "A4_thymus_lomo"),
    "liver":  ("A1", "A1_liver_lomo"),
    "skin":   ("A5", "A5_skin_lomo"),
}

# ── Statistics (same as run_baselines.py) ──────────────────────────────────────

def bootstrap_ci(y_true, y_score, n_boot=2000, seed=42):
    rng = np.random.RandomState(seed)
    n = len(y_true)
    boot_aurocs = []
    for _ in range(n_boot):
        idx = rng.randint(0, n, n)
        if len(np.unique(y_true[idx])) < 2:
            continue
        boot_aurocs.append(roc_auc_score(y_true[idx], y_score[idx]))
    boot_aurocs = np.array(boot_aurocs)
    return float(np.percentile(boot_aurocs, 2.5)), float(np.percentile(boot_aurocs, 97.5))


def permutation_pvalue(y_true, y_score, n_perm=1000, seed=42):
    rng = np.random.RandomState(seed)
    observed = roc_auc_score(y_true, y_score)
    null_scores = []
    for _ in range(n_perm):
        y_perm = rng.permutation(y_true)
        null_scores.append(roc_auc_score(y_perm, y_score))
    null_scores = np.array(null_scores)
    # Pseudocount to avoid p=0.000 (B2 fix)
    p = float((np.sum(null_scores >= observed) + 1) / (n_perm + 1))
    return p


# ── Models ─────────────────────────────────────────────────────────────────────

def build_lr(seed=42):
    return Pipeline([
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(
            solver="saga", l1_ratio=0.5, C=1.0,
            class_weight="balanced", max_iter=5000, random_state=seed,
        ))
    ])


def build_pca_lr(n_components, seed=42):
    return Pipeline([
        ("scaler", StandardScaler()),
        ("pca", PCA(n_components=n_components, random_state=seed)),
        ("clf", LogisticRegression(
            C=1.0, class_weight="balanced", max_iter=2000, random_state=seed,
        ))
    ])


# ── Mission mapping for composite missions ────────────────────────────────────
# Skin MHU-2 is split into dorsal/femoral subsites in GSVA files
MISSION_EXPAND = {
    "skin": {
        "MHU-2": ["MHU-2_dorsal", "MHU-2_femoral"],
    },
}


# ── Data loading ───────────────────────────────────────────────────────────────

def load_pathway_features(tissue, db, missions):
    """Concatenate GSVA scores for multiple missions (samples × pathways).
    Uses intersection of pathways across missions to avoid NaN from differing
    GSVA filtration (e.g. KEGG 218-231 pathways per mission)."""
    parts = []
    expand = MISSION_EXPAND.get(tissue, {})
    for mission in missions:
        sub_missions = expand.get(mission, [mission])
        for sub in sub_missions:
            fpath = PATHWAY_DIR / tissue / f"{sub}_gsva_{db}.csv"
            if not fpath.exists():
                raise FileNotFoundError(f"Missing: {fpath}")
            df = pd.read_csv(fpath, index_col=0)
            parts.append(df)
    # Intersect columns to avoid NaN from differing pathway sets
    if len(parts) > 1:
        common_cols = set(parts[0].columns)
        for p in parts[1:]:
            common_cols &= set(p.columns)
        common_cols = sorted(common_cols)
        parts = [p[common_cols] for p in parts]
    combined = pd.concat(parts, axis=0)
    # Drop duplicate samples (e.g. thymus MHU-2 GSVA contains MHU-1 samples)
    n_dup = combined.index.duplicated().sum()
    if n_dup > 0:
        print(f"    [WARN] Dropped {n_dup} duplicate sample(s) from {tissue}/{db}")
    combined = combined[~combined.index.duplicated(keep='first')]
    return combined


def load_fold(tissue, task_dir_name, db, fold):
    """Load train/test pathway features + labels for one LOMO fold."""
    task_dir = TASKS_DIR / task_dir_name
    fold_dir = task_dir / f"fold_{fold['test_mission']}_test"

    # Labels from task fold (same as gene-level LOMO)
    train_y = pd.read_csv(fold_dir / "train_y.csv", index_col=0).iloc[:, 0]
    test_y  = pd.read_csv(fold_dir / "test_y.csv",  index_col=0).iloc[:, 0]

    # Features from pathway scores
    train_X = load_pathway_features(tissue, db, fold["train_missions"])
    test_X  = load_pathway_features(tissue, db, [fold["test_mission"]])

    # Intersect pathways between train and test
    common_pathways = sorted(set(train_X.columns) & set(test_X.columns))
    train_X = train_X[common_pathways]
    test_X  = test_X[common_pathways]

    # Align samples — handle {mission}.{sample} prefix in fold labels
    # (gastrocnemius fold labels have "RR-5.GSM..." but GSVA has "GSM...")
    def align_samples(X_df, y_series):
        if y_series.index[0] in X_df.index:
            return X_df.loc[y_series.index]
        # Try stripping "{mission}." prefix from fold labels
        stripped = y_series.index.str.replace(r'^[^.]+\.', '', regex=True)
        if stripped[0] in X_df.index:
            X_aligned = X_df.loc[stripped]
            X_aligned.index = y_series.index  # restore original names
            return X_aligned
        raise KeyError(f"Cannot align samples. GSVA: {X_df.index[:3].tolist()}, "
                       f"Fold: {y_series.index[:3].tolist()}")

    train_X = align_samples(train_X, train_y)
    test_X  = align_samples(test_X, test_y)

    return (
        train_X.values.astype(float), train_y.values.astype(float),
        test_X.values.astype(float),  test_y.values.astype(float),
    )


# ── Main ───────────────────────────────────────────────────────────────────────

def run_lomo(tissue, db, n_boot=2000, n_perm=1000):
    if tissue not in TASK_MAP:
        raise ValueError(f"Unknown tissue: {tissue}. Choose from {list(TASK_MAP)}")

    task_id, task_dir_name = TASK_MAP[tissue]
    with open(TASKS_DIR / task_dir_name / "task_info.json") as f:
        task_info = json.load(f)
    folds = task_info["folds"]

    results = {"lr": [], "pca_lr": []}

    for fold in folds:
        test_mission = fold["test_mission"]
        print(f"\n  Fold: test={test_mission} | train={fold['train_missions']}", flush=True)

        train_X, train_y, test_X, test_y = load_fold(tissue, task_dir_name, db, fold)
        n_test = len(test_y)
        n_flt  = int((test_y == 1).sum())
        n_gnd  = int((test_y == 0).sum())
        n_comp = min(50, train_X.shape[0] - 1, train_X.shape[1])  # adaptive PCA components

        for model_name, model in [("lr", build_lr()), ("pca_lr", build_pca_lr(n_comp))]:
            t0 = time.time()
            model.fit(train_X, train_y)
            y_score = model.predict_proba(test_X)[:, 1]
            t1 = time.time()

            auroc = float(roc_auc_score(test_y, y_score))
            ci_lo, ci_hi = bootstrap_ci(test_y, y_score, n_boot=n_boot)
            perm_p = permutation_pvalue(test_y, y_score, n_perm=n_perm)

            print(f"    {model_name:8s} AUROC={auroc:.4f} CI=[{ci_lo:.3f},{ci_hi:.3f}] perm_p={perm_p:.3f} "
                  f"(n={n_test}: {n_flt}F+{n_gnd}G, t={t1-t0:.1f}s)", flush=True)

            results[model_name].append({
                "test_mission": test_mission,
                "train_missions": fold["train_missions"],
                "auroc": auroc,
                "ci_lower": ci_lo,
                "ci_upper": ci_hi,
                "perm_p": perm_p,
                "n_test": n_test,
                "n_flight_test": n_flt,
                "n_ground_test": n_gnd,
                "n_pathways": train_X.shape[1],
                "n_train": train_X.shape[0],
                "train_time_s": round(t1 - t0, 2),
            })

    # Summary
    summary = {}
    for model_name, fold_results in results.items():
        aurocs = [r["auroc"] for r in fold_results]
        ci_lowers = [r["ci_lower"] for r in fold_results]
        perm_ps = [r["perm_p"] for r in fold_results]
        mean_auroc = float(np.mean(aurocs))
        mean_ci_lower = float(np.mean(ci_lowers))
        mean_perm_p = float(np.mean(perm_ps))

        go_auroc = mean_auroc > 0.700
        go_ci    = mean_ci_lower > 0.500
        go_perm  = mean_perm_p < 0.050
        go       = go_auroc and go_ci and go_perm

        summary[model_name] = {
            "n_folds": len(fold_results),
            "mean_auroc": mean_auroc,
            "std_auroc": float(np.std(aurocs)),
            "mean_ci_lower": mean_ci_lower,
            "mean_perm_p": mean_perm_p,
            "go_auroc": go_auroc,
            "go_ci": go_ci,
            "go_perm": go_perm,
            "go": go,
            "folds": fold_results,
        }

        verdict = "✓ GO" if go else "✗ NO-GO"
        print(f"\n  [{model_name}] Mean AUROC={mean_auroc:.4f} ± {np.std(aurocs):.4f} | "
              f"CI lower={mean_ci_lower:.3f} | perm_p={mean_perm_p:.4f} → {verdict}")

    return summary


def main():
    parser = argparse.ArgumentParser(description="Pathway LOMO evaluation (GSVA features)")
    parser.add_argument("--tissue", required=True, choices=list(TASK_MAP), help="Tissue name")
    parser.add_argument("--db", default="hallmark",
                        choices=["hallmark", "kegg", "reactome", "mitocarta"],
                        help="Pathway database (default: hallmark)")
    parser.add_argument("--n-boot", type=int, default=2000, help="Bootstrap iterations (default: 2000)")
    parser.add_argument("--n-perm", type=int, default=1000, help="Permutation iterations (default: 1000)")
    args = parser.parse_args()

    print(f"\n=== Pathway LOMO: {args.tissue} ({args.db}) ===")
    print(f"n_boot={args.n_boot}, n_perm={args.n_perm}")

    summary = run_lomo(args.tissue, args.db, args.n_boot, args.n_perm)

    # Save output
    task_id = TASK_MAP[args.tissue][0]
    out_path = EVAL_DIR / f"{task_id}_pathway_{args.db}_results.json"
    EVAL_DIR.mkdir(exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nSaved → {out_path}")


if __name__ == "__main__":
    main()
