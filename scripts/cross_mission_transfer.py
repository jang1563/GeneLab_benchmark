#!/usr/bin/env python3
"""
cross_mission_transfer.py — GeneLab_benchmark: Category B Cross-Mission Transfer

Creates N×N pairwise transfer AUROC matrix for each tissue.
Tests H1: "Which tissue shows the most consistent cross-mission signatures?"

Two transfer methods:
  1. ML-based: Train classifier on mission i → test on mission j → AUROC
  2. LFC-signature: Compute LFC from mission i → project mission j → AUROC

Statistical features:
  - Bootstrap 95% CI for each pairwise AUROC
  - Bootstrap 95% CI for per-tissue mean AUROC
  - Permutation test for H1 cross-tissue comparison
  - Sensitivity analysis: liver with/without MHU-2 (small-sample artifact)

Output:
  processed/B_cross_mission/{tissue}/
    B_transfer_matrix_{method}.csv      # N×N AUROC
    B_transfer_ci_low_{method}.csv      # N×N CI lower bound
    B_transfer_ci_high_{method}.csv     # N×N CI upper bound
    B_transfer_details.json             # per-pair details with CIs
  evaluation/
    B_cross_mission_summary.json        # all tissues summary with H1 stats

Usage:
  python scripts/cross_mission_transfer.py --tissue liver
  python scripts/cross_mission_transfer.py --all
  python scripts/cross_mission_transfer.py --tissue liver --method lfc
  python scripts/cross_mission_transfer.py --all --compare-tissues
  python scripts/cross_mission_transfer.py --all --no-bootstrap   # fast mode
"""

import json
import argparse
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from itertools import combinations

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent.parent
PROCESSED_DIR = BASE_DIR / "processed" / "A_detection"
B_OUTPUT_DIR = BASE_DIR / "processed" / "B_cross_mission"
RESULTS_DIR = BASE_DIR / "evaluation"
TASKS_DIR = BASE_DIR / "tasks"

# ── Config ─────────────────────────────────────────────────────────────────────
VARIANCE_PERCENTILE = 0.25   # keep top 75% by variance (train-only)
FLIGHT_LABEL = "Flight"
GROUND_LABELS = {"GC", "VC"}
MIN_SAMPLES_PER_CLASS = 3    # minimum samples per class to attempt ML
N_BOOTSTRAP = 2000           # bootstrap resamples for CI
CI_ALPHA = 0.05              # 95% confidence interval
N_PERMUTATIONS = 10000       # permutation test iterations

# Tissues and their available missions
TISSUE_CONFIG = {
    "liver":         {"task_id": "B1", "min_missions": 3},
    "gastrocnemius": {"task_id": "B2", "min_missions": 3},
    "kidney":        {"task_id": "B3", "min_missions": 3},
    "thymus":        {"task_id": "B4", "min_missions": 3},
    "skin":          {"task_id": "B5", "min_missions": 3},
    "eye":           {"task_id": "B6", "min_missions": 3},
}


# ── Statistical Utilities ────────────────────────────────────────────────────

def bootstrap_auroc_ci(y_true, y_score, n_boot=None, alpha=CI_ALPHA,
                       random_state=42):
    """
    Bootstrap 95% CI for AUROC by resampling (y_true, y_score) pairs.
    Returns (ci_low, ci_high) or (NaN, NaN) if not computable.
    """
    from sklearn.metrics import roc_auc_score

    if n_boot is None:
        n_boot = N_BOOTSTRAP
    if n_boot <= 0:
        return np.nan, np.nan

    rng = np.random.RandomState(random_state)
    n = len(y_true)

    if n < 4 or len(np.unique(y_true)) < 2:
        return np.nan, np.nan

    boot_aurocs = []
    for _ in range(n_boot):
        idx = rng.randint(0, n, size=n)
        y_t = y_true[idx]
        y_s = y_score[idx]
        if len(np.unique(y_t)) < 2:
            continue
        try:
            boot_aurocs.append(roc_auc_score(y_t, y_s))
        except Exception:
            continue

    if len(boot_aurocs) < n_boot * 0.5:
        return np.nan, np.nan

    boot_aurocs = np.array(boot_aurocs)
    ci_low = float(np.percentile(boot_aurocs, 100 * alpha / 2))
    ci_high = float(np.percentile(boot_aurocs, 100 * (1 - alpha / 2)))
    return ci_low, ci_high


def bootstrap_mean_ci(values, n_boot=None, alpha=CI_ALPHA, random_state=42):
    """
    Bootstrap 95% CI for the mean of a set of values.
    Returns (mean, ci_low, ci_high).
    """
    if n_boot is None:
        n_boot = N_BOOTSTRAP

    values = np.array(values)
    mean_val = float(np.mean(values))
    n = len(values)

    if n < 2 or n_boot <= 0:
        return mean_val, np.nan, np.nan

    rng = np.random.RandomState(random_state)
    boot_means = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.randint(0, n, size=n)
        boot_means[i] = np.mean(values[idx])

    ci_low = float(np.percentile(boot_means, 100 * alpha / 2))
    ci_high = float(np.percentile(boot_means, 100 * (1 - alpha / 2)))
    return mean_val, ci_low, ci_high


def permutation_test_two_groups(values_a, values_b, n_perm=None,
                                 random_state=42):
    """
    Two-sided permutation test: H0 = mean(A) == mean(B).
    Returns p-value.
    """
    if n_perm is None:
        n_perm = N_PERMUTATIONS

    rng = np.random.RandomState(random_state)
    a = np.array(values_a)
    b = np.array(values_b)

    observed_diff = abs(np.mean(a) - np.mean(b))
    pooled = np.concatenate([a, b])
    n_a = len(a)

    count = 0
    for _ in range(n_perm):
        rng.shuffle(pooled)
        perm_diff = abs(np.mean(pooled[:n_a]) - np.mean(pooled[n_a:]))
        if perm_diff >= observed_diff:
            count += 1

    return (count + 1) / (n_perm + 1)  # +1 continuity correction


# ── Data Loading ─────────────────────────────────────────────────────────────

def load_tissue_data(tissue: str):
    """Load tissue-wide log2 normalized counts + metadata."""
    tissue_dir = PROCESSED_DIR / tissue
    counts_path = tissue_dir / f"{tissue}_all_missions_log2_norm.csv"
    meta_path = tissue_dir / f"{tissue}_all_missions_metadata.csv"

    if not counts_path.exists() or not meta_path.exists():
        print(f"  [ERROR] Data not found for {tissue}")
        return None

    counts = pd.read_csv(counts_path, index_col=0)
    meta = pd.read_csv(meta_path, index_col=0)

    # Drop non-gene columns
    non_gene_cols = [c for c in counts.columns if c in {"mission", "osd_id", "label"}]
    if non_gene_cols:
        counts = counts.drop(columns=non_gene_cols)
    counts = counts.select_dtypes(include=[np.number])

    # Align
    common = counts.index.intersection(meta.index)
    counts = counts.loc[common]
    meta = meta.loc[common]

    print(f"  Loaded {tissue}: {counts.shape[0]} samples x {counts.shape[1]} genes")
    return counts, meta


def get_binary_labels(meta: pd.DataFrame) -> pd.Series:
    """Extract binary labels: 1=Flight, 0=Ground."""
    label_col = None
    for col in ["label", "label_raw", "condition", "group"]:
        if col in meta.columns:
            label_col = col
            break
    if label_col is None:
        raise ValueError(f"No label column found. Columns: {list(meta.columns)}")

    labels_raw = meta[label_col]
    binary = pd.Series(np.nan, index=meta.index)
    binary[labels_raw == FLIGHT_LABEL] = 1
    binary[labels_raw.isin(GROUND_LABELS)] = 0
    return binary


# ── Transfer Methods ──────────────────────────────────────────────────────────

def ml_transfer(train_X, train_y, test_X, test_y, method="pca_lr"):
    """
    ML-based pairwise transfer.
    Returns dict: {auroc, y_true, y_score} for downstream bootstrap CI.
    """
    from sklearn.metrics import roc_auc_score
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    from sklearn.pipeline import Pipeline

    _nan = {"auroc": np.nan, "y_true": None, "y_score": None}

    X_train = train_X.values.astype(np.float32)
    X_test = test_X.values.astype(np.float32)
    y_train = train_y.values.astype(int)
    y_test = test_y.values.astype(int)

    if len(np.unique(y_test)) < 2 or len(np.unique(y_train)) < 2:
        return _nan

    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)

    if method == "pca_lr":
        n_comps = min(50, X_train.shape[0] - 1, X_train.shape[1])
        if n_comps < 2:
            return _nan
        model = Pipeline([
            ("pca", PCA(n_components=n_comps, random_state=42)),
            ("clf", LogisticRegression(C=1.0, class_weight="balanced",
                                       max_iter=1000, random_state=42))
        ])
    else:  # lr
        model = LogisticRegression(penalty="elasticnet", solver="saga",
                                    l1_ratio=0.5, C=1.0,
                                    class_weight="balanced",
                                    max_iter=2000, random_state=42)

    try:
        model.fit(X_train, y_train)
        y_score = model.predict_proba(X_test)[:, 1]
        auroc = float(roc_auc_score(y_test, y_score))
        return {"auroc": auroc, "y_true": y_test, "y_score": y_score}
    except Exception:
        return _nan


def lfc_signature_transfer(train_X, train_y, test_X, test_y):
    """
    LFC-signature transfer.
    Returns dict: {auroc, y_true, y_score} for downstream bootstrap CI.
    """
    from sklearn.metrics import roc_auc_score

    _nan = {"auroc": np.nan, "y_true": None, "y_score": None}

    y_train = train_y.values.astype(int)
    y_test = test_y.values.astype(int)

    if len(np.unique(y_test)) < 2:
        return _nan

    flight_mask = y_train == 1
    ground_mask = y_train == 0

    if flight_mask.sum() == 0 or ground_mask.sum() == 0:
        return _nan

    flight_mean = train_X.values[flight_mask].mean(axis=0)
    ground_mean = train_X.values[ground_mask].mean(axis=0)
    lfc = flight_mean - ground_mean

    lfc_norm = np.linalg.norm(lfc)
    if lfc_norm == 0:
        return _nan
    lfc_unit = lfc / lfc_norm

    test_centered = test_X.values - ground_mean
    projections = test_centered @ lfc_unit

    try:
        auroc = float(roc_auc_score(y_test, projections))
        return {"auroc": auroc, "y_true": y_test, "y_score": projections}
    except Exception:
        return _nan


# ── Pairwise Transfer Matrix ─────────────────────────────────────────────────

def compute_transfer_matrix(tissue: str, methods=("pca_lr", "lfc"),
                             verbose=True):
    """
    Compute N×N pairwise transfer AUROC matrix with bootstrap CIs.
    Rows = train mission, Columns = test mission.
    """
    result = load_tissue_data(tissue)
    if result is None:
        return None
    counts, meta = result

    binary_labels = get_binary_labels(meta)
    valid = ~binary_labels.isna()
    counts = counts[valid]
    meta = meta[valid]
    binary_labels = binary_labels[valid]

    missions = sorted(meta["mission"].unique())
    n_missions = len(missions)

    if n_missions < 2:
        print(f"  [SKIP] {tissue}: only {n_missions} mission(s)")
        return None

    print(f"  Missions ({n_missions}): {missions}")

    mission_sizes = {}
    for m in missions:
        m_mask = meta["mission"] == m
        n_flt = int((binary_labels[m_mask] == 1).sum())
        n_gnd = int((binary_labels[m_mask] == 0).sum())
        mission_sizes[m] = {"flight": n_flt, "ground": n_gnd,
                            "total": n_flt + n_gnd}
        print(f"    {m}: {n_flt} Flight + {n_gnd} Ground = {n_flt + n_gnd}")

    results = {}
    ci_matrices = {}
    details = []

    for method in methods:
        matrix = pd.DataFrame(np.nan, index=missions, columns=missions)
        ci_low_mat = pd.DataFrame(np.nan, index=missions, columns=missions)
        ci_high_mat = pd.DataFrame(np.nan, index=missions, columns=missions)
        print(f"\n  --- Method: {method} ---")

        for i, train_mission in enumerate(missions):
            for j, test_mission in enumerate(missions):
                if i == j:
                    continue

                train_mask = meta["mission"] == train_mission
                test_mask = meta["mission"] == test_mission

                train_X = counts[train_mask]
                test_X = counts[test_mask]
                train_y = binary_labels[train_mask]
                test_y = binary_labels[test_mask]

                n_flt_train = (train_y == 1).sum()
                n_gnd_train = (train_y == 0).sum()

                # Variance filter on train only (DD-03)
                gene_var = train_X.var(axis=0)
                threshold = gene_var.quantile(VARIANCE_PERCENTILE)
                selected_genes = gene_var[gene_var >= threshold].index.tolist()
                train_X = train_X[selected_genes]
                test_X = test_X[selected_genes]

                if method == "lfc":
                    res = lfc_signature_transfer(train_X, train_y,
                                                  test_X, test_y)
                else:
                    if n_flt_train < MIN_SAMPLES_PER_CLASS or \
                       n_gnd_train < MIN_SAMPLES_PER_CLASS:
                        res = {"auroc": np.nan, "y_true": None, "y_score": None}
                    else:
                        res = ml_transfer(train_X, train_y,
                                          test_X, test_y, method=method)

                auroc = res["auroc"]
                matrix.loc[train_mission, test_mission] = auroc

                # Bootstrap CI for this pair
                ci_low, ci_high = np.nan, np.nan
                if res["y_true"] is not None and not np.isnan(auroc):
                    ci_low, ci_high = bootstrap_auroc_ci(
                        res["y_true"], res["y_score"])
                ci_low_mat.loc[train_mission, test_mission] = ci_low
                ci_high_mat.loc[train_mission, test_mission] = ci_high

                detail = {
                    "train_mission": train_mission,
                    "test_mission": test_mission,
                    "method": method,
                    "auroc": auroc if not np.isnan(auroc) else None,
                    "ci_low": ci_low if not np.isnan(ci_low) else None,
                    "ci_high": ci_high if not np.isnan(ci_high) else None,
                    "n_train": int(train_mask.sum()),
                    "n_test": int(test_mask.sum()),
                    "n_genes": len(selected_genes),
                }
                details.append(detail)

                if verbose:
                    a_str = f"{auroc:.3f}" if not np.isnan(auroc) else "N/A"
                    ci_str = ""
                    if not np.isnan(ci_low):
                        ci_str = f" [{ci_low:.3f}, {ci_high:.3f}]"
                    print(f"    {train_mission} -> {test_mission}: "
                          f"AUROC={a_str}{ci_str} "
                          f"(train={train_mask.sum()}, test={test_mask.sum()})")

        # Summary
        valid_values = matrix.values[~np.isnan(matrix.values)]
        if len(valid_values) > 0:
            mean_val, mean_ci_lo, mean_ci_hi = bootstrap_mean_ci(valid_values)
            print(f"\n  {method} Transfer Matrix Summary:")
            ci_part = ""
            if not np.isnan(mean_ci_lo):
                ci_part = f" 95%CI [{mean_ci_lo:.3f}, {mean_ci_hi:.3f}]"
            print(f"    Mean AUROC = {mean_val:.3f}{ci_part}")
            print(f"    Range = [{np.min(valid_values):.3f}, "
                  f"{np.max(valid_values):.3f}]")
            print(f"    N valid pairs = {len(valid_values)}/"
                  f"{n_missions*(n_missions-1)}")

            # Asymmetry metric
            asymmetries = []
            for ii in range(n_missions):
                for jj in range(ii + 1, n_missions):
                    a = matrix.iloc[ii, jj]
                    b = matrix.iloc[jj, ii]
                    if not np.isnan(a) and not np.isnan(b):
                        asymmetries.append(abs(a - b))
            if asymmetries:
                print(f"    Mean asymmetry |A(i,j) - A(j,i)| = "
                      f"{np.mean(asymmetries):.3f}")

        results[method] = matrix
        ci_matrices[method] = {"low": ci_low_mat, "high": ci_high_mat}

    return {
        "tissue": tissue,
        "missions": missions,
        "n_missions": n_missions,
        "mission_sizes": mission_sizes,
        "matrices": results,
        "ci_matrices": ci_matrices,
        "details": details,
    }


# ── Cross-Tissue Comparison (H1) with Statistics ────────────────────────────

def compare_tissues(all_results: dict, method="pca_lr"):
    """
    H1 test: Compare mean transfer AUROC across tissues with:
    - Bootstrap 95% CI for each tissue mean
    - Pairwise permutation tests
    - Liver sensitivity: with vs without MHU-2 pairs
    """
    print(f"\n{'='*60}")
    print(f"H1 Cross-Tissue Comparison (method={method})")
    print(f"{'='*60}")

    tissue_aurocs = {}
    for tissue, res in all_results.items():
        if res is None or method not in res["matrices"]:
            continue
        mat = res["matrices"][method]
        valid = mat.values[~np.isnan(mat.values)]
        if len(valid) > 0:
            tissue_aurocs[tissue] = valid.tolist()

    # Bootstrap CI for each tissue mean
    tissue_stats = {}
    for tissue, aurocs in tissue_aurocs.items():
        mean_val, ci_low, ci_high = bootstrap_mean_ci(aurocs)
        tissue_stats[tissue] = {
            "mean": mean_val,
            "std": float(np.std(aurocs)),
            "ci_low": ci_low if not np.isnan(ci_low) else None,
            "ci_high": ci_high if not np.isnan(ci_high) else None,
            "n_pairs": len(aurocs),
            "min": float(np.min(aurocs)),
            "max": float(np.max(aurocs)),
        }

    sorted_tissues = sorted(tissue_stats.items(),
                            key=lambda x: x[1]["mean"], reverse=True)

    print(f"\n  Tissue ranking by mean transfer AUROC ({method}):")
    for rank, (tissue, stats) in enumerate(sorted_tissues, 1):
        ci_part = ""
        if stats["ci_low"] is not None:
            ci_part = f" 95%CI [{stats['ci_low']:.3f}, {stats['ci_high']:.3f}]"
        print(f"    {rank}. {tissue:20s} "
              f"mean={stats['mean']:.3f}{ci_part}  "
              f"(n_pairs={stats['n_pairs']})")

    # Pairwise permutation tests
    tissue_names = [t for t, _ in sorted_tissues]
    pairwise_pvalues = {}
    if len(tissue_names) >= 2:
        print(f"\n  Pairwise permutation tests "
              f"(two-sided, {N_PERMUTATIONS} permutations):")
        for t_a, t_b in combinations(tissue_names, 2):
            p_val = permutation_test_two_groups(
                tissue_aurocs[t_a], tissue_aurocs[t_b])
            pair_key = f"{t_a}_vs_{t_b}"
            pairwise_pvalues[pair_key] = p_val
            delta = abs(tissue_stats[t_a]["mean"] - tissue_stats[t_b]["mean"])
            sig = "*" if p_val < 0.05 else ("(*)" if p_val < 0.10 else "")
            print(f"    {t_a:15s} vs {t_b:15s}: "
                  f"delta={delta:.3f}, p={p_val:.4f} {sig}")

    # Liver sensitivity: with vs without MHU-2
    liver_sensitivity = None
    if "liver" in all_results and all_results["liver"] is not None:
        liver_res = all_results["liver"]
        if method in liver_res["matrices"]:
            mat = liver_res["matrices"][method]
            missions = liver_res["missions"]

            if "MHU-2" in missions:
                all_values = []
                excl_values = []
                for m_i in missions:
                    for m_j in missions:
                        if m_i == m_j:
                            continue
                        val = mat.loc[m_i, m_j]
                        if np.isnan(val):
                            continue
                        all_values.append(val)
                        if m_i != "MHU-2" and m_j != "MHU-2":
                            excl_values.append(val)

                if excl_values:
                    m_all, ci_lo_all, ci_hi_all = bootstrap_mean_ci(all_values)
                    m_excl, ci_lo_excl, ci_hi_excl = bootstrap_mean_ci(
                        excl_values)

                    liver_sensitivity = {
                        "with_mhu2": {
                            "mean": m_all,
                            "ci_low": ci_lo_all if not np.isnan(ci_lo_all) else None,
                            "ci_high": ci_hi_all if not np.isnan(ci_hi_all) else None,
                            "n_pairs": len(all_values),
                        },
                        "without_mhu2": {
                            "mean": m_excl,
                            "ci_low": ci_lo_excl if not np.isnan(ci_lo_excl) else None,
                            "ci_high": ci_hi_excl if not np.isnan(ci_hi_excl) else None,
                            "n_pairs": len(excl_values),
                        },
                    }

                    print(f"\n  Liver sensitivity (MHU-2 has only 6 samples):")
                    ci1 = ""
                    if liver_sensitivity["with_mhu2"]["ci_low"] is not None:
                        ci1 = (f" 95%CI [{liver_sensitivity['with_mhu2']['ci_low']:.3f},"
                               f" {liver_sensitivity['with_mhu2']['ci_high']:.3f}]")
                    ci2 = ""
                    if liver_sensitivity["without_mhu2"]["ci_low"] is not None:
                        ci2 = (f" 95%CI [{liver_sensitivity['without_mhu2']['ci_low']:.3f},"
                               f" {liver_sensitivity['without_mhu2']['ci_high']:.3f}]")
                    print(f"    With MHU-2:    mean={m_all:.3f}{ci1} "
                          f"(n={len(all_values)})")
                    print(f"    Without MHU-2: mean={m_excl:.3f}{ci2} "
                          f"(n={len(excl_values)})")

    return {
        "tissue_stats": tissue_stats,
        "pairwise_pvalues": pairwise_pvalues,
        "liver_sensitivity": liver_sensitivity,
    }


# ── Save Results ──────────────────────────────────────────────────────────────

def save_results(tissue_result: dict, verbose=True):
    """Save transfer matrices, CI matrices, and details for one tissue."""
    tissue = tissue_result["tissue"]
    out_dir = B_OUTPUT_DIR / tissue
    out_dir.mkdir(parents=True, exist_ok=True)

    for method, matrix in tissue_result["matrices"].items():
        csv_path = out_dir / f"B_transfer_matrix_{method}.csv"
        matrix.to_csv(csv_path)
        if verbose:
            print(f"  Saved: {csv_path.relative_to(BASE_DIR)}")

        # Save CI matrices
        if method in tissue_result.get("ci_matrices", {}):
            ci = tissue_result["ci_matrices"][method]
            ci_low_path = out_dir / f"B_transfer_ci_low_{method}.csv"
            ci_high_path = out_dir / f"B_transfer_ci_high_{method}.csv"
            ci["low"].to_csv(ci_low_path)
            ci["high"].to_csv(ci_high_path)

    # Save details JSON
    details_path = out_dir / "B_transfer_details.json"
    details = {
        "tissue": tissue,
        "missions": tissue_result["missions"],
        "n_missions": tissue_result["n_missions"],
        "mission_sizes": tissue_result.get("mission_sizes", {}),
        "generated_at": datetime.now().isoformat(),
        "config": {
            "variance_percentile": VARIANCE_PERCENTILE,
            "min_samples_per_class": MIN_SAMPLES_PER_CLASS,
            "n_bootstrap": N_BOOTSTRAP,
            "ci_alpha": CI_ALPHA,
        },
        "details": tissue_result["details"],
    }
    details_path.write_text(json.dumps(details, indent=2, default=str))
    if verbose:
        print(f"  Saved: {details_path.relative_to(BASE_DIR)}")


# ── B Task Directory Generator ────────────────────────────────────────────────

def create_b_task_dirs(tissue: str) -> None:
    """
    Create tasks/B{N}_{tissue}_cross_mission/ with pair_{TRAIN}_{TEST}/test_y.csv.
    Used to populate public benchmark task directories for Category B.
    Only binary samples (Flight vs GC/VC) are included; AG/BC excluded.
    """
    cfg = TISSUE_CONFIG.get(tissue)
    if cfg is None:
        print(f"  [ERROR] Unknown tissue: {tissue}")
        return

    result = load_tissue_data(tissue)
    if result is None:
        return
    _, meta = result

    task_id = cfg["task_id"]  # e.g., "B5"
    binary = get_binary_labels(meta)
    binary_mask = ~binary.isna()
    meta_b = meta[binary_mask]
    labels = binary[binary_mask]

    missions = sorted(meta_b["mission"].unique())
    if len(missions) < 2:
        print(f"  [ERROR] Need >= 2 missions, found: {missions}")
        return

    task_dir = TASKS_DIR / f"{task_id}_{tissue}_cross_mission"
    task_dir.mkdir(parents=True, exist_ok=True)

    n_created = 0
    for train_m in missions:
        for test_m in missions:
            if train_m == test_m:
                continue
            pair_dir = task_dir / f"pair_{train_m}_{test_m}"
            pair_dir.mkdir(exist_ok=True)
            test_mask = meta_b["mission"] == test_m
            test_y = labels[test_mask].astype(int).to_frame(name="label")
            out_path = pair_dir / "test_y.csv"
            test_y.to_csv(out_path)
            n_created += 1

    print(f"  Created {n_created} pairs → {task_dir.relative_to(BASE_DIR)}/")
    print(f"  Missions: {missions}")


# ── CLI ────────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="Category B: Cross-Mission Transfer Matrix"
    )
    parser.add_argument("--tissue", type=str,
                        help="Tissue to process (e.g., liver)")
    parser.add_argument("--all", action="store_true",
                        help="Process all tissues")
    parser.add_argument("--method", type=str, nargs="+",
                        default=["pca_lr", "lfc"],
                        choices=["lr", "pca_lr", "lfc"],
                        help="Transfer methods (default: pca_lr lfc)")
    parser.add_argument("--compare-tissues", action="store_true",
                        help="Run H1 cross-tissue comparison after all tissues")
    parser.add_argument("--quiet", action="store_true")
    parser.add_argument("--no-bootstrap", action="store_true",
                        help="Skip bootstrap CI computation (faster)")
    parser.add_argument("--create-task-dirs", action="store_true",
                        help="Create tasks/B{N}_{tissue}_cross_mission/ directories "
                             "with test_y.csv files (no transfer computation)")
    return parser.parse_args()


def main():
    args = parse_args()
    verbose = not args.quiet

    global N_BOOTSTRAP, N_PERMUTATIONS
    if args.no_bootstrap:
        N_BOOTSTRAP = 0
        N_PERMUTATIONS = 0

    if args.all:
        tissues = list(TISSUE_CONFIG.keys())
    elif args.tissue:
        tissues = [args.tissue]
    else:
        tissues = ["liver"]
        print("No tissue specified. Running liver by default.")

    # ── Create task directories mode (no transfer computation) ────────────────
    if args.create_task_dirs:
        for tissue in tissues:
            print(f"\n{'='*60}")
            print(f"Creating B task dirs: {tissue.upper()}")
            print(f"{'='*60}")
            create_b_task_dirs(tissue)
        print("\nDone.")
        return

    methods = tuple(args.method)
    all_results = {}

    for tissue in tissues:
        print(f"\n{'='*60}")
        print(f"Category B: {tissue.upper()} Cross-Mission Transfer")
        print(f"{'='*60}")

        result = compute_transfer_matrix(tissue, methods=methods,
                                          verbose=verbose)
        if result is not None:
            save_results(result, verbose=verbose)
            all_results[tissue] = result

    # Cross-tissue comparison (H1) with statistics
    h1_results = None
    if args.compare_tissues or args.all:
        if len(all_results) >= 2:
            h1_by_method = {}
            for method in methods:
                h1_by_method[method] = compare_tissues(all_results,
                                                        method=method)
            h1_results = h1_by_method
        else:
            print("\n  [SKIP] Need >= 2 tissues for H1 comparison")

    # Save overall summary
    if all_results:
        RESULTS_DIR.mkdir(parents=True, exist_ok=True)
        summary = {}
        for tissue, res in all_results.items():
            tissue_summary = {
                "tissue": tissue,
                "n_missions": res["n_missions"],
                "missions": res["missions"],
                "mission_sizes": res.get("mission_sizes", {}),
                "methods": {},
            }
            for method, matrix in res["matrices"].items():
                valid = matrix.values[~np.isnan(matrix.values)]
                if len(valid) > 0:
                    mean_val, ci_lo, ci_hi = bootstrap_mean_ci(valid)
                    tissue_summary["methods"][method] = {
                        "mean_auroc": mean_val,
                        "std_auroc": float(np.std(valid)),
                        "ci_low": ci_lo if not np.isnan(ci_lo) else None,
                        "ci_high": ci_hi if not np.isnan(ci_hi) else None,
                        "n_valid_pairs": len(valid),
                    }
            summary[tissue] = tissue_summary

        # Add H1 test results
        if h1_results is not None:
            summary["_h1_test"] = {}
            for method, h1_res in h1_results.items():
                summary["_h1_test"][method] = {
                    "pairwise_pvalues": h1_res["pairwise_pvalues"],
                    "liver_sensitivity": h1_res["liver_sensitivity"],
                }

        summary_path = RESULTS_DIR / "B_cross_mission_summary.json"
        # Merge with existing summary (don't overwrite other tissues)
        if summary_path.exists():
            existing = json.loads(summary_path.read_text())
            existing.update(summary)
            summary = existing
        summary_path.write_text(json.dumps(summary, indent=2))
        print(f"\n  Summary saved: {summary_path.relative_to(BASE_DIR)}")

    print("\nDone.")


if __name__ == "__main__":
    main()
