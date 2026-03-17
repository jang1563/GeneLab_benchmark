#!/usr/bin/env python3
"""
temporal_analysis.py — GeneLab_benchmark v2.0: Temporal & Covariate Analysis

T1: ISS-T vs LAR classification
  - Tests whether sacrifice timing (ISS-Terminal vs Live Animal Return)
    creates detectable transcriptomic differences
  - GC serves as preservation artifact baseline (DD-18, PMID 33376967)
  - Sub-tasks: T1a (RR-6 liver), T1b (RR-8 liver), T1c (RR-6 thymus),
               T1d (cross-mission RR-6↔RR-8)

T2: LAR recovery signature
  - Quantifies recovery toward baseline in LAR samples
  - Uses preservation-matched comparisons: FLT_ISS-T vs BSL_ISS-T, FLT_LAR vs BSL_LAR
  - PCA trajectory: BSL → FLT_ISS-T → FLT_LAR
  - Per-pathway recovery fraction (DD-20)

T3: Age × spaceflight interaction (RR-8 liver only)
  - Two-way ANOVA: pathway ~ flight + age + flight:age
  - Age classification within condition groups
  - Age-stratified spaceflight effect comparison

Output:
  v2/processed/T_temporal/
    T1_{tissue}_{mission}_results.json
    T1_cross_mission_results.json
    T2_{tissue}_{mission}_recovery.json
    T2_pca_trajectory.pdf
    T3_rr8_liver_results.json
  v2/evaluation/
    T_temporal_summary.json

Usage:
  python v2/scripts/temporal_analysis.py --task T1
  python v2/scripts/temporal_analysis.py --task T2
  python v2/scripts/temporal_analysis.py --task T3
  python v2/scripts/temporal_analysis.py --all
"""

import json
import sys
import argparse
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from collections import Counter

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")

# v2/scripts/ → import v1 utils from scripts/
V2_DIR = Path(__file__).resolve().parent.parent          # v2/
PROJECT_DIR = V2_DIR.parent                               # GeneLab_benchmark/
sys.path.insert(0, str(PROJECT_DIR / "scripts"))

from utils import (
    load_temporal_metadata, load_gene_features, load_pathway_features,
    align_features_with_meta,
)

# ── Paths ──────────────────────────────────────────────────────────────────────
# v1 data (read-only)
PROCESSED_DIR = PROJECT_DIR / "processed" / "A_detection"
PATHWAY_DIR = PROJECT_DIR / "processed" / "pathway_scores"
# v2 output
T_OUTPUT_DIR = V2_DIR / "processed" / "T_temporal"
RESULTS_DIR = V2_DIR / "evaluation"

# ── Config ─────────────────────────────────────────────────────────────────────
VARIANCE_PERCENTILE = 0.25   # keep top 75% by variance (train-only)
N_BOOTSTRAP = 2000
CI_ALPHA = 0.05
N_SPLITS = 5
N_REPEATS = 10
N_PERM = 10000
RECOVERY_MIN_DELTA = 0.1    # min |delta_flight| for recovery fraction (DD-20)


# ── Statistical Utilities ────────────────────────────────────────────────────
# Re-implemented here (not imported from condition_prediction.py) per plan decision:
# condition_prediction.py functions are tightly coupled to module constants.

def bootstrap_auroc_ci(y_true, y_score, n_boot=N_BOOTSTRAP, alpha=CI_ALPHA,
                       random_state=42):
    """Bootstrap 95% CI for AUROC."""
    from sklearn.metrics import roc_auc_score
    rng = np.random.RandomState(random_state)
    n = len(y_true)
    boot_aucs = []
    for _ in range(n_boot):
        idx = rng.randint(0, n, size=n)
        yt, ys = y_true[idx], y_score[idx]
        if len(np.unique(yt)) < 2:
            continue
        try:
            boot_aucs.append(roc_auc_score(yt, ys))
        except Exception:
            continue
    if len(boot_aucs) < 100:
        return np.nan, np.nan
    lo = np.percentile(boot_aucs, 100 * alpha / 2)
    hi = np.percentile(boot_aucs, 100 * (1 - alpha / 2))
    return lo, hi


def permutation_test_auroc(y_true, y_score, n_perm=N_PERM, random_state=42):
    """Permutation test: is AUROC significantly different from 0.5?"""
    from sklearn.metrics import roc_auc_score
    observed = roc_auc_score(y_true, y_score)
    rng = np.random.RandomState(random_state)
    count_ge = 0
    for _ in range(n_perm):
        perm = rng.permutation(y_true)
        try:
            pauc = roc_auc_score(perm, y_score)
        except ValueError:
            continue
        if pauc >= observed:
            count_ge += 1
    return (count_ge + 1) / (n_perm + 1)


# ── Core Classification ───────────────────────────────────────────────────────

def classify_binary_rskf(X, y, feature_mode="gene",
                         n_splits=N_SPLITS, n_repeats=N_REPEATS):
    """Binary classification with RepeatedStratifiedKFold.
    Returns AUROC, predictions, probabilities.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    from sklearn.model_selection import RepeatedStratifiedKFold
    from sklearn.metrics import roc_auc_score

    rskf = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats,
                                   random_state=42)
    n_samples = len(y)
    prob_accum = np.zeros(n_samples)
    count_accum = np.zeros(n_samples)

    for train_idx, test_idx in rskf.split(X, y):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train = y[train_idx]

        if feature_mode == "gene":
            var = np.var(X_train, axis=0)
            threshold = np.percentile(var, VARIANCE_PERCENTILE * 100)
            keep = var > threshold
            if keep.sum() < 50:
                keep = np.ones(X_train.shape[1], dtype=bool)
            X_train, X_test = X_train[:, keep], X_test[:, keep]
            n_comp = min(50, X_train.shape[0] - 1, X_train.shape[1])
            scaler = StandardScaler()
            X_train = scaler.fit_transform(X_train)
            X_test = scaler.transform(X_test)
            pca = PCA(n_components=n_comp, random_state=42)
            X_train = pca.fit_transform(X_train)
            X_test = pca.transform(X_test)
        else:
            scaler = StandardScaler()
            X_train = scaler.fit_transform(X_train)
            X_test = scaler.transform(X_test)

        clf = LogisticRegression(
            solver="lbfgs", class_weight="balanced",
            max_iter=2000, C=1.0, random_state=42,
        )
        clf.fit(X_train, y_train)
        proba = clf.predict_proba(X_test)
        # Probability of positive class (class 1)
        pos_idx = list(clf.classes_).index(1)
        prob_accum[test_idx] += proba[:, pos_idx]
        count_accum[test_idx] += 1

    prob_avg = prob_accum / np.maximum(count_accum, 1)
    y_pred = (prob_avg >= 0.5).astype(int)

    auroc = roc_auc_score(y, prob_avg)
    return auroc, y_pred, prob_avg


def classify_binary_loo(X, y, feature_mode="gene"):
    """Binary classification with Leave-One-Out CV."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    from sklearn.model_selection import LeaveOneOut
    from sklearn.metrics import roc_auc_score

    loo = LeaveOneOut()
    n_samples = len(y)
    y_pred = np.empty(n_samples, dtype=int)
    y_proba = np.zeros(n_samples)

    for train_idx, test_idx in loo.split(X, y):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train = y[train_idx]

        if feature_mode == "gene":
            # Variance filter BEFORE scaling (on raw data)
            var = np.var(X_train, axis=0)
            threshold = np.percentile(var, VARIANCE_PERCENTILE * 100)
            keep = var > threshold
            if keep.sum() < 20:
                keep = np.ones(X_train.shape[1], dtype=bool)
            X_train, X_test = X_train[:, keep], X_test[:, keep]

            scaler = StandardScaler()
            X_train = scaler.fit_transform(X_train)
            X_test = scaler.transform(X_test)

            n_comp = min(5, X_train.shape[0] - 1, X_train.shape[1])
            pca = PCA(n_components=n_comp, random_state=42)
            X_train = pca.fit_transform(X_train)
            X_test = pca.transform(X_test)
        else:
            scaler = StandardScaler()
            X_train = scaler.fit_transform(X_train)
            X_test = scaler.transform(X_test)

        clf = LogisticRegression(
            solver="lbfgs", class_weight="balanced",
            max_iter=2000, C=1.0, random_state=42,
        )
        clf.fit(X_train, y_train)
        y_pred[test_idx] = clf.predict(X_test)
        proba = clf.predict_proba(X_test)
        pos_idx = list(clf.classes_).index(1)
        y_proba[test_idx] = proba[:, pos_idx]

    auroc = roc_auc_score(y, y_proba)
    return auroc, y_pred, y_proba


def classify_train_test(X_train, y_train, X_test, y_test, feature_mode="gene"):
    """Train-all/test-all classification for cross-mission transfer."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    from sklearn.metrics import roc_auc_score

    if feature_mode == "gene":
        var = np.var(X_train, axis=0)
        threshold = np.percentile(var, VARIANCE_PERCENTILE * 100)
        keep = var > threshold
        if keep.sum() < 50:
            keep = np.ones(X_train.shape[1], dtype=bool)
        X_train, X_test = X_train[:, keep], X_test[:, keep]
        n_comp = min(50, X_train.shape[0] - 1, X_train.shape[1])
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)
        pca = PCA(n_components=n_comp, random_state=42)
        X_train = pca.fit_transform(X_train)
        X_test = pca.transform(X_test)
    else:
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)

    clf = LogisticRegression(
        solver="lbfgs", class_weight="balanced",
        max_iter=2000, C=1.0, random_state=42,
    )
    clf.fit(X_train, y_train)
    proba = clf.predict_proba(X_test)
    pos_idx = list(clf.classes_).index(1)
    y_proba = proba[:, pos_idx]
    y_pred = clf.predict(X_test)

    auroc = roc_auc_score(y_test, y_proba)
    return auroc, y_pred, y_proba


# ── JSON encoder ──────────────────────────────────────────────────────────────

class _NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.bool_):
            return bool(obj)
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj) if np.isfinite(obj) else None
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)


def save_json(data, path):
    """Save dict as JSON with numpy encoding."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, cls=_NumpyEncoder))
    print(f"  Saved: {path}")


# ── T1: ISS-T vs LAR Classification ─────────────────────────────────────────

def _run_t1_subtask(tissue, mission, condition, feature_mode, gene_feat, pathway_feat, meta):
    """Run a single T1 sub-task: ISS-T vs LAR within one condition.

    Args:
        condition: 'FLT', 'GC', or 'BSL'

    Returns dict with AUROC, CI, p-value.
    """
    from sklearn.metrics import roc_auc_score

    # Filter to condition + known timing
    mask_condition = meta["label"].str.upper().str.startswith(condition.upper())
    # Handle special labels: 'Flight' → 'FLT', 'BC' (baseline control)
    if condition == "FLT":
        mask_condition = meta["label"].isin(["Flight", "FLT"])
    elif condition == "GC":
        mask_condition = meta["label"] == "GC"
    elif condition == "BSL":
        mask_condition = meta["label"].isin(["BC", "BSL"])

    mask_timing = meta["sacrifice_timing"].isin(["ISS-T", "LAR"])
    mask = mask_condition & mask_timing

    meta_sub = meta[mask].copy()
    if len(meta_sub) < 6:
        return {"status": "SKIP", "reason": f"n={len(meta_sub)} too small",
                "n_total": len(meta_sub)}

    # Binary label: ISS-T=0, LAR=1
    y = (meta_sub["sacrifice_timing"] == "LAR").astype(int).values
    n_isst = (y == 0).sum()
    n_lar = (y == 1).sum()

    if n_isst < 2 or n_lar < 2:
        return {"status": "SKIP", "reason": f"ISS-T={n_isst}, LAR={n_lar}",
                "n_isst": n_isst, "n_lar": n_lar}

    feat_df = gene_feat if feature_mode == "gene" else pathway_feat
    if feat_df is None:
        return {"status": "SKIP", "reason": f"no {feature_mode} features"}

    try:
        feat_aligned, meta_aligned = align_features_with_meta(feat_df, meta_sub)
    except ValueError as e:
        return {"status": "SKIP", "reason": str(e)}

    X = feat_aligned.values.astype(float)
    X = np.nan_to_num(X, nan=0.0)
    y = (meta_aligned["sacrifice_timing"] == "LAR").astype(int).values

    n_total = len(y)
    n_isst = (y == 0).sum()
    n_lar = (y == 1).sum()

    # Choose CV method based on sample size
    min_class = min(n_isst, n_lar)
    if min_class >= 5:
        auroc, y_pred, y_proba = classify_binary_rskf(
            X, y, feature_mode=feature_mode,
            n_splits=min(N_SPLITS, min_class), n_repeats=N_REPEATS
        )
        cv_method = f"RepeatedStratifiedKFold({min(N_SPLITS, min_class)}x{N_REPEATS})"
    else:
        auroc, y_pred, y_proba = classify_binary_loo(X, y, feature_mode=feature_mode)
        cv_method = "LOO-CV"

    ci_lo, ci_hi = bootstrap_auroc_ci(y, y_proba)
    perm_p = permutation_test_auroc(y, y_proba)

    return {
        "status": "OK",
        "auroc": auroc,
        "ci_low": ci_lo,
        "ci_high": ci_hi,
        "perm_p": perm_p,
        "n_total": n_total,
        "n_isst": int(n_isst),
        "n_lar": int(n_lar),
        "cv_method": cv_method,
        "n_features": X.shape[1],
    }


def run_t1():
    """Run T1: ISS-T vs LAR classification across tissues/missions."""
    print("=" * 70)
    print("T1: ISS-T vs LAR Classification")
    print("  Confound warning: ISS-T vs LAR confounded with preservation method")
    print("  (RNAlater vs standard necropsy). PMID 33376967.")
    print("  GC results serve as preservation artifact baseline (DD-18).")
    print("=" * 70)

    T_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    all_results = {}

    # ── T1a/b: Liver (RR-6, RR-8) ──
    for mission in ["RR-6", "RR-8"]:
        task_label = "T1a" if mission == "RR-6" else "T1b"
        print(f"\n{'─' * 50}")
        print(f"  {task_label}: {mission} liver — ISS-T vs LAR")
        print(f"{'─' * 50}")

        meta = load_temporal_metadata("liver", mission=mission)
        gene_feat = load_gene_features("liver")
        pathway_feat = load_pathway_features("liver")

        results = {
            "task": task_label,
            "tissue": "liver",
            "mission": mission,
            "confound_warning": "PMID 33376967: ISS-T=RNAlater, LAR=necropsy",
            "conditions": {},
        }

        for condition in ["FLT", "GC", "BSL"]:
            for mode in ["gene", "pathway"]:
                key = f"{condition}_{mode}"
                print(f"\n  [{key}] ", end="")
                r = _run_t1_subtask(
                    "liver", mission, condition, mode,
                    gene_feat, pathway_feat, meta
                )
                if r["status"] == "OK":
                    print(f"AUROC={r['auroc']:.3f} [{r['ci_low']:.3f},{r['ci_high']:.3f}] "
                          f"p={r['perm_p']:.4f} (n={r['n_total']}: ISS-T={r['n_isst']}, LAR={r['n_lar']})")
                else:
                    print(f"SKIP: {r['reason']}")
                results["conditions"][key] = r

        # Interpretation
        flt_gene = results["conditions"].get("FLT_gene", {})
        gc_gene = results["conditions"].get("GC_gene", {})
        if flt_gene.get("status") == "OK" and gc_gene.get("status") == "OK":
            flt_auc = flt_gene["auroc"]
            gc_auc = gc_gene["auroc"]
            excess = flt_auc - gc_auc
            if excess > 0.05:
                note = ("excess > 0: FLT AUROC exceeds GC baseline → "
                        "potential biological ISS-T vs LAR difference beyond preservation artifact")
            elif excess < -0.05:
                note = ("excess < 0: GC AUROC exceeds FLT AUROC → "
                        "preservation artifact dominates; no evidence of biological timing effect")
            else:
                note = ("excess ≈ 0: FLT and GC AUROC comparable → "
                        "ISS-T vs LAR difference is primarily preservation artifact")
            results["interpretation"] = {
                "FLT_AUROC": flt_auc,
                "GC_AUROC_baseline": gc_auc,
                "excess_signal": excess,
                "note": note,
            }
            print(f"\n  Interpretation: FLT AUROC={flt_auc:.3f}, GC baseline={gc_auc:.3f}, excess={excess:+.3f}")

        save_json(results, T_OUTPUT_DIR / f"T1_{mission}_liver_results.json")
        all_results[f"{task_label}_{mission}_liver"] = results

    # ── T1c: Thymus (RR-6) ──
    print(f"\n{'─' * 50}")
    print(f"  T1c: RR-6 thymus — ISS-T vs LAR")
    print(f"{'─' * 50}")

    meta_thy = load_temporal_metadata("thymus", mission="RR-6")
    gene_feat_thy = load_gene_features("thymus")
    pathway_feat_thy = load_pathway_features("thymus")

    results_thy = {
        "task": "T1c",
        "tissue": "thymus",
        "mission": "RR-6",
        "confound_warning": "PMID 33376967",
        "conditions": {},
    }

    for condition in ["FLT", "GC"]:
        for mode in ["gene", "pathway"]:
            key = f"{condition}_{mode}"
            print(f"\n  [{key}] ", end="")
            r = _run_t1_subtask(
                "thymus", "RR-6", condition, mode,
                gene_feat_thy, pathway_feat_thy, meta_thy
            )
            if r["status"] == "OK":
                print(f"AUROC={r['auroc']:.3f} [{r['ci_low']:.3f},{r['ci_high']:.3f}] "
                      f"p={r['perm_p']:.4f} (n={r['n_total']})")
            else:
                print(f"SKIP: {r['reason']}")
            results_thy["conditions"][key] = r

    save_json(results_thy, T_OUTPUT_DIR / f"T1_RR-6_thymus_results.json")
    all_results["T1c_RR-6_thymus"] = results_thy

    # ── T1d: Cross-mission transfer (RR-6 ↔ RR-8 liver) ──
    print(f"\n{'─' * 50}")
    print(f"  T1d: Cross-mission ISS-T vs LAR transfer (RR-6 ↔ RR-8 liver)")
    print(f"{'─' * 50}")

    results_xm = {
        "task": "T1d",
        "tissue": "liver",
        "description": "Cross-mission temporal transfer",
        "transfers": {},
    }

    meta_rr6 = load_temporal_metadata("liver", mission="RR-6")
    meta_rr8 = load_temporal_metadata("liver", mission="RR-8")
    gene_feat_liver = load_gene_features("liver")

    for condition in ["FLT", "GC"]:
        for mode in ["gene", "pathway"]:
            feat_df = gene_feat_liver if mode == "gene" else load_pathway_features("liver")
            if feat_df is None:
                continue

            for train_m, test_m, train_meta, test_meta in [
                ("RR-6", "RR-8", meta_rr6, meta_rr8),
                ("RR-8", "RR-6", meta_rr8, meta_rr6),
            ]:
                # Filter condition + timing
                if condition == "FLT":
                    mask_tr = train_meta["label"].isin(["Flight", "FLT"]) & train_meta["sacrifice_timing"].isin(["ISS-T", "LAR"])
                    mask_te = test_meta["label"].isin(["Flight", "FLT"]) & test_meta["sacrifice_timing"].isin(["ISS-T", "LAR"])
                else:
                    mask_tr = (train_meta["label"] == condition) & train_meta["sacrifice_timing"].isin(["ISS-T", "LAR"])
                    mask_te = (test_meta["label"] == condition) & test_meta["sacrifice_timing"].isin(["ISS-T", "LAR"])

                meta_tr = train_meta[mask_tr]
                meta_te = test_meta[mask_te]

                if len(meta_tr) < 4 or len(meta_te) < 4:
                    continue

                try:
                    feat_tr, meta_tr_a = align_features_with_meta(feat_df, meta_tr)
                    feat_te, meta_te_a = align_features_with_meta(feat_df, meta_te)
                except ValueError:
                    continue

                # Ensure same features
                common_feats = sorted(set(feat_tr.columns) & set(feat_te.columns))
                X_train = feat_tr[common_feats].values.astype(float)
                X_test = feat_te[common_feats].values.astype(float)
                X_train = np.nan_to_num(X_train, nan=0.0)
                X_test = np.nan_to_num(X_test, nan=0.0)

                y_train = (meta_tr_a["sacrifice_timing"] == "LAR").astype(int).values
                y_test = (meta_te_a["sacrifice_timing"] == "LAR").astype(int).values

                if len(np.unique(y_train)) < 2 or len(np.unique(y_test)) < 2:
                    continue

                auroc, _, y_proba = classify_train_test(
                    X_train, y_train, X_test, y_test, feature_mode=mode
                )
                ci_lo, ci_hi = bootstrap_auroc_ci(y_test, y_proba)

                key = f"{condition}_{mode}_{train_m}→{test_m}"
                results_xm["transfers"][key] = {
                    "train_mission": train_m,
                    "test_mission": test_m,
                    "condition": condition,
                    "feature_mode": mode,
                    "auroc": auroc,
                    "ci_low": ci_lo,
                    "ci_high": ci_hi,
                    "n_train": len(y_train),
                    "n_test": len(y_test),
                }
                print(f"  [{key}] AUROC={auroc:.3f} [{ci_lo:.3f},{ci_hi:.3f}] "
                      f"(train={len(y_train)}, test={len(y_test)})")

    save_json(results_xm, T_OUTPUT_DIR / f"T1_cross_mission_results.json")
    all_results["T1d_cross_mission"] = results_xm

    return all_results


# ── T2: LAR Recovery Signature ───────────────────────────────────────────────

def run_t2():
    """Run T2: Recovery analysis for LAR samples."""
    print("=" * 70)
    print("T2: LAR Recovery Signature")
    print("  Preservation-matched analysis: FLT_X vs BSL_X (same method)")
    print("  Recovery fraction (DD-20): direction-aware 1 - (delta_return / delta_flight)")
    print("=" * 70)

    T_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    all_results = {}

    for mission in ["RR-6", "RR-8"]:
        print(f"\n{'─' * 50}")
        print(f"  T2: {mission} liver recovery analysis")
        print(f"{'─' * 50}")

        meta = load_temporal_metadata("liver", mission=mission)

        # Check BSL availability
        bsl_mask = meta["label"].isin(["BC", "BSL"])
        timing_mask = meta["sacrifice_timing"].isin(["ISS-T", "LAR"])
        n_bsl = (bsl_mask & timing_mask).sum()
        if n_bsl < 4:
            print(f"  SKIP: insufficient BSL samples with timing info (n={n_bsl})")
            continue

        # ── T2a: PCA trajectory ──
        print(f"\n  [PCA trajectory]")

        gene_feat = load_gene_features("liver")
        try:
            feat_aligned, meta_aligned = align_features_with_meta(gene_feat, meta)
        except ValueError as e:
            print(f"  SKIP: {e}")
            continue

        # Define groups
        groups = {}
        for label_val, condition_name in [("Flight", "FLT"), ("FLT", "FLT"),
                                           ("GC", "GC"), ("BC", "BSL"), ("BSL", "BSL")]:
            for timing in ["ISS-T", "LAR"]:
                mask = (meta_aligned["label"] == label_val) & (meta_aligned["sacrifice_timing"] == timing)
                grp_name = f"{condition_name}_{timing}"
                if mask.sum() > 0 and grp_name not in groups:
                    groups[grp_name] = mask

        print(f"  Groups found: {[f'{k}(n={v.sum()})' for k,v in groups.items()]}")

        # PCA
        from sklearn.preprocessing import StandardScaler
        from sklearn.decomposition import PCA

        X_all = feat_aligned.values.astype(float)
        X_all = np.nan_to_num(X_all, nan=0.0)

        # Variance filter
        var = np.var(X_all, axis=0)
        threshold = np.percentile(var, VARIANCE_PERCENTILE * 100)
        keep = var > threshold
        X_filt = X_all[:, keep]

        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X_filt)
        pca = PCA(n_components=min(10, X_scaled.shape[0] - 1), random_state=42)
        X_pca = pca.fit_transform(X_scaled)

        # Compute centroids
        centroids = {}
        for grp_name, mask in groups.items():
            if mask.sum() >= 2:
                centroids[grp_name] = X_pca[mask].mean(axis=0)

        # PCA distances for trajectory
        pca_distances = {}
        for pair in [("BSL_ISS-T", "FLT_ISS-T"), ("BSL_LAR", "FLT_LAR"),
                     ("FLT_ISS-T", "FLT_LAR"), ("BSL_ISS-T", "BSL_LAR"),
                     ("BSL_ISS-T", "FLT_LAR"), ("BSL_LAR", "FLT_ISS-T")]:
            if pair[0] in centroids and pair[1] in centroids:
                dist = np.linalg.norm(centroids[pair[0]] - centroids[pair[1]])
                pca_distances[f"{pair[0]}→{pair[1]}"] = dist

        print(f"  PCA distances:")
        for k, v in pca_distances.items():
            print(f"    {k}: {v:.3f}")

        # Recovery ratio (PCA-based)
        d_flight_isst = pca_distances.get("BSL_ISS-T→FLT_ISS-T", None)
        d_flight_lar = pca_distances.get("BSL_LAR→FLT_LAR", None)
        pca_recovery = None
        if d_flight_isst and d_flight_lar and d_flight_isst > 0.1:
            pca_recovery = d_flight_lar / d_flight_isst
            print(f"\n  PCA Recovery ratio R = d(FLT_LAR,BSL_LAR)/d(FLT_ISS-T,BSL_ISS-T) = {pca_recovery:.3f}")
            if pca_recovery < 1:
                print(f"    → R < 1: Recovery toward baseline")
            elif pca_recovery > 1:
                print(f"    → R > 1: Continued divergence")
            else:
                print(f"    → R ≈ 1: No recovery")

        # ── T2b: Per-pathway recovery fraction (GSVA) ──
        print(f"\n  [Per-pathway recovery (GSVA Hallmark)]")

        gsva_path = PATHWAY_DIR / "liver" / f"{mission}_gsva_hallmark.csv"
        if not gsva_path.exists():
            print(f"  SKIP: {gsva_path} not found")
            pathway_recovery = {}
        else:
            gsva = pd.read_csv(gsva_path, index_col=0)
            if "mission" in gsva.columns:
                gsva = gsva.drop(columns=["mission"])

            try:
                gsva_aligned, meta_gsva = align_features_with_meta(gsva, meta)
            except ValueError as e:
                print(f"  SKIP GSVA: {e}")
                gsva_aligned = None

            if gsva_aligned is not None:
                pathway_recovery = {}
                pathways = gsva_aligned.columns.tolist()

                for pw in pathways:
                    scores = gsva_aligned[pw]

                    # Preservation-matched deltas
                    bsl_isst = scores[(meta_gsva["label"].isin(["BC", "BSL"])) &
                                      (meta_gsva["sacrifice_timing"] == "ISS-T")]
                    flt_isst = scores[(meta_gsva["label"].isin(["Flight", "FLT"])) &
                                      (meta_gsva["sacrifice_timing"] == "ISS-T")]
                    bsl_lar = scores[(meta_gsva["label"].isin(["BC", "BSL"])) &
                                     (meta_gsva["sacrifice_timing"] == "LAR")]
                    flt_lar = scores[(meta_gsva["label"].isin(["Flight", "FLT"])) &
                                     (meta_gsva["sacrifice_timing"] == "LAR")]

                    if len(bsl_isst) < 2 or len(flt_isst) < 2 or len(bsl_lar) < 2 or len(flt_lar) < 2:
                        continue

                    delta_flight = flt_isst.mean() - bsl_isst.mean()
                    delta_return = flt_lar.mean() - bsl_lar.mean()

                    if abs(delta_flight) < RECOVERY_MIN_DELTA:
                        pathway_recovery[pw] = {
                            "delta_flight": delta_flight,
                            "delta_return": delta_return,
                            "recovery_fraction": None,
                            "note": f"|delta_flight|={abs(delta_flight):.3f} < {RECOVERY_MIN_DELTA}"
                        }
                        continue

                    # Direction-aware recovery fraction (DD-20):
                    # Same sign: recovery_frac = 1 - |delta_return/delta_flight|
                    #   1.0 = complete recovery, 0.0 = no recovery
                    # Opposite sign (direction reversal = overshoot past baseline):
                    #   recovery_frac = 1 + |delta_return/delta_flight| (but negative)
                    #   → results in recovery_frac > 1 or < 0 depending on magnitude
                    ratio = delta_return / delta_flight
                    if ratio >= 0:
                        # Same direction: partial recovery (0-1) or continued divergence (<0)
                        recovery_frac = 1.0 - ratio
                    else:
                        # Opposite direction: overshoot past baseline
                        recovery_frac = 1.0 - ratio  # = 1 + |ratio|, always > 1

                    pathway_recovery[pw] = {
                        "delta_flight": delta_flight,
                        "delta_return": delta_return,
                        "recovery_fraction": recovery_frac,
                        "direction_reversed": bool(ratio < 0),
                    }

                # Summary
                valid_recoveries = [v["recovery_fraction"] for v in pathway_recovery.values()
                                    if v["recovery_fraction"] is not None]
                if valid_recoveries:
                    mean_rec = np.mean(valid_recoveries)
                    n_recovering = sum(1 for r in valid_recoveries if r > 0)
                    n_total = len(valid_recoveries)
                    print(f"  Mean recovery fraction: {mean_rec:.3f}")
                    print(f"  Recovering pathways: {n_recovering}/{n_total}")

                    # Top 5 most recovered
                    sorted_pw = sorted(
                        [(k, v["recovery_fraction"]) for k, v in pathway_recovery.items()
                         if v["recovery_fraction"] is not None],
                        key=lambda x: x[1], reverse=True
                    )
                    print(f"\n  Top 5 most recovered:")
                    for pw, rf in sorted_pw[:5]:
                        print(f"    {pw}: {rf:.3f}")
                    print(f"\n  Bottom 5 (least recovered / overshoot):")
                    for pw, rf in sorted_pw[-5:]:
                        print(f"    {pw}: {rf:.3f}")
            else:
                pathway_recovery = {}

        # ── T2c: Classification-based recovery ──
        print(f"\n  [Classification recovery]")

        # Train BSL vs FLT classifier (ISS-T only, preservation-matched)
        bsl_isst_mask = meta_aligned["label"].isin(["BC", "BSL"]) & (meta_aligned["sacrifice_timing"] == "ISS-T")
        flt_isst_mask = meta_aligned["label"].isin(["Flight", "FLT"]) & (meta_aligned["sacrifice_timing"] == "ISS-T")
        flt_lar_mask = meta_aligned["label"].isin(["Flight", "FLT"]) & (meta_aligned["sacrifice_timing"] == "LAR")

        n_bsl_isst = bsl_isst_mask.sum()
        n_flt_isst = flt_isst_mask.sum()
        n_flt_lar = flt_lar_mask.sum()

        clf_shift = None
        if n_bsl_isst >= 3 and n_flt_isst >= 3 and n_flt_lar >= 2:
            train_mask = bsl_isst_mask | flt_isst_mask
            X_train_rec = X_all[train_mask]
            # Derive labels from actual metadata (not positional assumption)
            y_train_rec = meta_aligned.loc[train_mask, "label"].isin(
                ["Flight", "FLT"]
            ).astype(int).values
            X_test_rec = X_all[flt_lar_mask]

            from sklearn.preprocessing import StandardScaler
            from sklearn.decomposition import PCA as PCA_
            from sklearn.linear_model import LogisticRegression

            scaler2 = StandardScaler()
            X_tr_s = scaler2.fit_transform(X_train_rec[:, keep])
            X_te_s = scaler2.transform(X_test_rec[:, keep])
            n_comp2 = min(10, X_tr_s.shape[0] - 1, X_tr_s.shape[1])
            pca2 = PCA_(n_components=n_comp2, random_state=42)
            X_tr_p = pca2.fit_transform(X_tr_s)
            X_te_p = pca2.transform(X_te_s)

            clf = LogisticRegression(solver="lbfgs", max_iter=2000, C=1.0, random_state=42)
            clf.fit(X_tr_p, y_train_rec)

            # FLT_LAR flight probability
            flt_lar_proba = clf.predict_proba(X_te_p)[:, list(clf.classes_).index(1)]
            mean_flt_lar_prob = flt_lar_proba.mean()

            # FLT_ISS-T flight probability (should be ~1.0)
            flt_isst_proba = clf.predict_proba(pca2.transform(scaler2.transform(
                X_all[flt_isst_mask][:, keep]
            )))[:, list(clf.classes_).index(1)]
            mean_flt_isst_prob = flt_isst_proba.mean()

            clf_shift = {
                "analysis_type": "descriptive_projection",
                "reference_set": "BSL_ISS-T vs FLT_ISS-T training samples",
                "caveat": (
                    "FLT_ISS-T probabilities are in-sample training-set scores; "
                    "use as descriptive reference, not held-out validation evidence"
                ),
                "mean_flt_isst_flight_prob": mean_flt_isst_prob,
                "mean_flt_lar_flight_prob": mean_flt_lar_prob,
                "interpretation": "FLT_LAR projects closer to baseline than FLT_ISS-T"
                                  if mean_flt_lar_prob < mean_flt_isst_prob
                                  else "No baselineward projection shift detected"
            }
            print(f"  FLT_ISS-T mean flight_prob: {mean_flt_isst_prob:.3f}")
            print(f"  FLT_LAR mean flight_prob: {mean_flt_lar_prob:.3f}")
            if mean_flt_lar_prob < mean_flt_isst_prob:
                print(f"  → Descriptive shift: LAR prob {mean_flt_lar_prob:.3f} < ISS-T prob {mean_flt_isst_prob:.3f}")

        # Save results
        mission_results = {
            "task": "T2",
            "tissue": "liver",
            "mission": mission,
            "generated_at": datetime.now().isoformat(),
            "pca_distances": pca_distances,
            "pca_recovery_ratio": pca_recovery,
            "pca_variance_explained": pca.explained_variance_ratio_[:3].tolist(),
            "pathway_recovery": pathway_recovery,
            "classification_shift": clf_shift,
            "groups": {k: int(v.sum()) for k, v in groups.items()},
        }

        # Summary stats
        valid_recs = [v["recovery_fraction"] for v in pathway_recovery.values()
                      if v.get("recovery_fraction") is not None]
        if valid_recs:
            mission_results["pathway_recovery_summary"] = {
                "mean": np.mean(valid_recs),
                "median": np.median(valid_recs),
                "n_recovering": sum(1 for r in valid_recs if r > 0),
                "n_total": len(valid_recs),
            }

        save_json(mission_results, T_OUTPUT_DIR / f"T2_{mission}_liver_recovery.json")
        all_results[f"T2_{mission}"] = mission_results

    # ── PCA trajectory visualization ──
    _plot_t2_trajectory(all_results)

    return all_results


def _plot_t2_trajectory(results):
    """Generate PCA trajectory plot for T2."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  [SKIP] matplotlib not available for trajectory plot")
        return

    for key, res in results.items():
        mission = res.get("mission", "")
        meta = load_temporal_metadata("liver", mission=mission)
        gene_feat = load_gene_features("liver")

        try:
            feat_aligned, meta_aligned = align_features_with_meta(gene_feat, meta)
        except ValueError:
            continue

        from sklearn.preprocessing import StandardScaler
        from sklearn.decomposition import PCA

        X = feat_aligned.values.astype(float)
        X = np.nan_to_num(X, nan=0.0)
        var = np.var(X, axis=0)
        keep = var > np.percentile(var, VARIANCE_PERCENTILE * 100)
        X = X[:, keep]
        scaler = StandardScaler()
        X = scaler.fit_transform(X)
        pca = PCA(n_components=2, random_state=42)
        X_pca = pca.fit_transform(X)

        fig, ax = plt.subplots(1, 1, figsize=(8, 6))

        color_map = {
            "FLT_ISS-T": "#e74c3c", "FLT_LAR": "#f39c12",
            "GC_ISS-T": "#3498db", "GC_LAR": "#2ecc71",
            "BSL_ISS-T": "#9b59b6", "BSL_LAR": "#1abc9c",
        }
        marker_map = {
            "FLT_ISS-T": "^", "FLT_LAR": "v",
            "GC_ISS-T": "s", "GC_LAR": "D",
            "BSL_ISS-T": "o", "BSL_LAR": "p",
        }

        label_to_cond = {"Flight": "FLT", "FLT": "FLT", "GC": "GC", "BC": "BSL", "BSL": "BSL"}
        plotted = set()
        for i, (_, row) in enumerate(meta_aligned.iterrows()):
            cond = label_to_cond.get(row["label"], "?")
            timing = row["sacrifice_timing"]
            grp = f"{cond}_{timing}"
            if timing == "unknown":
                continue
            c = color_map.get(grp, "#7f8c8d")
            m = marker_map.get(grp, "x")
            label = grp if grp not in plotted else None
            plotted.add(grp)
            ax.scatter(X_pca[i, 0], X_pca[i, 1], c=c, marker=m, s=50,
                       alpha=0.7, label=label, edgecolors="white", linewidth=0.5)

        ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
        ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
        ax.set_title(f"T2 PCA Trajectory — {mission} Liver (ISS-T vs LAR)")
        ax.legend(fontsize=8, loc="best")
        ax.grid(alpha=0.3)

        outpath = T_OUTPUT_DIR / f"T2_{mission}_pca_trajectory.pdf"
        fig.tight_layout()
        fig.savefig(outpath, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {outpath}")


# ── T3: Age × Spaceflight Interaction ────────────────────────────────────────

def run_t3():
    """Run T3: Age × spaceflight interaction (RR-8 liver only)."""
    print("=" * 70)
    print("T3: Age × Spaceflight Interaction (RR-8 liver)")
    print("  OLD (32-week) vs YNG (10-12 week)")
    print("=" * 70)

    T_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    meta = load_temporal_metadata("liver", mission="RR-8")
    gene_feat = load_gene_features("liver")
    pathway_feat = load_pathway_features("liver")

    # Filter to known age
    meta_age = meta[meta["age_group"].isin(["OLD", "YNG"])].copy()
    print(f"\n  Total samples with age info: {len(meta_age)}")
    print(f"  Age distribution: {meta_age['age_group'].value_counts().to_dict()}")
    print(f"  Label distribution: {meta_age['label'].value_counts().to_dict()}")

    results = {
        "task": "T3",
        "tissue": "liver",
        "mission": "RR-8",
        "generated_at": datetime.now().isoformat(),
        "subtasks": {},
    }

    # ── T3a: Overall age classification ──
    print(f"\n  [T3a] Overall OLD vs YNG classification")

    for mode in ["gene", "pathway"]:
        feat_df = gene_feat if mode == "gene" else pathway_feat
        if feat_df is None:
            continue

        try:
            feat_aligned, meta_aligned = align_features_with_meta(feat_df, meta_age)
        except ValueError as e:
            print(f"    [{mode}] SKIP: {e}")
            continue

        X = feat_aligned.values.astype(float)
        X = np.nan_to_num(X, nan=0.0)
        y = (meta_aligned["age_group"] == "OLD").astype(int).values

        auroc, y_pred, y_proba = classify_binary_rskf(X, y, feature_mode=mode)
        ci_lo, ci_hi = bootstrap_auroc_ci(y, y_proba)
        perm_p = permutation_test_auroc(y, y_proba)

        print(f"    [{mode}] AUROC={auroc:.3f} [{ci_lo:.3f},{ci_hi:.3f}] p={perm_p:.4f} "
              f"(n={len(y)}: OLD={int((y==1).sum())}, YNG={int((y==0).sum())})")

        results["subtasks"][f"T3a_{mode}"] = {
            "description": "Overall OLD vs YNG classification",
            "auroc": auroc,
            "ci_low": ci_lo,
            "ci_high": ci_hi,
            "perm_p": perm_p,
            "n_old": int((y == 1).sum()),
            "n_yng": int((y == 0).sum()),
        }

    # ── T3b: Age within condition ──
    print(f"\n  [T3b] Age classification within condition groups")

    for condition, label_vals in [("FLT", ["Flight", "FLT"]),
                                   ("GC", ["GC"]),
                                   ("VIV", ["VC", "VIV"])]:
        mask_cond = meta_age["label"].isin(label_vals) & meta_age["age_group"].isin(["OLD", "YNG"])
        meta_sub = meta_age[mask_cond]

        if len(meta_sub) < 6:
            print(f"\n    [{condition}] SKIP: n={len(meta_sub)}")
            continue

        n_old = (meta_sub["age_group"] == "OLD").sum()
        n_yng = (meta_sub["age_group"] == "YNG").sum()
        print(f"\n    [{condition}] n={len(meta_sub)} (OLD={n_old}, YNG={n_yng})")

        for mode in ["gene", "pathway"]:
            feat_df = gene_feat if mode == "gene" else pathway_feat
            if feat_df is None:
                continue

            try:
                feat_aligned, meta_aligned = align_features_with_meta(feat_df, meta_sub)
            except ValueError:
                continue

            X = feat_aligned.values.astype(float)
            X = np.nan_to_num(X, nan=0.0)
            y = (meta_aligned["age_group"] == "OLD").astype(int).values

            min_class = min((y == 0).sum(), (y == 1).sum())
            if min_class >= 5:
                auroc, y_pred, y_proba = classify_binary_rskf(
                    X, y, feature_mode=mode,
                    n_splits=min(N_SPLITS, min_class), n_repeats=N_REPEATS
                )
                cv_method = "RSKF"
            else:
                auroc, y_pred, y_proba = classify_binary_loo(X, y, feature_mode=mode)
                cv_method = "LOO"

            ci_lo, ci_hi = bootstrap_auroc_ci(y, y_proba)
            perm_p = permutation_test_auroc(y, y_proba)

            print(f"      [{mode}] AUROC={auroc:.3f} [{ci_lo:.3f},{ci_hi:.3f}] p={perm_p:.4f} ({cv_method})")

            results["subtasks"][f"T3b_{condition}_{mode}"] = {
                "description": f"Age classification within {condition}",
                "condition": condition,
                "auroc": auroc,
                "ci_low": ci_lo,
                "ci_high": ci_hi,
                "perm_p": perm_p,
                "cv_method": cv_method,
                "n_old": int((y == 1).sum()),
                "n_yng": int((y == 0).sum()),
            }

    # ── T3c: Two-way ANOVA (pathway ~ flight + age + flight:age) ──
    print(f"\n  [T3c] Two-way ANOVA: pathway ~ flight + age + flight:age")
    print(f"        FLT + GC only, ISS-T only (preserve preservation control)")

    # Use ISS-T samples only to avoid preservation confound in ANOVA
    flt_gc_isst = meta_age[
        meta_age["label"].isin(["Flight", "FLT", "GC"]) &
        (meta_age["sacrifice_timing"] == "ISS-T") &
        meta_age["age_group"].isin(["OLD", "YNG"])
    ].copy()

    # If not enough ISS-T samples, use all timing
    if len(flt_gc_isst) < 20:
        print(f"    ISS-T only: n={len(flt_gc_isst)} (too few). Using all timing.")
        flt_gc_isst = meta_age[
            meta_age["label"].isin(["Flight", "FLT", "GC"]) &
            meta_age["age_group"].isin(["OLD", "YNG"])
        ].copy()

    # Create factors
    flt_gc_isst["flight"] = flt_gc_isst["label"].isin(["Flight", "FLT"]).astype(int)

    print(f"    ANOVA samples: n={len(flt_gc_isst)}")
    cross = pd.crosstab(
        flt_gc_isst["label"].map(lambda x: "FLT" if x in ["Flight", "FLT"] else "GC"),
        flt_gc_isst["age_group"]
    )
    print(f"    Cross-tabulation:\n{cross}")

    anova_results = {}
    if pathway_feat is not None:
        try:
            gsva_aligned, meta_anova = align_features_with_meta(pathway_feat, flt_gc_isst)
        except ValueError as e:
            print(f"    SKIP ANOVA: {e}")
            gsva_aligned = None

        if gsva_aligned is not None:
            import statsmodels.api as sm
            from statsmodels.formula.api import ols
            from statsmodels.stats.multitest import multipletests

            pathways = gsva_aligned.columns.tolist()
            anova_pvals = []

            for pw in pathways:
                df_anova = pd.DataFrame({
                    "score": gsva_aligned[pw].values,
                    "flight": meta_anova["flight"].values.astype(str),
                    "age": meta_anova["age_group"].values,
                })

                try:
                    model = ols("score ~ C(flight) + C(age) + C(flight):C(age)", data=df_anova).fit()
                    anova_table = sm.stats.anova_lm(model, typ=2)

                    p_flight = anova_table.loc["C(flight)", "PR(>F)"]
                    p_age = anova_table.loc["C(age)", "PR(>F)"]
                    p_interaction = anova_table.loc["C(flight):C(age)", "PR(>F)"]

                    anova_pvals.append({
                        "pathway": pw,
                        "p_flight": p_flight,
                        "p_age": p_age,
                        "p_interaction": p_interaction,
                    })
                except Exception as e:
                    print(f"    Warning: ANOVA failed for {pw}: {e}")

            if anova_pvals:
                # FDR correction (DD-19)
                interaction_ps = [r["p_interaction"] for r in anova_pvals]
                reject, fdr_pvals, _, _ = multipletests(interaction_ps, method="fdr_bh")

                for i, r in enumerate(anova_pvals):
                    r["fdr_interaction"] = fdr_pvals[i]
                    r["significant_interaction"] = bool(reject[i])

                n_sig = sum(reject)
                print(f"\n    ANOVA results: {len(anova_pvals)} pathways tested")
                print(f"    Significant interactions (FDR < 0.05): {n_sig}/{len(anova_pvals)}")

                if n_sig > 0:
                    sig_pws = [r for r in anova_pvals if r["significant_interaction"]]
                    sig_pws.sort(key=lambda x: x["fdr_interaction"])
                    print(f"\n    Significant interaction pathways:")
                    for r in sig_pws[:10]:
                        print(f"      {r['pathway']}: p_interaction={r['p_interaction']:.4f}, "
                              f"FDR={r['fdr_interaction']:.4f}")

                anova_results = {
                    "n_pathways": len(anova_pvals),
                    "n_significant_interaction": n_sig,
                    "n_samples": len(flt_gc_isst),
                    "timing_filter": "ISS-T" if "sacrifice_timing" in flt_gc_isst.columns and
                                     (flt_gc_isst["sacrifice_timing"] == "ISS-T").all() else "all",
                    "pathways": anova_pvals,
                }

    results["subtasks"]["T3c_anova"] = anova_results

    # ── T3d: Age-stratified spaceflight effect ──
    print(f"\n  [T3d] Age-stratified spaceflight AUROC")

    for age_grp in ["OLD", "YNG"]:
        mask_age = meta_age["age_group"] == age_grp
        mask_flt_gc = meta_age["label"].isin(["Flight", "FLT", "GC"])
        meta_sub = meta_age[mask_age & mask_flt_gc]

        if len(meta_sub) < 6:
            print(f"    [{age_grp}] SKIP: n={len(meta_sub)}")
            continue

        n_flt = meta_sub["label"].isin(["Flight", "FLT"]).sum()
        n_gc = (meta_sub["label"] == "GC").sum()
        print(f"\n    [{age_grp}] n={len(meta_sub)} (FLT={n_flt}, GC={n_gc})")

        for mode in ["gene", "pathway"]:
            feat_df = gene_feat if mode == "gene" else pathway_feat
            if feat_df is None:
                continue

            try:
                feat_aligned, meta_aligned = align_features_with_meta(feat_df, meta_sub)
            except ValueError:
                continue

            X = feat_aligned.values.astype(float)
            X = np.nan_to_num(X, nan=0.0)
            y = meta_aligned["label"].isin(["Flight", "FLT"]).astype(int).values

            min_class = min((y == 0).sum(), (y == 1).sum())
            if min_class >= 5:
                auroc, y_pred, y_proba = classify_binary_rskf(
                    X, y, feature_mode=mode,
                    n_splits=min(N_SPLITS, min_class), n_repeats=N_REPEATS
                )
            else:
                auroc, y_pred, y_proba = classify_binary_loo(X, y, feature_mode=mode)

            ci_lo, ci_hi = bootstrap_auroc_ci(y, y_proba)
            perm_p = permutation_test_auroc(y, y_proba)

            print(f"      [{mode}] AUROC={auroc:.3f} [{ci_lo:.3f},{ci_hi:.3f}] p={perm_p:.4f}")

            results["subtasks"][f"T3d_{age_grp}_{mode}"] = {
                "description": f"Spaceflight detection in {age_grp} mice",
                "age_group": age_grp,
                "auroc": auroc,
                "ci_low": ci_lo,
                "ci_high": ci_hi,
                "perm_p": perm_p,
                "n_flight": int((y == 1).sum()),
                "n_gc": int((y == 0).sum()),
            }

    # Interpretation: spaceflight amplifies aging?
    old_gene = results["subtasks"].get("T3d_OLD_gene", {})
    yng_gene = results["subtasks"].get("T3d_YNG_gene", {})
    if old_gene.get("auroc") and yng_gene.get("auroc"):
        delta = old_gene["auroc"] - yng_gene["auroc"]
        results["aging_hypothesis"] = {
            "OLD_auroc": old_gene["auroc"],
            "YNG_auroc": yng_gene["auroc"],
            "delta": delta,
            "interpretation": "Spaceflight amplifies age differences"
                              if delta > 0.05
                              else "No clear age-spaceflight amplification"
        }
        print(f"\n  Aging hypothesis: OLD AUROC={old_gene['auroc']:.3f}, "
              f"YNG AUROC={yng_gene['auroc']:.3f}, delta={delta:+.3f}")

    save_json(results, T_OUTPUT_DIR / f"T3_rr8_liver_results.json")
    return results


# ── Summary ──────────────────────────────────────────────────────────────────

def save_summary(t1_results=None, t2_results=None, t3_results=None):
    """Save combined T_temporal_summary.json, merging with existing data."""
    summary_path = RESULTS_DIR / "T_temporal_summary.json"

    # Load existing summary if present (merge, don't overwrite)
    existing = {}
    if summary_path.exists():
        try:
            existing = json.loads(summary_path.read_text())
        except (json.JSONDecodeError, OSError):
            existing = {}

    summary = {
        "generated_at": datetime.now().isoformat(),
        "version": "2.0",
    }

    # Preserve existing task results, update with new ones
    for task_key in ["T1", "T2", "T3"]:
        if task_key in existing:
            summary[task_key] = existing[task_key]

    if t1_results:
        summary["T1"] = t1_results
    if t2_results:
        summary["T2"] = t2_results
    if t3_results:
        summary["T3"] = t3_results

    save_json(summary, summary_path)


# ── CLI ────────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="Temporal & covariate analysis (v2.0)"
    )
    parser.add_argument(
        "--task", type=str, default=None,
        choices=["T1", "T2", "T3"],
        help="Task to run (T1/T2/T3)"
    )
    parser.add_argument(
        "--all", action="store_true",
        help="Run all tasks"
    )
    parser.add_argument(
        "--no-plot", action="store_true",
        help="Skip PCA trajectory plot"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    if not args.task and not args.all:
        print("Specify --task T1/T2/T3 or --all")
        return

    tasks = ["T1", "T2", "T3"] if args.all else [args.task]

    t1_results = None
    t2_results = None
    t3_results = None

    for task in tasks:
        if task == "T1":
            t1_results = run_t1()
        elif task == "T2":
            t2_results = run_t2()
        elif task == "T3":
            t3_results = run_t3()

    save_summary(t1_results, t2_results, t3_results)
    print("\n" + "=" * 70)
    print("Temporal analysis complete.")
    print(f"  Results: {T_OUTPUT_DIR}")
    print(f"  Summary: {RESULTS_DIR / 'T_temporal_summary.json'}")
    print("=" * 70)


if __name__ == "__main__":
    main()
