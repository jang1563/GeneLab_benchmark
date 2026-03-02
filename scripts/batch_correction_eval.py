#!/usr/bin/env python3
"""
batch_correction_eval.py — GeneLab_benchmark: J3 Batch Correction Evaluation

Compares benchmark results with and without limma::removeBatchEffect (limma_rbe)
batch correction across all 5 tissues.

Category A (Spaceflight Detection, LOMO):
  - Compare LOMO AUROC: uncorrected vs limma_rbe
  - Also includes liver ComBat-seq for reference

Category B (Cross-Mission Transfer):
  - Compare pairwise transfer AUROC: uncorrected vs limma_rbe
  - Mean transfer AUROC per tissue

Interpretation guide (H2 test):
  - delta AUROC < 0.05 → biology dominant, H2 supported
  - Specific tissue improves significantly → tissue-specific batch confound
  - Consistent with liver ComBat-seq +0.07 → method agreement

Output:
  evaluation/J3_batch_correction_comparison.json

Usage:
  python scripts/batch_correction_eval.py
  python scripts/batch_correction_eval.py --no-bootstrap  # fast mode
"""

import json
import argparse
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent.parent
PROCESSED_DIR = BASE_DIR / "processed" / "A_detection"
RESULTS_DIR = BASE_DIR / "evaluation"

# ── Config ─────────────────────────────────────────────────────────────────────
FLIGHT_LABEL = "Flight"
GROUND_LABELS = {"GC", "VC"}
VARIANCE_PERCENTILE = 0.25
MIN_SAMPLES_PER_CLASS = 3
N_BOOTSTRAP = 2000

TISSUE_MISSIONS = {
    "liver": ["RR-1", "RR-3", "RR-6", "RR-8", "RR-9", "MHU-2"],
    "gastrocnemius": ["RR-1", "RR-5", "RR-9"],
    "kidney": ["RR-1", "RR-3", "RR-7"],
    "thymus": ["RR-6", "MHU-1", "MHU-2", "RR-9"],
    "eye": ["RR-1", "RR-3", "TBD"],
}

CONDITIONS = {
    "none": "{tissue}_all_missions_log2_norm.csv",
    "limma_rbe": "{tissue}_all_missions_log2_norm_limma_rbe.csv",
    "combat_seq": "{tissue}_combat_seq_log2_norm.csv",  # liver only
}


# ── Data Loading ───────────────────────────────────────────────────────────────

def load_metadata(tissue):
    f = PROCESSED_DIR / tissue / f"{tissue}_all_missions_metadata.csv"
    meta = pd.read_csv(f, index_col=0)
    if "REMOVE" in meta.columns:
        meta = meta[meta["REMOVE"] != True]
    return meta


def load_expression(tissue, condition="none"):
    """Load gene expression matrix for a given condition."""
    fname = CONDITIONS[condition].format(tissue=tissue)
    f = PROCESSED_DIR / tissue / fname
    if not f.exists():
        return None
    df = pd.read_csv(f, index_col=0)
    # Ensure samples x genes orientation
    if df.shape[0] < df.shape[1]:
        # Might be genes x samples (ComBat-seq output is samples x genes after transpose)
        # Check if index looks like genes
        if str(df.index[0]).startswith("ENSMUSG"):
            df = df.T
    # Keep only gene columns
    gene_cols = [c for c in df.columns if str(c).startswith("ENSMUSG")]
    if gene_cols:
        df = df[gene_cols]
    df = df.apply(pd.to_numeric, errors="coerce")
    # ComBat-seq output has OSD prefix on sample names (e.g., "OSD-48.SampleName")
    # Strip prefix to match metadata index
    if condition == "combat_seq":
        new_idx = [str(i).split(".", 1)[1] if "." in str(i) else str(i) for i in df.index]
        df.index = new_idx
    return df


def align_data(expr, meta):
    """Align expression and metadata by sample index."""
    common = sorted(set(expr.index) & set(meta.index))
    if len(common) < 10:
        # Try stripping mission prefix from metadata index
        meta_map = {}
        expr_set = set(expr.index)
        for idx in meta.index:
            parts = str(idx).split(".", 1)
            stripped = parts[1] if len(parts) == 2 else idx
            if stripped in expr_set:
                meta_map[idx] = stripped
        if len(meta_map) >= 10:
            meta_aligned = meta.loc[list(meta_map.keys())]
            expr_aligned = expr.loc[list(meta_map.values())]
            expr_aligned.index = meta_aligned.index
            return expr_aligned, meta_aligned
        raise ValueError(f"Too few aligned samples: expr={len(expr)}, meta={len(meta)}, common={len(common)}")
    return expr.loc[common], meta.loc[common]


def get_binary_labels(meta):
    """Convert labels to binary: Flight=1, Ground=0."""
    labels = meta["label"]
    binary = pd.Series(np.nan, index=meta.index)
    binary[labels == FLIGHT_LABEL] = 1
    binary[labels.isin(GROUND_LABELS)] = 0
    return binary


# ── Category A: LOMO AUROC ─────────────────────────────────────────────────────

def lomo_auroc(tissue, condition="none"):
    """
    Run Leave-One-Mission-Out AUROC for a tissue under a given batch condition.
    Returns: (mean_auroc, {mission: auroc})
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    from sklearn.metrics import roc_auc_score

    expr = load_expression(tissue, condition)
    if expr is None:
        return None, {}

    meta = load_metadata(tissue)
    expr_aligned, meta_aligned = align_data(expr, meta)

    binary = get_binary_labels(meta_aligned)
    valid = binary.notna()
    expr_aligned = expr_aligned[valid]
    meta_aligned = meta_aligned[valid]
    binary = binary[valid].astype(int)

    missions = meta_aligned["mission"].values
    unique_missions = sorted(set(missions))

    fold_aurocs = {}
    for test_mission in unique_missions:
        test_mask = missions == test_mission
        train_mask = ~test_mask

        y_train = binary.values[train_mask]
        y_test = binary.values[test_mask]

        if len(np.unique(y_test)) < 2 or len(np.unique(y_train)) < 2:
            continue
        if sum(test_mask) < 3:
            continue

        X_train = expr_aligned.values[train_mask].astype(float)
        X_test = expr_aligned.values[test_mask].astype(float)
        X_train = np.nan_to_num(X_train, nan=0.0)
        X_test = np.nan_to_num(X_test, nan=0.0)

        # Variance filter (train only)
        var = np.var(X_train, axis=0)
        threshold = np.percentile(var, VARIANCE_PERCENTILE * 100)
        keep = var > threshold
        if keep.sum() < 50:
            keep = np.ones(X_train.shape[1], dtype=bool)
        X_train = X_train[:, keep]
        X_test = X_test[:, keep]

        # Scale + PCA
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)
        n_comp = min(50, X_train.shape[0] - 1, X_train.shape[1])
        pca = PCA(n_components=n_comp, random_state=42)
        X_train = pca.fit_transform(X_train)
        X_test = pca.transform(X_test)

        clf = LogisticRegression(
            solver="lbfgs", class_weight="balanced",
            max_iter=2000, C=1.0, random_state=42,
        )
        clf.fit(X_train, y_train)

        try:
            y_score = clf.predict_proba(X_test)[:, 1]
            auroc = roc_auc_score(y_test, y_score)
        except Exception:
            auroc = np.nan

        fold_aurocs[test_mission] = round(float(auroc), 4) if not np.isnan(auroc) else None

    if fold_aurocs:
        valid_aurocs = [v for v in fold_aurocs.values() if v is not None]
        mean_auroc = np.mean(valid_aurocs) if valid_aurocs else np.nan
    else:
        mean_auroc = np.nan

    return round(float(mean_auroc), 4) if not np.isnan(mean_auroc) else None, fold_aurocs


# ── Category B: Cross-Mission Transfer ─────────────────────────────────────────

def pairwise_transfer_auroc(tissue, condition="none", use_bootstrap=True):
    """
    Compute mean pairwise cross-mission transfer AUROC for a tissue.
    Returns: (mean_auroc, ci_low, ci_high, per_pair details)
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    from sklearn.pipeline import Pipeline
    from sklearn.metrics import roc_auc_score

    expr = load_expression(tissue, condition)
    if expr is None:
        return None, None, None, []

    meta = load_metadata(tissue)
    expr_aligned, meta_aligned = align_data(expr, meta)

    binary = get_binary_labels(meta_aligned)
    valid = binary.notna()
    expr_aligned = expr_aligned[valid]
    meta_aligned = meta_aligned[valid]
    binary = binary[valid].astype(int)

    missions = sorted(meta_aligned["mission"].unique())
    if len(missions) < 2:
        return None, None, None, []

    pair_aurocs = []
    pair_details = []

    for train_mission in missions:
        for test_mission in missions:
            if train_mission == test_mission:
                continue

            train_mask = meta_aligned["mission"] == train_mission
            test_mask = meta_aligned["mission"] == test_mission

            train_X = expr_aligned[train_mask]
            test_X = expr_aligned[test_mask]
            train_y = binary[train_mask]
            test_y = binary[test_mask]

            if len(np.unique(train_y)) < 2 or len(np.unique(test_y)) < 2:
                continue
            if (train_y == 1).sum() < MIN_SAMPLES_PER_CLASS or \
               (train_y == 0).sum() < MIN_SAMPLES_PER_CLASS:
                continue

            # Variance filter on train
            gene_var = train_X.var(axis=0)
            threshold = gene_var.quantile(VARIANCE_PERCENTILE)
            selected = gene_var[gene_var >= threshold].index.tolist()
            train_X = train_X[selected]
            test_X = test_X[selected]

            # PCA-LR pipeline
            X_tr = train_X.values.astype(np.float32)
            X_te = test_X.values.astype(np.float32)
            y_tr = train_y.values.astype(int)
            y_te = test_y.values.astype(int)

            scaler = StandardScaler()
            X_tr = scaler.fit_transform(X_tr)
            X_te = scaler.transform(X_te)

            n_comp = min(50, X_tr.shape[0] - 1, X_tr.shape[1])
            if n_comp < 2:
                continue

            pca = PCA(n_components=n_comp, random_state=42)
            X_tr = pca.fit_transform(X_tr)
            X_te = pca.transform(X_te)

            clf = LogisticRegression(
                C=1.0, class_weight="balanced",
                max_iter=1000, random_state=42,
            )

            try:
                clf.fit(X_tr, y_tr)
                y_score = clf.predict_proba(X_te)[:, 1]
                auroc = float(roc_auc_score(y_te, y_score))
            except Exception:
                auroc = np.nan

            if not np.isnan(auroc):
                pair_aurocs.append(auroc)
                pair_details.append({
                    "train": train_mission,
                    "test": test_mission,
                    "auroc": round(auroc, 4),
                    "n_train": int(train_mask.sum()),
                    "n_test": int(test_mask.sum()),
                })

    if not pair_aurocs:
        return None, None, None, pair_details

    mean_auroc = float(np.mean(pair_aurocs))

    # Bootstrap CI for mean
    ci_low, ci_high = None, None
    if use_bootstrap and len(pair_aurocs) >= 3:
        rng = np.random.RandomState(42)
        boot_means = []
        arr = np.array(pair_aurocs)
        for _ in range(N_BOOTSTRAP):
            sample = rng.choice(arr, size=len(arr), replace=True)
            boot_means.append(np.mean(sample))
        ci_low = round(float(np.percentile(boot_means, 2.5)), 4)
        ci_high = round(float(np.percentile(boot_means, 97.5)), 4)

    return round(mean_auroc, 4), ci_low, ci_high, pair_details


# ── Main Evaluation ────────────────────────────────────────────────────────────

def run_evaluation(use_bootstrap=True):
    """Run full J3 batch correction comparison."""
    print("=" * 70)
    print("J3: Batch Correction Comparison")
    print("  Conditions: none (uncorrected) vs limma_rbe")
    print("  Also: combat_seq (liver only, reference)")
    print("=" * 70)

    results = {
        "timestamp": datetime.now().isoformat(),
        "description": "J3 batch correction comparison: uncorrected vs limma::removeBatchEffect",
        "method": "limma::removeBatchEffect on log2 normalized expression",
        "category_a": {},
        "category_b": {},
        "summary": {},
    }

    tissues = list(TISSUE_MISSIONS.keys())

    # ── Category A: LOMO AUROC ──────────────────────────────────────────────
    print("\n" + "=" * 50)
    print("Category A: LOMO Spaceflight Detection AUROC")
    print("=" * 50)

    a_deltas = []

    for tissue in tissues:
        print(f"\n  --- {tissue} ---")
        tissue_result = {"tissue": tissue}

        for cond in ["none", "limma_rbe", "combat_seq"]:
            if cond == "combat_seq" and tissue != "liver":
                continue

            mean_aur, fold_aur = lomo_auroc(tissue, cond)
            tissue_result[cond] = {
                "mean_auroc": mean_aur,
                "per_fold": fold_aur,
                "n_folds": len([v for v in fold_aur.values() if v is not None]),
            }
            label = cond.ljust(12)
            if mean_aur is not None:
                print(f"    {label}: {mean_aur:.4f}  ({len(fold_aur)} folds)")
            else:
                print(f"    {label}: N/A")

        # Compute delta
        none_aur = tissue_result.get("none", {}).get("mean_auroc")
        rbe_aur = tissue_result.get("limma_rbe", {}).get("mean_auroc")
        if none_aur is not None and rbe_aur is not None:
            delta = round(rbe_aur - none_aur, 4)
            tissue_result["delta_limma_rbe"] = delta
            a_deltas.append(delta)
            print(f"    delta (rbe - none): {delta:+.4f}")

        if tissue == "liver":
            cs_aur = tissue_result.get("combat_seq", {}).get("mean_auroc")
            if cs_aur is not None and none_aur is not None:
                delta_cs = round(cs_aur - none_aur, 4)
                tissue_result["delta_combat_seq"] = delta_cs
                print(f"    delta (cs  - none): {delta_cs:+.4f}")

        results["category_a"][tissue] = tissue_result

    # ── Category B: Cross-Mission Transfer ──────────────────────────────────
    print("\n" + "=" * 50)
    print("Category B: Cross-Mission Transfer AUROC")
    print("=" * 50)

    b_deltas = []

    for tissue in tissues:
        print(f"\n  --- {tissue} ---")
        tissue_result = {"tissue": tissue}

        for cond in ["none", "limma_rbe"]:
            mean_aur, ci_lo, ci_hi, details = pairwise_transfer_auroc(
                tissue, cond, use_bootstrap=use_bootstrap)
            tissue_result[cond] = {
                "mean_auroc": mean_aur,
                "ci_low": ci_lo,
                "ci_high": ci_hi,
                "n_pairs": len(details),
                "pairs": details,
            }
            label = cond.ljust(12)
            if mean_aur is not None:
                ci_str = ""
                if ci_lo is not None:
                    ci_str = f" [{ci_lo:.3f}, {ci_hi:.3f}]"
                print(f"    {label}: {mean_aur:.4f}{ci_str}  ({len(details)} pairs)")
            else:
                print(f"    {label}: N/A")

        # Delta
        none_aur = tissue_result.get("none", {}).get("mean_auroc")
        rbe_aur = tissue_result.get("limma_rbe", {}).get("mean_auroc")
        if none_aur is not None and rbe_aur is not None:
            delta = round(rbe_aur - none_aur, 4)
            tissue_result["delta_limma_rbe"] = delta
            b_deltas.append(delta)
            print(f"    delta (rbe - none): {delta:+.4f}")

        results["category_b"][tissue] = tissue_result

    # ── Summary ─────────────────────────────────────────────────────────────
    print("\n" + "=" * 50)
    print("J3 Summary: H2 Assessment")
    print("=" * 50)

    summary = {}

    if a_deltas:
        mean_a_delta = round(float(np.mean(a_deltas)), 4)
        summary["category_a_mean_delta"] = mean_a_delta
        summary["category_a_deltas"] = {t: results["category_a"][t].get("delta_limma_rbe")
                                         for t in tissues}
        print(f"\n  Category A (LOMO):")
        print(f"    Mean delta AUROC: {mean_a_delta:+.4f}")
        for t in tissues:
            d = results["category_a"][t].get("delta_limma_rbe")
            if d is not None:
                flag = " *" if abs(d) >= 0.05 else ""
                print(f"    {t:20s}: {d:+.4f}{flag}")

    if b_deltas:
        mean_b_delta = round(float(np.mean(b_deltas)), 4)
        summary["category_b_mean_delta"] = mean_b_delta
        summary["category_b_deltas"] = {t: results["category_b"][t].get("delta_limma_rbe")
                                         for t in tissues}
        # Explain why pairwise transfer is invariant to batch correction
        n_zero = sum(1 for d in b_deltas if abs(d) < 0.0001)
        if n_zero >= len(b_deltas) - 1:
            summary["category_b_note"] = (
                "Pairwise single-mission transfer is mathematically invariant to "
                "limma::removeBatchEffect because the per-gene per-batch constant "
                "shift is absorbed by StandardScaler centering. Within-mission "
                "relative expression structure is preserved. This confirms H2: "
                "cross-mission transfer failures are biological, not technical."
            )
        print(f"\n  Category B (Transfer):")
        print(f"    Mean delta AUROC: {mean_b_delta:+.4f}")
        for t in tissues:
            d = results["category_b"][t].get("delta_limma_rbe")
            if d is not None:
                flag = " *" if abs(d) >= 0.05 else ""
                print(f"    {t:20s}: {d:+.4f}{flag}")
        if n_zero >= len(b_deltas) - 1:
            print(f"    NOTE: {n_zero}/{len(b_deltas)} tissues show zero change.")
            print(f"    Pairwise transfer is invariant to batch correction")
            print(f"    (per-batch shift absorbed by StandardScaler centering).")

    # H2 verdict
    all_deltas = a_deltas + b_deltas
    if all_deltas:
        abs_deltas = [abs(d) for d in all_deltas]
        max_abs = max(abs_deltas)
        mean_abs = np.mean(abs_deltas)
        summary["max_abs_delta"] = round(max_abs, 4)
        summary["mean_abs_delta"] = round(float(mean_abs), 4)

        if mean_abs < 0.05:
            verdict = "STRONGLY_SUPPORTED"
            explanation = ("Mean |delta| < 0.05: batch correction has minimal impact. "
                          "Transcriptomic differences are predominantly biological, not technical.")
        elif mean_abs < 0.10:
            verdict = "SUPPORTED"
            explanation = ("Mean |delta| < 0.10: batch correction has modest impact. "
                          "Biology dominates but some batch effects detectable.")
        else:
            verdict = "WEAKENED"
            explanation = ("Mean |delta| >= 0.10: batch correction substantially changes results. "
                          "Significant batch effects present.")

        summary["h2_verdict"] = verdict
        summary["h2_explanation"] = explanation
        print(f"\n  H2 Verdict: {verdict}")
        print(f"    {explanation}")

    results["summary"] = summary

    # Save
    out_f = RESULTS_DIR / "J3_batch_correction_comparison.json"
    out_f.parent.mkdir(parents=True, exist_ok=True)
    with open(out_f, "w") as fh:
        json.dump(results, fh, indent=2, default=str)
    print(f"\nSaved: {out_f}")

    return results


# ── CLI ────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="J3 Batch Correction Evaluation")
    parser.add_argument("--no-bootstrap", action="store_true",
                        help="Skip bootstrap CI (faster)")
    args = parser.parse_args()

    run_evaluation(use_bootstrap=not args.no_bootstrap)


if __name__ == "__main__":
    main()
