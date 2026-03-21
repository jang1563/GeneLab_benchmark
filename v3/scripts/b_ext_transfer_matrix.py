#!/usr/bin/env python3
"""
b_ext_transfer_matrix.py — v3 B_ext: Full cross-tissue transfer matrix.

Expands v1 Category C (4 predefined pairs) to full N×N matrix including
lung and colon from Phase 4.

Methods:
  A) Gene-level: Common genes → variance filter (train) → PCA → LR
  C) Pathway-level: fGSEA top-20 pathways → GSVA scores → LR
     (only for tissues with GSVA scores)

Usage:
    python v3/scripts/b_ext_transfer_matrix.py              # All pairs
    python v3/scripts/b_ext_transfer_matrix.py --fast        # Skip bootstrap/perm
    python v3/scripts/b_ext_transfer_matrix.py --pair liver kidney  # Single pair

Output:
    v3/evaluation/B_ext_transfer_matrix.json
"""

import argparse
import json
import sys
import warnings
from datetime import datetime
from itertools import permutations
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")

# ── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent.parent.parent  # GeneLab_benchmark/
PROCESSED_DIR = BASE_DIR / "processed" / "A_detection"
PATHWAY_DIR = BASE_DIR / "processed" / "pathway_scores"
FGSEA_DIR = BASE_DIR / "processed" / "fgsea"
OUTPUT_DIR = BASE_DIR / "v3" / "evaluation"

# ── Config ───────────────────────────────────────────────────────────────────
VARIANCE_PERCENTILE = 0.25   # DD-03: keep top 75% variance genes
N_BOOTSTRAP = 2000
N_PERMUTATIONS = 10000
CI_ALPHA = 0.05
TOP_N_PATHWAYS = 20

FLIGHT_LABELS = {"Flight"}
GROUND_LABELS = {"GC", "VC", "Ground"}  # v1 uses GC/VC, v3 new tissues use Ground

# ── Tissue definitions ───────────────────────────────────────────────────────
# Tissues with all_missions files (v1)
V1_TISSUES = {
    "liver":          {"missions": ["RR-1", "RR-3", "RR-6", "RR-8", "RR-9", "MHU-2"],
                       "has_gsva": True, "has_fgsea": True},
    "kidney":         {"missions": ["RR-1", "RR-3", "RR-7"],
                       "has_gsva": True, "has_fgsea": True},
    "thymus":         {"missions": ["RR-6", "MHU-1", "MHU-2", "RR-9"],
                       "has_gsva": True, "has_fgsea": True},
    "gastrocnemius":  {"missions": ["RR-1", "RR-5", "RR-9"],
                       "has_gsva": True, "has_fgsea": True},
    "eye":            {"missions": ["RR-1", "RR-3", "TBD"],
                       "has_gsva": True, "has_fgsea": True},
}

# New v3 tissues (per-mission files only)
V3_TISSUES = {
    "lung":  {"missions": ["RR-6"], "has_gsva": False, "has_fgsea": False},
    "colon": {"missions": ["RR-6"], "has_gsva": False, "has_fgsea": False},
}

ALL_TISSUES = {**V1_TISSUES, **V3_TISSUES}


# ── Data Loading ─────────────────────────────────────────────────────────────

def load_tissue_data(tissue):
    """Load log2 normalized counts + metadata for a tissue.

    For v1 tissues: use all_missions files.
    For v3 tissues: load per-mission files and concatenate.
    Returns (counts_df, meta_df) with only Flight/Ground samples.
    """
    tissue_dir = PROCESSED_DIR / tissue

    if tissue in V1_TISSUES:
        # Load all_missions files
        counts_path = tissue_dir / f"{tissue}_all_missions_log2_norm.csv"
        meta_path = tissue_dir / f"{tissue}_all_missions_metadata.csv"

        if not counts_path.exists() or not meta_path.exists():
            print(f"  [ERROR] all_missions data not found for {tissue}")
            return None, None

        counts = pd.read_csv(counts_path, index_col=0)
        meta = pd.read_csv(meta_path, index_col=0)

        # Drop non-gene columns from counts
        non_gene = [c for c in counts.columns
                    if c in {"mission", "osd_id", "label"} or
                    not str(c).startswith("ENSMUSG")]
        if non_gene:
            counts = counts.drop(columns=non_gene, errors="ignore")
        counts = counts.select_dtypes(include=[np.number])

        # Filter REMOVE rows
        if "REMOVE" in meta.columns:
            keep = meta["REMOVE"] != True
            meta = meta[keep]
            counts = counts.loc[counts.index.intersection(meta.index)]

    else:
        # v3 tissue: load per-mission files
        missions = ALL_TISSUES[tissue]["missions"]
        all_counts = []
        all_meta = []
        for mission in missions:
            cp = tissue_dir / f"{tissue}_{mission}_log2_norm.csv"
            mp = tissue_dir / f"{tissue}_{mission}_metadata.csv"
            if not cp.exists():
                print(f"  [WARN] Missing: {cp.name}")
                continue
            c = pd.read_csv(cp, index_col=0)
            m = pd.read_csv(mp, index_col=0)
            all_counts.append(c)
            all_meta.append(m)

        if not all_counts:
            print(f"  [ERROR] No data found for {tissue}")
            return None, None

        counts = pd.concat(all_counts)
        meta = pd.concat(all_meta)

        # Ensure counts has only gene columns
        gene_cols = [c for c in counts.columns if str(c).startswith("ENSMUSG")]
        non_gene = [c for c in counts.columns if not str(c).startswith("ENSMUSG")]
        if non_gene:
            counts = counts[gene_cols]
        counts = counts.apply(pd.to_numeric, errors="coerce")

    # Align counts and metadata
    common = counts.index.intersection(meta.index)
    counts = counts.loc[common]
    meta = meta.loc[common]

    # Extract binary labels
    labels = get_binary_labels(meta)
    valid = ~labels.isna()

    counts = counts[valid]
    meta = meta[valid]
    labels = labels[valid]

    n_flt = int((labels == 1).sum())
    n_gnd = int((labels == 0).sum())
    print(f"  {tissue}: {len(labels)} samples ({n_flt} Flight, {n_gnd} Ground)")

    return counts, labels


def get_binary_labels(meta):
    """Extract binary labels: 1=Flight, 0=Ground/GC/VC. NaN for others."""
    label_col = None
    for col in ["label", "label_raw", "condition", "group"]:
        if col in meta.columns:
            label_col = col
            break
    if label_col is None:
        raise ValueError(f"No label column. Columns: {list(meta.columns)}")

    labels_raw = meta[label_col]
    binary = pd.Series(np.nan, index=meta.index)
    binary[labels_raw.isin(FLIGHT_LABELS)] = 1
    binary[labels_raw.isin(GROUND_LABELS)] = 0
    return binary


def load_gsva_scores(tissue, db="hallmark"):
    """Load GSVA pathway scores for a tissue."""
    missions = ALL_TISSUES.get(tissue, {}).get("missions", [])
    all_scores = []

    for mission in missions:
        f = PATHWAY_DIR / tissue / f"{mission}_gsva_{db}.csv"
        if not f.exists():
            continue
        scores = pd.read_csv(f, index_col=0)
        scores["mission"] = mission
        all_scores.append(scores)

    if not all_scores:
        return None
    combined = pd.concat(all_scores)
    # Remove duplicates (e.g., MHU-1/MHU-2 may share sample names)
    combined = combined[~combined.index.duplicated(keep="first")]
    return combined


def load_fgsea_results(tissue, db="hallmark"):
    """Load fGSEA results across missions."""
    missions = ALL_TISSUES.get(tissue, {}).get("missions", [])
    dfs = []
    for mission in missions:
        f = FGSEA_DIR / tissue / f"{mission}_fgsea_{db}.csv"
        if not f.exists():
            continue
        dfs.append(pd.read_csv(f))
    if not dfs:
        return None
    return pd.concat(dfs, ignore_index=True)


# ── Statistical utilities ────────────────────────────────────────────────────

def bootstrap_auroc_ci(y_true, y_score, n_boot=None):
    """Bootstrap 95% CI for AUROC."""
    from sklearn.metrics import roc_auc_score

    if n_boot is None:
        n_boot = N_BOOTSTRAP
    if n_boot <= 0:
        return np.nan, np.nan

    rng = np.random.RandomState(42)
    n = len(y_true)
    if n < 4 or len(np.unique(y_true)) < 2:
        return np.nan, np.nan

    boots = []
    for _ in range(n_boot):
        idx = rng.randint(0, n, size=n)
        yt, ys = y_true[idx], y_score[idx]
        if len(np.unique(yt)) < 2:
            continue
        try:
            boots.append(roc_auc_score(yt, ys))
        except Exception:
            continue

    if len(boots) < n_boot * 0.5:
        return np.nan, np.nan
    return float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))


def permutation_test(y_true, y_score, n_perm=None):
    """Permutation p-value for AUROC."""
    from sklearn.metrics import roc_auc_score

    if n_perm is None:
        n_perm = N_PERMUTATIONS
    if n_perm <= 0:
        return np.nan

    rng = np.random.RandomState(42)
    y_true = np.array(y_true)
    y_score = np.array(y_score)

    if len(np.unique(y_true)) < 2:
        return np.nan
    try:
        observed = roc_auc_score(y_true, y_score)
    except Exception:
        return np.nan

    count = 0
    for _ in range(n_perm):
        perm_y = rng.permutation(y_true)
        if len(np.unique(perm_y)) < 2:
            continue
        try:
            if roc_auc_score(perm_y, y_score) >= observed:
                count += 1
        except Exception:
            continue

    return (count + 1) / (n_perm + 1)


# ── Transfer Methods ─────────────────────────────────────────────────────────

def method_a_transfer(train_counts, train_labels, test_counts, test_labels):
    """Method A: Gene-level direct transfer.

    Common genes → variance filter (train) → StandardScaler → PCA → LR.
    """
    from sklearn.decomposition import PCA
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    _nan = {"auroc": np.nan, "y_true": None, "y_score": None}

    common_genes = train_counts.columns.intersection(test_counts.columns)
    if len(common_genes) < 100:
        print(f"    [WARN] Only {len(common_genes)} common genes")
        if len(common_genes) < 10:
            return _nan

    X_train = train_counts[common_genes].copy()
    X_test = test_counts[common_genes].copy()

    # Variance filter on train (DD-03)
    gene_var = X_train.var(axis=0)
    threshold = gene_var.quantile(VARIANCE_PERCENTILE)
    keep = gene_var[gene_var >= threshold].index.tolist()
    X_train = X_train[keep]
    X_test = X_test[keep]

    y_train = train_labels.values.astype(int)
    y_test = test_labels.values.astype(int)

    if len(np.unique(y_train)) < 2 or len(np.unique(y_test)) < 2:
        return _nan

    scaler = StandardScaler()
    X_tr = scaler.fit_transform(X_train.values.astype(np.float32))
    X_te = scaler.transform(X_test.values.astype(np.float32))

    n_comps = min(50, X_tr.shape[0] - 1, X_tr.shape[1])
    if n_comps < 2:
        return _nan

    model = Pipeline([
        ("pca", PCA(n_components=n_comps, random_state=42)),
        ("clf", LogisticRegression(C=1.0, class_weight="balanced",
                                   max_iter=1000, random_state=42))
    ])

    try:
        model.fit(X_tr, y_train)
        y_score = model.predict_proba(X_te)[:, 1]
        auroc = float(roc_auc_score(y_test, y_score))
        return {"auroc": auroc, "y_true": y_test, "y_score": y_score,
                "n_features": len(keep), "n_common_genes": len(common_genes)}
    except Exception as e:
        print(f"    [ERROR] Method A: {e}")
        return _nan


def select_top_pathways(tissue, db="hallmark", top_n=TOP_N_PATHWAYS):
    """Select top-N pathways by |mean NES| from fGSEA results."""
    fgsea = load_fgsea_results(tissue, db)
    if fgsea is None:
        return [], None

    agg = (fgsea.groupby("pathway")
           .agg(mean_NES=("NES", "mean"),
                abs_mean_NES=("NES", lambda x: abs(x.mean())),
                n_missions=("mission", "nunique"),
                n_sig=("padj", lambda x: (x < 0.05).sum()))
           .sort_values("abs_mean_NES", ascending=False))

    selected = agg.head(top_n).index.tolist()
    return selected, agg.loc[selected]


def _align_samples(gsva_index, label_index):
    """Align GSVA sample names with metadata (handles mission-prefix mismatch)."""
    gsva_set = set(gsva_index)
    label_set = set(label_index)

    direct = gsva_set & label_set
    if len(direct) >= 5:
        direct = sorted(direct)
        return direct, direct

    # Strip mission prefix from label index ("RR-1.Sample" → "Sample")
    gsva_list, label_list = [], []
    for lab in label_index:
        parts = str(lab).split(".", 1)
        stripped = parts[1] if len(parts) == 2 else lab
        if stripped in gsva_set:
            gsva_list.append(stripped)
            label_list.append(lab)

    if len(gsva_list) >= 5:
        return gsva_list, label_list

    # Try reverse
    for g in gsva_index:
        parts = str(g).split(".", 1)
        if len(parts) == 2 and parts[1] in label_set:
            gsva_list.append(g)
            label_list.append(parts[1])

    return gsva_list, label_list


def method_c_transfer(train_tissue, test_tissue,
                      train_labels, test_labels, db="hallmark"):
    """Method C: Pathway-level transfer.

    Select top-N pathways from train tissue fGSEA → GSVA scores → LR.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.preprocessing import StandardScaler

    _nan = {"auroc": np.nan, "y_true": None, "y_score": None}

    # Check GSVA availability
    if not ALL_TISSUES[train_tissue].get("has_gsva"):
        return _nan
    if not ALL_TISSUES[test_tissue].get("has_gsva"):
        return _nan

    # Select pathways from train tissue
    selected, pw_stats = select_top_pathways(train_tissue, db)
    if not selected:
        print("    [ERROR] No pathways selected")
        return _nan

    # Load GSVA scores
    train_gsva = load_gsva_scores(train_tissue, db)
    test_gsva = load_gsva_scores(test_tissue, db)
    if train_gsva is None or test_gsva is None:
        return _nan

    train_feat = train_gsva.drop(columns=["mission"], errors="ignore")
    test_feat = test_gsva.drop(columns=["mission"], errors="ignore")

    common_pw = [p for p in selected
                 if p in train_feat.columns and p in test_feat.columns]
    if len(common_pw) < 3:
        print(f"    [ERROR] Only {len(common_pw)} common pathways")
        return _nan

    # Align samples
    tr_gsva_idx, tr_label_idx = _align_samples(train_feat.index, train_labels.index)
    te_gsva_idx, te_label_idx = _align_samples(test_feat.index, test_labels.index)

    if len(tr_gsva_idx) < 5 or len(te_gsva_idx) < 5:
        print(f"    [ERROR] Too few aligned: train={len(tr_gsva_idx)}, test={len(te_gsva_idx)}")
        return _nan

    X_train = train_feat.loc[tr_gsva_idx, common_pw].values.astype(np.float32)
    X_test = test_feat.loc[te_gsva_idx, common_pw].values.astype(np.float32)
    y_train = train_labels.loc[tr_label_idx].values.astype(int)
    y_test = test_labels.loc[te_label_idx].values.astype(int)

    # Safety: ensure X and y have same length after alignment
    if X_train.shape[0] != len(y_train) or X_test.shape[0] != len(y_test):
        print(f"    [ERROR] Length mismatch after alignment: "
              f"X_train={X_train.shape[0]} y_train={len(y_train)} "
              f"X_test={X_test.shape[0]} y_test={len(y_test)}")
        return _nan

    if len(np.unique(y_train)) < 2 or len(np.unique(y_test)) < 2:
        return _nan

    scaler = StandardScaler()
    X_tr = scaler.fit_transform(X_train)
    X_te = scaler.transform(X_test)

    model = LogisticRegression(C=1.0, class_weight="balanced",
                               max_iter=1000, random_state=42)

    try:
        model.fit(X_tr, y_train)
        y_score = model.predict_proba(X_te)[:, 1]
        auroc = float(roc_auc_score(y_test, y_score))

        coef = model.coef_[0]
        importance = sorted(zip(common_pw, coef.tolist()),
                           key=lambda x: abs(x[1]), reverse=True)

        return {"auroc": auroc, "y_true": y_test, "y_score": y_score,
                "n_features": len(common_pw),
                "selected_pathways": common_pw,
                "top_features": importance[:10]}
    except Exception as e:
        print(f"    [ERROR] Method C: {e}")
        return _nan


# ── Main evaluation ──────────────────────────────────────────────────────────

def evaluate_pair(train_tissue, test_tissue, fast=False):
    """Evaluate transfer from train_tissue → test_tissue using Methods A and C."""
    print(f"\n{'─'*60}")
    print(f"  {train_tissue} → {test_tissue}")
    print(f"{'─'*60}")

    train_counts, train_labels = load_tissue_data(train_tissue)
    test_counts, test_labels = load_tissue_data(test_tissue)

    if train_counts is None or test_counts is None:
        return None

    result = {
        "train_tissue": train_tissue,
        "test_tissue": test_tissue,
        "n_train": len(train_labels),
        "n_test": len(test_labels),
    }

    # Method A: Gene-level
    print(f"  Method A (gene-level):")
    res_a = method_a_transfer(train_counts, train_labels,
                              test_counts, test_labels)
    if not np.isnan(res_a.get("auroc", np.nan)):
        ci_low, ci_high = bootstrap_auroc_ci(
            res_a["y_true"], res_a["y_score"],
            n_boot=0 if fast else N_BOOTSTRAP)
        perm_p = permutation_test(
            res_a["y_true"], res_a["y_score"],
            n_perm=0 if fast else N_PERMUTATIONS)
        result["method_a"] = {
            "auroc": res_a["auroc"],
            "ci_low": None if np.isnan(ci_low) else ci_low,
            "ci_high": None if np.isnan(ci_high) else ci_high,
            "perm_p": None if np.isnan(perm_p) else perm_p,
            "n_features": res_a.get("n_features"),
            "n_common_genes": res_a.get("n_common_genes"),
        }
        print(f"    AUROC = {res_a['auroc']:.3f}")
    else:
        result["method_a"] = {"auroc": None}
        print(f"    AUROC = N/A")

    # Method C: Pathway-level
    print(f"  Method C (pathway-level):")
    res_c = method_c_transfer(train_tissue, test_tissue,
                              train_labels, test_labels)
    if not np.isnan(res_c.get("auroc", np.nan)):
        ci_low, ci_high = bootstrap_auroc_ci(
            res_c["y_true"], res_c["y_score"],
            n_boot=0 if fast else N_BOOTSTRAP)
        perm_p = permutation_test(
            res_c["y_true"], res_c["y_score"],
            n_perm=0 if fast else N_PERMUTATIONS)
        result["method_c"] = {
            "auroc": res_c["auroc"],
            "ci_low": None if np.isnan(ci_low) else ci_low,
            "ci_high": None if np.isnan(ci_high) else ci_high,
            "perm_p": None if np.isnan(perm_p) else perm_p,
            "n_features": res_c.get("n_features"),
            "top_features": res_c.get("top_features"),
        }
        print(f"    AUROC = {res_c['auroc']:.3f}")
    else:
        result["method_c"] = {"auroc": None}
        has_gsva = (ALL_TISSUES[train_tissue].get("has_gsva") and
                    ALL_TISSUES[test_tissue].get("has_gsva"))
        reason = "no GSVA" if not has_gsva else "failed"
        print(f"    AUROC = N/A ({reason})")

    return result


def build_matrix(tissue_list, all_pairs):
    """Build AUROC matrices from pair results."""
    n = len(tissue_list)
    idx = {t: i for i, t in enumerate(tissue_list)}
    mat_a = [[None] * n for _ in range(n)]
    mat_c = [[None] * n for _ in range(n)]

    for p in all_pairs:
        i = idx[p["train_tissue"]]
        j = idx[p["test_tissue"]]
        mat_a[i][j] = p.get("method_a", {}).get("auroc")
        mat_c[i][j] = p.get("method_c", {}).get("auroc")

    return mat_a, mat_c


def main():
    global N_BOOTSTRAP, N_PERMUTATIONS

    parser = argparse.ArgumentParser(
        description="v3 B_ext: Full cross-tissue transfer matrix")
    parser.add_argument("--fast", action="store_true",
                        help="Skip bootstrap/permutation (quick check)")
    parser.add_argument("--pair", nargs=2, metavar=("TRAIN", "TEST"),
                        help="Run single pair")
    parser.add_argument("--tissues", nargs="+",
                        help="Subset of tissues to evaluate")
    args = parser.parse_args()

    if args.fast:
        N_BOOTSTRAP = 0
        N_PERMUTATIONS = 0

    print("=" * 70)
    print("GeneLabBench v3 — B_ext: Cross-Tissue Transfer Matrix")
    print(f"Timestamp: {datetime.now().isoformat()}")
    print("=" * 70)

    # Determine tissues
    if args.tissues:
        tissues = [t for t in args.tissues if t in ALL_TISSUES]
    else:
        tissues = list(ALL_TISSUES.keys())
    print(f"\nTissues: {tissues}")

    # Generate pairs
    if args.pair:
        pairs_to_run = [(args.pair[0], args.pair[1])]
    else:
        pairs_to_run = list(permutations(tissues, 2))
    print(f"Total pairs: {len(pairs_to_run)}")

    # Run evaluations
    all_results = []
    for i, (train_t, test_t) in enumerate(pairs_to_run):
        print(f"\n[{i+1}/{len(pairs_to_run)}]", end="")
        result = evaluate_pair(train_t, test_t, fast=args.fast)
        if result:
            all_results.append(result)

    # Build matrices
    mat_a, mat_c = build_matrix(tissues, all_results)

    # Summary
    print(f"\n{'='*70}")
    print("TRANSFER MATRIX SUMMARY")
    print(f"{'='*70}")

    # Method A summary
    a_aurocs = [r["method_a"]["auroc"] for r in all_results
                if r["method_a"].get("auroc") is not None]
    if a_aurocs:
        print(f"\nMethod A (gene-level): {len(a_aurocs)} pairs evaluated")
        print(f"  Mean AUROC: {np.mean(a_aurocs):.3f}")
        print(f"  Range: [{min(a_aurocs):.3f}, {max(a_aurocs):.3f}]")

    # Method C summary
    c_aurocs = [r["method_c"]["auroc"] for r in all_results
                if r["method_c"].get("auroc") is not None]
    if c_aurocs:
        print(f"\nMethod C (pathway-level): {len(c_aurocs)} pairs evaluated")
        print(f"  Mean AUROC: {np.mean(c_aurocs):.3f}")
        print(f"  Range: [{min(c_aurocs):.3f}, {max(c_aurocs):.3f}]")

    # H3 comparison (pathway vs gene)
    both = [(r["method_a"]["auroc"], r["method_c"]["auroc"])
            for r in all_results
            if r["method_a"].get("auroc") is not None
            and r["method_c"].get("auroc") is not None]
    if both:
        c_wins = sum(1 for a, c in both if c > a)
        mean_diff = np.mean([c - a for a, c in both])
        print(f"\nH3 Test: Method C > A in {c_wins}/{len(both)} pairs")
        print(f"  Mean AUROC diff (C−A): {mean_diff:+.3f}")

    # Save
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output = {
        "tissues": tissues,
        "n_pairs": len(all_results),
        "method_a_matrix": mat_a,
        "method_c_matrix": mat_c,
        "pairs": all_results,
        "summary": {
            "method_a_mean": float(np.mean(a_aurocs)) if a_aurocs else None,
            "method_c_mean": float(np.mean(c_aurocs)) if c_aurocs else None,
            "h3_c_wins": c_wins if both else None,
            "h3_c_total": len(both) if both else None,
            "h3_mean_diff": float(mean_diff) if both else None,
        },
        "config": {
            "variance_percentile": VARIANCE_PERCENTILE,
            "n_bootstrap": N_BOOTSTRAP,
            "n_permutations": N_PERMUTATIONS,
            "top_n_pathways": TOP_N_PATHWAYS,
        },
        "timestamp": datetime.now().isoformat(),
    }

    out_path = OUTPUT_DIR / "B_ext_transfer_matrix.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, default=_json_default)
    print(f"\nSaved: {out_path}")


def _json_default(obj):
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, pd.DataFrame):
        return obj.to_dict()
    return str(obj)


if __name__ == "__main__":
    main()
