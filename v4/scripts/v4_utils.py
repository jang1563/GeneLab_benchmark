#!/usr/bin/env python3
"""
v4_utils.py — GeneLabBench v4: Extended data loading utilities

Extends scripts/utils.py with:
- 8-tissue TISSUE_MISSIONS (was 5 in v1)
- LOMO fold creation from raw data (for tissues without pre-split tasks/)
- Stratified k-fold for single-mission tissues (lung, colon)
- Bootstrap CI (n=2000) and permutation p-value (n=1000, DD-32)
"""

import json
import time
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent.parent.parent  # v4/scripts/ → project root
PROCESSED_DIR = BASE_DIR / "processed" / "A_detection"
TASKS_DIR = BASE_DIR / "tasks"
V4_EVAL_DIR = BASE_DIR / "v4" / "evaluation"

# ── Tissue Configuration (corrected from HPC metadata inspection) ─────────────
TISSUE_MISSIONS = {
    # v1 original (5 tissues)
    "liver":          ["RR-1", "RR-3", "RR-6", "RR-8", "RR-9", "MHU-2"],
    "gastrocnemius":  ["RR-1", "RR-5", "RR-9"],
    "kidney":         ["RR-1", "RR-3", "RR-7"],
    "thymus":         ["RR-6", "MHU-1", "MHU-2", "RR-9"],
    "eye":            ["RR-1", "RR-3", "TBD"],       # TBD = OSD-397 (16 samples, 9F+7GC)
    # v3/v4 additions (3 tissues)
    "skin":           ["MHU-2", "RR-6", "RR-7"],     # Corrected from plan's RR-1/RR-5/RR-6
    "lung":           ["RR-6"],                        # Single mission → 5-fold CV
    "colon":          ["RR-6"],                        # Single mission → 5-fold CV
}

# Tissues with pre-existing task/ fold splits (can reuse directly)
TASK_MAP = {
    "liver":         "A1_liver_lomo",
    "gastrocnemius": "A2_gastrocnemius_lomo",
    "kidney":        "A3_kidney_lomo",
    "thymus":        "A4_thymus_lomo",
    "skin":          "A5_skin_lomo",
    "eye":           "A6_eye_lomo",
    # lung, colon: no pre-existing task splits
}

# CV strategy per tissue
LOMO_TISSUES = ["liver", "gastrocnemius", "kidney", "thymus", "eye", "skin"]
KFOLD_TISSUES = ["lung", "colon"]  # Single-mission → 5-fold stratified CV (DD-30)

# Binary label encoding: Flight=1, all standard controls=0
# Excludes: AG (Artificial Gravity, MHU-2 third condition)
LABEL_MAP = {
    "Flight": 1, "FLT": 1,
    "GC": 0, "Ground Control": 0, "Ground": 0,
    "Vivarium": 0, "Basal": 0,
    "BC": 0, "VC": 0,  # Basal Control, Vivarium Control
}
# Labels to EXCLUDE from binary classification (not Flight, not standard control)
EXCLUDE_LABELS = {"AG"}  # Artificial Gravity (MHU-2)


# ── Label Encoding ────────────────────────────────────────────────────────────

def encode_labels(meta):
    """Encode labels to binary (Flight=1, Control=0), filtering excluded groups.

    Returns (y, valid_mask): y is integer array (0/1), valid_mask is boolean mask
    for samples with valid binary labels. Samples with labels in EXCLUDE_LABELS
    or unmapped labels are filtered out.
    """
    y = meta["label"].map(LABEL_MAP)
    # Warn about unmapped labels (not in LABEL_MAP and not in EXCLUDE_LABELS)
    unmapped = meta["label"][y.isna() & ~meta["label"].isin(EXCLUDE_LABELS)]
    if len(unmapped) > 0:
        unique_unmapped = unmapped.unique()
        warnings.warn(f"Unmapped labels (dropped): {list(unique_unmapped)}")
    valid_mask = y.notna()
    y = y.fillna(-1).astype(int)  # -1 for invalid, masked out by valid_mask
    return y, valid_mask


# ── Data Loading ──────────────────────────────────────────────────────────────

def load_metadata(tissue):
    """Load all-missions metadata for a tissue. Handles per-mission files if no all_missions."""
    # Try all-missions file first
    f = PROCESSED_DIR / tissue / f"{tissue}_all_missions_metadata.csv"
    if f.exists():
        meta = pd.read_csv(f, index_col=0)
        if "REMOVE" in meta.columns:
            meta = meta[meta["REMOVE"] != True]
        return meta

    # Fallback: single-mission file
    missions = TISSUE_MISSIONS.get(tissue, [])
    if len(missions) == 1:
        mission = missions[0]
        f = PROCESSED_DIR / tissue / f"{tissue}_{mission}_metadata.csv"
        if f.exists():
            meta = pd.read_csv(f, index_col=0)
            if "REMOVE" in meta.columns:
                meta = meta[meta["REMOVE"] != True]
            if "mission" not in meta.columns:
                meta["mission"] = mission
            return meta

    raise FileNotFoundError(f"No metadata found for tissue={tissue}")


def load_gene_features(tissue):
    """Load gene-level log2 normalized counts (samples × genes)."""
    # Try all-missions file first
    f = PROCESSED_DIR / tissue / f"{tissue}_all_missions_log2_norm.csv"
    if not f.exists():
        missions = TISSUE_MISSIONS.get(tissue, [])
        if len(missions) == 1:
            f = PROCESSED_DIR / tissue / f"{tissue}_{missions[0]}_log2_norm.csv"

    if not f.exists():
        raise FileNotFoundError(f"No gene features found for tissue={tissue}")

    df = pd.read_csv(f, index_col=0)
    # Auto-orient: genes should be columns (more columns than rows for gene data)
    if df.shape[0] > df.shape[1]:
        df = df.T
    # Drop non-gene columns
    meta_cols = [c for c in df.columns if not str(c).startswith("ENSMUSG")]
    if meta_cols:
        df = df.drop(columns=meta_cols)
    df = df.apply(pd.to_numeric, errors="coerce")
    return df


def align_features_with_meta(features, meta):
    """Align feature matrix with metadata by sample name."""
    feat_set = set(features.index)
    meta_set = set(meta.index)
    common = sorted(feat_set & meta_set)

    if len(common) >= 5:
        return features.loc[common], meta.loc[common]

    # Try stripping mission prefix
    meta_map = {}
    for idx in meta.index:
        parts = str(idx).split(".", 1)
        stripped = parts[1] if len(parts) == 2 else idx
        if stripped in feat_set:
            meta_map[idx] = stripped

    if len(meta_map) >= 5:
        meta_aligned = meta.loc[list(meta_map.keys())]
        feat_aligned = features.loc[list(meta_map.values())]
        feat_aligned.index = meta_aligned.index
        return feat_aligned, meta_aligned

    raise ValueError(
        f"Too few aligned samples: features={len(feat_set)}, "
        f"meta={len(meta_set)}, common={len(common)}"
    )


# ── Fold Generation ──────────────────────────────────────────────────────────

def get_lomo_folds(tissue, gene_features=None, meta=None):
    """Generate LOMO fold splits for multi-mission tissues.

    Returns list of dicts: [{train_X, train_y, test_X, test_y, test_mission, fold_name}, ...]
    """
    if gene_features is None:
        gene_features = load_gene_features(tissue)
    if meta is None:
        meta = load_metadata(tissue)
    gene_features, meta = align_features_with_meta(gene_features, meta)

    # Encode labels and filter to valid binary samples
    y_all, valid_mask = encode_labels(meta)

    # Apply valid mask to both features and metadata
    gene_features = gene_features[valid_mask]
    meta = meta[valid_mask]
    y_all = y_all[valid_mask]

    missions = TISSUE_MISSIONS[tissue]
    folds = []

    for test_mission in missions:
        test_mask = meta["mission"] == test_mission
        train_mask = ~test_mask

        if test_mask.sum() < 2 or train_mask.sum() < 5:
            continue

        # Skip if test set has only one class
        y_test = y_all[test_mask]
        if y_test.nunique() < 2:
            continue

        folds.append({
            "train_X": gene_features[train_mask].values.astype(np.float32),
            "train_y": y_all[train_mask].values.astype(int),
            "test_X": gene_features[test_mask].values.astype(np.float32),
            "test_y": y_test.values.astype(int),
            "test_mission": test_mission,
            "fold_name": f"fold_{test_mission}_test",
            "n_train": int(train_mask.sum()),
            "n_test": int(test_mask.sum()),
            "gene_names": list(gene_features.columns),
        })

    return folds


def get_kfold_splits(tissue, n_splits=5, seed=42, gene_features=None, meta=None):
    """Generate stratified k-fold splits for single-mission tissues.

    Returns list of dicts matching LOMO fold format.
    """
    if gene_features is None:
        gene_features = load_gene_features(tissue)
    if meta is None:
        meta = load_metadata(tissue)
    gene_features, meta = align_features_with_meta(gene_features, meta)

    # Encode labels and filter to valid binary samples
    y_encoded, valid_mask = encode_labels(meta)
    gene_features = gene_features[valid_mask]
    meta = meta[valid_mask]

    y_all = y_encoded[valid_mask].values.astype(int)
    X_all = gene_features.values.astype(np.float32)

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    folds = []

    for i, (train_idx, test_idx) in enumerate(skf.split(X_all, y_all)):
        folds.append({
            "train_X": X_all[train_idx],
            "train_y": y_all[train_idx],
            "test_X": X_all[test_idx],
            "test_y": y_all[test_idx],
            "test_mission": f"fold_{i+1}",
            "fold_name": f"kfold_{i+1}",
            "n_train": len(train_idx),
            "n_test": len(test_idx),
            "gene_names": list(gene_features.columns),
        })

    return folds


def get_folds(tissue, **kwargs):
    """Unified fold generation: LOMO for multi-mission, k-fold for single-mission."""
    if tissue in LOMO_TISSUES:
        return get_lomo_folds(tissue, **kwargs)
    elif tissue in KFOLD_TISSUES:
        return get_kfold_splits(tissue, **kwargs)
    else:
        raise ValueError(f"Unknown tissue: {tissue}")


def get_folds_from_task_dir(tissue):
    """Load pre-existing LOMO fold splits from tasks/ directory.

    Returns list of dicts matching the fold format.
    Only available for tissues in TASK_MAP (original v1 tissues + skin).
    """
    if tissue not in TASK_MAP:
        raise ValueError(f"No pre-existing task for tissue={tissue}. Use get_folds() instead.")

    task_dir = TASKS_DIR / TASK_MAP[tissue]
    if not task_dir.exists():
        raise FileNotFoundError(f"Task directory not found: {task_dir}")

    fold_dirs = sorted([d for d in task_dir.iterdir()
                        if d.is_dir() and d.name.startswith("fold_")])
    folds = []

    for fold_dir in fold_dirs:
        train_X = pd.read_csv(fold_dir / "train_X.csv", index_col=0)
        test_X = pd.read_csv(fold_dir / "test_X.csv", index_col=0)
        train_y = pd.read_csv(fold_dir / "train_y.csv", index_col=0).squeeze()
        test_y = pd.read_csv(fold_dir / "test_y.csv", index_col=0).squeeze()

        # Align features
        common_genes = train_X.columns.intersection(test_X.columns)
        train_X = train_X[common_genes]
        test_X = test_X[common_genes]

        # Read fold info
        fold_info_path = fold_dir / "fold_info.json"
        fold_info = json.loads(fold_info_path.read_text()) if fold_info_path.exists() else {}
        test_mission = fold_info.get("test_mission", fold_dir.name)

        folds.append({
            "train_X": train_X.values.astype(np.float32),
            "train_y": train_y.values.astype(int),
            "test_X": test_X.values.astype(np.float32),
            "test_y": test_y.values.astype(int),
            "test_mission": test_mission,
            "fold_name": fold_dir.name,
            "n_train": len(train_y),
            "n_test": len(test_y),
            "gene_names": list(common_genes),
        })

    return folds


# ── Metrics (DD-08 + DD-32) ──────────────────────────────────────────────────

def bootstrap_auroc(y_true, y_score, n_bootstrap=2000, alpha=0.05, seed=42):
    """Bootstrap 95% CI for AUROC. v4 default: n=2000 (v1 was 1000)."""
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


def permutation_pvalue(y_true, y_score, n_perm=1000, seed=42):
    """Permutation p-value for AUROC. H0: AUROC = random.
    DD-32: n_perm=1000 for all evaluations. Use n=5000 only for best-method report.
    """
    rng = np.random.default_rng(seed)
    observed = roc_auc_score(y_true, y_score)
    null_scores = []
    for _ in range(n_perm):
        y_perm = rng.permutation(y_true)
        null_scores.append(roc_auc_score(y_perm, y_score))
    p = float((np.sum(np.array(null_scores) >= observed) + 1) / (n_perm + 1))
    return p


# ── Results I/O (DD-29) ─────────────────────────────────────────────────────

def save_evaluation_result(tissue, feature_type, method_key, result_dict, output_dir=None):
    """Save evaluation result as JSON matching existing A1 schema."""
    if output_dir is None:
        output_dir = V4_EVAL_DIR
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fname = f"M1_{tissue}_{feature_type}_{method_key}.json"
    (output_dir / fname).write_text(json.dumps(result_dict, indent=2))
    return output_dir / fname


# ── Verification ──────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("v4_utils.py Verification")
    print("=" * 60)

    # Check tissue configurations
    for tissue, missions in TISSUE_MISSIONS.items():
        cv = "LOMO" if tissue in LOMO_TISSUES else "5-fold"
        print(f"  {tissue:15s}: {len(missions)} missions ({', '.join(missions)}) → {cv}")

    print(f"\n  LOMO tissues: {LOMO_TISSUES}")
    print(f"  K-fold tissues: {KFOLD_TISSUES}")
    print(f"  Task-mapped: {list(TASK_MAP.keys())}")

    # Test data loading for each tissue
    print(f"\n{'='*60}")
    print("Data loading test:")
    for tissue in TISSUE_MISSIONS:
        try:
            meta = load_metadata(tissue)
            genes = load_gene_features(tissue)
            genes, meta = align_features_with_meta(genes, meta)
            labels = meta["label"].value_counts()
            # Show binary encoding
            y, valid = encode_labels(meta)
            n_flight = int((y[valid] == 1).sum())
            n_ctrl = int((y[valid] == 0).sum())
            n_excluded = int((~valid).sum())
            print(f"  {tissue:15s}: {len(meta)} total, {n_flight} FLT + {n_ctrl} ctrl "
                  f"= {n_flight+n_ctrl} binary (excl {n_excluded}), "
                  f"{len(genes.columns)} genes")
        except Exception as e:
            print(f"  {tissue:15s}: FAILED — {e}")

    # Test fold generation
    print(f"\n{'='*60}")
    print("Fold generation test:")
    for tissue in TISSUE_MISSIONS:
        try:
            folds = get_folds(tissue)
            fold_info = [(f["test_mission"], f["n_train"], f["n_test"]) for f in folds]
            print(f"  {tissue:15s}: {len(folds)} folds — {fold_info}")
        except Exception as e:
            print(f"  {tissue:15s}: FAILED — {e}")
