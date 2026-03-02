#!/usr/bin/env python3
"""
cross_tissue_transfer.py — GeneLab_benchmark: Category C Cross-Tissue Transfer

Tests H3: "Pathway-level features transfer better across tissues than gene-level."

Three transfer methods:
  A) Gene-level:    Common genes → PCA → LR (train tissue → test tissue)
  B) DEG overlap:   Shared DEGs → PCA → LR
  C) Pathway-level: fGSEA top-N → GSVA scores → LR

Output:
  processed/C_cross_tissue/
    C{n}_{train}_{test}_results.json         # per-pair results
    C{n}_{train}_{test}_pathway_selection.csv # Method C pathway list
  evaluation/
    C_cross_tissue_summary.json              # all pairs summary + H3 test

Usage:
  python scripts/cross_tissue_transfer.py --pair C1
  python scripts/cross_tissue_transfer.py --all
  python scripts/cross_tissue_transfer.py --all --no-bootstrap
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
FGSEA_DIR = BASE_DIR / "processed" / "fgsea"
PATHWAY_DIR = BASE_DIR / "processed" / "pathway_scores"
C_OUTPUT_DIR = BASE_DIR / "processed" / "C_cross_tissue"
RESULTS_DIR = BASE_DIR / "evaluation"
DATA_DIR = BASE_DIR / "data" / "mouse"

# ── Config ─────────────────────────────────────────────────────────────────────
FLIGHT_LABEL = "Flight"
GROUND_LABELS = {"GC", "VC"}
VARIANCE_PERCENTILE = 0.25
N_BOOTSTRAP = 2000
CI_ALPHA = 0.05
N_PERMUTATIONS = 10000
TOP_N_PATHWAYS = 20
DEG_PADJ_THRESHOLD = 0.05
DEG_LFC_THRESHOLD = 0.5

# ── Tissue-Mission mapping ────────────────────────────────────────────────────
TISSUE_MISSIONS = {
    "liver": ["RR-1", "RR-3", "RR-6", "RR-8", "RR-9", "MHU-2"],
    "kidney": ["RR-1", "RR-3", "RR-7"],
    "thymus": ["RR-6", "MHU-2", "RR-9"],
    "gastrocnemius": ["RR-1", "RR-5", "RR-9"],
    "eye": ["RR-1", "RR-3", "TBD"],
}

# GLDS IDs for DGE file lookup
MISSION_GLDS = {
    "liver": {"RR-1": "GLDS-48", "RR-3": "GLDS-137", "RR-6": "GLDS-242",
              "RR-8": "GLDS-379", "RR-9": "GLDS-486", "MHU-2": "GLDS-245"},
    "kidney": {"RR-1": "GLDS-48", "RR-3": "GLDS-137", "RR-7": "GLDS-352"},
    "thymus": {"RR-6": "GLDS-242", "MHU-2": "GLDS-245", "RR-9": "GLDS-486"},
    "gastrocnemius": {"RR-1": "GLDS-48", "RR-9": "GLDS-486"},
    "eye": {"RR-1": "GLDS-48", "RR-3": "GLDS-137", "TBD": "GLDS-397"},
}

# Cross-tissue task definitions
CROSS_TISSUE_TASKS = {
    "C1": {"train": "liver", "test": "kidney",
            "biology": "metabolism, oxidative stress"},
    "C2": {"train": "liver", "test": "gastrocnemius",
            "biology": "energy metabolism"},
    "C3": {"train": "liver", "test": "thymus",
            "biology": "immune response (Kupffer cell)"},
    "C4": {"train": "thymus", "test": "kidney",
            "biology": "immune-renal axis"},
}


# ── Statistical Utilities ────────────────────────────────────────────────────

def bootstrap_auroc_ci(y_true, y_score, n_boot=None, alpha=CI_ALPHA,
                       random_state=42):
    """Bootstrap 95% CI for AUROC."""
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

    ci_low = float(np.percentile(boot_aurocs, 100 * alpha / 2))
    ci_high = float(np.percentile(boot_aurocs, 100 * (1 - alpha / 2)))
    return ci_low, ci_high


def permutation_test(y_true, y_score, n_perm=None, random_state=42):
    """Permutation test: H0 = AUROC == 0.5 (random). Returns p-value."""
    from sklearn.metrics import roc_auc_score

    if n_perm is None:
        n_perm = N_PERMUTATIONS
    if n_perm <= 0:
        return np.nan

    rng = np.random.RandomState(random_state)
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
            perm_auroc = roc_auc_score(perm_y, y_score)
            if perm_auroc >= observed:
                count += 1
        except Exception:
            continue

    return (count + 1) / (n_perm + 1)


# ── Data Loading ─────────────────────────────────────────────────────────────

def load_tissue_gene_data(tissue):
    """Load tissue-wide log2 normalized counts + metadata."""
    tissue_dir = PROCESSED_DIR / tissue
    counts_path = tissue_dir / f"{tissue}_all_missions_log2_norm.csv"
    meta_path = tissue_dir / f"{tissue}_all_missions_metadata.csv"

    if not counts_path.exists() or not meta_path.exists():
        print(f"  [ERROR] Data not found for {tissue}")
        return None, None

    counts = pd.read_csv(counts_path, index_col=0)
    meta = pd.read_csv(meta_path, index_col=0)

    # Drop non-gene columns
    non_gene_cols = [c for c in counts.columns
                     if c in {"mission", "osd_id", "label"}]
    if non_gene_cols:
        counts = counts.drop(columns=non_gene_cols)
    counts = counts.select_dtypes(include=[np.number])

    # Align
    common = counts.index.intersection(meta.index)
    counts = counts.loc[common]
    meta = meta.loc[common]

    return counts, meta


def load_gsva_scores(tissue, db="hallmark"):
    """Load GSVA scores across all missions for a tissue.
    Prefixes sample names with 'mission.' to match metadata conventions.
    """
    all_scores = []

    for mission in TISSUE_MISSIONS.get(tissue, []):
        f = PATHWAY_DIR / tissue / f"{mission}_gsva_{db}.csv"
        if not f.exists():
            continue
        scores = pd.read_csv(f, index_col=0)
        scores["mission"] = mission
        all_scores.append(scores)

    if not all_scores:
        return None
    return pd.concat(all_scores)


def load_fgsea_results(tissue, db="hallmark"):
    """Load and concatenate fGSEA results across all missions."""
    dfs = []
    for mission in TISSUE_MISSIONS.get(tissue, []):
        f = FGSEA_DIR / tissue / f"{mission}_fgsea_{db}.csv"
        if not f.exists():
            continue
        dfs.append(pd.read_csv(f))
    if not dfs:
        return None
    return pd.concat(dfs, ignore_index=True)


def get_binary_labels(meta):
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


def _align_samples(gsva_index, label_index):
    """
    Align GSVA sample names with metadata sample names.
    Handles mission-prefix mismatch: metadata may use 'RR-1.SampleName'
    while GSVA uses 'SampleName'.

    Returns (common_pairs, gsva_idx_list, label_idx_list) where each
    list preserves alignment.
    """
    gsva_set = set(gsva_index)
    label_set = set(label_index)

    # Try direct match first
    direct = gsva_set & label_set
    if len(direct) >= 5:
        direct = sorted(direct)
        return direct, direct, direct

    # Try stripping mission prefix from label index ("RR-1.Sample" → "Sample")
    label_to_stripped = {}
    for lab in label_index:
        parts = str(lab).split(".", 1)
        if len(parts) == 2:
            label_to_stripped[lab] = parts[1]
        else:
            label_to_stripped[lab] = lab

    gsva_idx_list = []
    label_idx_list = []
    for lab_orig, lab_stripped in label_to_stripped.items():
        if lab_stripped in gsva_set:
            gsva_idx_list.append(lab_stripped)
            label_idx_list.append(lab_orig)

    if len(gsva_idx_list) >= 5:
        return list(range(len(gsva_idx_list))), gsva_idx_list, label_idx_list

    # Try stripping prefix from GSVA index
    for gsva_name in gsva_index:
        parts = str(gsva_name).split(".", 1)
        if len(parts) == 2 and parts[1] in label_set:
            gsva_idx_list.append(gsva_name)
            label_idx_list.append(parts[1])

    return list(range(len(gsva_idx_list))), gsva_idx_list, label_idx_list


def load_tissue_degs(tissue, db="hallmark"):
    """
    Load DEG information from fGSEA leading edge genes.
    Returns set of significant gene symbols for the tissue.
    """
    fgsea = load_fgsea_results(tissue, db)
    if fgsea is None:
        return set()

    # Collect leading edge genes from significant pathways
    sig = fgsea[fgsea["padj"] < DEG_PADJ_THRESHOLD]
    all_genes = set()
    for _, row in sig.iterrows():
        if pd.notna(row.get("leadingEdge_str")):
            genes = str(row["leadingEdge_str"]).split("; ")
            all_genes.update(genes)
    return all_genes


def select_top_pathways(tissue, db="hallmark", top_n=None):
    """
    Select top-N pathways from a tissue's fGSEA results by |mean NES|.
    Returns list of pathway names and a DataFrame of pathway stats.
    """
    if top_n is None:
        top_n = TOP_N_PATHWAYS

    fgsea = load_fgsea_results(tissue, db)
    if fgsea is None:
        return [], None

    pathway_agg = (fgsea
                   .groupby("pathway")
                   .agg(
                       mean_NES=("NES", "mean"),
                       abs_mean_NES=("NES", lambda x: abs(x.mean())),
                       n_missions=("mission", "nunique"),
                       n_significant=("padj", lambda x: (x < 0.05).sum()),
                   )
                   .sort_values("abs_mean_NES", ascending=False))

    selected = pathway_agg.head(top_n).index.tolist()
    return selected, pathway_agg.loc[selected]


# ── Transfer Methods ──────────────────────────────────────────────────────────

def method_a_gene_transfer(train_counts, train_labels,
                           test_counts, test_labels):
    """
    Method A: Gene-level direct transfer.
    Find common genes → variance filter (train) → PCA → LR → test.
    """
    from sklearn.metrics import roc_auc_score
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import Pipeline

    _nan = {"auroc": np.nan, "y_true": None, "y_score": None}

    # Find common genes
    common_genes = train_counts.columns.intersection(test_counts.columns)
    if len(common_genes) < 100:
        print(f"    [WARN] Only {len(common_genes)} common genes")
        if len(common_genes) < 10:
            return _nan

    X_train = train_counts[common_genes].copy()
    X_test = test_counts[common_genes].copy()

    # Variance filter (train only, DD-03)
    gene_var = X_train.var(axis=0)
    threshold = gene_var.quantile(VARIANCE_PERCENTILE)
    keep_genes = gene_var[gene_var >= threshold].index.tolist()
    X_train = X_train[keep_genes]
    X_test = X_test[keep_genes]

    y_train = train_labels.values.astype(int)
    y_test = test_labels.values.astype(int)

    if len(np.unique(y_train)) < 2 or len(np.unique(y_test)) < 2:
        return _nan

    # StandardScaler + PCA + LR
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
                "n_features": len(keep_genes), "n_common_genes": len(common_genes)}
    except Exception as e:
        print(f"    [ERROR] Method A: {e}")
        return _nan


def method_b_deg_transfer(train_counts, train_labels,
                          test_counts, test_labels,
                          train_tissue, test_tissue, db="hallmark"):
    """
    Method B: DEG overlap transfer.
    Find DEGs in each tissue → intersect → PCA → LR.
    Uses fGSEA leading edge genes as DEG proxy (consistent with pipeline).
    """
    from sklearn.metrics import roc_auc_score
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import Pipeline

    _nan = {"auroc": np.nan, "y_true": None, "y_score": None}

    # Get DEGs from both tissues
    train_degs = load_tissue_degs(train_tissue, db)
    test_degs = load_tissue_degs(test_tissue, db)

    if not train_degs or not test_degs:
        print(f"    [WARN] No DEGs found: train={len(train_degs)}, test={len(test_degs)}")
        return _nan

    shared_degs = train_degs & test_degs
    print(f"    DEGs: train={len(train_degs)}, test={len(test_degs)}, "
          f"shared={len(shared_degs)}")

    if len(shared_degs) < 5:
        print(f"    [WARN] Too few shared DEGs")
        return _nan

    # Map gene symbols to ENSEMBL IDs in count data
    # Count data uses ENSMUSG IDs; DEGs are symbols
    # Need symbol-to-ensembl mapping
    # For now, check if count columns contain symbols or ENSMUSGs
    sample_col = train_counts.columns[0]
    if sample_col.startswith("ENSMUSG"):
        # Need to map symbols to ENSMUSG - use DGE file for mapping
        ensembl_map = _build_ensembl_map(train_tissue)
        if ensembl_map is None:
            print("    [WARN] Cannot map symbols to ENSEMBL")
            return _nan

        shared_ensembl = set()
        for sym in shared_degs:
            if sym in ensembl_map:
                ens = ensembl_map[sym]
                if ens in train_counts.columns and ens in test_counts.columns:
                    shared_ensembl.add(ens)

        if len(shared_ensembl) < 5:
            print(f"    [WARN] Only {len(shared_ensembl)} mapped shared DEGs")
            return _nan

        feature_genes = sorted(shared_ensembl)
    else:
        # Columns are already symbols
        feature_genes = sorted(shared_degs & set(train_counts.columns)
                               & set(test_counts.columns))

    print(f"    Using {len(feature_genes)} shared DEG features")

    X_train = train_counts[feature_genes].values.astype(np.float32)
    X_test = test_counts[feature_genes].values.astype(np.float32)
    y_train = train_labels.values.astype(int)
    y_test = test_labels.values.astype(int)

    if len(np.unique(y_train)) < 2 or len(np.unique(y_test)) < 2:
        return _nan

    scaler = StandardScaler()
    X_tr = scaler.fit_transform(X_train)
    X_te = scaler.transform(X_test)

    n_comps = min(min(20, len(feature_genes) // 2),
                  X_tr.shape[0] - 1, X_tr.shape[1])
    if n_comps < 2:
        # Few features: use LR directly without PCA
        model = LogisticRegression(C=1.0, class_weight="balanced",
                                   max_iter=1000, random_state=42)
    else:
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
                "n_features": len(feature_genes),
                "n_shared_degs": len(shared_degs)}
    except Exception as e:
        print(f"    [ERROR] Method B: {e}")
        return _nan


def _build_ensembl_map(tissue):
    """Build SYMBOL → ENSEMBL mapping from first available DGE file.
    Falls back to other tissues if the specified tissue has no DGE files.
    """
    # Try specified tissue first, then all others
    tissues_to_try = [tissue] + [t for t in TISSUE_MISSIONS if t != tissue]
    for t in tissues_to_try:
        for mission in TISSUE_MISSIONS.get(t, []):
            glds = MISSION_GLDS.get(t, {}).get(mission)
            if not glds:
                continue
            dge_path = DATA_DIR / t / mission / f"{glds}_rna_seq_differential_expression_GLbulkRNAseq.csv"
            if not dge_path.exists():
                dge_path = DATA_DIR / t / mission / f"{glds}_rna_seq_differential_expression.csv"
            if dge_path.exists():
                try:
                    dge = pd.read_csv(dge_path, usecols=["ENSEMBL", "SYMBOL"],
                                      nrows=50000)
                    dge = dge.dropna(subset=["ENSEMBL", "SYMBOL"])
                    if len(dge) > 100:
                        return dict(zip(dge["SYMBOL"], dge["ENSEMBL"]))
                except Exception:
                    continue
    return None


def method_c_pathway_transfer(train_tissue, test_tissue,
                               train_labels, test_labels,
                               train_meta, test_meta,
                               db="hallmark", top_n=None):
    """
    Method C: Pathway-level transfer.
    Select top-N pathways from train tissue fGSEA → GSVA scores → LR.
    """
    from sklearn.metrics import roc_auc_score
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression

    if top_n is None:
        top_n = TOP_N_PATHWAYS

    _nan = {"auroc": np.nan, "y_true": None, "y_score": None}

    # Step 1: Select pathways from train tissue
    selected, pathway_stats = select_top_pathways(train_tissue, db, top_n)
    if not selected:
        print("    [ERROR] No pathways selected")
        return _nan

    print(f"    Selected {len(selected)} pathways from {train_tissue}")
    for i, pw in enumerate(selected[:5]):
        row = pathway_stats.loc[pw]
        print(f"      {i+1}. {pw[:50]:50s} NES={row['mean_NES']:+.2f} "
              f"(sig {row['n_significant']:.0f}/{row['n_missions']:.0f})")

    # Step 2: Load GSVA scores
    train_gsva = load_gsva_scores(train_tissue, db)
    test_gsva = load_gsva_scores(test_tissue, db)

    if train_gsva is None or test_gsva is None:
        print("    [ERROR] Missing GSVA scores")
        return _nan

    # Drop mission column for features
    train_gsva_feat = train_gsva.drop(columns=["mission"], errors="ignore")
    test_gsva_feat = test_gsva.drop(columns=["mission"], errors="ignore")

    # Filter to selected pathways (available in both)
    common_pathways = [p for p in selected
                       if p in train_gsva_feat.columns
                       and p in test_gsva_feat.columns]
    if len(common_pathways) < 3:
        print(f"    [ERROR] Only {len(common_pathways)} common pathways")
        return _nan

    # Align samples with labels — handle mission-prefix mismatch
    # Metadata may have "RR-1.SampleName" while GSVA has "SampleName"
    train_common, train_gsva_idx, train_label_idx = _align_samples(
        train_gsva_feat.index, train_labels.index)
    test_common, test_gsva_idx, test_label_idx = _align_samples(
        test_gsva_feat.index, test_labels.index)

    if len(train_common) < 5 or len(test_common) < 5:
        print(f"    [ERROR] Too few aligned samples: "
              f"train={len(train_common)}, test={len(test_common)}")
        return _nan

    X_train = train_gsva_feat.loc[train_gsva_idx, common_pathways].values.astype(np.float32)
    X_test = test_gsva_feat.loc[test_gsva_idx, common_pathways].values.astype(np.float32)
    y_train = train_labels.loc[train_label_idx].values.astype(int)
    y_test = test_labels.loc[test_label_idx].values.astype(int)

    if len(np.unique(y_train)) < 2 or len(np.unique(y_test)) < 2:
        return _nan

    # LR (no PCA needed — only ~20 features)
    scaler = StandardScaler()
    X_tr = scaler.fit_transform(X_train)
    X_te = scaler.transform(X_test)

    model = LogisticRegression(C=1.0, class_weight="balanced",
                               max_iter=1000, random_state=42)

    try:
        model.fit(X_tr, y_train)
        y_score = model.predict_proba(X_te)[:, 1]
        auroc = float(roc_auc_score(y_test, y_score))

        # Feature importance
        coef = model.coef_[0]
        importance = sorted(zip(common_pathways, coef),
                           key=lambda x: abs(x[1]), reverse=True)

        return {"auroc": auroc, "y_true": y_test, "y_score": y_score,
                "n_features": len(common_pathways),
                "selected_pathways": common_pathways,
                "pathway_stats": pathway_stats,
                "top_features": importance[:10]}
    except Exception as e:
        print(f"    [ERROR] Method C: {e}")
        return _nan


# ── Main Transfer Evaluation ─────────────────────────────────────────────────

def evaluate_pair(task_id, train_tissue, test_tissue, db="hallmark",
                  verbose=True):
    """
    Run all three methods for a tissue pair and compare.
    """
    print(f"\n{'='*70}")
    print(f"Category C: {task_id} — {train_tissue} → {test_tissue}")
    print(f"{'='*70}")

    # Load gene-level data
    train_counts, train_meta = load_tissue_gene_data(train_tissue)
    test_counts, test_meta = load_tissue_gene_data(test_tissue)

    if train_counts is None or test_counts is None:
        print("[ERROR] Cannot load data")
        return None

    # Get labels
    train_labels = get_binary_labels(train_meta)
    test_labels = get_binary_labels(test_meta)

    # Remove NaN labels
    train_valid = ~train_labels.isna()
    test_valid = ~test_labels.isna()

    train_counts_v = train_counts[train_valid]
    train_labels_v = train_labels[train_valid]
    train_meta_v = train_meta[train_valid]
    test_counts_v = test_counts[test_valid]
    test_labels_v = test_labels[test_valid]
    test_meta_v = test_meta[test_valid]

    n_train = len(train_labels_v)
    n_test = len(test_labels_v)
    print(f"\n  Train ({train_tissue}): {n_train} samples "
          f"({int(train_labels_v.sum())} Flight, "
          f"{int((train_labels_v == 0).sum())} Ground)")
    print(f"  Test  ({test_tissue}):  {n_test} samples "
          f"({int(test_labels_v.sum())} Flight, "
          f"{int((test_labels_v == 0).sum())} Ground)")

    results = {"task_id": task_id, "train_tissue": train_tissue,
               "test_tissue": test_tissue, "db": db,
               "n_train": n_train, "n_test": n_test,
               "timestamp": datetime.now().isoformat()}

    # ── Method A: Gene-level ──
    print(f"\n  --- Method A: Gene-level transfer ---")
    res_a = method_a_gene_transfer(
        train_counts_v, train_labels_v, test_counts_v, test_labels_v)
    results["method_a"] = _format_result(res_a, "A")

    # ── Method B: DEG overlap ──
    print(f"\n  --- Method B: DEG overlap transfer ---")
    res_b = method_b_deg_transfer(
        train_counts_v, train_labels_v, test_counts_v, test_labels_v,
        train_tissue, test_tissue, db)
    results["method_b"] = _format_result(res_b, "B")

    # ── Method C: Pathway-level ──
    print(f"\n  --- Method C: Pathway-level transfer ---")
    res_c = method_c_pathway_transfer(
        train_tissue, test_tissue,
        train_labels_v, test_labels_v,
        train_meta_v, test_meta_v,
        db)
    results["method_c"] = _format_result(res_c, "C")

    # ── Per-mission breakdown (test tissue) ──
    print(f"\n  --- Per-test-mission breakdown ---")
    test_missions = sorted(test_meta_v["mission"].unique())
    mission_results = {}
    for tm in test_missions:
        tm_mask = test_meta_v["mission"] == tm
        n_tm = tm_mask.sum()
        if n_tm < 3:
            continue

        # Method C per-mission
        res_c_m = method_c_pathway_transfer(
            train_tissue, test_tissue,
            train_labels_v,
            test_labels_v[tm_mask],
            train_meta_v,
            test_meta_v[tm_mask],
            db)
        auroc_m = res_c_m["auroc"] if not np.isnan(res_c_m.get("auroc", np.nan)) else None

        mission_results[tm] = {
            "n_samples": int(n_tm),
            "method_c_auroc": auroc_m,
        }
        status = f"{auroc_m:.3f}" if auroc_m is not None else "N/A"
        print(f"    {tm} ({n_tm} samples): Method C AUROC = {status}")

    results["per_mission"] = mission_results

    # ── Save ──
    C_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = C_OUTPUT_DIR / f"{task_id}_{train_tissue}_{test_tissue}_results.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=_json_default)
    print(f"\n  Saved: {out_path.name}")

    # Save pathway selection if Method C succeeded
    if (res_c.get("pathway_stats") is not None and
            res_c.get("selected_pathways")):
        pw_path = C_OUTPUT_DIR / f"{task_id}_pathway_selection_{db}.csv"
        res_c["pathway_stats"].to_csv(pw_path)
        print(f"  Saved: {pw_path.name}")

    # ── Summary ──
    print(f"\n  {'Method':<12} {'AUROC':>8}  {'95% CI':>16}  {'p-value':>8}  {'Features':>8}")
    print(f"  {'-'*60}")
    for method_key, label in [("method_a", "A: Gene"),
                               ("method_b", "B: DEG"),
                               ("method_c", "C: Pathway")]:
        r = results[method_key]
        auroc_str = f"{r['auroc']:.3f}" if r["auroc"] is not None else "N/A"
        ci_str = (f"[{r['ci_low']:.3f}, {r['ci_high']:.3f}]"
                  if r.get("ci_low") is not None else "N/A")
        p_str = (f"{r['perm_p']:.4f}" if r.get("perm_p") is not None else "N/A")
        feat_str = str(r.get("n_features", "N/A"))
        print(f"  {label:<12} {auroc_str:>8}  {ci_str:>16}  {p_str:>8}  {feat_str:>8}")

    return results


def _format_result(res, method_label):
    """Format transfer result with bootstrap CI and permutation test."""
    if res.get("y_true") is None or np.isnan(res.get("auroc", np.nan)):
        return {"auroc": None, "ci_low": None, "ci_high": None,
                "perm_p": None, "n_features": res.get("n_features")}

    ci_low, ci_high = bootstrap_auroc_ci(res["y_true"], res["y_score"])
    perm_p = permutation_test(res["y_true"], res["y_score"])

    out = {
        "auroc": res["auroc"],
        "ci_low": ci_low if not np.isnan(ci_low) else None,
        "ci_high": ci_high if not np.isnan(ci_high) else None,
        "perm_p": perm_p if not np.isnan(perm_p) else None,
        "n_features": res.get("n_features"),
    }

    # Add extra info
    for key in ["n_common_genes", "n_shared_degs", "top_features"]:
        if key in res:
            val = res[key]
            if key == "top_features":
                val = [(name, float(coef)) for name, coef in val]
            out[key] = val

    return out


def _json_default(obj):
    """JSON serialization helper."""
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, pd.DataFrame):
        return obj.to_dict()
    if isinstance(obj, set):
        return sorted(obj)
    return str(obj)


# ── H3 Hypothesis Test ──────────────────────────────────────────────────────

def test_h3(all_results):
    """
    H3: Method C (pathway) > Method A (gene) for cross-tissue transfer.
    Paired comparison across tissue pairs.
    """
    print(f"\n{'='*70}")
    print("H3 Hypothesis Test: Pathway > Gene for Cross-Tissue Transfer")
    print(f"{'='*70}")

    method_a_aurocs = []
    method_b_aurocs = []
    method_c_aurocs = []
    pair_labels = []

    for task_id, res in all_results.items():
        a = res.get("method_a", {}).get("auroc")
        b = res.get("method_b", {}).get("auroc")
        c = res.get("method_c", {}).get("auroc")

        if a is not None:
            method_a_aurocs.append(a)
        if b is not None:
            method_b_aurocs.append(b)
        if c is not None:
            method_c_aurocs.append(c)
        pair_labels.append(task_id)

    h3_results = {
        "n_pairs": len(pair_labels),
        "pairs": pair_labels,
    }

    for label, vals in [("method_a", method_a_aurocs),
                         ("method_b", method_b_aurocs),
                         ("method_c", method_c_aurocs)]:
        if vals:
            mean_v = float(np.mean(vals))
            h3_results[f"{label}_mean"] = mean_v
            h3_results[f"{label}_values"] = [float(v) for v in vals]
            print(f"  {label}: mean={mean_v:.3f} ({vals})")

    # Paired sign test: C > A for each pair
    if method_a_aurocs and method_c_aurocs:
        n_pairs = min(len(method_a_aurocs), len(method_c_aurocs))
        c_wins = sum(1 for i in range(n_pairs)
                     if method_c_aurocs[i] > method_a_aurocs[i])
        h3_results["c_vs_a_wins"] = c_wins
        h3_results["c_vs_a_total"] = n_pairs
        print(f"\n  Method C > Method A: {c_wins}/{n_pairs} pairs")

        # Mean difference
        diffs = [method_c_aurocs[i] - method_a_aurocs[i]
                 for i in range(n_pairs)]
        h3_results["c_minus_a_mean"] = float(np.mean(diffs))
        h3_results["c_minus_a_diffs"] = [float(d) for d in diffs]
        print(f"  Mean AUROC difference (C-A): {np.mean(diffs):+.3f}")

    # B vs C
    if method_b_aurocs and method_c_aurocs:
        n_pairs = min(len(method_b_aurocs), len(method_c_aurocs))
        c_wins_b = sum(1 for i in range(n_pairs)
                       if method_c_aurocs[i] > method_b_aurocs[i])
        h3_results["c_vs_b_wins"] = c_wins_b
        h3_results["c_vs_b_total"] = n_pairs
        print(f"  Method C > Method B: {c_wins_b}/{n_pairs} pairs")

    return h3_results


# ── Sensitivity Analyses ────────────────────────────────────────────────────

def shared_mission_sensitivity(task_id, train_tissue, test_tissue, db="hallmark"):
    """
    C1 sensitivity: exclude shared missions between tissues to check
    batch confounding. If AUROC drops significantly, the original result
    was inflated by shared-mission batch effects.
    """
    from sklearn.metrics import roc_auc_score
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import Pipeline

    train_missions = set(TISSUE_MISSIONS.get(train_tissue, []))
    test_missions = set(TISSUE_MISSIONS.get(test_tissue, []))
    shared = train_missions & test_missions

    if not shared:
        print(f"  No shared missions between {train_tissue} and {test_tissue}")
        return None

    print(f"\n  === Shared Mission Sensitivity: {task_id} ===")
    print(f"  Shared missions: {sorted(shared)}")

    # Load data
    train_counts, train_meta = load_tissue_gene_data(train_tissue)
    test_counts, test_meta = load_tissue_gene_data(test_tissue)
    if train_counts is None or test_counts is None:
        return None

    train_labels = get_binary_labels(train_meta)
    test_labels = get_binary_labels(test_meta)

    train_valid = ~train_labels.isna()
    test_valid = ~test_labels.isna()
    train_counts = train_counts[train_valid]
    train_meta = train_meta[train_valid]
    train_labels = train_labels[train_valid]
    test_counts = test_counts[test_valid]
    test_meta = test_meta[test_valid]
    test_labels = test_labels[test_valid]

    results = {}

    # Test on non-shared missions only
    non_shared_mask = ~test_meta["mission"].isin(shared)
    if non_shared_mask.sum() >= 3:
        # Method A on non-shared test missions
        common_genes = train_counts.columns.intersection(test_counts.columns)
        gene_var = train_counts[common_genes].var(axis=0)
        threshold = gene_var.quantile(VARIANCE_PERCENTILE)
        keep = gene_var[gene_var >= threshold].index.tolist()

        X_tr = StandardScaler().fit_transform(
            train_counts[keep].values.astype(np.float32))
        scaler = StandardScaler().fit(train_counts[keep].values.astype(np.float32))
        X_tr = scaler.transform(train_counts[keep].values.astype(np.float32))
        X_te = scaler.transform(
            test_counts.loc[non_shared_mask, keep].values.astype(np.float32))
        y_tr = train_labels.values.astype(int)
        y_te = test_labels[non_shared_mask].values.astype(int)

        n_comps = min(50, X_tr.shape[0] - 1, X_tr.shape[1])
        model = Pipeline([
            ("pca", PCA(n_components=n_comps, random_state=42)),
            ("clf", LogisticRegression(C=1.0, class_weight="balanced",
                                       max_iter=1000, random_state=42))
        ])
        model.fit(X_tr, y_tr)
        y_score = model.predict_proba(X_te)[:, 1]
        auroc_excl = float(roc_auc_score(y_te, y_score))

        # Also test on shared missions only
        shared_mask = test_meta["mission"].isin(shared)
        X_te_sh = scaler.transform(
            test_counts.loc[shared_mask, keep].values.astype(np.float32))
        y_te_sh = test_labels[shared_mask].values.astype(int)
        y_score_sh = model.predict_proba(X_te_sh)[:, 1]
        auroc_shared = float(roc_auc_score(y_te_sh, y_score_sh))

        results = {
            "shared_missions": sorted(shared),
            "method_a_all": None,  # filled by caller
            "method_a_shared_only": auroc_shared,
            "method_a_excl_shared": auroc_excl,
            "n_shared_test": int(shared_mask.sum()),
            "n_non_shared_test": int(non_shared_mask.sum()),
        }

        print(f"  Method A on shared missions ({sorted(shared)}): {auroc_shared:.3f}")
        print(f"  Method A excluding shared: {auroc_excl:.3f}")
        print(f"  Batch confounding: {'LIKELY' if auroc_shared > auroc_excl + 0.1 else 'UNLIKELY'}")

    return results


# ── CLI ────────────────────────────────────────────────────────────────────────

def main():
    global N_BOOTSTRAP, N_PERMUTATIONS, TOP_N_PATHWAYS

    parser = argparse.ArgumentParser(
        description="Category C: Cross-Tissue Transfer (GeneLab_benchmark)")
    parser.add_argument("--pair", choices=list(CROSS_TISSUE_TASKS.keys()),
                        help="Run a specific tissue pair (e.g., C1)")
    parser.add_argument("--all", action="store_true",
                        help="Run all tissue pairs")
    parser.add_argument("--db", default="hallmark",
                        help="Gene set DB (default: hallmark)")
    parser.add_argument("--top-n", type=int, default=20,
                        help="Top pathways for Method C (default: 20)")
    parser.add_argument("--no-bootstrap", action="store_true",
                        help="Skip bootstrap CI and permutation tests")
    args = parser.parse_args()

    if args.no_bootstrap:
        N_BOOTSTRAP = 0
        N_PERMUTATIONS = 0

    TOP_N_PATHWAYS = args.top_n

    all_results = {}

    if args.pair:
        task = CROSS_TISSUE_TASKS[args.pair]
        result = evaluate_pair(args.pair, task["train"], task["test"], args.db)
        if result:
            all_results[args.pair] = result

    elif args.all:
        for task_id, task in CROSS_TISSUE_TASKS.items():
            result = evaluate_pair(task_id, task["train"], task["test"],
                                    args.db)
            if result:
                all_results[task_id] = result

        # Shared mission sensitivity (P1)
        sensitivity = {}
        for task_id, task in CROSS_TISSUE_TASKS.items():
            sens = shared_mission_sensitivity(
                task_id, task["train"], task["test"], args.db)
            if sens:
                a_auroc = all_results.get(task_id, {}).get("method_a", {}).get("auroc")
                sens["method_a_all"] = a_auroc
                sensitivity[task_id] = sens

        # H3 test
        if len(all_results) >= 2:
            h3 = test_h3(all_results)

            # Save summary
            RESULTS_DIR.mkdir(parents=True, exist_ok=True)
            summary = {
                "tasks": {},
                "h3_test": h3,
                "sensitivity": sensitivity,
                "config": {
                    "db": args.db,
                    "top_n_pathways": args.top_n,
                    "n_bootstrap": N_BOOTSTRAP,
                    "n_permutations": N_PERMUTATIONS,
                    "timestamp": datetime.now().isoformat(),
                },
            }
            for tid, res in all_results.items():
                summary["tasks"][tid] = {
                    "train": res["train_tissue"],
                    "test": res["test_tissue"],
                    "n_train": res["n_train"],
                    "n_test": res["n_test"],
                    "method_a": res["method_a"],
                    "method_b": res["method_b"],
                    "method_c": res["method_c"],
                }

            summary_path = RESULTS_DIR / "C_cross_tissue_summary.json"
            with open(summary_path, "w") as f:
                json.dump(summary, f, indent=2, default=_json_default)
            print(f"\n  Summary saved: {summary_path}")

    else:
        parser.error("Specify --pair or --all")


if __name__ == "__main__":
    main()
