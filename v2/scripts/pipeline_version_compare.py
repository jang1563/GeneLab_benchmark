#!/usr/bin/env python3
"""
pipeline_version_compare.py — GeneLab_benchmark v2.0: J1 Pipeline Version Comparison

Compares GeneLab pipeline versions by evaluating the same RR-1 samples
processed with two different pipelines:
  - GLDS-48:  Original pipeline (v1), 14 samples (FLT_I/C + GC_I/C)
  - GLDS-168: Reprocessed pipeline (v2, merged RR-1+RR-3), 20+ samples

J1 Question: Does pipeline version significantly affect ML performance?

Matching: GLDS-48 CAGES (C) subset ↔ GLDS-168 noERCC subset
  FLT: M25, M26, M28, M30 (GLDS-48 has M27; GLDS-168 has M29 — keep intersection)
  GC:  M36, M37, M38, M39, M40

Analyses:
  1. QC: Per-sample expression correlation between pipelines
  2. LOMO: Rebuild liver LOMO folds with GLDS-168 replacing GLDS-48
  3. fGSEA: NES correlation between pipeline versions
  4. Gene selection: Jaccard overlap of selected genes

Output:
  v2/processed/J1_pipeline/
  v2/evaluation/J1_pipeline_comparison.json

Usage:
  python v2/scripts/pipeline_version_compare.py
"""

import json
import sys
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from scipy import stats

warnings.filterwarnings("ignore", category=FutureWarning)

V2_DIR = Path(__file__).resolve().parent.parent
PROJECT_DIR = V2_DIR.parent
sys.path.insert(0, str(PROJECT_DIR / "scripts"))

from utils import (
    load_metadata, load_gene_features, load_pathway_features,
    align_features_with_meta,
)

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR = PROJECT_DIR / "data" / "mouse" / "liver"
PROCESSED_DIR = PROJECT_DIR / "processed" / "A_detection" / "liver"
J1_OUTPUT_DIR = V2_DIR / "processed" / "J1_pipeline"
RESULTS_DIR = V2_DIR / "evaluation"

# ── Sample mapping ────────────────────────────────────────────────────────────
# GLDS-48 Carcass(C) animal IDs that match GLDS-168 noERCC
# FLT_C: M25, M26, M27, M28, M30 (GLDS-48) ↔ M25, M26, M28, M29, M30 (GLDS-168)
#   Intersection: M25, M26, M28, M30 (M27 vs M29 mismatch — exclude both)
# GC_C: M36, M37, M38, M39, M40 (identical in both)

MATCHED_ANIMALS = {
    "FLT": ["M25", "M26", "M28", "M30"],  # M27/M29 excluded
    "GC": ["M36", "M37", "M38", "M39", "M40"],
}


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


def load_glds168_counts():
    """Load GLDS-168 normalized counts, filter to RR-1 noERCC FLT+GC."""
    path = DATA_DIR / "RR-1_3_(merged)" / "GLDS-168_rna_seq_Normalized_Counts.csv"
    print(f"  Loading GLDS-168: {path.name}")
    df = pd.read_csv(path, index_col=0)
    print(f"  Raw shape: {df.shape}")

    # Filter to RR-1 noERCC samples only
    rr1_noercc = [c for c in df.columns if "RR1" in c and "noERCC" in c]
    df_rr1 = df[rr1_noercc]
    print(f"  RR-1 noERCC samples: {len(rr1_noercc)}")

    # Filter to FLT + GC only (exclude BSL, VIV)
    flt_gc = [c for c in df_rr1.columns if ("_FLT_" in c or "_GC_" in c)]
    df_fltgc = df_rr1[flt_gc]
    print(f"  FLT+GC samples: {len(flt_gc)}")

    return df_fltgc


def load_glds48_counts():
    """Load GLDS-48 normalized counts."""
    path = DATA_DIR / "RR-1" / "GLDS-48_rna_seq_Normalized_Counts_GLbulkRNAseq.csv"
    print(f"  Loading GLDS-48: {path.name}")
    df = pd.read_csv(path, index_col=0)
    print(f"  Raw shape: {df.shape}")

    # Filter to Carcass (C) samples only (matching GLDS-168 noERCC)
    carcass = [c for c in df.columns if "_C_" in c]
    df_c = df[carcass]
    print(f"  Carcass(C) samples: {len(carcass)}")

    return df_c


def extract_animal_id(sample_name):
    """Extract animal ID (e.g., M25) from sample name."""
    parts = str(sample_name).split("_")
    for p in parts:
        if p.startswith("M") and p[1:].isdigit():
            return p
    return None


def match_samples(glds48_df, glds168_df):
    """Match samples between pipelines by animal ID."""
    matches = []

    for animal_id in MATCHED_ANIMALS["FLT"] + MATCHED_ANIMALS["GC"]:
        # Find in GLDS-48
        g48_col = None
        for c in glds48_df.columns:
            if animal_id == extract_animal_id(c):
                g48_col = c
                break

        # Find in GLDS-168
        g168_col = None
        for c in glds168_df.columns:
            if animal_id == extract_animal_id(c):
                g168_col = c
                break

        if g48_col and g168_col:
            condition = "FLT" if animal_id in MATCHED_ANIMALS["FLT"] else "GC"
            matches.append({
                "animal_id": animal_id,
                "condition": condition,
                "glds48_sample": g48_col,
                "glds168_sample": g168_col,
            })

    return matches


def compare_expression(glds48_df, glds168_df, matches):
    """Compare per-sample expression between pipelines."""
    print(f"\n  [Expression correlation]")
    print(f"  Matched samples: {len(matches)}")

    # Find common genes
    common_genes = sorted(set(glds48_df.index) & set(glds168_df.index))
    print(f"  Common genes: {len(common_genes)}")

    correlations = []
    for m in matches:
        g48_vals = glds48_df.loc[common_genes, m["glds48_sample"]].values.astype(float)
        g168_vals = glds168_df.loc[common_genes, m["glds168_sample"]].values.astype(float)

        # Remove NaN
        valid = ~(np.isnan(g48_vals) | np.isnan(g168_vals))
        g48_v = g48_vals[valid]
        g168_v = g168_vals[valid]

        # Log2 transform for correlation
        g48_log = np.log2(g48_v + 1)
        g168_log = np.log2(g168_v + 1)

        r, p = stats.pearsonr(g48_log, g168_log)
        correlations.append({
            "animal_id": m["animal_id"],
            "condition": m["condition"],
            "pearson_r": r,
            "p_value": p,
            "n_genes_valid": int(valid.sum()),
        })
        print(f"    {m['animal_id']} ({m['condition']}): r={r:.4f}")

    mean_r = np.mean([c["pearson_r"] for c in correlations])
    print(f"  Mean correlation: {mean_r:.4f}")

    return correlations, common_genes


def compare_lomo_performance(glds48_df, glds168_df, matches, common_genes):
    """Compare LOMO classification performance."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    from sklearn.metrics import roc_auc_score

    print(f"\n  [LOMO classification comparison]")

    # Load full liver data (all missions)
    full_gene = load_gene_features("liver")
    full_meta = load_metadata("liver")

    # Align
    full_gene_a, full_meta_a = align_features_with_meta(full_gene, full_meta)

    # For GLDS-48 (original): just use existing data
    rr1_mask_48 = full_meta_a["mission"] == "RR-1"
    other_mask = full_meta_a["mission"] != "RR-1"

    # Train on other missions, test on RR-1
    X_train = full_gene_a[other_mask].values.astype(float)
    X_test_48 = full_gene_a[rr1_mask_48].values.astype(float)
    y_train = (full_meta_a[other_mask]["label"] == "Flight").astype(int).values
    y_test_48 = (full_meta_a[rr1_mask_48]["label"] == "Flight").astype(int).values

    X_train = np.nan_to_num(X_train, nan=0.0)
    X_test_48 = np.nan_to_num(X_test_48, nan=0.0)

    # Pipeline
    scaler = StandardScaler()
    X_train_s = scaler.fit_transform(X_train)
    X_test_48_s = scaler.transform(X_test_48)

    n_comp = min(50, X_train_s.shape[0] - 1, X_train_s.shape[1])
    pca = PCA(n_components=n_comp, random_state=42)
    X_train_p = pca.fit_transform(X_train_s)
    X_test_48_p = pca.transform(X_test_48_s)

    clf = LogisticRegression(solver="lbfgs", max_iter=2000, C=1.0, random_state=42)
    clf.fit(X_train_p, y_train)

    if len(np.unique(y_test_48)) >= 2:
        proba_48 = clf.predict_proba(X_test_48_p)[:, list(clf.classes_).index(1)]
        auroc_48 = roc_auc_score(y_test_48, proba_48)
        print(f"  GLDS-48 (original) RR-1 fold AUROC: {auroc_48:.3f}")
    else:
        auroc_48 = None
        print(f"  GLDS-48: only one class in test set, cannot compute AUROC")

    # For GLDS-168: replace RR-1 data with GLDS-168 matched samples
    # Use log2(norm+1) like v1 pipeline
    g168_matched = glds168_df.loc[common_genes]
    g168_log2 = np.log2(g168_matched + 1)

    # Get matched sample names and conditions
    matched_168_samples = [m["glds168_sample"] for m in matches]
    matched_conditions = [m["condition"] for m in matches]

    X_test_168 = g168_log2[matched_168_samples].T.values.astype(float)
    y_test_168 = np.array([1 if c == "FLT" else 0 for c in matched_conditions])
    X_test_168 = np.nan_to_num(X_test_168, nan=0.0)

    # Need to align features: training used full_gene columns, GLDS-168 uses common_genes
    # Rebuild with common feature set
    train_gene_cols = list(full_gene_a.columns)
    # Map GLDS-168 features to training feature space
    feat_168 = pd.DataFrame(0.0, index=matched_168_samples, columns=train_gene_cols)
    for gene in common_genes:
        if gene in train_gene_cols:
            feat_168[gene] = g168_log2.loc[gene, matched_168_samples].values

    X_test_168_full = feat_168.values.astype(float)
    X_test_168_s = scaler.transform(X_test_168_full)
    X_test_168_p = pca.transform(X_test_168_s)

    if len(np.unique(y_test_168)) >= 2:
        proba_168 = clf.predict_proba(X_test_168_p)[:, list(clf.classes_).index(1)]
        auroc_168 = roc_auc_score(y_test_168, proba_168)
        print(f"  GLDS-168 (reprocessed) RR-1 fold AUROC: {auroc_168:.3f}")
    else:
        auroc_168 = None
        print(f"  GLDS-168: only one class in test set, cannot compute AUROC")

    delta = None
    if auroc_48 is not None and auroc_168 is not None:
        delta = auroc_168 - auroc_48
        print(f"  Delta (168 - 48): {delta:+.3f}")

    return {
        "glds48_auroc": auroc_48,
        "glds168_auroc": auroc_168,
        "delta": delta,
        "n_train": len(y_train),
        "n_test_48": len(y_test_48),
        "n_test_168": len(y_test_168),
        "n_pca_components": n_comp,
    }


def compare_gene_selection(glds48_df, glds168_df, common_genes):
    """Compare top variable genes between pipelines."""
    print(f"\n  [Gene selection comparison]")

    # Log2 transform
    g48_log = np.log2(glds48_df.loc[common_genes] + 1)
    g168_log = np.log2(glds168_df.loc[common_genes] + 1)

    # Variance per gene
    var_48 = g48_log.var(axis=1).sort_values(ascending=False)
    var_168 = g168_log.var(axis=1).sort_values(ascending=False)

    # Top N genes
    for n in [100, 500, 1000, 5000]:
        top48 = set(var_48.head(n).index)
        top168 = set(var_168.head(n).index)
        jaccard = len(top48 & top168) / len(top48 | top168)
        overlap = len(top48 & top168) / n
        print(f"    Top {n}: Jaccard={jaccard:.3f}, overlap={overlap:.3f}")

    # Variance rank correlation
    rank_48 = var_48.rank(ascending=False)
    rank_168 = var_168.rank(ascending=False)
    common_ranked = sorted(set(rank_48.index) & set(rank_168.index))
    rho, p = stats.spearmanr(
        rank_48.loc[common_ranked].values,
        rank_168.loc[common_ranked].values
    )
    print(f"  Variance rank correlation (Spearman): rho={rho:.4f}, p={p:.2e}")

    return {
        "jaccard_top100": len(set(var_48.head(100).index) & set(var_168.head(100).index)) / len(set(var_48.head(100).index) | set(var_168.head(100).index)),
        "jaccard_top500": len(set(var_48.head(500).index) & set(var_168.head(500).index)) / len(set(var_48.head(500).index) | set(var_168.head(500).index)),
        "jaccard_top1000": len(set(var_48.head(1000).index) & set(var_168.head(1000).index)) / len(set(var_48.head(1000).index) | set(var_168.head(1000).index)),
        "variance_rank_spearman_rho": rho,
        "variance_rank_spearman_p": p,
    }


def main():
    print("=" * 70)
    print("J1: Pipeline Version Comparison (GLDS-48 vs GLDS-168)")
    print("  RR-1 liver, C57BL/6J, CAGES/noERCC subset")
    print("=" * 70)

    J1_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Load data
    glds48 = load_glds48_counts()
    glds168 = load_glds168_counts()

    # Match samples
    matches = match_samples(glds48, glds168)
    print(f"\n  Matched samples: {len(matches)}")
    for m in matches:
        print(f"    {m['animal_id']} ({m['condition']}): "
              f"{m['glds48_sample']} ↔ {m['glds168_sample']}")

    # 1. Expression correlation
    correlations, common_genes = compare_expression(glds48, glds168, matches)

    # 2. Gene selection comparison
    gene_selection = compare_gene_selection(glds48, glds168, common_genes)

    # 3. LOMO performance comparison
    lomo_results = compare_lomo_performance(glds48, glds168, matches, common_genes)

    # Build results
    results = {
        "task": "J1",
        "description": "Pipeline version comparison: GLDS-48 (v1) vs GLDS-168 (v2)",
        "generated_at": datetime.now().isoformat(),
        "data": {
            "glds48": {
                "n_genes": glds48.shape[0],
                "n_samples": glds48.shape[1],
                "samples": list(glds48.columns),
            },
            "glds168": {
                "n_genes": glds168.shape[0],
                "n_samples": glds168.shape[1],
                "samples": list(glds168.columns),
            },
            "n_common_genes": len(common_genes),
            "n_matched_samples": len(matches),
            "matched_animals": matches,
        },
        "expression_correlation": {
            "per_sample": correlations,
            "mean_pearson_r": np.mean([c["pearson_r"] for c in correlations]),
        },
        "gene_selection": gene_selection,
        "ml_performance": lomo_results,
    }

    # Conclusion
    mean_r = results["expression_correlation"]["mean_pearson_r"]
    if mean_r > 0.99:
        conclusion = "STRONGLY_SUPPORTED"
        note = "Pipeline versions produce nearly identical expression values"
    elif mean_r > 0.95:
        conclusion = "PARTIALLY_SUPPORTED"
        note = "High correlation but some differences"
    else:
        conclusion = "NOT_SUPPORTED"
        note = "Significant differences between pipeline versions"

    results["conclusion"] = {
        "verdict": conclusion,
        "note": note,
        "criteria": {
            "expression_correlation": f"mean r = {mean_r:.4f} (threshold: 0.95)",
            "auroc_delta": f"{lomo_results.get('delta', 'N/A')}",
        }
    }

    print(f"\n  Conclusion: {conclusion} — {note}")

    outpath = RESULTS_DIR / "J1_pipeline_comparison.json"
    outpath.write_text(json.dumps(results, indent=2, cls=_NumpyEncoder))
    print(f"  Saved: {outpath}")


if __name__ == "__main__":
    main()
