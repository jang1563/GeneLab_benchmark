#!/usr/bin/env python3
"""
method_ranking.py — GeneLabBench v4 Phase 3: Method ranking + pairwise comparison

1. Friedman + Nemenyi post-hoc (6 LOMO + 8 all tissues)
2. DeLong pairwise (concatenated predictions per tissue)
3. Wilcoxon signed-rank (per-fold AUROCs)

Usage:
  python method_ranking.py
  python method_ranking.py --m1-summary M1_summary.json --predictions META_predictions.json
"""

import json
import sys
import argparse
import warnings
import numpy as np
from pathlib import Path
from itertools import combinations
from scipy.stats import friedmanchisquare, rankdata, wilcoxon, norm
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parent))

from v4_utils import TISSUE_MISSIONS, LOMO_TISSUES, V4_EVAL_DIR

ALL_TISSUES = list(TISSUE_MISSIONS.keys())

# Method keys (must match classifier_registry.py)
METHOD_KEYS = [
    "elasticnet_lr", "pca_lr", "rf", "xgb",
    "svm_rbf", "knn", "mlp", "tabnet"
]

# Nemenyi q_alpha critical values (alpha=0.05, from Demšar 2006)
NEMENYI_Q = {
    2: 1.960, 3: 2.343, 4: 2.569, 5: 2.728,
    6: 2.850, 7: 2.949, 8: 3.031, 9: 3.102, 10: 3.164
}


# ─────────────────────────────────────────────────────────────────────────────
# C1. Friedman + Nemenyi
# ─────────────────────────────────────────────────────────────────────────────

def friedman_nemenyi(auroc_matrix, tissue_names, method_names):
    """Friedman test + Nemenyi post-hoc on AUROC matrix.

    Args:
        auroc_matrix: shape (n_tissues, n_methods)
        tissue_names: list of tissue names
        method_names: list of method keys

    Returns:
        dict with chi2, p_value, avg_ranks, CD, significant_pairs, groups
    """
    n_tissues, k_methods = auroc_matrix.shape
    assert k_methods == len(method_names)

    # Friedman test
    try:
        stat, p_value = friedmanchisquare(
            *[auroc_matrix[:, j] for j in range(k_methods)]
        )
    except Exception as e:
        return {"error": str(e), "chi2": np.nan, "p_value": 1.0}

    # Average ranks per method (lower rank = better = higher AUROC)
    ranks = np.zeros_like(auroc_matrix)
    for i in range(n_tissues):
        ranks[i] = rankdata(-auroc_matrix[i])  # negative → higher AUROC gets rank 1
    avg_ranks = ranks.mean(axis=0)

    # Nemenyi critical difference
    q_alpha = NEMENYI_Q.get(k_methods, 3.031)
    cd = q_alpha * np.sqrt(k_methods * (k_methods + 1) / (6.0 * n_tissues))

    # Significant pairs (rank difference > CD)
    significant_pairs = []
    for i, j in combinations(range(k_methods), 2):
        if abs(avg_ranks[i] - avg_ranks[j]) > cd:
            significant_pairs.append([method_names[i], method_names[j]])

    # Group methods using maximal cliques (overlapping groups allowed)
    # Standard Nemenyi: a "group" is a maximal set of methods where ALL pairs
    # have rank difference <= CD. Methods can belong to multiple groups.
    sorted_indices = list(np.argsort(avg_ranks))

    # Find all maximal cliques of methods within CD of each other
    # Since methods are sorted by rank, cliques are contiguous runs where
    # first and last are within CD (rank-sorted property of 1D data)
    groups = []
    for start in range(len(sorted_indices)):
        clique = [sorted_indices[start]]
        for end in range(start + 1, len(sorted_indices)):
            idx = sorted_indices[end]
            # In 1D sorted data, if first-last <= CD, all pairs are within CD
            if avg_ranks[idx] - avg_ranks[sorted_indices[start]] <= cd:
                clique.append(idx)
            else:
                break
        # Keep only maximal cliques (not subsets of existing groups)
        if not any(set(clique).issubset(set(g)) for g in groups):
            groups.append(clique)

    # Remove non-maximal groups (subsets of larger groups)
    maximal_groups = []
    for g in groups:
        if not any(set(g) < set(other) for other in groups):
            maximal_groups.append(g)
    groups = maximal_groups

    nemenyi_groups = []
    for gi, group in enumerate(groups):
        nemenyi_groups.append({
            "group": chr(65 + gi),  # A, B, C, ...
            "methods": [method_names[i] for i in group],
            "avg_ranks": [round(float(avg_ranks[i]), 3) for i in group],
        })

    # Per-tissue ranks for detail
    per_tissue_ranks = {}
    for ti, tissue in enumerate(tissue_names):
        per_tissue_ranks[tissue] = {
            method_names[j]: int(ranks[ti, j]) for j in range(k_methods)
        }

    return {
        "chi2": round(float(stat), 3),
        "p_value": round(float(p_value), 6),
        "n_tissues": n_tissues,
        "k_methods": k_methods,
        "avg_ranks": {method_names[j]: round(float(avg_ranks[j]), 3) for j in range(k_methods)},
        "critical_difference": round(float(cd), 3),
        "q_alpha": q_alpha,
        "significant_pairs": significant_pairs,
        "nemenyi_groups": nemenyi_groups,
        "per_tissue_ranks": per_tissue_ranks,
    }


# ─────────────────────────────────────────────────────────────────────────────
# C2. DeLong Test
# ─────────────────────────────────────────────────────────────────────────────

def delong_test(y_true, y_score_a, y_score_b):
    """Paired DeLong test for two ROC curves on same samples.

    Uses placement values (structural components).
    O(n²) complexity — fine for n < 300.

    Returns: z_stat, p_value (two-sided), auc_a, auc_b
    """
    y_true = np.asarray(y_true, dtype=int)
    y_score_a = np.asarray(y_score_a, dtype=float)
    y_score_b = np.asarray(y_score_b, dtype=float)

    n1 = np.sum(y_true == 1)
    n0 = np.sum(y_true == 0)

    if n1 < 2 or n0 < 2:
        return 0.0, 1.0, np.nan, np.nan

    # Separate scores by class
    pos_a = y_score_a[y_true == 1]
    neg_a = y_score_a[y_true == 0]
    pos_b = y_score_b[y_true == 1]
    neg_b = y_score_b[y_true == 0]

    # Placement values V10 (positive placements)
    V10_a = np.array([np.mean(s > neg_a) + 0.5 * np.mean(s == neg_a) for s in pos_a])
    V10_b = np.array([np.mean(s > neg_b) + 0.5 * np.mean(s == neg_b) for s in pos_b])

    # Placement values V01 (negative placements)
    V01_a = np.array([np.mean(pos_a > s) + 0.5 * np.mean(pos_a == s) for s in neg_a])
    V01_b = np.array([np.mean(pos_b > s) + 0.5 * np.mean(pos_b == s) for s in neg_b])

    # AUCs
    auc_a = float(np.mean(V10_a))
    auc_b = float(np.mean(V10_b))

    # Covariance matrix of (AUC_a, AUC_b)
    # S10: covariance from positive class placements
    # S01: covariance from negative class placements
    V10 = np.stack([V10_a, V10_b])  # 2 × n1
    V01 = np.stack([V01_a, V01_b])  # 2 × n0

    if n1 < 2 or n0 < 2:
        return 0.0, 1.0, auc_a, auc_b

    S10 = np.cov(V10)  # 2×2
    S01 = np.cov(V01)  # 2×2

    S = S10 / n0 + S01 / n1

    # z-statistic for AUC difference
    diff = auc_a - auc_b
    var_diff = S[0, 0] + S[1, 1] - 2 * S[0, 1]

    if var_diff <= 0:
        return 0.0, 1.0, auc_a, auc_b

    z = diff / np.sqrt(var_diff)
    p = 2.0 * (1.0 - norm.cdf(abs(z)))

    return float(z), float(p), auc_a, auc_b


def run_delong_pairwise(predictions, tissues, method_keys):
    """Run DeLong test for all method pairs on each tissue.

    Concatenates fold predictions for each method pair.

    Returns:
        dict[tissue][pair_name] = {z, p_raw, p_fdr, auc_diff, ...}
    """
    results = {}

    for tissue in tissues:
        if tissue not in predictions:
            continue

        tissue_preds = predictions[tissue]
        pair_results = {}
        pair_names = []
        p_values_raw = []

        for m_a, m_b in combinations(method_keys, 2):
            if m_a not in tissue_preds or m_b not in tissue_preds:
                continue

            # Align folds — only use folds present in BOTH methods
            folds_a = {f["fold_name"]: f for f in tissue_preds[m_a]}
            folds_b = {f["fold_name"]: f for f in tissue_preds[m_b]}
            common_folds = sorted(set(folds_a.keys()) & set(folds_b.keys()))

            if len(common_folds) < 1:
                continue

            # Concatenate predictions from common folds
            y_true_all = np.concatenate([np.array(folds_a[fn]["y_true"]) for fn in common_folds])
            y_score_a_all = np.concatenate([np.array(folds_a[fn]["y_score"]) for fn in common_folds])
            y_score_b_all = np.concatenate([np.array(folds_b[fn]["y_score"]) for fn in common_folds])

            # Verify alignment
            y_true_b = np.concatenate([np.array(folds_b[fn]["y_true"]) for fn in common_folds])
            if not np.array_equal(y_true_all, y_true_b):
                warnings.warn(f"y_true mismatch for {tissue}/{m_a} vs {m_b}")
                continue

            z, p, auc_a, auc_b = delong_test(y_true_all, y_score_a_all, y_score_b_all)

            pair_name = f"{m_a}_vs_{m_b}"
            pair_results[pair_name] = {
                "z": round(z, 4),
                "p_raw": round(p, 6),
                "auc_a": round(auc_a, 4),
                "auc_b": round(auc_b, 4),
                "auc_diff": round(auc_a - auc_b, 4),
                "n_common_folds": len(common_folds),
                "n_samples": len(y_true_all),
            }
            pair_names.append(pair_name)
            p_values_raw.append(p)

        # BH-FDR correction across all pairs within this tissue
        if p_values_raw:
            _, p_fdr, _, _ = multipletests(p_values_raw, method="fdr_bh")
            for i, pn in enumerate(pair_names):
                pair_results[pn]["p_fdr"] = round(float(p_fdr[i]), 6)

        results[tissue] = pair_results

    return results


# ─────────────────────────────────────────────────────────────────────────────
# C3. Wilcoxon Signed-Rank
# ─────────────────────────────────────────────────────────────────────────────

def run_wilcoxon_pairwise(m1_data, tissues, method_keys, feature_type="gene"):
    """Wilcoxon signed-rank test on per-fold AUROCs for all method pairs.

    Uses M1 per-fold data (no raw predictions needed).

    LIMITATION: With n=3-6 folds, Wilcoxon minimum possible p-value is:
      n=3: min_p=0.250, n=4: min_p=0.125, n=5: min_p=0.0625, n=6: min_p=0.0625
    After BH-FDR correction over C(8,2)=28 pairs, NO comparison can reach
    significance at alpha=0.05. Results are included for completeness but
    should NOT be interpreted as evidence of method equivalence.
    DeLong test (on concatenated predictions) provides more power.

    Returns:
        dict[tissue][pair_name] = {stat, p_raw, p_fdr, n_folds}
    """
    results = {}

    for tissue in tissues:
        pair_results = {}
        pair_names = []
        p_values_raw = []

        for m_a, m_b in combinations(method_keys, 2):
            # Load per-fold AUROCs from M1 JSONs
            aurocs_a = _get_fold_aurocs(m1_data, tissue, feature_type, m_a)
            aurocs_b = _get_fold_aurocs(m1_data, tissue, feature_type, m_b)

            if aurocs_a is None or aurocs_b is None:
                continue

            # Align by fold name
            common = sorted(set(aurocs_a.keys()) & set(aurocs_b.keys()))
            if len(common) < 3:
                continue

            a_vals = np.array([aurocs_a[fn] for fn in common])
            b_vals = np.array([aurocs_b[fn] for fn in common])

            # Wilcoxon signed-rank test
            try:
                stat, p = wilcoxon(a_vals, b_vals, alternative="two-sided")
            except ValueError:
                # All differences are zero
                stat, p = 0.0, 1.0

            pair_name = f"{m_a}_vs_{m_b}"
            pair_results[pair_name] = {
                "stat": round(float(stat), 4),
                "p_raw": round(float(p), 6),
                "n_folds": len(common),
                "mean_diff": round(float(np.mean(a_vals - b_vals)), 4),
                "low_n_warning": len(common) < 6,
            }
            pair_names.append(pair_name)
            p_values_raw.append(p)

        # BH-FDR correction
        if p_values_raw:
            _, p_fdr, _, _ = multipletests(p_values_raw, method="fdr_bh")
            for i, pn in enumerate(pair_names):
                pair_results[pn]["p_fdr"] = round(float(p_fdr[i]), 6)

        results[tissue] = pair_results

    return results


def _get_fold_aurocs(m1_data, tissue, feature_type, method_key):
    """Extract {fold_name: auroc} from loaded M1 JSONs."""
    # Try loading individual M1 JSON
    m1_path = V4_EVAL_DIR / f"M1_{tissue}_{feature_type}_{method_key}.json"
    if m1_path.exists():
        with open(m1_path) as f:
            data = json.load(f)
        return {fold["fold_name"]: fold["auroc"] for fold in data.get("folds", [])}

    return None


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Method ranking + pairwise comparison")
    parser.add_argument("--m1-summary", type=str, default=None,
                        help="M1_summary.json path")
    parser.add_argument("--predictions", type=str, default=None,
                        help="META_predictions.json path")
    parser.add_argument("--output", type=str, default=None,
                        help="Output JSON path")
    args = parser.parse_args()

    m1_summary_path = Path(args.m1_summary) if args.m1_summary else V4_EVAL_DIR / "M1_summary.json"
    pred_path = Path(args.predictions) if args.predictions else V4_EVAL_DIR / "META_predictions.json"
    output_path = Path(args.output) if args.output else V4_EVAL_DIR / "META_method_ranking.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Load M1 summary
    with open(m1_summary_path) as f:
        m1_summary = json.load(f)

    # ── Friedman + Nemenyi ──

    # Build AUROC matrix from M1 summary (gene features)
    def build_auroc_matrix(tissues):
        matrix = np.full((len(tissues), len(METHOD_KEYS)), np.nan)
        for ti, tissue in enumerate(tissues):
            for mi, method in enumerate(METHOD_KEYS):
                try:
                    matrix[ti, mi] = m1_summary[tissue]["gene"][method]["auroc"]
                except (KeyError, TypeError):
                    pass
        return matrix

    # Primary: 6 LOMO tissues
    auroc_lomo = build_auroc_matrix(LOMO_TISSUES)
    # Drop any tissue/method with NaN
    valid_lomo = ~np.any(np.isnan(auroc_lomo), axis=1)
    lomo_tissues_valid = [t for t, v in zip(LOMO_TISSUES, valid_lomo) if v]
    auroc_lomo_valid = auroc_lomo[valid_lomo]

    print(f"Friedman (LOMO-6): {len(lomo_tissues_valid)} tissues × {len(METHOD_KEYS)} methods")
    friedman_lomo6 = friedman_nemenyi(auroc_lomo_valid, lomo_tissues_valid, METHOD_KEYS)

    # Supplementary: all 8 tissues
    auroc_all = build_auroc_matrix(ALL_TISSUES)
    valid_all = ~np.any(np.isnan(auroc_all), axis=1)
    all_tissues_valid = [t for t, v in zip(ALL_TISSUES, valid_all) if v]
    auroc_all_valid = auroc_all[valid_all]

    print(f"Friedman (All-8): {len(all_tissues_valid)} tissues × {len(METHOD_KEYS)} methods")
    friedman_all8 = friedman_nemenyi(auroc_all_valid, all_tissues_valid, METHOD_KEYS)

    # ── DeLong pairwise ──

    predictions = None
    delong_results = {}
    if pred_path.exists():
        with open(pred_path) as f:
            predictions = json.load(f)
        print(f"\nDeLong pairwise: {len(ALL_TISSUES)} tissues × C(8,2)=28 pairs")
        delong_results = run_delong_pairwise(predictions, ALL_TISSUES, METHOD_KEYS)
    else:
        print(f"\nWARNING: {pred_path} not found — skipping DeLong test")

    # ── Wilcoxon pairwise ──

    print(f"Wilcoxon pairwise: {len(ALL_TISSUES)} tissues")
    wilcoxon_results = run_wilcoxon_pairwise(None, ALL_TISSUES, METHOD_KEYS, "gene")

    # ── Assemble output ──

    output = {
        "friedman_lomo6": friedman_lomo6,
        "friedman_all8": friedman_all8,
        "delong_pairwise": delong_results,
        "wilcoxon_pairwise": wilcoxon_results,
        "wilcoxon_limitation": (
            "Wilcoxon signed-rank with n=3-6 folds cannot reach significance "
            "after BH-FDR correction over 28 pairs (min adjusted p >= 0.25). "
            "Included for completeness; DeLong test provides meaningful pairwise comparison."
        ),
        "methods": METHOD_KEYS,
        "feature_type": "gene",
    }

    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\nOutput: {output_path}")

    # Summary
    print(f"\n{'='*50}")
    print("Friedman (LOMO-6):")
    print(f"  chi2={friedman_lomo6.get('chi2', 'N/A')}, "
          f"p={friedman_lomo6.get('p_value', 'N/A')}")
    if "avg_ranks" in friedman_lomo6:
        sorted_ranks = sorted(friedman_lomo6["avg_ranks"].items(), key=lambda x: x[1])
        for method, rank in sorted_ranks:
            print(f"  {method:18s}: rank {rank:.2f}")
    print(f"  CD={friedman_lomo6.get('critical_difference', 'N/A')}")
    print(f"  Significant pairs: {len(friedman_lomo6.get('significant_pairs', []))}")

    if delong_results:
        n_sig = sum(
            1 for tissue_pairs in delong_results.values()
            for pair in tissue_pairs.values()
            if pair.get("p_fdr", 1.0) < 0.05
        )
        print(f"\nDeLong significant pairs (FDR<0.05): {n_sig}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
