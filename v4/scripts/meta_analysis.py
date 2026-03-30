#!/usr/bin/env python3
"""
meta_analysis.py — GeneLabBench v4 Phase 3: DerSimonian-Laird random-effects meta-analysis

Forest plots + I² heterogeneity + Cochran's Q for LOMO tissues.
Reads per-fold AUROC + n_flight_test + n_ground_test from M1 JSONs.

Design decisions:
  DD-26: Random-effects model, I²>75% = high, report both FE and RE
  DD-30: 6 LOMO tissues only (lung/colon excluded)

Usage:
  python meta_analysis.py                                          # Default
  python meta_analysis.py --input-dir v4/evaluation/ --output META_forest_plots.json
"""

import json
import sys
import argparse
import numpy as np
from pathlib import Path
from scipy.stats import chi2

sys.path.insert(0, str(Path(__file__).resolve().parent))

from v4_utils import LOMO_TISSUES, V4_EVAL_DIR

# Feature types to analyze
FEATURE_TYPES = ["gene", "pathway_hallmark", "pathway_kegg", "combined"]

# Methods (matching classifier_registry.py keys)
METHOD_KEYS = [
    "elasticnet_lr", "pca_lr", "rf", "xgb",
    "svm_rbf", "knn", "mlp", "tabnet"
]


def hanley_mcneil_var(auc, n1, n0):
    """Variance of AUC under Hanley-McNeil (1982) approximation.

    Args:
        auc: observed AUROC
        n1: number of positives (FLT)
        n0: number of negatives (GC)

    Returns:
        Estimated variance of AUC
    """
    if n1 < 1 or n0 < 1:
        return 1.0  # degenerate case
    # Clamp AUC away from exact 0 or 1 to avoid division issues
    auc = np.clip(auc, 0.001, 0.999)
    Q1 = auc / (2.0 - auc)
    Q2 = 2.0 * auc ** 2 / (1.0 + auc)
    var = (auc * (1 - auc)
           + (n1 - 1) * (Q1 - auc ** 2)
           + (n0 - 1) * (Q2 - auc ** 2)) / (n1 * n0)
    return max(var, 1e-10)


def dersimonian_laird(effects, variances):
    """DerSimonian-Laird random-effects meta-analysis.

    Args:
        effects: array of per-fold AUROCs (k values)
        variances: array of Hanley-McNeil variances (k values)

    Returns:
        dict with pooled estimates, heterogeneity stats
    """
    effects = np.asarray(effects, dtype=float)
    variances = np.asarray(variances, dtype=float)
    k = len(effects)

    if k < 2:
        return {
            "pooled_auroc_re": float(effects[0]) if k == 1 else np.nan,
            "se_re": np.nan, "ci_lower_re": np.nan, "ci_upper_re": np.nan,
            "pooled_auroc_fe": float(effects[0]) if k == 1 else np.nan,
            "tau2": 0.0, "I2": 0.0, "Q": 0.0, "Q_df": 0, "Q_pvalue": 1.0,
            "low_k_warning": True
        }

    # Fixed-effect weights
    w_fe = 1.0 / variances

    # Fixed-effect estimate
    mu_fe = np.sum(w_fe * effects) / np.sum(w_fe)
    se_fe = 1.0 / np.sqrt(np.sum(w_fe))

    # Cochran's Q
    Q = float(np.sum(w_fe * (effects - mu_fe) ** 2))
    Q_df = k - 1
    Q_pvalue = float(1.0 - chi2.cdf(Q, Q_df)) if Q_df > 0 else 1.0

    # Between-study variance (tau²)
    C = np.sum(w_fe) - np.sum(w_fe ** 2) / np.sum(w_fe)
    tau2 = max(0.0, (Q - Q_df) / C) if C > 0 else 0.0

    # Random-effects weights
    w_re = 1.0 / (variances + tau2)
    mu_re = np.sum(w_re * effects) / np.sum(w_re)
    se_re = 1.0 / np.sqrt(np.sum(w_re))

    # I²
    I2 = max(0.0, (Q - Q_df) / Q) * 100.0 if Q > 0 else 0.0

    # 95% CI
    ci_lower_re = mu_re - 1.96 * se_re
    ci_upper_re = mu_re + 1.96 * se_re
    ci_lower_fe = mu_fe - 1.96 * se_fe
    ci_upper_fe = mu_fe + 1.96 * se_fe

    return {
        "pooled_auroc_re": round(float(mu_re), 4),
        "se_re": round(float(se_re), 4),
        "ci_lower_re": round(float(ci_lower_re), 4),
        "ci_upper_re": round(float(ci_upper_re), 4),
        "pooled_auroc_fe": round(float(mu_fe), 4),
        "se_fe": round(float(se_fe), 4),
        "ci_lower_fe": round(float(ci_lower_fe), 4),
        "ci_upper_fe": round(float(ci_upper_fe), 4),
        "tau2": round(float(tau2), 6),
        "I2": round(float(I2), 1),
        "Q": round(float(Q), 3),
        "Q_df": Q_df,
        "Q_pvalue": round(float(Q_pvalue), 4),
        "low_k_warning": k < 5
    }


def build_forest_plot_data(tissue, method_key, feature_type, input_dir):
    """Build forest plot data for one tissue × method × feature_type.

    Returns dict with per-fold studies + pooled estimates + heterogeneity.
    """
    # Load M1 JSON
    m1_path = input_dir / f"M1_{tissue}_{feature_type}_{method_key}.json"
    if not m1_path.exists():
        return None

    with open(m1_path) as f:
        m1 = json.load(f)

    folds = m1.get("folds", [])
    if not folds:
        return None

    # Extract per-fold data
    effects = []
    variances = []
    studies = []

    for fold in folds:
        auroc = fold["auroc"]
        n1 = fold.get("n_flight_test", 0)
        n0 = fold.get("n_ground_test", 0)

        if n1 < 1 or n0 < 1:
            continue

        var = hanley_mcneil_var(auroc, n1, n0)
        se = np.sqrt(var)

        effects.append(auroc)
        variances.append(var)

        studies.append({
            "fold": fold.get("test_mission", fold.get("fold_name", "")),
            "fold_name": fold.get("fold_name", ""),
            "auroc": round(auroc, 4),
            "se": round(float(se), 4),
            "ci_lower": round(auroc - 1.96 * se, 4),
            "ci_upper": round(auroc + 1.96 * se, 4),
            "n_test": fold.get("n_test", n1 + n0),
            "n_flight": n1,
            "n_ground": n0,
        })

    if len(effects) < 1:
        return None

    # Run meta-analysis
    meta = dersimonian_laird(np.array(effects), np.array(variances))

    # Compute RE weights for display
    tau2 = meta["tau2"]
    for i, s in enumerate(studies):
        w_re = 1.0 / (variances[i] + tau2) if (variances[i] + tau2) > 0 else 0
        total_w = sum(1.0 / (v + tau2) for v in variances if (v + tau2) > 0)
        s["weight_re"] = round(w_re / total_w, 4) if total_w > 0 else 0

    result = {
        "tissue": tissue,
        "method": method_key,
        "feature_type": feature_type,
        "k": len(effects),
        "studies": studies,
        "pooled_re": {
            "auroc": meta["pooled_auroc_re"],
            "se": meta["se_re"],
            "ci_lower": meta["ci_lower_re"],
            "ci_upper": meta["ci_upper_re"],
        },
        "pooled_fe": {
            "auroc": meta["pooled_auroc_fe"],
            "se": meta["se_fe"],
            "ci_lower": meta["ci_lower_fe"],
            "ci_upper": meta["ci_upper_fe"],
        },
        "heterogeneity": {
            "I2": meta["I2"],
            "tau2": meta["tau2"],
            "Q": meta["Q"],
            "Q_df": meta["Q_df"],
            "Q_pvalue": meta["Q_pvalue"],
        },
        "low_k_warning": meta["low_k_warning"],
    }

    return result


def main():
    parser = argparse.ArgumentParser(
        description="DerSimonian-Laird meta-analysis for LOMO tissues"
    )
    parser.add_argument("--input-dir", type=str, default=None,
                        help="Directory with M1 JSONs (default: v4/evaluation/)")
    parser.add_argument("--output", type=str, default=None,
                        help="Output JSON path (default: META_forest_plots.json)")
    args = parser.parse_args()

    input_dir = Path(args.input_dir) if args.input_dir else V4_EVAL_DIR
    output_path = Path(args.output) if args.output else V4_EVAL_DIR / "META_forest_plots.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    all_results = []
    n_total = 0
    n_success = 0

    print(f"Meta-analysis: {len(LOMO_TISSUES)} LOMO tissues × "
          f"{len(METHOD_KEYS)} methods × {len(FEATURE_TYPES)} features")
    print(f"Input: {input_dir}")

    for tissue in LOMO_TISSUES:
        for feature_type in FEATURE_TYPES:
            for method_key in METHOD_KEYS:
                n_total += 1
                result = build_forest_plot_data(tissue, method_key, feature_type, input_dir)
                if result is not None:
                    all_results.append(result)
                    n_success += 1

    # Summary statistics
    summary = {
        "n_analyses": n_success,
        "n_attempted": n_total,
        "n_lomo_tissues": len(LOMO_TISSUES),
        "lomo_tissues": LOMO_TISSUES,
        "feature_types": FEATURE_TYPES,
        "methods": METHOD_KEYS,
    }

    # Organize by tissue for easy lookup
    by_tissue = {}
    for r in all_results:
        tissue = r["tissue"]
        if tissue not in by_tissue:
            by_tissue[tissue] = {}
        key = f"{r['feature_type']}_{r['method']}"
        by_tissue[tissue][key] = r

    output = {
        "summary": summary,
        "forest_plots": all_results,
        "by_tissue": by_tissue,
    }

    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\nCompleted: {n_success}/{n_total} analyses")
    print(f"Output: {output_path}")

    # Quick summary: I² per tissue for gene/PCA-LR
    print(f"\n{'='*50}")
    print("I² summary (gene features, PCA-LR):")
    for tissue in LOMO_TISSUES:
        key = f"gene_pca_lr"
        if tissue in by_tissue and key in by_tissue[tissue]:
            r = by_tissue[tissue][key]
            h = r["heterogeneity"]
            p = r["pooled_re"]
            warn = " [LOW K]" if r["low_k_warning"] else ""
            print(f"  {tissue:15s}: I²={h['I2']:5.1f}%  "
                  f"RE={p['auroc']:.3f} [{p['ci_lower']:.3f}, {p['ci_upper']:.3f}]  "
                  f"Q_p={h['Q_pvalue']:.3f}  k={r['k']}{warn}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
