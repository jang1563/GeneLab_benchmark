#!/usr/bin/env python3
"""
power_analysis.py — GeneLabBench v4 Phase 3: Power analysis for spaceflight detection

D1. Analytical power (Obuchowski-style) for required sample size
D2. Simulation-based power accounting for between-mission heterogeneity
D3. Sample recommendations per tissue

Reads pooled AUROC + tau² from META_forest_plots.json (meta_analysis.py output).

Usage:
  python power_analysis.py
  python power_analysis.py --meta-json META_forest_plots.json --output META_power_analysis.json
"""

import json
import sys
import argparse
import warnings
import numpy as np
from pathlib import Path
from scipy.stats import norm, wilcoxon, ttest_1samp
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).resolve().parent))

from v4_utils import LOMO_TISSUES, V4_EVAL_DIR


# ─────────────────────────────────────────────────────────────────────────────
# D1. Analytical Power
# ─────────────────────────────────────────────────────────────────────────────

def required_n_per_group(target_auroc, alpha=0.05, power=0.80):
    """Min samples per group (FLT/GC) for AUROC significantly > 0.5.

    Uses Obuchowski (1997) / Hanley-McNeil approximation for Wilcoxon test.
    Iterative solution (n appears on both sides of the equation).

    Args:
        target_auroc: expected true AUC (effect size)
        alpha: significance level (one-sided, H1: AUC > 0.5)
        power: target power (1 - beta)

    Returns:
        int: minimum n per group, or 10000 if not achievable
    """
    # One-sided test: H0: AUC=0.5 vs H1: AUC>0.5 (matches simulation)
    z_alpha = norm.ppf(1 - alpha)
    theta = np.clip(target_auroc, 0.501, 0.999)

    Q1 = theta / (2.0 - theta)
    Q2 = 2.0 * theta ** 2 / (1.0 + theta)

    for n in range(3, 10001):
        # Variance under H1 (Hanley-McNeil, equal groups n1=n0=n)
        # For analytical power estimation, equal groups is the standard assumption
        # since we're computing REQUIRED n. Real data group imbalance is handled
        # by the simulation-based power analysis below.
        n1, n0 = n, n
        var_h1 = (theta * (1 - theta)
                  + (n1 - 1) * (Q1 - theta ** 2)
                  + (n0 - 1) * (Q2 - theta ** 2)) / (n1 * n0)

        # z-score for the observed effect
        z_score = (theta - 0.5) / np.sqrt(var_h1)

        # Power = P(reject H0 | H1 true)
        actual_power = norm.cdf(z_score - z_alpha)

        if actual_power >= power:
            return n

    return 10000


def analytical_power_table(observed_aurocs):
    """Compute required n for each tissue at multiple power levels.

    Args:
        observed_aurocs: dict[tissue] → observed pooled AUROC

    Returns:
        dict[tissue] → required_n_per_group at 80/90/95% power
    """
    power_levels = [0.80, 0.90, 0.95]
    results = {}

    for tissue, auroc in observed_aurocs.items():
        required = {}
        for pwr in power_levels:
            n = required_n_per_group(auroc, power=pwr)
            required[f"power_{int(pwr*100)}"] = n
        results[tissue] = {
            "observed_pooled_auroc": round(auroc, 4),
            "required_n_per_group": required,
        }

    return results


# ─────────────────────────────────────────────────────────────────────────────
# D2. Simulation-Based Power
# ─────────────────────────────────────────────────────────────────────────────

def simulate_power(pooled_auroc, tau2, n_per_group, n_missions,
                   n_sim=5000, seed=42):
    """Simulate LOMO-CV power accounting for between-mission heterogeneity.

    Model:
        - Each mission has a true AUC ~ N(pooled_auroc, tau2)
        - Within each mission, FLT ~ N(mu, 1) and GC ~ N(0, 1)
          where mu = Phi^{-1}(AUC) * sqrt(2)
        - Significance: one-sample test on fold AUROCs vs 0.5

    Args:
        pooled_auroc: overall effect size
        tau2: between-mission heterogeneity variance
        n_per_group: number of samples per group per mission
        n_missions: number of missions (= number of LOMO folds)
        n_sim: number of simulations
        seed: random seed

    Returns:
        float: estimated power (fraction of significant simulations)
    """
    rng = np.random.default_rng(seed)
    significant = 0

    for _ in range(n_sim):
        # Sample mission-level true AUROCs
        mission_aurocs = rng.normal(pooled_auroc, np.sqrt(max(tau2, 1e-6)), n_missions)
        mission_aurocs = np.clip(mission_aurocs, 0.01, 0.99)

        # Per mission: simulate binary classification scores
        fold_aurocs = []
        for true_auc in mission_aurocs:
            mu = norm.ppf(true_auc) * np.sqrt(2)
            scores_flt = rng.normal(mu, 1, n_per_group)
            scores_gc = rng.normal(0, 1, n_per_group)
            y_true = np.array([1] * n_per_group + [0] * n_per_group)
            y_score = np.concatenate([scores_flt, scores_gc])

            try:
                auc = roc_auc_score(y_true, y_score)
            except ValueError:
                auc = 0.5
            fold_aurocs.append(auc)

        # Test significance
        fold_arr = np.array(fold_aurocs)
        try:
            if n_missions >= 5:
                _, p_val = wilcoxon(fold_arr - 0.5, alternative="greater")
            else:
                _, p_val = ttest_1samp(fold_arr, 0.5, alternative="greater")
        except (ValueError, ZeroDivisionError):
            p_val = 1.0

        if p_val < 0.05:
            significant += 1

    return significant / n_sim


def simulation_power_grid(tissue_params, n_per_group_grid=None,
                          n_missions_grid=None, n_sim=5000):
    """Run power simulation across parameter grid.

    Args:
        tissue_params: dict[tissue] → {pooled_auroc, tau2, current_k}
        n_per_group_grid: list of n_per_group values
        n_missions_grid: list of n_missions values

    Returns:
        dict[tissue] → power_curve list
    """
    if n_per_group_grid is None:
        n_per_group_grid = [3, 5, 8, 10, 15, 20, 30, 50]
    if n_missions_grid is None:
        n_missions_grid = [3, 4, 6, 8, 10, 15]

    results = {}

    for tissue, params in tissue_params.items():
        auroc = params["pooled_auroc"]
        tau2 = params["tau2"]
        current_k = params["current_k"]

        print(f"  {tissue}: AUC={auroc:.3f}, tau2={tau2:.4f}, k={current_k}")

        # Power curve at current number of missions
        power_at_current_k = []
        for n_pg in n_per_group_grid:
            pwr = simulate_power(auroc, tau2, n_pg, current_k, n_sim=n_sim)
            power_at_current_k.append({
                "n_per_group": n_pg,
                "n_missions": current_k,
                "power": round(pwr, 3),
            })

        # Power curve varying missions (fixed n_per_group=10)
        power_at_n10 = []
        for n_m in n_missions_grid:
            pwr = simulate_power(auroc, tau2, 10, n_m, n_sim=n_sim)
            power_at_n10.append({
                "n_per_group": 10,
                "n_missions": n_m,
                "power": round(pwr, 3),
            })

        results[tissue] = {
            "pooled_auroc": round(auroc, 4),
            "tau2": round(tau2, 6),
            "current_k": current_k,
            "power_curve_vary_n": power_at_current_k,
            "power_curve_vary_missions": power_at_n10,
        }

    return results


# ─────────────────────────────────────────────────────────────────────────────
# D3. Recommendations
# ─────────────────────────────────────────────────────────────────────────────

def generate_recommendations(analytical, simulation):
    """Classify tissues by power status and generate recommendations."""
    well_powered = []
    moderate = []
    underpowered = []

    for tissue in analytical:
        auroc = analytical[tissue]["observed_pooled_auroc"]
        n80 = analytical[tissue]["required_n_per_group"]["power_80"]

        if n80 <= 10:
            well_powered.append(tissue)
        elif n80 <= 25:
            moderate.append(tissue)
        else:
            underpowered.append(tissue)

    # Generate recommendation text
    recs = []
    if well_powered:
        recs.append(f"Well-powered tissues ({', '.join(well_powered)}): "
                    f"Current sample sizes sufficient for reliable detection.")
    if moderate:
        recs.append(f"Moderate power ({', '.join(moderate)}): "
                    f"Adding 1-2 missions would substantially improve detection reliability.")
    if underpowered:
        recs.append(f"Underpowered tissues ({', '.join(underpowered)}): "
                    f"Effect size is modest; require ≥20 samples/group or pathway features.")

    return {
        "well_powered": well_powered,
        "moderate": moderate,
        "underpowered": underpowered,
        "recommendation_text": " ".join(recs),
    }


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Power analysis for spaceflight detection")
    parser.add_argument("--meta-json", type=str, default=None,
                        help="META_forest_plots.json path")
    parser.add_argument("--output", type=str, default=None,
                        help="Output JSON path")
    parser.add_argument("--n-sim", type=int, default=5000,
                        help="Number of simulations per grid point (default: 5000)")
    args = parser.parse_args()

    meta_path = Path(args.meta_json) if args.meta_json else V4_EVAL_DIR / "META_forest_plots.json"
    output_path = Path(args.output) if args.output else V4_EVAL_DIR / "META_power_analysis.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Load meta-analysis results
    with open(meta_path) as f:
        meta = json.load(f)

    by_tissue = meta.get("by_tissue", {})

    # Extract tissue parameters (gene features, PCA-LR as reference method)
    tissue_params = {}
    observed_aurocs = {}

    for tissue in LOMO_TISSUES:
        key = "gene_pca_lr"
        if tissue in by_tissue and key in by_tissue[tissue]:
            r = by_tissue[tissue][key]
            auroc = r["pooled_re"]["auroc"]
            tau2 = r["heterogeneity"]["tau2"]
            k = r["k"]
            tissue_params[tissue] = {
                "pooled_auroc": auroc,
                "tau2": tau2,
                "current_k": k,
            }
            observed_aurocs[tissue] = auroc

    print(f"Power analysis for {len(tissue_params)} LOMO tissues")
    print(f"Simulations per grid point: {args.n_sim}")

    # D1: Analytical power
    print("\n─── Analytical Power ───")
    analytical = analytical_power_table(observed_aurocs)
    for tissue, data in sorted(analytical.items(), key=lambda x: x[1]["observed_pooled_auroc"], reverse=True):
        n80 = data["required_n_per_group"]["power_80"]
        print(f"  {tissue:15s}: AUC={data['observed_pooled_auroc']:.3f}  "
              f"n(80%)={n80}  n(90%)={data['required_n_per_group']['power_90']}  "
              f"n(95%)={data['required_n_per_group']['power_95']}")

    # D2: Simulation-based power
    print("\n─── Simulation Power ───")
    simulation = simulation_power_grid(tissue_params, n_sim=args.n_sim)

    # D3: Recommendations
    recommendations = generate_recommendations(analytical, simulation)

    # Assemble output
    output = {
        "analytical": analytical,
        "simulation": simulation,
        "recommendations": recommendations,
        "parameters": {
            "reference_method": "pca_lr",
            "feature_type": "gene",
            "n_sim": args.n_sim,
            "alpha": 0.05,
        },
    }

    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\nOutput: {output_path}")

    # Print recommendations
    print(f"\n{'='*50}")
    print("Recommendations:")
    print(f"  Well-powered: {recommendations['well_powered']}")
    print(f"  Moderate: {recommendations['moderate']}")
    print(f"  Underpowered: {recommendations['underpowered']}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
