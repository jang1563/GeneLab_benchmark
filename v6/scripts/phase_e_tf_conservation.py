#!/usr/bin/env python3
"""
Phase E: TF Activity Conservation — Mouse vs Human cfRNA TF activity

Compare mouse TF activity (decoupler, CollecTRI mouse) with human cfRNA
TF activity (decoupler, CollecTRI human) inferred from pre_vs_flight statistics.
"""

import sys
import os
import numpy as np
import pandas as pd
from scipy import stats
from datetime import datetime

sys.path.insert(0, os.path.dirname(__file__))
from v6_utils import (
    load_cfrna_de, load_tf_activity, load_ortholog_map,
    ALL_TISSUES, V6_EVAL, save_json
)


def infer_human_tf_activity(cfrna_de):
    """Infer TF activity from cfRNA gene-level statistics using decoupler.

    Uses edge_pre_vs_flight_diff as gene-level statistic.
    Runs ULM (Univariate Linear Model) with CollecTRI human regulons.
    """
    import decoupler as dc

    # Get CollecTRI human regulons
    print("  Loading CollecTRI human regulons...")
    net = dc.get_collectri(organism="human")
    print(f"  CollecTRI: {len(net)} interactions, "
          f"{net['source'].nunique()} TFs, {net['target'].nunique()} targets")

    # Prepare gene-level statistic as a 1-sample matrix
    gene_stat = cfrna_de["edge_pre_vs_flight_diff"].dropna()
    # Create a DataFrame with one "sample" (the group-level statistic)
    mat = pd.DataFrame({"pre_vs_flight": gene_stat}).T

    # Remove duplicate edges (CollecTRI may have repeated source-target pairs)
    n_before = len(net)
    net = net.drop_duplicates(subset=["source", "target"])
    n_after = len(net)
    if n_before != n_after:
        print(f"  Removed {n_before - n_after} duplicate edges ({n_before} → {n_after})")

    # Run ULM
    print("  Running ULM inference...")
    estimates, pvals = dc.run_ulm(mat=mat, net=net, verbose=False)

    # Extract results
    tf_results = {}
    for tf in estimates.columns:
        est = float(estimates.loc["pre_vs_flight", tf])
        pv = float(pvals.loc["pre_vs_flight", tf])
        tf_results[tf] = {
            "activity_score": round(est, 6),
            "p_value": round(pv, 6),
            "direction": "up_in_flight" if est > 0 else "down_in_flight",
        }

    return tf_results


def map_tf_names_mouse_to_human(mouse_tfs, orth_map):
    """Map mouse TF names to human (case conversion + ortholog fallback).
    Mouse: Stat1, Ahr → Human: STAT1, AHR
    """
    mapped = {}
    unmapped = []

    for mouse_tf in mouse_tfs:
        # Try direct uppercase (most TFs follow this pattern)
        human_tf = mouse_tf.upper()
        mapped[mouse_tf] = human_tf

    return mapped


def main():
    print("=" * 60)
    print("Phase E: TF Activity Conservation")
    print("=" * 60)

    # Step 1: Infer human TF activity from cfRNA
    print("\n1. Inferring human TF activity from cfRNA...")
    cfrna_de = load_cfrna_de()

    try:
        human_tfs = infer_human_tf_activity(cfrna_de)
        print(f"  Human TFs inferred: {len(human_tfs)}")
    except ImportError as e:
        print(f"  decoupler not available: {e}")
        print("  Using fallback: simple gene overlap scoring")
        human_tfs = _fallback_tf_activity(cfrna_de)
    except Exception as e:
        print(f"  Error: {e}")
        human_tfs = _fallback_tf_activity(cfrna_de)

    if not human_tfs:
        print("  ERROR: No human TF activity computed. Exiting.")
        return

    # Significant human TFs
    sig_human = {k: v for k, v in human_tfs.items() if v["p_value"] < 0.05}
    print(f"  Significant human TFs (p<0.05): {len(sig_human)}")

    # Step 2: Load mouse TF activity per tissue
    print("\n2. Loading mouse TF activity per tissue...")
    orth_map = load_ortholog_map()
    mouse_tf_all = {}

    for tissue in ALL_TISSUES:
        try:
            tf_data = load_tf_activity(tissue)
            tf_results = tf_data.get("tf_results", {})
            mouse_tf_all[tissue] = tf_results
            n_sig = tf_data.get("n_significant_fdr05", 0)
            print(f"  {tissue}: {len(tf_results)} TFs, {n_sig} significant")
        except FileNotFoundError:
            print(f"  {tissue}: SKIPPED (no TF data)")
        except Exception as e:
            print(f"  {tissue}: ERROR — {e}")

    # Step 3: Correlate mouse vs human TF activity per tissue
    print("\n3. Correlating mouse vs human TF activity...")
    correlation_results = {}

    human_tf_names = set(human_tfs.keys())

    for tissue, mouse_tfs in mouse_tf_all.items():
        # Map mouse TF names to human
        mouse_tf_names = set(mouse_tfs.keys())
        tf_name_map = map_tf_names_mouse_to_human(mouse_tf_names, orth_map)

        # Find shared TFs
        shared_tfs = []
        for mouse_name, human_name in tf_name_map.items():
            if human_name in human_tf_names and mouse_name in mouse_tfs:
                mouse_data = mouse_tfs[mouse_name]
                human_data = human_tfs[human_name]

                # Mouse: use t-statistic or mean_diff direction
                mouse_score = mouse_data.get("mean_flt", 0) - mouse_data.get("mean_gc", 0)
                human_score = human_data["activity_score"]

                shared_tfs.append({
                    "mouse_tf": mouse_name,
                    "human_tf": human_name,
                    "mouse_score": mouse_score,
                    "human_score": human_score,
                    "mouse_fdr": mouse_data.get("fdr_p", 1.0),
                    "human_p": human_data["p_value"],
                })

        n_shared = len(shared_tfs)
        if n_shared < 10:
            print(f"  {tissue:>15s}: only {n_shared} shared TFs, SKIPPING")
            continue

        mouse_scores = np.array([t["mouse_score"] for t in shared_tfs])
        human_scores = np.array([t["human_score"] for t in shared_tfs])

        # Spearman correlation
        rho, p_val = stats.spearmanr(mouse_scores, human_scores)

        # Bootstrap CI
        rng = np.random.default_rng(42)
        boot_rhos = []
        for _ in range(2000):
            idx = rng.integers(0, n_shared, size=n_shared)
            r, _ = stats.spearmanr(mouse_scores[idx], human_scores[idx])
            if not np.isnan(r):
                boot_rhos.append(r)
        ci_lower = np.percentile(boot_rhos, 2.5) if boot_rhos else np.nan
        ci_upper = np.percentile(boot_rhos, 97.5) if boot_rhos else np.nan

        # Direction concordance for jointly significant TFs
        sig_both = [t for t in shared_tfs
                    if t["mouse_fdr"] < 0.05 and t["human_p"] < 0.05]
        if sig_both:
            n_concordant = sum(1 for t in sig_both
                             if (t["mouse_score"] > 0) == (t["human_score"] > 0))
            concordance_rate = n_concordant / len(sig_both)
        else:
            n_concordant = 0
            concordance_rate = None

        # Top TFs by concordance (both significant, same direction)
        top_concordant = sorted(
            [t for t in sig_both if (t["mouse_score"] > 0) == (t["human_score"] > 0)],
            key=lambda x: abs(x["human_score"]), reverse=True)

        correlation_results[tissue] = {
            "n_shared_tfs": n_shared,
            "spearman_rho": round(float(rho), 4),
            "p_value": float(p_val),
            "ci_lower": round(float(ci_lower), 4),
            "ci_upper": round(float(ci_upper), 4),
            "significant": p_val < 0.05,
            "n_sig_both": len(sig_both),
            "n_concordant_sig": n_concordant,
            "concordance_rate": round(concordance_rate, 4) if concordance_rate is not None else None,
            "top_conserved_tfs": [
                {"mouse_tf": t["mouse_tf"], "human_tf": t["human_tf"],
                 "mouse_score": round(t["mouse_score"], 4),
                 "human_score": round(t["human_score"], 4)}
                for t in top_concordant[:10]
            ],
        }

        sig = "*" if p_val < 0.05 else ""
        conc_str = f"concordance={concordance_rate:.0%}" if concordance_rate is not None else "no sig overlap"
        print(f"  {tissue:>15s}: r={rho:+.3f} [{ci_lower:+.3f}, {ci_upper:+.3f}] "
              f"p={p_val:.4f}{sig} ({n_shared} TFs, {conc_str})")

    # Identify conserved TF regulators (significant in ≥2 mouse tissues + human)
    print("\n4. Identifying conserved TF regulators...")
    tf_tissue_count = {}
    for tissue, results in correlation_results.items():
        for tf_info in results.get("top_conserved_tfs", []):
            tf_name = tf_info["human_tf"]
            if tf_name not in tf_tissue_count:
                tf_tissue_count[tf_name] = {"tissues": [], "scores": []}
            tf_tissue_count[tf_name]["tissues"].append(tissue)
            tf_tissue_count[tf_name]["scores"].append(tf_info["mouse_score"])

    conserved_tfs = {k: v for k, v in tf_tissue_count.items() if len(v["tissues"]) >= 2}
    print(f"  Conserved TFs (sig in ≥2 mouse tissues + human): {len(conserved_tfs)}")
    for tf_name, info in sorted(conserved_tfs.items(),
                                 key=lambda x: len(x[1]["tissues"]), reverse=True)[:10]:
        tissues_str = ", ".join(info["tissues"])
        human_score = human_tfs.get(tf_name, {}).get("activity_score", 0)
        print(f"    {tf_name}: {len(info['tissues'])} tissues ({tissues_str}), "
              f"human={human_score:+.3f}")

    # Assemble output
    mean_rho = np.mean([v["spearman_rho"] for v in correlation_results.values()]) \
        if correlation_results else None

    output = {
        "n_human_tfs": len(human_tfs),
        "n_sig_human_tfs": len(sig_human),
        "n_tissues_analyzed": len(correlation_results),
        "per_tissue_correlations": correlation_results,
        "mean_rho": round(float(mean_rho), 4) if mean_rho is not None else None,
        "conserved_tfs": {k: {"n_tissues": len(v["tissues"]),
                              "tissues": v["tissues"],
                              "human_activity": round(
                                  human_tfs.get(k, {}).get("activity_score", 0), 4)}
                         for k, v in conserved_tfs.items()},
        "top_human_tfs": sorted(
            [{"tf": k, "activity": round(v["activity_score"], 4),
              "p_value": v["p_value"], "direction": v["direction"]}
             for k, v in human_tfs.items() if v["p_value"] < 0.05],
            key=lambda x: abs(x["activity"]), reverse=True)[:20],
        "timestamp": datetime.now().isoformat(),
    }

    save_json(output, V6_EVAL / "V6_E_tf_conservation.json")
    print("\nDone!")


def _fallback_tf_activity(cfrna_de):
    """Fallback TF activity using simple gene set scoring."""
    # Without decoupler, return empty
    print("  Warning: decoupler required for TF inference. Install with: pip install decoupler")
    return {}


if __name__ == "__main__":
    main()
