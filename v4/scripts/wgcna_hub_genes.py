#!/usr/bin/env python3
"""
wgcna_hub_genes.py — Extract hub genes from WGCNA module assignments.

Hub genes defined as: kME > 0.8 in assigned module AND in top-10 by kME per module.
Also computes intramodular hub score from kME values.

Input:  v4/wgcna_outputs/{tissue}/
          module_assignments.csv  (gene, module_color, kME_assigned, kME{color}...)
          eigengenes.csv          (samples × modules)
          trait_correlations_r.csv
          trait_correlations_p.csv
          summary.json

Output: v4/evaluation/WGCNA_hub_genes.json
        v4/evaluation/WGCNA_{tissue}_modules.json  (per-tissue, for downstream scripts)

Usage:
  python wgcna_hub_genes.py [--tissues liver kidney ...]
"""

import json
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR   = Path(__file__).resolve().parent.parent.parent
WGCNA_DIR  = BASE_DIR / "v4" / "wgcna_outputs"
EVAL_DIR   = BASE_DIR / "v4" / "evaluation"
EVAL_DIR.mkdir(parents=True, exist_ok=True)

LOMO_TISSUES = ["liver", "gastrocnemius", "kidney", "thymus", "eye", "skin"]

# Hub gene thresholds
KME_THRESHOLD   = 0.8   # minimum kME for hub classification
TOP_N_HUB       = 10    # top N per module


def load_tissue_wgcna(tissue):
    """Load WGCNA outputs for one tissue. Returns dict or None if missing."""
    tissue_dir = WGCNA_DIR / tissue

    mod_path    = tissue_dir / "module_assignments.csv"
    eigen_path  = tissue_dir / "eigengenes.csv"
    summary_path = tissue_dir / "summary.json"
    trait_r_path = tissue_dir / "trait_correlations_r.csv"
    trait_p_path = tissue_dir / "trait_correlations_p.csv"

    if not mod_path.exists():
        print(f"  [{tissue}] module_assignments.csv not found — skip")
        return None

    mod_df = pd.read_csv(mod_path)

    eigengenes = None
    if eigen_path.exists():
        eigengenes = pd.read_csv(eigen_path, index_col=0)

    summary = {}
    if summary_path.exists():
        with open(summary_path) as f:
            summary = json.load(f)

    trait_r = trait_p = None
    if trait_r_path.exists() and trait_p_path.exists():
        trait_r = pd.read_csv(trait_r_path)
        trait_p = pd.read_csv(trait_p_path)

    return {
        "mod_df": mod_df,
        "eigengenes": eigengenes,
        "summary": summary,
        "trait_r": trait_r,
        "trait_p": trait_p,
    }


def extract_hub_genes(mod_df, tissue):
    """Extract hub genes per module.

    Returns dict: {module_color: {hub_genes, kME_values, n_genes_total, top20_genes}}
    """
    modules = [m for m in mod_df["module_color"].unique() if m != "grey"]
    module_results = {}

    # Find all kME columns (pattern: kME{color})
    kme_cols = [c for c in mod_df.columns if c.startswith("kME") and c != "kME_assigned"]

    for mod_color in modules:
        mod_genes = mod_df[mod_df["module_color"] == mod_color]
        n_genes = len(mod_genes)

        # Identify the kME column for this module
        kme_col = f"kME{mod_color}"

        if kme_col not in mod_df.columns:
            # Fall back to kME_assigned if specific column missing
            if "kME_assigned" in mod_df.columns:
                kme_vals = mod_genes["kME_assigned"].fillna(0.0)
            else:
                module_results[mod_color] = {
                    "n_genes": n_genes,
                    "hub_genes": [],
                    "kME_values": {},
                    "top20_genes": [],
                    "error": f"kME column '{kme_col}' not found"
                }
                continue
        else:
            kme_vals = mod_genes[kme_col].fillna(0.0)

        # Sort by kME descending
        sorted_idx = kme_vals.sort_values(ascending=False).index
        sorted_genes = mod_genes.loc[sorted_idx, "gene"].values
        sorted_kme   = kme_vals.loc[sorted_idx].values

        # Top-10 hub genes
        top_n = min(TOP_N_HUB, len(sorted_genes))
        top_genes = sorted_genes[:top_n]
        top_kme   = sorted_kme[:top_n]

        # Hub criterion: kME > threshold AND in top-10
        hub_mask  = top_kme > KME_THRESHOLD
        hub_genes = top_genes[hub_mask].tolist()
        hub_kme   = top_kme[hub_mask].tolist()

        # Top-20 for broader reporting
        top20_n     = min(20, len(sorted_genes))
        top20_genes = sorted_genes[:top20_n].tolist()
        top20_kme   = sorted_kme[:top20_n].tolist()

        module_results[mod_color] = {
            "n_genes": n_genes,
            "hub_genes": hub_genes,
            "kME_values": {g: round(float(k), 4) for g, k in zip(hub_genes, hub_kme)},
            "top20_genes": top20_genes,
            "top20_kME": [round(float(k), 4) for k in top20_kme],
            "n_hub_genes": len(hub_genes),
            "mean_kME_top10": round(float(np.mean(top_kme)), 4) if len(top_kme) > 0 else None,
        }

    return module_results


def extract_trait_correlations(trait_r, trait_p):
    """Extract significant module-trait associations."""
    if trait_r is None or trait_p is None:
        return {}

    results = {}
    trait_cols = [c for c in trait_r.columns if c != "module"]

    for _, row_r in trait_r.iterrows():
        mod = row_r.get("module", str(row_r.name))
        if str(mod).startswith("ME"):
            mod = mod[2:]  # strip "ME" prefix

        row_p = trait_p[trait_p["module"] == row_r["module"]].squeeze() if "module" in trait_p.columns else None

        sig_traits = []
        for trait in trait_cols:
            r_val = float(row_r[trait]) if not pd.isna(row_r.get(trait, np.nan)) else None
            p_val = None
            if row_p is not None and not row_p.empty and trait in row_p.index:
                p_val = float(row_p[trait]) if not pd.isna(row_p[trait]) else None

            if r_val is not None and p_val is not None and p_val < 0.05:
                sig_traits.append({
                    "trait": trait,
                    "r": round(r_val, 4),
                    "p": round(p_val, 6),
                })

        results[mod] = sorted(sig_traits, key=lambda x: x["p"])

    return results


def process_tissue(tissue):
    """Process one tissue → returns tissue-level WGCNA summary dict."""
    data = load_tissue_wgcna(tissue)
    if data is None:
        return None

    mod_df    = data["mod_df"]
    summary   = data["summary"]
    eigengenes = data["eigengenes"]

    print(f"\n[{tissue}] Processing WGCNA outputs...")
    print(f"  Total genes: {len(mod_df)}, modules: {summary.get('n_modules_final', '?')}")

    # Extract hub genes
    hub_genes = extract_hub_genes(mod_df, tissue)

    # Trait correlations
    trait_cors = extract_trait_correlations(data["trait_r"], data["trait_p"])

    # Module-level summary
    modules_summary = {}
    for mod_color, hub_data in hub_genes.items():
        modules_summary[mod_color] = {
            **hub_data,
            "significant_trait_associations": trait_cors.get(mod_color, []),
        }

    # Eigengene statistics (if available)
    eigengene_stats = {}
    if eigengenes is not None:
        for col in eigengenes.columns:
            mod_col = col[2:] if col.startswith("ME") else col
            eigengene_stats[mod_col] = {
                "mean": round(float(eigengenes[col].mean()), 4),
                "std":  round(float(eigengenes[col].std()), 4),
            }

    # Most significant module-condition association (for reporting)
    condition_hits = []
    if data["trait_r"] is not None:
        for mod_color, sig_list in trait_cors.items():
            for item in sig_list:
                if item["trait"] == "condition":
                    condition_hits.append({
                        "module": mod_color,
                        "r": item["r"],
                        "p": item["p"],
                        "direction": "FLT_higher" if item["r"] > 0 else "GC_higher",
                    })
    condition_hits.sort(key=lambda x: x["p"])

    tissue_result = {
        "tissue": tissue,
        "n_samples": summary.get("n_samples"),
        "n_genes_input": summary.get("n_genes_input"),
        "soft_threshold_beta": summary.get("soft_threshold_beta"),
        "r2_achieved": summary.get("r2_achieved"),
        "r2_threshold_met": summary.get("r2_threshold_met", False),
        "n_modules": summary.get("n_modules_final"),
        "n_grey_genes": summary.get("n_grey_genes", 0),
        "modules": modules_summary,
        "condition_modules": condition_hits[:5],  # top-5 most significant
        "eigengene_stats": eigengene_stats,
        "network_type": "signed hybrid",
        "min_module_size": summary.get("min_module_size"),
    }

    # Save per-tissue module JSON for downstream use
    tissue_json_path = EVAL_DIR / f"WGCNA_{tissue}_modules.json"
    with open(tissue_json_path, "w") as f:
        json.dump(tissue_result, f, indent=2, default=str)
    print(f"  Wrote: {tissue_json_path.name}")

    # Report
    n_hub_total = sum(d.get("n_hub_genes", 0) for d in hub_genes.values())
    print(f"  Hub genes: {n_hub_total} across {len(hub_genes)} modules")
    print(f"  Condition-associated modules: {len(condition_hits)}")
    if condition_hits:
        best = condition_hits[0]
        print(f"  Best: {best['module']} (r={best['r']:.3f}, p={best['p']:.4f}, {best['direction']})")

    return tissue_result


def main():
    parser = argparse.ArgumentParser(description="Extract hub genes from WGCNA outputs")
    parser.add_argument("--tissues", nargs="+", default=LOMO_TISSUES,
                        help="Tissues to process (default: all 6 LOMO tissues)")
    args = parser.parse_args()

    all_results = {}
    summaries   = []

    for tissue in args.tissues:
        result = process_tissue(tissue)
        if result is not None:
            all_results[tissue] = result
            summaries.append({
                "tissue": tissue,
                "n_modules": result.get("n_modules"),
                "beta": result.get("soft_threshold_beta"),
                "r2_threshold_met": result.get("r2_threshold_met"),
                "n_condition_modules": len(result.get("condition_modules", [])),
            })

    # Save aggregate hub genes JSON
    hub_genes_path = EVAL_DIR / "WGCNA_hub_genes.json"
    output = {
        "tissues_processed": list(all_results.keys()),
        "kme_threshold": KME_THRESHOLD,
        "top_n_per_module": TOP_N_HUB,
        "tissues": all_results,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
    }
    with open(hub_genes_path, "w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nSaved: {hub_genes_path}")

    # Print summary table
    print("\n═══ Summary ═══")
    print(f"  {'Tissue':<15} {'Modules':<9} {'Beta':<6} {'R² OK':<7} {'Condition hits'}")
    for s in summaries:
        print(f"  {s['tissue']:<15} {str(s['n_modules']):<9} {str(s['beta']):<6} "
              f"{'Y' if s['r2_threshold_met'] else 'N':<7} {s['n_condition_modules']}")


if __name__ == "__main__":
    main()
