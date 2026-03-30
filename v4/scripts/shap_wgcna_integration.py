#!/usr/bin/env python3
"""
shap_wgcna_integration.py — Map SHAP top genes to WGCNA modules.

For each tissue:
  1. Load SHAP top-100 genes (averaged across 2 methods per tissue)
  2. Find which WGCNA module each SHAP gene belongs to
  3. Hypergeometric enrichment test: are SHAP genes concentrated in any module?
  4. Map consensus genes (≥3 tissues, ≥2 methods) to their modules

Input:
  v4/evaluation/SHAP_{tissue}_{method}.json    (16 files)
  v4/evaluation/SHAP_consensus.json
  v4/evaluation/WGCNA_{tissue}_modules.json    (6 files, from wgcna_hub_genes.py)

Output:
  v4/evaluation/SHAP_WGCNA_integration.json

Usage:
  python shap_wgcna_integration.py
"""

import json
import numpy as np
from pathlib import Path
from datetime import datetime
from scipy.stats import hypergeom

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR  = Path(__file__).resolve().parent.parent.parent
EVAL_DIR  = BASE_DIR / "v4" / "evaluation"

LOMO_TISSUES = ["liver", "gastrocnemius", "kidney", "thymus", "eye", "skin"]

# Top-2 methods per tissue (by Phase 1 gene AUROC)
TISSUE_TOP_METHODS = {
    "liver":         ["svm_rbf", "mlp"],
    "gastrocnemius": ["elasticnet_lr", "pca_lr"],
    "kidney":        ["xgb", "rf"],
    "thymus":        ["knn", "pca_lr"],
    "eye":           ["pca_lr", "xgb"],
    "skin":          ["elasticnet_lr", "pca_lr"],
    "lung":          ["elasticnet_lr", "pca_lr"],
    "colon":         ["pca_lr", "xgb"],
}

TOP_N_SHAP = 100
ALPHA      = 0.05


def load_shap_genes(tissue, method, top_n=TOP_N_SHAP):
    """Load top-N SHAP gene importances for a tissue/method pair.

    Returns list of (gene, importance) sorted descending.
    """
    shap_path = EVAL_DIR / f"SHAP_{tissue}_{method}.json"
    if not shap_path.exists():
        return []
    with open(shap_path) as f:
        data = json.load(f)

    gene_imp = data.get("gene_importance", {})
    if not gene_imp:
        return []

    sorted_genes = sorted(gene_imp.items(), key=lambda x: abs(float(x[1])), reverse=True)
    return sorted_genes[:top_n]


def get_top100_averaged(tissue):
    """Get top-100 genes averaged across both methods for a tissue.

    Returns: (top_genes_list, top_importances_array, all_gene_importances_dict)
    """
    methods = TISSUE_TOP_METHODS.get(tissue, [])
    if not methods:
        return [], np.array([]), {}

    # Collect per-method gene importances
    gene_scores = {}  # gene → sum of abs importance
    gene_counts  = {}  # gene → count

    all_genes_set = set()
    for method in methods:
        ranked = load_shap_genes(tissue, method, top_n=500)  # use top-500 for averaging
        if not ranked:
            continue
        for gene, imp in ranked:
            all_genes_set.add(gene)
            gene_scores[gene] = gene_scores.get(gene, 0.0) + abs(float(imp))
            gene_counts[gene] = gene_counts.get(gene, 0) + 1

    if not gene_scores:
        return [], np.array([]), {}

    # Average importance
    averaged = {g: gene_scores[g] / gene_counts[g] for g in gene_scores}
    sorted_genes = sorted(averaged.items(), key=lambda x: x[1], reverse=True)

    top_genes  = [g for g, _ in sorted_genes[:TOP_N_SHAP]]
    top_imps   = np.array([v for _, v in sorted_genes[:TOP_N_SHAP]])

    return top_genes, top_imps, averaged


def hypergeometric_enrichment(shap_genes, module_genes, all_wgcna_genes, top_n):
    """Test if SHAP top genes are enriched in a WGCNA module.

    H0: SHAP top genes distributed uniformly across WGCNA gene universe.
    N = total WGCNA genes, K = module size, n = top_n SHAP genes, k = overlap
    """
    N = len(all_wgcna_genes)
    K = len(module_genes)
    n = min(top_n, len(shap_genes))

    shap_set  = set(g.upper() for g in shap_genes)
    mod_set   = set(g.upper() for g in module_genes)
    overlap   = shap_set & mod_set
    k         = len(overlap)

    if N == 0 or K == 0 or n == 0:
        return k, 1.0, list(overlap)

    # P(X >= k) under hypergeometric
    p_val = hypergeom.sf(k - 1, N, K, n) if k > 0 else 1.0

    return k, float(p_val), list(overlap)


def load_wgcna_modules(tissue):
    """Load WGCNA module assignments from per-tissue JSON.

    Returns: {module_color: [gene1, gene2, ...], "all_genes": set()}
    or None if file missing.
    """
    wgcna_path = EVAL_DIR / f"WGCNA_{tissue}_modules.json"
    if not wgcna_path.exists():
        print(f"  [{tissue}] WGCNA modules file not found: {wgcna_path.name}")
        return None

    with open(wgcna_path) as f:
        data = json.load(f)

    modules_data = data.get("modules", {})

    # Rebuild gene→module map from hub + top20 genes (full gene list not stored in JSON)
    # We need to also load the original CSV for full membership
    mod_csv_path = BASE_DIR / "v4" / "wgcna_outputs" / tissue / "module_assignments.csv"
    if not mod_csv_path.exists():
        print(f"  [{tissue}] module_assignments.csv not found — using top20 genes only")
        # Fall back to top20 from JSON
        gene_to_module = {}
        all_genes = set()
        for mod_color, mod_data in modules_data.items():
            for gene in mod_data.get("top20_genes", []):
                gene_to_module[gene] = mod_color
                all_genes.add(gene)
        return {"gene_to_module": gene_to_module, "all_genes": all_genes, "modules": modules_data}

    import pandas as pd
    mod_df = pd.read_csv(mod_csv_path)
    gene_to_module = dict(zip(mod_df["gene"].astype(str), mod_df["module_color"].astype(str)))
    all_genes = set(mod_df["gene"].astype(str))

    # Build module→genes dict
    module_to_genes = {}
    for gene, mod in gene_to_module.items():
        module_to_genes.setdefault(mod, []).append(gene)

    return {
        "gene_to_module": gene_to_module,
        "module_to_genes": module_to_genes,
        "all_genes": all_genes,
        "modules": modules_data,
    }


def process_tissue(tissue, consensus_genes):
    """Integrate SHAP and WGCNA for one tissue."""
    print(f"\n[{tissue}] SHAP-WGCNA integration...")

    # Load SHAP top genes
    shap_top, shap_imps, all_shap = get_top100_averaged(tissue)
    if not shap_top:
        print(f"  No SHAP data found")
        return None

    # Load WGCNA modules
    wgcna = load_wgcna_modules(tissue)
    if wgcna is None:
        print(f"  No WGCNA data found")
        return None

    gene_to_module = wgcna["gene_to_module"]
    module_to_genes = wgcna.get("module_to_genes", {})
    all_wgcna_genes = wgcna["all_genes"]

    print(f"  SHAP top-100 genes, WGCNA: {len(all_wgcna_genes)} genes, "
          f"{len(module_to_genes)} modules")

    # Map SHAP top genes to modules
    shap_to_module = {}
    for gene in shap_top:
        mod = gene_to_module.get(gene)
        if mod is None:
            # Try uppercase / case-insensitive
            for wgcna_gene in all_wgcna_genes:
                if wgcna_gene.upper() == gene.upper():
                    mod = gene_to_module.get(wgcna_gene)
                    break
        shap_to_module[gene] = mod  # None if not in WGCNA

    n_mapped = sum(1 for m in shap_to_module.values() if m is not None)
    print(f"  SHAP→module mapping: {n_mapped}/{len(shap_top)} genes mapped")

    # Module counts
    from collections import Counter
    mod_counts = Counter(m for m in shap_to_module.values() if m is not None and m != "grey")

    # Hypergeometric enrichment per module
    enriched_modules = []
    for mod_color, mod_genes in module_to_genes.items():
        if mod_color == "grey":
            continue
        k, p_val, overlap = hypergeometric_enrichment(
            shap_top, mod_genes, all_wgcna_genes, TOP_N_SHAP
        )
        if p_val < ALPHA:
            enriched_modules.append({
                "module": mod_color,
                "n_shap_in_module": k,
                "n_module_genes": len(mod_genes),
                "overlap_genes": sorted(overlap)[:20],
                "p_value": round(p_val, 6),
                "enrichment_ratio": round(k / (len(mod_genes) / len(all_wgcna_genes) * TOP_N_SHAP), 3)
                    if len(mod_genes) > 0 else 0,
            })
    enriched_modules.sort(key=lambda x: x["p_value"])

    # Consensus gene module assignments
    consensus_assignments = {}
    for cg in consensus_genes:
        gene = cg.get("gene", "")
        mod = gene_to_module.get(gene)
        if mod is None:
            for wgcna_gene in all_wgcna_genes:
                if wgcna_gene.upper() == gene.upper():
                    mod = gene_to_module.get(wgcna_gene)
                    break
        if tissue in cg.get("tissues", []):
            consensus_assignments[gene] = mod

    print(f"  Enriched modules: {len(enriched_modules)} (p<{ALPHA})")
    print(f"  Consensus genes mapped: {sum(1 for m in consensus_assignments.values() if m)}"
          f"/{len(consensus_assignments)}")

    return {
        "tissue": tissue,
        "n_shap_top100": len(shap_top),
        "n_shap_in_wgcna": n_mapped,
        "shap_module_distribution": dict(mod_counts.most_common(10)),
        "enriched_modules": enriched_modules,
        "shap_module_map": shap_to_module,
        "consensus_gene_modules": consensus_assignments,
    }


def main():
    # Load consensus
    consensus_path = EVAL_DIR / "SHAP_consensus.json"
    consensus_genes = []
    if consensus_path.exists():
        with open(consensus_path) as f:
            cons_data = json.load(f)
        consensus_genes = cons_data.get("consensus_genes", [])
        print(f"Loaded {len(consensus_genes)} consensus genes")
    else:
        print("WARNING: SHAP_consensus.json not found")

    # Process all LOMO tissues
    all_results = {}
    for tissue in LOMO_TISSUES:
        result = process_tissue(tissue, consensus_genes)
        if result is not None:
            all_results[tissue] = result

    # Cross-tissue consensus module summary
    consensus_module_summary = {}
    for cg in consensus_genes:
        gene = cg["gene"]
        modules_across_tissues = {}
        for tissue, res in all_results.items():
            mod = res["consensus_gene_modules"].get(gene)
            if mod and mod != "grey":
                modules_across_tissues[tissue] = mod
        if modules_across_tissues:
            consensus_module_summary[gene] = {
                "n_tissues_with_module": len(modules_across_tissues),
                "tissue_modules": modules_across_tissues,
                "n_tissues_consensus": cg.get("n_tissues"),
                "n_methods_consensus": cg.get("n_methods"),
            }

    # Save
    output = {
        "tissues_processed": list(all_results.keys()),
        "top_n_shap": TOP_N_SHAP,
        "alpha": ALPHA,
        "tissues": all_results,
        "consensus_gene_modules": consensus_module_summary,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
    }

    out_path = EVAL_DIR / "SHAP_WGCNA_integration.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nSaved: {out_path}")

    # Summary
    print("\n═══ Summary ═══")
    for tissue, res in all_results.items():
        print(f"  {tissue:<15} mapped={res['n_shap_in_wgcna']}/{res['n_shap_top100']} "
              f"enriched_modules={len(res['enriched_modules'])}")

    if consensus_module_summary:
        print("\nConsensus gene module assignments:")
        for gene, info in consensus_module_summary.items():
            mods = ", ".join(f"{t}:{m}" for t, m in info["tissue_modules"].items())
            print(f"  {gene:<12} ({info['n_tissues_with_module']} tissues): {mods}")


if __name__ == "__main__":
    main()
