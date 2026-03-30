#!/usr/bin/env python3
"""
ppi_enrichment.py — PPI enrichment via STRING API v12 for SHAP top-100 genes.

For each tissue (all 8):
  1. Get top-100 SHAP genes (averaged across 2 methods)
  2. POST to STRING API → interaction network (edges, density, enrichment p)
  3. GET STRING functional enrichment (GO BP, KEGG, Reactome, p_fdr < 0.05)

Also queries WGCNA hub genes for 6 LOMO tissues.

Input:
  v4/evaluation/SHAP_{tissue}_{method}.json (16 files)
  v4/evaluation/WGCNA_hub_genes.json (optional, for hub gene PPI)

Output:
  v4/evaluation/PPI_enrichment.json

Usage:
  python ppi_enrichment.py [--tissues liver ...]
  python ppi_enrichment.py  # all 8 tissues

Notes:
  - STRING species 10090 = Mus musculus
  - Requires network access; rate-limit: 1 req/sec
  - required_score=400 (medium confidence)
"""

import json
import time
import argparse
import warnings
import numpy as np
import requests
from pathlib import Path
from datetime import datetime

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent.parent.parent
EVAL_DIR = BASE_DIR / "v4" / "evaluation"

# ── Constants ─────────────────────────────────────────────────────────────────
STRING_API_URL = "https://string-db.org/api/json"
STRING_SPECIES  = 10090   # Mus musculus
REQUIRED_SCORE  = 400     # medium confidence (0–1000 scale)
TOP_N_SHAP      = 100
ALPHA_FDR       = 0.05
REQUEST_DELAY   = 1.1     # seconds between requests (STRING rate limit)

ALL_TISSUES = ["liver", "gastrocnemius", "kidney", "thymus", "eye", "skin", "lung", "colon"]

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


def load_shap_top100(tissue):
    """Load top-100 genes averaged across 2 methods for a tissue.

    Returns: list of gene symbols (mouse, mixed case e.g. "Npas2")
    """
    methods = TISSUE_TOP_METHODS.get(tissue, [])
    gene_scores = {}
    gene_counts  = {}

    for method in methods:
        shap_path = EVAL_DIR / f"SHAP_{tissue}_{method}.json"
        if not shap_path.exists():
            continue
        with open(shap_path) as f:
            data = json.load(f)
        gene_imp = data.get("gene_importance", {})
        for gene, imp in gene_imp.items():
            gene_scores[gene] = gene_scores.get(gene, 0.0) + abs(float(imp))
            gene_counts[gene] = gene_counts.get(gene, 0) + 1

    if not gene_scores:
        return []

    averaged = {g: gene_scores[g] / gene_counts[g] for g in gene_scores}
    sorted_genes = sorted(averaged.keys(), key=lambda g: averaged[g], reverse=True)
    return sorted_genes[:TOP_N_SHAP]


def string_network(genes, label="query"):
    """Query STRING API for network interactions among a gene list.

    Returns dict with: edges, density, ppi_enrichment_pvalue, interaction_partners
    """
    if not genes:
        return {"error": "empty gene list"}

    url = f"{STRING_API_URL}/network"
    params = {
        "identifiers": "%0d".join(genes),
        "species": STRING_SPECIES,
        "required_score": REQUIRED_SCORE,
        "caller_identity": "GeneLabBench_v4",
    }

    try:
        resp = requests.post(url, data=params, timeout=30)
        resp.raise_for_status()
        interactions = resp.json()
    except Exception as e:
        return {"error": str(e)}

    if not isinstance(interactions, list) or len(interactions) == 0:
        return {
            "n_genes_queried": len(genes),
            "n_interactions": 0,
            "density": 0.0,
            "ppi_enrichment_pvalue": 1.0,
        }

    # Parse interactions
    n_edges = len(interactions)
    # Collect unique nodes that have interactions
    interacting_genes = set()
    for row in interactions:
        interacting_genes.add(row.get("preferredName_A", ""))
        interacting_genes.add(row.get("preferredName_B", ""))
    interacting_genes.discard("")

    n_nodes = len(interacting_genes)
    max_edges = n_nodes * (n_nodes - 1) / 2 if n_nodes > 1 else 1
    density = round(n_edges / max_edges, 4) if max_edges > 0 else 0.0

    # Look for STRING's own PPI enrichment p-value in results
    ppi_pvalue = None
    for row in interactions[:1]:
        if "enrichment" in row:
            try:
                ppi_pvalue = float(row["enrichment"])
            except:
                pass

    result = {
        "n_genes_queried": len(genes),
        "n_genes_with_interactions": n_nodes,
        "n_interactions": n_edges,
        "density": density,
        "ppi_enrichment_pvalue": ppi_pvalue,
        "top_interactors": sorted(list(interacting_genes))[:20],
    }
    return result


def string_enrichment(genes, label="query"):
    """Query STRING API for functional enrichment of a gene set.

    Returns list of significant terms (p_fdr < ALPHA_FDR) across
    GO Biological Process, KEGG, and Reactome.
    """
    if not genes:
        return []

    url = f"{STRING_API_URL}/enrichment"
    params = {
        "identifiers": "%0d".join(genes),
        "species": STRING_SPECIES,
        "caller_identity": "GeneLabBench_v4",
    }

    try:
        resp = requests.post(url, data=params, timeout=60)
        resp.raise_for_status()
        terms = resp.json()
    except Exception as e:
        return [{"error": str(e)}]

    if not isinstance(terms, list):
        return []

    # Filter to relevant categories and significant terms
    keep_categories = {"Process", "KEGG", "Reactome", "Component", "Function"}
    # STRING uses: "Process" = GO BP, "Component" = GO CC, "Function" = GO MF
    # Prioritize: Process (GO BP), KEGG, Reactome

    sig_terms = []
    for term in terms:
        cat  = term.get("category", "")
        fdr  = term.get("fdr", 1.0)
        desc = term.get("description", "")
        if cat in keep_categories and fdr < ALPHA_FDR:
            sig_terms.append({
                "category": cat,
                "term_id": term.get("term", ""),
                "description": desc,
                "p_fdr": round(float(fdr), 6),
                "n_genes_in_term": term.get("number_of_genes_in_background", None),
                "n_overlap": term.get("observed", None),
                "genes_in_term": (term.get("inputGenes") or [])[:10]
                    if isinstance(term.get("inputGenes"), list)
                    else str(term.get("inputGenes", "")).split(",")[:10],
            })

    sig_terms.sort(key=lambda x: x["p_fdr"])
    return sig_terms


def process_tissue(tissue, hub_genes_data=None):
    """Run STRING queries for one tissue's SHAP top-100 genes."""
    genes = load_shap_top100(tissue)
    if not genes:
        print(f"  [{tissue}] No SHAP genes found — skip")
        return None

    print(f"  [{tissue}] {len(genes)} SHAP top genes → STRING query...")

    # Network query
    network = string_network(genes, label=f"{tissue}_shap")
    time.sleep(REQUEST_DELAY)

    # Enrichment query
    enrichment = string_enrichment(genes, label=f"{tissue}_shap")
    time.sleep(REQUEST_DELAY)

    n_enrich = len([t for t in enrichment if "error" not in t])
    print(f"    Network: {network.get('n_interactions', 0)} interactions, "
          f"density={network.get('density', 0):.4f}")
    print(f"    Enrichment: {n_enrich} significant terms (p_fdr<{ALPHA_FDR})")

    result = {
        "tissue": tissue,
        "shap_genes_queried": genes,
        "n_shap_genes": len(genes),
        "network": network,
        "enrichment": enrichment[:50],  # cap at top-50 terms
        "n_significant_terms": n_enrich,
    }

    # Hub genes query (LOMO tissues only, if WGCNA available)
    if hub_genes_data and tissue in hub_genes_data.get("tissues", {}):
        tissue_hubs = hub_genes_data["tissues"][tissue]
        all_hub_genes = []
        for mod_data in tissue_hubs.get("modules", {}).values():
            all_hub_genes.extend(mod_data.get("hub_genes", []))
        all_hub_genes = list(dict.fromkeys(all_hub_genes))[:100]  # deduplicate, cap at 100

        if all_hub_genes:
            print(f"    Hub genes ({len(all_hub_genes)}) → STRING query...")
            hub_network = string_network(all_hub_genes, label=f"{tissue}_hub")
            time.sleep(REQUEST_DELAY)
            hub_enrichment = string_enrichment(all_hub_genes, label=f"{tissue}_hub")
            time.sleep(REQUEST_DELAY)

            result["hub_genes_network"] = hub_network
            result["hub_genes_enrichment"] = hub_enrichment[:30]
            print(f"    Hub network: {hub_network.get('n_interactions', 0)} interactions")

    return result


def main():
    parser = argparse.ArgumentParser(description="PPI enrichment via STRING API for SHAP top genes")
    parser.add_argument("--tissues", nargs="+", default=ALL_TISSUES,
                        help="Tissues to query (default: all 8)")
    parser.add_argument("--no-hub", action="store_true",
                        help="Skip hub gene STRING queries")
    args = parser.parse_args()

    # Load WGCNA hub genes if available
    hub_genes_data = None
    hub_path = EVAL_DIR / "WGCNA_hub_genes.json"
    if not args.no_hub and hub_path.exists():
        with open(hub_path) as f:
            hub_genes_data = json.load(f)
        print(f"Loaded hub genes for {len(hub_genes_data.get('tissues', {}))} tissues")
    else:
        print("WGCNA hub genes not available — skipping hub gene PPI queries")

    print(f"\nQuerying STRING API (species={STRING_SPECIES}, score≥{REQUIRED_SCORE})...")
    print("Note: ~1 sec delay between requests to respect rate limits\n")

    all_results = {}
    for tissue in args.tissues:
        result = process_tissue(tissue, hub_genes_data)
        if result is not None:
            all_results[tissue] = result

    # Save
    output = {
        "tissues_processed": list(all_results.keys()),
        "string_species": STRING_SPECIES,
        "required_score": REQUIRED_SCORE,
        "alpha_fdr": ALPHA_FDR,
        "top_n_shap": TOP_N_SHAP,
        "tissues": all_results,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
    }

    out_path = EVAL_DIR / "PPI_enrichment.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nSaved: {out_path}")

    # Summary
    print("\n═══ Summary ═══")
    print(f"  {'Tissue':<15} {'Interactions':<14} {'Density':<10} {'Enrich terms'}")
    for tissue, res in all_results.items():
        net = res.get("network", {})
        print(f"  {tissue:<15} {str(net.get('n_interactions', '?')):<14} "
              f"{str(net.get('density', '?')):<10} "
              f"{res.get('n_significant_terms', 0)}")

    # Tissues with significant enrichment
    sig = [t for t, r in all_results.items() if r.get("n_significant_terms", 0) > 0]
    print(f"\nTissues with ≥1 significant STRING term: {len(sig)}/8 — {', '.join(sig)}")


if __name__ == "__main__":
    main()
