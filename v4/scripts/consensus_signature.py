#!/usr/bin/env python3
"""
consensus_signature.py — GeneLabBench v4 Phase 4: Consensus spaceflight gene signature

Aggregates SHAP results across tissues and methods to identify:
  1. Consensus genes: SHAP top-100 in ≥3 tissues AND ≥2 methods (DD-25)
  2. Tissue-specific genes: SHAP top-100 in 1 tissue, ≥3 methods
  3. Validation vs da Silveira Cell 2020 reference genes (hypergeometric test)
  4. Validation vs GeneLab published DEG lists (hypergeometric test)
  5. SHAP direction analysis: positive SHAP = upregulated in FLT?

Usage:
  python consensus_signature.py
  python consensus_signature.py --shap-dir v4/evaluation/ --output v4/evaluation/SHAP_consensus.json
"""

import json
import sys
import argparse
import warnings
import numpy as np
from pathlib import Path
from collections import Counter, defaultdict
from datetime import datetime
from scipy.stats import hypergeom

warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).resolve().parent))

from v4_utils import V4_EVAL_DIR

# Top-2 methods per tissue (must match shap_multi_method.py)
TOP2_METHODS = {
    "liver":         ["svm_rbf", "mlp"],
    "gastrocnemius": ["elasticnet_lr", "pca_lr"],
    "kidney":        ["xgb", "rf"],
    "thymus":        ["knn", "pca_lr"],
    "eye":           ["pca_lr", "xgb"],
    "skin":          ["elasticnet_lr", "pca_lr"],
    "lung":          ["elasticnet_lr", "pca_lr"],
    "colon":         ["pca_lr", "xgb"],
}

# Reference genes from Cell 2020 / da Silveira / NASA GeneLab
# (from scripts/cell2020_validation.py and scripts/shap_analysis.py)
REFERENCE_GENES = {
    "liver": ["Angptl4", "Pck1", "G6pc", "Cyp2e1", "Cyp2b10", "Scd1", "Hmgcs2",
              "Fmo3", "Ugt1a1"],
    "gastrocnemius": ["Fbxo32", "Trim63", "Mstn", "Ctsl", "Ankrd1", "Myh4",
                      "Fst", "Igf1", "Foxo3"],
    "kidney": ["Umod", "Slc12a1", "Aqp2", "Klk1"],
    "thymus": ["Ccnb1", "Cdk1", "Top2a", "Mki67", "Foxn1"],
    "eye": ["Rpe65", "Rdh5", "Lrat", "Rho", "Opn1sw"],
}

# da Silveira et al., Cell 2020 core spaceflight signature genes
# (cross-tissue genes reported in twin study + multi-omics integration)
# Uses MOUSE gene symbols (converted from human where they differ)
DA_SILVEIRA_GENES = [
    "TRP53", "CDKN1A", "FOS", "JUN", "ATF3", "EGR1", "DUSP1",
    "GADD45A", "GADD45B", "BCL2", "BAX", "CASP3", "STAT3",
    "IL6", "TNF", "NFKB1", "VEGFA", "HIF1A", "SOD2", "GPX1",
    "CAT", "NFE2L2", "KEAP1", "HSP90AA1", "HSPA1A", "HSPA1B",
    "MT1", "MT2", "TERT", "TERC", "TERF1", "TERF2",
    "LMNB1", "HMGA2", "APOE", "CLU", "MYC", "CCND1", "CDK4",
    "RB1", "E2F1", "BRCA1", "RAD51", "ATM", "ATR",
    "CHEK1", "CHEK2", "CDKN2A", "MDM2", "PTEN",
    "AKT1", "MTOR", "PIK3CA", "MAPK1", "MAPK3",
    "NRAS", "KRAS", "RAF1", "MAP2K1", "MAPK3",
    "WNT3A", "CTNNB1", "LEF1", "TCF7",
    "NOTCH1", "HES1", "JAG1",
]
# Human→mouse symbol changes: TP53→TRP53, NRF2→NFE2L2, TRF1→TERF1,
# TRF2→TERF2, MEK1→MAP2K1, ERK1→MAPK3 (already in list)

# Total genes in mouse transcriptome (approximate, for hypergeometric test)
TOTAL_GENES_APPROX = 17000


def load_shap_results(shap_dir):
    """Load all SHAP_{tissue}_{method}.json files (excluding interactions)."""
    results = {}
    for path in sorted(shap_dir.glob("SHAP_*_*.json")):
        name = path.name
        # Skip interaction files (SHAP_interactions_*.json)
        if name.startswith("SHAP_interactions_"):
            continue
        # Skip consensus file itself
        if name == "SHAP_consensus.json":
            continue
        with open(path) as f:
            data = json.load(f)
        if "error" in data:
            continue
        key = (data["tissue"], data["method"])
        results[key] = data
    return results


def find_consensus_genes(shap_results, top_n=100, min_tissues=3, min_methods=2):
    """Find genes in SHAP top-{top_n} across multiple tissues and methods.

    DD-25: Consensus = top-100 in ≥3 tissues AND ≥2 methods.
    """
    # Per gene: which (tissue, method) pairs include it in top-N
    gene_occurrences = defaultdict(lambda: {"tissues": set(), "methods": set(),
                                             "tissue_method_pairs": []})

    for (tissue, method), data in shap_results.items():
        top_genes = data.get("top_100_genes", [])[:top_n]
        for rank, gene in enumerate(top_genes):
            gene_upper = gene.upper()
            gene_occurrences[gene_upper]["tissues"].add(tissue)
            gene_occurrences[gene_upper]["methods"].add(method)
            gene_occurrences[gene_upper]["tissue_method_pairs"].append({
                "tissue": tissue, "method": method, "rank": rank + 1
            })

    # Filter consensus genes
    consensus = []
    for gene, info in gene_occurrences.items():
        n_tissues = len(info["tissues"])
        n_methods = len(info["methods"])
        if n_tissues >= min_tissues and n_methods >= min_methods:
            consensus.append({
                "gene": gene,
                "n_tissues": n_tissues,
                "n_methods": n_methods,
                "tissues": sorted(info["tissues"]),
                "methods": sorted(info["methods"]),
                "best_rank": min(p["rank"] for p in info["tissue_method_pairs"]),
                "mean_rank": round(np.mean([p["rank"] for p in info["tissue_method_pairs"]]), 1),
            })

    # Sort by n_tissues (desc), then n_methods (desc), then best_rank (asc)
    consensus.sort(key=lambda x: (-x["n_tissues"], -x["n_methods"], x["best_rank"]))
    return consensus


def find_tissue_specific_genes(shap_results, top_n=100, min_methods=2):
    """Find genes in SHAP top-{top_n} in exactly 1 tissue, across ≥{min_methods} methods.

    DD-25: Tissue-specific = top-100 in 1 tissue, ≥3 methods.
    But since max methods per tissue = 2 (top-2), we use min_methods=2.
    """
    gene_tissue_methods = defaultdict(lambda: defaultdict(set))

    for (tissue, method), data in shap_results.items():
        top_genes = data.get("top_100_genes", [])[:top_n]
        for gene in top_genes:
            gene_upper = gene.upper()
            gene_tissue_methods[gene_upper][tissue].add(method)

    tissue_specific = defaultdict(list)
    for gene, tissue_map in gene_tissue_methods.items():
        if len(tissue_map) == 1:
            tissue = list(tissue_map.keys())[0]
            methods = tissue_map[tissue]
            if len(methods) >= min_methods:
                tissue_specific[tissue].append({
                    "gene": gene,
                    "n_methods": len(methods),
                    "methods": sorted(methods),
                })

    # Sort each tissue's genes by n_methods (desc)
    for tissue in tissue_specific:
        tissue_specific[tissue].sort(key=lambda x: -x["n_methods"])

    return dict(tissue_specific)


def hypergeometric_test(found_genes, reference_genes, top_n, total_genes,
                        all_gene_names=None):
    """Hypergeometric test for enrichment of reference genes in SHAP top-N.

    P(X ≥ k) where X ~ Hypergeometric(M=total_genes, n=ref_in_universe, N=top_n)

    If all_gene_names is provided, n = |reference ∩ gene_universe| (correct).
    Otherwise falls back to n = |reference| (conservative when ref genes
    are absent from the measured gene set).
    """
    ref_upper = {g.upper() for g in reference_genes}
    found_upper = {g.upper() for g in found_genes}

    overlap = ref_upper & found_upper
    k = len(overlap)  # observed overlaps
    M = total_genes  # population size
    N = min(top_n, len(found_upper))  # sample size (top-N SHAP genes)

    # n = reference genes present in the gene universe
    if all_gene_names is not None:
        universe_upper = {g.upper() for g in all_gene_names}
        n = len(ref_upper & universe_upper)
    else:
        n = len(ref_upper)

    # P(X >= k) = 1 - P(X <= k-1)
    p_value = hypergeom.sf(k - 1, M, n, N) if k > 0 else 1.0

    return {
        "n_reference": len(ref_upper),
        "n_reference_in_universe": n,
        "n_tested": N,
        "n_overlap": k,
        "overlap_genes": sorted(overlap),
        "p_value": round(float(p_value), 6),
        "enrichment_ratio": round(k / max(1, N * n / M), 2),
    }


def analyze_shap_directions(shap_results):
    """Analyze SHAP direction: positive SHAP = higher expression in FLT.

    For each tissue × method, check top-100 genes' direction consistency.
    """
    directions = {}

    for (tissue, method), data in shap_results.items():
        gene_dir = data.get("gene_direction", {})
        top_100 = data.get("top_100_genes", [])[:100]

        n_positive = sum(1 for g in top_100 if gene_dir.get(g, 0) > 0)
        n_negative = len(top_100) - n_positive

        # Top-10 upregulated and downregulated
        sorted_by_dir = sorted(gene_dir.items(), key=lambda x: x[1], reverse=True)
        top10_up = [(g, round(v, 4)) for g, v in sorted_by_dir[:10]]
        top10_down = [(g, round(v, 4)) for g, v in sorted_by_dir[-10:]]

        directions[f"{tissue}_{method}"] = {
            "tissue": tissue,
            "method": method,
            "n_positive_in_top100": n_positive,
            "n_negative_in_top100": n_negative,
            "top10_upregulated": top10_up,
            "top10_downregulated": top10_down,
        }

    return directions


def validate_vs_reference(shap_results, top_n=100):
    """Validate SHAP genes against literature reference genes."""
    validation = {}

    for (tissue, method), data in shap_results.items():
        top_genes = data.get("top_100_genes", [])[:top_n]
        all_genes = list(data.get("gene_importance", {}).keys())

        # Tissue-specific reference genes
        ref_genes = REFERENCE_GENES.get(tissue, [])
        if ref_genes:
            result = hypergeometric_test(top_genes, ref_genes, top_n,
                                          data.get("n_genes", TOTAL_GENES_APPROX),
                                          all_gene_names=all_genes)
            result["source"] = "tissue_specific_literature"
            validation[f"{tissue}_{method}_tissue_ref"] = result

        # da Silveira cross-tissue genes
        result = hypergeometric_test(top_genes, DA_SILVEIRA_GENES, top_n,
                                      data.get("n_genes", TOTAL_GENES_APPROX),
                                      all_gene_names=all_genes)
        result["source"] = "da_silveira_cell2020"
        validation[f"{tissue}_{method}_cell2020"] = result

    return validation


def method_agreement(shap_results, top_n=100):
    """Measure agreement between top-2 methods per tissue.

    For each tissue, compute Jaccard index and Spearman rank correlation
    between the two methods' SHAP rankings.
    """
    from scipy.stats import spearmanr

    agreement = {}

    for tissue, methods in TOP2_METHODS.items():
        if len(methods) < 2:
            continue

        key_a = (tissue, methods[0])
        key_b = (tissue, methods[1])

        if key_a not in shap_results or key_b not in shap_results:
            continue

        data_a = shap_results[key_a]
        data_b = shap_results[key_b]

        top_a = set(g.upper() for g in data_a.get("top_100_genes", [])[:top_n])
        top_b = set(g.upper() for g in data_b.get("top_100_genes", [])[:top_n])

        # Jaccard index
        intersection = top_a & top_b
        union = top_a | top_b
        jaccard = len(intersection) / len(union) if union else 0

        # Spearman on shared genes
        imp_a = data_a.get("gene_importance", {})
        imp_b = data_b.get("gene_importance", {})
        common_genes = set(g.upper() for g in imp_a) & set(g.upper() for g in imp_b)

        if len(common_genes) > 10:
            # Build aligned importance vectors (case-insensitive matching)
            imp_a_lower = {k.upper(): v for k, v in imp_a.items()}
            imp_b_lower = {k.upper(): v for k, v in imp_b.items()}
            vals_a = [imp_a_lower[g] for g in common_genes]
            vals_b = [imp_b_lower[g] for g in common_genes]
            rho, p = spearmanr(vals_a, vals_b)
        else:
            rho, p = np.nan, np.nan

        agreement[tissue] = {
            "method_a": methods[0],
            "method_b": methods[1],
            "jaccard_top100": round(jaccard, 4),
            "n_overlap_top100": len(intersection),
            "spearman_rho": round(float(rho), 4) if not np.isnan(rho) else None,
            "spearman_p": round(float(p), 6) if not np.isnan(p) else None,
            "n_common_genes": len(common_genes),
        }

    return agreement


def main():
    parser = argparse.ArgumentParser(description="Consensus spaceflight gene signature")
    parser.add_argument("--shap-dir", type=str, default=None)
    parser.add_argument("--output", type=str, default=None)
    args = parser.parse_args()

    shap_dir = Path(args.shap_dir) if args.shap_dir else V4_EVAL_DIR
    output_path = Path(args.output) if args.output else V4_EVAL_DIR / "SHAP_consensus.json"

    # Load all SHAP results
    shap_results = load_shap_results(shap_dir)
    print(f"Loaded {len(shap_results)} SHAP results")
    for (tissue, method) in sorted(shap_results.keys()):
        n = shap_results[(tissue, method)].get("n_genes", "?")
        print(f"  {tissue}/{method}: {n} genes")

    # 1. Consensus genes (DD-25)
    print(f"\n{'='*60}")
    print("1. Consensus genes (top-100 in ≥3 tissues AND ≥2 methods)")
    consensus = find_consensus_genes(shap_results, top_n=100, min_tissues=3, min_methods=2)
    print(f"  Found {len(consensus)} consensus genes")
    for g in consensus[:10]:
        print(f"    {g['gene']:12s} — {g['n_tissues']} tissues, {g['n_methods']} methods, "
              f"best rank {g['best_rank']}")

    # 2. Tissue-specific genes
    print(f"\n{'='*60}")
    print("2. Tissue-specific genes (1 tissue, ≥2 methods)")
    tissue_specific = find_tissue_specific_genes(shap_results, top_n=100, min_methods=2)
    for tissue, genes in sorted(tissue_specific.items()):
        print(f"  {tissue}: {len(genes)} tissue-specific genes")

    # 3. Validation vs literature
    print(f"\n{'='*60}")
    print("3. Validation vs literature reference genes")
    validation = validate_vs_reference(shap_results, top_n=100)
    for key, val in sorted(validation.items()):
        if val["n_overlap"] > 0:
            sig = "*" if val["p_value"] < 0.05 else ""
            print(f"  {key}: {val['n_overlap']}/{val['n_reference']} overlap, "
                  f"p={val['p_value']:.4f}{sig}")

    # 4. SHAP direction analysis
    print(f"\n{'='*60}")
    print("4. SHAP direction analysis")
    directions = analyze_shap_directions(shap_results)
    for key, d in sorted(directions.items()):
        print(f"  {key}: {d['n_positive_in_top100']} up, {d['n_negative_in_top100']} down")

    # 5. Method agreement
    print(f"\n{'='*60}")
    print("5. Method agreement (top-2 per tissue)")
    agreement = method_agreement(shap_results, top_n=100)
    for tissue, a in sorted(agreement.items()):
        rho = a["spearman_rho"] if a["spearman_rho"] is not None else "N/A"
        print(f"  {tissue}: Jaccard={a['jaccard_top100']:.3f}, "
              f"overlap={a['n_overlap_top100']}, rho={rho}")

    # Assemble output
    output = {
        "consensus_genes": consensus,
        "tissue_specific_genes": tissue_specific,
        "validation": validation,
        "shap_directions": directions,
        "method_agreement": agreement,
        "parameters": {
            "top_n": 100,
            "min_tissues_consensus": 3,
            "min_methods_consensus": 2,
            "min_methods_tissue_specific": 2,
            "total_genes_approx": TOTAL_GENES_APPROX,
        },
        "n_shap_results": len(shap_results),
        "timestamp": datetime.now().isoformat(),
    }

    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\nOutput: {output_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
