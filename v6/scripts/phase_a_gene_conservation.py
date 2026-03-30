#!/usr/bin/env python3
"""
Phase A: Gene-Level Conservation — Mouse SHAP genes vs Human cfRNA DRR genes

Tests whether mouse spaceflight-important genes (SHAP, WGCNA hubs, biomarker panel)
overlap with human differentially-regulated response (DRR) genes from cfRNA.
"""

import sys
import os
import numpy as np
from scipy import stats
from datetime import datetime

sys.path.insert(0, os.path.dirname(__file__))
from v6_utils import (
    load_ortholog_map, load_ensembl_to_symbol, map_mouse_to_human,
    load_cfrna_de, load_cfrna_drr, load_shap_consensus,
    load_all_shap_top_genes, load_wgcna_hub_genes, load_biomarker_panel,
    V6_EVAL, save_json
)


def hypergeometric_test(overlap_size, set_size, drr_size, universe_size):
    """One-sided hypergeometric test for enrichment.
    P(X >= overlap_size) where X ~ Hypergeometric(N, K, n).
    """
    # scipy: sf(k-1, M, n, N) = P(X >= k)
    p = stats.hypergeom.sf(overlap_size - 1, universe_size, drr_size, set_size)
    return p


def permutation_test(gene_set, drr_set, universe, n_perm=10000, seed=42):
    """Empirical permutation test: random gene sets of same size."""
    rng = np.random.default_rng(seed)
    universe_list = list(universe)
    observed = len(gene_set & drr_set)
    n_genes = len(gene_set)

    count_ge = 0
    for _ in range(n_perm):
        random_genes = set(rng.choice(universe_list, size=n_genes, replace=False))
        if len(random_genes & drr_set) >= observed:
            count_ge += 1

    return (count_ge + 1) / (n_perm + 1)  # add 1 for conservative estimate


def direction_concordance(mouse_genes_human, cfrna_de, shap_data_all):
    """Check if mouse SHAP direction matches human cfRNA direction.
    Returns list of {gene, mouse_direction, human_direction, concordant}.
    """
    results = []
    de_col = "edge_pre_vs_flight_diff"

    for mouse_gene, human_gene in mouse_genes_human.items():
        if human_gene not in cfrna_de.index:
            continue

        human_diff = cfrna_de.loc[human_gene, de_col]
        human_dir = "up" if human_diff > 0 else "down"

        # Get mouse direction from SHAP consensus or tissue-specific
        mouse_dir = None
        consensus = shap_data_all.get("shap_directions", {})
        if mouse_gene in consensus:
            d = consensus[mouse_gene]
            if isinstance(d, dict):
                # Take majority direction across tissues
                ups = sum(1 for v in d.values() if v == "up_in_flight")
                mouse_dir = "up" if ups > len(d) / 2 else "down"
            elif isinstance(d, str):
                mouse_dir = "up" if "up" in d else "down"

        if mouse_dir is None:
            continue

        results.append({
            "mouse_gene": mouse_gene,
            "human_gene": human_gene,
            "mouse_direction": mouse_dir,
            "human_direction": human_dir,
            "human_diff": round(float(human_diff), 4),
            "concordant": mouse_dir == human_dir,
        })

    return results


def main():
    print("=" * 60)
    print("Phase A: Gene-Level Conservation")
    print("=" * 60)

    # Load data
    print("\n1. Loading data...")
    orth_map = load_ortholog_map()
    ens2sym = load_ensembl_to_symbol()
    cfrna_de = load_cfrna_de()
    drr_genes = load_cfrna_drr()
    shap_consensus = load_shap_consensus()
    shap_all = load_all_shap_top_genes(top_n=100)
    wgcna_hubs = load_wgcna_hub_genes()
    biomarker = load_biomarker_panel()

    print(f"  Orthologs: {len(orth_map)} pairs")
    print(f"  cfRNA DE: {len(cfrna_de)} genes")
    print(f"  DRR genes: {len(drr_genes)}")

    # Build universe: human genes with mouse ortholog AND present in cfRNA
    human_genes_in_cfrna = set(cfrna_de.index)
    mouse_to_human = orth_map  # already mouse_symbol → human_symbol
    human_genes_with_orth = set(mouse_to_human.values())
    universe = human_genes_in_cfrna & human_genes_with_orth
    print(f"  Universe (orthologs ∩ cfRNA): {len(universe)}")
    drr_in_universe = drr_genes & universe
    print(f"  DRR in universe: {len(drr_in_universe)}")

    # Prepare gene sets
    gene_sets = {}

    # 1. SHAP consensus (5 genes)
    consensus_genes = [g["gene"] for g in shap_consensus.get("consensus_genes", [])]
    mapped_consensus, unmapped_consensus = map_mouse_to_human(
        consensus_genes, orth_map, ens2sym)
    gene_sets["shap_consensus"] = {
        "mouse_genes": consensus_genes,
        "human_genes": set(mapped_consensus.values()) & universe,
        "unmapped": unmapped_consensus,
        "description": "SHAP consensus (≥3 tissues, ≥2 methods)",
    }

    # 2. SHAP union top-100 (all tissues combined)
    all_shap = set()
    for tissue_genes in shap_all.values():
        all_shap |= tissue_genes
    mapped_shap_union, unmapped_shap_union = map_mouse_to_human(
        list(all_shap), orth_map, ens2sym)
    gene_sets["shap_union_top100"] = {
        "mouse_genes_n": len(all_shap),
        "human_genes": set(mapped_shap_union.values()) & universe,
        "unmapped_n": len(unmapped_shap_union),
        "description": "Union of SHAP top-100 across all 8 tissues",
    }

    # 3. Per-tissue SHAP top-100
    per_tissue_results = {}
    for tissue, tissue_genes in shap_all.items():
        mapped, unmapped = map_mouse_to_human(list(tissue_genes), orth_map, ens2sym)
        human_set = set(mapped.values()) & universe
        per_tissue_results[tissue] = {
            "mouse_n": len(tissue_genes),
            "human_n": len(human_set),
            "unmapped_n": len(unmapped),
            "human_genes": human_set,
        }

    # 4. WGCNA hub genes
    all_hubs = set()
    if isinstance(wgcna_hubs, dict):
        # Structure: {"tissues": {tissue: {"modules": {color: {"hub_genes": [...]}}}}}
        tissues_data = wgcna_hubs.get("tissues", wgcna_hubs)
        for tissue_name, tissue_data in tissues_data.items():
            if not isinstance(tissue_data, dict):
                continue
            modules = tissue_data.get("modules", {})
            for mod_color, mod_data in modules.items():
                if isinstance(mod_data, dict) and "hub_genes" in mod_data:
                    for g in mod_data["hub_genes"]:
                        # Handle pipe-separated gene names (e.g., "Rpl34|Rpl34-ps1")
                        all_hubs.add(g.split("|")[0])
    mapped_hubs, unmapped_hubs = map_mouse_to_human(list(all_hubs), orth_map, ens2sym)
    gene_sets["wgcna_hubs"] = {
        "mouse_genes_n": len(all_hubs),
        "human_genes": set(mapped_hubs.values()) & universe,
        "unmapped_n": len(unmapped_hubs),
        "description": "WGCNA hub genes (kME>0.8, 6 tissues)",
    }

    # 5. Biomarker panel (20 genes)
    panel_genes = []
    panel_human = []
    for g in biomarker.get("panel", []):
        mouse_name = g.get("gene", "")
        human_name = g.get("human_symbol", "")
        panel_genes.append(mouse_name)
        # Check if human_symbol is actually a valid human gene
        if human_name and human_name in universe:
            panel_human.append(human_name)
        else:
            # Try mapping through orthologs
            mapped, _ = map_mouse_to_human([mouse_name], orth_map, ens2sym)
            for h in mapped.values():
                if h in universe:
                    panel_human.append(h)
    gene_sets["biomarker_panel"] = {
        "mouse_genes": panel_genes,
        "human_genes": set(panel_human),
        "description": "v5 consensus biomarker panel (20 genes)",
    }

    # Run enrichment tests
    print("\n2. Running enrichment tests...")
    enrichment_results = {}

    for name, gs in gene_sets.items():
        human_set = gs["human_genes"]
        n_set = len(human_set)
        overlap = human_set & drr_in_universe
        n_overlap = len(overlap)

        if n_set == 0:
            enrichment_results[name] = {
                "n_human_genes": 0, "n_overlap": 0,
                "note": "No mappable genes in universe"
            }
            continue

        # Hypergeometric test
        p_hyper = hypergeometric_test(
            n_overlap, n_set, len(drr_in_universe), len(universe))

        # Fold enrichment
        expected = n_set * len(drr_in_universe) / len(universe)
        fold_enrich = n_overlap / expected if expected > 0 else 0

        # Permutation test
        p_perm = permutation_test(human_set, drr_in_universe, universe,
                                  n_perm=10000)

        enrichment_results[name] = {
            "description": gs["description"],
            "n_mouse_genes": gs.get("mouse_genes_n", len(gs.get("mouse_genes", []))),
            "n_human_genes": n_set,
            "n_overlap_drr": n_overlap,
            "overlap_genes": sorted(overlap),
            "expected_overlap": round(expected, 2),
            "fold_enrichment": round(fold_enrich, 3),
            "hypergeometric_p": float(p_hyper),
            "permutation_p": float(p_perm),
            "significant": p_hyper < 0.05,
        }

        print(f"  {name}: {n_set} genes, {n_overlap} DRR overlap, "
              f"FE={fold_enrich:.2f}, p_hyper={p_hyper:.4f}, p_perm={p_perm:.4f}")

    # Per-tissue enrichment
    print("\n3. Per-tissue SHAP-DRR overlap...")
    per_tissue_enrichment = {}
    for tissue, td in per_tissue_results.items():
        human_set = td["human_genes"]
        n_set = len(human_set)
        overlap = human_set & drr_in_universe
        n_overlap = len(overlap)

        if n_set == 0:
            per_tissue_enrichment[tissue] = {"n_human_genes": 0, "n_overlap": 0}
            continue

        p_hyper = hypergeometric_test(
            n_overlap, n_set, len(drr_in_universe), len(universe))
        expected = n_set * len(drr_in_universe) / len(universe)
        fold_enrich = n_overlap / expected if expected > 0 else 0

        per_tissue_enrichment[tissue] = {
            "n_mouse_genes": td["mouse_n"],
            "n_human_genes": n_set,
            "n_unmapped": td["unmapped_n"],
            "n_overlap_drr": n_overlap,
            "overlap_genes": sorted(overlap),
            "expected_overlap": round(expected, 2),
            "fold_enrichment": round(fold_enrich, 3),
            "hypergeometric_p": float(p_hyper),
            "significant": p_hyper < 0.05,
        }

        sig = "*" if p_hyper < 0.05 else ""
        print(f"  {tissue}: {n_set} human genes, {n_overlap} DRR overlap, "
              f"FE={fold_enrich:.2f}, p={p_hyper:.4f}{sig}")

    # Direction concordance
    print("\n4. Direction concordance...")
    # Use union of all SHAP mapped genes
    all_mapped_genes = {}
    for tissue_genes in shap_all.values():
        m, _ = map_mouse_to_human(list(tissue_genes), orth_map, ens2sym)
        all_mapped_genes.update(m)
    directions = direction_concordance(all_mapped_genes, cfrna_de, shap_consensus)
    n_concordant = sum(1 for d in directions if d["concordant"])
    n_total = len(directions)
    concordance_rate = n_concordant / n_total if n_total > 0 else 0
    print(f"  Direction concordance: {n_concordant}/{n_total} ({concordance_rate:.1%})")

    # Catalog untranslatable genes
    print("\n5. Cataloging untranslatable genes...")
    untranslatable = {
        "consensus_unmapped": unmapped_consensus,
        "notable_mouse_specific": [],
    }
    # Check known mouse-specific genes
    mouse_specific = ["MUP22", "Mup22", "GM24497", "Gm24497"]
    for g in mouse_specific:
        if g in all_shap or any(g in tissue_genes for tissue_genes in shap_all.values()):
            if g not in orth_map and g.capitalize() not in orth_map:
                untranslatable["notable_mouse_specific"].append(g)

    print(f"  Consensus unmapped: {unmapped_consensus}")
    print(f"  Notable mouse-specific: {untranslatable['notable_mouse_specific']}")

    # Assemble output
    output = {
        "universe_size": len(universe),
        "drr_in_universe": len(drr_in_universe),
        "enrichment_results": enrichment_results,
        "per_tissue_enrichment": per_tissue_enrichment,
        "direction_concordance": {
            "n_tested": n_total,
            "n_concordant": n_concordant,
            "concordance_rate": round(concordance_rate, 4),
            "top_concordant": sorted(
                [d for d in directions if d["concordant"]],
                key=lambda x: abs(x["human_diff"]), reverse=True)[:20],
            "top_discordant": sorted(
                [d for d in directions if not d["concordant"]],
                key=lambda x: abs(x["human_diff"]), reverse=True)[:20],
        },
        "untranslatable": untranslatable,
        "timestamp": datetime.now().isoformat(),
    }

    save_json(output, V6_EVAL / "V6_A_gene_conservation.json")
    print("\nDone!")


if __name__ == "__main__":
    main()
