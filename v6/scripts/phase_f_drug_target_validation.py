#!/usr/bin/env python3
"""
Phase F: Drug Target Human Validation

Check which mouse-derived drug targets are also affected in human spaceflight cfRNA.
Classify into tiers:
  Tier A (validated): Mouse drug target + human DRR gene
  Tier B (promising): Mouse drug target + human DE (FDR<0.05) but not DRR
  Tier C (detected): In cfRNA but not DE
  Tier D (untranslatable): Not detected in cfRNA
"""

import sys
import os
import numpy as np
from scipy import stats
from datetime import datetime

sys.path.insert(0, os.path.dirname(__file__))
from v6_utils import (
    load_ortholog_map, load_cfrna_de, load_cfrna_drr,
    load_drug_targets, V6_EVAL, save_json
)


def main():
    print("=" * 60)
    print("Phase F: Drug Target Human Validation")
    print("=" * 60)

    # Load data
    print("\n1. Loading data...")
    orth_map = load_ortholog_map()
    cfrna_de = load_cfrna_de()
    drr_genes = load_cfrna_drr()
    drug_data = load_drug_targets()

    # Extract all drug-gene interactions from tier1_approved
    tiers = drug_data.get("tiers", {})
    tier1 = tiers.get("tier1_approved", [])
    tier3 = tiers.get("tier3_preclinical", [])

    print(f"  Tier 1 (FDA-approved) interactions: {len(tier1)}")
    print(f"  Tier 3 (pre-clinical) interactions: {len(tier3)}")

    # Get unique drug target genes (human symbols from DGIdb)
    all_interactions = tier1 + tier3
    target_genes = {}  # human_symbol → list of drugs
    for interaction in all_interactions:
        gene = interaction.get("gene", "")
        drug = interaction.get("drug", "")
        tier = "tier1_approved" if interaction in tier1 else "tier3_preclinical"
        if gene not in target_genes:
            target_genes[gene] = []
        target_genes[gene].append({
            "drug": drug,
            "tier": tier,
            "interaction_types": interaction.get("interaction_types", []),
        })

    print(f"  Unique drug target genes: {len(target_genes)}")

    # Check each target gene in cfRNA
    print("\n2. Checking drug targets in human cfRNA...")
    gene_tiers = {"A": [], "B": [], "C": [], "D": []}
    tier_counts = {"A": 0, "B": 0, "C": 0, "D": 0}

    for human_gene, drugs in target_genes.items():
        if human_gene in cfrna_de.index:
            row = cfrna_de.loc[human_gene]
            fdr = float(row.get("edge_pre_vs_flight_fdr", 1.0))
            diff = float(row.get("edge_pre_vs_flight_diff", 0.0))
            is_drr = human_gene in drr_genes
            is_de = fdr < 0.05

            entry = {
                "gene": human_gene,
                "n_drugs": len(drugs),
                "top_drugs": [d["drug"] for d in drugs[:5]],
                "has_fda_approved": any(d["tier"] == "tier1_approved" for d in drugs),
                "de_fdr": round(fdr, 6),
                "de_diff": round(diff, 4),
                "direction": "up_in_flight" if diff > 0 else "down_in_flight",
            }

            if is_drr:
                gene_tiers["A"].append(entry)
                tier_counts["A"] += 1
            elif is_de:
                gene_tiers["B"].append(entry)
                tier_counts["B"] += 1
            else:
                gene_tiers["C"].append(entry)
                tier_counts["C"] += 1
        else:
            gene_tiers["D"].append({
                "gene": human_gene,
                "n_drugs": len(drugs),
                "top_drugs": [d["drug"] for d in drugs[:3]],
                "has_fda_approved": any(d["tier"] == "tier1_approved" for d in drugs),
            })
            tier_counts["D"] += 1

    print(f"\n  Tier A (validated — drug target + DRR): {tier_counts['A']}")
    print(f"  Tier B (promising — drug target + DE):  {tier_counts['B']}")
    print(f"  Tier C (detected — in cfRNA, not DE):   {tier_counts['C']}")
    print(f"  Tier D (untranslatable — not in cfRNA): {tier_counts['D']}")

    # Sort tiers by effect size
    gene_tiers["A"].sort(key=lambda x: abs(x.get("de_diff", 0)), reverse=True)
    gene_tiers["B"].sort(key=lambda x: abs(x.get("de_diff", 0)), reverse=True)

    # Print top Tier A candidates
    if gene_tiers["A"]:
        print("\n  Top Tier A countermeasure candidates:")
        for entry in gene_tiers["A"][:10]:
            dir_str = "↑" if entry["direction"] == "up_in_flight" else "↓"
            print(f"    {entry['gene']:>10s} {dir_str} (diff={entry['de_diff']:+.3f}) "
                  f"→ {', '.join(entry['top_drugs'][:3])}")

    if gene_tiers["B"]:
        print("\n  Top Tier B candidates:")
        for entry in gene_tiers["B"][:10]:
            dir_str = "↑" if entry["direction"] == "up_in_flight" else "↓"
            print(f"    {entry['gene']:>10s} {dir_str} (diff={entry['de_diff']:+.3f}, "
                  f"FDR={entry['de_fdr']:.4f}) → {', '.join(entry['top_drugs'][:3])}")

    # Fisher exact test: drug targets enriched among DRR?
    print("\n3. Fisher exact test: drug targets enriched among DRR...")
    human_genes_in_cfrna = set(cfrna_de.index)
    drug_target_set = set(target_genes.keys()) & human_genes_in_cfrna
    non_drug_target = human_genes_in_cfrna - drug_target_set
    drr_in_cfrna = drr_genes & human_genes_in_cfrna

    # Contingency table:
    #                  DRR    non-DRR
    # Drug target       a       b
    # Non-drug target   c       d
    a = len(drug_target_set & drr_in_cfrna)
    b = len(drug_target_set - drr_in_cfrna)
    c = len(non_drug_target & drr_in_cfrna)
    d = len(non_drug_target - drr_in_cfrna)

    odds_ratio, p_fisher = stats.fisher_exact([[a, b], [c, d]], alternative="greater")

    print(f"  Contingency table:")
    print(f"              DRR    non-DRR")
    print(f"  Drug target  {a:>4d}   {b:>6d}")
    print(f"  Non-target   {c:>4d}   {d:>6d}")
    print(f"  Odds ratio: {odds_ratio:.3f}")
    print(f"  Fisher p (one-sided): {p_fisher:.6f}")

    # Per-drug summary: which drugs target DRR genes?
    print("\n4. Building countermeasure table...")
    countermeasure_drugs = {}
    for entry in gene_tiers["A"]:
        gene = entry["gene"]
        for drug_info in target_genes.get(gene, []):
            drug_name = drug_info["drug"]
            if drug_name not in countermeasure_drugs:
                countermeasure_drugs[drug_name] = {
                    "drug": drug_name,
                    "tier": drug_info["tier"],
                    "target_genes_drr": [],
                }
            countermeasure_drugs[drug_name]["target_genes_drr"].append(gene)

    # Sort by number of DRR targets
    countermeasure_list = sorted(
        countermeasure_drugs.values(),
        key=lambda x: len(x["target_genes_drr"]), reverse=True)

    fda_countermeasures = [d for d in countermeasure_list
                          if d["tier"] == "tier1_approved"]
    print(f"  FDA-approved drugs targeting DRR genes: {len(fda_countermeasures)}")
    for d in fda_countermeasures[:10]:
        targets = ", ".join(d["target_genes_drr"][:5])
        print(f"    {d['drug']:>25s} → {targets}")

    # Assemble output
    output = {
        "n_drug_target_genes": len(target_genes),
        "n_in_cfrna": len(drug_target_set),
        "tier_counts": tier_counts,
        "tier_A_validated": gene_tiers["A"],
        "tier_B_promising": gene_tiers["B"],
        "tier_C_detected": len(gene_tiers["C"]),  # just count (large)
        "tier_D_untranslatable": len(gene_tiers["D"]),
        "fisher_test": {
            "contingency": [[a, b], [c, d]],
            "odds_ratio": round(float(odds_ratio), 4),
            "p_value": float(p_fisher),
            "significant": p_fisher < 0.05,
        },
        "countermeasure_drugs": countermeasure_list[:50],  # top 50
        "fda_countermeasures_n": len(fda_countermeasures),
        "timestamp": datetime.now().isoformat(),
    }

    save_json(output, V6_EVAL / "V6_F_drug_target_validation.json")
    print("\nDone!")


if __name__ == "__main__":
    main()
