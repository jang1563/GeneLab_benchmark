#!/usr/bin/env python
"""GeneLabBench v5 Phase 4: Drug Target Mapping.

Queries DGIdb v5 (GraphQL) and ChEMBL REST API to find FDA-approved
drugs and clinical compounds targeting spaceflight-affected genes.

Usage:
    python drug_target_mapping.py
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests


class SafeEncoder(json.JSONEncoder):
    """Handle numpy and other non-standard types for JSON serialization."""
    def default(self, obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, (np.bool_,)):
            return bool(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, set):
            return sorted(obj)
        return super().default(obj)

# ── paths ──
BASE_DIR = Path(__file__).resolve().parent.parent.parent
V4_EVAL_DIR = BASE_DIR / "v4" / "evaluation"
V5_EVAL_DIR = BASE_DIR / "v5" / "evaluation"
ORTHO_PATH = BASE_DIR / "v3" / "data" / "mouse_human_orthologs.tsv"

# DGIdb v5 GraphQL endpoint
DGIDB_URL = "https://dgidb.org/api/graphql"
# ChEMBL REST API
CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"


def load_orthologs():
    """Load mouse→human ortholog mapping."""
    df = pd.read_csv(ORTHO_PATH, sep="\t")
    return dict(zip(df["mouse_symbol"], df["human_symbol"]))


def collect_target_genes():
    """Collect all spaceflight-affected genes from SHAP + WGCNA."""
    genes = {}  # gene_symbol → {tissues, sources, max_shap_rank}

    # 1. SHAP consensus genes
    consensus_path = V4_EVAL_DIR / "SHAP_consensus.json"
    if consensus_path.exists():
        with open(consensus_path) as f:
            consensus = json.load(f)
        for entry in consensus.get("consensus_genes", []):
            g = entry["gene"]
            genes[g] = {
                "tissues": set(entry.get("tissues", [])),
                "sources": {"shap_consensus"},
                "shap_rank": entry.get("best_rank", 999),
            }
        for tissue, entries in consensus.get("tissue_specific_genes", {}).items():
            for entry in entries:
                g = entry["gene"]
                if g not in genes:
                    genes[g] = {"tissues": set(), "sources": set(), "shap_rank": 999}
                genes[g]["tissues"].add(tissue)
                genes[g]["sources"].add("shap_tissue_specific")

    # 2. SHAP per-tissue top-100
    for path in V4_EVAL_DIR.glob("SHAP_*_*.json"):
        if "consensus" in path.name or "interaction" in path.name or "WGCNA" in path.name:
            continue
        try:
            with open(path) as f:
                data = json.load(f)
            tissue = data.get("tissue", "")
            importance = data.get("mean_importance", data.get("gene_importance", {}))
            if not importance:
                continue
            ranked = sorted(importance.items(), key=lambda x: abs(x[1]), reverse=True)
            for rank, (g, _) in enumerate(ranked[:100]):
                if g not in genes:
                    genes[g] = {"tissues": set(), "sources": set(), "shap_rank": 999}
                genes[g]["tissues"].add(tissue)
                genes[g]["sources"].add("shap_top100")
                genes[g]["shap_rank"] = min(genes[g]["shap_rank"], rank)
        except (json.JSONDecodeError, KeyError):
            continue

    # 3. WGCNA hub genes
    hub_path = V4_EVAL_DIR / "WGCNA_hub_genes.json"
    if hub_path.exists():
        with open(hub_path) as f:
            hub_data = json.load(f)
        for tissue, tdata in hub_data.get("tissues", {}).items():
            for module, mdata in tdata.get("modules", {}).items():
                for g in mdata.get("hub_genes", []):
                    if g not in genes:
                        genes[g] = {"tissues": set(), "sources": set(), "shap_rank": 999}
                    genes[g]["tissues"].add(tissue)
                    genes[g]["sources"].add("wgcna_hub")

    # Convert sets to lists for JSON
    for g in genes:
        genes[g]["tissues"] = sorted(genes[g]["tissues"])
        genes[g]["sources"] = sorted(genes[g]["sources"])

    print(f"  Total unique target genes: {len(genes)}")
    return genes


def query_dgidb_graphql(human_genes, batch_size=100):
    """Query DGIdb v5 GraphQL API for drug-gene interactions."""
    all_interactions = {}

    gene_list = list(human_genes)
    n_batches = (len(gene_list) + batch_size - 1) // batch_size

    for i in range(n_batches):
        batch = gene_list[i*batch_size : (i+1)*batch_size]

        query = """
        query($names: [String!]!) {
          genes(names: $names) {
            nodes {
              name
              interactions {
                drug {
                  name
                  approved
                  conceptId
                }
                interactionScore
                interactionTypes {
                  type
                  directionality
                }
                publications {
                  pmid
                }
              }
            }
          }
        }
        """

        variables = {"names": batch}

        try:
            resp = requests.post(
                DGIDB_URL,
                json={"query": query, "variables": variables},
                headers={"Content-Type": "application/json"},
                timeout=30,
            )
            resp.raise_for_status()
            data = resp.json()

            nodes = data.get("data", {}).get("genes", {}).get("nodes", [])
            for node in nodes:
                gene_name = node.get("name", "")
                for interaction in node.get("interactions", []):
                    drug = interaction.get("drug", {})
                    drug_name = drug.get("name", "Unknown")
                    approved = drug.get("approved", False)

                    int_types = [t.get("type", "") for t in
                                 interaction.get("interactionTypes", [])]
                    pmids = [p.get("pmid", "") for p in
                             interaction.get("publications", [])[:5]]

                    if gene_name not in all_interactions:
                        all_interactions[gene_name] = []
                    all_interactions[gene_name].append({
                        "drug": drug_name,
                        "approved": approved,
                        "interaction_types": int_types,
                        "score": interaction.get("interactionScore"),
                        "pmids": pmids,
                    })

            print(f"    DGIdb batch {i+1}/{n_batches}: {len(nodes)} genes returned")

        except Exception as e:
            print(f"    DGIdb batch {i+1}/{n_batches} failed: {e}")

        time.sleep(1)  # Rate limit

    return all_interactions


def query_chembl_targets(human_genes, top_n=50):
    """Query ChEMBL REST API for top-N genes by SHAP importance."""
    results = {}
    genes_to_query = list(human_genes)[:top_n]

    for gene in genes_to_query:
        try:
            # Search for target by gene name
            url = f"{CHEMBL_BASE}/target/search.json"
            params = {"q": gene, "limit": 5, "format": "json"}
            resp = requests.get(url, params=params, timeout=10)

            if resp.status_code != 200:
                continue

            data = resp.json()
            targets = data.get("targets", [])
            if not targets:
                continue

            target = targets[0]  # Best match
            target_id = target.get("target_chembl_id", "")

            # Get approved drugs for this target
            drug_url = f"{CHEMBL_BASE}/mechanism.json"
            drug_params = {"target_chembl_id": target_id, "limit": 10, "format": "json"}
            drug_resp = requests.get(drug_url, params=drug_params, timeout=10)

            drugs = []
            if drug_resp.status_code == 200:
                mechs = drug_resp.json().get("mechanisms", [])
                for mech in mechs:
                    drugs.append({
                        "drug": mech.get("molecule_chembl_id", ""),
                        "drug_name": mech.get("molecule_name", ""),
                        "mechanism": mech.get("mechanism_of_action", ""),
                        "max_phase": mech.get("max_phase"),
                    })

            results[gene] = {
                "target_chembl_id": target_id,
                "target_type": target.get("target_type", ""),
                "organism": target.get("organism", ""),
                "n_drugs": len(drugs),
                "drugs": drugs,
            }

        except Exception as e:
            continue

        time.sleep(0.5)  # Rate limit

    print(f"  ChEMBL: queried {len(genes_to_query)} genes, {len(results)} with targets")
    return results


def categorize_drugs(dgidb_results, chembl_results):
    """Categorize drug hits into tiers. Deduplicates by (gene, drug) key."""
    # Collect all entries with tier assignment
    all_entries = {}  # (gene, drug_upper) → (tier, entry)

    # From DGIdb
    for gene, interactions in dgidb_results.items():
        for ix in interactions:
            drug_name = ix["drug"]
            key = (gene, drug_name.upper())
            tier = 1 if ix.get("approved") else 3
            entry = {
                "gene": gene,
                "drug": drug_name,
                "interaction_types": ix.get("interaction_types", []),
                "source": "DGIdb",
            }
            # Keep best (lowest) tier for duplicates
            if key not in all_entries or tier < all_entries[key][0]:
                all_entries[key] = (tier, entry)

    # From ChEMBL
    for gene, data in chembl_results.items():
        for drug in data.get("drugs", []):
            drug_name = drug.get("drug_name", drug.get("drug", ""))
            if not drug_name:
                continue
            key = (gene, drug_name.upper())
            phase = drug.get("max_phase")
            if phase and phase >= 4:
                tier = 1
            elif phase and phase >= 2:
                tier = 2
            else:
                tier = 3
            entry = {
                "gene": gene,
                "drug": drug_name,
                "mechanism": drug.get("mechanism", ""),
                "source": "ChEMBL",
            }
            if key not in all_entries or tier < all_entries[key][0]:
                all_entries[key] = (tier, entry)

    # Split by tier
    tier1, tier2, tier3 = [], [], []
    for (tier, entry) in all_entries.values():
        if tier == 1:
            tier1.append(entry)
        elif tier == 2:
            tier2.append(entry)
        else:
            tier3.append(entry)

    return {
        "tier1_approved": tier1,
        "tier2_clinical": tier2,
        "tier3_preclinical": tier3[:200],  # Cap to avoid huge output
    }


def wgcna_module_enrichment(target_genes_with_drugs, hub_data):
    """Test if drug targets are enriched in specific WGCNA modules.

    Uses all available gene lists from WGCNA hub data (hub_genes, top20_genes,
    all_genes) to maximize overlap detection.
    """
    from scipy.stats import fisher_exact

    druggable_genes = set(target_genes_with_drugs)
    enrichment_results = {}

    for tissue, tdata in hub_data.get("tissues", {}).items():
        n_total = tdata.get("n_genes_input", 5000)
        enrichment_results[tissue] = {}

        for module, mdata in tdata.get("modules", {}).items():
            # Collect ALL genes in this module from any available field
            module_genes = set()
            for field in ["hub_genes", "top20_genes", "all_genes", "genes"]:
                module_genes.update(mdata.get(field, []))
            # Also check kME dict keys (gene→kME score)
            for field in ["kME", "kme"]:
                if isinstance(mdata.get(field), dict):
                    module_genes.update(mdata[field].keys())

            n_module = mdata.get("n_genes", len(module_genes))
            if n_module < 1:
                continue

            overlap = druggable_genes & module_genes
            n_overlap = len(overlap)

            if n_overlap == 0:
                continue

            # Fisher exact test (2×2 table)
            n_drug = len(druggable_genes)
            table = [
                [n_overlap, n_drug - n_overlap],
                [n_module - n_overlap, max(0, n_total - n_module - n_drug + n_overlap)]
            ]
            table = [[max(0, x) for x in row] for row in table]

            try:
                odds, p = fisher_exact(table, alternative="greater")
            except ValueError:
                continue

            enrichment_results[tissue][module] = {
                "n_overlap": n_overlap,
                "n_module": n_module,
                "n_druggable": n_drug,
                "overlap_genes": sorted(overlap),
                "odds_ratio": round(float(odds), 2),
                "fisher_p": round(float(p), 6),
                "significant": bool(p < 0.05),
            }

    return enrichment_results


def main():
    V5_EVAL_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Phase 4: Drug Target Mapping")
    print("=" * 60)

    # Step 1: Collect target genes
    print("\n1. Collecting target genes from SHAP + WGCNA...")
    target_genes = collect_target_genes()

    # Step 2: Map mouse → human symbols
    print("\n2. Mapping mouse → human orthologs...")
    ortho = load_orthologs()
    mouse_to_human = {}
    for mouse_sym in target_genes:
        human_sym = ortho.get(mouse_sym)
        if human_sym:
            mouse_to_human[mouse_sym] = human_sym
    print(f"  Mapped {len(mouse_to_human)}/{len(target_genes)} genes to human orthologs")

    human_genes = set(mouse_to_human.values())
    human_to_mouse = {v: k for k, v in mouse_to_human.items()}

    # Step 3: Query DGIdb v5 (GraphQL)
    print("\n3. Querying DGIdb v5 (GraphQL)...")
    dgidb_results = query_dgidb_graphql(human_genes)
    n_with_drugs = len(dgidb_results)
    n_total_interactions = sum(len(v) for v in dgidb_results.values())
    print(f"  DGIdb: {n_with_drugs} genes with interactions, "
          f"{n_total_interactions} total interactions")

    # Step 4: Query ChEMBL (top-50 by SHAP rank)
    print("\n4. Querying ChEMBL (top-50 genes)...")
    # Sort by SHAP rank
    ranked_human = sorted(
        [(h, target_genes[m]["shap_rank"]) for m, h in mouse_to_human.items()],
        key=lambda x: x[1]
    )
    top50_human = [h for h, _ in ranked_human[:50]]
    chembl_results = query_chembl_targets(top50_human, top_n=50)

    # Step 5: Categorize into tiers
    print("\n5. Categorizing drug hits...")
    tiers = categorize_drugs(dgidb_results, chembl_results)
    print(f"  Tier 1 (FDA-approved): {len(tiers['tier1_approved'])}")
    print(f"  Tier 2 (Clinical): {len(tiers['tier2_clinical'])}")
    print(f"  Tier 3 (Pre-clinical): {len(tiers['tier3_preclinical'])}")

    # Step 6: WGCNA module enrichment
    print("\n6. Testing WGCNA module enrichment...")
    hub_path = V4_EVAL_DIR / "WGCNA_hub_genes.json"
    if hub_path.exists():
        with open(hub_path) as f:
            hub_data = json.load(f)
        druggable_mouse = {human_to_mouse.get(g, g) for g in dgidb_results.keys()}
        module_enrichment = wgcna_module_enrichment(druggable_mouse, hub_data)
        n_enriched = sum(len(v) for v in module_enrichment.values())
        print(f"  Enriched modules (p<0.05): {n_enriched}")
    else:
        module_enrichment = {}

    # Step 7: Per-tissue druggability summary
    print("\n7. Computing per-tissue druggability...")
    tissue_druggability = {}
    for mouse_sym, info in target_genes.items():
        human_sym = mouse_to_human.get(mouse_sym)
        has_drug = human_sym in dgidb_results or human_sym in chembl_results
        for tissue in info["tissues"]:
            if tissue not in tissue_druggability:
                tissue_druggability[tissue] = {"total": 0, "druggable": 0, "genes": []}
            tissue_druggability[tissue]["total"] += 1
            if has_drug:
                tissue_druggability[tissue]["druggable"] += 1
                tissue_druggability[tissue]["genes"].append(mouse_sym)

    for tissue in tissue_druggability:
        total = tissue_druggability[tissue]["total"]
        druggable = tissue_druggability[tissue]["druggable"]
        pct = druggable / total * 100 if total > 0 else 0
        tissue_druggability[tissue]["pct_druggable"] = round(pct, 1)
        tissue_druggability[tissue]["genes"] = tissue_druggability[tissue]["genes"][:20]
        print(f"  {tissue}: {druggable}/{total} ({pct:.1f}%) druggable")

    # Save results
    output = {
        "n_target_genes_mouse": len(target_genes),
        "n_mapped_to_human": len(mouse_to_human),
        "n_genes_with_dgidb_interactions": n_with_drugs,
        "n_total_dgidb_interactions": n_total_interactions,
        "n_chembl_targets": len(chembl_results),
        "tiers": tiers,
        "tissue_druggability": tissue_druggability,
        "wgcna_module_enrichment": module_enrichment,
        "dgidb_summary": {g: len(v) for g, v in dgidb_results.items()},
        "chembl_summary": {g: d.get("n_drugs", 0) for g, d in chembl_results.items()},
    }

    out_path = V5_EVAL_DIR / "drug_targets.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, cls=SafeEncoder)
    print(f"\nSaved: {out_path}")
    print("Done!")


if __name__ == "__main__":
    main()
