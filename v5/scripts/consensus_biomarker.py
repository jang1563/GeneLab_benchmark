#!/usr/bin/env python
"""GeneLabBench v5 Phase 5: Consensus Biomarker Panel.

Integrates all v4-v5 evidence sources to nominate a consensus spaceflight
biomarker panel. Scores each gene by multi-source evidence (max 12 points),
then validates the panel via PCA-LR on v4 LOMO folds.

Scoring: SHAP(0-3) + WGCNA(0-2) + multi-tissue(0-2) + network(0-1)
         + druggable(0-1) + literature(0-2) + cross-species(0-1) = max 12

Usage:
    python consensus_biomarker.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

# ── paths ──
BASE_DIR = Path(__file__).resolve().parent.parent.parent
V4_EVAL_DIR = BASE_DIR / "v4" / "evaluation"
V5_EVAL_DIR = BASE_DIR / "v5" / "evaluation"
ORTHO_PATH = BASE_DIR / "v3" / "data" / "mouse_human_orthologs.tsv"
SYMBOL_MAP_PATH = BASE_DIR / "processed" / "ensembl_symbol_map.csv"

sys.path.insert(0, str(BASE_DIR / "v4" / "scripts"))
from v4_utils import (
    TISSUE_MISSIONS, load_metadata, load_gene_features,
    align_features_with_meta, encode_labels, get_folds,
)

ALL_TISSUES = list(TISSUE_MISSIONS.keys())

# da Silveira Cell 2020 — 67 spaceflight consensus genes (human symbols, uppercased)
CELL2020_GENES = {
    "APOD", "ANGPTL4", "BCL6", "BNIP3", "BNIP3L", "C1QA", "C1QB", "C1QC",
    "CASP4", "CCL2", "CCL7", "CD14", "CD163", "CD68", "CDKN1A", "CEBPB",
    "CFB", "CTSS", "CYP1B1", "DDIT4", "DNAJB1", "EGR1", "EIF2AK2",
    "ERRFI1", "F13A1", "FKBP5", "FOS", "GADD45B", "GBP2", "GDF15",
    "HMOX1", "HP", "HSPA1A", "HSPA1B", "IDO1", "IFIT1", "IFIT3",
    "IFITM3", "IL1B", "IL6", "IRF7", "ISG15", "JUN", "KLF4", "LCN2",
    "LGALS3", "LYZ", "MMP9", "MT1A", "MT2A", "MX1", "NFKBIA", "NR4A1",
    "OAS2", "PLIN2", "PNRC1", "PPP1R15A", "S100A8", "S100A9", "SAA1",
    "SERPINA3", "SOCS3", "SPP1", "STAT1", "TIMP1", "TLR2", "TXNIP",
}


def load_orthologs():
    """Load mouse→human ortholog mapping."""
    df = pd.read_csv(ORTHO_PATH, sep="\t")
    return dict(zip(df["mouse_symbol"], df["human_symbol"]))


def load_symbol_map():
    """Load Ensembl → Symbol mapping."""
    df = pd.read_csv(SYMBOL_MAP_PATH)
    return dict(zip(df["ENSEMBL"], df["SYMBOL"]))


def collect_shap_genes():
    """Collect SHAP consensus + tissue-specific genes."""
    consensus_genes = {}  # gene → {n_tissues, tissues, sources}
    tissue_specific = {}  # gene → {tissues}

    # 1. SHAP consensus
    path = V4_EVAL_DIR / "SHAP_consensus.json"
    if path.exists():
        with open(path) as f:
            data = json.load(f)
        for entry in data.get("consensus_genes", []):
            g = entry["gene"]
            consensus_genes[g] = {
                "n_tissues": entry.get("n_tissues", 0),
                "n_methods": entry.get("n_methods", 0),
                "tissues": entry.get("tissues", []),
                "best_rank": entry.get("best_rank", 999),
            }
        for tissue, entries in data.get("tissue_specific_genes", {}).items():
            for entry in entries:
                g = entry["gene"]
                if g not in tissue_specific:
                    tissue_specific[g] = {"tissues": set()}
                tissue_specific[g]["tissues"].add(tissue)

    # 2. SHAP per-tissue top-100 (union across methods)
    shap_top100 = {}  # gene → set of tissues
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
            for g, _ in ranked[:100]:
                if g not in shap_top100:
                    shap_top100[g] = set()
                shap_top100[g].add(tissue)
        except (json.JSONDecodeError, KeyError):
            continue

    return consensus_genes, tissue_specific, shap_top100


def collect_wgcna_hub_genes():
    """Collect WGCNA hub genes and condition-associated modules."""
    hub_genes = {}  # gene → {tissues, modules, condition_associated}
    path = V4_EVAL_DIR / "WGCNA_hub_genes.json"
    if not path.exists():
        return hub_genes

    with open(path) as f:
        data = json.load(f)

    for tissue, tdata in data.get("tissues", {}).items():
        for module, mdata in tdata.get("modules", {}).items():
            is_sig = bool(mdata.get("significant_trait_associations"))
            for g in mdata.get("hub_genes", []):
                if g not in hub_genes:
                    hub_genes[g] = {"tissues": set(), "modules": set(),
                                    "condition_associated": False}
                hub_genes[g]["tissues"].add(tissue)
                hub_genes[g]["modules"].add(f"{tissue}:{module}")
                if is_sig:
                    hub_genes[g]["condition_associated"] = True

    return hub_genes


def collect_ppi_hubs():
    """Collect PPI hub genes from STRING enrichment."""
    ppi_hubs = set()
    path = V4_EVAL_DIR / "PPI_enrichment.json"
    if not path.exists():
        return ppi_hubs

    with open(path) as f:
        data = json.load(f)

    for tissue, tdata in data.get("tissues", {}).items():
        hub_network = tdata.get("hub_genes_network", {})
        for g in hub_network:
            ppi_hubs.add(g)

    return ppi_hubs


def collect_cross_organ_lr():
    """Collect genes involved in cross-organ L-R signaling (from v5 Phase 2)."""
    lr_genes = set()
    path = V5_EVAL_DIR / "cross_organ_signaling.json"
    if not path.exists():
        return lr_genes

    with open(path) as f:
        data = json.load(f)

    for key, pair_data in data.get("tissue_pairs", {}).items():
        for pair in pair_data.get("active_pairs", []):
            lr_genes.add(pair.get("ligand", ""))
            lr_genes.add(pair.get("receptor", ""))

    lr_genes.discard("")
    return lr_genes


def collect_drug_targets():
    """Collect druggable genes from v5 Phase 4."""
    fda_targets = set()
    clinical_targets = set()

    path = V5_EVAL_DIR / "drug_targets.json"
    if not path.exists():
        return fda_targets, clinical_targets

    with open(path) as f:
        data = json.load(f)

    tiers = data.get("tiers", {})
    for entry in tiers.get("tier1_approved", []):
        fda_targets.add(entry.get("gene", ""))
    for entry in tiers.get("tier2_clinical", []):
        clinical_targets.add(entry.get("gene", ""))

    fda_targets.discard("")
    clinical_targets.discard("")
    return fda_targets, clinical_targets


def collect_drosophila_orthologs():
    """Collect genes with significant Drosophila spaceflight response (from v3 E4)."""
    sig_dros_genes = set()

    # Try E4 multispecies NES results
    path = BASE_DIR / "v3" / "evaluation" / "E4_multispecies_nes.json"
    if not path.exists():
        return sig_dros_genes

    try:
        with open(path) as f:
            data = json.load(f)

        # Look for cross-species correlations or NES data
        dros_data = data.get("drosophila", data.get("cross_species", {}))
        if isinstance(dros_data, dict):
            for gene, info in dros_data.items():
                if isinstance(info, dict) and abs(info.get("nes", 0)) > 1.5:
                    sig_dros_genes.add(gene)
                elif isinstance(info, (int, float)) and abs(info) > 1.5:
                    sig_dros_genes.add(gene)
    except (json.JSONDecodeError, KeyError):
        pass

    return sig_dros_genes


def score_genes(consensus_genes, tissue_specific, shap_top100,
                wgcna_hubs, ppi_hubs, lr_genes,
                fda_targets, clinical_targets, cell2020_overlap,
                dros_genes, ortho_map):
    """Score each gene using 6 evidence sources (max 12 points).

    Scoring system (revised to reduce druggability bias):
      Source 1 — SHAP importance:        0-3 pts (consensus=3, tissue-specific=2, top-100=1)
      Source 2 — WGCNA hub:              0-2 pts (condition-associated=2, any hub=1)
      Source 3 — Multi-tissue breadth:   0-2 pts (≥4 tissues=2, ≥2 tissues=1)
      Source 4 — Network centrality:     0-1 pt  (PPI hub OR L-R pair)
      Source 5 — Druggability:           0-1 pt  (FDA-approved or clinical target)
      Source 6 — Literature validation:  0-2 pts (Cell 2020=2)
      Source 7 — Cross-species:          0-1 pt  (significant Drosophila ortholog)
      Max possible: 3+2+2+1+1+2+1 = 12
    """
    # Collect all candidate genes
    all_genes = set()
    all_genes.update(consensus_genes.keys())
    all_genes.update(tissue_specific.keys())
    all_genes.update(shap_top100.keys())
    all_genes.update(wgcna_hubs.keys())

    # Map mouse→human for drug/literature lookup
    mouse_to_human = {}
    for g in all_genes:
        h = ortho_map.get(g)
        if h:
            mouse_to_human[g] = h

    scored = {}
    for gene in all_genes:
        score = 0
        sources = {}
        human_sym = mouse_to_human.get(gene, gene.upper())

        # Collect tissue info first (needed for multi-tissue scoring)
        tissues = set()
        if gene in consensus_genes:
            tissues.update(consensus_genes[gene].get("tissues", []))
        if gene in tissue_specific:
            tissues.update(tissue_specific[gene].get("tissues", set()))
        if gene in shap_top100:
            tissues.update(shap_top100[gene])
        if gene in wgcna_hubs:
            tissues.update(wgcna_hubs[gene]["tissues"])

        # Source 1: SHAP importance (0-3 points) — primary evidence
        if gene in consensus_genes:
            score += 3
            sources["shap"] = 3
        elif gene in tissue_specific:
            score += 2
            sources["shap"] = 2
        elif gene in shap_top100:
            score += 1
            sources["shap"] = 1

        # Source 2: WGCNA hub gene (0-2 points)
        if gene in wgcna_hubs:
            if wgcna_hubs[gene]["condition_associated"]:
                score += 2
                sources["wgcna_hub"] = 2
            else:
                score += 1
                sources["wgcna_hub"] = 1

        # Source 3: Multi-tissue breadth (0-2 points) — NEW
        # Rewards genes detected across multiple tissues (robust signal)
        n_tiss = len(tissues)
        if n_tiss >= 4:
            score += 2
            sources["multi_tissue"] = 2
        elif n_tiss >= 2:
            score += 1
            sources["multi_tissue"] = 1

        # Source 4: Network centrality (0-1 point) — reduced from 0-2
        if gene in ppi_hubs or gene in lr_genes:
            score += 1
            sources["network"] = 1

        # Source 5: Druggability (0-1 point) — reduced from 0-2
        if human_sym in fda_targets or human_sym in clinical_targets:
            score += 1
            sources["druggable"] = 1

        # Source 6: Literature validation (0-2 points)
        if human_sym in cell2020_overlap:
            score += 2
            sources["literature"] = 2

        # Source 7: Cross-species conservation (0-1 point)
        if human_sym in dros_genes or gene in dros_genes:
            score += 1
            sources["cross_species"] = 1

        scored[gene] = {
            "score": score,
            "sources": sources,
            "tissues": sorted(tissues),
            "n_tissues": n_tiss,
            "human_symbol": human_sym,
        }

    return scored


def validate_panel(panel_genes, sym_map_ens_to_sym):
    """Validate panel by running PCA-LR on panel genes only using v4 folds."""
    # Build symbol → ensembl reverse map
    sym_to_ens = {}
    for ens, sym in sym_map_ens_to_sym.items():
        if sym not in sym_to_ens:
            sym_to_ens[sym] = ens

    results_per_tissue = {}

    for tissue in ALL_TISSUES:
        try:
            folds = get_folds(tissue)
        except Exception as e:
            print(f"  [WARN] {tissue}: get_folds failed — {e}")
            continue

        # Get gene names from first fold
        gene_names = folds[0].get("gene_names", [])
        if not gene_names:
            continue

        # Find panel gene indices (match by symbol)
        ens_to_sym = {g: sym_map_ens_to_sym.get(g, g) for g in gene_names}
        panel_indices = []
        panel_matched = []
        for i, ens_id in enumerate(gene_names):
            sym = ens_to_sym.get(ens_id, "")
            if sym in panel_genes:
                panel_indices.append(i)
                panel_matched.append(sym)

        if len(panel_indices) < 5:
            print(f"  [WARN] {tissue}: only {len(panel_indices)} panel genes found, skipping")
            results_per_tissue[tissue] = {
                "auroc": None, "n_panel_genes_found": len(panel_indices),
                "note": "insufficient panel genes"
            }
            continue

        fold_aurocs = []
        for fold in folds:
            X_train = fold["train_X"][:, panel_indices]
            y_train = fold["train_y"]
            X_test = fold["test_X"][:, panel_indices]
            y_test = fold["test_y"]

            if len(np.unique(y_test)) < 2:
                continue

            # Simple PCA-LR (adapt n_components for small panel)
            n_comp = min(5, X_train.shape[1] - 1, X_train.shape[0] - 1)
            if n_comp < 1:
                continue

            model = Pipeline([
                ("scaler", StandardScaler()),
                ("pca", PCA(n_components=n_comp)),
                ("lr", LogisticRegression(max_iter=1000, random_state=42)),
            ])

            try:
                model.fit(X_train, y_train)
                y_score = model.predict_proba(X_test)[:, 1]
                auroc = roc_auc_score(y_test, y_score)
                fold_aurocs.append(auroc)
            except Exception:
                continue

        if fold_aurocs:
            mean_auroc = float(np.mean(fold_aurocs))
            results_per_tissue[tissue] = {
                "auroc": round(mean_auroc, 4),
                "std_auroc": round(float(np.std(fold_aurocs)), 4),
                "n_folds": len(fold_aurocs),
                "n_panel_genes_found": len(panel_indices),
                "panel_genes_matched": panel_matched[:20],
            }
            print(f"  {tissue}: AUROC={mean_auroc:.3f} ({len(fold_aurocs)} folds, "
                  f"{len(panel_indices)} panel genes)")
        else:
            results_per_tissue[tissue] = {
                "auroc": None, "note": "no valid folds"
            }

    return results_per_tissue


def main():
    V5_EVAL_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Phase 5: Consensus Biomarker Panel")
    print("=" * 60)

    # Load mappings
    ortho_map = load_orthologs()
    sym_map = load_symbol_map()

    # Collect all evidence sources
    print("\n1. Collecting SHAP evidence...")
    consensus_genes, tissue_specific, shap_top100 = collect_shap_genes()
    print(f"  Consensus: {len(consensus_genes)}, "
          f"Tissue-specific: {len(tissue_specific)}, "
          f"Any top-100: {len(shap_top100)}")

    print("\n2. Collecting WGCNA hub genes...")
    wgcna_hubs = collect_wgcna_hub_genes()
    n_cond = sum(1 for v in wgcna_hubs.values() if v["condition_associated"])
    print(f"  Hub genes: {len(wgcna_hubs)} ({n_cond} in condition-associated modules)")

    print("\n3. Collecting PPI hub genes...")
    ppi_hubs = collect_ppi_hubs()
    print(f"  PPI hubs: {len(ppi_hubs)}")

    print("\n4. Collecting cross-organ L-R genes...")
    lr_genes = collect_cross_organ_lr()
    print(f"  L-R genes: {len(lr_genes)}")

    print("\n5. Collecting drug targets...")
    fda_targets, clinical_targets = collect_drug_targets()
    print(f"  FDA-approved targets: {len(fda_targets)}, "
          f"Clinical targets: {len(clinical_targets)}")

    print("\n6. Collecting cross-species (Drosophila) evidence...")
    dros_genes = collect_drosophila_orthologs()
    print(f"  Significant Drosophila orthologs: {len(dros_genes)}")

    # Cell 2020 overlap (mouse→human → check against Cell 2020 list)
    cell2020_set = CELL2020_GENES
    print(f"\n7. Cell 2020 reference: {len(cell2020_set)} genes")

    # Score all genes
    print("\n8. Scoring genes (7 sources, max 12 points)...")
    scored = score_genes(
        consensus_genes, tissue_specific, shap_top100,
        wgcna_hubs, ppi_hubs, lr_genes,
        fda_targets, clinical_targets, cell2020_set,
        dros_genes, ortho_map,
    )
    print(f"  Total candidate genes scored: {len(scored)}")

    # Rank and select top-20 panel
    ranked = sorted(scored.items(), key=lambda x: (-x[1]["score"], -x[1].get("n_tissues", 0)))
    panel_size = 20
    panel = ranked[:panel_size]

    print(f"\n  Top-{panel_size} panel:")
    for i, (gene, info) in enumerate(panel):
        src_str = ", ".join(f"{k}={v}" for k, v in info["sources"].items())
        print(f"    {i+1}. {gene} (score={info['score']}, "
              f"tissues={info['n_tissues']}, {src_str})")

    # Score distribution
    scores = [v["score"] for v in scored.values()]
    print(f"\n  Score distribution: min={min(scores)}, max={max(scores)}, "
          f"mean={np.mean(scores):.1f}, median={np.median(scores):.0f}")
    for threshold in [10, 8, 6, 4, 2]:
        n = sum(1 for s in scores if s >= threshold)
        print(f"    Score ≥ {threshold}: {n} genes")

    # Validate panel on v4 folds
    print("\n9. Validating panel via PCA-LR on v4 folds...")
    panel_gene_set = {gene for gene, _ in panel}
    panel_validation = validate_panel(panel_gene_set, sym_map)

    # Build drug connections for panel genes
    drug_connections = {}
    drug_path = V5_EVAL_DIR / "drug_targets.json"
    if drug_path.exists():
        with open(drug_path) as f:
            drug_data = json.load(f)
        tiers = drug_data.get("tiers", {})
        for tier_name in ["tier1_approved", "tier2_clinical"]:
            for entry in tiers.get(tier_name, []):
                g = entry.get("gene", "")
                if g in panel_gene_set or g in {scored[pg]["human_symbol"]
                                                  for pg, _ in panel if pg in scored}:
                    if g not in drug_connections:
                        drug_connections[g] = []
                    drug_connections[g].append({
                        "drug": entry.get("drug", ""),
                        "tier": tier_name,
                        "source": entry.get("source", ""),
                    })

    # Save results
    output = {
        "panel_size": panel_size,
        "n_candidates_scored": len(scored),
        "max_possible_score": 12,
        "panel": [
            {
                "rank": i + 1,
                "gene": gene,
                "human_symbol": info["human_symbol"],
                "score": info["score"],
                "sources": info["sources"],
                "tissues": info["tissues"],
                "n_tissues": info["n_tissues"],
                "drugs": drug_connections.get(info["human_symbol"], []),
            }
            for i, (gene, info) in enumerate(panel)
        ],
        "score_distribution": {
            "min": int(min(scores)),
            "max": int(max(scores)),
            "mean": round(float(np.mean(scores)), 2),
            "median": int(np.median(scores)),
            "n_ge10": sum(1 for s in scores if s >= 10),
            "n_ge8": sum(1 for s in scores if s >= 8),
            "n_ge6": sum(1 for s in scores if s >= 6),
            "n_ge4": sum(1 for s in scores if s >= 4),
        },
        "panel_validation": panel_validation,
        "evidence_sources": {
            "shap_consensus": len(consensus_genes),
            "shap_tissue_specific": len(tissue_specific),
            "shap_any_top100": len(shap_top100),
            "wgcna_hub_total": len(wgcna_hubs),
            "wgcna_hub_condition": n_cond,
            "ppi_hubs": len(ppi_hubs),
            "lr_genes": len(lr_genes),
            "fda_drug_targets": len(fda_targets),
            "clinical_drug_targets": len(clinical_targets),
            "drosophila_orthologs": len(dros_genes),
            "cell2020_reference": len(cell2020_set),
        },
    }

    out_path = V5_EVAL_DIR / "consensus_biomarker_panel.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved: {out_path}")
    print("Done!")


if __name__ == "__main__":
    main()
