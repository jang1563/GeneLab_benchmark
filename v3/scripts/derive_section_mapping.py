#!/usr/bin/env python3
"""
derive_section_mapping.py — Derive section→animal→condition mapping
from OSD-352 Visium data without Seurat metadata.

Strategy:
  1. Load pseudo-bulk for all 12 sections (from process_osd352_visium.py output)
  2. Hierarchical clustering to identify 6 same-animal pairs
  3. Cross-reference with companion bulk RNA-seq to assign conditions
     (map gene symbols → ENSMUSG via DEG file, correlate with bulk)
  4. Output: section_mapping.json

This script is needed because the Seurat RDS file on OSDR is corrupted.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
from scipy import stats

BASE_DIR = Path(__file__).resolve().parent.parent  # v3/
DATA_DIR = BASE_DIR / "data" / "spatial" / "osd352"
PROCESSED_DIR = DATA_DIR / "processed"


def load_pseudobulk():
    """Load pseudo-bulk matrix and section names."""
    X = np.load(PROCESSED_DIR / "pseudobulk_matrix.npy")
    gene_names = np.load(PROCESSED_DIR / "gene_names.npy", allow_pickle=True)
    with open(PROCESSED_DIR / "section_names.json") as f:
        section_names = json.load(f)
    return X, gene_names, section_names


def pair_sections_by_clustering(X, section_names, n_animals=6):
    """Pair sections into animals using hierarchical clustering.

    Same-animal sections should be the most correlated (technical replicates).

    Args:
        X: (12, n_genes) pseudo-bulk matrix
        section_names: list of 12 section names
        n_animals: expected number of animals (6)

    Returns:
        dict: animal_id → list of section names
    """
    # Pairwise correlation distance
    corr_dist = pdist(X, metric="correlation")  # 1 - pearson_r

    # Hierarchical clustering (complete linkage)
    Z = linkage(corr_dist, method="complete")

    # Cut tree to get n_animals clusters
    labels = fcluster(Z, t=n_animals, criterion="maxclust")

    # Group sections by cluster
    animals = {}
    for i, sec in enumerate(section_names):
        cluster_id = labels[i]
        if cluster_id not in animals:
            animals[cluster_id] = []
        animals[cluster_id].append(sec)

    # Verify: each animal should have exactly 2 sections
    print("\nSection pairs (from clustering):")
    for aid, secs in sorted(animals.items()):
        print(f"  Animal {aid}: {secs}")

    # Check for expected pairing
    pair_sizes = [len(secs) for secs in animals.values()]
    if all(s == 2 for s in pair_sizes):
        print("  ✓ All animals have exactly 2 sections")
    else:
        print(f"  ⚠ Unexpected pair sizes: {pair_sizes}")
        print("  Falling back to slide-based pairing")
        return pair_sections_by_slide(section_names)

    return animals


def pair_sections_by_slide(section_names):
    """Fallback: pair sections by slide (adjacent capture areas).

    Assumption: on each slide, A1+B1 = one animal, C1+D1 = another.
    """
    animals = {}
    aid = 1
    slides = sorted(set(s.rsplit("_", 1)[0] for s in section_names))

    for slide in slides:
        slide_sections = sorted(s for s in section_names
                                if s.startswith(slide))
        # Pair A1+B1 and C1+D1
        ab = [s for s in slide_sections if s.endswith(("A1", "B1"))]
        cd = [s for s in slide_sections if s.endswith(("C1", "D1"))]
        if len(ab) == 2:
            animals[aid] = ab
            aid += 1
        if len(cd) == 2:
            animals[aid] = cd
            aid += 1

    return animals


def assign_conditions_via_bulk(X, gene_names, section_names, animals):
    """Assign FLT/GC conditions by correlating with bulk RNA-seq.

    Steps:
    1. Load bulk normalized counts (ENSMUSG → gene_symbol via DEG file)
    2. Find common genes between Visium (symbols) and bulk (mapped to symbols)
    3. Correlate each animal's pseudo-bulk with each bulk sample
    4. Hungarian matching → assign animal IDs and conditions
    """
    # Load bulk data
    bulk_counts_path = DATA_DIR / "GLDS-352_rna_seq_Normalized_Counts.csv"
    bulk_deg_path = DATA_DIR / "GLDS-352_rna_seq_differential_expression.csv"

    if not bulk_counts_path.exists():
        print("  ERROR: Bulk counts not found")
        return None

    bulk_df = pd.read_csv(bulk_counts_path, index_col=0)  # ENSMUSG × samples
    print(f"\n  Bulk: {bulk_df.shape} (genes × samples)")

    # Get ENSMUSG → Symbol mapping from DEG file
    if bulk_deg_path.exists():
        deg = pd.read_csv(bulk_deg_path, index_col=0)
        # Find SYMBOL column
        sym_col = None
        for col in deg.columns:
            if col.upper() in ("SYMBOL", "GENE_SYMBOL", "GENENAME"):
                sym_col = col
                break
        if sym_col is None:
            # Try to find by content
            for col in deg.columns:
                sample = deg[col].dropna().astype(str).iloc[:5]
                if all(s.isalpha() or "-" in s for s in sample):
                    sym_col = col
                    break

        if sym_col:
            ensmusg_to_symbol = dict(zip(deg.index, deg[sym_col]))
            print(f"  Gene mapping: {len(ensmusg_to_symbol)} ENSMUSG → {sym_col}")
        else:
            print(f"  WARNING: No SYMBOL column in DEG. Columns: {list(deg.columns)[:10]}")
            return None
    else:
        print("  ERROR: DEG file not found for gene mapping")
        return None

    # Map bulk genes to symbols
    bulk_symbols = {}
    for ensmusg in bulk_df.index:
        sym = ensmusg_to_symbol.get(ensmusg)
        if sym and isinstance(sym, str) and sym != "nan":
            bulk_symbols[ensmusg] = sym

    # Find common genes (Visium symbol ∩ bulk symbol)
    visium_genes = set(gene_names)
    bulk_mapped_genes = set(bulk_symbols.values())
    common_symbols = sorted(visium_genes & bulk_mapped_genes)
    print(f"  Common genes: {len(common_symbols)} "
          f"(Visium: {len(visium_genes)}, Bulk mapped: {len(bulk_mapped_genes)})")

    if len(common_symbols) < 500:
        print("  WARNING: Too few common genes for reliable correlation")
        return None

    # Build aligned matrices
    # Visium: section × genes (indexed by symbol)
    visium_gene_idx = {g: i for i, g in enumerate(gene_names)}
    common_visium_idx = [visium_gene_idx[g] for g in common_symbols]
    X_visium = X[:, common_visium_idx]

    # Build animal-level pseudo-bulk
    animal_ids = sorted(animals.keys())
    X_animal = np.zeros((len(animal_ids), len(common_symbols)), dtype=np.float32)
    for i, aid in enumerate(animal_ids):
        sec_indices = [section_names.index(s) for s in animals[aid]]
        X_animal[i] = X_visium[sec_indices].mean(axis=0)

    # Bulk: samples × genes (indexed by symbol, log1p if needed)
    symbol_to_ensmusg = {}
    for e, s in bulk_symbols.items():
        if s in set(common_symbols):
            symbol_to_ensmusg.setdefault(s, e)

    bulk_common_ensmusg = [symbol_to_ensmusg[s] for s in common_symbols]
    X_bulk = bulk_df.loc[bulk_common_ensmusg].values.T.astype(np.float32)
    # log1p if values are large
    if X_bulk.max() > 100:
        X_bulk = np.log1p(X_bulk)

    bulk_sample_names = list(bulk_df.columns)

    # Correlate each Visium animal with each bulk sample
    print("\n  Correlation matrix (Visium animal × Bulk sample):")
    corr_matrix = np.zeros((len(animal_ids), len(bulk_sample_names)))
    for i in range(len(animal_ids)):
        for j in range(len(bulk_sample_names)):
            r, _ = stats.pearsonr(X_animal[i], X_bulk[j])
            corr_matrix[i, j] = r

    # Print correlation matrix
    print(f"  {'':>10}", end="")
    for bs in bulk_sample_names:
        short = bs.split("_")[-1]  # F1, F2, F7, G7, G8, G9
        print(f"  {short:>6}", end="")
    print()
    for i, aid in enumerate(animal_ids):
        secs = ",".join(s.replace("Sample_", "S") for s in animals[aid])
        print(f"  {secs:>10}", end="")
        for j in range(len(bulk_sample_names)):
            print(f"  {corr_matrix[i, j]:.3f}", end="")
        print()

    # Assign: each Visium animal → best-matching bulk sample
    # Greedy matching: iteratively pick the highest-correlation pair
    used_bulk = set()
    assigned_aids = set()
    assignments = {}  # animal_cluster_id → bulk_sample_name

    for _ in range(len(animal_ids)):
        best_r = -1
        best_i, best_j = -1, -1
        for i in range(len(animal_ids)):
            if animal_ids[i] in assigned_aids:
                continue
            for j in range(len(bulk_sample_names)):
                if j in used_bulk:
                    continue
                if corr_matrix[i, j] > best_r:
                    best_r = corr_matrix[i, j]
                    best_i, best_j = i, j
        if best_i >= 0:
            assignments[animal_ids[best_i]] = bulk_sample_names[best_j]
            assigned_aids.add(animal_ids[best_i])
            used_bulk.add(best_j)

    # Build final mapping
    mapping = {}
    print("\n  Final assignments:")
    for aid in animal_ids:
        bulk_name = assignments.get(aid, "unknown")
        parts = bulk_name.split("_")
        condition = "FLT" if "FLT" in parts else "GC" if "GC" in parts else "unknown"
        animal_id = parts[-1] if parts else "unknown"

        for sec in animals[aid]:
            mapping[sec] = {
                "animal": animal_id,
                "condition": condition,
                "bulk_sample": bulk_name,
                "cluster_id": int(aid)
            }
        r_val = corr_matrix[animal_ids.index(aid),
                             bulk_sample_names.index(bulk_name)]
        print(f"  Cluster {aid} ({animals[aid]}) → {bulk_name} "
              f"({condition}, r={r_val:.3f})")

    return mapping


def validate_mapping(mapping, section_names):
    """Validate the derived mapping."""
    print("\n=== Validation ===")

    # Check all sections mapped
    mapped = set(mapping.keys())
    expected = set(section_names)
    if mapped == expected:
        print("  ✓ All 12 sections mapped")
    else:
        print(f"  ✗ Missing: {expected - mapped}")

    # Check balanced design
    conditions = [mapping[s]["condition"] for s in section_names]
    n_flt = conditions.count("FLT")
    n_gc = conditions.count("GC")
    if n_flt == 6 and n_gc == 6:
        print(f"  ✓ Balanced: {n_flt} FLT, {n_gc} GC sections")
    else:
        print(f"  ⚠ Imbalanced: {n_flt} FLT, {n_gc} GC")

    # Check same-animal sections have same condition
    animals = {}
    for sec, info in mapping.items():
        a = info["animal"]
        if a not in animals:
            animals[a] = set()
        animals[a].add(info["condition"])
    for a, conds in animals.items():
        if len(conds) > 1:
            print(f"  ✗ Animal {a} has mixed conditions: {conds}")
        else:
            print(f"  ✓ Animal {a}: {list(conds)[0]}")

    # Check sections per animal
    sec_per_animal = {}
    for sec, info in mapping.items():
        a = info["animal"]
        sec_per_animal.setdefault(a, []).append(sec)
    for a, secs in sorted(sec_per_animal.items()):
        print(f"    {a}: {secs}")

    # Check slide distribution
    print("\n  Slide distribution:")
    for slide_num in ["158", "159", "304"]:
        slide_secs = [s for s in section_names if f"_{slide_num}_" in s]
        conds = [mapping[s]["condition"] for s in slide_secs]
        animals_on_slide = set(mapping[s]["animal"] for s in slide_secs)
        print(f"    Slide {slide_num}: {conds} — animals: {animals_on_slide}")


def main():
    print("=" * 60)
    print("OSD-352 — Deriving Section → Animal → Condition Mapping")
    print("=" * 60)

    # Step 1: Load pseudo-bulk
    print("\n[1] Loading pseudo-bulk data...")
    X, gene_names, section_names = load_pseudobulk()
    print(f"  Matrix: {X.shape} (sections × genes)")
    print(f"  Sections: {section_names}")

    # Step 2: Pair sections into animals
    print("\n[2] Clustering sections into animal pairs...")
    animals = pair_sections_by_clustering(X, section_names)

    # Step 3: Assign conditions via bulk correlation
    print("\n[3] Assigning conditions via bulk RNA-seq correlation...")
    mapping = assign_conditions_via_bulk(X, gene_names, section_names, animals)

    if mapping is None:
        print("\n  FAILED to derive mapping. Creating placeholder.")
        # Fallback: use slide-based pairing with unknown conditions
        mapping = {}
        for aid, secs in animals.items():
            for sec in secs:
                mapping[sec] = {
                    "animal": f"animal_{aid}",
                    "condition": "unknown",
                    "cluster_id": int(aid)
                }

    # Step 4: Validate
    validate_mapping(mapping, section_names)

    # Step 5: Save
    output_path = DATA_DIR / "section_mapping.json"
    with open(output_path, "w") as f:
        json.dump(mapping, f, indent=2)
    print(f"\n[DONE] Saved: {output_path}")

    return mapping


if __name__ == "__main__":
    main()
