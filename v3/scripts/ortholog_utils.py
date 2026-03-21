#!/usr/bin/env python3
"""
ortholog_utils.py — Cross-species ortholog mapping utilities.

Provides functions for:
  - Loading ortholog TSV files (mouse↔drosophila, mouse↔celegans)
  - Gene ID conversion between species
  - KEGG pathway intersection across species
  - Coverage statistics

Used by Phase 1 (E4/E5 multi-species NES concordance).
"""

import csv
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

BASE_DIR = Path(__file__).resolve().parent.parent.parent  # GeneLab_benchmark/
ORTHOLOG_DIR = BASE_DIR / "v3" / "processed" / "orthologs"

# Existing mouse-human ortholog from v1/v2
MOUSE_HUMAN_ORTHOLOG = BASE_DIR / "data" / "mouse" / "ensembl_mouse_human_orthologs_with_symbols.tsv"


def load_ortholog_map(
    species: str,
    orthology_type: Optional[str] = "ortholog_one2one",
) -> Dict[str, str]:
    """Load ortholog mapping: mouse_ensembl → target_ensembl.

    Args:
        species: "drosophila", "celegans", or "human"
        orthology_type: filter by type (None = all). Default "ortholog_one2one".

    Returns:
        Dict mapping mouse ENSMUSG ID → target Ensembl gene ID.
    """
    if species == "human":
        return _load_mouse_human()

    tsv_path = ORTHOLOG_DIR / f"mouse_{species}.tsv"
    if not tsv_path.exists():
        raise FileNotFoundError(f"Ortholog file not found: {tsv_path}")

    mapping = {}
    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if orthology_type and row["orthology_type"] != orthology_type:
                continue
            mouse_id = row["mouse_ensembl"]
            target_id = row[f"{species}_ensembl"]
            mapping[mouse_id] = target_id
    return mapping


def load_symbol_map(species: str) -> Dict[str, str]:
    """Load mouse_ensembl → target_symbol mapping.

    Useful for converting mouse gene IDs to target species gene symbols
    (e.g., for gene set analysis where gene sets use symbols).
    """
    if species == "human":
        return _load_mouse_human_symbols()

    tsv_path = ORTHOLOG_DIR / f"mouse_{species}.tsv"
    if not tsv_path.exists():
        raise FileNotFoundError(f"Ortholog file not found: {tsv_path}")

    mapping = {}
    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["orthology_type"] != "ortholog_one2one":
                continue
            mouse_id = row["mouse_ensembl"]
            target_sym = row[f"{species}_symbol"]
            if target_sym:
                mapping[mouse_id] = target_sym
    return mapping


def get_ortholog_coverage(
    species: str,
    gene_universe: Set[str],
) -> dict:
    """Calculate ortholog coverage for a gene universe.

    Args:
        species: target species
        gene_universe: set of mouse ENSMUSG gene IDs

    Returns:
        dict with coverage statistics
    """
    ortho_map = load_ortholog_map(species)
    covered = gene_universe & set(ortho_map.keys())
    return {
        "species": species,
        "total_genes": len(gene_universe),
        "genes_with_ortholog": len(covered),
        "coverage_pct": round(100 * len(covered) / max(len(gene_universe), 1), 1),
        "total_orthologs_available": len(ortho_map),
    }


def _load_mouse_human() -> Dict[str, str]:
    """Load existing v1 mouse-human ortholog map."""
    if not MOUSE_HUMAN_ORTHOLOG.exists():
        raise FileNotFoundError(f"Mouse-human ortholog not found: {MOUSE_HUMAN_ORTHOLOG}")
    mapping = {}
    with open(MOUSE_HUMAN_ORTHOLOG) as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)  # skip header
        for row in reader:
            if len(row) >= 2 and row[1]:  # has human ortholog
                mapping[row[0]] = row[1]
    return mapping


def _load_mouse_human_symbols() -> Dict[str, str]:
    """Load mouse_ensembl → human_symbol from v1 extended ortholog file."""
    if not MOUSE_HUMAN_ORTHOLOG.exists():
        raise FileNotFoundError(f"Mouse-human ortholog not found: {MOUSE_HUMAN_ORTHOLOG}")
    mapping = {}
    with open(MOUSE_HUMAN_ORTHOLOG) as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)
        for row in reader:
            if len(row) >= 3 and row[2]:
                mapping[row[0]] = row[2]  # col 2 = Human gene name
    return mapping


def list_available_species() -> List[str]:
    """List species with available ortholog mappings."""
    species = []
    if MOUSE_HUMAN_ORTHOLOG.exists():
        species.append("human")
    for tsv in ORTHOLOG_DIR.glob("mouse_*.tsv"):
        sp = tsv.stem.replace("mouse_", "")
        species.append(sp)
    return sorted(species)
