#!/usr/bin/env python3
"""
build_ortholog_db.py — Download cross-species orthologs from Ensembl BioMart.

Downloads ortholog mappings:
  - Mouse ↔ Drosophila (1:1 orthologs)
  - Mouse ↔ Arabidopsis (1:1 orthologs)
  - Mouse ↔ C. elegans (1:1 orthologs, if available)

Uses Ensembl REST API (BioMart XML query) for reliability.

Usage:
    python3 v3/scripts/build_ortholog_db.py

Output:
    v3/processed/orthologs/mouse_drosophila.tsv
    v3/processed/orthologs/mouse_arabidopsis.tsv
    v3/processed/orthologs/mouse_celegans.tsv
    v3/processed/orthologs/ortholog_summary.json
"""

import csv
import json
import sys
import time
from io import StringIO
from pathlib import Path
from urllib.request import urlopen, Request
from urllib.error import URLError, HTTPError
from urllib.parse import quote

BASE_DIR = Path(__file__).resolve().parent.parent.parent  # GeneLab_benchmark/
OUTPUT_DIR = BASE_DIR / "v3" / "processed" / "orthologs"

# Ensembl BioMart endpoints (primary + mirrors)
BIOMART_URLS = [
    "https://www.ensembl.org/biomart/martservice",
    "https://useast.ensembl.org/biomart/martservice",
    "https://asia.ensembl.org/biomart/martservice",
]

# Species configuration
# For BioMart: dataset = "{species}_gene_ensembl"
# Ortholog attributes: "{prefix}_homolog_ensembl_gene", "{prefix}_homolog_associated_gene_name"
SPECIES = {
    "drosophila": {
        "prefix": "dmelanogaster",
        "label": "Drosophila melanogaster",
        "output_file": "mouse_drosophila.tsv",
        "biomart": "ensembl",  # main Ensembl
    },
    "celegans": {
        "prefix": "celegans",
        "label": "Caenorhabditis elegans",
        "output_file": "mouse_celegans.tsv",
        "biomart": "ensembl",
    },
}
# NOTE: Arabidopsis is NOT included because:
# 1. Ensembl Plants BioMart does not have mouse homolog attributes
# 2. Phase 1 E4 uses KEGG pathway-level NES concordance (not gene-level orthologs)
# 3. KEGG gene sets exist per-species (mmu, dme, ath, cel) → pathway-level comparison

# BioMart XML template for mouse orthologs
# Query mouse genes and their orthologs in target species
BIOMART_XML_TEMPLATE = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6">
  <Dataset name="mmusculus_gene_ensembl" interface="default">
    <Attribute name="ensembl_gene_id"/>
    <Attribute name="external_gene_name"/>
    <Attribute name="{prefix}_homolog_ensembl_gene"/>
    <Attribute name="{prefix}_homolog_associated_gene_name"/>
    <Attribute name="{prefix}_homolog_orthology_type"/>
    <Attribute name="{prefix}_homolog_perc_id"/>
    <Attribute name="{prefix}_homolog_perc_id_r1"/>
  </Dataset>
</Query>"""



def fetch_biomart(xml_query, urls=None, timeout=300):
    """Fetch data from BioMart using XML query, trying mirrors on failure."""
    if urls is None:
        urls = BIOMART_URLS
    encoded = quote(xml_query.strip())
    last_error = None
    for url in urls:
        full_url = f"{url}?query={encoded}"
        req = Request(full_url, headers={"User-Agent": "GeneLabBench/3.0"})
        try:
            print(f"    Trying {url.split('//')[1].split('/')[0]}...")
            with urlopen(req, timeout=timeout) as resp:
                data = resp.read().decode("utf-8")
                if "ERROR" in data[:200] or "Query ERROR" in data[:200]:
                    last_error = f"BioMart error: {data[:500]}"
                    continue
                return data, None
        except (URLError, HTTPError) as e:
            last_error = str(e)
            continue
    return None, last_error


def parse_biomart_tsv(data):
    """Parse BioMart TSV output into ortholog records."""
    records = []
    reader = csv.reader(StringIO(data), delimiter="\t")
    header = next(reader, None)
    if not header:
        return records

    for row in reader:
        if len(row) < 7:
            continue
        mouse_gene_id = row[0]
        mouse_gene_name = row[1]
        target_gene_id = row[2]
        target_gene_name = row[3]
        orth_type = row[4]
        pct_id = row[5]
        pct_id_r1 = row[6]

        if not target_gene_id:
            continue
        records.append({
            "mouse_ensembl": mouse_gene_id,
            "mouse_symbol": mouse_gene_name,
            "target_ensembl": target_gene_id,
            "target_symbol": target_gene_name,
            "orthology_type": orth_type,
            "pct_id_mouse_to_target": pct_id,
            "pct_id_target_to_mouse": pct_id_r1,
        })

    return records


def save_tsv(records, output_path, species_label):
    """Save ortholog records as TSV."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "mouse_ensembl", "mouse_symbol",
            f"{species_label}_ensembl", f"{species_label}_symbol",
            "orthology_type", "pct_id_mouse_to_target", "pct_id_target_to_mouse",
        ])
        for r in records:
            writer.writerow([
                r["mouse_ensembl"], r["mouse_symbol"],
                r["target_ensembl"], r["target_symbol"],
                r["orthology_type"],
                r["pct_id_mouse_to_target"], r["pct_id_target_to_mouse"],
            ])


def main():
    print("=" * 70)
    print("GeneLabBench v3 — Ortholog Database Builder")
    print("=" * 70)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    summary = {}

    for species_key, cfg in SPECIES.items():
        label = cfg["label"]
        output_file = OUTPUT_DIR / cfg["output_file"]
        print(f"\n--- {label} ---")

        xml = BIOMART_XML_TEMPLATE.format(prefix=cfg["prefix"])
        print(f"  Querying Ensembl BioMart (with mirror fallback)...")
        data, error = fetch_biomart(xml)

        if error:
            print(f"  ERROR: {error}")
            summary[species_key] = {
                "species": label,
                "status": "error",
                "error": error,
            }
            continue

        records = parse_biomart_tsv(data)
        total = len(records)

        # Count by orthology type
        type_counts = {}
        for r in records:
            otype = r["orthology_type"]
            type_counts[otype] = type_counts.get(otype, 0) + 1

        # Filter to ortholog_one2one for clean mapping
        one2one = [r for r in records if r["orthology_type"] == "ortholog_one2one"]
        n_one2one = len(one2one)

        # Unique mouse genes with orthologs
        unique_mouse = len(set(r["mouse_ensembl"] for r in one2one))
        unique_target = len(set(r["target_ensembl"] for r in one2one))

        print(f"  Total ortholog pairs: {total}")
        print(f"  By type: {type_counts}")
        print(f"  1:1 orthologs: {n_one2one}")
        print(f"  Unique mouse genes: {unique_mouse}")
        print(f"  Unique {label} genes: {unique_target}")

        # Save all records (not just 1:1)
        save_tsv(records, output_file, species_key)
        print(f"  Saved to: {output_file}")

        summary[species_key] = {
            "species": label,
            "status": "ok",
            "total_pairs": total,
            "one2one_pairs": n_one2one,
            "unique_mouse_genes": unique_mouse,
            "unique_target_genes": unique_target,
            "type_counts": type_counts,
            "output_file": str(output_file),
        }

        time.sleep(1)  # Rate limit between queries

    # Save summary
    summary_path = OUTPUT_DIR / "ortholog_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nSummary saved to: {summary_path}")

    # Print final summary
    print("\n" + "=" * 70)
    print("ORTHOLOG SUMMARY")
    print("=" * 70)
    for key, s in summary.items():
        if s["status"] == "ok":
            print(f"  {s['species']}: {s['one2one_pairs']} 1:1 orthologs "
                  f"({s['unique_mouse_genes']} mouse genes)")
        else:
            print(f"  {s['species']}: ERROR — {s.get('error', 'unknown')}")


if __name__ == "__main__":
    main()
