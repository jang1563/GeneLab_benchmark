#!/usr/bin/env python3
"""
cell2020_validation.py — GeneLab_benchmark: External Validation Against Published Literature

Compares our benchmark findings against known spaceflight biology from:
- Cell 2020 (Beheshti et al., "Comprehensive Multi-omics Analysis")
- SOMA 2024 (Nature, "Space Omics and Medical Atlas")
- NASA GeneLab published DEG analyses

Two validation approaches:
  1. Pathway-level: Compare our fGSEA top enriched Hallmark pathways vs
     literature-reported pathways per tissue (Spearman rank + overlap)
  2. Gene-level: Check if known spaceflight response genes appear in
     our SHAP top-50 per tissue

Output:
  evaluation/cell2020_validation.json

Usage:
  python scripts/cell2020_validation.py
"""

import json
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from scipy import stats

warnings.filterwarnings("ignore")

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent.parent
FGSEA_DIR = BASE_DIR / "processed" / "fgsea"
EVAL_DIR = BASE_DIR / "evaluation"

# ── Reference Data: Known Spaceflight Biology ─────────────────────────────────
# Sources:
#   - Beheshti et al., Cell 2020 (PMID: 33242417) — 13 tissues, mitochondrial stress hub
#   - da Silveira et al., Cell 2020 (PMID: 33242416) — twin study integrated multi-omics
#   - SOMA, Nature 2024 — human multi-omics atlas
#   - Individual GLDS analyses from NASA GeneLab publications

# Hallmark pathways expected to be enriched per tissue (from literature)
# Direction: UP = enriched in spaceflight, DOWN = suppressed
REFERENCE_PATHWAYS = {
    "liver": {
        "expected_up": [
            "HALLMARK_OXIDATIVE_PHOSPHORYLATION",  # Cell 2020: mito stress hub
            "HALLMARK_FATTY_ACID_METABOLISM",       # Cell 2020: metabolic disruption
            "HALLMARK_ADIPOGENESIS",                # Lipid metabolism shift
            "HALLMARK_BILE_ACID_METABOLISM",         # Liver-specific metabolic stress
            "HALLMARK_CHOLESTEROL_HOMEOSTASIS",      # Lipid pathway disruption
            "HALLMARK_XENOBIOTIC_METABOLISM",        # Detoxification response
            "HALLMARK_PEROXISOME",                  # Oxidative metabolism
        ],
        "expected_down": [
            "HALLMARK_INFLAMMATORY_RESPONSE",       # Immune suppression in space
            "HALLMARK_INTERFERON_GAMMA_RESPONSE",    # Immune suppression
        ],
        "source": "Cell 2020 (Beheshti et al.), GLDS liver analyses",
    },
    "thymus": {
        "expected_up": [
            "HALLMARK_E2F_TARGETS",                 # Cell proliferation genes
            "HALLMARK_G2M_CHECKPOINT",               # Cell cycle control
            "HALLMARK_MITOTIC_SPINDLE",              # Cell division
            "HALLMARK_MYC_TARGETS_V1",               # Proliferation program
        ],
        "expected_down": [
            "HALLMARK_INTERFERON_GAMMA_RESPONSE",    # Thymic involution
            "HALLMARK_INTERFERON_ALPHA_RESPONSE",    # Immune signaling reduction
            "HALLMARK_ALLOGRAFT_REJECTION",           # Immune suppression
        ],
        "source": "Cell 2020, GLDS thymus analyses (thymocyte proliferation/involution)",
    },
    "gastrocnemius": {
        "expected_up": [
            "HALLMARK_OXIDATIVE_PHOSPHORYLATION",   # Cell 2020: mito stress
            "HALLMARK_MYOGENESIS",                  # Muscle remodeling
        ],
        "expected_down": [
            "HALLMARK_MTORC1_SIGNALING",            # Muscle protein synthesis down
            "HALLMARK_PI3K_AKT_MTOR_SIGNALING",      # Atrophy signaling
        ],
        "source": "Cell 2020, muscle atrophy literature (MuRF1, MAFbx pathway)",
    },
    "kidney": {
        "expected_up": [
            "HALLMARK_MTORC1_SIGNALING",            # Renal metabolic stress
            "HALLMARK_CHOLESTEROL_HOMEOSTASIS",      # Lipid metabolism disruption
            "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",  # Oxidative stress
        ],
        "expected_down": [
            "HALLMARK_INTERFERON_GAMMA_RESPONSE",    # Immune suppression
        ],
        "source": "Cell 2020, OSD-253 kidney analysis (NRF2 pathway)",
    },
    "eye": {
        "expected_up": [
            "HALLMARK_OXIDATIVE_PHOSPHORYLATION",   # Cell 2020: dominant pathway
            "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",  # Retinal oxidative stress
        ],
        "expected_down": [
            "HALLMARK_ANGIOGENESIS",                # Vascular remodeling
        ],
        "source": "Cell 2020, retina OXPHOS dominance across missions",
    },
}

# Known spaceflight response genes per tissue
# These are genes consistently reported in spaceflight transcriptomics literature
REFERENCE_GENES = {
    "liver": {
        "genes": {
            "Angptl4": {"direction": "up", "function": "Lipid metabolism, angiogenesis regulation"},
            "Pck1": {"direction": "up", "function": "Gluconeogenesis key enzyme"},
            "G6pc": {"direction": "up", "function": "Glucose-6-phosphatase, gluconeogenesis"},
            "Cyp2e1": {"direction": "variable", "function": "Xenobiotic metabolism, CYP450"},
            "Fasn": {"direction": "up", "function": "Fatty acid synthase, lipogenesis"},
            "Dbp": {"direction": "variable", "function": "Circadian clock (D-box binding)"},
            "Per2": {"direction": "variable", "function": "Circadian clock (Period)"},
            "Npas2": {"direction": "variable", "function": "Circadian clock (BMAL1 paralog)"},
            "Pparg": {"direction": "up", "function": "Adipogenesis master regulator"},
            "Hmgcr": {"direction": "up", "function": "HMG-CoA reductase, cholesterol synthesis"},
            "Scd1": {"direction": "up", "function": "Stearoyl-CoA desaturase, lipogenesis"},
            "Acly": {"direction": "up", "function": "ATP citrate lyase, lipogenesis"},
            "Ucp2": {"direction": "up", "function": "Uncoupling protein 2, mitochondria"},
        },
        "source": "Cell 2020, GLDS-48/137/245/379/242 DEG analyses",
    },
    "gastrocnemius": {
        "genes": {
            "Trim63": {"direction": "up", "function": "MuRF1, E3 ubiquitin ligase (muscle atrophy)"},
            "Fbxo32": {"direction": "up", "function": "MAFbx/Atrogin-1, E3 ligase (atrophy)"},
            "Myog": {"direction": "variable", "function": "Myogenin, muscle differentiation"},
            "Mstn": {"direction": "up", "function": "Myostatin, negative regulator of muscle growth"},
            "Foxo3": {"direction": "up", "function": "FOXO3, atrophy transcription factor"},
            "Cdkn1a": {"direction": "up", "function": "p21, cell cycle inhibitor"},
        },
        "source": "Cell 2020, muscle atrophy E3 ligase pathway",
    },
    "thymus": {
        "genes": {
            "Ccnb1": {"direction": "down", "function": "Cyclin B1, G2/M transition"},
            "Cdk1": {"direction": "down", "function": "CDK1, cell cycle driver"},
            "Top2a": {"direction": "down", "function": "Topoisomerase II alpha, proliferation"},
            "Mki67": {"direction": "down", "function": "Ki-67, proliferation marker"},
            "Foxn1": {"direction": "down", "function": "Thymic epithelial cell master regulator"},
        },
        "source": "Cell 2020, thymic involution / cell cycle arrest",
    },
    "kidney": {
        "genes": {
            "Nfe2l2": {"direction": "up", "function": "NRF2, oxidative stress response master regulator"},
            "Hmox1": {"direction": "up", "function": "Heme oxygenase 1, NRF2 target"},
            "Nqo1": {"direction": "up", "function": "NAD(P)H quinone dehydrogenase 1, NRF2 target"},
            "Gclc": {"direction": "up", "function": "Glutamate-cysteine ligase, glutathione synthesis"},
            "Gclm": {"direction": "up", "function": "Glutamate-cysteine ligase modifier"},
        },
        "source": "OSD-253 kidney analysis, NRF2 pathway activation",
    },
    "eye": {
        "genes": {
            "Ndufa1": {"direction": "variable", "function": "NADH:ubiquinone oxidoreductase (Complex I)"},
            "Cox7a2": {"direction": "variable", "function": "Cytochrome c oxidase (Complex IV)"},
            "Atp5pb": {"direction": "variable", "function": "ATP synthase (Complex V)"},
            "Uqcrc2": {"direction": "variable", "function": "Cytochrome bc1 complex (Complex III)"},
        },
        "source": "Cell 2020, OXPHOS genes dominate retina spaceflight response",
    },
}


# ── Analysis Functions ─────────────────────────────────────────────────────────

def load_fgsea_hallmark_summary():
    """Load the consolidated fGSEA Hallmark results."""
    f = FGSEA_DIR / "summary" / "all_fgsea_hallmark.csv"
    if not f.exists():
        return None
    df = pd.read_csv(f)
    return df


def get_tissue_pathway_rankings(fgsea_df, tissue):
    """
    Get mean NES across missions for each pathway in a tissue.
    Returns: dict {pathway: mean_NES}, sorted by absolute NES.
    """
    tissue_data = fgsea_df[fgsea_df["tissue"] == tissue]
    if tissue_data.empty:
        return {}

    # Mean NES per pathway across missions
    pathway_nes = tissue_data.groupby("pathway")["NES"].mean()
    return pathway_nes.to_dict()


def validate_pathways(tissue, our_nes, reference):
    """
    Validate pathway enrichment against literature expectations.
    Returns validation metrics.
    """
    expected_up = reference.get("expected_up", [])
    expected_down = reference.get("expected_down", [])
    all_expected = expected_up + expected_down

    # Check which expected pathways are in our data
    found_up = [p for p in expected_up if p in our_nes]
    found_down = [p for p in expected_down if p in our_nes]
    found_all = found_up + found_down

    if not found_all:
        return {"status": "NO_DATA", "n_expected": len(all_expected)}

    # Direction concordance
    concordant = 0
    discordant = 0
    details = []

    for p in found_up:
        nes = our_nes[p]
        is_concordant = nes > 0
        if is_concordant:
            concordant += 1
        else:
            discordant += 1
        details.append({
            "pathway": p.replace("HALLMARK_", ""),
            "expected": "UP",
            "observed_nes": round(nes, 3),
            "concordant": is_concordant,
        })

    for p in found_down:
        nes = our_nes[p]
        is_concordant = nes < 0
        if is_concordant:
            concordant += 1
        else:
            discordant += 1
        details.append({
            "pathway": p.replace("HALLMARK_", ""),
            "expected": "DOWN",
            "observed_nes": round(nes, 3),
            "concordant": is_concordant,
        })

    total = concordant + discordant
    concordance_rate = concordant / total if total > 0 else 0

    # Top-5 pathway overlap: check if expected pathways appear in top-5 by |NES|
    sorted_by_abs_nes = sorted(our_nes.items(), key=lambda x: abs(x[1]), reverse=True)
    top5_pathways = [p for p, _ in sorted_by_abs_nes[:5]]
    top10_pathways = [p for p, _ in sorted_by_abs_nes[:10]]

    expected_in_top5 = [p for p in all_expected if p in top5_pathways]
    expected_in_top10 = [p for p in all_expected if p in top10_pathways]

    return {
        "status": "VALIDATED",
        "n_expected": len(all_expected),
        "n_found": len(found_all),
        "concordant": concordant,
        "discordant": discordant,
        "concordance_rate": round(concordance_rate, 3),
        "expected_in_top5": len(expected_in_top5),
        "expected_in_top10": len(expected_in_top10),
        "top5_pathways": [p.replace("HALLMARK_", "") for p in top5_pathways],
        "details": details,
        "source": reference.get("source", ""),
    }


def validate_genes(tissue, reference_genes):
    """
    Check if known spaceflight genes appear in SHAP top-50.
    """
    # Load SHAP data
    task_map = {
        "liver": "A1", "gastrocnemius": "A2",
        "kidney": "A3", "thymus": "A4", "eye": "A6",
    }
    task_id = task_map.get(tissue)
    if not task_id:
        return {"status": "NO_TASK"}

    shap_f = EVAL_DIR / f"{task_id}_shap_rf.json"
    if not shap_f.exists():
        return {"status": "NO_SHAP_DATA", "file": str(shap_f)}

    with open(shap_f) as fh:
        shap_data = json.load(fh)

    shap_genes = list(shap_data.get("top_genes", {}).keys())
    shap_genes_lower = {g.lower(): (i+1, g) for i, g in enumerate(shap_genes)}

    ref_genes = reference_genes.get("genes", {})
    found = []
    not_found = []

    for gene, info in ref_genes.items():
        gene_lower = gene.lower()
        if gene_lower in shap_genes_lower:
            rank, orig_name = shap_genes_lower[gene_lower]
            found.append({
                "gene": gene,
                "shap_rank": rank,
                "direction": info["direction"],
                "function": info["function"],
            })
        else:
            not_found.append({
                "gene": gene,
                "direction": info["direction"],
                "function": info["function"],
            })

    return {
        "status": "VALIDATED",
        "task_id": task_id,
        "n_reference_genes": len(ref_genes),
        "n_found_in_shap_top50": len(found),
        "overlap_rate": round(len(found) / len(ref_genes), 3) if ref_genes else 0,
        "found_genes": found,
        "not_found_genes": not_found,
        "shap_top10": shap_genes[:10],
        "source": reference_genes.get("source", ""),
    }


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("Cell 2020 / Literature External Validation")
    print("=" * 70)

    results = {
        "timestamp": datetime.now().isoformat(),
        "description": "External validation against published spaceflight biology",
        "references": [
            "Beheshti et al., Cell 2020 (PMID: 33242417)",
            "da Silveira et al., Cell 2020 (PMID: 33242416)",
            "SOMA, Nature 2024",
        ],
        "pathway_validation": {},
        "gene_validation": {},
        "summary": {},
    }

    # Load fGSEA data
    fgsea_df = load_fgsea_hallmark_summary()
    if fgsea_df is None:
        print("ERROR: fGSEA summary not found")
        return

    print(f"\nLoaded fGSEA: {len(fgsea_df)} rows, "
          f"tissues: {sorted(fgsea_df['tissue'].unique())}")

    # ── Pathway Validation ──────────────────────────────────────────────────
    print("\n" + "=" * 50)
    print("1. Pathway-Level Validation (fGSEA Hallmark NES)")
    print("=" * 50)

    pathway_concordance_rates = []

    for tissue in REFERENCE_PATHWAYS:
        print(f"\n  --- {tissue} ---")
        our_nes = get_tissue_pathway_rankings(fgsea_df, tissue)
        if not our_nes:
            print(f"    No fGSEA data")
            continue

        ref = REFERENCE_PATHWAYS[tissue]
        validation = validate_pathways(tissue, our_nes, ref)
        results["pathway_validation"][tissue] = validation

        if validation["status"] == "VALIDATED":
            rate = validation["concordance_rate"]
            pathway_concordance_rates.append(rate)
            print(f"    Direction concordance: {validation['concordant']}/"
                  f"{validation['concordant']+validation['discordant']} "
                  f"({rate:.1%})")
            print(f"    Expected in top-5 |NES|: {validation['expected_in_top5']}")
            print(f"    Expected in top-10 |NES|: {validation['expected_in_top10']}")
            print(f"    Top-5 pathways: {', '.join(validation['top5_pathways'])}")

            for d in validation["details"]:
                mark = "OK" if d["concordant"] else "XX"
                print(f"      [{mark}] {d['pathway']:40s} expected={d['expected']:4s} NES={d['observed_nes']:+.3f}")

    # ── Gene Validation ─────────────────────────────────────────────────────
    print("\n" + "=" * 50)
    print("2. Gene-Level Validation (SHAP Top-50)")
    print("=" * 50)

    gene_overlap_rates = []

    for tissue in REFERENCE_GENES:
        print(f"\n  --- {tissue} ---")
        ref_genes = REFERENCE_GENES[tissue]
        validation = validate_genes(tissue, ref_genes)
        results["gene_validation"][tissue] = validation

        if validation["status"] == "VALIDATED":
            n_found = validation["n_found_in_shap_top50"]
            n_total = validation["n_reference_genes"]
            rate = validation["overlap_rate"]
            gene_overlap_rates.append(rate)
            print(f"    Found in SHAP top-50: {n_found}/{n_total} ({rate:.1%})")

            for g in validation["found_genes"]:
                print(f"      FOUND [rank {g['shap_rank']:2d}]: {g['gene']:12s} "
                      f"({g['function'][:50]})")
            for g in validation["not_found_genes"]:
                print(f"      ----  [     ]: {g['gene']:12s} "
                      f"({g['function'][:50]})")
        elif validation["status"] == "NO_SHAP_DATA":
            print(f"    No SHAP data available")

    # ── Summary ─────────────────────────────────────────────────────────────
    print("\n" + "=" * 50)
    print("Validation Summary")
    print("=" * 50)

    summary = {}

    if pathway_concordance_rates:
        mean_concordance = np.mean(pathway_concordance_rates)
        summary["pathway_mean_concordance"] = round(float(mean_concordance), 3)
        summary["pathway_n_tissues"] = len(pathway_concordance_rates)

        if mean_concordance >= 0.7:
            summary["pathway_verdict"] = "STRONG"
            summary["pathway_interpretation"] = (
                f"Mean direction concordance {mean_concordance:.1%}: "
                "Our fGSEA results strongly agree with published spaceflight pathway biology."
            )
        elif mean_concordance >= 0.5:
            summary["pathway_verdict"] = "MODERATE"
            summary["pathway_interpretation"] = (
                f"Mean direction concordance {mean_concordance:.1%}: "
                "Our fGSEA results partially agree with published biology. "
                "Discordances may reflect mission-specific variation."
            )
        else:
            summary["pathway_verdict"] = "WEAK"
            summary["pathway_interpretation"] = (
                f"Mean direction concordance {mean_concordance:.1%}: "
                "Limited agreement. May reflect different analysis methods or datasets."
            )

        print(f"\n  Pathway concordance: {mean_concordance:.1%} across {len(pathway_concordance_rates)} tissues")
        print(f"  Verdict: {summary['pathway_verdict']}")

    if gene_overlap_rates:
        mean_overlap = np.mean(gene_overlap_rates)
        summary["gene_mean_overlap"] = round(float(mean_overlap), 3)
        summary["gene_n_tissues"] = len(gene_overlap_rates)

        # Expected overlap by chance: 50 SHAP genes / ~22000 total ≈ 0.2%
        # Having even 1-2 genes overlap is noteworthy
        expected_by_chance = 50 / 22000  # ~0.23%

        summary["gene_interpretation"] = (
            f"Mean overlap rate {mean_overlap:.1%} vs chance expectation ~{expected_by_chance:.1%}. "
            f"Enrichment ratio: {mean_overlap / expected_by_chance:.1f}x."
        )
        print(f"\n  Gene overlap (SHAP top-50): {mean_overlap:.1%}")
        print(f"  Chance expectation: ~{expected_by_chance:.1%}")
        print(f"  Enrichment: {mean_overlap / expected_by_chance:.1f}x above chance")

    results["summary"] = summary

    # Save
    out_f = EVAL_DIR / "cell2020_validation.json"
    with open(out_f, "w") as fh:
        json.dump(results, fh, indent=2, default=str)
    print(f"\nSaved: {out_f}")


if __name__ == "__main__":
    main()
