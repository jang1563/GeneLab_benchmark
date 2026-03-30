#!/usr/bin/env python3
"""
Phase B: Pathway-Level Conservation — Mouse ssGSEA NES vs Human cfRNA NES

Correlate mouse tissue pathway NES with human cfRNA pathway NES across
50 Hallmark pathways. Test which mouse tissue best predicts human response.
"""

import sys
import os
import numpy as np
import pandas as pd
from scipy import stats
from datetime import datetime

sys.path.insert(0, os.path.dirname(__file__))
from v6_utils import (
    load_cfrna_de, load_pathway_scores_hallmark,
    ALL_TISSUES, V6_EVAL, GENELAB_BASE, PROCESSED, save_json
)


def _download_gene_sets():
    """Download Hallmark gene sets using gseapy's library download (cached)."""
    import gseapy as gp
    # Try multiple approaches to get gene sets
    try:
        # Method 1: gseapy get_library (downloads from Enrichr and caches)
        lib = gp.get_library("MSigDB_Hallmark_2020", organism="Human")
        if lib:
            return lib, "enrichr_download"
    except Exception as e:
        print(f"  get_library failed: {e}")

    try:
        # Method 2: gseapy get_library_name to verify, then manual download
        import json
        import urllib.request
        url = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020"
        response = urllib.request.urlopen(url, timeout=30)
        text = response.read().decode()
        gene_sets = {}
        for line in text.strip().split("\n"):
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                name = parts[0]
                genes = [g for g in parts[2:] if g.strip()]
                gene_sets[name] = genes
        if gene_sets:
            return gene_sets, "manual_download"
    except Exception as e:
        print(f"  Manual download failed: {e}")

    return None, None


def compute_human_pathway_nes(cfrna_de):
    """Compute human cfRNA pathway NES using gene-level pre_vs_flight statistics.

    Uses gseapy prerank with locally downloaded gene sets.
    Falls back to simple mean-of-ranks enrichment score.
    """
    import gseapy as gp

    # Prepare gene ranking: gene → pre_vs_flight_diff
    ranking = cfrna_de["edge_pre_vs_flight_diff"].dropna()
    ranking = ranking.sort_values(ascending=False)

    # Download gene sets first (cached by gseapy)
    print("  Downloading Hallmark gene sets...")
    gene_sets, dl_method = _download_gene_sets()

    if gene_sets:
        print(f"  Got {len(gene_sets)} gene sets via {dl_method}")
        try:
            res = gp.prerank(
                rnk=ranking,
                gene_sets=gene_sets,  # Pass dict directly
                min_size=10,
                max_size=500,
                permutation_num=1000,
                seed=42,
                no_plot=True,
                verbose=False,
            )

            res_df = res.res2d
            if "Term" in res_df.columns:
                term_col = "Term"
            elif "Name" in res_df.columns:
                term_col = "Name"
            else:
                term_col = res_df.columns[0]

            nes_dict = dict(zip(res_df[term_col], res_df["NES"].astype(float)))
            return nes_dict, "gseapy_prerank"
        except Exception as e:
            print(f"  gseapy prerank failed: {e}")
            print("  Using fallback: z-score enrichment")
            return _fallback_pathway_nes(ranking, gene_sets), "fallback_zscore"
    else:
        print("  Could not download gene sets")
        return {}, "failed"


def _fallback_pathway_nes(ranking, gene_sets):
    """Fallback: compute z-score enrichment from ranking data."""
    nes_dict = {}
    for name, genes in gene_sets.items():
        if isinstance(genes, list):
            overlap = [g for g in genes if g in ranking.index]
        else:
            overlap = [g for g in genes.split(",") if g.strip() in ranking.index]
        if len(overlap) >= 5:
            values = ranking.loc[overlap].values
            if np.std(values) > 0:
                nes = np.mean(values) / (np.std(values) / np.sqrt(len(overlap)))
            else:
                nes = 0
            nes_dict[name] = float(nes)
    return nes_dict


def compute_mouse_pathway_nes(tissue):
    """Compute mouse pathway NES as mean(FLT) - mean(GC) per Hallmark pathway.

    Returns dict {pathway_name: NES_diff} or None if data unavailable.
    """
    pw_df = load_pathway_scores_hallmark(tissue)
    if pw_df is None or pw_df.empty:
        return None

    # Load metadata to get FLT/GC labels
    sys.path.insert(0, str(GENELAB_BASE / "v4" / "scripts"))
    from v4_utils import load_metadata, encode_labels

    meta = load_metadata(tissue)
    y, valid_mask = encode_labels(meta)
    meta = meta[valid_mask]
    y = y[valid_mask]

    # Align pathway scores with metadata
    common_idx = pw_df.index.intersection(meta.index)
    if len(common_idx) < 5:
        return None

    pw_df = pw_df.loc[common_idx]
    y_aligned = y.loc[common_idx]

    # FLT mean - GC mean per pathway
    flt_mask = y_aligned == 1
    gc_mask = y_aligned == 0

    if flt_mask.sum() < 2 or gc_mask.sum() < 2:
        return None

    nes_dict = {}
    for col in pw_df.columns:
        flt_mean = pw_df.loc[flt_mask, col].mean()
        gc_mean = pw_df.loc[gc_mask, col].mean()
        pooled_std = pw_df[col].std()
        # Effect size (Cohen's d approximation)
        if pooled_std > 0:
            nes_dict[col] = float((flt_mean - gc_mean) / pooled_std)
        else:
            nes_dict[col] = 0.0

    return nes_dict


def normalize_pathway_names(name):
    """Normalize pathway names for matching across databases.
    Handles: HALLMARK_OXIDATIVE_PHOSPHORYLATION vs Oxidative Phosphorylation
    """
    name = name.upper().replace("HALLMARK_", "").replace(" ", "_")
    return name


def match_pathways(mouse_dict, human_dict):
    """Match pathways between mouse ssGSEA and human prerank by normalized name.
    Returns list of (pathway, mouse_nes, human_nes) tuples.
    """
    # Normalize both sets
    mouse_norm = {normalize_pathway_names(k): (k, v) for k, v in mouse_dict.items()}
    human_norm = {normalize_pathway_names(k): (k, v) for k, v in human_dict.items()}

    matched = []
    for norm_name in sorted(set(mouse_norm.keys()) & set(human_norm.keys())):
        mk, mv = mouse_norm[norm_name]
        hk, hv = human_norm[norm_name]
        matched.append((norm_name, mv, hv, mk, hk))

    return matched


def main():
    print("=" * 60)
    print("Phase B: Pathway-Level Conservation")
    print("=" * 60)

    # Step 1: Compute human cfRNA pathway NES
    print("\n1. Computing human cfRNA pathway NES...")
    cfrna_de = load_cfrna_de()
    human_nes, human_method = compute_human_pathway_nes(cfrna_de)
    print(f"  Method: {human_method}")
    print(f"  Human pathways: {len(human_nes)}")

    if not human_nes:
        print("  ERROR: No human pathway NES computed. Exiting.")
        return

    # Step 2: Compute mouse pathway NES per tissue
    print("\n2. Computing mouse pathway NES per tissue...")
    mouse_nes_all = {}
    for tissue in ALL_TISSUES:
        nes = compute_mouse_pathway_nes(tissue)
        if nes:
            mouse_nes_all[tissue] = nes
            print(f"  {tissue}: {len(nes)} pathways")
        else:
            print(f"  {tissue}: SKIPPED (no pathway data)")

    # Step 3: Correlate mouse vs human per tissue
    print("\n3. Correlating mouse vs human pathway NES...")
    correlation_results = {}

    for tissue, mouse_nes in mouse_nes_all.items():
        matched = match_pathways(mouse_nes, human_nes)
        n_matched = len(matched)

        if n_matched < 5:
            print(f"  {tissue}: only {n_matched} matched pathways, SKIPPING")
            continue

        mouse_vals = np.array([m[1] for m in matched])
        human_vals = np.array([m[2] for m in matched])

        # Z-score both for fair comparison
        mouse_z = (mouse_vals - np.mean(mouse_vals)) / np.std(mouse_vals) \
            if np.std(mouse_vals) > 0 else mouse_vals
        human_z = (human_vals - np.mean(human_vals)) / np.std(human_vals) \
            if np.std(human_vals) > 0 else human_vals

        # Spearman correlation
        rho, p_val = stats.spearmanr(mouse_z, human_z)

        # Bootstrap CI (n=2000)
        rng = np.random.default_rng(42)
        boot_rhos = []
        for _ in range(2000):
            idx = rng.integers(0, n_matched, size=n_matched)
            r, _ = stats.spearmanr(mouse_z[idx], human_z[idx])
            if not np.isnan(r):
                boot_rhos.append(r)
        ci_lower = np.percentile(boot_rhos, 2.5) if boot_rhos else np.nan
        ci_upper = np.percentile(boot_rhos, 97.5) if boot_rhos else np.nan

        # Top concordant and discordant pathways
        diffs = [(m[0], m[1], m[2], m[1] * m[2]) for m in matched]
        concordant = sorted([d for d in diffs if d[3] > 0],
                          key=lambda x: abs(x[3]), reverse=True)
        discordant = sorted([d for d in diffs if d[3] < 0],
                          key=lambda x: abs(x[3]), reverse=True)

        correlation_results[tissue] = {
            "n_matched_pathways": n_matched,
            "spearman_rho": round(float(rho), 4),
            "p_value": float(p_val),
            "ci_lower": round(float(ci_lower), 4),
            "ci_upper": round(float(ci_upper), 4),
            "significant": p_val < 0.05,
            "top_concordant": [
                {"pathway": d[0], "mouse_nes": round(d[1], 3),
                 "human_nes": round(d[2], 3)}
                for d in concordant[:5]
            ],
            "top_discordant": [
                {"pathway": d[0], "mouse_nes": round(d[1], 3),
                 "human_nes": round(d[2], 3)}
                for d in discordant[:5]
            ],
        }

        sig = "*" if p_val < 0.05 else ""
        print(f"  {tissue:>15s}: r={rho:+.3f} [{ci_lower:+.3f}, {ci_upper:+.3f}] "
              f"p={p_val:.4f}{sig} ({n_matched} pathways)")

    # Mean correlation across tissues
    if correlation_results:
        rhos = [v["spearman_rho"] for v in correlation_results.values()]
        mean_rho = np.mean(rhos)
        print(f"\n  Mean Spearman rho: {mean_rho:+.3f}")

    # Assemble output
    output = {
        "human_method": human_method,
        "n_human_pathways": len(human_nes),
        "human_nes": {k: round(v, 4) for k, v in human_nes.items()},
        "n_tissues_analyzed": len(correlation_results),
        "per_tissue_correlations": correlation_results,
        "mean_rho": round(float(mean_rho), 4) if correlation_results else None,
        "timestamp": datetime.now().isoformat(),
    }

    save_json(output, V6_EVAL / "V6_B_pathway_conservation.json")
    print("\nDone!")


if __name__ == "__main__":
    main()
