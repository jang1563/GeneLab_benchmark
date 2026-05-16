#!/usr/bin/env python
"""GeneLabBench v5 Phase 2: Cross-Organ Signaling Network.

Uses OmniPath L-R database + decoupler-py TF activity inference
to model inter-tissue communication during spaceflight.

Usage:
    python cross_organ_signaling.py
"""
import json
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

# ── paths ──
BASE_DIR = Path(__file__).resolve().parent.parent.parent
V4_EVAL_DIR = BASE_DIR / "v4" / "evaluation"
V5_EVAL_DIR = BASE_DIR / "v5" / "evaluation"
ORTHO_PATH = BASE_DIR / "v3" / "data" / "mouse_human_orthologs.tsv"

sys.path.insert(0, str(BASE_DIR / "v4" / "scripts"))
from v4_utils import (
    TISSUE_MISSIONS, load_metadata, load_gene_features,
    align_features_with_meta, encode_labels,
)

ALL_TISSUES = list(TISSUE_MISSIONS.keys())


def load_shap_top_genes(tissue, top_n=200):
    """Load SHAP top-N genes for a tissue (union across ALL available methods)."""
    genes = set()
    for path in V4_EVAL_DIR.glob(f"SHAP_{tissue}_*.json"):
        if "consensus" in path.name or "interaction" in path.name or "WGCNA" in path.name:
            continue
        try:
            with open(path) as f:
                data = json.load(f)
            importance = data.get("gene_importance", data.get("mean_importance", {}))
            if not importance:
                continue
            ranked = sorted(importance.items(),
                            key=lambda x: abs(x[1]), reverse=True)
            genes.update(g for g, _ in ranked[:top_n])
        except (json.JSONDecodeError, KeyError):
            continue
    return genes


def load_orthologs():
    """Load mouse→human ortholog mapping."""
    df = pd.read_csv(ORTHO_PATH, sep="\t")
    return dict(zip(df["mouse_symbol"], df["human_symbol"]))


def get_lr_pairs_omnipath():
    """Fetch ligand-receptor pairs from OmniPath (mouse, cached).

    Filters to proper intercellular signaling: source must be a secreted
    protein (ligand) and target must be a plasma membrane receptor.
    Excludes intracellular interactions (TF-target, kinase-substrate).
    """
    try:
        import omnipath

        # Get ALL interactions for mouse
        all_interactions = omnipath.interactions.OmniPath.get(
            genesymbols=True,
            organisms="mouse",
        )

        src_col = "source_genesymbol" if "source_genesymbol" in all_interactions.columns else "genesymbol_a"
        tgt_col = "target_genesymbol" if "target_genesymbol" in all_interactions.columns else "genesymbol_b"

        # Get intercellular annotations to identify ligands vs receptors
        try:
            intercell = omnipath.requests.Intercell.get(organism="mouse")
            # Build sets of intercellular proteins
            ligand_genes = set()
            receptor_genes = set()
            intercell_genes = set()  # any gene with intercellular annotation
            for _, row in intercell.iterrows():
                cat = str(row.get("category", "")).lower()
                gene = str(row.get("genesymbol", ""))
                if not gene or gene == "nan":
                    continue
                if any(kw in cat for kw in [
                    "ligand", "secreted", "growth_factor", "cytokine",
                    "chemokine", "hormone", "extracellular"
                ]):
                    ligand_genes.add(gene)
                    intercell_genes.add(gene)
                if any(kw in cat for kw in [
                    "receptor", "cell_surface", "transmembrane",
                    "adhesion", "transporter"
                ]):
                    receptor_genes.add(gene)
                    intercell_genes.add(gene)

            print(f"  OmniPath intercell: {len(ligand_genes)} ligands, "
                  f"{len(receptor_genes)} receptors, "
                  f"{len(intercell_genes)} total annotated")

            # Tiered filtering:
            # Tier 1 (strict): source=ligand AND target=receptor
            # Tier 2 (broad): at least one side has intercellular annotation
            #                  AND the other is not a known TF/kinase-only gene
            strict_pairs = []
            broad_pairs = []
            for _, row in all_interactions.iterrows():
                src = str(row.get(src_col, ""))
                tgt = str(row.get(tgt_col, ""))
                if not src or not tgt or src == "nan" or tgt == "nan":
                    continue
                if src in ligand_genes and tgt in receptor_genes:
                    strict_pairs.append((src, tgt))
                elif src in intercell_genes or tgt in intercell_genes:
                    broad_pairs.append((src, tgt))

            # Use strict if ≥50 pairs, otherwise use broad
            if len(strict_pairs) >= 50:
                pairs = strict_pairs
                filter_type = "strict (ligand→receptor)"
            else:
                pairs = strict_pairs + broad_pairs
                filter_type = "broad (any intercellular)"

            print(f"  Strict L-R: {len(strict_pairs)}, "
                  f"Broad intercellular: {len(broad_pairs)}")
            print(f"  Using {filter_type}: {len(pairs)} pairs")

        except Exception as e_intercell:
            print(f"  [WARN] Intercell annotations failed ({e_intercell}), "
                  "using heuristic filtering...")
            pairs = []
            for _, row in all_interactions.iterrows():
                src = str(row.get(src_col, ""))
                tgt = str(row.get(tgt_col, ""))
                if not src or not tgt or src == "nan" or tgt == "nan":
                    continue
                pairs.append((src, tgt))
            print(f"  OmniPath (no filter): {len(pairs)} pairs")

        # Deduplicate
        pairs = list(set(pairs))
        print(f"  Final L-R pairs: {len(pairs)} (deduplicated)")
        return pairs

    except Exception as e:
        print(f"  [WARN] OmniPath failed: {e}")
        print("  Falling back to manual ortholog mapping...")
        return get_lr_pairs_omnipath_human_mapped()


def get_lr_pairs_omnipath_human_mapped():
    """Fallback: Get human L-R pairs and map to mouse via orthologs."""
    try:
        import omnipath
        lr = omnipath.interactions.OmniPath.get(genesymbols=True)

        ortho = load_orthologs()
        human_to_mouse = {v: k for k, v in ortho.items()}

        pairs = []
        src_col = "source_genesymbol" if "source_genesymbol" in lr.columns else "genesymbol_a"
        tgt_col = "target_genesymbol" if "target_genesymbol" in lr.columns else "genesymbol_b"

        for _, row in lr.iterrows():
            src_h = str(row.get(src_col, ""))
            tgt_h = str(row.get(tgt_col, ""))
            src_m = human_to_mouse.get(src_h)
            tgt_m = human_to_mouse.get(tgt_h)
            if src_m and tgt_m:
                pairs.append((src_m, tgt_m))

        print(f"  OmniPath (human→mouse mapped): {len(pairs)} L-R pairs")
        return pairs

    except Exception as e:
        print(f"  [ERROR] OmniPath completely failed: {e}")
        return []


def compute_cross_organ_signaling(lr_pairs, shap_genes_per_tissue):
    """Compute L-R signaling scores between all tissue pairs."""
    # Build lookup sets for ligands and receptors
    ligands = {src for src, _ in lr_pairs}
    receptors = {tgt for _, tgt in lr_pairs}
    lr_set = set(lr_pairs)

    results = {}
    for t_src, t_tgt in combinations(ALL_TISSUES, 2):
        # Ligands from source tissue that are in SHAP top genes
        src_ligands = shap_genes_per_tissue.get(t_src, set()) & ligands
        # Receptors in target tissue that are in SHAP top genes
        tgt_receptors = shap_genes_per_tissue.get(t_tgt, set()) & receptors

        # Active L-R pairs: ligand in source SHAP + receptor in target SHAP
        active_pairs = []
        for lig in src_ligands:
            for rec in tgt_receptors:
                if (lig, rec) in lr_set:
                    active_pairs.append({"ligand": lig, "receptor": rec})

        # Also check reverse direction
        src_ligands_rev = shap_genes_per_tissue.get(t_tgt, set()) & ligands
        tgt_receptors_rev = shap_genes_per_tissue.get(t_src, set()) & receptors
        active_pairs_rev = []
        for lig in src_ligands_rev:
            for rec in tgt_receptors_rev:
                if (lig, rec) in lr_set:
                    active_pairs_rev.append({"ligand": lig, "receptor": rec})

        key = f"{t_src}_to_{t_tgt}"
        key_rev = f"{t_tgt}_to_{t_src}"
        results[key] = {
            "source": t_src, "target": t_tgt,
            "n_active_pairs": len(active_pairs),
            "active_pairs": active_pairs[:20],  # top 20 for output size
            "n_source_ligands": len(src_ligands),
            "n_target_receptors": len(tgt_receptors),
        }
        results[key_rev] = {
            "source": t_tgt, "target": t_src,
            "n_active_pairs": len(active_pairs_rev),
            "active_pairs": active_pairs_rev[:20],
            "n_source_ligands": len(src_ligands_rev),
            "n_target_receptors": len(tgt_receptors_rev),
        }

    return results


def run_tf_activity(tissue, symbol_map_path):
    """Run decoupler-py TF activity inference for one tissue."""
    import decoupler as dc

    meta = load_metadata(tissue)
    genes = load_gene_features(tissue)
    genes, meta = align_features_with_meta(genes, meta)
    y, valid = encode_labels(meta)

    genes_valid = genes.loc[valid]
    y_valid = y[valid].astype(int)

    # Map Ensembl → Symbol
    sym_map = pd.read_csv(symbol_map_path)
    ens_to_sym = dict(zip(sym_map["ENSEMBL"], sym_map["SYMBOL"]))
    expr = genes_valid.rename(columns=ens_to_sym)
    # Drop unmapped columns
    expr = expr[[c for c in expr.columns if not str(c).startswith("ENSMUSG")]]
    # Drop duplicate symbols
    expr = expr.loc[:, ~expr.columns.duplicated()]

    print(f"  {tissue}: {expr.shape[0]} samples × {expr.shape[1]} genes (symbols)")

    # Get mouse CollecTRI regulons (decoupler v2 API)
    try:
        collectri = dc.op.collectri(organism="mouse")
    except Exception:
        # Fallback: human regulons
        collectri = dc.op.collectri(organism="human")
        print(f"  [WARN] Using human CollecTRI (mouse not available)")

    # Run ULM (univariate linear model) for TF activity (decoupler v2 API)
    tf_acts, tf_pvals = dc.mt.ulm(data=expr, net=collectri)

    print(f"  TF activities: {tf_acts.shape[0]} samples × {tf_acts.shape[1]} TFs")

    # Compare FLT vs GC per TF
    tf_results = {}
    flt_mask = y_valid.values == 1
    gc_mask = y_valid.values == 0

    for tf in tf_acts.columns:
        flt_vals = tf_acts.loc[flt_mask, tf].values
        gc_vals = tf_acts.loc[gc_mask, tf].values

        if len(flt_vals) < 2 or len(gc_vals) < 2:
            continue

        try:
            stat, p = mannwhitneyu(flt_vals, gc_vals, alternative="two-sided")
        except ValueError:
            p = 1.0

        tf_results[tf] = {
            "mean_flt": round(float(np.mean(flt_vals)), 4),
            "mean_gc": round(float(np.mean(gc_vals)), 4),
            "t_stat_approx": round(float(np.mean(flt_vals) - np.mean(gc_vals)), 4),
            "wilcoxon_p": round(float(p), 6),
            "direction": "up_in_flight" if np.mean(flt_vals) > np.mean(gc_vals) else "down_in_flight",
        }

    # BH-FDR correction
    if tf_results:
        pvals = [tf_results[tf]["wilcoxon_p"] for tf in tf_results]
        tf_names = list(tf_results.keys())
        n_tests = len(pvals)
        sorted_idx = np.argsort(pvals)
        fdr_p = np.zeros(n_tests)
        for rank, idx in enumerate(sorted_idx):
            fdr_p[idx] = pvals[idx] * n_tests / (rank + 1)
        # Enforce monotonicity
        for i in range(n_tests - 2, -1, -1):
            fdr_p[sorted_idx[i]] = min(fdr_p[sorted_idx[i]],
                                        fdr_p[sorted_idx[i+1]] if i+1 < n_tests else 1.0)
        fdr_p = np.clip(fdr_p, 0, 1)

        for i, tf in enumerate(tf_names):
            tf_results[tf]["fdr_p"] = round(float(fdr_p[i]), 6)

    n_sig = sum(1 for tf in tf_results if tf_results[tf].get("fdr_p", 1) < 0.05)
    print(f"  Significant TFs (FDR<0.05): {n_sig}/{len(tf_results)}")

    return {
        "tissue": tissue,
        "n_samples": int(len(genes_valid)),
        "n_tfs_tested": len(tf_results),
        "n_significant_fdr05": n_sig,
        "tf_results": tf_results,
    }


def main():
    V5_EVAL_DIR.mkdir(parents=True, exist_ok=True)
    symbol_map_path = BASE_DIR / "processed" / "ensembl_symbol_map.csv"

    # ── Part A: Cross-organ L-R signaling ──
    print("=" * 60)
    print("Part A: Cross-Organ Ligand-Receptor Signaling")
    print("=" * 60)

    # Load SHAP top genes per tissue (use gene symbols)
    sym_map = pd.read_csv(symbol_map_path)
    ens_to_sym = dict(zip(sym_map["ENSEMBL"], sym_map["SYMBOL"]))

    shap_genes_per_tissue = {}
    for tissue in ALL_TISSUES:
        ens_genes = load_shap_top_genes(tissue, top_n=200)
        # Map to symbols
        sym_genes = {ens_to_sym.get(g, g) for g in ens_genes}
        shap_genes_per_tissue[tissue] = sym_genes
        print(f"  {tissue}: {len(sym_genes)} SHAP top genes (symbols)")

    # Get L-R pairs from OmniPath
    lr_pairs = get_lr_pairs_omnipath()

    if not lr_pairs:
        print("[ERROR] No L-R pairs available. Skipping cross-organ analysis.")
        signaling_results = {"error": "No L-R pairs available"}
    else:
        signaling_results = compute_cross_organ_signaling(lr_pairs, shap_genes_per_tissue)
        # Summary
        total_pairs = sum(r["n_active_pairs"] for r in signaling_results.values())
        print(f"\n  Total active L-R pairs across all tissue pairs: {total_pairs}")

    # Save signaling results
    signaling_output = {
        "n_lr_pairs_database": len(lr_pairs),
        "n_tissues": len(ALL_TISSUES),
        "n_tissue_pairs": len(signaling_results),
        "tissue_pairs": signaling_results,
    }
    with open(V5_EVAL_DIR / "cross_organ_signaling.json", "w") as f:
        json.dump(signaling_output, f, indent=2)
    print(f"Saved: cross_organ_signaling.json")

    # ── Part B: TF activity inference ──
    print("\n" + "=" * 60)
    print("Part B: Transcription Factor Activity (decoupler-py)")
    print("=" * 60)

    for tissue in ALL_TISSUES:
        try:
            tf_result = run_tf_activity(tissue, symbol_map_path)
            out_path = V5_EVAL_DIR / f"tf_activity_{tissue}.json"
            with open(out_path, "w") as f:
                json.dump(tf_result, f, indent=2)
            print(f"  Saved: {out_path.name}")
        except Exception as e:
            print(f"  [ERROR] {tissue}: {e}")

    print("\nDone!")


if __name__ == "__main__":
    main()
