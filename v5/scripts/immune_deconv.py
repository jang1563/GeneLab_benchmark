#!/usr/bin/env python
"""GeneLabBench v5 Phase 1: Immune Deconvolution via mMCP-counter.

Estimates immune/stromal cell type proportions from bulk RNA-seq,
then tests FLT vs GC differences per tissue.

Usage:
    python immune_deconv.py --tissue liver
    python immune_deconv.py  # all 8 tissues
"""
import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

# ── paths ──
BASE_DIR = Path(__file__).resolve().parent.parent.parent
PROCESSED_DIR = BASE_DIR / "processed" / "A_detection"
SYMBOL_MAP_PATH = BASE_DIR / "processed" / "ensembl_symbol_map.csv"
V5_EVAL_DIR = BASE_DIR / "v5" / "evaluation"
R_SCRIPT = Path(__file__).resolve().parent / "run_mmcp.R"

sys.path.insert(0, str(BASE_DIR / "v4" / "scripts"))
from v4_utils import (
    TISSUE_MISSIONS, LABEL_MAP, EXCLUDE_LABELS,
    load_metadata, load_gene_features, align_features_with_meta, encode_labels,
)

ALL_TISSUES = list(TISSUE_MISSIONS.keys())


def load_symbol_map():
    """Load Ensembl → Symbol mapping."""
    df = pd.read_csv(SYMBOL_MAP_PATH)
    return dict(zip(df["ENSEMBL"], df["SYMBOL"]))


def cliffs_delta(x, y):
    """Cliff's delta effect size (non-parametric)."""
    n_x, n_y = len(x), len(y)
    if n_x == 0 or n_y == 0:
        return 0.0
    more = sum(1 for xi in x for yi in y if xi > yi)
    less = sum(1 for xi in x for yi in y if xi < yi)
    return (more - less) / (n_x * n_y)


def run_deconv_for_tissue(tissue, symbol_map):
    """Run mMCP-counter deconvolution for one tissue."""
    print(f"\n{'='*60}")
    print(f"Tissue: {tissue}")
    print(f"{'='*60}")

    # Load expression + metadata
    meta = load_metadata(tissue)
    genes = load_gene_features(tissue)
    genes, meta = align_features_with_meta(genes, meta)

    # Encode labels → filter to FLT/GC only
    y, valid = encode_labels(meta)
    genes_valid = genes.loc[valid]
    meta_valid = meta.loc[valid]
    y_valid = y[valid].astype(int)

    print(f"  Samples: {len(meta_valid)} (FLT={sum(y_valid==1)}, GC={sum(y_valid==0)})")
    print(f"  Genes: {genes_valid.shape[1]}")

    # Map Ensembl → Symbol
    gene_cols = genes_valid.columns.tolist()
    mapped_symbols = []
    ensembl_to_sym = {}
    for ens_id in gene_cols:
        sym = symbol_map.get(ens_id, None)
        if sym and sym not in ensembl_to_sym.values():
            ensembl_to_sym[ens_id] = sym
            mapped_symbols.append(sym)

    # Build symbol-indexed expression (genes × samples for R)
    expr_sym = genes_valid[[e for e in ensembl_to_sym.keys()]].copy()
    expr_sym.columns = [ensembl_to_sym[e] for e in expr_sym.columns]
    expr_t = expr_sym.T  # genes × samples
    print(f"  Mapped to symbols: {len(mapped_symbols)} genes")

    # Write temp CSV, run R script, read result
    with tempfile.TemporaryDirectory() as tmpdir:
        input_csv = Path(tmpdir) / "expr.csv"
        output_csv = Path(tmpdir) / "deconv.csv"

        expr_t.to_csv(input_csv)
        print(f"  Running mMCP-counter via R...")

        result = subprocess.run(
            ["Rscript", str(R_SCRIPT), str(input_csv), str(output_csv)],
            capture_output=True, text=True, timeout=600,
        )
        if result.returncode != 0:
            print(f"  [ERROR] R script failed:\n{result.stderr[:500]}")
            return None

        print(f"  R stdout: {result.stdout.strip()[-200:]}")

        if not output_csv.exists():
            print(f"  [ERROR] Output CSV not created")
            return None

        deconv = pd.read_csv(output_csv, index_col=0)

    # deconv: cell_types × samples → transpose to samples × cell_types
    if deconv.shape[0] < deconv.shape[1]:
        deconv = deconv.T
    print(f"  Deconv matrix: {deconv.shape[0]} samples × {deconv.shape[1]} cell types")

    # Align deconv samples with metadata
    common_samples = deconv.index.intersection(meta_valid.index)
    if len(common_samples) < 5:
        # Try matching without mission prefix
        deconv_stripped = {s: s for s in deconv.index}
        meta_stripped = {s: s for s in meta_valid.index}
        # Fall back to positional alignment if sample names match count
        if deconv.shape[0] == len(meta_valid):
            print(f"  Using positional alignment ({deconv.shape[0]} samples)")
            deconv.index = meta_valid.index
            common_samples = meta_valid.index
        else:
            print(f"  [WARN] Only {len(common_samples)} common samples, skipping")
            return None

    deconv_aligned = deconv.loc[common_samples]
    y_aligned = y_valid.loc[common_samples]
    n_flt = int(sum(y_aligned == 1))
    n_gc = int(sum(y_aligned == 0))

    # Statistical tests per cell type
    cell_type_results = {}
    for ct in deconv_aligned.columns:
        scores = deconv_aligned[ct].values.astype(float)
        flt_scores = scores[y_aligned.values == 1]
        gc_scores = scores[y_aligned.values == 0]

        # Skip if all zeros or constant
        if np.std(scores) < 1e-10:
            cell_type_results[ct] = {
                "mean_flt": 0.0, "mean_gc": 0.0,
                "wilcoxon_p": 1.0, "cliffs_delta": 0.0,
                "direction": "none", "note": "constant scores"
            }
            continue

        try:
            stat, p = mannwhitneyu(flt_scores, gc_scores, alternative="two-sided")
        except ValueError:
            p = 1.0

        delta = cliffs_delta(flt_scores, gc_scores)
        # Direction follows Cliff's delta sign (rank-based, consistent with effect size)
        if abs(delta) < 1e-10:
            direction = "none"
        elif delta > 0:
            direction = "up_in_flight"
        else:
            direction = "down_in_flight"

        cell_type_results[ct] = {
            "mean_flt": round(float(np.mean(flt_scores)), 4),
            "mean_gc": round(float(np.mean(gc_scores)), 4),
            "std_flt": round(float(np.std(flt_scores)), 4),
            "std_gc": round(float(np.std(gc_scores)), 4),
            "wilcoxon_p": round(float(p), 6),
            "cliffs_delta": round(float(delta), 4),
            "direction": direction,
        }

    # BH-FDR correction
    pvals = [cell_type_results[ct]["wilcoxon_p"] for ct in cell_type_results
             if "note" not in cell_type_results[ct]]
    ct_names = [ct for ct in cell_type_results if "note" not in cell_type_results[ct]]

    if pvals:
        from scipy.stats import false_discovery_control
        try:
            fdr_p = false_discovery_control(pvals, method="bh")
        except AttributeError:
            # Older scipy: manual BH
            n_tests = len(pvals)
            sorted_idx = np.argsort(pvals)
            fdr_p = np.zeros(n_tests)
            for rank, idx in enumerate(sorted_idx):
                fdr_p[idx] = pvals[idx] * n_tests / (rank + 1)
            fdr_p = np.minimum.accumulate(fdr_p[np.argsort(sorted_idx)][::-1])[::-1]
            fdr_p = np.clip(fdr_p, 0, 1)

        for i, ct in enumerate(ct_names):
            cell_type_results[ct]["fdr_p"] = round(float(fdr_p[i]), 6)

    # Count significant
    n_sig = sum(1 for ct in cell_type_results
                if cell_type_results[ct].get("fdr_p", 1.0) < 0.05)
    print(f"  Significant cell types (FDR<0.05): {n_sig}/{len(cell_type_results)}")

    return {
        "tissue": tissue,
        "n_samples": len(common_samples),
        "n_flt": n_flt,
        "n_gc": n_gc,
        "n_cell_types": len(cell_type_results),
        "n_significant_fdr05": n_sig,
        "cell_types": cell_type_results,
    }


def main():
    parser = argparse.ArgumentParser(description="v5 Phase 1: Immune Deconvolution")
    parser.add_argument("--tissue", type=str, default=None,
                        help="Tissue to process (default: all 8)")
    args = parser.parse_args()

    V5_EVAL_DIR.mkdir(parents=True, exist_ok=True)
    symbol_map = load_symbol_map()
    print(f"Symbol map: {len(symbol_map)} entries")

    tissues = [args.tissue] if args.tissue else ALL_TISSUES

    for tissue in tissues:
        result = run_deconv_for_tissue(tissue, symbol_map)
        if result is None:
            print(f"  [SKIP] {tissue}: deconvolution failed")
            continue

        out_path = V5_EVAL_DIR / f"immune_deconv_{tissue}.json"
        with open(out_path, "w") as f:
            json.dump(result, f, indent=2)
        print(f"  Saved: {out_path}")

    print("\nDone!")


if __name__ == "__main__":
    main()
