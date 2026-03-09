#!/usr/bin/env python3
"""
dge_pipeline_comparison.py — GeneLab_benchmark: J2 DGE Pipeline Comparison Evaluation

Reads DGE results from 3 pipelines (DESeq2, edgeR, limma-voom) and computes:
  1. DEG Jaccard overlap at FDR < 0.05 and top-100/500
  2. Rank Spearman correlation of test statistics (all genes)
  3. Log2FC Pearson correlation
  4. GeneLab DGE concordance (our DESeq2 vs GeneLab's DESeq2)

Output: evaluation/J2_dge_pipeline_comparison.json

Usage:
    python scripts/dge_pipeline_comparison.py
    python scripts/dge_pipeline_comparison.py --dry-run
"""

import argparse
import json
import sys
import warnings
from datetime import datetime
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore", category=FutureWarning)

PROJECT_DIR = Path(__file__).resolve().parent.parent
DGE_DIR = PROJECT_DIR / "processed" / "J2_dge_comparison"
DATA_DIR = PROJECT_DIR / "data" / "mouse"
EVAL_DIR = PROJECT_DIR / "evaluation"

PIPELINES = ["deseq2", "edger", "limma_voom"]
PIPELINE_PAIRS = list(combinations(PIPELINES, 2))

TISSUE_MISSIONS = {
    "liver": [
        {"mission": "RR-1", "glds": "GLDS-48"},
        {"mission": "RR-3", "glds": "GLDS-137"},
        {"mission": "RR-6", "glds": "GLDS-245"},
        {"mission": "RR-8", "glds": "GLDS-379"},
        {"mission": "RR-9", "glds": "GLDS-242"},
        {"mission": "MHU-2", "glds": "GLDS-617"},
    ],
    "thymus": [
        {"mission": "RR-6", "glds": "GLDS-244"},
        {"mission": "MHU-2", "glds": "GLDS-289"},
        {"mission": "RR-9", "glds": "GLDS-421"},
    ],
}


class _NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.bool_):
            return bool(obj)
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj) if np.isfinite(obj) else None
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)


def load_dge(tissue, mission, pipeline):
    """Load a DGE result CSV."""
    path = DGE_DIR / tissue / f"{mission}_{pipeline}_dge.csv"
    if not path.exists():
        return None
    df = pd.read_csv(path)
    df = df.dropna(subset=["ENSEMBL"])
    df = df.set_index("ENSEMBL")
    return df


def jaccard(set_a, set_b):
    """Compute Jaccard similarity."""
    if len(set_a) == 0 and len(set_b) == 0:
        return 1.0
    intersection = len(set_a & set_b)
    union = len(set_a | set_b)
    return intersection / union if union > 0 else 0.0


def compare_pair(df_a, df_b, name_a, name_b):
    """Compare two DGE results."""
    # Align to common genes
    common = df_a.index.intersection(df_b.index)
    a = df_a.loc[common]
    b = df_b.loc[common]

    result = {
        "pair": f"{name_a}_vs_{name_b}",
        "n_common_genes": len(common),
    }

    # 1. DEG Jaccard at FDR < 0.05
    deg_a = set(a.index[a["adj_pvalue"] < 0.05])
    deg_b = set(b.index[b["adj_pvalue"] < 0.05])
    result["deg_jaccard_005"] = jaccard(deg_a, deg_b)
    result["n_deg_a"] = len(deg_a)
    result["n_deg_b"] = len(deg_b)
    result["n_deg_intersection"] = len(deg_a & deg_b)

    # 2. Top-k Jaccard (by p-value rank)
    for k in [100, 500]:
        top_a = set(a.nsmallest(k, "pvalue").index)
        top_b = set(b.nsmallest(k, "pvalue").index)
        result[f"top{k}_jaccard"] = jaccard(top_a, top_b)

    # 3. Log2FC Spearman rank correlation (universal across all pipelines)
    # Note: edgeR F-statistic is unsigned, so we use log2FC for rank comparison
    mask = a["log2FC"].notna() & b["log2FC"].notna()
    if mask.sum() > 10:
        rho, p = stats.spearmanr(a.loc[mask, "log2FC"], b.loc[mask, "log2FC"])
        result["log2fc_spearman_rho"] = float(rho)
        result["log2fc_spearman_p"] = float(p)

    # 4. Log2FC Pearson correlation
    mask_fc = a["log2FC"].notna() & b["log2FC"].notna()
    if mask_fc.sum() > 10:
        r, p = stats.pearsonr(a.loc[mask_fc, "log2FC"], b.loc[mask_fc, "log2FC"])
        result["log2fc_pearson_r"] = float(r)
        result["log2fc_pearson_p"] = float(p)

    # 5. Direction concordance (sign of log2FC among shared DEGs)
    shared_degs = deg_a & deg_b
    if len(shared_degs) > 0:
        signs_a = np.sign(a.loc[list(shared_degs), "log2FC"])
        signs_b = np.sign(b.loc[list(shared_degs), "log2FC"])
        result["direction_concordance"] = float((signs_a == signs_b).mean())
    else:
        result["direction_concordance"] = None

    return result


def compare_with_genelab(tissue, mission_info, our_deseq2):
    """Compare our DESeq2 replication vs GeneLab's original DGE."""
    mission = mission_info["mission"]
    glds = mission_info["glds"]

    # Find GeneLab DGE file
    mission_dir = mission.replace("MHU-2", "MHU-2")
    dge_patterns = [
        DATA_DIR / tissue / mission_dir / f"{glds}_rna_seq_differential_expression_GLbulkRNAseq.csv",
        DATA_DIR / tissue / mission_dir / f"{glds}_rna_seq_differential_expression.csv",
    ]

    gl_path = None
    for p in dge_patterns:
        if p.exists():
            gl_path = p
            break

    if gl_path is None:
        return None

    gl_dge = pd.read_csv(gl_path)

    # Find the right contrast column (Stat_ with Space Flight vs Ground Control)
    stat_cols = [c for c in gl_dge.columns if c.startswith("Stat_")]
    if not stat_cols:
        return None

    # Use detect_flight_contrast logic (simplified)
    target_col = None
    # Priority 1: GC-first with Space Flight
    for col in stat_cols:
        if "Ground" in col and "Space" in col and "Ground" in col.split("v")[0]:
            target_col = col
            break
    # Priority 2: Any with Flight and Control
    if target_col is None:
        for col in stat_cols:
            if ("Flight" in col or "FLT" in col) and ("Control" in col or "GC" in col):
                target_col = col
                break
    if target_col is None:
        return None

    # Extract corresponding log2fc and adj.p.value columns
    contrast_suffix = target_col.replace("Stat_", "")
    log2fc_col = f"Log2fc_{contrast_suffix}"
    adjp_col = f"Adj.p.value_{contrast_suffix}"

    if log2fc_col not in gl_dge.columns or adjp_col not in gl_dge.columns:
        return None

    # Build comparable dataframe
    gl_df = pd.DataFrame({
        "ENSEMBL": gl_dge["ENSEMBL"],
        "log2FC": gl_dge[log2fc_col],
        "stat": gl_dge[target_col],
        "adj_pvalue": gl_dge[adjp_col],
        "pvalue": gl_dge.get(f"P.value_{contrast_suffix}", np.nan),
    }).set_index("ENSEMBL").dropna(subset=["stat"])

    # Our DESeq2 result
    our_df = our_deseq2.dropna(subset=["stat"])

    common = gl_df.index.intersection(our_df.index)
    if len(common) < 100:
        return None

    gl = gl_df.loc[common]
    ours = our_df.loc[common]

    # GeneLab contrast is GC-first (positive = GC > FLT),
    # our DESeq2 is FLT vs GC (positive = FLT > GC).
    # Determine sign convention from contrast name.
    is_gc_first = "Ground" in target_col.split("v")[0]
    sign_flip = -1.0 if is_gc_first else 1.0

    # Spearman of log2FC (sign-corrected)
    mask = gl["log2FC"].notna() & ours["log2FC"].notna()
    rho, p = stats.spearmanr(gl.loc[mask, "log2FC"] * sign_flip, ours.loc[mask, "log2FC"])

    # Log2FC Pearson (sign-corrected)
    r_fc, _ = stats.pearsonr(gl.loc[mask, "log2FC"] * sign_flip, ours.loc[mask, "log2FC"])

    # DEG Jaccard
    deg_gl = set(gl.index[gl["adj_pvalue"] < 0.05])
    deg_ours = set(ours.index[ours["adj_pvalue"] < 0.05])
    jac = jaccard(deg_gl, deg_ours)

    # Top-100 Jaccard
    top_gl = set(gl.nsmallest(100, "pvalue").index)
    top_ours = set(ours.nsmallest(100, "pvalue").index)
    jac_top100 = jaccard(top_gl, top_ours)

    return {
        "genelab_contrast": target_col,
        "sign_corrected": is_gc_first,
        "n_common_genes": len(common),
        "log2fc_spearman_rho": float(rho),
        "log2fc_pearson_r": float(r_fc),
        "deg_jaccard_005": float(jac),
        "top100_jaccard": float(jac_top100),
        "n_deg_genelab": len(deg_gl),
        "n_deg_ours": len(deg_ours),
        "note": "Our DESeq2 (FLT vs GC binary) vs GeneLab DESeq2 (original multi-level contrast). Sign-corrected for direction convention.",
    }


def run_tissue(tissue):
    """Run all comparisons for a tissue."""
    results = {}

    for mi in TISSUE_MISSIONS[tissue]:
        mission = mi["mission"]
        print(f"\n--- {tissue} / {mission} ---")

        # Load all pipelines
        dge = {}
        for p in PIPELINES:
            df = load_dge(tissue, mission, p)
            if df is not None:
                dge[p] = df
                print(f"  {p}: {len(df)} genes, {(df['adj_pvalue'] < 0.05).sum()} DEGs")
            else:
                print(f"  {p}: NOT FOUND")

        if len(dge) < 2:
            print("  Skipping — fewer than 2 pipelines available")
            continue

        # Pairwise comparisons
        comparisons = []
        for p_a, p_b in PIPELINE_PAIRS:
            if p_a in dge and p_b in dge:
                comp = compare_pair(dge[p_a], dge[p_b], p_a, p_b)
                comparisons.append(comp)
                print(f"  {p_a} vs {p_b}: Jaccard(FDR<0.05)={comp['deg_jaccard_005']:.3f}, "
                      f"log2FC Spearman={comp.get('log2fc_spearman_rho', 'N/A'):.3f}, "
                      f"log2FC Pearson={comp.get('log2fc_pearson_r', 'N/A'):.3f}")

        # GeneLab concordance (our DESeq2 vs original)
        gl_concordance = None
        if "deseq2" in dge:
            gl_concordance = compare_with_genelab(tissue, mi, dge["deseq2"])
            if gl_concordance:
                print(f"  GeneLab concordance: log2FC Spearman={gl_concordance['log2fc_spearman_rho']:.3f}, "
                      f"Jaccard(FDR<0.05)={gl_concordance['deg_jaccard_005']:.3f}")

        results[mission] = {
            "tissue": tissue,
            "mission": mission,
            "n_pipelines": len(dge),
            "per_pipeline": {
                p: {"n_genes": len(df), "n_deg_005": int((df["adj_pvalue"] < 0.05).sum())}
                for p, df in dge.items()
            },
            "pairwise_comparisons": comparisons,
            "genelab_concordance": gl_concordance,
        }

    return results


def compute_summary(all_results):
    """Compute aggregate summary statistics."""
    all_spearman = []
    all_pearson = []
    all_jaccard = []
    all_top100 = []
    all_concordance = []
    gl_spearman = []

    for tissue_results in all_results.values():
        for mission_data in tissue_results.values():
            for comp in mission_data.get("pairwise_comparisons", []):
                if "log2fc_spearman_rho" in comp:
                    all_spearman.append(comp["log2fc_spearman_rho"])
                if "log2fc_pearson_r" in comp:
                    all_pearson.append(comp["log2fc_pearson_r"])
                if "deg_jaccard_005" in comp:
                    all_jaccard.append(comp["deg_jaccard_005"])
                if "top100_jaccard" in comp:
                    all_top100.append(comp["top100_jaccard"])
                if comp.get("direction_concordance") is not None:
                    all_concordance.append(comp["direction_concordance"])

            gl = mission_data.get("genelab_concordance")
            if gl and "log2fc_spearman_rho" in gl:
                gl_spearman.append(gl["log2fc_spearman_rho"])

    summary = {
        "n_missions": sum(len(v) for v in all_results.values()),
        "n_comparisons": len(all_spearman),
    }

    if all_spearman:
        summary["log2fc_spearman"] = {
            "mean": float(np.mean(all_spearman)),
            "min": float(np.min(all_spearman)),
            "max": float(np.max(all_spearman)),
        }
    if all_pearson:
        summary["log2fc_pearson"] = {
            "mean": float(np.mean(all_pearson)),
            "min": float(np.min(all_pearson)),
            "max": float(np.max(all_pearson)),
        }
    if all_jaccard:
        summary["deg_jaccard_005"] = {
            "mean": float(np.mean(all_jaccard)),
            "min": float(np.min(all_jaccard)),
            "max": float(np.max(all_jaccard)),
        }
    if all_top100:
        summary["top100_jaccard"] = {
            "mean": float(np.mean(all_top100)),
            "min": float(np.min(all_top100)),
            "max": float(np.max(all_top100)),
        }
    if all_concordance:
        summary["direction_concordance"] = {
            "mean": float(np.mean(all_concordance)),
            "min": float(np.min(all_concordance)),
        }
    if gl_spearman:
        summary["genelab_replication"] = {
            "mean_spearman": float(np.mean(gl_spearman)),
            "min_spearman": float(np.min(gl_spearman)),
            "n_missions": len(gl_spearman),
        }

    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true",
                        help="Print summary without writing files")
    args = parser.parse_args()

    print("=== J2 DGE Pipeline Comparison ===\n")

    all_results = {}
    for tissue in TISSUE_MISSIONS:
        print(f"\n{'='*50}")
        print(f"Tissue: {tissue}")
        print(f"{'='*50}")
        all_results[tissue] = run_tissue(tissue)

    summary = compute_summary(all_results)

    # Determine verdict
    mean_spearman = summary.get("log2fc_spearman", {}).get("mean", 0)
    mean_pearson = summary.get("log2fc_pearson", {}).get("mean", 0)
    mean_jaccard = summary.get("deg_jaccard_005", {}).get("mean", 0)

    if mean_spearman > 0.90 and mean_pearson > 0.90:
        verdict = "DGE pipeline choice has minimal impact: fold-change rankings highly conserved"
    elif mean_spearman > 0.70:
        verdict = "DGE pipeline choice has moderate impact: rankings consistent, DEG lists vary by stringency"
    else:
        verdict = "DGE pipeline choice significantly affects results"

    output = {
        "timestamp": datetime.now().isoformat(),
        "description": "J2: DGE pipeline comparison (DESeq2 vs edgeR vs limma-voom)",
        "design_decisions": ["DD-21"],
        "method": "DESeq2 Wald, edgeR QLF, limma-voom QW; FLT vs GC contrast",
        "scope": "liver (6 missions) + thymus (3 missions)",
        "pipelines": PIPELINES,
        "results": all_results,
        "summary": summary,
        "verdict": verdict,
        "evidence": (
            f"Log2FC Spearman: mean={mean_spearman:.3f}, "
            f"Log2FC Pearson: mean={mean_pearson:.3f}, "
            f"DEG Jaccard(FDR<0.05): mean={mean_jaccard:.3f}. "
            f"GeneLab replication: {summary.get('genelab_replication', {}).get('mean_spearman', 'N/A')}"
        ),
    }

    print(f"\n{'='*50}")
    print("SUMMARY")
    print(f"{'='*50}")
    print(f"  Missions: {summary['n_missions']}")
    print(f"  Comparisons: {summary['n_comparisons']}")
    if "log2fc_spearman" in summary:
        s = summary["log2fc_spearman"]
        print(f"  Log2FC Spearman: mean={s['mean']:.3f} [{s['min']:.3f}, {s['max']:.3f}]")
    if "log2fc_pearson" in summary:
        s = summary["log2fc_pearson"]
        print(f"  Log2FC Pearson: mean={s['mean']:.3f} [{s['min']:.3f}, {s['max']:.3f}]")
    if "deg_jaccard_005" in summary:
        s = summary["deg_jaccard_005"]
        print(f"  DEG Jaccard(FDR<0.05): mean={s['mean']:.3f} [{s['min']:.3f}, {s['max']:.3f}]")
    if "genelab_replication" in summary:
        s = summary["genelab_replication"]
        print(f"  GeneLab replication: mean log2FC Spearman={s['mean_spearman']:.3f} ({s['n_missions']} missions)")
    print(f"\n  Verdict: {verdict}")

    if args.dry_run:
        print("\n[DRY RUN] Would write to:")
        print(f"  {EVAL_DIR / 'J2_dge_pipeline_comparison.json'}")
        return

    EVAL_DIR.mkdir(parents=True, exist_ok=True)
    out_path = EVAL_DIR / "J2_dge_pipeline_comparison.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, cls=_NumpyEncoder)
    print(f"\nWrote {out_path} ({out_path.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
