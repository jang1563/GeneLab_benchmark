#!/usr/bin/env python3
"""
preprocess_pathways.py — GeneLab_benchmark: Pathway Feature Integration (DD-15)

Integrates fGSEA and GSVA results into the benchmark pipeline.
Provides gene-level vs pathway-level feature comparison across ALL categories.

Functions:
  1. Category C Method C: pathway-level cross-tissue transfer features
  2. Biological validation: SHAP top genes → pathway enrichment
  3. Cross-mission NES conservation analysis
  4. Gene vs pathway feature loading for generate_tasks.py

Usage:
  python scripts/preprocess_pathways.py --mode cross-tissue --train liver --test kidney
  python scripts/preprocess_pathways.py --mode conservation --tissue liver
  python scripts/preprocess_pathways.py --mode validate-shap --task A1
  python scripts/preprocess_pathways.py --mode merge-scores --tissue liver --db hallmark
"""

import json
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats as sp_stats

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent.parent
FGSEA_DIR = BASE_DIR / "processed" / "fgsea"
PATHWAY_SCORES_DIR = BASE_DIR / "processed" / "pathway_scores"
A_DETECTION_DIR = BASE_DIR / "processed" / "A_detection"
C_CROSS_TISSUE_DIR = BASE_DIR / "processed" / "C_cross_tissue"
EVALUATION_DIR = BASE_DIR / "evaluation"

# ── Tissue-Mission mapping ─────────────────────────────────────────────────────
TISSUE_MISSIONS = {
    "liver": ["RR-1", "RR-3", "RR-6", "RR-8", "RR-9", "MHU-2"],
    "kidney": ["RR-1", "RR-3", "RR-7"],
    "thymus": ["RR-6", "MHU-2", "RR-9"],
    "gastrocnemius": ["RR-1", "RR-5", "RR-9"],
    "eye": ["RR-1", "RR-3", "TBD"],
}

# ── Category C cross-tissue pairs ─────────────────────────────────────────────
CROSS_TISSUE_TASKS = {
    "C1": {"train": "liver", "test": "kidney", "biology": "metabolism, oxidative stress"},
    "C2": {"train": "liver", "test": "gastrocnemius", "biology": "energy metabolism"},
    "C3": {"train": "liver", "test": "thymus", "biology": "immune (Kupffer cell)"},
    "C4": {"train": "thymus", "test": "kidney", "biology": "immune-renal"},
    "C5": {"train": "all", "test": "holdout", "biology": "universal signature"},
}


# ── Utility Functions ──────────────────────────────────────────────────────────

def load_fgsea_results(tissue, mission, db="hallmark"):
    """Load fGSEA results for a single tissue/mission/db."""
    f = FGSEA_DIR / tissue / f"{mission}_fgsea_{db}.csv"
    if not f.exists():
        return None
    return pd.read_csv(f)


def load_gsva_scores(tissue, mission, db="hallmark", method="gsva"):
    """Load GSVA pathway scores for a single tissue/mission/db."""
    f = PATHWAY_SCORES_DIR / tissue / f"{mission}_{method}_{db}.csv"
    if not f.exists():
        return None
    return pd.read_csv(f, index_col=0)


def load_metadata(tissue, mission):
    """Load sample metadata for a tissue/mission."""
    f = A_DETECTION_DIR / tissue / f"{tissue}_{mission}_metadata.csv"
    if not f.exists():
        return None
    return pd.read_csv(f)


def load_all_fgsea(tissue, db="hallmark"):
    """Load and concatenate fGSEA results across all missions for a tissue."""
    dfs = []
    for mission in TISSUE_MISSIONS.get(tissue, []):
        df = load_fgsea_results(tissue, mission, db)
        if df is not None:
            dfs.append(df)
    if not dfs:
        return None
    return pd.concat(dfs, ignore_index=True)


def merge_gsva_scores(tissue, db="hallmark", method="gsva"):
    """Merge GSVA scores across missions for a tissue, adding mission/label columns."""
    all_scores = []
    all_meta = []

    for mission in TISSUE_MISSIONS.get(tissue, []):
        scores = load_gsva_scores(tissue, mission, db, method)
        meta = load_metadata(tissue, mission)

        if scores is None or meta is None:
            continue

        scores["mission"] = mission
        # Align metadata with scores by sample name
        if "sample" in meta.columns:
            meta = meta.set_index("sample")
        scores_with_meta = scores.copy()

        all_scores.append(scores_with_meta)

    if not all_scores:
        return None
    return pd.concat(all_scores, ignore_index=True)


# ── Mode 1: Cross-Tissue Transfer (Category C Method C) ───────────────────────

def build_cross_tissue_features(train_tissue, test_tissue, db="hallmark",
                                 top_n=20, method="gsva"):
    """
    Category C Method C: Pathway-level cross-tissue transfer.

    1. Select top-N pathways from train tissue fGSEA (by |NES|)
    2. Load GSVA scores for both tissues using those pathways
    3. Save as ML-ready feature matrices

    LOMO-aware: pathway selection from train missions only.
    """
    print(f"\n=== Category C Method C: {train_tissue} → {test_tissue} ===")

    # Step 1: Select top pathways from train tissue
    fgsea_all = load_all_fgsea(train_tissue, db)
    if fgsea_all is None or fgsea_all.empty:
        print(f"  [ERROR] No fGSEA results for {train_tissue}/{db}")
        return None

    # Aggregate NES across missions: mean |NES|, count significant
    pathway_agg = (fgsea_all
                   .groupby("pathway")
                   .agg(
                       mean_NES=("NES", "mean"),
                       abs_mean_NES=("NES", lambda x: abs(x.mean())),
                       std_NES=("NES", "std"),
                       n_missions=("mission", "nunique"),
                       n_significant=("padj", lambda x: (x < 0.05).sum()),
                       mean_padj=("padj", "mean"),
                   )
                   .sort_values("abs_mean_NES", ascending=False))

    selected = pathway_agg.head(top_n).index.tolist()
    print(f"  Selected {len(selected)} pathways (top by |NES|)")

    for i, pw in enumerate(selected[:5]):
        row = pathway_agg.loc[pw]
        print(f"    {i+1}. {pw[:50]:50s} NES={row['mean_NES']:+.2f} "
              f"(sig in {row['n_significant']:.0f}/{row['n_missions']:.0f} missions)")

    # Step 2-3: Load GSVA scores for selected pathways
    train_scores = merge_gsva_scores(train_tissue, db, method)
    test_scores = merge_gsva_scores(test_tissue, db, method)

    if train_scores is None or test_scores is None:
        print(f"  [ERROR] Missing GSVA scores")
        return None

    # Filter to selected pathways (columns)
    available_train = [p for p in selected if p in train_scores.columns]
    available_test = [p for p in selected if p in test_scores.columns]
    common = list(set(available_train) & set(available_test))

    if len(common) == 0:
        print("  [ERROR] No common pathways between train and test tissues")
        return None

    print(f"  Common pathways: {len(common)} / {top_n}")

    # Save
    C_CROSS_TISSUE_DIR.mkdir(parents=True, exist_ok=True)

    task_id = [k for k, v in CROSS_TISSUE_TASKS.items()
               if v["train"] == train_tissue and v["test"] == test_tissue]
    task_id = task_id[0] if task_id else f"{train_tissue}_{test_tissue}"

    # Save pathway selection
    pathway_agg.loc[selected].to_csv(
        C_CROSS_TISSUE_DIR / f"{task_id}_pathway_selection_{db}.csv")

    # Save feature matrices
    meta_cols = ["mission"]
    feature_cols = common

    train_out = train_scores[meta_cols + feature_cols].copy()
    test_out = test_scores[meta_cols + feature_cols].copy()

    train_out.to_csv(C_CROSS_TISSUE_DIR / f"{task_id}_train_{train_tissue}_{db}.csv",
                     index=False)
    test_out.to_csv(C_CROSS_TISSUE_DIR / f"{task_id}_test_{test_tissue}_{db}.csv",
                    index=False)

    print(f"  Saved: {task_id}_train/test_{db}.csv")
    print(f"  Train: {len(train_out)} samples × {len(feature_cols)} pathways")
    print(f"  Test:  {len(test_out)} samples × {len(feature_cols)} pathways")

    return {"selected_pathways": common, "train_shape": train_out.shape,
            "test_shape": test_out.shape}


# ── Mode 2: Cross-Mission NES Conservation ─────────────────────────────────────

def cross_mission_conservation(tissue, db="hallmark"):
    """
    Compare NES profiles across missions within a tissue.
    Output: Spearman correlation matrix (missions × missions).
    """
    print(f"\n=== NES Conservation: {tissue} / {db} ===")

    nes_profiles = {}
    for mission in TISSUE_MISSIONS.get(tissue, []):
        df = load_fgsea_results(tissue, mission, db)
        if df is not None:
            nes_profiles[mission] = df.set_index("pathway")["NES"]

    if len(nes_profiles) < 2:
        print(f"  [ERROR] Need ≥2 missions, got {len(nes_profiles)}")
        return None

    # Build NES matrix (pathways × missions)
    nes_matrix = pd.DataFrame(nes_profiles)
    nes_matrix = nes_matrix.dropna()
    print(f"  NES matrix: {nes_matrix.shape[0]} pathways × {nes_matrix.shape[1]} missions")

    # Spearman correlation
    corr_matrix = nes_matrix.corr(method="spearman")

    # P-values
    missions_list = list(nes_profiles.keys())
    n = len(missions_list)
    pval_matrix = pd.DataFrame(np.ones((n, n)), index=missions_list, columns=missions_list)
    for i in range(n):
        for j in range(i + 1, n):
            common = nes_matrix[[missions_list[i], missions_list[j]]].dropna()
            if len(common) >= 3:
                r, p = sp_stats.spearmanr(common.iloc[:, 0], common.iloc[:, 1])
                pval_matrix.iloc[i, j] = p
                pval_matrix.iloc[j, i] = p

    # Save
    out_dir = FGSEA_DIR / "summary"
    out_dir.mkdir(parents=True, exist_ok=True)

    corr_matrix.to_csv(out_dir / f"nes_correlation_{tissue}_{db}.csv")
    pval_matrix.to_csv(out_dir / f"nes_pvalue_{tissue}_{db}.csv")

    print(f"\n  Spearman correlation matrix:")
    print(corr_matrix.round(3).to_string())

    # Mean off-diagonal correlation
    mask = np.ones(corr_matrix.shape, dtype=bool)
    np.fill_diagonal(mask, False)
    mean_r = corr_matrix.values[mask].mean()
    print(f"\n  Mean off-diagonal Spearman r: {mean_r:.3f}")

    return corr_matrix


# ── Mode 3: SHAP Pathway Validation ───────────────────────────────────────────

def validate_shap_pathways(task_id, db="hallmark", top_n=100):
    """
    Check if SHAP top genes from Category A models are enriched in known pathways.
    Uses hypergeometric test (not fGSEA, since input is gene list not ranked).
    """
    from scipy.stats import hypergeom

    print(f"\n=== SHAP Pathway Validation: {task_id} / {db} ===")

    # Load SHAP results
    shap_file = EVALUATION_DIR / f"{task_id}_shap_importance.csv"
    if not shap_file.exists():
        print(f"  [SKIP] No SHAP file: {shap_file}")
        return None

    shap_df = pd.read_csv(shap_file)
    shap_genes = set(shap_df.head(top_n)["gene"].tolist())
    print(f"  SHAP top {top_n} genes loaded")

    # Load gene sets from fGSEA summary
    # Use msigdbr gene sets via the R-generated fGSEA results
    summary_file = FGSEA_DIR / "summary" / f"all_fgsea_{db}.csv"
    if not summary_file.exists():
        print(f"  [SKIP] No fGSEA summary: {summary_file}")
        return None

    fgsea_summary = pd.read_csv(summary_file)

    # For each pathway with leadingEdge, test enrichment in SHAP genes
    # Use universe = all genes in the model
    all_genes_file = list((A_DETECTION_DIR / task_id.split("_")[0].lower()).glob("*all_missions_log2_norm.csv"))
    if all_genes_file:
        all_genes_df = pd.read_csv(all_genes_file[0], index_col=0, nrows=1)
        universe_size = len(all_genes_df.columns)
    else:
        universe_size = 20000  # fallback

    results = []
    pathways = fgsea_summary.drop_duplicates("pathway")

    for _, row in pathways.iterrows():
        if pd.isna(row.get("leadingEdge_str", np.nan)):
            continue
        le_genes = set(str(row["leadingEdge_str"]).split("; "))
        overlap = shap_genes & le_genes

        if len(overlap) > 0:
            # Hypergeometric test
            # P(X >= |overlap|) where X ~ Hypergeom(N=universe, K=|pathway|, n=|shap|)
            p_val = hypergeom.sf(len(overlap) - 1, universe_size,
                                 len(le_genes), len(shap_genes))
            results.append({
                "pathway": row["pathway"],
                "n_overlap": len(overlap),
                "overlap_genes": "; ".join(sorted(overlap)),
                "pathway_size": len(le_genes),
                "shap_size": len(shap_genes),
                "hypergeom_p": p_val,
                "fgsea_NES": row.get("NES", np.nan),
            })

    if not results:
        print("  No overlapping pathways found")
        return None

    results_df = pd.DataFrame(results).sort_values("hypergeom_p")

    # BH FDR correction
    from scipy.stats import false_discovery_control
    try:
        results_df["hypergeom_q"] = false_discovery_control(results_df["hypergeom_p"])
    except (AttributeError, TypeError):
        # scipy < 1.12 fallback
        n = len(results_df)
        ranks = np.arange(1, n + 1)
        results_df["hypergeom_q"] = results_df["hypergeom_p"] * n / ranks

    # Save
    out_dir = FGSEA_DIR / "validation"
    out_dir.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(out_dir / f"{task_id}_shap_pathway_{db}.csv", index=False)

    print(f"  {len(results_df)} pathways with SHAP overlap")
    sig = results_df[results_df["hypergeom_q"] < 0.05]
    print(f"  {len(sig)} significant (FDR < 0.05)")
    if not sig.empty:
        for _, r in sig.head(5).iterrows():
            print(f"    {r['pathway'][:50]:50s}  overlap={r['n_overlap']}  q={r['hypergeom_q']:.1e}")

    return results_df


# ── Mode 4: Merge Scores for generate_tasks.py ────────────────────────────────

def merge_all_scores(tissue, db="hallmark", method="gsva"):
    """
    Merge GSVA scores across missions into a single file for generate_tasks.py.
    Output: processed/pathway_scores/{tissue}/{tissue}_all_missions_{method}_{db}.csv
    """
    print(f"\n=== Merge Scores: {tissue} / {db} / {method} ===")

    merged = merge_gsva_scores(tissue, db, method)
    if merged is None:
        print("  [ERROR] No scores to merge")
        return None

    out_file = PATHWAY_SCORES_DIR / tissue / f"{tissue}_all_missions_{method}_{db}.csv"
    merged.to_csv(out_file, index=False)
    print(f"  Saved: {out_file.name}")
    print(f"  Shape: {merged.shape}")

    return merged


# ── CLI ────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Pathway feature integration for GeneLab_benchmark")
    parser.add_argument("--mode", required=True,
                        choices=["cross-tissue", "conservation", "validate-shap",
                                 "merge-scores", "all"],
                        help="Processing mode")
    parser.add_argument("--tissue", help="Tissue name")
    parser.add_argument("--train", help="Training tissue (for cross-tissue)")
    parser.add_argument("--test", help="Test tissue (for cross-tissue)")
    parser.add_argument("--task", help="Task ID (for validate-shap, e.g., A1)")
    parser.add_argument("--db", default="hallmark",
                        help="Gene set DB (hallmark, kegg, reactome)")
    parser.add_argument("--top-n", type=int, default=20,
                        help="Number of top pathways to select (default: 20)")
    parser.add_argument("--method", default="gsva",
                        help="Scoring method: gsva or ssgsea")
    args = parser.parse_args()

    if args.mode == "cross-tissue":
        if not args.train or not args.test:
            parser.error("--train and --test required for cross-tissue mode")
        build_cross_tissue_features(args.train, args.test, args.db,
                                     args.top_n, args.method)

    elif args.mode == "conservation":
        if not args.tissue:
            parser.error("--tissue required for conservation mode")
        cross_mission_conservation(args.tissue, args.db)

    elif args.mode == "validate-shap":
        if not args.task:
            parser.error("--task required for validate-shap mode")
        validate_shap_pathways(args.task, args.db)

    elif args.mode == "merge-scores":
        if not args.tissue:
            parser.error("--tissue required for merge-scores mode")
        merge_all_scores(args.tissue, args.db, args.method)

    elif args.mode == "all":
        # Run conservation for all tissues
        for tissue in TISSUE_MISSIONS:
            cross_mission_conservation(tissue, args.db)

        # Run cross-tissue for defined pairs
        for task_id, task_def in CROSS_TISSUE_TASKS.items():
            if task_def["train"] != "all":
                build_cross_tissue_features(
                    task_def["train"], task_def["test"],
                    args.db, args.top_n, args.method)


if __name__ == "__main__":
    main()
