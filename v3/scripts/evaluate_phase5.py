#!/usr/bin/env python3
"""
evaluate_phase5.py — Phase 5 Radiation × Microgravity (HLU) analysis.

Tasks:
  R1: Pairwise classification — Radiation vs Ctrl, HLU vs Ctrl, Combined vs Ctrl
      → PCA-LR AUROC (LOO-CV due to small n)
  R2: Gene-level DGE concordance — ground stressor t-statistics vs ISS flight t-statistics
      (tissue-matched: skin_rad ↔ skin)
  R3: Combined vs single-stressor concordance with ISS
      → r(combined DGE, ISS) vs r(rad-only DGE, ISS) vs r(HLU-only DGE, ISS)

Datasets:
  brain (OSD-202): 2×2×2 (radiation × HLU × time), n=40
  spleen (OSD-211): 2×2, n=21
  skin_rad (OSD-237): 2×2, n=21

Output:
  v3/evaluation/R1_radiation_classification.json
  v3/evaluation/R2_radiation_dge_concordance.json
  v3/evaluation/R3_combined_vs_single.json
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler

BASE_DIR = Path(__file__).resolve().parent.parent.parent
RAD_DIR = BASE_DIR / "processed" / "R_radiation"
SKIN_DIR = BASE_DIR / "processed" / "A_detection" / "skin"
OUTPUT_DIR = BASE_DIR / "v3" / "evaluation"

VARIANCE_PERCENTILE = 0.25
N_PCA = 20
N_BOOTSTRAP = 2000
N_PERM = 1000
SEED = 42

TISSUES = ["brain", "spleen", "skin_rad"]

# v1 skin missions for ISS comparison
V1_SKIN_MISSIONS = ["MHU-2_(dorsal)", "MHU-2_(femoral)", "RR-6", "RR-7"]

# Ground labels that map to y=0
GROUND_LABELS = {"Ground", "GC", "Ground Control"}


def load_radiation_data(tissue):
    """Load preprocessed radiation data."""
    log2 = pd.read_csv(RAD_DIR / tissue / f"{tissue}_log2_norm.csv", index_col=0)
    meta = pd.read_csv(RAD_DIR / tissue / f"{tissue}_metadata.csv", index_col=0)

    if log2.shape[0] < log2.shape[1]:
        pass
    else:
        log2 = log2.T

    ensmusg_cols = [c for c in log2.columns if str(c).startswith("ENSMUSG")]
    log2 = log2[ensmusg_cols]
    log2 = log2.apply(pd.to_numeric, errors="coerce").fillna(0)

    common = sorted(set(log2.index) & set(meta.index))
    return log2.loc[common], meta.loc[common]


def load_v1_skin_mission(mission):
    """Load v1 skin data for a specific mission."""
    log2_path = SKIN_DIR / f"skin_{mission}_log2_norm.csv"
    meta_path = SKIN_DIR / f"skin_{mission}_metadata.csv"
    if not log2_path.exists():
        return None, None

    log2 = pd.read_csv(log2_path, index_col=0)
    meta = pd.read_csv(meta_path, index_col=0)

    if log2.shape[0] < log2.shape[1]:
        pass
    else:
        log2 = log2.T

    ensmusg_cols = [c for c in log2.columns if str(c).startswith("ENSMUSG")]
    log2 = log2[ensmusg_cols]
    log2 = log2.apply(pd.to_numeric, errors="coerce").fillna(0)

    common = sorted(set(log2.index) & set(meta.index))
    return log2.loc[common], meta.loc[common]


def bootstrap_auroc(y_true, y_score, n_boot=N_BOOTSTRAP, seed=SEED):
    """Bootstrap 95% CI for AUROC."""
    rng = np.random.default_rng(seed)
    n = len(y_true)
    scores = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if len(np.unique(y_true[idx])) < 2:
            continue
        scores.append(roc_auc_score(y_true[idx], y_score[idx]))
    if not scores:
        return (float("nan"), float("nan"))
    return float(np.percentile(scores, 2.5)), float(np.percentile(scores, 97.5))


def permutation_p(y_true, y_score, n_perm=N_PERM, seed=SEED):
    """Permutation test p-value."""
    rng = np.random.default_rng(seed)
    observed = roc_auc_score(y_true, y_score)
    null_scores = []
    for _ in range(n_perm):
        y_perm = rng.permutation(y_true)
        if len(np.unique(y_perm)) < 2:
            continue
        null_scores.append(roc_auc_score(y_perm, y_score))
    null_arr = np.array(null_scores)
    p = float((np.sum(null_arr >= observed) + 1) / (len(null_arr) + 1))
    return p


def classify_pairwise(X, y, name=""):
    """LOO-CV classification with PCA-LR. Returns AUROC + CI + p."""
    if len(np.unique(y)) < 2:
        print(f"    {name}: SKIP (single class)")
        return None
    if len(y) < 6:
        print(f"    {name}: SKIP (n={len(y)} < 6)")
        return None

    loo = LeaveOneOut()
    y_true_all = []
    y_score_all = []

    for train_idx, test_idx in loo.split(X, y):
        train_X, test_X = X.iloc[train_idx], X.iloc[test_idx]
        train_y, test_y = y[train_idx], y[test_idx]

        gene_var = train_X.var(axis=0)
        threshold = gene_var.quantile(VARIANCE_PERCENTILE)
        selected = gene_var[gene_var >= threshold].index.tolist()
        train_X_f = train_X[selected]
        test_X_f = test_X[selected]

        scaler = StandardScaler()
        train_X_s = scaler.fit_transform(train_X_f)
        test_X_s = scaler.transform(test_X_f)

        n_comp = min(N_PCA, train_X_s.shape[0] - 1, train_X_s.shape[1])
        if n_comp < 1:
            continue
        pca = PCA(n_components=n_comp, random_state=SEED)
        train_X_p = pca.fit_transform(train_X_s)
        test_X_p = pca.transform(test_X_s)

        lr = LogisticRegression(max_iter=1000, random_state=SEED)
        lr.fit(train_X_p, train_y)
        y_score_all.append(lr.predict_proba(test_X_p)[:, 1][0])
        y_true_all.append(test_y[0])

    y_true_arr = np.array(y_true_all)
    y_score_arr = np.array(y_score_all)

    if len(np.unique(y_true_arr)) < 2:
        print(f"    {name}: SKIP (single class in predictions)")
        return None

    auroc = roc_auc_score(y_true_arr, y_score_arr)
    ci_lo, ci_hi = bootstrap_auroc(y_true_arr, y_score_arr)
    perm_p = permutation_p(y_true_arr, y_score_arr)

    print(f"    {name}: AUROC={auroc:.3f} [{ci_lo:.3f}, {ci_hi:.3f}] "
          f"p={perm_p:.4f} (n={len(y)})")
    return {
        "contrast": name,
        "n": len(y),
        "auroc": round(auroc, 4),
        "ci_lower": round(ci_lo, 4),
        "ci_upper": round(ci_hi, 4),
        "perm_p": round(perm_p, 4),
    }


def compute_dge_tstat(X, y):
    """Compute Welch t-statistic per gene (group1 vs group0).
    Returns pd.Series indexed by ENSMUSG.
    """
    group0 = X.loc[y == 0]
    group1 = X.loc[y == 1]
    t_stats = {}
    for gene in X.columns:
        g0 = group0[gene].values
        g1 = group1[gene].values
        if g0.std() == 0 and g1.std() == 0:
            t_stats[gene] = 0.0
            continue
        t, _ = stats.ttest_ind(g1, g0, equal_var=False)
        t_stats[gene] = t if np.isfinite(t) else 0.0
    return pd.Series(t_stats)


def compute_iss_dge(mission):
    """Compute ISS Flight vs Ground DGE t-statistics for a v1 skin mission."""
    X, meta = load_v1_skin_mission(mission)
    if X is None:
        return None

    is_flight = meta["label"] == "Flight"
    is_ground = meta["label"].isin(GROUND_LABELS)
    mask = is_flight | is_ground
    if mask.sum() < 4:
        return None

    X_sub = X.loc[mask]
    y = is_flight.loc[mask].astype(int).values
    if len(np.unique(y)) < 2:
        return None

    return compute_dge_tstat(X_sub, y)


# --- R1: Classification ---

def evaluate_r1():
    """R1: Pairwise classification — Radiation/HLU/Combined vs Control."""
    print("\n=== R1: Radiation Pairwise Classification ===")
    results = {}

    for tissue in TISSUES:
        print(f"\n  [{tissue}]")
        X, meta = load_radiation_data(tissue)
        print(f"    Data: {X.shape[0]} samples × {X.shape[1]} genes")

        tissue_results = []

        # Contrast 1: Radiation effect (NL only, Rad vs Ctrl)
        mask = meta["hlu"] == 0
        if mask.sum() >= 6:
            X_sub = X.loc[mask]
            y_sub = meta.loc[mask, "radiation"].values
            r = classify_pairwise(X_sub, y_sub, "Rad_vs_Ctrl (NL only)")
            if r:
                tissue_results.append(r)

        # Contrast 2: HLU effect (non-irradiated only, HLU vs Ctrl)
        mask = meta["radiation"] == 0
        if mask.sum() >= 6:
            X_sub = X.loc[mask]
            y_sub = meta.loc[mask, "hlu"].values
            r = classify_pairwise(X_sub, y_sub, "HLU_vs_Ctrl (no radiation)")
            if r:
                tissue_results.append(r)

        # Contrast 3: Combined (HLU+Rad) vs Double-Control (NL+NoRad)
        mask_combined = (meta["hlu"] == 1) & (meta["radiation"] == 1)
        mask_ctrl = (meta["hlu"] == 0) & (meta["radiation"] == 0)
        mask = mask_combined | mask_ctrl
        if mask.sum() >= 6:
            X_sub = X.loc[mask]
            y_sub = mask_combined.loc[mask].astype(int).values
            r = classify_pairwise(X_sub, y_sub, "Combined_vs_Ctrl")
            if r:
                tissue_results.append(r)

        results[tissue] = tissue_results

    return results


# --- R2: Gene-Level DGE Concordance ---

def evaluate_r2():
    """R2: Gene-level DGE concordance (ground stressor vs ISS flight).
    Only skin_rad has tissue-matched v1 ISS data.
    """
    print("\n=== R2: Gene-Level DGE Concordance (Ground vs ISS) ===")

    results = {}

    # Compute radiation DGE for skin_rad
    print("\n  [skin_rad] Computing radiation DGE rankings...")
    X_rad, meta_rad = load_radiation_data("skin_rad")

    rad_contrasts = {}

    # Radiation only (within NL)
    mask_nl = meta_rad["hlu"] == 0
    X_nl = X_rad.loc[mask_nl]
    y_rad_only = meta_rad.loc[mask_nl, "radiation"].values
    if len(np.unique(y_rad_only)) == 2:
        rad_contrasts["radiation_only"] = compute_dge_tstat(X_nl, y_rad_only)
        print(f"    radiation_only: {len(rad_contrasts['radiation_only'])} genes")

    # HLU only (within no radiation)
    mask_norad = meta_rad["radiation"] == 0
    X_norad = X_rad.loc[mask_norad]
    y_hlu_only = meta_rad.loc[mask_norad, "hlu"].values
    if len(np.unique(y_hlu_only)) == 2:
        rad_contrasts["hlu_only"] = compute_dge_tstat(X_norad, y_hlu_only)
        print(f"    hlu_only: {len(rad_contrasts['hlu_only'])} genes")

    # Combined vs Control
    mask_comb = (meta_rad["hlu"] == 1) & (meta_rad["radiation"] == 1)
    mask_ctrl = (meta_rad["hlu"] == 0) & (meta_rad["radiation"] == 0)
    mask = mask_comb | mask_ctrl
    X_sub = X_rad.loc[mask]
    y_comb = mask_comb.loc[mask].astype(int).values
    if len(np.unique(y_comb)) == 2:
        rad_contrasts["combined"] = compute_dge_tstat(X_sub, y_comb)
        print(f"    combined: {len(rad_contrasts['combined'])} genes")

    # Compute ISS flight DGE for each v1 skin mission
    print("\n  Computing ISS flight DGE rankings...")
    iss_dge = {}
    for mission in V1_SKIN_MISSIONS:
        tstat = compute_iss_dge(mission)
        if tstat is not None:
            iss_dge[mission] = tstat
            print(f"    ISS {mission}: {len(tstat)} genes")

    # Correlate
    print("\n  Spearman correlations:")
    concordance_results = []
    for contrast_name, rad_tstat in rad_contrasts.items():
        for mission, iss_tstat in iss_dge.items():
            common_genes = sorted(set(rad_tstat.index) & set(iss_tstat.index))
            if len(common_genes) < 100:
                continue

            r_val, p_val = stats.spearmanr(
                rad_tstat.loc[common_genes].values,
                iss_tstat.loc[common_genes].values,
            )

            concordance_results.append({
                "ground_contrast": contrast_name,
                "iss_mission": mission,
                "n_common_genes": len(common_genes),
                "spearman_r": round(float(r_val), 4),
                "spearman_p": float(p_val),
            })
            print(f"    {contrast_name} vs ISS {mission}: "
                  f"r={r_val:.4f} (p={p_val:.2e}, n_genes={len(common_genes)})")

    results["skin_rad"] = concordance_results
    return results


# --- R3: Combined vs Single-Stressor ---

def evaluate_r3(r2_results):
    """R3: Is combined (HLU+Rad) more similar to ISS than single stressors?"""
    print("\n=== R3: Combined vs Single-Stressor Summary ===")

    if "skin_rad" not in r2_results:
        print("  No skin_rad data in R2 results")
        return {}

    concordances = r2_results["skin_rad"]
    if not concordances:
        print("  No concordance results")
        return {}

    # Group by ISS mission
    by_mission = {}
    for c in concordances:
        mission = c["iss_mission"]
        if mission not in by_mission:
            by_mission[mission] = {}
        by_mission[mission][c["ground_contrast"]] = c["spearman_r"]

    summary = []
    for mission, contrasts in by_mission.items():
        valid = {k: v for k, v in contrasts.items() if v is not None}
        if not valid:
            continue
        best = max(valid, key=valid.get)
        row = {
            "iss_mission": mission,
            "r_radiation_only": contrasts.get("radiation_only"),
            "r_hlu_only": contrasts.get("hlu_only"),
            "r_combined": contrasts.get("combined"),
            "best_match": best,
        }
        summary.append(row)
        print(f"  ISS {mission}: "
              f"rad={contrasts.get('radiation_only', '?')}, "
              f"hlu={contrasts.get('hlu_only', '?')}, "
              f"combined={contrasts.get('combined', '?')} "
              f"→ best={best}")

    return {"skin_rad": summary}


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # R1
    r1 = evaluate_r1()
    with open(OUTPUT_DIR / "R1_radiation_classification.json", "w") as f:
        json.dump(r1, f, indent=2)

    # R2
    r2 = evaluate_r2()
    with open(OUTPUT_DIR / "R2_radiation_dge_concordance.json", "w") as f:
        json.dump(r2, f, indent=2)

    # R3
    r3 = evaluate_r3(r2)
    with open(OUTPUT_DIR / "R3_combined_vs_single.json", "w") as f:
        json.dump(r3, f, indent=2)

    # Summary
    print("\n" + "=" * 70)
    print("PHASE 5 EVALUATION SUMMARY")
    print("=" * 70)

    print("\n  R1: Classification")
    for tissue, tissue_results in r1.items():
        for r in tissue_results:
            print(f"    {tissue} / {r['contrast']}: "
                  f"AUROC={r['auroc']:.3f} p={r['perm_p']:.4f}")

    if "skin_rad" in r2 and r2["skin_rad"]:
        print("\n  R2: Gene-Level DGE Concordance (skin_rad vs ISS skin)")
        for c in r2["skin_rad"]:
            print(f"    {c['ground_contrast']} vs ISS {c['iss_mission']}: "
                  f"r={c['spearman_r']:.4f}")


if __name__ == "__main__":
    main()
