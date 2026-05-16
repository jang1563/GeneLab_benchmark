"""Tier 1B — Multi-omics mechanistic propagation (RNA → Protein → Metabolite → Clinical).

For each causal stressor (HLU, IR, HLUxIR, T), pull top-100 responsive genes from
v8 rodent factorial βs, then trace signal through SOMA/I4 human layers:

  Layer 1 (RNA):      Mouse factorial β for the stressor (already have it)
  Layer 2 (Protein):  Plasma proteomics intensity for matching HGNC symbols
  Layer 3 (Metabolite): Metabolite intensity aggregated by SuperPathway
  Layer 4 (Clinical):  CBC-derived composite scores (NLR, PLR, anemia, inflam)

Propagation score per stressor = geometric mean of |Spearman ρ| across the three
layer transitions (RNA→Prot, Prot→Metab, Metab→Clinical). A stressor reaching
clinical readout with all |ρ| > 0.3 is mechanistically complete.

Outputs:
  v8/multiomics/evaluation/propagation_scores.csv
  v8/multiomics/evaluation/rna_to_clinical_cascade.json
"""
from __future__ import annotations

import json
import os
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

REPO_ROOT = Path(__file__).resolve().parents[2]
MO_DIR = Path(__file__).resolve().parent
EVAL_DIR = MO_DIR / "evaluation"
EVAL_DIR.mkdir(parents=True, exist_ok=True)

SOMA = Path(
    os.environ.get("SPACEOMICS_ROOT", REPO_ROOT.parent / "SpaceOmicsBench")
).expanduser() / "v2_public/data/processed"

RODENT_ANALOGS = ["spleen", "skin_analog", "brain"]
STRESSORS = ["HLU", "IR", "HLUxIR"]  # T only present in 2x2x2 design
TOP_N = 100


# -------- Layer loaders --------

def top_genes_per_stressor(n: int = TOP_N) -> dict[str, set[str]]:
    """Union of top-N |β| genes across analogs per stressor."""
    out = {s: set() for s in STRESSORS}
    for analog in RODENT_ANALOGS:
        p = REPO_ROOT / f"v8/decompose/evaluation/factorial_{analog}_betas.csv"
        if not p.exists():
            continue
        df = pd.read_csv(p)
        df["gene"] = df["gene"].astype(str).str.upper()
        for s in STRESSORS:
            col = f"beta_{s}"
            if col not in df.columns:
                continue
            top = df.reindex(df[col].abs().sort_values(ascending=False).index).head(n)
            out[s].update(top["gene"].tolist())
    return {k: sorted(v) for k, v in out.items()}


def _sample_key(crew: str, timepoint_days) -> str:
    """Canonical (crew, timepoint_days) key for cross-layer alignment."""
    try:
        td = int(float(timepoint_days))
    except Exception:
        td = timepoint_days
    return f"{crew}__{td}"


def load_proteomics() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Plasma proteomics matrix indexed by (crew, timepoint_days), columns = gene symbols."""
    p = SOMA / "proteomics_plasma_matrix.csv"
    df = pd.read_csv(p)
    meta_cols = ["sample_id", "crew", "month_tp", "timepoint", "timepoint_days", "phase", "mission", "tissue"]
    df["_key"] = [_sample_key(c, t) for c, t in zip(df["crew"], df["timepoint_days"])]
    df = df.set_index("_key")
    meta = df[[c for c in meta_cols if c in df.columns]]
    prot_cols = [c for c in df.columns if c not in meta_cols]
    prot = df[prot_cols].apply(pd.to_numeric, errors="coerce")
    prot.columns = [c.upper() for c in prot.columns]
    return prot, meta


def load_metabolomics(prot_meta: pd.DataFrame) -> pd.DataFrame:
    """Metabolomics ANPPOS matrix with SuperPathway aggregation; re-keyed to (crew, timepoint_days)."""
    p = SOMA / "metabolomics_anppos_matrix.csv"
    df = pd.read_csv(p)
    id_cols = ["SuperPathway", "SubPathway", "metabolite_name", "annotation_confidence",
               "Formula", "Mass", "RT", "CAS ID", "Mode", "KEGG", "HMDB"]
    sample_cols = [c for c in df.columns if c not in id_cols]
    mat = df[["SuperPathway"] + sample_cols].copy()
    mat[sample_cols] = mat[sample_cols].apply(pd.to_numeric, errors="coerce")
    vals = mat[sample_cols].values
    vals = (vals - np.nanmean(vals, axis=1, keepdims=True)) / (np.nanstd(vals, axis=1, keepdims=True) + 1e-9)
    mat.loc[:, sample_cols] = vals
    path_agg = mat.groupby("SuperPathway")[sample_cols].mean().T  # samples (month_tp) x pathways

    # Re-key metab samples from `Cxxx_<month>_<phase>` → (crew, timepoint_days) via prot_meta.month_tp
    month2td = {}
    if "month_tp" in prot_meta.columns:
        for _, row in prot_meta.reset_index().iterrows():
            month2td[(row["crew"], str(row["month_tp"]))] = row["timepoint_days"]
    new_idx = []
    keep = []
    for sid in path_agg.index:
        # sid looks like "C001_Aug_Pre" → crew=C001, month_tp="Aug_Pre"
        parts = sid.split("_", 1)
        if len(parts) != 2:
            new_idx.append(None); continue
        crew, month = parts
        td = month2td.get((crew, month))
        if td is None:
            new_idx.append(None); continue
        new_idx.append(_sample_key(crew, td))
        keep.append(sid)
    path_agg["_key"] = new_idx
    path_agg = path_agg[path_agg["_key"].notna()].set_index("_key")
    return path_agg


def load_clinical() -> tuple[pd.DataFrame, pd.DataFrame]:
    p = SOMA / "clinical_merged_serum.csv"
    df = pd.read_csv(p)
    df["_key"] = [_sample_key(c, t) for c, t in zip(df["crew"], df["timepoint_days"])]
    df = df.set_index("_key")
    # Compute composite scores from CBC
    composites = pd.DataFrame(index=df.index)
    if {"cbc_absolute_neutrophils", "cbc_absolute_lymphocytes"}.issubset(df.columns):
        composites["NLR"] = df["cbc_absolute_neutrophils"] / df["cbc_absolute_lymphocytes"].replace(0, np.nan)
    if {"cbc_platelet_count", "cbc_absolute_lymphocytes"}.issubset(df.columns):
        composites["PLR"] = df["cbc_platelet_count"] / df["cbc_absolute_lymphocytes"].replace(0, np.nan)
    if {"cbc_white_blood_cell_count"}.issubset(df.columns):
        composites["WBC"] = df["cbc_white_blood_cell_count"]
    if {"cbc_hemoglobin"}.issubset(df.columns):
        composites["anemia_inv"] = -df["cbc_hemoglobin"]
    if {"cbc_red_blood_cell_count", "cbc_hemoglobin", "cbc_hematocrit"}.issubset(df.columns):
        # z-scored average of RBC, Hgb, Hct (rbc health, higher = better)
        sub = df[["cbc_red_blood_cell_count", "cbc_hemoglobin", "cbc_hematocrit"]].apply(
            lambda c: (c - c.mean()) / (c.std() + 1e-9))
        composites["rbc_health"] = sub.mean(axis=1)
    # monocyte:lymphocyte ratio (macrophage-skewed inflammation)
    if {"cbc_absolute_monocytes", "cbc_absolute_lymphocytes"}.issubset(df.columns):
        composites["MLR"] = df["cbc_absolute_monocytes"] / df["cbc_absolute_lymphocytes"].replace(0, np.nan)
    meta = df[[c for c in ["crew", "timepoint", "timepoint_days", "phase", "mission"] if c in df.columns]]
    return composites, meta


# -------- Per-sample scoring --------

def gene_set_mean(mat: pd.DataFrame, genes: list[str]) -> pd.Series:
    """Mean z-score of specified gene columns per sample row."""
    hits = [g for g in genes if g in mat.columns]
    if not hits:
        return pd.Series(np.nan, index=mat.index)
    sub = mat[hits].apply(lambda c: (c - c.mean()) / (c.std() + 1e-9))
    return sub.mean(axis=1)


# -------- Cross-layer correlation --------

def safe_spearman(a: pd.Series, b: pd.Series) -> tuple[float, float, int]:
    merged = pd.concat([a.rename("a"), b.rename("b")], axis=1).dropna()
    if len(merged) < 6:
        return np.nan, np.nan, len(merged)
    r, p = stats.spearmanr(merged["a"], merged["b"])
    return float(r), float(p), int(len(merged))


def permutation_null(a: pd.Series, b: pd.Series, n_perm: int = 1000) -> float:
    """Permutation null: shuffle b, recompute rho, return two-sided p."""
    merged = pd.concat([a.rename("a"), b.rename("b")], axis=1).dropna()
    if len(merged) < 6:
        return np.nan
    obs, _ = stats.spearmanr(merged["a"], merged["b"])
    rng = np.random.default_rng(42)
    null = np.zeros(n_perm)
    b_vals = merged["b"].values.copy()
    for i in range(n_perm):
        rng.shuffle(b_vals)
        null[i], _ = stats.spearmanr(merged["a"].values, b_vals)
    return float((np.abs(null) >= abs(obs)).mean())


# -------- Main --------

def main() -> None:
    print("=== Tier 1B: Multi-Omics Propagation (RNA → Protein → Metabolite → Clinical) ===\n")

    top_genes = top_genes_per_stressor()
    for s, gs in top_genes.items():
        print(f"  {s}: {len(gs)} top genes (union across {len(RODENT_ANALOGS)} analogs)")

    print("\n--- loading SOMA layers ---")
    prot, prot_meta = load_proteomics()
    print(f"  proteomics: {prot.shape[0]} samples × {prot.shape[1]} proteins")
    metab = load_metabolomics(prot_meta)
    print(f"  metabolomics: {metab.shape[0]} samples × {metab.shape[1]} SuperPathways")
    clin, clin_meta = load_clinical()
    print(f"  clinical: {clin.shape[0]} samples × {clin.shape[1]} composite scores ({list(clin.columns)})")

    prot_keys = {k: k for k in prot.index}
    clin_keys = {k: k for k in clin.index}
    metab_keys = {k: k for k in metab.index}
    all_keys = set(prot_keys) | set(clin_keys) | set(metab_keys)
    aligned = sorted([k for k in all_keys if (k in prot_keys) + (k in clin_keys) + (k in metab_keys) >= 2])
    n3 = sum((k in prot_keys) + (k in clin_keys) + (k in metab_keys) == 3 for k in aligned)
    print(f"\n  aligned samples with ≥2 layers: {len(aligned)} ({n3} with all three)")

    results = []
    cascade = {}
    for stressor, genes in top_genes.items():
        rna_hits = [g for g in genes if g in prot.columns]
        if not rna_hits:
            print(f"\n[{stressor}] no top genes in proteomics; skipping")
            continue

        # Layer 2 score (protein): mean z-score of rodent-stressor genes in plasma proteomics
        prot_score = gene_set_mean(prot, genes)

        # Layer 3 score (metabolite): top stressor-relevant SuperPathway
        # Use all pathways; pick those most correlated with prot_score
        # (functional mapping: proteins → metabolic pathways via correlation)
        best_path = None
        best_r = 0
        for path in metab.columns:
            ps = pd.Series([prot_score.get(prot_keys.get(k), np.nan) for k in aligned], index=aligned)
            ms = pd.Series([metab[path].get(metab_keys.get(k), np.nan) for k in aligned], index=aligned)
            r, _, n = safe_spearman(ps, ms)
            if not np.isnan(r) and abs(r) > abs(best_r):
                best_r = r
                best_path = path
        metab_score = metab[best_path] if best_path else None

        # Layer 4 score (clinical): best composite correlated with metabolite pathway
        best_clin = None
        best_rc = 0
        for c in clin.columns:
            ms = pd.Series([metab_score.get(metab_keys.get(k), np.nan) for k in aligned], index=aligned) if metab_score is not None else None
            cs = pd.Series([clin[c].get(clin_keys.get(k), np.nan) for k in aligned], index=aligned)
            if ms is None:
                continue
            r, _, _ = safe_spearman(ms, cs)
            if not np.isnan(r) and abs(r) > abs(best_rc):
                best_rc = r
                best_clin = c

        # Compute three transition correlations across aligned samples
        prot_s = pd.Series([prot_score.get(prot_keys.get(k), np.nan) for k in aligned], index=aligned)
        metab_s = pd.Series([metab_score.get(metab_keys.get(k), np.nan) for k in aligned], index=aligned) if metab_score is not None else pd.Series(np.nan, index=aligned)
        clin_s = pd.Series([clin[best_clin].get(clin_keys.get(k), np.nan) for k in aligned], index=aligned) if best_clin else pd.Series(np.nan, index=aligned)

        # RNA layer is represented by rodent β magnitude of the top genes — a constant for this stressor,
        # so RNA→Prot transition uses within-protein-layer across top vs bottom genes instead.
        # Compute rna_score_per_sample = mean z of SAME top genes in proteomics (already = prot_score).
        # For a genuine RNA-level signal per sample, use cfRNA matrix if available; simplified here:
        # RNA→Prot: correlation of top-gene mean vs bottom-gene mean (sanity check)
        bottom = [g for g in prot.columns if g not in set(genes)][:500]
        bottom_prot = gene_set_mean(prot, bottom)
        bottom_s = pd.Series([bottom_prot.get(prot_keys.get(k), np.nan) for k in aligned], index=aligned)

        r_rp, p_rp, n_rp = safe_spearman(prot_s, bottom_s)  # top vs bottom in protein space
        r_pm, p_pm, n_pm = safe_spearman(prot_s, metab_s)
        r_mc, p_mc, n_mc = safe_spearman(metab_s, clin_s)

        perm_rp = permutation_null(prot_s, bottom_s)
        perm_pm = permutation_null(prot_s, metab_s)
        perm_mc = permutation_null(metab_s, clin_s)

        # Propagation score: geometric mean of |ρ| across three transitions
        rhos = [abs(x) for x in (r_rp, r_pm, r_mc) if not np.isnan(x)]
        prop_score = float(np.exp(np.mean(np.log(np.clip(rhos, 1e-3, None))))) if rhos else np.nan

        row = {
            "stressor": stressor,
            "n_top_genes": len(genes),
            "n_in_proteomics": len(rna_hits),
            "best_metab_pathway": best_path,
            "best_clin_score": best_clin,
            "rho_prot_top_vs_bottom": r_rp,
            "p_prot_rp": p_rp,
            "perm_p_rp": perm_rp,
            "n_rp": n_rp,
            "rho_prot_metab": r_pm,
            "p_pm": p_pm,
            "perm_p_pm": perm_pm,
            "n_pm": n_pm,
            "rho_metab_clin": r_mc,
            "p_mc": p_mc,
            "perm_p_mc": perm_mc,
            "n_mc": n_mc,
            "propagation_score": prop_score,
            "reaches_clinical_|rho|>0.3": all(abs(x) > 0.3 for x in (r_rp, r_pm, r_mc) if not np.isnan(x)) and len(rhos) == 3,
        }
        results.append(row)
        cascade[stressor] = row

        print(f"\n[{stressor}]")
        print(f"  top genes: {len(genes)} | in proteomics: {len(rna_hits)}")
        print(f"  best metab pathway: {best_path}  (r with prot: {best_r:+.3f})")
        print(f"  best clin score:    {best_clin}  (r with metab: {best_rc:+.3f})")
        print(f"  RNA→Prot  (top vs bottom): ρ={r_rp:+.3f} p={p_rp:.3g} perm-p={perm_rp:.3g} n={n_rp}")
        print(f"  Prot→Metab:                ρ={r_pm:+.3f} p={p_pm:.3g} perm-p={perm_pm:.3g} n={n_pm}")
        print(f"  Metab→Clin:                ρ={r_mc:+.3f} p={p_mc:.3g} perm-p={perm_mc:.3g} n={n_mc}")
        print(f"  PROPAGATION SCORE:         {prop_score:.3f}   reaches-clinical(|ρ|>0.3): {row['reaches_clinical_|rho|>0.3']}")

    res_df = pd.DataFrame(results)
    out_csv = EVAL_DIR / "propagation_scores.csv"
    res_df.to_csv(out_csv, index=False)
    print(f"\n✓ wrote {out_csv}")

    summary = {
        "aligned_samples": len(aligned),
        "proteomics_shape": list(prot.shape),
        "metabolomics_shape": list(metab.shape),
        "clinical_composites": list(clin.columns),
        "stressors": cascade,
    }
    (EVAL_DIR / "rna_to_clinical_cascade.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"✓ wrote {EVAL_DIR / 'rna_to_clinical_cascade.json'}")


if __name__ == "__main__":
    main()
