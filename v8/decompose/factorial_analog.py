"""v8 Pillar 2 (DECOMPOSE) — factorial HLU × radiation analog decomposition.

Fits per-gene linear model on ground-based 2×2 factorial studies:
  OSD-211 spleen  (HLU × 0.04 Gy 57Co; n=21)
  OSD-237 skin    (HLU × 0.04 Gy 57Co; n=21)

Model: log2(cpm+1) ~ HLU + IR + HLU:IR

Outputs per-gene coefficients (β_HLU, β_IR, β_int), from which we can compute:
  - µg component of gene's flight response
  - radiation component
  - synergy / antagonism

Then for each tissue we compare these analog coefficients against the ISS
flight signature (v8/intervene/signatures/{tissue}_ranked.csv) to estimate
what fraction of the flight signal is explained by µg alone vs GCR vs
interaction — the first-order answer needed for Mars-mission extrapolation.

Outputs
-------
  v8/decompose/evaluation/factorial_{tissue}_betas.csv
  v8/decompose/evaluation/factorial_flight_decomposition.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO_ROOT = Path(__file__).resolve().parents[2]
SIG_DIR = REPO_ROOT / "v8" / "intervene" / "signatures"
OUT_DIR = Path(__file__).resolve().parent / "evaluation"
SYMBOL_MAP = REPO_ROOT / "processed" / "ensembl_symbol_map.csv"


def _load_symbol_map() -> dict[str, str]:
    if not SYMBOL_MAP.exists():
        return {}
    m = pd.read_csv(SYMBOL_MAP)
    return dict(zip(m["ENSEMBL"].astype(str), m["SYMBOL"].astype(str).str.upper()))

DATASETS = {
    "hze_endocrine": {
        "counts": REPO_ROOT / "data/mouse/hze_analog/OSD-719/GLDS-641_rna_seq_Normalized_Counts_GLbulkRNAseq.csv",
        "samples": REPO_ROOT / "data/mouse/hze_analog/OSD-719/GLDS-641_rna_seq_SampleTable_GLbulkRNAseq.csv",
        "cond_map": {
            "Female...sham.irradiated": (0, 0),
            "Male...sham.irradiated":   (1, 0),
            "Female...mixed.radiation.field": (0, 1),
            "Male...mixed.radiation.field":   (1, 1),
        },
        "design": "2x2",  # Sex × HZE (mixed radiation field)
        "col_names_override": ["intercept", "Sex", "HZE", "SexxHZE"],
    },
    "spleen": {
        "counts": REPO_ROOT / "data/mouse/spleen/LDR_HLU/GLDS-211_rna_seq_Normalized_Counts.csv",
        "samples": REPO_ROOT / "data/mouse/spleen/LDR_HLU/GLDS-211_rna_seq_SampleTable.csv",
        "counts_candidates": [
            REPO_ROOT / "data/mouse/spleen/LDR_HLU/GLDS-211_rna_seq_Normalized_Counts.csv",
            REPO_ROOT / "data/mouse/spleen/LDR_HLU/GLDS-211_rna_seq_Normalized_Counts_GLbulkRNAseq.csv",
        ],
        "samples_candidates": [
            REPO_ROOT / "data/mouse/spleen/LDR_HLU/GLDS-211_rna_seq_SampleTable.csv",
            REPO_ROOT / "data/mouse/spleen/LDR_HLU/GLDS-211_rna_seq_SampleTable_GLbulkRNAseq.csv",
        ],
        "cond_map": {
            "HLLC_IRC": (0, 0),
            "HLLC_IR": (0, 1),
            "HLU_IRC": (1, 0),
            "HLU_IR": (1, 1),
            "Normally.Loaded.Control...non.irradiated": (0, 0),
            "Normally.Loaded.Control...cobalt.57.gamma.radiation": (0, 1),
            "Hindlimb.Unloaded...non.irradiated": (1, 0),
            "Hindlimb.Unloaded...cobalt.57.gamma.radiation": (1, 1),
        },
        "design": "2x2",  # HLU × IR
    },
    "skin_analog": {
        "counts": REPO_ROOT / "data/mouse/skin_rad/LDR_HLU/GLDS-237_rna_seq_Normalized_Counts.csv",
        "samples": REPO_ROOT / "data/mouse/skin_rad/LDR_HLU/GLDS-237_rna_seq_SampleTable.csv",
        "cond_map": {
            "Hindlimb.Loaded.Control...Non.irradiated.": (0, 0),
            "Hindlimb.Loaded.Control...Irradiated.with.0.04.Gy.57Co": (0, 1),
            "Hindlimb.Unloaded...Non.irradiated.": (1, 0),
            "Hindlimb.Unloaded...Irradiated.with.0.04.Gy.57Co": (1, 1),
        },
        "design": "2x2",
    },
    "brain": {
        "counts": REPO_ROOT / "data/mouse/retina_brain/LDR_HLU/GLDS-202_rna_seq_Normalized_Counts_GLbulkRNAseq.csv",
        "samples": REPO_ROOT / "data/mouse/retina_brain/LDR_HLU/GLDS-202_rna_seq_SampleTable_GLbulkRNAseq.csv",
        "cond_map": {
            # (HLU, IR, time4mo)
            "non.irradiated...Normally.Loaded.Control...1.month":        (0, 0, 0),
            "non.irradiated...Normally.Loaded.Control...4.month":        (0, 0, 1),
            "non.irradiated...Hindlimb.Unloaded...1.month":              (1, 0, 0),
            "non.irradiated...Hindlimb.Unloaded...4.month":              (1, 0, 1),
            "cobalt.57.gamma.radiation...Normally.Loaded.Control...1.month": (0, 1, 0),
            "cobalt.57.gamma.radiation...Normally.Loaded.Control...4.month": (0, 1, 1),
            "cobalt.57.gamma.radiation...Hindlimb.Unloaded...1.month":       (1, 1, 0),
            "cobalt.57.gamma.radiation...Hindlimb.Unloaded...4.month":       (1, 1, 1),
        },
        "design": "2x2x2",  # HLU × IR × Time(4mo vs 1mo)
    },
}


def _resolve_input_path(spec: dict, role: str) -> Path:
    candidates = spec.get(f"{role}_candidates", [spec[role]])
    for path in candidates:
        p = Path(path)
        if p.exists():
            return p
    return Path(candidates[0])


def fit_factorial_one(tissue_name: str, spec: dict) -> pd.DataFrame:
    counts = pd.read_csv(_resolve_input_path(spec, "counts"), index_col=0)
    samples = pd.read_csv(_resolve_input_path(spec, "samples"), index_col=0)
    # Normalize both index types to str for matching
    counts.columns = counts.columns.astype(str)
    samples.index = samples.index.astype(str)
    samples["condition"] = samples["condition"].astype(str)
    cond_map = spec["cond_map"]
    # Match sample index ↔ count column
    cols = [c for c in counts.columns if c in samples.index]
    if not cols:
        # strip quote/whitespace artifacts
        counts.columns = counts.columns.str.strip().str.strip('"')
        samples.index = samples.index.str.strip().str.strip('"')
        cols = [c for c in counts.columns if c in samples.index]
    if not cols:
        raise RuntimeError(f"No overlapping samples between counts and sampletable for {tissue_name}")
    counts = counts[cols]
    samples = samples.loc[cols]
    # Encode
    sample_cond = samples["condition"].map(cond_map)
    if sample_cond.isna().any():
        raise RuntimeError(f"{tissue_name}: unmapped conditions: {samples['condition'][sample_cond.isna()].unique().tolist()}")
    design = spec.get("design", "2x2")
    if design == "2x2":
        X_hlu = sample_cond.map(lambda t: t[0]).astype(float).values
        X_ir = sample_cond.map(lambda t: t[1]).astype(float).values
        X = np.column_stack([np.ones_like(X_hlu), X_hlu, X_ir, X_hlu * X_ir])
        col_names = spec.get("col_names_override", ["intercept", "HLU", "IR", "HLUxIR"])
    elif design == "2x2x2":
        X_hlu = sample_cond.map(lambda t: t[0]).astype(float).values
        X_ir = sample_cond.map(lambda t: t[1]).astype(float).values
        X_t = sample_cond.map(lambda t: t[2]).astype(float).values
        X = np.column_stack([
            np.ones_like(X_hlu), X_hlu, X_ir, X_t,
            X_hlu * X_ir, X_hlu * X_t, X_ir * X_t, X_hlu * X_ir * X_t,
        ])
        col_names = ["intercept", "HLU", "IR", "T", "HLUxIR", "HLUxT", "IRxT", "HLUxIRxT"]
    else:
        raise RuntimeError(f"Unknown design {design}")

    # Log-transform counts (normalized already)
    Y = np.log2(counts.values.astype(float) + 1.0)  # genes × samples
    n = X.shape[0]
    # OLS: β = (X'X)^-1 X'Y^T; residual df = n-4
    XtX_inv = np.linalg.pinv(X.T @ X)
    beta = Y @ X @ XtX_inv  # genes × 4
    resid = Y - beta @ X.T
    rss = (resid ** 2).sum(axis=1)
    dof = max(n - X.shape[1], 1)
    sigma2 = rss / dof
    se = np.sqrt(np.outer(sigma2, np.diag(XtX_inv)))  # genes × 4
    tstat = beta / np.where(se > 0, se, np.nan)
    # Two-sided p from t; use normal approximation for speed (dof=17 or 17)
    from scipy.stats import t as tdist
    pval = 2 * (1 - tdist.cdf(np.abs(tstat), df=dof))

    ens_to_sym = _load_symbol_map()
    raw_ids = counts.index.astype(str)
    symbols = [ens_to_sym.get(g.split(".")[0], g).upper() for g in raw_ids]
    df_cols: dict = {"gene": symbols, "ensembl_id": raw_ids}
    for i, name in enumerate(col_names):
        df_cols[f"beta_{name}"] = beta[:, i]
        if i > 0:
            df_cols[f"t_{name}"] = tstat[:, i]
            df_cols[f"p_{name}"] = pval[:, i]
    return pd.DataFrame(df_cols)


def load_flight_signature(tissue: str) -> pd.DataFrame | None:
    p = SIG_DIR / f"{tissue}_ranked.csv"
    if not p.exists():
        return None
    df = pd.read_csv(p)
    df["gene"] = df["gene"].astype(str).str.upper()
    return df[["gene", "mean_z"]].rename(columns={"mean_z": "flight_z"})


def decompose(tissue_name: str, flight_tissue: str, betas: pd.DataFrame, out_dir: Path) -> dict:
    """Compare flight Wald-z to all factorial beta columns via Spearman."""
    flight = load_flight_signature(flight_tissue)
    if flight is None:
        return {"flight_tissue": flight_tissue, "note": "no flight signature cached"}
    merged = betas.merge(flight, on="gene", how="inner")
    if len(merged) < 200:
        return {"flight_tissue": flight_tissue, "n_shared": int(len(merged)),
                "note": "insufficient overlap"}

    out = {"flight_tissue": flight_tissue, "analog_tissue": tissue_name,
           "n_shared_genes": int(len(merged))}
    beta_cols = [c for c in betas.columns if c.startswith("beta_") and c != "beta_intercept"]
    for col in beta_cols:
        comp = col.replace("beta_", "")
        rho, p = spearmanr(merged[col], merged["flight_z"])
        out[f"spearman_{comp}_vs_flight"] = {"r": float(rho), "p": float(p)}

    sig = merged[merged["flight_z"].abs() > 1]
    if len(sig) >= 50:
        for col in beta_cols:
            comp = col.replace("beta_", "")
            rho, p = spearmanr(sig[col], sig["flight_z"])
            out[f"spearman_{comp}_vs_flight_responsive"] = {"r": float(rho), "p": float(p), "n": int(len(sig))}

    top100_up = merged.nlargest(100, "flight_z")
    top100_dn = merged.nsmallest(100, "flight_z")
    top = pd.concat([top100_up, top100_dn])
    var_per = {col: float((top[col] ** 2).mean()) for col in beta_cols}
    total = sum(var_per.values()) or 1.0
    out["variance_attribution_top200"] = {col.replace("beta_", "") + "_frac": v / total for col, v in var_per.items()}
    return out


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary: dict = {}
    failures: list[str] = []

    for tissue_name, spec in DATASETS.items():
        print(f"\n=== {tissue_name} ===")
        try:
            betas = fit_factorial_one(tissue_name, spec)
        except Exception as e:
            print(f"  FAILED: {e}")
            summary[tissue_name] = {"error": str(e)}
            failures.append(tissue_name)
            continue
        out_path = OUT_DIR / f"factorial_{tissue_name}_betas.csv"
        betas.to_csv(out_path, index=False)
        print(f"  wrote {out_path}  ({len(betas)} genes)")

        # Reference significant genes for every factor
        sig_counts = {c.replace("p_", ""): int((betas[c] < 0.05).sum())
                      for c in betas.columns if c.startswith("p_")}
        print("  p<0.05:", " ".join(f"{k}={v}" for k, v in sig_counts.items()))

        # Compare to every flight tissue signature we have
        flight_tissues_all = ["thymus", "gastrocnemius", "skin", "eye", "liver", "kidney"]
        per_flight = {}
        for ft in flight_tissues_all:
            d = decompose(tissue_name, ft, betas, OUT_DIR)
            per_flight[ft] = d
        summary[tissue_name] = {"n_sig_p05": sig_counts,
                                "per_flight_tissue": per_flight}
        # Use same-lineage default for printed line
        flight_pairs = {"spleen": "thymus", "skin_analog": "skin", "brain": "eye"}
        flight_tissue = flight_pairs.get(tissue_name, tissue_name)
        decomp = per_flight.get(flight_tissue, {})

        def _fmt(k):
            v = decomp.get(k)
            if isinstance(v, dict) and isinstance(v.get("r"), (int, float)):
                return f"{v['r']:+.3f}"
            return "NA"
        print(f"  decomposition vs flight={flight_tissue}: "
              f"r_HLU={_fmt('spearman_HLU_vs_flight')}  "
              f"r_IR={_fmt('spearman_IR_vs_flight')}  "
              f"r_HLUxIR={_fmt('spearman_HLUxIR_vs_flight')}")
        print(f"    (flight-responsive |z|>1: "
              f"r_HLU={_fmt('spearman_HLU_vs_flight_responsive')}  "
              f"r_IR={_fmt('spearman_IR_vs_flight_responsive')}  "
              f"r_HLUxIR={_fmt('spearman_HLUxIR_vs_flight_responsive')})")
        va = decomp.get("variance_attribution_top200", {})
        if va:
            print(f"  top-200 var attribution: HLU={va['HLU_frac']:.2%} IR={va['IR_frac']:.2%} int={va['HLUxIR_frac']:.2%}")

    (OUT_DIR / "factorial_flight_decomposition.json").write_text(json.dumps(summary, indent=2))
    print(f"\nWrote {OUT_DIR}/factorial_flight_decomposition.json")
    if failures:
        raise SystemExit(f"DECOMPOSE factorial failed for: {', '.join(failures)}")


if __name__ == "__main__":
    main()
