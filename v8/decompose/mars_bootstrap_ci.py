"""Mars extrapolation with 95% bootstrap confidence intervals.

Resamples per-gene β ± SE from the factorial analogs, projects to Mars doses,
and computes credible intervals (2.5%, 97.5%) across 1000 resamples.

Outputs: mars_extrapolation_{analog}_with_ci.csv (delta_pred, delta_lo, delta_hi, ratio_lo, ratio_hi)
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DECOMPOSE_DIR = Path(__file__).resolve().parent / "evaluation"

MARS = {
    "HLU": 0.62,
    "IR": 12.9,
    "HLUxIR": 0.62 * 12.9,
    "T": 7.5,
    "HLUxT": 0.62 * 7.5,
    "IRxT": 12.9 * 7.5,
    "HLUxIRxT": 0.62 * 12.9 * 7.5,
}

ISS = {
    "HLU": 1.0,
    "IR": 1.0,
    "HLUxIR": 1.0,
    "T": 1.5,
    "HLUxT": 1.5,
    "IRxT": 1.5,
    "HLUxIRxT": 1.5,
}


def load_betas(analog: str) -> pd.DataFrame | None:
    p = DECOMPOSE_DIR / f"factorial_{analog}_betas.csv"
    if not p.exists():
        return None
    return pd.read_csv(p)


def project_with_bootstrap(betas: pd.DataFrame, dose: dict[str, float], n_resamples: int = 1000) -> pd.DataFrame:
    """Bootstrap Mars predictions: resample β from N(β_estimate, SE), project to dose."""
    beta_cols = [c for c in betas.columns if c.startswith("beta_") and c != "beta_intercept"]
    t_cols = [c for c in betas.columns if c.startswith("t_")]

    # Recover SE from β and t: SE = β / t
    se = pd.DataFrame({
        c.replace("beta_", ""): np.where(
            betas[c.replace("beta_", "t_")] != 0,
            np.abs(betas[c]) / np.abs(betas[c.replace("beta_", "t_")]),
            np.nan
        )
        for c in beta_cols if c.replace("beta_", "t_") in t_cols
    })

    # Bootstrap: resample each β from Normal(β_est, SE)
    predictions = np.zeros((len(betas), n_resamples))
    np.random.seed(42)
    for rep in range(n_resamples):
        delta = np.zeros(len(betas))
        for col in beta_cols:
            stressor = col.replace("beta_", "")
            d = dose.get(stressor, 0.0)
            # Resample β ± SE
            beta_sample = np.random.normal(betas[col].values, se[stressor].fillna(0).values)
            delta = delta + d * beta_sample
        predictions[:, rep] = delta

    # Compute quantiles
    delta_pred = predictions.mean(axis=1)
    delta_lo = np.percentile(predictions, 2.5, axis=1)
    delta_hi = np.percentile(predictions, 97.5, axis=1)

    out = pd.DataFrame({
        "gene": betas["gene"].astype(str).str.upper(),
        "delta_pred": delta_pred,
        "delta_lo": delta_lo,
        "delta_hi": delta_hi,
        "delta_se": predictions.std(axis=1),
    })
    return out


def main() -> None:
    analog_pairs = {
        "spleen": "thymus",
        "skin_analog": "skin",
        "brain": "eye",
    }

    for analog, flight_tissue in analog_pairs.items():
        betas = load_betas(analog)
        if betas is None:
            print(f"[{analog}] no betas cached, skipping")
            continue
        print(f"\n=== Mars bootstrap: {analog} ===")

        mars = project_with_bootstrap(betas, MARS, n_resamples=1000)
        iss = project_with_bootstrap(betas, ISS, n_resamples=1000)

        # Load original extrapolation for comparison
        orig_path = DECOMPOSE_DIR / f"mars_extrapolation_{analog}.csv"
        if orig_path.exists():
            orig = pd.read_csv(orig_path)
            # Merge with bootstrap CIs
            out = pd.DataFrame({
                "gene": mars["gene"],
                "mars_delta": mars["delta_pred"],
                "mars_delta_lo": mars["delta_lo"],
                "mars_delta_hi": mars["delta_hi"],
                "mars_se": mars["delta_se"],
                "iss_delta": iss["delta_pred"],
                "iss_delta_lo": iss["delta_lo"],
                "iss_delta_hi": iss["delta_hi"],
                "iss_se": iss["delta_se"],
                "mars_over_iss": mars["delta_pred"] / iss["delta_pred"].replace(0, np.nan),
                "mars_over_iss_lo": mars["delta_lo"] / iss["delta_hi"].replace(0, np.nan),  # inverted for lower bound
                "mars_over_iss_hi": mars["delta_hi"] / iss["delta_lo"].replace(0, np.nan),
            })
        else:
            out = pd.DataFrame({
                "gene": mars["gene"],
                "mars_delta": mars["delta_pred"],
                "mars_delta_lo": mars["delta_lo"],
                "mars_delta_hi": mars["delta_hi"],
                "mars_se": mars["delta_se"],
                "iss_delta": iss["delta_pred"],
                "iss_delta_lo": iss["delta_lo"],
                "iss_delta_hi": iss["delta_hi"],
                "iss_se": iss["delta_se"],
            })

        out_path = DECOMPOSE_DIR / f"mars_extrapolation_{analog}_with_ci.csv"
        out.to_csv(out_path, index=False)
        print(f"  wrote {out_path}  ({len(out)} genes with 95% CIs)")

        # Print top Mars-predicted genes with CIs
        valid = out.dropna(subset=["mars_delta"])
        valid = valid[np.isfinite(valid["mars_delta"])]
        valid = valid[valid["iss_delta"].abs() > 0.1]
        if len(valid):
            top_up = valid.nlargest(5, "mars_delta")
            top_dn = valid.nsmallest(5, "mars_delta")
            print(f"\n  Top Mars UP-regulated genes:")
            for _, row in top_up.iterrows():
                print(f"    {row['gene']:12s}  {row['mars_delta']:+7.2f} [{row['mars_delta_lo']:+7.2f}, {row['mars_delta_hi']:+7.2f}]")
            print(f"\n  Top Mars DOWN-regulated genes:")
            for _, row in top_dn.iterrows():
                print(f"    {row['gene']:12s}  {row['mars_delta']:+7.2f} [{row['mars_delta_lo']:+7.2f}, {row['mars_delta_hi']:+7.2f}]")


if __name__ == "__main__":
    main()
