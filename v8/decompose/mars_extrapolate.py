"""v8 Pillar 2 — Mars 900-day mission extrapolation from factorial betas.

Given per-gene factorial coefficients β_HLU, β_IR, β_HLUxIR, β_T (and 3-way
interactions where available), project the predicted transcriptomic response
to a Mars transit + surface mission:

  Mars dose setpoint
  ------------------
    µg        = 1.0 (full microgravity, transit) scaled by Mars gravity 0.38 on surface → 0.62 avg
    HZE       = 350 mGy / 900 d  vs ISS LEO ~70 mGy / 180 d → scale factor ≈ 12.9×
    time      = 900 d           vs analog 4 mo = 120 d      → scale factor ≈ 7.5×
    (the time axis is only defined for OSD-202 brain; elsewhere extrapolated
    from ISS mission-duration slope)

Predicted Mars response per gene:
  ΔExpr_Mars = β_HLU·µg_Mars + β_IR·HZE_Mars + β_HLUxIR·(µg·HZE) + β_T·t_Mars
             + β_HLUxT·(µg·t) + β_IRxT·(HZE·t) + β_HLUxIRxT·(µg·HZE·t)

with uncertainty propagated from β standard errors (delta-method 1-sigma CI).

Outputs
-------
  v8/decompose/evaluation/mars_extrapolation_{tissue}.csv    # per-gene predictions
  v8/decompose/evaluation/mars_summary.json                  # top-movers per tissue
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
DECOMPOSE_DIR = Path(__file__).resolve().parent / "evaluation"
SIG_DIR = REPO_ROOT / "v8" / "intervene" / "signatures"

# Mars mission dose vector (normalized to analog reference = 1.0)
MARS = {
    "HLU": 0.62,      # avg over transit (1.0) + Mars surface (0.38) = 0.62
    "IR": 12.9,       # 350 mGy/900d vs 0.04 Gy analog scaled for LEO
    "HLUxIR": 0.62 * 12.9,
    "T": 7.5,         # 900d vs 120d analog
    "HLUxT": 0.62 * 7.5,
    "IRxT": 12.9 * 7.5,
    "HLUxIRxT": 0.62 * 12.9 * 7.5,
}

ISS = {  # ISS reference for comparison (180 d nominal)
    "HLU": 1.0,
    "IR": 1.0,
    "HLUxIR": 1.0,
    "T": 1.5,   # 180 d / 120 d
    "HLUxT": 1.5,
    "IRxT": 1.5,
    "HLUxIRxT": 1.5,
}


def load_betas(analog: str) -> pd.DataFrame | None:
    p = DECOMPOSE_DIR / f"factorial_{analog}_betas.csv"
    if not p.exists():
        return None
    return pd.read_csv(p)


def project(betas: pd.DataFrame, dose: dict[str, float]) -> pd.DataFrame:
    """Return DataFrame with columns: gene, delta_pred, delta_se."""
    beta_cols = [c for c in betas.columns if c.startswith("beta_") and c != "beta_intercept"]
    t_cols = [c for c in betas.columns if c.startswith("t_")]

    # Recover SE from β and t
    se = pd.DataFrame({c.replace("beta_", ""): np.where(
        betas[c.replace("beta_", "t_")] != 0,
        betas[c] / betas[c.replace("beta_", "t_")],
        np.nan
    ) for c in beta_cols if c.replace("beta_", "t_") in t_cols})

    delta = pd.Series(0.0, index=betas.index)
    var = pd.Series(0.0, index=betas.index)
    for col in beta_cols:
        stressor = col.replace("beta_", "")
        d = dose.get(stressor, 0.0)
        delta = delta + d * betas[col].fillna(0)
        if stressor in se.columns:
            var = var + (d ** 2) * (se[stressor].fillna(0) ** 2)
    out = pd.DataFrame({
        "gene": betas["gene"].astype(str).str.upper(),
        "delta_pred": delta,
        "delta_se": np.sqrt(var),
    })
    return out


def compare_to_flight(mars: pd.DataFrame, tissue: str) -> dict:
    """Compare Mars-predicted Δ to cached flight Wald-z per gene."""
    p = SIG_DIR / f"{tissue}_ranked.csv"
    if not p.exists():
        return {"note": f"no flight signature for {tissue}"}
    sig = pd.read_csv(p)
    sig["gene"] = sig["gene"].astype(str).str.upper()
    m = mars.merge(sig[["gene", "mean_z"]], on="gene", how="inner")
    if len(m) < 200:
        return {"note": "insufficient overlap", "n": len(m)}
    from scipy.stats import spearmanr
    rho, p_val = spearmanr(m["delta_pred"], m["mean_z"])
    # Top Mars-predicted movers that also appear in flight signature
    m["mars_z"] = m["delta_pred"] / (m["delta_se"].replace(0, np.nan))
    m["combined_rank"] = m["delta_pred"].rank() + m["mean_z"].rank()
    top_up = m.nlargest(10, "delta_pred")[["gene", "delta_pred", "delta_se", "mean_z"]].to_dict("records")
    top_dn = m.nsmallest(10, "delta_pred")[["gene", "delta_pred", "delta_se", "mean_z"]].to_dict("records")
    return {
        "n_shared": int(len(m)),
        "spearman_r": round(float(rho), 4),
        "spearman_p": float(p_val),
        "top10_mars_up": top_up,
        "top10_mars_dn": top_dn,
    }


def main() -> None:
    summary: dict = {"mars_doses": MARS, "iss_reference": ISS, "analogs": {}}
    # Analog → flight-tissue pairing
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
        print(f"\n=== Mars extrapolation: {analog} → flight {flight_tissue} ===")

        mars = project(betas, MARS)
        iss = project(betas, ISS)

        out = pd.DataFrame({
            "gene": mars["gene"],
            "mars_delta": mars["delta_pred"],
            "mars_se": mars["delta_se"],
            "iss_delta": iss["delta_pred"],
            "iss_se": iss["delta_se"],
            "mars_over_iss": mars["delta_pred"] / iss["delta_pred"].replace(0, np.nan),
        })
        out_path = DECOMPOSE_DIR / f"mars_extrapolation_{analog}.csv"
        out.to_csv(out_path, index=False)
        print(f"  wrote {out_path}  ({len(out)} genes)")

        # Compare to flight signature
        comp = compare_to_flight(mars, flight_tissue)
        print(f"  vs flight[{flight_tissue}]: r={comp.get('spearman_r','NA')}  n={comp.get('n_shared','NA')}")

        # Top Mars vs ISS amplification
        valid = out.dropna(subset=["mars_over_iss"])
        valid = valid[np.isfinite(valid["mars_over_iss"])]
        valid = valid[valid["iss_delta"].abs() > 0.1]  # require meaningful ISS response
        if len(valid):
            top_amp = valid.nlargest(10, "mars_over_iss")[["gene", "mars_delta", "iss_delta", "mars_over_iss"]].to_dict("records")
            print(f"  Top Mars-amplified genes (|iss|>0.1): " +
                  ", ".join(f"{r['gene']}({r['mars_over_iss']:+.1f}x)" for r in top_amp[:5]))
        else:
            top_amp = []

        summary["analogs"][analog] = {
            "flight_tissue": flight_tissue,
            "n_genes": int(len(out)),
            "comparison_to_flight": comp,
            "top10_mars_amplification": top_amp,
        }

    (DECOMPOSE_DIR / "mars_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"\nWrote {DECOMPOSE_DIR}/mars_summary.json")


if __name__ == "__main__":
    main()
