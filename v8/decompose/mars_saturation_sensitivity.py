"""Conservative bounded-dose sensitivity check for Mars extrapolation.

The linear Mars projection is intentionally a stress test. This script keeps
that output, then asks which high-amplification flags survive simple bounded
dose transforms after 5x analog scale:

  - cap5: hard cap at 5x
  - sqrt_after5: 5 + sqrt(dose - 5)
  - log_after5: 5 + log1p(dose - 5)

These are not fitted nonlinear dose-response models. They are guardrails for
separating robust regime flags from linear-only amplification artifacts.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from mars_extrapolate import ISS, MARS, load_betas, project

DECOMPOSE_DIR = Path(__file__).resolve().parent / "evaluation"
ANALOG_PAIRS = {
    "spleen": "thymus",
    "skin_analog": "skin",
    "brain": "eye",
}
BOUNDED_SCENARIOS = ["cap5", "sqrt_after5", "log_after5"]
EXTREME_RATIO = 5.0
MIN_ABS_ISS_DELTA = 0.1


def _bounded_value(value: float, mode: str) -> float:
    if value <= EXTREME_RATIO:
        return value
    excess = value - EXTREME_RATIO
    if mode == "cap5":
        return EXTREME_RATIO
    if mode == "sqrt_after5":
        return EXTREME_RATIO + math.sqrt(excess)
    if mode == "log_after5":
        return EXTREME_RATIO + math.log1p(excess)
    raise ValueError(f"unknown bounded mode: {mode}")


def transform_dose(dose: dict[str, float], mode: str) -> dict[str, float]:
    """Apply a bounded high-dose transform without changing sub-5x terms."""
    if mode == "linear":
        return dict(dose)
    return {key: _bounded_value(float(value), mode) for key, value in dose.items()}


def _ratio(delta: pd.Series, iss_delta: pd.Series) -> pd.Series:
    denom = iss_delta.where(iss_delta.abs() >= MIN_ABS_ISS_DELTA)
    return delta / denom.replace(0, np.nan)


def summarize_one(analog: str) -> dict:
    betas = load_betas(analog)
    if betas is None:
        return {"error": "missing factorial betas"}

    iss = project(betas, ISS)
    out = pd.DataFrame({
        "gene": iss["gene"],
        "iss_delta": iss["delta_pred"],
    })

    scenario_doses = {"linear": MARS}
    for mode in BOUNDED_SCENARIOS:
        scenario_doses[mode] = transform_dose(MARS, mode)

    for mode, dose in scenario_doses.items():
        pred = project(betas, dose)
        out[f"{mode}_mars_delta"] = pred["delta_pred"]
        out[f"{mode}_mars_over_iss"] = _ratio(pred["delta_pred"], out["iss_delta"])

    bounded_ratio_cols = [f"{mode}_mars_over_iss" for mode in BOUNDED_SCENARIOS]
    bounded_delta_cols = [f"{mode}_mars_delta" for mode in BOUNDED_SCENARIOS]
    linear_ratio = out["linear_mars_over_iss"]
    linear_delta = out["linear_mars_delta"]

    out["linear_extreme_gt5"] = linear_ratio.abs() >= EXTREME_RATIO
    out["min_bounded_abs_over_iss"] = out[bounded_ratio_cols].abs().min(axis=1)
    out["max_bounded_abs_delta"] = out[bounded_delta_cols].abs().max(axis=1)
    out["max_abs_delta_shrink_frac"] = 1.0 - (
        out["max_bounded_abs_delta"] / linear_delta.abs().replace(0, np.nan)
    )
    sign_stable = np.logical_and.reduce([
        np.sign(out[f"{mode}_mars_delta"]) == np.sign(linear_delta)
        for mode in BOUNDED_SCENARIOS
    ])
    out["bounded_sign_stable"] = sign_stable
    out["robust_regime_flag"] = (
        out["linear_extreme_gt5"]
        & out["bounded_sign_stable"]
        & (out["min_bounded_abs_over_iss"] >= EXTREME_RATIO)
    )
    out["saturation_sensitive_flag"] = (
        out["linear_extreme_gt5"]
        & ((out["min_bounded_abs_over_iss"] < EXTREME_RATIO) | ~out["bounded_sign_stable"])
    )

    out_path = DECOMPOSE_DIR / f"mars_saturation_sensitivity_{analog}.csv"
    out.to_csv(out_path, index=False)

    ratio_evaluable = out["linear_mars_over_iss"].notna()
    robust = out[out["robust_regime_flag"]].copy()
    sensitive = out[out["saturation_sensitive_flag"]].copy()
    robust_top = robust.sort_values("min_bounded_abs_over_iss", ascending=False).head(10)
    sensitive["linear_minus_bounded_abs_over_iss"] = (
        sensitive["linear_mars_over_iss"].abs() - sensitive["min_bounded_abs_over_iss"]
    )
    sensitive_top = sensitive.sort_values("linear_minus_bounded_abs_over_iss", ascending=False).head(10)

    top_cols = [
        "gene",
        "linear_mars_over_iss",
        "cap5_mars_over_iss",
        "sqrt_after5_mars_over_iss",
        "log_after5_mars_over_iss",
        "min_bounded_abs_over_iss",
        "bounded_sign_stable",
    ]
    return {
        "flight_tissue": ANALOG_PAIRS[analog],
        "n_genes": int(len(out)),
        "n_ratio_evaluable": int(ratio_evaluable.sum()),
        "n_linear_extreme_gt5": int(out["linear_extreme_gt5"].sum()),
        "n_robust_regime_flags": int(out["robust_regime_flag"].sum()),
        "n_saturation_sensitive_flags": int(out["saturation_sensitive_flag"].sum()),
        "median_shrink_frac_among_linear_extreme": (
            float(out.loc[out["linear_extreme_gt5"], "max_abs_delta_shrink_frac"].median())
            if out["linear_extreme_gt5"].any() else None
        ),
        "output": str(out_path.relative_to(DECOMPOSE_DIR.parent.parent)),
        "top_robust_flags": robust_top[top_cols].to_dict("records"),
        "top_saturation_sensitive_flags": sensitive_top[top_cols].to_dict("records"),
    }


def main() -> None:
    DECOMPOSE_DIR.mkdir(parents=True, exist_ok=True)
    summary = {
        "schema_version": "0.1.0",
        "interpretation": (
            "Bounded-dose scenarios are conservative sensitivity checks. They "
            "do not fit nonlinear biology and should be interpreted as regime "
            "flag robustness screens."
        ),
        "linear_dose": MARS,
        "iss_reference": ISS,
        "bounded_doses": {mode: transform_dose(MARS, mode) for mode in BOUNDED_SCENARIOS},
        "thresholds": {
            "extreme_mars_over_iss": EXTREME_RATIO,
            "min_abs_iss_delta": MIN_ABS_ISS_DELTA,
        },
        "analogs": {},
    }

    for analog in ANALOG_PAIRS:
        result = summarize_one(analog)
        summary["analogs"][analog] = result
        print(
            f"{analog}: linear_extreme={result.get('n_linear_extreme_gt5')} "
            f"robust={result.get('n_robust_regime_flags')} "
            f"saturation_sensitive={result.get('n_saturation_sensitive_flags')}"
        )

    out_path = DECOMPOSE_DIR / "mars_saturation_summary.json"
    out_path.write_text(json.dumps(summary, indent=2))
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
