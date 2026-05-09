# DECOMPOSE Result Map

Pillar question: which stressor components explain spaceflight signatures, and
where do Mars-like regimes break linear extrapolation?

## Claims And Artifacts

| Claim | Primary output | Script | Status |
|---|---|---|---|
| Ground analog factorial models estimate HLU, radiation, time, and interaction coefficients per gene. | `factorial_*_betas.csv`, `factorial_flight_decomposition.json` | `v8/decompose/factorial_analog.py` | hpc_validated; raw-cache rerun still required before beta freeze |
| Mars projections are regime-change flags, not point predictions. | `mars_extrapolation_*.csv`, `mars_summary.json` | `v8/decompose/mars_extrapolate.py` | hpc_validated; exploratory interpretation |
| Bootstrap uncertainty must accompany any promoted Mars projection. | `mars_extrapolation_*_with_ci.csv` | `v8/decompose/mars_bootstrap_ci.py` | hpc_validated; exploratory interpretation |
| Radiation quality comparisons require low-LET and high-LET claims to be separated. | `factorial_hze_endocrine_betas.csv`, `factorial_flight_decomposition.json` | `v8/decompose/factorial_analog.py` | hpc_validated; mechanism claim remains scoped |

## Promotion Requirements

- Manifest must record OSDR/GLDS accessions, sample-table filters, encoded factor levels, and model formula.
- Mars figures must state dose/time extrapolation assumptions and uncertainty.
- Clinical or operational language is blocked until non-linear calibration and validation exist.
