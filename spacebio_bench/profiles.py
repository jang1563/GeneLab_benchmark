"""Metric profiles for SpaceBio-Bench v9 tasks."""

from __future__ import annotations

from copy import deepcopy
from typing import Any


METRIC_PROFILES: dict[str, dict[str, Any]] = {
    "genelab_minimal": {
        "description": "Primary bulk mission-held-out evaluation profile.",
        "metrics": [
            "macro_f1",
            "balanced_accuracy",
            "auroc",
            "calibration_error",
            "mission_discrimination",
        ],
    },
    "genelab_full": {
        "description": "Extended bulk profile for tissue, mission, and domain-shift analysis.",
        "metrics": [
            "macro_f1",
            "balanced_accuracy",
            "auroc",
            "calibration_error",
            "mission_discrimination",
            "per_tissue_f1",
            "cross_mission_transfer_entropy",
            "tissue_heldout_degradation",
            "species_domain_shift_delta",
            "bootstrap_ci",
        ],
    },
    "genelab_sc": {
        "description": "Single-cell profile for AnnData tasks where DE recovery is meaningful.",
        "metrics": [
            "de_overlap_at_n",
            "de_direction_match",
            "mission_discrimination",
            "state_label_auroc",
            "state_label_auprc",
            "expression_mae_when_applicable",
        ],
    },
    "genelab_organoid_pilot": {
        "description": "Draft profile for small public human organoid spaceflight pilots.",
        "metrics": [
            "balanced_accuracy",
            "auroc",
            "calibration_error",
            "de_direction_match",
            "signature_rank_correlation",
            "mission_discrimination",
        ],
    },
    "genelab_multispecies_pilot": {
        "description": "Draft profile for species-native non-mouse spaceflight pilots.",
        "metrics": [
            "balanced_accuracy",
            "auroc",
            "calibration_error",
            "condition_stratum_holdout_delta",
        ],
    },
    "genelab_multispecies_interaction_pilot": {
        "description": "Draft profile for non-mouse interaction-design spaceflight pilots.",
        "metrics": [
            "balanced_accuracy",
            "auroc",
            "calibration_error",
            "genotype_holdout_delta",
            "light_treatment_holdout_delta",
            "condition_stratum_holdout_delta",
        ],
    },
    "stressor_regime": {
        "description": "Radiation-quality and nonlinear-stressor regime profile.",
        "metrics": [
            "low_high_let_sign_consistency",
            "stressor_attribution_stability",
            "nonlinear_saturation_sensitivity",
            "heldout_analog_correlation",
            "uncertainty_width",
            "calibration_error",
        ],
    },
}


def get_metric_profile(name: str) -> dict[str, Any]:
    """Return a defensive copy of a metric profile."""

    try:
        return deepcopy(METRIC_PROFILES[name])
    except KeyError as exc:
        available = ", ".join(sorted(METRIC_PROFILES))
        raise KeyError(f"unknown metric profile {name!r}; available: {available}") from exc
