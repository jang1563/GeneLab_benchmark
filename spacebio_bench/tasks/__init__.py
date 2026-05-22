"""Task-family names for the v9 SpaceBio-Bench skeleton."""

from .legacy import legacy_task_info_to_manifest, load_legacy_task_info, write_manifest

TASK_FAMILIES = (
    "bulk_lomo",
    "bridge_cross_species",
    "sc_spaceflight",
    "stressor_radiation_quality",
    "intervention_hypothesis",
    "human_gated_protocol",
)

__all__ = [
    "TASK_FAMILIES",
    "legacy_task_info_to_manifest",
    "load_legacy_task_info",
    "write_manifest",
]
