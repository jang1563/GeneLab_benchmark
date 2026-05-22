"""SpaceBio-Bench v9 planning package.

This package starts as a thin benchmark API for v9 design work. It intentionally
does not promote v9 results; it provides contracts and metrics that can be used
to build reproducible task definitions.
"""

from .manifests import load_task_manifest, validate_task_manifest
from .profiles import METRIC_PROFILES, get_metric_profile
from .registry import TaskRegistry, summarize_task, write_task_index
from .source_audit import (
    audit_source_inventory,
    audit_source_row,
    parse_checksum_manifest,
    write_source_checksum_audit,
)
from .sample_table_audit import (
    audit_sample_table_inventory,
    audit_sample_table_row,
    build_sample_factor_rows,
    parse_condition_factors,
    parse_sample_table,
    write_sample_factor_table,
    write_sample_table_audit,
)
from .sources import build_source_inventory, write_source_inventory
from .evaluate import evaluate_submission
from .reports import (
    read_baseline_summary,
    write_baseline_summary,
    write_evaluation_report,
)
from .submissions import validate_submission
from .data import load_all_bulk_tasks, load_bulk_task, load_human_organoid_task
from .baselines import (
    NearestCentroidConfig,
    OrganoidBaselineConfig,
    SklearnBaselineConfig,
    run_human_organoid_baseline,
    run_human_organoid_donor_diagnostics,
    run_human_organoid_sensitivity_grid,
    run_nearest_centroid_all,
    run_nearest_centroid_task,
    run_organoid_nearest_centroid,
    run_sklearn_baseline_all,
    run_sklearn_baseline_task,
)
from .datapackage import build_v9_public_bulk_datapackage, write_datapackage
from .extension_sources import (
    HUMAN_ORGANOID_DRAFT_SOURCES,
    MULTISPECIES_DRAFT_SOURCES,
    read_extension_source_inventory,
    write_extension_source_inventory,
)
from .extension_tasks import (
    HUMAN_ORGANOID_TASK_ID,
    build_human_organoid_task_manifest,
)
from .expression_audit import (
    audit_expression_matrix_inventory,
    audit_expression_matrix_row,
    inspect_expression_matrix,
    write_expression_matrix_audit,
)
from .organoid_diagnostics import (
    build_organoid_donor_metadata_audit,
    build_organoid_group_diagnostics,
    build_organoid_sample_diagnostics,
    write_organoid_donor_metadata_audit,
    write_organoid_group_diagnostics,
    write_organoid_sample_diagnostics,
)
from .organoid_geo import (
    merge_geo_metadata_with_sample_factors,
    parse_geo_series_matrix,
    read_geo_series_matrix,
    write_geo_sample_metadata,
)
from .organoid_signature_audit import (
    audit_organoid_signature_reference_inventory,
    audit_organoid_signature_reference_row,
    parse_organoid_contrast_table,
    write_organoid_signature_reference_audit,
)

__all__ = [
    "METRIC_PROFILES",
    "NearestCentroidConfig",
    "HUMAN_ORGANOID_DRAFT_SOURCES",
    "HUMAN_ORGANOID_TASK_ID",
    "MULTISPECIES_DRAFT_SOURCES",
    "OrganoidBaselineConfig",
    "SklearnBaselineConfig",
    "TaskRegistry",
    "audit_source_inventory",
    "audit_sample_table_inventory",
    "audit_sample_table_row",
    "audit_expression_matrix_inventory",
    "audit_expression_matrix_row",
    "audit_organoid_signature_reference_inventory",
    "audit_organoid_signature_reference_row",
    "audit_source_row",
    "build_sample_factor_rows",
    "build_source_inventory",
    "build_human_organoid_task_manifest",
    "build_organoid_donor_metadata_audit",
    "build_organoid_group_diagnostics",
    "build_organoid_sample_diagnostics",
    "build_v9_public_bulk_datapackage",
    "evaluate_submission",
    "get_metric_profile",
    "inspect_expression_matrix",
    "load_all_bulk_tasks",
    "load_bulk_task",
    "load_human_organoid_task",
    "load_task_manifest",
    "merge_geo_metadata_with_sample_factors",
    "parse_checksum_manifest",
    "parse_condition_factors",
    "parse_geo_series_matrix",
    "parse_organoid_contrast_table",
    "parse_sample_table",
    "read_geo_series_matrix",
    "read_extension_source_inventory",
    "read_baseline_summary",
    "run_human_organoid_baseline",
    "run_human_organoid_donor_diagnostics",
    "run_human_organoid_sensitivity_grid",
    "run_nearest_centroid_all",
    "run_nearest_centroid_task",
    "run_organoid_nearest_centroid",
    "run_sklearn_baseline_all",
    "run_sklearn_baseline_task",
    "summarize_task",
    "validate_submission",
    "validate_task_manifest",
    "write_baseline_summary",
    "write_datapackage",
    "write_evaluation_report",
    "write_expression_matrix_audit",
    "write_extension_source_inventory",
    "write_geo_sample_metadata",
    "write_organoid_donor_metadata_audit",
    "write_organoid_group_diagnostics",
    "write_organoid_sample_diagnostics",
    "write_organoid_signature_reference_audit",
    "write_source_checksum_audit",
    "write_sample_table_audit",
    "write_sample_factor_table",
    "write_source_inventory",
    "write_task_index",
]
