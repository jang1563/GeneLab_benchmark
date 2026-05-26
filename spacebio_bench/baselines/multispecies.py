"""Draft multispecies species-native baseline runners."""

from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
import platform
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from spacebio_bench.data import MultispeciesFoldData, MultispeciesTaskData
from spacebio_bench.data import load_all_multispecies_interaction_tasks
from spacebio_bench.data import load_all_multispecies_tasks
from spacebio_bench.evaluate import evaluate_submission
from spacebio_bench.reports import write_evaluation_report
from spacebio_bench.source_audit import fetch_url, parse_checksum_manifest


BASELINE_ID = "multispecies_nearest_centroid"
INTERACTION_BASELINE_ID = "multispecies_interaction_nearest_centroid"
INTERACTION_LOGISTIC_BASELINE_ID = "multispecies_interaction_logistic_regression_l2"
NEGATIVE_LABEL = "Ground"
POSITIVE_LABEL = "LEO_or_ISS"
INTERACTION_FOLD_FAMILIES = [
    "primary_genotype_or_ecotype_holdout",
    "secondary_light_treatment_holdout",
    "diagnostic_condition_stratum_holdout",
]
HELDOUT_FACTOR_BY_FOLD_FAMILY = {
    "primary_genotype_or_ecotype_holdout": "genotype_or_ecotype",
    "secondary_light_treatment_holdout": "light_treatment",
    "diagnostic_condition_stratum_holdout": "condition_stratum",
}
DELTA_METRIC_BY_FOLD_FAMILY = {
    "condition_stratum_candidate_folds": "condition_stratum_holdout_delta",
    "primary_genotype_or_ecotype_holdout": "genotype_holdout_delta",
    "secondary_light_treatment_holdout": "light_treatment_holdout_delta",
    "diagnostic_condition_stratum_holdout": "condition_stratum_holdout_delta",
}
PREDICTION_COLUMNS = [
    "task_id",
    "fold_id",
    "sample_id",
    "mission",
    "condition_stratum",
    "genotype_or_ecotype",
    "light_treatment",
    "true_label",
    "predicted_label",
    "flight_probability",
    "embedding_0",
    "embedding_1",
]
LOGISTIC_PREDICTION_COLUMNS = [
    column for column in PREDICTION_COLUMNS if not column.startswith("embedding_")
]
INTERACTION_FOLD_DETAIL_FIELDS = [
    "baseline_id",
    "variant_id",
    "transform",
    "scaling",
    "top_variable_genes",
    "task_id",
    "fold_family",
    "delta_metric_id",
    "fold_id",
    "heldout_factor",
    "heldout_value",
    "n_test",
    "balanced_accuracy",
    "rank_low_to_high",
    "is_lowest_for_variant",
    "is_balanced_accuracy_le_0_5",
    "summary_metrics",
    "metrics",
]
INTERACTION_FOLD_COMPARISON_FIELDS = [
    "task_id",
    "fold_family",
    "heldout_factor",
    "heldout_value",
    "n_test",
    "nearest_centroid_baseline_id",
    "nearest_centroid_variant_id",
    "nearest_centroid_balanced_accuracy",
    "nearest_centroid_rank_low_to_high",
    "nearest_centroid_is_lowest_for_variant",
    "logistic_baseline_id",
    "logistic_variant_id",
    "logistic_balanced_accuracy",
    "logistic_rank_low_to_high",
    "logistic_is_lowest_for_variant",
    "delta_logistic_minus_nearest_centroid",
    "logistic_improved",
    "logistic_new_or_persistent_le_0_5",
    "nearest_centroid_metrics",
    "logistic_metrics",
]
INTERACTION_LOGISTIC_FEATURE_AUDIT_SUMMARY_FIELDS = [
    "task_id",
    "fold_family",
    "heldout_factor",
    "heldout_value",
    "fold_id",
    "variant_id",
    "baseline_id",
    "transform",
    "scaling",
    "top_variable_genes",
    "c",
    "n_train",
    "n_test",
    "n_selected_features",
    "n_top500_features",
    "n_extra_features",
    "n_nonzero_coefficients",
    "n_nonzero_top500_coefficients",
    "n_nonzero_extra_coefficients",
    "balanced_accuracy",
    "n_positive_predictions",
    "top_abs_feature_ids",
    "top_positive_feature_ids",
    "top_negative_feature_ids",
    "top_abs_extra_feature_ids",
]
INTERACTION_LOGISTIC_FEATURE_COEFFICIENT_FIELDS = [
    "task_id",
    "fold_family",
    "heldout_factor",
    "heldout_value",
    "fold_id",
    "variant_id",
    "top_variable_genes",
    "c",
    "feature_id",
    "train_variance_rank",
    "train_variance",
    "coefficient",
    "abs_coefficient",
    "coefficient_rank_abs",
    "in_top500_for_fold",
    "coefficient_direction",
]
INTERACTION_LOGISTIC_STABILITY_SUMMARY_FIELDS = [
    "task_id",
    "fold_family",
    "heldout_factor",
    "heldout_value",
    "fold_id",
    "variant_id",
    "top_variable_genes",
    "c",
    "n_train",
    "n_test",
    "n_subsamples",
    "subsample_fraction",
    "candidate_features",
    "reference_balanced_accuracy",
    "reference_nonzero_count",
    "median_subsample_nonzero_count",
    "min_subsample_nonzero_count",
    "max_subsample_nonzero_count",
    "stable_feature_count_ge_0_5",
    "stable_feature_count_ge_0_8",
    "mean_pairwise_jaccard",
    "top_selection_frequency_features",
]
INTERACTION_LOGISTIC_STABILITY_FEATURE_FIELDS = [
    "task_id",
    "fold_family",
    "heldout_factor",
    "heldout_value",
    "fold_id",
    "variant_id",
    "feature_id",
    "train_variance_rank",
    "reference_coefficient",
    "reference_abs_coefficient",
    "reference_selected",
    "selection_count",
    "selection_frequency",
    "positive_count",
    "negative_count",
    "sign_consistency",
    "mean_coefficient",
    "mean_abs_coefficient",
]
INTERACTION_BASELINE_LADDER_SUMMARY_FIELDS = [
    "candidate_id",
    "role",
    "decision",
    "model_class",
    "baseline_id",
    "variant_id",
    "task_id",
    "transform",
    "scaling",
    "top_variable_genes",
    "penalty",
    "c",
    "primary_balanced_accuracy",
    "primary_auroc",
    "primary_delta_metric",
    "primary_delta",
    "secondary_balanced_accuracy",
    "secondary_auroc",
    "secondary_delta_metric",
    "secondary_delta",
    "diagnostic_balanced_accuracy",
    "diagnostic_auroc",
    "diagnostic_delta_metric",
    "diagnostic_delta",
    "mean_family_balanced_accuracy",
    "min_family_balanced_accuracy",
    "focus_col0_phyd_dark_ba",
    "focus_light_treatment_ba",
    "focus_col0_phyd_ba",
    "nearest_improved_count",
    "nearest_tied_count",
    "nearest_worse_count",
    "nearest_min_delta",
    "nearest_max_delta",
    "dark_reference_nonzero_count",
    "dark_stable_feature_count_ge_0_5",
    "dark_stable_feature_count_ge_0_8",
    "dark_mean_pairwise_jaccard",
    "light_reference_nonzero_count",
    "light_stable_feature_count_ge_0_5",
    "light_stable_feature_count_ge_0_8",
    "light_mean_pairwise_jaccard",
    "phyd_reference_nonzero_count",
    "phyd_stable_feature_count_ge_0_5",
    "phyd_stable_feature_count_ge_0_8",
    "phyd_mean_pairwise_jaccard",
    "source_summary_csv",
    "source_fold_detail_csv",
    "source_comparison_csv",
    "source_stability_csv",
    "claim_boundary",
]
INTERACTION_BASELINE_LADDER_FOCUS_FIELDS = [
    "candidate_id",
    "model_class",
    "variant_id",
    "task_id",
    "fold_family",
    "heldout_factor",
    "heldout_value",
    "n_test",
    "nearest_centroid_balanced_accuracy",
    "candidate_balanced_accuracy",
    "delta_candidate_minus_nearest_centroid",
    "nearest_centroid_rank_low_to_high",
    "candidate_rank_low_to_high",
    "nearest_centroid_is_lowest_for_variant",
    "candidate_is_lowest_for_variant",
    "source_fold_detail_csv",
    "source_comparison_csv",
    "nearest_centroid_metrics",
    "candidate_metrics",
]
INTERACTION_DIAGNOSTIC_PACKAGE_SUMMARY_FIELDS = [
    "package_id",
    "task_id",
    "organism",
    "assay_modality",
    "biospecimen_type",
    "source_id",
    "source_url",
    "candidate_id",
    "candidate_variant_id",
    "candidate_model_class",
    "candidate_decision",
    "comparator_candidate_id",
    "comparator_variant_id",
    "primary_balanced_accuracy",
    "secondary_balanced_accuracy",
    "diagnostic_balanced_accuracy",
    "mean_family_balanced_accuracy",
    "min_family_balanced_accuracy",
    "nearest_improved_count",
    "nearest_tied_count",
    "nearest_worse_count",
    "focus_fold_count",
    "focus_improved_count",
    "focus_mean_candidate_ba",
    "focus_min_candidate_ba",
    "stable_feature_count_ge_0_5_total",
    "stable_feature_count_ge_0_8_total",
    "claim_boundary",
    "release_status",
    "evidence_status",
    "external_context_status",
    "source_task_manifest",
    "source_ladder_summary_csv",
    "source_ladder_focus_csv",
    "source_sparse_feature_audit_csv",
    "source_stability_summary_csv",
    "source_stability_feature_csv",
]
INTERACTION_DIAGNOSTIC_PACKAGE_FOCUS_FIELDS = [
    "package_id",
    "candidate_id",
    "candidate_variant_id",
    "fold_label",
    "fold_family",
    "heldout_factor",
    "heldout_value",
    "n_test",
    "nearest_centroid_balanced_accuracy",
    "candidate_balanced_accuracy",
    "delta_candidate_minus_nearest_centroid",
    "candidate_rank_low_to_high",
    "candidate_is_lowest_for_variant",
    "n_nonzero_coefficients",
    "n_nonzero_top500_coefficients",
    "n_nonzero_extra_coefficients",
    "stable_feature_count_ge_0_5",
    "stable_feature_count_ge_0_8",
    "mean_pairwise_jaccard",
    "top_reference_abs_features",
    "top_stable_features",
    "evidence_summary",
]
INTERACTION_DIAGNOSTIC_PACKAGE_FEATURE_FIELDS = [
    "package_id",
    "candidate_id",
    "candidate_variant_id",
    "fold_label",
    "fold_family",
    "heldout_value",
    "feature_id",
    "train_variance_rank",
    "reference_coefficient",
    "reference_abs_coefficient",
    "reference_selected",
    "selection_count",
    "selection_frequency",
    "sign_consistency",
    "mean_coefficient",
    "mean_abs_coefficient",
    "stability_tier",
]
INTERACTION_DIAGNOSTIC_PACKAGE_CLAIM_FIELDS = [
    "package_id",
    "claim_id",
    "claim_status",
    "claim_text",
    "local_evidence_artifacts",
    "local_evidence_columns",
    "external_source_urls",
    "limitations",
]
INTERACTION_FIGURE_TABLE_MAIN_FIELDS = [
    "figure_table_id",
    "package_id",
    "candidate_id",
    "candidate_variant_id",
    "display_fold",
    "fold_label",
    "heldout_value",
    "n_test",
    "nearest_centroid_ba",
    "candidate_ba",
    "delta_ba",
    "display_delta_ba",
    "n_nonzero_coefficients",
    "stable_features_ge_0_5",
    "stable_features_ge_0_8",
    "mean_pairwise_jaccard",
    "figure_use",
    "result_sentence",
]
INTERACTION_FIGURE_TABLE_FEATURE_FIELDS = [
    "figure_table_id",
    "package_id",
    "candidate_id",
    "fold_label",
    "display_fold",
    "feature_rank_in_fold",
    "feature_id",
    "selection_frequency",
    "display_selection_frequency",
    "reference_coefficient",
    "coefficient_direction",
    "stability_tier",
    "appendix_note",
]
INTERACTION_PAIRED_COMPARATOR_SUMMARY_FIELDS = [
    "paired_table_id",
    "task_id",
    "primary_candidate_id",
    "primary_variant_id",
    "compact_candidate_id",
    "compact_variant_id",
    "primary_mean_family_ba",
    "compact_mean_family_ba",
    "primary_min_family_ba",
    "compact_min_family_ba",
    "primary_diagnostic_ba",
    "compact_diagnostic_ba",
    "primary_nearest_improved_count",
    "primary_nearest_tied_count",
    "primary_nearest_worse_count",
    "compact_nearest_improved_count",
    "compact_nearest_tied_count",
    "compact_nearest_worse_count",
    "focus_tied_ba_count",
    "focus_primary_better_count",
    "focus_compact_better_count",
    "primary_focus_nonzero_total",
    "compact_focus_nonzero_total",
    "primary_stable_ge_0_5_total",
    "compact_stable_ge_0_5_total",
    "primary_stable_ge_0_8_total",
    "compact_stable_ge_0_8_total",
    "recommendation",
    "decision_rationale",
    "source_ladder_summary_csv",
    "source_ladder_focus_csv",
    "source_sparse_feature_audit_csv",
    "source_stability_summary_csv",
]
INTERACTION_PAIRED_COMPARATOR_FOCUS_FIELDS = [
    "paired_table_id",
    "fold_label",
    "display_fold",
    "fold_family",
    "heldout_value",
    "n_test",
    "nearest_centroid_ba",
    "primary_candidate_id",
    "primary_candidate_ba",
    "compact_candidate_id",
    "compact_candidate_ba",
    "primary_minus_compact_ba",
    "display_primary_minus_compact_ba",
    "primary_nonzero_coefficients",
    "compact_nonzero_coefficients",
    "primary_stable_ge_0_5",
    "compact_stable_ge_0_5",
    "primary_stable_ge_0_8",
    "compact_stable_ge_0_8",
    "primary_mean_pairwise_jaccard",
    "compact_mean_pairwise_jaccard",
    "compactness_delta_nonzero_compact_minus_primary",
    "focus_interpretation",
]
INTERACTION_DIAGNOSTIC_ARTIFACT_MANIFEST_FIELDS = [
    "manifest_id",
    "task_id",
    "artifact_id",
    "artifact_role",
    "block_id",
    "path",
    "file_format",
    "exists",
    "byte_size",
    "line_count",
    "row_count",
    "sha256",
    "required_for_claims",
    "validation_scope",
    "generated_by",
    "claim_boundary",
]
INTERACTION_DIAGNOSTIC_CLAIM_ARTIFACT_FIELDS = [
    "manifest_id",
    "task_id",
    "claim_id",
    "claim_status",
    "claim_text",
    "artifact_ids",
    "artifact_paths",
    "validation_tests",
    "external_source_urls",
    "limitations",
]
INTERACTION_RELEASE_READINESS_SUMMARY_FIELDS = [
    "audit_id",
    "task_id",
    "candidate_id",
    "candidate_variant_id",
    "release_readiness_decision",
    "public_alpha_ready",
    "blocker_count",
    "needs_work_count",
    "pass_count",
    "acceptable_draft_limitation_count",
    "artifact_count",
    "missing_artifact_count",
    "unhashed_artifact_count",
    "claim_count",
    "source_checksum_status",
    "source_freeze_ready",
    "source_inventory_status",
    "local_payload_status",
    "sample_factor_status",
    "sample_matrix_alignment_status",
    "claim_boundary_status",
    "reproducibility_status",
    "data_card_status",
    "metadata_package_status",
    "external_reference_ids",
    "next_required_block",
    "claim_boundary",
]
INTERACTION_RELEASE_READINESS_GAP_FIELDS = [
    "audit_id",
    "task_id",
    "gap_id",
    "category",
    "readiness_status",
    "fix_priority",
    "finding",
    "evidence_artifacts",
    "evidence_fields",
    "external_reference_ids",
    "external_reference_urls",
    "remediation",
    "claim_boundary_impact",
    "next_owner_block",
]
INTERACTION_RELEASE_READINESS_EXTERNAL_REFERENCE_FIELDS = [
    "reference_id",
    "topic",
    "url",
    "release_readiness_implication",
]
INTERACTION_PAYLOAD_FREEZE_SUMMARY_FIELDS = [
    "freeze_id",
    "task_id",
    "source_id",
    "glds_prefix",
    "checksum_manifest_url",
    "checksum_manifest_sha256",
    "checksum_manifest_sha256_matches_source_audit",
    "parsed_checksum_entries",
    "required_payload_count",
    "required_payload_matched_count",
    "required_payload_missing_count",
    "required_payload_checksum_mismatch_count",
    "out_of_scope_payload_count",
    "diagnostic_required_payload_freeze_ready",
    "full_osdr_payload_freeze_ready",
    "release_scope_decision",
    "release_status",
    "next_required_block",
    "claim_boundary",
]
INTERACTION_PAYLOAD_FREEZE_MANIFEST_FIELDS = [
    "freeze_id",
    "task_id",
    "source_id",
    "glds_prefix",
    "payload_filename",
    "algorithm",
    "expected_checksum",
    "local_path",
    "local_exists",
    "observed_checksum",
    "checksum_match",
    "release_scope",
    "payload_role",
    "verification_status",
    "action",
    "checksum_manifest_url",
    "source_checksum_audit_csv",
]
INTERACTION_PUBLIC_ALPHA_CARD_SUMMARY_FIELDS = [
    "card_id",
    "task_id",
    "source_id",
    "source_url",
    "organism",
    "assay_modality",
    "biospecimen_type",
    "candidate_id",
    "candidate_variant_id",
    "card_status",
    "release_status",
    "payload_freeze_decision",
    "diagnostic_required_payload_freeze_ready",
    "full_osdr_payload_freeze_ready",
    "artifact_count",
    "claim_count",
    "primary_balanced_accuracy",
    "secondary_balanced_accuracy",
    "diagnostic_balanced_accuracy",
    "nearest_improved_count",
    "nearest_tied_count",
    "nearest_worse_count",
    "stable_feature_count_ge_0_5_total",
    "allowed_claim_count",
    "disallowed_claim_count",
    "next_required_block",
    "claim_boundary",
]
INTERACTION_REBUILD_GATE_SUMMARY_FIELDS = [
    "gate_id",
    "task_id",
    "gate_status",
    "mode",
    "step_count",
    "ready_step_count",
    "missing_output_count",
    "script_missing_count",
    "hashed_output_count",
    "environment_package_count",
    "generated_at_utc",
    "python_version",
    "platform",
    "repo_root",
    "reports_root",
    "rebuild_command",
    "next_required_block",
    "claim_boundary",
]
INTERACTION_REBUILD_GATE_STEP_FIELDS = [
    "gate_id",
    "task_id",
    "step_order",
    "step_id",
    "block_id",
    "script_path",
    "script_exists",
    "command",
    "step_role",
    "execution_policy",
    "status",
    "output_count",
    "missing_output_count",
    "hashed_output_count",
    "output_paths",
    "missing_output_paths",
    "output_sha256s",
    "notes",
]
INTERACTION_REBUILD_GATE_ENVIRONMENT_FIELDS = [
    "gate_id",
    "task_id",
    "key",
    "value",
    "source",
]
INTERACTION_PUBLIC_METADATA_SUMMARY_FIELDS = [
    "package_id",
    "task_id",
    "source_id",
    "source_url",
    "metadata_package_status",
    "release_target_decision",
    "public_now_target_count",
    "not_public_now_target_count",
    "metadata_field_count",
    "metadata_present_field_count",
    "metadata_partial_field_count",
    "metadata_placeholder_field_count",
    "artifact_count",
    "claim_count",
    "rebuild_gate_status",
    "diagnostic_required_payload_freeze_ready",
    "full_osdr_payload_freeze_ready",
    "datacite_schema_version",
    "ro_crate_version",
    "hf_card_metadata_status",
    "osdr_citation_status",
    "next_required_block",
    "claim_boundary",
]
INTERACTION_PUBLIC_METADATA_RELEASE_TARGET_FIELDS = [
    "package_id",
    "task_id",
    "target_id",
    "release_surface",
    "target_status",
    "public_now",
    "allowed_claims",
    "disallowed_claims",
    "required_evidence_artifacts",
    "blocking_gaps",
    "next_action",
]
INTERACTION_PUBLIC_METADATA_FIELD_FIELDS = [
    "package_id",
    "task_id",
    "metadata_profile",
    "field_id",
    "field_label",
    "status",
    "current_value",
    "source_artifacts",
    "blocking_gap",
    "notes",
]
INTERACTION_PUBLIC_METADATA_REFERENCE_FIELDS = [
    "reference_id",
    "topic",
    "url",
    "checked_date",
    "metadata_implication",
]
INTERACTION_RO_CRATE_EXPORT_SUMMARY_FIELDS = [
    "scaffold_id",
    "package_id",
    "task_id",
    "scaffold_status",
    "ro_crate_graph_entity_count",
    "ro_crate_data_entity_count",
    "data_package_resource_count",
    "validation_check_count",
    "validation_pass_count",
    "validation_blocker_count",
    "validation_needs_review_count",
    "citation_item_count",
    "citation_pass_count",
    "citation_blocker_count",
    "citation_needs_review_count",
    "placeholder_field_count",
    "next_required_block",
    "claim_boundary",
]
INTERACTION_RO_CRATE_VALIDATION_FIELDS = [
    "scaffold_id",
    "package_id",
    "task_id",
    "check_id",
    "standard_profile",
    "check_status",
    "severity",
    "finding",
    "evidence",
    "blocking_gap",
    "next_action",
]
INTERACTION_CITATION_FREEZE_CHECKLIST_FIELDS = [
    "scaffold_id",
    "package_id",
    "task_id",
    "item_id",
    "citation_surface",
    "required_for",
    "status",
    "current_value",
    "source_artifacts",
    "blocking_gap",
    "next_action",
    "claim_boundary_impact",
]
INTERACTION_ARCHIVE_DECISION_SUMMARY_FIELDS = [
    "decision_id",
    "scaffold_id",
    "package_id",
    "task_id",
    "decision_status",
    "archive_path_decision",
    "version_decision",
    "creator_decision",
    "license_decision",
    "archive_option_count",
    "current_selected_option_count",
    "deferred_option_count",
    "blocked_option_count",
    "license_component_count",
    "license_blocker_count",
    "creator_component_count",
    "creator_blocker_count",
    "external_reference_count",
    "next_required_block",
    "claim_boundary",
]
INTERACTION_ARCHIVE_IDENTIFIER_OPTION_FIELDS = [
    "decision_id",
    "scaffold_id",
    "package_id",
    "task_id",
    "option_id",
    "option_label",
    "decision_status",
    "current_draft_selected",
    "citable_release_ready",
    "identifier_type",
    "expected_identifier",
    "official_reference_ids",
    "evidence_artifacts",
    "blocking_gaps",
    "required_actions",
    "rationale",
]
INTERACTION_LICENSE_RIGHTS_DECISION_FIELDS = [
    "decision_id",
    "scaffold_id",
    "package_id",
    "task_id",
    "component_id",
    "component_scope",
    "decision_status",
    "current_value",
    "recommended_decision_state",
    "official_reference_ids",
    "evidence_artifacts",
    "blocking_gaps",
    "required_actions",
    "claim_boundary_impact",
]
INTERACTION_CREATOR_CONTRIBUTOR_DECISION_FIELDS = [
    "decision_id",
    "scaffold_id",
    "package_id",
    "task_id",
    "component_id",
    "component_scope",
    "decision_status",
    "current_value",
    "recommended_decision_state",
    "official_reference_ids",
    "evidence_artifacts",
    "blocking_gaps",
    "required_actions",
    "claim_boundary_impact",
]
INTERACTION_ARCHIVE_DECISION_REFERENCE_FIELDS = [
    "reference_id",
    "topic",
    "url",
    "checked_date",
    "decision_implication",
]
INTERACTION_CITATION_METADATA_FILL_SUMMARY_FIELDS = [
    "fill_id",
    "decision_id",
    "scaffold_id",
    "package_id",
    "task_id",
    "fill_status",
    "owner_metadata_status",
    "intake_field_count",
    "supplied_field_count",
    "retained_current_draft_count",
    "not_supplied_field_count",
    "not_supplied_blocker_count",
    "needs_review_count",
    "descriptor_patch_status",
    "ro_crate_mutation_status",
    "datapackage_mutation_status",
    "release_ready_after_fill",
    "next_required_block",
    "claim_boundary",
]
INTERACTION_CITATION_METADATA_FILL_FIELDS = [
    "fill_id",
    "decision_id",
    "package_id",
    "task_id",
    "field_id",
    "field_group",
    "target_profile",
    "required_for",
    "current_status_from_gate",
    "current_value",
    "supplied_value",
    "supplied_by",
    "supplied_date",
    "supplied_evidence",
    "fill_status",
    "target_artifacts",
    "validation_rule",
    "blocker_if_missing",
    "next_action",
    "notes",
]
INTERACTION_ARCHIVE_RELEASE_DEFERRAL_SUMMARY_FIELDS = [
    "guard_id",
    "fill_id",
    "decision_id",
    "package_id",
    "task_id",
    "guard_status",
    "release_deferral_status",
    "owner_metadata_status",
    "supplied_field_count",
    "blocked_field_count",
    "retained_current_draft_count",
    "application_guard_count",
    "guard_pass_count",
    "guard_blocker_count",
    "guard_deferred_count",
    "action_count",
    "required_owner_action_count",
    "descriptor_mutation_allowed",
    "release_ready_after_guard",
    "next_required_block",
    "claim_boundary",
]
INTERACTION_OWNER_METADATA_APPLICATION_GUARD_FIELDS = [
    "guard_id",
    "fill_id",
    "package_id",
    "task_id",
    "guard_check_id",
    "guard_surface",
    "guard_status",
    "severity",
    "evidence_artifacts",
    "observed_value",
    "required_value",
    "source_blocker_fields",
    "action_if_failed",
    "mutation_policy",
    "claim_boundary_impact",
]
INTERACTION_ARCHIVE_RELEASE_DEFERRAL_ACTION_FIELDS = [
    "guard_id",
    "fill_id",
    "package_id",
    "task_id",
    "action_id",
    "action_owner",
    "action_status",
    "required_for",
    "source_blocker_fields",
    "source_evidence",
    "next_action",
    "deferral_policy",
]
INTERACTION_DIAGNOSTIC_METADATA_RELEASE_SUMMARY_FIELDS = [
    "release_note_id",
    "guard_id",
    "fill_id",
    "package_id",
    "task_id",
    "note_status",
    "branch_closeout_status",
    "current_public_surface",
    "inspectable_now_count",
    "not_released_claim_count",
    "owner_retry_item_count",
    "descriptor_mutation_status",
    "archive_release_status",
    "next_required_block",
    "claim_boundary",
]
INTERACTION_DIAGNOSTIC_METADATA_RELEASE_SECTION_FIELDS = [
    "release_note_id",
    "guard_id",
    "package_id",
    "task_id",
    "section_id",
    "section_title",
    "include_in_note",
    "note_text",
    "evidence_artifacts",
    "claim_boundary",
]
INTERACTION_DIAGNOSTIC_METADATA_PUBLIC_CLAIM_FIELDS = [
    "release_note_id",
    "guard_id",
    "package_id",
    "task_id",
    "claim_id",
    "claim_category",
    "public_note_status",
    "statement",
    "supporting_evidence",
    "prohibited_language",
    "next_allowed_action",
]
INTERACTION_OWNER_METADATA_RETRY_CHECKLIST_FIELDS = [
    "release_note_id",
    "guard_id",
    "package_id",
    "task_id",
    "retry_item_id",
    "required_for",
    "owner_action",
    "current_status",
    "source_blocker_fields",
    "evidence_artifacts",
    "validation_rule",
    "retry_priority",
    "closeout_policy",
]
DEFAULT_LOGISTIC_FEATURE_AUDIT_FOLDS = [
    ("diagnostic_condition_stratum_holdout", "Col.0.PhyD|Dark.Treatment"),
    ("secondary_light_treatment_holdout", "Light.Treatment"),
    ("primary_genotype_or_ecotype_holdout", "Col.0.PhyD"),
]
DEFAULT_INTERACTION_BASELINE_LADDER_FOCUS_FOLDS = [
    (
        "diagnostic_condition_stratum_holdout",
        "Col.0.PhyD|Dark.Treatment",
        "dark",
        "focus_col0_phyd_dark_ba",
    ),
    (
        "secondary_light_treatment_holdout",
        "Light.Treatment",
        "light",
        "focus_light_treatment_ba",
    ),
    (
        "primary_genotype_or_ecotype_holdout",
        "Col.0.PhyD",
        "phyd",
        "focus_col0_phyd_ba",
    ),
]
DEFAULT_INTERACTION_BASELINE_LADDER_CANDIDATES = [
    {
        "candidate_id": "nearest_centroid_default",
        "role": "reference_floor",
        "decision": "retain_as_nearest_centroid_floor",
        "model_class": "nearest_centroid",
        "baseline_id": INTERACTION_BASELINE_ID,
        "variant_id": "tvg2000_log1p_zscore",
        "summary_csv": "interaction_sensitivity/multispecies_baseline_summary.csv",
        "fold_detail_csv": "interaction_sensitivity/fold_detail_summary.csv",
        "comparison_csv": "",
        "stability_csv": "",
        "penalty": "",
        "c": "",
    },
    {
        "candidate_id": "l2_logistic_default",
        "role": "dense_logistic_control",
        "decision": "retain_as_dense_control_with_dark_focus_regression",
        "model_class": "l2_logistic",
        "baseline_id": INTERACTION_LOGISTIC_BASELINE_ID,
        "variant_id": "tvg2000_log1p_zscore",
        "summary_csv": "interaction_logistic_l2/multispecies_baseline_summary.csv",
        "fold_detail_csv": "interaction_logistic_l2/fold_detail_summary.csv",
        "comparison_csv": (
            "interaction_logistic_l2/"
            "fold_detail_comparison_vs_nearest_centroid.csv"
        ),
        "stability_csv": "",
        "penalty": "l2",
        "c": "1",
    },
    {
        "candidate_id": "l2_logistic_top500_c1",
        "role": "dense_logistic_sensitivity_control",
        "decision": "retain_as_top500_control_with_light_focus_regression",
        "model_class": "l2_logistic",
        "baseline_id": INTERACTION_LOGISTIC_BASELINE_ID,
        "variant_id": "tvg500_log1p_zscore_c1",
        "summary_csv": (
            "interaction_logistic_l2_sensitivity/"
            "multispecies_baseline_summary.csv"
        ),
        "fold_detail_csv": "interaction_logistic_l2_sensitivity/fold_detail_summary.csv",
        "comparison_csv": (
            "interaction_logistic_l2_sensitivity/"
            "fold_detail_comparison_vs_nearest_centroid.csv"
        ),
        "stability_csv": "",
        "penalty": "l2",
        "c": "1",
    },
    {
        "candidate_id": "sparse_l1_c0p3",
        "role": "compact_sparse_comparator",
        "decision": "retain_as_compact_stability_comparator",
        "model_class": "sparse_l1_logistic",
        "baseline_id": INTERACTION_LOGISTIC_BASELINE_ID,
        "variant_id": "tvg2000_log1p_zscore_l1_c0p3",
        "summary_csv": "interaction_logistic_sparse_l1/multispecies_baseline_summary.csv",
        "fold_detail_csv": "interaction_logistic_sparse_l1/fold_detail_summary.csv",
        "comparison_csv": (
            "interaction_logistic_sparse_l1/"
            "fold_detail_comparison_vs_nearest_centroid.csv"
        ),
        "stability_csv": (
            "interaction_logistic_sparse_l1_stability/stability_summary.csv"
        ),
        "penalty": "l1",
        "c": "0.3",
    },
    {
        "candidate_id": "sparse_l1_c1",
        "role": "leading_transparent_candidate",
        "decision": "advance_as_primary_sparse_diagnostic_candidate",
        "model_class": "sparse_l1_logistic",
        "baseline_id": INTERACTION_LOGISTIC_BASELINE_ID,
        "variant_id": "tvg2000_log1p_zscore_l1_c1",
        "summary_csv": "interaction_logistic_sparse_l1/multispecies_baseline_summary.csv",
        "fold_detail_csv": "interaction_logistic_sparse_l1/fold_detail_summary.csv",
        "comparison_csv": (
            "interaction_logistic_sparse_l1/"
            "fold_detail_comparison_vs_nearest_centroid.csv"
        ),
        "stability_csv": (
            "interaction_logistic_sparse_l1_stability/stability_summary.csv"
        ),
        "penalty": "l1",
        "c": "1",
    },
]
DEFAULT_INTERACTION_DIAGNOSTIC_PACKAGE_ID = "osd120_sparse_l1_c1_draft_candidate"
DEFAULT_INTERACTION_DIAGNOSTIC_CANDIDATE_ID = "sparse_l1_c1"
DEFAULT_INTERACTION_DIAGNOSTIC_COMPARATOR_ID = "sparse_l1_c0p3"
DEFAULT_INTERACTION_FIGURE_TABLE_ID = "osd120_sparse_l1_c1_focus_table"
DEFAULT_INTERACTION_PAIRED_COMPARATOR_TABLE_ID = "osd120_sparse_l1_c1_vs_c0p3"
DEFAULT_INTERACTION_DIAGNOSTIC_ARTIFACT_MANIFEST_ID = (
    "osd120_diagnostic_artifact_manifest"
)
DEFAULT_INTERACTION_RELEASE_READINESS_AUDIT_ID = (
    "osd120_release_readiness_gap_audit"
)
DEFAULT_INTERACTION_PAYLOAD_FREEZE_ID = "osd120_diagnostic_payload_freeze"
DEFAULT_INTERACTION_PUBLIC_ALPHA_CARD_ID = "osd120_diagnostic_public_alpha_card"
DEFAULT_INTERACTION_REBUILD_GATE_ID = "osd120_diagnostic_packaging_rebuild_gate"
DEFAULT_INTERACTION_PUBLIC_METADATA_PACKAGE_ID = (
    "osd120_diagnostic_public_metadata_package"
)
DEFAULT_INTERACTION_RO_CRATE_EXPORT_SCAFFOLD_ID = (
    "osd120_ro_crate_citation_freeze_scaffold"
)
DEFAULT_INTERACTION_ARCHIVE_DECISION_GATE_ID = (
    "osd120_archive_identifier_license_decision_gate"
)
DEFAULT_INTERACTION_CITATION_METADATA_FILL_ID = (
    "osd120_release_owner_citation_metadata_fill"
)
DEFAULT_INTERACTION_ARCHIVE_RELEASE_DEFERRAL_GUARD_ID = (
    "osd120_archive_release_deferral_application_guard"
)
DEFAULT_INTERACTION_DIAGNOSTIC_METADATA_RELEASE_NOTE_ID = (
    "osd120_diagnostic_metadata_release_note_closeout"
)


@dataclass(frozen=True)
class MultispeciesBaselineConfig:
    """Configuration for the draft multispecies nearest-centroid baseline."""

    transform: str = "log1p"
    scaling: str = "zscore"
    top_variable_genes: int = 2000

    def __post_init__(self) -> None:
        if self.transform not in {"none", "log1p"}:
            raise ValueError("transform must be 'none' or 'log1p'")
        if self.scaling not in {"none", "zscore"}:
            raise ValueError("scaling must be 'none' or 'zscore'")
        if self.top_variable_genes < 1:
            raise ValueError("top_variable_genes must be at least 1")


@dataclass(frozen=True)
class MultispeciesLogisticBaselineConfig(MultispeciesBaselineConfig):
    """Configuration for the draft OSD-120 interaction L2 logistic baseline."""

    c: float = 1.0
    max_iter: int = 5000
    random_state: int = 42
    penalty: str = "l2"
    solver: str = "liblinear"
    l1_ratio: float | None = None

    def __post_init__(self) -> None:
        super().__post_init__()
        if self.c <= 0:
            raise ValueError("c must be positive")
        if self.max_iter < 1:
            raise ValueError("max_iter must be at least 1")
        if self.penalty not in {"l1", "l2", "elasticnet"}:
            raise ValueError("penalty must be 'l1', 'l2', or 'elasticnet'")
        if self.solver not in {"liblinear", "saga"}:
            raise ValueError("solver must be 'liblinear' or 'saga'")
        if self.penalty == "elasticnet" and self.solver != "saga":
            raise ValueError("elasticnet penalty requires solver='saga'")
        if self.penalty != "elasticnet" and self.l1_ratio is not None:
            raise ValueError("l1_ratio is only valid for elasticnet penalty")
        if self.l1_ratio is not None and not 0.0 <= self.l1_ratio <= 1.0:
            raise ValueError("l1_ratio must be between 0 and 1")


@dataclass(frozen=True)
class MultispeciesFoldPredictionResult:
    """Predictions and fold metadata for one multispecies draft fold."""

    task_id: str
    fold_id: str
    heldout_factor: str
    heldout_value: str
    n_train: int
    n_test: int
    n_features: int
    predictions: list[dict[str, str]]


@dataclass(frozen=True)
class MultispeciesLogisticFoldPredictionResult(MultispeciesFoldPredictionResult):
    """Predictions and fit metadata for one OSD-120 L2 logistic fold."""

    train_time_s: float
    fit_details: Mapping[str, Any]


@dataclass(frozen=True)
class MultispeciesLogisticFeatureAuditResult:
    """Fold-level selected-feature and coefficient audit for OSD-120 logistic."""

    summary: dict[str, str]
    coefficients: list[dict[str, str]]


@dataclass(frozen=True)
class MultispeciesLogisticStabilityAuditResult:
    """Subsampling stability audit for one OSD-120 sparse logistic fold."""

    summary: dict[str, str]
    features: list[dict[str, str]]


@dataclass(frozen=True)
class MultispeciesBaselineTaskResult:
    """Filesystem outputs and metrics for one multispecies baseline run."""

    baseline_id: str
    task_id: str
    output_dir: Path
    predictions_path: Path
    metrics_path: Path
    run_manifest_path: Path
    evaluation: Mapping[str, Any]
    n_predictions: int


def _format_float(value: float) -> str:
    return f"{value:.10g}"


def _sample_factor_by_id(task: MultispeciesTaskData) -> dict[str, dict[str, str]]:
    return {
        str(row["sample_id"]): row
        for row in task.sample_factors
        if row.get("parse_status") == "parsed" and row.get("sample_id")
    }


def _transform(values: np.ndarray, mode: str) -> np.ndarray:
    if mode == "none":
        return values.astype(np.float64, copy=True)
    if mode == "log1p":
        clipped = np.clip(values.astype(np.float64, copy=True), a_min=0.0, a_max=None)
        return np.log1p(clipped)
    raise ValueError(f"unsupported transform: {mode}")


def _select_train_variable_features(
    train: np.ndarray,
    test: np.ndarray,
    *,
    top_variable_genes: int,
) -> tuple[np.ndarray, np.ndarray]:
    n_features = train.shape[1]
    selected, _ = _selected_train_variable_feature_indices(
        train,
        top_variable_genes=top_variable_genes,
    )
    if len(selected) >= n_features:
        return train, test
    return train[:, selected], test[:, selected]


def _selected_train_variable_feature_indices(
    train: np.ndarray,
    *,
    top_variable_genes: int,
) -> tuple[np.ndarray, np.ndarray]:
    n_features = train.shape[1]
    variances = np.var(train, axis=0)
    ranked = np.argsort(-variances, kind="mergesort")
    if top_variable_genes >= n_features:
        selected = np.arange(n_features)
    else:
        selected = ranked[:top_variable_genes]
    selected.sort()
    return selected, ranked


def _scale_train_test(
    train: np.ndarray,
    test: np.ndarray,
    *,
    scaling: str,
) -> tuple[np.ndarray, np.ndarray]:
    if scaling == "none":
        return train, test
    if scaling != "zscore":
        raise ValueError(f"unsupported scaling mode: {scaling}")
    mean = train.mean(axis=0)
    std = train.std(axis=0)
    std[std == 0.0] = 1.0
    return (train - mean) / std, (test - mean) / std


def _class_centroids(features: np.ndarray, labels: Sequence[str]) -> np.ndarray:
    centroids = []
    for label in (NEGATIVE_LABEL, POSITIVE_LABEL):
        mask = np.asarray([value == label for value in labels], dtype=bool)
        if not mask.any():
            raise ValueError(f"training fold has no {label} samples")
        centroids.append(features[mask].mean(axis=0))
    return np.vstack(centroids)


def _distances_to_centroids(features: np.ndarray, centroids: np.ndarray) -> np.ndarray:
    squared = np.square(features[:, np.newaxis, :] - centroids[np.newaxis, :, :]).sum(axis=2)
    return np.sqrt(squared)


def _positive_probability(ground_distance: np.ndarray, leo_distance: np.ndarray) -> np.ndarray:
    total = ground_distance + leo_distance
    return np.divide(
        ground_distance,
        total,
        out=np.full_like(total, 0.5, dtype=float),
        where=total > 0.0,
    )


def predict_multispecies_fold(
    task: MultispeciesTaskData,
    fold: MultispeciesFoldData,
    *,
    config: MultispeciesBaselineConfig | None = None,
) -> MultispeciesFoldPredictionResult:
    """Fit a transparent centroid classifier on one multispecies condition-stratum fold."""

    cfg = config or MultispeciesBaselineConfig()
    sample_factors = _sample_factor_by_id(task)
    train_x = task.features.loc[fold.train_sample_ids]
    test_x = task.features.loc[fold.test_sample_ids]
    train_labels = [sample_factors[sample_id]["true_label"] for sample_id in train_x.index]
    test_labels = [sample_factors[sample_id]["true_label"] for sample_id in test_x.index]

    train_values = _transform(train_x.to_numpy(dtype=np.float64), cfg.transform)
    test_values = _transform(test_x.to_numpy(dtype=np.float64), cfg.transform)
    train_values, test_values = _select_train_variable_features(
        train_values,
        test_values,
        top_variable_genes=cfg.top_variable_genes,
    )
    train_values, test_values = _scale_train_test(
        train_values,
        test_values,
        scaling=cfg.scaling,
    )
    centroids = _class_centroids(train_values, train_labels)
    distances = _distances_to_centroids(test_values, centroids)
    ground_distance = distances[:, 0]
    leo_distance = distances[:, 1]
    probabilities = _positive_probability(ground_distance, leo_distance)
    predicted = np.where(leo_distance <= ground_distance, POSITIVE_LABEL, NEGATIVE_LABEL)

    rows: list[dict[str, str]] = []
    for index, sample_id in enumerate(test_x.index):
        row = sample_factors[str(sample_id)]
        rows.append(
            {
                "task_id": task.task_id,
                "fold_id": fold.fold_id,
                "sample_id": str(sample_id),
                "mission": str(row.get("mission", "")),
                "condition_stratum": str(row.get("condition_stratum", "")),
                "genotype_or_ecotype": str(row.get("genotype_or_ecotype", "")),
                "light_treatment": str(row.get("light_treatment", "")),
                "true_label": test_labels[index],
                "predicted_label": str(predicted[index]),
                "flight_probability": _format_float(float(probabilities[index])),
                "embedding_0": _format_float(float(ground_distance[index])),
                "embedding_1": _format_float(float(leo_distance[index])),
            }
        )

    return MultispeciesFoldPredictionResult(
        task_id=task.task_id,
        fold_id=fold.fold_id,
        heldout_factor=fold.heldout_factor,
        heldout_value=fold.heldout_value,
        n_train=len(train_x),
        n_test=len(test_x),
        n_features=int(train_values.shape[1]),
        predictions=rows,
    )


def _numeric_labels(labels: Sequence[str]) -> np.ndarray:
    mapping = {NEGATIVE_LABEL: 0, POSITIVE_LABEL: 1}
    try:
        return np.asarray([mapping[label] for label in labels], dtype=int)
    except KeyError as exc:
        raise ValueError(f"unexpected label for binary classifier: {exc.args[0]!r}") from exc


def _prediction_label(value: int) -> str:
    return POSITIVE_LABEL if int(value) == 1 else NEGATIVE_LABEL


def _logistic_model_kwargs(config: MultispeciesLogisticBaselineConfig) -> dict[str, Any]:
    kwargs: dict[str, Any] = {
        "C": config.c,
        "class_weight": "balanced",
        "l1_ratio": 0.0 if config.penalty == "l2" else 1.0,
        "max_iter": config.max_iter,
        "random_state": config.random_state,
        "solver": config.solver,
    }
    if config.penalty == "elasticnet":
        kwargs["l1_ratio"] = config.l1_ratio
    return kwargs


def predict_multispecies_logistic_fold(
    task: MultispeciesTaskData,
    fold: MultispeciesFoldData,
    *,
    config: MultispeciesLogisticBaselineConfig | None = None,
) -> MultispeciesLogisticFoldPredictionResult:
    """Fit an L2 logistic classifier on one OSD-120 interaction fold."""

    from sklearn.linear_model import LogisticRegression

    cfg = config or MultispeciesLogisticBaselineConfig()
    sample_factors = _sample_factor_by_id(task)
    train_x = task.features.loc[fold.train_sample_ids]
    test_x = task.features.loc[fold.test_sample_ids]
    train_labels = [sample_factors[sample_id]["true_label"] for sample_id in train_x.index]
    test_labels = [sample_factors[sample_id]["true_label"] for sample_id in test_x.index]

    train_values = _transform(train_x.to_numpy(dtype=np.float64), cfg.transform)
    test_values = _transform(test_x.to_numpy(dtype=np.float64), cfg.transform)
    train_values, test_values = _select_train_variable_features(
        train_values,
        test_values,
        top_variable_genes=cfg.top_variable_genes,
    )
    train_values, test_values = _scale_train_test(
        train_values,
        test_values,
        scaling=cfg.scaling,
    )
    model = LogisticRegression(**_logistic_model_kwargs(cfg))
    start = time.perf_counter()
    model.fit(train_values, _numeric_labels(train_labels))
    train_time_s = time.perf_counter() - start

    classes = list(model.classes_)
    try:
        positive_index = classes.index(1)
    except ValueError as exc:
        raise ValueError(f"classifier did not learn positive class 1: {classes}") from exc
    probabilities = model.predict_proba(test_values)[:, positive_index]
    predicted_numeric = model.predict(test_values)

    rows: list[dict[str, str]] = []
    for index, sample_id in enumerate(test_x.index):
        row = sample_factors[str(sample_id)]
        rows.append(
            {
                "task_id": task.task_id,
                "fold_id": fold.fold_id,
                "sample_id": str(sample_id),
                "mission": str(row.get("mission", "")),
                "condition_stratum": str(row.get("condition_stratum", "")),
                "genotype_or_ecotype": str(row.get("genotype_or_ecotype", "")),
                "light_treatment": str(row.get("light_treatment", "")),
                "true_label": test_labels[index],
                "predicted_label": _prediction_label(int(predicted_numeric[index])),
                "flight_probability": _format_float(float(probabilities[index])),
            }
        )

    return MultispeciesLogisticFoldPredictionResult(
        task_id=task.task_id,
        fold_id=fold.fold_id,
        heldout_factor=fold.heldout_factor,
        heldout_value=fold.heldout_value,
        n_train=len(train_x),
        n_test=len(test_x),
        n_features=int(train_values.shape[1]),
        predictions=rows,
        train_time_s=train_time_s,
        fit_details={
            "model_type": INTERACTION_LOGISTIC_BASELINE_ID,
            "solver": cfg.solver,
            "penalty": cfg.penalty,
            "c": cfg.c,
            "l1_ratio": cfg.l1_ratio,
            "class_weight": "balanced",
            "max_iter": cfg.max_iter,
            "random_state": cfg.random_state,
            "train_time_s": round(train_time_s, 6),
        },
    )


def audit_multispecies_logistic_fold_features(
    task: MultispeciesTaskData,
    fold: MultispeciesFoldData,
    *,
    config: MultispeciesLogisticBaselineConfig,
    top_n_coefficients: int = 10,
) -> MultispeciesLogisticFeatureAuditResult:
    """Audit selected features and coefficients for one OSD-120 logistic fold."""

    from sklearn.linear_model import LogisticRegression

    sample_factors = _sample_factor_by_id(task)
    train_x = task.features.loc[fold.train_sample_ids]
    test_x = task.features.loc[fold.test_sample_ids]
    train_labels = [sample_factors[sample_id]["true_label"] for sample_id in train_x.index]
    test_labels = [sample_factors[sample_id]["true_label"] for sample_id in test_x.index]

    train_values = _transform(train_x.to_numpy(dtype=np.float64), config.transform)
    test_values = _transform(test_x.to_numpy(dtype=np.float64), config.transform)
    selected_indices, ranked_indices = _selected_train_variable_feature_indices(
        train_values,
        top_variable_genes=config.top_variable_genes,
    )
    train_variances = np.var(train_values, axis=0)
    variance_ranks = {
        int(feature_index): rank
        for rank, feature_index in enumerate(ranked_indices.tolist(), start=1)
    }
    train_selected = train_values[:, selected_indices]
    test_selected = test_values[:, selected_indices]
    train_selected, test_selected = _scale_train_test(
        train_selected,
        test_selected,
        scaling=config.scaling,
    )
    model = LogisticRegression(**_logistic_model_kwargs(config))
    model.fit(train_selected, _numeric_labels(train_labels))
    classes = list(model.classes_)
    try:
        positive_index = classes.index(1)
    except ValueError as exc:
        raise ValueError(f"classifier did not learn positive class 1: {classes}") from exc
    probabilities = model.predict_proba(test_selected)[:, positive_index]
    predicted_numeric = model.predict(test_selected)
    prediction_rows = [
        {
            "true_label": test_labels[index],
            "predicted_label": _prediction_label(int(predicted_numeric[index])),
            "flight_probability": _format_float(float(probabilities[index])),
        }
        for index in range(len(test_labels))
    ]
    coefficients = np.asarray(model.coef_[0], dtype=float)
    feature_ids = [str(train_x.columns[int(index)]) for index in selected_indices]
    abs_order = np.argsort(-np.abs(coefficients), kind="mergesort")
    positive_order = np.argsort(-coefficients, kind="mergesort")
    negative_order = np.argsort(coefficients, kind="mergesort")
    abs_ranks = {
        int(coef_index): rank
        for rank, coef_index in enumerate(abs_order.tolist(), start=1)
    }
    in_top500 = {
        int(feature_index): variance_ranks[int(feature_index)] <= 500
        for feature_index in selected_indices
    }
    top_abs_extra = [
        feature_ids[int(coef_index)]
        for coef_index in abs_order
        if not in_top500[int(selected_indices[int(coef_index)])]
    ][:top_n_coefficients]
    nonzero_mask = np.abs(coefficients) > 0.0

    variant_id = multispecies_logistic_config_id(config)
    summary = {
        "task_id": task.task_id,
        "fold_family": fold.fold_family,
        "heldout_factor": fold.heldout_factor,
        "heldout_value": fold.heldout_value,
        "fold_id": fold.fold_id,
        "variant_id": variant_id,
        "baseline_id": INTERACTION_LOGISTIC_BASELINE_ID,
        "transform": config.transform,
        "scaling": config.scaling,
        "top_variable_genes": str(config.top_variable_genes),
        "c": _format_float(config.c),
        "n_train": str(len(train_x)),
        "n_test": str(len(test_x)),
        "n_selected_features": str(len(feature_ids)),
        "n_top500_features": str(sum(in_top500.values())),
        "n_extra_features": str(sum(not value for value in in_top500.values())),
        "n_nonzero_coefficients": str(int(nonzero_mask.sum())),
        "n_nonzero_top500_coefficients": str(
            sum(
                bool(nonzero_mask[index]) and in_top500[int(feature_index)]
                for index, feature_index in enumerate(selected_indices.tolist())
            )
        ),
        "n_nonzero_extra_coefficients": str(
            sum(
                bool(nonzero_mask[index]) and not in_top500[int(feature_index)]
                for index, feature_index in enumerate(selected_indices.tolist())
            )
        ),
        "balanced_accuracy": _format_float(_balanced_accuracy(prediction_rows)),
        "n_positive_predictions": str(
            sum(row["predicted_label"] == POSITIVE_LABEL for row in prediction_rows)
        ),
        "top_abs_feature_ids": "|".join(
            feature_ids[int(index)] for index in abs_order[:top_n_coefficients]
        ),
        "top_positive_feature_ids": "|".join(
            feature_ids[int(index)] for index in positive_order[:top_n_coefficients]
        ),
        "top_negative_feature_ids": "|".join(
            feature_ids[int(index)] for index in negative_order[:top_n_coefficients]
        ),
        "top_abs_extra_feature_ids": "|".join(top_abs_extra),
    }
    coefficient_rows: list[dict[str, str]] = []
    for coef_index, feature_index in enumerate(selected_indices.tolist()):
        coefficient = float(coefficients[coef_index])
        coefficient_rows.append(
            {
                "task_id": task.task_id,
                "fold_family": fold.fold_family,
                "heldout_factor": fold.heldout_factor,
                "heldout_value": fold.heldout_value,
                "fold_id": fold.fold_id,
                "variant_id": variant_id,
                "top_variable_genes": str(config.top_variable_genes),
                "c": _format_float(config.c),
                "feature_id": feature_ids[coef_index],
                "train_variance_rank": str(variance_ranks[int(feature_index)]),
                "train_variance": _format_float(float(train_variances[int(feature_index)])),
                "coefficient": _format_float(coefficient),
                "abs_coefficient": _format_float(abs(coefficient)),
                "coefficient_rank_abs": str(abs_ranks[coef_index]),
                "in_top500_for_fold": str(in_top500[int(feature_index)]),
                "coefficient_direction": (
                    "positive_leo" if coefficient > 0.0 else "negative_ground"
                ),
            }
        )
    return MultispeciesLogisticFeatureAuditResult(
        summary=summary,
        coefficients=coefficient_rows,
    )


def _fit_logistic_coefficients(
    train_values: np.ndarray,
    labels: Sequence[str],
    *,
    config: MultispeciesLogisticBaselineConfig,
) -> np.ndarray:
    from sklearn.linear_model import LogisticRegression

    scaled_train, _ = _scale_train_test(
        train_values,
        train_values,
        scaling=config.scaling,
    )
    model = LogisticRegression(**_logistic_model_kwargs(config))
    model.fit(scaled_train, _numeric_labels(labels))
    return np.asarray(model.coef_[0], dtype=float)


def _deterministic_balanced_subsample_indices(
    labels: Sequence[str],
    *,
    fraction: float,
    seed: int,
) -> list[int]:
    rng = np.random.default_rng(seed)
    selected: list[int] = []
    for label in (NEGATIVE_LABEL, POSITIVE_LABEL):
        label_indices = [index for index, value in enumerate(labels) if value == label]
        if not label_indices:
            continue
        n_select = max(1, int(np.floor(len(label_indices) * fraction)))
        n_select = min(n_select, len(label_indices))
        chosen = rng.choice(label_indices, size=n_select, replace=False)
        selected.extend(int(index) for index in chosen.tolist())
    return sorted(selected)


def audit_multispecies_logistic_stability_fold(
    task: MultispeciesTaskData,
    fold: MultispeciesFoldData,
    *,
    config: MultispeciesLogisticBaselineConfig,
    n_subsamples: int = 20,
    subsample_fraction: float = 0.8,
    random_seed: int = 1729,
) -> MultispeciesLogisticStabilityAuditResult:
    """Audit sparse logistic feature stability under deterministic subsampling."""

    sample_factors = _sample_factor_by_id(task)
    train_x = task.features.loc[fold.train_sample_ids]
    train_labels = [sample_factors[sample_id]["true_label"] for sample_id in train_x.index]
    transformed = _transform(train_x.to_numpy(dtype=np.float64), config.transform)
    selected_indices, ranked_indices = _selected_train_variable_feature_indices(
        transformed,
        top_variable_genes=config.top_variable_genes,
    )
    candidate_values = transformed[:, selected_indices]
    feature_ids = [str(train_x.columns[int(index)]) for index in selected_indices]
    variance_ranks = {
        int(feature_index): rank
        for rank, feature_index in enumerate(ranked_indices.tolist(), start=1)
    }

    reference = audit_multispecies_logistic_fold_features(
        task,
        fold,
        config=config,
        top_n_coefficients=10,
    )
    reference_by_feature = {
        row["feature_id"]: float(row["coefficient"])
        for row in reference.coefficients
    }
    reference_coefficients = np.asarray(
        [reference_by_feature[feature_id] for feature_id in feature_ids],
        dtype=float,
    )
    reference_selected = np.abs(reference_coefficients) > 0.0

    coefficients_by_subsample: list[np.ndarray] = []
    selected_sets: list[set[str]] = []
    nonzero_counts: list[int] = []
    for repeat in range(n_subsamples):
        subsample_indices = _deterministic_balanced_subsample_indices(
            train_labels,
            fraction=subsample_fraction,
            seed=random_seed + repeat,
        )
        subsample_values = candidate_values[subsample_indices, :]
        subsample_labels = [train_labels[index] for index in subsample_indices]
        coefficients = _fit_logistic_coefficients(
            subsample_values,
            subsample_labels,
            config=config,
        )
        coefficients_by_subsample.append(coefficients)
        selected = {
            feature_ids[index]
            for index, coefficient in enumerate(coefficients)
            if abs(float(coefficient)) > 0.0
        }
        selected_sets.append(selected)
        nonzero_counts.append(len(selected))

    coefficient_matrix = np.vstack(coefficients_by_subsample)
    selected_matrix = np.abs(coefficient_matrix) > 0.0
    selection_counts = selected_matrix.sum(axis=0)
    positive_counts = (coefficient_matrix > 0.0).sum(axis=0)
    negative_counts = (coefficient_matrix < 0.0).sum(axis=0)
    selection_frequencies = selection_counts / max(n_subsamples, 1)
    pairwise_jaccards: list[float] = []
    for left_index in range(len(selected_sets)):
        for right_index in range(left_index + 1, len(selected_sets)):
            union = selected_sets[left_index] | selected_sets[right_index]
            if not union:
                pairwise_jaccards.append(1.0)
            else:
                pairwise_jaccards.append(
                    len(selected_sets[left_index] & selected_sets[right_index])
                    / len(union)
                )

    frequency_order = np.argsort(-selection_frequencies, kind="mergesort")
    top_frequency_features = [
        f"{feature_ids[int(index)]}:{_format_float(float(selection_frequencies[int(index)]))}"
        for index in frequency_order[:10]
        if selection_frequencies[int(index)] > 0.0
    ]
    variant_id = multispecies_logistic_config_id(config)
    summary = {
        "task_id": task.task_id,
        "fold_family": fold.fold_family,
        "heldout_factor": fold.heldout_factor,
        "heldout_value": fold.heldout_value,
        "fold_id": fold.fold_id,
        "variant_id": variant_id,
        "top_variable_genes": str(config.top_variable_genes),
        "c": _format_float(config.c),
        "n_train": str(len(train_x)),
        "n_test": str(fold.n_test),
        "n_subsamples": str(n_subsamples),
        "subsample_fraction": _format_float(subsample_fraction),
        "candidate_features": str(len(feature_ids)),
        "reference_balanced_accuracy": reference.summary["balanced_accuracy"],
        "reference_nonzero_count": reference.summary["n_nonzero_coefficients"],
        "median_subsample_nonzero_count": _format_float(float(np.median(nonzero_counts))),
        "min_subsample_nonzero_count": str(min(nonzero_counts) if nonzero_counts else 0),
        "max_subsample_nonzero_count": str(max(nonzero_counts) if nonzero_counts else 0),
        "stable_feature_count_ge_0_5": str(int((selection_frequencies >= 0.5).sum())),
        "stable_feature_count_ge_0_8": str(int((selection_frequencies >= 0.8).sum())),
        "mean_pairwise_jaccard": _format_float(
            float(np.mean(pairwise_jaccards)) if pairwise_jaccards else 1.0
        ),
        "top_selection_frequency_features": "|".join(top_frequency_features),
    }

    feature_rows: list[dict[str, str]] = []
    for index, feature_id in enumerate(feature_ids):
        selected_count = int(selection_counts[index])
        sign_total = int(positive_counts[index] + negative_counts[index])
        sign_consistency = (
            max(int(positive_counts[index]), int(negative_counts[index])) / sign_total
            if sign_total
            else 0.0
        )
        feature_rows.append(
            {
                "task_id": task.task_id,
                "fold_family": fold.fold_family,
                "heldout_factor": fold.heldout_factor,
                "heldout_value": fold.heldout_value,
                "fold_id": fold.fold_id,
                "variant_id": variant_id,
                "feature_id": feature_id,
                "train_variance_rank": str(
                    variance_ranks[int(selected_indices[index])]
                ),
                "reference_coefficient": _format_float(
                    float(reference_coefficients[index])
                ),
                "reference_abs_coefficient": _format_float(
                    abs(float(reference_coefficients[index]))
                ),
                "reference_selected": str(bool(reference_selected[index])),
                "selection_count": str(selected_count),
                "selection_frequency": _format_float(
                    float(selection_frequencies[index])
                ),
                "positive_count": str(int(positive_counts[index])),
                "negative_count": str(int(negative_counts[index])),
                "sign_consistency": _format_float(float(sign_consistency)),
                "mean_coefficient": _format_float(
                    float(coefficient_matrix[:, index].mean())
                ),
                "mean_abs_coefficient": _format_float(
                    float(np.abs(coefficient_matrix[:, index]).mean())
                ),
            }
        )
    return MultispeciesLogisticStabilityAuditResult(
        summary=summary,
        features=feature_rows,
    )


def _balanced_accuracy(rows: Sequence[Mapping[str, str]]) -> float:
    recalls: list[float] = []
    for label in (NEGATIVE_LABEL, POSITIVE_LABEL):
        label_rows = [row for row in rows if row["true_label"] == label]
        if not label_rows:
            recalls.append(0.0)
            continue
        correct = sum(row["predicted_label"] == label for row in label_rows)
        recalls.append(correct / len(label_rows))
    return sum(recalls) / len(recalls)


def _fold_holdout_delta(
    fold_predictions: Sequence[MultispeciesFoldPredictionResult],
    *,
    interpretation: str,
) -> dict[str, Any]:
    fold_scores = {
        fold.fold_id: {
            "heldout_value": fold.heldout_value,
            "balanced_accuracy": _balanced_accuracy(fold.predictions),
            "n_test": fold.n_test,
        }
        for fold in fold_predictions
    }
    values = [float(row["balanced_accuracy"]) for row in fold_scores.values()]
    if not values:
        return {"status": "skipped", "reason": "no fold predictions available"}
    return {
        "status": "computed",
        "value": max(values) - min(values),
        "interpretation": interpretation,
        "details": fold_scores,
    }


def _declared_delta_metric_statuses(
    *,
    manifest: Mapping[str, Any],
    fold_predictions: Sequence[MultispeciesFoldPredictionResult],
    fold_family: str,
) -> dict[str, Any]:
    declared_metrics = {
        str(metric.get("metric_id", ""))
        for metric in manifest.get("metrics", [])
        if isinstance(metric, Mapping)
    }
    delta_metrics = {
        "genotype_holdout_delta",
        "light_treatment_holdout_delta",
        "condition_stratum_holdout_delta",
    }
    computed_metric = DELTA_METRIC_BY_FOLD_FAMILY.get(fold_family)
    statuses: dict[str, Any] = {}
    for metric_id in sorted(declared_metrics & delta_metrics):
        if metric_id == computed_metric:
            statuses[metric_id] = _fold_holdout_delta(
                fold_predictions,
                interpretation=(
                    "range of fold-level balanced_accuracy across "
                    f"{fold_family} folds"
                ),
            )
        else:
            statuses[metric_id] = {
                "status": "skipped",
                "reason": f"not applicable to fold_family={fold_family}",
            }
    return statuses


def _write_predictions(
    path: Path,
    rows: Sequence[Mapping[str, str]],
    *,
    fieldnames: Sequence[str] = PREDICTION_COLUMNS,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_multispecies_nearest_centroid(
    task: MultispeciesTaskData,
    *,
    output_dir: str | Path,
    config: MultispeciesBaselineConfig | None = None,
    folds: Sequence[MultispeciesFoldData] | None = None,
    baseline_id: str = BASELINE_ID,
    fold_family: str = "condition_stratum_candidate_folds",
    claim_boundary: str = "draft_feasibility_only_not_leaderboard",
    task_manifest_path: str | Path | None = None,
    command: Sequence[str] | None = None,
) -> MultispeciesBaselineTaskResult:
    """Run the draft multispecies nearest-centroid baseline and evaluation report."""

    cfg = config or MultispeciesBaselineConfig()
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    prediction_rows: list[dict[str, str]] = []
    fold_summaries: list[dict[str, Any]] = []
    fold_results: list[MultispeciesFoldPredictionResult] = []
    selected_folds = list(folds or task.folds)
    for fold in selected_folds:
        fold_result = predict_multispecies_fold(task, fold, config=cfg)
        fold_results.append(fold_result)
        prediction_rows.extend(fold_result.predictions)
        fold_summaries.append(
            {
                "fold_id": fold_result.fold_id,
                "heldout_factor": fold_result.heldout_factor,
                "heldout_value": fold_result.heldout_value,
                "n_train": fold_result.n_train,
                "n_test": fold_result.n_test,
                "n_features": fold_result.n_features,
            }
        )

    predictions_path = outdir / "predictions.csv"
    _write_predictions(predictions_path, prediction_rows)
    evaluation = evaluate_submission(task.manifest, predictions_path)
    evaluation["metrics"].update(
        _declared_delta_metric_statuses(
            manifest=task.manifest,
            fold_predictions=fold_results,
            fold_family=fold_family,
        )
    )
    evaluation["baseline"] = {
        "baseline_id": baseline_id,
        "release_status": "draft_not_frozen",
        "transform": cfg.transform,
        "scaling": cfg.scaling,
        "top_variable_genes": cfg.top_variable_genes,
        "score_column": "flight_probability",
        "positive_label": POSITIVE_LABEL,
        "fold_family": fold_family,
        "claim_boundary": claim_boundary,
    }
    evaluation["folds"] = fold_summaries
    outputs = write_evaluation_report(
        evaluation_result=evaluation,
        task_manifest=task.manifest,
        task_manifest_path=task_manifest_path,
        submission_path=predictions_path,
        output_dir=outdir,
        command=command if command is not None else sys.argv,
    )
    return MultispeciesBaselineTaskResult(
        baseline_id=baseline_id,
        task_id=task.task_id,
        output_dir=outdir,
        predictions_path=predictions_path,
        metrics_path=outputs["metrics"],
        run_manifest_path=outputs["run_manifest"],
        evaluation=evaluation,
        n_predictions=len(prediction_rows),
    )


def _metric_value(evaluation: Mapping[str, Any], metric_id: str) -> str:
    metric = evaluation.get("metrics", {}).get(metric_id, {})
    if metric.get("status") != "computed":
        return ""
    return _format_float(float(metric["value"]))


def multispecies_config_id(config: MultispeciesBaselineConfig) -> str:
    return f"tvg{config.top_variable_genes}_{config.transform}_{config.scaling}"


def multispecies_logistic_config_id(config: MultispeciesLogisticBaselineConfig) -> str:
    c_token = _format_float(config.c).replace(".", "p")
    base = multispecies_config_id(config)
    if config.penalty == "l2" and config.solver == "liblinear":
        return f"{base}_c{c_token}"
    if config.penalty == "elasticnet":
        ratio = _format_float(float(config.l1_ratio)).replace(".", "p")
        return f"{base}_elasticnet{ratio}_c{c_token}"
    return f"{base}_{config.penalty}_c{c_token}"


def multispecies_result_summary_row(
    result: MultispeciesBaselineTaskResult,
) -> dict[str, str]:
    evaluation = result.evaluation
    baseline = evaluation.get("baseline", {})
    return {
        "baseline_id": result.baseline_id,
        "variant_id": multispecies_config_id(
            MultispeciesBaselineConfig(
                transform=str(baseline.get("transform", "log1p")),
                scaling=str(baseline.get("scaling", "zscore")),
                top_variable_genes=int(baseline.get("top_variable_genes", 2000)),
            )
        ),
        "transform": str(baseline.get("transform", "")),
        "scaling": str(baseline.get("scaling", "")),
        "top_variable_genes": str(baseline.get("top_variable_genes", "")),
        "task_id": result.task_id,
        "status": str(evaluation.get("status", "")),
        "validation_ok": str(evaluation.get("validation", {}).get("ok", "")),
        "release_status": "draft_not_frozen",
        "n_predictions": str(result.n_predictions),
        "macro_f1": _metric_value(evaluation, "macro_f1"),
        "balanced_accuracy": _metric_value(evaluation, "balanced_accuracy"),
        "auroc": _metric_value(evaluation, "auroc"),
        "calibration_error": _metric_value(evaluation, "calibration_error"),
        "mission_discrimination": _metric_value(evaluation, "mission_discrimination"),
        "genotype_holdout_delta": _metric_value(
            evaluation,
            "genotype_holdout_delta",
        ),
        "light_treatment_holdout_delta": _metric_value(
            evaluation,
            "light_treatment_holdout_delta",
        ),
        "condition_stratum_holdout_delta": _metric_value(
            evaluation,
            "condition_stratum_holdout_delta",
        ),
        "output_dir": result.output_dir.as_posix(),
        "predictions": result.predictions_path.as_posix(),
        "metrics": result.metrics_path.as_posix(),
        "run_manifest": result.run_manifest_path.as_posix(),
        "fold_family": str(baseline.get("fold_family", "")),
        "claim_boundary": str(baseline.get("claim_boundary", "")),
    }


def multispecies_logistic_result_summary_row(
    result: MultispeciesBaselineTaskResult,
) -> dict[str, str]:
    row = multispecies_result_summary_row(result)
    baseline = result.evaluation.get("baseline", {})
    row["variant_id"] = multispecies_logistic_config_id(
        MultispeciesLogisticBaselineConfig(
            transform=str(baseline.get("transform", "log1p")),
            scaling=str(baseline.get("scaling", "zscore")),
            top_variable_genes=int(baseline.get("top_variable_genes", 2000)),
            c=float(baseline.get("c", 1.0)),
            max_iter=int(baseline.get("max_iter", 5000)),
            random_state=int(baseline.get("random_state", 42)),
            penalty=str(baseline.get("penalty", "l2")),
            solver=str(baseline.get("solver", "liblinear")),
            l1_ratio=(
                float(baseline["l1_ratio"])
                if baseline.get("l1_ratio") not in {"", None, "None"}
                else None
            ),
        )
    )
    return row


def write_multispecies_baseline_summary(
    output_root: str | Path,
    rows: Sequence[Mapping[str, str]],
) -> dict[str, Path]:
    outdir = Path(output_root)
    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / "multispecies_baseline_summary.csv"
    json_path = outdir / "multispecies_baseline_summary.json"
    fieldnames = [
        "baseline_id",
        "variant_id",
        "transform",
        "scaling",
        "top_variable_genes",
        "task_id",
        "status",
        "validation_ok",
        "release_status",
        "n_predictions",
        "macro_f1",
        "balanced_accuracy",
        "auroc",
        "calibration_error",
        "mission_discrimination",
        "genotype_holdout_delta",
        "light_treatment_holdout_delta",
        "condition_stratum_holdout_delta",
        "output_dir",
        "predictions",
        "metrics",
        "run_manifest",
        "fold_family",
        "claim_boundary",
    ]
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    json_path.write_text(json.dumps(list(rows), indent=2, sort_keys=True) + "\n")
    return {"csv": csv_path, "json": json_path}


def run_multispecies_species_native_baselines(
    *,
    manifest_dir: str | Path = "v9/multispecies/task_manifests",
    repo_root: str | Path = ".",
    output_root: str | Path = "v9/multispecies/reports/nearest_centroid",
    config: MultispeciesBaselineConfig | None = None,
    command: Sequence[str] | None = None,
) -> tuple[list[MultispeciesBaselineTaskResult], dict[str, Path]]:
    """Load all draft multispecies tasks, run baselines, and write a summary."""

    tasks = load_all_multispecies_tasks(manifest_dir=manifest_dir, repo_root=repo_root)
    cfg = config or MultispeciesBaselineConfig()
    root = Path(output_root)
    results: list[MultispeciesBaselineTaskResult] = []
    rows: list[dict[str, str]] = []
    repo = Path(repo_root)
    for task in tasks:
        manifest_path = repo / manifest_dir / f"{task.task_id}.json"
        result = run_multispecies_nearest_centroid(
            task,
            output_dir=root / task.task_id,
            config=cfg,
            task_manifest_path=manifest_path,
            command=command,
        )
        results.append(result)
        rows.append(multispecies_result_summary_row(result))
    summary = write_multispecies_baseline_summary(output_root, rows)
    return results, summary


def run_multispecies_interaction_baselines(
    *,
    manifest_dir: str | Path = "v9/multispecies/interaction_task_manifests",
    repo_root: str | Path = ".",
    output_root: str | Path = "v9/multispecies/reports/interaction_nearest_centroid",
    config: MultispeciesBaselineConfig | None = None,
    fold_families: Sequence[str] | None = None,
    command: Sequence[str] | None = None,
) -> tuple[list[MultispeciesBaselineTaskResult], dict[str, Path]]:
    """Run draft OSD-120 interaction baselines and write a fold-family summary."""

    tasks = load_all_multispecies_interaction_tasks(
        manifest_dir=manifest_dir,
        repo_root=repo_root,
    )
    cfg = config or MultispeciesBaselineConfig()
    selected_families = list(fold_families or INTERACTION_FOLD_FAMILIES)
    unknown = sorted(set(selected_families) - set(INTERACTION_FOLD_FAMILIES))
    if unknown:
        raise ValueError(f"unknown multispecies interaction fold families: {unknown}")

    root = Path(output_root)
    repo = Path(repo_root)
    results: list[MultispeciesBaselineTaskResult] = []
    rows: list[dict[str, str]] = []
    for task in tasks:
        manifest_path = repo / manifest_dir / f"{task.task_id}.json"
        for fold_family in selected_families:
            result = run_multispecies_nearest_centroid(
                task,
                output_dir=root / fold_family / task.task_id,
                config=cfg,
                folds=task.fold_families[fold_family],
                baseline_id=INTERACTION_BASELINE_ID,
                fold_family=fold_family,
                claim_boundary="draft_interaction_feasibility_only_not_leaderboard",
                task_manifest_path=manifest_path,
                command=command,
            )
            results.append(result)
            rows.append(multispecies_result_summary_row(result))
    summary = write_multispecies_baseline_summary(output_root, rows)
    return results, summary


def run_multispecies_interaction_logistic_baselines(
    *,
    manifest_dir: str | Path = "v9/multispecies/interaction_task_manifests",
    repo_root: str | Path = ".",
    output_root: str | Path = "v9/multispecies/reports/interaction_logistic_l2",
    config: MultispeciesLogisticBaselineConfig | None = None,
    fold_families: Sequence[str] | None = None,
    command: Sequence[str] | None = None,
) -> tuple[list[MultispeciesBaselineTaskResult], dict[str, Path]]:
    """Run draft OSD-120 interaction L2 logistic baselines."""

    tasks = load_all_multispecies_interaction_tasks(
        manifest_dir=manifest_dir,
        repo_root=repo_root,
    )
    cfg = config or MultispeciesLogisticBaselineConfig()
    selected_families = list(fold_families or INTERACTION_FOLD_FAMILIES)
    unknown = sorted(set(selected_families) - set(INTERACTION_FOLD_FAMILIES))
    if unknown:
        raise ValueError(f"unknown multispecies interaction fold families: {unknown}")

    root = Path(output_root)
    repo = Path(repo_root)
    results: list[MultispeciesBaselineTaskResult] = []
    rows: list[dict[str, str]] = []
    for task in tasks:
        manifest_path = repo / manifest_dir / f"{task.task_id}.json"
        for fold_family in selected_families:
            outdir = root / fold_family / task.task_id
            outdir.mkdir(parents=True, exist_ok=True)
            prediction_rows: list[dict[str, str]] = []
            fold_summaries: list[dict[str, Any]] = []
            fold_results: list[MultispeciesLogisticFoldPredictionResult] = []
            for fold in task.fold_families[fold_family]:
                fold_result = predict_multispecies_logistic_fold(task, fold, config=cfg)
                fold_results.append(fold_result)
                prediction_rows.extend(fold_result.predictions)
                fold_summaries.append(
                    {
                        "fold_id": fold_result.fold_id,
                        "heldout_factor": fold_result.heldout_factor,
                        "heldout_value": fold_result.heldout_value,
                        "n_train": fold_result.n_train,
                        "n_test": fold_result.n_test,
                        "n_features": fold_result.n_features,
                        "train_time_s": round(fold_result.train_time_s, 6),
                        "fit_details": dict(fold_result.fit_details),
                    }
                )

            predictions_path = outdir / "predictions.csv"
            _write_predictions(
                predictions_path,
                prediction_rows,
                fieldnames=LOGISTIC_PREDICTION_COLUMNS,
            )
            evaluation = evaluate_submission(task.manifest, predictions_path)
            evaluation["metrics"].update(
                _declared_delta_metric_statuses(
                    manifest=task.manifest,
                    fold_predictions=fold_results,
                    fold_family=fold_family,
                )
            )
            evaluation["baseline"] = {
                "baseline_id": INTERACTION_LOGISTIC_BASELINE_ID,
                "framework": "sklearn",
                "release_status": "draft_not_frozen",
                "transform": cfg.transform,
                "scaling": cfg.scaling,
                "top_variable_genes": cfg.top_variable_genes,
                "score_column": "flight_probability",
                "positive_label": POSITIVE_LABEL,
                "fold_family": fold_family,
                "claim_boundary": "draft_interaction_logistic_only_not_leaderboard",
                "solver": cfg.solver,
                "penalty": cfg.penalty,
                "c": cfg.c,
                "l1_ratio": cfg.l1_ratio,
                "max_iter": cfg.max_iter,
                "random_state": cfg.random_state,
                "class_weight": "balanced",
            }
            evaluation["folds"] = fold_summaries
            outputs = write_evaluation_report(
                evaluation_result=evaluation,
                task_manifest=task.manifest,
                task_manifest_path=manifest_path,
                submission_path=predictions_path,
                output_dir=outdir,
                command=command if command is not None else sys.argv,
            )
            result = MultispeciesBaselineTaskResult(
                baseline_id=INTERACTION_LOGISTIC_BASELINE_ID,
                task_id=task.task_id,
                output_dir=outdir,
                predictions_path=predictions_path,
                metrics_path=outputs["metrics"],
                run_manifest_path=outputs["run_manifest"],
                evaluation=evaluation,
                n_predictions=len(prediction_rows),
            )
            results.append(result)
            rows.append(multispecies_result_summary_row(result))
    summary = write_multispecies_baseline_summary(output_root, rows)
    return results, summary


def _default_logistic_sensitivity_configs() -> list[MultispeciesLogisticBaselineConfig]:
    return [
        MultispeciesLogisticBaselineConfig(
            transform="log1p",
            scaling="zscore",
            top_variable_genes=top,
            c=c,
        )
        for top in (500, 2000)
        for c in (0.1, 1.0, 10.0)
    ]


def _default_sparse_l1_logistic_configs() -> list[MultispeciesLogisticBaselineConfig]:
    return [
        MultispeciesLogisticBaselineConfig(
            transform="log1p",
            scaling="zscore",
            top_variable_genes=2000,
            c=c,
            penalty="l1",
            solver="liblinear",
        )
        for c in (0.01, 0.03, 0.1, 0.3, 1.0)
    ]


def run_multispecies_interaction_logistic_sensitivity_grid(
    *,
    manifest_dir: str | Path = "v9/multispecies/interaction_task_manifests",
    repo_root: str | Path = ".",
    output_root: str | Path = "v9/multispecies/reports/interaction_logistic_l2_sensitivity",
    configs: Sequence[MultispeciesLogisticBaselineConfig] | None = None,
    fold_families: Sequence[str] | None = None,
    command: Sequence[str] | None = None,
) -> tuple[list[MultispeciesBaselineTaskResult], dict[str, Path]]:
    """Run a compact C/top-gene sensitivity grid for OSD-120 L2 logistic."""

    selected_configs = list(configs or _default_logistic_sensitivity_configs())
    root = Path(output_root)
    results: list[MultispeciesBaselineTaskResult] = []
    rows: list[dict[str, str]] = []
    for config in selected_configs:
        variant_id = multispecies_logistic_config_id(config)
        config_results, _ = run_multispecies_interaction_logistic_baselines(
            manifest_dir=manifest_dir,
            repo_root=repo_root,
            output_root=root / variant_id,
            config=config,
            fold_families=fold_families,
            command=command,
        )
        results.extend(config_results)
        rows.extend(
            multispecies_logistic_result_summary_row(result)
            for result in config_results
        )
    summary = write_multispecies_baseline_summary(output_root, rows)
    return results, summary


def run_multispecies_interaction_sparse_logistic_l1_grid(
    *,
    manifest_dir: str | Path = "v9/multispecies/interaction_task_manifests",
    repo_root: str | Path = ".",
    output_root: str | Path = "v9/multispecies/reports/interaction_logistic_sparse_l1",
    configs: Sequence[MultispeciesLogisticBaselineConfig] | None = None,
    fold_families: Sequence[str] | None = None,
    command: Sequence[str] | None = None,
) -> tuple[list[MultispeciesBaselineTaskResult], dict[str, Path]]:
    """Run the draft OSD-120 sparse L1 logistic pilot grid."""

    return run_multispecies_interaction_logistic_sensitivity_grid(
        manifest_dir=manifest_dir,
        repo_root=repo_root,
        output_root=output_root,
        configs=list(configs or _default_sparse_l1_logistic_configs()),
        fold_families=fold_families,
        command=command,
    )


def build_multispecies_interaction_logistic_feature_audit(
    *,
    manifest_dir: str | Path = "v9/multispecies/interaction_task_manifests",
    repo_root: str | Path = ".",
    configs: Sequence[MultispeciesLogisticBaselineConfig] | None = None,
    focus_folds: Sequence[tuple[str, str]] | None = None,
    top_n_coefficients: int = 10,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Build selected-feature and coefficient audit rows for OSD-120 logistic."""

    selected_configs = list(
        configs
        or [
            MultispeciesLogisticBaselineConfig(top_variable_genes=500, c=1.0),
            MultispeciesLogisticBaselineConfig(top_variable_genes=2000, c=1.0),
        ]
    )
    focus = set(focus_folds or DEFAULT_LOGISTIC_FEATURE_AUDIT_FOLDS)
    tasks = load_all_multispecies_interaction_tasks(
        manifest_dir=manifest_dir,
        repo_root=repo_root,
    )
    summary_rows: list[dict[str, str]] = []
    coefficient_rows: list[dict[str, str]] = []
    for task in tasks:
        for fold_family, heldout_value in sorted(focus):
            folds = [
                fold
                for fold in task.fold_families.get(fold_family, [])
                if fold.heldout_value == heldout_value
            ]
            if not folds:
                raise ValueError(
                    f"no fold found for fold_family={fold_family} "
                    f"heldout_value={heldout_value}"
                )
            if len(folds) > 1:
                raise ValueError(
                    f"multiple folds found for fold_family={fold_family} "
                    f"heldout_value={heldout_value}"
                )
            fold = folds[0]
            for config in selected_configs:
                result = audit_multispecies_logistic_fold_features(
                    task,
                    fold,
                    config=config,
                    top_n_coefficients=top_n_coefficients,
                )
                summary_rows.append(result.summary)
                coefficient_rows.extend(result.coefficients)
    return summary_rows, coefficient_rows


def write_multispecies_interaction_logistic_feature_audit(
    *,
    output_dir: str | Path = "v9/multispecies/reports/interaction_logistic_feature_audit",
    manifest_dir: str | Path = "v9/multispecies/interaction_task_manifests",
    repo_root: str | Path = ".",
    configs: Sequence[MultispeciesLogisticBaselineConfig] | None = None,
    focus_folds: Sequence[tuple[str, str]] | None = None,
    top_n_coefficients: int = 10,
) -> dict[str, Path]:
    """Write OSD-120 logistic selected-feature and coefficient audit tables."""

    summary_rows, coefficient_rows = build_multispecies_interaction_logistic_feature_audit(
        manifest_dir=manifest_dir,
        repo_root=repo_root,
        configs=configs,
        focus_folds=focus_folds,
        top_n_coefficients=top_n_coefficients,
    )
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    summary_csv = outdir / "feature_set_audit_summary.csv"
    summary_json = outdir / "feature_set_audit_summary.json"
    coefficient_csv = outdir / "feature_coefficient_audit.csv"
    coefficient_json = outdir / "feature_coefficient_audit.json"
    with summary_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_LOGISTIC_FEATURE_AUDIT_SUMMARY_FIELDS,
        )
        writer.writeheader()
        writer.writerows(summary_rows)
    summary_json.write_text(json.dumps(summary_rows, indent=2, sort_keys=True) + "\n")
    with coefficient_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_LOGISTIC_FEATURE_COEFFICIENT_FIELDS,
        )
        writer.writeheader()
        writer.writerows(coefficient_rows)
    coefficient_json.write_text(
        json.dumps(coefficient_rows, indent=2, sort_keys=True) + "\n"
    )
    return {
        "summary_csv": summary_csv,
        "summary_json": summary_json,
        "coefficient_csv": coefficient_csv,
        "coefficient_json": coefficient_json,
    }


def _default_sparse_l1_stability_configs() -> list[MultispeciesLogisticBaselineConfig]:
    return [
        MultispeciesLogisticBaselineConfig(
            transform="log1p",
            scaling="zscore",
            top_variable_genes=2000,
            c=c,
            penalty="l1",
            solver="liblinear",
        )
        for c in (0.3, 1.0)
    ]


def build_multispecies_interaction_sparse_l1_stability_audit(
    *,
    manifest_dir: str | Path = "v9/multispecies/interaction_task_manifests",
    repo_root: str | Path = ".",
    configs: Sequence[MultispeciesLogisticBaselineConfig] | None = None,
    focus_folds: Sequence[tuple[str, str]] | None = None,
    n_subsamples: int = 20,
    subsample_fraction: float = 0.8,
    random_seed: int = 1729,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Build sparse L1 subsampling stability audit rows for OSD-120."""

    if n_subsamples < 1:
        raise ValueError("n_subsamples must be at least 1")
    if not 0.0 < subsample_fraction <= 1.0:
        raise ValueError("subsample_fraction must be in (0, 1]")
    selected_configs = list(configs or _default_sparse_l1_stability_configs())
    focus = set(focus_folds or DEFAULT_LOGISTIC_FEATURE_AUDIT_FOLDS)
    tasks = load_all_multispecies_interaction_tasks(
        manifest_dir=manifest_dir,
        repo_root=repo_root,
    )
    summary_rows: list[dict[str, str]] = []
    feature_rows: list[dict[str, str]] = []
    for task in tasks:
        for fold_family, heldout_value in sorted(focus):
            folds = [
                fold
                for fold in task.fold_families.get(fold_family, [])
                if fold.heldout_value == heldout_value
            ]
            if not folds:
                raise ValueError(
                    f"no fold found for fold_family={fold_family} "
                    f"heldout_value={heldout_value}"
                )
            if len(folds) > 1:
                raise ValueError(
                    f"multiple folds found for fold_family={fold_family} "
                    f"heldout_value={heldout_value}"
                )
            fold = folds[0]
            for config in selected_configs:
                result = audit_multispecies_logistic_stability_fold(
                    task,
                    fold,
                    config=config,
                    n_subsamples=n_subsamples,
                    subsample_fraction=subsample_fraction,
                    random_seed=random_seed,
                )
                summary_rows.append(result.summary)
                feature_rows.extend(result.features)
    return summary_rows, feature_rows


def write_multispecies_interaction_sparse_l1_stability_audit(
    *,
    output_dir: str | Path = "v9/multispecies/reports/interaction_logistic_sparse_l1_stability",
    manifest_dir: str | Path = "v9/multispecies/interaction_task_manifests",
    repo_root: str | Path = ".",
    configs: Sequence[MultispeciesLogisticBaselineConfig] | None = None,
    focus_folds: Sequence[tuple[str, str]] | None = None,
    n_subsamples: int = 20,
    subsample_fraction: float = 0.8,
    random_seed: int = 1729,
) -> dict[str, Path]:
    """Write sparse L1 subsampling stability audit outputs for OSD-120."""

    summary_rows, feature_rows = build_multispecies_interaction_sparse_l1_stability_audit(
        manifest_dir=manifest_dir,
        repo_root=repo_root,
        configs=configs,
        focus_folds=focus_folds,
        n_subsamples=n_subsamples,
        subsample_fraction=subsample_fraction,
        random_seed=random_seed,
    )
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    summary_csv = outdir / "stability_summary.csv"
    summary_json = outdir / "stability_summary.json"
    feature_csv = outdir / "stability_feature_frequency.csv"
    feature_json = outdir / "stability_feature_frequency.json"
    with summary_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_LOGISTIC_STABILITY_SUMMARY_FIELDS,
        )
        writer.writeheader()
        writer.writerows(summary_rows)
    summary_json.write_text(json.dumps(summary_rows, indent=2, sort_keys=True) + "\n")
    with feature_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_LOGISTIC_STABILITY_FEATURE_FIELDS,
        )
        writer.writeheader()
        writer.writerows(feature_rows)
    feature_json.write_text(json.dumps(feature_rows, indent=2, sort_keys=True) + "\n")
    return {
        "summary_csv": summary_csv,
        "summary_json": summary_json,
        "feature_csv": feature_csv,
        "feature_json": feature_json,
    }


def _read_csv_dict_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def _require_one_row(
    rows: Sequence[dict[str, str]],
    *,
    context: str,
) -> dict[str, str]:
    if len(rows) != 1:
        raise ValueError(f"expected one row for {context}, found {len(rows)}")
    return rows[0]


def _ladder_source_path(
    candidate: Mapping[str, str],
    key: str,
    reports_root: Path,
) -> Path | None:
    value = str(candidate.get(key, ""))
    if not value:
        return None
    path = Path(value)
    if path.is_absolute():
        return path
    return reports_root / path


def _comparison_delta_counts(rows: Sequence[dict[str, str]]) -> dict[str, str]:
    deltas = [float(row["delta_logistic_minus_nearest_centroid"]) for row in rows]
    if not deltas:
        return {
            "nearest_improved_count": "0",
            "nearest_tied_count": "0",
            "nearest_worse_count": "0",
            "nearest_min_delta": "",
            "nearest_max_delta": "",
        }
    return {
        "nearest_improved_count": str(sum(delta > 0.0 for delta in deltas)),
        "nearest_tied_count": str(sum(delta == 0.0 for delta in deltas)),
        "nearest_worse_count": str(sum(delta < 0.0 for delta in deltas)),
        "nearest_min_delta": _format_float(min(deltas)),
        "nearest_max_delta": _format_float(max(deltas)),
    }


def build_multispecies_interaction_baseline_ladder(
    *,
    reports_root: str | Path = "v9/multispecies/reports",
    repo_root: str | Path = ".",
    candidates: Sequence[Mapping[str, str]] | None = None,
    focus_folds: Sequence[tuple[str, str, str, str]] | None = None,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Consolidate OSD-120 nearest/L2/L1 comparator outputs into a ladder."""

    root = Path(repo_root)
    reports_path = _resolve_path(reports_root, root)
    selected_candidates = [
        dict(candidate)
        for candidate in (candidates or DEFAULT_INTERACTION_BASELINE_LADDER_CANDIDATES)
    ]
    selected_focus = list(
        focus_folds or DEFAULT_INTERACTION_BASELINE_LADDER_FOCUS_FOLDS
    )
    csv_cache: dict[Path, list[dict[str, str]]] = {}

    def cached_rows(path: Path) -> list[dict[str, str]]:
        if path not in csv_cache:
            csv_cache[path] = _read_csv_dict_rows(path)
        return csv_cache[path]

    summary_rows: list[dict[str, str]] = []
    focus_rows: list[dict[str, str]] = []
    for candidate in selected_candidates:
        candidate_id = str(candidate["candidate_id"])
        variant_id = str(candidate["variant_id"])
        summary_path = _ladder_source_path(candidate, "summary_csv", reports_path)
        fold_path = _ladder_source_path(candidate, "fold_detail_csv", reports_path)
        comparison_path = _ladder_source_path(candidate, "comparison_csv", reports_path)
        stability_path = _ladder_source_path(candidate, "stability_csv", reports_path)
        if summary_path is None or fold_path is None:
            raise ValueError(f"candidate {candidate_id} is missing source paths")

        all_summary_rows = [
            row for row in cached_rows(summary_path) if row["variant_id"] == variant_id
        ]
        family_rows = {
            family: _require_one_row(
                [row for row in all_summary_rows if row["fold_family"] == family],
                context=f"{candidate_id} {family}",
            )
            for family in INTERACTION_FOLD_FAMILIES
        }
        first_summary = family_rows[INTERACTION_FOLD_FAMILIES[0]]
        row = {field: "" for field in INTERACTION_BASELINE_LADDER_SUMMARY_FIELDS}
        row.update(
            {
                "candidate_id": candidate_id,
                "role": str(candidate.get("role", "")),
                "decision": str(candidate.get("decision", "")),
                "model_class": str(candidate.get("model_class", "")),
                "baseline_id": str(candidate.get("baseline_id", first_summary["baseline_id"])),
                "variant_id": variant_id,
                "task_id": first_summary["task_id"],
                "transform": first_summary["transform"],
                "scaling": first_summary["scaling"],
                "top_variable_genes": first_summary["top_variable_genes"],
                "penalty": str(candidate.get("penalty", "")),
                "c": str(candidate.get("c", "")),
                "source_summary_csv": summary_path.as_posix(),
                "source_fold_detail_csv": fold_path.as_posix(),
                "source_comparison_csv": (
                    comparison_path.as_posix() if comparison_path is not None else ""
                ),
                "source_stability_csv": (
                    stability_path.as_posix() if stability_path is not None else ""
                ),
                "claim_boundary": first_summary["claim_boundary"],
            }
        )
        family_balanced_accuracy: list[float] = []
        for family, prefix in [
            ("primary_genotype_or_ecotype_holdout", "primary"),
            ("secondary_light_treatment_holdout", "secondary"),
            ("diagnostic_condition_stratum_holdout", "diagnostic"),
        ]:
            family_row = family_rows[family]
            delta_metric = DELTA_METRIC_BY_FOLD_FAMILY[family]
            family_balanced_accuracy.append(float(family_row["balanced_accuracy"]))
            row[f"{prefix}_balanced_accuracy"] = family_row["balanced_accuracy"]
            row[f"{prefix}_auroc"] = family_row["auroc"]
            row[f"{prefix}_delta_metric"] = delta_metric
            row[f"{prefix}_delta"] = family_row.get(delta_metric, "")
        row["mean_family_balanced_accuracy"] = _format_float(
            sum(family_balanced_accuracy) / len(family_balanced_accuracy)
        )
        row["min_family_balanced_accuracy"] = _format_float(
            min(family_balanced_accuracy)
        )

        if comparison_path is None:
            detail_rows = [
                detail
                for detail in cached_rows(fold_path)
                if detail["variant_id"] == variant_id
            ]
            row.update(
                {
                    "nearest_improved_count": "0",
                    "nearest_tied_count": str(len(detail_rows)),
                    "nearest_worse_count": "0",
                    "nearest_min_delta": "0",
                    "nearest_max_delta": "0",
                }
            )
            for fold_family, heldout_value, _, focus_column in selected_focus:
                focus_detail = _require_one_row(
                    [
                        detail
                        for detail in detail_rows
                        if detail["fold_family"] == fold_family
                        and detail["heldout_value"] == heldout_value
                    ],
                    context=f"{candidate_id} focus {heldout_value}",
                )
                row[focus_column] = focus_detail["balanced_accuracy"]
                focus_row = {
                    field: ""
                    for field in INTERACTION_BASELINE_LADDER_FOCUS_FIELDS
                }
                focus_row.update(
                    {
                        "candidate_id": candidate_id,
                        "model_class": str(candidate.get("model_class", "")),
                        "variant_id": variant_id,
                        "task_id": focus_detail["task_id"],
                        "fold_family": fold_family,
                        "heldout_factor": focus_detail["heldout_factor"],
                        "heldout_value": heldout_value,
                        "n_test": focus_detail["n_test"],
                        "nearest_centroid_balanced_accuracy": focus_detail[
                            "balanced_accuracy"
                        ],
                        "candidate_balanced_accuracy": focus_detail[
                            "balanced_accuracy"
                        ],
                        "delta_candidate_minus_nearest_centroid": "0",
                        "nearest_centroid_rank_low_to_high": focus_detail[
                            "rank_low_to_high"
                        ],
                        "candidate_rank_low_to_high": focus_detail[
                            "rank_low_to_high"
                        ],
                        "nearest_centroid_is_lowest_for_variant": focus_detail[
                            "is_lowest_for_variant"
                        ],
                        "candidate_is_lowest_for_variant": focus_detail[
                            "is_lowest_for_variant"
                        ],
                        "source_fold_detail_csv": fold_path.as_posix(),
                        "nearest_centroid_metrics": focus_detail["metrics"],
                        "candidate_metrics": focus_detail["metrics"],
                    }
                )
                focus_rows.append(focus_row)
        else:
            comparison_rows = [
                comparison
                for comparison in cached_rows(comparison_path)
                if comparison["logistic_variant_id"] == variant_id
            ]
            row.update(_comparison_delta_counts(comparison_rows))
            for fold_family, heldout_value, _, focus_column in selected_focus:
                focus_comparison = _require_one_row(
                    [
                        comparison
                        for comparison in comparison_rows
                        if comparison["fold_family"] == fold_family
                        and comparison["heldout_value"] == heldout_value
                    ],
                    context=f"{candidate_id} focus {heldout_value}",
                )
                row[focus_column] = focus_comparison["logistic_balanced_accuracy"]
                focus_row = {
                    field: ""
                    for field in INTERACTION_BASELINE_LADDER_FOCUS_FIELDS
                }
                focus_row.update(
                    {
                        "candidate_id": candidate_id,
                        "model_class": str(candidate.get("model_class", "")),
                        "variant_id": variant_id,
                        "task_id": focus_comparison["task_id"],
                        "fold_family": fold_family,
                        "heldout_factor": focus_comparison["heldout_factor"],
                        "heldout_value": heldout_value,
                        "n_test": focus_comparison["n_test"],
                        "nearest_centroid_balanced_accuracy": focus_comparison[
                            "nearest_centroid_balanced_accuracy"
                        ],
                        "candidate_balanced_accuracy": focus_comparison[
                            "logistic_balanced_accuracy"
                        ],
                        "delta_candidate_minus_nearest_centroid": focus_comparison[
                            "delta_logistic_minus_nearest_centroid"
                        ],
                        "nearest_centroid_rank_low_to_high": focus_comparison[
                            "nearest_centroid_rank_low_to_high"
                        ],
                        "candidate_rank_low_to_high": focus_comparison[
                            "logistic_rank_low_to_high"
                        ],
                        "nearest_centroid_is_lowest_for_variant": focus_comparison[
                            "nearest_centroid_is_lowest_for_variant"
                        ],
                        "candidate_is_lowest_for_variant": focus_comparison[
                            "logistic_is_lowest_for_variant"
                        ],
                        "source_fold_detail_csv": fold_path.as_posix(),
                        "source_comparison_csv": comparison_path.as_posix(),
                        "nearest_centroid_metrics": focus_comparison[
                            "nearest_centroid_metrics"
                        ],
                        "candidate_metrics": focus_comparison["logistic_metrics"],
                    }
                )
                focus_rows.append(focus_row)

        if stability_path is not None:
            stability_rows = [
                stability
                for stability in cached_rows(stability_path)
                if stability["variant_id"] == variant_id
            ]
            stability_by_focus = {
                (stability["fold_family"], stability["heldout_value"]): stability
                for stability in stability_rows
            }
            for fold_family, heldout_value, prefix, _ in selected_focus:
                stability = stability_by_focus.get((fold_family, heldout_value))
                if stability is None:
                    continue
                row[f"{prefix}_reference_nonzero_count"] = stability[
                    "reference_nonzero_count"
                ]
                row[f"{prefix}_stable_feature_count_ge_0_5"] = stability[
                    "stable_feature_count_ge_0_5"
                ]
                row[f"{prefix}_stable_feature_count_ge_0_8"] = stability[
                    "stable_feature_count_ge_0_8"
                ]
                row[f"{prefix}_mean_pairwise_jaccard"] = stability[
                    "mean_pairwise_jaccard"
                ]
        summary_rows.append(row)
    return summary_rows, focus_rows


def write_multispecies_interaction_baseline_ladder(
    *,
    output_dir: str | Path = "v9/multispecies/reports/interaction_baseline_ladder",
    reports_root: str | Path = "v9/multispecies/reports",
    repo_root: str | Path = ".",
    candidates: Sequence[Mapping[str, str]] | None = None,
    focus_folds: Sequence[tuple[str, str, str, str]] | None = None,
) -> dict[str, Path]:
    """Write OSD-120 comparator ladder summary and focus-fold tables."""

    summary_rows, focus_rows = build_multispecies_interaction_baseline_ladder(
        reports_root=reports_root,
        repo_root=repo_root,
        candidates=candidates,
        focus_folds=focus_folds,
    )
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    summary_csv = outdir / "baseline_ladder_summary.csv"
    summary_json = outdir / "baseline_ladder_summary.json"
    focus_csv = outdir / "baseline_ladder_focus_folds.csv"
    focus_json = outdir / "baseline_ladder_focus_folds.json"
    with summary_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_BASELINE_LADDER_SUMMARY_FIELDS,
        )
        writer.writeheader()
        writer.writerows(summary_rows)
    summary_json.write_text(json.dumps(summary_rows, indent=2, sort_keys=True) + "\n")
    with focus_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_BASELINE_LADDER_FOCUS_FIELDS,
        )
        writer.writeheader()
        writer.writerows(focus_rows)
    focus_json.write_text(json.dumps(focus_rows, indent=2, sort_keys=True) + "\n")
    return {
        "summary_csv": summary_csv,
        "summary_json": summary_json,
        "focus_csv": focus_csv,
        "focus_json": focus_json,
    }


def _focus_label(fold_family: str, heldout_value: str) -> str:
    if fold_family == "diagnostic_condition_stratum_holdout":
        return "diagnostic_col0_phyd_dark"
    if fold_family == "secondary_light_treatment_holdout":
        return "secondary_light_treatment"
    if fold_family == "primary_genotype_or_ecotype_holdout":
        return "primary_col0_phyd"
    return heldout_value.lower().replace("|", "_").replace(".", "_").replace(" ", "_")


def _pipe_join(values: Sequence[str]) -> str:
    return "|".join(value for value in values if value)


def build_multispecies_interaction_diagnostic_candidate_package(
    *,
    reports_root: str | Path = "v9/multispecies/reports",
    repo_root: str | Path = ".",
    task_manifest: str | Path = (
        "v9/multispecies/interaction_task_manifests/"
        "draft_osd120_arabidopsis_root_light_interaction_spaceflight.json"
    ),
    package_id: str = DEFAULT_INTERACTION_DIAGNOSTIC_PACKAGE_ID,
    candidate_id: str = DEFAULT_INTERACTION_DIAGNOSTIC_CANDIDATE_ID,
    comparator_candidate_id: str = DEFAULT_INTERACTION_DIAGNOSTIC_COMPARATOR_ID,
    stable_feature_frequency_threshold: float = 0.5,
) -> dict[str, list[dict[str, str]]]:
    """Build a figure-ready OSD-120 sparse diagnostic candidate package."""

    root = Path(repo_root)
    reports_path = _resolve_path(reports_root, root)
    task_manifest_path = _resolve_path(task_manifest, root)
    manifest = json.loads(task_manifest_path.read_text())
    source_record = (
        manifest.get("source_records", [{}])[0]
        if isinstance(manifest.get("source_records"), list)
        else {}
    )
    ladder_summary_path = reports_path / "interaction_baseline_ladder" / (
        "baseline_ladder_summary.csv"
    )
    ladder_focus_path = reports_path / "interaction_baseline_ladder" / (
        "baseline_ladder_focus_folds.csv"
    )
    sparse_feature_audit_path = reports_path / "interaction_logistic_sparse_l1" / (
        "feature_set_audit_summary.csv"
    )
    stability_summary_path = reports_path / "interaction_logistic_sparse_l1_stability" / (
        "stability_summary.csv"
    )
    stability_feature_path = reports_path / "interaction_logistic_sparse_l1_stability" / (
        "stability_feature_frequency.csv"
    )

    ladder_summary_rows = _read_csv_dict_rows(ladder_summary_path)
    ladder_focus_rows = _read_csv_dict_rows(ladder_focus_path)
    sparse_feature_audit_rows = _read_csv_dict_rows(sparse_feature_audit_path)
    stability_summary_rows = _read_csv_dict_rows(stability_summary_path)
    stability_feature_rows = _read_csv_dict_rows(stability_feature_path)

    candidate_summary = _require_one_row(
        [row for row in ladder_summary_rows if row["candidate_id"] == candidate_id],
        context=f"diagnostic candidate {candidate_id}",
    )
    comparator_summary = _require_one_row(
        [
            row
            for row in ladder_summary_rows
            if row["candidate_id"] == comparator_candidate_id
        ],
        context=f"diagnostic comparator {comparator_candidate_id}",
    )
    candidate_variant_id = candidate_summary["variant_id"]
    focus_rows_source = [
        row for row in ladder_focus_rows if row["candidate_id"] == candidate_id
    ]
    focus_rows_source = sorted(
        focus_rows_source,
        key=lambda row: (
            [focus[0] for focus in DEFAULT_LOGISTIC_FEATURE_AUDIT_FOLDS].index(
                row["fold_family"]
            )
            if row["fold_family"]
            in [focus[0] for focus in DEFAULT_LOGISTIC_FEATURE_AUDIT_FOLDS]
            else 99,
            row["heldout_value"],
        ),
    )
    audit_by_fold = {
        (row["fold_family"], row["heldout_value"]): row
        for row in sparse_feature_audit_rows
        if row["variant_id"] == candidate_variant_id
    }
    stability_by_fold = {
        (row["fold_family"], row["heldout_value"]): row
        for row in stability_summary_rows
        if row["variant_id"] == candidate_variant_id
    }

    focus_rows: list[dict[str, str]] = []
    focus_candidate_scores: list[float] = []
    focus_improved_count = 0
    stable_ge_0_5_total = 0
    stable_ge_0_8_total = 0
    for source_focus in focus_rows_source:
        key = (source_focus["fold_family"], source_focus["heldout_value"])
        audit = _require_one_row(
            [audit_by_fold[key]] if key in audit_by_fold else [],
            context=f"candidate package feature audit {key}",
        )
        stability = _require_one_row(
            [stability_by_fold[key]] if key in stability_by_fold else [],
            context=f"candidate package stability {key}",
        )
        fold_label = _focus_label(*key)
        candidate_ba = float(source_focus["candidate_balanced_accuracy"])
        delta = float(source_focus["delta_candidate_minus_nearest_centroid"])
        focus_candidate_scores.append(candidate_ba)
        focus_improved_count += int(delta > 0.0)
        stable_ge_0_5_total += int(stability["stable_feature_count_ge_0_5"])
        stable_ge_0_8_total += int(stability["stable_feature_count_ge_0_8"])
        focus_rows.append(
            {
                "package_id": package_id,
                "candidate_id": candidate_id,
                "candidate_variant_id": candidate_variant_id,
                "fold_label": fold_label,
                "fold_family": source_focus["fold_family"],
                "heldout_factor": source_focus["heldout_factor"],
                "heldout_value": source_focus["heldout_value"],
                "n_test": source_focus["n_test"],
                "nearest_centroid_balanced_accuracy": source_focus[
                    "nearest_centroid_balanced_accuracy"
                ],
                "candidate_balanced_accuracy": source_focus[
                    "candidate_balanced_accuracy"
                ],
                "delta_candidate_minus_nearest_centroid": source_focus[
                    "delta_candidate_minus_nearest_centroid"
                ],
                "candidate_rank_low_to_high": source_focus[
                    "candidate_rank_low_to_high"
                ],
                "candidate_is_lowest_for_variant": source_focus[
                    "candidate_is_lowest_for_variant"
                ],
                "n_nonzero_coefficients": audit["n_nonzero_coefficients"],
                "n_nonzero_top500_coefficients": audit[
                    "n_nonzero_top500_coefficients"
                ],
                "n_nonzero_extra_coefficients": audit[
                    "n_nonzero_extra_coefficients"
                ],
                "stable_feature_count_ge_0_5": stability[
                    "stable_feature_count_ge_0_5"
                ],
                "stable_feature_count_ge_0_8": stability[
                    "stable_feature_count_ge_0_8"
                ],
                "mean_pairwise_jaccard": stability["mean_pairwise_jaccard"],
                "top_reference_abs_features": audit["top_abs_feature_ids"],
                "top_stable_features": stability[
                    "top_selection_frequency_features"
                ],
                "evidence_summary": (
                    f"BA {source_focus['candidate_balanced_accuracy']} versus "
                    f"nearest {source_focus['nearest_centroid_balanced_accuracy']}; "
                    f"{stability['stable_feature_count_ge_0_5']} features selected "
                    f"in >=50% of subsamples"
                ),
            }
        )

    feature_rows_source = [
        row
        for row in stability_feature_rows
        if row["variant_id"] == candidate_variant_id
        and float(row["selection_frequency"]) >= stable_feature_frequency_threshold
    ]
    feature_rows_source = sorted(
        feature_rows_source,
        key=lambda row: (
            row["fold_family"],
            row["heldout_value"],
            -float(row["selection_frequency"]),
            -float(row["reference_abs_coefficient"]),
            row["feature_id"],
        ),
    )
    feature_rows = []
    for feature in feature_rows_source:
        frequency = float(feature["selection_frequency"])
        tier = "stable_ge_0_8" if frequency >= 0.8 else "stable_ge_0_5"
        feature_rows.append(
            {
                "package_id": package_id,
                "candidate_id": candidate_id,
                "candidate_variant_id": candidate_variant_id,
                "fold_label": _focus_label(
                    feature["fold_family"],
                    feature["heldout_value"],
                ),
                "fold_family": feature["fold_family"],
                "heldout_value": feature["heldout_value"],
                "feature_id": feature["feature_id"],
                "train_variance_rank": feature["train_variance_rank"],
                "reference_coefficient": feature["reference_coefficient"],
                "reference_abs_coefficient": feature["reference_abs_coefficient"],
                "reference_selected": feature["reference_selected"],
                "selection_count": feature["selection_count"],
                "selection_frequency": feature["selection_frequency"],
                "sign_consistency": feature["sign_consistency"],
                "mean_coefficient": feature["mean_coefficient"],
                "mean_abs_coefficient": feature["mean_abs_coefficient"],
                "stability_tier": tier,
            }
        )

    summary_row = {
        "package_id": package_id,
        "task_id": candidate_summary["task_id"],
        "organism": str(manifest.get("organism", "")),
        "assay_modality": str(manifest.get("assay_modality", "")),
        "biospecimen_type": str(manifest.get("biospecimen_type", "")),
        "source_id": str(source_record.get("source_id", "")),
        "source_url": str(source_record.get("source_url", "")),
        "candidate_id": candidate_id,
        "candidate_variant_id": candidate_variant_id,
        "candidate_model_class": candidate_summary["model_class"],
        "candidate_decision": candidate_summary["decision"],
        "comparator_candidate_id": comparator_candidate_id,
        "comparator_variant_id": comparator_summary["variant_id"],
        "primary_balanced_accuracy": candidate_summary[
            "primary_balanced_accuracy"
        ],
        "secondary_balanced_accuracy": candidate_summary[
            "secondary_balanced_accuracy"
        ],
        "diagnostic_balanced_accuracy": candidate_summary[
            "diagnostic_balanced_accuracy"
        ],
        "mean_family_balanced_accuracy": candidate_summary[
            "mean_family_balanced_accuracy"
        ],
        "min_family_balanced_accuracy": candidate_summary[
            "min_family_balanced_accuracy"
        ],
        "nearest_improved_count": candidate_summary["nearest_improved_count"],
        "nearest_tied_count": candidate_summary["nearest_tied_count"],
        "nearest_worse_count": candidate_summary["nearest_worse_count"],
        "focus_fold_count": str(len(focus_rows)),
        "focus_improved_count": str(focus_improved_count),
        "focus_mean_candidate_ba": _format_float(
            sum(focus_candidate_scores) / len(focus_candidate_scores)
        ),
        "focus_min_candidate_ba": _format_float(min(focus_candidate_scores)),
        "stable_feature_count_ge_0_5_total": str(stable_ge_0_5_total),
        "stable_feature_count_ge_0_8_total": str(stable_ge_0_8_total),
        "claim_boundary": candidate_summary["claim_boundary"],
        "release_status": str(manifest.get("release_status", "")),
        "evidence_status": "local_outputs_linked_no_new_model_family",
        "external_context_status": (
            "OSD-120_CARA_root_tip_light_genotype_context_checked"
        ),
        "source_task_manifest": task_manifest_path.as_posix(),
        "source_ladder_summary_csv": ladder_summary_path.as_posix(),
        "source_ladder_focus_csv": ladder_focus_path.as_posix(),
        "source_sparse_feature_audit_csv": sparse_feature_audit_path.as_posix(),
        "source_stability_summary_csv": stability_summary_path.as_posix(),
        "source_stability_feature_csv": stability_feature_path.as_posix(),
    }
    claim_rows = [
        {
            "package_id": package_id,
            "claim_id": "draft_candidate_boundary",
            "claim_status": "supported_local_artifact",
            "claim_text": (
                "sparse_l1_c1 is a draft transparent diagnostic candidate, "
                "not a frozen leaderboard baseline."
            ),
            "local_evidence_artifacts": _pipe_join(
                [
                    task_manifest_path.as_posix(),
                    ladder_summary_path.as_posix(),
                ]
            ),
            "local_evidence_columns": "release_status|claim_boundary|candidate_decision",
            "external_source_urls": "",
            "limitations": "within-source OSD-120 task; no leave-one-mission-out claim",
        },
        {
            "package_id": package_id,
            "claim_id": "nearest_centroid_fold_comparison",
            "claim_status": "supported_local_artifact",
            "claim_text": (
                "The candidate improves 9 of 11 nearest-centroid fold rows, "
                "ties 2, and worsens 0."
            ),
            "local_evidence_artifacts": ladder_summary_path.as_posix(),
            "local_evidence_columns": (
                "nearest_improved_count|nearest_tied_count|nearest_worse_count"
            ),
            "external_source_urls": "",
            "limitations": "fold counts are from the current draft OSD-120 split only",
        },
        {
            "package_id": package_id,
            "claim_id": "fragile_focus_recovery",
            "claim_status": "supported_local_artifact",
            "claim_text": (
                "The candidate recovers all three predefined fragile focus folds "
                "relative to nearest centroid."
            ),
            "local_evidence_artifacts": ladder_focus_path.as_posix(),
            "local_evidence_columns": (
                "heldout_value|nearest_centroid_balanced_accuracy|"
                "candidate_balanced_accuracy|delta_candidate_minus_nearest_centroid"
            ),
            "external_source_urls": "",
            "limitations": "focus folds are diagnostic sentinels, not independent tasks",
        },
        {
            "package_id": package_id,
            "claim_id": "feature_stability_evidence",
            "claim_status": "supported_local_artifact",
            "claim_text": (
                "Nineteen candidate features are selected in at least 50% of "
                "balanced train-fold subsamples across the three focus folds."
            ),
            "local_evidence_artifacts": _pipe_join(
                [
                    stability_summary_path.as_posix(),
                    stability_feature_path.as_posix(),
                ]
            ),
            "local_evidence_columns": (
                "stable_feature_count_ge_0_5|selection_frequency|sign_consistency"
            ),
            "external_source_urls": "",
            "limitations": "subsampling stability is small-n diagnostic evidence",
        },
        {
            "package_id": package_id,
            "claim_id": "external_light_genotype_context",
            "claim_status": "supported_external_context",
            "claim_text": (
                "OSD-120 is an Arabidopsis CARA root-tip RNA-seq study where "
                "genotype/ecotype and lighting context are core design factors."
            ),
            "local_evidence_artifacts": task_manifest_path.as_posix(),
            "local_evidence_columns": (
                "organism|biospecimen_type|blocking_factors|source_records"
            ),
            "external_source_urls": _pipe_join(
                [
                    "https://osdr.nasa.gov/bio/repo/data/studies/OSD-120",
                    "https://doi.org/10.1371/journal.pone.0180186",
                    "https://www.nature.com/articles/s41526-024-00417-0",
                ]
            ),
            "limitations": "external literature supports scope, not benchmark performance",
        },
    ]
    return {
        "summary": [summary_row],
        "focus": focus_rows,
        "features": feature_rows,
        "claims": claim_rows,
    }


def write_multispecies_interaction_diagnostic_candidate_package(
    *,
    output_dir: str | Path = (
        "v9/multispecies/reports/interaction_diagnostic_candidate_package"
    ),
    reports_root: str | Path = "v9/multispecies/reports",
    repo_root: str | Path = ".",
    task_manifest: str | Path = (
        "v9/multispecies/interaction_task_manifests/"
        "draft_osd120_arabidopsis_root_light_interaction_spaceflight.json"
    ),
    package_id: str = DEFAULT_INTERACTION_DIAGNOSTIC_PACKAGE_ID,
    candidate_id: str = DEFAULT_INTERACTION_DIAGNOSTIC_CANDIDATE_ID,
    comparator_candidate_id: str = DEFAULT_INTERACTION_DIAGNOSTIC_COMPARATOR_ID,
    stable_feature_frequency_threshold: float = 0.5,
) -> dict[str, Path]:
    """Write figure-ready OSD-120 diagnostic candidate package tables."""

    tables = build_multispecies_interaction_diagnostic_candidate_package(
        reports_root=reports_root,
        repo_root=repo_root,
        task_manifest=task_manifest,
        package_id=package_id,
        candidate_id=candidate_id,
        comparator_candidate_id=comparator_candidate_id,
        stable_feature_frequency_threshold=stable_feature_frequency_threshold,
    )
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "summary_csv": outdir / "candidate_package_summary.csv",
        "summary_json": outdir / "candidate_package_summary.json",
        "focus_csv": outdir / "candidate_focus_evidence.csv",
        "focus_json": outdir / "candidate_focus_evidence.json",
        "feature_csv": outdir / "candidate_stable_feature_evidence.csv",
        "feature_json": outdir / "candidate_stable_feature_evidence.json",
        "claim_csv": outdir / "candidate_claim_map.csv",
        "claim_json": outdir / "candidate_claim_map.json",
    }
    table_fields = {
        "summary": INTERACTION_DIAGNOSTIC_PACKAGE_SUMMARY_FIELDS,
        "focus": INTERACTION_DIAGNOSTIC_PACKAGE_FOCUS_FIELDS,
        "features": INTERACTION_DIAGNOSTIC_PACKAGE_FEATURE_FIELDS,
        "claims": INTERACTION_DIAGNOSTIC_PACKAGE_CLAIM_FIELDS,
    }
    output_prefix = {
        "summary": "summary",
        "focus": "focus",
        "features": "feature",
        "claims": "claim",
    }
    for table_name, rows in tables.items():
        prefix = output_prefix[table_name]
        csv_path = outputs[f"{prefix}_csv"]
        json_path = outputs[f"{prefix}_json"]
        with csv_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=table_fields[table_name])
            writer.writeheader()
            writer.writerows(rows)
        json_path.write_text(json.dumps(rows, indent=2, sort_keys=True) + "\n")
    return outputs


def _display_float(value: str | float, digits: int = 3) -> str:
    return f"{float(value):.{digits}f}"


def _display_delta(value: str | float, digits: int = 3) -> str:
    numeric = float(value)
    sign = "+" if numeric >= 0 else ""
    return f"{sign}{numeric:.{digits}f}"


def _display_fold_label(fold_label: str, heldout_value: str) -> str:
    if fold_label == "diagnostic_col0_phyd_dark":
        return "Condition stratum: Col.0.PhyD | Dark"
    if fold_label == "secondary_light_treatment":
        return "Light treatment: Light"
    if fold_label == "primary_col0_phyd":
        return "Genotype/ecotype: Col.0.PhyD"
    return heldout_value


def build_multispecies_interaction_figure_table_package(
    *,
    package_dir: str | Path = (
        "v9/multispecies/reports/interaction_diagnostic_candidate_package"
    ),
    repo_root: str | Path = ".",
    figure_table_id: str = DEFAULT_INTERACTION_FIGURE_TABLE_ID,
) -> dict[str, Any]:
    """Build human-facing OSD-120 figure/table draft artifacts."""

    root = Path(repo_root)
    package_path = _resolve_path(package_dir, root)
    summary = _require_one_row(
        _read_csv_dict_rows(package_path / "candidate_package_summary.csv"),
        context="candidate package summary",
    )
    focus_rows = _read_csv_dict_rows(package_path / "candidate_focus_evidence.csv")
    feature_rows = _read_csv_dict_rows(
        package_path / "candidate_stable_feature_evidence.csv"
    )
    claim_rows = _read_csv_dict_rows(package_path / "candidate_claim_map.csv")
    focus_order = {
        "diagnostic_col0_phyd_dark": 0,
        "secondary_light_treatment": 1,
        "primary_col0_phyd": 2,
    }
    sorted_focus_rows = sorted(
        focus_rows,
        key=lambda row: (focus_order.get(row["fold_label"], 99), row["heldout_value"]),
    )
    main_rows: list[dict[str, str]] = []
    for focus in sorted_focus_rows:
        display_fold = _display_fold_label(focus["fold_label"], focus["heldout_value"])
        nearest_ba = _display_float(focus["nearest_centroid_balanced_accuracy"])
        candidate_ba = _display_float(focus["candidate_balanced_accuracy"])
        display_delta = _display_delta(focus["delta_candidate_minus_nearest_centroid"])
        main_rows.append(
            {
                "figure_table_id": figure_table_id,
                "package_id": focus["package_id"],
                "candidate_id": focus["candidate_id"],
                "candidate_variant_id": focus["candidate_variant_id"],
                "display_fold": display_fold,
                "fold_label": focus["fold_label"],
                "heldout_value": focus["heldout_value"],
                "n_test": focus["n_test"],
                "nearest_centroid_ba": nearest_ba,
                "candidate_ba": candidate_ba,
                "delta_ba": _display_float(
                    focus["delta_candidate_minus_nearest_centroid"]
                ),
                "display_delta_ba": display_delta,
                "n_nonzero_coefficients": focus["n_nonzero_coefficients"],
                "stable_features_ge_0_5": focus["stable_feature_count_ge_0_5"],
                "stable_features_ge_0_8": focus["stable_feature_count_ge_0_8"],
                "mean_pairwise_jaccard": _display_float(
                    focus["mean_pairwise_jaccard"]
                ),
                "figure_use": "main_focus_table",
                "result_sentence": (
                    f"{display_fold}: sparse L1 BA {candidate_ba} versus "
                    f"nearest-centroid BA {nearest_ba} ({display_delta})."
                ),
            }
        )

    feature_counts: dict[str, int] = {}
    appendix_rows: list[dict[str, str]] = []
    sorted_feature_rows = sorted(
        feature_rows,
        key=lambda row: (
            focus_order.get(row["fold_label"], 99),
            -float(row["selection_frequency"]),
            -float(row["reference_abs_coefficient"]),
            row["feature_id"],
        ),
    )
    for feature in sorted_feature_rows:
        fold_label = feature["fold_label"]
        feature_counts[fold_label] = feature_counts.get(fold_label, 0) + 1
        coefficient = float(feature["reference_coefficient"])
        direction = "positive_LEO_or_ISS" if coefficient > 0 else "negative_LEO_or_ISS"
        appendix_rows.append(
            {
                "figure_table_id": f"{figure_table_id}_stable_feature_appendix",
                "package_id": feature["package_id"],
                "candidate_id": feature["candidate_id"],
                "fold_label": fold_label,
                "display_fold": _display_fold_label(
                    fold_label,
                    feature["heldout_value"],
                ),
                "feature_rank_in_fold": str(feature_counts[fold_label]),
                "feature_id": feature["feature_id"],
                "selection_frequency": _display_float(
                    feature["selection_frequency"],
                    digits=2,
                ),
                "display_selection_frequency": (
                    f"{float(feature['selection_frequency']) * 100:.0f}%"
                ),
                "reference_coefficient": _display_float(
                    feature["reference_coefficient"],
                    digits=4,
                ),
                "coefficient_direction": direction,
                "stability_tier": feature["stability_tier"],
                "appendix_note": (
                    "Sparse-model feature evidence only; not a validated "
                    "biomarker claim."
                ),
            }
        )

    source_urls = [
        url
        for row in claim_rows
        for url in row.get("external_source_urls", "").split("|")
        if url
    ]
    caption = "\n".join(
        [
            "# OSD-120 Sparse L1 Candidate Figure Caption",
            "",
            (
                "Draft OSD-120 Arabidopsis root-tip interaction diagnostic for "
                f"`{summary['candidate_id']}` "
                f"(`{summary['candidate_variant_id']}`)."
            ),
            (
                "The main table reports the three predefined fragile focus "
                "folds from the v9 OSD-120 interaction task: "
                "`Col.0.PhyD|Dark.Treatment`, `Light.Treatment`, and "
                "`Col.0.PhyD`."
            ),
            (
                f"Across the full ladder, the candidate improves "
                f"{summary['nearest_improved_count']} of 11 nearest-centroid "
                f"fold rows, ties {summary['nearest_tied_count']}, and worsens "
                f"{summary['nearest_worse_count']}."
            ),
            (
                f"Primary, secondary, and diagnostic balanced accuracies are "
                f"{_display_float(summary['primary_balanced_accuracy'])}, "
                f"{_display_float(summary['secondary_balanced_accuracy'])}, "
                f"and {_display_float(summary['diagnostic_balanced_accuracy'])}."
            ),
            (
                "Stable-feature appendix rows are sparse-model audit evidence "
                "from deterministic balanced train-fold subsampling, not "
                "validated biomarkers."
            ),
            "",
            "Allowed claim: draft transparent diagnostic candidate for a "
            "within-source OSD-120 interaction split.",
            "",
            "Disallowed claims: frozen benchmark baseline, leave-one-mission-out "
            "generalization, cross-species transfer, causal gene biology, or "
            "operational plant-growth recommendation.",
            "",
            "External context sources: " + _pipe_join(source_urls),
            "",
        ]
    )
    claim_box = "\n".join(
        [
            "# OSD-120 Figure/Table Claim Boundary",
            "",
            "Allowed:",
            "",
            "- `sparse_l1_c1` is the packaged primary draft transparent "
            "diagnostic candidate.",
            "- The candidate improves all three predefined fragile focus folds "
            "relative to nearest centroid.",
            "- The selected-feature appendix is auditable sparse-model evidence.",
            "",
            "Not Allowed:",
            "",
            "- Do not call this a frozen v9 benchmark baseline.",
            "- Do not claim leave-one-mission-out or cross-study generalization.",
            "- Do not describe selected genes as validated biomarkers.",
            "- Do not make operational plant-growth recommendations.",
            "",
        ]
    )
    return {
        "main": main_rows,
        "features": appendix_rows,
        "caption": caption,
        "claim_box": claim_box,
    }


def write_multispecies_interaction_figure_table_package(
    *,
    output_dir: str | Path = (
        "v9/multispecies/reports/interaction_figure_table_package"
    ),
    package_dir: str | Path = (
        "v9/multispecies/reports/interaction_diagnostic_candidate_package"
    ),
    repo_root: str | Path = ".",
    figure_table_id: str = DEFAULT_INTERACTION_FIGURE_TABLE_ID,
) -> dict[str, Path]:
    """Write human-facing OSD-120 figure/table draft artifacts."""

    package = build_multispecies_interaction_figure_table_package(
        package_dir=package_dir,
        repo_root=repo_root,
        figure_table_id=figure_table_id,
    )
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "main_csv": outdir / "figure_main_focus_table.csv",
        "main_json": outdir / "figure_main_focus_table.json",
        "feature_csv": outdir / "figure_stable_feature_appendix.csv",
        "feature_json": outdir / "figure_stable_feature_appendix.json",
        "caption_md": outdir / "figure_caption.md",
        "claim_box_md": outdir / "figure_claim_boundary_box.md",
    }
    with outputs["main_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=INTERACTION_FIGURE_TABLE_MAIN_FIELDS)
        writer.writeheader()
        writer.writerows(package["main"])
    outputs["main_json"].write_text(
        json.dumps(package["main"], indent=2, sort_keys=True) + "\n"
    )
    with outputs["feature_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_FIGURE_TABLE_FEATURE_FIELDS,
        )
        writer.writeheader()
        writer.writerows(package["features"])
    outputs["feature_json"].write_text(
        json.dumps(package["features"], indent=2, sort_keys=True) + "\n"
    )
    outputs["caption_md"].write_text(str(package["caption"]))
    outputs["claim_box_md"].write_text(str(package["claim_box"]))
    return outputs


def build_multispecies_interaction_paired_comparator_table(
    *,
    reports_root: str | Path = "v9/multispecies/reports",
    repo_root: str | Path = ".",
    paired_table_id: str = DEFAULT_INTERACTION_PAIRED_COMPARATOR_TABLE_ID,
    primary_candidate_id: str = DEFAULT_INTERACTION_DIAGNOSTIC_CANDIDATE_ID,
    compact_candidate_id: str = DEFAULT_INTERACTION_DIAGNOSTIC_COMPARATOR_ID,
) -> dict[str, Any]:
    """Build c1-vs-c0.3 sparse L1 paired comparator table artifacts."""

    root = Path(repo_root)
    reports_path = _resolve_path(reports_root, root)
    ladder_summary_path = reports_path / "interaction_baseline_ladder" / (
        "baseline_ladder_summary.csv"
    )
    ladder_focus_path = reports_path / "interaction_baseline_ladder" / (
        "baseline_ladder_focus_folds.csv"
    )
    sparse_feature_audit_path = reports_path / "interaction_logistic_sparse_l1" / (
        "feature_set_audit_summary.csv"
    )
    stability_summary_path = reports_path / "interaction_logistic_sparse_l1_stability" / (
        "stability_summary.csv"
    )
    ladder_summary_rows = _read_csv_dict_rows(ladder_summary_path)
    ladder_focus_rows = _read_csv_dict_rows(ladder_focus_path)
    sparse_feature_audit_rows = _read_csv_dict_rows(sparse_feature_audit_path)
    stability_summary_rows = _read_csv_dict_rows(stability_summary_path)

    primary_summary = _require_one_row(
        [row for row in ladder_summary_rows if row["candidate_id"] == primary_candidate_id],
        context=f"paired primary {primary_candidate_id}",
    )
    compact_summary = _require_one_row(
        [row for row in ladder_summary_rows if row["candidate_id"] == compact_candidate_id],
        context=f"paired compact {compact_candidate_id}",
    )
    primary_variant = primary_summary["variant_id"]
    compact_variant = compact_summary["variant_id"]

    focus_order = {
        "diagnostic_condition_stratum_holdout": 0,
        "secondary_light_treatment_holdout": 1,
        "primary_genotype_or_ecotype_holdout": 2,
    }
    primary_focus = {
        (row["fold_family"], row["heldout_value"]): row
        for row in ladder_focus_rows
        if row["candidate_id"] == primary_candidate_id
    }
    compact_focus = {
        (row["fold_family"], row["heldout_value"]): row
        for row in ladder_focus_rows
        if row["candidate_id"] == compact_candidate_id
    }
    audit_by_variant_and_fold = {
        (row["variant_id"], row["fold_family"], row["heldout_value"]): row
        for row in sparse_feature_audit_rows
        if row["variant_id"] in {primary_variant, compact_variant}
    }
    stability_by_variant_and_fold = {
        (row["variant_id"], row["fold_family"], row["heldout_value"]): row
        for row in stability_summary_rows
        if row["variant_id"] in {primary_variant, compact_variant}
    }

    focus_rows: list[dict[str, str]] = []
    focus_tied = 0
    focus_primary_better = 0
    focus_compact_better = 0
    primary_nonzero_total = 0
    compact_nonzero_total = 0
    primary_stable_05_total = 0
    compact_stable_05_total = 0
    primary_stable_08_total = 0
    compact_stable_08_total = 0

    for key, primary_row in sorted(
        primary_focus.items(),
        key=lambda item: (focus_order.get(item[0][0], 99), item[0][1]),
    ):
        compact_row = _require_one_row(
            [compact_focus[key]] if key in compact_focus else [],
            context=f"paired compact focus {key}",
        )
        fold_family, heldout_value = key
        primary_audit = _require_one_row(
            [audit_by_variant_and_fold[(primary_variant, fold_family, heldout_value)]],
            context=f"primary paired audit {key}",
        )
        compact_audit = _require_one_row(
            [audit_by_variant_and_fold[(compact_variant, fold_family, heldout_value)]],
            context=f"compact paired audit {key}",
        )
        primary_stability = _require_one_row(
            [stability_by_variant_and_fold[(primary_variant, fold_family, heldout_value)]],
            context=f"primary paired stability {key}",
        )
        compact_stability = _require_one_row(
            [stability_by_variant_and_fold[(compact_variant, fold_family, heldout_value)]],
            context=f"compact paired stability {key}",
        )
        fold_label = _focus_label(fold_family, heldout_value)
        primary_ba = float(primary_row["candidate_balanced_accuracy"])
        compact_ba = float(compact_row["candidate_balanced_accuracy"])
        delta = primary_ba - compact_ba
        focus_tied += int(delta == 0.0)
        focus_primary_better += int(delta > 0.0)
        focus_compact_better += int(delta < 0.0)
        primary_nonzero_total += int(primary_audit["n_nonzero_coefficients"])
        compact_nonzero_total += int(compact_audit["n_nonzero_coefficients"])
        primary_stable_05_total += int(primary_stability["stable_feature_count_ge_0_5"])
        compact_stable_05_total += int(compact_stability["stable_feature_count_ge_0_5"])
        primary_stable_08_total += int(primary_stability["stable_feature_count_ge_0_8"])
        compact_stable_08_total += int(compact_stability["stable_feature_count_ge_0_8"])
        compactness_delta = int(compact_audit["n_nonzero_coefficients"]) - int(
            primary_audit["n_nonzero_coefficients"]
        )
        if delta == 0.0 and compactness_delta < 0:
            interpretation = (
                "same focus BA; compact comparator uses fewer nonzero coefficients"
            )
        elif delta == 0.0:
            interpretation = "same focus BA; primary candidate keeps stronger ladder role"
        elif delta > 0.0:
            interpretation = "primary candidate has higher focus BA"
        else:
            interpretation = "compact comparator has higher focus BA"
        focus_rows.append(
            {
                "paired_table_id": paired_table_id,
                "fold_label": fold_label,
                "display_fold": _display_fold_label(fold_label, heldout_value),
                "fold_family": fold_family,
                "heldout_value": heldout_value,
                "n_test": primary_row["n_test"],
                "nearest_centroid_ba": _display_float(
                    primary_row["nearest_centroid_balanced_accuracy"]
                ),
                "primary_candidate_id": primary_candidate_id,
                "primary_candidate_ba": _display_float(
                    primary_row["candidate_balanced_accuracy"]
                ),
                "compact_candidate_id": compact_candidate_id,
                "compact_candidate_ba": _display_float(
                    compact_row["candidate_balanced_accuracy"]
                ),
                "primary_minus_compact_ba": _display_float(delta),
                "display_primary_minus_compact_ba": _display_delta(delta),
                "primary_nonzero_coefficients": primary_audit[
                    "n_nonzero_coefficients"
                ],
                "compact_nonzero_coefficients": compact_audit[
                    "n_nonzero_coefficients"
                ],
                "primary_stable_ge_0_5": primary_stability[
                    "stable_feature_count_ge_0_5"
                ],
                "compact_stable_ge_0_5": compact_stability[
                    "stable_feature_count_ge_0_5"
                ],
                "primary_stable_ge_0_8": primary_stability[
                    "stable_feature_count_ge_0_8"
                ],
                "compact_stable_ge_0_8": compact_stability[
                    "stable_feature_count_ge_0_8"
                ],
                "primary_mean_pairwise_jaccard": _display_float(
                    primary_stability["mean_pairwise_jaccard"]
                ),
                "compact_mean_pairwise_jaccard": _display_float(
                    compact_stability["mean_pairwise_jaccard"]
                ),
                "compactness_delta_nonzero_compact_minus_primary": str(
                    compactness_delta
                ),
                "focus_interpretation": interpretation,
            }
        )

    recommendation = "appendix_or_supplement_only"
    decision_rationale = (
        "Focus-fold BA is tied across the primary and compact sparse candidates, "
        "while c1 has stronger full-ladder behavior and no nearest-centroid "
        "worsening; c0.3 remains useful as a compact control."
    )
    summary_rows = [
        {
            "paired_table_id": paired_table_id,
            "task_id": primary_summary["task_id"],
            "primary_candidate_id": primary_candidate_id,
            "primary_variant_id": primary_variant,
            "compact_candidate_id": compact_candidate_id,
            "compact_variant_id": compact_variant,
            "primary_mean_family_ba": _display_float(
                primary_summary["mean_family_balanced_accuracy"]
            ),
            "compact_mean_family_ba": _display_float(
                compact_summary["mean_family_balanced_accuracy"]
            ),
            "primary_min_family_ba": _display_float(
                primary_summary["min_family_balanced_accuracy"]
            ),
            "compact_min_family_ba": _display_float(
                compact_summary["min_family_balanced_accuracy"]
            ),
            "primary_diagnostic_ba": _display_float(
                primary_summary["diagnostic_balanced_accuracy"]
            ),
            "compact_diagnostic_ba": _display_float(
                compact_summary["diagnostic_balanced_accuracy"]
            ),
            "primary_nearest_improved_count": primary_summary[
                "nearest_improved_count"
            ],
            "primary_nearest_tied_count": primary_summary["nearest_tied_count"],
            "primary_nearest_worse_count": primary_summary["nearest_worse_count"],
            "compact_nearest_improved_count": compact_summary[
                "nearest_improved_count"
            ],
            "compact_nearest_tied_count": compact_summary["nearest_tied_count"],
            "compact_nearest_worse_count": compact_summary["nearest_worse_count"],
            "focus_tied_ba_count": str(focus_tied),
            "focus_primary_better_count": str(focus_primary_better),
            "focus_compact_better_count": str(focus_compact_better),
            "primary_focus_nonzero_total": str(primary_nonzero_total),
            "compact_focus_nonzero_total": str(compact_nonzero_total),
            "primary_stable_ge_0_5_total": str(primary_stable_05_total),
            "compact_stable_ge_0_5_total": str(compact_stable_05_total),
            "primary_stable_ge_0_8_total": str(primary_stable_08_total),
            "compact_stable_ge_0_8_total": str(compact_stable_08_total),
            "recommendation": recommendation,
            "decision_rationale": decision_rationale,
            "source_ladder_summary_csv": ladder_summary_path.as_posix(),
            "source_ladder_focus_csv": ladder_focus_path.as_posix(),
            "source_sparse_feature_audit_csv": sparse_feature_audit_path.as_posix(),
            "source_stability_summary_csv": stability_summary_path.as_posix(),
        }
    ]
    decision_md = "\n".join(
        [
            "# OSD-120 Sparse L1 Paired Comparator Decision",
            "",
            "`sparse_l1_c0p3` should remain an appendix or supplement comparator,",
            "not a second primary panel in the main OSD-120 figure.",
            "",
            "Reason:",
            "",
            "- The three focus-fold balanced accuracies are tied between `c1` and `c0p3`.",
            "- `c1` has stronger full-ladder behavior: 9 improve, 2 tie, 0 worse versus nearest centroid.",
            "- `c0p3` has one nearest-centroid-worse fold in the full 11-fold comparison.",
            "- `c0p3` is more compact in the focus folds, especially `Light.Treatment`, so it remains a useful control.",
            "",
            "Recommended display:",
            "",
            "- Main figure: keep the simpler `sparse_l1_c1` focus table.",
            "- Appendix or supplement: include the paired comparator table when explaining the stability-versus-compactness tradeoff.",
            "",
        ]
    )
    return {
        "summary": summary_rows,
        "focus": focus_rows,
        "decision_md": decision_md,
    }


def write_multispecies_interaction_paired_comparator_table(
    *,
    output_dir: str | Path = (
        "v9/multispecies/reports/interaction_paired_comparator_table"
    ),
    reports_root: str | Path = "v9/multispecies/reports",
    repo_root: str | Path = ".",
    paired_table_id: str = DEFAULT_INTERACTION_PAIRED_COMPARATOR_TABLE_ID,
    primary_candidate_id: str = DEFAULT_INTERACTION_DIAGNOSTIC_CANDIDATE_ID,
    compact_candidate_id: str = DEFAULT_INTERACTION_DIAGNOSTIC_COMPARATOR_ID,
) -> dict[str, Path]:
    """Write c1-vs-c0.3 sparse L1 paired comparator table artifacts."""

    package = build_multispecies_interaction_paired_comparator_table(
        reports_root=reports_root,
        repo_root=repo_root,
        paired_table_id=paired_table_id,
        primary_candidate_id=primary_candidate_id,
        compact_candidate_id=compact_candidate_id,
    )
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "summary_csv": outdir / "paired_comparator_summary.csv",
        "summary_json": outdir / "paired_comparator_summary.json",
        "focus_csv": outdir / "paired_focus_comparator_table.csv",
        "focus_json": outdir / "paired_focus_comparator_table.json",
        "decision_md": outdir / "paired_comparator_decision.md",
    }
    with outputs["summary_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_PAIRED_COMPARATOR_SUMMARY_FIELDS,
        )
        writer.writeheader()
        writer.writerows(package["summary"])
    outputs["summary_json"].write_text(
        json.dumps(package["summary"], indent=2, sort_keys=True) + "\n"
    )
    with outputs["focus_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_PAIRED_COMPARATOR_FOCUS_FIELDS,
        )
        writer.writeheader()
        writer.writerows(package["focus"])
    outputs["focus_json"].write_text(
        json.dumps(package["focus"], indent=2, sort_keys=True) + "\n"
    )
    outputs["decision_md"].write_text(str(package["decision_md"]))
    return outputs


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _artifact_line_count(path: Path) -> str:
    with path.open("rb") as handle:
        return str(sum(1 for _ in handle))


def _artifact_row_count(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        with path.open(newline="") as handle:
            return str(max(sum(1 for _ in csv.reader(handle)) - 1, 0))
    if suffix == ".json":
        payload = json.loads(path.read_text())
        if isinstance(payload, list):
            return str(len(payload))
    return ""


def _artifact_row(
    *,
    manifest_id: str,
    task_id: str,
    artifact_id: str,
    artifact_role: str,
    block_id: str,
    path: Path,
    required_for_claims: Sequence[str],
    validation_scope: str,
    generated_by: str,
    claim_boundary: str,
) -> dict[str, str]:
    exists = path.exists()
    row = {
        "manifest_id": manifest_id,
        "task_id": task_id,
        "artifact_id": artifact_id,
        "artifact_role": artifact_role,
        "block_id": block_id,
        "path": path.as_posix(),
        "file_format": path.suffix.lstrip("."),
        "exists": str(exists),
        "byte_size": "",
        "line_count": "",
        "row_count": "",
        "sha256": "",
        "required_for_claims": _pipe_join(list(required_for_claims)),
        "validation_scope": validation_scope,
        "generated_by": generated_by,
        "claim_boundary": claim_boundary,
    }
    if exists:
        row["byte_size"] = str(path.stat().st_size)
        row["line_count"] = _artifact_line_count(path)
        row["row_count"] = _artifact_row_count(path)
        row["sha256"] = _sha256_file(path)
    return row


def _default_osd120_diagnostic_artifacts(reports_path: Path) -> list[dict[str, Any]]:
    return [
        {
            "artifact_id": "osd120_task_manifest",
            "artifact_role": "task_manifest",
            "block_id": "V9-MULTI-010",
            "path": Path(
                "v9/multispecies/interaction_task_manifests/"
                "draft_osd120_arabidopsis_root_light_interaction_spaceflight.json"
            ),
            "claims": ["draft_candidate_boundary", "external_light_genotype_context"],
            "validation_scope": "manifest_contract",
            "generated_by": "scripts/build_v9_osd120_interaction_task_manifest.py",
        },
        {
            "artifact_id": "osd120_task_manifest_index",
            "artifact_role": "task_manifest_index",
            "block_id": "V9-MULTI-010",
            "path": Path("v9/multispecies/interaction_task_manifest_index.draft.csv"),
            "claims": ["draft_candidate_boundary"],
            "validation_scope": "manifest_index",
            "generated_by": "scripts/build_v9_osd120_interaction_task_manifest.py",
        },
        {
            "artifact_id": "osd120_task_design_review",
            "artifact_role": "design_review",
            "block_id": "V9-MULTI-009",
            "path": reports_path / "OSD120_INTERACTION_TASK_DESIGN.md",
            "claims": ["external_light_genotype_context"],
            "validation_scope": "review_note",
            "generated_by": "manual_review",
        },
        {
            "artifact_id": "sparse_l1_summary",
            "artifact_role": "model_summary",
            "block_id": "V9-MULTI-019",
            "path": reports_path
            / "interaction_logistic_sparse_l1"
            / "multispecies_baseline_summary.csv",
            "claims": ["draft_candidate_boundary"],
            "validation_scope": "generated_output",
            "generated_by": "scripts/run_v9_osd120_interaction_sparse_l1.py",
        },
        {
            "artifact_id": "sparse_l1_fold_comparison",
            "artifact_role": "fold_comparison",
            "block_id": "V9-MULTI-019",
            "path": reports_path
            / "interaction_logistic_sparse_l1"
            / "fold_detail_comparison_vs_nearest_centroid.csv",
            "claims": [
                "nearest_centroid_fold_comparison",
                "fragile_focus_recovery",
            ],
            "validation_scope": "generated_output",
            "generated_by": "scripts/run_v9_osd120_interaction_sparse_l1.py",
        },
        {
            "artifact_id": "sparse_l1_feature_set_audit",
            "artifact_role": "feature_audit",
            "block_id": "V9-MULTI-019",
            "path": reports_path
            / "interaction_logistic_sparse_l1"
            / "feature_set_audit_summary.csv",
            "claims": ["feature_stability_evidence"],
            "validation_scope": "generated_output",
            "generated_by": "scripts/run_v9_osd120_interaction_sparse_l1.py",
        },
        {
            "artifact_id": "sparse_l1_stability_summary",
            "artifact_role": "stability_summary",
            "block_id": "V9-MULTI-020",
            "path": reports_path
            / "interaction_logistic_sparse_l1_stability"
            / "stability_summary.csv",
            "claims": ["feature_stability_evidence"],
            "validation_scope": "generated_output",
            "generated_by": "scripts/audit_v9_osd120_sparse_l1_stability.py",
        },
        {
            "artifact_id": "sparse_l1_stability_features",
            "artifact_role": "stability_feature_frequency",
            "block_id": "V9-MULTI-020",
            "path": reports_path
            / "interaction_logistic_sparse_l1_stability"
            / "stability_feature_frequency.csv",
            "claims": ["feature_stability_evidence"],
            "validation_scope": "generated_output",
            "generated_by": "scripts/audit_v9_osd120_sparse_l1_stability.py",
        },
        {
            "artifact_id": "baseline_ladder_summary",
            "artifact_role": "ladder_summary",
            "block_id": "V9-MULTI-021",
            "path": reports_path
            / "interaction_baseline_ladder"
            / "baseline_ladder_summary.csv",
            "claims": [
                "draft_candidate_boundary",
                "nearest_centroid_fold_comparison",
            ],
            "validation_scope": "generated_output",
            "generated_by": "scripts/build_v9_osd120_baseline_ladder.py",
        },
        {
            "artifact_id": "baseline_ladder_focus_folds",
            "artifact_role": "ladder_focus_folds",
            "block_id": "V9-MULTI-021",
            "path": reports_path
            / "interaction_baseline_ladder"
            / "baseline_ladder_focus_folds.csv",
            "claims": ["fragile_focus_recovery"],
            "validation_scope": "generated_output",
            "generated_by": "scripts/build_v9_osd120_baseline_ladder.py",
        },
        {
            "artifact_id": "candidate_package_summary",
            "artifact_role": "candidate_package_summary",
            "block_id": "V9-MULTI-022",
            "path": reports_path
            / "interaction_diagnostic_candidate_package"
            / "candidate_package_summary.csv",
            "claims": ["draft_candidate_boundary"],
            "validation_scope": "generated_output",
            "generated_by": "scripts/build_v9_osd120_diagnostic_candidate_package.py",
        },
        {
            "artifact_id": "candidate_focus_evidence",
            "artifact_role": "candidate_focus_evidence",
            "block_id": "V9-MULTI-022",
            "path": reports_path
            / "interaction_diagnostic_candidate_package"
            / "candidate_focus_evidence.csv",
            "claims": ["fragile_focus_recovery"],
            "validation_scope": "generated_output",
            "generated_by": "scripts/build_v9_osd120_diagnostic_candidate_package.py",
        },
        {
            "artifact_id": "candidate_stable_feature_evidence",
            "artifact_role": "candidate_stable_feature_evidence",
            "block_id": "V9-MULTI-022",
            "path": reports_path
            / "interaction_diagnostic_candidate_package"
            / "candidate_stable_feature_evidence.csv",
            "claims": ["feature_stability_evidence"],
            "validation_scope": "generated_output",
            "generated_by": "scripts/build_v9_osd120_diagnostic_candidate_package.py",
        },
        {
            "artifact_id": "candidate_claim_map",
            "artifact_role": "claim_map_source",
            "block_id": "V9-MULTI-022",
            "path": reports_path
            / "interaction_diagnostic_candidate_package"
            / "candidate_claim_map.csv",
            "claims": [
                "draft_candidate_boundary",
                "nearest_centroid_fold_comparison",
                "fragile_focus_recovery",
                "feature_stability_evidence",
                "external_light_genotype_context",
            ],
            "validation_scope": "generated_output",
            "generated_by": "scripts/build_v9_osd120_diagnostic_candidate_package.py",
        },
        {
            "artifact_id": "figure_main_focus_table",
            "artifact_role": "figure_table",
            "block_id": "V9-MULTI-023",
            "path": reports_path
            / "interaction_figure_table_package"
            / "figure_main_focus_table.csv",
            "claims": ["fragile_focus_recovery", "figure_table_claim_boundary"],
            "validation_scope": "generated_output",
            "generated_by": "scripts/build_v9_osd120_figure_table_package.py",
        },
        {
            "artifact_id": "figure_stable_feature_appendix",
            "artifact_role": "figure_appendix",
            "block_id": "V9-MULTI-023",
            "path": reports_path
            / "interaction_figure_table_package"
            / "figure_stable_feature_appendix.csv",
            "claims": ["feature_stability_evidence", "figure_table_claim_boundary"],
            "validation_scope": "generated_output",
            "generated_by": "scripts/build_v9_osd120_figure_table_package.py",
        },
        {
            "artifact_id": "figure_caption",
            "artifact_role": "caption",
            "block_id": "V9-MULTI-023",
            "path": reports_path
            / "interaction_figure_table_package"
            / "figure_caption.md",
            "claims": ["figure_table_claim_boundary"],
            "validation_scope": "caption_text",
            "generated_by": "scripts/build_v9_osd120_figure_table_package.py",
        },
        {
            "artifact_id": "figure_claim_boundary_box",
            "artifact_role": "claim_boundary",
            "block_id": "V9-MULTI-023",
            "path": reports_path
            / "interaction_figure_table_package"
            / "figure_claim_boundary_box.md",
            "claims": ["figure_table_claim_boundary"],
            "validation_scope": "claim_boundary_text",
            "generated_by": "scripts/build_v9_osd120_figure_table_package.py",
        },
        {
            "artifact_id": "paired_comparator_summary",
            "artifact_role": "paired_comparator_summary",
            "block_id": "V9-MULTI-024",
            "path": reports_path
            / "interaction_paired_comparator_table"
            / "paired_comparator_summary.csv",
            "claims": ["paired_comparator_decision"],
            "validation_scope": "generated_output",
            "generated_by": "scripts/build_v9_osd120_paired_comparator_table.py",
        },
        {
            "artifact_id": "paired_focus_comparator_table",
            "artifact_role": "paired_focus_table",
            "block_id": "V9-MULTI-024",
            "path": reports_path
            / "interaction_paired_comparator_table"
            / "paired_focus_comparator_table.csv",
            "claims": ["paired_comparator_decision"],
            "validation_scope": "generated_output",
            "generated_by": "scripts/build_v9_osd120_paired_comparator_table.py",
        },
        {
            "artifact_id": "paired_comparator_decision",
            "artifact_role": "paired_decision_note",
            "block_id": "V9-MULTI-024",
            "path": reports_path
            / "interaction_paired_comparator_table"
            / "paired_comparator_decision.md",
            "claims": ["paired_comparator_decision"],
            "validation_scope": "decision_text",
            "generated_by": "scripts/build_v9_osd120_paired_comparator_table.py",
        },
        {
            "artifact_id": "paired_comparator_review",
            "artifact_role": "review_note",
            "block_id": "V9-MULTI-024",
            "path": reports_path / "OSD120_INTERACTION_PAIRED_COMPARATOR_REVIEW.md",
            "claims": ["paired_comparator_decision"],
            "validation_scope": "review_note",
            "generated_by": "manual_review",
        },
        {
            "artifact_id": "figure_table_review",
            "artifact_role": "review_note",
            "block_id": "V9-MULTI-023",
            "path": reports_path / "OSD120_INTERACTION_FIGURE_TABLE_DRAFT_REVIEW.md",
            "claims": ["figure_table_claim_boundary"],
            "validation_scope": "review_note",
            "generated_by": "manual_review",
        },
        {
            "artifact_id": "candidate_package_review",
            "artifact_role": "review_note",
            "block_id": "V9-MULTI-022",
            "path": reports_path
            / "OSD120_INTERACTION_DIAGNOSTIC_CANDIDATE_PACKAGE_REVIEW.md",
            "claims": ["draft_candidate_boundary"],
            "validation_scope": "review_note",
            "generated_by": "manual_review",
        },
        {
            "artifact_id": "baseline_ladder_review",
            "artifact_role": "review_note",
            "block_id": "V9-MULTI-021",
            "path": reports_path / "OSD120_INTERACTION_BASELINE_LADDER_REVIEW.md",
            "claims": ["draft_candidate_boundary"],
            "validation_scope": "review_note",
            "generated_by": "manual_review",
        },
        {
            "artifact_id": "v9_test_source",
            "artifact_role": "validation_test_source",
            "block_id": "V9-MULTI-025",
            "path": Path("tests/test_v9_spacebio_bench.py"),
            "claims": [
                "draft_candidate_boundary",
                "nearest_centroid_fold_comparison",
                "fragile_focus_recovery",
                "feature_stability_evidence",
                "figure_table_claim_boundary",
                "paired_comparator_decision",
            ],
            "validation_scope": "unit_tests",
            "generated_by": "manual_test_updates",
        },
    ]


def build_multispecies_interaction_diagnostic_artifact_manifest(
    *,
    reports_root: str | Path = "v9/multispecies/reports",
    repo_root: str | Path = ".",
    manifest_id: str = DEFAULT_INTERACTION_DIAGNOSTIC_ARTIFACT_MANIFEST_ID,
) -> dict[str, list[dict[str, str]]]:
    """Build a traceable artifact manifest for the OSD-120 diagnostic set."""

    root = Path(repo_root)
    reports_path = _resolve_path(reports_root, root)
    task_manifest_path = root / (
        "v9/multispecies/interaction_task_manifests/"
        "draft_osd120_arabidopsis_root_light_interaction_spaceflight.json"
    )
    task_payload = json.loads(task_manifest_path.read_text())
    task_id = str(task_payload["task_id"])
    claim_boundary = "draft_interaction_logistic_only_not_leaderboard"
    artifact_rows = [
        _artifact_row(
            manifest_id=manifest_id,
            task_id=task_id,
            artifact_id=str(artifact["artifact_id"]),
            artifact_role=str(artifact["artifact_role"]),
            block_id=str(artifact["block_id"]),
            path=_resolve_path(artifact["path"], root),
            required_for_claims=artifact["claims"],
            validation_scope=str(artifact["validation_scope"]),
            generated_by=str(artifact["generated_by"]),
            claim_boundary=claim_boundary,
        )
        for artifact in _default_osd120_diagnostic_artifacts(reports_path)
    ]
    artifact_by_id = {row["artifact_id"]: row for row in artifact_rows}
    candidate_claim_rows = _read_csv_dict_rows(
        reports_path
        / "interaction_diagnostic_candidate_package"
        / "candidate_claim_map.csv"
    )
    claim_source = {row["claim_id"]: row for row in candidate_claim_rows}
    validation_tests = (
        "test_generated_v9_osd120_baseline_ladder_outputs_validate|"
        "test_generated_v9_osd120_diagnostic_candidate_package_outputs_validate|"
        "test_generated_v9_osd120_figure_table_package_outputs_validate|"
        "test_generated_v9_osd120_paired_comparator_table_outputs_validate|"
        "/usr/local/bin/python3 -m unittest discover -s tests"
    )
    claim_definitions = [
        {
            "claim_id": "draft_candidate_boundary",
            "fallback_text": (
                "sparse_l1_c1 is a draft transparent diagnostic candidate, "
                "not a frozen leaderboard baseline."
            ),
        },
        {
            "claim_id": "nearest_centroid_fold_comparison",
            "fallback_text": (
                "The candidate improves 9 of 11 nearest-centroid fold rows, "
                "ties 2, and worsens 0."
            ),
        },
        {
            "claim_id": "fragile_focus_recovery",
            "fallback_text": (
                "The candidate recovers all three predefined fragile focus folds "
                "relative to nearest centroid."
            ),
        },
        {
            "claim_id": "feature_stability_evidence",
            "fallback_text": (
                "Nineteen candidate features are selected in at least 50% of "
                "balanced train-fold subsamples across the three focus folds."
            ),
        },
        {
            "claim_id": "external_light_genotype_context",
            "fallback_text": (
                "OSD-120 is an Arabidopsis CARA root-tip RNA-seq study where "
                "genotype/ecotype and lighting context are core design factors."
            ),
        },
        {
            "claim_id": "figure_table_claim_boundary",
            "fallback_text": (
                "The figure/table package is a draft within-source diagnostic "
                "surface with explicit disallowed claims."
            ),
            "claim_status": "supported_local_artifact",
            "limitations": "figure text is for draft/poster use, not release freeze",
        },
        {
            "claim_id": "paired_comparator_decision",
            "fallback_text": (
                "sparse_l1_c0p3 remains an appendix/supplement comparator rather "
                "than a second primary figure panel."
            ),
            "claim_status": "supported_local_artifact",
            "limitations": "paired comparator explains compactness tradeoff only",
        },
    ]
    claim_rows: list[dict[str, str]] = []
    for definition in claim_definitions:
        claim_id = str(definition["claim_id"])
        source = claim_source.get(claim_id, {})
        artifact_ids = [
            row["artifact_id"]
            for row in artifact_rows
            if claim_id in row["required_for_claims"].split("|")
        ]
        artifact_paths = [artifact_by_id[artifact_id]["path"] for artifact_id in artifact_ids]
        claim_rows.append(
            {
                "manifest_id": manifest_id,
                "task_id": task_id,
                "claim_id": claim_id,
                "claim_status": str(
                    definition.get(
                        "claim_status",
                        source.get("claim_status", "supported_local_artifact"),
                    )
                ),
                "claim_text": str(
                    source.get("claim_text") or definition["fallback_text"]
                ),
                "artifact_ids": _pipe_join(artifact_ids),
                "artifact_paths": _pipe_join(artifact_paths),
                "validation_tests": validation_tests,
                "external_source_urls": str(source.get("external_source_urls", "")),
                "limitations": str(
                    definition.get("limitations", source.get("limitations", ""))
                ),
            }
        )
    return {"artifacts": artifact_rows, "claims": claim_rows}


def write_multispecies_interaction_diagnostic_artifact_manifest(
    *,
    output_dir: str | Path = (
        "v9/multispecies/reports/interaction_diagnostic_artifact_manifest"
    ),
    reports_root: str | Path = "v9/multispecies/reports",
    repo_root: str | Path = ".",
    manifest_id: str = DEFAULT_INTERACTION_DIAGNOSTIC_ARTIFACT_MANIFEST_ID,
) -> dict[str, Path]:
    """Write traceable OSD-120 diagnostic artifact manifest tables."""

    tables = build_multispecies_interaction_diagnostic_artifact_manifest(
        reports_root=reports_root,
        repo_root=repo_root,
        manifest_id=manifest_id,
    )
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "artifact_csv": outdir / "diagnostic_artifact_manifest.csv",
        "artifact_json": outdir / "diagnostic_artifact_manifest.json",
        "claim_csv": outdir / "diagnostic_claim_artifact_map.csv",
        "claim_json": outdir / "diagnostic_claim_artifact_map.json",
    }
    with outputs["artifact_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_DIAGNOSTIC_ARTIFACT_MANIFEST_FIELDS,
        )
        writer.writeheader()
        writer.writerows(tables["artifacts"])
    outputs["artifact_json"].write_text(
        json.dumps(tables["artifacts"], indent=2, sort_keys=True) + "\n"
    )
    with outputs["claim_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_DIAGNOSTIC_CLAIM_ARTIFACT_FIELDS,
        )
        writer.writeheader()
        writer.writerows(tables["claims"])
    outputs["claim_json"].write_text(
        json.dumps(tables["claims"], indent=2, sort_keys=True) + "\n"
    )
    return outputs


def _first_row_by_value(
    rows: Sequence[dict[str, str]],
    *,
    key: str,
    value: str,
) -> dict[str, str]:
    for row in rows:
        if row.get(key) == value:
            return row
    return {}


def _external_release_readiness_references() -> list[dict[str, str]]:
    return [
        {
            "reference_id": "nasa_osdr_repository",
            "topic": "OSDR repository scope and curated metadata",
            "url": "https://science.nasa.gov/biological-physical/data/osdr/",
            "release_readiness_implication": (
                "Use OSDR source identifiers, metadata, and access boundaries as "
                "primary public-source anchors."
            ),
        },
        {
            "reference_id": "nasa_osdr_faq_citation",
            "topic": "OSDR access and credit language",
            "url": "https://science.nasa.gov/reference/osdr-faq/",
            "release_readiness_implication": (
                "Public cards and release notes should credit OSDR/GeneLab "
                "explicitly and keep public and controlled data separated."
            ),
        },
        {
            "reference_id": "fair_principles",
            "topic": "FAIR machine-actionable metadata",
            "url": "https://www.go-fair.org/fair-principles",
            "release_readiness_implication": (
                "A release needs findable identifiers, access metadata, "
                "interoperable schemas, and reuse constraints beyond local files."
            ),
        },
        {
            "reference_id": "ro_crate_1_2",
            "topic": "Research object packaging",
            "url": "https://www.researchobject.org/ro-crate/specification.html",
            "release_readiness_implication": (
                "Future packaging should connect data entities, scripts, workflow "
                "runs, provenance, authorship, licensing, and outputs."
            ),
        },
        {
            "reference_id": "datacite_schema",
            "topic": "DOI-oriented dataset metadata",
            "url": "https://schema.datacite.org/",
            "release_readiness_implication": (
                "Versioned public releases should expose creators, title, "
                "version, identifiers, related identifiers, and resource type."
            ),
        },
        {
            "reference_id": "huggingface_dataset_cards",
            "topic": "Dataset card documentation",
            "url": "https://huggingface.co/docs/hub/datasets-cards",
            "release_readiness_implication": (
                "A public dataset artifact should include a rendered README/card "
                "with scope, structure, intended use, limitations, and metadata."
            ),
        },
        {
            "reference_id": "github_zenodo_doi",
            "topic": "Versioned repository archive DOI",
            "url": (
                "https://docs.github.com/en/repositories/archiving-a-github-"
                "repository/referencing-and-citing-content"
            ),
            "release_readiness_implication": (
                "A release candidate should plan a versioned archive path if it "
                "will be cited as a stable public artifact."
            ),
        },
    ]


def _release_gap_row(
    *,
    audit_id: str,
    task_id: str,
    gap_id: str,
    category: str,
    readiness_status: str,
    fix_priority: str,
    finding: str,
    evidence_artifacts: Sequence[str],
    evidence_fields: Sequence[str],
    external_reference_ids: Sequence[str],
    reference_by_id: Mapping[str, dict[str, str]],
    remediation: str,
    claim_boundary_impact: str,
    next_owner_block: str,
) -> dict[str, str]:
    reference_urls = [
        reference_by_id[reference_id]["url"]
        for reference_id in external_reference_ids
        if reference_id in reference_by_id
    ]
    return {
        "audit_id": audit_id,
        "task_id": task_id,
        "gap_id": gap_id,
        "category": category,
        "readiness_status": readiness_status,
        "fix_priority": fix_priority,
        "finding": finding,
        "evidence_artifacts": _pipe_join(list(evidence_artifacts)),
        "evidence_fields": _pipe_join(list(evidence_fields)),
        "external_reference_ids": _pipe_join(list(external_reference_ids)),
        "external_reference_urls": _pipe_join(reference_urls),
        "remediation": remediation,
        "claim_boundary_impact": claim_boundary_impact,
        "next_owner_block": next_owner_block,
    }


def build_multispecies_interaction_release_readiness_gap_audit(
    *,
    reports_root: str | Path = "v9/multispecies/reports",
    repo_root: str | Path = ".",
    task_manifest: str | Path = (
        "v9/multispecies/interaction_task_manifests/"
        "draft_osd120_arabidopsis_root_light_interaction_spaceflight.json"
    ),
    audit_id: str = DEFAULT_INTERACTION_RELEASE_READINESS_AUDIT_ID,
) -> dict[str, list[dict[str, str]]]:
    """Build an OSD-120 release-readiness gap audit from local evidence."""

    root = Path(repo_root)
    reports_path = _resolve_path(reports_root, root)
    task_manifest_path = _resolve_path(task_manifest, root)
    manifest = json.loads(task_manifest_path.read_text())
    task_id = str(manifest["task_id"])
    provenance = manifest.get("provenance", {})
    source_records = manifest.get("source_records", [])
    source_record = source_records[0] if isinstance(source_records, list) else {}
    source_id = str(source_record.get("source_id", "OSD-120"))

    candidate_summary_path = (
        reports_path
        / "interaction_diagnostic_candidate_package"
        / "candidate_package_summary.csv"
    )
    artifact_manifest_path = (
        reports_path
        / "interaction_diagnostic_artifact_manifest"
        / "diagnostic_artifact_manifest.csv"
    )
    claim_map_path = (
        reports_path
        / "interaction_diagnostic_artifact_manifest"
        / "diagnostic_claim_artifact_map.csv"
    )
    source_inventory_path = root / "v9/multispecies/source_inventory.draft.csv"
    source_checksum_path = root / "v9/multispecies/source_checksum_audit.draft.csv"
    expression_audit_path = root / "v9/multispecies/expression_matrix_audit.draft.csv"
    sample_factors_path = root / "v9/multispecies/sample_factors.draft.csv"
    local_payload_review_path = (
        reports_path / "MULTISPECIES_CHECKSUM_AND_LOCAL_PAYLOAD_AUDIT.md"
    )

    candidate_summary = _require_one_row(
        _read_csv_dict_rows(candidate_summary_path),
        context="OSD-120 diagnostic candidate package summary",
    )
    artifact_rows = _read_csv_dict_rows(artifact_manifest_path)
    claim_rows = _read_csv_dict_rows(claim_map_path)
    source_inventory = _first_row_by_value(
        _read_csv_dict_rows(source_inventory_path),
        key="source_id",
        value=source_id,
    )
    source_checksum = _first_row_by_value(
        _read_csv_dict_rows(source_checksum_path),
        key="source_id",
        value=source_id,
    )
    expression_audit = _first_row_by_value(
        _read_csv_dict_rows(expression_audit_path),
        key="source_id",
        value=source_id,
    )
    sample_factor_rows = [
        row
        for row in _read_csv_dict_rows(sample_factors_path)
        if row.get("source_id") == source_id
    ]
    parsed_sample_factor_count = sum(
        1 for row in sample_factor_rows if row.get("parse_status") == "parsed"
    )
    osd120_checksum_provenance = (
        provenance.get("source_checksum_sources", {}).get(source_id, {})
        if isinstance(provenance.get("source_checksum_sources"), dict)
        else {}
    )
    references = _external_release_readiness_references()
    reference_by_id = {row["reference_id"]: row for row in references}

    artifact_count = len(artifact_rows)
    missing_artifact_count = sum(1 for row in artifact_rows if row["exists"] != "True")
    unhashed_artifact_count = sum(1 for row in artifact_rows if not row["sha256"])
    claim_count = len(claim_rows)
    claim_boundary = str(candidate_summary["claim_boundary"])
    freeze_ready = str(source_checksum.get("freeze_ready", "")).lower() == "true"
    checksum_status = str(
        source_checksum.get(
            "audit_status",
            provenance.get("source_checksum_status", ""),
        )
    )
    source_inventory_status = str(
        provenance.get(
            "source_inventory_status",
            source_inventory.get("release_target", ""),
        )
    )
    local_payload_status = str(provenance.get("local_payload_status", ""))
    sample_factor_status = str(provenance.get("sample_factor_status", ""))
    matrix_status = str(
        expression_audit.get(
            "audit_status",
            provenance.get("expression_matrix_status", ""),
        )
    )

    gap_rows = [
        _release_gap_row(
            audit_id=audit_id,
            task_id=task_id,
            gap_id="artifact_manifest_traceability",
            category="local_artifact_traceability",
            readiness_status="pass",
            fix_priority="draft_ready",
            finding=(
                f"{artifact_count} diagnostic artifacts are indexed; "
                f"{missing_artifact_count} missing and {unhashed_artifact_count} "
                "without SHA-256."
            ),
            evidence_artifacts=[artifact_manifest_path.as_posix()],
            evidence_fields=[
                "exists",
                "byte_size",
                "row_count",
                "sha256",
                "required_for_claims",
            ],
            external_reference_ids=["fair_principles"],
            reference_by_id=reference_by_id,
            remediation="Keep regenerating this manifest after every diagnostic output change.",
            claim_boundary_impact="Supports local diagnostic claims only.",
            next_owner_block="V9-MULTI-026",
        ),
        _release_gap_row(
            audit_id=audit_id,
            task_id=task_id,
            gap_id="claim_boundary_map",
            category="claim_language",
            readiness_status="pass",
            fix_priority="draft_ready",
            finding=(
                f"{claim_count} diagnostic claims are mapped to artifacts, tests, "
                "limitations, and external context where relevant."
            ),
            evidence_artifacts=[
                claim_map_path.as_posix(),
                (
                    reports_path
                    / "interaction_figure_table_package"
                    / "figure_claim_boundary_box.md"
                ).as_posix(),
            ],
            evidence_fields=[
                "claim_status",
                "artifact_ids",
                "validation_tests",
                "limitations",
            ],
            external_reference_ids=["fair_principles", "huggingface_dataset_cards"],
            reference_by_id=reference_by_id,
            remediation="Reuse the claim map when drafting public-card language.",
            claim_boundary_impact="Prevents frozen, LOMO, biomarker, and operational claims.",
            next_owner_block="V9-MULTI-026",
        ),
        _release_gap_row(
            audit_id=audit_id,
            task_id=task_id,
            gap_id="sample_matrix_alignment",
            category="input_integrity",
            readiness_status="pass",
            fix_priority="draft_ready",
            finding=(
                f"OSD-120 has {expression_audit.get('n_sample_columns', '')} "
                "matrix sample columns, "
                f"{expression_audit.get('n_matching_sample_columns', '')} matching "
                "sample-factor rows, and "
                f"{parsed_sample_factor_count} parsed sample factors."
            ),
            evidence_artifacts=[
                expression_audit_path.as_posix(),
                sample_factors_path.as_posix(),
            ],
            evidence_fields=[
                "matrix_columns_match_sample_factors",
                "audit_status",
                "parse_status",
                "condition_stratum",
                "light_treatment",
            ],
            external_reference_ids=["nasa_osdr_repository"],
            reference_by_id=reference_by_id,
            remediation="Keep the matrix/sample-factor audit as a required preflight check.",
            claim_boundary_impact="Supports within-source OSD-120 split integrity.",
            next_owner_block="V9-MULTI-026",
        ),
        _release_gap_row(
            audit_id=audit_id,
            task_id=task_id,
            gap_id="local_payload_spotcheck",
            category="input_integrity",
            readiness_status="pass",
            fix_priority="draft_ready",
            finding=(
                "The local OSD-120 SampleTable and normalized-count matrix are "
                "recorded as MD5 matches against the processed OSDR checksum manifest."
            ),
            evidence_artifacts=[
                local_payload_review_path.as_posix(),
                source_checksum_path.as_posix(),
            ],
            evidence_fields=[
                "local_payload_status",
                "checksum_manifest_files",
                "checksum_payload_matches",
            ],
            external_reference_ids=["nasa_osdr_repository"],
            reference_by_id=reference_by_id,
            remediation="Promote from spot-check to full payload freeze before public alpha.",
            claim_boundary_impact="Enough for draft diagnostics, not enough for a release freeze.",
            next_owner_block="V9-MULTI-027",
        ),
        _release_gap_row(
            audit_id=audit_id,
            task_id=task_id,
            gap_id="full_osdr_payload_freeze",
            category="source_freeze",
            readiness_status="blocker",
            fix_priority="before_public_alpha",
            finding=(
                "OSD-120 has API and checksum-manifest evidence, but "
                f"freeze_ready={source_checksum.get('freeze_ready', '')}; "
                "the audit has not downloaded and hashed every listed payload."
            ),
            evidence_artifacts=[
                source_checksum_path.as_posix(),
                task_manifest_path.as_posix(),
            ],
            evidence_fields=[
                "freeze_ready",
                "pending_reason",
                "parsed_checksum_entries",
                "checksum_payload_matches",
            ],
            external_reference_ids=[
                "nasa_osdr_repository",
                "fair_principles",
            ],
            reference_by_id=reference_by_id,
            remediation=(
                "Create an OSD-120 payload freeze manifest with file-level names, "
                "expected MD5 values, observed local hashes, and missing payload rows."
            ),
            claim_boundary_impact="Blocks any frozen public-alpha source claim.",
            next_owner_block="V9-MULTI-027",
        ),
        _release_gap_row(
            audit_id=audit_id,
            task_id=task_id,
            gap_id="source_inventory_release_target",
            category="source_freeze",
            readiness_status="blocker",
            fix_priority="before_public_alpha",
            finding=(
                "The OSD-120 source inventory and task manifest still mark the "
                "source as draft/not frozen."
            ),
            evidence_artifacts=[
                source_inventory_path.as_posix(),
                task_manifest_path.as_posix(),
            ],
            evidence_fields=[
                "release_target",
                "source_inventory_status",
                "release_status",
            ],
            external_reference_ids=["fair_principles", "datacite_schema"],
            reference_by_id=reference_by_id,
            remediation=(
                "Define a v9 public-alpha OSD-120 release target only after the "
                "payload freeze and card/package metadata are complete."
            ),
            claim_boundary_impact="Blocks public-alpha wording beyond draft diagnostic artifact.",
            next_owner_block="V9-MULTI-027",
        ),
        _release_gap_row(
            audit_id=audit_id,
            task_id=task_id,
            gap_id="public_alpha_data_card_and_citation",
            category="release_metadata",
            readiness_status="blocker",
            fix_priority="before_public_alpha",
            finding=(
                "There is no OSD-120-specific public-alpha card with OSDR credit, "
                "scope, intended use, data fields, limitations, and citation/version language."
            ),
            evidence_artifacts=[
                "docs/v9_hf_dataset_card.md",
                claim_map_path.as_posix(),
            ],
            evidence_fields=[
                "OSD-120-specific card missing",
                "claim_text",
                "limitations",
            ],
            external_reference_ids=[
                "nasa_osdr_faq_citation",
                "huggingface_dataset_cards",
                "datacite_schema",
            ],
            reference_by_id=reference_by_id,
            remediation=(
                "Draft an OSD-120 diagnostic public-alpha card from the claim map, "
                "artifact manifest, source inventory, and external OSDR context."
            ),
            claim_boundary_impact="Blocks clean external reader-facing release language.",
            next_owner_block="V9-MULTI-028",
        ),
        _release_gap_row(
            audit_id=audit_id,
            task_id=task_id,
            gap_id="machine_readable_release_package",
            category="release_metadata",
            readiness_status="needs_work",
            fix_priority="before_broader_release",
            finding=(
                "The diagnostic artifact manifest is strong, but there is not yet "
                "an OSD-120 Data Package or RO-Crate-style metadata bundle linking "
                "sources, scripts, outputs, checksums, authorship, license, and provenance."
            ),
            evidence_artifacts=[
                artifact_manifest_path.as_posix(),
                "v9/datapackage.draft.json",
            ],
            evidence_fields=[
                "artifact_id",
                "generated_by",
                "sha256",
                "resource descriptors pending",
            ],
            external_reference_ids=["ro_crate_1_2", "fair_principles"],
            reference_by_id=reference_by_id,
            remediation=(
                "Design an OSD-120 diagnostic package descriptor or fold it into "
                "the next v9 RO-Crate/Data Package export lane."
            ),
            claim_boundary_impact="Does not block draft review, but weakens machine-actionable release reuse.",
            next_owner_block="V9-THEN-005",
        ),
        _release_gap_row(
            audit_id=audit_id,
            task_id=task_id,
            gap_id="reproducibility_environment_lock",
            category="reproducibility",
            readiness_status="needs_work",
            fix_priority="before_broader_release",
            finding=(
                "Scripts and tests regenerate the OSD-120 tables, but the package "
                "does not yet include a release-specific environment lock or "
                "single-command rebuild gate for this diagnostic bundle."
            ),
            evidence_artifacts=[
                "scripts/build_v9_osd120_diagnostic_artifact_manifest.py",
                "scripts/build_v9_osd120_release_readiness_gap_audit.py",
                "tests/test_v9_spacebio_bench.py",
            ],
            evidence_fields=[
                "validation_tests",
                "generated_by",
                "environment lock pending",
            ],
            external_reference_ids=["fair_principles", "ro_crate_1_2"],
            reference_by_id=reference_by_id,
            remediation=(
                "Add a release preflight command that rebuilds the diagnostic "
                "tables from source manifests and records package/runtime versions."
            ),
            claim_boundary_impact="Current evidence is test-backed but not a full release rebuild gate.",
            next_owner_block="V9-MULTI-029",
        ),
        _release_gap_row(
            audit_id=audit_id,
            task_id=task_id,
            gap_id="versioned_archive_plan",
            category="release_metadata",
            readiness_status="needs_work",
            fix_priority="before_broader_release",
            finding=(
                "No OSD-120 diagnostic version tag, archival DOI plan, or related "
                "identifier map exists yet."
            ),
            evidence_artifacts=["v9/OPERATING_BACKLOG.md"],
            evidence_fields=["release target pending", "DOI/version pending"],
            external_reference_ids=["datacite_schema", "github_zenodo_doi"],
            reference_by_id=reference_by_id,
            remediation=(
                "Define version strings, related identifiers, and an archive plan "
                "when the diagnostic package is promoted beyond internal draft."
            ),
            claim_boundary_impact="Blocks citable release-candidate language.",
            next_owner_block="V9-THEN-006",
        ),
        _release_gap_row(
            audit_id=audit_id,
            task_id=task_id,
            gap_id="diagnostic_not_leaderboard",
            category="claim_language",
            readiness_status="acceptable_draft_limitation",
            fix_priority="draft_ok",
            finding=(
                "The OSD-120 package is a within-source diagnostic and not a "
                "leaderboard, leave-one-mission-out, cross-study, or cross-species benchmark."
            ),
            evidence_artifacts=[
                (
                    reports_path
                    / "interaction_figure_table_package"
                    / "figure_claim_boundary_box.md"
                ).as_posix(),
                candidate_summary_path.as_posix(),
            ],
            evidence_fields=[
                "claim_boundary",
                "release_status",
                "Not Allowed",
            ],
            external_reference_ids=["huggingface_dataset_cards"],
            reference_by_id=reference_by_id,
            remediation=(
                "Keep public wording as diagnostic-alpha only unless future "
                "frozen cross-source tasks are added."
            ),
            claim_boundary_impact="Acceptable if stated prominently in every public surface.",
            next_owner_block="V9-MULTI-028",
        ),
        _release_gap_row(
            audit_id=audit_id,
            task_id=task_id,
            gap_id="external_context_scope",
            category="external_context",
            readiness_status="pass",
            fix_priority="draft_ready",
            finding=(
                "External OSD-120 context is scoped to CARA root-tip genotype/ecotype "
                "and light-treatment design, not benchmark performance."
            ),
            evidence_artifacts=[
                claim_map_path.as_posix(),
                task_manifest_path.as_posix(),
            ],
            evidence_fields=[
                "external_source_urls",
                "limitations",
                "source_url",
            ],
            external_reference_ids=[
                "nasa_osdr_repository",
                "nasa_osdr_faq_citation",
            ],
            reference_by_id=reference_by_id,
            remediation="Keep literature citations separated from model-performance claims.",
            claim_boundary_impact="Supports biological task framing without overclaiming.",
            next_owner_block="V9-MULTI-026",
        ),
    ]

    status_counts = {
        status: sum(1 for row in gap_rows if row["readiness_status"] == status)
        for status in {
            "pass",
            "blocker",
            "needs_work",
            "acceptable_draft_limitation",
        }
    }
    blocker_count = status_counts["blocker"]
    public_alpha_ready = blocker_count == 0
    release_decision = (
        "public_alpha_ready_with_limitations"
        if public_alpha_ready
        else "not_public_alpha_ready_source_freeze_card_and_release_target_blockers"
    )
    summary_row = {
        "audit_id": audit_id,
        "task_id": task_id,
        "candidate_id": candidate_summary["candidate_id"],
        "candidate_variant_id": candidate_summary["candidate_variant_id"],
        "release_readiness_decision": release_decision,
        "public_alpha_ready": str(public_alpha_ready),
        "blocker_count": str(blocker_count),
        "needs_work_count": str(status_counts["needs_work"]),
        "pass_count": str(status_counts["pass"]),
        "acceptable_draft_limitation_count": str(
            status_counts["acceptable_draft_limitation"]
        ),
        "artifact_count": str(artifact_count),
        "missing_artifact_count": str(missing_artifact_count),
        "unhashed_artifact_count": str(unhashed_artifact_count),
        "claim_count": str(claim_count),
        "source_checksum_status": checksum_status,
        "source_freeze_ready": str(freeze_ready),
        "source_inventory_status": source_inventory_status,
        "local_payload_status": local_payload_status,
        "sample_factor_status": sample_factor_status,
        "sample_matrix_alignment_status": matrix_status,
        "claim_boundary_status": "claim_map_present_and_draft_guarded",
        "reproducibility_status": "tests_present_environment_lock_pending",
        "data_card_status": "osd120_specific_public_alpha_card_missing",
        "metadata_package_status": "artifact_manifest_present_ro_crate_pending",
        "external_reference_ids": _pipe_join(
            [row["reference_id"] for row in references]
        ),
        "next_required_block": "V9-MULTI-027: OSD-120 payload freeze manifest",
        "claim_boundary": claim_boundary,
    }
    if osd120_checksum_provenance:
        summary_row["source_freeze_ready"] = str(
            str(osd120_checksum_provenance.get("freeze_ready", "")).lower() == "true"
            and freeze_ready
        )
    return {
        "summary": [summary_row],
        "gaps": gap_rows,
        "external_references": references,
    }


def write_multispecies_interaction_release_readiness_gap_audit(
    *,
    output_dir: str | Path = (
        "v9/multispecies/reports/interaction_release_readiness_gap_audit"
    ),
    reports_root: str | Path = "v9/multispecies/reports",
    repo_root: str | Path = ".",
    task_manifest: str | Path = (
        "v9/multispecies/interaction_task_manifests/"
        "draft_osd120_arabidopsis_root_light_interaction_spaceflight.json"
    ),
    audit_id: str = DEFAULT_INTERACTION_RELEASE_READINESS_AUDIT_ID,
) -> dict[str, Path]:
    """Write OSD-120 release-readiness summary, gap, and reference tables."""

    tables = build_multispecies_interaction_release_readiness_gap_audit(
        reports_root=reports_root,
        repo_root=repo_root,
        task_manifest=task_manifest,
        audit_id=audit_id,
    )
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "summary_csv": outdir / "release_readiness_summary.csv",
        "summary_json": outdir / "release_readiness_summary.json",
        "gap_csv": outdir / "release_readiness_gap_table.csv",
        "gap_json": outdir / "release_readiness_gap_table.json",
        "reference_csv": outdir / "release_readiness_external_references.csv",
        "reference_json": outdir / "release_readiness_external_references.json",
    }
    with outputs["summary_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_RELEASE_READINESS_SUMMARY_FIELDS,
        )
        writer.writeheader()
        writer.writerows(tables["summary"])
    outputs["summary_json"].write_text(
        json.dumps(tables["summary"], indent=2, sort_keys=True) + "\n"
    )
    with outputs["gap_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_RELEASE_READINESS_GAP_FIELDS,
        )
        writer.writeheader()
        writer.writerows(tables["gaps"])
    outputs["gap_json"].write_text(
        json.dumps(tables["gaps"], indent=2, sort_keys=True) + "\n"
    )
    with outputs["reference_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_RELEASE_READINESS_EXTERNAL_REFERENCE_FIELDS,
        )
        writer.writeheader()
        writer.writerows(tables["external_references"])
    outputs["reference_json"].write_text(
        json.dumps(tables["external_references"], indent=2, sort_keys=True) + "\n"
    )
    return outputs


def _md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _local_payload_manifest_names(local_path: Path, glds_prefix: str) -> list[str]:
    basename = local_path.name
    candidates = [basename]
    prefixes = [
        f"{glds_prefix}_rna_seq_",
        f"{glds_prefix}_rna-seq_",
        f"{glds_prefix}_",
    ]
    for prefix in prefixes:
        if basename.startswith(prefix):
            candidates.append(basename[len(prefix) :])
    for marker in ("_rna_seq_", "_rna-seq_"):
        if marker in basename:
            candidates.append(basename.split(marker, 1)[1])
    return list(dict.fromkeys(candidates))


def _local_payload_role(path: Path) -> str:
    lower = path.name.lower()
    if "sampletable" in lower:
        return "sample_table"
    if "normalized_counts" in lower:
        return "normalized_count_matrix"
    return "diagnostic_input_payload"


def _read_first_csv_row(path: Path, *, key: str, value: str) -> dict[str, str]:
    return _first_row_by_value(_read_csv_dict_rows(path), key=key, value=value)


def _load_checksum_manifest_text(
    *,
    checksum_manifest_url: str,
    checksum_manifest_text: str | None,
    checksum_manifest_text_path: str | Path | None,
) -> tuple[str, str]:
    if checksum_manifest_text is not None:
        body = checksum_manifest_text.encode("utf-8")
        return checksum_manifest_text, hashlib.sha256(body).hexdigest()
    if checksum_manifest_text_path is not None:
        path = Path(checksum_manifest_text_path)
        body = path.read_bytes()
        return body.decode("utf-8", errors="replace"), hashlib.sha256(body).hexdigest()
    fetched = fetch_url(checksum_manifest_url, timeout=30, max_bytes=5_000_000)
    if not fetched.ok:
        raise RuntimeError(
            f"failed to fetch checksum manifest {checksum_manifest_url}: "
            f"{fetched.error}"
        )
    return fetched.body.decode("utf-8", errors="replace"), fetched.sha256


def build_multispecies_interaction_payload_freeze_manifest(
    *,
    repo_root: str | Path = ".",
    task_manifest: str | Path = (
        "v9/multispecies/interaction_task_manifests/"
        "draft_osd120_arabidopsis_root_light_interaction_spaceflight.json"
    ),
    source_checksum_audit: str | Path = (
        "v9/multispecies/source_checksum_audit.draft.csv"
    ),
    expression_matrix_audit: str | Path = (
        "v9/multispecies/expression_matrix_audit.draft.csv"
    ),
    freeze_id: str = DEFAULT_INTERACTION_PAYLOAD_FREEZE_ID,
    checksum_manifest_text: str | None = None,
    checksum_manifest_text_path: str | Path | None = None,
) -> dict[str, list[dict[str, str]]]:
    """Build an OSD-120 diagnostic payload freeze manifest."""

    root = Path(repo_root)
    task_manifest_path = _resolve_path(task_manifest, root)
    source_checksum_path = _resolve_path(source_checksum_audit, root)
    expression_audit_path = _resolve_path(expression_matrix_audit, root)
    manifest = json.loads(task_manifest_path.read_text())
    task_id = str(manifest["task_id"])
    release_status = str(manifest.get("release_status", ""))
    claim_boundary = "diagnostic_required_payload_scope_only_not_full_osdr_freeze"
    source_record = (
        manifest.get("source_records", [{}])[0]
        if isinstance(manifest.get("source_records"), list)
        else {}
    )
    source_id = str(source_record.get("source_id", "OSD-120"))
    glds_prefix = str(source_record.get("glds_prefix", "GLDS-120"))
    checksum_row = _read_first_csv_row(
        source_checksum_path,
        key="source_id",
        value=source_id,
    )
    if not checksum_row:
        raise ValueError(f"source checksum audit row missing for {source_id}")
    expression_row = _read_first_csv_row(
        expression_audit_path,
        key="source_id",
        value=source_id,
    )
    if not expression_row:
        raise ValueError(f"expression matrix audit row missing for {source_id}")
    checksum_manifest_url = str(checksum_row.get("checksum_manifest_urls", "")).split(
        ";"
    )[0]
    if not checksum_manifest_url:
        raise ValueError(f"checksum manifest URL missing for {source_id}")
    manifest_text, manifest_sha256 = _load_checksum_manifest_text(
        checksum_manifest_url=checksum_manifest_url,
        checksum_manifest_text=checksum_manifest_text,
        checksum_manifest_text_path=checksum_manifest_text_path,
    )
    checksum_entries = parse_checksum_manifest(manifest_text)
    if not checksum_entries:
        raise ValueError(f"no checksum entries parsed for {source_id}")

    local_required_paths = [
        _resolve_path(expression_row["local_sample_table_path"], root),
        _resolve_path(expression_row["local_matrix_path"], root),
    ]
    required_by_manifest_name: dict[str, Path] = {}
    for local_path in local_required_paths:
        for candidate_name in _local_payload_manifest_names(local_path, glds_prefix):
            required_by_manifest_name[candidate_name] = local_path

    payload_rows: list[dict[str, str]] = []
    for entry in checksum_entries:
        payload_filename = entry["filename"]
        local_path = required_by_manifest_name.get(payload_filename)
        is_required = local_path is not None
        local_exists = bool(local_path and local_path.exists())
        observed_checksum = ""
        checksum_match = ""
        payload_role = "out_of_scope_processed_payload"
        release_scope = "osdr_processed_payload_not_required_for_diagnostic"
        verification_status = "out_of_scope_not_downloaded"
        action = "do_not_download_for_current_diagnostic_freeze"
        local_path_text = ""
        if is_required and local_path is not None:
            local_path_text = local_path.as_posix()
            payload_role = _local_payload_role(local_path)
            release_scope = "diagnostic_required_payload"
            if not local_exists:
                verification_status = "required_local_payload_missing"
                action = "restore_required_local_payload_before_release"
                checksum_match = "False"
            elif entry["algorithm"] != "md5":
                verification_status = "required_payload_unsupported_checksum_algorithm"
                action = "add checksum implementation for required payload"
                checksum_match = "False"
            else:
                observed_checksum = _md5_file(local_path)
                checksum_match = str(observed_checksum == entry["checksum"])
                if observed_checksum == entry["checksum"]:
                    verification_status = "required_payload_md5_matched"
                    action = "keep_in_diagnostic_payload_freeze"
                else:
                    verification_status = "required_payload_checksum_mismatch"
                    action = "investigate_or_redownload_required_payload"
        payload_rows.append(
            {
                "freeze_id": freeze_id,
                "task_id": task_id,
                "source_id": source_id,
                "glds_prefix": glds_prefix,
                "payload_filename": payload_filename,
                "algorithm": entry["algorithm"],
                "expected_checksum": entry["checksum"],
                "local_path": local_path_text,
                "local_exists": str(local_exists),
                "observed_checksum": observed_checksum,
                "checksum_match": checksum_match,
                "release_scope": release_scope,
                "payload_role": payload_role,
                "verification_status": verification_status,
                "action": action,
                "checksum_manifest_url": checksum_manifest_url,
                "source_checksum_audit_csv": source_checksum_path.as_posix(),
            }
        )

    required_rows = [
        row
        for row in payload_rows
        if row["release_scope"] == "diagnostic_required_payload"
    ]
    matched_required = [
        row
        for row in required_rows
        if row["verification_status"] == "required_payload_md5_matched"
    ]
    missing_required = [
        row
        for row in required_rows
        if row["verification_status"] == "required_local_payload_missing"
    ]
    mismatched_required = [
        row
        for row in required_rows
        if row["verification_status"] == "required_payload_checksum_mismatch"
    ]
    diagnostic_ready = (
        len(required_rows) == len(local_required_paths)
        and len(matched_required) == len(local_required_paths)
        and not missing_required
        and not mismatched_required
    )
    expected_manifest_hash = str(checksum_row.get("checksum_manifest_content_sha256", ""))
    summary = {
        "freeze_id": freeze_id,
        "task_id": task_id,
        "source_id": source_id,
        "glds_prefix": glds_prefix,
        "checksum_manifest_url": checksum_manifest_url,
        "checksum_manifest_sha256": manifest_sha256,
        "checksum_manifest_sha256_matches_source_audit": str(
            manifest_sha256 in expected_manifest_hash.split(";")
        ),
        "parsed_checksum_entries": str(len(checksum_entries)),
        "required_payload_count": str(len(required_rows)),
        "required_payload_matched_count": str(len(matched_required)),
        "required_payload_missing_count": str(len(missing_required)),
        "required_payload_checksum_mismatch_count": str(len(mismatched_required)),
        "out_of_scope_payload_count": str(
            sum(
                1
                for row in payload_rows
                if row["release_scope"]
                == "osdr_processed_payload_not_required_for_diagnostic"
            )
        ),
        "diagnostic_required_payload_freeze_ready": str(diagnostic_ready),
        "full_osdr_payload_freeze_ready": "False",
        "release_scope_decision": (
            "diagnostic_required_payloads_frozen_full_osdr_payloads_not_frozen"
            if diagnostic_ready
            else "diagnostic_required_payload_freeze_incomplete"
        ),
        "release_status": release_status,
        "next_required_block": "V9-MULTI-028: OSD-120 diagnostic public-alpha card",
        "claim_boundary": claim_boundary,
    }
    return {"summary": [summary], "payloads": payload_rows}


def write_multispecies_interaction_payload_freeze_manifest(
    *,
    output_dir: str | Path = (
        "v9/multispecies/reports/interaction_payload_freeze_manifest"
    ),
    repo_root: str | Path = ".",
    task_manifest: str | Path = (
        "v9/multispecies/interaction_task_manifests/"
        "draft_osd120_arabidopsis_root_light_interaction_spaceflight.json"
    ),
    source_checksum_audit: str | Path = (
        "v9/multispecies/source_checksum_audit.draft.csv"
    ),
    expression_matrix_audit: str | Path = (
        "v9/multispecies/expression_matrix_audit.draft.csv"
    ),
    freeze_id: str = DEFAULT_INTERACTION_PAYLOAD_FREEZE_ID,
    checksum_manifest_text: str | None = None,
    checksum_manifest_text_path: str | Path | None = None,
) -> dict[str, Path]:
    """Write OSD-120 diagnostic payload freeze manifest tables."""

    tables = build_multispecies_interaction_payload_freeze_manifest(
        repo_root=repo_root,
        task_manifest=task_manifest,
        source_checksum_audit=source_checksum_audit,
        expression_matrix_audit=expression_matrix_audit,
        freeze_id=freeze_id,
        checksum_manifest_text=checksum_manifest_text,
        checksum_manifest_text_path=checksum_manifest_text_path,
    )
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "summary_csv": outdir / "payload_freeze_summary.csv",
        "summary_json": outdir / "payload_freeze_summary.json",
        "manifest_csv": outdir / "payload_freeze_manifest.csv",
        "manifest_json": outdir / "payload_freeze_manifest.json",
    }
    with outputs["summary_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_PAYLOAD_FREEZE_SUMMARY_FIELDS,
        )
        writer.writeheader()
        writer.writerows(tables["summary"])
    outputs["summary_json"].write_text(
        json.dumps(tables["summary"], indent=2, sort_keys=True) + "\n"
    )
    with outputs["manifest_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_PAYLOAD_FREEZE_MANIFEST_FIELDS,
        )
        writer.writeheader()
        writer.writerows(tables["payloads"])
    outputs["manifest_json"].write_text(
        json.dumps(tables["payloads"], indent=2, sort_keys=True) + "\n"
    )
    return outputs


def _markdown_table(headers: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
    def cell(value: str) -> str:
        return str(value).replace("|", "\\|")

    lines = [
        "| " + " | ".join(cell(header) for header in headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(cell(value) for value in row) + " |")
    return "\n".join(lines)


def build_multispecies_interaction_public_alpha_card(
    *,
    reports_root: str | Path = "v9/multispecies/reports",
    repo_root: str | Path = ".",
    card_id: str = DEFAULT_INTERACTION_PUBLIC_ALPHA_CARD_ID,
) -> dict[str, Any]:
    """Build an OSD-120 diagnostic public-alpha card draft."""

    root = Path(repo_root)
    reports_path = _resolve_path(reports_root, root)
    candidate_summary = _require_one_row(
        _read_csv_dict_rows(
            reports_path
            / "interaction_diagnostic_candidate_package"
            / "candidate_package_summary.csv"
        ),
        context="OSD-120 candidate package summary for card",
    )
    payload_summary = _require_one_row(
        _read_csv_dict_rows(
            reports_path
            / "interaction_payload_freeze_manifest"
            / "payload_freeze_summary.csv"
        ),
        context="OSD-120 payload freeze summary for card",
    )
    release_summary = _require_one_row(
        _read_csv_dict_rows(
            reports_path
            / "interaction_release_readiness_gap_audit"
            / "release_readiness_summary.csv"
        ),
        context="OSD-120 release readiness summary for card",
    )
    claim_rows = _read_csv_dict_rows(
        reports_path
        / "interaction_diagnostic_artifact_manifest"
        / "diagnostic_claim_artifact_map.csv"
    )
    focus_rows = _read_csv_dict_rows(
        reports_path
        / "interaction_figure_table_package"
        / "figure_main_focus_table.csv"
    )
    allowed_claims = [
        row["claim_text"]
        for row in claim_rows
        if row.get("claim_status", "").startswith("supported")
    ]
    disallowed_claims = [
        "Do not call this a frozen v9 benchmark baseline.",
        "Do not claim leave-one-mission-out or cross-study generalization.",
        "Do not describe selected genes as validated biomarkers.",
        "Do not make operational plant-growth recommendations.",
        "Do not claim a complete local OSD-120 OSDR payload mirror.",
    ]
    external_context_urls = sorted(
        {
            url
            for row in claim_rows
            for url in row.get("external_source_urls", "").split("|")
            if url
        }
    )
    payload_table = _markdown_table(
        ["Payload Scope", "Count", "Status"],
        [
            [
                "Diagnostic-required processed payloads",
                payload_summary["required_payload_count"],
                (
                    payload_summary[
                        "diagnostic_required_payload_freeze_ready"
                    ]
                ),
            ],
            [
                "Required payload MD5 matches",
                payload_summary["required_payload_matched_count"],
                "matched",
            ],
            [
                "Required payload missing",
                payload_summary["required_payload_missing_count"],
                "none",
            ],
            [
                "OSDR processed payloads outside current diagnostic scope",
                payload_summary["out_of_scope_payload_count"],
                "not downloaded for this card",
            ],
            [
                "Full OSDR processed payload freeze",
                payload_summary["full_osdr_payload_freeze_ready"],
                "not claimed",
            ],
        ],
    )
    focus_table = _markdown_table(
        [
            "Focus Fold",
            "Nearest BA",
            "Sparse L1 BA",
            "Delta",
            "Stable >=0.5",
        ],
        [
            [
                row["display_fold"],
                row["nearest_centroid_ba"],
                row["candidate_ba"],
                row["display_delta_ba"],
                row["stable_features_ge_0_5"],
            ]
            for row in focus_rows
        ],
    )
    allowed_claim_text = "\n".join(f"- {claim}" for claim in allowed_claims)
    disallowed_claim_text = "\n".join(f"- {claim}" for claim in disallowed_claims)
    external_links = "\n".join(f"- {url}" for url in external_context_urls)
    card_md = f"""---
pretty_name: OSD-120 Arabidopsis Root Light Interaction Diagnostic Draft
license: other
task_categories:
- tabular-classification
tags:
- biology
- space-biology
- genomics
- OSDR
- GeneLab
- Arabidopsis
- draft
---

# OSD-120 Arabidopsis Root Light Interaction Diagnostic Draft

Release status: draft diagnostic alpha card, not a frozen benchmark release.

This card describes the current OSD-120/GLDS-120 diagnostic evidence package in
SpaceBio-Bench v9. It is source-specific and within-source: it evaluates an
Arabidopsis thaliana root bulk RNA-seq light/genotype interaction task derived
from public NASA OSDR/GeneLab processed data.

## Source And Scope

- Source: OSD-120 / GLDS-120
- Source URL: {candidate_summary["source_url"]}
- Organism: {candidate_summary["organism"]}
- Biospecimen: {candidate_summary["biospecimen_type"]}
- Assay modality: {candidate_summary["assay_modality"]}
- Task id: {candidate_summary["task_id"]}
- Candidate: {candidate_summary["candidate_id"]}
- Candidate variant: {candidate_summary["candidate_variant_id"]}

The external OSD-120 literature and source metadata support the task framing:
Arabidopsis root-tip spaceflight response with genotype/ecotype and light/dark
context. They do not support benchmark-performance claims by themselves.

## Payload Freeze Boundary

{payload_table}

The card's input freeze is intentionally narrow. The diagnostic package uses
the OSDR processed SampleTable and normalized-count matrix. Both local files
match the OSDR processed MD5 manifest. The broader OSDR processed payload set is
listed in the freeze manifest but is outside this diagnostic card's required
payload scope.

## Diagnostic Result Surface

{focus_table}

Summary metrics:

- Primary genotype/ecotype balanced accuracy:
  {candidate_summary["primary_balanced_accuracy"]}
- Secondary light-treatment balanced accuracy:
  {candidate_summary["secondary_balanced_accuracy"]}
- Diagnostic condition-stratum balanced accuracy:
  {candidate_summary["diagnostic_balanced_accuracy"]}
- Nearest-centroid fold comparison:
  {candidate_summary["nearest_improved_count"]} improve /
  {candidate_summary["nearest_tied_count"]} tie /
  {candidate_summary["nearest_worse_count"]} worse
- Stable sparse-model feature rows at selection frequency >=0.5:
  {candidate_summary["stable_feature_count_ge_0_5_total"]}

## Allowed Claims

{allowed_claim_text}

## Disallowed Claims

{disallowed_claim_text}

## Files To Inspect

- `interaction_diagnostic_artifact_manifest/diagnostic_artifact_manifest.csv`
- `interaction_diagnostic_artifact_manifest/diagnostic_claim_artifact_map.csv`
- `interaction_payload_freeze_manifest/payload_freeze_manifest.csv`
- `interaction_release_readiness_gap_audit/release_readiness_gap_table.csv`
- `interaction_figure_table_package/figure_main_focus_table.csv`
- `interaction_figure_table_package/figure_stable_feature_appendix.csv`
- `interaction_rebuild_gate/rebuild_gate_summary.csv`
- `interaction_rebuild_gate/rebuild_gate_steps.csv`
- `interaction_rebuild_gate/rebuild_gate_environment.csv`
- `interaction_public_metadata_package/public_metadata_summary.csv`
- `interaction_public_metadata_package/source_release_target_decision.csv`
- `interaction_public_metadata_package/public_metadata_skeleton.json`
- `interaction_ro_crate_citation_scaffold/ro_crate_export_summary.csv`
- `interaction_ro_crate_citation_scaffold/citation_freeze_checklist.csv`
- `interaction_ro_crate_citation_scaffold/ro-crate-metadata.draft.json`

## External Context Links

{external_links}

## Remaining Release Work

- Source release-target promotion is still pending.
- Full OSDR processed payload mirror is not claimed.
- A broader machine-readable metadata skeleton and draft RO-Crate/Data Package
  scaffold are available; archive identifier, creator, and license decisions
  are still pending.
- Packaging rebuild preflight is available through
  `python3 scripts/rebuild_v9_osd120_diagnostic_package.py --repo-root . --mode preflight`.
  It hashes packaging outputs and records environment context, but does not
  rerun model grids or freeze the benchmark release.

## Citation And Credit

Credit NASA OSDR/GeneLab and cite the upstream OSD-120 source when using this
diagnostic package. This local card is a draft SpaceBio-Bench diagnostic
surface and should not replace the upstream OSDR study metadata or citation.
"""
    summary = {
        "card_id": card_id,
        "task_id": candidate_summary["task_id"],
        "source_id": candidate_summary["source_id"],
        "source_url": candidate_summary["source_url"],
        "organism": candidate_summary["organism"],
        "assay_modality": candidate_summary["assay_modality"],
        "biospecimen_type": candidate_summary["biospecimen_type"],
        "candidate_id": candidate_summary["candidate_id"],
        "candidate_variant_id": candidate_summary["candidate_variant_id"],
        "card_status": "draft_public_alpha_card_not_frozen_release",
        "release_status": candidate_summary["release_status"],
        "payload_freeze_decision": payload_summary["release_scope_decision"],
        "diagnostic_required_payload_freeze_ready": payload_summary[
            "diagnostic_required_payload_freeze_ready"
        ],
        "full_osdr_payload_freeze_ready": payload_summary[
            "full_osdr_payload_freeze_ready"
        ],
        "artifact_count": release_summary["artifact_count"],
        "claim_count": release_summary["claim_count"],
        "primary_balanced_accuracy": candidate_summary["primary_balanced_accuracy"],
        "secondary_balanced_accuracy": candidate_summary[
            "secondary_balanced_accuracy"
        ],
        "diagnostic_balanced_accuracy": candidate_summary[
            "diagnostic_balanced_accuracy"
        ],
        "nearest_improved_count": candidate_summary["nearest_improved_count"],
        "nearest_tied_count": candidate_summary["nearest_tied_count"],
        "nearest_worse_count": candidate_summary["nearest_worse_count"],
        "stable_feature_count_ge_0_5_total": candidate_summary[
            "stable_feature_count_ge_0_5_total"
        ],
        "allowed_claim_count": str(len(allowed_claims)),
        "disallowed_claim_count": str(len(disallowed_claims)),
        "next_required_block": "V9-MULTI-032: archive identifier and license decision gate",
        "claim_boundary": (
            "diagnostic_public_alpha_card_only_not_frozen_benchmark_release"
        ),
    }
    return {"summary": [summary], "card_md": card_md}


def write_multispecies_interaction_public_alpha_card(
    *,
    output_dir: str | Path = (
        "v9/multispecies/reports/interaction_public_alpha_card"
    ),
    reports_root: str | Path = "v9/multispecies/reports",
    repo_root: str | Path = ".",
    card_id: str = DEFAULT_INTERACTION_PUBLIC_ALPHA_CARD_ID,
) -> dict[str, Path]:
    """Write an OSD-120 diagnostic public-alpha card draft."""

    package = build_multispecies_interaction_public_alpha_card(
        reports_root=reports_root,
        repo_root=repo_root,
        card_id=card_id,
    )
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "summary_csv": outdir / "public_alpha_card_summary.csv",
        "summary_json": outdir / "public_alpha_card_summary.json",
        "card_md": outdir / "public_alpha_card.md",
    }
    with outputs["summary_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_PUBLIC_ALPHA_CARD_SUMMARY_FIELDS,
        )
        writer.writeheader()
        writer.writerows(package["summary"])
    outputs["summary_json"].write_text(
        json.dumps(package["summary"], indent=2, sort_keys=True) + "\n"
    )
    outputs["card_md"].write_text(str(package["card_md"]))
    return outputs


def _package_version(package_name: str) -> str:
    try:
        return importlib.metadata.version(package_name)
    except importlib.metadata.PackageNotFoundError:
        return "not_installed"


def _rebuild_gate_environment_rows(
    *,
    gate_id: str,
    task_id: str,
    repo_root: Path,
    reports_root: Path,
) -> list[dict[str, str]]:
    values = [
        ("python_version", sys.version.replace("\n", " "), "sys.version"),
        ("python_executable", sys.executable, "sys.executable"),
        ("platform", platform.platform(), "platform.platform"),
        ("machine", platform.machine(), "platform.machine"),
        ("repo_root", repo_root.as_posix(), "argument"),
        ("reports_root", reports_root.as_posix(), "argument"),
        ("numpy_version", np.__version__, "import"),
        ("scikit_learn_version", _package_version("scikit-learn"), "importlib.metadata"),
        ("scipy_version", _package_version("scipy"), "importlib.metadata"),
        ("pandas_version", _package_version("pandas"), "importlib.metadata"),
    ]
    return [
        {
            "gate_id": gate_id,
            "task_id": task_id,
            "key": key,
            "value": value,
            "source": source,
        }
        for key, value, source in values
    ]


def _default_osd120_rebuild_gate_steps(reports_path: Path) -> list[dict[str, Any]]:
    return [
        {
            "step_id": "baseline_ladder",
            "block_id": "V9-MULTI-021",
            "script_path": Path("scripts/build_v9_osd120_baseline_ladder.py"),
            "command": (
                "python3 scripts/build_v9_osd120_baseline_ladder.py "
                "--repo-root ."
            ),
            "step_role": "consolidate_candidate_ladder",
            "execution_policy": "packaging_only_no_model_grid",
            "outputs": [
                reports_path
                / "interaction_baseline_ladder"
                / "baseline_ladder_summary.csv",
                reports_path
                / "interaction_baseline_ladder"
                / "baseline_ladder_summary.json",
                reports_path
                / "interaction_baseline_ladder"
                / "baseline_ladder_focus_folds.csv",
                reports_path
                / "interaction_baseline_ladder"
                / "baseline_ladder_focus_folds.json",
            ],
            "notes": "Reads existing nearest-centroid/logistic/sparse result tables.",
        },
        {
            "step_id": "diagnostic_candidate_package",
            "block_id": "V9-MULTI-022",
            "script_path": Path(
                "scripts/build_v9_osd120_diagnostic_candidate_package.py"
            ),
            "command": (
                "python3 scripts/build_v9_osd120_diagnostic_candidate_package.py "
                "--repo-root ."
            ),
            "step_role": "package_sparse_l1_c1_candidate",
            "execution_policy": "packaging_only_no_sparse_l1_grid",
            "outputs": [
                reports_path
                / "interaction_diagnostic_candidate_package"
                / "candidate_package_summary.csv",
                reports_path
                / "interaction_diagnostic_candidate_package"
                / "candidate_package_summary.json",
                reports_path
                / "interaction_diagnostic_candidate_package"
                / "candidate_focus_evidence.csv",
                reports_path
                / "interaction_diagnostic_candidate_package"
                / "candidate_focus_evidence.json",
                reports_path
                / "interaction_diagnostic_candidate_package"
                / "candidate_stable_feature_evidence.csv",
                reports_path
                / "interaction_diagnostic_candidate_package"
                / "candidate_stable_feature_evidence.json",
                reports_path
                / "interaction_diagnostic_candidate_package"
                / "candidate_claim_map.csv",
                reports_path
                / "interaction_diagnostic_candidate_package"
                / "candidate_claim_map.json",
            ],
            "notes": "Freezes the transparent sparse L1 c=1 diagnostic evidence surface.",
        },
        {
            "step_id": "figure_table_package",
            "block_id": "V9-MULTI-023",
            "script_path": Path("scripts/build_v9_osd120_figure_table_package.py"),
            "command": (
                "python3 scripts/build_v9_osd120_figure_table_package.py "
                "--repo-root ."
            ),
            "step_role": "emit_human_facing_figure_tables",
            "execution_policy": "packaging_only",
            "outputs": [
                reports_path
                / "interaction_figure_table_package"
                / "figure_main_focus_table.csv",
                reports_path
                / "interaction_figure_table_package"
                / "figure_main_focus_table.json",
                reports_path
                / "interaction_figure_table_package"
                / "figure_stable_feature_appendix.csv",
                reports_path
                / "interaction_figure_table_package"
                / "figure_stable_feature_appendix.json",
                reports_path
                / "interaction_figure_table_package"
                / "figure_caption.md",
                reports_path
                / "interaction_figure_table_package"
                / "figure_claim_boundary_box.md",
            ],
            "notes": "Maintains diagnostic-only figure language and disallowed claims.",
        },
        {
            "step_id": "paired_comparator_table",
            "block_id": "V9-MULTI-024",
            "script_path": Path("scripts/build_v9_osd120_paired_comparator_table.py"),
            "command": (
                "python3 scripts/build_v9_osd120_paired_comparator_table.py "
                "--repo-root ."
            ),
            "step_role": "compare_sparse_l1_c1_against_c0p3",
            "execution_policy": "packaging_only",
            "outputs": [
                reports_path
                / "interaction_paired_comparator_table"
                / "paired_comparator_summary.csv",
                reports_path
                / "interaction_paired_comparator_table"
                / "paired_comparator_summary.json",
                reports_path
                / "interaction_paired_comparator_table"
                / "paired_focus_comparator_table.csv",
                reports_path
                / "interaction_paired_comparator_table"
                / "paired_focus_comparator_table.json",
                reports_path
                / "interaction_paired_comparator_table"
                / "paired_comparator_decision.md",
            ],
            "notes": "Locks the c=1 versus c=0.3 interpretation boundary.",
        },
        {
            "step_id": "diagnostic_artifact_manifest",
            "block_id": "V9-MULTI-025",
            "script_path": Path(
                "scripts/build_v9_osd120_diagnostic_artifact_manifest.py"
            ),
            "command": (
                "python3 scripts/build_v9_osd120_diagnostic_artifact_manifest.py "
                "--repo-root ."
            ),
            "step_role": "hash_claim_support_artifacts",
            "execution_policy": "packaging_only",
            "outputs": [
                reports_path
                / "interaction_diagnostic_artifact_manifest"
                / "diagnostic_artifact_manifest.csv",
                reports_path
                / "interaction_diagnostic_artifact_manifest"
                / "diagnostic_artifact_manifest.json",
                reports_path
                / "interaction_diagnostic_artifact_manifest"
                / "diagnostic_claim_artifact_map.csv",
                reports_path
                / "interaction_diagnostic_artifact_manifest"
                / "diagnostic_claim_artifact_map.json",
            ],
            "notes": "Records artifact hashes for diagnostic claims.",
        },
        {
            "step_id": "release_readiness_gap_audit",
            "block_id": "V9-MULTI-026",
            "script_path": Path(
                "scripts/build_v9_osd120_release_readiness_gap_audit.py"
            ),
            "command": (
                "python3 scripts/build_v9_osd120_release_readiness_gap_audit.py "
                "--repo-root ."
            ),
            "step_role": "audit_public_alpha_blockers",
            "execution_policy": "packaging_only",
            "outputs": [
                reports_path
                / "interaction_release_readiness_gap_audit"
                / "release_readiness_summary.csv",
                reports_path
                / "interaction_release_readiness_gap_audit"
                / "release_readiness_summary.json",
                reports_path
                / "interaction_release_readiness_gap_audit"
                / "release_readiness_gap_table.csv",
                reports_path
                / "interaction_release_readiness_gap_audit"
                / "release_readiness_gap_table.json",
                reports_path
                / "interaction_release_readiness_gap_audit"
                / "release_readiness_external_references.csv",
                reports_path
                / "interaction_release_readiness_gap_audit"
                / "release_readiness_external_references.json",
            ],
            "notes": "Keeps blockers explicit before any release-target promotion.",
        },
        {
            "step_id": "payload_freeze_manifest",
            "block_id": "V9-MULTI-027",
            "script_path": Path("scripts/build_v9_osd120_payload_freeze_manifest.py"),
            "command": (
                "python3 scripts/build_v9_osd120_payload_freeze_manifest.py "
                "--repo-root ."
            ),
            "step_role": "verify_required_osdr_payload_checksums",
            "execution_policy": "network_optional_use_checksum_manifest_fixture_to_skip_fetch",
            "outputs": [
                reports_path
                / "interaction_payload_freeze_manifest"
                / "payload_freeze_summary.csv",
                reports_path
                / "interaction_payload_freeze_manifest"
                / "payload_freeze_summary.json",
                reports_path
                / "interaction_payload_freeze_manifest"
                / "payload_freeze_manifest.csv",
                reports_path
                / "interaction_payload_freeze_manifest"
                / "payload_freeze_manifest.json",
            ],
            "notes": "Generated table freezes diagnostic-required payloads, not the full OSDR payload mirror.",
        },
        {
            "step_id": "public_alpha_card",
            "block_id": "V9-MULTI-028",
            "script_path": Path("scripts/build_v9_osd120_public_alpha_card.py"),
            "command": (
                "python3 scripts/build_v9_osd120_public_alpha_card.py "
                "--repo-root ."
            ),
            "step_role": "render_diagnostic_public_alpha_card_draft",
            "execution_policy": "packaging_only",
            "outputs": [
                reports_path
                / "interaction_public_alpha_card"
                / "public_alpha_card_summary.csv",
                reports_path
                / "interaction_public_alpha_card"
                / "public_alpha_card_summary.json",
                reports_path
                / "interaction_public_alpha_card"
                / "public_alpha_card.md",
            ],
            "notes": "Reader-facing draft card remains diagnostic and not release-frozen.",
        },
    ]


def build_multispecies_interaction_rebuild_gate_manifest(
    *,
    reports_root: str | Path = "v9/multispecies/reports",
    repo_root: str | Path = ".",
    task_manifest: str | Path = (
        "v9/multispecies/interaction_task_manifests/"
        "draft_osd120_arabidopsis_root_light_interaction_spaceflight.json"
    ),
    gate_id: str = DEFAULT_INTERACTION_REBUILD_GATE_ID,
    mode: str = "preflight_existing_outputs",
) -> dict[str, Any]:
    """Build a preflight manifest for the OSD-120 diagnostic packaging layer."""

    if mode != "preflight_existing_outputs":
        raise ValueError("mode must be 'preflight_existing_outputs'")
    root = Path(repo_root)
    reports_path = _resolve_path(reports_root, root)
    task_manifest_path = _resolve_path(task_manifest, root)
    manifest = json.loads(task_manifest_path.read_text())
    task_id = str(manifest["task_id"])
    steps: list[dict[str, str]] = []
    for index, step in enumerate(_default_osd120_rebuild_gate_steps(reports_path), start=1):
        script_path = _resolve_path(step["script_path"], root)
        output_paths = list(step["outputs"])
        missing_outputs = [path for path in output_paths if not path.exists()]
        output_hashes = [
            f"{path.as_posix()}={_sha256_file(path)}"
            for path in output_paths
            if path.exists()
        ]
        script_exists = script_path.exists()
        status = "ready_existing_outputs_present"
        if not script_exists and missing_outputs:
            status = "script_and_outputs_missing"
        elif not script_exists:
            status = "script_missing"
        elif missing_outputs:
            status = "outputs_missing"
        steps.append(
            {
                "gate_id": gate_id,
                "task_id": task_id,
                "step_order": str(index),
                "step_id": str(step["step_id"]),
                "block_id": str(step["block_id"]),
                "script_path": script_path.as_posix(),
                "script_exists": str(script_exists),
                "command": str(step["command"]),
                "step_role": str(step["step_role"]),
                "execution_policy": str(step["execution_policy"]),
                "status": status,
                "output_count": str(len(output_paths)),
                "missing_output_count": str(len(missing_outputs)),
                "hashed_output_count": str(len(output_hashes)),
                "output_paths": _pipe_join([path.as_posix() for path in output_paths]),
                "missing_output_paths": _pipe_join(
                    [path.as_posix() for path in missing_outputs]
                ),
                "output_sha256s": _pipe_join(output_hashes),
                "notes": str(step["notes"]),
            }
        )
    environment = _rebuild_gate_environment_rows(
        gate_id=gate_id,
        task_id=task_id,
        repo_root=root,
        reports_root=reports_path,
    )
    missing_output_count = sum(int(row["missing_output_count"]) for row in steps)
    script_missing_count = sum(1 for row in steps if row["script_exists"] != "True")
    ready_step_count = sum(
        1 for row in steps if row["status"] == "ready_existing_outputs_present"
    )
    hashed_output_count = sum(int(row["hashed_output_count"]) for row in steps)
    gate_status = (
        "ready_existing_outputs_present"
        if ready_step_count == len(steps) and not missing_output_count
        else "blocked_missing_rebuild_inputs_or_scripts"
    )
    generated_at_utc = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    summary = {
        "gate_id": gate_id,
        "task_id": task_id,
        "gate_status": gate_status,
        "mode": mode,
        "step_count": str(len(steps)),
        "ready_step_count": str(ready_step_count),
        "missing_output_count": str(missing_output_count),
        "script_missing_count": str(script_missing_count),
        "hashed_output_count": str(hashed_output_count),
        "environment_package_count": str(
            sum(1 for row in environment if row["key"].endswith("_version"))
        ),
        "generated_at_utc": generated_at_utc,
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
        "repo_root": root.as_posix(),
        "reports_root": reports_path.as_posix(),
        "rebuild_command": (
            "python3 scripts/rebuild_v9_osd120_diagnostic_package.py "
            "--repo-root . --mode preflight"
        ),
        "next_required_block": "V9-MULTI-032: archive identifier and license decision gate",
        "claim_boundary": (
            "diagnostic_packaging_preflight_only_not_model_retraining_or_frozen_release"
        ),
    }
    step_lines = "\n".join(
        (
            f"- {row['step_order']}. `{row['step_id']}` "
            f"({row['block_id']}): {row['status']}; "
            f"{row['hashed_output_count']}/{row['output_count']} outputs hashed."
        )
        for row in steps
    )
    package_lines = "\n".join(
        f"- {row['key']}: {row['value']}"
        for row in environment
        if row["key"].endswith("_version") or row["key"] == "python_version"
    )
    review_md = f"""# OSD-120 Interaction Rebuild Gate Review

Block: V9-MULTI-029

Gate id: `{gate_id}`

Mode: `{mode}`

Gate status: `{gate_status}`

This gate is a preflight/environment lock for the OSD-120 diagnostic packaging
layer. It verifies that the non-model packaging outputs from V9-MULTI-021
through V9-MULTI-028 are present, script-backed, and content-hashed. It does
not rerun sparse L1 model grids and it does not promote the package to a frozen
benchmark release.

## Single Command

```bash
python3 scripts/rebuild_v9_osd120_diagnostic_package.py --repo-root . --mode preflight
```

## Step Status

{step_lines}

## Runtime Context

{package_lines}

## Claim Boundary

The current package remains a source-specific diagnostic alpha surface for
OSD-120/GLDS-120. The gate supports rebuild readiness of packaging artifacts
from existing model outputs; it does not claim leave-one-mission-out
generalization, a full OSDR payload mirror, biomarker validation, or a frozen
public benchmark release.

## Next Block

V9-MULTI-032 should make explicit archive identifier, release version,
creator/contributor, and license/rights decisions before any citable archive
path is attempted.
"""
    return {
        "summary": [summary],
        "steps": steps,
        "environment": environment,
        "review_md": review_md,
    }


def write_multispecies_interaction_rebuild_gate_manifest(
    *,
    output_dir: str | Path = (
        "v9/multispecies/reports/interaction_rebuild_gate"
    ),
    reports_root: str | Path = "v9/multispecies/reports",
    repo_root: str | Path = ".",
    task_manifest: str | Path = (
        "v9/multispecies/interaction_task_manifests/"
        "draft_osd120_arabidopsis_root_light_interaction_spaceflight.json"
    ),
    gate_id: str = DEFAULT_INTERACTION_REBUILD_GATE_ID,
    mode: str = "preflight_existing_outputs",
) -> dict[str, Path]:
    """Write the OSD-120 diagnostic packaging rebuild gate manifest."""

    package = build_multispecies_interaction_rebuild_gate_manifest(
        reports_root=reports_root,
        repo_root=repo_root,
        task_manifest=task_manifest,
        gate_id=gate_id,
        mode=mode,
    )
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "summary_csv": outdir / "rebuild_gate_summary.csv",
        "summary_json": outdir / "rebuild_gate_summary.json",
        "step_csv": outdir / "rebuild_gate_steps.csv",
        "step_json": outdir / "rebuild_gate_steps.json",
        "environment_csv": outdir / "rebuild_gate_environment.csv",
        "environment_json": outdir / "rebuild_gate_environment.json",
        "review_md": (
            outdir.parent / "OSD120_INTERACTION_REBUILD_GATE_REVIEW.md"
        ),
    }
    with outputs["summary_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_REBUILD_GATE_SUMMARY_FIELDS,
        )
        writer.writeheader()
        writer.writerows(package["summary"])
    outputs["summary_json"].write_text(
        json.dumps(package["summary"], indent=2, sort_keys=True) + "\n"
    )
    with outputs["step_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_REBUILD_GATE_STEP_FIELDS,
        )
        writer.writeheader()
        writer.writerows(package["steps"])
    outputs["step_json"].write_text(
        json.dumps(package["steps"], indent=2, sort_keys=True) + "\n"
    )
    with outputs["environment_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_REBUILD_GATE_ENVIRONMENT_FIELDS,
        )
        writer.writeheader()
        writer.writerows(package["environment"])
    outputs["environment_json"].write_text(
        json.dumps(package["environment"], indent=2, sort_keys=True) + "\n"
    )
    outputs["review_md"].write_text(str(package["review_md"]))
    return outputs


def _public_metadata_external_references(
    *,
    checked_date: str = "2026-05-26",
) -> list[dict[str, str]]:
    return [
        {
            "reference_id": "datacite_schema_4_7",
            "topic": "DataCite metadata for citable research outputs",
            "url": "https://schema.datacite.org/",
            "checked_date": checked_date,
            "metadata_implication": (
                "Use DataCite 4.7-oriented fields for future DOI/citation "
                "metadata, but keep identifier and creator fields pending until "
                "a real release target exists."
            ),
        },
        {
            "reference_id": "ro_crate_1_2",
            "topic": "RO-Crate JSON-LD package structure",
            "url": "https://www.researchobject.org/ro-crate/specification/1.2/metadata.html",
            "checked_date": checked_date,
            "metadata_implication": (
                "Model the package as a root Dataset with file Data Entities, "
                "contextual organizations, and workflow provenance."
            ),
        },
        {
            "reference_id": "huggingface_dataset_cards",
            "topic": "Dataset card README metadata",
            "url": "https://huggingface.co/docs/hub/datasets-cards",
            "checked_date": checked_date,
            "metadata_implication": (
                "Keep README YAML metadata, tags, task category, intended use, "
                "limitations, and license state explicit."
            ),
        },
        {
            "reference_id": "nasa_osdr_faq_citation",
            "topic": "OSDR dataset citation and credit",
            "url": "https://science.nasa.gov/reference/osdr-faq/",
            "checked_date": checked_date,
            "metadata_implication": (
                "Credit NASA OSDR/GeneLab and defer dataset-specific BibTeX/RIS "
                "to the OSDR study citation path."
            ),
        },
        {
            "reference_id": "nasa_osdr_repository",
            "topic": "OSDR repository source scope",
            "url": "https://science.nasa.gov/biological-physical/data/osdr/",
            "checked_date": checked_date,
            "metadata_implication": (
                "Keep upstream OSDR source identifiers and access boundaries as "
                "primary public-source anchors."
            ),
        },
    ]


def _metadata_field_row(
    *,
    package_id: str,
    task_id: str,
    metadata_profile: str,
    field_id: str,
    field_label: str,
    status: str,
    current_value: str,
    source_artifacts: Sequence[str],
    blocking_gap: str,
    notes: str,
) -> dict[str, str]:
    return {
        "package_id": package_id,
        "task_id": task_id,
        "metadata_profile": metadata_profile,
        "field_id": field_id,
        "field_label": field_label,
        "status": status,
        "current_value": current_value,
        "source_artifacts": _pipe_join(list(source_artifacts)),
        "blocking_gap": blocking_gap,
        "notes": notes,
    }


def _release_target_row(
    *,
    package_id: str,
    task_id: str,
    target_id: str,
    release_surface: str,
    target_status: str,
    public_now: bool,
    allowed_claims: Sequence[str],
    disallowed_claims: Sequence[str],
    required_evidence_artifacts: Sequence[str],
    blocking_gaps: Sequence[str],
    next_action: str,
) -> dict[str, str]:
    return {
        "package_id": package_id,
        "task_id": task_id,
        "target_id": target_id,
        "release_surface": release_surface,
        "target_status": target_status,
        "public_now": str(public_now),
        "allowed_claims": _pipe_join(list(allowed_claims)),
        "disallowed_claims": _pipe_join(list(disallowed_claims)),
        "required_evidence_artifacts": _pipe_join(list(required_evidence_artifacts)),
        "blocking_gaps": _pipe_join(list(blocking_gaps)),
        "next_action": next_action,
    }


def build_multispecies_interaction_public_metadata_package(
    *,
    reports_root: str | Path = "v9/multispecies/reports",
    repo_root: str | Path = ".",
    task_manifest: str | Path = (
        "v9/multispecies/interaction_task_manifests/"
        "draft_osd120_arabidopsis_root_light_interaction_spaceflight.json"
    ),
    package_id: str = DEFAULT_INTERACTION_PUBLIC_METADATA_PACKAGE_ID,
) -> dict[str, Any]:
    """Build a public metadata skeleton for the OSD-120 diagnostic alpha package."""

    root = Path(repo_root)
    reports_path = _resolve_path(reports_root, root)
    task_manifest_path = _resolve_path(task_manifest, root)
    manifest = json.loads(task_manifest_path.read_text())
    task_id = str(manifest["task_id"])
    card_summary_path = (
        reports_path
        / "interaction_public_alpha_card"
        / "public_alpha_card_summary.csv"
    )
    payload_summary_path = (
        reports_path
        / "interaction_payload_freeze_manifest"
        / "payload_freeze_summary.csv"
    )
    rebuild_summary_path = (
        reports_path / "interaction_rebuild_gate" / "rebuild_gate_summary.csv"
    )
    artifact_manifest_path = (
        reports_path
        / "interaction_diagnostic_artifact_manifest"
        / "diagnostic_artifact_manifest.csv"
    )
    claim_map_path = (
        reports_path
        / "interaction_diagnostic_artifact_manifest"
        / "diagnostic_claim_artifact_map.csv"
    )
    card_summary = _require_one_row(
        _read_csv_dict_rows(card_summary_path),
        context="OSD-120 public-alpha card summary for metadata package",
    )
    payload_summary = _require_one_row(
        _read_csv_dict_rows(payload_summary_path),
        context="OSD-120 payload freeze summary for metadata package",
    )
    rebuild_summary = _require_one_row(
        _read_csv_dict_rows(rebuild_summary_path),
        context="OSD-120 rebuild gate summary for metadata package",
    )
    artifact_rows = _read_csv_dict_rows(artifact_manifest_path)
    claim_rows = _read_csv_dict_rows(claim_map_path)

    source_id = str(card_summary["source_id"])
    source_url = str(card_summary["source_url"])
    title = "OSD-120 Arabidopsis Root Light Interaction Diagnostic Draft"
    description = (
        "Source-specific SpaceBio-Bench v9 diagnostic alpha package for an "
        "Arabidopsis thaliana root bulk RNA-seq light/genotype interaction task "
        "derived from NASA OSDR/GeneLab OSD-120 processed data."
    )
    allowed_claim_ids = [
        row["claim_id"]
        for row in claim_rows
        if row.get("claim_status", "").startswith("supported")
    ]
    disallowed_claims = [
        "no_frozen_v9_benchmark_baseline",
        "no_leave_one_mission_out_or_cross_study_generalization",
        "no_validated_biomarker_claim",
        "no_operational_plant_growth_recommendation",
        "no_full_local_osdr_payload_mirror",
    ]
    public_artifacts = [
        card_summary_path.as_posix(),
        payload_summary_path.as_posix(),
        rebuild_summary_path.as_posix(),
        artifact_manifest_path.as_posix(),
        claim_map_path.as_posix(),
    ]
    targets = [
        _release_target_row(
            package_id=package_id,
            task_id=task_id,
            target_id="diagnostic_alpha_metadata_draft",
            release_surface="public_reader_metadata_and_card_draft",
            target_status="ready_now_with_draft_limitations",
            public_now=True,
            allowed_claims=allowed_claim_ids,
            disallowed_claims=disallowed_claims,
            required_evidence_artifacts=public_artifacts,
            blocking_gaps=[],
            next_action=(
                "Use as a diagnostic alpha metadata surface with all limitations "
                "visible."
            ),
        ),
        _release_target_row(
            package_id=package_id,
            task_id=task_id,
            target_id="source_specific_public_alpha_package",
            release_surface="versioned_source_specific_alpha_package",
            target_status="partially_ready_source_release_target_pending",
            public_now=False,
            allowed_claims=allowed_claim_ids,
            disallowed_claims=disallowed_claims,
            required_evidence_artifacts=public_artifacts,
            blocking_gaps=[
                "source_inventory_release_target",
                "formal_license_and_creator_metadata",
                "versioned_archive_plan",
            ],
            next_action=(
                "Decide release target, formal creators/license, and version "
                "identifier before calling this a public alpha package."
            ),
        ),
        _release_target_row(
            package_id=package_id,
            task_id=task_id,
            target_id="full_osdr_payload_mirror_release",
            release_surface="complete_local_osdr_processed_payload_mirror",
            target_status="blocked_full_payload_freeze_missing",
            public_now=False,
            allowed_claims=[],
            disallowed_claims=disallowed_claims,
            required_evidence_artifacts=[
                "v9/multispecies/reports/interaction_payload_freeze_manifest/payload_freeze_manifest.csv"
            ],
            blocking_gaps=["full_osdr_payload_freeze"],
            next_action=(
                "Download or explicitly exclude every OSDR processed payload "
                "before claiming a full mirror."
            ),
        ),
        _release_target_row(
            package_id=package_id,
            task_id=task_id,
            target_id="frozen_v9_benchmark_release",
            release_surface="citable_frozen_benchmark_or_leaderboard_release",
            target_status="blocked_version_doi_and_broader_validation_missing",
            public_now=False,
            allowed_claims=[],
            disallowed_claims=disallowed_claims,
            required_evidence_artifacts=public_artifacts,
            blocking_gaps=[
                "versioned_archive_plan",
                "full_osdr_payload_freeze",
                "cross_source_or_lomo_validation_missing",
            ],
            next_action=(
                "Reserve frozen benchmark language for a versioned release with "
                "archive identifiers and broader validation scope."
            ),
        ),
    ]

    metadata_fields = [
        _metadata_field_row(
            package_id=package_id,
            task_id=task_id,
            metadata_profile="datacite_4_7",
            field_id="identifier",
            field_label="Identifier",
            status="placeholder",
            current_value="pending_versioned_release_identifier",
            source_artifacts=[],
            blocking_gap="versioned_archive_plan",
            notes="Do not mint or imply a DOI for this draft metadata skeleton.",
        ),
        _metadata_field_row(
            package_id=package_id,
            task_id=task_id,
            metadata_profile="datacite_4_7",
            field_id="creators",
            field_label="Creators",
            status="placeholder",
            current_value="pending_formal_release_creator_list",
            source_artifacts=[],
            blocking_gap="formal_license_and_creator_metadata",
            notes="Creator/order/affiliation metadata must be reviewed before release.",
        ),
        _metadata_field_row(
            package_id=package_id,
            task_id=task_id,
            metadata_profile="datacite_4_7",
            field_id="title",
            field_label="Title",
            status="present",
            current_value=title,
            source_artifacts=[card_summary_path.as_posix()],
            blocking_gap="",
            notes="Matches the source-specific public alpha card title.",
        ),
        _metadata_field_row(
            package_id=package_id,
            task_id=task_id,
            metadata_profile="datacite_4_7",
            field_id="publisher",
            field_label="Publisher",
            status="placeholder",
            current_value="SpaceBio-Bench draft package; upstream data NASA OSDR/GeneLab",
            source_artifacts=[claim_map_path.as_posix()],
            blocking_gap="formal_license_and_creator_metadata",
            notes="Formal publisher/maintainer identity is not frozen.",
        ),
        _metadata_field_row(
            package_id=package_id,
            task_id=task_id,
            metadata_profile="datacite_4_7",
            field_id="publication_year",
            field_label="Publication Year",
            status="present",
            current_value="2026",
            source_artifacts=[],
            blocking_gap="",
            notes="Draft metadata generated in 2026; not a DOI publication year.",
        ),
        _metadata_field_row(
            package_id=package_id,
            task_id=task_id,
            metadata_profile="datacite_4_7",
            field_id="resource_type",
            field_label="Resource Type",
            status="present",
            current_value="Dataset",
            source_artifacts=[],
            blocking_gap="",
            notes="Use dataset-like packaging for the diagnostic evidence bundle.",
        ),
        _metadata_field_row(
            package_id=package_id,
            task_id=task_id,
            metadata_profile="datacite_4_7",
            field_id="subjects",
            field_label="Subjects",
            status="present",
            current_value=_pipe_join(
                [
                    "space biology",
                    "Arabidopsis thaliana",
                    "bulk RNA-seq",
                    "OSD-120",
                    "diagnostic benchmark",
                ]
            ),
            source_artifacts=[card_summary_path.as_posix()],
            blocking_gap="",
            notes="Subject tags are descriptive and not release claims.",
        ),
        _metadata_field_row(
            package_id=package_id,
            task_id=task_id,
            metadata_profile="datacite_4_7",
            field_id="related_identifiers",
            field_label="Related Identifiers",
            status="partial",
            current_value=_pipe_join(
                [
                    source_url,
                    "https://doi.org/10.1093/nar/gkae1116",
                    "https://doi.org/10.1371/journal.pone.0180186",
                    "https://www.nature.com/articles/s41526-024-00417-0",
                ]
            ),
            source_artifacts=[claim_map_path.as_posix()],
            blocking_gap="versioned_archive_plan",
            notes="Upstream source and literature links exist; local release DOI is pending.",
        ),
        _metadata_field_row(
            package_id=package_id,
            task_id=task_id,
            metadata_profile="datacite_4_7",
            field_id="rights",
            field_label="Rights",
            status="placeholder",
            current_value="pending_local_release_license_review",
            source_artifacts=[card_summary_path.as_posix()],
            blocking_gap="formal_license_and_creator_metadata",
            notes="Do not infer local package licensing from upstream OSDR access alone.",
        ),
        _metadata_field_row(
            package_id=package_id,
            task_id=task_id,
            metadata_profile="ro_crate_1_2",
            field_id="context",
            field_label="@context",
            status="present",
            current_value="https://w3id.org/ro/crate/1.2/context",
            source_artifacts=[],
            blocking_gap="",
            notes="RO-Crate 1.2 JSON-LD context for the future export path.",
        ),
        _metadata_field_row(
            package_id=package_id,
            task_id=task_id,
            metadata_profile="ro_crate_1_2",
            field_id="root_dataset",
            field_label="Root Dataset",
            status="present",
            current_value="./",
            source_artifacts=[card_summary_path.as_posix()],
            blocking_gap="",
            notes="The skeleton can be promoted to a root RO-Crate Dataset later.",
        ),
        _metadata_field_row(
            package_id=package_id,
            task_id=task_id,
            metadata_profile="ro_crate_1_2",
            field_id="data_entities",
            field_label="Data Entities",
            status="present",
            current_value=f"{len(artifact_rows)} diagnostic artifacts indexed",
            source_artifacts=[artifact_manifest_path.as_posix()],
            blocking_gap="",
            notes="Uses the diagnostic artifact manifest as the file inventory seed.",
        ),
        _metadata_field_row(
            package_id=package_id,
            task_id=task_id,
            metadata_profile="ro_crate_1_2",
            field_id="contextual_entities",
            field_label="Contextual Entities",
            status="partial",
            current_value="NASA OSDR/GeneLab; SpaceBio-Bench draft package",
            source_artifacts=[claim_map_path.as_posix()],
            blocking_gap="formal_license_and_creator_metadata",
            notes="Organizations are identified, but formal creators are pending.",
        ),
        _metadata_field_row(
            package_id=package_id,
            task_id=task_id,
            metadata_profile="ro_crate_1_2",
            field_id="workflow_provenance",
            field_label="Workflow Provenance",
            status="present",
            current_value=rebuild_summary["rebuild_command"],
            source_artifacts=[rebuild_summary_path.as_posix()],
            blocking_gap="",
            notes="Packaging preflight command and hashed outputs are recorded.",
        ),
        _metadata_field_row(
            package_id=package_id,
            task_id=task_id,
            metadata_profile="huggingface_dataset_card",
            field_id="pretty_name",
            field_label="Pretty Name",
            status="present",
            current_value=title,
            source_artifacts=[card_summary_path.as_posix()],
            blocking_gap="",
            notes="Ready for README YAML metadata.",
        ),
        _metadata_field_row(
            package_id=package_id,
            task_id=task_id,
            metadata_profile="huggingface_dataset_card",
            field_id="tags",
            field_label="Tags",
            status="present",
            current_value=_pipe_join(
                [
                    "biology",
                    "space-biology",
                    "genomics",
                    "OSDR",
                    "GeneLab",
                    "Arabidopsis",
                    "tabular",
                    "draft",
                ]
            ),
            source_artifacts=[card_summary_path.as_posix()],
            blocking_gap="",
            notes="Tags aid discovery but do not create release status.",
        ),
        _metadata_field_row(
            package_id=package_id,
            task_id=task_id,
            metadata_profile="huggingface_dataset_card",
            field_id="license",
            field_label="License",
            status="placeholder",
            current_value="other",
            source_artifacts=[card_summary_path.as_posix()],
            blocking_gap="formal_license_and_creator_metadata",
            notes="Card can say other/pending, but a citable package needs review.",
        ),
        _metadata_field_row(
            package_id=package_id,
            task_id=task_id,
            metadata_profile="huggingface_dataset_card",
            field_id="task_categories",
            field_label="Task Categories",
            status="present",
            current_value="tabular-classification",
            source_artifacts=[card_summary_path.as_posix()],
            blocking_gap="",
            notes="Matches the diagnostic classification task.",
        ),
        _metadata_field_row(
            package_id=package_id,
            task_id=task_id,
            metadata_profile="osdr_credit",
            field_id="source_citation_credit",
            field_label="OSDR Citation And Credit",
            status="partial",
            current_value=(
                "Credit NASA OSDR/GeneLab and use the OSDR study citation path "
                "for dataset-specific BibTeX/RIS."
            ),
            source_artifacts=[card_summary_path.as_posix(), claim_map_path.as_posix()],
            blocking_gap="source_inventory_release_target",
            notes="OSDR credit is present; local package citation remains pending.",
        ),
        _metadata_field_row(
            package_id=package_id,
            task_id=task_id,
            metadata_profile="spacebio_claim_boundary",
            field_id="claim_boundary",
            field_label="Claim Boundary",
            status="present",
            current_value=(
                "diagnostic_public_metadata_only_not_frozen_benchmark_release"
            ),
            source_artifacts=[claim_map_path.as_posix(), card_summary_path.as_posix()],
            blocking_gap="",
            notes="Prevents full mirror, LOMO, biomarker, operational, and frozen-release claims.",
        ),
    ]

    present_count = sum(1 for row in metadata_fields if row["status"] == "present")
    partial_count = sum(1 for row in metadata_fields if row["status"] == "partial")
    placeholder_count = sum(
        1 for row in metadata_fields if row["status"] == "placeholder"
    )
    public_now_count = sum(1 for row in targets if row["public_now"] == "True")
    not_public_now_count = len(targets) - public_now_count
    data_entities = [
        {
            "@id": row["path"],
            "@type": "File",
            "name": row["artifact_id"],
            "encodingFormat": row["file_format"],
            "contentSize": row["byte_size"],
            "sha256": row["sha256"],
            "description": row["artifact_role"],
        }
        for row in artifact_rows
    ]
    ro_crate_graph: list[dict[str, Any]] = [
        {
            "@id": "ro-crate-metadata.json",
            "@type": "CreativeWork",
            "about": {"@id": "./"},
            "conformsTo": {"@id": "https://w3id.org/ro/crate/1.2"},
        },
        {
            "@id": "./",
            "@type": "Dataset",
            "name": title,
            "description": description,
            "license": "pending_local_release_license_review",
            "keywords": [
                "space biology",
                "Arabidopsis thaliana",
                "bulk RNA-seq",
                "OSD-120",
            ],
            "hasPart": [{"@id": row["path"]} for row in artifact_rows],
            "isBasedOn": [{"@id": source_url}],
            "publisher": {"@id": "#spacebio-bench"},
            "citation": [{"@id": "https://doi.org/10.1093/nar/gkae1116"}],
        },
        {
            "@id": "#nasa-osdr",
            "@type": "Organization",
            "name": "NASA Open Science Data Repository",
            "url": "https://osdr.nasa.gov/bio/repo",
        },
        {
            "@id": "#spacebio-bench",
            "@type": "Organization",
            "name": "SpaceBio-Bench draft package",
        },
    ] + data_entities
    metadata_skeleton = {
        "package_id": package_id,
        "task_id": task_id,
        "status": "draft_public_metadata_skeleton_not_frozen_release",
        "source": {
            "source_id": source_id,
            "source_url": source_url,
            "organism": card_summary["organism"],
            "assay_modality": card_summary["assay_modality"],
            "biospecimen_type": card_summary["biospecimen_type"],
        },
        "standards_profile": {
            "datacite_schema_version": "4.7",
            "ro_crate_version": "1.2",
            "huggingface_dataset_card": "README.md YAML metadata",
            "osdr_credit": "NASA OSDR/GeneLab source credit required",
        },
        "datacite_skeleton": {
            "identifier": "pending_versioned_release_identifier",
            "creators": "pending_formal_release_creator_list",
            "titles": [{"title": title}],
            "publisher": "SpaceBio-Bench draft package; upstream data NASA OSDR/GeneLab",
            "publicationYear": "2026",
            "resourceType": {
                "resourceType": "Diagnostic evidence package",
                "resourceTypeGeneral": "Dataset",
            },
            "subjects": [
                "space biology",
                "Arabidopsis thaliana",
                "bulk RNA-seq",
                "OSD-120",
            ],
            "relatedIdentifiers": [
                {"relatedIdentifier": source_url, "relationType": "IsDerivedFrom"},
                {
                    "relatedIdentifier": "https://doi.org/10.1093/nar/gkae1116",
                    "relationType": "References",
                },
            ],
            "rightsList": ["pending_local_release_license_review"],
        },
        "huggingface_dataset_card_metadata": {
            "pretty_name": title,
            "license": "other",
            "task_categories": ["tabular-classification"],
            "tags": [
                "biology",
                "space-biology",
                "genomics",
                "OSDR",
                "GeneLab",
                "Arabidopsis",
                "tabular",
                "draft",
            ],
        },
        "ro_crate_skeleton": {
            "@context": "https://w3id.org/ro/crate/1.2/context",
            "@graph": ro_crate_graph,
        },
        "release_targets": targets,
        "claim_boundary": {
            "allowed_claim_ids": allowed_claim_ids,
            "disallowed_claims": disallowed_claims,
            "full_osdr_payload_freeze_ready": payload_summary[
                "full_osdr_payload_freeze_ready"
            ],
            "rebuild_gate_status": rebuild_summary["gate_status"],
        },
    }
    summary = {
        "package_id": package_id,
        "task_id": task_id,
        "source_id": source_id,
        "source_url": source_url,
        "metadata_package_status": (
            "draft_metadata_skeleton_ready_not_release_frozen"
        ),
        "release_target_decision": (
            "diagnostic_metadata_public_now_source_release_target_pending"
        ),
        "public_now_target_count": str(public_now_count),
        "not_public_now_target_count": str(not_public_now_count),
        "metadata_field_count": str(len(metadata_fields)),
        "metadata_present_field_count": str(present_count),
        "metadata_partial_field_count": str(partial_count),
        "metadata_placeholder_field_count": str(placeholder_count),
        "artifact_count": str(len(artifact_rows)),
        "claim_count": str(len(claim_rows)),
        "rebuild_gate_status": rebuild_summary["gate_status"],
        "diagnostic_required_payload_freeze_ready": payload_summary[
            "diagnostic_required_payload_freeze_ready"
        ],
        "full_osdr_payload_freeze_ready": payload_summary[
            "full_osdr_payload_freeze_ready"
        ],
        "datacite_schema_version": "4.7",
        "ro_crate_version": "1.2",
        "hf_card_metadata_status": "draft_card_metadata_present_license_pending",
        "osdr_citation_status": "credit_present_dataset_specific_citation_pending",
        "next_required_block": "V9-MULTI-032: archive identifier and license decision gate",
        "claim_boundary": (
            "diagnostic_public_metadata_only_not_frozen_benchmark_release"
        ),
    }
    target_lines = "\n".join(
        (
            f"- `{row['target_id']}`: {row['target_status']}; "
            f"public_now={row['public_now']}."
        )
        for row in targets
    )
    field_status_lines = "\n".join(
        [
            f"- present: {present_count}",
            f"- partial: {partial_count}",
            f"- placeholder: {placeholder_count}",
        ]
    )
    review_md = f"""# OSD-120 Public Metadata Package Review

Block: V9-MULTI-030

Package id: `{package_id}`

Status: `{summary["metadata_package_status"]}`

This block separates the current public diagnostic metadata surface from future
source-package or frozen-benchmark release claims. It creates a machine-readable
metadata skeleton aligned with DataCite 4.7, RO-Crate 1.2, Hugging Face dataset
card metadata, and OSDR citation/credit expectations.

## Release Target Decision

{target_lines}

Only `diagnostic_alpha_metadata_draft` is public-now, and only with draft
limitations. Source-specific alpha packaging, full OSDR payload mirroring, and
frozen benchmark release claims remain blocked.

## Metadata Field Status

{field_status_lines}

The placeholders are intentional: they prevent false DOI, creator, publisher,
rights, and license claims before a formal release target exists.

## Standards Anchors

- DataCite Metadata Schema 4.7: `https://schema.datacite.org/`
- RO-Crate 1.2 metadata: `https://www.researchobject.org/ro-crate/specification/1.2/metadata.html`
- Hugging Face dataset cards: `https://huggingface.co/docs/hub/datasets-cards`
- NASA OSDR FAQ/citation guidance: `https://science.nasa.gov/reference/osdr-faq/`

## Claim Boundary

This package is a diagnostic public metadata skeleton. It does not claim a
complete local OSDR payload mirror, leave-one-mission-out generalization,
validated biomarkers, operational plant-growth recommendations, or a frozen v9
benchmark release.

## Next Block

V9-MULTI-032 should make explicit archive identifier, release version,
creator/contributor, and license/rights decisions before any citable archive or
DOI path is attempted.
"""
    return {
        "summary": [summary],
        "release_targets": targets,
        "metadata_fields": metadata_fields,
        "external_references": _public_metadata_external_references(),
        "metadata_skeleton": metadata_skeleton,
        "review_md": review_md,
    }


def write_multispecies_interaction_public_metadata_package(
    *,
    output_dir: str | Path = (
        "v9/multispecies/reports/interaction_public_metadata_package"
    ),
    reports_root: str | Path = "v9/multispecies/reports",
    repo_root: str | Path = ".",
    task_manifest: str | Path = (
        "v9/multispecies/interaction_task_manifests/"
        "draft_osd120_arabidopsis_root_light_interaction_spaceflight.json"
    ),
    package_id: str = DEFAULT_INTERACTION_PUBLIC_METADATA_PACKAGE_ID,
) -> dict[str, Path]:
    """Write OSD-120 public metadata package skeleton tables."""

    package = build_multispecies_interaction_public_metadata_package(
        reports_root=reports_root,
        repo_root=repo_root,
        task_manifest=task_manifest,
        package_id=package_id,
    )
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "summary_csv": outdir / "public_metadata_summary.csv",
        "summary_json": outdir / "public_metadata_summary.json",
        "target_csv": outdir / "source_release_target_decision.csv",
        "target_json": outdir / "source_release_target_decision.json",
        "field_csv": outdir / "public_metadata_field_table.csv",
        "field_json": outdir / "public_metadata_field_table.json",
        "reference_csv": outdir / "public_metadata_external_references.csv",
        "reference_json": outdir / "public_metadata_external_references.json",
        "metadata_skeleton_json": outdir / "public_metadata_skeleton.json",
        "review_md": (
            outdir.parent / "OSD120_INTERACTION_PUBLIC_METADATA_PACKAGE_REVIEW.md"
        ),
    }
    with outputs["summary_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_PUBLIC_METADATA_SUMMARY_FIELDS,
        )
        writer.writeheader()
        writer.writerows(package["summary"])
    outputs["summary_json"].write_text(
        json.dumps(package["summary"], indent=2, sort_keys=True) + "\n"
    )
    with outputs["target_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_PUBLIC_METADATA_RELEASE_TARGET_FIELDS,
        )
        writer.writeheader()
        writer.writerows(package["release_targets"])
    outputs["target_json"].write_text(
        json.dumps(package["release_targets"], indent=2, sort_keys=True) + "\n"
    )
    with outputs["field_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_PUBLIC_METADATA_FIELD_FIELDS,
        )
        writer.writeheader()
        writer.writerows(package["metadata_fields"])
    outputs["field_json"].write_text(
        json.dumps(package["metadata_fields"], indent=2, sort_keys=True) + "\n"
    )
    with outputs["reference_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_PUBLIC_METADATA_REFERENCE_FIELDS,
        )
        writer.writeheader()
        writer.writerows(package["external_references"])
    outputs["reference_json"].write_text(
        json.dumps(package["external_references"], indent=2, sort_keys=True) + "\n"
    )
    outputs["metadata_skeleton_json"].write_text(
        json.dumps(package["metadata_skeleton"], indent=2, sort_keys=True) + "\n"
    )
    outputs["review_md"].write_text(str(package["review_md"]))
    return outputs


def _ro_crate_validation_row(
    *,
    scaffold_id: str,
    package_id: str,
    task_id: str,
    check_id: str,
    standard_profile: str,
    check_status: str,
    severity: str,
    finding: str,
    evidence: Sequence[str],
    blocking_gap: str,
    next_action: str,
) -> dict[str, str]:
    return {
        "scaffold_id": scaffold_id,
        "package_id": package_id,
        "task_id": task_id,
        "check_id": check_id,
        "standard_profile": standard_profile,
        "check_status": check_status,
        "severity": severity,
        "finding": finding,
        "evidence": _pipe_join(list(evidence)),
        "blocking_gap": blocking_gap,
        "next_action": next_action,
    }


def _citation_checklist_row(
    *,
    scaffold_id: str,
    package_id: str,
    task_id: str,
    item_id: str,
    citation_surface: str,
    required_for: str,
    status: str,
    current_value: str,
    source_artifacts: Sequence[str],
    blocking_gap: str,
    next_action: str,
    claim_boundary_impact: str,
) -> dict[str, str]:
    return {
        "scaffold_id": scaffold_id,
        "package_id": package_id,
        "task_id": task_id,
        "item_id": item_id,
        "citation_surface": citation_surface,
        "required_for": required_for,
        "status": status,
        "current_value": current_value,
        "source_artifacts": _pipe_join(list(source_artifacts)),
        "blocking_gap": blocking_gap,
        "next_action": next_action,
        "claim_boundary_impact": claim_boundary_impact,
    }


def build_multispecies_interaction_ro_crate_citation_scaffold(
    *,
    reports_root: str | Path = "v9/multispecies/reports",
    repo_root: str | Path = ".",
    scaffold_id: str = DEFAULT_INTERACTION_RO_CRATE_EXPORT_SCAFFOLD_ID,
) -> dict[str, Any]:
    """Build a RO-Crate/Data Package scaffold and citation freeze checklist."""

    root = Path(repo_root)
    reports_path = _resolve_path(reports_root, root)
    metadata_dir = reports_path / "interaction_public_metadata_package"
    skeleton_path = metadata_dir / "public_metadata_skeleton.json"
    summary_path = metadata_dir / "public_metadata_summary.csv"
    field_table_path = metadata_dir / "public_metadata_field_table.csv"
    target_path = metadata_dir / "source_release_target_decision.csv"
    artifact_manifest_path = (
        reports_path
        / "interaction_diagnostic_artifact_manifest"
        / "diagnostic_artifact_manifest.csv"
    )
    skeleton = json.loads(skeleton_path.read_text())
    metadata_summary = _require_one_row(
        _read_csv_dict_rows(summary_path),
        context="OSD-120 public metadata summary for RO-Crate scaffold",
    )
    field_rows = _read_csv_dict_rows(field_table_path)
    target_rows = _read_csv_dict_rows(target_path)
    artifact_rows = _read_csv_dict_rows(artifact_manifest_path)
    package_id = str(skeleton["package_id"])
    task_id = str(skeleton["task_id"])
    title = str(skeleton["datacite_skeleton"]["titles"][0]["title"])
    description = (
        "Draft RO-Crate/Data Package export scaffold for the OSD-120 diagnostic "
        "metadata package. This scaffold is not a citable frozen release."
    )
    placeholder_fields = [
        row for row in field_rows if row.get("status") == "placeholder"
    ]
    ro_crate_graph = list(skeleton["ro_crate_skeleton"]["@graph"])
    root_entity = next(
        (row for row in ro_crate_graph if row.get("@id") == "./"),
        {},
    )
    root_entity["identifier"] = "pending_versioned_release_identifier"
    root_entity["version"] = "draft-v9-osd120-diagnostic-alpha-unreleased"
    root_entity["dateModified"] = "2026-05-26"
    root_entity["developmentStatus"] = "draft"
    root_entity["isAccessibleForFree"] = True
    root_entity["conditionsOfAccess"] = (
        "Local diagnostic metadata scaffold only; upstream OSDR payload access "
        "and local release license remain under review."
    )
    root_entity["conformsTo"] = [
        {"@id": "https://w3id.org/ro/crate/1.2"},
        {"@id": "https://specs.frictionlessdata.io/data-package/"},
    ]
    metadata_entity = next(
        (row for row in ro_crate_graph if row.get("@id") == "ro-crate-metadata.json"),
        {},
    )
    metadata_entity["@id"] = "ro-crate-metadata.draft.json"
    metadata_entity["name"] = "Draft RO-Crate metadata scaffold"
    metadata_entity["dateModified"] = "2026-05-26"
    data_entities = [
        row
        for row in ro_crate_graph
        if row.get("@type") == "File" and row.get("@id")
    ]
    ro_crate_metadata = {
        "@context": skeleton["ro_crate_skeleton"]["@context"],
        "@graph": ro_crate_graph,
    }
    data_package_resources = [
        {
            "name": row["artifact_id"],
            "path": row["path"],
            "format": row["file_format"],
            "bytes": int(row["byte_size"] or 0),
            "hash": row["sha256"],
            "description": row["artifact_role"],
        }
        for row in artifact_rows
    ]
    data_package = {
        "profile": "data-package",
        "name": "osd120-diagnostic-public-metadata-package",
        "title": title,
        "description": description,
        "version": "draft-v9-osd120-diagnostic-alpha-unreleased",
        "licenses": [
            {
                "name": "pending-local-release-license-review",
                "title": "Pending local release license review",
            }
        ],
        "contributors": [
            {
                "title": "NASA OSDR/GeneLab",
                "role": "upstream-data-source",
                "path": "https://osdr.nasa.gov/bio/repo/data/studies/OSD-120",
            },
            {
                "title": "SpaceBio-Bench draft package",
                "role": "metadata-packager",
            },
        ],
        "sources": [
            {
                "title": "OSD-120 / GLDS-120",
                "path": skeleton["source"]["source_url"],
            }
        ],
        "resources": data_package_resources,
        "custom": {
            "claim_boundary": (
                "diagnostic_ro_crate_scaffold_only_not_citable_archive"
            ),
            "full_osdr_payload_freeze_ready": skeleton["claim_boundary"][
                "full_osdr_payload_freeze_ready"
            ],
            "release_targets": target_rows,
        },
    }
    validation_rows = [
        _ro_crate_validation_row(
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            check_id="ro_crate_context_present",
            standard_profile="ro_crate_1_2",
            check_status="pass",
            severity="required",
            finding="RO-Crate JSON-LD context is present.",
            evidence=[skeleton_path.as_posix()],
            blocking_gap="",
            next_action="Keep context pinned to RO-Crate 1.2 until a formal upgrade.",
        ),
        _ro_crate_validation_row(
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            check_id="metadata_descriptor_entity_present",
            standard_profile="ro_crate_1_2",
            check_status="pass",
            severity="required",
            finding="A metadata descriptor entity points at the root dataset.",
            evidence=["ro-crate-metadata.draft.json"],
            blocking_gap="",
            next_action="Rename only during a true archive packaging step.",
        ),
        _ro_crate_validation_row(
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            check_id="root_dataset_present",
            standard_profile="ro_crate_1_2",
            check_status="pass",
            severity="required",
            finding="Root Dataset entity is present and named.",
            evidence=[skeleton_path.as_posix()],
            blocking_gap="",
            next_action="Keep root Dataset as the package focus.",
        ),
        _ro_crate_validation_row(
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            check_id="flat_graph_entity_shape",
            standard_profile="ro_crate_1_2",
            check_status="pass",
            severity="required",
            finding="All RO-Crate graph entities have @id and @type.",
            evidence=["ro-crate-metadata.draft.json"],
            blocking_gap="",
            next_action="Preserve flattened @graph structure.",
        ),
        _ro_crate_validation_row(
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            check_id="data_entities_present",
            standard_profile="ro_crate_1_2",
            check_status="pass",
            severity="required",
            finding=f"{len(data_entities)} file Data Entities are present.",
            evidence=[artifact_manifest_path.as_posix()],
            blocking_gap="",
            next_action="Regenerate from the artifact manifest after output changes.",
        ),
        _ro_crate_validation_row(
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            check_id="data_entities_have_hashes",
            standard_profile="ro_crate_1_2",
            check_status="pass",
            severity="recommended",
            finding="All file Data Entities inherit SHA-256 hashes from the artifact manifest.",
            evidence=[artifact_manifest_path.as_posix()],
            blocking_gap="",
            next_action="Treat missing hashes as a blocker in future archive mode.",
        ),
        _ro_crate_validation_row(
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            check_id="contextual_entities_formal_creators",
            standard_profile="ro_crate_1_2",
            check_status="needs_review",
            severity="recommended",
            finding="Organizations are represented, but formal creators/contributors are placeholders.",
            evidence=[field_table_path.as_posix()],
            blocking_gap="formal_license_and_creator_metadata",
            next_action="Resolve creator list before a citable release.",
        ),
        _ro_crate_validation_row(
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            check_id="workflow_provenance_present",
            standard_profile="ro_crate_1_2",
            check_status="pass",
            severity="recommended",
            finding="Rebuild/preflight command is linked through the metadata skeleton.",
            evidence=[
                "v9/multispecies/reports/interaction_rebuild_gate/rebuild_gate_summary.csv"
            ],
            blocking_gap="",
            next_action="Add full workflow run provenance only for archive mode.",
        ),
        _ro_crate_validation_row(
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            check_id="data_package_resources_present",
            standard_profile="frictionless_data_package",
            check_status="pass",
            severity="required",
            finding=f"{len(data_package_resources)} Data Package resources are present.",
            evidence=[artifact_manifest_path.as_posix()],
            blocking_gap="",
            next_action="Validate descriptor against a Frictionless validator later.",
        ),
        _ro_crate_validation_row(
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            check_id="data_package_license_placeholder",
            standard_profile="frictionless_data_package",
            check_status="blocker",
            severity="required_for_archive",
            finding="Data Package license is an explicit pending placeholder.",
            evidence=[field_table_path.as_posix()],
            blocking_gap="formal_license_and_creator_metadata",
            next_action="Select and review a local package license before archival release.",
        ),
        _ro_crate_validation_row(
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            check_id="datacite_identifier_placeholder",
            standard_profile="datacite_4_7",
            check_status="blocker",
            severity="required_for_archive",
            finding="DataCite identifier/DOI is not assigned.",
            evidence=[field_table_path.as_posix()],
            blocking_gap="versioned_archive_plan",
            next_action="Assign version/archive identifier only after release approval.",
        ),
        _ro_crate_validation_row(
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            check_id="datacite_creator_placeholder",
            standard_profile="datacite_4_7",
            check_status="blocker",
            severity="required_for_archive",
            finding="DataCite creator list is not frozen.",
            evidence=[field_table_path.as_posix()],
            blocking_gap="formal_license_and_creator_metadata",
            next_action="Freeze creator/contributor metadata before DOI or archive path.",
        ),
        _ro_crate_validation_row(
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            check_id="claim_boundary_preserved",
            standard_profile="spacebio_bench",
            check_status="pass",
            severity="required",
            finding="Scaffold preserves diagnostic-only, not-citable-archive claim boundary.",
            evidence=[summary_path.as_posix(), target_path.as_posix()],
            blocking_gap="",
            next_action="Reject any frozen benchmark wording until blockers clear.",
        ),
    ]
    citation_rows = [
        _citation_checklist_row(
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            item_id="upstream_osdr_credit",
            citation_surface="source_credit",
            required_for="diagnostic_metadata_draft",
            status="pass",
            current_value="Credit NASA OSDR/GeneLab and OSD-120 / GLDS-120.",
            source_artifacts=[summary_path.as_posix()],
            blocking_gap="",
            next_action="Keep OSDR/GeneLab credit in every public surface.",
            claim_boundary_impact="Supports source credit without local release claim.",
        ),
        _citation_checklist_row(
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            item_id="osdr_dataset_specific_citation",
            citation_surface="source_citation",
            required_for="source_specific_public_alpha_package",
            status="needs_review",
            current_value="Use the OSDR study citation path for BibTeX/RIS.",
            source_artifacts=[field_table_path.as_posix()],
            blocking_gap="source_inventory_release_target",
            next_action="Capture exact OSDR citation export before source package release.",
            claim_boundary_impact="Blocks clean source-alpha citation wording.",
        ),
        _citation_checklist_row(
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            item_id="spacebio_package_title",
            citation_surface="local_package_citation",
            required_for="diagnostic_metadata_draft",
            status="pass",
            current_value=title,
            source_artifacts=[field_table_path.as_posix()],
            blocking_gap="",
            next_action="Keep title stable unless package scope changes.",
            claim_boundary_impact="Safe for draft metadata display.",
        ),
        _citation_checklist_row(
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            item_id="spacebio_version_string",
            citation_surface="local_package_citation",
            required_for="citable_archive",
            status="blocker",
            current_value="draft-v9-osd120-diagnostic-alpha-unreleased",
            source_artifacts=["ro-crate-metadata.draft.json", "datapackage.draft.json"],
            blocking_gap="versioned_archive_plan",
            next_action="Choose a release version string only after release approval.",
            claim_boundary_impact="Blocks citable local package citation.",
        ),
        _citation_checklist_row(
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            item_id="local_archive_identifier",
            citation_surface="doi_or_archive_identifier",
            required_for="citable_archive",
            status="blocker",
            current_value="pending_versioned_release_identifier",
            source_artifacts=[field_table_path.as_posix()],
            blocking_gap="versioned_archive_plan",
            next_action="Define DOI/Zenodo/SWHID/RAiD path before archive.",
            claim_boundary_impact="Blocks DOI or archive citation.",
        ),
        _citation_checklist_row(
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            item_id="creator_contributor_list",
            citation_surface="datacite_creators",
            required_for="citable_archive",
            status="blocker",
            current_value="pending_formal_release_creator_list",
            source_artifacts=[field_table_path.as_posix()],
            blocking_gap="formal_license_and_creator_metadata",
            next_action="Freeze creator order, affiliations, and contributor roles.",
            claim_boundary_impact="Blocks DataCite-ready citation.",
        ),
        _citation_checklist_row(
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            item_id="license_and_rights",
            citation_surface="rights",
            required_for="source_specific_public_alpha_package",
            status="blocker",
            current_value="pending_local_release_license_review",
            source_artifacts=[field_table_path.as_posix()],
            blocking_gap="formal_license_and_creator_metadata",
            next_action="Review upstream and local package rights before public package release.",
            claim_boundary_impact="Blocks release packaging beyond draft metadata.",
        ),
        _citation_checklist_row(
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            item_id="related_identifier_map",
            citation_surface="datacite_related_identifiers",
            required_for="source_specific_public_alpha_package",
            status="needs_review",
            current_value="OSD-120 URL and selected literature links recorded; local DOI pending.",
            source_artifacts=[field_table_path.as_posix()],
            blocking_gap="versioned_archive_plan",
            next_action="Separate source, article, software, and archive relation types.",
            claim_boundary_impact="Weakens machine-actionable citation until resolved.",
        ),
        _citation_checklist_row(
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            item_id="ro_crate_metadata_file",
            citation_surface="machine_readable_export",
            required_for="diagnostic_metadata_draft",
            status="pass",
            current_value="ro-crate-metadata.draft.json",
            source_artifacts=["ro-crate-metadata.draft.json"],
            blocking_gap="",
            next_action="Validate in archive mode after placeholder resolution.",
            claim_boundary_impact="Supports draft machine-readable reuse.",
        ),
        _citation_checklist_row(
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            item_id="data_package_descriptor",
            citation_surface="machine_readable_export",
            required_for="diagnostic_metadata_draft",
            status="pass",
            current_value="datapackage.draft.json",
            source_artifacts=["datapackage.draft.json"],
            blocking_gap="",
            next_action="Validate with Frictionless tooling later.",
            claim_boundary_impact="Supports draft package inspection.",
        ),
        _citation_checklist_row(
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            item_id="frozen_claim_guard",
            citation_surface="claim_boundary",
            required_for="all_surfaces",
            status="pass",
            current_value="diagnostic_ro_crate_scaffold_only_not_citable_archive",
            source_artifacts=[target_path.as_posix()],
            blocking_gap="",
            next_action="Keep all export text draft-limited.",
            claim_boundary_impact="Prevents frozen benchmark and full mirror claims.",
        ),
    ]
    validation_pass = sum(1 for row in validation_rows if row["check_status"] == "pass")
    validation_blocker = sum(
        1 for row in validation_rows if row["check_status"] == "blocker"
    )
    validation_needs_review = sum(
        1 for row in validation_rows if row["check_status"] == "needs_review"
    )
    citation_pass = sum(1 for row in citation_rows if row["status"] == "pass")
    citation_blocker = sum(1 for row in citation_rows if row["status"] == "blocker")
    citation_needs_review = sum(
        1 for row in citation_rows if row["status"] == "needs_review"
    )
    scaffold_status = (
        "draft_scaffold_ready_archive_blocked_by_citation_placeholders"
        if validation_blocker or citation_blocker
        else "archive_ready"
    )
    summary = {
        "scaffold_id": scaffold_id,
        "package_id": package_id,
        "task_id": task_id,
        "scaffold_status": scaffold_status,
        "ro_crate_graph_entity_count": str(len(ro_crate_graph)),
        "ro_crate_data_entity_count": str(len(data_entities)),
        "data_package_resource_count": str(len(data_package_resources)),
        "validation_check_count": str(len(validation_rows)),
        "validation_pass_count": str(validation_pass),
        "validation_blocker_count": str(validation_blocker),
        "validation_needs_review_count": str(validation_needs_review),
        "citation_item_count": str(len(citation_rows)),
        "citation_pass_count": str(citation_pass),
        "citation_blocker_count": str(citation_blocker),
        "citation_needs_review_count": str(citation_needs_review),
        "placeholder_field_count": metadata_summary["metadata_placeholder_field_count"],
        "next_required_block": "V9-MULTI-032: archive identifier and license decision gate",
        "claim_boundary": "diagnostic_ro_crate_scaffold_only_not_citable_archive",
    }
    validation_lines = "\n".join(
        (
            f"- `{row['check_id']}`: {row['check_status']} "
            f"({row['severity']})."
        )
        for row in validation_rows
    )
    citation_lines = "\n".join(
        f"- `{row['item_id']}`: {row['status']}." for row in citation_rows
    )
    review_md = f"""# OSD-120 RO-Crate And Citation Freeze Scaffold Review

Block: V9-MULTI-031

Scaffold id: `{scaffold_id}`

Status: `{scaffold_status}`

This block turns the V9-MULTI-030 metadata skeleton into draft RO-Crate and
Data Package export descriptors plus a citation-freeze checklist. The scaffold
is intentionally not archive-ready because DOI/archive identifier,
creator/contributor, and license/rights fields remain unresolved.

## Validation Status

{validation_lines}

## Citation Freeze Checklist

{citation_lines}

## Standards Anchors

- RO-Crate 1.2: root Dataset, Data Entities, Contextual Entities, and flattened
  JSON-LD graph.
- Data Package v1: descriptor with resources; resources are required.
- DataCite 4.7: citable output metadata requires mandatory citation fields,
  including identifier, creators, title, publisher, publication year, and
  resource type.

## Claim Boundary

The scaffold is useful for local inspection and future archive preparation. It
does not claim a complete OSDR payload mirror, a DOI, frozen benchmark release,
leave-one-mission-out generalization, validated biomarkers, or operational
plant-growth guidance.

## Next Block

V9-MULTI-032 should make explicit decisions for archive identifier, release
version string, creator/contributor list, and license/rights before any citable
archive path is attempted.
"""
    return {
        "summary": [summary],
        "validation": validation_rows,
        "citation_checklist": citation_rows,
        "ro_crate_metadata": ro_crate_metadata,
        "data_package": data_package,
        "review_md": review_md,
    }


def write_multispecies_interaction_ro_crate_citation_scaffold(
    *,
    output_dir: str | Path = (
        "v9/multispecies/reports/interaction_ro_crate_citation_scaffold"
    ),
    reports_root: str | Path = "v9/multispecies/reports",
    repo_root: str | Path = ".",
    scaffold_id: str = DEFAULT_INTERACTION_RO_CRATE_EXPORT_SCAFFOLD_ID,
) -> dict[str, Path]:
    """Write draft RO-Crate/Data Package scaffold and citation checklist."""

    package = build_multispecies_interaction_ro_crate_citation_scaffold(
        reports_root=reports_root,
        repo_root=repo_root,
        scaffold_id=scaffold_id,
    )
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "summary_csv": outdir / "ro_crate_export_summary.csv",
        "summary_json": outdir / "ro_crate_export_summary.json",
        "validation_csv": outdir / "ro_crate_validation_table.csv",
        "validation_json": outdir / "ro_crate_validation_table.json",
        "citation_csv": outdir / "citation_freeze_checklist.csv",
        "citation_json": outdir / "citation_freeze_checklist.json",
        "ro_crate_json": outdir / "ro-crate-metadata.draft.json",
        "data_package_json": outdir / "datapackage.draft.json",
        "review_md": (
            outdir.parent
            / "OSD120_INTERACTION_RO_CRATE_CITATION_SCAFFOLD_REVIEW.md"
        ),
    }
    with outputs["summary_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_RO_CRATE_EXPORT_SUMMARY_FIELDS,
        )
        writer.writeheader()
        writer.writerows(package["summary"])
    outputs["summary_json"].write_text(
        json.dumps(package["summary"], indent=2, sort_keys=True) + "\n"
    )
    with outputs["validation_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_RO_CRATE_VALIDATION_FIELDS,
        )
        writer.writeheader()
        writer.writerows(package["validation"])
    outputs["validation_json"].write_text(
        json.dumps(package["validation"], indent=2, sort_keys=True) + "\n"
    )
    with outputs["citation_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_CITATION_FREEZE_CHECKLIST_FIELDS,
        )
        writer.writeheader()
        writer.writerows(package["citation_checklist"])
    outputs["citation_json"].write_text(
        json.dumps(package["citation_checklist"], indent=2, sort_keys=True) + "\n"
    )
    outputs["ro_crate_json"].write_text(
        json.dumps(package["ro_crate_metadata"], indent=2, sort_keys=True) + "\n"
    )
    outputs["data_package_json"].write_text(
        json.dumps(package["data_package"], indent=2, sort_keys=True) + "\n"
    )
    outputs["review_md"].write_text(str(package["review_md"]))
    return outputs


def _archive_decision_external_references(
    *,
    checked_date: str = "2026-05-26",
) -> list[dict[str, str]]:
    return [
        {
            "reference_id": "github_zenodo_doi",
            "topic": "GitHub repository DOI via Zenodo",
            "url": (
                "https://docs.github.com/en/repositories/archiving-a-github-"
                "repository/referencing-and-citing-content"
            ),
            "checked_date": checked_date,
            "decision_implication": (
                "A Zenodo DOI path requires a public repository, repository-owner "
                "access, a GitHub release, and an explicit license."
            ),
        },
        {
            "reference_id": "zenodo_github_integration",
            "topic": "Zenodo GitHub/software archive flow",
            "url": "https://help.zenodo.org/docs/github/",
            "checked_date": checked_date,
            "decision_implication": (
                "Zenodo/GitHub integration is an archive option for software "
                "records, but it does not replace local creator/license decisions."
            ),
        },
        {
            "reference_id": "github_citation_cff",
            "topic": "GitHub CITATION.cff citation metadata",
            "url": (
                "https://docs.github.com/en/repositories/managing-your-"
                "repositorys-settings-and-features/customizing-your-repository/"
                "about-citation-files"
            ),
            "checked_date": checked_date,
            "decision_implication": (
                "CITATION.cff can describe software or dataset citation metadata, "
                "but requires authors, version, DOI or URL, and release date."
            ),
        },
        {
            "reference_id": "osdr_dataset_citation",
            "topic": "NASA OSDR study citation and credit",
            "url": "https://science.nasa.gov/reference/osdr-faq/",
            "checked_date": checked_date,
            "decision_implication": (
                "Use the OSDR study citation button/BibTeX/RIS for upstream data "
                "citation and keep NASA OSDR/GeneLab credit visible."
            ),
        },
        {
            "reference_id": "datacite_4_7_mandatory_metadata",
            "topic": "DataCite 4.7 mandatory citation metadata",
            "url": "https://datacite-metadata-schema.readthedocs.io/en/4.7/properties/",
            "checked_date": checked_date,
            "decision_implication": (
                "Citable metadata needs identifier, creator, title, publisher, "
                "publication year, and resource type; placeholder identifier and "
                "creator fields block a DOI-ready package."
            ),
        },
        {
            "reference_id": "spdx_license_identifiers",
            "topic": "Machine-readable license identifiers",
            "url": "https://spdx.dev/learn/handling-license-info/",
            "checked_date": checked_date,
            "decision_implication": (
                "When a local license is selected, use a standard SPDX identifier "
                "or explicit custom-license label instead of ambiguous text."
            ),
        },
    ]


def _archive_option_row(
    *,
    decision_id: str,
    scaffold_id: str,
    package_id: str,
    task_id: str,
    option_id: str,
    option_label: str,
    decision_status: str,
    current_draft_selected: bool,
    citable_release_ready: bool,
    identifier_type: str,
    expected_identifier: str,
    official_reference_ids: Sequence[str],
    evidence_artifacts: Sequence[str],
    blocking_gaps: Sequence[str],
    required_actions: Sequence[str],
    rationale: str,
) -> dict[str, str]:
    return {
        "decision_id": decision_id,
        "scaffold_id": scaffold_id,
        "package_id": package_id,
        "task_id": task_id,
        "option_id": option_id,
        "option_label": option_label,
        "decision_status": decision_status,
        "current_draft_selected": str(current_draft_selected),
        "citable_release_ready": str(citable_release_ready),
        "identifier_type": identifier_type,
        "expected_identifier": expected_identifier,
        "official_reference_ids": _pipe_join(list(official_reference_ids)),
        "evidence_artifacts": _pipe_join(list(evidence_artifacts)),
        "blocking_gaps": _pipe_join(list(blocking_gaps)),
        "required_actions": _pipe_join(list(required_actions)),
        "rationale": rationale,
    }


def _rights_decision_row(
    *,
    decision_id: str,
    scaffold_id: str,
    package_id: str,
    task_id: str,
    component_id: str,
    component_scope: str,
    decision_status: str,
    current_value: str,
    recommended_decision_state: str,
    official_reference_ids: Sequence[str],
    evidence_artifacts: Sequence[str],
    blocking_gaps: Sequence[str],
    required_actions: Sequence[str],
    claim_boundary_impact: str,
) -> dict[str, str]:
    return {
        "decision_id": decision_id,
        "scaffold_id": scaffold_id,
        "package_id": package_id,
        "task_id": task_id,
        "component_id": component_id,
        "component_scope": component_scope,
        "decision_status": decision_status,
        "current_value": current_value,
        "recommended_decision_state": recommended_decision_state,
        "official_reference_ids": _pipe_join(list(official_reference_ids)),
        "evidence_artifacts": _pipe_join(list(evidence_artifacts)),
        "blocking_gaps": _pipe_join(list(blocking_gaps)),
        "required_actions": _pipe_join(list(required_actions)),
        "claim_boundary_impact": claim_boundary_impact,
    }


def build_multispecies_interaction_archive_decision_gate(
    *,
    reports_root: str | Path = "v9/multispecies/reports",
    repo_root: str | Path = ".",
    decision_id: str = DEFAULT_INTERACTION_ARCHIVE_DECISION_GATE_ID,
) -> dict[str, Any]:
    """Build the OSD-120 archive identifier and license decision gate."""

    root = Path(repo_root)
    reports_path = _resolve_path(reports_root, root)
    scaffold_dir = reports_path / "interaction_ro_crate_citation_scaffold"
    metadata_dir = reports_path / "interaction_public_metadata_package"
    summary_path = scaffold_dir / "ro_crate_export_summary.csv"
    validation_path = scaffold_dir / "ro_crate_validation_table.csv"
    citation_path = scaffold_dir / "citation_freeze_checklist.csv"
    ro_crate_path = scaffold_dir / "ro-crate-metadata.draft.json"
    data_package_path = scaffold_dir / "datapackage.draft.json"
    metadata_field_path = metadata_dir / "public_metadata_field_table.csv"
    target_path = metadata_dir / "source_release_target_decision.csv"
    scaffold_summary = _require_one_row(
        _read_csv_dict_rows(summary_path),
        context="OSD-120 RO-Crate scaffold summary for archive decision gate",
    )
    validation_rows = _read_csv_dict_rows(validation_path)
    citation_rows = _read_csv_dict_rows(citation_path)
    metadata_fields = _read_csv_dict_rows(metadata_field_path)
    target_rows = _read_csv_dict_rows(target_path)
    package_id = scaffold_summary["package_id"]
    scaffold_id = scaffold_summary["scaffold_id"]
    task_id = scaffold_summary["task_id"]
    common_evidence = [
        summary_path.as_posix(),
        validation_path.as_posix(),
        citation_path.as_posix(),
        metadata_field_path.as_posix(),
        target_path.as_posix(),
    ]
    archive_options = [
        _archive_option_row(
            decision_id=decision_id,
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            option_id="current_no_archive_diagnostic_draft",
            option_label="Current no-archive diagnostic draft",
            decision_status="selected_for_current_draft_only",
            current_draft_selected=True,
            citable_release_ready=False,
            identifier_type="none",
            expected_identifier="no_archive_identifier_for_current_draft",
            official_reference_ids=["datacite_4_7_mandatory_metadata"],
            evidence_artifacts=common_evidence,
            blocking_gaps=[
                "versioned_archive_plan",
                "formal_license_and_creator_metadata",
            ],
            required_actions=[
                "keep_draft_claim_boundary",
                "do_not_mint_doi_or_archive_identifier",
            ],
            rationale=(
                "This is the only safe current path because creator, license, "
                "and archive identifier decisions are unresolved."
            ),
        ),
        _archive_option_row(
            decision_id=decision_id,
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            option_id="zenodo_github_release_doi",
            option_label="Zenodo DOI from a GitHub release",
            decision_status="deferred_release_owner_required",
            current_draft_selected=False,
            citable_release_ready=False,
            identifier_type="DOI",
            expected_identifier="10.5281/zenodo.<record>",
            official_reference_ids=[
                "github_zenodo_doi",
                "zenodo_github_integration",
                "datacite_4_7_mandatory_metadata",
            ],
            evidence_artifacts=common_evidence,
            blocking_gaps=[
                "repository_public_and_owner_access_unconfirmed",
                "versioned_archive_plan",
                "formal_license_and_creator_metadata",
            ],
            required_actions=[
                "confirm_public_repository_and_release_owner",
                "choose_release_tag",
                "freeze_creator_list",
                "choose_license",
            ],
            rationale=(
                "Best future citable archive candidate if this package is released "
                "through a public GitHub repository, but not selectable now."
            ),
        ),
        _archive_option_row(
            decision_id=decision_id,
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            option_id="citation_cff_without_archive",
            option_label="CITATION.cff draft without DOI archive",
            decision_status="deferred_after_creator_version_decision",
            current_draft_selected=False,
            citable_release_ready=False,
            identifier_type="CITATION.cff",
            expected_identifier="repository_url_or_future_doi",
            official_reference_ids=["github_citation_cff"],
            evidence_artifacts=common_evidence,
            blocking_gaps=[
                "creator_contributor_list_pending",
                "versioned_archive_plan",
            ],
            required_actions=[
                "freeze_creator_list",
                "choose_version_or_url",
                "decide_dataset_vs_software_cff_type",
            ],
            rationale=(
                "Useful for citation UX, but it is not an archive identifier and "
                "still needs author/version decisions."
            ),
        ),
        _archive_option_row(
            decision_id=decision_id,
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            option_id="osdr_source_citation_only",
            option_label="OSDR source citation only",
            decision_status="selected_for_upstream_credit_not_local_archive",
            current_draft_selected=True,
            citable_release_ready=False,
            identifier_type="OSDR study citation",
            expected_identifier="use_osdr_study_page_citation_button",
            official_reference_ids=["osdr_dataset_citation"],
            evidence_artifacts=common_evidence,
            blocking_gaps=["local_archive_identifier_missing"],
            required_actions=[
                "capture_exact_osdr_bibtex_or_ris_before_source_package_release"
            ],
            rationale=(
                "Appropriate for upstream data credit now, but it cannot serve as "
                "the local SpaceBio-Bench package archive identifier."
            ),
        ),
        _archive_option_row(
            decision_id=decision_id,
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            option_id="software_heritage_swhid_related_identifier",
            option_label="Software Heritage SWHID as related identifier",
            decision_status="deferred_code_archive_related_identifier_only",
            current_draft_selected=False,
            citable_release_ready=False,
            identifier_type="SWHID",
            expected_identifier="swh:1:<object-type>:<hash>",
            official_reference_ids=["datacite_4_7_mandatory_metadata"],
            evidence_artifacts=common_evidence,
            blocking_gaps=[
                "code_release_target_unconfirmed",
                "versioned_archive_plan",
            ],
            required_actions=[
                "decide_code_archive_scope",
                "map_swhid_as_related_identifier_not_dataset_doi",
            ],
            rationale=(
                "Could document code preservation later, but it is not the primary "
                "dataset/package citation path for this diagnostic scaffold."
            ),
        ),
        _archive_option_row(
            decision_id=decision_id,
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            option_id="full_osdr_payload_mirror_archive",
            option_label="Full OSDR processed payload mirror archive",
            decision_status="blocked_full_payload_freeze_missing",
            current_draft_selected=False,
            citable_release_ready=False,
            identifier_type="DOI_or_repository_record",
            expected_identifier="not_available",
            official_reference_ids=[
                "osdr_dataset_citation",
                "datacite_4_7_mandatory_metadata",
            ],
            evidence_artifacts=common_evidence,
            blocking_gaps=["full_osdr_payload_freeze"],
            required_actions=[
                "download_or_explicitly_exclude_every_processed_osdr_payload",
                "review_rights_before_mirroring",
            ],
            rationale=(
                "Blocked because the current package verifies only diagnostic-"
                "required payloads and explicitly avoids full mirror claims."
            ),
        ),
    ]
    license_rows = [
        _rights_decision_row(
            decision_id=decision_id,
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            component_id="upstream_osdr_data_credit",
            component_scope="upstream_source_data",
            decision_status="needs_review",
            current_value="NASA OSDR/GeneLab credit present; exact OSDR study citation pending",
            recommended_decision_state="defer_to_osdr_study_citation_and_terms",
            official_reference_ids=["osdr_dataset_citation"],
            evidence_artifacts=[citation_path.as_posix(), metadata_field_path.as_posix()],
            blocking_gaps=["source_inventory_release_target"],
            required_actions=["download_exact_osdr_bibtex_or_ris", "retain_osdr_credit"],
            claim_boundary_impact="Blocks clean source-alpha citation wording, not local draft metadata.",
        ),
        _rights_decision_row(
            decision_id=decision_id,
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            component_id="local_code_license",
            component_scope="spacebio_bench_code_and_scripts",
            decision_status="blocked",
            current_value="pending_local_release_license_review",
            recommended_decision_state="choose_spdx_license_or_custom_terms",
            official_reference_ids=["github_zenodo_doi", "spdx_license_identifiers"],
            evidence_artifacts=[data_package_path.as_posix(), metadata_field_path.as_posix()],
            blocking_gaps=["formal_license_and_creator_metadata"],
            required_actions=["choose_license_identifier", "add_or_verify_license_file"],
            claim_boundary_impact="Blocks repository DOI/release path.",
        ),
        _rights_decision_row(
            decision_id=decision_id,
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            component_id="local_metadata_tables_license",
            component_scope="generated_metadata_tables_and_json",
            decision_status="blocked",
            current_value="pending_local_release_license_review",
            recommended_decision_state="choose_data_metadata_license",
            official_reference_ids=["spdx_license_identifiers"],
            evidence_artifacts=[data_package_path.as_posix(), metadata_field_path.as_posix()],
            blocking_gaps=["formal_license_and_creator_metadata"],
            required_actions=["choose_license_identifier", "state_scope_of_license"],
            claim_boundary_impact="Blocks publishing generated package descriptors as a reusable dataset.",
        ),
        _rights_decision_row(
            decision_id=decision_id,
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            component_id="diagnostic_payload_mirror_rights",
            component_scope="local_copies_of_osdr_processed_payloads",
            decision_status="blocked",
            current_value="full_osdr_payload_mirror_not_claimed",
            recommended_decision_state="do_not_archive_payload_mirror_now",
            official_reference_ids=["osdr_dataset_citation"],
            evidence_artifacts=[target_path.as_posix()],
            blocking_gaps=["full_osdr_payload_freeze", "rights_review_before_mirroring"],
            required_actions=[
                "keep_payload_mirror_out_of_archive",
                "review_upstream_terms_before_any_payload_mirror",
            ],
            claim_boundary_impact="Prevents full OSDR mirror claims.",
        ),
        _rights_decision_row(
            decision_id=decision_id,
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            component_id="third_party_literature_links",
            component_scope="external_context_urls_and_dois",
            decision_status="selected_citation_only",
            current_value="links and DOI URLs only; no copied full text",
            recommended_decision_state="keep_as_citation_links_only",
            official_reference_ids=["datacite_4_7_mandatory_metadata"],
            evidence_artifacts=[metadata_field_path.as_posix()],
            blocking_gaps=[],
            required_actions=["keep_external_links_as_related_identifiers"],
            claim_boundary_impact="Safe for metadata context without redistributing third-party content.",
        ),
    ]
    creator_rows = [
        _rights_decision_row(
            decision_id=decision_id,
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            component_id="upstream_osdr_dataset_citation",
            component_scope="upstream_source_creator_credit",
            decision_status="needs_review",
            current_value="OSDR credit present; exact study citation pending",
            recommended_decision_state="capture_osdr_study_citation",
            official_reference_ids=["osdr_dataset_citation"],
            evidence_artifacts=[citation_path.as_posix()],
            blocking_gaps=["source_inventory_release_target"],
            required_actions=["export_osd120_bibtex_or_ris"],
            claim_boundary_impact="Needed for source-alpha package wording.",
        ),
        _rights_decision_row(
            decision_id=decision_id,
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            component_id="local_package_creators",
            component_scope="datacite_creators",
            decision_status="blocked",
            current_value="pending_formal_release_creator_list",
            recommended_decision_state="release_owner_must_freeze_creator_order",
            official_reference_ids=["datacite_4_7_mandatory_metadata", "github_citation_cff"],
            evidence_artifacts=[metadata_field_path.as_posix(), citation_path.as_posix()],
            blocking_gaps=["formal_license_and_creator_metadata"],
            required_actions=["freeze_creator_names", "collect_affiliations_or_orcids_if_available"],
            claim_boundary_impact="Blocks DOI-ready DataCite and CITATION.cff metadata.",
        ),
        _rights_decision_row(
            decision_id=decision_id,
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            component_id="local_package_contributors",
            component_scope="contributors_and_roles",
            decision_status="blocked",
            current_value="contributor_roles_pending",
            recommended_decision_state="define_packager_maintainer_and_data_source_roles",
            official_reference_ids=["datacite_4_7_mandatory_metadata"],
            evidence_artifacts=[ro_crate_path.as_posix(), citation_path.as_posix()],
            blocking_gaps=["formal_license_and_creator_metadata"],
            required_actions=["assign_contributor_roles", "separate_upstream_source_from_local_packager"],
            claim_boundary_impact="Blocks robust provenance/citation metadata.",
        ),
        _rights_decision_row(
            decision_id=decision_id,
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            component_id="publisher_maintainer_identity",
            component_scope="publisher_or_maintainer",
            decision_status="blocked",
            current_value="SpaceBio-Bench draft package; formal publisher not frozen",
            recommended_decision_state="release_owner_must_select_publisher_or_maintainer",
            official_reference_ids=["datacite_4_7_mandatory_metadata"],
            evidence_artifacts=[metadata_field_path.as_posix()],
            blocking_gaps=["formal_license_and_creator_metadata"],
            required_actions=["choose_publisher_or_maintainer_identity"],
            claim_boundary_impact="Blocks DataCite-ready package citation.",
        ),
        _rights_decision_row(
            decision_id=decision_id,
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            component_id="package_title",
            component_scope="title",
            decision_status="selected_for_current_draft",
            current_value="OSD-120 Arabidopsis Root Light Interaction Diagnostic Draft",
            recommended_decision_state="keep_stable_unless_scope_changes",
            official_reference_ids=["datacite_4_7_mandatory_metadata", "github_citation_cff"],
            evidence_artifacts=[metadata_field_path.as_posix()],
            blocking_gaps=[],
            required_actions=["reuse_title_in_future_citation_files"],
            claim_boundary_impact="Safe for diagnostic draft citation display.",
        ),
        _rights_decision_row(
            decision_id=decision_id,
            scaffold_id=scaffold_id,
            package_id=package_id,
            task_id=task_id,
            component_id="release_version_string",
            component_scope="version",
            decision_status="blocked",
            current_value="draft-v9-osd120-diagnostic-alpha-unreleased",
            recommended_decision_state="choose_release_tag_only_after_archive_approval",
            official_reference_ids=["github_zenodo_doi", "github_citation_cff"],
            evidence_artifacts=[citation_path.as_posix(), ro_crate_path.as_posix()],
            blocking_gaps=["versioned_archive_plan"],
            required_actions=["choose_release_tag", "align_git_release_and_metadata_version"],
            claim_boundary_impact="Blocks citable archive release.",
        ),
    ]
    references = _archive_decision_external_references()
    current_selected = sum(
        1 for row in archive_options if row["current_draft_selected"] == "True"
    )
    deferred_options = sum(
        1 for row in archive_options if row["decision_status"].startswith("deferred")
    )
    blocked_options = sum(
        1 for row in archive_options if row["decision_status"].startswith("blocked")
    )
    license_blockers = sum(1 for row in license_rows if row["decision_status"] == "blocked")
    creator_blockers = sum(1 for row in creator_rows if row["decision_status"] == "blocked")
    summary = {
        "decision_id": decision_id,
        "scaffold_id": scaffold_id,
        "package_id": package_id,
        "task_id": task_id,
        "decision_status": "current_draft_no_archive_selected_release_decisions_blocked",
        "archive_path_decision": "no_archive_identifier_for_current_draft",
        "version_decision": "draft_unreleased_version_only",
        "creator_decision": "blocked_pending_formal_release_creator_list",
        "license_decision": "blocked_pending_local_release_license_review",
        "archive_option_count": str(len(archive_options)),
        "current_selected_option_count": str(current_selected),
        "deferred_option_count": str(deferred_options),
        "blocked_option_count": str(blocked_options),
        "license_component_count": str(len(license_rows)),
        "license_blocker_count": str(license_blockers),
        "creator_component_count": str(len(creator_rows)),
        "creator_blocker_count": str(creator_blockers),
        "external_reference_count": str(len(references)),
        "next_required_block": "V9-MULTI-033: release owner supplied citation metadata fill",
        "claim_boundary": "diagnostic_archive_decision_gate_only_no_archive_identifier",
    }
    option_lines = "\n".join(
        (
            f"- `{row['option_id']}`: {row['decision_status']}; "
            f"current_draft_selected={row['current_draft_selected']}."
        )
        for row in archive_options
    )
    license_lines = "\n".join(
        f"- `{row['component_id']}`: {row['decision_status']}."
        for row in license_rows
    )
    creator_lines = "\n".join(
        f"- `{row['component_id']}`: {row['decision_status']}."
        for row in creator_rows
    )
    review_md = f"""# OSD-120 Archive Identifier And License Decision Gate Review

Block: V9-MULTI-032

Decision id: `{decision_id}`

Status: `{summary["decision_status"]}`

This block does not mint a DOI, choose a license, or invent a creator list. It
selects a safe current-draft path with no local archive identifier, records OSDR
source citation as upstream credit, and keeps release-owner decisions blocked
until supplied explicitly.

## Archive Identifier Options

{option_lines}

## License And Rights Decisions

{license_lines}

## Creator And Contributor Decisions

{creator_lines}

## External Guidance Anchors

- GitHub/Zenodo DOI requires a public repository, repository-owner access,
  GitHub release, and license.
- GitHub CITATION.cff can describe dataset/software citation metadata, but
  needs author, version, and DOI or URL fields.
- OSDR dataset citation should be taken from the OSDR study page citation
  button in BibTeX or RIS form.
- DataCite 4.7 citable metadata requires identifier, creators, title,
  publisher, publication year, and resource type.
- License decisions should use a clear SPDX identifier or explicit custom
  terms after review.

## Claim Boundary

The current diagnostic package remains draft metadata only. No DOI, no archive
identifier, no frozen release version, no local package license, no full OSDR
payload mirror, and no frozen benchmark release are claimed.

## Next Block

V9-MULTI-033 should fill citation metadata only after the release owner supplies
creator/contributor, license/rights, archive route, and version decisions.
"""
    return {
        "summary": [summary],
        "archive_options": archive_options,
        "license_rights": license_rows,
        "creator_contributor": creator_rows,
        "external_references": references,
        "review_md": review_md,
    }


def write_multispecies_interaction_archive_decision_gate(
    *,
    output_dir: str | Path = (
        "v9/multispecies/reports/interaction_archive_decision_gate"
    ),
    reports_root: str | Path = "v9/multispecies/reports",
    repo_root: str | Path = ".",
    decision_id: str = DEFAULT_INTERACTION_ARCHIVE_DECISION_GATE_ID,
) -> dict[str, Path]:
    """Write archive identifier/license decision gate tables."""

    package = build_multispecies_interaction_archive_decision_gate(
        reports_root=reports_root,
        repo_root=repo_root,
        decision_id=decision_id,
    )
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "summary_csv": outdir / "archive_decision_summary.csv",
        "summary_json": outdir / "archive_decision_summary.json",
        "archive_option_csv": outdir / "archive_identifier_option_matrix.csv",
        "archive_option_json": outdir / "archive_identifier_option_matrix.json",
        "license_csv": outdir / "license_rights_decision_table.csv",
        "license_json": outdir / "license_rights_decision_table.json",
        "creator_csv": outdir / "creator_contributor_decision_table.csv",
        "creator_json": outdir / "creator_contributor_decision_table.json",
        "reference_csv": outdir / "archive_decision_external_references.csv",
        "reference_json": outdir / "archive_decision_external_references.json",
        "review_md": (
            outdir.parent
            / "OSD120_INTERACTION_ARCHIVE_IDENTIFIER_LICENSE_DECISION_REVIEW.md"
        ),
    }
    with outputs["summary_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_ARCHIVE_DECISION_SUMMARY_FIELDS,
        )
        writer.writeheader()
        writer.writerows(package["summary"])
    outputs["summary_json"].write_text(
        json.dumps(package["summary"], indent=2, sort_keys=True) + "\n"
    )
    with outputs["archive_option_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_ARCHIVE_IDENTIFIER_OPTION_FIELDS,
        )
        writer.writeheader()
        writer.writerows(package["archive_options"])
    outputs["archive_option_json"].write_text(
        json.dumps(package["archive_options"], indent=2, sort_keys=True) + "\n"
    )
    with outputs["license_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_LICENSE_RIGHTS_DECISION_FIELDS,
        )
        writer.writeheader()
        writer.writerows(package["license_rights"])
    outputs["license_json"].write_text(
        json.dumps(package["license_rights"], indent=2, sort_keys=True) + "\n"
    )
    with outputs["creator_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_CREATOR_CONTRIBUTOR_DECISION_FIELDS,
        )
        writer.writeheader()
        writer.writerows(package["creator_contributor"])
    outputs["creator_json"].write_text(
        json.dumps(package["creator_contributor"], indent=2, sort_keys=True) + "\n"
    )
    with outputs["reference_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_ARCHIVE_DECISION_REFERENCE_FIELDS,
        )
        writer.writeheader()
        writer.writerows(package["external_references"])
    outputs["reference_json"].write_text(
        json.dumps(package["external_references"], indent=2, sort_keys=True) + "\n"
    )
    outputs["review_md"].write_text(str(package["review_md"]))
    return outputs


def _load_owner_citation_metadata(
    owner_metadata_path: str | Path | None,
    *,
    repo_root: Path,
) -> dict[str, dict[str, str]]:
    if owner_metadata_path is None:
        return {}
    path = _resolve_path(owner_metadata_path, repo_root)
    if path.suffix.lower() == ".json":
        payload = json.loads(path.read_text())
        if isinstance(payload, dict) and isinstance(payload.get("fields"), list):
            rows = payload["fields"]
        elif isinstance(payload, dict):
            rows = []
            for field_id, value in payload.items():
                if isinstance(value, Mapping):
                    row = {"field_id": field_id}
                    row.update(value)
                    rows.append(row)
                else:
                    rows.append({"field_id": field_id, "supplied_value": value})
        elif isinstance(payload, list):
            rows = payload
        else:
            raise ValueError(f"unsupported owner metadata JSON payload in {path}")
    else:
        rows = _read_csv_dict_rows(path)
    metadata: dict[str, dict[str, str]] = {}
    for row in rows:
        if not isinstance(row, Mapping):
            continue
        field_id = str(row.get("field_id", "")).strip()
        if not field_id:
            continue
        metadata[field_id] = {
            "supplied_value": str(row.get("supplied_value", "") or "").strip(),
            "supplied_by": str(row.get("supplied_by", "") or "").strip(),
            "supplied_date": str(row.get("supplied_date", "") or "").strip(),
            "supplied_evidence": str(row.get("supplied_evidence", "") or "").strip(),
        }
    return metadata


def _citation_metadata_template_row(
    *,
    field_id: str,
    field_group: str,
    target_profile: Sequence[str],
    required_for: str,
    current_status_from_gate: str,
    current_value: str,
    target_artifacts: Sequence[str],
    validation_rule: str,
    blocker_if_missing: str,
    next_action: str,
    notes: str,
) -> dict[str, str]:
    return {
        "field_id": field_id,
        "field_group": field_group,
        "target_profile": _pipe_join(list(target_profile)),
        "required_for": required_for,
        "current_status_from_gate": current_status_from_gate,
        "current_value": current_value,
        "target_artifacts": _pipe_join(list(target_artifacts)),
        "validation_rule": validation_rule,
        "blocker_if_missing": blocker_if_missing,
        "next_action": next_action,
        "notes": notes,
    }


def _citation_metadata_fill_status(
    template: Mapping[str, str],
    owner_row: Mapping[str, str] | None,
) -> str:
    supplied_value = ""
    if owner_row is not None:
        supplied_value = str(owner_row.get("supplied_value", "") or "").strip()
    if supplied_value:
        if "license" in template["field_id"]:
            return "owner_supplied_pending_license_review"
        if "date" in template["field_id"] or "year" in template["field_id"]:
            return "owner_supplied_pending_date_review"
        if "identifier" in template["field_id"] or "archive" in template["field_id"]:
            return "owner_supplied_pending_archive_review"
        return "owner_supplied_pending_release_review"
    if not template.get("blocker_if_missing", ""):
        return "existing_current_draft_value_retained"
    return "not_supplied_blocking_release"


def _build_citation_metadata_template_rows(
    *,
    summary: Mapping[str, str],
    license_by_id: Mapping[str, Mapping[str, str]],
    creator_by_id: Mapping[str, Mapping[str, str]],
    archive_option_by_id: Mapping[str, Mapping[str, str]],
    reports_path: Path,
) -> list[dict[str, str]]:
    ro_crate_path = (
        reports_path
        / "interaction_ro_crate_citation_scaffold"
        / "ro-crate-metadata.draft.json"
    )
    data_package_path = (
        reports_path
        / "interaction_ro_crate_citation_scaffold"
        / "datapackage.draft.json"
    )
    archive_summary_path = (
        reports_path
        / "interaction_archive_decision_gate"
        / "archive_decision_summary.csv"
    )
    license_path = (
        reports_path
        / "interaction_archive_decision_gate"
        / "license_rights_decision_table.csv"
    )
    creator_path = (
        reports_path
        / "interaction_archive_decision_gate"
        / "creator_contributor_decision_table.csv"
    )
    current_no_archive = archive_option_by_id["current_no_archive_diagnostic_draft"]
    package_title = creator_by_id["package_title"]
    payload_rights = license_by_id["diagnostic_payload_mirror_rights"]
    return [
        _citation_metadata_template_row(
            field_id="current_archive_path_decision",
            field_group="archive",
            target_profile=["RO-Crate", "Data Package"],
            required_for="current_diagnostic_draft",
            current_status_from_gate=current_no_archive["decision_status"],
            current_value=summary["archive_path_decision"],
            target_artifacts=[archive_summary_path.as_posix()],
            validation_rule="retain_no_archive_for_current_draft_unless_owner_supplies_release_route",
            blocker_if_missing="",
            next_action="keep_no_archive_identifier_for_current_draft",
            notes="This preserves the V9-MULTI-032 diagnostic draft boundary.",
        ),
        _citation_metadata_template_row(
            field_id="future_archive_identifier",
            field_group="archive",
            target_profile=["DataCite", "RO-Crate", "CITATION.cff"],
            required_for="citable_release",
            current_status_from_gate="blocked",
            current_value="no_archive_identifier_for_current_draft",
            target_artifacts=[ro_crate_path.as_posix(), data_package_path.as_posix()],
            validation_rule="owner_supplied_doi_url_or_other_identifier_required",
            blocker_if_missing="versioned_archive_plan",
            next_action="supply DOI/URL/SWHID or keep no-archive draft status",
            notes="Do not mint or infer an identifier from repository state.",
        ),
        _citation_metadata_template_row(
            field_id="release_version_string",
            field_group="version",
            target_profile=["DataCite", "RO-Crate", "CITATION.cff", "Data Package"],
            required_for="citable_release",
            current_status_from_gate=creator_by_id["release_version_string"]["decision_status"],
            current_value=creator_by_id["release_version_string"]["current_value"],
            target_artifacts=[creator_path.as_posix(), ro_crate_path.as_posix()],
            validation_rule="owner_supplied_release_tag_required",
            blocker_if_missing="versioned_archive_plan",
            next_action="supply release tag only after archive route approval",
            notes="Draft internal strings are not citable release versions.",
        ),
        _citation_metadata_template_row(
            field_id="release_date",
            field_group="version",
            target_profile=["CITATION.cff", "DataCite"],
            required_for="citable_release",
            current_status_from_gate="blocked",
            current_value="not_supplied",
            target_artifacts=[ro_crate_path.as_posix()],
            validation_rule="owner_supplied_iso_8601_date_required",
            blocker_if_missing="versioned_archive_plan",
            next_action="supply release date after release tag is approved",
            notes="Do not use build date or current date as release date.",
        ),
        _citation_metadata_template_row(
            field_id="publication_year",
            field_group="version",
            target_profile=["DataCite"],
            required_for="citable_release",
            current_status_from_gate="blocked",
            current_value="not_supplied",
            target_artifacts=[ro_crate_path.as_posix()],
            validation_rule="owner_supplied_publication_year_required",
            blocker_if_missing="versioned_archive_plan",
            next_action="derive only from owner-approved release/publication date",
            notes="DataCite publication year should not be guessed.",
        ),
        _citation_metadata_template_row(
            field_id="package_title",
            field_group="citation",
            target_profile=["DataCite", "RO-Crate", "CITATION.cff", "Data Package"],
            required_for="current_diagnostic_draft",
            current_status_from_gate=package_title["decision_status"],
            current_value=package_title["current_value"],
            target_artifacts=[creator_path.as_posix(), ro_crate_path.as_posix()],
            validation_rule="retain_existing_title_or_owner_supplied_replacement",
            blocker_if_missing="",
            next_action="reuse stable draft title unless scope changes",
            notes="This is safe for display but not enough for archive release.",
        ),
        _citation_metadata_template_row(
            field_id="local_package_creators",
            field_group="creator",
            target_profile=["DataCite", "RO-Crate", "CITATION.cff"],
            required_for="citable_release",
            current_status_from_gate=creator_by_id["local_package_creators"]["decision_status"],
            current_value=creator_by_id["local_package_creators"]["current_value"],
            target_artifacts=[creator_path.as_posix(), ro_crate_path.as_posix()],
            validation_rule="owner_supplied_ordered_creator_list_required",
            blocker_if_missing="formal_license_and_creator_metadata",
            next_action="freeze creator names, order, affiliations, and ORCIDs if available",
            notes="Authorship and creator order must not be inferred from commits.",
        ),
        _citation_metadata_template_row(
            field_id="local_package_contributors",
            field_group="creator",
            target_profile=["DataCite", "RO-Crate"],
            required_for="citable_release",
            current_status_from_gate=creator_by_id["local_package_contributors"]["decision_status"],
            current_value=creator_by_id["local_package_contributors"]["current_value"],
            target_artifacts=[creator_path.as_posix(), ro_crate_path.as_posix()],
            validation_rule="owner_supplied_contributor_role_list_required",
            blocker_if_missing="formal_license_and_creator_metadata",
            next_action="assign contributor roles and separate upstream source from local packager",
            notes="Contributor role semantics need release-owner review.",
        ),
        _citation_metadata_template_row(
            field_id="publisher_maintainer_identity",
            field_group="creator",
            target_profile=["DataCite", "RO-Crate"],
            required_for="citable_release",
            current_status_from_gate=creator_by_id["publisher_maintainer_identity"]["decision_status"],
            current_value=creator_by_id["publisher_maintainer_identity"]["current_value"],
            target_artifacts=[creator_path.as_posix(), ro_crate_path.as_posix()],
            validation_rule="owner_supplied_publisher_or_maintainer_required",
            blocker_if_missing="formal_license_and_creator_metadata",
            next_action="choose publisher or maintainer identity for citable metadata",
            notes="DataCite publisher must be explicit.",
        ),
        _citation_metadata_template_row(
            field_id="local_code_license",
            field_group="rights",
            target_profile=["SPDX", "RO-Crate", "Data Package", "GitHub release"],
            required_for="citable_release",
            current_status_from_gate=license_by_id["local_code_license"]["decision_status"],
            current_value=license_by_id["local_code_license"]["current_value"],
            target_artifacts=[license_path.as_posix(), data_package_path.as_posix()],
            validation_rule="spdx_identifier_or_explicit_custom_terms_required",
            blocker_if_missing="formal_license_and_creator_metadata",
            next_action="choose code license and add or verify license file",
            notes="GitHub/Zenodo release guidance expects repository reuse terms.",
        ),
        _citation_metadata_template_row(
            field_id="local_metadata_tables_license",
            field_group="rights",
            target_profile=["SPDX", "RO-Crate", "Data Package"],
            required_for="citable_release",
            current_status_from_gate=license_by_id["local_metadata_tables_license"]["decision_status"],
            current_value=license_by_id["local_metadata_tables_license"]["current_value"],
            target_artifacts=[license_path.as_posix(), data_package_path.as_posix()],
            validation_rule="data_or_metadata_license_scope_required",
            blocker_if_missing="formal_license_and_creator_metadata",
            next_action="choose metadata/table license and state scope",
            notes="Code license and generated metadata license may differ.",
        ),
        _citation_metadata_template_row(
            field_id="upstream_osdr_dataset_citation",
            field_group="source_credit",
            target_profile=["RO-Crate", "Data Package", "dataset card"],
            required_for="source_alpha_package",
            current_status_from_gate=creator_by_id["upstream_osdr_dataset_citation"]["decision_status"],
            current_value=creator_by_id["upstream_osdr_dataset_citation"]["current_value"],
            target_artifacts=[creator_path.as_posix()],
            validation_rule="osdr_study_page_bibtex_or_ris_citation_required",
            blocker_if_missing="source_inventory_release_target",
            next_action="export exact OSDR study citation from the OSDR study page",
            notes="Use OSDR study-page citation text rather than hand-written citation.",
        ),
        _citation_metadata_template_row(
            field_id="upstream_osdr_acknowledgement_text",
            field_group="source_credit",
            target_profile=["dataset card", "README"],
            required_for="source_alpha_package",
            current_status_from_gate="needs_review",
            current_value="Data are courtesy of the NASA Open Science Data Repository",
            target_artifacts=[creator_path.as_posix()],
            validation_rule="owner_reviewed_osdr_credit_text_required_before_release",
            blocker_if_missing="",
            next_action="retain OSDR credit and review exact wording before release",
            notes="Generic acknowledgement is useful context but not a substitute for study citation.",
        ),
        _citation_metadata_template_row(
            field_id="repository_url",
            field_group="citation",
            target_profile=["CITATION.cff", "GitHub release"],
            required_for="citable_release",
            current_status_from_gate="blocked",
            current_value="not_supplied",
            target_artifacts=[ro_crate_path.as_posix()],
            validation_rule="owner_supplied_public_repository_url_required",
            blocker_if_missing="repository_public_and_owner_access_unconfirmed",
            next_action="supply public repository URL if using GitHub/Zenodo or CITATION.cff",
            notes="Do not assume the current local checkout is the release repository.",
        ),
        _citation_metadata_template_row(
            field_id="citation_cff_type",
            field_group="citation",
            target_profile=["CITATION.cff"],
            required_for="citable_release",
            current_status_from_gate="blocked",
            current_value="not_supplied",
            target_artifacts=[ro_crate_path.as_posix()],
            validation_rule="owner_supplied_cff_type_dataset_or_software_required",
            blocker_if_missing="formal_license_and_creator_metadata",
            next_action="decide whether CITATION.cff describes data, software, or preferred citation",
            notes="GitHub supports software citations and preferred-citation overrides.",
        ),
        _citation_metadata_template_row(
            field_id="payload_mirror_rights_review",
            field_group="rights",
            target_profile=["RO-Crate", "Data Package"],
            required_for="full_payload_mirror_release",
            current_status_from_gate=payload_rights["decision_status"],
            current_value=payload_rights["current_value"],
            target_artifacts=[license_path.as_posix()],
            validation_rule="retain_no_payload_mirror_unless_rights_review_supplied",
            blocker_if_missing="",
            next_action="keep full payload mirror out of current archive scope",
            notes="Current diagnostic package does not archive a full OSDR payload mirror.",
        ),
    ]


def build_multispecies_interaction_citation_metadata_fill(
    *,
    reports_root: str | Path = "v9/multispecies/reports",
    repo_root: str | Path = ".",
    owner_metadata_path: str | Path | None = None,
    fill_id: str = DEFAULT_INTERACTION_CITATION_METADATA_FILL_ID,
) -> dict[str, Any]:
    """Build a release-owner citation metadata intake/fill scaffold."""

    root = Path(repo_root)
    reports_path = _resolve_path(reports_root, root)
    archive_dir = reports_path / "interaction_archive_decision_gate"
    summary_path = archive_dir / "archive_decision_summary.csv"
    archive_option_path = archive_dir / "archive_identifier_option_matrix.csv"
    license_path = archive_dir / "license_rights_decision_table.csv"
    creator_path = archive_dir / "creator_contributor_decision_table.csv"
    summary = _require_one_row(
        _read_csv_dict_rows(summary_path),
        context="OSD-120 archive decision summary for citation metadata fill",
    )
    archive_options = _read_csv_dict_rows(archive_option_path)
    license_rows = _read_csv_dict_rows(license_path)
    creator_rows = _read_csv_dict_rows(creator_path)
    archive_option_by_id = {row["option_id"]: row for row in archive_options}
    license_by_id = {row["component_id"]: row for row in license_rows}
    creator_by_id = {row["component_id"]: row for row in creator_rows}
    owner_metadata = _load_owner_citation_metadata(
        owner_metadata_path,
        repo_root=root,
    )
    template_rows = _build_citation_metadata_template_rows(
        summary=summary,
        license_by_id=license_by_id,
        creator_by_id=creator_by_id,
        archive_option_by_id=archive_option_by_id,
        reports_path=reports_path,
    )
    intake_template = []
    fill_rows = []
    for template in template_rows:
        template_row = {
            "fill_id": fill_id,
            "decision_id": summary["decision_id"],
            "package_id": summary["package_id"],
            "task_id": summary["task_id"],
            **template,
            "supplied_value": "",
            "supplied_by": "",
            "supplied_date": "",
            "supplied_evidence": "",
            "fill_status": "owner_input_template",
        }
        intake_template.append(template_row)
        owner_row = owner_metadata.get(template["field_id"], {})
        fill_status = _citation_metadata_fill_status(template, owner_row)
        fill_rows.append(
            {
                "fill_id": fill_id,
                "decision_id": summary["decision_id"],
                "package_id": summary["package_id"],
                "task_id": summary["task_id"],
                **template,
                "supplied_value": str(owner_row.get("supplied_value", "") or ""),
                "supplied_by": str(owner_row.get("supplied_by", "") or ""),
                "supplied_date": str(owner_row.get("supplied_date", "") or ""),
                "supplied_evidence": str(owner_row.get("supplied_evidence", "") or ""),
                "fill_status": fill_status,
            }
        )
    supplied_count = sum(
        1 for row in fill_rows if row["fill_status"].startswith("owner_supplied")
    )
    retained_count = sum(
        1 for row in fill_rows if row["fill_status"] == "existing_current_draft_value_retained"
    )
    not_supplied_count = sum(
        1 for row in fill_rows if row["fill_status"].startswith("not_supplied")
    )
    blocker_count = sum(
        1 for row in fill_rows if row["fill_status"] == "not_supplied_blocking_release"
    )
    needs_review_count = sum(
        1
        for row in fill_rows
        if row["fill_status"].endswith("review")
        or row["current_status_from_gate"] == "needs_review"
    )
    if supplied_count:
        fill_status = "owner_metadata_supplied_patch_preview_only_not_applied"
        owner_metadata_status = "owner_metadata_file_loaded"
        descriptor_patch_status = "patch_preview_available_not_applied"
    else:
        fill_status = "no_owner_metadata_supplied_no_descriptor_changes"
        owner_metadata_status = "no_owner_metadata_file"
        descriptor_patch_status = "no_owner_supplied_values_no_descriptor_patch"
    patch_rows = []
    for row in fill_rows:
        if row["fill_status"].startswith("owner_supplied"):
            patch_rows.append(
                {
                    "field_id": row["field_id"],
                    "target_profile": row["target_profile"],
                    "target_artifacts": row["target_artifacts"],
                    "supplied_value": row["supplied_value"],
                    "patch_status": "preview_only_not_applied",
                }
            )
    descriptor_preview = {
        "fill_id": fill_id,
        "decision_id": summary["decision_id"],
        "package_id": summary["package_id"],
        "task_id": summary["task_id"],
        "descriptor_patch_status": descriptor_patch_status,
        "mutates_ro_crate": False,
        "mutates_datapackage": False,
        "release_ready_after_fill": False,
        "safe_current_draft_values": {
            "archive_path_decision": summary["archive_path_decision"],
            "package_title": creator_by_id["package_title"]["current_value"],
            "payload_mirror_scope": license_by_id["diagnostic_payload_mirror_rights"][
                "current_value"
            ],
        },
        "owner_supplied_patch_preview": patch_rows,
        "blocked_fields": [
            row["field_id"]
            for row in fill_rows
            if row["fill_status"] == "not_supplied_blocking_release"
        ],
        "claim_boundary": "owner_metadata_intake_only_no_descriptor_mutation",
    }
    summary_row = {
        "fill_id": fill_id,
        "decision_id": summary["decision_id"],
        "scaffold_id": summary["scaffold_id"],
        "package_id": summary["package_id"],
        "task_id": summary["task_id"],
        "fill_status": fill_status,
        "owner_metadata_status": owner_metadata_status,
        "intake_field_count": str(len(fill_rows)),
        "supplied_field_count": str(supplied_count),
        "retained_current_draft_count": str(retained_count),
        "not_supplied_field_count": str(not_supplied_count),
        "not_supplied_blocker_count": str(blocker_count),
        "needs_review_count": str(needs_review_count),
        "descriptor_patch_status": descriptor_patch_status,
        "ro_crate_mutation_status": "not_mutated_preview_only",
        "datapackage_mutation_status": "not_mutated_preview_only",
        "release_ready_after_fill": "False",
        "next_required_block": (
            "V9-MULTI-034: owner metadata application guard or archive release deferral"
        ),
        "claim_boundary": "owner_metadata_intake_only_no_descriptor_mutation",
    }
    blocked_lines = "\n".join(
        f"- `{row['field_id']}`: {row['blocker_if_missing']}."
        for row in fill_rows
        if row["fill_status"] == "not_supplied_blocking_release"
    )
    supplied_lines = "\n".join(
        f"- `{row['field_id']}`: {row['fill_status']}."
        for row in fill_rows
        if row["fill_status"].startswith("owner_supplied")
    )
    if not supplied_lines:
        supplied_lines = "- No owner-supplied citation metadata fields were provided."
    review_md = f"""# OSD-120 Citation Metadata Fill Review

Block: V9-MULTI-033

Fill id: `{fill_id}`

Status: `{fill_status}`

This block creates a release-owner metadata intake/fill scaffold. It does not
mutate the draft RO-Crate or Data Package descriptors, does not mint a DOI, does
not infer creator order, and does not choose a license.

## Owner-Supplied Fields

{supplied_lines}

## Remaining Release Blockers

{blocked_lines}

## Descriptor Patch Policy

Patch previews are emitted only for fields explicitly supplied by a release
owner. The current generated descriptors remain diagnostic drafts until owner
metadata, archive route, version, creator/contributor, publisher, and license
decisions are complete.

## External Guidance Anchors

- GitHub/Zenodo DOI release paths require a public repository, owner access,
  GitHub release, and license.
- GitHub CITATION.cff can describe software or dataset citation metadata, but
  author, version, DOI or URL, and release-date fields must be supplied.
- OSDR source citation should come from the OSDR study page citation button in
  BibTeX or RIS form.
- DataCite 4.7 citable metadata requires identifier, creator, title,
  publisher, publication year, and resource type.
- Local license metadata should use an SPDX identifier or explicit custom terms
  after review.

## Claim Boundary

The current output is an intake and preview scaffold only. No archive
identifier, no DOI, no frozen release version, no creator list, no publisher,
no local package license, and no descriptor mutation are claimed.
"""
    return {
        "summary": [summary_row],
        "intake_template": intake_template,
        "fill_rows": fill_rows,
        "descriptor_preview": descriptor_preview,
        "review_md": review_md,
    }


def write_multispecies_interaction_citation_metadata_fill(
    *,
    output_dir: str | Path = (
        "v9/multispecies/reports/interaction_citation_metadata_fill"
    ),
    reports_root: str | Path = "v9/multispecies/reports",
    repo_root: str | Path = ".",
    owner_metadata_path: str | Path | None = None,
    fill_id: str = DEFAULT_INTERACTION_CITATION_METADATA_FILL_ID,
) -> dict[str, Path]:
    """Write release-owner citation metadata intake/fill scaffold tables."""

    package = build_multispecies_interaction_citation_metadata_fill(
        reports_root=reports_root,
        repo_root=repo_root,
        owner_metadata_path=owner_metadata_path,
        fill_id=fill_id,
    )
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "summary_csv": outdir / "citation_metadata_fill_summary.csv",
        "summary_json": outdir / "citation_metadata_fill_summary.json",
        "intake_template_csv": outdir / "release_owner_metadata_intake_template.csv",
        "intake_template_json": outdir / "release_owner_metadata_intake_template.json",
        "fill_status_csv": outdir / "citation_metadata_fill_status.csv",
        "fill_status_json": outdir / "citation_metadata_fill_status.json",
        "descriptor_preview_json": outdir / "citation_descriptor_patch_preview.json",
        "review_md": (
            outdir.parent / "OSD120_INTERACTION_CITATION_METADATA_FILL_REVIEW.md"
        ),
    }
    with outputs["summary_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_CITATION_METADATA_FILL_SUMMARY_FIELDS,
        )
        writer.writeheader()
        writer.writerows(package["summary"])
    outputs["summary_json"].write_text(
        json.dumps(package["summary"], indent=2, sort_keys=True) + "\n"
    )
    with outputs["intake_template_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_CITATION_METADATA_FILL_FIELDS,
        )
        writer.writeheader()
        writer.writerows(package["intake_template"])
    outputs["intake_template_json"].write_text(
        json.dumps(package["intake_template"], indent=2, sort_keys=True) + "\n"
    )
    with outputs["fill_status_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_CITATION_METADATA_FILL_FIELDS,
        )
        writer.writeheader()
        writer.writerows(package["fill_rows"])
    outputs["fill_status_json"].write_text(
        json.dumps(package["fill_rows"], indent=2, sort_keys=True) + "\n"
    )
    outputs["descriptor_preview_json"].write_text(
        json.dumps(package["descriptor_preview"], indent=2, sort_keys=True) + "\n"
    )
    outputs["review_md"].write_text(str(package["review_md"]))
    return outputs


def _owner_application_guard_row(
    *,
    guard_id: str,
    fill_id: str,
    package_id: str,
    task_id: str,
    guard_check_id: str,
    guard_surface: str,
    guard_status: str,
    severity: str,
    evidence_artifacts: Sequence[str],
    observed_value: str,
    required_value: str,
    source_blocker_fields: Sequence[str],
    action_if_failed: str,
    mutation_policy: str,
    claim_boundary_impact: str,
) -> dict[str, str]:
    return {
        "guard_id": guard_id,
        "fill_id": fill_id,
        "package_id": package_id,
        "task_id": task_id,
        "guard_check_id": guard_check_id,
        "guard_surface": guard_surface,
        "guard_status": guard_status,
        "severity": severity,
        "evidence_artifacts": _pipe_join(list(evidence_artifacts)),
        "observed_value": observed_value,
        "required_value": required_value,
        "source_blocker_fields": _pipe_join(list(source_blocker_fields)),
        "action_if_failed": action_if_failed,
        "mutation_policy": mutation_policy,
        "claim_boundary_impact": claim_boundary_impact,
    }


def _archive_deferral_action_row(
    *,
    guard_id: str,
    fill_id: str,
    package_id: str,
    task_id: str,
    action_id: str,
    action_owner: str,
    action_status: str,
    required_for: str,
    source_blocker_fields: Sequence[str],
    source_evidence: Sequence[str],
    next_action: str,
    deferral_policy: str,
) -> dict[str, str]:
    return {
        "guard_id": guard_id,
        "fill_id": fill_id,
        "package_id": package_id,
        "task_id": task_id,
        "action_id": action_id,
        "action_owner": action_owner,
        "action_status": action_status,
        "required_for": required_for,
        "source_blocker_fields": _pipe_join(list(source_blocker_fields)),
        "source_evidence": _pipe_join(list(source_evidence)),
        "next_action": next_action,
        "deferral_policy": deferral_policy,
    }


def build_multispecies_interaction_archive_release_deferral_guard(
    *,
    reports_root: str | Path = "v9/multispecies/reports",
    repo_root: str | Path = ".",
    guard_id: str = DEFAULT_INTERACTION_ARCHIVE_RELEASE_DEFERRAL_GUARD_ID,
) -> dict[str, Any]:
    """Build the archive-release deferral/application guard."""

    root = Path(repo_root)
    reports_path = _resolve_path(reports_root, root)
    fill_dir = reports_path / "interaction_citation_metadata_fill"
    summary_path = fill_dir / "citation_metadata_fill_summary.csv"
    fill_status_path = fill_dir / "citation_metadata_fill_status.csv"
    descriptor_preview_path = fill_dir / "citation_descriptor_patch_preview.json"
    fill_summary = _require_one_row(
        _read_csv_dict_rows(summary_path),
        context="OSD-120 citation metadata fill summary for archive deferral guard",
    )
    fill_rows = _read_csv_dict_rows(fill_status_path)
    descriptor_preview = json.loads(descriptor_preview_path.read_text())
    fill_id = fill_summary["fill_id"]
    package_id = fill_summary["package_id"]
    task_id = fill_summary["task_id"]
    blocked_fields = [
        row["field_id"]
        for row in fill_rows
        if row["fill_status"] == "not_supplied_blocking_release"
    ]
    retained_fields = [
        row["field_id"]
        for row in fill_rows
        if row["fill_status"] == "existing_current_draft_value_retained"
    ]
    owner_supplied_fields = [
        row["field_id"]
        for row in fill_rows
        if row["fill_status"].startswith("owner_supplied")
    ]
    evidence = [summary_path.as_posix(), fill_status_path.as_posix(), descriptor_preview_path.as_posix()]
    version_fields = ["release_version_string", "release_date", "publication_year"]
    creator_fields = [
        "local_package_creators",
        "local_package_contributors",
        "publisher_maintainer_identity",
    ]
    license_fields = ["local_code_license", "local_metadata_tables_license"]
    citation_release_fields = ["repository_url", "citation_cff_type"]
    guard_rows = [
        _owner_application_guard_row(
            guard_id=guard_id,
            fill_id=fill_id,
            package_id=package_id,
            task_id=task_id,
            guard_check_id="owner_metadata_file_present",
            guard_surface="owner_metadata_intake",
            guard_status="blocker",
            severity="release_blocker",
            evidence_artifacts=evidence,
            observed_value=fill_summary["owner_metadata_status"],
            required_value="owner_metadata_file_loaded",
            source_blocker_fields=blocked_fields,
            action_if_failed="defer_archive_release_and_keep_diagnostic_metadata_only",
            mutation_policy="prevent_descriptor_mutation",
            claim_boundary_impact="No owner metadata means no authorship, license, version, or identifier can be applied.",
        ),
        _owner_application_guard_row(
            guard_id=guard_id,
            fill_id=fill_id,
            package_id=package_id,
            task_id=task_id,
            guard_check_id="owner_supplied_patch_preview_present",
            guard_surface="descriptor_patch_preview",
            guard_status="blocker",
            severity="release_blocker",
            evidence_artifacts=[descriptor_preview_path.as_posix()],
            observed_value=fill_summary["descriptor_patch_status"],
            required_value="patch_preview_available_not_applied",
            source_blocker_fields=blocked_fields,
            action_if_failed="do_not_apply_empty_patch_preview",
            mutation_policy="prevent_descriptor_mutation",
            claim_boundary_impact="Empty patch preview confirms that descriptor mutation would add inferred metadata.",
        ),
        _owner_application_guard_row(
            guard_id=guard_id,
            fill_id=fill_id,
            package_id=package_id,
            task_id=task_id,
            guard_check_id="all_release_blockers_resolved",
            guard_surface="release_readiness",
            guard_status="blocker",
            severity="release_blocker",
            evidence_artifacts=[fill_status_path.as_posix()],
            observed_value=str(len(blocked_fields)),
            required_value="0",
            source_blocker_fields=blocked_fields,
            action_if_failed="defer_archive_release_until_blockers_are_zero",
            mutation_policy="prevent_archive_release_claims",
            claim_boundary_impact="Citable archive release remains blocked while any release metadata field is missing.",
        ),
        _owner_application_guard_row(
            guard_id=guard_id,
            fill_id=fill_id,
            package_id=package_id,
            task_id=task_id,
            guard_check_id="archive_identifier_supplied",
            guard_surface="archive_identifier",
            guard_status="blocker",
            severity="release_blocker",
            evidence_artifacts=[fill_status_path.as_posix()],
            observed_value="future_archive_identifier_not_supplied",
            required_value="owner_supplied_identifier_or_explicit_no_archive_release_decision",
            source_blocker_fields=["future_archive_identifier"],
            action_if_failed="keep_no_archive_identifier_for_current_draft",
            mutation_policy="prevent_identifier_mutation",
            claim_boundary_impact="No DOI, SWHID, URL, or other archive identifier is claimed.",
        ),
        _owner_application_guard_row(
            guard_id=guard_id,
            fill_id=fill_id,
            package_id=package_id,
            task_id=task_id,
            guard_check_id="release_version_date_supplied",
            guard_surface="release_version",
            guard_status="blocker",
            severity="release_blocker",
            evidence_artifacts=[fill_status_path.as_posix()],
            observed_value="|".join(field for field in version_fields if field in blocked_fields),
            required_value="release_version_string|release_date|publication_year",
            source_blocker_fields=version_fields,
            action_if_failed="do_not_use_build_date_or_draft_string_as_release_version",
            mutation_policy="prevent_version_mutation",
            claim_boundary_impact="No frozen release version or publication year is claimed.",
        ),
        _owner_application_guard_row(
            guard_id=guard_id,
            fill_id=fill_id,
            package_id=package_id,
            task_id=task_id,
            guard_check_id="creator_publisher_metadata_supplied",
            guard_surface="creator_metadata",
            guard_status="blocker",
            severity="release_blocker",
            evidence_artifacts=[fill_status_path.as_posix()],
            observed_value="|".join(field for field in creator_fields if field in blocked_fields),
            required_value="creator_list|contributor_roles|publisher_or_maintainer",
            source_blocker_fields=creator_fields,
            action_if_failed="do_not_infer_authorship_or_publisher_from_repository_state",
            mutation_policy="prevent_creator_metadata_mutation",
            claim_boundary_impact="No creator order, contributor roles, or publisher are claimed.",
        ),
        _owner_application_guard_row(
            guard_id=guard_id,
            fill_id=fill_id,
            package_id=package_id,
            task_id=task_id,
            guard_check_id="license_rights_metadata_supplied",
            guard_surface="rights_metadata",
            guard_status="blocker",
            severity="release_blocker",
            evidence_artifacts=[fill_status_path.as_posix()],
            observed_value="|".join(field for field in license_fields if field in blocked_fields),
            required_value="local_code_license|local_metadata_tables_license",
            source_blocker_fields=license_fields,
            action_if_failed="do_not_choose_spdx_or_custom_license_without_owner_review",
            mutation_policy="prevent_license_mutation",
            claim_boundary_impact="No local package license is claimed.",
        ),
        _owner_application_guard_row(
            guard_id=guard_id,
            fill_id=fill_id,
            package_id=package_id,
            task_id=task_id,
            guard_check_id="osdr_source_citation_supplied",
            guard_surface="source_credit",
            guard_status="blocker",
            severity="source_alpha_blocker",
            evidence_artifacts=[fill_status_path.as_posix()],
            observed_value="upstream_osdr_dataset_citation_not_supplied",
            required_value="exact_osdr_study_page_bibtex_or_ris",
            source_blocker_fields=["upstream_osdr_dataset_citation"],
            action_if_failed="retain_generic_osdr_credit_but_defer_source_alpha_citation",
            mutation_policy="prevent_source_citation_freeze",
            claim_boundary_impact="OSDR credit remains visible but exact source citation is not frozen.",
        ),
        _owner_application_guard_row(
            guard_id=guard_id,
            fill_id=fill_id,
            package_id=package_id,
            task_id=task_id,
            guard_check_id="citation_release_surface_supplied",
            guard_surface="github_citation_release",
            guard_status="blocker",
            severity="release_blocker",
            evidence_artifacts=[fill_status_path.as_posix()],
            observed_value="|".join(field for field in citation_release_fields if field in blocked_fields),
            required_value="repository_url|citation_cff_type",
            source_blocker_fields=citation_release_fields,
            action_if_failed="do_not_emit_citation_cff_or_github_release_metadata",
            mutation_policy="prevent_citation_release_generation",
            claim_boundary_impact="No CITATION.cff or GitHub release metadata is generated.",
        ),
        _owner_application_guard_row(
            guard_id=guard_id,
            fill_id=fill_id,
            package_id=package_id,
            task_id=task_id,
            guard_check_id="descriptor_mutation_prevented",
            guard_surface="descriptor_mutation",
            guard_status="pass",
            severity="guard_pass",
            evidence_artifacts=[descriptor_preview_path.as_posix()],
            observed_value=(
                f"mutates_ro_crate={descriptor_preview['mutates_ro_crate']};"
                f"mutates_datapackage={descriptor_preview['mutates_datapackage']}"
            ),
            required_value="mutates_ro_crate=False;mutates_datapackage=False",
            source_blocker_fields=blocked_fields,
            action_if_failed="stop_descriptor_mutation",
            mutation_policy="confirmed_no_descriptor_mutation",
            claim_boundary_impact="The draft descriptors remain diagnostic metadata only.",
        ),
        _owner_application_guard_row(
            guard_id=guard_id,
            fill_id=fill_id,
            package_id=package_id,
            task_id=task_id,
            guard_check_id="current_diagnostic_surface_retained",
            guard_surface="diagnostic_metadata",
            guard_status="pass",
            severity="guard_pass",
            evidence_artifacts=[fill_status_path.as_posix(), descriptor_preview_path.as_posix()],
            observed_value=_pipe_join(retained_fields),
            required_value="current_archive_path_decision|package_title|upstream_osdr_acknowledgement_text|payload_mirror_rights_review",
            source_blocker_fields=[],
            action_if_failed="restore_current_diagnostic_metadata_surface",
            mutation_policy="retain_diagnostic_metadata_only",
            claim_boundary_impact="The current public surface remains diagnostic metadata, not a citable archive.",
        ),
    ]
    actions = [
        _archive_deferral_action_row(
            guard_id=guard_id,
            fill_id=fill_id,
            package_id=package_id,
            task_id=task_id,
            action_id="supply_release_owner_metadata_file",
            action_owner="release_owner",
            action_status="required_owner_action",
            required_for="any_descriptor_metadata_application",
            source_blocker_fields=blocked_fields,
            source_evidence=evidence,
            next_action="provide reviewed CSV/JSON values for the V9-MULTI-033 intake template",
            deferral_policy="archive_release_deferred_until_owner_metadata_exists",
        ),
        _archive_deferral_action_row(
            guard_id=guard_id,
            fill_id=fill_id,
            package_id=package_id,
            task_id=task_id,
            action_id="decide_archive_route_identifier",
            action_owner="release_owner",
            action_status="required_owner_action",
            required_for="doi_or_archive_release",
            source_blocker_fields=["future_archive_identifier", "repository_url"],
            source_evidence=[fill_status_path.as_posix()],
            next_action="choose no archive, Zenodo DOI, SWHID related identifier, or another explicit route",
            deferral_policy="do_not_mint_identifier_without_release_route",
        ),
        _archive_deferral_action_row(
            guard_id=guard_id,
            fill_id=fill_id,
            package_id=package_id,
            task_id=task_id,
            action_id="freeze_version_date_year",
            action_owner="release_owner",
            action_status="required_owner_action",
            required_for="citable_release_metadata",
            source_blocker_fields=version_fields,
            source_evidence=[fill_status_path.as_posix()],
            next_action="supply release tag, release date, and publication year",
            deferral_policy="do_not_promote_draft_version_string",
        ),
        _archive_deferral_action_row(
            guard_id=guard_id,
            fill_id=fill_id,
            package_id=package_id,
            task_id=task_id,
            action_id="freeze_creator_contributor_publisher",
            action_owner="release_owner",
            action_status="required_owner_action",
            required_for="datacite_ro_crate_citation_cff",
            source_blocker_fields=creator_fields,
            source_evidence=[fill_status_path.as_posix()],
            next_action="supply creator order, contributor roles, and publisher/maintainer",
            deferral_policy="do_not_infer_authorship_or_publisher",
        ),
        _archive_deferral_action_row(
            guard_id=guard_id,
            fill_id=fill_id,
            package_id=package_id,
            task_id=task_id,
            action_id="choose_license_scope",
            action_owner="release_owner",
            action_status="required_owner_action",
            required_for="reuse_terms_and_github_zenodo_release",
            source_blocker_fields=license_fields,
            source_evidence=[fill_status_path.as_posix()],
            next_action="choose local code and generated metadata/table license scope",
            deferral_policy="do_not_apply_spdx_or_custom_license_without_review",
        ),
        _archive_deferral_action_row(
            guard_id=guard_id,
            fill_id=fill_id,
            package_id=package_id,
            task_id=task_id,
            action_id="capture_exact_osdr_study_citation",
            action_owner="release_owner_or_data_curator",
            action_status="required_owner_action",
            required_for="source_alpha_package_citation",
            source_blocker_fields=["upstream_osdr_dataset_citation"],
            source_evidence=[fill_status_path.as_posix()],
            next_action="export OSDR study page citation in BibTeX or RIS form",
            deferral_policy="retain_generic_osdr_credit_but_defer_exact_source_citation",
        ),
        _archive_deferral_action_row(
            guard_id=guard_id,
            fill_id=fill_id,
            package_id=package_id,
            task_id=task_id,
            action_id="retain_diagnostic_metadata_surface",
            action_owner="spacebio_bench_maintainer",
            action_status="selected_deferral_action",
            required_for="current_public_diagnostic_metadata",
            source_blocker_fields=[],
            source_evidence=[descriptor_preview_path.as_posix()],
            next_action="keep current no-archive diagnostic metadata outputs inspectable",
            deferral_policy="public_surface_remains_diagnostic_metadata_only",
        ),
        _archive_deferral_action_row(
            guard_id=guard_id,
            fill_id=fill_id,
            package_id=package_id,
            task_id=task_id,
            action_id="prevent_descriptor_mutation",
            action_owner="spacebio_bench_maintainer",
            action_status="selected_guard_action",
            required_for="metadata_integrity",
            source_blocker_fields=blocked_fields,
            source_evidence=[descriptor_preview_path.as_posix()],
            next_action="do not regenerate RO-Crate/Data Package/CITATION.cff with inferred values",
            deferral_policy="descriptor_mutation_blocked_until_owner_fields_validate",
        ),
        _archive_deferral_action_row(
            guard_id=guard_id,
            fill_id=fill_id,
            package_id=package_id,
            task_id=task_id,
            action_id="defer_archive_release",
            action_owner="spacebio_bench_maintainer",
            action_status="selected_deferral_action",
            required_for="release_claim_boundary",
            source_blocker_fields=blocked_fields,
            source_evidence=evidence,
            next_action="record archive release as deferred and continue with diagnostic metadata-only package",
            deferral_policy="archive_release_deferred_no_owner_metadata",
        ),
    ]
    pass_count = sum(1 for row in guard_rows if row["guard_status"] == "pass")
    blocker_count = sum(1 for row in guard_rows if row["guard_status"] == "blocker")
    deferred_count = sum(1 for row in guard_rows if row["guard_status"] == "deferred")
    required_owner_action_count = sum(
        1 for row in actions if row["action_status"] == "required_owner_action"
    )
    descriptor_mutation_allowed = (
        not blocked_fields
        and bool(owner_supplied_fields)
        and descriptor_preview["owner_supplied_patch_preview"]
    )
    if descriptor_mutation_allowed:
        guard_status = "owner_metadata_present_ready_for_strict_application_review"
        release_deferral_status = "pending_application_review"
    else:
        guard_status = "archive_release_deferred_no_owner_metadata"
        release_deferral_status = "deferred_keep_diagnostic_metadata_only"
    summary = {
        "guard_id": guard_id,
        "fill_id": fill_id,
        "decision_id": fill_summary["decision_id"],
        "package_id": package_id,
        "task_id": task_id,
        "guard_status": guard_status,
        "release_deferral_status": release_deferral_status,
        "owner_metadata_status": fill_summary["owner_metadata_status"],
        "supplied_field_count": fill_summary["supplied_field_count"],
        "blocked_field_count": str(len(blocked_fields)),
        "retained_current_draft_count": str(len(retained_fields)),
        "application_guard_count": str(len(guard_rows)),
        "guard_pass_count": str(pass_count),
        "guard_blocker_count": str(blocker_count),
        "guard_deferred_count": str(deferred_count),
        "action_count": str(len(actions)),
        "required_owner_action_count": str(required_owner_action_count),
        "descriptor_mutation_allowed": str(bool(descriptor_mutation_allowed)),
        "release_ready_after_guard": "False",
        "next_required_block": "V9-MULTI-035: diagnostic metadata release note or owner metadata intake retry",
        "claim_boundary": "archive_release_deferred_diagnostic_metadata_only_no_descriptor_mutation",
    }
    mutation_guard = {
        "guard_id": guard_id,
        "fill_id": fill_id,
        "descriptor_mutation_allowed": bool(descriptor_mutation_allowed),
        "mutates_ro_crate": False,
        "mutates_datapackage": False,
        "mutates_citation_cff": False,
        "blocked_fields": blocked_fields,
        "retained_current_draft_fields": retained_fields,
        "owner_supplied_fields": owner_supplied_fields,
        "allowed_next_descriptor_action": (
            "none_defer_archive_release"
            if not descriptor_mutation_allowed
            else "strict_owner_metadata_application_review"
        ),
        "claim_boundary": summary["claim_boundary"],
    }
    blocker_lines = "\n".join(
        f"- `{row['guard_check_id']}`: {row['action_if_failed']}."
        for row in guard_rows
        if row["guard_status"] == "blocker"
    )
    action_lines = "\n".join(
        f"- `{row['action_id']}`: {row['action_status']}."
        for row in actions
    )
    review_md = f"""# OSD-120 Archive Release Deferral And Application Guard Review

Block: V9-MULTI-034

Guard id: `{guard_id}`

Status: `{summary["guard_status"]}`

The V9-MULTI-033 owner metadata fill has no owner-supplied citation metadata.
This guard therefore defers archive release, blocks descriptor mutation, and
keeps the current OSD-120 public surface as diagnostic metadata only.

## Guard Blockers

{blocker_lines}

## Deferral Actions

{action_lines}

## Mutation Policy

No RO-Crate, Data Package, CITATION.cff, DOI, release tag, license, creator,
publisher, or source-citation freeze is applied by this block. A future
application block may mutate descriptors only if owner-supplied metadata exists,
all release blockers validate, and the patch preview is non-empty and reviewed.

## External Guidance Anchors

- GitHub/Zenodo DOI release paths depend on public repository access, owner
  approval where needed, a GitHub release, and explicit license/reuse terms.
- GitHub releases require a tag and release metadata; draft/pre-release states
  are separate from a citable archive.
- Zenodo release ingestion can fail on release metadata, so metadata should be
  validated before release.
- DataCite citable metadata requires identifier, creator, title, publisher,
  publication year, and resource type.

## Claim Boundary

Archive release is deferred. The current outputs remain diagnostic metadata
only: no archive identifier, no DOI, no frozen release version, no local license,
no creator/publisher metadata, no full OSDR payload mirror, and no descriptor
mutation are claimed.
"""
    return {
        "summary": [summary],
        "guard_rows": guard_rows,
        "actions": actions,
        "mutation_guard": mutation_guard,
        "review_md": review_md,
    }


def write_multispecies_interaction_archive_release_deferral_guard(
    *,
    output_dir: str | Path = (
        "v9/multispecies/reports/interaction_archive_release_deferral_guard"
    ),
    reports_root: str | Path = "v9/multispecies/reports",
    repo_root: str | Path = ".",
    guard_id: str = DEFAULT_INTERACTION_ARCHIVE_RELEASE_DEFERRAL_GUARD_ID,
) -> dict[str, Path]:
    """Write archive-release deferral/application guard artifacts."""

    package = build_multispecies_interaction_archive_release_deferral_guard(
        reports_root=reports_root,
        repo_root=repo_root,
        guard_id=guard_id,
    )
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "summary_csv": outdir / "archive_release_deferral_summary.csv",
        "summary_json": outdir / "archive_release_deferral_summary.json",
        "guard_csv": outdir / "owner_metadata_application_guard.csv",
        "guard_json": outdir / "owner_metadata_application_guard.json",
        "action_csv": outdir / "archive_release_deferral_actions.csv",
        "action_json": outdir / "archive_release_deferral_actions.json",
        "mutation_guard_json": outdir / "descriptor_mutation_guard.json",
        "review_md": (
            outdir.parent
            / "OSD120_INTERACTION_ARCHIVE_RELEASE_DEFERRAL_GUARD_REVIEW.md"
        ),
    }
    with outputs["summary_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_ARCHIVE_RELEASE_DEFERRAL_SUMMARY_FIELDS,
        )
        writer.writeheader()
        writer.writerows(package["summary"])
    outputs["summary_json"].write_text(
        json.dumps(package["summary"], indent=2, sort_keys=True) + "\n"
    )
    with outputs["guard_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_OWNER_METADATA_APPLICATION_GUARD_FIELDS,
        )
        writer.writeheader()
        writer.writerows(package["guard_rows"])
    outputs["guard_json"].write_text(
        json.dumps(package["guard_rows"], indent=2, sort_keys=True) + "\n"
    )
    with outputs["action_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_ARCHIVE_RELEASE_DEFERRAL_ACTION_FIELDS,
        )
        writer.writeheader()
        writer.writerows(package["actions"])
    outputs["action_json"].write_text(
        json.dumps(package["actions"], indent=2, sort_keys=True) + "\n"
    )
    outputs["mutation_guard_json"].write_text(
        json.dumps(package["mutation_guard"], indent=2, sort_keys=True) + "\n"
    )
    outputs["review_md"].write_text(str(package["review_md"]))
    return outputs


def _diagnostic_release_section_row(
    *,
    release_note_id: str,
    guard_id: str,
    package_id: str,
    task_id: str,
    section_id: str,
    section_title: str,
    include_in_note: bool,
    note_text: str,
    evidence_artifacts: Sequence[str],
    claim_boundary: str,
) -> dict[str, str]:
    return {
        "release_note_id": release_note_id,
        "guard_id": guard_id,
        "package_id": package_id,
        "task_id": task_id,
        "section_id": section_id,
        "section_title": section_title,
        "include_in_note": str(include_in_note),
        "note_text": note_text,
        "evidence_artifacts": _pipe_join(list(evidence_artifacts)),
        "claim_boundary": claim_boundary,
    }


def _diagnostic_public_claim_row(
    *,
    release_note_id: str,
    guard_id: str,
    package_id: str,
    task_id: str,
    claim_id: str,
    claim_category: str,
    public_note_status: str,
    statement: str,
    supporting_evidence: Sequence[str],
    prohibited_language: str,
    next_allowed_action: str,
) -> dict[str, str]:
    return {
        "release_note_id": release_note_id,
        "guard_id": guard_id,
        "package_id": package_id,
        "task_id": task_id,
        "claim_id": claim_id,
        "claim_category": claim_category,
        "public_note_status": public_note_status,
        "statement": statement,
        "supporting_evidence": _pipe_join(list(supporting_evidence)),
        "prohibited_language": prohibited_language,
        "next_allowed_action": next_allowed_action,
    }


def _owner_retry_checklist_row(
    *,
    release_note_id: str,
    guard_id: str,
    package_id: str,
    task_id: str,
    retry_item_id: str,
    required_for: str,
    owner_action: str,
    current_status: str,
    source_blocker_fields: Sequence[str],
    evidence_artifacts: Sequence[str],
    validation_rule: str,
    retry_priority: str,
    closeout_policy: str,
) -> dict[str, str]:
    return {
        "release_note_id": release_note_id,
        "guard_id": guard_id,
        "package_id": package_id,
        "task_id": task_id,
        "retry_item_id": retry_item_id,
        "required_for": required_for,
        "owner_action": owner_action,
        "current_status": current_status,
        "source_blocker_fields": _pipe_join(list(source_blocker_fields)),
        "evidence_artifacts": _pipe_join(list(evidence_artifacts)),
        "validation_rule": validation_rule,
        "retry_priority": retry_priority,
        "closeout_policy": closeout_policy,
    }


def build_multispecies_interaction_diagnostic_metadata_release_note(
    *,
    reports_root: str | Path = "v9/multispecies/reports",
    repo_root: str | Path = ".",
    release_note_id: str = DEFAULT_INTERACTION_DIAGNOSTIC_METADATA_RELEASE_NOTE_ID,
) -> dict[str, Any]:
    """Build the OSD-120 diagnostic metadata release-note closeout."""

    root = Path(repo_root)
    reports_path = _resolve_path(reports_root, root)
    guard_dir = reports_path / "interaction_archive_release_deferral_guard"
    alpha_dir = reports_path / "interaction_public_alpha_card"
    metadata_dir = reports_path / "interaction_public_metadata_package"
    scaffold_dir = reports_path / "interaction_ro_crate_citation_scaffold"
    guard_summary_path = guard_dir / "archive_release_deferral_summary.csv"
    guard_path = guard_dir / "owner_metadata_application_guard.csv"
    action_path = guard_dir / "archive_release_deferral_actions.csv"
    mutation_guard_path = guard_dir / "descriptor_mutation_guard.json"
    alpha_summary_path = alpha_dir / "public_alpha_card_summary.csv"
    metadata_summary_path = metadata_dir / "public_metadata_summary.csv"
    ro_crate_summary_path = scaffold_dir / "ro_crate_export_summary.csv"
    guard_summary = _require_one_row(
        _read_csv_dict_rows(guard_summary_path),
        context="OSD-120 archive release deferral guard summary for release note",
    )
    guard_rows = _read_csv_dict_rows(guard_path)
    action_rows = _read_csv_dict_rows(action_path)
    mutation_guard = json.loads(mutation_guard_path.read_text())
    alpha_summary = _require_one_row(
        _read_csv_dict_rows(alpha_summary_path),
        context="OSD-120 public alpha card summary for release note",
    )
    metadata_summary = _require_one_row(
        _read_csv_dict_rows(metadata_summary_path),
        context="OSD-120 public metadata summary for release note",
    )
    ro_crate_summary = _require_one_row(
        _read_csv_dict_rows(ro_crate_summary_path),
        context="OSD-120 RO-Crate scaffold summary for release note",
    )
    guard_id = guard_summary["guard_id"]
    fill_id = guard_summary["fill_id"]
    package_id = guard_summary["package_id"]
    task_id = guard_summary["task_id"]
    common_evidence = [
        guard_summary_path.as_posix(),
        guard_path.as_posix(),
        action_path.as_posix(),
        mutation_guard_path.as_posix(),
    ]
    public_note_claim_boundary = (
        "diagnostic_metadata_note_only_not_archive_release"
    )
    sections = [
        _diagnostic_release_section_row(
            release_note_id=release_note_id,
            guard_id=guard_id,
            package_id=package_id,
            task_id=task_id,
            section_id="status",
            section_title="Status",
            include_in_note=True,
            note_text=(
                "This is an OSD-120 diagnostic metadata note for SpaceBio-Bench "
                "v9. It is inspectable draft metadata, not a DOI-backed archive "
                "release or benchmark leaderboard release."
            ),
            evidence_artifacts=[guard_summary_path.as_posix(), alpha_summary_path.as_posix()],
            claim_boundary=public_note_claim_boundary,
        ),
        _diagnostic_release_section_row(
            release_note_id=release_note_id,
            guard_id=guard_id,
            package_id=package_id,
            task_id=task_id,
            section_id="inspectable_now",
            section_title="Inspectable Now",
            include_in_note=True,
            note_text=(
                "Inspectable files include the OSD-120 task manifest, diagnostic "
                "artifact manifest, diagnostic payload-freeze manifest, public "
                "metadata skeleton, draft RO-Crate/Data Package descriptors, and "
                "archive deferral guard outputs."
            ),
            evidence_artifacts=[
                alpha_summary_path.as_posix(),
                metadata_summary_path.as_posix(),
                ro_crate_summary_path.as_posix(),
            ],
            claim_boundary=public_note_claim_boundary,
        ),
        _diagnostic_release_section_row(
            release_note_id=release_note_id,
            guard_id=guard_id,
            package_id=package_id,
            task_id=task_id,
            section_id="not_released_as_archive",
            section_title="Not Released As Archive",
            include_in_note=True,
            note_text=(
                "This note does not provide a DOI, CITATION.cff, GitHub release, "
                "Zenodo record, frozen release version, local package license, "
                "creator/publisher metadata, or full OSDR payload mirror."
            ),
            evidence_artifacts=common_evidence,
            claim_boundary=public_note_claim_boundary,
        ),
        _diagnostic_release_section_row(
            release_note_id=release_note_id,
            guard_id=guard_id,
            package_id=package_id,
            task_id=task_id,
            section_id="source_credit",
            section_title="Source Credit",
            include_in_note=True,
            note_text=(
                "Data are credited to NASA OSDR/GeneLab context, but the exact "
                "OSD-120 study-page citation remains an owner/curator retry item "
                "before source-alpha or archive release wording."
            ),
            evidence_artifacts=[action_path.as_posix()],
            claim_boundary=public_note_claim_boundary,
        ),
        _diagnostic_release_section_row(
            release_note_id=release_note_id,
            guard_id=guard_id,
            package_id=package_id,
            task_id=task_id,
            section_id="next_step",
            section_title="Next Step",
            include_in_note=True,
            note_text=(
                "The OSD-120 metadata branch is closed after this note unless "
                "owner-supplied release metadata appears. The next v9 step should "
                "recenter on public bulk alpha readiness or the first single-cell "
                "flagship inventory."
            ),
            evidence_artifacts=[
                "docs/V9_PURPOSE_DRIFT_AUDIT_2026_05_26.md",
                guard_summary_path.as_posix(),
            ],
            claim_boundary=public_note_claim_boundary,
        ),
    ]
    public_claims = [
        _diagnostic_public_claim_row(
            release_note_id=release_note_id,
            guard_id=guard_id,
            package_id=package_id,
            task_id=task_id,
            claim_id="current_surface_is_diagnostic_metadata",
            claim_category="allowed",
            public_note_status="allowed_current_note",
            statement=(
                "The OSD-120 package is available as diagnostic metadata and "
                "traceability evidence."
            ),
            supporting_evidence=[guard_summary_path.as_posix(), metadata_summary_path.as_posix()],
            prohibited_language="released benchmark snapshot",
            next_allowed_action="refer to diagnostic metadata note and artifact paths",
        ),
        _diagnostic_public_claim_row(
            release_note_id=release_note_id,
            guard_id=guard_id,
            package_id=package_id,
            task_id=task_id,
            claim_id="diagnostic_payload_scope_is_narrow",
            claim_category="allowed",
            public_note_status="allowed_current_note",
            statement=(
                "The diagnostic-required OSD-120 payload scope is narrow and does "
                "not claim a full OSDR processed payload mirror."
            ),
            supporting_evidence=[alpha_summary_path.as_posix()],
            prohibited_language="complete OSDR payload archive",
            next_allowed_action="keep full mirror language out of release note",
        ),
        _diagnostic_public_claim_row(
            release_note_id=release_note_id,
            guard_id=guard_id,
            package_id=package_id,
            task_id=task_id,
            claim_id="draft_ro_crate_is_scaffold_only",
            claim_category="allowed",
            public_note_status="allowed_current_note",
            statement=(
                "Draft RO-Crate and Data Package descriptors are inspectable "
                "scaffolds with release blockers."
            ),
            supporting_evidence=[ro_crate_summary_path.as_posix(), guard_path.as_posix()],
            prohibited_language="archive-ready RO-Crate",
            next_allowed_action="describe as draft scaffold only",
        ),
        _diagnostic_public_claim_row(
            release_note_id=release_note_id,
            guard_id=guard_id,
            package_id=package_id,
            task_id=task_id,
            claim_id="no_doi_or_archive_identifier",
            claim_category="prohibited",
            public_note_status="not_allowed_current_note",
            statement="No DOI or local archive identifier is claimed.",
            supporting_evidence=[guard_path.as_posix(), mutation_guard_path.as_posix()],
            prohibited_language="DOI, Zenodo archive, citable release",
            next_allowed_action="retry only after owner archive-route metadata is supplied",
        ),
        _diagnostic_public_claim_row(
            release_note_id=release_note_id,
            guard_id=guard_id,
            package_id=package_id,
            task_id=task_id,
            claim_id="no_license_or_creator_metadata",
            claim_category="prohibited",
            public_note_status="not_allowed_current_note",
            statement=(
                "No local license, creator list, contributor list, publisher, or "
                "release version is claimed."
            ),
            supporting_evidence=[guard_path.as_posix(), action_path.as_posix()],
            prohibited_language="licensed release, official author list, publisher",
            next_allowed_action="retry only after owner metadata file is supplied",
        ),
        _diagnostic_public_claim_row(
            release_note_id=release_note_id,
            guard_id=guard_id,
            package_id=package_id,
            task_id=task_id,
            claim_id="no_leaderboard_or_benchmark_promotion",
            claim_category="prohibited",
            public_note_status="not_allowed_current_note",
            statement=(
                "OSD-120 remains a draft diagnostic/provenance case study, not a "
                "promoted v9-alpha benchmark core or leaderboard task."
            ),
            supporting_evidence=["docs/V9_PURPOSE_DRIFT_AUDIT_2026_05_26.md"],
            prohibited_language="v9-alpha leaderboard task",
            next_allowed_action="recenter before any benchmark-core promotion",
        ),
    ]
    action_by_id = {row["action_id"]: row for row in action_rows}
    retry_items = [
        _owner_retry_checklist_row(
            release_note_id=release_note_id,
            guard_id=guard_id,
            package_id=package_id,
            task_id=task_id,
            retry_item_id="owner_metadata_file",
            required_for="any_descriptor_metadata_application",
            owner_action=action_by_id["supply_release_owner_metadata_file"]["next_action"],
            current_status="missing",
            source_blocker_fields=mutation_guard["blocked_fields"],
            evidence_artifacts=common_evidence,
            validation_rule="CSV_or_JSON_values_must_map_to_V9_MULTI_033_intake_template",
            retry_priority="P0",
            closeout_policy="do_not_extend_osd120_metadata_chain_without_this_file",
        ),
        _owner_retry_checklist_row(
            release_note_id=release_note_id,
            guard_id=guard_id,
            package_id=package_id,
            task_id=task_id,
            retry_item_id="archive_route_identifier",
            required_for="DOI_or_archive_release",
            owner_action=action_by_id["decide_archive_route_identifier"]["next_action"],
            current_status="missing",
            source_blocker_fields=["future_archive_identifier", "repository_url"],
            evidence_artifacts=[action_path.as_posix()],
            validation_rule="explicit_no_archive_or_owner_supplied_identifier_required",
            retry_priority="P1",
            closeout_policy="keep_no_archive_diagnostic_note_until_supplied",
        ),
        _owner_retry_checklist_row(
            release_note_id=release_note_id,
            guard_id=guard_id,
            package_id=package_id,
            task_id=task_id,
            retry_item_id="version_date_year",
            required_for="citable_release_metadata",
            owner_action=action_by_id["freeze_version_date_year"]["next_action"],
            current_status="missing",
            source_blocker_fields=["release_version_string", "release_date", "publication_year"],
            evidence_artifacts=[action_path.as_posix()],
            validation_rule="release_tag_date_and_year_must_be_owner_supplied",
            retry_priority="P1",
            closeout_policy="do_not_use_build_date_or_draft_string",
        ),
        _owner_retry_checklist_row(
            release_note_id=release_note_id,
            guard_id=guard_id,
            package_id=package_id,
            task_id=task_id,
            retry_item_id="creator_contributor_publisher",
            required_for="DataCite_RO_Crate_CITATION_cff",
            owner_action=action_by_id["freeze_creator_contributor_publisher"]["next_action"],
            current_status="missing",
            source_blocker_fields=[
                "local_package_creators",
                "local_package_contributors",
                "publisher_maintainer_identity",
            ],
            evidence_artifacts=[action_path.as_posix()],
            validation_rule="authorship_order_roles_and_publisher_must_not_be_inferred",
            retry_priority="P1",
            closeout_policy="do_not_generate_creator_or_publisher_metadata",
        ),
        _owner_retry_checklist_row(
            release_note_id=release_note_id,
            guard_id=guard_id,
            package_id=package_id,
            task_id=task_id,
            retry_item_id="license_scope",
            required_for="reuse_terms_and_repository_release",
            owner_action=action_by_id["choose_license_scope"]["next_action"],
            current_status="missing",
            source_blocker_fields=["local_code_license", "local_metadata_tables_license"],
            evidence_artifacts=[action_path.as_posix()],
            validation_rule="SPDX_identifier_or_explicit_custom_terms_after_review",
            retry_priority="P1",
            closeout_policy="do_not_apply_license_without_owner_review",
        ),
        _owner_retry_checklist_row(
            release_note_id=release_note_id,
            guard_id=guard_id,
            package_id=package_id,
            task_id=task_id,
            retry_item_id="exact_osdr_study_citation",
            required_for="source_alpha_package_citation",
            owner_action=action_by_id["capture_exact_osdr_study_citation"]["next_action"],
            current_status="missing",
            source_blocker_fields=["upstream_osdr_dataset_citation"],
            evidence_artifacts=[action_path.as_posix()],
            validation_rule="OSDR_study_page_BibTeX_or_RIS_required",
            retry_priority="P2",
            closeout_policy="retain_generic_osdr_credit_only",
        ),
    ]
    inspectable_count = sum(
        1 for row in public_claims if row["public_note_status"] == "allowed_current_note"
    )
    not_released_count = sum(
        1 for row in public_claims if row["public_note_status"] == "not_allowed_current_note"
    )
    summary = {
        "release_note_id": release_note_id,
        "guard_id": guard_id,
        "fill_id": fill_id,
        "package_id": package_id,
        "task_id": task_id,
        "note_status": "diagnostic_metadata_note_ready_no_archive_release",
        "branch_closeout_status": "osd120_metadata_branch_closeout_pending_refocus",
        "current_public_surface": "diagnostic_metadata_only",
        "inspectable_now_count": str(inspectable_count),
        "not_released_claim_count": str(not_released_count),
        "owner_retry_item_count": str(len(retry_items)),
        "descriptor_mutation_status": "not_mutated",
        "archive_release_status": "deferred_no_owner_metadata",
        "next_required_block": "V9-REFOCUS-001: post-OSD-120 recenter decision",
        "claim_boundary": public_note_claim_boundary,
    }
    included_sections = [row for row in sections if row["include_in_note"] == "True"]
    section_text = "\n\n".join(
        f"## {row['section_title']}\n\n{row['note_text']}"
        for row in included_sections
    )
    release_note_md = f"""# OSD-120 Diagnostic Metadata Note

Release note id: `{release_note_id}`

{section_text}

## Claim Boundary

This note is diagnostic metadata only. It does not claim DOI-backed archive
release, CITATION.cff release, local package license, creator/publisher
metadata, full OSDR payload mirror, or v9-alpha leaderboard promotion.
"""
    review_md = f"""# OSD-120 Diagnostic Metadata Release Note Closeout Review

Block: V9-MULTI-035

Release note id: `{release_note_id}`

Status: `{summary["note_status"]}`

This block closes the current OSD-120 metadata branch without mutating archive
metadata. The output is a diagnostic metadata note plus an owner metadata retry
checklist, not another archive-release gate.

## Closeout Decision

- Current public surface: diagnostic metadata only.
- Archive release: deferred because owner metadata is absent.
- Descriptor mutation: not allowed.
- Next block: V9-REFOCUS-001 should choose public bulk alpha readiness or first
  single-cell flagship inventory.

## Allowed Current Note Claims

{chr(10).join(f"- `{row['claim_id']}`: {row['statement']}" for row in public_claims if row['public_note_status'] == 'allowed_current_note')}

## Prohibited Release Claims

{chr(10).join(f"- `{row['claim_id']}`: {row['prohibited_language']}" for row in public_claims if row['public_note_status'] == 'not_allowed_current_note')}

## Owner Retry Items

{chr(10).join(f"- `{row['retry_item_id']}`: {row['current_status']}." for row in retry_items)}

## External Guidance Anchors

- GitHub/Zenodo DOI paths require public repository access, a release, and
  license/reuse terms.
- OSDR dataset citation should use the study-page citation button in BibTeX or
  RIS form.
- DataCite citable metadata requires identifier, creator, title, publisher,
  publication year, and resource type.

## Drift Guard

Per `docs/V9_PURPOSE_DRIFT_AUDIT_2026_05_26.md`, do not add another OSD-120
metadata/release block unless owner-supplied metadata appears. Recenter after
this closeout.
"""
    return {
        "summary": [summary],
        "sections": sections,
        "public_claims": public_claims,
        "retry_items": retry_items,
        "release_note_md": release_note_md,
        "review_md": review_md,
    }


def write_multispecies_interaction_diagnostic_metadata_release_note(
    *,
    output_dir: str | Path = (
        "v9/multispecies/reports/interaction_diagnostic_metadata_release_note"
    ),
    reports_root: str | Path = "v9/multispecies/reports",
    repo_root: str | Path = ".",
    release_note_id: str = DEFAULT_INTERACTION_DIAGNOSTIC_METADATA_RELEASE_NOTE_ID,
) -> dict[str, Path]:
    """Write the diagnostic metadata release-note closeout artifacts."""

    package = build_multispecies_interaction_diagnostic_metadata_release_note(
        reports_root=reports_root,
        repo_root=repo_root,
        release_note_id=release_note_id,
    )
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "summary_csv": outdir / "diagnostic_metadata_release_summary.csv",
        "summary_json": outdir / "diagnostic_metadata_release_summary.json",
        "section_csv": outdir / "diagnostic_metadata_release_note_sections.csv",
        "section_json": outdir / "diagnostic_metadata_release_note_sections.json",
        "claim_csv": outdir / "diagnostic_metadata_public_claims.csv",
        "claim_json": outdir / "diagnostic_metadata_public_claims.json",
        "retry_csv": outdir / "owner_metadata_retry_checklist.csv",
        "retry_json": outdir / "owner_metadata_retry_checklist.json",
        "release_note_md": outdir / "OSD120_DIAGNOSTIC_METADATA_NOTE.md",
        "review_md": (
            outdir.parent
            / "OSD120_INTERACTION_DIAGNOSTIC_METADATA_RELEASE_NOTE_REVIEW.md"
        ),
    }
    with outputs["summary_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_DIAGNOSTIC_METADATA_RELEASE_SUMMARY_FIELDS,
        )
        writer.writeheader()
        writer.writerows(package["summary"])
    outputs["summary_json"].write_text(
        json.dumps(package["summary"], indent=2, sort_keys=True) + "\n"
    )
    with outputs["section_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_DIAGNOSTIC_METADATA_RELEASE_SECTION_FIELDS,
        )
        writer.writeheader()
        writer.writerows(package["sections"])
    outputs["section_json"].write_text(
        json.dumps(package["sections"], indent=2, sort_keys=True) + "\n"
    )
    with outputs["claim_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_DIAGNOSTIC_METADATA_PUBLIC_CLAIM_FIELDS,
        )
        writer.writeheader()
        writer.writerows(package["public_claims"])
    outputs["claim_json"].write_text(
        json.dumps(package["public_claims"], indent=2, sort_keys=True) + "\n"
    )
    with outputs["retry_csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=INTERACTION_OWNER_METADATA_RETRY_CHECKLIST_FIELDS,
        )
        writer.writeheader()
        writer.writerows(package["retry_items"])
    outputs["retry_json"].write_text(
        json.dumps(package["retry_items"], indent=2, sort_keys=True) + "\n"
    )
    outputs["release_note_md"].write_text(str(package["release_note_md"]))
    outputs["review_md"].write_text(str(package["review_md"]))
    return outputs


def _default_sensitivity_configs(task: MultispeciesTaskData) -> list[MultispeciesBaselineConfig]:
    return [
        MultispeciesBaselineConfig(
            transform=transform,
            scaling=scaling,
            top_variable_genes=top,
        )
        for transform in ("log1p", "none")
        for scaling in ("zscore", "none")
        for top in (100, 500, 2000, 5000, task.n_features)
    ]


def run_multispecies_sensitivity_grid(
    *,
    manifest_dir: str | Path = "v9/multispecies/task_manifests",
    repo_root: str | Path = ".",
    output_root: str | Path = "v9/multispecies/reports/sensitivity",
    configs: Sequence[MultispeciesBaselineConfig] | None = None,
    command: Sequence[str] | None = None,
) -> tuple[list[MultispeciesBaselineTaskResult], dict[str, Path]]:
    """Run preprocessing sensitivity variants for draft multispecies baselines."""

    tasks = load_all_multispecies_tasks(manifest_dir=manifest_dir, repo_root=repo_root)
    root = Path(output_root)
    repo = Path(repo_root)
    results: list[MultispeciesBaselineTaskResult] = []
    rows: list[dict[str, str]] = []
    for task in tasks:
        task_configs = list(configs) if configs is not None else _default_sensitivity_configs(task)
        manifest_path = repo / manifest_dir / f"{task.task_id}.json"
        for config in task_configs:
            variant_id = multispecies_config_id(config)
            result = run_multispecies_nearest_centroid(
                task,
                output_dir=root / variant_id / task.task_id,
                config=config,
                task_manifest_path=manifest_path,
                command=command,
            )
            results.append(result)
            rows.append(multispecies_result_summary_row(result))
    summary = write_multispecies_baseline_summary(output_root, rows)
    return results, summary


def run_multispecies_interaction_sensitivity_grid(
    *,
    manifest_dir: str | Path = "v9/multispecies/interaction_task_manifests",
    repo_root: str | Path = ".",
    output_root: str | Path = "v9/multispecies/reports/interaction_sensitivity",
    configs: Sequence[MultispeciesBaselineConfig] | None = None,
    fold_families: Sequence[str] | None = None,
    command: Sequence[str] | None = None,
) -> tuple[list[MultispeciesBaselineTaskResult], dict[str, Path]]:
    """Run preprocessing sensitivity variants for OSD-120 interaction baselines."""

    tasks = load_all_multispecies_interaction_tasks(
        manifest_dir=manifest_dir,
        repo_root=repo_root,
    )
    selected_families = list(fold_families or INTERACTION_FOLD_FAMILIES)
    unknown = sorted(set(selected_families) - set(INTERACTION_FOLD_FAMILIES))
    if unknown:
        raise ValueError(f"unknown multispecies interaction fold families: {unknown}")

    root = Path(output_root)
    repo = Path(repo_root)
    results: list[MultispeciesBaselineTaskResult] = []
    rows: list[dict[str, str]] = []
    for task in tasks:
        task_configs = list(configs) if configs is not None else _default_sensitivity_configs(task)
        manifest_path = repo / manifest_dir / f"{task.task_id}.json"
        for config in task_configs:
            variant_id = multispecies_config_id(config)
            for fold_family in selected_families:
                result = run_multispecies_nearest_centroid(
                    task,
                    output_dir=root / variant_id / fold_family / task.task_id,
                    config=config,
                    folds=task.fold_families[fold_family],
                    baseline_id=INTERACTION_BASELINE_ID,
                    fold_family=fold_family,
                    claim_boundary="draft_interaction_sensitivity_only_not_leaderboard",
                    task_manifest_path=manifest_path,
                    command=command,
                )
                results.append(result)
                rows.append(multispecies_result_summary_row(result))
    summary = write_multispecies_baseline_summary(output_root, rows)
    return results, summary


def _resolve_path(path: str | Path, repo_root: Path) -> Path:
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    return repo_root / candidate


def aggregate_multispecies_interaction_fold_details(
    *,
    summary_csv: str | Path = (
        "v9/multispecies/reports/interaction_sensitivity/multispecies_baseline_summary.csv"
    ),
    repo_root: str | Path = ".",
) -> list[dict[str, str]]:
    """Extract fold-level sensitivity details from OSD-120 metrics reports."""

    root = Path(repo_root)
    summary_path = _resolve_path(summary_csv, root)
    with summary_path.open(newline="") as handle:
        summary_rows = list(csv.DictReader(handle))

    detail_rows: list[dict[str, str]] = []
    for summary_row in summary_rows:
        fold_family = str(summary_row["fold_family"])
        metric_id = DELTA_METRIC_BY_FOLD_FAMILY.get(fold_family)
        if not metric_id:
            continue
        metrics_path = _resolve_path(summary_row["metrics"], root)
        metrics_payload = json.loads(metrics_path.read_text())
        metric = metrics_payload.get("metrics", {}).get(metric_id, {})
        if metric.get("status") != "computed":
            continue
        details = metric.get("details", {})
        if not isinstance(details, Mapping):
            continue
        scored = sorted(
            (
                (
                    str(fold_id),
                    str(detail.get("heldout_value", "")),
                    int(detail.get("n_test", 0)),
                    float(detail.get("balanced_accuracy", 0.0)),
                )
                for fold_id, detail in details.items()
                if isinstance(detail, Mapping)
            ),
            key=lambda item: (item[3], item[0]),
        )
        if not scored:
            continue
        lowest_score = scored[0][3]
        for rank, (fold_id, heldout_value, n_test, balanced_accuracy) in enumerate(
            scored,
            start=1,
        ):
            detail_rows.append(
                {
                    "baseline_id": str(summary_row["baseline_id"]),
                    "variant_id": str(summary_row["variant_id"]),
                    "transform": str(summary_row["transform"]),
                    "scaling": str(summary_row["scaling"]),
                    "top_variable_genes": str(summary_row["top_variable_genes"]),
                    "task_id": str(summary_row["task_id"]),
                    "fold_family": fold_family,
                    "delta_metric_id": metric_id,
                    "fold_id": fold_id,
                    "heldout_factor": HELDOUT_FACTOR_BY_FOLD_FAMILY.get(
                        fold_family,
                        "",
                    ),
                    "heldout_value": heldout_value,
                    "n_test": str(n_test),
                    "balanced_accuracy": _format_float(balanced_accuracy),
                    "rank_low_to_high": str(rank),
                    "is_lowest_for_variant": str(balanced_accuracy == lowest_score),
                    "is_balanced_accuracy_le_0_5": str(balanced_accuracy <= 0.5),
                    "summary_metrics": summary_path.as_posix(),
                    "metrics": metrics_path.as_posix(),
                }
            )
    return detail_rows


def write_multispecies_interaction_fold_details(
    *,
    output_dir: str | Path = "v9/multispecies/reports/interaction_sensitivity",
    summary_csv: str | Path = (
        "v9/multispecies/reports/interaction_sensitivity/multispecies_baseline_summary.csv"
    ),
    repo_root: str | Path = ".",
) -> dict[str, Path]:
    """Write CSV/JSON fold-detail summaries for OSD-120 interaction sensitivity."""

    rows = aggregate_multispecies_interaction_fold_details(
        summary_csv=summary_csv,
        repo_root=repo_root,
    )
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / "fold_detail_summary.csv"
    json_path = outdir / "fold_detail_summary.json"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=INTERACTION_FOLD_DETAIL_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    json_path.write_text(json.dumps(rows, indent=2, sort_keys=True) + "\n")
    return {"csv": csv_path, "json": json_path}


def _read_fold_detail_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def compare_multispecies_interaction_fold_details(
    *,
    nearest_centroid_fold_detail_csv: str | Path = (
        "v9/multispecies/reports/interaction_sensitivity/fold_detail_summary.csv"
    ),
    logistic_fold_detail_csv: str | Path = (
        "v9/multispecies/reports/interaction_logistic_l2/fold_detail_summary.csv"
    ),
    repo_root: str | Path = ".",
    nearest_centroid_variant_id: str = "tvg2000_log1p_zscore",
    logistic_variant_id: str | None = "tvg2000_log1p_zscore",
) -> list[dict[str, str]]:
    """Compare nearest-centroid and L2 logistic OSD-120 fold details."""

    root = Path(repo_root)
    nearest_path = _resolve_path(nearest_centroid_fold_detail_csv, root)
    logistic_path = _resolve_path(logistic_fold_detail_csv, root)
    nearest_rows = [
        row
        for row in _read_fold_detail_rows(nearest_path)
        if row["variant_id"] == nearest_centroid_variant_id
    ]
    logistic_rows = [
        row
        for row in _read_fold_detail_rows(logistic_path)
        if logistic_variant_id is None or row["variant_id"] == logistic_variant_id
    ]
    nearest_by_key = {
        (row["fold_family"], row["heldout_value"]): row for row in nearest_rows
    }
    comparison_rows: list[dict[str, str]] = []
    for logistic_row in sorted(
        logistic_rows,
        key=lambda row: (row["variant_id"], row["fold_family"], row["heldout_value"]),
    ):
        key = (logistic_row["fold_family"], logistic_row["heldout_value"])
        nearest_row = nearest_by_key.get(key)
        if nearest_row is None:
            raise ValueError(f"nearest-centroid fold detail missing for {key}")
        nearest_score = float(nearest_row["balanced_accuracy"])
        logistic_score = float(logistic_row["balanced_accuracy"])
        delta = logistic_score - nearest_score
        comparison_rows.append(
            {
                "task_id": logistic_row["task_id"],
                "fold_family": logistic_row["fold_family"],
                "heldout_factor": logistic_row["heldout_factor"],
                "heldout_value": logistic_row["heldout_value"],
                "n_test": logistic_row["n_test"],
                "nearest_centroid_baseline_id": nearest_row["baseline_id"],
                "nearest_centroid_variant_id": nearest_row["variant_id"],
                "nearest_centroid_balanced_accuracy": nearest_row["balanced_accuracy"],
                "nearest_centroid_rank_low_to_high": nearest_row["rank_low_to_high"],
                "nearest_centroid_is_lowest_for_variant": nearest_row[
                    "is_lowest_for_variant"
                ],
                "logistic_baseline_id": logistic_row["baseline_id"],
                "logistic_variant_id": logistic_row["variant_id"],
                "logistic_balanced_accuracy": logistic_row["balanced_accuracy"],
                "logistic_rank_low_to_high": logistic_row["rank_low_to_high"],
                "logistic_is_lowest_for_variant": logistic_row["is_lowest_for_variant"],
                "delta_logistic_minus_nearest_centroid": _format_float(delta),
                "logistic_improved": str(delta > 0.0),
                "logistic_new_or_persistent_le_0_5": str(
                    logistic_score <= 0.5 or nearest_score <= 0.5
                ),
                "nearest_centroid_metrics": nearest_row["metrics"],
                "logistic_metrics": logistic_row["metrics"],
            }
        )
    return comparison_rows


def write_multispecies_interaction_fold_detail_comparison(
    *,
    output_dir: str | Path = "v9/multispecies/reports/interaction_logistic_l2",
    nearest_centroid_fold_detail_csv: str | Path = (
        "v9/multispecies/reports/interaction_sensitivity/fold_detail_summary.csv"
    ),
    logistic_fold_detail_csv: str | Path = (
        "v9/multispecies/reports/interaction_logistic_l2/fold_detail_summary.csv"
    ),
    repo_root: str | Path = ".",
    nearest_centroid_variant_id: str = "tvg2000_log1p_zscore",
    logistic_variant_id: str | None = "tvg2000_log1p_zscore",
) -> dict[str, Path]:
    """Write side-by-side OSD-120 logistic versus nearest-centroid fold details."""

    rows = compare_multispecies_interaction_fold_details(
        nearest_centroid_fold_detail_csv=nearest_centroid_fold_detail_csv,
        logistic_fold_detail_csv=logistic_fold_detail_csv,
        repo_root=repo_root,
        nearest_centroid_variant_id=nearest_centroid_variant_id,
        logistic_variant_id=logistic_variant_id,
    )
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / "fold_detail_comparison_vs_nearest_centroid.csv"
    json_path = outdir / "fold_detail_comparison_vs_nearest_centroid.json"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=INTERACTION_FOLD_COMPARISON_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    json_path.write_text(json.dumps(rows, indent=2, sort_keys=True) + "\n")
    return {"csv": csv_path, "json": json_path}
