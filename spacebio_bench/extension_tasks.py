"""Draft task manifests for v9 extension tracks."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Mapping, Sequence

from .extension_sources import HUMAN_ORGANOID_DRAFT_SOURCES, MULTISPECIES_DRAFT_SOURCES
from .manifests import validate_task_manifest
from .organoid_de_reference import (
    BENCHMARK_LOG2FC_ORIENTATION,
    RESPONSE_SIGNATURE_REQUIRED_COLUMNS,
)


HUMAN_ORGANOID_TASK_ID = "draft_human_organoid_spaceflight"
MULTISPECIES_ARABIDOPSIS_TASK_ID = "draft_osd37_arabidopsis_seedling_spaceflight"
MULTISPECIES_DROSOPHILA_TASK_ID = "draft_osd207_drosophila_whole_body_spaceflight"
MULTISPECIES_OSD120_INTERACTION_TASK_ID = (
    "draft_osd120_arabidopsis_root_light_interaction_spaceflight"
)
MULTISPECIES_TASK_IDS = [
    MULTISPECIES_ARABIDOPSIS_TASK_ID,
    MULTISPECIES_DROSOPHILA_TASK_ID,
]

MULTISPECIES_TASK_CONFIGS = {
    "OSD-37": {
        "task_id": MULTISPECIES_ARABIDOPSIS_TASK_ID,
        "title": "Arabidopsis seedling-pool species-native LEO vs ground draft task",
        "split_name": "arabidopsis_ecotype_blocked_leo_vs_ground_pilot",
        "primary_blocking_factor": "genotype_or_ecotype",
        "defer_note": "",
    },
    "OSD-207": {
        "task_id": MULTISPECIES_DROSOPHILA_TASK_ID,
        "title": "Drosophila whole-body species-native LEO vs ground draft task",
        "split_name": "drosophila_background_variant_blocked_leo_vs_ground_pilot",
        "primary_blocking_factor": "condition_stratum",
        "defer_note": "",
    },
}


TASK_SOURCE_FIELDS = [
    "source_id",
    "url_or_accession",
    "access_status",
    "checksum_status",
    "privacy_class",
    "glds_prefix",
    "mission",
    "notes",
    "organism",
    "taxon_id",
    "species_common_name",
    "material_type",
    "model_system",
    "biospecimen_type",
    "assay_modality",
    "platform",
    "spaceflight_environment",
    "ground_control_type",
    "donor_or_strain_block",
    "orthology_strategy",
    "feature_namespace",
    "external_accessions",
    "publication_urls",
    "organoid_type",
    "differentiation_protocol",
    "microglia_condition",
    "disease_context",
    "culture_hardware",
    "processed_file_status",
    "integration_status",
]


def _join_unique(rows: Sequence[Mapping[str, str]], key: str) -> str:
    values = sorted({str(row.get(key, "") or "") for row in rows if row.get(key)})
    return ";".join(values)


def _source_record(row: Mapping[str, str]) -> dict[str, str]:
    record = {field: str(row.get(field, "") or "") for field in TASK_SOURCE_FIELDS}
    record["source_url"] = str(row.get("osd_url", "") or "")
    return record


def _slug(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "_", value.strip().lower()).strip("_")
    return slug or "unknown"


def _count_by(rows: Sequence[Mapping[str, str]], key: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        value = str(row.get(key, "") or "")
        if not value:
            continue
        counts[value] = counts.get(value, 0) + 1
    return dict(sorted(counts.items()))


def _label_distribution(rows: Sequence[Mapping[str, str]]) -> dict[str, int]:
    return _count_by(rows, "true_label")


def _folds_for_factor(
    rows: Sequence[Mapping[str, str]],
    factor: str,
    *,
    prefix: str,
) -> list[dict[str, object]]:
    values = sorted({str(row.get(factor, "") or "") for row in rows if row.get(factor)})
    folds: list[dict[str, object]] = []
    for value in values:
        test_rows = [row for row in rows if row.get(factor) == value]
        train_rows = [row for row in rows if row.get(factor) != value]
        if not test_rows or not train_rows:
            continue
        folds.append(
            {
                "fold_id": f"{prefix}_{_slug(value)}",
                "heldout_factor": factor,
                "heldout_value": value,
                "n_train": len(train_rows),
                "n_test": len(test_rows),
                "train_label_distribution": _label_distribution(train_rows),
                "test_label_distribution": _label_distribution(test_rows),
                "status": "sample_count_backed_draft",
            }
        )
    return folds


def _sample_count_summary(rows: Sequence[Mapping[str, str]]) -> dict[str, object]:
    donor_distribution = _count_by(rows, "donor_or_line_id")
    ipsc_line_distribution = _count_by(rows, "ipsc_line_id")
    return {
        "n_samples": len(rows),
        "source_sample_rows": _count_by(rows, "source_id"),
        "label_distribution": _label_distribution(rows),
        "organoid_type_distribution": _count_by(rows, "organoid_type"),
        "microglia_condition_distribution": _count_by(rows, "microglia_condition"),
        "disease_context_distribution": _count_by(rows, "disease_context"),
        "donor_or_line_distribution": donor_distribution,
        "ipsc_line_distribution": ipsc_line_distribution,
    }


def _sample_backed_split(rows: Sequence[Mapping[str, str]]) -> dict[str, object]:
    parsed_rows = [row for row in rows if row.get("parse_status") == "parsed"]
    if not parsed_rows:
        return {}
    donor_distribution = _count_by(parsed_rows, "donor_or_line_id")
    donor_metadata_available = bool(donor_distribution)
    folds = [
        *_folds_for_factor(
            parsed_rows,
            "organoid_type",
            prefix="holdout_organoid_type",
        ),
        *_folds_for_factor(
            parsed_rows,
            "microglia_condition",
            prefix="holdout_microglia_condition",
        ),
    ]
    donor_diagnostic_folds = _folds_for_factor(
        parsed_rows,
        "donor_or_line_id",
        prefix="holdout_donor_or_line",
    )
    for fold in donor_diagnostic_folds:
        fold["status"] = "metadata_backed_pilot_only_not_default"
        fold["confounding_note"] = (
            "Donor/iPSC-line identifiers are parsed, but donor is not independently "
            "crossed with source, organoid fate, or disease context."
        )
    return {
        **_sample_count_summary(parsed_rows),
        "sample_factor_table": "v9/human_organoid/sample_factors.draft.csv",
        "sample_factor_status": (
            "condition_factors_and_geo_donor_line_metadata_parsed"
            if donor_metadata_available
            else "condition_factors_parsed"
        ),
        "donor_metadata_status": (
            "parsed_from_geo_series_matrix"
            if donor_metadata_available
            else "not_available_in_osdr_sample_table"
        ),
        "donor_blocking_status": (
            "available_as_metadata_but_not_default_fold_due_source_organoid_disease_confounding"
            if donor_metadata_available
            else "donor_block_named_but_unresolved"
        ),
        "candidate_folds": folds,
        "donor_diagnostic_folds": donor_diagnostic_folds,
    }


def _expression_matrix_summary(rows: Sequence[Mapping[str, str]]) -> dict[str, object]:
    aligned_rows = [
        row
        for row in rows
        if row.get("audit_status") == "matrix_downloaded_sample_aligned"
    ]
    if not aligned_rows:
        return {}
    return {
        "expression_matrix_audit": "v9/human_organoid/expression_matrix_audit.draft.csv",
        "expression_matrix_status": "matrix_downloaded_sample_aligned",
        "expression_matrix_sources": {
            str(row.get("source_id", "")): {
                "matrix_file": str(row.get("matrix_files", "")),
                "local_matrix_path": str(row.get("local_matrix_paths", "")),
                "n_feature_rows": int(row.get("n_feature_rows", "0") or 0),
                "n_sample_columns": int(row.get("n_sample_columns", "0") or 0),
                "matching_sample_count": int(row.get("matching_sample_count", "0") or 0),
                "sha256": str(row.get("matrix_sha256", "")),
            }
            for row in aligned_rows
            if row.get("source_id")
        },
    }


def _multispecies_sample_summary(rows: Sequence[Mapping[str, str]]) -> dict[str, object]:
    return {
        "n_samples": len(rows),
        "source_sample_rows": _count_by(rows, "source_id"),
        "label_distribution": _label_distribution(rows),
        "condition_stratum_distribution": _count_by(rows, "condition_stratum"),
        "genotype_or_ecotype_distribution": _count_by(rows, "genotype_or_ecotype"),
        "light_treatment_distribution": _count_by(rows, "light_treatment"),
    }


def _multispecies_expression_summary(
    rows: Sequence[Mapping[str, str]],
    source_id: str,
) -> dict[str, object]:
    source_rows = [row for row in rows if row.get("source_id") == source_id]
    aligned_rows = [
        row
        for row in source_rows
        if row.get("audit_status") == "matrix_local_sample_aligned"
        and row.get("matrix_columns_match_sample_factors") == "True"
    ]
    if not aligned_rows:
        return {}
    row = aligned_rows[0]
    return {
        "expression_matrix_audit": "v9/multispecies/expression_matrix_audit.draft.csv",
        "expression_matrix_status": "matrix_local_sample_aligned",
        "expression_matrix_sources": {
            source_id: {
                "local_matrix_path": str(row.get("local_matrix_path", "")),
                "n_feature_rows": int(row.get("n_feature_rows", "0") or 0),
                "n_sample_columns": int(row.get("n_sample_columns", "0") or 0),
                "n_matching_sample_columns": int(
                    row.get("n_matching_sample_columns", "0") or 0
                ),
                "matrix_columns_match_sample_factors": str(
                    row.get("matrix_columns_match_sample_factors", "")
                ),
            }
        },
    }


def _multispecies_checksum_summary(
    rows: Sequence[Mapping[str, str]],
    source_id: str,
) -> dict[str, object]:
    source_rows = [row for row in rows if row.get("source_id") == source_id]
    if not source_rows:
        return {}
    row = source_rows[0]
    return {
        "source_checksum_audit": "v9/multispecies/source_checksum_audit.draft.csv",
        "source_checksum_status": str(row.get("audit_status", "") or ""),
        "source_checksum_sources": {
            source_id: {
                "api_status": str(row.get("api_status", "") or ""),
                "n_files": int(row.get("n_files", "0") or 0),
                "checksum_manifest_count": int(
                    row.get("checksum_manifest_count", "0") or 0
                ),
                "parsed_checksum_entries": int(
                    row.get("parsed_checksum_entries", "0") or 0
                ),
                "checksum_payload_matches": int(
                    row.get("checksum_payload_matches", "0") or 0
                ),
                "freeze_ready": str(row.get("freeze_ready", "") or ""),
            }
        },
    }


def build_multispecies_task_manifests(
    rows: Sequence[Mapping[str, str]] | None = None,
    *,
    sample_factor_rows: Sequence[Mapping[str, str]] | None = None,
    expression_matrix_audit_rows: Sequence[Mapping[str, str]] | None = None,
    source_checksum_audit_rows: Sequence[Mapping[str, str]] | None = None,
) -> list[dict[str, object]]:
    """Build draft species-native multispecies task manifests for OSD-37 and OSD-207."""

    source_rows = list(rows or MULTISPECIES_DRAFT_SOURCES)
    source_by_id = {str(row.get("source_id", "") or ""): row for row in source_rows}
    sample_rows = list(sample_factor_rows or [])
    matrix_rows = list(expression_matrix_audit_rows or [])
    checksum_rows = list(source_checksum_audit_rows or [])
    manifests: list[dict[str, object]] = []

    for source_id in ("OSD-37", "OSD-207"):
        source_row = source_by_id.get(source_id)
        if source_row is None:
            continue
        config = MULTISPECIES_TASK_CONFIGS[source_id]
        task_sample_rows = [
            row
            for row in sample_rows
            if row.get("source_id") == source_id and row.get("parse_status") == "parsed"
        ]
        source_record = _source_record(source_row)
        checksum_summary = _multispecies_checksum_summary(checksum_rows, source_id)
        if checksum_summary:
            source_record["checksum_status"] = str(
                checksum_summary.get("source_checksum_status", "")
            )
            source_record["local_payload_status"] = "local_sample_table_and_matrix_md5_matched"
        sample_summary = _multispecies_sample_summary(task_sample_rows)
        expression_summary = _multispecies_expression_summary(matrix_rows, source_id)
        candidate_folds = _folds_for_factor(
            task_sample_rows,
            "condition_stratum",
            prefix="holdout_condition_stratum",
        )
        for fold in candidate_folds:
            fold["status"] = "sample_count_backed_candidate_not_default"
            fold["interpretation"] = (
                "Condition-stratum holdout checks robustness across genotype/ecotype "
                "context, not mission-held-out generalization."
            )
        blocking_factors = []
        for factor in (str(config["primary_blocking_factor"]), "condition_stratum"):
            if factor and factor not in blocking_factors:
                blocking_factors.append(factor)

        manifest: dict[str, object] = {
            "schema_version": "0.1.0",
            "task_id": str(config["task_id"]),
            "task_family": "multispecies_species_native_spaceflight",
            "title": str(config["title"]),
            "release_status": "draft_not_frozen",
            "organism": str(source_row.get("organism", "") or ""),
            "taxon_id": str(source_row.get("taxon_id", "") or ""),
            "species_common_name": str(source_row.get("species_common_name", "") or ""),
            "material_type": str(source_row.get("material_type", "") or ""),
            "model_system": str(source_row.get("model_system", "") or ""),
            "biospecimen_type": str(source_row.get("biospecimen_type", "") or ""),
            "assay_modality": str(source_row.get("assay_modality", "") or ""),
            "platform": str(source_row.get("platform", "") or ""),
            "spaceflight_environment": str(
                source_row.get("spaceflight_environment", "") or ""
            ),
            "ground_control_type": str(source_row.get("ground_control_type", "") or ""),
            "donor_or_strain_block": str(source_row.get("donor_or_strain_block", "") or ""),
            "orthology_strategy": "not_applicable_species_native_raw_gene_task",
            "feature_namespace": str(source_row.get("feature_namespace", "") or ""),
            "variant": str(source_row.get("variants", "") or ""),
            "source_records": [source_record],
            "split": {
                "name": str(config["split_name"]),
                "unit": "sample",
                "strategy": (
                    "Classify LEO/spaceflight versus ground-control samples within "
                    "one non-mouse OSDR source, with species-native features and "
                    "explicit genotype/ecotype condition strata."
                ),
                "status": "draft_sample_matrix_checksum_backed_pending_baseline",
                "target_labels": ["LEO_or_ISS", "Ground"],
                "n_sources": 1,
                "source_accessions": [source_id],
                "sample_factor_table": "v9/multispecies/sample_factors.draft.csv",
                "sample_factor_status": "condition_factors_parsed",
                "blocking_factors": blocking_factors,
                "candidate_folds": candidate_folds,
                "single_source_limitation": (
                    "This draft is a species-native within-source task, not a "
                    "leave-one-mission-out or cross-species generalization benchmark."
                ),
                **sample_summary,
                **expression_summary,
            },
            "metrics": [
                {
                    "metric_id": "balanced_accuracy",
                    "profile": "genelab_multispecies_pilot",
                    "interpretation": "Class-balance-aware LEO/spaceflight versus Ground score.",
                },
                {
                    "metric_id": "auroc",
                    "profile": "genelab_multispecies_pilot",
                    "interpretation": "Ranking quality for binary spaceflight probability outputs.",
                },
                {
                    "metric_id": "calibration_error",
                    "profile": "genelab_multispecies_pilot",
                    "interpretation": "Probability calibration check for small-sample predictions.",
                },
                {
                    "metric_id": "condition_stratum_holdout_delta",
                    "profile": "genelab_multispecies_pilot",
                    "interpretation": (
                        "Secondary robustness diagnostic across genotype/ecotype "
                        "condition-stratum candidate folds."
                    ),
                },
            ],
            "output": {
                "prediction_format": "csv",
                "primary_artifacts": ["predictions.csv", "metrics.json"],
                "required_columns": ["sample_id", "true_label", "predicted_label"],
                "label_domain": ["Ground", "LEO_or_ISS"],
                "positive_label": "LEO_or_ISS",
                "optional_columns": [
                    "leo_probability",
                    "flight_probability",
                    "condition_stratum",
                    "embedding_*",
                ],
            },
            "provenance": {
                "derived_from": "v9/multispecies/source_inventory.draft.json",
                "source_inventory_status": "draft_not_frozen",
                "sample_factor_table": "v9/multispecies/sample_factors.draft.csv",
                "sample_factor_status": "condition_factors_parsed",
                "expression_matrix_audit": "v9/multispecies/expression_matrix_audit.draft.csv",
                "expression_matrix_status": (
                    "matrix_local_sample_aligned"
                    if expression_summary
                    else "pending_expression_matrix_audit"
                ),
                "local_payload_status": "local_sample_table_and_matrix_md5_matched",
                "design_review": (
                    "v9/multispecies/reports/"
                    "MULTISPECIES_CHECKSUM_AND_LOCAL_PAYLOAD_AUDIT.md"
                ),
                **checksum_summary,
            },
            "limitations": [
                "Draft task manifest only; no baseline results are claimed.",
                "This is a species-native task and must not be compared as a raw-gene cross-species task.",
                "The split is within-source and cannot support leave-one-mission-out claims.",
                "Full OSDR payload freeze remains broader than the six local scaffold files hash-checked here.",
                "OSD-120 remains deferred to a separate light/genotype interaction-task design.",
            ],
        }
        validate_task_manifest(manifest)
        manifests.append(manifest)
    if not manifests:
        raise ValueError("no multispecies task manifests were generated")
    return manifests


def _interaction_folds(
    rows: Sequence[Mapping[str, str]],
    factor: str,
    *,
    prefix: str,
    status: str,
    interpretation: str,
) -> list[dict[str, object]]:
    folds = _folds_for_factor(rows, factor, prefix=prefix)
    for fold in folds:
        fold["status"] = status
        fold["interpretation"] = interpretation
    return folds


def build_osd120_interaction_task_manifest(
    rows: Sequence[Mapping[str, str]] | None = None,
    *,
    sample_factor_rows: Sequence[Mapping[str, str]] | None = None,
    expression_matrix_audit_rows: Sequence[Mapping[str, str]] | None = None,
    source_checksum_audit_rows: Sequence[Mapping[str, str]] | None = None,
) -> dict[str, object]:
    """Build the draft OSD-120 Arabidopsis light/genotype interaction manifest."""

    source_rows = list(rows or MULTISPECIES_DRAFT_SOURCES)
    source_by_id = {str(row.get("source_id", "") or ""): row for row in source_rows}
    source_id = "OSD-120"
    source_row = source_by_id.get(source_id)
    if source_row is None:
        raise ValueError("OSD-120 source row is required for the interaction manifest")

    sample_rows = [
        row
        for row in list(sample_factor_rows or [])
        if row.get("source_id") == source_id and row.get("parse_status") == "parsed"
    ]
    if not sample_rows:
        raise ValueError("OSD-120 interaction manifest requires parsed sample factors")

    source_record = _source_record(source_row)
    checksum_summary = _multispecies_checksum_summary(
        list(source_checksum_audit_rows or []),
        source_id,
    )
    if checksum_summary:
        source_record["checksum_status"] = str(
            checksum_summary.get("source_checksum_status", "")
        )
        source_record["local_payload_status"] = "local_sample_table_and_matrix_md5_matched"
    sample_summary = _multispecies_sample_summary(sample_rows)
    expression_summary = _multispecies_expression_summary(
        list(expression_matrix_audit_rows or []),
        source_id,
    )
    genotype_folds = _interaction_folds(
        sample_rows,
        "genotype_or_ecotype",
        prefix="holdout_genotype_or_ecotype",
        status="sample_count_backed_primary_candidate",
        interpretation=(
            "Primary candidate split testing transfer of LEO/Ground discrimination "
            "across Arabidopsis genotype/ecotype context while both light treatments "
            "remain represented in training."
        ),
    )
    light_folds = _interaction_folds(
        sample_rows,
        "light_treatment",
        prefix="holdout_light_treatment",
        status="sample_count_backed_secondary_candidate",
        interpretation=(
            "Secondary candidate split testing transfer across dark versus light "
            "treatment context."
        ),
    )
    condition_folds = _interaction_folds(
        sample_rows,
        "condition_stratum",
        prefix="holdout_condition_stratum",
        status="sample_count_backed_diagnostic_candidate",
        interpretation=(
            "Tertiary diagnostic split holding out one genotype/ecotype by "
            "light-treatment condition stratum at a time."
        ),
    )

    manifest: dict[str, object] = {
        "schema_version": "0.1.0",
        "task_id": MULTISPECIES_OSD120_INTERACTION_TASK_ID,
        "task_family": "multispecies_light_interaction_spaceflight",
        "title": "Arabidopsis root genotype-by-light LEO vs ground draft interaction task",
        "release_status": "draft_not_frozen",
        "organism": str(source_row.get("organism", "") or ""),
        "taxon_id": str(source_row.get("taxon_id", "") or ""),
        "species_common_name": str(source_row.get("species_common_name", "") or ""),
        "material_type": str(source_row.get("material_type", "") or ""),
        "model_system": str(source_row.get("model_system", "") or ""),
        "biospecimen_type": str(source_row.get("biospecimen_type", "") or ""),
        "assay_modality": str(source_row.get("assay_modality", "") or ""),
        "platform": str(source_row.get("platform", "") or ""),
        "spaceflight_environment": str(source_row.get("spaceflight_environment", "") or ""),
        "ground_control_type": str(source_row.get("ground_control_type", "") or ""),
        "donor_or_strain_block": str(source_row.get("donor_or_strain_block", "") or ""),
        "orthology_strategy": "not_applicable_species_native_interaction_task",
        "feature_namespace": str(source_row.get("feature_namespace", "") or ""),
        "variant": str(source_row.get("variants", "") or ""),
        "source_records": [source_record],
        "split": {
            "name": "arabidopsis_root_genotype_light_interaction_pilot",
            "unit": "sample",
            "strategy": (
                "Classify LEO/spaceflight versus ground-control root samples while "
                "explicitly exposing genotype/ecotype and light-treatment interaction "
                "structure through primary, secondary, and diagnostic holdout folds."
            ),
            "status": "draft_sample_matrix_checksum_backed_pending_baseline",
            "target_labels": ["LEO_or_ISS", "Ground"],
            "n_sources": 1,
            "source_accessions": [source_id],
            "sample_factor_table": "v9/multispecies/sample_factors.draft.csv",
            "sample_factor_status": "condition_factors_parsed",
            "blocking_factors": [
                "genotype_or_ecotype",
                "light_treatment",
                "condition_stratum",
            ],
            "primary_candidate_folds": genotype_folds,
            "candidate_folds": genotype_folds,
            "secondary_light_treatment_folds": light_folds,
            "condition_stratum_diagnostic_folds": condition_folds,
            "interaction_design_note": "v9/multispecies/reports/OSD120_INTERACTION_TASK_DESIGN.md",
            "single_source_limitation": (
                "This draft is an OSD-120 within-source interaction task, not a "
                "leave-one-mission-out or cross-species generalization benchmark."
            ),
            **sample_summary,
            **expression_summary,
        },
        "metrics": [
            {
                "metric_id": "balanced_accuracy",
                "profile": "genelab_multispecies_interaction_pilot",
                "interpretation": "Class-balance-aware LEO/spaceflight versus Ground score.",
            },
            {
                "metric_id": "auroc",
                "profile": "genelab_multispecies_interaction_pilot",
                "interpretation": "Ranking quality for binary spaceflight probability outputs.",
            },
            {
                "metric_id": "calibration_error",
                "profile": "genelab_multispecies_interaction_pilot",
                "interpretation": "Probability calibration check for small-sample predictions.",
            },
            {
                "metric_id": "genotype_holdout_delta",
                "profile": "genelab_multispecies_interaction_pilot",
                "interpretation": "Robustness range across genotype/ecotype holdout folds.",
            },
            {
                "metric_id": "light_treatment_holdout_delta",
                "profile": "genelab_multispecies_interaction_pilot",
                "interpretation": "Robustness range across dark/light treatment holdout folds.",
            },
            {
                "metric_id": "condition_stratum_holdout_delta",
                "profile": "genelab_multispecies_interaction_pilot",
                "interpretation": (
                    "Diagnostic robustness range across genotype/ecotype by light "
                    "condition-stratum holdout folds."
                ),
            },
        ],
        "output": {
            "prediction_format": "csv",
            "primary_artifacts": ["predictions.csv", "metrics.json"],
            "required_columns": ["sample_id", "true_label", "predicted_label"],
            "label_domain": ["Ground", "LEO_or_ISS"],
            "positive_label": "LEO_or_ISS",
            "optional_columns": [
                "leo_probability",
                "flight_probability",
                "genotype_or_ecotype",
                "light_treatment",
                "condition_stratum",
                "embedding_*",
            ],
        },
        "provenance": {
            "derived_from": "v9/multispecies/source_inventory.draft.json",
            "source_inventory_status": "draft_not_frozen",
            "sample_factor_table": "v9/multispecies/sample_factors.draft.csv",
            "sample_factor_status": "condition_factors_parsed",
            "expression_matrix_audit": "v9/multispecies/expression_matrix_audit.draft.csv",
            "expression_matrix_status": (
                "matrix_local_sample_aligned"
                if expression_summary
                else "pending_expression_matrix_audit"
            ),
            "local_payload_status": "local_sample_table_and_matrix_md5_matched",
            "design_review": "v9/multispecies/reports/OSD120_INTERACTION_TASK_DESIGN.md",
            **checksum_summary,
        },
        "limitations": [
            "Draft task manifest only; no baseline results are claimed.",
            "This is an Arabidopsis interaction task and should not be merged with the simpler OSD-37 species-native task.",
            "The split is within-source and cannot support leave-one-mission-out claims.",
            "Full OSDR payload freeze remains broader than the six local scaffold files hash-checked here.",
            "Interaction-specific delta metrics are diagnostics, not standalone leaderboard metrics.",
        ],
    }
    validate_task_manifest(manifest)
    return manifest


def _signature_reference_summary(rows: Sequence[Mapping[str, str]]) -> dict[str, object]:
    if not rows:
        return {}
    parsed_rows = [
        row
        for row in rows
        if row.get("audit_status") == "reference_tables_listed_contrast_definitions_parsed"
    ]
    policy_values = {
        str(row.get("recommended_metric_policy", "") or "")
        for row in rows
        if row.get("recommended_metric_policy")
    }
    summary_status = (
        "public_osdr_de_reference_tables_available_pending_contrast_freeze"
        if parsed_rows and len(parsed_rows) == len(rows)
        else _join_unique(rows, "metric_reference_status")
    )
    summary_policy = (
        "keep_classification_primary_enable_de_signature_after_frozen_contrast_subset"
        if policy_values
        == {"keep_classification_primary_enable_de_signature_after_frozen_contrast_subset"}
        else ";".join(sorted(policy_values))
    )
    return {
        "signature_reference_audit": "v9/human_organoid/signature_reference_audit.draft.csv",
        "signature_reference_status": summary_status,
        "signature_metric_policy": summary_policy,
        "signature_reference_sources": {
            str(row.get("source_id", "")): {
                "de_reference_files": str(row.get("de_reference_files", "")),
                "de_reference_file_sizes": str(row.get("de_reference_file_sizes", "")),
                "contrast_files": str(row.get("contrast_files", "")),
                "contrast_pair_count": str(row.get("contrast_pair_count", "")),
                "direct_spaceflight_contrast_count": str(
                    row.get("direct_spaceflight_contrast_count", "")
                ),
                "reversed_spaceflight_contrast_count": str(
                    row.get("reversed_spaceflight_contrast_count", "")
                ),
                "metric_reference_status": str(row.get("metric_reference_status", "")),
            }
            for row in rows
            if row.get("source_id")
        },
        "signature_reference_limitations": (
            "DE/signature metrics are reference-backed at the file-audit level, "
            "and diagnostic scoring is available when a valid response_signature.csv "
            "artifact is supplied. They remain non-primary for leaderboard claims."
        ),
        "de_reference_table": (
            "v9/human_organoid/de_references/human_organoid_de_reference.draft.csv.gz"
        ),
        "de_reference_manifest": (
            "v9/human_organoid/de_references/"
            "human_organoid_de_reference_manifest.draft.json"
        ),
        "response_signature_contract": {
            "artifact": "response_signature.csv",
            "required_columns": RESPONSE_SIGNATURE_REQUIRED_COLUMNS,
            "predicted_log2fc_orientation": BENCHMARK_LOG2FC_ORIENTATION,
            "status": "diagnostic_scorer_available_non_primary",
        },
    }


def build_human_organoid_task_manifest(
    rows: Sequence[Mapping[str, str]] | None = None,
    *,
    sample_factor_rows: Sequence[Mapping[str, str]] | None = None,
    expression_matrix_audit_rows: Sequence[Mapping[str, str]] | None = None,
    signature_reference_audit_rows: Sequence[Mapping[str, str]] | None = None,
) -> dict[str, object]:
    """Build a draft human organoid spaceflight task manifest.

    The manifest deliberately describes a pilot task contract, not a frozen
    leaderboard. Sample-level folds must wait for source payload and sample-table
    audit.
    """

    source_rows = list(rows or HUMAN_ORGANOID_DRAFT_SOURCES)
    if not source_rows:
        raise ValueError("human organoid task manifest requires source rows")
    source_records = [_source_record(row) for row in sorted(source_rows, key=lambda row: row["source_id"])]
    sample_backed = _sample_backed_split(list(sample_factor_rows or []))
    expression_backed = _expression_matrix_summary(list(expression_matrix_audit_rows or []))
    signature_backed = _signature_reference_summary(list(signature_reference_audit_rows or []))
    split_status = (
        "draft_sample_count_backed_pending_baseline"
        if sample_backed
        else "draft_pending_condition_factor_mapping"
    )
    candidate_folds: list[dict[str, object]] = list(
        sample_backed.get("candidate_folds", [])
        if sample_backed
        else [
            {
                "fold_id": "holdout_cortical_organoid_type",
                "test_organoid_type": "cortical_neural_organoid",
                "status": "candidate_only_after_sample_table_audit",
            },
            {
                "fold_id": "holdout_dopaminergic_organoid_type",
                "test_organoid_type": "dopaminergic_neural_organoid",
                "status": "candidate_only_after_sample_table_audit",
            },
            {
                "fold_id": "microglia_condition_blocked_split",
                "unit": "microglia_condition",
                "status": "candidate_only_after_sample_table_audit",
            },
        ]
    )

    manifest: dict[str, object] = {
        "schema_version": "0.1.0",
        "task_id": HUMAN_ORGANOID_TASK_ID,
        "task_family": "human_organoid_spaceflight",
        "title": "Human iPSC neural organoid LEO vs ground draft task",
        "release_status": "draft_not_frozen",
        "organism": "Homo sapiens",
        "taxon_id": "9606",
        "species_common_name": "human",
        "material_type": "cells_cultured",
        "model_system": "human_iPSC_neural_organoid",
        "biospecimen_type": "neural_organoid",
        "assay_modality": "bulk_rna_seq",
        "platform": "Illumina_RNA_seq",
        "spaceflight_environment": "LEO_ISS",
        "ground_control_type": "matched_ground_control",
        "donor_or_strain_block": "donor_block_required",
        "orthology_strategy": "not_applicable_within_human_task",
        "feature_namespace": "human_gene",
        "source_records": source_records,
        "split": {
            "name": "blocked_leo_vs_ground_pilot",
            "unit": "organoid_sample",
            "strategy": (
                "Classify LEO/ISS versus matched ground-control organoid samples "
                "with donor, organoid type, and microglia condition handled as "
                "blocking or stratification factors."
            ),
            "status": split_status,
            "target_labels": ["LEO_or_ISS", "Ground"],
            "n_sources": len(source_records),
            "reported_public_geo_sample_count": 42,
            "sample_table_audit": "v9/human_organoid/sample_table_audit.draft.csv",
            "sample_table_audit_status": "sample_table_parsed",
            "source_accessions": [record["source_id"] for record in source_records],
            "external_accessions": _join_unique(source_records, "external_accessions"),
            "organoid_types": _join_unique(source_records, "organoid_type"),
            "microglia_conditions": "with_iPSC_derived_microglia;without_iPSC_derived_microglia",
            "blocking_factors": [
                "donor_subject",
                "organoid_type",
                "microglia_condition",
            ],
            "candidate_folds": candidate_folds,
            "single_mission_limitation": (
                "OSD-863/OSD-871 derive from one ISS mission context, so this "
                "draft is not a leave-one-mission-out benchmark."
            ),
            **{key: value for key, value in sample_backed.items() if key != "candidate_folds"},
            **expression_backed,
        },
        "metrics": [
            {
                "metric_id": "balanced_accuracy",
                "profile": "genelab_organoid_pilot",
                "interpretation": "Class-balance-aware LEO/ISS versus Ground score.",
            },
            {
                "metric_id": "auroc",
                "profile": "genelab_organoid_pilot",
                "interpretation": "Ranking quality for binary LEO/ISS probability outputs.",
            },
            {
                "metric_id": "calibration_error",
                "profile": "genelab_organoid_pilot",
                "interpretation": "Probability calibration check for small-sample predictions.",
            },
            {
                "metric_id": "de_direction_match",
                "profile": "genelab_organoid_pilot",
                "interpretation": (
                    "Diagnostic direction agreement for response genes when a "
                    "valid response_signature artifact is supplied."
                ),
            },
            {
                "metric_id": "signature_rank_correlation",
                "profile": "genelab_organoid_pilot",
                "interpretation": (
                    "Diagnostic rank correlation between predicted and observed "
                    "organoid response signatures."
                ),
            },
            {
                "metric_id": "mission_discrimination",
                "profile": "genelab_organoid_pilot",
                "interpretation": "Representation diagnostic; not a mission-held-out claim for this single-mission pilot.",
            },
        ],
        "output": {
            "prediction_format": "csv",
            "primary_artifacts": ["predictions.csv", "metrics.json"],
            "required_columns": ["sample_id", "true_label", "predicted_label"],
            "label_domain": ["Ground", "LEO_or_ISS"],
            "positive_label": "LEO_or_ISS",
            "optional_columns": [
                "leo_probability",
                "flight_probability",
                "response_signature_score",
                "embedding_*",
            ],
        },
        "provenance": {
            "derived_from": "v9/human_organoid/source_inventory.draft.json",
            "source_inventory_status": "draft_not_frozen",
            "source_checksum_audit": "v9/human_organoid/source_checksum_audit.draft.csv",
            "source_checksum_status": "checksum_manifests_parsed_payloads_not_hashed",
            "source_payload_status": "not_downloaded_or_payload_hashed",
            "sample_table_audit": "v9/human_organoid/sample_table_audit.draft.csv",
            "sample_table_status": "sample_table_parsed",
            "sample_factor_status": (
                str(sample_backed.get("sample_factor_status"))
                if sample_backed
                else "pending_condition_factor_mapping"
            ),
            "geo_sample_metadata": "v9/human_organoid/geo_sample_metadata.draft.csv",
            "donor_or_line_status": (
                str(sample_backed.get("donor_metadata_status"))
                if sample_backed
                else "pending_geo_metadata_review"
            ),
            "expression_matrix_audit": "v9/human_organoid/expression_matrix_audit.draft.csv",
            "expression_matrix_status": (
                "matrix_downloaded_sample_aligned"
                if expression_backed
                else "pending_expression_matrix_audit"
            ),
            "signature_reference_audit": "v9/human_organoid/signature_reference_audit.draft.csv",
            "signature_reference_status": (
                str(signature_backed.get("signature_reference_status"))
                if signature_backed
                else "pending_signature_reference_audit"
            ),
            "signature_metric_policy": (
                str(signature_backed.get("signature_metric_policy"))
                if signature_backed
                else "classifier_primary_until_signature_reference_audit"
            ),
            "publication_urls": _join_unique(source_records, "publication_urls"),
        },
        "reference_signatures": signature_backed
        or {
            "signature_reference_audit": "v9/human_organoid/signature_reference_audit.draft.csv",
            "signature_reference_status": "pending_signature_reference_audit",
            "signature_metric_policy": "classifier_primary_until_signature_reference_audit",
        },
        "limitations": [
            "Draft task manifest only; no baseline results are claimed.",
            "Payload checksum audit is pending for OSD-863/OSD-871.",
            "Any organoid baseline output is draft-only and not a frozen leaderboard claim.",
            "single-mission organoid data cannot support leave-one-mission-out claims.",
            (
                "GEO-derived Subject/iPSC-line metadata are available, but donor, source, "
                "organoid fate, and disease context are not independently crossed; donor "
                "generalization remains pilot-only."
            ),
            (
                "DE/signature metrics are diagnostic-only and non-primary even when "
                "a response_signature artifact is supplied."
            ),
        ],
    }
    validate_task_manifest(manifest)
    return manifest


def write_task_manifest(manifest: Mapping[str, object], path: str | Path) -> Path:
    """Write a draft task manifest with stable JSON formatting."""

    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return output_path
