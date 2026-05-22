"""Draft task manifests for v9 extension tracks."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Mapping, Sequence

from .extension_sources import HUMAN_ORGANOID_DRAFT_SOURCES
from .manifests import validate_task_manifest


HUMAN_ORGANOID_TASK_ID = "draft_human_organoid_spaceflight"


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
            "but leaderboard scoring still requires a frozen contrast subset and "
            "explicit log2 fold-change orientation."
        ),
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
                    "Direction agreement for response genes after public OSDR DE "
                    "tables are pinned to a frozen contrast subset."
                ),
            },
            {
                "metric_id": "signature_rank_correlation",
                "profile": "genelab_organoid_pilot",
                "interpretation": (
                    "Rank correlation between predicted and observed organoid response "
                    "signatures after log2 fold-change orientation is frozen."
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
                "DE/signature metrics remain non-primary until public OSDR DE references "
                "are pinned to a frozen contrast subset."
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
