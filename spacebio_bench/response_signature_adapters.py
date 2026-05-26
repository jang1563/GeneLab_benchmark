"""Response-signature adapters for draft human organoid tasks."""

from __future__ import annotations

import csv
import gzip
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from .data import HumanOrganoidTaskData, OrganoidFoldData


SOURCE_TRANSFER_SIGNATURE_ID = "organoid_source_transfer_empirical_signature"
MICROGLIA_MATCHED_SOURCE_TRANSFER_SIGNATURE_ID = (
    "organoid_microglia_matched_source_transfer_empirical_signature"
)
SHARED_CONTROL_SOURCE_TRANSFER_SIGNATURE_ID = (
    "organoid_shared_control_disease_microglia_source_transfer_empirical_signature"
)
REFERENCE_USAGE_POLICY = "reference_not_used_for_signature_generation"
SHARED_CONTROL_DISEASE_CONTEXT = "no_known_diseases"
SOURCE_TRANSFER_RESPONSE_COLUMNS = [
    "task_id",
    "source_id",
    "contrast_id",
    "gene_symbol",
    "ensembl_id",
    "predicted_log2fc_leo_or_iss_minus_ground",
    "fold_id",
    "signature_model_id",
    "training_source_id",
    "target_source_id",
    "training_scope",
    "conditioning_strategy",
    "conditioning_factor",
    "conditioning_value",
    "target_disease_context",
    "target_microglia_condition",
    "n_train_ground",
    "n_train_leo_or_iss",
    "n_condition_train_ground",
    "n_condition_train_leo_or_iss",
    "signature_generation_method",
    "reference_usage_policy",
]


@dataclass(frozen=True)
class SourceTransferSignatureResult:
    """Generated response-signature rows and provenance for one adapter run."""

    signature_model_id: str
    response_rows: list[dict[str, str]]
    fold_summaries: list[dict[str, Any]]
    target_contrasts_by_source: dict[str, list[str]]
    contrast_summaries: list[dict[str, Any]] | None = None

    @property
    def n_response_rows(self) -> int:
        return len(self.response_rows)


def _format_float(value: float) -> str:
    return f"{value:.12g}"


def _sample_factor_by_id(task: HumanOrganoidTaskData) -> dict[str, dict[str, str]]:
    return {
        str(row["sample_id"]): row
        for row in task.sample_factors
        if row.get("parse_status") == "parsed" and row.get("sample_id")
    }


def _canonical_factor_value(value: object) -> str:
    text = str(value or "").strip().lower()
    chars = [char if char.isalnum() else "_" for char in text]
    return "_".join(part for part in "".join(chars).split("_") if part)


def _source_ids_for_samples(
    sample_factors_by_id: Mapping[str, Mapping[str, str]],
    sample_ids: Sequence[str],
) -> set[str]:
    return {
        str(sample_factors_by_id[sample_id].get("source_id", "") or "")
        for sample_id in sample_ids
        if sample_id in sample_factors_by_id
    } - {""}


def _source_transfer_folds(task: HumanOrganoidTaskData) -> list[tuple[str, str, OrganoidFoldData]]:
    sample_factors_by_id = _sample_factor_by_id(task)
    source_folds: list[tuple[str, str, OrganoidFoldData]] = []
    for fold in task.folds:
        train_sources = _source_ids_for_samples(sample_factors_by_id, fold.train_sample_ids)
        target_sources = _source_ids_for_samples(sample_factors_by_id, fold.test_sample_ids)
        if len(train_sources) != 1 or len(target_sources) != 1:
            continue
        training_source = next(iter(train_sources))
        target_source = next(iter(target_sources))
        if training_source == target_source:
            continue
        source_folds.append((target_source, training_source, fold))
    if not source_folds:
        raise ValueError("no single-source-to-single-source transfer folds found")
    return sorted(source_folds, key=lambda item: item[0])


def source_contrasts_from_de_manifest(path: str | Path) -> dict[str, list[str]]:
    """Read source-specific contrast ids from the derived DE-reference manifest."""

    metadata = source_contrast_metadata_from_de_manifest(path)
    return {
        source_id: [str(contrast["contrast_id"]) for contrast in contrasts]
        for source_id, contrasts in metadata.items()
    }


def source_contrast_metadata_from_de_manifest(path: str | Path) -> dict[str, list[dict[str, str]]]:
    """Read source-specific contrast metadata from the derived DE-reference manifest."""

    manifest = json.loads(Path(path).read_text())
    contrasts_by_source: dict[str, list[dict[str, str]]] = {}
    for source in manifest.get("sources", []):
        if not isinstance(source, Mapping):
            continue
        source_id = str(source.get("source_id", "") or "")
        contrasts = [
            {
                "source_id": source_id,
                "contrast_id": str(contrast.get("contrast_id", "") or ""),
                "disease_context": _canonical_factor_value(contrast.get("disease_context", "")),
                "microglia_condition": _canonical_factor_value(
                    contrast.get("microglia_condition", "")
                ),
            }
            for contrast in source.get("contrasts", [])
            if isinstance(contrast, Mapping) and contrast.get("contrast_id")
        ]
        if source_id and contrasts:
            contrasts_by_source[source_id] = contrasts
    if not contrasts_by_source:
        raise ValueError(f"DE-reference manifest has no source contrasts: {path}")
    return contrasts_by_source


def _label_counts_for_samples(
    sample_factors_by_id: Mapping[str, Mapping[str, str]],
    sample_ids: Sequence[str],
) -> dict[str, int]:
    ground = 0
    leo = 0
    for sample_id in sample_ids:
        label = str(sample_factors_by_id[sample_id].get("true_label", "") or "")
        if label == "Ground":
            ground += 1
        elif label == "LEO_or_ISS":
            leo += 1
    return {
        "n_train_ground": ground,
        "n_train_leo_or_iss": leo,
        "n_train_total": ground + leo,
    }


def _training_signature_for_sample_ids(
    task: HumanOrganoidTaskData,
    sample_ids: Sequence[str],
) -> tuple[np.ndarray, dict[str, int]]:
    sample_factors_by_id = _sample_factor_by_id(task)
    train_features = task.features.loc[list(sample_ids)]
    labels = [sample_factors_by_id[sample_id]["true_label"] for sample_id in train_features.index]
    ground_mask = np.asarray([label == "Ground" for label in labels], dtype=bool)
    leo_mask = np.asarray([label == "LEO_or_ISS" for label in labels], dtype=bool)
    if not ground_mask.any() or not leo_mask.any():
        raise ValueError("training signature requires both Ground and LEO_or_ISS samples")
    values = train_features.to_numpy(dtype=np.float64)
    ground_mean = values[ground_mask].mean(axis=0)
    leo_mean = values[leo_mask].mean(axis=0)
    signature = np.log2(np.clip(leo_mean, a_min=0.0, a_max=None) + 1.0) - np.log2(
        np.clip(ground_mean, a_min=0.0, a_max=None) + 1.0
    )
    return signature, {
        "n_train_ground": int(ground_mask.sum()),
        "n_train_leo_or_iss": int(leo_mask.sum()),
        "n_train_total": int(len(labels)),
    }


def _training_signature(
    task: HumanOrganoidTaskData,
    fold: OrganoidFoldData,
) -> tuple[np.ndarray, dict[str, int]]:
    return _training_signature_for_sample_ids(task, fold.train_sample_ids)


def build_source_transfer_response_signature(
    task: HumanOrganoidTaskData,
    *,
    de_reference_manifest_path: str | Path,
    max_features: int | None = None,
) -> SourceTransferSignatureResult:
    """Generate a non-leaky source-transfer empirical response signature."""

    contrasts_by_source = source_contrasts_from_de_manifest(de_reference_manifest_path)
    response_rows: list[dict[str, str]] = []
    fold_summaries: list[dict[str, Any]] = []
    features = [str(feature) for feature in task.features.columns]
    if max_features is not None:
        if max_features < 1:
            raise ValueError("max_features must be at least 1 when provided")
        features = features[:max_features]
    for target_source, training_source, fold in _source_transfer_folds(task):
        contrast_ids = contrasts_by_source.get(target_source, [])
        if not contrast_ids:
            continue
        signature, counts = _training_signature(task, fold)
        if max_features is not None:
            signature = signature[:max_features]
        fold_summaries.append(
            {
                "fold_id": fold.fold_id,
                "target_source_id": target_source,
                "training_source_id": training_source,
                "heldout_factor": fold.heldout_factor,
                "heldout_value": fold.heldout_value,
                **counts,
                "n_features": len(features),
                "n_target_contrasts": len(contrast_ids),
                "reference_usage_policy": REFERENCE_USAGE_POLICY,
            }
        )
        for contrast_id in contrast_ids:
            for feature, value in zip(features, signature):
                response_rows.append(
                    {
                        "task_id": task.task_id,
                        "source_id": target_source,
                        "contrast_id": contrast_id,
                        "gene_symbol": "",
                        "ensembl_id": feature,
                        "predicted_log2fc_leo_or_iss_minus_ground": _format_float(float(value)),
                        "fold_id": fold.fold_id,
                        "signature_model_id": SOURCE_TRANSFER_SIGNATURE_ID,
                        "training_source_id": training_source,
                        "target_source_id": target_source,
                        "training_scope": "source_transfer_organoid_type_holdout_train_samples",
                        "conditioning_strategy": "global_source_transfer",
                        "conditioning_factor": "",
                        "conditioning_value": "",
                        "target_disease_context": "",
                        "target_microglia_condition": "",
                        "n_train_ground": str(counts["n_train_ground"]),
                        "n_train_leo_or_iss": str(counts["n_train_leo_or_iss"]),
                        "n_condition_train_ground": "",
                        "n_condition_train_leo_or_iss": "",
                        "signature_generation_method": "log2_mean_leo_or_iss_plus1_minus_log2_mean_ground_plus1",
                        "reference_usage_policy": REFERENCE_USAGE_POLICY,
                    }
                )

    if not response_rows:
        raise ValueError("source-transfer adapter produced no response-signature rows")
    return SourceTransferSignatureResult(
        signature_model_id=SOURCE_TRANSFER_SIGNATURE_ID,
        response_rows=response_rows,
        fold_summaries=fold_summaries,
        target_contrasts_by_source={
            source_id: list(contrast_ids)
            for source_id, contrast_ids in sorted(contrasts_by_source.items())
        },
    )


def build_microglia_matched_source_transfer_response_signature(
    task: HumanOrganoidTaskData,
    *,
    de_reference_manifest_path: str | Path,
    max_features: int | None = None,
) -> SourceTransferSignatureResult:
    """Generate a source-transfer empirical response signature matched by microglia status."""

    contrasts_by_source = source_contrast_metadata_from_de_manifest(de_reference_manifest_path)
    sample_factors_by_id = _sample_factor_by_id(task)
    response_rows: list[dict[str, str]] = []
    fold_summaries: list[dict[str, Any]] = []
    contrast_summaries: list[dict[str, Any]] = []
    features = [str(feature) for feature in task.features.columns]
    if max_features is not None:
        if max_features < 1:
            raise ValueError("max_features must be at least 1 when provided")
        features = features[:max_features]

    for target_source, training_source, fold in _source_transfer_folds(task):
        contrast_metadata = contrasts_by_source.get(target_source, [])
        if not contrast_metadata:
            continue
        fold_counts = _label_counts_for_samples(sample_factors_by_id, fold.train_sample_ids)
        fold_summaries.append(
            {
                "fold_id": fold.fold_id,
                "target_source_id": target_source,
                "training_source_id": training_source,
                "heldout_factor": fold.heldout_factor,
                "heldout_value": fold.heldout_value,
                **fold_counts,
                "n_features": len(features),
                "n_target_contrasts": len(contrast_metadata),
                "conditioning_strategy": "microglia_matched_source_transfer",
                "conditioning_factor": "microglia_condition",
                "reference_usage_policy": REFERENCE_USAGE_POLICY,
            }
        )
        for contrast in contrast_metadata:
            microglia_condition = str(contrast["microglia_condition"])
            disease_context = str(contrast["disease_context"])
            condition_sample_ids = [
                sample_id
                for sample_id in fold.train_sample_ids
                if sample_id in sample_factors_by_id
                and sample_factors_by_id[sample_id].get("source_id") == training_source
                and sample_factors_by_id[sample_id].get("microglia_condition")
                == microglia_condition
            ]
            condition_counts = _label_counts_for_samples(sample_factors_by_id, condition_sample_ids)
            summary = {
                "fold_id": fold.fold_id,
                "contrast_id": str(contrast["contrast_id"]),
                "target_source_id": target_source,
                "training_source_id": training_source,
                "target_disease_context": disease_context,
                "target_microglia_condition": microglia_condition,
                "conditioning_strategy": "microglia_matched_source_transfer",
                "conditioning_factor": "microglia_condition",
                "conditioning_value": microglia_condition,
                "n_condition_train_ground": condition_counts["n_train_ground"],
                "n_condition_train_leo_or_iss": condition_counts["n_train_leo_or_iss"],
                "n_features": len(features),
                "reference_usage_policy": REFERENCE_USAGE_POLICY,
            }
            if not condition_counts["n_train_ground"] or not condition_counts["n_train_leo_or_iss"]:
                summary["status"] = "skipped_missing_condition_label_pair"
                contrast_summaries.append(summary)
                continue
            signature, _ = _training_signature_for_sample_ids(task, condition_sample_ids)
            if max_features is not None:
                signature = signature[:max_features]
            summary["status"] = "emitted"
            contrast_summaries.append(summary)
            for feature, value in zip(features, signature):
                response_rows.append(
                    {
                        "task_id": task.task_id,
                        "source_id": target_source,
                        "contrast_id": str(contrast["contrast_id"]),
                        "gene_symbol": "",
                        "ensembl_id": feature,
                        "predicted_log2fc_leo_or_iss_minus_ground": _format_float(float(value)),
                        "fold_id": fold.fold_id,
                        "signature_model_id": MICROGLIA_MATCHED_SOURCE_TRANSFER_SIGNATURE_ID,
                        "training_source_id": training_source,
                        "target_source_id": target_source,
                        "training_scope": (
                            "source_transfer_organoid_type_holdout_train_samples_"
                            "microglia_matched"
                        ),
                        "conditioning_strategy": "microglia_matched_source_transfer",
                        "conditioning_factor": "microglia_condition",
                        "conditioning_value": microglia_condition,
                        "target_disease_context": disease_context,
                        "target_microglia_condition": microglia_condition,
                        "n_train_ground": str(fold_counts["n_train_ground"]),
                        "n_train_leo_or_iss": str(fold_counts["n_train_leo_or_iss"]),
                        "n_condition_train_ground": str(condition_counts["n_train_ground"]),
                        "n_condition_train_leo_or_iss": str(
                            condition_counts["n_train_leo_or_iss"]
                        ),
                        "signature_generation_method": (
                            "log2_mean_microglia_matched_leo_or_iss_plus1_minus_"
                            "log2_mean_microglia_matched_ground_plus1"
                        ),
                        "reference_usage_policy": REFERENCE_USAGE_POLICY,
                    }
                )

    if not response_rows:
        raise ValueError("microglia-matched source-transfer adapter produced no rows")
    return SourceTransferSignatureResult(
        signature_model_id=MICROGLIA_MATCHED_SOURCE_TRANSFER_SIGNATURE_ID,
        response_rows=response_rows,
        fold_summaries=fold_summaries,
        target_contrasts_by_source={
            source_id: [str(contrast["contrast_id"]) for contrast in contrasts]
            for source_id, contrasts in sorted(contrasts_by_source.items())
        },
        contrast_summaries=contrast_summaries,
    )


def build_shared_control_source_transfer_response_signature(
    task: HumanOrganoidTaskData,
    *,
    de_reference_manifest_path: str | Path,
    max_features: int | None = None,
) -> SourceTransferSignatureResult:
    """Generate a partial shared-control disease+microglia source-transfer signature."""

    contrasts_by_source = source_contrast_metadata_from_de_manifest(de_reference_manifest_path)
    sample_factors_by_id = _sample_factor_by_id(task)
    response_rows: list[dict[str, str]] = []
    fold_summaries: list[dict[str, Any]] = []
    contrast_summaries: list[dict[str, Any]] = []
    features = [str(feature) for feature in task.features.columns]
    if max_features is not None:
        if max_features < 1:
            raise ValueError("max_features must be at least 1 when provided")
        features = features[:max_features]

    for target_source, training_source, fold in _source_transfer_folds(task):
        contrast_metadata = contrasts_by_source.get(target_source, [])
        if not contrast_metadata:
            continue
        fold_counts = _label_counts_for_samples(sample_factors_by_id, fold.train_sample_ids)
        fold_emitted = 0
        fold_skipped = 0
        for contrast in contrast_metadata:
            disease_context = str(contrast["disease_context"])
            microglia_condition = str(contrast["microglia_condition"])
            condition_sample_ids = [
                sample_id
                for sample_id in fold.train_sample_ids
                if sample_id in sample_factors_by_id
                and sample_factors_by_id[sample_id].get("source_id") == training_source
                and sample_factors_by_id[sample_id].get("disease_context") == disease_context
                and sample_factors_by_id[sample_id].get("microglia_condition")
                == microglia_condition
            ]
            condition_counts = _label_counts_for_samples(sample_factors_by_id, condition_sample_ids)
            summary = {
                "fold_id": fold.fold_id,
                "contrast_id": str(contrast["contrast_id"]),
                "target_source_id": target_source,
                "training_source_id": training_source,
                "target_disease_context": disease_context,
                "target_microglia_condition": microglia_condition,
                "conditioning_strategy": "shared_control_disease_microglia_source_transfer",
                "conditioning_factor": "disease_context+microglia_condition",
                "conditioning_value": f"{disease_context}|{microglia_condition}",
                "n_condition_train_ground": condition_counts["n_train_ground"],
                "n_condition_train_leo_or_iss": condition_counts["n_train_leo_or_iss"],
                "n_features": len(features),
                "reference_usage_policy": REFERENCE_USAGE_POLICY,
            }
            if disease_context != SHARED_CONTROL_DISEASE_CONTEXT:
                summary["status"] = "skipped_missing_shared_disease_context"
                summary["skip_reason"] = (
                    f"{disease_context} is not shared across OSD-863 and OSD-871"
                )
                contrast_summaries.append(summary)
                fold_skipped += 1
                continue
            if not condition_counts["n_train_ground"] or not condition_counts["n_train_leo_or_iss"]:
                summary["status"] = "skipped_missing_condition_label_pair"
                summary["skip_reason"] = "matched shared-control stratum lacks Ground or LEO/ISS"
                contrast_summaries.append(summary)
                fold_skipped += 1
                continue
            signature, _ = _training_signature_for_sample_ids(task, condition_sample_ids)
            if max_features is not None:
                signature = signature[:max_features]
            summary["status"] = "emitted"
            contrast_summaries.append(summary)
            fold_emitted += 1
            for feature, value in zip(features, signature):
                response_rows.append(
                    {
                        "task_id": task.task_id,
                        "source_id": target_source,
                        "contrast_id": str(contrast["contrast_id"]),
                        "gene_symbol": "",
                        "ensembl_id": feature,
                        "predicted_log2fc_leo_or_iss_minus_ground": _format_float(float(value)),
                        "fold_id": fold.fold_id,
                        "signature_model_id": SHARED_CONTROL_SOURCE_TRANSFER_SIGNATURE_ID,
                        "training_source_id": training_source,
                        "target_source_id": target_source,
                        "training_scope": (
                            "source_transfer_organoid_type_holdout_train_samples_"
                            "shared_control_disease_microglia_matched"
                        ),
                        "conditioning_strategy": (
                            "shared_control_disease_microglia_source_transfer"
                        ),
                        "conditioning_factor": "disease_context+microglia_condition",
                        "conditioning_value": f"{disease_context}|{microglia_condition}",
                        "target_disease_context": disease_context,
                        "target_microglia_condition": microglia_condition,
                        "n_train_ground": str(fold_counts["n_train_ground"]),
                        "n_train_leo_or_iss": str(fold_counts["n_train_leo_or_iss"]),
                        "n_condition_train_ground": str(condition_counts["n_train_ground"]),
                        "n_condition_train_leo_or_iss": str(
                            condition_counts["n_train_leo_or_iss"]
                        ),
                        "signature_generation_method": (
                            "log2_mean_shared_control_disease_microglia_matched_"
                            "leo_or_iss_plus1_minus_log2_mean_matched_ground_plus1"
                        ),
                        "reference_usage_policy": REFERENCE_USAGE_POLICY,
                    }
                )
        fold_summaries.append(
            {
                "fold_id": fold.fold_id,
                "target_source_id": target_source,
                "training_source_id": training_source,
                "heldout_factor": fold.heldout_factor,
                "heldout_value": fold.heldout_value,
                **fold_counts,
                "n_features": len(features),
                "n_target_contrasts": len(contrast_metadata),
                "n_emitted_contrasts": fold_emitted,
                "n_skipped_contrasts": fold_skipped,
                "conditioning_strategy": "shared_control_disease_microglia_source_transfer",
                "conditioning_factor": "disease_context+microglia_condition",
                "partial_coverage": True,
                "reference_usage_policy": REFERENCE_USAGE_POLICY,
            }
        )

    if not response_rows:
        raise ValueError("shared-control source-transfer adapter produced no rows")
    return SourceTransferSignatureResult(
        signature_model_id=SHARED_CONTROL_SOURCE_TRANSFER_SIGNATURE_ID,
        response_rows=response_rows,
        fold_summaries=fold_summaries,
        target_contrasts_by_source={
            source_id: [str(contrast["contrast_id"]) for contrast in contrasts]
            for source_id, contrasts in sorted(contrasts_by_source.items())
        },
        contrast_summaries=contrast_summaries,
    )


def write_source_transfer_response_signature(
    result: SourceTransferSignatureResult,
    *,
    response_signature_path: str | Path,
    metadata_path: str | Path,
) -> tuple[Path, Path]:
    """Write response-signature rows and adapter metadata."""

    response_path = Path(response_signature_path)
    meta_path = Path(metadata_path)
    response_path.parent.mkdir(parents=True, exist_ok=True)
    meta_path.parent.mkdir(parents=True, exist_ok=True)
    if response_path.name.endswith(".gz"):
        handle_context = gzip.open(response_path, "wt", newline="")
    else:
        handle_context = response_path.open("w", newline="")
    with handle_context as handle:
        writer = csv.DictWriter(handle, fieldnames=SOURCE_TRANSFER_RESPONSE_COLUMNS)
        writer.writeheader()
        for row in result.response_rows:
            writer.writerow({field: row.get(field, "") for field in SOURCE_TRANSFER_RESPONSE_COLUMNS})
    metadata = {
        "schema_version": "0.1.0",
        "artifact_type": "human_organoid_source_transfer_response_signature",
        "signature_model_id": result.signature_model_id,
        "release_status": "draft_not_frozen",
        "claim_boundary": "diagnostic_only_not_leaderboard",
        "reference_usage_policy": REFERENCE_USAGE_POLICY,
        "n_response_rows": result.n_response_rows,
        "fold_summaries": result.fold_summaries,
        "target_contrasts_by_source": result.target_contrasts_by_source,
        "notes": (
            "Response signatures are generated from source-transfer training samples only. "
            "The DE reference is not used during signature generation and should be used "
            "only for post hoc diagnostic scoring."
        ),
    }
    if result.contrast_summaries is not None:
        metadata["contrast_summaries"] = result.contrast_summaries
    meta_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    return response_path, meta_path


def write_microglia_matched_source_transfer_response_signature(
    result: SourceTransferSignatureResult,
    *,
    response_signature_path: str | Path,
    metadata_path: str | Path,
) -> tuple[Path, Path]:
    """Write microglia-matched response-signature rows and adapter metadata."""

    response_path, meta_path = write_source_transfer_response_signature(
        result,
        response_signature_path=response_signature_path,
        metadata_path=metadata_path,
    )
    metadata = json.loads(meta_path.read_text())
    metadata.update(
        {
            "artifact_type": "human_organoid_microglia_matched_source_transfer_response_signature",
            "conditioning_strategy": "microglia_matched_source_transfer",
            "conditioning_factor": "microglia_condition",
            "notes": (
                "Response signatures are generated from source-transfer training samples "
                "matched to the target contrast microglia condition. The DE reference is "
                "not used during signature generation and should be used only for post "
                "hoc diagnostic scoring."
            ),
        }
    )
    meta_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    return response_path, meta_path


def write_shared_control_source_transfer_response_signature(
    result: SourceTransferSignatureResult,
    *,
    response_signature_path: str | Path,
    metadata_path: str | Path,
) -> tuple[Path, Path]:
    """Write shared-control response-signature rows and adapter metadata."""

    response_path, meta_path = write_source_transfer_response_signature(
        result,
        response_signature_path=response_signature_path,
        metadata_path=metadata_path,
    )
    metadata = json.loads(meta_path.read_text())
    contrast_summaries = result.contrast_summaries or []
    n_emitted = sum(1 for summary in contrast_summaries if summary.get("status") == "emitted")
    n_skipped = sum(1 for summary in contrast_summaries if str(summary.get("status", "")).startswith("skipped"))
    metadata.update(
        {
            "artifact_type": "human_organoid_shared_control_source_transfer_response_signature",
            "conditioning_strategy": "shared_control_disease_microglia_source_transfer",
            "conditioning_factor": "disease_context+microglia_condition",
            "partial_coverage": True,
            "partial_coverage_policy": (
                "emit shared no_known_diseases target contrasts only; skip "
                "source-specific disease contexts without fallback"
            ),
            "n_emitted_contrasts": n_emitted,
            "n_skipped_contrasts": n_skipped,
            "notes": (
                "Response signatures are generated from source-transfer training samples "
                "matched to shared no_known_diseases plus target microglia condition. "
                "Source-specific disease contexts are skipped without fallback. The DE "
                "reference is not used during signature generation and should be used "
                "only for post hoc diagnostic scoring."
            ),
        }
    )
    meta_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    return response_path, meta_path
