"""Feature-effect artifacts and diagnostics for draft human organoid tasks."""

from __future__ import annotations

import csv
import gzip
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from .data import HumanOrganoidTaskData, OrganoidFoldData
from .response_signature_adapters import source_contrast_metadata_from_de_manifest


FEATURE_EFFECT_REQUIRED_COLUMNS = [
    "task_id",
    "fold_id",
    "source_id",
    "contrast_id",
    "feature_namespace",
    "gene_symbol",
    "ensembl_id",
    "feature_id",
    "effect_value",
    "effect_direction_positive_class",
    "effect_scale",
    "model_space",
    "classifier_model_id",
    "training_scope",
    "target_scope",
    "positive_label",
    "reference_usage_policy",
]
FEATURE_EFFECT_POLICY = "reference_not_used_for_effect_generation"
LOGISTIC_FEATURE_EFFECT_ID = "organoid_l2_logistic_gene_space_feature_effect"
PCA_LR_FEATURE_EFFECT_ID = "organoid_pca_lr_reconstructed_gene_space_feature_effect"
NEGATIVE_LABEL = "Ground"
POSITIVE_LABEL = "LEO_or_ISS"
EFFECT_VALUE_COLUMN = "effect_value"
REFERENCE_LOG2FC_COLUMN = "log2fc_leo_or_iss_minus_ground"


@dataclass(frozen=True)
class FeatureEffectResult:
    """Generated feature-effect rows and provenance for one run."""

    classifier_model_id: str
    effect_rows: list[dict[str, str]]
    fold_summaries: list[dict[str, Any]]
    target_contrasts_by_source: dict[str, list[str]]

    @property
    def n_effect_rows(self) -> int:
        return len(self.effect_rows)


def _format_float(value: float) -> str:
    return f"{value:.12g}"


def _open_text(path: str | Path):
    parsed = Path(path)
    if parsed.name.endswith(".gz"):
        return gzip.open(parsed, "rt", newline="")
    return parsed.open(newline="")


def _read_rows(path: str | Path) -> list[dict[str, str]]:
    with _open_text(path) as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def _as_float(value: str) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(parsed):
        return None
    return parsed


def _sign(value: float) -> int:
    if value > 0:
        return 1
    if value < 0:
        return -1
    return 0


def _rankdata(values: Sequence[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda index: values[index])
    ranks = [0.0] * len(values)
    cursor = 0
    while cursor < len(order):
        end = cursor + 1
        while end < len(order) and values[order[end]] == values[order[cursor]]:
            end += 1
        rank = (cursor + 1 + end) / 2.0
        for index in order[cursor:end]:
            ranks[index] = rank
        cursor = end
    return ranks


def _spearman(left: Sequence[float], right: Sequence[float]) -> float | None:
    if len(left) < 2 or len(right) < 2 or len(left) != len(right):
        return None
    left_ranks = np.asarray(_rankdata(left), dtype=float)
    right_ranks = np.asarray(_rankdata(right), dtype=float)
    if float(left_ranks.std()) == 0.0 or float(right_ranks.std()) == 0.0:
        return None
    return float(np.corrcoef(left_ranks, right_ranks)[0, 1])


def validate_feature_effect_artifact(
    path: str | Path,
    *,
    task_id: str | None = None,
) -> dict[str, Any]:
    """Validate the required feature-effect artifact surface."""

    parsed = Path(path)
    report: dict[str, Any] = {
        "ok": False,
        "path": parsed.as_posix(),
        "errors": [],
        "warnings": [],
        "columns": [],
        "n_rows": 0,
    }
    if not parsed.exists():
        report["errors"].append(f"feature effect artifact missing: {parsed}")
        return report
    try:
        with _open_text(parsed) as handle:
            reader = csv.DictReader(handle)
            columns = list(reader.fieldnames or [])
            report["columns"] = columns
            missing = [column for column in FEATURE_EFFECT_REQUIRED_COLUMNS if column not in columns]
            if missing:
                report["errors"].append("missing required columns: " + ", ".join(missing))
                return report
            for row_index, row in enumerate(reader, start=1):
                report["n_rows"] += 1
                if task_id is not None and row.get("task_id") != task_id:
                    report["errors"].append(
                        f"row {row_index} has task_id={row.get('task_id')!r}, expected {task_id!r}"
                    )
                    break
                if not row.get("gene_symbol") and not row.get("ensembl_id"):
                    report["errors"].append(f"row {row_index} is missing both gene_symbol and ensembl_id")
                    break
                if _as_float(str(row.get(EFFECT_VALUE_COLUMN, ""))) is None:
                    report["errors"].append(f"row {row_index} has non-numeric effect_value")
                    break
                if row.get("reference_usage_policy") != FEATURE_EFFECT_POLICY:
                    report["errors"].append(
                        f"row {row_index} has invalid reference_usage_policy"
                    )
                    break
    except OSError as exc:
        report["errors"].append(str(exc))
        return report

    if report["n_rows"] == 0:
        report["errors"].append("feature effect artifact has no rows")
    report["ok"] = not report["errors"]
    return report


def _reference_rows_by_key(rows: Iterable[Mapping[str, str]]) -> tuple[dict[tuple[str, str, str], Mapping[str, str]], dict[tuple[str, str, str], Mapping[str, str]], dict[str, Any]]:
    primary: dict[tuple[str, str, str], Mapping[str, str]] = {}
    fallback: dict[tuple[str, str, str], Mapping[str, str]] = {}
    duplicate_primary = 0
    duplicate_fallback = 0
    n_rows = 0
    for row in rows:
        n_rows += 1
        source_id = str(row.get("source_id", "") or "")
        contrast_id = str(row.get("contrast_id", "") or "")
        gene_symbol = str(row.get("gene_symbol", "") or "")
        ensembl_id = str(row.get("ensembl_id", "") or "")
        if source_id and contrast_id and gene_symbol:
            key = (source_id, contrast_id, gene_symbol)
            if key in primary:
                duplicate_primary += 1
            else:
                primary[key] = row
        if source_id and contrast_id and ensembl_id:
            key = (source_id, contrast_id, ensembl_id)
            if key in fallback:
                duplicate_fallback += 1
            else:
                fallback[key] = row
    return primary, fallback, {
        "n_reference_rows": n_rows,
        "n_primary_keys": len(primary),
        "n_fallback_keys": len(fallback),
        "duplicate_primary_keys": duplicate_primary,
        "duplicate_fallback_keys": duplicate_fallback,
    }


def _join_effects_to_reference(
    effect_rows: Sequence[Mapping[str, str]],
    reference_rows: Sequence[Mapping[str, str]],
) -> tuple[list[dict[str, Any]], dict[str, int], dict[str, Any]]:
    primary, fallback, reference_summary = _reference_rows_by_key(reference_rows)
    joined: list[dict[str, Any]] = []
    seen: set[tuple[str, str, str, str]] = set()
    counts = {
        "n_effect_rows": len(effect_rows),
        "n_joined_rows": 0,
        "n_primary_join_rows": 0,
        "n_fallback_join_rows": 0,
        "n_unmatched_rows": 0,
        "n_duplicate_effect_keys": 0,
        "n_non_numeric_reference_rows": 0,
        "n_non_numeric_effect_rows": 0,
    }
    for effect in effect_rows:
        effect_value = _as_float(str(effect.get(EFFECT_VALUE_COLUMN, "")))
        if effect_value is None:
            counts["n_non_numeric_effect_rows"] += 1
            continue
        source_id = str(effect.get("source_id", "") or "")
        contrast_id = str(effect.get("contrast_id", "") or "")
        gene_symbol = str(effect.get("gene_symbol", "") or "")
        ensembl_id = str(effect.get("ensembl_id", "") or "")
        dedupe_key = (source_id, contrast_id, gene_symbol, ensembl_id)
        if dedupe_key in seen:
            counts["n_duplicate_effect_keys"] += 1
            continue
        seen.add(dedupe_key)
        reference = None
        join_type = ""
        if gene_symbol:
            reference = primary.get((source_id, contrast_id, gene_symbol))
            join_type = "primary"
        if reference is None and ensembl_id:
            reference = fallback.get((source_id, contrast_id, ensembl_id))
            join_type = "fallback"
        if reference is None:
            counts["n_unmatched_rows"] += 1
            continue
        reference_value = _as_float(str(reference.get(REFERENCE_LOG2FC_COLUMN, "")))
        if reference_value is None:
            counts["n_non_numeric_reference_rows"] += 1
            continue
        significant = str(reference.get("significant_fdr_0_05", "")).lower() == "true"
        joined.append(
            {
                "source_id": source_id,
                "contrast_id": contrast_id,
                "gene_symbol": gene_symbol or str(reference.get("gene_symbol", "") or ""),
                "ensembl_id": ensembl_id or str(reference.get("ensembl_id", "") or ""),
                "effect_value": effect_value,
                "reference_log2fc": reference_value,
                "significant_fdr_0_05": significant,
                "join_type": join_type,
            }
        )
        counts["n_joined_rows"] += 1
        if join_type == "primary":
            counts["n_primary_join_rows"] += 1
        elif join_type == "fallback":
            counts["n_fallback_join_rows"] += 1
    return joined, counts, reference_summary


def _direction_match(joined: Sequence[Mapping[str, Any]]) -> tuple[float | None, dict[str, Any]]:
    significant = [row for row in joined if row.get("significant_fdr_0_05")]
    scored = [
        row
        for row in significant
        if _sign(float(row["effect_value"])) != 0
        and _sign(float(row["reference_log2fc"])) != 0
    ]
    matches = sum(
        _sign(float(row["effect_value"])) == _sign(float(row["reference_log2fc"]))
        for row in scored
    )
    details = {
        "n_joined_rows": len(joined),
        "n_significant_reference_rows": len(significant),
        "n_direction_scored": len(scored),
        "n_direction_matches": matches,
        "reference_filter": "significant_fdr_0_05",
    }
    if not scored:
        return None, details
    return matches / len(scored), details


def _rank_correlation(joined: Sequence[Mapping[str, Any]]) -> tuple[float | None, dict[str, Any]]:
    effect = [float(row["effect_value"]) for row in joined]
    reference = [float(row["reference_log2fc"]) for row in joined]
    value = _spearman(effect, reference)
    return value, {
        "correlation": "spearman",
        "n_joined_rows": len(joined),
        "n_rank_scored": len(effect),
        "reference_filter": "all_joined_rows",
    }


def _feature_key(row: Mapping[str, Any]) -> tuple[str, str, str]:
    return (
        str(row["source_id"]),
        str(row["contrast_id"]),
        str(row["ensembl_id"] or row["gene_symbol"]),
    )


def _log_comb(total: int, selected: int) -> float:
    if selected < 0 or selected > total:
        return float("-inf")
    return (
        math.lgamma(total + 1)
        - math.lgamma(selected + 1)
        - math.lgamma(total - selected + 1)
    )


def _hypergeometric_p_greater_equal(
    *,
    population_size: int,
    success_states: int,
    draws: int,
    observed: int,
) -> float | None:
    """Exact upper-tail probability for overlap under random top-k selection."""

    if population_size <= 0 or success_states < 0 or draws < 0:
        return None
    if success_states > population_size or draws > population_size:
        return None
    min_overlap = max(0, draws - (population_size - success_states))
    max_overlap = min(success_states, draws)
    if observed <= min_overlap:
        return 1.0
    if observed > max_overlap:
        return 0.0
    denominator = _log_comb(population_size, draws)
    log_terms = [
        _log_comb(success_states, overlap)
        + _log_comb(population_size - success_states, draws - overlap)
        - denominator
        for overlap in range(observed, max_overlap + 1)
    ]
    max_log = max(log_terms)
    probability = math.exp(max_log) * sum(
        math.exp(term - max_log) for term in log_terms
    )
    return min(1.0, max(0.0, probability))


def _top_k_overlap(joined: Sequence[Mapping[str, Any]], top_k: Sequence[int]) -> dict[str, Any]:
    significant_keys = {
        _feature_key(row)
        for row in joined
        if row.get("significant_fdr_0_05")
    }
    universe_keys = {_feature_key(row) for row in joined}
    ranked = sorted(
        joined,
        key=lambda row: abs(float(row["effect_value"])),
        reverse=True,
    )
    by_k: dict[str, Any] = {}
    for k in top_k:
        selected = ranked[: min(k, len(ranked))]
        selected_keys = {
            _feature_key(row)
            for row in selected
        }
        overlap = len(selected_keys & significant_keys)
        denominator = len(selected_keys) if selected_keys else 0
        universe_size = len(universe_keys)
        n_significant = len(significant_keys)
        expected = (
            denominator * n_significant / universe_size
            if universe_size and denominator
            else None
        )
        by_k[str(k)] = {
            "k": k,
            "n_selected": denominator,
            "n_feature_universe": universe_size,
            "n_significant_reference_rows": n_significant,
            "n_overlap": overlap,
            "overlap_fraction": (overlap / denominator) if denominator else None,
            "expected_overlap": expected,
            "enrichment": (overlap / expected) if expected else None,
            "hypergeometric_p_value_greater_equal": (
                _hypergeometric_p_greater_equal(
                    population_size=universe_size,
                    success_states=n_significant,
                    draws=denominator,
                    observed=overlap,
                )
                if universe_size and denominator
                else None
            ),
            "null_model": "hypergeometric_random_top_k_without_replacement",
        }
    return by_k


def compute_feature_effect_metrics(
    *,
    feature_effect_path: str | Path,
    reference_signature_path: str | Path,
    task_id: str | None = None,
    top_k: Sequence[int] = (50, 100, 250, 500),
) -> dict[str, Any]:
    """Compute diagnostic feature-effect metrics against the derived DE reference."""

    validation = validate_feature_effect_artifact(feature_effect_path, task_id=task_id)
    skipped_template = {
        "feature_effect_validation": validation,
        "feature_effect_path": Path(feature_effect_path).as_posix(),
        "reference_signature_path": Path(reference_signature_path).as_posix(),
    }
    if not validation["ok"]:
        reason = "; ".join(validation["errors"]) or "invalid feature effect artifact"
        return {
            **skipped_template,
            "metrics": {
                "feature_effect_direction_match": {"status": "skipped", "reason": reason},
                "feature_effect_rank_correlation": {"status": "skipped", "reason": reason},
                "feature_effect_top_k_de_overlap": {"status": "skipped", "reason": reason},
            },
        }

    effect_rows = _read_rows(feature_effect_path)
    reference_rows = _read_rows(reference_signature_path)
    joined, join_counts, reference_summary = _join_effects_to_reference(
        effect_rows,
        reference_rows,
    )
    per_contrast: dict[str, Any] = {}
    for contrast_id in sorted({str(row["contrast_id"]) for row in joined}):
        subset = [row for row in joined if row["contrast_id"] == contrast_id]
        direction, direction_details = _direction_match(subset)
        rank, rank_details = _rank_correlation(subset)
        per_contrast[contrast_id] = {
            "feature_effect_direction_match": direction,
            "feature_effect_direction_match_details": direction_details,
            "feature_effect_rank_correlation": rank,
            "feature_effect_rank_correlation_details": rank_details,
            "feature_effect_top_k_de_overlap": _top_k_overlap(subset, top_k),
            "n_joined_rows": len(subset),
        }

    direction, direction_details = _direction_match(joined)
    rank, rank_details = _rank_correlation(joined)
    common_details = {
        "join": join_counts,
        "reference": {
            **reference_summary,
            "reference_path": Path(reference_signature_path).as_posix(),
        },
        "validation": validation,
        "per_contrast": per_contrast,
    }
    return {
        **skipped_template,
        "metrics": {
            "feature_effect_direction_match": {
                "status": "computed" if direction is not None else "skipped",
                "value": direction,
                "details": {
                    **common_details,
                    "aggregate": direction_details,
                },
            },
            "feature_effect_rank_correlation": {
                "status": "computed" if rank is not None else "skipped",
                "value": rank,
                "details": {
                    **common_details,
                    "aggregate": rank_details,
                },
            },
            "feature_effect_top_k_de_overlap": {
                "status": "computed",
                "value": _top_k_overlap(joined, top_k),
                "details": {
                    **common_details,
                    "aggregate": {
                        "n_joined_rows": len(joined),
                        "reference_filter": "significant_fdr_0_05",
                    },
                },
            },
        },
    }


def _sample_factor_by_id(task: HumanOrganoidTaskData) -> dict[str, dict[str, str]]:
    return {
        str(row["sample_id"]): row
        for row in task.sample_factors
        if row.get("parse_status") == "parsed" and row.get("sample_id")
    }


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
        if training_source != target_source:
            source_folds.append((target_source, training_source, fold))
    if not source_folds:
        raise ValueError("no single-source-to-single-source transfer folds found")
    return sorted(source_folds, key=lambda item: item[0])


def _transform(values: np.ndarray, mode: str) -> np.ndarray:
    if mode == "none":
        return values.astype(np.float64, copy=True)
    if mode == "log1p":
        return np.log1p(np.clip(values.astype(np.float64, copy=True), a_min=0.0, a_max=None))
    raise ValueError(f"unsupported transform: {mode}")


def _select_train_variable_features(
    train: np.ndarray,
    feature_ids: Sequence[str],
    *,
    top_variable_genes: int,
) -> tuple[np.ndarray, list[str], np.ndarray]:
    if top_variable_genes < 1:
        raise ValueError("top_variable_genes must be at least 1")
    if top_variable_genes >= train.shape[1]:
        selected = np.arange(train.shape[1])
    else:
        variances = np.var(train, axis=0)
        selected = np.argsort(-variances, kind="mergesort")[:top_variable_genes]
        selected.sort()
    return train[:, selected], [str(feature_ids[index]) for index in selected], selected


def _scale_train(values: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mean = values.mean(axis=0)
    std = values.std(axis=0)
    std[std == 0.0] = 1.0
    return (values - mean) / std, mean, std


def _label_counts(labels: Sequence[str]) -> dict[str, int]:
    return {
        "n_train_ground": sum(label == NEGATIVE_LABEL for label in labels),
        "n_train_leo_or_iss": sum(label == POSITIVE_LABEL for label in labels),
    }


def _adaptive_pca_components(
    *,
    requested_components: int,
    n_train: int,
    n_features: int,
) -> int:
    if requested_components < 1:
        raise ValueError("pca_components must be at least 1")
    if n_train < 2:
        raise ValueError("PCA-LR feature effects require at least two training samples")
    return max(1, min(requested_components, n_train - 1, n_features))


def build_l2_logistic_feature_effect(
    task: HumanOrganoidTaskData,
    *,
    de_reference_manifest_path: str | Path,
    transform: str = "log1p",
    top_variable_genes: int = 2000,
    max_iter: int = 5000,
    random_state: int = 42,
    c: float = 1.0,
) -> FeatureEffectResult:
    """Fit source-transfer L2 logistic models and emit gene-space coefficients."""

    from sklearn.linear_model import LogisticRegression

    contrasts_by_source = source_contrast_metadata_from_de_manifest(de_reference_manifest_path)
    sample_factors_by_id = _sample_factor_by_id(task)
    feature_ids = [str(feature) for feature in task.features.columns]
    effect_rows: list[dict[str, str]] = []
    fold_summaries: list[dict[str, Any]] = []
    for target_source, training_source, fold in _source_transfer_folds(task):
        contrast_ids = [
            str(contrast["contrast_id"])
            for contrast in contrasts_by_source.get(target_source, [])
        ]
        if not contrast_ids:
            continue
        train_x = task.features.loc[fold.train_sample_ids]
        labels = [sample_factors_by_id[sample_id]["true_label"] for sample_id in train_x.index]
        counts = _label_counts(labels)
        if not counts["n_train_ground"] or not counts["n_train_leo_or_iss"]:
            raise ValueError(f"{fold.fold_id} requires both Ground and LEO_or_ISS training labels")
        y = np.asarray([1 if label == POSITIVE_LABEL else 0 for label in labels], dtype=int)
        transformed = _transform(train_x.to_numpy(dtype=np.float64), transform)
        selected_values, selected_features, selected_indices = _select_train_variable_features(
            transformed,
            feature_ids,
            top_variable_genes=top_variable_genes,
        )
        scaled_values, _, _ = _scale_train(selected_values)
        model = LogisticRegression(
            C=c,
            class_weight="balanced",
            max_iter=max_iter,
            random_state=random_state,
            solver="liblinear",
        )
        model.fit(scaled_values, y)
        classes = list(model.classes_)
        if classes != [0, 1]:
            raise ValueError(f"unexpected logistic classes for {fold.fold_id}: {classes}")
        coefficients = model.coef_[0]
        fold_summaries.append(
            {
                "fold_id": fold.fold_id,
                "target_source_id": target_source,
                "training_source_id": training_source,
                "heldout_factor": fold.heldout_factor,
                "heldout_value": fold.heldout_value,
                **counts,
                "n_features_model_input": int(len(selected_features)),
                "n_target_contrasts": len(contrast_ids),
                "transform": transform,
                "scaling": "train_zscore",
                "feature_selection": f"train_top_{len(selected_features)}_variable_genes",
                "classifier_model_id": LOGISTIC_FEATURE_EFFECT_ID,
                "regularization": "l2",
                "classifier_solver": "liblinear",
                "reference_usage_policy": FEATURE_EFFECT_POLICY,
            }
        )
        for contrast_id in contrast_ids:
            for feature_id, coefficient in zip(selected_features, coefficients):
                effect_rows.append(
                    {
                        "task_id": task.task_id,
                        "fold_id": fold.fold_id,
                        "source_id": target_source,
                        "contrast_id": contrast_id,
                        "feature_namespace": task.feature_namespace or "human_gene",
                        "gene_symbol": "",
                        "ensembl_id": feature_id,
                        "feature_id": feature_id,
                        "effect_value": _format_float(float(coefficient)),
                        "effect_direction_positive_class": POSITIVE_LABEL,
                        "effect_scale": "standardized_logistic_coefficient",
                        "model_space": "gene_space",
                        "classifier_model_id": LOGISTIC_FEATURE_EFFECT_ID,
                        "training_scope": "source_transfer_organoid_type_holdout_train_samples",
                        "target_scope": "target_source_contrast",
                        "positive_label": POSITIVE_LABEL,
                        "reference_usage_policy": FEATURE_EFFECT_POLICY,
                        "training_source_id": training_source,
                        "target_source_id": target_source,
                        "transform": transform,
                        "scaling": "train_zscore",
                        "feature_selection": f"train_top_{len(selected_features)}_variable_genes",
                        "n_train_ground": str(counts["n_train_ground"]),
                        "n_train_leo_or_iss": str(counts["n_train_leo_or_iss"]),
                        "n_features_model_input": str(len(selected_features)),
                        "n_features_emitted": str(len(selected_features)),
                        "regularization": "l2",
                        "classifier_solver": "liblinear",
                        "random_state": str(random_state),
                        "artifact_claim_boundary": "diagnostic_only_not_leaderboard",
                    }
                )

    if not effect_rows:
        raise ValueError("logistic feature-effect adapter produced no rows")
    return FeatureEffectResult(
        classifier_model_id=LOGISTIC_FEATURE_EFFECT_ID,
        effect_rows=effect_rows,
        fold_summaries=fold_summaries,
        target_contrasts_by_source={
            source_id: [str(contrast["contrast_id"]) for contrast in contrasts]
            for source_id, contrasts in sorted(contrasts_by_source.items())
        },
    )


def build_pca_lr_reconstructed_feature_effect(
    task: HumanOrganoidTaskData,
    *,
    de_reference_manifest_path: str | Path,
    transform: str = "log1p",
    top_variable_genes: int = 2000,
    pca_components: int = 50,
    max_iter: int = 5000,
    random_state: int = 42,
    c: float = 1.0,
) -> FeatureEffectResult:
    """Fit source-transfer PCA-LR models and emit reconstructed gene effects."""

    from sklearn.decomposition import PCA
    from sklearn.linear_model import LogisticRegression

    contrasts_by_source = source_contrast_metadata_from_de_manifest(de_reference_manifest_path)
    sample_factors_by_id = _sample_factor_by_id(task)
    feature_ids = [str(feature) for feature in task.features.columns]
    effect_rows: list[dict[str, str]] = []
    fold_summaries: list[dict[str, Any]] = []
    for target_source, training_source, fold in _source_transfer_folds(task):
        contrast_ids = [
            str(contrast["contrast_id"])
            for contrast in contrasts_by_source.get(target_source, [])
        ]
        if not contrast_ids:
            continue
        train_x = task.features.loc[fold.train_sample_ids]
        labels = [sample_factors_by_id[sample_id]["true_label"] for sample_id in train_x.index]
        counts = _label_counts(labels)
        if not counts["n_train_ground"] or not counts["n_train_leo_or_iss"]:
            raise ValueError(f"{fold.fold_id} requires both Ground and LEO_or_ISS training labels")
        y = np.asarray([1 if label == POSITIVE_LABEL else 0 for label in labels], dtype=int)
        transformed = _transform(train_x.to_numpy(dtype=np.float64), transform)
        selected_values, selected_features, _ = _select_train_variable_features(
            transformed,
            feature_ids,
            top_variable_genes=top_variable_genes,
        )
        scaled_values, _, _ = _scale_train(selected_values)
        n_components = _adaptive_pca_components(
            requested_components=pca_components,
            n_train=len(train_x),
            n_features=len(selected_features),
        )
        pca = PCA(n_components=n_components, random_state=random_state, whiten=False)
        pc_scores = pca.fit_transform(scaled_values)
        model = LogisticRegression(
            C=c,
            class_weight="balanced",
            max_iter=max_iter,
            random_state=random_state,
            solver="lbfgs",
        )
        model.fit(pc_scores, y)
        classes = list(model.classes_)
        if classes != [0, 1]:
            raise ValueError(f"unexpected logistic classes for {fold.fold_id}: {classes}")
        pc_coefficients = np.asarray(model.coef_[0], dtype=float)
        gene_coefficients = np.asarray(pca.components_.T @ pc_coefficients, dtype=float)
        if gene_coefficients.shape[0] != len(selected_features):
            raise ValueError(
                "PCA-LR reconstruction length mismatch for "
                f"{fold.fold_id}: {gene_coefficients.shape[0]} vs {len(selected_features)}"
            )
        fold_summaries.append(
            {
                "fold_id": fold.fold_id,
                "target_source_id": target_source,
                "training_source_id": training_source,
                "heldout_factor": fold.heldout_factor,
                "heldout_value": fold.heldout_value,
                **counts,
                "n_features_model_input": int(len(selected_features)),
                "n_target_contrasts": len(contrast_ids),
                "transform": transform,
                "scaling": "train_zscore",
                "feature_selection": f"train_top_{len(selected_features)}_variable_genes",
                "classifier_model_id": PCA_LR_FEATURE_EFFECT_ID,
                "regularization": "l2",
                "classifier_solver": "lbfgs",
                "reference_usage_policy": FEATURE_EFFECT_POLICY,
                "pca_components_requested": int(pca_components),
                "pca_components_used": int(n_components),
                "pca_explained_variance_ratio_sum": float(
                    np.sum(pca.explained_variance_ratio_)
                ),
                "pca_whiten": False,
                "reconstruction_formula": "pca.components_.T @ logistic.coef_[0]",
            }
        )
        for contrast_id in contrast_ids:
            for feature_id, coefficient in zip(selected_features, gene_coefficients):
                effect_rows.append(
                    {
                        "task_id": task.task_id,
                        "fold_id": fold.fold_id,
                        "source_id": target_source,
                        "contrast_id": contrast_id,
                        "feature_namespace": task.feature_namespace or "human_gene",
                        "gene_symbol": "",
                        "ensembl_id": feature_id,
                        "feature_id": feature_id,
                        "effect_value": _format_float(float(coefficient)),
                        "effect_direction_positive_class": POSITIVE_LABEL,
                        "effect_scale": "pca_lr_reconstructed_standardized_logistic_coefficient",
                        "model_space": "reconstructed_gene_space_from_pca",
                        "classifier_model_id": PCA_LR_FEATURE_EFFECT_ID,
                        "training_scope": "source_transfer_organoid_type_holdout_train_samples",
                        "target_scope": "target_source_contrast",
                        "positive_label": POSITIVE_LABEL,
                        "reference_usage_policy": FEATURE_EFFECT_POLICY,
                        "training_source_id": training_source,
                        "target_source_id": target_source,
                        "transform": transform,
                        "scaling": "train_zscore",
                        "feature_selection": f"train_top_{len(selected_features)}_variable_genes",
                        "n_train_ground": str(counts["n_train_ground"]),
                        "n_train_leo_or_iss": str(counts["n_train_leo_or_iss"]),
                        "n_features_model_input": str(len(selected_features)),
                        "n_features_emitted": str(len(selected_features)),
                        "regularization": "l2",
                        "classifier_solver": "lbfgs",
                        "random_state": str(random_state),
                        "artifact_claim_boundary": "diagnostic_only_not_leaderboard",
                        "pca_components_requested": str(pca_components),
                        "pca_components_used": str(n_components),
                        "pca_explained_variance_ratio_sum": _format_float(
                            float(np.sum(pca.explained_variance_ratio_))
                        ),
                        "pca_whiten": "False",
                        "reconstruction_formula": "pca.components_.T @ logistic.coef_[0]",
                    }
                )

    if not effect_rows:
        raise ValueError("PCA-LR feature-effect adapter produced no rows")
    return FeatureEffectResult(
        classifier_model_id=PCA_LR_FEATURE_EFFECT_ID,
        effect_rows=effect_rows,
        fold_summaries=fold_summaries,
        target_contrasts_by_source={
            source_id: [str(contrast["contrast_id"]) for contrast in contrasts]
            for source_id, contrasts in sorted(contrasts_by_source.items())
        },
    )


def write_feature_effect_artifact(
    result: FeatureEffectResult,
    *,
    feature_effect_path: str | Path,
    metadata_path: str | Path,
) -> tuple[Path, Path]:
    """Write feature-effect rows and metadata."""

    effect_path = Path(feature_effect_path)
    meta_path = Path(metadata_path)
    effect_path.parent.mkdir(parents=True, exist_ok=True)
    meta_path.parent.mkdir(parents=True, exist_ok=True)
    optional_columns = [
        "training_source_id",
        "target_source_id",
        "transform",
        "scaling",
        "feature_selection",
        "n_train_ground",
        "n_train_leo_or_iss",
        "n_features_model_input",
        "n_features_emitted",
        "regularization",
        "classifier_solver",
        "random_state",
        "artifact_claim_boundary",
        "pca_components_requested",
        "pca_components_used",
        "pca_explained_variance_ratio_sum",
        "pca_whiten",
        "reconstruction_formula",
    ]
    fieldnames = FEATURE_EFFECT_REQUIRED_COLUMNS + optional_columns
    if effect_path.name.endswith(".gz"):
        handle_context = gzip.open(effect_path, "wt", newline="")
    else:
        handle_context = effect_path.open("w", newline="")
    with handle_context as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in result.effect_rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})
    metadata = {
        "schema_version": "0.1.0",
        "artifact_type": "human_organoid_feature_effect",
        "classifier_model_id": result.classifier_model_id,
        "release_status": "draft_not_frozen",
        "claim_boundary": "diagnostic_only_not_leaderboard",
        "reference_usage_policy": FEATURE_EFFECT_POLICY,
        "n_effect_rows": result.n_effect_rows,
        "fold_summaries": result.fold_summaries,
        "target_contrasts_by_source": result.target_contrasts_by_source,
        "notes": (
            "Feature effects are discriminative classifier coefficients, not log2FC "
            "response signatures. The DE reference is not used during feature-effect "
            "generation and should be used only for post hoc diagnostic scoring."
        ),
    }
    meta_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    return effect_path, meta_path
