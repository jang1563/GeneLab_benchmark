"""Mission-discrimination metric for spaceflight domain-shift evaluation."""

from __future__ import annotations

from collections import defaultdict
from math import isfinite, sqrt
from typing import Hashable, Iterable, Sequence


def _as_matrix(embeddings: Iterable[Sequence[float]]) -> list[list[float]]:
    rows: list[list[float]] = []
    expected_dim: int | None = None
    for row_index, row in enumerate(embeddings):
        values = [float(value) for value in row]
        if not values:
            raise ValueError(f"embedding row {row_index} is empty")
        if any(not isfinite(value) for value in values):
            raise ValueError(f"embedding row {row_index} contains non-finite values")
        if expected_dim is None:
            expected_dim = len(values)
        elif len(values) != expected_dim:
            raise ValueError("all embedding rows must have the same dimension")
        rows.append(values)
    if not rows:
        raise ValueError("embeddings must contain at least one row")
    return rows


def _centroid(rows: list[list[float]], indices: list[int]) -> list[float]:
    if not indices:
        raise ValueError("cannot compute a centroid from no rows")
    dim = len(rows[0])
    return [sum(rows[index][j] for index in indices) / len(indices) for j in range(dim)]


def _euclidean(left: Sequence[float], right: Sequence[float]) -> float:
    return sqrt(sum((a - b) ** 2 for a, b in zip(left, right)))


def _cosine(left: Sequence[float], right: Sequence[float]) -> float:
    dot = sum(a * b for a, b in zip(left, right))
    left_norm = sqrt(sum(a * a for a in left))
    right_norm = sqrt(sum(b * b for b in right))
    if left_norm == 0.0 and right_norm == 0.0:
        return 0.0
    if left_norm == 0.0 or right_norm == 0.0:
        return 1.0
    return 1.0 - dot / (left_norm * right_norm)


def _distance(left: Sequence[float], right: Sequence[float], distance: str) -> float:
    if distance == "euclidean":
        return _euclidean(left, right)
    if distance == "cosine":
        return _cosine(left, right)
    raise ValueError("distance must be 'euclidean' or 'cosine'")


def mission_discrimination_score(
    embeddings: Iterable[Sequence[float]],
    mission_labels: Iterable[Hashable],
    *,
    sample_ids: Iterable[str] | None = None,
    distance: str = "euclidean",
    tie_tolerance: float = 1e-12,
) -> dict[str, object]:
    """Score whether samples are closest to their own mission centroid.

    For each sample, the metric ranks the distance from that sample to every
    mission centroid. The sample's own mission centroid is computed in
    leave-one-sample-out fashion to avoid self-nearest-neighbor leakage. The
    per-sample score is 1.0 when its own mission ranks first, 0.0 when it ranks
    last, and 0.5 under a complete tie. The returned aggregate is the mean of
    all scored samples.
    """

    rows = _as_matrix(embeddings)
    labels = list(mission_labels)
    if len(rows) != len(labels):
        raise ValueError("embeddings and mission_labels must have the same length")
    if len(set(labels)) < 2:
        raise ValueError("mission_discrimination requires at least two missions")

    ids = list(sample_ids) if sample_ids is not None else [str(i) for i in range(len(rows))]
    if len(ids) != len(rows):
        raise ValueError("sample_ids must have the same length as embeddings")

    label_to_indices: dict[Hashable, list[int]] = defaultdict(list)
    for index, label in enumerate(labels):
        label_to_indices[label].append(index)

    per_sample: list[dict[str, object]] = []
    skipped = 0
    for index, own_label in enumerate(labels):
        candidate_distances: dict[Hashable, float] = {}
        for label, indices in label_to_indices.items():
            centroid_indices = [i for i in indices if i != index] if label == own_label else indices
            if not centroid_indices:
                continue
            candidate_distances[label] = _distance(
                rows[index],
                _centroid(rows, centroid_indices),
                distance,
            )

        if own_label not in candidate_distances or len(candidate_distances) < 2:
            skipped += 1
            continue

        own_distance = candidate_distances[own_label]
        less = sum(
            dist < own_distance - tie_tolerance for dist in candidate_distances.values()
        )
        tied = sum(
            abs(dist - own_distance) <= tie_tolerance
            for dist in candidate_distances.values()
        )
        average_rank = less + (tied + 1) / 2.0
        denom = len(candidate_distances) - 1
        normalized = 1.0 - (average_rank - 1.0) / denom

        per_sample.append(
            {
                "sample_id": ids[index],
                "sample_index": index,
                "mission": str(own_label),
                "own_distance": own_distance,
                "rank": average_rank,
                "score": normalized,
                "n_candidate_missions": len(candidate_distances),
            }
        )

    if not per_sample:
        raise ValueError("no samples could be scored; each scored mission needs at least two samples")

    aggregate = sum(float(item["score"]) for item in per_sample) / len(per_sample)
    return {
        "metric_id": "mission_discrimination",
        "distance": distance,
        "score": aggregate,
        "n_scored": len(per_sample),
        "n_skipped": skipped,
        "per_sample": per_sample,
    }
