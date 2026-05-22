"""Reference-signature audit helpers for the draft human organoid task."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import re
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

from .source_audit import (
    BIODATA_API_BASE,
    FetchResult,
    extract_files_from_listing,
    fetch_url,
    source_file_api_url,
)


DE_REFERENCE_RE = re.compile(r"differential_expression.*\.csv$", re.I)
CONTRAST_RE = re.compile(r"contrasts.*\.csv$", re.I)
ENRICHMENT_RE = re.compile(
    r"(gsea|fgsea|enrich|gene_set|geneset|ontology|pathway|hallmark|signature)",
    re.I,
)
REFERENCE_SIGNATURE_AUDIT_FIELDS = [
    "source_id",
    "glds_prefix",
    "task_ids",
    "api_url",
    "api_status",
    "n_files",
    "candidate_reference_file_count",
    "candidate_reference_files",
    "candidate_reference_categories",
    "de_reference_file_count",
    "de_reference_files",
    "de_reference_urls",
    "de_reference_file_sizes",
    "de_reference_data_types",
    "contrast_file_count",
    "contrast_files",
    "contrast_urls",
    "contrast_sha256",
    "contrast_pair_count",
    "contrast_group_count",
    "contrast_groups",
    "direct_spaceflight_contrast_count",
    "direct_spaceflight_contrasts",
    "reversed_spaceflight_contrast_count",
    "reversed_spaceflight_contrasts",
    "disease_contexts",
    "microglia_conditions",
    "enrichment_reference_file_count",
    "enrichment_reference_files",
    "metric_reference_status",
    "recommended_metric_policy",
    "audit_status",
    "pending_reason",
]


def _join(values: Sequence[str]) -> str:
    return ";".join(sorted({value for value in values if value}))


def _download_url(metadata: Mapping[str, Any]) -> str:
    return str(
        metadata.get("URL")
        or metadata.get("url")
        or metadata.get("download_url")
        or metadata.get("REST_URL")
        or ""
    )


def _rest_url(metadata: Mapping[str, Any]) -> str:
    return str(metadata.get("REST_URL") or metadata.get("rest_url") or "")


def _metadata_value(metadata: Mapping[str, Any], key: str) -> str:
    nested = metadata.get("metadata")
    if isinstance(nested, Mapping) and key in nested:
        return str(nested.get(key, "") or "")
    return str(metadata.get(key, "") or "")


def _merge_file_api_metadata(
    source_id: str,
    filename: str,
    metadata: Mapping[str, Any],
    *,
    fetch: Callable[[str], FetchResult],
) -> Mapping[str, Any]:
    if _metadata_value(metadata, "file_size") and _metadata_value(metadata, "data_type"):
        return metadata
    url = _rest_url(metadata)
    if not url:
        return metadata
    result = fetch(url)
    if not result.ok:
        return metadata
    try:
        listing = json.loads(result.body.decode("utf-8"))
        file_record = extract_files_from_listing(source_id, listing).get(filename, {})
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError):
        return metadata
    merged: dict[str, Any] = dict(metadata)
    if isinstance(file_record, Mapping):
        nested = file_record.get("metadata")
        if isinstance(nested, Mapping):
            merged["metadata"] = {**dict(merged.get("metadata", {})), **dict(nested)}
    return merged


def _is_de_reference(filename: str, metadata: Mapping[str, Any]) -> bool:
    data_type = _metadata_value(metadata, "data_type").lower()
    subcategory = _metadata_value(metadata, "subcategory").lower()
    return bool(DE_REFERENCE_RE.search(filename)) or "differential expression" in data_type or (
        "differential gene expression" in subcategory
    )


def _is_contrast_definition(filename: str) -> bool:
    return bool(CONTRAST_RE.search(filename))


def _is_enrichment_reference(filename: str, metadata: Mapping[str, Any]) -> bool:
    data_type = _metadata_value(metadata, "data_type")
    subcategory = _metadata_value(metadata, "subcategory")
    haystack = " ".join([filename, data_type, subcategory])
    return bool(ENRICHMENT_RE.search(haystack))


def _condition_from_group(group: str) -> dict[str, str]:
    stripped = group.strip().strip("()")
    parts = [part.strip() for part in stripped.split(" & ")]
    if len(parts) != 3:
        return {"disease_context": "", "flight_condition": "", "microglia_condition": ""}
    disease_context, flight_condition, microglia_condition = parts
    return {
        "disease_context": disease_context,
        "flight_condition": flight_condition,
        "microglia_condition": microglia_condition,
    }


def parse_organoid_contrast_table(text: str) -> dict[str, object]:
    """Parse the OSDR bulk-RNA contrast-definition CSV for organoid studies."""

    reader = csv.reader(io.StringIO(text.lstrip("\ufeff")))
    rows = [row for row in reader if any(str(value or "").strip() for value in row)]
    if len(rows) < 3:
        return {
            "contrast_pair_count": 0,
            "contrast_group_count": 0,
            "contrast_groups": [],
            "direct_spaceflight_contrasts": [],
            "reversed_spaceflight_contrasts": [],
            "disease_contexts": [],
            "microglia_conditions": [],
            "parse_status": "missing_rows",
        }
    header, first_groups, second_groups = rows[:3]
    contrasts: list[dict[str, str]] = []
    for index in range(1, len(header)):
        label = str(header[index] or "").strip()
        first = str(first_groups[index] if index < len(first_groups) else "").strip()
        second = str(second_groups[index] if index < len(second_groups) else "").strip()
        if not label:
            continue
        first_display = label.split(")v(", 1)[0].strip("(")
        second_display = label.split(")v(", 1)[1].strip(")") if ")v(" in label else ""
        first_parsed = _condition_from_group(first_display)
        second_parsed = _condition_from_group(second_display)
        contrasts.append(
            {
                "label": label,
                "first_group": first,
                "second_group": second,
                "first_display": first_display,
                "second_display": second_display,
                "first_disease_context": first_parsed["disease_context"],
                "second_disease_context": second_parsed["disease_context"],
                "first_flight_condition": first_parsed["flight_condition"],
                "second_flight_condition": second_parsed["flight_condition"],
                "first_microglia_condition": first_parsed["microglia_condition"],
                "second_microglia_condition": second_parsed["microglia_condition"],
            }
        )

    contrast_groups = sorted(
        {
            value
            for contrast in contrasts
            for value in [contrast["first_display"], contrast["second_display"]]
            if value
        }
    )
    disease_contexts = sorted(
        {
            value
            for contrast in contrasts
            for value in [contrast["first_disease_context"], contrast["second_disease_context"]]
            if value
        }
    )
    microglia_conditions = sorted(
        {
            value
            for contrast in contrasts
            for value in [contrast["first_microglia_condition"], contrast["second_microglia_condition"]]
            if value
        }
    )
    direct: list[str] = []
    reversed_pairs: list[str] = []
    for contrast in contrasts:
        same_disease = contrast["first_disease_context"] == contrast["second_disease_context"]
        same_microglia = contrast["first_microglia_condition"] == contrast["second_microglia_condition"]
        condition_pair = {
            contrast["first_flight_condition"],
            contrast["second_flight_condition"],
        }
        if not same_disease or not same_microglia or condition_pair != {"Ground Control", "Space Flight"}:
            continue
        if (
            contrast["first_flight_condition"] == "Ground Control"
            and contrast["second_flight_condition"] == "Space Flight"
        ):
            direct.append(contrast["label"])
        elif (
            contrast["first_flight_condition"] == "Space Flight"
            and contrast["second_flight_condition"] == "Ground Control"
        ):
            reversed_pairs.append(contrast["label"])

    return {
        "contrast_pair_count": len(contrasts),
        "contrast_group_count": len(contrast_groups),
        "contrast_groups": contrast_groups,
        "direct_spaceflight_contrasts": sorted(direct),
        "reversed_spaceflight_contrasts": sorted(reversed_pairs),
        "disease_contexts": disease_contexts,
        "microglia_conditions": microglia_conditions,
        "parse_status": "parsed",
    }


def _classify_reference_file(filename: str, metadata: Mapping[str, Any]) -> str:
    categories: list[str] = []
    if _is_de_reference(filename, metadata):
        categories.append("differential_expression")
    if _is_contrast_definition(filename):
        categories.append("contrast_definition")
    if _is_enrichment_reference(filename, metadata):
        categories.append("enrichment_or_signature")
    return "+".join(categories)


def audit_organoid_signature_reference_row(
    row: Mapping[str, str],
    *,
    fetcher: Callable[[str], FetchResult] | None = None,
    api_base: str = BIODATA_API_BASE,
) -> dict[str, str]:
    source_id = str(row["source_id"])
    api_url = source_file_api_url(source_id, api_base=api_base)
    fetch = fetcher or (lambda url: fetch_url(url, timeout=60, max_bytes=None))
    listing_result = fetch(api_url)
    output = {field: "" for field in REFERENCE_SIGNATURE_AUDIT_FIELDS}
    output.update(
        {
            "source_id": source_id,
            "glds_prefix": str(row.get("glds_prefix", "")),
            "task_ids": str(row.get("task_ids", "")),
            "api_url": api_url,
        }
    )
    if not listing_result.ok:
        output.update(
            {
                "api_status": "error",
                "audit_status": "api_error",
                "metric_reference_status": "api_error",
                "recommended_metric_policy": "classifier_only_until_reference_audit_succeeds",
                "pending_reason": listing_result.error,
            }
        )
        return output

    try:
        listing = json.loads(listing_result.body.decode("utf-8"))
        files = extract_files_from_listing(source_id, listing)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        output.update(
            {
                "api_status": "invalid_json",
                "audit_status": "api_error",
                "metric_reference_status": "api_error",
                "recommended_metric_policy": "classifier_only_until_reference_audit_succeeds",
                "pending_reason": str(exc),
            }
        )
        return output

    de_files = [
        (filename, metadata)
        for filename, metadata in sorted(files.items())
        if _is_de_reference(filename, metadata)
    ]
    enriched_de_files = [
        (
            filename,
            _merge_file_api_metadata(
                source_id,
                filename,
                metadata,
                fetch=fetch,
            ),
        )
        for filename, metadata in de_files
    ]
    contrast_files = [
        (filename, metadata)
        for filename, metadata in sorted(files.items())
        if _is_contrast_definition(filename)
    ]
    enrichment_files = [
        (filename, metadata)
        for filename, metadata in sorted(files.items())
        if _is_enrichment_reference(filename, metadata)
        and not _is_de_reference(filename, metadata)
        and not _is_contrast_definition(filename)
    ]
    candidate_files = [
        (filename, metadata)
        for filename, metadata in sorted(files.items())
        if _classify_reference_file(filename, metadata)
    ]

    output["api_status"] = "ok"
    output["n_files"] = str(len(files))
    output["candidate_reference_file_count"] = str(len(candidate_files))
    output["candidate_reference_files"] = _join([filename for filename, _ in candidate_files])
    output["candidate_reference_categories"] = _join(
        [_classify_reference_file(filename, metadata) for filename, metadata in candidate_files]
    )
    output["de_reference_file_count"] = str(len(enriched_de_files))
    output["de_reference_files"] = _join([filename for filename, _ in enriched_de_files])
    output["de_reference_urls"] = _join([_download_url(metadata) for _, metadata in enriched_de_files])
    output["de_reference_file_sizes"] = _join(
        [_metadata_value(metadata, "file_size") for _, metadata in enriched_de_files]
    )
    output["de_reference_data_types"] = _join(
        [_metadata_value(metadata, "data_type") for _, metadata in enriched_de_files]
    )
    output["contrast_file_count"] = str(len(contrast_files))
    output["contrast_files"] = _join([filename for filename, _ in contrast_files])
    output["contrast_urls"] = _join([_download_url(metadata) for _, metadata in contrast_files])
    output["enrichment_reference_file_count"] = str(len(enrichment_files))
    output["enrichment_reference_files"] = _join([filename for filename, _ in enrichment_files])

    contrast_hashes: list[str] = []
    contrast_groups: list[str] = []
    direct_contrasts: list[str] = []
    reversed_contrasts: list[str] = []
    disease_contexts: list[str] = []
    microglia_conditions: list[str] = []
    contrast_counts: list[str] = []
    group_counts: list[str] = []
    contrast_errors: list[str] = []
    for filename, metadata in contrast_files:
        url = _download_url(metadata)
        if not url:
            contrast_errors.append(f"{filename}: missing download URL")
            continue
        contrast_result = fetch(url)
        if not contrast_result.ok:
            contrast_errors.append(f"{filename}: {contrast_result.error}")
            continue
        contrast_hashes.append(hashlib.sha256(contrast_result.body).hexdigest())
        parsed = parse_organoid_contrast_table(
            contrast_result.body.decode("utf-8", errors="replace")
        )
        if parsed["parse_status"] != "parsed":
            contrast_errors.append(f"{filename}: contrast table parse_status={parsed['parse_status']}")
            continue
        contrast_counts.append(str(parsed["contrast_pair_count"]))
        group_counts.append(str(parsed["contrast_group_count"]))
        contrast_groups.extend(str(value) for value in parsed["contrast_groups"])
        direct_contrasts.extend(str(value) for value in parsed["direct_spaceflight_contrasts"])
        reversed_contrasts.extend(str(value) for value in parsed["reversed_spaceflight_contrasts"])
        disease_contexts.extend(str(value) for value in parsed["disease_contexts"])
        microglia_conditions.extend(str(value) for value in parsed["microglia_conditions"])

    output["contrast_sha256"] = _join(contrast_hashes)
    output["contrast_pair_count"] = _join(contrast_counts)
    output["contrast_group_count"] = _join(group_counts)
    output["contrast_groups"] = _join(contrast_groups)
    output["direct_spaceflight_contrast_count"] = str(len(set(direct_contrasts)))
    output["direct_spaceflight_contrasts"] = _join(direct_contrasts)
    output["reversed_spaceflight_contrast_count"] = str(len(set(reversed_contrasts)))
    output["reversed_spaceflight_contrasts"] = _join(reversed_contrasts)
    output["disease_contexts"] = _join(disease_contexts)
    output["microglia_conditions"] = _join(microglia_conditions)

    if enriched_de_files and contrast_files and direct_contrasts:
        output.update(
            {
                "audit_status": "reference_tables_listed_contrast_definitions_parsed",
                "metric_reference_status": "public_osdr_de_reference_tables_available_pending_contrast_freeze",
                "recommended_metric_policy": (
                    "keep_classification_primary_enable_de_signature_after_frozen_contrast_subset"
                ),
                "pending_reason": (
                    "Public OSDR differential-expression tables and contrast definitions are listed. "
                    "DE/signature metrics should wait for a frozen contrast subset and sign-orientation "
                    "contract before leaderboard scoring."
                ),
            }
        )
    elif enriched_de_files and contrast_files:
        output.update(
            {
                "audit_status": "reference_tables_listed_contrast_parse_incomplete",
                "metric_reference_status": "de_reference_available_contrast_parse_pending",
                "recommended_metric_policy": "classifier_only_until_contrast_definitions_parse",
                "pending_reason": " | ".join(contrast_errors) or "No direct spaceflight contrast was parsed.",
            }
        )
    elif enriched_de_files:
        output.update(
            {
                "audit_status": "de_reference_listed_missing_contrasts",
                "metric_reference_status": "de_reference_available_contrast_definitions_missing",
                "recommended_metric_policy": "classifier_only_until_contrast_definitions_added",
                "pending_reason": "Differential-expression files are listed, but contrast files were not found.",
            }
        )
    else:
        output.update(
            {
                "audit_status": "no_de_reference_file_listed",
                "metric_reference_status": "no_public_de_reference_table_found",
                "recommended_metric_policy": "classifier_only_until_de_reference_table_added",
                "pending_reason": "OSDR file listing returned no differential-expression reference table.",
            }
        )
    return output


def audit_organoid_signature_reference_inventory(
    rows: Sequence[Mapping[str, str]],
    *,
    fetcher: Callable[[str], FetchResult] | None = None,
    api_base: str = BIODATA_API_BASE,
) -> list[dict[str, str]]:
    return [
        audit_organoid_signature_reference_row(
            row,
            fetcher=fetcher,
            api_base=api_base,
        )
        for row in rows
    ]


def write_organoid_signature_reference_audit(
    rows: Sequence[Mapping[str, str]],
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> tuple[Path, Path]:
    if not rows:
        raise ValueError("cannot write an empty signature-reference audit")
    output_csv = Path(csv_path)
    output_json = Path(json_path)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    normalized = [
        {field: str(row.get(field, "") or "") for field in REFERENCE_SIGNATURE_AUDIT_FIELDS}
        for row in rows
    ]
    with output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=REFERENCE_SIGNATURE_AUDIT_FIELDS)
        writer.writeheader()
        writer.writerows(normalized)
    output_json.write_text(json.dumps(normalized, indent=2, sort_keys=True) + "\n")
    return output_csv, output_json
