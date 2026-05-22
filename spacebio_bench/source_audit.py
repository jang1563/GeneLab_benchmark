"""Checksum-evidence audit helpers for SpaceBio-Bench source inventories."""

from __future__ import annotations

import csv
import hashlib
import json
import re
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence


BIODATA_API_BASE = "https://visualization.osdr.nasa.gov/biodata/api/v2"
MAX_CHECKSUM_MANIFEST_BYTES = 5_000_000
CHECKSUM_NAME_RE = re.compile(r"(md5sum|md5|sha1|sha256|sha512|checksum|checksums)", re.I)
CHECKSUM_VALUE_RE = re.compile(
    r"^(?:[A-Fa-f0-9]{32}|[A-Fa-f0-9]{40}|[A-Fa-f0-9]{64}|[A-Fa-f0-9]{128})$"
)
CHECKSUM_LINE_RE = re.compile(
    r"^\s*([A-Fa-f0-9]{32}|[A-Fa-f0-9]{40}|[A-Fa-f0-9]{64}|[A-Fa-f0-9]{128})\s+[* ]?(.+?)\s*$"
)
REVERSE_CHECKSUM_LINE_RE = re.compile(
    r"^\s*(.+?)\s+([A-Fa-f0-9]{32}|[A-Fa-f0-9]{40}|[A-Fa-f0-9]{64}|[A-Fa-f0-9]{128})\s*$"
)
BSD_CHECKSUM_LINE_RE = re.compile(
    r"^\s*(?:MD5|SHA1|SHA256|SHA512)\s*\((.+?)\)\s*=\s*"
    r"([A-Fa-f0-9]{32}|[A-Fa-f0-9]{40}|[A-Fa-f0-9]{64}|[A-Fa-f0-9]{128})\s*$",
    re.I,
)
SOURCE_CHECKSUM_AUDIT_FIELDS = [
    "source_id",
    "glds_prefix",
    "mission",
    "tissue",
    "task_ids",
    "api_url",
    "api_status",
    "api_response_sha256",
    "n_files",
    "checksum_manifest_count",
    "checksum_manifest_files",
    "checksum_manifest_urls",
    "parsed_checksum_entries",
    "checksum_algorithms",
    "checksum_payload_matches",
    "checksum_manifest_content_sha256",
    "audit_status",
    "freeze_ready",
    "pending_reason",
]


@dataclass(frozen=True)
class FetchResult:
    """Small wrapper around fetched HTTP content."""

    ok: bool
    url: str
    body: bytes = b""
    error: str = ""

    @property
    def sha256(self) -> str:
        if not self.body:
            return ""
        return hashlib.sha256(self.body).hexdigest()


def read_source_inventory(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def fetch_url(url: str, *, timeout: int = 30, max_bytes: int | None = None) -> FetchResult:
    try:
        with urllib.request.urlopen(url, timeout=timeout) as response:
            if max_bytes is None:
                body = response.read()
            else:
                body = response.read(max_bytes + 1)
                if len(body) > max_bytes:
                    return FetchResult(
                        ok=False,
                        url=url,
                        error=f"response exceeded max_bytes={max_bytes}",
                    )
        return FetchResult(ok=True, url=url, body=body)
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        return FetchResult(ok=False, url=url, error=str(exc))


def source_file_api_url(source_id: str, *, api_base: str = BIODATA_API_BASE) -> str:
    return f"{api_base.rstrip('/')}/dataset/{source_id}/files/"


def extract_files_from_listing(source_id: str, listing: Mapping[str, Any]) -> Mapping[str, Mapping[str, Any]]:
    branch = listing.get(source_id, listing)
    if isinstance(branch, Mapping) and isinstance(branch.get("files"), Mapping):
        files = branch["files"]
        return {
            str(filename): dict(metadata) if isinstance(metadata, Mapping) else {}
            for filename, metadata in files.items()
        }
    raise ValueError(f"OSDR file listing for {source_id} does not contain a files mapping")


def is_checksum_manifest(filename: str) -> bool:
    lower = filename.lower()
    if CHECKSUM_NAME_RE.search(lower):
        return True
    return lower.endswith((".md5", ".sha1", ".sha256", ".sha512"))


def checksum_algorithm(checksum: str) -> str:
    length = len(checksum)
    if length == 32:
        return "md5"
    if length == 40:
        return "sha1"
    if length == 64:
        return "sha256"
    if length == 128:
        return "sha512"
    return f"hex{length}"


def _parse_tabular_checksum_line(line: str) -> tuple[str, str] | None:
    fields = [field.strip() for field in line.split("\t")]
    if len(fields) < 3:
        return None
    for index, field in enumerate(fields):
        if not CHECKSUM_VALUE_RE.fullmatch(field):
            continue
        if index > 0 and fields[index - 1]:
            return field, fields[index - 1]
        if index + 1 < len(fields) and fields[index + 1]:
            return field, fields[index + 1]
    return None


def parse_checksum_manifest(text: str) -> list[dict[str, str]]:
    entries: list[dict[str, str]] = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        tabular_match = _parse_tabular_checksum_line(line)
        if tabular_match:
            checksum, filename = tabular_match
        else:
            match = CHECKSUM_LINE_RE.match(line)
            if match:
                checksum, filename = match.groups()
            else:
                bsd_match = BSD_CHECKSUM_LINE_RE.match(line)
                if bsd_match:
                    filename, checksum = bsd_match.groups()
                else:
                    reverse_match = REVERSE_CHECKSUM_LINE_RE.match(line)
                    if not reverse_match:
                        continue
                    filename, checksum = reverse_match.groups()
        entries.append(
            {
                "algorithm": checksum_algorithm(checksum),
                "checksum": checksum.lower(),
                "filename": filename.strip(),
            }
        )
    return entries


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


def _basename(filename: str) -> str:
    return filename.replace("\\", "/").rsplit("/", 1)[-1]


def _payload_name_matches(
    checksum_filename: str,
    payload_names: set[str],
    payload_basenames: set[str],
) -> bool:
    candidate = _basename(checksum_filename)
    if checksum_filename in payload_names or candidate in payload_basenames:
        return True
    return any(
        payload_name.endswith(candidate)
        for payload_name in payload_names | payload_basenames
    )


def audit_source_row(
    row: Mapping[str, str],
    *,
    fetcher: Callable[[str], FetchResult] | None = None,
    api_base: str = BIODATA_API_BASE,
    fetch_manifest_contents: bool = True,
) -> dict[str, str]:
    source_id = str(row["source_id"])
    api_url = source_file_api_url(source_id, api_base=api_base)
    if fetcher is None:
        fetch_listing = lambda url: fetch_url(url, timeout=30, max_bytes=None)
        fetch_manifest = lambda url: fetch_url(
            url,
            timeout=30,
            max_bytes=MAX_CHECKSUM_MANIFEST_BYTES,
        )
    else:
        fetch_listing = fetcher
        fetch_manifest = fetcher
    listing_result = fetch_listing(api_url)
    output = {field: "" for field in SOURCE_CHECKSUM_AUDIT_FIELDS}
    output.update(
        {
            "source_id": source_id,
            "glds_prefix": str(row.get("glds_prefix", "")),
            "mission": str(row.get("mission", "")),
            "tissue": str(row.get("tissue", "")),
            "task_ids": str(row.get("task_ids", "")),
            "api_url": api_url,
            "api_response_sha256": listing_result.sha256,
            "freeze_ready": "false",
        }
    )

    if not listing_result.ok:
        output.update(
            {
                "api_status": "error",
                "audit_status": "api_error",
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
                "pending_reason": str(exc),
            }
        )
        return output

    checksum_files = [
        (filename, metadata)
        for filename, metadata in sorted(files.items())
        if is_checksum_manifest(filename)
    ]
    output["api_status"] = "ok"
    output["n_files"] = str(len(files))
    output["checksum_manifest_count"] = str(len(checksum_files))
    output["checksum_manifest_files"] = _join([filename for filename, _ in checksum_files])
    output["checksum_manifest_urls"] = _join([_download_url(metadata) for _, metadata in checksum_files])

    if not checksum_files:
        output.update(
            {
                "audit_status": "no_checksum_manifest_listed",
                "pending_reason": "OSDR file listing returned no checksum manifest-like files.",
            }
        )
        return output

    parsed_entries: list[dict[str, str]] = []
    manifest_hashes: list[str] = []
    manifest_errors: list[str] = []
    if fetch_manifest_contents:
        for filename, metadata in checksum_files:
            url = _download_url(metadata)
            if not url:
                manifest_errors.append(f"{filename}: missing download URL")
                continue
            manifest_result = fetch_manifest(url)
            if not manifest_result.ok:
                manifest_errors.append(f"{filename}: {manifest_result.error}")
                continue
            manifest_hashes.append(manifest_result.sha256)
            try:
                parsed_entries.extend(
                    parse_checksum_manifest(manifest_result.body.decode("utf-8", errors="replace"))
                )
            except OSError as exc:
                manifest_errors.append(f"{filename}: {exc}")

    output["parsed_checksum_entries"] = str(len(parsed_entries))
    output["checksum_algorithms"] = _join([entry["algorithm"] for entry in parsed_entries])
    output["checksum_manifest_content_sha256"] = _join(manifest_hashes)
    payload_names = set(files)
    payload_basenames = {_basename(filename) for filename in files}
    output["checksum_payload_matches"] = str(
        sum(
            1
            for entry in parsed_entries
            if _payload_name_matches(entry["filename"], payload_names, payload_basenames)
        )
    )

    if parsed_entries:
        output.update(
            {
                "audit_status": "checksum_manifest_parsed",
                "pending_reason": (
                    "Checksum manifest entries were parsed, but payload files have not "
                    "been downloaded and verified in this audit."
                ),
            }
        )
    elif manifest_errors:
        output.update(
            {
                "audit_status": "checksum_manifest_listed_unparsed",
                "pending_reason": " | ".join(manifest_errors),
            }
        )
    else:
        pending_reason = (
            "Checksum manifest-like files were listed but not fetched."
            if not fetch_manifest_contents
            else (
                "Checksum manifest-like files were listed and fetched, but no "
                "checksum entries were parsed."
            )
        )
        output.update(
            {
                "audit_status": "checksum_manifest_listed_unparsed",
                "pending_reason": pending_reason,
            }
        )
    return output


def audit_source_inventory(
    rows: Sequence[Mapping[str, str]],
    *,
    fetcher: Callable[[str], FetchResult] | None = None,
    api_base: str = BIODATA_API_BASE,
    fetch_manifest_contents: bool = True,
) -> list[dict[str, str]]:
    return [
        audit_source_row(
            row,
            fetcher=fetcher,
            api_base=api_base,
            fetch_manifest_contents=fetch_manifest_contents,
        )
        for row in rows
    ]


def write_source_checksum_audit(
    rows: Sequence[Mapping[str, str]],
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> tuple[Path, Path]:
    if not rows:
        raise ValueError("cannot write an empty checksum audit")
    output_csv = Path(csv_path)
    output_json = Path(json_path)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    normalized = [
        {field: str(row.get(field, "") or "") for field in SOURCE_CHECKSUM_AUDIT_FIELDS}
        for row in rows
    ]
    with output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SOURCE_CHECKSUM_AUDIT_FIELDS)
        writer.writeheader()
        writer.writerows(normalized)
    output_json.write_text(json.dumps(normalized, indent=2, sort_keys=True) + "\n")
    return output_csv, output_json
