#!/usr/bin/env python3
"""Validate the public release manifest.

This checker intentionally uses only the Python standard library so it can run
inside lightweight release gates and clean checkouts.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = REPO_ROOT / "release" / "release_manifest.json"
HEX40 = re.compile(r"^[0-9a-f]{40}$")
HEX64 = re.compile(r"^[0-9a-f]{64}$")
LANE_STATUSES = {
    "canonical_historical_result_surface",
    "incubating_translational_extension",
    "metadata_only_public_bulk_alpha",
    "diagnostic_or_draft_lane",
}
CHECKSUM_STATUSES = {"not_recorded", "recorded", "external", "pending"}
GATE_STATUSES = {"planned", "implemented", "partial", "blocked"}


def load_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid JSON: {exc}") from exc


def repo_path(value: str) -> Path | None:
    if not value or value.startswith(("http://", "https://")):
        return None
    if "*" in value:
        return None
    return REPO_ROOT / value


def require_mapping(errors: list[str], where: str, value: Any) -> bool:
    if not isinstance(value, dict):
        errors.append(f"{where}: expected object, got {type(value).__name__}")
        return False
    return True


def require_list(errors: list[str], where: str, value: Any) -> bool:
    if not isinstance(value, list):
        errors.append(f"{where}: expected array, got {type(value).__name__}")
        return False
    return True


def require_fields(
    errors: list[str],
    where: str,
    mapping: dict[str, Any],
    fields: tuple[str, ...],
) -> None:
    for field in fields:
        if field not in mapping:
            errors.append(f"{where}: missing required key {field}")


def check_existing_path(errors: list[str], where: str, value: Any) -> None:
    if not isinstance(value, str):
        errors.append(f"{where}: expected string path")
        return
    path = repo_path(value)
    if path is None:
        return
    if not path.exists():
        errors.append(f"{where}: path does not exist: {value}")


def validate_release_manifest(data: Any) -> list[str]:
    errors: list[str] = []
    if not require_mapping(errors, "manifest", data):
        return errors

    require_fields(
        errors,
        "manifest",
        data,
        (
            "schema_version",
            "generated_at",
            "project",
            "source_control",
            "release_lanes",
            "artifacts",
            "schemas",
            "quality_gates",
            "open_decisions",
        ),
    )

    if data.get("schema_version") != "0.1.0":
        errors.append("manifest.schema_version: expected 0.1.0")

    project = data.get("project")
    if require_mapping(errors, "project", project):
        require_fields(
            errors,
            "project",
            project,
            ("name", "legacy_name", "repository", "huggingface_dataset"),
        )
        repository = project.get("repository")
        if not isinstance(repository, str) or not repository.startswith("https://github.com/"):
            errors.append("project.repository: expected GitHub URL")
        hf_dataset = project.get("huggingface_dataset")
        if not isinstance(hf_dataset, str) or "huggingface.co/datasets/" not in hf_dataset:
            errors.append("project.huggingface_dataset: expected Hugging Face dataset URL")

    source_control = data.get("source_control")
    if require_mapping(errors, "source_control", source_control):
        commit = source_control.get("git_commit_at_manifest_creation")
        if not isinstance(commit, str) or not HEX40.match(commit):
            errors.append("source_control.git_commit_at_manifest_creation: expected 40-char git SHA")

    lanes = data.get("release_lanes")
    lane_ids: set[str] = set()
    if require_list(errors, "release_lanes", lanes):
        for index, lane in enumerate(lanes):
            where = f"release_lanes[{index}]"
            if not require_mapping(errors, where, lane):
                continue
            require_fields(
                errors,
                where,
                lane,
                ("lane_id", "public_label", "status", "summary", "canonical_surfaces"),
            )
            lane_id = lane.get("lane_id")
            if not isinstance(lane_id, str) or not lane_id:
                errors.append(f"{where}.lane_id: expected non-empty string")
            elif lane_id in lane_ids:
                errors.append(f"{where}.lane_id: duplicate lane id {lane_id}")
            else:
                lane_ids.add(lane_id)
            if lane.get("status") not in LANE_STATUSES:
                errors.append(f"{where}.status: unsupported status {lane.get('status')!r}")
            surfaces = lane.get("canonical_surfaces")
            if require_list(errors, f"{where}.canonical_surfaces", surfaces):
                for surface_index, surface in enumerate(surfaces):
                    check_existing_path(
                        errors,
                        f"{where}.canonical_surfaces[{surface_index}]",
                        surface,
                    )

    artifacts = data.get("artifacts")
    artifact_ids: set[str] = set()
    if require_list(errors, "artifacts", artifacts):
        for index, artifact in enumerate(artifacts):
            where = f"artifacts[{index}]"
            if not require_mapping(errors, where, artifact):
                continue
            require_fields(
                errors,
                where,
                artifact,
                ("artifact_id", "lane_id", "role", "path", "kind", "checksum_status"),
            )
            artifact_id = artifact.get("artifact_id")
            if not isinstance(artifact_id, str) or not artifact_id:
                errors.append(f"{where}.artifact_id: expected non-empty string")
            elif artifact_id in artifact_ids:
                errors.append(f"{where}.artifact_id: duplicate artifact id {artifact_id}")
            else:
                artifact_ids.add(artifact_id)
            lane_id = artifact.get("lane_id")
            if isinstance(lane_id, str) and lane_ids and lane_id not in lane_ids:
                errors.append(f"{where}.lane_id: unknown lane id {lane_id}")
            if artifact.get("checksum_status") not in CHECKSUM_STATUSES:
                errors.append(
                    f"{where}.checksum_status: unsupported value "
                    f"{artifact.get('checksum_status')!r}"
                )
            check_existing_path(errors, f"{where}.path", artifact.get("path"))
            sha256 = artifact.get("sha256")
            if sha256 is not None and (not isinstance(sha256, str) or not HEX64.match(sha256)):
                errors.append(f"{where}.sha256: expected 64-char hex digest")

    schemas = data.get("schemas")
    if require_list(errors, "schemas", schemas):
        for index, schema in enumerate(schemas):
            where = f"schemas[{index}]"
            if not require_mapping(errors, where, schema):
                continue
            require_fields(errors, where, schema, ("schema_id", "path", "applies_to"))
            check_existing_path(errors, f"{where}.path", schema.get("path"))

    gates = data.get("quality_gates")
    if require_list(errors, "quality_gates", gates):
        for index, gate in enumerate(gates):
            where = f"quality_gates[{index}]"
            if not require_mapping(errors, where, gate):
                continue
            require_fields(errors, where, gate, ("gate_id", "status", "description"))
            if gate.get("status") not in GATE_STATUSES:
                errors.append(f"{where}.status: unsupported value {gate.get('status')!r}")

    decisions = data.get("open_decisions")
    if require_list(errors, "open_decisions", decisions):
        for index, decision in enumerate(decisions):
            where = f"open_decisions[{index}]"
            if not require_mapping(errors, where, decision):
                continue
            require_fields(
                errors,
                where,
                decision,
                ("decision_id", "question", "recommended_default"),
            )

    return errors


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Validate release/release_manifest.json")
    parser.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_MANIFEST,
        help=f"Release manifest path (default: {DEFAULT_MANIFEST.relative_to(REPO_ROOT)})",
    )
    args = parser.parse_args(argv)

    try:
        data = load_json(args.manifest)
    except ValueError as exc:
        print(f"[FAIL] {args.manifest}: {exc}", file=sys.stderr)
        return 1

    errors = validate_release_manifest(data)
    if errors:
        print(f"[FAIL] {args.manifest}: {len(errors)} issue(s)", file=sys.stderr)
        for error in errors:
            print(f"  - {error}", file=sys.stderr)
        return 1

    print(f"[OK] {args.manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
