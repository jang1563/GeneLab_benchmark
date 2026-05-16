#!/usr/bin/env python3
"""Validate v8 provenance manifests and beta-release metadata.

This checker intentionally uses only the Python standard library so it can run
inside the lightweight HPC release gate before optional dependencies are
installed.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
RUNS_DIR = REPO_ROOT / "v8" / "provenance" / "runs"
INPUT_FREEZE = REPO_ROOT / "v8" / "provenance" / "input_freeze.json"
ARTIFACT_MANIFEST = REPO_ROOT / "v8" / "release" / "v8_beta_artifact_manifest.json"

TOP_KEYS = {
    "schema_version",
    "result_id",
    "pillar",
    "status",
    "claim",
    "generated_at",
    "code",
    "inputs",
    "outputs",
    "environment",
    "metrics",
    "validation",
}
CLAIM_KEYS = {"short", "scope", "consumed_by"}
CODE_KEYS = {"script", "command", "git_commit", "parameters"}
VALIDATION_KEYS = {"local_tests_run", "hpc_tests_run", "checks", "notes"}
ARTIFACT_KEYS = {
    "role",
    "path_or_accession",
    "kind",
    "version",
    "checksum",
    "tracked_in_git",
    "notes",
}
PILLARS = {"bridge", "decompose", "intervene", "causal", "multiomics", "summary"}
STATUSES = {"draft", "hpc_validated", "release_candidate", "frozen"}
SCOPES = {"exploratory", "benchmark", "mechanistic", "release", "reproducibility"}
KINDS = {
    "local_file",
    "osdr_accession",
    "geo_accession",
    "soma_file",
    "url",
    "api",
    "hpc_local_cache",
    "other",
}
HEX64 = re.compile(r"^[0-9a-f]{64}$")


def load_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid JSON: {exc}") from exc


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def is_git_tracked(path: Path) -> bool:
    rel = path.relative_to(REPO_ROOT).as_posix()
    proc = subprocess.run(
        ["git", "ls-files", "--error-unmatch", rel],
        cwd=REPO_ROOT,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return proc.returncode == 0


def exact_repo_path(value: str) -> Path | None:
    if not value:
        return None
    if value.startswith(("$", "http://", "https://")):
        return None
    if any(token in value for token in ("*", "?", "[", "]", ";", "{", "}")):
        return None
    if value.startswith(("OSD-", "GLDS-", "GSE", "GSM")):
        return None
    path = REPO_ROOT / value
    return path if path.exists() else None


def parse_datetime(value: str) -> bool:
    try:
        dt.datetime.fromisoformat(value.replace("Z", "+00:00"))
        return True
    except ValueError:
        return False


def expect_type(errors: list[str], where: str, value: Any, typ: type | tuple[type, ...]) -> bool:
    if not isinstance(value, typ):
        if isinstance(typ, tuple):
            names = ", ".join(t.__name__ for t in typ)
        else:
            names = typ.__name__
        errors.append(f"{where}: expected {names}, got {type(value).__name__}")
        return False
    return True


def validate_artifact(
    errors: list[str],
    warnings: list[str],
    manifest: Path,
    where: str,
    artifact: Any,
    output: bool,
) -> None:
    if not expect_type(errors, where, artifact, dict):
        return
    extra = set(artifact) - ARTIFACT_KEYS
    if extra:
        errors.append(f"{where}: unexpected artifact keys {sorted(extra)}")
    for key in ("role", "path_or_accession"):
        if key not in artifact:
            errors.append(f"{where}: missing required key {key}")
        elif not isinstance(artifact[key], str):
            errors.append(f"{where}.{key}: expected string")
    if "kind" in artifact and artifact["kind"] not in KINDS:
        errors.append(f"{where}.kind: unsupported value {artifact['kind']!r}")
    if "tracked_in_git" in artifact and not isinstance(artifact["tracked_in_git"], bool):
        errors.append(f"{where}.tracked_in_git: expected boolean")

    path_value = artifact.get("path_or_accession")
    if not isinstance(path_value, str):
        return
    local_path = exact_repo_path(path_value)

    if output and artifact.get("tracked_in_git") is True and local_path is not None:
        if not is_git_tracked(local_path):
            errors.append(f"{where}: tracked_in_git=true but file is not tracked: {path_value}")
    if output and artifact.get("tracked_in_git") is True and local_path is None:
        if not any(token in path_value for token in ("*", ";", "{", "}", "$")):
            errors.append(f"{where}: tracked exact output path does not exist: {path_value}")

    checksum = artifact.get("checksum")
    if checksum is None:
        if output and artifact.get("tracked_in_git") is True and local_path is not None and local_path.is_file():
            warnings.append(f"{manifest.name} {where}: exact tracked output has no checksum")
        return
    if not isinstance(checksum, str):
        errors.append(f"{where}.checksum: expected string")
        return
    if not checksum.startswith("sha256:"):
        errors.append(f"{where}.checksum: expected sha256: prefix")
        return
    digest = checksum.removeprefix("sha256:")
    if HEX64.match(digest):
        if local_path is not None and local_path.is_file():
            actual = sha256_file(local_path)
            if actual != digest:
                errors.append(
                    f"{where}: checksum mismatch for {path_value}: expected {digest}, got {actual}"
                )
    elif ";" in digest or "=" in digest:
        # Compound checksums are used for grouped API-derived outputs.
        return
    else:
        errors.append(f"{where}.checksum: malformed sha256 digest")


def validate_run_manifest(path: Path) -> tuple[list[str], list[str]]:
    errors: list[str] = []
    warnings: list[str] = []
    try:
        data = load_json(path)
    except ValueError as exc:
        return [f"{path}: {exc}"], []

    if not expect_type(errors, path.name, data, dict):
        return errors, warnings
    extra = set(data) - TOP_KEYS
    if extra:
        errors.append(f"{path.name}: unexpected top-level keys {sorted(extra)}")

    required = {
        "schema_version",
        "result_id",
        "pillar",
        "status",
        "claim",
        "code",
        "inputs",
        "outputs",
        "validation",
    }
    for key in sorted(required):
        if key not in data:
            errors.append(f"{path.name}: missing required key {key}")

    if data.get("schema_version") != "0.1.0":
        errors.append(f"{path.name}: schema_version must be 0.1.0")
    if data.get("pillar") not in PILLARS:
        errors.append(f"{path.name}: unsupported pillar {data.get('pillar')!r}")
    if data.get("status") not in STATUSES:
        errors.append(f"{path.name}: unsupported status {data.get('status')!r}")
    if data.get("status") in {"hpc_validated", "release_candidate", "frozen"}:
        generated = data.get("generated_at")
        if not isinstance(generated, str) or not parse_datetime(generated):
            errors.append(f"{path.name}: validated manifests require ISO generated_at")

    claim = data.get("claim")
    if expect_type(errors, f"{path.name}.claim", claim, dict):
        extra = set(claim) - CLAIM_KEYS
        if extra:
            errors.append(f"{path.name}.claim: unexpected keys {sorted(extra)}")
        if not isinstance(claim.get("short"), str):
            errors.append(f"{path.name}.claim.short: expected string")
        if claim.get("scope") not in SCOPES:
            errors.append(f"{path.name}.claim.scope: unsupported value {claim.get('scope')!r}")
        if "consumed_by" in claim and not isinstance(claim["consumed_by"], list):
            errors.append(f"{path.name}.claim.consumed_by: expected list")

    code = data.get("code")
    if expect_type(errors, f"{path.name}.code", code, dict):
        extra = set(code) - CODE_KEYS
        if extra:
            errors.append(f"{path.name}.code: unexpected keys {sorted(extra)}")
        for key in ("script", "command"):
            if not isinstance(code.get(key), str):
                errors.append(f"{path.name}.code.{key}: expected string")

    validation = data.get("validation")
    if expect_type(errors, f"{path.name}.validation", validation, dict):
        extra = set(validation) - VALIDATION_KEYS
        if extra:
            errors.append(f"{path.name}.validation: unexpected keys {sorted(extra)}")
        for key in ("local_tests_run", "hpc_tests_run"):
            if not isinstance(validation.get(key), bool):
                errors.append(f"{path.name}.validation.{key}: expected boolean")
        if not isinstance(validation.get("notes"), str):
            errors.append(f"{path.name}.validation.notes: expected string")

    for section, output in (("inputs", False), ("outputs", True)):
        items = data.get(section)
        if not expect_type(errors, f"{path.name}.{section}", items, list):
            continue
        for idx, artifact in enumerate(items):
            validate_artifact(errors, warnings, path, f"{path.name}.{section}[{idx}]", artifact, output)

    text = path.read_text()
    if "TBD_HPC_RUN" in text:
        errors.append(f"{path.name}: contains TBD_HPC_RUN")
    if data.get("status") == "draft":
        errors.append(f"{path.name}: remains draft")
    return errors, warnings


def validate_input_freeze(require_frozen: bool) -> tuple[list[str], list[str]]:
    errors: list[str] = []
    warnings: list[str] = []
    if not INPUT_FREEZE.exists():
        return [f"{INPUT_FREEZE.relative_to(REPO_ROOT)}: missing"], []
    data = load_json(INPUT_FREEZE)
    required = {
        "schema_version",
        "freeze_id",
        "status",
        "generated_at",
        "repository_commit",
        "external_sources",
        "release_blockers",
    }
    for key in sorted(required):
        if key not in data:
            errors.append(f"input_freeze.json: missing {key}")
    if data.get("schema_version") != "0.1.0":
        errors.append("input_freeze.json: schema_version must be 0.1.0")
    if data.get("status") not in {"release_candidate", "frozen"}:
        errors.append("input_freeze.json: status must be release_candidate or frozen")
    if require_frozen and data.get("status") != "frozen":
        errors.append("input_freeze.json: --require-frozen requires status=frozen")
    if not isinstance(data.get("external_sources"), list) or not data.get("external_sources"):
        errors.append("input_freeze.json: external_sources must be a non-empty list")
    blockers = data.get("release_blockers", [])
    if not isinstance(blockers, list):
        errors.append("input_freeze.json: release_blockers must be a list")
    elif blockers:
        warnings.append(f"input_freeze.json: {len(blockers)} beta freeze blocker(s) remain")
        if require_frozen:
            errors.append("input_freeze.json: frozen release cannot have open blockers")
    return errors, warnings


def validate_artifact_manifest(require_frozen: bool) -> tuple[list[str], list[str]]:
    errors: list[str] = []
    warnings: list[str] = []
    if not ARTIFACT_MANIFEST.exists():
        return [f"{ARTIFACT_MANIFEST.relative_to(REPO_ROOT)}: missing"], []
    data = load_json(ARTIFACT_MANIFEST)
    required = {
        "schema_version",
        "bundle_id",
        "status",
        "generated_at",
        "git_commit",
        "targets",
        "files",
        "release_blockers",
    }
    for key in sorted(required):
        if key not in data:
            errors.append("v8_beta_artifact_manifest.json: missing " + key)
    if data.get("schema_version") != "0.1.0":
        errors.append("v8_beta_artifact_manifest.json: schema_version must be 0.1.0")
    if data.get("status") not in {"release_candidate", "frozen"}:
        errors.append("v8_beta_artifact_manifest.json: unsupported status")
    if require_frozen and data.get("status") != "frozen":
        errors.append("v8_beta_artifact_manifest.json: --require-frozen requires status=frozen")

    files = data.get("files")
    if isinstance(files, list):
        for idx, item in enumerate(files):
            if not isinstance(item, dict):
                errors.append(f"v8_beta_artifact_manifest.json.files[{idx}]: expected object")
                continue
            path_value = item.get("path")
            if not isinstance(path_value, str):
                errors.append(f"v8_beta_artifact_manifest.json.files[{idx}].path: expected string")
                continue
            local_path = exact_repo_path(path_value)
            if local_path is not None and item.get("tracked_in_git") is True and not is_git_tracked(local_path):
                errors.append(f"artifact manifest: {path_value} marked tracked but is not tracked")
            checksum = item.get("sha256")
            if checksum and local_path is not None and local_path.is_file():
                if sha256_file(local_path) != checksum:
                    errors.append(f"artifact manifest: checksum mismatch for {path_value}")
    else:
        errors.append("v8_beta_artifact_manifest.json: files must be a list")

    blockers = data.get("release_blockers", [])
    if not isinstance(blockers, list):
        errors.append("v8_beta_artifact_manifest.json: release_blockers must be a list")
    elif blockers:
        warnings.append(f"v8_beta_artifact_manifest.json: {len(blockers)} packaging blocker(s) remain")
        if require_frozen:
            errors.append("v8_beta_artifact_manifest.json: frozen release cannot have open blockers")
    return errors, warnings


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--require-frozen",
        action="store_true",
        help="Require input and artifact manifests to be frozen with no open blockers.",
    )
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    errors: list[str] = []
    warnings: list[str] = []

    manifests = sorted(RUNS_DIR.glob("*.json"))
    if not manifests:
        errors.append("v8/provenance/runs: no manifest JSON files found")
    for manifest in manifests:
        e, w = validate_run_manifest(manifest)
        errors.extend(e)
        warnings.extend(w)

    for validator in (validate_input_freeze, validate_artifact_manifest):
        e, w = validator(args.require_frozen)
        errors.extend(e)
        warnings.extend(w)

    if not args.quiet:
        for warning in warnings:
            print(f"WARNING: {warning}", file=sys.stderr)
    if errors:
        print("v8 provenance validation failed:", file=sys.stderr)
        for error in errors:
            print(f"  - {error}", file=sys.stderr)
        return 1
    if not args.quiet:
        print(f"v8 provenance validation OK ({len(manifests)} run manifests)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
