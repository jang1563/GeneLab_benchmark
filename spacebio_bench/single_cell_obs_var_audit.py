"""Single-cell AnnData obs/var audit helpers for SpaceBio-Bench v9."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any


DEFAULT_SC_OBS_VAR_AUDIT_ID = "v9_sc_005_rrrm1_blood_obs_var_audit"

SC_OBS_VAR_AUDIT_SUMMARY_FIELDS = [
    "audit_id",
    "decision_status",
    "task_id",
    "canonical_payload_path",
    "payload_path_status",
    "payload_sha256",
    "payload_size_bytes",
    "n_obs",
    "n_vars",
    "requirement_count",
    "pass_count",
    "fail_count",
    "skip_count",
    "blocker_count",
    "next_required_block",
    "claim_boundary",
]

SC_OBS_VAR_AUDIT_RESULT_FIELDS = [
    "audit_id",
    "field_id",
    "axis",
    "requirement_status",
    "audit_status",
    "observed_presence",
    "observed_completeness",
    "observed_summary",
    "blocker_if_missing",
    "audit_message",
]

SC_OBS_VAR_PAYLOAD_MANIFEST_FIELDS = [
    "audit_id",
    "payload_id",
    "payload_path",
    "path_status",
    "sha256",
    "size_bytes",
    "n_obs",
    "n_vars",
    "h5ad_read_status",
]


def _resolve_path(path: str | Path, repo_root: Path) -> Path:
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    return repo_root / candidate


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _requirement_rows(requirements_path: Path) -> list[dict[str, str]]:
    rows = _read_csv_rows(requirements_path)
    if not rows:
        raise ValueError(f"empty obs/var audit requirements table: {requirements_path}")
    return rows


def _missing_payload_results(
    *,
    audit_id: str,
    requirement_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    return [
        {
            "audit_id": audit_id,
            "field_id": row["field_id"],
            "axis": row["axis"],
            "requirement_status": row["requirement_status"],
            "audit_status": "skipped_no_local_payload",
            "observed_presence": "NA",
            "observed_completeness": "NA",
            "observed_summary": "NA",
            "blocker_if_missing": row["blocker_if_missing"],
            "audit_message": "Canonical AnnData payload is absent; field audit not run.",
        }
        for row in requirement_rows
    ]


def _load_anndata(path: Path) -> tuple[Any | None, str]:
    try:
        import anndata as ad  # type: ignore
    except ImportError:
        return None, "blocked_missing_anndata_dependency"
    try:
        return ad.read_h5ad(path), "read_ok"
    except Exception as exc:  # pragma: no cover - exercised only with bad h5ad files.
        return None, f"read_failed:{type(exc).__name__}"


def _presence_from_axis(adata: Any, row: dict[str, str]) -> tuple[str, str, str, str]:
    field = row["field_id"]
    axis = row["axis"]
    if axis == "obs":
        if field == "cell_id":
            count = len(adata.obs_names)
            unique = len(set(map(str, adata.obs_names))) == count
            status = "pass" if count and unique else "fail"
            return status, "present", "complete", f"obs_names={count};unique={unique}"
        if field in adata.obs.columns:
            series = adata.obs[field]
            completeness = float(series.notna().mean()) if len(series) else 0.0
            status = "pass" if completeness == 1.0 else "fail"
            return status, "present", f"{completeness:.6f}", f"n_values={len(series)}"
        aliases = [alias for alias in row["accepted_aliases"].split("|") if alias]
        present_aliases = [alias for alias in aliases if alias in adata.obs.columns]
        if present_aliases:
            return (
                "fail",
                "alias_present_requires_normalization",
                "NA",
                "|".join(present_aliases),
            )
        return "fail", "missing", "0.000000", "obs column absent"
    if axis == "var":
        if field in adata.var.columns:
            series = adata.var[field]
            completeness = float(series.notna().mean()) if len(series) else 0.0
            status = "pass" if completeness == 1.0 else "fail"
            return status, "present", f"{completeness:.6f}", f"n_values={len(series)}"
        if field in {"gene_symbol", "feature_id"} and len(adata.var_names):
            return (
                "fail",
                "alias_present_requires_normalization",
                "NA",
                f"var_names={len(adata.var_names)}",
            )
        return "fail", "missing", "0.000000", "var field absent"
    if axis == "uns":
        if field in adata.uns:
            return "pass", "present", "complete", "uns key present"
        return "fail", "missing", "0.000000", "uns key absent"
    if axis == "matrix" and field == "X":
        shape = getattr(adata.X, "shape", None)
        if shape == (adata.n_obs, adata.n_vars):
            return "pass", "present", "complete", f"shape={shape[0]}x{shape[1]}"
        return "fail", "missing_or_mismatched", "NA", "X shape missing or mismatched"
    if axis == "layers" and field == "layers.counts":
        if "counts" in adata.layers:
            layer = adata.layers["counts"]
            shape = getattr(layer, "shape", None)
            return "pass", "present", "complete", f"shape={shape[0]}x{shape[1]}"
        return "fail", "missing", "0.000000", "counts layer absent"
    if axis == "raw":
        if adata.raw is not None:
            return "pass", "present", "complete", f"shape={adata.raw.shape[0]}x{adata.raw.shape[1]}"
        return "pass", "missing_optional", "NA", "raw object absent"
    return "fail", "unsupported_axis", "NA", f"unsupported axis {axis}"


def _audit_payload_results(
    *,
    audit_id: str,
    requirement_rows: list[dict[str, str]],
    adata: Any,
) -> list[dict[str, str]]:
    results: list[dict[str, str]] = []
    for row in requirement_rows:
        status, presence, completeness, summary = _presence_from_axis(adata, row)
        results.append(
            {
                "audit_id": audit_id,
                "field_id": row["field_id"],
                "axis": row["axis"],
                "requirement_status": row["requirement_status"],
                "audit_status": status,
                "observed_presence": presence,
                "observed_completeness": completeness,
                "observed_summary": summary,
                "blocker_if_missing": row["blocker_if_missing"],
                "audit_message": "Field audit executed against local AnnData payload.",
            }
        )
    return results


def _is_blocker(row: dict[str, str]) -> bool:
    blocker = row["blocker_if_missing"]
    return blocker not in {
        "no_blocker",
        "warn_not_block",
        "warn_not_block_unless_needed_by_metric",
    }


def build_sc_obs_var_audit(
    *,
    repo_root: str | Path = ".",
    audit_id: str = DEFAULT_SC_OBS_VAR_AUDIT_ID,
    manifest_path: str | Path = (
        "v9/sc_spaceflight/task_manifests/"
        "draft_rrrm1_blood_single_cell_spaceflight.json"
    ),
    requirements_path: str | Path = "v9/sc_spaceflight/obs_var_audit_requirements.csv",
    canonical_payload_path: str = (
        "v9/sc_spaceflight/payloads/rrrm1_blood/"
        "OSD-918_blood_rrrm1_bench.h5ad"
    ),
) -> dict[str, Any]:
    """Audit a staged AnnData payload, or emit a skip report if it is absent."""

    root = Path(repo_root).resolve()
    resolved_manifest_path = _resolve_path(manifest_path, root)
    resolved_requirements_path = _resolve_path(requirements_path, root)
    payload_path = _resolve_path(canonical_payload_path, root)
    manifest = json.loads(resolved_manifest_path.read_text())
    requirement_rows = _requirement_rows(resolved_requirements_path)

    sha256 = "NA"
    size_bytes = "0"
    n_obs = "NA"
    n_vars = "NA"
    if not payload_path.exists():
        payload_status = "missing"
        h5ad_read_status = "not_attempted_no_local_payload"
        results = _missing_payload_results(
            audit_id=audit_id,
            requirement_rows=requirement_rows,
        )
        decision_status = "obs_var_audit_skipped_no_local_payload"
    else:
        payload_status = "local_exists"
        sha256 = _sha256_file(payload_path)
        size_bytes = str(payload_path.stat().st_size)
        adata, h5ad_read_status = _load_anndata(payload_path)
        if adata is None:
            results = _missing_payload_results(
                audit_id=audit_id,
                requirement_rows=requirement_rows,
            )
            decision_status = h5ad_read_status
        else:
            n_obs = str(adata.n_obs)
            n_vars = str(adata.n_vars)
            results = _audit_payload_results(
                audit_id=audit_id,
                requirement_rows=requirement_rows,
                adata=adata,
            )
            decision_status = (
                "obs_var_audit_passed"
                if all(row["audit_status"] == "pass" for row in results if _is_blocker(row))
                else "obs_var_audit_has_blockers"
            )

    pass_count = sum(1 for row in results if row["audit_status"] == "pass")
    fail_count = sum(1 for row in results if row["audit_status"] == "fail")
    skip_count = sum(1 for row in results if row["audit_status"].startswith("skipped"))
    blocker_count = sum(
        1
        for row in results
        if _is_blocker(row) and row["audit_status"] != "pass"
    )
    claim_boundary = "obs_var_audit_skip_only_no_payload_or_score_claim"
    summary = {
        "audit_id": audit_id,
        "decision_status": decision_status,
        "task_id": manifest["task_id"],
        "canonical_payload_path": canonical_payload_path,
        "payload_path_status": payload_status,
        "payload_sha256": sha256,
        "payload_size_bytes": size_bytes,
        "n_obs": n_obs,
        "n_vars": n_vars,
        "requirement_count": str(len(requirement_rows)),
        "pass_count": str(pass_count),
        "fail_count": str(fail_count),
        "skip_count": str(skip_count),
        "blocker_count": str(blocker_count),
        "next_required_block": "V9-SC-006: canonical payload staging or RRRM-1 h5ad regeneration",
        "claim_boundary": claim_boundary,
    }
    payload_manifest = [
        {
            "audit_id": audit_id,
            "payload_id": "rrrm1_blood_canonical_h5ad",
            "payload_path": canonical_payload_path,
            "path_status": payload_status,
            "sha256": sha256,
            "size_bytes": size_bytes,
            "n_obs": n_obs,
            "n_vars": n_vars,
            "h5ad_read_status": h5ad_read_status,
        }
    ]
    review_md = f"""# V9-SC-005 AnnData Obs/Var Audit

Status: `{decision_status}`

Task id: `{manifest["task_id"]}`

Claim boundary: `{claim_boundary}`

## Audit Result

- Canonical payload path: `{canonical_payload_path}`
- Payload path status: `{payload_status}`
- Requirements audited or skipped: {len(requirement_rows)}
- Pass rows: {pass_count}
- Fail rows: {fail_count}
- Skip rows: {skip_count}
- Blocker rows: {blocker_count}

## Interpretation

The canonical RRRM-1 blood AnnData payload is still absent, so this block emits
a machine-readable skip audit instead of failing or fabricating a pass. The
audit implementation is now in place; once the h5ad is staged at the canonical
path, the same API/CLI can hash it, read it, and evaluate the V9-SC-004
obs/var/uns/matrix requirements.

## Not Claimed

- No local h5ad payload is claimed when `payload_path_status` is `missing`.
- No payload checksum is claimed unless the file exists locally.
- No obs/var pass is claimed while all rows are skipped.
- No evaluator, benchmark result, or legacy RRRM score promotion is claimed.

## Next Block

Run `V9-SC-006: canonical payload staging or RRRM-1 h5ad regeneration`.
"""
    return {
        "summary": [summary],
        "results": results,
        "payload_manifest": payload_manifest,
        "review_md": review_md,
    }


def write_sc_obs_var_audit(
    *,
    output_dir: str | Path = "v9/sc_spaceflight",
    repo_root: str | Path = ".",
    audit_id: str = DEFAULT_SC_OBS_VAR_AUDIT_ID,
    manifest_path: str | Path = (
        "v9/sc_spaceflight/task_manifests/"
        "draft_rrrm1_blood_single_cell_spaceflight.json"
    ),
    requirements_path: str | Path = "v9/sc_spaceflight/obs_var_audit_requirements.csv",
    canonical_payload_path: str = (
        "v9/sc_spaceflight/payloads/rrrm1_blood/"
        "OSD-918_blood_rrrm1_bench.h5ad"
    ),
) -> dict[str, Path]:
    """Write the RRRM-1 blood obs/var audit outputs."""

    root = Path(repo_root).resolve()
    package = build_sc_obs_var_audit(
        repo_root=root,
        audit_id=audit_id,
        manifest_path=manifest_path,
        requirements_path=requirements_path,
        canonical_payload_path=canonical_payload_path,
    )
    outdir = _resolve_path(output_dir, root)
    outputs = {
        "summary_csv": outdir / "obs_var_audit_summary.csv",
        "summary_json": outdir / "obs_var_audit_summary.json",
        "results_csv": outdir / "obs_var_audit_results.csv",
        "results_json": outdir / "obs_var_audit_results.json",
        "payload_manifest_csv": outdir / "payload_manifest.csv",
        "payload_manifest_json": outdir / "payload_manifest.json",
        "review_md": root / "docs" / "V9_SC_OBS_VAR_AUDIT.md",
    }
    _write_csv(
        outputs["summary_csv"],
        package["summary"],
        SC_OBS_VAR_AUDIT_SUMMARY_FIELDS,
    )
    _write_json(outputs["summary_json"], package["summary"])
    _write_csv(
        outputs["results_csv"],
        package["results"],
        SC_OBS_VAR_AUDIT_RESULT_FIELDS,
    )
    _write_json(outputs["results_json"], package["results"])
    _write_csv(
        outputs["payload_manifest_csv"],
        package["payload_manifest"],
        SC_OBS_VAR_PAYLOAD_MANIFEST_FIELDS,
    )
    _write_json(outputs["payload_manifest_json"], package["payload_manifest"])
    outputs["review_md"].parent.mkdir(parents=True, exist_ok=True)
    outputs["review_md"].write_text(str(package["review_md"]))
    return outputs
