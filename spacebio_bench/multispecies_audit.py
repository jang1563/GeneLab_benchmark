"""Local audits for v9 multispecies extension inputs."""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path
from typing import Mapping, Sequence


MULTISPECIES_SAMPLE_FACTOR_FIELDS = [
    "source_id",
    "glds_prefix",
    "sample_id",
    "parse_status",
    "true_label",
    "condition_raw",
    "condition_stratum",
    "organism",
    "taxon_id",
    "species_common_name",
    "material_type",
    "model_system",
    "biospecimen_type",
    "assay_modality",
    "feature_namespace",
    "genotype_or_ecotype",
    "genotype_or_ecotype_type",
    "parental_strain_or_background",
    "allele_or_variant",
    "light_treatment",
    "replicate_hint",
    "spaceflight_environment",
    "ground_control_type",
    "task_family",
    "task_id",
    "variant",
    "processed_file_status",
    "local_sample_table_path",
    "local_matrix_path",
]

MULTISPECIES_EXPRESSION_MATRIX_AUDIT_FIELDS = [
    "source_id",
    "glds_prefix",
    "local_matrix_path",
    "local_sample_table_path",
    "n_feature_rows",
    "n_sample_columns",
    "n_sample_factor_rows",
    "n_matching_sample_columns",
    "n_missing_sample_factor_rows",
    "n_extra_matrix_columns",
    "matrix_columns_match_sample_factors",
    "missing_sample_factor_rows",
    "extra_matrix_columns",
    "feature_namespace",
    "organism",
    "assay_modality",
    "processed_file_status",
    "audit_status",
]

MULTISPECIES_LOCAL_FILES = {
    "OSD-207": {
        "matrix": "data/multispecies/drosophila/GLDS-207_rna_seq_Normalized_Counts_GLbulkRNAseq.csv",
        "sample_table": "data/multispecies/drosophila/GLDS-207_rna_seq_SampleTable_GLbulkRNAseq.csv",
    },
    "OSD-37": {
        "matrix": "data/multispecies/arabidopsis/GLDS-37_rna_seq_Normalized_Counts_GLbulkRNAseq.csv",
        "sample_table": "data/multispecies/arabidopsis/GLDS-37_rna_seq_SampleTable_GLbulkRNAseq.csv",
    },
    "OSD-120": {
        "matrix": "data/multispecies/arabidopsis/GLDS-120_rna_seq_Normalized_Counts_GLbulkRNAseq.csv",
        "sample_table": "data/multispecies/arabidopsis/GLDS-120_rna_seq_SampleTable_GLbulkRNAseq.csv",
    },
}


def _repo_path(repo_root: str | Path, relative_path: str) -> Path:
    return Path(repo_root) / relative_path


def _source_local_paths(
    source_row: Mapping[str, str],
    *,
    repo_root: str | Path,
) -> tuple[Path, Path]:
    source_id = str(source_row.get("source_id", "") or "")
    paths = MULTISPECIES_LOCAL_FILES.get(source_id)
    if paths is None:
        raise ValueError(f"no local multispecies files registered for {source_id}")
    return (
        _repo_path(repo_root, paths["matrix"]),
        _repo_path(repo_root, paths["sample_table"]),
    )


def _label_from_condition_token(token: str) -> str:
    if token == "Ground.Control":
        return "Ground"
    if token == "Space.Flight":
        return "LEO_or_ISS"
    return ""


def parse_multispecies_condition(source_id: str, condition: str) -> dict[str, str]:
    """Parse a local multispecies SampleTable condition into benchmark factors."""

    parts = [part for part in condition.split("...") if part]
    base = {
        "parse_status": "unparsed",
        "true_label": "",
        "condition_stratum": "",
        "genotype_or_ecotype": "",
        "genotype_or_ecotype_type": "",
        "parental_strain_or_background": "",
        "allele_or_variant": "",
        "light_treatment": "",
    }
    if source_id == "OSD-207" and len(parts) == 3:
        background, variant, label_token = parts
        genotype_or_ecotype = f"{background}_{variant}"
        base.update(
            {
                "parse_status": "parsed",
                "true_label": _label_from_condition_token(label_token),
                "condition_stratum": genotype_or_ecotype,
                "genotype_or_ecotype": genotype_or_ecotype,
                "genotype_or_ecotype_type": "drosophila_background_variant",
                "parental_strain_or_background": background,
                "allele_or_variant": variant,
            }
        )
    elif source_id == "OSD-37" and len(parts) == 2:
        label_token, ecotype = parts
        base.update(
            {
                "parse_status": "parsed",
                "true_label": _label_from_condition_token(label_token),
                "condition_stratum": ecotype,
                "genotype_or_ecotype": ecotype,
                "genotype_or_ecotype_type": "arabidopsis_ecotype",
            }
        )
    elif source_id == "OSD-120" and len(parts) == 3:
        genotype_or_ecotype, label_token, light_treatment = parts
        base.update(
            {
                "parse_status": "parsed",
                "true_label": _label_from_condition_token(label_token),
                "condition_stratum": f"{genotype_or_ecotype}|{light_treatment}",
                "genotype_or_ecotype": genotype_or_ecotype,
                "genotype_or_ecotype_type": "arabidopsis_genotype_or_ecotype",
                "light_treatment": light_treatment,
            }
        )
    if base["parse_status"] == "parsed" and not base["true_label"]:
        base["parse_status"] = "unparsed"
    return base


def _sample_id_field(fieldnames: Sequence[str] | None) -> str:
    for field in fieldnames or []:
        if field != "condition":
            return field
    raise ValueError("sample table must contain a sample-id column")


def _read_sample_table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        sample_field = _sample_id_field(reader.fieldnames)
        rows = []
        for row in reader:
            rows.append(
                {
                    "sample_id": str(row.get(sample_field, "") or ""),
                    "condition": str(row.get("condition", "") or ""),
                }
            )
    return rows


def _replicate_hint(sample_id: str) -> str:
    match = re.search(r"(Rep\d+)", sample_id)
    if match:
        return match.group(1)
    suffix = re.search(r"-(\d+)$", sample_id)
    if suffix:
        return f"sample_{suffix.group(1)}"
    return ""


def build_multispecies_sample_factor_rows(
    source_rows: Sequence[Mapping[str, str]],
    *,
    repo_root: str | Path = ".",
) -> list[dict[str, str]]:
    """Build sample-factor rows from local multispecies SampleTable files."""

    output: list[dict[str, str]] = []
    for source_row in source_rows:
        source_id = str(source_row.get("source_id", "") or "")
        matrix_path, sample_table_path = _source_local_paths(source_row, repo_root=repo_root)
        for sample_row in _read_sample_table(sample_table_path):
            condition = sample_row["condition"]
            parsed = parse_multispecies_condition(source_id, condition)
            output.append(
                {
                    "source_id": source_id,
                    "glds_prefix": str(source_row.get("glds_prefix", "") or ""),
                    "sample_id": sample_row["sample_id"],
                    "parse_status": parsed["parse_status"],
                    "true_label": parsed["true_label"],
                    "condition_raw": condition,
                    "condition_stratum": parsed["condition_stratum"],
                    "organism": str(source_row.get("organism", "") or ""),
                    "taxon_id": str(source_row.get("taxon_id", "") or ""),
                    "species_common_name": str(source_row.get("species_common_name", "") or ""),
                    "material_type": str(source_row.get("material_type", "") or ""),
                    "model_system": str(source_row.get("model_system", "") or ""),
                    "biospecimen_type": str(source_row.get("biospecimen_type", "") or ""),
                    "assay_modality": str(source_row.get("assay_modality", "") or ""),
                    "feature_namespace": str(source_row.get("feature_namespace", "") or ""),
                    "genotype_or_ecotype": parsed["genotype_or_ecotype"],
                    "genotype_or_ecotype_type": parsed["genotype_or_ecotype_type"],
                    "parental_strain_or_background": parsed["parental_strain_or_background"],
                    "allele_or_variant": parsed["allele_or_variant"],
                    "light_treatment": parsed["light_treatment"],
                    "replicate_hint": _replicate_hint(sample_row["sample_id"]),
                    "spaceflight_environment": str(source_row.get("spaceflight_environment", "") or ""),
                    "ground_control_type": str(source_row.get("ground_control_type", "") or ""),
                    "task_family": str(source_row.get("task_families", "") or ""),
                    "task_id": str(source_row.get("task_ids", "") or ""),
                    "variant": str(source_row.get("variants", "") or ""),
                    "processed_file_status": str(source_row.get("processed_file_status", "") or ""),
                    "local_sample_table_path": sample_table_path.as_posix(),
                    "local_matrix_path": matrix_path.as_posix(),
                }
            )
    return output


def _inspect_matrix(path: Path) -> tuple[int, list[str]]:
    with path.open(newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"{path} is empty") from exc
        sample_columns = [column for column in header[1:] if column]
        n_feature_rows = sum(1 for row in reader if row and any(cell != "" for cell in row))
    return n_feature_rows, sample_columns


def audit_multispecies_expression_matrices(
    source_rows: Sequence[Mapping[str, str]],
    *,
    sample_factor_rows: Sequence[Mapping[str, str]],
    repo_root: str | Path = ".",
) -> list[dict[str, str]]:
    """Check local multispecies matrices against parsed sample-factor rows."""

    factors_by_source: dict[str, list[Mapping[str, str]]] = {}
    for row in sample_factor_rows:
        source_id = str(row.get("source_id", "") or "")
        if source_id:
            factors_by_source.setdefault(source_id, []).append(row)

    output: list[dict[str, str]] = []
    for source_row in source_rows:
        source_id = str(source_row.get("source_id", "") or "")
        matrix_path, sample_table_path = _source_local_paths(source_row, repo_root=repo_root)
        source_factors = factors_by_source.get(source_id, [])
        factor_ids = {str(row.get("sample_id", "") or "") for row in source_factors}
        n_feature_rows, matrix_sample_columns = _inspect_matrix(matrix_path)
        matrix_ids = set(matrix_sample_columns)
        missing = sorted(factor_ids - matrix_ids)
        extra = sorted(matrix_ids - factor_ids)
        matching = sorted(factor_ids & matrix_ids)
        aligned = not missing and not extra and len(matrix_sample_columns) == len(source_factors)
        output.append(
            {
                "source_id": source_id,
                "glds_prefix": str(source_row.get("glds_prefix", "") or ""),
                "local_matrix_path": matrix_path.as_posix(),
                "local_sample_table_path": sample_table_path.as_posix(),
                "n_feature_rows": str(n_feature_rows),
                "n_sample_columns": str(len(matrix_sample_columns)),
                "n_sample_factor_rows": str(len(source_factors)),
                "n_matching_sample_columns": str(len(matching)),
                "n_missing_sample_factor_rows": str(len(missing)),
                "n_extra_matrix_columns": str(len(extra)),
                "matrix_columns_match_sample_factors": str(aligned),
                "missing_sample_factor_rows": ";".join(missing),
                "extra_matrix_columns": ";".join(extra),
                "feature_namespace": str(source_row.get("feature_namespace", "") or ""),
                "organism": str(source_row.get("organism", "") or ""),
                "assay_modality": str(source_row.get("assay_modality", "") or ""),
                "processed_file_status": str(source_row.get("processed_file_status", "") or ""),
                "audit_status": (
                    "matrix_local_sample_aligned" if aligned else "matrix_local_sample_mismatch"
                ),
            }
        )
    return output


def _write_rows(
    rows: Sequence[Mapping[str, str]],
    *,
    fields: Sequence[str],
    csv_path: str | Path,
    json_path: str | Path,
) -> tuple[Path, Path]:
    if not rows:
        raise ValueError("cannot write an empty multispecies audit table")
    output_csv = Path(csv_path)
    output_json = Path(json_path)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    normalized = [{field: str(row.get(field, "") or "") for field in fields} for row in rows]
    with output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(normalized)
    output_json.write_text(json.dumps(normalized, indent=2, sort_keys=True) + "\n")
    return output_csv, output_json


def write_multispecies_sample_factor_table(
    rows: Sequence[Mapping[str, str]],
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> tuple[Path, Path]:
    return _write_rows(
        rows,
        fields=MULTISPECIES_SAMPLE_FACTOR_FIELDS,
        csv_path=csv_path,
        json_path=json_path,
    )


def write_multispecies_expression_matrix_audit(
    rows: Sequence[Mapping[str, str]],
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> tuple[Path, Path]:
    return _write_rows(
        rows,
        fields=MULTISPECIES_EXPRESSION_MATRIX_AUDIT_FIELDS,
        csv_path=csv_path,
        json_path=json_path,
    )
