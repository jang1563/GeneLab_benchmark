"""Draft source inventories for v9 extension task families."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Mapping, Sequence

from .sources import OSDR_STUDY_URL, SOURCE_INVENTORY_FIELDS


EXTENSION_SOURCE_INVENTORY_FIELDS = [
    *SOURCE_INVENTORY_FIELDS,
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


def _osdr_url(source_id: str) -> str:
    return f"{OSDR_STUDY_URL}{source_id}"


HUMAN_ORGANOID_DRAFT_SOURCES = [
    {
        "source_id": "OSD-863",
        "glds_prefix": "GLDS-716",
        "osd_url": _osdr_url("OSD-863"),
        "url_or_accession": "OSD-863",
        "mission": "SpaceX-19",
        "tissue": "",
        "organism": "Homo sapiens",
        "taxon_id": "9606",
        "species_common_name": "human",
        "material_type": "cells_cultured",
        "model_system": "human_iPSC_neural_organoid",
        "biospecimen_type": "cortical_neural_organoid",
        "assay_modality": "bulk_rna_seq",
        "platform": "Illumina_RNA_seq",
        "spaceflight_environment": "LEO_ISS",
        "ground_control_type": "matched_ground_control",
        "donor_or_strain_block": "donor_block_required",
        "orthology_strategy": "not_applicable_within_human_task",
        "feature_namespace": "human_gene",
        "task_ids": "draft_human_organoid_spaceflight",
        "task_families": "human_organoid_spaceflight",
        "variants": "cortical",
        "access_status": "public",
        "privacy_class": "public",
        "checksum_status": "checksum_manifests_parsed_payloads_not_hashed",
        "release_target": "v9_extension_draft_not_frozen",
        "notes": "Cortical human iPSC-derived neural organoid RNA-seq; keep separate from mouse bulk LOMO.",
        "external_accessions": "GSE259421",
        "publication_urls": "https://academic.oup.com/stcltm/article/13/12/1186/7833382",
        "organoid_type": "cortical_neural_organoid",
        "differentiation_protocol": "published_iPSC_neural_organoid_protocol",
        "microglia_condition": "with_and_without_iPSC_derived_microglia",
        "disease_context": "none_reported_for_primary_RNA_seq_pilot",
        "culture_hardware": "ISS_culture_hardware_pending_manifest_detail",
        "processed_file_status": "OSDR_raw_and_processed_checksum_manifests_parsed",
        "integration_status": "draft_source_inventory_only",
    },
    {
        "source_id": "OSD-871",
        "glds_prefix": "GLDS-720",
        "osd_url": _osdr_url("OSD-871"),
        "url_or_accession": "OSD-871",
        "mission": "SpaceX-19",
        "tissue": "",
        "organism": "Homo sapiens",
        "taxon_id": "9606",
        "species_common_name": "human",
        "material_type": "cells_cultured",
        "model_system": "human_iPSC_neural_organoid",
        "biospecimen_type": "dopaminergic_neural_organoid",
        "assay_modality": "bulk_rna_seq",
        "platform": "Illumina_RNA_seq",
        "spaceflight_environment": "LEO_ISS",
        "ground_control_type": "matched_ground_control",
        "donor_or_strain_block": "donor_block_required",
        "orthology_strategy": "not_applicable_within_human_task",
        "feature_namespace": "human_gene",
        "task_ids": "draft_human_organoid_spaceflight",
        "task_families": "human_organoid_spaceflight",
        "variants": "dopaminergic",
        "access_status": "public",
        "privacy_class": "public",
        "checksum_status": "checksum_manifests_parsed_payloads_not_hashed",
        "release_target": "v9_extension_draft_not_frozen",
        "notes": "Dopaminergic human iPSC-derived neural organoid RNA-seq; pair with OSD-863 only through explicit blocking.",
        "external_accessions": "GSE259421",
        "publication_urls": "https://academic.oup.com/stcltm/article/13/12/1186/7833382",
        "organoid_type": "dopaminergic_neural_organoid",
        "differentiation_protocol": "published_iPSC_neural_organoid_protocol",
        "microglia_condition": "with_and_without_iPSC_derived_microglia",
        "disease_context": "none_reported_for_primary_RNA_seq_pilot",
        "culture_hardware": "ISS_culture_hardware_pending_manifest_detail",
        "processed_file_status": "OSDR_raw_and_processed_checksum_manifests_parsed",
        "integration_status": "draft_source_inventory_only",
    },
]


MULTISPECIES_DRAFT_SOURCES = [
    {
        "source_id": "OSD-207",
        "glds_prefix": "GLDS-207",
        "osd_url": _osdr_url("OSD-207"),
        "url_or_accession": "OSD-207",
        "mission": "spaceflight_vs_ground_control",
        "tissue": "",
        "organism": "Drosophila melanogaster",
        "taxon_id": "7227",
        "species_common_name": "fruit_fly",
        "material_type": "whole_organism_or_tissue",
        "model_system": "drosophila_spaceflight",
        "biospecimen_type": "female_whole_body",
        "assay_modality": "bulk_rna_seq",
        "platform": "GLbulkRNAseq_processed_counts",
        "spaceflight_environment": "spaceflight",
        "ground_control_type": "ground_control",
        "donor_or_strain_block": "genotype_block_required",
        "orthology_strategy": "pending_ortholog_or_pathway_namespace",
        "feature_namespace": "species_gene_ids_pending",
        "task_ids": "draft_multispecies_spaceflight",
        "task_families": "multispecies_spaceflight",
        "variants": "drosophila_whole_body",
        "access_status": "public",
        "privacy_class": "public",
        "checksum_status": "draft_pending_source_audit",
        "release_target": "v9_extension_draft_not_frozen",
        "notes": "Existing local pilot files under data/multispecies/drosophila; do not mix with mouse raw-gene leaderboard.",
        "processed_file_status": "local_processed_counts_and_sample_table_present",
        "integration_status": "draft_source_inventory_only",
    },
    {
        "source_id": "OSD-37",
        "glds_prefix": "GLDS-37",
        "osd_url": _osdr_url("OSD-37"),
        "url_or_accession": "OSD-37",
        "mission": "spaceflight_vs_ground_control",
        "tissue": "",
        "organism": "Arabidopsis thaliana",
        "taxon_id": "3702",
        "species_common_name": "thale_cress",
        "material_type": "plant_tissue",
        "model_system": "arabidopsis_spaceflight",
        "biospecimen_type": "seedling_pool",
        "assay_modality": "bulk_rna_seq",
        "platform": "GLbulkRNAseq_processed_counts",
        "spaceflight_environment": "spaceflight",
        "ground_control_type": "ground_control",
        "donor_or_strain_block": "ecotype_block_required",
        "orthology_strategy": "pending_ortholog_or_pathway_namespace",
        "feature_namespace": "species_gene_ids_pending",
        "task_ids": "draft_multispecies_spaceflight",
        "task_families": "multispecies_spaceflight",
        "variants": "arabidopsis_seedling_pool",
        "access_status": "public",
        "privacy_class": "public",
        "checksum_status": "draft_pending_source_audit",
        "release_target": "v9_extension_draft_not_frozen",
        "notes": "Existing local pilot files under data/multispecies/arabidopsis.",
        "processed_file_status": "local_processed_counts_and_sample_table_present",
        "integration_status": "draft_source_inventory_only",
    },
    {
        "source_id": "OSD-120",
        "glds_prefix": "GLDS-120",
        "osd_url": _osdr_url("OSD-120"),
        "url_or_accession": "OSD-120",
        "mission": "spaceflight_vs_ground_control_by_light_condition",
        "tissue": "",
        "organism": "Arabidopsis thaliana",
        "taxon_id": "3702",
        "species_common_name": "thale_cress",
        "material_type": "plant_tissue",
        "model_system": "arabidopsis_spaceflight",
        "biospecimen_type": "root",
        "assay_modality": "bulk_rna_seq",
        "platform": "GLbulkRNAseq_processed_counts",
        "spaceflight_environment": "spaceflight",
        "ground_control_type": "ground_control",
        "donor_or_strain_block": "genotype_and_light_block_required",
        "orthology_strategy": "pending_ortholog_or_pathway_namespace",
        "feature_namespace": "species_gene_ids_pending",
        "task_ids": "draft_multispecies_spaceflight",
        "task_families": "multispecies_spaceflight",
        "variants": "arabidopsis_root_light_interaction",
        "access_status": "public",
        "privacy_class": "public",
        "checksum_status": "draft_pending_source_audit",
        "release_target": "v9_extension_draft_not_frozen",
        "notes": "Existing local pilot files under data/multispecies/arabidopsis; split must handle light condition.",
        "processed_file_status": "local_processed_counts_and_sample_table_present",
        "integration_status": "draft_source_inventory_only",
    },
]


def normalize_extension_rows(
    rows: Sequence[Mapping[str, str]],
) -> list[dict[str, str]]:
    """Normalize extension source rows to the draft inventory field order."""

    if not rows:
        raise ValueError("cannot normalize an empty extension source inventory")
    return [
        {field: str(row.get(field, "") or "") for field in EXTENSION_SOURCE_INVENTORY_FIELDS}
        for row in rows
    ]


def write_extension_source_inventory(
    rows: Sequence[Mapping[str, str]],
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> tuple[Path, Path]:
    """Write a draft extension source inventory as CSV and JSON."""

    normalized = normalize_extension_rows(rows)
    output_csv = Path(csv_path)
    output_json = Path(json_path)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    output_json.parent.mkdir(parents=True, exist_ok=True)

    with output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=EXTENSION_SOURCE_INVENTORY_FIELDS)
        writer.writeheader()
        writer.writerows(normalized)
    output_json.write_text(json.dumps(normalized, indent=2, sort_keys=True) + "\n")
    return output_csv, output_json


def read_extension_source_inventory(path: str | Path) -> list[dict[str, str]]:
    """Read a draft extension source inventory from JSON or CSV."""

    inventory_path = Path(path)
    if inventory_path.suffix.lower() == ".json":
        payload = json.loads(inventory_path.read_text())
        if not isinstance(payload, list):
            raise ValueError(f"{inventory_path} must contain a list of source rows")
        return normalize_extension_rows(payload)
    with inventory_path.open(newline="") as handle:
        return normalize_extension_rows([dict(row) for row in csv.DictReader(handle)])
