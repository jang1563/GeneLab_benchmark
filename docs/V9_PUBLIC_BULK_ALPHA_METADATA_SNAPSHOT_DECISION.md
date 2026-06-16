# SpaceBio-Bench v9 Public Bulk Metadata Catalog Note

Status: `public_ready`

## Purpose

This note records the public v9 bulk metadata catalog scope. The catalog helps
readers inspect task manifests, source rows, fold indexes, checksum-audit
summaries, reference baselines, and package metadata from the GitHub
repository.

## Catalog Scope

| Area | Current public record |
|---|---|
| Task manifests | 8 public bulk LOMO task manifests |
| Fold definitions | 33 fold rows |
| Source records | 22 public NASA OSDR source rows |
| Checksum audit | 22 source rows with parsed OSDR API/checksum-manifest records |
| Baseline outputs | 24 reference baseline rows across 3 baseline families |
| Data package descriptor | `v9/datapackage.draft.json` |
| Reader-facing card | `docs/v9_hf_dataset_card.md` |

## Public Description

Use this description for the v9 public bulk surface:

> SpaceBio-Bench v9 public bulk is a metadata catalog for public mouse bulk
> RNA-seq mission-held-out tasks. It records task definitions, OSDR source
> coverage, fold indexes, checksum-audit summaries, and reference baseline
> outputs.

## Source And Dataset Notes

- NASA OSDR remains the upstream source for biological data.
- The GitHub catalog records metadata, audit summaries, and baseline outputs.
- The Hugging Face dataset card remains the entry point for processed public
  fold downloads in the v7.1 public package.
- Larger v9 payload bundles can be handled as separate release work when the
  package metadata and verification records are ready.

## Related Files

- `docs/v9_hf_dataset_card.md`
- `v9/README.md`
- `v9/task_manifest_index.csv`
- `v9/task_data_index.csv`
- `v9/source_inventory.csv`
- `v9/source_checksum_audit.csv`
- `v9/reports/bulk_lomo_baseline_summary.csv`
- `v9/datapackage.draft.json`
