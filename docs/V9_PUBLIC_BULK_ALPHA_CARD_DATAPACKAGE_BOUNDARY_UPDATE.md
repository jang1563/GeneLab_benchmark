# SpaceBio-Bench v9 Dataset Card And Data Package Note

Status: `public_ready`

## Purpose

This note records how the v9 public bulk metadata catalog is represented in the
reader-facing dataset card and the Data Package descriptor.

## Reader-Facing Card

The v9 card at `docs/v9_hf_dataset_card.md` presents the public bulk surface as
a metadata catalog. It summarizes:

- 8 public bulk leave-one-mission-out task manifests;
- 6 tissue contexts;
- 22 public NASA OSDR source rows;
- 33 fold definitions;
- 24 reference baseline runs across 3 baseline families;
- 21 catalog, audit, and output resources in `v9/datapackage.draft.json`.

The card is written for reviewers and method developers who want to inspect the
task catalog, source coverage, fold index, and reference outputs.

## Data Package Descriptor

The descriptor at `v9/datapackage.draft.json` lists the metadata and output
resources that make up the public bulk catalog. It includes task indexes,
source inventory files, checksum-audit summaries, baseline summaries, and
supporting report tables.

The descriptor is useful for machine reading and resource inventory checks. The
human-facing explanation remains `docs/v9_hf_dataset_card.md`.

## Related Files

- `docs/v9_hf_dataset_card.md`
- `v9/datapackage.draft.json`
- `v9/README.md`
- `v9/reports/README.md`
- `release/release_manifest.json`
