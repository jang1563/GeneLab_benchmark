---
pretty_name: SpaceBio-Bench v9 Public Bulk Metadata Catalog
license: other
task_categories:
- tabular-classification
tags:
- biology
- space-biology
- genomics
- benchmark
- tabular
- mission-held-out
- OSDR
- metadata-catalog
---

# SpaceBio-Bench v9 Public Bulk Metadata Catalog

This Hugging Face-style card documents the public bulk RNA-seq metadata catalog
for SpaceBio-Bench v9. It is designed for reviewers and method developers who
want to inspect task definitions, OSDR source coverage, fold indexes, and
baseline outputs directly from the GitHub repository.

The catalog focuses on metadata and reference outputs. Fold-level expression
matrices are indexed by the catalog when available, while packaged public fold
downloads remain documented in the v7.1 Hugging Face dataset card.

## Related Documents

- Main repository: <https://github.com/jang1563/GeneLab_benchmark>
- v7.1 public fold package card: `docs/hf_dataset_card.md`
- Release manifest: `release/release_manifest.json`
- Transparency card pack: `docs/SPACEBIOBENCH_TRANSPARENCY_CARD_PACK.md`

## Dataset Summary

SpaceBio-Bench v9 extends the mission-held-out benchmark design to a public
bulk RNA-seq task catalog. Each task is tissue-specific and asks whether a
method trained on one set of missions can generalize to a held-out mission.

Current catalog scope:

- 8 public bulk leave-one-mission-out task manifests.
- 6 tissue contexts: liver, gastrocnemius, kidney, thymus, skin, and eye.
- 22 public NASA OSDR source rows.
- 33 fold definitions.
- 24 reference baseline runs across 3 simple baseline families.
- 21 catalog, audit, and output resources in `v9/datapackage.draft.json`.

Scope:

- astronaut health-risk inference
- clinical, intervention, or countermeasure recommendations
- Mars-regime point prediction
- foundation-model leaderboard claims without adapter-specific validation
- controlled-access human sequence analysis

## Intended Uses

Use this catalog to:

- inspect public bulk mission-held-out task definitions;
- compare simple baseline outputs as reproducible reference rows;
- prepare method-development experiments for spaceflight biological domain
  shift;
- audit which public OSDR accessions support each task;
- keep per-task, per-mission evaluation visible instead of relying only on
  pooled averages.

Recommended reporting:

- Treat scores as benchmark evidence for the stated task definition.
- Report per-task results alongside pooled summaries.
- Compare stronger methods against the fixed task manifests and fold IDs.
- Cite the upstream OSDR datasets used by a specific analysis.

## Dataset Structure

### Core Catalog Files

Small Git-friendly files define the public bulk task catalog:

- `v9/task_manifests/*.json`
- `v9/task_manifest_index.csv`
- `v9/task_data_index.csv`
- `v9/source_inventory.csv`
- `v9/source_checksum_audit.csv`
- `v9/datapackage.draft.json`

### Indexed Fold Files

The catalog records fold-level files with the following contract:

- `train_X.csv` and `test_X.csv`: sample-by-gene expression matrices
- `train_y.csv` and `test_y.csv`: binary Flight/Ground labels
- `train_meta.csv` and `test_meta.csv`: sample metadata
- `selected_genes.txt`: fold-specific selected genes
- `fold_info.json`: held-out mission, training missions, and sample counts

### Baseline Output Files

Reference baseline outputs are stored under `v9/reports/`:

- `v9/reports/bulk_lomo_baseline_summary.csv`
- `v9/reports/nearest_centroid/bulk_lomo_summary.csv`
- `v9/reports/sklearn_baselines/bulk_lomo_summary.csv`
- per-task `predictions.csv`
- per-task `metrics.json`
- per-task `run_manifest.json`

## Task Table

| Task ID | Tissue | Variant | Missions | Folds | Sources | Mission labels |
|---|---|---:|---:|---:|---:|---|
| `A1_liver_bulk_lomo` | liver | canonical | 6 | 6 | 6 | MHU-2; RR-1; RR-3; RR-6; RR-8; RR-9 |
| `A1_liver_bulk_lomo_combat` | liver | combat | 6 | 6 | 6 | MHU-2; RR-1; RR-3; RR-6; RR-8; RR-9 |
| `A1_liver_bulk_lomo_iss_only` | liver | iss_only | 5 | 5 | 5 | RR-1; RR-3; RR-6; RR-8; RR-9 |
| `A2_gastrocnemius_bulk_lomo` | gastrocnemius | canonical | 3 | 3 | 3 | RR-1; RR-5; RR-9 |
| `A3_kidney_bulk_lomo` | kidney | canonical | 3 | 3 | 3 | RR-1; RR-3; RR-7 |
| `A4_thymus_bulk_lomo` | thymus | canonical | 4 | 4 | 3 | MHU-1; MHU-2; RR-6; RR-9 |
| `A5_skin_bulk_lomo` | skin | canonical | 3 | 3 | 4 | MHU-2; RR-6; RR-7 |
| `A6_eye_bulk_lomo` | eye | canonical | 3 | 3 | 3 | OSD-397; RR-1; RR-3 |

## Source Table

All source rows in this catalog are public OSDR accessions.

| OSDR source | GLDS prefix | Tissue | Mission | Task IDs |
|---|---|---|---|---|
| OSD-100 | GLDS-100 | eye | RR-1 | `A6_eye_bulk_lomo` |
| OSD-101 | GLDS-101 | gastrocnemius | RR-1 | `A2_gastrocnemius_bulk_lomo` |
| OSD-102 | GLDS-102 | kidney | RR-1 | `A3_kidney_bulk_lomo` |
| OSD-137 | GLDS-137 | liver | RR-3 | `A1_liver_bulk_lomo`; `A1_liver_bulk_lomo_combat`; `A1_liver_bulk_lomo_iss_only` |
| OSD-163 | GLDS-163 | kidney | RR-3 | `A3_kidney_bulk_lomo` |
| OSD-194 | GLDS-194 | eye | RR-3 | `A6_eye_bulk_lomo` |
| OSD-238 | GLDS-238 | skin | MHU-2 (dorsal) | `A5_skin_bulk_lomo` |
| OSD-239 | GLDS-239 | skin | MHU-2 (femoral) | `A5_skin_bulk_lomo` |
| OSD-242 | GLDS-242 | liver | RR-9 | `A1_liver_bulk_lomo`; `A1_liver_bulk_lomo_combat`; `A1_liver_bulk_lomo_iss_only` |
| OSD-243 | GLDS-243 | skin | RR-6 | `A5_skin_bulk_lomo` |
| OSD-244 | GLDS-244 | thymus | RR-6 | `A4_thymus_bulk_lomo` |
| OSD-245 | GLDS-245 | liver | RR-6 | `A1_liver_bulk_lomo`; `A1_liver_bulk_lomo_combat`; `A1_liver_bulk_lomo_iss_only` |
| OSD-253 | GLDS-253 | kidney | RR-7 | `A3_kidney_bulk_lomo` |
| OSD-254 | GLDS-254 | skin | RR-7 | `A5_skin_bulk_lomo` |
| OSD-289 | GLDS-289 | thymus | MHU-2 | `A4_thymus_bulk_lomo` |
| OSD-326 | GLDS-326 | gastrocnemius | RR-9 | `A2_gastrocnemius_bulk_lomo` |
| OSD-379 | GLDS-379 | liver | RR-8 | `A1_liver_bulk_lomo`; `A1_liver_bulk_lomo_combat`; `A1_liver_bulk_lomo_iss_only` |
| OSD-397 | GLDS-397 | eye | OSD-397 | `A6_eye_bulk_lomo` |
| OSD-401 | GLDS-401 | gastrocnemius | RR-5 | `A2_gastrocnemius_bulk_lomo` |
| OSD-421 | GLDS-421 | thymus | RR-9 | `A4_thymus_bulk_lomo` |
| OSD-48 | GLDS-48 | liver | RR-1 | `A1_liver_bulk_lomo`; `A1_liver_bulk_lomo_combat`; `A1_liver_bulk_lomo_iss_only` |
| OSD-686 | GLDS-617 | liver | MHU-2 | `A1_liver_bulk_lomo`; `A1_liver_bulk_lomo_combat` |

## Catalog Checks

The public catalog includes a compact audit table for OSDR source coverage:

- `v9/source_inventory.csv` records 22 public source rows.
- `v9/source_checksum_audit.csv` records 22 rows with `api_status=ok`.
- The audit parsed 39 checksum-manifest-like files and 8,439 MD5 entries.
- 8,275 parsed entries matched OSDR file-list payload names by exact,
  basename, or suffix matching.

These checks help users inspect source coverage before running or extending a
benchmark workflow.

## Baseline Results

The current catalog includes three simple baseline families across all eight
task rows. These are reproducible reference anchors for workflow comparison.

Mean metrics across eight task rows:

| Baseline | Macro-F1 | Balanced accuracy | AUROC | Calibration error | Mission discrimination |
|---|---:|---:|---:|---:|---:|
| Logistic regression L2 | 0.5532 | 0.5808 | 0.6870 | 0.3439 | NA |
| Nearest centroid | 0.5383 | 0.5685 | 0.6321 | 0.1132 | 0.8733 |
| PCA logistic regression | 0.5353 | 0.5619 | 0.6447 | 0.3747 | 0.9190 |

Report per-task rows when comparing methods. Pooled means are useful for
orientation, but tissue and mission variability are central to this benchmark.

## Data Fields

Key tabular resources:

- `task_manifest_index.csv`: task IDs, tissues, variants, missions, source
  counts, fold counts, and metric IDs.
- `task_data_index.csv`: fold IDs, held-out mission labels, train/test row
  counts, selected-gene counts, and fold paths.
- `source_inventory.csv`: OSDR accessions, GLDS prefixes, mission labels,
  tissues, access status, and release target.
- `source_checksum_audit.csv`: OSDR API status, file-list response hashes,
  checksum-manifest files, parsed entry counts, and payload-name match counts.
- `bulk_lomo_baseline_summary.csv`: baseline IDs, task IDs, validation status,
  prediction counts, metrics, and output paths.

## Privacy And Access

This catalog uses public OSDR source rows. It does not include controlled-access
human sequence data.

## Limitations

- Mouse bulk RNA-seq tasks are not a complete representation of space biology.
- Mission labels can combine hardware, time, vehicle, protocol, tissue, and
  processing effects.
- Some tasks include analog or special mission labels such as MHU-2 or OSD-397.
- Bulk RNA-seq tasks do not resolve cell-type-specific effects.
- Baseline rows serve as reference comparisons for workflow checks.
- The catalog describes task metadata and benchmark outputs rather than
  countermeasure, intervention, crew-health, or operational-readiness use cases.

## Licensing And Citation

The v9 metadata catalog uses `license: other` because source reuse follows the
terms of the underlying OSDR studies. Code in the GitHub repository is MIT
licensed. Cite OSDR study pages for the datasets used by a downstream analysis.

Recommended acknowledgment:

> Data are courtesy of the NASA Open Science Data Repository.

OSDR resource citation:

Gebre S G, Scott R T, Saravia-Butler A M, Lopez D K, Sanders L M, and Costes S
V. 2024. NASA Open Science Data Repository: Open Science for Life in Space.
Nucleic Acids Research 53(D1): D1697-D1710.
https://doi.org/10.1093/nar/gkae1116

## References

- Hugging Face Dataset Cards:
  https://huggingface.co/docs/hub/datasets-cards
- Frictionless Data Package:
  https://specs.frictionlessdata.io/data-package/
- NASA OSDR Biological Data API:
  https://visualization.osdr.nasa.gov/biodata/api/
- NASA OSDR FAQ:
  https://science.nasa.gov/reference/osdr-faq/
- NASA OSDR Terms and Conditions:
  https://science.nasa.gov/reference/osdr-help-terms-and-conditions/
- v9 Data Package descriptor: `v9/datapackage.draft.json`
