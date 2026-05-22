---
pretty_name: SpaceBio-Bench v9 Public Bulk Draft
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
- draft
---

# SpaceBio-Bench v9 Public Bulk Draft

This is a draft Hugging Face-style dataset card for the public bulk RNA-seq
portion of SpaceBio-Bench v9.

Release status: draft, not frozen.

This card describes the current local benchmark scaffold and should not yet be
used as final release language. The package has task manifests, source
inventory, OSDR checksum-manifest evidence, baseline outputs, and a draft
Frictionless Data Package descriptor. It does not yet have payload-level hash verification
for every distributed fold matrix.

## Dataset Summary

SpaceBio-Bench v9 is a provenance-first benchmark scaffold for testing
biological AI under spaceflight mission shift.

The current public bulk draft contains leave-one-mission-out classification
tasks built from public mouse bulk RNA-seq studies in NASA's Open Science Data
Repository (OSDR). Each task holds out one mission at a time and evaluates
whether a method can generalize across spaceflight or analog mission contexts
within a tissue-specific expression dataset.

Current scope:

- 8 generated public bulk LOMO task manifests.
- 6 tissue contexts: liver, gastrocnemius, kidney, thymus, skin, and eye.
- 22 deduplicated public OSDR source rows.
- 33 fold definitions.
- 24 baseline runs across 3 simple baseline families.
- 11 resources in `v9/datapackage.draft.json`.

Out of scope for this draft:

- gated human sequence data
- clinical or crew-health recommendations
- intervention or countermeasure claims
- Mars-regime point predictions
- foundation-model ranking claims
- frozen release or DOI claims

## Intended Uses

Intended uses:

- Benchmark simple and advanced methods on mission-held-out public bulk RNA-seq
  tasks.
- Test whether classifiers trained on some missions generalize to held-out
  missions.
- Compare baselines with explicit provenance, task manifests, and run manifests.
- Support method-development studies on spaceflight biological domain shift.
- Serve as the public bulk layer for future SpaceBio-Bench task families.

Responsible use:

- Treat scores as benchmark evidence, not biological mechanism proof.
- Inspect task manifests and source records before comparing methods.
- Report per-task results, not only pooled averages.
- Keep public bulk tasks separate from any future gated human track.
- Cite the upstream OSDR datasets used by a specific analysis.

Out-of-scope uses:

- inferring astronaut health risk
- recommending countermeasures
- claiming operational readiness for space missions
- treating analog missions as interchangeable with deep-space exposure
- ranking virtual-cell or foundation models without adapter-specific validation

## Dataset Structure

The current draft is organized around four artifact classes.

### Metadata Spine

Small Git-friendly files that define the benchmark:

- `v9/task_manifests/*.json`
- `v9/task_manifest_index.csv`
- `v9/task_data_index.csv`
- `v9/source_inventory.csv`
- `v9/source_checksum_audit.csv`
- `v9/datapackage.draft.json`

### Public Bulk Payload Bundle

Large fold-level benchmark inputs referenced by `v9/task_data_index.csv`:

- `train_X.csv`
- `test_X.csv`
- `train_y.csv`
- `test_y.csv`
- `train_meta.csv`
- `test_meta.csv`
- `selected_genes.txt`
- `fold_info.json`

These files are currently indexed but not payload-hash frozen.

### Benchmark Output Bundle

Generated baseline outputs:

- `v9/reports/bulk_lomo_baseline_summary.csv`
- `v9/reports/nearest_centroid/bulk_lomo_summary.csv`
- `v9/reports/sklearn_baselines/bulk_lomo_summary.csv`
- per-task `predictions.csv`
- per-task `metrics.json`
- per-task `run_manifest.json`

### Deferred Or Excluded Artifacts

Not part of the current public bulk draft:

- raw OSDR sequencing payload mirrors
- gated or controlled human data
- local caches or virtual environments
- user `submissions/`
- v8 intervention outputs
- single-cell and radiation/stressor flagship tasks

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

All current source rows are public OSDR accessions.

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

## Provenance And Integrity

Current provenance evidence:

- `v9/source_inventory.csv` has 22 public source rows.
- `v9/source_checksum_audit.csv` has 22 `api_status=ok` rows.
- All 22 source rows have `audit_status=checksum_manifest_parsed`.
- The audit found 39 checksum manifest-like files.
- The audit parsed 8,439 MD5 checksum entries.
- 8,275 parsed entries matched OSDR file-list payload names by exact,
  basename, or suffix matching.

Integrity status:

- `freeze_ready=false` for all current source rows.
- Checksum-manifest evidence exists.
- Payload files have not yet been downloaded and locally hashed as a distributed
  v9 package.
- The current descriptor is `v9/datapackage.draft.json`, not a release
  `datapackage.json`.

Draft Data Package status:

- `spacebio_bench:release_status = draft_not_frozen`
- `spacebio_bench:payload_verification_status = checksum_manifests_parsed_payloads_not_hashed`

## Baseline Results

The current scaffold includes three simple baseline families across all eight
tasks. These results are included to validate the benchmark workflow, not to
claim biological or model superiority.

Mean metrics across eight task rows:

| Baseline | Macro-F1 | Balanced accuracy | AUROC | Calibration error | Mission discrimination |
|---|---:|---:|---:|---:|---:|
| Logistic regression L2 | 0.5532 | 0.5808 | 0.6870 | 0.3439 | NA |
| Nearest centroid | 0.5383 | 0.5685 | 0.6321 | 0.1132 | 0.8733 |
| PCA logistic regression | 0.5353 | 0.5619 | 0.6447 | 0.3747 | 0.9190 |

Interpretation cautions:

- These are scaffold baselines.
- They are not tuned model-comparison endpoints.
- Per-task variability matters.
- Mission-discrimination scores depend on available embedding columns.
- Logistic regression L2 does not currently emit embeddings, so mission
  discrimination is not computed for that baseline.

## Data Fields

Key tabular resources:

- `task_manifest_index.csv`: task ids, tissues, variants, missions, source
  counts, fold counts, metric ids.
- `task_data_index.csv`: fold ids, held-out mission labels, train/test row
  counts, selected-gene counts, fold paths.
- `source_inventory.csv`: OSDR accessions, GLDS prefixes, mission labels,
  tissues, privacy/access status, release target.
- `source_checksum_audit.csv`: OSDR API status, file-list response hashes,
  checksum-manifest files, parsed entry counts, payload-name match counts,
  freeze status.
- `bulk_lomo_baseline_summary.csv`: baseline ids, task ids, validation status,
  prediction counts, metrics, and output paths.

## Privacy And Access

The current public bulk draft uses only source rows marked public in
`v9/source_inventory.csv`.

It does not include controlled-access human sequence data. OSDR documentation
notes that astronaut-derived studies can have public processed data while
sequence data may be controlled access through OSDR. This draft keeps the public
bulk benchmark independent of any such controlled data.

## Biases, Risks, And Limitations

Dataset limitations:

- Mouse tissues and missions are not a complete representation of space biology.
- Mission labels can conflate hardware, time, vehicle, protocol, tissue, and
  processing effects.
- Some tasks include analog or special mission labels such as MHU-2 or OSD-397.
- Bulk RNA-seq tasks do not evaluate cell-type-specific effects.
- Legacy processed fold matrices may reflect earlier preprocessing choices.
- The current package descriptor is draft-only.

Benchmark risks:

- A model can perform well by exploiting mission-correlated technical structure.
- Pooled averages can hide failure on a tissue or mission.
- Stronger methods should be compared only after payload verification and task
  documentation are finalized.
- Countermeasure, intervention, and crew-health language is not supported by
  this draft.

## Licensing And Citation

License status:

- The draft card uses `license: other` because final benchmark artifact
  licensing and OSDR source reuse language need release review.
- OSDR source datasets must be cited according to their OSDR study pages.

Recommended acknowledgment language for this draft:

> Data are courtesy of the NASA Open Science Data Repository.

OSDR citation guidance:

- Cite the NASA OSDR resource and the relevant individual OSDR datasets.
- Use the citation button on each OSDR study page for dataset-specific BibTeX or
  RIS where available.
- Include all OSDR-provided datasets that a downstream analysis uses.

OSDR resource citation from NASA FAQ:

Gebre S G, Scott R T, Saravia-Butler A M, Lopez D K, Sanders L M, and Costes S
V. 2024. NASA Open Science Data Repository: Open Science for Life in Space.
Nucleic Acids Research 53(D1): D1697-D1710.
https://doi.org/10.1093/nar/gkae1116

## Maintenance Notes

Before public release:

- replace draft package status with release status
- add payload-level SHA-256 manifest for distributed fold files
- verify distributed payload hashes
- finalize license field
- add dataset-specific OSDR citations
- decide whether to publish through Hugging Face, Zenodo, or both
- add RO-Crate export or release manifest

## References

- Hugging Face Dataset Cards:
  https://huggingface.co/docs/hub/datasets-cards
- NASA OSDR FAQ:
  https://science.nasa.gov/reference/osdr-faq/
- NASA OSDR Terms and Conditions:
  https://science.nasa.gov/reference/osdr-help-terms-and-conditions/
- Frictionless Data Package:
  https://specs.frictionlessdata.io/data-package/
- v9 public package design:
  `docs/V9_PUBLIC_BULK_PACKAGE_DESIGN.md`
- v9 draft Data Package descriptor:
  `v9/datapackage.draft.json`
