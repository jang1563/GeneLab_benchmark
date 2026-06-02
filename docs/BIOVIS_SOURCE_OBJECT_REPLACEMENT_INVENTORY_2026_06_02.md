# BioVis source-object replacement inventory

Date: 2026-06-02

## Anchor

This note continues the v0.4 source-proof module work. The goal is to decide
which real objects should replace the placeholder proof slots before building
production slides.

This is intentionally not a final slide pass. It is the evidence inventory that
prevents the next slide session from drifting into decoration or unsupported
claims.

## Generated Outputs

Builder:

- `scripts/build_biovis_source_object_replacement_inventory.py`

Output root:

- `output/biovis_source_object_replacement_inventory_v0_1/`

Generated files:

- `output/biovis_source_object_replacement_inventory_v0_1/source_object_replacement_inventory.json`
- `output/biovis_source_object_replacement_inventory_v0_1/source_object_replacement_inventory.csv`
- `output/biovis_source_object_replacement_inventory_v0_1/source_object_replacement_contact_sheet.png`
- `output/biovis_source_object_replacement_inventory_v0_1/source_object_replacement_inventory_qa.json`

## Local Inputs Read

Human organoid:

- `v9/human_organoid/source_inventory.draft.json`
- `v9/human_organoid/task_manifests/draft_human_organoid_spaceflight.json`
- `v9/human_organoid/expression_matrix_audit.draft.json`
- `v9/human_organoid/reports/source_transfer_signature/metrics.json`

Multispecies:

- `v9/multispecies/source_inventory.draft.json`
- `v9/multispecies/expression_matrix_audit.draft.json`
- `v9/multispecies/task_manifests/draft_osd37_arabidopsis_seedling_spaceflight.json`
- `v9/multispecies/task_manifests/draft_osd207_drosophila_whole_body_spaceflight.json`
- `v9/multispecies/interaction_task_manifests/draft_osd120_arabidopsis_root_light_interaction_spaceflight.json`

Single-cell:

- `v9/sc_spaceflight/task_manifests/draft_rrrm1_blood_single_cell_spaceflight.json`
- `v9/sc_spaceflight/obs_var_audit_summary.json`
- `v9/sc_spaceflight/asset_inventory_summary.json`

## External Source Surfaces

Official source-record targets:

- `https://osdr.nasa.gov/bio/repo/data/studies/OSD-863`
- `https://osdr.nasa.gov/bio/repo/data/studies/OSD-871`
- `https://osdr.nasa.gov/bio/repo/data/studies/OSD-37`
- `https://osdr.nasa.gov/bio/repo/data/studies/OSD-207`
- `https://osdr.nasa.gov/bio/repo/data/studies/OSD-120`
- `https://osdr.nasa.gov/bio/repo/data/studies/OSD-918`

High-value external proof details:

- OSD-863 is the strongest first source-record proof because the official OSDR
  page exposes the human cortical organoid record, GLDS-716, DOI, GSE259421, and
  the related OSD-871 study surface.
- OSD-120 is the strongest methods proof because the local manifest exposes a
  genotype-by-light interaction split with concrete sample counts.
- OSD-918 is useful as a source-record proof, but not as a single-cell embedding
  or result proof yet because the canonical local AnnData payload is missing.

## Inventory Verdict

Verdict: inventory ready for source-object replacement sprint.

Main-slide candidates:

- `organoid_osdr_source_records`
- `osd120_interaction_split`

Appendix or reviewer-backup candidates:

- `organoid_expression_matrix_qc`
- `organoid_diagnostic_result_claim`
- `multispecies_source_matrix_inventory`

Blocked or caveated candidate:

- `rrrm1_blood_sc_source_contract`

## Replacement Map

| Proof id | Module | Replacement object | Slide role | Boundary |
| --- | --- | --- | --- | --- |
| `organoid_osdr_source_records` | `source_dataset_record_plate` | OSD-863 and OSD-871 OSDR record screenshots | Main bridge or appendix source inventory | Source proof only |
| `organoid_expression_matrix_qc` | `expression_matrix_proof_plate` | Matrix audit rows and compact matrix/QC crop | Appendix or methods bridge | Matrix aligned, payload freeze still draft |
| `organoid_blocked_split` | `held_out_task_proof_plate` | Organoid type and microglia blocked-fold visual | Methods explainer | Single ISS mission, not mission-held-out |
| `organoid_diagnostic_result_claim` | `result_claim_plate` | Source-transfer signature metrics | Appendix result proof | Diagnostic only, no leaderboard claim |
| `multispecies_source_matrix_inventory` | `source_dataset_record_plate`; `expression_matrix_proof_plate` | OSD-37, OSD-207, OSD-120 source and matrix audit rows | Appendix or extension bridge | Species-native within-source pilots |
| `osd120_interaction_split` | `held_out_task_proof_plate` | Genotype-by-light split manifest visual | Main methods explainer candidate | Within-source interaction diagnostic |
| `rrrm1_blood_sc_source_contract` | `source_dataset_record_plate`; `single_cell_embedding_proof_plate` | OSD-918 record and AnnData blocker panel | Future single-cell lane proof | No local h5ad, no score claim |

## Design Implications

For a premium main slide, do not show all seven proof objects at once.

Recommended main-slide options:

1. Human organoid source proof:
   - Dark Z-stack.
   - One OSDR screenshot for OSD-863 as the dominant proof surface.
   - Smaller companion badge for OSD-871/GSE259421.
   - Minimal editable interpretation text: "human organoid extension is public
     and source-backed, but still draft/non-leaderboard."

2. OSD-120 method proof:
   - Dark Z-stack.
   - One held-out split module showing genotype/ecotype holdout.
   - One secondary lane showing dark/light treatment split.
   - Caveat footer: "within-source interaction diagnostic, not leave-one-mission-out."

3. Appendix evidence board:
   - Multiple proof modules allowed.
   - Include OSDR records, matrix row counts, sha256, fold counts, and caveats.

## Next Capture Sprint

P0 capture sequence:

1. Capture OSD-863 source-record screenshot from OSDR.
2. Capture OSD-871 source-record screenshot or at minimum the related-study
   proof visible from the OSD-863 page.
3. Render local organoid matrix audit proof crop:
   - OSD-863: 30,408 genes x 19 samples.
   - OSD-871: 30,269 genes x 23 samples.
4. Render OSD-120 split proof crop:
   - Col.0 train 24 test 12.
   - Col.0.PhyD train 24 test 12.
   - Wassilewskija.ecotype train 24 test 12.
5. Keep OSD-918 as a blocker/status proof until local AnnData exists.

Production-slide rule:

- One proof object dominates the slide.
- One method or claim-boundary module supports it.
- Source URL/date/status must remain visible.
- Placeholder objects must not be interpreted as measured evidence.
