# BioVis source-object proof crops

Date: 2026-06-02

## Anchor

This continues the v0.4 source-proof replacement sprint. The previous step
created an inventory of candidate real objects. This step creates the first
local proof crops from actual project evidence.

The goal is to replace schematic placeholders with proof objects that come from
real local matrix files and task manifests.

## Generated Outputs

Builder:

- `scripts/build_biovis_source_object_proof_crops.py`

Output root:

- `output/biovis_source_object_proof_crops_v0_1/`

Proof crops:

- `output/biovis_source_object_proof_crops_v0_1/proof_crops/organoid_matrix_audit_proof.png`
- `output/biovis_source_object_proof_crops_v0_1/proof_crops/osd120_interaction_split_proof.png`
- `output/biovis_source_object_proof_crops_v0_1/proof_crops/organoid_blocked_split_proof.png`

Production replacement mocks:

- `output/biovis_source_object_proof_crops_v0_1/production_mocks/01_dark_organoid_matrix_proof_replacement.png`
- `output/biovis_source_object_proof_crops_v0_1/production_mocks/02_dark_osd120_split_proof_replacement.png`
- `output/biovis_source_object_proof_crops_v0_1/production_mocks/03_proof_crop_contact_sheet.png`

Metadata:

- `output/biovis_source_object_proof_crops_v0_1/source_object_proof_crops_manifest.json`
- `output/biovis_source_object_proof_crops_v0_1/source_object_proof_crops_qa.json`

## Evidence Sources

Human organoid matrix proof:

- `v9/human_organoid/matrices/GLDS-716_rna_seq_Normalized_Counts_GLbulkRNAseq.csv`
- `v9/human_organoid/matrices/GLDS-720_rna_seq_Normalized_Counts_GLbulkRNAseq.csv`
- `v9/human_organoid/task_manifests/draft_human_organoid_spaceflight.json`

OSD-120 interaction split proof:

- `v9/multispecies/interaction_task_manifests/draft_osd120_arabidopsis_root_light_interaction_spaceflight.json`

Organoid blocked split proof:

- `v9/human_organoid/task_manifests/draft_human_organoid_spaceflight.json`

## Visual Review

Verdict: ready for production-slide replacement test.

What improved:

- The organoid matrix proof now uses actual normalized-count values, not a
  placeholder heatmap.
- Matrix metadata includes source IDs, gene/sample dimensions, and compact
  hashes.
- OSD-120 split proof makes train/test separation visible before score claims.
- Status/caveat text remains visible beside the proof object.

Best current main-slide candidate:

- `osd120_interaction_split_proof`

Best appendix or methods-support candidate:

- `organoid_matrix_audit_proof`

Useful but more caveated:

- `organoid_blocked_split_proof`

## Claim Boundaries

The generated crops prove local evidence presence and manifest-derived split
geometry. They are not final validated manuscript result figures.

Required caveats:

- Organoid matrix proof: matrix presence and sample alignment are supported;
  full payload freeze and benchmark claims remain draft.
- OSD-120 split proof: within-source interaction diagnostic only; not
  leave-one-mission-out and not cross-species generalization.
- Organoid split proof: single ISS mission pilot; useful for methods
  explanation, not for mission-held-out claims.

## Next Step

Capture official OSDR source-record screenshots and pair them with these local
proof crops:

1. OSD-863 source record.
2. OSD-871 source record or related-study proof from OSD-863.
3. OSD-120 source record.

After that, build one final dark production slide:

- one dominant proof crop,
- one OSDR source-record proof,
- one editable interpretation layer,
- one visible caveat/status footer.
