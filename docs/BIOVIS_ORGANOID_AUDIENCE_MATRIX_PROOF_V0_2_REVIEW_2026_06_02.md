# BioVis organoid audience matrix proof v0.2 review

Date: 2026-06-02

## Anchor

This continues the audience-facing cleanup pass started with the OSD-120 split
proof. The previous organoid matrix proof was evidence-backed, but it still read
like a reviewer/audit crop rather than a premium methods bridge.

v0.2 keeps the real normalized-count matrix evidence while making the visible
story easier for a first-time viewer:

> official OSDR source pages -> downloaded expression matrices -> draft/pilot
> claim boundary.

## Generated Outputs

Builder:

- `scripts/build_biovis_organoid_audience_matrix_proof_v2.py`

Output root:

- `output/biovis_organoid_audience_matrix_proof_v0_2/`

Generated proof crop:

- `output/biovis_organoid_audience_matrix_proof_v0_2/proof/human_organoid_downloaded_matrix_proof.png`

Generated slide-scale panel:

- `output/biovis_organoid_audience_matrix_proof_v0_2/panels/01_dark_organoid_clean_source_to_matrix.png`

Generated QA:

- `output/biovis_organoid_audience_matrix_proof_v0_2/qa/human_organoid_downloaded_matrix_proof_grayscale.png`
- `output/biovis_organoid_audience_matrix_proof_v0_2/qa/01_dark_organoid_clean_source_to_matrix_grayscale.png`
- `output/biovis_organoid_audience_matrix_proof_v0_2/organoid_audience_matrix_proof_manifest.json`
- `output/biovis_organoid_audience_matrix_proof_v0_2/organoid_audience_matrix_proof_qa.json`

## Evidence Source

Task manifest:

- `v9/human_organoid/task_manifests/draft_human_organoid_spaceflight.json`

Official source-record screenshots:

- `output/biovis_osdr_source_record_captures_v0_1/osd_863_viewport.png`
- `output/biovis_osdr_source_record_captures_v0_1/osd_871_viewport.png`

Official OSDR source records:

- `https://osdr.nasa.gov/bio/repo/data/studies/OSD-863`
- `https://osdr.nasa.gov/bio/repo/data/studies/OSD-871`

Local matrix files:

- `v9/human_organoid/matrices/GLDS-716_rna_seq_Normalized_Counts_GLbulkRNAseq.csv`
- `v9/human_organoid/matrices/GLDS-720_rna_seq_Normalized_Counts_GLbulkRNAseq.csv`

## Visible Terminology

Audience-facing terms:

- "Human organoid matrix proof"
- "Official OSDR pages"
- "Downloaded expression matrices"
- "Cortical organoids"
- "Dopaminergic organoids"
- "42 samples"
- "draft pilot"
- "not mission-held-out"

Internal/audit terms moved out of the visible layer:

- payload checksum status
- diagnostic scorer language
- detailed source-freeze status

The internal details still live in the manifest and project records. They do not
need to compete with the viewer path on a main-deck methods slide.

## Visual Review

Verdict: ready as the current organoid methods-bridge candidate.

What improved:

- The panel now reads as a source-to-matrix bridge instead of a standalone audit
  board.
- The two source records, OSD-863 and OSD-871, remain visible as official page
  proof objects.
- The matrix proof uses real normalized-count values from local matrices.
- Sample dimensions are prominent: 30,408 genes x 19 samples and 30,269 genes x
  23 samples.
- Hashes and file names are retained as provenance but visually subordinated.
- Claim boundary is visible in both the proof crop and 4K panel.
- Grayscale QA preserves the main hierarchy.

Fix made during visual QA:

- The first render placed gene/sample counts over the heatmap. The heatmap row
  height was reduced and metadata was moved into a separate lower band.

## Claim Boundary

This proof supports:

- official source-record continuity,
- local matrix presence,
- sample-column alignment,
- organoid extension readiness for methods discussion.

It does not support:

- final benchmark performance,
- leave-one-mission-out claims,
- completed payload freeze,
- final manuscript result claims.

## Recommended Use

Use `01_dark_organoid_clean_source_to_matrix.png` as the current organoid
source-to-matrix methods bridge.

Deck placement:

1. Place the 4K PNG as the Z2 evidence scene.
2. Add editable headline/callout text in the deck layer.
3. Keep OSDR source URLs and draft status visible in the footer.
4. Avoid presenting the organoid pilot as a finished benchmark result.

## Next Step

The two strongest methods bridge candidates are now:

- OSD-120: `output/biovis_osd120_audience_split_proof_v0_2/panels/01_dark_osd120_clean_source_to_task.png`
- Human organoids: `output/biovis_organoid_audience_matrix_proof_v0_2/panels/01_dark_organoid_clean_source_to_matrix.png`

Next, create a small "methods bridge shortlist" contact sheet and decide which
slides are main-deck, appendix, or reviewer-backup candidates before moving
into the actual slide-deck build.
