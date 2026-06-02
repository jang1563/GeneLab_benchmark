# BioVis OSD-120 audience split proof v0.2 review

Date: 2026-06-02

## Anchor

This continues the OSDR source-record proof-panel work. The previous source
panel successfully connected the official OSD-120 page to a local split proof,
but the embedded local proof still used internal analysis wording:

- "within-source diagnostic"
- "hidden"
- "interaction split"

v0.2 creates an audience-facing replacement that preserves the same
manifest-derived sample counts while lowering the visible terminology for
presentation use.

## Generated Outputs

Builder:

- `scripts/build_biovis_osd120_audience_split_proof_v2.py`

Output root:

- `output/biovis_osd120_audience_split_proof_v0_2/`

Generated proof crop:

- `output/biovis_osd120_audience_split_proof_v0_2/proof/osd120_same_study_task_check_proof.png`

Generated slide-scale panel:

- `output/biovis_osd120_audience_split_proof_v0_2/panels/01_dark_osd120_clean_source_to_task.png`

Generated QA:

- `output/biovis_osd120_audience_split_proof_v0_2/qa/osd120_same_study_task_check_proof_grayscale.png`
- `output/biovis_osd120_audience_split_proof_v0_2/qa/01_dark_osd120_clean_source_to_task_grayscale.png`
- `output/biovis_osd120_audience_split_proof_v0_2/osd120_audience_split_proof_manifest.json`
- `output/biovis_osd120_audience_split_proof_v0_2/osd120_audience_split_proof_qa.json`

## Evidence Source

Task manifest:

- `v9/multispecies/interaction_task_manifests/draft_osd120_arabidopsis_root_light_interaction_spaceflight.json`

Official source-record screenshot:

- `output/biovis_osdr_source_record_captures_v0_1/osd_120_viewport.png`

Official OSDR source record:

- `https://osdr.nasa.gov/bio/repo/data/studies/OSD-120`

## Terminology Cleanup

Old internal wording:

- "OSD-120 interaction split proof"
- "within-source diagnostic"
- "hidden"
- "not leave-one-mission-out"

New audience-facing wording:

- "OSD-120 same-study task check"
- "same OSDR study"
- "held out"
- "not mission-held-out"
- "not cross-species"

The goal is not to remove scientific precision. The goal is to make the method
legible before asking the viewer to understand benchmark design details.

## Visual Review

Verdict: stronger than the previous OSD-120 proof crop for a main-deck methods
bridge.

What improved:

- The title now tells the viewer the role of the object: a same-study task
  check.
- The train/test geometry is still visible and sample-count-backed.
- The audience can follow "train on these samples, test on the held-out group"
  without first parsing benchmark jargon.
- The claim boundary is visible in both the local proof crop and the 4K
  source-to-task panel.
- Grayscale QA preserves the main structure; color is not the only carrier of
  meaning.

Remaining caveats:

- The proof crop is still a methods bridge, not a scored result figure.
- "Mission-held-out" and "cross-species" remain technical terms, but they are
  necessary claim-boundary language and are paired with plain wording.
- The plant motif is intentionally secondary. It should not become the dominant
  evidence object.

## Recommended Use

Use `01_dark_osd120_clean_source_to_task.png` as the current OSD-120 methods
bridge candidate.

Deck placement:

1. Place the 4K PNG as the Z2 evidence scene.
2. Add editable slide text for the narrative headline and one sentence of
   interpretation.
3. Keep the source URL and draft status in the footer.
4. Do not present this as cross-species or leave-one-mission-out evidence.

## Next Step

Apply the same audience-facing cleanup to the organoid matrix proof:

- keep source IDs, dimensions, and hashes visible,
- reduce internal wording in the caption layer,
- make "official source page to downloaded matrix" the primary viewer path,
- preserve the draft/pilot claim boundary.
