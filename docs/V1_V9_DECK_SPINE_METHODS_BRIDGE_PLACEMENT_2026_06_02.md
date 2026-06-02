# v1-v9 deck spine methods-bridge placement

Date: 2026-06-02

## Anchor

This follows the methods-bridge shortlist commit and returns the work to the
original v1-v9 presentation/manuscript spine.

Sources reviewed:

- `docs/V1_V9_PRESENTATION_AND_MANUSCRIPT_MASTER_OUTLINE_2026_05_31.md`
- `docs/V1_V9_SLIDE_ASSET_MANIFEST_2026_05_31.md`
- `docs/VISUAL_TECHNICAL_PRODUCTION_PROTOCOL_2026_06_01.md`
- `docs/PREMIUM_BRIDGE_METHODS_NARRATION_PACK_B1_B4_2026_06_02.md`
- `docs/BIOVIS_METHODS_BRIDGE_SHORTLIST_REVIEW_2026_06_02.md`

## Verdict

No deck-spine drift if the new organoid and OSD-120 bridge assets stay in the
v9 extension section.

They should not replace:

- the B1-B4 core methods bridge,
- the v1-v7 benchmark result core,
- the v8 hypothesis-only boundary,
- the v9 public bulk alpha boundary.

The new bridge assets are useful because they make extension-track data
collection and task construction visible to first-time viewers. They become
drift if they start carrying the main scientific result story.

## Updated Placement Logic

### Opening and core methods

Use B1-B4 before performance results:

1. B1: public studies need an evaluation layer.
2. B2: source records become benchmark tasks.
3. B3: one whole mission is held out.
4. B4: feature processing and model fitting stay on the training side.

Role:

- explain data collection, task definition, held-out evaluation, and train-only
  processing;
- earn trust before any AUROC/model-comparison slide;
- keep organoid and OSD-120 out of the opening unless used as tiny examples in
  speaker notes.

### v1-v7 benchmark result core

Keep this as the center of gravity:

- tissue hierarchy,
- pathway rescue and confounder resistance,
- model comparison,
- v4-v7 hardening,
- temporal/single-cell lessons,
- negative results,
- biological interpretation,
- human translation,
- v7.1 public boundary.

Do not insert organoid or OSD-120 bridge slides here. They are extension
readiness examples, not completed v1-v7 benchmark result claims.

### v8 incubator

Keep v8 as hypothesis generation:

- SpaceMed bridge,
- stressor decomposition,
- intervention/causal boundary.

Do not use organoid or OSD-120 visuals to make v8 look more validated.

### v9 platform and extension section

This is where the selected bridge assets belong.

Recommended long-talk placement:

| Current master slide | Existing role | Update |
|---:|---|---|
| 18 | v9 Platform Architecture | Keep package/manifest/evaluator/run-manifest architecture first. |
| 19 | Public Bulk Alpha | Keep metadata-only alpha and payload-hash blocker visible. |
| 20 | Organoid Track | Use the v0.2 organoid source-to-matrix bridge as Z2 evidence surface if time allows. |
| 21 | Multispecies Track | Use the v0.2 OSD-120 source-to-task bridge as Z2 evidence surface only as a constrained same-study example. |
| 22 | Single-Cell Track | Use a blocker/status visual, not the organoid/OSD-120 proof panels. |
| 23 | Claimed / Not Claimed | Re-state draft/pilot and same-study boundaries. |
| 24 | Roadmap | Move extension receipts and superseded panels to appendix. |

Recommended short-talk placement:

- Keep slide 11 as one extension-track status slide.
- Use at most one selected bridge asset as a background/proof crop.
- Prefer the organoid source-to-matrix bridge if the talk needs human relevance.
- Prefer the OSD-120 source-to-task bridge if the talk needs task-construction
  rigor.
- Do not include both as separate short-talk slides unless the v9 extension
  section is the talk's main topic.

## Selected Asset Placement

Main-deck candidates:

- `output/biovis_organoid_audience_matrix_proof_v0_2/panels/01_dark_organoid_clean_source_to_matrix.png`
  - long talk: slide 20 evidence scene;
  - short talk: optional extension-track proof object;
  - manuscript: supplementary methods/provenance unless v9 organoid extension
    becomes a formal result section.
- `output/biovis_osd120_audience_split_proof_v0_2/panels/01_dark_osd120_clean_source_to_task.png`
  - long talk: slide 21 evidence scene;
  - short talk: optional task-construction proof object;
  - manuscript: supplementary methods/provenance unless v9 multispecies
    extension becomes a formal result section.

Appendix / reviewer backup:

- `output/biovis_osdr_source_record_proof_panel_v0_1/panels/01_osdr_source_record_contact_sheet.png`
- `output/biovis_osdr_source_record_proof_panel_v0_1/panels/04_dark_extension_source_lanes.png`
- raw OSDR captures under `output/biovis_osdr_source_record_captures_v0_1/`

Superseded:

- `output/biovis_osdr_source_record_proof_panel_v0_1/panels/02_dark_organoid_source_record_proof.png`
- `output/biovis_osdr_source_record_proof_panel_v0_1/panels/03_dark_osd120_source_record_proof.png`

## Caption And Overlay Rules

### Organoid bridge

Visible headline:

> Human organoid source records become local expression matrices.

Short caption:

> Official OSDR source records and downloaded normalized-count matrices support
> a draft human neural organoid extension lane. This supports source-to-matrix
> readiness only; it is not mission-held-out benchmark evidence.

Keep in speaker notes:

- OSD-863 and OSD-871 details;
- cortical and dopaminergic organoid distinctions;
- matrix row/sample dimensions;
- draft pilot status.

### OSD-120 bridge

Visible headline:

> One source study can define a constrained same-study task check.

Short caption:

> OSD-120 is used as a constrained Arabidopsis root RNA-seq task check with
> held-out groups inside the same source study. This is not mission-held-out or
> cross-species generalization evidence.

Keep in speaker notes:

- genotype/ecotype and light-treatment details;
- sample-balance explanation;
- why this is useful as a task-construction check but not as a headline result.

## Z-Stack Implementation Rule

For final slides, do not paste the contact sheet into the deck.

Use each selected bridge as:

- Z0: dark scientific canvas/atmosphere;
- Z1: subtle source-to-task or source-to-matrix measurement lane;
- Z2: selected PNG proof scene;
- Z3: editable headline, one focus ring, one transfer path;
- Z4: visible claim-boundary tag and source strip;
- Z5: optional motion/focus crop only if the slide is animated.

The current PNGs are strong evidence scenes, but the final deck should still
carry editable interpretation text. Do not bury the slide logic inside a
non-editable screenshot.

## Remaining Gaps Before Deck Assembly

High priority:

- feature-layer bridge explaining genes versus pathway summaries;
- model/metric caption explaining direct shared-row comparisons;
- v9 public bulk release-boundary slide;
- slide-level insertion QA for selected organoid and OSD-120 bridges.

Do not solve these by adding more organoid or OSD-120 proof variants. The next
design work should move laterally across the deck spine, not deeper into one
extension branch.

## Bottom Line

The methods-bridge sprint is now aligned and bounded.

Next deck work should be:

1. build a slide-order contact sheet around the full v1-v9 spine;
2. place selected bridge assets only at slides 20-21 or appendix;
3. return visual polish effort to the missing core bridge and result slides.
