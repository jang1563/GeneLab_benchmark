# BioVis OSDR source-record capture review

Date: 2026-06-02

## Anchor

This continues the v0.4 source-proof replacement sprint. The prior step created
local proof crops from real matrix files and task manifests. This step captures
official OSDR source-record pages and builds slide-scale proof panels that
connect source pages to local analysis objects.

Design rule used here:

> Layered PNG scene + editable scientific interpretation.

The goal is to make data collection and task construction visually legible for
first-time viewers without turning the figure into a low-quality card collage.

## Official Source Captures

Capture script:

- `scripts/capture_osdr_source_records.js`

Capture output root:

- `output/biovis_osdr_source_record_captures_v0_1/`

Captured official OSDR pages:

- OSD-863: `https://osdr.nasa.gov/bio/repo/data/studies/OSD-863`
- OSD-871: `https://osdr.nasa.gov/bio/repo/data/studies/OSD-871`
- OSD-120: `https://osdr.nasa.gov/bio/repo/data/studies/OSD-120`
- OSD-918: `https://osdr.nasa.gov/bio/repo/data/studies/OSD-918`

Generated capture files:

- `output/biovis_osdr_source_record_captures_v0_1/osd_863_viewport.png`
- `output/biovis_osdr_source_record_captures_v0_1/osd_863_fullpage.png`
- `output/biovis_osdr_source_record_captures_v0_1/osd_871_viewport.png`
- `output/biovis_osdr_source_record_captures_v0_1/osd_871_fullpage.png`
- `output/biovis_osdr_source_record_captures_v0_1/osd_120_viewport.png`
- `output/biovis_osdr_source_record_captures_v0_1/osd_120_fullpage.png`
- `output/biovis_osdr_source_record_captures_v0_1/osd_918_viewport.png`
- `output/biovis_osdr_source_record_captures_v0_1/osd_918_fullpage.png`
- `output/biovis_osdr_source_record_captures_v0_1/osdr_source_record_capture_manifest.json`
- `output/biovis_osdr_source_record_captures_v0_1/osdr_source_record_capture_summary.json`

Capture status:

- OSD-863: captured; 4/4 expected terms visible.
- OSD-871: captured; 4/4 expected terms visible.
- OSD-120: captured; 3/3 expected terms visible.
- OSD-918: captured; 3/3 expected terms visible.

## Proof Panel Outputs

Panel builder:

- `scripts/build_biovis_osdr_source_record_proof_panel.py`

Panel output root:

- `output/biovis_osdr_source_record_proof_panel_v0_1/`

Generated source crops:

- `output/biovis_osdr_source_record_proof_panel_v0_1/crops/osd_863_source_record_content_crop.png`
- `output/biovis_osdr_source_record_proof_panel_v0_1/crops/osd_871_source_record_content_crop.png`
- `output/biovis_osdr_source_record_proof_panel_v0_1/crops/osd_120_source_record_content_crop.png`
- `output/biovis_osdr_source_record_proof_panel_v0_1/crops/osd_918_source_record_content_crop.png`

Generated slide-scale panels:

- `output/biovis_osdr_source_record_proof_panel_v0_1/panels/01_osdr_source_record_contact_sheet.png`
- `output/biovis_osdr_source_record_proof_panel_v0_1/panels/02_dark_organoid_source_record_proof.png`
- `output/biovis_osdr_source_record_proof_panel_v0_1/panels/03_dark_osd120_source_record_proof.png`
- `output/biovis_osdr_source_record_proof_panel_v0_1/panels/04_dark_extension_source_lanes.png`

Generated QA assets:

- `output/biovis_osdr_source_record_proof_panel_v0_1/qa/01_osdr_source_record_contact_sheet_grayscale.png`
- `output/biovis_osdr_source_record_proof_panel_v0_1/qa/02_dark_organoid_source_record_proof_grayscale.png`
- `output/biovis_osdr_source_record_proof_panel_v0_1/qa/03_dark_osd120_source_record_proof_grayscale.png`
- `output/biovis_osdr_source_record_proof_panel_v0_1/qa/04_dark_extension_source_lanes_grayscale.png`
- `output/biovis_osdr_source_record_proof_panel_v0_1/osdr_source_record_proof_panel_manifest.json`
- `output/biovis_osdr_source_record_proof_panel_v0_1/osdr_source_record_proof_panel_qa.json`

## Visual Review

Verdict: ready as draft main-deck source-proof candidates, with manuscript-use
caveats.

Best main-deck candidates:

- `02_dark_organoid_source_record_proof.png`
- `03_dark_osd120_source_record_proof.png`

Best appendix or review-board candidates:

- `01_osdr_source_record_contact_sheet.png`
- `04_dark_extension_source_lanes.png`

What improved:

- OSDR title, GeneLab ID, DOI, source accession, and page date are visible in
  the proof objects.
- Raw browser chrome is mostly removed from the production panels.
- The viewer can follow the path from official source page to local matrix or
  split object.
- Z-stack is explicit: dark canvas, measurement grid/orbit, source proof
  surface, local proof surface, interpretation callout, claim boundary footer.
- Biological motifs are present but secondary: organoid/cell cluster for human
  organoid panels and root motif for OSD-120.
- Color is reinforced by labels; grayscale QA preserves the main hierarchy.

## Text And Jargon Audit

Changes already made:

- The shell label uses "human neural organoids" instead of making the viewer
  parse iPSC terminology first.
- The OSD-120 shell uses "same-study diagnostic" rather than "within-source
  diagnostic."
- The local evidence label uses "Downloaded matrix proof" rather than an
  abstract internal object label.
- Each main panel includes a plain-language caution or "What this proves"
  block.

Remaining caveats:

- The embedded local OSD-120 split proof still contains analytical wording such
  as "within-source diagnostic" because it was generated in the previous local
  proof-crop step. This is acceptable for a methods/reviewer version, but a
  cleaner audience-facing variant should be generated before final slides.
- "Matrix" remains necessary technical language. It is now paired with a visual
  heatmap proof and plain-language source-to-download framing.
- OSD-918 supports single-cell source availability, not benchmark readiness.

## Claim Boundary

These panels prove:

- Official OSDR source-record page availability.
- Visible accession context for organoid, plant, and single-cell extension
  lanes.
- Continuity between source records and local proof objects.

These panels do not prove:

- Full payload completeness.
- Final benchmark validity.
- Mission-held-out or cross-species generalization.
- That every extension lane is ready for scoring.

## Recommended Next Step

1. Build a cleaner audience-facing OSD-120 split proof variant with "same-study
   task check" wording.
2. Convert `02_dark_organoid_source_record_proof.png` and
   `03_dark_osd120_source_record_proof.png` into slide layouts with PNG proof
   scene plus editable callout text.
3. Add permanent citation/source footer treatment for manuscript-style use.
4. Run final slide-level QA after placing these panels into the actual deck:
   grayscale, projector contrast, 16:9 crop, and small-screen readability.
