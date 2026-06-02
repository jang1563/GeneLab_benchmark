# Biological Motif Asset Pack v0.2 Review

Date: 2026-06-02

Purpose:

- improve the first biological symbol pack, which was useful but too
  pictogram-like for premium slide use;
- create higher-end biological motif plates that can sit inside the
  SpaceBio-Bench Z-stack as restrained scientific context;
- verify the assets both as a contact sheet and as a full-canvas usage test.

## Output

Asset root:

- `assets/biovis_motif_pack_v0_2/`

Generated assets:

- `assets/biovis_motif_pack_v0_2/svg/`
- `assets/biovis_motif_pack_v0_2/png/`
- `assets/biovis_motif_pack_v0_2/manifest.json`
- `assets/biovis_motif_pack_v0_2/manifest.csv`
- `assets/biovis_motif_pack_v0_2/biovis_motif_pack_contact_sheet.png`
- `assets/biovis_motif_pack_v0_2/biovis_motif_pack_usage_sheet.png`

Generator:

- `scripts/build_biovis_motif_asset_pack_v2.py`

## Why v0.1 Was Not Enough

v0.1 solved availability, not premium feel.

Problems:

- the assets read as simple icons;
- several objects felt like generic clip-art when enlarged;
- the contact sheet looked usable but did not prove the assets could sit inside
  a high-end scientific slide;
- the biology specificity was present, but the z-depth and scientific texture
  were too light.

## v0.2 Design Change

v0.2 reframes the assets as biological motif plates rather than icons.

Changes:

- landscape transparent SVGs instead of square icon-first assets;
- subtle measurement-field lines and biological glow layers;
- more depth through soft shadows and translucent overlapping structures;
- more slide-native categories:
  - cell state;
  - mitochondrial stress;
  - organoid spheroid;
  - single-cell droplet;
  - tissue contour;
  - epithelial barrier;
  - vascular microsection;
  - spatial spot tissue;
  - RNA/DNA/protein/pathway motifs;
  - brain, kidney, liver contexts;
  - sample-to-matrix assay context.

## Manual Visual QA

Contact sheet:

- pass: v0.2 looks less like a generic icon set than v0.1;
- pass: visual identity colors are consistent with the current deck tokens;
- conditional pass: organ-context motifs still remain symbolic and should not
  be used as large hero objects.

Usage sheet:

- pass: motifs sit better on a full-canvas methods/provenance surface;
- pass: the assets behave as biological context rather than card-box
  decoration;
- pass: the feature-layer grouping is immediately more biological than a plain
  workflow rail;
- conditional pass: these should be used as Z2/Z3 context or callout objects,
  not as source-derived proof plates.

## Recommended Use

Use v0.2 for:

- B1/B2 biological context accents;
- the upcoming feature-layer bridge;
- organoid extension explanation;
- tissue/result slide anchors;
- assay/task construction signposts.

Do not use v0.2 for:

- replacing real plots, source-derived images, or manuscript figure panels;
- large decorative hero art;
- claiming evidence for a modality that has not been analyzed;
- implying cell type, organ, or assay specificity without nearby labels.

## Best Immediate Next Use

The feature-layer bridge should use:

- `rna_expression_trace`;
- `pathway_program_field`;
- `protein_complex_field`;
- optionally `sample_to_matrix_assay`.

Narrative:

> Expression signals are measured at the molecular level, summarized into
> biological programs, and evaluated under the same mission-held-out and
> train-only rules.

This is the first place where the new motif pack can make the deck feel more
biological without overcrowding B1-B4.

## Remaining Gap For v0.3

v0.2 is good enough for premium context motifs, but not yet enough for
source-like biological proof surfaces.

Recommended v0.3:

- add species/modality motifs: mouse, rat, fly, fish, plant, scRNA-seq droplet,
  spatial spot grid, proteomics spectrum;
- create grayscale/manuscript-safe variants;
- create one or two generated or source-derived biological texture plates for
  high-impact section openers;
- test motifs directly inside the next feature-layer slide rather than only in
  a QA sheet.

## Verdict

Pass as a premium context asset pack.

v0.2 is a meaningful upgrade over v0.1. It should be used as restrained
biological visual vocabulary, while source-derived figures and generated
texture plates remain responsible for high-evidence or hero-slide roles.
