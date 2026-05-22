# v9 Organoid and Species Extension Review

Status: evidence review and planning input
Date: 2026-05-21
Scope: SpaceBio-Bench v9 extension beyond current mouse public bulk LOMO tasks

## Bottom line

The short prior summary is directionally correct, but it should be made more
precise before it influences an abstract, manuscript, or v9 implementation.

Human spaceflight organoid data are real, public, and relevant. The strongest
near-term candidates are the NASA/OSDR cortical and dopaminergic human
iPSC-derived neural organoid RNA-seq datasets linked to GEO GSE259421. A 2026
Scientific Data descriptor also confirms a separate human brain organoid
proteomics dataset under ISS and ground conditions.

However, these should not be merged into the current mouse bulk LOMO leaderboard
as if they were the same task family. They should become a separate
`human_organoid_spaceflight` or `cell_model_spaceflight` track with its own
schema fields, split logic, metrics, and conservative claims.

Mouse-only is not a hard limit for v9. NASA GeneLab/OSDR explicitly spans many
biological systems, and this repository already contains older Drosophila and
Arabidopsis data assets plus historical cross-species outputs. The correct v9
move is to extend source/task schemas first, then add non-mouse species as
separate task families or pathway/ortholog-level transfer tasks.

## Evidence checked

### Human neural organoid RNA-seq

Confirmed sources:

- NASA Open Data Portal, cortical organoids:
  https://data.nasa.gov/dataset/effects-of-microgravity-on-human-ipsc-derived-neural-organoids-on-the-international-space-
- NASA Open Data Portal, dopaminergic organoids:
  https://data.nasa.gov/dataset/effects-of-microgravity-on-human-ipsc-derived-neural-organoids-on-the-international-space--c4013
- GEO GSE259421:
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE259421
- Marotta et al., Stem Cells Translational Medicine 2024:
  https://academic.oup.com/stcltm/article/13/12/1186/7833382

Key verified facts:

- The NASA portal has two public dataset pages:
  - cortical organoids, cross-linked to dopaminergic data under OSD-871
  - dopaminergic organoids, cross-linked to cortical data under OSD-863
- GEO GSE259421 is public and lists:
  - organism: Homo sapiens
  - experiment type: high-throughput expression profiling
  - 42 samples
  - cortical and dopaminergic organoids
  - with and without iPSC-derived microglia
  - LEO/ISS versus ground controls
  - bulk RNA-seq on individual organoids
- The GEO GSE259421 series matrix also exposes sample titles, Subject ids, cell
  line ids, cell type, treatment, BioSample accessions, and SRA accessions for
  the 42 samples. In v9 these are parsed into
  `v9/human_organoid/geo_sample_metadata.draft.csv` and merged into
  `sample_factors.draft.csv`.
- The 2024 paper reports human iPSC-derived neural organoids maintained on the
  ISS with matched Earth controls and links its underlying data to GSE259421.

Implication for v9:

These data can support a public human organoid RNA-seq task-card pilot. The task
should be framed as a small, human cell-model spaceflight benchmark or case
study, not as another mouse tissue LOMO fold. Donor/Subject holdouts should
remain diagnostic-only unless a separate task family explicitly accounts for
the donor-source-organoid-disease coupling.

Implementation update:

- v9 now parses GEO-derived donor/Subject and iPSC-line metadata for all 42
  samples.
- The donor-aware split decision keeps Subject holdouts under
  `v9/human_organoid/reports/donor_diagnostics/` with
  `claim_boundary=donor_diagnostic_only_not_leaderboard`.
- No separate human organoid donor-generalization leaderboard task is currently
  recommended.
- v9 now audits OSDR DE/signature reference evidence at
  `v9/human_organoid/signature_reference_audit.draft.csv`. OSDR lists public
  differential-expression tables and contrast-definition CSVs for both
  OSD-863/GLDS-716 and OSD-871/GLDS-720. Each source has 56 parsed contrast
  pairs, including four direct matched Ground Control versus Space Flight
  contrasts and four reversed matches.
- The DE/signature decision note keeps `de_direction_match` and
  `signature_rank_correlation` reference-backed but non-primary until a frozen
  contrast/signature contract is added.

### Human brain organoid proteomics

Confirmed source:

- Martins et al., Scientific Data 2026:
  https://www.nature.com/articles/s41597-026-06881-5.pdf

Key verified facts:

- The descriptor reports mass-spectrometry proteomic profiles from human brain
  organoids.
- The design compares ISS and ground conditions.
- The organoids include MECP2-deficient Rett syndrome model and wild-type
  control iPSC lines.
- The article points to PRIDE PXD069807 and a Zenodo record for processed
  identification/quantitation data.

Implication for v9:

This is valuable but should be later than the RNA-seq organoid pilot because it
is proteomics, small-N, and a different modality. Treat it as a multimodal or
disease-model organoid case study, not a first public leaderboard task.

### Non-mouse species

Confirmed online sources:

- NASA GeneLab overview:
  https://www.nasa.gov/centers-and-facilities/ames/what-is-nasas-genelab/
- NASA cross-species GeneLab article:
  https://www.nasa.gov/osdr-latest-news-multi-omics-cross-species-analysis-of-genelab-data-leads-to-new-nasa-investigation/

Key verified facts:

- NASA describes GeneLab/OSDR space biology studies across microorganisms,
  fruit flies, rodents, plants, fish, and human cells.
- NASA describes an X-Species Transcriptional Explorer covering processed
  GeneLab transcriptional data across human, mouse, rat, fly, worm, yeast, rice,
  tomato, Arabidopsis, Brassica, and Brachypodium.

Repository evidence:

- Current v9 public bulk source inventory is mouse-focused and has 22 public
  source rows in `v9/source_inventory.csv`.
- The current inventory does not yet have explicit `species`, `organism`,
  `material_type`, `model_system`, or `assay_modality` fields.
- Local non-mouse assets already exist:
  - `data/multispecies/drosophila/GLDS-207_*`
  - `data/multispecies/arabidopsis/GLDS-37_*`
  - `data/multispecies/arabidopsis/GLDS-120_*`
- Historical scripts and outputs exist for multispecies/cross-species analyses:
  - `v3/scripts/fetch_multispecies_osdr.py`
  - `v3/evaluation/E4_multispecies_nes.json`
  - `v6/evaluation/V6_C_cross_species_transfer.json`
  - `v8/bridge/evaluation/species_transfer_nes.csv`

Implication for v9:

Species expansion is feasible, and this repo already has early material for it.
But raw gene-level mouse, human, fly, worm, plant, and microbe tasks should not
be collapsed into one leaderboard. Cross-species tasks need ortholog, pathway,
or signature-level feature definitions and a separate metric profile.

## What should change in v9 design

### Add schema fields before adding data

Before organoids or non-mouse species enter v9 manifests, source/task metadata
should include:

- `organism`
- `taxon_id`
- `species_common_name`
- `material_type`
- `model_system`
- `biospecimen_type`
- `assay_modality`
- `platform`
- `spaceflight_environment`
- `ground_control_type`
- `donor_or_strain_block`
- `orthology_strategy`, for cross-species tasks
- `feature_namespace`, for raw genes versus ortholog groups versus pathways

For organoids specifically:

- `organoid_type`
- `differentiation_protocol`
- `microglia_condition`
- `donor_subject`
- `disease_context`
- `culture_hardware`

### Keep task families separate

Recommended task-family split:

| Family | First sources | Core question | Notes |
|---|---|---|---|
| `bulk_lomo` | current 22 public mouse sources | Mission-held-out mouse tissue flight/ground prediction | Current v9 alpha core |
| `human_organoid_spaceflight` | OSD-863, OSD-871, GSE259421 | Can models recover LEO-vs-ground or organoid response signatures under donor/type/microglia blocking? | Public RNA-seq pilot |
| `organoid_multimodal_case` | Scientific Data 2026 MECP2 proteomics | How do ISS/ground and genotype interact in a small human organoid proteomics dataset? | Case-study, not broad leaderboard |
| `multispecies_spaceflight` | OSD-207, OSD-37, OSD-120 first | How stable are flight signatures within non-mouse species? | Needs species-specific tasks |
| `bridge_cross_species` | mouse + human summaries + non-mouse pathway signatures | Which pathway-level signatures transfer across species/domains? | Ortholog/pathway level only |

### Avoid overclaiming in the current abstract

If the current ASGSR abstract is already near the word limit and currently
reports mouse bulk tasks, do not add organoids as if analyses are already done.
That would invite questions about integration, sample size, species alignment,
and whether organoid baselines were run.

Safe abstract language if one future-facing sentence is needed:

> The v9 framework is also designed to extend beyond mouse bulk tissues to
> public human iPSC-derived neural organoid RNA-seq studies and multispecies
> GeneLab resources.

Unsafe language:

- Do not claim organoid benchmark results before manifests, payload checks, and
  baselines exist.
- Do not imply that mouse tissue and human organoids are directly comparable on
  the same leaderboard.
- Do not describe the MECP2 proteomics descriptor as already integrated into v9.

## Recommended implementation order

1. Keep current v9 alpha focused on mouse public bulk LOMO until source freeze,
   RO-Crate, and release manifest work are complete.
2. Add schema support for organism, material, model system, assay modality, and
   feature namespace.
3. Create a source inventory extension draft for OSD-863 and OSD-871.
4. Draft one human organoid RNA-seq task manifest without running a leaderboard:
   - target: LEO versus ground or response-signature recovery
   - blocking: donor, organoid type, and microglia condition
   - status: pilot, not frozen
5. Audit processed files and checksums for OSD-863/OSD-871 using the same source
   audit pattern already added for v9 public bulk.
6. Run a feasibility baseline only after the manifest has clear split rules.
7. Separately inventory OSD-207, OSD-37, and OSD-120 as non-mouse pilot sources.
8. Decide whether multispecies v9 starts as:
   - species-specific flight/ground classifiers, or
   - ortholog/pathway-level cross-species transfer.

## Proposed backlog items

### V9-ORG-001: Human organoid schema extension

Goal:

Extend v9 manifest/source schemas for organism, material type, model system,
organoid metadata, and assay modality.

Done when:

- Existing mouse bulk manifests still validate.
- A draft human organoid source record can express OSD-863 and OSD-871 without
  overloading tissue fields.

### V9-ORG-002: OSD-863/OSD-871 source inventory

Goal:

Build a draft source inventory for public human neural organoid RNA-seq.

Done when:

- OSDR accessions, GEO accession, publication link, organism, material type,
  organoid type, microglia condition, and processed-file status are represented.

### V9-ORG-003: Human organoid task-card draft

Goal:

Draft one task manifest for GSE259421/OSD-863/OSD-871.

Done when:

- The split strategy explicitly handles donor, organoid type, microglia
  condition, and the single-mission limitation.

### V9-MULTI-001: Non-mouse source schema and inventory

Goal:

Inventory Drosophila, Arabidopsis, and additional candidate OSDR species without
mixing them into mouse bulk LOMO.

Done when:

- OSD-207, OSD-37, and OSD-120 have draft source rows with organism, feature
  namespace, assay, and task-family assignment.

### V9-MULTI-002: Ortholog/pathway feature strategy

Goal:

Decide how cross-species tasks represent features.

Done when:

- The design states whether each task uses raw gene ids, ortholog groups,
  pathways, or NES-style signatures.

## Decision recommendation

Use organoids and other species, but not by widening the current mouse task in
place. The stronger v9 story is:

1. Mouse public bulk LOMO is the frozen, reproducible alpha core.
2. Human organoids are the first public human cell-model extension.
3. Drosophila/Arabidopsis and other species become a cross-species or
   multispecies branch after schema expansion.
4. Each branch has its own metric profile and claim boundary.
