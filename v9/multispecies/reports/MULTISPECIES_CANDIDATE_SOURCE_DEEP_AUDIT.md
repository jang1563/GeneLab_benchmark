# Multispecies Candidate Source Deep Audit

Status: source deep audit complete
Date: 2026-05-23
Block: `V9-MULTI-003`

## Scope

This note audits the first three non-mouse v9 multispecies candidates:

- OSD-207 / GLDS-207, `Drosophila melanogaster`
- OSD-37 / GLDS-37, `Arabidopsis thaliana`
- OSD-120 / GLDS-120, `Arabidopsis thaliana`

The goal is to decide which candidates should move toward draft task manifests,
which should be deferred, and what evidence is still missing.

## External Checks

Primary source checks used for this audit:

- NASA OSDR overview and data/API context:
  https://science.nasa.gov/biological-physical/data/osdr/
- NASA OSDR API overview:
  https://www.nasa.gov/osdr-api/
- NASA GLDS-37 article page:
  https://www.nasa.gov/osdr-latest-news-variation-in-the-transcriptome-of-different-ecotypes-of-arabidopsis-thaliana-reveals-signatures-of-oxidative-stress-in-plant-responses-to-spaceflight/
- NASA GLDS-120 plant-growth article page:
  https://www.nasa.gov/osdr-latest-news-genelab-and-awgs-enable-knowledge-of-plant-growth-in-spaceflight/
- PubMed GLDS-120 root-tip transcriptome article:
  https://pubmed.ncbi.nlm.nih.gov/38935184/

OSDR study URLs for the three candidates are recorded in the source inventory,
but the web-rendered OSDR study pages are dynamic. The next implementation
block should query the OSDR file-list API directly and record checksums, just
as the mouse and organoid audits do.

## Local Evidence

Local processed count matrices and sample tables exist:

| Source | Count Matrix | Sample Table |
|---|---|---|
| OSD-207 | `data/multispecies/drosophila/GLDS-207_rna_seq_Normalized_Counts_GLbulkRNAseq.csv` | `data/multispecies/drosophila/GLDS-207_rna_seq_SampleTable_GLbulkRNAseq.csv` |
| OSD-37 | `data/multispecies/arabidopsis/GLDS-37_rna_seq_Normalized_Counts_GLbulkRNAseq.csv` | `data/multispecies/arabidopsis/GLDS-37_rna_seq_SampleTable_GLbulkRNAseq.csv` |
| OSD-120 | `data/multispecies/arabidopsis/GLDS-120_rna_seq_Normalized_Counts_GLbulkRNAseq.csv` | `data/multispecies/arabidopsis/GLDS-120_rna_seq_SampleTable_GLbulkRNAseq.csv` |

Observed local dimensions:

| Source | Organism | Feature Rows | Sample Columns | Sample Rows | Condition Strata |
|---|---|---:|---:|---:|---:|
| OSD-207 | `Drosophila melanogaster` | 15,999 | 32 | 32 | 8 |
| OSD-37 | `Arabidopsis thaliana` | 28,067 | 56 | 56 | 8 |
| OSD-120 | `Arabidopsis thaliana` | 24,740 | 36 | 36 | 12 |

The local fetch script describes the intended source roles:

- OSD-207: fly whole-body flight versus ground, four genotype blocks.
- OSD-37: Arabidopsis seedling-pool flight versus ground, four ecotype blocks.
- OSD-120: Arabidopsis root flight versus ground by light condition, three
  genotype groups.

Historical cross-species output also exists, but should not be promoted yet.
`v3/evaluation/E4_multispecies_nes.json` contains Drosophila GLDS-207 pathway
NES evidence with 222 pathways. Arabidopsis bridge evidence was limited by
cross-kingdom mapping, so species-native tasks should come first.

## Candidate: OSD-207

Local condition counts:

| Stratum | Ground | Flight |
|---|---:|---:|
| Canton.S / Sei.ts1 | 4 | 4 |
| Canton.S / Wild.Type | 4 | 4 |
| w1118 / KCNQ370 | 4 | 4 |
| w1118 / Wild.Type | 4 | 4 |

Decision: `go_after_source_audit`.

Recommended task type:

- `multispecies_spaceflight`
- species-native Drosophila flight-versus-ground classification
- feature namespace: `species_gene_ids`
- split/blocking: genotype-aware folds

Why:

- local matrix and sample table are complete and aligned by sample count;
- each genotype stratum has balanced 4 ground and 4 flight samples;
- Drosophila is a distinct non-mouse animal branch and useful for species
  breadth.

Do not use as:

- raw-gene cross-species bridge with mouse or human;
- primary pathway/NES bridge before source integrity is checked.

Missing evidence:

- OSDR API file-list and checksum audit;
- sample-factor parser that extracts genotype and flight status;
- matrix column alignment audit against parsed sample factors;
- decision on whether all genotypes are one task with genotype blocking or four
  within-genotype micro tasks.

## Candidate: OSD-37

Local condition counts:

| Ecotype | Ground | Flight |
|---|---:|---:|
| Col.0 | 8 | 8 |
| Cvi.0 | 6 | 6 |
| Ler.0 | 6 | 6 |
| Ws.2 | 8 | 8 |

Decision: `first_plant_go_after_source_audit`.

Recommended task type:

- `multispecies_spaceflight`
- species-native Arabidopsis seedling-pool flight-versus-ground classification
- feature namespace: `species_gene_ids`
- split/blocking: ecotype-aware folds

Why:

- largest and cleanest of the three local non-mouse candidates;
- balanced flight/ground design within four ecotypes;
- NASA page explicitly links GLDS-37 to Arabidopsis ecotype transcriptome
  variation in spaceflight and oxidative stress response.

Missing evidence:

- OSDR API file-list and checksum audit;
- sample-factor parser for ecotype and flight status;
- matrix column alignment audit;
- response-signature/reference audit to see whether public DE tables exist and
  are stable enough for a diagnostic scorer.

## Candidate: OSD-120

Local condition counts:

| Genotype / Ecotype | Light | Ground | Flight |
|---|---|---:|---:|
| Col.0 | Dark | 3 | 3 |
| Col.0 | Light | 3 | 3 |
| Col.0.PhyD | Dark | 3 | 3 |
| Col.0.PhyD | Light | 3 | 3 |
| Wassilewskija | Dark | 3 | 3 |
| Wassilewskija | Light | 3 | 3 |

Decision: `defer_to_second_plant_or_interaction_task`.

Recommended task type:

- not first species-native classifier;
- later Arabidopsis light/genotype interaction diagnostic;
- possible response-signature or pathway/NES task after OSD-37 parser exists.

Why:

- sample table is complete but has a 3 genotype/ecotype x 2 light x 2 flight
  design;
- flight versus ground is entangled with light-condition interpretation unless
  the split explicitly stratifies by light;
- NASA and PubMed pages support GLDS-120 as a root-tip/light-condition
  transcriptome dataset, which is scientifically valuable but more complex than
  OSD-37.

Missing evidence:

- OSDR API file-list and checksum audit;
- parser for genotype/ecotype, light treatment, and flight status;
- decision on whether the first task is light-stratified flight/ground,
  genotype-stratified flight/ground, or explicit flight x light interaction;
- response-signature/reference audit.

## Cross-Species Bridge Decision

Do not start with a cross-species raw-gene task.

Do not start with an ortholog-group bridge.

Recommended sequence:

1. Build species-native source-factor and matrix audits for OSD-37 and OSD-207.
2. Draft species-native task cards using `species_gene_ids`.
3. Add OSD-120 as a second plant interaction task after the parser exists.
4. Revisit pathway/NES bridge tasks after source-native integrity is checked.

Reason:

- the repository already saw Drosophila pathway/NES evidence, but cross-species
  concordance was not strong enough to anchor a benchmark claim;
- Arabidopsis-to-human or Arabidopsis-to-mouse ortholog coverage is expected to
  be sparse and should not be the first bridge surface;
- pathway/NES remains the best later bridge namespace, consistent with
  `docs/V9_MULTISPECIES_FEATURE_STRATEGY.md`.

## Go/Defer Summary

| Source | Decision | First Task Candidate | Priority |
|---|---|---|---|
| OSD-37 | go after source audit | Arabidopsis ecotype-blocked flight/ground | 1 |
| OSD-207 | go after source audit | Drosophila genotype-blocked flight/ground | 2 |
| OSD-120 | defer | Arabidopsis flight x light/genotype interaction | 3 |

## Next Block

`V9-MULTI-004: Multispecies sample-factor and matrix audit scaffold`

Expected work:

- add parser logic for OSD-37, OSD-207, and OSD-120 sample tables;
- write `v9/multispecies/sample_factors.draft.csv` and `.json`;
- write `v9/multispecies/expression_matrix_audit.draft.csv` and `.json`;
- verify matrix columns match parsed sample ids;
- keep OSD-37 and OSD-207 as go candidates, OSD-120 as deferred until split
  design.
