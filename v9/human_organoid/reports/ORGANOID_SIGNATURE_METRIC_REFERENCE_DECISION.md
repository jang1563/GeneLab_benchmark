# Human Organoid DE/Signature Metric Reference Decision

Status: draft reference-policy decision; V9-ORG-015 diagnostic scorer added
Date: 2026-05-21
Task: `draft_human_organoid_spaceflight`

## Decision

Public, machine-readable differential-expression references are available for
the OSD-863/GLDS-716 cortical and OSD-871/GLDS-720 dopaminergic organoid
sources. They are strong enough to keep `de_direction_match` and
`signature_rank_correlation` in the draft organoid metric profile as
reference-backed future metrics.

They are not yet strong enough to make DE/signature scoring a primary
leaderboard metric in the current pilot baseline. The next implementation step
must freeze a contrast subset, log2 fold-change orientation, gene-id mapping,
and submission artifact format before those metrics are computed.

For now:

- Classification metrics remain primary for the current draft pilot.
- `de_direction_match` and `signature_rank_correlation` are declared as
  reference-backed and now have a draft contrast/signature input contract plus
  diagnostic scorer. They still remain non-primary until response-signature
  artifacts are produced by real model adapters and reviewed.
- Baseline outputs under `v9/human_organoid/reports/nearest_centroid/`,
  `sensitivity/`, and `donor_diagnostics/` should not claim DE/signature
  performance yet.

## Evidence

Local audit artifacts:

- `v9/human_organoid/signature_reference_audit.draft.csv`
- `v9/human_organoid/signature_reference_audit.draft.json`
- `v9/human_organoid/de_references/human_organoid_de_reference.draft.csv.gz`
- `v9/human_organoid/de_references/human_organoid_de_reference_manifest.draft.json`
- `spacebio_bench/organoid_signature_audit.py`
- `spacebio_bench/organoid_de_reference.py`
- `spacebio_bench/signature_metrics.py`
- `scripts/audit_v9_human_organoid_signature_references.py`
- `scripts/build_v9_human_organoid_de_reference.py`

External sources checked:

- OSDR file-list API for OSD-863:
  https://visualization.osdr.nasa.gov/biodata/api/v2/dataset/OSD-863/files/
- OSDR file-list API for OSD-871:
  https://visualization.osdr.nasa.gov/biodata/api/v2/dataset/OSD-871/files/
- OSDR file metadata for OSD-863 differential expression:
  https://visualization.osdr.nasa.gov/biodata/api/v2/dataset/OSD-863/file/GLDS-716_rna_seq_differential_expression_GLbulkRNAseq.csv/
- OSDR file metadata for OSD-871 differential expression:
  https://visualization.osdr.nasa.gov/biodata/api/v2/dataset/OSD-871/file/GLDS-720_rna_seq_differential_expression_GLbulkRNAseq.csv/
- GSE259421 metadata via OmicsDI:
  https://www.omicsdi.org/dataset/geo/GSE259421
- Marotta et al. 2024:
  https://academic.oup.com/stcltm/article/13/12/1186/7833382
- NASA OSDR overview:
  https://science.nasa.gov/biological-physical/data/osdr/
- NASA OSDR FAQ:
  https://science.nasa.gov/reference/osdr-faq/

## Audit Findings

OSD-863 / GLDS-716:

- OSDR file listing returned 338 files.
- Reference candidates: 3 files.
- DE reference files: 2.
  - `GLDS-716_rna_seq_differential_expression_GLbulkRNAseq.csv`
  - `GLDS-716_rna_seq_differential_expression_rRNArm_GLbulkRNAseq.csv`
- Contrast definition files: 1.
  - `GLDS-716_rna_seq_contrasts_GLbulkRNAseq.csv`
- Contrast definitions parsed: 56 contrast pairs across 8 groups.
- Direct matched Ground Control versus Space Flight contrasts: 4.
- Reversed matched Space Flight versus Ground Control contrasts: 4.
- Disease contexts: `no known diseases`, `primary progressive multiple sclerosis`.
- Microglia conditions: `with microglia`, `without microglia`.

OSD-871 / GLDS-720:

- OSDR file listing returned 402 files.
- Reference candidates: 3 files.
- DE reference files: 2.
  - `GLDS-720_rna_seq_differential_expression_GLbulkRNAseq.csv`
  - `GLDS-720_rna_seq_differential_expression_rRNArm_GLbulkRNAseq.csv`
- Contrast definition files: 1.
  - `GLDS-720_rna_seq_contrasts_GLbulkRNAseq.csv`
- Contrast definitions parsed: 56 contrast pairs across 8 groups.
- Direct matched Ground Control versus Space Flight contrasts: 4.
- Reversed matched Space Flight versus Ground Control contrasts: 4.
- Disease contexts: `no known diseases`, `Sporadic Parkinson disease`.
- Microglia conditions: `with microglia`, `without microglia`.

The peer-reviewed paper reports RNA-seq processing with featureCounts and
analysis using DESeq2, GSEA, gProfiler/GO, and TFEA. The OSDR file audit found
packaged DE tables and contrast definitions, but did not find separate packaged
GSEA/GO/TFEA result tables in the OSDR file listing. Therefore the first
metric-ready reference should be DE-table based, not enrichment-table based.

## Metric Policy

Use the following policy in v9:

1. Keep `balanced_accuracy`, `auroc`, and `calibration_error` as the only
   primary computed organoid metrics for the current pilot baseline.
2. Keep `de_direction_match` and `signature_rank_correlation` declared in the
   `genelab_organoid_pilot` metric profile.
3. Mark their reference state as
   `public_osdr_de_reference_tables_available_pending_contrast_freeze`.
4. Use the V9-ORG-014/V9-ORG-015 contract for diagnostic scorer implementation:
   - canonical non-rRNArm `*_differential_expression_GLbulkRNAseq.csv` tables;
   - eight direct matched source-specific Ground Control versus Space Flight
     contrasts;
   - log2 fold-change orientation normalized to `LEO_or_ISS - Ground`;
   - human gene symbol as the primary feature key with Ensembl fallback;
   - `response_signature.csv` as a separate gene-level artifact rather than
     trying to infer response signatures from `predictions.csv`.

## Candidate Frozen Contrasts

For each source, the audit found four direct matched spaceflight contrasts:

- healthy/no-known-disease, with microglia: Ground Control versus Space Flight;
- healthy/no-known-disease, without microglia: Ground Control versus Space Flight;
- disease context, with microglia: Ground Control versus Space Flight;
- disease context, without microglia: Ground Control versus Space Flight.

The contrast table also contains reversed versions. V9-ORG-014 uses only the
direct Ground Control versus Space Flight labels and negates source log2FC
values so the benchmark orientation is `LEO_or_ISS - Ground`.

## V9-ORG-014 Output

The frozen-input draft reference now exists at
`v9/human_organoid/de_references/human_organoid_de_reference.draft.csv.gz`.
Its manifest records two public DE source files, eight direct contrasts, 242,708
gene/contrast rows, and 2,368 rows with `adj_p_value <= 0.05`. The contract note
is `v9/human_organoid/reports/ORGANOID_DE_REFERENCE_CONTRACT.md`.

## V9-ORG-015 Output

The evaluator now computes diagnostic `de_direction_match` and
`signature_rank_correlation` when a valid `response_signature.csv` is supplied.
It validates required columns, joins to the DE reference by
`source_id + contrast_id + gene_symbol` with Ensembl fallback, reports aggregate
and per-contrast details, and still emits precise skip reasons when the artifact
is missing or invalid. The scorer note is
`v9/human_organoid/reports/ORGANOID_RESPONSE_SIGNATURE_SCORER.md`.

## V9-ORG-016 Output

The response-signature smoke report now exists at
`v9/human_organoid/reports/response_signature_smoke/`. It uses the real derived
DE reference table and a small mirrored fixture to verify that response
validation, reference joining, run-manifest provenance, and diagnostic metric
calculation work end to end. Because the response signature is derived from the
reference itself, this report is explicitly not a model-performance claim.

## V9-ORG-017 Output

The first model-produced adapter path is specified in
`v9/human_organoid/reports/ORGANOID_RESPONSE_SIGNATURE_ADAPTER_DESIGN.md`. The
recommended baseline is a source-transfer empirical response signature: predict
all OSD-863 target contrasts from OSD-871 training samples and all OSD-871
target contrasts from OSD-863 training samples, without using target-source
expression or DE references during signature generation.

## V9-ORG-018 Output

The source-transfer adapter report now exists at
`v9/human_organoid/reports/source_transfer_signature/`. It emits compressed
`response_signature.csv.gz`, uses `reference_not_used_for_signature_generation`,
and computes diagnostic signature metrics against the derived DE reference. The
draft diagnostic scores are `de_direction_match=0.7706734867860188` and
`signature_rank_correlation=0.1760078660242601`; these are source-transfer
diagnostics only, not leaderboard claims.

## Risk Boundary

The available reference tables are good enough for metric design, but not enough
to claim a stable benchmark score without another freeze step. Main risks:

- The DE tables are large source-level files, so a compact derived reference
  table should be generated and checksummed before scoring.
- The current task has source, organoid fate, disease context, and donor
  confounding. A DE metric can evaluate response-signature alignment, but it
  cannot by itself create donor-generalization or mission-held-out evidence.
- The nearest-centroid pilot does not emit gene-level response signatures, so
  its existing prediction outputs cannot be retrofitted into DE/signature
  scoring without an additional model-output artifact.

## Next Step

Open V9-ORG-019 as `Human organoid source-transfer diagnostic review`. It
should interpret the source-transfer signature report, inspect per-contrast
details, and decide whether to keep this adapter as the first real diagnostic
baseline.
