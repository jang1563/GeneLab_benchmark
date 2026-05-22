# Human Organoid DE/Signature Metric Reference Decision

Status: draft reference-policy decision
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
  reference-backed but pending a frozen contrast/signature contract.
- Baseline outputs under `v9/human_organoid/reports/nearest_centroid/`,
  `sensitivity/`, and `donor_diagnostics/` should not claim DE/signature
  performance yet.

## Evidence

Local audit artifacts:

- `v9/human_organoid/signature_reference_audit.draft.csv`
- `v9/human_organoid/signature_reference_audit.draft.json`
- `spacebio_bench/organoid_signature_audit.py`
- `scripts/audit_v9_human_organoid_signature_references.py`

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
4. Do not score DE/signature metrics until V9-ORG-014 defines:
   - canonical DE table choice, likely non-rRNArm `*_differential_expression_GLbulkRNAseq.csv`;
   - frozen direct contrast subset;
   - log2 fold-change orientation normalized to `LEO_or_ISS - Ground`;
   - feature key, likely human gene symbol plus Ensembl fallback;
   - minimum reference-gene and significance thresholds;
   - submission artifact format, likely a separate `response_signature.csv`
     rather than trying to infer gene-level signatures from `predictions.csv`.

## Candidate Frozen Contrasts

For each source, the audit found four direct matched spaceflight contrasts:

- healthy/no-known-disease, with microglia: Ground Control versus Space Flight;
- healthy/no-known-disease, without microglia: Ground Control versus Space Flight;
- disease context, with microglia: Ground Control versus Space Flight;
- disease context, without microglia: Ground Control versus Space Flight.

The contrast table also contains reversed versions. V9-ORG-014 must define a
single sign convention and invert log2 fold-change values when needed.

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

Open V9-ORG-014 as `Human organoid frozen DE contrast extraction and signature
metric contract`. It should generate a small derived reference table from the
public OSDR DE references, define the submission schema for gene-level response
signatures, and add skip-aware evaluator plumbing.
