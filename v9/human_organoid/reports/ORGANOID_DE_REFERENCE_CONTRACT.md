# Human Organoid DE Reference Contract

Status: draft metric-input contract
Date: 2026-05-21
Task: `draft_human_organoid_spaceflight`
Block: `V9-ORG-014`

## Decision

The draft human organoid task now has a derived, checksummed DE reference
artifact for future response-signature metrics:

- `v9/human_organoid/de_references/human_organoid_de_reference.draft.csv.gz`
- `v9/human_organoid/de_references/human_organoid_de_reference_manifest.draft.json`
- `scripts/build_v9_human_organoid_de_reference.py`
- `spacebio_bench/organoid_de_reference.py`

The reference remains `draft_not_frozen`. It is metric-ready as an input
contract, but DE/signature metrics remain non-primary until a scorer is promoted
and benchmark submissions include a valid `response_signature.csv` artifact.

## Source Freeze Rule

Use the public OSDR non-rRNArm bulk RNA-seq DE tables:

- OSD-863 / GLDS-716:
  `GLDS-716_rna_seq_differential_expression_GLbulkRNAseq.csv`
- OSD-871 / GLDS-720:
  `GLDS-720_rna_seq_differential_expression_GLbulkRNAseq.csv`

Exclude the `rRNArm` supplementary DE tables from the benchmark reference.
The generated manifest records the source download URLs, OSDR-reported file
sizes, and SHA-256 hashes of the downloaded DE CSV payloads.

## Contrast Freeze Rule

For each source, keep the four direct matched contrast labels where disease
context and microglia condition match and the contrast is:

```text
Ground Control v Space Flight
```

The resulting artifact has:

- 2 sources.
- 8 direct source-specific contrasts.
- 242,708 reference rows.
- 2,368 rows with `adj_p_value <= 0.05`.

## Orientation Contract

The row-level metric orientation is:

```text
LEO_or_ISS_minus_Ground
```

OSDR direct contrast columns are labeled as `Log2fc_(Ground Control)v(Space
Flight)`. The derived table treats the source label as `group_a - group_b` and
negates the source value so positive values mean higher expression in LEO/ISS
than matched ground control.

This sign convention is explicit in the manifest:

- source: `group_a_minus_group_b_assumed_from_Log2fc_A_v_B_header`
- transform: `negated_ground_control_v_space_flight_to_leo_or_iss_minus_ground`
- benchmark: `LEO_or_ISS_minus_Ground`

## Row Schema

The compressed reference CSV keeps only row-level fields needed for metric
joins and scoring:

```text
source_id,glds_prefix,task_id,organoid_type,disease_context,
microglia_condition,contrast_id,ensembl_id,gene_symbol,entrez_id,
log2fc_leo_or_iss_minus_ground,source_log2fc_value,stat,p_value,
adj_p_value,significant_fdr_0_05,source_de_file
```

Contrast labels, source groups, source URLs, SHA-256 hashes, and orientation
metadata are stored once in the JSON manifest rather than repeated for every
gene row.

## Response Signature Contract

Future DE/signature scoring expects a separate submission-side artifact:

```text
response_signature.csv
```

Required columns:

```text
task_id
source_id
contrast_id
gene_symbol
ensembl_id
predicted_log2fc_leo_or_iss_minus_ground
```

Primary join key:

```text
source_id + contrast_id + gene_symbol
```

Fallback join key:

```text
source_id + contrast_id + ensembl_id
```

Optional future columns may include `fold_id`, `method_id`, `prediction_rank`,
`prediction_uncertainty`, and `feature_namespace`, but they are not required for
the first contract.

## Evaluator Policy

The evaluator now records declared organoid DE/signature metrics as skipped
with a precise reason when `response_signature.csv` is missing:

- `de_direction_match`
- `signature_rank_correlation`

This prevents the metrics from silently disappearing from `metrics.json` while
keeping current classification baselines valid.

V9-ORG-015 adds diagnostic scoring when a valid `response_signature.csv` is
supplied. See:

- `v9/human_organoid/reports/ORGANOID_RESPONSE_SIGNATURE_SCORER.md`

## Boundary

This contract does not make DE/signature scoring a primary leaderboard metric.
The organoid pilot still has coupled source, organoid fate, disease, microglia,
and donor/iPSC-line factors. A future scorer should report these metrics as
response-signature diagnostics first, not as standalone evidence of
donor-generalization or mission-held-out performance.

## Next Step

`V9-ORG-016` added the dry-run/example response-signature report under
`v9/human_organoid/reports/response_signature_smoke/`. The next implementation
step is `V9-ORG-017`: design the first real model-produced response-signature
adapter without leaking the DE reference into predictions.
