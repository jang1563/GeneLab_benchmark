# Human Organoid Response-Signature Scorer

Status: draft diagnostic scorer implemented
Date: 2026-05-22
Task: `draft_human_organoid_spaceflight`
Block: `V9-ORG-015`

## Scope

The evaluator can now compute organoid response-signature diagnostics when a
submission includes:

```text
response_signature.csv
```

The scorer joins that artifact to the derived DE reference table from
V9-ORG-014:

```text
v9/human_organoid/de_references/human_organoid_de_reference.draft.csv.gz
```

Classification metrics remain primary. These response-signature metrics are
diagnostic because the organoid task is still small and confounded by source,
organoid fate, disease context, microglia condition, and donor/iPSC-line.

## Required Response Columns

```text
task_id
source_id
contrast_id
gene_symbol
ensembl_id
predicted_log2fc_leo_or_iss_minus_ground
```

Validation rules:

- `task_id` must match the task manifest when present.
- `predicted_log2fc_leo_or_iss_minus_ground` must be numeric.
- each row must include `gene_symbol` or `ensembl_id`.
- missing, unreadable, empty, or invalid artifacts skip DE/signature metrics
  with explicit reasons.

## Join Contract

Primary join:

```text
source_id + contrast_id + gene_symbol
```

Fallback join:

```text
source_id + contrast_id + ensembl_id
```

Duplicate response keys are skipped after the first joined row so a model cannot
inflate a metric by repeating the same feature.

## Metrics

`de_direction_match`:

- uses joined reference rows where `significant_fdr_0_05 == true`;
- skips zero-direction rows;
- reports the fraction where predicted and reference log2FC signs match.

`signature_rank_correlation`:

- uses all joined rows;
- computes Spearman correlation between predicted and reference log2FC;
- skips when fewer than two joined rows or no rank variation is available.

Both metrics include aggregate details, join counts, reference table summary,
validation details, and per-contrast diagnostics in `metrics.json`.

## CLI

```bash
python scripts/evaluate_v9_submission.py \
  --task-manifest v9/human_organoid/task_manifests/draft_human_organoid_spaceflight.json \
  --submission path/to/predictions.csv \
  --response-signature path/to/response_signature.csv \
  --output-dir path/to/report
```

The evaluator reads the DE reference path from the task manifest. A custom
reference can be supplied for tests or dry runs:

```bash
python scripts/evaluate_v9_submission.py \
  --task-manifest path/to/task.json \
  --submission path/to/predictions.csv \
  --response-signature path/to/response_signature.csv \
  --reference-signature-table path/to/reference.csv.gz \
  --output-dir path/to/report
```

## Implementation

- `spacebio_bench/signature_metrics.py`
- `spacebio_bench/evaluate.py`
- `scripts/evaluate_v9_submission.py`
- `tests/test_v9_spacebio_bench.py`

## Boundary

The current nearest-centroid organoid baseline does not emit gene-level
response signatures, so its existing outputs still skip these diagnostics. A
future model adapter can opt in by writing `response_signature.csv`.

## V9-ORG-016 Smoke Test

The scorer has an end-to-end smoke report at:

```text
v9/human_organoid/reports/response_signature_smoke/
```

That report uses the real derived DE reference table, builds a small
`response_signature.csv` by mirroring selected reference rows, and records
computed `de_direction_match` and `signature_rank_correlation` values. Because
the response signature is constructed from the reference itself, the values are
only a join/scoring plumbing check and are not model-performance evidence.

## Source-Transfer Baseline

`V9-ORG-017` selected a conservative source-transfer empirical signature
adapter design. See:

- `v9/human_organoid/reports/ORGANOID_RESPONSE_SIGNATURE_ADAPTER_DESIGN.md`

`V9-ORG-018` implemented that adapter and generated a diagnostic-only
source-transfer response-signature report at:

- `v9/human_organoid/reports/source_transfer_signature/`

## Next Step

`V9-ORG-019` should review the source-transfer diagnostic details and decide
whether the next adapter should add per-condition training signatures,
classifier-coefficient signatures, or keep the source-transfer baseline as the
first conservative diagnostic.
