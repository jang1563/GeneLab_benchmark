# V2 Fix Follow-Up Log

Created: 2026-03-10

Purpose: track post-review fixes, rationale, regeneration status, and follow-up work for `v2`.

## Priority Order

1. J1 pipeline comparison
2. T2 temporal recovery interpretation
3. Tier 3 LLM evaluation separation
4. T3 age metadata and documentation consistency

## Status Snapshot

| Item | Status | Notes |
|---|---|---|
| J1 matched-sample evaluation | Completed | Recomputed on the same 9 RR-1 animals for both pipelines |
| T2 recovery interpretation | Completed | Renamed classifier block to descriptive shift and aligned DD-20/docs |
| Tier 3 parsed-only vs end-to-end | Completed | Parser hardened, submissions regenerated, `parse_metrics_*` emitted for all providers/tasks |
| T3 age consistency | Completed | RR-8 runsheet confirms 32-week / 10 to 12 week labels; docs patched to match |

## Change Log

### 2026-03-10

- Created follow-up log to keep implementation and regeneration history in one place.
- Started J1 fix design:
  - Current issue: ML comparison uses different test sets (`RR-1` full fold vs matched subset).
  - Target: evaluate both pipelines on the same matched animals only.
- Implemented J1 fix in `v2/scripts/pipeline_version_compare.py`:
  - GLDS-48 test set now restricted to the matched RR-1 carcass samples.
  - GLDS-168 test set remains the matched reprocessed subset.
  - `ml_performance` now records comparison scope and matched sample IDs.
- Regenerated `v2/evaluation/J1_pipeline_comparison.json`.
- Updated `v2/evaluation/V2_RESULTS_SUMMARY.md` J1 section.

- Implemented T2 interpretation fix:
  - Renamed `classification_recovery` to `classification_shift`.
  - Added explicit caveat that `FLT_ISS-T` scores are in-sample reference scores.
  - Updated DD-20 wording to match the code's direction-aware formula.
  - Aligned `V2_RESULTS_SUMMARY.md`, `V2_PAPER_CONTENT.md`, `T_temporal_summary.json`, and `evaluation/T2_recovery_summary.json`.
- Started Tier 3 fix:
  - `llm_call_api.py` no longer treats cached responses with errors or truncation as complete.
  - `llm_parse_responses.py` now writes companion `parse_metrics_*` JSONs with pooled end-to-end AUROC and parsed-only AUROC.
  - Added parse summary metadata into submission JSONs.
- Hardened Tier 3 parser to remove scientific Python dependencies during evaluation:
  - Added canonical task-directory mapping for `A1`-`A6`
  - Replaced pandas/sklearn AUROC evaluation with stdlib CSV parsing and rank-based AUROC
- Refined Tier 3 parser heuristics:
  - Added explicit `final answer` / `prediction` phrase detection
  - Added first-sentence narrative verdict recovery for truncated outputs such as `consistent with spaceflight`
  - Added narrow task-specific residual fallbacks for recurring archived `A1`/`A2` truncation patterns
  - Kept Groq-style analysis-only truncations as parse failures to avoid label guessing
- Regenerated Tier 3 outputs:
  - Updated all `v2/evaluation/submission_*_zeroshot_A*.json`
  - Created `v2/evaluation/parse_metrics_*_{A1..A6}.json`
  - Updated `v2/evaluation/V2_RESULTS_SUMMARY.md` to separate end-to-end vs parsed-only metrics
  - Formatted `v2/docs/V2_SUPPLEMENT_TABLES.md` Supplementary Table S1 in journal-style caption/legend/footnote form for parse-aware Tier 3 reporting
- Confirmed T3 age source of truth from `data/mouse/liver/RR-8/GLDS-379_rna_seq_bulkRNASeq_v2_runsheet.csv`:
  - OLD samples are labeled `32 week`
  - YNG samples are labeled `10 to 12 week`
- Patched RR-8 age wording in `v2/docs/DATA_CATALOG_V2.md` and `v2/docs/V2_PAPER_CONTENT.md` to match the runsheet metadata.

### J1 verification

- Matched comparison scope: same 9 RR-1 animals only.
- Updated result:
  - `glds48_auroc = 0.700`
  - `glds168_auroc = 0.600`
  - `delta = -0.100`
  - `n_test_48 = 9`
  - `n_test_168 = 9`
- Interpretation change:
  - Previous summary claim of "minimal downstream ML effect" is no longer supported.
  - Current matched-animal comparison indicates a meaningful performance shift.

## Regeneration Checklist

| Output | Needed | Status |
|---|---|---|
| `v2/evaluation/J1_pipeline_comparison.json` | Yes | Completed |
| `v2/evaluation/V2_RESULTS_SUMMARY.md` | Yes | Completed |
| `v2/evaluation/T_temporal_summary.json` | Yes | Completed (schema/interpretation patch) |
| `evaluation/T2_recovery_summary.json` | Yes | Completed |
| `v2/evaluation/submission_*_zeroshot_A*.json` | Yes | Completed |
| `v2/evaluation/parse_metrics_*_{A1..A6}.json` | Yes | Completed |

## Verification Notes

- J1 regenerated successfully after patching the comparison scope.
- Runtime note:
  - Local default Python lacked scientific packages.
  - J1 regeneration required `uv run --with numpy --with pandas --with scipy --with scikit-learn`.
- T2 note:
  - Metric values were unchanged.
  - Downstream JSON summaries were patched to reflect the new schema/interpretation without recomputing the underlying scores.

## Next Actions

1. Decide whether to stop at the current conservative parser or pursue broader tissue-marker heuristics for the remaining Gemini/Groq failures.
2. If needed, promote `V2_SUPPLEMENT_TABLES.md` Supplementary Table S1 into the manuscript submission package format.
3. Continue with any remaining non-Tier-3 follow-up items outside this log.

## Tier 3 verification

- Parser regeneration completed with base `python3`; no external scientific packages required.
- Provider-level parse and mean AUROC summary:
  - DeepSeek-V3: 549/549 parsed, mean end-to-end AUROC 0.507, mean parsed-only AUROC 0.507
  - Gemini-2.5-Flash: 529/549 parsed, mean end-to-end AUROC 0.493, mean parsed-only AUROC 0.493
  - Llama-3.3-70B (Groq): 489/549 parsed, mean end-to-end AUROC 0.485, mean parsed-only AUROC 0.443
- Interpretation change:
  - Provider comparisons now explicitly distinguish model discrimination from output-format robustness.
  - Decision: parsed-only Tier 3 metrics belong in supplement/reviewer material, not the main text.

## T3 metadata verification

- No remaining `56-57 wk` / `16-17 wk` wording was found in active `v2/docs/` content or `v2/evaluation/V2_RESULTS_SUMMARY.md` after patching.
