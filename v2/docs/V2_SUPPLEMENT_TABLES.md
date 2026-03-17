# Supplementary Tables

## Supplementary Table S1. Parse-aware decomposition of Tier 3 zero-shot LLM performance

Legend: Tier 3 evaluates zero-shot text classification of spaceflight versus ground control from top-50 gene z-score prompts. End-to-end AUROC reflects the benchmark submission exactly, with parse failures assigned the neutral fallback score used in the submission JSON (`0.5`). Parsed-only AUROC was computed only on samples for which the parser recovered an explicit answer. Main-text reporting should use end-to-end AUROC; the parsed-only metric is provided here to separate model discrimination from output-format robustness.

Abbreviations: AUROC, area under the receiver operating characteristic curve; LLM, large language model.

### S1A. Provider-level summary

| Provider | Parsed | Total | Parse rate | Mean end-to-end AUROC | Mean parsed-only AUROC |
|---|---:|---:|---:|---:|---:|
| DeepSeek-V3 | 549 | 549 | 100.0% | 0.507 | 0.507 |
| Gemini-2.5-Flash | 529 | 549 | 96.4% | 0.493 | 0.493 |
| Llama-3.3-70B (Groq) | 489 | 549 | 89.1% | 0.485 | 0.443 |

### S1B. Per-task detail

| Task | Tissue | Metric | DeepSeek-V3 | Gemini-2.5-Flash | Llama-3.3-70B (Groq) |
|---|---|---|---:|---:|---:|
| A1 | Liver | End-to-end AUROC | 0.451 | 0.475 | 0.470 |
| A1 | Liver | Parsed-only AUROC | 0.451 | 0.475 | 0.459 |
| A1 | Liver | Parse rate | 100.0% | 100.0% | 92.7% |
| A2 | Gastrocnemius | End-to-end AUROC | 0.535 | 0.461 | 0.480 |
| A2 | Gastrocnemius | Parsed-only AUROC | 0.535 | 0.461 | 0.214 |
| A2 | Gastrocnemius | Parse rate | 100.0% | 100.0% | 53.1% |
| A3 | Kidney | End-to-end AUROC | 0.621 | 0.539 | 0.531 |
| A3 | Kidney | Parsed-only AUROC | 0.621 | 0.539 | 0.533 |
| A3 | Kidney | Parse rate | 100.0% | 100.0% | 91.5% |
| A4 | Thymus | End-to-end AUROC | 0.459 | 0.521 | 0.479 |
| A4 | Thymus | Parsed-only AUROC | 0.459 | 0.538 | 0.460 |
| A4 | Thymus | Parse rate | 100.0% | 95.5% | 92.5% |
| A5 | Skin | End-to-end AUROC | 0.482 | 0.560 | 0.446 |
| A5 | Skin | Parsed-only AUROC | 0.482 | 0.568 | 0.463 |
| A5 | Skin | Parse rate | 100.0% | 88.2% | 97.1% |
| A6 | Eye | End-to-end AUROC | 0.496 | 0.404 | 0.504 |
| A6 | Eye | Parsed-only AUROC | 0.496 | 0.379 | 0.528 |
| A6 | Eye | Parse rate | 100.0% | 86.5% | 64.9% |

### Footnotes

1. DeepSeek-V3 required no parse fallback recovery; end-to-end and parsed-only metrics are therefore identical.
2. Gemini-2.5-Flash parse rate improved after conservative recovery of truncated narrative verdicts and a narrow set of recurring archived `A1` and `A2` truncation patterns.
3. Llama-3.3-70B (Groq) shows the largest divergence between end-to-end and parsed-only metrics because truncation frequently removed an explicit answer token before the response terminated.
4. Parse rates are reported as `parsed/total` over all held-out samples pooled across folds within each task.
