# Text LLM Evaluation Format — GeneLab Benchmark (Tier 3)

**Version**: 1.0
**Design Decision**: DD-16
**Status**: Phase 2 specification (implementation pending)

---

## Overview

Tier 3 evaluates **text-based LLMs** (GPT-4o, Claude, Llama 3) on spaceflight detection without access to raw expression values. The model receives a natural-language description of a sample's top variable genes and must classify it as Flight or Ground.

This format tests whether LLMs can leverage **domain prior knowledge** about spaceflight biology to make predictions — complementing Tier 1 (classical ML on tabular data) and Tier 2 (Geneformer on rank-encoded sequences).

---

## Input Prompt Format

Each sample is encoded as a prompt using the top-N most variably expressed genes (selected from the LOMO training set variance filter). Gene names are given as `symbol (ENSMUSG...)` pairs.

### System Prompt

```
You are a bioinformatics expert specializing in spaceflight transcriptomics and mouse gene expression.
```

### User Prompt Template

```
A mouse {tissue} RNA-seq sample shows the highest expression variability in the following {n_genes} genes
(listed in descending order of variance across the training cohort):

{gene_1} ({ensembl_1}), {gene_2} ({ensembl_2}), ..., {gene_N} ({ensembl_N})

Based on your knowledge of spaceflight biology, predict whether this sample is from:
(A) Flight — collected from mice aboard the International Space Station or similar microgravity environment
(B) Ground — collected from ground control mice under normal gravity

Answer with the letter A or B, followed by a confidence score between 0.0 (uncertain) and 1.0 (certain).
Format: "<A or B> <confidence>"
Example: "A 0.82"
```

**Variables**:

| Variable | Description |
|----------|-------------|
| `{tissue}` | `thymus` or `gastrocnemius` |
| `{n_genes}` | Number of genes included (default: 50) |
| `{gene_i}` | Gene symbol (e.g., `Fbxo5`) |
| `{ensembl_i}` | Ensembl mouse gene ID (e.g., `ENSMUSG00000019773`) |

---

## Gene Selection Protocol

Genes are selected per-fold using the **LOMO-aware variance filter** applied to training samples only:

1. Compute per-gene variance across all training samples in the fold
2. Rank genes by variance (descending)
3. Select top `n_genes` (default: 50) from the same `selected_genes.txt` used for classical ML

**Critical rule**: Genes are selected from the *training set only* — the LLM must not receive information about test-set-specific expression patterns. This preserves the LOMO independence guarantee (DD-04).

Gene order in the prompt reflects variance rank: gene_1 = highest variance.

---

## Evaluation Variants

Three prompt strategies are evaluated, all using the same JSON submission format:

### 1. Zero-Shot (default)

No labeled examples provided. Only the gene list and task description.

```
{user_prompt_template}
```

### 2. Few-Shot (k = 3, 5, 10)

k labeled examples prepended before the test prompt. Examples are drawn from the *training set* of the same fold.

```
Here are {k} example samples with known labels:

Example 1 (Flight):
Genes: {gene_A1} ({ensembl_A1}), {gene_A2} ({ensembl_A2}), ...

Example 2 (Ground):
Genes: {gene_B1} ({ensembl_B1}), ...

...

Now predict the following sample:
{user_prompt_template}
```

### 3. Chain-of-Thought (CoT)

Request step-by-step reasoning before the final answer.

```
{user_prompt_template}

Before answering, briefly explain which genes (if any) you recognize as relevant to spaceflight biology,
and how they inform your prediction. Then provide your final answer on the last line as "<A or B> <confidence>".
```

---

## Output Parsing

The LLM response must contain:
1. A letter `A` (Flight) or `B` (Ground)
2. A confidence score in `[0.0, 1.0]`

**Parsing logic** (implemented in evaluation script):

```python
import re

def parse_llm_response(response: str) -> tuple[int, float]:
    """
    Returns (label, flight_probability).
    label: 1 = Flight (A), 0 = Ground (B)
    flight_probability: confidence if Flight, 1-confidence if Ground
    """
    # Match "A 0.82" or "B 0.35" anywhere in response
    pattern = r'\b([AB])\s+(0\.\d+|1\.0|0|1)\b'
    matches = re.findall(pattern, response.upper())
    if not matches:
        raise ValueError(f"Could not parse LLM response: {response[:100]}")
    letter, conf_str = matches[-1]   # take last match (final answer in CoT)
    conf = float(conf_str)
    label = 1 if letter == "A" else 0
    flight_prob = conf if label == 1 else (1.0 - conf)
    return label, flight_prob
```

**Flight probability** (for AUROC computation):
- If answer is `A` (Flight) with confidence `c` → `flight_prob = c`
- If answer is `B` (Ground) with confidence `c` → `flight_prob = 1 - c`

---

## Submission Format

Text LLM submissions use the **same JSON format** as Tier 1/2 submissions (see `docs/submission_format.md`), with additional metadata fields:

```json
{
  "task_id": "A4",
  "model_name": "GPT-4o_zeroshot",
  "model_description": "GPT-4o (gpt-4o-2024-11-20), zero-shot, top-50 genes",
  "tier": "3",
  "llm_variant": "zero_shot",
  "n_genes_in_prompt": 50,
  "submission_date": "2026-03-01",
  "predictions": {
    "fold_MHU-1_test": {
      "Mmus_C57-6J_TMS_MHU1_FLT_uG_Rep1": 0.82,
      "Mmus_C57-6J_TMS_MHU1_FLT_uG_Rep2": 0.78,
      "Mmus_C57-6J_TMS_MHU1_FLT_uG_Rep3": 0.91,
      "Mmus_C57-6CR_TMS_MHU1_GC_1G_Rep1": 0.21,
      "Mmus_C57-6CR_TMS_MHU1_GC_1G_Rep2": 0.15,
      "Mmus_C57-6CR_TMS_MHU1_GC_1G_Rep3": 0.09
    },
    "fold_MHU-2_test": { "...": "..." },
    "fold_RR-6_test":  { "...": "..." },
    "fold_RR-9_test":  { "...": "..." }
  }
}
```

**Additional fields for Tier 3**:

| Field | Type | Description |
|-------|------|-------------|
| `tier` | string | Must be `"3"` |
| `llm_variant` | string | `zero_shot` / `few_shot_k3` / `few_shot_k5` / `few_shot_k10` / `chain_of_thought` |
| `n_genes_in_prompt` | int | Number of genes per prompt (default: 50) |

---

## Generating Prompts from Task Files

```python
import pandas as pd
from pathlib import Path

def make_llm_prompt(
    fold_dir: Path,
    sample_id: str,
    tissue: str = "thymus",
    n_genes: int = 50,
    gene_symbols: dict = None,   # ENSMUSG → symbol mapping
) -> str:
    """Generate a zero-shot prompt for a single test sample."""
    selected_genes = (fold_dir / "selected_genes.txt").read_text().splitlines()
    selected_genes = selected_genes[:n_genes]  # top-N by variance (train-only order)

    gene_strs = []
    for ensembl in selected_genes:
        symbol = gene_symbols.get(ensembl, ensembl) if gene_symbols else ensembl
        gene_strs.append(f"{symbol} ({ensembl})")

    genes_formatted = ", ".join(gene_strs)

    return (
        f"A mouse {tissue} RNA-seq sample shows the highest expression variability "
        f"in the following {n_genes} genes (listed in descending order of variance "
        f"across the training cohort):\n\n"
        f"{genes_formatted}\n\n"
        f"Based on your knowledge of spaceflight biology, predict whether this sample is from:\n"
        f"(A) Flight — collected from mice aboard the International Space Station or similar "
        f"microgravity environment\n"
        f"(B) Ground — collected from ground control mice under normal gravity\n\n"
        f"Answer with the letter A or B, followed by a confidence score between "
        f"0.0 (uncertain) and 1.0 (certain).\n"
        f"Format: \"<A or B> <confidence>\"\n"
        f"Example: \"A 0.82\""
    )
```

**Gene symbol mapping**: Obtain from `data/mouse/ensembl_mouse_symbols.tsv` or via Ensembl BioMart query (same approach as Geneformer ortholog mapping in `scripts/geneformer_tokenize.py`).

---

## Target Models (Phase 2)

| Model | Provider | Tier | Variant |
|-------|----------|------|---------|
| GPT-4o (`gpt-4o-2024-11-20`) | OpenAI | 3 | Zero-shot, Few-shot k=3/5/10, CoT |
| Claude Opus (`claude-opus-4-6`) | Anthropic | 3 | Zero-shot, Few-shot k=3/5/10, CoT |
| Llama 3.1-70B | Meta (open) | 3 | Zero-shot, Few-shot k=3 |

---

## Evaluation

Text LLM submissions are evaluated identically to Tier 1/2:

```bash
python scripts/evaluate_submission.py \
    --submission GPT-4o_zeroshot_A4.json \
    --task A4
```

Output stored to:
```
evaluation/llm_baseline/{model_name}_A{n}_{tissue}.json
```

---

## Design Principles (DD-16)

1. **LOMO integrity**: Gene list derived from training variance only. LLM sees no test-set-specific patterns.
2. **No data leakage**: Training labels are excluded from zero-shot prompts. Few-shot examples use *training* samples only.
3. **Calibrated probabilities**: Parse confidence from LLM response → convert to `flight_probability` → AUROC evaluation (same pipeline as all tiers).
4. **Reproducibility**: Fix model version (API snapshot date), temperature=0, seed where supported.
5. **Comparability**: All three tiers evaluated on identical fold splits, same AUROC/CI/perm_p pipeline.
