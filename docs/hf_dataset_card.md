---
license: cc-by-4.0
task_categories:
  - tabular-classification
tags:
  - genomics
  - transcriptomics
  - spaceflight
  - space-biology
  - benchmarking
  - nasa-osdr
  - bulk-rna-seq
  - mouse
  - rna-seq
  - mission-held-out
  - foundation-model
  - pathway-analysis
  - spacebio-bench
size_categories:
  - 100M<n<1GB
language:
  - en
pretty_name: "SpaceBio-Bench / GeneLab Benchmark v7.1.2 Public Fold Package"
viewer: false
---

# SpaceBio-Bench / GeneLab Benchmark v7.1.2 Public Fold Package

**Processed mission-held-out transcriptomics folds for evaluating whether
machine-learning and foundation-model methods generalize spaceflight biological
signatures across missions.**

Public status: **v7.1.2 public-card/metadata patch over canonical v7.1 results**

Dataset freeze: **2026-03-01**

Patch scope: documentation, public metadata, and access guidance. It does not
introduce new benchmark result generation.

Code and full documentation: <https://github.com/jang1563/GeneLab_benchmark>

Maintainer / citation author: JangKeun Kim, Weill Cornell Medicine.

![SpaceBio-Bench benchmark at a glance](assets/hf_benchmark_summary.png)

## What Is In This Dataset

This Hugging Face dataset contains self-contained public fold packages from the
GeneLab Benchmark v1-v7 surface. Each fold holds out one mission as the test set
and provides all files needed to train on the remaining missions and evaluate on
the held-out mission.

| Public package item | Description |
|---|---|
| `train_X.csv`, `test_X.csv` | Sample-by-gene expression matrices |
| `train_y.csv`, `test_y.csv` | Binary labels: `1 = Flight`, `0 = Ground` |
| `train_meta.csv`, `test_meta.csv` | Sample-level metadata used for fold auditing |
| `fold_info.json` | Held-out mission, train missions, and sample-count audit metadata |
| `selected_genes.txt` | Fold-specific genes selected from training missions only |
| `task_info.json` | Task-level metadata and source summary |

The web Dataset Viewer is disabled because these are high-dimensional
sample-by-gene matrices plus JSON artifacts. Use direct downloads for reliable
access.

## Public Fold Layout

```text
genelab-benchmark/
├── A2_gastrocnemius_lomo/
│   ├── task_info.json
│   ├── fold_RR-1_test/
│   ├── fold_RR-5_test/
│   └── fold_RR-9_test/
├── A4_thymus_lomo/
│   ├── task_info.json
│   ├── fold_MHU-1_test/
│   ├── fold_MHU-2_test/
│   ├── fold_RR-6_test/
│   ├── fold_RR-9_test/
│   └── fold_RR-23_holdout/
├── A5_skin_lomo/
│   ├── task_info.json
│   ├── fold_MHU-2_test/
│   ├── fold_RR-6_test/
│   ├── fold_RR-7_test/
│   └── fold_RR-7_holdout/
├── A6_eye_lomo/
│   ├── task_info.json
│   ├── fold_RR-1_test/
│   ├── fold_RR-3_test/
│   └── fold_OSD-397_test/
├── v4/evaluation/
├── v5/evaluation/
└── v6/evaluation/
```

`fold_OSD-397_test` is the stable public label for the third A6 eye fold.

## Scope

| Dimension | Coverage |
|---|---|
| Full v1-v7 benchmark surface | 8 tissues |
| Public source catalog | 24+ NASA OSDR accessions |
| Processed sample scope | 600+ binary/control samples across release layers |
| v4 multi-method evaluation | 8 tissues x 8 classifiers x 4 feature types = 256 evaluations |
| Public HF fold package | 4 reviewer-facing LOMO tasks plus selected result artifacts |

The full GitHub benchmark also includes historical v2-v7 extensions for
temporal dynamics, cross-species analysis, single-cell and spatial pilots,
foundation-model comparisons, graph/network baselines, and biological
interpretation layers. This HF repository is optimized for processed dataset
access; GitHub is the complete methods, code, and release-documentation surface.

## Download Example

```python
from huggingface_hub import hf_hub_download
import pandas as pd

repo_id = "jang1563/genelab-benchmark"
fold = "A5_skin_lomo/fold_RR-7_test"

def hf_csv(name):
    return pd.read_csv(
        hf_hub_download(
            repo_id=repo_id,
            filename=f"{fold}/{name}",
            repo_type="dataset",
        ),
        index_col=0,
    )

train_X = hf_csv("train_X.csv")
train_y = hf_csv("train_y.csv").iloc[:, 0]
test_X = hf_csv("test_X.csv")
test_y = hf_csv("test_y.csv").iloc[:, 0]
test_meta = hf_csv("test_meta.csv")

print(train_X.shape, train_y.shape, test_X.shape, test_y.shape)
```

Download a complete task:

```python
from huggingface_hub import snapshot_download

snapshot_download(
    repo_id="jang1563/genelab-benchmark",
    repo_type="dataset",
    allow_patterns="A5_skin_lomo/**",
    local_dir="./data/genelab-benchmark",
)
```

## File Contract

| File | Contract |
|---|---|
| Feature matrices | Rows are sample IDs; columns are Ensembl mouse gene IDs |
| Expression values | Log2(DESeq2 size-factor normalized counts + 1) |
| Gene selection | Top 75th percentile variance, computed on training missions only |
| Labels | Binary Flight/Ground labels |
| Metadata | Sample and fold metadata for auditability |

The fold design prevents test-mission leakage by applying variance filtering
inside each training split.

## Evaluation Summary

The canonical result source is `docs/CANONICAL_RESULTS_V7_1.md` in the GitHub
repository.

| Result surface | Takeaway |
|---|---|
| Multi-method benchmark | PCA-LR is the strongest 8-tissue gene-level baseline in v4, with mean AUROC 0.776. |
| Best tissue rows | Thymus 0.948, colon 0.921, lung 0.901, kidney 0.829 across best method-feature combinations. |
| Cross-mission transfer | Thymus and gastrocnemius show the strongest mission-transfer signal; liver and kidney are harder. |
| Pathway features | Pathway representations rescue some weaker gene-level tissues, especially kidney and eye. |
| Foundation models | Tested gene-expression foundation models underperform tuned classical baselines on small-n bulk RNA-seq mission shift. |
| Held-out validation | Thymus RR-23 AUROC 0.905; skin RR-7 AUROC 0.885. |

## Intended Use

Use this dataset to:

- evaluate spaceflight transcriptomics classifiers under mission-held-out shift;
- compare classical ML, foundation-model, and adapter methods on fixed folds;
- test preprocessing or feature representations without changing test missions;
- reproduce public benchmark summaries from the GitHub repository.

For full methods and release status, use the GitHub documentation and release
manifest.

## Release Labels

| Surface | Public label |
|---|---|
| v7.1 GeneLab Benchmark | Canonical historical result surface and citation target |
| v7.1.2 public-card patch | Documentation and metadata patch over v7.1 results |

This HF dataset card describes the v7.1 public fold package with the v7.1.2
public-card patch.

## Citation

Please cite the software and benchmark using the GitHub `CITATION.cff` metadata.

```bibtex
@dataset{kim2026genelab,
  title = {SpaceBio-Bench / GeneLab Benchmark: Mission-Held-Out Spaceflight Transcriptomics Benchmark},
  author = {Kim, JangKeun},
  year = {2026},
  url = {https://huggingface.co/datasets/jang1563/genelab-benchmark},
  note = {v7.1.2 documentation, public-card, and metadata patch over canonical v7.1 results; data freeze 2026-03-01}
}
```

Source data: NASA Open Science Data Repository (OSDR), <https://osdr.nasa.gov/bio/repo/>.

## License

- Processed dataset package: CC-BY-4.0
- Code: MIT, in the GitHub repository
- Source data: NASA OSDR public data; follow individual source-dataset terms
