# Literature Query Log For SpaceBio-Bench Positioning

Date: 2026-05-31

Purpose: formal query log supporting the positioning of a possible v1-v9
GeneLab Benchmark / SpaceBio-Bench presentation and manuscript.

This log follows the correction that `SOMA` and `SpaceOmicsBench` are
self/sibling projects, not external competitors. They are included only as
ecosystem anchors, not as external novelty threats.

## Search Scope

Search mode:

- Web search plus targeted source opening.
- Exact phrase searches were prioritized.
- PubMed/Nature/NASA/NTRS/AWS/GitHub/Arc/OpenProblems sources were preferred.

Limitations:

- This is not yet a complete systematic review.
- Google Scholar manual export was not performed in this pass.
- Citation counts were not audited.
- Search-engine ranking can drift; repeat this query log before submission.

## Inclusion And Exclusion Rules

Include as close external neighbor if the result has at least two of:

- spaceflight or space biology;
- GeneLab/OSDR/GLDS data;
- transcriptomics or omics;
- AI/ML/foundation model/representation learning;
- benchmark, standardized task, or evaluation framework.

Exclude or down-rank if:

- it is JK/self/sibling ecosystem work (`SOMA`, `SpaceOmicsBench`, public
  `genelab-benchmark` artifact);
- it is generic biomedical AI without spaceflight/GeneLab relevance;
- it is a data repository without benchmark/evaluation layer;
- it is general single-cell benchmarking without spaceflight relevance.

## Query Results

| Query | Key hits | Classification | Decision |
|---|---|---|---|
| `"spaceflight transcriptomics benchmark"` | public `genelab-benchmark`; NASA NTRS BPS benchmark; GLARE | self artifact + close external | Include NASA BPS and GLARE; treat `genelab-benchmark` as self artifact |
| `"GeneLab" "machine learning" "benchmark"` | NASA BPS RNA-seq benchmark; GLARE; AI4LS materials | close external | Include NASA BPS, GLARE, AI4LS |
| `"space biology AI benchmark"` | NASA NTRS benchmark; AI4LS; unrelated clinical/therapeutic benchmarks | mixed | Include NASA BPS/AI4LS; exclude generic biomedical AI |
| `"leave-one-mission-out" "spaceflight"` | mostly this project's public artifact / weak hits | self/low external | Supports novelty of LOMO phrasing |
| `"GeneLab" "representation learning" "GLARE"` | GLARE npj Microgravity / PMC | close external | Include as closest GeneLab ML pipeline |
| `"NASA BPS RNA Sequencing Benchmark Dataset"` | AWS Open Data Registry, NTRS poster, NASA AI4LS | close external | Include as closest benchmark-dataset precedent |
| `"Guided Transfer Learning" "space biology" RNA-seq` | arXiv GTL, NASA AI4LS | external method ecosystem | Include as method context, not benchmark duplicate |
| `"causal inference machine learning ensemble" "space-flown mice"` | Scientific Reports CRISP paper, PubMed/PMC | external method/discovery | Include as method ecosystem |
| `"OpenProblems" "Nature Biotechnology" 2025 benchmark single-cell` | OpenProblems Nature Biotechnology paper/site | benchmark architecture | Include as framework precedent |
| `"cell-eval" "Arc Institute" "Virtual Cell Challenge" metrics` | Arc cell-eval GitHub, VCC metrics article | benchmark architecture | Include for v9 single-cell metric profile context |
| `"spaceflight" "foundation model" transcriptomics GeneLab` | weak direct hits; AI4LS and general scFM literature | context only | Use to motivate careful FM evaluation, not close duplicate |

## Included External References

### NASA BPS RNA Sequencing Benchmark Training Dataset / AI4LS

Sources:

- NASA AI4LS:
  https://www.nasa.gov/ames/space-biosciences/research-branch/artificial-intelligence-for-life-in-space/
- NTRS ASGSR 2023 poster:
  https://ntrs.nasa.gov/citations/20230015992
- AWS Open Data Registry:
  https://registry.opendata.aws/bps_rnaseq/

Why include:

- Closest external benchmark-dataset precedent.
- Explicitly AI/ML-ready and space biology RNA-seq.
- Mouse liver data from NASA GeneLab plus synthetic/GAN augmentation.

How it differs from SpaceBio-Bench:

- Liver-centric, not multi-tissue.
- Synthetic training-data orientation, not empirical mission-held-out
  generalization.
- Does not provide v1-v7 style tissue hierarchy, pathway-rescue analysis,
  FM/LLM/GNN comparison, or v9 provenance-first task manifests.

Recommended citation role:

- "Prior NASA AI-ready benchmark dataset for space biology RNA-seq."

### GLARE

Sources:

- npj Microgravity:
  https://www.nature.com/articles/s41526-025-00525-5
- PMC:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC12569054/

Why include:

- Closest external GeneLab/GLDS machine-learning analysis pipeline.
- Uses representation learning and clustering to discover hidden patterns.
- Demonstrated on CARA / OSD-120 and other GLDS datasets.

How it differs from SpaceBio-Bench:

- Discovery pipeline, not benchmark platform.
- No fixed mission-held-out split/metric surface.
- Does not evaluate model tiers across tissues under LOMO domain shift.

Recommended citation role:

- "Prior ML representation-learning pipeline for GeneLab transcriptomics."

### NASA AI4LS Method/Discovery Ecosystem

Sources:

- NASA AI4LS overview:
  https://www.nasa.gov/ames/space-biosciences/research-branch/artificial-intelligence-for-life-in-space/
- Guided Transfer Learning:
  https://arxiv.org/abs/2311.12045
- Causal inference ML ensemble:
  https://pubmed.ncbi.nlm.nih.gov/39824847/
- Scientific Reports CRISP full article:
  https://www.nature.com/articles/s41598-024-81394-y
- Age-response ML paper:
  https://xbio.jmir.org/2026/1/e73041/

Why include:

- Shows NASA space-biology AI/ML ecosystem is active.
- Covers transfer learning, causal inference, knowledge graphs, synthetic data,
  federated learning, and spaceflight mouse transcriptomics.

How it differs from SpaceBio-Bench:

- Method/discovery orientation.
- Not a unified benchmark platform.
- Generally focused on specific tissues, phenotypes, or methods rather than
  multi-tissue mission-held-out model generalization.

Recommended citation role:

- "Relevant space-biology AI method ecosystem that motivates a standardized
  evaluation layer."

### NASA OSDR / GeneLab

Sources:

- NASA GeneLab:
  https://www.nasa.gov/centers-and-facilities/ames/what-is-nasas-genelab/
- OSDR Biological Data API:
  https://visualization.osdr.nasa.gov/biodata/api/

Why include:

- Source data infrastructure for this project.

How it differs from SpaceBio-Bench:

- Data repository and API, not benchmark/evaluation layer.

Recommended citation role:

- "Source data spine."

### OpenProblems

Sources:

- Nature Biotechnology 2025:
  https://www.nature.com/articles/s41587-025-02694-w
- OpenProblems:
  https://openproblems.bio/

Why include:

- Strong benchmark-platform precedent.
- Community-guided living benchmark architecture.

How it differs from SpaceBio-Bench:

- Single-cell analysis platform.
- Not spaceflight-specific.
- Does not target GeneLab/OSDR mission-held-out bulk transcriptomics.

Recommended citation role:

- "Benchmark architecture precedent."

### Arc cell-eval / Virtual Cell Challenge

Sources:

- cell-eval:
  https://github.com/ArcInstitute/cell-eval
- VCC metrics/data article:
  https://arcinstitute.org/news/behind-the-data-virtual-cell-challenge
- VCC wrap-up:
  https://arcinstitute.org/news/virtual-cell-challenge-2025-wrap-up

Why include:

- Perturbation-prediction metric and challenge ecosystem.
- Relevant to v9 single-cell/scaffold metric design.

How it differs from SpaceBio-Bench:

- Perturbation prediction in single-cell systems, not spaceflight exposure
  generalization.
- Should inform `genelab_sc` metric design but not replace bulk LOMO metrics.

Recommended citation role:

- "Single-cell perturbation metric precedent."

## Self / Sibling References

### SOMA

Sources:

- SOMA Nature 2024:
  https://www.nature.com/articles/s41586-024-07639-y
- SOMA Browser:
  https://soma.weill.cornell.edu/apps/SOMA_Browser/

Use:

- Human spaceflight atlas and validation anchor.
- Sibling ecosystem, not external competitor.

### SpaceOmicsBench

Source:

- https://huggingface.co/datasets/jang1563/SpaceOmicsBench

Use:

- Human/civilian-astronaut multi-omics and LLM benchmark branch.
- Sibling ecosystem, not external competitor.

### Public GeneLab Benchmark Artifact

Source:

- https://huggingface.co/datasets/jang1563/genelab-benchmark

Use:

- Current public release surface / companion artifact for this project.
- Not external prior work.

## Excluded Or Down-Ranked Hits

| Hit type | Reason |
|---|---|
| Generic therapeutic/clinical AI benchmarks such as CURE-Bench | Not spaceflight/GeneLab/OSDR-specific |
| Generic genomic sequence benchmarks | Not spaceflight transcriptomics or GeneLab benchmark |
| Space AI/geospatial foundation-model results | Not biological omics |
| General single-cell foundation-model papers without benchmark platform or spaceflight data | Useful background only |
| Search-result mirrors, blogs, and duplicate PDFs | Use primary paper/source instead |

## Novelty Assessment From Query Log

The external landscape supports these statements:

- There are AI-ready benchmark datasets for space biology RNA-seq, especially
  NASA BPS mouse liver/synthetic-data resources.
- There are GeneLab ML analysis pipelines such as GLARE.
- There are NASA AI4LS method/discovery papers using GeneLab/OSDR data.
- There are general benchmark frameworks for single-cell analysis and
  perturbation prediction.

The scan did not identify an external project that combines:

- multi-tissue GeneLab/OSDR mouse transcriptomics;
- leave-one-mission-out / mission-held-out evaluation as the core task;
- gene/pathway feature abstraction and transfer analysis;
- classical ML vs FM/LLM/GNN comparison under spaceflight domain shift;
- held-out validation and biological interpretation layers;
- provenance-first v9 task manifests and extension tracks.

Therefore, the paper should not claim broad firstness, but can credibly claim a
distinct mission-held-out GeneLab/OSDR benchmark/platform niche.

## Manuscript Language Draft

Safe introduction sentence:

> Prior space-biology AI efforts include NASA AI-ready RNA-seq benchmark
> datasets, GeneLab representation-learning pipelines, and method-specific
> studies in transfer learning, causal inference, and knowledge graphs. These
> resources motivate, but do not replace, a benchmark centered on
> mission-held-out generalization across GeneLab/OSDR tissues.

Safe novelty sentence:

> SpaceBio-Bench addresses this gap by defining mission-held-out
> transcriptomic tasks, standardized metrics, and provenance-backed task
> manifests for evaluating biological AI under spaceflight domain shift.

Unsafe sentence:

> This is the first AI benchmark for space biology.

## Next Search Pass

Before submission, repeat or extend with:

- PubMed exact query export.
- Google Scholar exact query export.
- Semantic Scholar citation-neighborhood scan around NASA BPS benchmark, GLARE,
  CRISP, and Guided Transfer Learning.
- Manual review of NTRS records linked to BPS RNA-seq and microscopy benchmark
  datasets.

