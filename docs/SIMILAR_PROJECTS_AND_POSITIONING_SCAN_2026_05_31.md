# Similar Projects And Positioning Scan

Date: 2026-05-31

Purpose: external positioning scan for a possible v1-v9 GeneLab Benchmark /
SpaceBio-Bench presentation and paper.

Companion query log:

- `docs/LITERATURE_QUERY_LOG_SPACEBIO_BENCH_2026_05_31.md`

## Correction After JK Note

`SpaceOmicsBench` and `SOMA` are JK/sibling projects, not external
competitors. They should be treated as the broader self/sibling ecosystem:

- `SOMA` = human spaceflight multi-omics atlas / biobank anchor.
- `SpaceOmicsBench` = human/civilian-astronaut multi-omics and LLM benchmark
  branch.
- `GeneLab Benchmark / SpaceBio-Bench` = GeneLab/OSDR mouse-first
  mission-held-out and provenance-first benchmark/platform branch.

This means the external competitive scan should focus on other NASA/OSDR AI
benchmark efforts, GeneLab ML pipelines, and general benchmark frameworks.

## Bottom Line

No clean one-to-one external duplicate was found in this scan.

The closest external neighbors are:

1. NASA BPS RNA Sequencing Benchmark Training Dataset / AI4LS benchmark work.
2. GLARE: GeneLab Representation learning pipelinE.
3. NASA AI4LS method/discovery papers and projects.
4. NASA OSDR/GeneLab infrastructure and AWG discovery papers.
5. OpenProblems / cell-eval / Virtual Cell Challenge as benchmark-framework
   precedents.

The strongest positioning remains:

> SpaceBio-Bench is a mission-held-out, provenance-first benchmark/platform for
> evaluating biological AI under GeneLab/OSDR spaceflight omics domain shift.

The key differentiator is not "first spaceflight AI benchmark." It is:

- mission-held-out generalization;
- multi-tissue mouse GeneLab/OSDR transcriptomics;
- pathway and feature abstraction under domain shift;
- classical ML vs FM/LLM/GNN comparison under small-n shift;
- provenance-first task manifests and release boundaries;
- extension lanes for human organoids, multispecies, and single-cell tasks.

## Self / Sibling Ecosystem

### SpaceOmicsBench

Source checked:

- https://huggingface.co/datasets/jang1563/SpaceOmicsBench

Relationship:

- Self/sibling project, not an external competitor.
- Human/civilian-astronaut spaceflight biomedical multi-omics and LLM benchmark
  branch.

How to use in this paper:

- Cite as part of the broader JK/Mason spaceflight-AI benchmark ecosystem.
- Use to show that the ecosystem already covers human multi-omics benchmark
  tasks, while this paper contributes the GeneLab/OSDR mission-held-out
  transcriptomics benchmark branch.
- Avoid broad "first spaceflight omics benchmark" claims because
  SpaceOmicsBench already occupies that broader phrase.

### SOMA

Sources checked:

- Nature 2024 SOMA paper:
  https://www.nature.com/articles/s41586-024-07639-y
- SOMA Browser:
  https://soma.weill.cornell.edu/apps/SOMA_Browser/

Relationship:

- Self/sibling ecosystem anchor, not an external competitor.
- Human spaceflight atlas/biobank and validation source.

How to use in this paper:

- Cite as the human-spaceflight atlas anchor.
- Use as motivation for cross-species and human-bridge tasks.
- Do not frame SpaceBio-Bench as competing with SOMA. The clean relationship is
  atlas/resource layer versus benchmark/evaluation layer.

## Closest External Neighbors

### 1. NASA BPS RNA Sequencing Benchmark Training Dataset / AI4LS

Sources checked:

- NASA AI4LS:
  https://www.nasa.gov/ames/space-biosciences/research-branch/artificial-intelligence-for-life-in-space/
- NTRS ASGSR 2023 poster:
  https://ntrs.nasa.gov/citations/20230015992
- AWS Open Data Registry NASA SMD AI tag:
  https://registry.opendata.aws/tag/NASA-SMD-AI/

What it is:

- NASA AI4LS describes standardized AI-ready RNA-seq and microscopy benchmark
  datasets for AI/ML benchmarking, training, and data challenges.
- The NTRS poster says the RNA-seq benchmark used space-flown mouse liver
  samples from NASA GeneLab and GAN augmentation to synthesize a larger
  benchmark/training dataset.

Overlap:

- NASA/GeneLab RNA-seq.
- AI/ML benchmark framing.
- Mouse spaceflight data.

Difference:

- Liver-centric.
- Synthetic-data and AI-ready training dataset orientation.
- Does not appear to define a multi-tissue leave-one-mission-out benchmark,
  pathway-abstraction comparison, multi-tier FM/LLM/GNN evaluation, or v9
  provenance-first task platform.

Positioning:

- Cite as the closest external "space biology AI-ready benchmark dataset"
  precedent.
- Contrast with this project's empirical mission-held-out, multi-tissue,
  multi-method evaluation.

### 2. GLARE: GeneLab Representation Learning Pipeline

Sources checked:

- Nature / npj Microgravity:
  https://www.nature.com/articles/s41526-025-00525-5
- PMC mirror:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC12569054/

What it is:

- An open-source pipeline for GeneLab/GLDS transcriptomic datasets.
- Uses projection, sparse autoencoders, self-supervised representation
  learning, ensemble clustering, and downstream biological interpretation.
- Demonstrated on CARA / OSD-120 and additional GLDS datasets OSD-217,
  OSD-406, and OSD-427.

Overlap:

- GeneLab/OSDR transcriptomics.
- ML applied to spaceflight omics.
- OSD-120 is also part of this repository's v9 multispecies branch.

Difference:

- GLARE is an analysis/discovery pipeline.
- It is not a fixed mission-held-out benchmark with standardized splits,
  metrics, model tiers, and provenance manifests.

Positioning:

- Cite as closest external ML-on-GeneLab analytical pipeline.
- Contrast "representation learning for hidden pattern discovery" with
  "mission-held-out benchmark for model generalization."

### 3. NASA AI4LS Method And Discovery Ecosystem

Sources checked:

- NASA AI4LS:
  https://www.nasa.gov/ames/space-biosciences/research-branch/artificial-intelligence-for-life-in-space/
- Guided Transfer Learning:
  https://arxiv.org/abs/2311.12045
- Causal inference ML ensemble:
  https://pubmed.ncbi.nlm.nih.gov/39824847/
- Age-response ML in murine mammary tissue:
  https://xbio.jmir.org/2026/1/e73041/

What it is:

- NASA AI4LS covers LLMs/foundation models, knowledge graphs, digital twins,
  causal inference, multimodal integration, federated learning, guided transfer
  learning, synthetic RNA-seq data, and benchmark datasets.
- Relevant papers include explainable multi-omics ML in space-flown mouse
  muscle, CRISP causal inference in space-flown mouse liver, guided transfer
  learning for small RNA-seq, and age-response ML in murine mammary tissue.

Overlap:

- Space biology AI/ML.
- GeneLab/OSDR-derived data.
- HDLSS RNA-seq modeling and small-sample ML methods.

Difference:

- Mostly method/discovery projects and data resources.
- Not a unified multi-tissue, mission-held-out benchmark/platform.

Positioning:

- Cite as the NASA space-biology AI/ML ecosystem.
- Frame SpaceBio-Bench as the standardized evaluation layer that can test such
  methods under cross-mission domain shift.

### 4. NASA OSDR / GeneLab

Sources checked:

- NASA GeneLab overview:
  https://www.nasa.gov/centers-and-facilities/ames/what-is-nasas-genelab/
- OSDR Biological Data API:
  https://visualization.osdr.nasa.gov/biodata/api/
- OSDR visualization portal:
  https://visualization.osdr.nasa.gov/

What it is:

- Public data repository and analysis ecosystem for space biology.
- Source spine for GeneLab Benchmark / SpaceBio-Bench.

Overlap:

- This project depends on OSDR/GeneLab.

Difference:

- OSDR/GeneLab is data infrastructure.
- It is not itself a model-evaluation benchmark with frozen splits and metrics.

Positioning:

- SpaceBio-Bench is a benchmark layer built on top of OSDR/GeneLab.

### 5. The Biology Of Spaceflight / GeneLab AWG Papers

Sources checked:

- NASA Cell Press package summary:
  https://www.nasa.gov/osdr-latest-news-nov-25-release-the-biology-of-spaceflight-a-package-of-29-scientific-papers-published-in-five-cell-press-journals/
- NASA AWG page:
  https://science.nasa.gov/biological-physical/data/awg

What it is:

- Large GeneLab/spaceflight biology discovery-paper ecosystem.
- Includes multi-omics analyses, radiation/countermeasure studies,
  tissue-specific biology, and GeneLab-derived secondary analyses.

Overlap:

- Spaceflight omics.
- GeneLab reuse and systems biology.

Difference:

- Discovery-paper corpus, not standardized AI benchmark infrastructure.

Positioning:

- Cite as biological motivation and proof that GeneLab secondary analysis is
  impactful.
- Present SpaceBio-Bench as a benchmark/evaluation complement to discovery
  papers.

### 6. X-Species GeneLab / TOAST / Cross-Species Explorer

Sources checked:

- NASA cross-species GeneLab article:
  https://www.nasa.gov/osdr-latest-news-multi-omics-cross-species-analysis-of-genelab-data-leads-to-new-nasa-investigation/
- Gilroy Lab X-Species GeneLab page:
  https://astrobiology.botany.wisc.edu/astrobotany-toast/x-species-genelab

What it is:

- Cross-species GeneLab exploratory tools.
- NASA describes transcriptional data across human, mouse, rat, fly, worm,
  yeast, rice, tomato, Arabidopsis, Brassica, and Brachypodium.

Overlap:

- Cross-species GeneLab analysis.
- Motivation for v9 multispecies and bridge tracks.

Difference:

- Exploratory visualization/data-mining ecosystem.
- Not a formal benchmark with fixed train/test splits and model comparisons.

Positioning:

- Cite as cross-species precedent.
- SpaceBio-Bench can convert cross-species exploration into auditable benchmark
  contracts.

## General Benchmark And Foundation-Model Evaluation Neighbors

### 7. OpenProblems

Sources checked:

- Nature Biotechnology 2025 paper:
  https://www.nature.com/articles/s41587-025-02694-w
- OpenProblems:
  https://openproblems.bio/

What it is:

- Community benchmark framework for single-cell analysis.

Relevance:

- Architectural precedent: tasks, datasets, metrics, methods, reproducible
  runs, and leaderboards.

Difference:

- Not spaceflight-specific.
- Does not address bulk GeneLab mission-held-out generalization.

### 8. Arc cell-eval / Virtual Cell Challenge / STATE

Sources checked:

- cell-eval:
  https://github.com/ArcInstitute/cell-eval
- VCC metrics:
  https://arcinstitute.org/news/behind-the-data-virtual-cell-challenge
- STATE:
  https://github.com/ArcInstitute/state

What it is:

- Perturbation-prediction model evaluation suite and challenge ecosystem.

Relevance:

- Useful for future single-cell and perturbation-like metric profiles.

Difference:

- Not spaceflight exposure generalization.
- Should not replace bulk LOMO metrics.

### 9. General Single-Cell Foundation Model Benchmarks

Sources checked:

- Nature Methods 2025 research highlight:
  https://www.nature.com/articles/s41592-025-02735-x
- PubMed zero-shot evaluation:
  https://pubmed.ncbi.nlm.nih.gov/40251685/

What it is:

- General single-cell foundation model evaluation literature.

Relevance:

- Supports the need for careful FM evaluation and the observation that FM
  performance is not automatically superior under distribution shift.

Difference:

- Not spaceflight-specific and not mission-held-out.

## Claim Guidance

Avoid:

- "First spaceflight omics benchmark."
- "First AI benchmark for space biology."
- "First use of GeneLab for ML."
- "First human spaceflight omics atlas."

Safer:

- "A mission-held-out benchmark for GeneLab/OSDR spaceflight
  transcriptomics."
- "A multi-tissue benchmark showing how biological AI models generalize across
  unseen spaceflight missions."
- "A provenance-first task-manifest platform for spaceflight omics model
  evaluation."
- "The GeneLab/OSDR mission-held-out branch of a broader JK/Mason
  spaceflight-AI ecosystem that includes SOMA and SpaceOmicsBench."

Best paper framing:

> Existing spaceflight omics resources and AI efforts include human multi-omics
> atlases, AI-ready RNA-seq benchmark datasets, representation-learning
> pipelines, and discovery-oriented ML studies. However, they do not center the
> specific problem of mission-held-out generalization across GeneLab/OSDR mouse
> tissues with standardized splits, pathway abstractions, model tiers, and
> provenance-backed task manifests. SpaceBio-Bench fills this gap.

## Recommended Citation Groups

Self/sibling ecosystem:

- SOMA.
- SpaceOmicsBench.
- Public GeneLab Benchmark dataset card, if treated as a companion artifact.

External space-biology AI/data context:

- NASA OSDR/GeneLab.
- NASA BPS RNA Sequencing Benchmark Training Dataset / AI4LS.
- GLARE.
- NASA AI4LS method papers: explainable multi-omics ML, causal inference,
  guided transfer learning, synthetic data, federated learning.
- The Biology of Spaceflight / AWG papers.
- X-Species GeneLab / TOAST.

Benchmark architecture context:

- OpenProblems.
- Arc cell-eval / Virtual Cell Challenge.
- General single-cell foundation model benchmark literature.

## Next Work Needed

Before manuscript submission:

1. Repeat the formal query log with PubMed export, Google Scholar export, and
   Semantic Scholar citation-neighborhood checks around NASA BPS, GLARE, CRISP,
   and Guided Transfer Learning.
2. Decide how to cite `SpaceOmicsBench` and `SOMA`: self/sibling prior work,
   companion ecosystem, or introduction context.
3. Decide how to cite the public `genelab-benchmark` artifact: data
   availability, companion dataset, or preliminary release.
4. Resolve the v4 canonical/raw table issue before presenting comparative
   method results.
