# v9 Source and Competitor Matrix

Status: planning support for v9 design
Companion docs:

- `docs/V9_EXTERNAL_DEEP_RESEARCH_2026_05_21.md`
- `docs/V9_DESIGN_OPTIONS.md`

## Data and platform sources

| Source | URL | Class | Access status | v9 use | Main risk |
|---|---|---|---|---|---|
| NASA OSDR | https://science.nasa.gov/biological-physical/data/osdr/ | Space biology repository | Public | Source-of-truth for public datasets and accessions | API/data drift without frozen manifests |
| OSDR Biological Data API | https://visualization.osdr.nasa.gov/biodata/api/ | Metadata/data API | Public | Programmatic source manifests, dataset refresh, task construction | Query contract changes |
| SOMA Nature 2024 | https://www.nature.com/articles/s41586-024-07639-y | Human spaceflight multi-omics atlas | Public paper; data split across repositories | Human bridge and validation anchor | Small crew count and privacy constraints |
| SOMA Browser | https://soma.weill.cornell.edu/apps/SOMA_Browser/ | Human spaceflight browser | Public browser | Manual/source validation and exploratory bridge context | Not a benchmark API by itself |
| SOMA data/code collection | https://www.nature.com/collections/ebdbcahdgc/data-code | Linked data/code collection | Mixed public artifacts | Source traceability for human-spaceflight claims | Heterogeneous packaging |
| TRISH EXPAND | https://cdn.bcm.edu/academic-centers/space-medicine/expand | Commercial astronaut health database | Proposal/access-controlled | Future gated human track protocol | Access, IRB, consent, privacy |
| SpaceOmicsBench | https://huggingface.co/datasets/jang1563/SpaceOmicsBench | Processed benchmark dataset | Public HF dataset | Bridge/human summaries and packaging precedent | Upstream version must be pinned |
| L1000CDS2 | https://maayanlab.cloud/L1000CDS2/help/ | LINCS signature query tool | Public API/tool | Hypothesis-only reversal triage | API defaults and db-version ambiguity |
| Enrichr | https://maayanlab.cloud/Enrichr/ | Gene-set enrichment tool | Public API/tool | CRISPR/gene-set orthogonal support | Library versioning and API snapshots |
| CLUE/CMap | https://clue.io/ | LINCS/CMap resource | Account/research-use resources | Optional advanced perturbation resource | Licensing, redistribution, key management |

## Benchmark and model ecosystem

| Resource | URL | What it is | v9 lesson | Use directly? |
|---|---|---|---|---|
| OpenProblems | https://openproblems.bio/ | Continuous single-cell benchmark platform | Task registry, metrics, methods, leaderboard architecture | Architectural model, not direct dependency |
| OpenProblems paper | https://www.nature.com/articles/s41587-025-02694-w | 2025 benchmark framework paper | Community benchmark framing | Cite and emulate structure |
| Arc cell-eval | https://github.com/ArcInstitute/cell-eval | Perturbation-prediction metric package | Metric profiles and perturbation-specific scoring | Yes, only for scRNA-seq/AnnData tasks |
| Arc VCC metrics | https://arcinstitute.org/news/behind-the-data-virtual-cell-challenge | Explanation of VCC metrics | DES/PDS/MAE interpretation | Conceptual alignment |
| Arc STATE | https://arcinstitute.org/news/virtual-cell-model-state | Virtual-cell perturbation model | Model adapter target and field signal | Adapter only if reproducible |
| X-Cell | https://huggingface.co/Xaira-Therapeutics/X-Cell | Xaira virtual-cell model card | Emerging model target, licensing caveats | Not a hard dependency |
| scArchon | https://link.springer.com/article/10.1186/s13059-026-04104-z | 2026 virtual-cell benchmark paper | Ranking instability and metric sensitivity | Cite as warning and motivation |

## Standards and release infrastructure

| Standard/tool | URL | v9 role |
|---|---|---|
| FAIR principles | https://www.nature.com/articles/sdata201618 | Justification for findable, accessible, interoperable, reusable benchmark artifacts |
| RO-Crate | https://www.researchobject.org/ro-crate/ | Metadata export target for frozen releases |
| Hugging Face dataset cards | https://huggingface.co/docs/hub/datasets-cards | Public dataset documentation and artifact packaging |
| GitHub citation/Zenodo guidance | https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content | DOI and release archive planning |

## Competitive gap summary

| External ecosystem | What it does well | What it does not cover | v9 opportunity |
|---|---|---|---|
| OSDR/GeneLab | Public space biology data | Standard AI benchmark API and leaderboard | Convert curated tasks into runnable benchmark contracts |
| SOMA | Human astronaut multi-omics | Large-N model ranking | Human validation/bridge layer |
| TRISH EXPAND | Rich commercial astronaut data | Open public benchmark core | Gated protocol and schema readiness |
| OpenProblems | Benchmark architecture | Spaceflight domain shift | Reuse architecture ideas for space biology |
| cell-eval/VCC | Perturbation prediction metrics | Bulk mission-held-out tasks | Apply only to scRNA-seq subtracks |
| Virtual-cell models | Fast model innovation | Robust spaceflight OOD evaluation | SpaceBio-Bench as stress test |
| LINCS/Enrichr/CMap | Perturbation signature lookup | Validated countermeasures | Hypothesis-only intervention triage |

## Design decisions supported by this matrix

1. Public v9 cannot depend on restricted human data.
2. v9 should expose task contracts before adding model adapters.
3. cell-eval alignment should be scoped to single-cell tasks.
4. Radiation quality deserves its own benchmark family.
5. Intervention claims remain hypothesis triage.
6. Source and run manifests are part of the benchmark, not release paperwork.

