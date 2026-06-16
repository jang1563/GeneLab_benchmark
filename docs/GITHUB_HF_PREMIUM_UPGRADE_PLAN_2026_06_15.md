# GitHub and Hugging Face Premium Upgrade Plan

Status: strategic execution plan
Date: 2026-06-15
Scope: GitHub, Hugging Face, citation metadata, release manifests, and public
documentation for GeneLab Benchmark / SpaceBio-Bench.

## Purpose

This plan upgrades the public project surface so that both human readers and
machines can understand the benchmark without reconstructing context from many
separate documents.

The goal is not simply a better README. The goal is a coherent public package:

- GitHub explains the project, methods, status, and release history.
- Hugging Face distributes benchmark artifacts with clear dataset-card metadata.
- Machine-readable manifests define versions, tasks, files, checksums, schemas,
  release status, and citation targets.
- Public-facing prose stays concise and confident, while detailed boundary and
  reproducibility information lives in cards, schemas, and manifests.

## Current Diagnosis

The repository already has strong raw material:

- `README.md` introduces the benchmark and summarizes v1-v7 results.
- `docs/hf_dataset_card.md` is a strong Hugging Face card for the public fold
  package.
- `docs/v9_hf_dataset_card.md` clearly labels the v9 public bulk track as a
  metadata-only alpha.
- `CITATION.cff` records the v7.1 software citation.
- `docs/SPACEBIOBENCH_TRANSPARENCY_CARD_PACK.md` maps the system, evaluation,
  release-readiness, and claim-register documents.
- `spacebio_bench/schemas/` already contains task, run, and metric-profile
  schemas.
- `spacebio_bench/datapackage.py` and `v9/datapackage.draft.json` begin the
  machine-readable data-package spine.

The main gap is not content volume. The main gap is a single public source of
truth. Today, README, Hugging Face card, citation metadata, v8 release notes,
and v9 alpha documents each carry part of the truth. A new reader or automated
agent has to infer which surface is current, which release lane is frozen, and
which artifact is authoritative.

## North Star

SpaceBio-Bench should feel like a premium benchmark platform:

- A scientist can quickly understand the benchmark question, data sources,
  task design, main results, and valid uses.
- A reviewer can trace every public claim to a release lane, result file,
  method description, and validation status.
- A machine agent can parse the release manifest and discover tasks, schemas,
  files, checksums, licenses, source accessions, and citation metadata.
- A contributor can run validation, rebuild public indices, and prepare a
  release without reading the entire documentation history.

## Public Status Model

Use one controlled status vocabulary across README, Hugging Face cards,
release manifests, and citation metadata.

| Lane | Public label | Meaning | Primary surfaces |
|---|---|---|---|
| v7.1 | canonical historical result surface | Public v1-v7 benchmark result and documentation surface | `README.md`, `docs/CANONICAL_RESULTS_V7_1.md`, `docs/hf_dataset_card.md`, `CITATION.cff` |
| v8 | incubating translational extension | Reproducibility hardening and artifact split in progress | `docs/V8_BETA_RELEASE_PLAN_2026_05_10.md`, `v8/release/` |
| v9 public bulk | metadata-only public bulk alpha | Task/source/provenance metadata and scaffold baselines are present; frozen payload release is not yet declared | `docs/v9_hf_dataset_card.md`, `v9/datapackage.draft.json`, `v9/task_manifest_index.csv` |
| v9 extension lanes | diagnostic or draft lanes | Single-cell, organoid, multispecies, and stressor extensions under design or diagnostic evaluation | `v9/README.md`, lane-specific reports |

Public pages should not expose too much internal decision language. Detailed
release boundaries should remain available, but the landing experience should
read as a clean product surface: what this is, what is available, how to use it,
and where deeper validation lives.

## Target Information Architecture

### GitHub

GitHub should be the canonical project operating system.

Recommended top-level public structure:

```text
README.md
CITATION.cff
LICENSE
CHANGELOG.md
CONTRIBUTING.md
docs/
  START_HERE.md
  RELEASE_INDEX.md
  cards/
  results/
  methods/
  releases/
  api/
spacebio_bench/
  schemas/
release/
  release_manifest.json
  release_manifest.schema.json
  checksums.sha256
examples/
scripts/
  release/
```

The root README should become a compact public landing page:

1. One-sentence benchmark identity.
2. Current public release table.
3. "Use the benchmark" quickstart.
4. "What is measured" task map.
5. "What the main results say" concise result panel.
6. Links to dataset card, transparency card pack, canonical results, v9 alpha,
   and citation.

Long result tables, historical detail, and lane-specific caveats should move
to linked docs.

### Hugging Face

Hugging Face should be the canonical dataset distribution surface.

Recommended dataset package shape:

```text
README.md
manifest.json
datapackage.json
checksums.sha256
task_manifest_index.csv
task_data_index.csv
source_inventory.csv
examples/
assets/
v7_public_folds/
v9_public_bulk_alpha/
```

The HF `README.md` should keep a strong YAML metadata block because Hugging Face
uses dataset-card metadata for discoverability and display. The prose should be
shorter than GitHub docs and optimized for dataset users:

1. What the package contains.
2. Which release lane it represents.
3. How to download/load it.
4. File layout.
5. Citation.
6. Links back to GitHub for full methods and transparency documentation.

If v7 public folds and v9 metadata alpha continue to share one HF repository,
the lane separation must be obvious in both paths and metadata. A separate
HF dataset repo for v9 is cleaner if public users are likely to confuse v7
benchmark folds with v9 alpha metadata.

## Machine-Readable Spine

Add a canonical `release/release_manifest.json` and validate it in CI.

Minimum fields:

```json
{
  "schema_version": "1.0.0",
  "project": {
    "name": "SpaceBio-Bench",
    "legacy_name": "GeneLab Benchmark",
    "repository": "https://github.com/jang1563/GeneLab_benchmark"
  },
  "release_lanes": [],
  "artifacts": [],
  "schemas": [],
  "huggingface": {},
  "citation": {},
  "checksums": {},
  "validation": {}
}
```

The manifest should connect:

- release lane and status,
- Git commit and optional tag,
- Hugging Face repo id and revision,
- Zenodo or DOI target when available,
- artifact paths and SHA-256 hashes,
- task manifests and schemas,
- source accessions and licenses,
- citation version and preferred citation.

This file should become the input for generated public surfaces:

- README release table,
- HF card release table,
- release index,
- citation consistency check,
- upload dry-run report,
- post-upload remote validation.

## Human-Facing Premium Assets

The public project should have a small number of polished, reusable visual
assets:

- benchmark task map,
- release-lane status map,
- artifact flow from OSDR source to task manifest to evaluation output,
- model-family comparison panel,
- HF package layout diagram.

These should be stored under `docs/assets/` and reused in README, HF card,
slide decks, and manuscript support material. Each asset should have a source
script and visual QA record.

## Automation Plan

Add release commands that make premium quality repeatable:

```bash
make release-qa
make build-public-manifest
make build-hf-package
make upload-hf-dry-run
make validate-hf-remote
```

Recommended checks:

- JSON Schema validation for release, task, run, and metric manifests.
- README, HF card, and `CITATION.cff` version/status consistency.
- Broken-link check for public docs.
- SHA-256 verification for release artifacts.
- HF dry-run diff before upload.
- HF remote smoke download after upload.
- Status-language lint so alpha/draft surfaces do not accidentally claim a
  frozen release.
- Large-file and local-cache guard before release commits.

## 30/60/90 Day Execution Board

### First 30 days: public coherence

- Create `release/release_manifest.schema.json`.
- Create first `release/release_manifest.json` covering v7.1, v8, and v9
  public bulk alpha.
- Add a generated `docs/RELEASE_INDEX.md`.
- Rewrite root README into a compact landing page.
- Keep historical result detail in `docs/CANONICAL_RESULTS_V7_1.md`.
- Add a README/HF/CITATION consistency test.

Deliverable: a new reader can understand the public state in under five
minutes, and an automated check can confirm that status labels agree.

### Days 31-60: HF package hardening

- Decide whether v9 alpha stays in the existing HF repo or moves to a separate
  dataset repo.
- Generate HF card sections from the release manifest.
- Add `manifest.json`, `datapackage.json`, and `checksums.sha256` to the HF
  package.
- Upgrade `scripts/upload_to_hf.py` into a release package builder with dry-run
  diff and post-upload validation.
- Add one polished HF summary graphic generated from source data.

Deliverable: HF becomes a deterministic, auditable distribution package rather
than a manually curated upload target.

### Days 61-90: CI and release discipline

- Add GitHub Actions for schema validation, docs consistency, and link checks.
- Add a release checklist template.
- Add contributor and issue templates.
- Add release-note generation from `release_manifest.json`.
- Record HF revision after every upload.
- Prepare Zenodo/DOI metadata if a frozen release is ready.

Deliverable: public releases can be prepared from a clean checkout with a
repeatable validation path.

### 3-6 months: platform maturity

- Add RO-Crate or DataCite-compatible metadata for archival release.
- Add richer examples for classical ML, foundation-model adapters, and
  submission validation.
- Add `docs/AI_AGENT_README.md` or `llms.txt` for machine-agent orientation.
- Publish stable API docs if `spacebio_bench` becomes an installable package.
- Add benchmark cards for each promoted task family.

Deliverable: SpaceBio-Bench reads as a professional benchmark platform, not a
single project repository.

## Decision Points

### Public name

Recommended approach:

- Use "SpaceBio-Bench" as the forward-looking platform name.
- Preserve "GeneLab Benchmark" as the historical v1-v7 name and repository
  alias.
- Make the transition explicit in README, HF card, and citation metadata.

### HF repository split

Option A: one HF repo.

- Pros: one public dataset link, simpler historical continuity.
- Cons: v7 public folds and v9 metadata alpha can be confused.

Option B: separate HF repos for stable/frozen and alpha/living tracks.

- Pros: cleaner user mental model and safer release language.
- Cons: more metadata maintenance.

Recommendation: keep the existing `jang1563/genelab-benchmark` for v7 public
folds and create a separate v9 alpha repo when public distribution expands
beyond metadata indices.

### Package identity

Decide whether the repo should remain documentation/scripts only or become an
installable Python package. If the latter, add `pyproject.toml`, console
entrypoints, and API docs.

## Quality Bar

A public update is premium quality only when all of the following are true:

- A human can identify the current public release state from the first screen
  of the README.
- The HF card says exactly what is packaged and what release lane it belongs to.
- `CITATION.cff` matches the release lane being cited.
- Every public artifact has a path, role, version, and checksum in a manifest.
- Every promoted task has a schema-valid task manifest.
- Every promoted result has a run manifest or equivalent provenance record.
- Public docs do not require internal project history to interpret.
- Release boundaries are clear without sounding apologetic or defensive.

## Immediate Next Actions

1. Create `release/release_manifest.schema.json`.
2. Create `release/release_manifest.json` with current v7.1, v8, and v9 public
   lane status.
3. Add a small validator script for release-manifest schema and required public
   files.
4. Draft the compact README structure without deleting existing result detail.
5. Draft the HF package layout decision: one repo now, separate v9 repo later.

## External Platform References

- Hugging Face Dataset Cards: https://huggingface.co/docs/hub/datasets-cards
- Hugging Face Model Cards: https://huggingface.co/docs/hub/model-cards
- GitHub README documentation: https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-readmes
- GitHub CITATION files: https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files
