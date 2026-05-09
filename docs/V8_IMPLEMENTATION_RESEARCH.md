# v8 Implementation Research Plan

Date: 2026-05-08

This note closes the v7 release surface and turns the v8 SpaceMed incubator into
an implementation plan that can be executed on HPC without relying on local
tests or laptop-scale computation.

## Current Decision

v7 should be closed as a public benchmark cleanup release. v8 should remain an
incubating analysis layer until every manuscript claim has a machine-readable
link to its input data, code command, output artifact, and validation status.

The local workstation is for static review, small metadata edits, and release
hygiene only. Regression tests, v8 summary regeneration, fGSEA, bootstrap, API
queries, and figure regeneration should run on HPC.

## v7 Closure Gate

Close v7 only after these checks are satisfied on HPC:

- `OSD-397` is the public eye mission label across README, task metadata,
  Hugging Face card, submission docs, and result summaries.
- `TBD` appears only as a legacy storage alias where old processed directories
  still require it.
- The fold rename is staged as a true rename:
  `tasks/A6_eye_lomo/fold_TBD_test` to
  `tasks/A6_eye_lomo/fold_OSD-397_test`.
- Portable paths are used for SpaceOmicsBench-dependent code through
  `SPACEOMICS_ROOT` or sibling discovery.
- No staged `.claude/`, `v8/bridge/geo_cache/`, `__pycache__/`, personal local
  paths, or single file larger than 50 MB.
- HPC validation command:
  `bash scripts/hpc_release_validate.sh --v8-summary`

Current static review shows the main remaining v7 closure action is not code
logic but staging discipline: the old `fold_TBD_test` deletions and the new
`fold_OSD-397_test` files must be staged together and inspected with
`git diff --cached -M --summary`.

## Research-Backed Implementation Decisions

### 1. Data access and provenance

NASA OSDR exposes biological data through REST and query interfaces, with JSON
and table outputs for metadata and data interrogation. v8 should therefore
record OSDR accessions, assay names, sample filters, file URLs, and checksums in
manifest files rather than relying on informal path notes.

SOMA is the right human-spaceflight anchor for BRIDGE because the 2024 Nature
SOMA paper describes an integrated astronaut multi-omics repository spanning
NASA Twins, JAXA CFE, Inspiration4, Axiom, and Polaris, with OSDR deposition and
SOMA browser access.

Implementation consequence:

- Add one manifest per pillar under `v8/provenance/runs/`.
- Require every table in `v8/RESULTS_SUMMARY.md` to point to one manifest.
- Keep `SPACEOMICS_ROOT` as the local pointer to SpaceOmicsBench and store the
  resolved dataset release hash or tag in the manifest.

### 2. Artifact split

GitHub documents repository health limits, recommends storing generated files
outside Git, and enforces a 100 MB single-object limit. Git LFS can hold larger
files by pointer, but per-file limits still apply by account plan. Hugging Face
Hub is suitable for public dataset-style artifacts, while Zenodo should archive
versioned release snapshots and mint DOIs. `CITATION.cff` should remain the
software citation entrypoint.

Implementation consequence:

- Git: code, small CSV/JSON summaries, signatures, manifests, final figures.
- Hugging Face dataset repo: benchmark tables, compact public result bundles,
  and regenerated pathway summaries.
- Zenodo: immutable release DOI for v7.1 and v8.0-beta snapshots.
- HPC/object storage: GEO cache downloads, SpaceOmicsBench raw inputs, and
  regenerable intermediate matrices.

### 3. FAIR and RO-Crate alignment

The FAIR principles require data and metadata to be findable, accessible,
interoperable, and reusable. RO-Crate provides a lightweight metadata packaging
pattern for research data, code, inputs, outputs, provenance, and authorship.
v8 does not need a full RO-Crate export immediately, but its manifest should be
RO-Crate-compatible.

Implementation consequence:

- Use `v8/provenance/manifest.schema.json` as the internal contract.
- Add a later `ro-crate-metadata.json` exporter once the beta result set is
  frozen.
- Treat missing provenance as a release blocker, not a documentation nicety.

### 4. Evaluation strategy

GeneLab v1-v7 should keep LOMO AUROC, macro-F1, balanced accuracy, and
per-tissue transfer as primary bulk-transcriptomics metrics. cell-eval/VCC
metrics are useful for v8 only where the biology matches: single-cell
perturbation prediction and DE recovery. Arc's VCC framing uses DES, PDS, and
MAE to capture DE recovery, perturbation distinguishability, and transcriptome
wide expression error.

Implementation consequence:

- Do not replace bulk LOMO metrics with PDS.
- Add a PDS-inspired `mission_discrimination` metric for embeddings or mission
  signatures.
- Use `cell-eval` directly only for future scRNA-seq/AnnData subtracks.
- Add DE overlap and direction-match metrics for perturbation signatures.

### 5. BRIDGE pillar

Current result surface is strong but needs reproducibility hardening:

- Mouse NES plus SpaceOmics features reports RF AUROC 0.888 vs 0.712 baseline
  after the 2026-05-09 HPC full-MSigDB rerun.
- Full MSigDB expansion increases the supervised lift and changes which tissues
  carry predictive information.
- The scientific claim should be phrased as pathway-level transfer, not raw
  gene-level human prediction.

Implementation tasks:

- BRIDGE manifests now cover `link_spaceomicsbench.py`,
  `tissue_nes_bridge.py`, `supervised_conservation.py`, and the supervised
  leakage audit.
- The 2026-05-09 HPC leakage audit records that the supervised model excludes
  the target label from all 14 model features, uses unique pathway merge keys,
  hashes deterministic stratified folds, and finds no near-perfect
  single-feature label proxy. Upstream SpaceOmicsBench feature-builder freezing
  remains a beta requirement.
- Save bootstrap seed, fold assignment, model hyperparameters, and feature list.
- Add a small `bridge/evaluation/README.md` mapping each result file to a claim.

### 6. DECOMPOSE pillar

The first-pass Mars extrapolation is scientifically useful as a regime-change
detector, but it should not be sold as a point predictor. Linear extrapolation
breaks at high dose/time amplification, so manuscript language must describe
Mars outputs as stressor-sensitivity flags until non-linear calibration exists.

Implementation tasks:

- Add manifests for every factorial analog dataset.
- Keep `raw_cache_audit.json` with the current HPC bundle so full-rerun
  readiness is explicit rather than inferred from failed wrapper attempts. The
  2026-05-09 HPC rerun restored OSD-211/237/202 processed caches and records
  all eight required count/sample-table files as present.
- Store the exact design formula, encoded factor levels, dose units, and
  contrast definitions.
- Separate low-LET and high-LET radiation quality claims.
- Add saturating or piecewise dose-response models as a v8-beta requirement,
  not a v8-alpha blocker.
- Report uncertainty propagation before any Mars-risk figure is promoted.

### 7. INTERVENE pillar

L1000CDS2 is appropriate for first-pass signature reversal because its API
supports reverse-mode gene-set and cosine searches over LINCS small molecule
signatures. Enrichr is appropriate for CRISPR KO orthogonal checks through its
submit-then-enrich API pattern. CLUE/CMap is valuable but access keys and data
files are individual, research-use resources and should not be redistributed.

Implementation tasks:

- Store request payload hashes, parsed-output hashes, query timestamp, gene
  symbol normalization rule, and top-N signature size. The current
  `api_snapshot_manifest.json` covers deterministic payloads and tracked parsed
  outputs without re-calling external APIs.
- Pin L1000CDS2 `db-version` in every query manifest instead of using untracked
  defaults; the current historical snapshot still records `latest`, so a
  concrete upstream version remains a beta-freeze item if the API exposes one.
- Keep API responses or large raw query dumps outside Git if they exceed the
  artifact policy.
- Label drug results as hypothesis-generating countermeasure candidates.
- Maintain `safety_triage.csv`: target class, known class-level toxicity
  liability, broad mechanism, tissue count, min reversal, and orthogonal
  support. This table is for pathway triage only, not clinical prioritization.

### 8. Causal/Digital Twin claims

The DAG and ICP layer should remain an integration scaffold until identifiability
and intervention assumptions are explicit. "Digital twin" is acceptable as a
project direction; manuscript claims should be careful unless the intervention
model has validated counterfactual predictions.

Implementation tasks:

- Add a DAG manifest with edge source, effect type, evidence file, and
  directionality assumption.
- Separate observational stability from causal intervention language.
- Promote do-calculus examples only after the assumptions are listed and
  falsification checks exist.

## Immediate v8 Build Order

1. Freeze v7.1 patch scope and stage the `OSD-397` rename as one coherent diff.
2. Add v8 provenance manifests for existing results without recomputing locally.
3. Add result-to-claim maps for BRIDGE, DECOMPOSE, INTERVENE, and Causal.
4. Add HPC entrypoints per pillar:
   `scripts/hpc_v8_bridge.sh`, `scripts/hpc_v8_decompose.sh`,
   `scripts/hpc_v8_intervene.sh`, `scripts/hpc_v8_causal.sh`, and
   `scripts/hpc_v8_summary.sh`.
5. Run the validation and computation only on HPC.
6. Rewrite v8 manuscript language where the current draft overstates clinical
   translation or point prediction.

## v8.0-alpha Definition of Done

- All current v8 scripts are portable through repo-relative paths or documented
  environment variables.
- Every existing v8 result file has a manifest stub.
- `v8/RESULTS_SUMMARY.md` has no empty key-result fields caused by schema drift.
- Large GEO caches and raw SpaceOmicsBench/SOMA files are not tracked.
- Drug and Mars claims are explicitly marked as exploratory unless supported by
  independent validation.

## v8.0-beta Definition of Done

- Full HPC recomputation from clean checkout.
- Frozen OSDR/SOMA/SpaceOmicsBench versions.
- One-command pillar regeneration on HPC.
- Manifest validation in CI or HPC gate.
- Final figures generated from manifest-linked result files.
- Hugging Face dataset bundle and Zenodo DOI prepared.

## Sources

- NASA OSDR Biological Data API:
  https://visualization.osdr.nasa.gov/biodata/api/
- NASA OSDR Developer API:
  https://www.nasa.gov/reference/osdr-developer-api/
- SOMA Nature paper:
  https://www.nature.com/articles/s41586-024-07639-y
- GitHub repository limits:
  https://docs.github.com/en/repositories/creating-and-managing-repositories/repository-limits
- GitHub Large File Storage:
  https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-git-large-file-storage
- Hugging Face Hub storage limits:
  https://huggingface.co/docs/hub/storage-limits
- Zenodo GitHub release archiving:
  https://help.zenodo.org/docs/github/archive-software/github-upload/
- GitHub CITATION files:
  https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files
- FAIR Guiding Principles:
  https://www.nature.com/articles/sdata201618
- RO-Crate:
  https://www.researchobject.org/ro-crate/
- Arc cell-eval:
  https://github.com/ArcInstitute/cell-eval
- Arc Virtual Cell Challenge metrics:
  https://arcinstitute.org/news/virtual-cell-challenge-2025-wrap-up
- L1000CDS2 API:
  https://maayanlab.cloud/L1000CDS2/help/
- Enrichr API protocol:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC8152575/
- CLUE developer resources and terms:
  https://clue.io/developer-resources
