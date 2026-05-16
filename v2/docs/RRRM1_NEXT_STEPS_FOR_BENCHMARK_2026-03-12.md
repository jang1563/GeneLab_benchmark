# RRRM-1 Next Steps for GeneLab Benchmark

Date: 2026-03-12
Scope: what is still worth doing for RRRM-1 scRNA-seq in the context of `GeneLab_benchmark`

## Current state

- Four benchmark-aligned tissues are now processed end-to-end:
  - `OSD-918` blood
  - `OSD-920` eye
  - `OSD-924` muscle
  - `OSD-934` skin
- Available outputs already cover the minimum benchmark-ready path:
  - `Solo.out/GeneFull/filtered`
  - per-sample `.h5ad`
  - merged `.h5ad`
  - tissue-specific processed `.h5ad`
  - broad annotation tables and summary proportions

## What counts as "enough" for GeneLab_benchmark

For this project, RRRM-1 does not need a publication-grade atlas before it becomes useful.
The practical threshold is:

- stable input matrices
- reproducible processing provenance
- coarse but defensible cell-type labels
- machine-readable manifests pointing to final benchmark inputs

That threshold has now been met.

## Recommended next steps

### Priority 1: Freeze benchmark metadata

- Treat [`RRRM1_BENCHMARK_READY_MANIFEST_2026-03-12.csv`](./RRRM1_BENCHMARK_READY_MANIFEST_2026-03-12.csv) as the dataset freeze record
- Update project-level summaries to mark RRRM-1 scRNA as benchmark-ready
- Keep the download, STARsolo, conversion, merge, exploratory, annotation, and summary scripts under version control

### Priority 2: Do only the single-cell work that materially reduces future rework

These are worth doing because they make the dataset cleaner and more defensible without changing the benchmark scope.

1. Tissue subtype refinement
- Review `*_cluster_markers_top20.csv` and `*_cluster_annotation.csv`
- Split obvious coarse labels where the evidence is strong
- Do not force subtype labels where marker support is weak

2. Doublet review
- Run one pass of doublet detection per tissue
- Flag suspicious clusters rather than aggressively deleting cells on the first pass

3. Ambient RNA review
- Check for obvious lineage leakage
- Prioritize blood and skin because they have strong dominant populations and may contaminate minority clusters

4. Keep integration modest
- Use tissue-specific objects first
- Avoid forcing a single cross-tissue biological integration unless the benchmark directly needs it

## Raw FASTQ cleanup review

### RRRM-1 selected tissues only

Current raw FASTQ footprint on Cayuga scratch:

- `OSD-918`: about `61G`
- `OSD-920`: about `69G`
- `OSD-924`: about `56G`
- `OSD-934`: about `122G`
- Total selected RRRM-1 raw FASTQ footprint: about `308G`

### Recommendation

- Zero-byte `Aligned.sortedByCoord.out.bam` placeholders can be deleted immediately
- Raw FASTQ directories for the four selected RRRM-1 tissues can be deleted from scratch **after**:
  - manifest and status docs are frozen
  - all final `.h5ad` and annotation outputs remain readable
  - the OSDR download script is preserved for re-download

### Rationale

- RRRM-1 benchmark outputs are already materialized in processed form
- raw FASTQs remain re-downloadable from OSDR if needed
- keeping 308G of selected raw FASTQs on scratch has low value once processed objects are frozen

### Caution

- Do not delete raw FASTQs if you expect to rerun alignment with different chemistry assumptions soon
- Do not delete unrelated `OSD-904` to `OSD-916` raw directories unless they are separately reviewed

## Practical recommendation

- Freeze metadata now
- Keep the current broad annotation as the benchmark baseline
- Perform targeted subtype/doublet/ambient cleanup next
- Then delete only the four selected RRRM-1 raw FASTQ directories from scratch if no rerun is planned
