# RRRM-1 scRNA-seq Downstream Plan

Date: 2026-03-11
Scope: RRRM-1 mouse scRNA-seq (`OSD-918`, `OSD-920`, `OSD-924`, `OSD-934`)

## 1. Current recovery state

- `OSD-918`: STARsolo counting completed; `Solo.out/Gene` and `Solo.out/GeneFull` usable
- `OSD-920`: STARsolo counting completed; `Solo.out/Gene` and `Solo.out/GeneFull` usable
- `OSD-924`: successful completion
- `OSD-934`: rerun in progress after adding `--limitSjdbInsertNsj 1500000`

Interpretation:
- `918/920` do not need realignment unless BAM files are specifically required.
- Primary downstream input should be treated as STARsolo matrix output, not BAM.

## 2. Immediate goals

1. Freeze per-OSD count matrices
2. Convert each OSD to `.h5ad`
3. Run sample-level QC
4. Build tissue-level annotated objects
5. Decide whether cross-tissue integration is needed or whether tissue-specific analyses are sufficient

## 3. Input conventions

Preferred matrix source:
- `Solo.out/GeneFull/filtered`

Fallback:
- `Solo.out/Gene/filtered`

Rationale:
- `GeneFull` is generally preferable for 10x single-cell when intronic reads may be informative.
- Keep `Gene` outputs for sensitivity checks and compatibility comparisons.

Expected sample-to-tissue mapping:
- `OSD-918` → blood
- `OSD-920` → eye
- `OSD-924` → muscle
- `OSD-934` → skin

## 4. Blocking item

`rrrm1_h5ad_convert.py` currently depends on `scanpy`, which was not available in the default local or Cayuga shell Python tested during recovery.

Action:
- identify or create a Python environment with:
  - `scanpy`
  - `anndata`
  - `scipy`
  - `pandas`
  - `numpy`

If Python setup is slow, use R/Seurat first for QC and return to `.h5ad` export later.

## 5. Phase 1: data freeze and conversion

### 5.1 Freeze matrix outputs

For each OSD:
- verify presence of:
  - `barcodes.tsv`
  - `features.tsv`
  - `matrix.mtx`
- record matrix mode used:
  - `GeneFull/filtered`
  - `Gene/filtered`

Recommended manifest fields:
- `osd`
- `tissue`
- `status`
- `matrix_mode`
- `matrix_path`
- `raw_path`
- `log_path`
- `job_id`

### 5.2 Convert to `.h5ad`

Target outputs:
- `OSD-918_blood_rrrm1.h5ad`
- `OSD-920_eye_rrrm1.h5ad`
- `OSD-924_muscle_rrrm1.h5ad`
- `OSD-934_skin_rrrm1.h5ad`

Required obs metadata:
- `cell_barcode`
- `osd`
- `tissue`
- `study = RRRM-1`
- `mission = RR-8`
- `feature_mode = GeneFull`

Required var metadata:
- `gene_symbol`
- `mt`

## 6. Phase 2: sample-level QC

Use the current baseline filter first:
- `min_genes >= 200`
- `pct_counts_mt < 25`
- `min_cells per gene >= 3`

Generate for each sample:
- cell count before and after QC
- median UMI
- median genes
- mitochondrial fraction distribution
- top marker genes by abundance

Plots:
- violin: `n_genes_by_counts`, `total_counts`, `pct_counts_mt`
- scatter: `total_counts` vs `pct_counts_mt`
- scatter: `total_counts` vs `n_genes_by_counts`
- knee/barcode rank if useful from raw matrices

Decision point:
- if one tissue shows clearly shifted distributions, use tissue-specific QC thresholds instead of forcing one global rule.

## 7. Phase 3: doublets and ambient RNA

Minimal requirement:
- run one doublet detection pass per OSD

Preferred tools:
- Python: `scrublet`
- R: `DoubletFinder` if Seurat workflow is used first

Ambient RNA:
- quick check only at first pass
- compare known marker leakage across major lineages
- inspect raw vs filtered matrix scale

Do not block the first atlas on aggressive ambient correction unless obvious contamination is seen.

## 8. Phase 4: normalization and embeddings

Default first pass:
- library-size normalization
- log1p transform
- HVG selection per sample or per tissue
- PCA
- neighborhood graph
- UMAP
- Leiden clustering

Recommended analysis unit:
- tissue-specific objects first

Rationale:
- blood, eye, muscle, and skin are biologically distinct enough that a single global embedding will mostly recapitulate tissue identity.
- tissue-first analysis is lower risk and gives usable results faster.

## 9. Phase 5: annotation

### 9.1 Broad annotation targets

Blood:
- erythroid
- lymphoid
- myeloid
- platelet/megakaryocyte-like

Eye:
- epithelial
- stromal
- neuronal/glial
- immune

Muscle:
- myofiber-associated / myonuclei-related
- satellite/progenitor
- fibro-adipogenic / stromal
- endothelial
- immune

Skin:
- keratinocyte
- fibroblast
- endothelial
- immune
- appendage-associated populations if present

### 9.2 Annotation method

First pass:
- manual marker-based annotation

Second pass:
- reference transfer if a suitable mouse reference is available

Annotation outputs:
- cluster-to-celltype table
- marker table for each cluster
- confidence notes for ambiguous clusters

## 10. Phase 6: first biological readouts

Per tissue:
- total cells retained after QC
- cluster composition
- top marker genes by cluster
- proportion table

Cross-sample summary:
- one page summary of all four tissues
- whether cell recovery looks plausible
- obvious technical outliers

Potential first biological questions:
- blood: immune composition and erythroid signal
- muscle: stress/regeneration programs
- skin: keratinocyte vs stromal composition
- eye: major compartment identity and immune infiltration

## 11. Deliverables

Minimum deliverables for first downstream milestone:
- four per-OSD filtered matrix manifests
- four `.h5ad` objects
- four QC summaries
- four tissue-level UMAPs
- cluster marker tables
- annotation table
- brief downstream memo summarizing data quality and next biological questions

## 12. Suggested execution order

1. Confirm `OSD-934` rerun success
2. Build/locate Python environment with `scanpy`
3. Convert `918`, `920`, `924`, `934` to `.h5ad`
4. Run sample-level QC
5. Perform tissue-specific clustering and annotation
6. Generate first summary figures and tables

## 13. Nice-to-have improvements

- add a manifest-writing step to `rrrm1_h5ad_convert.py`
- optionally allow `--gene-mode GeneFull|Gene`
- emit a per-sample QC JSON/CSV alongside `.h5ad`
- preserve job/log provenance inside `adata.uns`
