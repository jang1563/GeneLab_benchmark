# GeneLabBench v3 HPC / Blocker Runbook

Last updated: 2026-03-20

## 1. UCE

### Goal
- Run UCE embeddings on GeneLabBench tissues on Cayuga without failing during env/bootstrap.

### Preconditions
- Project synced to: `/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark`
- Conda env exists:
  `/home/fs01/jak4013/miniconda3/miniconda3/envs/uce_env`
- `v3/logs/` exists before `sbatch` submission
- `uce-eval-single-anndata` is on `PATH` inside `uce_env`

### Batch script
- Script: `v3/scripts/slurm_uce_eval.sh`
- Default args: `--all --nlayers 4 --batch-size 25`
- Override path:
  `UCE_ARGS="--tissue liver --fast --embeddings-only" /opt/ohpc/pub/software/slurm/24.05.2/bin/sbatch v3/scripts/slurm_uce_eval.sh`

### Smoke test sequence
1. Sync patched scripts to Cayuga.
2. Confirm env activation from user Miniconda init, not system Anaconda.
3. Submit smoke job with:
   `UCE_ARGS="--tissue liver --fast --embeddings-only"`
4. Check `v3/logs/uce_eval_<JOBID>.out` and `.err`.
5. Only after smoke success, submit full job with default args.

### Success criteria
- Job stays alive past startup.
- `.err` does not contain env/CLI bootstrap failures.
- `FM_uce.json` or tissue embeddings are written under `v3/evaluation/`.

### Failure triage
- `EnvironmentNameNotFound`:
  source wrong conda init; use user Miniconda path.
- `uce-eval-single-anndata` missing:
  env incomplete; reinstall UCE package or fix PATH.
- Output `.h5ad` missing:
  CLI ran but did not emit embeddings; inspect CLI args/package version.

## 2. scFoundation

### Goal
- Run scFoundation evaluator on Cayuga with the real upstream directory layout and installed model weights.

### Preconditions
- Repo exists at:
  `/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/scFoundation`
- Env exists:
  `/home/fs01/jak4013/miniconda3/miniconda3/envs/scfoundation_env`
- Required runtime files:
  `scFoundation/model/get_embedding.py`
  `scFoundation/model/OS_scRNA_gene_index.19264.tsv`
  `scFoundation/model/models/models.ckpt`

### Batch script
- Script: `v3/scripts/slurm_scfoundation_eval.sh`
- Default args: `--tissue liver`
- Override example:
  `SCF_ARGS="--tissue liver --fast" /opt/ohpc/pub/software/slurm/24.05.2/bin/sbatch v3/scripts/slurm_scfoundation_eval.sh`

### Smoke test sequence
1. Sync patched evaluator + SLURM script to Cayuga.
2. Confirm `models/models.ckpt` exists.
3. Submit smoke job with:
   `SCF_ARGS="--tissue liver --fast"`
4. Inspect `v3/logs/scf_eval_<JOBID>.out` and `.err`.
5. If smoke succeeds, expand to additional tissues.

### Success criteria
- No `get_embedding.py` path error.
- No missing-weight error.
- `FM_scfoundation.json` is written under `v3/evaluation/`.

### Failure triage
- `get_embedding.py` not found:
  evaluator pointed at repo root only; must resolve nested `model/` runtime dir.
- `models.ckpt` missing:
  weights not downloaded yet; use SharePoint link in `model/models/download.txt`.
  As of 2026-03-20, the public link still requires authentication from our side.
- CUDA/runtime import failure:
  fix `scfoundation_env` packages before resubmission.

## 3. F3e (GLDS-270 Heart Visium)

### Goal
- Decide quickly whether F3e can move from blocked to exploratory execution.

### Minimum required inputs
- One of:
  processed Visium matrix outputs
  `filtered_feature_bc_matrix.h5` plus `spatial/` metadata
  a complete raw FASTQ + image + reference bundle for Space Ranger

### Decision path
1. Search Cayuga/local cache for processed GLDS-270 spatial outputs.
2. If found, stage a lightweight exploratory pipeline.
3. If only raw FASTQs exist, estimate preprocessing cost and decide explicitly whether to proceed.
4. If neither processed outputs nor preprocessing inputs are complete, keep blocked.

### Exit criteria
- `Go`: processed spatial inputs available.
- `No-Go`: raw FASTQ only with no preprocessing commitment.

## 4. I1 (Gene + Body/Organ Weight)

### Goal
- Determine whether multimodal integration is actually feasible with available metadata.

### Minimum required table schema
- `sample_id` or stable join key
- `animal_id` if sample-level join is indirect
- `mission`
- `label`
- `body_weight`
- at least one organ-weight column

### Decision path
1. Search OSDR/local metadata for phenotype tables.
2. Attempt one join for a single tissue, preferably liver.
3. Count matched samples per class.
4. If matched set is too small or ambiguous, keep deferred.

### Exit criteria
- `Go`: stable join + enough matched samples for a pilot model.
- `No-Go`: no stable join key or no usable phenotype matrix.

## 5. Recommended Operating Order

1. Patch locally.
2. Sync `v3/scripts/` and `v3/docs/` to Cayuga.
3. Run UCE smoke job.
4. Run scFoundation smoke job.
5. If both smoke tests pass, launch full runs.
6. Sync logs and evaluation JSONs back to Dropbox.
7. Reassess F3e and I1 with explicit go/no-go outcomes.
