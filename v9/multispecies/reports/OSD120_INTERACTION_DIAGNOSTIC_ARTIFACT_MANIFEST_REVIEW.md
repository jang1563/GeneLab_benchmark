# OSD-120 Diagnostic Artifact Manifest Review

Status: V9-MULTI-025 complete as of 2026-05-26.

This block indexes the OSD-120 diagnostic evidence set as one traceable artifact
manifest. The goal is not to add new model results, but to make every local
claim auditable through files, hashes, row counts, and validation tests.

Primary outputs:

- `interaction_diagnostic_artifact_manifest/diagnostic_artifact_manifest.csv`
- `interaction_diagnostic_artifact_manifest/diagnostic_artifact_manifest.json`
- `interaction_diagnostic_artifact_manifest/diagnostic_claim_artifact_map.csv`
- `interaction_diagnostic_artifact_manifest/diagnostic_claim_artifact_map.json`

## Manifest Scope

The manifest indexes 26 artifacts across V9-MULTI-009 through V9-MULTI-025:

- OSD-120 task manifest and task-manifest index
- task-design review
- sparse L1 model summary, fold comparison, and feature audit
- sparse L1 stability summary and feature-frequency table
- baseline ladder summary and focus-fold table
- diagnostic candidate package
- figure/table package
- paired comparator package
- review notes
- validation test source

Each artifact row records file format, existence, byte size, line count, row
count when applicable, SHA-256, claim ids, validation scope, generator, and the
draft-only claim boundary.

## Claim Map

The claim map contains seven traceable claim rows:

- `draft_candidate_boundary`
- `nearest_centroid_fold_comparison`
- `fragile_focus_recovery`
- `feature_stability_evidence`
- `external_light_genotype_context`
- `figure_table_claim_boundary`
- `paired_comparator_decision`

Each claim row links back to artifact ids, artifact paths, validation tests,
external context URLs when relevant, and limitations.

## Decision

The OSD-120 diagnostic evidence set is now internally traceable. The current
draft still remains explicitly non-frozen: the manifest supports local
diagnostic interpretation, not a release snapshot, leaderboard baseline, or
leave-one-mission-out generalization claim.

Recommended next block: run an OSD-120 release/readiness gap audit that checks
which gaps remain before this diagnostic package could be promoted from
draft evidence to a cleaner public alpha artifact.
