# V9 Purpose Drift Audit

Date: 2026-05-26

Status: strategic alignment checkpoint updated after V9-REFOCUS-001

## Verdict

V9 has not lost its original purpose, but it has accumulated a moderate drift
risk.

The current OSD-120 chain is well aligned with the provenance-first and
claim-boundary parts of the original plan. It is not yet aligned with the
flagship-track timing in the original v9 design, because many consecutive
blocks have stayed inside one multispecies diagnostic/release-metadata branch
instead of returning to the public benchmark core or the first intended
single-cell flagship.

Decision:

- Keep the OSD-120 work as a bounded Phase 4A diagnostic/provenance case study.
- Treat V9-MULTI-035 as the completed closeout/retry note, not another
  expanding archive-metadata chain.
- Recenter on the public bulk alpha freeze path first, then return to the
  first single-cell flagship inventory after the bulk alpha gap matrix.

## Original V9 Contract

The original documents define v9 as a benchmark platform build, not a single
analysis sprint.

Core purpose:

- Public, living, provenance-first benchmark for biological AI under
  spaceflight domain shift.
- Test generalization across missions, tissues, species, modalities,
  perturbations, and stressor regimes.
- Every promoted claim links to source, code, output, metric, and validation
  manifests.

Primary order of operations:

- Build the platform spine first: package, task registry, metric profiles,
  manifests, dataset bundle.
- Choose a first flagship public benchmark, originally mission-held-out
  scRNA-seq/virtual-cell evaluation.
- Keep radiation/stressor tasks as the second flagship or v9.1 track.
- Add human organoid and multispecies extension only after schema support
  prevents conflating mouse tissue, human organoid, and non-mouse species.

Guardrails:

- Freeze claims before expanding scope.
- Use simple baselines before foundation-model adapters.
- Keep public benchmark tasks independent of gated human data.
- Keep intervention and Mars-regime outputs hypothesis-only unless separately
  validated.
- Maintain frozen snapshots and living refreshes as separate tracks.
- Add one flagship scientific track at a time.

## Current Alignment

Strong alignment:

- The package/API spine is real: manifests, registry, evaluator, submissions,
  run manifests, baseline runners, and tests exist.
- Public bulk LOMO tasks are still the stable core surface and remain separate
  from organoid/multispecies draft tracks.
- Organoid and multispecies work added schema fields and source boundaries
  before making baseline claims.
- OSD-120 outputs repeatedly preserve draft-only claim boundaries.
- The release-metadata chain prevents overclaiming: no DOI, no archive
  identifier, no frozen release version, no owner-inferred creator list, no
  license, and no descriptor mutation are claimed.
- Gated human data, intervention claims, and Mars-style projection claims remain
  outside the public benchmark core.

Partial alignment:

- OSD-120 is a useful Phase 4A example of non-mouse task construction,
  fragile-fold diagnostics, model transparency, and release-boundary discipline.
- The OSD-120 RO-Crate/Data Package scaffolds support the provenance-first plan,
  but they are not the same as finishing a public v9-alpha benchmark snapshot.
- The organoid response-signature and feature-effect work is scientifically
  useful, but it remains diagnostic rather than primary leaderboard work.

Misalignment risk:

- The original first flagship direction was mission-held-out scRNA-seq /
  virtual-cell evaluation. That track has not yet been re-entered.
- Phase 4A was supposed to come after the v9 alpha public bulk spine was stable.
  We have done substantial Phase 4A work before fully resolving public bulk
  payload freeze and release metadata.
- OSD-120 has become a long branch. It is still bounded, but one more
  open-ended metadata chain would start to feel like release paperwork taking
  the place of benchmark-platform progress.
- The number of draft diagnostic reports is high. Without a recentering memo,
  the narrative could become "many analyses" rather than "one benchmark
  platform with clear task families."

## Drift Risk Table

| Area | Status | Drift Risk | Decision |
|---|---|---:|---|
| Platform spine | Strongly aligned | Low | Keep using tests/manifests as the backbone. |
| Public bulk alpha | Aligned but incomplete | Medium | Revisit payload freeze, packaging, and baseline summary after OSD-120 closeout. |
| Single-cell flagship | Under-entered | Medium-high | Restore as a near-term candidate after V9-MULTI-035. |
| Organoid diagnostics | Aligned as diagnostic extension | Medium | Keep as secondary/draft-alpha, not primary benchmark narrative. |
| Multispecies OSD-120 | Aligned as Phase 4A case study | Medium | Stop expansion after diagnostic release note unless owner metadata arrives. |
| Release metadata | Aligned with provenance-first guardrail | Low-medium | Keep guard, but do not let archive metadata work displace task-family progress. |
| Claims/release language | Strongly aligned | Low | Continue explicit no-release/no-DOI/no-leaderboard boundaries. |
| Foundation-model adapters | Not over-expanded | Low | Still defer until simple baselines and reproducibility are strong. |

## Recenter Decision

V9-MULTI-035 is the OSD-120 release-metadata closeout block. It produced a
diagnostic metadata release note plus owner-metadata retry checklist, then
stopped the OSD-120 metadata chain unless new owner-supplied release metadata
appears.

After V9-MULTI-035, choose one of two recenter paths:

1. Public bulk alpha recenter:
   - Confirm current public bulk source/payload freeze gaps.
   - Decide whether a minimal alpha snapshot can be released without full
     payload mirroring.
   - Update dataset-card and package boundary around the public bulk core.

2. Single-cell flagship recenter:
   - Inventory RRRM/scRNA-seq assets.
   - Draft one AnnData-style task manifest.
   - Define the first `genelab_sc` metric contract before model adapters.

Preferred next path:

- V9-REFOCUS-001 selected public bulk alpha readiness first.
- V9-BULK-ALPHA-001 produced the public bulk alpha freeze-gap matrix.
- Run V9-BULK-ALPHA-002 to decide metadata-only alpha wording versus
  payload-mirror-first.
- Return to V9-SC-001 after the bulk alpha claim/payload boundary is explicit.

## Hard Guardrails From This Audit

- Do not promote OSD-120 to the v9-alpha benchmark core.
- Do not add another OSD-120 archive/release metadata block after
  V9-MULTI-035 unless owner-supplied metadata arrives.
- Do not generate DOI, CITATION.cff, license, release tag, publisher, or
  creator metadata from inference.
- Do not let organoid/multispecies diagnostics become the first flagship unless
  the plan is explicitly revised.
- Keep all current OSD-120 language diagnostic-only and not-leaderboard.

## Bottom Line

The project is still on mission: the recent work strengthens provenance,
reproducibility, and claim discipline. The risk is not scientific wrongness; it
is narrative and sequencing drift. The fix is to close the OSD-120 metadata
branch cleanly and then return to the benchmark-platform path.
