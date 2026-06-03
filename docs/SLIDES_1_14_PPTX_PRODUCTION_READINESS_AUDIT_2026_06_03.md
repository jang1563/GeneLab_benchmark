# Slides 1-14 PPTX Production Readiness Audit

Date: 2026-06-03

## Scope

This audit checks whether slides 1-14 are ready to move from visual-asset production into actual PPTX assembly. It does not declare the deck final. It asks a narrower production question: can the current assets support a premium layered slide build without content drift, missing overlays, or hidden claim-boundary problems?

## Decision

Proceed to PPTX assembly as a production draft.

The required build rule remains: layered PNG scene plus editable scientific interpretation. The rendered previews should not be pasted as final all-in-one slides. Each slide should use its scene plate as the Z0-Z2 proof surface, then rebuild title, bridge text, source/scope line, caveat, and presenter-facing interpretation as editable PPTX objects.

## Audit Outputs

- Manifest: `output/slides_1_14_pptx_readiness_audit_v0_1/slides_1_14_pptx_readiness_manifest.json`
- JSON audit: `output/slides_1_14_pptx_readiness_audit_v0_1/slides_1_14_pptx_readiness_audit.json`
- CSV matrix: `output/slides_1_14_pptx_readiness_audit_v0_1/slides_1_14_pptx_readiness_matrix.csv`
- Summary: `output/slides_1_14_pptx_readiness_audit_v0_1/slides_1_14_pptx_readiness_summary.md`
- Status board: `output/slides_1_14_pptx_readiness_audit_v0_1/slides_1_14_pptx_readiness_status_board.png`
- Audit grayscale QA for early result slides 7-9: `output/slides_1_14_pptx_readiness_audit_v0_1/qa/`

## Automatic Gate Result

All 14 slides pass the core readiness gate.

- Rendered assets exist: pass
- Scene plates exist: pass
- Editable overlay or caption sources exist: pass
- Grayscale QA exists: pass
- Forbidden visible terms absent: pass
- Pre-slide-20 extension terms absent: pass
- Blockers: none

The audit also generated individual grayscale QA files for slides 7-9, because the early result family had source scene plates and layered scenes but no individual grayscale QA files at the slide level.

## Content Drift Check

The first 14 slides remain aligned with the original deck purpose: opening thesis, methods bridge, feature-layer explanation, early v1-v7 result logic, and core result spine.

No organoid, OSD-120, or v8/v9 extension evidence appears in the visible text scanned for slides 1-14. That material should stay out of the main deck until after slide 20 or until a later extension block is intentionally introduced.

## Slide Block Readiness

Slides 1-3 are ready as the opening block. They establish the benchmark thesis, external positioning, and project evolution. The PPTX build should keep firstness claims out and preserve the separation between completed results and v8/v9 extension/platform layers.

Slides 4-5 are ready as compressed dark methods slides. Slide 4 must carry the task-unit definition in notes or backup. Slide 5 must carry the train-only mission-held-out guardrail in notes or backup. Do not expand these into dense implementation diagrams unless the talk is explicitly methods-heavy.

Slide 6 is ready as the feature-layer bridge. It should explain the difference between gene-level and pathway-summary views before tissue/model results appear.

Slides 7-9 are ready as the first result family. Their purpose is to establish tissue-specific transfer, selected pathway rescue, and model-comparison limits before the hardened result spine.

Slides 10-14 are ready as the core result spine. Keep every claim boundary intact: hardening is coverage evidence, temporal/RRRM is guardrail evidence, negative results define current task limits, biological interpretation is hypothesis generation, and human translation is partial alignment only.

## Viewer-Language Watchpoints

These terms are acceptable in a technical talk, but they need a spoken bridge for first-time viewers:

- Slide 1: GeneLab, OSDR, held-out
- Slide 2: held-out, v1-v9
- Slide 3: v1-v9, v8/v9
- Slide 5: held-out
- Slide 6: held-out
- Slide 7: AUROC, held-out
- Slide 9: PCA-LR, foundation-model
- Slide 11: RRRM
- Slide 12: foundation-model

The recommendation is not to remove every technical term. Instead, define them once, keep the visible slide text plain, and put the heavier definitions in notes or backup.

## Required PPTX Build Rules

1. Use scene plates as visual proof surfaces, not as finished slides.
2. Rebuild all major slide text as editable objects.
3. Keep source/scope and caveat lines visible but subordinate.
4. Insert speaker notes for slides 4, 5, 6, and 10-14 before rehearsal.
5. Keep table-like evidence out of main figure cards; move it to backup or manuscript tables.
6. Keep v8/v9 extension visuals out of slides 1-14.
7. After PPTX assembly, rerender the full deck and repeat color, grayscale, and overlap QA on the actual exported slides.

## Residual Risks

Slide 10 has a tight word budget and should not gain extra label text during PPTX assembly.

Slides 7-9 are good main-figure candidates, but their final source/caveat lines should be rebuilt as explicit editable text rather than relying on tiny embedded labels.

Slide 13 should continue to use current JSON-derived values and should not reuse legacy HTML narrative statements that conflict with the current source summaries.

The readiness board is an audit artifact, not a slide design target.

## Next Anchor

Build the actual slides 1-14 PPTX skeleton using the audited scene plates, editable overlay text, and speaker-note requirements from this audit.
