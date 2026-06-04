---
title: SpaceBio-Bench Transparency Card Pack
page_type: transparency_card_pack
status: draft_public_ready
last_reviewed: 2026-06-04
claim_boundary: transparency_card_pack_no_new_release_claim
---

# SpaceBio-Bench Transparency Card Pack

This card pack is the public-facing map for SpaceBio-Bench scope, evaluation
interpretation, release readiness, and claim boundaries. It is meant to help a
reader understand what the benchmark currently supports before relying on a
headline score or release label.

This pack does not introduce new benchmark results or approve a new release.

## Current Public Summary

- v1-v7 / v7.1 is the canonical historical result surface for the original
  GeneLab Benchmark.
- v8 is an incubating translational extension and should not be mixed into
  v7.1 benchmark claims.
- v9 public bulk is a metadata-only alpha: task/source/provenance metadata and
  scaffold baselines are present, but payload-level hash verification remains
  pending.
- v9 single-cell, organoid, and multispecies lanes remain diagnostic or draft
  surfaces unless a lane-specific release boundary says otherwise.
- Current scores should be read as benchmark evidence, not biological mechanism
  proof, clinical guidance, countermeasure evidence, or mission-readiness
  evidence.

## Card Index

| Card | What it answers | File |
|---|---|---|
| Portfolio Brief | Why does this project matter as a research and portfolio artifact? | `docs/SPACEBIOBENCH_PORTFOLIO_BRIEF.md` |
| System Card | What is SpaceBio-Bench, what surfaces exist, and what is out of scope? | `docs/SPACEBIOBENCH_SYSTEM_CARD.md` |
| Evaluation Card | How should task, fold, baseline, metric, and pooled results be interpreted? | `docs/SPACEBIOBENCH_EVALUATION_CARD.md` |
| Release Readiness Card | Which release tier is each surface in, and what gates block stronger wording? | `docs/SPACEBIOBENCH_RELEASE_READINESS_CARD.md` |
| Claim Register | Which claims are supported, blocked, or future-only? | `docs/SPACEBIOBENCH_CLAIM_REGISTER.md` |
| v9 HF Dataset Card Draft | What does the v9 public bulk metadata alpha contain? | `docs/v9_hf_dataset_card.md` |
| v7.1 Canonical Results | What is the locked public result and scope source for v1-v7? | `docs/CANONICAL_RESULTS_V7_1.md` |

## Publication Surface Guidance

For GitHub:

- Link this pack from the root README and v9 README.
- Keep the four cards in `docs/`.
- Present the pack as transparency documentation, not as a new release.

For Hugging Face:

- Keep the dataset card concise.
- Add a short "Transparency and Release Boundary" section linking back to this
  GitHub card pack.
- Do not paste the full card pack into the Hugging Face README unless the full
  docs directory is also published there.

## Current Readiness Labels

| Surface | Recommended label | Stronger wording currently blocked |
|---|---|---|
| v1-v7 / v7.1 | canonical historical result surface | New v7.1 result generation |
| v8 SpaceMed | incubating translational extension | Countermeasure or intervention recommendation |
| v9 public bulk | metadata-only public bulk alpha | Frozen payload release; DOI/archive-ready release |
| v9 extension lanes | diagnostic or draft lane | Public leaderboard; frozen payload release |

## Minimum Public-Update Checklist

- Root README links to this card pack.
- Root README links to the portfolio brief.
- v9 README links to this card pack.
- v9 Hugging Face-style dataset card has a concise transparency section.
- `SPACEBIOBENCH_RELEASE_READINESS_CARD.md` is used as the release-readiness
  surface.
- No public-facing text claims a frozen v9 payload release, DOI/archive release,
  state-of-the-art leaderboard, clinical use, crew-health guidance,
  countermeasure recommendation, or Mars-regime prediction.
