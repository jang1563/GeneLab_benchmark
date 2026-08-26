# SpaceBio-Bench Conference Deck v0.5 Rehearsal Timing

Date: 2026-06-03

## Scope

This rehearsal pass turns the 24-slide v0.5 conference deck into a timed talk
track. It does not change the PPTX. The goal is to stabilize pacing before an
institutional-template or branding pass.

Primary deck:

`output/spacebiobench_conference_deck_v0_5/spacebiobench_conference_deck_v0_5.pptx`

## Speaker Thesis

SpaceBio-Bench is not a firstness claim and not a clinical translation claim.
It is a mission-held-out benchmark for testing whether biological AI signals
survive spaceflight domain shift, with a strict separation between completed
benchmark evidence, hypothesis-only translation, platform readiness, and blocked
extension tracks.

## 20-Minute Timing Plan

Target runtime: 19:45, leaving 15 seconds of buffer.

| Slide | Time | Cumulative | Purpose | Talk Track |
| --- | ---: | ---: | --- | --- |
| 1 | 0:40 | 0:40 | Thesis | "The question is not whether AI can describe biology in-distribution. The question is whether a model trained on known space biology missions can survive a new mission context." |
| 2 | 0:45 | 1:25 | Positioning | "This is a distinct niche, not a firstness claim. Existing resources matter, but the gap here is mission-held-out biological evaluation under spaceflight domain shift." |
| 3 | 0:50 | 2:15 | Story map | "The talk has three layers: completed benchmark results, hypothesis-only translation, and platformization. I will keep those layers separate throughout." |
| 4 | 1:05 | 3:20 | Data contract | "A benchmark task starts as a public record, but it becomes useful only when study, mission, sample, matrix, split, feature view, and metric stay attached as a contract." |
| 5 | 0:55 | 4:15 | Held-out protocol | "Held-out here means an entire mission is hidden. Preprocessing, feature choices, and model selection are made on train-side missions before the hidden mission is scored." |
| 6 | 1:00 | 5:15 | Feature and metric primer | "The same task can be scored from genes, pathways, or transformed model inputs. AUROC is the score grammar: 0.5 is chance, the dot is the score, and the interval is uncertainty." |
| 7 | 0:55 | 6:10 | First result | "The first result is tissue-specific transfer. Some tissues show stronger held-out separation, while others stay near chance. The claim is tissue-specific, not universal." |
| 8 | 0:45 | 6:55 | Pathway rescue | "Pathway summaries can suppress unwanted labels in selected settings. This supports a selected rescue claim, not universal pathway superiority." |
| 9 | 0:45 | 7:40 | Model comparison | "Scale alone does not solve transfer. In these tested settings, larger or more complex representations do not automatically beat simpler baselines." |
| 10 | 0:55 | 8:35 | Hardening | "The hardening grid reduces cherry-pick risk: 8 tissues, 8 classifiers, and 4 feature views. This is coverage evidence, not a new leaderboard." |
| 11 | 0:55 | 9:30 | Guardrails | "Timepoint and RRRM labels need guardrails. Preservation can dominate, recovery is descriptive, and the RRRM single-cell pilot is underpowered." |
| 12 | 0:45 | 10:15 | Negative boundary | "Negative results are useful because they define task limits. They do not prove absence of biology; they tell us where the current benchmark evidence stops." |
| 13 | 0:55 | 11:10 | Biology interpretation | "Benchmark hits become hypotheses. They can prioritize pathway or target follow-up, but they are not mechanism proof and not treatment evidence." |
| 14 | 0:55 | 12:05 | Human translation | "Human translation is partial, not direct. We can see some pathway or target-tier alignment, but weak direct gene transfer prevents a clinical claim." |
| 15 | 0:40 | 12:45 | v7 boundary | "v7.1 is a canonical result surface. It freezes documentation and claim discipline; it is not a new benchmark run." |
| 16 | 0:50 | 13:35 | v8 transition | "v8 uses benchmark signals as a hypothesis incubator. That is useful, but it remains downstream of the completed benchmark evidence." |
| 17 | 0:45 | 14:20 | Countermeasure boundary | "Stressors are not countermeasure recommendations. Provenance machinery can make the evidence traceable, but it does not make the biology operational." |
| 18 | 0:55 | 15:15 | Platform readiness | "The platform terms matter because metadata readiness is not payload readiness. Release claims wait for payload mirroring, hash verification, and evaluator runs." |
| 19 | 0:40 | 15:55 | v9 platform | "v9 platformizes the benchmark with manifests, evaluator outputs, and run records. This is reproducibility infrastructure, not a new biological claim." |
| 20 | 0:45 | 16:40 | Public bulk alpha | "Public bulk alpha is useful because 22 source records are parsed, but 0 payloads are locally hash-verified. That is metadata alpha, not frozen release evidence." |
| 21 | 0:35 | 17:15 | Organoid extension | "The organoid extension shows source records becoming local matrices, but it remains a draft pilot track outside the public bulk core." |
| 22 | 0:35 | 17:50 | OSD-120 | "OSD-120 is a same-study diagnostic, not the flagship mission-held-out transfer task. Compactness is not generalization." |
| 23 | 0:45 | 18:35 | Single-cell blocker | "The single-cell metric spec exists, but the processed payload is missing. Raw FASTQ is not enough for a leaderboard score." |
| 24 | 1:10 | 19:45 | Close | "The closing rule is simple: completed core, metadata alpha, hypothesis layer, and blocked scores stay separate. The next work is payload freeze, QA, and release-plus-paper alignment." |

## 20-Minute Section Rhythm

- Opening and positioning: 2:15
- Methods scaffold: 3:00
- Result spine: 6:50
- v7/v8 boundary: 2:15
- v9 platform/readiness: 2:25
- Extensions: 1:55
- Close: 1:10

The highest-value speaking time is slides 4-7 and 18-20. If a live audience
looks confused, slow down there and compress slides 11-14.

## 15-Minute Cut

Target runtime: 14:15, leaving 45 seconds of buffer.

| Slide | Time | Cut Instruction |
| --- | ---: | --- |
| 1 | 0:25 | State only the core question. |
| 2 | 0:30 | Say "distinct niche, not firstness"; skip named landscape details. |
| 3 | 0:35 | Name the three evidence layers and move on. |
| 4 | 0:55 | Explain the task contract; do not read all fields. |
| 5 | 0:45 | Define mission-held-out in one sentence. |
| 6 | 0:45 | Define feature views and AUROC grammar; skip transformed-input detail. |
| 7 | 0:40 | Use as the only full result-reading example. |
| 8 | 0:30 | One sentence: selected pathway rescue. |
| 9 | 0:30 | One sentence: scale alone does not transfer. |
| 10 | 0:45 | Keep the 8 x 8 x 4 = 256 hardening claim. |
| 11 | 0:35 | Say guardrails only; skip RRRM composition detail. |
| 12 | 0:30 | Say task limits, not absence of biology. |
| 13 | 0:35 | Say hypotheses, not mechanisms. |
| 14 | 0:40 | Say partial translation, not clinical transfer. |
| 15 | 0:25 | v7.1 is documentation boundary. |
| 16 | 0:35 | v8 is hypothesis incubator. |
| 17 | 0:25 | No countermeasure claim. |
| 18 | 0:50 | Keep metadata versus payload readiness. |
| 19 | 0:25 | Platform is reproducibility infrastructure. |
| 20 | 0:35 | 22/22 parsed, 0/22 hash-verified. |
| 21 | 0:25 | Organoid is draft extension. |
| 22 | 0:25 | OSD-120 is same-study diagnostic. |
| 23 | 0:35 | Single-cell is blocked by missing processed payload. |
| 24 | 0:55 | Close with separated claims and next work. |

## Fast-Forward Rules

Use these when the chair gives a warning or the audience needs more time on
methods:

- Merge slides 8-9 verbally: "Pathways help in selected settings; scale alone
  does not solve transfer."
- Merge slides 11-12 verbally: "Guardrail and negative results define the
  limits of the current benchmark."
- Merge slides 16-17 verbally: "v8 is hypothesis generation, not operational
  recommendation."
- Merge slides 21-22 verbally: "These are diagnostic extension tracks, not core
  benchmark claims."
- Never fast-forward slides 4, 5, 6, 7, 18, 20, or 23; those are where the
  audience learns how to read the deck.

## Lines To Avoid

- Do not say "first AI benchmark for space biology."
- Do not say "validated countermeasure" or "treatment implication."
- Do not say "Mars risk prediction" or imply operational astronaut-health
  recommendations.
- Do not say "v9 public bulk is released" or "payloads are frozen."
- Do not say "single-cell leaderboard" before the processed payload audit
  passes.
- Do not describe negative results as absence of biology.
- Do not describe pathway results as universally superior to gene-level models.

## Rehearsal Notes

- Slide 4 should feel like the first true teaching slide, not a methods
  footnote.
- Slide 5 should be delivered slowly enough that "hidden mission" lands.
- Slide 6 should make AUROC readable before slide 7 appears.
- Slide 7 is the first proof slide; treat it as a guided example.
- Slides 10-14 should feel like a hardened result spine, not five unrelated
  findings.
- Slide 18 is the key platform transition; it should reset the audience before
  v9 details.
- Slide 23 should be crisp and slightly slower: raw data exists, processed
  payload does not, therefore no score.

## Next Anchor

Continue from:

`SpaceBio-Bench v0.5 rehearsal timing -> speaker notes enrichment / institutional template decision`

Recommended next branch: either enrich the PPTX speaker notes with the 20-minute
talk track, or apply an institutional template after confirming whether the
target slot is 15 or 20 minutes.
