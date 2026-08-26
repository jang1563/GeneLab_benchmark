# SpaceBio-Bench Conference Deck v0.5 - Platform Slide 02 Visual Review

Date: 2026-06-09

## Scope

- Deck: `output/spacebiobench_conference_deck_v0_5/spacebiobench_conference_deck_v0_5.pptx`
- Slide reviewed: slide 20, `Manifests, evaluator, and run records`
- Review anchor: make the v9 platform layer read as reproducibility infrastructure, not a new biological claim or leaderboard.

## Change Made

Replaced the original document-scene screenshot with an editable native object-chain diagram:

- `1 MANIFEST` - `Declare inputs`
- `2 EVALUATOR` - `Run metric code`
- `3 RUN RECORD` - `Store the trace`
- bottom rule: `Reproducibility surface`

Speaker notes were updated so the talk track starts from the left-to-right object chain.

## Review

The original screenshot asset had sufficient resolution but weak slide-level readability. The document scene placed the useful evidence in small, partially visible panels over a bright background, so first-time viewers could not quickly infer how manifests, evaluator outputs, and run records connect.

The replacement diagram makes the platform contract explicit:

- manifests declare the task inputs,
- the evaluator runs metric code,
- run records store the trace,
- together they create an audit/recompute surface.

This is more useful for a mixed audience than preserving a low-readability screenshot.

## Zoom / Audience QA

Additional zoom crops were inspected from the rendered slide:

- `outputs/019e8bd1-3b42-7821-9e1c-beebfa8e2ece/presentations/spacebiobench-conference-deck/qa/slide20_zoom/01_header_claim.png`
- `outputs/019e8bd1-3b42-7821-9e1c-beebfa8e2ece/presentations/spacebiobench-conference-deck/qa/slide20_zoom/02_object_chain.png`
- `outputs/019e8bd1-3b42-7821-9e1c-beebfa8e2ece/presentations/spacebiobench-conference-deck/qa/slide20_zoom/03_repro_surface.png`
- `outputs/019e8bd1-3b42-7821-9e1c-beebfa8e2ece/presentations/spacebiobench-conference-deck/qa/slide20_zoom/04_bottom_strip.png`

Audience-readability judgment:

- Title, subtitle, top-right claim, three object cards, arrows, reproducibility-surface rule, bottom strip, and footer guardrail are readable.
- No overlaps or cramped regions were observed.
- The object-chain diagram is more legible than the original screenshot scene.
- The slide now carries the claim through editable text and shapes rather than small raster annotations.

Resolution judgment:

- Rendered QA preview: 1280 x 720.
- Slide 20 has no embedded raster image asset in the PPTX after the edit.
- Projection risk is low because the slide is vector/text based and uses large type.

## QA

- PPTX build completed.
- Rendered slide 20 inspected visually.
- Zoom crops inspected.
- Contact sheet inspected for slide 19-21 platform flow.
- PPTX XML: 25 slides.
- Speaker notes XML: 25 notes.
- Slide 20 object-chain text present: yes.
- Slide 20 reproducibility-surface text present: yes.
- Slide 20 no-new-score/leaderboard rule present: yes.
- Slide 20 note mentions object chain: yes.
- Forbidden visible phrase `single-cell leaderboard`: none.
- Extension terms remain confined to later extension slides:
  - slide 22: `organoid`
  - slide 23: `OSD-120`
- Extension-term hits in slides 1-21: none.
- PDF export completed.
- PDF page count: 25.
- Output sizes:
  - PPTX: 31,239,067 bytes
  - PDF: 6,580,935 bytes
  - contact sheet: 3,831,074 bytes
  - speaker notes: 14,816 bytes

## Decision

Keep the slide 20 edit. The native object-chain diagram is clearer, more accessible, and better aligned with the reproducibility-infrastructure claim.
