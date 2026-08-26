# MLCB 2026 Submission Draft

Status: working LaTeX draft
Target: MLCB 2026 8-page paper
Formatting basis checked 2026-05-10:

- 8-page paper, references excluded from page limit
- 11 pt font
- 1 inch margins
- single-column preferred
- single-blind review; anonymization is not required by the current MLCB
  submit page

## Build

From this directory:

```bash
tectonic main.tex
```

or, from the repository root:

```bash
tectonic submissions/mlcb_2026/main.tex
```

## Claim Source

All numerical claims in this draft should be checked against:

- `docs/CANONICAL_RESULTS_V7_1.md`
- `evaluation/RESULTS_SUMMARY.md`
- `docs/CONFERENCE_SUBMISSION_PACKET_2026_05_10.md`

Do not mix v8 SpaceMed intervention, countermeasure, or Mars-extrapolation
claims into this v7.1 benchmark paper.

## Next Edits

1. Replace placeholder figure boxes with static exports from existing benchmark
   artifacts.
2. Replace TODO bibliography entries with exact citation metadata.
3. Compress to 8 pages after figures are inserted.
4. Decide whether to opt into PMLR publication before final submission.
