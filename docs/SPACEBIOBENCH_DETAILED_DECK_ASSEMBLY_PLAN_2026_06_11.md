# SpaceBio-Bench Detailed Deck Assembly Plan

Date: 2026-06-11

Purpose: convert the detailed-deck strategy into an assembly-ready slide order.
This plan treats the detailed deck as a guided technical walkthrough for a mixed
computational biology, space biology, and ML/LLM audience.

## Assembly Decision

Build the deck in three passes:

1. Core walkthrough, slides 1-33: task, leakage control, metric/model primers,
   core transfer result, pathway conservation, held-out validation, model
   comparison, DGE robustness, and external biology support.
2. Biology and translation layer, slides 34-49: temporal, single-cell, spatial,
   systems biology, mouse-human, and v8 translational incubation.
3. Platform and release layer, slides 50-60: reproducibility objects, public
   bulk status, payload readiness, organoid/OSD-120/single-cell extensions, and
   final claim-status close.

The first build should assemble slides 1-33 before expanding to the full
60-slide technical deck. This keeps the highest-value reasoning chain tight:

`task contract -> hidden mission -> leakage guard -> metric/model primer ->
tissue transfer -> biology explains transfer -> held-out validation ->
robustness -> external biology`.

## Ready-To-Place Proof Assets

These 33 assets are QA-ready at `3840x2160` and should anchor the first assembly
pass:

| Planned slide | Asset | Role |
|---:|---|---|
| 1 | SpaceBio-Bench opening title | Introduces the mission-held-out benchmark thesis and reading contract. |
| 2 | Why this needs a benchmark | Explains how public study records become transfer-evaluation evidence. |
| 3 | Evidence layers set claim strength | Shows v1-v9 as evidence layers and separates benchmark, translation, and platform-readiness claim strength. |
| 4 | What counts as a task | Defines the five fixed task fields using an A4 thymus task-record example. |
| 5 | From source record to matrix | Shows how source records, metadata, normalized expression, task filters, and matrix orientation create the benchmark input surface. |
| 6 | What the model actually sees | Defines sample rows, feature views, classifier input, and score-row fields for mixed ML / biology audiences. |
| 7 | Mission-held-out protocol | Explains the evaluation unit. |
| 8 | Train-only leakage guard | Makes the leakage-control rule visible. |
| 9 | Metric primer | Teaches AUROC, uncertainty, chance, and GO gating. |
| 10 | Model-family bridge | Shows one fixed task contract across model families. |
| 11 | Model-family primer | Explains the three input surfaces. |
| 12 | Evidence scope ladder | Teaches how to separate score, robustness, biology, and follow-up evidence. |
| 13 | Worked tissue score | Teaches how to read one tissue row. |
| 14 | Tissue transfer hierarchy | Shows the ranked core tissue result. |
| 15 | Liver heterogeneity explainer | Explains why more missions can reveal lower transfer consistency. |
| 16 | Transfer matrix behind ranking | Shows the directed mission-pair evidence behind the tissue means. |
| 17 | Pathway nuisance-signal check | Shows how pathway features reduce selected coupled-label signals. |
| 18 | NES conservation predicts transfer | Connects pathway conservation to transfer performance. |
| 19 | Screen transfer feasibility before training | Turns NES conservation into a preflight model-triage workflow. |
| 20 | Held-out validation | Shows RR-23 thymus and RR-7 skin independent mission checks. |
| 21 | Negative controls | Shows shuffled-label, random-feature, permutation, and housekeeping audit gates. |
| 22 | Strong baselines matter | Shows the matched classical floor that model claims must clear. |
| 23 | Classical ML result surface | Shows tissue, model, and feature-view structure across the v4 grid. |
| 24 | What a foundation model adds | Shows pretraining, input translation, encoder, adapter, and shared score gate. |
| 25 | Text LLM checks are prompt diagnostics | Shows the prompt input surface, parser, parse rates, and prompt-diagnostic readout. |
| 26 | Task fit beats scale | Shows shared six-task model-tier readout and tissue-local deltas versus matched classical floors. |
| 27 | Bulk RNA-seq is a hard surface for cell FMs | Shows single-cell pretraining, adaptation pinch points, and bulk mission-held-out scoring surface. |
| 28 | Method hardening preserves the main readout | Shows the v4 8 x 8 x 4 grid, tissue-by-feature readout, and top-row stability. |
| 29 | Newer model ideas preserve the lesson | Shows scPRINT-2 and GNN probes under the same paired benchmark readout. |
| 30 | DGE robustness | Shows stable log2FC ranks across DGE callers. |
| 31 | External biology validation | Shows Cell 2020 pathway concordance and SHAP gene enrichment support. |
| 32 | Evidence stack | Shows how split, leakage, metric, model breadth, robustness, and biology layers set claim status. |
| 33 | Core benchmark takeaway | Closes Acts 1-3 and bridges the reader into biology interpretation. |

## Flow And Rhythm

Use a recurring rhythm:

1. Teach the reading rule.
2. Show the proof object.
3. Add scope/readout.
4. Move to the next layer.

This matters because the audience includes people who may not be fluent in ML,
LOMO evaluation, foundation models, or pathway enrichment. The deck should make
the method visible on the slide, with spoken explanation reinforcing what the
audience can already see.

## Key Sequencing Choices

- Put mission-held-out and train-only guard before any score.
- Put AUROC and model-family primers before model comparison.
- Put the worked tissue score before the full tissue ranking.
- Put liver heterogeneity immediately after the tissue hierarchy so sample count
  stays separate from transfer strength.
- Put the pathway feature check before NES conservation: first show what the
  feature view changes, then show why conserved pathways predict transfer.
- Put NES conservation before the preflight feasibility slide: first establish
  the conservation-performance link, then translate it into a model-triage
  workflow.
- Put the preflight feasibility slide before held-out validation: the audience
  sees how fGSEA can prioritize modeling, then held-out missions test the
  generalization pattern.
- Put DGE robustness before external biology validation: first verify that the
  pathway/DGE layer is stable across caller choices, then show biological
  concordance with published literature.
- Put the strong-baseline floor before the model result section: first define
  the matched benchmark floor, then compare classical ML, foundation models, and
  prompt-only diagnostics.
- Put the classical ML result surface immediately after the baseline floor:
  first show that the floor is structured by tissue, model family, and feature
  view, then ask what pretrained models add.
- Put the foundation-model primer before text LLM and model-tier result slides:
  first separate pretrained expression encoders from prompt-only language checks,
  then compare tested-setting outcomes.
- Put the model-family result before DGE/external validation so the audience
  understands what kind of benchmark signal is being supported.

## Slide Spine V2

| # | Act | Slide Title | Main Question | Proof Object / Asset Status |
|---:|---|---|---|---|
| 1 | Open | SpaceBio-Bench | What is the thesis? | Native title scene. |
| 2 | Open | Why This Needs A Benchmark | How do public study records become transfer-evaluation evidence? | Ready asset: OSDR/study/mission benchmark-need map. |
| 3 | Open | Evidence Layers Set Claim Strength | How do v1-v9 layers relate? | Ready asset: evidence-layer claim-strength map. |
| 4 | Method | What Counts As A Task | What must be fixed before scoring? | Ready asset: five-field task-record anatomy. |
| 5 | Method | From Source Record To Matrix | How does public data become benchmark input? | Ready asset: source-to-matrix construction map. |
| 6 | Method | What The Model Actually Sees | What is the sample-by-feature input surface? | Ready asset: gene/pathway/compressed input bridge. |
| 7 | Method | The Test Set Is An Entire Mission | What is the hidden unit? | Ready asset: mission-held-out protocol. |
| 8 | Method | Training Choices Stop Before The Hidden Mission | How is leakage controlled? | Ready asset: train-only leakage guard. |
| 9 | Method | How To Read AUROC | What does a score mean? | Ready asset: metric primer. |
| 10 | Method | One Task Contract, Three Model Families | How do model families stay comparable? | Ready asset: model-family bridge. |
| 11 | Method | Same Benchmark, Three Input Surfaces | What does each model family receive? | Ready asset: model-family primer. |
| 12 | Method | What Counts As Evidence | What is benchmark evidence, biological support, and scope? | Ready asset: evidence scope ladder. |
| 13 | Core Result | Read One Tissue Score | How should one row be interpreted? | Ready asset: worked tissue score. |
| 14 | Core Result | Thymus Leads The Transfer Hierarchy | Which tissues generalize best? | Ready asset: tissue transfer hierarchy. |
| 15 | Core Result | More Liver Missions Expose Heterogeneity | Why does liver rank lower? | Ready asset: liver heterogeneity explainer. |
| 16 | Core Result | The Transfer Matrix Behind The Ranking | Which mission-pair evidence creates the mean? | Ready asset: transfer matrix behind ranking. |
| 17 | Core Result | Pathway Features Reduce Selected Nuisance Signals | How do pathway summaries change the task? | Ready asset: pathway nuisance-signal check. |
| 18 | Core Result | Conserved Pathways Predict Transfer | Why might some tissues generalize? | Ready asset: NES conservation predicts transfer. |
| 19 | Core Result | Screen Transfer Feasibility Before Training | How can fGSEA guide model triage? | Ready asset: transfer-feasibility preflight workflow. |
| 20 | Core Result | Held-Out Missions Confirm The Signal | Does the pattern survive reserved missions? | Ready asset: held-out validation. |
| 21 | Core Result | Negative Controls Anchor The Readout | Where should signal drop toward chance? | Ready asset: negative-control gate summary. |
| 22 | Models | Strong Baselines Matter | What is the benchmark floor? | Ready asset: strong-baseline floor. |
| 23 | Models | Classical ML Result Surface | How strong are tuned classical models? | Ready asset: classical ML result surface. |
| 24 | Models | What A Foundation Model Adds | What does pretraining contribute? | Ready asset: foundation-model pretrain-to-adapt diagram. |
| 25 | Models | Text LLM Checks Are Prompt Diagnostics | What is the LLM input surface? | Ready asset: text LLM prompt diagnostic diagram. |
| 26 | Models | Task Fit Beats Scale | What happens across model tiers? | Ready asset: model-tier shared readout and tissue-local deltas. |
| 27 | Models | Bulk RNA-seq Is A Hard Surface For Cell FMs | Why is transfer hard for single-cell-pretrained models? | Ready asset: bulk adaptation-surface explainer. |
| 28 | Models | Method Hardening Preserves The Main Readout | Does v4 broaden the evidence surface? | Ready asset: v4 method-hardening grid. |
| 29 | Models | Newer Model Ideas Preserve The Lesson | What do scPRINT/GNN checks add? | Ready asset: v7 scPRINT/GNN paired comparison. |
| 30 | Robustness | DGE Ranks Hold; DEG Lists Move | Is the DGE layer stable across callers? | Ready asset: DGE robustness. |
| 31 | Robustness | Published Biology Supports The Pathway Signal | Does literature biology align? | Ready asset: external biology validation. |
| 32 | Robustness | Evidence Stack Turns Scores Into Claim Status | How do split, leakage, metric, DGE, and biology checks work together? | Ready asset: integrated evidence stack. |
| 33 | Robustness | Core Benchmark Takeaway | What can the first three acts claim? | Ready asset: recap bridge into biology layer. |
| 34 | Biology | Temporal Labels Need Guardrails | How do preservation/recovery labels behave? | Ready asset: temporal context guardrail. |
| 35 | Biology | Single-Cell Pilots Provide Context | What do RRRM pilots contribute? | Ready asset: RRRM single-cell pilot context. |
| 36 | Biology | Spatial And Weak-Signal Cases Define Scope | How do weak-signal cases inform the benchmark? | Ready asset: spatial weak-signal scope boundary. |
| 37 | Biology | Systems Biology Adds Interpretation | How do immune, TF, metabolic, target, and biomarker layers connect? | Ready asset: five-lane systems-biology interpretation map. |
| 38 | Biology | Immune And TF Activity Prioritize Follow-Up | Which biology axes appear most interpretable? | Ready asset: immune/TF follow-up priority map. |
| 39 | Biology | Target And Biomarker Layers Are Triage | How should downstream target evidence be used? | Ready asset: target funnel and biomarker panel triage map. |
| 40 | Translation | Mouse-To-Human Transfer Is Partial | Which signals cross species? | Source-remake from v6 translation. |
| 41 | Translation | Translation Details Matter | Which ortholog/cfRNA/TF/target constraints matter? | Native compact evidence table. |
| 42 | Translation | Prediction To Interpretation To Hypothesis | How should biology claims be layered? | Native triage summary. |
| 43 | v8 | Why Translation Appears After Benchmarking | Why place v8 after the evidence stack? | Native transition. |
| 44 | v8 | Four Incubator Pillars | What are BRIDGE, DECOMPOSE, INTERVENE, CAUSAL? | Native v8 pillar map. |
| 45 | v8 | Mouse Pathways Improve Human Prediction | What does BRIDGE add? | Reuse/remake v8 Figure1. |
| 46 | v8 | Stressor Effects Interact | What does DECOMPOSE show? | Reuse/remake v8 Figure2. |
| 47 | v8 | Perturbation Hits Prioritize Hypotheses | What does INTERVENE add? | Native prioritization map from v8 Figure4. |
| 48 | v8 | Causal Maps Organize Assumptions | What does CAUSAL add? | Simplified DAG from v8 Figure5. |
| 49 | v8 | Incubator Summary | What claim layer does v8 occupy? | Native scope summary. |
| 50 | Platform | Platform Turn | How does the work become reproducible infrastructure? | Native transition. |
| 51 | Platform | Manifest, Evaluator, Run Record | What objects make runs reproducible? | Reuse/remake v9 architecture. |
| 52 | Platform | Task Registry And Metric Profiles | How do typed tasks and metrics keep runs comparable? | Native registry/metric diagram. |
| 53 | Platform | Public Bulk Metadata Alpha | What exists in the public bulk alpha? | Reuse v9 public bulk scope scene. |
| 54 | Platform | Payload Readiness Ladder | What separates parsed source from frozen release? | Native readiness ladder. |
| 55 | Platform | Data Package Scope | What does the Data Package prove? | Reuse v9 document scene. |
| 56 | Extension | Organoid Extension Is A Biology Check | Where does organoid work fit? | Reuse/remake organoid diagnostic surface. |
| 57 | Extension | OSD-120 Is Same-Study Diagnostic | How should OSD-120 be read? | Reuse/remake source/split proof. |
| 58 | Extension | Single-Cell Scoring Needs Payload Freeze | What must be frozen for public single-cell scoring? | Native payload-freeze diagram. |
| 59 | Close | Roadmap To Release | What remains before release/paper alignment? | Native roadmap. |
| 60 | Close | What We Can Claim | What is the final claim-status map? | Native close slide. |

## Immediate Build Order

1. Slides 7-15 have been assembled into a first Act 1-2 contact sheet and pass
   sequence QA.
2. Slides 7-33 ready proof assets have been assembled into a core proof-asset
   contact sheet and pass sequence, grayscale, and copy-tone QA.
3. Slides 1-33 have been assembled into a canonical walkthrough contact sheet
   and pass sequence, grayscale, thumbnail-distinctness, and flow QA.

## Open Build Risks

- Slide 34 is now ready as the first Act 4 biology-context proof asset.
- Slides 35-60 remain proof-asset build/remake work before a full 60-slide
  walkthrough can be assembled.
- Slide 12 is now ready and should be used as the reader contract before the
  result section.
- The older 60-slide spine uses legacy scope wording in places; visible slide
  copy should use `Scope`, `Readout`, and claim-status language.
- Acts 4-6 should be assembled after either the full first 33-slide walkthrough
  passes contact-sheet QA or the user explicitly chooses to continue beyond the
  proof-asset gate.
