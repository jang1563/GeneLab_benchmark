# Visual Identity Research And Depth Strategy

Date: 2026-06-01

Purpose: pause slide rendering and define the visual identity before rebuilding
more bridge slides. The previous bridge style emerged from production
constraints rather than a deliberate brand/design decision. This document
converts external research, internal critique, expert-role simulation, and
technical implementation constraints into a usable design direction.

## Executive Decision

Do not continue with the current "light methods/provenance grammar" as the
brand.

It is a useful production grammar, but it is not a visual identity. The target
identity should be:

> Mission-grade scientific benchmark system

This means:

- NASA/mission-systems seriousness without copying NASA marks or becoming
  decorative space wallpaper;
- high-profile scientific editorial restraint for figures and manuscript
  compatibility;
- premium strategy-deck clarity for bridge slides and first-time viewers;
- explicit z-depth so viewers can distinguish source evidence, measurement,
  interpretation, and caveat layers.

## What Went Wrong

The style became what it is because negative constraints arrived before
positive identity decisions.

We decided:

- no dashboard cards;
- no cheap card boxes;
- no decorative space background;
- no dense tables inside figures;
- use thin rules;
- keep caveats/source-traceability visible;
- separate method slides from dark result slides.

Those are good constraints, but they do not by themselves create a premium
brand. They produced careful, thin, off-white schematics. The result is
intellectually safe but visually underpowered.

## External Research Anchors

### NASA / Mission Systems

Research signals:

- NASA's public design language emphasizes strong identity control, clear
  typography, contrast, accessibility, and disciplined use of official marks.
- Modern NASA web/design-system references use systematized typography, color,
  components, and accessibility checks rather than ad hoc visual treatment.

Implication for SpaceBio-Bench:

- borrow mission-grade hierarchy, data-room restraint, and technical trust;
- do not imitate NASA logos, patches, or official insignia;
- do not use decorative starfield/space wallpaper;
- use orbital/mission cues only when they explain evaluation structure.

Useful sources:

- NASA Brand Center: https://www.nasa.gov/nasa-brand-center/
- NASA Graphics Standards Manual archive: https://www.nasa.gov/image-article/nasa-graphics-standards-manual/
- NASA Horizon Design System: https://hds.nasa.gov/

### High-Profile Scientific Editorial

Research signals:

- Journal figure guidance prioritizes legibility, clear labels, accessibility,
  and accurate separation of figure panels, legends, and tables.
- Graphical-abstract guidance favors a single clear story, minimal clutter, and
  self-contained interpretation.

Implication for SpaceBio-Bench:

- use figure-like discipline for result slides;
- keep exact rows and dense metadata in tables or appendix artifacts;
- make each main slide answer one question without relying on speaker repair;
- maintain manuscript variants that can remove atmospheric depth.

Useful sources:

- Nature final figure guidance: https://www.nature.com/nature/for-authors/final-submission
- Nature figure preparation: https://research-figure-guide.nature.com/figures/preparing-figures-our-specifications/
- Cell Press graphical abstracts: https://www.cell.com/pb/assets/raw/shared/figureguidelines/GA_guide.pdf

### Depth, Elevation, And Layer Semantics

Research signals:

- Interface systems treat elevation/depth as hierarchy and interaction meaning,
  not decoration.
- Shadows, overlays, and material surfaces should communicate what is above,
  below, active, or constrained.
- Depth becomes confusing when every object casts the same shadow or when
  nested panels mimic cards.

Implication for SpaceBio-Bench:

- z-depth must encode scientific meaning:
  evidence below interpretation, boundaries above evidence, caveats above all;
- one dominant depth event per slide is better than many small shadows;
- use crop windows, backplates, rails, and focus rings before heavy shadows.

Useful sources:

- Material Design elevation: https://m3.material.io/styles/elevation/overview
- Apple Human Interface Guidelines, materials: https://developer.apple.com/design/human-interface-guidelines/materials
- Adobe Spectrum tokens: https://spectrum.adobe.com/page/design-tokens/

## Professional Simulation

The following is an internal expert-role simulation, not an external human
review. It is designed to pressure-test the direction before more rendering.

### Role 1. Scientific Art Director

Assessment:

- Current bridge slides are too schematic and thin.
- They do not yet have the optical weight of an invited-talk or high-profile
  conference deck.
- The strongest path is not more labels; it is one large evidence object and a
  disciplined interpretation layer.

Recommendation:

- treat every bridge slide as an assertion-evidence slide;
- enlarge the scientific object until it carries the slide;
- keep the headline short enough to be remembered after the slide disappears.

### Role 2. NASA Mission Interface Designer

Assessment:

- The current style uses pale fields but lacks mission-system specificity.
- Mission-held-out evaluation should feel like a control room split, not a
  generic flowchart.

Recommendation:

- use mission rails, boundary gates, fold orbits, source ticks, and assay lanes
  as measurement infrastructure;
- reserve red for trust boundaries or hidden-test status only;
- avoid NASA-logo imitation and decorative aerospace motifs.

### Role 3. Journal Figure Editor

Assessment:

- The figure/table separation rule is correct.
- Some bridge text is too deck-native for manuscript reuse.

Recommendation:

- keep two export modes:
  deck mode with atmosphere/depth;
  manuscript mode with flatter background and tighter captions;
- remove informal phrases like `loose matrix` from manuscript variants;
- never put source inventory rows inside a main figure.

### Role 4. Scientific Strategy Deck Reviewer

Assessment:

- The current slides explain the method, but they do not yet sell the system.
- Premium strategy decks make the viewer feel the operating model immediately.

Recommendation:

- start with the operating model: public studies become auditable tasks;
- use fewer, larger objects;
- make before/after structure visible before reading any text.

### Role 5. Technical Production Lead

Assessment:

- The renderers are reproducible but tokenization is insufficient.
- Colors, elevation, layer order, and object semantics are still hard-coded in
  slide scripts.

Recommendation:

- introduce visual identity tokens;
- add z-depth tokens with semantic names;
- make every renderer declare which z-layers it uses;
- audit for card-box anti-patterns, excessive shadows, and missing evidence
  layers.

## Chosen Brand Territory

### Name

Mission-grade scientific benchmark system.

### Personality

- precise;
- source-traceable;
- mission-aware;
- calm but not weak;
- technical without becoming software UI;
- premium through restraint, spacing, and evidence, not ornament.

### Rejected Territories

- Generic academic beige schematic.
- NASA cosplay with patches, stars, orbit wallpaper, or agency-like branding.
- Consulting dashboard cards.
- Dark sci-fi control room.
- Startup SaaS product UI.
- Dense manuscript figure pasted into a slide without hierarchy.

## Visual System

### Typography

Deck mode:

- display/headline: geometric grotesk or Inter-like sans;
- body/callout: clean sans;
- data/code/source tick: mono only when it signals provenance.

Implementation fallback:

- Matplotlib should use available sans defaults unless a bundled font is
  explicitly available.
- Do not let font choice block reproducibility.

### Color

Use color by semantic role, not decoration.

Core roles:

- `ink`: headline and primary statements;
- `mission_navy`: deep system anchor;
- `source_blue`: public/source/provenance;
- `bio_green`: training/biological evidence;
- `assay_teal`: assay/tissue/modality context;
- `label_amber`: label contrast;
- `test_red`: held-out/test/trust boundary only;
- `model_purple`: model/fit operation;
- `paper`: evidence surface.

Rules:

- red is never a general highlight;
- purple should stay secondary and should not dominate the deck;
- avoid one-note beige by anchoring warm paper with navy/blue/green/teal;
- use color plus position/shape so grayscale remains interpretable.

### Layout

Slide types:

| Type | Primary Object | Depth Use | Notes |
|---|---|---|---|
| Thesis / B1 | fragmented sources compressing into evaluation surface | large atmospheric evidence field | must not look like agenda |
| Method bridge | one operating-model scene | measurement rail + trust boundary | large enough for thumbnail reading |
| Result slide | source-derived result proof object | proof crop + interpretation overlay | no decorative workflow |
| Provenance/release | document/evidence surface | document stack + focus windows | clear caveat/status layer |
| Manuscript figure | data figure | minimal or no atmosphere | captions carry caveats |

## Z-Stack And Depth Strategy

Depth must answer: what is evidence, what is measurement, what is
interpretation, and what is a caveat?

### Layer Contract

| Layer | Name | Meaning | Implementation |
|---|---|---|---|
| Z0 | Canvas atmosphere | project environment, not evidence | subtle warm/cool field, faint grain, no wallpaper |
| Z1 | Measurement infrastructure | rails, ticks, lanes, fold arcs | thin rules, low alpha, never dominant |
| Z2 | Evidence surface | source-derived or source-like object | large plate/band/crop, minimal border |
| Z3 | Proof objects | task record, mission nodes, result crop | sharper edges, slight elevation, semantic color |
| Z4 | Interpretation overlay | headline, callouts, focus labels | editable text, no paragraph blocks |
| Z5 | Trust boundary | held-out mission, train-only guard | restrained red line, focus window, gate |
| Z6 | Status/caveat | claim boundary, source-traceability | compact text, low visual weight |
| Z7 | Motion/focus | one directional path, crop window, transfer line | one dominant movement only |

### Depth Rules

1. Do not use depth as decoration.
2. One slide should have one dominant elevation event.
3. Avoid large bordered panels that read as card boxes.
4. Prefer full-width evidence bands, crop windows, and material shifts over
   floating cards.
5. Shadows should be directional and restrained:
   small offset, low opacity, tied to proof objects.
6. Text should not cast shadows.
7. If an object is not source evidence or a focus object, it should probably be
   flat.
8. Manuscript export should flatten Z0-Z2 and keep Z4-Z6 as captions/labels.

## Technical Implementation Update

### New Token File

Use:

- `config/visual_identity/spacebio_bench_visual_identity_v0_1.json`

This file defines:

- brand territory;
- color roles;
- typography roles;
- z-depth tokens;
- slide-type contracts;
- anti-patterns;
- export modes.

### Renderer Requirements

Future visual renderers should:

1. load visual identity tokens instead of hard-coding all colors;
2. declare slide type: `thesis`, `method_bridge`, `result`, `provenance`, or
   `manuscript`;
3. declare used z-layers in the manifest;
4. output both `scene_plate.png` and `overlay_spec.json`;
5. include a `depth_audit` field in `qa.json`;
6. produce a contact sheet for family review;
7. fail QA if the slide becomes a card-box layout.

### Suggested Renderer API

```python
tokens = load_visual_identity("config/visual_identity/spacebio_bench_visual_identity_v0_1.json")
canvas = draw_canvas(ax, tokens, slide_type="method_bridge")
draw_measurement_layer(ax, tokens, rail="mission")
draw_evidence_surface(ax, tokens, surface="source_trace")
draw_proof_object(ax, tokens, object_type="task_contract")
draw_interpretation_overlay(ax, overlay_spec, tokens)
write_layer_manifest(slide_id, used_layers=["Z0", "Z1", "Z2", "Z3", "Z4", "Z5", "Z6"])
```

### Next Implementation Step

Before rebuilding B1 or revising B2-B4 again:

1. create one style-frame contact sheet with three directions:
   `editorial`, `mission_system`, `hybrid`;
2. score each direction against:
   scientific trust,
   premium impression,
   thumbnail recognizability,
   manuscript portability,
   technical reproducibility;
3. choose one direction and lock visual tokens;
4. only then render B1-B4 as a family.

## Immediate Decision For Current Bridge Assets

Current status:

- `output/premium_bridge_scenes/`: logic wireframes only;
- `output/premium_bridge_rebuild_scenes/`: better candidates, but still not a
  locked brand system.

Decision:

- do not create more final-looking bridge slides until a style-frame choice is
  made;
- use the new token file as the production source of truth;
- treat B1 as the first real test of the chosen identity.
