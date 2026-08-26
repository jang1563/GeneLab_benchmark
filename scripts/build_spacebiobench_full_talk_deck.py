#!/usr/bin/env python3
"""Build the SpaceBio-Bench 24-slide full-talk editable PPTX draft.

The script writes temporary slide modules inside a Presentations skill workspace,
then uses the bundled artifact-tool builder to export PPTX, rendered previews,
layout JSON, and a contact sheet. Final deliverables are copied into the
project's output directory for the SpaceBio-Bench deck lineage.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CREATED = "2026-06-03"
THREAD_ID = os.environ.get("CODEX_THREAD_ID", "manual-20260603-spacebiobench-full-talk-v0-3")
TASK_SLUG = "spacebiobench-full-talk-deck"
WORKSPACE = ROOT / "outputs" / THREAD_ID / "presentations" / TASK_SLUG
SLIDES_DIR = WORKSPACE / "slides"
PREVIEW_DIR = WORKSPACE / "preview"
LAYOUT_DIR = WORKSPACE / "layout" / "final"
QA_DIR = WORKSPACE / "qa"
WORKSPACE_OUTPUT_DIR = WORKSPACE / "output"
PROJECT_OUTPUT_DIR = ROOT / "output" / "spacebiobench_full_talk_deck_v0_3"

SKILL_DIR = Path(
    "<PATH>"
)
NODE = Path(
    os.environ.get(
        "PRESENTATIONS_NODE",
        "<PATH>",
    )
)
BUILD_ARTIFACT_DECK = SKILL_DIR / "scripts" / "build_artifact_deck.mjs"

AUDIT_JSON = ROOT / "output" / "slides_1_14_pptx_readiness_audit_v0_1" / "slides_1_14_pptx_readiness_audit.json"
FINAL_PPTX = WORKSPACE_OUTPUT_DIR / "spacebiobench_full_talk_deck_v0_3.pptx"
CONTACT_SHEET = PREVIEW_DIR / "spacebiobench_full_talk_deck_contact_sheet.png"
BUILD_MANIFEST = WORKSPACE_OUTPUT_DIR / "spacebiobench_full_talk_deck_artifact_manifest.json"


COLORS = {
    "void": "#070A0E",
    "panel": "#0B1117",
    "panel2": "#101923",
    "ink": "#F4F7F8",
    "soft": "#B9C7D2",
    "muted": "#788898",
    "rule": "#33465A",
    "teal": "#1AA090",
    "sky": "#6BAED6",
    "amber": "#B4832F",
    "rose": "#B45A7E",
    "green": "#178B63",
    "violet": "#7B68A8",
}

SECTION_ACCENT = {
    "Opening": "teal",
    "Core methods": "sky",
    "v1-v7 results": "amber",
    "core result spine": "green",
    "v8 incubator": "rose",
    "v9 platform": "sky",
    "v9 extension": "violet",
    "Close": "teal",
}

EXTENSION_SLIDES = [
    {
        "slide": 15,
        "section": "v1-v7 results",
        "title": "v7.1 boundary",
        "role": "release boundary",
        "headline": "Canonical result surface, not a new result",
        "asset": "",
        "native": "v7_boundary",
        "overlay": {
            "eyebrow": "V1-V7 RESULT BOUNDARY",
            "title": "Canonical result surface, not a new result",
            "subtitle": "v7.1 freezes public-facing counts, subset labels, and claim discipline.",
            "bridge": "Separate completed benchmark evidence from translation and platform layers.",
            "caveat": "Documentation consistency patch only; no v8 intervention or Mars claim.",
            "source": "docs/CANONICAL_RESULTS_V7_1.md",
            "layout": "core",
        },
        "notes": (
            "Slide 15: v7.1 boundary\n"
            "Claim: Canonical result surface, not a new benchmark run.\n"
            "Claim boundary: v7.1 is a documentation consistency patch; v8 SpaceMed claims do not enter the v7.1 benchmark paper."
        ),
    },
    {
        "slide": 16,
        "section": "v8 incubator",
        "title": "v8 SpaceMed",
        "role": "hypothesis incubator",
        "headline": "Translation is hypothesis generation only",
        "asset": "v8/figures/Figure1_Species_Transfer.png",
        "assetFit": "contain",
        "overlay": {
            "eyebrow": "V8 INCUBATOR",
            "title": "Translation is hypothesis generation only",
            "subtitle": "Species-transfer views are downstream hypotheses, not benchmark proof.",
            "bridge": "After the completed benchmark spine, v8 becomes a bounded translational incubator.",
            "caveat": "No clinical, countermeasure, or operational recommendation claim.",
            "source": "v8 Figure 1 + docs/V8_BETA_RELEASE_PLAN_2026_05_10.md",
            "layout": "extension",
        },
        "notes": (
            "Slide 16: v8 SpaceMed\n"
            "Claim: translation is hypothesis generation only.\n"
            "Presenter bridge: this is downstream use of benchmark signals, not a new validated intervention surface."
        ),
    },
    {
        "slide": 17,
        "section": "v8 incubator",
        "title": "v8 claim boundary",
        "role": "hypothesis boundary",
        "headline": "Stressors are not countermeasure recommendations",
        "asset": "v8/figures/Figure2_Stressor_Decomposition.png",
        "assetFit": "contain",
        "overlay": {
            "eyebrow": "V8 CLAIM BOUNDARY",
            "title": "Stressors are not countermeasure recommendations",
            "subtitle": "The beta machinery hardens provenance before any release claim.",
            "bridge": "Use v8 as a controlled hypothesis layer, not an astronaut-health recommendation.",
            "caveat": "No Mars point estimate, treatment, or countermeasure recommendation claim.",
            "source": "v8 Figure 2 + v8 beta provenance/release plan",
            "layout": "extension",
        },
        "notes": (
            "Slide 17: v8 claim boundary\n"
            "Claim: no countermeasure or Mars-risk claim.\n"
            "Presenter bridge: release machinery validates provenance; it does not make the biology operational."
        ),
    },
    {
        "slide": 18,
        "section": "v9 platform",
        "title": "v9 platform",
        "role": "platform architecture",
        "headline": "Manifests, evaluator, and run records",
        "asset": "output/premium_v9_document_scenes/v9_platform_provenance_document_scene_plate.png",
        "assetFit": "cover",
        "overlay": {
            "eyebrow": "V9 PLATFORM",
            "title": "Manifests, evaluator, and run records",
            "subtitle": "v9 turns result artifacts into auditable benchmark infrastructure.",
            "bridge": "The platform layer is about reproducibility surfaces, not new biological claims.",
            "caveat": "Platform provenance surface; not a frozen release or leaderboard claim.",
            "source": "v9 platform provenance document scene + v9 manifests/evaluator outputs",
            "layout": "extension",
        },
        "notes": (
            "Slide 18: v9 platform\n"
            "Claim: v9 platformizes the benchmark with manifests, evaluator, run records, and provenance."
        ),
    },
    {
        "slide": 19,
        "section": "v9 platform",
        "title": "Public bulk alpha",
        "role": "release boundary",
        "headline": "Metadata alpha, payload hashes blocked",
        "asset": "output/premium_v9_document_scenes/v9_public_bulk_boundary_document_scene_plate.png",
        "assetFit": "cover",
        "overlay": {
            "eyebrow": "PUBLIC BULK ALPHA",
            "title": "Metadata alpha, payload hashes blocked",
            "subtitle": "22/22 sources have parsed checksum-manifest evidence; 0/22 payloads are locally hash-verified.",
            "bridge": "A metadata snapshot can be useful without pretending to be a frozen payload release.",
            "caveat": "Metadata-only alpha; no frozen payload mirror, DOI release, or leaderboard claim.",
            "source": "docs/V9_PUBLIC_BULK_ALPHA_METADATA_SNAPSHOT_DECISION.md",
            "layout": "extension",
        },
        "notes": (
            "Slide 19: public bulk alpha\n"
            "Claim: metadata alpha is allowed with explicit payload blockers.\n"
            "Do not call it a frozen public benchmark release."
        ),
    },
    {
        "slide": 20,
        "section": "v9 extension",
        "title": "Organoid extension",
        "role": "extension proof",
        "headline": "Source records become local matrices",
        "asset": "output/biovis_organoid_audience_matrix_proof_v0_2/panels/01_dark_organoid_clean_source_to_matrix.png",
        "assetFit": "cover",
        "overlay": {
            "eyebrow": "V9 EXTENSION",
            "title": "Source records become local matrices",
            "subtitle": "Human organoid evidence stays a draft extension track, not the public bulk core.",
            "bridge": "Use one organoid bridge here, after the completed result and platform sections.",
            "caveat": "Draft pilot extension only; not mission-held-out benchmark evidence.",
            "source": "organoid source-to-matrix proof panel + v9 organoid draft manifests",
            "layout": "extension",
        },
        "notes": (
            "Slide 20: organoid extension\n"
            "Claim: source records become local expression matrices.\n"
            "Boundary: draft pilot extension, not public bulk alpha core."
        ),
    },
    {
        "slide": 21,
        "section": "v9 extension",
        "title": "Multispecies task check",
        "role": "extension proof",
        "headline": "Same-study, not mission-held-out",
        "asset": "output/biovis_osd120_audience_split_proof_v0_2/panels/01_dark_osd120_clean_source_to_task.png",
        "assetFit": "cover",
        "overlay": {
            "eyebrow": "V9 EXTENSION",
            "title": "Same-study, not mission-held-out",
            "subtitle": "OSD-120 is useful as an interaction diagnostic, but it is not the flagship transfer task.",
            "bridge": "Keep the OSD-120 bridge honest: compactness is not cross-mission generalization.",
            "caveat": "Diagnostic same-study task only; no cross-species or leaderboard claim.",
            "source": "OSD-120 source-to-task proof panel + multispecies diagnostic reports",
            "layout": "extension",
        },
        "notes": (
            "Slide 21: multispecies task check\n"
            "Claim: OSD-120 is same-study, not mission-held-out.\n"
            "Boundary: diagnostic extension only."
        ),
    },
    {
        "slide": 22,
        "section": "v9 extension",
        "title": "Single-cell extension",
        "role": "blocker status",
        "headline": "Metric spec exists; payload blocker remains",
        "asset": "",
        "native": "single_cell_blocker",
        "overlay": {
            "eyebrow": "V9 EXTENSION",
            "title": "Metric spec exists; payload blocker remains",
            "subtitle": "OSD-918 blood has manifest and metric scaffolding, but no processed h5ad or STARsolo payload is available.",
            "bridge": "Do not promote RRRM single-cell scores until the payload passes audit.",
            "caveat": "No canonical h5ad copy, regeneration, evaluator, or legacy-score promotion claim.",
            "source": "docs/V9_SC_OSDR_PROCESSED_PAYLOAD_DISCOVERY.md",
            "layout": "extension",
        },
        "notes": (
            "Slide 22: single-cell extension\n"
            "Claim: inventory and metric spec exist; payload blocker remains.\n"
            "OSDR lists raw FASTQ pairs for 8/8 blood SRX, but no processed h5ad or STARsolo bundle."
        ),
    },
    {
        "slide": 23,
        "section": "Close",
        "title": "Claim boundary",
        "role": "claim matrix",
        "headline": "Keep claims separated",
        "asset": "",
        "native": "claim_matrix",
        "overlay": {
            "eyebrow": "CLAIM BOUNDARY",
            "title": "Keep claims separated",
            "subtitle": "The deck is strongest when completed evidence, alpha metadata, and blocked extensions stay separate.",
            "bridge": "Use the closing matrix to say what is claimed, what is exploratory, and what is explicitly not claimed.",
            "caveat": "No frozen v9 release, no v8 countermeasure, no single-cell leaderboard claim.",
            "source": "v7.1 canonical boundary + v8/v9 release-boundary docs",
            "layout": "extension",
        },
        "notes": (
            "Slide 23: claim boundary\n"
            "Use this as the closing guardrail before roadmap.\n"
            "Keep completed results separate from alpha metadata, hypothesis-only translation, and blocked tracks."
        ),
    },
    {
        "slide": 24,
        "section": "Close",
        "title": "Roadmap",
        "role": "roadmap",
        "headline": "Freeze payloads, QA, release, manuscript",
        "asset": "",
        "native": "roadmap",
        "overlay": {
            "eyebrow": "ROADMAP",
            "title": "Freeze payloads, QA, release, manuscript",
            "subtitle": "The next work turns polished evidence into a release-ready benchmark package.",
            "bridge": "Close with disciplined next steps, not inflated claims.",
            "caveat": "Roadmap only; no claim that all release blockers are resolved.",
            "source": "docs/RELEASE_ROADMAP.md + v8/v9 operating backlog",
            "layout": "extension",
        },
        "notes": (
            "Slide 24: roadmap\n"
            "Claim: next work is payload freeze, QA, release packaging, and manuscript assembly.\n"
            "Boundary: roadmap only."
        ),
    },
]


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def clean_text(text: str) -> str:
    return re.sub(r"\s+", " ", text or "").strip()


def wrap_source(text: str, max_len: int = 92) -> str:
    text = clean_text(text)
    if len(text) <= max_len:
        return text
    return text[: max_len - 1].rstrip() + "."


POLISHED_SOURCE_LINES = {
    1: "public GeneLab/OSDR study records + v1-v7 benchmark core",
    2: "external landscape scan + v1-v9 positioning notes",
    3: "V1-V9 master outline + deck assembly bridge review",
    4: "B1/B2 narration pack + source/task inventory",
    5: "B3/B4 narration pack + mission-held-out split rule",
    6: "feature-layer bridge overlay spec + v1-v7 benchmark core",
    7: "v1-v7 mouse bulk LOMO AUROC + tissue-transfer scene",
    8: "v1-v7 pathway-summary diagnostics + feature-layer QA",
    9: "shared 6-task model comparison + tissue-delta summaries",
    10: "v4 hardening summary: 8 tissues x 8 classifiers x 4 feature views",
    11: "v2 temporal/RRRM review + RRRM single-cell pilot summaries",
    12: "v3 negative-boundary review + spatial/pretrained controls",
    13: "v5 biological-interpretation summaries + target-triage evidence",
    14: "v6 human-translation ladder + target-tier alignment summaries",
}

POLISHED_CAVEATS = {
    10: "Hardening and coverage evidence only; not a new leaderboard.",
    11: "Guardrail evidence only; recovery is descriptive and RRRM is underpowered.",
    12: "Task-limit evidence only; no universal absence-of-biology claim.",
    13: "Hypothesis and target-triage evidence only; no mechanism or treatment claim.",
    14: "Partial pathway and target-tier alignment only; no clinical claim.",
}


def js(value) -> str:
    return json.dumps(value, ensure_ascii=True)


def choose_overlay(row: dict) -> dict:
    visible = [clean_text(item) for item in row.get("visible_text_sample", []) if clean_text(item)]
    section = row["section"]
    slide_number = row["slide"]

    if slide_number in {1, 2, 3, 4, 5}:
        eyebrow = visible[0] if visible else section.upper()
        title = visible[1] if len(visible) > 1 else row["headline"]
        subtitle = visible[2] if len(visible) > 2 else ""
        bridge = visible[3] if len(visible) > 3 else ""
        caveat = visible[4] if len(visible) > 4 else row.get("claim_boundary", "")
        source = visible[5] if len(visible) > 5 else row.get("editable_overlay_source", "")
        if slide_number == 1:
            title = "SpaceBio-Bench"
            subtitle = "Testing biological AI under spaceflight domain shift."
        if slide_number == 3:
            title = "From benchmark results to platform"
            subtitle = "v1-v9 separates completed evidence, hypothesis-only translation, and platformization."
        return {
            "eyebrow": eyebrow,
            "title": title,
            "subtitle": subtitle,
            "bridge": bridge,
            "caveat": caveat,
            "source": POLISHED_SOURCE_LINES.get(slide_number, source),
            "layout": "opening" if slide_number <= 3 else "method",
        }

    if slide_number == 6:
        return {
            "eyebrow": "FEATURE LAYERS",
            "title": row["headline"],
            "subtitle": visible[1] if len(visible) > 1 else "Same samples, different feature views.",
            "bridge": "Define the feature view before reading score surfaces.",
            "caveat": row.get("claim_boundary", ""),
            "source": POLISHED_SOURCE_LINES[6],
            "layout": "method",
        }

    if slide_number in {7, 8, 9}:
        return {
            "eyebrow": "EARLY RESULT FAMILY",
            "title": visible[0] if visible else row["headline"],
            "subtitle": visible[1] if len(visible) > 1 else row["headline"],
            "bridge": row.get("headline", ""),
            "caveat": row.get("claim_boundary", ""),
            "source": POLISHED_SOURCE_LINES[slide_number],
            "layout": "result",
        }

    return {
        "eyebrow": "CORE RESULT SPINE",
        "title": visible[0] if visible else row["headline"],
        "subtitle": visible[1] if len(visible) > 1 else row["headline"],
        "bridge": "",
        "caveat": POLISHED_CAVEATS.get(
            slide_number,
            visible[2] if len(visible) > 2 else row.get("claim_boundary", ""),
        ),
        "source": POLISHED_SOURCE_LINES[slide_number],
        "layout": "core",
    }


def speaker_notes(row: dict, overlay: dict) -> str:
    lines = [
        f"Slide {row['slide']}: {row['title']}",
        f"Claim: {overlay['title']}",
        f"PPTX instruction: {row.get('pptx_instruction', '')}",
        f"Claim boundary: {row.get('claim_boundary', '')}",
    ]
    if row.get("viewer_language_watch_hits"):
        lines.append("Viewer-language watch: " + ", ".join(row["viewer_language_watch_hits"]))
    if row.get("speaker_notes_required"):
        lines.append("Required presenter bridge:")
        lines.extend(f"- {item}" for item in row["speaker_notes_required"])
    return "\n".join(line for line in lines if clean_text(line))


def build_deck_spec() -> list[dict]:
    audit = json.loads(AUDIT_JSON.read_text(encoding="utf-8"))
    deck_spec: list[dict] = []
    for row in audit["slides"]:
        overlay = choose_overlay(row)
        accent = COLORS[SECTION_ACCENT.get(row["section"], "teal")]
        deck_spec.append(
            {
                "slide": row["slide"],
                "section": row["section"],
                "title": row["title"],
                "role": row["role"],
                "asset": str((ROOT / row["scene_plate"]).resolve()),
                "accent": accent,
                "overlay": overlay,
                "notes": speaker_notes(row, overlay),
                "watch": row.get("viewer_language_watch_hits", []),
            }
        )
    for row in EXTENSION_SLIDES:
        slide_spec = {
            "slide": row["slide"],
            "section": row["section"],
            "title": row["title"],
            "role": row["role"],
            "asset": str((ROOT / row["asset"]).resolve()) if row.get("asset") else "",
            "assetFit": row.get("assetFit", "cover"),
            "native": row.get("native", ""),
            "accent": COLORS[SECTION_ACCENT.get(row["section"], "teal")],
            "overlay": row["overlay"],
            "notes": row["notes"],
            "watch": [],
        }
        deck_spec.append(slide_spec)
    return deck_spec


def write_workspace_notes(deck_spec: list[dict]) -> None:
    WORKSPACE.mkdir(parents=True, exist_ok=True)
    QA_DIR.mkdir(parents=True, exist_ok=True)
    profile = [
        "task mode: targeted-edit",
        "primary deck-profile: engineering-platform",
        "secondary profile gates: appendix-heavy source discipline; scientific figure claim-boundary discipline",
        "required proof objects: workflow/method bridges, feature-layer bridge, metric result proof scenes, translation-boundary slide",
        "source/asset requirements: use audited scene plates from local output; rebuild interpretation text as editable PPTX objects",
        "brand authenticity constraints: no external logos or invented identity marks",
        "profile-specific QA gates: preserve technical precision, keep claim attached to proof object, render contact sheet, inspect footer readability",
        "known missing inputs: no final institutional template supplied; this is a polished production draft",
    ]
    (WORKSPACE / "profile-plan.txt").write_text("\n".join(profile) + "\n", encoding="utf-8")

    source_lines = [
        "Source notes",
        f"Date: {CREATED}",
        "",
        "Primary local source: output/slides_1_14_pptx_readiness_audit_v0_1/slides_1_14_pptx_readiness_audit.json",
        "All slide proof surfaces use audited local scene plates. No external logo or identity asset is embedded.",
        "",
    ]
    for slide in deck_spec:
        asset = slide.get("asset")
        source_lines.append(
            f"Slide {slide['slide']}: {rel(Path(asset)) if asset else 'native editable scene'}"
        )
    (WORKSPACE / "source-notes.txt").write_text("\n".join(source_lines) + "\n", encoding="utf-8")

    claim_lines = ["Claim spine", f"Date: {CREATED}", ""]
    for slide in deck_spec:
        claim_lines.append(f"{slide['slide']}. {slide['overlay']['title']} | {slide['overlay']['caveat']}")
    (WORKSPACE / "claim-spine.txt").write_text("\n".join(claim_lines) + "\n", encoding="utf-8")

    contact_lines = [
        "Contact-sheet plan",
        f"Date: {CREATED}",
        "",
        "Macro-layout rhythm:",
        "1-3 opening scenes with atmospheric proof surfaces",
        "4-6 methods bridge scenes with workflow proof surfaces",
        "7-9 early result proof scenes",
        "10-14 hardened core result proof scenes",
        "",
        "Hard rule: no rendered preview should become the final all-in-one slide; text shell remains editable.",
        "Polish rule: footer caveats and source lines must remain readable at full-slide view and not drift into scene labels.",
    ]
    (WORKSPACE / "contact-sheet-plan.txt").write_text("\n".join(contact_lines) + "\n", encoding="utf-8")


SHARED_MODULE = r'''
const COLORS = {
  void: "#070A0E",
  panel: "#0B1117",
  panel2: "#101923",
  ink: "#F4F7F8",
  soft: "#B9C7D2",
  muted: "#788898",
  rule: "#33465A",
  sky: "#6BAED6",
  green: "#178B63",
  amber: "#B4832F",
};

function text(slide, ctx, opts) {
  return ctx.addText(slide, {
    ...opts,
    typeface: opts.typeface || "Aptos",
    insets: opts.insets || { left: 0, right: 0, top: 0, bottom: 0 },
  });
}

function rule(slide, ctx, x, y, w, color) {
  ctx.addShape(slide, {
    x,
    y,
    w,
    h: 2,
    geometry: "rect",
    fill: color,
    line: ctx.line("#00000000", 0),
  });
}

function matte(slide, ctx, y, h) {
  ctx.addShape(slide, {
    x: 0,
    y,
    w: ctx.W,
    h,
    geometry: "rect",
    fill: COLORS.void,
    line: ctx.line("#00000000", 0),
  });
}

function panel(slide, ctx, x, y, w, h, color = COLORS.panel, line = COLORS.rule) {
  ctx.addShape(slide, {
    x,
    y,
    w,
    h,
    geometry: "rect",
    fill: color,
    line: ctx.line(line, 1),
  });
}

function label(slide, ctx, value, x, y, w, color, size = 13, bold = false) {
  text(slide, ctx, {
    text: value,
    x,
    y,
    w,
    h: 24,
    fontSize: size,
    bold,
    color,
  });
}

function card(slide, ctx, item, x, y, w, h, accent) {
  panel(slide, ctx, x, y, w, h, COLORS.panel, accent);
  text(slide, ctx, {
    text: item.kicker,
    x: x + 18,
    y: y + 16,
    w: w - 36,
    h: 18,
    fontSize: 10,
    bold: true,
    color: accent,
  });
  text(slide, ctx, {
    text: item.title,
    x: x + 18,
    y: y + 42,
    w: w - 36,
    h: 28,
    fontSize: 17,
    bold: true,
    color: COLORS.ink,
  });
  text(slide, ctx, {
    text: item.body,
    x: x + 18,
    y: y + 78,
    w: w - 36,
    h: h - 94,
    fontSize: 11,
    color: COLORS.soft,
  });
}

function fontForTitle(title, layout) {
  if (layout === "opening") return title.length > 54 ? 35 : 42;
  if (layout === "core") return title.length > 42 ? 27 : 32;
  return title.length > 48 ? 29 : 34;
}

function drawNativeScene(slide, ctx, spec) {
  if (spec.native === "v7_boundary") {
    const items = [
      { kicker: "COMPLETED", title: "v1-v7 benchmark", body: "8 tissues, public OSDR sources, multi-method evaluation, foundation-model comparisons." },
      { kicker: "PATCH", title: "v7.1 boundary", body: "Public counts, subset labels, license/citation language, and safe claim wording." },
      { kicker: "SEPARATE", title: "v8/v9 tracks", body: "Translation and platform layers remain explicitly downstream of completed benchmark evidence." },
    ];
    items.forEach((item, i) => card(slide, ctx, item, 90 + i * 380, 214, 320, 214, spec.accent));
    rule(slide, ctx, 150, 500, 978, spec.accent);
    label(slide, ctx, "public-facing claim set stays stable before extension material appears", 312, 522, 656, COLORS.soft, 14, true);
    return;
  }

  if (spec.native === "single_cell_blocker") {
    const items = [
      { kicker: "READY", title: "Manifest", body: "OSD-918 blood task contract exists with 8 samples, 4 Flight, 4 Ground." },
      { kicker: "READY", title: "Metric spec", body: "genelab_sc metric formulas, required inputs, and skip policies are defined." },
      { kicker: "BLOCKED", title: "Payload", body: "No processed h5ad, STARsolo bundle, or processed checksum manifest is listed." },
    ];
    items.forEach((item, i) => card(slide, ctx, item, 76 + i * 282, 202, 248, 236, i === 2 ? "#B45A7E" : spec.accent));
    panel(slide, ctx, 930, 202, 270, 236, COLORS.panel2, spec.accent);
    label(slide, ctx, "OSDR OSD-918 file list", 954, 226, 222, spec.accent, 12, true);
    label(slide, ctx, "19 files listed", 954, 268, 210, COLORS.ink, 20, true);
    label(slide, ctx, "16 raw FASTQ", 954, 310, 210, COLORS.sky, 17, true);
    label(slide, ctx, "8/8 blood SRX pairs", 954, 348, 210, COLORS.green, 17, true);
    label(slide, ctx, "0 processed payloads", 954, 386, 210, "#B45A7E", 17, true);
    return;
  }

  if (spec.native === "claim_matrix") {
    const items = [
      { kicker: "CLAIM", title: "Completed result core", body: "v1-v7 benchmark evidence, bounded model comparisons, and hardened result spine." },
      { kicker: "ALLOW", title: "Metadata alpha", body: "v9 task/source/provenance metadata can be reviewed with payload blockers visible." },
      { kicker: "BOUND", title: "Hypothesis layer", body: "v8 translation is hypothesis generation, not a countermeasure recommendation." },
      { kicker: "BLOCK", title: "Single-cell score", body: "No RRRM leaderboard or legacy-score promotion before payload audit passes." },
    ];
    items.forEach((item, i) => card(slide, ctx, item, 62 + i * 304, 196, 260, 246, [COLORS.green, COLORS.sky, "#B45A7E", COLORS.amber][i]));
    return;
  }

  if (spec.native === "roadmap") {
    const items = [
      { kicker: "1", title: "Payload freeze", body: "Hash-verify public bulk payloads and resolve RRRM single-cell source path." },
      { kicker: "2", title: "QA loop", body: "Render deck, verify text/source boundaries, and rerun figure accessibility checks." },
      { kicker: "3", title: "Release package", body: "Finalize Data Package, dataset card, provenance records, and archive split." },
      { kicker: "4", title: "Manuscript", body: "Align figures, claims, methods text, and supplement tables." },
    ];
    items.forEach((item, i) => {
      card(slide, ctx, item, 62 + i * 304, 214, 260, 214, spec.accent);
      if (i < 3) rule(slide, ctx, 324 + i * 304, 318, 36, spec.accent);
    });
  }
}

export async function addShellSlide(presentation, ctx, spec) {
  const slide = presentation.slides.add();
  ctx.addShape(slide, {
    x: 0,
    y: 0,
    w: ctx.W,
    h: ctx.H,
    geometry: "rect",
    fill: COLORS.void,
    line: ctx.line("#00000000", 0),
  });
  if (spec.asset) {
    if (spec.assetFit === "contain") {
      panel(slide, ctx, 94, 178, 1092, 414, COLORS.panel, spec.accent);
      await ctx.addImage(slide, {
        path: spec.asset,
        x: 116,
        y: 198,
        w: 1048,
        h: 374,
        fit: "contain",
        alt: `${spec.title} proof image`,
        name: `proof_image_${spec.slide}`,
      });
    } else {
      await ctx.addImage(slide, {
        path: spec.asset,
        x: 0,
        y: 0,
        w: ctx.W,
        h: ctx.H,
        fit: "cover",
        alt: `${spec.title} audited scene plate`,
        name: `scene_plate_${spec.slide}`,
      });
    }
  } else {
    drawNativeScene(slide, ctx, spec);
  }

  matte(slide, ctx, 0, spec.overlay.layout === "opening" ? 178 : spec.overlay.layout === "core" ? 154 : 170);
  matte(slide, ctx, 624, 96);

  rule(slide, ctx, 58, spec.overlay.layout === "opening" ? 158 : spec.overlay.layout === "core" ? 136 : 150, 388, spec.accent);
  text(slide, ctx, {
    text: `${String(spec.slide).padStart(2, "0")} / ${spec.overlay.eyebrow}`,
    x: 58,
    y: 30,
    w: 440,
    h: 22,
    fontSize: 11,
    bold: true,
    color: spec.accent,
    name: `eyebrow_${spec.slide}`,
  });

  text(slide, ctx, {
    text: spec.overlay.title,
    x: 58,
    y: spec.overlay.layout === "opening" ? 54 : spec.overlay.layout === "core" ? 52 : 58,
    w: spec.overlay.layout === "opening" ? 790 : 720,
    h: spec.overlay.layout === "opening" ? 56 : spec.overlay.layout === "core" ? 38 : 42,
    fontSize: fontForTitle(spec.overlay.title, spec.overlay.layout),
    bold: true,
    color: COLORS.ink,
    typeface: "Aptos Display",
    name: `title_${spec.slide}`,
  });

  text(slide, ctx, {
    text: spec.overlay.subtitle,
    x: 60,
    y: spec.overlay.layout === "opening" ? 121 : spec.overlay.layout === "core" ? 96 : 108,
    w: spec.overlay.layout === "opening" ? 740 : 690,
    h: spec.overlay.layout === "core" ? 38 : 40,
    fontSize: spec.overlay.layout === "core" ? 13 : 15,
    color: COLORS.soft,
    name: `subtitle_${spec.slide}`,
  });

  if (spec.overlay.bridge) {
    text(slide, ctx, {
      text: spec.overlay.bridge,
      x: 858,
      y: spec.overlay.layout === "opening" ? 52 : 45,
      w: 318,
      h: 58,
      fontSize: spec.overlay.layout === "opening" ? 17 : 15,
      color: COLORS.ink,
      bold: spec.overlay.layout === "opening",
      name: `bridge_${spec.slide}`,
    });
    rule(slide, ctx, 858, spec.overlay.layout === "opening" ? 114 : 104, 230, spec.accent);
  }

  text(slide, ctx, {
    text: spec.overlay.caveat,
    x: 58,
    y: 641,
    w: 922,
    h: 30,
    fontSize: 11.5,
    color: COLORS.ink,
    name: `claim_boundary_${spec.slide}`,
  });

  text(slide, ctx, {
    text: spec.overlay.source,
    x: 58,
    y: 674,
    w: 1058,
    h: 20,
    fontSize: 9.2,
    color: COLORS.soft,
    name: `source_${spec.slide}`,
  });

  if (slide.speakerNotes?.setText) {
    slide.speakerNotes.setText(spec.notes || "");
  }

  return slide;
}
'''


def write_slide_modules(deck_spec: list[dict]) -> None:
    if SLIDES_DIR.exists():
        shutil.rmtree(SLIDES_DIR)
    SLIDES_DIR.mkdir(parents=True, exist_ok=True)
    (SLIDES_DIR / "shared.mjs").write_text(SHARED_MODULE.strip() + "\n", encoding="utf-8")
    for slide in deck_spec:
        module = f"""import {{ addShellSlide }} from "./shared.mjs";

const spec = {js(slide)};

export async function slide{slide['slide']:02d}(presentation, ctx) {{
  return addShellSlide(presentation, ctx, spec);
}}
"""
        (SLIDES_DIR / f"slide-{slide['slide']:02d}.mjs").write_text(module, encoding="utf-8")


def run_build() -> dict:
    WORKSPACE_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PREVIEW_DIR.mkdir(parents=True, exist_ok=True)
    LAYOUT_DIR.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(NODE),
        str(BUILD_ARTIFACT_DECK),
        "--workspace",
        str(WORKSPACE),
        "--slides-dir",
        str(SLIDES_DIR),
        "--out",
        str(FINAL_PPTX),
        "--preview-dir",
        str(PREVIEW_DIR),
        "--layout-dir",
        str(LAYOUT_DIR),
        "--contact-sheet",
        str(CONTACT_SHEET),
        "--manifest",
        str(BUILD_MANIFEST),
        "--slide-count",
        "24",
        "--slide-size",
        "1280x720",
        "--scale",
        "1",
    ]
    result = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, check=True)
    return json.loads(result.stdout)


def copy_final_deliverables(build_manifest: dict, deck_spec: list[dict]) -> dict:
    PROJECT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    copied = {
        "pptx": PROJECT_OUTPUT_DIR / "spacebiobench_full_talk_deck_v0_3.pptx",
        "contact_sheet": PROJECT_OUTPUT_DIR / "spacebiobench_full_talk_deck_contact_sheet.png",
        "artifact_manifest": PROJECT_OUTPUT_DIR / "spacebiobench_full_talk_deck_artifact_manifest.json",
        "speaker_notes": PROJECT_OUTPUT_DIR / "spacebiobench_full_talk_deck_speaker_notes.md",
    }
    shutil.copy2(FINAL_PPTX, copied["pptx"])
    shutil.copy2(CONTACT_SHEET, copied["contact_sheet"])
    shutil.copy2(BUILD_MANIFEST, copied["artifact_manifest"])

    note_lines = ["# SpaceBio-Bench Full-Talk Deck Speaker Notes", "", f"Date: {CREATED}", ""]
    for slide in deck_spec:
        note_lines.extend([f"## Slide {slide['slide']}: {slide['title']}", "", slide["notes"], ""])
    copied["speaker_notes"].write_text("\n".join(note_lines), encoding="utf-8")

    production_manifest = {
        "created": CREATED,
        "artifact_role": "24-slide SpaceBio-Bench full-talk editable PPTX draft",
        "workspace": rel(WORKSPACE),
        "presentation_skill": "artifact-tool presentation build",
        "decision_source": "output/slides_1_14_pptx_readiness_audit_v0_1/slides_1_14_pptx_readiness_audit.json",
        "outputs": {key: rel(path) for key, path in copied.items()},
        "workspace_outputs": {
            "pptx": rel(FINAL_PPTX),
            "contact_sheet": rel(CONTACT_SHEET),
            "artifact_manifest": rel(BUILD_MANIFEST),
            "preview_dir": rel(PREVIEW_DIR),
            "layout_dir": rel(LAYOUT_DIR),
        },
        "qa_summary": {
            "slide_count": build_manifest.get("slideCount"),
            "pptx_bytes": build_manifest.get("outputBytes"),
            "editable_shell": "title, subtitle, bridge, source, claim boundary, and speaker notes are built as artifact-tool text objects",
            "scene_plate_rule": "audited scene plates are used instead of rendered previews",
            "full_talk_extension": "slides 15-24 add v7.1 boundary, v8 hypothesis incubator, v9 platform/public-bulk alpha, organoid/OSD-120/single-cell extensions, claim boundary, and roadmap",
            "known_status": "full-talk production draft, not final institutional-template deck",
        },
    }
    production_path = PROJECT_OUTPUT_DIR / "spacebiobench_full_talk_deck_manifest.json"
    production_path.write_text(json.dumps(production_manifest, indent=2), encoding="utf-8")
    copied["production_manifest"] = production_path
    return production_manifest


def main() -> None:
    deck_spec = build_deck_spec()
    write_workspace_notes(deck_spec)
    write_slide_modules(deck_spec)
    build_manifest = run_build()
    production_manifest = copy_final_deliverables(build_manifest, deck_spec)
    print(json.dumps(production_manifest, indent=2))


if __name__ == "__main__":
    main()
