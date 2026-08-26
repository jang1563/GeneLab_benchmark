#!/usr/bin/env python3
"""Build the SpaceBio-Bench v0.4 depth-master seminar deck.

This script preserves the v0.3 visual system but turns the methods and platform
sections into a teaching-first 30-slide deck. It imports the v0.3 builder for
shared paths, source parsing, and existing extension specs, then writes a new
artifact-tool workspace and output package.
"""

from __future__ import annotations

from copy import deepcopy
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import build_spacebiobench_full_talk_deck as base


ROOT = base.ROOT
CREATED = "2026-06-03"
THREAD_ID = os.environ.get("CODEX_THREAD_ID", "manual-20260603-spacebiobench-depth-master-v0-4")
TASK_SLUG = "spacebiobench-depth-master-deck"
WORKSPACE = ROOT / "outputs" / THREAD_ID / "presentations" / TASK_SLUG
SLIDES_DIR = WORKSPACE / "slides"
PREVIEW_DIR = WORKSPACE / "preview"
LAYOUT_DIR = WORKSPACE / "layout" / "final"
QA_DIR = WORKSPACE / "qa"
WORKSPACE_OUTPUT_DIR = WORKSPACE / "output"
PROJECT_OUTPUT_DIR = ROOT / "output" / "spacebiobench_depth_master_deck_v0_4"

FINAL_PPTX = WORKSPACE_OUTPUT_DIR / "spacebiobench_depth_master_deck_v0_4.pptx"
CONTACT_SHEET = PREVIEW_DIR / "spacebiobench_depth_master_deck_contact_sheet.png"
BUILD_MANIFEST = WORKSPACE_OUTPUT_DIR / "spacebiobench_depth_master_deck_artifact_manifest.json"


DEPTH_SCENE_HELPERS = r'''
function flowNode(slide, ctx, item, x, y, w, h, accent) {
  panel(slide, ctx, x, y, w, h, COLORS.panel, accent);
  text(slide, ctx, {
    text: item.kicker,
    x: x + 12,
    y: y + 12,
    w: w - 24,
    h: 16,
    fontSize: 9,
    bold: true,
    color: accent,
  });
  text(slide, ctx, {
    text: item.title,
    x: x + 12,
    y: y + 34,
    w: w - 24,
    h: 28,
    fontSize: 15,
    bold: true,
    color: COLORS.ink,
  });
  text(slide, ctx, {
    text: item.body,
    x: x + 12,
    y: y + 68,
    w: w - 24,
    h: h - 78,
    fontSize: 9.8,
    color: COLORS.soft,
  });
}

function connector(slide, ctx, x1, y, x2, color) {
  rule(slide, ctx, x1, y, x2 - x1, color);
  text(slide, ctx, { text: ">", x: x2 - 8, y: y - 11, w: 16, h: 18, fontSize: 12, bold: true, color });
}

function smallTerm(slide, ctx, item, x, y, w, accent) {
  panel(slide, ctx, x, y, w, 88, COLORS.panel, accent);
  text(slide, ctx, {
    text: item.term,
    x: x + 14,
    y: y + 13,
    w: w - 28,
    h: 22,
    fontSize: 15,
    bold: true,
    color: COLORS.ink,
  });
  text(slide, ctx, {
    text: item.def,
    x: x + 14,
    y: y + 42,
    w: w - 28,
    h: 38,
    fontSize: 9.8,
    color: COLORS.soft,
  });
}

function drawDataStrip(slide, ctx, spec) {
  const fields = (spec.dataCard || []).slice(0, 5);
  if (!fields.length) return;
  const x = 58;
  const y = 586;
  const gap = 10;
  const w = (1164 - gap * (fields.length - 1)) / fields.length;
  fields.forEach((item, i) => {
    const px = x + i * (w + gap);
    panel(slide, ctx, px, y, w, 32, COLORS.void, spec.accent);
    text(slide, ctx, {
      text: item.k,
      x: px + 9,
      y: y + 6,
      w: 48,
      h: 14,
      fontSize: 8.5,
      bold: true,
      color: spec.accent,
    });
    text(slide, ctx, {
      text: item.v,
      x: px + 58,
      y: y + 6,
      w: w - 66,
      h: 20,
      fontSize: 8.8,
      color: COLORS.ink,
    });
  });
}

function drawDataOrientation(slide, ctx, spec) {
  const nodes = [
    { kicker: "1", title: "Repository", body: "GeneLab / OSDR public space biology records." },
    { kicker: "2", title: "Study", body: "Experiment context, assay, organism, and metadata." },
    { kicker: "3", title: "Mission", body: "Flight or analog context used for split boundaries." },
    { kicker: "4", title: "Tissue", body: "Biological context; transfer can differ by organ." },
    { kicker: "5", title: "Sample", body: "One measured biospecimen with labels attached." },
    { kicker: "6", title: "Matrix", body: "Samples by molecular features after processing." },
    { kicker: "7", title: "Task", body: "Auditable prediction unit with metric and provenance." },
  ];
  nodes.forEach((node, i) => {
    const x = 44 + i * 172;
    flowNode(slide, ctx, node, x, 218, 142, 150, spec.accent);
    if (i < nodes.length - 1) connector(slide, ctx, x + 144, 292, x + 166, spec.accent);
  });
  panel(slide, ctx, 148, 430, 984, 82, COLORS.panel2, spec.accent);
  label(slide, ctx, "Teaching rule", 174, 452, 130, spec.accent, 11, true);
  text(slide, ctx, {
    text: "Every score shown later should be traceable back to this chain: source context, biological context, sample matrix, task definition, and metric.",
    x: 302,
    y: 448,
    w: 790,
    h: 44,
    fontSize: 15,
    color: COLORS.ink,
  });
}

function drawTaskRecord(slide, ctx, spec) {
  panel(slide, ctx, 74, 196, 500, 350, COLORS.panel, spec.accent);
  label(slide, ctx, "Example task contract", 100, 220, 260, spec.accent, 12, true);
  const rows = [
    ["Source", "OSDR / GeneLab study record"],
    ["Biology", "mission, organism, tissue, assay"],
    ["Samples", "sample IDs plus condition labels"],
    ["Split", "mission-held-out train/test rule"],
    ["Features", "gene matrix, pathway scores, or model input"],
    ["Metric", "AUROC or task-specific benchmark score"],
    ["Provenance", "task ID, data version, run record, caveat"],
  ];
  rows.forEach((row, i) => {
    const y = 260 + i * 36;
    text(slide, ctx, { text: row[0], x: 104, y, w: 92, h: 20, fontSize: 11, bold: true, color: spec.accent });
    text(slide, ctx, { text: row[1], x: 206, y, w: 330, h: 20, fontSize: 11.5, color: COLORS.ink });
  });
  const cards = [
    { kicker: "WHY", title: "A task is not a row in a table", body: "It is a source-bound contract. The model score is interpretable only when the split, label, feature view, and metric travel together." },
    { kicker: "AUDIT", title: "The contract prevents drift", body: "Later platform slides preserve this contract with manifests, evaluator outputs, checksums, and run records." },
  ];
  cards.forEach((item, i) => card(slide, ctx, item, 660, 220 + i * 168, 430, 136, i === 0 ? COLORS.sky : COLORS.amber));
}

function drawMissionProtocol(slide, ctx, spec) {
  panel(slide, ctx, 70, 214, 558, 230, COLORS.panel, COLORS.green);
  label(slide, ctx, "Train-side missions", 100, 238, 210, COLORS.green, 12, true);
  ["M1", "M2", "M3"].forEach((m, i) => {
    panel(slide, ctx, 110 + i * 146, 292, 100, 70, COLORS.panel2, COLORS.green);
    label(slide, ctx, m, 142 + i * 146, 314, 48, COLORS.ink, 20, true);
  });
  text(slide, ctx, { text: "preprocessing, feature choices, model selection, and thresholds are fixed here", x: 112, y: 382, w: 430, h: 34, fontSize: 13, color: COLORS.soft });
  ctx.addShape(slide, { x: 650, y: 190, w: 4, h: 306, geometry: "rect", fill: COLORS.amber, line: ctx.line("#00000000", 0) });
  text(slide, ctx, { text: "freeze boundary", x: 608, y: 504, w: 120, h: 20, fontSize: 12, bold: true, color: COLORS.amber });
  panel(slide, ctx, 712, 214, 392, 230, COLORS.panel, "#B45A7E");
  label(slide, ctx, "Hidden mission", 742, 238, 210, "#B45A7E", 12, true);
  panel(slide, ctx, 842, 292, 100, 70, COLORS.panel2, "#B45A7E");
  label(slide, ctx, "M4", 874, 314, 48, COLORS.ink, 20, true);
  text(slide, ctx, { text: "score once after the boundary; no tuning against the hidden mission", x: 760, y: 382, w: 306, h: 34, fontSize: 13, color: COLORS.soft });
  panel(slide, ctx, 250, 510, 780, 56, COLORS.panel2, spec.accent);
  text(slide, ctx, { text: "Method claim: the hidden test set is one whole mission context, with train-side choices frozen before scoring.", x: 274, y: 526, w: 730, h: 28, fontSize: 13, bold: true, color: COLORS.ink });
}

function drawFeaturePrimer(slide, ctx, spec) {
  const lanes = [
    { kicker: "VIEW A", title: "Gene matrix", body: "sample x gene expression values; high dimensional and noisy", accent: COLORS.sky },
    { kicker: "VIEW B", title: "Pathway scores", body: "genes summarized into biological programs before classification", accent: COLORS.green },
    { kicker: "VIEW C", title: "Model input", body: "embedding or transformed representation used by a classifier", accent: COLORS.amber },
  ];
  lanes.forEach((lane, i) => {
    const y = 202 + i * 120;
    panel(slide, ctx, 82, y, 280, 86, COLORS.panel, lane.accent);
    label(slide, ctx, lane.kicker, 104, y + 14, 90, lane.accent, 10, true);
    label(slide, ctx, lane.title, 104, y + 40, 200, COLORS.ink, 17, true);
    text(slide, ctx, { text: lane.body, x: 384, y: y + 20, w: 322, h: 42, fontSize: 12, color: COLORS.soft });
    connector(slide, ctx, 720, y + 43, 814, lane.accent);
  });
  panel(slide, ctx, 828, 274, 158, 132, COLORS.panel2, spec.accent);
  label(slide, ctx, "Same evaluator", 852, 310, 112, spec.accent, 12, true);
  label(slide, ctx, "one held-out score", 852, 350, 118, COLORS.ink, 16, true);
  connector(slide, ctx, 994, 340, 1078, spec.accent);
  panel(slide, ctx, 1094, 286, 110, 104, COLORS.panel, spec.accent);
  label(slide, ctx, "AUROC", 1116, 320, 80, COLORS.ink, 21, true);
  text(slide, ctx, { text: "Same samples, different feature views.", x: 414, y: 516, w: 456, h: 28, fontSize: 16, bold: true, color: COLORS.ink });
}

function drawScoreRecord(slide, ctx, spec) {
  const fields = [
    { kicker: "TASK", title: "task_id", body: "Which source-bound prediction unit?" },
    { kicker: "SPLIT", title: "held-out mission", body: "Which mission was hidden?" },
    { kicker: "MODEL", title: "classifier", body: "Which method was trained?" },
    { kicker: "FEATURE", title: "view", body: "Genes, pathways, or model input?" },
    { kicker: "METRIC", title: "AUROC", body: "How was separation scored?" },
    { kicker: "BOUND", title: "caveat", body: "What is not being claimed?" },
  ];
  fields.forEach((item, i) => flowNode(slide, ctx, item, 58 + i * 202, 228, 174, 154, spec.accent));
  panel(slide, ctx, 210, 454, 860, 72, COLORS.panel2, spec.accent);
  text(slide, ctx, {
    text: "A benchmark result is a score record, not just a number. If any field is missing, the audience cannot tell what the number means.",
    x: 242,
    y: 476,
    w: 800,
    h: 34,
    fontSize: 16,
    bold: true,
    color: COLORS.ink,
  });
}

function drawPlotPrimer(slide, ctx, spec) {
  panel(slide, ctx, 124, 202, 646, 342, COLORS.panel, spec.accent);
  text(slide, ctx, { text: "AUROC", x: 164, y: 224, w: 100, h: 24, fontSize: 15, bold: true, color: spec.accent });
  rule(slide, ctx, 230, 492, 440, COLORS.rule);
  [0.5, 0.7, 0.9].forEach((v, i) => {
    const x = 230 + i * 220;
    ctx.addShape(slide, { x, y: 260, w: 2, h: 232, geometry: "rect", fill: i === 0 ? COLORS.amber : COLORS.rule, line: ctx.line("#00000000", 0) });
    text(slide, ctx, { text: String(v), x: x - 12, y: 502, w: 40, h: 18, fontSize: 10, color: i === 0 ? COLORS.amber : COLORS.soft });
  });
  const rows = [
    ["Tissue A", 0.84, COLORS.green],
    ["Tissue B", 0.68, COLORS.sky],
    ["Tissue C", 0.51, COLORS.amber],
  ];
  rows.forEach((row, i) => {
    const y = 304 + i * 62;
    text(slide, ctx, { text: row[0], x: 164, y: y - 8, w: 80, h: 18, fontSize: 11, color: COLORS.ink });
    rule(slide, ctx, 316, y, 168, COLORS.muted);
    const dotX = 230 + (row[1] - 0.5) / 0.4 * 440;
    ctx.addShape(slide, { x: dotX - 5, y: y - 5, w: 10, h: 10, geometry: "ellipse", fill: row[2], line: ctx.line("#00000000", 0) });
    text(slide, ctx, { text: row[1].toFixed(2), x: dotX + 12, y: y - 10, w: 60, h: 18, fontSize: 10, bold: true, color: row[2] });
  });
  const legend = [
    { term: "0.5 line", def: "chance-level separation for a binary task" },
    { term: "Dot", def: "score for one tissue or task family" },
    { term: "Line", def: "uncertainty interval or variability cue" },
  ];
  legend.forEach((item, i) => smallTerm(slide, ctx, item, 828, 216 + i * 106, 314, i === 0 ? COLORS.amber : spec.accent));
}

function drawV8Transition(slide, ctx, spec) {
  const items = [
    { kicker: "DONE", title: "Benchmark evidence", body: "v1-v7 asks whether models transfer across space-biology tasks." },
    { kicker: "THEN", title: "Hypothesis incubator", body: "v8 uses signals to prioritize translational questions, not interventions." },
    { kicker: "BOUNDARY", title: "No operational claim", body: "No Mars-risk estimate, treatment, clinical, or countermeasure recommendation." },
  ];
  items.forEach((item, i) => {
    const accent = [COLORS.green, COLORS.sky, "#B45A7E"][i];
    card(slide, ctx, item, 104 + i * 360, 240, 292, 204, accent);
    if (i < 2) connector(slide, ctx, 398 + i * 360, 342, 456 + i * 360, accent);
  });
  panel(slide, ctx, 286, 502, 708, 42, COLORS.panel2, spec.accent);
  text(slide, ctx, { text: "Transition rule: translation appears only after the benchmark evidence and keeps a separate claim status.", x: 310, y: 514, w: 660, h: 18, fontSize: 13, bold: true, color: COLORS.ink });
}

function drawPlatformGlossary(slide, ctx, spec) {
  const terms = [
    { term: "Manifest", def: "expected files, task records, and source paths" },
    { term: "Evaluator", def: "code path that computes benchmark metrics" },
    { term: "Run record", def: "stored output from one completed evaluation" },
    { term: "Checksum", def: "file identity check for a local payload" },
    { term: "Payload", def: "analysis-ready data object used for scoring" },
    { term: "Blocker", def: "missing or unverifiable object that prevents release claims" },
  ];
  terms.forEach((item, i) => smallTerm(slide, ctx, item, 82 + (i % 3) * 384, 214 + Math.floor(i / 3) * 130, 318, i < 3 ? COLORS.sky : COLORS.amber));
  panel(slide, ctx, 226, 508, 828, 42, COLORS.panel2, spec.accent);
  text(slide, ctx, { text: "Key distinction: metadata readiness is not the same as payload readiness.", x: 252, y: 520, w: 780, h: 18, fontSize: 14, bold: true, color: COLORS.ink });
}

function drawPayloadLadder(slide, ctx, spec) {
  const steps = [
    { kicker: "1", title: "Metadata parsed", body: "source records and expected checksums are known", accent: COLORS.green },
    { kicker: "2", title: "Payload mirrored", body: "analysis-ready files exist locally", accent: COLORS.amber },
    { kicker: "3", title: "Hash verified", body: "local files match expected identity", accent: COLORS.amber },
    { kicker: "4", title: "Evaluator run", body: "metrics are recomputed from frozen inputs", accent: COLORS.sky },
    { kicker: "5", title: "Release frozen", body: "DOI/data package and leaderboard claims allowed", accent: "#B45A7E" },
  ];
  steps.forEach((item, i) => {
    flowNode(slide, ctx, item, 66 + i * 238, 236, 188, 168, item.accent);
    if (i < steps.length - 1) connector(slide, ctx, 254 + i * 238, 320, 296 + i * 238, item.accent);
  });
  panel(slide, ctx, 276, 474, 728, 48, COLORS.panel2, spec.accent);
  text(slide, ctx, { text: "v9 public bulk alpha is useful at the metadata layer, but release evidence requires payload and hash verification.", x: 302, y: 488, w: 682, h: 22, fontSize: 13.5, bold: true, color: COLORS.ink });
}
'''


DEPTH_SCENE_BRANCHES = r'''
  if (spec.native === "data_orientation") {
    drawDataOrientation(slide, ctx, spec);
    return;
  }
  if (spec.native === "task_record") {
    drawTaskRecord(slide, ctx, spec);
    return;
  }
  if (spec.native === "mission_protocol") {
    drawMissionProtocol(slide, ctx, spec);
    return;
  }
  if (spec.native === "feature_primer") {
    drawFeaturePrimer(slide, ctx, spec);
    return;
  }
  if (spec.native === "score_record") {
    drawScoreRecord(slide, ctx, spec);
    return;
  }
  if (spec.native === "plot_primer") {
    drawPlotPrimer(slide, ctx, spec);
    return;
  }
  if (spec.native === "v8_transition") {
    drawV8Transition(slide, ctx, spec);
    return;
  }
  if (spec.native === "platform_glossary") {
    drawPlatformGlossary(slide, ctx, spec);
    return;
  }
  if (spec.native === "payload_ladder") {
    drawPayloadLadder(slide, ctx, spec);
    return;
  }
'''


SHARED_MODULE = base.SHARED_MODULE.replace(
    "function drawNativeScene(slide, ctx, spec) {\n",
    DEPTH_SCENE_HELPERS + "\nfunction drawNativeScene(slide, ctx, spec) {\n" + DEPTH_SCENE_BRANCHES,
).replace(
    "  matte(slide, ctx, 0,",
    "  if (spec.dataCard) drawDataStrip(slide, ctx, spec);\n\n  matte(slide, ctx, 0,",
)


def native_slide(
    slide: int,
    section: str,
    title: str,
    role: str,
    native: str,
    overlay: dict,
    notes: str,
    accent_key: str = "Core methods",
) -> dict:
    return {
        "slide": slide,
        "section": section,
        "title": title,
        "role": role,
        "asset": "",
        "assetFit": "cover",
        "native": native,
        "accent": base.COLORS[base.SECTION_ACCENT.get(accent_key, "sky")],
        "overlay": overlay,
        "notes": notes,
        "watch": [],
    }


def audit_slide(row: dict, new_slide: int, data_card: list[dict] | None = None, bridge: str | None = None) -> dict:
    overlay = deepcopy(base.choose_overlay(row))
    if bridge is not None:
        overlay["bridge"] = bridge
    accent = base.COLORS[base.SECTION_ACCENT.get(row["section"], "teal")]
    spec = {
        "slide": new_slide,
        "section": row["section"],
        "title": row["title"],
        "role": row["role"],
        "asset": str((ROOT / row["scene_plate"]).resolve()),
        "assetFit": "cover",
        "native": "",
        "accent": accent,
        "overlay": overlay,
        "notes": base.speaker_notes(row, overlay),
        "watch": row.get("viewer_language_watch_hits", []),
    }
    if data_card:
        spec["dataCard"] = data_card
        spec["notes"] += "\nDepth pass: add on-slide data-card fields so first-time viewers can read the result contract."
    return spec


def extension_slide(row: dict, new_slide: int, data_card: list[dict] | None = None) -> dict:
    spec = {
        "slide": new_slide,
        "section": row["section"],
        "title": row["title"],
        "role": row["role"],
        "asset": str((ROOT / row["asset"]).resolve()) if row.get("asset") else "",
        "assetFit": row.get("assetFit", "cover"),
        "native": row.get("native", ""),
        "accent": base.COLORS[base.SECTION_ACCENT.get(row["section"], "teal")],
        "overlay": deepcopy(row["overlay"]),
        "notes": re.sub(r"Slide \d+", f"Slide {new_slide}", row["notes"], count=1),
        "watch": [],
    }
    if data_card:
        spec["dataCard"] = data_card
        spec["notes"] += "\nDepth pass: visible data-card fields clarify the release or extension status."
    return spec


def depth_slides() -> list[dict]:
    method_accent = "Core methods"
    return [
        native_slide(
            4,
            "Core methods",
            "Data orientation",
            "teaching bridge",
            "data_orientation",
            {
                "eyebrow": "DATA ORIENTATION",
                "title": "Start with the data unit",
                "subtitle": "SpaceBio-Bench keeps repository, mission, tissue, sample, matrix, and task context attached.",
                "bridge": "A benchmark task is not just a matrix; it is a source-bound data contract.",
                "caveat": "Orientation slide only; no result or novelty claim.",
                "source": "GeneLab/OSDR source-task inventory + depth audit",
                "layout": "method",
            },
            "Slide 4: data orientation\nDefine the basic data chain before methods or results: repository, study, mission, tissue, sample, matrix, task.",
            method_accent,
        ),
        native_slide(
            5,
            "Core methods",
            "Task record anatomy",
            "teaching bridge",
            "task_record",
            {
                "eyebrow": "TASK ANATOMY",
                "title": "One task is a contract",
                "subtitle": "Every score needs the source, biology, samples, split, features, metric, and provenance.",
                "bridge": "This is what makes later results auditable instead of anecdotal.",
                "caveat": "Task-contract explanation only; not a performance result.",
                "source": "source/task inventory + benchmark run-record design",
                "layout": "method",
            },
            "Slide 5: task record anatomy\nUse a concrete task contract so first-time viewers can see what one benchmark record contains.",
            method_accent,
        ),
        native_slide(
            6,
            "Core methods",
            "Mission-held-out protocol",
            "teaching bridge",
            "mission_protocol",
            {
                "eyebrow": "METHODS GUARD",
                "title": "Held-out means a hidden mission",
                "subtitle": "Model and preprocessing choices are frozen before the hidden mission is scored.",
                "bridge": "This is stricter than random sample validation because the hidden context is an entire mission.",
                "caveat": "Mission-held-out protocol; do not call it random cross-validation.",
                "source": "B3/B4 narration pack + mission-held-out split rule",
                "layout": "method",
            },
            "Slide 6: mission-held-out protocol\nExplain train-side choices, freeze boundary, hidden mission, and score-once rule.",
            method_accent,
        ),
        native_slide(
            7,
            "Core methods",
            "Feature-view primer",
            "teaching bridge",
            "feature_primer",
            {
                "eyebrow": "FEATURE VIEWS",
                "title": "Same samples, different views",
                "subtitle": "Genes, pathways, and transformed model inputs can be evaluated under the same held-out rule.",
                "bridge": "Define the feature view before reading any score surface.",
                "caveat": "Feature-view bridge only; not a quantitative result.",
                "source": "feature-layer bridge overlay spec + v1-v7 benchmark core",
                "layout": "method",
            },
            "Slide 7: feature-view primer\nBridge gene-level matrices, pathway summaries, and model inputs before showing tissue and model comparisons.",
            method_accent,
        ),
        native_slide(
            8,
            "Core methods",
            "Score record anatomy",
            "teaching bridge",
            "score_record",
            {
                "eyebrow": "SCORE RECORD",
                "title": "A score is only readable with its contract",
                "subtitle": "Task ID, split, model, feature view, metric, and caveat travel together.",
                "bridge": "This is the unit that connects methods to result slides.",
                "caveat": "Score-record explanation only; no leaderboard claim.",
                "source": "benchmark run-record schema + depth audit",
                "layout": "method",
            },
            "Slide 8: score record anatomy\nTeach that a result is a complete score record, not just a number.",
            method_accent,
        ),
        native_slide(
            9,
            "Core methods",
            "Plot-reading primer",
            "teaching bridge",
            "plot_primer",
            {
                "eyebrow": "PLOT PRIMER",
                "title": "Read the first result before judging it",
                "subtitle": "AUROC, chance line, score dot, uncertainty line, and tissue row are the basic grammar.",
                "bridge": "A strong result in one tissue is not a universal transfer claim.",
                "caveat": "Plot-reading primer; schematic values are illustrative.",
                "source": "v1-v7 AUROC result grammar + depth audit",
                "layout": "method",
            },
            "Slide 9: plot-reading primer\nDefine AUROC 0.5 as chance, dot as score, line as uncertainty, and row as tissue or task family.",
            method_accent,
        ),
    ]


RESULT_DATA_CARDS = {
    7: [
        {"k": "DATA", "v": "mouse bulk"},
        {"k": "SPLIT", "v": "mission-held-out"},
        {"k": "UNIT", "v": "tissue task"},
        {"k": "METRIC", "v": "AUROC"},
    ],
    8: [
        {"k": "DATA", "v": "selected tasks"},
        {"k": "VIEW", "v": "pathway summaries"},
        {"k": "TEST", "v": "held-out score"},
        {"k": "CLAIM", "v": "selected rescue"},
    ],
    9: [
        {"k": "DATA", "v": "shared 6-task panel"},
        {"k": "MODEL", "v": "tested settings"},
        {"k": "METRIC", "v": "AUROC / delta"},
        {"k": "CLAIM", "v": "no scale-only win"},
    ],
    10: [
        {"k": "GRID", "v": "8 tissues"},
        {"k": "MODELS", "v": "8 classifiers"},
        {"k": "VIEWS", "v": "4 feature sets"},
        {"k": "TOTAL", "v": "256 evaluations"},
    ],
    11: [
        {"k": "DATA", "v": "RR temporal / RRRM"},
        {"k": "READ", "v": "guardrail"},
        {"k": "STATUS", "v": "underpowered pilot"},
        {"k": "CLAIM", "v": "descriptive only"},
    ],
    12: [
        {"k": "DATA", "v": "negative controls"},
        {"k": "TEST", "v": "task limits"},
        {"k": "MODELS", "v": "tested settings"},
        {"k": "CLAIM", "v": "boundary, not absence"},
    ],
    13: [
        {"k": "INPUT", "v": "benchmark hits"},
        {"k": "OUTPUT", "v": "hypothesis triage"},
        {"k": "BIOLOGY", "v": "pathway / target tier"},
        {"k": "CLAIM", "v": "no mechanism proof"},
    ],
    14: [
        {"k": "SOURCE", "v": "mouse signal"},
        {"k": "COMPARE", "v": "human perturbation"},
        {"k": "READ", "v": "partial alignment"},
        {"k": "CLAIM", "v": "no clinical transfer"},
    ],
}


EXTENSION_DATA_CARDS = {
    18: [
        {"k": "ROLE", "v": "result boundary"},
        {"k": "STATUS", "v": "documentation patch"},
        {"k": "CORE", "v": "v1-v7 evidence"},
        {"k": "EXCLUDES", "v": "v8/v9 claims"},
    ],
    20: [
        {"k": "ROLE", "v": "hypothesis layer"},
        {"k": "INPUT", "v": "benchmark signals"},
        {"k": "STATUS", "v": "incubator"},
        {"k": "EXCLUDES", "v": "intervention claim"},
    ],
    21: [
        {"k": "ROLE", "v": "boundary"},
        {"k": "STATUS", "v": "provenance beta"},
        {"k": "CLAIM", "v": "not operational"},
        {"k": "EXCLUDES", "v": "countermeasure"},
    ],
    23: [
        {"k": "ROLE", "v": "platform"},
        {"k": "OBJECTS", "v": "manifests / runs"},
        {"k": "CLAIM", "v": "reproducibility"},
        {"k": "EXCLUDES", "v": "leaderboard"},
    ],
    24: [
        {"k": "STATUS", "v": "metadata alpha"},
        {"k": "SOURCES", "v": "22/22 parsed"},
        {"k": "PAYLOAD", "v": "0/22 hash-verified"},
        {"k": "CLAIM", "v": "not frozen"},
    ],
    26: [
        {"k": "TRACK", "v": "organoid"},
        {"k": "OBJECT", "v": "local matrices"},
        {"k": "STATUS", "v": "draft extension"},
        {"k": "EXCLUDES", "v": "core benchmark"},
    ],
    27: [
        {"k": "TRACK", "v": "OSD-120"},
        {"k": "SPLIT", "v": "same-study"},
        {"k": "ROLE", "v": "diagnostic"},
        {"k": "EXCLUDES", "v": "mission-held-out"},
    ],
}


def build_deck_spec() -> list[dict]:
    audit = json.loads(base.AUDIT_JSON.read_text(encoding="utf-8"))
    rows = {row["slide"]: row for row in audit["slides"]}
    extensions = {row["slide"]: row for row in base.EXTENSION_SLIDES}

    deck: list[dict] = [
        audit_slide(rows[1], 1),
        audit_slide(rows[2], 2, bridge="First define the niche, then define the data unit."),
        audit_slide(rows[3], 3, bridge="The talk separates benchmark evidence from translation and platform layers."),
        *depth_slides(),
    ]

    result_bridges = {
        7: "Read as tissue-specific transfer, not universal transfer.",
        8: "Read as selected pathway rescue, not universal superiority.",
        9: "Read as tested-setting comparison, not a universal model ranking.",
        10: "Read as coverage and robustness, not a new leaderboard.",
        11: "Read as guardrail evidence; recovery and RRRM remain cautious.",
        12: "Read as task-limit evidence, not absence of biology.",
        13: "Read as hypothesis triage, not mechanism proof.",
        14: "Read as partial translation, not clinical transfer.",
    }
    for old_slide in range(7, 15):
        deck.append(
            audit_slide(
                rows[old_slide],
                old_slide + 3,
                data_card=RESULT_DATA_CARDS[old_slide],
                bridge=result_bridges[old_slide],
            )
        )

    deck.append(extension_slide(extensions[15], 18, EXTENSION_DATA_CARDS[18]))
    deck.append(
        native_slide(
            19,
            "v8 incubator",
            "Translation transition",
            "teaching bridge",
            "v8_transition",
            {
                "eyebrow": "V8 TRANSITION",
                "title": "Translation appears after evidence",
                "subtitle": "The completed benchmark can motivate hypotheses without turning them into recommendations.",
                "bridge": "This transition protects the completed result spine from downstream overclaiming.",
                "caveat": "Transition slide only; no clinical, Mars-risk, or countermeasure claim.",
                "source": "docs/CANONICAL_RESULTS_V7_1.md + docs/V8_BETA_RELEASE_PLAN_2026_05_10.md",
                "layout": "extension",
            },
            "Slide 19: v8 transition\nExplain why translation appears after benchmark evidence and why it has a separate claim status.",
            "v8 incubator",
        )
    )
    deck.append(extension_slide(extensions[16], 20, EXTENSION_DATA_CARDS[20]))
    deck.append(extension_slide(extensions[17], 21, EXTENSION_DATA_CARDS[21]))
    deck.append(
        native_slide(
            22,
            "v9 platform",
            "Platform glossary",
            "teaching bridge",
            "platform_glossary",
            {
                "eyebrow": "V9 GLOSSARY",
                "title": "Platform terms before platform claims",
                "subtitle": "Manifest, evaluator, run record, checksum, payload, and blocker mean different things.",
                "bridge": "Define infrastructure terms before asking the audience to judge release readiness.",
                "caveat": "Glossary slide only; not a frozen-release claim.",
                "source": "v9 platform provenance notes + depth audit",
                "layout": "extension",
            },
            "Slide 22: platform glossary\nDefine manifest, evaluator, run record, checksum, payload, and blocker before the v9 platform slides.",
            "v9 platform",
        )
    )
    deck.append(extension_slide(extensions[18], 23, EXTENSION_DATA_CARDS[23]))
    deck.append(extension_slide(extensions[19], 24, EXTENSION_DATA_CARDS[24]))
    deck.append(
        native_slide(
            25,
            "v9 platform",
            "Payload readiness ladder",
            "teaching bridge",
            "payload_ladder",
            {
                "eyebrow": "RELEASE READINESS",
                "title": "Metadata is not payload readiness",
                "subtitle": "Parsed source records are useful, but frozen release evidence needs local payload and hash verification.",
                "bridge": "This ladder explains why alpha sections stay alpha.",
                "caveat": "Readiness framework only; no claim that blockers are resolved.",
                "source": "docs/V9_PUBLIC_BULK_ALPHA_METADATA_SNAPSHOT_DECISION.md + depth audit",
                "layout": "extension",
            },
            "Slide 25: payload readiness ladder\nTeach the difference between metadata parsed, payload mirrored, hash verified, evaluator run, and release frozen.",
            "v9 platform",
        )
    )
    deck.append(extension_slide(extensions[20], 26, EXTENSION_DATA_CARDS[26]))
    deck.append(extension_slide(extensions[21], 27, EXTENSION_DATA_CARDS[27]))
    deck.append(extension_slide(extensions[22], 28))
    deck.append(extension_slide(extensions[23], 29))
    deck.append(extension_slide(extensions[24], 30))
    return deck


def write_workspace_notes(deck_spec: list[dict]) -> None:
    WORKSPACE.mkdir(parents=True, exist_ok=True)
    QA_DIR.mkdir(parents=True, exist_ok=True)
    profile = [
        "task mode: create / depth-master extension of v0.3 deck",
        "primary deck-profile: engineering-platform",
        "secondary profile gates: scientific pedagogy, appendix-heavy source discipline, claim-boundary discipline",
        "required proof objects: data orientation, task anatomy, held-out protocol, feature-view primer, score-record anatomy, plot primer, platform glossary",
        "source/asset requirements: preserve audited scene plates for result proof surfaces; build teaching slides as native editable scenes",
        "brand authenticity constraints: no external logos or invented identity marks",
        "profile-specific QA gates: first-time viewer can identify data unit, split rule, feature view, metric, and release blocker status",
        "known missing inputs: no final institutional template supplied; this is a seminar depth master",
    ]
    (WORKSPACE / "profile-plan.txt").write_text("\n".join(profile) + "\n", encoding="utf-8")

    source_lines = [
        "Source notes",
        f"Date: {CREATED}",
        "",
        "Primary local source: output/slides_1_14_pptx_readiness_audit_v0_1/slides_1_14_pptx_readiness_audit.json",
        "Depth audit source: docs/SPACEBIOBENCH_FULL_TALK_DECK_DEPTH_AUDIT_2026_06_03.md",
        "Native teaching slides are editable artifact-tool scenes. Existing result proof surfaces use audited local scene plates.",
        "",
    ]
    for slide in deck_spec:
        asset = slide.get("asset")
        source_lines.append(f"Slide {slide['slide']}: {base.rel(Path(asset)) if asset else 'native editable teaching scene'}")
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
        "1-3 opening and positioning",
        "4-9 teaching-first methods and metric scaffolding",
        "10-17 result spine with data-card strips",
        "18-21 v7/v8 boundary and translation transition",
        "22-25 v9 platform glossary and readiness ladder",
        "26-28 extension tracks and blocker status",
        "29-30 claim matrix and roadmap",
        "",
        "Depth rule: every non-obvious term should either be defined before use or carried in a visible data/status card.",
    ]
    (WORKSPACE / "contact-sheet-plan.txt").write_text("\n".join(contact_lines) + "\n", encoding="utf-8")


def write_slide_modules(deck_spec: list[dict]) -> None:
    if SLIDES_DIR.exists():
        shutil.rmtree(SLIDES_DIR)
    SLIDES_DIR.mkdir(parents=True, exist_ok=True)
    (SLIDES_DIR / "shared.mjs").write_text(SHARED_MODULE.strip() + "\n", encoding="utf-8")
    for slide in deck_spec:
        module = f"""import {{ addShellSlide }} from "./shared.mjs";

const spec = {base.js(slide)};

export async function slide{slide['slide']:02d}(presentation, ctx) {{
  return addShellSlide(presentation, ctx, spec);
}}
"""
        (SLIDES_DIR / f"slide-{slide['slide']:02d}.mjs").write_text(module, encoding="utf-8")


def run_build(deck_spec: list[dict]) -> dict:
    WORKSPACE_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PREVIEW_DIR.mkdir(parents=True, exist_ok=True)
    LAYOUT_DIR.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(base.NODE),
        str(base.BUILD_ARTIFACT_DECK),
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
        str(len(deck_spec)),
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
        "pptx": PROJECT_OUTPUT_DIR / "spacebiobench_depth_master_deck_v0_4.pptx",
        "contact_sheet": PROJECT_OUTPUT_DIR / "spacebiobench_depth_master_deck_contact_sheet.png",
        "artifact_manifest": PROJECT_OUTPUT_DIR / "spacebiobench_depth_master_deck_artifact_manifest.json",
        "speaker_notes": PROJECT_OUTPUT_DIR / "spacebiobench_depth_master_deck_speaker_notes.md",
    }
    shutil.copy2(FINAL_PPTX, copied["pptx"])
    shutil.copy2(CONTACT_SHEET, copied["contact_sheet"])
    shutil.copy2(BUILD_MANIFEST, copied["artifact_manifest"])

    note_lines = ["# SpaceBio-Bench Depth-Master Deck Speaker Notes", "", f"Date: {CREATED}", ""]
    for slide in deck_spec:
        note_lines.extend([f"## Slide {slide['slide']}: {slide['title']}", "", slide["notes"], ""])
    copied["speaker_notes"].write_text("\n".join(note_lines), encoding="utf-8")

    production_manifest = {
        "created": CREATED,
        "artifact_role": "30-slide SpaceBio-Bench v0.4 depth-master seminar deck",
        "workspace": base.rel(WORKSPACE),
        "presentation_skill": "artifact-tool presentation build",
        "decision_source": "docs/SPACEBIOBENCH_FULL_TALK_DECK_DEPTH_AUDIT_2026_06_03.md",
        "outputs": {key: base.rel(path) for key, path in copied.items()},
        "workspace_outputs": {
            "pptx": base.rel(FINAL_PPTX),
            "contact_sheet": base.rel(CONTACT_SHEET),
            "artifact_manifest": base.rel(BUILD_MANIFEST),
            "preview_dir": base.rel(PREVIEW_DIR),
            "layout_dir": base.rel(LAYOUT_DIR),
        },
        "qa_summary": {
            "slide_count": build_manifest.get("slideCount"),
            "pptx_bytes": build_manifest.get("outputBytes"),
            "depth_additions": "six teaching methods slides plus v8 transition, v9 glossary, and payload readiness ladder",
            "data_card_rule": "result and platform proof slides carry compact visible data/status cards",
            "known_status": "seminar depth master, not final 24-slide conference cut-down or institutional-template deck",
        },
    }
    production_path = PROJECT_OUTPUT_DIR / "spacebiobench_depth_master_deck_manifest.json"
    production_path.write_text(json.dumps(production_manifest, indent=2), encoding="utf-8")
    copied["production_manifest"] = production_path
    return production_manifest


def main() -> None:
    deck_spec = build_deck_spec()
    write_workspace_notes(deck_spec)
    write_slide_modules(deck_spec)
    build_manifest = run_build(deck_spec)
    production_manifest = copy_final_deliverables(build_manifest, deck_spec)
    print(json.dumps(production_manifest, indent=2))


if __name__ == "__main__":
    main()
