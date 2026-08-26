#!/usr/bin/env python3
"""Build the SpaceBio-Bench conference cut-down deck.

This derives a shorter conference deck from the v0.4 depth master. It preserves
the new teaching scaffolding, but compresses the methods and platform readiness
sections so the talk can fit a shorter slot without dropping the data-unit,
held-out, feature-view, metric, and release-boundary explanations.
"""

from __future__ import annotations

from copy import deepcopy
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import build_spacebiobench_depth_master_deck as depth
import build_spacebiobench_full_talk_deck as base


ROOT = base.ROOT
CREATED = "2026-06-03"
THREAD_ID = os.environ.get("CODEX_THREAD_ID", "manual-20260603-spacebiobench-conference-v0-5")
TASK_SLUG = "spacebiobench-conference-deck"
WORKSPACE = ROOT / "outputs" / THREAD_ID / "presentations" / TASK_SLUG
SLIDES_DIR = WORKSPACE / "slides"
PREVIEW_DIR = WORKSPACE / "preview"
LAYOUT_DIR = WORKSPACE / "layout" / "final"
QA_DIR = WORKSPACE / "qa"
WORKSPACE_OUTPUT_DIR = WORKSPACE / "output"
PROJECT_OUTPUT_DIR = ROOT / "output" / "spacebiobench_conference_deck_v0_5"

FINAL_PPTX = WORKSPACE_OUTPUT_DIR / "spacebiobench_conference_deck_v0_5.pptx"
CONTACT_SHEET = PREVIEW_DIR / "spacebiobench_conference_deck_contact_sheet.png"
BUILD_MANIFEST = WORKSPACE_OUTPUT_DIR / "spacebiobench_conference_deck_artifact_manifest.json"


CONFERENCE_SCENE_HELPERS = r'''
function drawTaskContractCompressed(slide, ctx, spec) {
  const nodes = [
    { kicker: "1", title: "Repository", body: "GeneLab / OSDR public records." },
    { kicker: "2", title: "Study", body: "Experiment, organism, assay, metadata." },
    { kicker: "3", title: "Mission", body: "Context used for split boundaries." },
    { kicker: "4", title: "Sample", body: "Measured biospecimen with labels." },
    { kicker: "5", title: "Matrix", body: "Samples by molecular features." },
    { kicker: "6", title: "Task", body: "Prediction contract with metric." },
  ];
  nodes.forEach((node, i) => {
    const x = 52 + i * 190;
    flowNode(slide, ctx, node, x, 206, 156, 134, spec.accent);
    if (i < nodes.length - 1) connector(slide, ctx, x + 158, 270, x + 184, spec.accent);
  });

  panel(slide, ctx, 76, 404, 514, 126, COLORS.panel, spec.accent);
  label(slide, ctx, "Task contract fields", 100, 424, 220, spec.accent, 12, true);
  const left = [["Source", "study record"], ["Biology", "mission / tissue / assay"], ["Samples", "IDs + condition labels"]];
  const right = [["Split", "held-out rule"], ["Features", "gene / pathway / model view"], ["Metric", "AUROC or task score"]];
  left.forEach((row, i) => {
    text(slide, ctx, { text: row[0], x: 104, y: 458 + i * 28, w: 82, h: 18, fontSize: 10, bold: true, color: spec.accent });
    text(slide, ctx, { text: row[1], x: 190, y: 458 + i * 28, w: 170, h: 18, fontSize: 10.5, color: COLORS.ink });
  });
  right.forEach((row, i) => {
    text(slide, ctx, { text: row[0], x: 366, y: 458 + i * 28, w: 82, h: 18, fontSize: 10, bold: true, color: spec.accent });
    text(slide, ctx, { text: row[1], x: 450, y: 458 + i * 28, w: 128, h: 18, fontSize: 10.5, color: COLORS.ink });
  });

  panel(slide, ctx, 660, 404, 460, 126, COLORS.panel2, COLORS.amber);
  label(slide, ctx, "Conference shortcut", 684, 424, 220, COLORS.amber, 12, true);
  text(slide, ctx, {
    text: "When later result slides show a data strip, read it as the compact version of this contract.",
    x: 684,
    y: 458,
    w: 396,
    h: 44,
    fontSize: 15,
    bold: true,
    color: COLORS.ink,
  });
}

function drawFeatureScoreCompressed(slide, ctx, spec) {
  const lanes = [
    { label: "Gene matrix", detail: "sample x gene expression", color: COLORS.sky },
    { label: "Pathway scores", detail: "biological program summaries", color: COLORS.green },
    { label: "Compressed input", detail: "numeric model input", color: COLORS.amber },
  ];
  lanes.forEach((lane, i) => {
    const y = 214 + i * 82;
    panel(slide, ctx, 78, y, 246, 58, COLORS.panel, lane.color);
    label(slide, ctx, lane.label, 98, y + 16, 190, COLORS.ink, 15, true);
    text(slide, ctx, { text: lane.detail, x: 350, y: y + 16, w: 214, h: 24, fontSize: 10.8, color: COLORS.soft });
    connector(slide, ctx, 586, y + 30, 656, lane.color);
  });
  panel(slide, ctx, 676, 252, 160, 116, COLORS.panel2, spec.accent);
  label(slide, ctx, "Evaluator", 708, 284, 108, spec.accent, 12, true);
  label(slide, ctx, "one held-out score", 708, 324, 126, COLORS.ink, 15, true);
  connector(slide, ctx, 850, 310, 914, spec.accent);
  panel(slide, ctx, 930, 214, 220, 208, COLORS.panel, spec.accent);
  label(slide, ctx, "AUROC primer", 954, 236, 130, spec.accent, 12, true);
  rule(slide, ctx, 974, 342, 108, COLORS.rule);
  ctx.addShape(slide, { x: 990, y: 278, w: 2, h: 64, geometry: "rect", fill: COLORS.amber, line: ctx.line("#00000000", 0) });
  ctx.addShape(slide, { x: 1060, y: 304, w: 10, h: 10, geometry: "ellipse", fill: COLORS.green, line: ctx.line("#00000000", 0) });
  text(slide, ctx, { text: "0.5 = chance", x: 954, y: 362, w: 112, h: 16, fontSize: 9.6, color: COLORS.amber });
  text(slide, ctx, { text: "dot = score", x: 954, y: 384, w: 100, h: 16, fontSize: 9.6, color: COLORS.ink });
  text(slide, ctx, { text: "line = uncertainty", x: 954, y: 406, w: 126, h: 16, fontSize: 9.6, color: COLORS.soft });
  panel(slide, ctx, 238, 486, 804, 44, COLORS.panel2, spec.accent);
  text(slide, ctx, {
    text: "Read every result as: task contract + feature view + held-out score + caveat.",
    x: 270,
    y: 500,
    w: 748,
    h: 18,
    fontSize: 14,
    bold: true,
    color: COLORS.ink,
  });
}

function drawModelViewBridge(slide, ctx, spec) {
  panel(slide, ctx, 70, 214, 250, 176, COLORS.panel, COLORS.green);
  label(slide, ctx, "BIOLOGY", 94, 236, 120, COLORS.green, 11, true);
  label(slide, ctx, "Sample context", 94, 272, 180, COLORS.ink, 19, true);
  text(slide, ctx, {
    text: "mission, tissue, assay, condition labels",
    x: 94,
    y: 312,
    w: 182,
    h: 46,
    fontSize: 12.4,
    color: COLORS.soft,
  });
  connector(slide, ctx, 326, 302, 382, COLORS.green);
  text(slide, ctx, { text: "processed into", x: 328, y: 320, w: 92, h: 18, fontSize: 9.6, color: COLORS.muted });

  panel(slide, ctx, 404, 184, 426, 266, COLORS.panel2, spec.accent);
  label(slide, ctx, "NUMERICAL VIEWS", 432, 206, 180, spec.accent, 11, true);
  const views = [
    { title: "Gene matrix", body: "sample x gene expression values", color: COLORS.sky },
    { title: "Pathway scores", body: "genes summarized into biological programs", color: COLORS.green },
    { title: "Compressed input", body: "numeric summary used by the model", color: COLORS.amber },
  ];
  views.forEach((view, i) => {
    const y = 246 + i * 62;
    panel(slide, ctx, 432, y, 350, 48, COLORS.panel, view.color);
    text(slide, ctx, { text: view.title, x: 450, y: y + 9, w: 130, h: 18, fontSize: 13.5, bold: true, color: COLORS.ink });
    text(slide, ctx, { text: view.body, x: 590, y: y + 10, w: 172, h: 18, fontSize: 9.8, color: COLORS.soft });
  });
  connector(slide, ctx, 848, 302, 904, spec.accent);
  text(slide, ctx, { text: "fed to", x: 862, y: 320, w: 58, h: 18, fontSize: 9.6, color: COLORS.muted });

  panel(slide, ctx, 926, 214, 286, 176, COLORS.panel, COLORS.amber);
  label(slide, ctx, "EVALUATOR", 950, 236, 132, COLORS.amber, 11, true);
  label(slide, ctx, "Train, then score", 950, 272, 204, COLORS.ink, 19, true);
  text(slide, ctx, {
    text: "learn on known missions; score the hidden mission once",
    x: 950,
    y: 312,
    w: 214,
    h: 46,
    fontSize: 12.4,
    color: COLORS.soft,
  });

  panel(slide, ctx, 150, 500, 980, 56, COLORS.panel2, spec.accent);
  text(slide, ctx, {
    text: "Reading rule: a score is not a general AI/biology claim; it is one numerical view scored on one hidden task with a caveat.",
    x: 180,
    y: 516,
    w: 920,
    h: 24,
    fontSize: 14,
    bold: true,
    color: COLORS.ink,
  });
}

function drawModelComparisonSetup(slide, ctx, spec) {
  panel(slide, ctx, 70, 192, 226, 286, COLORS.panel, COLORS.sky);
  label(slide, ctx, "FIXED FIRST", 94, 216, 120, COLORS.sky, 11, true);
  const fixed = [
    { k: "split", v: "hidden mission" },
    { k: "view", v: "gene / pathway / compressed" },
    { k: "metric", v: "AUROC + delta" },
  ];
  fixed.forEach((row, i) => {
    const y = 264 + i * 58;
    text(slide, ctx, { text: row.k, x: 94, y, w: 62, h: 16, fontSize: 9.8, bold: true, color: COLORS.sky });
    text(slide, ctx, { text: row.v, x: 160, y, w: 150, h: 18, fontSize: 11.4, color: COLORS.ink });
    if (i < fixed.length - 1) rule(slide, ctx, 94, y + 36, 164, COLORS.rule);
  });
  connector(slide, ctx, 302, 334, 356, COLORS.sky);
  text(slide, ctx, { text: "then compare", x: 304, y: 352, w: 90, h: 18, fontSize: 9.4, color: COLORS.muted });

  const tiers = [
    {
      kicker: "CLASSICAL ML",
      title: "task-trained models",
      body: "PCA-LR / RF / LR use gene or pathway features",
      color: COLORS.sky,
    },
    {
      kicker: "FOUNDATION MODEL",
      title: "pretrained, then adapted",
      body: "scGPT / Geneformer start from expression pretraining",
      color: COLORS.green,
    },
    {
      kicker: "TEXT LLM CHECK",
      title: "prompt-only language model",
      body: "Gemini / Llama / DeepSeek see prompts, not matrices",
      color: COLORS.amber,
    },
  ];
  tiers.forEach((tier, i) => {
    const y = 174 + i * 122;
    panel(slide, ctx, 374, y, 382, 92, COLORS.panel2, tier.color);
    label(slide, ctx, tier.kicker, 398, y + 18, 172, tier.color, 10.2, true);
    text(slide, ctx, { text: tier.title, x: 398, y: y + 42, w: 260, h: 22, fontSize: 16, bold: true, color: COLORS.ink });
    text(slide, ctx, { text: tier.body, x: 398, y: y + 66, w: 300, h: 18, fontSize: 10.4, color: COLORS.soft });
    connector(slide, ctx, 762, y + 46, 838, tier.color);
  });

  panel(slide, ctx, 858, 222, 310, 216, COLORS.panel, spec.accent);
  label(slide, ctx, "ONE RESULT SURFACE", 884, 246, 190, spec.accent, 11, true);
  text(slide, ctx, {
    text: "Same task, same hidden mission, same metric",
    x: 884,
    y: 286,
    w: 226,
    h: 40,
    fontSize: 18,
    bold: true,
    color: COLORS.ink,
  });
  const read = [
    { k: "score", v: "AUROC" },
    { k: "baseline", v: "matched classical row" },
    { k: "claim", v: "tested-setting comparison" },
  ];
  read.forEach((row, i) => {
    text(slide, ctx, { text: row.k, x: 884, y: 350 + i * 26, w: 76, h: 16, fontSize: 9.6, bold: true, color: spec.accent });
    text(slide, ctx, { text: row.v, x: 970, y: 350 + i * 26, w: 172, h: 16, fontSize: 10.8, color: COLORS.soft });
  });

  panel(slide, ctx, 188, 552, 904, 44, COLORS.panel2, spec.accent);
  text(slide, ctx, {
    text: "Reading rule: this is a controlled model-family stress test, not a universal intelligence ranking.",
    x: 220,
    y: 564,
    w: 840,
    h: 20,
    fontSize: 13.8,
    bold: true,
    color: COLORS.ink,
  });
}

function drawResultReadingGuide(slide, ctx, spec) {
  const guide = spec.resultGuide;
  if (!guide) return;
  panel(slide, ctx, 76, 218, 286, 250, COLORS.panel, spec.accent);
  label(slide, ctx, guide.kicker || "READ", 100, 238, 120, spec.accent, 11, true);
  text(slide, ctx, {
    text: guide.title,
    x: 100,
    y: 266,
    w: 238,
    h: 30,
    fontSize: 17,
    bold: true,
    color: COLORS.ink,
  });
  (guide.rows || []).slice(0, 3).forEach((row, i) => {
    const y = 318 + i * 38;
    text(slide, ctx, { text: row.k, x: 100, y, w: 82, h: 16, fontSize: 9.8, bold: true, color: row.color || spec.accent });
    text(slide, ctx, { text: row.v, x: 184, y, w: 144, h: 20, fontSize: 10.4, color: COLORS.soft });
  });
  connector(slide, ctx, 368, 326, 514, spec.accent);
  text(slide, ctx, {
    text: guide.footer,
    x: 100,
    y: 424,
    w: 238,
    h: 30,
    fontSize: 10.6,
    bold: true,
    color: COLORS.ink,
  });
}

function drawGuardrailBadges(slide, ctx, spec) {
  const badges = spec.guardrailBadges || [];
  badges.forEach((badge) => {
    const y = badge.y || 164;
    panel(slide, ctx, badge.x, y, badge.w, 44, COLORS.panel2, badge.color || spec.accent);
    text(slide, ctx, {
      text: badge.kicker,
      x: badge.x + 14,
      y: y + 8,
      w: badge.w - 28,
      h: 14,
      fontSize: 8.8,
      bold: true,
      color: badge.color || spec.accent,
    });
    text(slide, ctx, {
      text: badge.title,
      x: badge.x + 14,
      y: y + 24,
      w: badge.w - 28,
      h: 16,
      fontSize: 10.6,
      bold: true,
      color: COLORS.ink,
    });
  });
}

function drawPlatformReadinessCompressed(slide, ctx, spec) {
  const terms = [
    { term: "Manifest", def: "expected files and task records" },
    { term: "Evaluator", def: "code path that computes metrics" },
    { term: "Run record", def: "stored evaluation output" },
    { term: "Checksum", def: "file identity check" },
  ];
  terms.forEach((item, i) => smallTerm(slide, ctx, item, 70 + i * 300, 206, 250, i < 2 ? COLORS.sky : COLORS.amber));

  const steps = [
    { title: "Metadata parsed", color: COLORS.green },
    { title: "Payload mirrored", color: COLORS.amber },
    { title: "Hash verified", color: COLORS.amber },
    { title: "Evaluator run", color: COLORS.sky },
    { title: "Release frozen", color: "#B45A7E" },
  ];
  steps.forEach((step, i) => {
    const x = 82 + i * 234;
    panel(slide, ctx, x, 388, 172, 92, COLORS.panel, step.color);
    text(slide, ctx, { text: String(i + 1), x: x + 14, y: 402, w: 22, h: 18, fontSize: 10, bold: true, color: step.color });
    text(slide, ctx, { text: step.title, x: x + 14, y: 432, w: 132, h: 28, fontSize: 14, bold: true, color: COLORS.ink });
    if (i < steps.length - 1) connector(slide, ctx, x + 174, 434, x + 220, step.color);
  });
  panel(slide, ctx, 250, 522, 780, 38, COLORS.panel2, spec.accent);
  text(slide, ctx, {
    text: "The public-bulk alpha is useful at metadata readiness; release claims wait for payload and hash verification.",
    x: 274,
    y: 532,
    w: 730,
    h: 18,
    fontSize: 12.8,
    bold: true,
    color: COLORS.ink,
  });
}

function drawPlatformObjectsCompressed(slide, ctx, spec) {
  const objects = [
    {
      kicker: "1 MANIFEST",
      title: "Declare inputs",
      body: "expected files, task IDs, checksums, and source records",
      color: COLORS.sky,
    },
    {
      kicker: "2 EVALUATOR",
      title: "Run metric code",
      body: "the scoring path computes metrics from frozen task inputs",
      color: COLORS.green,
    },
    {
      kicker: "3 RUN RECORD",
      title: "Store the trace",
      body: "parameters, outputs, and provenance stay available for audit",
      color: COLORS.amber,
    },
  ];
  objects.forEach((item, i) => {
    const x = 82 + i * 356;
    panel(slide, ctx, x, 226, 276, 156, COLORS.panel, item.color);
    text(slide, ctx, {
      text: item.kicker,
      x: x + 18,
      y: 244,
      w: 130,
      h: 16,
      fontSize: 9.5,
      bold: true,
      color: item.color,
    });
    text(slide, ctx, {
      text: item.title,
      x: x + 18,
      y: 278,
      w: 218,
      h: 24,
      fontSize: 18,
      bold: true,
      color: COLORS.ink,
    });
    text(slide, ctx, {
      text: item.body,
      x: x + 18,
      y: 320,
      w: 226,
      h: 38,
      fontSize: 11,
      color: COLORS.soft,
    });
    if (i < objects.length - 1) {
      connector(slide, ctx, x + 278, 304, x + 346, item.color);
    }
  });

  panel(slide, ctx, 240, 438, 800, 68, COLORS.panel2, spec.accent);
  text(slide, ctx, {
    text: "Reproducibility surface",
    x: 270,
    y: 452,
    w: 230,
    h: 22,
    fontSize: 17,
    bold: true,
    color: COLORS.ink,
  });
  text(slide, ctx, {
    text: "These objects let others recompute and audit the benchmark; they do not add a new biological score or leaderboard.",
    x: 512,
    y: 454,
    w: 480,
    h: 28,
    fontSize: 12.4,
    bold: true,
    color: COLORS.soft,
  });
}

function drawMetadataAlphaStatus(slide, ctx, spec) {
  const states = [
    {
      kicker: "METADATA",
      value: "22 / 22",
      title: "sources parsed",
      body: "checksum-manifest evidence is present for every source record",
      color: COLORS.green,
    },
    {
      kicker: "PAYLOAD",
      value: "0 / 22",
      title: "hash-verified",
      body: "no local payload mirror has passed hash verification yet",
      color: "#B45A7E",
    },
    {
      kicker: "CLAIM",
      value: "NOT",
      title: "frozen release",
      body: "metadata snapshot only; no DOI release, score, or leaderboard claim",
      color: COLORS.amber,
    },
  ];
  states.forEach((item, i) => {
    const x = 94 + i * 364;
    panel(slide, ctx, x, 218, 300, 184, COLORS.panel, item.color);
    text(slide, ctx, {
      text: item.kicker,
      x: x + 20,
      y: 238,
      w: 120,
      h: 16,
      fontSize: 9.5,
      bold: true,
      color: item.color,
    });
    text(slide, ctx, {
      text: item.value,
      x: x + 20,
      y: 272,
      w: 138,
      h: 42,
      fontSize: 30,
      bold: true,
      color: item.color,
    });
    text(slide, ctx, {
      text: item.title,
      x: x + 168,
      y: 284,
      w: 108,
      h: 26,
      fontSize: 15,
      bold: true,
      color: COLORS.ink,
    });
    text(slide, ctx, {
      text: item.body,
      x: x + 20,
      y: 336,
      w: 246,
      h: 38,
      fontSize: 10.8,
      color: COLORS.soft,
    });
    if (i < states.length - 1) {
      connector(slide, ctx, x + 302, 310, x + 356, item.color);
    }
  });

  panel(slide, ctx, 230, 456, 820, 66, COLORS.panel2, spec.accent);
  text(slide, ctx, {
    text: "Allowed alpha",
    x: 260,
    y: 470,
    w: 164,
    h: 22,
    fontSize: 16,
    bold: true,
    color: COLORS.ink,
  });
  text(slide, ctx, {
    text: "Metadata can be reviewed now; release claims wait until payload hashes pass.",
    x: 440,
    y: 472,
    w: 560,
    h: 26,
    fontSize: 12.8,
    bold: true,
    color: COLORS.soft,
  });
}

function drawSingleCellBlockerReadable(slide, ctx, spec) {
  const items = [
    { kicker: "READY", title: "Manifest", body: "OSD-918 blood task contract exists with 8 samples, 4 Flight, 4 Ground.", color: spec.accent },
    { kicker: "READY", title: "Metric spec", body: "genelab_sc metric formulas, required inputs, and skip policies are defined.", color: spec.accent },
    { kicker: "BLOCKED", title: "Payload", body: "No processed h5ad, STARsolo bundle, or processed checksum manifest is listed.", color: "#B45A7E" },
  ];
  items.forEach((item, i) => card(slide, ctx, item, 76 + i * 282, 202, 248, 216, item.color));

  panel(slide, ctx, 930, 202, 270, 216, COLORS.panel2, spec.accent);
  label(slide, ctx, "OSDR OSD-918 file list", 954, 226, 222, spec.accent, 12, true);
  label(slide, ctx, "19 files listed", 954, 266, 210, COLORS.ink, 20, true);
  label(slide, ctx, "16 raw FASTQ", 954, 306, 210, COLORS.sky, 17, true);
  label(slide, ctx, "8/8 blood SRX pairs", 954, 344, 210, COLORS.green, 17, true);
  label(slide, ctx, "0 processed payloads", 954, 382, 210, "#B45A7E", 17, true);

  const steps = [
    { kicker: "RAW EXISTS", title: "source files", color: COLORS.sky },
    { kicker: "PAYLOAD MISSING", title: "no h5ad / STARsolo", color: "#B45A7E" },
    { kicker: "SCORE BLOCKED", title: "no public score", color: COLORS.amber },
  ];
  steps.forEach((step, i) => {
    const x = 164 + i * 340;
    panel(slide, ctx, x, 474, 276, 68, COLORS.panel2, step.color);
    label(slide, ctx, step.kicker, x + 22, 490, 220, step.color, 10, true);
    label(slide, ctx, step.title, x + 22, 514, 220, COLORS.ink, 16, true);
    if (i < steps.length - 1) connector(slide, ctx, x + 278, 508, x + 332, step.color);
  });
  text(slide, ctx, {
    text: "Raw FASTQ is source material; scoring waits for an analysis-ready processed payload.",
    x: 276,
    y: 558,
    w: 728,
    h: 22,
    fontSize: 13,
    bold: true,
    color: COLORS.soft,
  });
}

function drawCloseRoadmap(slide, ctx, spec) {
  const claims = [
    { kicker: "CLAIM", title: "Completed core", body: "v1-v7 benchmark evidence and bounded model comparisons." },
    { kicker: "ALLOW", title: "Metadata alpha", body: "v9 provenance can be reviewed while payload blockers stay visible." },
    { kicker: "BOUND", title: "Hypothesis layer", body: "v8 translation is not an intervention or countermeasure claim." },
    { kicker: "BLOCK", title: "Single-cell score", body: "No RRRM score promotion before processed payload audit passes." },
  ];
  claims.forEach((item, i) => card(slide, ctx, item, 58 + i * 304, 198, 260, 150, [COLORS.green, COLORS.sky, "#B45A7E", COLORS.amber][i]));

  const road = [
    { kicker: "1", title: "Freeze payloads", body: "hash-verify public bulk and resolve RRRM single-cell source path" },
    { kicker: "2", title: "QA package", body: "rerender deck, verify sources, rerun accessibility checks" },
    { kicker: "3", title: "Release + paper", body: "align Data Package, dataset card, claims, methods, supplement" },
  ];
  road.forEach((item, i) => {
    card(slide, ctx, item, 148 + i * 340, 404, 284, 126, spec.accent);
    if (i < 2) connector(slide, ctx, 434 + i * 340, 466, 486 + i * 340, spec.accent);
  });
}
'''


CONFERENCE_SCENE_BRANCHES = r'''
  if (spec.native === "task_contract_compressed") {
    drawTaskContractCompressed(slide, ctx, spec);
    return;
  }
  if (spec.native === "feature_score_compressed") {
    drawFeatureScoreCompressed(slide, ctx, spec);
    return;
  }
  if (spec.native === "model_view_bridge") {
    drawModelViewBridge(slide, ctx, spec);
    return;
  }
  if (spec.native === "model_comparison_setup") {
    drawModelComparisonSetup(slide, ctx, spec);
    return;
  }
  if (spec.native === "platform_readiness_compressed") {
    drawPlatformReadinessCompressed(slide, ctx, spec);
    return;
  }
  if (spec.native === "platform_objects_compressed") {
    drawPlatformObjectsCompressed(slide, ctx, spec);
    return;
  }
  if (spec.native === "metadata_alpha_status") {
    drawMetadataAlphaStatus(slide, ctx, spec);
    return;
  }
  if (spec.native === "single_cell_blocker") {
    drawSingleCellBlockerReadable(slide, ctx, spec);
    return;
  }
  if (spec.native === "close_roadmap") {
    drawCloseRoadmap(slide, ctx, spec);
    return;
  }
'''


SHARED_MODULE = depth.SHARED_MODULE.replace(
    "function drawNativeScene(slide, ctx, spec) {\n",
    CONFERENCE_SCENE_HELPERS + "\nfunction drawNativeScene(slide, ctx, spec) {\n" + CONFERENCE_SCENE_BRANCHES,
).replace(
    "  if (spec.dataCard) drawDataStrip(slide, ctx, spec);\n\n  matte(slide, ctx, 0,",
    "  if (spec.dataCard) drawDataStrip(slide, ctx, spec);\n  if (spec.resultGuide) drawResultReadingGuide(slide, ctx, spec);\n  if (spec.guardrailBadges) drawGuardrailBadges(slide, ctx, spec);\n\n  matte(slide, ctx, 0,",
).replace(
    """      await ctx.addImage(slide, {
        path: spec.asset,
        x: 0,
        y: 0,
        w: ctx.W,
        h: ctx.H,
        fit: "cover",""",
    """      await ctx.addImage(slide, {
        path: spec.asset,
        x: spec.assetXOffset || 0,
        y: spec.assetYOffset || 0,
        w: ctx.W,
        h: ctx.H,
        fit: "cover",""",
)


def native_slide(
    slide: int,
    section: str,
    title: str,
    role: str,
    native: str,
    overlay: dict,
    notes: str,
    accent_key: str,
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
    spec = depth.audit_slide(row, new_slide, data_card=data_card, bridge=bridge)
    spec["slide"] = new_slide
    return spec


def extension_slide(row: dict, new_slide: int, data_card: list[dict] | None = None, bridge: str | None = None) -> dict:
    spec = depth.extension_slide(row, new_slide, data_card=data_card)
    if bridge is not None:
        spec["overlay"]["bridge"] = bridge
    return spec


def with_track(spec: dict, track_key: int | str) -> dict:
    spec["track_key"] = track_key
    return spec


CONFERENCE_TALK_TRACK = {
    1: {
        "time20": "0:40",
        "time15": "0:25",
        "purpose": "Thesis",
        "talk": "The question is not whether AI can describe biology in-distribution. The question is whether a model trained on known space biology missions can survive a new mission context.",
        "cut": "State only the core question.",
        "guardrail": "Do not claim first AI benchmark for space biology.",
    },
    2: {
        "time20": "0:45",
        "time15": "0:30",
        "purpose": "Positioning",
        "talk": "This is a distinct niche, not a firstness claim. Existing resources matter, but the gap here is mission-held-out biological evaluation under spaceflight domain shift.",
        "cut": "Say distinct niche, not firstness; skip named landscape details.",
        "guardrail": "Do not overstate novelty against NASA BPS, GLARE, OpenProblems, or cell-eval.",
    },
    3: {
        "time20": "0:50",
        "time15": "0:35",
        "purpose": "Story map",
        "talk": "The talk has three layers: completed benchmark results, hypothesis-only translation, and platformization. Keep those layers separate throughout.",
        "cut": "Name the three evidence layers and move on.",
        "guardrail": "Do not present v1-v9 outputs as equal-strength scientific discoveries.",
    },
    4: {
        "time20": "1:05",
        "time15": "0:55",
        "purpose": "Data contract",
        "talk": "A benchmark task starts as a public record, but it becomes useful only when study, mission, sample, matrix, split, feature view, and metric stay attached as a contract.",
        "cut": "Explain the task contract; do not read all fields.",
        "guardrail": "Methods scaffold only; no result or novelty claim.",
    },
    5: {
        "time20": "0:55",
        "time15": "0:45",
        "purpose": "Held-out protocol",
        "talk": "Held-out here means an entire mission is hidden. Preprocessing, feature choices, and model selection are made on train-side missions before the hidden mission is scored.",
        "cut": "Define mission-held-out in one sentence.",
        "guardrail": "Do not call this random cross-validation.",
    },
    "model_view_bridge": {
        "time20": "0:50",
        "time15": "0:35",
        "purpose": "Plain-language ML bridge",
        "talk": "Before the first result, translate what the model actually receives. It does not see space biology directly; it sees numerical views of samples, such as a gene matrix, pathway scores, or a compressed input representation.",
        "cut": "Say that the model sees numerical views, not biology itself.",
        "guardrail": "Conceptual method bridge only; no model-performance or biological interpretation claim.",
    },
    6: {
        "time20": "1:00",
        "time15": "0:45",
        "purpose": "Feature and metric primer",
        "talk": "The same task can be scored from genes, pathways, or compressed model inputs. AUROC is the score grammar: 0.5 is chance, the dot is the score, and the interval is uncertainty.",
        "cut": "Define feature views and AUROC grammar; skip compressed-input detail.",
        "guardrail": "Illustrative score grammar only; do not imply a result on this slide.",
    },
    "model_compare_setup": {
        "time20": "0:55",
        "time15": "0:35",
        "purpose": "Model-family setup",
        "talk": "Before the model comparison result, fix the benchmark surface. Classical ML means task-trained models using gene or pathway features, foundation models start from expression pretraining and are adapted here, and text LLM checks are prompt-only language-model diagnostics. They are compared only after the hidden mission split, input view, and metric are fixed.",
        "cut": "Define the three model families, then say same split and same metric.",
        "guardrail": "Method comparison scaffold only; not a universal AI ranking.",
    },
    7: {
        "time20": "0:55",
        "time15": "0:40",
        "purpose": "First result",
        "talk": "Use this first result as the worked example. Read the guide first: 0.5 is chance, the dot is the tissue score, and the horizontal line is uncertainty. Some tissues sit to the right of chance, while others remain near chance.",
        "cut": "Use as the only full result-reading example: chance line, dot, uncertainty, tissue-specific claim.",
        "guardrail": "Do not imply one tissue result generalizes to all tissues.",
    },
    8: {
        "time20": "0:45",
        "time15": "0:30",
        "purpose": "Pathway rescue",
        "talk": "Use the visual guide first: the unwanted labels are mission, hardware, and gravity labels; pathway summaries sometimes make those labels weaker. That supports a selected rescue claim, not universal pathway superiority.",
        "cut": "Read the guide only: pathway summaries weaken selected unwanted labels.",
        "guardrail": "Do not claim pathway features always outperform gene-level models.",
    },
    9: {
        "time20": "0:45",
        "time15": "0:30",
        "purpose": "Model comparison",
        "talk": "Use the model guide first: the three families from the methods slide now appear on the same result surface. Classical PCA-LR remains a strong baseline, foundation-model gains are local rather than automatic, and text LLM checks stay near the tested-setting boundary. The claim is no scale-only win in these tested settings.",
        "cut": "Read the guide only: model families, strong baseline, local gains, no scale-only win.",
        "guardrail": "Do not universalize foundation-model failure beyond tested settings.",
    },
    10: {
        "time20": "0:55",
        "time15": "0:45",
        "purpose": "Hardening",
        "talk": "The hardening grid reduces cherry-pick risk: 8 tissues, 8 classifiers, and 4 feature views. This is coverage evidence, not a new ranking surface.",
        "cut": "Keep the 8 x 8 x 4 = 256 hardening claim.",
        "guardrail": "Hardening and coverage evidence only; not a new ranking claim.",
    },
    11: {
        "time20": "0:55",
        "time15": "0:35",
        "purpose": "Guardrails",
        "talk": "Read the guardrail badges left to right: preservation can dominate, recovery is descriptive only, and the RRRM pilot is ready but underpowered. These are useful signals, but none should be over-read.",
        "cut": "Read the badges only: preservation dominates, recovery is descriptive, RRRM is underpowered.",
        "guardrail": "Recovery projection is descriptive; RRRM composition is underpowered.",
    },
    12: {
        "time20": "0:45",
        "time15": "0:30",
        "purpose": "Negative boundary",
        "talk": "Read the boundary badges first: one control shows no score signal, the failure-mode panel marks the evidence boundary, and the model check shows no automatic gain. Negative results are useful because they define task limits; they do not prove absence of biology.",
        "cut": "Read the badges only: no score signal, boundary evidence, no automatic gain.",
        "guardrail": "No universal absence-of-biology or universal model-failure claim.",
    },
    13: {
        "time20": "0:55",
        "time15": "0:35",
        "purpose": "Biology interpretation",
        "talk": "Read the triage badges first: benchmark hits are signals to inspect, biological context helps prioritize follow-up, and the claim status remains hypothesis only. This is useful biology triage, not mechanism proof or treatment evidence.",
        "cut": "Read the badges only: hit, context, hypothesis only.",
        "guardrail": "Interpretation and target triage only; no treatment or mechanism proof.",
    },
    14: {
        "time20": "0:55",
        "time15": "0:40",
        "purpose": "Human translation",
        "talk": "Read the translation badges first: mouse evidence starts the comparison, pathway and target-tier alignment are partial, and the claim limit remains no clinical transfer. Weak direct gene transfer is why this stays a translation hypothesis rather than a clinical claim.",
        "cut": "Read the badges only: mouse signal, partial alignment, no clinical transfer.",
        "guardrail": "No direct gene-transfer, clinical, or countermeasure claim.",
    },
    15: {
        "time20": "0:40",
        "time15": "0:25",
        "purpose": "v7 boundary",
        "talk": "v7.1 is a canonical result surface. It freezes documentation and claim discipline; it is not a new benchmark run.",
        "cut": "v7.1 is documentation boundary.",
        "guardrail": "Do not mix v8 SpaceMed claims into the v7.1 benchmark paper claim.",
    },
    16: {
        "time20": "0:50",
        "time15": "0:35",
        "purpose": "v8 transition",
        "talk": "v8 uses benchmark signals as a hypothesis incubator. That is useful, but it remains downstream of the completed benchmark evidence.",
        "cut": "v8 is hypothesis incubator.",
        "guardrail": "Do not present v8 as a validated intervention surface.",
    },
    17: {
        "time20": "0:45",
        "time15": "0:25",
        "purpose": "Countermeasure boundary",
        "talk": "Read this as a claim boundary: stressor decomposition can make hypothesis evidence traceable, but it does not make the biology operational. There is no countermeasure recommendation here.",
        "cut": "No countermeasure recommendation.",
        "guardrail": "No Mars point estimate, treatment, or operational recommendation claim.",
    },
    18: {
        "time20": "0:55",
        "time15": "0:50",
        "purpose": "Platform readiness",
        "talk": "The platform terms matter because metadata readiness is not payload readiness. Release claims wait for payload mirroring, hash verification, and evaluator runs.",
        "cut": "Keep metadata versus payload readiness.",
        "guardrail": "Readiness framework only; do not claim blockers are resolved.",
    },
    19: {
        "time20": "0:40",
        "time15": "0:25",
        "purpose": "v9 platform",
        "talk": "Read the object chain left to right: manifests declare inputs, the evaluator runs metric code, and run records store the trace. Together they create a reproducibility surface, not a new biological score.",
        "cut": "Read the object chain only: manifest, evaluator, run record, reproducibility surface.",
        "guardrail": "Not a frozen release or ranking claim.",
    },
    20: {
        "time20": "0:45",
        "time15": "0:35",
        "purpose": "Public bulk alpha",
        "talk": "Read the three status cards: 22 of 22 source records are parsed, 0 of 22 payloads are hash-verified, and the release claim is not frozen. That is a useful metadata alpha, not frozen release evidence.",
        "cut": "Read the three cards only: 22 parsed, 0 hash-verified, not frozen.",
        "guardrail": "No frozen payload mirror, DOI release, or score claim.",
    },
    21: {
        "time20": "0:35",
        "time15": "0:25",
        "purpose": "Organoid extension",
        "talk": "Read the badges first: public OSDR source records become local sample-by-gene matrices, but this remains a draft organoid extension outside the public bulk core.",
        "cut": "Read badges only: source record, local matrix, draft boundary.",
        "guardrail": "Draft pilot extension only; not mission-held-out benchmark evidence.",
    },
    22: {
        "time20": "0:35",
        "time15": "0:25",
        "purpose": "OSD-120",
        "talk": "Read the badges first: OSD-120 starts from one public source study, the train/test split stays within that same study, and the claim limit is not mission-held-out generalization.",
        "cut": "Read badges only: source study, same-study split, not mission-held-out.",
        "guardrail": "Diagnostic extension only; not the core held-out transfer claim.",
    },
    23: {
        "time20": "0:45",
        "time15": "0:35",
        "purpose": "Single-cell blocker",
        "talk": "Read the bottom rule first: raw FASTQ source files exist, but the processed h5ad or STARsolo payload is missing, so the public metric remains blocked.",
        "cut": "Read bottom rule only: raw exists, payload missing, score blocked.",
        "guardrail": "Do not promote a score before the processed payload audit passes.",
    },
    24: {
        "time20": "1:10",
        "time15": "0:55",
        "purpose": "Close",
        "talk": "The closing rule is simple: completed core, metadata alpha, hypothesis layer, and blocked scores stay separate. The next work is payload freeze, QA, and release-plus-paper alignment.",
        "cut": "Close with separated claims and next work.",
        "guardrail": "No frozen v9 release, no v8 countermeasure, no RRRM score claim.",
    },
}


def enrich_speaker_notes(deck_spec: list[dict]) -> None:
    for slide in deck_spec:
        track = CONFERENCE_TALK_TRACK[slide.get("track_key", slide["slide"])]
        slide["notes"] = "\n".join(
            [
                f"Slide {slide['slide']}: {slide['overlay']['title']}",
                f"Purpose: {track['purpose']}",
                f"20-minute timing: {track['time20']}",
                f"15-minute cut timing: {track['time15']}",
                "",
                "20-minute talk track:",
                track["talk"],
                "",
                "15-minute cut:",
                track["cut"],
                "",
                "Claim guardrail:",
                track["guardrail"],
                "",
                "Presenter cue:",
                slide["overlay"].get("bridge", ""),
            ]
        )


def build_deck_spec() -> list[dict]:
    audit = json.loads(base.AUDIT_JSON.read_text(encoding="utf-8"))
    rows = {row["slide"]: row for row in audit["slides"]}
    extensions = {row["slide"]: row for row in base.EXTENSION_SLIDES}

    deck: list[dict] = [
        audit_slide(rows[1], 1),
        audit_slide(rows[2], 2, bridge="First define the niche, then define the evaluation unit."),
        audit_slide(rows[3], 3, bridge="The task separates benchmark evidence from translation and platform layers."),
        native_slide(
            4,
            "Core methods",
            "Data-to-task contract",
            "compressed teaching bridge",
            "task_contract_compressed",
            {
                "eyebrow": "DATA + TASK",
                "title": "Public data becomes a task contract",
                "subtitle": "Repository records stay attached to study, mission, sample, matrix, split, feature, and metric context.",
                "bridge": "This is the compressed version of the seminar data orientation.",
                "caveat": "Methods scaffold only; no result or novelty claim.",
                "source": "GeneLab/OSDR source-task inventory + v0.4 depth master",
                "layout": "method",
            },
            "Slide 4: data-to-task contract\nConference cut-down: combine data orientation and task anatomy into one visible contract.",
            "Core methods",
        ),
        native_slide(
            5,
            "Core methods",
            "Mission-held-out protocol",
            "compressed teaching bridge",
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
            "Slide 5: mission-held-out protocol\nKeep the train-side choice, freeze boundary, hidden mission, score-once rule.",
            "Core methods",
        ),
        with_track(
            native_slide(
                6,
                "Core methods",
                "Model-input bridge",
                "plain-language methods bridge",
                "model_view_bridge",
                {
                    "eyebrow": "MODEL VIEW",
                    "title": "What the model actually sees",
                    "subtitle": "AI models do not see space biology directly; the benchmark gives them numerical views of samples.",
                    "bridge": "Translate ML terms before the score appears.",
                    "caveat": "Conceptual method bridge only; no performance or biology claim.",
                    "source": "v0.4 feature-view primer + score-record grammar",
                    "layout": "method",
                },
                "Slide 6: model-input bridge\nTranslate ML terms into plain benchmark-reading language.",
                "Core methods",
            ),
            "model_view_bridge",
        ),
        with_track(
            native_slide(
                7,
                "Core methods",
                "Feature and score primer",
                "compressed teaching bridge",
                "feature_score_compressed",
                {
                    "eyebrow": "FEATURE + SCORE",
                    "title": "Feature views feed one held-out score",
                    "subtitle": "The same task can be scored from genes, pathway summaries, or compressed model inputs.",
                    "bridge": "Define the score grammar before the first result.",
                    "caveat": "Methods and plot-reading primer only; illustrative AUROC grammar.",
                    "source": "feature-layer bridge + v1-v7 AUROC result grammar + v0.4 depth master",
                    "layout": "method",
                },
                "Slide 7: feature and score primer\nCompress feature views, score record, and AUROC primer into one conference slide.",
                "Core methods",
            ),
            6,
        ),
        with_track(
            native_slide(
                8,
                "Core methods",
                "Model comparison setup",
                "foundation-model bridge",
                "model_comparison_setup",
                {
                    "eyebrow": "MODEL TIERS",
                    "title": "Model families share one task",
                    "subtitle": "Classical ML, foundation models, and text LLM checks are compared only after the split, input view, and metric are fixed.",
                    "bridge": "Define model tiers before the first model-comparison result.",
                    "caveat": "Method comparison scaffold only; no universal AI ranking claim.",
                    "source": "V1 model-tier methods + foundation-model comparison audit + v0.4 depth master",
                    "layout": "method",
                },
                "Slide 8: model comparison setup\nBridge feature views to classical ML, foundation models, and text LLM diagnostic checks.",
                "Core methods",
            ),
            "model_compare_setup",
        ),
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
        spec = audit_slide(
            rows[old_slide],
            old_slide + 2,
            data_card=depth.RESULT_DATA_CARDS[old_slide],
            bridge=result_bridges[old_slide],
        )
        if old_slide == 7:
            spec["resultGuide"] = {
                "kicker": "PLOT GUIDE",
                "title": "Read this first",
                "rows": [
                    {"k": "0.5", "v": "chance baseline", "color": base.COLORS["amber"]},
                    {"k": "dot", "v": "tissue score", "color": base.COLORS["green"]},
                    {"k": "line", "v": "uncertainty", "color": base.COLORS["sky"]},
                ],
                "footer": "Right of 0.5 suggests transfer; claim stays tissue-specific.",
            }
        if old_slide == 8:
            spec["assetYOffset"] = 16
            spec["resultGuide"] = {
                "kicker": "RESCUE GUIDE",
                "title": "Lower unwanted labels",
                "rows": [
                    {"k": "labels", "v": "mission / hardware / gravity", "color": base.COLORS["amber"]},
                    {"k": "view", "v": "pathway summaries", "color": base.COLORS["green"]},
                    {"k": "effect", "v": "lower nuisance scores", "color": base.COLORS["sky"]},
                ],
                "footer": "Rescue means weaker confounding labels, not better biology everywhere.",
            }
        if old_slide == 9:
            spec["resultGuide"] = {
                "kicker": "MODEL GUIDE",
                "title": "Compare locally",
                "rows": [
                    {"k": "families", "v": "classical / FM / text LLM", "color": base.COLORS["amber"]},
                    {"k": "baseline", "v": "PCA-LR remains strong", "color": base.COLORS["green"]},
                    {"k": "result", "v": "gains local, often negative", "color": base.COLORS["sky"]},
                ],
                "footer": "No scale-only win; not a universal ranking.",
            }
        if old_slide == 11:
            spec["guardrailBadges"] = [
                {
                    "x": 78,
                    "w": 304,
                    "kicker": "1 PRESERVATION",
                    "title": "can dominate",
                    "color": base.COLORS["amber"],
                },
                {
                    "x": 472,
                    "w": 304,
                    "kicker": "2 RECOVERY",
                    "title": "descriptive only",
                    "color": base.COLORS["sky"],
                },
                {
                    "x": 866,
                    "w": 304,
                    "kicker": "3 RRRM PILOT",
                    "title": "ready, underpowered",
                    "color": base.COLORS["rose"],
                },
            ]
        if old_slide == 12:
            spec["guardrailBadges"] = [
                {
                    "x": 78,
                    "w": 304,
                    "kicker": "1 SPATIAL CONTROL",
                    "title": "no score signal",
                    "color": base.COLORS["rose"],
                },
                {
                    "x": 504,
                    "w": 304,
                    "kicker": "2 FAILURE MODE",
                    "title": "boundary evidence",
                    "color": base.COLORS["green"],
                },
                {
                    "x": 866,
                    "w": 304,
                    "kicker": "3 MODEL CHECK",
                    "title": "no automatic gain",
                    "color": base.COLORS["sky"],
                },
            ]
        if old_slide == 13:
            spec["guardrailBadges"] = [
                {
                    "x": 78,
                    "w": 304,
                    "kicker": "1 BENCHMARK HIT",
                    "title": "signal to inspect",
                    "color": base.COLORS["rose"],
                },
                {
                    "x": 504,
                    "w": 304,
                    "kicker": "2 BIOLOGY CONTEXT",
                    "title": "prioritize follow-up",
                    "color": base.COLORS["green"],
                },
                {
                    "x": 866,
                    "w": 304,
                    "kicker": "3 CLAIM STATUS",
                    "title": "hypothesis only",
                    "color": base.COLORS["amber"],
                },
            ]
        if old_slide == 14:
            spec["guardrailBadges"] = [
                {
                    "x": 78,
                    "w": 214,
                    "kicker": "1 MOUSE SIGNAL",
                    "title": "starting evidence",
                    "color": base.COLORS["sky"],
                },
                {
                    "x": 504,
                    "w": 304,
                    "kicker": "2 PARTIAL ALIGNMENT",
                    "title": "pathway + tier only",
                    "color": base.COLORS["green"],
                },
                {
                    "x": 866,
                    "w": 304,
                    "kicker": "3 CLAIM LIMIT",
                    "title": "no clinical transfer",
                    "color": base.COLORS["amber"],
                },
            ]
        deck.append(
            with_track(
                spec,
                old_slide,
            )
        )

    deck.append(with_track(extension_slide(extensions[15], 17, depth.EXTENSION_DATA_CARDS[18]), 15))
    deck.append(
        with_track(
            extension_slide(
                extensions[16],
                18,
                depth.EXTENSION_DATA_CARDS[20],
                bridge="Translation is a bounded hypothesis layer downstream of benchmark evidence.",
            ),
            16,
        )
    )
    countermeasure_boundary = deepcopy(extensions[17])
    countermeasure_boundary["headline"] = "No countermeasure recommendation"
    countermeasure_boundary["asset"] = "output/spacebiobench_conference_deck_v0_5/assets/Figure2_Stressor_Decomposition_AB_only.png"
    countermeasure_boundary["overlay"] = deepcopy(countermeasure_boundary["overlay"])
    countermeasure_boundary["overlay"]["title"] = "No countermeasure recommendation"
    countermeasure_boundary["overlay"]["subtitle"] = "Stressors remain hypothesis signals before any release claim."
    deck.append(with_track(extension_slide(countermeasure_boundary, 19, depth.EXTENSION_DATA_CARDS[21]), 17))
    deck.append(
        with_track(
            native_slide(
                20,
                "v9 platform",
                "Platform readiness",
                "compressed teaching bridge",
                "platform_readiness_compressed",
                {
                    "eyebrow": "V9 READINESS",
                    "title": "Metadata is not payload readiness",
                    "subtitle": "Manifest and run-record terms matter because release claims require payload and hash verification.",
                    "bridge": "This compresses the glossary and readiness ladder into one release-boundary slide.",
                    "caveat": "Readiness framework only; no claim that blockers are resolved.",
                    "source": "v9 platform provenance + public bulk alpha decision + v0.4 depth master",
                    "layout": "extension",
                },
                "Slide 20: platform readiness\nCompress v9 glossary and payload-readiness ladder into one conference slide.",
                "v9 platform",
            ),
            18,
        )
    )
    platform_objects = native_slide(
        21,
        "v9 platform",
        "Platform object chain",
        "platform architecture",
        "platform_objects_compressed",
        {
            "eyebrow": "V9 PLATFORM",
            "title": "Manifests, evaluator, and run records",
            "subtitle": "v9 turns result artifacts into auditable benchmark infrastructure.",
            "bridge": "The platform layer is about reproducibility surfaces, not new biological claims.",
            "caveat": "Platform provenance surface; not a frozen release or leaderboard claim.",
            "source": "v9 platform provenance + manifests/evaluator/run-record design",
            "layout": "extension",
        },
        "Slide 21: platform object chain\nTranslate manifests, evaluator outputs, and run records into a reproducibility-surface diagram.",
        "v9 platform",
    )
    platform_objects["dataCard"] = depth.EXTENSION_DATA_CARDS[23]
    deck.append(with_track(platform_objects, 19))
    metadata_alpha = native_slide(
        22,
        "v9 platform",
        "Public bulk alpha",
        "release boundary",
        "metadata_alpha_status",
        {
            "eyebrow": "PUBLIC BULK ALPHA",
            "title": "Metadata alpha, payload hashes blocked",
            "subtitle": "22/22 sources have parsed checksum-manifest evidence; 0/22 payloads are locally hash-verified.",
            "bridge": "A metadata snapshot can be useful without pretending to be a frozen payload release.",
            "caveat": "Metadata-only alpha; no frozen payload mirror, DOI release, or leaderboard claim.",
            "source": "docs/V9_PUBLIC_BULK_ALPHA_METADATA_SNAPSHOT_DECISION.md",
            "layout": "extension",
        },
        "Slide 22: metadata alpha status\nShow parsed metadata, blocked payload hashes, and release boundary as large status cards.",
        "v9 platform",
    )
    metadata_alpha["dataCard"] = depth.EXTENSION_DATA_CARDS[24]
    deck.append(with_track(metadata_alpha, 20))
    organoid_bridge = extension_slide(extensions[20], 23, depth.EXTENSION_DATA_CARDS[26])
    organoid_bridge["assetYOffset"] = 26
    organoid_bridge["guardrailBadges"] = [
        {
            "x": 40,
            "y": 532,
            "w": 384,
            "kicker": "1 SOURCE RECORD",
            "title": "public OSDR page",
            "color": base.COLORS["sky"],
        },
        {
            "x": 424,
            "y": 532,
            "w": 410,
            "kicker": "2 LOCAL MATRIX",
            "title": "sample-by-gene table",
            "color": base.COLORS["green"],
        },
        {
            "x": 834,
            "y": 532,
            "w": 352,
            "kicker": "3 DRAFT BOUNDARY",
            "title": "outside core benchmark",
            "color": base.COLORS["amber"],
        },
    ]
    deck.append(with_track(organoid_bridge, 21))
    osd120_bridge = extension_slide(extensions[21], 24, depth.EXTENSION_DATA_CARDS[27])
    osd120_bridge["assetYOffset"] = 26
    osd120_bridge["guardrailBadges"] = [
        {
            "x": 40,
            "y": 532,
            "w": 384,
            "kicker": "1 SOURCE STUDY",
            "title": "public OSDR record",
            "color": base.COLORS["sky"],
        },
        {
            "x": 424,
            "y": 532,
            "w": 410,
            "kicker": "2 SAME-STUDY SPLIT",
            "title": "train/test inside study",
            "color": base.COLORS["green"],
        },
        {
            "x": 834,
            "y": 532,
            "w": 352,
            "kicker": "3 CLAIM LIMIT",
            "title": "not mission-held-out",
            "color": base.COLORS["amber"],
        },
    ]
    deck.append(with_track(osd120_bridge, 22))
    deck.append(with_track(extension_slide(extensions[22], 25), 23))
    deck.append(
        with_track(
            native_slide(
                26,
                "Close",
                "Claim boundary and roadmap",
                "conference close",
                "close_roadmap",
                {
                    "eyebrow": "CLOSE",
                    "title": "Separate claims, then freeze",
                    "subtitle": "Completed evidence, alpha metadata, hypothesis-only translation, and blocked scores stay separate.",
                    "bridge": "Close with next steps that match each claim status.",
                    "caveat": "No frozen v9 release, no v8 countermeasure, no RRRM score claim.",
                    "source": "v7.1 canonical boundary + v8/v9 release-boundary docs + release roadmap",
                    "layout": "extension",
                },
                "Slide 26: claim boundary and roadmap\nConference close: combine the claim matrix with concrete next steps.",
                "Close",
            ),
            24,
        )
    )
    enrich_speaker_notes(deck)
    return deck


def write_workspace_notes(deck_spec: list[dict]) -> None:
    WORKSPACE.mkdir(parents=True, exist_ok=True)
    QA_DIR.mkdir(parents=True, exist_ok=True)
    profile = [
        "task mode: targeted-edit / conference cut-down from v0.4 depth master",
        "primary deck-profile: engineering-platform",
        "secondary profile gates: scientific pedagogy; conference time-boxing; claim-boundary discipline",
        "required proof objects: compressed data-task contract, mission-held-out protocol, model-view bridge, feature/score primer, result data cards, platform readiness boundary",
        "source/asset requirements: preserve audited scene plates for result proof surfaces; use native editable scenes for compressed explainers",
        "brand authenticity constraints: no external logos or invented identity marks",
        "profile-specific QA gates: 26 slides, teaching anchors remain visible, extension terms begin only in extension section, PDF smoke succeeds",
        "known missing inputs: no final institutional template supplied; this is a conference cut-down draft",
    ]
    (WORKSPACE / "profile-plan.txt").write_text("\n".join(profile) + "\n", encoding="utf-8")

    source_lines = [
        "Source notes",
        f"Date: {CREATED}",
        "",
        "Primary local source: docs/SPACEBIOBENCH_DEPTH_MASTER_DECK_V0_4_REVIEW_2026_06_03.md",
        "Base depth master: scripts/build_spacebiobench_depth_master_deck.py",
        "Conference cut-down keeps core result proof surfaces and compresses teaching/native scenes.",
        "",
    ]
    for slide in deck_spec:
        asset = slide.get("asset")
        source_lines.append(f"Slide {slide['slide']}: {base.rel(Path(asset)) if asset else 'native editable conference scene'}")
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
        "4-8 compressed teaching methods, model-view bridge, and model-family setup",
        "9-16 result spine with visible data-card strips",
        "17-19 v7/v8 boundary and hypothesis layer",
        "20-22 v9 platform and public-bulk readiness",
        "23-25 extension tracks and blocker status",
        "26 claim boundary plus roadmap close",
        "",
        "Conference rule: no teaching anchor from v0.4 disappears; each is compressed into an on-slide cue.",
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
        "pptx": PROJECT_OUTPUT_DIR / "spacebiobench_conference_deck_v0_5.pptx",
        "contact_sheet": PROJECT_OUTPUT_DIR / "spacebiobench_conference_deck_contact_sheet.png",
        "artifact_manifest": PROJECT_OUTPUT_DIR / "spacebiobench_conference_deck_artifact_manifest.json",
        "speaker_notes": PROJECT_OUTPUT_DIR / "spacebiobench_conference_deck_speaker_notes.md",
    }
    shutil.copy2(FINAL_PPTX, copied["pptx"])
    shutil.copy2(CONTACT_SHEET, copied["contact_sheet"])
    shutil.copy2(BUILD_MANIFEST, copied["artifact_manifest"])

    note_lines = ["# SpaceBio-Bench Conference Deck Speaker Notes", "", f"Date: {CREATED}", ""]
    for slide in deck_spec:
        note_lines.extend([f"## Slide {slide['slide']}: {slide['title']}", "", slide["notes"], ""])
    copied["speaker_notes"].write_text("\n".join(note_lines), encoding="utf-8")

    production_manifest = {
        "created": CREATED,
        "artifact_role": "26-slide SpaceBio-Bench v0.5 methodology-insert conference draft",
        "workspace": base.rel(WORKSPACE),
        "presentation_skill": "artifact-tool presentation build",
        "decision_source": "docs/SPACEBIOBENCH_DEPTH_MASTER_DECK_V0_4_REVIEW_2026_06_03.md",
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
            "conference_compression": "v0.4 methods slides 4-9 compressed into 4-8 with model-view and model-family bridges; v9 glossary/readiness compressed into slide 20; close matrix and roadmap combined into slide 26",
            "data_card_rule": "result and platform proof slides retain compact visible data/status cards",
            "known_status": "conference cut-down draft, not final institutional-template deck",
        },
    }
    production_path = PROJECT_OUTPUT_DIR / "spacebiobench_conference_deck_manifest.json"
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
