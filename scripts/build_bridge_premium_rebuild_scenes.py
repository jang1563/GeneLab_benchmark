#!/usr/bin/env python3
"""Render premium rebuild bridge-method slide scenes."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "output" / ".matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle


OUT_ROOT = ROOT / "output" / "premium_bridge_rebuild_scenes"
SLIDE_DPI = 240
CANVAS = {"width_px": 3840, "height_px": 2160, "aspect_ratio": "16:9"}
CREATED = "2026-06-01"

COLORS = {
    "bg": "#F4F0E8",
    "paper": "#FCFCFA",
    "ink": "#17212B",
    "muted": "#5D6978",
    "rule": "#AEB8C5",
    "blue": "#286EA6",
    "green": "#15815E",
    "amber": "#A36F13",
    "red": "#B91C1C",
    "teal": "#159A8A",
    "purple": "#6D3EDB",
    "shadow": "#9C9487",
    "slate": "#283544",
}


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def write_json(path: Path, data: Any) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def y_from_slide(value: float) -> float:
    return 1.0 - value


def count_words(items: list[dict[str, Any]]) -> int:
    words = 0
    for item in items:
        content = str(item.get("content", ""))
        words += len([token for token in content.replace("/", " ").split() if token])
    return words


def make_figure() -> tuple[plt.Figure, plt.Axes]:
    fig = plt.figure(figsize=(16, 9), dpi=SLIDE_DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    return fig, ax


def draw_canvas(ax: plt.Axes, seed: int) -> None:
    h, w = 720, 1280
    yy, xx = np.mgrid[0:h, 0:w]
    rng = np.random.default_rng(seed)
    base = np.zeros((h, w, 3), dtype=float)
    base[..., 0] = 0.955
    base[..., 1] = 0.938
    base[..., 2] = 0.904
    grain = rng.normal(0, 0.0045, size=(h, w, 1))
    vignette = ((xx / w - 0.55) ** 2 + (yy / h - 0.45) ** 2)[..., None] * np.array([0.045, 0.035, 0.020])
    cool_left = (1 - xx / w)[..., None] * np.array([0.000, 0.006, 0.010])
    texture = np.clip(base + grain - vignette + cool_left, 0, 1)
    ax.imshow(texture, extent=(0, 1, 0, 1), origin="lower", zorder=0)
    ax.set_aspect("auto")

    for y in [0.190, 0.328, 0.500, 0.672, 0.835]:
        ax.plot([0.055, 0.945], [y, y], color=COLORS["rule"], alpha=0.20, linewidth=0.85, transform=ax.transAxes, zorder=1)
    for x in [0.075, 0.275, 0.500, 0.725, 0.925]:
        ax.plot([x, x], [0.095, 0.870], color=COLORS["rule"], alpha=0.10, linewidth=0.75, transform=ax.transAxes, zorder=1)


def rect(
    ax: plt.Axes,
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    face: str = "paper",
    edge: str | None = "rule",
    alpha: float = 1.0,
    lw: float = 0.8,
    z: float = 3,
) -> None:
    ax.add_patch(
        Rectangle(
            (x, y),
            w,
            h,
            transform=ax.transAxes,
            facecolor=COLORS.get(face, face),
            edgecolor=COLORS.get(edge, edge) if edge else "none",
            linewidth=lw,
            alpha=alpha,
            zorder=z,
        )
    )


def shadow_rect(ax: plt.Axes, x: float, y: float, w: float, h: float, *, alpha: float = 0.16, z: float = 2) -> None:
    rect(ax, x + 0.010, y - 0.012, w, h, face="shadow", edge=None, alpha=alpha, z=z)


def arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float], color: str = "muted", *, alpha: float = 0.65, lw: float = 1.2, z: float = 6) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=12,
            linewidth=lw,
            color=COLORS[color],
            alpha=alpha,
            transform=ax.transAxes,
            zorder=z,
        )
    )


def node(ax: plt.Axes, x: float, y: float, label: str, color: str, *, radius: float = 0.024, z: float = 7) -> None:
    ax.add_patch(Circle((x + 0.004, y - 0.006), radius, transform=ax.transAxes, facecolor=COLORS["shadow"], edgecolor="none", alpha=0.16, zorder=z - 1))
    ax.add_patch(Circle((x, y), radius, transform=ax.transAxes, facecolor=COLORS["paper"], edgecolor=COLORS[color], linewidth=1.45, alpha=0.98, zorder=z))
    ax.text(x, y, label, color=COLORS[color], fontsize=7.0, fontweight="bold", ha="center", va="center", transform=ax.transAxes, zorder=z + 1)


def draw_b2(ax: plt.Axes) -> None:
    draw_canvas(ax, 20260621)
    rect(ax, 0.055, 0.320, 0.890, 0.365, face="#FBFAF6", edge=None, alpha=0.50, z=2)
    rect(ax, 0.065, 0.320, 0.190, 0.365, face="#F4F6F5", edge=None, alpha=0.44, z=2.2)
    rect(ax, 0.745, 0.320, 0.190, 0.365, face="#F7F4FD", edge=None, alpha=0.36, z=2.2)
    ax.plot([0.105, 0.865], [0.500, 0.500], color=COLORS["rule"], linewidth=1.15, alpha=0.48, transform=ax.transAxes, zorder=4)

    for idx, (dx, dy, a) in enumerate([(0.000, 0.000, 0.95), (0.018, 0.018, 0.55), (0.036, 0.036, 0.28)]):
        shadow_rect(ax, 0.105 + dx, 0.405 + dy, 0.150, 0.205, alpha=0.08 * a, z=4 + idx * 0.05)
        rect(ax, 0.100 + dx, 0.412 + dy, 0.150, 0.205, face="paper", edge="rule", alpha=a, lw=0.70, z=5 + idx * 0.05)
    for y in [0.572, 0.542, 0.512, 0.482, 0.452]:
        ax.plot([0.126, 0.222], [y, y], color=COLORS["rule"], alpha=0.64, linewidth=0.82, transform=ax.transAxes, zorder=7)
    for x, y in [(0.218, 0.435), (0.198, 0.435), (0.178, 0.435)]:
        ax.add_patch(Circle((x, y), 0.006, transform=ax.transAxes, facecolor=COLORS["muted"], edgecolor="none", alpha=0.38, zorder=7))

    ax.add_patch(Circle((0.330, 0.505), 0.092, transform=ax.transAxes, facecolor="#FDFDFB", edgecolor=COLORS["green"], linewidth=1.35, alpha=0.96, zorder=5))
    ax.add_patch(Circle((0.330, 0.505), 0.054, transform=ax.transAxes, facecolor="none", edgecolor=COLORS["green"], linewidth=1.05, alpha=0.34, zorder=6))
    theta = np.linspace(0.15, 2.85 * np.pi, 160)
    ax.plot(0.330 + np.cos(theta) * 0.088, 0.505 + np.sin(theta) * 0.040, color=COLORS["green"], alpha=0.20, linewidth=1.0, transform=ax.transAxes, zorder=5)

    sample_x = [0.455, 0.485, 0.515, 0.455, 0.485, 0.515, 0.455, 0.485, 0.515]
    sample_y = [0.565, 0.565, 0.565, 0.505, 0.505, 0.505, 0.445, 0.445, 0.445]
    for idx, (x, y) in enumerate(zip(sample_x, sample_y)):
        color = "blue" if idx % 2 == 0 else "teal"
        ax.add_patch(Circle((x, y), 0.013, transform=ax.transAxes, facecolor=COLORS["paper"], edgecolor=COLORS[color], linewidth=1.05, alpha=0.98, zorder=7))
    for x in [0.455, 0.485, 0.515]:
        ax.plot([x, x], [0.430, 0.582], color=COLORS["rule"], alpha=0.14, linewidth=0.8, transform=ax.transAxes, zorder=4)

    for y, color in [(0.535, "amber"), (0.465, "slate")]:
        shadow_rect(ax, 0.590, y - 0.028, 0.108, 0.056, alpha=0.07, z=5)
        rect(ax, 0.584, y - 0.025, 0.108, 0.056, face="paper", edge=color, alpha=0.96, lw=0.95, z=6)
        ax.plot([0.606, 0.670], [y, y], color=COLORS[color], alpha=0.50, linewidth=0.9, transform=ax.transAxes, zorder=7)

    rect(ax, 0.715, 0.430, 0.105, 0.150, face="paper", edge="teal", alpha=0.95, lw=0.95, z=6)
    for y in [0.552, 0.530, 0.508, 0.486, 0.464]:
        ax.plot([0.736, 0.799], [y, y], color=COLORS["teal"], alpha=0.28, linewidth=1.0, transform=ax.transAxes, zorder=7)

    shadow_rect(ax, 0.846, 0.370, 0.088, 0.265, alpha=0.17, z=5)
    rect(ax, 0.835, 0.386, 0.088, 0.265, face="paper", edge="rule", alpha=0.98, lw=0.85, z=7)
    ax.plot([0.850, 0.908], [0.603, 0.603], color=COLORS["purple"], alpha=0.55, linewidth=1.1, transform=ax.transAxes, zorder=8)
    for y in [0.568, 0.535, 0.502, 0.469, 0.436]:
        ax.plot([0.850, 0.904], [y, y], color=COLORS["rule"], alpha=0.56, linewidth=0.75, transform=ax.transAxes, zorder=8)
    for x, color in [(0.852, "green"), (0.870, "blue"), (0.888, "amber"), (0.906, "teal")]:
        ax.add_patch(Circle((x, 0.410), 0.0065, transform=ax.transAxes, facecolor=COLORS[color], edgecolor="none", alpha=0.80, zorder=8))

    for start, end in [((0.255, 0.500), (0.300, 0.500)), ((0.395, 0.500), (0.435, 0.500)), ((0.532, 0.500), (0.575, 0.500)), ((0.695, 0.500), (0.708, 0.500)), ((0.820, 0.500), (0.832, 0.500))]:
        arrow(ax, start, end, "muted", alpha=0.42, lw=1.1, z=8)


def draw_b3(ax: plt.Axes) -> None:
    draw_canvas(ax, 20260622)
    rect(ax, 0.055, 0.320, 0.890, 0.365, face="#FBFAF6", edge=None, alpha=0.50, z=2)
    rect(ax, 0.080, 0.335, 0.565, 0.335, face="#F2F8F5", edge=None, alpha=0.50, z=2.2)
    rect(ax, 0.690, 0.335, 0.205, 0.335, face="#FFF5F5", edge=None, alpha=0.54, z=2.2)

    theta = np.linspace(0.0, np.pi, 120)
    ax.plot(0.355 + np.cos(theta) * 0.255, 0.505 + np.sin(theta) * 0.120, color=COLORS["green"], alpha=0.18, linewidth=1.6, transform=ax.transAxes, zorder=4)
    ax.plot([0.150, 0.610], [0.505, 0.505], color=COLORS["green"], alpha=0.30, linewidth=5.0, transform=ax.transAxes, zorder=4)
    for x, label in [(0.165, "M1"), (0.300, "M2"), (0.435, "M3")]:
        node(ax, x, 0.505, label, "green", radius=0.036, z=7)
        rect(ax, x - 0.045, 0.408, 0.090, 0.050, face="paper", edge="rule", alpha=0.72, lw=0.55, z=5)
        for y in [0.443, 0.430, 0.417]:
            ax.plot([x - 0.030, x + 0.030], [y, y], color=COLORS["rule"], alpha=0.40, linewidth=0.65, transform=ax.transAxes, zorder=6)

    rect(ax, 0.545, 0.430, 0.090, 0.150, face="paper", edge="rule", alpha=0.94, lw=0.80, z=6)
    for y in [0.555, 0.530, 0.505, 0.480, 0.455]:
        ax.plot([0.565, 0.616], [y, y], color=COLORS["rule"], alpha=0.62, linewidth=0.78, transform=ax.transAxes, zorder=7)

    ax.plot([0.670, 0.670], [0.318, 0.715], color=COLORS["red"], alpha=0.78, linewidth=1.55, transform=ax.transAxes, zorder=9)
    ax.plot([0.682, 0.682], [0.340, 0.694], color=COLORS["red"], alpha=0.22, linewidth=0.90, transform=ax.transAxes, zorder=8)
    node(ax, 0.775, 0.505, "M4", "red", radius=0.040, z=8)
    rect(ax, 0.720, 0.395, 0.110, 0.070, face="paper", edge="red", alpha=0.70, lw=0.70, z=6)
    for y in [0.440, 0.425, 0.410]:
        ax.plot([0.740, 0.812], [y, y], color=COLORS["red"], alpha=0.22, linewidth=0.70, transform=ax.transAxes, zorder=7)
    ax.add_patch(Circle((0.850, 0.630), 0.040, transform=ax.transAxes, facecolor=COLORS["paper"], edgecolor=COLORS["red"], linewidth=1.15, alpha=0.96, zorder=8))
    ax.plot([0.832, 0.868], [0.630, 0.630], color=COLORS["red"], alpha=0.72, linewidth=1.0, transform=ax.transAxes, zorder=9)
    ax.plot([0.850, 0.850], [0.612, 0.648], color=COLORS["red"], alpha=0.72, linewidth=1.0, transform=ax.transAxes, zorder=9)
    arrow(ax, (0.815, 0.535), (0.842, 0.590), "red", alpha=0.50, lw=1.0, z=9)


def draw_b4(ax: plt.Axes) -> None:
    draw_canvas(ax, 20260623)
    rect(ax, 0.055, 0.320, 0.890, 0.365, face="#FBFAF6", edge=None, alpha=0.50, z=2)
    rect(ax, 0.105, 0.430, 0.610, 0.165, face="#F3F8F8", edge=None, alpha=0.60, z=3.2)
    rect(ax, 0.105, 0.330, 0.610, 0.070, face="#FFF7F7", edge=None, alpha=0.54, z=3.2)
    ax.plot([0.125, 0.705], [0.515, 0.515], color=COLORS["rule"], alpha=0.66, linewidth=1.35, transform=ax.transAxes, zorder=4)
    ax.plot([0.125, 0.755], [0.365, 0.365], color=COLORS["red"], alpha=0.34, linewidth=1.15, transform=ax.transAxes, zorder=4)

    for x, label in [(0.145, "M1"), (0.202, "M2"), (0.259, "M3")]:
        node(ax, x, 0.515, label, "green", radius=0.023, z=7)
    node(ax, 0.145, 0.365, "M4", "red", radius=0.025, z=7)

    for x, color, h in [(0.360, "blue", 0.170), (0.515, "teal", 0.170), (0.665, "purple", 0.170)]:
        shadow_rect(ax, x - 0.038, 0.430, 0.076, h, alpha=0.08, z=5)
        rect(ax, x - 0.034, 0.442, 0.068, h, face="paper", edge=color, alpha=0.94, lw=0.95, z=6)
        ax.plot([x, x], [0.457, 0.582], color=COLORS[color], alpha=0.65, linewidth=1.0, transform=ax.transAxes, zorder=7)
        for y in [0.559, 0.537, 0.493]:
            ax.plot([x - 0.018, x + 0.018], [y, y], color=COLORS[color], alpha=0.28, linewidth=0.80, transform=ax.transAxes, zorder=7)

    ax.plot([0.742, 0.742], [0.335, 0.695], color=COLORS["red"], alpha=0.78, linewidth=1.55, transform=ax.transAxes, zorder=9)
    ax.plot([0.758, 0.758], [0.355, 0.675], color=COLORS["red"], alpha=0.20, linewidth=0.90, transform=ax.transAxes, zorder=8)
    arrow(ax, (0.700, 0.515), (0.810, 0.580), "muted", alpha=0.55, lw=1.1, z=8)
    arrow(ax, (0.178, 0.365), (0.815, 0.365), "red", alpha=0.34, lw=1.0, z=8)
    ax.add_patch(Circle((0.842, 0.590), 0.040, transform=ax.transAxes, facecolor=COLORS["paper"], edgecolor=COLORS["red"], linewidth=1.10, alpha=0.96, zorder=9))
    for y in [0.590, 0.610, 0.570]:
        ax.plot([0.824, 0.860], [y, y], color=COLORS["red"], alpha=0.52, linewidth=0.80, transform=ax.transAxes, zorder=10)


DRAWERS = {
    "b2_study_to_task_premium": draw_b2,
    "b3_mission_held_out_premium": draw_b3,
    "b4_train_only_guard_premium": draw_b4,
}


SLIDES: list[dict[str, Any]] = [
    {
        "slide_id": "b2_study_to_task_premium",
        "decision_headline": "A task is a traceable contract, not a loose matrix",
        "audience_question": "What is the unit of the benchmark?",
        "claim_boundary": "a task contract ties public source records to mission context, samples, labels, tissue, assay, and evaluation scope",
        "visual_move": "source evidence consolidates into one task contract surface",
        "content_brief": "docs/VISUAL_BRIDGE_CONTENT_BRIEFS_B1_B4_2026_06_01.md",
        "technical_preflight": "docs/VISUAL_BRIDGE_PREMIUM_REBUILD_CRITIQUE_2026_06_01.md",
        "evidence_sources": [
            {"path": "docs/VISUAL_BRIDGE_PREMIUM_REBUILD_CRITIQUE_2026_06_01.md", "role": "premium rebuild criteria"},
            {"path": "v9/source_inventory.csv", "role": "public source inventory"},
            {"path": "v9/task_manifest_index.csv", "role": "task contract index"},
            {"path": "docs/VISUAL_BRIDGE_CONTENT_BRIEFS_B1_B4_2026_06_01.md", "role": "B2 brief"},
        ],
        "forbidden_visible_terms": ["raw accession", "payload", "artifact", "RRRM", "alpha", "LOMO", "/Users/", "function", "class"],
        "overlay": {
            "text": [
                {"id": "headline", "role": "decision_headline", "content": "A task is a traceable contract, not a loose matrix", "x": 0.065, "y": 0.128, "font_pt": 26, "color": "ink", "z": "Z3"},
                {"id": "source", "role": "primary_callout", "content": "public source", "x": 0.108, "y": 0.338, "font_pt": 9.5, "color": "muted", "z": "Z3"},
                {"id": "mission", "role": "primary_callout", "content": "mission context", "x": 0.282, "y": 0.338, "font_pt": 9.5, "color": "green", "z": "Z3"},
                {"id": "samples", "role": "primary_callout", "content": "sample set", "x": 0.448, "y": 0.338, "font_pt": 9.5, "color": "blue", "z": "Z3"},
                {"id": "labels", "role": "primary_callout", "content": "label contrast", "x": 0.590, "y": 0.338, "font_pt": 9.5, "color": "amber", "z": "Z3"},
                {"id": "tissue_assay", "role": "primary_callout", "content": "tissue + assay", "x": 0.715, "y": 0.338, "font_pt": 9.5, "color": "teal", "z": "Z3"},
                {"id": "contract", "role": "primary_callout", "content": "task contract", "x": 0.832, "y": 0.338, "font_pt": 9.5, "color": "purple", "z": "Z3"},
            ],
            "status_labels": [
                {"id": "scope", "role": "claim_boundary", "content": "source rows remain auditable", "x": 0.065, "y": 0.865, "font_pt": 7.8, "color": "muted", "z": "Z4"},
                {"id": "caveat", "role": "trust_caveat", "content": "Exact rows live outside this figure.", "x": 0.065, "y": 0.900, "font_pt": 7.4, "color": "muted", "z": "Z4"},
            ],
            "focus_marks": [{"id": "contract_path", "role": "flow_path", "shape": "evidence_rail", "x0": 0.105, "x1": 0.865, "y": 0.500, "color": "muted", "z": "Z5"}],
        },
    },
    {
        "slide_id": "b3_mission_held_out_premium",
        "decision_headline": "The test set is an entire mission",
        "audience_question": "What is held out during evaluation?",
        "claim_boundary": "evaluation hides a mission-level group rather than a random sample",
        "visual_move": "training missions feed the model while one mission remains behind the boundary",
        "content_brief": "docs/VISUAL_BRIDGE_CONTENT_BRIEFS_B1_B4_2026_06_01.md",
        "technical_preflight": "docs/VISUAL_BRIDGE_PREMIUM_REBUILD_CRITIQUE_2026_06_01.md",
        "evidence_sources": [
            {"path": "docs/VISUAL_BRIDGE_PREMIUM_REBUILD_CRITIQUE_2026_06_01.md", "role": "premium rebuild criteria"},
            {"path": "docs/VISUAL_METHODS_EXPLANATION_GAP_MAP_2026_06_01.md", "role": "held-out mission explanation gap"},
            {"path": "docs/VISUAL_BRIDGE_CONTENT_BRIEFS_B1_B4_2026_06_01.md", "role": "B3 brief"},
        ],
        "forbidden_visible_terms": ["LOMO", "random CV", "cross-validation", "payload", "RRRM", "alpha", "NES", "macro-F1", "/Users/"],
        "overlay": {
            "text": [
                {"id": "headline", "role": "decision_headline", "content": "The test set is an entire mission", "x": 0.065, "y": 0.128, "font_pt": 27, "color": "ink", "z": "Z3"},
                {"id": "train", "role": "primary_callout", "content": "training missions", "x": 0.145, "y": 0.338, "font_pt": 10.2, "color": "green", "z": "Z3"},
                {"id": "model", "role": "primary_callout", "content": "fit before the boundary", "x": 0.535, "y": 0.338, "font_pt": 9.5, "color": "muted", "z": "Z3"},
                {"id": "hidden", "role": "primary_callout", "content": "held-out mission", "x": 0.725, "y": 0.338, "font_pt": 10.2, "color": "red", "z": "Z3"},
                {"id": "score", "role": "secondary_label", "content": "score only after training", "x": 0.805, "y": 0.285, "font_pt": 8.4, "color": "red", "z": "Z3"},
            ],
            "status_labels": [
                {"id": "scope", "role": "claim_boundary", "content": "mission-held-out evaluation", "x": 0.065, "y": 0.865, "font_pt": 7.8, "color": "muted", "z": "Z4"},
                {"id": "caveat", "role": "trust_caveat", "content": "Hidden-mission samples stay outside training.", "x": 0.065, "y": 0.900, "font_pt": 7.4, "color": "muted", "z": "Z4"},
            ],
            "focus_marks": [{"id": "boundary", "role": "boundary", "shape": "vertical_rule", "x": 0.670, "y0": 0.285, "y1": 0.720, "color": "red", "z": "Z5"}],
        },
    },
    {
        "slide_id": "b4_train_only_guard_premium",
        "decision_headline": "Training choices stop before the hidden mission",
        "audience_question": "How does the benchmark avoid leakage?",
        "claim_boundary": "feature selection, scaling, and model fitting are learned from training missions only",
        "visual_move": "hidden mission bypasses the processing workspace and joins only at scoring",
        "content_brief": "docs/VISUAL_BRIDGE_CONTENT_BRIEFS_B1_B4_2026_06_01.md",
        "technical_preflight": "docs/VISUAL_BRIDGE_PREMIUM_REBUILD_CRITIQUE_2026_06_01.md",
        "evidence_sources": [
            {"path": "docs/VISUAL_BRIDGE_PREMIUM_REBUILD_CRITIQUE_2026_06_01.md", "role": "premium rebuild criteria"},
            {"path": "docs/VISUAL_METHODS_EXPLANATION_GAP_MAP_2026_06_01.md", "role": "leakage prevention explanation gap"},
            {"path": "docs/VISUAL_BRIDGE_CONTENT_BRIEFS_B1_B4_2026_06_01.md", "role": "B4 brief"},
        ],
        "forbidden_visible_terms": ["LOMO", "payload", "artifact", "macro-F1", "/Users/", "sklearn", "StandardScaler", "fit_transform"],
        "overlay": {
            "text": [
                {"id": "headline", "role": "decision_headline", "content": "Training choices stop before the hidden mission", "x": 0.065, "y": 0.128, "font_pt": 27, "color": "ink", "z": "Z3"},
                {"id": "train", "role": "primary_callout", "content": "training missions only", "x": 0.125, "y": 0.345, "font_pt": 9.5, "color": "green", "z": "Z3"},
                {"id": "features", "role": "process_verb", "content": "choose features", "x": 0.325, "y": 0.365, "font_pt": 9.0, "color": "blue", "z": "Z3"},
                {"id": "scale", "role": "process_verb", "content": "scale", "x": 0.498, "y": 0.365, "font_pt": 9.0, "color": "teal", "z": "Z3"},
                {"id": "fit", "role": "process_verb", "content": "fit model", "x": 0.640, "y": 0.365, "font_pt": 9.0, "color": "purple", "z": "Z3"},
                {"id": "hidden", "role": "primary_callout", "content": "hidden mission bypass", "x": 0.205, "y": 0.655, "font_pt": 8.8, "color": "red", "z": "Z3"},
                {"id": "score", "role": "process_verb", "content": "score", "x": 0.823, "y": 0.540, "font_pt": 9.0, "color": "red", "z": "Z3"},
            ],
            "status_labels": [
                {"id": "scope", "role": "claim_boundary", "content": "train-only processing", "x": 0.065, "y": 0.865, "font_pt": 7.8, "color": "muted", "z": "Z4"},
                {"id": "caveat", "role": "trust_caveat", "content": "No hidden mission data shapes the feature space.", "x": 0.065, "y": 0.900, "font_pt": 7.4, "color": "muted", "z": "Z4"},
            ],
            "focus_marks": [{"id": "guard", "role": "guard", "shape": "vertical_rule", "x": 0.742, "y0": 0.330, "y1": 0.700, "color": "red", "z": "Z5"}],
        },
    },
]


def build_contract(slide: dict[str, Any]) -> dict[str, Any]:
    slide_id = slide["slide_id"]
    slide_dir = OUT_ROOT / slide_id
    slide_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "scene_plate": rel(slide_dir / "scene_plate.png"),
        "rendered_preview_png": rel(slide_dir / "rendered_preview.png"),
        "rendered_preview_pdf": rel(slide_dir / "rendered_preview.pdf"),
        "overlay_spec": rel(slide_dir / "overlay_spec.json"),
        "manifest": rel(slide_dir / "manifest.json"),
        "qa": rel(slide_dir / "qa.json"),
    }
    text_items = slide["overlay"]["text"] + slide["overlay"]["status_labels"]
    overlay = {
        "slide_id": slide_id,
        "stage": "pre_render",
        "canvas": CANVAS,
        "coordinate_system": "normalized_0_1",
        "text": slide["overlay"]["text"],
        "status_labels": slide["overlay"]["status_labels"],
        "focus_marks": slide["overlay"]["focus_marks"],
        "visible_word_count": count_words(text_items),
        "visible_word_budget": 45,
        "forbidden_visible_terms": slide["forbidden_visible_terms"],
    }
    manifest = {
        "slide_id": slide_id,
        "created": CREATED,
        "stage": "pre_render",
        "content_brief": slide["content_brief"],
        "technical_preflight": slide["technical_preflight"],
        "audience_question": slide["audience_question"],
        "decision_headline": slide["decision_headline"],
        "visual_move": slide["visual_move"],
        "claim_boundary": slide["claim_boundary"],
        "evidence_sources": slide["evidence_sources"],
        "generator": "scripts/build_bridge_premium_rebuild_scenes.py",
        "outputs": outputs,
        "qa": {"stage": "pre_render", "pre_render_gate_expected": True, "post_render_required_before_use": True},
    }
    qa = {
        "slide_id": slide_id,
        "stage": "pre_render",
        "created": CREATED,
        "pre_render_gate": {
            "content_brief_declared": True,
            "technical_preflight_declared": True,
            "evidence_sources_declared": True,
            "claim_boundary_declared": True,
            "visible_text_word_count": overlay["visible_word_count"],
            "visible_text_budget": overlay["visible_word_budget"],
            "output_paths_declared": True,
            "overlay_spec_declared": True,
            "manifest_declared": True,
        },
        "manual_review_pending": [
            "full_size_render_inspection",
            "thumbnail_contact_sheet_inspection",
            "text_overlap_check",
            "premium_quality_gate",
        ],
    }
    write_json(slide_dir / "overlay_spec.json", overlay)
    write_json(slide_dir / "manifest.json", manifest)
    write_json(slide_dir / "qa.json", qa)
    return {"manifest": manifest, "overlay": overlay, "qa": qa}


def render_overlay(ax: plt.Axes, overlay: dict[str, Any]) -> None:
    for item in list(overlay["text"]) + list(overlay["status_labels"]):
        color = COLORS.get(str(item.get("color", "ink")), COLORS["ink"])
        role = str(item.get("role", ""))
        weight = "bold" if role in {"decision_headline", "primary_callout", "process_verb", "claim_boundary"} else "normal"
        va = "top" if role == "decision_headline" else "center"
        ax.text(
            float(item["x"]),
            y_from_slide(float(item["y"])),
            str(item["content"]),
            color=color,
            fontsize=float(item.get("font_pt", 8.0)),
            fontweight=weight,
            ha="left",
            va=va,
            transform=ax.transAxes,
            zorder=20,
        )


def render_slide(slide: dict[str, Any]) -> dict[str, str]:
    contract = build_contract(slide)
    slide_id = slide["slide_id"]
    slide_dir = OUT_ROOT / slide_id
    manifest = contract["manifest"]
    overlay = contract["overlay"]
    outputs = manifest["outputs"]
    scene_plate = ROOT / outputs["scene_plate"]
    preview_png = ROOT / outputs["rendered_preview_png"]
    preview_pdf = ROOT / outputs["rendered_preview_pdf"]

    fig, ax = make_figure()
    DRAWERS[slide_id](ax)
    fig.savefig(scene_plate, dpi=SLIDE_DPI, facecolor=COLORS["bg"])
    plt.close(fig)

    fig, ax = make_figure()
    DRAWERS[slide_id](ax)
    render_overlay(ax, overlay)
    fig.savefig(preview_png, dpi=SLIDE_DPI, facecolor=COLORS["bg"])
    fig.savefig(preview_pdf, facecolor=COLORS["bg"])
    plt.close(fig)

    manifest["stage"] = "post_render"
    manifest["renderer"] = "scripts/build_bridge_premium_rebuild_scenes.py"
    manifest["qa"]["stage"] = "post_render"
    manifest["qa"]["render_outputs_exist"] = True
    overlay["stage"] = "post_render"
    qa = contract["qa"]
    qa["stage"] = "post_render"
    qa["post_render_gate"] = {
        "rendered_outputs": {
            "scene_plate": rel(scene_plate),
            "rendered_preview_png": rel(preview_png),
            "rendered_preview_pdf": rel(preview_pdf),
        },
        "manual_full_size_inspection": "pending",
        "manual_thumbnail_inspection": "pending",
        "text_overlap_check": "pending",
        "premium_quality_gate": "pending",
        "visual_verdict": "premium rebuild render; awaiting manual inspection",
    }
    write_json(slide_dir / "manifest.json", manifest)
    write_json(slide_dir / "overlay_spec.json", overlay)
    write_json(slide_dir / "qa.json", qa)
    return {"scene_plate": rel(scene_plate), "rendered_preview_png": rel(preview_png), "rendered_preview_pdf": rel(preview_pdf)}


def build_contact_sheet(slide_ids: list[str]) -> dict[str, str]:
    output = OUT_ROOT / "bridge_methods_premium_rebuild_contact_sheet.png"
    fig = plt.figure(figsize=(19.0, 5.7), dpi=220)
    grid = fig.add_gridspec(1, len(slide_ids), left=0.018, right=0.988, top=0.790, bottom=0.060, wspace=0.026)
    fig.suptitle("Premium bridge-method rebuild candidates: B2-B4", x=0.018, y=0.946, ha="left", fontsize=12.8, fontweight="bold")
    fig.text(
        0.018,
        0.882,
        "Rebuild criterion: visible thesis at thumbnail scale, larger evidence surface, stronger depth, sparse interpretation overlay.",
        fontsize=7.8,
        color=COLORS["muted"],
        ha="left",
    )
    slides_by_id = {slide["slide_id"]: slide for slide in SLIDES}
    for idx, slide_id in enumerate(slide_ids):
        ax = fig.add_subplot(grid[0, idx])
        source = OUT_ROOT / slide_id / "rendered_preview.png"
        ax.imshow(mpimg.imread(source))
        claim = slides_by_id[slide_id]["decision_headline"]
        ax.set_title(f"{idx + 1}. {slide_id} | {claim}", loc="left", fontsize=6.7, pad=3)
        ax.axis("off")
    fig.savefig(output, dpi=220, facecolor="white")
    plt.close(fig)
    manifest = {
        "contact_sheet": rel(output),
        "slides": [
            {
                "slide_id": slide_id,
                "claim": slides_by_id[slide_id]["decision_headline"],
                "source": rel(OUT_ROOT / slide_id / "rendered_preview.png"),
            }
            for slide_id in slide_ids
        ],
    }
    manifest_path = OUT_ROOT / "bridge_methods_premium_rebuild_contact_sheet_manifest.json"
    write_json(manifest_path, manifest)
    return {"contact_sheet": rel(output), "manifest": rel(manifest_path)}


def selected_slides(slide_id: str) -> list[dict[str, Any]]:
    if slide_id == "all":
        return SLIDES
    for slide in SLIDES:
        if slide["slide_id"] == slide_id:
            return [slide]
    known = ", ".join(slide["slide_id"] for slide in SLIDES)
    raise SystemExit(f"Unknown slide_id {slide_id!r}. Known slide IDs: {known}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--slide", default="all", help="Slide ID to render, or 'all'.")
    args = parser.parse_args()
    slides = selected_slides(args.slide)
    rendered = {slide["slide_id"]: render_slide(slide) for slide in slides}
    contact_sheet = build_contact_sheet([slide["slide_id"] for slide in slides]) if len(slides) > 1 else None
    print(json.dumps({"rendered": rendered, "contact_sheet": contact_sheet}, indent=2))


if __name__ == "__main__":
    main()
