#!/usr/bin/env python3
"""Build layered premium slide-scene pilots from audited proof figures."""

from __future__ import annotations

import json
import os
import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "output" / ".matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Arc, Ellipse, FancyArrowPatch, Rectangle


PREMIUM_FIGURE_DIR = ROOT / "output" / "premium_figures"
OUT_DIR = ROOT / "output" / "premium_slide_scenes"

FIG1_PROOF = PREMIUM_FIGURE_DIR / "premium_fig1_tissue_transfer_hierarchy.png"
FIG1_MANIFEST = PREMIUM_FIGURE_DIR / "manifests" / "premium_fig1_tissue_transfer_hierarchy_manifest.json"
FIG2_PROOF = PREMIUM_FIGURE_DIR / "manuscript_variants" / "premium_fig2_pathway_rescue_manuscript.png"
FIG2_MANIFEST = PREMIUM_FIGURE_DIR / "manifests" / "premium_fig2_pathway_rescue_manuscript_manifest.json"
FIG3_PROOF = PREMIUM_FIGURE_DIR / "manuscript_variants" / "premium_fig3_model_tier_comparison_manuscript.png"
FIG3_MANIFEST = PREMIUM_FIGURE_DIR / "manifests" / "premium_fig3_model_tier_comparison_manuscript_manifest.json"
FIG6_PROOF = PREMIUM_FIGURE_DIR / "manuscript_variants" / "premium_fig6_organoid_biology_check_manuscript.png"
FIG6_MANIFEST = PREMIUM_FIGURE_DIR / "manifests" / "premium_fig6_organoid_biology_check_manuscript_manifest.json"
NUMERIC_AUDIT = ROOT / "output" / "premium_audit" / "premium_visual_numeric_audit.csv"

SLIDE_DPI = 240

COLORS = {
    "bg": "#071019",
    "bg2": "#0A1722",
    "ink": "#F8FAFC",
    "muted": "#AAB5C4",
    "grid": "#8FA3B5",
    "blue": "#2F6C9F",
    "amber": "#D99A22",
    "green": "#1F8F5F",
    "teal": "#1F9D8A",
    "purple": "#7C3AED",
    "violet": "#9B7AF7",
    "shadow": "#02060A",
    "paper": "#F9FAFB",
}


def load_rgb(path: Path) -> np.ndarray:
    image = mpimg.imread(path).astype(float)
    if image.max() > 1.0:
        image = image / 255.0
    if image.shape[-1] == 4:
        alpha = image[..., 3:4]
        image = image[..., :3] * alpha + (1.0 - alpha)
    return image[..., :3]


def crop_fig1_plot(image: np.ndarray) -> np.ndarray:
    h, w = image.shape[:2]
    y0 = int(h * 0.235)
    y1 = int(h * 0.930)
    x0 = int(w * 0.020)
    x1 = int(w * 0.975)
    return image[y0:y1, x0:x1, :]


def crop_fraction(image: np.ndarray, x0: float, x1: float, y0: float, y1: float) -> np.ndarray:
    h, w = image.shape[:2]
    return image[int(h * y0) : int(h * y1), int(w * x0) : int(w * x1), :]


def draw_canvas(ax: plt.Axes, rng: np.random.Generator) -> None:
    ax.add_patch(Rectangle((0, 0), 1, 1, color=COLORS["bg"], transform=ax.transAxes, zorder=0))

    h, w = 540, 960
    yy, xx = np.mgrid[0:h, 0:w]
    base = np.zeros((h, w, 3), dtype=float)
    base[..., 0] = 0.026
    base[..., 1] = 0.055
    base[..., 2] = 0.078
    fine_noise = rng.normal(0, 0.012, size=(h, w, 1))
    left_depth = np.clip(1.25 - xx / w, 0.0, 1.0)[..., None] * np.array([0.020, 0.045, 0.070])
    vignette = 0.90 + 0.10 * np.cos((xx / w - 0.15) * np.pi)
    texture = np.clip((base + fine_noise + left_depth) * vignette[..., None], 0, 1)
    ax.imshow(texture, extent=(0, 1, 0, 1), origin="lower", zorder=0)
    ax.set_aspect("auto")

    for x in np.linspace(0.08, 0.94, 9):
        ax.plot([x, x], [0.08, 0.92], color=COLORS["grid"], alpha=0.045, linewidth=0.7, zorder=1)
    for y in np.linspace(0.12, 0.88, 7):
        ax.plot([0.05, 0.96], [y, y], color=COLORS["grid"], alpha=0.035, linewidth=0.7, zorder=1)

    for _ in range(180):
        x = rng.uniform(0.02, 0.98)
        y = rng.uniform(0.04, 0.96)
        length = rng.uniform(0.006, 0.020)
        alpha = rng.uniform(0.035, 0.090)
        ax.plot(
            [x, min(x + length, 0.99)],
            [y, y + rng.uniform(-0.002, 0.002)],
            color=COLORS["grid"],
            alpha=alpha,
            linewidth=0.45,
            zorder=1,
        )


def draw_measurement_layer(ax: plt.Axes) -> None:
    arc = Arc(
        (0.33, 0.70),
        width=0.82,
        height=0.70,
        angle=0,
        theta1=198,
        theta2=352,
        color=COLORS["blue"],
        linewidth=1.1,
        alpha=0.28,
        zorder=2,
    )
    ax.add_patch(arc)
    for x, label in [(0.10, "train"), (0.20, "train"), (0.30, "train"), (0.40, "held-out")]:
        ax.plot([x, x], [0.167, 0.197], color=COLORS["grid"], alpha=0.45, linewidth=0.9, zorder=2)
        ax.text(x, 0.135, label, color=COLORS["muted"], fontsize=6.2, ha="center", alpha=0.68, zorder=2)
    ax.plot([0.08, 0.42], [0.182, 0.182], color=COLORS["grid"], alpha=0.22, linewidth=1.0, zorder=2)
    ax.text(0.08, 0.215, "one mission held out", color=COLORS["muted"], fontsize=7.0, alpha=0.74, zorder=2)


def draw_pathway_measurement_layer(ax: plt.Axes) -> None:
    ax.plot([0.070, 0.365], [0.185, 0.185], color=COLORS["grid"], alpha=0.24, linewidth=1.0, zorder=2)
    for x, label, color in [
        (0.080, "gene input", COLORS["blue"]),
        (0.215, "pathway summary", COLORS["teal"]),
        (0.355, "held-out task", COLORS["muted"]),
    ]:
        ax.plot([x, x], [0.166, 0.205], color=color, alpha=0.55, linewidth=0.9, zorder=2)
        ax.text(x, 0.135, label, color=COLORS["muted"], fontsize=6.3, ha="center", alpha=0.72, zorder=2)
    for x, y, r, color in [
        (0.105, 0.812, 0.0045, COLORS["blue"]),
        (0.140, 0.842, 0.0034, COLORS["blue"]),
        (0.188, 0.817, 0.0038, COLORS["teal"]),
        (0.238, 0.848, 0.0042, COLORS["teal"]),
        (0.302, 0.812, 0.0032, COLORS["muted"]),
    ]:
        ax.add_patch(Ellipse((x, y), r, r * 1.25, transform=ax.transAxes, color=color, alpha=0.18, zorder=2))


def draw_model_measurement_layer(ax: plt.Axes) -> None:
    ax.plot([0.070, 0.372], [0.185, 0.185], color=COLORS["grid"], alpha=0.24, linewidth=1.0, zorder=2)
    for x, label, color in [
        (0.083, "classical", COLORS["blue"]),
        (0.220, "single-cell", COLORS["purple"]),
        (0.360, "text LLM", "#C026D3"),
    ]:
        ax.plot([x, x], [0.166, 0.205], color=color, alpha=0.55, linewidth=0.9, zorder=2)
        ax.text(x, 0.135, label, color=COLORS["muted"], fontsize=6.3, ha="center", alpha=0.72, zorder=2)
    for x, y, w, color in [
        (0.130, 0.817, 0.034, COLORS["blue"]),
        (0.205, 0.842, 0.026, COLORS["purple"]),
        (0.286, 0.818, 0.031, "#C026D3"),
        (0.345, 0.846, 0.020, COLORS["muted"]),
    ]:
        ax.plot([x, x + w], [y, y + 0.004], color=color, alpha=0.14, linewidth=1.0, zorder=2)


def draw_organoid_measurement_layer(ax: plt.Axes) -> None:
    ax.plot([0.070, 0.370], [0.185, 0.185], color=COLORS["grid"], alpha=0.24, linewidth=1.0, zorder=2)
    for x, label, color in [
        (0.083, "2 sources", COLORS["amber"]),
        (0.220, "42 samples", COLORS["amber"]),
        (0.360, "gene checks", COLORS["green"]),
    ]:
        ax.plot([x, x], [0.166, 0.205], color=color, alpha=0.55, linewidth=0.9, zorder=2)
        ax.text(x, 0.135, label, color=COLORS["muted"], fontsize=6.3, ha="center", alpha=0.72, zorder=2)
    for x, y, r, color in [
        (0.116, 0.820, 0.009, COLORS["green"]),
        (0.143, 0.844, 0.006, COLORS["green"]),
        (0.178, 0.813, 0.004, COLORS["amber"]),
        (0.246, 0.843, 0.006, COLORS["blue"]),
        (0.294, 0.816, 0.004, COLORS["green"]),
    ]:
        ax.add_patch(Ellipse((x, y), r, r * 1.2, transform=ax.transAxes, fill=False, color=color, alpha=0.16, linewidth=1.0, zorder=2))


def draw_proof_plate(ax: plt.Axes, image: np.ndarray, extent: tuple[float, float, float, float]) -> None:
    x0, x1, y0, y1 = extent
    ax.add_patch(
        Rectangle(
            (x0 + 0.014, y0 - 0.018),
            x1 - x0,
            y1 - y0,
            transform=ax.transAxes,
            color=COLORS["shadow"],
            alpha=0.36,
            linewidth=0,
            zorder=3,
        )
    )
    ax.add_patch(
        Rectangle(
            (x0 - 0.006, y0 - 0.006),
            (x1 - x0) + 0.012,
            (y1 - y0) + 0.012,
            transform=ax.transAxes,
            facecolor=COLORS["paper"],
            edgecolor="#D8DEE6",
            linewidth=0.8,
            zorder=4,
        )
    )
    ax.imshow(image, extent=extent, origin="upper", zorder=5)


def plate_point(extent: tuple[float, float, float, float], x: float, y: float) -> tuple[float, float]:
    x0, x1, y0, y1 = extent
    return x0 + x * (x1 - x0), y0 + y * (y1 - y0)


def draw_focus_window(
    ax: plt.Axes,
    extent: tuple[float, float, float, float],
    x0: float,
    x1: float,
    y0: float,
    y1: float,
    color: str,
    *,
    alpha: float = 0.78,
) -> None:
    ax0, ay0 = plate_point(extent, x0, y0)
    ax1, ay1 = plate_point(extent, x1, y1)
    left, right = sorted([ax0, ax1])
    bottom, top = sorted([ay0, ay1])
    width = right - left
    height = top - bottom
    corner = min(width, height) * 0.22
    ax.add_patch(
        Rectangle(
            (left, bottom),
            width,
            height,
            transform=ax.transAxes,
            facecolor=color,
            edgecolor="none",
            alpha=0.045,
            zorder=8,
        )
    )
    segments = [
        ([left, left + corner], [top, top]),
        ([left, left], [top - corner, top]),
        ([right - corner, right], [top, top]),
        ([right, right], [top - corner, top]),
        ([left, left + corner], [bottom, bottom]),
        ([left, left], [bottom, bottom + corner]),
        ([right - corner, right], [bottom, bottom]),
        ([right, right], [bottom, bottom + corner]),
    ]
    for xs, ys in segments:
        ax.plot(xs, ys, color=color, alpha=alpha, linewidth=1.25, transform=ax.transAxes, zorder=9)


def draw_interpretation_overlay(ax: plt.Axes, proof_extent: tuple[float, float, float, float]) -> None:
    ax.text(
        0.066,
        0.805,
        "Some tissues\ntransfer",
        color=COLORS["ink"],
        fontsize=24,
        fontweight="bold",
        ha="left",
        va="top",
        linespacing=1.03,
        zorder=8,
    )
    ax.text(
        0.068,
        0.635,
        "Held-out missions reveal\nstrong vs near-chance transfer.",
        color=COLORS["muted"],
        fontsize=10.5,
        ha="left",
        va="top",
        linespacing=1.22,
        zorder=8,
    )

    ax.text(
        0.071,
        0.468,
        "Thymus and gastrocnemius",
        color=COLORS["amber"],
        fontsize=10.5,
        fontweight="bold",
        ha="left",
        zorder=8,
    )
    ax.text(0.071, 0.435, "retain transfer signal", color=COLORS["ink"], fontsize=9.2, ha="left", zorder=8)
    ax.text(
        0.071,
        0.358,
        "Liver and kidney",
        color=COLORS["muted"],
        fontsize=10.5,
        fontweight="bold",
        ha="left",
        zorder=8,
    )
    ax.text(0.071, 0.325, "sit close to chance", color=COLORS["ink"], fontsize=9.2, ha="left", zorder=8)

    high = plate_point(proof_extent, 0.735, 0.760)
    low = plate_point(proof_extent, 0.430, 0.255)
    ax.add_patch(
        Ellipse(high, 0.082, 0.168, transform=ax.transAxes, fill=False, color=COLORS["amber"], linewidth=1.7, alpha=0.88, zorder=9)
    )
    ax.add_patch(
        Ellipse(high, 0.108, 0.212, transform=ax.transAxes, fill=False, color=COLORS["amber"], linewidth=0.75, alpha=0.25, zorder=9)
    )
    ax.add_patch(
        Ellipse(low, 0.092, 0.145, transform=ax.transAxes, fill=False, color=COLORS["muted"], linewidth=1.7, alpha=0.82, zorder=9)
    )
    ax.add_patch(
        Ellipse(low, 0.120, 0.188, transform=ax.transAxes, fill=False, color=COLORS["muted"], linewidth=0.75, alpha=0.22, zorder=9)
    )

    arrow = FancyArrowPatch(
        plate_point(proof_extent, 0.690, 0.650),
        plate_point(proof_extent, 0.485, 0.360),
        arrowstyle="-|>",
        mutation_scale=12,
        linewidth=1.2,
        color=COLORS["amber"],
        alpha=0.62,
        transform=ax.transAxes,
        zorder=9,
        connectionstyle="arc3,rad=-0.18",
    )
    ax.add_patch(arrow)
    ax.text(
        0.402,
        0.852,
        "audited proof figure",
        color=COLORS["muted"],
        fontsize=7.2,
        ha="left",
        alpha=0.78,
        zorder=8,
    )
    ax.text(
        0.066,
        0.075,
        "Mouse bulk RNA-seq | one-mission-held-out AUROC",
        color=COLORS["muted"],
        fontsize=8.0,
        ha="left",
        alpha=0.82,
        zorder=8,
    )


def fig1_overlay_spec() -> dict[str, object]:
    return {
        "slide_id": "fig1_tissue_transfer_layered_scene_pilot",
        "editable_overlay_elements": [
            {
                "layer": "Z3",
                "type": "title",
                "text": "Some tissues transfer",
                "purpose": "main interpretive claim",
            },
            {
                "layer": "Z3",
                "type": "callout",
                "text": "Held-out missions reveal strong vs near-chance transfer.",
                "purpose": "plain-language guide",
            },
            {
                "layer": "Z3",
                "type": "focus_ring",
                "targets": ["high-transfer tissue group", "near-chance tissue group"],
                "purpose": "direct attention to transfer break",
            },
            {
                "layer": "Z4",
                "type": "scope_label",
                "text": "Mouse bulk RNA-seq | one-mission-held-out AUROC",
                "purpose": "trust and scope boundary",
            },
            {
                "layer": "Z5",
                "type": "trajectory",
                "text": "one transfer-break arrow",
                "purpose": "single visual movement",
            },
        ],
    }


def write_fig1_manifest(outputs: dict[str, Path]) -> Path:
    manifest = {
        "slide_id": "fig1_tissue_transfer_layered_scene_pilot",
        "created": "2026-06-01",
        "design_rule": "Layered PNG scene + editable scientific interpretation",
        "source_proof_object": str(FIG1_PROOF.relative_to(ROOT)),
        "source_manifest": str(FIG1_MANIFEST.relative_to(ROOT)),
        "numeric_audit": str(NUMERIC_AUDIT.relative_to(ROOT)),
        "layers": [
            {"layer": "Z0", "name": "canvas_atmosphere", "implementation": "procedural dark transcriptomic texture"},
            {"layer": "Z1", "name": "measurement_layer", "implementation": "held-out mission arc and train/held-out ticks"},
            {"layer": "Z2", "name": "evidence_surface", "implementation": "source-verified Fig1 proof object with backplate"},
            {"layer": "Z3", "name": "interpretation_layer", "implementation": "editable-equivalent title, callout, and focus rings"},
            {"layer": "Z4", "name": "trust_status_layer", "implementation": "short scope label; no raw paths in visible slide"},
            {"layer": "Z5", "name": "motion_focus_layer", "implementation": "single transfer-break arrow"},
        ],
        "claim_boundary": "Mouse bulk GeneLab/OSDR one-mission-held-out tissue-transfer benchmark; not a universal tissue or model claim.",
        "outputs": {key: str(path.relative_to(ROOT)) for key, path in outputs.items()},
    }
    path = OUT_DIR / "fig1_tissue_transfer_layered_scene_manifest.json"
    path.write_text(json.dumps(manifest, indent=2) + "\n")
    return path


def draw_fig2_interpretation_overlay(
    ax: plt.Axes,
    panel_a_extent: tuple[float, float, float, float],
    panel_b_extent: tuple[float, float, float, float],
) -> None:
    ax.text(
        0.066,
        0.810,
        "Pathways suppress\nunwanted labels",
        color=COLORS["ink"],
        fontsize=21.5,
        fontweight="bold",
        ha="left",
        va="top",
        linespacing=1.02,
        zorder=8,
    )
    ax.text(
        0.068,
        0.645,
        "Mission, hardware, and gravity labels\nbecome weaker after pathway summarization.",
        color=COLORS["muted"],
        fontsize=9.7,
        ha="left",
        va="top",
        linespacing=1.22,
        zorder=8,
    )
    ax.text(
        0.071,
        0.488,
        "Signal reduced",
        color=COLORS["teal"],
        fontsize=10.6,
        fontweight="bold",
        ha="left",
        zorder=8,
    )
    ax.text(0.071, 0.456, "mission, hardware, gravity labels", color=COLORS["ink"], fontsize=9.1, ha="left", zorder=8)
    ax.text(
        0.071,
        0.360,
        "Selective gains",
        color=COLORS["amber"],
        fontsize=10.6,
        fontweight="bold",
        ha="left",
        zorder=8,
    )
    ax.text(0.071, 0.328, "kidney and eye improve", color=COLORS["ink"], fontsize=9.1, ha="left", zorder=8)

    # Panel A focus: paired drops from gene-level to pathway-level inputs.
    draw_focus_window(ax, panel_a_extent, 0.855, 0.982, 0.178, 0.855, COLORS["blue"], alpha=0.66)
    draw_focus_window(ax, panel_a_extent, 0.145, 0.425, 0.178, 0.855, COLORS["teal"], alpha=0.64)
    ax.add_patch(
        FancyArrowPatch(
            plate_point(panel_a_extent, 0.810, 0.205),
            plate_point(panel_a_extent, 0.405, 0.205),
            arrowstyle="-|>",
            mutation_scale=10,
            linewidth=1.05,
            color=COLORS["teal"],
            alpha=0.58,
            transform=ax.transAxes,
            zorder=9,
            connectionstyle="arc3,rad=0.08",
        )
    )

    # Panel B focus: selected pathway-minus-gene improvements.
    draw_focus_window(ax, panel_b_extent, 0.545, 0.980, 0.515, 0.760, COLORS["amber"], alpha=0.78)
    ax.text(
        0.452,
        0.852,
        "audited proof crops",
        color=COLORS["muted"],
        fontsize=7.2,
        ha="left",
        alpha=0.78,
        zorder=8,
    )
    ax.text(
        0.066,
        0.075,
        "Diagnostic check: mission, hardware, gravity labels",
        color=COLORS["muted"],
        fontsize=8.0,
        ha="left",
        alpha=0.82,
        zorder=8,
    )


def fig2_overlay_spec() -> dict[str, object]:
    return {
        "slide_id": "fig2_pathway_layered_scene_pilot",
        "editable_overlay_elements": [
            {
                "layer": "Z3",
                "type": "title",
                "text": "Pathways suppress unwanted labels",
                "purpose": "main interpretive claim",
            },
            {
                "layer": "Z3",
                "type": "callout",
                "text": "Mission, hardware, and gravity labels become weaker after pathway summarization.",
                "purpose": "plain-language guide",
            },
            {
                "layer": "Z3",
                "type": "focus_window",
                "targets": ["gene endpoints", "pathway endpoints", "selected detection gains"],
                "purpose": "direct attention to suppression and selected improvement",
            },
            {
                "layer": "Z4",
                "type": "scope_label",
                "text": "Diagnostic check: mission, hardware, gravity labels",
                "purpose": "trust and scope boundary",
            },
            {
                "layer": "Z5",
                "type": "trajectory",
                "text": "gene-level labels to pathway summary",
                "purpose": "single visual movement",
            },
        ],
    }


def write_fig2_manifest(outputs: dict[str, Path]) -> Path:
    manifest = {
        "slide_id": "fig2_pathway_layered_scene_pilot",
        "created": "2026-06-01",
        "design_rule": "Layered PNG scene + editable scientific interpretation",
        "source_proof_object": str(FIG2_PROOF.relative_to(ROOT)),
        "source_manifest": str(FIG2_MANIFEST.relative_to(ROOT)),
        "numeric_audit": str(NUMERIC_AUDIT.relative_to(ROOT)),
        "layers": [
            {"layer": "Z0", "name": "canvas_atmosphere", "implementation": "procedural dark transcriptomic texture"},
            {"layer": "Z1", "name": "measurement_layer", "implementation": "gene-to-pathway assay lane"},
            {"layer": "Z2", "name": "evidence_surface", "implementation": "source-derived Fig2 Panel A and Panel B crops with backplates"},
            {"layer": "Z3", "name": "interpretation_layer", "implementation": "editable-equivalent title, callout, and focus rings"},
            {"layer": "Z4", "name": "trust_status_layer", "implementation": "short diagnostic-scope label; no raw paths in visible slide"},
            {"layer": "Z5", "name": "motion_focus_layer", "implementation": "single gene-to-pathway transition arrow"},
        ],
        "claim_boundary": "Diagnostic checks on coupled mission, hardware, and gravity labels; selected pathway gains are tissue/task-specific.",
        "outputs": {key: str(path.relative_to(ROOT)) for key, path in outputs.items()},
    }
    path = OUT_DIR / "fig2_pathway_layered_scene_manifest.json"
    path.write_text(json.dumps(manifest, indent=2) + "\n")
    return path


def draw_fig3_interpretation_overlay(
    ax: plt.Axes,
    panel_a_extent: tuple[float, float, float, float],
    panel_b_extent: tuple[float, float, float, float],
) -> None:
    ax.text(
        0.066,
        0.810,
        "Scale alone\ndoes not transfer",
        color=COLORS["ink"],
        fontsize=22.2,
        fontweight="bold",
        ha="left",
        va="top",
        linespacing=1.02,
        zorder=8,
    )
    ax.text(
        0.068,
        0.645,
        "On shared tasks, PCA-LR remains highest;\nmodel gains are local and tissue-specific.",
        color=COLORS["muted"],
        fontsize=9.7,
        ha="left",
        va="top",
        linespacing=1.22,
        zorder=8,
    )
    ax.text(
        0.071,
        0.488,
        "Shared task rows",
        color=COLORS["blue"],
        fontsize=10.6,
        fontweight="bold",
        ha="left",
        zorder=8,
    )
    ax.text(0.071, 0.456, "direct comparison only", color=COLORS["ink"], fontsize=9.1, ha="left", zorder=8)
    ax.text(
        0.071,
        0.360,
        "Tissue deltas",
        color=COLORS["violet"],
        fontsize=10.6,
        fontweight="bold",
        ha="left",
        zorder=8,
    )
    ax.text(0.071, 0.328, "mostly below matched baseline", color=COLORS["ink"], fontsize=9.1, ha="left", zorder=8)

    # Panel A focus: keep attention on the shared-task comparison rather than a leaderboard.
    draw_focus_window(ax, panel_a_extent, 0.805, 0.970, 0.660, 0.900, COLORS["blue"], alpha=0.70)
    draw_focus_window(ax, panel_a_extent, 0.575, 0.760, 0.500, 0.730, COLORS["purple"], alpha=0.62)

    # Panel B focus: foundation-model deltas are usually below their matched classical baseline.
    draw_focus_window(ax, panel_b_extent, 0.170, 0.735, 0.145, 0.785, COLORS["violet"], alpha=0.58)
    ax.add_patch(
        FancyArrowPatch(
            (0.666, 0.512),
            (0.666, 0.480),
            arrowstyle="-|>",
            mutation_scale=10,
            linewidth=1.05,
            color=COLORS["violet"],
            alpha=0.50,
            transform=ax.transAxes,
            zorder=9,
            connectionstyle="arc3,rad=-0.08",
        )
    )
    ax.text(
        0.452,
        0.852,
        "audited proof crops",
        color=COLORS["muted"],
        fontsize=7.2,
        ha="left",
        alpha=0.78,
        zorder=8,
    )
    ax.text(
        0.066,
        0.075,
        "Shared 6-task rows are direct comparisons",
        color=COLORS["muted"],
        fontsize=8.0,
        ha="left",
        alpha=0.82,
        zorder=8,
    )


def fig3_overlay_spec() -> dict[str, object]:
    return {
        "slide_id": "fig3_model_tier_layered_scene_pilot",
        "editable_overlay_elements": [
            {
                "layer": "Z3",
                "type": "title",
                "text": "Scale alone does not transfer",
                "purpose": "main interpretive claim",
            },
            {
                "layer": "Z3",
                "type": "callout",
                "text": "On shared tasks, PCA-LR remains highest; model gains are local and tissue-specific.",
                "purpose": "plain-language guide",
            },
            {
                "layer": "Z3",
                "type": "focus_window",
                "targets": ["shared-task PCA-LR row", "single-cell rows", "negative tissue deltas"],
                "purpose": "direct attention to comparison boundary and tissue-specific failures",
            },
            {
                "layer": "Z4",
                "type": "scope_label",
                "text": "Shared 6-task rows are direct comparisons",
                "purpose": "trust and scope boundary",
            },
            {
                "layer": "Z5",
                "type": "trajectory",
                "text": "aggregate shared-task comparison to tissue deltas",
                "purpose": "single visual movement",
            },
        ],
    }


def write_fig3_manifest(outputs: dict[str, Path]) -> Path:
    manifest = {
        "slide_id": "fig3_model_tier_layered_scene_pilot",
        "created": "2026-06-01",
        "design_rule": "Layered PNG scene + editable scientific interpretation",
        "source_proof_object": str(FIG3_PROOF.relative_to(ROOT)),
        "source_manifest": str(FIG3_MANIFEST.relative_to(ROOT)),
        "numeric_audit": str(NUMERIC_AUDIT.relative_to(ROOT)),
        "layers": [
            {"layer": "Z0", "name": "canvas_atmosphere", "implementation": "procedural dark benchmark texture"},
            {"layer": "Z1", "name": "measurement_layer", "implementation": "quiet model-family lane"},
            {"layer": "Z2", "name": "evidence_surface", "implementation": "source-derived Fig3 Panel A and Panel B crops with backplates"},
            {"layer": "Z3", "name": "interpretation_layer", "implementation": "editable-equivalent title, callout, and corner focus windows"},
            {"layer": "Z4", "name": "trust_status_layer", "implementation": "short shared-row comparability label; no raw paths in visible slide"},
            {"layer": "Z5", "name": "motion_focus_layer", "implementation": "single transition from aggregate comparison to tissue deltas"},
        ],
        "claim_boundary": "Only shared six-task rows are direct model comparisons; tissue deltas are matched local comparisons, not a universal model ranking.",
        "outputs": {key: str(path.relative_to(ROOT)) for key, path in outputs.items()},
    }
    path = OUT_DIR / "fig3_model_tier_layered_scene_manifest.json"
    path.write_text(json.dumps(manifest, indent=2) + "\n")
    return path


def draw_fig6_interpretation_overlay(
    ax: plt.Axes,
    footprint_extent: tuple[float, float, float, float],
    checks_extent: tuple[float, float, float, float],
) -> None:
    ax.text(
        0.066,
        0.810,
        "Organoids add\nbiology checks",
        color=COLORS["ink"],
        fontsize=22.2,
        fontweight="bold",
        ha="left",
        va="top",
        linespacing=1.02,
        zorder=8,
    )
    ax.text(
        0.068,
        0.645,
        "Two public sources and 42 samples support\ndraft diagnostics, not validation.",
        color=COLORS["muted"],
        fontsize=9.7,
        ha="left",
        va="top",
        linespacing=1.22,
        zorder=8,
    )
    ax.text(
        0.071,
        0.488,
        "Small public set",
        color=COLORS["amber"],
        fontsize=10.6,
        fontweight="bold",
        ha="left",
        zorder=8,
    )
    ax.text(0.071, 0.456, "2 sources, 42 samples", color=COLORS["ink"], fontsize=9.1, ha="left", zorder=8)
    ax.text(
        0.071,
        0.360,
        "Separate evidence",
        color=COLORS["green"],
        fontsize=10.6,
        fontweight="bold",
        ha="left",
        zorder=8,
    )
    ax.text(0.071, 0.328, "prediction plus gene-response checks", color=COLORS["ink"], fontsize=9.1, ha="left", zorder=8)

    draw_focus_window(ax, footprint_extent, 0.045, 0.405, 0.330, 0.850, COLORS["amber"], alpha=0.70)
    draw_focus_window(ax, checks_extent, 0.060, 0.950, 0.075, 0.470, COLORS["green"], alpha=0.60)
    ax.add_patch(
        FancyArrowPatch(
            (0.666, 0.648),
            (0.666, 0.616),
            arrowstyle="-|>",
            mutation_scale=10,
            linewidth=1.05,
            color=COLORS["green"],
            alpha=0.50,
            transform=ax.transAxes,
            zorder=9,
            connectionstyle="arc3,rad=-0.08",
        )
    )
    ax.text(
        0.430,
        0.855,
        "audited proof crops",
        color=COLORS["muted"],
        fontsize=7.2,
        ha="left",
        alpha=0.78,
        zorder=8,
    )
    ax.text(
        0.066,
        0.075,
        "Draft extension: source factors remain coupled",
        color=COLORS["muted"],
        fontsize=8.0,
        ha="left",
        alpha=0.82,
        zorder=8,
    )


def fig6_overlay_spec() -> dict[str, object]:
    return {
        "slide_id": "fig6_organoid_layered_scene_pilot",
        "editable_overlay_elements": [
            {
                "layer": "Z3",
                "type": "title",
                "text": "Organoids add biology checks",
                "purpose": "main interpretive claim",
            },
            {
                "layer": "Z3",
                "type": "callout",
                "text": "Two public sources and 42 samples support draft diagnostics, not validation.",
                "purpose": "plain-language guide",
            },
            {
                "layer": "Z3",
                "type": "focus_window",
                "targets": ["small public dataset footprint", "gene-response overlap checks"],
                "purpose": "direct attention to scope and secondary biology evidence",
            },
            {
                "layer": "Z4",
                "type": "scope_label",
                "text": "Draft extension: source factors remain coupled",
                "purpose": "trust and scope boundary",
            },
            {
                "layer": "Z5",
                "type": "trajectory",
                "text": "source footprint to biology checks",
                "purpose": "single visual movement",
            },
        ],
    }


def write_fig6_manifest(outputs: dict[str, Path]) -> Path:
    manifest = {
        "slide_id": "fig6_organoid_layered_scene_pilot",
        "created": "2026-06-01",
        "design_rule": "Layered PNG scene + editable scientific interpretation",
        "source_proof_object": str(FIG6_PROOF.relative_to(ROOT)),
        "source_manifest": str(FIG6_MANIFEST.relative_to(ROOT)),
        "numeric_audit": str(NUMERIC_AUDIT.relative_to(ROOT)),
        "layers": [
            {"layer": "Z0", "name": "canvas_atmosphere", "implementation": "procedural dark cellular texture"},
            {"layer": "Z1", "name": "measurement_layer", "implementation": "quiet source/sample/gene-check lane"},
            {"layer": "Z2", "name": "evidence_surface", "implementation": "source-derived Fig6 footprint and biology-check crops with backplates"},
            {"layer": "Z3", "name": "interpretation_layer", "implementation": "editable-equivalent title, callout, and corner focus windows"},
            {"layer": "Z4", "name": "trust_status_layer", "implementation": "draft-extension caution label; no raw paths in visible slide"},
            {"layer": "Z5", "name": "motion_focus_layer", "implementation": "single transition from source footprint to biology checks"},
        ],
        "claim_boundary": "Human organoid evidence is a small draft diagnostic extension; source, disease, donor, and microglia factors remain coupled.",
        "outputs": {key: str(path.relative_to(ROOT)) for key, path in outputs.items()},
    }
    path = OUT_DIR / "fig6_organoid_layered_scene_manifest.json"
    path.write_text(json.dumps(manifest, indent=2) + "\n")
    return path


def build_fig1_scene() -> dict[str, Path]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    proof = crop_fig1_plot(load_rgb(FIG1_PROOF))
    proof_extent = (0.435, 0.965, 0.272, 0.685)

    rng = np.random.default_rng(20260601)
    fig = plt.figure(figsize=(16, 9), dpi=SLIDE_DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    draw_canvas(ax, rng)
    draw_measurement_layer(ax)
    draw_proof_plate(ax, proof, proof_extent)

    plate_path = OUT_DIR / "fig1_tissue_transfer_scene_plate.png"
    fig.savefig(plate_path, dpi=SLIDE_DPI, facecolor=COLORS["bg"])

    draw_interpretation_overlay(ax, proof_extent)
    final_png = OUT_DIR / "fig1_tissue_transfer_layered_scene.png"
    final_pdf = OUT_DIR / "fig1_tissue_transfer_layered_scene.pdf"
    fig.savefig(final_png, dpi=SLIDE_DPI, facecolor=COLORS["bg"])
    fig.savefig(final_pdf, facecolor=COLORS["bg"])
    plt.close(fig)

    overlay_path = OUT_DIR / "fig1_tissue_transfer_editable_overlay_spec.json"
    overlay_path.write_text(json.dumps(fig1_overlay_spec(), indent=2) + "\n")

    outputs = {
        "scene_plate_png": plate_path,
        "final_composite_png": final_png,
        "final_composite_pdf": final_pdf,
        "editable_overlay_spec": overlay_path,
    }
    outputs["layer_manifest"] = write_fig1_manifest(outputs)
    return outputs


def build_fig2_scene() -> dict[str, Path]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    proof = load_rgb(FIG2_PROOF)
    panel_a = crop_fraction(proof, 0.035, 0.990, 0.158, 0.555)
    panel_b = crop_fraction(proof, 0.055, 0.510, 0.575, 0.936)
    panel_a_extent = (0.452, 0.960, 0.535, 0.790)
    panel_b_extent = (0.552, 0.818, 0.218, 0.430)

    rng = np.random.default_rng(20260602)
    fig = plt.figure(figsize=(16, 9), dpi=SLIDE_DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    draw_canvas(ax, rng)
    draw_pathway_measurement_layer(ax)
    draw_proof_plate(ax, panel_a, panel_a_extent)
    draw_proof_plate(ax, panel_b, panel_b_extent)

    plate_path = OUT_DIR / "fig2_pathway_scene_plate.png"
    fig.savefig(plate_path, dpi=SLIDE_DPI, facecolor=COLORS["bg"])

    draw_fig2_interpretation_overlay(ax, panel_a_extent, panel_b_extent)
    final_png = OUT_DIR / "fig2_pathway_layered_scene.png"
    final_pdf = OUT_DIR / "fig2_pathway_layered_scene.pdf"
    fig.savefig(final_png, dpi=SLIDE_DPI, facecolor=COLORS["bg"])
    fig.savefig(final_pdf, facecolor=COLORS["bg"])
    plt.close(fig)

    overlay_path = OUT_DIR / "fig2_pathway_editable_overlay_spec.json"
    overlay_path.write_text(json.dumps(fig2_overlay_spec(), indent=2) + "\n")

    outputs = {
        "scene_plate_png": plate_path,
        "final_composite_png": final_png,
        "final_composite_pdf": final_pdf,
        "editable_overlay_spec": overlay_path,
    }
    outputs["layer_manifest"] = write_fig2_manifest(outputs)
    return outputs


def build_fig3_scene() -> dict[str, Path]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    proof = load_rgb(FIG3_PROOF)
    panel_a = crop_fraction(proof, 0.000, 0.998, 0.165, 0.555)
    panel_b = crop_fraction(proof, 0.000, 0.998, 0.565, 0.955)
    panel_a_extent = (0.425, 0.962, 0.525, 0.805)
    panel_b_extent = (0.425, 0.962, 0.180, 0.468)

    rng = np.random.default_rng(20260603)
    fig = plt.figure(figsize=(16, 9), dpi=SLIDE_DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    draw_canvas(ax, rng)
    draw_model_measurement_layer(ax)
    draw_proof_plate(ax, panel_a, panel_a_extent)
    draw_proof_plate(ax, panel_b, panel_b_extent)

    plate_path = OUT_DIR / "fig3_model_tier_scene_plate.png"
    fig.savefig(plate_path, dpi=SLIDE_DPI, facecolor=COLORS["bg"])

    draw_fig3_interpretation_overlay(ax, panel_a_extent, panel_b_extent)
    final_png = OUT_DIR / "fig3_model_tier_layered_scene.png"
    final_pdf = OUT_DIR / "fig3_model_tier_layered_scene.pdf"
    fig.savefig(final_png, dpi=SLIDE_DPI, facecolor=COLORS["bg"])
    fig.savefig(final_pdf, facecolor=COLORS["bg"])
    plt.close(fig)

    overlay_path = OUT_DIR / "fig3_model_tier_editable_overlay_spec.json"
    overlay_path.write_text(json.dumps(fig3_overlay_spec(), indent=2) + "\n")

    outputs = {
        "scene_plate_png": plate_path,
        "final_composite_png": final_png,
        "final_composite_pdf": final_pdf,
        "editable_overlay_spec": overlay_path,
    }
    outputs["layer_manifest"] = write_fig3_manifest(outputs)
    return outputs


def build_fig6_scene() -> dict[str, Path]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    proof = load_rgb(FIG6_PROOF)
    footprint = crop_fraction(proof, 0.205, 0.980, 0.180, 0.335)
    checks = crop_fraction(proof, 0.050, 0.998, 0.390, 0.980)
    footprint_extent = (0.430, 0.962, 0.680, 0.830)
    checks_extent = (0.430, 0.962, 0.182, 0.622)

    rng = np.random.default_rng(20260606)
    fig = plt.figure(figsize=(16, 9), dpi=SLIDE_DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    draw_canvas(ax, rng)
    draw_organoid_measurement_layer(ax)
    draw_proof_plate(ax, footprint, footprint_extent)
    draw_proof_plate(ax, checks, checks_extent)

    plate_path = OUT_DIR / "fig6_organoid_scene_plate.png"
    fig.savefig(plate_path, dpi=SLIDE_DPI, facecolor=COLORS["bg"])

    draw_fig6_interpretation_overlay(ax, footprint_extent, checks_extent)
    final_png = OUT_DIR / "fig6_organoid_layered_scene.png"
    final_pdf = OUT_DIR / "fig6_organoid_layered_scene.pdf"
    fig.savefig(final_png, dpi=SLIDE_DPI, facecolor=COLORS["bg"])
    fig.savefig(final_pdf, facecolor=COLORS["bg"])
    plt.close(fig)

    overlay_path = OUT_DIR / "fig6_organoid_editable_overlay_spec.json"
    overlay_path.write_text(json.dumps(fig6_overlay_spec(), indent=2) + "\n")

    outputs = {
        "scene_plate_png": plate_path,
        "final_composite_png": final_png,
        "final_composite_pdf": final_pdf,
        "editable_overlay_spec": overlay_path,
    }
    outputs["layer_manifest"] = write_fig6_manifest(outputs)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scene", choices=["fig1", "fig2", "fig3", "fig6", "all"], default="all")
    args = parser.parse_args()
    built: dict[str, dict[str, str]] = {}
    if args.scene in {"fig1", "all"}:
        built["fig1"] = {key: str(path.relative_to(ROOT)) for key, path in build_fig1_scene().items()}
    if args.scene in {"fig2", "all"}:
        built["fig2"] = {key: str(path.relative_to(ROOT)) for key, path in build_fig2_scene().items()}
    if args.scene in {"fig3", "all"}:
        built["fig3"] = {key: str(path.relative_to(ROOT)) for key, path in build_fig3_scene().items()}
    if args.scene in {"fig6", "all"}:
        built["fig6"] = {key: str(path.relative_to(ROOT)) for key, path in build_fig6_scene().items()}
    print(json.dumps(built, indent=2))


if __name__ == "__main__":
    main()
