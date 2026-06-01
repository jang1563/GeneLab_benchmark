#!/usr/bin/env python3
"""Build light provenance-document slide scenes for v9 platform/resource slides."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "output" / ".matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, Rectangle


PREMIUM_FIGURE_DIR = ROOT / "output" / "premium_figures"
OUT_DIR = ROOT / "output" / "premium_v9_document_scenes"

FIG4_PROOF = PREMIUM_FIGURE_DIR / "premium_fig4_v9_platform_architecture.png"
FIG4_MANIFEST = PREMIUM_FIGURE_DIR / "manifests" / "premium_fig4_v9_platform_status_manifest.json"
FIG5_PROOF = PREMIUM_FIGURE_DIR / "premium_fig5_public_bulk_release_boundary_schematic.png"
FIG5_MANIFEST = PREMIUM_FIGURE_DIR / "manifests" / "premium_fig5_public_bulk_alpha_boundary_manifest.json"
NUMERIC_AUDIT = ROOT / "output" / "premium_audit" / "premium_visual_numeric_audit.csv"

SLIDE_DPI = 240

COLORS = {
    "paper": "#F7F4EF",
    "paper2": "#ECE7DF",
    "ink": "#1F2933",
    "muted": "#687385",
    "rule": "#B9C1CC",
    "green": "#1F8F5F",
    "amber": "#D99A22",
    "blue": "#2F6C9F",
    "red": "#C81E1E",
    "shadow": "#B7B0A6",
    "white": "#FCFCFB",
}


def load_rgb(path: Path) -> np.ndarray:
    image = mpimg.imread(path).astype(float)
    if image.max() > 1.0:
        image = image / 255.0
    if image.shape[-1] == 4:
        alpha = image[..., 3:4]
        image = image[..., :3] * alpha + (1.0 - alpha)
    return image[..., :3]


def crop_fraction(image: np.ndarray, x0: float, x1: float, y0: float, y1: float) -> np.ndarray:
    h, w = image.shape[:2]
    return image[int(h * y0) : int(h * y1), int(w * x0) : int(w * x1), :]


def draw_document_canvas(ax: plt.Axes) -> None:
    ax.add_patch(Rectangle((0, 0), 1, 1, color=COLORS["paper"], transform=ax.transAxes, zorder=0))
    h, w = 540, 960
    yy, xx = np.mgrid[0:h, 0:w]
    base = np.zeros((h, w, 3), dtype=float)
    base[..., 0] = 0.965
    base[..., 1] = 0.948
    base[..., 2] = 0.918
    grain = np.random.default_rng(20260607).normal(0, 0.006, size=(h, w, 1))
    diagonal = ((xx + yy) / (w + h))[..., None] * np.array([0.012, 0.010, 0.006])
    texture = np.clip(base + grain - diagonal, 0, 1)
    ax.imshow(texture, extent=(0, 1, 0, 1), origin="lower", zorder=0)
    ax.set_aspect("auto")

    for y in np.linspace(0.115, 0.885, 8):
        ax.plot([0.055, 0.945], [y, y], color=COLORS["rule"], alpha=0.18, linewidth=0.75, zorder=1)
    for x in [0.055, 0.395, 0.945]:
        ax.plot([x, x], [0.085, 0.905], color=COLORS["rule"], alpha=0.14, linewidth=0.75, zorder=1)


def draw_provenance_rail(ax: plt.Axes, items: list[tuple[float, str, str]]) -> None:
    y = 0.160
    ax.plot([0.070, 0.360], [y, y], color=COLORS["rule"], alpha=0.85, linewidth=1.0, zorder=2)
    for x, label, color in items:
        ax.plot([x, x], [y - 0.026, y + 0.026], color=color, alpha=0.90, linewidth=1.1, zorder=2)
        ax.text(x, y - 0.062, label, color=COLORS["muted"], fontsize=6.4, ha="center", va="top", zorder=2)


def draw_proof_plate(ax: plt.Axes, image: np.ndarray, extent: tuple[float, float, float, float]) -> None:
    x0, x1, y0, y1 = extent
    ax.add_patch(
        Rectangle(
            (x0 + 0.012, y0 - 0.014),
            x1 - x0,
            y1 - y0,
            transform=ax.transAxes,
            color=COLORS["shadow"],
            alpha=0.25,
            linewidth=0,
            zorder=3,
        )
    )
    ax.add_patch(
        Rectangle(
            (x0 - 0.005, y0 - 0.005),
            (x1 - x0) + 0.010,
            (y1 - y0) + 0.010,
            transform=ax.transAxes,
            facecolor=COLORS["white"],
            edgecolor="#D4D9DF",
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
    alpha: float = 0.76,
) -> None:
    ax0, ay0 = plate_point(extent, x0, y0)
    ax1, ay1 = plate_point(extent, x1, y1)
    left, right = sorted([ax0, ax1])
    bottom, top = sorted([ay0, ay1])
    width = right - left
    height = top - bottom
    corner = min(width, height) * 0.24
    ax.add_patch(
        Rectangle(
            (left, bottom),
            width,
            height,
            transform=ax.transAxes,
            facecolor=color,
            edgecolor="none",
            alpha=0.050,
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


def draw_left_text(
    ax: plt.Axes,
    title: str,
    callout: str,
    blocks: list[tuple[str, str, str]],
    scope: str,
) -> None:
    ax.text(
        0.066,
        0.820,
        title,
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
        callout,
        color=COLORS["muted"],
        fontsize=9.6,
        ha="left",
        va="top",
        linespacing=1.22,
        zorder=8,
    )
    y = 0.490
    for heading, body, color in blocks:
        ax.text(0.071, y, heading, color=color, fontsize=10.6, fontweight="bold", ha="left", zorder=8)
        ax.text(0.071, y - 0.032, body, color=COLORS["ink"], fontsize=9.0, ha="left", zorder=8)
        y -= 0.126
    ax.text(0.066, 0.075, scope, color=COLORS["muted"], fontsize=8.0, ha="left", alpha=0.88, zorder=8)


def platform_overlay_spec() -> dict[str, object]:
    return {
        "slide_id": "v9_platform_provenance_document_scene",
        "grammar": "light provenance-document scene",
        "editable_overlay_elements": [
            {"layer": "Z3", "type": "title", "text": "V9 is a staged evidence resource"},
            {
                "layer": "Z3",
                "type": "callout",
                "text": "Public bulk is strongest; extension lanes remain diagnostic.",
            },
            {"layer": "Z4", "type": "scope_label", "text": "Platform status; not a final benchmark result"},
        ],
    }


def boundary_overlay_spec() -> dict[str, object]:
    return {
        "slide_id": "v9_public_bulk_boundary_document_scene",
        "grammar": "light provenance-document scene",
        "editable_overlay_elements": [
            {"layer": "Z3", "type": "title", "text": "Public bulk is metadata-ready"},
            {
                "layer": "Z3",
                "type": "callout",
                "text": "Data-file mirroring remains a separate release gate.",
            },
            {"layer": "Z4", "type": "scope_label", "text": "Metadata preview only; not payload-frozen"},
        ],
    }


def write_manifest(
    slide_id: str,
    outputs: dict[str, Path],
    source_proof: Path,
    source_manifest: Path,
    claim_boundary: str,
) -> Path:
    manifest = {
        "slide_id": slide_id,
        "created": "2026-06-01",
        "design_rule": "Light provenance-document scene + editable release-boundary interpretation",
        "source_proof_object": str(source_proof.relative_to(ROOT)),
        "source_manifest": str(source_manifest.relative_to(ROOT)),
        "numeric_audit": str(NUMERIC_AUDIT.relative_to(ROOT)),
        "layers": [
            {"layer": "Z0", "name": "canvas_document", "implementation": "off-white evidence-paper texture"},
            {"layer": "Z1", "name": "provenance_measurement", "implementation": "thin provenance rail and document rules"},
            {"layer": "Z2", "name": "evidence_surface", "implementation": "source-derived proof crops on paper evidence plates"},
            {"layer": "Z3", "name": "interpretation_layer", "implementation": "editable-equivalent title, callout, and focus windows"},
            {"layer": "Z4", "name": "trust_status_layer", "implementation": "short release-boundary label"},
            {"layer": "Z5", "name": "motion_focus_layer", "implementation": "single provenance or release-boundary arrow"},
        ],
        "claim_boundary": claim_boundary,
        "outputs": {key: str(path.relative_to(ROOT)) for key, path in outputs.items()},
    }
    path = OUT_DIR / f"{slide_id}_manifest.json"
    path.write_text(json.dumps(manifest, indent=2) + "\n")
    return path


def build_platform_scene() -> dict[str, Path]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    proof = load_rgb(FIG4_PROOF)
    architecture = crop_fraction(proof, 0.060, 0.985, 0.205, 0.565)
    evidence_depth = crop_fraction(proof, 0.060, 0.985, 0.625, 0.930)
    arch_extent = (0.430, 0.960, 0.580, 0.810)
    depth_extent = (0.430, 0.960, 0.205, 0.455)

    fig = plt.figure(figsize=(16, 9), dpi=SLIDE_DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    draw_document_canvas(ax)
    draw_provenance_rail(
        ax,
        [
            (0.085, "sources", COLORS["green"]),
            (0.197, "tasks", COLORS["blue"]),
            (0.305, "checks", COLORS["amber"]),
        ],
    )
    draw_proof_plate(ax, architecture, arch_extent)
    draw_proof_plate(ax, evidence_depth, depth_extent)
    plate_path = OUT_DIR / "v9_platform_provenance_document_scene_plate.png"
    fig.savefig(plate_path, dpi=SLIDE_DPI, facecolor=COLORS["paper"])

    draw_left_text(
        ax,
        "V9 is a staged\nevidence resource",
        "Public bulk is strongest; extension lanes\nremain diagnostic or blocked.",
        [
            ("Metadata indexed", "source, task, checksum, baseline records", COLORS["green"]),
            ("Extension lanes", "organoid and multispecies remain draft", COLORS["amber"]),
            ("Payload gaps", "single-cell and local file copies need staging", COLORS["red"]),
        ],
        "Platform status; not a final benchmark result",
    )
    draw_focus_window(ax, depth_extent, 0.015, 0.455, 0.075, 0.310, COLORS["green"], alpha=0.74)
    draw_focus_window(ax, depth_extent, 0.760, 0.985, 0.045, 0.950, COLORS["red"], alpha=0.68)
    ax.add_patch(
        FancyArrowPatch(
            (0.662, 0.550),
            (0.662, 0.485),
            arrowstyle="-|>",
            mutation_scale=10,
            linewidth=1.05,
            color=COLORS["blue"],
            alpha=0.55,
            transform=ax.transAxes,
            zorder=9,
        )
    )

    final_png = OUT_DIR / "v9_platform_provenance_document_scene.png"
    final_pdf = OUT_DIR / "v9_platform_provenance_document_scene.pdf"
    fig.savefig(final_png, dpi=SLIDE_DPI, facecolor=COLORS["paper"])
    fig.savefig(final_pdf, facecolor=COLORS["paper"])
    plt.close(fig)

    overlay_path = OUT_DIR / "v9_platform_provenance_document_overlay_spec.json"
    overlay_path.write_text(json.dumps(platform_overlay_spec(), indent=2) + "\n")
    outputs = {
        "scene_plate_png": plate_path,
        "final_composite_png": final_png,
        "final_composite_pdf": final_pdf,
        "editable_overlay_spec": overlay_path,
    }
    outputs["layer_manifest"] = write_manifest(
        "v9_platform_provenance_document_scene",
        outputs,
        FIG4_PROOF,
        FIG4_MANIFEST,
        "V9 platform figure is a resource/status overview; extension tracks are not final benchmark results.",
    )
    return outputs


def build_boundary_scene() -> dict[str, Path]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    proof = load_rgb(FIG5_PROOF)
    boundary = crop_fraction(proof, 0.050, 0.985, 0.150, 0.705)
    boundary_extent = (0.415, 0.960, 0.365, 0.740)

    fig = plt.figure(figsize=(16, 9), dpi=SLIDE_DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    draw_document_canvas(ax)
    draw_provenance_rail(
        ax,
        [
            (0.085, "metadata", COLORS["green"]),
            (0.220, "mirror", COLORS["red"]),
            (0.355, "package", COLORS["blue"]),
        ],
    )
    draw_proof_plate(ax, boundary, boundary_extent)
    plate_path = OUT_DIR / "v9_public_bulk_boundary_document_scene_plate.png"
    fig.savefig(plate_path, dpi=SLIDE_DPI, facecolor=COLORS["paper"])

    draw_left_text(
        ax,
        "Public bulk is\nmetadata-ready",
        "Data-file mirroring remains a separate\nrelease gate before payload freeze.",
        [
            ("Available now", "source, task, checksum, baseline metadata", COLORS["green"]),
            ("Still pending", "local data-file copy and hash checks", COLORS["red"]),
            ("Boundary intact", "metadata evidence stays separate", COLORS["blue"]),
        ],
        "Metadata preview only; not payload-frozen",
    )
    draw_focus_window(ax, boundary_extent, 0.035, 0.320, 0.285, 0.735, COLORS["green"], alpha=0.74)
    draw_focus_window(ax, boundary_extent, 0.395, 0.670, 0.285, 0.735, COLORS["red"], alpha=0.72)
    ax.add_patch(
        FancyArrowPatch(
            plate_point(boundary_extent, 0.310, 0.455),
            plate_point(boundary_extent, 0.425, 0.455),
            arrowstyle="-|>",
            mutation_scale=10,
            linewidth=1.05,
            color=COLORS["red"],
            alpha=0.55,
            transform=ax.transAxes,
            zorder=9,
        )
    )

    final_png = OUT_DIR / "v9_public_bulk_boundary_document_scene.png"
    final_pdf = OUT_DIR / "v9_public_bulk_boundary_document_scene.pdf"
    fig.savefig(final_png, dpi=SLIDE_DPI, facecolor=COLORS["paper"])
    fig.savefig(final_pdf, facecolor=COLORS["paper"])
    plt.close(fig)

    overlay_path = OUT_DIR / "v9_public_bulk_boundary_document_overlay_spec.json"
    overlay_path.write_text(json.dumps(boundary_overlay_spec(), indent=2) + "\n")
    outputs = {
        "scene_plate_png": plate_path,
        "final_composite_png": final_png,
        "final_composite_pdf": final_pdf,
        "editable_overlay_spec": overlay_path,
    }
    outputs["layer_manifest"] = write_manifest(
        "v9_public_bulk_boundary_document_scene",
        outputs,
        FIG5_PROOF,
        FIG5_MANIFEST,
        "Public-bulk v9 is metadata-preview ready; local data-file mirroring and hash verification remain separate gates.",
    )
    return outputs


def write_contact_sheet(built: dict[str, dict[str, Path]]) -> Path:
    fig = plt.figure(figsize=(14.5, 8.2), dpi=210)
    grid = fig.add_gridspec(2, 2, left=0.035, right=0.985, top=0.875, bottom=0.060, wspace=0.045, hspace=0.115)
    fig.suptitle(
        "V9 provenance-document scene review",
        x=0.035,
        y=0.955,
        ha="left",
        fontsize=14,
        fontweight="bold",
    )
    order = [
        ("platform", "final_composite_png", "Platform scene"),
        ("boundary", "final_composite_png", "Public-bulk boundary scene"),
        ("platform", "scene_plate_png", "Platform scene plate"),
        ("boundary", "scene_plate_png", "Boundary scene plate"),
    ]
    for idx, (scene_key, output_key, title) in enumerate(order):
        ax = fig.add_subplot(grid[idx // 2, idx % 2])
        ax.imshow(mpimg.imread(built[scene_key][output_key]))
        ax.set_title(title, loc="left", fontsize=8.2, pad=3)
        ax.axis("off")
    output = OUT_DIR / "v9_provenance_document_contact_sheet.png"
    fig.savefig(output, dpi=210, facecolor="white")
    plt.close(fig)
    manifest = {
        "contact_sheet": str(output.relative_to(ROOT)),
        "scenes": {
            scene_key: {key: str(path.relative_to(ROOT)) for key, path in outputs.items()}
            for scene_key, outputs in built.items()
        },
    }
    (OUT_DIR / "v9_provenance_document_contact_sheet_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scene", choices=["platform", "boundary", "all"], default="all")
    args = parser.parse_args()
    built: dict[str, dict[str, Path]] = {}
    if args.scene in {"platform", "all"}:
        built["platform"] = build_platform_scene()
    if args.scene in {"boundary", "all"}:
        built["boundary"] = build_boundary_scene()
    if args.scene == "all":
        built["contact_sheet"] = {"contact_sheet": write_contact_sheet(built)}
    print(json.dumps({scene: {key: str(path.relative_to(ROOT)) for key, path in outputs.items()} for scene, outputs in built.items()}, indent=2))


if __name__ == "__main__":
    main()
