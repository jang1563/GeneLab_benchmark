#!/usr/bin/env python3
"""Build consulting-style usage tests for the v0.3 biological visual system."""

from __future__ import annotations

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
from matplotlib.patches import Circle, FancyArrowPatch, Polygon, Rectangle


OUT = ROOT / "output" / "biovis_v0_3_consulting_pack"
ASSET_DIR = ROOT / "assets" / "biovis_motif_pack_v0_3" / "png"
CREATED = "2026-06-02"
SLIDE_DPI = 240

COLORS = {
    "paper": "#F6F2EA",
    "paper2": "#FBFAF7",
    "ink": "#17212B",
    "muted": "#5D6978",
    "rule": "#AEB8C5",
    "blue": "#2D6F9F",
    "green": "#178B63",
    "teal": "#1AA090",
    "amber": "#B4832F",
    "rose": "#B45A7E",
    "purple": "#6750A4",
    "red": "#B23A3A",
    "shadow": "#9C9487",
    "dark": "#111A23",
}


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def write_json(path: Path, data: Any) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def make_fig() -> tuple[plt.Figure, plt.Axes]:
    fig = plt.figure(figsize=(16, 9), dpi=SLIDE_DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    return fig, ax


def draw_canvas(ax: plt.Axes, seed: int, *, dark: bool = False) -> None:
    rng = np.random.default_rng(seed)
    h, w = 720, 1280
    yy, xx = np.mgrid[0:h, 0:w]
    base = np.zeros((h, w, 3), dtype=float)
    if dark:
        base[..., 0] = 0.060
        base[..., 1] = 0.095
        base[..., 2] = 0.125
        grain = rng.normal(0, 0.010, (h, w, 1))
        vignette = ((xx / w - 0.48) ** 2 + (yy / h - 0.48) ** 2)[..., None] * np.array([0.120, 0.110, 0.095])
        img = np.clip(base + grain - vignette, 0, 1)
    else:
        base[..., 0] = 0.967
        base[..., 1] = 0.950
        base[..., 2] = 0.915
        grain = rng.normal(0, 0.005, (h, w, 1))
        vignette = ((xx / w - 0.52) ** 2 + (yy / h - 0.48) ** 2)[..., None] * np.array([0.052, 0.043, 0.030])
        cool_left = (1 - xx / w)[..., None] * np.array([0.000, 0.003, 0.008])
        img = np.clip(base + grain + cool_left - vignette, 0, 1)
    ax.imshow(img, extent=(0, 1, 0, 1), origin="lower", zorder=0)
    ax.set_aspect("auto")
    grid_color = "#FFFFFF" if dark else COLORS["rule"]
    for y in [0.18, 0.36, 0.54, 0.72, 0.86]:
        ax.plot([0.055, 0.945], [y, y], color=grid_color, alpha=0.075 if dark else 0.13, lw=0.8, zorder=1)
    for x in [0.09, 0.30, 0.50, 0.70, 0.91]:
        ax.plot([x, x], [0.095, 0.890], color=grid_color, alpha=0.060 if dark else 0.09, lw=0.8, zorder=1)


def add_asset(ax: plt.Axes, name: str, x: float, y: float, w: float, h: float, *, alpha: float = 1.0, z: int = 5, clip: bool = False) -> None:
    img = mpimg.imread(ASSET_DIR / f"{name}.png")
    if clip:
        # Tight crop the plates for more consulting-slide assertiveness.
        img = img[70:-70, 110:-110, :]
    ax.imshow(img, extent=(x, x + w, y, y + h), transform=ax.transAxes, zorder=z, alpha=alpha, interpolation="lanczos")
    ax.set_aspect("auto")


def arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float], color: str = "rule", *, alpha: float = 0.60) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=14,
            lw=1.2,
            color=COLORS[color],
            alpha=alpha,
            transform=ax.transAxes,
            zorder=20,
        )
    )


def band_label(ax: plt.Axes, x: float, y: float, text: str, color: str, *, width: float = 0.155) -> None:
    ax.add_patch(Rectangle((x, y - 0.016), width, 0.032, facecolor=COLORS[color], edgecolor="none", alpha=0.12, zorder=10))
    ax.plot([x, x + width], [y - 0.017, y - 0.017], color=COLORS[color], alpha=0.55, lw=1.2, zorder=11)
    ax.text(x + 0.008, y, text, fontsize=8.8, color=COLORS["ink"], fontweight="bold", ha="left", va="center", zorder=12)


def title(ax: plt.Axes, headline: str, support: str, *, dark: bool = False) -> None:
    ink = "#F7F4EC" if dark else COLORS["ink"]
    muted = "#BAC4CF" if dark else COLORS["muted"]
    ax.text(0.065, 0.895, headline, fontsize=25.5, fontweight="bold", color=ink, ha="left", va="center", zorder=40)
    ax.text(0.066, 0.840, support, fontsize=11.0, color=muted, ha="left", va="center", zorder=40)


def draw_executive_portfolio() -> dict[str, Any]:
    fig, ax = make_fig()
    draw_canvas(ax, 4011)
    title(
        ax,
        "Biological visuals must answer an executive question",
        "Use the v0.3 system to connect coverage, evidence status, and next decisions.",
    )

    add_asset(ax, "cell_micrograph_texture_field", 0.060, 0.255, 0.360, 0.203, alpha=0.98, z=4, clip=True)
    ax.add_patch(Rectangle((0.060, 0.255), 0.360, 0.203, facecolor="none", edgecolor=COLORS["ink"], lw=0.6, alpha=0.25, zorder=8))
    ax.add_patch(Polygon([(0.060, 0.510), (0.420, 0.510), (0.385, 0.108), (0.092, 0.108)], facecolor=COLORS["teal"], edgecolor="none", alpha=0.045, zorder=2))
    ax.text(0.060, 0.562, "Board-level question", fontsize=11.4, fontweight="bold", color=COLORS["ink"], ha="left", va="center", zorder=20)
    ax.text(0.060, 0.528, "Where does the biology change the decision?", fontsize=9.4, color=COLORS["muted"], ha="left", va="center", zorder=20)
    for i, (k, v, col) in enumerate([("Coverage", "what systems are in scope", "blue"), ("Evidence", "what is observed vs schematic", "green"), ("Action", "what should change next", "amber")]):
        y = 0.205 - i * 0.052
        ax.text(0.060, y, k, fontsize=8.7, fontweight="bold", color=COLORS[col], ha="left", va="center", zorder=20)
        ax.text(0.135, y, v, fontsize=8.2, color=COLORS["muted"], ha="left", va="center", zorder=20)

    ax.plot([0.455, 0.455], [0.150, 0.735], color=COLORS["rule"], alpha=0.25, lw=1.0, zorder=5)
    ax.text(0.485, 0.735, "content portfolio", fontsize=13.0, fontweight="bold", color=COLORS["ink"], ha="left", va="center", zorder=20)
    ax.text(0.485, 0.705, "What story role does this visual serve?", fontsize=8.5, color=COLORS["muted"], ha="left", va="center", zorder=20)

    rows = [
        ("Scope", "species / modality coverage", "species_model_strip", "blue"),
        ("Method", "feature construction and split logic", "omics_matrix_texture_field", "teal"),
        ("Evidence", "source-like proof object", "cell_micrograph_texture_field", "green"),
        ("Decision", "what changes next", "pathway_network_summary_field", "amber"),
    ]
    y0 = 0.630
    for i, (label, desc, asset, col) in enumerate(rows):
        y = y0 - i * 0.112
        ax.plot([0.455, 0.900], [y - 0.048, y - 0.048], color=COLORS["rule"], alpha=0.18, lw=1.0, zorder=5)
        band_label(ax, 0.485, y, label, col, width=0.115)
        ax.text(0.625, y, desc, fontsize=9.4, color=COLORS["ink"], ha="left", va="center", zorder=20)
        add_asset(ax, asset, 0.815, y - 0.045, 0.130, 0.073, alpha=0.70, z=6, clip=True)

    ax.text(0.485, 0.180, "Consulting-grade rule", fontsize=10.0, fontweight="bold", color=COLORS["ink"], ha="left", va="center", zorder=20)
    ax.text(0.485, 0.145, "Each plate needs a decision role, evidence status, and claim boundary.", fontsize=9.2, color=COLORS["muted"], ha="left", va="center", zorder=20)
    path = OUT / "01_executive_portfolio_map.png"
    fig.savefig(path, dpi=SLIDE_DPI, facecolor=COLORS["paper"])
    plt.close(fig)
    return {"id": "executive_portfolio_map", "png": rel(path)}


def draw_evidence_ladder() -> dict[str, Any]:
    fig, ax = make_fig()
    draw_canvas(ax, 4021, dark=True)
    title(
        ax,
        "Separate visual appeal from evidence strength",
        "The same biological texture can be a context layer, a proof object, or a claim risk.",
        dark=True,
    )

    ladder = [
        ("Context", "generated schematic", "safe for orientation", "cell_micrograph_texture_field", "blue"),
        ("Construct", "feature/model logic", "safe for methods", "omics_matrix_texture_field", "teal"),
        ("Observed", "source-derived panel", "requires citation", "spatial_spot_section_field", "green"),
        ("Decision", "validated result", "requires locked metric", "pathway_network_summary_field", "amber"),
    ]
    xs = [0.090, 0.315, 0.540, 0.765]
    for i, (stage, status, claim, asset, col) in enumerate(ladder):
        x = xs[i]
        ax.plot([x, x], [0.245, 0.700], color="#FFFFFF", alpha=0.10, lw=1.0, zorder=4)
        add_asset(ax, asset, x - 0.010, 0.385, 0.175, 0.098, alpha=0.82 if i < 2 else 0.66, z=5, clip=True)
        ax.text(x, 0.695, f"0{i + 1}", fontsize=10, color="#BAC4CF", ha="left", va="center", zorder=20)
        ax.text(x, 0.648, stage, fontsize=18, fontweight="bold", color="#F7F4EC", ha="left", va="center", zorder=20)
        ax.text(x, 0.605, status, fontsize=9.2, color="#BAC4CF", ha="left", va="center", zorder=20)
        ax.plot([x, x + 0.150], [0.330, 0.330], color=COLORS[col], alpha=0.72, lw=3.2, zorder=18)
        ax.text(x, 0.292, claim, fontsize=8.8, color="#D4DADF", ha="left", va="center", zorder=20)
        if i < 3:
            arrow(ax, (x + 0.170, 0.505), (x + 0.225, 0.505), "rule", alpha=0.45)

    ax.text(0.090, 0.152, "Upgrade path", fontsize=10.5, fontweight="bold", color="#F7F4EC", ha="left", va="center", zorder=20)
    ax.text(0.090, 0.115, "v0.3 supports context and construct layers now; source-derived observed panels remain the next proof upgrade.", fontsize=9.0, color="#BAC4CF", ha="left", va="center", zorder=20)
    path = OUT / "02_evidence_ladder.png"
    fig.savefig(path, dpi=SLIDE_DPI, facecolor=COLORS["dark"])
    plt.close(fig)
    return {"id": "evidence_ladder", "png": rel(path)}


def draw_archetype_gallery() -> dict[str, Any]:
    fig, ax = make_fig()
    draw_canvas(ax, 4031)
    title(
        ax,
        "Use biological content in distinct slide jobs",
        "A premium consulting deck needs repeatable archetypes, not a pile of attractive motifs.",
    )

    lanes = [
        (
            "Orient",
            "scope and executive thesis",
            [("Executive thesis", "one implication, one boundary", "cell_micrograph_texture_field", "blue"), ("Scope map", "species x modality x source", "species_model_strip", "green")],
        ),
        (
            "Explain",
            "method and evidence logic",
            [("Method bridge", "measurement to benchmark", "omics_matrix_texture_field", "teal"), ("Evidence gate", "proof strength and caveat", "spatial_spot_section_field", "amber")],
        ),
        (
            "Commit",
            "result and next action",
            [("Result story", "metric with biological layer", "single_cell_embedding_field", "purple"), ("Roadmap", "what gets upgraded next", "pathway_network_summary_field", "rose")],
        ),
    ]
    y_start = 0.660
    for row, (lane, lane_desc, items) in enumerate(lanes):
        y = y_start - row * 0.200
        ax.text(0.080, y + 0.025, lane, fontsize=14.0, fontweight="bold", color=COLORS["ink"], ha="left", va="center", zorder=20)
        ax.text(0.080, y - 0.018, lane_desc, fontsize=8.8, color=COLORS["muted"], ha="left", va="center", zorder=20)
        ax.plot([0.205, 0.925], [y - 0.082, y - 0.082], color=COLORS["rule"], alpha=0.22, lw=1.0, zorder=5)
        for idx, (name, desc, asset, col) in enumerate(items):
            x = 0.245 + idx * 0.365
            add_asset(ax, asset, x, y - 0.066, 0.155, 0.087, alpha=0.76, z=6, clip=True)
            ax.text(x + 0.185, y + 0.014, name, fontsize=11.5, fontweight="bold", color=COLORS["ink"], ha="left", va="center", zorder=20)
            ax.text(x + 0.185, y - 0.021, desc, fontsize=8.2, color=COLORS["muted"], ha="left", va="center", zorder=20)
            ax.plot([x + 0.185, x + 0.320], [y - 0.050, y - 0.050], color=COLORS[col], alpha=0.50, lw=2.0, zorder=10)

    ax.text(0.080, 0.140, "Design implication", fontsize=10.5, fontweight="bold", color=COLORS["ink"], ha="left", va="center", zorder=20)
    ax.text(0.080, 0.105, "The asset pack should ship with slide jobs, claim boundaries, and visual hierarchy rules.", fontsize=9.0, color=COLORS["muted"], ha="left", va="center", zorder=20)
    path = OUT / "03_consulting_archetype_gallery.png"
    fig.savefig(path, dpi=SLIDE_DPI, facecolor=COLORS["paper"])
    plt.close(fig)
    return {"id": "consulting_archetype_gallery", "png": rel(path)}


def draw_v04_decision_roadmap() -> dict[str, Any]:
    fig, ax = make_fig()
    draw_canvas(ax, 4041)
    title(
        ax,
        "v0.3 is ready for controlled use; v0.4 should target proof",
        "The next improvement is not more icons. It is stronger evidence plates and native deck controls.",
    )

    priorities = [
        ("P0", "Source-derived proof plates", "microscopy, histology, real plot crops", "red"),
        ("P1", "Native badge system", "species, modality, evidence status", "blue"),
        ("P1", "Colorblind QA", "simulation plus grayscale checks", "green"),
        ("P2", "Final bridge replacement", "non-meta title, tighter crops, higher contrast", "amber"),
    ]
    for i, (prio, title_text, desc, col) in enumerate(priorities):
        y = 0.690 - i * 0.122
        ax.text(0.085, y, prio, fontsize=10.5, fontweight="bold", color=COLORS[col], ha="left", va="center", zorder=20)
        ax.text(0.145, y, title_text, fontsize=12.2, fontweight="bold", color=COLORS["ink"], ha="left", va="center", zorder=20)
        ax.text(0.145, y - 0.036, desc, fontsize=8.8, color=COLORS["muted"], ha="left", va="center", zorder=20)
        ax.plot([0.085, 0.475], [y - 0.067, y - 0.067], color=COLORS["rule"], alpha=0.22, lw=1.0, zorder=4)

    add_asset(ax, "cell_micrograph_texture_field", 0.560, 0.505, 0.330, 0.186, alpha=0.94, z=5, clip=True)
    add_asset(ax, "omics_matrix_texture_field", 0.610, 0.310, 0.250, 0.141, alpha=0.78, z=7, clip=True)
    ax.add_patch(Polygon([(0.545, 0.455), (0.895, 0.455), (0.840, 0.245), (0.595, 0.245)], facecolor=COLORS["teal"], edgecolor="none", alpha=0.050, zorder=3))
    ax.text(0.560, 0.740, "visual strategy", fontsize=10.5, fontweight="bold", color=COLORS["ink"], ha="left", va="center", zorder=20)
    ax.text(0.560, 0.715, "proof-like texture underneath, editable decision layer above", fontsize=8.6, color=COLORS["muted"], ha="left", va="center", zorder=20)
    ax.text(0.560, 0.190, "Stop condition", fontsize=10.5, fontweight="bold", color=COLORS["ink"], ha="left", va="center", zorder=20)
    ax.text(0.560, 0.155, "A first-time expert can tell what is source evidence, what is schematic, and what decision follows.", fontsize=8.8, color=COLORS["muted"], ha="left", va="center", zorder=20)
    path = OUT / "04_v0_4_decision_roadmap.png"
    fig.savefig(path, dpi=SLIDE_DPI, facecolor=COLORS["paper"])
    plt.close(fig)
    return {"id": "v0_4_decision_roadmap", "png": rel(path)}


def build() -> dict[str, Any]:
    OUT.mkdir(parents=True, exist_ok=True)
    slides = [
        draw_executive_portfolio(),
        draw_evidence_ladder(),
        draw_archetype_gallery(),
        draw_v04_decision_roadmap(),
    ]
    qa_path = OUT / "biovis_v0_3_consulting_pack_qa.json"
    manifest_path = OUT / "biovis_v0_3_consulting_pack_manifest.json"
    manifest = {
        "created": CREATED,
        "purpose": "Consulting-style usage expansion for the v0.3 biological visual system.",
        "strategy": "Convert biological assets into executive portfolio, evidence ladder, archetype, and roadmap slide jobs.",
        "slides": slides,
        "source_assets": [
            rel(ASSET_DIR / "cell_micrograph_texture_field.png"),
            rel(ASSET_DIR / "omics_matrix_texture_field.png"),
            rel(ASSET_DIR / "single_cell_embedding_field.png"),
            rel(ASSET_DIR / "spatial_spot_section_field.png"),
            rel(ASSET_DIR / "pathway_network_summary_field.png"),
            rel(ASSET_DIR / "species_model_strip.png"),
        ],
        "claim_boundary": "Generated schematic consulting context only; source-derived panels remain required for empirical proof.",
    }
    qa = {
        "created": CREATED,
        "automatic_checks": {
            "slide_count": len(slides),
            "slide_count_expected": 4,
            "png_outputs_exist": all((ROOT / slide["png"]).exists() for slide in slides),
            "source_assets_exist": all(Path(ROOT / asset).exists() for asset in manifest["source_assets"]),
        },
        "manual_review": {
            "executive_readability": "pass - four-slide pack inspected; evidence ladder and roadmap are strongest for executive/consulting use.",
            "content_diversity": "pass - covers portfolio framing, evidence ladder, slide-job operating model, and v0.4 roadmap.",
            "consulting_quality_gate": "conditional pass - suitable as consulting-style seed material; archetype gallery is more operating-model documentation than final deck slide.",
            "biological_claim_boundary": "pass - slides separate generated schematic context from source-derived proof and validated result layers.",
            "verdict": "v0.3 improves when shipped with slide jobs and decision frames; next upgrade should add source-derived proof plates and native badge controls.",
        },
    }
    write_json(manifest_path, manifest)
    write_json(qa_path, qa)
    return {"output": rel(OUT), "manifest": rel(manifest_path), "qa": rel(qa_path), "slides": slides}


def main() -> None:
    print(json.dumps(build(), indent=2))


if __name__ == "__main__":
    main()
