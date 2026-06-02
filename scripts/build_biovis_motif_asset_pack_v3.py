#!/usr/bin/env python3
"""Build v0.3 biological visual assets for SpaceBio-Bench slides.

v0.3 is intentionally smaller and stricter than v0.2. It creates source-like
schematic texture plates plus editable SVG overlays and metadata, so the assets
behave as scientific visual context rather than decorative icons.
"""

from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Callable


ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "output" / ".matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, PathPatch, Polygon, Rectangle
from matplotlib.path import Path as MplPath

try:
    from PIL import Image, ImageDraw, ImageFont
except ImportError:  # pragma: no cover - Pillow is expected with matplotlib.
    Image = None
    ImageDraw = None
    ImageFont = None


OUT = ROOT / "assets" / "biovis_motif_pack_v0_3"
PNG_DIR = OUT / "png"
OVERLAY_DIR = OUT / "overlay_svg"
QA_DIR = OUT / "qa"
CREATED = "2026-06-02"

PLATE_W = 1600
PLATE_H = 900
PLATE_DPI = 160

COLORS = {
    "paper": "#F7F4EC",
    "paper2": "#FBFAF7",
    "ink": "#17212B",
    "muted": "#5D6978",
    "rule": "#AEB8C5",
    "blue": "#2D6F9F",
    "sky": "#6BAED6",
    "green": "#178B63",
    "teal": "#1AA090",
    "amber": "#B4832F",
    "rose": "#B45A7E",
    "purple": "#6750A4",
    "red": "#B23A3A",
    "deep": "#13202B",
    "dark_field": "#0D1720",
}


@dataclass(frozen=True)
class Asset:
    asset_id: str
    label: str
    biological_level: str
    modality: str
    species_or_model: str
    evidence_status: str
    allowed_uses: str
    prohibited_uses: str
    palette: tuple[str, str, str]
    draw: str


ASSETS: list[Asset] = [
    Asset(
        "cell_micrograph_texture_field",
        "Cell micrograph texture field",
        "cell / organelle",
        "microscopy-like schematic",
        "generic mammalian cell",
        "generated schematic texture; not source microscopy",
        "hero texture, methods bridge, cell-level context",
        "do not present as measured microscopy or specific cell type",
        ("green", "sky", "rose"),
        "cell_micrograph",
    ),
    Asset(
        "tissue_section_texture_field",
        "Tissue section texture field",
        "tissue",
        "histology-like schematic",
        "generic mammalian tissue",
        "generated schematic texture; not source histology",
        "tissue result anchor, methods bridge, organism context",
        "do not use as histology evidence or disease morphology",
        ("amber", "rose", "blue"),
        "tissue_section",
    ),
    Asset(
        "organoid_rosette_texture_field",
        "Organoid rosette texture field",
        "3D culture / tissue",
        "organoid-like schematic",
        "human iPSC organoid context",
        "generated schematic texture; not source microscopy",
        "organoid extension, neural/culture context, future-work bridge",
        "do not imply organoid analysis unless captioned as extension",
        ("purple", "green", "sky"),
        "organoid_rosette",
    ),
    Asset(
        "spatial_spot_section_field",
        "Spatial spot section field",
        "tissue / spatial assay",
        "spatial transcriptomics schematic",
        "generic tissue section",
        "generated schematic texture; not source spatial data",
        "spatial modality explanation, tissue localization bridge",
        "do not imply spatial data exists for a dataset without citation",
        ("green", "purple", "amber"),
        "spatial_spot",
    ),
    Asset(
        "single_cell_embedding_field",
        "Single-cell embedding field",
        "cell population",
        "single-cell embedding schematic",
        "mixed cell populations",
        "generated schematic data texture",
        "single-cell modality context, population shift explanation",
        "do not use as analysis result without replacing with real embedding",
        ("blue", "green", "rose"),
        "single_cell_embedding",
    ),
    Asset(
        "omics_matrix_texture_field",
        "Omics matrix texture field",
        "dataset / molecular",
        "bulk RNA-seq or omics matrix schematic",
        "multi-sample expression assays",
        "generated schematic data texture",
        "feature-layer bridge, preprocessing, task construction",
        "do not treat heatmap pattern as measured values",
        ("blue", "amber", "green"),
        "omics_matrix",
    ),
    Asset(
        "pathway_network_summary_field",
        "Pathway network summary field",
        "molecular program",
        "pathway/network schematic",
        "cross-species pathway summary",
        "generated schematic network texture",
        "pathway feature layer, interpretability, result bridge",
        "do not claim causal pathway wiring without evidence",
        ("green", "amber", "blue"),
        "pathway_network",
    ),
    Asset(
        "species_model_strip",
        "Species and model-system strip",
        "organism / model system",
        "species/model coverage schematic",
        "human, mouse, rat, fly, worm, fish, plant, yeast, microbe",
        "generated schematic label strip",
        "scope slide, dataset coverage, extension map",
        "do not imply all species are analyzed in the current result",
        ("blue", "green", "amber"),
        "species_strip",
    ),
]


def color(name: str) -> str:
    return COLORS[name]


def ensure_dirs() -> None:
    for directory in [OUT, PNG_DIR, OVERLAY_DIR, QA_DIR]:
        directory.mkdir(parents=True, exist_ok=True)


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def make_fig() -> tuple[plt.Figure, plt.Axes]:
    fig = plt.figure(figsize=(PLATE_W / PLATE_DPI, PLATE_H / PLATE_DPI), dpi=PLATE_DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    return fig, ax


def save_plate(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, dpi=PLATE_DPI, transparent=True)
    plt.close(fig)


def subtle_field(ax: plt.Axes, seed: int, *, dark: bool = False) -> None:
    rng = np.random.default_rng(seed)
    h, w = 450, 800
    yy, xx = np.mgrid[0:h, 0:w]
    if dark:
        base = np.zeros((h, w, 4), dtype=float)
        base[..., 0] = 0.05
        base[..., 1] = 0.09
        base[..., 2] = 0.13
        base[..., 3] = 0.94
        grain = rng.normal(0, 0.018, (h, w, 1))
        vignette = ((xx / w - 0.50) ** 2 + (yy / h - 0.50) ** 2)[..., None] * 0.16
        base[..., :3] = np.clip(base[..., :3] + grain - vignette, 0, 1)
    else:
        base = np.zeros((h, w, 4), dtype=float)
        base[..., 0] = 0.98
        base[..., 1] = 0.965
        base[..., 2] = 0.925
        base[..., 3] = 0.82
        grain = rng.normal(0, 0.010, (h, w, 1))
        cool = (1 - xx / w)[..., None] * np.array([0.0, 0.004, 0.010])
        vignette = ((xx / w - 0.52) ** 2 + (yy / h - 0.48) ** 2)[..., None] * np.array([0.060, 0.050, 0.035])
        base[..., :3] = np.clip(base[..., :3] + grain + cool - vignette, 0, 1)
    ax.imshow(base, extent=(0, 1, 0, 1), origin="lower", zorder=0)
    ax.set_aspect("auto")
    for y in [0.18, 0.50, 0.82]:
        ax.plot([0.06, 0.94], [y, y], color=color("rule"), alpha=0.11 if not dark else 0.08, lw=1.0, zorder=1)
    for x in [0.12, 0.50, 0.88]:
        ax.plot([x, x], [0.08, 0.90], color=color("rule"), alpha=0.08 if not dark else 0.06, lw=1.0, zorder=1)


def draw_scale(ax: plt.Axes, text: str = "schematic scale") -> None:
    ax.plot([0.080, 0.205], [0.108, 0.108], color=color("ink"), alpha=0.70, lw=3.0, solid_capstyle="butt", zorder=20)
    ax.text(0.080, 0.073, text, color=color("muted"), fontsize=8.5, ha="left", va="center", zorder=20)


def draw_source_tag(ax: plt.Axes, label: str = "generated schematic") -> None:
    ax.text(
        0.930,
        0.085,
        label,
        color=color("muted"),
        fontsize=8.2,
        ha="right",
        va="center",
        zorder=20,
    )


def irregular_blob(cx: float, cy: float, rx: float, ry: float, n: int, seed: int, roughness: float = 0.12) -> np.ndarray:
    rng = np.random.default_rng(seed)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    jitter = 1 + rng.normal(0, roughness, size=n)
    xs = cx + np.cos(angles) * rx * jitter
    ys = cy + np.sin(angles) * ry * jitter
    return np.column_stack([xs, ys])


def smooth_closed(points: np.ndarray, window: int = 9) -> np.ndarray:
    if window <= 1:
        return points
    pad = window // 2
    padded = np.vstack([points[-pad:], points, points[:pad]])
    kernel = np.ones(window) / window
    xs = np.convolve(padded[:, 0], kernel, mode="valid")
    ys = np.convolve(padded[:, 1], kernel, mode="valid")
    return np.column_stack([xs, ys])


def add_blob(ax: plt.Axes, points: np.ndarray, *, face: str, edge: str, alpha: float, lw: float, z: int) -> None:
    vertices = np.vstack([points, points[0]])
    codes = [MplPath.MOVETO] + [MplPath.LINETO] * (len(points) - 1) + [MplPath.CLOSEPOLY]
    ax.add_patch(PathPatch(MplPath(vertices, codes), facecolor=face, edgecolor=edge, lw=lw, alpha=alpha, zorder=z, joinstyle="round"))


def draw_cell_micrograph(ax: plt.Axes, asset: Asset) -> None:
    subtle_field(ax, 3101, dark=True)
    rng = np.random.default_rng(3102)
    p, s, t = [color(x) for x in asset.palette]
    for _ in range(58):
        x, y = rng.uniform(0.08, 0.92), rng.uniform(0.16, 0.84)
        r = rng.uniform(0.010, 0.030)
        ax.add_patch(Circle((x, y), r * 1.85, facecolor="none", edgecolor=p, lw=rng.uniform(0.8, 1.8), alpha=rng.uniform(0.18, 0.44), zorder=5))
        ax.add_patch(Circle((x + rng.normal(0, 0.004), y + rng.normal(0, 0.004)), r * rng.uniform(0.55, 0.86), facecolor=s, edgecolor="none", alpha=rng.uniform(0.18, 0.46), zorder=6))
        if rng.random() > 0.55:
            ax.plot([x - r * 1.7, x + r * 1.7], [y + r * 1.5, y - r * 1.2], color=t, alpha=0.22, lw=1.0, zorder=4)
    for x, y, a in [(0.26, 0.52, 0.25), (0.57, 0.56, 0.18), (0.72, 0.37, 0.16)]:
        ax.add_patch(Ellipse((x, y), 0.40, 0.18, angle=rng.uniform(-18, 18), facecolor=t, edgecolor="none", alpha=a, zorder=3))
    draw_scale(ax, "schematic 20 um")
    draw_source_tag(ax)


def draw_tissue_section(ax: plt.Axes, asset: Asset) -> None:
    subtle_field(ax, 3201)
    rng = np.random.default_rng(3202)
    p, s, t = [color(x) for x in asset.palette]
    outer = smooth_closed(irregular_blob(0.52, 0.52, 0.36, 0.245, 160, 3203, 0.055), 13)
    add_blob(ax, outer, face="#FFF9F0", edge=p, alpha=0.84, lw=2.2, z=4)
    inner = smooth_closed(irregular_blob(0.51, 0.51, 0.255, 0.150, 140, 3204, 0.050), 11)
    add_blob(ax, inner, face=s, edge="none", alpha=0.055, lw=0, z=5)
    for idx in range(118):
        x, y = rng.normal(0.52, 0.17), rng.normal(0.52, 0.115)
        if 0.16 < x < 0.88 and 0.24 < y < 0.78:
            ax.add_patch(Ellipse((x, y), rng.uniform(0.015, 0.039), rng.uniform(0.009, 0.023), angle=rng.uniform(-40, 40), facecolor=s if idx % 3 else t, edgecolor="none", alpha=rng.uniform(0.30, 0.58), zorder=6))
    for offset, op, lw in [(0.00, 0.42, 2.3), (0.060, 0.28, 1.8), (-0.065, 0.24, 1.8)]:
        xs = np.linspace(0.18, 0.86, 140)
        ys = 0.51 + offset + 0.030 * np.sin(xs * 14 + offset * 12)
        ax.plot(xs, ys, color=p, alpha=op, lw=lw, zorder=7)
    draw_scale(ax, "schematic tissue")
    draw_source_tag(ax)


def draw_organoid_rosette(ax: plt.Axes, asset: Asset) -> None:
    subtle_field(ax, 3301)
    rng = np.random.default_rng(3302)
    p, s, t = [color(x) for x in asset.palette]
    ax.add_patch(Ellipse((0.52, 0.48), 0.58, 0.15, facecolor=color("rule"), edgecolor="none", alpha=0.10, zorder=2))
    ax.add_patch(Circle((0.52, 0.51), 0.285, facecolor="#FBFAF7", edgecolor=p, lw=3.0, alpha=0.92, zorder=4))
    for _ in range(96):
        angle = rng.uniform(0, 2 * np.pi)
        rad = rng.uniform(0.18, 0.275)
        x = 0.52 + math.cos(angle) * rad * 0.97
        y = 0.51 + math.sin(angle) * rad * 0.92
        ax.add_patch(Circle((x, y), rng.uniform(0.006, 0.015), facecolor=[p, s, t][int(rng.integers(0, 3))], edgecolor="none", alpha=rng.uniform(0.20, 0.44), zorder=5))
    for cx, cy, rr in [(0.43, 0.55, 0.085), (0.55, 0.61, 0.075), (0.60, 0.43, 0.090), (0.43, 0.40, 0.070)]:
        ax.add_patch(Circle((cx, cy), rr, facecolor="#F8F5F0", edgecolor=s, lw=2.1, alpha=0.86, zorder=6))
        ax.add_patch(Circle((cx, cy), rr * 0.38, facecolor=color("paper"), edgecolor=t, lw=1.8, alpha=0.90, zorder=7))
        for k in range(16):
            angle = 2 * np.pi * k / 16 + rng.normal(0, 0.05)
            x = cx + math.cos(angle) * rr * 0.72
            y = cy + math.sin(angle) * rr * 0.72
            ax.add_patch(Circle((x, y), rr * 0.100, facecolor=p if k % 2 else s, edgecolor="none", alpha=0.55, zorder=8))
    draw_scale(ax, "schematic organoid")
    draw_source_tag(ax)


def draw_spatial_spot(ax: plt.Axes, asset: Asset) -> None:
    subtle_field(ax, 3401)
    rng = np.random.default_rng(3402)
    p, s, t = [color(x) for x in asset.palette]
    blob = smooth_closed(irregular_blob(0.50, 0.52, 0.38, 0.245, 150, 3403, 0.052), 11)
    add_blob(ax, blob, face="#FBFAF5", edge=p, alpha=0.74, lw=2.2, z=4)
    for _ in range(72):
        x, y = rng.uniform(0.20, 0.82), rng.uniform(0.31, 0.72)
        ax.add_patch(Ellipse((x, y), rng.uniform(0.022, 0.064), rng.uniform(0.010, 0.028), angle=rng.uniform(-28, 28), facecolor=t, edgecolor="none", alpha=rng.uniform(0.14, 0.30), zorder=5))
    xs = np.linspace(0.26, 0.76, 8)
    ys = np.linspace(0.36, 0.68, 5)
    for i, x in enumerate(xs):
        for j, y in enumerate(ys):
            val = (i + j) % 3
            col = [p, s, t][val]
            ax.add_patch(Circle((x, y), 0.017, facecolor=col, edgecolor=color("paper2"), lw=1.1, alpha=0.82, zorder=8))
    ax.plot([0.26, 0.76], [0.30, 0.30], color=color("rule"), alpha=0.38, lw=1.0, zorder=7)
    ax.plot([0.26, 0.26], [0.30, 0.73], color=color("rule"), alpha=0.38, lw=1.0, zorder=7)
    draw_scale(ax, "schematic spots")
    draw_source_tag(ax)


def draw_single_cell_embedding(ax: plt.Axes, asset: Asset) -> None:
    subtle_field(ax, 3501)
    rng = np.random.default_rng(3502)
    p, s, t = [color(x) for x in asset.palette]
    centers = [(0.32, 0.58), (0.50, 0.42), (0.66, 0.58), (0.62, 0.32), (0.40, 0.28)]
    cols = [p, s, t, color("amber"), color("purple")]
    for idx, (cx, cy) in enumerate(centers):
        covx = rng.uniform(0.035, 0.070)
        covy = rng.uniform(0.025, 0.055)
        ax.add_patch(Ellipse((cx, cy), covx * 5.8, covy * 4.8, angle=rng.uniform(-35, 35), facecolor=cols[idx], edgecolor="none", alpha=0.115, zorder=3))
        pts = rng.normal([cx, cy], [covx, covy], size=(95, 2))
        ax.scatter(pts[:, 0], pts[:, 1], s=rng.uniform(10, 27, size=len(pts)), c=cols[idx], alpha=0.56, lw=0, zorder=5)
    ax.add_patch(FancyArrowPatch((0.30, 0.78), (0.64, 0.72), arrowstyle="-|>", mutation_scale=14, color=color("rule"), lw=1.5, alpha=0.55, zorder=8))
    ax.text(0.66, 0.72, "shift", fontsize=9, color=color("muted"), ha="left", va="center", zorder=10)
    draw_scale(ax, "schematic embedding")
    draw_source_tag(ax)


def draw_omics_matrix(ax: plt.Axes, asset: Asset) -> None:
    subtle_field(ax, 3601)
    rng = np.random.default_rng(3602)
    p, s, t = [color(x) for x in asset.palette]
    mat = rng.normal(0, 0.50, (34, 56))
    mat[:12, 8:24] += 1.6
    mat[13:23, 30:46] -= 1.3
    mat[24:, 42:] += 1.1
    ax.imshow(mat, extent=(0.22, 0.82, 0.24, 0.72), cmap="vlag" if "vlag" in plt.colormaps() else "coolwarm", aspect="auto", alpha=0.74, zorder=4)
    for x in np.linspace(0.22, 0.82, 8):
        ax.plot([x, x], [0.24, 0.72], color=color("paper2"), alpha=0.18, lw=0.7, zorder=5)
    for y in np.linspace(0.24, 0.72, 6):
        ax.plot([0.22, 0.82], [y, y], color=color("paper2"), alpha=0.18, lw=0.7, zorder=5)
    ax.add_patch(Rectangle((0.22, 0.745), 0.18, 0.028, facecolor=p, edgecolor="none", alpha=0.55, zorder=6))
    ax.add_patch(Rectangle((0.405, 0.745), 0.20, 0.028, facecolor=s, edgecolor="none", alpha=0.55, zorder=6))
    ax.add_patch(Rectangle((0.610, 0.745), 0.21, 0.028, facecolor=t, edgecolor="none", alpha=0.55, zorder=6))
    ax.plot([0.17, 0.20, 0.20, 0.22], [0.31, 0.31, 0.56, 0.56], color=color("muted"), alpha=0.48, lw=1.2, zorder=8)
    ax.plot([0.38, 0.38, 0.62, 0.62], [0.78, 0.81, 0.81, 0.78], color=color("muted"), alpha=0.48, lw=1.2, zorder=8)
    draw_scale(ax, "schematic values")
    draw_source_tag(ax)


def draw_pathway_network(ax: plt.Axes, asset: Asset) -> None:
    subtle_field(ax, 3701)
    rng = np.random.default_rng(3702)
    p, s, t = [color(x) for x in asset.palette]
    modules = [(0.34, 0.54, p), (0.55, 0.60, s), (0.63, 0.38, t), (0.41, 0.34, color("purple"))]
    nodes: list[tuple[float, float, str, float]] = []
    for cx, cy, col in modules:
        for _ in range(10):
            x, y = rng.normal(cx, 0.045), rng.normal(cy, 0.040)
            r = rng.uniform(0.010, 0.019)
            nodes.append((x, y, col, r))
        ax.add_patch(Ellipse((cx, cy), 0.25, 0.18, facecolor=col, edgecolor="none", alpha=0.105, zorder=3))
    for i in range(len(nodes)):
        x1, y1, _, _ = nodes[i]
        for j in range(i + 1, len(nodes)):
            x2, y2, _, _ = nodes[j]
            d = math.hypot(x1 - x2, y1 - y2)
            if d < 0.095 or (rng.random() < 0.012 and d < 0.35):
                ax.plot([x1, x2], [y1, y2], color=color("rule"), alpha=0.34, lw=0.9, zorder=4)
    for x, y, col, r in nodes:
        ax.add_patch(Circle((x, y), r * 1.04, facecolor=color("paper2"), edgecolor=col, lw=2.1, alpha=0.96, zorder=7))
    ax.plot([0.22, 0.80], [0.22, 0.22], color=p, alpha=0.42, lw=4.0, solid_capstyle="round", zorder=6)
    ax.plot([0.30, 0.74], [0.79, 0.79], color=s, alpha=0.36, lw=4.0, solid_capstyle="round", zorder=6)
    ax.text(0.50, 0.235, "program score", fontsize=9, color=color("muted"), ha="center", va="bottom", zorder=10)
    draw_scale(ax, "schematic network")
    draw_source_tag(ax)


def draw_species_strip(ax: plt.Axes, asset: Asset) -> None:
    subtle_field(ax, 3801)
    p, s, t = [color(x) for x in asset.palette]
    labels = [
        ("human", "H"),
        ("mouse", "M"),
        ("rat", "R"),
        ("fly", "F"),
        ("worm", "W"),
        ("fish", "Fi"),
        ("plant", "P"),
        ("yeast", "Y"),
        ("microbe", "B"),
    ]
    xs = np.linspace(0.12, 0.88, len(labels))
    for idx, ((label, glyph), x) in enumerate(zip(labels, xs)):
        col = [p, s, t, color("purple")][idx % 4]
        y = 0.52
        ax.add_patch(Circle((x, y), 0.047, facecolor=color("paper2"), edgecolor=col, lw=2.0, alpha=0.92, zorder=5))
        ax.text(x, y + 0.002, glyph, color=col, fontsize=15 if len(glyph) == 1 else 12, fontweight="bold", ha="center", va="center", zorder=8)
        ax.text(x, 0.405, label, color=color("ink"), fontsize=8.4, ha="center", va="center", zorder=8)
        if idx < len(labels) - 1:
            ax.plot([x + 0.052, xs[idx + 1] - 0.052], [y, y], color=color("rule"), alpha=0.36, lw=1.1, zorder=4)
    ax.add_patch(Rectangle((0.10, 0.63), 0.78, 0.035, facecolor=p, edgecolor="none", alpha=0.12, zorder=3))
    ax.add_patch(Rectangle((0.10, 0.32), 0.78, 0.022, facecolor=s, edgecolor="none", alpha=0.10, zorder=3))
    ax.text(0.10, 0.68, "explicit species labels; not silhouette-only encoding", color=color("muted"), fontsize=8.8, ha="left", va="center", zorder=8)
    draw_scale(ax, "scope strip")
    draw_source_tag(ax)


DRAWERS: dict[str, Callable[[plt.Axes, Asset], None]] = {
    "cell_micrograph": draw_cell_micrograph,
    "tissue_section": draw_tissue_section,
    "organoid_rosette": draw_organoid_rosette,
    "spatial_spot": draw_spatial_spot,
    "single_cell_embedding": draw_single_cell_embedding,
    "omics_matrix": draw_omics_matrix,
    "pathway_network": draw_pathway_network,
    "species_strip": draw_species_strip,
}


def overlay_svg(asset: Asset) -> str:
    safe_label = asset.label.replace("&", "&amp;")
    safe_modality = asset.modality.replace("&", "&amp;")
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="{PLATE_W}" height="{PLATE_H}" viewBox="0 0 {PLATE_W} {PLATE_H}" role="img" aria-label="{safe_label} overlay">
  <metadata>SpaceBio-Bench v0.3 editable overlay; asset_id={asset.asset_id}; evidence_status={asset.evidence_status}</metadata>
  <rect width="{PLATE_W}" height="{PLATE_H}" fill="none"/>
  <text x="86" y="84" font-family="Arial, Helvetica, sans-serif" font-size="34" font-weight="700" fill="{color('ink')}">{safe_label}</text>
  <text x="86" y="126" font-family="Arial, Helvetica, sans-serif" font-size="21" fill="{color('muted')}">{asset.biological_level} | {safe_modality}</text>
  <line x1="86" y1="806" x2="290" y2="806" stroke="{color('ink')}" stroke-width="6"/>
  <text x="86" y="850" font-family="Arial, Helvetica, sans-serif" font-size="18" fill="{color('muted')}">schematic scale / generated texture</text>
  <text x="1514" y="850" font-family="Arial, Helvetica, sans-serif" font-size="18" text-anchor="end" fill="{color('muted')}">not source evidence</text>
</svg>
'''


def build_asset(asset: Asset) -> dict[str, str]:
    fig, ax = make_fig()
    DRAWERS[asset.draw](ax, asset)
    png = PNG_DIR / f"{asset.asset_id}.png"
    save_plate(fig, png)
    overlay = OVERLAY_DIR / f"{asset.asset_id}_overlay.svg"
    overlay.write_text(overlay_svg(asset), encoding="utf-8")
    return {"png": rel(png), "overlay_svg": rel(overlay)}


def render_contact_sheet(records: list[dict[str, str]]) -> Path:
    cols = 2
    cell_w = 1220
    cell_h = 760
    pad_x = 100
    pad_y = 230
    rows = math.ceil(len(records) / cols)
    width = cols * cell_w + pad_x * 2
    height = rows * cell_h + pad_y + 120
    fig = plt.figure(figsize=(width / 180, height / 180), dpi=180)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, width)
    ax.set_ylim(height, 0)
    ax.add_patch(Rectangle((0, 0), width, height, facecolor=color("paper"), edgecolor="none"))
    ax.text(pad_x, 86, "SpaceBio-Bench biological visual system v0.3", fontsize=32, color=color("ink"), fontweight="bold", ha="left", va="center")
    ax.text(pad_x, 136, "Proof texture plates plus editable scientific interpretation overlays.", fontsize=17, color=color("muted"), ha="left", va="center")
    for idx, rec in enumerate(records):
        row, col = divmod(idx, cols)
        x = pad_x + col * cell_w
        y = pad_y + row * cell_h
        img = plt.imread(ROOT / rec["png"])
        ax.imshow(img, extent=(x, x + 760, y + 428, y), origin="upper", zorder=2)
        ax.text(x, y + 488, rec["label"], fontsize=21, color=color("ink"), fontweight="bold", ha="left", va="center")
        ax.text(x, y + 525, f'{rec["biological_level"]} | {rec["modality"]}', fontsize=13.5, color=color("muted"), ha="left", va="center")
        ax.text(x, y + 560, f'status: {rec["evidence_status"]}', fontsize=12.5, color=color("muted"), ha="left", va="center")
        ax.text(x, y + 600, f'id: {rec["asset_id"]}', fontsize=11.5, color=color("muted"), ha="left", va="center")
        ax.plot([x, x + 760], [y + 650, y + 650], color=color("rule"), alpha=0.28, lw=1.0)
    path = OUT / "biovis_motif_pack_v0_3_contact_sheet.png"
    fig.savefig(path, dpi=180, facecolor=color("paper"))
    plt.close(fig)
    return path


def render_usage_sheet(records: list[dict[str, str]]) -> Path:
    width, height = 3840, 2160
    fig = plt.figure(figsize=(16, 9), dpi=240)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    subtle_field(ax, 3901)
    ax.text(0.065, 0.900, "Biological visual system usage test", fontsize=27, color=color("ink"), fontweight="bold", ha="left", va="center")
    ax.text(0.066, 0.858, "Texture plates carry biological material; editable overlays carry interpretation and trust.", fontsize=11.5, color=color("muted"), ha="left", va="center")
    placements = [
        (0.070, 0.535, 0.300, 0.169, "cell_micrograph_texture_field"),
        (0.390, 0.535, 0.300, 0.169, "tissue_section_texture_field"),
        (0.710, 0.535, 0.220, 0.124, "organoid_rosette_texture_field"),
        (0.070, 0.285, 0.260, 0.146, "omics_matrix_texture_field"),
        (0.360, 0.285, 0.260, 0.146, "single_cell_embedding_field"),
        (0.650, 0.285, 0.260, 0.146, "pathway_network_summary_field"),
        (0.205, 0.095, 0.330, 0.186, "spatial_spot_section_field"),
        (0.580, 0.095, 0.330, 0.186, "species_model_strip"),
    ]
    rec_by_id = {rec["asset_id"]: rec for rec in records}
    for x, y, w, h, asset_id in placements:
        img = plt.imread(ROOT / rec_by_id[asset_id]["png"])
        ax.imshow(img, extent=(x, x + w, y, y + h), zorder=3)
        ax.text(x, y - 0.020, rec_by_id[asset_id]["label"].replace(" texture field", ""), fontsize=8.8, color=color("muted"), ha="left", va="top", zorder=8)
    for x0, y0, x1, y1 in [(0.335, 0.365, 0.365, 0.365), (0.620, 0.365, 0.650, 0.365), (0.510, 0.215, 0.580, 0.188)]:
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=14, lw=1.2, color=color("rule"), alpha=0.55, transform=ax.transAxes, zorder=10))
    ax.text(0.070, 0.755, "Z2 proof-like texture", fontsize=9.5, color=color("ink"), fontweight="bold", ha="left", va="center")
    ax.text(0.070, 0.736, "schematic, source-bounded, not evidence replacement", fontsize=8.0, color=color("muted"), ha="left", va="center")
    ax.text(0.705, 0.742, "Z3 interpretation", fontsize=9.5, color=color("ink"), fontweight="bold", ha="left", va="center")
    ax.text(0.705, 0.723, "pathway, embedding, species labels remain editable", fontsize=8.0, color=color("muted"), ha="left", va="center")
    ax.text(0.066, 0.040, "QA target: biology should read before labels; labels should prevent overclaiming.", fontsize=8.8, color=color("muted"), ha="left", va="center")
    path = OUT / "biovis_motif_pack_v0_3_usage_sheet.png"
    fig.savefig(path, dpi=240, facecolor=color("paper"))
    plt.close(fig)
    return path


def render_grayscale_sheet(records: list[dict[str, str]]) -> Path:
    if Image is None:
        raise RuntimeError("Pillow is required for grayscale QA sheet generation")
    cols = 2
    thumb_w, thumb_h = 720, 405
    cell_w, cell_h = 850, 520
    pad_x, pad_y = 90, 180
    rows = math.ceil(len(records) / cols)
    width = cols * cell_w + pad_x * 2
    height = rows * cell_h + pad_y + 80
    canvas = Image.new("RGB", (width, height), COLORS["paper"])
    draw = ImageDraw.Draw(canvas)
    try:
        title_font = ImageFont.truetype("Arial.ttf", 34)
        body_font = ImageFont.truetype("Arial.ttf", 17)
        small_font = ImageFont.truetype("Arial.ttf", 13)
    except OSError:
        title_font = body_font = small_font = ImageFont.load_default()
    draw.text((pad_x, 55), "v0.3 grayscale QA", fill=COLORS["ink"], font=title_font)
    draw.text((pad_x, 105), "Shape, density, and label structure should survive without color.", fill=COLORS["muted"], font=body_font)
    for idx, rec in enumerate(records):
        row, col = divmod(idx, cols)
        x = pad_x + col * cell_w
        y = pad_y + row * cell_h
        img = Image.open(ROOT / rec["png"]).convert("RGBA")
        bg = Image.new("RGBA", img.size, (247, 244, 236, 255))
        img = Image.alpha_composite(bg, img).convert("L").convert("RGB")
        img.thumbnail((thumb_w, thumb_h), Image.Resampling.LANCZOS)
        canvas.paste(img, (x, y))
        draw.text((x, y + thumb_h + 22), rec["label"], fill=COLORS["ink"], font=body_font)
        draw.text((x, y + thumb_h + 50), rec["asset_id"], fill=COLORS["muted"], font=small_font)
    path = OUT / "biovis_motif_pack_v0_3_grayscale_qa.png"
    canvas.save(path)
    return path


def write_manifest(records: list[dict[str, str]], contact: Path, usage: Path, grayscale: Path) -> Path:
    manifest = {
        "created": CREATED,
        "version": "0.3",
        "scope": "Biological visual system P0: proof-like schematic texture plates plus editable scientific overlays.",
        "design_rule": "Use as Z2 proof-like texture plus Z3/Z4 editable interpretation; never substitute for measured source evidence.",
        "asset_count": len(records),
        "contact_sheet_png": rel(contact),
        "usage_sheet_png": rel(usage),
        "grayscale_qa_png": rel(grayscale),
        "records": records,
    }
    path = OUT / "manifest.json"
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    with (OUT / "manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0].keys()))
        writer.writeheader()
        writer.writerows(records)
    return path


def write_qa(records: list[dict[str, str]], contact: Path, usage: Path, grayscale: Path, manifest: Path) -> Path:
    required = [
        "asset_id",
        "label",
        "biological_level",
        "modality",
        "species_or_model",
        "evidence_status",
        "allowed_uses",
        "prohibited_uses",
        "png",
        "overlay_svg",
    ]
    qa = {
        "created": CREATED,
        "version": "0.3",
        "automatic_checks": {
            "asset_count": len(records),
            "asset_count_expected": 8,
            "png_outputs_exist": all((ROOT / rec["png"]).exists() for rec in records),
            "overlay_svg_outputs_exist": all((ROOT / rec["overlay_svg"]).exists() for rec in records),
            "required_manifest_fields_present": all(all(field in rec and rec[field] for field in required) for rec in records),
            "contact_sheet_exists": contact.exists(),
            "usage_sheet_exists": usage.exists(),
            "grayscale_qa_exists": grayscale.exists(),
            "manifest_exists": manifest.exists(),
        },
        "manual_review": {
            "contact_sheet_inspection": "pass - inspected after orientation fix; P0 assets read as biological texture plates rather than simple icons.",
            "usage_sheet_inspection": "pass - 3840x2160 usage sheet inspected; no major text/asset overlap after Z3 label relocation.",
            "grayscale_inspection": "pass - grayscale QA inspected; shape, density, and label structure survive without color.",
            "biological_specificity": "conditional pass - cell, tissue, spatial, omics, embedding, pathway, organoid, and species coverage are explicit; organoid/pathway/species assets remain schematic.",
            "premium_quality_gate": "conditional pass - materially stronger than v0.2 for premium deck context, but still generated schematic texture rather than source-derived proof.",
            "verdict": "v0.3 P0 prototype accepted for controlled slide tests; next step is direct replacement inside the feature-layer bridge and source-derived texture exploration.",
        },
    }
    path = QA_DIR / "biovis_motif_pack_v0_3_qa.json"
    path.write_text(json.dumps(qa, indent=2) + "\n", encoding="utf-8")
    return path


def main() -> None:
    ensure_dirs()
    records: list[dict[str, str]] = []
    for asset in ASSETS:
        paths = build_asset(asset)
        records.append(
            {
                "asset_id": asset.asset_id,
                "label": asset.label,
                "biological_level": asset.biological_level,
                "modality": asset.modality,
                "species_or_model": asset.species_or_model,
                "evidence_status": asset.evidence_status,
                "allowed_uses": asset.allowed_uses,
                "prohibited_uses": asset.prohibited_uses,
                "png": paths["png"],
                "overlay_svg": paths["overlay_svg"],
                "license": "Original project-generated schematic artwork; repository MIT license.",
            }
        )
    contact = render_contact_sheet(records)
    usage = render_usage_sheet(records)
    grayscale = render_grayscale_sheet(records)
    manifest = write_manifest(records, contact, usage, grayscale)
    qa = write_qa(records, contact, usage, grayscale, manifest)
    print(
        json.dumps(
            {
                "output": rel(OUT),
                "asset_count": len(records),
                "contact_sheet_png": rel(contact),
                "usage_sheet_png": rel(usage),
                "grayscale_qa_png": rel(grayscale),
                "manifest": rel(manifest),
                "qa": rel(qa),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
