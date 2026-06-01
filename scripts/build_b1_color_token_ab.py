#!/usr/bin/env python3
"""Render B1 color-token A/B candidates without changing global brand tokens."""

from __future__ import annotations

import copy
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
from PIL import Image, ImageDraw, ImageFont

from build_b1_brand_stress_test_scene import (
    CANVAS,
    CREATED,
    SLIDE_DPI,
    build_overlay,
    color,
    draw_scene,
    load_tokens,
    render_overlay,
    write_json,
)
from build_premium_bridge_color_qa import (
    contrast_ratio,
    hex_to_rgb01,
    simulate_image,
    simulate_rgb,
    rgb01_to_hex,
    rgb_distance,
)


OUT_ROOT = ROOT / "output" / "premium_bridge_color_ab"
SLIDE_ID = "b1_evaluation_layer"

VARIANTS = [
    {
        "variant_id": "baseline_current",
        "label": "baseline",
        "description": "current visual-identity v0.1 color tokens",
        "overrides": {},
    },
    {
        "variant_id": "a_amber_dark",
        "label": "A",
        "description": "darken amber only for better small-text contrast",
        "overrides": {"label_amber": "#A36F13"},
    },
    {
        "variant_id": "b_amber_dark_purple_calm",
        "label": "B",
        "description": "darken amber and calm purple for a more editorial task accent",
        "overrides": {"label_amber": "#A36F13", "model_purple": "#5D4BA8"},
    },
]


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def load_font(size: int, *, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    font_paths = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/Library/Fonts/Arial Bold.ttf" if bold else "/Library/Fonts/Arial.ttf",
    ]
    for path in font_paths:
        if Path(path).exists():
            return ImageFont.truetype(path, size=size)
    return ImageFont.load_default()


def y_from_slide(value: float) -> float:
    return 1.0 - value


def make_figure() -> tuple[plt.Figure, plt.Axes]:
    fig = plt.figure(figsize=(16, 9), dpi=SLIDE_DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    return fig, ax


def tokens_for_variant(base_tokens: dict[str, Any], variant: dict[str, Any]) -> dict[str, Any]:
    tokens = copy.deepcopy(base_tokens)
    for key, value in variant["overrides"].items():
        tokens["colors"][key] = value
    return tokens


def render_variant(base_tokens: dict[str, Any], variant: dict[str, Any]) -> dict[str, Any]:
    variant_id = variant["variant_id"]
    variant_dir = OUT_ROOT / variant_id
    variant_dir.mkdir(parents=True, exist_ok=True)
    tokens = tokens_for_variant(base_tokens, variant)
    overlay = build_overlay()
    outputs = {
        "rendered_preview_png": rel(variant_dir / "rendered_preview.png"),
        "rendered_preview_pdf": rel(variant_dir / "rendered_preview.pdf"),
        "grayscale_png": rel(variant_dir / "rendered_preview_grayscale.png"),
        "tritanopia_png": rel(variant_dir / "rendered_preview_tritanopia.png"),
        "manifest": rel(variant_dir / "manifest.json"),
        "qa": rel(variant_dir / "qa.json"),
    }
    preview_png = ROOT / outputs["rendered_preview_png"]
    preview_pdf = ROOT / outputs["rendered_preview_pdf"]

    fig, ax = make_figure()
    draw_scene(ax, tokens)
    render_overlay(ax, tokens, overlay)
    fig.savefig(preview_png, dpi=SLIDE_DPI, facecolor=color(tokens, "canvas_warm"))
    fig.savefig(preview_pdf, facecolor=color(tokens, "canvas_warm"))
    plt.close(fig)

    image = Image.open(preview_png).convert("RGB")
    simulate_image(image, "grayscale").save(ROOT / outputs["grayscale_png"])
    simulate_image(image, "tritanopia").save(ROOT / outputs["tritanopia_png"])

    metrics = token_metrics(tokens)
    qa = {
        "created": CREATED,
        "slide_id": SLIDE_ID,
        "variant_id": variant_id,
        "description": variant["description"],
        "overrides": variant["overrides"],
        "outputs": outputs,
        "image_dimensions": list(image.size),
        "expected_dimensions_match": list(image.size) == [CANVAS["width_px"], CANVAS["height_px"]],
        "token_metrics": metrics,
        "manual_review": {
            "original_readability": "pending",
            "grayscale_readability": "pending",
            "tritanopia_readability": "pending",
            "premium_verdict": "pending",
        },
    }
    manifest = {
        "created": CREATED,
        "slide_id": SLIDE_ID,
        "variant_id": variant_id,
        "description": variant["description"],
        "base_token_version": "0.1",
        "overrides": variant["overrides"],
        "generator": "scripts/build_b1_color_token_ab.py",
        "outputs": outputs,
        "qa": {"rendered": True, "manual_review_required": True},
    }
    write_json(ROOT / outputs["qa"], qa)
    write_json(ROOT / outputs["manifest"], manifest)
    return {"variant": variant, "outputs": outputs, "qa": qa}


def token_metrics(tokens: dict[str, Any]) -> dict[str, Any]:
    colors = tokens["colors"]
    paper = hex_to_rgb01(colors["paper"])
    names = ["source_blue", "bio_green", "assay_teal", "label_amber", "model_purple", "muted"]
    records: dict[str, Any] = {}
    for name in names:
        rgb = hex_to_rgb01(colors[name])
        records[name] = {
            "hex": colors[name],
            "contrast_vs_paper": round(contrast_ratio(rgb, paper), 2),
            "tritanopia_hex": rgb01_to_hex(simulate_rgb(rgb, "tritanopia")),
        }
    pair_names = [("source_blue", "bio_green"), ("source_blue", "assay_teal"), ("bio_green", "assay_teal"), ("model_purple", "muted")]
    pairwise = {}
    for left, right in pair_names:
        left_rgb = hex_to_rgb01(colors[left])
        right_rgb = hex_to_rgb01(colors[right])
        pairwise[f"{left}__{right}"] = {
            "original_distance": round(rgb_distance(left_rgb, right_rgb), 1),
            "tritanopia_distance": round(rgb_distance(simulate_rgb(left_rgb, "tritanopia"), simulate_rgb(right_rgb, "tritanopia")), 1),
        }
    return {"records": records, "pairwise": pairwise}


def panel_label(width: int, title: str, subtitle: str) -> Image.Image:
    bar = Image.new("RGB", (width, 58), "#FFFFFF")
    draw = ImageDraw.Draw(bar)
    draw.text((16, 9), title, fill="#17212B", font=load_font(21, bold=True))
    draw.text((16, 34), subtitle, fill="#5D6978", font=load_font(14))
    return bar


def build_variant_sheet(rendered: list[dict[str, Any]], mode: str, output: Path, title: str) -> None:
    panels = []
    for item in rendered:
        source = ROOT / item["outputs"][mode]
        image = Image.open(source).convert("RGB")
        scale = 760 / image.width
        thumb = image.resize((760, int(image.height * scale)), Image.Resampling.LANCZOS)
        variant = item["variant"]
        subtitle = variant["description"]
        panel = Image.new("RGB", (thumb.width, thumb.height + 58), "#FFFFFF")
        panel.paste(panel_label(thumb.width, f"{variant['label']} | {variant['variant_id']}", subtitle), (0, 0))
        panel.paste(thumb, (0, 58))
        panels.append(panel)
    gap = 22
    header_h = 92
    sheet = Image.new("RGB", (len(panels) * panels[0].width + (len(panels) - 1) * gap, header_h + panels[0].height), "#FFFFFF")
    draw = ImageDraw.Draw(sheet)
    draw.text((0, 0), title, fill="#17212B", font=load_font(30, bold=True))
    draw.text((0, 48), "B1 color-token candidates. Global visual identity config is not changed.", fill="#5D6978", font=load_font(17))
    for idx, panel in enumerate(panels):
        sheet.paste(panel, (idx * (panel.width + gap), header_h))
    sheet.save(output)


def build_palette_delta_sheet(rendered: list[dict[str, Any]], output: Path) -> None:
    fig = plt.figure(figsize=(12.8, 5.8), dpi=220)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.text(0.035, 0.915, "B1 color-token A/B palette deltas", fontsize=15.5, fontweight="bold", color="#17212B", transform=ax.transAxes)
    ax.text(0.035, 0.862, "Compare amber and purple against the current token set before final deck polish.", fontsize=8.2, color="#5D6978", transform=ax.transAxes)
    row_y = [0.685, 0.470, 0.255]
    for y, item in zip(row_y, rendered):
        variant = item["variant"]
        qa = item["qa"]
        metrics = qa["token_metrics"]["records"]
        ax.text(0.035, y + 0.055, f"{variant['label']} | {variant['variant_id']}", fontsize=9.2, fontweight="bold", color="#17212B", transform=ax.transAxes)
        ax.text(0.035, y + 0.015, variant["description"], fontsize=7.4, color="#5D6978", transform=ax.transAxes)
        for idx, name in enumerate(["label_amber", "model_purple", "source_blue", "bio_green", "assay_teal"]):
            x = 0.335 + idx * 0.118
            hex_value = metrics[name]["hex"]
            ax.add_patch(plt.Rectangle((x, y), 0.082, 0.072, facecolor=hex_value, edgecolor="#FFFFFF", linewidth=1.0, transform=ax.transAxes))
            ax.text(x, y - 0.027, name.replace("_", " "), fontsize=6.3, color="#5D6978", ha="left", transform=ax.transAxes)
            ax.text(x, y - 0.052, f"cr {metrics[name]['contrast_vs_paper']}", fontsize=6.1, color="#5D6978", ha="left", transform=ax.transAxes)
    fig.savefig(output, dpi=220, facecolor="white")
    plt.close(fig)


def build() -> dict[str, Any]:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    base_tokens = load_tokens()
    rendered = [render_variant(base_tokens, variant) for variant in VARIANTS]

    sheets = {
        "original_comparison": rel(OUT_ROOT / "b1_color_token_ab_original.png"),
        "grayscale_comparison": rel(OUT_ROOT / "b1_color_token_ab_grayscale.png"),
        "tritanopia_comparison": rel(OUT_ROOT / "b1_color_token_ab_tritanopia.png"),
        "palette_delta": rel(OUT_ROOT / "b1_color_token_ab_palette_delta.png"),
        "manifest": rel(OUT_ROOT / "manifest.json"),
        "qa": rel(OUT_ROOT / "qa.json"),
    }
    build_variant_sheet(rendered, "rendered_preview_png", ROOT / sheets["original_comparison"], "B1 color-token A/B - original")
    build_variant_sheet(rendered, "grayscale_png", ROOT / sheets["grayscale_comparison"], "B1 color-token A/B - grayscale")
    build_variant_sheet(rendered, "tritanopia_png", ROOT / sheets["tritanopia_comparison"], "B1 color-token A/B - tritanopia")
    build_palette_delta_sheet(rendered, ROOT / sheets["palette_delta"])

    qa = {
        "created": CREATED,
        "slide_id": SLIDE_ID,
        "variants": [
            {
                "variant_id": item["variant"]["variant_id"],
                "description": item["variant"]["description"],
                "overrides": item["variant"]["overrides"],
                "outputs": item["outputs"],
                "metrics": item["qa"]["token_metrics"],
            }
            for item in rendered
        ],
        "comparison_sheets": sheets,
        "manual_review": {
            "original_preference": "pending",
            "grayscale_preference": "pending",
            "tritanopia_preference": "pending",
            "premium_preference": "pending",
            "recommended_token_change": "pending",
        },
    }
    manifest = {
        "created": CREATED,
        "slide_id": SLIDE_ID,
        "purpose": "A/B comparison for final-polish color token candidates",
        "variants": VARIANTS,
        "comparison_sheets": sheets,
        "generator": "scripts/build_b1_color_token_ab.py",
    }
    write_json(ROOT / sheets["qa"], qa)
    write_json(ROOT / sheets["manifest"], manifest)
    return {"comparison_sheets": sheets, "qa": qa}


def main() -> None:
    print(json.dumps(build(), indent=2))


if __name__ == "__main__":
    main()
