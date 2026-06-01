#!/usr/bin/env python3
"""Build color, contrast, grayscale, and CVD QA assets for premium bridge slides."""

from __future__ import annotations

import json
import os
from itertools import combinations
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "output" / ".matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw, ImageFont


TOKEN_PATH = ROOT / "config" / "visual_identity" / "spacebio_bench_visual_identity_v0_1.json"
OUT_DIR = ROOT / "output" / "premium_bridge_color_qa"
CREATED = "2026-06-01"

SOURCES = {
    "b1_evaluation_layer": ROOT / "output" / "premium_bridge_scenes" / "b1_evaluation_layer" / "rendered_preview.png",
    "b1_b4_family_contact_sheet": ROOT / "output" / "premium_bridge_family_review" / "b1_b4_premium_family_contact_sheet.png",
}

SEMANTIC_COLORS = [
    "source_blue",
    "bio_green",
    "assay_teal",
    "label_amber",
    "test_red",
    "model_purple",
    "muted",
    "rule",
]

CVD_MATRICES = {
    "protanopia": np.array(
        [
            [0.152286, 1.052583, -0.204868],
            [0.114503, 0.786281, 0.099216],
            [-0.003882, -0.048116, 1.051998],
        ],
        dtype=float,
    ),
    "deuteranopia": np.array(
        [
            [0.367322, 0.860646, -0.227968],
            [0.280085, 0.672501, 0.047413],
            [-0.011820, 0.042940, 0.968881],
        ],
        dtype=float,
    ),
    "tritanopia": np.array(
        [
            [1.255528, -0.076749, -0.178779],
            [-0.078411, 0.930809, 0.147602],
            [0.004733, 0.691367, 0.303900],
        ],
        dtype=float,
    ),
}


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def write_json(path: Path, data: Any) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def load_tokens() -> dict[str, Any]:
    return json.loads(TOKEN_PATH.read_text(encoding="utf-8"))


def hex_to_rgb01(value: str) -> np.ndarray:
    value = value.strip().lstrip("#")
    return np.array([int(value[idx : idx + 2], 16) / 255.0 for idx in (0, 2, 4)], dtype=float)


def rgb01_to_hex(rgb: np.ndarray) -> str:
    rgb8 = np.clip(np.rint(rgb * 255.0), 0, 255).astype(int)
    return "#" + "".join(f"{channel:02X}" for channel in rgb8)


def srgb_to_linear(rgb: np.ndarray) -> np.ndarray:
    return np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)


def linear_to_srgb(rgb: np.ndarray) -> np.ndarray:
    rgb = np.clip(rgb, 0.0, 1.0)
    return np.where(rgb <= 0.0031308, rgb * 12.92, 1.055 * np.power(rgb, 1 / 2.4) - 0.055)


def relative_luminance(rgb: np.ndarray) -> float:
    linear = srgb_to_linear(rgb)
    return float(linear[0] * 0.2126 + linear[1] * 0.7152 + linear[2] * 0.0722)


def contrast_ratio(foreground: np.ndarray, background: np.ndarray) -> float:
    lum_a = relative_luminance(foreground)
    lum_b = relative_luminance(background)
    lighter = max(lum_a, lum_b)
    darker = min(lum_a, lum_b)
    return float((lighter + 0.05) / (darker + 0.05))


def simulate_rgb(rgb: np.ndarray, mode: str) -> np.ndarray:
    if mode == "grayscale":
        lum = relative_luminance(rgb)
        return np.array([lum, lum, lum], dtype=float)
    matrix = CVD_MATRICES[mode]
    linear = srgb_to_linear(rgb)
    simulated = linear @ matrix.T
    return linear_to_srgb(simulated)


def simulate_image(image: Image.Image, mode: str) -> Image.Image:
    arr = np.asarray(image.convert("RGB"), dtype=float) / 255.0
    if mode == "grayscale":
        linear = srgb_to_linear(arr)
        lum = linear[..., 0] * 0.2126 + linear[..., 1] * 0.7152 + linear[..., 2] * 0.0722
        out = linear_to_srgb(np.repeat(lum[..., None], 3, axis=2))
    else:
        linear = srgb_to_linear(arr)
        out = linear_to_srgb(np.einsum("...c,dc->...d", linear, CVD_MATRICES[mode]))
    return Image.fromarray(np.clip(np.rint(out * 255.0), 0, 255).astype(np.uint8), mode="RGB")


def rgb_distance(rgb_a: np.ndarray, rgb_b: np.ndarray) -> float:
    return float(np.linalg.norm((rgb_a - rgb_b) * 255.0))


def load_font(size: int, *, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    font_paths = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/Library/Fonts/Arial Bold.ttf" if bold else "/Library/Fonts/Arial.ttf",
    ]
    for path in font_paths:
        if Path(path).exists():
            return ImageFont.truetype(path, size=size)
    return ImageFont.load_default()


def label_bar(width: int, label: str) -> Image.Image:
    bar = Image.new("RGB", (width, 46), "#FFFFFF")
    draw = ImageDraw.Draw(bar)
    draw.text((18, 13), label, fill="#17212B", font=load_font(20, bold=True))
    return bar


def build_mode_sheet(source: Path, output: Path, title: str, *, thumb_width: int = 900) -> dict[str, str]:
    original = Image.open(source).convert("RGB")
    variants = {
        "original": original,
        "grayscale": simulate_image(original, "grayscale"),
        "protanopia": simulate_image(original, "protanopia"),
        "deuteranopia": simulate_image(original, "deuteranopia"),
        "tritanopia": simulate_image(original, "tritanopia"),
    }
    variant_paths: dict[str, str] = {}
    resized_panels = []
    for mode, image in variants.items():
        variant_path = output.with_name(f"{output.stem}_{mode}.png")
        image.save(variant_path)
        variant_paths[mode] = rel(variant_path)
        scale = thumb_width / image.width
        thumb = image.resize((thumb_width, int(image.height * scale)), Image.Resampling.LANCZOS)
        panel = Image.new("RGB", (thumb.width, thumb.height + 46), "#FFFFFF")
        panel.paste(label_bar(thumb.width, mode), (0, 0))
        panel.paste(thumb, (0, 46))
        resized_panels.append(panel)

    gap = 24
    columns = 2
    rows = 3
    panel_w = resized_panels[0].width
    panel_h = resized_panels[0].height
    header_h = 82
    sheet = Image.new("RGB", (columns * panel_w + gap, header_h + rows * panel_h + (rows - 1) * gap), "#FFFFFF")
    draw = ImageDraw.Draw(sheet)
    draw.text((0, 0), title, fill="#17212B", font=load_font(28, bold=True))
    draw.text((0, 44), "Color QA sheet: original, grayscale, and color-vision simulations.", fill="#5D6978", font=load_font(18))
    for idx, panel in enumerate(resized_panels):
        x = (idx % columns) * (panel_w + gap)
        y = header_h + (idx // columns) * (panel_h + gap)
        sheet.paste(panel, (x, y))
    output.parent.mkdir(parents=True, exist_ok=True)
    sheet.save(output)
    variant_paths["sheet"] = rel(output)
    return variant_paths


def palette_metrics(tokens: dict[str, Any]) -> dict[str, Any]:
    colors = tokens["colors"]
    backgrounds = {
        "paper": hex_to_rgb01(colors["paper"]),
        "canvas_warm": hex_to_rgb01(colors["canvas_warm"]),
        "canvas_cool": hex_to_rgb01(colors["canvas_cool"]),
    }
    records: list[dict[str, Any]] = []
    for name in SEMANTIC_COLORS:
        rgb = hex_to_rgb01(colors[name])
        records.append(
            {
                "name": name,
                "hex": colors[name],
                "luminance": round(relative_luminance(rgb), 4),
                "contrast_vs_paper": round(contrast_ratio(rgb, backgrounds["paper"]), 2),
                "contrast_vs_canvas_warm": round(contrast_ratio(rgb, backgrounds["canvas_warm"]), 2),
                "contrast_vs_canvas_cool": round(contrast_ratio(rgb, backgrounds["canvas_cool"]), 2),
                "simulated_hex": {mode: rgb01_to_hex(simulate_rgb(rgb, mode)) for mode in ["grayscale", "protanopia", "deuteranopia", "tritanopia"]},
            }
        )

    pairwise: list[dict[str, Any]] = []
    for left, right in combinations(SEMANTIC_COLORS, 2):
        rgb_left = hex_to_rgb01(colors[left])
        rgb_right = hex_to_rgb01(colors[right])
        distances = {
            "original": round(rgb_distance(rgb_left, rgb_right), 1),
            "grayscale": round(rgb_distance(simulate_rgb(rgb_left, "grayscale"), simulate_rgb(rgb_right, "grayscale")), 1),
            "protanopia": round(rgb_distance(simulate_rgb(rgb_left, "protanopia"), simulate_rgb(rgb_right, "protanopia")), 1),
            "deuteranopia": round(rgb_distance(simulate_rgb(rgb_left, "deuteranopia"), simulate_rgb(rgb_right, "deuteranopia")), 1),
            "tritanopia": round(rgb_distance(simulate_rgb(rgb_left, "tritanopia"), simulate_rgb(rgb_right, "tritanopia")), 1),
        }
        min_cvd = min(distances["protanopia"], distances["deuteranopia"], distances["tritanopia"])
        pairwise.append(
            {
                "pair": [left, right],
                "distances": distances,
                "min_colorblind_distance": min_cvd,
                "flag": "warn" if min_cvd < 55 else "pass",
            }
        )
    return {"records": records, "pairwise": pairwise}


def build_palette_sheet(tokens: dict[str, Any], metrics: dict[str, Any], output: Path) -> None:
    colors = tokens["colors"]
    rows = metrics["records"]
    modes = ["original", "grayscale", "protanopia", "deuteranopia", "tritanopia"]
    fig = plt.figure(figsize=(13.5, 7.2), dpi=220)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.text(0.035, 0.940, "SpaceBio-Bench semantic palette QA", fontsize=16, fontweight="bold", color=colors["ink"], transform=ax.transAxes)
    ax.text(
        0.035,
        0.900,
        "Swatches show original and color-vision simulations; contrast ratios are against the paper surface.",
        fontsize=8.8,
        color=colors["muted"],
        transform=ax.transAxes,
    )
    col_x = [0.300, 0.410, 0.520, 0.630, 0.740]
    for idx, mode in enumerate(modes):
        ax.text(col_x[idx], 0.842, mode, fontsize=7.8, color=colors["muted"], ha="center", transform=ax.transAxes)
    y = 0.785
    row_h = 0.080
    for record in rows:
        name = record["name"]
        ax.text(0.035, y + 0.017, name, fontsize=8.8, color=colors["ink"], fontweight="bold", transform=ax.transAxes)
        ax.text(
            0.160,
            y + 0.017,
            f"{record['hex']} | contrast {record['contrast_vs_paper']}",
            fontsize=7.7,
            color=colors["muted"],
            transform=ax.transAxes,
        )
        for idx, mode in enumerate(modes):
            hex_value = colors[name] if mode == "original" else record["simulated_hex"][mode]
            ax.add_patch(plt.Rectangle((col_x[idx] - 0.045, y), 0.090, 0.048, facecolor=hex_value, edgecolor="#FFFFFF", linewidth=1.0, transform=ax.transAxes))
        y -= row_h
    flagged = [item for item in metrics["pairwise"] if item["flag"] == "warn"]
    ax.text(0.035, 0.115, "Pairwise warning threshold: min simulated RGB distance < 55.", fontsize=8.0, color=colors["muted"], transform=ax.transAxes)
    if flagged:
        summary = "; ".join(f"{left}/{right}: {item['min_colorblind_distance']}" for item in flagged[:6] for left, right in [item["pair"]])
    else:
        summary = "No pairwise semantic-color warning under the current threshold."
    ax.text(0.035, 0.080, summary, fontsize=7.7, color=colors["muted"], transform=ax.transAxes)
    fig.savefig(output, dpi=220, facecolor="white")
    plt.close(fig)


def build() -> dict[str, Any]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    tokens = load_tokens()
    metrics = palette_metrics(tokens)
    outputs = {
        "b1_modes": build_mode_sheet(SOURCES["b1_evaluation_layer"], OUT_DIR / "b1_evaluation_layer_color_modes.png", "B1 evaluation-layer color QA", thumb_width=900),
        "family_modes": build_mode_sheet(SOURCES["b1_b4_family_contact_sheet"], OUT_DIR / "b1_b4_family_color_modes.png", "B1-B4 bridge-family color QA", thumb_width=840),
        "palette_sheet": rel(OUT_DIR / "semantic_palette_color_qa.png"),
        "metrics_json": rel(OUT_DIR / "semantic_palette_color_metrics.json"),
        "qa_json": rel(OUT_DIR / "premium_bridge_color_qa.json"),
    }
    build_palette_sheet(tokens, metrics, OUT_DIR / "semantic_palette_color_qa.png")

    low_contrast = [
        record
        for record in metrics["records"]
        if record["name"] not in {"rule"} and record["contrast_vs_paper"] < 3.0
    ]
    pairwise_warnings = [item for item in metrics["pairwise"] if item["flag"] == "warn"]
    qa = {
        "created": CREATED,
        "scope": "B1 evaluation-layer slide plus B1-B4 bridge-family contact sheet",
        "source_images": {name: rel(path) for name, path in SOURCES.items()},
        "outputs": outputs,
        "automatic_checks": {
            "source_images_exist": all(path.exists() for path in SOURCES.values()),
            "semantic_color_records": len(metrics["records"]),
            "text_color_contrast_warnings_vs_paper": low_contrast,
            "pairwise_colorblind_distance_warnings": pairwise_warnings,
        },
        "manual_review": {
            "b1_color_distinction": "pending",
            "b1_grayscale_distinction": "pending",
            "b1_colorblind_distinction": "pending",
            "family_color_consistency": "pending",
            "premium_palette_verdict": "pending",
        },
    }
    write_json(OUT_DIR / "semantic_palette_color_metrics.json", metrics)
    write_json(OUT_DIR / "premium_bridge_color_qa.json", qa)
    return {"outputs": outputs, "qa": qa}


def main() -> None:
    print(json.dumps(build(), indent=2))


if __name__ == "__main__":
    main()
