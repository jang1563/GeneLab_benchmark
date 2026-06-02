#!/usr/bin/env python3
"""Build a 24-slide v1-v9 deck-spine contact sheet.

This is a planning board for deck assembly. It keeps the benchmark result core
visible while showing where v9 organoid and OSD-120 methods bridges are allowed
to enter.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont, ImageOps


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "v1_v9_deck_spine_contact_sheet_v0_1"
QA = OUT / "qa"
CREATED = "2026-06-02"

COLORS = {
    "void": "#070A0E",
    "deep": "#0B1117",
    "panel": "#111A24",
    "panel2": "#172330",
    "ink": "#F4F7F8",
    "soft": "#B8C6D1",
    "muted": "#778694",
    "rule": "#33465A",
    "teal": "#1AA090",
    "sky": "#6BAED6",
    "blue": "#2D6F9F",
    "green": "#178B63",
    "amber": "#B4832F",
    "rose": "#B45A7E",
    "red": "#B23A3A",
}


def rgb(token: str) -> tuple[int, int, int]:
    value = COLORS[token].lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


def rgba(token: str, alpha: int) -> tuple[int, int, int, int]:
    return rgb(token) + (alpha,)


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def font(size: int, *, bold: bool = False):
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size=size)
        except OSError:
            continue
    return ImageFont.load_default()


F = {
    "title": font(58, bold=True),
    "subtitle": font(26),
    "section": font(20, bold=True),
    "num": font(30, bold=True),
    "h": font(21, bold=True),
    "body": font(15),
    "tiny": font(12),
}


def p(path: str) -> Path:
    return ROOT / path


SLIDES = [
    {
        "n": 1,
        "title": "SpaceBio-Bench",
        "section": "Opening",
        "claim": "Mission-held-out benchmark platform thesis.",
        "asset_status": "new title visual",
        "tone": "teal",
    },
    {
        "n": 2,
        "title": "External gap",
        "section": "Opening",
        "claim": "Distinct niche without claiming first AI benchmark.",
        "asset_status": "new positioning matrix",
        "tone": "teal",
    },
    {
        "n": 3,
        "title": "Project evolution",
        "section": "Opening",
        "claim": "v1-v9 timeline from analysis to platform.",
        "asset_status": "new timeline",
        "tone": "teal",
    },
    {
        "n": 4,
        "title": "Evaluation layer",
        "section": "Core methods",
        "claim": "Public studies become auditable tasks.",
        "asset_status": "B1/B2 compressed",
        "tone": "sky",
        "asset": p("output/premium_bridge_scenes/b1_evaluation_layer/rendered_preview.png"),
    },
    {
        "n": 5,
        "title": "Mission-held-out guard",
        "section": "Core methods",
        "claim": "A whole mission is held out; processing stays train-only.",
        "asset_status": "B3/B4 compressed",
        "tone": "sky",
        "asset": p("output/premium_bridge_rebuild_scenes/b3_mission_held_out_premium/rendered_preview.png"),
    },
    {
        "n": 6,
        "title": "Feature layers",
        "section": "Core methods",
        "claim": "Genes versus pathway summaries before results.",
        "asset_status": "needed bridge",
        "tone": "sky",
    },
    {
        "n": 7,
        "title": "Tissue hierarchy",
        "section": "v1-v7 results",
        "claim": "Cross-mission transfer is tissue-specific.",
        "asset_status": "existing premium scene",
        "tone": "blue",
        "asset": p("output/premium_slide_scenes/fig1_tissue_transfer_layered_scene.png"),
    },
    {
        "n": 8,
        "title": "Pathway rescue",
        "section": "v1-v7 results",
        "claim": "Pathway abstraction can rescue selected weak tasks.",
        "asset_status": "existing premium scene",
        "tone": "blue",
        "asset": p("output/premium_slide_scenes/fig2_pathway_layered_scene.png"),
    },
    {
        "n": 9,
        "title": "Model comparison",
        "section": "v1-v7 results",
        "claim": "Scale alone does not solve transfer.",
        "asset_status": "existing premium scene",
        "tone": "blue",
        "asset": p("output/premium_slide_scenes/fig3_model_tier_layered_scene.png"),
    },
    {
        "n": 10,
        "title": "v4 hardening",
        "section": "v1-v7 results",
        "claim": "Broader tissue/method grid reduces cherry-pick risk.",
        "asset_status": "needs static export",
        "tone": "blue",
    },
    {
        "n": 11,
        "title": "Temporal and RRRM lessons",
        "section": "v1-v7 results",
        "claim": "Recovery, preservation, and single-cell lessons.",
        "asset_status": "needs static export",
        "tone": "blue",
    },
    {
        "n": 12,
        "title": "Negative results",
        "section": "v1-v7 results",
        "claim": "Failure modes are part of benchmark value.",
        "asset_status": "needs static export",
        "tone": "blue",
    },
    {
        "n": 13,
        "title": "Biological interpretation",
        "section": "v1-v7 results",
        "claim": "Immune, TF, metabolism, and biomarker interpretation.",
        "asset_status": "needs panel selection",
        "tone": "blue",
    },
    {
        "n": 14,
        "title": "Human translation",
        "section": "v1-v7 results",
        "claim": "Pathway/target evidence, not clean gene transfer.",
        "asset_status": "needs static export",
        "tone": "blue",
    },
    {
        "n": 15,
        "title": "v7.1 boundary",
        "section": "v1-v7 results",
        "claim": "Canonical public-safe counts and claim discipline.",
        "asset_status": "release-boundary slide",
        "tone": "rose",
    },
    {
        "n": 16,
        "title": "v8 SpaceMed",
        "section": "v8 incubator",
        "claim": "Translation is hypothesis generation only.",
        "asset_status": "existing PNG",
        "tone": "amber",
        "asset": p("v8/figures/Figure1_Species_Transfer.png"),
    },
    {
        "n": 17,
        "title": "v8 claim boundary",
        "section": "v8 incubator",
        "claim": "No countermeasure recommendation or Mars point claim.",
        "asset_status": "boundary overlay needed",
        "tone": "amber",
        "asset": p("v8/figures/Figure2_Stressor_Decomposition.png"),
    },
    {
        "n": 18,
        "title": "v9 platform",
        "section": "v9 platform",
        "claim": "Package, manifests, evaluator, and run records.",
        "asset_status": "existing document scene",
        "tone": "green",
        "asset": p("output/premium_v9_document_scenes/v9_platform_provenance_document_scene.png"),
    },
    {
        "n": 19,
        "title": "Public bulk alpha",
        "section": "v9 platform",
        "claim": "Metadata-ready, payload hashes blocked.",
        "asset_status": "existing document scene",
        "tone": "green",
        "asset": p("output/premium_v9_document_scenes/v9_public_bulk_boundary_document_scene.png"),
    },
    {
        "n": 20,
        "title": "Organoid extension",
        "section": "v9 extension",
        "claim": "Source records become local expression matrices.",
        "asset_status": "selected bridge",
        "tone": "teal",
        "asset": p("output/biovis_organoid_audience_matrix_proof_v0_2/panels/01_dark_organoid_clean_source_to_matrix.png"),
    },
    {
        "n": 21,
        "title": "Multispecies task check",
        "section": "v9 extension",
        "claim": "OSD-120 is same-study, not mission-held-out.",
        "asset_status": "selected bridge",
        "tone": "green",
        "asset": p("output/biovis_osd120_audience_split_proof_v0_2/panels/01_dark_osd120_clean_source_to_task.png"),
    },
    {
        "n": 22,
        "title": "Single-cell extension",
        "section": "v9 extension",
        "claim": "Inventory and metric spec exist; payload blocker remains.",
        "asset_status": "needs blocker visual",
        "tone": "green",
    },
    {
        "n": 23,
        "title": "Claim boundary",
        "section": "Close",
        "claim": "Separate completed results from draft tracks.",
        "asset_status": "needed matrix",
        "tone": "rose",
    },
    {
        "n": 24,
        "title": "Roadmap",
        "section": "Close",
        "claim": "Payload freeze, QA, release, and manuscript plan.",
        "asset_status": "needed roadmap",
        "tone": "rose",
    },
]


def ensure() -> None:
    QA.mkdir(parents=True, exist_ok=True)


def draw_wrapped(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    text: str,
    font_obj,
    fill,
    *,
    width: int,
    line_height: int,
    max_lines: int,
) -> int:
    words = text.split()
    lines: list[str] = []
    current = ""
    for word in words:
        candidate = f"{current} {word}".strip()
        if draw.textlength(candidate, font=font_obj) <= width:
            current = candidate
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)
    if len(lines) > max_lines:
        lines = lines[:max_lines]
        lines[-1] = lines[-1].rstrip(".") + "..."
    x, y = xy
    for idx, line in enumerate(lines):
        draw.text((x, y + idx * line_height), line, font=font_obj, fill=fill)
    return y + len(lines) * line_height


def fit_cover(image: Image.Image, size: tuple[int, int]) -> Image.Image:
    return ImageOps.fit(image.convert("RGB"), size, method=Image.Resampling.LANCZOS, centering=(0.5, 0.48))


def draw_background(canvas: Image.Image) -> None:
    draw = ImageDraw.Draw(canvas)
    w, h = canvas.size
    draw.rectangle((0, 0, w, h), fill=rgb("void"))
    for y in range(0, h, 2):
        t = y / max(1, h - 1)
        color = tuple(int(rgb("void")[i] * (1 - t) + rgb("panel")[i] * t) for i in range(3))
        draw.line((0, y, w, y), fill=color + (255,), width=2)
    for x in range(180, w, 180):
        draw.line((x, 0, x, h), fill=rgba("rule", 24), width=1)
    for y in range(160, h, 160):
        draw.line((0, y, w, y), fill=rgba("rule", 22), width=1)
    for idx, radius in enumerate([920, 1180, 1440]):
        center = (int(w * 0.72), int(h * 0.16))
        bbox = (
            center[0] - radius,
            center[1] - int(radius * 0.33),
            center[0] + radius,
            center[1] + int(radius * 0.33),
        )
        draw.arc(bbox, 198, 350, fill=rgba("sky", 44 - idx * 10), width=3)


def draw_placeholder(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], tone: str, text: str) -> None:
    x, y, w, h = box
    draw.rounded_rectangle((x, y, x + w, y + h), radius=10, fill=rgba("panel2", 220), outline=rgba(tone, 120), width=2)
    mid_y = y + h // 2
    for idx in range(5):
        cx = x + 70 + idx * ((w - 140) // 4)
        draw.ellipse((cx - 9, mid_y - 9, cx + 9, mid_y + 9), fill=rgba(tone, 185))
        if idx:
            px = x + 70 + (idx - 1) * ((w - 140) // 4)
            draw.line((px + 12, mid_y, cx - 12, mid_y), fill=rgba(tone, 120), width=3)
    draw_wrapped(draw, (x + 22, y + 24), text, F["body"], rgb("soft"), width=w - 44, line_height=20, max_lines=2)


def draw_thumb(canvas: Image.Image, slide: dict, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(canvas)
    x, y, w, h = box
    tone = slide["tone"]
    asset = slide.get("asset")
    if asset and asset.exists():
        shadow = Image.new("RGBA", canvas.size, (0, 0, 0, 0))
        sd = ImageDraw.Draw(shadow)
        sd.rounded_rectangle((x + 10, y + 12, x + w + 10, y + h + 12), radius=13, fill=(0, 0, 0, 118))
        shadow = shadow.filter(ImageFilter.GaussianBlur(12))
        canvas.alpha_composite(shadow)
        thumb = fit_cover(Image.open(asset), (w, h)).convert("RGBA")
        mask = Image.new("L", (w, h), 0)
        md = ImageDraw.Draw(mask)
        md.rounded_rectangle((0, 0, w, h), radius=10, fill=255)
        canvas.paste(thumb, (x, y), mask)
        draw.rounded_rectangle((x, y, x + w, y + h), radius=10, outline=rgba(tone, 220), width=3)
    else:
        draw_placeholder(draw, box, tone, slide["asset_status"])


def draw_tile(canvas: Image.Image, slide: dict, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(canvas)
    x, y, w, h = box
    tone = slide["tone"]
    draw.rounded_rectangle((x, y, x + w, y + h), radius=12, fill=rgba("deep", 232), outline=rgba("rule", 170), width=2)
    draw.rectangle((x, y, x + w, y + 7), fill=rgba(tone, 235))

    draw.text((x + 18, y + 20), f"{slide['n']:02d}", font=F["num"], fill=rgb(tone))
    draw.text((x + 80, y + 24), slide["section"], font=F["tiny"], fill=rgb("muted"))
    draw_wrapped(draw, (x + 80, y + 46), slide["title"], F["h"], rgb("ink"), width=w - 102, line_height=24, max_lines=2)

    thumb_box = (x + 18, y + 100, w - 36, 152)
    draw_thumb(canvas, slide, thumb_box)

    y2 = y + 270
    draw_wrapped(draw, (x + 18, y2), slide["claim"], F["body"], rgb("soft"), width=w - 36, line_height=20, max_lines=2)
    status = slide["asset_status"]
    status_w = min(w - 36, int(draw.textlength(status, font=F["tiny"])) + 28)
    draw.rounded_rectangle((x + 18, y + h - 42, x + 18 + status_w, y + h - 18), radius=9, fill=rgba("panel2", 235), outline=rgba(tone, 180), width=1)
    draw.text((x + 32, y + h - 36), status, font=F["tiny"], fill=rgb("soft"))


def render_contact_sheet() -> Path:
    canvas = Image.new("RGBA", (3840, 2160), (0, 0, 0, 255))
    draw_background(canvas)
    draw = ImageDraw.Draw(canvas)

    draw.text((105, 74), "v1-v9 deck spine contact sheet", font=F["title"], fill=rgb("ink"))
    draw_wrapped(
        draw,
        (108, 148),
        "24-slide full-talk order. B1-B4 are compressed into two early methods slides; organoid and OSD-120 stay in the v9 extension section.",
        F["subtitle"],
        rgb("soft"),
        width=2520,
        line_height=34,
        max_lines=2,
    )
    draw.rounded_rectangle((2860, 76, 3715, 180), radius=16, fill=rgba("panel", 220), outline=rgba("rule", 170), width=2)
    draw.text((2888, 100), "Guardrail", font=F["section"], fill=rgb("ink"))
    draw.text((2888, 132), "v1-v7 result core first; extension proofs at slides 20-21 only.", font=F["body"], fill=rgb("soft"))

    margin_x = 105
    margin_y = 250
    gutter_x = 24
    gutter_y = 28
    tile_w = 590
    tile_h = 405

    for idx, slide in enumerate(SLIDES):
        row = idx // 6
        col = idx % 6
        box = (margin_x + col * (tile_w + gutter_x), margin_y + row * (tile_h + gutter_y), tile_w, tile_h)
        draw_tile(canvas, slide, box)

    footer_y = 2024
    draw.line((105, footer_y - 28, 3715, footer_y - 28), fill=rgba("rule", 160), width=1)
    draw.text((105, footer_y), "Next work: feature-layer bridge, model/metric caption bridge, v9 public bulk boundary QA, slide-level insertion QA.", font=F["body"], fill=rgb("soft"))
    draw.text((2918, footer_y), "review artifact, not final slide", font=F["body"], fill=rgb("muted"))

    out = OUT / "v1_v9_deck_spine_contact_sheet.png"
    canvas.convert("RGB").save(out, quality=95)
    return out


def render_grayscale(contact_sheet: Path) -> Path:
    out = QA / "v1_v9_deck_spine_contact_sheet_grayscale.png"
    Image.open(contact_sheet).convert("L").convert("RGB").save(out, quality=95)
    return out


def write_manifest(contact_sheet: Path, grayscale: Path) -> tuple[Path, Path]:
    slide_rows = []
    for slide in SLIDES:
        slide_rows.append(
            {
                "slide_number": slide["n"],
                "title": slide["title"],
                "section": slide["section"],
                "claim": slide["claim"],
                "asset_status": slide["asset_status"],
                "source_asset": rel(slide["asset"]) if slide.get("asset") else "",
            }
        )

    order_csv = OUT / "v1_v9_deck_spine_order.csv"
    with order_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(slide_rows[0].keys()))
        writer.writeheader()
        writer.writerows(slide_rows)

    manifest = {
        "created": CREATED,
        "purpose": "24-slide v1-v9 deck-spine planning contact sheet before deck assembly.",
        "outputs": {
            "contact_sheet": rel(contact_sheet),
            "grayscale_qa": rel(grayscale),
            "slide_order_csv": rel(order_csv),
        },
        "slide_count": len(SLIDES),
        "spine_decision": "B1-B4 compressed into two early methods slides for a 24-slide full-talk version; expand to 26 slides only if methods explanation becomes the talk focus.",
        "extension_guardrail": "Organoid and OSD-120 bridge assets appear only in the v9 extension section, slides 20 and 21.",
        "source_docs": [
            "docs/V1_V9_DECK_SPINE_METHODS_BRIDGE_PLACEMENT_2026_06_02.md",
            "docs/V1_V9_PRESENTATION_AND_MANUSCRIPT_MASTER_OUTLINE_2026_05_31.md",
            "docs/V1_V9_SLIDE_ASSET_MANIFEST_2026_05_31.md",
            "docs/PREMIUM_BRIDGE_METHODS_NARRATION_PACK_B1_B4_2026_06_02.md",
        ],
        "slides": slide_rows,
    }

    manifest_path = OUT / "v1_v9_deck_spine_contact_sheet_manifest.json"
    with manifest_path.open("w") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")
    return manifest_path, order_csv


def main() -> None:
    ensure()
    contact_sheet = render_contact_sheet()
    grayscale = render_grayscale(contact_sheet)
    manifest, order_csv = write_manifest(contact_sheet, grayscale)
    print(
        json.dumps(
            {
                "contact_sheet": rel(contact_sheet),
                "grayscale": rel(grayscale),
                "manifest": rel(manifest),
                "order_csv": rel(order_csv),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
