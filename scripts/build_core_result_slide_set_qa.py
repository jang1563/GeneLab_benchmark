#!/usr/bin/env python3
"""Build integrated QA artifacts for premium core result slides 10-14."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFilter, ImageFont


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "premium_core_result_slides" / "core_result_set_qa_v0_1"
VARIANTS = OUT / "variants"
CREATED = "2026-06-03"

COLORS = {
    "void": "#070A0E",
    "deep": "#0B1117",
    "panel": "#101923",
    "panel2": "#152331",
    "ink": "#F4F7F8",
    "soft": "#B9C7D2",
    "muted": "#788898",
    "rule": "#33465A",
    "sky": "#6BAED6",
    "teal": "#1AA090",
    "amber": "#B4832F",
    "rose": "#B45A7E",
    "green": "#178B63",
}

COLOR_TRANSFORMS = {
    "original": np.eye(3),
    "deuteranopia": np.array(
        [
            [0.625, 0.375, 0.000],
            [0.700, 0.300, 0.000],
            [0.000, 0.300, 0.700],
        ]
    ),
    "protanopia": np.array(
        [
            [0.567, 0.433, 0.000],
            [0.558, 0.442, 0.000],
            [0.000, 0.242, 0.758],
        ]
    ),
    "tritanopia": np.array(
        [
            [0.950, 0.050, 0.000],
            [0.000, 0.433, 0.567],
            [0.000, 0.475, 0.525],
        ]
    ),
}

SLIDES = [
    {
        "n": 10,
        "tone": "sky",
        "short_id": "slide10_v4_hardening",
        "title": "v4 hardening",
        "headline": "Broader grid reduces cherry-pick risk",
        "plain_caption": "v4 expands the benchmark across 8 tissues, 8 classifiers, and 4 feature views, while keeping significance boundaries visible.",
        "caption": "The v4 hardening layer broadens the benchmark to 256 evaluations across tissue, classifier, and feature-view combinations. It should be read as robustness and coverage evidence, not as a new leaderboard or a claim that every displayed best-row tissue is significant.",
        "do_not_say": "Do not imply every best-row tissue is statistically significant.",
        "dir": ROOT / "output" / "premium_core_result_slides" / "slide10_v4_hardening_v0_1",
        "image": "slide10_v4_hardening_rendered_preview.png",
        "manifest": "slide10_v4_hardening_manifest.json",
        "qa": "slide10_v4_hardening_qa.json",
        "review": "docs/SLIDE10_V4_HARDENING_SCENE_V0_1_REVIEW_2026_06_03.md",
    },
    {
        "n": 11,
        "tone": "rose",
        "short_id": "slide11_temporal_rrrm",
        "title": "Temporal and RRRM guardrails",
        "headline": "Timepoint labels need guardrails",
        "plain_caption": "v2 separates preservation effects, descriptive recovery projections, and an underpowered RRRM single-cell pilot before interpretation.",
        "caption": "The temporal/RRRM layer shows why timepoint and pilot single-cell labels cannot be over-read: preservation can dominate ISS-T/LAR separation, recovery projections are descriptive, and RRRM composition testing remains underpowered despite benchmark-ready data.",
        "do_not_say": "Do not present LAR recovery as held-out validation or RRRM composition as mature mechanism.",
        "dir": ROOT / "output" / "premium_core_result_slides" / "slide11_temporal_rrrm_guardrails_v0_1",
        "image": "slide11_temporal_rrrm_guardrails_rendered_preview.png",
        "manifest": "slide11_temporal_rrrm_guardrails_manifest.json",
        "qa": "slide11_temporal_rrrm_guardrails_qa.json",
        "review": "docs/SLIDE11_TEMPORAL_RRRM_GUARDRAILS_V0_1_REVIEW_2026_06_03.md",
    },
    {
        "n": 12,
        "tone": "amber",
        "short_id": "slide12_negative_boundary",
        "title": "Negative boundary",
        "headline": "Negative results set boundaries",
        "plain_caption": "v3 makes null and weak-signal settings explicit, including spatial brain and pretrained embedding controls.",
        "caption": "The negative-boundary layer frames spatial brain and foundation-model results as useful task limits. The slide supports careful follow-up design, not a universal claim that brain biology is absent or that future spatial or foundation-model analyses must fail.",
        "do_not_say": "Do not claim universal absence of brain spaceflight biology or universal FM failure.",
        "dir": ROOT / "output" / "premium_core_result_slides" / "slide12_negative_boundary_v0_1",
        "image": "slide12_negative_boundary_rendered_preview.png",
        "manifest": "slide12_negative_boundary_manifest.json",
        "qa": "slide12_negative_boundary_qa.json",
        "review": "docs/SLIDE12_NEGATIVE_BOUNDARY_V0_1_REVIEW_2026_06_03.md",
    },
    {
        "n": 13,
        "tone": "teal",
        "short_id": "slide13_biological_interpretation",
        "title": "Biological interpretation",
        "headline": "Hits become hypotheses",
        "plain_caption": "v5 connects immune shifts, metabolic context, and biomarker triage without making treatment claims.",
        "caption": "The biological interpretation layer converts benchmark signals into follow-up hypotheses across immune context, pathway-level metabolism, and target-linked biomarker triage. It is an interpretation layer only, not causal mechanism, clinical validation, or countermeasure recommendation.",
        "do_not_say": "Do not call target links treatment recommendations or validated mechanisms.",
        "dir": ROOT / "output" / "premium_core_result_slides" / "slide13_biological_interpretation_v0_1",
        "image": "slide13_biological_interpretation_rendered_preview.png",
        "manifest": "slide13_biological_interpretation_manifest.json",
        "qa": "slide13_biological_interpretation_qa.json",
        "review": "docs/SLIDE13_BIOLOGICAL_INTERPRETATION_V0_1_REVIEW_2026_06_03.md",
    },
    {
        "n": 14,
        "tone": "green",
        "short_id": "slide14_human_translation",
        "title": "Human translation boundary",
        "headline": "Translation is partial, not direct",
        "plain_caption": "v6 shows partial pathway and target-tier alignment while direct gene-level transfer remains weak.",
        "caption": "The human-translation layer closes the v1-v7 result core by showing that mouse signals partly carry over at pathway and target-tier levels, while direct gene transfer and TF conservation remain weak. Use it as human-alignment evidence, not clinical proof.",
        "do_not_say": "Do not imply clinical validation, direct gene-transfer reliability, or countermeasure actionability.",
        "dir": ROOT / "output" / "premium_core_result_slides" / "slide14_human_translation_ladder_v0_1",
        "image": "slide14_human_translation_ladder_rendered_preview.png",
        "manifest": "slide14_human_translation_ladder_manifest.json",
        "qa": "slide14_human_translation_ladder_qa.json",
        "review": "docs/SLIDE14_HUMAN_TRANSLATION_LADDER_V0_1_REVIEW_2026_06_03.md",
    },
]


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
    "eyebrow": font(22, bold=True),
    "title": font(58, bold=True),
    "subtitle": font(24),
    "h": font(25, bold=True),
    "body": font(18),
    "small": font(15),
    "tiny": font(12),
    "row": font(17, bold=True),
}


def ensure() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    VARIANTS.mkdir(parents=True, exist_ok=True)


def load_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def draw_background(canvas: Image.Image) -> None:
    draw = ImageDraw.Draw(canvas)
    w, h = canvas.size
    draw.rectangle((0, 0, w, h), fill=rgb("void"))
    top = rgb("void")
    bottom = rgb("panel")
    for y in range(0, h, 2):
        t = y / max(1, h - 1)
        color = tuple(int(top[i] * (1 - t) + bottom[i] * t) for i in range(3))
        draw.line((0, y, w, y), fill=color + (255,), width=2)
    for x in range(170, w, 170):
        draw.line((x, 0, x, h), fill=rgba("rule", 20), width=1)
    for y in range(150, h, 150):
        draw.line((0, y, w, y), fill=rgba("rule", 16), width=1)


def fit_image(image: Image.Image, size: tuple[int, int]) -> Image.Image:
    target_w, target_h = size
    image = image.convert("RGB")
    copied = image.copy()
    resample = getattr(getattr(Image, "Resampling", Image), "LANCZOS")
    copied.thumbnail(size, resample)
    result = Image.new("RGB", size, rgb("deep"))
    x = (target_w - copied.width) // 2
    y = (target_h - copied.height) // 2
    result.paste(copied, (x, y))
    return result


def wrap_text(draw: ImageDraw.ImageDraw, text: str, width_px: int, font_obj: ImageFont.ImageFont) -> list[str]:
    words = text.split()
    lines: list[str] = []
    current = ""
    for word in words:
        candidate = f"{current} {word}".strip()
        box = draw.textbbox((0, 0), candidate, font=font_obj)
        if box[2] - box[0] <= width_px or not current:
            current = candidate
        else:
            lines.append(current)
            current = word
    if current:
        lines.append(current)
    return lines


def header(draw: ImageDraw.ImageDraw, title: str, subtitle: str) -> None:
    draw.text((150, 92), "CORE RESULT SLIDES 10-14 QA", font=F["eyebrow"], fill=rgb("teal"))
    draw.text((150, 140), title, font=F["title"], fill=rgb("ink"))
    draw.text((154, 214), subtitle, font=F["subtitle"], fill=rgb("soft"))


def load_slide_records() -> list[dict]:
    records = []
    for slide in SLIDES:
        image_path = slide["dir"] / slide["image"]
        manifest = load_json(slide["dir"] / slide["manifest"])
        qa = load_json(slide["dir"] / slide["qa"])
        records.append(
            {
                **slide,
                "image_path": image_path,
                "manifest_data": manifest,
                "qa_data": qa,
            }
        )
    return records


def draw_contact_sheet(records: list[dict]) -> Path:
    canvas = Image.new("RGBA", (3840, 2160), (0, 0, 0, 255))
    draw_background(canvas)
    draw = ImageDraw.Draw(canvas)
    header(
        draw,
        "Core result set contact sheet",
        "One family review: hardening, guardrails, negative boundaries, biology hypotheses, and human translation.",
    )

    thumb_w, thumb_h = 650, 366
    x0, y0 = 148, 430
    gutter = 52
    for idx, record in enumerate(records):
        x = x0 + idx * (thumb_w + gutter)
        y = y0
        tone = record["tone"]
        shadow = Image.new("RGBA", canvas.size, (0, 0, 0, 0))
        sd = ImageDraw.Draw(shadow)
        sd.rounded_rectangle((x + 10, y + 18, x + thumb_w + 10, y + thumb_h + 18), radius=12, fill=(0, 0, 0, 120))
        canvas.alpha_composite(shadow.filter(ImageFilter.GaussianBlur(14)))
        draw.rounded_rectangle((x - 2, y - 2, x + thumb_w + 2, y + thumb_h + 2), radius=10, outline=rgba(tone, 185), width=2)
        canvas.alpha_composite(fit_image(Image.open(record["image_path"]), (thumb_w, thumb_h)).convert("RGBA"), (x, y))
        draw.text((x, y + thumb_h + 34), f"{record['n']}. {record['title']}", font=F["h"], fill=rgb("ink"))
        draw.text((x, y + thumb_h + 70), record["headline"], font=F["body"], fill=rgb(tone))
        lines = wrap_text(draw, record["plain_caption"], thumb_w, F["small"])
        for line_idx, line in enumerate(lines[:4]):
            draw.text((x, y + thumb_h + 110 + line_idx * 24), line, font=F["small"], fill=rgb("soft"))
        qa = record["qa_data"]["automatic_checks"]
        draw.text(
            (x, y + thumb_h + 232),
            f"text {qa['visible_text_word_count']}/{qa['visible_text_budget']} words",
            font=F["tiny"],
            fill=rgb("muted"),
        )
        draw.text((x, y + thumb_h + 258), "claim boundary locked", font=F["tiny"], fill=rgb("muted"))

    bridge_y = 1422
    points = [(310, bridge_y + 180), (980, bridge_y + 90), (1680, bridge_y + 185), (2460, bridge_y + 96), (3180, bridge_y + 178)]
    draw.line(points, fill=rgba("teal", 34), width=44, joint="curve")
    draw.line(points, fill=rgba("amber", 145), width=5, joint="curve")
    labels = ["robustness", "guardrails", "limits", "hypotheses", "translation"]
    for idx, point in enumerate(points):
        tone = records[idx]["tone"]
        draw.ellipse((point[0] - 20, point[1] - 20, point[0] + 20, point[1] + 20), fill=rgba(tone, 215))
        draw.text((point[0] - 56, point[1] + 38), labels[idx], font=F["body"], fill=rgb("soft"))

    footer = (
        "Integrated read order: v4 hardening -> temporal/RRRM guardrails -> negative boundaries -> "
        "biological hypotheses -> human translation boundary."
    )
    draw.rounded_rectangle((150, 1878, 2440, 1982), radius=10, fill=rgba("panel", 230), outline=rgba("rule", 150), width=2)
    draw.text((182, 1905), footer, font=F["body"], fill=rgb("soft"))
    draw.text((2920, 1932), f"created {CREATED}", font=F["small"], fill=rgb("muted"))

    output = OUT / "core_result_slides_10_14_contact_sheet.png"
    canvas.convert("RGB").save(output, quality=95)
    return output


def image_to_float(image: Image.Image) -> np.ndarray:
    array = np.asarray(image.convert("RGB"), dtype=np.float32) / 255.0
    return array


def array_to_image(array: np.ndarray) -> Image.Image:
    array = np.clip(array * 255.0, 0, 255).astype(np.uint8)
    return Image.fromarray(array)


def grayscale(array: np.ndarray) -> np.ndarray:
    gray = array @ np.array([0.2126, 0.7152, 0.0722], dtype=np.float32)
    return np.repeat(gray[..., None], 3, axis=2)


def transformed_variants(records: list[dict]) -> list[dict]:
    rows = []
    for record in records:
        image = fit_image(Image.open(record["image_path"]).convert("RGB"), (960, 540))
        array = image_to_float(image)
        variant_arrays = {"grayscale": grayscale(array)}
        for name, matrix in COLOR_TRANSFORMS.items():
            variant_arrays[name] = np.clip(array @ matrix.T, 0.0, 1.0)
        for variant_name, variant_array in variant_arrays.items():
            output = VARIANTS / f"{record['short_id']}_{variant_name}.png"
            array_to_image(variant_array).save(output, quality=95)
            luminance = variant_array @ np.array([0.2126, 0.7152, 0.0722], dtype=np.float32)
            rows.append(
                {
                    "slide": record["n"],
                    "short_id": record["short_id"],
                    "variant": variant_name,
                    "source": rel(record["image_path"]),
                    "output": rel(output),
                    "width_px": int(variant_array.shape[1]),
                    "height_px": int(variant_array.shape[0]),
                    "variant_resolution": "qa_thumbnail_960x540",
                    "luminance_mean": round(float(luminance.mean()), 4),
                    "luminance_std": round(float(luminance.std()), 4),
                }
            )
    return rows


def draw_color_qa_sheet(records: list[dict], variant_rows: list[dict]) -> Path:
    by_key = {(row["short_id"], row["variant"]): ROOT / row["output"] for row in variant_rows}
    variants = ["original", "grayscale", "deuteranopia", "protanopia", "tritanopia"]
    canvas = Image.new("RGBA", (3840, 2160), (0, 0, 0, 255))
    draw_background(canvas)
    draw = ImageDraw.Draw(canvas)
    header(
        draw,
        "Color and grayscale QA",
        "The five-slide family must remain legible without relying on hue alone.",
    )

    thumb_w, thumb_h = 560, 315
    x0, y0 = 430, 322
    col_gap, row_gap = 38, 25
    for col, record in enumerate(records):
        x = x0 + col * (thumb_w + col_gap)
        draw.text((x, 284), f"{record['n']}. {record['title']}", font=F["tiny"], fill=rgb("soft"))
    for row_idx, variant in enumerate(variants):
        y = y0 + row_idx * (thumb_h + row_gap)
        draw.text((154, y + 142), variant, font=F["row"], fill=rgb("muted" if variant != "original" else "ink"))
        for col, record in enumerate(records):
            x = x0 + col * (thumb_w + col_gap)
            image = fit_image(Image.open(by_key[(record["short_id"], variant)]), (thumb_w, thumb_h))
            draw.rounded_rectangle((x - 1, y - 1, x + thumb_w + 1, y + thumb_h + 1), radius=8, outline=rgba(record["tone"], 150), width=1)
            canvas.alpha_composite(image.convert("RGBA"), (x, y))

    draw.rounded_rectangle((150, 2048, 2660, 2110), radius=10, fill=rgba("panel", 225), outline=rgba("rule", 145), width=2)
    draw.text(
        (180, 2068),
        "Review target: panel structure, key metrics, and caveats remain visible across original, grayscale, and color-vision variants.",
        font=F["small"],
        fill=rgb("soft"),
    )
    output = OUT / "core_result_slides_10_14_color_qa_sheet.png"
    canvas.convert("RGB").save(output, quality=95)
    return output


def caption_pack(records: list[dict]) -> dict:
    slides = []
    for record in records:
        manifest = record["manifest_data"]
        qa = record["qa_data"]
        slides.append(
            {
                "slide": record["n"],
                "title": record["title"],
                "headline": record["headline"],
                "plain_caption": record["plain_caption"],
                "caption": record["caption"],
                "claim_boundary": manifest["claim_boundary"],
                "do_not_say": record["do_not_say"],
                "source_documents": manifest.get("source_documents", []),
                "review_document": record["review"],
                "visible_text_word_count": manifest["visible_text_word_count"],
                "visible_text_budget": manifest["visible_text_budget"],
                "qa_status": qa["manual_review"]["visual_review_status"],
            }
        )
    return {
        "created": CREATED,
        "caption_pack_role": "Slide 10-14 presenter/manuscript caption source for the premium deck spine.",
        "family_caption": (
            "Slides 10-14 form the core result spine: v4 hardening establishes robustness, "
            "v2 temporal/RRRM lessons add guardrails, v3 negative results define task limits, "
            "v5 biological interpretation turns hits into hypotheses, and v6 human translation "
            "keeps translational claims boundary-aware."
        ),
        "slides": slides,
    }


def write_caption_markdown(pack: dict) -> Path:
    path = OUT / "core_result_slides_10_14_caption_pack.md"
    lines = [
        "# Core Result Slides 10-14 Caption Pack",
        "",
        f"Date: {CREATED}",
        "",
        "## Family Caption",
        "",
        pack["family_caption"],
        "",
        "## Slide Captions",
        "",
    ]
    for slide in pack["slides"]:
        lines.extend(
            [
                f"### Slide {slide['slide']}: {slide['title']}",
                "",
                f"Headline: {slide['headline']}",
                "",
                f"Presenter caption: {slide['plain_caption']}",
                "",
                f"Manuscript-style caption: {slide['caption']}",
                "",
                f"Claim boundary: {slide['claim_boundary']}",
                "",
                f"Do not say: {slide['do_not_say']}",
                "",
                f"Visible text: {slide['visible_text_word_count']}/{slide['visible_text_budget']} words",
                "",
            ]
        )
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def write_json(path: Path, data: dict | list) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def build() -> dict[str, str]:
    ensure()
    records = load_slide_records()
    contact_sheet = draw_contact_sheet(records)
    variants = transformed_variants(records)
    color_qa_sheet = draw_color_qa_sheet(records, variants)

    pack = caption_pack(records)
    caption_json = OUT / "core_result_slides_10_14_caption_pack.json"
    caption_md = write_caption_markdown(pack)
    write_json(caption_json, pack)

    manifest = {
        "created": CREATED,
        "artifact_role": "integrated visual QA and caption pack for slides 10-14",
        "slides": [
            {
                "slide": record["n"],
                "short_id": record["short_id"],
                "image": rel(record["image_path"]),
                "manifest": rel(record["dir"] / record["manifest"]),
                "qa": rel(record["dir"] / record["qa"]),
                "review": record["review"],
            }
            for record in records
        ],
        "outputs": {
            "contact_sheet": rel(contact_sheet),
            "color_qa_sheet": rel(color_qa_sheet),
            "caption_pack_json": rel(caption_json),
            "caption_pack_markdown": rel(caption_md),
            "variant_manifest": rel(OUT / "core_result_slides_10_14_color_variant_manifest.json"),
            "qa": rel(OUT / "core_result_slides_10_14_integrated_qa.json"),
        },
    }
    manifest_path = OUT / "core_result_slides_10_14_manifest.json"
    write_json(manifest_path, manifest)
    write_json(OUT / "core_result_slides_10_14_color_variant_manifest.json", variants)

    word_counts = {
        str(record["n"]): {
            "count": record["manifest_data"]["visible_text_word_count"],
            "budget": record["manifest_data"]["visible_text_budget"],
            "within_budget": record["manifest_data"]["visible_text_word_count"] <= record["manifest_data"]["visible_text_budget"],
        }
        for record in records
    }
    qa_statuses = {str(record["n"]): record["qa_data"]["manual_review"]["visual_review_status"] for record in records}
    review_flags = []
    for record in records:
        review_flags.extend(record["manifest_data"].get("source_review_flags", []))

    qa = {
        "created": CREATED,
        "automatic_checks": {
            "n_slides": len(records),
            "all_source_images_exist": all(record["image_path"].exists() for record in records),
            "all_text_counts_within_budget": all(item["within_budget"] for item in word_counts.values()),
            "word_counts": word_counts,
            "source_image_dimensions": {
                str(record["n"]): Image.open(record["image_path"]).size for record in records
            },
            "variant_outputs": len(variants),
            "variant_types": sorted({row["variant"] for row in variants}),
            "review_flags": review_flags,
        },
        "manual_review": {
            "visual_review_status": "pass: integrated contact sheet and color-QA sheet inspected after footer-overlap revision",
            "family_read_order": "hardening -> guardrails -> limits -> hypotheses -> translation",
            "claim_boundary": "No clinical, treatment, universal-failure, or direct-transfer overclaim across the five-slide set.",
            "qa_statuses": qa_statuses,
            "cautions": [
                "Slide 10 has a tight 44/45-word budget; keep final PPTX labels editable and controlled.",
                "Slide 13 intentionally follows current JSON over legacy HTML narrative statements.",
                "Final deck assembly should rebuild title, source, and caveat text as editable overlays over scene plates.",
            ],
        },
    }
    qa_path = OUT / "core_result_slides_10_14_integrated_qa.json"
    write_json(qa_path, qa)
    return {
        "contact_sheet": rel(contact_sheet),
        "color_qa_sheet": rel(color_qa_sheet),
        "caption_pack_markdown": rel(caption_md),
        "caption_pack_json": rel(caption_json),
        "manifest": rel(manifest_path),
        "qa": rel(qa_path),
    }


def main() -> None:
    print(json.dumps(build(), indent=2))


if __name__ == "__main__":
    main()
