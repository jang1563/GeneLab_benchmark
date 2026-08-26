#!/usr/bin/env python3
"""Build visual QA crops for slide 5 source-to-matrix asset."""

from __future__ import annotations

import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent.parent
ASSET_DIR = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
    / "assets"
    / "source_to_matrix"
)
SOURCE = ASSET_DIR / "from_source_record_to_matrix_premium.png"
QA_DIR = ASSET_DIR / "qa_crops"
STRICT_DIR = ASSET_DIR / "strict_qa"


def font(size: int, *, bold: bool = False) -> ImageFont.ImageFont:
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Supplemental/Helvetica Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Helvetica.ttf",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size)
        except OSError:
            continue
    return ImageFont.load_default()


F = {"label": font(28, bold=True)}


def fit(image: Image.Image, size: tuple[int, int]) -> Image.Image:
    image = image.convert("RGB").copy()
    resample = getattr(getattr(Image, "Resampling", Image), "LANCZOS")
    image.thumbnail(size, resample)
    result = Image.new("RGB", size, "#0C111A")
    result.paste(image, ((size[0] - image.width) // 2, (size[1] - image.height) // 2))
    return result


def write_contact_sheet(crops: list[dict[str, object]], directory: Path, sheet_name: str) -> Path:
    width = 1700
    thumb_w = 1560
    y = 42
    blocks: list[tuple[str, Image.Image, int]] = []
    for crop in crops:
        name = str(crop["name"])
        image = Image.open(directory / f"source_to_matrix_{name}.png")
        target_h = 420 if image.width > image.height * 1.8 else 620
        fitted = fit(image, (thumb_w, target_h))
        blocks.append((name, fitted, y))
        y += fitted.height + 98

    canvas = Image.new("RGB", (width, y + 20), "#0C111A")
    draw = ImageDraw.Draw(canvas)
    for name, image, top in blocks:
        draw.text((70, top - 32), name.replace("_", " "), font=F["label"], fill="#F3F7FC")
        draw.rounded_rectangle((66, top - 2, 1634, top + image.height + 2), radius=18, outline="#2A394D", width=2)
        canvas.paste(image, (70, top))
    out = directory / sheet_name
    canvas.save(out, quality=95)
    return out


def crop_set(target_dir: Path, specs: list[dict[str, object]], sheet_name: str) -> Path:
    target_dir.mkdir(parents=True, exist_ok=True)
    source = Image.open(SOURCE).convert("RGB")
    manifest = []
    for spec in specs:
        name = str(spec["name"])
        box = tuple(int(v) for v in spec["box"])
        out = target_dir / f"source_to_matrix_{name}.png"
        source.crop(box).save(out, quality=95)
        manifest.append({"name": name, "box": box, "path": str(out.relative_to(ROOT))})
    sheet = write_contact_sheet(specs, target_dir, sheet_name)
    (target_dir / "source_to_matrix_crop_manifest.json").write_text(
        json.dumps({"source": str(SOURCE.relative_to(ROOT)), "crops": manifest, "contact_sheet": str(sheet.relative_to(ROOT))}, indent=2) + "\n",
        encoding="utf-8",
    )
    return sheet


def main() -> None:
    qa_specs = [
        {"name": "01_header_badges", "box": (100, 40, 3740, 235)},
        {"name": "02_title_subtitle", "box": (120, 350, 2280, 620)},
        {"name": "03_source_panel", "box": (120, 650, 1065, 1625)},
        {"name": "04_matrix_alignment_panel", "box": (1050, 650, 2705, 1625)},
        {"name": "05_output_surface_panel", "box": (2700, 650, 3725, 1625)},
        {"name": "06_stage_pipeline", "box": (130, 1615, 3710, 1840)},
        {"name": "07_footer_scope", "box": (130, 1870, 3710, 2095)},
        {"name": "08_full_downscaled", "box": (0, 0, 3840, 2160)},
    ]
    strict_specs = [
        {"name": "01_full_thumbnail", "box": (0, 0, 3840, 2160)},
        {"name": "02_header_badges", "box": (100, 40, 3740, 235)},
        {"name": "03_title_subtitle", "box": (120, 350, 2280, 620)},
        {"name": "04_source_cards", "box": (185, 880, 980, 1390)},
        {"name": "05_sample_metadata", "box": (185, 1360, 980, 1545)},
        {"name": "06_actual_matrix_crop", "box": (1110, 865, 1905, 1425)},
        {"name": "07_filter_bars", "box": (1900, 900, 2585, 1358)},
        {"name": "08_orientation_callout", "box": (1900, 1350, 2585, 1525)},
        {"name": "09_output_matrix", "box": (2765, 885, 3665, 1168)},
        {"name": "10_next_slide_reader_use", "box": (2765, 1200, 3665, 1548)},
        {"name": "11_stage_pipeline", "box": (130, 1615, 3710, 1840)},
        {"name": "12_footer_scope", "box": (130, 1870, 3710, 2095)},
    ]
    qa_sheet = crop_set(QA_DIR, qa_specs, "source_to_matrix_qa_crop_contact_sheet.png")
    strict_sheet = crop_set(STRICT_DIR, strict_specs, "source_to_matrix_STRICT_QA_contact_sheet.png")
    print(json.dumps({"qa_contact_sheet": str(qa_sheet.relative_to(ROOT)), "strict_contact_sheet": str(strict_sheet.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()
