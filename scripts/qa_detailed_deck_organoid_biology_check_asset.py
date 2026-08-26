#!/usr/bin/env python3
"""Build strict visual QA crops for the organoid biology-check asset."""

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
    / "organoid_biology_check"
)
SLIDE = ASSET_DIR / "organoid_extension_is_a_biology_check_premium.png"
QA_DIR = ASSET_DIR / "qa_crops"
STRICT_DIR = ASSET_DIR / "strict_qa"
QA_DIR.mkdir(parents=True, exist_ok=True)
STRICT_DIR.mkdir(parents=True, exist_ok=True)


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


F_TITLE = font(44, bold=True)
F_LABEL = font(24, bold=True)


def fit(image: Image.Image, size: tuple[int, int]) -> Image.Image:
    image = image.convert("RGB").copy()
    resample = getattr(getattr(Image, "Resampling", Image), "LANCZOS")
    image.thumbnail(size, resample)
    canvas = Image.new("RGB", size, "#0C111A")
    canvas.paste(image, ((size[0] - image.width) // 2, (size[1] - image.height) // 2))
    return canvas


def crop_set(target_dir: Path, specs: list[tuple[str, tuple[int, int, int, int]]], filename: str) -> Path:
    slide = Image.open(SLIDE).convert("RGB")
    entries = []
    for name, box in specs:
        crop = slide.crop(box)
        out = target_dir / f"organoid_biology_check_{name}.png"
        crop.save(out, quality=95)
        entries.append({"name": name, "box": list(box), "path": str(out.relative_to(ROOT))})

    thumb_w, thumb_h = 1480, 620
    pad_x, pad_y = 70, 88
    width = thumb_w + pad_x * 2
    height = 56 + len(specs) * (thumb_h + pad_y) + 40
    canvas = Image.new("RGB", (width, height), "#0C111A")
    draw = ImageDraw.Draw(canvas)
    y = 28
    for name, _box in specs:
        crop = Image.open(target_dir / f"organoid_biology_check_{name}.png")
        draw.text((pad_x, y), name.replace("_", " "), font=F_TITLE, fill="#F3F7FC")
        y += 54
        draw.rounded_rectangle((pad_x - 6, y - 6, pad_x + thumb_w + 6, y + thumb_h + 6), radius=18, outline="#2A394D", width=2)
        canvas.paste(fit(crop, (thumb_w, thumb_h)), (pad_x, y))
        y += thumb_h + pad_y
    sheet = target_dir / filename
    canvas.save(sheet, quality=95)
    (target_dir / "organoid_biology_check_crop_manifest.json").write_text(json.dumps(entries, indent=2) + "\n")
    return sheet


def main() -> None:
    qa_specs = [
        ("full_thumbnail", (0, 0, 3840, 2160)),
        ("title_badges", (90, 42, 2850, 520)),
        ("footprint_panel", (90, 580, 1810, 1420)),
        ("readout_panel", (1810, 580, 3745, 1420)),
        ("bar_chart", (90, 1420, 3745, 1890)),
        ("footer", (90, 1880, 3745, 2105)),
    ]
    strict_specs = [
        ("01_full_thumbnail", (0, 0, 3840, 2160)),
        ("02_header_title", (90, 42, 2940, 520)),
        ("03_footprint_nodes", (440, 770, 1795, 1315)),
        ("04_rosette_card", (145, 770, 540, 1325)),
        ("05_readout_sliders", (1880, 790, 3680, 1320)),
        ("06_bar_values", (120, 1450, 3720, 1856)),
        ("07_footer", (120, 1910, 3720, 2076)),
    ]
    qa_sheet = crop_set(QA_DIR, qa_specs, "organoid_biology_check_qa_crop_contact_sheet.png")
    strict_sheet = crop_set(STRICT_DIR, strict_specs, "organoid_biology_check_STRICT_QA_contact_sheet.png")
    print(json.dumps({"qa_contact_sheet": str(qa_sheet.relative_to(ROOT)), "strict_contact_sheet": str(strict_sheet.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()
