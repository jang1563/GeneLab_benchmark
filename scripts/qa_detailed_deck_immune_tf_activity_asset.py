#!/usr/bin/env python3
"""Build visual QA crops for slide 38 immune/TF activity asset."""

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
    / "immune_tf_activity"
)
SOURCE = ASSET_DIR / "immune_tf_activity_prioritize_followup_premium.png"
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
    canvas = Image.new("RGB", size, "#0C111A")
    canvas.paste(image, ((size[0] - image.width) // 2, (size[1] - image.height) // 2))
    return canvas


def write_contact_sheet(crops: list[dict[str, object]], directory: Path, sheet_name: str) -> Path:
    width = 1700
    thumb_w = 1560
    y = 42
    blocks: list[tuple[str, Image.Image, int]] = []
    for crop in crops:
        name = str(crop["name"])
        image = Image.open(directory / f"immune_tf_activity_{name}.png")
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
        out = target_dir / f"immune_tf_activity_{name}.png"
        source.crop(box).save(out, quality=95)
        manifest.append({"name": name, "box": box, "path": str(out.relative_to(ROOT))})
    sheet = write_contact_sheet(specs, target_dir, sheet_name)
    (target_dir / "immune_tf_activity_crop_manifest.json").write_text(
        json.dumps({"source": str(SOURCE.relative_to(ROOT)), "crops": manifest, "contact_sheet": str(sheet.relative_to(ROOT))}, indent=2) + "\n",
        encoding="utf-8",
    )
    return sheet


def main() -> None:
    qa_specs = [
        {"name": "01_header_badges", "box": (100, 40, 3740, 235)},
        {"name": "02_title_subtitle", "box": (120, 350, 3120, 625)},
        {"name": "03_immune_panel", "box": (120, 610, 1260, 1685)},
        {"name": "04_tf_panel", "box": (1265, 610, 2410, 1685)},
        {"name": "05_priority_panel", "box": (2415, 610, 3745, 1685)},
        {"name": "06_reader_rule", "box": (120, 1680, 3720, 1845)},
        {"name": "07_footer_scope", "box": (120, 1865, 3720, 2095)},
        {"name": "08_full_downscaled", "box": (0, 0, 3840, 2160)},
    ]
    strict_specs = [
        {"name": "01_full_thumbnail", "box": (0, 0, 3840, 2160)},
        {"name": "02_title_subtitle", "box": (120, 350, 3120, 625)},
        {"name": "03_immune_ladder", "box": (190, 1025, 1195, 1395)},
        {"name": "04_immune_cell_chips", "box": (190, 1410, 1195, 1605)},
        {"name": "05_tf_ladder", "box": (1335, 1025, 2340, 1395)},
        {"name": "06_tf_examples", "box": (1335, 1410, 2340, 1605)},
        {"name": "07_priority_scatter", "box": (2480, 1015, 3685, 1485)},
        {"name": "08_priority_queue", "box": (2480, 1490, 3685, 1675)},
        {"name": "09_reader_rule_arrows", "box": (150, 1695, 3690, 1828)},
        {"name": "10_footer", "box": (120, 1865, 3720, 2095)},
    ]
    qa_sheet = crop_set(QA_DIR, qa_specs, "immune_tf_activity_qa_crop_contact_sheet.png")
    strict_sheet = crop_set(STRICT_DIR, strict_specs, "immune_tf_activity_STRICT_QA_contact_sheet.png")
    print(json.dumps({"qa_contact_sheet": str(qa_sheet.relative_to(ROOT)), "strict_contact_sheet": str(strict_sheet.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()
