#!/usr/bin/env python3
"""Build visual QA crops for slide 3 evidence-layer asset."""

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
    / "evidence_layers"
)
SOURCE = ASSET_DIR / "evidence_layers_set_claim_strength_premium.png"
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
        image = Image.open(directory / f"evidence_layers_{name}.png")
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
        out = target_dir / f"evidence_layers_{name}.png"
        source.crop(box).save(out, quality=95)
        manifest.append({"name": name, "box": box, "path": str(out.relative_to(ROOT))})
    sheet = write_contact_sheet(specs, target_dir, sheet_name)
    (target_dir / "evidence_layers_crop_manifest.json").write_text(
        json.dumps({"source": str(SOURCE.relative_to(ROOT)), "crops": manifest, "contact_sheet": str(sheet.relative_to(ROOT))}, indent=2) + "\n",
        encoding="utf-8",
    )
    return sheet


def main() -> None:
    qa_specs = [
        {"name": "01_header_badges", "box": (100, 40, 3740, 235)},
        {"name": "02_title_subtitle", "box": (120, 350, 1980, 620)},
        {"name": "03_claim_ladder_top", "box": (120, 650, 2140, 1125)},
        {"name": "04_claim_ladder_bottom", "box": (120, 1080, 2140, 1645)},
        {"name": "05_version_layer_map", "box": (2160, 650, 3720, 1645)},
        {"name": "06_metric_strip", "box": (130, 1628, 3710, 1825)},
        {"name": "07_footer_scope", "box": (130, 1870, 3710, 2095)},
        {"name": "08_full_downscaled", "box": (0, 0, 3840, 2160)},
    ]
    strict_specs = [
        {"name": "01_full_thumbnail", "box": (0, 0, 3840, 2160)},
        {"name": "02_header_badges", "box": (100, 40, 3740, 235)},
        {"name": "03_title_subtitle", "box": (120, 350, 1980, 620)},
        {"name": "04_ladder_layers_1_2", "box": (175, 865, 2080, 1160)},
        {"name": "05_ladder_layers_3_5", "box": (175, 1120, 2080, 1588)},
        {"name": "06_v1_v7_card", "box": (2200, 870, 3660, 1175)},
        {"name": "07_v8_v9_cards", "box": (2200, 1160, 3660, 1495)},
        {"name": "08_reader_rule", "box": (2200, 1470, 3660, 1595)},
        {"name": "09_metric_strip", "box": (130, 1628, 3710, 1825)},
        {"name": "10_footer_scope", "box": (130, 1870, 3710, 2095)},
    ]
    qa_sheet = crop_set(QA_DIR, qa_specs, "evidence_layers_qa_crop_contact_sheet.png")
    strict_sheet = crop_set(STRICT_DIR, strict_specs, "evidence_layers_STRICT_QA_contact_sheet.png")
    print(json.dumps({"qa_contact_sheet": str(qa_sheet.relative_to(ROOT)), "strict_contact_sheet": str(strict_sheet.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()
