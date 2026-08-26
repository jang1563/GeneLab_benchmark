#!/usr/bin/env python3
"""Build strict visual QA crops for the bulk adaptation-surface asset."""

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
    / "bulk_surface_hard"
)
SOURCE = ASSET_DIR / "bulk_rnaseq_hard_surface_for_cell_fms_premium.png"
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
        image = Image.open(directory / f"bulk_surface_hard_{name}.png")
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
        out = target_dir / f"bulk_surface_hard_{name}.png"
        source.crop(box).save(out, quality=95)
        manifest.append({"name": name, "box": box, "path": str(out.relative_to(ROOT))})
    sheet = write_contact_sheet(specs, target_dir, sheet_name)
    (target_dir / "bulk_surface_hard_crop_manifest.json").write_text(
        json.dumps({"source": str(SOURCE.relative_to(ROOT)), "crops": manifest, "contact_sheet": str(sheet.relative_to(ROOT))}, indent=2) + "\n"
    )
    return sheet


def main() -> None:
    qa_specs = [
        {"name": "01_header_badges", "box": (120, 40, 3720, 285)},
        {"name": "02_pretraining_surface", "box": (120, 330, 1115, 1840)},
        {"name": "03_adaptation_pinch_points", "box": (1110, 330, 2505, 1840)},
        {"name": "04_benchmark_surface", "box": (2490, 330, 3730, 1840)},
        {"name": "05_footer", "box": (120, 1850, 3720, 2125)},
        {"name": "06_full_downscaled", "box": (0, 0, 3840, 2160)},
    ]
    strict_specs = [
        {"name": "01_full_thumbnail", "box": (0, 0, 3840, 2160)},
        {"name": "02_header_title_badges", "box": (120, 40, 3720, 285)},
        {"name": "03_cell_cloud", "box": (170, 360, 1060, 960)},
        {"name": "04_pretraining_cards", "box": (170, 920, 1060, 1810)},
        {"name": "05_bottleneck_top", "box": (1170, 360, 2435, 835)},
        {"name": "06_bottleneck_middle", "box": (1170, 780, 2435, 1350)},
        {"name": "07_bottleneck_bottom", "box": (1170, 1300, 2435, 1810)},
        {"name": "08_bulk_vector", "box": (2550, 370, 3655, 1015)},
        {"name": "09_mission_and_metrics", "box": (2550, 970, 3655, 1525)},
        {"name": "10_right_readout", "box": (2550, 1475, 3655, 1810)},
        {"name": "11_footer_source", "box": (130, 1860, 3710, 2130)},
    ]
    qa_sheet = crop_set(QA_DIR, qa_specs, "bulk_surface_hard_qa_crop_contact_sheet.png")
    strict_sheet = crop_set(STRICT_DIR, strict_specs, "bulk_surface_hard_STRICT_QA_contact_sheet.png")
    print(json.dumps({"qa_contact_sheet": str(qa_sheet.relative_to(ROOT)), "strict_contact_sheet": str(strict_sheet.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()
