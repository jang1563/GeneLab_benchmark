#!/usr/bin/env python3
"""Build strict visual QA crops for the NES conservation asset."""

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
    / "nes_conservation"
)
SOURCE = ASSET_DIR / "nes_conservation_predicts_transfer_premium.png"
QA_DIR = ASSET_DIR / "qa_crops_v2"
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
        image = Image.open(directory / f"nes_conservation_{name}.png")
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
        out = target_dir / f"nes_conservation_{name}.png"
        source.crop(box).save(out, quality=95)
        manifest.append({"name": name, "box": box, "path": str(out.relative_to(ROOT))})
    sheet = write_contact_sheet(specs, target_dir, sheet_name)
    (target_dir / "nes_conservation_crop_manifest.json").write_text(
        json.dumps({"source": str(SOURCE.relative_to(ROOT)), "crops": manifest, "contact_sheet": str(sheet.relative_to(ROOT))}, indent=2) + "\n"
    )
    return sheet


def main() -> None:
    qa_specs = [
        {"name": "01_header_badges", "box": (120, 40, 3720, 260)},
        {"name": "02_scatter_panel", "box": (130, 330, 2315, 1820)},
        {"name": "03_rank_panel", "box": (2360, 330, 3715, 1820)},
        {"name": "04_footer", "box": (130, 1860, 3710, 2065)},
        {"name": "05_full_downscaled", "box": (0, 0, 3840, 2160)},
    ]
    strict_specs = [
        {"name": "01_full_thumbnail", "box": (0, 0, 3840, 2160)},
        {"name": "02_header_title_badges", "box": (120, 40, 3720, 260)},
        {"name": "03_scatter_top_labels", "box": (260, 500, 2245, 900)},
        {"name": "04_scatter_mid_labels", "box": (260, 820, 1730, 1280)},
        {"name": "05_scatter_low_labels", "box": (260, 1240, 1160, 1690)},
        {"name": "06_axis_labels", "box": (150, 1570, 2260, 1800)},
        {"name": "07_rank_score_cards", "box": (2410, 520, 3680, 860)},
        {"name": "08_rank_bridge", "box": (2410, 890, 3680, 1375)},
        {"name": "09_scope_card", "box": (2410, 1360, 3680, 1665)},
        {"name": "10_footer", "box": (150, 1880, 3690, 2045)},
    ]
    qa_sheet = crop_set(QA_DIR, qa_specs, "nes_conservation_qa_crop_contact_sheet.png")
    strict_sheet = crop_set(STRICT_DIR, strict_specs, "nes_conservation_STRICT_QA_contact_sheet.png")
    print(json.dumps({"qa_contact_sheet": str(qa_sheet.relative_to(ROOT)), "strict_contact_sheet": str(strict_sheet.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()
