#!/usr/bin/env python3
"""Build visual QA crops for slide 46 DECOMPOSE interactions asset."""

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
    / "decompose_interactions"
)
SOURCE = ASSET_DIR / "stressor_interactions_dominate_combined_effects_premium.png"
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
        image = Image.open(directory / f"decompose_interactions_{name}.png")
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
        out = target_dir / f"decompose_interactions_{name}.png"
        source.crop(box).save(out, quality=95)
        manifest.append({"name": name, "box": box, "path": str(out.relative_to(ROOT))})
    sheet = write_contact_sheet(specs, target_dir, sheet_name)
    (target_dir / "decompose_interactions_crop_manifest.json").write_text(
        json.dumps({"source": str(SOURCE.relative_to(ROOT)), "crops": manifest, "contact_sheet": str(sheet.relative_to(ROOT))}, indent=2) + "\n",
        encoding="utf-8",
    )
    return sheet


def main() -> None:
    qa_specs = [
        {"name": "01_header_badges", "box": (100, 40, 3740, 235)},
        {"name": "02_title_subtitle", "box": (120, 310, 3400, 520)},
        {"name": "03_method_panel", "box": (120, 610, 1210, 1462)},
        {"name": "04_surface_panel", "box": (1265, 610, 2575, 1462)},
        {"name": "05_icp_panel", "box": (2630, 610, 3720, 1462)},
        {"name": "06_decompose_flow_strip", "box": (120, 1530, 3720, 1848)},
        {"name": "07_footer_scope", "box": (120, 1910, 3720, 2076)},
        {"name": "08_full_downscaled", "box": (0, 0, 3840, 2160)},
    ]
    strict_specs = [
        {"name": "01_full_thumbnail", "box": (0, 0, 3840, 2160)},
        {"name": "02_header_title_badges", "box": (100, 50, 3740, 520)},
        {"name": "03_method_factorial_grid", "box": (120, 720, 1210, 1170)},
        {"name": "04_method_formula_rule", "box": (120, 1160, 1210, 1462)},
        {"name": "05_surface_bars", "box": (1265, 830, 2575, 1295)},
        {"name": "06_surface_callout", "box": (1265, 1268, 2575, 1462)},
        {"name": "07_icp_bars", "box": (2630, 760, 3720, 1255)},
        {"name": "08_icp_lower_cards", "box": (2630, 1240, 3720, 1462)},
        {"name": "09_flow_nodes", "box": (120, 1530, 3720, 1848)},
        {"name": "10_footer", "box": (120, 1910, 3720, 2076)},
    ]
    qa_sheet = crop_set(QA_DIR, qa_specs, "decompose_interactions_qa_crop_contact_sheet.png")
    strict_sheet = crop_set(STRICT_DIR, strict_specs, "decompose_interactions_STRICT_QA_contact_sheet.png")
    print(json.dumps({"qa_contact_sheet": str(qa_sheet.relative_to(ROOT)), "strict_contact_sheet": str(strict_sheet.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()
