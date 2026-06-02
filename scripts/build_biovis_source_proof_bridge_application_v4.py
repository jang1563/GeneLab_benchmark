#!/usr/bin/env python3
"""Stress-test v0.4 source-proof modules inside bridge-style layouts."""

from __future__ import annotations

import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "biovis_source_proof_bridge_application_v0_4"
PROOF = ROOT / "assets" / "biovis_source_proof_module_pack_v0_4"
SYMBOL = ROOT / "assets" / "biovis_symbol_module_pack_v0_3"
CREATED = "2026-06-02"
W, H = 3840, 2160

COLORS = {
    "paper": "#F7F4EC",
    "paper2": "#FBFAF7",
    "ink": "#17212B",
    "muted": "#5D6978",
    "rule": "#AEB8C5",
    "blue": "#2D6F9F",
    "sky": "#6BAED6",
    "green": "#178B63",
    "teal": "#1AA090",
    "amber": "#B4832F",
    "red": "#B23A3A",
    "deep": "#0D1720",
    "deep2": "#111D27",
    "deep3": "#172636",
    "white": "#F2F6F8",
    "soft": "#B8C4CF",
}


def rgb(token: str) -> tuple[int, int, int]:
    value = COLORS[token].lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


def rgba(token: str, alpha: int) -> tuple[int, int, int, int]:
    return rgb(token) + (alpha,)


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
    "eyebrow": font(28, bold=True),
    "title": font(72, bold=True),
    "subtitle": font(34),
    "h2": font(38, bold=True),
    "h3": font(28, bold=True),
    "body": font(26),
    "small": font(19),
    "tiny": font(15),
}


def ensure() -> None:
    OUT.mkdir(parents=True, exist_ok=True)


def draw_text(draw: ImageDraw.ImageDraw, xy: tuple[int, int], text: str, style: str, color: str) -> None:
    draw.text(xy, text, font=F[style], fill=rgb(color))


def paste_fit(canvas: Image.Image, path: Path, box: tuple[int, int, int, int], *, opacity: float = 1.0) -> None:
    image = Image.open(path).convert("RGBA")
    x, y, w, h = box
    scale = min(w / image.width, h / image.height)
    image = image.resize((max(1, int(image.width * scale)), max(1, int(image.height * scale))), Image.Resampling.LANCZOS)
    if opacity < 1.0:
        alpha = image.getchannel("A")
        alpha = alpha.point(lambda px: int(px * opacity))
        image.putalpha(alpha)
    canvas.alpha_composite(image, (x + (w - image.width) // 2, y + (h - image.height) // 2))


def card(draw: ImageDraw.ImageDraw, xywh: tuple[int, int, int, int], *, dark: bool = False) -> None:
    x, y, w, h = xywh
    fill = "deep2" if dark else "paper2"
    outline = "deep3" if dark else "rule"
    draw.rounded_rectangle((x, y, x + w, y + h), radius=22, fill=rgba(fill, 238), outline=rgba(outline, 210), width=2)


def badge(draw: ImageDraw.ImageDraw, xy: tuple[int, int], text: str, tone: str, *, dark: bool = False) -> None:
    x, y = xy
    width = max(130, 48 + int(draw.textlength(text, font=F["small"])))
    fill = "deep2" if dark else "paper2"
    ink = "white" if dark else "ink"
    draw.rounded_rectangle((x, y, x + width, y + 44), radius=18, fill=rgba(fill, 235), outline=rgb(tone), width=3)
    draw.ellipse((x + 18, y + 15, x + 30, y + 27), fill=rgb(tone))
    draw.text((x + 42, y + 12), text, font=F["small"], fill=rgb(ink))


def proof_ready_dark() -> Image.Image:
    canvas = Image.new("RGBA", (W, H), rgb("deep") + (255,))
    draw = ImageDraw.Draw(canvas)
    for x in [420, 1920, 3320]:
        draw.line((x, 160, x, H - 170), fill=rgba("deep3", 180), width=2)
    for y in [310, 1050, 1810]:
        draw.line((160, y, W - 160, y), fill=rgba("deep3", 210), width=2)

    draw_text(draw, (160, 108), "SOURCE-PROOF PLACEHOLDER TEST", "eyebrow", "sky")
    draw_text(draw, (160, 164), "Design the evidence frame before the final figure exists", "title", "white")
    draw_text(draw, (160, 256), "v0.4 modules reserve proof object, status, source line, and caveat", "subtitle", "soft")

    paste_fit(canvas, SYMBOL / "modules_dark" / "png" / "space_bio_assay_lane.png", (160, 374, 2020, 515))
    card(draw, (2260, 360, 1480, 570), dark=True)
    paste_fit(canvas, PROOF / "modules_dark" / "png" / "source_dataset_record_plate.png", (2295, 390, 700, 390))
    paste_fit(canvas, PROOF / "modules_dark" / "png" / "expression_matrix_proof_plate.png", (3020, 390, 700, 390))
    badge(draw, (2295, 810), "source slot", "blue", dark=True)
    badge(draw, (2470, 810), "processed", "teal", dark=True)
    badge(draw, (2645, 810), "caveat", "red", dark=True)

    paste_fit(canvas, PROOF / "modules_dark" / "png" / "held_out_task_proof_plate.png", (160, 1110, 1040, 610))
    paste_fit(canvas, PROOF / "modules_dark" / "png" / "single_cell_embedding_proof_plate.png", (1265, 1110, 1040, 610))
    paste_fit(canvas, PROOF / "modules_dark" / "png" / "result_claim_plate.png", (2370, 1110, 1040, 610))

    draw_text(draw, (160, 1860), "Verdict", "h2", "white")
    notes = [
        ("Pass: evidence slot is explicit", "final screenshot/plot can replace placeholder without redesign", "green"),
        ("Pass: status travels with proof", "source, processed, held-out, caveat, and validation states remain visible", "blue"),
        ("Caution: one major proof module per slide", "multiple large modules should be appendix or review-board material", "amber"),
    ]
    x = 160
    for head, body, tone in notes:
        card(draw, (x, 1930, 1080, 112), dark=True)
        draw.ellipse((x + 28, 1968, x + 52, 1992), fill=rgb(tone))
        draw_text(draw, (x + 74, 1948), head, "h3", "white")
        draw_text(draw, (x + 74, 1986), body, "small", "soft")
        x += 1160
    return canvas.convert("RGB")


def proof_ready_light() -> Image.Image:
    canvas = Image.new("RGBA", (W, H), rgb("paper") + (255,))
    draw = ImageDraw.Draw(canvas)
    draw_text(draw, (160, 108), "SOURCE-PROOF PLACEHOLDER TEST", "eyebrow", "blue")
    draw_text(draw, (160, 164), "A final slide should have proof slots, not just motifs", "title", "ink")
    draw_text(draw, (160, 256), "Light-canvas test: source record, matrix proof, held-out split, and claim footer", "subtitle", "muted")
    draw.line((160, 320, W - 160, 320), fill=rgba("rule", 210), width=2)

    paste_fit(canvas, PROOF / "modules" / "png" / "source_dataset_record_plate.png", (150, 395, 1100, 640))
    paste_fit(canvas, PROOF / "modules" / "png" / "expression_matrix_proof_plate.png", (1350, 395, 1100, 640))
    paste_fit(canvas, PROOF / "modules" / "png" / "held_out_task_proof_plate.png", (2550, 395, 1100, 640))
    paste_fit(canvas, PROOF / "modules" / "png" / "result_claim_plate.png", (150, 1190, 1180, 680))
    paste_fit(canvas, SYMBOL / "modules" / "png" / "claim_boundary_bar.png", (1500, 1330, 1460, 360))

    card(draw, (3080, 1265, 510, 390), dark=False)
    draw_text(draw, (3115, 1315), "Production use", "h2", "ink")
    bullets = [
        "1. Replace proof slot first",
        "2. Keep status badges attached",
        "3. Add source/date footer",
        "4. Then write the claim",
    ]
    y = 1375
    for item in bullets:
        draw_text(draw, (3115, y), item, "body", "muted")
        y += 58

    draw.line((160, 1968, W - 160, 1968), fill=rgba("rule", 190), width=2)
    draw_text(draw, (160, 2028), "Design rule: proof placeholders prevent evidence from being retrofitted after the slide is already composed.", "body", "muted")
    return canvas.convert("RGB")


def contact_sheet(light: Image.Image, dark: Image.Image) -> Image.Image:
    canvas = Image.new("RGB", (2400, 1500), rgb("paper"))
    draw = ImageDraw.Draw(canvas)
    draw_text(draw, (72, 48), "v0.4 source-proof bridge application", "h2", "ink")
    draw_text(draw, (72, 96), "Stress test: evidence-ready modules replacing empty bridge space.", "body", "muted")
    thumb_w, thumb_h = 1040, 585
    canvas.paste(light.resize((thumb_w, thumb_h), Image.Resampling.LANCZOS), (72, 190))
    canvas.paste(dark.resize((thumb_w, thumb_h), Image.Resampling.LANCZOS), (72 + thumb_w + 120, 190))
    draw_text(draw, (72, 832), "01 light proof-ready bridge", "h3", "ink")
    draw_text(draw, (72 + thumb_w + 120, 832), "02 dark proof-ready bridge", "h3", "ink")
    notes = [
        ("Best use", "one proof module plus one method module on a main slide"),
        ("Appendix use", "multiple proof modules are useful for review boards and QA packs"),
        ("Claim safety", "source and caveat state travel with the proof object"),
        ("Next gap", "replace placeholders with actual OSDR/GEO/result screenshots"),
    ]
    for idx, (head, body) in enumerate(notes):
        x = 72 + (idx % 2) * 1120
        y = 940 + (idx // 2) * 190
        card(draw, (x, y, 1020, 132), dark=False)
        draw_text(draw, (x + 30, y + 30), head, "h3", "ink")
        draw_text(draw, (x + 30, y + 76), body, "small", "muted")
    return canvas


def write_manifest() -> None:
    manifest = {
        "created": CREATED,
        "purpose": "Apply v0.4 source-proof placeholder modules inside bridge-style layouts.",
        "inputs": {
            "source_proof_pack": "assets/biovis_source_proof_module_pack_v0_4",
            "symbol_module_pack": "assets/biovis_symbol_module_pack_v0_3",
        },
        "outputs": {
            "light": "output/biovis_source_proof_bridge_application_v0_4/01_light_source_proof_bridge.png",
            "dark": "output/biovis_source_proof_bridge_application_v0_4/02_dark_source_proof_bridge.png",
            "contact_sheet": "output/biovis_source_proof_bridge_application_v0_4/03_source_proof_bridge_contact_sheet.png",
        },
    }
    (OUT / "biovis_source_proof_bridge_application_v0_4_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    qa = {
        "created": CREATED,
        "automatic_checks": {
            "light_png_exists": (OUT / "01_light_source_proof_bridge.png").exists(),
            "dark_png_exists": (OUT / "02_dark_source_proof_bridge.png").exists(),
            "contact_sheet_exists": (OUT / "03_source_proof_bridge_contact_sheet.png").exists(),
        },
        "manual_review": {
            "initial_verdict": "pass with conditions",
            "light_bridge": "conditional pass: clear source-proof logic, but density is closer to appendix or review-board material than a final main slide",
            "dark_bridge": "pass: strongest premium direction; proof objects, status badges, and caveat/status surfaces stay readable on a layered dark field",
            "proof_slot_clarity": "pass: real source screenshots and result plots can replace placeholders without redesigning the slide",
            "claim_boundary": "pass: source, processed, held-out, caveat, and validation states travel with the proof object",
            "main_slide_rule": "use one dominant proof module plus one method or claim-boundary module for a final main-deck slide",
            "next_step": "replace placeholder surfaces with actual OSDR/GeneLab/GEO records and real analysis outputs, then re-test a production slide",
        },
    }
    (OUT / "biovis_source_proof_bridge_application_v0_4_qa.json").write_text(json.dumps(qa, indent=2) + "\n")


def main() -> None:
    ensure()
    dark = proof_ready_dark()
    light = proof_ready_light()
    sheet = contact_sheet(light, dark)
    light.save(OUT / "01_light_source_proof_bridge.png")
    dark.save(OUT / "02_dark_source_proof_bridge.png")
    sheet.save(OUT / "03_source_proof_bridge_contact_sheet.png")
    write_manifest()
    print(json.dumps({"output": str(OUT.relative_to(ROOT)), "slides": 2, "contact_sheet": True}, indent=2))


if __name__ == "__main__":
    main()
