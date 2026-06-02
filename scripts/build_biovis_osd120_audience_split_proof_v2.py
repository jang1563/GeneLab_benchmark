#!/usr/bin/env python3
"""Build audience-facing OSD-120 split proof visuals.

v0.2 keeps the same manifest-derived sample counts as the earlier proof crop,
but lowers internal terminology for presentation use.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont, ImageOps


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = ROOT / "v9" / "multispecies" / "interaction_task_manifests" / "draft_osd120_arabidopsis_root_light_interaction_spaceflight.json"
CAPTURE_PATH = ROOT / "output" / "biovis_osdr_source_record_captures_v0_1" / "osd_120_viewport.png"
OUT = ROOT / "output" / "biovis_osd120_audience_split_proof_v0_2"
PROOF = OUT / "proof"
PANELS = OUT / "panels"
QA = OUT / "qa"
CREATED = "2026-06-02"

COLORS = {
    "void": "#080B0F",
    "deep": "#0C1218",
    "deep2": "#101922",
    "deep3": "#172231",
    "panel": "#13202B",
    "ink": "#F4F7F8",
    "soft": "#B7C5D1",
    "muted": "#718293",
    "rule": "#33465A",
    "blue": "#2D6F9F",
    "sky": "#6BAED6",
    "teal": "#1AA090",
    "green": "#178B63",
    "amber": "#B4832F",
    "red": "#B23A3A",
    "violet": "#6750A4",
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
    "eyebrow": font(20, bold=True),
    "title": font(48, bold=True),
    "subtitle": font(25),
    "h1": font(34, bold=True),
    "h2": font(26, bold=True),
    "body": font(20),
    "small": font(16),
    "tiny": font(13),
    "slide_title": font(66, bold=True),
    "slide_subtitle": font(32),
}


VISIBLE_TEXT = [
    "OSD-120 same-study task check",
    "Arabidopsis root RNA-seq",
    "Hold out one genotype",
    "Training samples",
    "Held-out samples",
    "Same OSDR study",
    "Not mission-held-out",
    "Not cross-species",
]


def ensure() -> None:
    for directory in [OUT, PROOF, PANELS, QA]:
        directory.mkdir(parents=True, exist_ok=True)


def load_manifest() -> dict:
    with MANIFEST_PATH.open() as handle:
        return json.load(handle)


def fit_cover_top(image: Image.Image, size: tuple[int, int]) -> Image.Image:
    return ImageOps.fit(image.convert("RGB"), size, method=Image.Resampling.LANCZOS, centering=(0.5, 0.0))


def fit_contain(image: Image.Image, size: tuple[int, int]) -> Image.Image:
    copy = image.convert("RGB")
    copy.thumbnail(size, Image.Resampling.LANCZOS)
    return copy


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


def card(draw: ImageDraw.ImageDraw, xywh: tuple[int, int, int, int], *, alpha: int = 226, outline: str = "rule", radius: int = 24) -> None:
    x, y, w, h = xywh
    draw.rounded_rectangle((x, y, x + w, y + h), radius=radius, fill=rgba("deep2", alpha), outline=rgba(outline, 230), width=2)


def pill(draw: ImageDraw.ImageDraw, xy: tuple[int, int], text: str, tone: str, *, fill_token: str = "deep2") -> None:
    x, y = xy
    width = max(112, int(draw.textlength(text, font=F["small"])) + 42)
    draw.rounded_rectangle((x, y, x + width, y + 34), radius=14, fill=rgba(fill_token, 236), outline=rgb(tone), width=2)
    draw.ellipse((x + 14, y + 12, x + 24, y + 22), fill=rgb(tone))
    draw.text((x + 34, y + 9), text, font=F["small"], fill=rgb("ink"))


def draw_background(canvas: Image.Image, *, accent: str = "green") -> None:
    draw = ImageDraw.Draw(canvas)
    w, h = canvas.size
    draw.rectangle((0, 0, w, h), fill=rgb("void"))
    for y in range(0, h, 2):
        t = y / max(1, h - 1)
        color = tuple(int(rgb("void")[i] * (1 - t) + rgb("deep2")[i] * t) for i in range(3))
        draw.line((0, y, w, y), fill=color + (255,), width=2)
    for x in range(160, w, 160):
        draw.line((x, 0, x, h), fill=rgba("rule", 36), width=1)
    for y in range(140, h, 140):
        draw.line((0, y, w, y), fill=rgba("rule", 30), width=1)
    center = (int(w * 0.65), int(h * 0.34))
    for idx, radius in enumerate([700, 920, 1140]):
        bbox = (center[0] - radius, center[1] - int(radius * 0.36), center[0] + radius, center[1] + int(radius * 0.36))
        draw.arc(bbox, 198, 348, fill=rgba(accent, 58 - idx * 12), width=3)
    for idx in range(36):
        x = int((idx * 257) % w)
        y = int((idx * 421) % h)
        r = 2 + idx % 3
        draw.ellipse((x - r, y - r, x + r, y + r), fill=rgba("soft", 36))


def draw_root_motif(draw: ImageDraw.ImageDraw, xy: tuple[int, int], *, scale: float = 1.0) -> None:
    x, y = xy
    draw.line((x, y, x, y + int(128 * scale)), fill=rgba("green", 130), width=max(2, int(5 * scale)))
    for idx, (dy, direction) in enumerate([(30, -1), (54, 1), (84, -1), (108, 1)]):
        x2 = x + direction * int((42 + idx * 9) * scale)
        y2 = y + int((dy + 22) * scale)
        draw.line((x, y + int(dy * scale), x2, y2), fill=rgba("green", 120), width=max(1, int(3 * scale)))
    draw.ellipse((x - int(34 * scale), y - int(42 * scale), x + int(6 * scale), y + int(2 * scale)), outline=rgba("green", 115), width=max(2, int(3 * scale)))
    draw.ellipse((x - int(2 * scale), y - int(44 * scale), x + int(48 * scale), y + int(2 * scale)), outline=rgba("green", 115), width=max(2, int(3 * scale)))


def label_distribution_text(dist: dict) -> str:
    ground = dist.get("Ground", 0)
    leo = dist.get("LEO_or_ISS", 0)
    return f"Ground {ground} / ISS {leo}"


def clean_value(value: str) -> str:
    return value.replace(".Treatment", "").replace(".ecotype", "").replace(".PhyD", " PhyD").replace(".", " ")


def draw_split_lane(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    fold: dict,
    *,
    width: int,
    tone: str,
    compact: bool = False,
) -> None:
    x, y = xy
    train_w = int(width * 0.55)
    test_w = int(width * 0.30)
    draw.text((x, y), clean_value(str(fold["heldout_value"])), font=F["h2"], fill=rgb("ink"))
    draw.text((x, y + 32), "held-out group", font=F["small"], fill=rgb("soft"))
    y0 = y + 68
    draw.rounded_rectangle((x, y0, x + train_w, y0 + 76), radius=16, fill=rgba("deep3", 238), outline=rgb("blue"), width=2)
    draw.rounded_rectangle((x + train_w + 96, y0, x + train_w + 96 + test_w, y0 + 76), radius=16, fill=rgba("deep3", 238), outline=rgb(tone), width=2)
    draw.text((x + 28, y0 + 15), f"TRAIN {fold['n_train']}", font=F["h2"], fill=rgb("ink"))
    draw.text((x + train_w + 124, y0 + 15), f"TEST {fold['n_test']}", font=F["h2"], fill=rgb("ink"))
    dash_y = y0 + 38
    for xx in range(x + train_w + 22, x + train_w + 78, 14):
        draw.line((xx, dash_y, xx + 8, dash_y), fill=rgb("amber"), width=3)
    draw.text((x + train_w + 22, y0 + 8), "held out", font=F["tiny"], fill=rgb("amber"))
    if not compact:
        draw.text((x + 28, y0 + 104), f"train balance: {label_distribution_text(fold['train_label_distribution'])}", font=F["small"], fill=rgb("soft"))
        draw.text((x + train_w + 124, y0 + 104), f"test balance: {label_distribution_text(fold['test_label_distribution'])}", font=F["small"], fill=rgb("soft"))


def render_clean_split_proof(manifest: dict) -> Path:
    split = manifest["split"]
    primary = split["primary_candidate_folds"]
    secondary = split["secondary_light_treatment_folds"]
    canvas = Image.new("RGBA", (2200, 1250), rgba("deep", 255))
    draw = ImageDraw.Draw(canvas)
    draw_root_motif(draw, (1980, 92), scale=1.1)
    draw.text((70, 54), "OSD-120 same-study task check", font=F["title"], fill=rgb("ink"))
    draw_wrapped(
        draw,
        (70, 110),
        "Train on Arabidopsis root RNA-seq samples, then test on a held-out genotype or light condition.",
        F["subtitle"],
        rgb("soft"),
        width=1560,
        line_height=32,
        max_lines=2,
    )
    draw.line((70, 172, 2130, 172), fill=rgba("rule", 255), width=2)
    pill(draw, (70, 202), "OSD-120", "green")
    pill(draw, (190, 202), "36 samples", "sky")
    pill(draw, (334, 202), "same OSDR study", "amber")
    pill(draw, (548, 202), "not mission-held-out", "red")
    pill(draw, (786, 202), "not cross-species", "red")

    card(draw, (70, 286, 1280, 762), alpha=216, outline="rule")
    draw.text((110, 326), "Primary split: hold out one genotype", font=F["h1"], fill=rgb("ink"))
    draw.text((110, 366), "Each test group keeps Ground and ISS labels balanced.", font=F["body"], fill=rgb("soft"))
    y = 440
    for fold in primary:
        draw_split_lane(draw, (110, y), fold, width=1040, tone="green")
        y += 190

    card(draw, (1410, 286, 720, 762), alpha=216, outline="rule")
    draw.text((1450, 326), "Secondary check", font=F["h1"], fill=rgb("ink"))
    draw.text((1450, 366), "Hold out light condition.", font=F["body"], fill=rgb("soft"))
    y = 446
    for fold in secondary:
        draw_split_lane(draw, (1450, y), fold, width=580, tone="amber", compact=True)
        y += 236

    card(draw, (70, 1084, 2060, 96), alpha=190, outline="rule", radius=18)
    boundary_label = "Claim boundary"
    draw.text((104, 1118), boundary_label, font=F["h2"], fill=rgb("ink"))
    body_x = 104 + int(draw.textlength(boundary_label, font=F["h2"])) + 34
    draw.text((body_x, 1122), "This shows a same-study task check for OSD-120. It is not evidence of mission-held-out or cross-species generalization.", font=F["body"], fill=rgb("soft"))
    out = PROOF / "osd120_same_study_task_check_proof.png"
    canvas.convert("RGB").save(out)
    return out


def crop_osdr_source() -> Image.Image:
    source = Image.open(CAPTURE_PATH).convert("RGB")
    return source.crop((292, 145, 1580, 1120))


def shadowed_surface(
    canvas: Image.Image,
    image: Image.Image,
    xy: tuple[int, int],
    *,
    width: int,
    height: int,
    radius: int,
    outline: tuple[int, int, int, int],
    contain: bool = False,
) -> None:
    x, y = xy
    surface = (fit_contain(image, (width - 70, height - 70)) if contain else fit_cover_top(image, (width, height))).convert("RGBA")
    shadow = Image.new("RGBA", canvas.size, (0, 0, 0, 0))
    sd = ImageDraw.Draw(shadow)
    sd.rounded_rectangle((x + 18, y + 26, x + width + 18, y + height + 26), radius=radius, fill=(0, 0, 0, 155))
    shadow = shadow.filter(ImageFilter.GaussianBlur(34))
    canvas.alpha_composite(shadow)
    draw = ImageDraw.Draw(canvas)
    draw.rounded_rectangle((x, y, x + width, y + height), radius=radius, fill=rgba("deep2", 238), outline=outline, width=3)
    if contain:
        px = x + (width - surface.width) // 2
        py = y + (height - surface.height) // 2
        canvas.alpha_composite(surface, (px, py))
    else:
        mask = Image.new("L", (width, height), 0)
        md = ImageDraw.Draw(mask)
        md.rounded_rectangle((0, 0, width, height), radius=radius, fill=255)
        canvas.paste(surface, (x, y), mask)
        draw.rounded_rectangle((x, y, x + width, y + height), radius=radius, outline=outline, width=3)


def arrow(draw: ImageDraw.ImageDraw, start: tuple[int, int], end: tuple[int, int], *, tone: str = "green") -> None:
    x1, y1 = start
    x2, y2 = end
    mid = ((x1 + x2) // 2, min(y1, y2) - 80)
    points = []
    for i in range(44):
        t = i / 43
        x = int((1 - t) ** 2 * x1 + 2 * (1 - t) * t * mid[0] + t**2 * x2)
        y = int((1 - t) ** 2 * y1 + 2 * (1 - t) * t * mid[1] + t**2 * y2)
        points.append((x, y))
    draw.line(points, fill=rgba(tone, 205), width=5)
    angle = math.atan2(y2 - points[-3][1], x2 - points[-3][0])
    size = 23
    p1 = (x2, y2)
    p2 = (int(x2 - size * math.cos(angle - 0.45)), int(y2 - size * math.sin(angle - 0.45)))
    p3 = (int(x2 - size * math.cos(angle + 0.45)), int(y2 - size * math.sin(angle + 0.45)))
    draw.polygon([p1, p2, p3], fill=rgba(tone, 220))


def caption(draw: ImageDraw.ImageDraw, xywh: tuple[int, int, int, int], title: str, body: str, *, tone: str) -> None:
    x, y, w, h = xywh
    draw.rounded_rectangle((x, y, x + w, y + h), radius=24, fill=rgba("deep2", 228), outline=rgba(tone, 160), width=2)
    draw.text((x + 30, y + 28), title, font=F["h2"], fill=rgb("ink"))
    draw_wrapped(draw, (x + 30, y + 76), body, F["small"], rgb("soft"), width=w - 60, line_height=27, max_lines=4)


def render_source_to_task_panel(proof_path: Path) -> Path:
    canvas = Image.new("RGBA", (3840, 2160), rgba("void", 255))
    draw_background(canvas, accent="green")
    draw = ImageDraw.Draw(canvas)
    draw_root_motif(draw, (3360, 250), scale=2.7)
    draw.text((120, 92), "OSD-120 source to same-study task", font=F["slide_title"], fill=rgb("ink"))
    draw_wrapped(
        draw,
        (122, 176),
        "Official Arabidopsis source page on the left; manifest-derived train/test task check on the right.",
        F["slide_subtitle"],
        rgb("soft"),
        width=2280,
        line_height=42,
        max_lines=2,
    )
    pill(draw, (122, 260), "official source page", "green")
    pill(draw, (376, 260), "local task proof", "sky")
    pill(draw, (608, 260), "same-study check", "amber")
    pill(draw, (842, 260), "claim boundary visible", "red")

    source_crop = crop_osdr_source()
    proof = Image.open(proof_path).convert("RGB")
    draw.text((120, 350), "Official OSDR source page", font=F["h2"], fill=rgb("ink"))
    draw.text((1780, 350), "Audience-facing local proof", font=F["h2"], fill=rgb("ink"))
    shadowed_surface(canvas, source_crop, (120, 405), width=1460, height=980, radius=34, outline=rgba("green", 230))
    shadowed_surface(canvas, proof, (1780, 405), width=1900, height=1080, radius=38, outline=rgba("amber", 230), contain=True)
    arrow(draw, (1580, 890), (1780, 930), tone="green")

    caption(
        draw,
        (120, 1524, 1460, 238),
        "Viewer path",
        "Start at the official OSDR study record, then move to the local split object. The labels avoid internal shorthand so non-specialists can follow the method.",
        tone="green",
    )
    caption(
        draw,
        (1780, 1524, 1420, 238),
        "What this proves",
        "The local task has sample-count-backed train/test folds for one OSDR study. It does not prove mission-held-out or cross-species generalization.",
        tone="amber",
    )
    draw.rounded_rectangle((120, 1948, 3720, 2040), radius=22, fill=rgba("deep2", 230), outline=rgba("rule", 220), width=2)
    draw.text((158, 1981), "Claim boundary", font=F["h2"], fill=rgb("ink"))
    draw_wrapped(
        draw,
        (410, 1984),
        "Use this as a methods bridge, not as a final scored result figure. Keep the OSDR source URL and draft status in the deck footer.",
        F["small"],
        rgb("soft"),
        width=3140,
        line_height=27,
        max_lines=2,
    )
    out = PANELS / "01_dark_osd120_clean_source_to_task.png"
    canvas.convert("RGB").save(out)
    return out


def grayscale_exports(paths: list[Path]) -> list[Path]:
    outs = []
    for path in paths:
        out = QA / f"{path.stem}_grayscale.png"
        Image.open(path).convert("L").save(out)
        outs.append(out)
    return outs


def luminance(color: tuple[int, int, int]) -> float:
    values = [v / 255 for v in color]
    linear = [v / 12.92 if v <= 0.03928 else ((v + 0.055) / 1.055) ** 2.4 for v in values]
    return 0.2126 * linear[0] + 0.7152 * linear[1] + 0.0722 * linear[2]


def contrast(c1: tuple[int, int, int], c2: tuple[int, int, int]) -> float:
    l1 = luminance(c1)
    l2 = luminance(c2)
    return (max(l1, l2) + 0.05) / (min(l1, l2) + 0.05)


def write_metadata(proof_path: Path, panel_path: Path, grayscale: list[Path], manifest: dict) -> None:
    blocked_terms = ["within-source", "hidden", "interaction split"]
    visible = " ".join(VISIBLE_TEXT).lower()
    data = {
        "created": CREATED,
        "source_manifest": rel(MANIFEST_PATH),
        "source_record": "https://osdr.nasa.gov/bio/repo/data/studies/OSD-120",
        "outputs": {
            "proof": rel(proof_path),
            "panel": rel(panel_path),
            "grayscale_qa": [rel(path) for path in grayscale],
        },
        "sample_counts": {
            "n_samples": manifest["split"]["n_samples"],
            "label_distribution": manifest["split"]["label_distribution"],
            "genotype_or_ecotype_distribution": manifest["split"]["genotype_or_ecotype_distribution"],
            "light_treatment_distribution": manifest["split"]["light_treatment_distribution"],
        },
        "visible_text_policy": {
            "preferred_terms": VISIBLE_TEXT,
            "blocked_terms_absent_from_new_visible_text": {term: term.lower() not in visible for term in blocked_terms},
        },
        "claim_boundary": "Same-study OSD-120 task check only; not mission-held-out and not cross-species generalization.",
    }
    (OUT / "osd120_audience_split_proof_manifest.json").write_text(json.dumps(data, indent=2) + "\n")
    qa = {
        "created": CREATED,
        "contrast_checks": {
            "ink_on_deep": round(contrast(rgb("ink"), rgb("deep")), 2),
            "soft_on_deep": round(contrast(rgb("soft"), rgb("deep")), 2),
            "ink_on_void": round(contrast(rgb("ink"), rgb("void")), 2),
        },
        "visual_checks": [
            "Audience-facing proof removes internal terms from visible labels.",
            "Train/test separation remains sample-count-backed from the task manifest.",
            "Claim boundary is visible in both the proof crop and slide-scale panel.",
            "Color is reinforced by text labels and grayscale QA outputs.",
        ],
        "review_required": [
            "Final deck should use editable text for callouts, with this PNG as the evidence surface.",
            "Do not use as final result or generalization evidence.",
        ],
    }
    (OUT / "osd120_audience_split_proof_qa.json").write_text(json.dumps(qa, indent=2) + "\n")


def main() -> None:
    ensure()
    manifest = load_manifest()
    proof_path = render_clean_split_proof(manifest)
    panel_path = render_source_to_task_panel(proof_path)
    grayscale = grayscale_exports([proof_path, panel_path])
    write_metadata(proof_path, panel_path, grayscale, manifest)
    print(
        json.dumps(
            {
                "output": rel(OUT),
                "proof": rel(proof_path),
                "panel": rel(panel_path),
                "grayscale_qa": [rel(path) for path in grayscale],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
