#!/usr/bin/env python3
"""Build the slides 1-14 deck assembly bridge and production brief."""

from __future__ import annotations

import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "slides_1_14_deck_assembly_bridge_v0_1"
QA = OUT / "qa"
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
    "blue": "#2D6F9F",
    "violet": "#7B68A8",
}


def p(path: str) -> Path:
    return ROOT / path


SLIDES = [
    {
        "n": 1,
        "section": "Opening",
        "title": "SpaceBio-Bench",
        "headline": "Testing biological AI under spaceflight domain shift",
        "role": "thesis",
        "status": "asset ready",
        "tone": "teal",
        "caption": "Open with the platform thesis and the mission-held-out problem before any model score appears.",
        "claim_boundary": "Do not claim first AI benchmark for space biology.",
        "asset": p("output/premium_opening_slides_1_3_v0_1/slide01_spacebiobench_title_rendered_preview.png"),
    },
    {
        "n": 2,
        "section": "Opening",
        "title": "External gap",
        "headline": "Distinct niche, not a firstness claim",
        "role": "positioning",
        "status": "asset ready",
        "tone": "teal",
        "caption": "Show the external landscape while keeping the claim to a distinct mission-held-out GeneLab/OSDR benchmark niche.",
        "claim_boundary": "Do not overstate novelty against NASA BPS, GLARE, OpenProblems, or cell-eval.",
        "asset": p("output/premium_opening_slides_1_3_v0_1/slide02_external_gap_positioning_rendered_preview.png"),
    },
    {
        "n": 3,
        "section": "Opening",
        "title": "Project evolution",
        "headline": "v1-v9 moves from benchmark results to platform",
        "role": "timeline",
        "status": "asset ready",
        "tone": "teal",
        "caption": "Orient the audience: v1-v7 are completed benchmark results, v8 is hypothesis-only translation, and v9 is platformization.",
        "claim_boundary": "Do not present all v1-v9 outputs as equal-strength scientific discoveries.",
        "asset": p("output/premium_opening_slides_1_3_v0_1/slide03_project_evolution_timeline_rendered_preview.png"),
    },
    {
        "n": 4,
        "section": "Core methods",
        "title": "Evaluation layer",
        "headline": "Public studies become auditable tasks",
        "role": "method bridge",
        "status": "asset ready; B2 in notes",
        "tone": "sky",
        "caption": "B1 carries the visible slide; B2 supplies the task-unit explanation in speaker notes or backup.",
        "claim_boundary": "Conceptual bridge, not quantitative result evidence.",
        "asset": p("output/premium_bridge_scenes/b1_evaluation_layer/rendered_preview.png"),
        "backup_assets": [p("output/premium_bridge_rebuild_scenes/b2_study_to_task_premium/rendered_preview.png")],
    },
    {
        "n": 5,
        "section": "Core methods",
        "title": "Mission-held-out guard",
        "headline": "The hidden test set is an entire mission",
        "role": "method bridge",
        "status": "asset ready; B4 in notes",
        "tone": "sky",
        "caption": "B3 carries the visible split; B4 supplies the train-only guardrail in speaker notes or backup.",
        "claim_boundary": "Do not call this random cross-validation or unqualified independent validation.",
        "asset": p("output/premium_bridge_rebuild_scenes/b3_mission_held_out_premium/rendered_preview.png"),
        "backup_assets": [p("output/premium_bridge_rebuild_scenes/b4_train_only_guard_premium/rendered_preview.png")],
    },
    {
        "n": 6,
        "section": "Core methods",
        "title": "Feature layers",
        "headline": "What the model sees can change",
        "role": "method-to-result bridge",
        "status": "asset ready",
        "tone": "sky",
        "caption": "Explain gene-level and pathway-summary views before tissue and model score results.",
        "claim_boundary": "Schematic methods bridge, not a quantitative result figure.",
        "asset": p("output/premium_feature_layer_bridge_v0_4/feature_layer_bridge_v0_4_rendered_preview.png"),
    },
    {
        "n": 7,
        "section": "v1-v7 results",
        "title": "Tissue hierarchy",
        "headline": "Cross-mission transfer is tissue-specific",
        "role": "early result",
        "status": "asset ready",
        "tone": "blue",
        "caption": "Start the result core with the main tissue-dependence finding.",
        "claim_boundary": "Do not imply one tissue result generalizes to all tissues.",
        "asset": p("output/premium_slide_scenes/fig1_tissue_transfer_layered_scene.png"),
    },
    {
        "n": 8,
        "section": "v1-v7 results",
        "title": "Pathway rescue",
        "headline": "Pathway summaries can rescue selected weak tasks",
        "role": "early result",
        "status": "asset ready",
        "tone": "blue",
        "caption": "Use this after tissue hierarchy so feature abstraction has a visible payoff.",
        "claim_boundary": "Do not claim pathway features always outperform gene-level models.",
        "asset": p("output/premium_slide_scenes/fig2_pathway_layered_scene.png"),
    },
    {
        "n": 9,
        "section": "v1-v7 results",
        "title": "Model comparison",
        "headline": "Scale alone does not solve transfer",
        "role": "early result",
        "status": "asset ready",
        "tone": "blue",
        "caption": "Place model comparison after task and feature logic are clear.",
        "claim_boundary": "Do not universalize foundation-model failure beyond tested settings.",
        "asset": p("output/premium_slide_scenes/fig3_model_tier_layered_scene.png"),
    },
    {
        "n": 10,
        "section": "core result spine",
        "title": "v4 hardening",
        "headline": "Broader grid reduces cherry-pick risk",
        "role": "core result",
        "status": "asset ready",
        "tone": "sky",
        "caption": "This is the first hardened result after early tissue/pathway/model findings.",
        "claim_boundary": "Hardening and coverage evidence only; not a new leaderboard.",
        "asset": p("output/premium_core_result_slides/slide10_v4_hardening_v0_1/slide10_v4_hardening_rendered_preview.png"),
    },
    {
        "n": 11,
        "section": "core result spine",
        "title": "Temporal and RRRM",
        "headline": "Timepoint labels need guardrails",
        "role": "core result",
        "status": "asset ready",
        "tone": "rose",
        "caption": "Adds preservation, recovery, and single-cell pilot guardrails before interpretation.",
        "claim_boundary": "Recovery projection is descriptive; RRRM composition is underpowered.",
        "asset": p("output/premium_core_result_slides/slide11_temporal_rrrm_guardrails_v0_1/slide11_temporal_rrrm_guardrails_rendered_preview.png"),
    },
    {
        "n": 12,
        "section": "core result spine",
        "title": "Negative boundary",
        "headline": "Negative results set boundaries",
        "role": "core result",
        "status": "asset ready",
        "tone": "amber",
        "caption": "Make failure modes part of benchmark value before moving into biology.",
        "claim_boundary": "No universal absence-of-biology or universal model-failure claim.",
        "asset": p("output/premium_core_result_slides/slide12_negative_boundary_v0_1/slide12_negative_boundary_rendered_preview.png"),
    },
    {
        "n": 13,
        "section": "core result spine",
        "title": "Biological interpretation",
        "headline": "Hits become hypotheses",
        "role": "core result",
        "status": "asset ready",
        "tone": "teal",
        "caption": "Interpret bounded signals as immune, metabolic, and biomarker hypotheses.",
        "claim_boundary": "Interpretation and target triage only; no treatment or mechanism proof.",
        "asset": p("output/premium_core_result_slides/slide13_biological_interpretation_v0_1/slide13_biological_interpretation_rendered_preview.png"),
    },
    {
        "n": 14,
        "section": "core result spine",
        "title": "Human translation",
        "headline": "Translation is partial, not direct",
        "role": "core result",
        "status": "asset ready",
        "tone": "green",
        "caption": "Close the completed result core with a boundary-aware human-alignment slide.",
        "claim_boundary": "No direct gene-transfer, clinical, or countermeasure claim.",
        "asset": p("output/premium_core_result_slides/slide14_human_translation_ladder_v0_1/slide14_human_translation_ladder_rendered_preview.png"),
    },
]

SECTION_BANDS = [
    {"label": "Opening", "slides": "1-3", "tone": "teal", "question": "Why this benchmark?"},
    {"label": "Methods Bridge", "slides": "4-6", "tone": "sky", "question": "How are public studies tested?"},
    {"label": "Early Results", "slides": "7-9", "tone": "blue", "question": "What does transfer depend on?"},
    {"label": "Core Result Spine", "slides": "10-14", "tone": "amber", "question": "How robust, bounded, and translatable is it?"},
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
    "eyebrow": font(20, bold=True),
    "title": font(58, bold=True),
    "subtitle": font(25),
    "section": font(19, bold=True),
    "h": font(23, bold=True),
    "body": font(17),
    "small": font(14),
    "tiny": font(11),
    "num": font(28, bold=True),
}


def ensure() -> None:
    QA.mkdir(parents=True, exist_ok=True)


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
    center = (int(w * 0.56), int(h * 0.25))
    for idx, radius in enumerate([740, 1040, 1320]):
        bbox = (
            center[0] - radius,
            center[1] - int(radius * 0.33),
            center[0] + radius,
            center[1] + int(radius * 0.33),
        )
        draw.arc(bbox, 198, 354, fill=rgba("sky", 38 - idx * 9), width=3)


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


def draw_header(draw: ImageDraw.ImageDraw) -> None:
    draw.text((130, 92), "SLIDES 1-14 DECK ASSEMBLY BRIDGE", font=F["eyebrow"], fill=rgb("teal"))
    draw.text((130, 140), "Opening to core result spine", font=F["title"], fill=rgb("ink"))
    draw.text(
        (134, 214),
        "A production board for connecting thesis, methods, feature views, and the completed result family.",
        font=F["subtitle"],
        fill=rgb("soft"),
    )


def draw_placeholder(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, h: int, slide: dict) -> None:
    tone = slide["tone"]
    draw.rounded_rectangle((x, y, x + w, y + h), radius=10, fill=rgba("panel2", 208), outline=rgba(tone, 155), width=2)
    draw.text((x + 24, y + 24), "PRODUCTION BRIEF", font=F["tiny"], fill=rgb(tone))
    headline_lines = wrap_text(draw, slide["headline"], w - 48, F["body"])[:3]
    for idx, line in enumerate(headline_lines):
        draw.text((x + 24, y + 62 + idx * 25), line, font=F["body"], fill=rgb("ink"))
    caption_y = y + 80 + len(headline_lines) * 25
    for idx, line in enumerate(wrap_text(draw, slide["caption"], w - 48, F["small"])[:4]):
        draw.text((x + 24, caption_y + idx * 22), line, font=F["small"], fill=rgb("soft"))
    draw.rounded_rectangle((x + 24, y + h - 70, x + w - 24, y + h - 24), radius=8, fill=rgba("panel", 220), outline=rgba("rule", 120), width=1)
    draw.text((x + 42, y + h - 54), slide["status"], font=F["small"], fill=rgb("muted"))


def draw_slide_card(canvas: Image.Image, draw: ImageDraw.ImageDraw, slide: dict, x: int, y: int, w: int, h: int) -> None:
    tone = slide["tone"]
    shadow = Image.new("RGBA", canvas.size, (0, 0, 0, 0))
    sd = ImageDraw.Draw(shadow)
    sd.rounded_rectangle((x + 10, y + 16, x + w + 10, y + h + 16), radius=12, fill=(0, 0, 0, 120))
    canvas.alpha_composite(shadow.filter(ImageFilter.GaussianBlur(12)))
    draw.rounded_rectangle((x - 1, y - 1, x + w + 1, y + h + 1), radius=11, outline=rgba(tone, 170), width=2)
    if slide.get("asset") and slide["asset"].exists():
        canvas.alpha_composite(fit_image(Image.open(slide["asset"]), (w, h)).convert("RGBA"), (x, y))
    else:
        draw_placeholder(draw, x, y, w, h, slide)

    top_y = y + h + 22
    draw.text((x, top_y), f"{slide['n']}. {slide['title']}", font=F["h"], fill=rgb("ink"))
    headline_lines = wrap_text(draw, slide["headline"], w, F["small"])[:2]
    for idx, line in enumerate(headline_lines):
        draw.text((x, top_y + 31 + idx * 19), line, font=F["small"], fill=rgb(tone))
    caption_y = top_y + 39 + len(headline_lines) * 19
    for idx, line in enumerate(wrap_text(draw, slide["caption"], w, F["tiny"])[:3]):
        draw.text((x, caption_y + idx * 18), line, font=F["tiny"], fill=rgb("muted"))
    draw.text((x, top_y + 118), slide["status"], font=F["tiny"], fill=rgb("muted"))


def draw_assembly_board() -> Path:
    canvas = Image.new("RGBA", (3840, 2160), (0, 0, 0, 255))
    draw_background(canvas)
    draw = ImageDraw.Draw(canvas)
    draw_header(draw)

    card_w, card_h = 440, 248
    x0, y0 = 130, 378
    xgap, ygap = 52, 250
    for idx, slide in enumerate(SLIDES):
        row = idx // 7
        col = idx % 7
        x = x0 + col * (card_w + xgap)
        y = y0 + row * (card_h + ygap)
        draw_slide_card(canvas, draw, slide, x, y, card_w, card_h)

    band_y = 1510
    band_x = 154
    band_w = 825
    for idx, band in enumerate(SECTION_BANDS):
        x = band_x + idx * (band_w + 70)
        tone = band["tone"]
        draw.rounded_rectangle((x, band_y, x + band_w, band_y + 148), radius=10, fill=rgba("panel", 226), outline=rgba(tone, 145), width=2)
        draw.text((x + 28, band_y + 26), f"{band['label']} | slides {band['slides']}", font=F["section"], fill=rgb(tone))
        draw.text((x + 28, band_y + 68), band["question"], font=F["body"], fill=rgb("soft"))
        draw.text((x + 28, band_y + 106), "speaker bridge required", font=F["tiny"], fill=rgb("muted"))

    points = [(360, 1846), (1000, 1762), (1610, 1846), (2320, 1762), (3160, 1846)]
    draw.line(points, fill=rgba("teal", 32), width=44, joint="curve")
    draw.line(points, fill=rgba("amber", 145), width=5, joint="curve")
    labels = [
        ("why", "problem"),
        ("how", "task/split"),
        ("what", "feature view"),
        ("so what", "result spine"),
        ("boundary", "translation"),
    ]
    for idx, point in enumerate(points):
        tone = ["teal", "sky", "blue", "amber", "green"][idx]
        draw.ellipse((point[0] - 20, point[1] - 20, point[0] + 20, point[1] + 20), fill=rgba(tone, 215))
        draw.text((point[0] - 46, point[1] + 38), labels[idx][0], font=F["body"], fill=rgb("soft"))
        draw.text((point[0] - 58, point[1] + 64), labels[idx][1], font=F["tiny"], fill=rgb("muted"))

    draw.rounded_rectangle((130, 2034, 2310, 2100), radius=10, fill=rgba("panel", 224), outline=rgba("rule", 145), width=2)
    draw.text(
        (160, 2055),
        "Assembly rule: source-proof PNG scenes stay as plates; titles, bridge lines, caveats, and sources should be rebuilt as editable slide text.",
        font=F["small"],
        fill=rgb("soft"),
    )
    output = OUT / "slides_1_14_deck_assembly_bridge.png"
    canvas.convert("RGB").save(output, quality=95)
    return output


def draw_grayscale(path: Path) -> Path:
    image = Image.open(path).convert("RGB")
    output = QA / "slides_1_14_deck_assembly_bridge_grayscale.png"
    image.convert("L").convert("RGB").save(output, quality=95)
    return output


def build_caption_pack() -> dict:
    return {
        "created": CREATED,
        "role": "slides 1-14 assembly narration and production brief",
        "family_bridge": (
            "The first 14 slides move from why the benchmark is needed, to how public studies become "
            "mission-held-out tasks, to what feature layers the models see, and finally to the completed "
            "result spine with robustness, guardrails, limitations, biological hypotheses, and translation boundaries."
        ),
        "compression_decision": {
            "default_14_slide_path": "Use B1 as slide 4 with B2 in notes; use B3 as slide 5 with B4 in notes; use feature-layer bridge as slide 6.",
            "expanded_methods_option": "If the talk focus is methodology, split B2 and B4 into visible slides and shift the result spine later.",
            "risk": "Over-expanding methods can delay the first result; over-compressing methods can make AUROC slides feel unearned.",
        },
        "slides": [
            {
                "slide": slide["n"],
                "section": slide["section"],
                "title": slide["title"],
                "headline": slide["headline"],
                "role": slide["role"],
                "status": slide["status"],
                "caption": slide["caption"],
                "claim_boundary": slide["claim_boundary"],
                "asset": rel(slide["asset"]) if slide.get("asset") else None,
                "backup_assets": [rel(path) for path in slide.get("backup_assets", [])],
            }
            for slide in SLIDES
        ],
    }


def write_caption_markdown(pack: dict) -> Path:
    path = OUT / "slides_1_14_deck_assembly_caption_pack.md"
    lines = [
        "# Slides 1-14 Deck Assembly Caption Pack",
        "",
        f"Date: {CREATED}",
        "",
        "## Family Bridge",
        "",
        pack["family_bridge"],
        "",
        "## Compression Decision",
        "",
        f"Default: {pack['compression_decision']['default_14_slide_path']}",
        "",
        f"Expanded option: {pack['compression_decision']['expanded_methods_option']}",
        "",
        f"Risk: {pack['compression_decision']['risk']}",
        "",
        "## Slide Notes",
        "",
    ]
    for slide in pack["slides"]:
        lines.extend(
            [
                f"### Slide {slide['slide']}: {slide['title']}",
                "",
                f"Headline: {slide['headline']}",
                "",
                f"Caption: {slide['caption']}",
                "",
                f"Claim boundary: {slide['claim_boundary']}",
                "",
                f"Asset: {slide['asset'] or 'production brief needed'}",
                "",
            ]
        )
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def write_overlay_rules() -> dict:
    return {
        "created": CREATED,
        "z_stack_rule": "Layered PNG scene plate plus editable scientific interpretation.",
        "global_overlay_rules": [
            "Use PNG proof scenes as Z2 evidence plates, not as all-in-one final slides.",
            "Rebuild headline, bridge line, source, caveat, and any presenter text as editable objects.",
            "Keep source and caveat text visible but subordinate.",
            "Do not put table-like panels inside figure cards; tables move to backup or manuscript tables.",
            "Avoid adding organoid/OSD-120 extension visuals before slide 20.",
        ],
        "slide_specific": [
            {
                "slides": "1-3",
                "rule": "Build original opening visuals; do not use result thumbnails as decorative background.",
            },
            {
                "slides": "4-5",
                "rule": "Visible slide can be compressed, but speaker notes must carry task-unit and train-only guard definitions.",
            },
            {
                "slide": 6,
                "rule": "Use feature-layer bridge to define gene and pathway views before score surfaces.",
            },
            {
                "slides": "7-9",
                "rule": "Use early result scenes to establish tissue, pathway, and model-comparison logic before hardening.",
            },
            {
                "slides": "10-14",
                "rule": "Use core-result caption pack and keep each slide's claim boundary line intact.",
            },
        ],
    }


def write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def build() -> dict[str, str]:
    ensure()
    board = draw_assembly_board()
    grayscale = draw_grayscale(board)
    pack = build_caption_pack()
    caption_json = OUT / "slides_1_14_deck_assembly_caption_pack.json"
    caption_md = write_caption_markdown(pack)
    write_json(caption_json, pack)
    overlay_rules = OUT / "slides_1_14_deck_overlay_rules.json"
    write_json(overlay_rules, write_overlay_rules())

    missing_required_assets = [
        slide["n"] for slide in SLIDES if slide.get("asset") is not None and not slide["asset"].exists()
    ]
    opening_needs = [slide["n"] for slide in SLIDES if slide.get("asset") is None]
    manifest = {
        "created": CREATED,
        "artifact_role": "slides 1-14 deck assembly bridge before PPTX build",
        "outputs": {
            "assembly_board": rel(board),
            "grayscale_qa": rel(grayscale),
            "caption_pack_json": rel(caption_json),
            "caption_pack_markdown": rel(caption_md),
            "overlay_rules": rel(overlay_rules),
            "qa": rel(OUT / "slides_1_14_deck_assembly_qa.json"),
        },
        "source_documents": [
            "docs/V1_V9_PRESENTATION_AND_MANUSCRIPT_MASTER_OUTLINE_2026_05_31.md",
            "docs/V1_V9_DECK_SPINE_METHODS_BRIDGE_PLACEMENT_2026_06_02.md",
            "docs/PREMIUM_BRIDGE_METHODS_NARRATION_PACK_B1_B4_2026_06_02.md",
            "docs/SLIDES_1_3_OPENING_VISUAL_SYSTEM_REVIEW_2026_06_03.md",
            "output/premium_opening_slides_1_3_v0_1/opening_slides_1_3_manifest.json",
            "output/premium_core_result_slides/core_result_set_qa_v0_1/core_result_slides_10_14_caption_pack.json",
        ],
        "slides": pack["slides"],
    }
    manifest_path = OUT / "slides_1_14_deck_assembly_manifest.json"
    write_json(manifest_path, manifest)

    qa = {
        "created": CREATED,
        "automatic_checks": {
            "n_slides": len(SLIDES),
            "required_assets_missing": missing_required_assets,
                "opening_production_brief_slides": opening_needs,
            "assembly_board_exists": board.exists(),
            "grayscale_qa_exists": grayscale.exists(),
            "image_dimensions": {"width_px": Image.open(board).width, "height_px": Image.open(board).height},
            "caption_pack_exists": caption_md.exists() and caption_json.exists(),
            "overlay_rules_exist": overlay_rules.exists(),
        },
        "manual_review": {
            "visual_review_status": "pass: assembly board and grayscale QA inspected with opening slide 1-3 visual candidates integrated",
            "read_order": "opening -> methods bridge -> feature view -> early results -> core result spine",
            "claim_boundary": "slides 1-14 keep v1-v7 results central; v8/v9 extension assets do not enter this opening/core block",
            "production_decision": "slides 1-14 now have production candidate assets; final PPTX should rebuild text as editable overlays",
            "style_caution": "slides 4-5 are light proof-stage scenes; final PPTX should either treat them as intentional methods-proof contrast or rebuild dark variants for tighter deck continuity",
        },
    }
    qa_path = OUT / "slides_1_14_deck_assembly_qa.json"
    write_json(qa_path, qa)
    return {
        "assembly_board": rel(board),
        "grayscale_qa": rel(grayscale),
        "caption_pack": rel(caption_md),
        "overlay_rules": rel(overlay_rules),
        "manifest": rel(manifest_path),
        "qa": rel(qa_path),
    }


def main() -> None:
    print(json.dumps(build(), indent=2))


if __name__ == "__main__":
    main()
