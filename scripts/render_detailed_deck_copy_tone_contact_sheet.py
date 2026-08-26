#!/usr/bin/env python3
"""Render a contact sheet for detailed-deck assets after copy-tone edits."""

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
)
OUT_DIR = ASSET_DIR / "copy_tone_qa"
OUT_DIR.mkdir(parents=True, exist_ok=True)

ASSETS = [
    ("SpaceBio-Bench opening title", ASSET_DIR / "opening_title" / "spacebiobench_opening_title_premium.png"),
    ("Why this needs a benchmark", ASSET_DIR / "benchmark_need" / "why_this_needs_a_benchmark_premium.png"),
    ("Evidence layers organize the story", ASSET_DIR / "evidence_layers" / "evidence_layers_set_claim_strength_premium.png"),
    ("What counts as a task", ASSET_DIR / "task_record" / "what_counts_as_a_task_premium.png"),
    ("From source record to matrix", ASSET_DIR / "source_to_matrix" / "from_source_record_to_matrix_premium.png"),
    ("What the model actually sees", ASSET_DIR / "model_input_surface" / "what_the_model_actually_sees_premium.png"),
    ("Mission-held-out protocol", ASSET_DIR / "mission_holdout" / "mission_heldout_protocol_premium.png"),
    ("Train-only leakage guard", ASSET_DIR / "train_only_guard" / "train_only_feature_guard_premium.png"),
    ("Metric primer", ASSET_DIR / "metric_primer" / "metric_primer_auroc_uncertainty_premium.png"),
    ("Model-family bridge", ASSET_DIR / "model_family_bridge" / "model_family_bridge_premium.png"),
    ("Model-family primer", ASSET_DIR / "model_family" / "model_family_input_surface_premium.png"),
    ("Evidence scope ladder", ASSET_DIR / "evidence_scope_ladder" / "evidence_scope_ladder_premium.png"),
    ("Worked tissue score", ASSET_DIR / "worked_tissue_score" / "worked_tissue_score_thymus_transfer_premium.png"),
    ("Tissue transfer hierarchy", ASSET_DIR / "tissue_hierarchy" / "tissue_transfer_hierarchy_premium.png"),
    ("Liver heterogeneity explainer", ASSET_DIR / "liver_heterogeneity" / "liver_mission_heterogeneity_premium.png"),
    ("Transfer matrix behind ranking", ASSET_DIR / "transfer_matrix" / "transfer_matrix_behind_ranking_premium.png"),
    ("Pathway nuisance-signal check", ASSET_DIR / "pathway_nuisance" / "pathway_features_reduce_selected_nuisance_premium.png"),
    ("NES conservation predicts transfer", ASSET_DIR / "nes_conservation" / "nes_conservation_predicts_transfer_premium.png"),
    ("Screen transfer feasibility", ASSET_DIR / "transfer_feasibility" / "screen_transfer_feasibility_before_training_premium.png"),
    ("Held-out validation", ASSET_DIR / "heldout_validation" / "heldout_missions_confirm_signal_premium.png"),
    ("Negative controls", ASSET_DIR / "negative_controls" / "negative_controls_anchor_readout_premium.png"),
    ("Strong baselines", ASSET_DIR / "strong_baselines" / "strong_baselines_make_model_claims_meaningful_premium.png"),
    ("Classical ML result surface", ASSET_DIR / "classical_result_surface" / "classical_ml_result_surface_premium.png"),
    ("Foundation model primer", ASSET_DIR / "foundation_model_adds" / "what_a_foundation_model_adds_premium.png"),
    ("Text LLM prompt diagnostics", ASSET_DIR / "text_llm_prompt_diagnostics" / "text_llm_prompt_diagnostics_premium.png"),
    ("Model tier result", ASSET_DIR / "model_tier_result" / "model_tier_result_task_fit_beats_scale_premium.png"),
    ("Bulk adaptation surface", ASSET_DIR / "bulk_surface_hard" / "bulk_rnaseq_hard_surface_for_cell_fms_premium.png"),
    ("Method hardening grid", ASSET_DIR / "method_hardening" / "method_hardening_preserves_main_readout_premium.png"),
    ("New model ideas", ASSET_DIR / "new_model_ideas" / "newer_model_ideas_preserve_the_lesson_premium.png"),
    ("DGE robustness", ASSET_DIR / "dge_robustness" / "dge_pipeline_robustness_premium.png"),
    ("External biology validation", ASSET_DIR / "external_validation" / "external_biology_validation_premium.png"),
    ("Evidence stack makes scores readable", ASSET_DIR / "evidence_stack" / "evidence_stack_turns_scores_into_claim_status_premium.png"),
    ("Core benchmark takeaway", ASSET_DIR / "core_takeaway" / "core_benchmark_takeaway_premium.png"),
    ("Temporal labels need context", ASSET_DIR / "temporal_guardrails" / "temporal_labels_need_guardrails_premium.png"),
    ("Single-cell pilots provide context", ASSET_DIR / "single_cell_pilots" / "single_cell_pilots_provide_context_premium.png"),
    ("Spatial weak-signal cases sharpen the readout", ASSET_DIR / "spatial_scope" / "spatial_weak_signal_cases_define_scope_premium.png"),
    ("Systems biology adds interpretation", ASSET_DIR / "systems_biology" / "systems_biology_adds_interpretation_premium.png"),
    ("Immune and TF activity prioritize follow-up", ASSET_DIR / "immune_tf_activity" / "immune_tf_activity_prioritize_followup_premium.png"),
    ("Target and biomarker layers are triage", ASSET_DIR / "target_biomarker" / "target_biomarker_layers_are_triage_premium.png"),
    ("Pathways carry mouse-to-human transfer", ASSET_DIR / "mouse_human_transfer" / "pathways_carry_mouse_to_human_transfer_premium.png"),
    ("Translation details shape the readout", ASSET_DIR / "translation_details" / "translation_details_define_claim_scope_premium.png"),
    ("Prediction becomes hypothesis through layered readouts", ASSET_DIR / "prediction_hypothesis" / "prediction_becomes_hypothesis_through_evidence_gates_premium.png"),
    ("v8 starts after the core readout", ASSET_DIR / "v8_transition" / "v8_starts_after_the_evidence_gates_premium.png"),
    ("v8 has four incubator pillars", ASSET_DIR / "v8_pillars" / "v8_has_four_incubator_pillars_premium.png"),
    ("Mouse pathways improve human prediction", ASSET_DIR / "bridge_human_prediction" / "mouse_pathways_improve_human_prediction_premium.png"),
    ("Stressor interactions dominate combined effects", ASSET_DIR / "decompose_interactions" / "stressor_interactions_dominate_combined_effects_premium.png"),
    ("Perturbation hits prioritize follow-up axes", ASSET_DIR / "intervene_prioritization" / "perturbation_hits_prioritize_followup_axes_premium.png"),
    ("Causal maps connect stressor axes to tissue readouts", ASSET_DIR / "causal_map" / "causal_maps_connect_stressor_axes_to_tissue_readouts_premium.png"),
    ("v8 turns benchmark evidence into follow-up workstreams", ASSET_DIR / "v8_summary" / "v8_turns_benchmark_evidence_into_followup_workstreams_premium.png"),
    ("Platform turn packages results into a reproducible benchmark", ASSET_DIR / "platform_turn" / "platform_turn_packages_results_into_reproducible_benchmark_premium.png"),
    ("Manifest evaluator and run record make scores rerunnable", ASSET_DIR / "platform_object_chain" / "manifest_evaluator_run_record_make_scores_rerunnable_premium.png"),
    ("Task registry and metric profiles keep runs comparable", ASSET_DIR / "task_registry_metric_profiles" / "task_registry_and_metric_profiles_keep_runs_comparable_premium.png"),
    ("Public bulk metadata alpha makes the bulk core inspectable", ASSET_DIR / "public_bulk_metadata_alpha" / "public_bulk_metadata_alpha_makes_the_bulk_core_inspectable_premium.png"),
    ("File readiness ladder turns metadata into a release package", ASSET_DIR / "file_readiness_ladder" / "file_readiness_ladder_turns_metadata_into_release_package_premium.png"),
    ("Data package makes the metadata alpha portable", ASSET_DIR / "data_package_catalog" / "data_package_makes_the_metadata_alpha_portable_premium.png"),
    ("Organoid extension is a biology check", ASSET_DIR / "organoid_biology_check" / "organoid_extension_is_a_biology_check_premium.png"),
    ("OSD-120 is a same-study interaction check", ASSET_DIR / "osd120_interaction_check" / "osd120_is_a_same_study_interaction_check_premium.png"),
    ("Single-cell scoring needs a fixed input set", ASSET_DIR / "single_cell_fixed_inputs" / "single_cell_scoring_needs_fixed_inputs_premium.png"),
    ("Roadmap to release", ASSET_DIR / "release_roadmap" / "roadmap_to_release_premium.png"),
    ("What this deck establishes", ASSET_DIR / "final_establishes" / "what_this_deck_establishes_premium.png"),
]


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
F_SMALL = font(18)


def fit(image: Image.Image, size: tuple[int, int]) -> Image.Image:
    image = image.convert("RGB").copy()
    resample = getattr(getattr(Image, "Resampling", Image), "LANCZOS")
    image.thumbnail(size, resample)
    canvas = Image.new("RGB", size, "#0C111A")
    canvas.paste(image, ((size[0] - image.width) // 2, (size[1] - image.height) // 2))
    return canvas


def main() -> None:
    thumb_w, thumb_h = 720, 405
    pad_x, pad_y = 70, 86
    cols = 2
    rows = (len(ASSETS) + cols - 1) // cols
    width = pad_x * 2 + cols * thumb_w + (cols - 1) * 60
    height = 150 + rows * (thumb_h + 76) + pad_y
    canvas = Image.new("RGB", (width, height), "#0C111A")
    draw = ImageDraw.Draw(canvas)
    draw.text((pad_x, 42), "Detailed deck copy-tone QA", font=F_TITLE, fill="#F3F7FC")
    draw.text((pad_x, 96), "Visible copy pass: scope language without defensive framing.", font=F_SMALL, fill="#A8B4C4")

    entries = []
    for i, (label, path) in enumerate(ASSETS):
        row, col = divmod(i, cols)
        x = pad_x + col * (thumb_w + 60)
        y = 150 + row * (thumb_h + 76)
        image = fit(Image.open(path), (thumb_w, thumb_h))
        draw.rounded_rectangle((x - 4, y - 4, x + thumb_w + 4, y + thumb_h + 4), radius=18, outline="#2A394D", width=2)
        canvas.paste(image, (x, y))
        draw.text((x, y + thumb_h + 18), label, font=F_LABEL, fill="#F3F7FC")
        entries.append({"label": label, "path": str(path.relative_to(ROOT))})

    out = OUT_DIR / "detailed_deck_copy_tone_contact_sheet.png"
    canvas.save(out, quality=95)
    manifest = OUT_DIR / "detailed_deck_copy_tone_contact_sheet_manifest.json"
    manifest.write_text(json.dumps({"contact_sheet": str(out.relative_to(ROOT)), "assets": entries}, indent=2) + "\n")
    print(json.dumps({"contact_sheet": str(out.relative_to(ROOT)), "manifest": str(manifest.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()
