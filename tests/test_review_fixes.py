import importlib
import json
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


class ReviewFixTests(unittest.TestCase):
    def read_repo_text(self, relative_path: str) -> str:
        return (REPO_ROOT / relative_path).read_text()

    def load_module(self, dotted_path: str):
        return importlib.import_module(dotted_path)

    def test_run_baselines_lr_uses_elasticnet_penalty(self):
        run_baselines = self.load_module("scripts.run_baselines")
        model = run_baselines.build_lr()
        self.assertEqual(model.named_steps["clf"].get_params()["penalty"], "elasticnet")

    def test_run_pathway_lomo_lr_uses_elasticnet_penalty(self):
        run_pathway_lomo = self.load_module("scripts.run_pathway_lomo")
        model = run_pathway_lomo.build_lr()
        self.assertEqual(model.named_steps["clf"].get_params()["penalty"], "elasticnet")

    def test_phase1_ci_threshold_is_shared(self):
        benchmark_common = self.load_module("scripts.benchmark_common")
        evaluate_submission = self.load_module("scripts.evaluate_submission")
        run_pathway_lomo = self.load_module("scripts.run_pathway_lomo")
        self.assertEqual(benchmark_common.PHASE1_CI_LOWER_THRESHOLD, 0.600)
        self.assertEqual(
            evaluate_submission.CI_LOWER_THRESHOLD,
            benchmark_common.PHASE1_CI_LOWER_THRESHOLD,
        )
        self.assertEqual(
            run_pathway_lomo.PHASE1_CI_LOWER_THRESHOLD,
            benchmark_common.PHASE1_CI_LOWER_THRESHOLD,
        )

    def test_task_variant_suffix_preserves_existing_variant_names(self):
        benchmark_common = self.load_module("scripts.benchmark_common")
        with tempfile.TemporaryDirectory() as tmpdir:
            tasks_dir = Path(tmpdir)
            canonical = tasks_dir / "A1_liver_lomo"
            combat = tasks_dir / "A1_liver_lomo_combat"
            iss_only = tasks_dir / "A1_liver_lomo_iss_only"
            canonical.mkdir()
            combat.mkdir()
            iss_only.mkdir()

            self.assertEqual(
                benchmark_common.task_variant_suffix("A1", canonical, tasks_dir),
                "",
            )
            self.assertEqual(
                benchmark_common.task_variant_suffix("A1", combat, tasks_dir),
                "_combat",
            )
            self.assertEqual(
                benchmark_common.task_variant_suffix("A1", iss_only, tasks_dir),
                "_iss_only",
            )

    def test_utils_include_skin_and_skin_mission_expansion(self):
        utils = self.load_module("scripts.utils")
        self.assertIn("skin", utils.TISSUE_MISSIONS)
        self.assertEqual(
            utils.MISSION_EXPAND["skin"]["MHU-2"],
            ["MHU-2_dorsal", "MHU-2_femoral"],
        )

    def test_generate_submission_supports_all_published_v1_a_and_b_tasks(self):
        generate_submission = self.load_module("scripts.generate_submission")
        for task_id in ["A3", "B1", "B3"]:
            self.assertIn(task_id, generate_submission.ALL_TASKS)

    def test_quality_filter_harmonizes_skin_variants(self):
        import pandas as pd
        quality_filter = self.load_module("scripts.quality_filter")

        self.assertEqual(
            quality_filter.canonical_mission_name("OSD-238", "MHU-2 (dorsal)"),
            "MHU-2",
        )
        self.assertEqual(
            quality_filter.canonical_mission_name("OSD-239", "MHU-2 (femoral)"),
            "MHU-2",
        )

        counts = pd.DataFrame(
            [[1, 2, 3]],
            index=["ENSMUSG0000001"],
            columns=[
                "Mmus_C57-6J_SKN_FLT_25days_Rep1_F1",
                "Mmus_C3H-HeJ_SKN_VIV_25days_Rep1_V2",
                "Mmus_C57-6J_SKN_VIV_25days_Rep1_V1",
            ],
        )
        labels = {
            "Mmus_C57-6J_SKN_FLT_25days_Rep1_F1": "Flight",
            "Mmus_C3H-HeJ_SKN_VIV_25days_Rep1_V2": "VC",
            "Mmus_C57-6J_SKN_VIV_25days_Rep1_V1": "VC",
        }

        filtered_counts, filtered_labels, canonical_mission, notes = quality_filter.apply_study_harmonization(
            "OSD-254",
            counts,
            labels,
        )

        self.assertEqual(canonical_mission, "RR-7")
        self.assertEqual(
            list(filtered_counts.columns),
            [
                "Mmus_C57-6J_SKN_FLT_25days_Rep1_F1",
                "Mmus_C57-6J_SKN_VIV_25days_Rep1_V1",
            ],
        )
        self.assertNotIn("Mmus_C3H-HeJ_SKN_VIV_25days_Rep1_V2", filtered_labels)
        self.assertTrue(any("C57BL/6J subset" in note for note in notes))

    def test_analysis_registries_match_shared_v1_surface(self):
        preprocess_pathways = self.load_module("scripts.preprocess_pathways")
        cross_tissue_transfer = self.load_module("scripts.cross_tissue_transfer")
        batch_correction_eval = self.load_module("scripts.batch_correction_eval")
        self.assertIn("MHU-1", preprocess_pathways.TISSUE_MISSIONS["thymus"])
        self.assertEqual(cross_tissue_transfer.MISSION_GLDS["liver"]["RR-6"], "GLDS-245")
        self.assertEqual(cross_tissue_transfer.MISSION_GLDS["liver"]["RR-9"], "GLDS-242")
        self.assertEqual(cross_tissue_transfer.MISSION_GLDS["liver"]["MHU-2"], "GLDS-617")
        self.assertEqual(cross_tissue_transfer.MISSION_GLDS["thymus"]["RR-6"], "GLDS-244")
        self.assertEqual(cross_tissue_transfer.MISSION_GLDS["thymus"]["MHU-2"], "GLDS-289")
        self.assertEqual(cross_tissue_transfer.MISSION_GLDS["thymus"]["RR-9"], "GLDS-421")
        self.assertIn("skin", batch_correction_eval.TISSUE_MISSIONS)

    def test_root_eye_registry_uses_stable_osd_identifier_with_legacy_storage_alias(self):
        utils = self.load_module("scripts.utils")
        self.assertIn("OSD-397", utils.TISSUE_MISSIONS["eye"])
        self.assertNotIn("TBD", utils.TISSUE_MISSIONS["eye"])
        self.assertEqual(utils.MISSION_FILE_ALIASES["OSD-397"], "TBD")
        self.assertEqual(utils.MISSION_ALIASES["TBD"], "OSD-397")

    def test_root_eye_scripts_use_stable_public_label(self):
        generate_tasks = self.read_repo_text("scripts/generate_tasks.py")
        fetch_osdr = self.read_repo_text("scripts/fetch_osdr.py")
        catalog = self.read_repo_text("scripts/catalog_datasets.py")
        compute_pathways = self.read_repo_text("scripts/compute_pathway_scores.R")
        run_fgsea = self.read_repo_text("scripts/run_fgsea.R")
        hpc_submit = self.read_repo_text("scripts/hpc_submit_all_tissues.sh")
        root_readme = self.read_repo_text("README.md")

        self.assertIn('MISSION_ALIASES = {"TBD": "OSD-397"}', generate_tasks)
        self.assertIn('"mission": "OSD-397"', fetch_osdr)
        self.assertIn('"mission": "OSD-397"', catalog)
        self.assertIn('list(mission = "OSD-397", dir = "TBD")', compute_pathways)
        self.assertIn('list(mission = "MHU-1", dir = "MHU-2", glds = "GLDS-289", osd = "OSD-289")', run_fgsea)
        self.assertIn('list(mission = "OSD-397", dir = "TBD", glds = "GLDS-397", osd = "OSD-397")', run_fgsea)
        self.assertIn('"A6|A6_eye_lomo|RR-1 RR-3 OSD-397|3"', hpc_submit)
        self.assertIn("| Eye | OSD-100, 194, 397 | RR-1, RR-3, OSD-397 | 37 |", root_readme)

    def test_eye_docs_use_stable_osd397_label(self):
        hf_card = self.read_repo_text("docs/hf_dataset_card.md")
        paper = self.read_repo_text("docs/V1_PAPER_CONTENT.md")
        tasks_readme = self.read_repo_text("tasks/README.md")

        self.assertIn("fold_OSD-397_test", hf_card)
        self.assertNotIn("fold_TBD_test", hf_card)
        self.assertIn("| Eye | RR-1, RR-3, OSD-397 | 3 | 37 | 2a |", paper)
        self.assertIn("| `B6` | Eye | RR-1, RR-3, OSD-397 | 6 | 0.754 [0.688, 0.838] |", tasks_readme)

    def test_scgpt_task_lists_match_v1_task_surface(self):
        finetune = self.read_repo_text("scripts/scgpt_finetune.py")
        tokenize = self.read_repo_text("scripts/scgpt_tokenize.py")

        for src in [finetune, tokenize]:
            self.assertIn('"A2": ("A2_gastrocnemius_lomo", ["RR-1", "RR-5", "RR-9"])', src)
            self.assertIn('"A3": ("A3_kidney_lomo", ["RR-1", "RR-3", "RR-7"])', src)
            self.assertIn('"A6": ("A6_eye_lomo", ["RR-1", "RR-3", "OSD-397"])', src)

    def test_scgpt_aggregator_discovers_scgpt_results_dir_from_eval_root(self):
        aggregate = self.load_module("scripts.aggregate_scgpt_results")
        aggregate_src = self.read_repo_text("scripts/aggregate_scgpt_results.py")
        with tempfile.TemporaryDirectory() as tmpdir:
            eval_root = Path(tmpdir) / "evaluation"
            scgpt_dir = eval_root / "scgpt"
            scgpt_dir.mkdir(parents=True)
            (scgpt_dir / "scgpt_whole_human_A1_RR-1_result.json").write_text("{}")

            resolved = aggregate.resolve_results_dir(eval_root, "whole_human")

            self.assertEqual(resolved, scgpt_dir)
        self.assertIn('"reference_summary": repo_relative_path(reference_summary)', aggregate_src)

    def test_scgpt_aggregator_uses_tracked_reference_summary_instead_of_hardcoded_values(self):
        aggregate = self.load_module("scripts.aggregate_scgpt_results")
        with tempfile.TemporaryDirectory() as tmpdir:
            eval_root = Path(tmpdir) / "evaluation"
            scgpt_dir = eval_root / "scgpt"
            scgpt_dir.mkdir(parents=True)

            reference_payload = {
                "tissues": {
                    "A1": {
                        "baseline_mean_auroc": 0.91234,
                        "geneformer_mean_auroc": 0.12345,
                    }
                },
                "overall": {
                    "baseline_mean": 0.81234,
                    "geneformer_mean": 0.22345,
                },
            }
            reference_path = eval_root / "geneformer_mouse_gf_all_tissues_summary.json"
            reference_path.write_text(json.dumps(reference_payload))

            resolved_reference = aggregate.resolve_reference_summary(scgpt_dir, eval_root, None)
            self.assertEqual(resolved_reference, reference_path)

            metrics = aggregate.load_reference_metrics(resolved_reference)
            summary = aggregate.summarize_task(
                "A1",
                [
                    {
                        "test_mission": "RR-1",
                        "auroc": 0.7,
                        "best_epoch": 1,
                        "n_train": 10,
                        "n_test": 5,
                    }
                ],
                "whole_human",
                metrics["tasks"],
            )

            self.assertEqual(summary["baseline_auroc"], 0.9123)
            self.assertEqual(summary["geneformer_auroc"], 0.1235)
            self.assertEqual(summary["delta_vs_baseline"], -0.2123)
            self.assertEqual(summary["delta_vs_geneformer"], 0.5765)
            self.assertEqual(metrics["overall"]["baseline"], 0.81234)
            self.assertEqual(metrics["overall"]["geneformer"], 0.22345)

    def test_scgpt_docs_and_summary_match_data_driven_paths_and_current_win_counts(self):
        aggregate_src = self.read_repo_text("scripts/aggregate_scgpt_results.py")
        results_summary = self.read_repo_text("evaluation/RESULTS_SUMMARY.md")

        self.assertIn("writes the combined summary alongside those per-fold result JSONs by default", aggregate_src)
        self.assertIn("Classical ML wins 4/6 tissues vs scGPT", results_summary)
        self.assertNotIn("Classical ML wins 5/6 tissues vs scGPT", results_summary)
        self.assertIn("*Results file: `evaluation/scgpt/scgpt_whole_human_all_tissues_summary.json`*", results_summary)

    def test_fm_docs_use_current_scgpt_and_model_count_language(self):
        readme = self.read_repo_text("README.md")
        hf_card = self.read_repo_text("docs/hf_dataset_card.md")
        paper = self.read_repo_text("docs/V1_PAPER_CONTENT.md")
        outline = self.read_repo_text("docs/PAPER_OUTLINE.md")
        results_summary = self.read_repo_text("evaluation/RESULTS_SUMMARY.md")

        self.assertIn("4 gene-expression foundation models + 3 text LLMs evaluated", readme)
        self.assertIn("| scGPT | 12L Transformer, 33M human cells | 0.666 | -0.093 |", readme)
        self.assertIn("Foundation / Language Models | 4 gene-expression FMs + 3 text LLMs", hf_card)
        self.assertIn("scGPT | 0.666 (6-tissue mean, v1)", hf_card)
        self.assertIn("Classical ML wins 4/6 tissues vs scGPT and 6/6 tissues vs Geneformer.", paper)
        self.assertIn("| scGPT vs Baseline mean delta | -0.093 | 4/6 Baseline wins |", paper)
        self.assertIn("https://github.com/jang1563/GeneLab_benchmark", paper)
        self.assertNotIn("[GitHub repository URL]", paper)
        self.assertIn("scGPT: mean 0.666 vs Classical ML 0.758 (delta=-0.093)", outline)
        self.assertIn("gene-expression FMs, text LLMs", outline)
        self.assertIn("Geneformer, scGPT, + 3 text LLMs", outline)
        self.assertIn("https://github.com/jang1563/GeneLab_benchmark", outline)
        self.assertIn("| scGPT | 0.628 | 0.685 | **0.556** | 0.782 | 0.650 | — | — | 0.666 |", results_summary)

    def test_submission_docs_match_current_v1_task_surface_and_thresholds(self):
        readme = self.read_repo_text("README.md")
        submission = self.read_repo_text("docs/submission_format.md")

        self.assertIn("| 95% CI lower | Bootstrap CI (N=2000) lower bound | > 0.600 |", readme)
        self.assertIn("| **Go/No-Go** | AUROC > 0.700 AND CI lower > 0.600 |", submission)
        self.assertIn("| `A6` | Eye | 3 | 37 | ~21k genes (log2 normalized) |", submission)
        self.assertIn("| `B6` | Eye | RR-1, RR-3, OSD-397 | 6 | `tasks/B6_eye_cross_mission/` |", submission)
        self.assertIn("Submit separate JSON files for each supported task (`A1`-`A6`, `B1`-`B6`).", submission)
        self.assertIn("https://github.com/jang1563/GeneLab_benchmark/issues", submission)
        self.assertNotIn("[TBD when public]", submission)

    def test_root_readme_repository_map_avoids_stale_hard_coded_counts(self):
        readme = self.read_repo_text("README.md")

        self.assertIn("tasks/                          <- Public task inputs (benchmark + sensitivity tasks)", readme)
        self.assertIn("v2 analysis and runner scripts", readme)
        self.assertIn("v4 result JSONs + SHAP/WGCNA outputs", readme)
        self.assertNotIn("Public task inputs (17 directories)", readme)
        self.assertNotIn("v1 pipeline scripts (35 Python/R/shell)", readme)
        self.assertNotIn("19 Python scripts", readme)
        self.assertNotIn("18 scripts (classifiers, SHAP, WGCNA, PPI)", readme)

    def test_root_readme_includes_v7_in_status_structure_and_changelog(self):
        readme = self.read_repo_text("README.md")

        self.assertIn("Version: v7.0 (2026-04-12) | Dataset freeze: 2026-03-01", readme)
        self.assertIn("Status: **v1–v7 Complete**", readme)
        self.assertIn("| **v7.0** | Unified/foundation-model benchmarking: scPRINT2, GNN/WGCNA graph baselines, cross-method synthesis, and signal hierarchy analysis | **Complete** | `v7/` |", readme)
        self.assertIn("├── v7/                             <- Unified/foundation-model benchmarking", readme)
        self.assertIn("| v7.0 | 2026-04-12 | Unified benchmark layer complete:", readme)
        self.assertNotIn("v7 Graph Neural Networks In Progress", readme)

    def test_public_release_metadata_uses_v7_consistently(self):
        readme = self.read_repo_text("README.md")
        hf_card = self.read_repo_text("docs/hf_dataset_card.md")
        citation = self.read_repo_text("CITATION.cff")

        self.assertIn("Version: v7.0 (2026-04-12) | Dataset freeze: 2026-03-01", readme)
        self.assertIn("Canonical v7.1 documentation source:", readme)
        self.assertIn("note    = {v7.0 with v7.1.1 documentation, public-card, and metadata consistency patch; data freeze 2026-03-01}", readme)
        self.assertIn("Version: v7.0 with v7.1 documentation consistency patch | Dataset freeze: 2026-03-01", hf_card)
        self.assertIn("note    = {v7.0 with v7.1.1 documentation, public-card, and metadata consistency patch; data freeze 2026-03-01}", hf_card)
        self.assertIn('version: "7.1.1"', citation)
        self.assertIn('date-released: "2026-06-05"', citation)
        self.assertIn('notes: "Manuscript in preparation; v7.1.1 documentation, public-card, and metadata consistency patch."', citation)
        self.assertIn("documentation, public-card, and metadata consistency patch", citation)
        self.assertIn('family-names: "Kim"', citation)
        self.assertIn('given-names: "Jihoon"', citation)
        self.assertIn('affiliation: "Weill Cornell Medicine"', citation)
        self.assertNotIn("Kang", citation)
        self.assertNotIn("Jaeyoung", citation)
        self.assertNotIn("JangKeun", readme + hf_card + citation)
        self.assertNotIn("blob/v3/docs/SPACEBIOBENCH", hf_card)
        self.assertNotIn('version: "5.0.0"', citation)
        self.assertNotIn("Target journal:", citation)

    def test_v2_rrrm1_wrapper_rebuilds_merged_input_and_uses_tissue_scoped_steps(self):
        wrapper = self.read_repo_text("v2/scripts/rrrm1_f2_pipeline_wrapper.sh")
        self.assertIn("rrrm1_merge_h5ad.py", wrapper)
        self.assertIn('rrrm1_initial_scanpy.py \\', wrapper)
        self.assertIn('rrrm1_broad_annotate.py \\', wrapper)
        self.assertIn('rrrm1_singlecell_hardening.py \\', wrapper)
        self.assertGreaterEqual(wrapper.count('--tissue "${TISSUE}"'), 3)

    def test_v2_rrrm1_single_tissue_entrypoints_exist(self):
        initial_scanpy = self.read_repo_text("v2/scripts/rrrm1_initial_scanpy.py")
        broad_annotate = self.read_repo_text("v2/scripts/rrrm1_broad_annotate.py")

        self.assertIn('parser.add_argument(\n        "--tissue"', initial_scanpy)
        self.assertIn("if args.tissue is None:", initial_scanpy)
        self.assertIn("save_tissue_outputs(adata, [args.tissue])", initial_scanpy)

        self.assertIn('parser.add_argument(\n        "--tissue"', broad_annotate)
        self.assertIn("tissues = [args.tissue] if args.tissue else sorted(MARKERS)", broad_annotate)

    def test_v2_f2_paths_are_repo_local_or_override_aware(self):
        f2a = self.read_repo_text("v2/scripts/rrrm1_f2a_composition.py")
        f2b = self.read_repo_text("v2/scripts/rrrm1_f2b_pseudobulk_fgsea.py")
        f2c = self.read_repo_text("v2/scripts/rrrm1_f2c_loao_classifier.py")
        f2d = self.read_repo_text("v2/scripts/rrrm1_f2d_crossspecies.py")

        for src in [f2a, f2b, f2c]:
            self.assertIn('RRRM1_BASE_DIR', src)
            self.assertIn('REPO_ROOT / "v2" / "evaluation"', src)

        self.assertIn('out_csv_dir = PROCESSED_ROOT / f"F2B_{tissue}"', f2b)
        self.assertIn('BASE_DIR_ENV = os.environ.get("RRRM1_BASE_DIR")', f2d)
        self.assertIn('base_dir / "processed" / "F2B_blood"', f2d)
        self.assertIn('REPO_ROOT / "v2" / "processed" / "F2B_blood"', f2d)

    def test_v2_f2a_requires_broad_celltype_aware_inputs(self):
        f2a = self.read_repo_text("v2/scripts/rrrm1_f2a_composition.py")
        self.assertIn('if "broad_celltype" in adata.obs.columns:', f2a)
        self.assertIn("skipping: missing broad_celltype labels", f2a)
        self.assertIn("No annotated h5ad with broad_celltype found", f2a)

    def test_v2_e2_recomputes_jaxa_on_shared_pathway_subset(self):
        e2 = self.read_repo_text("v2/scripts/mission_conservation_comparison.py")
        self.assertIn(
            "common_pathways = sorted(set(i4_nes.keys()) & set(jaxa_nes.keys()) & set(mouse_nes.keys()))",
            e2,
        )
        self.assertIn("jaxa_stats = spearman_with_ci(jaxa_vec, mouse_vec)", e2)
        self.assertIn('"JAXA_vs_mouse_full_e1"', e2)
        self.assertIn('"comparison_basis"', e2)

    def test_v2_e3_bootstrap_resamples_paired_pathway_indices(self):
        e3 = self.read_repo_text("v2/scripts/e3_cfrna_celltype_origin.py")
        self.assertIn('ct_df = ct_df[ct_df["comparison"] == "FP1_R1_vs_preflight"]', e3)
        self.assertIn("idx = rng.integers(0, n, size=n)", e3)
        self.assertIn("stats.spearmanr(x[idx], y[idx]).statistic", e3)
        self.assertNotIn("rng.choice(x", e3)
        self.assertNotIn("rng.choice(y", e3)

    def test_v2_docs_drop_stale_paths_and_missing_script_references(self):
        readme = self.read_repo_text("v2/README.md")
        summary = self.read_repo_text("v2/evaluation/V2_RESULTS_SUMMARY.md")

        self.assertNotIn("rrrm1_benchmark.py", readme)
        self.assertNotIn("v3/evaluation", readme)
        self.assertNotIn("RRRM1_SC_BENCHMARK_PLAN_V3", readme)
        self.assertIn("RRRM1_DOWNSTREAM_PLAN_2026-03-11.md", readme)

        self.assertNotIn("v2/processed/E_crossspecies/", summary)
        self.assertIn("v2/evaluation/E1_crossspecies_nes.json", summary)
        self.assertIn("v2/evaluation/E2_mission_conservation.json", summary)
        self.assertIn("v2/evaluation/E3_cfrna_origin.json", summary)

    def test_v2_eye_llm_and_radiation_surfaces_use_osd397(self):
        radiation = self.read_repo_text("v2/scripts/annotate_radiation.py")
        task_meta = self.read_repo_text("v2/processed/llm_prompts/A6/task_meta.json")
        prompts = self.read_repo_text("v2/processed/llm_prompts/A6/fold_OSD-397_test_prompts.json")
        submission = self.read_repo_text("v2/evaluation/submission_gemini-2.5-flash_zeroshot_A6.json")

        self.assertIn('"OSD-397": 35', radiation)
        self.assertNotIn('"TBD": 35', radiation)
        self.assertIn("fold_OSD-397_test", task_meta)
        self.assertIn("fold_OSD-397_test", prompts)
        self.assertIn("fold_OSD-397_test", submission)

    def test_v4_eye_registry_uses_stable_osd_identifier_with_legacy_storage_alias(self):
        v4_utils = self.read_repo_text("v4/scripts/v4_utils.py")
        multi_method = self.read_repo_text("v4/scripts/multi_method_eval.py")
        cond_v4 = self.read_repo_text("v4/scripts/condition_prediction_v4.py")
        pathway_gen = self.read_repo_text("v4/scripts/generate_pathway_scores.py")

        self.assertIn('"eye":            ["RR-1", "RR-3", "OSD-397"]', v4_utils)
        self.assertIn('MISSION_FILE_ALIASES = {"OSD-397": "TBD"}', v4_utils)
        self.assertIn('storage_mission = MISSION_FILE_ALIASES.get(mission, mission)', multi_method)
        self.assertIn('storage_mission = MISSION_FILE_ALIASES.get(mission, mission)', cond_v4)
        self.assertIn('storage_mission = MISSION_FILE_ALIASES.get(mission, mission)', pathway_gen)

    def test_v4_and_v7_generated_figures_are_refreshed(self):
        fig1 = self.read_repo_text("v4/figures/html/Fig1_benchmark.html")
        fig7 = self.read_repo_text("v7/figures/html/v7_methods_comparison.html")

        self.assertIn('"fold": "OSD-397"', fig1)
        self.assertNotIn('"fold": "TBD"', fig1)
        self.assertIn('"test_mission": "OSD-397"', fig7)
        self.assertNotIn('"test_mission": "TBD"', fig7)

    def test_gnn_existing_result_must_match_invocation(self):
        gnn_wgcna = self.load_module("v7.unified.gnn_wgcna")
        existing = {
            "tissue": "liver",
            "graph_type": "wgcna",
            "topology_scope": "train_fold",
            "n_edges_per_gene": 10,
            "n_bootstrap": 1000,
            "n_perm": 100,
            "seed": 42,
        }
        self.assertTrue(
            gnn_wgcna.result_matches_request(
                existing,
                tissue="liver",
                graph_type="wgcna",
                topology_scope="train_fold",
                n_edges_per_gene=10,
                n_bootstrap=1000,
                n_perm=100,
                seed=42,
            )
        )
        self.assertFalse(
            gnn_wgcna.result_matches_request(
                existing,
                tissue="liver",
                graph_type="wgcna",
                topology_scope="full_dataset",
                n_edges_per_gene=10,
                n_bootstrap=1000,
                n_perm=100,
                seed=42,
            )
        )

    def test_scprint_cache_and_results_are_checkpoint_scoped(self):
        scprint2_benchmark = self.load_module("v7.unified.scprint2_benchmark")
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = Path(tmpdir)
            ckpt_a = Path("/tmp/medium-v1.5_fixed.ckpt")
            ckpt_b = Path("/tmp/alternate.ckpt")
            cache_a = scprint2_benchmark.embedding_cache_file(cache_dir, "liver", ckpt_a)
            cache_b = scprint2_benchmark.embedding_cache_file(cache_dir, "liver", ckpt_b)

            self.assertNotEqual(cache_a, cache_b)

            existing = {
                "method": "scPRINT-2",
                "tissue": "liver",
                "ckpt_path": str(ckpt_a),
                "n_bootstrap": 1000,
                "n_perm": 100,
                "seed": 42,
            }
            self.assertTrue(
                scprint2_benchmark.result_matches_request(
                    existing,
                    tissue="liver",
                    ckpt_path=ckpt_a,
                    n_bootstrap=1000,
                    n_perm=100,
                    seed=42,
                )
            )
            self.assertFalse(
                scprint2_benchmark.result_matches_request(
                    existing,
                    tissue="liver",
                    ckpt_path=ckpt_b,
                    n_bootstrap=1000,
                    n_perm=100,
                    seed=42,
                )
            )

    def test_v3_eye_registry_uses_stable_osd_identifier(self):
        scfoundation = self.read_repo_text("v3/scripts/scfoundation_eval.py")
        uce = self.read_repo_text("v3/scripts/uce_eval.py")
        b_ext = self.read_repo_text("v3/scripts/b_ext_transfer_matrix.py")

        for src in [scfoundation, uce, b_ext]:
            self.assertIn("OSD-397", src)

        self.assertIn('MISSION_FILE_ALIASES = {"OSD-397": "TBD"}', b_ext)

    def test_v3_readme_drops_missing_doc_references(self):
        readme = self.read_repo_text("v3/README.md")
        self.assertNotIn("docs/PLAN_V3.md", readme)
        self.assertNotIn("docs/HPC_FM_RUNBOOK.md", readme)
        self.assertIn("docs/DATA_CATALOG_V3.md", readme)
        self.assertIn("Tracked example SLURM logs", readme)

    def test_v6_utils_requires_explicit_spaceomics_root(self):
        v6_utils = self.read_repo_text("v6/scripts/v6_utils.py")
        self.assertIn('SPACEOMICS_ROOT = os.environ.get("SPACEOMICS_ROOT")', v6_utils)
        self.assertNotIn('Path(os.environ.get("SPACEOMICS_ROOT", "."))', v6_utils)
        self.assertIn("Set SPACEOMICS_ROOT to a local SpaceOmicsBench checkout", v6_utils)

    def test_v7_scprint_hpc_and_cli_share_checkpoint_contract(self):
        hpc_scprint = self.read_repo_text("v7/scripts/hpc_scprint.sh")
        hpc_setup = self.read_repo_text("v7/scripts/hpc_setup_v7.sh")
        scprint = self.read_repo_text("v7/unified/scprint2_benchmark.py")

        self.assertIn('CKPT_PATH="${SCPRINT_CKPT_PATH:-$BASE_DIR/v7/models/scprint2/medium-v1.5.ckpt}"', hpc_scprint)
        self.assertIn("medium-v1.5.ckpt", hpc_setup)
        self.assertIn('parser.add_argument("--ckpt-path", default="medium-v1.5.ckpt"', scprint)
        self.assertIn('add_candidate(model_dir / "medium-v1.5_fixed.ckpt")', scprint)
        self.assertNotIn("/home/fs01/jak4013/Dropbox", hpc_scprint)

    def test_v7_scprint_artifacts_store_portable_checkpoint_paths(self):
        scprint = self.read_repo_text("v7/unified/scprint2_benchmark.py")
        eval_json = self.read_repo_text("v7/evaluation/SCPRINT2_eye.json")
        fig_html = self.read_repo_text("v7/figures/html/v7_methods_comparison.html")

        self.assertIn("if not ckpt_path:", scprint)
        self.assertIn('return str(Path(*parts[v7_idx:]).with_name(canonical_name))', scprint)
        self.assertIn('"ckpt_path": "v7/models/scprint2/medium-v1.5.ckpt"', eval_json)
        self.assertIn('"ckpt_path": "v7/models/scprint2/medium-v1.5.ckpt"', fig_html)
        self.assertNotIn("/home/fs01/jak4013/Dropbox", eval_json)
        self.assertNotIn("/home/fs01/jak4013/Dropbox", fig_html)

    def test_v7_signal_hierarchy_uses_real_fm_outputs_and_explicit_spaceomics_path(self):
        signal_hierarchy = self.read_repo_text("v7/unified/signal_hierarchy.py")
        self.assertIn("SPACEOMICS_RESULTS_JSON", signal_hierarchy)
        self.assertIn('V3_EVAL_DIR / "FM_scfoundation.json"', signal_hierarchy)
        self.assertIn('V3_EVAL_DIR / "FM_uce.json"', signal_hierarchy)
        self.assertIn('ROOT_EVAL_DIR / "geneformer_mouse_gf_all_tissues_summary.json"', signal_hierarchy)
        self.assertIn('ROOT_EVAL_DIR / "scgpt" / "scgpt_whole_human_all_tissues_summary.json"', signal_hierarchy)
        self.assertNotIn("MEMORY.md", signal_hierarchy)

    def test_v7_methods_figure_uses_geneformer_summary_instead_of_memory_constant(self):
        fig = self.read_repo_text("v7/unified/fig_v7_methods.py")
        self.assertIn('geneformer_mouse_gf_all_tissues_summary.json', fig)
        self.assertIn("def load_geneformer_summary()", fig)
        self.assertNotIn("MEMORY.md", fig)
        self.assertNotIn("gf_mean = 0.476", fig)

    def test_release_hygiene_ignores_local_and_large_v8_artifacts(self):
        gitignore = self.read_repo_text(".gitignore")
        self.assertIn(".claude/", gitignore)
        self.assertIn("v8/bridge/geo_cache/", gitignore)
        self.assertIn("v8/**/__pycache__/", gitignore)

    def test_release_hygiene_no_personal_paths_in_public_sources(self):
        checked_files = [
            "v8/README.md",
            "v8/bridge/tissue_nes_bridge.py",
            "v8/bridge/supervised_conservation.py",
            "v8/bridge/link_spaceomicsbench.py",
            "v8/bridge/leakage_audit.py",
            "v8/multiomics/propagation.py",
            "v8/RESULTS_SUMMARY.py",
        ]
        forbidden = ["/Users/jak4013", "~/.claude"]
        offenders = [
            rel
            for rel in checked_files
            if any(token in self.read_repo_text(rel) for token in forbidden)
        ]
        self.assertEqual(offenders, [])

    def test_v8_bridge_leakage_audit_is_part_of_hpc_bridge_gate(self):
        bridge_sh = self.read_repo_text("scripts/hpc_v8_bridge.sh")
        audit = self.read_repo_text("v8/bridge/leakage_audit.py")
        readme = self.read_repo_text("v8/bridge/evaluation/README.md")

        self.assertIn("v8/bridge/leakage_audit.py", bridge_sh)
        self.assertIn("label_excluded_from_features", audit)
        self.assertIn("StratifiedKFold", audit)
        self.assertIn("bridge_leakage_audit.json", readme)

    def test_v8_intervene_api_snapshot_records_payload_hashes_without_requery(self):
        snapshot = self.read_repo_text("v8/intervene/api_snapshot_manifest.py")
        readme = self.read_repo_text("v8/intervene/evaluation/README.md")

        self.assertIn("not_recalled_by_this_script", snapshot)
        self.assertIn("payload_sha256", snapshot)
        self.assertIn("parsed_output_sha256", snapshot)
        self.assertIn("api_snapshot_manifest.json", readme)

    def test_v8_decompose_raw_cache_audit_records_full_rerun_readiness(self):
        audit = self.read_repo_text("v8/decompose/raw_cache_audit.py")
        factorial = self.read_repo_text("v8/decompose/factorial_analog.py")
        readme = self.read_repo_text("v8/decompose/evaluation/README.md")
        hpc = self.read_repo_text("scripts/hpc_v8_decompose.sh")

        self.assertIn("counts_candidates", factorial)
        self.assertIn("DECOMPOSE factorial failed for", factorial)
        self.assertIn("full_rerun_ready", audit)
        self.assertIn("missing_files", audit)
        self.assertIn("raw_cache_audit.json", readme)
        self.assertIn("v8/decompose/raw_cache_audit.py", hpc)

    def test_v8_summary_uses_decompose_variance_not_sig_count_fraction(self):
        summary = self.read_repo_text("v8/RESULTS_SUMMARY.py")

        self.assertIn("variance_attribution_top200", summary)
        self.assertNotIn("int_count / max(total_sig, 1)", summary)

    def test_v8_decompose_mars_saturation_is_bounded_sensitivity_not_fit(self):
        saturation = self.read_repo_text("v8/decompose/mars_saturation_sensitivity.py")
        hpc = self.read_repo_text("scripts/hpc_v8_decompose.sh")
        readme = self.read_repo_text("v8/decompose/evaluation/README.md")

        self.assertIn("cap5", saturation)
        self.assertIn("sqrt_after5", saturation)
        self.assertIn("log_after5", saturation)
        self.assertIn("not fitted nonlinear dose-response models", saturation)
        self.assertIn("v8/decompose/mars_saturation_sensitivity.py", hpc)
        self.assertIn("mars_saturation_summary.json", readme)

    def test_v8_intervene_safety_triage_keeps_candidates_hypothesis_only(self):
        triage = self.read_repo_text("v8/intervene/safety_triage.py")
        readme = self.read_repo_text("v8/intervene/evaluation/README.md")
        hpc = self.read_repo_text("scripts/hpc_v8_intervene.sh")

        self.assertIn("known_toxicity_class", triage)
        self.assertIn("hypothesis-generating target/pathway triage only", triage)
        self.assertIn("safety_triage.csv", readme)
        self.assertIn("v8/intervene/safety_triage.py", hpc)

    def test_v8_beta_gate_validates_provenance_and_freeze_metadata(self):
        validator = self.read_repo_text("scripts/validate_v8_provenance.py")
        release_gate = self.read_repo_text("scripts/hpc_release_validate.sh")
        rebuild = self.read_repo_text("scripts/hpc_v8_beta_rebuild.sh")
        input_freeze = self.read_repo_text("v8/provenance/input_freeze.json")
        artifact_manifest = self.read_repo_text("v8/release/v8_beta_artifact_manifest.json")
        beta_plan = self.read_repo_text("docs/V8_BETA_RELEASE_PLAN_2026_05_10.md")

        self.assertIn("validate_run_manifest", validator)
        self.assertIn("validate_v8_provenance.py", release_gate)
        self.assertIn("hpc_v8_bridge.sh", rebuild)
        self.assertIn("hpc_v8_decompose.sh", rebuild)
        self.assertIn("hpc_v8_intervene.sh", rebuild)
        self.assertIn("v8/figures/generate_main_figures.py", rebuild)
        self.assertIn('"status": "release_candidate"', input_freeze)
        self.assertIn("spaceomicsbench.v2_public", input_freeze)
        self.assertIn("l1000cds2.api_snapshot", input_freeze)
        self.assertIn("hugging_face_dataset", artifact_manifest)
        self.assertIn("zenodo", artifact_manifest)
        self.assertIn("Still open before declaring v8.0-beta frozen", beta_plan)


if __name__ == "__main__":
    unittest.main()
