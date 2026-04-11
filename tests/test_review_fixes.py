import importlib
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

    def test_v7_signal_hierarchy_uses_real_fm_outputs_and_explicit_spaceomics_path(self):
        signal_hierarchy = self.read_repo_text("v7/unified/signal_hierarchy.py")
        self.assertIn("SPACEOMICS_RESULTS_JSON", signal_hierarchy)
        self.assertIn('V3_EVAL_DIR / "FM_scfoundation.json"', signal_hierarchy)
        self.assertIn('V3_EVAL_DIR / "FM_uce.json"', signal_hierarchy)
        self.assertIn('ROOT_EVAL_DIR / "geneformer_mouse_gf_all_tissues_summary.json"', signal_hierarchy)
        self.assertIn('ROOT_EVAL_DIR / "scgpt" / "scgpt_whole_human_all_tissues_summary.json"', signal_hierarchy)
        self.assertNotIn("MEMORY.md", signal_hierarchy)


if __name__ == "__main__":
    unittest.main()
