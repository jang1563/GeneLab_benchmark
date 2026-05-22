import csv
import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


class V9SpaceBioBenchTests(unittest.TestCase):
    def test_mission_discrimination_scores_perfect_cluster_layout(self):
        from spacebio_bench.metrics import mission_discrimination_score

        result = mission_discrimination_score(
            embeddings=[
                [0.0, 0.0],
                [0.0, 0.2],
                [10.0, 0.0],
                [10.0, 0.2],
                [0.0, 10.0],
                [0.2, 10.0],
            ],
            mission_labels=["RR-1", "RR-1", "RR-3", "RR-3", "RR-9", "RR-9"],
        )

        self.assertEqual(result["metric_id"], "mission_discrimination")
        self.assertEqual(result["n_scored"], 6)
        self.assertEqual(result["n_skipped"], 0)
        self.assertAlmostEqual(result["score"], 1.0)

    def test_mission_discrimination_complete_tie_scores_random_expectation(self):
        from spacebio_bench.metrics import mission_discrimination_score

        result = mission_discrimination_score(
            embeddings=[
                [1.0, 1.0],
                [1.0, 1.0],
                [1.0, 1.0],
                [1.0, 1.0],
            ],
            mission_labels=["RR-1", "RR-1", "RR-3", "RR-3"],
        )

        self.assertAlmostEqual(result["score"], 0.5)

    def test_mission_discrimination_skips_singleton_own_centroids(self):
        from spacebio_bench.metrics import mission_discrimination_score

        result = mission_discrimination_score(
            embeddings=[
                [0.0, 0.0],
                [0.0, 0.2],
                [10.0, 0.0],
            ],
            mission_labels=["RR-1", "RR-1", "RR-3"],
        )

        self.assertEqual(result["n_scored"], 2)
        self.assertEqual(result["n_skipped"], 1)

    def test_task_manifest_validator_accepts_minimal_v9_contract(self):
        from spacebio_bench import validate_task_manifest

        manifest = {
            "schema_version": "0.1.0",
            "task_id": "v9_demo_bulk_lomo",
            "task_family": "bulk_lomo",
            "title": "Demo bulk mission-held-out task",
            "source_records": [
                {
                    "source_id": "osdr.demo",
                    "url_or_accession": "OSD-000",
                    "access_status": "public",
                    "checksum_status": "pending",
                }
            ],
            "split": {
                "name": "leave_one_mission_out",
                "unit": "mission",
                "strategy": "train on all missions except held-out mission",
            },
            "metrics": [
                {
                    "metric_id": "mission_discrimination",
                    "profile": "genelab_minimal",
                    "interpretation": "Ranks whether samples are nearest to their own mission centroid.",
                }
            ],
            "output": {
                "prediction_format": "csv",
                "primary_artifacts": ["predictions.csv", "metrics.json"],
            },
        }

        self.assertTrue(validate_task_manifest(manifest))

    def test_metric_profiles_define_v9_core_profiles(self):
        from spacebio_bench import METRIC_PROFILES, get_metric_profile

        self.assertIn("genelab_minimal", METRIC_PROFILES)
        self.assertIn("genelab_sc", METRIC_PROFILES)
        self.assertIn(
            "mission_discrimination",
            get_metric_profile("genelab_minimal")["metrics"],
        )
        self.assertIn("de_overlap_at_n", get_metric_profile("genelab_sc")["metrics"])
        self.assertIn(
            "signature_rank_correlation",
            get_metric_profile("genelab_organoid_pilot")["metrics"],
        )

    def test_legacy_task_info_converter_exports_bulk_lomo_manifest(self):
        from spacebio_bench.tasks import legacy_task_info_to_manifest, load_legacy_task_info

        task_dir = REPO_ROOT / "tasks" / "A2_gastrocnemius_lomo"
        manifest = legacy_task_info_to_manifest(
            load_legacy_task_info(task_dir / "task_info.json"),
            task_dir=task_dir,
        )

        self.assertEqual(manifest["task_id"], "A2_gastrocnemius_bulk_lomo")
        self.assertEqual(manifest["task_family"], "bulk_lomo")
        self.assertEqual(manifest["organism"], "Mus musculus")
        self.assertEqual(manifest["assay_modality"], "bulk_rna_seq")
        self.assertEqual(manifest["split"]["missions"], ["RR-1", "RR-5", "RR-9"])
        self.assertEqual(
            [source["source_id"] for source in manifest["source_records"]],
            ["OSD-101", "OSD-326", "OSD-401"],
        )
        self.assertEqual(manifest["source_records"][0]["feature_namespace"], "mouse_gene")

    def test_legacy_variant_converter_filters_sources_to_split_missions(self):
        from spacebio_bench.tasks import legacy_task_info_to_manifest, load_legacy_task_info

        task_dir = REPO_ROOT / "tasks" / "A1_liver_lomo_iss_only"
        manifest = legacy_task_info_to_manifest(
            load_legacy_task_info(task_dir / "task_info.json"),
            task_dir=task_dir,
        )

        self.assertEqual(manifest["task_id"], "A1_liver_bulk_lomo_iss_only")
        self.assertEqual(manifest["variant"], "iss_only")
        self.assertEqual(manifest["legacy_task_id"], "A1_iss_only")
        self.assertNotIn(
            "OSD-686",
            [source["source_id"] for source in manifest["source_records"]],
        )

    def test_generated_v9_task_manifests_validate(self):
        from spacebio_bench import load_task_manifest

        manifest_dir = REPO_ROOT / "v9" / "task_manifests"
        manifests = sorted(manifest_dir.glob("*.json"))
        self.assertEqual(len(manifests), 8)
        for manifest_path in manifests:
            manifest = load_task_manifest(manifest_path)
            self.assertEqual(manifest["task_family"], "bulk_lomo")
            self.assertTrue(manifest["source_records"])

    def test_task_registry_loads_generated_manifests(self):
        from spacebio_bench import TaskRegistry

        registry = TaskRegistry.from_dir(REPO_ROOT / "v9" / "task_manifests")

        self.assertEqual(len(registry), 8)
        self.assertEqual(registry.get("A2_gastrocnemius_bulk_lomo")["tissue"], "gastrocnemius")
        self.assertIn("A1_liver_bulk_lomo_iss_only", registry.task_ids())

    def test_task_index_files_cover_generated_manifests(self):
        index_path = REPO_ROOT / "v9" / "task_manifest_index.csv"
        with index_path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))

        self.assertEqual(len(rows), 8)
        iss_only = next(row for row in rows if row["task_id"] == "A1_liver_bulk_lomo_iss_only")
        self.assertEqual(iss_only["n_sources"], "5")
        self.assertNotIn("OSD-686", iss_only["source_ids"])

    def test_task_index_writer_round_trips_registry_to_temp_outputs(self):
        from spacebio_bench import TaskRegistry, write_task_index

        registry = TaskRegistry.from_dir(REPO_ROOT / "v9" / "task_manifests")
        with tempfile.TemporaryDirectory() as tmpdir:
            csv_path, json_path = write_task_index(
                registry,
                csv_path=Path(tmpdir) / "index.csv",
                json_path=Path(tmpdir) / "index.json",
            )

            self.assertTrue(csv_path.exists())
            self.assertTrue(json_path.exists())

    def write_submission(self, path: Path, rows: list[dict[str, object]]) -> None:
        fieldnames = list(rows[0])
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

    def synthetic_submission_rows(self) -> list[dict[str, object]]:
        return [
            {
                "task_id": "A2_gastrocnemius_bulk_lomo",
                "sample_id": "rr1_f1",
                "mission": "RR-1",
                "true_label": "Flight",
                "predicted_label": "Flight",
                "flight_probability": 0.91,
                "embedding_0": 0.0,
                "embedding_1": 0.0,
            },
            {
                "task_id": "A2_gastrocnemius_bulk_lomo",
                "sample_id": "rr1_g1",
                "mission": "RR-1",
                "true_label": "Ground",
                "predicted_label": "Ground",
                "flight_probability": 0.08,
                "embedding_0": 0.0,
                "embedding_1": 0.2,
            },
            {
                "task_id": "A2_gastrocnemius_bulk_lomo",
                "sample_id": "rr5_f1",
                "mission": "RR-5",
                "true_label": "Flight",
                "predicted_label": "Flight",
                "flight_probability": 0.88,
                "embedding_0": 10.0,
                "embedding_1": 0.0,
            },
            {
                "task_id": "A2_gastrocnemius_bulk_lomo",
                "sample_id": "rr5_g1",
                "mission": "RR-5",
                "true_label": "Ground",
                "predicted_label": "Ground",
                "flight_probability": 0.12,
                "embedding_0": 10.0,
                "embedding_1": 0.2,
            },
        ]

    def test_submission_validator_accepts_valid_prediction_csv(self):
        from spacebio_bench import load_task_manifest
        from spacebio_bench.submissions import validate_submission

        manifest = load_task_manifest(
            REPO_ROOT / "v9" / "task_manifests" / "A2_gastrocnemius_bulk_lomo.json"
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            submission = Path(tmpdir) / "submission.csv"
            self.write_submission(submission, self.synthetic_submission_rows())
            report = validate_submission(manifest, submission)

        self.assertTrue(report.ok)
        self.assertEqual(report.n_rows, 4)
        self.assertEqual(report.embedding_columns, ["embedding_0", "embedding_1"])

    def test_submission_validator_rejects_missing_required_column(self):
        from spacebio_bench import load_task_manifest
        from spacebio_bench.submissions import validate_submission

        manifest = load_task_manifest(
            REPO_ROOT / "v9" / "task_manifests" / "A2_gastrocnemius_bulk_lomo.json"
        )
        rows = self.synthetic_submission_rows()
        for row in rows:
            row.pop("predicted_label")

        with tempfile.TemporaryDirectory() as tmpdir:
            submission = Path(tmpdir) / "submission.csv"
            self.write_submission(submission, rows)
            report = validate_submission(manifest, submission)

        self.assertFalse(report.ok)
        self.assertTrue(any("missing required columns" in error for error in report.errors))

    def test_submission_validator_rejects_empty_file(self):
        from spacebio_bench import load_task_manifest
        from spacebio_bench.submissions import validate_submission

        manifest = load_task_manifest(
            REPO_ROOT / "v9" / "task_manifests" / "A2_gastrocnemius_bulk_lomo.json"
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            submission = Path(tmpdir) / "submission.csv"
            submission.write_text("sample_id,true_label,predicted_label\n")
            report = validate_submission(manifest, submission)

        self.assertFalse(report.ok)
        self.assertIn("submission file has no prediction rows", report.errors)

    def test_submission_validator_rejects_task_id_mismatch(self):
        from spacebio_bench import load_task_manifest
        from spacebio_bench.submissions import validate_submission

        manifest = load_task_manifest(
            REPO_ROOT / "v9" / "task_manifests" / "A2_gastrocnemius_bulk_lomo.json"
        )
        rows = self.synthetic_submission_rows()
        rows[0]["task_id"] = "A1_liver_bulk_lomo"
        with tempfile.TemporaryDirectory() as tmpdir:
            submission = Path(tmpdir) / "submission.csv"
            self.write_submission(submission, rows)
            report = validate_submission(manifest, submission)

        self.assertFalse(report.ok)
        self.assertTrue(any("does not match manifest task_id" in error for error in report.errors))

    def test_evaluator_computes_genelab_minimal_metrics(self):
        from spacebio_bench import evaluate_submission, load_task_manifest

        manifest = load_task_manifest(
            REPO_ROOT / "v9" / "task_manifests" / "A2_gastrocnemius_bulk_lomo.json"
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            submission = Path(tmpdir) / "submission.csv"
            self.write_submission(submission, self.synthetic_submission_rows())
            result = evaluate_submission(manifest, submission)

        self.assertEqual(result["status"], "evaluated")
        self.assertAlmostEqual(result["metrics"]["macro_f1"]["value"], 1.0)
        self.assertAlmostEqual(result["metrics"]["balanced_accuracy"]["value"], 1.0)
        self.assertAlmostEqual(result["metrics"]["auroc"]["value"], 1.0)
        self.assertEqual(
            result["metrics"]["mission_discrimination"]["status"],
            "computed",
        )

    def test_evaluator_skips_optional_metrics_when_optional_columns_missing(self):
        from spacebio_bench import evaluate_submission, load_task_manifest

        manifest = load_task_manifest(
            REPO_ROOT / "v9" / "task_manifests" / "A2_gastrocnemius_bulk_lomo.json"
        )
        rows = []
        for row in self.synthetic_submission_rows():
            rows.append(
                {
                    "sample_id": row["sample_id"],
                    "true_label": row["true_label"],
                    "predicted_label": row["predicted_label"],
                }
            )

        with tempfile.TemporaryDirectory() as tmpdir:
            submission = Path(tmpdir) / "submission.csv"
            self.write_submission(submission, rows)
            result = evaluate_submission(manifest, submission)

        self.assertEqual(result["status"], "evaluated")
        self.assertEqual(result["metrics"]["auroc"]["status"], "skipped")
        self.assertEqual(result["metrics"]["mission_discrimination"]["status"], "skipped")

    def test_evaluator_skips_mission_discrimination_for_single_mission_embeddings(self):
        from spacebio_bench import evaluate_submission

        manifest = {
            "task_id": "single_mission_demo",
            "output": {
                "required_columns": ["sample_id", "true_label", "predicted_label"],
                "label_domain": ["Ground", "LEO_or_ISS"],
                "positive_label": "LEO_or_ISS",
            },
        }
        rows = [
            {
                "sample_id": "s1",
                "mission": "SpaceX-19",
                "true_label": "Ground",
                "predicted_label": "Ground",
                "flight_probability": 0.2,
                "embedding_0": 0.0,
            },
            {
                "sample_id": "s2",
                "mission": "SpaceX-19",
                "true_label": "LEO_or_ISS",
                "predicted_label": "LEO_or_ISS",
                "flight_probability": 0.8,
                "embedding_0": 1.0,
            },
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            submission = Path(tmpdir) / "submission.csv"
            self.write_submission(submission, rows)
            result = evaluate_submission(manifest, submission)

        self.assertEqual(result["status"], "evaluated")
        self.assertEqual(result["positive_label"], "LEO_or_ISS")
        self.assertEqual(result["metrics"]["mission_discrimination"]["status"], "skipped")
        self.assertIn("at least two missions", result["metrics"]["mission_discrimination"]["reason"])

    def test_evaluation_report_writes_metrics_and_run_manifest(self):
        from spacebio_bench import (
            evaluate_submission,
            load_task_manifest,
            write_evaluation_report,
        )

        manifest_path = REPO_ROOT / "v9" / "task_manifests" / "A2_gastrocnemius_bulk_lomo.json"
        manifest = load_task_manifest(manifest_path)
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            submission = tmp / "submission.csv"
            self.write_submission(submission, self.synthetic_submission_rows())
            result = evaluate_submission(manifest, submission)
            outputs = write_evaluation_report(
                evaluation_result=result,
                task_manifest=manifest,
                task_manifest_path=manifest_path,
                submission_path=submission,
                output_dir=tmp / "report",
                command=["spacebio-bench-test"],
            )

            metrics = json.loads(outputs["metrics"].read_text())
            run_manifest = json.loads(outputs["run_manifest"].read_text())

        self.assertEqual(metrics["status"], "evaluated")
        self.assertEqual(run_manifest["run_type"], "spacebio_bench_evaluation")
        self.assertEqual(run_manifest["task_id"], "A2_gastrocnemius_bulk_lomo")
        self.assertEqual(len(run_manifest["input_files"]), 2)

    def test_evaluate_v9_submission_cli_writes_report_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            submission = tmp / "submission.csv"
            self.write_submission(submission, self.synthetic_submission_rows())
            output_dir = tmp / "out"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "evaluate_v9_submission.py"),
                    "--task-manifest",
                    str(
                        REPO_ROOT
                        / "v9"
                        / "task_manifests"
                        / "A2_gastrocnemius_bulk_lomo.json"
                    ),
                    "--submission",
                    str(submission),
                    "--output-dir",
                    str(output_dir),
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )

            self.assertIn("metrics.json", result.stdout)
            self.assertTrue((output_dir / "metrics.json").exists())
            self.assertTrue((output_dir / "run_manifest.json").exists())

    def test_bulk_task_loader_exposes_fold_paths_and_counts(self):
        from spacebio_bench import load_bulk_task

        task = load_bulk_task(
            "A2_gastrocnemius_bulk_lomo",
            repo_root=REPO_ROOT,
        )

        self.assertEqual(task.task_id, "A2_gastrocnemius_bulk_lomo")
        self.assertEqual(len(task.folds), 3)
        first = task.folds[0]
        self.assertEqual(first.test_mission, "RR-1")
        self.assertEqual(first.train_row_count, 20)
        self.assertEqual(first.test_row_count, 12)
        self.assertEqual(first.selected_gene_count, 15760)
        self.assertTrue(first.paths["train_X.csv"].exists())

    def test_all_bulk_task_loader_covers_generated_manifest_set(self):
        from spacebio_bench import load_all_bulk_tasks

        tasks = load_all_bulk_tasks(repo_root=REPO_ROOT)

        self.assertEqual(len(tasks), 8)
        self.assertEqual(sum(len(task.folds) for task in tasks), 33)

    def test_task_data_index_files_cover_all_bulk_folds(self):
        index_path = REPO_ROOT / "v9" / "task_data_index.csv"
        with index_path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))

        self.assertEqual(len(rows), 33)
        first = rows[0]
        self.assertEqual(first["task_id"], "A1_liver_bulk_lomo")
        self.assertEqual(first["train_y_rows"], first["n_train"])
        self.assertFalse(first["train_X"].startswith("/"))

    def test_build_v9_task_data_index_cli_writes_temp_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "build_v9_task_data_index.py"),
                    "--repo-root",
                    str(REPO_ROOT),
                    "--csv",
                    str(tmp / "task_data_index.csv"),
                    "--json",
                    str(tmp / "task_data_index.json"),
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )

            self.assertIn("task_data_index.csv", result.stdout)
            self.assertTrue((tmp / "task_data_index.csv").exists())
            self.assertTrue((tmp / "task_data_index.json").exists())

    def test_source_inventory_deduplicates_manifest_sources(self):
        from spacebio_bench import TaskRegistry, build_source_inventory

        registry = TaskRegistry.from_dir(REPO_ROOT / "v9" / "task_manifests")
        rows = build_source_inventory(registry)
        by_source = {row["source_id"]: row for row in rows}

        self.assertEqual(len(rows), 22)
        self.assertEqual(by_source["OSD-101"]["mission"], "RR-1")
        self.assertEqual(by_source["OSD-101"]["tissue"], "gastrocnemius")
        self.assertEqual(by_source["OSD-101"]["organism"], "Mus musculus")
        self.assertEqual(by_source["OSD-101"]["taxon_id"], "10090")
        self.assertEqual(by_source["OSD-101"]["material_type"], "animal_tissue")
        self.assertEqual(by_source["OSD-101"]["feature_namespace"], "mouse_gene")
        self.assertIn("A1_liver_bulk_lomo", by_source["OSD-686"]["task_ids"])
        self.assertIn("A1_liver_bulk_lomo_combat", by_source["OSD-686"]["task_ids"])
        self.assertNotIn("A1_liver_bulk_lomo_iss_only", by_source["OSD-686"]["task_ids"])
        self.assertEqual(by_source["OSD-686"]["checksum_status"], "legacy_task_source_unfrozen")

    def test_build_v9_source_inventory_cli_writes_temp_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "build_v9_source_inventory.py"),
                    "--manifest-dir",
                    str(REPO_ROOT / "v9" / "task_manifests"),
                    "--csv",
                    str(tmp / "source_inventory.csv"),
                    "--json",
                    str(tmp / "source_inventory.json"),
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )

            with (tmp / "source_inventory.csv").open(newline="") as handle:
                rows = list(csv.DictReader(handle))

        self.assertIn("source_inventory.csv", result.stdout)
        self.assertEqual(len(rows), 22)
        self.assertEqual(rows[0]["release_target"], "v9_alpha_public_bulk_candidate")
        self.assertIn("organism", rows[0])
        self.assertIn("feature_namespace", rows[0])

    def test_build_v9_extension_source_inventories_cli_writes_temp_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "build_v9_extension_source_inventories.py"),
                    "--human-organoid-csv",
                    str(tmp / "human_organoid.csv"),
                    "--human-organoid-json",
                    str(tmp / "human_organoid.json"),
                    "--multispecies-csv",
                    str(tmp / "multispecies.csv"),
                    "--multispecies-json",
                    str(tmp / "multispecies.json"),
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )

            with (tmp / "human_organoid.csv").open(newline="") as handle:
                human_rows = list(csv.DictReader(handle))
            with (tmp / "multispecies.csv").open(newline="") as handle:
                multispecies_rows = list(csv.DictReader(handle))

        self.assertIn("human_organoid.csv", result.stdout)
        self.assertEqual({row["source_id"] for row in human_rows}, {"OSD-863", "OSD-871"})
        self.assertEqual(human_rows[0]["organism"], "Homo sapiens")
        self.assertEqual(human_rows[0]["task_families"], "human_organoid_spaceflight")
        self.assertEqual(len(multispecies_rows), 3)
        self.assertIn("OSD-207", {row["source_id"] for row in multispecies_rows})
        self.assertIn("Arabidopsis thaliana", {row["organism"] for row in multispecies_rows})
        self.assertEqual(multispecies_rows[0]["release_target"], "v9_extension_draft_not_frozen")

    def test_human_organoid_task_manifest_builder_preserves_draft_boundaries(self):
        from spacebio_bench import (
            HUMAN_ORGANOID_DRAFT_SOURCES,
            build_human_organoid_task_manifest,
        )

        with (REPO_ROOT / "v9" / "human_organoid" / "sample_factors.draft.csv").open(
            newline=""
        ) as handle:
            sample_factor_rows = list(csv.DictReader(handle))
        with (REPO_ROOT / "v9" / "human_organoid" / "expression_matrix_audit.draft.csv").open(
            newline=""
        ) as handle:
            expression_matrix_audit_rows = list(csv.DictReader(handle))
        manifest = build_human_organoid_task_manifest(
            HUMAN_ORGANOID_DRAFT_SOURCES,
            sample_factor_rows=sample_factor_rows,
            expression_matrix_audit_rows=expression_matrix_audit_rows,
        )

        self.assertEqual(manifest["task_id"], "draft_human_organoid_spaceflight")
        self.assertEqual(manifest["task_family"], "human_organoid_spaceflight")
        self.assertEqual(manifest["release_status"], "draft_not_frozen")
        self.assertEqual(manifest["organism"], "Homo sapiens")
        self.assertEqual(manifest["feature_namespace"], "human_gene")
        self.assertEqual(
            [source["source_id"] for source in manifest["source_records"]],
            ["OSD-863", "OSD-871"],
        )
        self.assertEqual(
            manifest["split"]["status"],
            "draft_sample_count_backed_pending_baseline",
        )
        self.assertEqual(
            manifest["split"]["sample_factor_status"],
            "condition_factors_and_geo_donor_line_metadata_parsed",
        )
        self.assertEqual(manifest["split"]["n_samples"], 42)
        self.assertEqual(manifest["split"]["source_sample_rows"]["OSD-863"], 19)
        self.assertEqual(
            manifest["split"]["donor_or_line_distribution"],
            {"Subject1": 11, "Subject2": 8, "Subject3": 10, "Subject4": 13},
        )
        self.assertEqual(
            manifest["split"]["donor_metadata_status"],
            "parsed_from_geo_series_matrix",
        )
        self.assertEqual(len(manifest["split"]["candidate_folds"]), 4)
        self.assertEqual(len(manifest["split"]["donor_diagnostic_folds"]), 4)
        self.assertEqual(
            {fold["status"] for fold in manifest["split"]["donor_diagnostic_folds"]},
            {"metadata_backed_pilot_only_not_default"},
        )
        self.assertIn("donor_subject", manifest["split"]["blocking_factors"])
        self.assertEqual(
            manifest["split"]["expression_matrix_status"],
            "matrix_downloaded_sample_aligned",
        )
        expression_sources = manifest["split"]["expression_matrix_sources"]
        self.assertEqual(set(expression_sources), {"OSD-863", "OSD-871"})
        self.assertEqual(expression_sources["OSD-863"]["n_sample_columns"], 19)
        self.assertEqual(expression_sources["OSD-871"]["n_feature_rows"], 30269)
        self.assertEqual(
            manifest["provenance"]["expression_matrix_status"],
            "matrix_downloaded_sample_aligned",
        )
        self.assertIn("single-mission", manifest["limitations"][3])

    def test_build_v9_human_organoid_task_manifest_cli_writes_temp_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            source_json = tmp / "source_inventory.draft.json"
            source_json.write_text(
                (REPO_ROOT / "v9" / "human_organoid" / "source_inventory.draft.json").read_text()
            )
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "build_v9_human_organoid_task_manifest.py"),
                    "--source-inventory",
                    str(source_json),
                    "--sample-factor-table",
                    str(REPO_ROOT / "v9" / "human_organoid" / "sample_factors.draft.csv"),
                    "--expression-matrix-audit",
                    str(
                        REPO_ROOT
                        / "v9"
                        / "human_organoid"
                        / "expression_matrix_audit.draft.csv"
                    ),
                    "--output-dir",
                    str(tmp / "task_manifests"),
                    "--index-csv",
                    str(tmp / "task_manifest_index.draft.csv"),
                    "--index-json",
                    str(tmp / "task_manifest_index.draft.json"),
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            manifest_path = tmp / "task_manifests" / "draft_human_organoid_spaceflight.json"
            payload = json.loads(manifest_path.read_text())
            with (tmp / "task_manifest_index.draft.csv").open(newline="") as handle:
                index_rows = list(csv.DictReader(handle))

        self.assertIn("draft_human_organoid_spaceflight.json", result.stdout)
        self.assertEqual(payload["split"]["reported_public_geo_sample_count"], 42)
        self.assertEqual(payload["split"]["status"], "draft_sample_count_backed_pending_baseline")
        self.assertEqual(
            payload["provenance"]["sample_table_status"],
            "sample_table_parsed",
        )
        self.assertEqual(
            payload["provenance"]["sample_factor_status"],
            "condition_factors_and_geo_donor_line_metadata_parsed",
        )
        self.assertEqual(
            payload["provenance"]["donor_or_line_status"],
            "parsed_from_geo_series_matrix",
        )
        self.assertEqual(
            payload["provenance"]["expression_matrix_status"],
            "matrix_downloaded_sample_aligned",
        )
        self.assertEqual(
            payload["split"]["expression_matrix_status"],
            "matrix_downloaded_sample_aligned",
        )
        self.assertEqual(
            payload["provenance"]["source_checksum_status"],
            "checksum_manifests_parsed_payloads_not_hashed",
        )
        self.assertEqual(index_rows[0]["task_family"], "human_organoid_spaceflight")

    def test_human_organoid_checksum_audit_records_two_draft_sources(self):
        audit_path = REPO_ROOT / "v9" / "human_organoid" / "source_checksum_audit.draft.csv"
        with audit_path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))

        self.assertEqual({row["source_id"] for row in rows}, {"OSD-863", "OSD-871"})
        self.assertEqual({row["api_status"] for row in rows}, {"ok"})
        self.assertEqual({row["audit_status"] for row in rows}, {"checksum_manifest_parsed"})
        self.assertEqual({row["freeze_ready"] for row in rows}, {"false"})
        self.assertTrue(all(int(row["parsed_checksum_entries"]) > 0 for row in rows))

    def test_human_organoid_sample_table_audit_records_two_draft_sources(self):
        audit_path = REPO_ROOT / "v9" / "human_organoid" / "sample_table_audit.draft.csv"
        with audit_path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))

        self.assertEqual({row["source_id"] for row in rows}, {"OSD-863", "OSD-871"})
        self.assertEqual({row["api_status"] for row in rows}, {"ok"})
        self.assertEqual({row["audit_status"] for row in rows}, {"sample_table_parsed"})
        self.assertTrue(all(int(row["sample_rows"]) > 0 for row in rows))

    def test_human_organoid_sample_factors_record_all_public_samples(self):
        factor_path = REPO_ROOT / "v9" / "human_organoid" / "sample_factors.draft.csv"
        with factor_path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))

        self.assertEqual(len(rows), 42)
        self.assertEqual({row["parse_status"] for row in rows}, {"parsed"})
        self.assertEqual({row["true_label"] for row in rows}, {"Ground", "LEO_or_ISS"})
        self.assertEqual({row["microglia_condition"] for row in rows}, {"with_microglia", "without_microglia"})
        self.assertEqual(
            {row["organoid_type"] for row in rows},
            {"cortical_neural_organoid", "dopaminergic_neural_organoid"},
        )
        self.assertEqual({row["donor_or_line_id"] for row in rows}, {"Subject1", "Subject2", "Subject3", "Subject4"})
        self.assertEqual(
            {row["ipsc_line_id"] for row in rows},
            {"051121-01-MR-017", "AK003-01-MR-008", "UEC741iPS517", "HDF410iPS504"},
        )
        self.assertEqual(
            {row["geo_metadata_status"] for row in rows},
            {"parsed_from_geo_series_matrix"},
        )

    def test_human_organoid_geo_metadata_records_public_geo_series_factors(self):
        geo_path = REPO_ROOT / "v9" / "human_organoid" / "geo_sample_metadata.draft.csv"
        with geo_path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))

        self.assertEqual(len(rows), 42)
        self.assertEqual({row["geo_series_accession"] for row in rows}, {"GSE259421"})
        self.assertEqual({row["geo_pubmed_id"] for row in rows}, {"39441987"})
        self.assertEqual({row["donor_or_line_id"] for row in rows}, {"Subject1", "Subject2", "Subject3", "Subject4"})
        self.assertEqual(
            {row["geo_title_spaceflight_condition"] for row in rows},
            {"Ground", "LEO_or_ISS"},
        )
        self.assertTrue(all(row["geo_biosample_accession"].startswith("SAMN") for row in rows))

    def test_expression_matrix_inspector_counts_rows_and_samples(self):
        from spacebio_bench import inspect_expression_matrix

        inspected = inspect_expression_matrix("gene_id,s1,s2\nENSG1,1,2\n\nENSG2,3,4\n")

        self.assertEqual(inspected["n_feature_rows"], 2)
        self.assertEqual(inspected["n_matrix_columns"], 3)
        self.assertEqual(inspected["sample_columns"], ["s1", "s2"])

    def test_human_organoid_expression_matrix_audit_records_aligned_matrices(self):
        audit_path = REPO_ROOT / "v9" / "human_organoid" / "expression_matrix_audit.draft.csv"
        with audit_path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))

        self.assertEqual({row["source_id"] for row in rows}, {"OSD-863", "OSD-871"})
        self.assertEqual({row["api_status"] for row in rows}, {"ok"})
        self.assertEqual({row["audit_status"] for row in rows}, {"matrix_downloaded_sample_aligned"})
        for row in rows:
            self.assertEqual(row["downloaded_matrix_count"], "1")
            self.assertEqual(row["matching_sample_count"], row["expected_sample_count"])
            self.assertEqual(row["missing_samples"], "")
            self.assertEqual(row["extra_samples"], "")
            self.assertTrue(row["matrix_files"].endswith("_rna_seq_Normalized_Counts_GLbulkRNAseq.csv"))
            self.assertNotIn("rRNArm", row["matrix_files"])
            self.assertNotIn("Unnormalized", row["matrix_files"])
            self.assertTrue(int(row["n_feature_rows"]) > 30000)

    def test_human_organoid_loader_aligns_matrices_and_draft_folds(self):
        from spacebio_bench import load_human_organoid_task

        task = load_human_organoid_task(repo_root=REPO_ROOT)

        self.assertEqual(task.task_id, "draft_human_organoid_spaceflight")
        self.assertEqual(task.n_samples, 42)
        self.assertEqual(task.n_features, 27986)
        self.assertEqual(set(task.matrix_paths), {"OSD-863", "OSD-871"})
        self.assertEqual(len(task.folds), 4)
        self.assertEqual(len(task.diagnostic_folds), 4)
        self.assertEqual(task.folds[0].fold_id, "holdout_organoid_type_cortical_neural_organoid")
        self.assertEqual(
            task.diagnostic_folds[0].fold_id,
            "holdout_donor_or_line_subject1",
        )
        self.assertEqual(len(task.folds[0].train_sample_ids), 23)
        self.assertEqual(len(task.folds[0].test_sample_ids), 19)
        self.assertEqual(set(task.features.index), {row["sample_id"] for row in task.sample_factors})

    def test_human_organoid_nearest_centroid_baseline_writes_draft_report(self):
        from spacebio_bench import (
            OrganoidBaselineConfig,
            load_human_organoid_task,
            run_organoid_nearest_centroid,
        )

        task = load_human_organoid_task(repo_root=REPO_ROOT)
        with tempfile.TemporaryDirectory() as tmpdir:
            result = run_organoid_nearest_centroid(
                task,
                output_dir=Path(tmpdir) / "organoid",
                config=OrganoidBaselineConfig(top_variable_genes=100),
                task_manifest_path=(
                    REPO_ROOT
                    / "v9"
                    / "human_organoid"
                    / "task_manifests"
                    / "draft_human_organoid_spaceflight.json"
                ),
                command=["test-organoid-baseline"],
            )
            with result.predictions_path.open(newline="") as handle:
                predictions = list(csv.DictReader(handle))

        self.assertEqual(result.n_predictions, 84)
        self.assertEqual(len(predictions), 84)
        self.assertEqual(result.evaluation["status"], "evaluated")
        self.assertEqual(result.evaluation["positive_label"], "LEO_or_ISS")
        self.assertEqual(result.evaluation["baseline"]["release_status"], "draft_not_frozen")
        self.assertEqual(
            result.evaluation["metrics"]["mission_discrimination"]["status"],
            "skipped",
        )
        self.assertEqual(
            result.evaluation["baseline"]["fold_family"],
            "default_sample_count_backed_folds",
        )

    def test_human_organoid_donor_diagnostic_baseline_writes_draft_report(self):
        from spacebio_bench import (
            OrganoidBaselineConfig,
            run_human_organoid_donor_diagnostics,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            result, summary = run_human_organoid_donor_diagnostics(
                repo_root=REPO_ROOT,
                output_root=Path(tmpdir) / "donor_diagnostics",
                config=OrganoidBaselineConfig(top_variable_genes=100),
                command=["test-organoid-donor-diagnostics"],
            )
            with result.predictions_path.open(newline="") as handle:
                predictions = list(csv.DictReader(handle))
            with summary["csv"].open(newline="") as handle:
                summary_rows = list(csv.DictReader(handle))

        self.assertEqual(result.n_predictions, 42)
        self.assertEqual(len(predictions), 42)
        self.assertEqual(result.evaluation["status"], "evaluated")
        self.assertEqual(
            result.evaluation["baseline"]["fold_family"],
            "geo_donor_or_line_diagnostic_folds",
        )
        self.assertEqual(
            result.evaluation["baseline"]["claim_boundary"],
            "donor_diagnostic_only_not_leaderboard",
        )
        self.assertEqual(
            {row["fold_id"] for row in predictions},
            {
                "holdout_donor_or_line_subject1",
                "holdout_donor_or_line_subject2",
                "holdout_donor_or_line_subject3",
                "holdout_donor_or_line_subject4",
            },
        )
        self.assertEqual(summary_rows[0]["fold_family"], "geo_donor_or_line_diagnostic_folds")
        self.assertEqual(summary_rows[0]["claim_boundary"], "donor_diagnostic_only_not_leaderboard")

    def test_run_v9_human_organoid_donor_diagnostics_cli_writes_summary(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "donor_diagnostics"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "run_v9_human_organoid_donor_diagnostics.py"),
                    "--repo-root",
                    str(REPO_ROOT),
                    "--output-dir",
                    str(output_dir),
                    "--top-variable-genes",
                    "50",
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            with (output_dir / "human_organoid_baseline_summary.csv").open(
                newline=""
            ) as handle:
                rows = list(csv.DictReader(handle))

        self.assertIn("human_organoid_baseline_summary.csv", result.stdout)
        self.assertEqual(rows[0]["fold_family"], "geo_donor_or_line_diagnostic_folds")
        self.assertEqual(rows[0]["claim_boundary"], "donor_diagnostic_only_not_leaderboard")

    def test_human_organoid_sensitivity_grid_writes_variant_summary(self):
        from spacebio_bench import (
            OrganoidBaselineConfig,
            run_human_organoid_sensitivity_grid,
        )

        configs = [
            OrganoidBaselineConfig(transform="log1p", scaling="zscore", top_variable_genes=50),
            OrganoidBaselineConfig(transform="none", scaling="none", top_variable_genes=50),
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            _, summary = run_human_organoid_sensitivity_grid(
                repo_root=REPO_ROOT,
                output_root=Path(tmpdir) / "sensitivity",
                configs=configs,
                command=["test-organoid-sensitivity"],
            )
            with summary["csv"].open(newline="") as handle:
                rows = list(csv.DictReader(handle))

        self.assertEqual(len(rows), 2)
        self.assertEqual(
            {row["variant_id"] for row in rows},
            {"tvg50_log1p_zscore", "tvg50_none_none"},
        )
        self.assertTrue(all(row["status"] == "evaluated" for row in rows))
        self.assertTrue(all(row["release_status"] == "draft_not_frozen" for row in rows))

    def test_human_organoid_diagnostics_record_scale_and_donor_gaps(self):
        from spacebio_bench import (
            build_organoid_donor_metadata_audit,
            build_organoid_group_diagnostics,
            build_organoid_sample_diagnostics,
            load_human_organoid_task,
        )

        task = load_human_organoid_task(repo_root=REPO_ROOT)
        sample_rows = build_organoid_sample_diagnostics(task)
        group_rows = build_organoid_group_diagnostics(sample_rows)
        with (REPO_ROOT / "v9" / "human_organoid" / "sample_table_audit.draft.csv").open(
            newline=""
        ) as handle:
            donor_rows = build_organoid_donor_metadata_audit(
                list(csv.DictReader(handle)),
                sample_factor_rows=task.sample_factors,
            )

        self.assertEqual(len(sample_rows), 42)
        self.assertTrue(all(row["n_features"] == "27986" for row in sample_rows))
        self.assertEqual(
            {row["donor_or_line_id"] for row in sample_rows},
            {"Subject1", "Subject2", "Subject3", "Subject4"},
        )
        self.assertTrue(all(0.0 <= float(row["zero_fraction"]) <= 1.0 for row in sample_rows))
        label_groups = {
            row["group_value"]: row
            for row in group_rows
            if row["group_field"] == "true_label"
        }
        self.assertEqual(set(label_groups), {"Ground", "LEO_or_ISS"})
        donor_groups = {
            row["group_value"]: row
            for row in group_rows
            if row["group_field"] == "donor_or_line_id"
        }
        self.assertEqual(set(donor_groups), {"Subject1", "Subject2", "Subject3", "Subject4"})
        self.assertEqual({row["donor_field_status"] for row in donor_rows}, {"recovered_from_geo_series_matrix"})
        self.assertEqual({row["blocking_status"] for row in donor_rows}, {"available_but_confounded_pilot_only"})

    def test_parse_checksum_manifest_extracts_common_line_formats(self):
        from spacebio_bench import parse_checksum_manifest

        entries = parse_checksum_manifest(
            "0123456789abcdef0123456789abcdef  sample.fastq.gz\n"
            "MD5 (other.fastq.gz) = fedcba9876543210fedcba9876543210\n"
            "third.fastq.gz\t00112233445566778899aabbccddeeff\n"
            "Merged Sequence Data\tfourth.fastq.gz\tffeeddccbbaa99887766554433221100\n"
            "# ignored comment\n"
            "not a checksum line\n"
        )

        self.assertEqual(len(entries), 4)
        self.assertEqual(entries[0]["algorithm"], "md5")
        self.assertEqual(entries[0]["filename"], "sample.fastq.gz")
        self.assertEqual(entries[1]["filename"], "other.fastq.gz")
        self.assertEqual(entries[2]["filename"], "third.fastq.gz")
        self.assertEqual(entries[3]["filename"], "fourth.fastq.gz")

    def test_source_checksum_audit_parses_mock_osdr_listing_and_manifest(self):
        from spacebio_bench.source_audit import FetchResult, audit_source_row

        api_url = "https://example.test/api/dataset/OSD-TEST/files/"
        manifest_url = "https://example.test/download/check.md5sum.txt"
        listing = {
            "OSD-TEST": {
                "files": {
                    "GLDS-TEST.md5sum.txt": {"REST_URL": manifest_url},
                    "GLDS-TEST_rna_seq_sample.fastq.gz": {
                        "REST_URL": "https://example.test/download/sample.fastq.gz"
                    },
                }
            }
        }

        def fake_fetcher(url: str) -> FetchResult:
            if url == api_url:
                return FetchResult(ok=True, url=url, body=json.dumps(listing).encode())
            if url == manifest_url:
                return FetchResult(
                    ok=True,
                    url=url,
                    body=b"0123456789abcdef0123456789abcdef  sample.fastq.gz\n",
                )
            return FetchResult(ok=False, url=url, error="unexpected url")

        row = {
            "source_id": "OSD-TEST",
            "glds_prefix": "GLDS-TEST",
            "mission": "RR-X",
            "tissue": "liver",
            "task_ids": "demo_task",
        }
        audit = audit_source_row(
            row,
            fetcher=fake_fetcher,
            api_base="https://example.test/api",
        )

        self.assertEqual(audit["api_status"], "ok")
        self.assertEqual(audit["n_files"], "2")
        self.assertEqual(audit["checksum_manifest_count"], "1")
        self.assertEqual(audit["parsed_checksum_entries"], "1")
        self.assertEqual(audit["checksum_algorithms"], "md5")
        self.assertEqual(audit["checksum_payload_matches"], "1")
        self.assertEqual(audit["audit_status"], "checksum_manifest_parsed")
        self.assertEqual(audit["freeze_ready"], "false")

    def test_source_checksum_audit_records_missing_manifest_status(self):
        from spacebio_bench.source_audit import FetchResult, audit_source_row

        listing = {
            "OSD-NOCHECK": {
                "files": {
                    "sample.fastq.gz": {"REST_URL": "https://example.test/download/sample.fastq.gz"}
                }
            }
        }

        def fake_fetcher(url: str) -> FetchResult:
            return FetchResult(ok=True, url=url, body=json.dumps(listing).encode())

        audit = audit_source_row(
            {"source_id": "OSD-NOCHECK"},
            fetcher=fake_fetcher,
            api_base="https://example.test/api",
        )

        self.assertEqual(audit["api_status"], "ok")
        self.assertEqual(audit["checksum_manifest_count"], "0")
        self.assertEqual(audit["audit_status"], "no_checksum_manifest_listed")

    def test_parse_organoid_contrast_table_identifies_spaceflight_pairs(self):
        from spacebio_bench import parse_organoid_contrast_table

        contrast_text = (
            '"","(control & Ground Control & with microglia)v'
            '(control & Space Flight & with microglia)",'
            '"(control & Space Flight & with microglia)v'
            '(control & Ground Control & with microglia)",'
            '"(control & Ground Control & with microglia)v'
            '(disease & Ground Control & with microglia)"\n'
            '"1","control...Ground.Control...with.microglia",'
            '"control...Space.Flight...with.microglia",'
            '"control...Ground.Control...with.microglia"\n'
            '"2","control...Space.Flight...with.microglia",'
            '"control...Ground.Control...with.microglia",'
            '"disease...Ground.Control...with.microglia"\n'
        )

        parsed = parse_organoid_contrast_table(contrast_text)

        self.assertEqual(parsed["parse_status"], "parsed")
        self.assertEqual(parsed["contrast_pair_count"], 3)
        self.assertEqual(len(parsed["direct_spaceflight_contrasts"]), 1)
        self.assertEqual(len(parsed["reversed_spaceflight_contrasts"]), 1)
        self.assertEqual(set(parsed["disease_contexts"]), {"control", "disease"})

    def test_organoid_signature_reference_audit_parses_mock_osdr_listing(self):
        from spacebio_bench.organoid_signature_audit import (
            audit_organoid_signature_reference_row,
        )
        from spacebio_bench.source_audit import FetchResult

        api_url = "https://example.test/api/dataset/OSD-SIG/files/"
        contrast_url = "https://example.test/download/contrasts.csv"
        listing = {
            "OSD-SIG": {
                "files": {
                    "GLDS-SIG_rna_seq_differential_expression_GLbulkRNAseq.csv": {
                        "REST_URL": "https://example.test/download/de.csv",
                        "metadata": {
                            "data_type": "differential expression table",
                            "file_size": 12345,
                        },
                    },
                    "GLDS-SIG_rna_seq_contrasts_GLbulkRNAseq.csv": {
                        "REST_URL": contrast_url,
                    },
                    "GLDS-SIG_rna_seq_Normalized_Counts_GLbulkRNAseq.csv": {
                        "REST_URL": "https://example.test/download/counts.csv",
                    },
                }
            }
        }
        contrast_text = (
            '"","(control & Ground Control & with microglia)v'
            '(control & Space Flight & with microglia)",'
            '"(control & Space Flight & with microglia)v'
            '(control & Ground Control & with microglia)"\n'
            '"1","control...Ground.Control...with.microglia",'
            '"control...Space.Flight...with.microglia"\n'
            '"2","control...Space.Flight...with.microglia",'
            '"control...Ground.Control...with.microglia"\n'
        )

        def fake_fetcher(url: str) -> FetchResult:
            if url == api_url:
                return FetchResult(ok=True, url=url, body=json.dumps(listing).encode())
            if url == contrast_url:
                return FetchResult(ok=True, url=url, body=contrast_text.encode())
            return FetchResult(ok=False, url=url, error="unexpected url")

        audit = audit_organoid_signature_reference_row(
            {"source_id": "OSD-SIG", "glds_prefix": "GLDS-SIG"},
            fetcher=fake_fetcher,
            api_base="https://example.test/api",
        )

        self.assertEqual(audit["api_status"], "ok")
        self.assertEqual(audit["de_reference_file_count"], "1")
        self.assertEqual(audit["contrast_file_count"], "1")
        self.assertEqual(audit["contrast_pair_count"], "2")
        self.assertEqual(audit["direct_spaceflight_contrast_count"], "1")
        self.assertEqual(audit["reversed_spaceflight_contrast_count"], "1")
        self.assertEqual(
            audit["metric_reference_status"],
            "public_osdr_de_reference_tables_available_pending_contrast_freeze",
        )
        self.assertEqual(
            audit["recommended_metric_policy"],
            "keep_classification_primary_enable_de_signature_after_frozen_contrast_subset",
        )

    def test_human_organoid_task_manifest_records_signature_reference_policy(self):
        from spacebio_bench import (
            HUMAN_ORGANOID_DRAFT_SOURCES,
            build_human_organoid_task_manifest,
        )

        signature_rows = [
            {
                "source_id": "OSD-863",
                "audit_status": "reference_tables_listed_contrast_definitions_parsed",
                "de_reference_files": "GLDS-716_rna_seq_differential_expression_GLbulkRNAseq.csv",
                "de_reference_file_sizes": "141841931",
                "contrast_files": "GLDS-716_rna_seq_contrasts_GLbulkRNAseq.csv",
                "contrast_pair_count": "56",
                "direct_spaceflight_contrast_count": "4",
                "reversed_spaceflight_contrast_count": "4",
                "metric_reference_status": (
                    "public_osdr_de_reference_tables_available_pending_contrast_freeze"
                ),
                "recommended_metric_policy": (
                    "keep_classification_primary_enable_de_signature_after_frozen_contrast_subset"
                ),
            },
            {
                "source_id": "OSD-871",
                "audit_status": "reference_tables_listed_contrast_definitions_parsed",
                "de_reference_files": "GLDS-720_rna_seq_differential_expression_GLbulkRNAseq.csv",
                "de_reference_file_sizes": "143992071",
                "contrast_files": "GLDS-720_rna_seq_contrasts_GLbulkRNAseq.csv",
                "contrast_pair_count": "56",
                "direct_spaceflight_contrast_count": "4",
                "reversed_spaceflight_contrast_count": "4",
                "metric_reference_status": (
                    "public_osdr_de_reference_tables_available_pending_contrast_freeze"
                ),
                "recommended_metric_policy": (
                    "keep_classification_primary_enable_de_signature_after_frozen_contrast_subset"
                ),
            },
        ]

        manifest = build_human_organoid_task_manifest(
            HUMAN_ORGANOID_DRAFT_SOURCES,
            signature_reference_audit_rows=signature_rows,
        )

        self.assertEqual(
            manifest["provenance"]["signature_reference_status"],
            "public_osdr_de_reference_tables_available_pending_contrast_freeze",
        )
        self.assertEqual(
            manifest["reference_signatures"]["signature_reference_sources"]["OSD-863"][
                "direct_spaceflight_contrast_count"
            ],
            "4",
        )

    def test_sample_table_audit_parses_mock_osdr_listing_and_table(self):
        from spacebio_bench.sample_table_audit import audit_sample_table_row
        from spacebio_bench.source_audit import FetchResult

        api_url = "https://example.test/api/dataset/OSD-SAMPLE/files/"
        table_url = "https://example.test/download/SampleTable.csv"
        listing = {
            "OSD-SAMPLE": {
                "files": {
                    "GLDS-SAMPLE_rna_seq_SampleTable_GLbulkRNAseq.csv": {
                        "REST_URL": table_url
                    }
                }
            }
        }

        def fake_fetcher(url: str) -> FetchResult:
            if url == api_url:
                return FetchResult(ok=True, url=url, body=json.dumps(listing).encode())
            if url == table_url:
                return FetchResult(
                    ok=True,
                    url=url,
                    body=(
                        "Sample Name,Spaceflight,Organoid Type,Microglia Condition\n"
                        "s1,LEO,cortical,with\n"
                        "s2,Ground,cortical,without\n"
                    ).encode(),
                )
            return FetchResult(ok=False, url=url, error="unexpected url")

        audit = audit_sample_table_row(
            {"source_id": "OSD-SAMPLE", "glds_prefix": "GLDS-SAMPLE"},
            fetcher=fake_fetcher,
            api_base="https://example.test/api",
        )

        self.assertEqual(audit["api_status"], "ok")
        self.assertEqual(audit["sample_table_file_count"], "1")
        self.assertEqual(audit["sample_rows"], "2")
        self.assertEqual(audit["audit_status"], "sample_table_parsed")
        self.assertIn("Microglia Condition", audit["candidate_condition_columns"])

    def test_parse_condition_factors_maps_organoid_condition_tokens(self):
        from spacebio_bench import parse_condition_factors

        factors = parse_condition_factors(
            "Sporadic.Parkinson.disease...Space.Flight...with.microglia"
        )

        self.assertEqual(factors["disease_context"], "sporadic_parkinson_disease")
        self.assertEqual(factors["spaceflight_condition"], "LEO_or_ISS")
        self.assertEqual(factors["true_label"], "LEO_or_ISS")
        self.assertEqual(factors["microglia_condition"], "with_microglia")
        self.assertEqual(factors["parse_status"], "parsed")

    def test_parse_geo_series_matrix_recovers_subject_and_ipsc_line(self):
        from spacebio_bench import parse_geo_series_matrix

        text = (
            '!Series_geo_accession\t"GSE259421"\n'
            '!Series_pubmed_id\t"39441987"\n'
            '!Sample_title\t'
            '"Subject1, Cortical organoid, without microglia, Ground1"\t'
            '"Subject4, Dopaminergic organoid,  with microglia , LEO2"\n'
            '!Sample_geo_accession\t"GSM1"\t"GSM2"\n'
            '!Sample_source_name_ch1\t"051121-01-MR-017"\t"HDF410iPS504"\n'
            '!Sample_characteristics_ch1\t'
            '"cell line: 051121-01-MR-017"\t"cell line: HDF410iPS504"\n'
            '!Sample_characteristics_ch1\t'
            '"cell type: iPSC-derived cortical organoid without microglia"\t'
            '"cell type: iPSC-derived dopaminergic organoid with microglia"\n'
            '!Sample_characteristics_ch1\t'
            '"treatment: ground control, no microgravity exposure"\t'
            '"treatment: 38 days microgravity exposure"\n'
            '!Sample_relation\t'
            '"BioSample: https://www.ncbi.nlm.nih.gov/biosample/SAMN1"\t'
            '"BioSample: https://www.ncbi.nlm.nih.gov/biosample/SAMN2"\n'
            '!Sample_relation\t'
            '"SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRX1"\t'
            '"SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRX2"\n'
            "!series_matrix_table_begin\n"
        )

        rows = parse_geo_series_matrix(text)

        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["donor_or_line_id"], "Subject1")
        self.assertEqual(rows[0]["ipsc_line_id"], "051121-01-MR-017")
        self.assertEqual(rows[0]["geo_title_spaceflight_condition"], "Ground")
        self.assertEqual(rows[1]["donor_or_line_id"], "Subject4")
        self.assertEqual(rows[1]["geo_title_microglia_condition"], "with_microglia")
        self.assertEqual(rows[1]["geo_title_spaceflight_condition"], "LEO_or_ISS")
        self.assertEqual(rows[1]["geo_biosample_accession"], "SAMN2")
        self.assertEqual(rows[1]["geo_sra_accession"], "SRX2")

    def test_build_v9_human_organoid_geo_metadata_cli_writes_temp_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            matrix = tmp / "GSE259421_series_matrix.txt.gz"
            matrix_text = (
                '!Series_geo_accession\t"GSE259421"\n'
                '!Sample_title\t'
                '"Subject1, Cortical organoid, without microglia, Ground1"\t'
                '"Subject1, Cortical organoid, with microglia, LEO1"\n'
                '!Sample_geo_accession\t"GSM1"\t"GSM2"\n'
                '!Sample_source_name_ch1\t"line-a"\t"line-a"\n'
                '!Sample_characteristics_ch1\t"cell line: line-a"\t"cell line: line-a"\n'
                '!Sample_characteristics_ch1\t'
                '"cell type: iPSC-derived cortical organoid without microglia"\t'
                '"cell type: iPSC-derived cortical organoid with microglia"\n'
                '!Sample_characteristics_ch1\t'
                '"treatment: ground control, no microgravity exposure"\t'
                '"treatment: 38 days microgravity exposure"\n'
                "!series_matrix_table_begin\n"
            )
            with gzip.open(matrix, "wt", encoding="utf-8") as handle:
                handle.write(matrix_text)
            sample_factors = tmp / "sample_factors.csv"
            sample_factors.write_text(
                "source_id,glds_prefix,sample_id,raw_condition,disease_context,"
                "spaceflight_condition,true_label,microglia_condition,organoid_type,"
                "sample_table_file,parse_status\n"
                "OSD-X,GLDS-X,GSM1,a,b,Ground,Ground,without_microglia,"
                "cortical_neural_organoid,table.csv,parsed\n"
                "OSD-X,GLDS-X,GSM2,a,b,LEO_or_ISS,LEO_or_ISS,with_microglia,"
                "cortical_neural_organoid,table.csv,parsed\n"
            )
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "build_v9_human_organoid_geo_metadata.py"),
                    "--series-matrix-gz",
                    str(matrix),
                    "--sample-factors",
                    str(sample_factors),
                    "--geo-csv",
                    str(tmp / "geo.csv"),
                    "--geo-json",
                    str(tmp / "geo.json"),
                    "--merged-csv",
                    str(tmp / "merged.csv"),
                    "--merged-json",
                    str(tmp / "merged.json"),
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            with (tmp / "merged.csv").open(newline="") as handle:
                merged = list(csv.DictReader(handle))

        self.assertIn("geo.csv", result.stdout)
        self.assertEqual({row["donor_or_line_id"] for row in merged}, {"Subject1"})
        self.assertEqual({row["ipsc_line_id"] for row in merged}, {"line-a"})
        self.assertEqual(
            {row["geo_metadata_status"] for row in merged},
            {"parsed_from_geo_series_matrix"},
        )

    def test_write_source_checksum_audit_writes_temp_outputs(self):
        from spacebio_bench import write_source_checksum_audit

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            csv_path, json_path = write_source_checksum_audit(
                [
                    {
                        "source_id": "OSD-TEST",
                        "api_status": "ok",
                        "audit_status": "checksum_manifest_parsed",
                    }
                ],
                csv_path=tmp / "source_checksum_audit.csv",
                json_path=tmp / "source_checksum_audit.json",
            )
            with csv_path.open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            payload = json.loads(json_path.read_text())

        self.assertEqual(rows[0]["source_id"], "OSD-TEST")
        self.assertEqual(payload[0]["audit_status"], "checksum_manifest_parsed")

    def test_v9_datapackage_draft_describes_current_public_bulk_resources(self):
        from spacebio_bench import build_v9_public_bulk_datapackage

        package = build_v9_public_bulk_datapackage(repo_root=REPO_ROOT)
        resources = {resource["name"]: resource for resource in package["resources"]}

        self.assertEqual(package["profile"], "data-package")
        self.assertEqual(package["spacebio_bench:release_status"], "draft_not_frozen")
        self.assertIn("source_checksum_audit", resources)
        self.assertIn("baseline_predictions", resources)
        self.assertEqual(resources["task_manifests"]["spacebio_bench:file_count"], 8)
        self.assertEqual(resources["baseline_predictions"]["spacebio_bench:file_count"], 24)
        self.assertEqual(
            resources["source_checksum_audit"]["schema"]["primaryKey"],
            ["source_id"],
        )
        for resource in package["resources"]:
            paths = resource["path"] if isinstance(resource["path"], list) else [resource["path"]]
            self.assertTrue(paths)
            for path in paths:
                self.assertFalse(path.startswith("/"))
                self.assertFalse(path.startswith("../"))
                self.assertNotIn("submissions", path)

    def test_build_v9_datapackage_draft_cli_writes_temp_descriptor(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "build_v9_datapackage_draft.py"),
                    "--repo-root",
                    str(REPO_ROOT),
                    "--json",
                    str(tmp / "datapackage.draft.json"),
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            payload = json.loads((tmp / "datapackage.draft.json").read_text())

        self.assertIn("datapackage.draft.json", result.stdout)
        self.assertEqual(payload["name"], "spacebio-bench-v9-public-bulk-draft")
        self.assertEqual(len(payload["resources"]), 11)

    def test_v9_hf_dataset_card_records_draft_boundaries(self):
        card = (REPO_ROOT / "docs" / "v9_hf_dataset_card.md").read_text()

        required_phrases = [
            "Release status: draft, not frozen.",
            "payload-level hash verification",
            "22 deduplicated public OSDR source rows",
            "8 generated public bulk LOMO task manifests",
            "draft_not_frozen",
            "checksum_manifests_parsed_payloads_not_hashed",
            "controlled-access human sequence data",
            "NASA Open Science Data Repository",
        ]
        for phrase in required_phrases:
            self.assertIn(phrase, card)
        self.assertIn("Out-of-scope uses", card)

    def test_nearest_centroid_predicts_real_a2_fold(self):
        from spacebio_bench import NearestCentroidConfig, load_bulk_task
        from spacebio_bench.baselines import predict_fold

        task = load_bulk_task("A2_gastrocnemius_bulk_lomo", repo_root=REPO_ROOT)
        result = predict_fold(task.folds[0], config=NearestCentroidConfig())

        self.assertEqual(result.fold_id, "fold_RR-1_test")
        self.assertEqual(result.n_test, 12)
        self.assertEqual(len(result.predictions), 12)
        self.assertEqual(
            set(result.predictions[0]),
            {
                "task_id",
                "fold_id",
                "sample_id",
                "mission",
                "true_label",
                "predicted_label",
                "flight_probability",
                "embedding_0",
                "embedding_1",
            },
        )
        probabilities = [float(row["flight_probability"]) for row in result.predictions]
        self.assertTrue(all(0.0 <= value <= 1.0 for value in probabilities))
        self.assertTrue(
            set(row["predicted_label"] for row in result.predictions)
            <= {"Flight", "Ground"}
        )

    def test_nearest_centroid_task_runner_writes_evaluable_outputs(self):
        from spacebio_bench import NearestCentroidConfig, load_bulk_task, run_nearest_centroid_task

        task = load_bulk_task("A2_gastrocnemius_bulk_lomo", repo_root=REPO_ROOT)
        with tempfile.TemporaryDirectory() as tmpdir:
            result = run_nearest_centroid_task(
                task,
                output_dir=Path(tmpdir) / "nearest_centroid",
                config=NearestCentroidConfig(),
                task_manifest_path=(
                    REPO_ROOT
                    / "v9"
                    / "task_manifests"
                    / "A2_gastrocnemius_bulk_lomo.json"
                ),
                command=["spacebio-bench-test", "nearest-centroid"],
            )

            with result.predictions_path.open(newline="") as handle:
                predictions = list(csv.DictReader(handle))
            metrics = json.loads(result.metrics_path.read_text())
            run_manifest = json.loads(result.run_manifest_path.read_text())

        self.assertEqual(result.n_predictions, 32)
        self.assertEqual(len(predictions), 32)
        self.assertEqual(metrics["status"], "evaluated")
        self.assertEqual(metrics["baseline"]["baseline_id"], "nearest_centroid")
        self.assertEqual(run_manifest["status"], "evaluated")

    def test_run_v9_nearest_centroid_cli_writes_task_summary(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "reports"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "run_v9_nearest_centroid.py"),
                    "--task-id",
                    "A2_gastrocnemius_bulk_lomo",
                    "--repo-root",
                    str(REPO_ROOT),
                    "--output-dir",
                    str(output_dir),
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )

            with (output_dir / "bulk_lomo_summary.csv").open(newline="") as handle:
                rows = list(csv.DictReader(handle))

        self.assertIn("bulk_lomo_summary.csv", result.stdout)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["task_id"], "A2_gastrocnemius_bulk_lomo")
        self.assertEqual(rows[0]["status"], "evaluated")

    def test_sklearn_pca_lr_predicts_real_a2_fold(self):
        from spacebio_bench import SklearnBaselineConfig, load_bulk_task
        from spacebio_bench.baselines import predict_sklearn_fold

        task = load_bulk_task("A2_gastrocnemius_bulk_lomo", repo_root=REPO_ROOT)
        result = predict_sklearn_fold(
            task.folds[0],
            config=SklearnBaselineConfig(baseline_id="pca_lr"),
        )

        self.assertEqual(result.fold_id, "fold_RR-1_test")
        self.assertEqual(result.n_test, 12)
        self.assertEqual(len(result.predictions), 12)
        self.assertEqual(result.fit_details["pca_components"], 19)
        self.assertIn("embedding_0", result.predictions[0])
        probabilities = [float(row["flight_probability"]) for row in result.predictions]
        self.assertTrue(all(0.0 <= value <= 1.0 for value in probabilities))

    def test_sklearn_logistic_runner_writes_evaluable_outputs(self):
        from spacebio_bench import (
            SklearnBaselineConfig,
            load_bulk_task,
            run_sklearn_baseline_task,
        )

        task = load_bulk_task("A2_gastrocnemius_bulk_lomo", repo_root=REPO_ROOT)
        with tempfile.TemporaryDirectory() as tmpdir:
            result = run_sklearn_baseline_task(
                task,
                output_dir=Path(tmpdir) / "logistic",
                config=SklearnBaselineConfig(baseline_id="logistic_l2"),
                task_manifest_path=(
                    REPO_ROOT
                    / "v9"
                    / "task_manifests"
                    / "A2_gastrocnemius_bulk_lomo.json"
                ),
                command=["spacebio-bench-test", "logistic-l2"],
            )

            with result.predictions_path.open(newline="") as handle:
                predictions = list(csv.DictReader(handle))
            metrics = json.loads(result.metrics_path.read_text())

        self.assertEqual(result.n_predictions, 32)
        self.assertEqual(len(predictions), 32)
        self.assertEqual(metrics["status"], "evaluated")
        self.assertEqual(metrics["baseline"]["baseline_id"], "logistic_regression_l2")
        self.assertEqual(metrics["metrics"]["mission_discrimination"]["status"], "skipped")

    def test_run_v9_sklearn_baselines_cli_writes_task_summary(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "reports"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "run_v9_sklearn_baselines.py"),
                    "--baseline",
                    "pca_lr",
                    "--task-id",
                    "A2_gastrocnemius_bulk_lomo",
                    "--repo-root",
                    str(REPO_ROOT),
                    "--output-dir",
                    str(output_dir),
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )

            with (output_dir / "bulk_lomo_summary.csv").open(newline="") as handle:
                rows = list(csv.DictReader(handle))

        self.assertIn("bulk_lomo_summary.csv", result.stdout)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["baseline_id"], "pca_logistic_regression")
        self.assertEqual(rows[0]["task_id"], "A2_gastrocnemius_bulk_lomo")
        self.assertEqual(rows[0]["status"], "evaluated")

    def test_build_v9_baseline_summary_cli_combines_baseline_tables(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            nearest = tmp / "nearest.csv"
            sklearn = tmp / "sklearn.csv"
            nearest.write_text(
                "task_id,status,validation_ok,n_predictions,macro_f1,balanced_accuracy,"
                "auroc,calibration_error,mission_discrimination,output_dir,predictions,"
                "metrics,run_manifest,error\n"
                "A2_gastrocnemius_bulk_lomo,evaluated,True,32,0.33,0.5,0.7,0.1,1,"
                "nearest/out,nearest/pred.csv,nearest/metrics.json,nearest/run.json,\n"
            )
            sklearn.write_text(
                "baseline_id,task_id,status,validation_ok,n_predictions,macro_f1,"
                "balanced_accuracy,auroc,calibration_error,mission_discrimination,"
                "output_dir,predictions,metrics,run_manifest,error\n"
                "pca_logistic_regression,A2_gastrocnemius_bulk_lomo,evaluated,True,"
                "32,0.44,0.55,0.8,0.2,0.9,pca/out,pca/pred.csv,pca/metrics.json,"
                "pca/run.json,\n"
            )
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "build_v9_baseline_summary.py"),
                    "--nearest-centroid-summary",
                    str(nearest),
                    "--sklearn-summary",
                    str(sklearn),
                    "--csv",
                    str(tmp / "combined.csv"),
                    "--json",
                    str(tmp / "combined.json"),
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )

            with (tmp / "combined.csv").open(newline="") as handle:
                rows = list(csv.DictReader(handle))

        self.assertIn("combined.csv", result.stdout)
        self.assertEqual(len(rows), 2)
        self.assertEqual(
            {row["baseline_id"] for row in rows},
            {"nearest_centroid", "pca_logistic_regression"},
        )


if __name__ == "__main__":
    unittest.main()
