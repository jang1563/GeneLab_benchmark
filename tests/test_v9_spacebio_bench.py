import csv
import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def write_osd918_osdr_file_listing_fixture(path: Path) -> None:
    srxs = [
        "SRX28491856",
        "SRX28491916",
        "SRX28491975",
        "SRX28491923",
        "SRX28491934",
        "SRX28491993",
        "SRX28491866",
        "SRX28491877",
    ]
    filenames = [
        "OSD-918_metadata_OSD-918-ISA.zip",
        "GLDS-746_scRNA-Seq_raw_md5sum_GLscRNAseq.txt",
        "GLDS-746_scRNA_Seq_raw_multiqc_GLscRNAseq_report.zip",
    ]
    for srx in srxs:
        filenames.append(f"GLDS-746_scRNA-Seq_{srx}_R1_raw.fastq.gz")
        filenames.append(f"GLDS-746_scRNA-Seq_{srx}_R2_raw.fastq.gz")
    files = {
        filename: {
            "REST_URL": (
                "https://osdr.nasa.gov/geode-py/ws/studies/OSD-918/download?"
                f"source=datamanager&file={filename}"
            ),
            "URL": (
                "https://visualization.osdr.nasa.gov/biodata/files/OSD-918/"
                f"{filename}"
            ),
        }
        for filename in filenames
    }
    path.write_text(json.dumps({"OSD-918": {"files": files}}, indent=2) + "\n")


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
        self.assertIn("genelab_multispecies_pilot", METRIC_PROFILES)
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

    def test_run_human_organoid_response_signature_smoke_cli_writes_report(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            reference = tmp / "reference.csv.gz"
            with gzip.open(reference, "wt", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "source_id",
                        "contrast_id",
                        "gene_symbol",
                        "ensembl_id",
                        "log2fc_leo_or_iss_minus_ground",
                        "significant_fdr_0_05",
                    ],
                )
                writer.writeheader()
                writer.writerows(
                    [
                        {
                            "source_id": "OSD-X",
                            "contrast_id": "contrast_a",
                            "gene_symbol": "GENE1",
                            "ensembl_id": "ENSG1",
                            "log2fc_leo_or_iss_minus_ground": "1.0",
                            "significant_fdr_0_05": "true",
                        },
                        {
                            "source_id": "OSD-X",
                            "contrast_id": "contrast_a",
                            "gene_symbol": "GENE2",
                            "ensembl_id": "ENSG2",
                            "log2fc_leo_or_iss_minus_ground": "-1.0",
                            "significant_fdr_0_05": "true",
                        },
                    ]
                )
            output_dir = tmp / "response_signature_smoke"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "run_v9_human_organoid_response_signature_smoke.py"),
                    "--task-manifest",
                    str(
                        REPO_ROOT
                        / "v9"
                        / "human_organoid"
                        / "task_manifests"
                        / "draft_human_organoid_spaceflight.json"
                    ),
                    "--reference-signature-table",
                    str(reference),
                    "--output-dir",
                    str(output_dir),
                    "--significant-per-contrast",
                    "2",
                    "--background-per-contrast",
                    "0",
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            metrics = json.loads((output_dir / "metrics.json").read_text())
            run_manifest = json.loads((output_dir / "run_manifest.json").read_text())
            smoke_note = (output_dir / "README.md").read_text()

        self.assertIn("response_signature.csv", result.stdout)
        self.assertEqual(metrics["metrics"]["de_direction_match"]["status"], "computed")
        self.assertAlmostEqual(metrics["metrics"]["de_direction_match"]["value"], 1.0)
        self.assertEqual(
            metrics["metrics"]["signature_rank_correlation"]["status"],
            "computed",
        )
        self.assertEqual(
            {
                input_file["role"]
                for input_file in run_manifest["input_files"]
            },
            {
                "submission",
                "task_manifest",
                "response_signature",
                "reference_signature_table",
            },
        )
        self.assertIn("not a biological model result", smoke_note)

    def test_source_transfer_response_signature_adapter_excludes_target_source(self):
        from spacebio_bench import (
            REFERENCE_USAGE_POLICY,
            build_source_transfer_response_signature,
            load_human_organoid_task,
        )

        task = load_human_organoid_task(repo_root=REPO_ROOT)
        with tempfile.TemporaryDirectory() as tmpdir:
            manifest_path = Path(tmpdir) / "de_reference_manifest.json"
            manifest_path.write_text(
                json.dumps(
                    {
                        "sources": [
                            {
                                "source_id": "OSD-863",
                                "contrasts": [{"contrast_id": "osd_863_demo"}],
                            },
                            {
                                "source_id": "OSD-871",
                                "contrasts": [{"contrast_id": "osd_871_demo"}],
                            },
                        ]
                    }
                )
            )
            result = build_source_transfer_response_signature(
                task,
                de_reference_manifest_path=manifest_path,
            )

        self.assertEqual(result.signature_model_id, "organoid_source_transfer_empirical_signature")
        self.assertEqual(result.n_response_rows, task.n_features * 2)
        rows_by_source = {}
        for row in result.response_rows:
            rows_by_source.setdefault(row["source_id"], row)
        self.assertEqual(rows_by_source["OSD-863"]["training_source_id"], "OSD-871")
        self.assertEqual(rows_by_source["OSD-871"]["training_source_id"], "OSD-863")
        self.assertTrue(
            all(row["reference_usage_policy"] == REFERENCE_USAGE_POLICY for row in result.response_rows)
        )
        for summary in result.fold_summaries:
            self.assertNotEqual(summary["target_source_id"], summary["training_source_id"])

    def test_run_human_organoid_source_transfer_signature_cli_writes_report(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            output_dir = tmp / "source_transfer_signature"
            de_manifest = tmp / "de_reference_manifest.json"
            de_manifest.write_text(
                json.dumps(
                    {
                        "sources": [
                            {
                                "source_id": "OSD-863",
                                "contrasts": [{"contrast_id": "osd_863_demo"}],
                            },
                            {
                                "source_id": "OSD-871",
                                "contrasts": [{"contrast_id": "osd_871_demo"}],
                            },
                        ]
                    }
                )
            )
            reference = tmp / "reference.csv.gz"
            with gzip.open(reference, "wt", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "source_id",
                        "contrast_id",
                        "gene_symbol",
                        "ensembl_id",
                        "log2fc_leo_or_iss_minus_ground",
                        "significant_fdr_0_05",
                    ],
                )
                writer.writeheader()
                for source_id, contrast_id in [
                    ("OSD-863", "osd_863_demo"),
                    ("OSD-871", "osd_871_demo"),
                ]:
                    writer.writerows(
                        [
                            {
                                "source_id": source_id,
                                "contrast_id": contrast_id,
                                "gene_symbol": "",
                                "ensembl_id": "ENSG00000000003",
                                "log2fc_leo_or_iss_minus_ground": "1.0",
                                "significant_fdr_0_05": "true",
                            },
                            {
                                "source_id": source_id,
                                "contrast_id": contrast_id,
                                "gene_symbol": "",
                                "ensembl_id": "ENSG00000000005",
                                "log2fc_leo_or_iss_minus_ground": "-1.0",
                                "significant_fdr_0_05": "true",
                            },
                        ]
                    )
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "run_v9_human_organoid_source_transfer_signature.py"),
                    "--repo-root",
                    str(REPO_ROOT),
                    "--task-manifest",
                    str(
                        REPO_ROOT
                        / "v9"
                        / "human_organoid"
                        / "task_manifests"
                        / "draft_human_organoid_spaceflight.json"
                    ),
                    "--de-reference-manifest",
                    str(de_manifest),
                    "--reference-signature-table",
                    str(reference),
                    "--output-dir",
                    str(output_dir),
                    "--max-features",
                    "5",
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            metrics = json.loads((output_dir / "metrics.json").read_text())
            metadata = json.loads((output_dir / "response_signature_metadata.json").read_text())
            run_manifest = json.loads((output_dir / "run_manifest.json").read_text())
            readme = (output_dir / "README.md").read_text()

        self.assertIn("response_signature_metadata.json", result.stdout)
        self.assertEqual(metadata["reference_usage_policy"], "reference_not_used_for_signature_generation")
        self.assertEqual(metadata["n_response_rows"], 10)
        self.assertEqual(metrics["metrics"]["de_direction_match"]["status"], "computed")
        self.assertEqual(metrics["metrics"]["signature_rank_correlation"]["status"], "computed")
        self.assertIn("not a leaderboard result", readme)
        self.assertIn(
            "response_signature",
            {input_file["role"] for input_file in run_manifest["input_files"]},
        )

    def test_microglia_matched_source_transfer_adapter_matches_target_microglia(self):
        from spacebio_bench import (
            MICROGLIA_MATCHED_SOURCE_TRANSFER_SIGNATURE_ID,
            REFERENCE_USAGE_POLICY,
            build_microglia_matched_source_transfer_response_signature,
            load_human_organoid_task,
        )

        task = load_human_organoid_task(repo_root=REPO_ROOT)
        with tempfile.TemporaryDirectory() as tmpdir:
            manifest_path = Path(tmpdir) / "de_reference_manifest.json"
            manifest_path.write_text(
                json.dumps(
                    {
                        "sources": [
                            {
                                "source_id": "OSD-863",
                                "contrasts": [
                                    {
                                        "contrast_id": "osd_863_with",
                                        "disease_context": "no known diseases",
                                        "microglia_condition": "with microglia",
                                    },
                                    {
                                        "contrast_id": "osd_863_without",
                                        "disease_context": "primary progressive multiple sclerosis",
                                        "microglia_condition": "without microglia",
                                    },
                                ],
                            },
                            {
                                "source_id": "OSD-871",
                                "contrasts": [
                                    {
                                        "contrast_id": "osd_871_with",
                                        "disease_context": "no known diseases",
                                        "microglia_condition": "with microglia",
                                    },
                                    {
                                        "contrast_id": "osd_871_without",
                                        "disease_context": "Sporadic Parkinson disease",
                                        "microglia_condition": "without microglia",
                                    },
                                ],
                            },
                        ]
                    }
                )
            )
            result = build_microglia_matched_source_transfer_response_signature(
                task,
                de_reference_manifest_path=manifest_path,
                max_features=3,
            )

        self.assertEqual(
            result.signature_model_id,
            MICROGLIA_MATCHED_SOURCE_TRANSFER_SIGNATURE_ID,
        )
        self.assertEqual(result.n_response_rows, 12)
        self.assertTrue(result.contrast_summaries)
        emitted = [summary for summary in result.contrast_summaries or [] if summary["status"] == "emitted"]
        self.assertEqual(len(emitted), 4)
        for row in result.response_rows:
            self.assertNotEqual(row["source_id"], row["training_source_id"])
            self.assertEqual(row["conditioning_factor"], "microglia_condition")
            self.assertEqual(row["conditioning_value"], row["target_microglia_condition"])
            self.assertGreater(int(row["n_condition_train_ground"]), 0)
            self.assertGreater(int(row["n_condition_train_leo_or_iss"]), 0)
            self.assertEqual(row["reference_usage_policy"], REFERENCE_USAGE_POLICY)

    def test_run_human_organoid_microglia_source_transfer_signature_cli_writes_report(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            output_dir = tmp / "microglia_source_transfer_signature"
            de_manifest = tmp / "de_reference_manifest.json"
            contrasts = [
                ("OSD-863", "osd_863_with", "no known diseases", "with microglia"),
                (
                    "OSD-863",
                    "osd_863_without",
                    "primary progressive multiple sclerosis",
                    "without microglia",
                ),
                ("OSD-871", "osd_871_with", "no known diseases", "with microglia"),
                (
                    "OSD-871",
                    "osd_871_without",
                    "Sporadic Parkinson disease",
                    "without microglia",
                ),
            ]
            de_manifest.write_text(
                json.dumps(
                    {
                        "sources": [
                            {
                                "source_id": "OSD-863",
                                "contrasts": [
                                    {
                                        "contrast_id": contrast_id,
                                        "disease_context": disease,
                                        "microglia_condition": microglia,
                                    }
                                    for source_id, contrast_id, disease, microglia in contrasts
                                    if source_id == "OSD-863"
                                ],
                            },
                            {
                                "source_id": "OSD-871",
                                "contrasts": [
                                    {
                                        "contrast_id": contrast_id,
                                        "disease_context": disease,
                                        "microglia_condition": microglia,
                                    }
                                    for source_id, contrast_id, disease, microglia in contrasts
                                    if source_id == "OSD-871"
                                ],
                            },
                        ]
                    }
                )
            )
            reference = tmp / "reference.csv.gz"
            with gzip.open(reference, "wt", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "source_id",
                        "contrast_id",
                        "gene_symbol",
                        "ensembl_id",
                        "log2fc_leo_or_iss_minus_ground",
                        "significant_fdr_0_05",
                    ],
                )
                writer.writeheader()
                for source_id, contrast_id, _, _ in contrasts:
                    writer.writerows(
                        [
                            {
                                "source_id": source_id,
                                "contrast_id": contrast_id,
                                "gene_symbol": "",
                                "ensembl_id": "ENSG00000000003",
                                "log2fc_leo_or_iss_minus_ground": "1.0",
                                "significant_fdr_0_05": "true",
                            },
                            {
                                "source_id": source_id,
                                "contrast_id": contrast_id,
                                "gene_symbol": "",
                                "ensembl_id": "ENSG00000000005",
                                "log2fc_leo_or_iss_minus_ground": "-1.0",
                                "significant_fdr_0_05": "true",
                            },
                        ]
                    )
            result = subprocess.run(
                [
                    sys.executable,
                    str(
                        REPO_ROOT
                        / "scripts"
                        / "run_v9_human_organoid_microglia_source_transfer_signature.py"
                    ),
                    "--repo-root",
                    str(REPO_ROOT),
                    "--task-manifest",
                    str(
                        REPO_ROOT
                        / "v9"
                        / "human_organoid"
                        / "task_manifests"
                        / "draft_human_organoid_spaceflight.json"
                    ),
                    "--de-reference-manifest",
                    str(de_manifest),
                    "--reference-signature-table",
                    str(reference),
                    "--global-source-transfer-metrics",
                    str(
                        REPO_ROOT
                        / "v9"
                        / "human_organoid"
                        / "reports"
                        / "source_transfer_signature"
                        / "metrics.json"
                    ),
                    "--output-dir",
                    str(output_dir),
                    "--max-features",
                    "5",
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            metrics = json.loads((output_dir / "metrics.json").read_text())
            metadata = json.loads((output_dir / "response_signature_metadata.json").read_text())
            readme = (output_dir / "README.md").read_text()

        self.assertIn("response_signature_metadata.json", result.stdout)
        self.assertEqual(
            metadata["signature_model_id"],
            "organoid_microglia_matched_source_transfer_empirical_signature",
        )
        self.assertEqual(metadata["conditioning_strategy"], "microglia_matched_source_transfer")
        self.assertEqual(metadata["n_response_rows"], 20)
        self.assertEqual(metrics["metrics"]["de_direction_match"]["status"], "computed")
        self.assertEqual(metrics["metrics"]["signature_rank_correlation"]["status"], "computed")
        self.assertIn("Microglia-Matched", readme)
        self.assertIn("plumbing", readme)

    def test_shared_control_source_transfer_adapter_skips_source_specific_diseases(self):
        from spacebio_bench import (
            REFERENCE_USAGE_POLICY,
            SHARED_CONTROL_SOURCE_TRANSFER_SIGNATURE_ID,
            build_shared_control_source_transfer_response_signature,
            load_human_organoid_task,
        )

        task = load_human_organoid_task(repo_root=REPO_ROOT)
        with tempfile.TemporaryDirectory() as tmpdir:
            manifest_path = Path(tmpdir) / "de_reference_manifest.json"
            manifest_path.write_text(
                json.dumps(
                    {
                        "sources": [
                            {
                                "source_id": "OSD-863",
                                "contrasts": [
                                    {
                                        "contrast_id": "osd_863_control",
                                        "disease_context": "no known diseases",
                                        "microglia_condition": "with microglia",
                                    },
                                    {
                                        "contrast_id": "osd_863_ppms",
                                        "disease_context": "primary progressive multiple sclerosis",
                                        "microglia_condition": "without microglia",
                                    },
                                ],
                            },
                            {
                                "source_id": "OSD-871",
                                "contrasts": [
                                    {
                                        "contrast_id": "osd_871_control",
                                        "disease_context": "no known diseases",
                                        "microglia_condition": "without microglia",
                                    },
                                    {
                                        "contrast_id": "osd_871_parkinson",
                                        "disease_context": "Sporadic Parkinson disease",
                                        "microglia_condition": "with microglia",
                                    },
                                ],
                            },
                        ]
                    }
                )
            )
            result = build_shared_control_source_transfer_response_signature(
                task,
                de_reference_manifest_path=manifest_path,
                max_features=3,
            )

        self.assertEqual(result.signature_model_id, SHARED_CONTROL_SOURCE_TRANSFER_SIGNATURE_ID)
        self.assertEqual(result.n_response_rows, 6)
        self.assertTrue(result.contrast_summaries)
        emitted = [summary for summary in result.contrast_summaries or [] if summary["status"] == "emitted"]
        skipped = [
            summary
            for summary in result.contrast_summaries or []
            if summary["status"] == "skipped_missing_shared_disease_context"
        ]
        self.assertEqual(len(emitted), 2)
        self.assertEqual(len(skipped), 2)
        for row in result.response_rows:
            self.assertEqual(row["target_disease_context"], "no_known_diseases")
            self.assertNotEqual(row["source_id"], row["training_source_id"])
            self.assertEqual(
                row["conditioning_factor"],
                "disease_context+microglia_condition",
            )
            self.assertGreater(int(row["n_condition_train_ground"]), 0)
            self.assertGreater(int(row["n_condition_train_leo_or_iss"]), 0)
            self.assertEqual(row["reference_usage_policy"], REFERENCE_USAGE_POLICY)

    def test_run_human_organoid_shared_control_source_transfer_signature_cli_writes_report(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            output_dir = tmp / "shared_control_source_transfer_signature"
            de_manifest = tmp / "de_reference_manifest.json"
            contrasts = [
                ("OSD-863", "osd_863_control", "no known diseases", "with microglia"),
                (
                    "OSD-863",
                    "osd_863_ppms",
                    "primary progressive multiple sclerosis",
                    "without microglia",
                ),
                ("OSD-871", "osd_871_control", "no known diseases", "without microglia"),
                (
                    "OSD-871",
                    "osd_871_parkinson",
                    "Sporadic Parkinson disease",
                    "with microglia",
                ),
            ]
            de_manifest.write_text(
                json.dumps(
                    {
                        "sources": [
                            {
                                "source_id": "OSD-863",
                                "contrasts": [
                                    {
                                        "contrast_id": contrast_id,
                                        "disease_context": disease,
                                        "microglia_condition": microglia,
                                    }
                                    for source_id, contrast_id, disease, microglia in contrasts
                                    if source_id == "OSD-863"
                                ],
                            },
                            {
                                "source_id": "OSD-871",
                                "contrasts": [
                                    {
                                        "contrast_id": contrast_id,
                                        "disease_context": disease,
                                        "microglia_condition": microglia,
                                    }
                                    for source_id, contrast_id, disease, microglia in contrasts
                                    if source_id == "OSD-871"
                                ],
                            },
                        ]
                    }
                )
            )
            reference = tmp / "reference.csv.gz"
            emitted_contrasts = [
                ("OSD-863", "osd_863_control"),
                ("OSD-871", "osd_871_control"),
            ]
            with gzip.open(reference, "wt", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "source_id",
                        "contrast_id",
                        "gene_symbol",
                        "ensembl_id",
                        "log2fc_leo_or_iss_minus_ground",
                        "significant_fdr_0_05",
                    ],
                )
                writer.writeheader()
                for source_id, contrast_id in emitted_contrasts:
                    writer.writerows(
                        [
                            {
                                "source_id": source_id,
                                "contrast_id": contrast_id,
                                "gene_symbol": "",
                                "ensembl_id": "ENSG00000000003",
                                "log2fc_leo_or_iss_minus_ground": "1.0",
                                "significant_fdr_0_05": "true",
                            },
                            {
                                "source_id": source_id,
                                "contrast_id": contrast_id,
                                "gene_symbol": "",
                                "ensembl_id": "ENSG00000000005",
                                "log2fc_leo_or_iss_minus_ground": "-1.0",
                                "significant_fdr_0_05": "true",
                            },
                        ]
                    )
            result = subprocess.run(
                [
                    sys.executable,
                    str(
                        REPO_ROOT
                        / "scripts"
                        / "run_v9_human_organoid_shared_control_source_transfer_signature.py"
                    ),
                    "--repo-root",
                    str(REPO_ROOT),
                    "--task-manifest",
                    str(
                        REPO_ROOT
                        / "v9"
                        / "human_organoid"
                        / "task_manifests"
                        / "draft_human_organoid_spaceflight.json"
                    ),
                    "--de-reference-manifest",
                    str(de_manifest),
                    "--reference-signature-table",
                    str(reference),
                    "--output-dir",
                    str(output_dir),
                    "--max-features",
                    "5",
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            metrics = json.loads((output_dir / "metrics.json").read_text())
            metadata = json.loads((output_dir / "response_signature_metadata.json").read_text())
            readme = (output_dir / "README.md").read_text()

        self.assertIn("response_signature_metadata.json", result.stdout)
        self.assertEqual(
            metadata["signature_model_id"],
            "organoid_shared_control_disease_microglia_source_transfer_empirical_signature",
        )
        self.assertEqual(metadata["n_response_rows"], 10)
        self.assertEqual(metadata["n_emitted_contrasts"], 2)
        self.assertEqual(metadata["n_skipped_contrasts"], 2)
        self.assertTrue(metadata["partial_coverage"])
        self.assertEqual(metrics["metrics"]["de_direction_match"]["status"], "computed")
        self.assertEqual(metrics["metrics"]["signature_rank_correlation"]["status"], "computed")
        self.assertIn("partial-coverage", readme)
        self.assertIn("skipped target contrasts", readme)

    def test_feature_effect_metrics_compute_rank_direction_and_top_k(self):
        from spacebio_bench import compute_feature_effect_metrics

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            reference = tmp / "reference.csv.gz"
            with gzip.open(reference, "wt", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "source_id",
                        "contrast_id",
                        "gene_symbol",
                        "ensembl_id",
                        "log2fc_leo_or_iss_minus_ground",
                        "significant_fdr_0_05",
                    ],
                )
                writer.writeheader()
                writer.writerows(
                    [
                        {
                            "source_id": "OSD-863",
                            "contrast_id": "c1",
                            "gene_symbol": "",
                            "ensembl_id": "ENSG1",
                            "log2fc_leo_or_iss_minus_ground": "2.0",
                            "significant_fdr_0_05": "true",
                        },
                        {
                            "source_id": "OSD-863",
                            "contrast_id": "c1",
                            "gene_symbol": "",
                            "ensembl_id": "ENSG2",
                            "log2fc_leo_or_iss_minus_ground": "-1.0",
                            "significant_fdr_0_05": "true",
                        },
                        {
                            "source_id": "OSD-863",
                            "contrast_id": "c1",
                            "gene_symbol": "",
                            "ensembl_id": "ENSG3",
                            "log2fc_leo_or_iss_minus_ground": "0.5",
                            "significant_fdr_0_05": "false",
                        },
                    ]
                )
            effects = tmp / "feature_effect.csv"
            with effects.open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "task_id",
                        "fold_id",
                        "source_id",
                        "contrast_id",
                        "feature_namespace",
                        "gene_symbol",
                        "ensembl_id",
                        "feature_id",
                        "effect_value",
                        "effect_direction_positive_class",
                        "effect_scale",
                        "model_space",
                        "classifier_model_id",
                        "training_scope",
                        "target_scope",
                        "positive_label",
                        "reference_usage_policy",
                    ],
                )
                writer.writeheader()
                for ensembl_id, effect_value in [
                    ("ENSG1", "3.0"),
                    ("ENSG2", "-2.0"),
                    ("ENSG3", "0.1"),
                ]:
                    writer.writerow(
                        {
                            "task_id": "draft_human_organoid_spaceflight",
                            "fold_id": "fold",
                            "source_id": "OSD-863",
                            "contrast_id": "c1",
                            "feature_namespace": "human_gene",
                            "gene_symbol": "",
                            "ensembl_id": ensembl_id,
                            "feature_id": ensembl_id,
                            "effect_value": effect_value,
                            "effect_direction_positive_class": "LEO_or_ISS",
                            "effect_scale": "standardized_logistic_coefficient",
                            "model_space": "gene_space",
                            "classifier_model_id": "demo",
                            "training_scope": "source_transfer",
                            "target_scope": "target_source_contrast",
                            "positive_label": "LEO_or_ISS",
                            "reference_usage_policy": "reference_not_used_for_effect_generation",
                        }
                    )
            result = compute_feature_effect_metrics(
                feature_effect_path=effects,
                reference_signature_path=reference,
                task_id="draft_human_organoid_spaceflight",
                top_k=(2,),
            )

        self.assertTrue(result["feature_effect_validation"]["ok"])
        self.assertEqual(
            result["metrics"]["feature_effect_direction_match"]["value"],
            1.0,
        )
        self.assertEqual(
            result["metrics"]["feature_effect_top_k_de_overlap"]["value"]["2"]["n_overlap"],
            2,
        )
        top_k_metric = result["metrics"]["feature_effect_top_k_de_overlap"]["value"]["2"]
        self.assertEqual(top_k_metric["n_feature_universe"], 3)
        self.assertAlmostEqual(top_k_metric["expected_overlap"], 4 / 3)
        self.assertAlmostEqual(top_k_metric["enrichment"], 1.5)
        self.assertAlmostEqual(
            top_k_metric["hypergeometric_p_value_greater_equal"],
            1 / 3,
        )
        self.assertEqual(
            top_k_metric["null_model"],
            "hypergeometric_random_top_k_without_replacement",
        )
        self.assertEqual(
            result["metrics"]["feature_effect_rank_correlation"]["status"],
            "computed",
        )

    def test_l2_logistic_feature_effect_adapter_excludes_target_source(self):
        from spacebio_bench import (
            FEATURE_EFFECT_POLICY,
            LOGISTIC_FEATURE_EFFECT_ID,
            build_l2_logistic_feature_effect,
            load_human_organoid_task,
        )

        task = load_human_organoid_task(repo_root=REPO_ROOT)
        with tempfile.TemporaryDirectory() as tmpdir:
            manifest_path = Path(tmpdir) / "de_reference_manifest.json"
            manifest_path.write_text(
                json.dumps(
                    {
                        "sources": [
                            {
                                "source_id": "OSD-863",
                                "contrasts": [{"contrast_id": "osd_863_demo"}],
                            },
                            {
                                "source_id": "OSD-871",
                                "contrasts": [{"contrast_id": "osd_871_demo"}],
                            },
                        ]
                    }
                )
            )
            result = build_l2_logistic_feature_effect(
                task,
                de_reference_manifest_path=manifest_path,
                top_variable_genes=5,
            )

        self.assertEqual(result.classifier_model_id, LOGISTIC_FEATURE_EFFECT_ID)
        self.assertEqual(result.n_effect_rows, 10)
        for row in result.effect_rows:
            self.assertNotEqual(row["training_source_id"], row["target_source_id"])
            self.assertEqual(row["effect_scale"], "standardized_logistic_coefficient")
            self.assertEqual(row["model_space"], "gene_space")
            self.assertEqual(row["reference_usage_policy"], FEATURE_EFFECT_POLICY)

    def test_run_human_organoid_logistic_feature_effect_cli_writes_report(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "logistic_feature_effect"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "run_v9_human_organoid_logistic_feature_effect.py"),
                    "--repo-root",
                    str(REPO_ROOT),
                    "--task-manifest",
                    str(
                        REPO_ROOT
                        / "v9"
                        / "human_organoid"
                        / "task_manifests"
                        / "draft_human_organoid_spaceflight.json"
                    ),
                    "--reference-signature-table",
                    str(
                        REPO_ROOT
                        / "v9"
                        / "human_organoid"
                        / "de_references"
                        / "human_organoid_de_reference.draft.csv.gz"
                    ),
                    "--de-reference-manifest",
                    str(
                        REPO_ROOT
                        / "v9"
                        / "human_organoid"
                        / "de_references"
                        / "human_organoid_de_reference_manifest.draft.json"
                    ),
                    "--output-dir",
                    str(output_dir),
                    "--top-variable-genes",
                    "5",
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            metrics = json.loads((output_dir / "metrics.json").read_text())
            metadata = json.loads((output_dir / "feature_effect_metadata.json").read_text())
            run_manifest = json.loads((output_dir / "run_manifest.json").read_text())
            readme = (output_dir / "README.md").read_text()

        self.assertIn("feature_effect_metadata.json", result.stdout)
        self.assertEqual(metadata["classifier_model_id"], "organoid_l2_logistic_gene_space_feature_effect")
        self.assertEqual(metadata["n_effect_rows"], 40)
        self.assertEqual(metrics["metrics"]["feature_effect_direction_match"]["status"], "computed")
        self.assertEqual(metrics["metrics"]["feature_effect_rank_correlation"]["status"], "computed")
        top_k_value = metrics["metrics"]["feature_effect_top_k_de_overlap"]["value"]["50"]
        self.assertIn("expected_overlap", top_k_value)
        self.assertIn("enrichment", top_k_value)
        self.assertIn("hypergeometric_p_value_greater_equal", top_k_value)
        self.assertIn("not log2FC", readme)
        self.assertIn("Hypergeometric p>=x", readme)
        self.assertIn("feature_effect", {item["role"] for item in run_manifest["input_files"]})

    def test_pca_lr_reconstructed_feature_effect_matches_reconstruction_formula(self):
        import numpy as np
        import pandas as pd
        from sklearn.decomposition import PCA
        from sklearn.linear_model import LogisticRegression

        from spacebio_bench import (
            PCA_LR_FEATURE_EFFECT_ID,
            build_pca_lr_reconstructed_feature_effect,
        )
        from spacebio_bench.data import HumanOrganoidTaskData, OrganoidFoldData

        train_ids = ["b_ground_1", "b_flight_1", "b_ground_2", "b_flight_2"]
        test_ids = ["a_ground_1", "a_flight_1"]
        features = pd.DataFrame(
            [
                [1.0, 2.0, 0.5, 4.0],
                [3.0, 1.0, 2.5, 5.5],
                [1.5, 2.2, 0.7, 3.8],
                [3.2, 0.8, 2.8, 5.8],
                [0.2, 1.0, 0.1, 2.0],
                [2.0, 0.3, 1.7, 3.3],
            ],
            index=train_ids + test_ids,
            columns=["ENSG1", "ENSG2", "ENSG3", "ENSG4"],
        )
        sample_factors = [
            {
                "parse_status": "parsed",
                "sample_id": sample_id,
                "source_id": "OSD-B",
                "true_label": "LEO_or_ISS" if "flight" in sample_id else "Ground",
            }
            for sample_id in train_ids
        ] + [
            {
                "parse_status": "parsed",
                "sample_id": sample_id,
                "source_id": "OSD-A",
                "true_label": "LEO_or_ISS" if "flight" in sample_id else "Ground",
            }
            for sample_id in test_ids
        ]
        task = HumanOrganoidTaskData(
            manifest={"task_id": "draft_human_organoid_spaceflight"},
            features=features,
            sample_factors=sample_factors,
            folds=[
                OrganoidFoldData(
                    task_id="draft_human_organoid_spaceflight",
                    fold_id="fold_source_transfer",
                    heldout_factor="source_id",
                    heldout_value="OSD-A",
                    train_sample_ids=train_ids,
                    test_sample_ids=test_ids,
                    train_label_distribution={"Ground": 2, "LEO_or_ISS": 2},
                    test_label_distribution={"Ground": 1, "LEO_or_ISS": 1},
                    status="active",
                )
            ],
            diagnostic_folds=[],
            matrix_paths={},
            feature_namespace="human_gene",
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            manifest_path = Path(tmpdir) / "de_reference_manifest.json"
            manifest_path.write_text(
                json.dumps(
                    {
                        "sources": [
                            {
                                "source_id": "OSD-A",
                                "contrasts": [{"contrast_id": "contrast_a"}],
                            }
                        ]
                    }
                )
            )
            result = build_pca_lr_reconstructed_feature_effect(
                task,
                de_reference_manifest_path=manifest_path,
                transform="none",
                top_variable_genes=4,
                pca_components=2,
            )

        train_values = features.loc[train_ids].to_numpy(dtype=np.float64)
        scaled = (train_values - train_values.mean(axis=0)) / train_values.std(axis=0)
        pca = PCA(n_components=2, random_state=42, whiten=False)
        pc_scores = pca.fit_transform(scaled)
        model = LogisticRegression(
            C=1.0,
            class_weight="balanced",
            max_iter=5000,
            random_state=42,
            solver="lbfgs",
        )
        model.fit(pc_scores, np.asarray([0, 1, 0, 1], dtype=int))
        expected = pca.components_.T @ model.coef_[0]

        self.assertEqual(result.classifier_model_id, PCA_LR_FEATURE_EFFECT_ID)
        self.assertEqual(result.n_effect_rows, 4)
        observed = {
            row["ensembl_id"]: float(row["effect_value"])
            for row in result.effect_rows
        }
        for feature_id, expected_value in zip(features.columns, expected):
            self.assertAlmostEqual(observed[feature_id], float(expected_value))
        for row in result.effect_rows:
            self.assertEqual(
                row["effect_scale"],
                "pca_lr_reconstructed_standardized_logistic_coefficient",
            )
            self.assertEqual(row["model_space"], "reconstructed_gene_space_from_pca")
            self.assertNotEqual(row["training_source_id"], row["target_source_id"])

    def test_run_human_organoid_pca_lr_feature_effect_cli_writes_report(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "pca_lr_feature_effect"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "run_v9_human_organoid_pca_lr_feature_effect.py"),
                    "--repo-root",
                    str(REPO_ROOT),
                    "--task-manifest",
                    str(
                        REPO_ROOT
                        / "v9"
                        / "human_organoid"
                        / "task_manifests"
                        / "draft_human_organoid_spaceflight.json"
                    ),
                    "--reference-signature-table",
                    str(
                        REPO_ROOT
                        / "v9"
                        / "human_organoid"
                        / "de_references"
                        / "human_organoid_de_reference.draft.csv.gz"
                    ),
                    "--de-reference-manifest",
                    str(
                        REPO_ROOT
                        / "v9"
                        / "human_organoid"
                        / "de_references"
                        / "human_organoid_de_reference_manifest.draft.json"
                    ),
                    "--output-dir",
                    str(output_dir),
                    "--top-variable-genes",
                    "5",
                    "--pca-components",
                    "2",
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            metrics = json.loads((output_dir / "metrics.json").read_text())
            metadata = json.loads((output_dir / "feature_effect_metadata.json").read_text())
            run_manifest = json.loads((output_dir / "run_manifest.json").read_text())
            readme = (output_dir / "README.md").read_text()

        self.assertIn("feature_effect_metadata.json", result.stdout)
        self.assertEqual(
            metadata["classifier_model_id"],
            "organoid_pca_lr_reconstructed_gene_space_feature_effect",
        )
        self.assertEqual(metadata["n_effect_rows"], 40)
        self.assertEqual(metrics["metrics"]["feature_effect_direction_match"]["status"], "computed")
        self.assertEqual(metrics["metrics"]["feature_effect_rank_correlation"]["status"], "computed")
        top_k_value = metrics["metrics"]["feature_effect_top_k_de_overlap"]["value"]["50"]
        self.assertIn("expected_overlap", top_k_value)
        self.assertIn("hypergeometric_p_value_greater_equal", top_k_value)
        self.assertIn("pca.components_.T @ logistic.coef_[0]", readme)
        self.assertIn("feature_effect", {item["role"] for item in run_manifest["input_files"]})

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

    def test_multispecies_sample_factor_builder_parses_local_condition_strata(self):
        from spacebio_bench import (
            build_multispecies_sample_factor_rows,
            read_extension_source_inventory,
        )

        source_rows = read_extension_source_inventory(
            REPO_ROOT / "v9" / "multispecies" / "source_inventory.draft.csv"
        )
        rows = build_multispecies_sample_factor_rows(source_rows, repo_root=REPO_ROOT)
        rows_by_source = Counter(row["source_id"] for row in rows)
        labels_by_source = Counter((row["source_id"], row["true_label"]) for row in rows)

        self.assertEqual(len(rows), 124)
        self.assertEqual(rows_by_source, {"OSD-207": 32, "OSD-37": 56, "OSD-120": 36})
        self.assertEqual({row["parse_status"] for row in rows}, {"parsed"})
        self.assertEqual(labels_by_source[("OSD-207", "Ground")], 16)
        self.assertEqual(labels_by_source[("OSD-207", "LEO_or_ISS")], 16)
        self.assertEqual(labels_by_source[("OSD-37", "Ground")], 28)
        self.assertEqual(labels_by_source[("OSD-37", "LEO_or_ISS")], 28)
        self.assertEqual(labels_by_source[("OSD-120", "Ground")], 18)
        self.assertEqual(labels_by_source[("OSD-120", "LEO_or_ISS")], 18)
        self.assertEqual(
            {
                row["genotype_or_ecotype"]
                for row in rows
                if row["source_id"] == "OSD-37"
            },
            {"Col.0", "Cvi.0", "Ler.0", "Ws.2"},
        )
        self.assertEqual(
            {
                row["light_treatment"]
                for row in rows
                if row["source_id"] == "OSD-120"
            },
            {"Dark.Treatment", "Light.Treatment"},
        )

    def test_multispecies_expression_matrix_audit_aligns_local_columns(self):
        from spacebio_bench import (
            audit_multispecies_expression_matrices,
            build_multispecies_sample_factor_rows,
            read_extension_source_inventory,
        )

        source_rows = read_extension_source_inventory(
            REPO_ROOT / "v9" / "multispecies" / "source_inventory.draft.csv"
        )
        sample_factor_rows = build_multispecies_sample_factor_rows(
            source_rows,
            repo_root=REPO_ROOT,
        )
        rows = audit_multispecies_expression_matrices(
            source_rows,
            sample_factor_rows=sample_factor_rows,
            repo_root=REPO_ROOT,
        )
        by_source = {row["source_id"]: row for row in rows}

        self.assertEqual(set(by_source), {"OSD-207", "OSD-37", "OSD-120"})
        self.assertEqual({row["audit_status"] for row in rows}, {"matrix_local_sample_aligned"})
        self.assertEqual({row["matrix_columns_match_sample_factors"] for row in rows}, {"True"})
        self.assertEqual(by_source["OSD-207"]["n_feature_rows"], "15999")
        self.assertEqual(by_source["OSD-207"]["n_sample_columns"], "32")
        self.assertEqual(by_source["OSD-37"]["n_feature_rows"], "28067")
        self.assertEqual(by_source["OSD-37"]["n_sample_columns"], "56")
        self.assertEqual(by_source["OSD-120"]["n_feature_rows"], "24740")
        self.assertEqual(by_source["OSD-120"]["n_sample_columns"], "36")
        self.assertTrue(all(row["missing_sample_factor_rows"] == "" for row in rows))
        self.assertTrue(all(row["extra_matrix_columns"] == "" for row in rows))

    def test_audit_v9_multispecies_inputs_cli_writes_temp_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "audit_v9_multispecies_inputs.py"),
                    "--source-inventory",
                    str(REPO_ROOT / "v9" / "multispecies" / "source_inventory.draft.csv"),
                    "--repo-root",
                    str(REPO_ROOT),
                    "--sample-factors-csv",
                    str(tmp / "sample_factors.csv"),
                    "--sample-factors-json",
                    str(tmp / "sample_factors.json"),
                    "--matrix-audit-csv",
                    str(tmp / "matrix_audit.csv"),
                    "--matrix-audit-json",
                    str(tmp / "matrix_audit.json"),
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )

            with (tmp / "sample_factors.csv").open(newline="") as handle:
                sample_rows = list(csv.DictReader(handle))
            with (tmp / "matrix_audit.csv").open(newline="") as handle:
                matrix_rows = list(csv.DictReader(handle))

        self.assertIn("sample_factors.csv", result.stdout)
        self.assertIn("matrix_audit.csv", result.stdout)
        self.assertEqual(len(sample_rows), 124)
        self.assertEqual(len(matrix_rows), 3)
        self.assertEqual({row["matrix_columns_match_sample_factors"] for row in matrix_rows}, {"True"})

    def test_multispecies_task_manifest_builder_creates_species_native_manifests(self):
        from spacebio_bench import (
            build_multispecies_task_manifests,
            read_extension_source_inventory,
        )

        source_rows = read_extension_source_inventory(
            REPO_ROOT / "v9" / "multispecies" / "source_inventory.draft.csv"
        )
        with (REPO_ROOT / "v9" / "multispecies" / "sample_factors.draft.csv").open(
            newline=""
        ) as handle:
            sample_factor_rows = list(csv.DictReader(handle))
        with (
            REPO_ROOT / "v9" / "multispecies" / "expression_matrix_audit.draft.csv"
        ).open(newline="") as handle:
            expression_matrix_audit_rows = list(csv.DictReader(handle))
        with (
            REPO_ROOT / "v9" / "multispecies" / "source_checksum_audit.draft.csv"
        ).open(newline="") as handle:
            source_checksum_audit_rows = list(csv.DictReader(handle))

        manifests = build_multispecies_task_manifests(
            source_rows,
            sample_factor_rows=sample_factor_rows,
            expression_matrix_audit_rows=expression_matrix_audit_rows,
            source_checksum_audit_rows=source_checksum_audit_rows,
        )
        by_task = {manifest["task_id"]: manifest for manifest in manifests}

        self.assertEqual(
            set(by_task),
            {
                "draft_osd37_arabidopsis_seedling_spaceflight",
                "draft_osd207_drosophila_whole_body_spaceflight",
            },
        )
        arabidopsis = by_task["draft_osd37_arabidopsis_seedling_spaceflight"]
        drosophila = by_task["draft_osd207_drosophila_whole_body_spaceflight"]
        self.assertEqual(arabidopsis["organism"], "Arabidopsis thaliana")
        self.assertEqual(arabidopsis["split"]["n_samples"], 56)
        self.assertEqual(arabidopsis["split"]["label_distribution"], {"Ground": 28, "LEO_or_ISS": 28})
        self.assertEqual(len(arabidopsis["split"]["candidate_folds"]), 4)
        self.assertEqual(
            arabidopsis["split"]["expression_matrix_sources"]["OSD-37"]["n_feature_rows"],
            28067,
        )
        self.assertEqual(
            arabidopsis["provenance"]["source_checksum_status"],
            "checksum_manifest_parsed",
        )
        self.assertEqual(drosophila["organism"], "Drosophila melanogaster")
        self.assertEqual(drosophila["split"]["n_samples"], 32)
        self.assertEqual(drosophila["split"]["label_distribution"], {"Ground": 16, "LEO_or_ISS": 16})
        self.assertEqual(len(drosophila["split"]["candidate_folds"]), 4)
        self.assertEqual(
            {source["source_id"] for source in arabidopsis["source_records"] + drosophila["source_records"]},
            {"OSD-37", "OSD-207"},
        )

    def test_build_v9_multispecies_task_manifests_cli_writes_temp_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "build_v9_multispecies_task_manifests.py"),
                    "--source-inventory",
                    str(REPO_ROOT / "v9" / "multispecies" / "source_inventory.draft.csv"),
                    "--sample-factor-table",
                    str(REPO_ROOT / "v9" / "multispecies" / "sample_factors.draft.csv"),
                    "--expression-matrix-audit",
                    str(REPO_ROOT / "v9" / "multispecies" / "expression_matrix_audit.draft.csv"),
                    "--source-checksum-audit",
                    str(REPO_ROOT / "v9" / "multispecies" / "source_checksum_audit.draft.csv"),
                    "--output-dir",
                    str(tmp / "task_manifests"),
                    "--index-csv",
                    str(tmp / "task_manifest_index.csv"),
                    "--index-json",
                    str(tmp / "task_manifest_index.json"),
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )

            manifests = sorted((tmp / "task_manifests").glob("*.json"))
            with (tmp / "task_manifest_index.csv").open(newline="") as handle:
                index_rows = list(csv.DictReader(handle))

        self.assertIn("draft_osd37_arabidopsis_seedling_spaceflight.json", result.stdout)
        self.assertEqual(len(manifests), 2)
        self.assertEqual(len(index_rows), 2)
        self.assertEqual(
            {row["task_family"] for row in index_rows},
            {"multispecies_species_native_spaceflight"},
        )

    def test_osd120_interaction_task_manifest_builder_records_all_fold_families(self):
        from spacebio_bench import (
            build_osd120_interaction_task_manifest,
            read_extension_source_inventory,
        )

        source_rows = read_extension_source_inventory(
            REPO_ROOT / "v9" / "multispecies" / "source_inventory.draft.csv"
        )
        with (REPO_ROOT / "v9" / "multispecies" / "sample_factors.draft.csv").open(
            newline=""
        ) as handle:
            sample_factor_rows = list(csv.DictReader(handle))
        with (
            REPO_ROOT / "v9" / "multispecies" / "expression_matrix_audit.draft.csv"
        ).open(newline="") as handle:
            expression_matrix_audit_rows = list(csv.DictReader(handle))
        with (
            REPO_ROOT / "v9" / "multispecies" / "source_checksum_audit.draft.csv"
        ).open(newline="") as handle:
            source_checksum_audit_rows = list(csv.DictReader(handle))

        manifest = build_osd120_interaction_task_manifest(
            source_rows,
            sample_factor_rows=sample_factor_rows,
            expression_matrix_audit_rows=expression_matrix_audit_rows,
            source_checksum_audit_rows=source_checksum_audit_rows,
        )

        self.assertEqual(
            manifest["task_id"],
            "draft_osd120_arabidopsis_root_light_interaction_spaceflight",
        )
        self.assertEqual(manifest["task_family"], "multispecies_light_interaction_spaceflight")
        self.assertEqual(manifest["split"]["n_samples"], 36)
        self.assertEqual(manifest["split"]["label_distribution"], {"Ground": 18, "LEO_or_ISS": 18})
        self.assertEqual(
            manifest["split"]["genotype_or_ecotype_distribution"],
            {"Col.0": 12, "Col.0.PhyD": 12, "Wassilewskija.ecotype": 12},
        )
        self.assertEqual(
            manifest["split"]["light_treatment_distribution"],
            {"Dark.Treatment": 18, "Light.Treatment": 18},
        )
        self.assertEqual(len(manifest["split"]["primary_candidate_folds"]), 3)
        self.assertEqual(len(manifest["split"]["secondary_light_treatment_folds"]), 2)
        self.assertEqual(len(manifest["split"]["condition_stratum_diagnostic_folds"]), 6)
        self.assertEqual(
            {fold["n_test"] for fold in manifest["split"]["primary_candidate_folds"]},
            {12},
        )
        self.assertEqual(
            {fold["n_test"] for fold in manifest["split"]["condition_stratum_diagnostic_folds"]},
            {6},
        )
        self.assertEqual(
            manifest["split"]["expression_matrix_sources"]["OSD-120"]["n_feature_rows"],
            24740,
        )
        self.assertEqual(
            manifest["provenance"]["source_checksum_status"],
            "checksum_manifest_parsed",
        )

    def test_build_v9_osd120_interaction_task_manifest_cli_writes_temp_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "build_v9_osd120_interaction_task_manifest.py"),
                    "--source-inventory",
                    str(REPO_ROOT / "v9" / "multispecies" / "source_inventory.draft.csv"),
                    "--sample-factor-table",
                    str(REPO_ROOT / "v9" / "multispecies" / "sample_factors.draft.csv"),
                    "--expression-matrix-audit",
                    str(REPO_ROOT / "v9" / "multispecies" / "expression_matrix_audit.draft.csv"),
                    "--source-checksum-audit",
                    str(REPO_ROOT / "v9" / "multispecies" / "source_checksum_audit.draft.csv"),
                    "--output-dir",
                    str(tmp / "interaction_task_manifests"),
                    "--index-csv",
                    str(tmp / "interaction_task_manifest_index.csv"),
                    "--index-json",
                    str(tmp / "interaction_task_manifest_index.json"),
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )

            manifests = sorted((tmp / "interaction_task_manifests").glob("*.json"))
            with (tmp / "interaction_task_manifest_index.csv").open(newline="") as handle:
                index_rows = list(csv.DictReader(handle))

        self.assertIn("draft_osd120_arabidopsis_root_light_interaction_spaceflight.json", result.stdout)
        self.assertEqual(len(manifests), 1)
        self.assertEqual(len(index_rows), 1)
        self.assertEqual(index_rows[0]["task_family"], "multispecies_light_interaction_spaceflight")

    def test_generated_v9_osd120_interaction_task_manifest_validates(self):
        from spacebio_bench import TaskRegistry, load_task_manifest

        manifest_dir = REPO_ROOT / "v9" / "multispecies" / "interaction_task_manifests"
        manifest_path = (
            manifest_dir / "draft_osd120_arabidopsis_root_light_interaction_spaceflight.json"
        )
        manifest = load_task_manifest(manifest_path)

        self.assertEqual(
            manifest["task_id"],
            "draft_osd120_arabidopsis_root_light_interaction_spaceflight",
        )
        self.assertEqual(manifest["task_family"], "multispecies_light_interaction_spaceflight")
        self.assertEqual(manifest["split"]["n_samples"], 36)
        self.assertEqual(len(manifest["split"]["primary_candidate_folds"]), 3)
        self.assertEqual(len(manifest["split"]["secondary_light_treatment_folds"]), 2)
        self.assertEqual(len(manifest["split"]["condition_stratum_diagnostic_folds"]), 6)
        self.assertEqual(
            manifest["provenance"]["source_checksum_status"],
            "checksum_manifest_parsed",
        )
        registry = TaskRegistry.from_dir(manifest_dir)
        self.assertEqual(len(registry), 1)
        self.assertEqual(registry.task_ids(), [manifest["task_id"]])

    def test_multispecies_interaction_loader_aligns_osd120_matrix_and_folds(self):
        from spacebio_bench import load_multispecies_interaction_task

        task = load_multispecies_interaction_task(
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "interaction_task_manifests"
            / "draft_osd120_arabidopsis_root_light_interaction_spaceflight.json",
            repo_root=REPO_ROOT,
        )

        self.assertEqual(task.n_samples, 36)
        self.assertEqual(task.n_features, 24740)
        self.assertEqual(set(task.matrix_paths), {"OSD-120"})
        self.assertEqual(len(task.folds), 3)
        self.assertEqual(
            {key: len(value) for key, value in task.fold_families.items()},
            {
                "primary_genotype_or_ecotype_holdout": 3,
                "secondary_light_treatment_holdout": 2,
                "diagnostic_condition_stratum_holdout": 6,
            },
        )
        self.assertEqual(
            {fold.n_test for fold in task.fold_families["primary_genotype_or_ecotype_holdout"]},
            {12},
        )
        self.assertEqual(
            {fold.n_test for fold in task.fold_families["secondary_light_treatment_holdout"]},
            {18},
        )
        self.assertEqual(
            {fold.n_test for fold in task.fold_families["diagnostic_condition_stratum_holdout"]},
            {6},
        )

    def test_multispecies_interaction_baseline_writes_fold_family_reports(self):
        from spacebio_bench import (
            MultispeciesBaselineConfig,
            run_multispecies_interaction_baselines,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            _, summary = run_multispecies_interaction_baselines(
                repo_root=REPO_ROOT,
                manifest_dir=REPO_ROOT / "v9" / "multispecies" / "interaction_task_manifests",
                output_root=Path(tmpdir) / "interaction_nearest_centroid",
                config=MultispeciesBaselineConfig(top_variable_genes=100),
                command=["test-osd120-interaction-baseline"],
            )
            with summary["csv"].open(newline="") as handle:
                rows = list(csv.DictReader(handle))

        self.assertEqual(len(rows), 3)
        self.assertEqual(
            {row["fold_family"] for row in rows},
            {
                "primary_genotype_or_ecotype_holdout",
                "secondary_light_treatment_holdout",
                "diagnostic_condition_stratum_holdout",
            },
        )
        self.assertEqual(
            {row["baseline_id"] for row in rows},
            {"multispecies_interaction_nearest_centroid"},
        )
        self.assertEqual({row["status"] for row in rows}, {"evaluated"})
        self.assertEqual({row["n_predictions"] for row in rows}, {"36"})
        primary = next(
            row for row in rows if row["fold_family"] == "primary_genotype_or_ecotype_holdout"
        )
        secondary = next(
            row for row in rows if row["fold_family"] == "secondary_light_treatment_holdout"
        )
        diagnostic = next(
            row for row in rows if row["fold_family"] == "diagnostic_condition_stratum_holdout"
        )
        self.assertTrue(primary["genotype_holdout_delta"])
        self.assertTrue(secondary["light_treatment_holdout_delta"])
        self.assertTrue(diagnostic["condition_stratum_holdout_delta"])

    def test_run_v9_osd120_interaction_baseline_cli_writes_summary(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_nearest_centroid"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "run_v9_osd120_interaction_baseline.py"),
                    "--repo-root",
                    str(REPO_ROOT),
                    "--manifest-dir",
                    str(REPO_ROOT / "v9" / "multispecies" / "interaction_task_manifests"),
                    "--output-dir",
                    str(output_dir),
                    "--top-variable-genes",
                    "100",
                    "--fold-family",
                    "primary_genotype_or_ecotype_holdout",
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            with (output_dir / "multispecies_baseline_summary.csv").open(
                newline=""
            ) as handle:
                rows = list(csv.DictReader(handle))

        self.assertIn("multispecies_baseline_summary.csv", result.stdout)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["fold_family"], "primary_genotype_or_ecotype_holdout")
        self.assertEqual(rows[0]["status"], "evaluated")

    def test_generated_v9_osd120_interaction_baseline_outputs_validate(self):
        summary_path = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "interaction_nearest_centroid"
            / "multispecies_baseline_summary.csv"
        )
        with summary_path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))

        self.assertEqual(len(rows), 3)
        self.assertEqual(
            {row["fold_family"] for row in rows},
            {
                "primary_genotype_or_ecotype_holdout",
                "secondary_light_treatment_holdout",
                "diagnostic_condition_stratum_holdout",
            },
        )
        self.assertEqual(
            {row["claim_boundary"] for row in rows},
            {"draft_interaction_feasibility_only_not_leaderboard"},
        )
        for row in rows:
            self.assertTrue((REPO_ROOT / row["predictions"]).exists())
            self.assertTrue((REPO_ROOT / row["metrics"]).exists())
            self.assertTrue((REPO_ROOT / row["run_manifest"]).exists())

    def test_multispecies_interaction_sensitivity_grid_writes_variant_summary(self):
        from spacebio_bench import (
            MultispeciesBaselineConfig,
            run_multispecies_interaction_sensitivity_grid,
        )

        configs = [
            MultispeciesBaselineConfig(transform="log1p", scaling="zscore", top_variable_genes=100),
            MultispeciesBaselineConfig(transform="none", scaling="none", top_variable_genes=100),
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            _, summary = run_multispecies_interaction_sensitivity_grid(
                repo_root=REPO_ROOT,
                manifest_dir=REPO_ROOT / "v9" / "multispecies" / "interaction_task_manifests",
                output_root=Path(tmpdir) / "interaction_sensitivity",
                configs=configs,
                command=["test-osd120-interaction-sensitivity"],
            )
            with summary["csv"].open(newline="") as handle:
                rows = list(csv.DictReader(handle))

        self.assertEqual(len(rows), 6)
        self.assertEqual(
            {row["variant_id"] for row in rows},
            {"tvg100_log1p_zscore", "tvg100_none_none"},
        )
        self.assertEqual(
            {row["fold_family"] for row in rows},
            {
                "primary_genotype_or_ecotype_holdout",
                "secondary_light_treatment_holdout",
                "diagnostic_condition_stratum_holdout",
            },
        )
        self.assertTrue(all(row["status"] == "evaluated" for row in rows))
        self.assertEqual(
            {row["claim_boundary"] for row in rows},
            {"draft_interaction_sensitivity_only_not_leaderboard"},
        )

    def test_run_v9_osd120_interaction_sensitivity_cli_writes_summary(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_sensitivity"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "run_v9_osd120_interaction_sensitivity.py"),
                    "--repo-root",
                    str(REPO_ROOT),
                    "--manifest-dir",
                    str(REPO_ROOT / "v9" / "multispecies" / "interaction_task_manifests"),
                    "--output-dir",
                    str(output_dir),
                    "--top-variable-genes",
                    "100",
                    "--transform",
                    "log1p",
                    "--scaling",
                    "zscore",
                    "--fold-family",
                    "secondary_light_treatment_holdout",
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            with (output_dir / "multispecies_baseline_summary.csv").open(
                newline=""
            ) as handle:
                rows = list(csv.DictReader(handle))

        self.assertIn("multispecies_baseline_summary.csv", result.stdout)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["variant_id"], "tvg100_log1p_zscore")
        self.assertEqual(rows[0]["fold_family"], "secondary_light_treatment_holdout")
        self.assertEqual(rows[0]["status"], "evaluated")

    def test_generated_v9_osd120_interaction_sensitivity_outputs_validate(self):
        summary_path = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "interaction_sensitivity"
            / "multispecies_baseline_summary.csv"
        )
        with summary_path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))

        self.assertEqual(len(rows), 60)
        self.assertEqual(len({row["variant_id"] for row in rows}), 20)
        self.assertEqual(
            Counter(row["fold_family"] for row in rows),
            Counter(
                {
                    "primary_genotype_or_ecotype_holdout": 20,
                    "secondary_light_treatment_holdout": 20,
                    "diagnostic_condition_stratum_holdout": 20,
                }
            ),
        )
        self.assertEqual({row["status"] for row in rows}, {"evaluated"})
        self.assertEqual(
            {row["claim_boundary"] for row in rows},
            {"draft_interaction_sensitivity_only_not_leaderboard"},
        )
        for row in rows:
            self.assertTrue((REPO_ROOT / row["predictions"]).exists())
            self.assertTrue((REPO_ROOT / row["metrics"]).exists())
            self.assertTrue((REPO_ROOT / row["run_manifest"]).exists())

    def test_multispecies_interaction_fold_detail_aggregation_writes_rows(self):
        from spacebio_bench import (
            MultispeciesBaselineConfig,
            run_multispecies_interaction_sensitivity_grid,
            write_multispecies_interaction_fold_details,
        )

        configs = [
            MultispeciesBaselineConfig(transform="log1p", scaling="zscore", top_variable_genes=100),
            MultispeciesBaselineConfig(transform="none", scaling="none", top_variable_genes=100),
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            _, summary = run_multispecies_interaction_sensitivity_grid(
                repo_root=REPO_ROOT,
                manifest_dir=REPO_ROOT / "v9" / "multispecies" / "interaction_task_manifests",
                output_root=tmp / "interaction_sensitivity",
                configs=configs,
                command=["test-osd120-fold-detail-aggregation"],
            )
            outputs = write_multispecies_interaction_fold_details(
                summary_csv=summary["csv"],
                repo_root=REPO_ROOT,
                output_dir=tmp / "interaction_sensitivity",
            )
            with outputs["csv"].open(newline="") as handle:
                rows = list(csv.DictReader(handle))

        self.assertEqual(len(rows), 22)
        self.assertEqual(
            Counter(row["fold_family"] for row in rows),
            Counter(
                {
                    "primary_genotype_or_ecotype_holdout": 6,
                    "secondary_light_treatment_holdout": 4,
                    "diagnostic_condition_stratum_holdout": 12,
                }
            ),
        )
        self.assertIn("is_lowest_for_variant", rows[0])
        self.assertIn("rank_low_to_high", rows[0])

    def test_build_v9_osd120_interaction_fold_details_cli_writes_summary(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "fold_details"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "build_v9_osd120_interaction_fold_details.py"),
                    "--repo-root",
                    str(REPO_ROOT),
                    "--summary-csv",
                    str(
                        REPO_ROOT
                        / "v9"
                        / "multispecies"
                        / "reports"
                        / "interaction_sensitivity"
                        / "multispecies_baseline_summary.csv"
                    ),
                    "--output-dir",
                    str(output_dir),
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            with (output_dir / "fold_detail_summary.csv").open(newline="") as handle:
                rows = list(csv.DictReader(handle))

        self.assertIn("fold_detail_summary.csv", result.stdout)
        self.assertEqual(len(rows), 220)
        self.assertEqual(
            {row["delta_metric_id"] for row in rows},
            {
                "genotype_holdout_delta",
                "light_treatment_holdout_delta",
                "condition_stratum_holdout_delta",
            },
        )

    def test_generated_v9_osd120_interaction_fold_detail_outputs_validate(self):
        path = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "interaction_sensitivity"
            / "fold_detail_summary.csv"
        )
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))

        self.assertEqual(len(rows), 220)
        self.assertEqual(
            Counter(row["fold_family"] for row in rows),
            Counter(
                {
                    "primary_genotype_or_ecotype_holdout": 60,
                    "secondary_light_treatment_holdout": 40,
                    "diagnostic_condition_stratum_holdout": 120,
                }
            ),
        )
        lowest = Counter(
            (row["fold_family"], row["heldout_value"])
            for row in rows
            if row["is_lowest_for_variant"] == "True"
        )
        self.assertEqual(
            lowest[("secondary_light_treatment_holdout", "Light.Treatment")],
            19,
        )
        self.assertEqual(
            lowest[
                (
                    "diagnostic_condition_stratum_holdout",
                    "Wassilewskija.ecotype|Dark.Treatment",
                )
            ],
            16,
        )
        self.assertEqual(
            lowest[("primary_genotype_or_ecotype_holdout", "Wassilewskija.ecotype")],
            12,
        )

    def test_multispecies_interaction_logistic_baseline_writes_fold_family_reports(self):
        from spacebio_bench import (
            MultispeciesLogisticBaselineConfig,
            run_multispecies_interaction_logistic_baselines,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            _, summary = run_multispecies_interaction_logistic_baselines(
                repo_root=REPO_ROOT,
                manifest_dir=REPO_ROOT / "v9" / "multispecies" / "interaction_task_manifests",
                output_root=Path(tmpdir) / "interaction_logistic_l2",
                config=MultispeciesLogisticBaselineConfig(top_variable_genes=100),
                command=["test-osd120-interaction-logistic"],
            )
            with summary["csv"].open(newline="") as handle:
                rows = list(csv.DictReader(handle))

        self.assertEqual(len(rows), 3)
        self.assertEqual(
            {row["baseline_id"] for row in rows},
            {"multispecies_interaction_logistic_regression_l2"},
        )
        self.assertEqual(
            {row["fold_family"] for row in rows},
            {
                "primary_genotype_or_ecotype_holdout",
                "secondary_light_treatment_holdout",
                "diagnostic_condition_stratum_holdout",
            },
        )
        self.assertEqual({row["status"] for row in rows}, {"evaluated"})
        self.assertEqual({row["n_predictions"] for row in rows}, {"36"})
        self.assertEqual(
            {row["claim_boundary"] for row in rows},
            {"draft_interaction_logistic_only_not_leaderboard"},
        )

    def test_run_v9_osd120_interaction_logistic_cli_writes_summary(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_logistic_l2"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "run_v9_osd120_interaction_logistic.py"),
                    "--repo-root",
                    str(REPO_ROOT),
                    "--manifest-dir",
                    str(REPO_ROOT / "v9" / "multispecies" / "interaction_task_manifests"),
                    "--output-dir",
                    str(output_dir),
                    "--top-variable-genes",
                    "100",
                    "--fold-family",
                    "diagnostic_condition_stratum_holdout",
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            with (output_dir / "multispecies_baseline_summary.csv").open(
                newline=""
            ) as handle:
                rows = list(csv.DictReader(handle))

        self.assertIn("multispecies_baseline_summary.csv", result.stdout)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["fold_family"], "diagnostic_condition_stratum_holdout")
        self.assertEqual(rows[0]["baseline_id"], "multispecies_interaction_logistic_regression_l2")
        self.assertEqual(rows[0]["status"], "evaluated")

    def test_generated_v9_osd120_interaction_logistic_outputs_validate(self):
        summary_path = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "interaction_logistic_l2"
            / "multispecies_baseline_summary.csv"
        )
        with summary_path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))

        self.assertEqual(len(rows), 3)
        self.assertEqual(
            {row["baseline_id"] for row in rows},
            {"multispecies_interaction_logistic_regression_l2"},
        )
        self.assertEqual(
            {row["fold_family"] for row in rows},
            {
                "primary_genotype_or_ecotype_holdout",
                "secondary_light_treatment_holdout",
                "diagnostic_condition_stratum_holdout",
            },
        )
        self.assertEqual({row["status"] for row in rows}, {"evaluated"})
        self.assertEqual(
            {row["claim_boundary"] for row in rows},
            {"draft_interaction_logistic_only_not_leaderboard"},
        )
        for row in rows:
            self.assertTrue((REPO_ROOT / row["predictions"]).exists())
            self.assertTrue((REPO_ROOT / row["metrics"]).exists())
            self.assertTrue((REPO_ROOT / row["run_manifest"]).exists())

    def test_v9_osd120_logistic_fold_detail_comparison_api(self):
        from spacebio_bench import (
            write_multispecies_interaction_fold_detail_comparison,
            write_multispecies_interaction_fold_details,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_logistic_l2"
            logistic_outputs = write_multispecies_interaction_fold_details(
                summary_csv=(
                    REPO_ROOT
                    / "v9"
                    / "multispecies"
                    / "reports"
                    / "interaction_logistic_l2"
                    / "multispecies_baseline_summary.csv"
                ),
                repo_root=REPO_ROOT,
                output_dir=output_dir,
            )
            comparison_outputs = write_multispecies_interaction_fold_detail_comparison(
                output_dir=output_dir,
                nearest_centroid_fold_detail_csv=(
                    REPO_ROOT
                    / "v9"
                    / "multispecies"
                    / "reports"
                    / "interaction_sensitivity"
                    / "fold_detail_summary.csv"
                ),
                logistic_fold_detail_csv=logistic_outputs["csv"],
                repo_root=REPO_ROOT,
            )
            with logistic_outputs["csv"].open(newline="") as handle:
                fold_rows = list(csv.DictReader(handle))
            with comparison_outputs["csv"].open(newline="") as handle:
                comparison_rows = list(csv.DictReader(handle))

        self.assertEqual(len(fold_rows), 11)
        self.assertEqual(len(comparison_rows), 11)
        rows_by_key = {
            (row["fold_family"], row["heldout_value"]): row
            for row in comparison_rows
        }
        weak_dark = rows_by_key[
            ("diagnostic_condition_stratum_holdout", "Col.0.PhyD|Dark.Treatment")
        ]
        self.assertAlmostEqual(
            float(weak_dark["nearest_centroid_balanced_accuracy"]),
            0.5,
        )
        self.assertAlmostEqual(float(weak_dark["logistic_balanced_accuracy"]), 1.0 / 3.0)
        self.assertEqual(weak_dark["logistic_improved"], "False")
        self.assertEqual(weak_dark["logistic_new_or_persistent_le_0_5"], "True")

        light = rows_by_key[
            ("secondary_light_treatment_holdout", "Light.Treatment")
        ]
        self.assertAlmostEqual(
            float(light["nearest_centroid_balanced_accuracy"]),
            0.5555555556,
        )
        self.assertAlmostEqual(float(light["logistic_balanced_accuracy"]), 0.7222222222)
        self.assertEqual(light["logistic_improved"], "True")

    def test_build_v9_osd120_logistic_fold_comparison_cli(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_logistic_l2"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "build_v9_osd120_logistic_fold_comparison.py"),
                    "--repo-root",
                    str(REPO_ROOT),
                    "--nearest-centroid-fold-detail-csv",
                    str(
                        REPO_ROOT
                        / "v9"
                        / "multispecies"
                        / "reports"
                        / "interaction_sensitivity"
                        / "fold_detail_summary.csv"
                    ),
                    "--logistic-summary-csv",
                    str(
                        REPO_ROOT
                        / "v9"
                        / "multispecies"
                        / "reports"
                        / "interaction_logistic_l2"
                        / "multispecies_baseline_summary.csv"
                    ),
                    "--logistic-output-dir",
                    str(output_dir),
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            with (output_dir / "fold_detail_summary.csv").open(newline="") as handle:
                fold_rows = list(csv.DictReader(handle))
            with (
                output_dir / "fold_detail_comparison_vs_nearest_centroid.csv"
            ).open(newline="") as handle:
                comparison_rows = list(csv.DictReader(handle))

        self.assertIn("fold_detail_comparison_vs_nearest_centroid.csv", result.stdout)
        self.assertEqual(len(fold_rows), 11)
        self.assertEqual(len(comparison_rows), 11)

    def test_generated_v9_osd120_logistic_fold_detail_comparison_outputs_validate(self):
        output_dir = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "interaction_logistic_l2"
        )
        fold_csv = output_dir / "fold_detail_summary.csv"
        comparison_csv = output_dir / "fold_detail_comparison_vs_nearest_centroid.csv"
        comparison_json = output_dir / "fold_detail_comparison_vs_nearest_centroid.json"
        with fold_csv.open(newline="") as handle:
            fold_rows = list(csv.DictReader(handle))
        with comparison_csv.open(newline="") as handle:
            comparison_rows = list(csv.DictReader(handle))
        comparison_payload = json.loads(comparison_json.read_text())

        self.assertEqual(len(fold_rows), 11)
        self.assertEqual(len(comparison_rows), 11)
        self.assertEqual(len(comparison_payload), 11)
        self.assertEqual(
            sum(row["logistic_improved"] == "True" for row in comparison_rows),
            8,
        )
        rows_by_key = {
            (row["fold_family"], row["heldout_value"]): row
            for row in comparison_rows
        }
        self.assertAlmostEqual(
            float(
                rows_by_key[
                    (
                        "diagnostic_condition_stratum_holdout",
                        "Col.0.PhyD|Dark.Treatment",
                    )
                ]["logistic_balanced_accuracy"]
            ),
            1.0 / 3.0,
        )

    def test_v9_osd120_logistic_sensitivity_grid_api_writes_comparison(self):
        from spacebio_bench import (
            MultispeciesLogisticBaselineConfig,
            multispecies_logistic_config_id,
            run_multispecies_interaction_logistic_sensitivity_grid,
            write_multispecies_interaction_fold_detail_comparison,
            write_multispecies_interaction_fold_details,
        )

        configs = [
            MultispeciesLogisticBaselineConfig(top_variable_genes=100, c=0.1),
            MultispeciesLogisticBaselineConfig(top_variable_genes=100, c=1.0),
        ]
        self.assertEqual(
            multispecies_logistic_config_id(configs[0]),
            "tvg100_log1p_zscore_c0p1",
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_logistic_l2_sensitivity"
            _, summary = run_multispecies_interaction_logistic_sensitivity_grid(
                manifest_dir=(
                    REPO_ROOT
                    / "v9"
                    / "multispecies"
                    / "interaction_task_manifests"
                ),
                repo_root=REPO_ROOT,
                output_root=output_dir,
                configs=configs,
                fold_families=["diagnostic_condition_stratum_holdout"],
            )
            fold_outputs = write_multispecies_interaction_fold_details(
                summary_csv=summary["csv"],
                repo_root=REPO_ROOT,
                output_dir=output_dir,
            )
            comparison_outputs = write_multispecies_interaction_fold_detail_comparison(
                output_dir=output_dir,
                nearest_centroid_fold_detail_csv=(
                    REPO_ROOT
                    / "v9"
                    / "multispecies"
                    / "reports"
                    / "interaction_sensitivity"
                    / "fold_detail_summary.csv"
                ),
                logistic_fold_detail_csv=fold_outputs["csv"],
                repo_root=REPO_ROOT,
                logistic_variant_id=None,
            )
            with summary["csv"].open(newline="") as handle:
                summary_rows = list(csv.DictReader(handle))
            with fold_outputs["csv"].open(newline="") as handle:
                fold_rows = list(csv.DictReader(handle))
            with comparison_outputs["csv"].open(newline="") as handle:
                comparison_rows = list(csv.DictReader(handle))

        self.assertEqual(len(summary_rows), 2)
        self.assertEqual(len(fold_rows), 12)
        self.assertEqual(len(comparison_rows), 12)
        self.assertEqual(
            {row["logistic_variant_id"] for row in comparison_rows},
            {"tvg100_log1p_zscore_c0p1", "tvg100_log1p_zscore_c1"},
        )
        self.assertEqual(
            sum(
                row["heldout_value"] == "Col.0.PhyD|Dark.Treatment"
                for row in comparison_rows
            ),
            2,
        )

    def test_run_v9_osd120_interaction_logistic_sensitivity_cli(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_logistic_l2_sensitivity"
            result = subprocess.run(
                [
                    sys.executable,
                    str(
                        REPO_ROOT
                        / "scripts"
                        / "run_v9_osd120_interaction_logistic_sensitivity.py"
                    ),
                    "--repo-root",
                    str(REPO_ROOT),
                    "--manifest-dir",
                    str(
                        REPO_ROOT
                        / "v9"
                        / "multispecies"
                        / "interaction_task_manifests"
                    ),
                    "--nearest-centroid-fold-detail-csv",
                    str(
                        REPO_ROOT
                        / "v9"
                        / "multispecies"
                        / "reports"
                        / "interaction_sensitivity"
                        / "fold_detail_summary.csv"
                    ),
                    "--output-dir",
                    str(output_dir),
                    "--fold-family",
                    "secondary_light_treatment_holdout",
                    "--top-variable-genes",
                    "100",
                    "--c",
                    "0.1",
                    "--c",
                    "1.0",
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            with (output_dir / "multispecies_baseline_summary.csv").open(
                newline=""
            ) as handle:
                summary_rows = list(csv.DictReader(handle))
            with (output_dir / "fold_detail_summary.csv").open(newline="") as handle:
                fold_rows = list(csv.DictReader(handle))
            with (
                output_dir / "fold_detail_comparison_vs_nearest_centroid.csv"
            ).open(newline="") as handle:
                comparison_rows = list(csv.DictReader(handle))

        self.assertIn("fold_detail_comparison_vs_nearest_centroid.csv", result.stdout)
        self.assertEqual(len(summary_rows), 2)
        self.assertEqual(len(fold_rows), 4)
        self.assertEqual(len(comparison_rows), 4)

    def test_generated_v9_osd120_logistic_sensitivity_outputs_validate(self):
        output_dir = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "interaction_logistic_l2_sensitivity"
        )
        with (output_dir / "multispecies_baseline_summary.csv").open(
            newline=""
        ) as handle:
            summary_rows = list(csv.DictReader(handle))
        with (output_dir / "fold_detail_summary.csv").open(newline="") as handle:
            fold_rows = list(csv.DictReader(handle))
        with (
            output_dir / "fold_detail_comparison_vs_nearest_centroid.csv"
        ).open(newline="") as handle:
            comparison_rows = list(csv.DictReader(handle))

        self.assertEqual(len(summary_rows), 18)
        self.assertEqual(len(fold_rows), 66)
        self.assertEqual(len(comparison_rows), 66)
        self.assertEqual(
            {row["variant_id"] for row in summary_rows},
            {
                "tvg500_log1p_zscore_c0p1",
                "tvg500_log1p_zscore_c1",
                "tvg500_log1p_zscore_c10",
                "tvg2000_log1p_zscore_c0p1",
                "tvg2000_log1p_zscore_c1",
                "tvg2000_log1p_zscore_c10",
            },
        )
        self.assertEqual(
            sum(
                row["heldout_value"] == "Col.0.PhyD|Dark.Treatment"
                for row in comparison_rows
            ),
            6,
        )
        rows_by_variant_and_fold = {
            (row["logistic_variant_id"], row["heldout_value"]): row
            for row in comparison_rows
        }
        for c_token in ("c0p1", "c1", "c10"):
            self.assertAlmostEqual(
                float(
                    rows_by_variant_and_fold[
                        (f"tvg2000_log1p_zscore_{c_token}", "Col.0.PhyD|Dark.Treatment")
                    ]["logistic_balanced_accuracy"]
                ),
                1.0 / 3.0,
            )
            self.assertAlmostEqual(
                float(
                    rows_by_variant_and_fold[
                        (f"tvg500_log1p_zscore_{c_token}", "Col.0.PhyD|Dark.Treatment")
                    ]["logistic_balanced_accuracy"]
                ),
                2.0 / 3.0,
            )

    def test_v9_osd120_logistic_feature_audit_api(self):
        from spacebio_bench import (
            MultispeciesLogisticBaselineConfig,
            write_multispecies_interaction_logistic_feature_audit,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_logistic_feature_audit"
            outputs = write_multispecies_interaction_logistic_feature_audit(
                manifest_dir=(
                    REPO_ROOT
                    / "v9"
                    / "multispecies"
                    / "interaction_task_manifests"
                ),
                repo_root=REPO_ROOT,
                output_dir=output_dir,
                configs=[
                    MultispeciesLogisticBaselineConfig(top_variable_genes=25, c=1.0),
                    MultispeciesLogisticBaselineConfig(top_variable_genes=50, c=1.0),
                ],
                focus_folds=[
                    (
                        "diagnostic_condition_stratum_holdout",
                        "Col.0.PhyD|Dark.Treatment",
                    )
                ],
                top_n_coefficients=5,
            )
            with outputs["summary_csv"].open(newline="") as handle:
                summary_rows = list(csv.DictReader(handle))
            with outputs["coefficient_csv"].open(newline="") as handle:
                coefficient_rows = list(csv.DictReader(handle))

        self.assertEqual(len(summary_rows), 2)
        self.assertEqual(len(coefficient_rows), 75)
        self.assertEqual(
            {row["variant_id"] for row in summary_rows},
            {"tvg25_log1p_zscore_c1", "tvg50_log1p_zscore_c1"},
        )
        self.assertEqual({row["n_extra_features"] for row in summary_rows}, {"0"})
        self.assertTrue(summary_rows[0]["top_abs_feature_ids"])

    def test_audit_v9_osd120_logistic_feature_sets_cli(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_logistic_feature_audit"
            result = subprocess.run(
                [
                    sys.executable,
                    str(
                        REPO_ROOT
                        / "scripts"
                        / "audit_v9_osd120_logistic_feature_sets.py"
                    ),
                    "--repo-root",
                    str(REPO_ROOT),
                    "--manifest-dir",
                    str(
                        REPO_ROOT
                        / "v9"
                        / "multispecies"
                        / "interaction_task_manifests"
                    ),
                    "--output-dir",
                    str(output_dir),
                    "--focus-fold",
                    "secondary_light_treatment_holdout=Light.Treatment",
                    "--top-variable-genes",
                    "25",
                    "--top-variable-genes",
                    "50",
                    "--top-n-coefficients",
                    "5",
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            with (output_dir / "feature_set_audit_summary.csv").open(
                newline=""
            ) as handle:
                summary_rows = list(csv.DictReader(handle))
            with (output_dir / "feature_coefficient_audit.csv").open(
                newline=""
            ) as handle:
                coefficient_rows = list(csv.DictReader(handle))

        self.assertIn("feature_coefficient_audit.csv", result.stdout)
        self.assertEqual(len(summary_rows), 2)
        self.assertEqual(len(coefficient_rows), 75)

    def test_generated_v9_osd120_logistic_feature_audit_outputs_validate(self):
        output_dir = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "interaction_logistic_feature_audit"
        )
        with (output_dir / "feature_set_audit_summary.csv").open(
            newline=""
        ) as handle:
            summary_rows = list(csv.DictReader(handle))
        with (output_dir / "feature_coefficient_audit.csv").open(
            newline=""
        ) as handle:
            coefficient_rows = list(csv.DictReader(handle))
        summary_payload = json.loads(
            (output_dir / "feature_set_audit_summary.json").read_text()
        )

        self.assertEqual(len(summary_rows), 6)
        self.assertEqual(len(summary_payload), 6)
        self.assertEqual(len(coefficient_rows), 7500)
        rows_by_variant_and_fold = {
            (row["variant_id"], row["heldout_value"]): row for row in summary_rows
        }
        dark_top500 = rows_by_variant_and_fold[
            ("tvg500_log1p_zscore_c1", "Col.0.PhyD|Dark.Treatment")
        ]
        dark_top2000 = rows_by_variant_and_fold[
            ("tvg2000_log1p_zscore_c1", "Col.0.PhyD|Dark.Treatment")
        ]
        self.assertAlmostEqual(float(dark_top500["balanced_accuracy"]), 2.0 / 3.0)
        self.assertAlmostEqual(float(dark_top2000["balanced_accuracy"]), 1.0 / 3.0)
        self.assertEqual(dark_top500["n_extra_features"], "0")
        self.assertEqual(dark_top2000["n_extra_features"], "1500")
        self.assertTrue(dark_top2000["top_abs_extra_feature_ids"])

    def test_run_v9_osd120_interaction_sparse_l1_cli(self):
        from spacebio_bench import (
            MultispeciesLogisticBaselineConfig,
            multispecies_logistic_config_id,
        )

        self.assertEqual(
            multispecies_logistic_config_id(
                MultispeciesLogisticBaselineConfig(
                    top_variable_genes=100,
                    c=0.1,
                    penalty="l1",
                    solver="liblinear",
                )
            ),
            "tvg100_log1p_zscore_l1_c0p1",
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_logistic_sparse_l1"
            result = subprocess.run(
                [
                    sys.executable,
                    str(
                        REPO_ROOT
                        / "scripts"
                        / "run_v9_osd120_interaction_sparse_l1.py"
                    ),
                    "--repo-root",
                    str(REPO_ROOT),
                    "--manifest-dir",
                    str(
                        REPO_ROOT
                        / "v9"
                        / "multispecies"
                        / "interaction_task_manifests"
                    ),
                    "--nearest-centroid-fold-detail-csv",
                    str(
                        REPO_ROOT
                        / "v9"
                        / "multispecies"
                        / "reports"
                        / "interaction_sensitivity"
                        / "fold_detail_summary.csv"
                    ),
                    "--output-dir",
                    str(output_dir),
                    "--fold-family",
                    "secondary_light_treatment_holdout",
                    "--focus-fold",
                    "secondary_light_treatment_holdout=Light.Treatment",
                    "--top-variable-genes",
                    "100",
                    "--c",
                    "0.1",
                    "--top-n-coefficients",
                    "5",
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            with (output_dir / "multispecies_baseline_summary.csv").open(
                newline=""
            ) as handle:
                summary_rows = list(csv.DictReader(handle))
            with (output_dir / "fold_detail_summary.csv").open(newline="") as handle:
                fold_rows = list(csv.DictReader(handle))
            with (
                output_dir / "fold_detail_comparison_vs_nearest_centroid.csv"
            ).open(newline="") as handle:
                comparison_rows = list(csv.DictReader(handle))
            with (output_dir / "feature_set_audit_summary.csv").open(
                newline=""
            ) as handle:
                audit_rows = list(csv.DictReader(handle))
            with (output_dir / "feature_coefficient_audit.csv").open(
                newline=""
            ) as handle:
                coefficient_rows = list(csv.DictReader(handle))

        self.assertIn("feature_coefficient_audit.csv", result.stdout)
        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(len(fold_rows), 2)
        self.assertEqual(len(comparison_rows), 2)
        self.assertEqual(len(audit_rows), 1)
        self.assertEqual(len(coefficient_rows), 100)
        self.assertEqual(audit_rows[0]["variant_id"], "tvg100_log1p_zscore_l1_c0p1")

    def test_generated_v9_osd120_sparse_l1_outputs_validate(self):
        output_dir = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "interaction_logistic_sparse_l1"
        )
        with (output_dir / "multispecies_baseline_summary.csv").open(
            newline=""
        ) as handle:
            summary_rows = list(csv.DictReader(handle))
        with (output_dir / "fold_detail_summary.csv").open(newline="") as handle:
            fold_rows = list(csv.DictReader(handle))
        with (
            output_dir / "fold_detail_comparison_vs_nearest_centroid.csv"
        ).open(newline="") as handle:
            comparison_rows = list(csv.DictReader(handle))
        with (output_dir / "feature_set_audit_summary.csv").open(
            newline=""
        ) as handle:
            audit_rows = list(csv.DictReader(handle))
        with (output_dir / "feature_coefficient_audit.csv").open(
            newline=""
        ) as handle:
            coefficient_rows = list(csv.DictReader(handle))

        self.assertEqual(len(summary_rows), 15)
        self.assertEqual(len(fold_rows), 55)
        self.assertEqual(len(comparison_rows), 55)
        self.assertEqual(len(audit_rows), 15)
        self.assertEqual(len(coefficient_rows), 30000)
        self.assertEqual(
            {row["variant_id"] for row in summary_rows},
            {
                "tvg2000_log1p_zscore_l1_c0p01",
                "tvg2000_log1p_zscore_l1_c0p03",
                "tvg2000_log1p_zscore_l1_c0p1",
                "tvg2000_log1p_zscore_l1_c0p3",
                "tvg2000_log1p_zscore_l1_c1",
            },
        )
        self.assertTrue(
            all(int(row["n_nonzero_coefficients"]) <= 2000 for row in audit_rows)
        )
        audit_by_variant_and_fold = {
            (row["variant_id"], row["heldout_value"]): row for row in audit_rows
        }
        c1_dark = audit_by_variant_and_fold[
            ("tvg2000_log1p_zscore_l1_c1", "Col.0.PhyD|Dark.Treatment")
        ]
        c1_light = audit_by_variant_and_fold[
            ("tvg2000_log1p_zscore_l1_c1", "Light.Treatment")
        ]
        c1_phyd = audit_by_variant_and_fold[
            ("tvg2000_log1p_zscore_l1_c1", "Col.0.PhyD")
        ]
        self.assertAlmostEqual(float(c1_dark["balanced_accuracy"]), 2.0 / 3.0)
        self.assertAlmostEqual(float(c1_light["balanced_accuracy"]), 0.8333333333)
        self.assertAlmostEqual(float(c1_phyd["balanced_accuracy"]), 0.9166666667)
        self.assertGreater(int(c1_dark["n_nonzero_coefficients"]), 0)

    def test_v9_osd120_sparse_l1_stability_audit_api(self):
        from spacebio_bench import (
            MultispeciesLogisticBaselineConfig,
            write_multispecies_interaction_sparse_l1_stability_audit,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_logistic_sparse_l1_stability"
            outputs = write_multispecies_interaction_sparse_l1_stability_audit(
                manifest_dir=(
                    REPO_ROOT
                    / "v9"
                    / "multispecies"
                    / "interaction_task_manifests"
                ),
                repo_root=REPO_ROOT,
                output_dir=output_dir,
                configs=[
                    MultispeciesLogisticBaselineConfig(
                        top_variable_genes=100,
                        c=0.3,
                        penalty="l1",
                        solver="liblinear",
                    )
                ],
                focus_folds=[
                    (
                        "diagnostic_condition_stratum_holdout",
                        "Col.0.PhyD|Dark.Treatment",
                    )
                ],
                n_subsamples=4,
                subsample_fraction=0.75,
            )
            with outputs["summary_csv"].open(newline="") as handle:
                summary_rows = list(csv.DictReader(handle))
            with outputs["feature_csv"].open(newline="") as handle:
                feature_rows = list(csv.DictReader(handle))

        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(len(feature_rows), 100)
        self.assertEqual(summary_rows[0]["n_subsamples"], "4")
        self.assertEqual(summary_rows[0]["variant_id"], "tvg100_log1p_zscore_l1_c0p3")
        self.assertTrue(0.0 <= float(summary_rows[0]["mean_pairwise_jaccard"]) <= 1.0)

    def test_audit_v9_osd120_sparse_l1_stability_cli(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_logistic_sparse_l1_stability"
            result = subprocess.run(
                [
                    sys.executable,
                    str(
                        REPO_ROOT
                        / "scripts"
                        / "audit_v9_osd120_sparse_l1_stability.py"
                    ),
                    "--repo-root",
                    str(REPO_ROOT),
                    "--manifest-dir",
                    str(
                        REPO_ROOT
                        / "v9"
                        / "multispecies"
                        / "interaction_task_manifests"
                    ),
                    "--output-dir",
                    str(output_dir),
                    "--focus-fold",
                    "secondary_light_treatment_holdout=Light.Treatment",
                    "--top-variable-genes",
                    "100",
                    "--c",
                    "0.3",
                    "--n-subsamples",
                    "4",
                    "--subsample-fraction",
                    "0.75",
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            with (output_dir / "stability_summary.csv").open(newline="") as handle:
                summary_rows = list(csv.DictReader(handle))
            with (
                output_dir / "stability_feature_frequency.csv"
            ).open(newline="") as handle:
                feature_rows = list(csv.DictReader(handle))

        self.assertIn("stability_feature_frequency.csv", result.stdout)
        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(len(feature_rows), 100)

    def test_generated_v9_osd120_sparse_l1_stability_outputs_validate(self):
        output_dir = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "interaction_logistic_sparse_l1_stability"
        )
        with (output_dir / "stability_summary.csv").open(newline="") as handle:
            summary_rows = list(csv.DictReader(handle))
        with (output_dir / "stability_feature_frequency.csv").open(
            newline=""
        ) as handle:
            feature_rows = list(csv.DictReader(handle))

        self.assertEqual(len(summary_rows), 6)
        self.assertEqual(len(feature_rows), 12000)
        self.assertEqual({row["n_subsamples"] for row in summary_rows}, {"20"})
        self.assertEqual(
            {row["variant_id"] for row in summary_rows},
            {
                "tvg2000_log1p_zscore_l1_c0p3",
                "tvg2000_log1p_zscore_l1_c1",
            },
        )
        self.assertTrue(
            all(0.0 <= float(row["mean_pairwise_jaccard"]) <= 1.0 for row in summary_rows)
        )
        rows_by_variant_and_fold = {
            (row["variant_id"], row["heldout_value"]): row for row in summary_rows
        }
        dark_c1 = rows_by_variant_and_fold[
            ("tvg2000_log1p_zscore_l1_c1", "Col.0.PhyD|Dark.Treatment")
        ]
        light_c03 = rows_by_variant_and_fold[
            ("tvg2000_log1p_zscore_l1_c0p3", "Light.Treatment")
        ]
        self.assertEqual(dark_c1["stable_feature_count_ge_0_5"], "10")
        self.assertEqual(dark_c1["stable_feature_count_ge_0_8"], "7")
        self.assertIn("AT1G54970:1", dark_c1["top_selection_frequency_features"])
        self.assertEqual(light_c03["reference_nonzero_count"], "3")
        self.assertIn("AT5G04700:1", light_c03["top_selection_frequency_features"])

    def test_v9_osd120_baseline_ladder_api(self):
        from spacebio_bench import write_multispecies_interaction_baseline_ladder

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_baseline_ladder"
            outputs = write_multispecies_interaction_baseline_ladder(
                repo_root=REPO_ROOT,
                output_dir=output_dir,
            )
            with outputs["summary_csv"].open(newline="") as handle:
                summary_rows = list(csv.DictReader(handle))
            with outputs["focus_csv"].open(newline="") as handle:
                focus_rows = list(csv.DictReader(handle))
            summary_payload = json.loads(outputs["summary_json"].read_text())
            focus_payload = json.loads(outputs["focus_json"].read_text())

        self.assertEqual(len(summary_rows), 5)
        self.assertEqual(len(focus_rows), 15)
        self.assertEqual(len(summary_payload), 5)
        self.assertEqual(len(focus_payload), 15)
        rows_by_candidate = {row["candidate_id"]: row for row in summary_rows}
        c1 = rows_by_candidate["sparse_l1_c1"]
        self.assertEqual(c1["decision"], "advance_as_primary_sparse_diagnostic_candidate")
        self.assertEqual(c1["nearest_worse_count"], "0")
        self.assertAlmostEqual(float(c1["primary_balanced_accuracy"]), 0.9166666667)
        self.assertEqual(c1["dark_stable_feature_count_ge_0_5"], "10")

    def test_build_v9_osd120_baseline_ladder_cli(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_baseline_ladder"
            result = subprocess.run(
                [
                    sys.executable,
                    str(
                        REPO_ROOT
                        / "scripts"
                        / "build_v9_osd120_baseline_ladder.py"
                    ),
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
            with (output_dir / "baseline_ladder_summary.csv").open(
                newline=""
            ) as handle:
                summary_rows = list(csv.DictReader(handle))
            with (output_dir / "baseline_ladder_focus_folds.csv").open(
                newline=""
            ) as handle:
                focus_rows = list(csv.DictReader(handle))

        self.assertIn("baseline_ladder_focus_folds.csv", result.stdout)
        self.assertEqual(len(summary_rows), 5)
        self.assertEqual(len(focus_rows), 15)

    def test_generated_v9_osd120_baseline_ladder_outputs_validate(self):
        output_dir = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "interaction_baseline_ladder"
        )
        with (output_dir / "baseline_ladder_summary.csv").open(newline="") as handle:
            summary_rows = list(csv.DictReader(handle))
        with (output_dir / "baseline_ladder_focus_folds.csv").open(
            newline=""
        ) as handle:
            focus_rows = list(csv.DictReader(handle))
        summary_payload = json.loads(
            (output_dir / "baseline_ladder_summary.json").read_text()
        )
        focus_payload = json.loads(
            (output_dir / "baseline_ladder_focus_folds.json").read_text()
        )

        self.assertEqual(len(summary_rows), 5)
        self.assertEqual(len(focus_rows), 15)
        self.assertEqual(len(summary_payload), 5)
        self.assertEqual(len(focus_payload), 15)
        rows_by_candidate = {row["candidate_id"]: row for row in summary_rows}
        nearest = rows_by_candidate["nearest_centroid_default"]
        c03 = rows_by_candidate["sparse_l1_c0p3"]
        c1 = rows_by_candidate["sparse_l1_c1"]
        self.assertAlmostEqual(float(nearest["mean_family_balanced_accuracy"]), 2.0 / 3.0)
        self.assertEqual(nearest["nearest_tied_count"], "11")
        self.assertAlmostEqual(float(c03["diagnostic_balanced_accuracy"]), 0.9166666667)
        self.assertAlmostEqual(float(c1["primary_balanced_accuracy"]), 0.9166666667)
        self.assertAlmostEqual(float(c1["secondary_balanced_accuracy"]), 0.8333333333)
        self.assertAlmostEqual(float(c1["diagnostic_balanced_accuracy"]), 0.8888888889)
        self.assertEqual(c1["nearest_worse_count"], "0")
        self.assertEqual(c1["dark_stable_feature_count_ge_0_5"], "10")
        focus_by_candidate_and_value = {
            (row["candidate_id"], row["heldout_value"]): row for row in focus_rows
        }
        c1_dark = focus_by_candidate_and_value[
            ("sparse_l1_c1", "Col.0.PhyD|Dark.Treatment")
        ]
        c1_light = focus_by_candidate_and_value[
            ("sparse_l1_c1", "Light.Treatment")
        ]
        self.assertAlmostEqual(float(c1_dark["candidate_balanced_accuracy"]), 2.0 / 3.0)
        self.assertAlmostEqual(float(c1_light["delta_candidate_minus_nearest_centroid"]), 0.2777777777)

    def test_v9_osd120_diagnostic_candidate_package_api(self):
        from spacebio_bench import (
            write_multispecies_interaction_diagnostic_candidate_package,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_diagnostic_candidate_package"
            outputs = write_multispecies_interaction_diagnostic_candidate_package(
                repo_root=REPO_ROOT,
                output_dir=output_dir,
            )
            with outputs["summary_csv"].open(newline="") as handle:
                summary_rows = list(csv.DictReader(handle))
            with outputs["focus_csv"].open(newline="") as handle:
                focus_rows = list(csv.DictReader(handle))
            with outputs["feature_csv"].open(newline="") as handle:
                feature_rows = list(csv.DictReader(handle))
            with outputs["claim_csv"].open(newline="") as handle:
                claim_rows = list(csv.DictReader(handle))

        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(len(focus_rows), 3)
        self.assertEqual(len(feature_rows), 19)
        self.assertEqual(len(claim_rows), 5)
        summary = summary_rows[0]
        self.assertEqual(summary["candidate_id"], "sparse_l1_c1")
        self.assertEqual(summary["stable_feature_count_ge_0_5_total"], "19")
        self.assertEqual(summary["nearest_worse_count"], "0")
        self.assertIn(
            "https://doi.org/10.1371/journal.pone.0180186",
            next(row for row in claim_rows if row["claim_id"] == "external_light_genotype_context")[
                "external_source_urls"
            ],
        )

    def test_build_v9_osd120_diagnostic_candidate_package_cli(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_diagnostic_candidate_package"
            result = subprocess.run(
                [
                    sys.executable,
                    str(
                        REPO_ROOT
                        / "scripts"
                        / "build_v9_osd120_diagnostic_candidate_package.py"
                    ),
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
            with (output_dir / "candidate_package_summary.csv").open(
                newline=""
            ) as handle:
                summary_rows = list(csv.DictReader(handle))
            with (output_dir / "candidate_stable_feature_evidence.csv").open(
                newline=""
            ) as handle:
                feature_rows = list(csv.DictReader(handle))

        self.assertIn("candidate_claim_map.csv", result.stdout)
        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(len(feature_rows), 19)

    def test_generated_v9_osd120_diagnostic_candidate_package_outputs_validate(self):
        output_dir = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "interaction_diagnostic_candidate_package"
        )
        with (output_dir / "candidate_package_summary.csv").open(newline="") as handle:
            summary_rows = list(csv.DictReader(handle))
        with (output_dir / "candidate_focus_evidence.csv").open(newline="") as handle:
            focus_rows = list(csv.DictReader(handle))
        with (output_dir / "candidate_stable_feature_evidence.csv").open(
            newline=""
        ) as handle:
            feature_rows = list(csv.DictReader(handle))
        with (output_dir / "candidate_claim_map.csv").open(newline="") as handle:
            claim_rows = list(csv.DictReader(handle))
        summary_payload = json.loads(
            (output_dir / "candidate_package_summary.json").read_text()
        )

        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(len(summary_payload), 1)
        self.assertEqual(len(focus_rows), 3)
        self.assertEqual(len(feature_rows), 19)
        self.assertEqual(len(claim_rows), 5)
        summary = summary_rows[0]
        self.assertEqual(summary["package_id"], "osd120_sparse_l1_c1_draft_candidate")
        self.assertEqual(summary["candidate_variant_id"], "tvg2000_log1p_zscore_l1_c1")
        self.assertEqual(summary["focus_improved_count"], "3")
        self.assertAlmostEqual(float(summary["focus_min_candidate_ba"]), 2.0 / 3.0)
        focus_by_value = {row["heldout_value"]: row for row in focus_rows}
        self.assertEqual(
            focus_by_value["Col.0.PhyD|Dark.Treatment"]["stable_feature_count_ge_0_8"],
            "7",
        )
        self.assertEqual(
            focus_by_value["Light.Treatment"]["n_nonzero_extra_coefficients"],
            "8",
        )
        self.assertEqual(
            Counter(row["stability_tier"] for row in feature_rows),
            Counter({"stable_ge_0_8": 10, "stable_ge_0_5": 9}),
        )

    def test_v9_osd120_figure_table_package_api(self):
        from spacebio_bench import write_multispecies_interaction_figure_table_package

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_figure_table_package"
            outputs = write_multispecies_interaction_figure_table_package(
                repo_root=REPO_ROOT,
                output_dir=output_dir,
            )
            with outputs["main_csv"].open(newline="") as handle:
                main_rows = list(csv.DictReader(handle))
            with outputs["feature_csv"].open(newline="") as handle:
                feature_rows = list(csv.DictReader(handle))
            caption = outputs["caption_md"].read_text()
            claim_box = outputs["claim_box_md"].read_text()

        self.assertEqual(len(main_rows), 3)
        self.assertEqual(len(feature_rows), 19)
        self.assertEqual(main_rows[0]["display_delta_ba"], "+0.167")
        self.assertIn("not validated biomarkers", caption)
        self.assertIn("Not Allowed", claim_box)

    def test_build_v9_osd120_figure_table_package_cli(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_figure_table_package"
            result = subprocess.run(
                [
                    sys.executable,
                    str(
                        REPO_ROOT
                        / "scripts"
                        / "build_v9_osd120_figure_table_package.py"
                    ),
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
            with (output_dir / "figure_main_focus_table.csv").open(
                newline=""
            ) as handle:
                main_rows = list(csv.DictReader(handle))
            caption = (output_dir / "figure_caption.md").read_text()

        self.assertIn("figure_claim_boundary_box.md", result.stdout)
        self.assertEqual(len(main_rows), 3)
        self.assertIn("within-source OSD-120 interaction split", caption)

    def test_generated_v9_osd120_figure_table_package_outputs_validate(self):
        output_dir = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "interaction_figure_table_package"
        )
        with (output_dir / "figure_main_focus_table.csv").open(newline="") as handle:
            main_rows = list(csv.DictReader(handle))
        with (output_dir / "figure_stable_feature_appendix.csv").open(
            newline=""
        ) as handle:
            feature_rows = list(csv.DictReader(handle))
        main_payload = json.loads(
            (output_dir / "figure_main_focus_table.json").read_text()
        )
        feature_payload = json.loads(
            (output_dir / "figure_stable_feature_appendix.json").read_text()
        )
        caption = (output_dir / "figure_caption.md").read_text()
        claim_box = (output_dir / "figure_claim_boundary_box.md").read_text()

        self.assertEqual(len(main_rows), 3)
        self.assertEqual(len(main_payload), 3)
        self.assertEqual(len(feature_rows), 19)
        self.assertEqual(len(feature_payload), 19)
        rows_by_fold = {row["fold_label"]: row for row in main_rows}
        self.assertEqual(
            rows_by_fold["primary_col0_phyd"]["candidate_ba"],
            "0.917",
        )
        self.assertEqual(
            rows_by_fold["diagnostic_col0_phyd_dark"]["stable_features_ge_0_8"],
            "7",
        )
        self.assertEqual(
            Counter(row["coefficient_direction"] for row in feature_rows),
            Counter({"positive_LEO_or_ISS": 12, "negative_LEO_or_ISS": 7}),
        )
        self.assertIn("Disallowed claims", caption)
        self.assertIn("Do not call this a frozen v9 benchmark baseline", claim_box)

    def test_v9_osd120_paired_comparator_table_api(self):
        from spacebio_bench import write_multispecies_interaction_paired_comparator_table

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_paired_comparator_table"
            outputs = write_multispecies_interaction_paired_comparator_table(
                repo_root=REPO_ROOT,
                output_dir=output_dir,
            )
            with outputs["summary_csv"].open(newline="") as handle:
                summary_rows = list(csv.DictReader(handle))
            with outputs["focus_csv"].open(newline="") as handle:
                focus_rows = list(csv.DictReader(handle))
            decision = outputs["decision_md"].read_text()

        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(len(focus_rows), 3)
        summary = summary_rows[0]
        self.assertEqual(summary["focus_tied_ba_count"], "3")
        self.assertEqual(summary["recommendation"], "appendix_or_supplement_only")
        self.assertEqual(summary["primary_nearest_worse_count"], "0")
        self.assertEqual(summary["compact_nearest_worse_count"], "1")
        self.assertIn("not a second primary panel", decision)

    def test_build_v9_osd120_paired_comparator_table_cli(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_paired_comparator_table"
            result = subprocess.run(
                [
                    sys.executable,
                    str(
                        REPO_ROOT
                        / "scripts"
                        / "build_v9_osd120_paired_comparator_table.py"
                    ),
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
            with (output_dir / "paired_focus_comparator_table.csv").open(
                newline=""
            ) as handle:
                focus_rows = list(csv.DictReader(handle))

        self.assertIn("paired_comparator_decision.md", result.stdout)
        self.assertEqual(len(focus_rows), 3)
        self.assertTrue(
            all(row["display_primary_minus_compact_ba"] == "+0.000" for row in focus_rows)
        )

    def test_generated_v9_osd120_paired_comparator_table_outputs_validate(self):
        output_dir = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "interaction_paired_comparator_table"
        )
        with (output_dir / "paired_comparator_summary.csv").open(newline="") as handle:
            summary_rows = list(csv.DictReader(handle))
        with (output_dir / "paired_focus_comparator_table.csv").open(
            newline=""
        ) as handle:
            focus_rows = list(csv.DictReader(handle))
        summary_payload = json.loads(
            (output_dir / "paired_comparator_summary.json").read_text()
        )
        focus_payload = json.loads(
            (output_dir / "paired_focus_comparator_table.json").read_text()
        )
        decision = (output_dir / "paired_comparator_decision.md").read_text()

        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(len(summary_payload), 1)
        self.assertEqual(len(focus_rows), 3)
        self.assertEqual(len(focus_payload), 3)
        summary = summary_rows[0]
        self.assertEqual(summary["primary_focus_nonzero_total"], "32")
        self.assertEqual(summary["compact_focus_nonzero_total"], "19")
        self.assertEqual(summary["primary_stable_ge_0_5_total"], "19")
        self.assertEqual(summary["compact_stable_ge_0_5_total"], "8")
        self.assertEqual(summary["focus_tied_ba_count"], "3")
        rows_by_fold = {row["fold_label"]: row for row in focus_rows}
        self.assertEqual(
            rows_by_fold["secondary_light_treatment"][
                "compactness_delta_nonzero_compact_minus_primary"
            ],
            "-7",
        )
        self.assertIn("Appendix or supplement", decision)

    def test_v9_osd120_diagnostic_artifact_manifest_api(self):
        from spacebio_bench import (
            write_multispecies_interaction_diagnostic_artifact_manifest,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_diagnostic_artifact_manifest"
            outputs = write_multispecies_interaction_diagnostic_artifact_manifest(
                repo_root=REPO_ROOT,
                output_dir=output_dir,
            )
            with outputs["artifact_csv"].open(newline="") as handle:
                artifact_rows = list(csv.DictReader(handle))
            with outputs["claim_csv"].open(newline="") as handle:
                claim_rows = list(csv.DictReader(handle))

        self.assertEqual(len(artifact_rows), 26)
        self.assertEqual(len(claim_rows), 7)
        self.assertTrue(all(row["exists"] == "True" for row in artifact_rows))
        self.assertTrue(all(row["sha256"] for row in artifact_rows))
        claim_by_id = {row["claim_id"]: row for row in claim_rows}
        self.assertIn(
            "figure_main_focus_table",
            claim_by_id["fragile_focus_recovery"]["artifact_ids"],
        )
        self.assertIn(
            "paired_comparator_summary",
            claim_by_id["paired_comparator_decision"]["artifact_ids"],
        )

    def test_build_v9_osd120_diagnostic_artifact_manifest_cli(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_diagnostic_artifact_manifest"
            result = subprocess.run(
                [
                    sys.executable,
                    str(
                        REPO_ROOT
                        / "scripts"
                        / "build_v9_osd120_diagnostic_artifact_manifest.py"
                    ),
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
            with (output_dir / "diagnostic_artifact_manifest.csv").open(
                newline=""
            ) as handle:
                artifact_rows = list(csv.DictReader(handle))
            with (output_dir / "diagnostic_claim_artifact_map.csv").open(
                newline=""
            ) as handle:
                claim_rows = list(csv.DictReader(handle))

        self.assertIn("diagnostic_claim_artifact_map.json", result.stdout)
        self.assertEqual(len(artifact_rows), 26)
        self.assertEqual(len(claim_rows), 7)

    def test_generated_v9_osd120_diagnostic_artifact_manifest_outputs_validate(self):
        output_dir = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "interaction_diagnostic_artifact_manifest"
        )
        with (output_dir / "diagnostic_artifact_manifest.csv").open(
            newline=""
        ) as handle:
            artifact_rows = list(csv.DictReader(handle))
        with (output_dir / "diagnostic_claim_artifact_map.csv").open(
            newline=""
        ) as handle:
            claim_rows = list(csv.DictReader(handle))
        artifact_payload = json.loads(
            (output_dir / "diagnostic_artifact_manifest.json").read_text()
        )
        claim_payload = json.loads(
            (output_dir / "diagnostic_claim_artifact_map.json").read_text()
        )

        self.assertEqual(len(artifact_rows), 26)
        self.assertEqual(len(artifact_payload), 26)
        self.assertEqual(len(claim_rows), 7)
        self.assertEqual(len(claim_payload), 7)
        self.assertTrue(all(row["exists"] == "True" for row in artifact_rows))
        artifacts_by_id = {row["artifact_id"]: row for row in artifact_rows}
        self.assertEqual(
            artifacts_by_id["candidate_stable_feature_evidence"]["row_count"],
            "19",
        )
        self.assertEqual(
            artifacts_by_id["figure_main_focus_table"]["row_count"],
            "3",
        )
        self.assertEqual(
            artifacts_by_id["paired_focus_comparator_table"]["row_count"],
            "3",
        )
        self.assertEqual(
            artifacts_by_id["v9_test_source"]["validation_scope"],
            "unit_tests",
        )
        claim_by_id = {row["claim_id"]: row for row in claim_rows}
        self.assertIn(
            "https://doi.org/10.1371/journal.pone.0180186",
            claim_by_id["external_light_genotype_context"][
                "external_source_urls"
            ],
        )
        self.assertIn(
            "/usr/local/bin/python3 -m unittest discover -s tests",
            claim_by_id["draft_candidate_boundary"]["validation_tests"],
        )

    def test_v9_osd120_release_readiness_gap_audit_api(self):
        from spacebio_bench import (
            write_multispecies_interaction_release_readiness_gap_audit,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_release_readiness_gap_audit"
            outputs = write_multispecies_interaction_release_readiness_gap_audit(
                repo_root=REPO_ROOT,
                output_dir=output_dir,
            )
            with outputs["summary_csv"].open(newline="") as handle:
                summary_rows = list(csv.DictReader(handle))
            with outputs["gap_csv"].open(newline="") as handle:
                gap_rows = list(csv.DictReader(handle))
            with outputs["reference_csv"].open(newline="") as handle:
                reference_rows = list(csv.DictReader(handle))

        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(len(gap_rows), 12)
        self.assertEqual(len(reference_rows), 7)
        summary = summary_rows[0]
        self.assertEqual(summary["public_alpha_ready"], "False")
        self.assertEqual(summary["blocker_count"], "3")
        self.assertEqual(summary["artifact_count"], "26")
        self.assertEqual(summary["missing_artifact_count"], "0")
        gap_by_id = {row["gap_id"]: row for row in gap_rows}
        self.assertEqual(
            gap_by_id["full_osdr_payload_freeze"]["readiness_status"],
            "blocker",
        )
        self.assertEqual(
            gap_by_id["artifact_manifest_traceability"]["readiness_status"],
            "pass",
        )
        self.assertIn(
            "ro_crate_1_2",
            {row["reference_id"] for row in reference_rows},
        )

    def test_build_v9_osd120_release_readiness_gap_audit_cli(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_release_readiness_gap_audit"
            result = subprocess.run(
                [
                    sys.executable,
                    str(
                        REPO_ROOT
                        / "scripts"
                        / "build_v9_osd120_release_readiness_gap_audit.py"
                    ),
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
            with (output_dir / "release_readiness_summary.csv").open(
                newline=""
            ) as handle:
                summary_rows = list(csv.DictReader(handle))
            with (output_dir / "release_readiness_gap_table.csv").open(
                newline=""
            ) as handle:
                gap_rows = list(csv.DictReader(handle))

        self.assertIn("release_readiness_external_references.json", result.stdout)
        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(len(gap_rows), 12)

    def test_generated_v9_osd120_release_readiness_gap_audit_outputs_validate(self):
        output_dir = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "interaction_release_readiness_gap_audit"
        )
        with (output_dir / "release_readiness_summary.csv").open(
            newline=""
        ) as handle:
            summary_rows = list(csv.DictReader(handle))
        with (output_dir / "release_readiness_gap_table.csv").open(
            newline=""
        ) as handle:
            gap_rows = list(csv.DictReader(handle))
        with (output_dir / "release_readiness_external_references.csv").open(
            newline=""
        ) as handle:
            reference_rows = list(csv.DictReader(handle))
        summary_payload = json.loads(
            (output_dir / "release_readiness_summary.json").read_text()
        )
        gap_payload = json.loads(
            (output_dir / "release_readiness_gap_table.json").read_text()
        )

        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(len(summary_payload), 1)
        self.assertEqual(len(gap_rows), 12)
        self.assertEqual(len(gap_payload), 12)
        self.assertEqual(len(reference_rows), 7)
        summary = summary_rows[0]
        self.assertEqual(
            summary["release_readiness_decision"],
            "not_public_alpha_ready_source_freeze_card_and_release_target_blockers",
        )
        self.assertEqual(summary["source_checksum_status"], "checksum_manifest_parsed")
        self.assertEqual(summary["source_freeze_ready"], "False")
        gap_by_id = {row["gap_id"]: row for row in gap_rows}
        self.assertEqual(
            gap_by_id["public_alpha_data_card_and_citation"]["readiness_status"],
            "blocker",
        )
        self.assertIn(
            "https://huggingface.co/docs/hub/datasets-cards",
            gap_by_id["public_alpha_data_card_and_citation"][
                "external_reference_urls"
            ],
        )
        self.assertEqual(
            gap_by_id["diagnostic_not_leaderboard"]["readiness_status"],
            "acceptable_draft_limitation",
        )

    def test_v9_osd120_payload_freeze_manifest_api(self):
        from spacebio_bench import (
            write_multispecies_interaction_payload_freeze_manifest,
        )

        checksum_manifest_text = "\n".join(
            [
                "7e4a87241e4272909c8d9a47b3d45e1f SampleTable_GLbulkRNAseq.csv",
                "d2e9e835f60ef0752a47ecdcbca1b9af Normalized_Counts_GLbulkRNAseq.csv",
                "0123456789abcdef0123456789abcdef out_of_scope.txt",
            ]
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_payload_freeze_manifest"
            outputs = write_multispecies_interaction_payload_freeze_manifest(
                repo_root=REPO_ROOT,
                output_dir=output_dir,
                checksum_manifest_text=checksum_manifest_text,
            )
            with outputs["summary_csv"].open(newline="") as handle:
                summary_rows = list(csv.DictReader(handle))
            with outputs["manifest_csv"].open(newline="") as handle:
                payload_rows = list(csv.DictReader(handle))

        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(len(payload_rows), 3)
        summary = summary_rows[0]
        self.assertEqual(summary["required_payload_count"], "2")
        self.assertEqual(summary["required_payload_matched_count"], "2")
        self.assertEqual(summary["out_of_scope_payload_count"], "1")
        self.assertEqual(summary["diagnostic_required_payload_freeze_ready"], "True")
        payload_by_name = {row["payload_filename"]: row for row in payload_rows}
        self.assertEqual(
            payload_by_name["SampleTable_GLbulkRNAseq.csv"]["payload_role"],
            "sample_table",
        )
        self.assertEqual(
            payload_by_name["Normalized_Counts_GLbulkRNAseq.csv"][
                "verification_status"
            ],
            "required_payload_md5_matched",
        )

    def test_build_v9_osd120_payload_freeze_manifest_cli(self):
        checksum_manifest_text = "\n".join(
            [
                "7e4a87241e4272909c8d9a47b3d45e1f SampleTable_GLbulkRNAseq.csv",
                "d2e9e835f60ef0752a47ecdcbca1b9af Normalized_Counts_GLbulkRNAseq.csv",
                "0123456789abcdef0123456789abcdef out_of_scope.txt",
            ]
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_payload_freeze_manifest"
            fixture = Path(tmpdir) / "osd120_md5sum.tsv"
            fixture.write_text(checksum_manifest_text + "\n")
            result = subprocess.run(
                [
                    sys.executable,
                    str(
                        REPO_ROOT
                        / "scripts"
                        / "build_v9_osd120_payload_freeze_manifest.py"
                    ),
                    "--repo-root",
                    str(REPO_ROOT),
                    "--checksum-manifest-text-path",
                    str(fixture),
                    "--output-dir",
                    str(output_dir),
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            with (output_dir / "payload_freeze_summary.csv").open(
                newline=""
            ) as handle:
                summary_rows = list(csv.DictReader(handle))
            with (output_dir / "payload_freeze_manifest.csv").open(
                newline=""
            ) as handle:
                payload_rows = list(csv.DictReader(handle))

        self.assertIn("payload_freeze_manifest.json", result.stdout)
        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(len(payload_rows), 3)
        self.assertEqual(summary_rows[0]["required_payload_missing_count"], "0")

    def test_generated_v9_osd120_payload_freeze_manifest_outputs_validate(self):
        output_dir = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "interaction_payload_freeze_manifest"
        )
        with (output_dir / "payload_freeze_summary.csv").open(newline="") as handle:
            summary_rows = list(csv.DictReader(handle))
        with (output_dir / "payload_freeze_manifest.csv").open(newline="") as handle:
            payload_rows = list(csv.DictReader(handle))
        summary_payload = json.loads(
            (output_dir / "payload_freeze_summary.json").read_text()
        )
        manifest_payload = json.loads(
            (output_dir / "payload_freeze_manifest.json").read_text()
        )

        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(len(summary_payload), 1)
        self.assertEqual(len(payload_rows), 533)
        self.assertEqual(len(manifest_payload), 533)
        summary = summary_rows[0]
        self.assertEqual(summary["required_payload_count"], "2")
        self.assertEqual(summary["required_payload_matched_count"], "2")
        self.assertEqual(summary["out_of_scope_payload_count"], "531")
        self.assertEqual(
            summary["release_scope_decision"],
            "diagnostic_required_payloads_frozen_full_osdr_payloads_not_frozen",
        )
        payload_by_name = {row["payload_filename"]: row for row in payload_rows}
        self.assertEqual(
            payload_by_name["SampleTable_GLbulkRNAseq.csv"]["checksum_match"],
            "True",
        )
        self.assertEqual(
            payload_by_name["Normalized_Counts_GLbulkRNAseq.csv"][
                "observed_checksum"
            ],
            "d2e9e835f60ef0752a47ecdcbca1b9af",
        )
        self.assertEqual(
            payload_by_name["Normalized_Counts_rRNArm_GLbulkRNAseq.csv"][
                "verification_status"
            ],
            "out_of_scope_not_downloaded",
        )

    def test_v9_osd120_public_alpha_card_api(self):
        from spacebio_bench import write_multispecies_interaction_public_alpha_card

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_public_alpha_card"
            outputs = write_multispecies_interaction_public_alpha_card(
                repo_root=REPO_ROOT,
                output_dir=output_dir,
            )
            with outputs["summary_csv"].open(newline="") as handle:
                summary_rows = list(csv.DictReader(handle))
            card = outputs["card_md"].read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(summary["card_status"], "draft_public_alpha_card_not_frozen_release")
        self.assertEqual(summary["diagnostic_required_payload_freeze_ready"], "True")
        self.assertEqual(summary["full_osdr_payload_freeze_ready"], "False")
        self.assertIn("Payload Freeze Boundary", card)
        self.assertIn("Do not claim a complete local OSD-120 OSDR payload mirror", card)

    def test_build_v9_osd120_public_alpha_card_cli(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_public_alpha_card"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "build_v9_osd120_public_alpha_card.py"),
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
            with (output_dir / "public_alpha_card_summary.csv").open(
                newline=""
            ) as handle:
                summary_rows = list(csv.DictReader(handle))
            card = (output_dir / "public_alpha_card.md").read_text()

        self.assertIn("public_alpha_card.md", result.stdout)
        self.assertEqual(len(summary_rows), 1)
        self.assertIn("Release status: draft diagnostic alpha card", card)

    def test_generated_v9_osd120_public_alpha_card_outputs_validate(self):
        output_dir = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "interaction_public_alpha_card"
        )
        with (output_dir / "public_alpha_card_summary.csv").open(
            newline=""
        ) as handle:
            summary_rows = list(csv.DictReader(handle))
        summary_payload = json.loads(
            (output_dir / "public_alpha_card_summary.json").read_text()
        )
        card = (output_dir / "public_alpha_card.md").read_text()

        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(len(summary_payload), 1)
        summary = summary_rows[0]
        self.assertEqual(summary["allowed_claim_count"], "7")
        self.assertEqual(summary["disallowed_claim_count"], "5")
        self.assertEqual(summary["artifact_count"], "26")
        self.assertEqual(summary["claim_count"], "7")
        self.assertIn("Condition stratum: Col.0.PhyD \\| Dark", card)
        self.assertIn("Full OSDR processed payload mirror is not claimed.", card)
        self.assertIn("https://osdr.nasa.gov/bio/repo/data/studies/OSD-120", card)

    def test_v9_osd120_rebuild_gate_manifest_api(self):
        from spacebio_bench import write_multispecies_interaction_rebuild_gate_manifest

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_rebuild_gate"
            outputs = write_multispecies_interaction_rebuild_gate_manifest(
                repo_root=REPO_ROOT,
                output_dir=output_dir,
            )
            with outputs["summary_csv"].open(newline="") as handle:
                summary_rows = list(csv.DictReader(handle))
            with outputs["step_csv"].open(newline="") as handle:
                step_rows = list(csv.DictReader(handle))
            environment_rows = json.loads(outputs["environment_json"].read_text())
            review = outputs["review_md"].read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(summary["gate_status"], "ready_existing_outputs_present")
        self.assertEqual(summary["step_count"], "8")
        self.assertEqual(summary["missing_output_count"], "0")
        self.assertEqual(summary["script_missing_count"], "0")
        self.assertEqual(len(step_rows), 8)
        self.assertEqual(
            {row["status"] for row in step_rows},
            {"ready_existing_outputs_present"},
        )
        self.assertEqual(
            step_rows[6]["execution_policy"],
            "network_optional_use_checksum_manifest_fixture_to_skip_fetch",
        )
        self.assertIn("numpy_version", {row["key"] for row in environment_rows})
        self.assertIn("--mode preflight", review)

    def test_rebuild_v9_osd120_diagnostic_package_cli(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_rebuild_gate"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "rebuild_v9_osd120_diagnostic_package.py"),
                    "--repo-root",
                    str(REPO_ROOT),
                    "--output-dir",
                    str(output_dir),
                    "--mode",
                    "preflight",
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            with (output_dir / "rebuild_gate_summary.csv").open(
                newline=""
            ) as handle:
                summary_rows = list(csv.DictReader(handle))
            with (output_dir / "rebuild_gate_steps.csv").open(newline="") as handle:
                step_rows = list(csv.DictReader(handle))

        self.assertIn("rebuild_gate_environment.json", result.stdout)
        self.assertIn("OSD120_INTERACTION_REBUILD_GATE_REVIEW.md", result.stdout)
        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(summary_rows[0]["gate_status"], "ready_existing_outputs_present")
        self.assertEqual(len(step_rows), 8)

    def test_generated_v9_osd120_rebuild_gate_outputs_validate(self):
        output_dir = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "interaction_rebuild_gate"
        )
        with (output_dir / "rebuild_gate_summary.csv").open(newline="") as handle:
            summary_rows = list(csv.DictReader(handle))
        with (output_dir / "rebuild_gate_steps.csv").open(newline="") as handle:
            step_rows = list(csv.DictReader(handle))
        environment_payload = json.loads(
            (output_dir / "rebuild_gate_environment.json").read_text()
        )
        review = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "OSD120_INTERACTION_REBUILD_GATE_REVIEW.md"
        ).read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(summary["gate_status"], "ready_existing_outputs_present")
        self.assertEqual(summary["ready_step_count"], "8")
        self.assertEqual(summary["hashed_output_count"], "40")
        self.assertEqual(summary["claim_boundary"], "diagnostic_packaging_preflight_only_not_model_retraining_or_frozen_release")
        self.assertEqual(len(step_rows), 8)
        self.assertTrue(all(row["missing_output_count"] == "0" for row in step_rows))
        self.assertIn("scikit_learn_version", {row["key"] for row in environment_payload})
        self.assertIn("python3 scripts/rebuild_v9_osd120_diagnostic_package.py", review)
        self.assertIn("not rerun sparse L1 model grids", review)

    def test_v9_osd120_public_metadata_package_api(self):
        from spacebio_bench import write_multispecies_interaction_public_metadata_package

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_public_metadata_package"
            outputs = write_multispecies_interaction_public_metadata_package(
                repo_root=REPO_ROOT,
                output_dir=output_dir,
            )
            with outputs["summary_csv"].open(newline="") as handle:
                summary_rows = list(csv.DictReader(handle))
            with outputs["target_csv"].open(newline="") as handle:
                target_rows = list(csv.DictReader(handle))
            with outputs["field_csv"].open(newline="") as handle:
                field_rows = list(csv.DictReader(handle))
            skeleton = json.loads(outputs["metadata_skeleton_json"].read_text())
            review = outputs["review_md"].read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(
            summary["metadata_package_status"],
            "draft_metadata_skeleton_ready_not_release_frozen",
        )
        self.assertEqual(summary["public_now_target_count"], "1")
        self.assertEqual(summary["not_public_now_target_count"], "3")
        self.assertEqual(summary["metadata_field_count"], "20")
        self.assertEqual(summary["metadata_present_field_count"], "12")
        self.assertEqual(summary["metadata_partial_field_count"], "3")
        self.assertEqual(summary["metadata_placeholder_field_count"], "5")
        self.assertEqual(len(target_rows), 4)
        self.assertEqual(target_rows[0]["target_id"], "diagnostic_alpha_metadata_draft")
        self.assertEqual(target_rows[0]["public_now"], "True")
        self.assertEqual(
            {row["status"] for row in field_rows},
            {"present", "partial", "placeholder"},
        )
        self.assertEqual(skeleton["standards_profile"]["datacite_schema_version"], "4.7")
        self.assertEqual(skeleton["standards_profile"]["ro_crate_version"], "1.2")
        self.assertIn("DataCite 4.7", review)

    def test_build_v9_osd120_public_metadata_package_cli(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_public_metadata_package"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "build_v9_osd120_public_metadata_package.py"),
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
            with (output_dir / "public_metadata_summary.csv").open(
                newline=""
            ) as handle:
                summary_rows = list(csv.DictReader(handle))
            skeleton = json.loads(
                (output_dir / "public_metadata_skeleton.json").read_text()
            )

        self.assertIn("public_metadata_skeleton.json", result.stdout)
        self.assertIn("OSD120_INTERACTION_PUBLIC_METADATA_PACKAGE_REVIEW.md", result.stdout)
        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(summary_rows[0]["datacite_schema_version"], "4.7")
        self.assertEqual(skeleton["release_targets"][0]["public_now"], "True")

    def test_generated_v9_osd120_public_metadata_package_outputs_validate(self):
        output_dir = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "interaction_public_metadata_package"
        )
        with (output_dir / "public_metadata_summary.csv").open(newline="") as handle:
            summary_rows = list(csv.DictReader(handle))
        with (output_dir / "source_release_target_decision.csv").open(
            newline=""
        ) as handle:
            target_rows = list(csv.DictReader(handle))
        field_payload = json.loads(
            (output_dir / "public_metadata_field_table.json").read_text()
        )
        with (output_dir / "public_metadata_external_references.csv").open(
            newline=""
        ) as handle:
            reference_rows = list(csv.DictReader(handle))
        skeleton = json.loads(
            (output_dir / "public_metadata_skeleton.json").read_text()
        )
        review = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "OSD120_INTERACTION_PUBLIC_METADATA_PACKAGE_REVIEW.md"
        ).read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(summary["release_target_decision"], "diagnostic_metadata_public_now_source_release_target_pending")
        self.assertEqual(summary["full_osdr_payload_freeze_ready"], "False")
        self.assertEqual(summary["rebuild_gate_status"], "ready_existing_outputs_present")
        self.assertEqual(summary["claim_boundary"], "diagnostic_public_metadata_only_not_frozen_benchmark_release")
        self.assertEqual(len(target_rows), 4)
        self.assertEqual(
            target_rows[2]["target_status"],
            "blocked_full_payload_freeze_missing",
        )
        self.assertEqual(len(field_payload), 20)
        self.assertIn("datacite_schema_4_7", {row["reference_id"] for row in reference_rows})
        self.assertEqual(
            skeleton["ro_crate_skeleton"]["@context"],
            "https://w3id.org/ro/crate/1.2/context",
        )
        self.assertEqual(len(skeleton["ro_crate_skeleton"]["@graph"]), 30)
        self.assertIn("does not claim a", review)

    def test_v9_osd120_ro_crate_citation_scaffold_api(self):
        from spacebio_bench import write_multispecies_interaction_ro_crate_citation_scaffold

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_ro_crate_citation_scaffold"
            outputs = write_multispecies_interaction_ro_crate_citation_scaffold(
                repo_root=REPO_ROOT,
                output_dir=output_dir,
            )
            with outputs["summary_csv"].open(newline="") as handle:
                summary_rows = list(csv.DictReader(handle))
            with outputs["validation_csv"].open(newline="") as handle:
                validation_rows = list(csv.DictReader(handle))
            with outputs["citation_csv"].open(newline="") as handle:
                citation_rows = list(csv.DictReader(handle))
            ro_crate = json.loads(outputs["ro_crate_json"].read_text())
            data_package = json.loads(outputs["data_package_json"].read_text())
            review = outputs["review_md"].read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(
            summary["scaffold_status"],
            "draft_scaffold_ready_archive_blocked_by_citation_placeholders",
        )
        self.assertEqual(summary["ro_crate_graph_entity_count"], "30")
        self.assertEqual(summary["ro_crate_data_entity_count"], "26")
        self.assertEqual(summary["data_package_resource_count"], "26")
        self.assertEqual(summary["validation_blocker_count"], "3")
        self.assertEqual(summary["citation_blocker_count"], "4")
        self.assertEqual(len(validation_rows), 13)
        self.assertEqual(len(citation_rows), 11)
        self.assertEqual(ro_crate["@context"], "https://w3id.org/ro/crate/1.2/context")
        root = next(row for row in ro_crate["@graph"] if row["@id"] == "./")
        self.assertEqual(root["identifier"], "pending_versioned_release_identifier")
        self.assertEqual(data_package["profile"], "data-package")
        self.assertEqual(data_package["licenses"][0]["name"], "pending-local-release-license-review")
        self.assertIn("not archive-ready", review)

    def test_build_v9_osd120_ro_crate_citation_scaffold_cli(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_ro_crate_citation_scaffold"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "build_v9_osd120_ro_crate_citation_scaffold.py"),
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
            with (output_dir / "ro_crate_export_summary.csv").open(
                newline=""
            ) as handle:
                summary_rows = list(csv.DictReader(handle))
            data_package = json.loads((output_dir / "datapackage.draft.json").read_text())

        self.assertIn("ro-crate-metadata.draft.json", result.stdout)
        self.assertIn("datapackage.draft.json", result.stdout)
        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(summary_rows[0]["citation_needs_review_count"], "2")
        self.assertEqual(len(data_package["resources"]), 26)

    def test_generated_v9_osd120_ro_crate_citation_scaffold_outputs_validate(self):
        output_dir = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "interaction_ro_crate_citation_scaffold"
        )
        with (output_dir / "ro_crate_export_summary.csv").open(newline="") as handle:
            summary_rows = list(csv.DictReader(handle))
        with (output_dir / "ro_crate_validation_table.csv").open(
            newline=""
        ) as handle:
            validation_rows = list(csv.DictReader(handle))
        with (output_dir / "citation_freeze_checklist.csv").open(
            newline=""
        ) as handle:
            citation_rows = list(csv.DictReader(handle))
        ro_crate = json.loads((output_dir / "ro-crate-metadata.draft.json").read_text())
        data_package = json.loads((output_dir / "datapackage.draft.json").read_text())
        review = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "OSD120_INTERACTION_RO_CRATE_CITATION_SCAFFOLD_REVIEW.md"
        ).read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(summary["placeholder_field_count"], "5")
        self.assertEqual(summary["next_required_block"], "V9-MULTI-032: archive identifier and license decision gate")
        self.assertEqual(summary["claim_boundary"], "diagnostic_ro_crate_scaffold_only_not_citable_archive")
        self.assertEqual(
            {row["check_status"] for row in validation_rows},
            {"pass", "needs_review", "blocker"},
        )
        self.assertIn("datacite_identifier_placeholder", {row["check_id"] for row in validation_rows})
        self.assertIn("license_and_rights", {row["item_id"] for row in citation_rows})
        self.assertEqual(len(ro_crate["@graph"]), 30)
        self.assertEqual(len(data_package["resources"]), 26)
        self.assertEqual(data_package["custom"]["claim_boundary"], "diagnostic_ro_crate_scaffold_only_not_citable_archive")
        self.assertIn("V9-MULTI-032", review)

    def test_v9_osd120_archive_decision_gate_api(self):
        from spacebio_bench import write_multispecies_interaction_archive_decision_gate

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_archive_decision_gate"
            outputs = write_multispecies_interaction_archive_decision_gate(
                repo_root=REPO_ROOT,
                output_dir=output_dir,
            )
            with outputs["summary_csv"].open(newline="") as handle:
                summary_rows = list(csv.DictReader(handle))
            options = json.loads(outputs["archive_option_json"].read_text())
            license_rows = json.loads(outputs["license_json"].read_text())
            creator_rows = json.loads(outputs["creator_json"].read_text())
            reference_rows = json.loads(outputs["reference_json"].read_text())
            review = outputs["review_md"].read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(
            summary["decision_status"],
            "current_draft_no_archive_selected_release_decisions_blocked",
        )
        self.assertEqual(summary["archive_option_count"], "6")
        self.assertEqual(summary["current_selected_option_count"], "2")
        self.assertEqual(summary["deferred_option_count"], "3")
        self.assertEqual(summary["blocked_option_count"], "1")
        self.assertEqual(summary["license_blocker_count"], "3")
        self.assertEqual(summary["creator_blocker_count"], "4")
        self.assertEqual(summary["external_reference_count"], "6")

        option_by_id = {row["option_id"]: row for row in options}
        self.assertEqual(
            option_by_id["current_no_archive_diagnostic_draft"]["current_draft_selected"],
            "True",
        )
        self.assertEqual(
            option_by_id["zenodo_github_release_doi"]["decision_status"],
            "deferred_release_owner_required",
        )
        license_by_id = {row["component_id"]: row for row in license_rows}
        creator_by_id = {row["component_id"]: row for row in creator_rows}
        self.assertEqual(license_by_id["local_code_license"]["decision_status"], "blocked")
        self.assertEqual(creator_by_id["local_package_creators"]["decision_status"], "blocked")
        self.assertIn("github_zenodo_doi", {row["reference_id"] for row in reference_rows})
        self.assertIn("does not mint a DOI", review)

    def test_build_v9_osd120_archive_decision_gate_cli(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_archive_decision_gate"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "build_v9_osd120_archive_decision_gate.py"),
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
            with (output_dir / "archive_decision_summary.csv").open(
                newline=""
            ) as handle:
                summary_rows = list(csv.DictReader(handle))
            with (output_dir / "archive_decision_external_references.csv").open(
                newline=""
            ) as handle:
                reference_rows = list(csv.DictReader(handle))

        self.assertIn("archive_decision_summary.csv", result.stdout)
        self.assertIn("OSD120_INTERACTION_ARCHIVE_IDENTIFIER_LICENSE_DECISION_REVIEW.md", result.stdout)
        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(
            summary_rows[0]["archive_path_decision"],
            "no_archive_identifier_for_current_draft",
        )
        self.assertEqual(len(reference_rows), 6)

    def test_generated_v9_osd120_archive_decision_gate_outputs_validate(self):
        output_dir = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "interaction_archive_decision_gate"
        )
        with (output_dir / "archive_decision_summary.csv").open(newline="") as handle:
            summary_rows = list(csv.DictReader(handle))
        with (output_dir / "archive_identifier_option_matrix.csv").open(
            newline=""
        ) as handle:
            option_rows = list(csv.DictReader(handle))
        with (output_dir / "license_rights_decision_table.csv").open(
            newline=""
        ) as handle:
            license_rows = list(csv.DictReader(handle))
        with (output_dir / "creator_contributor_decision_table.csv").open(
            newline=""
        ) as handle:
            creator_rows = list(csv.DictReader(handle))
        with (output_dir / "archive_decision_external_references.csv").open(
            newline=""
        ) as handle:
            reference_rows = list(csv.DictReader(handle))
        review = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "OSD120_INTERACTION_ARCHIVE_IDENTIFIER_LICENSE_DECISION_REVIEW.md"
        ).read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(
            summary["next_required_block"],
            "V9-MULTI-033: release owner supplied citation metadata fill",
        )
        self.assertEqual(
            summary["claim_boundary"],
            "diagnostic_archive_decision_gate_only_no_archive_identifier",
        )
        self.assertEqual(summary["archive_option_count"], "6")
        self.assertEqual(summary["license_blocker_count"], "3")
        self.assertEqual(summary["creator_blocker_count"], "4")
        self.assertEqual(
            {row["decision_status"] for row in option_rows},
            {
                "selected_for_current_draft_only",
                "deferred_release_owner_required",
                "deferred_after_creator_version_decision",
                "selected_for_upstream_credit_not_local_archive",
                "deferred_code_archive_related_identifier_only",
                "blocked_full_payload_freeze_missing",
            },
        )
        self.assertEqual(
            sum(1 for row in license_rows if row["decision_status"] == "blocked"),
            3,
        )
        self.assertEqual(
            sum(1 for row in creator_rows if row["decision_status"] == "blocked"),
            4,
        )
        self.assertIn("spdx_license_identifiers", {row["reference_id"] for row in reference_rows})
        self.assertIn("No DOI", review)
        self.assertIn("V9-MULTI-033", review)

    def test_v9_osd120_citation_metadata_fill_api_no_owner_metadata(self):
        from spacebio_bench import write_multispecies_interaction_citation_metadata_fill

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_citation_metadata_fill"
            outputs = write_multispecies_interaction_citation_metadata_fill(
                repo_root=REPO_ROOT,
                output_dir=output_dir,
            )
            with outputs["summary_csv"].open(newline="") as handle:
                summary_rows = list(csv.DictReader(handle))
            fill_rows = json.loads(outputs["fill_status_json"].read_text())
            descriptor_preview = json.loads(outputs["descriptor_preview_json"].read_text())
            review = outputs["review_md"].read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(
            summary["fill_status"],
            "no_owner_metadata_supplied_no_descriptor_changes",
        )
        self.assertEqual(summary["intake_field_count"], "16")
        self.assertEqual(summary["supplied_field_count"], "0")
        self.assertEqual(summary["retained_current_draft_count"], "4")
        self.assertEqual(summary["not_supplied_blocker_count"], "12")
        self.assertEqual(summary["needs_review_count"], "2")
        self.assertEqual(summary["ro_crate_mutation_status"], "not_mutated_preview_only")
        self.assertEqual(summary["datapackage_mutation_status"], "not_mutated_preview_only")
        fill_by_id = {row["field_id"]: row for row in fill_rows}
        self.assertEqual(
            fill_by_id["future_archive_identifier"]["fill_status"],
            "not_supplied_blocking_release",
        )
        self.assertEqual(
            fill_by_id["package_title"]["fill_status"],
            "existing_current_draft_value_retained",
        )
        self.assertEqual(
            fill_by_id["local_code_license"]["validation_rule"],
            "spdx_identifier_or_explicit_custom_terms_required",
        )
        self.assertFalse(descriptor_preview["mutates_ro_crate"])
        self.assertFalse(descriptor_preview["mutates_datapackage"])
        self.assertEqual(descriptor_preview["owner_supplied_patch_preview"], [])
        self.assertIn("does not mint a DOI", review)

    def test_build_v9_osd120_citation_metadata_fill_cli_with_owner_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            owner_metadata = tmp / "owner_metadata.csv"
            owner_metadata.write_text(
                "\n".join(
                    [
                        "field_id,supplied_value,supplied_by,supplied_date,supplied_evidence",
                        "package_title,Owner Reviewed OSD-120 Draft,release-owner,2026-05-26,email-record",
                        "local_code_license,MIT,release-owner,2026-05-26,LICENSE-review",
                    ]
                )
                + "\n"
            )
            output_dir = tmp / "interaction_citation_metadata_fill"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "build_v9_osd120_citation_metadata_fill.py"),
                    "--repo-root",
                    str(REPO_ROOT),
                    "--output-dir",
                    str(output_dir),
                    "--owner-metadata",
                    str(owner_metadata),
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            with (output_dir / "citation_metadata_fill_summary.csv").open(
                newline=""
            ) as handle:
                summary_rows = list(csv.DictReader(handle))
            with (output_dir / "citation_metadata_fill_status.csv").open(
                newline=""
            ) as handle:
                fill_rows = list(csv.DictReader(handle))
            descriptor_preview = json.loads(
                (output_dir / "citation_descriptor_patch_preview.json").read_text()
            )

        self.assertIn("citation_metadata_fill_summary.csv", result.stdout)
        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(
            summary["fill_status"],
            "owner_metadata_supplied_patch_preview_only_not_applied",
        )
        self.assertEqual(summary["supplied_field_count"], "2")
        self.assertEqual(summary["descriptor_patch_status"], "patch_preview_available_not_applied")
        fill_by_id = {row["field_id"]: row for row in fill_rows}
        self.assertEqual(
            fill_by_id["local_code_license"]["fill_status"],
            "owner_supplied_pending_license_review",
        )
        self.assertEqual(
            fill_by_id["package_title"]["fill_status"],
            "owner_supplied_pending_release_review",
        )
        self.assertEqual(len(descriptor_preview["owner_supplied_patch_preview"]), 2)
        self.assertFalse(descriptor_preview["release_ready_after_fill"])

    def test_generated_v9_osd120_citation_metadata_fill_outputs_validate(self):
        output_dir = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "interaction_citation_metadata_fill"
        )
        with (output_dir / "citation_metadata_fill_summary.csv").open(
            newline=""
        ) as handle:
            summary_rows = list(csv.DictReader(handle))
        with (output_dir / "citation_metadata_fill_status.csv").open(
            newline=""
        ) as handle:
            fill_rows = list(csv.DictReader(handle))
        descriptor_preview = json.loads(
            (output_dir / "citation_descriptor_patch_preview.json").read_text()
        )
        review = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "OSD120_INTERACTION_CITATION_METADATA_FILL_REVIEW.md"
        ).read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(
            summary["claim_boundary"],
            "owner_metadata_intake_only_no_descriptor_mutation",
        )
        self.assertEqual(summary["release_ready_after_fill"], "False")
        self.assertEqual(summary["not_supplied_blocker_count"], "12")
        self.assertEqual({row["fill_status"] for row in fill_rows}, {
            "existing_current_draft_value_retained",
            "not_supplied_blocking_release",
        })
        self.assertIn("local_package_creators", descriptor_preview["blocked_fields"])
        self.assertEqual(descriptor_preview["safe_current_draft_values"]["archive_path_decision"], "no_archive_identifier_for_current_draft")
        self.assertIn("No archive", review)
        self.assertIn("no descriptor mutation", review)

    def test_v9_osd120_archive_release_deferral_guard_api(self):
        from spacebio_bench import (
            write_multispecies_interaction_archive_release_deferral_guard,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_archive_release_deferral_guard"
            outputs = write_multispecies_interaction_archive_release_deferral_guard(
                repo_root=REPO_ROOT,
                output_dir=output_dir,
            )
            with outputs["summary_csv"].open(newline="") as handle:
                summary_rows = list(csv.DictReader(handle))
            guard_rows = json.loads(outputs["guard_json"].read_text())
            action_rows = json.loads(outputs["action_json"].read_text())
            mutation_guard = json.loads(outputs["mutation_guard_json"].read_text())
            review = outputs["review_md"].read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(summary["guard_status"], "archive_release_deferred_no_owner_metadata")
        self.assertEqual(summary["release_deferral_status"], "deferred_keep_diagnostic_metadata_only")
        self.assertEqual(summary["blocked_field_count"], "12")
        self.assertEqual(summary["application_guard_count"], "11")
        self.assertEqual(summary["guard_pass_count"], "2")
        self.assertEqual(summary["guard_blocker_count"], "9")
        self.assertEqual(summary["action_count"], "9")
        self.assertEqual(summary["required_owner_action_count"], "6")
        self.assertEqual(summary["descriptor_mutation_allowed"], "False")
        guard_by_id = {row["guard_check_id"]: row for row in guard_rows}
        self.assertEqual(guard_by_id["descriptor_mutation_prevented"]["guard_status"], "pass")
        self.assertEqual(guard_by_id["owner_metadata_file_present"]["guard_status"], "blocker")
        self.assertIn("future_archive_identifier", mutation_guard["blocked_fields"])
        self.assertFalse(mutation_guard["descriptor_mutation_allowed"])
        self.assertIn("defer_archive_release", {row["action_id"] for row in action_rows})
        self.assertIn("defers archive release", review)

    def test_build_v9_osd120_archive_release_deferral_guard_cli(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_archive_release_deferral_guard"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "build_v9_osd120_archive_release_deferral_guard.py"),
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
            with (output_dir / "archive_release_deferral_summary.csv").open(
                newline=""
            ) as handle:
                summary_rows = list(csv.DictReader(handle))
            with (output_dir / "owner_metadata_application_guard.csv").open(
                newline=""
            ) as handle:
                guard_rows = list(csv.DictReader(handle))

        self.assertIn("archive_release_deferral_summary.csv", result.stdout)
        self.assertIn("descriptor_mutation_guard.json", result.stdout)
        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(summary_rows[0]["release_ready_after_guard"], "False")
        self.assertEqual(len(guard_rows), 11)

    def test_generated_v9_osd120_archive_release_deferral_guard_outputs_validate(self):
        output_dir = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "interaction_archive_release_deferral_guard"
        )
        with (output_dir / "archive_release_deferral_summary.csv").open(
            newline=""
        ) as handle:
            summary_rows = list(csv.DictReader(handle))
        with (output_dir / "owner_metadata_application_guard.csv").open(
            newline=""
        ) as handle:
            guard_rows = list(csv.DictReader(handle))
        with (output_dir / "archive_release_deferral_actions.csv").open(
            newline=""
        ) as handle:
            action_rows = list(csv.DictReader(handle))
        mutation_guard = json.loads((output_dir / "descriptor_mutation_guard.json").read_text())
        review = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "OSD120_INTERACTION_ARCHIVE_RELEASE_DEFERRAL_GUARD_REVIEW.md"
        ).read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(
            summary["claim_boundary"],
            "archive_release_deferred_diagnostic_metadata_only_no_descriptor_mutation",
        )
        self.assertEqual(
            summary["next_required_block"],
            "V9-MULTI-035: diagnostic metadata release note or owner metadata intake retry",
        )
        self.assertEqual(
            {row["guard_status"] for row in guard_rows},
            {"blocker", "pass"},
        )
        self.assertEqual(
            sum(1 for row in action_rows if row["action_status"] == "required_owner_action"),
            6,
        )
        self.assertEqual(mutation_guard["allowed_next_descriptor_action"], "none_defer_archive_release")
        self.assertFalse(mutation_guard["mutates_citation_cff"])
        self.assertIn("Archive release is deferred", review)
        self.assertIn("no descriptor", review)

    def test_v9_osd120_diagnostic_metadata_release_note_api(self):
        from spacebio_bench import (
            write_multispecies_interaction_diagnostic_metadata_release_note,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_diagnostic_metadata_release_note"
            outputs = write_multispecies_interaction_diagnostic_metadata_release_note(
                repo_root=REPO_ROOT,
                output_dir=output_dir,
            )
            with outputs["summary_csv"].open(newline="") as handle:
                summary_rows = list(csv.DictReader(handle))
            sections = json.loads(outputs["section_json"].read_text())
            claims = json.loads(outputs["claim_json"].read_text())
            retry_items = json.loads(outputs["retry_json"].read_text())
            release_note = outputs["release_note_md"].read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(
            summary["note_status"],
            "diagnostic_metadata_note_ready_no_archive_release",
        )
        self.assertEqual(
            summary["branch_closeout_status"],
            "osd120_metadata_branch_closeout_pending_refocus",
        )
        self.assertEqual(summary["current_public_surface"], "diagnostic_metadata_only")
        self.assertEqual(summary["inspectable_now_count"], "3")
        self.assertEqual(summary["not_released_claim_count"], "3")
        self.assertEqual(summary["owner_retry_item_count"], "6")
        self.assertEqual(summary["descriptor_mutation_status"], "not_mutated")
        self.assertEqual(summary["archive_release_status"], "deferred_no_owner_metadata")
        self.assertEqual(
            summary["next_required_block"],
            "V9-REFOCUS-001: post-OSD-120 recenter decision",
        )
        self.assertEqual(
            summary["claim_boundary"],
            "diagnostic_metadata_note_only_not_archive_release",
        )
        self.assertEqual(len(sections), 5)
        self.assertEqual({row["include_in_note"] for row in sections}, {"True"})
        self.assertEqual(
            {row["public_note_status"] for row in claims},
            {"allowed_current_note", "not_allowed_current_note"},
        )
        self.assertEqual(
            {row["retry_item_id"] for row in retry_items},
            {
                "owner_metadata_file",
                "archive_route_identifier",
                "version_date_year",
                "creator_contributor_publisher",
                "license_scope",
                "exact_osdr_study_citation",
            },
        )
        self.assertIn("diagnostic metadata only", release_note)

    def test_build_v9_osd120_diagnostic_metadata_release_note_cli(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "interaction_diagnostic_metadata_release_note"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "build_v9_osd120_diagnostic_metadata_release_note.py"),
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
            with (output_dir / "diagnostic_metadata_release_summary.csv").open(
                newline=""
            ) as handle:
                summary_rows = list(csv.DictReader(handle))
            with (output_dir / "owner_metadata_retry_checklist.csv").open(
                newline=""
            ) as handle:
                retry_rows = list(csv.DictReader(handle))

        self.assertIn("diagnostic_metadata_release_summary.csv", result.stdout)
        self.assertIn("OSD120_DIAGNOSTIC_METADATA_NOTE.md", result.stdout)
        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(len(retry_rows), 6)

    def test_generated_v9_osd120_diagnostic_metadata_release_note_outputs_validate(self):
        output_dir = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "interaction_diagnostic_metadata_release_note"
        )
        with (output_dir / "diagnostic_metadata_release_summary.csv").open(
            newline=""
        ) as handle:
            summary_rows = list(csv.DictReader(handle))
        with (output_dir / "diagnostic_metadata_public_claims.csv").open(
            newline=""
        ) as handle:
            claim_rows = list(csv.DictReader(handle))
        with (output_dir / "owner_metadata_retry_checklist.csv").open(
            newline=""
        ) as handle:
            retry_rows = list(csv.DictReader(handle))
        review = (
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "reports"
            / "OSD120_INTERACTION_DIAGNOSTIC_METADATA_RELEASE_NOTE_REVIEW.md"
        ).read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(
            summary["next_required_block"],
            "V9-REFOCUS-001: post-OSD-120 recenter decision",
        )
        self.assertEqual(summary["inspectable_now_count"], "3")
        self.assertEqual(summary["not_released_claim_count"], "3")
        self.assertEqual(summary["owner_retry_item_count"], "6")
        self.assertEqual(
            summary["claim_boundary"],
            "diagnostic_metadata_note_only_not_archive_release",
        )
        self.assertEqual(
            sum(1 for row in claim_rows if row["public_note_status"] == "allowed_current_note"),
            3,
        )
        self.assertEqual(
            sum(1 for row in claim_rows if row["public_note_status"] == "not_allowed_current_note"),
            3,
        )
        self.assertEqual({row["retry_priority"] for row in retry_rows}, {"P0", "P1", "P2"})
        self.assertIn("Recenter after", review)
        self.assertIn("not another archive-release gate", review)

    def test_v9_recenter_decision_api_selects_public_bulk_alpha(self):
        from spacebio_bench import write_v9_recenter_decision

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "recenter_decision"
            outputs = write_v9_recenter_decision(
                repo_root=REPO_ROOT,
                output_dir=output_dir,
            )
            with outputs["summary_csv"].open(newline="") as handle:
                summary_rows = list(csv.DictReader(handle))
            candidates = json.loads(outputs["candidate_json"].read_text())
            actions = json.loads(outputs["action_json"].read_text())
            assets = json.loads(outputs["asset_json"].read_text())
            review = outputs["review_md"].read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(summary["decision_status"], "selected_public_bulk_alpha_recenter")
        self.assertEqual(summary["selected_next_lane"], "public_bulk_alpha")
        self.assertEqual(
            summary["selected_next_block"],
            "V9-BULK-ALPHA-001: public bulk alpha freeze-gap matrix",
        )
        self.assertEqual(summary["bulk_task_count"], "8")
        self.assertEqual(summary["bulk_fold_count"], "33")
        self.assertEqual(summary["bulk_source_count"], "22")
        self.assertEqual(summary["bulk_freeze_ready_source_count"], "0")
        self.assertEqual(summary["bulk_baseline_row_count"], "24")
        self.assertEqual(summary["single_cell_h5ad_count"], "0")
        self.assertGreater(
            int(summary["public_bulk_readiness_score"]),
            int(summary["single_cell_readiness_score"]),
        )
        decisions = {row["candidate_id"]: row["decision"] for row in candidates}
        self.assertEqual(decisions["public_bulk_alpha"], "selected_next")
        self.assertEqual(
            decisions["single_cell_flagship"],
            "defer_until_after_bulk_alpha_gap_matrix",
        )
        self.assertIn("bulk_alpha_gap_matrix", {row["action_id"] for row in actions})
        self.assertTrue(any("rrrm" in path.lower() for path in assets))
        self.assertIn("sequencing correction", review)

    def test_build_v9_recenter_decision_cli(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "recenter_decision"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "build_v9_recenter_decision.py"),
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
            with (output_dir / "recenter_decision_summary.csv").open(
                newline=""
            ) as handle:
                summary_rows = list(csv.DictReader(handle))
            with (output_dir / "recenter_candidate_matrix.csv").open(
                newline=""
            ) as handle:
                candidate_rows = list(csv.DictReader(handle))

        self.assertIn("recenter_decision_summary.csv", result.stdout)
        self.assertIn("V9_REFOCUS_001_POST_OSD120_RECENTER_DECISION.md", result.stdout)
        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(len(candidate_rows), 2)

    def test_generated_v9_recenter_decision_outputs_validate(self):
        output_dir = REPO_ROOT / "v9" / "reports" / "recenter_decision"
        with (output_dir / "recenter_decision_summary.csv").open(
            newline=""
        ) as handle:
            summary_rows = list(csv.DictReader(handle))
        with (output_dir / "recenter_candidate_matrix.csv").open(
            newline=""
        ) as handle:
            candidate_rows = list(csv.DictReader(handle))
        with (output_dir / "recenter_next_block_actions.csv").open(
            newline=""
        ) as handle:
            action_rows = list(csv.DictReader(handle))
        review = (
            REPO_ROOT
            / "docs"
            / "V9_REFOCUS_001_POST_OSD120_RECENTER_DECISION.md"
        ).read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(summary["selected_next_lane"], "public_bulk_alpha")
        self.assertEqual(summary["bulk_checksum_parsed_source_count"], "22")
        self.assertEqual(summary["bulk_freeze_ready_source_count"], "0")
        self.assertEqual(summary["single_cell_v9_manifest_count"], "0")
        self.assertEqual(
            summary["claim_boundary"],
            "recenter_decision_only_no_new_benchmark_or_release_claim",
        )
        self.assertEqual(
            {row["decision"] for row in candidate_rows},
            {"selected_next", "defer_until_after_bulk_alpha_gap_matrix"},
        )
        self.assertEqual({row["lane_id"] for row in action_rows}, {"public_bulk_alpha"})
        self.assertIn("does not promote OSD-120", review)

    def test_v9_public_bulk_alpha_gap_matrix_api(self):
        from spacebio_bench import write_public_bulk_alpha_gap_matrix

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "public_bulk_alpha_gap_matrix"
            outputs = write_public_bulk_alpha_gap_matrix(
                repo_root=REPO_ROOT,
                output_dir=output_dir,
            )
            with outputs["summary_csv"].open(newline="") as handle:
                summary_rows = list(csv.DictReader(handle))
            gap_rows = json.loads(outputs["gap_json"].read_text())
            payload_rows = json.loads(outputs["payload_json"].read_text())
            claim_rows = json.loads(outputs["claim_json"].read_text())
            update_rows = json.loads(outputs["update_json"].read_text())
            review = outputs["review_md"].read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(
            summary["decision_status"],
            "metadata_alpha_gap_matrix_ready_payload_hash_blocked",
        )
        self.assertEqual(summary["bulk_task_count"], "8")
        self.assertEqual(summary["bulk_fold_count"], "33")
        self.assertEqual(summary["bulk_source_count"], "22")
        self.assertEqual(summary["checksum_parsed_source_count"], "22")
        self.assertEqual(summary["freeze_ready_source_count"], "0")
        self.assertEqual(summary["baseline_row_count"], "24")
        self.assertEqual(summary["pass_count"], "6")
        self.assertEqual(summary["blocker_count"], "2")
        self.assertEqual(summary["needs_update_count"], "2")
        self.assertEqual(
            summary["next_required_block"],
            "V9-BULK-ALPHA-002: metadata-only alpha snapshot decision",
        )
        by_gap = {row["gap_id"]: row for row in gap_rows}
        self.assertEqual(by_gap["payload_hash_verification"]["readiness_status"], "blocker")
        self.assertEqual(by_gap["minimal_alpha_snapshot_decision"]["readiness_status"], "blocker")
        self.assertEqual(len(payload_rows), 22)
        self.assertEqual(
            {row["current_boundary"] for row in payload_rows},
            {"checksum_manifest_evidence_only"},
        )
        self.assertEqual(
            {row["claim_status"] for row in claim_rows},
            {"allowed_current_alpha_gap", "prohibited_current_alpha_gap"},
        )
        self.assertEqual(len(update_rows), 3)
        self.assertIn("0/22", review)

    def test_build_v9_public_bulk_alpha_gap_matrix_cli(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "public_bulk_alpha_gap_matrix"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "build_v9_public_bulk_alpha_gap_matrix.py"),
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
            with (output_dir / "public_bulk_alpha_gap_summary.csv").open(
                newline=""
            ) as handle:
                summary_rows = list(csv.DictReader(handle))
            with (output_dir / "public_bulk_alpha_gap_matrix.csv").open(
                newline=""
            ) as handle:
                gap_rows = list(csv.DictReader(handle))

        self.assertIn("public_bulk_alpha_gap_summary.csv", result.stdout)
        self.assertIn("V9_PUBLIC_BULK_ALPHA_FREEZE_GAP_MATRIX_REVIEW.md", result.stdout)
        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(len(gap_rows), 10)

    def test_generated_v9_public_bulk_alpha_gap_matrix_outputs_validate(self):
        output_dir = REPO_ROOT / "v9" / "reports" / "public_bulk_alpha_gap_matrix"
        with (output_dir / "public_bulk_alpha_gap_summary.csv").open(
            newline=""
        ) as handle:
            summary_rows = list(csv.DictReader(handle))
        with (output_dir / "public_bulk_alpha_gap_matrix.csv").open(
            newline=""
        ) as handle:
            gap_rows = list(csv.DictReader(handle))
        with (output_dir / "payload_hash_boundary.csv").open(newline="") as handle:
            payload_rows = list(csv.DictReader(handle))
        with (output_dir / "public_bulk_alpha_claim_boundary.csv").open(
            newline=""
        ) as handle:
            claim_rows = list(csv.DictReader(handle))
        review = (
            REPO_ROOT
            / "docs"
            / "V9_PUBLIC_BULK_ALPHA_FREEZE_GAP_MATRIX_REVIEW.md"
        ).read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(
            summary["claim_boundary"],
            "public_bulk_alpha_gap_matrix_no_release_approval",
        )
        self.assertEqual(summary["pass_count"], "6")
        self.assertEqual(summary["blocker_count"], "2")
        self.assertEqual(summary["needs_update_count"], "2")
        self.assertEqual({row["readiness_status"] for row in gap_rows}, {
            "pass",
            "blocker",
            "needs_update",
        })
        self.assertEqual(
            {row["release_blocker"] for row in payload_rows},
            {"local_payload_hash_verification_pending"},
        )
        self.assertEqual(
            sum(1 for row in claim_rows if row["claim_status"] == "allowed_current_alpha_gap"),
            3,
        )
        self.assertEqual(
            sum(1 for row in claim_rows if row["claim_status"] == "prohibited_current_alpha_gap"),
            3,
        )
        self.assertIn("not ready for frozen release language", review)

    def test_v9_public_bulk_alpha_snapshot_decision_api(self):
        from spacebio_bench import write_public_bulk_alpha_snapshot_decision

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "public_bulk_alpha_snapshot_decision"
            outputs = write_public_bulk_alpha_snapshot_decision(
                repo_root=REPO_ROOT,
                output_dir=output_dir,
            )
            with outputs["summary_csv"].open(newline="") as handle:
                summary_rows = list(csv.DictReader(handle))
            option_rows = json.loads(outputs["option_json"].read_text())
            claim_rows = json.loads(outputs["claim_json"].read_text())
            language_rows = json.loads(outputs["language_json"].read_text())
            action_rows = json.loads(outputs["action_json"].read_text())
            review = outputs["review_md"].read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(
            summary["decision_status"],
            "metadata_only_alpha_snapshot_allowed_with_payload_blockers",
        )
        self.assertEqual(summary["selected_path"], "metadata_only_alpha_snapshot")
        self.assertEqual(summary["deferred_path"], "payload_mirror_first")
        self.assertEqual(summary["metadata_only_allowed"], "true")
        self.assertEqual(summary["payload_release_allowed"], "false")
        self.assertEqual(summary["bulk_source_count"], "22")
        self.assertEqual(summary["freeze_ready_source_count"], "0")
        self.assertEqual(summary["gap_blocker_count"], "2")
        self.assertEqual(
            summary["next_required_block"],
            "V9-BULK-ALPHA-003: dataset card and Data Package alpha boundary update",
        )
        self.assertEqual(
            summary["claim_boundary"],
            "metadata_only_public_bulk_alpha_no_payload_release",
        )
        self.assertEqual(
            {row["decision"] for row in option_rows},
            {"selected", "deferred", "rejected_for_current_sequence"},
        )
        self.assertEqual(
            {row["claim_status"] for row in claim_rows},
            {"allowed_metadata_alpha", "prohibited_payload_release"},
        )
        self.assertTrue(
            any("metadata-only alpha snapshot" in row["snippet_text"] for row in language_rows)
        )
        self.assertEqual(len(action_rows), 4)
        self.assertIn("not a frozen payload release", review)

    def test_build_v9_public_bulk_alpha_snapshot_decision_cli(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "public_bulk_alpha_snapshot_decision"
            result = subprocess.run(
                [
                    sys.executable,
                    str(
                        REPO_ROOT
                        / "scripts"
                        / "build_v9_public_bulk_alpha_snapshot_decision.py"
                    ),
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
            with (output_dir / "snapshot_decision_summary.csv").open(
                newline=""
            ) as handle:
                summary_rows = list(csv.DictReader(handle))
            with (output_dir / "snapshot_option_matrix.csv").open(
                newline=""
            ) as handle:
                option_rows = list(csv.DictReader(handle))

        self.assertIn("snapshot_decision_summary.csv", result.stdout)
        self.assertIn(
            "V9_PUBLIC_BULK_ALPHA_METADATA_SNAPSHOT_DECISION.md",
            result.stdout,
        )
        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(len(option_rows), 3)

    def test_generated_v9_public_bulk_alpha_snapshot_decision_outputs_validate(self):
        output_dir = REPO_ROOT / "v9" / "reports" / "public_bulk_alpha_snapshot_decision"
        with (output_dir / "snapshot_decision_summary.csv").open(
            newline=""
        ) as handle:
            summary_rows = list(csv.DictReader(handle))
        with (output_dir / "snapshot_option_matrix.csv").open(
            newline=""
        ) as handle:
            option_rows = list(csv.DictReader(handle))
        with (output_dir / "snapshot_claim_boundary.csv").open(
            newline=""
        ) as handle:
            claim_rows = list(csv.DictReader(handle))
        with (output_dir / "snapshot_next_actions.csv").open(
            newline=""
        ) as handle:
            action_rows = list(csv.DictReader(handle))
        review = (
            REPO_ROOT
            / "docs"
            / "V9_PUBLIC_BULK_ALPHA_METADATA_SNAPSHOT_DECISION.md"
        ).read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(summary["selected_path"], "metadata_only_alpha_snapshot")
        self.assertEqual(summary["metadata_only_allowed"], "true")
        self.assertEqual(summary["payload_release_allowed"], "false")
        self.assertEqual(
            {row["status"] for row in option_rows},
            {
                "allowed_with_explicit_blockers",
                "valid_but_not_required_before_metadata_alpha",
                "too_conservative_for_metadata_scaffold",
            },
        )
        self.assertEqual(
            sum(1 for row in claim_rows if row["claim_status"] == "allowed_metadata_alpha"),
            4,
        )
        self.assertEqual(
            sum(1 for row in claim_rows if row["claim_status"] == "prohibited_payload_release"),
            4,
        )
        self.assertEqual(
            {row["action_status"] for row in action_rows},
            {"pending_next_block", "deferred_after_metadata_alpha", "guarded_no_change"},
        )
        self.assertIn("metadata-only alpha snapshot", review)
        self.assertIn("not a frozen payload release", review)

    def test_v9_single_cell_asset_inventory_api(self):
        from spacebio_bench import write_single_cell_asset_inventory

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "sc_spaceflight"
            outputs = write_single_cell_asset_inventory(
                repo_root=REPO_ROOT,
                output_dir=output_dir,
            )
            with outputs["summary_csv"].open(newline="") as handle:
                summary_rows = list(csv.DictReader(handle))
            asset_rows = json.loads(outputs["asset_json"].read_text())
            payload_rows = json.loads(outputs["payload_json"].read_text())
            review = outputs["review_md"].read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(
            summary["decision_status"],
            "legacy_rrrm_assets_indexed_no_local_anndata_payload",
        )
        self.assertEqual(summary["total_asset_count"], "54")
        self.assertEqual(summary["rrrm1_asset_count"], "41")
        self.assertEqual(summary["rrrm2_asset_count"], "13")
        self.assertEqual(summary["documentation_count"], "8")
        self.assertEqual(summary["script_count"], "31")
        self.assertEqual(summary["generated_cache_count"], "6")
        self.assertEqual(summary["local_anndata_payload_count"], "0")
        current_sc_manifest_count = len(
            list((REPO_ROOT / "v9" / "sc_spaceflight" / "task_manifests").glob("*.json"))
        )
        self.assertEqual(summary["v9_sc_manifest_count"], str(current_sc_manifest_count))
        self.assertEqual(summary["metric_profile_status"], "genelab_sc_profile_present")
        self.assertEqual(
            summary["next_required_block"],
            "V9-SC-002: AnnData task manifest draft",
        )
        self.assertEqual(
            summary["claim_boundary"],
            "single_cell_asset_inventory_only_no_v9_sc_task_or_payload_claim",
        )
        paths = {row["asset_path"] for row in asset_rows}
        self.assertIn("v2/docs/RRRM1_SAMPLE_MANIFEST_2026-03-12.csv", paths)
        self.assertIn("v3/evaluation/F5C_rrrm2_loao.json", paths)
        self.assertIn("exclude_from_v9", {row["promotion_status"] for row in asset_rows})
        self.assertEqual(payload_rows, [])
        self.assertIn("does not claim local AnnData payload availability", review)

    def test_build_v9_sc_rrrm_asset_inventory_cli(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "sc_spaceflight"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "build_v9_sc_rrrm_asset_inventory.py"),
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
            with (output_dir / "asset_inventory_summary.csv").open(
                newline=""
            ) as handle:
                summary_rows = list(csv.DictReader(handle))
            with (output_dir / "asset_inventory.csv").open(newline="") as handle:
                asset_rows = list(csv.DictReader(handle))

        self.assertIn("asset_inventory_summary.csv", result.stdout)
        self.assertIn("asset_inventory.md", result.stdout)
        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(summary_rows[0]["total_asset_count"], "54")
        self.assertEqual(len(asset_rows), 54)

    def test_generated_v9_single_cell_asset_inventory_outputs_validate(self):
        output_dir = REPO_ROOT / "v9" / "sc_spaceflight"
        with (output_dir / "asset_inventory_summary.csv").open(
            newline=""
        ) as handle:
            summary_rows = list(csv.DictReader(handle))
        with (output_dir / "asset_inventory.csv").open(newline="") as handle:
            asset_rows = list(csv.DictReader(handle))
        payload_rows = json.loads((output_dir / "local_payload_scan.json").read_text())
        review = (output_dir / "asset_inventory.md").read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(summary["total_asset_count"], "54")
        self.assertEqual(summary["local_h5ad_count"], "0")
        self.assertEqual(summary["local_loom_count"], "0")
        self.assertEqual(summary["local_mtx_count"], "0")
        self.assertEqual(summary["v9_sc_manifest_count"], "0")
        self.assertEqual(
            {row["asset_family"] for row in asset_rows},
            {"RRRM-1", "RRRM-2"},
        )
        self.assertEqual(payload_rows, [])
        self.assertIn("V9-SC-002: AnnData task manifest draft", review)
        self.assertIn("promote legacy RRRM scores as v9 benchmark results", review)

    def test_v9_sc_anndata_manifest_draft_api(self):
        from spacebio_bench import (
            load_task_manifest,
            write_sc_anndata_manifest_draft,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "sc_spaceflight"
            outputs = write_sc_anndata_manifest_draft(
                repo_root=REPO_ROOT,
                output_dir=output_dir,
            )
            with outputs["summary_csv"].open(newline="") as handle:
                summary_rows = list(csv.DictReader(handle))
            manifest = load_task_manifest(outputs["manifest_json"])
            blockers = json.loads(outputs["blockers_json"].read_text())
            review = outputs["review_md"].read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(
            summary["decision_status"],
            "draft_anndata_manifest_contract_ready_payload_blocked",
        )
        self.assertEqual(summary["task_id"], "draft_rrrm1_blood_single_cell_spaceflight")
        self.assertEqual(summary["selected_source_id"], "OSD-918")
        self.assertEqual(summary["selected_tissue"], "blood")
        self.assertEqual(summary["sample_count"], "8")
        self.assertEqual(summary["flight_sample_count"], "4")
        self.assertEqual(summary["ground_sample_count"], "4")
        self.assertEqual(summary["qc_cell_count"], "4395")
        self.assertEqual(summary["local_h5ad_count"], "0")
        self.assertEqual(summary["local_payload_status"], "blocked_no_local_anndata_payload")
        self.assertEqual(
            summary["next_required_block"],
            "V9-SC-003: genelab_sc metric specification",
        )
        self.assertEqual(manifest["task_family"], "sc_spaceflight")
        self.assertEqual(manifest["runnable_status"], "blocked_no_local_anndata_payload")
        self.assertEqual(manifest["metrics"][0]["profile"], "genelab_sc")
        self.assertIn("broad_celltype", manifest["anndata_contract"]["required_obs_fields"])
        fold_ids = [fold["fold_id"] for fold in manifest["split"]["candidate_folds"]]
        self.assertEqual(len(fold_ids), 8)
        self.assertEqual(len(set(fold_ids)), 8)
        self.assertEqual(
            {row["blocker_status"] for row in blockers},
            {"blocking_runnable_task", "blocking_evaluator", "blocking_leaderboard", "claim_guard"},
        )
        self.assertIn("not a runnable v9 single-cell benchmark task", review)

    def test_build_v9_sc_anndata_manifest_draft_cli(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "sc_spaceflight"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "build_v9_sc_anndata_manifest_draft.py"),
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
            with (output_dir / "anndata_manifest_draft_summary.csv").open(
                newline=""
            ) as handle:
                summary_rows = list(csv.DictReader(handle))
            manifest_path = (
                output_dir
                / "task_manifests"
                / "draft_rrrm1_blood_single_cell_spaceflight.json"
            )
            manifest = json.loads(manifest_path.read_text())

        self.assertIn("anndata_manifest_draft_summary.csv", result.stdout)
        self.assertIn("draft_rrrm1_blood_single_cell_spaceflight.json", result.stdout)
        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(manifest["source_records"][0]["source_id"], "OSD-918")
        self.assertEqual(manifest["split"]["label_distribution"], {"Flight": 4, "Ground": 4})

    def test_generated_v9_sc_anndata_manifest_draft_outputs_validate(self):
        from spacebio_bench import load_task_manifest

        output_dir = REPO_ROOT / "v9" / "sc_spaceflight"
        manifest_path = (
            output_dir
            / "task_manifests"
            / "draft_rrrm1_blood_single_cell_spaceflight.json"
        )
        with (output_dir / "anndata_manifest_draft_summary.csv").open(
            newline=""
        ) as handle:
            summary_rows = list(csv.DictReader(handle))
        with (output_dir / "anndata_manifest_blockers.csv").open(
            newline=""
        ) as handle:
            blocker_rows = list(csv.DictReader(handle))
        manifest = load_task_manifest(manifest_path)
        review = (output_dir / "ANNDATA_MANIFEST_DRAFT_REVIEW.md").read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(
            summary["claim_boundary"],
            "anndata_manifest_contract_only_no_local_payload_or_score_claim",
        )
        self.assertEqual(summary["v9_manifest_validation_status"], "passes_minimal_manifest_validator")
        self.assertEqual(manifest["release_status"], "draft_non_runnable_contract")
        self.assertEqual(
            manifest["provenance"]["legacy_score_status"],
            "reference_only_not_v9_benchmark_result",
        )
        self.assertEqual(len(manifest["split"]["candidate_folds"]), 8)
        self.assertEqual(len(blocker_rows), 4)
        self.assertIn("local_anndata_payload_absent", {row["blocker_id"] for row in blocker_rows})
        self.assertIn("V9-SC-003: genelab_sc metric specification", review)

    def test_v9_sc_metric_spec_api(self):
        from spacebio_bench import build_sc_metric_spec, get_metric_profile, write_sc_metric_spec

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "sc_spaceflight"
            outputs = write_sc_metric_spec(
                repo_root=REPO_ROOT,
                output_dir=output_dir,
            )
            with outputs["summary_csv"].open(newline="") as handle:
                summary_rows = list(csv.DictReader(handle))
            metrics = json.loads(outputs["metrics_json"].read_text())
            inputs = json.loads(outputs["inputs_json"].read_text())
            skip_rows = json.loads(outputs["skip_json"].read_text())
            review = outputs["review_md"].read_text()

        package = build_sc_metric_spec(repo_root=REPO_ROOT)
        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(
            summary["decision_status"],
            "genelab_sc_metric_spec_ready_no_evaluator",
        )
        self.assertEqual(summary["profile"], "genelab_sc")
        self.assertEqual(summary["profile_metric_count"], "6")
        self.assertEqual(summary["primary_metric_count"], "2")
        self.assertEqual(summary["diagnostic_metric_count"], "3")
        self.assertEqual(summary["optional_metric_count"], "1")
        self.assertEqual(summary["manifest_task_id"], "draft_rrrm1_blood_single_cell_spaceflight")
        self.assertEqual(summary["local_payload_status"], "blocked_no_local_anndata_payload")
        self.assertEqual(
            summary["next_required_block"],
            "V9-SC-004: AnnData payload staging and obs/var audit plan",
        )
        self.assertEqual(
            {row["metric_id"] for row in metrics},
            set(get_metric_profile("genelab_sc")["metrics"]),
        )
        self.assertEqual(
            {row["metric_role"] for row in metrics},
            {
                "primary_after_payload_freeze",
                "diagnostic_representation",
                "diagnostic_de_recovery",
                "optional_reconstruction",
            },
        )
        self.assertIn("de_reference_table", {row["input_id"] for row in inputs})
        self.assertEqual({row["reported_value"] for row in skip_rows}, {"NA"})
        self.assertEqual(len(package["metrics"]), 6)
        self.assertIn("No evaluator", review)

    def test_build_v9_sc_metric_spec_cli(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "sc_spaceflight"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "build_v9_sc_metric_spec.py"),
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
            with (output_dir / "metric_spec_summary.csv").open(
                newline=""
            ) as handle:
                summary_rows = list(csv.DictReader(handle))
            with (output_dir / "metric_spec_metrics.csv").open(newline="") as handle:
                metric_rows = list(csv.DictReader(handle))

        self.assertIn("metric_spec_summary.csv", result.stdout)
        self.assertIn("V9_SC_METRIC_SPEC.md", result.stdout)
        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(len(metric_rows), 6)
        self.assertEqual(
            {row["metric_id"] for row in metric_rows},
            {
                "de_overlap_at_n",
                "de_direction_match",
                "mission_discrimination",
                "state_label_auroc",
                "state_label_auprc",
                "expression_mae_when_applicable",
            },
        )

    def test_generated_v9_sc_metric_spec_outputs_validate(self):
        output_dir = REPO_ROOT / "v9" / "sc_spaceflight"
        with (output_dir / "metric_spec_summary.csv").open(newline="") as handle:
            summary_rows = list(csv.DictReader(handle))
        with (output_dir / "metric_spec_metrics.csv").open(newline="") as handle:
            metric_rows = list(csv.DictReader(handle))
        with (output_dir / "metric_spec_required_inputs.csv").open(
            newline=""
        ) as handle:
            input_rows = list(csv.DictReader(handle))
        with (output_dir / "metric_spec_skip_policy.csv").open(newline="") as handle:
            skip_rows = list(csv.DictReader(handle))
        review = (REPO_ROOT / "docs" / "V9_SC_METRIC_SPEC.md").read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(
            summary["claim_boundary"],
            "genelab_sc_metric_spec_only_no_evaluator_or_score_claim",
        )
        self.assertEqual(summary["profile_metric_count"], "6")
        self.assertEqual(len(metric_rows), 6)
        self.assertEqual(len(input_rows), 7)
        self.assertEqual(len(skip_rows), 6)
        by_metric = {row["metric_id"]: row for row in metric_rows}
        self.assertEqual(
            by_metric["state_label_auroc"]["metric_role"],
            "primary_after_payload_freeze",
        )
        self.assertEqual(
            by_metric["de_overlap_at_n"]["claim_status"],
            "diagnostic_only_until_de_reference_frozen",
        )
        self.assertIn("flight_probability", review)
        self.assertIn("V9-SC-004: AnnData payload staging and obs/var audit plan", review)

    def test_v9_sc_payload_staging_plan_api(self):
        from spacebio_bench import build_sc_payload_staging_plan, write_sc_payload_staging_plan

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "sc_spaceflight"
            outputs = write_sc_payload_staging_plan(
                repo_root=REPO_ROOT,
                output_dir=output_dir,
            )
            with outputs["summary_csv"].open(newline="") as handle:
                summary_rows = list(csv.DictReader(handle))
            candidates = json.loads(outputs["candidates_json"].read_text())
            requirements = json.loads(outputs["audit_json"].read_text())
            actions = json.loads(outputs["actions_json"].read_text())
            review = outputs["review_md"].read_text()

        package = build_sc_payload_staging_plan(repo_root=REPO_ROOT)
        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(
            summary["decision_status"],
            "payload_staging_plan_ready_no_local_payload",
        )
        self.assertEqual(summary["task_id"], "draft_rrrm1_blood_single_cell_spaceflight")
        self.assertEqual(summary["selected_source_id"], "OSD-918")
        self.assertEqual(summary["selected_tissue"], "blood")
        self.assertEqual(summary["local_payload_present"], "false")
        self.assertEqual(summary["required_obs_field_count"], "9")
        self.assertEqual(summary["required_var_field_count"], "2")
        self.assertEqual(summary["required_uns_field_count"], "4")
        self.assertEqual(summary["staging_action_count"], "7")
        self.assertEqual(
            summary["next_required_block"],
            "V9-SC-005: AnnData obs/var audit implementation",
        )
        self.assertIn("canonical_v9_payload_target", {row["candidate_id"] for row in candidates})
        self.assertIn("legacy_annotated_h5ad", {row["candidate_id"] for row in candidates})
        self.assertIn("broad_celltype", {row["field_id"] for row in requirements})
        self.assertIn("gene_symbol", {row["field_id"] for row in requirements})
        self.assertEqual(len(actions), 7)
        self.assertEqual(len(package["summary"]), 1)
        self.assertIn("No local h5ad payload is claimed", review)

    def test_build_v9_sc_payload_staging_plan_cli(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "sc_spaceflight"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "build_v9_sc_payload_staging_plan.py"),
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
            with (output_dir / "payload_staging_plan_summary.csv").open(
                newline=""
            ) as handle:
                summary_rows = list(csv.DictReader(handle))
            with (output_dir / "obs_var_audit_requirements.csv").open(
                newline=""
            ) as handle:
                requirement_rows = list(csv.DictReader(handle))

        self.assertIn("payload_staging_plan_summary.csv", result.stdout)
        self.assertIn("V9_SC_PAYLOAD_STAGING_PLAN.md", result.stdout)
        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(summary_rows[0]["local_payload_present"], "false")
        self.assertIn("flight_ground_label", {row["field_id"] for row in requirement_rows})

    def test_generated_v9_sc_payload_staging_plan_outputs_validate(self):
        output_dir = REPO_ROOT / "v9" / "sc_spaceflight"
        with (output_dir / "payload_staging_plan_summary.csv").open(
            newline=""
        ) as handle:
            summary_rows = list(csv.DictReader(handle))
        with (output_dir / "payload_staging_candidates.csv").open(
            newline=""
        ) as handle:
            candidate_rows = list(csv.DictReader(handle))
        with (output_dir / "obs_var_audit_requirements.csv").open(
            newline=""
        ) as handle:
            requirement_rows = list(csv.DictReader(handle))
        with (output_dir / "payload_staging_actions.csv").open(
            newline=""
        ) as handle:
            action_rows = list(csv.DictReader(handle))
        review = (REPO_ROOT / "docs" / "V9_SC_PAYLOAD_STAGING_PLAN.md").read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(
            summary["claim_boundary"],
            "payload_staging_plan_only_no_local_payload_or_score_claim",
        )
        self.assertEqual(summary["canonical_payload_path"], "v9/sc_spaceflight/payloads/rrrm1_blood/OSD-918_blood_rrrm1_bench.h5ad")
        self.assertEqual(len(candidate_rows), 4)
        self.assertEqual(len(action_rows), 7)
        by_candidate = {row["candidate_id"]: row for row in candidate_rows}
        self.assertEqual(
            by_candidate["canonical_v9_payload_target"]["path_status"],
            "planned_repo_path_absent",
        )
        self.assertIn(
            "payload_not_runnable",
            {row["blocker_if_missing"] for row in requirement_rows},
        )
        self.assertIn("V9-SC-005: AnnData obs/var audit implementation", review)

    def test_v9_sc_obs_var_audit_api_skips_without_payload(self):
        from spacebio_bench import build_sc_obs_var_audit, write_sc_obs_var_audit

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "sc_spaceflight"
            outputs = write_sc_obs_var_audit(
                repo_root=REPO_ROOT,
                output_dir=output_dir,
            )
            with outputs["summary_csv"].open(newline="") as handle:
                summary_rows = list(csv.DictReader(handle))
            results = json.loads(outputs["results_json"].read_text())
            payload_manifest = json.loads(outputs["payload_manifest_json"].read_text())
            review = outputs["review_md"].read_text()

        package = build_sc_obs_var_audit(repo_root=REPO_ROOT)
        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(
            summary["decision_status"],
            "obs_var_audit_skipped_no_local_payload",
        )
        self.assertEqual(summary["task_id"], "draft_rrrm1_blood_single_cell_spaceflight")
        self.assertEqual(summary["payload_path_status"], "missing")
        self.assertEqual(summary["requirement_count"], "27")
        self.assertEqual(summary["pass_count"], "0")
        self.assertEqual(summary["fail_count"], "0")
        self.assertEqual(summary["skip_count"], "27")
        self.assertEqual(summary["blocker_count"], "17")
        self.assertEqual(
            summary["next_required_block"],
            "V9-SC-006: canonical payload staging or RRRM-1 h5ad regeneration",
        )
        self.assertEqual({row["audit_status"] for row in results}, {"skipped_no_local_payload"})
        self.assertEqual(payload_manifest[0]["path_status"], "missing")
        self.assertEqual(len(package["results"]), 27)
        self.assertIn("machine-readable skip audit", review)

    def test_audit_v9_sc_obs_var_cli(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "sc_spaceflight"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "audit_v9_sc_obs_var.py"),
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
            with (output_dir / "obs_var_audit_summary.csv").open(
                newline=""
            ) as handle:
                summary_rows = list(csv.DictReader(handle))
            with (output_dir / "payload_manifest.csv").open(newline="") as handle:
                payload_rows = list(csv.DictReader(handle))

        self.assertIn("obs_var_audit_summary.csv", result.stdout)
        self.assertIn("V9_SC_OBS_VAR_AUDIT.md", result.stdout)
        self.assertEqual(summary_rows[0]["payload_path_status"], "missing")
        self.assertEqual(payload_rows[0]["h5ad_read_status"], "not_attempted_no_local_payload")

    def test_generated_v9_sc_obs_var_audit_outputs_validate(self):
        output_dir = REPO_ROOT / "v9" / "sc_spaceflight"
        with (output_dir / "obs_var_audit_summary.csv").open(newline="") as handle:
            summary_rows = list(csv.DictReader(handle))
        with (output_dir / "obs_var_audit_results.csv").open(newline="") as handle:
            result_rows = list(csv.DictReader(handle))
        with (output_dir / "payload_manifest.csv").open(newline="") as handle:
            payload_rows = list(csv.DictReader(handle))
        review = (REPO_ROOT / "docs" / "V9_SC_OBS_VAR_AUDIT.md").read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(
            summary["claim_boundary"],
            "obs_var_audit_skip_only_no_payload_or_score_claim",
        )
        self.assertEqual(summary["payload_sha256"], "NA")
        self.assertEqual(summary["n_obs"], "NA")
        self.assertEqual(len(result_rows), 27)
        self.assertEqual({row["audit_status"] for row in result_rows}, {"skipped_no_local_payload"})
        self.assertEqual(payload_rows[0]["path_status"], "missing")
        self.assertIn("No obs/var pass is claimed", review)

    def test_v9_sc_payload_staging_execution_api_blocks_without_candidate_payload(self):
        from spacebio_bench import (
            build_sc_payload_staging_execution,
            write_sc_payload_staging_execution,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "sc_spaceflight"
            outputs = write_sc_payload_staging_execution(
                repo_root=REPO_ROOT,
                output_dir=output_dir,
            )
            with outputs["summary_csv"].open(newline="") as handle:
                summary_rows = list(csv.DictReader(handle))
            candidates = json.loads(outputs["candidates_json"].read_text())
            regeneration_steps = json.loads(outputs["regeneration_json"].read_text())
            review = outputs["review_md"].read_text()

        package = build_sc_payload_staging_execution(repo_root=REPO_ROOT)
        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(
            summary["decision_status"],
            "payload_staging_execution_blocked_no_candidate_payload",
        )
        self.assertEqual(summary["task_id"], "draft_rrrm1_blood_single_cell_spaceflight")
        self.assertEqual(summary["canonical_payload_status"], "planned_repo_path_absent")
        self.assertEqual(summary["selected_route"], "prepare_regeneration_from_starsolo_per_srx")
        self.assertEqual(summary["payload_manifest_status"], "missing")
        self.assertEqual(summary["obs_var_audit_status"], "obs_var_audit_skipped_no_local_payload")
        self.assertEqual(summary["local_payload_staged"], "false")
        self.assertEqual(summary["candidate_count"], "4")
        self.assertEqual(summary["regeneration_step_count"], "5")
        self.assertEqual(
            summary["claim_boundary"],
            "payload_staging_execution_no_payload_or_score_claim",
        )
        by_candidate = {row["candidate_id"]: row for row in candidates}
        self.assertEqual(
            by_candidate["legacy_annotated_h5ad"]["action_decision"],
            "preferred_route_blocked",
        )
        self.assertEqual(
            by_candidate["starsolo_per_srx_regeneration"]["action_decision"],
            "regeneration_route_requires_per_srx_matrices",
        )
        self.assertIn("confirm_per_srx_starsolo_matrices", {row["step_id"] for row in regeneration_steps})
        self.assertEqual(len(package["regeneration_steps"]), 5)
        self.assertIn("No local h5ad payload is claimed", review)

    def test_build_v9_sc_payload_staging_execution_cli(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "sc_spaceflight"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "build_v9_sc_payload_staging_execution.py"),
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
            with (output_dir / "payload_staging_execution_summary.csv").open(
                newline=""
            ) as handle:
                summary_rows = list(csv.DictReader(handle))
            with (output_dir / "payload_regeneration_steps.csv").open(
                newline=""
            ) as handle:
                regeneration_rows = list(csv.DictReader(handle))

        self.assertIn("payload_staging_execution_summary.csv", result.stdout)
        self.assertIn("V9_SC_PAYLOAD_STAGING_EXECUTION.md", result.stdout)
        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(summary_rows[0]["local_payload_staged"], "false")
        self.assertEqual(len(regeneration_rows), 5)

    def test_generated_v9_sc_payload_staging_execution_outputs_validate(self):
        output_dir = REPO_ROOT / "v9" / "sc_spaceflight"
        with (output_dir / "payload_staging_execution_summary.csv").open(
            newline=""
        ) as handle:
            summary_rows = list(csv.DictReader(handle))
        with (output_dir / "payload_staging_execution_candidates.csv").open(
            newline=""
        ) as handle:
            candidate_rows = list(csv.DictReader(handle))
        with (output_dir / "payload_regeneration_steps.csv").open(
            newline=""
        ) as handle:
            regeneration_rows = list(csv.DictReader(handle))
        review = (REPO_ROOT / "docs" / "V9_SC_PAYLOAD_STAGING_EXECUTION.md").read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(
            summary["claim_boundary"],
            "payload_staging_execution_no_payload_or_score_claim",
        )
        self.assertEqual(summary["canonical_payload_status"], "planned_repo_path_absent")
        self.assertEqual(summary["payload_manifest_status"], "missing")
        self.assertEqual(summary["obs_var_audit_status"], "obs_var_audit_skipped_no_local_payload")
        self.assertEqual(len(candidate_rows), 4)
        self.assertEqual(len(regeneration_rows), 5)
        self.assertIn(
            "blocked_until_external_payload_available",
            {row["gate_status"] for row in regeneration_rows},
        )
        self.assertIn("No evaluator, leaderboard result, or legacy RRRM score promotion is claimed", review)

    def test_v9_sc_external_payload_availability_api_blocks_without_source_payload(self):
        from spacebio_bench import (
            build_sc_external_payload_availability,
            write_sc_external_payload_availability,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "sc_spaceflight"
            checked_bases = [
                str(Path(tmpdir) / "checked_base_a"),
                str(Path(tmpdir) / "checked_base_b"),
            ]
            outputs = write_sc_external_payload_availability(
                repo_root=REPO_ROOT,
                output_dir=output_dir,
                checked_bases=checked_bases,
            )
            with outputs["summary_csv"].open(newline="") as handle:
                summary_rows = list(csv.DictReader(handle))
            candidates = json.loads(outputs["candidates_json"].read_text())
            matrix_rows = json.loads(outputs["matrix_json"].read_text())
            copy_decision = json.loads(outputs["copy_decision_json"].read_text())
            review = outputs["review_md"].read_text()
            package = build_sc_external_payload_availability(
                repo_root=REPO_ROOT,
                checked_bases=checked_bases,
            )

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(
            summary["decision_status"],
            "external_payload_availability_blocked_no_h5ad_or_starsolo_matrices",
        )
        self.assertEqual(summary["task_id"], "draft_rrrm1_blood_single_cell_spaceflight")
        self.assertEqual(summary["source_id"], "OSD-918")
        self.assertEqual(summary["tissue"], "blood")
        self.assertEqual(summary["checked_base_count"], "2")
        self.assertEqual(summary["expected_blood_srx_count"], "8")
        self.assertEqual(summary["checked_matrix_row_count"], "16")
        self.assertEqual(summary["annotated_h5ad_found"], "false")
        self.assertEqual(summary["labeled_h5ad_found"], "false")
        self.assertEqual(summary["complete_starsolo_srx_count"], "0")
        self.assertEqual(summary["canonical_copy_allowed"], "false")
        self.assertEqual(summary["regeneration_allowed"], "false")
        self.assertEqual(
            summary["next_required_block"],
            "V9-SC-006c: OSDR processed payload discovery or owner scratch path request",
        )
        self.assertEqual(
            summary["claim_boundary"],
            "external_payload_availability_no_payload_or_score_claim",
        )
        self.assertEqual(len(candidates), 5)
        self.assertEqual(
            {row["candidate_type"] for row in candidates},
            {"canonical_target", "preferred_annotated_h5ad", "metadata_labeled_h5ad"},
        )
        self.assertEqual(len(matrix_rows), 16)
        self.assertEqual({row["complete_matrix_bundle"] for row in matrix_rows}, {"false"})
        self.assertEqual(len({row["srx"] for row in matrix_rows}), 8)
        self.assertEqual(copy_decision[0]["source_payload_status"], "no_source_payload_found")
        self.assertEqual(copy_decision[0]["copy_allowed"], "false")
        self.assertEqual(len(package["matrix_availability"]), 16)
        self.assertIn("No canonical or external RRRM-1 blood source payload", review)

    def test_build_v9_sc_external_payload_availability_cli(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "sc_spaceflight"
            checked_bases = [
                str(Path(tmpdir) / "checked_base_a"),
                str(Path(tmpdir) / "checked_base_b"),
            ]
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "build_v9_sc_external_payload_availability.py"),
                    "--repo-root",
                    str(REPO_ROOT),
                    "--output-dir",
                    str(output_dir),
                    "--checked-base",
                    checked_bases[0],
                    "--checked-base",
                    checked_bases[1],
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            with (output_dir / "external_payload_availability_summary.csv").open(
                newline=""
            ) as handle:
                summary_rows = list(csv.DictReader(handle))
            with (output_dir / "external_payload_candidates.csv").open(
                newline=""
            ) as handle:
                candidate_rows = list(csv.DictReader(handle))
            with (output_dir / "external_starsolo_matrix_availability.csv").open(
                newline=""
            ) as handle:
                matrix_rows = list(csv.DictReader(handle))

        self.assertIn("external_payload_availability_summary.csv", result.stdout)
        self.assertIn("V9_SC_EXTERNAL_PAYLOAD_AVAILABILITY.md", result.stdout)
        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(summary_rows[0]["checked_matrix_row_count"], "16")
        self.assertEqual(summary_rows[0]["canonical_copy_allowed"], "false")
        self.assertEqual(len(candidate_rows), 5)
        self.assertEqual(len(matrix_rows), 16)

    def test_generated_v9_sc_external_payload_availability_outputs_validate(self):
        output_dir = REPO_ROOT / "v9" / "sc_spaceflight"
        with (output_dir / "external_payload_availability_summary.csv").open(
            newline=""
        ) as handle:
            summary_rows = list(csv.DictReader(handle))
        with (output_dir / "external_payload_candidates.csv").open(
            newline=""
        ) as handle:
            candidate_rows = list(csv.DictReader(handle))
        with (output_dir / "external_starsolo_matrix_availability.csv").open(
            newline=""
        ) as handle:
            matrix_rows = list(csv.DictReader(handle))
        with (output_dir / "canonical_payload_copy_decision.csv").open(
            newline=""
        ) as handle:
            copy_rows = list(csv.DictReader(handle))
        review = (REPO_ROOT / "docs" / "V9_SC_EXTERNAL_PAYLOAD_AVAILABILITY.md").read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(
            summary["decision_status"],
            "external_payload_availability_blocked_no_h5ad_or_starsolo_matrices",
        )
        self.assertEqual(
            summary["claim_boundary"],
            "external_payload_availability_no_payload_or_score_claim",
        )
        self.assertEqual(summary["expected_blood_srx_count"], "8")
        self.assertEqual(summary["checked_matrix_row_count"], "16")
        self.assertEqual(summary["canonical_copy_allowed"], "false")
        self.assertEqual(summary["regeneration_allowed"], "false")
        self.assertEqual(len(candidate_rows), 5)
        self.assertEqual(len(matrix_rows), 16)
        self.assertEqual({row["complete_matrix_bundle"] for row in matrix_rows}, {"false"})
        self.assertEqual(copy_rows[0]["copy_allowed"], "false")
        self.assertIn("No h5ad was copied into the canonical v9 payload path", review)

    def test_v9_sc_processed_payload_discovery_api_blocks_without_processed_payload(self):
        from spacebio_bench import (
            build_sc_processed_payload_discovery,
            write_sc_processed_payload_discovery,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fixture = tmp / "osd918_files.json"
            write_osd918_osdr_file_listing_fixture(fixture)
            output_dir = tmp / "sc_spaceflight"
            outputs = write_sc_processed_payload_discovery(
                repo_root=REPO_ROOT,
                output_dir=output_dir,
                api_listing_json=fixture,
            )
            with outputs["summary_csv"].open(newline="") as handle:
                summary_rows = list(csv.DictReader(handle))
            files = json.loads(outputs["files_json"].read_text())
            coverage = json.loads(outputs["coverage_json"].read_text())
            owner_requests = json.loads(outputs["owner_request_json"].read_text())
            deferral = json.loads(outputs["deferral_json"].read_text())
            review = outputs["review_md"].read_text()
            package = build_sc_processed_payload_discovery(
                repo_root=REPO_ROOT,
                api_listing_json=fixture,
            )

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(
            summary["decision_status"],
            "osdr_processed_payload_unavailable_owner_scratch_request_required",
        )
        self.assertEqual(summary["task_id"], "draft_rrrm1_blood_single_cell_spaceflight")
        self.assertEqual(summary["source_id"], "OSD-918")
        self.assertEqual(summary["glds_prefix"], "GLDS-746")
        self.assertEqual(summary["tissue"], "blood")
        self.assertEqual(summary["api_status"], "ok")
        self.assertEqual(summary["osdr_file_count"], "19")
        self.assertEqual(summary["metadata_file_count"], "1")
        self.assertEqual(summary["raw_fastq_count"], "16")
        self.assertEqual(summary["raw_checksum_count"], "1")
        self.assertEqual(summary["raw_multiqc_count"], "1")
        self.assertEqual(summary["processed_h5ad_count"], "0")
        self.assertEqual(summary["processed_starsolo_count"], "0")
        self.assertEqual(summary["processed_checksum_manifest_count"], "0")
        self.assertEqual(summary["expected_blood_srx_count"], "8")
        self.assertEqual(summary["complete_expected_fastq_pair_count"], "8")
        self.assertEqual(summary["missing_expected_fastq_pair_count"], "0")
        self.assertEqual(summary["canonical_copy_allowed"], "false")
        self.assertEqual(summary["regeneration_allowed"], "false")
        self.assertEqual(summary["owner_scratch_request_required"], "true")
        self.assertEqual(
            summary["next_required_block"],
            "V9-SC-006d: owner scratch path intake or raw FASTQ regeneration feasibility decision",
        )
        self.assertEqual(
            summary["claim_boundary"],
            "osdr_processed_payload_discovery_only_no_payload_copy_or_score_claim",
        )
        self.assertEqual(len(files), 19)
        self.assertEqual(len(coverage), 8)
        self.assertEqual({row["fastq_pair_complete"] for row in coverage}, {"true"})
        self.assertEqual(
            {row["processed_matrix_status"] for row in coverage},
            {"not_listed_in_osdr_file_api"},
        )
        self.assertEqual(len(owner_requests), 3)
        self.assertEqual(len(deferral), 1)
        self.assertEqual(deferral[0]["copy_allowed"], "false")
        self.assertEqual(len(package["files"]), 19)
        self.assertIn("No OSDR payload was downloaded", review)

    def test_build_v9_sc_osdr_processed_payload_discovery_cli(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fixture = tmp / "osd918_files.json"
            write_osd918_osdr_file_listing_fixture(fixture)
            output_dir = tmp / "sc_spaceflight"
            result = subprocess.run(
                [
                    sys.executable,
                    str(
                        REPO_ROOT
                        / "scripts"
                        / "build_v9_sc_osdr_processed_payload_discovery.py"
                    ),
                    "--repo-root",
                    str(REPO_ROOT),
                    "--output-dir",
                    str(output_dir),
                    "--api-listing-json",
                    str(fixture),
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            with (output_dir / "osdr_processed_payload_discovery_summary.csv").open(
                newline=""
            ) as handle:
                summary_rows = list(csv.DictReader(handle))
            with (output_dir / "osdr_file_discovery.csv").open(newline="") as handle:
                file_rows = list(csv.DictReader(handle))
            with (output_dir / "owner_scratch_request.csv").open(newline="") as handle:
                request_rows = list(csv.DictReader(handle))

        self.assertIn("osdr_processed_payload_discovery_summary.csv", result.stdout)
        self.assertIn("V9_SC_OSDR_PROCESSED_PAYLOAD_DISCOVERY.md", result.stdout)
        self.assertEqual(len(summary_rows), 1)
        self.assertEqual(summary_rows[0]["osdr_file_count"], "19")
        self.assertEqual(summary_rows[0]["owner_scratch_request_required"], "true")
        self.assertEqual(len(file_rows), 19)
        self.assertEqual(len(request_rows), 3)

    def test_generated_v9_sc_osdr_processed_payload_discovery_outputs_validate(self):
        output_dir = REPO_ROOT / "v9" / "sc_spaceflight"
        with (output_dir / "osdr_processed_payload_discovery_summary.csv").open(
            newline=""
        ) as handle:
            summary_rows = list(csv.DictReader(handle))
        with (output_dir / "osdr_file_discovery.csv").open(newline="") as handle:
            file_rows = list(csv.DictReader(handle))
        with (output_dir / "osdr_expected_srx_coverage.csv").open(newline="") as handle:
            coverage_rows = list(csv.DictReader(handle))
        with (output_dir / "owner_scratch_request.csv").open(newline="") as handle:
            request_rows = list(csv.DictReader(handle))
        with (output_dir / "processed_payload_deferral_decision.csv").open(
            newline=""
        ) as handle:
            deferral_rows = list(csv.DictReader(handle))
        review = (
            REPO_ROOT / "docs" / "V9_SC_OSDR_PROCESSED_PAYLOAD_DISCOVERY.md"
        ).read_text()

        self.assertEqual(len(summary_rows), 1)
        summary = summary_rows[0]
        self.assertEqual(
            summary["decision_status"],
            "osdr_processed_payload_unavailable_owner_scratch_request_required",
        )
        self.assertEqual(summary["osdr_file_count"], "19")
        self.assertEqual(summary["raw_fastq_count"], "16")
        self.assertEqual(summary["processed_h5ad_count"], "0")
        self.assertEqual(summary["processed_starsolo_count"], "0")
        self.assertEqual(summary["complete_expected_fastq_pair_count"], "8")
        self.assertEqual(summary["canonical_copy_allowed"], "false")
        self.assertEqual(summary["regeneration_allowed"], "false")
        self.assertEqual(summary["owner_scratch_request_required"], "true")
        self.assertEqual(
            summary["claim_boundary"],
            "osdr_processed_payload_discovery_only_no_payload_copy_or_score_claim",
        )
        self.assertEqual(len(file_rows), 19)
        self.assertEqual(len(coverage_rows), 8)
        self.assertEqual(len(request_rows), 3)
        self.assertEqual(len(deferral_rows), 1)
        self.assertEqual(deferral_rows[0]["copy_allowed"], "false")
        self.assertIn("Complete expected raw FASTQ pairs listed by OSDR: 8/8", review)

    def test_generated_v9_multispecies_task_manifests_validate(self):
        from spacebio_bench import TaskRegistry, load_task_manifest

        manifest_dir = REPO_ROOT / "v9" / "multispecies" / "task_manifests"
        manifests = sorted(manifest_dir.glob("*.json"))
        self.assertEqual(len(manifests), 2)
        payloads = [load_task_manifest(path) for path in manifests]
        self.assertEqual(
            {payload["task_id"] for payload in payloads},
            {
                "draft_osd37_arabidopsis_seedling_spaceflight",
                "draft_osd207_drosophila_whole_body_spaceflight",
            },
        )
        self.assertEqual(
            {payload["provenance"]["source_checksum_status"] for payload in payloads},
            {"checksum_manifest_parsed"},
        )
        registry = TaskRegistry.from_dir(manifest_dir)
        self.assertEqual(len(registry), 2)

    def test_multispecies_loader_aligns_species_native_matrices_and_folds(self):
        from spacebio_bench import load_multispecies_task

        arabidopsis = load_multispecies_task(
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "task_manifests"
            / "draft_osd37_arabidopsis_seedling_spaceflight.json",
            repo_root=REPO_ROOT,
        )
        drosophila = load_multispecies_task(
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "task_manifests"
            / "draft_osd207_drosophila_whole_body_spaceflight.json",
            repo_root=REPO_ROOT,
        )

        self.assertEqual(arabidopsis.n_samples, 56)
        self.assertEqual(arabidopsis.n_features, 28067)
        self.assertEqual(len(arabidopsis.folds), 4)
        self.assertEqual(set(arabidopsis.matrix_paths), {"OSD-37"})
        self.assertEqual(drosophila.n_samples, 32)
        self.assertEqual(drosophila.n_features, 15999)
        self.assertEqual(len(drosophila.folds), 4)
        self.assertEqual(set(drosophila.matrix_paths), {"OSD-207"})
        self.assertEqual(
            set(arabidopsis.features.index),
            {row["sample_id"] for row in arabidopsis.sample_factors},
        )

    def test_multispecies_nearest_centroid_baseline_writes_draft_report(self):
        from spacebio_bench import (
            MultispeciesBaselineConfig,
            load_multispecies_task,
            run_multispecies_nearest_centroid,
        )

        task = load_multispecies_task(
            REPO_ROOT
            / "v9"
            / "multispecies"
            / "task_manifests"
            / "draft_osd207_drosophila_whole_body_spaceflight.json",
            repo_root=REPO_ROOT,
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            result = run_multispecies_nearest_centroid(
                task,
                output_dir=Path(tmpdir) / "multispecies",
                config=MultispeciesBaselineConfig(top_variable_genes=100),
                task_manifest_path=(
                    REPO_ROOT
                    / "v9"
                    / "multispecies"
                    / "task_manifests"
                    / "draft_osd207_drosophila_whole_body_spaceflight.json"
                ),
                command=["test-multispecies-baseline"],
            )
            with result.predictions_path.open(newline="") as handle:
                predictions = list(csv.DictReader(handle))

        self.assertEqual(result.n_predictions, 32)
        self.assertEqual(len(predictions), 32)
        self.assertEqual(result.evaluation["status"], "evaluated")
        self.assertEqual(result.evaluation["positive_label"], "LEO_or_ISS")
        self.assertEqual(result.evaluation["baseline"]["release_status"], "draft_not_frozen")
        self.assertEqual(
            result.evaluation["baseline"]["fold_family"],
            "condition_stratum_candidate_folds",
        )
        self.assertEqual(
            result.evaluation["metrics"]["condition_stratum_holdout_delta"]["status"],
            "computed",
        )

    def test_run_v9_multispecies_baseline_cli_writes_summary(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "nearest_centroid"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "run_v9_multispecies_baseline.py"),
                    "--repo-root",
                    str(REPO_ROOT),
                    "--manifest-dir",
                    str(REPO_ROOT / "v9" / "multispecies" / "task_manifests"),
                    "--output-dir",
                    str(output_dir),
                    "--top-variable-genes",
                    "100",
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            with (output_dir / "multispecies_baseline_summary.csv").open(
                newline=""
            ) as handle:
                rows = list(csv.DictReader(handle))

        self.assertIn("multispecies_baseline_summary.csv", result.stdout)
        self.assertEqual(len(rows), 2)
        self.assertEqual({row["baseline_id"] for row in rows}, {"multispecies_nearest_centroid"})
        self.assertEqual({row["status"] for row in rows}, {"evaluated"})
        self.assertEqual({row["claim_boundary"] for row in rows}, {"draft_feasibility_only_not_leaderboard"})

    def test_multispecies_sensitivity_grid_writes_variant_summary(self):
        from spacebio_bench import (
            MultispeciesBaselineConfig,
            run_multispecies_sensitivity_grid,
        )

        configs = [
            MultispeciesBaselineConfig(transform="log1p", scaling="zscore", top_variable_genes=100),
            MultispeciesBaselineConfig(transform="none", scaling="none", top_variable_genes=100),
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            _, summary = run_multispecies_sensitivity_grid(
                repo_root=REPO_ROOT,
                manifest_dir=REPO_ROOT / "v9" / "multispecies" / "task_manifests",
                output_root=Path(tmpdir) / "sensitivity",
                configs=configs,
                command=["test-multispecies-sensitivity"],
            )
            with summary["csv"].open(newline="") as handle:
                rows = list(csv.DictReader(handle))

        self.assertEqual(len(rows), 4)
        self.assertEqual(
            {row["variant_id"] for row in rows},
            {"tvg100_log1p_zscore", "tvg100_none_none"},
        )
        self.assertEqual(
            {row["task_id"] for row in rows},
            {
                "draft_osd37_arabidopsis_seedling_spaceflight",
                "draft_osd207_drosophila_whole_body_spaceflight",
            },
        )
        self.assertTrue(all(row["status"] == "evaluated" for row in rows))

    def test_run_v9_multispecies_sensitivity_cli_writes_summary(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "sensitivity"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "run_v9_multispecies_sensitivity.py"),
                    "--repo-root",
                    str(REPO_ROOT),
                    "--manifest-dir",
                    str(REPO_ROOT / "v9" / "multispecies" / "task_manifests"),
                    "--output-dir",
                    str(output_dir),
                    "--top-variable-genes",
                    "100",
                    "--transform",
                    "log1p",
                    "--scaling",
                    "zscore",
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            with (output_dir / "multispecies_baseline_summary.csv").open(
                newline=""
            ) as handle:
                rows = list(csv.DictReader(handle))

        self.assertIn("multispecies_baseline_summary.csv", result.stdout)
        self.assertEqual(len(rows), 2)
        self.assertEqual({row["variant_id"] for row in rows}, {"tvg100_log1p_zscore"})
        self.assertEqual({row["status"] for row in rows}, {"evaluated"})

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

    def test_organoid_de_reference_extractor_normalizes_direct_spaceflight_sign(self):
        from spacebio_bench.organoid_de_reference import build_organoid_de_reference
        from spacebio_bench.source_audit import FetchResult

        direct_contrast = (
            "(control & Ground Control & with microglia)v"
            "(control & Space Flight & with microglia)"
        )
        de_url = "https://example.test/de.csv"
        de_csv = (
            '"ENSEMBL","SYMBOL","GENENAME","REFSEQ","ENTREZID","STRING_id",'
            '"GOSLIM_IDS",'
            f'"Log2fc_{direct_contrast}","Stat_{direct_contrast}",'
            f'"P.value_{direct_contrast}","Adj.p.value_{direct_contrast}"\n'
            '"string","string","string","string","number","string","string",'
            '"number","number","number","number"\n'
            '"ENSG000001","GENE1","Gene 1","NM_1","101","","",'
            '"2.0","5.0","0.01","0.02"\n'
            '"ENSG000002","GENE2","Gene 2","NM_2","102","","",'
            '"-1.5","-3.0","0.2","0.5"\n'
        ).encode()
        source_rows = [
            {
                "source_id": "OSD-TEST",
                "glds_prefix": "GLDS-TEST",
                "organoid_type": "mock_neural_organoid",
            }
        ]
        audit_rows = [
            {
                "source_id": "OSD-TEST",
                "glds_prefix": "GLDS-TEST",
                "task_ids": "draft_human_organoid_spaceflight",
                "de_reference_files": (
                    "GLDS-TEST_rna_seq_differential_expression_GLbulkRNAseq.csv;"
                    "GLDS-TEST_rna_seq_differential_expression_rRNArm_GLbulkRNAseq.csv"
                ),
                "de_reference_urls": f"{de_url};https://example.test/rrna.csv",
                "de_reference_file_sizes": "123;456",
                "de_reference_data_types": (
                    "differential expression table;"
                    "differential expression table - supplementary"
                ),
                "direct_spaceflight_contrasts": direct_contrast,
            }
        ]

        def fake_fetcher(url: str) -> FetchResult:
            if url == de_url:
                return FetchResult(ok=True, url=url, body=de_csv)
            return FetchResult(ok=False, url=url, error="unexpected url")

        rows, manifest = build_organoid_de_reference(
            source_rows=source_rows,
            audit_rows=audit_rows,
            fetcher=fake_fetcher,
        )

        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["source_de_file"], "GLDS-TEST_rna_seq_differential_expression_GLbulkRNAseq.csv")
        self.assertEqual(rows[0]["source_log2fc_value"], "2")
        self.assertEqual(rows[0]["log2fc_leo_or_iss_minus_ground"], "-2")
        self.assertEqual(rows[0]["significant_fdr_0_05"], "true")
        self.assertEqual(rows[1]["log2fc_leo_or_iss_minus_ground"], "1.5")
        self.assertEqual(rows[1]["significant_fdr_0_05"], "false")
        self.assertEqual(manifest["totals"]["n_rows"], 2)
        self.assertEqual(manifest["totals"]["n_significant_fdr_0_05"], 1)
        self.assertIn("response_signature_contract", manifest)

    def test_write_organoid_de_reference_writes_csv_and_manifest(self):
        from spacebio_bench.organoid_de_reference import write_organoid_de_reference
        from spacebio_bench.source_audit import FetchResult

        direct_contrast = (
            "(control & Ground Control & without microglia)v"
            "(control & Space Flight & without microglia)"
        )
        de_csv = (
            "ENSEMBL,SYMBOL,GENENAME,REFSEQ,ENTREZID,"
            f"Log2fc_{direct_contrast},Stat_{direct_contrast},"
            f"P.value_{direct_contrast},Adj.p.value_{direct_contrast}\n"
            "ENSG000001,GENE1,Gene 1,NM_1,101,0.25,2.0,0.03,0.04\n"
        ).encode()
        source_rows = [{"source_id": "OSD-WRITE", "glds_prefix": "GLDS-WRITE"}]
        audit_rows = [
            {
                "source_id": "OSD-WRITE",
                "glds_prefix": "GLDS-WRITE",
                "de_reference_files": "GLDS-WRITE_rna_seq_differential_expression_GLbulkRNAseq.csv",
                "de_reference_urls": "https://example.test/write-de.csv",
                "direct_spaceflight_contrasts": direct_contrast,
            }
        ]

        def fake_fetcher(url: str) -> FetchResult:
            return FetchResult(ok=True, url=url, body=de_csv)

        with tempfile.TemporaryDirectory() as tmpdir:
            csv_path, json_path = write_organoid_de_reference(
                source_rows=source_rows,
                audit_rows=audit_rows,
                csv_path=Path(tmpdir) / "de_reference.csv",
                json_path=Path(tmpdir) / "de_reference_manifest.json",
                fetcher=fake_fetcher,
            )
            with csv_path.open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            manifest = json.loads(json_path.read_text())

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["contrast_id"], "osd_write_control_without_microglia_leo_or_iss_vs_ground")
        self.assertEqual(manifest["totals"]["n_sources"], 1)
        self.assertEqual(
            manifest["log2fc_orientation_contract"]["benchmark_orientation"],
            "LEO_or_ISS_minus_Ground",
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
        self.assertEqual(
            manifest["reference_signatures"]["response_signature_contract"][
                "predicted_log2fc_orientation"
            ],
            "LEO_or_ISS_minus_Ground",
        )

    def test_organoid_evaluator_skips_signature_metrics_without_response_artifact(self):
        from spacebio_bench import evaluate_submission, load_task_manifest

        manifest = load_task_manifest(
            REPO_ROOT
            / "v9"
            / "human_organoid"
            / "task_manifests"
            / "draft_human_organoid_spaceflight.json"
        )
        rows = [
            {
                "sample_id": "s1",
                "true_label": "Ground",
                "predicted_label": "Ground",
                "flight_probability": 0.2,
            },
            {
                "sample_id": "s2",
                "true_label": "LEO_or_ISS",
                "predicted_label": "LEO_or_ISS",
                "flight_probability": 0.8,
            },
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            submission = Path(tmpdir) / "submission.csv"
            self.write_submission(submission, rows)
            result = evaluate_submission(manifest, submission)

        self.assertEqual(result["status"], "evaluated")
        self.assertEqual(result["metrics"]["de_direction_match"]["status"], "skipped")
        self.assertEqual(result["metrics"]["signature_rank_correlation"]["status"], "skipped")
        self.assertIn(
            "response_signature.csv artifact missing",
            result["metrics"]["de_direction_match"]["reason"],
        )

    def test_response_signature_metrics_compute_direction_and_rank(self):
        from spacebio_bench import compute_response_signature_metrics

        manifest = {
            "task_id": "draft_human_organoid_spaceflight",
            "reference_signatures": {},
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            reference_path = tmp / "reference.csv.gz"
            with gzip.open(reference_path, "wt", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "source_id",
                        "contrast_id",
                        "gene_symbol",
                        "ensembl_id",
                        "log2fc_leo_or_iss_minus_ground",
                        "significant_fdr_0_05",
                    ],
                )
                writer.writeheader()
                writer.writerows(
                    [
                        {
                            "source_id": "OSD-X",
                            "contrast_id": "contrast_a",
                            "gene_symbol": "GENE1",
                            "ensembl_id": "ENSG1",
                            "log2fc_leo_or_iss_minus_ground": "1.0",
                            "significant_fdr_0_05": "true",
                        },
                        {
                            "source_id": "OSD-X",
                            "contrast_id": "contrast_a",
                            "gene_symbol": "GENE2",
                            "ensembl_id": "ENSG2",
                            "log2fc_leo_or_iss_minus_ground": "-1.0",
                            "significant_fdr_0_05": "true",
                        },
                        {
                            "source_id": "OSD-X",
                            "contrast_id": "contrast_a",
                            "gene_symbol": "GENE3",
                            "ensembl_id": "ENSG3",
                            "log2fc_leo_or_iss_minus_ground": "0.2",
                            "significant_fdr_0_05": "false",
                        },
                    ]
                )
            response_path = tmp / "response_signature.csv"
            with response_path.open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "task_id",
                        "source_id",
                        "contrast_id",
                        "gene_symbol",
                        "ensembl_id",
                        "predicted_log2fc_leo_or_iss_minus_ground",
                    ],
                )
                writer.writeheader()
                writer.writerows(
                    [
                        {
                            "task_id": "draft_human_organoid_spaceflight",
                            "source_id": "OSD-X",
                            "contrast_id": "contrast_a",
                            "gene_symbol": "GENE1",
                            "ensembl_id": "ENSG1",
                            "predicted_log2fc_leo_or_iss_minus_ground": "2.0",
                        },
                        {
                            "task_id": "draft_human_organoid_spaceflight",
                            "source_id": "OSD-X",
                            "contrast_id": "contrast_a",
                            "gene_symbol": "GENE2",
                            "ensembl_id": "ENSG2",
                            "predicted_log2fc_leo_or_iss_minus_ground": "0.5",
                        },
                        {
                            "task_id": "draft_human_organoid_spaceflight",
                            "source_id": "OSD-X",
                            "contrast_id": "contrast_a",
                            "gene_symbol": "GENE3",
                            "ensembl_id": "ENSG3",
                            "predicted_log2fc_leo_or_iss_minus_ground": "0.1",
                        },
                    ]
                )
            result = compute_response_signature_metrics(
                manifest=manifest,
                response_signature_path=response_path,
                reference_signature_path=reference_path,
            )

        self.assertTrue(result["validation"]["ok"])
        self.assertEqual(result["metrics"]["de_direction_match"]["status"], "computed")
        self.assertAlmostEqual(result["metrics"]["de_direction_match"]["value"], 0.5)
        self.assertEqual(
            result["metrics"]["de_direction_match"]["details"]["aggregate"][
                "n_direction_scored"
            ],
            2,
        )
        self.assertEqual(result["metrics"]["signature_rank_correlation"]["status"], "computed")
        self.assertAlmostEqual(result["metrics"]["signature_rank_correlation"]["value"], 0.5)

    def test_evaluator_computes_signature_metrics_when_response_artifact_supplied(self):
        from spacebio_bench import evaluate_submission

        manifest = {
            "task_id": "draft_human_organoid_spaceflight",
            "metrics": [
                {"metric_id": "de_direction_match"},
                {"metric_id": "signature_rank_correlation"},
            ],
            "output": {
                "required_columns": ["sample_id", "true_label", "predicted_label"],
                "label_domain": ["Ground", "LEO_or_ISS"],
                "positive_label": "LEO_or_ISS",
            },
        }
        prediction_rows = [
            {"sample_id": "s1", "true_label": "Ground", "predicted_label": "Ground"},
            {
                "sample_id": "s2",
                "true_label": "LEO_or_ISS",
                "predicted_label": "LEO_or_ISS",
            },
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            submission = tmp / "predictions.csv"
            self.write_submission(submission, prediction_rows)
            reference = tmp / "reference.csv"
            self.write_submission(
                reference,
                [
                    {
                        "source_id": "OSD-X",
                        "contrast_id": "contrast_a",
                        "gene_symbol": "GENE1",
                        "ensembl_id": "ENSG1",
                        "log2fc_leo_or_iss_minus_ground": "1.0",
                        "significant_fdr_0_05": "true",
                    },
                    {
                        "source_id": "OSD-X",
                        "contrast_id": "contrast_a",
                        "gene_symbol": "GENE2",
                        "ensembl_id": "ENSG2",
                        "log2fc_leo_or_iss_minus_ground": "-1.0",
                        "significant_fdr_0_05": "true",
                    },
                ],
            )
            response = tmp / "response_signature.csv"
            self.write_submission(
                response,
                [
                    {
                        "task_id": "draft_human_organoid_spaceflight",
                        "source_id": "OSD-X",
                        "contrast_id": "contrast_a",
                        "gene_symbol": "GENE1",
                        "ensembl_id": "ENSG1",
                        "predicted_log2fc_leo_or_iss_minus_ground": "1.0",
                    },
                    {
                        "task_id": "draft_human_organoid_spaceflight",
                        "source_id": "OSD-X",
                        "contrast_id": "contrast_a",
                        "gene_symbol": "GENE2",
                        "ensembl_id": "ENSG2",
                        "predicted_log2fc_leo_or_iss_minus_ground": "-1.0",
                    },
                ],
            )
            result = evaluate_submission(
                manifest,
                submission,
                response_signature_path=response,
                reference_signature_path=reference,
            )

        self.assertEqual(result["status"], "evaluated")
        self.assertEqual(result["response_signature_validation"]["ok"], True)
        self.assertEqual(result["metrics"]["de_direction_match"]["status"], "computed")
        self.assertAlmostEqual(result["metrics"]["de_direction_match"]["value"], 1.0)
        self.assertEqual(
            result["metrics"]["signature_rank_correlation"]["status"],
            "computed",
        )

    def test_response_signature_metrics_skip_invalid_artifact(self):
        from spacebio_bench import compute_response_signature_metrics

        manifest = {"task_id": "draft_human_organoid_spaceflight"}
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            reference = tmp / "reference.csv"
            self.write_submission(
                reference,
                [
                    {
                        "source_id": "OSD-X",
                        "contrast_id": "contrast_a",
                        "gene_symbol": "GENE1",
                        "ensembl_id": "ENSG1",
                        "log2fc_leo_or_iss_minus_ground": "1.0",
                        "significant_fdr_0_05": "true",
                    }
                ],
            )
            response = tmp / "response_signature.csv"
            response.write_text("task_id,source_id,contrast_id,gene_symbol,ensembl_id\n")
            result = compute_response_signature_metrics(
                manifest=manifest,
                response_signature_path=response,
                reference_signature_path=reference,
            )

        self.assertFalse(result["validation"]["ok"])
        self.assertEqual(result["metrics"]["de_direction_match"]["status"], "skipped")
        self.assertIn("missing required", result["metrics"]["de_direction_match"]["reason"])

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
        self.assertEqual(
            package["spacebio_bench:release_status"],
            "metadata_alpha_not_frozen",
        )
        self.assertEqual(
            package["spacebio_bench:alpha_snapshot_status"],
            "metadata_only_alpha_snapshot",
        )
        self.assertEqual(
            package["spacebio_bench:claim_boundary"],
            "metadata_only_public_bulk_alpha_no_payload_release",
        )
        self.assertFalse(package["spacebio_bench:payload_release_allowed"])
        self.assertIn("source_checksum_audit", resources)
        self.assertIn("baseline_predictions", resources)
        self.assertIn("public_bulk_alpha_gap_matrix", resources)
        self.assertIn("public_bulk_alpha_snapshot_claim_boundary", resources)
        self.assertEqual(resources["task_manifests"]["spacebio_bench:file_count"], 8)
        self.assertEqual(resources["baseline_predictions"]["spacebio_bench:file_count"], 24)
        self.assertEqual(
            resources["public_bulk_alpha_gap_matrix"]["spacebio_bench:bundle_part"],
            "alpha_boundary_metadata",
        )
        self.assertEqual(
            resources["public_bulk_alpha_snapshot_claim_boundary"]["schema"]["primaryKey"],
            ["snapshot_decision_id", "claim_id"],
        )
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
        self.assertEqual(len(payload["resources"]), 21)
        self.assertEqual(
            payload["spacebio_bench:claim_boundary"],
            "metadata_only_public_bulk_alpha_no_payload_release",
        )

    def test_generated_v9_datapackage_draft_records_metadata_alpha_boundary(self):
        payload = json.loads((REPO_ROOT / "v9" / "datapackage.draft.json").read_text())
        resources = {resource["name"]: resource for resource in payload["resources"]}

        self.assertEqual(len(resources), 21)
        self.assertEqual(
            payload["spacebio_bench:release_status"],
            "metadata_alpha_not_frozen",
        )
        self.assertEqual(
            payload["spacebio_bench:alpha_snapshot_status"],
            "metadata_only_alpha_snapshot",
        )
        self.assertFalse(payload["spacebio_bench:payload_release_allowed"])
        self.assertEqual(
            payload["spacebio_bench:payload_verification_status"],
            "checksum_manifests_parsed_payloads_not_hashed",
        )
        alpha_resources = {
            name
            for name, resource in resources.items()
            if resource.get("spacebio_bench:bundle_part") == "alpha_boundary_metadata"
        }
        self.assertEqual(
            alpha_resources,
            {
                "public_bulk_alpha_gap_summary",
                "public_bulk_alpha_gap_matrix",
                "public_bulk_alpha_payload_hash_boundary",
                "public_bulk_alpha_claim_boundary",
                "public_bulk_alpha_package_update_plan",
                "public_bulk_alpha_snapshot_decision_summary",
                "public_bulk_alpha_snapshot_option_matrix",
                "public_bulk_alpha_snapshot_claim_boundary",
                "public_bulk_alpha_snapshot_language_snippets",
                "public_bulk_alpha_snapshot_next_actions",
            },
        )
        self.assertNotIn("payload_mirror", resources)

    def test_v9_hf_dataset_card_records_draft_boundaries(self):
        card = (REPO_ROOT / "docs" / "v9_hf_dataset_card.md").read_text()

        required_phrases = [
            "Release status: metadata-only alpha snapshot, not frozen.",
            "metadata-only alpha snapshot",
            "not a frozen payload release",
            "does not include a local payload mirror",
            "payload-level hash verification",
            "22 deduplicated public OSDR source rows",
            "8 generated public bulk LOMO task manifests",
            "metadata_alpha_not_frozen",
            "metadata_only_public_bulk_alpha_no_payload_release",
            "payload_release_allowed = false",
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
