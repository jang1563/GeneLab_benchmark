import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
UPLOAD_SCRIPT = REPO_ROOT / "scripts" / "upload_to_hf.py"


class HuggingFaceUploadWorkflowTests(unittest.TestCase):
    def run_upload_plan(self, *args: str) -> tuple[subprocess.CompletedProcess[str], dict]:
        with tempfile.TemporaryDirectory() as tmpdir:
            plan_path = Path(tmpdir) / "upload_plan.json"
            proc = subprocess.run(
                [
                    sys.executable,
                    str(UPLOAD_SCRIPT),
                    *args,
                    "--dry-run",
                    "--write-upload-plan",
                    str(plan_path),
                ],
                cwd=REPO_ROOT,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
            )
            payload = json.loads(plan_path.read_text()) if plan_path.exists() else {}
            return proc, payload

    def test_card_only_dry_run_plan_has_card_files(self):
        proc, payload = self.run_upload_plan("--card-only")

        self.assertEqual(proc.returncode, 0, proc.stderr)
        repo_paths = {row["repo_path"] for row in payload["files"]}
        self.assertIn("README.md", repo_paths)
        self.assertIn("assets/hf_benchmark_summary.png", repo_paths)
        self.assertNotIn("manifest.json", repo_paths)

    def test_manifest_only_dry_run_plan_has_release_metadata(self):
        proc, payload = self.run_upload_plan("--manifest-only")

        self.assertEqual(proc.returncode, 0, proc.stderr)
        repo_paths = {row["repo_path"] for row in payload["files"]}
        self.assertEqual(
            repo_paths,
            {"manifest.json", "metadata/release_manifest.schema.json"},
        )

    def test_task_dry_run_plan_includes_task_card_and_manifest(self):
        proc, payload = self.run_upload_plan("--task", "A5", "--skip-prune")

        self.assertEqual(proc.returncode, 0, proc.stderr)
        repo_paths = {row["repo_path"] for row in payload["files"]}
        self.assertIn("A5_skin_lomo/task_info.json", repo_paths)
        self.assertIn("A5_skin_lomo/fold_RR-7_test/fold_info.json", repo_paths)
        self.assertIn("README.md", repo_paths)
        self.assertIn("manifest.json", repo_paths)

        first = payload["files"][0]
        self.assertIn("sha256", first)
        self.assertIn("size_bytes", first)


if __name__ == "__main__":
    unittest.main()
