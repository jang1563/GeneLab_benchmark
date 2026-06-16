import json
import subprocess
import sys
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


class ReleaseManifestTests(unittest.TestCase):
    def test_release_manifest_validator_passes_current_manifest(self):
        proc = subprocess.run(
            [
                sys.executable,
                str(REPO_ROOT / "scripts" / "validate_release_manifest.py"),
                "--manifest",
                str(REPO_ROOT / "release" / "release_manifest.json"),
            ],
            cwd=REPO_ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )

        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("[OK]", proc.stdout)

    def test_release_manifest_declares_public_release_lanes(self):
        manifest = json.loads((REPO_ROOT / "release" / "release_manifest.json").read_text())
        lanes = {lane["lane_id"]: lane for lane in manifest["release_lanes"]}

        self.assertEqual(
            lanes["v7.1"]["status"],
            "canonical_historical_result_surface",
        )
        self.assertEqual(
            lanes["v9_public_bulk"]["status"],
            "metadata_only_public_bulk_alpha",
        )
        self.assertIn("docs/hf_dataset_card.md", lanes["v7.1"]["canonical_surfaces"])
        self.assertIn("docs/v9_hf_dataset_card.md", lanes["v9_public_bulk"]["canonical_surfaces"])

    def test_release_manifest_artifacts_reference_known_lanes(self):
        manifest = json.loads((REPO_ROOT / "release" / "release_manifest.json").read_text())
        lane_ids = {lane["lane_id"] for lane in manifest["release_lanes"]}

        for artifact in manifest["artifacts"]:
            self.assertIn(artifact["lane_id"], lane_ids)
            artifact_path = REPO_ROOT / artifact["path"]
            self.assertTrue(artifact_path.exists(), artifact["path"])


if __name__ == "__main__":
    unittest.main()
