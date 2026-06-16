import subprocess
import sys
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


class PublicDocsConsistencyTests(unittest.TestCase):
    def test_public_docs_consistency_script_passes(self):
        proc = subprocess.run(
            [
                sys.executable,
                str(REPO_ROOT / "scripts" / "validate_public_docs_consistency.py"),
            ],
            cwd=REPO_ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )

        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("[OK]", proc.stdout)


if __name__ == "__main__":
    unittest.main()
