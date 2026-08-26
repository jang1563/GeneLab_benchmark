#!/usr/bin/env python3
"""Build the V9-SC-003 genelab_sc metric specification."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import DEFAULT_SC_METRIC_SPEC_ID, write_sc_metric_spec  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root used to resolve the single-cell manifest draft.",
    )
    parser.add_argument(
        "--output-dir",
        default="v9/sc_spaceflight",
        help="Directory for v9 single-cell metric-specification artifacts.",
    )
    parser.add_argument(
        "--metric-spec-id",
        default=DEFAULT_SC_METRIC_SPEC_ID,
        help="Identifier for the metric-specification artifacts.",
    )
    parser.add_argument(
        "--manifest-path",
        default=(
            "v9/sc_spaceflight/task_manifests/"
            "draft_rrrm1_blood_single_cell_spaceflight.json"
        ),
        help="Draft sc_spaceflight manifest used to anchor the metric spec.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = write_sc_metric_spec(
        repo_root=args.repo_root,
        output_dir=args.output_dir,
        metric_spec_id=args.metric_spec_id,
        manifest_path=args.manifest_path,
    )
    print(outputs["summary_csv"])
    print(outputs["summary_json"])
    print(outputs["metrics_csv"])
    print(outputs["metrics_json"])
    print(outputs["inputs_csv"])
    print(outputs["inputs_json"])
    print(outputs["skip_csv"])
    print(outputs["skip_json"])
    print(outputs["review_md"])


if __name__ == "__main__":
    main()
