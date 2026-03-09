#!/usr/bin/env python3
"""
annotate_radiation.py — GeneLab_benchmark v2.0: T4 Radiation Metadata

Annotates metadata with estimated radiation dose based on mission duration.
ISS dose rate ≈ 200-230 µGy/day (nearly constant across LEO missions).

Mission-specific data (from NASA dosimetry):
  RR-1:  7.4 mGy / 37d → 0.200 mGy/d
  MHU-1: 8.05 mGy / 35d → 0.230 mGy/d
  General: ~0.215 mGy/d (average)

Output:
  v2/processed/T4_radiation/
    radiation_metadata.json

Usage:
  python v2/scripts/annotate_radiation.py
"""

import json
import sys
from pathlib import Path
from datetime import datetime

V2_DIR = Path(__file__).resolve().parent.parent
PROJECT_DIR = V2_DIR.parent
sys.path.insert(0, str(PROJECT_DIR / "scripts"))

from utils import load_temporal_metadata

# ── Mission dosimetry ──────────────────────────────────────────────────────────

# Duration in days for each mission (from OSDR metadata)
MISSION_DURATION = {
    "RR-1": 37,
    "RR-3": 30,    # ~30 days (SpaceX CRS-4, ISS 30d)
    "RR-5": 35,    # ~35 days
    "RR-6": 35,    # ~35 days (similar to MHU-1)
    "RR-7": 35,    # ~35 days
    "RR-8": 35,    # ~35 days (similar duration, includes LAR ~40d total)
    "RR-9": 35,    # ~35 days
    "MHU-1": 35,
    "MHU-2": 35,
    "TBD": 35,     # eye mission, ~35d estimate
}

# Dose rates (mGy/day) from NASA dosimetry
# Nearly constant for LEO ISS missions
DOSE_RATE_MGY_D = {
    "RR-1": 0.200,   # 7.4 mGy / 37d
    "MHU-1": 0.230,  # 8.05 mGy / 35d
}
DEFAULT_DOSE_RATE = 0.215  # average

# GC/BSL/VIV labels (ground controls, no radiation)
GROUND_LABELS = {"GC", "BC", "BSL", "VC", "VIV"}


def annotate_radiation():
    """Generate radiation metadata for all missions."""
    print("=" * 60)
    print("T4: Radiation Metadata Annotation")
    print("  ISS LEO dose rate ≈ 0.200-0.230 mGy/day")
    print("=" * 60)

    output_dir = V2_DIR / "processed" / "T4_radiation"
    output_dir.mkdir(parents=True, exist_ok=True)

    radiation_data = {
        "generated_at": datetime.now().isoformat(),
        "task": "T4",
        "description": "Radiation dose annotation based on ISS LEO dosimetry",
        "sources": [
            "NASA Environmental Data for ISS",
            "RR-1: 7.4 mGy total / 37d flight",
            "MHU-1: 8.05 mGy total / 35d flight",
        ],
        "note": "ISS LEO dose rate is nearly constant across missions (~0.215 mGy/d). "
                "Dose-response analysis not meaningful with current data (requires variable doses). "
                "This annotation is for metadata completeness and future use.",
        "missions": {},
    }

    for mission, duration in sorted(MISSION_DURATION.items()):
        rate = DOSE_RATE_MGY_D.get(mission, DEFAULT_DOSE_RATE)
        total_dose = rate * duration

        radiation_data["missions"][mission] = {
            "duration_days": duration,
            "dose_rate_mgy_per_day": rate,
            "total_dose_mgy": round(total_dose, 2),
            "dose_source": "measured" if mission in DOSE_RATE_MGY_D else "estimated",
            "ground_control_dose_mgy": 0.0,
            "note": f"Flight samples: {total_dose:.2f} mGy, GC/BSL/VIV: 0 mGy"
        }

        print(f"  {mission}: {duration}d × {rate:.3f} mGy/d = {total_dose:.2f} mGy "
              f"({'measured' if mission in DOSE_RATE_MGY_D else 'estimated'})")

    # Summary
    doses = [v["total_dose_mgy"] for v in radiation_data["missions"].values()]
    radiation_data["summary"] = {
        "min_dose_mgy": min(doses),
        "max_dose_mgy": max(doses),
        "dose_range_mgy": round(max(doses) - min(doses), 2),
        "conclusion": "Dose range too narrow for dose-response analysis "
                      "(6.0-8.1 mGy, all LEO). Ground-based variable-dose experiments needed."
    }

    outpath = output_dir / "radiation_metadata.json"
    outpath.write_text(json.dumps(radiation_data, indent=2))
    print(f"\n  Saved: {outpath}")
    print(f"\n  Conclusion: All missions 6-8 mGy range (LEO). "
          f"Dose-response not feasible with current data.")


if __name__ == "__main__":
    annotate_radiation()
