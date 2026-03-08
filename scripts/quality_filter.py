#!/usr/bin/env python3
"""
quality_filter.py — GeneLab_benchmark: RNA-seq Quality Filtering

Applies sample-level and mission-level QC filters to downloaded normalized counts.
Outputs clean data to processed/ directory and a QC report.

QC criteria (from PLAN.md v0.5 Section 7):

  Sample-level (hard cutoffs):
    - Total counts (library size) >= 1,000,000
    - Detected genes (count > 0) >= 5,000
    - PCA outlier <= 4 SD from centroid (computed per tissue × condition)
    - Pearson r with group median <= 0.99 per group (same-group correlation check)

  Mission-level (eligibility gates):
    - n_samples >= 3 per group (flight and ground control)
    - At least one ground control type present (GC or VC)

  Normalization for Category A input (DESIGN_DECISIONS DD-01):
    - Input: DESeq2 normalized counts (downloaded from GeneLab pipeline)
    - Output: log2(normalized counts + 1) per sample

  Feature selection (DESIGN_DECISIONS DD-03):
    - Low-expression filter (applied globally, once):
        Keep genes present (count > 1) in >= 20% of ALL samples within tissue
    - Variance filter: applied INSIDE LOMO loop in generate_tasks.py (NOT here)

Usage:
  python quality_filter.py --tissue liver         # process liver
  python quality_filter.py --all                  # process all tissues
  python quality_filter.py --tissue liver --report # generate QC report only
  python quality_filter.py --check-input liver     # verify input files exist

Output structure:
  processed/A_detection/{tissue}/
    {tissue}_normalized_log2.csv       # log2(norm+1), samples × genes
    {tissue}_sample_metadata.csv       # sample labels + mission + QC flags
    {tissue}_qc_report.json            # per-sample QC metrics
    {tissue}_removed_samples.csv       # removed samples with reason

Design decisions applied:
  - DD-01: Category A uses log2(normalized counts), NOT LFC
  - DD-03: Variance filter NOT applied here (must be inside LOMO loop)
  - DD-04: Mission integrity preserved (no cross-mission sample mixing)
  - DD-05: GLDS-168 excluded from Category A output
  - DD-06: Strain tracked (C57BL/6J = Track 2a, others = Track 2b)
"""

import json
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from scipy.stats import pearsonr

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data" / "mouse"
PROCESSED_DIR = BASE_DIR / "processed" / "A_detection"
QC_DIR = BASE_DIR / "processed" / "qc_reports"
GLDS_VERIFIED_JSON = BASE_DIR / "GLDS_verified.json"

# ── QC thresholds (PLAN.md v0.5, Section 7.1) ─────────────────────────────────
QC_MIN_LIBRARY_SIZE = 1_000_000    # total counts
QC_MIN_DETECTED_GENES = 5_000      # genes with count > 0
QC_PCA_SD_CUTOFF = 4.0             # SD from centroid
QC_GROUP_CORR_CUTOFF = 0.9999      # Pairwise Pearson r — true technical dups only
QC_MIN_EXPR_PCT = 0.20             # low-expression filter: gene present in ≥20% samples
QC_MIN_SAMPLES_PER_GROUP = 3       # minimum biological replicates

# ── Mission metadata ───────────────────────────────────────────────────────────
# Maps OSD → mission info. Used to label processed outputs.
OSD_TO_MISSION = {
    "OSD-48":  {"mission": "RR-1",  "duration_days": 37,  "strain": "C57BL/6J", "track": "2a"},
    "OSD-137": {"mission": "RR-3",  "duration_days": 40,  "strain": "C57BL/6J", "track": "2a"},
    "OSD-245": {"mission": "RR-6",  "duration_days": 35,  "strain": "C57BL/6J", "track": "2a"},
    "OSD-379": {"mission": "RR-8",  "duration_days": 30,  "strain": "C57BL/6J", "track": "2a"},
    "OSD-242": {"mission": "RR-9",  "duration_days": 33,  "strain": "C57BL/6J", "track": "2a"},
    "OSD-686": {"mission": "MHU-2", "duration_days": 35,  "strain": "C57BL/6J", "track": "2a",
                "note": "3 groups: uG, 1G-centrifuge, GC. Use uG vs GC for A1."},
    "OSD-101": {"mission": "RR-1",  "duration_days": 37,  "strain": "C57BL/6J", "track": "2a"},
    "OSD-401": {"mission": "RR-5",  "duration_days": 35,  "strain": "C57BL/6J", "track": "2a"},
    "OSD-326": {"mission": "RR-9",  "duration_days": 33,  "strain": "C57BL/6J", "track": "2a"},
    "OSD-102": {"mission": "RR-1",  "duration_days": 37,  "strain": "C57BL/6J", "track": "2a"},
    "OSD-163": {"mission": "RR-3",  "duration_days": 40,  "strain": "C57BL/6J", "track": "2a"},
    "OSD-253": {"mission": "RR-7",  "duration_days": 75,  "strain": "C57BL/6J", "track": "2a"},
    "OSD-244": {"mission": "RR-6",  "duration_days": 35,  "strain": "C57BL/6J", "track": "2a"},
    "OSD-289": {"mission": "MHU-2", "duration_days": 35,  "strain": "C57BL/6J", "track": "2a"},
    "OSD-421": {"mission": "RR-9",  "duration_days": 33,  "strain": "C57BL/6J", "track": "2a"},
    "OSD-100": {"mission": "RR-1",  "duration_days": 37,  "strain": "C57BL/6J", "track": "2a"},
    "OSD-194": {"mission": "RR-3",  "duration_days": 40,  "strain": "C57BL/6J", "track": "2a"},
    "OSD-397": {"mission": "TBD",   "duration_days": None, "strain": "Mus musculus", "track": "2a"},
    "OSD-238": {"mission": "MHU-2 (dorsal)", "duration_days": 35, "strain": "C57BL/6J", "track": "2a"},
    "OSD-239": {"mission": "MHU-2 (femoral)","duration_days": 35, "strain": "C57BL/6J", "track": "2a"},
    "OSD-243": {"mission": "RR-6",  "duration_days": 35,  "strain": "C57BL/6J", "track": "2a"},
    "OSD-254": {"mission": "RR-7",  "duration_days": 75,  "strain": "C57BL/6J", "track": "2a",
                "note": "Mixed strain (C57BL/6J + C3H/HeJ). C57BL/6J subset only for Track 2a."},
    "OSD-295": {"mission": "HU",    "duration_days": None, "strain": "?", "track": "HU"},
}

# ── Tissue → OSD mapping ───────────────────────────────────────────────────────
TISSUE_OSD_MAP = {
    "liver":          ["OSD-48", "OSD-137", "OSD-245", "OSD-379", "OSD-242", "OSD-686"],
    "gastrocnemius":  ["OSD-101", "OSD-401", "OSD-326"],
    "kidney":         ["OSD-102", "OSD-163", "OSD-253"],
    "thymus":         ["OSD-244", "OSD-289", "OSD-421"],
    "eye":            ["OSD-100", "OSD-194", "OSD-397"],
    "skin":           ["OSD-238", "OSD-239", "OSD-243", "OSD-254"],
    "soleus_HU":      ["OSD-295"],
}


# ── File loading ───────────────────────────────────────────────────────────────

def find_normalized_counts_file(osd_dir: Path, glds_prefix: str) -> Path | None:
    """
    Find normalized counts CSV for a study.
    Tries GeneLab v2 naming first, falls back to v1.
    """
    candidates = [
        f"{glds_prefix}_rna_seq_Normalized_Counts_GLbulkRNAseq.csv",
        f"{glds_prefix}_rna_seq_Normalized_Counts_rRNArm_GLbulkRNAseq.csv",
        f"{glds_prefix}_rna_seq_Normalized_Counts.csv",
        f"{glds_prefix}_rna_seq_ERCC_Normalized_Counts.csv",
    ]
    for cname in candidates:
        fpath = osd_dir / cname
        if fpath.exists():
            return fpath
    return None


def find_sample_table_file(osd_dir: Path, glds_prefix: str) -> Path | None:
    """Find sample table (metadata) CSV."""
    candidates = [
        f"{glds_prefix}_rna_seq_SampleTable_GLbulkRNAseq.csv",
        f"{glds_prefix}_rna_seq_SampleTable.csv",
        f"{glds_prefix}_rna_seq_bulkRNASeq_v2_runsheet.csv",
        f"{glds_prefix}_rna_seq_bulkRNASeq_v1_runsheet.csv",
    ]
    for cname in candidates:
        fpath = osd_dir / cname
        if fpath.exists():
            return fpath
    return None


def load_normalized_counts(filepath: Path) -> pd.DataFrame:
    """
    Load GeneLab normalized counts CSV.
    Returns DataFrame: rows=genes, columns=samples.
    Handles both row-name and no-row-name formats.
    """
    df = pd.read_csv(filepath, index_col=0)
    # Drop non-numeric columns if any
    df = df.select_dtypes(include=[np.number])
    return df


def load_sample_table(filepath: Path) -> pd.DataFrame:
    """Load sample table / runsheet CSV."""
    return pd.read_csv(filepath, index_col=0)


def infer_sample_labels(sample_table: pd.DataFrame,
                        osd_id: str) -> dict[str, str]:
    """
    Infer Flight/Ground labels from sample table.
    Returns {sample_name: label} where label ∈ {'Flight', 'GC', 'VC', 'BC', 'AG'}.

    Handles multiple GeneLab naming conventions.
    """
    labels = {}

    # Try 'condition' column first (standard SampleTable format)
    condition_col = None
    for col in ["condition", "Condition", "Factor Value[Spaceflight]",
                "Factor.Value.Spaceflight.", "group"]:
        if col in sample_table.columns:
            condition_col = col
            break

    if condition_col is None:
        print(f"    [WARN] No condition column found for {osd_id}. Columns: {list(sample_table.columns)}")
        return {}

    # IMPORTANT: more-specific patterns must come BEFORE less-specific ones
    # because the first match wins. E.g., "1g by centrifugation" must precede
    # "flight" so that "Space.Flight...1G.by.centrifugation" → AG, not Flight.
    keyword_map = {
        # ── Specific patterns first ──────────────────────────────────────────
        "1g by centrifugation": "AG",  # OSD-686 MHU-2 artificial gravity (before "flight")
        "1g on earth": "GC",           # OSD-686 MHU-2 GC (before "gc")
        "ground control": "GC",
        "ground.control": "GC",
        "vivarium control": "VC",
        "vivarium.control": "VC",
        # ── Generic patterns ─────────────────────────────────────────────────
        "flight": "Flight",
        "spaceflight": "Flight",
        "flt": "Flight",
        "microgravity": "Flight",
        "ug": "Flight",              # microgravity (μg) shorthand
        "gc": "GC",
        "cc": "GC",                  # Cage Control (OSD-242 RR-9)
        "vivarium": "VC",
        "viv": "VC",
        "basal": "BC",
        "baseline": "BC",
        "bsl": "BC",                 # Basal shorthand (OSD-242 RR-9)
    }

    for sample, row in sample_table.iterrows():
        cond_raw = str(row[condition_col]).lower().strip().replace(".", " ")
        label = "Unknown"
        for key, lbl in keyword_map.items():
            if key in cond_raw:
                label = lbl
                break
        labels[str(sample)] = label

    return labels


# ── Temporal metadata (v2.0) ─────────────────────────────────────────────────

def infer_sacrifice_timing(sample_name: str, mission: str) -> str:
    """
    Infer sacrifice timing from sample name.
    Returns: 'ISS-T' | 'LAR' | 'unknown'

    ISS-T (ISS Terminal) = sacrificed on-orbit, preserved with RNAlater
    LAR (Live Animal Return) = returned alive, standard necropsy on ground

    Parsing rules verified against OSDR runsheets:
      - RR-6, RR-8:  '_ISS-T_' or '_LAR_' in sample name
      - RR-1 liver:   '_FLT_C_' / '_GC_C_' → ISS-T (Carcass = RNAlater)
                       '_FLT_I_' / '_GC_I_' → LAR (Euthanasia = live return)
                       Verified via GLDS-48 runsheet 'Dissection Condition' column
      - MHU-2:        all LAR (JAXA protocol, mice returned live)
      - RR-3, RR-9:   no temporal split → 'unknown'
    """
    name = sample_name.upper()

    # RR-6, RR-8: explicit ISS-T / LAR in name
    if '_ISS-T_' in name or '_ISS-T' == name[-5:]:
        return 'ISS-T'
    if '_LAR_' in name or '_LAR' == name[-4:]:
        return 'LAR'

    # RR-1 liver: _C_ = Carcass (ISS-T), _I_ = Euthanasia (LAR)
    # Pattern: _{FLT|GC|BSL}_{C|I}_ — must follow group label
    mission_upper = mission.upper().replace(' ', '')
    if mission_upper in ('RR-1',):
        import re
        # Match _FLT_C_, _GC_C_, _BSL_C_ etc. (Carcass → ISS-T)
        if re.search(r'_(FLT|GC|BSL|VIV|VC)_C_', name):
            return 'ISS-T'
        # Match _FLT_I_, _GC_I_, etc. (Euthanasia → LAR)
        if re.search(r'_(FLT|GC|BSL|VIV|VC)_I_', name):
            return 'LAR'

    # MHU-2: all returned live (JAXA protocol)
    if 'MHU' in mission_upper:
        return 'LAR'

    return 'unknown'


def infer_age_group(sample_name: str) -> str:
    """
    Infer age group from sample name (RR-8 only).
    Returns: 'OLD' | 'YNG' | 'unknown'

    RR-8 mice: OLD = 32-week, YNG = 10-12 week
    Pattern: '_OLD_' or '_YNG_' in sample name
    """
    name = sample_name.upper()
    if '_OLD_' in name or name.endswith('_OLD'):
        return 'OLD'
    if '_YNG_' in name or name.endswith('_YNG'):
        return 'YNG'
    return 'unknown'


def enrich_temporal_metadata(tissue: str, verbose: bool = True) -> None:
    """
    Post-hoc enrichment: add sacrifice_timing and age_group columns
    to existing metadata CSVs in processed/A_detection/{tissue}/.
    Does NOT re-run QC — only reads and updates metadata files.
    """
    outdir = PROCESSED_DIR / tissue

    # Find all metadata CSVs
    meta_files = sorted(outdir.glob(f"{tissue}_*_metadata.csv"))
    if not meta_files:
        print(f"  [SKIP] No metadata files found in {outdir}")
        return

    print(f"\n  Enriching temporal metadata for {tissue}...")

    for meta_path in meta_files:
        meta = pd.read_csv(meta_path, index_col=0)

        # Determine mission from metadata
        if 'mission' in meta.columns:
            missions = meta['mission'].unique()
        else:
            missions = ['unknown']

        # Add sacrifice_timing
        timing_values = []
        age_values = []
        for sample_name in meta.index:
            mission = meta.loc[sample_name, 'mission'] if 'mission' in meta.columns else 'unknown'
            timing_values.append(infer_sacrifice_timing(str(sample_name), str(mission)))
            age_values.append(infer_age_group(str(sample_name)))

        meta['sacrifice_timing'] = timing_values
        meta['age_group'] = age_values

        # Report
        timing_counts = pd.Series(timing_values).value_counts().to_dict()
        age_counts = pd.Series(age_values).value_counts().to_dict()

        if verbose:
            mission_str = ', '.join(str(m) for m in missions)
            print(f"    {meta_path.name}: {mission_str}")
            print(f"      sacrifice_timing: {timing_counts}")
            if any(v != 'unknown' for v in age_values):
                print(f"      age_group: {age_counts}")

        # Save
        meta.to_csv(meta_path)

    print(f"  Done: {len(meta_files)} metadata files enriched.")


# ── QC functions ───────────────────────────────────────────────────────────────

def compute_sample_qc_metrics(counts: pd.DataFrame) -> pd.DataFrame:
    """
    Compute per-sample QC metrics.
    counts: rows=genes, columns=samples
    Returns DataFrame with QC metrics per sample.
    """
    metrics = pd.DataFrame(index=counts.columns)
    metrics["library_size"] = counts.sum(axis=0)
    metrics["n_detected_genes"] = (counts > 0).sum(axis=0)
    metrics["pct_top100_genes"] = (
        counts.apply(lambda col: col.nlargest(100).sum()) / metrics["library_size"] * 100
    )
    return metrics


def pca_outlier_detection(log2_counts: pd.DataFrame,
                          sd_cutoff: float = QC_PCA_SD_CUTOFF) -> pd.Series:
    """
    Detect PCA outliers: samples > sd_cutoff SDs from centroid on PC1+PC2.
    log2_counts: rows=samples, columns=genes
    Returns boolean Series (True = outlier).
    """
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA

    if log2_counts.shape[0] < 3:
        return pd.Series(False, index=log2_counts.index)

    X = log2_counts.values
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    n_components = min(2, X.shape[0] - 1, X.shape[1] - 1)
    pca = PCA(n_components=n_components)
    coords = pca.fit_transform(X_scaled)

    centroid = coords.mean(axis=0)
    dists = np.linalg.norm(coords - centroid, axis=1)
    mean_d, std_d = dists.mean(), dists.std()
    outlier = (dists - mean_d) > (sd_cutoff * std_d)

    return pd.Series(outlier, index=log2_counts.index)


def within_group_correlation_check(log2_counts: pd.DataFrame,
                                   labels: pd.Series,
                                   corr_cutoff: float = QC_GROUP_CORR_CUTOFF) -> pd.Series:
    """
    Flag samples that are near-exact duplicates of another sample (technical duplicates).
    Uses PAIRWISE correlation between samples (NOT correlation to group median).
    Only removes one sample from each pair with r > corr_cutoff (default 0.9999).

    Biological replicates typically show r ~ 0.97-0.99. True technical duplicates
    (same sample sequenced twice) show r ~ 0.9999-1.000.
    """
    suspicious = pd.Series(False, index=log2_counts.index)
    already_flagged = set()

    for group in labels.unique():
        members = labels[labels == group].index
        if len(members) < 2:
            continue
        group_data = log2_counts.loc[members]

        for i, s1 in enumerate(members):
            if s1 in already_flagged:
                continue
            for s2 in members[i+1:]:
                if s2 in already_flagged:
                    continue
                v1, v2 = group_data.loc[s1].values, group_data.loc[s2].values
                if np.std(v1) < 1e-10 or np.std(v2) < 1e-10:
                    continue
                r, _ = pearsonr(v1, v2)
                if r >= corr_cutoff:
                    # Flag the second of the pair as duplicate
                    suspicious[s2] = True
                    already_flagged.add(s2)

    return suspicious


def low_expression_filter(counts: pd.DataFrame,
                           min_pct: float = QC_MIN_EXPR_PCT) -> pd.Index:
    """
    Remove genes present (count > 1) in fewer than min_pct of samples.
    Applied GLOBALLY within a tissue (across all missions).
    Returns index of genes to KEEP.

    NOTE: Variance filter is NOT applied here.
    See DESIGN_DECISIONS DD-03 — must be inside LOMO loop in generate_tasks.py.
    """
    n_samples = counts.shape[1]
    min_samples = int(np.ceil(min_pct * n_samples))
    expressed = (counts > 1).sum(axis=1) >= min_samples
    return counts.index[expressed]


# ── Per-study processing ───────────────────────────────────────────────────────

def process_study(osd_id: str, tissue: str, verbose: bool = True) -> dict | None:
    """
    Process a single study: load → QC → log2 transform.
    Returns dict with processed data and QC results, or None on failure.
    """
    # Find data directory
    mission_info = OSD_TO_MISSION.get(osd_id, {})
    mission = mission_info.get("mission", "unknown")
    mission_clean = mission.replace(" ", "_").replace("/", "_").replace("+", "_")
    osd_dir = DATA_DIR / tissue / mission_clean
    glds_prefix = f"GLDS-{osd_id.replace('OSD-', '')}"

    # Special case: OSD-686 uses GLDS-617 prefix
    if osd_id == "OSD-686":
        glds_prefix = "GLDS-617"

    if verbose:
        print(f"\n  [{osd_id}] {tissue} / {mission}")

    # Load files
    counts_file = find_normalized_counts_file(osd_dir, glds_prefix)
    if counts_file is None:
        print(f"    [SKIP] No normalized counts file found in {osd_dir}")
        print(f"    Run: python fetch_osdr.py --osd {osd_id}")
        return None

    meta_file = find_sample_table_file(osd_dir, glds_prefix)

    counts = load_normalized_counts(counts_file)
    if verbose:
        print(f"    Loaded: {counts.shape[1]} samples × {counts.shape[0]} genes")

    # Load sample metadata
    sample_labels = {}
    if meta_file:
        meta = load_sample_table(meta_file)
        sample_labels = infer_sample_labels(meta, osd_id)
    else:
        print(f"    [WARN] No sample table found. Labels will be Unknown.")

    # ── Sample-level QC ───────────────────────────────────────────────────────
    qc_metrics = compute_sample_qc_metrics(counts)

    # Flag samples by criteria
    low_library = qc_metrics["library_size"] < QC_MIN_LIBRARY_SIZE
    low_genes = qc_metrics["n_detected_genes"] < QC_MIN_DETECTED_GENES

    # Log2 transform for PCA and correlation checks
    log2_counts = np.log2(counts + 1).T  # → rows=samples, cols=genes

    pca_outliers = pca_outlier_detection(log2_counts)

    # Only check within-group correlation if we have labels
    labels_series = pd.Series(
        {s: sample_labels.get(s, "Unknown") for s in log2_counts.index}
    )
    dup_flags = within_group_correlation_check(log2_counts, labels_series)

    # Combine QC flags
    remove_mask = low_library | low_genes | pca_outliers | dup_flags
    qc_metrics["remove_low_library"] = low_library
    qc_metrics["remove_low_genes"] = low_genes
    qc_metrics["remove_pca_outlier"] = pca_outliers
    qc_metrics["remove_duplicate"] = dup_flags
    qc_metrics["REMOVE"] = remove_mask
    qc_metrics["label"] = labels_series
    qc_metrics["mission"] = mission
    qc_metrics["osd_id"] = osd_id
    qc_metrics["strain"] = mission_info.get("strain", "?")
    qc_metrics["track"] = mission_info.get("track", "?")
    qc_metrics["duration_days"] = mission_info.get("duration_days", None)

    n_removed = remove_mask.sum()
    n_kept = (~remove_mask).sum()

    if verbose:
        print(f"    QC: {n_kept} kept, {n_removed} removed "
              f"(low_lib={low_library.sum()}, low_genes={low_genes.sum()}, "
              f"pca_out={pca_outliers.sum()}, dup={dup_flags.sum()})")

    # Filter counts
    counts_clean = counts.loc[:, ~remove_mask]

    # Check label distribution
    label_counts = labels_series[~remove_mask].value_counts()
    if verbose:
        print(f"    Groups: {dict(label_counts)}")

    # Mission-level eligibility check
    flight_n = label_counts.get("Flight", 0)
    gc_n = label_counts.get("GC", 0)
    vc_n = label_counts.get("VC", 0)
    eligible = (flight_n >= QC_MIN_SAMPLES_PER_GROUP and
                (gc_n >= QC_MIN_SAMPLES_PER_GROUP or vc_n >= QC_MIN_SAMPLES_PER_GROUP))

    if not eligible:
        print(f"    [WARN] Mission {mission} fails eligibility: "
              f"Flight={flight_n}, GC={gc_n}, VC={vc_n} (need ≥{QC_MIN_SAMPLES_PER_GROUP}/group)")

    return {
        "osd_id": osd_id,
        "tissue": tissue,
        "mission": mission,
        "counts_clean": counts_clean,
        "log2_counts_clean": np.log2(counts_clean + 1),
        "qc_metrics": qc_metrics,
        "label_series": labels_series[~remove_mask],
        "mission_info": mission_info,
        "eligible": eligible,
        "n_samples_kept": n_kept,
        "n_samples_removed": n_removed,
        "label_distribution": dict(label_counts),
    }


# ── Tissue-level aggregation ───────────────────────────────────────────────────

def process_tissue(tissue: str, verbose: bool = True) -> None:
    """
    Process all studies for a tissue:
    1. Load and QC each study
    2. Apply global low-expression filter across all missions
    3. Save clean normalized counts + metadata
    """
    osd_ids = TISSUE_OSD_MAP.get(tissue)
    if not osd_ids:
        print(f"[ERROR] Unknown tissue: {tissue}. Known: {list(TISSUE_OSD_MAP.keys())}")
        return

    print(f"\n{'='*60}")
    print(f"Processing: {tissue.upper()} ({len(osd_ids)} studies)")
    print(f"{'='*60}")

    # Process each study
    study_results = []
    for osd_id in osd_ids:
        result = process_study(osd_id, tissue, verbose=verbose)
        if result is not None:
            study_results.append(result)

    if not study_results:
        print(f"  [SKIP] No data loaded for {tissue}. Download first.")
        return

    # ── Global low-expression filter (across all missions) ────────────────────
    # Collect all samples across missions into one matrix
    all_counts = pd.concat(
        [r["counts_clean"] for r in study_results], axis=1
    ).fillna(0)

    kept_genes = low_expression_filter(all_counts, min_pct=QC_MIN_EXPR_PCT)
    n_total = len(all_counts)
    n_kept = len(kept_genes)
    print(f"\n  Low-expression filter: {n_kept} / {n_total} genes kept "
          f"(≥{QC_MIN_EXPR_PCT:.0%} samples with count > 1)")
    print(f"  ⚠️  Variance filter NOT applied here — must be inside LOMO loop")

    # ── Save per-study processed files ────────────────────────────────────────
    outdir = PROCESSED_DIR / tissue
    outdir.mkdir(parents=True, exist_ok=True)
    qc_dir = QC_DIR / tissue
    qc_dir.mkdir(parents=True, exist_ok=True)

    all_log2_frames = []
    all_meta_frames = []
    qc_summary = []

    for r in study_results:
        osd_id = r["osd_id"]
        mission = r["mission"]
        mission_clean = mission.replace(" ", "_").replace("/", "_").replace("+", "_")

        # Apply global gene filter (intersect with genes present in this study)
        study_genes = r["log2_counts_clean"].index
        study_kept = kept_genes.intersection(study_genes)
        log2_filtered = r["log2_counts_clean"].loc[study_kept]

        # Save per-mission files
        mission_log2_path = outdir / f"{tissue}_{mission_clean}_log2_norm.csv"
        mission_meta_path = outdir / f"{tissue}_{mission_clean}_metadata.csv"
        log2_filtered.to_csv(mission_log2_path)

        # Build metadata frame
        meta_frame = r["qc_metrics"][~r["qc_metrics"]["REMOVE"]].copy()
        meta_frame["tissue"] = tissue
        meta_frame.to_csv(mission_meta_path)

        # Collect for tissue-wide concatenation
        log2_t = log2_filtered.T  # rows=samples
        log2_t["mission"] = mission
        log2_t["osd_id"] = osd_id
        all_log2_frames.append(log2_t)
        all_meta_frames.append(meta_frame)

        qc_summary.append({
            "osd_id": osd_id,
            "mission": mission,
            "tissue": tissue,
            "n_samples_kept": r["n_samples_kept"],
            "n_samples_removed": r["n_samples_removed"],
            "label_distribution": r["label_distribution"],
            "eligible_for_task": r["eligible"],
            "strain": r["mission_info"].get("strain", "?"),
            "track": r["mission_info"].get("track", "?"),
            "duration_days": r["mission_info"].get("duration_days"),
        })

    # ── Save tissue-wide combined files ───────────────────────────────────────
    if all_log2_frames:
        combined_log2 = pd.concat(all_log2_frames)
        # Fill NaN with 0: genes absent in a mission → log2(0+1) = 0 (not detected)
        # NaN arises when missions have different gene sets after pipeline-version differences
        gene_cols = [c for c in combined_log2.columns if c not in {"mission", "osd_id", "label"}]
        combined_log2[gene_cols] = combined_log2[gene_cols].fillna(0.0)
        combined_meta = pd.concat(all_meta_frames)

        combined_log2.to_csv(outdir / f"{tissue}_all_missions_log2_norm.csv")
        combined_meta.to_csv(outdir / f"{tissue}_all_missions_metadata.csv")

        eligible_missions = [q["mission"] for q in qc_summary if q["eligible_for_task"]]
        n_total_samples = combined_log2.shape[0] - len(all_log2_frames)  # subtract mission/osd cols

        print(f"\n  ✓ Tissue-wide file: {combined_log2.shape[0]} samples × {len(kept_genes)} genes")
        print(f"  ✓ Eligible missions: {eligible_missions}")
        print(f"  ✓ Saved to: {outdir}")

    # ── Save QC report ────────────────────────────────────────────────────────
    qc_report = {
        "generated_at": datetime.now().isoformat(),
        "tissue": tissue,
        "plan_version": "0.5",
        "n_genes_after_low_expr_filter": int(n_kept),
        "n_genes_before_filter": int(n_total),
        "low_expr_threshold": f"count > 1 in >= {QC_MIN_EXPR_PCT:.0%} of samples",
        "variance_filter": "NOT applied here — must be inside LOMO loop (DD-03)",
        "qc_thresholds": {
            "min_library_size": QC_MIN_LIBRARY_SIZE,
            "min_detected_genes": QC_MIN_DETECTED_GENES,
            "pca_sd_cutoff": QC_PCA_SD_CUTOFF,
            "within_group_corr_cutoff": QC_GROUP_CORR_CUTOFF,
            "min_samples_per_group": QC_MIN_SAMPLES_PER_GROUP,
        },
        "studies": qc_summary,
    }

    qc_path = qc_dir / f"{tissue}_qc_report.json"

    class _NumpyEncoder(json.JSONEncoder):
        def default(self, obj):
            import numpy as np
            if isinstance(obj, np.bool_):
                return bool(obj)
            if isinstance(obj, np.integer):
                return int(obj)
            if isinstance(obj, np.floating):
                return float(obj)
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            return super().default(obj)

    qc_path.write_text(json.dumps(qc_report, indent=2, cls=_NumpyEncoder))
    print(f"  ✓ QC report: {qc_path}")


# ── CLI ────────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="Apply QC filters to downloaded RNA-seq normalized counts"
    )
    parser.add_argument(
        "--tissue", type=str, default=None,
        help=f"Tissue to process. Options: {list(TISSUE_OSD_MAP.keys())}"
    )
    parser.add_argument(
        "--all", action="store_true",
        help="Process all tissues"
    )
    parser.add_argument(
        "--check-input", type=str, default=None, metavar="TISSUE",
        help="Check if input files exist for a tissue (no processing)"
    )
    parser.add_argument(
        "--report", action="store_true",
        help="Print summary of what QC would do (no output written)"
    )
    parser.add_argument(
        "--enrich-temporal", action="store_true",
        help="Add sacrifice_timing and age_group columns to existing metadata CSVs (v2.0)"
    )
    parser.add_argument(
        "--verbose", action="store_true", default=True,
        help="Verbose output (default: True)"
    )
    return parser.parse_args()


def check_input_files(tissue: str) -> None:
    """Check which input files exist for a tissue."""
    osd_ids = TISSUE_OSD_MAP.get(tissue, [])
    print(f"\n=== Input file check: {tissue} ===")
    for osd_id in osd_ids:
        mission_info = OSD_TO_MISSION.get(osd_id, {})
        mission = mission_info.get("mission", "unknown")
        mission_clean = mission.replace(" ", "_").replace("/", "_").replace("+", "_")
        osd_dir = DATA_DIR / tissue / mission_clean
        glds_prefix = f"GLDS-{osd_id.replace('OSD-', '')}"
        if osd_id == "OSD-686":
            glds_prefix = "GLDS-617"

        counts_file = find_normalized_counts_file(osd_dir, glds_prefix)
        meta_file = find_sample_table_file(osd_dir, glds_prefix)

        counts_ok = "✅" if counts_file else "❌"
        meta_ok = "✅" if meta_file else "⚠️"
        print(f"  {osd_id} ({mission}): counts={counts_ok} meta={meta_ok}")
        if counts_file:
            print(f"    counts: {counts_file.name}")
        else:
            print(f"    counts: NOT FOUND in {osd_dir}")


def main():
    args = parse_args()

    if args.check_input:
        check_input_files(args.check_input)
        return

    # v2.0: Enrich existing metadata with temporal columns
    if args.enrich_temporal:
        tissues = list(TISSUE_OSD_MAP.keys()) if args.all else ([args.tissue] if args.tissue else list(TISSUE_OSD_MAP.keys()))
        print("=" * 60)
        print("GeneLab_benchmark — Temporal Metadata Enrichment (v2.0)")
        print(f"Tissues: {tissues}")
        print("=" * 60)
        for tissue in tissues:
            enrich_temporal_metadata(tissue, verbose=args.verbose)
        print("\n✓ Temporal enrichment complete.")
        return

    tissues_to_process = []
    if args.all:
        tissues_to_process = list(TISSUE_OSD_MAP.keys())
    elif args.tissue:
        tissues_to_process = [args.tissue]
    else:
        print("Specify --tissue TISSUE or --all")
        print(f"Available tissues: {list(TISSUE_OSD_MAP.keys())}")
        return

    print("=" * 60)
    print("GeneLab_benchmark — Quality Filter")
    print(f"Tissues: {tissues_to_process}")
    print(f"Output: {PROCESSED_DIR}")
    print("=" * 60)

    for tissue in tissues_to_process:
        process_tissue(tissue, verbose=args.verbose)

    print("\n✓ QC filtering complete.")
    print(f"  Processed data: {PROCESSED_DIR}")
    print(f"  QC reports:     {QC_DIR}")
    print()
    print("Next: python scripts/generate_tasks.py --tissue liver")
    print("      (applies LOMO-aware variance filter + generates task splits)")


if __name__ == "__main__":
    main()
