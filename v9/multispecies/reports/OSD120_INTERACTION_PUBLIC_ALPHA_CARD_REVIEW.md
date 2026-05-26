# OSD-120 Public-Alpha Card Draft Review

Status: V9-MULTI-028 complete as of 2026-05-26; wording refreshed after
V9-MULTI-029 and V9-MULTI-030.

This block converts the OSD-120 diagnostic evidence set into a source-specific
public-alpha card draft. The card is not a release announcement and does not
promote OSD-120 to a frozen benchmark. It gives an external reader enough
context to inspect the diagnostic package without confusing it with a complete
OSDR mirror, leave-one-mission-out task, or validated biological claim.

Primary outputs:

- `interaction_public_alpha_card/public_alpha_card.md`
- `interaction_public_alpha_card/public_alpha_card_summary.csv`
- `interaction_public_alpha_card/public_alpha_card_summary.json`

## Card Scope

The card records:

- source id and source URL for OSD-120 / GLDS-120;
- organism, biospecimen, assay modality, task id, candidate id, and variant id;
- payload-freeze boundary from V9-MULTI-027;
- the three focus-fold diagnostic table rows;
- summary balanced-accuracy and nearest-centroid comparison values;
- seven allowed claims from the claim-to-artifact map;
- five disallowed claims;
- files an external reader should inspect;
- external context links;
- remaining release work.

## Boundary Decision

The card uses this release language:

- draft diagnostic alpha card;
- not a frozen benchmark release;
- diagnostic-required payloads are MD5 matched;
- full OSDR processed payload mirror is not claimed;
- within-source OSD-120 diagnostic only.

This resolves the immediate "OSD-120-specific public card" documentation gap
for the local diagnostic package, while preserving the remaining release
boundary: source release-target promotion, archive identifier, and license
decisions are still pending. The V9-MULTI-029 rebuild/environment gate is
available as a packaging preflight, V9-MULTI-030 provides a public metadata
skeleton, and V9-MULTI-031 provides draft RO-Crate/Data Package export
descriptors, but none of these artifacts rerun model grids or promote the
package to a frozen benchmark release.

## Safe External Summary

The safe one-paragraph summary is:

OSD-120 is a draft SpaceBio-Bench v9 diagnostic card for an Arabidopsis root
bulk RNA-seq light/genotype interaction task derived from public NASA
OSDR/GeneLab processed data. The packaged sparse L1 candidate recovers three
predefined fragile focus folds relative to nearest centroid, and the two
diagnostic-required processed payloads match the OSDR MD5 manifest. The package
is not a frozen benchmark release, not a full OSDR payload mirror, and not a
leave-one-mission-out or biomarker-validation claim.

## Rebuild Gate

Use:

```bash
python3 scripts/rebuild_v9_osd120_diagnostic_package.py --repo-root . --mode preflight
```

This writes `interaction_rebuild_gate/rebuild_gate_summary.csv`,
`interaction_rebuild_gate/rebuild_gate_steps.csv`, and
`interaction_rebuild_gate/rebuild_gate_environment.csv`.

## Public Metadata Package

Inspect:

- `interaction_public_metadata_package/public_metadata_summary.csv`
- `interaction_public_metadata_package/source_release_target_decision.csv`
- `interaction_public_metadata_package/public_metadata_skeleton.json`

The metadata package separates a public-now diagnostic metadata draft from
blocked source-alpha, full-mirror, and frozen-benchmark release targets.

## RO-Crate And Citation Scaffold

Inspect:

- `interaction_ro_crate_citation_scaffold/ro_crate_export_summary.csv`
- `interaction_ro_crate_citation_scaffold/ro_crate_validation_table.csv`
- `interaction_ro_crate_citation_scaffold/citation_freeze_checklist.csv`
- `interaction_ro_crate_citation_scaffold/ro-crate-metadata.draft.json`
- `interaction_ro_crate_citation_scaffold/datapackage.draft.json`

The scaffold is machine-readable but intentionally not archive-ready because
DOI/archive identifier, creator/contributor, and license/rights fields remain
unresolved.

## Next Block

Proceed to `V9-MULTI-032: archive identifier and license decision gate`. The
next block should decide the archive identifier/version, creator/contributor,
and license/rights path before any citable archive is attempted.
