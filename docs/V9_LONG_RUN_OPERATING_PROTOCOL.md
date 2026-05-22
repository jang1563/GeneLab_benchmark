# v9 Long-Run Operating Protocol

Status: active operating protocol
Date: 2026-05-21

This protocol exists because v9 work should not proceed as a string of tiny
1-3 minute sessions. When the user says "continue", the default operating mode
is a 30-60 minute uninterrupted work block with local implementation, external
research, documentation, and verification kept in one continuous arc.

## Default Block Contract

Use this contract unless the user explicitly asks for a short answer, a quick
status report, or a pause.

- Target duration: 45 minutes.
- Minimum meaningful duration: 30 minutes.
- Maximum ordinary duration: 60 minutes before a checkpoint.
- Do not send a final answer after a small subtask if the active block still has
  useful work remaining.
- Send short progress updates while working, but reserve the final answer for a
  real checkpoint.
- If a command or data run is still active, keep the turn alive and continue
  polling until the run finishes or clearly fails.
- If the user sends a continue request, resume the active block from the latest
  backlog state, not from a fresh planning reset.

## Three-Lane Execution Model

Every long block should run three lanes in parallel where possible.

### Lane A: Local Build

Concrete repository work:

- scripts
- package APIs
- manifests
- generated reports
- tests
- docs/backlog updates

### Lane B: External Deep Research

Primary-source research that changes the next implementation choices:

- OSDR API and OSDR citation/access guidance.
- FAIR data and metadata principles.
- BagIt checksum/package integrity rules.
- RO-Crate metadata/provenance structure.
- DataCite metadata and versioning expectations.
- Frictionless Data Package / Tabular Data Package conventions for CSV-heavy
  artifacts.

Research should end in local decisions, not just notes. Each source scan should
answer:

- What field or artifact should v9 add?
- What should v9 avoid claiming?
- What evidence is needed before release language is safe?
- Which script/doc/schema should change next?

### Lane C: Verification And Continuity

Every long block should leave the next block easy to start:

- tests or explicit reason tests were not run
- regenerated artifact paths
- no local absolute-path leaks in public v9 artifacts
- backlog status update
- next active block named in `v9/OPERATING_BACKLOG.md`

## Current External Research Anchors

These are the primary references for the next provenance/release blocks.

| Topic | Source | Operating Implication |
|---|---|---|
| OSDR metadata | NASA OSDR Biological Data API: `https://visualization.osdr.nasa.gov/biodata/api/` | Use dataset, assay, sample, and file-list endpoints to derive source evidence instead of hand-maintaining accessions only. |
| OSDR citation/access | NASA OSDR FAQ: `https://science.nasa.gov/reference/osdr-faq/` | Dataset citations and OSDR credit language should be explicit in dataset-card and release docs. Public vs controlled human data must remain separated. |
| FAIR | Wilkinson et al. 2016 FAIR Guiding Principles: `https://www.nature.com/articles/sdata201618` | v9 source freeze needs identifiers, access metadata, interoperability metadata, and reuse constraints, not only files. |
| BagIt | RFC 8493 BagIt: `https://www.rfc-editor.org/rfc/rfc8493` | Checksum audit should produce a file-level payload manifest; a source is not "frozen" until listed payload checksums verify. |
| RO-Crate | RO-Crate technical overview: `https://www.researchobject.org/ro-crate/technical_overview` | A future `ro-crate-metadata.json` should describe files, authors/license/identifiers/provenance, and workflow outputs together. |
| DataCite | DataCite Metadata Schema: `https://schema.datacite.org/` | Release metadata should align with DOI-oriented dataset metadata: creators, titles, identifiers, related identifiers, versions, and resource type. |
| Data Package | Frictionless Data Package: `https://specs.frictionlessdata.io/guides/data-package/` | CSV-heavy v9 tables should have machine-readable schemas and resource descriptors, even before RO-Crate export. |

## Next Five Long Blocks

### Block 1: Source Checksum Audit

Target: 60 minutes.

Status: complete as of 2026-05-21 for OSDR API and checksum-manifest evidence.

Goal:

Promote `v9/source_inventory.csv` from accession inventory to checksum-audit
input.

Work:

- Inspect OSDR API file-list shape for representative sources such as OSD-48,
  OSD-101, OSD-686, and OSD-397.
- Implement `scripts/audit_v9_source_checksums.py`.
- Generate `v9/source_checksum_audit.csv` and `.json`.
- Record whether each source has file-level checksum evidence, API evidence, or
  an explicit pending reason.
- Update backlog and long-horizon plan.

Done when:

- All 22 source rows have an audit status.
- Failures are captured as data, not as crashed runs.
- Tests cover parsing and failure reporting.

Actual result:

- `scripts/audit_v9_source_checksums.py` generated
  `v9/source_checksum_audit.csv` and `.json`.
- All 22 public bulk source rows returned `api_status=ok` and
  `audit_status=checksum_manifest_parsed`.
- The audit found 39 checksum manifests, parsed 8,439 MD5 entries, and matched
  8,275 entries back to OSDR file-list payload names.
- All rows keep `freeze_ready=false` because payload files have not been
  downloaded and locally hashed.

### Block 2: Public Data Package Design

Target: 45-60 minutes.

Status: complete as of 2026-05-21 for draft package boundary and descriptor.

Goal:

Design how v9 public bulk artifacts should be packaged outside Git.

Work:

- Draft `docs/V9_PUBLIC_BULK_PACKAGE_DESIGN.md`.
- Decide which artifacts belong in a lightweight Git repo versus downloadable
  data bundle.
- Define a future `datapackage.json` shape for source inventory, task manifests,
  fold data index, baseline summaries, and reports.
- Add schema/test stubs only if the package boundaries are clear.

Done when:

- Public bundle boundaries are documented.
- No release language claims a frozen dataset yet.

Actual result:

- `docs/V9_PUBLIC_BULK_PACKAGE_DESIGN.md` defines metadata spine, public bulk
  payload bundle, benchmark output bundle, and deferred/excluded artifacts.
- `scripts/build_v9_datapackage_draft.py` generates
  `v9/datapackage.draft.json`.
- The draft descriptor has 11 resources, including task/source/provenance
  tables, eight task manifests, 24 baseline prediction CSVs, 24 metrics JSON
  files, and 24 run-manifest JSON files.
- Large fold payload files are indexed by `v9/task_data_index.csv` but remain
  outside payload-level hash freeze.

### Block 3: RO-Crate Export Design

Target: 45-60 minutes.

Goal:

Turn the benchmark scaffold into a citable research object plan.

Work:

- Draft `docs/V9_RO_CRATE_EXPORT_DESIGN.md`.
- Map source inventory, task manifests, reports, scripts, and benchmark outputs
  into RO-Crate entities.
- Specify provenance links from source accession to task manifest to predictions
  to metrics.

Done when:

- Export design names entities, identifiers, files, licenses, authorship,
  workflow runs, and checksum dependencies.

### Block 4: Dataset Card And Citation Draft

Target: 45-60 minutes.

Status: complete as of 2026-05-21 for draft dataset-card language.

Goal:

Prepare release-facing language without overclaiming.

Work:

- Draft `docs/v9_hf_dataset_card.md`.
- Include OSDR credit/citation language, source table, task table, metrics,
  limitations, privacy boundaries, and release status.
- Separate public bulk alpha from future single-cell/human/gated tracks.

Done when:

- An external reader can understand what v9 alpha is and is not.

Actual result:

- `docs/v9_hf_dataset_card.md` contains Hugging Face-style YAML metadata and
  release-facing sections for summary, use, structure, tasks, sources,
  provenance, baselines, privacy, limitations, license/citation, and maintenance.
- The card explicitly says the package is draft/not frozen and that checksum
  manifests have been parsed but payloads have not been locally hashed.
- OSDR citation and access guidance is linked to NASA OSDR FAQ and terms pages.

### Block 5: First Flagship Decision Research

Target: 60 minutes.

Goal:

Choose whether the next flagship implementation should be single-cell
spaceflight or radiation/stressor tasks.

Work:

- Inventory local RRRM/scRNA-seq assets.
- Inventory local v8 DECOMPOSE/radiation assets.
- Compare readiness, scientific novelty, implementation risk, and release
  safety.
- Update `docs/V9_DESIGN_OPTIONS.md` and `v9/OPERATING_BACKLOG.md`.

Done when:

- A concrete next flagship block is selected with files and tests named.

## Stop Conditions

Stop and checkpoint only when one of these is true:

- A 30-60 minute block reaches a coherent artifact boundary.
- A destructive or network-escalated action requires user approval.
- A dependency or data access issue blocks all useful local work.
- The user explicitly asks to stop, pause, or only report status.

## Next Active Block

`V9-ORG-014: Human organoid frozen DE contrast extraction and signature metric
contract` is the next active implementation block. It should start from
`v9/human_organoid/signature_reference_audit.draft.csv`,
`v9/human_organoid/reports/ORGANOID_SIGNATURE_METRIC_REFERENCE_DECISION.md`,
`v9/human_organoid/task_manifests/draft_human_organoid_spaceflight.json`, and
the OSDR DE/contrast files for OSD-863 and OSD-871. V9-ORG-013 established that
public OSDR DE tables and contrast definitions exist. The next block should
derive a compact, checksummed DE reference table, normalize log2 fold-change
orientation to `LEO_or_ISS - Ground`, define a `response_signature.csv`
submission artifact, and keep DE/signature metrics skip-aware until that
artifact is supplied.

`V9-THEN-005: RO-Crate Export Design` remains the next packaging/release block
when the active lane returns to release metadata. It should start from
`docs/V9_PUBLIC_BULK_PACKAGE_DESIGN.md`, `docs/v9_hf_dataset_card.md`,
`v9/datapackage.draft.json`, `v9/source_inventory.csv`, and
`v9/source_checksum_audit.csv`.
