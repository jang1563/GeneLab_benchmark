# Contributing To SpaceBio-Bench

Thank you for helping improve SpaceBio-Bench / GeneLab Benchmark. The most
useful contributions are small, reproducible, and tied to the public release
surface.

## Good First Contributions

- Report a broken download path, stale link, or unclear dataset-card wording.
- Reproduce a public fold or release QA command and share the exact command and
  output.
- Improve documentation around task layout, submission format, or citation.
- Submit benchmark predictions using the public JSON format.

## Benchmark Submissions

Prediction submissions should follow
[docs/submission_format.md](docs/submission_format.md). Include:

- task ID and fold IDs;
- model name and short model description;
- exact preprocessing or feature representation;
- the command used to validate or evaluate the submission.

Use a GitHub issue for discussion before opening a pull request with benchmark
results or new public artifacts.

## Pull Requests

Keep pull requests focused. For documentation, metadata, or public-card changes,
run:

```bash
make release-qa
make hpc-public-qa
```

For code changes that touch benchmark behavior, also run the relevant unit tests
or scripts and describe what changed in the PR body.

## Public Release Boundaries

Use the release labels in [README.md](README.md) and
[release/release_manifest.json](release/release_manifest.json). In short:

- v7.1 is the canonical historical result surface.
- v7.1.2 is a public-card and metadata patch over v7.1 results.
- v9 public bulk is a metadata catalog and baseline-summary surface.

Please avoid mixing future-lane claims into the v7.1 result surface.
