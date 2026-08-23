# Changelog

Notable changes to `genomic-benchmarks-qc`. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); this project aims at
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Flags are thresholds on a computed score, so **any change to how a score is
computed or where a boundary sits is a breaking change** for anyone gating a
build on the output, and is called out as such here.

## [Unreleased]

### Added

- A documentation site at
  <https://genomic-benchmarks.github.io/genomic-benchmarks-qc/>, with eight
  worked examples, a task-shaped guide, and CLI and API references generated
  from the code.
- `examples/`: eight example datasets chosen so that every check is exercised at
  every severity it can reach, each with its provenance and expected flags in a
  `meta.toml`, plus `examples/build.py` to regenerate their reports and assert
  the flags still match.

### Changed

- `example_datasets/` moved to `examples/<name>/data/`. Commands in the README
  and CI were updated; anything that referenced the old paths needs the same.
- `example_outputs/` removed. Reports are built from the committed data and
  published to the docs site rather than committed, so a published report is
  always the one the current code produces. Their URLs moved from
  `/example_outputs/<dataset>/...` to `/reports/<example>/...`.
- `utils.testing.direct_feature_model` and
  `utils.testing.flag_significant_differences` gained type annotations. No
  behaviour change.

## [0.9.0] and earlier

No changelog was kept. See the
[commit history](https://github.com/genomic-benchmarks/genomic-benchmarks-qc/commits/main).
