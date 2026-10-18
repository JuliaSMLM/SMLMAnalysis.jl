# Changelog

All notable changes to SMLMAnalysis.jl are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Before 1.0, a release that raises
the last number (0.1.Z) is non-breaking; one that raises the middle number
(0.Y.0) may break.

## [Unreleased]

## [0.1.1] - 2026-10-18

### Fixed

- `driftcorrect`: `max_drift_nm` is now measured from the reference frame that
  `DriftConfig.dataset_mode` implies (#69). The step's log line, `stats.md` and per-dataset
  table use the same definition.
  - `:registered` (the default): from each dataset's own first frame, so the value no longer
    includes that dataset's registration shift, which is still reported as
    `max_intershift_nm`. In 0.1.0 a 5 µm misregistration read as 4650 nm of drift.
  - `:continuous`: from the acquisition's first frame, as in 0.1.0.
  - Registered-mode values are not comparable with 0.1.0: the inter shift is no longer
    included, and a steady drift now reads its full first-to-last size.
  - Anything that gates on `max_drift_nm` must also check `max_intershift_nm`.

### Changed

- Requires SMLMDriftCorrection 0.2.9 (was 0.2.8). 0.2.9 fixes a registered-mode regression
  in 0.2.7 and 0.2.8 that could shift whole datasets by microns while reporting convergence.
  Registered-mode results from SMLMAnalysis 0.1.0 with SMLMDriftCorrection 0.2.8 should be
  checked (`max_intershift_nm`) or rerun.
- From SMLMDriftCorrection 0.2.9, `DriftConfig.shift_scale` is unused; a non-default value
  logs a one-time warning.

## [0.1.0] - 2026-10-11

First registered release. The merged pull requests are listed in the
[v0.1.0 release notes](https://github.com/JuliaSMLM/SMLMAnalysis.jl/releases/tag/v0.1.0).

[Unreleased]: https://github.com/JuliaSMLM/SMLMAnalysis.jl/compare/v0.1.1...HEAD
[0.1.1]: https://github.com/JuliaSMLM/SMLMAnalysis.jl/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/JuliaSMLM/SMLMAnalysis.jl/releases/tag/v0.1.0
