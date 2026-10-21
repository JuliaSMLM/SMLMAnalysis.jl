# Changelog

All notable changes to SMLMAnalysis.jl are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Before 1.0, a release that raises
the last number (0.1.Z) is non-breaking; one that raises the middle number
(0.Y.0) may break.

## [Unreleased]

### Changed

- The SMART and MIC H5 readers follow the lab image layout (decision 0040): images load as
  `[y, x]` then frame and camera maps as `(ny, nx)`. A dataset's `dimension_order`
  attribute decides its stored layout; untagged SMART files and untagged MIC files whose
  `Metadata` `Writer` does not start with `MicroscopeSeqSR` (MATLAB `MIC_SEQ_SRcollect`)
  are legacy `[x, y]`, converted once at load together with their calibration maps.
  - Legacy MATLAB MIC data and maps now load transposed relative to 0.1.1, so
    localizations from those files swap x and y relative to 0.1.1 runs.
  - Tagged SMART files (`dimension_order` present) are no longer transposed twice; this
    also differs from 0.1.1.
  - `load_smart_h5_info` reports `width`/`height` as `nx`/`ny`; an unrecognized
    `dimension_order` throws `ArgumentError`.

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
  - `04_driftcorrect/stats.md` now reads `**Max drift**: <x> nm from <reference>` where
    0.1.0 read `**Max intra-dataset drift**: <x> nm`. Scripts should read
    `stepinfo(info, :driftcorrect).summary[:max_drift_nm]` rather than parse the file.

### Changed

- Requires SMLMDriftCorrection 0.2.9 (was 0.2.8). 0.2.9 fixes a registered-mode regression
  in 0.2.7 and 0.2.8 that could shift whole datasets by microns while reporting convergence.
  Registered-mode results from SMLMAnalysis 0.1.0 with SMLMDriftCorrection 0.2.8 should be
  checked (`max_intershift_nm`) or rerun.
- `DriftConfig.shift_scale` has had no effect since SMLMDriftCorrection 0.2.7. From 0.2.9 a
  non-default value logs a one-time warning; drop the setting to silence it.

## [0.1.0] - 2026-10-11

First registered release. The merged pull requests are listed in the
[v0.1.0 release notes](https://github.com/JuliaSMLM/SMLMAnalysis.jl/releases/tag/v0.1.0).

[Unreleased]: https://github.com/JuliaSMLM/SMLMAnalysis.jl/compare/v0.1.1...HEAD
[0.1.1]: https://github.com/JuliaSMLM/SMLMAnalysis.jl/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/JuliaSMLM/SMLMAnalysis.jl/releases/tag/v0.1.0
