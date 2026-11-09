# Changelog

All notable changes to SMLMAnalysis.jl are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Before 1.0, a release that raises
the last number (0.1.Z) is non-breaking; one that raises the middle number
(0.Y.0) may break.

## [Unreleased]

### Changed

- The `DetectFitConfig` default detection box is now 9 pixels (was 11). This changes
  results for callers that rely on the default; pass `BoxerConfig(boxsize = 11, ...)` to
  keep the old behavior.

## [0.2.0] - 2026-10-21

### Breaking

The SMART and MIC H5 readers follow the lab image layout (decision 0040): images load as
`[y, x]` then frame (first index = row = y) and camera maps as `(ny, nx)`. Some files now
load transposed relative to 0.1.1, so x and y swap in results computed from them. Positions,
renders, overlays and anything else that depends on orientation change; distances,
densities, cluster and track statistics do not. Do not mix such results with 0.1.1 output
without swapping x and y.

- MATLAB-era SEQ files (MIC, written by `MIC_SEQ_SRcollect`): data and calibration maps now
  load as `[y, x]`; 0.1.1 loaded them `[x, y]`.
- SMART files written after 2026-01-29 (tagged `dimension_order`): `smart_h5_to_array` and
  file-based detectfit now return `[y, x]`; 0.1.1 transposed them to `[x, y]`.
- SMART files written before 2026-01-29 (untagged): `load_smart_h5` now converts them to
  `[y, x]`; 0.1.1 returned the stored `[x, y]`. `smart_h5_to_array` and file-based detectfit
  already returned `[y, x]` for these files and are unchanged.
- SEQv2 files from MicroscopeSeqSR load as stored, as in 0.1.1.
- Camera maps must be `(ny, nx)` of the data: a calibration map of another size, or an
  unrecognized `dimension_order`, throws `ArgumentError`.

### Changed

- How a file's stored layout is decided: a dataset's `dimension_order` attribute wins,
  matched on exact spellings (`"y,x,frame"`, `"y,x"`, `"y,x,z"`, `"HWN"`, `"HW"` are
  `[y, x]`; `"x,y,frame"`, `"x,y"`, `"x,y,z"` are `[x, y]`). Untagged SMART files are
  legacy `[x, y]`. An untagged MIC file is `[y, x]` only when its `Metadata` `Writer` starts
  with `MicroscopeSeqSR`; otherwise it is legacy MATLAB `[x, y]`. Legacy data and maps are
  converted together, once, at load. An untagged map takes its file's data layout.
- The tag of every MIC block is checked before any block is read, including with
  `max_blocks` and in the calibration and camera loaders.
- `load_smart_h5_info` and `load_mic_h5_info` report `width`/`height` as `nx`/`ny` of the
  array the loaders return.
- No change to `max_drift_nm` since 0.1.1. Users upgrading from 0.1.0 also get the 0.1.1
  change below (#69): registered-mode values are not comparable with 0.1.0.

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

[Unreleased]: https://github.com/JuliaSMLM/SMLMAnalysis.jl/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/JuliaSMLM/SMLMAnalysis.jl/compare/v0.1.1...v0.2.0
[0.1.1]: https://github.com/JuliaSMLM/SMLMAnalysis.jl/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/JuliaSMLM/SMLMAnalysis.jl/releases/tag/v0.1.0
