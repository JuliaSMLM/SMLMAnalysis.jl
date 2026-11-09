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
- DetectFit drops fits whose position lies outside the camera. GaussMLE lets a fit move
  1.5 pixels past its box, so a box on the camera edge could return a position outside the
  image, with an uncertainty of hundreds of nanometers. The count is `n_out_of_field` in
  the DetectFit info, its log line and `stats.md`.
- `DetectFitConfig` now throws `ArgumentError` when `camera` is given together with
  `pixel_size` or `qe != 1.0`. Those values were ignored before; build the camera with
  `build_camera_from_mic_h5(path; pixel_size, qe)` instead.

### Fixed

- Automatic density filtering preserves sparse populations with low neighbor-count peaks, rejecting only localizations without neighbors.
- DetectFit skips datasets with no detections while retaining their dataset indices and fit info; all-empty runs report the detection settings error.
- Quality filtering and diagnostics support emitters without p-values; requesting `pvalue` filtering warns and skips that criterion.
- DetectFitConfig rejects conflicting data-source keywords (`path` with `paths`, or `dataset_frames` without `path`) instead of silently ignoring them.
- `Checkpoint.END` saves the last localization checkpoint when rendering or statistics steps follow it. A user-defined step still counts as a checkpoint writer, as documented in Extending the Pipeline.
- FFTA 0.3.2 is excluded. That release, loaded through Makie, conflicts with FFTW and breaks drift's cross-correlation: `method = :fft` alignment crashes, and drift silently loses its cross-correlation start. FFTA 0.3.3 fixes the conflict (JuliaMath/FFTA.jl#146) and is allowed, as are versions before 0.3.2.
- `save_smld` records the emitters' own number type, taken from the coordinate field, instead of the camera's, so a result reloads as it was saved.
- `FilterConfig.psf_sigma`, `CompositeRenderConfig.clip_percentile` and `DensityFilterConfig.min_neighbors` accept whole-number input, as the other ranges already did.
- A MIC frame range reads only the blocks it overlaps; a range late in a file read every earlier frame into memory. One-frame blocks stored as 2-D datasets still load.
- The fit overlay colors fits rejected by the `z` (blue) and `sigma_z` (magenta) filters; they were drawn as passing.
- With `outdir`, the log-file tee forwards to the caller's logger only what that logger accepts: a run under `NullLogger()` no longer throws, an Error-level logger no longer receives Info lines, and `log.txt` keeps every message it recorded before.
- `build_camera_from_mic_h5` on a MIC file without a Calibration group throws `ArgumentError` naming the missing entries and the camera to pass instead, not a `KeyError` from HDF5.
- The docstring, README and documentation examples set `psf_sigma` with `boxsize`; without it, detection had no threshold and boxed every local maximum.
- DetectFit warns when cleanup removes every fit, with the counts removed as non-finite and as outside the camera.
- Diagnostic plots no longer throw when their data has a single value (a flat sample frame in the detection or fit overlay, or constant photon counts in the intensity-filter figures) with Makie older than 0.24.14.
- The color range of a diagnostic plot stays finite and distinct when the data's value is too large for Makie's single-precision range (above about 3.4e38).
- Filter and intensity-filter figures no longer throw on axis limits that collapse to one value or go infinite: one localization, zero background, zero photons, or an open filter bound such as `FilterConfig(precision = (0, Inf))`. An open bound draws no rejected region.
- The filter p-value panel no longer throws when every p-value is 1; it spans one decade, and its bars, density axis and limits use the same bin edges.
- The filter precision panel no longer throws when every precision is above its 1 µm display cutoff.
- The detection and fit overlays no longer throw on sample frames that contain NaN pixels; the contrast stretch uses the finite pixels.
- The fit overlay applies an explicit `FilterConfig(psf_sigma = (0, hi))`: fits the filter rejects for PSF width are drawn cyan, not green. Only `:auto` is skipped when no PSF-width mode can be found.
- The filter figures draw every bound that rejects a fit, and the fits it rejects: an upper photon bound, a lower `sigma_z` bound, an upper p-value bound and the anisotropic PSF bounds were applied but not shown, and a shaded region could hold none of its rejected fits because the histogram stopped at the data's 2-98% range. Each fit is now placed by the filter's own test on its original value: a shaded region holds exactly the fits that side of the range rejects, every bound is a bar edge, and bars are drawn on the exact edges. Rejected fits beyond the axis are counted in the end bar; passing fits outside the panel's view are clipped, as before. The panel title names every other fit not drawn where it sits: rejected fits with a NaN value (a NaN p-value, which DetectFit keeps), rejected fits a transform leaves no place for (p = 0 under an upper p-value bound of 0), and passing fits at coincident bounds (panels that overlay σ_x and σ_y count each component's values by name); and it names a rejecting bound too large to draw, including a lower precision bound. A positive lower `z` bound or negative upper one no longer falls outside the view. A lower precision bound is a dashed line, not a shaded region, because the filter tests the larger of σ_x and σ_y.
- The fit overlay legend prints each range exactly as the filter tests it, in the filter's units (μm for precision and PSF width), including the p-value range (which a NaN p-value also fails), the PSF-width bounds the filter applies, and ranges that reject everything. It and the photon panel no longer throw on a photon bound or count too large for an integer, such as `FilterConfig(photons = (200, 1e20))`. The legend wraps to the figure's width, and the figure widens when a frame is too narrow for it, so a long legend (the tutorial's filter with `psf_sigma = :auto`) is no longer cut off.

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
