```@meta
CurrentModule = SMLMAnalysis
```

# Troubleshooting

A practical guide to the failure modes you are most likely to hit. Each entry is
**symptom → likely cause → fix**, with a link to the step page that covers it in
depth. Most problems are either a step-ordering mistake (caught by Julia's
dispatch as a `MethodError`) or a parameter mismatched to the data or optics.

## `MethodError` from `analyze` (wrong step order)

**Symptom.** An error like
`MethodError: no method matching analyze(::Vector{...}, ::FilterConfig)`.

**Cause.** A step that consumes a `BasicSMLD` was placed before any
`DetectFitConfig`, or the steps are otherwise in an impossible order. Only
`DetectFitConfig` turns raw image data into an `smld`; every other step expects
an existing `smld`, so dispatching one on a raw image stack finds no method.

**Fix.** Make `DetectFitConfig` the first step:

```julia
steps = [
    DetectFitConfig(boxer = BoxerConfig(boxsize = 9, psf_sigma = 0.130)),
    FilterConfig(photons = (500.0, Inf)),   # now receives a BasicSMLD
    RenderConfig(zoom = 20),
]
```

This is dispatch working as designed — a wrong order is a clear error, not a
silent wrong answer. See [The Pipeline Model](@ref).

## `Pipeline produced no SMLD`

**Symptom.** `Pipeline produced no SMLD. Did you include a DetectFitConfig step?`

**Cause.** The `steps` vector contains no `DetectFitConfig`, so the orchestrator
never produced any localizations to thread through the pipeline.

**Fix.** Add a `DetectFitConfig` as the first step (see
[Detection & Fitting](@ref "Detection & Fitting")). If you are resuming from a
saved `smld` (`smld = load_smld("after_detectfit.h5")`), run the later steps
directly on it rather than through a pipeline that starts from images.

## Too many or too few detections

**Symptom.** The detection/fit overlay boxes noise instead of real spots, or
misses obvious molecules; `DetectFitInfo.n_fits` is far below `n_rois`.

**Cause.** `BoxerConfig.boxsize` and `psf_sigma` are mismatched to your optics:
a wrong `psf_sigma` (expected PSF width, in microns) mis-thresholds detection,
and a `boxsize` too small for the PSF (or too large, merging neighbors) makes the
fit reject candidates. A large `n_rois`–`n_fits` gap is the tell.

**Fix.** Tune `boxsize` and `psf_sigma` to your pixel size and PSF, and choose a
`fitter` PSF model that matches (`GaussianXYNBS` for free width, `AstigmaticXYZNB`
for 3D); inspect the overlay figure after each adjustment. See
[Detection & Fitting](@ref "Detection & Fitting").

## Large inter-dataset shift warnings / bad drift

**Symptom.** A `PROGRESS`-level warning that inter-dataset shifts exceed ~500 nm,
a jagged `drift_trajectory.png`, or a smeared reconstruction after drift
correction.

**Cause.** The wrong `DriftConfig.dataset_mode`. Correcting one continuous movie
as `:registered` (or independent overlapping acquisitions as `:continuous`)
produces large spurious inter-dataset shifts.

**Fix.** Match the mode to how the data was acquired — `:continuous` for one long
movie (e.g. `DriftConfig(degree = 3, dataset_mode = :continuous,
chunk_frames = 4000)` for long runs), `:registered` for multiple overlapping
acquisitions of the same field. See [Drift Correction](@ref). Note that
per-dataset frame numbering is by design (required by the Legendre
normalization), not a bug.

## Crowded data / suspiciously bright localizations

**Symptom.** A few localizations with abnormally high photon counts, or a
reconstruction that looks denser/blurrier than expected.

**Cause.** When two fluorophores are active inside one diffraction-limited spot,
the fit returns a single localization with a corrupted position and an inflated
photon count. Separately, frame connection can over- or under-link in dense data.

**Fix.** Add an [Intensity Filter](@ref) (e.g.
`IntensityFilterConfig(cutoff = 0.01, field_mode = :gaussian)`) to reject
multi-emitter events with its Poisson upper-tail test; check its `p2_estimate` to
gauge how crowded the data is. For [Frame Connection](@ref), tune `max_frame_gap`:
too large over-links distinct emitters into one track, too small fragments one
emitter into many. A `compression` near 1 in the summary means almost nothing
connected.

## Out-of-memory on large acquisitions

**Symptom.** Julia runs out of memory while loading or detecting on a big
acquisition.

**Cause.** Holding the entire image stack in memory at once.

**Fix.** Process the acquisition as multiple datasets — detection and fitting
loop over datasets individually, so memory use stays bounded. Use file-based
detection (a single file split into blocks, or one file per dataset) rather than
passing all images in memory:

```julia
# MIC blocks auto-detected as datasets, loaded block by block
DetectFitConfig(path = "data.h5", h5_format = :mic,
                boxer = BoxerConfig(boxsize = 9))
```

See [Multi-Dataset](@ref "Multi-Dataset Acquisitions").

## GPU / CUDA errors during fitting

**Symptom.** A CUDA-related error raised inside the [Detection & Fitting](@ref
"Detection & Fitting") step.

**Cause.** Fitting is GPU-accelerated through GaussMLE, which needs a working
CUDA setup. The error surfaces when CUDA is unavailable or misconfigured.

**Fix.** Confirm your CUDA installation (driver, `CUDA.jl` functional); see
[Installation & Setup](@ref) for requirements and the GaussMLE docs for its
GPU/CPU options.

## Empty output / no figures written

**Symptom.** The pipeline runs cleanly but no diagnostic figures or `stats.md`
appear in the step output folders.

**Cause.** Figures are gated on two things: an `outdir` must be set, and
`verbose` must be high enough. At `Verbosity.SILENT`/`PROGRESS` only counts and
timing are emitted — basic figures (overlays, fit-quality, drift plots) start at
`Verbosity.STANDARD`, diagnostics at `DETAILED`, movies at `DEBUG`.

**Fix.** Set an output directory and raise the verbosity:

```julia
config = AnalysisConfig(camera = cam, steps = [...],
                        outdir = "output/", verbose = Verbosity.DETAILED)
```

See [Running a Pipeline](@ref) for the verbosity levels and what each one writes.
