# [Guide](@id Guide)

Conceptual reference for topics that need more explanation than a docstring provides.

## Pipeline Architecture

### Design philosophy

SMLMAnalysis uses **Julia's method dispatch as the pipeline routing mechanism**. There is no step registry, no string-based lookup, and no `if/elseif` chain in the pipeline loop. Each analysis step is a config type (`<: AbstractSMLMConfig`), and calling `analyze(state, config)` dispatches to the correct implementation based on the types of both arguments. This means:

- **Composability**: Steps can be freely reordered, repeated, or omitted -- the pipeline loop doesn't care what the steps are, only that `analyze()` has a method for `(state_type, config_type)`.
- **Extensibility**: Adding a new step requires no changes to the pipeline orchestrator. Define a config type and an `analyze()` method, and the step works automatically.
- **Error clarity**: Wrong step ordering produces a `MethodError` (e.g., `FilterConfig` before `DetectFitConfig` gives "no method matching `analyze(::Vector{...}, ::FilterConfig)`"), not a silent failure.

Every `analyze()` call returns a `(result, info)` tuple following the JuliaSMLM tuple-pattern. The result becomes the input state for the next step; the info is collected for provenance.

### The `steps` vector

`AnalysisConfig.steps` is a composable list of `AbstractSMLMConfig` subtypes. The orchestrator in `analysis.jl` iterates through them sequentially, threading state between steps:

```
DetectFitConfig (required first)
        │ produces smld
        ▼
  ┌─────────────────────────────────┐
  │ FilterConfig         (0+ times) │
  │ IntensityFilterConfig(0+ times) │
  │ FrameConnectConfig   (0-1)      │  ← Any order,
  │ DriftConfig          (0-1)      │    any combination
  │ DensityFilterConfig  (0+ times) │
  │ BaGoLConfig          (0-1)      │
  │ RenderConfig         (0+ times) │
  └─────────────────────────────────┘
```

`CalibrationConfig` is **not** a top-level pipeline step — it is a sub-config of `FrameConnectConfig` (set via the `calibration=` field). When provided, the calibration runs inside frame connection (link → calibrate → combine) and its results are exposed on `FrameConnectInfo.calibration`.

The only ordering constraint: `DetectFitConfig` must be first because it is the only step that produces a `BasicSMLD` from raw image data. All other steps receive and return an existing `smld`.

### Two-layer step architecture

Each step has two layers:

1. **Internal function** (not exported) -- does the actual work, returns `(result, YourStepInfo)`:
   - `filter_step(smld, cfg; outdir, step_number, verbose)` → `(filtered_smld, FilterInfo)`
   - `densityfilter_step(smld, cfg; ...)` → `(filtered_smld, DensityFilterInfo)`
   - Or delegates to an upstream package: `SMLMDriftCorrection.driftcorrect(smld, cfg)` → `(corrected, DriftInfo)`

2. **`analyze()` dispatch** -- thin wrapper that times the call and creates a `StepInfo`:
   ```julia
   function analyze(smld::BasicSMLD, cfg::FilterConfig;
                    outdir=nothing, step_number::Int=0, verbose::Int=Verbosity.STANDARD, kwargs...)
       t = @elapsed (filtered, filter_info) = filter_step(smld, cfg;
           outdir=outdir, step_number=step_number, verbose=verbose)
       (filtered, StepInfo(step_number, cfg, t, _step_summary(filter_info); info=filter_info))
   end
   ```

The `kwargs...` in the `analyze()` signature is required -- the pipeline loop passes context keywords (like `smld_raw`) and each step must accept (and ignore) keywords it doesn't use.

### Pipeline orchestrator

The orchestrator (`_run_pipeline` in `analysis.jl`) is a simple fold:

```julia
for (i, cfg) in enumerate(steps)
    cfg = _prepare_step(cfg, camera)       # inject camera into DetectFitConfig
    (state, step_info) = analyze(state, cfg; outdir=outdir, step_number=i, verbose=v)
    push!(step_infos, step_info)
end
```

No `isa` checks, no step registry. Dispatch handles everything. After the loop, the orchestrator collects `StepInfo`s into `AnalysisInfo` and extracts special state (`smld_connected`, `drift_model`) from the collected infos.

### State threading

The orchestrator maintains `smld` as the working state. Each step receives the current `smld` and returns an updated one. Additional state captured from specific steps:

- `smld_connected` — from `FrameConnectConfig` (tracks with multi-frame info)
- `drift_model` — from `DriftConfig` (fitted drift polynomial)

### Repeatable and optional steps

- **`FilterConfig`**: Can appear multiple times (coarse filter early, tighter filter after connection)
- **`RenderConfig`**: Can appear multiple times (different zooms, colormaps, or strategies)
- **`DensityFilterConfig`**, **`IntensityFilterConfig`**: Can appear multiple times
- **`FrameConnectConfig`** and **`DriftConfig`**: Typically used once, but not enforced. Uncertainty calibration is configured inside `FrameConnectConfig` via the `calibration=` field; it is not a separate pipeline step.

### Config provenance

| Config | Defined in | Notes |
|--------|-----------|-------|
| `DetectFitConfig` | SMLMAnalysis | Wraps SMLMBoxer + GaussMLE internally |
| `FilterConfig` | SMLMAnalysis | Pure SMLMAnalysis logic |
| `IntensityFilterConfig` | SMLMAnalysis | Poisson upper-tail test against excitation field |
| `FrameConnectConfig` | SMLMFrameConnection | Re-exported via const alias |
| `CalibrationConfig` | SMLMFrameConnection | Sub-config of `FrameConnectConfig.calibration` (not a pipeline step) |
| `DriftConfig` | SMLMDriftCorrection | Re-exported via const alias |
| `DensityFilterConfig` | SMLMAnalysis | Pure SMLMAnalysis logic |
| `BaGoLConfig` | SMLMAnalysis | Wraps `SMLMBaGoL.run_bagol` |
| `RenderConfig` | SMLMRender | Re-exported via const alias |

SMLMAnalysis defines some step configs locally (`DetectFitConfig`, `FilterConfig`, `IntensityFilterConfig`, `DensityFilterConfig`, `BaGoLConfig`) and re-exports others from upstream packages (`FrameConnectConfig`, `CalibrationConfig`, `DriftConfig`, `RenderConfig`). Extending the pipeline with a new upstream package follows the same re-export pattern.

## Adding a Custom Step

### Overview

Adding a new step to SMLMAnalysis requires **no changes to the pipeline orchestrator**. The dispatch-based architecture means you only need to:

1. Define a config type (`<: AbstractSMLMConfig`)
2. Define an info type (`<: AbstractSMLMInfo`)
3. Implement the internal work function
4. Implement the `analyze()` dispatch wrapper
5. Implement `_step_summary()` for the info type
6. Include and export from the module

### Worked example: adding a spatial ROI filter

Here is a complete, minimal step that filters localizations to a spatial region:

**Step 1: Config and info types**

```julia
# src/steps/spatialfilter.jl

"""
    SpatialFilterConfig <: AbstractSMLMConfig

Filter localizations to a rectangular spatial region.

# Keywords
- `x_range`: (min, max) in microns
- `y_range`: (min, max) in microns
"""
@kwdef struct SpatialFilterConfig <: SMLMData.AbstractSMLMConfig
    x_range::Tuple{Float64, Float64}
    y_range::Tuple{Float64, Float64}
end

"""Info from spatial filtering step."""
struct SpatialFilterInfo <: SMLMData.AbstractSMLMInfo
    n_before::Int
    n_after::Int
end
```

**Step 2: Internal work function**

The internal function does the actual computation. It receives the standard keyword arguments (`outdir`, `step_number`, `verbose`) and returns `(result, info)`:

```julia
function spatialfilter_step(smld::BasicSMLD, cfg::SpatialFilterConfig;
                            outdir::Union{String,Nothing}=nothing,
                            step_number::Int=0,
                            verbose::Int=Verbosity.STANDARD)
    v = verbose
    dir = step_outdir(outdir, step_number, cfg)

    v >= Verbosity.PROGRESS && @info "[$step_number] $(step_name(cfg))"

    n_before = length(smld.emitters)
    keep = [cfg.x_range[1] <= e.x <= cfg.x_range[2] &&
            cfg.y_range[1] <= e.y <= cfg.y_range[2] for e in smld.emitters]
    filtered = BasicSMLD(smld.emitters[keep], smld.camera, smld.n_frames, smld.n_datasets, smld.metadata)
    n_after = length(filtered.emitters)

    if dir !== nothing
        mkpath(dir)
        _save_config!(dir, cfg)        # writes config.toml
    end

    v >= Verbosity.PROGRESS && @info "  → $n_after / $n_before"
    (filtered, SpatialFilterInfo(n_before, n_after))
end
```

**Step 3: Summary dispatch**

`_step_summary` converts the info struct to a `Dict{Symbol, Any}` for the pipeline summary table:

```julia
_step_summary(info::SpatialFilterInfo) = Dict{Symbol,Any}(
    :n_before => info.n_before,
    :n_after => info.n_after,
    :acceptance => round(info.n_after / max(1, info.n_before), digits=3)
)
```

**Step 4: `analyze()` dispatch**

The thin wrapper times the call and creates a `StepInfo`. The `kwargs...` is required for pipeline compatibility:

```julia
function analyze(smld::BasicSMLD, cfg::SpatialFilterConfig;
                 outdir=nothing, step_number::Int=0, verbose::Int=Verbosity.STANDARD, kwargs...)
    t = @elapsed (filtered, sf_info) = spatialfilter_step(smld, cfg;
        outdir=outdir, step_number=step_number, verbose=verbose)
    (filtered, StepInfo(step_number, cfg, t, _step_summary(sf_info); info=sf_info))
end
```

**Step 5: Include and export**

In `src/SMLMAnalysis.jl`:

```julia
include("steps/spatialfilter.jl")
export SpatialFilterConfig, SpatialFilterInfo
```

**Usage** -- works immediately in both workflows:

```julia
# Config-driven
config = AnalysisConfig(
    camera = cam,
    steps = [
        DetectFitConfig(boxer=BoxerConfig(boxsize=9)),
        SpatialFilterConfig(x_range=(5.0, 15.0), y_range=(5.0, 15.0)),
        RenderConfig(zoom=20),
    ],
)
(result, info) = analyze(image_stacks, config)

# Step-by-step
(smld, _) = analyze(smld, SpatialFilterConfig(x_range=(5.0, 15.0), y_range=(5.0, 15.0)))
```

### Key patterns

**StepInfo creation**: Every `analyze()` dispatch must return `(result, StepInfo(...))`. The `StepInfo` constructor takes `(step_number, config, elapsed_s, summary_dict; info=typed_info)`.

**`_step_summary` dispatch**: Implement `_step_summary(info::YourInfo)` returning `Dict{Symbol, Any}`. This is used for the pipeline summary table (`summary.md`). The fallback returns an empty dict.

**Shared helpers** (from `src/steps/common.jl`):
- `step_outdir(outdir, step_number, cfg)` -- computes `outdir/02_spatialfilter/`
- `_save_config!(dir, cfg)` -- writes `config.toml` with all config fields
- `_save_info!(dir, info)` -- writes `info.toml` with scalar fields from the info struct
- `step_name(cfg)` -- derives name from type (e.g., `SpatialFilterConfig` → `"spatialfilter"`)

**Pipeline cache** (for inter-step data passing via filesystem):
- `save_cache(outdir, "filename.jld2"; key=value)` -- save data for downstream steps
- `load_cache(outdir, "filename.jld2")` -- load cached data (returns `nothing` if missing)

### Re-exporting from an upstream package

If the step is implemented in a separate package (like `DriftConfig` from SMLMDriftCorrection), use a const alias and delegate:

```julia
# In SMLMAnalysis.jl
const MyStepConfig = MyUpstreamPackage.MyStepConfig

# In steps/mystep.jl
function analyze(smld::BasicSMLD, cfg::MyUpstreamPackage.MyStepConfig;
                 outdir=nothing, step_number::Int=0, verbose::Int=Verbosity.STANDARD, kwargs...)
    t = @elapsed (result, upstream_info) = MyUpstreamPackage.mystep(smld, cfg)
    (result, StepInfo(step_number, cfg, t, _step_summary(upstream_info); info=upstream_info))
end
```

This keeps the upstream package independent while making its config a first-class pipeline step.

## Multi-Dataset Architecture

SMLM acquisitions are often split into multiple datasets -- either multiple files from the same sample or a single long movie divided into segments. SMLMAnalysis handles this natively.

### Why multi-dataset?

- **Memory efficiency**: Process one dataset at a time without loading all images
- **Drift correction**: Each dataset is corrected independently, then aligned via inter-dataset shifts
- **Frame numbering**: Frames are 1 to `n_frames_per_dataset` within each dataset, NOT global indices. This is required by SMLMDriftCorrection's Legendre polynomial normalization to [-1, 1]
- **Dataset tracking**: Every emitter carries a `dataset` field so downstream steps know which dataset it belongs to

### Configuration

Dataset boundaries are encoded in the data structure, not user-specified integers:

```julia
# In-memory: Vector{Array} defines dataset boundaries
image_stacks = [images1, images2, images3, images4]  # 4 datasets
(result, info) = analyze(image_stacks, config)

# Step-by-step: same Vector{Array} input
(smld, info) = analyze(image_stacks, DetectFitConfig(
    camera=camera, boxer=BoxerConfig(boxsize=9)))

# File-based: MIC format auto-detects blocks as datasets
config = AnalysisConfig(
    camera = cam,
    steps = [DetectFitConfig(path="data.h5", h5_format=:mic, boxer=BoxerConfig(boxsize=9)), ...],
)
(result, info) = analyze(config)  # No data argument needed

# Multiple files: one per dataset
config = AnalysisConfig(
    camera = cam,
    steps = [DetectFitConfig(paths=["d1.h5", "d2.h5", "d3.h5", "d4.h5"], boxer=BoxerConfig(boxsize=9)), ...],
)
(result, info) = analyze(config)
```

### SMLD structure

After detection and fitting, the combined SMLD has:
- `smld.n_frames` = frames per dataset (not total)
- `smld.n_datasets` = number of datasets
- Each `emitter.dataset` identifies its source dataset
- Each `emitter.frame` is relative to its dataset (1-based)

## Drift Correction Modes

`DriftConfig` supports two primary modes, controlled by the `dataset_mode` field.

### Registered mode (`dataset_mode=:registered`, default)

For multi-dataset acquisitions where datasets are spatially registered (e.g., the stage returns to approximately the same position between datasets):

```julia
DriftConfig(
    degree = 2,
    dataset_mode = :registered,
    quality = :singlepass
)
```

- Each dataset gets its own intra-dataset drift polynomial
- Inter-dataset alignment via entropy optimization over spatial overlap
- Warns if inter-dataset shifts exceed `intershift_threshold_nm` (default 500nm)

### Continuous mode (`dataset_mode=:continuous`)

For a single long acquisition split into chunks for processing:

```julia
# Short acquisition (< 4000 frames): single polynomial
DriftConfig(
    degree = 5,
    dataset_mode = :continuous,
    n_chunks = 0
)

# Long acquisition: chunked
DriftConfig(
    degree = 3,
    dataset_mode = :continuous,
    chunk_frames = 4000
)
```

- Drift accumulates continuously across chunks/datasets
- Chunking recommended for > 4000 frames (`chunk_frames=4000` is a reasonable maximum per chunk)
- Higher polynomial degree may be needed for complex drift patterns

### Quality settings

- `:singlepass` (default): Single entropy optimization pass. Fast and usually sufficient.
- `:iterative`: Iterates until convergence. Slower but may give better results for difficult drift patterns.

### Auto ROI

`auto_roi=true` (default) automatically selects a dense subregion for entropy estimation, which is faster and gives better signal for datasets with uneven emitter density.

## Saving and Resuming

The functional pipeline supports saving intermediate state for later resume:

### SMLD checkpoints

Save/load SMLD after expensive steps (detectfit) using HDF5:

```julia
# Save after expensive detectfit
(smld, info) = analyze(image_stacks, DetectFitConfig(
    camera=camera, boxer=BoxerConfig(boxsize=9)))
save_smld("output/after_detectfit.h5", smld)

# Resume later - try different filter parameters
smld = load_smld("output/after_detectfit.h5")
(smld, _) = analyze(smld, FilterConfig(photons=(300.0, Inf)))
```

### Full pipeline state

Save/load the complete pipeline state (SMLD, drift model, connected SMLD, etc.) using JLD2:

```julia
# After running full pipeline
(result, info) = analyze(image_stacks, config)
save_pipeline_state("output/pipeline.jld2", result)

# Resume later
state = load_pipeline_state("output/pipeline.jld2")
state.smld              # Final SMLD
state.drift_model       # Drift model
state.smld_connected    # Connected SMLD
```

## Verbosity Levels

Control output detail via the `verbose` field:

| Level | Constant | Output |
|-------|----------|--------|
| 0 | `Verbosity.SILENT` | Errors only |
| 1 | `Verbosity.PROGRESS` | Step names, counts, timing |
| 2 | `Verbosity.STANDARD` | + stats.md, basic figures (fit_quality, overlays, drift plots) |
| 3 | `Verbosity.DETAILED` | + diagnostic plots, per-filter breakdowns, localizations_per_frame |
| 4 | `Verbosity.DEBUG` | + MP4 animations, frame-by-frame analysis |

Set via `AnalysisConfig` or step function keyword:

```julia
config = AnalysisConfig(camera=cam, steps=[...], verbose=Verbosity.DETAILED)

# Or per-step
(smld, info) = analyze(data, DetectFitConfig(camera=cam, boxer=BoxerConfig(boxsize=9)); verbose=Verbosity.DEBUG)
```

## Uncertainty Calibration

Uncertainty calibration compares reported CRLB uncertainties against observed frame-to-frame scatter from linked emitters. It is **not** a separate pipeline step; it is enabled by setting the `calibration=` field on `FrameConnectConfig`. When set, frame connection runs as **link → calibrate → combine**, so the per-track combine uses the corrected weights in a single pass:

```julia
FrameConnectConfig(
    max_frame_gap = 5,
    calibration = CalibrationConfig(clamp_k_to_one=true),
)
```

Both `FrameConnectConfig` and `CalibrationConfig` are re-exported from SMLMFrameConnection via const aliases. `CalibrationConfig` and `CalibrationResult` are not standalone steps — they participate only via the `FrameConnectConfig.calibration` field, and results are surfaced on `FrameConnectInfo.calibration`.

`FrameConnectConfig` also accepts a `track_length=(min, max)` filter (localizations per track) that drops tracks whose linked-localization count falls outside the range — e.g. `(2.0, Inf)` removes single-frame blinks, while a finite upper bound drops sticky/fiducial tracks. The number of tracks removed is reported as `FrameConnectInfo.n_filtered`.

### Model

The calibration fits: `observed_variance = A + B * CRLB_variance`

Where:
- **A** (nm^2): Additive motion/vibration variance (`sigma_motion = sqrt(A)`)
- **B** (dimensionless): CRLB scale factor (`k = sqrt(B)`)

### Corrected uncertainties

After calibration, each localization's uncertainty is corrected:

```
sigma_corrected = sqrt(sigma_motion^2 + k^2 * sigma_CRLB^2)
```

Track recombination then uses weighted averaging with corrected uncertainties, giving properly calibrated combined positions.

### Configuration

```julia
FrameConnectConfig(
    max_frame_gap = 5,
    calibration = CalibrationConfig(
        clamp_k_to_one = true,        # k >= 1 (CRLB is theoretical lower bound)
        filter_high_chi2 = false,     # Optional: remove tracks with high chi-squared pairs
        chi2_filter_threshold = 6.0,
    ),
)
```

### Accessing calibration results

```julia
(smld, fc_info) = analyze(smld, FrameConnectConfig(
    max_frame_gap=5, calibration=CalibrationConfig()))

cal = fc_info.info.calibration   # CalibrationResult, or nothing if disabled
cal.k_scale                       # k = sqrt(B)
cal.sigma_motion_nm               # sqrt(A) in nm
cal.mean_chi2                     # ~2.0 means well-calibrated
cal.r_squared                     # WLS fit R²
```

### Interpreting results

- **Mean chi-squared ~ 2.0**: Uncertainties well-calibrated (expected for chi-squared(2) distribution)
- **Mean chi-squared > 2.5**: Uncertainties underestimated (motion, vibration, double-emitter fits)
- **Mean chi-squared < 1.5**: Uncertainties overestimated (conservative CRLB)

## I/O Formats

### Native SMLD format (HDF5)

Save and load localization results with full metadata:

```julia
save_smld("results.h5", smld; source_file="data.h5", drift_model=dm)
smld = load_smld("results.h5")
smld_info("results.h5")  # Quick summary without loading data
```

The HDF5 file stores emitter positions, uncertainties, camera calibration, drift model coefficients, and provenance information. Supports all emitter types including GaussMLE types with fitted PSF widths.

### SMART microscope H5

Import data from SMART microscope acquisitions:

```julia
images = load_smart_h5("acquisition.h5")          # returns the image array
info = load_smart_h5_info("acquisition.h5")        # metadata only
images, info = smart_h5_to_array("acquisition.h5") # array (rows, cols, frames) + metadata
```

### MIC format

Import data from MIC (MATLAB Instrument Control):

```julia
images, dataset_indices = load_mic_h5("experiment.h5")  # array + per-frame block index
info = load_mic_h5_info("experiment.h5")                # metadata only
block = load_mic_h5_block("experiment.h5", 1)           # single block (memory efficient)
```

Block-based loading is automatic when using `DetectFitConfig(path=..., h5_format=:mic)`.

### File-based detection

`DetectFitConfig` supports loading directly from H5 files, avoiding the need to hold all images in memory:

```julia
# Auto-detect format, MIC blocks become datasets
config = AnalysisConfig(
    camera = cam,
    steps = [DetectFitConfig(path="data.h5", h5_format=:mic, boxer=BoxerConfig(boxsize=9)), ...],
)
(result, info) = analyze(config)

# Multiple files (one per dataset)
config = AnalysisConfig(
    camera = cam,
    steps = [DetectFitConfig(paths=["d1.h5", "d2.h5"], boxer=BoxerConfig(boxsize=9)), ...],
)
(result, info) = analyze(config)
```

## Multi-Target Analysis

`MultiTargetConfig` orchestrates independent analysis pipelines for multiple color channels and produces composite multi-channel renders.

### Configuration

Each channel gets its own `AnalysisConfig` and image data. The `MultiTargetConfig` ties them together with labels, default colors, and an ordered list of multi-target steps that operate on the resulting `Vector{BasicSMLD}`:

```julia
mt = MultiTargetConfig(
    labels = [:IgG, :C1q],
    colors = [:cyan, :magenta],   # default: cyan/magenta for 2, CMY for 3
    steps = [
        CompositeRenderConfig(zoom=20.0, strategy=GaussianRender()),
        CrossAlignConfig(method=:entropy),
        CompositeRenderConfig(zoom=20.0, strategy=GaussianRender()),  # post-alignment render
        CrossCorrConfig(r_max=0.5, dr=0.005),
    ],
    outdir = "output/cell1/",
)

(result, info) = analyze([
    (image_stacks_647, config_647),
    (image_stacks_568, config_568),
], mt)
```

Multi-target steps are dispatched on `analyze(smlds::Vector{BasicSMLD}, cfg::AbstractMultiTargetStep)`:

- **`CompositeRenderConfig`** — render multi-channel composite image (pass-through; SMLDs not modified).
- **`CrossAlignConfig`** — align channels via entropy or FFT cross-correlation (state-modifying; returns aligned SMLDs).
- **`CrossCorrConfig`** — pair cross-correlation g(r).

### Accessing results

```julia
result.smlds              # Vector{BasicSMLD}, one per channel (may be aligned)
result[:IgG]              # Per-channel AnalysisResult
result[:IgG].smld         # Channel SMLD
info.channels[:IgG]       # Per-channel AnalysisInfo
info.step_infos           # Vector{StepInfo} for the multi-target steps
info.steps                # Dict{Symbol, AbstractSMLMInfo} keyed by step name
```

### Output structure

Multi-target analysis writes per-channel outputs and composite renders:

```
output/cell1/
├── IgG/                  # Per-channel pipeline output
│   ├── 01_detectfit/
│   ├── 02_filter/
│   └── ...
├── C1q/
│   └── ...
├── composite/            # Multi-channel overlay renders + alignment / crosscorr outputs
│   ├── 01_compositerender/
│   ├── 02_crossalign/
│   ├── 03_compositerender/
│   ├── 04_crosscorr/
│   └── README.md
├── smld_IgG.h5           # Per-channel SMLD files
├── smld_C1q.h5
└── multi_target_config.toml
```

### Composite rendering

`CompositeRenderConfig` accepts the same render strategies as `RenderConfig` (`GaussianRender`, `HistogramRender`, `CircleRender`, ...). Per-step `colors` override the `MultiTargetConfig.colors` defaults. For non-histogram strategies, intensities are clipped (default 99th percentile) before per-channel normalization to improve contrast. Histogram overlays use saturate mode (count=1 = full brightness).
