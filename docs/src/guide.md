# [Guide](@id Guide)

Conceptual reference for topics that need more explanation than a docstring provides.

## Pipeline Architecture

### The `steps` vector

`AnalysisConfig.steps` is a composable list of `AbstractSMLMConfig` subtypes. The orchestrator in `analysis.jl` iterates through them sequentially, threading state between steps:

```
DetectFitConfig (required first)
        │ produces smld
        ▼
  ┌─────────────────────────────┐
  │ FilterConfig       (0+ times)│
  │ FrameConnectConfig (0-1)     │  ← Any order,
  │ CalibrationConfig  (0-1)     │    any combination
  │ DriftConfig        (0-1)     │
  │ DensityFilterConfig(0+ times)│
  │ RenderConfig       (0+ times)│
  └─────────────────────────────┘
```

The only ordering constraint: `DetectFitConfig` must be first because it is the only step that produces a `BasicSMLD` from raw image data. All other steps receive and return an existing `smld`.

### Two-level dispatch

The pipeline operates at two levels:

1. **Orchestrator** (`analyze(data, AnalysisConfig)`): Loops over `steps`, calls each step, threads `smld` state, and collects `StepRecord`s into `AnalysisInfo`
2. **Step dispatch** (`analyze(smld, StepConfig)`): Each step config type has its own `analyze()` method that does the actual work

The step-by-step workflow (`analyze(smld, cfg)`) already supports any config in any order — the orchestrator just automates the loop and collects metadata.

### State threading

The orchestrator maintains `smld` as the working state. Each step receives the current `smld` and returns an updated one. Additional state captured from specific steps:

- `smld_raw` — from `DetectFitConfig` (pre-filter SMLD)
- `smld_connected` — from `FrameConnectConfig` (tracks with multi-frame info)
- `drift_model` — from `DriftConfig` (fitted drift polynomial)

### Repeatable and optional steps

- **`FilterConfig`**: Can appear multiple times (coarse filter early, tighter filter after connection)
- **`RenderConfig`**: Can appear multiple times (different zooms, colormaps, or strategies)
- **`DensityFilterConfig`**: Can appear multiple times
- **`FrameConnectConfig`**, **`CalibrationConfig`**, and **`DriftConfig`**: Typically used once, but not enforced

### Config provenance

| Config | Defined in | Notes |
|--------|-----------|-------|
| `DetectFitConfig` | SMLMAnalysis | Wraps SMLMBoxer + GaussMLE internally |
| `FilterConfig` | SMLMAnalysis | Pure SMLMAnalysis logic |
| `FrameConnectConfig` | SMLMFrameConnection | Re-exported via const alias |
| `CalibrationConfig` | SMLMAnalysis | Pure SMLMAnalysis logic |
| `DriftConfig` | SMLMDriftCorrection | Re-exported via const alias |
| `DensityFilterConfig` | SMLMAnalysis | Pure SMLMAnalysis logic |
| `RenderConfig` | SMLMRender | Re-exported via const alias |

SMLMAnalysis defines some step configs locally (`DetectFitConfig`, `FilterConfig`, `CalibrationConfig`, `DensityFilterConfig`) and re-exports others from upstream packages (`FrameConnectConfig`, `DriftConfig`, `RenderConfig`). Extending the pipeline with a new upstream package follows the same re-export pattern.

### Adding a new step

1. Define `@kwdef struct YourConfig <: AbstractSMLMConfig` in `src/steps/yourstep.jl`
2. Implement internal function `yourstep(smld, cfg; outdir, step_number, verbose)` returning `(result, info_namedtuple)`
3. Add dispatch: `analyze(smld::BasicSMLD, cfg::YourConfig; kw...) = yourstep(smld, cfg; kw...)`
4. Add `elseif cfg isa YourConfig` in **both** `analyze()` methods in `analysis.jl` (data-based and file-based)
5. Export `YourConfig` from `SMLMAnalysis.jl`

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

Uncertainty calibration is a separate pipeline step (`CalibrationConfig`) that compares reported CRLB uncertainties against observed frame-to-frame scatter from frame connection. It must follow `FrameConnectConfig` in the pipeline.

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
CalibrationConfig(
    clamp_k_to_one = true,        # k >= 1 (CRLB is theoretical lower bound)
    filter_high_chi2 = false,     # Optional: remove tracks with high chi-squared pairs
    chi2_filter_threshold = 6.0
)
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
images, info = load_smart_h5("acquisition.h5")
info = load_smart_h5_info("acquisition.h5")  # Metadata only
```

### MIC format

Import data from MIC (MATLAB Instrument Control):

```julia
images, metadata = load_mic_h5("experiment.h5")
info = load_mic_h5_info("experiment.h5")  # Metadata only
block = load_mic_h5_block("experiment.h5", 1)  # Single block (memory efficient)
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

Each channel gets its own `AnalysisConfig` and image data. The `MultiTargetConfig` ties them together with labels, colors, and composite rendering settings:

```julia
mt = MultiTargetConfig(
    labels = [:IgG, :C1q],
    colors = [:red, :green],
    render_zoom = 20,
    render_strategies = [GaussianRender(), CircleRender()],
    clip_percentile = 0.99,
    outdir = "output/cell1/",
)

(result, info) = analyze([
    (image_stacks_647, config_647),
    (image_stacks_568, config_568),
], mt)
```

### Accessing results

```julia
result.smlds              # Vector{BasicSMLD}, one per channel
result[:IgG]              # Per-channel AnalysisResult
result[:IgG].smld         # Channel SMLD
info.channels[:IgG]       # Per-channel AnalysisInfo
info.composite_renders    # Vector{RenderInfo} from composite renders
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
├── composite/            # Multi-channel overlay renders
│   ├── gaussianrender_20x.png
│   └── circlerender_20x.png
├── smld_IgG.h5           # Per-channel SMLD files
├── smld_C1q.h5
└── multi_target_config.toml
```

### Composite rendering

For non-histogram strategies, `clip_percentile` (default 0.99) clips outlier intensities before normalization to improve contrast. Histogram overlays use saturate mode instead (count=1 = full brightness).
