# SMLMAnalysis API Overview

AI-parseable API reference for SMLMAnalysis.jl.

## Core Functions

### analyze

```julia
(result, info) = analyze(data, config::AnalysisConfig) -> (AnalysisResult, AnalysisInfo)
(result, info) = analyze(config::AnalysisConfig) -> (AnalysisResult, AnalysisInfo)  # file-based
(result, info) = analyze(data, steps...; camera, kwargs...) -> (AnalysisResult, AnalysisInfo)  # varargs
```

Run complete SMLM analysis pipeline. Returns tuple of (AnalysisResult, AnalysisInfo).

**Arguments:**
- `data`: Image data - `Vector{AbstractArray{<:Real,3}}` (multi-dataset) or `AbstractArray{<:Real,3}` (single dataset)
- `config`: `AnalysisConfig` with camera, steps, and output settings

For file-based workflows (MIC/SMART H5), use `analyze(config)` with `DetectFitConfig(path=...)`.

### Pure Step Functions

```julia
(smld, info) = detectfit(data, camera, cfg::DetectFitConfig; outdir, step_number, verbose)
(smld, info) = detectfit(camera, cfg::DetectFitConfig; ...)  # file-based
(smld, info) = filter_step(smld, cfg::FilterConfig; smld_raw, outdir, step_number, verbose)
(smld, info) = frameconnect_step(smld, cfg::FrameConnectConfig; outdir, step_number, verbose)
(smld, info) = driftcorrect_step(smld, cfg::DriftCorrectConfig; outdir, step_number, verbose)
(smld, info) = densityfilter_step(smld, cfg::DensityFilterConfig; outdir, step_number, verbose)
(image, info) = render_step(smld, cfg::RenderConfig; outdir, step_number, verbose)
```

Each step function returns `(result, NamedTuple)` where the NamedTuple includes:
- `step_record::StepRecord` - timing, config, summary stats
- Step-specific fields (e.g., `smld_raw` from detectfit, `smld_connected` from frameconnect, `drift_model` from driftcorrect)

## Types

### AnalysisResult

Immutable result from `analyze()`.

**Fields:**
- `smld::BasicSMLD` - Final localizations
- `smld_connected::Union{BasicSMLD, Nothing}` - Connected localizations (if frameconnect was run)
- `drift_model` - Drift correction model (if driftcorrect was run)

### AnalysisConfig

Complete pipeline description.

**Fields:**
- `camera::AbstractCamera` - Camera model (required)
- `steps::Vector{AbstractSMLMConfig}` - Ordered pipeline steps
- `roi::Union{NamedTuple, Nothing}` - Optional ROI as `(x=100:300, y=50:200)`
- `outdir::Union{String, Nothing}` - Output directory
- `verbose::Int` - Verbosity level (default: STANDARD)

### AnalysisInfo

Aggregated metadata from pipeline run.

**Fields:**
- `elapsed_s::Float64` - Total elapsed time in seconds
- `steps::Dict{Symbol, Any}` - Per-step info keyed by step name
- `step_records::Vector{StepRecord}` - Full step history

**Step info types:**
- `:detectfit` -> `(boxes=BoxesInfo, fit=FitInfo)`
- `:filter` -> `nothing`
- `:frameconnect` -> `FrameConnectInfo`
- `:driftcorrect` -> `DriftInfo`
- `:render` -> `Vector{RenderInfo}`

### StepRecord

Logged after each step execution.

**Fields:**
- `number::Int` - Step number
- `name::String` - Step name
- `config::StepConfig` - Configuration used
- `timestamp::DateTime` - Execution time
- `timing::Float64` - Duration in seconds
- `summary::Dict{Symbol, Any}` - Step statistics
- `info::Any` - Upstream package info

### DataSource

Lazy loading wrapper for image data.

**Fields:**
- `images::Union{AbstractArray{<:Real,3}, Nothing}` - Single dataset
- `images_vec::Union{Vector{<:AbstractArray{<:Real,3}}, Nothing}` - Multiple datasets
- `path::Union{String, Nothing}` - File path for deferred loading
- `frame_range::Union{UnitRange{Int}, Nothing}` - Frame subset

**Constructors:**
```julia
DataSource(images)           # Single 3D array (1 dataset)
DataSource(image_stacks)     # Vector{Array} (N datasets)
DataSource(path)             # File path (lazy loading)
DataSource()                 # Empty (file-based DetectFitConfig)
```

### MultiTargetConfig

Configuration for multi-channel analysis.

**Fields:**
- `labels::Vector{Symbol}` - Channel labels (e.g., `[:IgG, :C1q]`)
- `colors::Vector{Symbol}` - Colors per channel
- `render_zoom::Float64` - Zoom for composite renders
- `render_strategies::Vector{RenderingStrategy}` - Rendering strategies
- `outdir::String` - Output directory

### MultiTargetResult

Result of multi-channel analysis. Access per-channel results via `result[:label]`.

**Fields:**
- `labels::Vector{Symbol}` - Channel labels
- `smlds::Vector{BasicSMLD}` - Per-channel SMLDs
- `channels::Dict{Symbol, AnalysisResult}` - Per-channel results

## Step Configs

### DetectFitConfig

Combined detection + fitting step.

```julia
DetectFitConfig(;
    boxsize=11,
    overlap=2.0,
    min_photons=500.0,
    psf_sigma=0.135,
    backend=:auto,            # :auto, :gpu, :cpu
    psf_model=:variable,      # :fixed, :variable, :anisotropic
    psf_sigma_fit=0.135f0,    # For :fixed only
    iterations=20,
    path=nothing,             # File-based loading
    paths=nothing,            # Multiple files (one per dataset)
    dataset_frames=nothing,   # Explicit frame ranges
    h5_format=:auto,          # :auto, :smart, :mic
    filter_min_photons=500.0,
    filter_max_precision=0.007,
    filter_min_pvalue=1e-6,
)
```

Dataset boundaries are inferred from data structure (not a user integer):
- `Vector{Array}` data -> N datasets
- `path` with `:mic` format -> blocks auto-detected as datasets
- `paths` -> one file per dataset

### FilterConfig

Filter localizations by quality metrics.

```julia
FilterConfig(;
    photons=(0.0, Inf),       # (min, max) photon range
    precision=(0.0, Inf),     # (min, max) precision range (um)
    pvalue=(0.0, 1.0),        # (min, max) pvalue range
    psf_sigma=nothing,        # (min, max) PSF sigma range (um), or :auto
)
```

### FrameConnectConfig

Link localizations across frames.

```julia
FrameConnectConfig(;
    max_frame_gap=5,
    max_sigma_dist=5.0,
    n_density_neighbors=2,
    max_neighbors=2,
    calibrate=true,
    clamp_k_to_one=true,
    filter_high_chi2=false,
    chi2_filter_threshold=6.0,
)
```

### DriftCorrectConfig

Correct sample drift.

```julia
DriftCorrectConfig(;
    degree=2,
    continuous=false,         # true: continuous, false: registered
    n_chunks=0,
    chunk_frames=0,
    maxn=200,
    quality=:singlepass,      # :singlepass or :iterative
    warn_large_intershift=true,
    intershift_threshold_nm=500.0,
)
```

### DensityFilterConfig

Filter by local neighbor density.

```julia
DensityFilterConfig(;
    n_sigma=2.0,              # Neighbor search radius (sigma units)
    min_neighbors=:auto,      # :auto uses valley method
)
```

### RenderConfig (from SMLMRender)

Generate super-resolution images. Uses SMLMRender.RenderConfig directly.

```julia
RenderConfig(;
    strategy=GaussianRender(),  # GaussianRender(), HistogramRender(), CircleRender(), EllipseRender()
    zoom=nothing,               # Zoom factor for camera FOV mode
    pixel_size=nothing,         # Pixel size in nm (data bounds mode)
    colormap=nothing,           # :inferno, :viridis, :turbo, etc.
    color_by=nothing,           # nil=density, or :absolute_frame, :photons, :z
    clip_percentile=0.99,
    filename=nothing,           # Auto-generated if outdir set
)
```

## Info Types (Re-exported)

### FitInfo (from GaussMLE)

```julia
struct FitInfo
    elapsed_ns::UInt64
    backend::Symbol       # :cpu or :gpu
    device_id::Int        # GPU device, -1 for CPU
    n_fits::Int
    n_converged::Int
end
```

### BoxesInfo (from SMLMBoxer)

```julia
struct BoxesInfo
    backend::Symbol
    elapsed_ns::UInt64
    device_id::Int
end
```

### FrameConnectInfo (from SMLMFrameConnection)

```julia
struct FrameConnectInfo{T}
    connected::BasicSMLD{T}
    n_input::Int
    n_tracks::Int
    n_combined::Int
    k_on::Float64
    k_off::Float64
    k_bleach::Float64
    p_miss::Float64
    initial_density::Vector{Float64}
    elapsed_ns::UInt64
    algorithm::Symbol
    n_preclusters::Int
end
```

### DriftInfo (from SMLMDriftCorrection)

```julia
struct DriftInfo
    model::AbstractIntraInter
    elapsed_ns::UInt64
    backend::Symbol
    iterations::Int
    converged::Bool
    entropy::Float64
    history::Vector{Float64}
end
```

### RenderInfo (from SMLMRender)

```julia
struct RenderInfo
    elapsed_ns::UInt64
    backend::Symbol
    device_id::Int
    n_emitters_rendered::Int
    output_size::Tuple{Int,Int}
    pixel_size_nm::Float64
    strategy::Symbol
    color_mode::Symbol
    field_range::Union{Nothing, Tuple{Float64,Float64}}
end
```

### SimInfo (from SMLMSim)

```julia
struct SimInfo
    elapsed_ns::UInt64
    backend::Symbol
    device_id::Int
    seed::Union{UInt64, Nothing}
    smld_true::Union{SMLD, Nothing}
    smld_model::Union{SMLD, Nothing}
    n_patterns::Int
    n_emitters::Int
    n_localizations::Int
    n_frames::Int
end
```

### ImageInfo (from SMLMSim)

```julia
struct ImageInfo
    elapsed_ns::UInt64
    backend::Symbol
    device_id::Int
    frames_generated::Int
    n_photons_total::Float64
    output_size::Tuple{Int,Int,Int}
end
```

## Verbosity Levels

```julia
module Verbosity
    const SILENT = 0    # Errors only
    const PROGRESS = 1  # Step names + counts
    const STANDARD = 2  # + stats.md, basic figures
    const DETAILED = 3  # + diagnostic plots
    const DEBUG = 4     # + MP4, frame-by-frame
end
```

## I/O Functions

### save_smld / load_smld

```julia
save_smld(path::String, smld::BasicSMLD)
smld = load_smld(path::String) -> BasicSMLD
```

HDF5 serialization of SMLD data.

### save_pipeline_state / load_pipeline_state

```julia
save_pipeline_state(path, result::AnalysisResult; smld_raw, step_records, camera)
state = load_pipeline_state(path)  # Returns NamedTuple with smld, smld_raw, etc.
```

JLD2-based full pipeline state save/restore.

### H5 Loading

```julia
# SMART microscope format
data, info = load_smart_h5(path)
info = load_smart_h5_info(path)

# LidkeLab MIC format
images, metadata = load_lidkelab_h5(path)
info = load_lidkelab_h5_info(path)
block = load_lidkelab_h5_block(path, block_index)
```
