# SMLMAnalysis API Overview

AI-parseable API reference for SMLMAnalysis.jl.

## Core Functions

### analyze

```julia
(result, info) = analyze(data, camera; kwargs...) -> (Analysis, AnalysisInfo)
```

Run complete SMLM analysis pipeline. Returns tuple of (result, info).

**Arguments:**
- `data`: Image stack (H×W×N array) or path to data file
- `camera`: `IdealCamera` or `SCMOSCamera`

**Keywords:**
- `outdir=nothing`: Output directory
- `verbose=Verbosity.STANDARD`: Verbosity level
- `n_datasets=1`: Number of datasets
- Detection: `boxsize=11`, `detect_min_photons=500.0`, `psf_sigma=0.135`, `use_gpu=true`
- Fitting: `psf_model=:variable`, `iterations=20`
- Filtering: `filter=true`, `min_photons=500.0`, `max_precision=0.007`, `min_pvalue=1e-3`
- Frame connection: `frameconnect=false`, `maxframegap=5`
- Drift: `drift=true`, `degree=2`
- Isolated: `isolated=false`, `n_sigma=2.0`
- Render: `render=true`, `render_zoom=20`

### run_step!

```julia
run_step!(a::Analysis, cfg::StepConfig) -> Analysis
```

Execute a pipeline step. Mutates `Analysis` in place.

### reset!

```julia
reset!(a::Analysis, step::Int) -> Analysis
reset!(a::Analysis) -> Analysis
```

Reset to checkpoint at step N, or to initial state.

### get_analysis_info

```julia
get_analysis_info(a::Analysis) -> AnalysisInfo
```

Extract AnalysisInfo from an Analysis object. Useful when running steps interactively with `run_step!` and wanting to get the aggregated info at the end.

## Types

### Analysis

Mutable state container for pipeline execution.

**Fields:**
- `data::DataSource` - Input images
- `camera::AbstractCamera` - Camera model
- `n_datasets::Int` - Number of datasets
- `n_frames_per_dataset::Int` - Frames per dataset
- `smld_raw::BasicSMLD` - Raw localizations (after detectfit)
- `smld::BasicSMLD` - Current localizations
- `smld_connected::BasicSMLD` - Connected localizations (after frameconnect)
- `drift_model` - Drift correction model
- `steps::Vector{StepRecord}` - Step history

**Constructor:**
```julia
Analysis(data, camera; roi=nothing, n_datasets=1, outdir=nothing,
         verbose=Verbosity.STANDARD, checkpoint=false)
```

### AnalysisInfo

Aggregated metadata from pipeline run.

**Fields:**
- `elapsed_ns::UInt64` - Total wall time in nanoseconds
- `steps::Dict{Symbol, Any}` - Per-step info keyed by step name

**Step info types:**
- `:detectfit` → `(boxes=BoxesInfo, fit=FitInfo)`
- `:filter` → `nothing`
- `:frameconnect` → `ConnectInfo`
- `:driftcorrect` → `DriftInfo`
- `:render` → `Vector{RenderInfo}`

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

## Step Configs

### DetectFitConfig

Combined detection + fitting step.

```julia
DetectFitConfig(;
    boxsize=11,
    overlap=2.0,
    min_photons=500.0,
    psf_sigma=0.135,
    use_gpu=true,
    psf_model=:variable,      # :fixed, :variable, :anisotropic
    psf_sigma_fit=0.135f0,    # For :fixed only
    iterations=20,
    path=nothing,             # File-based loading
    n_datasets=1,
    h5_format=:auto,          # :auto, :smart, :mic
    verbose=Verbosity.STANDARD
)
```

### FilterConfig

Filter localizations by quality metrics.

```julia
FilterConfig(;
    photons=(0.0, Inf),       # (min, max) photon range
    precision=(0.0, Inf),     # (min, max) precision range (μm)
    pvalue=(0.0, 1.0),        # (min, max) pvalue range
    psf_sigma=nothing,        # (min, max) PSF sigma range (μm)
    verbose=Verbosity.STANDARD
)
```

### FrameConnectConfig

Link localizations across frames.

```julia
FrameConnectConfig(;
    maxframegap=5,
    nsigmadev=5.0,
    nnearestclusters=2,
    nmaxnn=2,
    calibrate=true,
    clamp_k_to_one=true,
    filter_high_chi2=false,
    chi2_filter_threshold=6.0,
    verbose=Verbosity.STANDARD
)
```

### DriftCorrectConfig

Correct sample drift.

```julia
DriftCorrectConfig(;
    degree=2,
    continuous=false,         # true: TYPE 1 continuous, false: TYPE 2 registered
    n_chunks=0,
    chunk_frames=0,
    maxn=200,
    quality=:singlepass,      # :singlepass or :iterative
    warn_large_intershift=true,
    intershift_threshold_nm=500.0,
    verbose=Verbosity.STANDARD
)
```

### IsolatedConfig

Filter isolated emitters.

```julia
IsolatedConfig(;
    n_sigma=2.0,              # Neighbor search radius (σ units)
    verbose=Verbosity.STANDARD
)
```

### RenderConfig (from SMLMRender)

Generate super-resolution images. Uses SMLMRender.RenderConfig directly.
Each render is one step call with its own output folder.

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

### ConnectInfo (from SMLMFrameConnection)

```julia
struct ConnectInfo{T}
    connected::BasicSMLD{T}
    n_input::Int
    n_tracks::Int
    n_combined::Int
    k_on::Float64
    k_off::Float64
    k_bleach::Float64
    p_miss::Float64
    initialdensity::Vector{Float64}
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

### resume_analysis

```julia
a = resume_analysis(outdir::String; images=nothing) -> Analysis
```

Resume from disk checkpoint.

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
