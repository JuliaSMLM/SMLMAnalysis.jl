# SMLMAnalysis API Overview

AI-parseable API reference for SMLMAnalysis.jl.

## Design Concepts

**Dispatch-based pipeline**: The pipeline loop calls `analyze(state, config)` for each step. Routing is pure Julia method dispatch on `(state_type, config_type)` -- no step registry, no `isa` checks. Adding a new step requires only a new config type and `analyze()` method; the orchestrator needs no changes.

**Tuple-pattern returns**: Every `analyze()` call returns `(result, info)`. Pipeline-level: `(AnalysisResult, AnalysisInfo)`. Step-level: `(smld_or_image, StepInfo)`.

**Two-layer steps**: Each step has an internal function that does the work (e.g., `filter_step()` → `(filtered, FilterInfo)`) and a thin `analyze()` wrapper that times it and creates a `StepInfo`. The `StepInfo` wraps the typed info with step number, timing, config, and summary dict.

**Composability**: Steps can be reordered, repeated, or omitted freely after `DetectFitConfig`. Wrong ordering gives a `MethodError`.

**Config provenance**: Some configs are defined locally (`DetectFitConfig`, `FilterConfig`), others are re-exported from upstream packages via const aliases (`DriftConfig = SMLMDriftCorrection.DriftConfig`).

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

### Step Dispatch

Individual steps use `analyze()` with typed configs:

```julia
(smld, info)  = analyze(data, cfg::DetectFitConfig; outdir, step_number, verbose)
(smld, info)  = analyze(cfg::DetectFitConfig; ...)  # file-based (requires camera in config)
(smld, info)  = analyze(smld, cfg::FilterConfig; outdir, step_number, verbose)
(smld, info)  = analyze(smld, cfg::FrameConnectConfig; outdir, step_number, verbose)
(smld, info)  = analyze(smld, cfg::DriftConfig; outdir, step_number, verbose)
(smld, info)  = analyze(smld, cfg::DensityFilterConfig; outdir, step_number, verbose)
(smld, info)  = analyze(smld, cfg::IntensityFilterConfig; outdir, step_number, verbose)
(smld, info)  = analyze(smld, cfg::BaGoLConfig; outdir, step_number, verbose)
(smld, info)  = analyze(smld, cfg::RenderConfig; outdir, step_number, verbose)  # image written to outdir; smld passes through
```

Calibration is not a standalone step: set `FrameConnectConfig(calibration=CalibrationConfig(...))`
and it runs inside frame connection. Clustering, edge classification, and the
multi-target steps (`CompositeRenderConfig`, `CrossAlignConfig`, `CrossCorrConfig`)
dispatch on their own config types — see [Step Configs](#step-configs) below.

Each returns `(result, StepInfo)` where StepInfo wraps:
- Timing (`elapsed_s`), config, step number, summary dict
- Typed `info` field (e.g., `FilterInfo`, `DriftInfo`, `FrameConnectInfo`)

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
- `step_infos::Vector{StepInfo}` - Full step history

**Step info types:**
- `:detectfit` -> `(boxes=BoxesInfo, fit=FitInfo)`
- `:filter` -> `nothing`
- `:frameconnect` -> `FrameConnectInfo`
- `:driftcorrect` -> `DriftInfo`
- `:render` -> `Vector{RenderInfo}`

### StepInfo <: AbstractSMLMInfo

Logged after each step execution.

**Fields:**
- `number::Int` - Step number in the pipeline
- `name::String` - Step name (derived from config type, e.g., `"filter"`)
- `config::AbstractSMLMConfig` - Configuration used
- `timestamp::DateTime` - When the step completed
- `elapsed_s::Float64` - Duration in seconds
- `summary::Dict{Symbol, Any}` - Step statistics (counts, rates, etc.)
- `info::Union{AbstractSMLMInfo, Nothing}` - Typed upstream info struct (FilterInfo, DriftInfo, etc.)

**Constructor:**
```julia
StepInfo(number, cfg, elapsed_s, summary_dict; info=typed_info)
```

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
- `colors::Vector{Symbol}` - Colors per channel (default: cyan/magenta for 2, CMY for 3)
- `steps::Vector{AbstractMultiTargetStep}` - Multi-target steps run after the
  per-channel pipelines (`CompositeRenderConfig`, `CrossAlignConfig`, `CrossCorrConfig`)
- `outdir::String` - Output directory

### MultiTargetResult

Result of multi-channel analysis. Access per-channel results via `result[:label]`.

**Fields:**
- `labels::Vector{Symbol}` - Channel labels
- `smlds::Vector{BasicSMLD}` - Per-channel SMLDs
- `channels::Dict{Symbol, AnalysisResult}` - Per-channel results

## Step Configs

### DetectFitConfig

Combined detection + fitting step. Embeds native upstream configs (`BoxerConfig` from SMLMBoxer, `GaussMLEConfig` from GaussMLE).

```julia
DetectFitConfig(;
    boxer=BoxerConfig(boxsize=11, psf_sigma=0.135),  # Detection config
    fitter=GaussMLEConfig(psf_model=GaussianXYNBS(), iterations=20),  # Fitting config
    camera=nothing,           # Required for analyze() dispatch; injected by AnalysisConfig pipeline
    path=nothing,             # File-based loading
    paths=nothing,            # Multiple files (one per dataset)
    dataset_frames=nothing,   # Explicit frame ranges
    h5_format=:auto,          # :auto, :smart, :mic
)
```

**BoxerConfig** key fields: `boxsize` (Int, default 7), `psf_sigma` (Float64, microns), `min_photons` (Float64), `overlap` (Float64), `backend` (Symbol), `sigma_small`/`sigma_large` (advanced).

**GaussMLEConfig** key fields: `psf_model` (PSFModel), `iterations` (Int), `backend` (Symbol), `constraints` (ParameterConstraints), `batch_size` (Int).

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

### FrameConnectConfig (from SMLMFrameConnection)

Link localizations across frames.

```julia
FrameConnectConfig(;
    max_frame_gap=5,
    max_sigma_dist=5.0,
    n_density_neighbors=2,
    max_neighbors=2,
)
```

### CalibrationConfig (nested in FrameConnectConfig)

Calibrate localization uncertainties using frame-to-frame scatter. Not a standalone
step: pass it as `FrameConnectConfig(calibration=CalibrationConfig(...))` and it runs
inside frame connection (link → calibrate → combine). Results land in
`FrameConnectInfo.calibration`.

```julia
CalibrationConfig(;
    clamp_k_to_one=true,
    filter_high_chi2=false,
    chi2_filter_threshold=6.0,
)
```

### DriftConfig (from SMLMDriftCorrection)

Correct sample drift.

```julia
DriftConfig(;
    degree=2,
    dataset_mode=:registered,     # :registered or :continuous
    n_chunks=0,
    chunk_frames=0,
    maxn=200,
    quality=:singlepass,          # :singlepass, :iterative, or :fft
    auto_roi=false,
    max_iterations=10,
    convergence_tol=0.001,
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

### IntensityFilterConfig

Reject multi-emitter events via a Poisson upper-tail test against a spatially-varying
excitation field.

```julia
IntensityFilterConfig(;
    cutoff=0.01,                # p-value cutoff
    field_mode=:gaussian,       # :uniform or :gaussian (2D beam fit for λ(x,y))
    rate_percentile=0.95,       # percentile for single-emitter rate estimation
    estimate_p2=true,           # estimate the double-emitter fraction
)
```

### BaGoLConfig (from SMLMBaGoL)

Bayesian grouping of localizations (RJMCMC). Run after FrameConnect, before Render.
Fields pass directly to `SMLMBaGoL.run_bagol`.

```julia
BaGoLConfig(;
    μ=10.0,                     # expected locs per emitter
    shape=2.0,                  # NegBin shape (1=dSTORM, >1=DNA-PAINT)
    learn_distribution=true,    # true/false/:mu/:shape
    n_iterations=10000, burn_in=2000,
    partition_sigma=3.0,        # DBSCAN threshold in sigma units
    posterior_pixel_size=0.002, # 2 nm posterior image (0.0 disables)
)
```

### Clustering & spatial statistics (from SMLMClustering)

`analyze(smld, cfg)` where `cfg` is a cluster or statistics config. Config types are
re-exported: `DBSCANConfig`, `HDBSCANConfig`, `HierarchicalConfig`, `VoronoiConfig`
(clustering); `HopkinsConfig`, `VoronoiDensityConfig` (spatial-tendency statistics).

```julia
(smld, info) = analyze(smld, DBSCANConfig(eps_nm=50.0, min_points=5))
```

### Edge classification (from SMLMClustering)

Non-destructive interior/edge labelling of a cell's localizations. Config types
re-exported: `OuterPolygonConfig`, `KdeValleyConfig`. The class is written to
`info.class` / `interior_mask`; emitters are not removed.

```julia
(smld, info) = analyze(smld, OuterPolygonConfig())
```

### Multi-target steps

Dispatch on `Vector{BasicSMLD}` inside a `MultiTargetConfig` pipeline:

- `CompositeRenderConfig(; strategy, zoom, colors, clip_percentile=:auto, scalebar)` —
  multi-channel composite render (pass-through).
- `CrossAlignConfig(; method=:entropy, maxn, histbinsize)` — cross-channel alignment
  (state-modifying; returns aligned SMLDs).
- `CrossCorrConfig(; r_max, dr)` — pairwise cross-correlation g(r).

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

# MIC (MATLAB Instrument Control) format
images, metadata = load_mic_h5(path)
info = load_mic_h5_info(path)
block = load_mic_h5_block(path, block_index)
```

## AI Assistant Guide

### install_agent_guide

```julia
install_agent_guide(; tool=:claude, scope=:project, track=false,
                      overwrite=false, dir=pwd()) -> String
```

Writes a hierarchical, version-stamped ecosystem guide (the `analyze()` pipeline plus
every sub-package's `api_overview.md`, read from each `pkgdir`) for an AI coding
assistant. Returns the installed skill/bundle directory. Follows the lab
skills-installer convention (namespaced dir, `x-` provenance stamp, own-install
idempotent refresh).

- `tool`: `:claude` → Claude Code skill (`.claude/skills/smlma-ecosystem/SKILL.md` + `reference/*.md`); `:codex` → `smlm-agent-guide/` bundle + a managed block in `AGENTS.md`.
- `scope`: `:project` (into `dir`) or `:user` (into `~/.claude` / `~/.codex`).
- `track`: project scope only — when `false` (default) the guide is added to `.gitignore`; `track=true` commits it.
- `overwrite`: replace a target **not** installed by this installer (hand-made, or another package's). Re-running refreshes our own stamped install without it.

### uninstall_agent_guide / agent_guide_status

```julia
uninstall_agent_guide(; tool=:claude, scope=:project, dir=pwd()) -> Vector{String}
agent_guide_status(; tool=:claude, scope=:project, dir=pwd()) -> NamedTuple
```

`uninstall_agent_guide` removes a previously installed guide, but only targets carrying
this installer's provenance stamp (a hand-made skill or another package's install is
left untouched); returns the paths removed. `agent_guide_status` is a doctor: returns
`(; installed, path, source_version, source_commit, current_version, stale)`, with
`stale=true` when the installed version differs from the resolved one.
