"""
Core types for the SMLMAnalysis pipeline.
"""

using Dates

# ============================================================
# Verbosity Levels
# ============================================================
module Verbosity
    const SILENT = 0    # Errors only
    const PROGRESS = 1  # Step names + counts
    const STANDARD = 2  # + stats.md, basic figures
    const DETAILED = 3  # + diagnostic plots
    const DEBUG = 4     # + MP4, frame-by-frame, heavy viz
end

# ============================================================
# Data Source - lazy loading wrapper
# ============================================================

"""
    DataSource

Lazy loading wrapper for SMLM image data. Supports:
- Single 3D array (one dataset)
- Vector of 3D arrays (multiple datasets, boundaries encoded in data structure)
- File path for deferred loading

# Constructors
```julia
DataSource(images)                          # From single 3D array (1 dataset)
DataSource(image_stacks)                    # From Vector{Array} (N datasets)
DataSource(path; frame_range=nothing)       # From file path (lazy)
```
"""
struct DataSource
    images::Union{AbstractArray{<:Real,3}, Nothing}
    images_vec::Union{Vector{<:AbstractArray{<:Real,3}}, Nothing}
    path::Union{String, Nothing}
    frame_range::Union{UnitRange{Int}, Nothing}
end

DataSource(images::AbstractArray{<:Real,3}) = DataSource(images, nothing, nothing, nothing)
DataSource(vec::Vector{<:AbstractArray{<:Real,3}}) = DataSource(nothing, vec, nothing, nothing)
DataSource(path::String; frame_range=nothing) = DataSource(nothing, nothing, path, frame_range)
# Empty data source (for file-based DetectFitConfig workflows)
DataSource() = DataSource(nothing, nothing, nothing, nothing)

function get_images(ds::DataSource)
    ds.images !== nothing && return ds.images
    ds.images_vec !== nothing && error("DataSource holds multiple datasets. Access via ds.images_vec.")
    ds.path !== nothing || error("No data source specified")
    data, _ = smart_h5_to_array(ds.path; max_frames=ds.frame_range === nothing ? nothing : last(ds.frame_range))
    if ds.frame_range !== nothing
        return data[:, :, ds.frame_range]
    end
    data
end

"""
    n_datasets(ds::DataSource) -> Int

Number of datasets in this data source.
"""
function n_datasets(ds::DataSource)
    ds.images_vec !== nothing && return length(ds.images_vec)
    ds.images !== nothing && return 1
    1  # file-based: determined at load time
end

"""
    n_frames_per_dataset(ds::DataSource) -> Int

Frames per dataset. For Vector{Array}, uses first element.
"""
function n_frames_per_dataset(ds::DataSource)
    ds.images_vec !== nothing && return size(ds.images_vec[1], 3)
    ds.images !== nothing && return size(ds.images, 3)
    0  # file-based: determined at load time
end

# ============================================================
# ROI Cropping - preserves physical coordinates
# ============================================================
"""
    crop_camera(camera::AbstractCamera, roi_x::UnitRange, roi_y::UnitRange)

Crop camera to specified pixel ROI while preserving physical coordinates.

Returns a new camera with:
- pixel_edges_x/y sliced to cover ROI (preserves original μm positions)
- For SCMOSCamera: offset/gain/readnoise/qe matrices cropped accordingly

# Arguments
- `camera`: Original camera
- `roi_x`: Pixel range in x (columns), e.g., 100:300
- `roi_y`: Pixel range in y (rows), e.g., 50:200

# Example
```julia
cam = IdealCamera(512, 512, 0.1)
cam_cropped = crop_camera(cam, 100:300, 50:200)
# cam_cropped.pixel_edges_x[1] = 10.0 μm (not 0!)
```
"""
function crop_camera(camera::SMLMData.IdealCamera, roi_x::UnitRange{Int}, roi_y::UnitRange{Int})
    # Slice pixel edges (need +1 for edges)
    new_edges_x = camera.pixel_edges_x[roi_x.start:roi_x.stop+1]
    new_edges_y = camera.pixel_edges_y[roi_y.start:roi_y.stop+1]
    SMLMData.IdealCamera(new_edges_x, new_edges_y)
end

function crop_camera(camera::SMLMData.SCMOSCamera, roi_x::UnitRange{Int}, roi_y::UnitRange{Int})
    # Slice pixel edges (need +1 for edges)
    new_edges_x = camera.pixel_edges_x[roi_x.start:roi_x.stop+1]
    new_edges_y = camera.pixel_edges_y[roi_y.start:roi_y.stop+1]

    # Helper to crop matrix or keep scalar
    function crop_param(param, roi_y, roi_x)
        param isa Matrix ? param[roi_y, roi_x] : param
    end

    SMLMData.SCMOSCamera(
        new_edges_x,
        new_edges_y,
        crop_param(camera.offset, roi_y, roi_x),
        crop_param(camera.gain, roi_y, roi_x),
        crop_param(camera.readnoise, roi_y, roi_x),
        crop_param(camera.qe, roi_y, roi_x)
    )
end

"""
    crop_images(images::AbstractArray{T,3}, roi_x::UnitRange, roi_y::UnitRange) where T

Crop image stack to specified pixel ROI.

# Arguments
- `images`: 3D array (y, x, frames) - Julia convention: [row, col]
- `roi_x`: Pixel range in x (columns)
- `roi_y`: Pixel range in y (rows)
"""
function crop_images(images::AbstractArray{T,3}, roi_x::UnitRange{Int}, roi_y::UnitRange{Int}) where T
    images[roi_y, roi_x, :]
end

# ============================================================
# Step Configs - use AbstractSMLMConfig from SMLMData
# ============================================================

# Alias for backward compatibility within this package
const StepConfig = SMLMData.AbstractSMLMConfig

"""
    AbstractMultiTargetStep <: AbstractSMLMConfig

Abstract supertype for steps that operate on `Vector{BasicSMLD}` in the
multi-target pipeline (composite rendering, cross-channel alignment, etc.).
"""
abstract type AbstractMultiTargetStep <: SMLMData.AbstractSMLMConfig end

"""
    step_name(cfg::AbstractSMLMConfig) -> String

Derive step name from config type (e.g., `FilterConfig` → `"filter"`, `DriftConfig` → `"drift"`).
"""
step_name(cfg::SMLMData.AbstractSMLMConfig) = lowercase(replace(string(nameof(typeof(cfg))), r"Config|Options" => ""))

# ============================================================
# StepInfo - typed step record (replaces StepRecord)
# ============================================================
"""
    StepInfo <: AbstractSMLMInfo

Record of a completed pipeline step, stored in AnalysisInfo.

# Fields
- `number::Int`: Step number in the pipeline
- `name::String`: Step name (derived from config type, e.g. `"filter"`)
- `config::AbstractSMLMConfig`: The config used for this step
- `timestamp::DateTime`: When the step completed
- `elapsed_s::Float64`: Elapsed time in seconds
- `summary::Dict{Symbol, Any}`: Summary statistics (counts, acceptance rates, etc.)
- `info::Union{AbstractSMLMInfo, Nothing}`: Typed upstream info struct (DriftInfo, FrameConnectInfo, etc.)
"""
struct StepInfo <: SMLMData.AbstractSMLMInfo
    number::Int
    name::String
    config::SMLMData.AbstractSMLMConfig
    timestamp::DateTime
    elapsed_s::Float64
    summary::Dict{Symbol, Any}
    info::Union{SMLMData.AbstractSMLMInfo, Nothing}
end

function StepInfo(number::Int, cfg::SMLMData.AbstractSMLMConfig, elapsed_s::Float64, summary::Dict{Symbol,Any}; info=nothing)
    StepInfo(number, step_name(cfg), cfg, now(), elapsed_s, summary, info)
end

# ============================================================
# StepRecord - deprecated, kept for JLD2 checkpoint backward compat
# ============================================================
struct StepRecord
    number::Int
    name::String
    config::SMLMData.AbstractSMLMConfig
    timestamp::DateTime
    timing::Float64
    summary::Dict{Symbol, Any}
    info::Any
end

"""Convert legacy StepRecord to StepInfo."""
function StepInfo(r::StepRecord)
    StepInfo(r.number, r.name, r.config, r.timestamp, r.timing, r.summary, nothing)
end

# ============================================================
# Native info structs for SMLMAnalysis-owned steps
# ============================================================

"""
    DetectFitInfo <: AbstractSMLMInfo

Info from combined detection and fitting step.
"""
struct DetectFitInfo <: SMLMData.AbstractSMLMInfo
    boxes_info::Vector{Any}
    fit_info::Vector{Any}
    n_datasets::Int
    n_rois::Int
    n_fits::Int
    n_frames_per_dataset::Int
    elapsed_s::Float64
end

"""
    FilterInfo <: AbstractSMLMInfo

Info from quality filtering step.
"""
struct FilterInfo <: SMLMData.AbstractSMLMInfo
    n_before::Int
    n_after::Int
    elapsed_s::Float64
end

"""
    DensityFilterInfo <: AbstractSMLMInfo

Info from density-based filtering step.
"""
struct DensityFilterInfo <: SMLMData.AbstractSMLMInfo
    n_before::Int
    n_after::Int
    threshold::Int
    elapsed_s::Float64
end

"""
    CompositeRenderInfo <: AbstractSMLMInfo

Info from a multi-channel composite render step.
"""
struct CompositeRenderInfo <: SMLMData.AbstractSMLMInfo
    render_info::SMLMRender.RenderInfo
    strategy::Symbol
    zoom::Float64
    n_channels::Int
    elapsed_s::Float64
end

"""
    CrossAlignInfo <: AbstractSMLMInfo

Info from cross-channel alignment step.
"""
struct CrossAlignInfo <: SMLMData.AbstractSMLMInfo
    align_info::SMLMDriftCorrection.AlignInfo
    shifts::Vector{Vector{Float64}}
    max_shift_nm::Float64
    elapsed_s::Float64
end

# ============================================================
# AnalysisInfo - aggregated info from all steps (tuple-pattern)
# ============================================================
"""
    AnalysisInfo <: AbstractSMLMInfo

Aggregated metadata from all analysis steps, following the tuple-pattern.

# Fields
- `elapsed_s::Float64`: Total elapsed time in seconds
- `steps::Dict{Symbol, Any}`: Step name → upstream info struct mapping
- `step_infos::Vector{StepInfo}`: Full step history with timing and config
"""
struct AnalysisInfo <: SMLMData.AbstractSMLMInfo
    elapsed_s::Float64
    steps::Dict{Symbol, Any}
    step_infos::Vector{StepInfo}
end

AnalysisInfo() = AnalysisInfo(0.0, Dict{Symbol, Any}(), StepInfo[])

# ============================================================
# AnalysisResult - immutable result from analyze()
# ============================================================
"""
    AnalysisResult

Immutable result from `analyze()`. Replaces the old mutable `Analysis` struct.

# Fields
- `smld::BasicSMLD`: Final SMLD after all steps
- `smld_connected::Union{BasicSMLD, Nothing}`: Connected SMLD (if frameconnect was run)
- `drift_model::Any`: Drift model (if driftcorrect was run)

# Access
```julia
(result, info) = analyze(image_stacks, config)
result.smld               # Final SMLD
result.drift_model        # Drift model for plotting
result.smld_connected     # Connected SMLD for track analysis
info.steps[:driftcorrect] # Step info from upstream packages
```
"""
struct AnalysisResult
    smld::SMLMData.BasicSMLD
    smld_connected::Union{SMLMData.BasicSMLD, Nothing}
    drift_model::Any
end

function Base.show(io::IO, r::AnalysisResult)
    n = length(r.smld.emitters)
    print(io, "AnalysisResult: $n localizations")
    r.drift_model !== nothing && print(io, ", drift corrected")
end

# ============================================================
# AnalysisConfig - pipeline configuration (uniform interface)
# ============================================================
"""
    AnalysisConfig <: AbstractSMLMConfig

Complete description of an SMLM analysis pipeline.

The `steps` vector contains upstream package configs (BoxerConfig, GaussMLEConfig, etc.)
and SMLMAnalysis-specific configs (FilterConfig, DensityFilterConfig). The pipeline executes
steps in order.

# Fields
- `camera::AbstractCamera`: Camera model (required, no default)
- `steps::Vector{SMLMData.AbstractSMLMConfig}`: Ordered pipeline steps
- `roi::Union{NamedTuple, Nothing}`: Optional ROI as `(x=100:300, y=50:200)` to crop images/camera
- `outdir::Union{String, Nothing}`: Output directory for results
- `verbose::Int`: Verbosity level (default: STANDARD)

# Example
```julia
config = AnalysisConfig(
    camera = cam,
    steps = [
        DetectFitConfig(
            boxer=BoxerConfig(boxsize=9, psf_sigma=0.130),
            fitter=GaussMLEConfig(psf_model=GaussianXYNBS(), iterations=20)),
        FilterConfig(photons=(500.0, Inf)),
        DriftConfig(degree=2, dataset_mode=:registered),
        RenderConfig(zoom=20, colormap=:inferno),
    ],
    outdir = "output/",
)
(result, info) = analyze(image_stacks, config)
```
"""
@kwdef struct AnalysisConfig <: SMLMData.AbstractSMLMConfig
    camera::SMLMData.AbstractCamera
    steps::Vector{SMLMData.AbstractSMLMConfig} = SMLMData.AbstractSMLMConfig[]
    roi::Union{@NamedTuple{x::UnitRange{Int}, y::UnitRange{Int}}, Nothing} = nothing
    outdir::Union{String, Nothing} = nothing
    verbose::Int = Verbosity.STANDARD
end

# Varargs constructor: AnalysisConfig(step1, step2, ...; camera=cam, outdir="out/")
function AnalysisConfig(steps::SMLMData.AbstractSMLMConfig...; camera::SMLMData.AbstractCamera, kwargs...)
    AnalysisConfig(; camera=camera, steps=collect(SMLMData.AbstractSMLMConfig, steps), kwargs...)
end

# ============================================================
# Pretty printing
# ============================================================

function Base.show(io::IO, cfg::SMLMData.AbstractSMLMConfig)
    T = typeof(cfg)
    fields = fieldnames(T)
    vals = [string(f, "=", getfield(cfg, f)) for f in fields]
    print(io, "$(nameof(T))($(join(vals, ", ")))")
end

function Base.show(io::IO, si::StepInfo)
    print(io, "Step $(si.number): $(si.name) ($(round(si.elapsed_s, digits=2))s)")
end

function Base.show(io::IO, info::AnalysisInfo)
    n = length(info.step_infos)
    print(io, "AnalysisInfo: $n steps, $(round(info.elapsed_s, digits=2))s")
end

function Base.show(io::IO, ::MIME"text/plain", info::AnalysisInfo)
    println(io, "AnalysisInfo: $(length(info.step_infos)) steps, $(round(info.elapsed_s, digits=2))s")
    for s in info.step_infos
        println(io, "  $(s.number). $(s.name) ($(round(s.elapsed_s, digits=2))s)")
        for (k, v) in s.summary
            println(io, "      $k: $v")
        end
    end
end

# ============================================================
# Multi-Target Types
# ============================================================

"""
    _default_colors(n::Int) -> Vector{Symbol}

Default color palette for multi-target rendering.
2 channels: cyan/magenta. 3 channels: cyan/magenta/yellow.
Supports up to 6 channels; provide explicit colors for more.
"""
function _default_colors(n::Int)
    defaults = [:cyan, :magenta, :yellow, :red, :green, :blue]
    n <= length(defaults) || error("Provide explicit colors for >$(length(defaults)) channels")
    defaults[1:n]
end

"""
    MultiTargetConfig <: AbstractSMLMConfig

Configuration for multi-target (multi-color) SMLM analysis.

Each channel runs an independent `analyze(data, AnalysisConfig)` pipeline,
then cross-channel steps (composite rendering, alignment) are executed
via dispatch on `AbstractMultiTargetStep` configs.

# Fields
- `labels::Vector{Symbol}`: Channel labels (e.g., `[:IgG, :C1q]`)
- `colors::Vector{Symbol}`: Colors per channel (default: cyan/magenta for 2, CMY for 3)
- `steps::Vector{AbstractMultiTargetStep}`: Ordered multi-target steps (renders, alignment, etc.)
- `outdir::String`: Output directory
- `verbose::Int`: Verbosity level

# Example
```julia
mt = MultiTargetConfig(
    labels = [:IgG, :C1q],
    steps = [
        CompositeRenderConfig(zoom=20.0, strategy=GaussianRender()),
        CrossAlignConfig(method=:entropy),
        CompositeRenderConfig(zoom=20.0, strategy=GaussianRender()),
    ],
    outdir = "output/cell1/",
)
```
"""
@kwdef struct MultiTargetConfig <: SMLMData.AbstractSMLMConfig
    labels::Vector{Symbol}
    colors::Vector{Symbol} = _default_colors(length(labels))
    steps::Vector{AbstractMultiTargetStep} = AbstractMultiTargetStep[]
    outdir::String
    verbose::Int = Verbosity.STANDARD
end

"""
    MultiTargetResult

Result of a multi-target analysis. Holds per-channel `AnalysisResult` objects and
the final SMLD vectors for composite rendering.

# Fields
- `labels::Vector{Symbol}`: Channel labels in order
- `smlds::Vector{SMLMData.BasicSMLD}`: Per-channel SMLD results (may be aligned)
- `channels::Dict{Symbol, AnalysisResult}`: Per-channel results
- `step_infos::Vector{StepInfo}`: Multi-target step history
- `outdir::String`: Output directory

# Indexing
```julia
result[:IgG]         # Access per-channel AnalysisResult
keys(result)         # Channel labels
result.smlds         # Vector of all SMLDs
```
"""
struct MultiTargetResult
    labels::Vector{Symbol}
    smlds::Vector{SMLMData.BasicSMLD}
    channels::Dict{Symbol, AnalysisResult}
    step_infos::Vector{StepInfo}
    outdir::String
end

Base.getindex(mtr::MultiTargetResult, key::Symbol) = mtr.channels[key]
Base.keys(mtr::MultiTargetResult) = mtr.labels

"""
    MultiTargetInfo <: AbstractSMLMInfo

Aggregated metadata from a multi-target analysis.

# Fields
- `elapsed_s::Float64`: Total elapsed time in seconds
- `channels::Dict{Symbol, AnalysisInfo}`: Per-channel analysis info
- `step_infos::Vector{StepInfo}`: Multi-target step history (composite renders, alignment, etc.)
- `steps::Dict{Symbol, Any}`: Step name → info mapping for convenience
"""
struct MultiTargetInfo <: SMLMData.AbstractSMLMInfo
    elapsed_s::Float64
    channels::Dict{Symbol, AnalysisInfo}
    step_infos::Vector{StepInfo}
    steps::Dict{Symbol, Any}
end

function Base.show(io::IO, mtr::MultiTargetResult)
    n = sum(length(s.emitters) for s in mtr.smlds)
    print(io, "MultiTargetResult: $(length(mtr.labels)) channels, $n total localizations")
end

function Base.show(io::IO, ::MIME"text/plain", mtr::MultiTargetResult)
    println(io, "MultiTargetResult: $(length(mtr.labels)) channels")
    for (i, label) in enumerate(mtr.labels)
        println(io, "  $label: $(length(mtr.smlds[i].emitters)) localizations")
    end
    print(io, "Output: $(mtr.outdir)")
end
