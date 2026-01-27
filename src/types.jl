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
struct DataSource
    images::Union{AbstractArray{<:Real,3}, Nothing}
    path::Union{String, Nothing}
    frame_range::Union{UnitRange{Int}, Nothing}
end

DataSource(images::AbstractArray{<:Real,3}) = DataSource(images, nothing, nothing)
DataSource(path::String; frame_range=nothing) = DataSource(nothing, path, frame_range)

function get_images(ds::DataSource)
    ds.images !== nothing && return ds.images
    ds.path !== nothing || error("No data source specified")
    data, _ = smart_h5_to_array(ds.path; max_frames=ds.frame_range === nothing ? nothing : last(ds.frame_range))
    if ds.frame_range !== nothing
        return data[:, :, ds.frame_range]
    end
    data
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
# Step Configs - abstract type and common interface
# ============================================================
abstract type StepConfig end

# Get the name field from any config
step_name(cfg::StepConfig) = cfg.name
step_verbose(cfg::StepConfig) = cfg.verbose

# ============================================================
# Step Record - what gets logged after each step
# ============================================================
struct StepRecord
    number::Int
    name::String
    config::StepConfig
    timestamp::DateTime
    timing::Float64
    summary::Dict{Symbol, Any}
end

function StepRecord(number::Int, cfg::StepConfig, timing::Float64, summary::Dict{Symbol,Any})
    StepRecord(number, step_name(cfg), cfg, now(), timing, summary)
end

# ============================================================
# Analysis Checkpoint - snapshot of state for reset
# ============================================================
struct AnalysisCheckpoint
    roi_batch::Union{SMLMData.ROIBatch, Nothing}
    roi_datasets::Union{Vector{Int}, Nothing}  # Dataset index for each ROI
    smld_raw::Union{SMLMData.BasicSMLD, Nothing}
    smld::Union{SMLMData.BasicSMLD, Nothing}
    smld_connected::Union{SMLMData.BasicSMLD, Nothing}
    drift_model::Any
    bagol_result::Any  # BaGoLDiagnostics from SMLMBaGoL
    bagol_smld::Union{SMLMData.BasicSMLD, Nothing}  # Grouped emitters from BaGoL
end

# ============================================================
# Analysis - main state container
# ============================================================
mutable struct Analysis
    # Input (set once)
    data::DataSource
    camera::SMLMData.AbstractCamera

    # Multi-dataset info
    n_datasets::Int
    n_frames_per_dataset::Int

    # Pipeline products (updated by steps)
    roi_batch::Union{SMLMData.ROIBatch, Nothing}
    roi_datasets::Union{Vector{Int}, Nothing}  # Dataset index for each ROI
    smld_raw::Union{SMLMData.BasicSMLD, Nothing}
    smld::Union{SMLMData.BasicSMLD, Nothing}
    smld_connected::Union{SMLMData.BasicSMLD, Nothing}
    drift_model::Any
    bagol_result::Any  # BaGoLDiagnostics from SMLMBaGoL
    bagol_smld::Union{SMLMData.BasicSMLD, Nothing}  # Grouped emitters from BaGoL

    # Checkpoints for reset
    checkpoints::Dict{Int, AnalysisCheckpoint}

    # History
    steps::Vector{StepRecord}
    step_counter::Int

    # Output control
    outdir::Union{String, Nothing}
    verbose::Int
    checkpoint::Bool  # Persist checkpoints to outdir/checkpoints/
end

"""
    Analysis(data, camera; roi=nothing, n_datasets=1, outdir=nothing, verbose=Verbosity.STANDARD, checkpoint=false)

Create an Analysis from image data and camera.

# Arguments
- `data`: Image array (y, x, frames) or DataSource
- `camera`: AbstractCamera defining pixel geometry

# Keywords
- `roi`: Optional ROI as NamedTuple `(x=100:300, y=50:200)` to crop images/camera.
         Physical coordinates are preserved (cropped camera retains original μm positions).
- `n_datasets`: Number of datasets (for multi-dataset acquisitions)
- `outdir`: Output directory for results
- `verbose`: Verbosity level (default: Verbosity.STANDARD)
- `checkpoint`: Enable disk persistence of checkpoints

# Example
```julia
# Full FOV
a = Analysis(images, camera)

# Crop to ROI while preserving coordinates
a = Analysis(images, camera; roi=(x=100:300, y=50:200))
```
"""
function Analysis(data, camera::SMLMData.AbstractCamera; roi=nothing, n_datasets=1, outdir=nothing, verbose=Verbosity.STANDARD, checkpoint=false)
    ds = data isa DataSource ? data : DataSource(data)
    images = get_images(ds)

    # Apply ROI cropping if specified
    if roi !== nothing
        roi_x = roi.x
        roi_y = roi.y
        images = crop_images(images, roi_x, roi_y)
        camera = crop_camera(camera, roi_x, roi_y)
        # Wrap cropped images as new DataSource
        ds = DataSource(images)
    end

    # Compute n_frames_per_dataset from total frames
    total_frames = size(images, 3)
    n_frames_per_dataset = div(total_frames, n_datasets)
    if n_frames_per_dataset * n_datasets != total_frames
        @warn "Total frames ($total_frames) not evenly divisible by n_datasets ($n_datasets)"
    end
    Analysis(
        ds,
        camera,
        n_datasets,
        n_frames_per_dataset,
        nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing,  # roi_batch, roi_datasets, smld_raw, smld, smld_connected, drift_model, bagol_result, bagol_smld
        Dict{Int, AnalysisCheckpoint}(),
        StepRecord[],
        0,
        outdir,
        verbose,
        checkpoint
    )
end

# Convenience: Analysis from path
function Analysis(path::String, camera::SMLMData.AbstractCamera; frame_range=nothing, roi=nothing, n_datasets=1, outdir=nothing, verbose=Verbosity.STANDARD, checkpoint=false)
    Analysis(DataSource(path; frame_range), camera; roi, n_datasets, outdir, verbose, checkpoint)
end

# Convenience: Analysis without data (for DetectFitConfig workflow where data comes from config)
function Analysis(camera::SMLMData.AbstractCamera; outdir=nothing, verbose=Verbosity.STANDARD, checkpoint=false)
    Analysis(
        DataSource(nothing, nothing, nothing),  # Empty data source
        camera,
        1,   # n_datasets - will be set by DetectFitConfig
        0,   # n_frames_per_dataset - will be set by DetectFitConfig
        nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing,  # roi_batch, roi_datasets, smld_raw, smld, smld_connected, drift_model, bagol_result, bagol_smld
        Dict{Int, AnalysisCheckpoint}(),
        StepRecord[],
        0,
        outdir,
        verbose,
        checkpoint
    )
end

# ============================================================
# Pretty printing
# ============================================================
function Base.show(io::IO, a::Analysis)
    n_steps = length(a.steps)
    n_locs = a.smld === nothing ? 0 : length(a.smld.emitters)
    print(io, "Analysis: $n_steps steps, $n_locs localizations")
end

function Base.show(io::IO, ::MIME"text/plain", a::Analysis)
    println(io, "Analysis: $(length(a.steps)) steps")
    println(io)

    if isempty(a.steps)
        println(io, "  (no steps run yet)")
    else
        for s in a.steps
            cp_marker = haskey(a.checkpoints, s.number) ? " [checkpoint]" : ""
            println(io, "  $(s.number). $(s.name)$cp_marker ($(round(s.timing, digits=2))s)")
            for (k, v) in s.summary
                println(io, "      $k: $v")
            end
        end
    end

    println(io)
    if !isempty(a.checkpoints)
        println(io, "Checkpoints (memory): $(sort(collect(keys(a.checkpoints))))")
    end
    if a.smld !== nothing
        println(io, "Current: $(length(a.smld.emitters)) localizations")
    end
    if a.outdir !== nothing
        println(io, "Output: $(a.outdir)")
        if a.checkpoint
            println(io, "Checkpoint persistence: enabled")
        end
    end
end

function Base.show(io::IO, cfg::StepConfig)
    T = typeof(cfg)
    fields = fieldnames(T)
    vals = [string(f, "=", getfield(cfg, f)) for f in fields if f != :verbose]
    print(io, "$(nameof(T))($(join(vals, ", ")))")
end

function Base.show(io::IO, r::StepRecord)
    print(io, "Step $(r.number): $(r.name) ($(round(r.timing, digits=2))s)")
end
