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

function Analysis(data, camera::SMLMData.AbstractCamera; n_datasets=1, outdir=nothing, verbose=Verbosity.STANDARD, checkpoint=false)
    ds = data isa DataSource ? data : DataSource(data)
    # Compute n_frames_per_dataset from total frames
    images = get_images(ds)
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
        nothing, nothing, nothing, nothing, nothing, nothing,  # roi_batch, roi_datasets, smld_raw, smld, smld_connected, drift_model
        Dict{Int, AnalysisCheckpoint}(),
        StepRecord[],
        0,
        outdir,
        verbose,
        checkpoint
    )
end

# Convenience: Analysis from path
function Analysis(path::String, camera::SMLMData.AbstractCamera; frame_range=nothing, n_datasets=1, outdir=nothing, verbose=Verbosity.STANDARD, checkpoint=false)
    Analysis(DataSource(path; frame_range), camera; n_datasets, outdir, verbose, checkpoint)
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
