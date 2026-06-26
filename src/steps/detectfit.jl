"""
Combined Detect+Fit step as a pure function.

Supports three calling modes:
1. **Vector of image stacks**: `detectfit(data::Vector{<:AbstractArray{<:Real,3}}, camera, cfg)`
   Each element is a dataset. Primary path.
2. **Single image stack**: `detectfit(images::AbstractArray{<:Real,3}, camera, cfg)`
   Wraps into a 1-element vector.
3. **File-based**: `detectfit(camera, cfg)` with `cfg.path`/`cfg.paths` set.
   Loads each dataset from file, processes, frees memory.

For each dataset:
1. Detect ROIs via SMLMBoxer.getboxes
2. Fit localizations via GaussMLE.fit
3. Append to combined SMLD
"""

"""
    DetectFitConfig <: AbstractSMLMConfig

Combined detection and fitting step. Detects ROIs via SMLMBoxer and fits
localizations via GaussMLE in a single step, with per-dataset processing.

# Embedded Configs
- `boxer::BoxerConfig`: Detection parameters (boxsize, min_photons, psf_sigma, etc.)
- `fitter::GaussMLEConfig`: Fitting parameters (psf_model, iterations, constraints, etc.)

# Data Source Keywords (for file-based workflows)
- `path`: Single H5 file path
- `paths`: Vector of H5 file paths (one per dataset)
- `dataset_frames`: Explicit frame ranges per dataset
- `h5_format`: `:auto`, `:smart`, or `:mic`

# Dataset Selection
- `datasets`: Optional `AbstractVector{Int}` (e.g. `1:19` or `[1,2,3,5]`) selecting a
  subset of the resolved source slots. Applies uniformly across all source modes:
  MIC auto-blocks, multi-path `paths`, explicit `dataset_frames`, and in-memory
  `Vector{Array}` input. Selected slots are reindexed to contiguous `1:length(datasets)`
  in the output SMLD (original source indices preserved in `DetectFitInfo.selected_source_indices`).
  Default `nothing` uses all resolved slots.
"""
@kwdef struct DetectFitConfig <: SMLMData.AbstractSMLMConfig
    # Embedded upstream configs
    boxer::SMLMBoxer.BoxerConfig = SMLMBoxer.BoxerConfig(boxsize=11, psf_sigma=0.135)
    fitter::GaussMLEConfig = GaussMLEConfig(psf_model=GaussianXYNBS(), iterations=20)

    # Camera (optional - injected by AnalysisConfig pipeline, required for standalone analyze())
    camera::Union{SMLMData.AbstractCamera, Nothing} = nothing

    # Data source (for file-based workflows)
    path::Union{String, Nothing} = nothing
    paths::Union{Vector{String}, Nothing} = nothing
    dataset_frames::Union{Vector{UnitRange{Int}}, Nothing} = nothing

    # Dataset selection: subset of resolved source slots to include
    datasets::Union{AbstractVector{Int}, Nothing} = nothing

    # H5 format: :auto (detect), :smart (SMART microscope), :mic (MATLAB Instrument Control)
    h5_format::Symbol = :auto

    # Auto-camera from MIC H5: when set with camera=nothing and h5_format=:mic,
    # SCMOSCamera is built from H5 calibration at detectfit time
    pixel_size::Union{Float64, Nothing} = nothing
    qe::Float64 = 1.0

    # DEBUG-verbosity detection frame movies: realtime playback fps (= 1/exposure_s) for the
    # per-gallery-frame MP4s written under <step>/frame_movies/. nothing -> 20 fps fallback.
    movie_fps::Union{Float64, Nothing} = nothing
end

"""
    detectfit(data, camera, cfg; kwargs...) -> (smld, DetectFitInfo)

Run combined detection and fitting on image data.

# Arguments
- `data::Vector{<:AbstractArray{<:Real,3}}`: Vector of image stacks (one per dataset)
- `camera::AbstractCamera`: Camera model
- `cfg::DetectFitConfig`: Detection and fitting configuration

# Keywords
- `outdir::Union{String,Nothing}=nothing`: Output directory for results
- `step_number::Int=1`: Step number in the pipeline
- `verbose::Int=Verbosity.STANDARD`: Verbosity level

# Returns
`(smld::BasicSMLD, info::DetectFitInfo)`
"""
function detectfit(data::Vector{<:AbstractArray{<:Real,3}}, camera::SMLMData.AbstractCamera, cfg::DetectFitConfig;
                   outdir::Union{String,Nothing}=nothing,
                   step_number::Int=1,
                   verbose::Int=Verbosity.STANDARD,
                   checkpoint::Int=Checkpoint.EXPENSIVE)
    v = verbose
    dir = step_outdir(outdir, step_number, cfg)

    # Apply dataset selection (Option C: uniform across all source modes)
    selected_indices = cfg.datasets === nothing ? nothing : collect(cfg.datasets)
    data = _select_sources(data, cfg.datasets)
    n_datasets_val = length(data)

    v >= Verbosity.PROGRESS && @info "[$step_number] $(step_name(cfg))" n_datasets=n_datasets_val psf_model=typeof(cfg.fitter.psf_model)

    # Process each dataset
    all_emitters = AbstractEmitter[]
    n_frames_per_dataset = 0
    total_rois = 0
    total_fits = 0

    # Collect info from each dataset (tuple-pattern)
    all_boxes_info = []
    all_fit_info = []

    # Sample data for overlay plots (capture 12 frames spread across all datasets for 3x4 grid)
    n_sample_frames = 12
    sample_image_slices = []
    sample_roi_data = []
    sample_roi_x = Int[]
    sample_roi_y = Int[]
    sample_roi_frames = Int[]       # remapped to 1:n_total_samples for indexing sample_images
    sample_abs_frames = Int[]       # absolute frame labels for display
    samples_collected = 0
    sample_movie_stacks = Vector{Any}()   # DEBUG: ~100-frame raw window per gallery frame
    sample_movie_starts = Int[]           # absolute frame # of each window's first frame
    frame_offset = 0

    t = @elapsed begin
        for (ds, images) in enumerate(data)
            n_frames_ds = size(images, 3)
            n_frames_per_dataset = max(n_frames_per_dataset, n_frames_ds)

            v >= Verbosity.PROGRESS && @info "  Dataset $ds: $(size(images)) images"

            # Detect (tuple-pattern: returns (ROIBatch, BoxesInfo))
            (roi_batch, boxes_info) = SMLMBoxer.getboxes(images, camera, cfg.boxer)
            push!(all_boxes_info, boxes_info)
            n_rois = length(roi_batch)
            total_rois += n_rois

            # Fit (tuple-pattern: returns (BasicSMLD, FitInfo))
            (smld_ds, fit_info) = GaussMLE.fit(roi_batch, cfg.fitter)
            push!(all_fit_info, fit_info)
            n_fits = length(smld_ds.emitters)
            total_fits += n_fits

            v >= Verbosity.PROGRESS && @info "    $n_rois ROIs -> $n_fits fits (detect: $(round(boxes_info.elapsed_s, digits=2))s/$(boxes_info.backend), fit: $(round(fit_info.elapsed_s, digits=2))s/$(fit_info.backend))"

            # Capture sample data spread across all datasets for overlay plots
            if dir !== nothing && samples_collected < n_sample_frames
                remaining = n_sample_frames - samples_collected
                remaining_ds = n_datasets_val - ds + 1
                n_this = clamp(remaining ÷ remaining_ds, 1, min(n_frames_ds, remaining))
                idxs = n_this == 1 ? [cld(n_frames_ds, 2)] : [round(Int, x) for x in range(1, n_frames_ds, length=n_this)]

                for idx in idxs
                    push!(sample_image_slices, collect(images[:, :, idx]))
                    push!(sample_abs_frames, frame_offset + idx)
                    if v >= Verbosity.DEBUG
                        # 100-frame window kept in-block (shifted if idx near the block end);
                        # the gallery frame idx is always within it. Order matches sample_abs_frames.
                        i1 = clamp(idx, 1, max(1, n_frames_ds - 99))
                        push!(sample_movie_stacks, collect(images[:, :, i1:min(i1 + 99, n_frames_ds)]))
                        push!(sample_movie_starts, frame_offset + i1)
                    end
                end

                frame_to_sample = Dict(f => samples_collected + i for (i, f) in enumerate(idxs))
                for (ri, f) in enumerate(roi_batch.frame_indices)
                    if haskey(frame_to_sample, f)
                        push!(sample_roi_data, roi_batch.data[:, :, ri])
                        push!(sample_roi_x, roi_batch.x_corners[ri])
                        push!(sample_roi_y, roi_batch.y_corners[ri])
                        push!(sample_roi_frames, frame_to_sample[f])
                    end
                end
                samples_collected += n_this
            end
            frame_offset += n_frames_ds

            # Set dataset field and append
            for e in smld_ds.emitters
                push!(all_emitters, _with_dataset(e, ds))
            end
        end

        # Assemble sample data for overlay plots
        sample_images = nothing
        sample_roi_batch = nothing
        sample_original_frames = nothing
        if !isempty(sample_image_slices)
            sample_images = cat(sample_image_slices..., dims=3)
            sample_original_frames = sample_abs_frames
            if !isempty(sample_roi_data)
                sample_roi_batch = ROIBatch(cat(sample_roi_data..., dims=3),
                    sample_roi_x, sample_roi_y, sample_roi_frames, camera)
            end
        end
    end

    # Create combined SMLD
    if isempty(all_emitters)
        error("No localizations found across all datasets")
    end

    smld = BasicSMLD(all_emitters, camera, n_frames_per_dataset, n_datasets_val, Dict{String,Any}())

    detect_info = DetectFitInfo(all_boxes_info, all_fit_info,
        n_datasets_val, total_rois, total_fits, n_frames_per_dataset, t, selected_indices)

    if dir !== nothing
        _save_detectfit_outputs!(dir, outdir, smld, camera, cfg, v, t, total_rois, total_fits,
                                 n_datasets_val, n_frames_per_dataset,
                                 sample_images, sample_roi_batch, sample_original_frames,
                                 all_boxes_info, all_fit_info)

        if v >= Verbosity.DEBUG && !isempty(sample_movie_stacks)
            _write_detection_frame_movies(dir, sample_movie_stacks, sample_abs_frames, sample_movie_starts,
                                          something(cfg.movie_fps, 20.0), v)
        end

        if checkpoint >= Checkpoint.EXPENSIVE
            _save_step_smld(dir, smld; filename="smld_raw.jld2")
        end
    end

    v >= Verbosity.PROGRESS && @info "  -> $total_fits fits from $total_rois ROIs across $n_datasets_val datasets ($(round(t, digits=2))s)"

    return (smld, detect_info)
end

"""
    detectfit(images::AbstractArray{<:Real,3}, camera, cfg; kwargs...) -> (smld, info)

Convenience method for a single image stack. Wraps into a 1-element vector.
"""
function detectfit(images::AbstractArray{<:Real,3}, camera::SMLMData.AbstractCamera, cfg::DetectFitConfig; kwargs...)
    detectfit([images], camera, cfg; kwargs...)
end

"""
    detectfit(camera, cfg; kwargs...) -> (smld, info)

File-based detectfit. Loads image data from `cfg.path` or `cfg.paths`.
Each data source is loaded one at a time for memory efficiency.
"""
function detectfit(camera::SMLMData.AbstractCamera, cfg::DetectFitConfig;
                   outdir::Union{String,Nothing}=nothing,
                   step_number::Int=1,
                   verbose::Int=Verbosity.STANDARD,
                   checkpoint::Int=Checkpoint.EXPENSIVE)
    (cfg.path !== nothing || cfg.paths !== nothing) || error("File-based detectfit requires path or paths in config")
    sources = _resolve_file_sources(cfg)

    # Apply dataset selection (Option C: uniform across all source modes)
    selected_indices = cfg.datasets === nothing ? nothing : collect(cfg.datasets)
    sources = _select_sources(sources, cfg.datasets)

    v = verbose
    dir = step_outdir(outdir, step_number, cfg)
    n_datasets_val = length(sources)

    v >= Verbosity.PROGRESS && @info "[$step_number] $(step_name(cfg)) [file-based]" n_datasets=n_datasets_val psf_model=typeof(cfg.fitter.psf_model)

    # Process each dataset
    all_emitters = AbstractEmitter[]
    n_frames_per_dataset = 0
    total_rois = 0
    total_fits = 0

    all_boxes_info = []
    all_fit_info = []

    # Sample data for overlay plots (spread across all datasets)
    n_sample_frames = 12
    sample_image_slices = []
    sample_roi_data = []
    sample_roi_x = Int[]
    sample_roi_y = Int[]
    sample_roi_frames = Int[]
    sample_abs_frames = Int[]
    samples_collected = 0
    sample_movie_stacks = Vector{Any}()   # DEBUG: ~100-frame raw window per gallery frame
    sample_movie_starts = Int[]           # absolute frame # of each window's first frame
    frame_offset = 0

    t = @elapsed begin
        for (ds, source) in enumerate(sources)
            images = _load_source(source, v)
            n_frames_ds = size(images, 3)
            n_frames_per_dataset = max(n_frames_per_dataset, n_frames_ds)

            v >= Verbosity.PROGRESS && @info "  Dataset $ds: $(size(images)) images"

            (roi_batch, boxes_info) = SMLMBoxer.getboxes(images, camera, cfg.boxer)
            push!(all_boxes_info, boxes_info)
            n_rois = length(roi_batch)
            total_rois += n_rois

            (smld_ds, fit_info) = GaussMLE.fit(roi_batch, cfg.fitter)
            push!(all_fit_info, fit_info)
            n_fits = length(smld_ds.emitters)
            total_fits += n_fits

            v >= Verbosity.PROGRESS && @info "    $n_rois ROIs -> $n_fits fits (detect: $(round(boxes_info.elapsed_s, digits=2))s/$(boxes_info.backend), fit: $(round(fit_info.elapsed_s, digits=2))s/$(fit_info.backend))"

            # Capture sample data spread across all datasets for overlay plots
            if dir !== nothing && samples_collected < n_sample_frames
                remaining = n_sample_frames - samples_collected
                remaining_ds = n_datasets_val - ds + 1
                n_this = clamp(remaining ÷ remaining_ds, 1, min(n_frames_ds, remaining))
                idxs = n_this == 1 ? [cld(n_frames_ds, 2)] : [round(Int, x) for x in range(1, n_frames_ds, length=n_this)]

                for idx in idxs
                    push!(sample_image_slices, collect(images[:, :, idx]))
                    push!(sample_abs_frames, frame_offset + idx)
                    if v >= Verbosity.DEBUG
                        # 100-frame window kept in-block (shifted if idx near the block end);
                        # the gallery frame idx is always within it. Order matches sample_abs_frames.
                        i1 = clamp(idx, 1, max(1, n_frames_ds - 99))
                        push!(sample_movie_stacks, collect(images[:, :, i1:min(i1 + 99, n_frames_ds)]))
                        push!(sample_movie_starts, frame_offset + i1)
                    end
                end

                frame_to_sample = Dict(f => samples_collected + i for (i, f) in enumerate(idxs))
                for (ri, f) in enumerate(roi_batch.frame_indices)
                    if haskey(frame_to_sample, f)
                        push!(sample_roi_data, roi_batch.data[:, :, ri])
                        push!(sample_roi_x, roi_batch.x_corners[ri])
                        push!(sample_roi_y, roi_batch.y_corners[ri])
                        push!(sample_roi_frames, frame_to_sample[f])
                    end
                end
                samples_collected += n_this
            end
            frame_offset += n_frames_ds

            for e in smld_ds.emitters
                push!(all_emitters, _with_dataset(e, ds))
            end

            # Images freed when loop iteration ends (GC)
        end

        # Assemble sample data for overlay plots
        sample_images = nothing
        sample_roi_batch = nothing
        sample_original_frames = nothing
        if !isempty(sample_image_slices)
            sample_images = cat(sample_image_slices..., dims=3)
            sample_original_frames = sample_abs_frames
            if !isempty(sample_roi_data)
                sample_roi_batch = ROIBatch(cat(sample_roi_data..., dims=3),
                    sample_roi_x, sample_roi_y, sample_roi_frames, camera)
            end
        end
    end

    if isempty(all_emitters)
        error("No localizations found across all datasets")
    end

    smld = BasicSMLD(all_emitters, camera, n_frames_per_dataset, n_datasets_val, Dict{String,Any}())

    detect_info = DetectFitInfo(all_boxes_info, all_fit_info,
        n_datasets_val, total_rois, total_fits, n_frames_per_dataset, t, selected_indices)

    if dir !== nothing
        _save_detectfit_outputs!(dir, outdir, smld, camera, cfg, v, t, total_rois, total_fits,
                                 n_datasets_val, n_frames_per_dataset,
                                 sample_images, sample_roi_batch, sample_original_frames,
                                 all_boxes_info, all_fit_info)

        if v >= Verbosity.DEBUG && !isempty(sample_movie_stacks)
            _write_detection_frame_movies(dir, sample_movie_stacks, sample_abs_frames, sample_movie_starts,
                                          something(cfg.movie_fps, 20.0), v)
        end

        if checkpoint >= Checkpoint.EXPENSIVE
            _save_step_smld(dir, smld; filename="smld_raw.jld2")
        end
    end

    v >= Verbosity.PROGRESS && @info "  -> $total_fits fits from $total_rois ROIs across $n_datasets_val datasets ($(round(t, digits=2))s)"

    return (smld, detect_info)
end

# ============================================================
# Camera injection helper
# ============================================================

"""
    _inject_camera(cfg::DetectFitConfig, camera::AbstractCamera) -> DetectFitConfig

Inject camera into DetectFitConfig if not already set. Used by AnalysisConfig pipeline.
"""
function _inject_camera(cfg::DetectFitConfig, camera::SMLMData.AbstractCamera)
    cfg.camera !== nothing && return cfg
    DetectFitConfig(; camera=camera, [f => getfield(cfg, f) for f in fieldnames(DetectFitConfig) if f != :camera]...)
end

"""
    _auto_camera(cfg::DetectFitConfig) -> AbstractCamera

Resolve camera from DetectFitConfig. Uses cfg.camera if set, otherwise auto-builds
from MIC H5 calibration when pixel_size is provided.
"""
function _auto_camera(cfg::DetectFitConfig)
    cfg.camera !== nothing && return cfg.camera
    if cfg.pixel_size !== nothing
        h5_path = cfg.path !== nothing ? cfg.path :
                  cfg.paths !== nothing ? cfg.paths[1] :
                  error("Auto-camera requires path or paths in DetectFitConfig")
        format = cfg.h5_format == :auto ? _detect_h5_format(h5_path) : cfg.h5_format
        format == :mic || error("Auto-camera from H5 only supported for :mic format, got :$format")
        return build_camera_from_mic_h5(h5_path; pixel_size=cfg.pixel_size, qe=cfg.qe)
    end
    error("DetectFitConfig requires camera or pixel_size for auto-camera from MIC H5")
end

# ============================================================
# analyze() dispatch methods
# ============================================================

_step_summary(info::DetectFitInfo) = Dict{Symbol,Any}(
    :n_datasets => info.n_datasets,
    :n_rois => info.n_rois,
    :n_fits => info.n_fits,
    :n_frames_per_dataset => info.n_frames_per_dataset
)

"""
    analyze(data, cfg::DetectFitConfig; kwargs...) -> (smld, StepInfo)

Run combined detection and fitting. Camera must be set in `cfg.camera`.
"""
function analyze(data::Vector{<:AbstractArray{<:Real,3}}, cfg::DetectFitConfig;
                 outdir=nothing, step_number::Int=1, verbose::Int=Verbosity.STANDARD,
                 checkpoint::Int=Checkpoint.EXPENSIVE, kwargs...)
    camera = _auto_camera(cfg)
    t = @elapsed (smld, detect_info) = detectfit(data, camera, cfg;
        outdir=outdir, step_number=step_number, verbose=verbose, checkpoint=checkpoint)
    (smld, StepInfo(step_number, cfg, t, _step_summary(detect_info); info=detect_info))
end

function analyze(images::AbstractArray{<:Real,3}, cfg::DetectFitConfig; kwargs...)
    analyze([images], cfg; kwargs...)
end

"""File-based dispatch for pipeline use: `analyze(nothing, DetectFitConfig(path=...))`.
Routes to file-based `analyze(cfg::DetectFitConfig)` when no data is provided."""
analyze(::Nothing, cfg::DetectFitConfig; kwargs...) = analyze(cfg; kwargs...)

"""
    analyze(cfg::DetectFitConfig; kwargs...) -> (smld, StepInfo)

File-based detection and fitting. Requires `cfg.path` or `cfg.paths` and `cfg.camera`.
"""
function analyze(cfg::DetectFitConfig;
                 outdir=nothing, step_number::Int=1, verbose::Int=Verbosity.STANDARD,
                 checkpoint::Int=Checkpoint.EXPENSIVE, kwargs...)
    camera = _auto_camera(cfg)
    t = @elapsed (smld, detect_info) = detectfit(camera, cfg;
        outdir=outdir, step_number=step_number, verbose=verbose, checkpoint=checkpoint)
    (smld, StepInfo(step_number, cfg, t, _step_summary(detect_info); info=detect_info))
end

# ============================================================
# File source resolution
# ============================================================

"""Detect H5 file format by checking internal structure"""
function _detect_h5_format(filepath::String)
    HDF5.h5open(filepath, "r") do f
        if haskey(f, "Main/data")
            return :smart
        elseif haskey(f, "Channel01/Zposition001")
            return :mic
        else
            error("Unknown H5 format: $filepath")
        end
    end
end

"""
Resolve file-based sources from config options.

Returns a vector of NamedTuples describing each dataset to load.
Dataset count is determined from the data structure, not from a config field.
"""
function _resolve_file_sources(cfg::DetectFitConfig)
    # Multiple files: one per dataset
    if cfg.paths !== nothing
        format = cfg.h5_format == :auto ? _detect_h5_format(cfg.paths[1]) : cfg.h5_format
        return [(path=p, frame_range=nothing, format=format) for p in cfg.paths]
    end

    # Must have single path
    cfg.path === nothing && error("Must specify path or paths")

    # Detect format
    format = cfg.h5_format == :auto ? _detect_h5_format(cfg.path) : cfg.h5_format

    # Explicit frame ranges
    if cfg.dataset_frames !== nothing
        return [(path=cfg.path, frame_range=r, format=format) for r in cfg.dataset_frames]
    end

    # MIC format: auto-detect blocks, each block is a dataset
    if format == :mic
        info = load_mic_h5_info(cfg.path)
        return [(path=cfg.path, block=ds, format=format) for ds in 1:info.n_blocks]
    end

    # Single dataset (SMART format or MIC with single block)
    return [(path=cfg.path, frame_range=nothing, format=format)]
end

"""
Apply `cfg.datasets` selection to a resolved source list. Bounds-checks and
returns the selected subset in the order given. Pass-through when `cfg.datasets`
is `nothing`.
"""
function _select_sources(sources::AbstractVector, sel::Union{AbstractVector{Int}, Nothing})
    sel === nothing && return sources
    n = length(sources)
    for i in sel
        (1 <= i <= n) || error("DetectFitConfig.datasets contains index $i, valid range is 1:$n")
    end
    sources[sel]
end

"""Load images from a data source"""
function _load_source(source, v)
    if source.format == :smart
        if source.frame_range === nothing
            data, _ = smart_h5_to_array(source.path)
            return data
        else
            # Load specific frame range from SMART H5
            HDF5.h5open(source.path, "r") do file
                raw = file["Main/data"][:, :, source.frame_range]
                # Transpose to (height, width, frames)
                permutedims(raw, (2, 1, 3))
            end
        end
    elseif source.format == :mic
        if haskey(source, :block)
            # Block-based loading (efficient for MIC format)
            return load_mic_h5_block(source.path, source.block)
        elseif source.frame_range === nothing
            # Load all frames
            images, _ = load_mic_h5(source.path)
            return images
        else
            # Load specific frame range (less efficient, but needed for custom ranges)
            images, _ = load_mic_h5(source.path; max_frames=last(source.frame_range))
            return images[:, :, source.frame_range]
        end
    else
        error("Unknown H5 format: $(source.format)")
    end
end

# ============================================================
# Output functions
# ============================================================

"""
    _write_detection_frame_movies(dir, stacks, abs_frames, fps, v)

DEBUG-verbosity output. For each detection-overlay gallery frame, write a short realtime
MP4 of the raw camera frames. `stacks[k]` is a (H, W, ~100)-frame window containing gallery
frame `abs_frames[k]` (same order as the overlay panel). Contrast is ONE global linear
stretch over the whole 100-frame window: black point = 0.1st percentile (0.001), white point =
99.9th percentile (0.999) of all pixels in the window. Played at `fps` (realtime = 1/exposure).
A label (frame / fps / stretch percentiles / black-white values) is burned into each frame.
Files: `dir/frame_movies/detection_frame_<gallery_frame>.mp4`. Requires `ffmpeg` on PATH.
"""
function _write_detection_frame_movies(dir, stacks, abs_frames, start_frames, fps, v)
    ffmpeg = Sys.which("ffmpeg")
    if ffmpeg === nothing
        v >= Verbosity.PROGRESS && @warn "ffmpeg not found on PATH; skipping detection frame movies"
        return
    end
    mdir = joinpath(dir, "frame_movies"); mkpath(mdir)
    lo_pct, hi_pct = 0.001, 0.999             # global black / white percentile points
    font = ""                                  # monospace font for the burned-in label (optional)
    for fc in ("/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf",
               "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf")
        isfile(fc) && (font = fc; break)
    end
    for (k, stack) in enumerate(stacks)
        absf = abs_frames[k]
        h, w, nf = size(stack)
        flat = vec(Float32.(stack))
        vmin = Float32(quantile(flat, lo_pct))   # black point
        vmax = Float32(quantile(flat, hi_pct))   # white point
        vmax <= vmin && (vmax = vmin + 1f0)
        fn = joinpath(mdir, "detection_frame_$(absf).mp4")
        # burned-in label: per-frame ABSOLUTE frame # (eif: start + output-frame n) / realtime fps /
        # stretch percentiles / black-white values
        wstart = start_frames[k]
        label = "frame %{eif\\:$(wstart)+n\\:d}  $(round(Int,fps))fps  stretch $(lo_pct)-$(hi_pct)  black $(round(Int,vmin)) white $(round(Int,vmax))"
        vf = isempty(font) ? "format=yuv420p" :
             "format=yuv420p,drawtext=fontfile=$font:text='$label':x=10:y=10:fontsize=14:fontcolor=white:box=1:boxcolor=black@0.5:boxborderw=5"
        # raw gray frames piped to ffmpeg; -framerate sets realtime playback; even dims (500x500) ok for yuv420p
        cmd = `$ffmpeg -y -loglevel error -f rawvideo -pixel_format gray -video_size $(w)x$(h) -framerate $(fps) -i pipe:0 -an -vf $vf -c:v libx264 -crf 18 -movflags +faststart $fn`
        try
            open(cmd, "w") do io
                for fr in 1:nf
                    sl = @view stack[:, :, fr]
                    scaled = clamp.((Float32.(sl) .- vmin) ./ (vmax - vmin), 0f0, 1f0)
                    bytes = round.(UInt8, scaled .* 255f0)        # H×W
                    write(io, vec(permutedims(bytes)))            # row-major (W per row) for ffmpeg rawvideo
                end
            end
            v >= Verbosity.DETAILED && @info "    movie: $(basename(fn)) ($nf frames @ $(round(fps, digits=1)) fps)"
        catch err
            v >= Verbosity.PROGRESS && @warn "    failed to write $(basename(fn))" exception=err
        end
    end
    v >= Verbosity.PROGRESS && @info "  wrote $(length(stacks)) detection frame movie(s) -> $mdir"
end

function _save_detectfit_outputs!(dir, outdir, smld, camera, cfg, v, t, n_rois, n_fits,
                                  n_datasets, n_frames_per_dataset,
                                  sample_images, sample_roi_batch, sample_original_frames,
                                  all_boxes_info, all_fit_info)
    mkpath(dir)
    _save_config!(dir, cfg)

    # Write upstream info structs to info.toml
    # Write header first, then append sections
    open(joinpath(dir, "info.toml"), "w") do io
        println(io, "# Upstream package info")
    end
    n_ds = length(all_boxes_info)
    if n_ds == 1
        _save_info!(dir, all_boxes_info[1]; section="boxes_info")
        _save_info!(dir, all_fit_info[1]; section="fit_info")
    else
        for i in 1:n_ds
            _save_info!(dir, all_boxes_info[i]; section="boxes_info_$i")
            _save_info!(dir, all_fit_info[i]; section="fit_info_$i")
        end
    end

    if v >= Verbosity.STANDARD
        _write_detectfit_stats(dir, smld, cfg, t, n_rois, n_fits, n_datasets, n_frames_per_dataset)

        # Generate detection overlay and save sample cache for filter step
        if sample_images !== nothing && sample_roi_batch !== nothing
            box_colors = fill(:yellow, length(sample_roi_batch))
            _save_box_overlay(dir, "detection_overlay.png", sample_images, sample_roi_batch, box_colors;
                              title_prefix="Detection Frame", frame_labels=sample_original_frames)

            _save_detectfit_sample_cache(outdir, smld, sample_images, sample_roi_batch, sample_original_frames)
        end

        _save_loc_per_frame(dir, smld; title="Localizations per Frame (raw fits)")
    end
end

function _write_detectfit_stats(dir, smld, cfg, t, n_rois, n_fits, n_datasets, n_frames_per_dataset)
    emitters = smld.emitters
    n = length(emitters)

    photons = [e.photons for e in emitters]
    bg = [e.bg for e in emitters]
    σ_x = [e.σ_x for e in emitters]
    σ_y = [e.σ_y for e in emitters]
    pvalue = [e.pvalue for e in emitters]

    # Calculate photobleaching rate from localizations per frame (absolute frames)
    n_frames = smld.n_frames
    n_total = n_frames * smld.n_datasets
    frame_counts = zeros(Int, n_total)
    for e in emitters
        abs_frame = (e.dataset - 1) * n_frames + e.frame
        if abs_frame >= 1 && abs_frame <= n_total
            frame_counts[abs_frame] += 1
        end
    end
    bleach_result = _estimate_bleaching_rate(frame_counts)

    filepath = joinpath(dir, "stats.md")
    open(filepath, "w") do io
        println(io, "# DetectFit Statistics\n")
        println(io, "## Summary")
        println(io, "- **Datasets**: $n_datasets")
        println(io, "- **Frames/dataset**: $n_frames_per_dataset")
        println(io, "- **ROIs detected**: $n_rois")
        println(io, "- **Fits**: $n_fits")
        println(io, "- **Time**: $(round(t, digits=2))s ($(round(n_fits/t/1000, digits=1))k fits/s)")

        # Photobleaching rate (observed)
        if bleach_result !== nothing && bleach_result.r_squared > 0.5
            println(io, "")
            println(io, "## Photobleaching (from loc/frame decay)")
            println(io, "- **Model**: N(t) = a + b*exp(-k*t)")
            println(io, "- **k_observed**: $(round(bleach_result.k_bleach, sigdigits=3)) /frame")
            println(io, "- **Half-life**: $(round(bleach_result.half_life, digits=0)) frames")
            println(io, "- **a (offset)**: $(round(Int, bleach_result.offset)) loc/frame")
            println(io, "- **b (amplitude)**: $(round(Int, bleach_result.N_0)) loc/frame")
            println(io, "- **R^2**: $(round(bleach_result.r_squared, digits=3))")
            println(io, "")
            println(io, "*Note: k_observed = k_bleach * P_on. For GenericFluor, divide by duty cycle.*")
        end

        println(io, "")
        println(io, "## Detection Parameters")
        println(io, "- boxsize: $(cfg.boxer.boxsize)")
        println(io, "- min_photons: $(cfg.boxer.min_photons)")
        println(io, "- psf_sigma: $(cfg.boxer.psf_sigma)")
        println(io, "")
        println(io, "## Fit Parameters")
        println(io, "- psf_model: $(typeof(cfg.fitter.psf_model))")
        println(io, "- iterations: $(cfg.fitter.iterations)")
        println(io, "")
        println(io, "## Distributions\n")
        println(io, "| Parameter | Median | 5% | 95% |")
        println(io, "|-----------|--------|-----|-----|")
        println(io, "| Photons | $(round(median(photons), digits=0)) | $(round(quantile(photons, 0.05), digits=0)) | $(round(quantile(photons, 0.95), digits=0)) |")
        println(io, "| Background | $(round(median(bg), digits=1)) | $(round(quantile(bg, 0.05), digits=1)) | $(round(quantile(bg, 0.95), digits=1)) |")
        println(io, "| σ_x (nm) | $(round(median(σ_x)*1000, digits=1)) | $(round(quantile(σ_x, 0.05)*1000, digits=1)) | $(round(quantile(σ_x, 0.95)*1000, digits=1)) |")
        println(io, "| σ_y (nm) | $(round(median(σ_y)*1000, digits=1)) | $(round(quantile(σ_y, 0.05)*1000, digits=1)) | $(round(quantile(σ_y, 0.95)*1000, digits=1)) |")
        println(io, "")
        println(io, "## P-value")
        pval_pass = sum(pvalue .> 0.001) / n
        println(io, "- pvalue > 0.001: $(round(100*pval_pass, digits=1))%")
        println(io, "- pvalue > 0.01: $(round(100*sum(pvalue .> 0.01)/n, digits=1))%")
    end
end

"""
Save sample frame data to pipeline cache for the filter step to generate fit_overlay.png.

Decomposes ROIBatch into plain arrays (avoids camera serialization issues with JLD2).
"""
function _save_detectfit_sample_cache(outdir, smld, sample_images, sample_roi_batch, sample_original_frames)
    save_cache(outdir, "detectfit_samples.jld2";
        sample_images = sample_images,
        sample_roi_data = sample_roi_batch.data,
        sample_roi_x = sample_roi_batch.x_corners,
        sample_roi_y = sample_roi_batch.y_corners,
        sample_roi_frames = sample_roi_batch.frame_indices,
        sample_roi_size = sample_roi_batch.roi_size,
        sample_original_frames = sample_original_frames,
        n_frames = smld.n_frames,
        n_datasets = smld.n_datasets,
    )
end
