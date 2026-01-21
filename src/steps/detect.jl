"""
Detection step - wraps SMLMBoxer.getboxes
"""

@kwdef struct DetectConfig <: StepConfig
    name::String = "detect"
    # SMLMBoxer.getboxes kwargs
    boxsize::Int = 11
    overlap::Float64 = 2.0
    min_photons::Float64 = 500.0
    psf_sigma::Float64 = 0.135
    use_gpu::Bool = true
    # Extra
    verbose::Int = Verbosity.STANDARD
end

function run_step!(a::Analysis, cfg::DetectConfig)
    a.step_counter += 1
    v = _get_verbose(a, cfg)
    dir = _stepdir(a, cfg)

    v >= Verbosity.PROGRESS && @info "[$(a.step_counter)] $(cfg.name)" boxsize=cfg.boxsize min_photons=cfg.min_photons

    images = get_images(a.data)

    # Loop over datasets for memory efficiency and natural per-dataset frame numbers
    t = @elapsed begin
        all_roi_batches = ROIBatch[]
        all_datasets = Int[]  # Track dataset for each ROI

        for ds in 1:a.n_datasets
            # Extract this dataset's frames
            frame_start = (ds - 1) * a.n_frames_per_dataset + 1
            frame_end = ds * a.n_frames_per_dataset
            dataset_images = @view images[:, :, frame_start:frame_end]

            # Detect ROIs for this dataset (frames will be 1:n_frames_per_dataset)
            roi_batch_ds = SMLMBoxer.getboxes(dataset_images, a.camera;
                boxsize = cfg.boxsize,
                overlap = cfg.overlap,
                min_photons = cfg.min_photons,
                psf_sigma = cfg.psf_sigma,
                use_gpu = cfg.use_gpu
            )

            push!(all_roi_batches, roi_batch_ds)
            append!(all_datasets, fill(ds, length(roi_batch_ds)))

            v >= Verbosity.DETAILED && @info "  Dataset $ds: $(length(roi_batch_ds)) ROIs"
        end

        # Concatenate ROI batches and store dataset info
        a.roi_batch = _concat_roi_batches(all_roi_batches)
        a.roi_datasets = all_datasets  # Store dataset indices separately
    end

    n_rois = length(a.roi_batch)
    if n_rois == 0
        error("No particles detected. Try lowering min_photons (currently $(cfg.min_photons))")
    end

    summary = Dict{Symbol,Any}(:n_rois => n_rois, :n_frames => size(images, 3), :n_datasets => a.n_datasets)
    _record!(a, cfg, t, summary)
    _checkpoint!(a)  # Auto-checkpoint after expensive step

    if dir !== nothing
        _save_step_outputs!(dir, a, cfg, v, t, images)
    end

    v >= Verbosity.PROGRESS && @info "  → $n_rois ROIs ($(round(t, digits=2))s)"
    a
end

"""Concatenate multiple ROI batches into one"""
function _concat_roi_batches(batches::Vector{ROIBatch})
    isempty(batches) && error("No ROI batches to concatenate")
    length(batches) == 1 && return batches[1]

    # Concatenate arrays
    # ROIBatch.data is 3D: (roi_size, roi_size, n_rois) - concatenate along 3rd dim
    all_data = cat([b.data for b in batches]..., dims=3)
    all_x_corners = vcat([b.x_corners for b in batches]...)
    all_y_corners = vcat([b.y_corners for b in batches]...)
    all_frame_indices = vcat([b.frame_indices for b in batches]...)

    # ROIBatch constructor computes roi_size from data dimensions
    ROIBatch(all_data, all_x_corners, all_y_corners, all_frame_indices,
             batches[1].camera)
end

function _save_step_outputs!(dir::String, a::Analysis, cfg::DetectConfig, v::Int, t::Float64, images)
    mkpath(dir)
    _save_config!(dir, cfg)

    if v >= Verbosity.STANDARD
        _write_detect_stats(dir, a.roi_batch, images, cfg, t)
        _save_detect_figures(dir, a.roi_batch, images, a.camera, cfg)
    end

    if v >= Verbosity.DETAILED
        _save_detect_detailed(dir, a.roi_batch, images, cfg)
    end

    if v >= Verbosity.DEBUG
        _save_detect_debug(dir, a.roi_batch, images, cfg)
    end
end

function _write_detect_stats(dir, roi_batch, images, cfg, t)
    n_rois = length(roi_batch)
    n_frames = size(images, 3)
    rois_per_frame = n_rois / n_frames

    filepath = joinpath(dir, "stats.md")
    open(filepath, "w") do io
        println(io, "# Detection Statistics\n")
        println(io, "## Summary")
        println(io, "- **ROIs detected**: $n_rois")
        println(io, "- **Frames**: $n_frames")
        println(io, "- **ROIs/frame**: $(round(rois_per_frame, digits=1))")
        println(io, "- **Time**: $(round(t, digits=2))s")
        println(io, "")
        println(io, "## Parameters")
        println(io, "- boxsize: $(cfg.boxsize)")
        println(io, "- min_photons: $(cfg.min_photons)")
        println(io, "- psf_sigma: $(cfg.psf_sigma)")
        println(io, "- overlap: $(cfg.overlap)")
    end
end

function _save_detect_figures(dir, roi_batch, images, camera, cfg)
    # All detection boxes are yellow (no fit status yet)
    box_colors = fill(:yellow, length(roi_batch))
    _save_box_overlay(dir, "detection_overlay.png", images, roi_batch, box_colors)
end

function _save_detect_detailed(dir, roi_batch, images, cfg)
    # Per-frame ROI counts
    n_frames = size(images, 3)
    counts = zeros(Int, n_frames)
    for f in roi_batch.frame_indices
        counts[f] += 1
    end

    fig = Figure(size=(800, 400))
    ax = Axis(fig[1, 1], xlabel="Frame", ylabel="ROIs", title="ROIs per Frame")
    lines!(ax, 1:n_frames, counts)
    hlines!(ax, [mean(counts)], color=:red, linestyle=:dash, label="mean")
    save(joinpath(dir, "rois_per_frame.png"), fig)
end

function _save_detect_debug(dir, roi_batch, images, cfg)
    # TODO: MP4 of detections over time
    # For now, save ROI montage
    n_rois = min(25, length(roi_batch))
    if n_rois == 0
        return
    end

    roi_size = roi_batch.roi_size
    ncols = 5
    nrows = ceil(Int, n_rois / ncols)

    fig = Figure(size=(ncols * 100, nrows * 100))
    for i in 1:n_rois
        row = div(i - 1, ncols) + 1
        col = mod(i - 1, ncols) + 1
        ax = Axis(fig[row, col], aspect=DataAspect())
        heatmap!(ax, roi_batch.data[:, :, i], colormap=:grays)
        hidedecorations!(ax)
        hidespines!(ax)
    end
    save(joinpath(dir, "roi_montage.png"), fig)
end

