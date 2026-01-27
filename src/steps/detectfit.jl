"""
Combined Detect+Fit step with per-dataset processing.

Supports two modes:
1. **In-memory images**: When Analysis has loaded images (via constructor), processes
   datasets by slicing the image array.
2. **File-based loading**: When paths are specified, loads one dataset at a time
   for memory efficiency with arbitrarily large acquisitions.

For each dataset:
1. Get images (slice from memory or load from file)
2. Detect ROIs
3. Fit localizations
4. Append to combined SMLD
5. (File mode only) Free images from memory
"""

@kwdef struct DetectFitConfig <: StepConfig
    name::String = "detectfit"

    # Data source (optional - use Analysis.data if not specified)
    # Option 1: Single file with n_datasets (splits frames evenly)
    path::Union{String, Nothing} = nothing
    n_datasets::Int = 1  # Number of datasets (for file-based: splits frames evenly)
    # Option 2: Single file with explicit frame ranges
    dataset_frames::Union{Vector{UnitRange{Int}}, Nothing} = nothing
    # Option 3: Multiple files (one per dataset)
    paths::Union{Vector{String}, Nothing} = nothing

    # H5 format: :auto (detect), :smart (SMART microscope), :mic (MATLAB Instrument Control)
    h5_format::Symbol = :auto

    # Detection params (passed to SMLMBoxer.getboxes)
    boxsize::Int = 11
    overlap::Float64 = 2.0
    min_photons::Float64 = 500.0
    psf_sigma::Float64 = 0.135
    use_gpu::Bool = true

    # Fit params (passed to GaussMLE.fit)
    psf_model::Symbol = :variable  # :fixed, :variable, :anisotropic
    psf_sigma_fit::Float32 = 0.135f0  # For :fixed only
    iterations::Int = 20

    # Filter preview thresholds (for fit_quality plot visualization)
    # Set these to match your intended FilterConfig settings
    filter_min_photons::Float64 = 500.0
    filter_max_precision::Float64 = 0.007  # 7nm default
    filter_min_pvalue::Float64 = 1e-6

    # Extra
    verbose::Int = Verbosity.STANDARD
end

function run_step!(a::Analysis, cfg::DetectFitConfig)
    a.step_counter += 1
    v = _get_verbose(a, cfg)
    dir = _stepdir(a, cfg)

    # Determine data sources: use Analysis.data if available, else use file paths from config
    sources, source_mode = _resolve_data_sources(a, cfg)
    n_datasets = length(sources)

    v >= Verbosity.PROGRESS && @info "[$(a.step_counter)] $(cfg.name)" n_datasets=n_datasets psf_model=cfg.psf_model

    # Setup fitter
    psf = if cfg.psf_model == :fixed
        GaussianXYNB(cfg.psf_sigma_fit)
    elseif cfg.psf_model == :variable
        GaussianXYNBS()
    elseif cfg.psf_model == :anisotropic
        GaussianXYNBSXSY()
    else
        error("Unknown psf_model: $(cfg.psf_model)")
    end
    fitter = GaussMLEFitter(psf_model=psf, iterations=cfg.iterations)

    # Process each dataset
    all_emitters = AbstractEmitter[]
    n_frames_per_dataset = 0
    total_rois = 0
    total_fits = 0

    # Sample data for overlay plots (capture 12 frames from first dataset for 3x4 grid)
    sample_images = nothing
    sample_roi_batch = nothing
    sample_original_frames = nothing  # Original frame numbers for display
    n_sample_frames = 12

    t = @elapsed begin
        for (ds, source) in enumerate(sources)
            # Get images for this dataset (slice from memory or load from file)
            images = _get_dataset_images(source, source_mode, v)
            n_frames_ds = size(images, 3)
            n_frames_per_dataset = max(n_frames_per_dataset, n_frames_ds)

            v >= Verbosity.DETAILED && @info "  Dataset $ds: $(size(images)) images"

            # Detect
            roi_batch = SMLMBoxer.getboxes(images, a.camera;
                boxsize = cfg.boxsize,
                overlap = cfg.overlap,
                min_photons = cfg.min_photons,
                psf_sigma = cfg.psf_sigma,
                use_gpu = cfg.use_gpu
            )
            n_rois = length(roi_batch)
            total_rois += n_rois

            v >= Verbosity.DETAILED && @info "    Detected $n_rois ROIs"

            # Fit
            smld_ds = GaussMLE.fit(fitter, roi_batch)
            n_fits = length(smld_ds.emitters)
            total_fits += n_fits

            v >= Verbosity.DETAILED && @info "    Fitted $n_fits localizations"

            # Capture sample data from first dataset for overlay plots
            # Sample frames spread across the whole dataset, not just the beginning
            if ds == 1 && dir !== nothing
                n_sample = min(n_sample_frames, n_frames_ds)
                # Pick frames evenly spread across the dataset
                sample_frame_indices = [round(Int, x) for x in range(1, n_frames_ds, length=n_sample)]
                sample_images = collect(images[:, :, sample_frame_indices])
                sample_original_frames = sample_frame_indices  # Store for display titles

                # Create mapping from original frame index to sample index
                frame_to_sample = Dict(f => i for (i, f) in enumerate(sample_frame_indices))

                # Filter ROIs to only those in sample frames and remap frame indices
                sample_mask = [f in keys(frame_to_sample) for f in roi_batch.frame_indices]
                if any(sample_mask)
                    # Remap frame indices to 1:n_sample for the sample_images array
                    remapped_frames = [frame_to_sample[f] for f in roi_batch.frame_indices[sample_mask]]
                    sample_roi_batch = ROIBatch(
                        roi_batch.data[:, :, sample_mask],
                        roi_batch.x_corners[sample_mask],
                        roi_batch.y_corners[sample_mask],
                        remapped_frames,
                        roi_batch.camera
                    )
                end
            end

            # Set dataset field and append
            for e in smld_ds.emitters
                push!(all_emitters, _with_dataset(e, ds))
            end

            # Images freed when loop iteration ends (GC)
        end
    end

    # Create combined SMLD
    if isempty(all_emitters)
        error("No localizations found across all datasets")
    end

    a.smld_raw = BasicSMLD(all_emitters, a.camera, n_frames_per_dataset, n_datasets, Dict{String,Any}())
    a.smld = a.smld_raw
    a.n_datasets = n_datasets
    a.n_frames_per_dataset = n_frames_per_dataset

    summary = Dict{Symbol,Any}(
        :n_datasets => n_datasets,
        :n_rois => total_rois,
        :n_fits => total_fits,
        :n_frames_per_dataset => n_frames_per_dataset
    )
    _record!(a, cfg, t, summary)
    _checkpoint!(a)

    if dir !== nothing
        _save_detectfit_outputs!(dir, a, cfg, v, t, total_rois, total_fits, sample_images, sample_roi_batch, sample_original_frames)
    end

    v >= Verbosity.PROGRESS && @info "  → $total_fits fits from $total_rois ROIs across $n_datasets datasets ($(round(t, digits=2))s)"
    a
end

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
Resolve data sources from Analysis and/or config.

Returns (sources, mode) where:
- mode = :memory → sources are frame ranges into Analysis.data images
- mode = :file → sources are file paths/ranges for loading

Priority:
1. If config specifies paths/path, use file-based loading
2. If Analysis.data has images, use in-memory slicing
3. Error if no data source available
"""
function _resolve_data_sources(a::Analysis, cfg::DetectFitConfig)
    # Check if config specifies file-based loading
    if cfg.paths !== nothing || cfg.path !== nothing
        return _resolve_file_sources(cfg), :file
    end

    # Check if Analysis has in-memory images
    if a.data.images !== nothing
        return _resolve_memory_sources(a), :memory
    end

    # Check if Analysis has a path in DataSource
    if a.data.path !== nothing
        # Create a config-like source from Analysis.data
        return _resolve_datasource_path(a), :file
    end

    error("No data source: specify path/paths in DetectFitConfig or provide images to Analysis constructor")
end

"""Resolve file-based sources from config options"""
function _resolve_file_sources(cfg::DetectFitConfig)
    # Option 3: Multiple files
    if cfg.paths !== nothing
        format = cfg.h5_format == :auto ? _detect_h5_format(cfg.paths[1]) : cfg.h5_format
        return [(path=p, frame_range=nothing, format=format) for p in cfg.paths]
    end

    # Must have single path for options 1 & 2
    cfg.path === nothing && error("Must specify path or paths")

    # Detect format
    format = cfg.h5_format == :auto ? _detect_h5_format(cfg.path) : cfg.h5_format

    # Option 2: Explicit frame ranges
    if cfg.dataset_frames !== nothing
        return [(path=cfg.path, frame_range=r, format=format) for r in cfg.dataset_frames]
    end

    # Option 1: Single file with n_datasets
    if cfg.n_datasets == 1
        return [(path=cfg.path, frame_range=nothing, format=format)]
    end

    # Multiple datasets from single file - determine frame ranges
    if format == :smart
        info = load_smart_h5_info(cfg.path)
        n_frames = info.nframes
    else  # :mic
        info = load_lidkelab_h5_info(cfg.path)
        n_frames = info.n_frames
        # For MIC format: use block-based loading when n_datasets matches n_blocks
        if cfg.n_datasets == info.n_blocks
            return [(path=cfg.path, block=ds, format=format) for ds in 1:cfg.n_datasets]
        end
    end

    # Split frames evenly across datasets
    frames_per_ds = div(n_frames, cfg.n_datasets)
    sources = []
    for ds in 1:cfg.n_datasets
        start_frame = (ds - 1) * frames_per_ds + 1
        end_frame = ds == cfg.n_datasets ? n_frames : ds * frames_per_ds
        push!(sources, (path=cfg.path, frame_range=start_frame:end_frame, format=format))
    end
    return sources
end

"""Resolve in-memory sources by slicing Analysis.data into datasets"""
function _resolve_memory_sources(a::Analysis)
    n_datasets = a.n_datasets
    n_frames_per_dataset = a.n_frames_per_dataset
    images = a.data.images

    sources = []
    for ds in 1:n_datasets
        frame_start = (ds - 1) * n_frames_per_dataset + 1
        frame_end = ds * n_frames_per_dataset
        push!(sources, (images=images, frame_range=frame_start:frame_end))
    end
    return sources
end

"""Resolve sources when Analysis.data has a path but config doesn't"""
function _resolve_datasource_path(a::Analysis)
    path = a.data.path
    format = _detect_h5_format(path)

    # Use frame_range from DataSource if specified
    frame_range = a.data.frame_range

    # For now, single dataset from the path
    return [(path=path, frame_range=frame_range, format=format)]
end

"""Get images for a dataset from either memory or file"""
function _get_dataset_images(source, mode::Symbol, v::Int)
    if mode == :memory
        # Slice from in-memory images
        return @view source.images[:, :, source.frame_range]
    else
        # Load from file
        return _load_source(source, v)
    end
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
            return load_lidkelab_h5_block(source.path, source.block)
        elseif source.frame_range === nothing
            # Load all frames
            images, _ = load_lidkelab_h5(source.path)
            return images
        else
            # Load specific frame range (less efficient, but needed for custom ranges)
            images, _ = load_lidkelab_h5(source.path; max_frames=last(source.frame_range))
            return images[:, :, source.frame_range]
        end
    else
        error("Unknown H5 format: $(source.format)")
    end
end

function _save_detectfit_outputs!(dir, a, cfg, v, t, n_rois, n_fits, sample_images, sample_roi_batch, sample_original_frames)
    mkpath(dir)
    _save_config!(dir, cfg)

    if v >= Verbosity.STANDARD
        _write_detectfit_stats(dir, a, cfg, t, n_rois, n_fits)
        _save_detectfit_figures(dir, a.smld, cfg)

        # Generate overlay plots if we have sample data
        if sample_images !== nothing && sample_roi_batch !== nothing
            _save_detectfit_overlays(dir, a.smld, sample_roi_batch, sample_images, cfg, sample_original_frames)
        end
    end

    if v >= Verbosity.DETAILED
        _save_detectfit_detailed(dir, a.smld)
    end
end

function _write_detectfit_stats(dir, a, cfg, t, n_rois, n_fits)
    emitters = a.smld.emitters
    n = length(emitters)

    photons = [e.photons for e in emitters]
    bg = [e.bg for e in emitters]
    σ_x = [e.σ_x for e in emitters]
    σ_y = [e.σ_y for e in emitters]
    pvalue = [e.pvalue for e in emitters]

    # Calculate photobleaching rate from localizations per frame
    n_frames = a.smld.n_frames
    frame_counts = zeros(Int, n_frames)
    for e in emitters
        if e.frame >= 1 && e.frame <= n_frames
            frame_counts[e.frame] += 1
        end
    end
    bleach_result = _estimate_bleaching_rate(frame_counts)

    filepath = joinpath(dir, "stats.md")
    open(filepath, "w") do io
        println(io, "# DetectFit Statistics\n")
        println(io, "## Summary")
        println(io, "- **Datasets**: $(a.n_datasets)")
        println(io, "- **Frames/dataset**: $(a.n_frames_per_dataset)")
        println(io, "- **ROIs detected**: $n_rois")
        println(io, "- **Fits**: $n_fits")
        println(io, "- **Time**: $(round(t, digits=2))s ($(round(n_fits/t/1000, digits=1))k fits/s)")

        # Photobleaching rate (observed)
        if bleach_result !== nothing && bleach_result.r_squared > 0.5
            println(io, "")
            println(io, "## Photobleaching (from loc/frame decay)")
            println(io, "- **k_observed**: $(round(bleach_result.k_bleach, sigdigits=3)) /frame")
            println(io, "- **Half-life**: $(round(bleach_result.half_life, digits=0)) frames")
            println(io, "- **R²**: $(round(bleach_result.r_squared, digits=3))")
            println(io, "")
            println(io, "*Note: k_observed = k_bleach × P_on. For GenericFluor, divide by duty cycle.*")
        end

        println(io, "")
        println(io, "## Detection Parameters")
        println(io, "- boxsize: $(cfg.boxsize)")
        println(io, "- min_photons: $(cfg.min_photons)")
        println(io, "- psf_sigma: $(cfg.psf_sigma)")
        println(io, "")
        println(io, "## Fit Parameters")
        println(io, "- psf_model: $(cfg.psf_model)")
        println(io, "- iterations: $(cfg.iterations)")
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

function _save_detectfit_figures(dir, smld, cfg)
    emitters = smld.emitters
    isempty(emitters) && return

    photons = [e.photons for e in emitters]
    bg = [e.bg for e in emitters]
    σ_x = [e.σ_x for e in emitters] .* 1000  # precision in nm
    σ_y = [e.σ_y for e in emitters] .* 1000
    pvalue = [e.pvalue for e in emitters]

    # Check PSF model type
    has_psf_iso = hasproperty(emitters[1], :σ)
    has_psf_aniso = hasproperty(emitters[1], :σx)

    # Colors for consistent styling
    REJECTED_COLOR = (:gray30, 0.5)
    MEAN_COLOR = :blue
    MEDIAN_COLOR = :red
    THRESHOLD_COLOR = :black

    fig = Figure(size=(1200, 900))

    # Row 1: Photons and Background
    p98 = quantile(photons, 0.98)
    ax1 = Axis(fig[1, 1], xlabel="Photons", ylabel="Count", title="Photon Distribution")
    photons_lo = cfg.filter_min_photons
    vspan!(ax1, 0, photons_lo, color=REJECTED_COLOR)
    hist!(ax1, photons[photons .<= p98], bins=50)
    vlines!(ax1, [mean(photons)], color=MEAN_COLOR, linestyle=:solid, linewidth=2)
    vlines!(ax1, [median(photons)], color=MEDIAN_COLOR, linestyle=:dash, linewidth=2)
    vlines!(ax1, [photons_lo], color=THRESHOLD_COLOR, linestyle=:dot, linewidth=2)
    xlims!(ax1, 0, p98)
    text!(ax1, 0.97, 0.95, text="mean: $(round(Int, mean(photons)))\nmedian: $(round(Int, median(photons)))",
          align=(:right, :top), space=:relative, fontsize=10)

    bg98 = quantile(bg, 0.98)
    ax2 = Axis(fig[1, 2], xlabel="Background", ylabel="Count", title="Background Distribution")
    hist!(ax2, bg[bg .<= bg98], bins=50)
    vlines!(ax2, [mean(bg)], color=MEAN_COLOR, linestyle=:solid, linewidth=2)
    vlines!(ax2, [median(bg)], color=MEDIAN_COLOR, linestyle=:dash, linewidth=2)
    xlims!(ax2, 0, bg98)
    text!(ax2, 0.97, 0.95, text="mean: $(round(mean(bg), digits=1))\nmedian: $(round(median(bg), digits=1))",
          align=(:right, :top), space=:relative, fontsize=10)

    # Row 2: Precision and P-value
    prec_hi = cfg.filter_max_precision * 1000  # convert μm to nm
    prec_data = vcat(σ_x, σ_y)
    prec98 = quantile(prec_data, 0.98)
    # Ensure rejected region is at least 30% of plot width for visibility
    prec_xlim = max(prec98, prec_hi * 1.5)

    ax3 = Axis(fig[2, 1], xlabel="Localization Precision (nm)", ylabel="Count", title="Precision Distribution")
    vspan!(ax3, prec_hi, prec_xlim, color=REJECTED_COLOR)
    hist!(ax3, σ_x[σ_x .<= prec98], bins=50, color=(:blue, 0.5), label="σ_x")
    hist!(ax3, σ_y[σ_y .<= prec98], bins=50, color=(:red, 0.5), label="σ_y")
    vlines!(ax3, [prec_hi], color=THRESHOLD_COLOR, linestyle=:dot, linewidth=2)
    xlims!(ax3, 0, prec_xlim)
    axislegend(ax3, position=:rt, framevisible=false, labelsize=9)
    text!(ax3, 0.97, 0.70, text="σ_x: $(round(median(σ_x), digits=1)) nm\nσ_y: $(round(median(σ_y), digits=1)) nm",
          align=(:right, :top), space=:relative, fontsize=10)

    ax4 = Axis(fig[2, 2], xlabel="log₁₀(p-value)", ylabel="Density", title="P-value Distribution")
    pval_nonzero = pvalue[pvalue .> 0]
    pval_thresh = cfg.filter_min_pvalue

    if !isempty(pval_nonzero)
        log_pval = log10.(pval_nonzero)
        pval_lo = quantile(log_pval, 0.02)
        log_pval_filtered = log_pval[log_pval .>= pval_lo]
        vspan!(ax4, pval_lo - 1, log10(pval_thresh), color=REJECTED_COLOR)
        hist!(ax4, log_pval_filtered, bins=50, normalization=:pdf, color=(:steelblue, 0.7))
        # Compute histogram max for y limits (based on data, not theory)
        nbins = 50
        bin_edges = range(pval_lo, 0, length=nbins+1)
        bin_width = step(bin_edges)
        counts = zeros(Int, nbins)
        for v in log_pval_filtered
            idx = clamp(floor(Int, (v - pval_lo) / bin_width) + 1, 1, nbins)
            counts[idx] += 1
        end
        max_density = maximum(counts) / (length(log_pval_filtered) * bin_width)
        # Theory curve (uniform p-values → exponential in log space)
        u_range = range(0, -pval_lo, length=100)
        theory_pdf = log(10) .* (10.0 .^ (-u_range))
        lines!(ax4, -u_range, theory_pdf, color=:red, linewidth=2, label="Uniform theory")
        vlines!(ax4, [mean(log_pval)], color=MEAN_COLOR, linestyle=:solid, linewidth=1.5)
        vlines!(ax4, [median(log_pval)], color=MEDIAN_COLOR, linestyle=:dash, linewidth=1.5)
        vlines!(ax4, [log10(pval_thresh)], color=THRESHOLD_COLOR, linestyle=:dot, linewidth=2)
        xlims!(ax4, pval_lo, 0)
        ylims!(ax4, 0, max_density * 1.1)
        pval_pass_pct = round(100 * sum(pvalue .> pval_thresh) / length(pvalue), digits=1)
        text!(ax4, 0.03, 0.95, text="pass: $(pval_pass_pct)%\nthreshold: $(pval_thresh)",
              align=(:left, :top), space=:relative, fontsize=10)
    end

    # Row 3: PSF Sigma
    if has_psf_aniso
        psf_σx = [e.σx for e in emitters] .* 1000
        psf_σy = [e.σy for e in emitters] .* 1000
        psf_data = vcat(psf_σx, psf_σy)
        psf98 = quantile(psf_data, 0.98)
        psf02 = quantile(psf_data, 0.02)
        mode_x = _calculate_mode([e.σx for e in emitters]) * 1000
        mode_y = _calculate_mode([e.σy for e in emitters]) * 1000
        mode_avg = (mode_x + mode_y) / 2
        tol = 0.10

        ax5 = Axis(fig[3, 1], xlabel="Fitted PSF σ (nm)", ylabel="Count", title="PSF Width Distribution")
        vspan!(ax5, psf02 * 0.9, mode_avg * (1 - tol), color=REJECTED_COLOR)
        vspan!(ax5, mode_avg * (1 + tol), psf98 * 1.1, color=REJECTED_COLOR)
        hist!(ax5, psf_σx[(psf_σx .>= psf02) .& (psf_σx .<= psf98)], bins=50, color=(:blue, 0.5), label="σx")
        hist!(ax5, psf_σy[(psf_σy .>= psf02) .& (psf_σy .<= psf98)], bins=50, color=(:red, 0.5), label="σy")
        vlines!(ax5, [mode_avg * (1 - tol), mode_avg * (1 + tol)], color=THRESHOLD_COLOR, linestyle=:dot, linewidth=2)
        xlims!(ax5, psf02, psf98)
        axislegend(ax5, position=:rt, framevisible=false, labelsize=9)
        text!(ax5, 0.97, 0.70, text="mode σx: $(round(mode_x, digits=1)) nm\nmode σy: $(round(mode_y, digits=1)) nm\ntol: ±$(round(Int, tol*100))%",
              align=(:right, :top), space=:relative, fontsize=10)

    elseif has_psf_iso
        psf_σ = [e.σ for e in emitters] .* 1000
        psf98 = quantile(psf_σ, 0.98)
        psf02 = quantile(psf_σ, 0.02)
        mode_σ = _calculate_mode([e.σ for e in emitters]) * 1000
        tol = 0.10
        lo_bound = mode_σ * (1 - tol)
        hi_bound = mode_σ * (1 + tol)
        psf_xmin = min(psf02, lo_bound * 0.95)
        psf_xmax = max(psf98, hi_bound * 1.05)

        ax5 = Axis(fig[3, 1], xlabel="Fitted PSF σ (nm)", ylabel="Count", title="PSF Width Distribution")
        vspan!(ax5, psf_xmin, lo_bound, color=REJECTED_COLOR)
        vspan!(ax5, hi_bound, psf_xmax, color=REJECTED_COLOR)
        hist!(ax5, psf_σ[(psf_σ .>= psf02) .& (psf_σ .<= psf98)], bins=50)
        vlines!(ax5, [mean(psf_σ)], color=MEAN_COLOR, linestyle=:solid, linewidth=2)
        vlines!(ax5, [median(psf_σ)], color=MEDIAN_COLOR, linestyle=:dash, linewidth=2)
        vlines!(ax5, [lo_bound, hi_bound], color=THRESHOLD_COLOR, linestyle=:dot, linewidth=2)
        xlims!(ax5, psf_xmin, psf_xmax)
        text!(ax5, 0.97, 0.95, text="mode: $(round(mode_σ, digits=1)) nm\nmean: $(round(mean(psf_σ), digits=1)) nm\ntol: ±$(round(Int, tol*100))%",
              align=(:right, :top), space=:relative, fontsize=10)
    else
        ax5 = Axis(fig[3, 1], xlabel="PSF σ (nm)", ylabel="", title="PSF Width (Fixed)")
        fixed_σ = cfg.psf_sigma_fit * 1000
        vlines!(ax5, [fixed_σ], color=:blue, linewidth=3)
        text!(ax5, 0.5, 0.5, text="Fixed: $(round(fixed_σ, digits=1)) nm",
              align=(:center, :center), space=:relative, fontsize=14)
        hideydecorations!(ax5)
    end

    # Row 3, Col 2: Legend
    ax6 = Axis(fig[3, 2], title="Legend")
    hidedecorations!(ax6)
    hidespines!(ax6)
    lines!(ax6, [0.1, 0.25], [0.8, 0.8], color=MEAN_COLOR, linewidth=2)
    text!(ax6, 0.3, 0.8, text="Mean", fontsize=12)
    lines!(ax6, [0.1, 0.25], [0.6, 0.6], color=MEDIAN_COLOR, linewidth=2, linestyle=:dash)
    text!(ax6, 0.3, 0.6, text="Median", fontsize=12)
    lines!(ax6, [0.1, 0.25], [0.4, 0.4], color=THRESHOLD_COLOR, linewidth=2, linestyle=:dot)
    text!(ax6, 0.3, 0.4, text="Filter Threshold", fontsize=12)
    poly!(ax6, Point2f[(0.1, 0.15), (0.25, 0.15), (0.25, 0.25), (0.1, 0.25)], color=REJECTED_COLOR)
    text!(ax6, 0.3, 0.2, text="Rejected Region", fontsize=12)
    xlims!(ax6, 0, 1)
    ylims!(ax6, 0, 1)

    save(joinpath(dir, "fit_quality.png"), fig)
end

"""
Generate overlay plots showing detection boxes colored by fit status.
Uses sample frames from first dataset to avoid loading all data again.
"""
function _save_detectfit_overlays(dir, smld, sample_roi_batch, sample_images, cfg, sample_original_frames)
    isempty(sample_roi_batch) && return

    # Create mapping from sample index (1:N) to original frame number
    sample_to_original = Dict(i => f for (i, f) in enumerate(sample_original_frames))
    original_frame_set = Set(sample_original_frames)

    # Get emitters only from the sampled frames (using original frame numbers)
    sample_emitters = filter(e -> e.frame in original_frame_set, smld.emitters)

    # Detection overlay: all boxes yellow (detection view)
    box_colors = fill(:yellow, length(sample_roi_batch))
    _save_box_overlay(dir, "detection_overlay.png", sample_images, sample_roi_batch, box_colors;
                      title_prefix="Detection Frame", frame_labels=sample_original_frames)

    # Fit overlay: boxes colored by fit status
    # Match emitters to ROIs by position (approximate)
    if !isempty(sample_emitters)
        fit_colors = Symbol[]
        # Get pixel size from camera edges (assumes uniform pixels)
        pix_size = smld.camera.pixel_edges_x[2] - smld.camera.pixel_edges_x[1]
        for i in 1:length(sample_roi_batch)
            roi_frame_idx = sample_roi_batch.frame_indices[i]  # Remapped index (1:N)
            original_frame = sample_to_original[roi_frame_idx]  # Original frame number
            roi_x = sample_roi_batch.x_corners[i] + sample_roi_batch.roi_size ÷ 2
            roi_y = sample_roi_batch.y_corners[i] + sample_roi_batch.roi_size ÷ 2

            # Find matching emitter using original frame number
            best_emitter = nothing
            best_dist = Inf
            for e in sample_emitters
                if e.frame == original_frame
                    # Convert emitter position (microns) to pixels
                    ex_px = e.x / pix_size
                    ey_px = e.y / pix_size
                    dist = sqrt((ex_px - roi_x)^2 + (ey_px - roi_y)^2)
                    if dist < best_dist
                        best_dist = dist
                        best_emitter = e
                    end
                end
            end

            if best_emitter !== nothing && best_dist < sample_roi_batch.roi_size
                push!(fit_colors, _detectfit_box_color(best_emitter;
                    min_photons=cfg.filter_min_photons,
                    max_precision=cfg.filter_max_precision,
                    min_pvalue=cfg.filter_min_pvalue))
            else
                push!(fit_colors, :gray)  # No matching emitter found
            end
        end
        _save_fit_overlay_with_legend(dir, sample_images, sample_roi_batch, fit_colors, cfg, sample_original_frames)
    end
end

"""Save fit overlay with color info in figure title"""
function _save_fit_overlay_with_legend(dir, images, roi_batch, box_colors, cfg, frame_labels)
    n_frames = size(images, 3)
    frame_indices = [round(Int, x) for x in range(1, n_frames, length=min(12, n_frames))]
    # Use provided frame_labels for display titles
    display_labels = frame_labels !== nothing ? frame_labels : frame_indices

    # Contrast stretch
    sample_frames = frame_indices[1:min(4, length(frame_indices))]
    sample_data = vec(images[:, :, sample_frames])
    bg_level = median(sample_data)
    pmin = Float64(bg_level)
    pmax = Float64(quantile(sample_data, 0.995))

    # Build title with color legend info
    prec_nm = round(cfg.filter_max_precision * 1000, digits=1)
    title_str = "Fit Status: green=pass, red=photons<$(round(Int, cfg.filter_min_photons)), orange=prec>$(prec_nm)nm, purple=pval<$(cfg.filter_min_pvalue), gray=no match"

    n_rows = ceil(Int, length(frame_indices) / 4)
    fig = Figure(size=(1200, 50 + 250 * n_rows))
    box_size = roi_batch.roi_size

    # Add title at top
    Label(fig[0, 1:4], title_str, fontsize=12)

    for (idx, frame_num) in enumerate(frame_indices)
        row = div(idx - 1, 4) + 1
        col = mod(idx - 1, 4) + 1
        display_frame = display_labels[idx]

        ax = Axis(fig[row, col], title="Frame $display_frame", aspect=DataAspect(), yreversed=true)
        frame_data = images[:, :, frame_num]'
        heatmap!(ax, frame_data, colormap=:grays, colorrange=(pmin, pmax))

        frame_mask = roi_batch.frame_indices .== frame_num
        if any(frame_mask)
            det_x = roi_batch.x_corners[frame_mask]
            det_y = roi_batch.y_corners[frame_mask]
            colors = box_colors[frame_mask]
            for (x, y, c) in zip(det_x, det_y, colors)
                lines!(ax, [x, x+box_size, x+box_size, x, x],
                          [y, y, y+box_size, y+box_size, y],
                    color=c, linewidth=0.5)
            end
        end
        hidedecorations!(ax)
    end

    save(joinpath(dir, "fit_overlay.png"), fig)
end

"""Determine box color based on fit status (like fit.jl)"""
function _detectfit_box_color(e; min_photons=500.0, max_precision=0.007, min_pvalue=1e-6)
    e.photons < min_photons && return :red           # failed photons
    prec = sqrt(e.σ_x^2 + e.σ_y^2) / sqrt(2)
    prec > max_precision && return :orange           # failed precision
    e.pvalue < min_pvalue && return :purple          # failed pvalue
    return :green                                     # accepted
end

"""
Estimate observed bleaching rate from localizations per frame decay.
Fits N(t) = N_0 * exp(-k_observed * t) using linear regression on log(N).

Note: This gives the OBSERVED decay rate, not the true k_bleach for GenericFluor.
Since bleaching only occurs from the On state:
    k_observed = k_bleach * P_on
    where P_on = k_on/(k_on+k_off) is the duty cycle

To get true k_bleach: k_bleach = k_observed / P_on

Returns (k_bleach, N_0, half_life, r_squared) or nothing if fit fails.
"""
function _estimate_bleaching_rate(frame_counts::Vector{Int})
    # Filter out zero counts and use frames with sufficient data
    valid_mask = frame_counts .> 0
    valid_frames = findall(valid_mask)
    valid_counts = frame_counts[valid_mask]

    length(valid_counts) < 10 && return nothing

    # Smooth the data (rolling average) to reduce noise
    window = min(50, length(valid_counts) ÷ 10)
    if window > 1
        smoothed = [mean(valid_counts[max(1, i-window):min(end, i+window)]) for i in 1:length(valid_counts)]
    else
        smoothed = Float64.(valid_counts)
    end

    # Linear regression on log(N) vs frame
    # log(N) = log(N_0) - k_bleach * t
    x = Float64.(valid_frames)
    y = log.(max.(smoothed, 1.0))  # Avoid log(0)

    n = length(x)
    sum_x = sum(x)
    sum_y = sum(y)
    sum_xy = sum(x .* y)
    sum_x2 = sum(x .^ 2)

    denom = n * sum_x2 - sum_x^2
    abs(denom) < 1e-10 && return nothing

    slope = (n * sum_xy - sum_x * sum_y) / denom
    intercept = (sum_y - slope * sum_x) / n

    # k_bleach should be positive (slope is negative for decay)
    k_bleach = -slope
    k_bleach <= 0 && return nothing  # No decay detected

    N_0 = exp(intercept)
    half_life = log(2) / k_bleach

    # R² calculation
    y_pred = intercept .+ slope .* x
    ss_res = sum((y .- y_pred) .^ 2)
    ss_tot = sum((y .- mean(y)) .^ 2)
    r_squared = ss_tot > 0 ? 1 - ss_res / ss_tot : 0.0

    (k_bleach=k_bleach, N_0=N_0, half_life=half_life, r_squared=r_squared,
     valid_frames=valid_frames, smoothed=smoothed)
end

"""Generate detailed plots for DETAILED verbosity"""
function _save_detectfit_detailed(dir, smld)
    emitters = smld.emitters
    isempty(emitters) && return nothing

    # ROIs per frame plot
    n_frames = smld.n_frames
    frame_counts = zeros(Int, n_frames)
    for e in emitters
        if e.frame >= 1 && e.frame <= n_frames
            frame_counts[e.frame] += 1
        end
    end

    # Estimate photobleaching rate
    bleach_result = _estimate_bleaching_rate(frame_counts)

    fig = Figure(size=(900, 400))
    ax = Axis(fig[1, 1], xlabel="Frame", ylabel="Localizations", title="Localizations per Frame")
    lines!(ax, 1:n_frames, frame_counts, color=(:blue, 0.5), linewidth=0.5, label="Raw")

    # Add exponential decay fit if successful
    if bleach_result !== nothing
        # Plot smoothed data
        lines!(ax, bleach_result.valid_frames, bleach_result.smoothed,
               color=:blue, linewidth=1.5, label="Smoothed")

        # Plot fitted exponential
        fit_frames = 1:n_frames
        fit_counts = bleach_result.N_0 .* exp.(-bleach_result.k_bleach .* fit_frames)
        lines!(ax, fit_frames, fit_counts, color=:red, linewidth=2, linestyle=:dash,
               label="Fit: k=$(round(bleach_result.k_bleach, sigdigits=3))/frame")

        # Add annotation
        text!(ax, 0.95, 0.95,
            text="k_observed = $(round(bleach_result.k_bleach, sigdigits=3)) /frame\nτ_1/2 = $(round(bleach_result.half_life, digits=0)) frames\nR² = $(round(bleach_result.r_squared, digits=3))",
            align=(:right, :top),
            space=:relative,
            fontsize=11)
    else
        hlines!(ax, [mean(frame_counts)], color=:red, linestyle=:dash,
                label="mean ($(round(mean(frame_counts), digits=1)))")
    end

    axislegend(ax, position=:rt, framevisible=false)
    save(joinpath(dir, "localizations_per_frame.png"), fig)

    return bleach_result
end
