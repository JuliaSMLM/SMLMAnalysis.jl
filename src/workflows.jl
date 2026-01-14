"""
    workflows.jl

High-level workflow functions that chain multiple SMLM analysis steps together.
Each workflow function handles data transformations and optionally tracks
processing steps in an SMLMWorkflow object.
"""

using SMLMData
using SMLMSim
using SMLMBoxer
using GaussMLE
using MicroscopePSFs
using SMLMRender
using SMLMDriftCorrection
using SMLMFrameConnection
using NearestNeighbors
using CairoMakie
using Statistics
using StatsBase: countmap
using TOML
using LinearAlgebra: det, inv

# =============================================================================
# AnalysisConfig - Configuration for SMLM analysis workflow
# =============================================================================

"""
    AnalysisConfig

Configuration struct for SMLM analysis workflow. Fields are ordered by pipeline flow.
Use `true` to enable a step, `false` to skip it.

# Example
```julia
# All defaults
config = AnalysisConfig()

# Custom configuration
config = AnalysisConfig(
    detect_min_photons = 500.0,  # Detection threshold
    min_photons = 500.0,         # Filter threshold
    render = true,
    outdir = "output/"
)

# Save/load config
save_config("myconfig.toml", config)
config = load_config("myconfig.toml")
```
"""
Base.@kwdef struct AnalysisConfig
    # === Detection (always runs) ===
    boxsize::Int = 11
    overlap::Float64 = 2.0
    detect_min_photons::Float64 = 500.0  # Detection threshold in photons (auto-converted to DoG threshold)
    detect_gpu::Bool = true

    # === PSF sigma (used for both detection and fixed-sigma fitting) ===
    # Detection: DoG uses sigma_small = psf_sigma, sigma_large = 2*psf_sigma (in pixels)
    # Fitting: Only used when fit_model=:fixed (GaussianXYNB)
    psf_sigma::Float32 = 0.135f0    # PSF sigma in MICRONS (~135nm typical for TIRF)

    # === Fitting (always runs) ===
    # PSF model: :fixed (GaussianXYNB), :variable (GaussianXYNBS), :anisotropic (GaussianXYNBSXSY)
    fit_model::Symbol = :variable   # Default to variable sigma for better pvalue
    iterations::Int = 20
    fit_device::Symbol = :auto      # :cpu, :gpu, :auto

    # === Filtering (runs after fitting) ===
    filter::Bool = true
    min_photons::Union{Float64, Nothing} = 500.0
    max_precision::Union{Float64, Nothing} = 0.015  # Max localization precision in μm (15nm default)
    psf_sigma_mode_tolerance::Union{Float64, Nothing} = 0.10  # Keep PSF sigma within ±10% of mode (for variable sigma fits)
    min_pvalue::Union{Float64, Nothing} = 1e-3    # p-value threshold

    # === Frame Connection (runs after filtering) ===
    frameconnect::Bool = false      # OFF by default
    fc_maxframegap::Int = 5         # Max frame gap between connected locs
    fc_nsigmadev::Float64 = 5.0     # Sigma multiplier for preclustering distance
    fc_nnearestclusters::Int = 2    # Nearest clusters for density estimation
    fc_nmaxnn::Int = 2              # Max nearest neighbors for preclustering

    # === Uncertainty Calibration (runs after frame connection) ===
    calibrate_uncertainties::Bool = true  # Adjust uncertainties using k and σ_motion from frame connection

    # === Drift Correction (runs after uncertainty calibration) ===
    drift::Bool = true              # ON by default for DNA-PAINT
    drift_degree::Int = 2           # Polynomial degree (2 usually sufficient)
    drift_cost_fun::String = "Kdtree"  # "Kdtree" (fast) or "Entropy"
    drift_model::String = "Polynomial" # "Polynomial" or "LegendrePoly"
    dataset_indices::Union{Nothing, Vector{Int}} = nothing  # Map frame → dataset ID for multi-dataset drift correction

    # === Isolated Emitter Filter (runs after drift correction) ===
    filter_isolated::Bool = false   # OFF by default
    isolated_n_sigma::Float64 = 2.0 # Neighbor if dist < n_sigma * sqrt(σ_i² + σ_j²)
    isolated_min_neighbors::Union{Int, Symbol} = :auto  # :auto uses triangle method, or set Int manually

    # === Rendering ===
    render::Bool = true             # ON by default - primary output
    render_gaussian::Bool = true    # Gaussian blur @ 20x + inferno
    render_histogram::Bool = true   # Histogram @ 10x + time coloring
    render_circles::Bool = true     # Circles @ 50x + time coloring
    render_gaussian_zoom::Int = 20
    render_histogram_zoom::Int = 10
    render_circles_zoom::Int = 50
    render_time_colormap::Symbol = :turbo  # colormap for time-colored renders
    render_clip_percentile::Union{Float64, Symbol} = :auto  # :auto adapts to n_locs, or set 0.0-1.0

    # === BAGOL (Deep Learning) ===
    bagol::Bool = false             # OFF by default
    bagol_model::String = "default"

    # === Output ===
    outdir::Union{String, Nothing} = nothing
    save_figures::Bool = true
    save_smld::Bool = true
end

"""
    save_config(filepath, config::AnalysisConfig)

Save an AnalysisConfig to a TOML file.
"""
function save_config(filepath::String, config::AnalysisConfig)
    d = Dict{String, Any}()
    for field in fieldnames(AnalysisConfig)
        val = getfield(config, field)
        if val !== nothing
            # Convert Symbol to String for TOML compatibility
            d[string(field)] = val isa Symbol ? string(val) : val
        end
    end
    open(filepath, "w") do io
        TOML.print(io, d)
    end
    return filepath
end

"""
    load_config(filepath) -> AnalysisConfig

Load an AnalysisConfig from a TOML file.
"""
function load_config(filepath::String)
    d = TOML.parsefile(filepath)
    kwargs = Dict{Symbol, Any}()
    # Fields that should be Symbols
    symbol_fields = (:fit_device, :render_colormap, :render_strategy)
    for (k, v) in d
        sym = Symbol(k)
        # Convert string back to Symbol for appropriate fields
        if sym in symbol_fields && v isa String
            kwargs[sym] = Symbol(v)
        else
            kwargs[sym] = v
        end
    end
    return AnalysisConfig(; kwargs...)
end

# =============================================================================
# AnalysisResult - Return type from analyze()
# =============================================================================

"""
    AnalysisResult

Result container from `analyze()`. Contains the final SMLD plus intermediate results.

# Fields
- `smld`: Final BasicSMLD (after all processing)
- `smld_raw`: Raw fitted SMLD (before filtering/frameconnect)
- `roi_batch`: Detection results
- `timings`: Dict of step timings in seconds
- `workflow`: SMLMWorkflow provenance tracking
"""
struct AnalysisResult{T}
    smld::BasicSMLD{T}
    smld_raw::BasicSMLD{T}
    roi_batch::ROIBatch
    timings::Dict{String, Float64}
    workflow::SMLMWorkflow
end

# Convenience accessors
Base.getproperty(r::AnalysisResult, s::Symbol) =
    s == :emitters ? r.smld.emitters : getfield(r, s)

function Base.show(io::IO, r::AnalysisResult)
    n_raw = length(r.smld_raw.emitters)
    n_final = length(r.smld.emitters)
    total_time = sum(values(r.timings))
    print(io, "AnalysisResult: $n_final localizations ($n_raw raw) in $(round(total_time, digits=2))s")
end

# =============================================================================
# analyze() - Main workflow function
# =============================================================================

"""
    analyze(data, camera) -> AnalysisResult
    analyze(data, camera; kwargs...) -> AnalysisResult
    analyze(data, camera, config::AnalysisConfig) -> AnalysisResult

Run SMLM analysis workflow on image data.

# Arguments
- `data::Array{<:Real, 3}`: Image stack (height, width, nframes)
- `camera::AbstractCamera`: Camera calibration

# Keyword Arguments
Any field from `AnalysisConfig` can be passed as a keyword argument.
See `?AnalysisConfig` for full list.

# Returns
`AnalysisResult` containing:
- `smld`: Final processed localizations
- `smld_raw`: Raw fitted localizations (before filtering)
- `roi_batch`: Detection results
- `timings`: Dict of step timings
- `workflow`: Provenance tracking

# Examples
```julia
# All defaults
result = analyze(data, camera)

# With custom parameters (detection and filter thresholds)
result = analyze(data, camera; detect_min_photons=500.0, min_photons=500.0)

# Skip frame connection, enable rendering
result = analyze(data, camera; frameconnect=false, render=true, outdir="output/")

# From config file
config = load_config("myconfig.toml")
result = analyze(data, camera, config)
```
"""
# Kwargs version - creates config from kwargs (handles zero-arg case too)
function analyze(data::AbstractArray{<:Real, 3}, camera::SMLMData.AbstractCamera; kwargs...)
    return analyze(data, camera, AnalysisConfig(; kwargs...))
end

function analyze(data::AbstractArray{<:Real, 3}, camera::SMLMData.AbstractCamera, config::AnalysisConfig)
    timings = Dict{String, Float64}()
    workflow = SMLMWorkflow("SMLM Analysis")

    # Setup output directory if specified
    # Directory numbering matches pipeline order:
    # detect → fit → filter → frameconnect → calibrate → drift → isolated → render
    if config.outdir !== nothing
        mkpath(config.outdir)
        mkpath(joinpath(config.outdir, "01_detection"))
        mkpath(joinpath(config.outdir, "02_fitting"))
        mkpath(joinpath(config.outdir, "03_filtered"))
        if config.frameconnect
            mkpath(joinpath(config.outdir, "04_frameconnect"))
            if config.calibrate_uncertainties
                mkpath(joinpath(config.outdir, "05_calibration"))
            end
        end
        if config.drift
            mkpath(joinpath(config.outdir, "06_drift"))
        end
        if config.filter_isolated
            mkpath(joinpath(config.outdir, "07_isolated"))
        end
        if config.render
            mkpath(joinpath(config.outdir, "08_superres"))
        end
        mkpath(joinpath(config.outdir, "results"))

        # Save config for reproducibility
        save_config(joinpath(config.outdir, "config.toml"), config)
    end

    println("="^60)
    println("SMLM Analysis")
    println("="^60)
    println("  Data: $(size(data, 1))×$(size(data, 2)) × $(size(data, 3)) frames")

    # =========================================================================
    # Step 1: Detection (PSF-aware DoG)
    # =========================================================================
    print("  Detection... ")
    t = @elapsed roi_batch = getboxes(data, camera;
        boxsize = config.boxsize,
        overlap = config.overlap,
        psf_sigma = config.psf_sigma,  # DoG derives sigma_small=psf_sigma, sigma_large=2*psf_sigma
        min_photons = config.detect_min_photons,  # SMLMBoxer converts to DoG threshold automatically
        use_gpu = config.detect_gpu
    )
    timings["detection"] = t
    println("$(length(roi_batch)) ROIs ($(round(t, digits=2))s)")

    add_step!(workflow, "Detection", :getboxes,
        Dict{Symbol,Any}(:boxsize => config.boxsize, :psf_sigma => config.psf_sigma, :min_photons => config.detect_min_photons),
        :SMLMBoxer, "$(length(roi_batch)) ROIs")

    if length(roi_batch) == 0
        error("No particles detected! Check detect_min_photons threshold (currently $(config.detect_min_photons))")
    end

    # Save detection figures and stats
    if config.outdir !== nothing && config.save_figures
        _save_detection_figures(data, roi_batch, camera, config)
        _write_detection_stats(roi_batch, data, config, t)
    end

    # =========================================================================
    # Step 2: Fitting
    # =========================================================================
    print("  Fitting... ")
    device = config.fit_device == :auto ? nothing : config.fit_device

    # Select PSF model based on fit_model option
    psf_model = if config.fit_model == :fixed
        GaussianXYNB(config.psf_sigma)
    elseif config.fit_model == :variable
        GaussianXYNBS()  # Fits isotropic sigma
    elseif config.fit_model == :anisotropic
        GaussianXYNBSXSY()  # Fits sigma_x, sigma_y
    else
        error("Unknown fit_model: $(config.fit_model). Use :fixed, :variable, or :anisotropic")
    end

    fitter = GaussMLEFitter(
        psf_model = psf_model,
        iterations = config.iterations,
        device = device
    )
    t = @elapsed smld_raw = GaussMLE.fit(fitter, roi_batch)
    timings["fitting"] = t
    println("$(length(smld_raw.emitters)) fits ($(round(t, digits=2))s)")

    add_step!(workflow, "Fitting", :fit,
        Dict{Symbol,Any}(:psf_sigma => config.psf_sigma, :iterations => config.iterations),
        :GaussMLE, "$(length(smld_raw.emitters)) fits")

    smld = smld_raw  # Will be modified by subsequent steps

    # Save fitting figures and stats
    if config.outdir !== nothing && config.save_figures
        _save_fitting_figures(smld_raw, roi_batch, data, camera, config)
        _write_fitting_stats(smld_raw, config, t)
    end

    # =========================================================================
    # Step 3: Filtering (optional) - before frame connection
    # =========================================================================
    if config.filter
        print("  Filtering... ")
        n_before = length(smld.emitters)
        t = @elapsed smld = _filter_smld(smld, config)
        timings["filtering"] = t
        n_after = length(smld.emitters)
        pct = round(100 * n_after / n_before, digits=1)
        println("$n_after / $n_before ($pct%) ($(round(t, digits=2))s)")

        add_step!(workflow, "Filtering", :filter,
            Dict{Symbol,Any}(:min_photons => config.min_photons, :max_precision => config.max_precision, :psf_sigma_mode_tolerance => config.psf_sigma_mode_tolerance),
            :SMLMAnalysis, "$n_after / $n_before accepted")

        # Save filter stats
        if config.outdir !== nothing
            _write_filter_stats(smld_raw, smld, config, t)
        end
    end

    # =========================================================================
    # Step 3.5: Set dataset indices (for multi-dataset drift correction)
    # =========================================================================
    if config.dataset_indices !== nothing
        print("  Setting dataset indices... ")
        for e in smld.emitters
            if 1 <= e.frame <= length(config.dataset_indices)
                e.dataset = config.dataset_indices[e.frame]
            end
        end
        n_datasets = length(unique(config.dataset_indices))
        # Reconstruct SMLD with updated n_datasets (BasicSMLD is immutable)
        smld = BasicSMLD(smld.emitters, smld.camera, smld.n_frames, n_datasets, smld.metadata)
        println("$n_datasets datasets")
    end

    # =========================================================================
    # Step 4: Frame Connection (optional) - connects repeated localizations
    # Runs BEFORE drift correction so combined locs give higher precision for drift
    # =========================================================================
    if config.frameconnect
        print("  Frame connection... ")
        n_before = length(smld.emitters)

        t = @elapsed fc_result = SMLMFrameConnection.frameconnect(smld;
            maxframegap = config.fc_maxframegap,
            nsigmadev = config.fc_nsigmadev,
            nnearestclusters = config.fc_nnearestclusters,
            nmaxnn = config.fc_nmaxnn
        )
        timings["frameconnect"] = t

        smld_connected = fc_result.connected  # Original locs with track_id
        smld = fc_result.combined             # Combined high-precision locs
        fc_params = fc_result.params

        n_after = length(smld.emitters)
        n_tracks = n_after
        compression = round(n_before / n_after, digits=1)
        println("$n_before → $n_tracks tracks ($(compression)x) ($(round(t, digits=2))s)")

        add_step!(workflow, "Frame Connection", :frameconnect,
            Dict{Symbol,Any}(:maxframegap => config.fc_maxframegap, :nsigmadev => config.fc_nsigmadev),
            :SMLMFrameConnection, "$n_tracks tracks from $n_before locs")

        # Analyze drift from linked emitters
        drift_analysis = _analyze_frameconnect_drift(smld_connected)

        # Save frame connection figures and stats
        if config.outdir !== nothing
            _save_frameconnect_figures(smld_connected, config)
            _save_frameconnect_drift_figures(drift_analysis, config)
            _save_uncertainty_calibration_figure(drift_analysis, config)
            _write_frameconnect_stats(n_before, n_after, smld_connected, fc_params, config, t;
                                      drift_analysis=drift_analysis)
        end

        # =====================================================================
        # Step 4b: Uncertainty Calibration (optional) - adjusts uncertainties
        # Uses k and σ_motion from frame connection analysis
        # Regenerates frame-connected results with corrected uncertainties
        # =====================================================================
        if config.calibrate_uncertainties && !isnan(drift_analysis.calibration.A) && !isnan(drift_analysis.calibration.B)
            print("  Uncertainty calibration... ")
            cal = drift_analysis.calibration

            # Extract calibration parameters (stored as variance in nm²)
            # σ_motion² = A (nm²), k² = B
            σ_motion_nm = sqrt(max(0.0, cal.A))  # nm
            k_scale = sqrt(max(1.0, cal.B))      # dimensionless, minimum 1.0
            σ_motion = σ_motion_nm / 1000.0      # convert to μm

            t_cal = @elapsed begin
                # Apply correction to connected localizations and recombine tracks
                smld_connected_corrected, smld_calibrated = _apply_uncertainty_calibration(
                    smld_connected, σ_motion, k_scale)
                smld = smld_calibrated
            end
            timings["calibration"] = t_cal

            println("k=$(round(k_scale, digits=2)), σ_motion=$(round(σ_motion_nm, digits=1))nm ($(round(t_cal, digits=2))s)")

            add_step!(workflow, "Uncertainty Calibration", :calibrate,
                Dict{Symbol,Any}(:k_scale => k_scale, :sigma_motion => σ_motion_nm),
                :SMLMAnalysis, "k=$(round(k_scale, digits=2)), σ_motion=$(round(σ_motion_nm, digits=1))nm")

            # Save calibration stats
            if config.outdir !== nothing
                _write_calibration_stats(σ_motion_nm, k_scale, cal, config, t_cal)
            end
        end
    end

    # =========================================================================
    # Step 5: Drift Correction (optional) - operates on combined locs from FC
    # =========================================================================
    if config.drift
        print("  Drift correction... ")

        t = @elapsed drift_result = driftcorrect(smld;
            degree = config.drift_degree,
            cost_fun = config.drift_cost_fun,
            intramodel = config.drift_model,
            verbose = 0
        )
        timings["drift"] = t

        smld_corrected = drift_result.smld
        drift_model = drift_result.model
        n_datasets = length(drift_model.intra)

        # Extract drift info for reporting
        DC = SMLMDriftCorrection
        n_frames = smld.n_frames
        frames = collect(1:n_frames)

        # Get frame ranges per dataset from emitters (polynomials only valid within their dataset's frames)
        frame_ranges = Dict{Int, Tuple{Int,Int}}()
        for ds in 1:n_datasets
            ds_frames_list = [e.frame for e in smld.emitters if e.dataset == ds]
            if !isempty(ds_frames_list)
                frame_ranges[ds] = (minimum(ds_frames_list), maximum(ds_frames_list))
            end
        end

        # Calculate max drift across all datasets
        max_drift_x = 0.0
        max_drift_y = 0.0
        for ds in 1:n_datasets
            # Only evaluate polynomial over frames that belong to this dataset
            ds_frames = haskey(frame_ranges, ds) ? collect(frame_ranges[ds][1]:frame_ranges[ds][2]) : frames
            drift_x = [DC.applydrift(0.0, f, drift_model.intra[ds].dm[1]) for f in ds_frames]
            drift_y = [DC.applydrift(0.0, f, drift_model.intra[ds].dm[2]) for f in ds_frames]
            max_drift_x = max(max_drift_x, maximum(abs.(drift_x)) * 1000)
            max_drift_y = max(max_drift_y, maximum(abs.(drift_y)) * 1000)
        end

        if n_datasets == 1
            println("max $(round(max_drift_x, digits=1))nm X, $(round(max_drift_y, digits=1))nm Y ($(round(t, digits=2))s)")
        else
            println("$n_datasets datasets, max $(round(max_drift_x, digits=1))nm X, $(round(max_drift_y, digits=1))nm Y ($(round(t, digits=2))s)")
        end

        add_step!(workflow, "Drift Correction", :driftcorrect,
            Dict{Symbol,Any}(:degree => config.drift_degree, :cost_fun => config.drift_cost_fun, :n_datasets => n_datasets),
            :SMLMDriftCorrection, "max $(round(max(max_drift_x, max_drift_y), digits=0))nm")

        # Save drift figures and stats
        if config.outdir !== nothing
            _save_drift_figures(drift_model, smld, config)
            _write_drift_stats(drift_model, smld, config, t)
        end

        smld = smld_corrected
    end

    # =========================================================================
    # Step 6: Isolated Emitter Filter (optional) - runs after drift correction
    # =========================================================================
    if config.filter_isolated
        print("  Isolated filter... ")
        n_before = length(smld.emitters)
        t = @elapsed smld, neighbor_counts, threshold_used = _filter_isolated(smld, config)
        timings["isolated_filter"] = t
        n_after = length(smld.emitters)
        n_rejected = n_before - n_after
        pct_rejected = round(100 * n_rejected / n_before, digits=1)
        auto_str = config.isolated_min_neighbors == :auto ? " (auto)" : ""
        println("$n_rejected rejected ($pct_rejected%) threshold=$threshold_used$auto_str ($(round(t, digits=2))s)")

        add_step!(workflow, "Isolated Filter", :filter_isolated,
            Dict{Symbol,Any}(:n_sigma => config.isolated_n_sigma, :min_neighbors => threshold_used),
            :SMLMAnalysis, "$n_rejected isolated emitters rejected")

        # Save isolated filter figures and stats
        if config.outdir !== nothing
            _save_isolated_figures(neighbor_counts, threshold_used, config)
            _write_isolated_stats(n_before, n_after, neighbor_counts, threshold_used, config, t)
        end
    end

    # =========================================================================
    # Step 7: Rendering (optional) - generates multiple render types
    # =========================================================================
    if config.render
        print("  Rendering... ")
        t_total = 0.0
        render_count = 0

        # Determine clip percentile
        clip_pct = if config.render_clip_percentile === :auto
            _adaptive_clip_percentile(length(smld.emitters))
        else
            config.render_clip_percentile
        end

        if config.outdir !== nothing
            # 1. Gaussian render @ 20x + inferno (primary super-res image)
            if config.render_gaussian
                t = @elapsed render(smld;
                    strategy = GaussianRender(),
                    zoom = config.render_gaussian_zoom,
                    colormap = :inferno,
                    clip_percentile = clip_pct,
                    filename = joinpath(config.outdir, "08_superres", "gaussian_inferno.png")
                )
                t_total += t
                render_count += 1
            end

            # 2. Histogram render @ 10x + time coloring (temporal coverage)
            if config.render_histogram
                t = @elapsed render(smld;
                    strategy = HistogramRender(),
                    zoom = config.render_histogram_zoom,
                    color_by = :frame,
                    colormap = config.render_time_colormap,
                    filename = joinpath(config.outdir, "08_superres", "histogram_time.png")
                )
                t_total += t
                render_count += 1
            end

            # 3. Circles render @ 50x + time coloring (individual loc inspection)
            if config.render_circles
                n_locs = length(smld.emitters)
                if n_locs > 100_000 && config.render_circles_zoom > 30
                    @warn "Large circles render: $(n_locs) locs @ $(config.render_circles_zoom)x may produce large image"
                end
                t = @elapsed render(smld;
                    strategy = CircleRender(),
                    zoom = config.render_circles_zoom,
                    color_by = :frame,
                    colormap = config.render_time_colormap,
                    filename = joinpath(config.outdir, "08_superres", "circles_time.png")
                )
                t_total += t
                render_count += 1
            end

            timings["rendering"] = t_total
            println("$render_count images ($(round(t_total, digits=2))s)")
        else
            # No output directory - just render gaussian to display
            t = @elapsed render(smld;
                strategy = GaussianRender(),
                zoom = config.render_gaussian_zoom,
                colormap = :inferno,
                clip_percentile = clip_pct
            )
            timings["rendering"] = t
            println("($(round(t, digits=2))s)")
        end

        add_step!(workflow, "Rendering", :render,
            Dict{Symbol,Any}(:gaussian => config.render_gaussian, :histogram => config.render_histogram, :circles => config.render_circles),
            :SMLMRender, "$render_count renders")

        # Save render stats
        if config.outdir !== nothing
            _write_render_stats(smld, config, t_total)
        end
    end

    # =========================================================================
    # Step 6: BAGOL (optional, default off)
    # =========================================================================
    if config.bagol
        print("  BAGOL... ")
        # TODO: Integrate SMLMDeepFit
        println("(skipped - package not yet integrated)")
        timings["bagol"] = 0.0
    end

    # =========================================================================
    # Save results
    # =========================================================================
    if config.outdir !== nothing && config.save_smld
        save_smld(joinpath(config.outdir, "results", "smld_final.h5"), smld)
        save_smld(joinpath(config.outdir, "results", "smld_raw.h5"), smld_raw)
        println("  Saved SMLD files to $(config.outdir)/results/")
    end

    # Write summary with health check
    if config.outdir !== nothing
        _write_summary(roi_batch, smld_raw, smld, data, config, timings)
    end

    # Summary
    total = sum(values(timings))
    println("-"^60)
    println("  Total: $(round(total, digits=2))s")
    println("="^60)

    return AnalysisResult(smld, smld_raw, roi_batch, timings, workflow)
end

# =============================================================================
# Helper functions for analyze()
# =============================================================================

"""Adaptive clip percentile based on number of localizations.
Higher percentile for dense data = less clipping = less saturation."""
function _adaptive_clip_percentile(n_locs::Int)
    if n_locs < 50_000
        return 0.99
    elseif n_locs < 200_000
        return 0.995
    elseif n_locs < 500_000
        return 0.999
    else
        return 0.9999  # Very dense - minimal clipping
    end
end

"""Calculate figure size for grid overlay plots based on data aspect ratio."""
function _grid_figure_size(data; n_cols=4, n_rows=3, panel_height=200)
    data_height, data_width = size(data, 1), size(data, 2)
    data_aspect = data_width / data_height
    panel_width = round(Int, panel_height * data_aspect)
    fig_width = panel_width * n_cols + 100   # columns + margins
    fig_height = panel_height * n_rows + 150  # rows + titles/margins
    return (fig_width, fig_height)
end

"""Calculate mode of values using histogram binning."""
function _calculate_mode(values::Vector{T}; n_bins=100) where T<:Real
    if isempty(values)
        return zero(T)
    end

    # Filter out NaN and Inf
    valid = filter(x -> isfinite(x) && x > 0, values)
    if isempty(valid)
        return zero(T)
    end

    # Create histogram bins
    lo, hi = quantile(valid, [0.01, 0.99])  # Robust range
    if lo >= hi
        return median(valid)  # Fall back to median
    end

    edges = range(lo, hi, length=n_bins+1)
    counts = zeros(Int, n_bins)

    for v in valid
        if lo <= v <= hi
            bin_idx = clamp(floor(Int, (v - lo) / (hi - lo) * n_bins) + 1, 1, n_bins)
            counts[bin_idx] += 1
        end
    end

    # Find mode (bin center with highest count)
    mode_idx = argmax(counts)
    mode_value = (edges[mode_idx] + edges[mode_idx+1]) / 2

    return T(mode_value)
end

"""Filter SMLD based on config criteria."""
function _filter_smld(smld::BasicSMLD, config::AnalysisConfig)
    emitters = smld.emitters
    mask = trues(length(emitters))

    if config.min_photons !== nothing
        mask .&= [e.photons > config.min_photons for e in emitters]
    end

    # Precision filter (localization precision from CRLB)
    if config.max_precision !== nothing
        mask .&= [max(e.σ_x, e.σ_y) < config.max_precision for e in emitters]
    end

    # PSF sigma mode filter (for variable sigma fits - filters on fitted PSF width)
    if config.psf_sigma_mode_tolerance !== nothing && length(emitters) > 0
        tol = config.psf_sigma_mode_tolerance

        # Check for isotropic sigma field (GaussianXYNBS)
        if hasfield(typeof(emitters[1]), :σ)
            psf_sigmas = [e.σ for e in emitters]
            psf_sigma_mode = _calculate_mode(psf_sigmas)

            if psf_sigma_mode > 0
                lo = psf_sigma_mode * (1 - tol)
                hi = psf_sigma_mode * (1 + tol)
                mask .&= [lo <= e.σ <= hi for e in emitters]
            end
        # Check for anisotropic sigma fields (GaussianXYNBSXSY)
        elseif hasfield(typeof(emitters[1]), :σx) && hasfield(typeof(emitters[1]), :σy)
            psf_sigmas_x = [e.σx for e in emitters]
            psf_sigmas_y = [e.σy for e in emitters]
            psf_sigma_mode_x = _calculate_mode(psf_sigmas_x)
            psf_sigma_mode_y = _calculate_mode(psf_sigmas_y)

            if psf_sigma_mode_x > 0 && psf_sigma_mode_y > 0
                lo_x = psf_sigma_mode_x * (1 - tol)
                hi_x = psf_sigma_mode_x * (1 + tol)
                lo_y = psf_sigma_mode_y * (1 - tol)
                hi_y = psf_sigma_mode_y * (1 + tol)
                mask .&= [lo_x <= e.σx <= hi_x && lo_y <= e.σy <= hi_y for e in emitters]
            end
        end
    end

    if config.min_pvalue !== nothing
        mask .&= [e.pvalue > config.min_pvalue for e in emitters]
    end

    filtered = emitters[mask]
    return BasicSMLD(filtered, smld.camera, smld.n_frames, smld.n_datasets, smld.metadata)
end

"""
    _filter_isolated(smld, config) -> (filtered_smld, neighbor_counts, threshold)

Filter out isolated emitters based on precision-weighted neighbor counting.
Two emitters are neighbors if their distance < n_sigma × sqrt(σ_i² + σ_j²).

When `isolated_min_neighbors = :auto`, uses triangle method to find optimal threshold.

Returns filtered SMLD, neighbor counts array, and the threshold used.
"""
function _filter_isolated(smld::BasicSMLD, config::AnalysisConfig)
    emitters = smld.emitters
    n = length(emitters)

    if n == 0
        return smld, Int[], 0
    end

    n_sigma = config.isolated_n_sigma

    # Get precision for each emitter (localization uncertainty)
    σ = [sqrt(e.σ_x^2 + e.σ_y^2) for e in emitters]
    max_σ = maximum(σ)
    max_radius = n_sigma * 2 * max_σ  # Conservative search radius for KD-tree

    # Build KD-tree for fast neighbor search
    coords = zeros(2, n)
    for i in 1:n
        coords[1, i] = emitters[i].x
        coords[2, i] = emitters[i].y
    end
    tree = KDTree(coords)

    # Count precision-weighted neighbors for each emitter
    neighbor_counts = zeros(Int, n)
    for i in 1:n
        # Find candidates within max radius
        point = [emitters[i].x, emitters[i].y]
        candidates = inrange(tree, point, max_radius)

        for j in candidates
            j == i && continue
            dist = sqrt((emitters[i].x - emitters[j].x)^2 +
                       (emitters[i].y - emitters[j].y)^2)
            σ_combined = sqrt(σ[i]^2 + σ[j]^2)
            if dist < n_sigma * σ_combined
                neighbor_counts[i] += 1
            end
        end
    end

    # Determine threshold
    if config.isolated_min_neighbors == :auto
        min_neighbors = _triangle_threshold(neighbor_counts)
    else
        min_neighbors = config.isolated_min_neighbors
    end

    # Filter - keep emitters with enough neighbors
    keep = neighbor_counts .>= min_neighbors
    filtered = emitters[keep]

    return BasicSMLD(filtered, smld.camera, smld.n_frames, smld.n_datasets, smld.metadata), neighbor_counts, min_neighbors
end

"""
    _triangle_threshold(counts) -> threshold

Triangle method for automatic threshold selection.

Finds the threshold that maximizes the perpendicular distance from a line
connecting the histogram peak to the histogram tail.

Works well for distributions with a peak at low values and a long tail
(like neighbor count histograms where most noise has 0-1 neighbors).
"""
function _triangle_threshold(counts::Vector{Int})
    if isempty(counts)
        return 1
    end

    max_count = maximum(counts)
    if max_count == 0
        return 1
    end

    # Build histogram (bins 0, 1, 2, ..., max_count)
    hist = zeros(Int, max_count + 1)
    for c in counts
        hist[c + 1] += 1  # +1 for 1-based indexing
    end

    # Find peak (mode) - usually at 0 or low values for noise
    peak_idx = argmax(hist)
    peak_val = hist[peak_idx]

    # Find last non-zero bin (tail end)
    tail_idx = findlast(x -> x > 0, hist)
    if tail_idx === nothing || tail_idx <= peak_idx
        return 1
    end

    # Line from peak to tail: (peak_idx, peak_val) to (tail_idx, hist[tail_idx])
    # Normalized line direction
    dx = tail_idx - peak_idx
    dy = hist[tail_idx] - peak_val
    line_len = sqrt(dx^2 + dy^2)
    if line_len == 0
        return 1
    end

    # Find point with maximum perpendicular distance from line
    max_dist = 0.0
    best_idx = peak_idx

    for i in peak_idx:tail_idx
        # Vector from peak to point i
        px = i - peak_idx
        py = hist[i] - peak_val

        # Perpendicular distance = |cross product| / line_length
        # cross = dx*py - dy*px (z-component of 3D cross product)
        cross = dx * py - dy * px
        dist = abs(cross) / line_len

        if dist > max_dist
            max_dist = dist
            best_idx = i
        end
    end

    # Return threshold (convert from 1-based histogram index to neighbor count)
    threshold = best_idx - 1  # Convert back to 0-based count

    # Ensure minimum threshold of 1 (at least 1 neighbor required)
    return max(1, threshold)
end

"""
    _apply_uncertainty_calibration(smld_connected, σ_motion, k_scale)

Apply uncertainty calibration to connected localizations and recombine tracks.

The calibration model is: σ²_corrected = σ²_motion + k² × σ²_CRLB

For each localization:
- σ_x_corrected = √(σ_motion² + k² × σ_x²)
- σ_y_corrected = √(σ_motion² + k² × σ_y²)

Then recombines tracks using weighted averaging with corrected uncertainties.

Returns:
- smld_connected_corrected: Connected localizations with corrected uncertainties
- smld_combined: Recombined tracks with corrected uncertainties
"""
function _apply_uncertainty_calibration(smld_connected::BasicSMLD, σ_motion::Float64, k_scale::Float64)
    emitters = smld_connected.emitters

    # Apply correction to each localization
    σ_motion_sq = σ_motion^2
    k_sq = k_scale^2

    corrected_emitters = map(emitters) do e
        # Corrected uncertainties: σ_corrected = √(σ_motion² + k² × σ_CRLB²)
        σ_x_corrected = Float32(sqrt(σ_motion_sq + k_sq * e.σ_x^2))
        σ_y_corrected = Float32(sqrt(σ_motion_sq + k_sq * e.σ_y^2))

        # Create new emitter with corrected uncertainties using setproperties pattern
        _copy_emitter_with_uncertainty(e, σ_x_corrected, σ_y_corrected)
    end

    smld_connected_corrected = BasicSMLD(corrected_emitters, smld_connected.camera,
                                         smld_connected.n_frames, smld_connected.n_datasets,
                                         smld_connected.metadata)

    # Recombine tracks with corrected uncertainties
    smld_combined = _recombine_tracks(smld_connected_corrected)

    return smld_connected_corrected, smld_combined
end

"""Copy emitter with new σ_x and σ_y values. Works for any GaussMLE emitter type."""
function _copy_emitter_with_uncertainty(e::GaussMLE.Emitter2DFitGaussMLE, σ_x::Float32, σ_y::Float32)
    GaussMLE.Emitter2DFitGaussMLE(
        e.x, e.y, e.photons, e.bg, σ_x, σ_y, e.σ_photons, e.σ_bg,
        e.pvalue, e.frame, e.dataset, e.track_id, e.id
    )
end

function _copy_emitter_with_uncertainty(e::GaussMLE.Emitter2DFitSigma, σ_x::Float32, σ_y::Float32)
    GaussMLE.Emitter2DFitSigma(
        e.x, e.y, e.photons, e.bg, e.σ, σ_x, σ_y, e.σ_photons, e.σ_bg, e.σ_σ,
        e.pvalue, e.frame, e.dataset, e.track_id, e.id
    )
end

function _copy_emitter_with_uncertainty(e::GaussMLE.Emitter2DFitSigmaXY, σ_x::Float32, σ_y::Float32)
    GaussMLE.Emitter2DFitSigmaXY(
        e.x, e.y, e.photons, e.bg, e.σx, e.σy, σ_x, σ_y, e.σ_photons, e.σ_bg, e.σ_σx, e.σ_σy,
        e.pvalue, e.frame, e.dataset, e.track_id, e.id
    )
end

# Fallback for SMLMData Emitter2DFit (no pvalue field)
function _copy_emitter_with_uncertainty(e::SMLMData.Emitter2DFit, σ_x::Float32, σ_y::Float32)
    SMLMData.Emitter2DFit(
        e.x, e.y, e.photons, e.bg, σ_x, σ_y, e.σ_photons, e.σ_bg,
        e.frame, e.dataset, e.track_id, e.id
    )
end

"""
    _recombine_tracks(smld_connected)

Recombine localizations into tracks using weighted averaging.

Groups localizations by track_id and computes weighted average position
and combined uncertainty for each track.

The combined emitter uses:
- Position: weighted average (weight = 1/σ²)
- Uncertainty: σ_combined = 1/√(Σ 1/σ²)
- Frame: middle frame of track
- Photons: sum of photons
- Background: mean background
- pvalue: geometric mean of pvalues
"""
function _recombine_tracks(smld_connected::BasicSMLD)
    emitters = smld_connected.emitters
    EmitterType = eltype(emitters)

    # Group by track_id
    track_dict = Dict{Int, Vector{EmitterType}}()
    for e in emitters
        if e.track_id > 0
            if !haskey(track_dict, e.track_id)
                track_dict[e.track_id] = EmitterType[]
            end
            push!(track_dict[e.track_id], e)
        end
    end

    # Combine each track
    combined_emitters = EmitterType[]

    for (track_id, track_locs) in track_dict
        n = length(track_locs)

        # Weighted average position
        sum_wx, sum_wy = 0.0, 0.0
        sum_w_x, sum_w_y = 0.0, 0.0
        sum_photons = 0.0
        sum_bg = 0.0
        log_pvalue_sum = 0.0
        min_frame, max_frame = typemax(Int), 0
        dataset = track_locs[1].dataset

        for e in track_locs
            # Weights = 1/variance
            w_x = 1.0 / e.σ_x^2
            w_y = 1.0 / e.σ_y^2

            sum_wx += w_x * e.x
            sum_wy += w_y * e.y
            sum_w_x += w_x
            sum_w_y += w_y

            sum_photons += e.photons
            sum_bg += e.bg
            if hasproperty(e, :pvalue)
                log_pvalue_sum += log(max(1e-300, e.pvalue))
            end

            min_frame = min(min_frame, e.frame)
            max_frame = max(max_frame, e.frame)
        end

        # Combined position and uncertainty
        x_combined = Float32(sum_wx / sum_w_x)
        y_combined = Float32(sum_wy / sum_w_y)
        σ_x_combined = Float32(1.0 / sqrt(sum_w_x))
        σ_y_combined = Float32(1.0 / sqrt(sum_w_y))

        # Other combined properties
        frame_combined = div(min_frame + max_frame, 2)
        bg_combined = Float32(sum_bg / n)
        pvalue_combined = Float32(exp(log_pvalue_sum / n))  # geometric mean

        # Create combined emitter using first emitter as template
        template = track_locs[1]
        combined_e = _create_combined_emitter(template, x_combined, y_combined,
                                               Float32(sum_photons), bg_combined,
                                               σ_x_combined, σ_y_combined,
                                               frame_combined, dataset, track_id, pvalue_combined)
        push!(combined_emitters, combined_e)
    end

    return BasicSMLD(combined_emitters, smld_connected.camera,
                     smld_connected.n_frames, smld_connected.n_datasets,
                     smld_connected.metadata)
end

"""Create combined emitter from weighted average. Works for any GaussMLE emitter type."""
function _create_combined_emitter(template::GaussMLE.Emitter2DFitGaussMLE,
                                   x, y, photons, bg, σ_x, σ_y, frame, dataset, track_id, pvalue)
    GaussMLE.Emitter2DFitGaussMLE(x, y, photons, bg, σ_x, σ_y, σ_x, σ_y,  # σ_photons/σ_bg approximated
                                  pvalue, frame, dataset, track_id, 0)
end

function _create_combined_emitter(template::GaussMLE.Emitter2DFitSigma,
                                   x, y, photons, bg, σ_x, σ_y, frame, dataset, track_id, pvalue)
    σ_mean = (σ_x + σ_y) / 2
    GaussMLE.Emitter2DFitSigma(x, y, photons, bg, σ_mean, σ_x, σ_y, σ_x, σ_y, σ_mean,
                               pvalue, frame, dataset, track_id, 0)
end

function _create_combined_emitter(template::GaussMLE.Emitter2DFitSigmaXY,
                                   x, y, photons, bg, σ_x, σ_y, frame, dataset, track_id, pvalue)
    # For anisotropic, use mean of template's PSF sigmas as representative
    σx_psf = (template.σx + σ_x) / 2  # Approximate combined PSF sigma
    σy_psf = (template.σy + σ_y) / 2
    GaussMLE.Emitter2DFitSigmaXY(x, y, photons, bg, σx_psf, σy_psf, σ_x, σ_y,
                                  σ_x, σ_y, σ_x, σ_y,  # Approximations for σ_photons, etc.
                                  pvalue, frame, dataset, track_id, 0)
end

# Fallback for SMLMData Emitter2DFit
function _create_combined_emitter(template::SMLMData.Emitter2DFit,
                                   x, y, photons, bg, σ_x, σ_y, frame, dataset, track_id, pvalue)
    SMLMData.Emitter2DFit(x, y, photons, bg, σ_x, σ_y, σ_x, σ_y,
                          frame, dataset, track_id, 0)
end

"""Write calibration statistics to markdown file."""
function _write_calibration_stats(σ_motion_nm::Float64, k_scale::Float64, cal, config::AnalysisConfig, elapsed::Float64)
    filepath = joinpath(config.outdir, "05_calibration", "calibration_stats.md")

    open(filepath, "w") do io
        println(io, "# Uncertainty Calibration Applied")
        println(io)
        println(io, "## Calibration Parameters")
        println(io, "| Parameter | Value | Description |")
        println(io, "|-----------|-------|-------------|")
        println(io, "| σ_motion | $(round(σ_motion_nm, digits=1)) nm | Frame-to-frame motion/vibration |")
        println(io, "| k (CRLB scale) | $(round(k_scale, digits=2)) | Multiply σ_CRLB by this |")
        println(io)
        println(io, "## Correction Formula")
        println(io, "```")
        println(io, "σ_corrected = √(σ_motion² + k² × σ_CRLB²)")
        println(io, "```")
        println(io)
        println(io, "## Effect")
        println(io, "- All localization uncertainties have been adjusted using this model")
        println(io, "- Frame-connected track uncertainties recalculated with corrected weights")
        println(io, "- Time: $(round(elapsed, digits=2))s")
    end
end

"""Save detection overlay figures."""
function _save_detection_figures(data, roi_batch, camera, config)
    nframes = size(data, 3)
    frame_indices = [round(Int, x) for x in range(1, nframes, length=12)]

    # Intensity range
    pmin = Float64(quantile(vec(data[:,:,1]), 0.01))
    pmax = Float64(quantile(vec(data[:,:,1]), 0.99))

    fig = Figure(size=_grid_figure_size(data))
    box_size = roi_batch.roi_size

    for (idx, frame_num) in enumerate(frame_indices)
        row = div(idx - 1, 4) + 1
        col = mod(idx - 1, 4) + 1

        ax = Axis(fig[row, col],
            title = "Frame $frame_num",
            aspect = DataAspect(),
            yreversed = true
        )

        frame_data = data[:, :, frame_num]'
        heatmap!(ax, frame_data, colormap=:grays, colorrange=(pmin, pmax))

        frame_mask = roi_batch.frame_indices .== frame_num
        if any(frame_mask)
            det_x = roi_batch.x_corners[frame_mask]
            det_y = roi_batch.y_corners[frame_mask]
            for (x, y) in zip(det_x, det_y)
                lines!(ax, [x, x+box_size, x+box_size, x, x],
                          [y, y, y+box_size, y+box_size, y],
                    color=:yellow, linewidth=0.5)
            end
        end
        hidedecorations!(ax)
    end

    save(joinpath(config.outdir, "01_detection", "detection_overlay.png"), fig)
end

"""Save fitting quality figures."""
function _save_fitting_figures(smld, roi_batch, data, camera, config)
    emitters = smld.emitters
    photons = [e.photons for e in emitters]
    bg = [e.bg for e in emitters]
    σ_x = [e.σ_x for e in emitters]
    σ_y = [e.σ_y for e in emitters]
    pvalue = [e.pvalue for e in emitters]

    # Fit quality histograms (6 panels, 3x2 grid)
    fig = Figure(size=(1600, 1200))

    # Photons histogram (0 to min(1e5, max_photons))
    photon_max = min(1e5, maximum(photons))
    ax1 = Axis(fig[1, 1], xlabel="Photons", ylabel="Count", title="Photon Distribution",
               limits=(0, photon_max, nothing, nothing))
    hist!(ax1, clamp.(photons, 0, photon_max), bins=range(0, photon_max, length=51))
    vlines!(ax1, [median(photons)], color=:red, linestyle=:dash)

    # Background histogram
    ax2 = Axis(fig[1, 2], xlabel="Background (ADU)", ylabel="Count", title="Background Distribution")
    hist!(ax2, bg, bins=50)
    vlines!(ax2, [median(bg)], color=:red, linestyle=:dash)

    # Precision histograms (σ_x, σ_y) - fixed range 0-50nm
    σ_x_nm = σ_x .* 1000
    σ_y_nm = σ_y .* 1000

    ax3 = Axis(fig[2, 1], xlabel="σ_x (nm)", ylabel="Count", title="X Precision Distribution",
               limits=(0, 50, nothing, nothing))
    hist!(ax3, clamp.(σ_x_nm, 0, 50), bins=range(0, 50, length=51))
    vlines!(ax3, [median(σ_x_nm)], color=:red, linestyle=:dash, label="Median")
    # Show precision threshold
    if config.max_precision !== nothing
        vlines!(ax3, [config.max_precision * 1000], color=:orange, linestyle=:solid, label="Threshold")
    end

    ax4 = Axis(fig[2, 2], xlabel="σ_y (nm)", ylabel="Count", title="Y Precision Distribution",
               limits=(0, 50, nothing, nothing))
    hist!(ax4, clamp.(σ_y_nm, 0, 50), bins=range(0, 50, length=51))
    vlines!(ax4, [median(σ_y_nm)], color=:red, linestyle=:dash)
    if config.max_precision !== nothing
        vlines!(ax4, [config.max_precision * 1000], color=:orange, linestyle=:solid)
    end

    # p-value histogram - LOG SCALE x-axis
    ax5 = Axis(fig[3, 1], xlabel="log₁₀(p-value)", ylabel="Count", title="Fit Quality (p-value)")
    pval_nonzero = pvalue[pvalue .> 0]
    n_zero = sum(pvalue .== 0)
    if isempty(pval_nonzero)
        text!(ax5, 0.5, 0.5, text="All p-values = 0\n(check PSF model)", align=(:center, :center),
              space=:relative)
    else
        log_pval = log10.(pval_nonzero)
        hist!(ax5, log_pval, bins=50)
        # Show threshold on log scale
        if config.min_pvalue !== nothing && config.min_pvalue > 0
            vlines!(ax5, [log10(config.min_pvalue)], color=:orange, linestyle=:solid)
        end
        # Annotate how many are zero
        if n_zero > 0
            text!(ax5, 0.02, 0.98, text="$(n_zero) fits with pvalue=0 ($(round(100*n_zero/length(pvalue), digits=1))%)",
                  align=(:left, :top), space=:relative, fontsize=10)
        end
    end

    # PSF Sigma histogram (for variable sigma fits)
    has_sigma = hasfield(typeof(emitters[1]), :σ)
    has_sigma_xy = hasfield(typeof(emitters[1]), :σx) && hasfield(typeof(emitters[1]), :σy)

    # Initialize variables for acceptance calculation
    psf_sigma_mode = 0.0
    psf_sigma_mode_x = 0.0
    psf_sigma_mode_y = 0.0

    if has_sigma
        psf_sigmas = [e.σ for e in emitters]
        psf_sigma_nm = psf_sigmas .* 1000
        psf_sigma_mode = _calculate_mode(psf_sigmas)
        psf_sigma_mode_nm = psf_sigma_mode * 1000

        ax6 = Axis(fig[3, 2], xlabel="PSF σ (nm)", ylabel="Count", title="PSF Sigma Distribution",
                   limits=(50, 300, nothing, nothing))
        hist!(ax6, clamp.(psf_sigma_nm, 50, 300), bins=range(50, 300, length=51))
        vlines!(ax6, [median(psf_sigma_nm)], color=:red, linestyle=:dash, label="Median")
        if config.psf_sigma_mode_tolerance !== nothing && psf_sigma_mode > 0
            tol = config.psf_sigma_mode_tolerance
            lo_nm = psf_sigma_mode_nm * (1 - tol)
            hi_nm = psf_sigma_mode_nm * (1 + tol)
            vlines!(ax6, [psf_sigma_mode_nm], color=:blue, linestyle=:solid, label="Mode")
            vlines!(ax6, [lo_nm, hi_nm], color=:orange, linestyle=:solid, label="Bounds")
        end
    elseif has_sigma_xy
        # Anisotropic PSF sigma - show both σx and σy on same plot
        psf_sigmas_x = [e.σx for e in emitters]
        psf_sigmas_y = [e.σy for e in emitters]
        psf_sigma_mode_x = _calculate_mode(psf_sigmas_x)
        psf_sigma_mode_y = _calculate_mode(psf_sigmas_y)

        ax6 = Axis(fig[3, 2], xlabel="PSF σ (nm)", ylabel="Count", title="PSF Sigma Distribution (σx=blue, σy=red)",
                   limits=(50, 300, nothing, nothing))
        hist!(ax6, clamp.(psf_sigmas_x .* 1000, 50, 300), bins=range(50, 300, length=51), color=(:blue, 0.5), label="σx")
        hist!(ax6, clamp.(psf_sigmas_y .* 1000, 50, 300), bins=range(50, 300, length=51), color=(:red, 0.5), label="σy")
        if config.psf_sigma_mode_tolerance !== nothing
            tol = config.psf_sigma_mode_tolerance
            if psf_sigma_mode_x > 0
                lo_x = psf_sigma_mode_x * 1000 * (1 - tol)
                hi_x = psf_sigma_mode_x * 1000 * (1 + tol)
                vlines!(ax6, [psf_sigma_mode_x * 1000], color=:blue, linestyle=:solid)
                vlines!(ax6, [lo_x, hi_x], color=:blue, linestyle=:dash)
            end
            if psf_sigma_mode_y > 0
                lo_y = psf_sigma_mode_y * 1000 * (1 - tol)
                hi_y = psf_sigma_mode_y * 1000 * (1 + tol)
                vlines!(ax6, [psf_sigma_mode_y * 1000], color=:red, linestyle=:solid)
                vlines!(ax6, [lo_y, hi_y], color=:red, linestyle=:dash)
            end
        end
    else
        # For fixed sigma fits, show photons vs background
        ax6 = Axis(fig[3, 2], xlabel="Photons", ylabel="Background (ADU)", title="Photons vs Background")
        scatter!(ax6, photons, bg, markersize=2, alpha=0.3)
    end

    save(joinpath(config.outdir, "02_fitting", "fit_quality.png"), fig)

    # Fit acceptance panel
    precision_values = [sqrt(e.σ_x^2 + e.σ_y^2)/sqrt(2) for e in emitters]

    photon_ok = config.min_photons === nothing ? trues(length(emitters)) : photons .> config.min_photons

    # Precision filter
    precision_ok = config.max_precision === nothing ? trues(length(emitters)) :
                   [max(e.σ_x, e.σ_y) < config.max_precision for e in emitters]

    # PSF sigma mode filter (for variable sigma fits)
    if has_sigma && config.psf_sigma_mode_tolerance !== nothing && psf_sigma_mode > 0
        tol = config.psf_sigma_mode_tolerance
        lo = psf_sigma_mode * (1 - tol)
        hi = psf_sigma_mode * (1 + tol)
        psf_sigma_ok = [lo <= e.σ <= hi for e in emitters]
    elseif has_sigma_xy && config.psf_sigma_mode_tolerance !== nothing && psf_sigma_mode_x > 0 && psf_sigma_mode_y > 0
        tol = config.psf_sigma_mode_tolerance
        lo_x = psf_sigma_mode_x * (1 - tol)
        hi_x = psf_sigma_mode_x * (1 + tol)
        lo_y = psf_sigma_mode_y * (1 - tol)
        hi_y = psf_sigma_mode_y * (1 + tol)
        psf_sigma_ok = [lo_x <= e.σx <= hi_x && lo_y <= e.σy <= hi_y for e in emitters]
    else
        psf_sigma_ok = trues(length(emitters))
    end

    pvalue_ok = config.min_pvalue === nothing ? trues(length(emitters)) : pvalue .> config.min_pvalue
    accepted = photon_ok .& precision_ok .& psf_sigma_ok .& pvalue_ok

    n_total = length(emitters)
    n_accepted = sum(accepted)
    accept_pct = round(100 * n_accepted / n_total, digits=1)

    nframes = size(data, 3)
    frame_indices = [round(Int, x) for x in range(1, nframes, length=12)]
    pmin = Float64(quantile(vec(data[:,:,1]), 0.01))
    pmax = Float64(quantile(vec(data[:,:,1]), 0.99))

    fig = Figure(size=_grid_figure_size(data))
    Label(fig[0, 1:4],
        "Fit Acceptance: $n_accepted/$n_total ($accept_pct%) — Green=Accepted, Red=Photons, Orange=σ, Purple=pvalue",
        fontsize=14, tellwidth=false)

    box_size = roi_batch.roi_size

    for (idx, frame_num) in enumerate(frame_indices)
        row = div(idx - 1, 4) + 1
        col = mod(idx - 1, 4) + 1

        frame_mask = roi_batch.frame_indices .== frame_num
        n_in_frame = sum(frame_mask)
        n_acc_frame = sum(accepted[frame_mask])

        ax = Axis(fig[row, col],
            title = "Frame $frame_num ($n_acc_frame/$n_in_frame)",
            aspect = DataAspect(),
            yreversed = true
        )

        frame_data = data[:, :, frame_num]'
        heatmap!(ax, frame_data, colormap=:grays, colorrange=(pmin, pmax))

        frame_locs = findall(frame_mask)
        if !isempty(frame_locs)
            det_x = roi_batch.x_corners[frame_mask]
            det_y = roi_batch.y_corners[frame_mask]
            frame_accepted = accepted[frame_mask]
            frame_photon_ok = photon_ok[frame_mask]
            frame_precision_ok = precision_ok[frame_mask]
            frame_psf_sigma_ok = psf_sigma_ok[frame_mask]
            frame_pvalue_ok = pvalue_ok[frame_mask]

            for pass in [false, true]
                for j in eachindex(frame_locs)
                    if frame_accepted[j] == pass
                        bx, by = det_x[j], det_y[j]
                        # Color: green=accepted, red=photons, orange=precision/psf_sigma, purple=pvalue
                        c = if frame_accepted[j]
                            :green
                        elseif !frame_photon_ok[j]
                            :red
                        elseif !frame_precision_ok[j] || !frame_psf_sigma_ok[j]
                            :orange
                        else
                            :purple  # pvalue failed
                        end
                        lines!(ax, [bx, bx+box_size, bx+box_size, bx, bx],
                                  [by, by, by+box_size, by+box_size, by],
                            color = (c, pass ? 1.0 : 0.7), linewidth = 0.5)
                    end
                end
            end
        end
        hidedecorations!(ax)
    end

    save(joinpath(config.outdir, "02_fitting", "fit_acceptance.png"), fig)
end

# =============================================================================
# Stats Output Functions - Markdown files for human and LLM diagnostics
# =============================================================================

"""Write detection statistics markdown file."""
function _write_detection_stats(roi_batch, data, config, elapsed_time)
    nframes = size(data, 3)
    n_rois = length(roi_batch)

    # ROIs per frame
    rois_per_frame = [sum(roi_batch.frame_indices .== f) for f in 1:nframes]

    # Intensity stats from ROIs
    roi_intensities = [maximum(roi_batch.data[:,:,i]) for i in 1:min(n_rois, 10000)]
    bg_estimates = [minimum(roi_batch.data[:,:,i]) for i in 1:min(n_rois, 10000)]
    signal_above_bg = roi_intensities .- bg_estimates

    filepath = joinpath(config.outdir, "01_detection", "detection_stats.md")
    open(filepath, "w") do io
        println(io, "# Detection Statistics\n")
        println(io, "## Summary")
        println(io, "- **Total ROIs detected**: $(n_rois)")
        println(io, "- **Frames analyzed**: $(nframes)")
        println(io, "- **Time**: $(round(elapsed_time, digits=2))s ($(round(nframes/elapsed_time, digits=0)) frames/s)")
        println(io, "- **ROIs per frame**: mean=$(round(mean(rois_per_frame), digits=1)), std=$(round(std(rois_per_frame), digits=1))")
        println(io, "- **ROIs/frame range**: $(minimum(rois_per_frame)) - $(maximum(rois_per_frame))")
        println(io, "")
        println(io, "## Detection Parameters")
        println(io, "- Box size: $(config.boxsize)")
        println(io, "- Min photons threshold: $(config.detect_min_photons)")
        println(io, "- PSF sigma: $(round(config.psf_sigma * 1000, digits=0)) nm")
        println(io, "- Overlap: $(config.overlap)")
        println(io, "")
        println(io, "## ROI Intensity (sampled)")
        println(io, "- Peak intensity: median=$(round(median(roi_intensities), digits=0)) ADU")
        println(io, "- Background estimate: median=$(round(median(bg_estimates), digits=0)) ADU")
        println(io, "- Signal above background: median=$(round(median(signal_above_bg), digits=0)) ADU")
        println(io, "")
        println(io, "## Health Check")
        rpf_mean = mean(rois_per_frame)
        rpf_cv = std(rois_per_frame) / rpf_mean
        println(io, "- ROIs/frame CV: $(round(rpf_cv, digits=2)) ", rpf_cv < 0.3 ? "✓" : "⚠ (high variation)")
        println(io, "- Signal/background: $(round(median(signal_above_bg)/median(bg_estimates), digits=1))x ",
                median(signal_above_bg) > median(bg_estimates) ? "✓" : "⚠")
    end
end

"""Write fitting statistics markdown file."""
function _write_fitting_stats(smld_raw, config, elapsed_time)
    emitters = smld_raw.emitters
    n = length(emitters)

    photons = [e.photons for e in emitters]
    bg = [e.bg for e in emitters]
    σ_x = [e.σ_x for e in emitters]
    σ_y = [e.σ_y for e in emitters]
    pvalue = [e.pvalue for e in emitters]

    # Check if we have fitted PSF sigma (GaussianXYNBS or GaussianXYNBSXSY)
    has_sigma = hasfield(typeof(emitters[1]), :σ)
    has_sigma_xy = hasfield(typeof(emitters[1]), :σx)

    filepath = joinpath(config.outdir, "02_fitting", "fitting_stats.md")
    open(filepath, "w") do io
        println(io, "# Fitting Statistics\n")
        println(io, "## Summary")
        println(io, "- **Total fits**: $(n)")
        println(io, "- **Time**: $(round(elapsed_time, digits=2))s ($(round(n/elapsed_time/1000, digits=0))k fits/s)")
        println(io, "- **Model**: $(config.fit_model)")
        println(io, "- **Iterations**: $(config.iterations)")
        println(io, "")

        println(io, "## Parameter Distributions\n")
        println(io, "| Parameter | 5% | 25% | 50% (median) | 75% | 95% |")
        println(io, "|-----------|-----|-----|--------------|-----|-----|")

        for (name, vals, unit, scale) in [
            ("Photons", photons, "", 1),
            ("Background", bg, "e⁻", 1),
            ("Precision σ_x", σ_x, "nm", 1000),
            ("Precision σ_y", σ_y, "nm", 1000),
        ]
            p = quantile(vals .* scale, [0.05, 0.25, 0.50, 0.75, 0.95])
            println(io, "| $name ($unit) | $(round(p[1], digits=1)) | $(round(p[2], digits=1)) | $(round(p[3], digits=1)) | $(round(p[4], digits=1)) | $(round(p[5], digits=1)) |")
        end

        if has_sigma
            σ = [e.σ for e in emitters]
            p = quantile(σ .* 1000, [0.05, 0.25, 0.50, 0.75, 0.95])
            println(io, "| PSF σ (nm) | $(round(p[1], digits=1)) | $(round(p[2], digits=1)) | $(round(p[3], digits=1)) | $(round(p[4], digits=1)) | $(round(p[5], digits=1)) |")
        end

        if has_sigma_xy
            σx = [e.σx for e in emitters]
            σy = [e.σy for e in emitters]
            px = quantile(σx .* 1000, [0.05, 0.25, 0.50, 0.75, 0.95])
            py = quantile(σy .* 1000, [0.05, 0.25, 0.50, 0.75, 0.95])
            println(io, "| PSF σx (nm) | $(round(px[1], digits=1)) | $(round(px[2], digits=1)) | $(round(px[3], digits=1)) | $(round(px[4], digits=1)) | $(round(px[5], digits=1)) |")
            println(io, "| PSF σy (nm) | $(round(py[1], digits=1)) | $(round(py[2], digits=1)) | $(round(py[3], digits=1)) | $(round(py[4], digits=1)) | $(round(py[5], digits=1)) |")
        end

        println(io, "")
        println(io, "## P-value Distribution")
        n_zero = sum(pvalue .== 0)
        println(io, "- pvalue = 0: $(n_zero) ($(round(100*n_zero/n, digits=1))%)")
        println(io, "- pvalue > 0: $(sum(pvalue .> 0)) ($(round(100*sum(pvalue .> 0)/n, digits=1))%)")
        println(io, "- pvalue > 0.001: $(sum(pvalue .> 0.001)) ($(round(100*sum(pvalue .> 0.001)/n, digits=1))%)")
        println(io, "- pvalue > 0.01: $(sum(pvalue .> 0.01)) ($(round(100*sum(pvalue .> 0.01)/n, digits=1))%)")
        println(io, "- pvalue > 0.05: $(sum(pvalue .> 0.05)) ($(round(100*sum(pvalue .> 0.05)/n, digits=1))%)")
        println(io, "")
        println(io, "## Health Check")
        pval_pass = sum(pvalue .> 0.001) / n
        prec_med = median(σ_x) * 1000
        println(io, "- Precision σ_x median: $(round(prec_med, digits=1)) nm ", prec_med < 20 ? "✓" : "⚠")
        println(io, "- pvalue > 0.001: $(round(100*pval_pass, digits=1))% ", pval_pass > 0.05 ? "✓" : "⚠ (low)")
        println(io, "- Photons median: $(round(median(photons), digits=0)) ", 1000 < median(photons) < 100000 ? "✓" : "⚠")
    end
end

"""Write filter statistics markdown file."""
function _write_filter_stats(smld_raw, smld_filtered, config, elapsed_time)
    n_raw = length(smld_raw.emitters)
    n_filtered = length(smld_filtered.emitters)

    emitters = smld_raw.emitters
    photons = [e.photons for e in emitters]
    σ_x = [e.σ_x for e in emitters]
    σ_y = [e.σ_y for e in emitters]
    pvalue = [e.pvalue for e in emitters]

    # Check for PSF sigma fields
    has_sigma = hasfield(typeof(emitters[1]), :σ)
    has_sigma_xy = hasfield(typeof(emitters[1]), :σx) && hasfield(typeof(emitters[1]), :σy)

    # Calculate modes for sigma fields
    psf_sigma_mode = 0.0
    psf_sigma_mode_x = 0.0
    psf_sigma_mode_y = 0.0

    if has_sigma
        psf_sigmas = [e.σ for e in emitters]
        psf_sigma_mode = _calculate_mode(psf_sigmas)
    elseif has_sigma_xy
        psf_sigmas_x = [e.σx for e in emitters]
        psf_sigmas_y = [e.σy for e in emitters]
        psf_sigma_mode_x = _calculate_mode(psf_sigmas_x)
        psf_sigma_mode_y = _calculate_mode(psf_sigmas_y)
    end

    # Calculate per-filter pass rates
    photon_ok = config.min_photons === nothing ? trues(n_raw) : photons .> config.min_photons

    # Precision filter
    precision_ok = config.max_precision === nothing ? trues(n_raw) :
                   [max(e.σ_x, e.σ_y) < config.max_precision for e in emitters]

    # PSF sigma mode filter
    if has_sigma && config.psf_sigma_mode_tolerance !== nothing && psf_sigma_mode > 0
        tol = config.psf_sigma_mode_tolerance
        lo = psf_sigma_mode * (1 - tol)
        hi = psf_sigma_mode * (1 + tol)
        psf_sigma_ok = [lo <= e.σ <= hi for e in emitters]
    elseif has_sigma_xy && config.psf_sigma_mode_tolerance !== nothing && psf_sigma_mode_x > 0 && psf_sigma_mode_y > 0
        tol = config.psf_sigma_mode_tolerance
        lo_x = psf_sigma_mode_x * (1 - tol)
        hi_x = psf_sigma_mode_x * (1 + tol)
        lo_y = psf_sigma_mode_y * (1 - tol)
        hi_y = psf_sigma_mode_y * (1 + tol)
        psf_sigma_ok = [lo_x <= e.σx <= hi_x && lo_y <= e.σy <= hi_y for e in emitters]
    else
        psf_sigma_ok = trues(n_raw)
    end

    pvalue_ok = config.min_pvalue === nothing ? trues(n_raw) : pvalue .> config.min_pvalue

    n_photon_pass = sum(photon_ok)
    n_precision_pass = sum(precision_ok)
    n_psf_sigma_pass = sum(psf_sigma_ok)
    n_pvalue_pass = sum(pvalue_ok)

    filepath = joinpath(config.outdir, "03_filtered", "filter_stats.md")
    open(filepath, "w") do io
        println(io, "# Filter Statistics\n")
        println(io, "## Summary")
        println(io, "- **Input**: $(n_raw) localizations")
        println(io, "- **Output**: $(n_filtered) localizations")
        println(io, "- **Acceptance rate**: $(round(100*n_filtered/n_raw, digits=1))%")
        println(io, "- **Time**: $(round(elapsed_time*1000, digits=1))ms")
        println(io, "")

        # Report PSF sigma mode for variable sigma fits
        if has_sigma && config.psf_sigma_mode_tolerance !== nothing && psf_sigma_mode > 0
            psf_sigma_mode_nm = psf_sigma_mode * 1000
            println(io, "## PSF Sigma Mode Analysis")
            println(io, "- **PSF sigma mode**: $(round(psf_sigma_mode_nm, digits=1)) nm")
            tol_pct = round(config.psf_sigma_mode_tolerance * 100, digits=0)
            lo_nm = round(psf_sigma_mode_nm * (1 - config.psf_sigma_mode_tolerance), digits=1)
            hi_nm = round(psf_sigma_mode_nm * (1 + config.psf_sigma_mode_tolerance), digits=1)
            println(io, "- **Acceptance range**: $(lo_nm) - $(hi_nm) nm (mode ± $(tol_pct)%)")
            println(io, "")
        elseif has_sigma_xy && config.psf_sigma_mode_tolerance !== nothing && psf_sigma_mode_x > 0 && psf_sigma_mode_y > 0
            psf_sigma_mode_x_nm = psf_sigma_mode_x * 1000
            psf_sigma_mode_y_nm = psf_sigma_mode_y * 1000
            println(io, "## PSF Sigma Mode Analysis (Anisotropic)")
            println(io, "- **PSF σx mode**: $(round(psf_sigma_mode_x_nm, digits=1)) nm")
            println(io, "- **PSF σy mode**: $(round(psf_sigma_mode_y_nm, digits=1)) nm")
            tol_pct = round(config.psf_sigma_mode_tolerance * 100, digits=0)
            lo_x_nm = round(psf_sigma_mode_x_nm * (1 - config.psf_sigma_mode_tolerance), digits=1)
            hi_x_nm = round(psf_sigma_mode_x_nm * (1 + config.psf_sigma_mode_tolerance), digits=1)
            lo_y_nm = round(psf_sigma_mode_y_nm * (1 - config.psf_sigma_mode_tolerance), digits=1)
            hi_y_nm = round(psf_sigma_mode_y_nm * (1 + config.psf_sigma_mode_tolerance), digits=1)
            println(io, "- **σx acceptance range**: $(lo_x_nm) - $(hi_x_nm) nm (mode ± $(tol_pct)%)")
            println(io, "- **σy acceptance range**: $(lo_y_nm) - $(hi_y_nm) nm (mode ± $(tol_pct)%)")
            println(io, "")
        end

        println(io, "## Per-Filter Results\n")
        println(io, "| Filter | Threshold | Pass | Fail | % Pass |")
        println(io, "|--------|-----------|------|------|--------|")

        if config.min_photons !== nothing
            println(io, "| Photons | >$(config.min_photons) | $(n_photon_pass) | $(n_raw - n_photon_pass) | $(round(100*n_photon_pass/n_raw, digits=1))% |")
        end
        if config.max_precision !== nothing
            println(io, "| Precision | <$(config.max_precision*1000)nm | $(n_precision_pass) | $(n_raw - n_precision_pass) | $(round(100*n_precision_pass/n_raw, digits=1))% |")
        end
        if has_sigma && config.psf_sigma_mode_tolerance !== nothing && psf_sigma_mode > 0
            tol_pct = round(config.psf_sigma_mode_tolerance * 100, digits=0)
            println(io, "| PSF σ | mode ±$(tol_pct)% | $(n_psf_sigma_pass) | $(n_raw - n_psf_sigma_pass) | $(round(100*n_psf_sigma_pass/n_raw, digits=1))% |")
        elseif has_sigma_xy && config.psf_sigma_mode_tolerance !== nothing && psf_sigma_mode_x > 0 && psf_sigma_mode_y > 0
            tol_pct = round(config.psf_sigma_mode_tolerance * 100, digits=0)
            println(io, "| PSF σx,σy | mode ±$(tol_pct)% each | $(n_psf_sigma_pass) | $(n_raw - n_psf_sigma_pass) | $(round(100*n_psf_sigma_pass/n_raw, digits=1))% |")
        end
        if config.min_pvalue !== nothing
            println(io, "| P-value | >$(config.min_pvalue) | $(n_pvalue_pass) | $(n_raw - n_pvalue_pass) | $(round(100*n_pvalue_pass/n_raw, digits=1))% |")
        end

        println(io, "")
        println(io, "## Limiting Factor")
        min_pass = min(n_photon_pass, n_precision_pass, n_psf_sigma_pass, n_pvalue_pass)
        if min_pass == n_pvalue_pass && config.min_pvalue !== nothing
            println(io, "**P-value filter** is the most restrictive ($(round(100*n_pvalue_pass/n_raw, digits=1))% pass)")
        elseif min_pass == n_precision_pass && config.max_precision !== nothing
            println(io, "**Precision filter** is the most restrictive ($(round(100*n_precision_pass/n_raw, digits=1))% pass)")
        elseif min_pass == n_psf_sigma_pass && config.psf_sigma_mode_tolerance !== nothing
            println(io, "**PSF σ filter** is the most restrictive ($(round(100*n_psf_sigma_pass/n_raw, digits=1))% pass)")
        elseif min_pass == n_photon_pass && config.min_photons !== nothing
            println(io, "**Photon filter** is the most restrictive ($(round(100*n_photon_pass/n_raw, digits=1))% pass)")
        end

        println(io, "")
        println(io, "## Health Check")
        acc_rate = n_filtered / n_raw
        println(io, "- Acceptance rate: $(round(100*acc_rate, digits=1))% ",
                0.01 < acc_rate < 0.5 ? "✓" : (acc_rate < 0.01 ? "⚠ (too strict)" : "⚠ (too loose)"))
    end
end

"""Save drift correction figures - handles single and multi-dataset modes."""
function _save_drift_figures(drift_model, smld, config)
    DC = SMLMDriftCorrection
    n_datasets = length(drift_model.intra)
    n_frames = smld.n_frames
    frames = collect(1:n_frames)

    # Color palette for datasets
    colors = [:blue, :red, :green, :orange, :purple, :cyan, :magenta, :brown]

    if n_datasets == 1
        # Single dataset - simple 3-panel plot
        drift_x = [DC.applydrift(0.0, f, drift_model.intra[1].dm[1]) for f in frames]
        drift_y = [DC.applydrift(0.0, f, drift_model.intra[1].dm[2]) for f in frames]
        drift_x_nm = drift_x .* 1000
        drift_y_nm = drift_y .* 1000

        fig = Figure(size=(1400, 400))

        ax1 = Axis(fig[1, 1], xlabel="Frame", ylabel="X Drift (nm)", title="X Drift vs Frame")
        lines!(ax1, frames, drift_x_nm, color=:blue, linewidth=1.5)
        hlines!(ax1, [0], color=:gray, linestyle=:dash)

        ax2 = Axis(fig[1, 2], xlabel="Frame", ylabel="Y Drift (nm)", title="Y Drift vs Frame")
        lines!(ax2, frames, drift_y_nm, color=:red, linewidth=1.5)
        hlines!(ax2, [0], color=:gray, linestyle=:dash)

        ax3 = Axis(fig[1, 3], xlabel="X Drift (nm)", ylabel="Y Drift (nm)",
                   title="XY Drift Path", aspect=DataAspect())
        lines!(ax3, drift_x_nm, drift_y_nm, color=:black, linewidth=1.5)
        scatter!(ax3, [drift_x_nm[1]], [drift_y_nm[1]], color=:green, markersize=12, label="Start")
        scatter!(ax3, [drift_x_nm[end]], [drift_y_nm[end]], color=:red, markersize=12, label="End")
        axislegend(ax3, position=:lt)

        save(joinpath(config.outdir, "06_drift", "drift_trajectory.png"), fig)
    else
        # Multi-dataset - show per-dataset trajectories + inter-dataset shifts
        fig = Figure(size=(1600, 800))

        # Get frame ranges per dataset from emitters
        emitters = smld.emitters
        frame_ranges = Dict{Int, Tuple{Int,Int}}()
        for ds in 1:n_datasets
            ds_frames = [e.frame for e in emitters if e.dataset == ds]
            if !isempty(ds_frames)
                frame_ranges[ds] = (minimum(ds_frames), maximum(ds_frames))
            end
        end

        # Top row: X and Y drift by dataset
        ax1 = Axis(fig[1, 1], xlabel="Frame", ylabel="X Drift (nm)",
                   title="X Drift per Dataset (intra-dataset)")
        ax2 = Axis(fig[1, 2], xlabel="Frame", ylabel="Y Drift (nm)",
                   title="Y Drift per Dataset (intra-dataset)")

        for (i, ds) in enumerate(1:n_datasets)
            color = colors[mod1(i, length(colors))]
            # Only evaluate polynomial over frames that belong to this dataset
            ds_frames = haskey(frame_ranges, ds) ? collect(frame_ranges[ds][1]:frame_ranges[ds][2]) : frames
            drift_x = [DC.applydrift(0.0, f, drift_model.intra[ds].dm[1]) for f in ds_frames]
            drift_y = [DC.applydrift(0.0, f, drift_model.intra[ds].dm[2]) for f in ds_frames]
            lines!(ax1, ds_frames, drift_x .* 1000, color=color, linewidth=1.5, label="Dataset $ds")
            lines!(ax2, ds_frames, drift_y .* 1000, color=color, linewidth=1.5, label="Dataset $ds")
        end
        hlines!(ax1, [0], color=:gray, linestyle=:dash)
        hlines!(ax2, [0], color=:gray, linestyle=:dash)
        axislegend(ax1, position=:lt)
        axislegend(ax2, position=:lt)

        # Bottom left: XY paths per dataset
        ax3 = Axis(fig[2, 1], xlabel="X Drift (nm)", ylabel="Y Drift (nm)",
                   title="XY Drift Paths per Dataset", aspect=DataAspect())
        for (i, ds) in enumerate(1:n_datasets)
            color = colors[mod1(i, length(colors))]
            # Only evaluate polynomial over frames that belong to this dataset
            ds_frames = haskey(frame_ranges, ds) ? collect(frame_ranges[ds][1]:frame_ranges[ds][2]) : frames
            drift_x = [DC.applydrift(0.0, f, drift_model.intra[ds].dm[1]) for f in ds_frames] .* 1000
            drift_y = [DC.applydrift(0.0, f, drift_model.intra[ds].dm[2]) for f in ds_frames] .* 1000
            lines!(ax3, drift_x, drift_y, color=color, linewidth=1.5, label="Dataset $ds")
        end
        axislegend(ax3, position=:lt)

        # Bottom right: Inter-dataset shifts table
        ax4 = Axis(fig[2, 2], title="Inter-Dataset Alignment Shifts")
        hidedecorations!(ax4)
        hidespines!(ax4)

        shift_text = "Dataset → Shift X (nm) → Shift Y (nm)\n" * "─"^40 * "\n"
        for ds in 1:n_datasets
            shift_x = drift_model.inter[ds].dm[1] * 1000
            shift_y = drift_model.inter[ds].dm[2] * 1000
            shift_text *= "   $ds    →  $(round(shift_x, digits=1))  →  $(round(shift_y, digits=1))\n"
        end
        text!(ax4, 0.5, 0.5, text=shift_text, align=(:center, :center),
              fontsize=14)

        save(joinpath(config.outdir, "06_drift", "drift_trajectory.png"), fig)
    end
end

"""Write drift correction statistics markdown file - handles single and multi-dataset modes."""
function _write_drift_stats(drift_model, smld, config, elapsed_time)
    DC = SMLMDriftCorrection
    n_datasets = length(drift_model.intra)
    n_frames = smld.n_frames
    frames = collect(1:n_frames)

    filepath = joinpath(config.outdir, "06_drift", "drift_stats.md")
    open(filepath, "w") do io
        println(io, "# Drift Correction Statistics\n")
        println(io, "## Summary")
        println(io, "- **Time**: $(round(elapsed_time, digits=2))s")
        println(io, "- **Frames**: $n_frames")
        println(io, "- **Datasets**: $n_datasets")

        if n_datasets == 1
            # Single dataset - simple stats
            drift_x = [DC.applydrift(0.0, f, drift_model.intra[1].dm[1]) for f in frames]
            drift_y = [DC.applydrift(0.0, f, drift_model.intra[1].dm[2]) for f in frames]
            drift_x_nm = drift_x .* 1000
            drift_y_nm = drift_y .* 1000

            max_x = maximum(abs.(drift_x_nm))
            max_y = maximum(abs.(drift_y_nm))
            total_x = drift_x_nm[end] - drift_x_nm[1]
            total_y = drift_y_nm[end] - drift_y_nm[1]
            total_dist = sqrt(total_x^2 + total_y^2)

            println(io, "- **Max X drift**: $(round(max_x, digits=1)) nm")
            println(io, "- **Max Y drift**: $(round(max_y, digits=1)) nm")
            println(io, "- **Total displacement**: $(round(total_dist, digits=1)) nm")
            println(io, "")
            println(io, "## Parameters")
            println(io, "- Model: $(config.drift_model)")
            println(io, "- Degree: $(config.drift_degree)")
            println(io, "- Cost function: $(config.drift_cost_fun)")
            println(io, "")
            println(io, "## Health Check")
            println(io, "- Max drift: $(round(max(max_x, max_y), digits=0)) nm ",
                    max(max_x, max_y) < 500 ? "✓" : "⚠ (large drift)")
            println(io, "- Drift rate: $(round(total_dist / n_frames * 1000, digits=2)) nm/kframe")
        else
            # Multi-dataset - per-dataset stats + inter-dataset alignment
            println(io, "")
            println(io, "## Parameters")
            println(io, "- Model: $(config.drift_model)")
            println(io, "- Degree: $(config.drift_degree)")
            println(io, "- Cost function: $(config.drift_cost_fun)")
            println(io, "")
            println(io, "## Per-Dataset Intra-Dataset Drift")
            println(io, "")
            println(io, "| Dataset | Max X (nm) | Max Y (nm) | Total (nm) |")
            println(io, "|---------|------------|------------|------------|")

            # Get frame ranges per dataset from emitters (polynomials only valid within their dataset's frames)
            emitters = smld.emitters
            frame_ranges = Dict{Int, Tuple{Int,Int}}()
            for ds in 1:n_datasets
                ds_frames_list = [e.frame for e in emitters if e.dataset == ds]
                if !isempty(ds_frames_list)
                    frame_ranges[ds] = (minimum(ds_frames_list), maximum(ds_frames_list))
                end
            end

            overall_max = 0.0
            for ds in 1:n_datasets
                # Only evaluate polynomial over frames that belong to this dataset
                ds_frames = haskey(frame_ranges, ds) ? collect(frame_ranges[ds][1]:frame_ranges[ds][2]) : frames
                drift_x = [DC.applydrift(0.0, f, drift_model.intra[ds].dm[1]) for f in ds_frames]
                drift_y = [DC.applydrift(0.0, f, drift_model.intra[ds].dm[2]) for f in ds_frames]
                drift_x_nm = drift_x .* 1000
                drift_y_nm = drift_y .* 1000

                max_x = maximum(abs.(drift_x_nm))
                max_y = maximum(abs.(drift_y_nm))
                total_x = drift_x_nm[end] - drift_x_nm[1]
                total_y = drift_y_nm[end] - drift_y_nm[1]
                total_dist = sqrt(total_x^2 + total_y^2)

                overall_max = max(overall_max, max_x, max_y)
                println(io, "| $ds | $(round(max_x, digits=1)) | $(round(max_y, digits=1)) | $(round(total_dist, digits=1)) |")
            end

            println(io, "")
            println(io, "## Inter-Dataset Alignment Shifts")
            println(io, "")
            println(io, "| Dataset | Shift X (nm) | Shift Y (nm) |")
            println(io, "|---------|--------------|--------------|")

            max_shift = 0.0
            for ds in 1:n_datasets
                shift_x = drift_model.inter[ds].dm[1] * 1000
                shift_y = drift_model.inter[ds].dm[2] * 1000
                max_shift = max(max_shift, abs(shift_x), abs(shift_y))
                println(io, "| $ds | $(round(shift_x, digits=1)) | $(round(shift_y, digits=1)) |")
            end

            println(io, "")
            println(io, "## Health Check")
            println(io, "- Max intra-dataset drift: $(round(overall_max, digits=0)) nm ",
                    overall_max < 500 ? "✓" : "⚠ (large drift)")
            println(io, "- Max inter-dataset shift: $(round(max_shift, digits=0)) nm ",
                    max_shift < 1000 ? "✓" : "⚠ (large shift - check brightfield registration)")
        end
    end
end

"""Save isolated emitter filter figures - neighbor count histogram with triangle method visualization."""
function _save_isolated_figures(neighbor_counts, threshold, config)
    if isempty(neighbor_counts)
        return
    end

    auto_mode = config.isolated_min_neighbors == :auto
    title_suffix = auto_mode ? " (auto: triangle method)" : ""

    fig = Figure(size=(800, 500))

    ax = Axis(fig[1, 1],
        xlabel = "Number of Neighbors (within $(config.isolated_n_sigma)σ)",
        ylabel = "Count",
        title = "Neighbor Count Distribution (threshold = $threshold$title_suffix)"
    )

    # Histogram of neighbor counts
    max_count = min(maximum(neighbor_counts), 50)  # Cap at 50 for visualization
    hist!(ax, clamp.(neighbor_counts, 0, max_count), bins=range(0, max_count, length=max_count+1))

    # Mark threshold
    vlines!(ax, [threshold - 0.5], color=:red, linestyle=:dash, linewidth=2,
            label="Threshold ($threshold)")

    # Annotate rejection stats
    n_rejected = sum(neighbor_counts .< threshold)
    n_total = length(neighbor_counts)
    pct_rejected = round(100 * n_rejected / n_total, digits=1)

    method_str = auto_mode ? "Triangle method" : "Manual"
    text!(ax, 0.95, 0.95,
        text = "Method: $method_str\nThreshold: $threshold neighbors\nRejected: $n_rejected ($pct_rejected%)\nKept: $(n_total - n_rejected)",
        align = (:right, :top),
        space = :relative,
        fontsize = 11
    )

    axislegend(ax, position=:rt)

    save(joinpath(config.outdir, "07_isolated", "neighbor_histogram.png"), fig)
end

"""Write isolated emitter filter statistics markdown file."""
function _write_isolated_stats(n_before, n_after, neighbor_counts, threshold, config, elapsed_time)
    n_rejected = n_before - n_after
    auto_mode = config.isolated_min_neighbors == :auto

    filepath = joinpath(config.outdir, "07_isolated", "isolated_stats.md")
    open(filepath, "w") do io
        println(io, "# Isolated Emitter Filter Statistics\n")
        println(io, "## Summary")
        println(io, "- **Input**: $(n_before) localizations")
        println(io, "- **Output**: $(n_after) localizations")
        println(io, "- **Rejected**: $(n_rejected) ($(round(100*n_rejected/n_before, digits=1))%)")
        println(io, "- **Time**: $(round(elapsed_time, digits=2))s")
        println(io, "")
        println(io, "## Parameters")
        println(io, "- n_sigma: $(config.isolated_n_sigma) (neighbor if dist < n_sigma × σ_combined)")
        if auto_mode
            println(io, "- min_neighbors: **$threshold** (auto: triangle method)")
        else
            println(io, "- min_neighbors: $threshold (manual)")
        end
        println(io, "")
        println(io, "## Neighbor Count Distribution")
        if !isempty(neighbor_counts)
            p = quantile(neighbor_counts, [0.0, 0.05, 0.25, 0.50, 0.75, 0.95, 1.0])
            println(io, "- Min: $(Int(p[1]))")
            println(io, "- 5th percentile: $(Int(p[2]))")
            println(io, "- 25th percentile: $(Int(p[3]))")
            println(io, "- Median: $(Int(p[4]))")
            println(io, "- 75th percentile: $(Int(p[5]))")
            println(io, "- 95th percentile: $(Int(p[6]))")
            println(io, "- Max: $(Int(p[7]))")
        end
        println(io, "")
        println(io, "## Health Check")
        rejection_rate = n_rejected / n_before
        println(io, "- Rejection rate: $(round(100*rejection_rate, digits=1))% ",
                rejection_rate < 0.2 ? "✓" : "⚠ (high rejection)")
    end
end

"""Save frame connection figures - histogram of locs per track."""
function _save_frameconnect_figures(smld_connected, config)
    emitters = smld_connected.emitters
    if isempty(emitters)
        return
    end

    # Count locs per track
    track_ids = [e.track_id for e in emitters]
    counts = countmap(track_ids)
    locs_per_track = collect(values(counts))

    fig = Figure(size=(800, 500))

    ax = Axis(fig[1, 1],
        xlabel = "Localizations per Track",
        ylabel = "Count (tracks)",
        title = "Frame Connection: Track Size Distribution"
    )

    # Histogram of track sizes (cap at 50 for viz)
    max_locs = min(maximum(locs_per_track), 50)
    hist!(ax, clamp.(locs_per_track, 1, max_locs),
          bins=range(0.5, max_locs + 0.5, length=max_locs + 1),
          color=:steelblue)

    # Stats annotation
    n_tracks = length(locs_per_track)
    n_singles = sum(locs_per_track .== 1)
    pct_singles = round(100 * n_singles / n_tracks, digits=1)
    med_size = median(locs_per_track)
    max_size = maximum(locs_per_track)

    text!(ax, 0.95, 0.95,
        text = "Tracks: $n_tracks\nSingles: $n_singles ($pct_singles%)\nMedian: $(Int(med_size))\nMax: $max_size",
        align = (:right, :top),
        space = :relative,
        fontsize = 12
    )

    save(joinpath(config.outdir, "04_frameconnect", "track_histogram.png"), fig)
end

"""
Analyze frame-to-frame drift from linked localizations.

For multi-dataset analysis, drift is calculated per dataset since frame numbers
reset between datasets and datasets may be acquired at different times.

Returns NamedTuple with:
- frame_shifts: Dict{Int, Vector} mapping dataset_id => shifts per frame transition
- chi2_values: Vector of χ² values for each pair (should follow χ²(2) if uncertainties correct)
- mean_chi2: Mean χ² (expected = 2 for correct uncertainties)
- n_pairs_total: Total number of frame-to-frame pairs analyzed
- n_datasets: Number of datasets with tracked emitters
"""
function _analyze_frameconnect_drift(smld_connected)
    emitters = smld_connected.emitters
    n_datasets = smld_connected.n_datasets

    # Group emitters by track_id
    track_dict = Dict{Int, Vector{eltype(emitters)}}()
    for e in emitters
        if e.track_id > 0
            if !haskey(track_dict, e.track_id)
                track_dict[e.track_id] = eltype(emitters)[]
            end
            push!(track_dict[e.track_id], e)
        end
    end

    # Collect all frame-to-frame pairs, keyed by (dataset_id, frame)
    # Structure: (dataset, frame) => Vector of (Δx, Δy, var_x, var_y)
    frame_pairs = Dict{Tuple{Int, Int}, Vector{NTuple{4, Float64}}}()
    chi2_values = Float64[]

    for (track_id, track_emitters) in track_dict
        # Sort by dataset first, then frame
        sort!(track_emitters, by = e -> (e.dataset, e.frame))

        # Find consecutive frame pairs within same dataset
        for i in 1:(length(track_emitters) - 1)
            e1, e2 = track_emitters[i], track_emitters[i + 1]

            # Only consecutive frames within same dataset
            if e2.dataset == e1.dataset && e2.frame == e1.frame + 1
                Δx = Float64(e2.x - e1.x)
                Δy = Float64(e2.y - e1.y)
                var_x = Float64(e1.σ_x^2 + e2.σ_x^2)
                var_y = Float64(e1.σ_y^2 + e2.σ_y^2)

                key = (e1.dataset, e1.frame)
                if !haskey(frame_pairs, key)
                    frame_pairs[key] = NTuple{4, Float64}[]
                end
                push!(frame_pairs[key], (Δx, Δy, var_x, var_y))

                # Chi-squared for this pair
                if var_x > 0 && var_y > 0
                    χ2 = Δx^2 / var_x + Δy^2 / var_y
                    push!(chi2_values, χ2)
                end
            end
        end
    end

    # Calculate weighted mean shift for each (dataset, frame) transition
    # Output: Dict{dataset_id => Vector of (frame, Δx, Δy, σ_Δx, σ_Δy, n_pairs)}
    frame_shifts = Dict{Int, Vector{Tuple{Int, Float64, Float64, Float64, Float64, Int}}}()

    for (dataset_id, frame) in sort(collect(keys(frame_pairs)))
        pairs = frame_pairs[(dataset_id, frame)]
        n_pairs = length(pairs)

        if n_pairs > 0
            # Weighted average: weight = 1/variance
            sum_wx, sum_wy = 0.0, 0.0
            sum_w_x, sum_w_y = 0.0, 0.0

            for (Δx, Δy, var_x, var_y) in pairs
                if var_x > 0 && var_y > 0
                    w_x, w_y = 1.0 / var_x, 1.0 / var_y
                    sum_wx += w_x * Δx
                    sum_wy += w_y * Δy
                    sum_w_x += w_x
                    sum_w_y += w_y
                end
            end

            if sum_w_x > 0 && sum_w_y > 0
                Δx_mean = sum_wx / sum_w_x
                Δy_mean = sum_wy / sum_w_y
                σ_Δx = 1.0 / sqrt(sum_w_x)
                σ_Δy = 1.0 / sqrt(sum_w_y)

                if !haskey(frame_shifts, dataset_id)
                    frame_shifts[dataset_id] = Tuple{Int, Float64, Float64, Float64, Float64, Int}[]
                end
                push!(frame_shifts[dataset_id], (frame, Δx_mean, Δy_mean, σ_Δx, σ_Δy, n_pairs))
            end
        end
    end

    mean_chi2 = isempty(chi2_values) ? NaN : mean(chi2_values)

    # Collect all pair data for uncertainty calibration analysis
    # Each entry: (Δx², Δy², var_x, var_y)
    pair_data = NTuple{4, Float64}[]
    for pairs in values(frame_pairs)
        for (Δx, Δy, var_x, var_y) in pairs
            if var_x > 0 && var_y > 0
                push!(pair_data, (Δx^2, Δy^2, var_x, var_y))
            end
        end
    end

    # Fit uncertainty calibration model: observed_var = A + B * reported_var
    # Using simple linear regression on binned data
    calibration = _fit_uncertainty_calibration(pair_data)

    # Estimate stage motion from variance of mean shifts
    # If stage vibrates, all emitters move together → mean shift variance > expected from uncertainties
    # Collect all mean shifts across datasets
    all_Δx_mean = Float64[]
    all_Δy_mean = Float64[]
    all_σ_Δx = Float64[]
    all_σ_Δy = Float64[]
    for (dataset_id, shifts) in frame_shifts
        for (frame, Δx, Δy, σx, σy, n) in shifts
            push!(all_Δx_mean, Δx * 1000)  # Convert to nm
            push!(all_Δy_mean, Δy * 1000)
            push!(all_σ_Δx, σx * 1000)
            push!(all_σ_Δy, σy * 1000)
        end
    end

    # Observed variance of mean shifts
    var_Δx_observed = length(all_Δx_mean) > 1 ? var(all_Δx_mean) : 0.0
    var_Δy_observed = length(all_Δy_mean) > 1 ? var(all_Δy_mean) : 0.0

    # Expected variance from localization uncertainty alone
    var_Δx_expected = length(all_σ_Δx) > 0 ? mean(all_σ_Δx.^2) : 0.0
    var_Δy_expected = length(all_σ_Δy) > 0 ? mean(all_σ_Δy.^2) : 0.0

    # Excess variance = stage motion
    excess_var_x = max(0.0, var_Δx_observed - var_Δx_expected)
    excess_var_y = max(0.0, var_Δy_observed - var_Δy_expected)
    σ_stage_x = sqrt(excess_var_x)
    σ_stage_y = sqrt(excess_var_y)
    σ_stage = sqrt((excess_var_x + excess_var_y) / 2)  # Combined estimate

    stage_motion = (
        σ_stage_x = σ_stage_x,
        σ_stage_y = σ_stage_y,
        σ_stage = σ_stage,
        var_observed_x = var_Δx_observed,
        var_observed_y = var_Δy_observed,
        var_expected_x = var_Δx_expected,
        var_expected_y = var_Δy_expected,
        n_frames = length(all_Δx_mean)
    )

    return (
        frame_shifts = frame_shifts,
        chi2_values = chi2_values,
        mean_chi2 = mean_chi2,
        n_pairs_total = length(chi2_values),
        n_datasets = length(frame_shifts),
        pair_data = pair_data,
        calibration = calibration,
        stage_motion = stage_motion
    )
end

"""
Fit uncertainty calibration model: observed_variance = A + B * CRLB_variance

Bins pairs by reported variance, computes mean observed variance per bin,
fits linear model to determine if discrepancy is additive (A), multiplicative (B), or both.

Returns NamedTuple with:
- A: additive term (nm²) - represents motion/vibration variance
- B: multiplicative factor - CRLB scale correction
- A_σ, B_σ: uncertainties on A and B
- r_squared: goodness of fit
- bin_centers, bin_observed, bin_expected: binned data for plotting
"""
function _fit_uncertainty_calibration(pair_data)
    if length(pair_data) < 100
        return (A = NaN, B = NaN, A_σ = NaN, B_σ = NaN, r_squared = NaN,
                bin_centers = Float64[], bin_observed = Float64[], bin_expected = Float64[], n_per_bin = Int[])
    end

    # Combine x and y data: use average of x and y for each pair
    reported_var = [(p[3] + p[4]) / 2 for p in pair_data]  # Average of var_x and var_y
    observed_var = [(p[1] + p[2]) / 2 for p in pair_data]   # Average of Δx² and Δy²

    # Convert to nm² for interpretability
    reported_var_nm2 = reported_var .* 1e6
    observed_var_nm2 = observed_var .* 1e6

    # Bin by reported variance (use quantiles for equal-count bins)
    n_bins = 20
    sorted_idx = sortperm(reported_var_nm2)
    bin_size = length(sorted_idx) ÷ n_bins

    bin_centers = Float64[]
    bin_observed = Float64[]
    bin_expected = Float64[]
    n_per_bin = Int[]

    for i in 1:n_bins
        start_idx = (i - 1) * bin_size + 1
        end_idx = i == n_bins ? length(sorted_idx) : i * bin_size
        bin_indices = sorted_idx[start_idx:end_idx]

        push!(bin_centers, mean(reported_var_nm2[bin_indices]))
        push!(bin_observed, mean(observed_var_nm2[bin_indices]))
        push!(bin_expected, mean(reported_var_nm2[bin_indices]))  # Expected if perfectly calibrated
        push!(n_per_bin, length(bin_indices))
    end

    # Linear regression: observed = A + B * reported
    # Using normal equations: [A; B] = (X'X)^-1 X'y
    X = hcat(ones(n_bins), bin_centers)
    y = bin_observed
    XtX = X' * X
    Xty = X' * y

    if det(XtX) ≈ 0
        return (A = NaN, B = NaN, A_σ = NaN, B_σ = NaN, r_squared = NaN,
                bin_centers = bin_centers, bin_observed = bin_observed, bin_expected = bin_expected, n_per_bin = n_per_bin)
    end

    coeffs = XtX \ Xty
    A, B = coeffs[1], coeffs[2]

    # Compute R² and parameter uncertainties
    y_pred = X * coeffs
    ss_res = sum((y .- y_pred).^2)
    ss_tot = sum((y .- mean(y)).^2)
    r_squared = 1 - ss_res / ss_tot

    # Standard errors
    mse = ss_res / (n_bins - 2)
    var_coeffs = mse * inv(XtX)
    A_σ = sqrt(var_coeffs[1, 1])
    B_σ = sqrt(var_coeffs[2, 2])

    return (A = A, B = B, A_σ = A_σ, B_σ = B_σ, r_squared = r_squared,
            bin_centers = bin_centers, bin_observed = bin_observed, bin_expected = bin_expected, n_per_bin = n_per_bin)
end

"""Save frame connection drift analysis figures.

For multi-dataset analysis, plots separate drift trajectories per dataset with
distinct colors. Cumulative drift is calculated separately per dataset since
datasets may be acquired at different times.
"""
function _save_frameconnect_drift_figures(drift_analysis, config)
    frame_shifts = drift_analysis.frame_shifts  # Dict{dataset_id => Vector of shifts}
    chi2_values = drift_analysis.chi2_values
    n_datasets = drift_analysis.n_datasets

    if isempty(frame_shifts)
        return
    end

    fig = Figure(size=(1600, 800))

    # Color palette for datasets
    colors_x = [:blue, :teal, :purple, :navy, :cyan]
    colors_y = [:red, :orange, :magenta, :brown, :coral]

    # Top left: Frame-to-frame X shift with error bars
    ax1 = Axis(fig[1, 1], xlabel="Frame", ylabel="ΔX (nm)",
               title="Frame-to-Frame X Shift")
    hlines!(ax1, [0], color=:gray, linestyle=:dash)

    # Top right: Frame-to-frame Y shift with error bars
    ax2 = Axis(fig[1, 2], xlabel="Frame", ylabel="ΔY (nm)",
               title="Frame-to-Frame Y Shift")
    hlines!(ax2, [0], color=:gray, linestyle=:dash)

    # Bottom left: Cumulative drift with error envelope
    title_suffix = n_datasets > 1 ? " (per dataset)" : ""
    ax3 = Axis(fig[2, 1], xlabel="Frame", ylabel="Cumulative Drift (nm)",
               title="Cumulative Drift from Linked Emitters$title_suffix")

    # Plot each dataset
    for (i, (dataset_id, shifts)) in enumerate(sort(collect(frame_shifts)))
        if isempty(shifts)
            continue
        end

        frames = [s[1] for s in shifts]
        Δx = [s[2] for s in shifts] .* 1000  # Convert to nm
        Δy = [s[3] for s in shifts] .* 1000
        σ_Δx = [s[4] for s in shifts] .* 1000
        σ_Δy = [s[5] for s in shifts] .* 1000

        color_x = colors_x[mod1(i, length(colors_x))]
        color_y = colors_y[mod1(i, length(colors_y))]

        # Frame-to-frame shifts
        alpha = n_datasets > 1 ? 0.15 : 0.3
        band!(ax1, frames, Δx .- σ_Δx, Δx .+ σ_Δx, color=(color_x, alpha))
        lines!(ax1, frames, Δx, color=color_x, linewidth=1)

        band!(ax2, frames, Δy .- σ_Δy, Δy .+ σ_Δy, color=(color_y, alpha))
        lines!(ax2, frames, Δy, color=color_y, linewidth=1)

        # Cumulative drift per dataset
        cum_x = cumsum(Δx)
        cum_y = cumsum(Δy)
        cum_σx = sqrt.(cumsum(σ_Δx .^ 2))
        cum_σy = sqrt.(cumsum(σ_Δy .^ 2))

        if n_datasets > 1
            label_x = "X (ds $dataset_id)"
            label_y = "Y (ds $dataset_id)"
        else
            label_x = "X"
            label_y = "Y"
        end

        band!(ax3, frames, cum_x .- cum_σx, cum_x .+ cum_σx, color=(color_x, 0.15))
        band!(ax3, frames, cum_y .- cum_σy, cum_y .+ cum_σy, color=(color_y, 0.15))
        lines!(ax3, frames, cum_x, color=color_x, linewidth=1.5, label=label_x)
        lines!(ax3, frames, cum_y, color=color_y, linewidth=1.5, label=label_y)
    end

    axislegend(ax3, position=:lt)

    # Bottom right: Chi-squared histogram (pooled across all datasets)
    ax4 = Axis(fig[2, 2], xlabel="χ²", ylabel="Count",
               title="Uncertainty Validation (expected χ²=2)")
    if !isempty(chi2_values)
        # Clip extreme values for visualization
        chi2_clipped = clamp.(chi2_values, 0, 20)
        hist!(ax4, chi2_clipped, bins=range(0, 20, length=41), color=:steelblue)
        vlines!(ax4, [2.0], color=:red, linewidth=2, linestyle=:dash, label="Expected")
        vlines!(ax4, [mean(chi2_values)], color=:orange, linewidth=2, label="Observed mean")
        axislegend(ax4, position=:rt)

        # Annotation
        mean_chi2 = drift_analysis.mean_chi2
        text!(ax4, 0.95, 0.95,
            text = "Mean χ² = $(round(mean_chi2, digits=2))\nExpected = 2.0\nn = $(length(chi2_values))",
            align = (:right, :top),
            space = :relative,
            fontsize = 11
        )
    end

    save(joinpath(config.outdir, "04_frameconnect", "drift_from_tracks.png"), fig)
end

"""Save uncertainty calibration figure: observed vs reported variance."""
function _save_uncertainty_calibration_figure(drift_analysis, config)
    cal = drift_analysis.calibration

    if isnan(cal.A) || isempty(cal.bin_centers)
        return
    end

    fig = Figure(size=(800, 600))

    ax = Axis(fig[1, 1],
        xlabel = "Reported Variance σ²_CRLB (nm²)",
        ylabel = "Observed Variance ⟨Δ²⟩ (nm²)",
        title = "Uncertainty Calibration: Observed = A + B × Reported"
    )

    # Plot binned data
    scatter!(ax, cal.bin_centers, cal.bin_observed, markersize=8, color=:steelblue, label="Binned data")

    # Plot 1:1 line (perfect calibration)
    x_range = range(minimum(cal.bin_centers), maximum(cal.bin_centers), length=100)
    lines!(ax, collect(x_range), collect(x_range), color=:gray, linestyle=:dash, linewidth=2, label="Perfect (1:1)")

    # Plot fitted line
    y_fit = cal.A .+ cal.B .* x_range
    lines!(ax, collect(x_range), collect(y_fit), color=:red, linewidth=2, label="Fit: A + B×σ²")

    axislegend(ax, position=:lt)

    # Annotation with fit results (show precision, not variance)
    σ_motion = cal.A > 0 ? sqrt(cal.A) : 0.0
    σ_motion_err = cal.A > 0 ? cal.A_σ / (2 * σ_motion) : 0.0
    k_scale = sqrt(cal.B)
    k_scale_err = cal.B_σ / (2 * k_scale)

    interpretation = if σ_motion < 2.5 && abs(k_scale - 1) < 0.15
        "Well calibrated"
    elseif σ_motion < 2.5 && k_scale > 1.15
        "CRLB underestimates by $(round(k_scale, digits=2))×"
    elseif σ_motion > 2.5 && abs(k_scale - 1) < 0.15
        "Motion noise: $(round(σ_motion, digits=1)) nm"
    else
        "Both: motion $(round(σ_motion, digits=1)) nm, CRLB ×$(round(k_scale, digits=2))"
    end

    text!(ax, 0.95, 0.05,
        text = "σ_motion = $(round(σ_motion, digits=1)) ± $(round(σ_motion_err, digits=1)) nm\nk = $(round(k_scale, digits=2)) ± $(round(k_scale_err, digits=2))\nR² = $(round(cal.r_squared, digits=3))\n\n$interpretation",
        align = (:right, :bottom),
        space = :relative,
        fontsize = 11
    )

    save(joinpath(config.outdir, "04_frameconnect", "uncertainty_calibration.png"), fig)
end

"""Write frame connection statistics markdown file."""
function _write_frameconnect_stats(n_before, n_after, smld_connected, fc_params, config, elapsed_time; drift_analysis=nothing)
    # Count locs per track
    track_ids = [e.track_id for e in smld_connected.emitters]
    counts = countmap(track_ids)
    locs_per_track = collect(values(counts))

    compression = n_before / n_after
    n_tracks = n_after
    n_singles = sum(locs_per_track .== 1)

    filepath = joinpath(config.outdir, "04_frameconnect", "frameconnect_stats.md")
    open(filepath, "w") do io
        println(io, "# Frame Connection Statistics\n")
        println(io, "## Summary")
        println(io, "- **Input**: $(n_before) localizations")
        println(io, "- **Output**: $(n_after) tracks")
        println(io, "- **Compression**: $(round(compression, digits=2))x")
        println(io, "- **Time**: $(round(elapsed_time, digits=2))s")
        println(io, "")
        println(io, "## Track Statistics")
        println(io, "- Total tracks: $(n_tracks)")
        println(io, "- Single-loc tracks: $(n_singles) ($(round(100*n_singles/n_tracks, digits=1))%)")
        println(io, "- Mean locs/track: $(round(mean(locs_per_track), digits=2))")
        println(io, "- Median locs/track: $(Int(median(locs_per_track)))")
        println(io, "- Max locs/track: $(maximum(locs_per_track))")
        println(io, "")
        println(io, "## Parameters")
        println(io, "- maxframegap: $(config.fc_maxframegap)")
        println(io, "- nsigmadev: $(config.fc_nsigmadev)")
        println(io, "- nnearestclusters: $(config.fc_nnearestclusters)")
        println(io, "- nmaxnn: $(config.fc_nmaxnn)")
        println(io, "")

        # Drift analysis from linked emitters
        if drift_analysis !== nothing && drift_analysis.n_pairs_total > 0
            println(io, "## Drift Analysis from Linked Emitters")
            println(io, "")
            println(io, "Frame-to-frame shifts estimated from $(drift_analysis.n_pairs_total) consecutive-frame pairs.")
            if drift_analysis.n_datasets > 1
                println(io, "Drift calculated separately for $(drift_analysis.n_datasets) datasets.")
            end
            println(io, "")

            # Total drift per dataset
            frame_shifts = drift_analysis.frame_shifts
            if !isempty(frame_shifts)
                println(io, "| Dataset | Total ΔX (nm) | Total ΔY (nm) |")
                println(io, "|---------|---------------|---------------|")

                for (dataset_id, shifts) in sort(collect(frame_shifts))
                    if !isempty(shifts)
                        Δx_total = sum(s[2] for s in shifts) * 1000  # nm
                        Δy_total = sum(s[3] for s in shifts) * 1000
                        σ_Δx_total = sqrt(sum(s[4]^2 for s in shifts)) * 1000
                        σ_Δy_total = sqrt(sum(s[5]^2 for s in shifts)) * 1000

                        println(io, "| $dataset_id | $(round(Δx_total, digits=1)) ± $(round(σ_Δx_total, digits=1)) | $(round(Δy_total, digits=1)) ± $(round(σ_Δy_total, digits=1)) |")
                    end
                end
                println(io, "")
            end

            # Uncertainty validation
            println(io, "## Uncertainty Validation (χ² test)")
            println(io, "")
            println(io, "If localization uncertainties are correct, χ² should follow χ²(2) with mean=2.")
            println(io, "")
            mean_chi2 = drift_analysis.mean_chi2
            println(io, "- **Observed mean χ²**: $(round(mean_chi2, digits=2))")
            println(io, "- **Expected mean χ²**: 2.0")
            println(io, "- **Ratio**: $(round(mean_chi2 / 2.0, digits=2))")
            println(io, "")

            if mean_chi2 > 3.0
                println(io, "⚠ **Uncertainties underestimated** - observed motion exceeds expected from uncertainties")
            elseif mean_chi2 < 1.0
                println(io, "⚠ **Uncertainties overestimated** - emitters more stable than uncertainties predict")
            else
                println(io, "✓ Uncertainties appear well-calibrated")
            end
            println(io, "")

            # Uncertainty calibration model
            cal = drift_analysis.calibration
            if !isnan(cal.A)
                println(io, "## Uncertainty Calibration Model")
                println(io, "")
                println(io, "Fit model: `σ²_observed = σ²_motion + k² × σ²_CRLB`")
                println(io, "")

                # Convert variance to precision with error propagation
                # σ = √A, so σ_err = A_err / (2√A)
                σ_motion = cal.A > 0 ? sqrt(cal.A) : 0.0
                σ_motion_err = cal.A > 0 ? cal.A_σ / (2 * σ_motion) : 0.0
                k_scale = sqrt(cal.B)
                k_scale_err = cal.B_σ / (2 * k_scale)

                println(io, "| Parameter | Value | Description |")
                println(io, "|-----------|-------|-------------|")
                println(io, "| σ_motion | $(round(σ_motion, digits=1)) ± $(round(σ_motion_err, digits=1)) nm | Frame-to-frame motion/vibration |")
                println(io, "| k (CRLB scale) | $(round(k_scale, digits=2)) ± $(round(k_scale_err, digits=2)) | Multiply σ_CRLB by this |")
                println(io, "| R² | $(round(cal.r_squared, digits=3)) | Fit quality |")
                println(io, "")

                # Interpretation
                if σ_motion < 2.5 && abs(k_scale - 1) < 0.15
                    println(io, "✓ **Well calibrated** - CRLB matches observed variance")
                elseif σ_motion < 2.5 && k_scale > 1.15
                    println(io, "⚠ **CRLB underestimates** - multiply σ_CRLB by $(round(k_scale, digits=2))")
                elseif σ_motion > 2.5 && abs(k_scale - 1) < 0.15
                    println(io, "⚠ **Motion noise** - $(round(σ_motion, digits=1)) nm per frame (stage vibration, thermal)")
                else
                    println(io, "⚠ **Both effects** - motion $(round(σ_motion, digits=1)) nm + CRLB scale $(round(k_scale, digits=2))×")
                end
                println(io, "")
            end

            # Stage motion analysis - compare to calibration model
            sm = drift_analysis.stage_motion
            if sm.n_frames > 10
                println(io, "## Stage Motion Analysis")
                println(io, "")
                println(io, "Variance of frame-to-frame mean shifts vs expected from localization uncertainty:")
                println(io, "")
                println(io, "| Direction | Var(observed) | Var(expected) | Excess → σ_stage |")
                println(io, "|-----------|---------------|---------------|------------------|")
                println(io, "| X | $(round(sm.var_observed_x, digits=1)) nm² | $(round(sm.var_expected_x, digits=1)) nm² | $(round(sm.σ_stage_x, digits=1)) nm |")
                println(io, "| Y | $(round(sm.var_observed_y, digits=1)) nm² | $(round(sm.var_expected_y, digits=1)) nm² | $(round(sm.σ_stage_y, digits=1)) nm |")
                println(io, "| **Combined** | | | **$(round(sm.σ_stage, digits=1)) nm** |")
                println(io, "")

                # Compare to calibration model
                if !isnan(cal.A)
                    println(io, "**Comparison:**")
                    println(io, "- From calibration model: σ_motion = $(round(σ_motion, digits=1)) nm")
                    println(io, "- From mean shift variance: σ_stage = $(round(sm.σ_stage, digits=1)) nm")
                    println(io, "")
                    if abs(σ_motion - sm.σ_stage) < 1.5
                        println(io, "✓ **Consistent** - motion is coherent stage vibration/drift")
                    elseif σ_motion > sm.σ_stage + 1.5
                        println(io, "⚠ **Calibration σ > stage σ** - some motion is per-emitter (sample drift, diffusion)")
                    else
                        println(io, "⚠ **Stage σ > calibration σ** - unexpected, check data quality")
                    end
                    println(io, "")
                end
            end
        end

        println(io, "## Health Check")
        single_rate = n_singles / n_tracks
        println(io, "- Single track rate: $(round(100*single_rate, digits=1))% ",
                single_rate > 0.5 ? "⚠ (many isolated locs)" : "✓")
        println(io, "- Compression ratio: $(round(compression, digits=1))x ",
                compression < 1.5 ? "⚠ (low compression)" : "✓")
    end
end

"""Write render statistics markdown file."""
function _write_render_stats(smld, config, elapsed_time)
    emitters = smld.emitters
    n = length(emitters)

    x = [e.x for e in emitters]
    y = [e.y for e in emitters]

    x_range = maximum(x) - minimum(x)
    y_range = maximum(y) - minimum(y)
    area_um2 = x_range * y_range
    density = n / area_um2

    # Estimate Nyquist resolution
    σ_x = [e.σ_x for e in emitters]
    mean_precision = mean(σ_x) * 1000  # nm
    nyquist_resolution = 2 * mean_precision  # nm

    filepath = joinpath(config.outdir, "08_superres", "render_stats.md")
    open(filepath, "w") do io
        println(io, "# Render Statistics\n")
        println(io, "## Summary")
        println(io, "- **Localizations rendered**: $(n)")
        println(io, "- **Time**: $(round(elapsed_time, digits=2))s")
        println(io, "- **Field of view**: $(round(x_range, digits=1)) × $(round(y_range, digits=1)) μm")
        println(io, "- **Area**: $(round(area_um2, digits=1)) μm²")
        println(io, "- **Localization density**: $(round(density, digits=1)) /μm²")
        println(io, "")
        println(io, "## Renders Generated\n")
        println(io, "| Render | Zoom | Pixel Size | Purpose |")
        println(io, "|--------|------|------------|---------|")
        if config.render_gaussian
            ps = round(1000 / config.render_gaussian_zoom * 0.078, digits=1)
            println(io, "| Gaussian (inferno) | $(config.render_gaussian_zoom)x | $(ps) nm | Primary super-res |")
        end
        if config.render_histogram
            ps = round(1000 / config.render_histogram_zoom * 0.078, digits=1)
            println(io, "| Histogram (time) | $(config.render_histogram_zoom)x | $(ps) nm | Temporal coverage |")
        end
        if config.render_circles
            ps = round(1000 / config.render_circles_zoom * 0.078, digits=1)
            println(io, "| Circles (time) | $(config.render_circles_zoom)x | $(ps) nm | Individual locs |")
        end
        println(io, "")
        println(io, "## Resolution Estimate")
        println(io, "- Mean precision: $(round(mean_precision, digits=1)) nm")
        println(io, "- Nyquist resolution (2×precision): $(round(nyquist_resolution, digits=1)) nm")
        println(io, "")
        println(io, "## Health Check")
        println(io, "- Localization density: $(round(density, digits=1)) /μm² ", density > 10 ? "✓" : "⚠ (sparse)")
    end
end

"""Write summary markdown file with health check."""
function _write_summary(roi_batch, smld_raw, smld_filtered, data, config, timings)
    nframes = size(data, 3)
    n_rois = length(roi_batch)
    n_raw = length(smld_raw.emitters)
    n_filtered = length(smld_filtered.emitters)

    emitters_raw = smld_raw.emitters
    photons = [e.photons for e in emitters_raw]
    σ_x = [e.σ_x for e in emitters_raw]
    pvalue = [e.pvalue for e in emitters_raw]

    has_sigma = hasfield(typeof(emitters_raw[1]), :σ)
    has_sigma_xy = hasfield(typeof(emitters_raw[1]), :σx) && hasfield(typeof(emitters_raw[1]), :σy)
    psf_sigma = has_sigma ? median([e.σ for e in emitters_raw]) * 1000 : nothing
    psf_sigma_x = has_sigma_xy ? median([e.σx for e in emitters_raw]) * 1000 : nothing
    psf_sigma_y = has_sigma_xy ? median([e.σy for e in emitters_raw]) * 1000 : nothing

    rois_per_frame = n_rois / nframes
    pval_pass_rate = sum(pvalue .> 0.001) / n_raw
    acceptance_rate = n_filtered / n_raw

    filepath = joinpath(config.outdir, "summary.md")
    open(filepath, "w") do io
        println(io, "# SMLM Analysis Summary\n")
        println(io, "## Quick Health Check\n")
        println(io, "```")
        println(io, "Detection:    $(round(rois_per_frame, digits=1)) ROIs/frame     ", 50 < rois_per_frame < 500 ? "✓" : "⚠")
        println(io, "Fitting:      $(round(100*pval_pass_rate, digits=1))% pvalue>0.001  ", pval_pass_rate > 0.05 ? "✓" : "⚠")
        if psf_sigma !== nothing
            println(io, "PSF sigma:    $(round(psf_sigma, digits=0)) nm            ", 80 < psf_sigma < 250 ? "✓" : "⚠")
        elseif psf_sigma_x !== nothing && psf_sigma_y !== nothing
            avg_sigma = (psf_sigma_x + psf_sigma_y) / 2
            println(io, "PSF σx,σy:    $(round(psf_sigma_x, digits=0)),$(round(psf_sigma_y, digits=0)) nm     ", 80 < avg_sigma < 250 ? "✓" : "⚠")
        end
        println(io, "Precision:    $(round(median(σ_x)*1000, digits=1)) nm median      ", median(σ_x)*1000 < 20 ? "✓" : "⚠")
        println(io, "Photons:      $(round(median(photons), digits=0)) median      ", 1000 < median(photons) < 100000 ? "✓" : "⚠")
        println(io, "Filtering:    $(round(100*acceptance_rate, digits=1))% accepted      ", 0.01 < acceptance_rate < 0.5 ? "✓" : "⚠")
        println(io, "```")

        # Overall status
        psf_sigma_ok = if psf_sigma !== nothing
            80 < psf_sigma < 250
        elseif psf_sigma_x !== nothing && psf_sigma_y !== nothing
            avg = (psf_sigma_x + psf_sigma_y) / 2
            80 < avg < 250
        else
            true  # No PSF sigma to check (fixed model)
        end

        all_ok = (50 < rois_per_frame < 500) &&
                 (pval_pass_rate > 0.05) &&
                 psf_sigma_ok &&
                 (median(σ_x)*1000 < 20) &&
                 (1000 < median(photons) < 100000) &&
                 (0.01 < acceptance_rate < 0.5)
        println(io, "\n**STATUS: $(all_ok ? "HEALTHY" : "NEEDS ATTENTION")**\n")

        println(io, "## Pipeline Results\n")
        println(io, "| Step | Count | Time |")
        println(io, "|------|-------|------|")
        println(io, "| Detection | $(n_rois) ROIs | $(round(get(timings, "detection", 0.0), digits=1))s |")
        println(io, "| Fitting | $(n_raw) fits | $(round(get(timings, "fitting", 0.0), digits=1))s |")
        println(io, "| Filtering | $(n_filtered) kept ($(round(100*acceptance_rate, digits=1))%) | $(round(get(timings, "filtering", 0.0), digits=1))s |")
        if config.frameconnect
            println(io, "| Frame Connection | - | $(round(get(timings, "frameconnect", 0.0), digits=1))s |")
        end
        if config.calibrate_uncertainties && config.frameconnect
            println(io, "| Uncertainty Calibration | - | $(round(get(timings, "calibration", 0.0), digits=1))s |")
        end
        if config.drift
            println(io, "| Drift Correction | - | $(round(get(timings, "drift", 0.0), digits=1))s |")
        end
        if config.filter_isolated
            println(io, "| Isolated Filter | - | $(round(get(timings, "isolated_filter", 0.0), digits=1))s |")
        end
        if config.render
            println(io, "| Rendering | - | $(round(get(timings, "rendering", 0.0), digits=1))s |")
        end
        println(io, "| **Total** | **$(n_filtered) localizations** | **$(round(sum(values(timings)), digits=1))s** |")

        println(io, "\n## Configuration\n")
        println(io, "- Frames: $(nframes)")
        println(io, "- Fit model: $(config.fit_model)")
        # Build filter description
        filter_parts = String[]
        if config.min_photons !== nothing
            push!(filter_parts, "photons>$(config.min_photons)")
        end
        if config.max_precision !== nothing
            push!(filter_parts, "precision<$(config.max_precision*1000)nm")
        end
        if config.psf_sigma_mode_tolerance !== nothing
            push!(filter_parts, "PSFσ=mode±$(round(config.psf_sigma_mode_tolerance*100))%")
        end
        if config.min_pvalue !== nothing
            push!(filter_parts, "pvalue>$(config.min_pvalue)")
        end
        println(io, "- Filter: $(join(filter_parts, ", "))")

        println(io, "\n## Output Files\n")
        println(io, "- `01_detection/detection_overlay.png` - Detection visualization")
        println(io, "- `01_detection/detection_stats.md` - Detection statistics")
        println(io, "- `02_fitting/fit_quality.png` - Fit parameter histograms")
        println(io, "- `02_fitting/fit_acceptance.png` - Acceptance visualization")
        println(io, "- `02_fitting/fitting_stats.md` - Fitting statistics")
        println(io, "- `03_filtered/filter_stats.md` - Filter statistics")
        if config.frameconnect
            println(io, "- `04_frameconnect/track_histogram.png` - Track size histogram")
            println(io, "- `04_frameconnect/frameconnect_stats.md` - Frame connection statistics")
        end
        if config.calibrate_uncertainties && config.frameconnect
            println(io, "- `05_calibration/calibration_stats.md` - Uncertainty calibration parameters")
        end
        if config.drift
            println(io, "- `06_drift/drift_trajectory.png` - Drift trajectory plots")
            println(io, "- `06_drift/drift_stats.md` - Drift correction statistics")
        end
        if config.filter_isolated
            println(io, "- `07_isolated/neighbor_histogram.png` - Neighbor count histogram")
            println(io, "- `07_isolated/isolated_stats.md` - Isolated filter statistics")
        end
        if config.render
            if config.render_gaussian
                println(io, "- `08_superres/gaussian_inferno.png` - Gaussian render ($(config.render_gaussian_zoom)x)")
            end
            if config.render_histogram
                println(io, "- `08_superres/histogram_time.png` - Histogram by time ($(config.render_histogram_zoom)x)")
            end
            if config.render_circles
                println(io, "- `08_superres/circles_time.png` - Circles by time ($(config.render_circles_zoom)x)")
            end
            println(io, "- `08_superres/render_stats.md` - Render statistics")
        end
        println(io, "- `config.toml` - Analysis configuration")
        println(io, "- `results/smld_raw.h5` - Raw fitted localizations")
        println(io, "- `results/smld_final.h5` - Filtered localizations")
    end
end

"""
    simulate_detect_fit_workflow(sim_params, camera; detect_params=NamedTuple(),
                                  fit_params=NamedTuple(), workflow=nothing)
                                  → (smld_fitted, smld_ground_truth, workflow)

Complete SMLM localization workflow: simulation → detection → fitting.

This workflow:
1. Simulates SMLM data with ground truth localizations
2. Generates realistic camera images
3. Detects particles using difference of Gaussians
4. Fits Gaussian PSFs to get precise localizations

# Arguments
- `sim_params`: SMLMSim parameters (StaticSMLMParams or DiffusionSMLMParams)
- `camera`: Camera object (IdealCamera or SCMOSCamera) for image generation and detection

## Keyword Arguments
- `detect_params`: NamedTuple of parameters for `getboxes()` (boxsize, overlap, psf_sigma, min_photons, use_gpu)
- `fit_params`: NamedTuple of parameters for `GaussMLEFitter()` (psf_model, iterations, device, batch_size)
- `workflow`: Optional SMLMWorkflow object to track processing steps

# Returns
Named tuple with:
- `smld_fitted`: BasicSMLD with fitted localizations
- `smld_ground_truth`: BasicSMLD with ground truth from simulation
- `images`: Generated camera images (for inspection)
- `detections`: Raw detection results from boxer
- `fit_results`: Raw GaussMLE fit results
- `workflow`: Updated workflow object (or new one if none provided)

# Example
```julia
using SMLMAnalysis, SMLMData, SMLMSim

# Setup simulation
camera = IdealCamera(1:256, 1:256, 0.1f0)  # 256×256, 100nm pixels
sim_params = StaticSMLMParams(
    density = 1.0,      # molecules per μm²
    σ_psf = 0.13,       # PSF width in μm
    photons = 1000,     # photons per emitter
    background = 10     # background photons
)

# Run complete workflow
result = simulate_detect_fit_workflow(sim_params, camera)

# Access results
println("Ground truth: ", length(result.smld_ground_truth.emitters), " emitters")
println("Detected: ", size(result.detections.coords_pixels, 1), " particles")
println("Fitted: ", length(result.smld_fitted.emitters), " localizations")

# Inspect workflow
println(result.workflow)
```
"""
function simulate_detect_fit_workflow(
    sim_params::SMLMSim.SMLMSimParams,
    camera::SMLMData.AbstractCamera;
    pattern = nothing,
    molecule = SMLMSim.Core.GenericFluor(photons=1e4, k_off=50.0, k_on=1e-2),
    detect_params::NamedTuple = NamedTuple(),
    fit_params::NamedTuple = NamedTuple(),
    workflow::Union{SMLMWorkflow,Nothing} = nothing
)
    # Create workflow tracker if not provided
    if workflow === nothing
        workflow = SMLMWorkflow("Simulate → Detect → Fit")
    end

    # Step 1: Simulate SMLM data
    simulated = simulate(sim_params; pattern=pattern, molecule=molecule, camera=camera)

    # Handle different simulation return types
    if simulated isa Tuple && length(simulated) == 3
        # StaticSMLM returns (pattern, smld_true, smld_noisy)
        smld_ground_truth = simulated[3]  # Use noisy version as "ground truth" with applied localization noise
    else
        # DiffusionSMLM returns just SMLD
        smld_ground_truth = simulated
    end

    add_step!(workflow, "Simulate SMLM Data", :simulate,
             Dict{Symbol,Any}(:sim_params => string(typeof(sim_params))),
             :SMLMSim, summarize_smld(smld_ground_truth))

    # Step 2: Generate camera images with noise
    psf = MicroscopePSFs.GaussianPSF(sim_params.σ_psf)
    images = gen_images(smld_ground_truth, psf;
                       bg=10.0,              # Background photons/pixel
                       poisson_noise=true,   # Add shot noise
                       camera_noise=false)   # IdealCamera (no read noise)

    add_step!(workflow, "Generate Camera Images", :gen_images,
             Dict{Symbol,Any}(:psf_sigma => sim_params.σ_psf, :n_frames => size(images, 3)),
             :SMLMSim, "Image stack size $(size(images))")

    # Step 3: Detect particles (PSF-aware DoG)
    # Merge default params with user params
    default_detect = (
        boxsize = 11,
        overlap = 2.0,
        psf_sigma = sim_params.σ_psf,  # Use simulation PSF sigma
        min_photons = 500.0,           # Detection threshold in photons
        use_gpu = true
    )
    detect_params_full = merge(default_detect, detect_params)

    detections = getboxes(images, camera;
                         boxsize = detect_params_full.boxsize,
                         overlap = detect_params_full.overlap,
                         psf_sigma = detect_params_full.psf_sigma,
                         min_photons = detect_params_full.min_photons,
                         use_gpu = detect_params_full.use_gpu)

    # getboxes now returns ROIBatch directly
    roi_batch = detections

    add_step!(workflow, "Detect Particles", :getboxes,
             Dict{Symbol,Any}(pairs(detect_params_full)...),
             :SMLMBoxer, summarize_boxer_result(roi_batch))

    # Check if any detections were made
    if length(roi_batch) == 0
        @warn "No particles detected! Returning empty results. Check detection threshold (min_photons) or simulation parameters."
        empty_smld = SMLMData.BasicSMLD(
            SMLMData.Emitter2DFit{Float32}[],
            camera,
            sim_params.nframes,
            1,
            Dict{String,Any}("workflow" => "simulate_detect_fit", "warning" => "no detections")
        )
        return (
            smld_fitted = empty_smld,
            smld_ground_truth = smld_ground_truth,
            images = images,
            detections = roi_batch,
            fit_results = nothing,
            workflow = workflow
        )
    end

    # Step 5: Fit Gaussian PSFs
    # Merge default params with user params
    default_fit = (
        psf_model = GaussianXYNB(Float32(sim_params.σ_psf)),
        iterations = 20,
        device = nothing,  # auto-detect
        batch_size = 10_000
    )
    fit_params_full = merge(default_fit, fit_params)

    fitter = GaussMLEFitter(;
        psf_model = fit_params_full.psf_model,
        iterations = fit_params_full.iterations,
        device = fit_params_full.device,
        batch_size = fit_params_full.batch_size
    )

    # fit() now returns BasicSMLD directly
    smld_fitted = GaussMLE.fit(fitter, roi_batch)

    add_step!(workflow, "Fit Gaussian PSFs", :fit,
             Dict{Symbol,Any}(:psf_model => string(typeof(fit_params_full.psf_model)),
                  :iterations => fit_params_full.iterations,
                  :device => string(fitter.device)),
             :GaussMLE, summarize_smld(smld_fitted))

    # Return complete results
    return (
        smld_fitted = smld_fitted,
        smld_ground_truth = smld_ground_truth,
        images = images,
        detections = roi_batch,
        workflow = workflow
    )
end

"""
    standard_localization_workflow(images, camera; kwargs...) → (smld, workflow)

Standard localization workflow for experimental data: detection → fitting.

Similar to simulate_detect_fit_workflow but starts with existing images
rather than simulation.

# Arguments
- `images`: Camera image stack (2D or 3D array)
- `camera`: Camera object (IdealCamera or SCMOSCamera)

## Keyword Arguments
- `detect_params`: NamedTuple of detection parameters
- `fit_params`: NamedTuple of fitting parameters
- `workflow`: Optional SMLMWorkflow object

# Returns
Named tuple with:
- `smld`: BasicSMLD with fitted localizations
- `detections`: Raw detection results
- `fit_results`: Raw GaussMLE results
- `workflow`: Updated workflow object

# Example
```julia
# Load experimental data
images = load("experiment.tif")
camera = IdealCamera(1:512, 1:512, 0.065f0)

# Run workflow
result = standard_localization_workflow(images, camera;
    detect_params = (min_photons=500.0, use_gpu=true),
    fit_params = (iterations=30,)
)

# Save results
save_smld(result.smld, "localizations.h5")
```
"""
function standard_localization_workflow(
    images::AbstractArray,
    camera::SMLMData.AbstractCamera;
    detect_params::NamedTuple = NamedTuple(),
    fit_params::NamedTuple = NamedTuple(),
    workflow::Union{SMLMWorkflow,Nothing} = nothing
)
    # Create workflow tracker if not provided
    if workflow === nothing
        workflow = SMLMWorkflow("Detect → Fit")
    end

    # Step 1: Detect particles (PSF-aware DoG)
    default_detect = (
        boxsize = 11,
        overlap = 2.0,
        psf_sigma = 0.135,     # Default PSF sigma in microns
        min_photons = 500.0,   # Detection threshold in photons
        use_gpu = true
    )
    detect_params_full = merge(default_detect, detect_params)

    detections = getboxes(images, camera;
                         boxsize = detect_params_full.boxsize,
                         overlap = detect_params_full.overlap,
                         psf_sigma = detect_params_full.psf_sigma,
                         min_photons = detect_params_full.min_photons,
                         use_gpu = detect_params_full.use_gpu)

    # getboxes now returns ROIBatch directly
    roi_batch = detections

    add_step!(workflow, "Detect Particles", :getboxes,
             Dict{Symbol,Any}(pairs(detect_params_full)...),
             :SMLMBoxer, summarize_boxer_result(roi_batch))

    # Step 3: Fit Gaussian PSFs
    default_fit = (
        psf_model = GaussianXYNB(1.3f0),  # Default PSF sigma
        iterations = 20,
        device = nothing,
        batch_size = 10_000
    )
    fit_params_full = merge(default_fit, fit_params)

    fitter = GaussMLEFitter(;
        psf_model = fit_params_full.psf_model,
        iterations = fit_params_full.iterations,
        device = fit_params_full.device,
        batch_size = fit_params_full.batch_size
    )

    # fit() now returns BasicSMLD directly
    smld = GaussMLE.fit(fitter, roi_batch)

    add_step!(workflow, "Fit Gaussian PSFs", :fit,
             Dict{Symbol,Any}(:psf_model => string(typeof(fit_params_full.psf_model)),
                  :iterations => fit_params_full.iterations,
                  :device => string(fitter.device)),
             :GaussMLE, summarize_smld(smld))

    return (
        smld = smld,
        detections = roi_batch,
        workflow = workflow
    )
end
