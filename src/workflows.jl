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
using CairoMakie
using Statistics
using TOML

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
    minval = 100.0,
    min_photons = 500.0,
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
    minval::Float64 = 10.0
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

    # === Filtering (runs before drift correction) ===
    filter::Bool = true
    min_photons::Union{Float64, Nothing} = 500.0
    max_sigma::Union{Float64, Nothing} = 0.015    # 15nm precision threshold
    min_pvalue::Union{Float64, Nothing} = 1e-3    # p-value threshold

    # === Drift Correction (runs on filtered data) ===
    drift::Bool = true              # ON by default for DNA-PAINT
    drift_degree::Int = 2           # Polynomial degree (2 usually sufficient)
    drift_cost_fun::String = "Kdtree"  # "Kdtree" (fast) or "Entropy"
    drift_model::String = "Polynomial" # "Polynomial" or "LegendrePoly"

    # === Frame Connection (runs on drift-corrected data) ===
    frameconnect::Bool = false      # OFF by default (skip for most workflows)
    fc_max_distance::Float64 = 0.1  # microns
    fc_max_gap::Int = 1             # frames
    fc_max_blinks::Int = 1

    # === Rendering ===
    render::Bool = true             # ON by default - primary output
    render_gaussian::Bool = true    # Gaussian blur @ 20x + inferno
    render_histogram::Bool = true   # Histogram @ 10x + time coloring
    render_circles::Bool = true     # Circles @ 50x + time coloring
    render_gaussian_zoom::Int = 20
    render_histogram_zoom::Int = 10
    render_circles_zoom::Int = 50
    render_time_colormap::Symbol = :turbo  # colormap for time-colored renders

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

# With custom parameters
result = analyze(data, camera; minval=100.0, min_photons=500.0)

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
    if config.outdir !== nothing
        mkpath(config.outdir)
        mkpath(joinpath(config.outdir, "01_detection"))
        mkpath(joinpath(config.outdir, "02_fitting"))
        mkpath(joinpath(config.outdir, "03_filtered"))
        if config.frameconnect
            mkpath(joinpath(config.outdir, "04_frameconnect"))
        end
        if config.render
            mkpath(joinpath(config.outdir, "05_superres"))
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
        minval = config.minval,
        use_gpu = config.detect_gpu
    )
    timings["detection"] = t
    println("$(length(roi_batch)) ROIs ($(round(t, digits=2))s)")

    add_step!(workflow, "Detection", :getboxes,
        Dict{Symbol,Any}(:boxsize => config.boxsize, :psf_sigma => config.psf_sigma, :minval => config.minval),
        :SMLMBoxer, "$(length(roi_batch)) ROIs")

    if length(roi_batch) == 0
        error("No particles detected! Check minval threshold (currently $(config.minval))")
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
    t = @elapsed smld_raw = fit(fitter, roi_batch)
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
            Dict{Symbol,Any}(:min_photons => config.min_photons, :max_sigma => config.max_sigma),
            :SMLMAnalysis, "$n_after / $n_before accepted")

        # Save filter stats
        if config.outdir !== nothing
            _write_filter_stats(smld_raw, smld, config, t)
        end
    end

    # =========================================================================
    # Step 4: Drift Correction (optional) - operates on filtered data
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

        # Extract drift curves from model (more accurate than coord differencing)
        n_frames = smld.n_frames
        frames = collect(1:n_frames)
        DC = SMLMDriftCorrection
        drift_x = [DC.applydrift(0.0, f, drift_model.intra[1].dm[1]) for f in frames]
        drift_y = [DC.applydrift(0.0, f, drift_model.intra[1].dm[2]) for f in frames]

        max_drift_x = maximum(abs.(drift_x)) * 1000  # nm
        max_drift_y = maximum(abs.(drift_y)) * 1000  # nm
        println("max $(round(max_drift_x, digits=1))nm X, $(round(max_drift_y, digits=1))nm Y ($(round(t, digits=2))s)")

        add_step!(workflow, "Drift Correction", :driftcorrect,
            Dict{Symbol,Any}(:degree => config.drift_degree, :cost_fun => config.drift_cost_fun),
            :SMLMDriftCorrection, "max $(round(max(max_drift_x, max_drift_y), digits=0))nm")

        # Save drift figures and stats
        if config.outdir !== nothing
            mkpath(joinpath(config.outdir, "04_drift"))
            _save_drift_figures(frames, drift_x, drift_y, config)
            _write_drift_stats(frames, drift_x, drift_y, config, t)
        end

        smld = smld_corrected
    end

    # =========================================================================
    # Step 5: Frame Connection (optional) - operates on drift-corrected data
    # =========================================================================
    if config.frameconnect
        print("  Frame connection... ")
        # TODO: Integrate SMLMFrameConnection when dependency resolved
        # For now, skip with warning
        println("(skipped - package not yet integrated)")
        timings["frameconnect"] = 0.0
        # t = @elapsed smld = frameconnect(smld;
        #     max_distance=config.fc_max_distance,
        #     max_gap=config.fc_max_gap)
        # timings["frameconnect"] = t
        # add_step!(workflow, "Frame Connection", :frameconnect, ...)
    end

    # =========================================================================
    # Step 6: Rendering (optional) - generates multiple render types
    # =========================================================================
    if config.render
        print("  Rendering... ")
        t_total = 0.0
        render_count = 0

        if config.outdir !== nothing
            # 1. Gaussian render @ 20x + inferno (primary super-res image)
            if config.render_gaussian
                t = @elapsed render(smld;
                    strategy = GaussianRender(),
                    zoom = config.render_gaussian_zoom,
                    colormap = :inferno,
                    filename = joinpath(config.outdir, "05_superres", "gaussian_inferno.png")
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
                    filename = joinpath(config.outdir, "05_superres", "histogram_time.png")
                )
                t_total += t
                render_count += 1
            end

            # 3. Circles render @ 50x + time coloring (individual loc inspection)
            if config.render_circles
                t = @elapsed render(smld;
                    strategy = CircleRender(),
                    zoom = config.render_circles_zoom,
                    color_by = :frame,
                    colormap = config.render_time_colormap,
                    filename = joinpath(config.outdir, "05_superres", "circles_time.png")
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
                colormap = :inferno
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

"""Filter SMLD based on config criteria."""
function _filter_smld(smld::BasicSMLD, config::AnalysisConfig)
    emitters = smld.emitters
    mask = trues(length(emitters))

    if config.min_photons !== nothing
        mask .&= [e.photons > config.min_photons for e in emitters]
    end

    if config.max_sigma !== nothing
        mask .&= [max(e.σ_x, e.σ_y) < config.max_sigma for e in emitters]
    end

    if config.min_pvalue !== nothing
        mask .&= [e.pvalue > config.min_pvalue for e in emitters]
    end

    filtered = emitters[mask]
    return BasicSMLD(filtered, smld.camera, smld.n_frames, smld.n_datasets, smld.metadata)
end

"""Save detection overlay figures."""
function _save_detection_figures(data, roi_batch, camera, config)
    nframes = size(data, 3)
    frame_indices = [round(Int, x) for x in range(1, nframes, length=12)]

    # Intensity range
    pmin = Float64(quantile(vec(data[:,:,1]), 0.01))
    pmax = Float64(quantile(vec(data[:,:,1]), 0.99))

    # Wide figure for data aspect ratio
    fig = Figure(size=(2400, 700))
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

    # Photons histogram
    ax1 = Axis(fig[1, 1], xlabel="Photons", ylabel="Count", title="Photon Distribution")
    hist!(ax1, photons, bins=50)
    vlines!(ax1, [median(photons)], color=:red, linestyle=:dash)

    # Background histogram
    ax2 = Axis(fig[1, 2], xlabel="Background (ADU)", ylabel="Count", title="Background Distribution")
    hist!(ax2, bg, bins=50)
    vlines!(ax2, [median(bg)], color=:red, linestyle=:dash)

    # σ_x histogram (in nm) - fixed range 0-50nm for precision
    σ_x_nm = σ_x .* 1000
    ax3 = Axis(fig[2, 1], xlabel="σ_x (nm)", ylabel="Count", title="X Precision Distribution",
               limits=(0, 50, nothing, nothing))
    hist!(ax3, clamp.(σ_x_nm, 0, 50), bins=range(0, 50, length=51))
    vlines!(ax3, [median(σ_x_nm)], color=:red, linestyle=:dash)
    if config.max_sigma !== nothing
        vlines!(ax3, [config.max_sigma * 1000], color=:orange, linestyle=:solid, label="Threshold")
    end

    # σ_y histogram (in nm) - fixed range 0-50nm for precision
    σ_y_nm = σ_y .* 1000
    ax4 = Axis(fig[2, 2], xlabel="σ_y (nm)", ylabel="Count", title="Y Precision Distribution",
               limits=(0, 50, nothing, nothing))
    hist!(ax4, clamp.(σ_y_nm, 0, 50), bins=range(0, 50, length=51))
    vlines!(ax4, [median(σ_y_nm)], color=:red, linestyle=:dash)
    if config.max_sigma !== nothing
        vlines!(ax4, [config.max_sigma * 1000], color=:orange, linestyle=:solid)
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

    # Photons vs Background scatter
    ax6 = Axis(fig[3, 2], xlabel="Photons", ylabel="Background (ADU)", title="Photons vs Background")
    scatter!(ax6, photons, bg, markersize=2, alpha=0.3)

    save(joinpath(config.outdir, "02_fitting", "fit_quality.png"), fig)

    # Fit acceptance panel
    precision_values = [sqrt(e.σ_x^2 + e.σ_y^2)/sqrt(2) for e in emitters]

    photon_ok = config.min_photons === nothing ? trues(length(emitters)) : photons .> config.min_photons
    sigma_ok = config.max_sigma === nothing ? trues(length(emitters)) :
               [max(e.σ_x, e.σ_y) < config.max_sigma for e in emitters]
    pvalue_ok = config.min_pvalue === nothing ? trues(length(emitters)) : pvalue .> config.min_pvalue
    accepted = photon_ok .& sigma_ok .& pvalue_ok

    n_total = length(emitters)
    n_accepted = sum(accepted)
    accept_pct = round(100 * n_accepted / n_total, digits=1)

    nframes = size(data, 3)
    frame_indices = [round(Int, x) for x in range(1, nframes, length=12)]
    pmin = Float64(quantile(vec(data[:,:,1]), 0.01))
    pmax = Float64(quantile(vec(data[:,:,1]), 0.99))

    fig = Figure(size=(2400, 700))
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
            frame_sigma_ok = sigma_ok[frame_mask]
            frame_pvalue_ok = pvalue_ok[frame_mask]

            for pass in [false, true]
                for j in eachindex(frame_locs)
                    if frame_accepted[j] == pass
                        bx, by = det_x[j], det_y[j]
                        # Color: green=accepted, red=photons, orange=sigma, purple=pvalue
                        c = if frame_accepted[j]
                            :green
                        elseif !frame_photon_ok[j]
                            :red
                        elseif !frame_sigma_ok[j]
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
        println(io, "- Min value threshold: $(config.minval)")
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

    # Calculate per-filter pass rates
    photon_ok = config.min_photons === nothing ? trues(n_raw) : photons .> config.min_photons
    sigma_ok = config.max_sigma === nothing ? trues(n_raw) :
               [max(e.σ_x, e.σ_y) < config.max_sigma for e in emitters]
    pvalue_ok = config.min_pvalue === nothing ? trues(n_raw) : pvalue .> config.min_pvalue

    n_photon_pass = sum(photon_ok)
    n_sigma_pass = sum(sigma_ok)
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
        println(io, "## Per-Filter Results\n")
        println(io, "| Filter | Threshold | Pass | Fail | % Pass |")
        println(io, "|--------|-----------|------|------|--------|")

        if config.min_photons !== nothing
            println(io, "| Photons | >$(config.min_photons) | $(n_photon_pass) | $(n_raw - n_photon_pass) | $(round(100*n_photon_pass/n_raw, digits=1))% |")
        end
        if config.max_sigma !== nothing
            println(io, "| Precision | <$(config.max_sigma*1000)nm | $(n_sigma_pass) | $(n_raw - n_sigma_pass) | $(round(100*n_sigma_pass/n_raw, digits=1))% |")
        end
        if config.min_pvalue !== nothing
            println(io, "| P-value | >$(config.min_pvalue) | $(n_pvalue_pass) | $(n_raw - n_pvalue_pass) | $(round(100*n_pvalue_pass/n_raw, digits=1))% |")
        end

        println(io, "")
        println(io, "## Limiting Factor")
        min_pass = min(n_photon_pass, n_sigma_pass, n_pvalue_pass)
        if min_pass == n_pvalue_pass && config.min_pvalue !== nothing
            println(io, "**P-value filter** is the most restrictive ($(round(100*n_pvalue_pass/n_raw, digits=1))% pass)")
        elseif min_pass == n_sigma_pass && config.max_sigma !== nothing
            println(io, "**Precision filter** is the most restrictive ($(round(100*n_sigma_pass/n_raw, digits=1))% pass)")
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

"""Save drift correction figures."""
function _save_drift_figures(frames, drift_x, drift_y, config)
    # Convert to nm
    drift_x_nm = drift_x .* 1000
    drift_y_nm = drift_y .* 1000

    fig = Figure(size=(1400, 400))

    # X drift vs frame
    ax1 = Axis(fig[1, 1], xlabel="Frame", ylabel="X Drift (nm)", title="X Drift vs Frame")
    lines!(ax1, frames, drift_x_nm, color=:blue, linewidth=1.5)
    hlines!(ax1, [0], color=:gray, linestyle=:dash)

    # Y drift vs frame
    ax2 = Axis(fig[1, 2], xlabel="Frame", ylabel="Y Drift (nm)", title="Y Drift vs Frame")
    lines!(ax2, frames, drift_y_nm, color=:red, linewidth=1.5)
    hlines!(ax2, [0], color=:gray, linestyle=:dash)

    # XY drift path
    ax3 = Axis(fig[1, 3], xlabel="X Drift (nm)", ylabel="Y Drift (nm)",
               title="XY Drift Path", aspect=DataAspect())
    lines!(ax3, drift_x_nm, drift_y_nm, color=:black, linewidth=1.5)
    scatter!(ax3, [drift_x_nm[1]], [drift_y_nm[1]], color=:green, markersize=12, label="Start")
    scatter!(ax3, [drift_x_nm[end]], [drift_y_nm[end]], color=:red, markersize=12, label="End")
    axislegend(ax3, position=:lt)

    save(joinpath(config.outdir, "04_drift", "drift_trajectory.png"), fig)
end

"""Write drift correction statistics markdown file."""
function _write_drift_stats(frames, drift_x, drift_y, config, elapsed_time)
    # Convert to nm
    drift_x_nm = drift_x .* 1000
    drift_y_nm = drift_y .* 1000

    max_x = maximum(abs.(drift_x_nm))
    max_y = maximum(abs.(drift_y_nm))
    total_x = drift_x_nm[end] - drift_x_nm[1]
    total_y = drift_y_nm[end] - drift_y_nm[1]
    total_dist = sqrt(total_x^2 + total_y^2)

    filepath = joinpath(config.outdir, "04_drift", "drift_stats.md")
    open(filepath, "w") do io
        println(io, "# Drift Correction Statistics\n")
        println(io, "## Summary")
        println(io, "- **Time**: $(round(elapsed_time, digits=2))s")
        println(io, "- **Frames**: $(length(frames))")
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
        println(io, "- Drift rate: $(round(total_dist / length(frames) * 1000, digits=2)) nm/kframe")
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

    filepath = joinpath(config.outdir, "05_superres", "render_stats.md")
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
    psf_sigma = has_sigma ? median([e.σ for e in emitters_raw]) * 1000 : nothing

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
        end
        println(io, "Precision:    $(round(median(σ_x)*1000, digits=1)) nm median      ", median(σ_x)*1000 < 20 ? "✓" : "⚠")
        println(io, "Photons:      $(round(median(photons), digits=0)) median      ", 1000 < median(photons) < 100000 ? "✓" : "⚠")
        println(io, "Filtering:    $(round(100*acceptance_rate, digits=1))% accepted      ", 0.01 < acceptance_rate < 0.5 ? "✓" : "⚠")
        println(io, "```")

        # Overall status
        all_ok = (50 < rois_per_frame < 500) &&
                 (pval_pass_rate > 0.05) &&
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
        if config.drift
            println(io, "| Drift Correction | - | $(round(get(timings, "drift", 0.0), digits=1))s |")
        end
        if config.render
            println(io, "| Rendering | - | $(round(get(timings, "rendering", 0.0), digits=1))s |")
        end
        println(io, "| **Total** | **$(n_filtered) localizations** | **$(round(sum(values(timings)), digits=1))s** |")

        println(io, "\n## Configuration\n")
        println(io, "- Frames: $(nframes)")
        println(io, "- Fit model: $(config.fit_model)")
        println(io, "- Filter: photons>$(something(config.min_photons, "none")), σ<$(something(config.max_sigma, "none") === nothing ? "none" : "$(config.max_sigma*1000)nm"), pvalue>$(something(config.min_pvalue, "none"))")

        println(io, "\n## Output Files\n")
        println(io, "- `01_detection/detection_overlay.png` - Detection visualization")
        println(io, "- `01_detection/detection_stats.md` - Detection statistics")
        println(io, "- `02_fitting/fit_quality.png` - Fit parameter histograms")
        println(io, "- `02_fitting/fit_acceptance.png` - Acceptance visualization")
        println(io, "- `02_fitting/fitting_stats.md` - Fitting statistics")
        println(io, "- `03_filtered/filter_stats.md` - Filter statistics")
        if config.drift
            println(io, "- `04_drift/drift_trajectory.png` - Drift trajectory plots")
            println(io, "- `04_drift/drift_stats.md` - Drift correction statistics")
        end
        if config.render
            if config.render_gaussian
                println(io, "- `05_superres/gaussian_inferno.png` - Gaussian render ($(config.render_gaussian_zoom)x)")
            end
            if config.render_histogram
                println(io, "- `05_superres/histogram_time.png` - Histogram by time ($(config.render_histogram_zoom)x)")
            end
            if config.render_circles
                println(io, "- `05_superres/circles_time.png` - Circles by time ($(config.render_circles_zoom)x)")
            end
            println(io, "- `05_superres/render_stats.md` - Render statistics")
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
- `detect_params`: NamedTuple of parameters for `getboxes()` (boxsize, overlap, sigma_small, sigma_large, minval, use_gpu)
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
        minval = 10.0,
        use_gpu = true
    )
    detect_params_full = merge(default_detect, detect_params)

    detections = getboxes(images, camera;
                         boxsize = detect_params_full.boxsize,
                         overlap = detect_params_full.overlap,
                         psf_sigma = detect_params_full.psf_sigma,
                         minval = detect_params_full.minval,
                         use_gpu = detect_params_full.use_gpu)

    # getboxes now returns ROIBatch directly
    roi_batch = detections

    add_step!(workflow, "Detect Particles", :getboxes,
             Dict{Symbol,Any}(pairs(detect_params_full)...),
             :SMLMBoxer, summarize_boxer_result(roi_batch))

    # Check if any detections were made
    if length(roi_batch) == 0
        @warn "No particles detected! Returning empty results. Check detection threshold (minval) or simulation parameters."
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
    smld_fitted = fit(fitter, roi_batch)

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
    detect_params = (minval=50.0, use_gpu=true),
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
        psf_sigma = 0.135,  # Default PSF sigma in microns
        minval = 10.0,
        use_gpu = true
    )
    detect_params_full = merge(default_detect, detect_params)

    detections = getboxes(images, camera;
                         boxsize = detect_params_full.boxsize,
                         overlap = detect_params_full.overlap,
                         psf_sigma = detect_params_full.psf_sigma,
                         minval = detect_params_full.minval,
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
    smld = fit(fitter, roi_batch)

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
