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
    sigma_small::Float64 = 1.0
    sigma_large::Float64 = 2.0
    minval::Float64 = 10.0
    detect_gpu::Bool = true

    # === Fitting (always runs) ===
    psf_sigma::Float32 = 1.3f0
    iterations::Int = 20
    fit_device::Symbol = :auto      # :cpu, :gpu, :auto

    # === Filtering (runs before frame connection) ===
    filter::Bool = true
    min_photons::Union{Float64, Nothing} = 500.0
    max_sigma::Union{Float64, Nothing} = nothing  # microns, e.g., 0.03
    min_pvalue::Union{Float64, Nothing} = nothing

    # === Frame Connection (runs on filtered data) ===
    frameconnect::Bool = true
    fc_max_distance::Float64 = 0.1  # microns
    fc_max_gap::Int = 1             # frames
    fc_max_blinks::Int = 1

    # === Rendering ===
    render::Bool = true             # ON by default - primary output
    render_zoom::Int = 20
    render_colormap::Symbol = :inferno
    render_strategy::Symbol = :gaussian  # :gaussian or :histogram

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
    # Step 1: Detection
    # =========================================================================
    print("  Detection... ")
    t = @elapsed roi_batch = getboxes(data, camera;
        boxsize = config.boxsize,
        overlap = config.overlap,
        sigma_small = config.sigma_small,
        sigma_large = config.sigma_large,
        minval = config.minval,
        use_gpu = config.detect_gpu
    )
    timings["detection"] = t
    println("$(length(roi_batch)) ROIs ($(round(t, digits=2))s)")

    add_step!(workflow, "Detection", :getboxes,
        Dict{Symbol,Any}(:boxsize => config.boxsize, :minval => config.minval),
        :SMLMBoxer, "$(length(roi_batch)) ROIs")

    if length(roi_batch) == 0
        error("No particles detected! Check minval threshold (currently $(config.minval))")
    end

    # Save detection figures
    if config.outdir !== nothing && config.save_figures
        _save_detection_figures(data, roi_batch, camera, config)
    end

    # =========================================================================
    # Step 2: Fitting
    # =========================================================================
    print("  Fitting... ")
    device = config.fit_device == :auto ? nothing : config.fit_device
    fitter = GaussMLEFitter(
        psf_model = GaussianXYNB(config.psf_sigma),
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

    # Save fitting figures
    if config.outdir !== nothing && config.save_figures
        _save_fitting_figures(smld_raw, roi_batch, data, camera, config)
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
    end

    # =========================================================================
    # Step 4: Frame Connection (optional) - operates on filtered data
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
    # Step 5: Rendering (optional)
    # =========================================================================
    if config.render
        print("  Rendering... ")
        strategy = config.render_strategy == :gaussian ? GaussianRender() : HistogramRender()

        if config.outdir !== nothing
            t = @elapsed render(smld;
                strategy = strategy,
                zoom = config.render_zoom,
                colormap = config.render_colormap,
                filename = joinpath(config.outdir, "05_superres", "superres_$(config.render_strategy).png")
            )
            timings["rendering"] = t
            println("saved ($(round(t, digits=2))s)")
        else
            t = @elapsed render(smld;
                strategy = strategy,
                zoom = config.render_zoom,
                colormap = config.render_colormap
            )
            timings["rendering"] = t
            println("($(round(t, digits=2))s)")
        end

        add_step!(workflow, "Rendering", :render,
            Dict{Symbol,Any}(:zoom => config.render_zoom, :strategy => config.render_strategy),
            :SMLMRender, "$(config.render_zoom)x zoom")
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
    σ_x = [e.σ_x for e in emitters]
    σ_y = [e.σ_y for e in emitters]

    # Fit quality histograms
    fig = Figure(size=(1200, 400))

    ax1 = Axis(fig[1, 1], xlabel="Photons", ylabel="Count", title="Photon Distribution")
    hist!(ax1, photons, bins=50)

    ax2 = Axis(fig[1, 2], xlabel="σ_x (μm)", ylabel="Count", title="PSF Width X")
    hist!(ax2, σ_x, bins=50)

    ax3 = Axis(fig[1, 3], xlabel="σ_y (μm)", ylabel="Count", title="PSF Width Y")
    hist!(ax3, σ_y, bins=50)

    save(joinpath(config.outdir, "02_fitting", "fit_quality.png"), fig)

    # Fit acceptance panel
    precision_values = [sqrt(e.σ_x^2 + e.σ_y^2)/sqrt(2) for e in emitters]

    photon_ok = config.min_photons === nothing ? trues(length(emitters)) : photons .> config.min_photons
    sigma_ok = config.max_sigma === nothing ? trues(length(emitters)) :
               [max(e.σ_x, e.σ_y) < config.max_sigma for e in emitters]
    accepted = photon_ok .& sigma_ok

    n_total = length(emitters)
    n_accepted = sum(accepted)
    accept_pct = round(100 * n_accepted / n_total, digits=1)

    nframes = size(data, 3)
    frame_indices = [round(Int, x) for x in range(1, nframes, length=12)]
    pmin = Float64(quantile(vec(data[:,:,1]), 0.01))
    pmax = Float64(quantile(vec(data[:,:,1]), 0.99))

    fig = Figure(size=(2400, 700))
    Label(fig[0, 1:4],
        "Fit Acceptance: $n_accepted/$n_total ($accept_pct%) — Green=Accepted, Orange=σ, Red=Photons",
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

            for pass in [false, true]
                for j in eachindex(frame_locs)
                    if frame_accepted[j] == pass
                        bx, by = det_x[j], det_y[j]
                        c = frame_accepted[j] ? :green : (!frame_photon_ok[j] ? :red : :orange)
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

    # Step 3: Detect particles
    # Merge default params with user params
    default_detect = (
        boxsize = 11,
        overlap = 2.0,
        sigma_small = 1.0,
        sigma_large = 2.0,
        minval = 10.0,
        use_gpu = true
    )
    detect_params_full = merge(default_detect, detect_params)

    detections = getboxes(images, camera;
                         boxsize = detect_params_full.boxsize,
                         overlap = detect_params_full.overlap,
                         sigma_small = detect_params_full.sigma_small,
                         sigma_large = detect_params_full.sigma_large,
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

    # Step 1: Detect particles
    default_detect = (
        boxsize = 11,
        overlap = 2.0,
        sigma_small = 1.0,
        sigma_large = 2.0,
        minval = 10.0,
        use_gpu = true
    )
    detect_params_full = merge(default_detect, detect_params)

    detections = getboxes(images, camera;
                         boxsize = detect_params_full.boxsize,
                         overlap = detect_params_full.overlap,
                         sigma_small = detect_params_full.sigma_small,
                         sigma_large = detect_params_full.sigma_large,
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
