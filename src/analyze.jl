"""
    analyze.jl

Main SMLM analysis pipeline function.
"""

using SMLMData
using SMLMBoxer
using GaussMLE
using SMLMRender
using SMLMDriftCorrection
using SMLMFrameConnection

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
result = analyze(data, camera; detect_min_photons=500.0, min_photons=500.0)

# From config file
config = load_config("myconfig.toml")
result = analyze(data, camera, config)
```
"""
function analyze(data::AbstractArray{<:Real, 3}, camera::SMLMData.AbstractCamera; kwargs...)
    return analyze(data, camera, AnalysisConfig(; kwargs...))
end

function analyze(data::AbstractArray{<:Real, 3}, camera::SMLMData.AbstractCamera, config::AnalysisConfig)
    timings = Dict{String, Float64}()
    workflow = SMLMWorkflow("SMLM Analysis")
    DC = SMLMDriftCorrection

    # Setup output directory if specified
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
        psf_sigma = config.psf_sigma,
        min_photons = config.detect_min_photons,
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

    if config.outdir !== nothing && config.save_figures
        Figures.save_detection_figures(data, roi_batch, camera, config)
        Stats.write_detection_stats(roi_batch, data, config, t)
    end

    # =========================================================================
    # Step 2: Fitting
    # =========================================================================
    print("  Fitting... ")
    device = config.fit_device == :auto ? nothing : config.fit_device

    psf_model = if config.fit_model == :fixed
        GaussianXYNB(config.psf_sigma)
    elseif config.fit_model == :variable
        GaussianXYNBS()
    elseif config.fit_model == :anisotropic
        GaussianXYNBSXSY()
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

    smld = smld_raw

    if config.outdir !== nothing && config.save_figures
        Figures.save_fitting_figures(smld_raw, roi_batch, data, camera, config, _calculate_mode)
        Stats.write_fitting_stats(smld_raw, config, t)
    end

    # =========================================================================
    # Step 3: Filtering
    # =========================================================================
    if config.filter
        print("  Filtering... ")
        n_before = length(smld.emitters)
        t = @elapsed smld = filter_smld(smld, config)
        timings["filtering"] = t
        n_after = length(smld.emitters)
        pct = round(100 * n_after / n_before, digits=1)
        println("$n_after / $n_before ($pct%) ($(round(t, digits=2))s)")

        add_step!(workflow, "Filtering", :filter,
            Dict{Symbol,Any}(:min_photons => config.min_photons, :max_precision => config.max_precision),
            :SMLMAnalysis, "$n_after / $n_before accepted")

        if config.outdir !== nothing
            Stats.write_filter_stats(smld_raw, smld, config, t, _calculate_mode)
        end
    end

    # =========================================================================
    # Step 3.5: Set dataset indices
    # =========================================================================
    if config.dataset_indices !== nothing
        print("  Setting dataset indices... ")
        for e in smld.emitters
            if 1 <= e.frame <= length(config.dataset_indices)
                e.dataset = config.dataset_indices[e.frame]
            end
        end
        n_datasets = length(unique(config.dataset_indices))
        smld = BasicSMLD(smld.emitters, smld.camera, smld.n_frames, n_datasets, smld.metadata)
        println("$n_datasets datasets")
    end

    # =========================================================================
    # Step 4: Frame Connection
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

        smld_connected = fc_result.connected
        smld = fc_result.combined
        fc_params = fc_result.params

        n_after = length(smld.emitters)
        compression = round(n_before / n_after, digits=1)
        println("$n_before → $n_after tracks ($(compression)x) ($(round(t, digits=2))s)")

        add_step!(workflow, "Frame Connection", :frameconnect,
            Dict{Symbol,Any}(:maxframegap => config.fc_maxframegap),
            :SMLMFrameConnection, "$n_after tracks from $n_before locs")

        drift_analysis = analyze_frameconnect_drift(smld_connected)

        if config.outdir !== nothing
            Figures.save_frameconnect_figures(smld_connected, config)
            Figures.save_frameconnect_drift_figures(drift_analysis, config)
            Figures.save_uncertainty_calibration_figure(drift_analysis, config)
            Stats.write_frameconnect_stats(n_before, n_after, smld_connected, fc_params, config, t;
                                           drift_analysis=drift_analysis)
        end

        # =====================================================================
        # Step 4b: Uncertainty Calibration
        # =====================================================================
        if config.calibrate_uncertainties && !isnan(drift_analysis.calibration.A) && !isnan(drift_analysis.calibration.B)
            print("  Uncertainty calibration... ")
            cal = drift_analysis.calibration

            σ_motion_nm = sqrt(max(0.0, cal.A))
            k_scale = sqrt(max(1.0, cal.B))
            σ_motion = σ_motion_nm / 1000.0

            t_cal = @elapsed begin
                smld_connected_corrected, smld_calibrated = apply_uncertainty_calibration(
                    smld_connected, σ_motion, k_scale)
                smld = smld_calibrated
            end
            timings["calibration"] = t_cal

            println("k=$(round(k_scale, digits=2)), σ_motion=$(round(σ_motion_nm, digits=1))nm ($(round(t_cal, digits=2))s)")

            add_step!(workflow, "Uncertainty Calibration", :calibrate,
                Dict{Symbol,Any}(:k_scale => k_scale, :sigma_motion => σ_motion_nm),
                :SMLMAnalysis, "k=$(round(k_scale, digits=2)), σ_motion=$(round(σ_motion_nm, digits=1))nm")

            if config.outdir !== nothing
                Stats.write_calibration_stats(σ_motion_nm, k_scale, cal, config, t_cal)
            end
        end
    end

    # =========================================================================
    # Step 5: Drift Correction
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

        # Get frame ranges per dataset
        frame_ranges = Dict{Int, Tuple{Int,Int}}()
        for ds in 1:n_datasets
            ds_frames = [e.frame for e in smld.emitters if e.dataset == ds]
            if !isempty(ds_frames)
                frame_ranges[ds] = (minimum(ds_frames), maximum(ds_frames))
            end
        end

        # Calculate max drift
        max_drift_x = 0.0
        max_drift_y = 0.0
        frames = collect(1:smld.n_frames)
        for ds in 1:n_datasets
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
            Dict{Symbol,Any}(:degree => config.drift_degree, :n_datasets => n_datasets),
            :SMLMDriftCorrection, "max $(round(max(max_drift_x, max_drift_y), digits=0))nm")

        if config.outdir !== nothing
            Figures.save_drift_figures(drift_model, smld, config, DC)
            Stats.write_drift_stats(drift_model, smld, config, t, DC)
        end

        smld = smld_corrected
    end

    # =========================================================================
    # Step 6: Isolated Emitter Filter
    # =========================================================================
    if config.filter_isolated
        print("  Isolated filter... ")
        n_before = length(smld.emitters)
        t = @elapsed smld, neighbor_counts, threshold_used = filter_isolated(smld, config)
        timings["isolated_filter"] = t
        n_after = length(smld.emitters)
        n_rejected = n_before - n_after
        pct_rejected = round(100 * n_rejected / n_before, digits=1)
        auto_str = config.isolated_min_neighbors == :auto ? " (auto)" : ""
        println("$n_rejected rejected ($pct_rejected%) threshold=$threshold_used$auto_str ($(round(t, digits=2))s)")

        add_step!(workflow, "Isolated Filter", :filter_isolated,
            Dict{Symbol,Any}(:n_sigma => config.isolated_n_sigma, :min_neighbors => threshold_used),
            :SMLMAnalysis, "$n_rejected isolated emitters rejected")

        if config.outdir !== nothing
            Figures.save_isolated_figures(neighbor_counts, threshold_used, config)
            Stats.write_isolated_stats(n_before, n_after, neighbor_counts, threshold_used, config, t)
        end
    end

    # =========================================================================
    # Step 7: Rendering
    # =========================================================================
    if config.render
        print("  Rendering... ")
        t_total = 0.0
        render_count = 0

        clip_pct = if config.render_clip_percentile === :auto
            adaptive_clip_percentile(length(smld.emitters))
        else
            config.render_clip_percentile
        end

        if config.outdir !== nothing
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

            Stats.write_render_stats(smld, config, t_total)
        else
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
            Dict{Symbol,Any}(:gaussian => config.render_gaussian, :histogram => config.render_histogram),
            :SMLMRender, "$render_count renders")
    end

    # =========================================================================
    # Save results
    # =========================================================================
    if config.outdir !== nothing && config.save_smld
        save_smld(joinpath(config.outdir, "results", "smld_final.h5"), smld)
        save_smld(joinpath(config.outdir, "results", "smld_raw.h5"), smld_raw)
        println("  Saved SMLD files to $(config.outdir)/results/")
    end

    if config.outdir !== nothing
        Stats.write_summary(roi_batch, smld_raw, smld, data, config, timings)
    end

    # Summary
    total = sum(values(timings))
    println("-"^60)
    println("  Total: $(round(total, digits=2))s")
    println("="^60)

    return AnalysisResult(smld, smld_raw, roi_batch, timings, workflow)
end
