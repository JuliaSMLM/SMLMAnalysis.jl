"""
    figures.jl

Figure generation functions for SMLM analysis pipeline.
All functions save PNG files to the appropriate output subdirectory.
"""

using CairoMakie
using Statistics
using StatsBase: countmap

"""Calculate figure size for grid overlay plots based on data aspect ratio."""
function grid_figure_size(data; n_cols=4, n_rows=3, panel_height=200)
    data_height, data_width = size(data, 1), size(data, 2)
    data_aspect = data_width / data_height
    panel_width = round(Int, panel_height * data_aspect)
    fig_width = panel_width * n_cols + 100   # columns + margins
    fig_height = panel_height * n_rows + 150  # rows + titles/margins
    return (fig_width, fig_height)
end

# =============================================================================
# Detection Figures
# =============================================================================

"""Save detection overlay figures."""
function save_detection_figures(data, roi_batch, camera, config)
    nframes = size(data, 3)
    frame_indices = [round(Int, x) for x in range(1, nframes, length=12)]

    # Intensity range
    pmin = Float64(quantile(vec(data[:,:,1]), 0.01))
    pmax = Float64(quantile(vec(data[:,:,1]), 0.99))

    fig = Figure(size=grid_figure_size(data))
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

# =============================================================================
# Fitting Figures
# =============================================================================

"""Save fitting quality figures."""
function save_fitting_figures(smld, roi_batch, data, camera, config, calculate_mode_fn)
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
        if config.min_pvalue !== nothing && config.min_pvalue > 0
            vlines!(ax5, [log10(config.min_pvalue)], color=:orange, linestyle=:solid)
        end
        if n_zero > 0
            text!(ax5, 0.02, 0.98, text="$(n_zero) fits with pvalue=0 ($(round(100*n_zero/length(pvalue), digits=1))%)",
                  align=(:left, :top), space=:relative, fontsize=10)
        end
    end

    # PSF Sigma histogram (for variable sigma fits)
    has_sigma = hasfield(typeof(emitters[1]), :σ)
    has_sigma_xy = hasfield(typeof(emitters[1]), :σx) && hasfield(typeof(emitters[1]), :σy)

    if has_sigma
        psf_sigmas = [e.σ for e in emitters]
        psf_sigma_nm = psf_sigmas .* 1000
        psf_sigma_mode = calculate_mode_fn(psf_sigmas)
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
        psf_sigmas_x = [e.σx for e in emitters]
        psf_sigmas_y = [e.σy for e in emitters]

        ax6 = Axis(fig[3, 2], xlabel="PSF σ (nm)", ylabel="Count", title="PSF Sigma Distribution (σx=blue, σy=red)",
                   limits=(50, 300, nothing, nothing))
        hist!(ax6, clamp.(psf_sigmas_x .* 1000, 50, 300), bins=range(50, 300, length=51), color=(:blue, 0.5), label="σx")
        hist!(ax6, clamp.(psf_sigmas_y .* 1000, 50, 300), bins=range(50, 300, length=51), color=(:red, 0.5), label="σy")
    else
        ax6 = Axis(fig[3, 2], xlabel="Photons", ylabel="Background (ADU)", title="Photons vs Background")
        scatter!(ax6, photons, bg, markersize=2, alpha=0.3)
    end

    save(joinpath(config.outdir, "02_fitting", "fit_quality.png"), fig)
end

# =============================================================================
# Drift Correction Figures
# =============================================================================

"""Save drift correction figures - handles single and multi-dataset modes."""
function save_drift_figures(drift_model, smld, config, DC)
    n_datasets = length(drift_model.intra)
    n_frames = smld.n_frames
    frames = collect(1:n_frames)

    # Color palette for datasets
    colors = [:blue, :red, :green, :orange, :purple, :cyan, :magenta, :brown]

    # Get frame ranges per dataset from emitters
    emitters = smld.emitters
    frame_ranges = Dict{Int, Tuple{Int,Int}}()
    for ds in 1:n_datasets
        ds_frames = [e.frame for e in emitters if e.dataset == ds]
        if !isempty(ds_frames)
            frame_ranges[ds] = (minimum(ds_frames), maximum(ds_frames))
        end
    end

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

        ax1 = Axis(fig[1, 1], xlabel="Frame", ylabel="X Drift (nm)",
                   title="X Drift per Dataset (intra-dataset)")
        ax2 = Axis(fig[1, 2], xlabel="Frame", ylabel="Y Drift (nm)",
                   title="Y Drift per Dataset (intra-dataset)")

        for (i, ds) in enumerate(1:n_datasets)
            color = colors[mod1(i, length(colors))]
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

        ax3 = Axis(fig[2, 1], xlabel="X Drift (nm)", ylabel="Y Drift (nm)",
                   title="XY Drift Paths per Dataset", aspect=DataAspect())
        for (i, ds) in enumerate(1:n_datasets)
            color = colors[mod1(i, length(colors))]
            ds_frames = haskey(frame_ranges, ds) ? collect(frame_ranges[ds][1]:frame_ranges[ds][2]) : frames
            drift_x = [DC.applydrift(0.0, f, drift_model.intra[ds].dm[1]) for f in ds_frames] .* 1000
            drift_y = [DC.applydrift(0.0, f, drift_model.intra[ds].dm[2]) for f in ds_frames] .* 1000
            lines!(ax3, drift_x, drift_y, color=color, linewidth=1.5, label="Dataset $ds")
        end
        axislegend(ax3, position=:lt)

        ax4 = Axis(fig[2, 2], title="Inter-Dataset Alignment Shifts")
        hidedecorations!(ax4)
        hidespines!(ax4)

        shift_text = "Dataset → Shift X (nm) → Shift Y (nm)\n" * "─"^40 * "\n"
        for ds in 1:n_datasets
            shift_x = drift_model.inter[ds].dm[1] * 1000
            shift_y = drift_model.inter[ds].dm[2] * 1000
            shift_text *= "   $ds    →  $(round(shift_x, digits=1))  →  $(round(shift_y, digits=1))\n"
        end
        text!(ax4, 0.5, 0.5, text=shift_text, align=(:center, :center), fontsize=14)

        save(joinpath(config.outdir, "06_drift", "drift_trajectory.png"), fig)
    end
end

# =============================================================================
# Frame Connection Figures
# =============================================================================

"""Save frame connection figures - histogram of locs per track."""
function save_frameconnect_figures(smld_connected, config)
    emitters = smld_connected.emitters
    if isempty(emitters)
        return
    end

    track_ids = [e.track_id for e in emitters]
    counts = countmap(track_ids)
    locs_per_track = collect(values(counts))

    fig = Figure(size=(800, 500))

    ax = Axis(fig[1, 1],
        xlabel = "Localizations per Track",
        ylabel = "Count (tracks)",
        title = "Frame Connection: Track Size Distribution"
    )

    max_locs = min(maximum(locs_per_track), 50)
    hist!(ax, clamp.(locs_per_track, 1, max_locs),
          bins=range(0.5, max_locs + 0.5, length=max_locs + 1),
          color=:steelblue)

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

"""Save frame connection drift analysis figures."""
function save_frameconnect_drift_figures(drift_analysis, config)
    frame_shifts = drift_analysis.frame_shifts
    chi2_values = drift_analysis.chi2_values
    n_datasets = drift_analysis.n_datasets

    if isempty(frame_shifts)
        return
    end

    fig = Figure(size=(1600, 800))

    colors_x = [:blue, :teal, :purple, :navy, :cyan]
    colors_y = [:red, :orange, :magenta, :brown, :coral]

    ax1 = Axis(fig[1, 1], xlabel="Frame", ylabel="ΔX (nm)", title="Frame-to-Frame X Shift")
    hlines!(ax1, [0], color=:gray, linestyle=:dash)

    ax2 = Axis(fig[1, 2], xlabel="Frame", ylabel="ΔY (nm)", title="Frame-to-Frame Y Shift")
    hlines!(ax2, [0], color=:gray, linestyle=:dash)

    title_suffix = n_datasets > 1 ? " (per dataset)" : ""
    ax3 = Axis(fig[2, 1], xlabel="Frame", ylabel="Cumulative Drift (nm)",
               title="Cumulative Drift from Linked Emitters$title_suffix")

    for (i, (dataset_id, shifts)) in enumerate(sort(collect(frame_shifts)))
        if isempty(shifts)
            continue
        end

        frames = [s[1] for s in shifts]
        Δx = [s[2] for s in shifts] .* 1000
        Δy = [s[3] for s in shifts] .* 1000
        σ_Δx = [s[4] for s in shifts] .* 1000
        σ_Δy = [s[5] for s in shifts] .* 1000

        color_x = colors_x[mod1(i, length(colors_x))]
        color_y = colors_y[mod1(i, length(colors_y))]

        alpha = n_datasets > 1 ? 0.15 : 0.3
        band!(ax1, frames, Δx .- σ_Δx, Δx .+ σ_Δx, color=(color_x, alpha))
        lines!(ax1, frames, Δx, color=color_x, linewidth=1)

        band!(ax2, frames, Δy .- σ_Δy, Δy .+ σ_Δy, color=(color_y, alpha))
        lines!(ax2, frames, Δy, color=color_y, linewidth=1)

        cum_x = cumsum(Δx)
        cum_y = cumsum(Δy)
        cum_σx = sqrt.(cumsum(σ_Δx .^ 2))
        cum_σy = sqrt.(cumsum(σ_Δy .^ 2))

        label_x = n_datasets > 1 ? "X (ds $dataset_id)" : "X"
        label_y = n_datasets > 1 ? "Y (ds $dataset_id)" : "Y"

        band!(ax3, frames, cum_x .- cum_σx, cum_x .+ cum_σx, color=(color_x, 0.15))
        band!(ax3, frames, cum_y .- cum_σy, cum_y .+ cum_σy, color=(color_y, 0.15))
        lines!(ax3, frames, cum_x, color=color_x, linewidth=1.5, label=label_x)
        lines!(ax3, frames, cum_y, color=color_y, linewidth=1.5, label=label_y)
    end

    axislegend(ax3, position=:lt)

    ax4 = Axis(fig[2, 2], xlabel="χ²", ylabel="Count",
               title="Uncertainty Validation (expected χ²=2)")
    if !isempty(chi2_values)
        chi2_clipped = clamp.(chi2_values, 0, 20)
        hist!(ax4, chi2_clipped, bins=range(0, 20, length=41), color=:steelblue)
        vlines!(ax4, [2.0], color=:red, linewidth=2, linestyle=:dash, label="Expected")
        vlines!(ax4, [mean(chi2_values)], color=:orange, linewidth=2, label="Observed mean")
        axislegend(ax4, position=:rt)

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
function save_uncertainty_calibration_figure(drift_analysis, config)
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

    scatter!(ax, cal.bin_centers, cal.bin_observed, markersize=8, color=:steelblue, label="Binned data")

    x_range = range(minimum(cal.bin_centers), maximum(cal.bin_centers), length=100)
    lines!(ax, collect(x_range), collect(x_range), color=:gray, linestyle=:dash, linewidth=2, label="Perfect (1:1)")

    y_fit = cal.A .+ cal.B .* x_range
    lines!(ax, collect(x_range), collect(y_fit), color=:red, linewidth=2, label="Fit: A + B×σ²")

    axislegend(ax, position=:lt)

    σ_motion = cal.A > 0 ? sqrt(cal.A) : 0.0
    k_scale = sqrt(cal.B)

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
        text = "A (motion²) = $(round(cal.A, digits=1)) nm²\nB (k²) = $(round(cal.B, digits=2))\nR² = $(round(cal.r_squared, digits=3))\n\n$(interpretation)",
        align = (:right, :bottom),
        space = :relative,
        fontsize = 11
    )

    save(joinpath(config.outdir, "05_calibration", "uncertainty_calibration.png"), fig)
end

# =============================================================================
# Isolated Filter Figures
# =============================================================================

"""Save isolated emitter filter figures - neighbor count histogram with triangle method visualization."""
function save_isolated_figures(neighbor_counts, threshold, config)
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

    max_count = min(maximum(neighbor_counts), 50)
    hist!(ax, clamp.(neighbor_counts, 0, max_count), bins=range(0, max_count, length=max_count+1))

    vlines!(ax, [threshold - 0.5], color=:red, linestyle=:dash, linewidth=2,
            label="Threshold ($threshold)")

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
