"""
Frame connection step - wraps SMLMFrameConnection.frameconnect
"""

"""
    FrameConnectConfig <: AbstractSMLMConfig

Frame connection step. Links localizations of the same emitter across
consecutive frames using SMLMFrameConnection, then optionally calibrates
uncertainties.

# Keywords
- `max_frame_gap`: Maximum dark frames allowed in a track (default: 5)
- `max_sigma_dist`: Spatial matching threshold in sigma (default: 5.0)
- `calibrate`: Run uncertainty calibration after connection (default: `true`)
- `clamp_k_to_one`: Prevent CRLB scale factor k < 1 (default: `true`)
- `filter_high_chi2`: Remove tracks with high chi-squared pairs (default: `false`)
- `chi2_filter_threshold`: Chi-squared threshold for track removal (default: 6.0)
"""
@kwdef struct FrameConnectConfig <: SMLMData.AbstractSMLMConfig
    # SMLMFrameConnection.frameconnect kwargs
    max_frame_gap::Int = 5
    max_sigma_dist::Float64 = 5.0
    n_density_neighbors::Int = 2
    max_neighbors::Int = 2
    # Extra options
    calibrate::Bool = true  # Run uncertainty calibration after
    clamp_k_to_one::Bool = true  # Don't allow k < 1 (CRLB is theoretical lower bound)
    # Chi² filtering - removes tracks with high chi² pairs (likely double-emitter fits)
    filter_high_chi2::Bool = false
    chi2_filter_threshold::Float64 = 6.0  # ~99.7% of chi²(2) is below this
end

"""
    frameconnect_step(smld, cfg; outdir, step_number, verbose) -> (combined_smld, info_nt)

Run frame connection on `smld`, returning the combined (recombined) SMLD as the
primary result, plus a NamedTuple with step metadata.

# Returns
`(combined_smld, (step_record, smld_connected, calibration_result, connect_info))`
"""
function frameconnect_step(smld::BasicSMLD, cfg::FrameConnectConfig;
                           outdir::Union{String,Nothing}=nothing,
                           step_number::Int=0,
                           verbose::Int=Verbosity.STANDARD)
    v = verbose
    dir = step_outdir(outdir, step_number, cfg)

    v >= Verbosity.PROGRESS && @info "[$step_number] $(step_name(cfg))" max_frame_gap=cfg.max_frame_gap

    n_before = length(smld.emitters)

    # Tuple-pattern: returns (combined, FrameConnectInfo) where FrameConnectInfo contains .connected
    t = @elapsed (combined, connect_info) = SMLMFrameConnection.frameconnect(smld;
        max_frame_gap = cfg.max_frame_gap,
        max_sigma_dist = cfg.max_sigma_dist,
        n_density_neighbors = cfg.n_density_neighbors,
        max_neighbors = cfg.max_neighbors
    )

    smld_connected = connect_info.connected

    # Optional chi² filtering to remove tracks with high chi² pairs (likely double-emitter fits)
    n_tracks_filtered = 0
    n_locs_filtered = 0
    if cfg.filter_high_chi2
        n_locs_before = length(smld_connected.emitters)
        smld_connected, combined, n_tracks_filtered = _filter_high_chi2_tracks(
            smld_connected, cfg.chi2_filter_threshold)
        n_locs_filtered = n_locs_before - length(smld_connected.emitters)
        v >= Verbosity.PROGRESS && @info "  Chi² filter: removed $n_tracks_filtered tracks ($n_locs_filtered localizations)"
    end

    n_after = length(combined.emitters)
    compression = n_before / n_after

    summary = Dict{Symbol,Any}(
        :n_before => n_before,
        :n_after => n_after,
        :compression => round(compression, digits=1),
        :n_tracks_filtered => n_tracks_filtered,
        :n_locs_filtered => n_locs_filtered
    )

    # Optional uncertainty calibration
    calibration_result = nothing
    if cfg.calibrate
        combined, calibration_result = _analyze_and_calibrate(smld_connected, combined, cfg, summary)
    end

    record = StepRecord(step_number, cfg, t, summary; info=connect_info)

    if dir !== nothing
        _save_frameconnect_outputs!(dir, cfg, v, t, n_before, n_after,
                                    smld_connected, smld.n_frames,
                                    calibration_result, summary, connect_info)
    end

    v >= Verbosity.PROGRESS && @info "  → $n_after tracks ($(round(compression, digits=1))x) ($(round(t, digits=2))s)"
    (combined, (step_record=record, smld_connected=smld_connected, calibration_result=calibration_result, connect_info=connect_info))
end

"""
    _filter_high_chi2_tracks(smld_connected, threshold) -> (filtered_connected, filtered_combined, n_removed)

Remove tracks that have any frame-to-frame pair with chi² > threshold.
These are likely double-emitter fits where position shifts between frames.

Returns filtered smld_connected, recombined smld, and count of removed tracks.
"""
function _filter_high_chi2_tracks(smld_connected, threshold::Float64)
    emitters = smld_connected.emitters
    EmitterType = eltype(emitters)

    # Group emitters by track_id
    track_dict = Dict{Int, Vector{EmitterType}}()
    for e in emitters
        if e.track_id > 0
            if !haskey(track_dict, e.track_id)
                track_dict[e.track_id] = EmitterType[]
            end
            push!(track_dict[e.track_id], e)
        end
    end

    # Find tracks with any high chi² pair
    bad_tracks = Set{Int}()

    for (track_id, track_emitters) in track_dict
        # Sort by dataset then frame
        sort!(track_emitters, by = e -> (e.dataset, e.frame))

        # Check consecutive frame pairs within same dataset
        for i in 1:(length(track_emitters) - 1)
            e1, e2 = track_emitters[i], track_emitters[i + 1]

            # Only consecutive frames within same dataset
            if e2.dataset == e1.dataset && e2.frame == e1.frame + 1
                Δx = Float64(e2.x - e1.x)
                Δy = Float64(e2.y - e1.y)
                var_x = Float64(e1.σ_x^2 + e2.σ_x^2)
                var_y = Float64(e1.σ_y^2 + e2.σ_y^2)

                if var_x > 0 && var_y > 0
                    χ2 = Δx^2 / var_x + Δy^2 / var_y
                    if χ2 > threshold
                        push!(bad_tracks, track_id)
                        break  # No need to check more pairs in this track
                    end
                end
            end
        end
    end

    n_removed = length(bad_tracks)

    # Filter out bad tracks
    good_emitters = filter(e -> e.track_id == 0 || !(e.track_id in bad_tracks), emitters)

    filtered_connected = BasicSMLD(good_emitters, smld_connected.camera,
                                    smld_connected.n_frames, smld_connected.n_datasets,
                                    smld_connected.metadata)

    # Recombine remaining tracks
    filtered_combined = recombine_tracks(filtered_connected)

    return filtered_connected, filtered_combined, n_removed
end

"""
    _analyze_and_calibrate(smld_connected, smld_combined, cfg, summary) -> (calibrated_smld, drift_analysis)

Run uncertainty calibration on connected tracks. Pure function -- returns the
calibrated SMLD and drift analysis result without mutating anything except the
`summary` dict (for logging).
"""
function _analyze_and_calibrate(smld_connected::BasicSMLD, smld_combined::BasicSMLD,
                                 cfg::FrameConnectConfig, summary::Dict{Symbol,Any})
    drift_analysis = analyze_frameconnect_drift(smld_connected)
    summary[:mean_chi2] = round(drift_analysis.mean_chi2, digits=2)

    cal = drift_analysis.calibration
    calibrated_smld = smld_combined
    if !isnan(cal.A) && !isnan(cal.B)
        # σ_motion: always clamp A >= 0 (negative motion variance is unphysical)
        A_clamped = cal.A < 0.0
        σ_motion_nm = sqrt(max(0.0, cal.A))

        # k_scale: sqrt(B) but optionally clamped to >= 1.0
        if cfg.clamp_k_to_one
            k_scale = sqrt(max(1.0, cal.B))
        else
            k_scale = sqrt(max(0.0, cal.B))  # Still prevent negative/NaN
        end

        σ_motion = σ_motion_nm / 1000.0

        _, calibrated_smld = apply_uncertainty_calibration(smld_connected, σ_motion, k_scale)

        summary[:k_scale] = round(k_scale, digits=2)
        summary[:sigma_motion_nm] = round(σ_motion_nm, digits=1)
        summary[:k_clamped] = cfg.clamp_k_to_one && cal.B < 1.0
        summary[:A_clamped] = A_clamped
        summary[:A_raw] = cal.A
    end

    (calibrated_smld, drift_analysis)
end

function _save_frameconnect_outputs!(dir::String, cfg::FrameConnectConfig, v::Int, t::Float64,
                             n_before::Int, n_after::Int,
                             smld_connected::BasicSMLD, n_frames_per_dataset::Int,
                             cal_result, summary, connect_info)
    mkpath(dir)
    _save_config!(dir, cfg)
    _save_info!(dir, connect_info)

    if v >= Verbosity.STANDARD
        koff_stats = _save_frameconnect_figures(dir, smld_connected, cal_result)
        _write_frameconnect_stats(dir, cfg, n_before, n_after, t, cal_result, summary, koff_stats)
        # Calibration and drift plots at STANDARD level
        if cal_result !== nothing
            _save_calibration_plot(dir, cal_result, summary)
            _save_drift_jitter_plot(dir, cal_result, n_frames_per_dataset)
        end
    end

    if v >= Verbosity.DETAILED && cal_result !== nothing
        _save_frameconnect_detailed(dir, cal_result)
    end
end

function _write_frameconnect_stats(dir, cfg, n_before, n_after, t, cal_result, summary, koff_stats=nothing)
    compression = n_before / n_after

    filepath = joinpath(dir, "stats.md")
    open(filepath, "w") do io
        println(io, "# Frame Connection Statistics\n")
        println(io, "## Summary")
        println(io, "- **Input localizations**: $n_before")
        println(io, "- **Output tracks**: $n_after")
        println(io, "- **Compression**: $(round(compression, digits=1))x")
        println(io, "- **Time**: $(round(t, digits=2))s")

        # k_off estimation from track lengths
        if koff_stats !== nothing && !isnan(koff_stats.k_off)
            println(io, "")
            println(io, "## Blinking Kinetics (from track lengths)")
            println(io, "- **Mean track length**: $(round(koff_stats.mean_length, digits=2)) frames")
            println(io, "- **Estimated k_off**: $(round(koff_stats.k_off, digits=3)) /frame")
        end
        println(io, "")
        println(io, "## Parameters")
        println(io, "- max_frame_gap: $(cfg.max_frame_gap)")
        println(io, "- max_sigma_dist: $(cfg.max_sigma_dist)")

        # Chi² filtering info
        if cfg.filter_high_chi2
            println(io, "")
            println(io, "## Chi² Track Filtering")
            println(io, "- Threshold: χ² > $(cfg.chi2_filter_threshold)")
            println(io, "- Tracks removed: $(summary[:n_tracks_filtered])")
            println(io, "- Localizations removed: $(summary[:n_locs_filtered])")
        end

        if cal_result !== nothing && !isnan(cal_result.mean_chi2)
            println(io, "")
            println(io, "## Uncertainty Calibration")
            println(io, "- Mean χ²: $(round(cal_result.mean_chi2, digits=2)) (expected: 2.0)")
            cal = cal_result.calibration
            if !isnan(cal.A)
                k_raw = sqrt(max(0.0, cal.B))
                k_used = get(summary, :k_scale, k_raw)
                k_clamped = get(summary, :k_clamped, false)
                k_str = k_clamped ? "$(round(k_used, digits=2)) (clamped from $(round(k_raw, digits=2)))" : "$(round(k_used, digits=2))"
                println(io, "- k (CRLB scale): $k_str")

                A_clamped = get(summary, :A_clamped, false)
                σ_motion_nm = sqrt(max(0.0, cal.A))
                if A_clamped
                    A_raw = get(summary, :A_raw, cal.A)
                    println(io, "- σ_motion: $(round(σ_motion_nm, digits=1)) nm (A clamped from $(round(A_raw, digits=1)) nm²)")
                else
                    println(io, "- σ_motion: $(round(σ_motion_nm, digits=1)) nm")
                end
                println(io, "")
                println(io, "### Calibration Fit Details")
                println(io, "- Total frame pairs: $(cal.n_total)")
                println(io, "- Filtered (χ² > $(cal.chi2_threshold)): $(cal.n_filtered) ($(round(100*cal.n_filtered/cal.n_total, digits=1))%)")
                println(io, "- Used for fit: $(cal.n_total - cal.n_filtered)")
                println(io, "- Fit: A=$(round(cal.A, digits=1)) nm², B=$(round(cal.B, digits=3))")
                println(io, "- R²: $(round(cal.r_squared, digits=3))")
            end
        end
    end
end

"""
Estimate k_off from track length distribution.
For geometric distribution: mean = 1/k_off, so k_off = 1/mean
"""
function _estimate_koff(locs_per_track::Vector{Int})
    isempty(locs_per_track) && return (mean_length=NaN, k_off=NaN)

    mean_length = mean(locs_per_track)
    # k_off = probability of turning off per frame
    # For geometric distribution: E[n] = 1/p, so p = 1/E[n]
    k_off = 1.0 / mean_length

    (mean_length=mean_length, k_off=k_off)
end

function _save_frameconnect_figures(dir, smld_connected, cal_result)
    emitters = smld_connected.emitters
    isempty(emitters) && return nothing

    # Track size histogram
    track_ids = [e.track_id for e in emitters]
    counts = Dict{Int,Int}()
    for tid in track_ids
        counts[tid] = get(counts, tid, 0) + 1
    end
    locs_per_track = collect(values(counts))

    # Estimate k_off from track lengths
    koff_stats = _estimate_koff(locs_per_track)

    fig = Figure(size=(800, 500))
    ax = Axis(fig[1, 1], xlabel="Localizations per Track", ylabel="Count", title="Track Size Distribution")
    max_locs = min(maximum(locs_per_track), 50)
    hist!(ax, clamp.(locs_per_track, 1, max_locs), bins=max_locs)
    vlines!(ax, [median(locs_per_track)], color=:red, linestyle=:dash, label="Median")
    vlines!(ax, [koff_stats.mean_length], color=:blue, linestyle=:dot, label="Mean")

    # Add k_off annotation
    text!(ax, 0.95, 0.95,
        text="Mean: $(round(koff_stats.mean_length, digits=2)) frames\nk_off: $(round(koff_stats.k_off, digits=3)) /frame",
        align=(:right, :top),
        space=:relative,
        fontsize=12)

    axislegend(ax, position=:rt)
    save(joinpath(dir, "track_histogram.png"), fig)

    return koff_stats

    # Chi2 histogram if available
    if cal_result !== nothing && !isempty(cal_result.chi2_values)
        fig2 = Figure(size=(700, 500))
        ax2 = Axis(fig2[1, 1], xlabel="χ²", ylabel="Density", title="Uncertainty Validation (χ² Distribution)")

        chi2_vals = cal_result.chi2_values
        chi2_clipped = clamp.(chi2_vals, 0, 15)

        # Histogram normalized to density
        hist!(ax2, chi2_clipped, bins=30, normalization=:pdf, color=(:steelblue, 0.7), label="Observed")

        # Theoretical χ²(2) distribution: f(x) = 0.5 * exp(-x/2)
        x_theory = range(0.01, 15, length=200)
        y_theory = 0.5 .* exp.(-x_theory ./ 2)
        lines!(ax2, collect(x_theory), y_theory, color=:red, linewidth=2, label="χ²(2) theory")

        # Vertical lines for means
        obs_mean = mean(chi2_vals)
        vlines!(ax2, [2.0], color=:red, linestyle=:dash, linewidth=1.5, label="Expected mean (2.0)")
        vlines!(ax2, [obs_mean], color=:orange, linewidth=2, label="Observed mean ($(round(obs_mean, digits=2)))")

        axislegend(ax2, position=:rt)

        # Interpretation text
        if obs_mean > 2.5
            interp = "Uncertainties underestimated"
        elseif obs_mean < 1.5
            interp = "Uncertainties overestimated"
        else
            interp = "Uncertainties well-calibrated"
        end
        text!(ax2, 0.97, 0.55, text=interp, align=(:right, :top), space=:relative, fontsize=12)

        save(joinpath(dir, "chi2_histogram.png"), fig2)
    end
end

function _save_calibration_plot(dir, cal_result, summary)
    # Uncertainty calibration plot with chi² fit
    cal = cal_result.calibration
    isempty(cal.bin_centers) && return

    fig = Figure(size=(700, 500))
    ax = Axis(fig[1, 1],
        xlabel = "Reported Variance (nm²)",
        ylabel = "Observed Variance (nm²)",
        title = "Uncertainty Calibration: χ² = A + B·σ²"
    )

    # Data points
    scatter!(ax, cal.bin_centers, cal.bin_observed, markersize=8, label="Data")

    # 1:1 line
    x_range = range(minimum(cal.bin_centers), maximum(cal.bin_centers), length=100)
    lines!(ax, collect(x_range), collect(x_range), color=:gray, linestyle=:dash, label="1:1 (ideal)")

    # Fit line
    A, B = cal.A, cal.B
    k_raw = sqrt(max(0.0, B))
    k_used = get(summary, :k_scale, k_raw)
    k_clamped = get(summary, :k_clamped, false)
    σ_motion = sqrt(max(0.0, A))
    lines!(ax, collect(x_range), A .+ B .* collect(x_range), color=:red, linewidth=2,
        label="Fit: A=$(round(A, digits=1)), B=$(round(B, digits=2))")

    axislegend(ax, position=:lt)

    # Add text annotation with calibration results
    n_used = cal.n_total - cal.n_filtered
    filter_pct = round(100 * cal.n_filtered / cal.n_total, digits=1)
    k_str = k_clamped ? "$(round(k_used, digits=2)) (clamped)" : "$(round(k_used, digits=2))"

    A_clamped = get(summary, :A_clamped, false)
    σ_str = A_clamped ? "$(round(σ_motion, digits=1)) nm (A clamped)" : "$(round(σ_motion, digits=1)) nm"

    text!(ax, 0.95, 0.05,
        text="k = $k_str\nσ_motion = $σ_str\n$(n_used) pairs ($(filter_pct)% filtered)",
        align=(:right, :bottom),
        space=:relative,
        fontsize=12)

    save(joinpath(dir, "uncertainty_calibration.png"), fig)
end

function _save_frameconnect_detailed(dir, cal_result)
    # Additional detailed outputs (kept for DETAILED verbosity)
end

"""
Plot frame-to-frame drift/jitter estimated from linked emitters.
Shows both instantaneous shifts and cumulative drift using global frame index.
"""
function _save_drift_jitter_plot(dir, cal_result, n_frames_per_dataset::Int)
    frame_shifts = cal_result.frame_shifts
    isempty(frame_shifts) && return

    n_datasets = length(frame_shifts)

    fig = Figure(size=(1200, 800))

    colors_x = [:blue, :teal, :purple, :navy]
    colors_y = [:red, :orange, :magenta, :brown]

    # Frame-to-frame X shift
    ax1 = Axis(fig[1, 1], xlabel="Global Frame", ylabel="ΔX (nm)",
               title="Frame-to-Frame X Shift (jitter)")
    hlines!(ax1, [0], color=:gray, linestyle=:dash)

    # Frame-to-frame Y shift
    ax2 = Axis(fig[1, 2], xlabel="Global Frame", ylabel="ΔY (nm)",
               title="Frame-to-Frame Y Shift (jitter)")
    hlines!(ax2, [0], color=:gray, linestyle=:dash)

    # Cumulative drift
    ax3 = Axis(fig[2, 1:2], xlabel="Global Frame", ylabel="Cumulative Drift (nm)",
               title="Cumulative Drift from Linked Emitters")
    hlines!(ax3, [0], color=:gray, linestyle=:dash)

    # Collect all data for continuous cumulative plot
    all_global_frames = Int[]
    all_Δx = Float64[]
    all_Δy = Float64[]
    all_σ_Δx = Float64[]
    all_σ_Δy = Float64[]

    for (i, (dataset_id, shifts)) in enumerate(sort(collect(frame_shifts)))
        isempty(shifts) && continue

        # Convert per-dataset frames to global frames
        local_frames = [s[1] for s in shifts]
        global_frames = local_frames .+ (dataset_id - 1) * n_frames_per_dataset

        Δx = [s[2] for s in shifts] .* 1000  # Convert to nm
        Δy = [s[3] for s in shifts] .* 1000
        σ_Δx = [s[4] for s in shifts] .* 1000
        σ_Δy = [s[5] for s in shifts] .* 1000

        color_x = colors_x[mod1(i, length(colors_x))]
        color_y = colors_y[mod1(i, length(colors_y))]

        alpha = n_datasets > 1 ? 0.15 : 0.3

        # Plot instantaneous shifts with uncertainty bands
        band!(ax1, global_frames, Δx .- σ_Δx, Δx .+ σ_Δx, color=(color_x, alpha))
        lines!(ax1, global_frames, Δx, color=color_x, linewidth=0.5)

        band!(ax2, global_frames, Δy .- σ_Δy, Δy .+ σ_Δy, color=(color_y, alpha))
        lines!(ax2, global_frames, Δy, color=color_y, linewidth=0.5)

        # Collect for cumulative
        append!(all_global_frames, global_frames)
        append!(all_Δx, Δx)
        append!(all_Δy, Δy)
        append!(all_σ_Δx, σ_Δx)
        append!(all_σ_Δy, σ_Δy)
    end

    # Mark dataset boundaries
    for ds in 2:n_datasets
        boundary = (ds - 1) * n_frames_per_dataset
        vlines!(ax1, [boundary], color=:lightgray, linestyle=:dot)
        vlines!(ax2, [boundary], color=:lightgray, linestyle=:dot)
        vlines!(ax3, [boundary], color=:lightgray, linestyle=:dot)
    end

    # Sort by global frame for cumulative calculation
    if !isempty(all_global_frames)
        sort_idx = sortperm(all_global_frames)
        sorted_frames = all_global_frames[sort_idx]
        sorted_Δx = all_Δx[sort_idx]
        sorted_Δy = all_Δy[sort_idx]
        sorted_σ_Δx = all_σ_Δx[sort_idx]
        sorted_σ_Δy = all_σ_Δy[sort_idx]

        # Cumulative drift (continuous across all datasets)
        cum_x = cumsum(sorted_Δx)
        cum_y = cumsum(sorted_Δy)
        cum_σx = sqrt.(cumsum(sorted_σ_Δx .^ 2))
        cum_σy = sqrt.(cumsum(sorted_σ_Δy .^ 2))

        band!(ax3, sorted_frames, cum_x .- cum_σx, cum_x .+ cum_σx, color=(:blue, 0.15))
        band!(ax3, sorted_frames, cum_y .- cum_σy, cum_y .+ cum_σy, color=(:red, 0.15))
        lines!(ax3, sorted_frames, cum_x, color=:blue, linewidth=1.5, label="X")
        lines!(ax3, sorted_frames, cum_y, color=:red, linewidth=1.5, label="Y")
    end

    axislegend(ax3, position=:lt)

    save(joinpath(dir, "drift_jitter.png"), fig)
end
