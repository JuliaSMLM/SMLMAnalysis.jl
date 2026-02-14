"""
Frame connection step - uses SMLMFrameConnection.FrameConnectConfig directly.

No wrapper config -- the upstream FrameConnectConfig is used as a pipeline step.
Calibration is optionally enabled via FrameConnectConfig.calibration field.
"""

"""
    frameconnect_step(smld, cfg; outdir, step_number, verbose) -> (combined_smld, FrameConnectInfo)

Run frame connection on `smld`, returning the combined (recombined) SMLD as the
primary result, plus the upstream FrameConnectInfo.

When `cfg.calibration` is set, uncertainty calibration is applied before combination:
link -> calibrate -> combine (single pass with correct weights).

# Returns
`(combined_smld, FrameConnectInfo)`
"""
function frameconnect_step(smld::BasicSMLD, cfg::SMLMFrameConnection.FrameConnectConfig;
                           outdir::Union{String,Nothing}=nothing,
                           step_number::Int=0,
                           verbose::Int=Verbosity.STANDARD)
    v = verbose
    dir = step_outdir(outdir, step_number, cfg)

    v >= Verbosity.PROGRESS && @info "[$step_number] $(step_name(cfg))" max_frame_gap=cfg.max_frame_gap

    n_before = length(smld.emitters)

    # Config dispatch: returns (combined, FrameConnectInfo) where FrameConnectInfo contains .connected
    t = @elapsed (combined, connect_info) = SMLMFrameConnection.frameconnect(smld, cfg)

    smld_connected = connect_info.connected

    n_after = length(combined.emitters)
    compression = n_before / n_after

    if dir !== nothing
        _save_frameconnect_outputs!(dir, cfg, v, t, n_before, n_after,
                                    smld_connected, connect_info)
    end

    cal = connect_info.calibration
    if cal !== nothing && cal.calibration_applied && v >= Verbosity.PROGRESS
        @info "  Calibration: k=$(round(cal.k_scale, digits=2)), sigma_motion=$(round(cal.sigma_motion_nm, digits=1))nm, chi2=$(round(cal.mean_chi2, digits=2))"
    end

    v >= Verbosity.PROGRESS && @info "  -> $n_after tracks ($(round(compression, digits=1))x) ($(round(t, digits=2))s)"
    (combined, connect_info)
end

function _step_summary(info::SMLMFrameConnection.FrameConnectInfo)
    d = Dict{Symbol,Any}(
        :n_before => info.n_input,
        :n_after => info.n_combined,
        :compression => round(info.n_input / max(1, info.n_combined), digits=1),
    )
    cal = info.calibration
    if cal !== nothing && cal.calibration_applied
        d[:k_scale] = round(cal.k_scale, digits=2)
        d[:sigma_motion_nm] = round(cal.sigma_motion_nm, digits=1)
        d[:mean_chi2] = round(cal.mean_chi2, digits=2)
    end
    d
end

"""
    analyze(smld, cfg::SMLMFrameConnection.FrameConnectConfig; kwargs...) -> (combined_smld, StepInfo)

Run frame connection on localizations.
"""
function analyze(smld::BasicSMLD, cfg::SMLMFrameConnection.FrameConnectConfig;
                 outdir=nothing, step_number::Int=0, verbose::Int=Verbosity.STANDARD, kwargs...)
    t = @elapsed (combined, connect_info) = frameconnect_step(smld, cfg;
        outdir=outdir, step_number=step_number, verbose=verbose)
    (combined, StepInfo(step_number, cfg, t, _step_summary(connect_info); info=connect_info))
end

function _save_frameconnect_outputs!(dir::String, cfg::SMLMFrameConnection.FrameConnectConfig,
                             v::Int, t::Float64,
                             n_before::Int, n_after::Int,
                             smld_connected::BasicSMLD, connect_info)
    mkpath(dir)
    _save_config!(dir, cfg)
    _save_info!(dir, connect_info)
    if connect_info.calibration !== nothing
        _save_info!(dir, connect_info.calibration; section="calibration")
    end

    if v >= Verbosity.STANDARD
        koff_stats = _save_frameconnect_figures(dir, smld_connected)
        _write_frameconnect_stats(dir, cfg, n_before, n_after, t, koff_stats, connect_info.calibration)
        if connect_info.calibration !== nothing
            _save_calibration_figures(dir, connect_info.calibration, smld_connected)
        end
    end
end

function _write_frameconnect_stats(dir, cfg, n_before, n_after, t, koff_stats=nothing, cal=nothing)
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

        # Calibration results
        if cal !== nothing && cal.calibration_applied
            println(io, "")
            println(io, "## Uncertainty Calibration")
            println(io, "- **k (CRLB scale)**: $(round(cal.k_scale, digits=2))")
            println(io, "- **sigma_motion**: $(round(cal.sigma_motion_nm, digits=1)) nm")
            println(io, "- **Mean chi2**: $(round(cal.mean_chi2, digits=2)) (expected: 2.0)")
            println(io, "- **R2**: $(round(cal.r_squared, digits=3))")
            println(io, "- **Frame pairs**: $(cal.n_pairs)")
            if cal.n_tracks_filtered > 0
                println(io, "- **Tracks filtered (chi2)**: $(cal.n_tracks_filtered)")
            end
        elseif cal !== nothing && !cal.calibration_applied
            println(io, "")
            println(io, "## Uncertainty Calibration")
            println(io, "- **Status**: skipped ($(cal.warning))")
        end

        println(io, "")
        println(io, "## Parameters")
        println(io, "- max_frame_gap: $(cfg.max_frame_gap)")
        println(io, "- max_sigma_dist: $(cfg.max_sigma_dist)")
        if cfg.calibration !== nothing
            println(io, "- calibration.clamp_k_to_one: $(cfg.calibration.clamp_k_to_one)")
            if cfg.calibration.filter_high_chi2
                println(io, "- calibration.filter_high_chi2: true")
                println(io, "- calibration.chi2_filter_threshold: $(cfg.calibration.chi2_filter_threshold)")
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
    k_off = 1.0 / mean_length

    (mean_length=mean_length, k_off=k_off)
end

function _save_frameconnect_figures(dir, smld_connected)
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
end

# ============================================================
# Calibration diagnostic plots (data sourced from CalibrationResult)
# ============================================================

"""
Save calibration diagnostic plots from FrameConnectInfo.calibration.
Called when calibration was enabled and outdir is set.
"""
function _save_calibration_figures(dir, cal::SMLMFrameConnection.CalibrationResult,
                                   smld_connected::BasicSMLD)
    !cal.calibration_applied && return
    isempty(cal.bin_centers) && return

    _save_calibration_plot(dir, cal)
    _save_shift_histogram(dir, cal)
    _save_drift_jitter_plot(dir, smld_connected)
end

"""
Plot uncertainty calibration: binned observed vs CRLB variance with WLS fit line.
FC stores bin data in μm²; convert to nm² for readable axes.
"""
function _save_calibration_plot(dir, cal::SMLMFrameConnection.CalibrationResult)
    # Convert from μm² to nm² (×1e6)
    bc_nm2 = cal.bin_centers .* 1e6
    bo_nm2 = cal.bin_observed .* 1e6
    A_nm2 = cal.A * 1e6
    B = cal.B  # dimensionless ratio

    fig = Figure(size=(700, 500))
    ax = Axis(fig[1, 1],
        xlabel = "CRLB Variance (nm², per-axis average)",
        ylabel = "Observed Variance (nm², per-axis average)",
        title = "Uncertainty Calibration: obs = A + B * CRLB"
    )

    # Data points
    scatter!(ax, bc_nm2, bo_nm2, markersize=8, label="Binned data")

    # 1:1 line
    x_range = range(minimum(bc_nm2), maximum(bc_nm2), length=100)
    lines!(ax, collect(x_range), collect(x_range), color=:gray, linestyle=:dash, label="1:1 (ideal)")

    # Fit line
    lines!(ax, collect(x_range), A_nm2 .+ B .* collect(x_range), color=:red, linewidth=2,
        label="Fit: A=$(round(A_nm2, sigdigits=3)) nm², B=$(round(B, digits=2))")

    axislegend(ax, position=:lt)

    # Annotation with clamped/filtered info
    k_str = "k = $(round(cal.k_scale, digits=2))"
    if cal.k_scale == 1.0 && B < 1.0
        k_str *= " (clamped)"
    end
    σ_str = "sigma_motion = $(round(cal.sigma_motion_nm, digits=1)) nm"
    if cal.k_scale == 1.0 && cal.sigma_motion_nm == 0.0
        σ_str *= " (A clamped)"
    end
    pairs_str = "$(cal.n_pairs) pairs"
    if cal.n_tracks_filtered > 0
        filt_pct = round(100 * cal.n_tracks_filtered / (cal.n_tracks_used + cal.n_tracks_filtered), digits=1)
        pairs_str *= " ($(filt_pct)% filtered)"
    end
    pairs_str *= ", R2=$(round(cal.r_squared, digits=3))"

    text!(ax, 0.95, 0.05,
        text="$k_str\n$σ_str\n$pairs_str",
        align=(:right, :bottom),
        space=:relative,
        fontsize=12)

    save(joinpath(dir, "uncertainty_calibration.png"), fig)
end

"""
Plot histogram of frame-to-frame shifts from calibration analysis.
"""
function _save_shift_histogram(dir, cal::SMLMFrameConnection.CalibrationResult)
    isempty(cal.frame_shifts) && return

    # Collect all shifts across datasets
    all_dx = Float64[]
    all_dy = Float64[]
    for (_, shifts) in cal.frame_shifts
        for (dx, dy) in shifts
            push!(all_dx, dx * 1000)  # Convert to nm
            push!(all_dy, dy * 1000)
        end
    end

    isempty(all_dx) && return

    fig = Figure(size=(900, 400))

    ax1 = Axis(fig[1, 1], xlabel="DX (nm)", ylabel="Count", title="Frame-to-Frame X Shifts")
    hist!(ax1, all_dx, bins=50, color=(:blue, 0.6))
    vlines!(ax1, [0], color=:red, linestyle=:dash)

    ax2 = Axis(fig[1, 2], xlabel="DY (nm)", ylabel="Count", title="Frame-to-Frame Y Shifts")
    hist!(ax2, all_dy, bins=50, color=(:red, 0.6))
    vlines!(ax2, [0], color=:blue, linestyle=:dash)

    # Annotate with std dev
    σ_dx = std(all_dx)
    σ_dy = std(all_dy)
    text!(ax1, 0.95, 0.95, text="sigma = $(round(σ_dx, digits=1)) nm",
        align=(:right, :top), space=:relative, fontsize=11)
    text!(ax2, 0.95, 0.95, text="sigma = $(round(σ_dy, digits=1)) nm",
        align=(:right, :top), space=:relative, fontsize=11)

    save(joinpath(dir, "shift_histogram.png"), fig)
end

"""
Plot frame-to-frame drift jitter and cumulative drift from linked emitters.
Three panels: jitter X, jitter Y (top), cumulative drift (bottom).
Uses global frame numbering (absolute_frame across datasets).
"""
function _save_drift_jitter_plot(dir, smld_connected::BasicSMLD)
    emitters = smld_connected.emitters
    isempty(emitters) && return
    n_frames = smld_connected.n_frames

    # Group emitters by track_id
    track_dict = Dict{Int, Vector{eltype(emitters)}}()
    for e in emitters
        e.track_id > 0 || continue
        if haskey(track_dict, e.track_id)
            push!(track_dict[e.track_id], e)
        else
            track_dict[e.track_id] = [e]
        end
    end

    # Collect weighted mean shift per (dataset, frame) transition
    # Accumulate: key = (dataset, frame) => (sum_wx, sum_wy, sum_w_x, sum_w_y, count)
    shift_accum = Dict{Tuple{Int,Int}, NTuple{5,Float64}}()

    for (_, track_locs) in track_dict
        sort!(track_locs, by = e -> (e.dataset, e.frame))
        for i in 1:(length(track_locs) - 1)
            e1, e2 = track_locs[i], track_locs[i + 1]
            e1.dataset == e2.dataset && e2.frame == e1.frame + 1 || continue

            dx = Float64(e2.x - e1.x)
            dy = Float64(e2.y - e1.y)
            vx = Float64(e1.σ_x^2 + e2.σ_x^2)
            vy = Float64(e1.σ_y^2 + e2.σ_y^2)
            (vx > 0 && vy > 0) || continue

            w_x, w_y = 1.0 / vx, 1.0 / vy
            key = (e1.dataset, e1.frame)
            if haskey(shift_accum, key)
                s = shift_accum[key]
                shift_accum[key] = (s[1] + w_x * dx, s[2] + w_y * dy,
                                    s[3] + w_x, s[4] + w_y, s[5] + 1.0)
            else
                shift_accum[key] = (w_x * dx, w_y * dy, w_x, w_y, 1.0)
            end
        end
    end

    isempty(shift_accum) && return

    # Build per-dataset arrays with global frame numbering
    datasets = sort(unique(first.(keys(shift_accum))))

    # Per-dataset data: (global_frames, dx_nm, dy_nm, σ_dx_nm, σ_dy_nm)
    dataset_data = Dict{Int, NamedTuple{(:gframes, :dx, :dy, :σ_dx, :σ_dy),
                         Tuple{Vector{Int}, Vector{Float64}, Vector{Float64},
                               Vector{Float64}, Vector{Float64}}}}()

    for ds in datasets
        frames = Int[]
        dx_nm = Float64[]
        dy_nm = Float64[]
        σ_dx_nm = Float64[]
        σ_dy_nm = Float64[]
        for (key, s) in shift_accum
            key[1] == ds || continue
            # Global frame = (dataset-1) * n_frames + frame
            gframe = (ds - 1) * n_frames + key[2]
            push!(frames, gframe)
            push!(dx_nm, s[1] / s[3] * 1000)   # weighted mean, μm -> nm
            push!(dy_nm, s[2] / s[4] * 1000)
            push!(σ_dx_nm, 1000.0 / sqrt(s[3])) # σ of weighted mean
            push!(σ_dy_nm, 1000.0 / sqrt(s[4]))
        end
        perm = sortperm(frames)
        dataset_data[ds] = (gframes=frames[perm], dx=dx_nm[perm], dy=dy_nm[perm],
                            σ_dx=σ_dx_nm[perm], σ_dy=σ_dy_nm[perm])
    end

    # Build global sorted arrays for cumulative drift
    all_gframes = Int[]
    all_dx = Float64[]
    all_dy = Float64[]
    all_σ_dx = Float64[]
    all_σ_dy = Float64[]
    for (_, dd) in dataset_data
        append!(all_gframes, dd.gframes)
        append!(all_dx, dd.dx)
        append!(all_dy, dd.dy)
        append!(all_σ_dx, dd.σ_dx)
        append!(all_σ_dy, dd.σ_dy)
    end
    gperm = sortperm(all_gframes)
    all_gframes = all_gframes[gperm]
    all_dx = all_dx[gperm]
    all_dy = all_dy[gperm]
    all_σ_dx = all_σ_dx[gperm]
    all_σ_dy = all_σ_dy[gperm]

    # Cumulative drift and uncertainty
    cum_dx = cumsum(all_dx)
    cum_dy = cumsum(all_dy)
    cum_σ_dx = sqrt.(cumsum(all_σ_dx .^ 2))
    cum_σ_dy = sqrt.(cumsum(all_σ_dy .^ 2))

    colors = [:blue, :red, :green, :orange, :purple, :cyan, :magenta, :brown]

    fig = Figure(size=(900, 700))

    # Top row: jitter
    ax1 = Axis(fig[1, 1], xlabel="Global Frame", ylabel="DX (nm)",
               title="Frame-to-Frame X Shift (jitter)")
    ax2 = Axis(fig[1, 2], xlabel="Global Frame", ylabel="DY (nm)",
               title="Frame-to-Frame Y Shift (jitter)")

    for (di, ds) in enumerate(datasets)
        dd = dataset_data[ds]
        c = colors[mod1(di, length(colors))]
        lines!(ax1, dd.gframes, dd.dx, color=(c, 0.5), linewidth=0.5)
        lines!(ax2, dd.gframes, dd.dy, color=(c, 0.5), linewidth=0.5)
    end
    hlines!(ax1, [0], color=:black, linestyle=:dash, linewidth=0.5)
    hlines!(ax2, [0], color=:black, linestyle=:dash, linewidth=0.5)

    # Bottom row: cumulative drift with confidence band
    ax3 = Axis(fig[2, 1:2], xlabel="Global Frame", ylabel="Cumulative Drift (nm)",
               title="Cumulative Drift from Linked Emitters")

    band!(ax3, all_gframes, cum_dx .- cum_σ_dx, cum_dx .+ cum_σ_dx, color=(:red, 0.2))
    lines!(ax3, all_gframes, cum_dx, color=:red, linewidth=1, label="X")
    band!(ax3, all_gframes, cum_dy .- cum_σ_dy, cum_dy .+ cum_σ_dy, color=(:blue, 0.2))
    lines!(ax3, all_gframes, cum_dy, color=:blue, linewidth=1, label="Y")
    hlines!(ax3, [0], color=:black, linestyle=:dash, linewidth=0.5)
    axislegend(ax3, position=:lt)

    save(joinpath(dir, "drift_jitter.png"), fig)
end
