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
                 outdir=nothing, step_number::Int=0, verbose::Int=Verbosity.STANDARD)
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

    if v >= Verbosity.STANDARD
        koff_stats = _save_frameconnect_figures(dir, smld_connected)
        _write_frameconnect_stats(dir, cfg, n_before, n_after, t, koff_stats, connect_info.calibration)
        if connect_info.calibration !== nothing
            _save_calibration_figures(dir, connect_info.calibration)
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
function _save_calibration_figures(dir, cal::SMLMFrameConnection.CalibrationResult)
    !cal.calibration_applied && return
    isempty(cal.bin_centers) && return

    _save_calibration_plot(dir, cal)
    _save_shift_histogram(dir, cal)
end

"""
Plot uncertainty calibration: binned observed vs CRLB variance with WLS fit line.
"""
function _save_calibration_plot(dir, cal::SMLMFrameConnection.CalibrationResult)
    fig = Figure(size=(700, 500))
    ax = Axis(fig[1, 1],
        xlabel = "CRLB Variance (combined x+y)",
        ylabel = "Observed Variance (combined x+y)",
        title = "Uncertainty Calibration: obs = A + B * CRLB"
    )

    # Data points
    scatter!(ax, cal.bin_centers, cal.bin_observed, markersize=8, label="Binned data")

    # 1:1 line
    x_range = range(minimum(cal.bin_centers), maximum(cal.bin_centers), length=100)
    lines!(ax, collect(x_range), collect(x_range), color=:gray, linestyle=:dash, label="1:1 (ideal)")

    # Fit line
    lines!(ax, collect(x_range), cal.A .+ cal.B .* collect(x_range), color=:red, linewidth=2,
        label="Fit: A=$(round(cal.A, sigdigits=3)), B=$(round(cal.B, digits=2))")

    axislegend(ax, position=:lt)

    # Annotation
    text!(ax, 0.95, 0.05,
        text="k = $(round(cal.k_scale, digits=2))\nsigma_motion = $(round(cal.sigma_motion_nm, digits=1)) nm\n$(cal.n_pairs) pairs, R2=$(round(cal.r_squared, digits=3))",
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
