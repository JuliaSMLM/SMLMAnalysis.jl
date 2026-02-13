"""
Frame connection step - uses SMLMFrameConnection.FrameConnectConfig directly.

No wrapper config -- the upstream FrameConnectConfig is used as a pipeline step.
Calibration and chi2 filtering are now in CalibrationConfig (calibration_step.jl).
"""

"""
    frameconnect_step(smld, cfg; outdir, step_number, verbose) -> (combined_smld, info_nt)

Run frame connection on `smld`, returning the combined (recombined) SMLD as the
primary result, plus a NamedTuple with step metadata.

# Returns
`(combined_smld, (step_record, smld_connected, connect_info))`
"""
function frameconnect_step(smld::BasicSMLD, cfg::SMLMFrameConnection.FrameConnectConfig;
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

    n_after = length(combined.emitters)
    compression = n_before / n_after

    summary = Dict{Symbol,Any}(
        :n_before => n_before,
        :n_after => n_after,
        :compression => round(compression, digits=1),
    )

    record = StepRecord(step_number, cfg, t, summary; info=connect_info)

    if dir !== nothing
        _save_frameconnect_outputs!(dir, cfg, v, t, n_before, n_after,
                                    smld_connected, connect_info)
    end

    v >= Verbosity.PROGRESS && @info "  -> $n_after tracks ($(round(compression, digits=1))x) ($(round(t, digits=2))s)"
    (combined, (step_record=record, smld_connected=smld_connected, connect_info=connect_info))
end

"""
    analyze(smld, cfg::SMLMFrameConnection.FrameConnectConfig; kwargs...) -> (combined_smld, info)

Run frame connection on localizations.
"""
analyze(smld::BasicSMLD, cfg::SMLMFrameConnection.FrameConnectConfig; kwargs...) = frameconnect_step(smld, cfg; kwargs...)

function _save_frameconnect_outputs!(dir::String, cfg::SMLMFrameConnection.FrameConnectConfig,
                             v::Int, t::Float64,
                             n_before::Int, n_after::Int,
                             smld_connected::BasicSMLD, connect_info)
    mkpath(dir)
    _save_config!(dir, cfg)
    _save_info!(dir, connect_info)

    if v >= Verbosity.STANDARD
        koff_stats = _save_frameconnect_figures(dir, smld_connected)
        _write_frameconnect_stats(dir, cfg, n_before, n_after, t, koff_stats)
    end
end

function _write_frameconnect_stats(dir, cfg, n_before, n_after, t, koff_stats=nothing)
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
