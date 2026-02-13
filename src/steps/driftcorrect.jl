"""
Drift correction step - uses SMLMDriftCorrection.DriftConfig directly.

No wrapper config -- the upstream DriftConfig is used as a pipeline step.
Intershift warnings are always checked and emitted at PROGRESS verbosity.
"""

# Alias for convenience (like RenderConfig pattern)
const DriftConfig = SMLMDriftCorrection.DriftConfig

# Override step_name to keep backward-compatible output directory naming
step_name(::SMLMDriftCorrection.DriftConfig) = "driftcorrect"

"""
    driftcorrect_step(smld, cfg; outdir, step_number, verbose) -> (corrected_smld, DriftInfo)

Run drift correction on `smld`, returning the corrected SMLD as the primary
result, plus the upstream DriftInfo.

# Returns
`(corrected_smld, DriftInfo)`
"""
function driftcorrect_step(smld::BasicSMLD, cfg::SMLMDriftCorrection.DriftConfig;
                           outdir::Union{String,Nothing}=nothing,
                           step_number::Int=0,
                           verbose::Int=Verbosity.STANDARD)
    v = verbose
    dir = step_outdir(outdir, step_number, cfg)

    # Log info
    if cfg.n_chunks > 0 || cfg.chunk_frames > 0
        chunks_info = cfg.n_chunks > 0 ? "$(cfg.n_chunks) chunks" : "$(cfg.chunk_frames) frames/chunk"
        v >= Verbosity.PROGRESS && @info "[$step_number] $(step_name(cfg))" degree=cfg.degree dataset_mode=cfg.dataset_mode chunks=chunks_info
    else
        v >= Verbosity.PROGRESS && @info "[$step_number] $(step_name(cfg))" degree=cfg.degree dataset_mode=cfg.dataset_mode
    end

    # Config dispatch with verbose injected from pipeline verbosity
    drift_cfg = SMLMDriftCorrection.DriftConfig(;
        (f => getfield(cfg, f) for f in fieldnames(SMLMDriftCorrection.DriftConfig) if f != :verbose)...,
        verbose = v >= Verbosity.DETAILED ? 1 : 0
    )
    t = @elapsed (corrected_smld, drift_info) = SMLMDriftCorrection.driftcorrect(smld, drift_cfg)

    drift_model = drift_info.model
    n_datasets_val = drift_model.ndatasets
    n_frames = smld.n_frames

    # Calculate inter-shift magnitudes (in nm)
    inter_shifts = _calc_inter_shifts(drift_model)
    max_intershift = n_datasets_val > 1 ? maximum(inter_shifts[2:end]) : 0.0

    # Calculate max intra-dataset drift
    max_drift = _calc_max_drift(drift_model, n_frames; n_chunks=cfg.n_chunks)

    # Diagnostic warnings (always check, warn at PROGRESS level)
    intershift_threshold_nm = 500.0
    if cfg.dataset_mode == :registered && max_intershift > intershift_threshold_nm
        v >= Verbosity.PROGRESS && @warn "  Large inter-dataset shifts detected ($(round(max_intershift, digits=1))nm). " *
            "If data was acquired without registration, consider dataset_mode=:continuous"
    end
    if cfg.dataset_mode == :continuous && max_intershift > intershift_threshold_nm
        v >= Verbosity.PROGRESS && @warn "  Large inter-dataset shifts ($(round(max_intershift, digits=1))nm) " *
            "unexpected for continuous acquisition - check data alignment"
    end

    if dir !== nothing
        converged = hasproperty(drift_info, :converged) ? drift_info.converged : nothing
        iterations = hasproperty(drift_info, :iterations) ? drift_info.iterations : nothing
        _save_driftcorrect_outputs!(dir, drift_model, cfg, v, t, max_drift, inter_shifts, n_frames, converged, iterations, drift_info)
    end

    v >= Verbosity.PROGRESS && @info "  -> max drift $(round(max_drift, digits=1))nm, inter-shift $(round(max_intershift, digits=1))nm ($(round(t, digits=2))s)"
    (corrected_smld, drift_info)
end

"""
    analyze(smld, cfg::SMLMDriftCorrection.DriftConfig; kwargs...) -> (corrected_smld, StepInfo)

Run drift correction on localizations.
"""
function analyze(smld::BasicSMLD, cfg::SMLMDriftCorrection.DriftConfig;
                 outdir=nothing, step_number::Int=0, verbose::Int=Verbosity.STANDARD)
    n_frames = smld.n_frames
    t = @elapsed (corrected, drift_info) = driftcorrect_step(smld, cfg;
        outdir=outdir, step_number=step_number, verbose=verbose)

    # Build summary (needs smld.n_frames for max_drift calculation)
    drift_model = drift_info.model
    n_datasets_val = drift_model.ndatasets
    inter_shifts = _calc_inter_shifts(drift_model)
    max_intershift = n_datasets_val > 1 ? maximum(inter_shifts[2:end]) : 0.0
    max_drift = _calc_max_drift(drift_model, n_frames; n_chunks=cfg.n_chunks)
    converged = hasproperty(drift_info, :converged) ? drift_info.converged : nothing
    iterations = hasproperty(drift_info, :iterations) ? drift_info.iterations : nothing

    summary = Dict{Symbol,Any}(
        :max_drift_nm => round(max_drift, digits=1),
        :max_intershift_nm => round(max_intershift, digits=1),
        :n_datasets => n_datasets_val,
        :n_frames => n_frames,
        :dataset_mode => cfg.dataset_mode,
        :quality => cfg.quality,
        :converged => converged,
        :iterations => iterations
    )

    (corrected, StepInfo(step_number, cfg, t, summary; info=drift_info))
end

"""Calculate inter-dataset shift magnitudes in nm (Euclidean distance)"""
function _calc_inter_shifts(drift_model)
    n_datasets = drift_model.ndatasets
    shifts = Float64[]
    for ds in 1:n_datasets
        dx = drift_model.inter[ds].dm[1] * 1000  # um to nm
        dy = drift_model.inter[ds].dm[2] * 1000
        push!(shifts, sqrt(dx^2 + dy^2))
    end
    shifts
end

"""Calculate max intra-dataset drift using drift_trajectory"""
function _calc_max_drift(drift_model, n_frames; n_chunks::Int=0)
    DC = SMLMDriftCorrection
    n_datasets = drift_model.ndatasets

    # Use drift_trajectory to get all drift values
    max_drift = 0.0
    for ds in 1:n_datasets
        traj = DC.drift_trajectory(drift_model; dataset=ds)
        # Note: traj.x and traj.y are in um, convert to nm
        max_drift = max(max_drift, maximum(abs.(traj.x)) * 1000, maximum(abs.(traj.y)) * 1000)
    end
    max_drift
end

function _save_driftcorrect_outputs!(dir::String, drift_model, cfg::SMLMDriftCorrection.DriftConfig, v::Int, t::Float64,
                             max_drift::Float64, inter_shifts::Vector{Float64}, n_frames::Int,
                             converged::Union{Bool,Nothing}, iterations::Union{Int,Nothing}, drift_info)
    mkpath(dir)
    _save_config!(dir, cfg)
    _save_info!(dir, drift_info)

    if v >= Verbosity.STANDARD
        _write_drift_stats(dir, cfg, drift_model, t, max_drift, inter_shifts, n_frames, converged, iterations)
        _save_drift_figures(dir, drift_model, n_frames, cfg.dataset_mode; n_chunks=cfg.n_chunks)
    end

    if v >= Verbosity.DETAILED
        _save_drift_detailed(dir, drift_model, n_frames, inter_shifts; n_chunks=cfg.n_chunks)
    end
end

function _write_drift_stats(dir, cfg, drift_model, t, max_drift, inter_shifts, n_frames,
                            converged::Union{Bool,Nothing}, iterations::Union{Int,Nothing})
    n_datasets = drift_model.ndatasets
    max_intershift = n_datasets > 1 ? maximum(inter_shifts[2:end]) : 0.0

    mode_str = cfg.dataset_mode == :continuous ? "Continuous" : "Registered"

    filepath = joinpath(dir, "stats.md")
    open(filepath, "w") do io
        println(io, "# Drift Correction Statistics\n")
        println(io, "## Summary")
        println(io, "- **Mode**: $mode_str ($(cfg.dataset_mode))")
        println(io, "- **Quality**: $(cfg.quality)")
        if cfg.quality == :iterative && converged !== nothing
            println(io, "- **Converged**: $(converged)")
            println(io, "- **Iterations**: $(iterations)")
        end
        println(io, "- **Max intra-dataset drift**: $(round(max_drift, digits=1)) nm")
        if n_datasets > 1
            println(io, "- **Max inter-dataset shift**: $(round(max_intershift, digits=1)) nm")
        end
        println(io, "- **Datasets**: $n_datasets")
        println(io, "- **Frames per dataset**: $n_frames")
        println(io, "- **Time**: $(round(t, digits=2))s")
        println(io, "")
        println(io, "## Parameters")
        println(io, "- degree: $(cfg.degree)")
        println(io, "- dataset_mode: $(cfg.dataset_mode)")
        println(io, "- quality: $(cfg.quality)")
        if cfg.n_chunks > 0
            println(io, "- n_chunks: $(cfg.n_chunks)")
        end

        if n_datasets > 1
            println(io, "")
            println(io, "## Inter-Dataset Shifts")
            println(io, "| Dataset | Shift (nm) |")
            println(io, "|---------|------------|")
            for (ds, shift) in enumerate(inter_shifts)
                println(io, "| $ds | $(round(shift, digits=1)) |")
            end
        end
    end
end

function _save_drift_figures(dir, drift_model, n_frames, dataset_mode::Symbol; n_chunks::Int=0)
    DC = SMLMDriftCorrection
    n_datasets = drift_model.ndatasets
    continuous = dataset_mode == :continuous

    if n_datasets == 1
        # Single dataset - use drift_trajectory
        traj = DC.drift_trajectory(drift_model; dataset=1)
        drift_x = traj.x .* 1000  # um to nm
        drift_y = traj.y .* 1000

        fig = Figure(size=(1200, 400))

        ax1 = Axis(fig[1, 1], xlabel="Frame", ylabel="X Drift (nm)", title="X Drift")
        lines!(ax1, traj.frames, drift_x, color=:blue)
        hlines!(ax1, [0], color=:gray, linestyle=:dash)

        ax2 = Axis(fig[1, 2], xlabel="Frame", ylabel="Y Drift (nm)", title="Y Drift")
        lines!(ax2, traj.frames, drift_y, color=:red)
        hlines!(ax2, [0], color=:gray, linestyle=:dash)

        ax3 = Axis(fig[1, 3], xlabel="X (nm)", ylabel="Y (nm)", title="XY Path", aspect=DataAspect())
        lines!(ax3, drift_x, drift_y, color=:black)
        scatter!(ax3, [drift_x[1]], [drift_y[1]], color=:green, markersize=10)
        scatter!(ax3, [drift_x[end]], [drift_y[end]], color=:red, markersize=10)

        save(joinpath(dir, "drift_trajectory.png"), fig)
    else
        # Multi-dataset: always show global frame continuous view
        _save_continuous_drift_figure(dir, drift_model, n_frames; n_chunks=n_chunks, continuous=continuous)
    end
end

"""Plot drift trajectory with global frame view for multi-dataset acquisitions.
Uses drift_trajectory default mode (inter + intra per chunk) so boundary
discontinuities between chunks are visible as gaps."""
function _save_continuous_drift_figure(dir, drift_model, n_frames; n_chunks::Int=0, continuous::Bool=true)
    DC = SMLMDriftCorrection
    n_datasets = drift_model.ndatasets

    # Default mode: inter[ds] + intra[ds] per chunk, preserving boundary gaps
    traj = DC.drift_trajectory(drift_model)

    # Anchor at origin (subtract first point)
    drift_x = (traj.x .- traj.x[1]) .* 1000  # um to nm
    drift_y = (traj.y .- traj.y[1]) .* 1000

    mode_str = continuous ? "Continuous" : "Registered"
    label = n_chunks > 0 ? "Chunk" : "Dataset"
    fig = Figure(size=(1200, 600))

    total_frames = maximum(traj.frames)

    colors = [:blue, :red, :green, :orange, :purple, :cyan, :magenta, :brown,
              :darkblue, :darkred, :darkgreen, :darkorange, :violet, :teal, :pink, :chocolate,
              :navy, :crimson, :forestgreen, :coral]

    # Top row: X and Y drift vs frame, per-chunk segments
    top = fig[1, 1:2] = GridLayout()
    ax1 = Axis(top[1, 1], xlabel="Global Frame", ylabel="X Drift (nm)",
               title="X Drift ($(n_datasets) $(lowercase(label))s, $mode_str)")
    hlines!(ax1, [0], color=:gray, linestyle=:dash)

    ax2 = Axis(top[1, 2], xlabel="Global Frame", ylabel="Y Drift (nm)",
               title="Y Drift")
    hlines!(ax2, [0], color=:gray, linestyle=:dash)

    for ds in 1:n_datasets
        mask = traj.dataset .== ds
        c = colors[mod1(ds, length(colors))]
        lines!(ax1, traj.frames[mask], drift_x[mask], color=c)
        lines!(ax2, traj.frames[mask], drift_y[mask], color=c)
    end
    xlims!(ax1, 1, total_frames)
    xlims!(ax2, 1, total_frames)

    # Bottom row: XY trajectory per-chunk
    bottom = fig[2, 1:2] = GridLayout()
    ax3 = Axis(bottom[1, 1], xlabel="X (nm)", ylabel="Y (nm)",
               title="XY Drift Trajectory", aspect=DataAspect())

    for ds in 1:n_datasets
        mask = traj.dataset .== ds
        c = colors[mod1(ds, length(colors))]
        lines!(ax3, drift_x[mask], drift_y[mask], color=c, linewidth=1.0, label="$label $ds")
    end

    # Mark start and end
    first_idx = argmin(traj.frames)
    last_idx = argmax(traj.frames)
    scatter!(ax3, [drift_x[first_idx]], [drift_y[first_idx]], color=:green, markersize=12, label="Start")
    scatter!(ax3, [drift_x[last_idx]], [drift_y[last_idx]], color=:red, markersize=12, label="End")

    if n_datasets <= 8
        Legend(bottom[1, 2], ax3, framevisible=false)
    end

    save(joinpath(dir, "drift_trajectory.png"), fig)
end

function _save_drift_detailed(dir, drift_model, n_frames, inter_shifts; n_chunks::Int=0)
    DC = SMLMDriftCorrection
    n_datasets = drift_model.ndatasets

    filepath = joinpath(dir, "per_dataset.md")
    open(filepath, "w") do io
        label = n_chunks > 0 ? "Chunk" : "Dataset"
        println(io, "# Per-$label Drift Details\n")
        println(io, "| $label | Max X (nm) | Max Y (nm) | Inter-Shift (nm) |")
        println(io, "|---------|------------|------------|------------------|")

        for ds in 1:n_datasets
            traj = DC.drift_trajectory(drift_model; dataset=ds)
            drift_x = traj.x .* 1000  # um to nm
            drift_y = traj.y .* 1000
            inter = inter_shifts[ds]
            println(io, "| $ds | $(round(maximum(abs.(drift_x)), digits=1)) | $(round(maximum(abs.(drift_y)), digits=1)) | $(round(inter, digits=1)) |")
        end
    end
end
