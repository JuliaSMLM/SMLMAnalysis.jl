"""
Drift correction step - wraps SMLMDriftCorrection.driftcorrect
"""

@kwdef struct DriftCorrectConfig <: SMLMData.AbstractSMLMConfig
    # SMLMDriftCorrection.driftcorrect kwargs
    degree::Int = 2
    # Acquisition type
    continuous::Bool = false  # true: TYPE 1 continuous acquisition (drift accumulates across datasets)
                              # false: TYPE 2 registered acquisition (each dataset starts near zero)
    # Chunking for continuous datasets
    n_chunks::Int = 0         # Split continuous data into N chunks (0 = no chunking)
    chunk_frames::Int = 0     # Alternative: frames per chunk (0 = use n_chunks instead)
    maxn::Int = 200           # Maximum neighbors for entropy calculation
    # Quality mode
    quality::Symbol = :singlepass  # :singlepass (fast) or :iterative (slower, converges)
    # ROI selection
    auto_roi::Bool = true  # Automatically select dense ROI for faster estimation
    # Diagnostics
    warn_large_intershift::Bool = true  # Warn if TYPE 2 has large inter-dataset shifts
    intershift_threshold_nm::Float64 = 500.0  # nm threshold for warning
end

function run_step!(a::Analysis, cfg::DriftCorrectConfig)
    a.smld === nothing && error("Must run Fit first")
    a.step_counter += 1
    v = a.verbose
    dir = _stepdir(a, cfg)

    # Determine dataset mode
    dataset_mode = cfg.continuous ? :continuous : :registered

    # Log info
    if cfg.n_chunks > 0 || cfg.chunk_frames > 0
        chunks_info = cfg.n_chunks > 0 ? "$(cfg.n_chunks) chunks" : "$(cfg.chunk_frames) frames/chunk"
        v >= Verbosity.PROGRESS && @info "[$(a.step_counter)] $(step_name(cfg))" degree=cfg.degree continuous=cfg.continuous chunks=chunks_info
    else
        v >= Verbosity.PROGRESS && @info "[$(a.step_counter)] $(step_name(cfg))" degree=cfg.degree continuous=cfg.continuous
    end

    # Tuple-pattern: returns (corrected_smld, DriftInfo) where DriftInfo contains .model
    t = @elapsed (corrected_smld, drift_info) = SMLMDriftCorrection.driftcorrect(a.smld;
        degree = cfg.degree,
        dataset_mode = dataset_mode,
        n_chunks = cfg.n_chunks,
        chunk_frames = cfg.chunk_frames,
        maxn = cfg.maxn,
        quality = cfg.quality,
        auto_roi = cfg.auto_roi,
        verbose = v >= Verbosity.DETAILED ? 1 : 0
    )

    a.smld = corrected_smld
    a.drift_model = drift_info.model

    n_datasets = a.drift_model.ndatasets
    n_frames = a.smld.n_frames

    # Calculate inter-shift magnitudes (in nm)
    inter_shifts = _calc_inter_shifts(a.drift_model)
    max_intershift = n_datasets > 1 ? maximum(inter_shifts[2:end]) : 0.0

    # Calculate max intra-dataset drift
    max_drift = _calc_max_drift(a.drift_model, n_frames; n_chunks=cfg.n_chunks)

    # Diagnostic warnings
    if cfg.warn_large_intershift && !cfg.continuous && max_intershift > cfg.intershift_threshold_nm
        v >= Verbosity.PROGRESS && @warn "  Large inter-dataset shifts detected ($(round(max_intershift, digits=1))nm). " *
            "If data was acquired without registration, consider continuous=true"
    end
    if cfg.continuous && max_intershift > cfg.intershift_threshold_nm
        v >= Verbosity.PROGRESS && @warn "  Large inter-dataset shifts ($(round(max_intershift, digits=1))nm) " *
            "unexpected for continuous acquisition - check data alignment"
    end

    # Capture convergence info from drift_info (tuple-pattern)
    converged = hasproperty(drift_info, :converged) ? drift_info.converged : nothing
    iterations = hasproperty(drift_info, :iterations) ? drift_info.iterations : nothing

    summary = Dict{Symbol,Any}(
        :max_drift_nm => round(max_drift, digits=1),
        :max_intershift_nm => round(max_intershift, digits=1),
        :n_datasets => n_datasets,
        :n_frames => n_frames,
        :continuous => cfg.continuous,
        :quality => cfg.quality,
        :converged => converged,
        :iterations => iterations
    )
    # Include drift_info in step record (tuple-pattern)
    _record!(a, cfg, t, summary; info=drift_info)
    _checkpoint!(a)

    if dir !== nothing
        _save_step_outputs!(dir, a, cfg, v, t, max_drift, inter_shifts, n_frames, converged, iterations, drift_info)
    end

    v >= Verbosity.PROGRESS && @info "  → max drift $(round(max_drift, digits=1))nm, inter-shift $(round(max_intershift, digits=1))nm ($(round(t, digits=2))s)"
    a
end

"""Calculate inter-dataset shift magnitudes in nm (Euclidean distance)"""
function _calc_inter_shifts(drift_model)
    n_datasets = drift_model.ndatasets
    shifts = Float64[]
    for ds in 1:n_datasets
        dx = drift_model.inter[ds].dm[1] * 1000  # μm to nm
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
        # Note: traj.x and traj.y are in μm, convert to nm
        max_drift = max(max_drift, maximum(abs.(traj.x)) * 1000, maximum(abs.(traj.y)) * 1000)
    end
    max_drift
end

function _save_step_outputs!(dir::String, a::Analysis, cfg::DriftCorrectConfig, v::Int, t::Float64,
                             max_drift::Float64, inter_shifts::Vector{Float64}, n_frames::Int,
                             converged::Union{Bool,Nothing}, iterations::Union{Int,Nothing}, drift_info)
    mkpath(dir)
    _save_config!(dir, cfg)
    _save_info!(dir, drift_info)

    if v >= Verbosity.STANDARD
        _write_drift_stats(dir, cfg, a.drift_model, t, max_drift, inter_shifts, n_frames, converged, iterations)
        _save_drift_figures(dir, a.drift_model, n_frames, cfg.continuous; n_chunks=cfg.n_chunks)
    end

    if v >= Verbosity.DETAILED
        _save_drift_detailed(dir, a.drift_model, n_frames, inter_shifts; n_chunks=cfg.n_chunks)
    end
end

function _write_drift_stats(dir, cfg, drift_model, t, max_drift, inter_shifts, n_frames,
                            converged::Union{Bool,Nothing}, iterations::Union{Int,Nothing})
    n_datasets = drift_model.ndatasets
    max_intershift = n_datasets > 1 ? maximum(inter_shifts[2:end]) : 0.0

    filepath = joinpath(dir, "stats.md")
    open(filepath, "w") do io
        println(io, "# Drift Correction Statistics\n")
        println(io, "## Summary")
        println(io, "- **Mode**: $(cfg.continuous ? "Continuous (TYPE 1)" : "Registered (TYPE 2)")")
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
        println(io, "- continuous: $(cfg.continuous)")
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

function _save_drift_figures(dir, drift_model, n_frames, continuous::Bool; n_chunks::Int=0)
    DC = SMLMDriftCorrection
    n_datasets = drift_model.ndatasets

    if n_datasets == 1
        # Single dataset - use drift_trajectory
        traj = DC.drift_trajectory(drift_model; dataset=1)
        drift_x = traj.x .* 1000  # μm to nm
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
        # This gives a clear picture of drift evolution across entire acquisition
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
    drift_x = (traj.x .- traj.x[1]) .* 1000  # μm to nm
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
            drift_x = traj.x .* 1000  # μm to nm
            drift_y = traj.y .* 1000
            inter = inter_shifts[ds]
            println(io, "| $ds | $(round(maximum(abs.(drift_x)), digits=1)) | $(round(maximum(abs.(drift_y)), digits=1)) | $(round(inter, digits=1)) |")
        end
    end
end
