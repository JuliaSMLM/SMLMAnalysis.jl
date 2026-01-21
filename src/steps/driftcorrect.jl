"""
Drift correction step - wraps SMLMDriftCorrection.driftcorrect
"""

@kwdef struct DriftCorrectConfig <: StepConfig
    name::String = "driftcorrect"
    # SMLMDriftCorrection.driftcorrect kwargs
    degree::Int = 2
    intramodel::String = "Polynomial"  # "Polynomial" or "LegendrePoly"
    cost_fun::String = "Kdtree"        # "Kdtree" or "Entropy"
    # Acquisition type
    continuous::Bool = false  # true: TYPE 1 continuous acquisition (drift accumulates across datasets)
                              # false: TYPE 2 registered acquisition (each dataset starts near zero)
    # Diagnostics
    warn_large_intershift::Bool = true  # Warn if TYPE 2 has large inter-dataset shifts
    intershift_threshold_nm::Float64 = 500.0  # nm threshold for warning
    verbose::Int = Verbosity.STANDARD
end

function run_step!(a::Analysis, cfg::DriftCorrectConfig)
    a.smld === nothing && error("Must run Fit first")
    a.step_counter += 1
    v = _get_verbose(a, cfg)
    dir = _stepdir(a, cfg)

    v >= Verbosity.PROGRESS && @info "[$(a.step_counter)] $(cfg.name)" degree=cfg.degree model=cfg.intramodel continuous=cfg.continuous

    t = @elapsed drift_result = SMLMDriftCorrection.driftcorrect(a.smld;
        degree = cfg.degree,
        intramodel = cfg.intramodel,
        cost_fun = cfg.cost_fun,
        verbose = 0
    )

    a.smld = drift_result.smld
    a.drift_model = drift_result.model

    n_datasets = length(a.drift_model.intra)
    n_frames = a.smld.n_frames

    # Calculate inter-shift magnitudes (in nm)
    inter_shifts = _calc_inter_shifts(a.drift_model)
    max_intershift = n_datasets > 1 ? maximum(inter_shifts[2:end]) : 0.0

    # Calculate max intra-dataset drift
    max_drift = _calc_max_drift(a.drift_model, n_frames)

    # Diagnostic warnings
    if cfg.warn_large_intershift && !cfg.continuous && max_intershift > cfg.intershift_threshold_nm
        v >= Verbosity.PROGRESS && @warn "  Large inter-dataset shifts detected ($(round(max_intershift, digits=1))nm). " *
            "If data was acquired without registration, consider continuous=true"
    end
    if cfg.continuous && max_intershift > cfg.intershift_threshold_nm
        v >= Verbosity.PROGRESS && @warn "  Large inter-dataset shifts ($(round(max_intershift, digits=1))nm) " *
            "unexpected for continuous acquisition - check data alignment"
    end

    summary = Dict{Symbol,Any}(
        :max_drift_nm => round(max_drift, digits=1),
        :max_intershift_nm => round(max_intershift, digits=1),
        :n_datasets => n_datasets,
        :n_frames => n_frames,
        :continuous => cfg.continuous
    )
    _record!(a, cfg, t, summary)
    _checkpoint!(a)

    if dir !== nothing
        _save_step_outputs!(dir, a, cfg, v, t, max_drift, inter_shifts, n_frames)
    end

    v >= Verbosity.PROGRESS && @info "  → max drift $(round(max_drift, digits=1))nm, inter-shift $(round(max_intershift, digits=1))nm ($(round(t, digits=2))s)"
    a
end

"""Calculate inter-dataset shift magnitudes in nm (Euclidean distance)"""
function _calc_inter_shifts(drift_model)
    n_datasets = length(drift_model.inter)
    shifts = Float64[]
    for ds in 1:n_datasets
        dx = drift_model.inter[ds].dm[1] * 1000  # μm to nm
        dy = drift_model.inter[ds].dm[2] * 1000
        push!(shifts, sqrt(dx^2 + dy^2))
    end
    shifts
end

function _calc_max_drift(drift_model, n_frames)
    DC = SMLMDriftCorrection
    n_datasets = length(drift_model.intra)
    frames = 1:n_frames

    max_drift = 0.0
    for ds in 1:n_datasets
        drift_x = [DC.applydrift(0.0, f, drift_model.intra[ds].dm[1]) for f in frames]
        drift_y = [DC.applydrift(0.0, f, drift_model.intra[ds].dm[2]) for f in frames]
        max_drift = max(max_drift, maximum(abs.(drift_x)) * 1000, maximum(abs.(drift_y)) * 1000)
    end
    max_drift
end

function _save_step_outputs!(dir::String, a::Analysis, cfg::DriftCorrectConfig, v::Int, t::Float64,
                             max_drift::Float64, inter_shifts::Vector{Float64}, n_frames::Int)
    mkpath(dir)
    _save_config!(dir, cfg)

    if v >= Verbosity.STANDARD
        _write_drift_stats(dir, cfg, a.drift_model, t, max_drift, inter_shifts, n_frames)
        _save_drift_figures(dir, a.drift_model, n_frames, cfg.continuous)
    end

    if v >= Verbosity.DETAILED
        _save_drift_detailed(dir, a.drift_model, n_frames, inter_shifts)
    end
end

function _write_drift_stats(dir, cfg, drift_model, t, max_drift, inter_shifts, n_frames)
    n_datasets = length(drift_model.intra)
    max_intershift = n_datasets > 1 ? maximum(inter_shifts[2:end]) : 0.0

    filepath = joinpath(dir, "stats.md")
    open(filepath, "w") do io
        println(io, "# Drift Correction Statistics\n")
        println(io, "## Summary")
        println(io, "- **Mode**: $(cfg.continuous ? "Continuous (TYPE 1)" : "Registered (TYPE 2)")")
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
        println(io, "- model: $(cfg.intramodel)")
        println(io, "- cost_fun: $(cfg.cost_fun)")
        println(io, "- continuous: $(cfg.continuous)")

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

function _save_drift_figures(dir, drift_model, n_frames, continuous::Bool)
    DC = SMLMDriftCorrection
    n_datasets = length(drift_model.intra)
    frames = collect(1:n_frames)

    if n_datasets == 1
        drift_x = [DC.applydrift(0.0, f, drift_model.intra[1].dm[1]) for f in frames] .* 1000
        drift_y = [DC.applydrift(0.0, f, drift_model.intra[1].dm[2]) for f in frames] .* 1000

        fig = Figure(size=(1200, 400))

        ax1 = Axis(fig[1, 1], xlabel="Frame", ylabel="X Drift (nm)", title="X Drift")
        lines!(ax1, frames, drift_x, color=:blue)
        hlines!(ax1, [0], color=:gray, linestyle=:dash)

        ax2 = Axis(fig[1, 2], xlabel="Frame", ylabel="Y Drift (nm)", title="Y Drift")
        lines!(ax2, frames, drift_y, color=:red)
        hlines!(ax2, [0], color=:gray, linestyle=:dash)

        ax3 = Axis(fig[1, 3], xlabel="X (nm)", ylabel="Y (nm)", title="XY Path", aspect=DataAspect())
        lines!(ax3, drift_x, drift_y, color=:black)
        scatter!(ax3, [drift_x[1]], [drift_y[1]], color=:green, markersize=10)
        scatter!(ax3, [drift_x[end]], [drift_y[end]], color=:red, markersize=10)

        save(joinpath(dir, "drift_trajectory.png"), fig)
    elseif continuous
        # Continuous mode: show cumulative drift trajectory
        _save_continuous_drift_figure(dir, drift_model, n_frames)
    else
        # Registered mode: show per-dataset drift
        _save_perdataset_drift_figure(dir, drift_model, n_frames)
    end
end

"""Plot cumulative drift trajectory for continuous acquisition (TYPE 1)"""
function _save_continuous_drift_figure(dir, drift_model, n_frames)
    DC = SMLMDriftCorrection
    n_datasets = length(drift_model.intra)
    frames = collect(1:n_frames)

    # Anchor at dataset 1's inter-shift (subtract to start at origin)
    inter1_x = drift_model.inter[1].dm[1] * 1000
    inter1_y = drift_model.inter[1].dm[2] * 1000

    # Build continuous trajectory
    all_global_frames = Int[]
    all_drift_x = Float64[]
    all_drift_y = Float64[]

    for ds in 1:n_datasets
        # Total drift = inter[ds] + intra[ds](frame), anchored at origin
        inter_x = (drift_model.inter[ds].dm[1] * 1000) - inter1_x
        inter_y = (drift_model.inter[ds].dm[2] * 1000) - inter1_y

        for f in frames
            global_frame = (ds - 1) * n_frames + f
            intra_x = DC.applydrift(0.0, f, drift_model.intra[ds].dm[1]) * 1000
            intra_y = DC.applydrift(0.0, f, drift_model.intra[ds].dm[2]) * 1000

            push!(all_global_frames, global_frame)
            push!(all_drift_x, inter_x + intra_x)
            push!(all_drift_y, inter_y + intra_y)
        end
    end

    total_frames = n_datasets * n_frames

    fig = Figure(size=(1200, 600))

    ax1 = Axis(fig[1, 1], xlabel="Global Frame", ylabel="X Drift (nm)",
               title="Continuous X Drift ($(n_datasets) datasets)")
    lines!(ax1, all_global_frames, all_drift_x, color=:blue)
    hlines!(ax1, [0], color=:gray, linestyle=:dash)
    # Mark dataset boundaries
    for ds in 2:n_datasets
        vlines!(ax1, [(ds-1) * n_frames], color=:lightgray, linestyle=:dot)
    end

    ax2 = Axis(fig[1, 2], xlabel="Global Frame", ylabel="Y Drift (nm)",
               title="Continuous Y Drift")
    lines!(ax2, all_global_frames, all_drift_y, color=:red)
    hlines!(ax2, [0], color=:gray, linestyle=:dash)
    for ds in 2:n_datasets
        vlines!(ax2, [(ds-1) * n_frames], color=:lightgray, linestyle=:dot)
    end

    ax3 = Axis(fig[2, 1:2], xlabel="X (nm)", ylabel="Y (nm)",
               title="Continuous XY Trajectory", aspect=DataAspect())
    lines!(ax3, all_drift_x, all_drift_y, color=:black, linewidth=0.5)
    scatter!(ax3, [all_drift_x[1]], [all_drift_y[1]], color=:green, markersize=12, label="Start")
    scatter!(ax3, [all_drift_x[end]], [all_drift_y[end]], color=:red, markersize=12, label="End")
    axislegend(ax3, position=:lt)

    save(joinpath(dir, "drift_trajectory.png"), fig)
end

"""Plot per-dataset drift for registered acquisition (TYPE 2)"""
function _save_perdataset_drift_figure(dir, drift_model, n_frames)
    DC = SMLMDriftCorrection
    n_datasets = length(drift_model.intra)
    frames = collect(1:n_frames)
    colors = [:blue, :red, :green, :orange, :purple, :cyan, :magenta, :brown]

    fig = Figure(size=(1200, 600))

    ax1 = Axis(fig[1, 1], xlabel="Frame", ylabel="X Drift (nm)", title="X Drift per Dataset")
    ax2 = Axis(fig[1, 2], xlabel="Frame", ylabel="Y Drift (nm)", title="Y Drift per Dataset")
    ax3 = Axis(fig[2, 1:2], xlabel="X (nm)", ylabel="Y (nm)", title="XY Paths per Dataset", aspect=DataAspect())

    for ds in 1:n_datasets
        c = colors[mod1(ds, length(colors))]
        drift_x = [DC.applydrift(0.0, f, drift_model.intra[ds].dm[1]) for f in frames] .* 1000
        drift_y = [DC.applydrift(0.0, f, drift_model.intra[ds].dm[2]) for f in frames] .* 1000

        lines!(ax1, frames, drift_x, color=c, label="DS $ds")
        lines!(ax2, frames, drift_y, color=c, label="DS $ds")
        lines!(ax3, drift_x, drift_y, color=c, label="DS $ds")
    end

    if n_datasets <= 6
        axislegend(ax1, position=:lt)
    end
    save(joinpath(dir, "drift_trajectory.png"), fig)
end

function _save_drift_detailed(dir, drift_model, n_frames, inter_shifts)
    DC = SMLMDriftCorrection
    n_datasets = length(drift_model.intra)
    frames = 1:n_frames

    filepath = joinpath(dir, "per_dataset.md")
    open(filepath, "w") do io
        println(io, "# Per-Dataset Drift Details\n")
        println(io, "| Dataset | Max X (nm) | Max Y (nm) | Inter-Shift (nm) |")
        println(io, "|---------|------------|------------|------------------|")

        for ds in 1:n_datasets
            drift_x = [DC.applydrift(0.0, f, drift_model.intra[ds].dm[1]) for f in frames] .* 1000
            drift_y = [DC.applydrift(0.0, f, drift_model.intra[ds].dm[2]) for f in frames] .* 1000
            inter = inter_shifts[ds]
            println(io, "| $ds | $(round(maximum(abs.(drift_x)), digits=1)) | $(round(maximum(abs.(drift_y)), digits=1)) | $(round(inter, digits=1)) |")
        end
    end
end
