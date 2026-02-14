"""
Analysis orchestrator: analyze() function and helpers.

The analyze() function is a fold over pipeline steps, threading state
through analyze() dispatch and collecting StepInfos. Step routing is
pure dispatch on config type — no isa checks in the loop.
"""

# ============================================================
# Log file support
# ============================================================

"""
    TeeLogger(loggers)

Logger that forwards messages to multiple loggers.
Used internally to tee @info output to both console and log file.
"""
struct TeeLogger <: Logging.AbstractLogger
    loggers::Vector{Logging.AbstractLogger}
end

Logging.min_enabled_level(tl::TeeLogger) = minimum(Logging.min_enabled_level(l) for l in tl.loggers)

function Logging.shouldlog(tl::TeeLogger, level, _module, group, id)
    any(Logging.shouldlog(l, level, _module, group, id) for l in tl.loggers)
end

function Logging.handle_message(tl::TeeLogger, args...; kwargs...)
    for l in tl.loggers
        Logging.handle_message(l, args...; kwargs...)
    end
end

Logging.catch_exceptions(::TeeLogger) = true

"""
    _with_log_file(f, outdir)

Run `f()` with @info output teed to `outdir/log.txt`.
If outdir is nothing, just runs f() directly.
"""
function _with_log_file(f, outdir)
    outdir === nothing && return f()

    logpath = joinpath(outdir, "log.txt")
    open(logpath, "a") do io
        println(io, "=== SMLMAnalysis $(Dates.now()) ===")
        file_logger = Logging.SimpleLogger(io, Logging.Info)
        tee = TeeLogger(Logging.AbstractLogger[Logging.current_logger(), file_logger])
        Logging.with_logger(tee) do
            f()
        end
        flush(io)
    end
end

# ============================================================
# Data normalization helpers
# ============================================================

"""Normalize input data to Vector{AbstractArray{<:Real,3}} for uniform processing."""
function _normalize_data(data::Vector{<:AbstractArray{<:Real,3}})
    data
end

function _normalize_data(data::AbstractArray{<:Real,3})
    [data]
end

"""Apply ROI cropping to data."""
function _apply_roi(data::Vector{<:AbstractArray{<:Real,3}}, roi)
    [crop_images(img, roi.x, roi.y) for img in data]
end

function _apply_roi(data::AbstractArray{<:Real,3}, roi)
    crop_images(data, roi.x, roi.y)
end

# ============================================================
# Step preparation (dispatch-based camera injection)
# ============================================================

"""Pre-dispatch hook: inject pipeline-level config into step configs."""
_prepare_step(cfg::SMLMData.AbstractSMLMConfig, ::SMLMData.AbstractCamera) = cfg
_prepare_step(cfg::DetectFitConfig, camera::SMLMData.AbstractCamera) = _inject_camera(cfg, camera)

# ============================================================
# Primary interface: analyze(data, config::AnalysisConfig)
# ============================================================

"""
    analyze(data, config::AnalysisConfig) -> (AnalysisResult, AnalysisInfo)

Run SMLM analysis pipeline defined by config.

Returns a tuple of (AnalysisResult, AnalysisInfo) following the JuliaSMLM tuple-pattern.

# Arguments
- `data`: Image data - one of:
  - `Vector{AbstractArray{<:Real,3}}`: Multiple datasets (primary path)
  - `AbstractArray{<:Real,3}`: Single dataset
- `config`: AnalysisConfig with camera, steps, and output settings

# Example
```julia
config = AnalysisConfig(
    camera = cam,
    steps = [
        DetectFitConfig(boxer=BoxerConfig(boxsize=9)),
        FilterConfig(photons=(500.0, Inf)),
        DriftConfig(degree=2),
        RenderConfig(zoom=20, colormap=:inferno),
    ],
    outdir = "output/",
)
(result, info) = analyze(image_stacks, config)
result.smld               # Final SMLD
info.steps[:detectfit]    # DetectFit step info
info.steps[:driftcorrect] # Drift step info
```
"""
function analyze(data, config::AnalysisConfig)
    camera = config.camera

    # Apply ROI if specified
    if config.roi !== nothing
        data = _apply_roi(data, config.roi)
        camera = crop_camera(camera, config.roi.x, config.roi.y)
    end

    _run_pipeline(_normalize_data(data), config.steps, camera, config.outdir, config.verbose)
end

"""
    analyze(data, steps::AbstractSMLMConfig...; camera, kwargs...) -> (AnalysisResult, AnalysisInfo)

Convenience varargs form. Builds AnalysisConfig from positional step configs and keyword arguments.

# Example
```julia
(result, info) = analyze(image_stacks,
    DetectFitConfig(boxer=BoxerConfig(boxsize=9)),
    FilterConfig(photons=(500.0, Inf)),
    DriftConfig(degree=2);
    camera=cam, outdir="output/")
```
"""
function analyze(data, steps::SMLMData.AbstractSMLMConfig...;
                 camera::SMLMData.AbstractCamera,
                 roi=nothing,
                 kwargs...)
    config = AnalysisConfig(steps...; camera=camera, roi=roi, kwargs...)
    analyze(data, config)
end

# ============================================================
# File-based analyze (DetectFitConfig with path/paths)
# ============================================================

"""
    analyze(config::AnalysisConfig) -> (AnalysisResult, AnalysisInfo)

Run analysis from file paths specified in DetectFitConfig.
No image data argument needed - data is loaded from files.

# Example
```julia
config = AnalysisConfig(
    camera = cam,
    steps = [
        DetectFitConfig(path="data.h5"),
        FilterConfig(photons=(500.0, Inf)),
    ],
    outdir = "output/",
)
(result, info) = analyze(config)
```
"""
function analyze(config::AnalysisConfig)
    _run_pipeline(nothing, config.steps, config.camera, config.outdir, config.verbose)
end

# ============================================================
# Pipeline loop (pure dispatch, no isa routing)
# ============================================================

"""
    _run_pipeline(initial_state, steps, camera, outdir, verbose)

Internal pipeline executor. Folds `analyze()` dispatch over the steps vector,
threading state through each step. Step routing is determined entirely by
Julia's method dispatch on `(state_type, config_type)`.

Wrong step ordering produces a MethodError — e.g., FilterConfig before
DetectFitConfig gives `no method matching analyze(::Vector{...}, ::FilterConfig)`.
"""
function _run_pipeline(initial_state, steps::Vector{SMLMData.AbstractSMLMConfig},
                       camera::SMLMData.AbstractCamera,
                       outdir::Union{String,Nothing}, v::Int)
    t_start = time_ns()
    outdir !== nothing && mkpath(outdir)

    step_infos = StepInfo[]
    state = initial_state
    last_smld = nothing

    _with_log_file(outdir) do
        for (i, cfg) in enumerate(steps)
            cfg = _prepare_step(cfg, camera)

            (state, step_info) = analyze(state, cfg;
                outdir=outdir, step_number=i, verbose=v)

            push!(step_infos, step_info)

            # Track last SMLD for result construction
            if state isa SMLMData.BasicSMLD
                last_smld = state
            end
        end
    end

    # Build result
    last_smld === nothing && error("Pipeline produced no SMLD. Did you include a DetectFitConfig step?")

    smld_connected = nothing
    drift_model = nothing
    for si in step_infos
        if si.info isa SMLMFrameConnection.FrameConnectInfo
            smld_connected = si.info.connected
        end
        if si.config isa SMLMDriftCorrection.DriftConfig && si.info !== nothing
            drift_model = si.info.model
        end
    end
    result = AnalysisResult(last_smld, smld_connected, drift_model)

    elapsed_s = (time_ns() - t_start) / 1e9
    steps_dict = Dict{Symbol, Any}()
    for si in step_infos
        if si.info !== nothing
            steps_dict[Symbol(si.name)] = si.info
        end
    end
    info = AnalysisInfo(elapsed_s, steps_dict, step_infos)

    if outdir !== nothing
        _write_summary(outdir, step_infos, last_smld)
    end

    v >= Verbosity.PROGRESS && @info "Pipeline complete: $(length(last_smld.emitters)) localizations ($(round(elapsed_s, digits=2))s)"

    (result, info)
end

# ============================================================
# Pipeline summary
# ============================================================

function _write_summary(outdir::String, step_infos::Vector{StepInfo}, smld)
    filepath = joinpath(outdir, "summary.md")
    open(filepath, "w") do io
        println(io, "# Analysis Summary\n")
        println(io, "Generated: $(Dates.now())")
        println(io, "")

        total_time = sum(s.elapsed_s for s in step_infos)
        println(io, "## Pipeline")
        println(io, "")
        println(io, "| Step | Name | Time | Result |")
        println(io, "|------|------|------|--------|")

        for s in step_infos
            result_str = if haskey(s.summary, :n_rois)
                "$(s.summary[:n_rois]) ROIs"
            elseif haskey(s.summary, :n_fits)
                "$(s.summary[:n_fits]) fits"
            elseif haskey(s.summary, :n_after) && haskey(s.summary, :n_before)
                "$(s.summary[:n_after])/$(s.summary[:n_before])"
            elseif haskey(s.summary, :compression)
                "$(s.summary[:compression])x compression"
            elseif haskey(s.summary, :max_drift_nm)
                "$(s.summary[:max_drift_nm])nm max drift"
            elseif haskey(s.summary, :strategy)
                "$(s.summary[:strategy])"
            else
                ""
            end
            println(io, "| $(s.number) | $(s.name) | $(round(s.elapsed_s, digits=2))s | $result_str |")
        end

        println(io, "")
        println(io, "**Total time**: $(round(total_time, digits=2))s")

        if smld !== nothing
            println(io, "")
            println(io, "**Final**: $(length(smld.emitters)) localizations")
        end
    end
end
