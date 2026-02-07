"""
Core analysis functions: reset!, checkpoint!, run_recipe, helpers
"""

# ============================================================
# Checkpointing
# ============================================================

"""
    _checkpoint!(a::Analysis; save::Union{Bool,Nothing}=nothing)

Save checkpoint to memory, and optionally to disk.

The `save` argument controls disk persistence:
- `nothing`: use `a.checkpoint` setting (default)
- `true`: force save to disk
- `false`: skip disk save (memory only)
"""
function _checkpoint!(a::Analysis; save::Union{Bool,Nothing}=nothing)
    # Always save to memory
    a.checkpoints[a.step_counter] = AnalysisCheckpoint(
        a.roi_batch,
        a.roi_datasets,
        a.smld_raw,
        a.smld,
        a.smld_connected,
        a.drift_model,
        a.bagol_result,
        a.bagol_smld
    )

    # Determine if we should save to disk
    should_save = save === nothing ? a.checkpoint : save

    if should_save && a.outdir !== nothing
        path = _save_checkpoint!(a)
        if path !== nothing && a.verbose >= Verbosity.PROGRESS
            @info "Checkpoint saved to disk: $(basename(path))"
        end
    end
end

function checkpoint!(a::Analysis)
    _checkpoint!(a)
    @info "Checkpoint saved at step $(a.step_counter)"
    a
end

# ============================================================
# Reset
# ============================================================

function reset!(a::Analysis)
    a.roi_batch = nothing
    a.roi_datasets = nothing
    a.smld_raw = nothing
    a.smld = nothing
    a.smld_connected = nothing
    a.drift_model = nothing
    a.bagol_result = nothing
    a.bagol_smld = nothing
    empty!(a.checkpoints)
    empty!(a.steps)
    a.step_counter = 0
    @info "Reset to initial state"
    a
end

function reset!(a::Analysis, step::Int)
    step < 0 && error("Step must be >= 0 (0 = initial state)")
    step == 0 && return reset!(a)

    # Try memory first, then disk
    if haskey(a.checkpoints, step)
        cp = a.checkpoints[step]
        a.roi_batch = cp.roi_batch
        a.roi_datasets = cp.roi_datasets
        a.smld_raw = cp.smld_raw
        a.smld = cp.smld
        a.smld_connected = cp.smld_connected
        a.drift_model = cp.drift_model
        a.bagol_result = cp.bagol_result
        a.bagol_smld = cp.bagol_smld

        # Trim history to step
        a.steps = a.steps[1:step]
        a.step_counter = step

        # Remove later checkpoints (memory)
        for k in collect(keys(a.checkpoints))
            k > step && delete!(a.checkpoints, k)
        end

        @info "Reset to after step $step: $(a.steps[end].name)"
    else
        # Try loading from disk
        path = _find_checkpoint(a, step)
        if path === nothing
            available_mem = sort(collect(keys(a.checkpoints)))
            available_disk = Int[]
            if a.outdir !== nothing
                available_disk = first.(_list_checkpoints(a.outdir))
            end
            error("No checkpoint at step $step. Memory: $available_mem, Disk: $available_disk")
        end

        # Load from disk
        data = _load_checkpoint(path)
        a.smld_raw = data.smld_raw
        a.smld = data.smld
        a.smld_connected = data.smld_connected
        a.drift_model = data.drift_model
        a.bagol_result = data.bagol_result
        a.bagol_smld = data.bagol_smld
        a.roi_batch = nothing  # Not saved, would need to regenerate
        a.roi_datasets = data.roi_datasets

        # Restore history up to this step
        a.steps = data.steps[1:min(step, length(data.steps))]
        a.step_counter = step

        # Clear memory checkpoints beyond this step
        for k in collect(keys(a.checkpoints))
            k > step && delete!(a.checkpoints, k)
        end

        # Store loaded checkpoint in memory for future resets
        a.checkpoints[step] = AnalysisCheckpoint(
            nothing,  # roi_batch not saved
            data.roi_datasets,
            data.smld_raw,
            data.smld,
            data.smld_connected,
            data.drift_model,
            get(data, :bagol_result, nothing),
            get(data, :bagol_smld, nothing)
        )

        @info "Reset to after step $step (loaded from disk): $(a.steps[end].name)"
    end
    a
end

# ============================================================
# Debug helper
# ============================================================

function debug!(a::Analysis, cfg::SMLMData.AbstractSMLMConfig)
    # Temporarily set analysis to DEBUG verbosity, run step, restore
    old_verbose = a.verbose
    a.verbose = Verbosity.DEBUG
    try
        run_step!(a, cfg)
    finally
        a.verbose = old_verbose
    end
end

# ============================================================
# Helpers used by all steps
# ============================================================

function _stepdir(a::Analysis, cfg::SMLMData.AbstractSMLMConfig)
    a.outdir === nothing && return nothing
    joinpath(a.outdir, "$(lpad(a.step_counter, 2, '0'))_$(step_name(cfg))")
end

function _record!(a::Analysis, cfg::SMLMData.AbstractSMLMConfig, t::Float64, summary::Dict{Symbol,Any}; info=nothing)
    push!(a.steps, StepRecord(a.step_counter, cfg, t, summary; info=info))
end

function _save_config!(dir::String, cfg::SMLMData.AbstractSMLMConfig)
    filepath = joinpath(dir, "config.toml")
    open(filepath, "w") do io
        println(io, "# $(nameof(typeof(cfg)))")
        println(io, "type = \"$(nameof(typeof(cfg)))\"")
        for f in fieldnames(typeof(cfg))
            v = getfield(cfg, f)
            if v isa String
                println(io, "$f = \"$v\"")
            elseif v isa Symbol
                println(io, "$f = \"$v\"")
            elseif v !== nothing
                println(io, "$f = $v")
            end
        end
    end
end

"""
    _save_info!(dir::String, info; section::String="")

Write upstream Info struct fields to `info.toml` in TOML format.

Writes scalar fields (numbers, bools, strings, symbols, tuples of scalars).
Skips complex fields (arrays, dicts, structs like BasicSMLD, models).

When `section` is empty, writes a fresh file with type header.
When `section` is provided, appends a `[section]` block.
"""
function _save_info!(dir::String, info; section::String="")
    filepath = joinpath(dir, "info.toml")
    open(filepath, section == "" ? "w" : "a") do io
        if section == ""
            println(io, "# Upstream package info")
            println(io, "type = \"$(nameof(typeof(info)))\"")
        else
            println(io, "\n[$section]")
        end
        for f in fieldnames(typeof(info))
            v = getfield(info, f)
            _write_info_field!(io, f, v)
        end
    end
end

"""Write a single field to info.toml, skipping complex types."""
function _write_info_field!(io::IO, name::Symbol, v::Number)
    println(io, "$name = $v")
end
function _write_info_field!(io::IO, name::Symbol, v::Bool)
    println(io, "$name = $v")
end
function _write_info_field!(io::IO, name::Symbol, v::String)
    println(io, "$name = \"$v\"")
end
function _write_info_field!(io::IO, name::Symbol, v::Symbol)
    println(io, "$name = \"$v\"")
end
function _write_info_field!(io::IO, name::Symbol, v::Nothing)
    println(io, "$name = \"nothing\"")
end
function _write_info_field!(io::IO, name::Symbol, v::Tuple)
    # Only write tuples of scalars
    if all(x -> x isa Union{Number, Bool, String, Symbol}, v)
        vals = join([x isa String || x isa Symbol ? "\"$x\"" : "$x" for x in v], ", ")
        println(io, "$name = [$vals]")
    end
    # Skip tuples containing complex types
end
function _write_info_field!(io::IO, ::Symbol, ::Any)
    # Skip: AbstractVector, AbstractArray, AbstractDict, complex structs
end

function _write_summary(a::Analysis)
    a.outdir === nothing && return

    filepath = joinpath(a.outdir, "summary.md")
    open(filepath, "w") do io
        println(io, "# Analysis Summary\n")
        println(io, "Generated: $(Dates.now())")
        println(io, "")

        total_time = sum(s.timing for s in a.steps)
        println(io, "## Pipeline")
        println(io, "")
        println(io, "| Step | Name | Time | Result |")
        println(io, "|------|------|------|--------|")

        for s in a.steps
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
            println(io, "| $(s.number) | $(s.name) | $(round(s.timing, digits=2))s | $result_str |")
        end

        println(io, "")
        println(io, "**Total time**: $(round(total_time, digits=2))s")

        if a.smld !== nothing
            println(io, "")
            println(io, "**Final**: $(length(a.smld.emitters)) localizations")
        end
    end
end

# ============================================================
# Convenience: one-shot analyze
# ============================================================

# ============================================================
# Primary interface: analyze(data, config::AnalysisConfig)
# ============================================================

"""
    analyze(data, config::AnalysisConfig) -> (Analysis, AnalysisInfo)

Run SMLM analysis pipeline defined by config.

Returns a tuple of (Analysis, AnalysisInfo) following the JuliaSMLM tuple-pattern.

# Arguments
- `data`: Image stack (H×W×N array) or path to data file
- `config`: AnalysisConfig with camera, steps, and output settings

# Example
```julia
config = AnalysisConfig(
    camera = cam,
    steps = [
        DetectFitConfig(boxsize=9),
        FilterConfig(photons=(500.0, Inf)),
        DriftCorrectConfig(degree=2),
        RenderConfig(zoom=20, colormap=:inferno),
    ],
    outdir = "output/",
)
(result, info) = analyze(images, config)
result.smld               # Final SMLD
info.steps[:detectfit]    # DetectFit step info
info.steps[:driftcorrect] # Drift step info
```
"""
function analyze(data, config::AnalysisConfig)
    t_start = time_ns()

    # Create analysis object
    a = Analysis(data, config.camera;
        roi=config.roi, outdir=config.outdir, verbose=config.verbose, checkpoint=config.checkpoint)

    # Execute pipeline
    _execute_pipeline!(a, config.steps)

    # Write summary if output directory specified
    if config.outdir !== nothing
        _write_summary(a)
    end

    # Build AnalysisInfo from step records (tuple-pattern)
    elapsed_s = (time_ns() - t_start) / 1e9
    info = _build_analysis_info(a, elapsed_s)

    (a, info)
end

"""
    analyze(data, steps::AbstractSMLMConfig...; camera, kwargs...) -> (Analysis, AnalysisInfo)

Convenience varargs form. Builds AnalysisConfig from positional step configs and keyword arguments.

# Example
```julia
(result, info) = analyze(images,
    DetectFitConfig(boxsize=9),
    FilterConfig(photons=(500.0, Inf)),
    DriftCorrectConfig(degree=2);
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

"""
    _execute_pipeline!(a::Analysis, steps)

Execute pipeline steps with automatic fusion of adjacent boxer+fitter configs.
"""
function _execute_pipeline!(a::Analysis, steps::Vector{SMLMData.AbstractSMLMConfig})
    i = 1
    while i <= length(steps)
        # TODO: Detect adjacent BoxerConfig + FitConfig for fusion when upstream configs are used directly
        # For now, execute each step sequentially
        run_step!(a, steps[i])
        i += 1
    end
end

"""
    get_config(a::Analysis) -> AnalysisConfig

Extract an AnalysisConfig from a completed Analysis, enabling reproducibility.

After interactive exploration with run_step!, extract the config that produced the result:

```julia
a = Analysis(images, camera; outdir="output/")
run_step!(a, DetectFitConfig(boxsize=9))
run_step!(a, FilterConfig(photons=(500.0, Inf)))

config = get_config(a)  # Reproducible config from step history
(result, info) = analyze(new_images, config)
```
"""
function get_config(a::Analysis)
    AnalysisConfig(
        camera = a.camera,
        steps = SMLMData.AbstractSMLMConfig[s.config for s in a.steps],
        roi = nothing,  # ROI consumed at construction (camera/images already cropped)
        outdir = a.outdir,
        verbose = a.verbose,
        checkpoint = a.checkpoint
    )
end

# ============================================================
# Legacy kwargs interface (backward compatible)
# ============================================================

"""
    analyze(data, camera::AbstractCamera; kwargs...) -> (Analysis, AnalysisInfo)

Legacy convenience interface with flat keyword arguments.
Builds a default pipeline from kwargs and runs it.

See `analyze(data, config::AnalysisConfig)` for the primary interface.
"""
function analyze(data, camera::SMLMData.AbstractCamera;
                 outdir=nothing,
                 verbose=Verbosity.STANDARD,
                 n_datasets::Int=1,
                 # Detection
                 boxsize=11, detect_min_photons=500.0, psf_sigma=0.135, backend=:auto,
                 # Fitting
                 psf_model=:variable, iterations=20,
                 # Filtering
                 filter=true, min_photons=500.0, max_precision=0.007, min_pvalue=1e-3,
                 # Frame connection
                 frameconnect=false, maxframegap=5,
                 # Drift
                 drift=true, degree=2,
                 # Isolated
                 isolated=false, n_sigma=2.0,
                 # Render
                 render=true, render_zoom=20)

    # Build steps from kwargs
    steps = SMLMData.AbstractSMLMConfig[]

    push!(steps, DetectFitConfig(
        boxsize=boxsize,
        min_photons=detect_min_photons,
        psf_sigma=psf_sigma,
        backend=backend,
        psf_model=psf_model,
        iterations=iterations,
        n_datasets=n_datasets,
        filter_min_photons=min_photons,
        filter_max_precision=max_precision,
        filter_min_pvalue=min_pvalue
    ))

    if filter
        push!(steps, FilterConfig(
            photons=(min_photons, Inf),
            precision=(0.0, max_precision),
            pvalue=(min_pvalue, 1.0)
        ))
    end

    if frameconnect
        push!(steps, FrameConnectConfig(maxframegap=maxframegap))
    end

    if drift
        push!(steps, DriftCorrectConfig(degree=degree))
    end

    if isolated
        push!(steps, IsolatedConfig(n_sigma=n_sigma))
    end

    if render
        push!(steps, SMLMRender.RenderConfig(zoom=render_zoom, colormap=:inferno))
    end

    config = AnalysisConfig(
        camera=camera, steps=steps,
        outdir=outdir, verbose=verbose
    )
    analyze(data, config)
end

"""
    _build_analysis_info(a::Analysis, elapsed_s::Float64) -> AnalysisInfo

Build AnalysisInfo from step records, aggregating per-step info structs.
"""
function _build_analysis_info(a::Analysis, elapsed_s::Float64)
    steps = Dict{Symbol, Any}()
    for step in a.steps
        step_name_sym = Symbol(step.name)
        if step.info !== nothing
            steps[step_name_sym] = step.info
        end
    end
    AnalysisInfo(elapsed_s, steps)
end

"""
    get_analysis_info(a::Analysis) -> AnalysisInfo

Extract AnalysisInfo from an Analysis object.

Useful when running steps interactively with run_step! and wanting
to get the aggregated info at the end.
"""
function get_analysis_info(a::Analysis)
    # Sum up timing from all steps
    elapsed_s = sum(s.timing for s in a.steps; init=0.0)
    _build_analysis_info(a, elapsed_s)
end
