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
        a.drift_model
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
            data.drift_model
        )

        @info "Reset to after step $step (loaded from disk): $(a.steps[end].name)"
    end
    a
end

# ============================================================
# Recipe execution
# ============================================================

function run_recipe(recipe::Vector{<:StepConfig}, data, camera::SMLMData.AbstractCamera;
                    outdir=nothing, verbose=Verbosity.STANDARD, checkpoint=false)
    a = Analysis(data, camera; outdir, verbose, checkpoint)

    if outdir !== nothing
        mkpath(outdir)
        _save_recipe(outdir, recipe)
    end

    for cfg in recipe
        run_step!(a, cfg)
    end

    if outdir !== nothing
        _write_summary(a)
    end

    a
end

function _save_recipe(outdir, recipe)
    filepath = joinpath(outdir, "recipe.toml")
    open(filepath, "w") do io
        println(io, "# Analysis Recipe")
        println(io, "# Generated: $(Dates.now())")
        println(io, "")

        for (i, cfg) in enumerate(recipe)
            println(io, "[[step]]")
            println(io, "number = $i")
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
            println(io, "")
        end
    end
end

# ============================================================
# Debug helper
# ============================================================

function debug!(a::Analysis, cfg::StepConfig)
    # Create modified config with DEBUG verbosity
    T = typeof(cfg)
    fields = fieldnames(T)
    vals = Dict{Symbol,Any}()
    for f in fields
        if f == :verbose
            vals[f] = Verbosity.DEBUG
        else
            vals[f] = getfield(cfg, f)
        end
    end
    new_cfg = T(; vals...)
    run_step!(a, new_cfg)
end

# ============================================================
# Helpers used by all steps
# ============================================================

function _get_verbose(a::Analysis, cfg::StepConfig)
    # Use step's verbose if not STANDARD, else use analysis default
    cfg.verbose != Verbosity.STANDARD ? cfg.verbose : a.verbose
end

function _stepdir(a::Analysis, cfg::StepConfig)
    a.outdir === nothing && return nothing
    joinpath(a.outdir, "$(lpad(a.step_counter, 2, '0'))_$(step_name(cfg))")
end

function _record!(a::Analysis, cfg::StepConfig, t::Float64, summary::Dict{Symbol,Any})
    push!(a.steps, StepRecord(a.step_counter, cfg, t, summary))
end

function _save_config!(dir::String, cfg::StepConfig)
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
# Convenience: one-shot analyze (backwards compat wrapper)
# ============================================================

function analyze(data, camera::SMLMData.AbstractCamera;
                 outdir=nothing,
                 verbose=Verbosity.STANDARD,
                 # Detection
                 boxsize=11, detect_min_photons=500.0, psf_sigma=0.135, use_gpu=true,
                 # Fitting
                 psf_model=:variable, iterations=20,
                 # Filtering
                 filter=true, min_photons=500.0, max_precision=0.015, min_pvalue=1e-3,
                 # Frame connection
                 frameconnect=false, maxframegap=5,
                 # Drift
                 drift=true, degree=2,
                 # Isolated
                 isolated=false, n_sigma=2.0,
                 # Render
                 render=true, render_zoom=20)

    recipe = StepConfig[]

    push!(recipe, DetectConfig(
        boxsize=boxsize,
        min_photons=detect_min_photons,
        psf_sigma=psf_sigma,
        use_gpu=use_gpu,
        verbose=verbose
    ))

    push!(recipe, FitConfig(
        psf_model=psf_model,
        iterations=iterations,
        verbose=verbose
    ))

    if filter
        push!(recipe, FilterConfig(
            min_photons=min_photons,
            max_precision=max_precision,
            min_pvalue=min_pvalue,
            verbose=verbose
        ))
    end

    if frameconnect
        push!(recipe, FrameConnectConfig(
            maxframegap=maxframegap,
            verbose=verbose
        ))
    end

    if drift
        push!(recipe, DriftCorrectConfig(
            degree=degree,
            verbose=verbose
        ))
    end

    if isolated
        push!(recipe, IsolatedConfig(
            n_sigma=n_sigma,
            verbose=verbose
        ))
    end

    if render
        push!(recipe, RenderConfig(
            zoom=render_zoom,
            verbose=verbose
        ))
    end

    run_recipe(recipe, data, camera; outdir, verbose)
end
