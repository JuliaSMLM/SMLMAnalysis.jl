"""
Multi-target (multi-color) analysis orchestration.

Loops over per-channel `analyze(data, config)` calls, saves per-channel SMLDs,
then dispatches multi-target steps (composite renders, cross-channel alignment, etc.)
on the resulting `Vector{BasicSMLD}`.
"""

"""
    _resolve_colors(cfg::CompositeRenderConfig, defaults::Vector{Symbol}) -> Vector{Symbol}

Use per-step colors if specified, otherwise fall back to MultiTargetConfig defaults.
"""
_resolve_colors(cfg::CompositeRenderConfig, defaults::Vector{Symbol}) = cfg.colors !== nothing ? cfg.colors : defaults
_resolve_colors(::AbstractMultiTargetStep, defaults::Vector{Symbol}) = defaults

"""
    analyze(channels::Vector{<:Tuple}, config::MultiTargetConfig) -> (MultiTargetResult, MultiTargetInfo)

Run independent analysis pipelines for each channel, then execute multi-target
steps (composite rendering, cross-channel alignment, etc.) via dispatch.

Each element of `channels` is a `(data, AnalysisConfig)` tuple where `data` is an
image stack (or Vector of stacks) or file path. The `config.labels` must match the
number of channels.

# Arguments
- `channels`: Vector of `(data, AnalysisConfig)` tuples, one per target/color
- `config`: MultiTargetConfig with labels, colors, steps, and output settings

# Returns
`(MultiTargetResult, MultiTargetInfo)` tuple following the JuliaSMLM convention.

# Example
```julia
mt = MultiTargetConfig(
    labels = [:IgG, :C1q],
    steps = [
        CompositeRenderConfig(zoom=20.0, strategy=GaussianRender()),
        CrossAlignConfig(method=:entropy),
        CompositeRenderConfig(zoom=20.0, strategy=GaussianRender()),
    ],
    outdir = "output/cell1/",
)

(result, info) = analyze([
    (image_stacks_647, config_647),
    (image_stacks_568, config_568),
], mt)

result.smlds              # Vector{BasicSMLD}
result[:IgG].smld         # Per-channel access
info.channels[:IgG]       # Per-channel AnalysisInfo
```
"""
function analyze(channels::Vector{<:Tuple}, config::MultiTargetConfig)
    t_start = time_ns()
    v = config.verbose

    # Validate inputs
    n = length(channels)
    length(config.labels) == n || error("Number of labels ($(length(config.labels))) must match number of channels ($n)")
    length(config.colors) == n || error("Number of colors ($(length(config.colors))) must match number of channels ($n)")
    length(unique(config.labels)) == n || error("Channel labels must be unique: $(config.labels)")

    mkpath(config.outdir)

    v >= Verbosity.PROGRESS && @info "Multi-target analysis: $(n) channels $(config.labels)"

    # Phase 1: Per-channel pipelines
    channel_results = Dict{Symbol, AnalysisResult}()
    channel_infos = Dict{Symbol, AnalysisInfo}()
    smlds = SMLMData.BasicSMLD[]

    for (i, (data, acfg)) in enumerate(channels)
        label = config.labels[i]
        v >= Verbosity.PROGRESS && @info "Channel $i/$n: $label"

        # Reconstruct AnalysisConfig with per-channel outdir
        channel_outdir = joinpath(config.outdir, string(label))
        mkpath(channel_outdir)
        channel_cfg = AnalysisConfig(
            camera = acfg.camera,
            steps = acfg.steps,
            roi = acfg.roi,
            outdir = channel_outdir,
            verbose = acfg.verbose,
            checkpoint = acfg.checkpoint,
        )

        (result, info) = analyze(data, channel_cfg)
        channel_results[label] = result
        channel_infos[label] = info
        push!(smlds, result.smld)

        # Save per-channel SMLD
        smld_path = joinpath(config.outdir, "smld_$(label).h5")
        save_smld(smld_path, result.smld; drift_model=result.drift_model)
        v >= Verbosity.PROGRESS && @info "  Saved $smld_path ($(length(result.smld.emitters)) localizations)"
    end

    # Phase 2: Multi-target step dispatch
    composite_dir = joinpath(config.outdir, "composite")
    mkpath(composite_dir)
    step_infos = StepInfo[]
    steps_dict = Dict{Symbol, Any}()

    state = smlds
    for (i, step_cfg) in enumerate(config.steps)
        colors = _resolve_colors(step_cfg, config.colors)
        (state, step_info) = analyze(state, step_cfg;
            outdir=composite_dir, step_number=i, verbose=v, colors=colors, labels=config.labels)
        push!(step_infos, step_info)

        # Store in steps dict (use step_name with index for duplicates)
        sname = step_name(step_cfg)
        key = haskey(steps_dict, Symbol(sname)) ? Symbol("$(sname)_$i") : Symbol(sname)
        steps_dict[key] = step_info.info
    end

    # Write composite readme
    _write_composite_readme!(composite_dir, config, state, step_infos)

    # Save config
    _save_multitarget_config!(config)

    # Build result
    elapsed_s = (time_ns() - t_start) / 1e9
    result = MultiTargetResult(config.labels, state, channel_results, step_infos, config.outdir)
    info = MultiTargetInfo(elapsed_s, channel_infos, step_infos, steps_dict)

    v >= Verbosity.PROGRESS && @info "Multi-target complete: $(sum(length(s.emitters) for s in state)) total localizations ($(round(elapsed_s, digits=1))s)"

    (result, info)
end

"""
    _write_composite_readme!(composite_dir, config, smlds, step_infos)

Write a README.md in the composite directory documenting the color scheme,
channel labels, multi-target steps, and per-channel localization counts.
"""
function _write_composite_readme!(composite_dir::String, config::MultiTargetConfig,
                                  smlds::Vector{<:SMLMData.BasicSMLD},
                                  step_infos::Vector{StepInfo})
    filepath = joinpath(composite_dir, "README.md")
    open(filepath, "w") do io
        println(io, "# Composite Output")
        println(io)
        println(io, "## Color Scheme")
        println(io)
        println(io, "| Channel | Label | Color | Localizations |")
        println(io, "|---------|-------|-------|---------------|")
        for (i, label) in enumerate(config.labels)
            color = config.colors[i]
            n = length(smlds[i].emitters)
            println(io, "| $i | $label | $color | $n |")
        end
        println(io)

        # Steps
        println(io, "## Steps")
        println(io)
        for si in step_infos
            println(io, "- **$(si.number). $(si.name)** ($(round(si.elapsed_s, digits=2))s)")
            for (k, v) in si.summary
                println(io, "  - $k: $v")
            end
        end
    end
end

"""
    _save_multitarget_config!(config::MultiTargetConfig)

Serialize MultiTargetConfig to TOML file in the output directory.
"""
function _save_multitarget_config!(config::MultiTargetConfig)
    filepath = joinpath(config.outdir, "multi_target_config.toml")
    open(filepath, "w") do io
        println(io, "# MultiTargetConfig")
        println(io, "type = \"MultiTargetConfig\"")
        println(io, "labels = [$(join(["\"$l\"" for l in config.labels], ", "))]")
        println(io, "colors = [$(join(["\"$c\"" for c in config.colors], ", "))]")
        println(io, "verbose = $(config.verbose)")
        println(io, "outdir = \"$(config.outdir)\"")
        println(io, "")
        println(io, "# Steps")
        for (i, s) in enumerate(config.steps)
            println(io, "[[steps]]")
            println(io, "type = \"$(nameof(typeof(s)))\"")
            _write_config_fields!(io, s)
            println(io)
        end
    end
end
