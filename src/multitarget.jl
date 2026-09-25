"""
Multi-target (multi-color) analysis orchestration.

Loops over per-channel `analyze(data, config)` calls, dispatches multi-target steps
(composite renders, cross-channel alignment, etc.) on the resulting
`Vector{BasicSMLD}`, then saves the final (aligned) per-channel SMLDs.
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
        CrossAlignConfig(),   # entropy alignment (upstream AlignConfig defaults)
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
    end

    # Phase 2: Multi-target step dispatch
    composite_dir = joinpath(config.outdir, "composite")
    mkpath(composite_dir)
    step_infos = StepInfo[]

    state = smlds
    for (i, step_cfg) in enumerate(config.steps)
        colors = _resolve_colors(step_cfg, config.colors)
        (state, step_info) = analyze(state, step_cfg;
            outdir=composite_dir, step_number=i, verbose=v, colors=colors, labels=config.labels)
        push!(step_infos, step_info)
    end

    # Phase 3: save the final (aligned) per-channel SMLDs and point the channel results at them
    _finalize_channels!(channel_results, state, config.labels, config.outdir; verbose=v)

    # Write composite readme
    _write_composite_readme!(composite_dir, config, state, step_infos)

    # Save config
    _save_multitarget_config!(config)

    # Build result
    elapsed_s = (time_ns() - t_start) / 1e9
    result = MultiTargetResult(config.labels, state, channel_results, step_infos, config.outdir)
    info = MultiTargetInfo(elapsed_s, channel_infos, step_infos)

    v >= Verbosity.PROGRESS && @info "Multi-target complete: $(sum(length(s.emitters) for s in state)) total localizations ($(round(elapsed_s, digits=1))s)"

    (result, info)
end

"""
    _finalize_channels!(channel_results, state, labels, outdir; verbose) -> channel_results

Save `smld_<label>.h5` from the post-multi-target-step `state` and rebuild each
channel's `AnalysisResult` around it, so `result[label].smld` agrees with
`result.smlds`. `smld_connected` and `drift_model` stay from the channel run.
Requires that `state` is a `Vector{<:SMLMData.BasicSMLD}` with one entry per
label, in label order — the exact type `_write_composite_readme!` requires of
it right after this call. A custom multi-target step that returns anything
else (wrong length; an untyped container like `Any[...]` even if its actual
elements happen to be `BasicSMLD`s; a `view`/`SubArray` rather than a
`Vector`) leaves the label mapping undefined or breaks that type contract,
and is a contract violation either way, so this throws `ArgumentError` rather
than silently skipping the save (a caller relying on `_write_composite_readme!`
right after this would otherwise hit a `BoundsError` or `MethodError` instead
of a clear error).
"""
function _finalize_channels!(channel_results::Dict{Symbol,AnalysisResult}, state,
                             labels::Vector{Symbol}, outdir::String;
                             verbose::Int=Verbosity.STANDARD)
    if !(state isa Vector{<:SMLMData.BasicSMLD} && length(state) == length(labels))
        throw(ArgumentError("Multi-target steps must return a Vector with one BasicSMLD per channel, in label order; got $(typeof(state)) for labels $labels"))
    end
    for (i, label) in enumerate(labels)
        smld = state[i]
        cr = channel_results[label]
        smld_path = joinpath(outdir, "smld_$(label).h5")
        save_smld(smld_path, smld; drift_model=cr.drift_model)
        verbose >= Verbosity.PROGRESS && @info "  Saved $smld_path ($(length(smld.emitters)) localizations)"
        channel_results[label] = AnalysisResult(smld, cr.smld_connected, cr.drift_model)
    end
    channel_results
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
            # table_prefix="steps." so a nested config field (e.g. CompositeRenderConfig's
            # strategy, CrossAlignConfig's align) writes as `[steps.strategy]`, which TOML
            # attaches to this array-of-tables element, not a document-root `[strategy]`
            # that the next `[[steps]]` entry would collide with.
            _write_config_fields!(io, s; table_prefix="steps.")
            println(io)
        end
    end
end
