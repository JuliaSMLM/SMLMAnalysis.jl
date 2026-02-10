"""
Multi-target (multi-color) analysis orchestration.

Loops over per-channel `analyze(data, config)` calls, saves per-channel SMLDs,
and produces composite multi-channel renders.
"""

"""
    analyze(channels::Vector{<:Tuple}, config::MultiTargetConfig) -> (MultiTargetResult, MultiTargetInfo)

Run independent analysis pipelines for each channel and produce composite renders.

Each element of `channels` is a `(data, AnalysisConfig)` tuple where `data` is an
image stack (or Vector of stacks) or file path. The `config.labels` must match the
number of channels.

# Arguments
- `channels`: Vector of `(data, AnalysisConfig)` tuples, one per target/color
- `config`: MultiTargetConfig with labels, colors, rendering, and output settings

# Returns
`(MultiTargetResult, MultiTargetInfo)` tuple following the JuliaSMLM convention.

# Example
```julia
mt = MultiTargetConfig(
    labels = [:IgG, :C1q],
    colors = [:red, :green],
    render_zoom = 20,
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

    # Run per-channel pipelines
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

    # Composite renders
    composite_dir = joinpath(config.outdir, "composite")
    mkpath(composite_dir)
    composite_render_infos = SMLMRender.RenderInfo[]

    for strategy in config.render_strategies
        strategy_name = lowercase(string(nameof(typeof(strategy))))
        zoom_str = "$(round(Int, config.render_zoom))x"
        filename = joinpath(composite_dir, "$(strategy_name)_$(zoom_str).png")

        v >= Verbosity.PROGRESS && @info "Composite render: $strategy_name $(zoom_str)"

        # Histogram overlays: saturate mode (count=1 = full brightness, clamp handles >1)
        # Other strategies: clip+normalize for smooth intensity scaling
        is_histogram = strategy isa HistogramRender
        cp = is_histogram ? nothing : config.clip_percentile
        ne = is_histogram ? false : true

        (_, render_info) = SMLMRender.render(smlds;
            colors = config.colors,
            strategy = strategy,
            zoom = config.render_zoom,
            clip_percentile = cp,
            normalize_each = ne,
            filename = filename,
        )
        push!(composite_render_infos, render_info)
    end

    # Save config
    _save_multitarget_config!(config)

    # Build result
    elapsed_s = (time_ns() - t_start) / 1e9
    result = MultiTargetResult(config.labels, smlds, channel_results, nothing, config.outdir)
    info = MultiTargetInfo(elapsed_s, channel_infos, composite_render_infos)

    v >= Verbosity.PROGRESS && @info "Multi-target complete: $(sum(length(s.emitters) for s in smlds)) total localizations ($(round(elapsed_s, digits=1))s)"

    (result, info)
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
        println(io, "render_zoom = $(config.render_zoom)")
        println(io, "verbose = $(config.verbose)")
        println(io, "outdir = \"$(config.outdir)\"")
        println(io, "")
        println(io, "# Rendering strategies")
        for (i, s) in enumerate(config.render_strategies)
            println(io, "[[render_strategies]]")
            println(io, "type = \"$(nameof(typeof(s)))\"")
        end
    end
end
