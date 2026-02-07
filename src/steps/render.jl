"""
Render step - dispatches directly on SMLMRender.RenderConfig.

Each render is one step call with its own numbered output folder:
```julia
run_step!(a, RenderConfig(zoom=20, colormap=:inferno))
run_step!(a, RenderConfig(zoom=10, colormap=:turbo, color_by=:absolute_frame))
```
"""

# Override step_name so dirs are "07_render" not "07_renderconfig"
step_name(::SMLMRender.RenderConfig) = "render"

function run_step!(a::Analysis, cfg::SMLMRender.RenderConfig)
    a.smld === nothing && error("Must run Fit first")
    a.step_counter += 1
    v = a.verbose
    dir = _stepdir(a, cfg)

    v >= Verbosity.PROGRESS && @info "[$(a.step_counter)] render" zoom=cfg.zoom colormap=cfg.colormap

    # Set filename into config if output dir exists and user didn't set one
    render_cfg = if dir !== nothing && cfg.filename === nothing
        mkpath(dir)
        # Build descriptive filename
        strategy_name = lowercase(string(nameof(typeof(cfg.strategy))))
        color_suffix = cfg.color_by === nothing ? "" : "_$(cfg.color_by)"
        cmap = cfg.colormap === nothing ? "" : "_$(cfg.colormap)"
        zoom_str = cfg.zoom === nothing ? "" : "_$(round(Int, cfg.zoom))x"
        filename = joinpath(dir, "$(strategy_name)$(cmap)$(color_suffix)$(zoom_str).png")
        # Reconstruct with filename set (immutable struct)
        SMLMRender.RenderConfig(
            cfg.strategy, cfg.pixel_size, cfg.zoom, cfg.roi, cfg.target,
            cfg.colormap, cfg.color_by, cfg.color, cfg.categorical,
            cfg.clip_percentile, cfg.field_range, cfg.field_clip_percentiles,
            cfg.backend, filename
        )
    else
        cfg
    end

    # Tuple-pattern: returns (image, RenderInfo)
    local render_info
    t = @elapsed begin
        (_, render_info) = SMLMRender.render(a.smld, render_cfg)
    end

    n_locs = length(a.smld.emitters)
    summary = Dict{Symbol,Any}(
        :n_locs => n_locs,
        :strategy => render_info.strategy,
        :output_size => render_info.output_size
    )

    _record!(a, cfg, t, summary; info=render_info)

    if dir !== nothing
        mkpath(dir)
        _save_config!(dir, cfg)
        _save_info!(dir, render_info)
        if v >= Verbosity.STANDARD
            _write_render_stats(dir, cfg, render_info, n_locs, t)
        end
    end

    v >= Verbosity.PROGRESS && @info "  → render $(render_info.output_size) ($(round(t, digits=2))s)"
    a
end

function _write_render_stats(dir, cfg::SMLMRender.RenderConfig, render_info, n_locs, t)
    filepath = joinpath(dir, "stats.md")
    open(filepath, "w") do io
        println(io, "# Render Statistics\n")
        println(io, "## Summary")
        println(io, "- **Localizations**: $n_locs")
        println(io, "- **Emitters rendered**: $(render_info.n_emitters_rendered)")
        println(io, "- **Output size**: $(render_info.output_size)")
        println(io, "- **Pixel size**: $(round(render_info.pixel_size_nm, digits=1)) nm")
        println(io, "- **Strategy**: $(render_info.strategy)")
        println(io, "- **Color mode**: $(render_info.color_mode)")
        println(io, "- **Time**: $(round(t, digits=2))s")
    end
end
