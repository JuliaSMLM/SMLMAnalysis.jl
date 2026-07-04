"""
Render step - dispatches directly on SMLMRender.RenderConfig.

Each render is one step call with its own numbered output folder:
```julia
render_step(smld, RenderConfig(zoom=20, colormap=:inferno))
render_step(smld, RenderConfig(zoom=10, colormap=:turbo, color_by=:absolute_frame))
```
"""

# Override step_name so dirs are "07_render" not "07_renderconfig"
step_name(::SMLMRender.RenderConfig) = "render"

"""
    render_step(smld, cfg; outdir=nothing, step_number=0, verbose=Verbosity.STANDARD)

Render localizations to a super-resolution image. Returns `(render_image, RenderInfo)`.

# Arguments
- `smld::BasicSMLD`: Input localizations
- `cfg::SMLMRender.RenderConfig`: Render configuration

# Keyword Arguments
- `outdir`: Output directory (nothing to skip file output)
- `step_number`: Step number for output directory naming
- `verbose`: Verbosity level

# Returns
`(render_image, RenderInfo)`
"""
function render_step(smld::BasicSMLD, cfg::SMLMRender.RenderConfig;
                     outdir::Union{String,Nothing}=nothing,
                     step_number::Int=0,
                     verbose::Int=Verbosity.STANDARD)
    v = verbose
    dir = step_outdir(outdir, step_number, cfg)

    v >= Verbosity.PROGRESS && @info "[$step_number] render" zoom=cfg.zoom colormap=cfg.colormap

    # Set filename into config if output dir exists and user didn't set one
    render_cfg = if dir !== nothing && cfg.filename === nothing
        mkpath(dir)
        # Build descriptive filename
        strategy_name = lowercase(string(nameof(typeof(cfg.strategy))))
        color_suffix = cfg.color_by === nothing ? "" : "_$(cfg.color_by)"
        cmap = cfg.colormap === nothing ? "" : "_$(cfg.colormap)"
        zoom_str = cfg.zoom === nothing ? "" : "_$(round(Int, cfg.zoom))x"
        filename = joinpath(dir, "$(strategy_name)$(cmap)$(color_suffix)$(zoom_str).png")
        # Reconstruct with filename set (immutable struct) - use kwargs to be robust to field additions
        SMLMRender.RenderConfig(;
            strategy=cfg.strategy, pixel_size=cfg.pixel_size, zoom=cfg.zoom,
            roi=cfg.roi, target=cfg.target, colormap=cfg.colormap,
            color_by=cfg.color_by, color=cfg.color, categorical=cfg.categorical,
            clip_percentile=cfg.clip_percentile, field_range=cfg.field_range,
            field_clip_percentiles=cfg.field_clip_percentiles,
            backend=cfg.backend, filename=filename,
            scalebar=cfg.scalebar, scalebar_length=cfg.scalebar_length,
            scalebar_position=cfg.scalebar_position, scalebar_color=cfg.scalebar_color,
        )
    else
        cfg
    end

    # Tuple-pattern: returns (image, RenderInfo)
    local render_info, render_image
    t = @elapsed begin
        (render_image, render_info) = SMLMRender.render(smld, render_cfg)
    end

    n_locs = length(smld.emitters)

    if dir !== nothing
        mkpath(dir)
        _save_config!(dir, cfg)
        _save_info!(dir, render_info)
        if v >= Verbosity.STANDARD
            _write_render_stats(dir, cfg, render_info, n_locs, t)
        end
    end

    v >= Verbosity.PROGRESS && @info "  → render $(render_info.output_size) ($(round(t, digits=2))s)"
    (render_image, render_info)
end

_step_summary(info::SMLMRender.RenderInfo) = Dict{Symbol,Any}(
    :n_locs => info.n_emitters_rendered,
    :strategy => info.strategy,
    :output_size => info.output_size
)

"""
    analyze(smld, cfg::RenderConfig; kwargs...) -> (smld, StepInfo)

Render localizations to a super-resolution image. The image is saved to disk
(via render_step); the smld passes through so subsequent pipeline steps can
operate on it. Use `render_step` or `SMLMRender.render` directly to get the image.
"""
function analyze(smld::BasicSMLD, cfg::SMLMRender.RenderConfig;
                 outdir=nothing, step_number::Int=0, verbose::Int=Verbosity.STANDARD, kwargs...)
    t = @elapsed (render_image, render_info) = render_step(smld, cfg;
        outdir=outdir, step_number=step_number, verbose=verbose)
    (smld, StepInfo(step_number, cfg, t, _step_summary(render_info); info=render_info))
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
