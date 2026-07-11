"""
Composite render step — multi-channel rendering for the multi-target pipeline.

Dispatches on `CompositeRenderConfig <: AbstractMultiTargetStep` operating on
`Vector{BasicSMLD}`. Pass-through: SMLDs are not modified.
"""

"""
    CompositeRenderConfig <: AbstractMultiTargetStep

Configuration for a composite multi-channel render in the multi-target pipeline.

# Fields
- `strategy`: Rendering strategy (default: GaussianRender())
- `zoom`: Zoom factor (default: 20.0)
- `colors`: Per-channel colors (nothing = inherit from MultiTargetConfig)
- `clip_percentile`: Intensity clipping. `:auto` (default) picks per strategy —
  saturate (no clip) for histogram, 0.99 for others. A `Float64` clips at that
  percentile; `nothing` forces saturate mode.
- `normalize_each`: Per-channel normalization (nothing = auto: false for histogram, true for others)
- `scalebar`: Enable scale bar (default: true)
- `scalebar_length`: Scale bar length in μm (nothing = auto)
- `scalebar_position`: Scale bar corner (default: :br)
- `scalebar_color`: Scale bar color (default: :white)
"""
@kwdef struct CompositeRenderConfig <: AbstractMultiTargetStep
    strategy::SMLMRender.RenderingStrategy = GaussianRender()
    zoom::Float64 = 20.0
    colors::Union{Vector{Symbol}, Nothing} = nothing
    clip_percentile::Union{Float64, Nothing, Symbol} = :auto
    normalize_each::Union{Bool, Nothing} = nothing
    scalebar::Bool = true
    scalebar_length::Union{Float64, Nothing} = nothing
    scalebar_position::Symbol = :br
    scalebar_color::Symbol = :white
end

step_name(::CompositeRenderConfig) = "compositerender"

"""
    composite_render_step(smlds, cfg; outdir, step_number, verbose, colors) -> (smlds, CompositeRenderInfo)

Render a multi-channel composite image. SMLDs pass through unmodified.
"""
function composite_render_step(smlds::Vector{<:SMLMData.BasicSMLD}, cfg::CompositeRenderConfig;
                               outdir::Union{String,Nothing}=nothing,
                               step_number::Int=0,
                               verbose::Int=Verbosity.STANDARD,
                               colors::Vector{Symbol}=Symbol[])
    v = verbose
    dir = step_outdir(outdir, step_number, cfg)

    strategy_name = lowercase(string(nameof(typeof(cfg.strategy))))
    zoom_str = "$(round(Int, cfg.zoom))x"
    v >= Verbosity.PROGRESS && @info "[$step_number] compositerender: $strategy_name $zoom_str"

    # Resolve normalize_each: histogram=false (saturate), others=true (clip+normalize)
    is_histogram = cfg.strategy isa HistogramRender
    ne = cfg.normalize_each !== nothing ? cfg.normalize_each : !is_histogram
    # clip_percentile :auto → per-strategy default (histogram saturates, others clip at
    # 0.99); an explicit Float64/nothing is honored as-is (a user 0.99 is no longer
    # silently reinterpreted as saturate).
    cp = cfg.clip_percentile === :auto ? (is_histogram ? nothing : 0.99) : cfg.clip_percentile

    # Build filename
    filename = if dir !== nothing
        mkpath(dir)
        joinpath(dir, "$(strategy_name)_$(zoom_str).png")
    else
        nothing
    end

    local render_info
    t = @elapsed begin
        (_, render_info) = SMLMRender.render(smlds;
            colors=colors,
            strategy=cfg.strategy,
            zoom=cfg.zoom,
            clip_percentile=cp,
            normalize_each=ne,
            filename=filename,
            scalebar=cfg.scalebar,
            scalebar_length=cfg.scalebar_length,
            scalebar_position=cfg.scalebar_position,
            scalebar_color=cfg.scalebar_color,
        )
    end

    if dir !== nothing
        _save_config!(dir, cfg)
        _save_info!(dir, render_info)
        if v >= Verbosity.STANDARD
            _write_composite_render_stats(dir, cfg, render_info, smlds, t)
        end
    end

    v >= Verbosity.PROGRESS && @info "  -> composite $(render_info.output_size) ($(round(t, digits=2))s)"

    info = CompositeRenderInfo(render_info, Symbol(strategy_name), cfg.zoom, length(smlds), t)
    (smlds, info)
end

_step_summary(info::CompositeRenderInfo) = Dict{Symbol,Any}(
    :strategy => info.strategy,
    :zoom => info.zoom,
    :n_channels => info.n_channels,
    :output_size => info.render_info.output_size,
)

"""
    analyze(smlds::Vector{BasicSMLD}, cfg::CompositeRenderConfig; kwargs...) -> (smlds, StepInfo)

Multi-target dispatch: composite render. SMLDs pass through.
"""
function analyze(smlds::Vector{<:SMLMData.BasicSMLD}, cfg::CompositeRenderConfig;
                 outdir=nothing, step_number::Int=0, verbose::Int=Verbosity.STANDARD,
                 colors::Vector{Symbol}=Symbol[], kwargs...)
    t = @elapsed (smlds, cr_info) = composite_render_step(smlds, cfg;
        outdir=outdir, step_number=step_number, verbose=verbose, colors=colors)
    (smlds, StepInfo(step_number, cfg, t, _step_summary(cr_info); info=cr_info))
end

function _write_composite_render_stats(dir, cfg::CompositeRenderConfig, render_info, smlds, t)
    filepath = joinpath(dir, "stats.md")
    open(filepath, "w") do io
        println(io, "# Composite Render Statistics\n")
        println(io, "## Summary")
        println(io, "- **Strategy**: $(nameof(typeof(cfg.strategy)))")
        println(io, "- **Zoom**: $(cfg.zoom)x")
        println(io, "- **Channels**: $(length(smlds))")
        for (i, smld) in enumerate(smlds)
            println(io, "  - Channel $i: $(length(smld.emitters)) localizations")
        end
        println(io, "- **Output size**: $(render_info.output_size)")
        println(io, "- **Pixel size**: $(round(render_info.pixel_size_nm, digits=1)) nm")
        println(io, "- **Time**: $(round(t, digits=2))s")
    end
end
