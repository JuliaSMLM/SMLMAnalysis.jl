"""
Render step - wraps SMLMRender.render
"""

"""
Single render specification for multi-render output.
"""
@kwdef struct RenderSpec
    strategy::Symbol = :gaussian    # :gaussian, :histogram, :circles
    zoom::Int = 20
    colormap::Symbol = :inferno
    color_by::Union{Symbol, Nothing} = nothing  # nothing=density, or :frame, :photons, :z, etc.
    clip_percentile::Float64 = 0.999  # Clip at this percentile for intensity scaling
end

# Default renders: gaussian density, histogram time, circles time
const DEFAULT_RENDERS = [
    RenderSpec(strategy=:gaussian, zoom=20, colormap=:inferno, color_by=nothing, clip_percentile=0.999),
    RenderSpec(strategy=:histogram, zoom=10, colormap=:turbo, color_by=:frame, clip_percentile=0.999),
    RenderSpec(strategy=:circles, zoom=50, colormap=:turbo, color_by=:frame, clip_percentile=0.999),
]

@kwdef struct RenderConfig <: StepConfig
    name::String = "render"
    # Render specifications - vector of RenderSpec or single spec via legacy fields
    renders::Vector{RenderSpec} = DEFAULT_RENDERS
    # Legacy single-render fields (used if renders is empty or for backward compat)
    strategy::Symbol = :gaussian
    zoom::Int = 20
    colormap::Symbol = :inferno
    color_by::Union{Symbol, Nothing} = nothing
    clip_percentile::Union{Float64, Symbol} = :auto
    # Extra
    verbose::Int = Verbosity.STANDARD
end

function run_step!(a::Analysis, cfg::RenderConfig)
    a.smld === nothing && error("Must run Fit first")
    a.step_counter += 1
    v = _get_verbose(a, cfg)
    dir = _stepdir(a, cfg)

    # Use renders list if non-empty, otherwise fall back to single spec from legacy fields
    specs = if !isempty(cfg.renders)
        cfg.renders
    else
        [RenderSpec(strategy=cfg.strategy, zoom=cfg.zoom, colormap=cfg.colormap, color_by=cfg.color_by)]
    end

    v >= Verbosity.PROGRESS && @info "[$(a.step_counter)] $(cfg.name)" n_renders=length(specs)

    t = @elapsed for spec in specs
        _render_single(a.smld, spec, dir, v)
    end

    n_locs = length(a.smld.emitters)
    summary = Dict{Symbol,Any}(
        :n_locs => n_locs,
        :n_renders => length(specs),
        :renders => [(s.strategy, s.zoom, s.colormap, s.color_by) for s in specs]
    )
    _record!(a, cfg, t, summary)

    if dir !== nothing
        _save_step_outputs!(dir, a, cfg, v, t, n_locs, specs)
    end

    v >= Verbosity.PROGRESS && @info "  → $(length(specs)) renders ($(round(t, digits=2))s)"
    a
end

function _render_single(smld, spec::RenderSpec, dir::Union{String,Nothing}, v::Int)
    strategy_obj = if spec.strategy == :gaussian
        GaussianRender()
    elseif spec.strategy == :histogram
        HistogramRender()
    elseif spec.strategy == :circles
        CircleRender()
    else
        error("Unknown strategy: $(spec.strategy). Use :gaussian, :histogram, or :circles")
    end

    # Build filename
    color_suffix = spec.color_by === nothing ? "" : "_$(spec.color_by)"
    filename = dir === nothing ? nothing : joinpath(dir, "$(spec.strategy)_$(spec.colormap)$(color_suffix)_$(spec.zoom)x.png")

    # Build render kwargs
    render_kwargs = Dict{Symbol, Any}(
        :strategy => strategy_obj,
        :zoom => spec.zoom,
        :colormap => spec.colormap,
        :clip_percentile => spec.clip_percentile,
        :filename => filename
    )
    if spec.color_by !== nothing
        render_kwargs[:color_by] = spec.color_by
    end

    SMLMRender.render(smld; render_kwargs...)

    v >= Verbosity.DETAILED && @info "    $(spec.strategy) $(spec.colormap) $(spec.zoom)x" color_by=spec.color_by
end

function _adaptive_clip_percentile(n_locs::Int)
    if n_locs < 50_000
        0.99
    elseif n_locs < 200_000
        0.995
    elseif n_locs < 500_000
        0.999
    else
        0.9999
    end
end

function _save_step_outputs!(dir::String, a::Analysis, cfg::RenderConfig, v::Int, t::Float64, n_locs::Int, specs::Vector{RenderSpec})
    mkpath(dir)
    _save_config!(dir, cfg)

    if v >= Verbosity.STANDARD
        _write_render_stats(dir, cfg, n_locs, t, specs)
    end
end

function _write_render_stats(dir, cfg, n_locs, t, specs)
    filepath = joinpath(dir, "stats.md")
    open(filepath, "w") do io
        println(io, "# Render Statistics\n")
        println(io, "## Summary")
        println(io, "- **Localizations**: $n_locs")
        println(io, "- **Renders**: $(length(specs))")
        println(io, "- **Time**: $(round(t, digits=2))s")
        println(io, "")
        println(io, "## Renders")
        println(io, "| Strategy | Zoom | Colormap | Color By | Clip |")
        println(io, "|----------|------|----------|----------|------|")
        for spec in specs
            color_by_str = spec.color_by === nothing ? "density" : string(spec.color_by)
            println(io, "| $(spec.strategy) | $(spec.zoom)x | $(spec.colormap) | $color_by_str | $(spec.clip_percentile) |")
        end
    end
end

