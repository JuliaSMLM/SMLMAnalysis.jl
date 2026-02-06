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
    RenderSpec(strategy=:histogram, zoom=10, colormap=:turbo, color_by=:absolute_frame, clip_percentile=0.999),
    RenderSpec(strategy=:circles, zoom=50, colormap=:turbo, color_by=:absolute_frame, clip_percentile=0.999),
]

@kwdef struct RenderConfig <: SMLMData.AbstractSMLMConfig
    # Render specifications - vector of RenderSpec or single spec via legacy fields
    renders::Vector{RenderSpec} = DEFAULT_RENDERS
    # Legacy single-render fields (used if renders is empty or for backward compat)
    strategy::Symbol = :gaussian
    zoom::Int = 20
    colormap::Symbol = :inferno
    color_by::Union{Symbol, Nothing} = nothing
    clip_percentile::Union{Float64, Symbol} = :auto
end

function run_step!(a::Analysis, cfg::RenderConfig)
    a.smld === nothing && error("Must run Fit first")
    a.step_counter += 1
    v = a.verbose
    dir = _stepdir(a, cfg)

    # Use renders list if non-empty, otherwise fall back to single spec from legacy fields
    specs = if !isempty(cfg.renders)
        cfg.renders
    else
        [RenderSpec(strategy=cfg.strategy, zoom=cfg.zoom, colormap=cfg.colormap, color_by=cfg.color_by)]
    end

    v >= Verbosity.PROGRESS && @info "[$(a.step_counter)] $(step_name(cfg))" n_renders=length(specs)

    # Collect render_info from each render (tuple-pattern)
    all_render_info = []
    t = @elapsed for spec in specs
        render_info = _render_single(a.smld, spec, dir, v)
        push!(all_render_info, render_info)
    end

    n_locs = length(a.smld.emitters)
    summary = Dict{Symbol,Any}(
        :n_locs => n_locs,
        :n_renders => length(specs),
        :renders => [(s.strategy, s.zoom, s.colormap, s.color_by) for s in specs]
    )

    # Aggregate render info (tuple-pattern)
    step_info = (
        render_info = all_render_info,
        elapsed_s = t
    )
    _record!(a, cfg, t, summary; info=step_info)

    if dir !== nothing
        _save_step_outputs!(dir, a, cfg, v, t, n_locs, specs, all_render_info)
    end

    v >= Verbosity.PROGRESS && @info "  → $(length(specs)) renders ($(round(t, digits=2))s)"
    a
end

"""
Render single specification, returns RenderInfo (tuple-pattern).
"""
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

    # Tuple-pattern: returns (image, RenderInfo)
    (image, render_info) = SMLMRender.render(smld; render_kwargs...)

    v >= Verbosity.DETAILED && @info "    $(spec.strategy) $(spec.colormap) $(spec.zoom)x" color_by=spec.color_by

    return render_info
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

function _save_step_outputs!(dir::String, a::Analysis, cfg::RenderConfig, v::Int, t::Float64, n_locs::Int, specs::Vector{RenderSpec}, all_render_info)
    mkpath(dir)
    _save_config!(dir, cfg)

    # Write upstream info structs to info.toml
    # Write header first, then append sections
    open(joinpath(dir, "info.toml"), "w") do io
        println(io, "# Upstream package info")
    end
    n_renders = length(all_render_info)
    if n_renders == 1
        _save_info!(dir, all_render_info[1]; section="render_info")
    else
        for i in 1:n_renders
            _save_info!(dir, all_render_info[i]; section="render_info_$i")
        end
    end

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

