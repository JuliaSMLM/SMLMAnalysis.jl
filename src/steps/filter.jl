"""
Filtering step - filters localizations by various criteria
"""

"""
    FilterConfig <: AbstractSMLMConfig

Quality-based filtering of localizations. All criteria use `(min, max)` tuples.

# Keywords
- `photons`: Photon count range, e.g. `(500.0, Inf)`
- `precision`: Lateral localization precision range in microns, `max(σ_x, σ_y)`, e.g. `(0.0, 0.007)`
- `pvalue`: Goodness-of-fit p-value range, e.g. `(1e-3, 1.0)`. Warns and has no
  effect on emitters without p-values (frame-connected, BaGoL, or other output).
- `psf_sigma`: PSF width filter. `:auto` uses mode ± 10%, or explicit `(min, max)` in microns
- `z`: Axial position range in microns (3D fits only), e.g. `(-0.4, 0.4)`. Useful to
  drop localizations railed to the z-grid edge. No-op on 2D SMLDs.
- `sigma_z`: Axial precision range in microns (3D fits only), e.g. `(0.0, 0.030)`. The
  axial analog of `precision`. No-op on 2D SMLDs.

All filters default to `nothing` (disabled).
"""
@kwdef struct FilterConfig <: SMLMData.AbstractSMLMConfig
    # All filters use (min, max) tuples. Use -Inf/Inf for unbounded.
    photons::Union{Tuple{Float64, Float64}, Nothing} = nothing      # (min, max)
    precision::Union{Tuple{Float64, Float64}, Nothing} = nothing    # (min, max) in microns
    pvalue::Union{Tuple{Float64, Float64}, Nothing} = nothing       # (min, max)
    # PSF sigma: :auto (mode ± 10%), or (min, max) tuple in microns
    psf_sigma::Union{Symbol, Tuple{Float64, Float64}, Nothing} = nothing
    # 3D-only (Emitter3DFit/…PV); no-op when emitters lack z / σ_z.
    z::Union{Tuple{Float64, Float64}, Nothing} = nothing            # (min, max) axial position, microns
    sigma_z::Union{Tuple{Float64, Float64}, Nothing} = nothing      # (min, max) axial precision, microns

    function FilterConfig(photons, precision, pvalue, psf_sigma, z, sigma_z)
        # A Union with two non-Nothing alternatives does not convert tuple elements.
        bounds = psf_sigma isa Tuple{Real, Real} ? convert(Tuple{Float64, Float64}, psf_sigma) : psf_sigma
        return new(photons, precision, pvalue, bounds, z, sigma_z)
    end
end

"""
    filter_step(smld, cfg; outdir=nothing, step_number=0, verbose=Verbosity.STANDARD)

Filter localizations by quality criteria. Returns `(filtered_smld, FilterInfo)`.

Diagnostic plots use `smld` (the input) for pre-filter distributions.

# Arguments
- `smld::BasicSMLD`: Input localizations
- `cfg::FilterConfig`: Filter criteria

# Keyword Arguments
- `outdir`: Output directory (nothing to skip file output)
- `step_number`: Step number for output directory naming
- `verbose`: Verbosity level

# Returns
`(filtered_smld, FilterInfo)`
"""
function filter_step(
        smld::BasicSMLD, cfg::FilterConfig;
        outdir::Union{String, Nothing} = nothing,
        step_number::Int = 0,
        verbose::Int = Verbosity.STANDARD
    )
    v = verbose
    dir = step_outdir(outdir, step_number, cfg)

    v >= Verbosity.PROGRESS && @info "[$step_number] $(step_name(cfg))" photons = cfg.photons precision = cfg.precision

    n_before = length(smld.emitters)
    t = @elapsed filtered = _filter_smld(smld, cfg)
    n_after = length(filtered.emitters)

    if dir !== nothing
        _save_filter_outputs!(dir, outdir, cfg, v, t, n_before, n_after, smld, filtered)
    end

    v >= Verbosity.PROGRESS && @info "  → $n_after / $n_before ($(round(t, digits = 2))s)"
    return (filtered, FilterInfo(n_before, n_after, t))
end

_step_summary(info::FilterInfo) = Dict{Symbol, Any}(
    :n_before => info.n_before,
    :n_after => info.n_after,
    :acceptance => round(info.n_after / max(1, info.n_before), digits = 3)
)

"""
    analyze(smld, cfg::FilterConfig; outdir, step_number, verbose, checkpoint) -> (filtered_smld, StepInfo)

Filter localizations by quality criteria.
"""
function analyze(
        smld::BasicSMLD, cfg::FilterConfig;
        outdir = nothing, step_number::Int = 0, verbose::Int = Verbosity.STANDARD,
        checkpoint::Int = Checkpoint.EXPENSIVE
    )
    t = @elapsed (filtered, filter_info) = filter_step(
        smld, cfg;
        outdir = outdir, step_number = step_number, verbose = verbose
    )

    if checkpoint >= Checkpoint.ALL
        dir = step_outdir(outdir, step_number, cfg)
        _save_step_smld(dir, filtered; filename = "smld_filtered.h5")
    end

    return (filtered, StepInfo(step_number, cfg, t, _step_summary(filter_info); info = filter_info))
end

# The (lo, hi) range test FilterConfig applies. A NaN fails it, so the fit overlay
# (_fit_box_color) uses the same predicate to colour exactly what the filter rejects.
_in_range(v, (lo, hi)) = lo <= v <= hi

# Percentage of `values` the filter keeps for `range`, by the same test the filter uses.
_pass_percent(values, range) = round(100 * count(v -> _in_range(v, range), values) / length(values), digits = 1)

# Smallest and largest non-NaN value (either may be infinite), or (NaN, NaN) if there is none.
function _value_extent(values)
    lo, hi = Inf, -Inf
    for v in values
        isnan(v) && continue
        lo = min(lo, v)
        hi = max(hi, v)
    end
    return lo <= hi ? (Float64(lo), Float64(hi)) : (NaN, NaN)
end

# Bounds whose magnitude cannot be drawn on a Float32 axis do not set the view.
const _PANEL_MAX_BOUND = Float64(floatmax(Float32)) / 4

"""
    _panel_geometry(window, extent, bounds; nbins = 50) -> (; view, edges, lines)

The x view, bin edges and threshold lines of one fit_quality.png histogram panel.

- `window = (lo, hi)`: what the panel shows without a filter (data quantiles).
- `extent = (min, max)`: the data's full range, to tell whether a bound rejects any value.
- `bounds`: the filter ranges `(lo, hi)` drawn on this panel, in the panel's units.

A finite bound is drawn when it rejects a value (`lo > min` or `hi < max`) or lies within one
window width of the window; a bound beyond that rejects nothing and is left off, as is one at or
below `floor` that rejects nothing (a lower bound of 0 on photons or a precision). The view
covers the window and every drawn bound; a side that a bound extends (or touches) gets 5% of the
view's width as margin, so each drawn bound lies strictly inside the view. `edges` are about
`nbins` equal bins over the view with every drawn bound inserted as an edge, so no bar straddles
a threshold. `lines` are the drawn bounds.
"""
function _panel_geometry(window, extent, bounds; nbins::Int = 50, floor = -Inf)
    dlo, dhi = _plot_limits(Float64(window[1]), Float64(window[2]))
    dmin, dmax = extent
    reach = dhi - dlo
    lines = Float64[]
    for (lo, hi) in bounds, (b, rejects) in ((lo, lo > dmin), (hi, hi < dmax))
        (isfinite(b) && abs(b) <= _PANEL_MAX_BOUND) || continue
        (rejects || (b > floor && dlo - reach <= b <= dhi + reach)) && push!(lines, Float64(b))
    end
    lo_ext = minimum(lines; init = dlo)
    hi_ext = maximum(lines; init = dhi)
    pad = 0.05 * (hi_ext - lo_ext)
    vlo = any(<=(dlo), lines) ? lo_ext - pad : dlo
    vhi = any(>=(dhi), lines) ? hi_ext + pad : dhi
    view = _plot_limits(vlo, vhi)
    lines = sort!(unique!(filter(t -> view[1] < t < view[2], lines)))
    w = (view[2] - view[1]) / nbins
    uniform = range(view[1], view[2], length = nbins + 1)
    interior = filter(e -> all(t -> abs(e - t) >= w / 4, lines), uniform[2:(end - 1)])
    edges = sort!(unique!(vcat(view[1], interior, lines, view[2])))
    return (; view, edges, lines)
end

"""
    _panel_spans(view, bounds, sides) -> Vector{Tuple{Float64, Float64}}

The parts of `view` that the filter range `bounds = (lo, hi)` rejects, on the `sides` it shades:
`(view_lo, lo)` below a lower bound and `(hi, view_hi)` above an upper bound, clipped to the view.
A range that rejects everything (`lo = Inf`, `hi = -Inf`, a NaN end) shades the whole view.
"""
function _panel_spans(view, bounds, sides)
    vlo, vhi = view
    lo, hi = bounds
    spans = Tuple{Float64, Float64}[]
    if :lo in sides && (isnan(lo) || lo > vlo)
        push!(spans, (vlo, isnan(lo) ? vhi : min(Float64(lo), vhi)))
    end
    if :hi in sides && (isnan(hi) || hi < vhi)
        push!(spans, (isnan(hi) ? vlo : max(Float64(hi), vlo), vhi))
    end
    return spans
end

"""
    _panel_counts(values, edges, bounds, spans) -> (; counts, beyond)

Bar counts over `edges` (bins `[edges[i], edges[i+1])`, the top edge in the last bin), compared
with the stored edges rather than divided by the bin width. A value on an upper bound passes the
filter, so it goes in the bar below that bound, outside the shaded span. Values beyond the view on
a side whose shaded span reaches the view's end are rejected, so they are drawn in that end's bar
(`beyond` counts them); on any other side they are not drawn, as the window clips them.
"""
function _panel_counts(values, edges, bounds, spans)
    nb = length(edges) - 1
    counts = zeros(Int, nb)
    vlo, vhi = first(edges), last(edges)
    pile_lo = any(s -> s[1] == vlo, spans)
    pile_hi = any(s -> s[2] == vhi, spans)
    hi = bounds === nothing ? NaN : Float64(bounds[2])
    rejected(v) = bounds !== nothing && !_in_range(v, bounds)
    beyond = 0
    for v in values
        isnan(v) && continue
        if v < vlo
            pile_lo && rejected(v) && (counts[1] += 1; beyond += 1)
        elseif v > vhi
            pile_hi && rejected(v) && (counts[nb] += 1; beyond += 1)
        else
            i = min(searchsortedlast(edges, v), nb)
            v == hi && i > 1 && edges[i] == v && (i -= 1)
            counts[i] += 1
        end
    end
    return (; counts, beyond)
end

"""
    _histogram_panel(window, series; density = false)

Everything one fit_quality.png histogram panel draws, as data. `series` holds one
`(; name, values, bounds, sides)` per overlaid histogram: values in the panel's units, the
filter range applied to them (or `nothing`) and the sides of it that reject on their own
(`(:lo, :hi)`; the precision panel shades only `(:hi,)`, since the filter tests max(σx, σy)).
Returns `(; view, edges, density, series)`, each series with its `spans`, `lines` (bounds on
shaded sides), `dashed` (bounds on unshaded sides), `counts` and `beyond`.
"""
function _histogram_panel(window, series; density::Bool = false, floor = -Inf)
    extent = _value_extent(reduce(vcat, [s.values for s in series]; init = Float64[]))
    g = _panel_geometry(window, extent, [s.bounds for s in series if s.bounds !== nothing]; floor)
    out = map(series) do s
        if s.bounds === nothing
            c = _panel_counts(s.values, g.edges, nothing, Tuple{Float64, Float64}[])
            return (;
                s.name, s.values, s.bounds, s.sides, spans = Tuple{Float64, Float64}[],
                lines = Float64[], dashed = Float64[], c.counts, c.beyond,
            )
        end
        spans = _panel_spans(g.view, s.bounds, s.sides)
        on_line(b) = any(==(b), g.lines)
        lo, hi = s.bounds
        lines = Float64[b for (b, side) in ((lo, :lo), (hi, :hi)) if side in s.sides && on_line(b)]
        dashed = Float64[b for (b, side) in ((lo, :lo), (hi, :hi)) if !(side in s.sides) && on_line(b)]
        c = _panel_counts(s.values, g.edges, s.bounds, spans)
        return (; s.name, s.values, s.bounds, s.sides, spans, lines, dashed, c.counts, c.beyond)
    end
    return (; g.view, g.edges, density, series = out)
end

_um_to_nm(b) = b === nothing ? nothing : (b[1] * 1000, b[2] * 1000)

# p-value bounds in log10 space, where p = 0 is -Inf: a lower bound at or below 0 and an upper
# bound at or above 1 reject nothing; an upper bound of 0 keeps only p = 0 (-Inf); a negative
# upper bound keeps nothing (NaN fails every comparison). No logarithm of a non-positive number.
function _log10_pvalue_bounds(bounds)
    bounds === nothing && return nothing
    p_lo, p_hi = bounds
    lo = p_lo > 0 ? log10(p_lo) : -Inf
    hi = p_hi >= 1 ? Inf : (p_hi > 0 ? log10(p_hi) : (p_hi == 0 ? -Inf : NaN))
    return (lo, hi)
end

"""
    _quality_panels(emitters, cfg) -> Dict{Symbol, NamedTuple}

The histogram panels of fit_quality.png as data (see `_histogram_panel`), keyed `:photons`,
`:precision`, `:pvalue`, `:psf`, `:z` and `:sigma_z` (each only when the emitters have the
field). Bounds are the ones `_filter_smld` applies: PSF widths come from `_resolve_psf_bounds`,
so a degenerate `:auto` that the filter skips draws nothing. The figure draws exactly these
spans, lines and counts, so a test of this function tests what the figure shows.
"""
function _quality_panels(emitters, cfg::FilterConfig)
    panels = Dict{Symbol, Any}()
    isempty(emitters) && return panels
    e1 = first(emitters)
    both = (:lo, :hi)

    photons = Float64[e.photons for e in emitters]
    panels[:photons] = _histogram_panel(
        (0.0, quantile(photons, 0.98)),
        [(; name = "photons", values = photons, bounds = cfg.photons, sides = both)];
        floor = 0.0
    )

    σ_x = Float64[e.σ_x for e in emitters] .* 1000
    σ_y = Float64[e.σ_y for e in emitters] .* 1000
    # Truncate at 1 μm before taking the window's quantile; with nothing below it, show 0-1 μm.
    prec_data = filter(<=(1000.0), vcat(σ_x, σ_y))
    prec98 = isempty(prec_data) ? 1000.0 : quantile(prec_data, 0.98)
    pb = _um_to_nm(cfg.precision)
    panels[:precision] = _histogram_panel(
        (0.0, prec98),
        [
            (; name = "σ_x", values = σ_x, bounds = pb, sides = (:hi,)),
            (; name = "σ_y", values = σ_y, bounds = pb, sides = (:hi,)),
        ];
        floor = 0.0
    )

    if _has_pvalue(emitters)
        p = Float64[e.pvalue for e in emitters]
        logp = log10.(filter(>=(0), p))  # p = 0 gives -Inf, below every view
        finite_logp = filter(isfinite, logp)
        pval_lo = isempty(finite_logp) ? -1.0 : quantile(finite_logp, 0.02)
        # p-values at (or rounding to) 1 put the 2nd percentile at 0; span one decade instead.
        pval_lo >= -1.0e-3 && (pval_lo = -1.0)
        panels[:pvalue] = _histogram_panel(
            (pval_lo, 0.0),
            [(; name = "log10(p)", values = logp, bounds = _log10_pvalue_bounds(cfg.pvalue), sides = both)];
            density = true
        )
    end

    psf = _resolve_psf_bounds(emitters, cfg)
    if hasproperty(e1, :σx) && hasproperty(e1, :σy)
        psf_σx = Float64[e.σx for e in emitters] .* 1000
        psf_σy = Float64[e.σy for e in emitters] .* 1000
        pooled = vcat(psf_σx, psf_σy)
        bx = psf === nothing ? nothing : (psf.lo_x * 1000, psf.hi_x * 1000)
        by = psf === nothing ? nothing : (psf.lo_y * 1000, psf.hi_y * 1000)
        panels[:psf] = _histogram_panel(
            (quantile(pooled, 0.02), quantile(pooled, 0.98)),
            [
                (; name = "σx", values = psf_σx, bounds = bx, sides = both),
                (; name = "σy", values = psf_σy, bounds = by, sides = both),
            ];
            floor = 0.0
        )
    elseif hasproperty(e1, :σ)
        psf_σ = Float64[e.σ for e in emitters] .* 1000
        b = psf === nothing ? nothing : (psf.lo * 1000, psf.hi * 1000)
        panels[:psf] = _histogram_panel(
            (quantile(psf_σ, 0.02), quantile(psf_σ, 0.98)),
            [(; name = "σ", values = psf_σ, bounds = b, sides = both)];
            floor = 0.0
        )
    end

    if hasproperty(e1, :z) && hasproperty(e1, :σ_z)
        z_nm = Float64[e.z for e in emitters] .* 1000
        # Keep the full grid extent (quantiles, not tight) so a spike at the z-grid edge shows.
        panels[:z] = _histogram_panel(
            (quantile(z_nm, 0.002), quantile(z_nm, 0.998)),
            [(; name = "z", values = z_nm, bounds = _um_to_nm(cfg.z), sides = both)]
        )
        σz_nm = Float64[e.σ_z for e in emitters] .* 1000
        panels[:sigma_z] = _histogram_panel(
            (0.0, quantile(σz_nm, 0.98)),
            [(; name = "σ_z", values = σz_nm, bounds = _um_to_nm(cfg.sigma_z), sides = both)];
            floor = 0.0
        )
    end
    return panels
end

"""
    _draw_histogram_panel!(ax, p, colors, alpha; rejected_color, threshold_color) -> (; plotted, beyond, ymax)

Draw a `_histogram_panel` on `ax`: shaded spans under the bars (gray for one series or when all
series share the same spans, else each series' own color, faint), bars from the panel's counts
(density for the p-value panel) in `colors[k]` at `alpha`, threshold lines on top (dotted; dashed
for a bound on an unshaded side), and the view as the x limits. `plotted` says whether any bar
was drawn, `beyond` how many values sit in end bars from beyond the view, `ymax` the tallest bar.
"""
function _draw_histogram_panel!(ax, p, colors, alpha; rejected_color, threshold_color)
    shared = length(p.series) == 1 || allequal([(s.spans, s.lines, s.dashed) for s in p.series])
    for (k, s) in enumerate(p.series)
        shared && k > 1 && break
        c = shared ? rejected_color : (colors[k], 0.2)
        for (a, b) in s.spans
            vspan!(ax, a, b, color = c)
        end
    end
    widths = diff(p.edges)
    centers = p.edges[1:(end - 1)] .+ widths ./ 2
    plotted = false
    ymax = 0.0
    for (k, s) in enumerate(p.series)
        n = sum(s.counts)
        n == 0 && continue
        y = p.density ? s.counts ./ (n .* widths) : Float64.(s.counts)
        barplot!(ax, centers, y, width = widths, gap = 0, color = (colors[k], alpha), label = s.name)
        plotted = true
        ymax = max(ymax, maximum(y))
    end
    for (k, s) in enumerate(p.series)
        shared && k > 1 && break
        lc = shared ? threshold_color : colors[k]
        isempty(s.lines) || vlines!(ax, s.lines, color = lc, linestyle = :dot, linewidth = 2)
        isempty(s.dashed) || vlines!(ax, s.dashed, color = lc, linestyle = :dash, linewidth = 2)
    end
    xlims!(ax, p.view...)
    return (; plotted, beyond = sum(s.beyond for s in p.series), ymax)
end

_beyond_note(n) = n == 0 ? "" : " ($n rejected beyond the axis, in the end bar)"

# P-values are optional on post-fit output; empty populations have no panel or criterion.
_has_pvalue(emitters::AbstractVector) = !isempty(emitters) && _has_pvalue(first(emitters))
_has_pvalue(emitter::AbstractEmitter) = hasproperty(emitter, :pvalue)

function _filter_smld(smld::BasicSMLD, cfg::FilterConfig)
    emitters = smld.emitters
    mask = trues(length(emitters))

    if cfg.photons !== nothing
        mask .&= [_in_range(e.photons, cfg.photons) for e in emitters]
    end

    if cfg.precision !== nothing
        mask .&= [_in_range(max(e.σ_x, e.σ_y), cfg.precision) for e in emitters]
    end

    if cfg.pvalue !== nothing
        if _has_pvalue(emitters)
            mask .&= [_in_range(e.pvalue, cfg.pvalue) for e in emitters]
        elseif !isempty(emitters)
            @warn "FilterConfig.pvalue has no effect: emitters have no p-value field (frame-connected, BaGoL and other output may lack p-values)"
        end
    end

    # 3D-only axial filters; no-op on 2D SMLDs (emitters without z / σ_z). Warn once
    # (not per-emitter) when the bound is set but has nothing to act on, rather than
    # silently doing nothing -- the empty-emitters case is skipped since we can't
    # tell 2D from 3D without an emitter to inspect.
    if cfg.z !== nothing && length(emitters) > 0
        if hasproperty(emitters[1], :z)
            lo, hi = cfg.z
            mask .&= [lo <= e.z <= hi for e in emitters]
        else
            @warn "FilterConfig.z has no effect: emitters have no z field (2D SMLD)"
        end
    end

    if cfg.sigma_z !== nothing && length(emitters) > 0
        if hasproperty(emitters[1], :σ_z)
            lo, hi = cfg.sigma_z
            mask .&= [lo <= e.σ_z <= hi for e in emitters]
        else
            @warn "FilterConfig.sigma_z has no effect: emitters have no σ_z field (2D SMLD)"
        end
    end

    if cfg.psf_sigma !== nothing && length(emitters) > 0
        # Determine bounds: :auto calculates mode ± 10%, or use explicit (min, max).
        # :auto returns the (0.0, 0.0) sentinel when mode estimation is degenerate
        # (e.g. all-zero σ) — skip rather than filter everything out. An explicit
        # (lo, hi), including a 0.0 lower bound, is a deliberate choice and is
        # always applied as given.
        is_auto = cfg.psf_sigma === :auto
        if hasproperty(emitters[1], :σ)
            lo, hi = _get_psf_sigma_bounds(cfg.psf_sigma, [e.σ for e in emitters])
            if !is_auto || (lo > 0 && hi > 0)
                mask .&= [lo <= e.σ <= hi for e in emitters]
            end
        elseif hasproperty(emitters[1], :σx) && hasproperty(emitters[1], :σy)
            lo_x, hi_x = _get_psf_sigma_bounds(cfg.psf_sigma, [e.σx for e in emitters])
            lo_y, hi_y = _get_psf_sigma_bounds(cfg.psf_sigma, [e.σy for e in emitters])
            if !is_auto || (lo_x > 0 && hi_x > 0 && lo_y > 0 && hi_y > 0)
                mask .&= [lo_x <= e.σx <= hi_x && lo_y <= e.σy <= hi_y for e in emitters]
            end
        else
            @warn "FilterConfig.psf_sigma has no effect: emitter type $(typeof(emitters[1])) has no PSF-width field"
        end
    end

    filtered = emitters[mask]
    return BasicSMLD(filtered, smld.camera, smld.n_frames, smld.n_datasets, smld.metadata)
end

"""
    _get_psf_sigma_bounds(range_spec, values) -> (lo, hi)

Calculate PSF sigma filter bounds.
- `:auto` → mode ± 10%
- `(min, max)` → explicit bounds in microns
"""
function _get_psf_sigma_bounds(range_spec, values::Vector)
    if range_spec === :auto
        mode = _calculate_mode(values)
        mode > 0 || return (0.0, 0.0)
        return (mode * 0.9, mode * 1.1)
    elseif range_spec isa Tuple{Float64, Float64}
        return range_spec
    else
        throw(ArgumentError("psf_sigma_range must be :auto or (min, max) tuple, got: $range_spec"))
    end
end

function _save_filter_outputs!(
        dir::String, outdir::Union{String, Nothing}, cfg::FilterConfig, v::Int, t::Float64,
        n_before::Int, n_after::Int,
        smld_input::BasicSMLD, smld_filtered::BasicSMLD
    )
    mkpath(dir)
    _save_config!(dir, cfg)

    if v >= Verbosity.STANDARD
        _write_filter_stats(dir, cfg, n_before, n_after, t, smld_input.emitters)
        _save_filter_quality_figures(dir, smld_input, cfg)
        _save_fit_overlay_from_cache(dir, outdir, smld_input, cfg)
        _save_loc_per_frame(dir, smld_filtered; title = "Localizations per Frame (post-filter)")
    end

    return if v >= Verbosity.DETAILED
        _save_filter_detailed(dir, smld_input, smld_filtered, cfg)
    end
end

function _write_filter_stats(dir, cfg, n_before, n_after, t, emitters)
    acceptance = n_before == 0 ? 0.0 : n_after / n_before

    filepath = joinpath(dir, "stats.md")
    return open(filepath, "w") do io
        println(io, "# Filter Statistics\n")
        println(io, "## Summary")
        println(io, "- **Input**: $n_before")
        println(io, "- **Output**: $n_after")
        println(io, "- **Acceptance**: $(round(100 * acceptance, digits = 1))%")
        println(io, "- **Time**: $(round(t * 1000, digits = 1))ms")
        println(io, "")
        println(io, "## Criteria Applied")
        if cfg.photons !== nothing
            lo, hi = cfg.photons
            println(io, "- photons: $(lo) - $(hi == Inf ? "∞" : hi)")
        end
        if cfg.precision !== nothing
            lo, hi = cfg.precision
            println(io, "- precision: $(round(lo * 1000, digits = 1)) - $(hi == Inf ? "∞" : round(hi * 1000, digits = 1)) nm")
        end
        if cfg.pvalue !== nothing && _has_pvalue(emitters)
            lo, hi = cfg.pvalue
            println(io, "- pvalue: $(lo) - $(hi)")
        end
        if cfg.psf_sigma !== nothing
            if cfg.psf_sigma === :auto
                println(io, "- psf_sigma: :auto (mode ± 10%)")
            else
                lo, hi = cfg.psf_sigma
                println(io, "- psf_sigma: $(round(lo * 1000, digits = 1)) - $(round(hi * 1000, digits = 1)) nm")
            end
        end
    end
end

function _save_filter_detailed(dir, smld_raw, smld_filtered, cfg)
    # Show which criteria rejected what
    emitters = smld_raw.emitters
    n = length(emitters)

    filepath = joinpath(dir, "detailed_stats.md")
    return open(filepath, "w") do io
        println(io, "# Filter Breakdown\n")
        println(io, "| Criterion | Pass | Fail | % Pass |")
        println(io, "|-----------|------|------|--------|")

        if cfg.photons !== nothing
            lo, hi = cfg.photons
            pass = sum(lo <= e.photons <= hi for e in emitters)
            hi_str = hi == Inf ? "∞" : string(hi)
            println(io, "| Photons ∈ [$lo, $hi_str] | $pass | $(n - pass) | $(round(100 * pass / n, digits = 1))% |")
        end

        if cfg.precision !== nothing
            lo, hi = cfg.precision
            pass = sum(lo <= max(e.σ_x, e.σ_y) <= hi for e in emitters)
            hi_str = hi == Inf ? "∞" : "$(round(hi * 1000, digits = 1))nm"
            println(io, "| Precision ∈ [$(round(lo * 1000, digits = 1))nm, $hi_str] | $pass | $(n - pass) | $(round(100 * pass / n, digits = 1))% |")
        end

        if cfg.pvalue !== nothing && _has_pvalue(emitters)
            lo, hi = cfg.pvalue
            pass = sum(lo <= e.pvalue <= hi for e in emitters)
            println(io, "| P-value ∈ [$lo, $hi] | $pass | $(n - pass) | $(round(100 * pass / n, digits = 1))% |")
        end
    end
end

"""
    _save_filter_quality_figures(dir, smld_raw, cfg)

Generate fit quality distribution plots with filter thresholds highlighted.
Shows photons, background, precision, p-value (when present), and PSF sigma distributions
from the raw (unfiltered) SMLD, with rejected regions marked per FilterConfig.
"""
function _save_filter_quality_figures(dir, smld_raw, cfg::FilterConfig)
    emitters = smld_raw.emitters
    isempty(emitters) && return

    # Every histogram panel is drawn from _quality_panels: its view, shaded spans, threshold
    # lines and bar counts, which the tests check against the filter.
    panels = _quality_panels(emitters, cfg)
    is_3d = haskey(panels, :z)

    REJECTED_COLOR = (:gray30, 0.5)
    MEAN_COLOR = :blue
    MEDIAN_COLOR = :red
    THRESHOLD_COLOR = :black
    draw!(ax, p, colors, alpha) = _draw_histogram_panel!(
        ax, p, colors, alpha; rejected_color = REJECTED_COLOR, threshold_color = THRESHOLD_COLOR
    )
    center_lines!(ax, v) = (
        vlines!(ax, [mean(v)], color = MEAN_COLOR, linestyle = :solid, linewidth = 2);
        vlines!(ax, [median(v)], color = MEDIAN_COLOR, linestyle = :dash, linewidth = 2)
    )

    fig = Figure(size = (1200, is_3d ? 1200 : 900))

    # Row 1: Photons and Background
    photons = [e.photons for e in emitters]
    ax1 = Axis(fig[1, 1], xlabel = "Photons", ylabel = "Count")
    d = draw!(ax1, panels[:photons], [:steelblue], 0.8)
    ax1.title = "Photon Distribution" * _beyond_note(d.beyond)
    center_lines!(ax1, photons)
    text!(
        ax1, 0.97, 0.95, text = "mean: $(round(Int, mean(photons)))\nmedian: $(round(Int, median(photons)))",
        align = (:right, :top), space = :relative, fontsize = 10
    )

    # The filter has no background criterion: a plain 0-98% histogram.
    bg = [e.bg for e in emitters]
    bg98 = quantile(bg, 0.98)
    ax2 = Axis(fig[1, 2], xlabel = "Background", ylabel = "Count", title = "Background Distribution")
    hist!(ax2, bg[bg .<= bg98], bins = 50)
    center_lines!(ax2, bg)
    xlims!(ax2, _plot_limits(0, bg98)...)
    text!(
        ax2, 0.97, 0.95, text = "mean: $(round(mean(bg), digits = 1))\nmedian: $(round(median(bg), digits = 1))",
        align = (:right, :top), space = :relative, fontsize = 10
    )

    # Row 2: Precision and P-value. The filter tests max(σx, σy): a component above the upper
    # bound is rejected, so that region is shaded; a component below the lower bound may still
    # pass (the other axis decides), so the lower bound is a dashed line with a note.
    σ_x = [e.σ_x for e in emitters] .* 1000
    σ_y = [e.σ_y for e in emitters] .* 1000
    ax3 = Axis(fig[2, 1], xlabel = "Localization Precision (nm)", ylabel = "Count")
    d = draw!(ax3, panels[:precision], [:blue, :red], 0.5)
    ax3.title = "Precision Distribution" * _beyond_note(d.beyond)
    d.plotted && axislegend(ax3, position = :rt, framevisible = false, labelsize = 9)
    prec_note = isempty(panels[:precision].series[1].dashed) ? "" : "\nlower bound on max(σ_x, σ_y)"
    text!(
        ax3, 0.97, 0.7, text = "σ_x: $(round(median(σ_x), digits = 1)) nm\nσ_y: $(round(median(σ_y), digits = 1)) nm$(prec_note)",
        align = (:right, :top), space = :relative, fontsize = 10
    )

    if haskey(panels, :pvalue)
        p = panels[:pvalue]
        pvalue = [e.pvalue for e in emitters]
        ax4 = Axis(fig[2, 2], xlabel = "log10(p-value)", ylabel = "Density")
        d = draw!(ax4, p, [:steelblue], 0.7)
        ax4.title = "P-value Distribution" * _beyond_note(d.beyond)
        # Uniform p-values are exponential in log10 space.
        u_hi = min(p.view[2], 0.0)
        if p.view[1] < u_hi
            u = range(p.view[1], u_hi, length = 100)
            lines!(ax4, u, log(10) .* (10.0 .^ u), color = :red, linewidth = 2, label = "Uniform theory")
        end
        logp = filter(isfinite, p.series[1].values)
        if !isempty(logp)
            vlines!(ax4, [mean(logp)], color = MEAN_COLOR, linestyle = :solid, linewidth = 1.5)
            vlines!(ax4, [median(logp)], color = MEDIAN_COLOR, linestyle = :dash, linewidth = 1.5)
        end
        if cfg.pvalue !== nothing
            bound_txt = cfg.pvalue[2] < 1 ? "range: $(cfg.pvalue)" : "threshold: $(cfg.pvalue[1])"
            text!(
                ax4, 0.03, 0.95, text = "pass: $(_pass_percent(pvalue, cfg.pvalue))%\n$(bound_txt)",
                align = (:left, :top), space = :relative, fontsize = 10
            )
        end
        ylims!(ax4, _plot_limits(0, d.ymax * 1.1)...)
    end

    # Row 3: PSF Sigma
    if haskey(panels, :psf)
        p = panels[:psf]
        ax5 = Axis(fig[3, 1], xlabel = "Fitted PSF σ (nm)", ylabel = "Count")
        if length(p.series) == 2
            d = draw!(ax5, p, [:blue, :red], 0.5)
            d.plotted && axislegend(ax5, position = :rt, framevisible = false, labelsize = 9)
            mode_x = _calculate_mode([e.σx for e in emitters]) * 1000
            mode_y = _calculate_mode([e.σy for e in emitters]) * 1000
            txt = "mode σx: $(round(mode_x, digits = 1)) nm\nmode σy: $(round(mode_y, digits = 1)) nm"
            text!(ax5, 0.97, 0.7, text = txt, align = (:right, :top), space = :relative, fontsize = 10)
        else
            d = draw!(ax5, p, [:steelblue], 0.8)
            psf_σ = p.series[1].values
            center_lines!(ax5, psf_σ)
            mode_σ = _calculate_mode([e.σ for e in emitters]) * 1000
            txt = "mode: $(round(mode_σ, digits = 1)) nm\nmean: $(round(mean(psf_σ), digits = 1)) nm"
            text!(ax5, 0.97, 0.95, text = txt, align = (:right, :top), space = :relative, fontsize = 10)
        end
        ax5.title = "PSF Width Distribution" * _beyond_note(d.beyond)
    else
        ax5 = Axis(fig[3, 1], xlabel = "PSF σ (nm)", ylabel = "", title = "PSF Width (Fixed)")
        text!(
            ax5, 0.5, 0.5, text = "Fixed PSF model",
            align = (:center, :center), space = :relative, fontsize = 14
        )
        hideydecorations!(ax5)
    end

    # Row 3, Col 2: Legend
    ax6 = Axis(fig[3, 2], title = "Legend")
    hidedecorations!(ax6)
    hidespines!(ax6)
    lines!(ax6, [0.1, 0.25], [0.8, 0.8], color = MEAN_COLOR, linewidth = 2)
    text!(ax6, 0.3, 0.8, text = "Mean", fontsize = 12)
    lines!(ax6, [0.1, 0.25], [0.6, 0.6], color = MEDIAN_COLOR, linewidth = 2, linestyle = :dash)
    text!(ax6, 0.3, 0.6, text = "Median", fontsize = 12)
    lines!(ax6, [0.1, 0.25], [0.4, 0.4], color = THRESHOLD_COLOR, linewidth = 2, linestyle = :dot)
    text!(ax6, 0.3, 0.4, text = "Filter Threshold", fontsize = 12)
    poly!(ax6, Point2f[(0.1, 0.15), (0.25, 0.15), (0.25, 0.25), (0.1, 0.25)], color = REJECTED_COLOR)
    text!(ax6, 0.3, 0.2, text = "Rejected Region", fontsize = 12)
    xlims!(ax6, 0, 1)
    ylims!(ax6, 0, 1)

    # Row 4 (3D only): axial position z and axial precision σ_z. The z panel makes
    # z-grid-edge-railed degenerate localizations obvious as a spike at the boundary.
    if is_3d
        z_nm = panels[:z].series[1].values
        ax7 = Axis(fig[4, 1], xlabel = "Axial position z (nm)", ylabel = "Count")
        d = draw!(ax7, panels[:z], [:purple], 0.6)
        ax7.title = "z Distribution" * _beyond_note(d.beyond)
        center_lines!(ax7, z_nm)
        text!(
            ax7, 0.97, 0.95, text = "median: $(round(median(z_nm), digits = 1)) nm\nrange: [$(round(minimum(z_nm), digits = 0)), $(round(maximum(z_nm), digits = 0))] nm",
            align = (:right, :top), space = :relative, fontsize = 10
        )

        σz_nm = panels[:sigma_z].series[1].values
        ax8 = Axis(fig[4, 2], xlabel = "Axial precision σ_z (nm)", ylabel = "Count")
        d = draw!(ax8, panels[:sigma_z], [:green], 0.5)
        ax8.title = "Axial Precision Distribution" * _beyond_note(d.beyond)
        center_lines!(ax8, σz_nm)
        text!(
            ax8, 0.97, 0.95, text = "median: $(round(median(σz_nm), digits = 1)) nm\nmean: $(round(mean(σz_nm), digits = 1)) nm",
            align = (:right, :top), space = :relative, fontsize = 10
        )
    end

    return save(joinpath(dir, "fit_quality.png"), fig)
end


# ============================================================
# Fit overlay from detectfit cache (colored by FilterConfig thresholds)
# ============================================================

"""
Load detectfit sample cache and generate fit_overlay.png with FilterConfig thresholds.
Gracefully skips if cache is missing (e.g., standalone filter without prior detectfit).
"""
function _save_fit_overlay_from_cache(dir, outdir, smld_raw, cfg::FilterConfig)
    cache = load_cache(outdir, "detectfit_samples.jld2")
    cache === nothing && return

    # Reconstruct ROIBatch from cached components using smld_raw's camera
    roi_batch = ROIBatch(
        cache["sample_roi_data"],
        cache["sample_roi_x"],
        cache["sample_roi_y"],
        cache["sample_roi_frames"],
        smld_raw.camera,
    )

    return _save_fit_overlay(
        dir, smld_raw, roi_batch,
        cache["sample_images"], cache["sample_original_frames"],
        cache["n_frames"], cache["n_datasets"], cfg
    )
end

"""
Generate fit overlay: boxes colored by fit quality using FilterConfig thresholds.

Colors: green=pass, red=photons fail, orange=precision fail, purple=pvalue fail, gray=no match.
"""
function _save_fit_overlay(
        dir, smld, sample_roi_batch, sample_images, sample_original_frames,
        n_frames, n_datasets, cfg::FilterConfig
    )
    isempty(sample_roi_batch) && return

    # Map sample index (1:N) to absolute frame number
    sample_to_original = Dict(i => f for (i, f) in enumerate(sample_original_frames))
    original_frame_set = Set(sample_original_frames)

    # Build absolute frame for each emitter and filter to sampled frames
    emitter_abs_frames = [(e, (e.dataset - 1) * n_frames + e.frame) for e in smld.emitters]
    sample_emitters = [(e, af) for (e, af) in emitter_abs_frames if af in original_frame_set]
    isempty(sample_emitters) && return

    # Pre-resolve PSF sigma bounds (needed for box coloring)
    psf_bounds = _resolve_psf_bounds(smld.emitters, cfg)

    # Color each ROI box by matching emitter quality
    fit_colors = Symbol[]
    pix_size = smld.camera.pixel_edges_x[2] - smld.camera.pixel_edges_x[1]
    x_origin = smld.camera.pixel_edges_x[1]
    y_origin = smld.camera.pixel_edges_y[1]

    for i in 1:length(sample_roi_batch)
        roi_frame_idx = sample_roi_batch.frame_indices[i]
        original_frame = sample_to_original[roi_frame_idx]
        roi_x = sample_roi_batch.x_corners[i] + sample_roi_batch.roi_size ÷ 2
        roi_y = sample_roi_batch.y_corners[i] + sample_roi_batch.roi_size ÷ 2

        # Find closest matching emitter by position in this frame
        best_emitter = nothing
        best_dist = Inf
        for (e, af) in sample_emitters
            if af == original_frame
                ex_px = (e.x - x_origin) / pix_size
                ey_px = (e.y - y_origin) / pix_size
                dist = sqrt((ex_px - roi_x)^2 + (ey_px - roi_y)^2)
                if dist < best_dist
                    best_dist = dist
                    best_emitter = e
                end
            end
        end

        if best_emitter !== nothing && best_dist < sample_roi_batch.roi_size
            push!(fit_colors, _fit_box_color(best_emitter, cfg, psf_bounds))
        else
            push!(fit_colors, :gray)
        end
    end

    title = _fit_overlay_title(cfg, smld.emitters)
    return _save_box_overlay(
        dir, "fit_overlay.png", sample_images, sample_roi_batch, fit_colors;
        title_prefix = "Frame", frame_labels = sample_original_frames, suptitle = title
    )
end

"""
    _resolve_psf_bounds(emitters, cfg) -> Union{NamedTuple, Nothing}

Pre-resolve PSF sigma filter bounds from the emitter population.
Returns named tuple with :iso or :aniso bounds, or nothing if psf_sigma filter is off.
"""
function _resolve_psf_bounds(emitters, cfg::FilterConfig)
    cfg.psf_sigma === nothing && return nothing
    isempty(emitters) && return nothing

    # Same rule as _filter_smld: the positive-bound sentinel (no mode found) disables only
    # :auto; explicit bounds, including a zero lower bound, always apply.
    is_auto = cfg.psf_sigma === :auto
    if hasproperty(emitters[1], :σ)
        lo, hi = _get_psf_sigma_bounds(cfg.psf_sigma, [e.σ for e in emitters])
        !is_auto || (lo > 0 && hi > 0) || return nothing
        return (kind = :iso, lo = lo, hi = hi)
    elseif hasproperty(emitters[1], :σx) && hasproperty(emitters[1], :σy)
        lo_x, hi_x = _get_psf_sigma_bounds(cfg.psf_sigma, [e.σx for e in emitters])
        lo_y, hi_y = _get_psf_sigma_bounds(cfg.psf_sigma, [e.σy for e in emitters])
        !is_auto || (lo_x > 0 && hi_x > 0 && lo_y > 0 && hi_y > 0) || return nothing
        return (kind = :aniso, lo_x = lo_x, hi_x = hi_x, lo_y = lo_y, hi_y = hi_y)
    end
    return nothing
end

"""Determine box color based on FilterConfig thresholds. All-nothing filters = all green."""
function _fit_box_color(e, cfg::FilterConfig, psf_bounds)
    if cfg.photons !== nothing
        _in_range(e.photons, cfg.photons) || return :red
    end
    if cfg.precision !== nothing
        _in_range(max(e.σ_x, e.σ_y), cfg.precision) || return :orange
    end
    if cfg.pvalue !== nothing && _has_pvalue(e)
        _in_range(e.pvalue, cfg.pvalue) || return :purple
    end
    if cfg.z !== nothing && hasproperty(e, :z)
        _in_range(e.z, cfg.z) || return :blue
    end
    if cfg.sigma_z !== nothing && hasproperty(e, :σ_z)
        _in_range(e.σ_z, cfg.sigma_z) || return :magenta
    end
    if psf_bounds !== nothing
        if psf_bounds.kind === :iso
            (psf_bounds.lo <= e.σ <= psf_bounds.hi) || return :cyan
        elseif psf_bounds.kind === :aniso
            (psf_bounds.lo_x <= e.σx <= psf_bounds.hi_x && psf_bounds.lo_y <= e.σy <= psf_bounds.hi_y) || return :cyan
        end
    end
    return :green
end

# A photon bound as a label: whole photons when it fits in an Int, else 3 significant digits
# (round(Int, 1e20) would throw, failing the filter step while it writes figures).
_photon_label(x) = abs(x) < 1.0e15 ? string(round(Int, x)) : string(round(x, sigdigits = 3))

"""Build title string from actual FilterConfig thresholds."""
function _fit_overlay_title(cfg::FilterConfig, emitters)
    parts = String["green=pass"]
    # Open (infinite) bounds reject nothing and get no label.
    if cfg.photons !== nothing
        lo, hi = cfg.photons
        conds = String[]
        isfinite(lo) && push!(conds, "<$(_photon_label(lo))")
        isfinite(hi) && push!(conds, ">$(_photon_label(hi))")
        isempty(conds) || push!(parts, "red=photons" * join(conds, "/"))
    end
    if cfg.precision !== nothing
        lo, hi = cfg.precision
        conds = String[]
        # Precision is never negative, so a lower bound of 0 rejects nothing either.
        isfinite(lo) && lo > 0 && push!(conds, "<$(round(lo * 1000, digits = 1))nm")
        isfinite(hi) && push!(conds, ">$(round(hi * 1000, digits = 1))nm")
        isempty(conds) || push!(parts, "orange=prec" * join(conds, "/"))
    end
    if cfg.pvalue !== nothing && _has_pvalue(emitters)
        # The filter rejects any p-value that fails its range test, NaN included (DetectFit
        # keeps fits whose p-value alone is NaN), so the label names the range, not an inequality.
        push!(parts, "purple=pval not in [$(cfg.pvalue[1]), $(cfg.pvalue[2])]")
    end
    if cfg.z !== nothing && !isempty(emitters) && hasproperty(first(emitters), :z)
        push!(parts, "blue=z outside $(cfg.z)μm")
    end
    if cfg.sigma_z !== nothing && !isempty(emitters) && hasproperty(first(emitters), :σ_z)
        push!(parts, "magenta=σ_z outside $(cfg.sigma_z)μm")
    end
    # The PSF bounds the filter applies (none when :auto finds no mode and the filter skips it).
    psf = _resolve_psf_bounds(emitters, cfg)
    nm(lo, hi) = "[$(round(lo * 1000, digits = 1)), $(round(hi * 1000, digits = 1))]nm"
    if psf !== nothing && psf.kind === :iso
        push!(parts, "cyan=psf_σ not in $(nm(psf.lo, psf.hi))")
    elseif psf !== nothing
        push!(parts, "cyan=psf_σx not in $(nm(psf.lo_x, psf.hi_x)) or σy not in $(nm(psf.lo_y, psf.hi_y))")
    end
    push!(parts, "gray=no match")
    return "Fit: " * join(parts, "  ")
end
