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

"""
    _rejected_spans(lo, hi, view_lo, view_hi) -> (; spans, thresholds)

What a filter range `(lo, hi)` rejects inside a figure's x window `[view_lo, view_hi]`:
`[view_lo, lo]` when `lo` is finite and above `view_lo`, and `[hi, view_hi]` when `hi` is
finite and below `view_hi`, with a threshold line at each such bound. An open (infinite)
bound, or one at the window's edge, rejects nothing shown and gets neither.
"""
function _rejected_spans(lo, hi, view_lo, view_hi)
    spans = Tuple{Float64, Float64}[]
    thresholds = Float64[]
    if isfinite(lo) && lo > view_lo
        push!(spans, (Float64(view_lo), Float64(lo)))
        push!(thresholds, Float64(lo))
    end
    if isfinite(hi) && hi < view_hi
        push!(spans, (Float64(hi), Float64(view_hi)))
        push!(thresholds, Float64(hi))
    end
    return (; spans, thresholds)
end

function _draw_rejected!(ax, r, region_color, line_color)
    for (a, b) in r.spans
        vspan!(ax, a, b, color = region_color)
    end
    isempty(r.thresholds) || vlines!(ax, r.thresholds, color = line_color, linestyle = :dot, linewidth = 2)
    return nothing
end

# Percentage of `values` the filter keeps for `range`, by the same test the filter uses.
_pass_percent(values, range) = round(100 * count(v -> _in_range(v, range), values) / length(values), digits = 1)

"""
    _bin_counts(values, edges) -> Vector{Int}

Count `values` into the bins `[edges[i], edges[i+1])` by comparing with the stored edge values
(dividing by the bin width can put a value that sits on an edge into the bin below). The top
edge belongs to the last bin; values outside `[first(edges), last(edges)]` are not counted.
"""
function _bin_counts(values, edges::AbstractVector)
    nb = length(edges) - 1
    counts = zeros(Int, nb)
    for v in values
        first(edges) <= v <= last(edges) || continue
        counts[min(searchsortedlast(edges, v), nb)] += 1
    end
    return counts
end

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

    photons = [e.photons for e in emitters]
    bg = [e.bg for e in emitters]
    σ_x = [e.σ_x for e in emitters] .* 1000  # precision in nm
    σ_y = [e.σ_y for e in emitters] .* 1000

    # Check PSF model type
    has_psf_iso = hasproperty(emitters[1], :σ)
    has_psf_aniso = hasproperty(emitters[1], :σx)
    # 3D fits get an extra row of z-position + axial-precision panels.
    is_3d = hasproperty(emitters[1], :z) && hasproperty(emitters[1], :σ_z)

    # Colors for consistent styling
    REJECTED_COLOR = (:gray30, 0.5)
    MEAN_COLOR = :blue
    MEDIAN_COLOR = :red
    THRESHOLD_COLOR = :black

    fig = Figure(size = (1200, is_3d ? 1200 : 900))

    # Row 1: Photons and Background
    p98 = quantile(photons, 0.98)
    photons_xlim = p98
    ax1 = Axis(fig[1, 1], xlabel = "Photons", ylabel = "Count", title = "Photon Distribution")
    if cfg.photons !== nothing
        photons_lo, photons_hi = cfg.photons
        # A finite upper bound widens the view so its rejected region shows; open
        # (infinite) bounds reject nothing and draw nothing.
        isfinite(photons_hi) && (photons_xlim = max(p98, photons_hi * 1.05))
        _draw_rejected!(ax1, _rejected_spans(photons_lo, photons_hi, 0, photons_xlim), REJECTED_COLOR, THRESHOLD_COLOR)
    end
    hist!(ax1, photons[photons .<= p98], bins = 50)
    vlines!(ax1, [mean(photons)], color = MEAN_COLOR, linestyle = :solid, linewidth = 2)
    vlines!(ax1, [median(photons)], color = MEDIAN_COLOR, linestyle = :dash, linewidth = 2)
    xlims!(ax1, _plot_limits(0, photons_xlim)...)
    text!(
        ax1, 0.97, 0.95, text = "mean: $(round(Int, mean(photons)))\nmedian: $(round(Int, median(photons)))",
        align = (:right, :top), space = :relative, fontsize = 10
    )

    bg98 = quantile(bg, 0.98)
    ax2 = Axis(fig[1, 2], xlabel = "Background", ylabel = "Count", title = "Background Distribution")
    hist!(ax2, bg[bg .<= bg98], bins = 50)
    vlines!(ax2, [mean(bg)], color = MEAN_COLOR, linestyle = :solid, linewidth = 2)
    vlines!(ax2, [median(bg)], color = MEDIAN_COLOR, linestyle = :dash, linewidth = 2)
    xlims!(ax2, _plot_limits(0, bg98)...)
    text!(
        ax2, 0.97, 0.95, text = "mean: $(round(mean(bg), digits = 1))\nmedian: $(round(median(bg), digits = 1))",
        align = (:right, :top), space = :relative, fontsize = 10
    )

    # Row 2: Precision and P-value
    prec_data = vcat(σ_x, σ_y)
    prec_data = prec_data[prec_data .<= 1000.0]  # truncate at 1 μm before percentile stretch
    # Every precision above the cutoff leaves nothing to stretch: show the empty 0-1 μm window.
    prec98 = isempty(prec_data) ? 1000.0 : quantile(prec_data, 0.98)
    prec_xlim = prec98

    ax3 = Axis(fig[2, 1], xlabel = "Localization Precision (nm)", ylabel = "Count", title = "Precision Distribution")
    if cfg.precision !== nothing
        prec_lo, prec_hi = cfg.precision .* 1000  # convert um to nm
        # Finite bounds widen the view so their rejected regions show; open bounds draw nothing.
        isfinite(prec_hi) && (prec_xlim = max(prec98, prec_hi * 1.2))
        isfinite(prec_lo) && prec_lo > 0 && (prec_xlim = max(prec_xlim, prec_lo * 1.2))
        _draw_rejected!(ax3, _rejected_spans(prec_lo, prec_hi, 0, prec_xlim), REJECTED_COLOR, THRESHOLD_COLOR)
    end
    σx_view = σ_x[σ_x .<= prec_xlim]
    σy_view = σ_y[σ_y .<= prec_xlim]
    !isempty(σx_view) && hist!(ax3, σx_view, bins = 50, color = (:blue, 0.5), label = "σ_x")
    !isempty(σy_view) && hist!(ax3, σy_view, bins = 50, color = (:red, 0.5), label = "σ_y")
    xlims!(ax3, _plot_limits(0, prec_xlim)...)
    (isempty(σx_view) && isempty(σy_view)) || axislegend(ax3, position = :rt, framevisible = false, labelsize = 9)
    text!(
        ax3, 0.97, 0.7, text = "σ_x: $(round(median(σ_x), digits = 1)) nm\nσ_y: $(round(median(σ_y), digits = 1)) nm",
        align = (:right, :top), space = :relative, fontsize = 10
    )

    if _has_pvalue(emitters)
        pvalue = [e.pvalue for e in emitters]
        ax4 = Axis(fig[2, 2], xlabel = "log10(p-value)", ylabel = "Density", title = "P-value Distribution")
        pval_nonzero = pvalue[pvalue .> 0]
        pval_thresh = cfg.pvalue !== nothing ? cfg.pvalue[1] : nothing

        if !isempty(pval_nonzero)
            log_pval = log10.(pval_nonzero)
            pval_lo = quantile(log_pval, 0.02)
            # p-values at (or rounding to) 1 put pval_lo at 0, which would collapse every
            # edge to 0; span one decade instead. The same edges draw the histogram, size the
            # y axis and set the x limits.
            if pval_lo >= -1.0e-3
                pval_lo = -1.0
            end
            log_pval_filtered = log_pval[pval_lo .<= log_pval .<= 0]
            nbins = 50
            bin_edges = range(pval_lo, 0, length = nbins + 1)
            bin_width = step(bin_edges)
            if cfg.pvalue !== nothing
                # In log10 space: p = 0 and p = 1 bounds reject nothing.
                p_lo, p_hi = cfg.pvalue
                lo_log = p_lo > 0 ? log10(p_lo) : -Inf
                hi_log = p_hi < 1 ? log10(p_hi) : Inf
                _draw_rejected!(ax4, _rejected_spans(lo_log, hi_log, pval_lo - 1, 0.0), REJECTED_COLOR, THRESHOLD_COLOR)
            end
            # Count once and draw the bars from these counts, so the bars and the y limit
            # come from the same numbers (a separate hist! can bin a value on an edge
            # differently and clip its peak). p = 1 (log 0, the top edge) goes in the last bin.
            edges = collect(bin_edges)
            counts = _bin_counts(log_pval_filtered, edges)
            density = counts ./ (length(log_pval_filtered) * bin_width)
            bin_centers = (edges[1:nbins] .+ edges[2:end]) ./ 2
            barplot!(ax4, bin_centers, density, width = bin_width, gap = 0, color = (:steelblue, 0.7))
            max_density = maximum(density)
            # Theory curve (uniform p-values -> exponential in log space)
            u_range = range(0, -pval_lo, length = 100)
            theory_pdf = log(10) .* (10.0 .^ (-u_range))
            lines!(ax4, -u_range, theory_pdf, color = :red, linewidth = 2, label = "Uniform theory")
            vlines!(ax4, [mean(log_pval)], color = MEAN_COLOR, linestyle = :solid, linewidth = 1.5)
            vlines!(ax4, [median(log_pval)], color = MEDIAN_COLOR, linestyle = :dash, linewidth = 1.5)
            if pval_thresh !== nothing
                pval_pass_pct = _pass_percent(pvalue, cfg.pvalue)
                bound_txt = cfg.pvalue[2] < 1 ? "range: $(cfg.pvalue)" : "threshold: $(pval_thresh)"
                text!(
                    ax4, 0.03, 0.95, text = "pass: $(pval_pass_pct)%\n$(bound_txt)",
                    align = (:left, :top), space = :relative, fontsize = 10
                )
            end
            xlims!(ax4, _plot_limits(pval_lo, 0)...)
            ylims!(ax4, _plot_limits(0, max_density * 1.1)...)
        end
    end

    # Row 3: PSF Sigma
    if has_psf_aniso
        psf_σx = [e.σx for e in emitters] .* 1000
        psf_σy = [e.σy for e in emitters] .* 1000
        psf_data = vcat(psf_σx, psf_σy)
        psf98 = quantile(psf_data, 0.98)
        psf02 = quantile(psf_data, 0.02)
        mode_x = _calculate_mode([e.σx for e in emitters]) * 1000
        mode_y = _calculate_mode([e.σy for e in emitters]) * 1000

        # Compute per-axis bounds matching _filter_smld logic
        if cfg.psf_sigma !== nothing
            if cfg.psf_sigma === :auto
                tol = 0.1
                lo_x_nm, hi_x_nm = mode_x * (1 - tol), mode_x * (1 + tol)
                lo_y_nm, hi_y_nm = mode_y * (1 - tol), mode_y * (1 + tol)
            else
                lo_x_nm, hi_x_nm = cfg.psf_sigma[1] * 1000, cfg.psf_sigma[2] * 1000
                lo_y_nm, hi_y_nm = lo_x_nm, hi_x_nm
            end
        end

        ax5 = Axis(fig[3, 1], xlabel = "Fitted PSF σ (nm)", ylabel = "Count", title = "PSF Width Distribution")
        # Show per-axis filter bounds (matching actual filter behavior)
        if cfg.psf_sigma !== nothing
            # Open (infinite) bounds reject nothing and have no line to draw; Makie
            # cannot take limits from an empty vlines!.
            bounds_x = filter(isfinite, [lo_x_nm, hi_x_nm])
            bounds_y = filter(isfinite, [lo_y_nm, hi_y_nm])
            isempty(bounds_x) || vlines!(ax5, bounds_x, color = :blue, linestyle = :dot, linewidth = 2, label = "σx bounds")
            isempty(bounds_y) || vlines!(ax5, bounds_y, color = :red, linestyle = :dot, linewidth = 2, label = "σy bounds")
        end
        # The 2-98% window is pooled over both axes, so one axis (or both, for a few
        # emitters) can fall entirely outside it; Makie before 0.24.17 throws on an empty hist!.
        psf_σx_view = psf_σx[(psf_σx .>= psf02) .& (psf_σx .<= psf98)]
        psf_σy_view = psf_σy[(psf_σy .>= psf02) .& (psf_σy .<= psf98)]
        isempty(psf_σx_view) || hist!(ax5, psf_σx_view, bins = 50, color = (:blue, 0.5), label = "σx")
        isempty(psf_σy_view) || hist!(ax5, psf_σy_view, bins = 50, color = (:red, 0.5), label = "σy")
        xlims!(ax5, _plot_limits(psf02, psf98)...)
        (isempty(psf_σx_view) && isempty(psf_σy_view)) || axislegend(ax5, position = :rt, framevisible = false, labelsize = 9)
        text!(
            ax5, 0.97, 0.7, text = "mode σx: $(round(mode_x, digits = 1)) nm\nmode σy: $(round(mode_y, digits = 1)) nm",
            align = (:right, :top), space = :relative, fontsize = 10
        )

    elseif has_psf_iso
        psf_σ = [e.σ for e in emitters] .* 1000
        psf98 = quantile(psf_σ, 0.98)
        psf02 = quantile(psf_σ, 0.02)
        mode_σ = _calculate_mode([e.σ for e in emitters]) * 1000

        ax5 = Axis(fig[3, 1], xlabel = "Fitted PSF σ (nm)", ylabel = "Count", title = "PSF Width Distribution")
        if cfg.psf_sigma !== nothing
            if cfg.psf_sigma === :auto
                tol = 0.1
                lo_bound = mode_σ * (1 - tol)
                hi_bound = mode_σ * (1 + tol)
            else
                lo_bound = cfg.psf_sigma[1] * 1000
                hi_bound = cfg.psf_sigma[2] * 1000
            end
            # The view extends to finite bounds only; an open bound rejects nothing on its side.
            psf_xmin = isfinite(lo_bound) ? min(psf02, lo_bound * 0.95) : psf02
            psf_xmax = isfinite(hi_bound) ? max(psf98, hi_bound * 1.05) : psf98
            isfinite(lo_bound) && vspan!(ax5, psf_xmin, lo_bound, color = REJECTED_COLOR)
            isfinite(hi_bound) && vspan!(ax5, hi_bound, psf_xmax, color = REJECTED_COLOR)
            bounds = filter(isfinite, [lo_bound, hi_bound])
            isempty(bounds) || vlines!(ax5, bounds, color = THRESHOLD_COLOR, linestyle = :dot, linewidth = 2)
        else
            psf_xmin = psf02
            psf_xmax = psf98
        end
        # Two distinct values both fall outside their own 2-98% window.
        psf_σ_view = psf_σ[(psf_σ .>= psf02) .& (psf_σ .<= psf98)]
        isempty(psf_σ_view) || hist!(ax5, psf_σ_view, bins = 50)
        vlines!(ax5, [mean(psf_σ)], color = MEAN_COLOR, linestyle = :solid, linewidth = 2)
        vlines!(ax5, [median(psf_σ)], color = MEDIAN_COLOR, linestyle = :dash, linewidth = 2)
        xlims!(ax5, _plot_limits(psf_xmin, psf_xmax)...)
        text!(
            ax5, 0.97, 0.95, text = "mode: $(round(mode_σ, digits = 1)) nm\nmean: $(round(mean(psf_σ), digits = 1)) nm",
            align = (:right, :top), space = :relative, fontsize = 10
        )
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
        z_nm = [e.z for e in emitters] .* 1000      # axial position in nm
        σz_nm = [e.σ_z for e in emitters] .* 1000   # axial precision in nm

        # --- z position panel ---
        # Keep the full grid extent visible (quantiles, not tight) so a boundary
        # spike is not clipped; extend to include any filter bounds.
        z_lo_view = quantile(z_nm, 0.002)
        z_hi_view = quantile(z_nm, 0.998)
        ax7 = Axis(fig[4, 1], xlabel = "Axial position z (nm)", ylabel = "Count", title = "z Distribution")
        if cfg.z !== nothing
            z_lo, z_hi = cfg.z[1] * 1000, cfg.z[2] * 1000
            # Only finite bounds widen the view or draw a rejected region.
            isfinite(z_lo) && (z_lo_view = min(z_lo_view, z_lo * 1.1))
            isfinite(z_hi) && (z_hi_view = max(z_hi_view, z_hi * 1.1))
            isfinite(z_lo) && vspan!(ax7, z_lo_view, z_lo, color = REJECTED_COLOR)
            isfinite(z_hi) && vspan!(ax7, z_hi, z_hi_view, color = REJECTED_COLOR)
            z_bounds = filter(isfinite, [z_lo, z_hi])
            isempty(z_bounds) || vlines!(ax7, z_bounds, color = THRESHOLD_COLOR, linestyle = :dot, linewidth = 2)
        end
        z_view = z_nm[(z_nm .>= z_lo_view) .& (z_nm .<= z_hi_view)]
        !isempty(z_view) && hist!(ax7, z_view, bins = 50, color = (:purple, 0.6))
        vlines!(ax7, [mean(z_nm)], color = MEAN_COLOR, linestyle = :solid, linewidth = 2)
        vlines!(ax7, [median(z_nm)], color = MEDIAN_COLOR, linestyle = :dash, linewidth = 2)
        xlims!(ax7, _plot_limits(z_lo_view, z_hi_view)...)
        text!(
            ax7, 0.97, 0.95, text = "median: $(round(median(z_nm), digits = 1)) nm\nrange: [$(round(minimum(z_nm), digits = 0)), $(round(maximum(z_nm), digits = 0))] nm",
            align = (:right, :top), space = :relative, fontsize = 10
        )

        # --- σ_z (axial precision) panel ---
        σz98 = quantile(σz_nm, 0.98)
        σz_xlim = σz98
        ax8 = Axis(fig[4, 2], xlabel = "Axial precision σ_z (nm)", ylabel = "Count", title = "Axial Precision Distribution")
        if cfg.sigma_z !== nothing
            σz_lo, σz_hi = cfg.sigma_z .* 1000  # convert μm to nm
            # Finite bounds widen the view so their rejected regions show; open bounds draw nothing.
            isfinite(σz_hi) && (σz_xlim = max(σz98, σz_hi * 1.2))
            isfinite(σz_lo) && σz_lo > 0 && (σz_xlim = max(σz_xlim, σz_lo * 1.2))
            _draw_rejected!(ax8, _rejected_spans(σz_lo, σz_hi, 0, σz_xlim), REJECTED_COLOR, THRESHOLD_COLOR)
        end
        σz_view = σz_nm[σz_nm .<= σz_xlim]
        !isempty(σz_view) && hist!(ax8, σz_view, bins = 50, color = (:green, 0.5))
        vlines!(ax8, [mean(σz_nm)], color = MEAN_COLOR, linestyle = :solid, linewidth = 2)
        vlines!(ax8, [median(σz_nm)], color = MEDIAN_COLOR, linestyle = :dash, linewidth = 2)
        xlims!(ax8, _plot_limits(0, σz_xlim)...)
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

"""Build title string from actual FilterConfig thresholds."""
function _fit_overlay_title(cfg::FilterConfig, emitters)
    parts = String["green=pass"]
    # Open (infinite) bounds reject nothing and get no label; round(Int, ±Inf) would throw.
    if cfg.photons !== nothing
        lo, hi = cfg.photons
        conds = String[]
        isfinite(lo) && push!(conds, "<$(round(Int, lo))")
        isfinite(hi) && push!(conds, ">$(round(Int, hi))")
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
        lo, hi = cfg.pvalue
        conds = String[]
        # p-values lie in [0, 1], so bounds of 0 and 1 reject nothing.
        lo > 0 && push!(conds, "<$(lo)")
        hi < 1 && push!(conds, ">$(hi)")
        isempty(conds) || push!(parts, "purple=pval" * join(conds, "/"))
    end
    if cfg.z !== nothing && !isempty(emitters) && hasproperty(first(emitters), :z)
        push!(parts, "blue=z outside $(cfg.z)μm")
    end
    if cfg.sigma_z !== nothing && !isempty(emitters) && hasproperty(first(emitters), :σ_z)
        push!(parts, "magenta=σ_z outside $(cfg.sigma_z)μm")
    end
    if cfg.psf_sigma !== nothing
        push!(parts, "cyan=psf_σ fail")
    end
    push!(parts, "gray=no match")
    return "Fit: " * join(parts, "  ")
end
