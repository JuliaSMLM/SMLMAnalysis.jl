# Shared checks of the filter figures, included by colorrange.jl (Core) and
# long/filter_figures.jl (Long). Run on its own, as a Core test file, it only defines them.
using SMLMAnalysis
using GaussMLE
using Random
using Test

# Every axis of fit_quality.png that can collapse or go infinite. Each case must save the
# figure; the filter result must be the same with and without the figures.
function _quality_figure_saves(smld, cfg)
    return mktempdir() do dir
        SMLMAnalysis._save_filter_quality_figures(dir, smld, cfg)
        isfile(joinpath(dir, "fit_quality.png"))
    end
end

_sigma_emitter(; σ = 0.13, bg = 5.0, prec = 0.01, pvalue = 0.4, frame = 1, photons = 1000.0) =
    GaussMLE.Emitter2DFitSigma{Float64}(
    0.8, 0.8, photons, bg, σ, prec, prec, 0.0, 30.0, 0.5, 0.002, pvalue, frame, 1, 0, frame
)

# The original values and range of each panel series, read straight from the emitters, as
# the filter tests them.
function _panel_oracle(emitters, cfg)
    e1 = first(emitters)
    get(f) = Float64[f(e) for e in emitters]
    both = (:lo, :hi)
    o = Dict{Symbol, Vector}()
    o[:photons] = [(get(e -> e.photons), cfg.photons, both)]
    o[:precision] = [(get(e -> e.σ_x), cfg.precision, (:hi,)), (get(e -> e.σ_y), cfg.precision, (:hi,))]
    hasproperty(e1, :pvalue) && (o[:pvalue] = [(get(e -> e.pvalue), cfg.pvalue, both)])
    psf = SMLMAnalysis._resolve_psf_bounds(emitters, cfg)
    if hasproperty(e1, :σx)
        bx = psf === nothing ? nothing : (psf.lo_x, psf.hi_x)
        by = psf === nothing ? nothing : (psf.lo_y, psf.hi_y)
        o[:psf] = [(get(e -> e.σx), bx, both), (get(e -> e.σy), by, both)]
    elseif hasproperty(e1, :σ)
        o[:psf] = [(get(e -> e.σ), psf === nothing ? nothing : (psf.lo, psf.hi), both)]
    end
    if hasproperty(e1, :z)
        o[:z] = [(get(e -> e.z), cfg.z, both)]
        o[:sigma_z] = [(get(e -> e.σ_z), cfg.sigma_z, both)]
    end
    return o
end

_in_span(edges, i, s) = s !== nothing && edges[i] >= s[1] && edges[i + 1] <= s[2]
_drawn_in(s, edges, span) = sum((s.counts[i] for i in eachindex(s.counts) if _in_span(edges, i, span)); init = 0)

# What each fit_quality.png panel draws, checked against the filter's own test on the original
# values. Every bound that rejects a fit lies strictly inside the view, on a bar edge. Each
# shaded span holds exactly the fits that side of the range rejects (those beyond the view on
# that side included), every rejected fit with no value is counted, and every passing fit in
# the view is drawn once outside the spans (or counted as unplaced at coincident bounds).
function _check_quality_panels(emitters, cfg)
    panels = SMLMAnalysis._quality_panels(emitters, cfg)
    oracle = _panel_oracle(emitters, cfg)
    @test sort(collect(keys(panels))) == sort(collect(keys(oracle)))
    for (key, p) in panels
        vlo, vhi = p.view
        @test vlo < vhi
        @test issorted(p.edges) && allunique(p.edges) && first(p.edges) == vlo && last(p.edges) == vhi
        for (s, (raw, b, sides)) in zip(p.series, oracle[key])
            @test s.sides == sides
            ok = filter(k -> !isnan(raw[k]) && !isnan(s.values[k]), eachindex(raw))
            if b === nothing
                @test s.spans == (; lo = nothing, hi = nothing)
                @test sum(s.counts) == count(k -> vlo <= s.values[k] <= vhi, ok)
                continue
            end
            lo, hi = b
            lowc(k) = :lo in sides && !(lo <= raw[k])
            highc(k) = :hi in sides && !(raw[k] <= hi)
            # Rejected fits with no agreeing bar are counted by side, never as passing.
            @test _drawn_in(s, p.edges, s.spans.lo) + s.lost_lo == count(lowc, ok)
            @test _drawn_in(s, p.edges, s.spans.hi) + s.lost_hi == count(highc, ok)
            @test s.invalid == count(k -> isnan(raw[k]) || isnan(s.values[k]), eachindex(raw))
            passing_in_view = count(k -> !lowc(k) && !highc(k) && vlo <= s.values[k] <= vhi, ok)
            @test sum(s.counts) + s.unplaced + s.lost == passing_in_view + count(k -> lowc(k) || highc(k), ok)
            @test s.unplaced <= passing_in_view
            # Placement, checked independently: when every fit's own bin (or end bar, for a
            # rejected fit beyond the view) agrees with its class, the bars are that histogram.
            nb = length(p.edges) - 1
            natural = zeros(Int, nb)
            agree = true
            for k in ok
                t = s.values[k]
                rej = lowc(k) || highc(k)
                i = vlo <= t <= vhi ? min(searchsortedlast(p.edges, t), nb) : (rej ? (t < vlo ? 1 : nb) : 0)
                i == 0 && continue
                agree &= _in_span(p.edges, i, s.spans.lo) == lowc(k) && _in_span(p.edges, i, s.spans.hi) == highc(k)
                natural[i] += 1
            end
            agree && @test s.counts == natural
            for (x, rejects) in ((s.pbounds[1], any(k -> !(lo <= raw[k]), ok)), (s.pbounds[2], any(k -> !(raw[k] <= hi), ok)))
                if rejects && SMLMAnalysis._panel_drawable(x)
                    @test vlo < x < vhi
                    @test x in p.edges
                    @test x in s.lines || x in s.dashed
                end
            end
        end
    end
    return panels
end

# What the figure itself draws on each panel's Axis, read back from Makie: the shaded spans,
# the bar rectangles (corners on stored edges, heights equal to the counts) and the threshold
# lines all match the panel data, and the x limits are the view.
function _check_quality_figure(emitters, cfg)
    panels = SMLMAnalysis._quality_panels(emitters, cfg)
    r = SMLMAnalysis._filter_quality_figure(emitters, cfg)
    Mk = SMLMAnalysis.CairoMakie.Makie
    for (key, p) in panels
        ax = r.axes[key]
        # The title carries the panel's note, and the note names every fit not drawn where it sits.
        note = SMLMAnalysis._panel_note(p)
        @test endswith(ax.title[], note)
        tot(f) = sum(f(s) for s in p.series)
        # One series counts fits; overlaid series count each series' values by name.
        for (f, phrase) in (
                (s -> s.beyond, "rejected beyond the axis"), (s -> s.invalid, "with NaN rejected"),
                (s -> s.lost, "rejected with no place"), (s -> s.unplaced, "passing at the bound"),
            )
            tot(f) > 0 || continue
            if length(p.series) == 1
                @test occursin("$(tot(f)) $phrase", note)
            else
                @test all(s -> f(s) == 0 || occursin("$(f(s)) $(s.name)", note), p.series)
                @test occursin("values $phrase", note)
            end
        end
        @test all(n -> occursin(n, note), reduce(vcat, [s.notes for s in p.series]; init = String[]))
        plots = ax.scene.plots
        kind(pl) = nameof(Mk.plotfunc(pl))
        spans = sort([(Float64(pl[1][]), Float64(pl[2][])) for pl in plots if kind(pl) === :vspan])
        shared = length(p.series) == 1 || allequal([(s.spans, s.lines, s.dashed) for s in p.series])
        want = [sp for (k, s) in enumerate(p.series) if !(shared && k > 1) for sp in (s.spans.lo, s.spans.hi) if sp !== nothing]
        @test spans == sort(want)
        bars = [pl[1][] for pl in plots if kind(pl) === :poly]
        drawn = [s for s in p.series if sum(s.counts) > 0]
        @test length(bars) == length(drawn)
        widths = diff(p.edges)
        for (polys, s) in zip(bars, drawn)
            y = p.density ? s.counts ./ (sum(s.counts) .* widths) : Float64.(s.counts)
            # Each bar is exactly the rectangle on its two stored edges, at its count's height.
            want_polys = [
                [(p.edges[i], 0.0), (p.edges[i + 1], 0.0), (p.edges[i + 1], y[i]), (p.edges[i], y[i])]
                    for i in eachindex(y) if y[i] > 0
            ]
            @test [[(Float64(q[1]), Float64(q[2])) for q in poly] for poly in polys] == want_polys
        end
        vl = reduce(vcat, [collect(pl[1][]) for pl in plots if kind(pl) === :vlines]; init = Float64[])
        @test all(x -> x in vl, reduce(vcat, [vcat(s.lines, s.dashed) for s in p.series]; init = Float64[]))
        @test ax.limits[][1] == p.view
    end
    return r
end

_fit2d(; photons = 1000.0, bg = 5.0, σ = 0.13, px = 0.01, py = 0.01, pvalue = 0.5, i = 1) =
    GaussMLE.Emitter2DFitSigma{Float64}(0.8, 0.8, photons, bg, σ, px, py, 0.0, 30.0, 0.5, 0.002, pvalue, i, 1, 0, i)

# Mostly ordinary fits; a bright tail above the 98th percentile, a few wide PSFs, poor
# precisions, p-values far below every view, p = 0 and NaN p-values (DetectFit keeps them).
_sweep_fits(rng; n = 2000) = [
    _fit2d(
        photons = i <= 60 ? 6000 + 3000 * rand(rng) : 1000 + 3000 * rand(rng),
        σ = i % 50 == 0 ? 0.25 : 0.12 + 0.03 * rand(rng),
        px = i % 97 == 0 ? 0.5 : 0.005 + 0.02 * rand(rng),
        py = 0.005 + 0.02 * rand(rng),
        pvalue = i % 40 == 0 ? 1.0e-12 : (i % 333 == 0 ? 0.0 : (i % 401 == 0 ? NaN : rand(rng))),
        i = i
    ) for i in 1:n
]

# Bars inside a series' shaded spans.
_shaded(s, edges) = sum(
    (s.counts[i] for i in eachindex(s.counts) if _in_span(edges, i, s.spans.lo) || _in_span(edges, i, s.spans.hi));
    init = 0
)
