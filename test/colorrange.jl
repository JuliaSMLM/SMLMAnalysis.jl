using SMLMAnalysis
using GaussMLE
using Random
using Test

# Makie's interpolated colormaps throw on a colorrange with equal ends (CairoMakie 0.15.6 /
# Makie 0.24.6, the lowest-compat CI job). Each diagnostic plot must go through
# _plot_colorrange, so a degenerate range from data cannot abort a save.

@testset "Plot color range helper" begin
    f = SMLMAnalysis._plot_colorrange

    # A nondegenerate range is returned unchanged.
    @test f(0.0, 1.0) == (0.0, 1.0)
    @test f(-2.5, 7.0) == (-2.5, 7.0)

    # Equal ends widen to a finite range around the shared value.
    for c in (20.0, 0.0, -3.0, 1.0e6)
        lo, hi = f(c, c)
        @test isfinite(lo) && isfinite(hi)
        @test lo < c < hi
    end

    # Ends distinct in Float64 but equal in Float32 (Makie's Vec2f) also widen.
    lo, hi = f(1000.0, 1000.00001)
    @test Float32(hi) > Float32(lo)

    # NaN or Inf ends: the finite end centers the range; no finite end gives (0, 1).
    @test f(NaN, NaN) == (0.0, 1.0)
    @test f(Inf, Inf) == (0.0, 1.0)
    @test f(-Inf, NaN) == (0.0, 1.0)
    lo, hi = f(NaN, 5.0)
    @test lo < 5.0 < hi
    lo, hi = f(2.0, Inf)
    @test lo < 2.0 < hi

    # Ends equal and beyond Float32 (Inf32 in Makie's Vec2f) give a bounded finite range; a
    # value at the Float32 limit widens without rounding an end to Inf32.
    for g in (f, SMLMAnalysis._plot_limits), c in (1.0e40, -1.0e40, Float64(floatmax(Float32)), -Float64(floatmax(Float32)))
        lo, hi = g(c, c)
        @test isfinite(Float32(lo)) && isfinite(Float32(hi))
        @test Float32(hi) > Float32(lo)
    end
end

@testset "Plot axis limit helper" begin
    f = SMLMAnalysis._plot_limits
    @test f(0.0, 250.0) == (0.0, 250.0)
    # Equal ends (one localization, a constant column) widen to a readable window.
    lo, hi = f(130.0, 130.0)
    @test lo < 130.0 < hi && hi - lo > 1
    lo, hi = f(0, 0)
    @test lo < 0 < hi
    # An open filter bound has no finite end to show: the finite end centers the window.
    lo, hi = f(0, Inf)
    @test isfinite(lo) && isfinite(hi) && lo < 0 < hi
    # Reversed finite ends are put in order.
    @test f(0.0, -5.0) == (-5.0, 0.0)
end

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

@testset "Filter quality figure limits on degenerate data" begin
    cam = IdealCamera(16, 16, 0.1)
    smld_of(es) = SMLMAnalysis.BasicSMLD(es, cam, 20, 1, Dict{String, Any}())
    varied = [
        _sigma_emitter(σ = 0.12 + 0.002i, bg = 5.0 + i, prec = 0.005 + 0.0005i, pvalue = i / 21, frame = i, photons = 500.0 + 50i)
            for i in 1:20
    ]

    # One localization: PSF width limits (130, 130), and every other panel from one value.
    one = smld_of([_sigma_emitter()])
    for cfg in (FilterConfig(), FilterConfig(psf_sigma = :auto), FilterConfig(psf_sigma = (0.1, 0.2)))
        @test _quality_figure_saves(one, cfg)
    end
    xy = GaussMLE.Emitter2DFitSigmaXY{Float64}(
        0.8, 0.8, 1000.0, 5.0, 0.13, 0.14, 0.01, 0.01, 0.0, 30.0, 0.5, 0.002, 0.002, 0.4, 1, 1, 0, 1
    )
    for cfg in (FilterConfig(), FilterConfig(psf_sigma = :auto), FilterConfig(psf_sigma = (-Inf, Inf)))
        @test _quality_figure_saves(smld_of([xy]), cfg)
    end
    @test _quality_figure_saves(one, FilterConfig(psf_sigma = (-Inf, Inf)))
    # Two distinct PSF widths both lie outside their 2-98% display window (empty histogram).
    @test _quality_figure_saves(smld_of([_sigma_emitter(σ = 0.12), _sigma_emitter(σ = 0.14, frame = 2)]), FilterConfig())

    # Zero background gives background limits (0, 0).
    @test _quality_figure_saves(smld_of([_sigma_emitter(bg = 0.0, frame = i) for i in 1:5]), FilterConfig())

    # Open filter bounds (an ordinary user setting) give infinite limits.
    for cfg in (
            FilterConfig(precision = (0, Inf)), FilterConfig(photons = (-Inf, Inf)),
            FilterConfig(psf_sigma = (0, Inf)), FilterConfig(pvalue = (0, 1)),
        )
        @test _quality_figure_saves(smld_of(varied), cfg)
        plain, _ = analyze(smld_of(varied), cfg; verbose = Verbosity.SILENT)
        mktempdir() do outdir
            drawn, _ = analyze(smld_of(varied), cfg; outdir = outdir, checkpoint = Checkpoint.NONE)
            @test drawn.emitters == plain.emitters
            @test isfile(joinpath(outdir, "00_filter", "fit_quality.png"))
        end
    end
    e3 = [
        SMLMAnalysis.Emitter3DFit{Float64}(
            0.8, 0.8, 0.1, 1000.0, 10.0, 0.01, 0.01, 0.03, 30.0, 1.0; frame = 1
        ),
    ]
    @test _quality_figure_saves(smld_of(e3), FilterConfig(z = (-Inf, Inf), sigma_z = (0, Inf)))

    # p-values equal to 1: the 2nd percentile of log10(p) is 0.
    @test _quality_figure_saves(smld_of([_sigma_emitter(pvalue = 1.0, frame = i) for i in 1:5]), FilterConfig())
    @test _quality_figure_saves(smld_of([_sigma_emitter(pvalue = 1.0)]), FilterConfig(pvalue = (1.0e-3, 1.0)))

    # Every precision above the 1 μm display cutoff leaves nothing to take a quantile of.
    far = smld_of([_sigma_emitter(prec = 2.0)])
    for cfg in (FilterConfig(), FilterConfig(precision = (0, 0.02)))
        @test _quality_figure_saves(far, cfg)
    end
end

@testset "Box overlay with non-finite pixels" begin
    cam = IdealCamera(16, 16, 0.1)
    batch = SMLMAnalysis.ROIBatch(fill(20.0f0, 7, 7, 1), [5], [5], [1], cam)
    some_nan = fill(20.0f0, 16, 16, 1)
    some_nan[3, 4, 1] = NaN32
    some_nan[8:10, 8:10, 1] .= 200.0f0
    for images in (some_nan, fill(NaN32, 16, 16, 1))
        mktempdir() do dir
            SMLMAnalysis._save_box_overlay(dir, "overlay.png", images, batch, [:green])
            @test isfile(joinpath(dir, "overlay.png"))
        end
    end
end

@testset "Intensity filter diagnostics with constant photons" begin
    # Constant photons make the photon scatter range (p01, p99), the binned rate
    # heatmap range (min, max over bins) and the uniform field heatmap range all degenerate.
    cam = IdealCamera(16, 16, 0.1)
    emitters = [
        SMLMAnalysis.Emitter2DFit{Float64}(
            0.1 * ((i - 1) % 14), 0.1 * ((i - 1) ÷ 14), 1000.0, 10.0,
            0.01, 0.01, 30.0, 1.0; frame = i
        ) for i in 1:200
    ]
    smld = SMLMAnalysis.BasicSMLD(emitters, cam, 200, 1, Dict{String, Any}())
    cfg = IntensityFilterConfig(field_mode = :uniform, estimate_p2 = false, n_bins = 2)

    # Without an outdir no figures are drawn; the analysis result is the reference.
    _, info_plain = SMLMAnalysis.intensityfilter_step(smld, cfg)

    mktempdir() do outdir
        filtered, info = SMLMAnalysis.intensityfilter_step(smld, cfg; outdir = outdir)

        # The figures must not change the analysis: same emitters kept with and without them.
        @test info.n_after == info_plain.n_after
        @test length(filtered.emitters) == info_plain.n_after
        @test info.n_before == 200

        field_png = String[]
        for (root, _, files) in walkdir(outdir), file in files
            file == "excitation_field.png" && push!(field_png, joinpath(root, file))
        end
        @test length(field_png) == 1
        @test all(p -> filesize(p) > 0, field_png)
    end
end

@testset "Intensity filter diagnostics with zero photons" begin
    # Zero photons put the raw and field-normalized photon axes at (0, 0).
    cam = IdealCamera(16, 16, 0.1)
    emitters = [
        SMLMAnalysis.Emitter2DFit{Float64}(
            0.1 * ((i - 1) % 14), 0.1 * ((i - 1) ÷ 14), 0.0, 10.0,
            0.01, 0.01, 30.0, 1.0; frame = i
        ) for i in 1:200
    ]
    smld = SMLMAnalysis.BasicSMLD(emitters, cam, 200, 1, Dict{String, Any}())
    cfg = IntensityFilterConfig(field_mode = :uniform, estimate_p2 = false, n_bins = 2)
    _, info_plain = SMLMAnalysis.intensityfilter_step(smld, cfg)
    mktempdir() do outdir
        filtered, info = SMLMAnalysis.intensityfilter_step(smld, cfg; outdir = outdir)
        @test info.n_after == info_plain.n_after
        @test length(filtered.emitters) == info_plain.n_after
        pngs = [f for (_, _, files) in walkdir(outdir) for f in files if f == "photon_distribution.png"]
        @test length(pngs) == 1
    end
end

@testset "Panel bar counts and pass percentage" begin
    pc = SMLMAnalysis._panel_counts
    nospan = (; lo = nothing, hi = nothing)
    both = (:lo, :hi)
    edges = collect(range(-2.0, 0.0, length = 51))
    # A value on an edge goes in the bar that starts there: -1.8 and -1.79 both belong in
    # [-1.8, -1.76). Dividing by the bin width put -1.8 in the bar below.
    vals = vcat(fill(log10(0.01), 10), fill(log10(10^-1.8), 50), fill(log10(10^-1.79), 50))
    c = pc(vals, vals, nothing, both, edges, nospan)
    @test c.counts[1] == 10 && c.counts[6] == 100 && sum(c.counts) == 110 && c.beyond == 0
    # The top edge goes in the last bar; with no filter, values beyond the view are not drawn.
    v = [0.0, 0.5, -3.0, NaN]
    c = pc(v, v, nothing, both, edges, nospan)
    @test c.counts == [zeros(Int, 49); 1] && c.beyond == 0 && c.invalid == 0
    e = [0.0, 1.0, 2.0, 3.0]
    # A value on an upper bound passes, so it goes below the bound, outside the shaded span.
    c = pc([2.0], [2.0], (-Inf, 2.0), both, e, (; lo = nothing, hi = (2.0, 3.0)))
    @test c.counts == [0, 1, 0]
    # The class comes from the original value: a p-value just above 0.01 is rejected by
    # (0, 0.01) although its log10 equals the bound's.
    p = nextfloat(0.01)
    @test log10(p) == log10(0.01)
    c = pc([log10(p)], [p], (0.0, 0.01), both, [-3.0, -2.0, -1.0], (; lo = nothing, hi = (-2.0, -1.0)))
    @test c.counts == [0, 1]
    # Coincident bounds: a passing value has no bar outside both spans.
    c = pc([1.0], [1.0], (1.0, 1.0), both, e, (; lo = (0.0, 1.0), hi = (1.0, 3.0)))
    @test c.counts == [0, 0, 0] && c.unplaced == 1
    # Beyond the view: rejected values land in the end bar; a NaN is counted, not drawn.
    c = pc([-5.0, -Inf, 9.0, NaN], [-5.0, -Inf, 9.0, NaN], (1.0, 2.0), both, e, (; lo = (0.0, 1.0), hi = (2.0, 3.0)))
    @test c.counts == [2, 0, 1] && c.beyond == 3 && c.invalid == 1
    pp = SMLMAnalysis._pass_percent
    @test pp([0.4, 0.4, 0.4], (0.4, 1.0)) == 100.0
    @test pp([0.1, 0.5, 0.9], (0.2, 0.8)) == round(100 / 3, digits = 1)
end

@testset "Panel geometry" begin
    g(window, bounds; kw...) = SMLMAnalysis._panel_geometry(
        window, [(; lo, hi, rejects_lo = rl, rejects_hi = rh) for (lo, hi, rl, rh) in bounds]; kw...
    )
    # A rejecting bound beyond the window extends the view, with margin, and is an edge.
    r = g((0.0, 4463.0), [(-Inf, 6000.0, false, true)])
    @test r.view[1] == 0.0 && 6000 < r.view[2] && r.lines == [6000.0] && 6000.0 in r.edges
    # Sign-safe margin: positive bounds below and above the data stay inside the view.
    r = g((200.0, 300.0), [(150.0, 350.0, false, false)])
    @test r.view[1] < 150 && r.view[2] > 350 && r.lines == [150.0, 350.0]
    # A bound that rejects nothing and lies more than a window width away does not widen it.
    r = g((0.0, 4000.0), [(200.0, 1.0e20, false, false)])
    @test r.lines == [200.0] && r.view[2] == 4000.0
    # A rejecting bound as large as Float32 allows is still drawn inside the view.
    r = g((0.0, 1000.0), [(1.0e38, Inf, true, false)])
    @test r.view[1] < 1.0e38 < r.view[2] && 1.0e38 in r.edges
    # A non-rejecting bound at the panel's floor is not drawn (photons or precision of 0).
    r = g((0.0, 10.0), [(0.0, Inf, false, false)]; floor = 0.0)
    @test r.view == (0.0, 10.0) && isempty(r.lines)
    # Every bound inside the final view is an edge, even one that did not widen it.
    r = g((0.0, 1000.0), [(1.0e6, 5.0e5, true, false)])
    @test 5.0e5 in r.edges && 1.0e6 in r.edges
    @test issorted(r.edges) && allunique(r.edges)
    s = SMLMAnalysis._panel_spans
    @test s((0.0, 10.0), (Inf, Inf), (:lo, :hi)) == (; lo = (0.0, 10.0), hi = nothing)
    @test s((0.0, 10.0), (NaN, 1.0), (:lo, :hi)) == (; lo = (0.0, 10.0), hi = (1.0, 10.0))
    @test s((0.0, 10.0), (2.0, 8.0), (:hi,)) == (; lo = nothing, hi = (8.0, 10.0))
    @test s((0.0, 10.0), (-Inf, Inf), (:lo, :hi)) == (; lo = nothing, hi = nothing)
end

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
            @test _drawn_in(s, p.edges, s.spans.lo) == count(lowc, ok)
            @test _drawn_in(s, p.edges, s.spans.hi) == count(highc, ok)
            @test s.invalid == count(k -> isnan(raw[k]) || isnan(s.values[k]), eachindex(raw))
            passing_in_view = count(k -> !lowc(k) && !highc(k) && vlo <= s.values[k] <= vhi, ok)
            @test sum(s.counts) + s.unplaced == passing_in_view + count(k -> lowc(k) || highc(k), ok)
            for (x, rejects) in ((s.pbounds[1], any(k -> raw[k] < lo, ok)), (s.pbounds[2], any(k -> raw[k] > hi, ok)))
                if isfinite(x) && rejects && abs(x) <= SMLMAnalysis._PANEL_MAX_BOUND
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
            got = zeros(length(y))
            for poly in polys
                xs = unique(first.(poly))
                @test length(xs) == 2 && all(x -> x in p.edges, xs)
                i = findfirst(==(minimum(xs)), p.edges)
                @test p.edges[i + 1] == maximum(xs)
                got[i] = maximum(last.(poly))
            end
            @test got == y
        end
        vl = reduce(vcat, [collect(pl[1][]) for pl in plots if kind(pl) === :vlines]; init = Float64[])
        @test all(x -> x in vl, reduce(vcat, [vcat(s.lines, s.dashed) for s in p.series]; init = Float64[]))
        @test ax.limits[][1] == p.view
    end
    return r
end

_fit2d(; photons = 1000.0, bg = 5.0, σ = 0.13, px = 0.01, py = 0.01, pvalue = 0.5, i = 1) =
    GaussMLE.Emitter2DFitSigma{Float64}(0.8, 0.8, photons, bg, σ, px, py, 0.0, 30.0, 0.5, 0.002, pvalue, i, 1, 0, i)

@testset "Filter figures draw what the filter rejects" begin
    cam = IdealCamera(16, 16, 0.1)
    smld_of(es) = SMLMAnalysis.BasicSMLD(es, cam, length(es), 1, Dict{String, Any}())
    rng = Random.Xoshiro(7)
    n = 2000
    # Mostly ordinary fits; a bright tail above the 98th percentile, a few wide PSFs, poor
    # precisions, p-values far below every view, p = 0 and NaN p-values (DetectFit keeps them).
    es = [
        _fit2d(
            photons = i <= 60 ? 6000 + 3000 * rand(rng) : 1000 + 3000 * rand(rng),
            σ = i % 50 == 0 ? 0.25 : 0.12 + 0.03 * rand(rng),
            px = i % 97 == 0 ? 0.5 : 0.005 + 0.02 * rand(rng),
            py = 0.005 + 0.02 * rand(rng),
            pvalue = i % 40 == 0 ? 1.0e-12 : (i % 333 == 0 ? 0.0 : (i % 401 == 0 ? NaN : rand(rng))),
            i = i
        ) for i in 1:n
    ]
    smld = smld_of(es)
    removed(cfg) = n - length(SMLMAnalysis._filter_smld(smld, cfg).emitters)
    shaded(s, edges) = sum(
        (s.counts[i] for i in eachindex(s.counts) if _in_span(edges, i, s.spans.lo) || _in_span(edges, i, s.spans.hi));
        init = 0
    )

    # Review case: an upper photon bound above the 98th percentile. Every fit it rejects is
    # drawn in the shaded band, and the passing fits between the percentile and the bound show.
    cfg = FilterConfig(photons = (-Inf, 6000))
    p = _check_quality_panels(es, cfg)[:photons]
    @test shaded(p.series[1], p.edges) == removed(cfg) == 60
    @test sum(p.series[1].counts) == n
    _check_quality_figure(es, cfg)
    # Review case: a lower p-value bound far below the data's 2% quantile. p = 0 is drawn in
    # the band; NaN p-values are counted and named in the panel title.
    cfg = FilterConfig(pvalue = (1.0e-9, 1))
    p = _check_quality_panels(es, cfg)[:pvalue]
    @test shaded(p.series[1], p.edges) + p.series[1].invalid == removed(cfg) == 50 + 6 + 4
    @test occursin("4 with NaN rejected", SMLMAnalysis._panel_note(p))
    _check_quality_figure(es, cfg)
    for cfg in (
            FilterConfig(photons = (2000, Inf)), FilterConfig(photons = (1500, 3500)),
            FilterConfig(photons = (1.0e6, Inf)), FilterConfig(photons = (-Inf, -1)),
            FilterConfig(photons = (200, 1.0e20)), FilterConfig(photons = (0, Inf)),
            FilterConfig(photons = (3000, 2000)), FilterConfig(photons = (1.0e38, Inf)),
            FilterConfig(photons = (Inf, Inf)), FilterConfig(photons = (NaN, Inf)),
            FilterConfig(pvalue = (0.2, 0.8)), FilterConfig(pvalue = (0, 0)), FilterConfig(pvalue = (-Inf, -0.1)),
            FilterConfig(pvalue = (1.0e-3, 1)), FilterConfig(pvalue = (1.0, 1.0)), FilterConfig(pvalue = (NaN, 1)),
            FilterConfig(psf_sigma = :auto), FilterConfig(psf_sigma = (0.13, 0.14)),
            FilterConfig(psf_sigma = (0.3, Inf)),
            FilterConfig(precision = (0.01, 0.02)), FilterConfig(precision = (0.02, Inf)), FilterConfig(precision = (0, 0.015)),
        )
        p = _check_quality_panels(es, cfg)
        key = cfg.photons !== nothing ? :photons : cfg.pvalue !== nothing ? :pvalue : cfg.psf_sigma !== nothing ? :psf : :precision
        s = p[key].series[1]
        # One series and a range whose two sides do not overlap: the shaded bars plus the
        # rejected fits with no value are exactly what the filter removes.
        if length(p[key].series) == 1 && s.pbounds !== nothing && !(s.pbounds[1] > s.pbounds[2])
            @test shaded(s, p[key].edges) + s.invalid == removed(cfg)
        end
        _check_quality_figure(es, cfg)
        @test _quality_figure_saves(smld, cfg)
    end
    # The precision lower bound (on max(σx, σy)) is a dashed line, the upper bound shaded.
    p = _check_quality_panels(es, FilterConfig(precision = (0.01, 0.02)))[:precision]
    @test all(s -> s.dashed == [10.0] && s.lines == [20.0], p.series)

    # Coincident bounds: the passing fit is named in the title, not drawn in a shaded bar.
    one = [_fit2d(photons = 1000.0)]
    p = _check_quality_panels(one, FilterConfig(photons = (1000, 1000)))[:photons]
    @test p.series[1].unplaced == 1 && shaded(p.series[1], p.edges) == 0
    @test occursin("1 passing at the bound", SMLMAnalysis._panel_note(p))
    # A fit on the upper bound passes: no bar crosses the threshold, and it is not shaded.
    p = _check_quality_panels([_fit2d(photons = 1000.01)], FilterConfig(photons = (-Inf, 1000.01)))[:photons]
    @test shaded(p.series[1], p.edges) == 0
    _check_quality_figure([_fit2d(photons = 1000.01)], FilterConfig(photons = (-Inf, 1000.01)))
    # log10 rounding: p just above 0.01 is rejected by (0, 0.01) and drawn shaded.
    p = _check_quality_panels([_fit2d(pvalue = nextfloat(0.01)), _fit2d(pvalue = 0.005, i = 2)], FilterConfig(pvalue = (0, 0.01)))[:pvalue]
    @test shaded(p.series[1], p.edges) == 1
    # A photon count too large for an Int does not stop the figure.
    @test _quality_figure_saves(smld_of([_fit2d(photons = 1.0e20)]), FilterConfig())

    # Anisotropic PSF: per-axis :auto bounds differ, explicit bounds are shared.
    xy = [
        GaussMLE.Emitter2DFitSigmaXY{Float64}(
            0.8, 0.8, 1000.0, 5.0, 0.12 + 0.02 * rand(rng) + (i % 30 == 0 ? 0.1 : 0.0),
            0.15 + 0.02 * rand(rng), 0.01, 0.01, 0.0, 30.0, 0.5, 0.002, 0.002, rand(rng), i, 1, 0, i
        ) for i in 1:500
    ]
    for cfg in (FilterConfig(psf_sigma = :auto), FilterConfig(psf_sigma = (0.13, 0.16)), FilterConfig(psf_sigma = (0.2, 0.4)))
        p = _check_quality_panels(xy, cfg)[:psf]
        @test all(s -> s.pbounds !== nothing, p.series)
        _check_quality_figure(xy, cfg)
        @test _quality_figure_saves(smld_of(xy), cfg)
    end

    # 3D, all z in 0.2-0.3 µm. Positive bounds below the data stay inside the view (the old
    # 1.1x widening moved them toward zero). z = (0.1, 0.4) rejects nothing and lies more than
    # a window width away, so it widens nothing and nothing is shaded.
    e3 = [
        SMLMAnalysis.Emitter3DFit{Float64}(
            0.8, 0.8, 0.2 + 0.1 * rand(rng), 1000.0, 10.0, 0.01, 0.01, 0.02 + 0.03 * rand(rng), 30.0, 1.0; frame = i
        ) for i in 1:500
    ]
    p = _check_quality_panels(e3, FilterConfig(z = (0.15, 0.35)))[:z]
    @test p.view[1] < 150 && p.view[2] > 350 && p.series[1].lines ≈ [150.0, 350.0]
    p = _check_quality_panels(e3, FilterConfig(z = (0.1, 0.4)))[:z]
    @test isempty(p.series[1].lines) && p.series[1].spans == (; lo = nothing, hi = nothing)
    for cfg in (FilterConfig(z = (0.22, 0.28)), FilterConfig(z = (-0.4, -0.1)), FilterConfig(sigma_z = (0.03, 0.04)), FilterConfig(sigma_z = (0, 0.03)))
        _check_quality_panels(e3, cfg)
        _check_quality_figure(e3, cfg)
        @test _quality_figure_saves(smld_of(e3), cfg)
    end
end

@testset "Fit overlay title" begin
    t = SMLMAnalysis._fit_overlay_title
    es = [_sigma_emitter()]
    # Each range is printed exactly as the filter tests it.
    @test occursin("red=photons not in [500, 5000]", t(FilterConfig(photons = (500, 5000)), es))
    @test occursin("red=photons not in [500, Inf]", t(FilterConfig(photons = (500, Inf)), es))
    @test occursin("red=photons not in [-Inf, 5000]", t(FilterConfig(photons = (-Inf, 5000)), es))
    @test occursin("red=photons not in [200, 1.23456e20]", t(FilterConfig(photons = (200, 1.23456e20)), es))
    @test occursin("red=photons not in [500.5, 5000]", t(FilterConfig(photons = (500.5, 5000)), es))
    # Ranges that reject everything are labeled too; a range that rejects nothing is not.
    @test occursin("red=photons not in [Inf, Inf]", t(FilterConfig(photons = (Inf, Inf)), es))
    @test occursin("red=photons not in [NaN, Inf]", t(FilterConfig(photons = (NaN, Inf)), es))
    @test !occursin("photons", t(FilterConfig(photons = (-Inf, Inf)), es))
    @test !occursin("prec", t(FilterConfig(precision = (0, Inf)), es))
    @test occursin("orange=prec not in [0, 20]nm", t(FilterConfig(precision = (0, 0.02)), es))
    @test occursin("orange=prec not in [20, Inf]nm", t(FilterConfig(precision = (0.02, Inf)), es))
    @test occursin("orange=prec not in [5, 20]nm", t(FilterConfig(precision = (0.005, 0.02)), es))
    # The p-value label names the range, which a NaN p-value also fails.
    @test occursin("purple=pval not in [0.2, 0.8]", t(FilterConfig(pvalue = (0.2, 0.8)), es))
    nan_p = _sigma_emitter(pvalue = NaN)
    cfg = FilterConfig(pvalue = (0.001, 1.0))
    @test !SMLMAnalysis._in_range(nan_p.pvalue, cfg.pvalue)
    @test occursin("purple=pval not in [0.001, 1.0]", t(cfg, [nan_p]))
    # PSF bounds as the filter applies them, without rounding two distinct bounds together;
    # a degenerate :auto that the filter skips has no label.
    @test occursin("cyan=psf_σ not in [100, 200]nm", t(FilterConfig(psf_sigma = (0.1, 0.2)), es))
    @test occursin("cyan=psf_σ not in [130.01, 130.02]nm", t(FilterConfig(psf_sigma = (0.13001, 0.13002)), es))
    @test occursin("cyan=psf_σ not in [117, 143]nm", t(FilterConfig(psf_sigma = :auto), es))
    @test !occursin("cyan", t(FilterConfig(psf_sigma = :auto), [_sigma_emitter(σ = 0.0)]))
    xy = GaussMLE.Emitter2DFitSigmaXY{Float64}(0.8, 0.8, 1000.0, 5.0, 0.13, 0.14, 0.01, 0.01, 0.0, 30.0, 0.5, 0.002, 0.002, 0.4, 1, 1, 0, 1)
    @test occursin("cyan=psf_σx not in [100, 200]nm or σy not in [100, 200]nm", t(FilterConfig(psf_sigma = (0.1, 0.2)), [xy]))
end

@testset "Filter figures save for each side of a range" begin
    cam = IdealCamera(16, 16, 0.1)
    smld_of(es) = SMLMAnalysis.BasicSMLD(es, cam, 20, 1, Dict{String, Any}())
    @test _quality_figure_saves(smld_of([_sigma_emitter(photons = 1000.0)]), FilterConfig(photons = (-Inf, 800)))
    @test _quality_figure_saves(smld_of([_sigma_emitter(prec = 0.01)]), FilterConfig(precision = (0.02, Inf)))
    @test _quality_figure_saves(smld_of([_sigma_emitter(pvalue = 0.9)]), FilterConfig(pvalue = (0.2, 0.8)))
    e3 = [
        SMLMAnalysis.Emitter3DFit{Float64}(
            0.8, 0.8, 0.1, 1000.0, 10.0, 0.01, 0.01, 0.01, 30.0, 1.0; frame = 1
        ),
    ]
    @test _quality_figure_saves(smld_of(e3), FilterConfig(sigma_z = (0.02, Inf)))
    # Unequal axes: the filter tests max(σx, σy) = 30 nm, so a 20 nm lower bound keeps this fit.
    unequal = GaussMLE.Emitter2DFitSigma{Float64}(0.8, 0.8, 1000.0, 5.0, 0.13, 0.01, 0.03, 0.0, 30.0, 0.5, 0.002, 0.4, 1, 1, 0, 1)
    @test _quality_figure_saves(smld_of([unequal]), FilterConfig(precision = (0.02, Inf)))
    es = [_sigma_emitter(pvalue = p, frame = i) for (i, p) in enumerate((0.1, 0.5, 0.9))]
    for cfg in (FilterConfig(pvalue = (-Inf, -0.1)), FilterConfig(pvalue = (0, 0)), FilterConfig(pvalue = (0.2, 0.8)))
        @test _quality_figure_saves(smld_of(es), cfg)
    end
end
