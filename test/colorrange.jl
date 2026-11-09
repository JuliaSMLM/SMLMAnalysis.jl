using SMLMAnalysis
using GaussMLE
using Random
using Test

# The shared checks of the filter figures (also used by long/filter_figures.jl).
include("filter_checks.jl")

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


@testset "Filter figures draw what the filter rejects" begin
    cam = IdealCamera(16, 16, 0.1)
    smld_of(es) = SMLMAnalysis.BasicSMLD(es, cam, length(es), 1, Dict{String, Any}())
    rng = Random.Xoshiro(7)
    es = _sweep_fits(rng)
    n = length(es)
    smld = smld_of(es)
    removed(cfg) = n - length(SMLMAnalysis._filter_smld(smld, cfg).emitters)
    shaded = _shaded

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
    # Rejected fits with no agreeing bar are counted as rejected, never as passing: p = 0 is
    # rejected by the lower bound of (0.1, 0) but sits at the upper bound's -Inf in log10.
    for (ps, b) in (([0.0, 0.5], (0.1, 0.0)), ([0.0], (NaN, 0.0)))
        es2 = [_fit2d(pvalue = q, i = i) for (i, q) in enumerate(ps)]
        p = _check_quality_panels(es2, FilterConfig(pvalue = b))[:pvalue]
        @test p.series[1].unplaced == 0 && p.series[1].lost >= 1
        @test occursin("rejected with no place on the axis", SMLMAnalysis._panel_note(p))
        _check_quality_figure(es2, FilterConfig(pvalue = b))
    end
    # Reversed bounds that log10 collapses: both sides reject the fit.
    lb = 1.0e-12
    _check_quality_panels([_fit2d(pvalue = nextfloat(lb))], FilterConfig(pvalue = (nextfloat(nextfloat(lb)), lb)))
    # A rejecting bound as large as Float32 allows is drawn inside the view.
    p = _check_quality_panels([_fit2d(photons = 1000.0)], FilterConfig(photons = (3.2e38, Inf)))[:photons]
    @test p.view[1] < 3.2e38 < p.view[2] && 3.2e38 in p.series[1].lines
    # Precision lower bounds that cannot be drawn are named, and nothing claims to be shaded.
    for lo in (1.0e40, Inf, NaN)
        p = _check_quality_panels([_fit2d(px = 0.01, py = 0.01)], FilterConfig(precision = (lo, Inf)))[:precision]
        note = SMLMAnalysis._panel_note(p)
        @test occursin("lower limit $(lo) on max(σ_x, σ_y) not drawn", note) && !occursin("all shaded", note)
        @test all(s -> s.spans == (; lo = nothing, hi = nothing), p.series)
    end
    # One fit with both precision components beyond the axis is one rejected fit: the title
    # counts values per component, not two fits.
    pr = vcat([_fit2d(px = 0.01, py = 0.01, i = i) for i in 1:99], [_fit2d(px = 0.5, py = 0.5, i = 100)])
    p = _check_quality_panels(pr, FilterConfig(precision = (0, 0.02)))[:precision]
    note = SMLMAnalysis._panel_note(p)
    @test occursin("1 σ_x and 1 σ_y values rejected beyond the axis", note) && !occursin("2 rejected", note)
    @test 100 - length(SMLMAnalysis._filter_smld(smld_of(pr), FilterConfig(precision = (0, 0.02))).emitters) == 1
    _check_quality_figure(pr, FilterConfig(precision = (0, 0.02)))
    # A photon count too large for an Int does not stop the figure.
    @test _quality_figure_saves(smld_of([_fit2d(photons = 1.0e20)]), FilterConfig())

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
    @test occursin("orange=prec not in [0, 0.02]μm", t(FilterConfig(precision = (0, 0.02)), es))
    @test occursin("orange=prec not in [0.02, Inf]μm", t(FilterConfig(precision = (0.02, Inf)), es))
    @test occursin("orange=prec not in [0.005, 0.02]μm", t(FilterConfig(precision = (0.005, 0.02)), es))
    # A bound just above a fit's value prints so that the fit reads as rejected.
    @test occursin("red=photons not in [1000.000000001, Inf]", t(FilterConfig(photons = (1000.000000001, Inf)), es))
    # The p-value label names the range, which a NaN p-value also fails.
    @test occursin("purple=pval not in [0.2, 0.8]", t(FilterConfig(pvalue = (0.2, 0.8)), es))
    nan_p = _sigma_emitter(pvalue = NaN)
    cfg = FilterConfig(pvalue = (0.001, 1.0))
    @test !SMLMAnalysis._in_range(nan_p.pvalue, cfg.pvalue)
    @test occursin("purple=pval not in [0.001, 1.0]", t(cfg, [nan_p]))
    # PSF bounds as the filter applies them, without rounding two distinct bounds together;
    # a degenerate :auto that the filter skips has no label.
    @test occursin("cyan=psf_σ not in [0.1, 0.2]μm", t(FilterConfig(psf_sigma = (0.1, 0.2)), es))
    @test occursin("cyan=psf_σ not in [0.13001, 0.13002]μm", t(FilterConfig(psf_sigma = (0.13001, 0.13002)), es))
    # :auto bounds print as computed: a fit of 0.126 μm below a bound of 0.12600000000000003
    # reads as outside it, as its cyan box says.
    auto_es = [_sigma_emitter(σ = σ, frame = i) for (i, σ) in enumerate((0.14, 0.14, 0.126))]
    b = SMLMAnalysis._resolve_psf_bounds(auto_es, FilterConfig(psf_sigma = :auto))
    @test occursin("cyan=psf_σ not in [$(b.lo), $(b.hi)]μm", t(FilterConfig(psf_sigma = :auto), auto_es))
    @test SMLMAnalysis._fit_box_color(auto_es[3], FilterConfig(psf_sigma = :auto), b) === (0.126 < b.lo ? :cyan : :green)
    @test !occursin("cyan", t(FilterConfig(psf_sigma = :auto), [_sigma_emitter(σ = 0.0)]))
    xy = GaussMLE.Emitter2DFitSigmaXY{Float64}(0.8, 0.8, 1000.0, 5.0, 0.13, 0.14, 0.01, 0.01, 0.0, 30.0, 0.5, 0.002, 0.002, 0.4, 1, 1, 0, 1)
    @test occursin("cyan=psf_σx not in [0.1, 0.2]μm or σy not in [0.1, 0.2]μm", t(FilterConfig(psf_sigma = (0.1, 0.2)), [xy]))
end

@testset "Fit overlay legend fits the figure" begin
    # The tutorial's FilterConfig, and anisotropic :auto with a photon range, give legends far
    # wider than one line; wrapped to the figure, no line is wider than it and nothing is lost.
    iso = [_sigma_emitter(σ = 0.12 + 0.0005i, frame = i) for i in 1:60]
    xy = [
        GaussMLE.Emitter2DFitSigmaXY{Float64}(
            0.8, 0.8, 1000.0, 5.0, 0.12 + 0.0007i, 0.15 + 0.0003i, 0.01, 0.01, 0.0, 30.0, 0.5, 0.002, 0.002, 0.4, i, 1, 0, i
        ) for i in 1:60
    ]
    tutorial = FilterConfig(photons = (500.0, Inf), precision = (0.0, 0.007), pvalue = (1.0e-3, 1.0), psf_sigma = :auto)
    aniso = FilterConfig(photons = (200.0, 5000.0), precision = (0.0, 0.007), pvalue = (1.0e-3, 1.0), psf_sigma = :auto)
    # Very narrow frames (a 124 px figure) widen the figure; so does a word wider than it.
    cases = (
        (SMLMAnalysis._fit_overlay_title(tutorial, iso), (128, 128)),
        (SMLMAnalysis._fit_overlay_title(aniso, xy), (128, 128)),
        (SMLMAnalysis._fit_overlay_title(aniso, xy), (256, 64)),
        (SMLMAnalysis._fit_overlay_title(aniso, xy), (512, 16)),
        ("Fit: " * repeat("W", 100), (128, 128)),
    )
    for (title, (h, w)) in cases
        r = SMLMAnalysis._box_overlay_figure(
            zeros(Float32, h, w, 12), Float64[], Float64[], Int[], 7, Symbol[]; suptitle = title
        )
        width = r.fig.scene.viewport[].widths[1]
        @test width >= SMLMAnalysis._grid_figure_size(zeros(h, w, 1))[1]
        @test r.label.layoutobservables.autosize[][1] <= width - 20
        @test split(r.label.text[]) == split(title)
    end
    @test all(t -> length(t[1]) > 100, cases)
    # A theme that turns on Makie's Label word wrap does not break the overlay or its fit.
    Mk = SMLMAnalysis.CairoMakie
    Mk.with_theme(Mk.Theme(Label = (word_wrap = true,))) do
        for (title, (h, w)) in cases[[1, 4]]
            r = SMLMAnalysis._box_overlay_figure(
                zeros(Float32, h, w, 12), Float64[], Float64[], Int[], 7, Symbol[]; suptitle = title
            )
            @test r.label.layoutobservables.autosize[][1] <= r.fig.scene.viewport[].widths[1] - 20
        end
        mktempdir() do dir
            SMLMAnalysis._save_box_overlay(dir, "o.png", zeros(Float32, 128, 128, 12), Float64[], Float64[], Int[], 7, Symbol[]; suptitle = cases[2][1])
            @test isfile(joinpath(dir, "o.png"))
        end
    end
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
