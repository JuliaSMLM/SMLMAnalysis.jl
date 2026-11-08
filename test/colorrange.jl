using SMLMAnalysis
using GaussMLE
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

@testset "P-value panel with values on a bin edge" begin
    # log10 p-values -2, -1.8 and -1.79 with edges from -2 in steps of 0.04: -1.8 and -1.79
    # share a bin by exact edges but not by floating-point division. The bars and the y
    # limit now come from one count, so the panel saves and its peak is not clipped.
    cam = IdealCamera(16, 16, 0.1)
    ps = vcat(fill(0.01, 10), fill(10^-1.8, 50), fill(10^-1.79, 50))
    es = [_sigma_emitter(pvalue = p, frame = i) for (i, p) in enumerate(ps)]
    @test _quality_figure_saves(SMLMAnalysis.BasicSMLD(es, cam, length(es), 1, Dict{String, Any}()), FilterConfig())
end

@testset "Fit overlay title with open photon and precision bounds" begin
    t = SMLMAnalysis._fit_overlay_title
    es = [_sigma_emitter()]
    @test occursin("red=photons<500/>5000", t(FilterConfig(photons = (500, 5000)), es))
    @test occursin("red=photons<500  ", t(FilterConfig(photons = (500, Inf)), es))
    # Open bounds reject nothing: no label, and no round(Int, ±Inf).
    @test !occursin("photons", t(FilterConfig(photons = (-Inf, Inf)), es))
    @test occursin("red=photons>5000", t(FilterConfig(photons = (-Inf, 5000)), es))
    @test !occursin("prec", t(FilterConfig(precision = (0, Inf)), es))
    @test occursin("orange=prec>20.0nm", t(FilterConfig(precision = (0, 0.02)), es))
    # An active lower precision bound is labeled too.
    @test occursin("orange=prec<20.0nm", t(FilterConfig(precision = (0.02, Inf)), es))
    @test occursin("orange=prec<5.0nm/>20.0nm", t(FilterConfig(precision = (0.005, 0.02)), es))
end

@testset "P-value panel counts and pass percentage" begin
    # Values on an edge go in the bin that starts there: log10 p = -1.8 and -1.79 both
    # belong in [-1.8, -1.76) for edges from -2 in steps of 0.04. Dividing by the bin
    # width put -1.8 in the bin below.
    edges = collect(range(-2.0, 0.0, length = 51))
    vals = vcat(fill(log10(0.01), 10), fill(log10(10^-1.8), 50), fill(log10(10^-1.79), 50))
    c = SMLMAnalysis._bin_counts(vals, edges)
    @test c[1] == 10 && c[6] == 100 && sum(c) == 110
    # p = 1 (the top edge) goes in the last bin; values outside the edges are not counted.
    @test SMLMAnalysis._bin_counts([0.0, 0.5, -3.0], edges) == [zeros(Int, 49); 1]
    # The pass percentage uses the filter's own test: inclusive ends and the upper bound.
    pp = SMLMAnalysis._pass_percent
    @test pp([0.4, 0.4, 0.4], (0.4, 1.0)) == 100.0
    @test pp([0.1, 0.5, 0.9], (0.2, 0.8)) == round(100 / 3, digits = 1)
end

@testset "Filter figures show both sides of a filter range" begin
    rs = SMLMAnalysis._rejected_spans
    # Lower bound only, upper only, both, open, and a bound at the window edge.
    @test rs(500, Inf, 0, 1000) == (spans = [(0.0, 500.0)], thresholds = [500.0])
    @test rs(-Inf, 800, 0, 840) == (spans = [(800.0, 840.0)], thresholds = [800.0])
    @test rs(5, 20, 0, 24) == (spans = [(0.0, 5.0), (20.0, 24.0)], thresholds = [5.0, 20.0])
    @test rs(-Inf, Inf, 0, 10) == (spans = Tuple{Float64, Float64}[], thresholds = Float64[])
    @test rs(0, Inf, 0, 10).spans == Tuple{Float64, Float64}[]
    # log10 p window: p from 0.2 to 0.8 rejects on both sides.
    r = rs(log10(0.2), log10(0.8), -4.0, 0.0)
    @test r.thresholds ≈ [log10(0.2), log10(0.8)]

    # The overlay legend names an active upper p-value bound.
    t = SMLMAnalysis._fit_overlay_title
    es = [_sigma_emitter()]
    @test occursin("purple=pval<0.2/>0.8", t(FilterConfig(pvalue = (0.2, 0.8)), es))
    @test occursin("purple=pval<0.001  ", t(FilterConfig(pvalue = (0.001, 1.0)), es))

    # Each input the filter rejects on the side that was not drawn before saves its figure.
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
end

@testset "P-value panel window and bounds" begin
    pw = SMLMAnalysis._pvalue_window
    rs = SMLMAnalysis._rejected_spans
    # Both thresholds of (0.2, 0.8) lie inside the window even when the data sit near p = 1.
    w = pw(-0.0457575, (0.2, 0.8))
    @test w.lo_log ≈ log10(0.2) && w.hi_log ≈ log10(0.8)
    @test w.x_lo <= w.lo_log < w.hi_log < 0
    @test all(t -> w.x_lo <= t <= 0, rs(w.lo_log, w.hi_log, w.x_lo, 0.0).thresholds)
    # An ordinary lower bound inside the data window leaves the window as it was.
    w = pw(-5.0, (0.001, 1.0))
    @test w.x_lo == -5.0 && w.lo_log ≈ -3.0 && w.hi_log == Inf
    # Upper bounds at or below 0 pass no positive p-value: the whole window is rejected and
    # no logarithm of a non-positive number is taken.
    for b in ((-Inf, -0.1), (0.0, 0.0))
        w = pw(-2.0, b)
        @test w.hi_log == w.x_lo == -2.0
        @test rs(w.lo_log, w.hi_log, w.x_lo, 0.0).spans == [(-2.0, 0.0)]
    end
    @test pw(-2.0, nothing) == (x_lo = -2.0, lo_log = -Inf, hi_log = Inf)

    cam = IdealCamera(16, 16, 0.1)
    es = [_sigma_emitter(pvalue = p, frame = i) for (i, p) in enumerate((0.1, 0.5, 0.9))]
    smld = SMLMAnalysis.BasicSMLD(es, cam, 3, 1, Dict{String, Any}())
    for cfg in (FilterConfig(pvalue = (-Inf, -0.1)), FilterConfig(pvalue = (0, 0)), FilterConfig(pvalue = (0.2, 0.8)))
        @test _quality_figure_saves(smld, cfg)
    end
end
