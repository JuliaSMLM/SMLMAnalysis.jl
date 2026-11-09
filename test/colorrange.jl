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
    edges = collect(range(-2.0, 0.0, length = 51))
    none = Tuple{Float64, Float64}[]
    # A value on an edge goes in the bar that starts there: -1.8 and -1.79 both belong in
    # [-1.8, -1.76). Dividing by the bin width put -1.8 in the bar below.
    vals = vcat(fill(log10(0.01), 10), fill(log10(10^-1.8), 50), fill(log10(10^-1.79), 50))
    c = SMLMAnalysis._panel_counts(vals, edges, nothing, none)
    @test c.counts[1] == 10 && c.counts[6] == 100 && sum(c.counts) == 110 && c.beyond == 0
    # The top edge goes in the last bar; values beyond the view on an unshaded side are not drawn.
    c = SMLMAnalysis._panel_counts([0.0, 0.5, -3.0, NaN], edges, nothing, none)
    @test c.counts == [zeros(Int, 49); 1] && c.beyond == 0
    # A value on an upper bound passes the filter, so it goes in the bar below the bound.
    e = [0.0, 1.0, 2.0, 3.0]
    c = SMLMAnalysis._panel_counts([2.0], e, (-Inf, 2.0), [(2.0, 3.0)])
    @test c.counts == [0, 1, 0]
    # A value on a lower bound passes and goes in the bar that starts at the bound.
    c = SMLMAnalysis._panel_counts([1.0], e, (1.0, Inf), [(0.0, 1.0)])
    @test c.counts == [0, 1, 0]
    # Beyond the view on a shaded side: rejected, drawn in that end's bar.
    c = SMLMAnalysis._panel_counts([-5.0, -Inf, 9.0], e, (1.0, 2.0), [(0.0, 1.0), (2.0, 3.0)])
    @test c.counts == [2, 0, 1] && c.beyond == 3
    # The pass percentage uses the filter's own test: inclusive ends and the upper bound.
    pp = SMLMAnalysis._pass_percent
    @test pp([0.4, 0.4, 0.4], (0.4, 1.0)) == 100.0
    @test pp([0.1, 0.5, 0.9], (0.2, 0.8)) == round(100 / 3, digits = 1)
end

@testset "Panel geometry" begin
    g = SMLMAnalysis._panel_geometry
    # A bound beyond the window that rejects values extends the view, with margin, and is an edge.
    r = g((0.0, 4463.0), (10.0, 9000.0), [(-Inf, 6000.0)])
    @test r.view[1] == 0.0 && 6000 < r.view[2] && r.lines == [6000.0] && 6000.0 in r.edges
    # Sign-safe margin: a positive lower z bound below the data stays inside the view.
    r = g((200.0, 300.0), (200.0, 300.0), [(100.0, 400.0)])
    @test r.view[1] < 100 && r.view[2] > 400 && r.lines == [100.0, 400.0]
    # A bound that rejects nothing and lies far from the data is left off.
    r = g((0.0, 4000.0), (10.0, 5000.0), [(200.0, 1.0e20)])
    @test r.lines == [200.0] && r.view[2] == 4000.0
    # A bound on the window edge gets margin, so it is strictly inside.
    r = g((0.0, 10.0), (0.0, 10.0), [(0.0, Inf)])
    @test r.view[1] < 0.0 && r.lines == [0.0]
    # ...unless it rejects nothing and sits at the panel's floor (photons or precision of 0).
    r = g((0.0, 10.0), (0.5, 10.0), [(0.0, Inf)]; floor = 0.0)
    @test r.view == (0.0, 10.0) && isempty(r.lines)
    # Edges increase strictly and no bar straddles a bound.
    r = g((0.0, 1.0), (0.0, 1.0), [(0.3001, 0.7)])
    @test issorted(r.edges) && allunique(r.edges) && 0.3001 in r.edges && 0.7 in r.edges
    # Spans: whole view for a range that rejects everything; only the shaded sides.
    s = SMLMAnalysis._panel_spans
    @test s((0.0, 10.0), (Inf, Inf), (:lo, :hi)) == [(0.0, 10.0)]
    @test s((0.0, 10.0), (-Inf, -Inf), (:lo, :hi)) == [(0.0, 10.0)]
    @test s((0.0, 10.0), (2.0, 8.0), (:hi,)) == [(8.0, 10.0)]
    @test s((0.0, 10.0), (-Inf, Inf), (:lo, :hi)) == Tuple{Float64, Float64}[]
end

# What each fit_quality.png panel draws, checked against the filter. Every bound that rejects a
# value lies strictly inside the view, on a bar edge. In each shaded span the bars hold exactly
# the values that side of the range rejects, those beyond the view on that side included. Every
# value inside the view is drawn once.
function _check_quality_panels(emitters, cfg)
    panels = SMLMAnalysis._quality_panels(emitters, cfg)
    for (key, p) in panels
        vlo, vhi = p.view
        @test vlo < vhi
        @test issorted(p.edges) && allunique(p.edges) && first(p.edges) == vlo && last(p.edges) == vhi
        for s in p.series
            vals = filter(!isnan, s.values)
            @test sum(s.counts) == count(v -> vlo <= v <= vhi, vals) + s.beyond
            s.bounds === nothing && continue
            lo, hi = s.bounds
            for (b, rejects) in ((lo, any(<(lo), vals)), (hi, any(>(hi), vals)))
                if isfinite(b) && rejects && abs(b) < 1.0e30
                    @test vlo < b < vhi
                    @test b in p.edges
                    @test b in s.lines || b in s.dashed
                end
            end
            expected = Tuple{Float64, Float64}[]
            rejected = Int[]
            # A NaN end rejects everything (it fails every comparison).
            if :lo in s.sides && (isnan(lo) || lo > vlo)
                push!(expected, (vlo, isnan(lo) ? vhi : min(lo, vhi)))
                push!(rejected, isnan(lo) ? length(vals) : count(<(lo), vals))
            end
            if :hi in s.sides && (isnan(hi) || hi < vhi)
                push!(expected, (isnan(hi) ? vlo : max(hi, vlo), vhi))
                push!(rejected, isnan(hi) ? length(vals) : count(>(hi), vals))
            end
            @test s.spans == expected
            for ((a, b), n) in zip(s.spans, rejected)
                drawn = sum((s.counts[i] for i in eachindex(s.counts) if p.edges[i] >= a && p.edges[i + 1] <= b); init = 0)
                @test drawn == n
            end
        end
    end
    return panels
end

_fit2d(; photons, bg = 5.0, σ = 0.13, px = 0.01, py = 0.01, pvalue = 0.5, i = 1) =
    GaussMLE.Emitter2DFitSigma{Float64}(0.8, 0.8, photons, bg, σ, px, py, 0.0, 30.0, 0.5, 0.002, pvalue, i, 1, 0, i)

@testset "Filter figures draw what the filter rejects" begin
    cam = IdealCamera(16, 16, 0.1)
    smld_of(es) = SMLMAnalysis.BasicSMLD(es, cam, length(es), 1, Dict{String, Any}())
    rng = Random.Xoshiro(7)
    n = 2000
    # Mostly ordinary fits; a bright tail above the 98th percentile, a few wide PSFs, poor
    # precisions and p-values far below every view.
    es = [
        _fit2d(
            photons = i <= 60 ? 6000 + 3000 * rand(rng) : 1000 + 3000 * rand(rng),
            σ = i % 50 == 0 ? 0.25 : 0.12 + 0.03 * rand(rng),
            px = i % 97 == 0 ? 0.5 : 0.005 + 0.02 * rand(rng),
            py = 0.005 + 0.02 * rand(rng),
            pvalue = i % 40 == 0 ? 1.0e-12 : (i % 333 == 0 ? 0.0 : rand(rng)),
            i = i
        ) for i in 1:n
    ]
    smld = smld_of(es)
    removed(cfg) = n - length(SMLMAnalysis._filter_smld(smld, cfg).emitters)
    spans_total(s, edges) = sum(
        (s.counts[i] for i in eachindex(s.counts) if any(sp -> edges[i] >= sp[1] && edges[i + 1] <= sp[2], s.spans));
        init = 0
    )

    # Review case: an upper photon bound above the 98th percentile. Every fit it rejects is
    # drawn in the shaded band, and the passing fits between the percentile and the bound show.
    cfg = FilterConfig(photons = (-Inf, 6000))
    p = _check_quality_panels(es, cfg)[:photons]
    @test spans_total(p.series[1], p.edges) == removed(cfg) == 60
    @test sum(p.series[1].counts) == n
    # Review case: a lower p-value bound far below the data's 2% quantile; p = 0 included.
    cfg = FilterConfig(pvalue = (1.0e-9, 1))
    p = _check_quality_panels(es, cfg)[:pvalue]
    @test spans_total(p.series[1], p.edges) == removed(cfg) == 50 + 6
    for cfg in (
            FilterConfig(photons = (2000, Inf)), FilterConfig(photons = (1500, 3500)),
            FilterConfig(photons = (1.0e6, Inf)), FilterConfig(photons = (-Inf, -1)),
            FilterConfig(photons = (200, 1.0e20)), FilterConfig(photons = (0, Inf)),
            FilterConfig(photons = (3000, 2000)),
            FilterConfig(pvalue = (0.2, 0.8)), FilterConfig(pvalue = (0, 0)), FilterConfig(pvalue = (-Inf, -0.1)),
            FilterConfig(pvalue = (1.0e-3, 1)), FilterConfig(pvalue = (1.0, 1.0)),
            FilterConfig(psf_sigma = :auto), FilterConfig(psf_sigma = (0.13, 0.14)),
            FilterConfig(psf_sigma = (0.3, Inf)),
        )
        p = _check_quality_panels(es, cfg)
        key = cfg.photons !== nothing ? :photons : cfg.pvalue !== nothing ? :pvalue : :psf
        lo, hi = p[key].series[1].bounds === nothing ? (-Inf, Inf) : p[key].series[1].bounds
        lo <= hi && @test spans_total(p[key].series[1], p[key].edges) == removed(cfg)
        @test _quality_figure_saves(smld, cfg)
    end
    # Precision: the upper bound is shaded and its rejected components drawn there; the lower
    # bound (on max(σx, σy)) is a dashed line.
    for cfg in (FilterConfig(precision = (0.01, 0.02)), FilterConfig(precision = (0.02, Inf)), FilterConfig(precision = (0, 0.015)))
        p = _check_quality_panels(es, cfg)[:precision]
        lo = cfg.precision[1] * 1000
        @test all(s -> (lo > 0) == (lo in s.dashed), p.series)
        @test _quality_figure_saves(smld, cfg)
    end

    # Anisotropic PSF: per-axis :auto bounds differ, explicit bounds are shared.
    xy = [
        GaussMLE.Emitter2DFitSigmaXY{Float64}(
            0.8, 0.8, 1000.0, 5.0, 0.12 + 0.02 * rand(rng) + (i % 30 == 0 ? 0.1 : 0.0),
            0.15 + 0.02 * rand(rng), 0.01, 0.01, 0.0, 30.0, 0.5, 0.002, 0.002, rand(rng), i, 1, 0, i
        ) for i in 1:500
    ]
    for cfg in (FilterConfig(psf_sigma = :auto), FilterConfig(psf_sigma = (0.13, 0.16)), FilterConfig(psf_sigma = (0.2, 0.4)))
        p = _check_quality_panels(xy, cfg)[:psf]
        @test all(s -> s.bounds !== nothing, p.series)
        @test _quality_figure_saves(smld_of(xy), cfg)
    end

    # 3D, all z in 0.2-0.3 µm. Positive bounds below the data stay inside the view (the old
    # 1.1x widening moved them toward zero). z = (0.1, 0.4) rejects nothing and lies more than
    # a window width away, so it is not drawn and nothing is shaded.
    e3 = [
        SMLMAnalysis.Emitter3DFit{Float64}(
            0.8, 0.8, 0.2 + 0.1 * rand(rng), 1000.0, 10.0, 0.01, 0.01, 0.02 + 0.03 * rand(rng), 30.0, 1.0; frame = i
        ) for i in 1:500
    ]
    p = _check_quality_panels(e3, FilterConfig(z = (0.15, 0.35)))[:z]
    @test p.view[1] < 150 && p.view[2] > 350 && p.series[1].lines ≈ [150.0, 350.0]
    p = _check_quality_panels(e3, FilterConfig(z = (0.1, 0.4)))[:z]
    @test isempty(p.series[1].lines) && isempty(p.series[1].spans)
    for cfg in (FilterConfig(z = (0.22, 0.28)), FilterConfig(z = (-0.4, -0.1)), FilterConfig(sigma_z = (0.03, 0.04)), FilterConfig(sigma_z = (0, 0.03)))
        _check_quality_panels(e3, cfg)
        @test _quality_figure_saves(smld_of(e3), cfg)
    end
end

@testset "Fit overlay title" begin
    t = SMLMAnalysis._fit_overlay_title
    es = [_sigma_emitter()]
    @test occursin("red=photons<500/>5000", t(FilterConfig(photons = (500, 5000)), es))
    @test occursin("red=photons<500  ", t(FilterConfig(photons = (500, Inf)), es))
    # Open bounds reject nothing: no label.
    @test !occursin("photons", t(FilterConfig(photons = (-Inf, Inf)), es))
    @test occursin("red=photons>5000", t(FilterConfig(photons = (-Inf, 5000)), es))
    # A finite bound too large for an Int is labeled, not an InexactError.
    @test occursin("red=photons<200/>1.0e20", t(FilterConfig(photons = (200, 1.0e20)), es))
    @test !occursin("prec", t(FilterConfig(precision = (0, Inf)), es))
    @test occursin("orange=prec>20.0nm", t(FilterConfig(precision = (0, 0.02)), es))
    @test occursin("orange=prec<20.0nm", t(FilterConfig(precision = (0.02, Inf)), es))
    @test occursin("orange=prec<5.0nm/>20.0nm", t(FilterConfig(precision = (0.005, 0.02)), es))
    # The p-value label names the range, which a NaN p-value also fails.
    @test occursin("purple=pval not in [0.2, 0.8]", t(FilterConfig(pvalue = (0.2, 0.8)), es))
    nan_p = _sigma_emitter(pvalue = NaN)
    cfg = FilterConfig(pvalue = (0.001, 1.0))
    @test !SMLMAnalysis._in_range(nan_p.pvalue, cfg.pvalue)
    @test occursin("purple=pval not in [0.001, 1.0]", t(cfg, [nan_p]))
    # PSF bounds are named as the filter applies them; a skipped degenerate :auto has no label.
    @test occursin("cyan=psf_σ not in [100.0, 200.0]nm", t(FilterConfig(psf_sigma = (0.1, 0.2)), es))
    @test occursin("cyan=psf_σ not in [117.0, 143.0]nm", t(FilterConfig(psf_sigma = :auto), es))
    @test !occursin("cyan", t(FilterConfig(psf_sigma = :auto), [_sigma_emitter(σ = 0.0)]))
    xy = GaussMLE.Emitter2DFitSigmaXY{Float64}(0.8, 0.8, 1000.0, 5.0, 0.13, 0.14, 0.01, 0.01, 0.0, 30.0, 0.5, 0.002, 0.002, 0.4, 1, 1, 0, 1)
    @test occursin("cyan=psf_σx not in [100.0, 200.0]nm or σy not in [100.0, 200.0]nm", t(FilterConfig(psf_sigma = (0.1, 0.2)), [xy]))
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
