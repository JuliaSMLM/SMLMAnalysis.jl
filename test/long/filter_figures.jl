# The filter figures over many bound settings (the review cases stay in Core, colorrange.jl).
using SMLMAnalysis
using GaussMLE
using Random
using Test

include(joinpath(dirname(@__DIR__), "filter_checks.jl"))

@testset "Filter figures across bound settings" begin
    cam = IdealCamera(16, 16, 0.1)
    smld_of(es) = SMLMAnalysis.BasicSMLD(es, cam, length(es), 1, Dict{String, Any}())
    rng = Random.Xoshiro(7)
    es = _sweep_fits(rng)
    n = length(es)
    smld = smld_of(es)
    removed(cfg) = n - length(SMLMAnalysis._filter_smld(smld, cfg).emitters)
    shaded = _shaded

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

    e3 = [
        SMLMAnalysis.Emitter3DFit{Float64}(
            0.8, 0.8, 0.2 + 0.1 * rand(rng), 1000.0, 10.0, 0.01, 0.01, 0.02 + 0.03 * rand(rng), 30.0, 1.0; frame = i
        ) for i in 1:500
    ]
    for cfg in (FilterConfig(z = (0.22, 0.28)), FilterConfig(z = (-0.4, -0.1)), FilterConfig(sigma_z = (0.03, 0.04)), FilterConfig(sigma_z = (0, 0.03)))
        _check_quality_panels(e3, cfg)
        _check_quality_figure(e3, cfg)
        @test _quality_figure_saves(smld_of(e3), cfg)
    end
end

@testset "Filter figures read back from Makie: review cases" begin
    es = _sweep_fits(Random.Xoshiro(7))
    for cfg in (FilterConfig(photons = (-Inf, 6000)), FilterConfig(pvalue = (1.0e-9, 1)))
        _check_quality_figure(es, cfg)
    end
    _check_quality_figure([_fit2d(photons = 1000.01)], FilterConfig(photons = (-Inf, 1000.01)))
    for (ps, b) in (([0.0, 0.5], (0.1, 0.0)), ([0.0], (NaN, 0.0)))
        _check_quality_figure([_fit2d(pvalue = q, i = i) for (i, q) in enumerate(ps)], FilterConfig(pvalue = b))
    end
    pr = vcat([_fit2d(px = 0.01, py = 0.01, i = i) for i in 1:99], [_fit2d(px = 0.5, py = 0.5, i = 100)])
    _check_quality_figure(pr, FilterConfig(precision = (0, 0.02)))
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
