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
