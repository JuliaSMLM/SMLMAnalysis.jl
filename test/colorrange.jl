using SMLMAnalysis
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
